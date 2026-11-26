#ifndef THCALLBACK_SURF_EAMCC_H_
#define THCALLBACK_SURF_EAMCC_H_

#include "particledatabase.hpp"
#include "random.hpp"
#include "globals.h"
#include "funct.h"
#include "constants.hpp"
#include "geometry.hpp"
#include "error.hpp"
#include "CallbackThreading.h"
#include "SurfaceEventLedger.h"
#include <atomic>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

// Safety constants (same as THCallback.h)
#define MAX_GENERATION_DEPTH 5
#define MAX_PARTICLE_COUNT 10000000
#ifndef MAX_SECONDARY_PARTICLES_PER_PASS
#define MAX_SECONDARY_PARTICLES_PER_PASS 1000000
#endif
#define MIN_VELOCITY_THRESHOLD 1e3
#define MAX_VELOCITY_SCALE SPEED_C
// Per-species cap on secondary electrons from one impact. The Fortran passes different
// values for the two projectile types: 26 for ion impacts (ion_surf_coll.f:101) and 3 for
// electron impacts (e_surf_coll.f:91). The electron cap is effectively never binding, since
// etaSEC = ETA1*exp(CTsec*(1-cos theta)) <= 1.13*exp(1.005) = 3.09.
// Defaults, overridable per run through physics.surfaceCollisions. The ion cap and the
// secondary energy are configurable (SimulationParameters::getSurfaceMaxSecondaryElectrons
// and ::getSurfaceSecondaryElectronEnergyEV); the ELECTRON-impact cap is not, because
// etaSEC = ETA1*exp(CTsec*(1-cos theta)) <= 1.13*exp(1.005) = 3.09 makes 3 unreachable, so
// a control for it could never change a result.
#define DEFAULT_MAX_SECONDARY_ELECTRONS 26
#define MAX_SECONDARY_ELECTRONS_E 3
#define DEFAULT_SECONDARY_ELECTRON_ENERGY 10.0  // eV

// Backscattered-electron energy spectrum, Fubiani/de Esch/Simonin/Hemsworth,
// Phys. Rev. ST Accel. Beams 11, 014202 (2008), Sec. II A.
// ALPHA_BACK: exponent alpha of Eq. (1); the paper notes alpha = 2.2 for all energies.
// OMEGA_BACK: Omega of Eq. (6), B_theta(E0, pi/2, pi/2) = Omega, from Matsukawa et al.
#define ALPHA_BACK 2.2
#define OMEGA_BACK 0.55

// Generation number offset to distinguish surface-generated from volume-generated particles
// Volume-generated: gen + 1 (e.g., 0->1, 1->2, 2->3, ...)
// Surface-generated: gen + SURFACE_GENERATION_OFFSET (e.g., 0->101, 1->102, 2->103, ...)
// 
// To distinguish in code:
//   - Check if gen >= SURFACE_GENERATION_OFFSET for surface-generated
//   - Use gen % 100 to get effective generation (0-4) for diagnostics
//   - Volume: gen % 100 = 1,2,3,4,5...
//   - Surface: gen % 100 = 1,2,3,4,5... (same effective gen, but gen >= 101)
#define SURFACE_GENERATION_OFFSET 101

// Default number of macroparticles used to carry the sampled secondary-electron yield.
// See THCallback_surf_EAMCC::generateSecondaryElectrons for what this does and why it is
// mean-preserving. 0 restores the original one-macroparticle-per-electron behaviour.
#define DEFAULT_SECONDARY_ELECTRON_BUNDLE 4

namespace thcallback_surf_detail {

/*! \brief Solid index of the impact this thread is currently processing.
 *
 *  Scratch state shared between operator() and the emission records made further down the
 *  same call, rather than a parameter, because it would otherwise have to be threaded
 *  through four handler signatures that already take eight or nine arguments. It is
 *  thread-local because every worker thread shares one callback object, so a plain member
 *  would be overwritten by whichever thread hit a surface most recently and emissions would
 *  be booked against the wrong electrode.
 */
inline uint32_t& current_hit_solid() {
    static thread_local uint32_t solid = 0;
    return solid;
}

} // namespace thcallback_surf_detail

/**
 * @class THCallback_surf_EAMCC
 * @brief Handles surface collisions with EAMCC secondary particle generation
 * 
 * This class implements the EAMCC (Electron and Ion Monte Carlo Collision) 
 * surface interaction model for secondary particle generation.
 * Based on the Fortran routines in EAMCCsecondaries folder.
 * 
 * Features:
 * - Electron backscattering and secondary emission
 * - Ion backscattering and secondary electron emission
 * - Angular-dependent emission probabilities
 * - Proper energy and charge conservation
 */
class THCallback_surf_EAMCC : public TrajectoryEndCallback {
private:
    Geometry &_geom;
    ParticleDataBase3D* _pdb;
    double _mass;
    double _minimum_z;
    bool _debugprint;

    // Signed surface energy/charge balance. Not owned; may be null, in which case net
    // accounting is simply unavailable and the gross columns are unaffected.
    SurfaceEventLedger* _ledger;

    // Warn-once flags and the running secondary count. Atomic because every worker thread
    // shares this one callback object.
    std::atomic<bool> _secondary_limit_warned;
    std::atomic<bool> _electron_limit_warned;

    // Particle-database size is tracked here rather than read back from _pdb->size():
    // std::vector::size() raced against the push_back() inside add_particle() is a data
    // race, and this counter answers the same question without one.
    size_t _initial_particle_count;
    std::atomic<size_t> _secondaries_created;

    // Max macroparticles used to represent one impact's secondary-electron yield; 0 means
    // emit one macroparticle per electron, as the original code did.
    int _secondary_electron_bundle;

    // Per-impact secondary-electron cap for ion/neutral impacts, and the birth energy of a
    // true secondary [eV]. Both come from the run configuration.
    int _max_secondary_electrons;
    double _secondary_electron_energy;

    /*! \brief Have we hit a particle cap?
     *
     *  Kept identical to THCallback.h: the two callbacks share the same ceilings and must
     *  agree, or the gas and surface cascades stop at different points.
     *
     *   - MAX_SECONDARY_PARTICLES_PER_PASS bounds the secondaries CREATED in this pass.
     *     Primaries are not part of it: they are an input, not cascade product.
     *   - MAX_PARTICLE_COUNT bounds the TOTAL population, primaries included.
     *
     *  The previous form compared primaries plus secondaries against the SECONDARY cap, so
     *  a run configured with as many primaries as that cap generated no secondaries at all.
     */
    bool secondaryBudgetExhausted() const {
        const size_t created = _secondaries_created.load(std::memory_order_relaxed);
        return created >= MAX_SECONDARY_PARTICLES_PER_PASS ||
               _initial_particle_count + created >= MAX_PARTICLE_COUNT;
    }

    /*! \brief Emit a warn-once message, at most one thread at a time. */
    void warnOnce(std::atomic<bool>& flag, const char* message) {
        bool expected = false;
        if (flag.compare_exchange_strong(expected, true)) {
            CB_LOG_LOCK();
            logfile << "WARNING: " << message << endl;
        }
    }


    // Helper function for safe velocity scaling
    bool validateAndScaleVelocity(double vel, double minvel, const double originalVel[3], double scaledVel[3]) {
        if (vel < MIN_VELOCITY_THRESHOLD || !std::isfinite(vel)) {
            if (debug) logfile << "WARNING: Invalid velocity in secondary generation: " << vel << endl;
            return false;
        }
        
        double scale_factor = std::min(1.0, minvel/vel);
        
        for (int i = 0; i < 3; i++) {
            scaledVel[i] = scale_factor * originalVel[i];
            
            if (!std::isfinite(scaledVel[i]) || abs(scaledVel[i]) > MAX_VELOCITY_SCALE) {
                if (debug) logfile << "WARNING: Invalid scaled velocity component " << i << ": " << scaledVel[i] << endl;
                scaledVel[i] = 0.0;
            }
        }
        
        return true;
    }
    
    // Position validation
    bool validatePosition(const Vec3D& loc) {
        if (!std::isfinite(loc[0]) || !std::isfinite(loc[1]) || !std::isfinite(loc[2])) {
            logfile << "ERROR: Invalid particle position in secondary generation" << endl;
            return false;
        }
        return true;
    }

    /**
     * @brief Calculate backscattering probability for electrons (ETA0)
     * Based on ORNL Redbook data and experimental fits
     * @param E0 Incident electron energy in eV
     * @return Backscattering probability (0-1)
     */
    double ETA0(double E0) {
        // Data from ORNL Redbook Vol 1, D-4 for Cu
        // Energy points in eV
        const double E[] = {0., 100., 200., 500., 1000., 2000., 4000., 7000.,
                           10000., 20000., 40000., 70000., 100000.,
                           200000., 400000., 700000., 1.0e6, 4.0e6};
        // Backscattering probabilities
        const double eta[] = {0.0, 0.124, 0.161, 0.262, 0.266, 0.282, 0.315, 0.322,
                             0.323, 0.318, 0.309, 0.299, 0.285,
                             0.275, 0.260, 0.240, 0.210, 0.11};
        const int n = 18;
        
        if (E0 <= E[0]) return eta[0];
        if (E0 >= E[n-1]) return eta[n-1];
        
        // Linear interpolation
        for (int i = 1; i < n; i++) {
            if (E0 > E[i-1] && E0 <= E[i]) {
                return eta[i-1] + (E0 - E[i-1]) * (eta[i] - eta[i-1]) / (E[i] - E[i-1]);
            }
        }
        
        return eta[n-1];
    }

    /**
     * @brief Calculate CT factor for backscattered electrons
     * Based on P-F Staub et al. [J. of Phys. D 27, 1533 (1994)]
     * @param E0 Incident electron energy in eV
     * @param eta0 Backscattering probability at normal incidence
     * @return CT factor for angular dependence
     */
    double CTback(double E0, double eta0) {
        double kappa = 1.0 - exp(-1.83 * pow(E0 * 1e-3, 0.25));
        return kappa * log(1.0 / eta0);
    }

    /**
     * @brief Calculate true secondary electron emission yield (ETA1)
     * Based on CERN LHC Project Report 472 (fully conditioned copper surfaces)
     * @param Ep Incident particle energy in eV
     * @return Secondary electron yield
     */
    double ETA1(double Ep) {
        // Parameters from CERN report
        const double s = 1.35;
        const double Emax = 318.0;  // eV
        const double DeltaMax = 1.13;
        const double A0 = 20.69989;
        const double A1 = -7.07605;
        const double A2 = 0.483547;
        const double E0 = 56.914686;  // eV
        
        double EpEm = Ep / Emax;
        double deltaS = DeltaMax * s * EpEm / (s - 1.0 + pow(EpEm, s));
        
        if (EpEm >= 1.0) {
            return deltaS;
        }
        
        double LnEpE0 = log(Ep + E0);
        double f = A0 + A1 * LnEpE0 + A2 * LnEpE0 * LnEpE0;
        f = exp(f);
        double deltaT = deltaS / (1.0 - f);
        return deltaT;
    }

    /**
     * @brief Calculate CT factor for true secondary electron emission
     * Based on Koshikawa et al. [J. of Phys. D 6, 1360 (1973)]
     * @param E0 Incident electron energy in eV
     * @return CT factor for angular dependence
     */
    double CTsec(double E0) {
        const double E[] = {0.0, 500.0, 1000.0, 2000.0, 3000.0, 6000.0, 10000.0};
        const double CT[] = {0.692, 0.692, 0.731, 0.835, 0.956, 0.971, 1.005};
        const int n = 7;
        
        if (E0 <= E[0]) return CT[0];
        if (E0 >= E[n-1]) return CT[n-1];
        
        // Linear interpolation
        for (int i = 1; i < n; i++) {
            if (E0 > E[i-1] && E0 <= E[i]) {
                return CT[i-1] + (E0 - E[i-1]) * (CT[i] - CT[i-1]) / (E[i] - E[i-1]);
            }
        }
        
        return CT[n-1];
    }

    /**
     * @brief Table I fit parameters for the backscattered-electron energy spectrum
     *
     * Fubiani et al., Phys. Rev. ST Accel. Beams 11, 014202 (2008), Table I. Fitted to
     * Matsukawa et al. (10, 20 keV) and Sternglass et al. (2, 370 keV). The paper assumes
     * negligible variation of eta(Ehat) below 2 keV and above 370 keV, with linear
     * interpolation in between.
     *
     * The tabulated tau satisfies Eq. (6), B_theta(E0, pi/2, pi/2) = OMEGA_BACK, to better
     * than 1% at all four points (0.2*e^1.02 = 0.555, 0.24*e^0.824 = 0.547,
     * 0.265*e^0.73 = 0.550, 0.273*e^0.70 = 0.550), which is a useful check on the table.
     *
     * @param E0_keV Incident electron energy in keV
     * @param B0 Fit parameter B0 (out)
     * @param p Fit parameter p (out)
     * @param tau Fit parameter tau (out)
     */
    void backscatterFitParams(double E0_keV, double& B0, double& p, double& tau) {
        const double E[]    = {2.0,  10.0,  20.0,  370.0};
        const double vB0[]  = {0.20, 0.24,  0.265, 0.273};
        const double vp[]   = {0.32, 0.27,  0.27,  0.27};
        const double vtau[] = {0.51, 0.412, 0.365, 0.35};
        const int n = 4;

        if (E0_keV <= E[0])   { B0 = vB0[0];   p = vp[0];   tau = vtau[0];   return; }
        if (E0_keV >= E[n-1]) { B0 = vB0[n-1]; p = vp[n-1]; tau = vtau[n-1]; return; }

        for (int i = 1; i < n; i++) {
            if (E0_keV <= E[i]) {
                double f = (E0_keV - E[i-1]) / (E[i] - E[i-1]);
                B0  = vB0[i-1]  + f * (vB0[i]  - vB0[i-1]);
                p   = vp[i-1]   + f * (vp[i]   - vp[i-1]);
                tau = vtau[i-1] + f * (vtau[i] - vtau[i-1]);
                return;
            }
        }
        B0 = vB0[n-1]; p = vp[n-1]; tau = vtau[n-1];
    }

    /**
     * @brief Sample energy and emission direction of a backscattered electron
     *
     * Reconstruction of the Fortran routine ebackdist_anal, which is called at
     * e_surf_coll.f:65 (the live path, since option=2 is hardcoded at e_surf_coll.f:55)
     * but is itself absent from EAMCCsecondaries/. Implements Eqs. (1)-(7) of Fubiani
     * et al., Phys. Rev. ST Accel. Beams 11, 014202 (2008), Sec. II A:
     *
     *   eta(Ehat) = S exp[ -( K / (1 - gamma*Ehat^alpha) )^p ]              (1)
     *   S         = eta_b0 exp(K^p)                                         (2)
     *   gamma     = 1 - exp[ -6 |ln B_theta|^(-3/2) ]                       (3)
     *   K         = 70 |ln B_theta|^4                                       (4)
     *   B_theta   = B0 exp[tau(1-cos theta1)] exp[tau(1-cos theta2)]        (5)
     *   Ehat      = { (1/gamma)[ 1 - K/ln^(1/p)(S/P) ] }^(1/alpha)          (7)
     *
     * eta(Ehat) is the INTEGRATED (cumulative) spectrum, not the differential one:
     * eta(0) = S exp(-K^p) = eta_b0 is the total backscatter probability, and Fig. 2
     * plots g = -d(eta)/d(Ehat).
     *
     * Because this routine runs only after backscattering has already been decided by
     * Eq. (8), the energy must come from the distribution CONDITIONED on backscattering,
     * i.e. eta normalized by eta(0) -- which is what the paper means by "normalizing and
     * inverting Eq. (1)". Normalizing replaces S by S/eta_b0 = exp(K^p), so eta_b0
     * cancels out of Eq. (7) and P is uniform on (0,1]. Using Eq. (7) with the
     * un-normalized S instead would drive the bracket negative for every P > eta_b0
     * (roughly 70% of draws, since eta_b0 ~ 0.3) and return NaN, so the normalized
     * reading is the intended one.
     *
     * The emission direction is drawn FIRST because the scattering angle theta2 feeds
     * back into the spectrum through B_theta (Eq. 5). That coupling is precisely why the
     * original is a self-contained routine that builds the particle itself instead of
     * deferring to make_sec_wall_part: energy and angle are correlated and cannot be
     * drawn independently.
     *
     * @param E0 Incident electron energy in eV
     * @param theta1 Angle of incidence from the surface normal (radians)
     * @param Eout Backscattered electron energy in eV (out)
     * @param theta2 Emission angle from the surface normal (radians) (out)
     * @param phi Azimuthal emission angle (radians) (out)
     * @return true if a valid energy was sampled
     */
    bool sampleBackscatteredElectron(double E0, double theta1,
                                     double& Eout, double& theta2, double& phi) {
        // Emission direction. The paper states the backscattered electron is reemitted
        // "in an arbitrary direction {theta2, phi}", theta2 in [0, pi/2], phi in [0, 2pi],
        // explicitly assuming isotropic scattering with no preferred direction as a
        // function of the incoming angle theta1. The paper does not state the measure in
        // theta2, so this follows the same convention the rest of this code uses to draw
        // an emission direction (make_sec_wall_part.f:41-45).
        // Two independent deviates. The previous code asked a one-dimensional MTRandom for
        // two numbers, so rnd[1] -- and therefore phi -- was whatever happened to be on the
        // stack; see ThreadRng::fill in CallbackThreading.h.
        theta2 = asin(cbthread::rng().uniform());
        phi = 2.0 * M_PI * cbthread::rng().uniform();

        double eta_b0 = ETA0(E0);
        if (eta_b0 <= 0.0) return false;

        double B0, p, tau;
        backscatterFitParams(E0 * 1e-3, B0, p, tau);

        // Eq. (5). Both angles enter symmetrically, which is how the spectrum peak moves
        // toward high Ehat at grazing incidence and grazing scattering (Fig. 2).
        double B_theta = B0 * exp(tau * (1.0 - cos(theta1)))
                            * exp(tau * (1.0 - cos(theta2)));
        if (!(B_theta > 0.0) || B_theta >= 1.0) return false;

        double lnB = fabs(log(B_theta));
        if (lnB < 1e-12) return false;

        double gamma = 1.0 - exp(-6.0 * pow(lnB, -1.5));   // Eq. (3)
        double K     = 70.0 * pow(lnB, 4.0);               // Eq. (4)
        if (!std::isfinite(gamma) || !std::isfinite(K) || gamma <= 0.0 || K <= 0.0) {
            return false;
        }

        // Eq. (7), with the normalized amplitude S/eta_b0 = exp(K^p) so that
        // ln(S/P) becomes K^p - ln(P).
        double P = cbthread::rng().uniform();
        if (P <= 0.0) P = 1e-12;      // guard the logarithm
        if (P > 1.0)  P = 1.0;

        double Kp = pow(K, p);
        double L  = pow(Kp - log(P), 1.0/p);               // ln^(1/p)( exp(K^p)/P )
        if (!std::isfinite(L) || L <= 0.0) return false;

        double bracket = (1.0 - K/L) / gamma;
        if (bracket < 0.0) bracket = 0.0;                  // P -> 1 corresponds to Ehat = 0

        double Ehat = pow(bracket, 1.0/ALPHA_BACK);
        if (!std::isfinite(Ehat)) return false;

        // Ehat can marginally exceed 1 for very small P, because Eq. (7) has
        // Ehat_max = (1/gamma)^(1/alpha) and gamma < 1 by construction (about 1.009 at
        // 20 keV normal incidence, reached for P < 3e-4). The Fortran caller treats
        // Eout > Ek as an error condition (e_surf_coll.f:77-80); clamping is the benign
        // equivalent and keeps energy conservation intact.
        if (Ehat > 1.0) Ehat = 1.0;

        Eout = Ehat * E0;
        return true;
    }

    /**
     * @brief Calculate backscattering coefficient for ions (RN)
     * Based on ORNL Redbook Vol 3, E22
     * @param E0 Incident ion energy in keV
     * @param mass_ratio Mass ratio (for D, H2, He scaling)
     * @return Backscattering coefficient
     */
    double RN(double E0, double mass_ratio = 1.0) {
        // Adjust energy for mass scaling (velocity scaling)
        E0 = E0 / mass_ratio;
        
        // Fit parameters: log(RN) = sum a_i * log(E0)^i
        const double a[] = {-1.55626, -0.2385, -0.0673113, -0.00899714, 3.23125e-6};
        const int n = 5;
        
        if (E0 <= 0.0) return 0.0;
        
        double logRN = 0.0;
        double logE = log(E0);
        for (int i = 0; i < n; i++) {
            logRN += a[i] * pow(logE, i);
        }
        
        return exp(logRN);
    }

    /**
     * @brief Calculate energy ratio for backscattered ions (RE)
     * Based on ORNL Redbook Vol 3, E23
     * @param E0 Incident ion energy in keV
     * @param mass_ratio Mass ratio (for D, H2, He scaling)
     * @return Energy ratio RE/RN
     */
    double RE(double E0, double mass_ratio = 1.0) {
        // Adjust energy for mass scaling
        E0 = E0 / mass_ratio;
        
        // Fit parameters: log(RE) = sum a_i * log(E0)^i
        const double a[] = {-2.33586, -0.450527, -0.0732621, -0.00687502, -4.56457e-6};
        const int n = 5;
        
        if (E0 <= 0.0) return 0.0;
        
        double logRE = 0.0;
        double logE = log(E0);
        for (int i = 0; i < n; i++) {
            logRE += a[i] * pow(logE, i);
        }
        
        return exp(logRE);
    }

    /**
     * @brief Sample the charge state of a backscattered heavy particle
     *
     * Fubiani et al. 2008, Sec. II C: "a backscattered ion may suffer a change of charge
     * state. It is typically found that for proton impacts the backscattered particles are
     * predominantly neutrals (~100%-85% for backscattered energy ratios Ehat = E_kb/E0
     * ranging between 0 and 1), followed by positive ions (~0%-13%), and lastly negative
     * ions (~0%-5.5%)." EAMCC uses an average profile from Eckstein & Matschke,
     * Phys. Rev. B 14, 3231 (1976) (the paper's ref [32]).
     *
     * APPROXIMATION: ref [32] is not available here, so the branching ratios are linearly
     * interpolated in Ehat between the endpoints the paper quotes, then normalized -- the
     * quoted upper values sum to 103.5%, which reflects the precision of the quoted ranges
     * rather than a real excess. Replace this with the tabulated profile from ref [32]
     * when it can be obtained.
     *
     * The paper states the same treatment is applied to all heavy particles; note that for
     * molecular projectiles a negative charge state is physically dubious, but that is what
     * the described model does and the branch is rare (<=5.5%).
     *
     * @param Ehat Backscattered energy ratio E_kb/E0, clamped to [0,1]
     * @return Charge state in units of e: -1, 0 or +1
     */
    double sampleBackscatteredChargeState(double Ehat) {
        if (Ehat < 0.0) Ehat = 0.0;
        if (Ehat > 1.0) Ehat = 1.0;

        double f_neutral = 1.00 + Ehat * (0.850 - 1.00);   // 100% -> 85%
        double f_pos     =        Ehat *  0.130;           //   0% -> 13%
        double f_neg     =        Ehat *  0.055;           //   0% -> 5.5%

        double tot = f_neutral + f_pos + f_neg;
        if (!(tot > 0.0)) return 0.0;

        double r = cbthread::rng().uniform() * tot;

        if (r < f_neutral)         return  0.0;
        if (r < f_neutral + f_pos) return  1.0;
        return -1.0;
    }

    /**
     * @brief Calculate secondary electron yield for ion impacts (eta_ion)
     * Based on ORNL Redbook Vol 3, C10 and Svensson et al PRB 25
     * @param E0 Incident ion energy in keV
     * @param mass_ratio Mass ratio (for D, H2, He scaling)
     * @return Secondary electron yield
     */
    double eta_ion(double E0, double mass_ratio = 1.0) {
        // Adjust energy for mass scaling
        E0 = E0 / mass_ratio;
        
        // Fit parameters: log(eta_ion) = sum a_i * log(E0)^i
        const double a[] = {-1.3854, 0.544935, -0.0124299, -0.00307087, -0.000558421};
        const int n = 5;
        
        if (E0 <= 0.0) return 0.0;
        
        double log_eta = 0.0;
        double logE = log(E0);
        for (int i = 0; i < n; i++) {
            log_eta += a[i] * pow(logE, i);
        }
        
        return exp(log_eta);
    }

    /**
     * @brief Handle electron surface collision
     * @param particle Incident particle
     * @param loc Impact location
     * @param vel Incident velocity
     * @param energy Incident energy in eV
     * @param theta Impact angle with surface normal (radians)
     * @param normal Surface normal vector
     * @param tang1 First tangent vector
     * @param tang2 Second tangent vector
     */
    void handleElectronSurfaceCollision(ParticleBase *particle, const Vec3D& loc, 
                                         const Vec3D& vel, double energy, double theta,
                                         const Vec3D& normal, const Vec3D& tang1, const Vec3D& tang2) {
        // Calculate backscattering probability
        double etaBACK = ETA0(energy);
        double CT = CTback(energy, etaBACK);
        etaBACK = etaBACK * exp(CT * (1.0 - cos(theta)));
        
        // Calculate secondary electron emission probability
        double etaSEC = ETA1(energy);
        CT = CTsec(energy);
        etaSEC = etaSEC * exp(CT * (1.0 - cos(theta)));
        
        if (_debugprint) {
            CB_LOG_LOCK();
            logfile << "Electron surface collision @ E=" << energy*1e-3 << " keV, "
                    << "theta=" << theta*180.0/M_PI << " deg" << endl;
            logfile << "  etaBACK=" << etaBACK << ", etaSEC=" << etaSEC << endl;
        }

        double Elost = energy;

        // Check for backscattering
        if (cbthread::rng().uniform() <= etaBACK && Elost >= 1.0) {
            // Backscattered electron. Energy and emission direction come from the
            // reconstructed ebackdist_anal (Eqs. 1-7 of Fubiani et al. 2008). The spectrum
            // is defined relative to the INCIDENT energy, Ehat = E_kb/E0, so it is sampled
            // from `energy` rather than from the running Elost.
            double Eout, theta2, phi;
            if (sampleBackscatteredElectron(energy, theta, Eout, theta2, phi)) {
                if (Eout > Elost) Eout = Elost;
                if (Eout < 0.0) Eout = 0.0;

                // Debit Elost ONLY if the particle was really created, otherwise the energy
                // simply disappears from the books.
                if (emitSecondaryParticle(particle, loc, normal, tang1, tang2, Eout,
                                          -1.0, MASS_E/MASS_U,
                                          particle->gen() + SURFACE_GENERATION_OFFSET,
                                          theta2, phi)) {
                    Elost -= Eout;

                    if (_debugprint) {
                        CB_LOG_LOCK();
                        logfile << "  Electron backscattered with E=" << Eout*1e-3 << " keV"
                                << " (Ehat=" << Eout/energy << "), emitted at "
                                << theta2*180.0/M_PI << " deg from normal" << endl;
                    }
                } else if (_debugprint) {
                    CB_LOG_LOCK();
                    logfile << "  Backscattered electron not created; Elost unchanged" << endl;
                }
            } else if (_debugprint) {
                CB_LOG_LOCK();
                logfile << "  Backscatter energy sampling failed at E=" << energy
                        << " eV, no backscattered electron emitted" << endl;
            }
        }
        
        // Generate secondary electrons (Fortran cap for electron impacts: e_surf_coll.f:91)
        generateSecondaryElectrons(particle, loc, normal, tang1, tang2, etaSEC, Elost,
                                   MAX_SECONDARY_ELECTRONS_E);
    }

    /**
     * @brief Handle ion surface collision
     * @param particle Incident particle
     * @param loc Impact location
     * @param vel Incident velocity
     * @param energy Incident energy in eV
     * @param theta Impact angle with surface normal (radians)
     * @param normal Surface normal vector
     * @param tang1 First tangent vector
     * @param tang2 Second tangent vector
     * @param kind Projectile species, from identify_particle_species()
     */
    void handleIonSurfaceCollision(ParticleBase *particle, const Vec3D& loc,
                                    const Vec3D& vel, double energy, double theta,
                                    const Vec3D& normal, const Vec3D& tang1, const Vec3D& tang2,
                                    particle_kind kind) {
        // Velocity-equivalence scaling for the RN, RE and eta_ion fits.
        //
        // All three fits are tabulated for atomic hydrogen. The Fortran implements the
        // rule documented in its own headers -- "a projectile with the same velocity as
        // an H atom has the same cross section" -- by evaluating the fit at E0/A, where
        // A is the projectile mass in u: D and H2 at E0/2, D2 and He at E0/4
        // (ion_surf_coll.f:132-134, :176-178, :220-222).
        //
        // Taking A straight from the particle mass reproduces every species the Fortran
        // enumerates (H->1, D->2, H2->2, D2->4, He->4) and extends consistently to H3+
        // (->3), so no species windows are needed for this part.
        double mass_ratio = particle->m() / MASS_U;
        if (mass_ratio < 1.0) mass_ratio = 1.0;  // never extrapolate below atomic hydrogen

        double E0_keV = energy * 1e-3;  // Convert to keV
        
        // Calculate backscattering probability
        double mu_i = 0.5;  // Free parameter
        double etaBACK = RN(E0_keV, mass_ratio);
        etaBACK = etaBACK / ((1.0 - mu_i) * cos(theta) + mu_i);
        if (etaBACK > 1.0) etaBACK = 1.0;
        
        // Calculate secondary electron emission probability
        double mu_e = 1.45;
        double etaSEC = eta_ion(E0_keV, mass_ratio);

        // Enhanced electron yield for molecular projectiles (Redbook Vol 1, C-10).
        // The Fortran applies this to H2, D2 and He but not to atomic H or D
        // (ion_surf_coll.f:47-49). D and H2 share A=2, so the mass alone cannot
        // separate them -- the species kind is genuinely required here.
        // H3+ is included as a consistent extension; the Fortran has no H3 species.
        // Note He is not part of the particle_kind enum, so it arrives as
        // PARTICLE_WRONG and does not pick up this factor.
        if (kind == PARTICLE_H2P || kind == PARTICLE_H20 || kind == PARTICLE_H3P) {
            etaSEC = etaSEC * 2.90 / 1.32;
        }

        etaSEC = etaSEC * exp(mu_e * (1.0 - cos(theta)));
        
        if (_debugprint) {
            CB_LOG_LOCK();
            logfile << "Ion surface collision @ E=" << E0_keV << " keV, "
                    << "theta=" << theta*180.0/M_PI << " deg" << endl;
            logfile << "  etaBACK=" << etaBACK << ", etaSEC=" << etaSEC << endl;
        }

        double Elost = energy;

        // Check for backscattering
        if (cbthread::rng().uniform() <= etaBACK && Elost >= 1.0) {
            // Backscattered ion
            double RE_val = RE(E0_keV, mass_ratio);
            double RN_val = RN(E0_keV, mass_ratio);
            double Eout = (RE_val / RN_val) * Elost;
            
            if (Eout >= 0.0 && Eout <= Elost) {
                
                // Create backscattered heavy particle. Same mass as the incident particle,
                // but the charge state is resampled: backscattered heavies emerge mostly
                // neutral (Fubiani et al. 2008, Sec. II C). Keeping the incident charge
                // would be qualitatively wrong -- a neutral is field-blind and flies
                // ballistically, while a positive one is accelerated back toward the source.
                //
                // add_particle() expects the charge in units of e and the mass in u
                // (particledatabaseimp.hpp:641 applies CHARGE_E*q and MASS_U*m), so the SI
                // value returned by m() has to be converted.
                double q_back = sampleBackscatteredChargeState(Eout / energy);

                if (_debugprint) {
                    CB_LOG_LOCK();
                    logfile << "  Backscattered heavy charge state: " << q_back << " e"
                            << " (was " << particle->q()/CHARGE_E << " e)" << endl;
                }

                // Debit Elost ONLY if the particle was really created.
                if (createSecondaryParticle(particle, loc, normal, tang1, tang2, Eout,
                                            q_back, particle->m()/MASS_U,
                                            particle->gen() + SURFACE_GENERATION_OFFSET)) {
                    Elost -= Eout;

                    if (_debugprint) {
                        CB_LOG_LOCK();
                        logfile << "  Ion backscattered with E=" << Eout*1e-3 << " keV" << endl;
                    }
                } else if (_debugprint) {
                    CB_LOG_LOCK();
                    logfile << "  Backscattered ion not created; Elost unchanged" << endl;
                }
            }
        }
        
        // Generate secondary electrons (Fortran cap for ion impacts: ion_surf_coll.f:101)
        generateSecondaryElectrons(particle, loc, normal, tang1, tang2, etaSEC, Elost,
                                   _max_secondary_electrons);
    }

    /**
     * @brief Generate secondary electrons after surface impact
     *
     * Two equivalent samplers, selected by _secondary_electron_bundle:
     *
     *   bundle <= 0  Analog. One unit-weight macroparticle per emitted electron, exactly as
     *                the Fortran does. Up to lsec_max (26 for ion impacts) macroparticles
     *                from a single impact.
     *
     *   bundle >= 1  Bundled. The same yield is carried by at most `bundle` macroparticles,
     *                each weighted so that the represented electron count is unchanged.
     *
     * WHY BUNDLING IS MEAN-PRESERVING
     *
     * The analog loop emits a first electron with probability etaSEC and then, for
     * etaSEC > 1, keeps testing eta_remaining -= 1 against a fresh deviate. That is a
     * Bernoulli decomposition of etaSEC: for etaSEC = 3.4 it emits 3 electrons with
     * certainty plus a fourth with probability 0.4, so the expected count is exactly
     * etaSEC. It is additionally truncated by the energy budget (Elost >= 10 eV per
     * electron) and by lsec_max. Applying those same two caps to the yield up front,
     *
     *     eta_eff = min( etaSEC, Elost/E_sec, lsec_max )
     *
     * and then emitting n = min(ceil(eta_eff), bundle) macroparticles of weight
     * eta_eff/n gives the same expected electron count, with zero sampling variance in
     * the count rather than Bernoulli variance. Nothing about the emission physics --
     * energy, angular law, charge -- changes; only how many macroparticles carry it.
     *
     * The eta_eff <= 1 case is left analog on purpose. Bundling there would emit a
     * fractional-weight macroparticle on EVERY impact, where the analog sampler emits one
     * on only a fraction of them, so it would increase the particle count rather than
     * reduce it. Bundling only pays off in the high-yield tail it was written for.
     *
     * @param particle Incident particle
     * @param loc Impact location
     * @param normal Surface normal
     * @param tang1 First tangent
     * @param tang2 Second tangent
     * @param etaSEC Secondary electron yield
     * @param Elost Energy available for secondary emission (in eV)
     * @param lsec_max Per-impact cap on the represented electron count
     */
    void generateSecondaryElectrons(ParticleBase *particle, const Vec3D& loc,
                                     const Vec3D& normal, const Vec3D& tang1, const Vec3D& tang2,
                                     double etaSEC, double& Elost, int lsec_max) {
        if (Elost < _secondary_electron_energy) return;
        if (!std::isfinite(etaSEC) || etaSEC <= 0.0) return;

        if (_secondary_electron_bundle < 1) {
            generateSecondaryElectronsAnalog(particle, loc, normal, tang1, tang2,
                                             etaSEC, Elost, lsec_max);
            return;
        }

        // Yield actually affordable at this impact, after the energy budget and the
        // per-impact cap the analog sampler would have applied anyway.
        double eta_eff = etaSEC;
        const double budget = Elost / _secondary_electron_energy;
        if (eta_eff > budget) eta_eff = budget;
        if (eta_eff > static_cast<double>(lsec_max)) {
            eta_eff = static_cast<double>(lsec_max);
            warnOnce(_electron_limit_warned,
                     "Maximum secondary electron count reached; suppressing repeated warnings");
        }
        if (eta_eff <= 0.0) return;

        // Bundle only where it actually reduces the macroparticle count. The analog
        // sampler emits ceil(eta_eff) macroparticles at most, so when that already fits
        // inside the bundle there is nothing to gain -- and rounding UP to n = ceil(eta)
        // would emit slightly MORE macroparticles than the analog sampler does on average
        // (eta = 3.4 would cost 4 instead of 3.4). Deferring to the analog sampler there
        // keeps the low- and mid-yield range bit-for-bit unchanged and confines bundling to
        // the high-yield tail it was written for.
        const int analog_macroparticles = static_cast<int>(ceil(eta_eff));
        if (analog_macroparticles <= _secondary_electron_bundle) {
            generateSecondaryElectronsAnalog(particle, loc, normal, tang1, tang2,
                                             etaSEC, Elost, lsec_max);
            return;
        }

        const int n = _secondary_electron_bundle;
        const double weight = eta_eff / static_cast<double>(n);

        const double energy_per_macroparticle = weight * _secondary_electron_energy;

        for (int k = 0; k < n; k++) {
            if (Elost < energy_per_macroparticle) break;

            // Debit Elost ONLY when the particle is really created (createSecondaryParticle
            // can refuse on generation depth, particle-count cap or invalid velocity) --
            // otherwise the energy vanishes from the books.
            if (!createSecondaryParticle(particle, loc, normal, tang1, tang2,
                                         _secondary_electron_energy,
                                         -1.0, MASS_E/MASS_U,
                                         particle->gen() + SURFACE_GENERATION_OFFSET,
                                         weight)) {
                if (_debugprint) {
                    CB_LOG_LOCK();
                    logfile << "  Secondary electron " << (k+1) << "/" << n
                            << " not created; Elost unchanged, stopping" << endl;
                }
                break;
            }

            Elost -= energy_per_macroparticle;

            if (_debugprint) {
                CB_LOG_LOCK();
                logfile << "  Secondary electron macroparticle " << (k+1) << "/" << n
                        << ", E=" << _secondary_electron_energy << " eV, weight=" << weight
                        << " (yield etaSEC=" << etaSEC << ", effective=" << eta_eff << ")" << endl;
            }
        }
    }

    /**
     * @brief Original one-macroparticle-per-electron sampler.
     *
     * Retained verbatim (bar thread-safety) so that setting
     * physics.surfaceCollisions.secondaryElectronBundle to 0 reproduces pre-bundling runs.
     */
    void generateSecondaryElectronsAnalog(ParticleBase *particle, const Vec3D& loc,
                                          const Vec3D& normal, const Vec3D& tang1, const Vec3D& tang2,
                                          double etaSEC, double& Elost, int lsec_max) {
        // Check for first secondary electron.
        if (cbthread::rng().uniform() <= etaSEC && Elost >= _secondary_electron_energy) {
            if (!createSecondaryParticle(particle, loc, normal, tang1, tang2,
                                         _secondary_electron_energy,
                                         -1.0, MASS_E/MASS_U,
                                         particle->gen() + SURFACE_GENERATION_OFFSET)) {
                if (_debugprint) {
                    CB_LOG_LOCK();
                    logfile << "  Secondary electron not created; Elost unchanged" << endl;
                }
                return;
            }

            Elost -= _secondary_electron_energy;

            if (_debugprint) {
                CB_LOG_LOCK();
                logfile << "  Secondary electron emitted, E=" << _secondary_electron_energy << " eV" << endl;
            }

            // Check for additional secondary electrons if yield > 1
            if (etaSEC > 1.0) {
                double eta_remaining = etaSEC - 1.0;
                int cnt_lsec = 1;

                while (eta_remaining > 0.0 && Elost >= _secondary_electron_energy &&
                       cnt_lsec < lsec_max) {
                    if (cbthread::rng().uniform() <= eta_remaining) {
                        if (!createSecondaryParticle(particle, loc, normal, tang1, tang2,
                                                     _secondary_electron_energy,
                                                     -1.0, MASS_E/MASS_U,
                                                     particle->gen() + SURFACE_GENERATION_OFFSET)) {
                            if (_debugprint) {
                                CB_LOG_LOCK();
                                logfile << "  Additional secondary electron not created; "
                                        << "Elost unchanged, stopping" << endl;
                            }
                            break;
                        }

                        Elost -= _secondary_electron_energy;
                        cnt_lsec++;

                        if (_debugprint) {
                            CB_LOG_LOCK();
                            logfile << "  Additional secondary electron #" << cnt_lsec
                                    << ", E=" << _secondary_electron_energy << " eV" << endl;
                        }
                    }

                    eta_remaining -= 1.0;
                }

                if (cnt_lsec >= lsec_max) {
                    warnOnce(_electron_limit_warned,
                             "Maximum secondary electron count reached; suppressing repeated warnings");
                }
            }
        }
    }

    /**
     * @brief Create a secondary particle with proper angular distribution
     * Based on make_sec_wall_part Fortran routine
     * @param particle Incident particle (for current/charge info)
     * @param loc Impact location
     * @param normal Surface normal
     * @param tang1 First tangent
     * @param tang2 Second tangent
     * @param Ek Secondary particle energy in eV
     * @param charge Secondary particle charge
     * @param mass Secondary particle mass (in atomic mass units)
     * @param generation Generation number
     * @param weight How many real secondaries this one macroparticle stands for, relative
     *        to the incident macroparticle. 1.0 reproduces the original behaviour.
     */
    bool createSecondaryParticle(ParticleBase *particle, const Vec3D& loc,
                                 const Vec3D& normal, const Vec3D& tang1, const Vec3D& tang2,
                                 double Ek, double charge, double mass, int generation,
                                 double weight = 1.0) {
        // Generate random emission angles (make_sec_wall_part.f:41-45).
        // theta = asin(rnd), phi = 2*pi*rnd. Two independent deviates: the previous code
        // read two numbers out of a one-dimensional MTRandom, so phi was uninitialized
        // stack memory rather than a uniform azimuth.
        const double theta = asin(cbthread::rng().uniform());
        const double phi = 2.0 * M_PI * cbthread::rng().uniform();
        return emitSecondaryParticle(particle, loc, normal, tang1, tang2, Ek, charge, mass,
                                     generation, theta, phi, weight);
    }

    /**
     * @brief Create a secondary particle emitted in a given direction
     *
     * Same as createSecondaryParticle but with the emission direction supplied by the
     * caller, for cases where the direction is correlated with the sampled energy and
     * therefore cannot be drawn independently here (see sampleBackscatteredElectron).
     *
     * @param theta Emission angle from the surface normal (radians)
     * @param phi Azimuthal emission angle (radians)
     * @param weight How many real secondaries this macroparticle stands for, relative to
     *        the incident macroparticle. Scales the child's IQ and, with it, every
     *        current, power and space-charge contribution the child makes downstream, plus
     *        the emitting surface's ledger debit. 1.0 reproduces the original behaviour.
     * @return true only if the particle was actually added to the database. Callers MUST
     *         gate the Elost debit on this: debiting energy for a particle that was never
     *         created makes it vanish from the books, which the Phase-2 net accounting
     *         would then report as a surface debited with no matching credit.
     */
    bool emitSecondaryParticle(ParticleBase *particle, const Vec3D& loc,
                               const Vec3D& normal, const Vec3D& tang1, const Vec3D& tang2,
                               double Ek, double charge, double mass, int generation,
                               double theta, double phi, double weight = 1.0) {
        if (!std::isfinite(weight) || weight <= 0.0) {
            return false;
        }
        // Check generation depth using modulo 100 to get effective generation
        // This works for both volume-generated (gen+1) and surface-generated (gen+101) particles
        int effective_gen = generation % 100;
        if (effective_gen >= MAX_GENERATION_DEPTH) {
            if (_debugprint) {
                CB_LOG_LOCK();
                logfile << "WARNING: Maximum generation depth reached (effective gen="
                        << effective_gen << ", full gen=" << generation << ")" << endl;
            }
            return false;
        }

        // Check particle count
        if (secondaryBudgetExhausted()) {
            warnOnce(_secondary_limit_warned,
                     "Secondary particle pass limit reached; stopping further surface "
                     "secondary creation");
            return false;
        }

        // Calculate velocity from energy (non-relativistic approximation for low energies)
        // For high energies, use relativistic formula
        //
        // Ek is in eV, so the rest energy must be in eV as well: mass*MASS_U*SPEED_C2
        // is in joules and has to be divided by CHARGE_E. Matches the Fortran, which
        // carries mc2 in eV throughout (make_sec_wall_part.f:33-38, mec2/mpc2).
        double mc2 = mass * MASS_U * SPEED_C2 / CHARGE_E;  // Rest energy in eV
        double gamma = 1.0 + Ek / mc2;
        double Vtot;
        
        if (gamma > 1.01) {
            // Relativistic case
            Vtot = sqrt(gamma * gamma - 1.0) / gamma * SPEED_C;
        } else {
            // Non-relativistic case
            Vtot = sqrt(2.0 * Ek * CHARGE_E / (mass * MASS_U));
        }
        
        // Validate velocity
        if (!std::isfinite(Vtot) || Vtot < 0.0 || Vtot > MAX_VELOCITY_SCALE) {
            if (debug) {
                CB_LOG_LOCK();
                logfile << "WARNING: Invalid secondary particle velocity: " << Vtot << endl;
            }
            return false;
        }
        
        // Calculate velocity components in local coordinate system (u, t, n)
        double V11 = sin(theta) * cos(phi);
        double V22 = sin(theta) * sin(phi);
        double V33 = cos(theta);
        
        // Transform to global coordinates.
        //
        // The local frame is (u, t, n) = (tang1, tang2, normal), so the cos(theta)
        // component must go along the surface normal (make_sec_wall_part.f:58-60).
        // Since theta = asin(rnd) is in [0, pi/2], cos(theta) >= 0 and every
        // secondary is emitted into the outward hemisphere. Mapping sin(theta)*cos(phi)
        // onto the normal instead would send half of them into the solid, because
        // cos(phi) is negative over half of the sampled phi range.
        Vec3D vel_global = V11 * tang1 + V22 * tang2 + V33 * normal;
        vel_global *= Vtot;
        
        // Adjust location off the surface, verifying the offset actually lands in vacuum.
        //
        // A fixed 0.01*h offset (h = mesh spacing) is only ~1% of a cell, and
        // Geometry::surface_normal() is a local, finite-difference estimate: at a thin
        // grid web or a concave/re-entrant corner (routine for a multi-aperture
        // accelerator grid), that offset can still land inside the emitting solid or a
        // neighbouring one. A particle whose very first position is invalid there is a
        // DT-INDEPENDENT failure for the trajectory integrator -- ParticleP3D::get_derivatives
        // (particles.cpp) returns IBSIMU_DERIV_ERROR for every step-size retry from the same
        // starting point, so gsl_odeiv_evolve_apply's step-halving loop in particleiterator.hpp
        // never succeeds and dt underflows to 0.0, throwing exactly the observed
        // Error("too small step size") -- deterministically, on whichever iteration first
        // creates a surface secondary there (only the last one, since THCallback_surf_EAMCC
        // is only attached on the final iteration; see ManageSimulation_New.cpp).
        //
        // Escalate the offset until Geometry::inside() confirms vacuum (return value 0;
        // see geometry.hpp), and give up rather than launch a doomed particle if it never
        // clears the solid within a reasonable multiple of the cell size.
        static const double kOffsetFactors[] = { 0.01, 0.05, 0.25, 1.0, 2.0, 4.0 };
        const double mesh_h = _geom.h();
        Vec3D loc_adj = loc + kOffsetFactors[0] * mesh_h * normal;
        bool placed_in_vacuum = ( _geom.inside( loc_adj ) == 0 );
        if( !placed_in_vacuum ) {
            for( size_t fi = 1; fi < sizeof(kOffsetFactors)/sizeof(kOffsetFactors[0]); fi++ ) {
                loc_adj = loc + kOffsetFactors[fi] * mesh_h * normal;
                if( _geom.inside( loc_adj ) == 0 ) {
                    placed_in_vacuum = true;
                    if( _debugprint ) {
                        CB_LOG_LOCK();
                        logfile << "  Secondary offset escalated to " << kOffsetFactors[fi]
                                << "*h to clear solid" << endl;
                    }
                    break;
                }
            }
        }
        if( !placed_in_vacuum ) {
            CB_LOG_LOCK();
            logfile << "WARNING: Secondary particle start point still inside a solid after "
                     << "offset escalation (loc=" << loc << ", normal=" << normal
                     << "); not created" << endl;
            return false;
        }

        if (_debugprint) {
            CB_LOG_LOCK();
            double emission_angle = acos(V33) * 180.0 / M_PI;
            logfile << "  Secondary particle: E=" << Ek << " eV, "
                    << "angle=" << emission_angle << " deg, "
                    << "charge=" << charge << ", mass=" << mass
                    << ", weight=" << weight << endl;
        }

        // Create ParticleP3D: (t, x, vx, y, vy, z, vz) - velocity components directly
        ParticleP3D sec_part(0.0,  // time
                            loc_adj[0], vel_global[0],
                            loc_adj[1], vel_global[1],
                            loc_adj[2], vel_global[2]);
        
        // Secondary macroparticle current.
        //
        // IQ is sign-carrying in this codebase: negative for negative and neutral species,
        // positive for positive ions (which is why the power column comes out positive,
        // wdata = IQ * Ekin/q). The gas-phase callback enforces this via
        // thcallback_detail::child_particle_current (THCallback.h:148-153), flipping the
        // sign for positive-ion products. Passing the parent's IQ through unchanged, as
        // this did, is only correct when the parent is negative: a secondary electron off a
        // positive-ion impact would inherit IQ > 0 while carrying q < 0, and then register
        // NEGATIVE power and a wrong-sign current on the grid it lands on.
        //
        // Setting the sign from the child's own charge reproduces child_particle_current
        // for negative parents and stays correct for positive ones.
        //
        // `weight` scales the macroparticle: a weight-w child stands for w times as many
        // real particles per second as the incident macroparticle did. Because everything
        // downstream -- space charge (particleiterator.hpp deposits particle.IQ()), the
        // grid current column, and the power column wdata = IQ*Ekin/q -- is built from IQ
        // and never from a particle count, scaling IQ is all that is needed for a bundled
        // secondary to carry exactly the physics of the electrons it replaces.
        double child_IQ = (charge > 0.0 ? 1.0 : -1.0) * fabs(particle->IQ()) * weight;

        // Add to particle database. Catch ibsimu's own Error (e.g. ErrorNoMem on
        // allocation failure) as well as std::exception: Error does NOT derive from
        // std::exception (error.hpp), so a plain std::exception handler cannot see it
        // and it would otherwise unwind uncaught all the way out of the trajectory
        // scheduler -> std::terminate() -> abort() (SIGABRT), instead of degrading to
        // "this one secondary was not created" as intended here.
        //
        // add_particle() is safe to call from a worker thread: it takes the scheduler's
        // problem-vector mutex (particledatabaseimp.hpp:626-634), which is exactly the
        // protocol Scheduler documents for adding problems from any thread.
        try {
            _pdb->add_particle(child_IQ, charge, mass, generation, sec_part);
            _secondaries_created.fetch_add(1, std::memory_order_relaxed);
        } catch (Error& e) {
            CB_LOG_LOCK();
            logfile << "ERROR: Failed to add secondary particle: " << e.get_error_message() << endl;
            return false;
        } catch (const std::exception& e) {
            CB_LOG_LOCK();
            logfile << "ERROR: Failed to add secondary particle: " << e.what() << endl;
            return false;
        }

        // Debit the emitting surface. Only reached once the particle really exists, so the
        // ledger can never claim energy left a surface that nothing carried away.
        //
        // The rate is taken from the PARENT, never from the child: a neutral child has
        // q == 0, so |IQ_child|/|q_child| is a division by zero. The parent's own rate is
        // the right one anyway -- one incident macroparticle emits one secondary
        // macroparticle, so they represent the same number of real particles per second.
        // (A neutral parent is handled with the same CHARGE_E fallback PowerStruct::add
        // uses, so that the rate stays a particle rate rather than becoming infinite.)
        //
        // With bundling that one-for-one correspondence no longer holds, so the parent rate
        // is scaled by `weight` -- the same factor applied to child_IQ above. Without this
        // the surface would be debited for one electron while the child carried away w of
        // them, and the net power column would no longer balance.
        if (_ledger != NULL) {
            // Unit care: particle->q()/m() are SI (coulombs, kilograms), whereas this
            // function's `charge` and `mass` parameters are in elementary charges and
            // atomic mass units -- that is what ParticleDataBase3D::add_particle expects.
            // identify_particle_species compares against absolute kg/C thresholds, so the
            // child's values have to be converted before either use.
            double q_parent  = particle->q();                    // [C]
            double q_scale   = (q_parent == 0.0) ? CHARGE_E : fabs(q_parent);
            double rate      = weight * fabs(particle->IQ()) / q_scale;  // particles per second

            double q_child_C = charge * CHARGE_E;                // signed [C]; 0 for neutrals
            double m_child_kg = mass * MASS_U;
            double Ek_joules = Ek * CHARGE_E;                    // Ek arrives in eV

            // Negative: energy and charge LEAVE this surface. current_A is formed from the
            // child's signed physical charge, so a backscattered neutral removes no charge
            // while still removing energy -- the same rule Phase 1 applied to arrivals.
            _ledger->record(thcallback_surf_detail::current_hit_solid(), loc,
                            static_cast<int>(identify_particle_species(m_child_kg, q_child_C,
                                                                       _mass)),
                            generation,
                            -(rate * Ek_joules),
                            charge == 0.0 ? 0.0 : -(rate * q_child_C));
        }

        return true;
    }

public:
    /**
     * @brief Constructor
     * @param geom Geometry object for surface normal calculation
     * @param pdb Particle database for adding secondary particles
     * @param mass Reference mass for particle type identification
     * @param debugprint Enable debug output
     * @param minimum_z Impacts below this z produce no secondaries
     * @param secondary_electron_bundle Max macroparticles per impact used to carry the
     *        secondary-electron yield; 0 selects the original per-electron sampler.
     *        See generateSecondaryElectrons.
     * @param max_secondary_electrons Per-impact cap on secondary electrons from an ION or
     *        NEUTRAL impact (lsec_max). Electron impacts keep MAX_SECONDARY_ELECTRONS_E.
     * @param secondary_electron_energy Birth energy of a true secondary electron [eV].
     */
    THCallback_surf_EAMCC(Geometry &geom,
                          ParticleDataBase3D* pdb,
                          double mass = 1.0,
                          bool debugprint = false,
                          double minimum_z = 7.0e-3,
                          int secondary_electron_bundle = DEFAULT_SECONDARY_ELECTRON_BUNDLE,
                          int max_secondary_electrons = DEFAULT_MAX_SECONDARY_ELECTRONS,
                          double secondary_electron_energy = DEFAULT_SECONDARY_ELECTRON_ENERGY)
        : _geom(geom), _pdb(pdb), _mass(mass), _minimum_z(minimum_z), _debugprint(debugprint),
                    _ledger(NULL), _secondary_limit_warned(false), _electron_limit_warned(false),
                    _initial_particle_count(pdb ? pdb->size() : 0), _secondaries_created(0),
                    _secondary_electron_bundle(secondary_electron_bundle),
                    _max_secondary_electrons(max_secondary_electrons),
                    _secondary_electron_energy(secondary_electron_energy) {
        // Last line of defence. A cap below 1 would switch true secondary emission off
        // entirely, and a birth energy of 0 would create electrons at rest that never
        // leave the surface and stall the trajectory integrator. Neither is a plausible
        // intent, so fall back rather than honour them.
        if (_max_secondary_electrons < 1) {
            CB_LOG_LOCK();
            logfile << "WARNING: maximumSecondaryElectrons=" << _max_secondary_electrons
                    << " is not usable; falling back to " << DEFAULT_MAX_SECONDARY_ELECTRONS << endl;
            _max_secondary_electrons = DEFAULT_MAX_SECONDARY_ELECTRONS;
        }
        if (!(_secondary_electron_energy > 0.0) || !std::isfinite(_secondary_electron_energy)) {
            CB_LOG_LOCK();
            logfile << "WARNING: secondaryElectronEnergyEV=" << _secondary_electron_energy
                    << " is not usable; falling back to " << DEFAULT_SECONDARY_ELECTRON_ENERGY
                    << " eV" << endl;
            _secondary_electron_energy = DEFAULT_SECONDARY_ELECTRON_ENERGY;
        }

        {
            CB_LOG_LOCK();
            logfile << "THCallback_surf_EAMCC: max secondary electrons per ion impact="
                    << _max_secondary_electrons
                    << " (electron impacts fixed at " << MAX_SECONDARY_ELECTRONS_E << ")"
                    << ", secondary electron energy=" << _secondary_electron_energy << " eV"
                    << ", bundle="
                    << (_secondary_electron_bundle < 1
                            ? std::string("off, one macroparticle per electron")
                            : std::to_string(_secondary_electron_bundle))
                    << endl;
        }
    }

    /*! \brief Secondary macroparticles created by this callback in the current pass. */
    size_t secondariesCreated() const {
        return _secondaries_created.load(std::memory_order_relaxed);
    }

    /**
     * @brief Attach the signed surface energy/charge ledger.
     *
     * Optional: with no ledger attached the callback behaves exactly as before and only
     * the gross (arrival-side) power and current columns are available. The ledger is not
     * owned, and must outlive the tracking pass.
     */
    void setEventLedger(SurfaceEventLedger* ledger) { _ledger = ledger; }

    virtual ~THCallback_surf_EAMCC() {}

    /**
     * @brief Main callback function called when particle trajectory ends
     * @param particle Particle that hit the surface
     * @param pdb Particle database
     */
    virtual void operator()(ParticleBase *particle, class ParticleDataBase *pdb) {
        Particle3D *p3d = (Particle3D *)(particle);
        Vec3D loc = p3d->location();
        Vec3D vel = p3d->velocity();
        
        // Validate position
        if (!validatePosition(loc)) {
            return;
        }
        
        // Only generate secondary particles beyond the configured impact plane.
        if (loc[2] <= _minimum_z) {
            return;
        }
        
        // Limit generation depth to prevent infinite chains
        // Check effective generation (works for both volume and surface-generated particles)
        int generation = particle->gen();
        int effective_gen = generation % 100;
        if (effective_gen >= MAX_GENERATION_DEPTH) {
            if (_debugprint) {
                CB_LOG_LOCK();
                logfile << "WARNING: Maximum generation depth reached (gen=" << generation
                        << ", effective=" << effective_gen << "), stopping surface secondary generation" << endl;
            }
            return;
        }

        // Check particle count to prevent memory overflow
        if (secondaryBudgetExhausted()) {
            warnOnce(_secondary_limit_warned,
                     "Maximum particle count reached; stopping surface secondary generation");
            return;
        }

        // Calculate energy
        double vel_mag = vel.ssqr();
        if (vel_mag < 0.0 || !std::isfinite(vel_mag)) {
            if (debug) {
                CB_LOG_LOCK();
                logfile << "WARNING: Invalid velocity magnitude" << endl;
            }
            return;
        }
        vel_mag = sqrt(vel_mag);
        double energy = particle->m() * SPEED_C2 * (1.0 / sqrt(1.0 - (vel_mag * vel_mag) / (SPEED_C2)) - 1.0) / CHARGE_E;  // in eV

        if (!std::isfinite(energy) || energy < 0.0) {
            if (debug) {
                CB_LOG_LOCK();
                logfile << "WARNING: Invalid particle energy: " << energy << " eV" << endl;
            }
            return;
        }

        // Get surface normal
        Vec3D normal = _geom.surface_normal(loc);
        if (normal.ssqr() < 1e-10) {
            if (debug) {
                CB_LOG_LOCK();
                logfile << "WARNING: Could not determine surface normal" << endl;
            }
            return;
        }
        normal.normalize();

        // Only impacts on a real solid (a grid) produce surface secondaries.
        //
        // IBSimu calls this callback at the end of EVERY trajectory, whatever the reason
        // (particleiterator.hpp:1475 -- after PARTICLE_COLL, PARTICLE_OUT, PARTICLE_NSTP or
        // PARTICLE_TIME alike), and Geometry::inside() returns boundary indices 1..6 for
        // out-of-mesh coordinates (geometry.cpp:343-354). A particle that merely left the
        // domain therefore still gets a valid surface normal, so without this check the
        // entire transmitted beam would emit secondaries as if it had struck a grid.
        // Defined solids are indexed from 7 upward.
        uint32_t hit_solid = _geom.inside(loc);
        if (hit_solid < 7) {
            hit_solid = _geom.inside(loc - 0.1 * _geom.h() * normal);
        }
        if (hit_solid < 7) {
            if (_debugprint) {
                CB_LOG_LOCK();
                logfile << "  Trajectory ended without hitting a defined solid (index="
                        << hit_solid << "), no surface secondaries" << endl;
            }
            return;
        }

        // Publish the resolved solid for the emission records made further down this call.
        // Thread-local: every worker thread runs this same callback object concurrently.
        thcallback_surf_detail::current_hit_solid() = hit_solid;

        // Calculate impact angle
        Vec3D vel_normalized = vel;
        vel_normalized.normalize();
        double cos_theta = -(vel_normalized[0]*normal[0] + vel_normalized[1]*normal[1] + vel_normalized[2]*normal[2]);  // Negative because particle is coming in
        if (cos_theta > 1.0) cos_theta = 1.0;
        double theta = acos(cos_theta);

        // Unphysical impact angle (velocity and outward normal disagree). The Fortran does
        // not discard the particle: it declares the angle "good" and substitutes 45 deg
        // (checkcollisions.f:128-221, counting these in cnt_badth). Clamping to pi/2
        // instead would pick grazing incidence, which maximizes both etaBACK and etaSEC
        // through the exp(CT*(1-cos theta)) factors and roughly doubles the yield on every
        // such impact.
        if (cos_theta < 0.0) {
            theta = M_PI / 4.0;
            if (_debugprint) {
                CB_LOG_LOCK();
                logfile << "  Bad impact angle (cos_theta=" << cos_theta
                        << "), substituting 45 deg" << endl;
            }
        }
        
        // Create tangent vectors
        Vec3D tang1 = normal.arb_perpendicular();
        Vec3D tang2 = cross(normal, tang1);
        tang1.normalize();
        tang2.normalize();
        
        // Determine particle type using the same species classification as the rest of
        // the code (funct.cpp identify_particle_species, cf. THCallback.h:201). _mass is
        // the working-species ion mass in u, as passed by ManageSimulation_New.cpp.
        particle_kind kind = identify_particle_species(particle->m(), particle->q(), _mass);

        if (kind == PARTICLE_E) {
            handleElectronSurfaceCollision(particle, loc, vel, energy, theta, normal, tang1, tang2);
        } else if (particle->m() > 0.5 * MASS_U) {
            // Everything heavier than an electron goes through the ion model. Projectiles
            // outside the H/D family windows come back as PARTICLE_WRONG (He while running
            // hydrogen, for instance); they are still handled, because the RN/RE/eta_ion
            // scaling depends only on the mass and only the molecular yield boost needs
            // the kind.
            if (kind == PARTICLE_WRONG && _debugprint) {
                CB_LOG_LOCK();
                logfile << "  Unclassified projectile (m=" << particle->m()/MASS_U
                        << " u, q=" << particle->q()/CHARGE_E
                        << " e), using ion surface model" << endl;
            }
            handleIonSurfaceCollision(particle, loc, vel, energy, theta, normal, tang1, tang2, kind);
        } else if (_debugprint) {
            CB_LOG_LOCK();
            logfile << "  Projectile is neither an electron nor heavy enough for the ion "
                    << "model (m=" << particle->m()/MASS_U << " u), no surface secondaries" << endl;
        }
        
        // Mark original particle as absorbed (or could be backscattered, handled above)
        // The particle status will be set by the trajectory handler
    }
};

#endif /* THCALLBACK_SURF_EAMCC_H_ */

