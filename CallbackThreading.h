#ifndef CALLBACKTHREADING_H_
#define CALLBACKTHREADING_H_

/*
 * CallbackThreading.h
 *
 * Per-thread support for the trajectory callbacks, so that ibsimu can run
 * pdb->iterate_trajectories() on more than one thread while stripping and/or surface
 * collisions are active.
 *
 * WHY THIS EXISTS
 *
 * ibsimu builds one ParticleIterator per thread but hands every one of them the SAME
 * callback object (particledatabaseimp.hpp:687-700). Anything mutable inside a callback is
 * therefore shared across worker threads. The library itself is already prepared for this:
 * Scheduler documents that "all threads can add problems" as long as the problem vector is
 * protected, and ParticleDataBasePPImp::add_particle already takes that lock
 * (particledatabaseimp.hpp:626-634). What was missing was thread-safety inside our own
 * callbacks, which is what this header provides.
 *
 * RANDOM NUMBERS
 *
 * ibsimu's MTRandom is mutex-protected internally, so sharing one across threads is safe
 * but serializes every draw -- and the gas callback draws on every integration step of
 * every particle, so that mutex would become the bottleneck the threading was meant to
 * remove. MTRandom also offers no way to seed an instance (MTRandom(N) takes the number of
 * DIMENSIONS, not a seed, and gsl_rng_alloc always starts from the default seed), so
 * handing each thread its own MTRandom would give every thread the IDENTICAL stream --
 * perfectly correlated sampling, which is statistically much worse than the contention.
 *
 * Hence a small per-thread generator seeded through std::seed_seq from a global base seed
 * plus a distinct per-thread stream index. Each thread gets an independent, well-separated
 * MT19937-64 stream and never touches a shared lock.
 *
 * REPRODUCIBILITY: with threadCount == 1 there is exactly one stream and runs are
 * bit-reproducible. With N > 1 the assignment of particles to streams depends on OS thread
 * scheduling, so runs are statistically equivalent but not bit-identical -- the usual
 * property of parallel Monte Carlo. Set run.threadCount to 1 when a byte-exact rerun is
 * needed.
 */

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <mutex>
#include <random>

namespace cbthread {

/*! \brief Uniform (0,1) generator owned by exactly one thread.
 *
 *  Values are drawn from the half-open interval [0,1). Callers that feed the value into a
 *  logarithm must guard against 0 exactly as they did with MTRandom, whose
 *  gsl_rng_uniform() has the same range.
 */
class ThreadRng {
public:
    explicit ThreadRng(std::seed_seq& seq) : _engine(seq), _dist(0.0, 1.0) {}

    /*! \brief One uniform deviate in [0,1). */
    double uniform() { return _dist(_engine); }

    /*! \brief Fill \a n uniform deviates into \a x.
     *
     *  Unlike Random::get(), which silently fills only as many values as the generator was
     *  constructed with, this takes the count explicitly. The old code asked a
     *  one-dimensional MTRandom for two numbers in several places and used the second,
     *  uninitialized, one; an explicit count makes that class of mistake impossible.
     */
    void fill(double* x, std::size_t n) {
        for (std::size_t i = 0; i < n; ++i) {
            x[i] = uniform();
        }
    }

private:
    ThreadRng(const ThreadRng&);
    ThreadRng& operator=(const ThreadRng&);

    std::mt19937_64 _engine;
    std::uniform_real_distribution<double> _dist;
};

namespace detail {

inline std::atomic<std::uint64_t>& base_seed_storage() {
    // Arbitrary but fixed, so that a default-configured run is repeatable at threadCount 1.
    static std::atomic<std::uint64_t> seed(20240917ULL);
    return seed;
}

inline std::atomic<std::uint64_t>& stream_counter() {
    static std::atomic<std::uint64_t> counter(0);
    return counter;
}

} // namespace detail

/*! \brief Set the base seed for every thread stream created from now on.
 *
 *  Must be called before tracking starts; threads that have already created their
 *  generator keep the stream they have.
 */
inline void set_base_seed(std::uint64_t seed) {
    detail::base_seed_storage().store(seed, std::memory_order_relaxed);
}

/*! \brief This thread's generator, created on first use. */
inline ThreadRng& rng() {
    static thread_local ThreadRng* instance = NULL;
    if (instance == NULL) {
        const std::uint64_t base = detail::base_seed_storage().load(std::memory_order_relaxed);
        const std::uint64_t stream = detail::stream_counter().fetch_add(1, std::memory_order_relaxed);
        std::seed_seq seq{ static_cast<std::uint32_t>(base & 0xffffffffULL),
                           static_cast<std::uint32_t>(base >> 32),
                           static_cast<std::uint32_t>(stream & 0xffffffffULL),
                           static_cast<std::uint32_t>(stream >> 32) };
        // Leaked deliberately: thread_local destruction order relative to the trajectory
        // scheduler's thread teardown is not worth the risk for a few hundred bytes.
        instance = new ThreadRng(seq);
    }
    return *instance;
}

/*! \brief How many distinct thread streams have been created so far.
 *
 *  One per thread that has actually drawn a random number from inside a callback, so this
 *  reports how many worker threads really entered the Monte Carlo path -- the direct check
 *  that raising simulation.threadCount had the intended effect.
 */
inline std::uint64_t stream_count() {
    return detail::stream_counter().load(std::memory_order_relaxed);
}

/*! \brief Serializes writes to the global logfile from callback code.
 *
 *  std::ofstream gives no atomicity for a chain of operator<< calls, so multi-line warning
 *  blocks from two threads would interleave into unreadable output.
 */
inline std::mutex& log_mutex() {
    static std::mutex m;
    return m;
}

} // namespace cbthread

/*! \brief Guard one logfile statement.
 *
 *  Usage: CB_LOG_LOCK(); logfile << ... << endl;  -- the lock lives to the end of the
 *  enclosing scope, so wrap multi-statement blocks in braces.
 */
#define CB_LOG_LOCK() std::lock_guard<std::mutex> _cb_log_guard(cbthread::log_mutex())

#endif /* CALLBACKTHREADING_H_ */
