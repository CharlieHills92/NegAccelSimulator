      subroutine e_surf_coll(Ek,Elost,theta,iseed,u,t,n,
     &        x_g,y_g,z_g,Ysec)
c     ==============================================================
c     VERSION:         0.1
c     LAST MOD:       SEP/06
c     MOD AUTHOR:    G. Fubiani
c     COMMENTS:      Calculates electron collision with surface
c                    Theta corresponds to the impact angle between
c                    particle and wall.
c     NOTE:          Functions ETA0, ETA1 and have been
c                    written by H. P. L. de Esch
c                    Info. on the emission probability coeff. used
c                    may be found in the EFDA report "Electrons in SINGAP",
c                    H. P. L. De Esch, April 2006.
c     --------------------------------------------------------------
      implicit none
      include 'particle_info.h'
      include 'ode_arrays.h'
      integer iseed,ptype_info,option
      double precision Ysec(nmax_sec),u(3),t(3),n(3),x_g,y_g,z_g
      double precision Ek,theta,etaBACK,CTsec,CTback,CT,etaSEC,eta0, 
     &     eta1,rnd,ran2,Eout,Elost

c     Note: Ek is already in eV
c     Initialization of Elost. It represent the energy lost by particle 
c     after collision. For now it is set equal to Ek
      Elost=Ek
      etaBACK = eta0(Ek)
      CT = CTback(Ek,etaBACK)
      etaBACK=etaBACK*dexp( CT*( 1.d0 - dcos(theta) ) )
      
      etaSEC = eta1(Ek)
      CT = CTsec(Ek)
      
      etaSEC=etaSEC*dexp( CT*( 1.d0 - dcos(theta) ) )
       
      !print*,'x(mm)=',Ysec(4),'y(mm),',Ysec(5),'z(mm),',Ysec(6) 
    
      if(flag_info.eq.1)then 
         write(6,10)etaBACK,etaSEC,Ek*1.d-3
 10      format('Probabilities: backscatter=',f5.2,
     &          ', secondary=',f5.2,' @ Ek=',f8.2,' keV')
      endif

C
C     Backscattering. 
C
      rnd=ran2(iseed)
      if(rnd.le.etaBACK) then ! Backscattered

         ptype_info=1
         n_sec=n_sec+1
         flag_sec=1

         option=2
         if(option.eq.1) then  ! Back. energy spectra option #1

            call ebackdist_data(Eout,Ek,iseed) ! Direct interp. of exp. data

            call make_sec_wall_part(Ysec,Eout,iseed,u,t,n,x_g,y_g,z_g,
     &           ptype_info)

         else ! Back. energy spectra option #2

            call ebackdist_anal(Eout,Ek,theta,Ysec,iseed,u,t,n,x_g,y_g,
     &           z_g,ptype_info)  ! Use analytical fit of exp. data

         endif

         Elost=Elost-Eout ! Energy lost by incoming particle (eV)

         if(flag_info.eq.1)then         
            write(6,30) pname(ptype),Eout*1.d-3
 30         format(a2,' is backscattered. Eout=',f8.2,' keV')
         endif

         if(Eout.lt.0.d0.or.Eout.gt.Ek) then ! WARNING
            print*, 'Warning: output error from back. elec. dist.'
            flag_err=1
         endif
       
      else ! Not backscattered

      call checknetcharge(charge(1)) ! Cnt net charge absorbed by grid  
  
      endif

C
C     Secondary electrons. 
C    
      call generate_sec_e(etaSEC,Elost,Ysec,iseed,u,t,n,x_g,y_g,z_g,3)

      return
      end

      function ETA0(E0)
C     =================
C     BACKSCATTERING
C     Re-emission probability for an incident electron with energy E0 (eV).
C     Originally hacked from A. Simonin ELSTOP which has energies
C     below 500 eV and above 30 kV undefined (ie: set to zero!)
C
C     Changed 25 May 2005.
C     Redbook data (D-4) for Cu from 500 eV to 100 keV
C                   Scale using Mo to extrapolate down to 200 eV
C                   Scale using Ag and Au to extrapolate down to 100 eV
C                   PJ Ebert, phys rev 183(1969)422 fig. 18 for 1 MeV point
C                   AMD Assa'd, Scanning Microscopy 12(1998)1 checks OK.
C
C                   L. Wang, Phys. Rev. ST Accel. Beams 8(2005)094201
C                        gives data between 10 keV and 5 MeV
C
C     Some conflict with Darlington 1972.
C
      implicit none
      double precision    E0, ETA0
      double precision    E(18),eta(18)
      integer i

      data   E /                            0.,   100.,   200.,   
     >                    500.,  1000.,  2000.,  4000.,  7000.,
     >          10000., 20000., 40000., 70000.,100000.,
     >                 200000.,400000.,700000., 1.0e6 , 4.0e6  /
      data eta /                         0.0,    0.124,  0.161,  
     >                   0.262,  0.266,  0.282,  0.315,  0.322,
     >           0.323,  0.318,  0.309,  0.299,  0.285,
     >                   0.275,  0.260,  0.240,  0.210,  0.11  /

      eta0 = 0.

      do 100 i=2,18,1
      if( E0.gt.E(i-1) .and. E0.le.E(i) )then
         eta0=eta(i-1) + (E0-E(i-1))*(eta(i)-eta(i-1))/(E(i)-E(i-1))
      end if
100   continue

      return
      end


      function CTback(E0,eta0)
c     ==============================================================
c     VERSION:         0.1
c     LAST MOD:      FEB/07
c     MOD AUTHOR:    G. Fubiani
c     COMMENTS:      CT factor for backscattered electrons 
c                    with an incident electron E0 in keV.
c     NOTE:          Replaces data from A. Simonin ELSTOP Functions
c                    New Expression for CTback comes from 
c                    P-F Staub et al. [J. of Phys. D 27, 1533 (1994)], 
c                    Eqs (4) and (5).
c     --------------------------------------------------------------
      implicit none
      double precision E0, CTback,kappa,eta0

      kappa = 1.d0 - dexp( -1.83d0*(E0*1.d-3)**0.25d0 )
      CTback = kappa*log(1.d0/eta0) 

      return
      end


      function ETA1(Ep)
C     =================
C     True Secondary electrons
C     Re-emission probability for an incident electron with energy E0 (eV).
C     Toshikawa J. Phys. D (Appl. Phys.)  6(1973)1369   0.5 - 10 keV
C     Mc Allister ??? (1921)122   5 eV - 750 eV
C     Henrist, Hilleret et al: CERN LHC Project Report 472 (24 June 2002)
C
C     Adopted the fit for FULLY CONDITIONED copper surfaces from CERN.
C     need 2-10 mA sec / mm2 electrons to achieve this state. 
C
C     No redbook data on this.

      implicit none
C- input/output
      double precision    Ep, eta1
C- for deltaS
      double precision    s, Emax, DeltaMax
      double precision    EpEm
C- for f
      double precision    A0, A1, A2, E0, LnEpE0
C-
      double precision    deltaS, f, deltaT
C- data
      data    s, Emax, DeltaMax  /1.35, 318.0, 1.13/
      data    A0, A1, A2         /20.69989, -7.07605, 0.483547/
      data    E0                 /56.914686/

      EpEm   = Ep / Emax
      deltaS = DeltaMax * s * EpEm / ( s-1.d0 + EpEm**s ) 
      ! formula di Furman
      if( EpEm.ge.1.d0 )then
         eta1 = deltaS
         return
      end if
      LnEpE0 = dlog( Ep+E0 )
      f      = A0 + A1*LnEpE0 + A2*LnEpE0*LnEpE0 
      f      = dexp(f)
      deltaT = deltaS / (1.0-f)
      eta1   = deltaT
      
      return
      end

      function CTsec(E0)
c     ==============================================================
c     VERSION:         0.1
c     LAST MOD:      FEB/07
c     MOD AUTHOR:    G. Fubiani
c     COMMENTS:      CT factor for true secondary electron emission 
c                    with an incident electron E0 in eV.
c     NOTE:          Replaces data from A. Simonin ELSTOP Functions
c                    New set of data comes from Koshikawa et al. 
c                    [J. of Phys. D 6, 1360 (1973)] Fig 7.
c                    CT is calculated as follow:
c                    CT(E0)=log(eta_th/eta0)/(1.-cos(theta))
c                    An average over 3 or 4 points in taken (20, 30, 40 
c                    and 60 degrees)
c     --------------------------------------------------------------
      implicit none
      double precision    E0, CTsec
      double precision    E(7), CT(7)
      integer i


      data   E /    0.d0,  500.d0, 1000.d0, 2000.d0, 3000.d0, 6000.d0,
     >            10000.d0  /
      data CT /  0.692d0, 0.692d0, 0.731d0, 0.835d0 ,0.956d0, 0.971d0,
     >            1.005d0   /
     
   !   data   E /    0.d0, 499.d0,  500.d0, 2000.d0, 3000.d0, 6000.d0,
   !  >            10000.d0  /
   !   data CT /  0.d0, 0.d0, 0.692d0, 0.835d0 ,0.956d0, 0.971d0,
   !  >            1.005d0   /
      CTsec = CT(1)

      do 100 i=2,7,1
         if( E0.gt.E(i-1) .and. E0.le.E(i) ) then
            CTsec = CT(i-1) + ( E0-E(i-1) )*( CT(i)-CT(i-1) )/
     &           ( E(i)-E(i-1) )
         endif
 100  continue

      if( E0.ge.E(7) ) CTsec = CT(7)  
      return
      end
