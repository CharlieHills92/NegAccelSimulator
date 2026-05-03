      subroutine ion_surf_coll(Ek,Elost,theta,iseed,u,t,n,
     &        x_g,y_g,z_g,Ysec)
c     ==============================================================
c     VERSION:         0.1
c     LAST MOD:       SEP/06
c     MOD AUTHOR:    G. Fubiani
c     COMMENTS:      Calculates ion collision with surface
c                    Theta corresponds to the impact angle between
c                    particle and wall.
c     NOTE:          RN, RE comes from ORNL Redbook Vol 3, E22 and E24
c                    eta_ion is a combinaison of Redbook data Vol 3 
c                    C10 (energy range 2 keV to 150 keV) and renormalized
c                    data from Svensson et al PRB 25 (range 150 keV to 400
c                    keV). mu_e is taken from proton impacts on Ni (Redbook
c                    Vol 3).
c                    
c     --------------------------------------------------------------
      implicit none
      include 'particle_info.h'
      include 'ode_arrays.h'
      integer iseed,ptype_info
      double precision Ysec(nmax_sec),u(3),t(3),n(3),x_g,y_g,z_g
      double precision Ek,theta,etaBACK,mu_e,mu_i,etaSEC,RE,RN,eta_ion,
     &     rnd,ran2,Eout,Elost, sput_norm, eta_sput_Cu_plus
c     Note: Ek in already in eV
c     Initialization of Elost. It represent the energy lost by particle 
c     after collision. For now it is set equal to Ek
      Elost=Ek

c
c     Set mu_i
c     Right now enforce etaBACK twice bigger at grazing incidence 
c
      mu_i=0.5d0 ! Warning: this is a free parameter

      etaBACK = RN(Ek*1.d-3)
      etaBACK = etaBACK/( (1.d0-mu_i)*dcos(theta) + mu_i  )
      if(etaBACK.gt.1.d0) etaBACK=1.d0 ! Obviously

c
c     Set mu_e
c     Enforce eta_ion = eta{max}_ion*exp[mu_e*(1-cos(theta))]
c
      mu_e=1.45d0

      etaSEC = eta_ion(Ek*1.d-3)
      if(ptype.eq.5.or.ptype.eq.9.or.ptype.eq.10) then  ! For H2, D2 or He
         etaSEC=etaSEC*2.90d0/1.32d0 ! See Redbook Vol 1, C-10
      endif
     
      etaSEC  = etaSEC*dexp( mu_e*(1.d0-dcos(theta)) ) ! Angular correction

     

      if(flag_info.eq.1)then
      !if (1.eq.1)then ! modified AP
         write(6,10)etaBACK,etaSEC,Ek*1.d-3
 10      format('Probabilities: backscatter=',f5.2,
     &          ', secondary=',f5.2,' @ Ek=',f8.2,' keV')
      endif

C
C     Backscattered ions
C     Ek in RN and RE must be in keV.
C
      rnd=ran2(iseed)
      if(rnd.le.etaBACK.and.Elost.ge.1.d0) then ! Assume Ek > 1eV   

         Eout=RE(Ek*1.d-3)/RN(Ek*1.d-3)*Ek ! ORNL Redbook Vol 3, E2

         if( (Elost-Eout).ge.0.d0 ) then ! Ensure that Ek > Eout

            Elost=Elost-Eout ! Energy lost by incoming particle (eV)

            if(flag_info.eq.1)then         
               write(6,30) pname(ptype),Eout*1.d-3
 30            format(a2,' is backscattered. Eout=',f8.2,' keV')
            endif

            n_sec=n_sec+1
            flag_sec=1
            ptype_info=ptype
            call make_sec_wall_part(Ysec,Eout,iseed,u,t,n,x_g,y_g,z_g,
     &           ptype_info)

         else

            call checknetcharge(charge(ptype)) ! Cnt net charge absorbed by grid   

         endif

      else

      call checknetcharge(charge(ptype)) ! Cnt net charge absorbed by grid  

      endif

C
C     Secondary electrons. 
C     
      call generate_sec_e(etaSEC,Elost,Ysec,iseed,u,t,n,x_g,y_g,z_g,26)
    !  call generate_Cu_plus(eta_sput_Cu_plus,Elost,Ysec,iseed,u,t,n,x_g,
    ! & y_g,z_g,26)
      return
      end

      function eta_ion(E0)
c     ==============================================================
c     VERSION:         0.1
c     LAST MOD:       SEP/06
c     MOD AUTHOR:    G. Fubiani
c     COMMENTS:      Redbook data Vol 3 C10 (energy range 2 keV to 150 keV) 
c                    and renormalized data from Svensson et al PRB 25 (range 
c                    150 keV to 400 keV).
c     NOTE:          (1) Fit log(eta_ion) = sum^5_{i=1} a_i log(Ek)**i was
c                    obtained using gnuplot.
c                    Ek must be in keV
c                    (2) For Deuterium and H2 we assume that a velocity 
c                    identical to the one of an Hydrogen atom gives the 
c                    same cross section, i.e., in terms of kinetic energies, 
c                    for a given Deuterium energy, the cross section is 
c                    identical to Hydrogen with half that energy.
c                    For Di-deuterium, we take a fouth of Hydrogen energy,
c                    same for Helium 
c     --------------------------------------------------------------
      implicit none
      include 'particle_info.h'
      integer i,ntot 
      parameter (ntot=5)
      double precision a(ntot),eta_ion,E0

      if(ptype.ge.5) E0=E0/2.d0 ! For D or H2 Ek_H == Ek_D(H2)/2
      if(ptype.eq.9) E0=E0/4.d0 ! For Di-deuterium Ek_H == Ek_D2/4
      if(ptype.eq.10) E0=E0/4.d0 ! For Helium Ek_H == Ek_He/4

      a(1)=-1.3854d0
      a(2)=0.544935d0
      a(3)=-0.0124299d0
      a(4)=-0.00307087d0
      a(5)=-0.000558421

      eta_ion=0.d0

      do i=1,ntot
         eta_ion = eta_ion + a(i)*dlog(E0)**(i-1)
      enddo

      eta_ion=dexp(eta_ion)

      return
      end

      function RN(E0)
c     ==============================================================
c     VERSION:         0.1
c     LAST MOD:       SEP/06
c     MOD AUTHOR:    G. Fubiani
c     COMMENTS:      Data from ORNL Redbook Vol 3, E22
c     NOTE:          (1) Fit log(RN) = sum^5_{i=1} a_i log(Ek)**i was
c                    obtained using gnuplot.
c                    Ek must be in keV
c                    (2) For Deuterium and H2 we assume that a velocity 
c                    identical to the one of an Hydrogen atom gives the 
c                    same cross section, i.e., in terms of kinetic energies, 
c                    for a given Deuterium energy, the cross section is 
c                    identical to Hydrogen with half that energy. 
c                    For Di-deuterium, we take a fouth of Hydrogen energy,
c                    same for Helium 
c     --------------------------------------------------------------
      implicit none
      include 'particle_info.h'
      integer i,ntot 
      parameter (ntot=5)
      double precision a(ntot),RN,E0

      if(ptype.ge.5) E0=E0/2.d0 ! For D or H2 Ek_H == Ek_D(H2)/2
      if(ptype.eq.9) E0=E0/4.d0 ! For Di-deuterium Ek_H == Ek_D2/4
      if(ptype.eq.10) E0=E0/4.d0 ! For Helium Ek_H == Ek_He/4

      a(1)=-1.55626d0
      a(2)=-0.2385d0
      a(3)=-0.0673113d0
      a(4)=-0.00899714d0
      a(5)=3.23125d-6

      RN=0.d0

      do i=1,ntot
         RN= RN + a(i)*dlog(E0)**(i-1)
      enddo

      RN=dexp(RN)

      return
      end

      function RE(E0)
c     ==============================================================
c     VERSION:         0.1
c     LAST MOD:       SEP/06
c     MOD AUTHOR:    G. Fubiani
c     COMMENTS:      Data from ORNL Redbook Vol 3, E23
c     NOTE:          (1) Fit log(RE) = sum^5_{i=1} a_i log(Ek)**i was
c                    obtained using gnuplot.
c                    Ek must be in keV
c                    (2) For Deuterium and H2 we assume that a velocity 
c                    identical to the one of an Hydrogen atom gives the 
c                    same cross section, i.e., in terms of kinetic energies, 
c                    for a given Deuterium energy, the cross section is 
c                    identical to Hydrogen with half that energy. 
c                    For Di-deuterium, we take a fouth of Hydrogen energy,
c                    same for Helium
c     --------------------------------------------------------------
      implicit none
      include 'particle_info.h'
      integer i,ntot 
      parameter (ntot=5)
      double precision a(ntot),RE,E0

      if(ptype.ge.5) E0=E0/2.d0 ! For D or H2 Ek_H == Ek_D(H2)/2
      if(ptype.eq.9) E0=E0/4.d0 ! For Di-deuterium Ek_H == Ek_D2/4
      if(ptype.eq.10) E0=E0/4.d0 ! For Helium Ek_H == Ek_He/4


      a(1)=-2.33586d0
      a(2)=-0.450527d0
      a(3)=-0.0732621d0
      a(4)=-0.00687502d0
      a(5)=-4.56457d-6

      RE=0.d0

      do i=1,ntot
         RE= RE + a(i)*dlog(E0)**(i-1)
      enddo

      RE=dexp(RE)

      return
      end

      

