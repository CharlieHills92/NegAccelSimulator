      subroutine make_sec_wall_part(Ysec,Ek,iseed,u,t,n,x_g,y_g,z_g,
     &     ptype_info)
c     ==============================================================
c     VERSION:         0.1
c     LAST MOD:       SEP/06
c     MOD AUTHOR:    G. Fubiani
c     COMMENTS:      generate secondary particle's 6D coordinates
c     NOTE:                   /
c     --------------------------------------------------------------

      implicit none
      include 'ode_arrays.h'
      include 'constants.h'
      include 'particle_info.h'
      include 'magfield_arrays.h' !added (2013/01)
      integer iseed,ptype_info
      double precision Ysec(nmax_sec),rnd,ran2,Ek,gamma,theta,phi
      double precision u(3),t(3),n(3),x_g,y_g,z_g,V11,V22,V33,
     &     Vx,Vy,Vz,Vtot,prod,mc2

C     Ysec(i)=ux_i
C     Ysec(i+1)=uy_i
C     Ysec(i+2)=uz_i
C     Ysec(i+3)=x_i
C     Ysec(i+4)=y_i
C     Ysec(i+5)=z_i
C     Ysec(i+6)=ptype_info
C     Ysec(i+7)=sym


      if( flag_err.eq.1 )RETURN

      if(ptype_info.eq.1) mc2=mec2 ! m_e c^2 (eV)
      if(ptype_info.gt.1.and.ptype_info.le.4) mc2=mpc2 ! m_H c^2 (eV)
      if(ptype_info.eq.5) mc2=2.d0*mpc2 ! m_H2 c^2 (eV)
      if(ptype_info.gt.5.and.ptype_info.le.8) mc2=2.d0*mpc2 ! m_D c^2 (eV)
      if(ptype_info.eq.9.or.ptype_info.eq.10) mc2=4.d0*mpc2 ! m_{D2,He} c^2 (eV)
      gamma = 1.d0 + Ek/mc2
      Vtot= dsqrt( gamma**2-1.d0  )/gamma*c

      rnd=ran2(iseed)
      theta = dasin(rnd)

      rnd=ran2(iseed)
      phi = 2*pi*rnd

C
C     Emission angle secondary particle relative to surface normal
C     ie: theta and phi are in the (x',y',z') = (u,t,n) system
C
      V11 = dsin(theta)*dcos(phi)
      V22 = dsin(theta)*dsin(phi)
      V33 = dcos(theta)

C
C     Transform this back to (x,y,z)
C
      Vx = V11*u(1) + V22*t(1) + V33*n(1)
      Vy = V11*u(2) + V22*t(2) + V33*n(2)
      Vz = V11*u(3) + V22*t(3) + V33*n(3)

      if(flag_info.eq.1) then
         prod = Vx*n(1) + Vy*n(2) + Vz*n(3)
         prod = dacos(prod)
         write(6,20) prod*180/pi
 20         format('Secondary emission angle:',f6.2,1x,'degree')
      endif

      Vx=Vx*Vtot
      Vy=Vy*Vtot
      Vz=Vz*Vtot

      if( ind_sec+8 .gt. nmax_sec )then
C        write(6,*)'MAKE_SEC_WALL_PART - ind_sec out of range,',
C    >             ' increase nmax_sec'
C        write(6,*)'# primary particle: ',i_info
C        write(6,*)'ind_sec, nmax_sec = ',ind_sec, nmax_sec
         flag = 1
         flag_lost = 1
         flag_err = 1
         RETURN
      end if

      Ysec(ind_sec) = Vx/c*gamma
      Ysec(ind_sec+1) = Vy/c*gamma
      Ysec(ind_sec+2) = Vz/c*gamma
      Ysec(ind_sec+3) = x_g
      Ysec(ind_sec+4) = y_g
      Ysec(ind_sec+5) = z_g
      Ysec(ind_sec+6) = ptype_info
      Ysec(ind_sec+7) = sym

      ind_sec = ind_sec + 8
      if (ptype_info.eq.11) print*,Vtot
      return
      end
