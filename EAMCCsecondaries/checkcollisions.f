      subroutine checkcollisions(Y,Ysec,Ek,iseed,i)
c     ==============================================================
c     VERSION:         0.1
c     LAST MOD:       SEP/06
c     MOD AUTHOR:    G. Fubiani
c     COMMENTS:      Checks if particle encounter a collision
c     NOTE:               /
c     --------------------------------------------------------------
      implicit none
      include 'slac_arrays.h'
      include 'ode_arrays.h'
      include 'particle_info.h'
      include 'constants.h'
      integer iseed !,side nf side variabile globale 
      double precision Y(nmax),Ysec(nmax_sec),Ek
      integer i

      if( flag_err.eq.1 )RETURN
C
C     Collisions with grids
C     avoid mistaken interception by freshly emitted particles (electrons)
C      
      if(NSTEP.eq.0)flag_g=0
      if(flag_g.eq.1) call gridcollisions(Y,Ysec,Ek,iseed,i)
      if( flag_err.eq.1 )RETURN

C
C     Collisions with background gas (H2/D2)
C 
     
      if(ptype.gt.1.and.flag_src.eq.0) call gascollisions(Y,Ysec,Ek,
     &     iseed)               ! flag_src=1 == inside ion source!

      return
      end

      subroutine gridcollisions(Y,Ysec,Ek,iseed,i)
c     ==============================================================
c     VERSION:         0.1
c     LAST MOD:       SEP/06
c     MOD AUTHOR:    G. Fubiani
c     COMMENTS:      Checks if particle encounter a collision with wall
c     NOTE:          n(3) is the normal vector to the surface, pointing
c                    outward, t(3) the tangent vector and u(3)=t^n the
c                    3rd vector defining an orthonormal basis.
c     --------------------------------------------------------------
      implicit none
      include 'slac_arrays.h'
      include 'ode_arrays.h'
      include 'magfield_arrays.h'
      include 'particle_info.h'
      include 'constants.h'
      integer iseed,ierr !,side
      double precision Y(nmax),Ysec(nmax_sec),u(3),t(3),n(3),
     &     partvect(3),norm,prod,theta,Ek,Elost,x_g,y_g,z_g
      double precision vx,vy,vz,gamma,delta1
      double precision xCentre !added by Kisaki
      integer iz, icount,i
      integer ix,iy,ii,jj,bound_l,bound_r,bound_b,bound_u !added (2013/01)
      data icount /0/

      ierr=0
C      xCentre=(dxmax3D+dxmin3D)/2.d0 !added by Kisaki

C
C     Find vector normal (pointing outward) to the wall. x_g, y_g
C     z_g are the location of impact in grid unit
C
      
C         call surfacevector(Y,u,t,n,x_g,y_g,z_g,side)  ! 2D routine   nf esclusa è di slac lavora in coordinate cilindriche
C      else
C         if( neFiles.le.1 )then    !neFiles=1 se c'è solo mappa SLAC, ma ora non esiste più, tolta la definizione            
C            write(6,*)'gridcoll - part outside domain, XYZ=',
C     >      1000*Y(4),1000*Y(5),1000*Y(6),' mm'
C            flag_err = 1
C            return
C         else
C            call surf3Dvector(Y,u,t,n,x_g,y_g,z_g,side)   ! 3D routine
C         end if
C      end if
C***********************************************************************************
      if( Y(6).le.dZmin3D+deltaz_E(1) )then   
            write(6,*)'gridcoll - part outside domain, XYZ=',
     &                1000*Y(4),1000*Y(5),1000*Y(6),' mm'
             flag_err = 1
             return
      end if                 !NF se la z è inferiore a 0, limite inferiore del dominio, allora errore e return
             
      call surf3Dvector(Y,u,t,n,x_g,y_g,z_g,i)   ! 3D routine
      
      if( flag_err.eq.1 )return

      if(flag_info.eq.1) then 
         write(*,10) pname(ptype),x_g*1.d3,y_g*1.d3,z_g*1.d3
 10      format(a2,' has collided with grid @ (x,y,z) =',
     &        '(',f6.2,',',f6.2,',',f6.2,') [mm]')
         write(*,11) grid_ind,sym
 11      format('Grid segment # ',i2,', hole #',i2)
      endif


C
C     Calculate particle incoming angle
C
      gamma = dsqrt( 1.d0 + Y(1)**2 + Y(2)**2 + Y(3)**2 )

      vx = Y(1)/gamma*c
      vy = Y(2)/gamma*c
      vz = Y(3)/gamma*c

      norm=dsqrt( vx**2 + vy**2 + vz**2)
      if( norm.gt.0 )then !||partvect|| = 1
         partvect(1)=vx/norm
         partvect(2)=vy/norm
         partvect(3)=vz/norm
      else
         partvect(1)=0
         partvect(2)=0
         partvect(3)=0
      end if

c     a.b = ||a||*|b||*cos(a,b)
c     here ||a||=|b||=1
      prod= partvect(1)*n(1) + partvect(2)*n(2) + partvect(3)*n(3)
      prod=-prod ! correct for particle vector pointing inward
      theta= dacos(prod)

      if( ABS(theta).gt.pi/2.d0 .or. norm.eq.0.d0 ) then ! Found bad angle 
C
C Rather then deleting the particle, I declare the impact angle "good"
C and allow it to possibly produce secondary particles !
         cnt_badth=cnt_badth+1

c        theta = pi/4
c        goto 15
C
C - Disabled -
c        ierr=1
c        cnt_badth=cnt_badth+1
c        call checknetcharge(charge(ptype)) ! Collect charge on grid

c ---------------NF quanto segue è il calcolo del nodo iz ed è stato sostituito dalle due righe successive
c         if( Y(6)/deltaZ_E .ge. 0 )then            è uguale alla espressione in get3Dpot 
c            iz = int( Y(6)/deltaZ_E + 0.0 ) + 1    NF iz è indice lungo z precedente la particella
c         else
c            iz = int( Y(6)/deltaZ_E - 1.0 ) + 1
c         end if
c---------------------------------------------------
		 
		numel_zinf=floor((Y(6)-dzmin3D)/delta)+1   !NF Y(6) è la zi, coordinata della particella i lungo z
		iz=index_z(numel_zinf)           !NF indice z del nodo precedente la particella
C     
C     -----    added (2013/01)    -----   NF stesso discorso per x e y, la modifica è uguale a quella in get3Dpot
C     
C         if( Y(4)/deltaX_E .ge. 0 )then
C            ix = int( Y(4)/deltaX_E + 0.0 )
C         else
C            ix = int( Y(4)/deltaX_E - 1.0 )
C         end if
C         if( Y(5)/deltaY_E .ge. 0 )then
C            iy = int( Y(5)/deltaY_E + 0.0 )
C         else
C            iy = int( Y(5)/deltaY_E - 1.0 )
C         end if

        numel_xinf=floor((Y(4)-dxmin3D)/delta)+1    !NF definizione nodi ix e iy come in get3Dpot
        numel_yinf=floor((Y(5)-dymin3D)/delta)+1
      
        ix=index_x(numel_xinf)  
        iy=index_y(numel_yinf)
		 
         do jj=1,naper_y
         do ii=1,naper_x
            bound_b = ( jjmax(ii,jj-1,iz) + jjmin(ii,jj,iz) )/2 + 1
            bound_u = ( jjmax(ii,jj,iz) + jjmin(ii,jj+1,iz) )/2
            if(jj.eq.1) bound_b = iyEmin
            if(naper_y.eq.1 .or. jj.eq.naper_y) bound_u = iyEmax

            bound_l = ( iimax(ii-1,jj,iz) + iimin(ii,jj,iz) )/2 + 1
            bound_r = ( iimax(ii,jj,iz) + iimin(ii+1,jj,iz) )/2
            if(ii.eq.1) bound_l = ixEmin
            if(naper_x.eq.1 .or. ii.eq.naper_x) bound_r = ixEmax

            if(ix.ge.bound_l .and. ix.le.bound_r .and.
     >           iy.ge.bound_b .and. iy.le.bound_u) then
               Rpart_old = dsqrt( (old_part_pos(1)-Xgrid(ii,jj,iz))**2 +
     >              (old_part_pos(2)-Ygrid(ii,jj,iz))**2 )
            end if
         end do
         end do
C
C     --------------------------------
C
C         Rpart_old = dsqrt( (old_part_pos(1)-Xgrid(iz))**2 +
C     >                      (old_part_pos(2)-Ygrid(iz))**2 ) !comment out by Kisaki
C         Rpart_old = dsqrt( ((old_part_pos(1)-xCentre)-Xgrid(iz))**2 +
C     >                      (old_part_pos(2)-Ygrid(iz))**2 ) !added by Kisaki
         icount = icount+1
         if(icount.le.-1)then
         write(6,*)'GRIDCOLLISIONS -- Bad impact angle surf normal: ',
     >              180*theta/pi,' degrees.'
         write(6,*)'grid  X Y Z=',1000*x_g,1000*y_g,1000*z_g,' mm'
         write(6,*)'oldP  X Y Z=',1000*old_part_pos(1),
     >    1000*old_part_pos(2),1000*old_part_pos(3),' mm',
     >             ' Rpart_OLD=',1000*Rpart_old,' mm.'
         write(6,*)'part  X Y Z=',1000*Y(4),1000*Y(5),1000*Y(6),' mm',
     >             ' Rpart=',1000*Rpart,' mm.'
C        write(6,*)'   IX IY IZ=',ix,iy,iz    
         write(6,*)'surf normal=',n(1),n(2),n(3)
         write(6,*)'part vector=',partvect(1),partvect(2),partvect(3)
         delta1 = dsqrt( (old_part_pos(1)-Y(4))**2 + 
     &                   (old_part_pos(2)-Y(5))**2 +
     &                   (old_part_pos(3)-Y(6))**2 )
         write(6,*)'Stepsize was',1000*Delta1,' mm.  NSTEP=',nstep,
     >             ' part=',pname(ptype)
         end if
         if( icount.ge.100 .and. icount.le.-1 )then
            if(icount.eq.100)write(6,*)'LAST SUCH MESSAGE !'
            write(6,*)'NO MORE SUCH MESSAGES !'
         end if
         theta = pi/4

      endif

 15   continue
      if(flag_info.eq.1) then 
         write(*,20) theta/pi*180.d0 
 20      format('Impact angle with surface normal:',f6.2,1x,'degree')
      endif

C
C     Check wall collisions for electrons
C
      if(ptype.eq.1.and.ierr.eq.0) call e_surf_coll(Ek,Elost,
     &     theta,iseed,u,t,n,x_g,y_g,z_g,Ysec)

C
C     Check wall collisions for ions
C
      if(ptype.gt.1.and.ierr.eq.0) call ion_surf_coll(Ek,Elost,
     &     theta,iseed,u,t,n,x_g,y_g,z_g,Ysec)

C
C     Save impact information
C     
      if(ierr.eq.0) then ! Collect energy lost on grid
         call savimpact(Y,x_g,y_g,z_g,Elost) 
      else
         call savimpact(Y,x_g,y_g,z_g,Ek)
      endif

      return
      end


      subroutine gascollisions(Y,Ysec,Ek,iseed)
c     ==============================================================
c     VERSION:         0.1
c     LAST MOD:       SEP/06
c     MOD AUTHOR:    G. Fubiani
c     COMMENTS:      Checks if particle encounter a collision with gas
c     WARNING:       Having Ek*1.d-3 as an input to the stripping and 
c                    ionization functions does not induce any change of 
c                    the true value of Ek outside the function even when 
c                    Ek is set to Ek/2 for heavy particles. This has been 
c                    checked!
c     --------------------------------------------------------------
      implicit none
      include 'slac_arrays.h'
      include 'ode_arrays.h'
      include 'particle_info.h'
      include 'constants.h'
      include 'density_info.h'
      include 'magfield_arrays.h'
      integer nmax_nu
      parameter (nmax_nu=7)
      integer iseed,coll_type,i,cnt_coll,nsort,indx(nmax_nu)
      double precision Y(nmax),Ysec(nmax_sec),Ek,
     &     Deltal,ng,sigma_stripping_H2,P_tot,nu(nmax_nu),nu_sum,
     &     nu_tot,rnd,ran2,sigma_dblstripping_H2,sigma_ionization_H0_H2,
     &     sort_arr(nmax_nu),sigma_stripping_He,sigma_dblstripping_He,
     &     sigma_ionization_H0_He,sigma_ioni_H0_loss_He,
     &     sigma_stripping_H0_He

C
C     old_part_pos(1) == x_old
C     old_part_pos(2) == y_old
C     old_part_pos(3) == z_old
C     Y(4) == x
C     Y(5) == y
C     Y(6) == z
C     Ek is already in eV
C

      cnt_coll=0 ! initialize counter
      coll_type=0

C
C     Find Deltal (distance advanced by particle during time interval DT)
C
      Deltal = dsqrt( (old_part_pos(1)-Y(4))**2 + 
     &                (old_part_pos(2)-Y(5))**2 +
     &                (old_part_pos(3)-Y(6))**2 )


C
C     Grab local background gas density profile (in m-3)
C
      
      do i=2,cnt_mag,1 ! Perform linear interpolation

         if( Y(6).gt.density(i-1,1).and.Y(6).le.density(i,1) ) then

            ng= density(i-1,2) + ( Y(6) - density(i-1,1) )*
     &           ( density(i,2) - density(i-1,2) )/
     &           ( density(i,1) - density(i-1,1) )

         endif

      enddo

      if ( Y(6) .gt. density(cnt_mag,1) ) ng=0.0
      
      if (flag_ext_x .eq. 1 .and. abs(Y(4)) .gt. x_uni(ixEmax-1)) then    !nf disattiva collisioni mettendo ng=0 serve solo per tomografia
      ng=0.0
      end if
      
      if (flag_ext_y .eq. 1 .and. abs(Y(5)) .gt. y_uni(iyEmax-1)) then
      ng=0.0
      end if
      
c      if (flag_ext_z .eq. 1 .and. Y(6) .gt. z_uni(izEmax2-1)) then
c      ng=0.0
c      end if
C
C     =============================================================
C     STRIPPING LOSSES
C
C     Cumulative stripping loss prob. (Ek must be in keV)
C     for H- + H2 -> H0 + H2 + e- (or D- + D2 -> D0 + D2 + e-)
C     Redbook data Vol 1, F8
C     for H- + He -> H0 + He + e- (or D- + He -> D0 + He + e-)
C     Redbook data Vol 1, F13
C     for H0 + He -> H+ + He + e- (or D0 + He -> D+ + He + e-)
C     Redbook data Vol 1, E11
C
C     Double stripping loss prob. (Ek must be in keV)
C     for H- + H2 -> H+ + H2 + 2e- (D- + D2 -> D+ + D2 + 2e-)
C     Redbook data Vol 1, F10
C     for H- + He -> H+ + He + 2e- (or D- + He -> D+ + He + 2e-)
C     Redbook data Vol 1, F15
C
C
C     =============================================================
C     IONIZATION OF BACKGROUND GAS
C
C     H2 (D2) background gas ionization probability (Ek must be in keV)
C     for H- + H2 -> H- + H2+ + e- (D- + D2 -> D- + D2+ + e-) and
C     for H0 + H2 -> H0 + H2+ + e- (D0 + D2 -> D0 + D2+ + e-)
C     Redbook data Vol 1, D4

C     He background gas ionization probability (Ek must be in keV)
C     for H- + He -> H- + He+ + e- (D- + He -> D- + He+ + e-) and
C     for H0 + He -> H0 + He+ + e- (D0 + He -> D0 + He+ + e-)
C     DuBois et al., PRA 40, 3605 (1989) Fig 3(a)
C
C     Double ionization (background gas  + projectile)
C     for H- + He -> H0 + He+ + 2e- (D- + He -> D0 + He+ + 2e-)
C     for H0 + He -> H+ + He+ + 2e- (D0 + He -> D+ + He+ + 2e-)
C     DuBois et al., PRA 40, 3605 (1989) Fig 3(c)
C
C     ptype=2 == H- (ptype=6 == D-), ptype=3 == H0 (ptype=7 == D0)
C     ptype=10 == He

C
C     Check for H- (D-) collision type
C

      if( ptype.eq.2.or.ptype.eq.6 ) then 


C     Background gas same as primary particle
            
         nu(1) = ng*sigma_stripping_H2(Ek*1.d-3) ! Single stripping
         nu(2) = ng*sigma_dblstripping_H2(Ek*1.d-3) ! Double stripping
         nu(3) = ng*sigma_ionization_H0_H2(Ek*1.d-3) ! Ionization
         
C     Background gas is Helium
         nu(4) = ng_He*sigma_stripping_He(Ek*1.d-3) ! Single stripping
         nu(5) = ng_He*sigma_dblstripping_He(Ek*1.d-3) ! Double stripping
         nu(6) = ng_He*sigma_ionization_H0_He(Ek*1.d-3) ! Pure ionization
         nu(7) = ng_He*sigma_ioni_H0_loss_He(Ek*1.d-3) ! Ionization + H- (D-) ele. loss
         

         if(ionization.eq.'n'.or.flag_g.eq.1) then ! Excluding ionization
            nu(3) = 0.d0
            nu(6) = 0.d0
            nu(7) = 0.d0            
         endif

         if(gtype.eq.10) then ! If Helium background gas is added
            nsort=7
         else
            nsort=3               
         endif
            
         nu_tot = 0.d0 ! Initialize nu_tot
         do i=1,nsort
            nu_tot = nu_tot + nu(i)
         enddo


         P_tot = 1.d0 - dexp( -nu_tot*Deltal ) ! Calculate coll. prob.

         rnd=ran2(iseed)
         if( rnd.le.P_tot ) then ! A reaction occured

            do i=1,nmax_nu ! Initialize sort_arr 
               sort_arr(i) = 0.d0 
               indx(i) = 0
            enddo

            do i=1,nsort ! Sort reactions
               sort_arr(i) = nu(i)/nu_tot
            enddo

            call indexx(nsort,sort_arr,indx)

            rnd=ran2(iseed)
            nu_sum=0.d0
            do i=1,nsort ! Select one reaction [nu_1/nu_tot < (nu_1+nu_2)/n_tot < ...]

               if( rnd.le.( nu_sum + sort_arr(indx(i)) ) ) then
                  coll_type=indx(i) ! coll_type = 1-7
                  goto 55
               endif

               nu_sum = nu_sum + sort_arr(indx(i))

            enddo
         
 55      endif

      endif !  End Check for H- (D-) collision type

C
C     Check for H0 (D0) collision type
C

      if( ptype.eq.3.or.ptype.eq.7 ) then

         if(ionization.eq.'y'.and.flag_g.eq.0) then ! Ionization by H0 (D0)

C     Background gas same as primary particle

            nu(1) = ng*sigma_ionization_H0_H2(Ek*1.d-3) ! Pure ionization

C     Background gas is helium

            nu(2) = ng_He*sigma_ionization_H0_He(Ek*1.d-3) ! Pure ionization
            nu(3) = ng_He*sigma_ioni_H0_loss_He(Ek*1.d-3) ! Ionization + H0 (D0) electron loss
            nu(4) = ng_He*sigma_stripping_H0_He(Ek*1.d-3) ! H0 (D0) single stripping
               
            if(gtype.eq.10) then  ! If He background gas is added
               nsort=4
            else
               nsort=1
            endif
            
            nu_tot = 0.d0 ! Initialize nu_tot
            do i=1,nsort
               nu_tot = nu_tot + nu(i)
            enddo

            P_tot = 1.d0 - dexp( -nu_tot*Deltal ) ! Calculate coll. prob.

            rnd=ran2(iseed)
            if( rnd.le.P_tot ) then ! A reaction occured

               do i=1,nmax_nu    ! Initialize sort_arr 
                  sort_arr(i) = 0.d0 
                  indx(i) = 0
               enddo

               do i=1,nsort ! Sort reactions
                  sort_arr(i) = nu(i)/nu_tot
               enddo

               call indexx(nsort,sort_arr,indx)

               rnd=ran2(iseed)
               nu_sum=0.d0
               do i=1,nsort ! Select one reaction [nu_1/nu_tot < (nu_1+nu_2)/n_tot < ...]

                  if( rnd.le.( nu_sum + sort_arr(indx(i)) ) ) then
                     coll_type = 7 + indx(i) ! coll_type = 8 to 11
                     goto 56
                  endif

                  nu_sum = nu_sum + sort_arr(indx(i))

               enddo
               
 56         endif

         endif

      endif ! End check for H0 (D0) collision type

C
C     Generate secondary particles
C     

! if single stripping of H- (D-) or H0 (D0)
      if(coll_type.eq.1.or.coll_type.eq.4.or.coll_type.eq.11) then

         if(flag_info.eq.1) then         
            write(6,10) Y(4)*1.d3,Y(5)*1.d3,Y(6)*1.d3
 10         format ('Stripping reaction occurred @ (x, y, z)=',
     &           '(',f6.2,',',f6.2,',',f6.2,') [mm]')
         endif

         n_sec=n_sec+2 ! Generate two sec. particles (H/D or H+/D+,e)
         flag_sec=1 

         call make_sec_part(Y,Ysec,coll_type,iseed)

         cnt_sec(ptype+1) = cnt_sec(ptype+1) + 1 ! Counter for H0(D0)/H+(D+)
         cnt_sec(1) = cnt_sec(1) + 1  ! Counter for sec. elec.
         flag=1 ! kill current (H- or H0) particle
         flag_des=1

      endif

! if double stripping of H- (D-)
      if(coll_type.eq.2.or.coll_type.eq.5) then
         
         if(flag_info.eq.1) then         
            write(6,20) Y(4)*1.d3,Y(5)*1.d3,Y(6)*1.d3
 20         format ('Double stripping reaction occurred @ (x, y, z)=',
     &           '(',f6.2,',',f6.2,',',f6.2,') [mm]')
         endif

         n_sec=n_sec+3 ! Generate 3 sec. particles (H/D,2e)
         flag_sec=1 

         call make_sec_part(Y,Ysec,coll_type,iseed)

         cnt_sec(ptype+2) = cnt_sec(ptype+2) + 1 ! Counter for H+/D+
         cnt_sec(1) = cnt_sec(1) + 2  ! Counter for sec. elec.
         flag=1 ! kill current (H-) particle
         flag_des=1

      endif


! if ionization of background gas by H- (D-), H0 (D0)
      if(coll_type.eq.3.or.coll_type.eq.6.or.coll_type.eq.8.or.
     &     coll_type.eq.9) then

         if(coll_type.eq.3) i=ptype+3
         if(coll_type.eq.8) i=ptype+2
         if(coll_type.eq.6.or.coll_type.eq.9) i=gtype 
            
         if(flag_info.eq.1) then         
            write(6,30) pname(i),pname(ptype),Y(4)*1.d3,
     &           Y(5)*1.d3,Y(6)*1.d3
 30         format ('Ionization of ',1x,a2,1x,'by',1x,a2,1x,
     &           'occurred @ (x, y, z)=','(',f6.2,',',f6.2,',',
     &           f6.2,') [mm]')
         endif

         n_sec=n_sec+2          ! Generate 2 sec. particles (H2+/D2+,He+,e)
         flag_sec=1             

         call make_sec_part(Y,Ysec,coll_type,iseed)

         cnt_sec(i) = cnt_sec(i) + 1 ! Counter for H2+/D2+/He+
         cnt_sec(1) = cnt_sec(1) + 1 ! Counter for sec. elec.
!     Note: the incident particle is NOT killed (by using flag=1)

      endif

! if ionization of He background gas and primary particle H*(D*)
      if(coll_type.eq.7.or.coll_type.eq.10) then
            
         if(flag_info.eq.1) then         
            write(6,31) pname(gtype),pname(ptype),Y(4)*1.d3,
     &           Y(5)*1.d3,Y(6)*1.d3
 31         format ('Ionization of both',1x,a2,1x,'and',1x,a2,1x,
     &           'occurred @ (x, y, z)=','(',f6.2,',',f6.2,',',
     &           f6.2,') [mm]')
         endif

         n_sec=n_sec+4          ! Generate 2 sec. particles (H0/D0,He+,2e) or
         flag_sec=1             ! (H+/D+,He+,2e)

         call make_sec_part(Y,Ysec,coll_type,iseed)

         cnt_sec(ptype+1) = cnt_sec(ptype+1) + 1 ! Counter for H+/D+ or H0/D0
         cnt_sec(gtype) = cnt_sec(gtype) + 1 ! Counter for He+
         cnt_sec(1) = cnt_sec(1) + 2 ! Counter for sec. elec.
         flag=1 ! kill current (H-/D- or H0/D0) particle
         flag_des=1

      endif

      return
      end

      function sigma_stripping_H2(Ek)
c     ==============================================================
c     VERSION:         0.1
c     LAST MOD:       SEP/06
c     MOD AUTHOR:    G. Fubiani
c     COMMENTS:      Calculate cross section for electron detachment in H-(D-) 
c                    i.e., reaction H- + H2 -> H + H2 + e- (D- + D2 -> D + D2 
c                    + e-)
c                    Data from ORNL redbook vol. 1, F-8.
c     NOTE:          For Deuterium we assume that a velocity identical to
c                    the one of an Hydrogen atom gives the same cross 
c                    section, i.e., in terms of kinetic energies, for a 
c                    given Deuterium energy, the cross section is identical
c                    to Hydrogen with half that energy. 
c     --------------------------------------------------------------
      implicit none
      include 'particle_info.h'
      integer i,ntot
      parameter(ntot=36)
      double precision Ehyd(ntot)  ! H- energy in keV
      double precision sigH2(ntot) ! Cross-section in 10**(-20) m2
      double precision Ek,sigma_stripping_H2

      if(ptype.ge.5) Ek=Ek/2.d0 ! For Deuterium Ek_H == Ek_D/2

      data Ehyd   /             0.0023, 0.004,  0.007,
     >          0.010,  0.015,  0.020,  0.040,  0.070,
     >          0.1,    0.15,   0.2,    0.4,    0.7,
     >          1.0,    1.4,    2.0,    4.0,    7.0,
     >         10.0,   15.0,   20.0,   40.0,   70.0,
     >        100.0,  150.0,  200.0,  400.0,  700.0,
     >       1000.0, 1500.0, 2000.0, 4000.0, 7000.0,
     >      10000.0,15000.0,17000.0 /

      data sigH2  /             0.891,  1.92,  3.18,
     >          3.68,   4.14,   4.28,   4.15,  4.36,
     >          4.80,   5.60,   6.29,   8.41, 10.00,
     >         10.60,  11.00,  11.30,  11.40, 10.90,
     >         10.20,   9.24,   8.36,   6.33,  4.82,
     >          3.95,   3.06,   2.51,   1.43,  0.902,
     >          0.643,  0.441,  0.329,  0.155, 0.0813,
     >          0.0516, 0.0317, 0.0272 /


      sigma_stripping_H2 = 0.d0

C
C     Linear interpolation
C
      do i=2,ntot,1
         if( Ek.gt.Ehyd(i-1) .and. Ek.le.Ehyd(i) ) then
            sigma_stripping_H2 = sigH2(i-1) + ( Ek - Ehyd(i-1) )*
     &           ( sigH2(i) - sigH2(i-1) )/( Ehyd(i) - Ehyd(i-1) )
         end if
      enddo

      sigma_stripping_H2 = sigma_stripping_H2*1.d-20

      return
      end

      function sigma_dblstripping_H2(Ek)
c     ==============================================================
c     VERSION:         0.1
c     LAST MOD:       OCT/06
c     MOD AUTHOR:    G. Fubiani
c     COMMENTS:      Calculate cross section for double stripping. 
c                    H- + H2 -> H+ + H2 + 2e- ( D- + D2 -> D+ + D2 + 2e-)
c                    Redbook data Vol 1 F10 (energy range 1 keV to 1 MeV)
c     NOTE:          (1) Fit log(sigma) = sum^4_{i=1} a_i log(Ek)**i was
c                    obtained using gnuplot.
c                    Ek must be in keV, sigma originally in cm2
c                    (2) For Deuterium we assume that a velocity identical to
c                    the one of an Hydrogen atom gives the same cross 
c                    section, i.e., in terms of kinetic energies, for a 
c                    given Deuterium energy, the cross section is identical
c                    to Hydrogen with half that energy. 
c     --------------------------------------------------------------
      implicit none
      include 'particle_info.h'
      integer i,ntot 
      parameter (ntot=4)
      double precision a(ntot),sigma_dblstripping_H2,Ek

      if(ptype.ge.5) Ek=Ek/2.d0 ! For Deuterium Ek_H == Ek_D/2

      a(1)=-38.8333d0
      a(2)=0.859947d0
      a(3)=-0.160914d0
      a(4)=0.00113803d0

      sigma_dblstripping_H2=0.d0

      do i=1,ntot
         sigma_dblstripping_H2 = sigma_dblstripping_H2 + 
     &        a(i)*dlog(Ek)**(i-1)
      enddo

      sigma_dblstripping_H2=dexp(sigma_dblstripping_H2)*1.d-4 ! convert cm2 to m2

      return
      end

      function sigma_ionization_H0_H2(Ek)
c     ==============================================================
c     VERSION:         0.1
c     LAST MOD:       OCT/06
c     MOD AUTHOR:    G. Fubiani
c     COMMENTS:      Calculate cross section for ionization of H2 (D2). 
c                    H0 + H2 -> H0 + H2+ + e- ( D0 + D2 -> D0 + D2+ + e-)
c                    Redbook data Vol 1 D4 (energy range 2 keV to 400 keV)
c     NOTE:          (1) Fit log(sigma) = sum^9_{i=1} a_i log(Ek)**i was
c                    obtained using gnuplot.
c                    Ek must be in keV, sigma originally in cm2
c                    (2) For Deuterium we assume that a velocity identical to
c                    the one of an Hydrogen atom gives the same cross 
c                    section, i.e., in terms of kinetic energies, for a 
c                    given Deuterium energy, the cross section is identical
c                    to Hydrogen with half that energy. 
c                    (3) For ionization with H- (D-) we use the same cross
c                    section as for H0 (D0)
c     --------------------------------------------------------------
      implicit none
      include 'particle_info.h'
      integer i,ntot 
      parameter (ntot=9)
      double precision a(ntot),sigma_ionization_H0_H2,Ek

      if(ptype.ge.5) Ek=Ek/2.d0 ! For Deuterium Ek_H == Ek_D/2

      a(1)=-39.7161d0
      a(2)=1.47946d0
      a(3)=-0.0522287d0
      a(4)=-0.0236038d0
      a(5)=-0.0096701d0
      a(6)=0.00691263d0
      a(7)=-0.00234029d0
      a(8)=0.000365356d0
      a(9)=-2.02354d-5

      sigma_ionization_H0_H2=0.d0

      do i=1,ntot
         sigma_ionization_H0_H2 = sigma_ionization_H0_H2 + 
     &        a(i)*dlog(Ek)**(i-1)
      enddo

      sigma_ionization_H0_H2=dexp(sigma_ionization_H0_H2)*1.d-4 ! convert cm2 to m2
      if(Ek.le.1.d-2) sigma_ionization_H0_H2=0.d0 ! Ek < 0.01 keV 

      return
      end

      function sigma_ionization_Hpl_H2(Ek)
c     ==============================================================
c     VERSION:         0.1
c     LAST MOD:       NOV/06
c     MOD AUTHOR:    G. Fubiani
c     COMMENTS:      Calculate cross section for ionization of H2 (D2). 
c                    H+ + H2 -> H+ + H2+ + e- ( D+ + D2 -> D+ + D2+ + e-)
c                    Redbook data Vol 1 D4 (energy range 2 keV to 400 keV)
c     NOTE:          (1) Fit log(sigma) = sum^9_{i=1} a_i log(Ek)**i was
c                    obtained using gnuplot.
c                    Ek must be in keV, sigma originally in cm2
c                    (2) For Deuterium we assume that a velocity identical to
c                    the one of an Hydrogen atom gives the same cross 
c                    section, i.e., in terms of kinetic energies, for a 
c                    given Deuterium energy, the cross section is identical
c                    to Hydrogen with half that energy. 
c     --------------------------------------------------------------
      implicit none
      include 'particle_info.h'
      integer i,ntot 
      parameter (ntot=4)
      double precision a(ntot),sigma_ionization_Hpl_H2,Ek

      if(ptype.ge.5) Ek=Ek/2.d0 ! For Deuterium Ek_H == Ek_D/2

      a(1)=-45.5559d0
      a(2)=4.7871d0
      a(3)=-0.728689d0
      a(4)=0.0284248d0

      sigma_ionization_Hpl_H2=0.d0

      do i=1,ntot
         sigma_ionization_Hpl_H2 = sigma_ionization_Hpl_H2 + 
     &        a(i)*dlog(Ek)**(i-1)
      enddo

      sigma_ionization_Hpl_H2=dexp(sigma_ionization_Hpl_H2)*1.d-4 ! convert cm2 to m2
      if(Ek.le.1.d-2) sigma_ionization_Hpl_H2=0.d0 ! Ek < 0.01 keV 

      return
      end

      function sigma_stripping_He(Ek)
c     ==============================================================
c     VERSION:         0.1
c     LAST MOD:       OCT/07
c     MOD AUTHOR:    G. Fubiani
c     COMMENTS:      Calculate cross section for electron detachment in H-(D-) 
c                    i.e., reaction H- + He -> H + He + e- (D- + He -> D + He 
c                    + e-)
c                    Data from ORNL redbook vol. 1, F-13.
c     NOTE:          For Deuterium we assume that a velocity identical to
c                    the one of an Hydrogen atom gives the same cross 
c                    section, i.e., in terms of kinetic energies, for a 
c                    given Deuterium energy, the cross section is identical
c                    to Hydrogen with half that energy. 
c     --------------------------------------------------------------
      implicit none
      include 'particle_info.h'
      integer i,ntot
      parameter(ntot=36)
      double precision Ehyd(ntot)  ! H- energy in keV
      double precision sigHe(ntot) ! Cross-section in 10**(-20) m2
      double precision Ek,sigma_stripping_He

      if(ptype.ge.5) Ek=Ek/2.d0 ! For Deuterium Ek_H == Ek_D/2

      data Ehyd   /             0.0023, 0.004,  0.007,
     >          0.010,  0.015,  0.020,  0.040,  0.070,
     >          0.1,    0.15,   0.2,    0.4,    0.7,
     >          1.0,    1.4,    2.0,    4.0,    7.0,
     >         10.0,   15.0,   20.0,   40.0,   70.0,
     >        100.0,  150.0,  200.0,  400.0,  700.0,
     >       1000.0, 1500.0, 2000.0, 4000.0, 7000.0,
     >      10000.0,15000.0,17000.0 /

      data sigHe  /             1.412,  2.33,  2.94,
     >          3.03,   3.14,   3.22,   3.42,  3.70,
     >          3.91,   4.21,   4.48,   5.10,  5.60,
     >          6.02,   6.31,   6.53,   6.73,  6.34,
     >          6.08,   5.50,   5.01,   3.83,  2.80,
     >          2.24,   1.72,   1.42,   0.833, 0.527,
     >          0.400,  0.292,  0.227,  0.123, 0.0698,
     >          0.0463, 0.0282, 0.021  /


      sigma_stripping_He = 0.d0

C
C     Linear interpolation
C
      do i=2,ntot,1
         if( Ek.gt.Ehyd(i-1) .and. Ek.le.Ehyd(i) ) then
            sigma_stripping_He = sigHe(i-1) + ( Ek - Ehyd(i-1) )*
     &           ( sigHe(i) - sigHe(i-1) )/( Ehyd(i) - Ehyd(i-1) )
         end if
      enddo

      sigma_stripping_He = sigma_stripping_He*1.d-20

      return
      end


      function sigma_dblstripping_He(Ek)
c     ==============================================================
c     VERSION:         0.1
c     LAST MOD:       OCT/06
c     MOD AUTHOR:    G. Fubiani
c     COMMENTS:      Calculate cross section for double stripping. 
c                    H- + He -> H+ + He + 2e- ( D- + He -> D+ + He + 2e-)
c                    Redbook data Vol 1 F15 (energy range 0.4 keV to 1 MeV)
c     NOTE:          (1) Fit log(sigma) = sum^5_{i=1} a_i log(Ek)**i was
c                    obtained using gnuplot.
c                    Ek must be in keV, sigma originally in cm2
c                    (2) For Deuterium we assume that a velocity identical to
c                    the one of an Hydrogen atom gives the same cross 
c                    section, i.e., in terms of kinetic energies, for a 
c                    given Deuterium energy, the cross section is identical
c                    to Hydrogen with half that energy. 
c     --------------------------------------------------------------
      implicit none
      include 'particle_info.h'
      integer i,ntot 
      parameter (ntot=5)
      double precision a(ntot),sigma_dblstripping_He,Ek

      if(ptype.ge.5) Ek=Ek/2.d0 ! For Deuterium Ek_H == Ek_D/2

      a(1)=-38.8363d0
      a(2)=1.57652d0
      a(3)=-0.521506d0
      a(4)=0.0467808d0
      a(5)=-0.00168294d0

      sigma_dblstripping_He=0.d0

      do i=1,ntot
         sigma_dblstripping_He = sigma_dblstripping_He + 
     &        a(i)*dlog(Ek)**(i-1)
      enddo

      sigma_dblstripping_He=dexp(sigma_dblstripping_He)*1.d-4 ! convert cm2 to m2

      return
      end

      function sigma_ionization_H0_He(Ek)
c     ==============================================================
c     VERSION:         0.1
c     LAST MOD:       OCT/06
c     MOD AUTHOR:    G. Fubiani
c     COMMENTS:      Calculate cross section for ionization of H2 (D2). 
c                    H0 + He -> H0 + He+ + e- ( D0 + He -> D0 + He+ + e-)
c                    DuBois PRA 40, 3605 (1989) Fig 3(a) (energy range 25 keV to 1 MeV)
c     NOTE:          (1) Fit log(sigma) = sum^5_{i=1} a_i log(Ek)**i was
c                    obtained using gnuplot.
c                    Ek must be in keV, sigma originally in cm2
c                    (2) For Deuterium we assume that a velocity identical to
c                    the one of an Hydrogen atom gives the same cross 
c                    section, i.e., in terms of kinetic energies, for a 
c                    given Deuterium energy, the cross section is identical
c                    to Hydrogen with half that energy. 
c                    (3) For ionization with H- (D-) we use the same cross
c                    section as for H0 (D0)
c     --------------------------------------------------------------
      implicit none
      include 'particle_info.h'
      integer i,ntot 
      parameter (ntot=5)
      double precision a(ntot),sigma_ionization_H0_He,Ek

      if(ptype.ge.5) Ek=Ek/2.d0 ! For Deuterium Ek_H == Ek_D/2

      a(1)=-42.458d0
      a(2)=2.74073d0
      a(3)=-0.323338d0
      a(4)=-0.026335d0
      a(5)=0.0036327d0

      sigma_ionization_H0_He=0.d0

      do i=1,ntot
         sigma_ionization_H0_He = sigma_ionization_H0_He + 
     &        a(i)*dlog(Ek)**(i-1)
      enddo

      sigma_ionization_H0_He=dexp(sigma_ionization_H0_He)*1.d-4 ! convert cm2 to m2
      if(Ek.le.1.d-2) sigma_ionization_H0_He=0.d0 ! Ek < 0.01 keV 

      return
      end

      function sigma_ioni_H0_loss_He(Ek)
c     ==============================================================
c     VERSION:         0.1
c     LAST MOD:       OCT/06
c     MOD AUTHOR:    G. Fubiani
c     COMMENTS:      Calculate cross section for ionization of H2 (D2). 
c                    H0 + He -> H+ + He+ + 2e- ( D0 + He -> D+ + He+ + 2e-)
c                    DuBois PRA 40, 3605 (1989) Fig 3(c) (energy range 25 keV to 1 MeV)
c     NOTE:          (1) Fit log(sigma) = sum^6_{i=1} a_i log(Ek)**i was
c                    obtained using gnuplot.
c                    Ek must be in keV, sigma originally in cm2
c                    (2) For Deuterium we assume that a velocity identical to
c                    the one of an Hydrogen atom gives the same cross 
c                    section, i.e., in terms of kinetic energies, for a 
c                    given Deuterium energy, the cross section is identical
c                    to Hydrogen with half that energy. 
c                    (3) For ionization with H- (D-) we use the same cross
c                    section as for H0 (D0)
c     --------------------------------------------------------------
      implicit none
      include 'particle_info.h'
      integer i,ntot 
      parameter (ntot=6)
      double precision a(ntot),sigma_ioni_H0_loss_He,Ek

      if(ptype.ge.5) Ek=Ek/2.d0 ! For Deuterium Ek_H == Ek_D/2

      a(1)=-81.6381d0
      a(2)=25.9177d0
      a(3)=-4.30974d0
      a(4)=-0.172419d0
      a(5)=0.104279d0
      a(6)=-0.00713259d0 

      sigma_ioni_H0_loss_He=0.d0

      do i=1,ntot
         sigma_ioni_H0_loss_He = sigma_ioni_H0_loss_He + 
     &        a(i)*dlog(Ek)**(i-1)
      enddo

      sigma_ioni_H0_loss_He=dexp(sigma_ioni_H0_loss_He)*1.d-4 ! convert cm2 to m2
      if(Ek.le.1.d-2) sigma_ioni_H0_loss_He=0.d0 ! Ek < 0.01 keV 

      return
      end

      function sigma_stripping_H0_He(Ek)
c     ==============================================================
c     VERSION:         0.1
c     LAST MOD:       OCT/06
c     MOD AUTHOR:    G. Fubiani
c     COMMENTS:      Calculate cross section for single neutral stripping. 
c                    H + He -> H+ + He + e- ( D + He -> D+ + He + e-)
c                    Redbook data Vol 1 E11 (energy range 0.05 keV to 20 MeV)
c     NOTE:          (1) Fit log(sigma) = sum^5_{i=1} a_i log(Ek)**i was
c                    obtained using gnuplot.
c                    Ek must be in keV, sigma originally in cm2
c                    (2) For Deuterium we assume that a velocity identical to
c                    the one of an Hydrogen atom gives the same cross 
c                    section, i.e., in terms of kinetic energies, for a 
c                    given Deuterium energy, the cross section is identical
c                    to Hydrogen with half that energy. 
c     --------------------------------------------------------------
      implicit none
      include 'particle_info.h'
      integer i,ntot 
      parameter (ntot=5)
      double precision a(ntot),sigma_stripping_H0_He,Ek

      if(ptype.ge.5) Ek=Ek/2.d0 ! For Deuterium Ek_H == Ek_D/2

      a(1)=-37.9671d0
      a(2)=1.06995d0
      a(3)=-0.194566d0
      a(4)=-0.000891849d0
      a(5)=0.00058525d0

      sigma_stripping_H0_He = 0.d0

      do i=1,ntot
         sigma_stripping_H0_He = sigma_stripping_H0_He + 
     &        a(i)*dlog(Ek)**(i-1)
      enddo

      sigma_stripping_H0_He = dexp(sigma_stripping_H0_He)*1.d-4 ! convert cm2 to m2

      return
      end
