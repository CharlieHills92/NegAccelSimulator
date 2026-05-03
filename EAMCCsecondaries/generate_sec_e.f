      subroutine generate_sec_e(etaSEC,Elost,Ysec,iseed,u,t,n,x_g,y_g,
     &     z_g,lsec_max)
c     ==============================================================
c     VERSION:         0.1
c     LAST MOD:       FEB/07
c     MOD AUTHOR:    G. Fubiani
c     COMMENTS:      Generate true secondary electrons after
c                    particle impact on grids
c     NOTE:          Secondary electrons are assumed to be emitted at low 
c                    energy. For now we suppose Eout = 10. eV, following 
c                    H. de Esch assumption.
C                    For ions see for instance ORNL REdbook Vol 3, C17
c     --------------------------------------------------------------
      implicit none
      include 'particle_info.h'
      include 'ode_arrays.h'
      integer iseed,ptype_info,cnt_lsec,lsec_max
      double precision Ysec(nmax_sec),u(3),t(3),n(3),x_g,y_g,z_g
      double precision etaSEC,rnd,ran2,Eout,Elost

C      open(68,file='./DATA/secondary_emission.dat',form='formatted',
C     >     access='append') !added (2013/01)

      rnd=ran2(iseed) ! Start checking for a single secondary electron
      !rnd=0    !nf da togliere
      if(rnd.le.etaSEC.and.Elost.gt.10.d0) then 

         Eout=10.d0
         Elost=Elost-Eout       ! Energy lost by incoming particle (eV)

         call checknetcharge(-charge(1)) ! Cnt net charge absorbed by grid

C         write(68,*) grid_ind, ptype, charge(ptype),-1 !added (2013/01)
C         numsec = numsec + 1
         if(flag_info.eq.1)then
            write(6,40) Eout
 40         format('Secondary electron emission. Eout=',f6.2,' eV') 
         endif

         n_sec=n_sec+1
         cnt_sec(1) = cnt_sec(1) + 1 ! Counter for sec. elec. 
         flag_sec=1
         ptype_info=1
         call make_sec_wall_part(Ysec,Eout,iseed,u,t,n,x_g,y_g,z_g,
     &        ptype_info)

C
C     More than one secondary electron emmited. 
C 
         if( etaSEC.gt.1.d0 ) then ! Start checking for more than 1 sec. elec.

            etaSEC=etaSEC-1.d0
            cnt_lsec=1

            do while( etaSEC.gt.0.d0.and.Elost.gt.10.d0 )
               
               rnd=ran2(iseed)
               if( rnd.le.etaSEC ) then

                  call checknetcharge(-charge(1)) ! Cnt net charge absorbed by grid  
C                  write(68,*) grid_ind, ptype, charge(ptype), -1 !added (2013/01)

                  Eout=10.d0
                  Elost=Elost-Eout ! Energy lost by incoming particle (eV)
                  
                  if(flag_info.eq.1) then
                     write(6,50) cnt_lsec+1,Eout
 50                  format('Emission of',1x,i2,1x,
     &                    'secondary electrons. Eout=',f6.2,' eV') 
                  endif

                  n_sec=n_sec+1 ! Total # of sec. part. generated
                  cnt_sec(1) = cnt_sec(1) + 1 ! Counter for sec. elec.
                  call make_sec_wall_part(Ysec,Eout,iseed,u,t,n,
     &                 x_g,y_g,z_g,ptype_info)

                  cnt_lsec=cnt_lsec+1 ! Update cnt_lsec
                  nmax_nsec=MAX(nmax_nsec,cnt_lsec)

                  if(cnt_lsec.gt.lsec_max) then ! WARNING
                     write(*,60) cnt_lsec
 60                  format(' WARNING: cnt_lsec reached max. value',
     &                    1x,i2)
                     flag_err=1
                  endif

               endif

               etaSEC=etaSEC-1.d0 ! update etaSEC

            enddo

         endif ! End checking for more than 1 sec. electron
         
      endif  ! End checking for a single secondary electron
      
      close(68)

      return
      end
