       program default_model
c      for generate default model
c      in a Gaussian Form
       
       implicit none
       integer*4 i,nf,nk,L
       double precision wmin,wmax,stepw,tk
       double precision PI,gamma,f,J,wl,w_tmp
       double precision muller_func
       double precision,pointer :: w(:),model(:)
       double precision normaler
       parameter(PI = 3.1415926535897932453626D0)
       parameter(gamma=0.1d0)
       parameter(nf=10000)
       parameter(L=128)
       allocate(w(nf))
       allocate(model(nf))
       J    = 1d0
       wmin = 0.d0
       wmax = PI*J*10.d0
       stepw = (wmax-wmin)/DBLE(nf)
!
!       write(*,*) 'In default nk=?'
!       read(*,9) nk
!9      format(I2)
       write(*,*) 'wmax=',wmax,' stepw=',stepw
       tk = 2d0*PI/DBLE(L)*DBLE(nk)
c****approximated S=1 HAF dispersion relation
c       wl = DSQRT(0.4105**2+(2.d0*J)**2*(tk-PI)**2)
c       wl = DSQRT(0.4105**2 + (2.d0*J)**2*DSIN(tk)**2)
       wl = 0.41d0
c      
       normaler = 0.d0 
       w_tmp = 0.d0
       DO i = 1, 1000
          w(i) = w_tmp
          model(i) = DEXP(-((w(i)-wl)/gamma)*((w(i)-wl)/gamma))
     &/DSQRT(PI)/gamma
          if(model(i).le.1d-6) model(i)=1.0d-6
c          model(i) = 1d0/PI/J/2.d0
          w_tmp = w_tmp + stepw
          normaler = normaler + model(i)*stepw
       ENDDO

       OPEN(9 ,FILE='deflt.dat', STATUS='unknown')
       OPEN(10,FILE='f0.dat'   , STATUS='unknown')
       write(10,*) '#normalor:',tk
       write(9,*)  '#normalor:',tk
       DO i = 1, 560
          write(9, *) w(i),model(i)/normaler,stepw
          WRITE(10,*) w(i),model(i)/normaler,stepw
       ENDDO
       CLOSE(10)
       CLOSE(9)

       deallocate(w,model)

       END

      function muller_func(tk,w,jb)
      implicit none
      real*8 muller_func
      real*8 stepfunc,factor
      real*8 tk,w,jb,pi
      real*8 wl,wu
      pi = 3.141592653589793d0

      wl = pi*jb*DSIN(tk)/2d0
      wu = pi*jb*DSIN(tk/2d0)

      factor    =1d0/DSQRT(DABS(w*w-wl*wl))/pi
      muller_func=stepfunc(w,wl)*stepfunc(wu,w)*factor
c      write(*,11) muller_fun
c11    format(1x,f10.5)
      return
      end

      function stepfunc(x,y)
      implicit none
      real*8 stepfunc
      real*8 x, y

      if(x.gt.y) then
        stepfunc= 1d0
      else if(x.eq.y) then
        stepfunc=0.5d0
      else
        stepfunc=0d0
      endif

      return
      end

