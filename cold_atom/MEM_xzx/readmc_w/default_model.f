       program default_model
c      for generate default model
c      in a Gaussian Form
       
       implicit none
       integer*4 i,nf,nk,L
       double precision wmin,wmax,stepw,tk
       double precision PI,w,m,gamma,f,J,wl
       double precision muller_func
       parameter(PI = 3.1415926535897932453626D0)
       parameter(gamma=0.2d0)
       parameter(nf=2000)
       parameter(L=128)
       J    = 1d0
       wmin = 0.d0
       wmax = PI*J*8.d0
       stepw = (wmax-wmin) / DBLE(nf)
c       write(*,*) 'In default nk=?'
c       read(*,9) nk
c9      format(I2)
c       write(*,*) 'wmax=',wmax,' stepw=',stepw
c       pause
       w = 0.d0
       tk = 2d0*PI/DBLE(L)*DBLE(nk+1)
c       wl = pi*J*DSIN(tk)/2d0
       wl = 0.4105*J
       OPEN(9,FILE='deflt.dat',STATUS='unknown')
       DO i = 1, 560
          m = DEXP(-((w-wl)/gamma)*((w-wl)/gamma))/DSQRT(PI)/gamma
c           m = 1d0/PI/J/2.d0
c          WRITE(*,*) w,m,stepw
          WRITE(9,*) w,m,stepw
          w = w+stepw
c          pause
       ENDDO
       CLOSE(9)

       OPEN(10,FILE='f0.dat',STATUS='unknown')
       w = 0.d0
       DO i = 1, 560
c          f = DEXP(-w)/(1d0-DEXP(-PI*J))
c          f = muller_func(tk,w,J)+0.15d0
          f = DEXP(-((w-wl)/gamma)*((w-wl)/gamma))/DSQRT(PI)/gamma
c          f = 1.d0/PI/J/2.d0
c          WRITE(*,*) w,f
          WRITE(10,*) w,f,stepw
          w = w+stepw
c          pause
       ENDDO
       CLOSE(10)

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

