       program default_model
c
c      generate default model for density of state
c      in a Gaussian Form 
c      N_k=1 for Fermi
c      N_k=0 for Boson
c      N_k=2 for Boson + Fermi
c       
       implicit none
       integer*4 i,nf,N_k
       double precision wmin,wmax,stepw,tk
       double precision PI,w,m,gamma,f,J,wl
       double precision muller_func
       parameter(PI = 3.1415926535897932453626D0)
       parameter(gamma=1.2d0)
       parameter(nf=2000)
c
c      range of frequency (wmin,wmax)=(0,31.4)
       wmin = 0.d0
       wmax = PI*10.d0
c      resolution of frequency
       stepw = (wmax-wmin) / DBLE(nf)
       write(*,*) 'wmax=',wmax,' stepw=',stepw
c
       write(*,*) 'Choose default model N_k=?'
       read(*,9) N_k
 9     format(I2)
c
c      specify central value for Gaussian form
       wl = 0.d0
c      normalization factor tk
       tk = 1.d0
c       
       w = stepw*DBLE(560)*(-1.d0)
       OPEN(9,FILE='deflt.dat',STATUS='unknown')
       write(9,*) '#normalor:',tk
       DO i = -560, 560
          m = DEXP(-((w-wl)/gamma)*((w-wl)/gamma))/DSQRT(PI)/gamma
c           m = 1.d0/PI/2.d0
c          WRITE(*,*) w,m,stepw
          WRITE(9,*) w,m,stepw
          w = w+stepw
       ENDDO
       CLOSE(9)

       OPEN(10,FILE='f0.dat',STATUS='unknown')
       write(10,*) '#normalor:',tk
       w = stepw*DBLE(560)*(-1.d0)
       DO i = -560, 560
          f = DEXP(-((w-wl)/gamma)*((w-wl)/gamma))/DSQRT(PI)/gamma
c          f = 1.d0/PI/2.d0
c          WRITE(*,*) w,f
          WRITE(10,*) w,f,stepw
          w = w+stepw
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

