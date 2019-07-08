       program default_model
c      for generate default model
c      in a Gaussian Form
       
       implicit none
       integer*4 i,nf,nk,L
       double precision wmin,wmax
       double precision PI,gamma,f,J,wl,w_tmp
       double precision muller_func
       double precision,pointer :: w(:),model(:),stepw(:)
       double precision normaler
       parameter(PI = 3.1415926535897932453626D0)
       parameter(gamma=2d0)
       parameter(nf=1120)
       parameter(L=32)
       allocate(w(nf))
       allocate(model(nf))
       allocate(stepw(nf))
       J    = 1d0
       wmin = 0.d0
       wmax = PI*J*5.d0
!       stepw = (wmax-wmin)/DBLE(nf)
!
       wl = 0.1d0
c      
       normaler = 0.d0 
       w_tmp = 0.d0
       DO i = 1, nf
          w(i) = w_tmp
          model(i) = DEXP(-(w(i)-wl)**2/gamma**2)
     &/DSQRT(PI)/gamma
!          stepw(i) = 0.25d-6*DBLE(i)**2
          stepw = (wmax-wmin)/DBLE(nf)
          w_tmp = w_tmp + stepw(i)
          normaler = normaler + model(i)*stepw(i)
       ENDDO

       OPEN(9 ,FILE='deflt.dat', STATUS='unknown')
       OPEN(10,FILE='f0.dat'   , STATUS='unknown')
       write(10,*) '#normalor:',normaler
       write(9,*)  '#normalor:',normaler
       DO i = 1, nf
          write(9, *) w(i),model(i)/normaler,stepw(i)
          WRITE(10,*) w(i),model(i)/normaler,stepw(i)
       ENDDO
       CLOSE(10)
       CLOSE(9)

       deallocate(w,model,stepw)

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

