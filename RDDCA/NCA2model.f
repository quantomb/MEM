        program NCA2model
c*******************************************************************************
c       This block is designed to take the NCA-DCA Hubbard model single-particle
c       A_bar(K,w) and generate A_bar(K,w) on a Lorentzian grid, suitable
c       as a MEM default model.
c*******************************************************************************

        implicit none
        integer, parameter :: nf=340, nwmax=1000
        real(8), parameter :: wmax=8.0_8, delta=0.5_8, pi=3.1515927_8

        integer :: nw, i, n
        real(8) :: A_nca(nwmax),w_nca(nwmax),A_model(nf),w(nf),dw(nf),
     &             scoef(nwmax), y, r1

        character*72 linemc
c*******************************************************************************
        
c       Open a datafile for the input spectra
        write(6,1) 
 1      format('enter the NCA spectra filename:  ',$)
        read(5,1000) linemc
        open(unit=9,file=linemc,status='old')

c       Open a datafile for the parameters
        write(6,2) 
 2      format('enter the MEM sepctra filename:  ',$)
        read(5,1000) linemc
        open(unit=10,file=linemc,status='unknown')
 1000   format(a72)
 
        
c       Read in the NCA spectra A_nca(w) and the frequency grid
        do i=1,1000
          read(9,*,end=15) w_nca(i),A_nca(i)
        end do
 15     continue
        nw=i-1


c       Now Akima spline the spectra
        call spline(w_nca,A_nca,nw,0.0_8,0.0_8,scoef)

c       To evaluate at w(l)
c       call splint(w_nca,A_nca,scoef,nw,w(l),A_model(l))       

        do i=-nf/2+1,nf/2
          n=i+nf/2
          y=(i-0.5_8)/real(nf)
          w(n)=delta*tan(pi*y)
          if(abs(w(n)).le.wmax) then
            dw(n)=pi*(delta**2+(w(n))**2)/(delta*nf)
            call splint(w_nca,A_nca,scoef,nw,w(n),r1)
            write(10,*) w(n),r1,dw(n)
          end if
        end do

        stop
        end program NCA2model
        
        
        
        SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
        implicit none
        integer, parameter :: NMAX=1000
        integer :: N,I,K
        real(8) :: YP1,YPN,SIG,P,QN,UN
        real(8) :: X(N),Y(N),Y2(N),Uspline(NMAX)
        IF (YP1.GT..99E30) THEN
          Y2(1)=0.0
          Uspline(1)=0.0
        ELSE
          Y2(1)=-0.5
          Uspline(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
        ENDIF
        DO 11 I=2,N-1
          SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
          P=SIG*Y2(I-1)+2.0
          Y2(I)=(SIG-1.)/P
          Uspline(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     &      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*Uspline(I-1))/P
11      CONTINUE
        IF (YPN.GT..99E30) THEN
          QN=0.0
          UN=0.0
        ELSE
          QN=0.5
          UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
        ENDIF
        Y2(N)=(UN-QN*Uspline(N-1))/(QN*Y2(N-1)+1.)
        DO 12 K=N-1,1,-1
          Y2(K)=Y2(K)*Y2(K+1)+Uspline(K)
12      CONTINUE
        RETURN
        END


        SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
        implicit none
        integer :: N,K,KLO,KHI
        real(8) :: XA(N),YA(N),Y2A(N),X,Y,H,A,B
        KLO=1
        KHI=N
1       IF (KHI-KLO.GT.1) THEN
          K=(KHI+KLO)/2
          IF(XA(K).GT.X)THEN
            KHI=K
          ELSE
            KLO=K
          ENDIF
        GOTO 1
        ENDIF
        H=XA(KHI)-XA(KLO)
        IF (H.EQ.0.) THEN
          write(6,*) 'Bad XA input, Aborting...'
          CALL EXIT()
        ENDIF
        A=(XA(KHI)-X)/H
        B=(X-XA(KLO))/H
        Y=A*YA(KLO)+B*YA(KHI)+
     &    ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
        RETURN
        END
