        program energy_surface
c*******************************************************************************
c       to compile type make spectra (or make all) using the 
c       associated Makefile.  This code was designed to read the
c       magority of its inputs from a qmc sigma file.
c*******************************************************************************
        use Global
        use module_spectra
        implicit none
c*******************************************************************************
        integer :: filenum, istat, i, j, n, run, meas, ick,ios, 
     &             nkpoints, info, infot, ickx, icky, i_energy

        real(kind) :: r1, r2, kx, ky,  Epsk, kx1, ky1, kx2, 
     &                ky2, offset, Energy
        
        complex(kind) :: c1, c2, blint
        
        character*72 linemc,fmt1,string1
c*******************************************************************************


c*******************************************************************************
c       Read in some parameters
c       Open a datafile for the parameters
c*******************************************************************************
        write(6,2) 
 2      format('       enter the QMC sigma filenumber:  ',$)
        read(5,*) filenum


        if(filenum.lt.10) then  ! P. 402, Chapman
          fmt1="('sigma',i1,'.dat')"
        else if(filenum.lt.100) then
          fmt1="('sigma',i2,'.dat')"
        else if(filenum.lt.1000) then
          fmt1="('sigma',i3,'.dat')"
        else if(filenum.lt.10000) then
          fmt1="('sigma',i4,'.dat')"
        else
          write(6,*) 'filenum to large'
          stop
        end if  
        write(string1,fmt1) filenum
                
c       Now open the QMC self energy 
        open(unit=20,file=string1,status='old')


 1000   format(a72)
c       Read some parameters from the QMC sigma file
 401    read(20,1000) linemc
        if(index(linemc(1:72),'Results').ne.0 ) then
          read(20,*) nl,nwn,fill,cluster, ntot,ndim
        else
         goto 401
        endif
        read(cluster,"(i3)",iostat=ios) Nc ! should barf with characters
          if (ios.gt.0) then
            read(cluster,"(i2)",iostat=ios) Nc   ! ignore trailing charac.
            if (ios.gt.0) then
              read(cluster,"(i1)",iostat=ios) Nc ! ignore trailing charac.
              if (ios.gt.0) then
                write(6,*) "Cannot read Nc from string ",cluster
                stop                   
              endif
            endif
          endif

402    read(20,"(a72)") linemc
        if(index(linemc(1:72),'tprime').ne.0 ) then
           if(index(linemc(1:72),'tdprime').ne.0 ) then
              if(index(linemc(1:72),'t3prime').ne.0 ) then
                 read(20,*) ed,U,beta,i,run,meas,tprime,tdprime,t3prime
              else
                 read(20,*) ed,U,beta,i,run,meas,tprime,tdprime
                 t3prime=0.d0
              endif
           else 
              read(20,*) ed,U,beta,i,run,meas,tprime
              tdprime=0.d0
              t3prime=0.d0
           endif
        else
          goto 402
        endif
        temp=1.0/beta   



c*******************************************************************************        
c       Form the symmetry tables and allocate the arrays.
c*******************************************************************************        
        call allocate_1p
        call tables
        nf=1000  ! for declaration only.  It will be reset below
        call allocate_Spec0         
        call allocate_Spec1         

        
        
        do ick=1,Ncw  ! do all K in the irreducible wedge
c         Open some files.
          if(ick.lt.10) then
            if(filenum.lt.10) then
              fmt1="('Sigma',i1,'_',i1,'.dat')"
            else if(filenum.lt.100) then
              fmt1="('Sigma',i2,'_',i1,'.dat')"
            else if(filenum.lt.1000) then
              fmt1="('Sigma',i3,'_',i1,'.dat')"
            else if(filenum.lt.10000) then
              fmt1="('Sigma',i4,'_',i1,'.dat')"
            end if
          else if(ick.lt.100) then
            if(filenum.lt.10) then
              fmt1="('Sigma',i1,'_',i2,'.dat')"
            else if(filenum.lt.100) then
              fmt1="('Sigma',i2,'_',i2,'.dat')"
            else if(filenum.lt.1000) then
              fmt1="('Sigma',i3,'_',i2,'.dat')"
            else if(filenum.lt.10000) then
              fmt1="('Sigma',i4,'_',i2,'.dat')"
            end if
          else
            write(6,*) 'ick too large'
            stop
          end if
            
          write(string1,fmt1) filenum,ick
          open(unit=21,file=string1,status='old')
          
          read(21,1000) linemc    
          do i=1,1000
            read(21,*,IOSTAT=istat) w(i),r1,r2
            if(istat.ne.0) goto 998
            Sigma(i,ick)=cmplx(r1,r2)
          end do
 998      continue
          if(ick.eq.1) then
            nf=i-1
          else if(nf.ne.i-1) then       ! trap for errors
            write(6,*) 'Sigma files are not compatible'
            stop
          end if
        end do 
        
        if(nf.gt.1000) then
          write(6,*) 'nf in the data file is too large'
          stop
        end if
        
c       Now fill in Sigma(w,k), for k outside the IW.
        do ick=Ncw+1,Nc
        do i=1,nf
          Sigma(i,ick)=Sigma(i,ickmap(ick))
        end do
        end do

        if(ispline.eq.1) then
c          call allocate_Spec1
c         Form the 2D interpolation arrays needed to interpolate
c         the self energy with a bicubic spline.
c         1. form the grid in the space of g1 and g2
          do i=-N_sp,N_sp
            g1_sp(i)=i
            g2_sp(i)=i
          end do

          do n=1,nf
            do i=-N_sp,N_sp
            do j=-N_sp,N_sp
              ick=ict(i,j,0)
              Sig_sp_r(i,j,n)=real(Sigma(n,ick))
              Sig_sp_i(i,j,n)=aimag(Sigma(n,ick))
            end do
            end do
            call splie2(g1_sp(-N_sp),g2_sp(-N_sp),
     &                  Sig_sp_r(-N_sp,-N_sp,n),
     &                  2*N_sp+1,2*N_sp+1,
     &                  Derivs_sp_r(-N_sp,-N_sp,n))
            call splie2(g1_sp(-N_sp),g2_sp(-N_sp),
     &                  Sig_sp_i(-N_sp,-N_sp,n),
     &                  2*N_sp+1,2*N_sp+1,
     &                  Derivs_sp_i(-N_sp,-N_sp,n))
          end do
        end if

c       Now do the interpolation, and calculate the spectrum at
c       a sequence of locations.

        write(6,3) 
 3      format('Enter Energy,nkpoints:  ',$)
        read(5,*) Energy,nkpoints
        do i=1,nf
          if(w(i).ge.Energy) exit
        end do
        i_energy=i
        
        open(unit=22,file='energy_surface.dat',status='unknown')
        open(unit=21,file='energy_surface0.dat',status='unknown')
        write(6,*) 'output in energy_surface.dat'
        write(6,*) 'gnuplot output in energy_surface0.dat'
        write(22,"('#  kx            ky            w      ',   
     &   '    ReAkw        ImAkw       ',
     & ' ReSigkw       ImSigkw')")
        write(21,"('#  kx    ky     -1/pi ImAkw       ')")
        kx1=0.0
        ky1=0.0
        kx2=3.1415927
        ky2=3.1415927
        do ickx=1,nkpoints
        do icky=1,nkpoints
          kx=kx1 + real(ickx-1)*(kx2-kx1)/real(nkpoints-1,kind)
          ky=ky1 + real(icky-1)*(ky2-ky1)/real(nkpoints-1,kind)
c          c2=blint(kx,ky,i_energy)                      ! Sigma(k,w)
          c2=0.5*(blint(kx,ky,i_energy)+blint(ky,kx,i_energy)) 
          c1=1.0/(w(i_energy)-ed-Epsk(kx,ky)-c2)  ! G(k,w) no ed => none in mem
          write(22,"(7(g12.6,1x))") kx,ky,w(i_energy),real(c1),
     &                              aimag(c1),real(c2),aimag(c2)
          write(21,"(3(g12.6,1x))") kx,ky,-aimag(c1)/pi
        end do
          write(22,*) ' '
          write(21,*) ' '
        end do

c*******************************************************************************
c       Wrap things up
c*******************************************************************************
        call deallocate_1p
        call deallocate_Spec0         
        if(ispline.eq.1) call deallocate_Spec1        
        
        stop
        end
                  
        
        
        complex(8) function blint(kx,ky,n)
c       
c       This subroutine performs a bilinear interpolation of
c       the self energy, at any location (kx,ky) in the zone.
c
        use Global
        use module_spectra
        implicit none
        integer :: n,i1,j1,i2,j2,ic1,ic2,ic3,ic4
        real(kind) :: Kx,Ky,dx,dy,Kx3,Ky3,g1,g2,r7,r8
        complex(kind) :: z1,z2,z3,z4
        
        if(ispline.eq.0) then
c         Use a two-dimensional bilinear interpolation within the square
c
c
c         1*        *2
c
c
c
c
c         3*        *4
c
c

c         First find the coordinates of the cell which contains (kx,ky)
          i1=int((kx*g(1,1)+ky*g(1,2))/(g(1,1)*g(1,1)+g(1,2)*g(1,2))+Nc)-Nc
          i2=i1+1
          j1=int((kx*g(2,1)+ky*g(2,2))/(g(2,1)*g(2,1)+g(2,2)*g(2,2))+Nc)-Nc
          j2=j1+1
c         now get the K-point on the lower left.
          Kx3=i1*g(1,1)+j1*g(2,1)
          Ky3=i1*g(1,2)+j1*g(2,2)
c         index the four corners
          ic1=ict(i1,j2,0)
          ic2=ict(i2,j2,0)
          ic3=ict(i1,j1,0)
          ic4=ict(i2,j1,0)
          z1=Sigma(n,ic1)
          z2=Sigma(n,ic2)
          z3=Sigma(n,ic3)
          z4=Sigma(n,ic4)
c         calculate the slopes
          dx=((kx-Kx3)*g(1,1)
     &      +(ky-Ky3)*g(1,2))/(g(1,1)*g(1,1)+g(1,2)*g(1,2))
          dy=((kx-Kx3)*g(2,1)
     &      +(ky-Ky3)*g(2,2))/(g(2,1)*g(2,1)+g(2,2)*g(2,2))
c         interpolate.
          blint=z3+(z4-z3)*dx+(z1-z3)*dy+(z2-z1-z4+z3)*dx*dy
        else if(ispline.eq.1) then
c         Use a bicubic spline.
c         Find the location on the grid of principle translation
c         vectors in K-space.
          g1=(kx*g(2,2)-ky*g(2,1))/(g(1,1)*g(2,2)-g(1,2)*g(2,1))
          g2=(kx*g(1,2)-ky*g(1,1))/(g(2,1)*g(1,2)-g(2,2)*g(1,1))
c         Interpolate for the real part of Sigma(kx,ky,n)
          call splin2(g1_sp(-N_sp),g2_sp(-N_sp),
     &                Sig_sp_r(-N_sp,-N_sp,n),
     &                Derivs_sp_r(-N_sp,-N_sp,n),
     &                2*N_sp+1,2*N_sp+1,
     &                g1,g2,r7)
c         Interpolate for the imaginary part of Sigma(kx,ky,n)
          call splin2(g1_sp(-N_sp),g2_sp(-N_sp),
     &                Sig_sp_i(-N_sp,-N_sp,n),
     &                Derivs_sp_i(-N_sp,-N_sp,n),
     &                2*N_sp+1,2*N_sp+1,
     &                g1,g2,r8)
     
          blint=cmplx(r7,r8)
        else
          write(6,*) 'BAD ispline'
          stop
        end if
        return
        end
                
        
        real*8 function Epsk(kx,ky)
        use Global
        implicit none
        real(kind) :: kx, ky    
  

            Epsk=-0.5_kind*(cos(kx)+cos(ky)) -
     &       tprime*(cos(kx)*cos(ky)-1.0_kind)
     &             - 0.5*tdprime*(cos(2.d0*kx)+cos(2.d0*ky))
     &             - t3prime*(cos(2.d0*kx)*cos(ky)
     &                         +cos(kx)*cos(2.d0*ky))  
        end function Epsk

c       NUMERICAL RECIPES 2D SPLINE

      SUBROUTINE SPLIE2(X1A,X2A,YA,M,N,Y2A) 
      PARAMETER (NN=100) 
      IMPLICIT real(8) (A - H, O - Z) 
      real(8) X1A(M),X2A(N),YA(M,N),Y2A(M,N),YTMP(NN),
     &           Y2TMP(NN) 
      DO 13 J=1,M 
        DO 11 K=1,N 
          YTMP(K)=YA(J,K) 
11      CONTINUE 
        CALL SPLINE(X2A,YTMP,N,1.d30,1.d30,Y2TMP) 
        DO 12 K=1,N 
          Y2A(J,K)=Y2TMP(K) 
12      CONTINUE 
13    CONTINUE 
      RETURN 
      END 


      SUBROUTINE SPLIN2(X1A,X2A,YA,Y2A,M,N,X1,X2,Y) 
      PARAMETER (NN=100) 
      IMPLICIT real(8) (A - H, O - Z) 
      real(8) X1A(M),X2A(N),YA(M,N),Y2A(M,N),YTMP(NN),Y2TMP(NN),
     &        YYTMP(NN) 
      DO 12 J=1,M 
        DO 11 K=1,N 
          YTMP(K)=YA(J,K) 
          Y2TMP(K)=Y2A(J,K) 
11      CONTINUE 
        CALL SPLINT(X2A,YTMP,Y2TMP,N,X2,YYTMP(J)) 
12    CONTINUE 
      CALL SPLINE(X1A,YYTMP,M,1.d30,1.d30,Y2TMP) 
      CALL SPLINT(X1A,YYTMP,Y2TMP,M,X1,Y) 
      RETURN 
      END 

      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2) 
      IMPLICIT real(8) (A - H, O - Z) 
      PARAMETER (NMAX=100) 
      real(8) X(N),Y(N),Y2(N),U(NMAX) 
      IF (YP1.GT..99E30) THEN 
        Y2(1)=0. 
        U(1)=0. 
      ELSE 
        Y2(1)=-0.5 
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1) 
      ENDIF 
      DO 11 I=2,N-1 
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1)) 
        P=SIG*Y2(I-1)+2. 
        Y2(I)=(SIG-1.)/P 
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1)) 
     *      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P 
11    CONTINUE 
      IF (YPN.GT..99E30) THEN 
        QN=0. 
        UN=0. 
      ELSE 
        QN=0.5 
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1))) 
      ENDIF 
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.) 
      DO 12 K=N-1,1,-1 
        Y2(K)=Y2(K)*Y2(K+1)+U(K) 
12    CONTINUE 
      RETURN 
      END 


      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y) 
      IMPLICIT real(8) (A - H, O - Z) 
      real(8) XA(N),YA(N),Y2A(N) 
      KLO=1 
      KHI=N 
1     IF (KHI-KLO.GT.1) THEN 
        K=(KHI+KLO)/2 
        IF(XA(K).GT.X)THEN 
          KHI=K 
        ELSE 
          KLO=K 
        ENDIF 
      GOTO 1 
      ENDIF 
      H=XA(KHI)-XA(KLO) 
      IF (H.EQ.0.) PAUSE 'Bad XA input.' 
      A=(XA(KHI)-X)/H 
      B=(X-XA(KLO))/H 
      Y=A*YA(KLO)+B*YA(KHI)+ 
     *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6. 
      RETURN 
      END 


