        program default
c******************************************************************************
c       to compile type make 2default (or make all) using the 
c       associated Makefile.  This code was designed to read the
c       magority of its inputs from a qmc sigma file.  It uses
c       the coarse-grained A(K,w) to generate a two-particle 
c       default model.
c******************************************************************************
c******************************************************************************
        use Global
        use module_spectra
        implicit none
c******************************************************************************
c       some local declarations 
        integer :: filenum, istat, i, j, n, run, meas, ick, jck,
     &             ik, isk, nkpoints, info, infot, ikx, iky, deg, nfb

        real(kind) :: r1, r2, r3, kx, ky, fill, kx1, ky1, kx2, 
     &                ky2, offset, res, hall, vx, vyy, fermi, delta,
     &                y
        
        complex(kind) :: c1, c2
        
        character*72 linemc,fmt1,string1
c******************************************************************************


c******************************************************************************
c       Read in some parameters
c******************************************************************************
c       Open a datafile for the parameters
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
                
c       Now open the QMC self energy sigma.dat
        open(unit=20,file=string1,status='old')

 1000   format(a72)
c       Read some parameters from the QMC sigma file
 401    read(20,1000) linemc
        if(index(linemc(1:72),'Results').ne.0 ) then
          read(20,*) nl,nwn,fill,Nc,nover
        else
         goto 401
        endif

 402    read(20,1000) linemc
        if(index(linemc(1:72),'tprime').ne.0 ) then
          read(20,*) ed,U,beta,i,run,meas,tprime
        else
         goto 402
        endif
        temp=1.0/beta   
        if(nover.le.40) nover=140
        

c******************************************************************************
c       This is all the information we need to allocate the arrays
c******************************************************************************
        call allocate_1p
        nf=1000      !set preliminary value just for allocation
        call allocate_Spec0

c******************************************************************************
c       calculate the symmetry tables
c******************************************************************************
        call tables
        
c******************************************************************************
c       Open the MEM dos files, readin the spectra, and then spline it.
c******************************************************************************
        do ick=1,Ncw  ! do all K in the irreducible wedge
c         Open some files.
          if(ick.lt.10) then
            if(filenum.lt.10) then
              fmt1="('dos',i1,'_',i1,'QMC')"
            else if(filenum.lt.100) then
              fmt1="('dos',i2,'_',i1,'QMC')"
            else if(filenum.lt.1000) then
              fmt1="('dos',i3,'_',i1,'QMC')"
            else if(filenum.lt.10000) then
              fmt1="('dos',i4,'_',i1,'QMC')"
            end if
          else if(ick.lt.100) then
            if(filenum.lt.10) then
              fmt1="('dos',i1,'_',i2,'QMC')"
            else if(filenum.lt.100) then
              fmt1="('dos',i2,'_',i2,'QMC')"
            else if(filenum.lt.1000) then
              fmt1="('dos',i3,'_',i2,'QMC')"
            else if(filenum.lt.10000) then
              fmt1="('dos',i4,'_',i2,'QMC')"
            end if
          else
            write(6,*) 'ick too large'
            stop
          end if
            
          write(string1,fmt1) filenum,ick        !open dos_ickQMC
          open(unit=21,file=string1,status='old')
          
          read(21,1000) linemc    
          do i=1,1000
            read(21,*,IOSTAT=istat) w(i),r1,dw(i)
            if(istat.ne.0) goto 998
            Akw(i,ick)=r1
          end do
 998      continue
          if(ick.eq.1) then
            nf=i-1
          else if(nf.ne.i-1) then       ! trap for errors
            write(6,*) 'dos files are not compatible'
            stop
          end if
          
          call spline(w,Akw(1,ick),nf,0.0_8,0.0_8,scAkw(1,ick))

c         To evaluate at w(l)
c         call splint(w,Akw(1,ick),scAkw(1,ick),nf,w(l),A_spline(l))    

          
        end do 
 
c******************************************************************************
c       Now generate an inhomogeneous frequency grid, and calculate the
c       polarization bubble on this grid (both real and imag. parts).
c******************************************************************************

        delta=0.7 
        nfb=300     
        do i=1,nfb
          y=(i-0.5)/float(2*nfb)  !positive freqs. only.
          wb(i)=delta*tan(pi*y)
          if(wb(i).gt.w(nf)) exit
          dwb(i)=pi*(delta**2+wb(i)**2)/(delta*2*nfb)   
        end do 
        nfb=i-1
c                  ---
c             2 pi \   /
c chi2(Q,w) = ---- |   | dw' A(K+Q,w'+w/2) A(K,w'-w/2) [f(w'-w/2)-f(w'+w/2)]
c              Nc  /   /
c                  ---
c                   K
        chi2=0.0
        do ick=1,Ncw                   ! Q
        do jck=1,Nc                    ! K sum over full zone
          isk=ickmap(ickplus(ick,jck)) ! K+Q
          ik=ickmap(jck)               ! K mapped into FBZ
        do i=1,nfb                     ! w
        do j=1,nf                      ! w' sum over frequency
          call splint(w,Akw(1,isk),scAkw(1,isk),nf,w(j)+0.5*wb(i),r1) !A(K+Q,w'+w/2)
          call splint(w,Akw(1,ik),scAkw(1,ik),nf,w(j)-0.5*wb(i),r2)   !A(K,w'-w/2)
          chi2(i,ick)=chi2(i,ick)+dw(j)*r1*r2*
     &         (fermi(w(j)-0.5*wb(i),beta)-fermi(w(j)+0.5*wb(i),beta))
        end do
        end do
        end do
        end do

        chi2=pi*chi2/float(Nc)       ! normalize

c       Check that chi2/w is postive (i.e. like the conductivity)
        do ick=1,Ncw
          if(chi2(1,ick).lt.0.0.or.chi2(nfb/2,ick).lt.0.0) then
            write(6,*) 'sign error in chi2(i,ick),ick=',ick
            stop
          end if
        end do

c       Form the real part, chi1, using Kramers Kronig
c
c                 /    -chi2(x)/pi
c       chi1(w) = | dx -----------
c                 /      w - x
c
        do ick=1,Ncw
        do i=1,nfb
          chi1(i,ick)=0.0
          do j=1,i-1          ! positive x < w
            chi1(i,ick)=chi1(i,ick)+dwb(j)*chi2(j,ick)/(wb(i)-wb(j))
          end do
          do j=i+1,nfb        ! positive x > w
            chi1(i,ick)=chi1(i,ick)+dwb(j)*chi2(j,ick)/(wb(i)-wb(j))
          end do
          do j=1,nfb          ! negative x 
            chi1(i,ick)=chi1(i,ick)-dwb(j)*chi2(j,ick)/(wb(i)+wb(j))
          end do          
        end do
        end do
        chi1=-chi1/pi

c       Check that chi1(w-->0) is postive
        do ick=1,Ncw
          if(chi1(1,ick).lt.0.0) then
            write(6,*) 'sign error in chi1(1,ick),ick=',ick
            stop
          end if
        end do

c******************************************************************************
c       Now chi=chi1+ii*chi2 is the bare bubble.  To form the PT
c       bubble, we will use
c
c                chi                    chi
c       chi_c= --------        chi_s= -------
c               1+U*chi               1-U*chi
c
c******************************************************************************

        chi_c=2.0*(chi1+ii*chi2)/(1+U*(chi1+ii*chi2))
        chi_s=2.0*(chi1+ii*chi2)/(1-U*(chi1+ii*chi2))

c       Form the normalized default model for -cha2(w)/w
        chi2=aimag(chi_c)
        do ick=1,Ncw
        r1=0.0  
        do i=1,nfb
          chi2(i,ick)=chi2(i,ick)/(pi*wb(i))
          r1=r1+dwb(i)*chi2(i,ick)
        end do
        chi2(:,ick)=0.5*chi2(:,ick)/r1
        end do


c******************************************************************************
c       Now that we have the 2-particle spectra, open the appropriate files
c       and write it out. 
c******************************************************************************

        do ick=1,Ncw  ! do all K in the irreducible wedge
c         Open some files.
          if(ick.lt.10) then
            fmt1="('model_chak',i1)"
          else if(ick.lt.100) then
            fmt1="('model_chak',i1)"
          else
            write(6,*) 'ick too large'
            stop
          end if
            
          write(string1,fmt1)        ick        !open modelchakick
          open(unit=21,file=string1,status='unknown')
          do i=1,nfb
            write(21,*) wb(i),chi2(i,ick),dwb(i)
          end do
          close(21)
        end do
        
        open(unit=21,file='model_cha',status='unknown')
        do i=1,nfb
        r1=0.0
        do ick=1,Ncw
          r1=r1+ickdeg(ick)*chi2(i,ick)
        end do
        write(21,*) wb(i),r1/float(Nc),dwb(i)
        end do

c       Form the normalized default model for -chi2(w)/w
        chi2=aimag(chi_s)
        do ick=1,Ncw
        r1=0.0  
        do i=1,nfb
          chi2(i,ick)=chi2(i,ick)/(pi*wb(i))
          r1=r1+dwb(i)*chi2(i,ick)
        end do
        chi2(:,ick)=0.5*chi2(:,ick)/r1
        end do

        do ick=1,Ncw  ! do all K in the irreducible wedge
c         Open some files.
          if(ick.lt.10) then
            fmt1="('model_chik',i1)"
          else if(ick.lt.100) then
            fmt1="('model_chik',i1)"
          else
            write(6,*) 'ick too large'
            stop
          end if
            
          write(string1,fmt1)        ick        !open modelchakick
          open(unit=21,file=string1,status='unknown')
          do i=1,nfb
            write(21,*) wb(i),chi2(i,ick),dwb(i)
          end do
          close(21)
        end do
        
        open(unit=21,file='model_chi',status='unknown')
        do i=1,nfb
        r1=0.0
        do ick=1,Ncw
          r1=r1+ickdeg(ick)*chi2(i,ick)
        end do
        write(21,*) wb(i),r1/float(Nc),dwb(i)
        end do
        
c******************************************************************************
c       Wrap things up 
c******************************************************************************
        call deallocate_1p
        call deallocate_Spec0
        stop
        end program default
                  
        double precision function fermi(x,beta)
c******************************************************************************
c       A fermi function function
c******************************************************************************
        real(8) :: r1,beta,x
        real(8), parameter :: one=1.0_8, half=0.5_8
        r1=exp(-beta*abs(x))
        fermi=half*((r1+one)+(r1-one)*sign(one,x))/(one+r1)
        return
        end
                
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


