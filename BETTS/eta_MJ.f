        program eta_MJ
c*******************************************************************************
c       to compile type make spectra (or make all) using the 
c       associated Makefile.  This code was designed to read the
c       magority of its inputs from a qmc sigma file.
c*******************************************************************************
        use Global
        use module_spectra
        implicit none
c*******************************************************************************
        character*72 linemc,fmt1,string1

        integer, parameter :: nover=50, nfp=4000
        integer :: filenum, istat, i, j, n, run, meas, ick, ibofip(0:nfp),  iaofip(0:nfp), 
     &             nkpoints, info, infot, ikx, iky, deg((2*nover+3)**2/8), itunnel,
     &             ios, izero, im, in, nfx, ifx, iplow, iphigh, ip, ikdex, ipx, ib, ia

        real(kind), parameter :: delta=0.8_kind
        real(kind) :: r1, r2, r3, kx, ky, Epsk, kx1, ky1, kx2, sumgd2, ek,
     &                ky2, offset, dfermi,kplot, fermi,
     &                L11, L12, L21, L22,
     &                t_start, t_stop,
     &                eta(0:nfp), dwm ,Akwp((2*nover+3)**2/8,0:nfp), gd((2*nover+3)**2/8),
     &                fermip(0:nfp), wp(0:nfp), dwp, sgd, eta0
        
        complex(kind) :: c1, c2, c3, blint
c*******************************************************************************
c       Open some files and read in some parameters
c*******************************************************************************
c       Open a datafile for the parameters (if the filenumber
c       is negative, then the tunneling data will be calculated
c       and the sign changed to positive).
        write(6,"('If filenumber<0, Gtunnel is calculated')")
        write(6,"('       Enter the QMC sigma filenumber:  ',$)") 
        read(5,*) filenum

        if(filenum.lt.0) then
          filenum=-filenum
          itunnel=1
        else
          itunnel=0
        end if

c       Now open the QMC self energy sigma#.dat
        if(filenum.lt.10) then  ! P. 402, Chapman
          fmt1="('sigma',i1,'.dat')"
        else if(filenum.lt.100) then
          fmt1="('sigma',i2,'.dat')"
        else if(filenum.lt.1000) then
          fmt1="('sigma',i3,'.dat')"
        else if(filenum.lt.10000) then
          fmt1="('sigma',i4,'.dat')"
        else
          write(6,*) 'filenum too large'
          stop
        end if  
        write(string1,fmt1) filenum
        open(unit=20,file=string1,status='old')
             
c       Read some parameters from the QMC sigma#.dat file
 401    read(20,"(a72)") linemc
        if(index(linemc(1:72),'Results').ne.0 ) then
          read(20,*) nl,nwn,fill,cluster,ntot,ndim
        else
          goto 401
        endif
        if(ndim.ne.2) stop 'bad dimensionality in spectra'

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

c         Open a file to store w(ifx),weta(ifx)+weta(nf-ifx+1),ifx
          open(unit=444,file='weta2.dat',status='unknown')
c         Open a file to store w(ifx),weta(ifx)
          open(unit=448,file='weta.dat',status='unknown')
c*******************************************************************************
c       allocate some arrays and get some symmetry tables
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
          
          read(21,"(a72)") linemc    
          do i=1,1000
            read(21,*,IOSTAT=istat) w(i),r1,r2,r3,dw(i)
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

c       Now fill in Sigma(w,k), for k in the first Brillouin zone but
c       outside the irreducible wedge.
        do ick=Ncw+1,Nc
        do i=1,nf
          Sigma(i,ick)=Sigma(i,ickmap(ick))
        end do
        end do

        if(ispline.eq.1) then
c         Form the 2D interpolation arrays needed to interpolate
c         the self energy with a bicubic spline.  First, form the 
c         grid in the space of the reciprocal lattice vectors g1 and g2
          do i=-N_sp,N_sp
            g1_sp(i)=i
            g2_sp(i)=i
          end do

c         Now form the data to be splined over an extended zone scheme
          do n=1,nf
            do i=-N_sp,N_sp
            do j=-N_sp,N_sp
              ick=ict(i,j,0)
              Sig_sp_r(i,j,n)=(1-icum)*real(Sigma(n,ick)) + icum*real(oner/(w(n)-Sigma(n,ick)))
              Sig_sp_i(i,j,n)=(1-icum)*aimag(Sigma(n,ick)) + icum*aimag(oner/(w(n)-Sigma(n,ick)))
            end do
            end do
c           Form the interpolation data.
            call splie2(g1_sp(-N_sp),g2_sp(-N_sp),
     &                  Sig_sp_r(-N_sp,-N_sp,n),
     &                  2*N_sp+1,2*N_sp+1,
     &                  Derivs_sp_r(-N_sp,-N_sp,n))
            call splie2(g1_sp(-N_sp),g2_sp(-N_sp),
     &                  Sig_sp_i(-N_sp,-N_sp,n),
     &                  2*N_sp+1,2*N_sp+1,
     &                  Derivs_sp_i(-N_sp,-N_sp,n))
          end do
        end if ! (ispline.eq.1)
      
        write(6,*) '********************************************'
        call CPU_TIME ( t_start)
        write ( *, * ) 'start timer',nf


c       Now we must introduce a homogeneous grid wp, and interpolate A(k,w) onto this grid
        dwp=(w(nf)-w(1))/real(nfp)   
        do ip = 0 , nfp
          wp(ip) = w(1) + ip*dwp  ! Then for a general frequency x, the closest ip = int((x-w(1)+eps)/dwp)
          fermip(ip) = fermi(wp(ip),beta)
        end do

c       Form ibofip a table which give the w(i) closest to but below wp(io)
        do ip=1,nfp-1
        do i=1,nf
          if(w(i).ge.wp(ip)) then
            ibofip(ip)=i-1
            exit
          end if
        end do
        end do
c       Form iaofip a table which give the w(i) closest to but above wp(io)
        do ip=1,nfp-1
        do i=nf,1,-1
          if(w(i).le.wp(ip)) then
            iaofip(ip)=i+1 
            exit
          end if
        end do
        end do
c       Now interpolate A(k,w) onto the homogeneous wedge
        ikdex = 0; sgd=0.0d0; eta0=0.0
        do ikx=0,nover  ! restrict the calculation to the IRW.
        do iky=0,ikx
          ikdex = ikdex +1 
          if (ikx.eq.iky.and.ikx.eq.0
     &             .or.
     &        ikx.eq.iky.and.ikx.eq.nover) then
            deg(ikdex)=1
          else if(ikx.eq.iky.or.iky.eq.0.or.ikx.eq.nover) then
            deg(ikdex)=4          
          else
            deg(ikdex)=8
          end if
          kx=ikx*pi/float(nover)
          ky=iky*pi/float(nover)
          ek = Epsk(kx,ky)
          gd(ikdex) = cos(kx) - cos(ky)

          sgd=sgd+deg(ikdex)*(gd(ikdex))**2

          c2=0.125*(blint(kx,ky,1)  + blint(-kx,-ky,1)     ! Sigma( k,w(1))
     &               +blint(-kx,ky,1) + blint(kx,-ky,1)    ! averaged over
     &               +blint(ky,kx,1)  + blint(-ky,-kx,1)   ! 8 cubic point
     &               +blint(-ky,kx,1) + blint(ky,-kx,1))   ! group operations
          Akwp(ikdex,0)= aimag(1.0/(w(1)-ed-ek-c2))        ! -pi A(k,w(1))
          c2=0.125*(blint(kx,ky,nf)  + blint(-kx,-ky,nf)     ! Sigma( k,w(1))
     &               +blint(-kx,ky,nf) + blint(kx,-ky,nf)    ! averaged over
     &               +blint(ky,kx,nf)  + blint(-ky,-kx,nf)   ! 8 cubic point
     &               +blint(-ky,kx,nf) + blint(ky,-kx,nf))   ! group operations
          Akwp(ikdex,nfp)= aimag(1.0/(w(nf)-ed-ek-c2))        ! -pi A(k,w(1))
          
          do ip=1,nfp-1
            ib=ibofip(ip)                                      ! index of w <~ wp(ip)
            c2=0.125*(blint(kx,ky,ib)  + blint(-kx,-ky,ib)     ! Sigma( k,w(1))
     &               +blint(-kx,ky,ib) + blint(kx,-ky,ib)    ! averaged over
     &               +blint(ky,kx,ib)  + blint(-ky,-kx,ib)   ! 8 cubic point
     &               +blint(-ky,kx,ib) + blint(ky,-kx,ib))   ! group operations
            r1= aimag(1.0/(w(ib)-ed-ek-c2))                  ! -pi A(k,w(i))
            ia=iaofip(ip)                                      ! index of w >~ wp(ip)
            c3=0.125*(blint(kx,ky,ia)  + blint(-kx,-ky,ia)   ! Sigma( k,x/2+w)
     &               +blint(-kx,ky,ia) + blint(kx,-ky,ia)    ! averaged over
     &               +blint(ky,kx,ia)  + blint(-ky,-kx,ia)   ! 8 cubic point
     &               +blint(-ky,kx,ia) + blint(ky,-kx,ia))   ! group operations
            r2     = aimag(1.0/(w(ia)-ed-ek-c3))             ! -pi A(k,w(i+1))
            Akwp(ikdex,ip) = r1 + (r2-r1)*(wp(ip)-w(ib))/(w(ia)-w(ib))
          end do
          do i=1,nf                                        ! index of w
            im = nf - i + 1                                ! index of -w
            r2     = exp(-beta*abs(w(i)))
            dfermi = -beta*r2/((1.0+r2)**2)                ! df/dw
            c2=0.125*(blint(kx,ky,i)  + blint(-kx,-ky,i)   ! Sigma(k,w)
     &               +blint(-kx,ky,i) + blint(kx,-ky,i)    ! averaged over
     &               +blint(ky,kx,i)  + blint(-ky,-kx,i)   ! 8 cubic point
     &               +blint(-ky,kx,i) + blint(ky,-kx,i))   ! group operations
            r1     = aimag(1.0/(w(i)-ed-ek-c2))            ! A(k,w)
            c3=0.125*(blint(kx,ky,im)  + blint(-kx,-ky,im)   ! Sigma(k,w)
     &               +blint(-kx,ky,im) + blint(kx,-ky,im)    ! averaged over
     &               +blint(ky,kx,im)  + blint(-ky,-kx,im)   ! 8 cubic point
     &               +blint(-ky,kx,im) + blint(ky,-kx,im))   ! group operations
            r3     = aimag(1.0/(w(im)-ed-ek-c3))            ! A(k,-w)
            eta0=eta0 -                                    ! eta = SUM_k INT g(k)^2 dw A(k,w) A(-k,-w) -df/dw
     &         dw(i)*deg(ikdex)*dfermi*(gd(ikdex)**2)*r1*r3
          end do
        end do
        end do
        eta0=2.0*eta0/(pi*sgd)
 
c       Test the interpolant at k=0
        do ip=0,nfp
          write(455,*) wp(ip),Akwp(1,ip)
        end do
        write(455,*) '& '
        do i=1,nf
          write(455,*) w(i),aimag(1.0/(w(i)-ed- Epsk(0,0)-blint(0,0,i)))
        end do
        write(455,*) '& '
c       Now use Akw(ikdex,ip), wp(ip), fermip(ip), to calculate
c
c       chi''(w)/w = SUM_k INT g(k)^2 dw A(k,w+x/2) A(-k,-w+x/2) [f(w+x/2)-f(w-x/2)] /x   
c       or        
c       chi''(w)/w = SUM_k INT dw g(k)^2  A(k,w+x) A(-k,-w) ( f(w+x) - f(w)) /x
c
c       ip = int((x-w(1)+eps)/dwp)
c
c       ip=nfp/2 <=> wp(ip)=0         -wp(ip) <=> wp(nfp-ip)
c
c       Then, since w(1) = -dwp*nfp/2
c       wp(ip) + wp(ipp) = 2*w(1) + (ip+ipp)*dwp = w(1) + (ip + ipp -nfp/2)*dwp = wp(ip+ipp-nfp/2)
c       wp(ip) - wp(ipp) = (ip-ipp)*dwp = -dwp*nfp/2 + (ip-ipp+nfp/2)*dwp = w(1) + (ip-ipp+nfp/2)*dwp = wp(ip-ipp+nfp/2)
        
        do ipx= int((-2.0-w(1)+eps)/dwp),int((2.0-w(1)+eps)/dwp)
          eta(ipx)=0.0
          ikdex=0
          do ikx=0,nover  ! restrict the sum on k to the IRW.
          do iky=0,ikx
            ikdex = ikdex + 1 
            do ip = max(ipx-nfp/2,-ipx+nfp/2),min(nfp+ipx-nfp/2,nfp-ipx+nfp/2)      ! max(w(1) +/-x) < w < min(w(nf) +/-x)                              
              eta(ipx) = eta(ipx) + 
     &        deg(ikdex)*(gd(ikdex)**2)*Akwp(ikdex,ip+ipx-nfp/2)*Akwp(ikdex,nfp-ip)*(fermip(ip+ipx-nfp/2)-fermip(ip))
            end do
          end do
          end do

          eta(ipx)=eta(ipx)*dwp/(pi*sgd)
          write(448,"(1x,f9.4,2x,f9.4,2x,i4)") wp(ipx),eta(ipx),ipx
        enddo ! ipx

        write(444,"(1x,f9.4,2x,f9.4,2x,f9.4)") wp(nfp/2),eta0,dwp
        do ipx= nfp/2+1,int((2.0-w(1)+eps)/dwp)
           write(444,"(1x,f9.4,2x,f9.4,2x,f9.4)") wp(ipx),-(eta(ipx)-eta(nfp-ipx))/wp(ipx),dwp
        enddo ! ipx

        call CPU_TIME ( t_stop )
        write ( *, * ) 'Elapsed CPU time = ', t_stop - t_start
c*******************************************************************************        
c       wrap up
c*******************************************************************************        
        call deallocate_1p
        call deallocate_Spec0
        call deallocate_Spec1      
        stop
        end program eta_MJ
       
        complex(8) function blint(kx,ky,n)
c       
c       This subroutine performs a bilinear interpolation of
c       the self energy, at any location (kx,ky) in the zone.
c
        use Global
        use module_spectra
        implicit none
        integer :: n,i1,j1,i2,j2,ic1,ic2,ic3,ic4,ic,ik
        real(kind) :: Kx,Ky,dx,dy,g1,g2,r7,r8,Kx3,Ky3
        complex(kind) :: z1,z2,z3,z4,Sapprox(Nc)
        
        if(ispline.eq.0) then ! perform a bilinear interpolation
          if(Nc.ne.4) stop'doesnt work for betts'
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
c         First find the coordinates of the cell which contains (kx,ky)

          i1=int((ky*g(2,1)-kx*g(2,2))/(g(2,1)*g(1,2)-g(1,1)*g(2,2))+Nc)
     &      -Nc
          j1=int((kx*g(1,2)-ky*g(1,1))/(g(2,1)*g(1,2)-g(1,1)*g(2,2))+Nc)
     &      -Nc

          i2=i1+1
          j2=j1+1
c          index the four corners
          ic1=ict(i1,j2,0)
          ic2=ict(i2,j2,0)
          ic3=ict(i1,j1,0)
          ic4=ict(i2,j1,0)
c         now get the K-point on the lower left.
          Kx3=i1*g(1,1)+j1*g(2,1)
          Ky3=i1*g(1,2)+j1*g(2,2)

          z1=(1-icum)*Sigma(n,ic1) + icum/(w(n)-Sigma(n,ic1))
          z2=(1-icum)*Sigma(n,ic2) + icum/(w(n)-Sigma(n,ic2))
          z3=(1-icum)*Sigma(n,ic3) + icum/(w(n)-Sigma(n,ic3))
          z4=(1-icum)*Sigma(n,ic4) + icum/(w(n)-Sigma(n,ic4))
 
c         calculate the slopes
          dx=((kx-Kx3)*g(1,1)+(ky-Ky3)*g(1,2))
     &      /(g(1,1)*g(1,1)+g(1,2)*g(1,2))
          dy=((kx-Kx3)*g(2,1)+(ky-Ky3)*g(2,2))
     &      /(g(2,1)*g(2,1)+g(2,2)*g(2,2))
c         interpolate.
          blint=z3+(z4-z3)*dx+(z1-z3)*dy+(z2-z1-z4+z3)*dx*dy
          blint=(1-icum)*blint + icum*(w(n)-oner/blint)          ! use Sigma or cumulant, depending on icum

        else if(ispline.eq.1) then ! perform a bicubic interpolation
c         Use a bicubic spline.
c         Find the location on the grid of principle translation
c         vectors in K-space.
          g1=(kx*g(2,2)-ky*g(2,1))/(g(1,1)*g(2,2)-g(1,2)*g(2,1))
          g2=(kx*g(1,2)-ky*g(1,1))/(g(2,1)*g(1,2)-g(2,2)*g(1,1))
c         Interpolate for the real part of Sigma(kx,ky,n) or the cumulant
          call splin2(g1_sp(-N_sp),g2_sp(-N_sp),
     &                Sig_sp_r(-N_sp,-N_sp,n),
     &                Derivs_sp_r(-N_sp,-N_sp,n),
     &                2*N_sp+1,2*N_sp+1,
     &                g1,g2,r7)
c         Interpolate for the imaginary part of Sigma(kx,ky,n) or the cumulant
          call splin2(g1_sp(-N_sp),g2_sp(-N_sp),
     &                Sig_sp_i(-N_sp,-N_sp,n),
     &                Derivs_sp_i(-N_sp,-N_sp,n),
     &                2*N_sp+1,2*N_sp+1,
     &                g1,g2,r8)
     

          blint=cmplx(r7,r8)  !return sigma or the cumulant
          blint=(1-icum)*blint + icum*(w(n)-oner/blint)          ! use Sigma or cumulant, depending on icum
       else if(ispline.eq.2) then ! use Fourier interpolation
c interpolation by introducing a factor exp(i kt X) in sigma
          Sapprox=(0.d0,0.d0)
          do ic=1,Nc
          do ik=1,Nc        
             Sapprox(ic)= Sapprox(ic)
     &             +FTCoefs_K_to_R(ik,ic)*((1-icum)*Sigma(n,ik)+icum/(w(n)-Sigma(n,ik)))
          enddo
          enddo

          blint=(0.d0,0.d0)
          do ic=1,Nc                              
             z1=exp(+ii*(kx*Rc(1,ic)+ky*Rc(2,ic)))
             blint=blint+z1*Sapprox(ic)
          enddo
          blint=(1-icum)*blint + icum*(w(n)-oner/blint)          ! use Sigma or cumulant, depending on icum
        else
          write(6,*) 'BAD ispline'
          stop
        end if
        return
        end
                
        
        real*8 function Epsk(kx,ky)
c******************************************************************************
c       2D dispersion function
c******************************************************************************
        use Global
        implicit none
        real(kind) :: kx, ky    
       
         Epsk=-0.5_kind*(cos(kx)+cos(ky)) -
     &       tprime*(cos(kx)*cos(ky)-1.0_kind)
     &             - 0.5*tdprime*(cos(2.d0*kx)+cos(2.d0*ky))
     &             - t3prime*(cos(2.d0*kx)*cos(ky)
     &                         +cos(kx)*cos(2.d0*ky))       
        end function Epsk
        

        double precision function fermi(x,beta)
c******************************************************************************
c       A fermi function
c******************************************************************************
        real(8) :: r1,beta,x
        real(8), parameter :: one=1.0_8, half=0.5_8
        r1=exp(-beta*abs(x))
        fermi=half*((r1+one)+(r1-one)*sign(one,x))/(one+r1)
        return
        end
        
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
