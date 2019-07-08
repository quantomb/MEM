        program etawl
c*******************************************************************************
c       to compile type make eta (or make all) using the 
c       associated Makefile.  This code was designed to read the
c       magority of its inputs from a qmc sigma file.
c
c       in this version Sig_sp_i and Derivs_sp_i sit on a much densor homogenous 
c       grid for later eta pick up a inhomogenous grid up to 3 digital
c*******************************************************************************
        use Global
        use module_spectra
        implicit none
c*******************************************************************************
        character*72 linemc,fmt1,string1

        integer, parameter :: nover=100
        integer :: filenum, istat, i, j, k, n, run, meas, ick, 
     &             nkpoints, info, infot, ikx, iky, deg, itunnel,
     &             ios, izero, im, in, nfx, ifx,ifx2,is2,
     &             im1,nfd

        real(kind), parameter :: delta=0.5_kind
        real(kind) :: r1, r2, r3, kx, ky, Epsk, kx1, ky1, kx2, gd, sumgd2, ek,
     &                ky2, offset, dfermi, kplot, fermi,
     &                L11, L12, L21, L22,
     &                t_start, t_stop,
     &                eta, weta(1000), dwm, 
     &                wmax,dw3,y,
     &                r11,r31,eta1,dfermi1
        
        complex(kind) :: c1, c2, c3, blint,
     &                   c21,c31

c********* local allocatable arrays:
        real(kind), allocatable :: ReSigma(:,:), ImSigma(:,:), 
     &                             ReSigmadv(:,:), ImSigmadv(:,:)


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

c    set up frequency range
c    frequency range 2*wmax is set in sopt_dca.f or you could read from Sigmaxxxx_x.dat
        wmax=7.2d0
        dw3=0.004d0   ! It should smaller than the smallest dw in the inhomogenous grid, if not, think dw3i again.
        nfh=int(wmax*2.0/dw3)*2+1
        nf=10000  ! for declaration only.  It will be reset below

        call allocate_Spec0         
        call allocate_Spec1         
       
c*******************************************************************************
c       allocate some local arrays just for etawl code
c*******************************************************************************


        allocate (ReSigma(nf,Nc),stat=info)

        allocate (ImSigma(nf,Nc),stat=info)

        allocate (ReSigmadv(nf,Nc),stat=info)

        allocate (ImSigmadv(nf,Nc),stat=info)


 
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

c 
c       interpolate \Sig(K,\omega) to a 3 digital dense frequency grid
c
c  set up frequency grid, wh: 3 digital homogenous grid.  w3i: 3 digital round off inhomogenous grid
        do i=1,nfh
          wh(i)=-wmax-dw3+real(i)*dw3*0.5
          dwh(i)=dw3*0.5
        enddo
        do n=1,nfh
          write(822,"(i,2(f11.6))")n,wh(n),dwh(n)
        enddo

         write(*,*) 'line220'

        do i=1,nf/2
          n=i+nf/2+1
          y=(i-0.5_kind)/real(nf)
c          w3i(n)=real( dw3* int(delta*tan(pi*y)/dw3-0.5) )+0.5*dw3
          w3i(n)=real( dw3* int(w(n-1)/dw3+0.5) )
c          dw3i(n)=pi*(delta**2+(w3i(n))**2)/(delta*nf)
        enddo
        w3i(nf/2+1)=0.0d0
c        dw3i(n)=pi*(delta**2+0.0d0**2)/(delta*nf)
        do i=1, nf/2
          w3i(i)=-w3i(nf-i+2)
c          dw3i(i)=dw3i(nf-i+2)
        enddo
c  form dw3i by dw3i(n)=(w3i(n+1)-w3i(n-1))/2
        do n=2,nf
          dw3i(n)=(w3i(n+1)-w3i(n-1))/2.0d0
        enddo 
          dw3i(1)=w3i(2)-w3i(1); dw3i(nf+1)=w3i(nf+1)-w3i(nf)

        do n=1,nf+1
          write(811,"(i,4(f11.6))")n,w(n),w3i(n),dw(n),dw3i(n)
        enddo        

        do i=1,nf
        do j=1,Nc
          ReSigma(i,j)=real(Sigma(i,j))
          ImSigma(i,j)=aimag(Sigma(i,j))
        enddo
        enddo
        do i=1,Nc
          call spline(w,ReSigma(:,i),nf,0.0_8,0.0_8,ReSigmadv(:,i))
          call spline(w,ImSigma(:,i),nf,0.0_8,0.0_8,ImSigmadv(:,i))
        enddo

           write(*,*) 'line241, nf=',nf

        if(ispline.eq.1) then

c         Form the 2D interpolation arrays needed to interpolate
c         the self energy with a bicubic spline.  First, form the 
c         grid in the space of the reciprocal lattice vectors g1 and g2
          do i=-N_sp,N_sp
            g1_sp(i)=i
            g2_sp(i)=i
          end do
c         Now form the data to be splined over an extended zone scheme
          do n=1,nfh
            do i=-N_sp,N_sp
            do j=-N_sp,N_sp
              ick=ict(i,j,0)

c              if(i.eq.0 .and. j.eq.0 .)write(*,*) 'line258,n=',n

              call splint(w,ReSigma(:,ick),ReSigmadv(:,ick),nf,wh(n),r1)  !r1=Re(Sig(K,wh))
              call splint(w,ImSigma(:,ick),ImSigmadv(:,ick),nf,wh(n),r2)  !r2=Im(Sig(K,wh))

              Sig_sp_r(i,j,n)=(1-icum)*r1 + icum* real( oner/(w(n)-(r1*(1.0_kind,0.0_kind)+ii*r2)) )
              Sig_sp_i(i,j,n)=(1-icum)*r2 + icum*aimag( oner/(w(n)-(r1*(1.0_kind,0.0_kind)+ii*r2)) )


c              Sig_sp_r(i,j,n)=(1-icum)*real(Sigma(n,ick)) + icum*real(oner/(w(n)-Sigma(n,ick)))
c              Sig_sp_i(i,j,n)=(1-icum)*aimag(Sigma(n,ick)) + icum*aimag(oner/(w(n)-Sigma(n,ick)))
            end do
            end do
c          write(*,*) 'line270,n=',n

c           Form the interpolation data.
            call splie2(g1_sp(-N_sp),g2_sp(-N_sp),
     &                  Sig_sp_r(-N_sp,-N_sp,n),
     &                  2*N_sp+1,2*N_sp+1,
     &                  Derivs_sp_r(-N_sp,-N_sp,n))
            call splie2(g1_sp(-N_sp),g2_sp(-N_sp),
     &                  Sig_sp_i(-N_sp,-N_sp,n),
     &                  2*N_sp+1,2*N_sp+1,
     &                  Derivs_sp_i(-N_sp,-N_sp,n))
          enddo
c          write(*,*) 'line280'
        end if ! (ispline.eq.1)

c check 3 frequency grid w wh w3i,  checked, w3i(i) falls on w(i) exactly
        do i=1,nf
          write(611,*) real(i)/real(nf), w(i)
        enddo
        do i=1,nf+1
          write(622,*) real(i)/real(nf+1), w3i(i)
        enddo
        do i=1,nfh+1
          write(633,*) real(i)/real(nfh+1), wh(i)
        enddo

        do i=1,nf
          write(644,*) real(i)/real(nf), dw(i)
        enddo
        do i=1,nf+1
          write(655,*) real(i)/real(nf+1), dw3i(i)
        enddo



c   interpolate Sigma to denser homogenous grid , checked successfully

        do i=0,2
        do j=0,i
          do n=1,nf
            write(711,*)w(n), ReSigma(n,ict(i,j,0))
            write(712,*)w(n), ImSigma(n,ict(i,j,0))
          enddo
          do n=1,nfh
            write(721,*)wh(n), Sig_sp_r(i,j,n)
            write(722,*)wh(n), Sig_sp_i(i,j,n)
          enddo
          write(711,*)'&'
          write(712,*)'&'
          write(721,*)'&'
          write(722,*)'&'          
        enddo
        enddo      

c  check -w+x/2 and w+x/2 calculation : typical error within 0.005
            k=0
            do ifx=nf/2+20,nf/2+30
            do i=1,nf+1                                        ! index of w
              k=k+1 
 
              im = int( ( w3i(i)+w3i(ifx)*0.5d0)/(dw3*0.5d0) )+(nfh-1)/2+1                            ! calculate index of x/2+w in wh(nfh)
              in = int( (-w3i(i)+w3i(ifx)*0.5d0)/(dw3*0.5d0) )+(nfh-1)/2+1                            ! calculate index of x/2-w in wh(nfh)
c      protect the boundary
c              im = max(im,1); im = min(im,nfh)
c              in = max(in,1); in = min(in,nfh)
              write(733,*) k, ( w3i(i)+w3i(ifx)*0.5d0)-wh(im)
              write(734,*) k, (-w3i(i)+w3i(ifx)*0.5d0)-wh(in)
            enddo
            enddo

        write(6,*) '********************************************'
        call CPU_TIME ( t_start)
        write ( *, * ) 'start timer',nf

cc********************************************
cc form integer table xo2pw(i,ifx) to take care of x/2+w
cc form integer table xo2mw(i,ifx) to take care of x/2-w
cc********************************************
c        do ifx=nf/2-140+1,nf/2+140,3    
c        do i=1,nf
cc          xo2pw(i,ifx)=int(ATAN( (w(ifx)/2+w(i)) /delta)*nf/pi+0.5+nf/2)     
cc          xo2mw(i,ifx)=int(ATAN( (w(ifx)/2-w(i)) /delta)*nf/pi+0.5+nf/2)      
c          dwm=20000.0d0
c          do j=1,nf
c          if( abs(w(ifx)/2.0d0+w(i)-w(j)) .lt. dwm )then
c             dwm=abs(w(ifx)/2.0d0+w(i)-w(j))
c             xo2pw(i,ifx)=j
c          endif
c          enddo
c          xo2pw(i,ifx)=min(xo2pw(i,ifx),nf)
c          xo2pw(i,ifx)=max(xo2pw(i,ifx),1) 
c          dwm=20000.0d0
c          do j=1,nf
c          if( abs(w(ifx)/2.0d0-w(i)-w(j)) .lt. dwm )then
c             dwm=abs(w(ifx)/2.0d0-w(i)-w(j))
c             xo2mw(i,ifx)=j
c          endif
c          enddo
cc          xo2mw(i,ifx)=min(xo2mw(i,ifx),nf) 
cc          xo2mw(i,ifx)=max(xo2mw(i,ifx),1)   
cc          write(*,*) ifx, i, xo2pw(i,ifx), xo2mw(i,ifx)
cc          write(*,*) w(ifx)/2.0d0+w(i), w(xo2pw(i,ifx)), w(ifx)/2.0d0-w(i), w(xo2mw(i,ifx))
c        enddo
c        enddo
cc***************************************************
cc       calculate the normalize factor for the form factor gdwave as: sumgd2=\Sig_k g(k)^2=sum_ick^(nover**2) gdwave(ick)^2
cc***************************************************
c        sumgd2=0.0d0  
c        do ikx=0,nover  ! restrict the calculation to the IRW.
c        do iky=0,ikx
c          if (ikx.eq.iky.and.ikx.eq.0
c     &               .or.
c     &          ikx.eq.iky.and.ikx.eq.nover) then
c              deg=1
c          else if(ikx.eq.iky.or.iky.eq.0.or.ikx.eq.nover) then
c              deg=4          
c          else
c              deg=8
c          end if
c          kx=ikx*pi/float(nover)
c          ky=iky*pi/float(nover)
c          gd = cos(kx) - cos(ky)
c          sumgd2=sumgd2 + gd*gd
c        enddo
c        enddo 
c        write(911,*) '# sumgd2=',sumgd2, ' w, eta(w) i(w)'

        do ifx2 =1,100,3
        do is2 = 1,2
        ifx=(nf/2+1)+ifx2*(3-2*is2)
        eta=0.0
        eta1=0.0
          do ikx=0,nover  ! restrict the calculation to the IRW.
          do iky=0,ikx

c          do ikx=43,43
c          do iky=28,28

            if (ikx.eq.iky.and.ikx.eq.0
     &               .or.
     &          ikx.eq.iky.and.ikx.eq.nover) then
              deg=1
            else if(ikx.eq.iky.or.iky.eq.0.or.ikx.eq.nover) then
              deg=4          
            else
              deg=8
            end if
            kx=ikx*pi/float(nover)
            ky=iky*pi/float(nover)
            gd = cos(kx) - cos(ky)
            ek = Epsk(kx,ky)
            do i=1,nf+1                                        ! index of w
      
              im = int( ( w3i(i)+w3i(ifx)*0.5d0)/(dw3*0.5d0) )+(nfh-1)/2+1                            ! calculate index of x/2+w in wh(nfh)
              in = int( (-w3i(i)+w3i(ifx)*0.5d0)/(dw3*0.5d0) )+(nfh-1)/2+1                            ! calculate index of x/2-w in wh(nfh)
c      protect the boundary
              im = max(im,1); im = min(im,nfh)
              in = max(in,1); in = min(in,nfh)


c              write(*,*)'397 wh(im),in=',im,in             

              dfermi = ( fermi(0.5d0*w3i(ifx)+w3i(i),beta) + fermi(0.5d0*w3i(ifx)-w3i(i),beta) -1 )/ w3i(ifx) 
c              r2     = exp(-beta*abs(w3i(i)))
c              dfermi = -beta*r2/((1.0+r2)**2)                ! df/dw

 
c       [f(x/2+w)+f(x/2-w)-1]/x;  use w3i 3 digital round off inhomogenous w grid.
c              dfermi = ( fermi(0.5*w(ifx)+w(i),beta) + fermi(0.5*w(ifx)-w(i),beta) -1 )/ w(ifx)   ! [f(x/2+w)+f(x/2-w)-1]/x
c              dfermi = ( fermi(w(im),beta)+fermi(w(in),beta)-1 )/ w(ifx)   ! [f(x/2+w)+f(x/2-w)-1]/x 

              c2=0.125d0*(blint(kx,ky,im)  + blint(-kx,-ky,im)   ! Sigma( k,x/2+w)
     &                 +blint(-kx,ky,im) + blint(kx,-ky,im)    ! averaged over
     &                 +blint(ky,kx,im)  + blint(-ky,-kx,im)   ! 8 cubic point
     &                 +blint(-ky,kx,im) + blint(ky,-kx,im))   ! group operations
c              r1     = aimag(1.0d0/(w3i(i)+w3i(ifx)*0.5d0-ed-ek-c2))             ! -pi A(k,x/2+w)
              r1     = aimag(1.0d0/(wh(im)-ed-ek-c2))             ! -pi A(k,x/2+w)

              c3=0.125d0*(blint(kx,ky,in)  + blint(-kx,-ky,in)   ! Sigma(-k,x/2-w)
     &                 +blint(-kx,ky,in) + blint(kx,-ky,in)    ! averaged over
     &                 +blint(ky,kx,in)  + blint(-ky,-kx,in)   ! 8 cubic point
     &                 +blint(-ky,kx,in) + blint(ky,-kx,in))   ! group operations
              r3     = aimag(1.0d0/(wh(in)-ed-ek-c3))             ! -pi A(k,x/2-w)

              eta=eta-dw3i(i)*deg*dfermi*(gd**2)*r1*r3                                    
c eta(x) = SUM_k INT g(k)^2 dw A(k,w+x/2) A(-k,-w+x/2) [f(w+x/2)-f(w-x/2)] /x 

c              write(944,"(f11.6,2x,f11.6,2x,i4)") w3i(i),eta,i
c              write(945,"(f11.6,2x,f11.6,2x,i4)") w3i(i),dfermi,i  !pass
c              write(946,"(f11.6,2x,f11.6,2x,i4)") w3i(i),r1*r3,i
c              write(947,"(f11.6,2x,f11.6,2x,i4)") w3i(i),dw3i(i),i !pass
c              write(948,"(f11.6,2x,f11.6,2x,i4)") w3i(i),real(c2),i
c              write(949,"(f11.6,2x,f11.6,2x,i4)") w3i(i),real(c3),i
c              write(950,"(f11.6,2x,f11.6,2x,i4)") w3i(i),aimag(c2),i
c              write(951,"(f11.6,2x,f11.6,2x,i4)") w3i(i),aimag(c3),i
c              write(952,"(f11.6,2x,f11.6,2x,i4)") w3i(i),r1,i
c              write(953,"(f11.6,2x,f11.6,2x,i4)") w3i(i),r3,i
c              write(954,"(f11.6,2x,f11.6,2x,i4)") ed,ek,i
    
            end do
          end do
          end do

          eta=eta/(float(nover**2)*pi**2)
c          eta1=eta1/(float(nover**2)*pi**2)
          weta(ifx)=eta
          write(448,*) w3i(ifx),eta,ifx

        enddo! use ifx2 and is2 together
          write(922,"(f,1x,f,2x,i4)") w3i(ifx),
     &         weta(nf/2+1+ifx2)/2.0d0+weta(nf/2+1-ifx2)/2.0d0,ifx
        enddo! to determine ifx

        do ifx=nf/2+2,nf/2+2+99,3
          write(444,"(1x,f9.4,2x,f9.4,2x,i4)") w3i(ifx),weta(ifx)/2.0d0+weta(nf-ifx+2)/2.0d0,ifx
          write(911,"(1x,f9.4,2x,f9.4,2x,i4)") w3i(ifx),weta(ifx)/2.0d0+weta(nf-ifx+2)/2.0d0,ifx
        enddo

        call CPU_TIME ( t_stop )
        write ( *, * ) 'Elapsed CPU time = ', t_stop - t_start
c*******************************************************************************        
c       wrap up
c*******************************************************************************        
        call deallocate_1p
        call deallocate_Spec0
        call deallocate_Spec1
        deallocate (ReSigma  )
        deallocate (ImSigma  )
        deallocate (ReSigmadv)
        deallocate (IMSigmadv)
      
        stop
        end program etawl
       
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
     

          if(r8.gt.-gamma) r8=-gamma
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
      PARAMETER (NN=10000) 
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
      PARAMETER (NN=10000) 
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
c  Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., y_i=f(x_i), with 
c  x_1<x_2<...<x_N, and given values yp1 and ypn for the first derivative of the interpolating 
c  function at points 1 and n, respectively, this routine returns an array y2(1:n) of length n 
c  which contains the second derivatives of the interpolating function at the tabulated points
c  x_i. 
      IMPLICIT real(8) (A - H, O - Z) 
      PARAMETER (NMAX=10000) 
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
c  Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the xa_i's 
c  in order), and given the array y2a(1:n), which is the output from spline above, and given a 
c  value of x, this routine returns a cubic-spline interpolated value y.
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
