        program spectra
c*******************************************************************************
c       to compile type make spectra (or make all) using the 
c       associated Makefile.  This code was designed to read the
c       magority of its inputs from a qmc sigma file.
c*******************************************************************************
c*******************************************************************************
        use Global
        use module_spectra
        implicit none
c*******************************************************************************
        character*72 linemc,fmt1,string1

        integer, parameter :: nover=100
        integer :: filenum, istat, i, j, n, run, meas, ick, 
     &             nkpoints, info, infot, ikx, iky, deg, itunnel,
     &             ios, izero, im

        real(kind) :: r1, r2, r3, kx, ky, Epsk, kx1, ky1, kx2, gd, ek,
     &                ky2, offset, res, hall, eta, vx, vyy, dfermi,kplot,
     &                L11, L12, L21, L22, thermopower, peltier, kappa,
     &                sgd, kfx, kfy
        
        complex(kind) :: c1, c2, c3, blint
        
c*******************************************************************************


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
        
        if (iconductivity.eq.1) then
        
c         Open a file for the lifetime tau(w) 
          if(filenum.lt.10) then  ! P. 402, Chapman
            fmt1="('tau',i1,'.dat')"
          else if(filenum.lt.100) then
            fmt1="('tau',i2,'.dat')"
          else if(filenum.lt.1000) then
            fmt1="('tau',i3,'.dat')"
          else if(filenum.lt.10000) then
            fmt1="('tau',i4,'.dat')"
          else
            write(6,*) 'filenum to large'
            stop
          end if  
          write(string1,fmt1) filenum
          open(unit=45,file=string1,status='unknown')
          
c         Open a file for the resistivity, L_ij, etc.
          open(unit=44,file='resis.dat',status='unknown')
        end if



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
          do i=1,2000
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
        do i=1,nf
          kx=pi/2
          ky=pi/2
          c1=0.125*(blint(kx,ky,i)  + blint(-kx,-ky,i)  ! Sigma(k,w) with C_4
     &                +blint(-kx,ky,i) + blint(kx,-ky,i)   ! symmetry imposed
     &                +blint(ky,kx,i)  + blint(-ky,-kx,i)
     &                +blint(-ky,kx,i) + blint(ky,-kx,i)  )
          write(666,*) w(i),real(c1),aimag(c1)
        end do 


c****************************************************************************************************
c        find the vectors k1 and k2 
c****************************************************************************************************

        write(6,"('Enter kx1,ky1,kx2,ky2,offset,nkpoints:  ',$)") 
c 3      format('Enter kx1,ky1,kx2,ky2,offset,nkpoints:  ',$)
        read(5,*) kx1,ky1,kx2,ky2,offset,nkpoints
             
        write(6,*) '********************************************'
        write(6,*) '*** G(k,w) and Sigma(k,w) are in Spectra ***'
        write(6,*) '*** A(k,w) is in Akw_spectra             ***'
        write(6,*) '*** Sigma(k_F,w) is in Sigmak_F          ***'

c****************************************************************************************************
c       Find the Fermi surface along the direction k1 --> k2 using a very fine k-mesh
c****************************************************************************************************
        r2=zeror
        do ick=1,10*nkpoints 
          kx=kx1 + real(ick-1)*(kx2-kx1)/max(10*nkpoints-1,1)
          ky=ky1 + real(ick-1)*(ky2-ky1)/max(10*nkpoints-1,1)
c         w=0 is half way between nf/2 and nf/2+1
          i=nf/2
          c1=0.125*(blint(kx,ky,i)  + blint(-kx,-ky,i)  ! Sigma(k,w) with C_4
     &                +blint(-kx,ky,i) + blint(kx,-ky,i)   ! symmetry imposed
     &                +blint(ky,kx,i)  + blint(-ky,-kx,i)
     &                +blint(-ky,kx,i) + blint(ky,-kx,i)  )
          i=nf/2+1
          c2=0.125*(blint(kx,ky,i)  + blint(-kx,-ky,i)  ! Sigma(k,w) with C_4
     &                +blint(-kx,ky,i) + blint(kx,-ky,i)   ! symmetry imposed
     &                +blint(ky,kx,i)  + blint(-ky,-kx,i)
     &                +blint(-ky,kx,i) + blint(ky,-kx,i)  )
          r1=-aimag(1.0/(w(nf/2)-ed-Epsk(kx,ky)-c1) + 1.0/(w(nf/2+1)-ed-Epsk(kx,ky)-c2))
          if(r1.gt.r2) then
            kfx=kx
            kfy=ky
            r2=r1
          end if
        end do
        write(6,*) 'k_F=',kfx,kfy
c****************************************************************************************************
c       Calculate and export A(k_f,w) and Sigma(k_F,w)
c****************************************************************************************************
        open(unit=27,file='Sigmak_F.dat',status='unknown')
        write(27,"('#   w        Re(Sig(k_F,w)  Im(Sig(k_F,w)  A(k_F,w)')")
        do i=1,nf
          c1=0.125*(blint(kfx,kfy,i)  + blint(-kfx,-kfy,i)  ! Sigma(k,w) with C_4
     &                +blint(-kfx,kfy,i) + blint(kfx,-kfy,i)   ! symmetry imposed
     &                +blint(kfy,kfx,i)  + blint(-kfy,-kfx,i)
     &                +blint(-kfy,kfx,i) + blint(kfy,-kfx,i)  )
          write(27,"(4(g12.6,2x))") w(i),real(c1),aimag(c1),-aimag(1.0/(w(i)-ed-Epsk(kfx,kfy)-c1))/pi           
        end do

          
c****************************************************************************************************
c       Calculate and export A(k,w) along the prescribed direction
c****************************************************************************************************
        open(unit=22,file='Spectra',status='unknown')
        write(22,"('#   w          ReGkw        ImGkw       ReSigkw  ',
     &             '    ImSigkw')")
        open(unit=21,file='Akw_spectra',status='unknown')
        write(21,"('#   w      A(k,w)=-Im G(k,w)/pi ')")

        open(unit=24,file='False_color.dat',status='unknown')
        write(24,"('# k     w      A(k,w)=-Im G(k,w)/pi ')")


        do ick=1,nkpoints ! nkpoints may be 1
          kx=kx1 + real(ick-1)*(kx2-kx1)/max(nkpoints-1,1)
          ky=ky1 + real(ick-1)*(ky2-ky1)/max(nkpoints-1,1)

          do i=1,nf
             c2=0.125*(blint(kx,ky,i)  + blint(-kx,-ky,i)  ! Sigma(k,w) with C_4
     &                +blint(-kx,ky,i) + blint(kx,-ky,i)   ! symmetry imposed
     &                +blint(ky,kx,i)  + blint(-ky,-kx,i)
     &                +blint(-ky,kx,i) + blint(ky,-kx,i)  )
c            c2=0.5*(blint(kx,ky,i)+blint(ky,kx,i))
c            c2=blint(kx,ky,i)
            c1=1.0/(w(i)-ed-Epsk(kx,ky)-c2)            !G(k,w)
            write(22,"(5(g12.6,1x))") w(i),real(c1),aimag(c1),
     &                           real(c2),aimag(c2)
            write(21,"(5(g12.6,1x))") w(i),
     &                            -aimag(c1)/pi+(ick-1)*offset
        
            write(24,"(5(g12.6,1x))") ick,w(i),
     &                            -aimag(c1)/pi

          end do
          write(21,*) '& '
          write(22,*) '& '
          write(24,*) ' '
        end do

 
c****************************************************************************************************
c     Calculate and export the MDC curves
c****************************************************************************************************

      open(unit=25,file='MDC',status='unknown')
      write(6,*) '*** MDC distribution of A(k,w) is in MDC***'
      write(25,"('#   kplot      A(k,w)          w          kx ',
     &      '        ky')")
 
        do i=1,nf,5
        if ((w(i).le.0.d0).and.(w(i).gt.-4.d0)) then
        do ick=1,nkpoints
          kx=kx1 + real(ick-1)*(kx2-kx1)/max(nkpoints-1,1)
          ky=ky1 + real(ick-1)*(ky2-ky1)/max(nkpoints-1,1)

          kplot=ky
         if (abs(kx2-kx1).gt.eps) kplot=kx
        

           c2=0.125*(blint(kx,ky,nf/2)  + blint(-kx,-ky,nf/2)
     &            +blint(-kx,ky,nf/2) + blint(kx,-ky,nf/2)
     &            +blint(ky,kx,nf/2)  + blint(-ky,-kx,nf/2)
     &            +blint(-ky,kx,nf/2) + blint(ky,-kx,nf/2)   )
c            c2=0.5*(blint(kx,ky,i)+blint(ky,kx,i))
c            c2=blint(kx,ky,i)
            c1=1.0/(w(nf/2)-ed-Epsk(kx,ky)-c2)


c       if ((w(i).le.0.d0).and.(w(i).gt.-4.d0)) then
        write(25,"(5(g12.6,1x))") kplot, -aimag(c1)/pi+(i-1-nf/4)*offset*0.25, 
     &   w(i), kx, ky
c       endif   
 
        enddo   !ick
        write(25,*) '& '
        endif
        enddo   !i

c****************************************************************************************************
c       Now interpolate the self energy to calculate A(k,w=0) in the
c       upper right quadrant of the first Brillouin zone.
c****************************************************************************************************
        open(unit=26,file='Akw0',status='unknown')
        write(6,*) '***  distribution of A(k,w=0) is in Akw0***'
        write(25,"('#   kx     ky      A(k,w=0)     ')")
c       Find the zero frequency index
        do i=1,nf
          if(w(i).ge.0.0d0) then
            izero=i
            exit
          end if
        end do
 
        do ikx=0,60
        do iky=0,60
          kx=ikx*pi/real(60)
          ky=iky*pi/real(60)
          c2=0.125*(blint(kx,ky,izero)  + blint(-kx,-ky,izero)
     &            +blint(-kx,ky,izero) + blint(kx,-ky,izero)
     &            +blint(ky,kx,izero)  + blint(-ky,-kx,izero)
     &            +blint(-ky,kx,izero) + blint(ky,-kx,izero)   )
          c1=1.0/(w(izero)-ed-Epsk(kx,ky)-c2)
          r1=(-1.0d0/pi)*aimag(c1)
          write(26,*) kx,ky,r1
 
        enddo   !ikx
        write(26,*) '  '
        enddo   !iky

      
        if(itunnel.eq.1) then        
c****************************************************************************************************
c         Now calculate the Green function near the FS, to be used to 
c         calculate the STM data.
c****************************************************************************************************

          write(6,*) '*** G(k,w) is tabulated in Gtunnel       ***'
          open(unit=23,file='Gtunnel',status='unknown')
          write(23,"('#     w          kx          ky         ReGkw ',
     &               '       ImGkw ')")

          do ikx=0,256
          do iky=0,ikx
            kx=ikx*pi/real(256)
            ky=iky*pi/real(256)
            do i=1,nf
              if(abs(w(i)).lt.0.2) then
                 c2=0.125*(blint(kx,ky,i)  + blint(-kx,-ky,i)
     &                    +blint(-kx,ky,i) + blint(kx,-ky,i)
     &                    +blint(ky,kx,i)  + blint(-ky,-kx,i)
     &                    +blint(-ky,kx,i) + blint(ky,-kx,i)  )
c                 c2=0.5*(blint(kx,ky,i)+blint(ky,kx,i))
c                c2=blint(kx,ky,i)
                c1=1.0/(w(i)-ed-Epsk(kx,ky)-c2)
                write(23,"(5(g12.6,1x))") w(i),kx,ky,real(c1),aimag(c1)
              end if
            end do
          end do
          end do  
        end if                
          
 
      
        write(6,*) '********************************************'
        if (iconductivity.eq.1) then
c****************************************************************************************************
c         Now that we have the interpolated spectra, calculate the 
c         transport using the lifetime approximation (which is only
c         exact when D=oo, or Nc=1 or 4).
c****************************************************************************************************

          res=0.0; hall=0.0; eta=0.0; tau(:)=0.0; sgd = 0.0
          do ikx=0,nover  ! restrict the calculation to the IRW.
          do iky=0,ikx
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
            vx = sin(kx)*(0.5 + tprime*cos(ky))
            vyy = cos(ky)*(0.5 + tprime*cos(kx))
            gd = cos(kx) - cos(ky)
            sgd=sgd+deg*(gd)**2
            ek = Epsk(kx,ky)
            do i=1,nf                                        ! index of w
              im = nf - i + 1                                ! index of -w
              r2     = exp(-beta*abs(w(i)))
              dfermi = -beta*r2/((1.0+r2)**2)                ! df/dw
              c2=0.125*(blint(kx,ky,i)  + blint(-kx,-ky,i)   ! Sigma(k,w)
     &                 +blint(-kx,ky,i) + blint(kx,-ky,i)    ! averaged over
     &                 +blint(ky,kx,i)  + blint(-ky,-kx,i)   ! 8 cubic point
     &                 +blint(-ky,kx,i) + blint(ky,-kx,i))   ! group operations
              r1     = aimag(1.0/(w(i)-ed-ek-c2))            ! -pi A(k,w)
              c3=0.125*(blint(kx,ky,im)  + blint(-kx,-ky,im) ! Sigma(-k,-w)
     &                 +blint(-kx,ky,im) + blint(kx,-ky,im)  ! averaged over
     &                 +blint(ky,kx,im)  + blint(-ky,-kx,im) ! 8 cubic point
     &                 +blint(-ky,kx,im) + blint(ky,-kx,im)) ! group operations
              r3     = aimag(1.0/(w(im)-ed-ek-c3))           ! -pi A(+/-k,-w)
              tau(i)=tau(i) + deg*(r1*vx)**2                 ! tau(w)
              res=res -
     &           dw(i)*deg*dfermi*(r1*vx)**2
              hall=hall +
     &           dw(i)*deg*dfermi*(vx**2)*vyy*(r1)**3
              eta=eta -                                      ! eta = SUM_k INT g(k)^2 dw A(k,w) A(-k,-w) -df/dw
     &           dw(i)*deg*dfermi*(gd**2)*r1*r3
            end do
          end do
          end do
          res=(float(nover**2)*pi**2)/res         ! 1/L11 units e^2 pi/(2 hbar a)
          hall=hall/(float(nover**2)*pi**3)       ! 2 pi^2 e^3 a B /(3 hbar^2)
          hall=hall/res**2
          tau=tau/(float(nover**2)*pi**2)
          eta=2.0*eta/(sgd*pi)
          L11=0.0;L21=0.0;L12=0.0;L22=0.0
          do i=1,nf
            r2     = exp(-beta*abs(w(i)))
            dfermi = -beta*r2/((1.0+r2)**2)     !
            L11=L11+dw(i)*dfermi*tau(i)         !         /    df         i+j-2
            L12=L12+dw(i)*dfermi*tau(i)*w(i)    ! L_ij = -| dw -- tau(w) w
            L21=L21+dw(i)*dfermi*tau(i)*w(i)    !         /    dw
            L22=L22+dw(i)*dfermi*tau(i)*w(i)**2 ! 
            write(45,"(1x,f9.4,2x,f9.4,2x,f9.4)") w(i),tau(i),
     &        (-1.0/pi)*aimag(1.0/(w(i) - ed - Epsk(0.0,pi)-blint(0.0,pi,i)))
          end do
          thermopower=-beta*L12/L11    ! S = (−1/eT )(L12/L11)
          kappa= beta*(L22-L12**2/L11) ! Thermal conductivity
          peltier=L21/L11              ! Peltier Coef.
          write(44,"(f7.4,10(1x,f10.5))") 
     &          temp,res,thermopower,kappa,peltier,hall,eta,L11,L12,L21,L22
       
        endif ! (iconductivity.eq.1)
        
c*******************************************************************************        
c       wrap up
c*******************************************************************************        
        call deallocate_1p
        call deallocate_Spec0
        call deallocate_Spec1      
        stop
        end program spectra
                  
        
        
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


