        program spectra_star
c*******************************************************************************
c       This code reads in the real-frequency self energy data
c       extracted by dsigcluster from the MEM A(k,w), and then
c       calculates the interpolated A(k,w) on the lattice.  The
c       interpolation of the self enrgy is done using a star
c       interpolation.  
c
c       To compile type make spectra_star (or make all) using the 
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
     &             nkpoints, info, infot, ikx, iky, deg, ios

        real(kind) :: r1, r2, r3, kx, ky, kz, Epsk, kx1, ky1, kz1,
     &                kx2, ky2, kz2, offset, res, hall, vx, vyy, dfermi
        
        complex(kind) :: c1, c2, blint
        
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


c       Read some parameters from the QMC sigma file
 401    read(20,"(a72)") linemc
        if(index(linemc(1:72),'Results').ne.0 ) then
          read(20,*) nl,nwn,fill,cluster,ntot,ndim
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
          read(20,*) ed,U,beta,i,run,meas,tprime
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
        
        do ick=1,Ncw  ! do all K in the irreducible wedge
c         Open the self energy files.
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
          
c         Read in the real-frequency self energy
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

c       Make sure that nf is le the limit set above        
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


c       Now do the interpolation, and calculate the spectrum at
c       a sequence of locations.

        if(ndim.eq.2) then
          write(6,"('Enter kx1,ky1,kx2,ky2,offset,nkpoints:  ',$)") 
          read(5,*) kx1,ky1,kx2,ky2,offset,nkpoints
        else if(ndim.eq.3) then
          write(6,"('Enter kx1,ky1,kz1,kx2,ky2,kz2,offset,nkpoints:  ',$)") 
          read(5,*) kx1,ky1,kz1,kx2,ky2,kz2,offset,nkpoints
        end if
        
        open(unit=22,file='Spectra',status='unknown')
        open(unit=21,file='Akw_spectra',status='unknown')
        open(unit=24,file='False_color.dat',status='unknown')
        write(6,*) '********************************************'
        write(6,*) '*** G(k,w) and Sigma(k,w) are in Spectra ***'
       write(6,*) '*** A(k,w) is in Akw_spectra             ***'
        write(22,"('#   w          ReGkw        ImGkw       ReSigkw  ',
     &             '    ImSigkw')")
        write(21,"('#   w      A(k,w)=-Im G(k,w)/pi ')")
        write(24,"('# k     w      A(k,w)=-Im G(k,w)/pi ')")
        if(ndim.eq.2) then
          do ick=1,nkpoints ! nkpoints may be 1
            kx=kx1 + real(ick-1)*(kx2-kx1)/max(nkpoints-1,1)
            ky=ky1 + real(ick-1)*(ky2-ky1)/max(nkpoints-1,1)

            do i=1,nf
              c2=0.5*(blint(kx,ky,zeror,i)+blint(ky,kx,zeror,i))
              c1=1.0/(w(i)-ed-Epsk(kx,ky,zeror)-c2)
              write(22,"(5(g12.6,1x))") w(i),real(c1),aimag(c1),
     &                             real(c2),aimag(c2)
              write(21,"(5(g12.6,1x))") w(i),
     &                              -aimag(c1)/pi+(ick-1)*offset
        
              write(24,"(5(g12.6,1x))") ick,w(i),
     &                              -aimag(c1)/pi

            end do
            write(21,*) '& '
            write(22,*) '& '
            write(24,*) ' '
          end do
        else if(ndim.eq.3) then
          do ick=1,nkpoints ! nkpoints may be 1
            kx=kx1 + real(ick-1)*(kx2-kx1)/max(nkpoints-1,1)
            ky=ky1 + real(ick-1)*(ky2-ky1)/max(nkpoints-1,1)
            kz=kz1 + real(ick-1)*(kz2-kz1)/max(nkpoints-1,1)

            do i=1,nf
              c2=blint(kx,ky,kz,i)
              c1=1.0/(w(i)-ed-Epsk(kx,ky,kz)-c2)
              write(22,"(5(g12.6,1x))") w(i),real(c1),aimag(c1),
     &                             real(c2),aimag(c2)
              write(21,"(5(g12.6,1x))") w(i),-aimag(c1)/pi+(ick-1)*offset
        
              write(24,"(5(g12.6,1x))") ick,w(i),-aimag(c1)/pi

            end do
            write(21,*) '& '
            write(22,*) '& '
            write(24,*) ' '
          end do
        end if
      
 
      
        write(6,*) '********************************************'
        if (ndim.eq.2) then
c         Now that we have the interpolated spectra, calculate the 
c         transport using the lifetime approximation (which is only
c         exact when D=oo, or Nc=1 or 4).
          open(unit=44,file='resis.dat',status='unknown')

          res=0.0
          hall=0.0
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
            r1 = Epsk(kx,ky,zeror)
            do i=1,nf
              r2     = exp(-beta*abs(w(i)))
              dfermi = -beta*r2/((1.0+r2)**2)
              c2     = 0.5*(blint(kx,ky,zeror,i)+ blint(ky,kx,zeror,i))
              r1     = aimag(1.0/(w(i)-ed-r1-c2))
              res=res -
     &           dw(i)*deg*dfermi*(r1*vx)**2
              hall=hall +
     &           dw(i)*deg*dfermi*(vx**2)*vyy*(r1)**3
            end do
          end do
          end do
          res=(float(nover**2)*pi**2)/res         ! units e^2 pi/(2 hbar a)
          hall=hall/(float(nover**2)*pi**3)       ! 2 pi^2 e^3 a B /(3 hbar^2)
          hall=hall/res**2
          write(44,*) temp,res,hall
        endif
        
c*******************************************************************************
c       wrap up
c*******************************************************************************        
        call deallocate_1p
        call deallocate_Spec0
        stop
        end program spectra_star
                  
        
        
        complex(8) function blint(kx,ky,kz,n)
c       
c       This subroutine performs a bilinear interpolation of
c       the self energy, at any location (kx,ky) in the zone.
c
        use Global
        use module_spectra
        implicit none
        integer :: n,i1,j1,i2,j2,ic1,ic2,ic3,ic4,ic,ik
        real(kind) :: Kx,Ky,Kz,dx,dy,g1,g2,r7,r8,Kx3,Ky3
        complex(kind) :: z1,z2,z3,z4,Sapprox(Nc)
        
c       interpolation by introducing a factor exp(i kt X) in sigma
        Sapprox=(0.d0,0.d0)
        do ic=1,Nc
        do ik=1,Nc        
          Sapprox(ic)= Sapprox(ic)
     &          +FTCoefs_K_to_R(ik,ic)*Sigma(n,ik)
        enddo
        enddo

        blint=(0.d0,0.d0)
        do ic=1,Nc                              
           z1=exp(+ii*(kx*Rc(1,ic)+ky*Rc(2,ic)+kz*Rc(3,ic)))
           blint=blint+z1*Sapprox(ic)
        enddo
        return
        end
                
        
        real*8 function Epsk(kx,ky,kz)
        use Global
        implicit none
        real(kind) :: kx, ky, kz  
        
        if(ndim.eq.2) then
          Epsk=-0.5_kind*(cos(kx)+cos(ky)) -
     &       tprime*(cos(kx)*cos(ky)-1.0_kind)
        else if(ndim.eq.3) then
          Epsk=-0.5_kind*(cos(kx)+cos(ky)+cos(ky)) 
     &        - tprime*(cos(kx)*cos(ky)
     &                 +cos(kx)*cos(kz)
     &                 +cos(ky)*cos(kz))
        
        end if
        
        end function Epsk
        
        
