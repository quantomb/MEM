        program spectra_star
c*******************************************************************************
c       This code reads in the real-frequency self energy data
c       extracted by dsigcluster from the MEM A(k,w), and then
c       calculates the interpolated A(k,w) on the lattice.  The
c       interpolation os the self enrgy is done using a star
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

        integer, parameter :: nover=128
        integer :: filenum, istat, i, j, n, run, meas, ick, ic, ik, 
     &             nkpoints, info, infot, ikx, iky, deg, ios,
     &             npad1,npad2,npad3   

        real(kind) :: r1, r2, r3, kx, ky, kz, Epsk, kx1, ky1, kz1, resclk,
     &                kx2, ky2, kz2, kp1(3), offset, res, hall, vx, vyy, dfermi,
     &                ksym(3,8)
        
        complex(kind) :: c1, c2, c22, blint, inter_sigmaf
        complex(kind) :: cn0(24),cn1,cn2,cn3,ccn
	complex(kind) :: blintG, blintlng, blintLS, blint2, blintcmlnt
        
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
     

c*******************************************************************************
c       Write out the K information
c*******************************************************************************
        
        if(ndim.eq.2) then
          write(6,*) ' ick  Kcx(ick)   Kcy(ick)'
          do ick=1,Nc
            write(6,"(1x,i3,1x,f9.6,2x,f9.6)") ick,Kc(1,ick),Kc(2,ick)
          end do
        else if(ndim.eq.3) then
          write(6,*) ' ick  Kcx(ick)   Kcy(ick)   Kcz(ick)'
          do ick=1,Nc
             write(6,"(1x,i3,1x,f9.6,2x,f9.6,2x,f9.6)") 
     &		   ick,Kc(1,ick),Kc(2,ick),Kc(3,ick)
          end do
        end if

   
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
        write(24,"('# k  rescale-k   w      A(k,w)=-Im G(k,w)/pi ')")
        if(ndim.eq.2) then
          do ick=1,nkpoints ! nkpoints may be 1
            kx=kx1 + real(ick-1)*(kx2-kx1)/max(nkpoints-1,1)
            ky=ky1 + real(ick-1)*(ky2-ky1)/max(nkpoints-1,1)

            do i=1,nf
c              c2=0.5*(blintcmlnt(kx,ky,zeror,i)+blintcmlnt(ky,kx,zeror,i))
              c2=0.125*(blintcmlnt(kx,ky,zeror,i)  + blintcmlnt(-kx,-ky,zeror,i)  ! Sigma(k,w) with C_4
     &                +blintcmlnt(-kx,ky,zeror,i) + blintcmlnt(kx,-ky,zeror,i)   ! symmetry imposed
     &                +blintcmlnt(ky,kx,zeror,i)  + blintcmlnt(-ky,-kx,zeror,i)
     &                +blintcmlnt(-ky,kx,zeror,i) + blintcmlnt(ky,-kx,zeror,i)  )

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

c         set symmetrize table for kx,ky,kz
          ksym(1,1)= 1.0; ksym(2,1)= 1.0; ksym(3,1)= 1.0
          ksym(1,2)=-1.0; ksym(2,2)= 1.0; ksym(3,2)= 1.0
          ksym(1,3)=-1.0; ksym(2,3)=-1.0; ksym(3,3)= 1.0
          ksym(1,4)=-1.0; ksym(2,4)=-1.0; ksym(3,4)=-1.0
          ksym(1,5)=-1.0; ksym(2,5)= 1.0; ksym(3,5)=-1.0
          ksym(1,6)= 1.0; ksym(2,6)= 1.0; ksym(3,6)=-1.0
          ksym(1,7)= 1.0; ksym(2,7)=-1.0; ksym(3,7)=-1.0
          ksym(1,8)= 1.0; ksym(2,8)=-1.0; ksym(3,8)= 1.0
         

          do ick=1,nkpoints ! nkpoints may be 1
            kx=kx1 + real(ick-1)*(kx2-kx1)/max(nkpoints-1,1)
            ky=ky1 + real(ick-1)*(ky2-ky1)/max(nkpoints-1,1)
            kz=kz1 + real(ick-1)*(kz2-kz1)/max(nkpoints-1,1)



            do i=1,nf
             c1=(0.0d0,0.0d0)
             c22=(0.0d0,0.0d0)
             do j=1,8

c              c22=c22+blint(kx*ksym(1,j),ky*ksym(2,j),kz*ksym(3,j),i) 
              c22=c22+blintcmlnt(kx*ksym(1,j),ky*ksym(2,j),kz*ksym(3,j),i)         

c              c2=blint(kx,ky,kz,i)

c              if(w(i).lt.0.0)then
c                c2=blintcmlnt(kx,ky,kz,i)
c              else
c                c2=blint(kx,ky,kz,i)
c              endif
c              c2=blintLS(kx,ky,kz,i)
c              c2=blintG(kx,ky,kz,i)

c           if((aimag(c2).gt.-0.18).and.(w(i).lt.0.8).and.(w(i).gt.-1.2)) then
c              if ( (aimag(c2).gt.-0.1) .and. 
c     &         abs(w(i)-ed-Epsk(kx,ky,kz)-real(c2)).lt.0.1 ) then
c            check subroutine interpol_sigma to see wether it use Sigma or G. 
c                call interpol_sigma(kp1,inter_sigmaf,i)
c                c2=inter_sigmaf
c                c2=real(c2)+ii*aimag(inter_sigmaf)
c                c2=real(c2)+ii*(0.0-0.1*exp((-aimag(c2)-0.1)*4.0))
c                 c2=real(c2)+ii*(-0.18+(aimag(c2)+0.18)*(0.18-0.01)/(0.18+0.36))
c              endif
              
              enddo ! j=1,8
             c2=c22/8.0d0
c              c1=1.0/(w(i)-ed-Epsk(kx,ky,kz)-c2)
               c1=1.0/(1.0/c2-ed-Epsk(kx,ky,kz))       
         
c              if(w(i).lt.0.0)then
c                 c1=1.0/(1.0/c2-ed-Epsk(kx,ky,kz))
c              else            
c                c1=1.0/(w(i)-ed-Epsk(kx,ky,kz)-c2)
c              endif
               
c              c1=c2
             

c              If (aimag(c1).lt.0.0) then
              write(22,"(5(g12.6,1x))") w(i),real(c1),aimag(c1),
     &                             real(c2),aimag(c2)
              write(21,"(5(g12.6,1x))") w(i),-aimag(c1)/pi+(ick-1)*offset
              

c         column 2nd is for rescale ick to k_F (here k_F=1.12353)
              resclk=dble(ick-1)/dble(nkpoints)*pi/1.12353        
              write(24,"(4(g12.6,1x))") ick,resclk,w(i),-aimag(c1)/pi
c              endif

            end do
            write(21,*) '& '
            write(22,*) '& '
            write(24,*) ' '
          end do
        end if
      
 
      
c        write(6,*) '********************************************'
c        if (ndim.eq.2) then
cc         Now that we have the interpolated spectra, calculate the 
cc         transport using the lifetime approximation (which is only
cc         exact when D=oo, or Nc=1 or 4).
c          open(unit=44,file='resis.dat',status='unknown')
c          res=0.0
c          hall=0.0
c          do ikx=0,nover  ! restrict the calculation to the IRW.
c          do iky=0,ikx
c            if (ikx.eq.iky.and.ikx.eq.0
c     &               .or.
c     &          ikx.eq.iky.and.ikx.eq.nover) then
c              deg=1
c            else if(ikx.eq.iky.or.iky.eq.0.or.ikx.eq.nover) then
c              deg=4          
c            else
c              deg=8
c            end if
c            kx=ikx*pi/float(nover)
c            ky=iky*pi/float(nover)
c            vx = sin(kx)*(0.5 + tprime*cos(ky))
c            vyy = cos(ky)*(0.5 + tprime*cos(kx))
c            r1 = Epsk(kx,ky,zeror)
c            do i=1,nf
c              r2     = exp(-beta*abs(w(i)))
c              dfermi = -beta*r2/((1.0+r2)**2)
c              c2     = 0.5*(blint(kx,ky,zeror,i)+ blint(ky,kx,zeror,i))
c              r1     = aimag(1.0/(w(i)-ed-r1-c2))
c              res=res -
c     &           dw(i)*deg*dfermi*(r1*vx)**2
c              hall=hall +
c     &           dw(i)*deg*dfermi*(vx**2)*vyy*(r1)**3
c            end do
c          end do
c          end do
c          res=(float(nover**2)*pi**2)/res         ! units e^2 pi/(2 hbar a)
c          hall=hall/(float(nover**2)*pi**3)       ! 2 pi^2 e^3 a B /(3 hbar^2)
c          hall=hall/res**2
c          write(44,*) temp,res,hall
c        endif
        
c*******************************************************************************
c       wrap up
c*******************************************************************************        
        call deallocate_1p
        call deallocate_Spec0
        stop
        end program spectra_star                    
        
        complex(8) function blint(kx,ky,kz,n)
c       
c       This subroutine performs a Fourier interpolation of
c       the self energy, at any location (kx,ky,kz) in the zone.
c
        use Global
        use module_spectra
        implicit none
        integer :: n,i1,j1,i2,j2,ic1,ic2,ic3,ic4,ic,ik
        real(kind) :: Kx,Ky,Kz,dx,dy,g1,g2,r7,r8,Kx3,Ky3
        real(kind) :: resig,imsig
        complex(kind) :: z1,z2,z3,z4,Sapprox(Nc)
                
c       interpolation by introducing a factor exp(i kt X) in sigma
        Sapprox=(0.d0,0.d0)
        do ic=1,Nc
        do ik=1,Nc        
          Sapprox(ic)= Sapprox(ic)
     &          +FTCoefs_K_to_R(ik,ic)*Sigma(n,ik)
        enddo
        enddo

        blint=(0.0d0,0.0d0)
        do ic=1,Nc                              
           z1=exp(+ii*(kx*Rc(1,ic)+ky*Rc(2,ic)+kz*Rc(3,ic)))
           blint=blint+z1*Sapprox(ic)
        enddo
        return
        end



        complex(8) function blintcmlnt(kx,ky,kz,n)
c       
c       This subroutine performs a Fourier interpolation of
c       the self energy, at any location (kx,ky,kz) in the zone.
c
        use Global
        use module_spectra
        implicit none
        integer :: n,i1,j1,i2,j2,ic1,ic2,ic3,ic4,ic,ik
        real(kind) :: Kx,Ky,Kz,dx,dy,g1,g2,r7,r8,Kx3,Ky3
        real(kind) :: resig,imsig
        complex(kind) :: z1,cmlnt,Sapprox(Nc)
                
c       interpolation by introducing a factor exp(i kt X) in sigma
        Sapprox=(0.d0,0.d0)
        do ic=1,Nc
        do ik=1,Nc
          cmlnt=1.0/(w(n)-Sigma(n,ik))        
          Sapprox(ic)= Sapprox(ic)
     &          +FTCoefs_K_to_R(ik,ic)*cmlnt
        enddo
        enddo

        blintcmlnt=(0.0d0,0.0d0)
        do ic=1,Nc       
             z1=exp(+ii*(kx*Rc(1,ic)+ky*Rc(2,ic)+kz*Rc(3,ic)))
             blintcmlnt=blintcmlnt+z1*Sapprox(ic)
        enddo
        return
        end




c===============================================================================

        complex(8) function blintspt(kx,ky,kz,n)
c       
c       This subroutine performs a Fourier interpolation of
c       the self energy, at any location (kx,ky,kz) in the zone.
c
        use Global
        use module_spectra
        implicit none
        integer :: n,i1,j1,i2,j2,ic1,ic2,ic3,ic4,ic,ik
        real(kind) :: Kx,Ky,Kz,dx,dy,g1,g2,r7,r8,Kx3,Ky3
        real(kind) :: resig,imsig
        complex(kind) :: z1,z2,z3,z4,Sapprox(Nc)
        
        
c       interpolation by introducing a factor exp(i kt X) in sigma
        Sapprox=(0.d0,0.d0)
        do ic=1,Nc
        resig=0.0d0
        imsig=0.0d0
        do ik=1,Nc        
c          Sapprox(ic)= Sapprox(ic)
c     &          +FTCoefs_K_to_R(ik,ic)*Sigma(n,ik)
        resig=resig
     &      +real(FTCoefs_K_to_R(ik,ic))*real(Sigma(n,ik))
        imsig=imsig
     &      +real(FTCoefs_K_to_R(ik,ic))*aimag(Sigma(n,ik))
c        write(234,*)'ic=',ic,'ik=',ik
c        write(234,*)'real(FTCoefs_K_to_R(ik,ic))=',real(FTCoefs_K_to_R(ik,ic))
c        write(234,*)'aimag(FTCoefs_K_to_R(ik,ic))=',aimag(FTCoefs_K_to_R(ik,ic))
c        write(234,*)'real(Sigma(n,ik))=',real(Sigma(n,ik))
c        write(234,*)'aimag(Sigma(n,ik)=',aimag(Sigma(n,ik))
        enddo
         Sapprox(ic)=resig+ii*imsig
        enddo

        blintspt=(0.d0,0.d0)
        do ic=1,Nc                              
c           z1=exp(+ii*(kx*Rc(1,ic)+ky*Rc(2,ic)+kz*Rc(3,ic)))
c           blintspt=blintspt+z1*Sapprox(ic)
c           z1=exp(+ii*(kx*Rc(1,ic)+ky*Rc(2,ic)+kz*Rc(3,ic)))
           z2=cos(kx*Rc(1,ic)+ky*Rc(2,ic)+kz*Rc(3,ic))
           z3=cos(kx*Rc(1,ic)+ky*Rc(2,ic)+kz*Rc(3,ic))
           blintspt=blintspt+z2*real(Sapprox(ic))
     &          +ii*z3*aimag(Sapprox(ic))
        enddo
        return
        end

        complex(8) function blintLS(kx,ky,kz,n)
c       
c       This subroutine performs a Fourier interpolation of
c       the self energy, using the basis has the lattice group symmetry
c       at any location (kx,ky,kz) in the zone.
c   
        use Global
        use module_spectra
        implicit none
        integer :: n,i1,j1,i2,j2,ic1,ic2,ic3,ic4,ic,ik,igroup,igo
        real(kind) :: Kx,Ky,Kz,dx,dy,g1,g2,r7,r8,Kx3,Ky3
        complex(kind) :: z1,z2,z3
        complex(kind) :: Sapprxls(Rgonum)
        
c       interpolation by introducing a factor exp(i kt X) in sigma
        Sapprxls=(0.d0,0.d0)
        
        do ik=1,Rgonum
          z2=(0.d0,0.d0) 
          igo=0
          do ic=1,Ncw
          do igroup=1,48
            igo=igo+1
            z3=exp( -ii*(Kgo(1,igo)*Rgo(1,ik)+Kgo(2,igo)*Rgo(2,ik)
     & +Kgo(3,igo)*Rgo(3,ik)) )     
            z2=z2 + z3*Sigma(n,ic)/real(Kgodeg(igo),kind)         
          enddo
          enddo
          Sapprxls(ik)=z2
        enddo

        blintLS=(0.d0,0.d0)
        do ic=1,Rgonum                              
           z1=exp(+ii*(kx*Rgo(1,ic)+ky*Rgo(2,ic)+kz*Rgo(3,ic)))
           blintLS=blintLS+z1*Sapprxls(ic)
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
          Epsk=-0.5_kind*(cos(kx)+cos(ky)+cos(kz)) 
     &        - tprime*(cos(kx)*cos(ky)
     &                 +cos(kx)*cos(kz)
     &                 +cos(ky)*cos(kz))
        
        end if
        
        end function Epsk
        
        
!===================================================================

      subroutine interpol_sigma(kp1,inter_sigmaf,nwi)
!  find the volume of tetrahedral
      use Global
      use module_spectra
      implicit none

      real(kind), intent(in)  :: kp1(3)
	integer, intent(in)  :: nwi
      complex(kind), intent(out) :: inter_sigmaf

      real(kind) :: vol(4)
      real(kind)     :: tetra(3,Nc)   ! tetrahedron : tetra(4,:) is polarization
      integer :: tetra_Nc(Nc), j
      

       call find_tetra(kp1,tetra,tetra_Nc)

c#ifdef USE_TETRA_INTERPOL
    
       call cal_vol(kp1,tetra,vol)

c      inter_sigmaf(:,:,:) = vol(1)*sigmaf(:,:,tetra_Nc(1),:) + vol(2)*sigmaf(:,:,tetra_Nc(2),:) &
c                       + vol(3)*sigmaf(:,:,tetra_Nc(3),:) + vol(4)*sigmaf(:,:,tetra_Nc(4),:)

c       inter_sigmaf=  vol(1)*GKw(nwi,tetra_Nc(1)) 
c     &               +vol(2)*GKw(nwi,tetra_Nc(2)) 
c     &               +vol(3)*GKw(nwi,tetra_Nc(3))  
c     &		      +vol(4)*GKw(nwi,tetra_Nc(4))

       inter_sigmaf=  vol(1)*Sigma(nwi,tetra_Nc(1)) 
     &               +vol(2)*Sigma(nwi,tetra_Nc(2)) 
     &               +vol(3)*Sigma(nwi,tetra_Nc(3))  
     &  	     +vol(4)*Sigma(nwi,tetra_Nc(4))


c#else
      inter_sigmaf= Sigma(nwi,tetra_Nc(1))
c#endif

      RETURN
      END subroutine interpol_sigma

!-------------------------------------------------------

      subroutine find_tetra(kp1,tetra,tetra_Nc)
!  find the appropriate tetrahedron for (Kx,Ky,Kz)
      use Global
      use module_spectra
      implicit none

      real(kind)    , intent(in)  :: kp1(3)
      real(kind)    , intent(out) :: tetra(3,Nc)
      integer, intent(out) :: tetra_Nc(Nc)
 
      real(kind) :: area, volume

      integer :: leng_ic(7*Nc), min_i, leng_i
      real(kind)     :: leng(7*Nc), leng_K(3,7*Nc), min_leng
      integer :: ic, i, j, k, l, itest

!******************************************
! Principle lattice vectors in K-space
!******************************************
      real(kind)     :: K0(3,4)

      K0(:,1)=(/ 0.d0   , 0.d0   , 0.d0   /)
      K0(:,2)=(/-2.*pi  , 0.d0   , 0.d0   /)
      K0(:,3)=(/ 0.d0   , -2.*pi , 0.d0   /)
      K0(:,4)=(/ 0.d0   , 0.d0   , -2.*pi /)


      leng_i = 0

      DO ic=1,Nc
      do i=1,4
         leng_i = leng_i + 1
         leng_ic(leng_i) = ic
         leng_K(:,leng_i) = Kc(:,ic)+K0(:,i)
         leng(leng_i) = sum((leng_K(:,leng_i)-kp1(:))**2)
      enddo
      do i=2,4
         leng_i = leng_i + 1
         leng_ic(leng_i) = ic
         leng_K(:,leng_i) = Kc(:,ic)-K0(:,i)
         leng(leng_i) = sum((leng_K(:,leng_i)-kp1(:))**2)
      enddo
      ENDDO

      j = 0
      tetra1 : DO
      min_i    = 1
      min_leng = leng(1)
      DO i=2,7*Nc
         if ( leng(i) < min_leng ) then
            min_i    = i
            min_leng = leng(i)
         endif
      ENDDO
      leng(min_i) = 1.d30
      itest = 0
      tetra2 : do k=1,j-1
      do l=k+1,j
         if ( area(leng_K(1:3,min_i),tetra,k,l) < eps ) then
            itest = 1
            exit tetra2
         endif
      enddo
      enddo tetra2
      if (itest == 0 .and. 
     & (j<3 .or. volume(leng_K(1:3,min_i),tetra,1,2,3)>eps )) then
         j = j + 1
         tetra(1:3,j) = leng_K(1:3,min_i)
         tetra_Nc(j) = leng_ic(min_i)
         if (j >= 4) exit tetra1
      endif
      ENDDO tetra1

      RETURN
      END subroutine find_tetra

!-----------------------------------------

      subroutine cal_vol(kp1,tetra,vol)
!  find the appropriate cluster for (Kx,Ky,Kz)

      use Global
      use module_spectra
      implicit none

      real(kind) :: volume      

      real(kind), intent(in)  :: kp1(3)
      real(kind), intent(in)  :: tetra(3,Nc)
      real(kind), intent(out) :: vol(4)

      real(kind)     :: vol0

      vol(1) = volume(kp1,tetra,2,3,4)
      vol(2) = volume(kp1,tetra,3,4,1)
      vol(3) = volume(kp1,tetra,4,1,2)
      vol(4) = volume(kp1,tetra,1,2,3)

      vol0 = sum(vol(:))
      vol(:) = vol(:)/vol0

      RETURN
      END subroutine cal_vol

!-----------------------------------------

      real(kind) FUNCTION volume(kp1,tetra,ic1,ic2,ic3)
!  find the volume of tetrahedral

      use Global
      use module_spectra
      implicit none

      real(kind)    , intent(in) :: tetra(3,Nc), kp1(3)
      integer, intent(in) :: ic1, ic2, ic3
      real(kind) :: Pa(3), Pb(3), Pc(3)

      Pa(1:3) = tetra(1:3,ic1)-kp1(1:3)
      Pb(1:3) = tetra(1:3,ic2)-kp1(1:3)
      Pc(1:3) = tetra(1:3,ic3)-kp1(1:3)

      volume = Pa(1)*(Pb(2)*Pc(3)-Pb(3)*Pc(2)) 
     &     - Pa(2)*(Pb(1)*Pc(3)-Pb(3)*Pc(1)) 
     &     + Pa(3)*(Pb(1)*Pc(2)-Pb(2)*Pc(1))

      volume = abs(volume) / 6.d0

      RETURN
      END FUNCTION volume

!-------------------------------------------------------

      real(kind) FUNCTION area(kp1,tetra,ic1,ic2)
!  find the volume of tetrahedral

      use Global
      use module_spectra
      implicit none

      real(kind)    , intent(in) :: tetra(3,Nc), kp1(3)
      integer, intent(in) :: ic1, ic2
      real(kind) :: leng_side(3), sum_leng

      leng_side(1) = sqrt(sum((tetra(1:3,ic1)-kp1(1:3))**2))
      leng_side(2) = sqrt(sum((tetra(1:3,ic2)-kp1(1:3))**2))
      leng_side(3) = sqrt(sum((tetra(1:3,ic2)-tetra(1:3,ic1))**2))

      sum_leng = 0.5d0*sum(leng_side(:))

      area = sqrt(sum_leng*(sum_leng-leng_side(1)) 
     &                  *(sum_leng-leng_side(2)) 
     &                  *(sum_leng-leng_side(3)))

      RETURN
      END FUNCTION area


