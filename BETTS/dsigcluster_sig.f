        program dsigcluster
c*******************************************************************************
c       This block is designed to take the DCA Hubbard model single-particle
c       A_bar(K,w) and extract from it the self-energy Sigma(K,w) by 
c       inverting the coarse-graining step.  
c       This block was designed to run on Linux boxes.
c*******************************************************************************
        use Global
        implicit none
c*******************************************************************************
c       Some local definitions
c
c       Characters
        character*72 linemc

c       Integers
        integer, parameter :: itermax=1000, npolish=5, nfmax=800
        integer :: i,j,ick,run,meas,iter,nw,ised,nguess,iguess,
     &             ipolish,ios
     
c       Reals        
        real(kind), parameter :: delta=0.001_kind
        real(kind) :: Epsk,w(nfmax),dw(nfmax),Arf(nfmax),
     &                Gr(nfmax),Gi(nfmax),Kx,Ky,Kz,tolerance,
     &                damping,shift,chcrt0,fill2

c       Complex
        complex(kind) :: q,Az,Bz,Cz,z1,z2,z3,znew,denom,P1,P2,P3,
     &                   sigma(nfmax),Grf(nfmax),lastsig
c*******************************************************************************
        
c*******************************************************************************
c       Open some files
c*******************************************************************************

c       Open a datafile for the input spectra
        write(6,"('enter the QMC spectra filename:  ',$)") 
        read(5,"(a72)") linemc
        open(unit=9,file=linemc,status='old')

c       Open a datafile for the parameters
        write(6,"('  enter the QMC sigma filename:  ',$)") 
        read(5,"(a72)") linemc
        open(unit=92,file=linemc,status='old')
c 1000   format(a72)
 
c       Open a datafile for the output Sigma(w)
        open(unit=67,file='Sigma.dat',status='unknown')
        write(67,"('#      w   Re[Sigma(K,w)] Im[Sigma(K,w)]   A(K,w)        dw  ')")

c*******************************************************************************
c       Read in some parameters
c*******************************************************************************
c       Read some parameters from the QMC sigma file
 401    read(92,"(a72)") linemc
        if(index(linemc(1:72),'Results').ne.0 ) then
         read(92,*) nl,nwn,fill,cluster,ntot,ndim,nw,chcrt0
         fill2=halfr*fill
        else
          goto 401
        endif

c  in order to have better resolution for the sharp quasi-particle peak in the FL
c  case, temporally increase ntot from 100 to 200.
c        ntot=200

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


 402    read(92,"(a72)") linemc
        if(index(linemc(1:72),'tprime').ne.0 ) then
           if(index(linemc(1:72),'tdprime').ne.0 ) then
              if(index(linemc(1:72),'t3prime').ne.0 ) then
                 read(92,*) ed,U,beta,i,run,meas,tprime,tdprime,t3prime
              else
                 read(92,*) ed,U,beta,i,run,meas,tprime,tdprime
                 t3prime=0.d0
              endif
           else 
              read(92,*) ed,U,beta,i,run,meas,tprime
              tdprime=0.d0
              t3prime=0.d0
           endif
        else
          goto 402
        endif
        temp=1.0/beta   
        
c*******************************************************************************
c       Allocate the arrays and Create the tables
c*******************************************************************************
        call allocate_1p
        call tables
        
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
            write(6,"(1x,i3,1x,f9.6,2x,f9.6,2x,f9.6)") ick,Kc(1,ick),Kc(2,ick),Kc(3,ick)
          end do
        end if

c*******************************************************************************
c       Read in a few more things
c*******************************************************************************
        write(6,"('             Ncw= ',i2,' enter ick:  ',$)") Ncw
        read(5,*) ick

        
        
c       Read in the sigma(w) and the frequency grid, and form Sigma(i)
        do i=1,1000
          read(9,*,end=15) w(i),Arf(i),dw(i)
          Arf(i)=Arf(i)*U**2*fill2*(oner-fill2)
        end do
 15     continue
        nw=i-1




c*******************************************************************************
c       Now integrate to find the real part of the self energy Sigma(ick,w), G means Sig in the following.
c*******************************************************************************
        Gi(:)=-pi*Arf(:)
        do i=2,nw-1
          Gr(i)=0.0_kind
          do j=1,i-1
            Gr(i)=Gr(i) + dw(j)*Arf(j)/(w(i)-w(j))
          end do
          do j=i+1,nw
            Gr(i)=Gr(i) + dw(j)*Arf(j)/(w(i)-w(j))
          end do
        end do

c*******************************************************************************
c       Add the Hartree term to the real part of the self energy Sigma(ick,w)
c*******************************************************************************
        Gr=Gr+U*(fill2-halfr)
        sigma(:)=cmplx(Gr(:),Gi(:))

c*******************************************************************************
c       Calculate A(ick,w) as a sanity check
c*******************************************************************************

c        write(5100,*) "ed,Epsk(zeror,zeror,zeror)=",ed,Epsk(zeror,zeror,zeror)
c        write(5100,*) "U,fill2=",U,fill2
        if(Nc.eq.1) then
          Arf(:)=-aimag(
     &            1.0_kind/(w(:)-ed-Epsk(zeror,zeror,zeror)-sigma(:))
     &            )/pi
        else
          Arf(:)=-aimag(
     &            1.0_kind/(w(:)-ed-Epsk(Kc(1,ick),Kc(2,ick),Kc(3,ick))-sigma(:))
     &            )/pi

        end if


c*******************************************************************************
c       Write everything out to Sigma.dat
c*******************************************************************************

        do i=1,nw
          write(67,"(f11.6,2x,f11.6,2x,f11.6,2x,f11.6,2x,f11.6)") 
     &      w(i),Gr(i),Gi(i),Arf(i),dw(i)
        end do

        stop
        end program dsigcluster


        real(8) function Epsk(kx,ky,kz)
c*******************************************************************************
        use Global
        implicit none
c*******************************************************************************
        real(8) kx,ky,kz
c*******************************************************************************
        if(ndim.eq.2) then
          Epsk=-0.5_kind*(cos(kx)+cos(ky)) -
     &       tprime*(cos(kx)*cos(ky)-1.0_kind)
     &             - 0.5*tdprime*(cos(2.d0*kx)+cos(2.d0*ky))
     &             - t3prime*(cos(2.d0*kx)*cos(ky)
     &                         +cos(kx)*cos(2.d0*ky))   


        else if(ndim.eq.3) then
          Epsk=-0.5_kind*(cos(kx)+cos(ky)+cos(kz))
     &        - tprime*(cos(kx)*cos(ky)
     &                 +cos(kx)*cos(kz)
     &                 +cos(ky)*cos(kz))
        
        end if

        return
        end
