        program dsigcluster
c*******************************************************************************
c       This block is designed to take the DCA Hubbard model single-particle
c       A_bar(K,w) and extract from it the self-energy Sigma_bar(K,w) by 
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
        integer, parameter :: itermax=2000, npolish=10, nfmax=800
        integer :: i,j,ick,run,meas,iter,nw,ised,nguess,iguess,
     &             ipolish
     
c       Reals        
        real(kind), parameter :: delta=0.01_kind
        real(kind) :: Epsk,fill,w(nfmax),dw(nfmax),Arf(nfmax),
     &                Gr(nfmax),Gi(nfmax),Kx,Ky,tolerance,
     &                damping,shift

c       Complex
        complex(kind) :: q,Az,Bz,Cz,z1,z2,z3,znew,denom,P1,P2,P3,
     &                   sigma(nfmax),Grf(nfmax),func,lastsig
c*******************************************************************************
        
c*******************************************************************************
c       Open some files
c*******************************************************************************

c       Open a datafile for the input spectra
        write(6,1) 
 1      format('enter the QMC spectra filename:  ',$)
        read(5,1000) linemc
        open(unit=9,file=linemc,status='old')

c       Open a datafile for the parameters
        write(6,2) 
 2      format('  enter the QMC sigma filename:  ',$)
        read(5,1000) linemc
        open(unit=92,file=linemc,status='old')
 1000   format(a72)
 
c       Open a datafile for the output Sigma(w)
        open(unit=67,file='Sigma.dat',status='unknown')
 
c*******************************************************************************
c       Read in some parameters
c*******************************************************************************
c       Read some parameters from the QMC sigma file
 401    read(92,1000) linemc
        if(index(linemc(1:72),'Results').ne.0 ) then
          read(92,*) nl,nwn,fill,Nc,nover
          nover=2*nover ! increase the precision
        else
          goto 401
        endif
 402    read(92,1000) linemc
        if(index(linemc(1:72),'tprime').ne.0 ) then
          read(92,*) ed,U,beta,i,run,meas,tprime
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
        
        write(6,*) ' ick  Kcx(ick)   Kcy(ick)'
        do ick=1,Nc
          write(6,6) ick,Kc(1,ick),Kc(2,ick)
        end do

c*******************************************************************************
c       Read in a few more things
c*******************************************************************************
 6      format(1x,i3,1x,f9.6,2x,f9.6)
        write(6,3) Nc
 3      format('              Nc= ',i2,' enter ick:  ',$)
        read(5,*) ick
        write(6,4)
 4      format('  enter tolerance for the root:  ',$)
        read(5,*) tolerance
        write(6,5)
 5      format('    enter damping for the root:  ',$)
        read(5,*) damping

        
        
c       Read in the local spectra A(w) and the frequency grid
        do i=1,1000
          read(9,*,end=15) w(i),Arf(i),dw(i)
        end do
 15     continue
        nw=i-1


c*******************************************************************************
c       Now spline the spectra
c*******************************************************************************
c       call spline(w,A,nw+1,0.0,0.0,scoef)

c       To evaluate at x(l)
c       call splint(w,A,scoef,nf+1,x(l),A_spline(l))    

c*******************************************************************************
c       Now integrate to find the real part of the spectra
c*******************************************************************************

        do i=2,nw-1
          Gi(i)=-pi*Arf(i)
          Gr(i)=0.0_kind
          do j=1,i-1
            Gr(i)=Gr(i) + dw(j)*Arf(j)/(w(i)-w(j))
          end do
          do j=i+1,nw
            Gr(i)=Gr(i) + dw(j)*Arf(j)/(w(i)-w(j))
          end do
        end do
        Gi(1)=-pi*Arf(1)
        Gi(nw)=-pi*Arf(nw)
c       Gr(1)=Gr(2) - (w(2)-w(1))*(Gr(3)-Gr(2))/(w(3)-w(2))
c       Gr(nw)=Gr(nw-1) + 
c     &     (w(nw)-w(nw-1))*(Gr(nw-1)-Gr(nw-2))/(w(nw-1)-w(nw-2))

        do i=1,nw
          Grf(i)=cmplx(Gr(i),Gi(i))
          sigma(i)=(0.0_kind,0.0_kind)
          write(48,*) w(i),Gr(i),Gi(i)
        end do

c       Now attempt to solve for Sigma(K,w) in
c
c                     ---
c                 Nc  \              1
c       G(K,w) = ----      ------------------------------
c                 N   /    w - ed - Epsk(K+k) - Sigma(K,w)
c                     ---
c                      k
c
c       using a Muller complex root finding routine
c       First make three guesses for Sigma
        write(67,7)
c               123456789012*123456789012*123456789012*12345678901* 
 7      format('#      w   Re[Sigma(K,w)] Im[Sigma(K,w)]  A(K,w)',
     &        '   dw  ')
        Kx=Kc(1,ick)
        Ky=Kc(2,ick)
        nguess=3
        sigma(:)=(0.0,0.0)
                  
        do i=2,nw-1      ! For each frequency search for a solution
          do iguess=1,nguess
c           For each frequency i allow nguess possible initial conditions
            if(iguess.eq.1) then
c             just use the last values of z1, z2,z3 as IC.
              z1= sigma(i-1)
              z2=z1 + delta
              z3=z1 - delta*ii
            else if(iguess.eq.2) then
c             use the coarse grained result
              z1= w(i)-ed-Epsk(Kx,Ky,tprime)-1.0_kind/Grf(i)
              z2=z1 + delta
              z3=z1 - delta*ii
            else if(iguess.eq.3) then
c             Perhaps there is a gap so that the real part has abruptly
c             switched sign.  A corresponding guess would then be:
              z1= -conjg(sigma(i-1))
              z2=z1 + delta
              z3=z1 - delta*ii
            else
c             Apparently everything has failed
              write(6,"('no solution found for i=',i3)") i
              stop   
            end if
            P1 = func(w(i),z1,ick)-Grf(i)   !i
            P2 = func(w(i),z2,ick)-Grf(i)   !i-1
            P3 = func(w(i),z3,ick)-Grf(i)   !i-2
            do iter=1,itermax
              q = (z1-z2)/(z2-z3)
              Az = q*P1 - q*(1+q)*P2 + q**2*P3
              Bz = (2*q + 1)*P1 - (1+q)**2*P2 + q**2*P3
              Cz = (1+q)*P1
              denom= Bz + sqrt(Bz**2 - 4*Az*Cz)
              if(abs(denom).lt.abs(Bz - sqrt(Bz**2 - 4*Az*Cz)))
     &          denom = Bz - sqrt(Bz**2 - 4*Az*Cz)
              znew = z1 - (z1-z2)*2*Cz/denom
              z3=z2
              z2=z1
              z1=damping*z1 +(1.0-damping)*znew
              P3=P2
              P2=P1
              P1=func(w(i),z1,ick)-Grf(i)
              if(abs(P1).lt.tolerance.and.iter.gt.2) goto 100
            end do
          end do

 100      sigma(i)=z1
        end do
        write(6,"('done first RS 0  Sigma(0)=',2g12.5)") 
     &        sigma(nw/2)
        
c*******************************************************************************
c       Now polish the roots a maximum of npolish times
c*******************************************************************************
        do ipolish=1,npolish

c         Now integrate to find the real part of Sigma again,
c         using the Kramers-Kronig relations.
c
c                 -1 /     Gi(w')
c         Gr(w) = ---| dw' ------
c                 pi /      w-w'
c
c
c                  1 /     Gr(w')
c         Gi(w) = ---| dw' ------
c                 pi /      w-w'
c

          Arf(:) = - aimag(sigma(:))/pi
          Gr(:)=0.0_kind
          do i=2,nw-1
            do j=1,i-1
              Gr(i)=Gr(i) + dw(j)*Arf(j)/(w(i)-w(j))
            end do
            do j=i+1,nw
              Gr(i)=Gr(i) + dw(j)*Arf(j)/(w(i)-w(j))
            end do
            write(44,*) w(i),Gr(i)
          end do
          close(44)
          

          Arf(:) = + Gr(:)/pi
          Gi(:)=0.0_kind
          do i=2,nw-1
            do j=1,i-1
              Gi(i)=Gi(i) + dw(j)*Arf(j)/(w(i)-w(j))
            end do
            do j=i+1,nw
              Gi(i)=Gi(i) + dw(j)*Arf(j)/(w(i)-w(j))
            end do
            write(45,*) w(i),Gi(i)
          end do
          close(45)
        
c*******************************************************************************
c         Find the "Hartree" shifts (the shift on the imaginary part
c         is needed since the range does not include the region 
c         where Re[Sigma(w)]~a/w.  Perhaps it would be better to 
c         extend the frequency range of the default model).
c*******************************************************************************
          shift=0.0
          do i=2,nw-1
            shift=shift + (Gr(i)-real(sigma(i)))
          end do
          shift=shift/real(nw-2)
        
          Gr(:)=Gr(:)-shift

          shift=0.0
          do i=2,nw-1
            shift=shift + (Gi(i)-aimag(sigma(i)))
          end do
          shift=shift/real(nw-2)
        
          Gi(:)=Gi(:)-shift
        
          sigma(:)=cmplx(Gr(:),Gi(:))
          write(6,"('done polish KK',i2,' Sigma(0)=',2g12.5)") 
     &                              ipolish,sigma(nw/2)

        
c*******************************************************************************
c         Try to find the roots yet again!
c*******************************************************************************
        
          do i=2,nw-1      ! For each frequency search for a solution
            z1= sigma(i)
            z2=z1 + delta
            z3=z1 - delta*ii
            P1 = func(w(i),z1,ick)-Grf(i)   !i
            P2 = func(w(i),z2,ick)-Grf(i)   !i-1
            P3 = func(w(i),z3,ick)-Grf(i)   !i-2
            do iter=1,itermax
              q = (z1-z2)/(z2-z3)
              Az = q*P1 - q*(1+q)*P2 + q**2*P3
              Bz = (2*q + 1)*P1 - (1+q)**2*P2 + q**2*P3
              Cz = (1+q)*P1
              denom= Bz + sqrt(Bz**2 - 4*Az*Cz)
              if(abs(denom).lt.abs(Bz - sqrt(Bz**2 - 4*Az*Cz)))
     &          denom = Bz - sqrt(Bz**2 - 4*Az*Cz)
              znew = z1 - (z1-z2)*2*Cz/denom
              z3=z2
              z2=z1
              z1=damping*z1 +(1.0-damping)*znew
              P3=P2
              P2=P1
              P1=func(w(i),z1,ick)-Grf(i)
              if(abs(P1).lt.tolerance.and.iter.gt.2) goto 101
            end do
 101        sigma(i)=z1
          end do
        write(6,"('done polish RS',i2,' Sigma(0)=',2g12.5)") 
     &                              ipolish,sigma(nw/2)
        if(abs(sigma(nw/2)-lastsig).lt.0.00001) goto 102
        lastsig=sigma(nw/2)
        end do
        
 102    if(Nc.eq.1) then
          Arf(:)=-aimag(
     &            1.0_kind/(w(:)-ed-Epsk(pi,0.0_kind,tprime)-sigma(:))
     &            )/pi
        else
          Arf(:)=-aimag(
     &            1.0_kind/(w(:)-ed-Epsk(Kx,Ky,tprime)-sigma(:))
     &            )/pi
        end if

        do i=1,nw
          write(67,8) w(i),real(sigma(i)),aimag(sigma(i)),Arf(i),dw(i)
        end do

 8      format(f11.6,2x,f11.6,2x,f11.6,2x,f11.6,2x,f11.6)

        stop
        end program dsigcluster
        
        complex(8) function func(w,z,ick)
c*******************************************************************************
        use Global
        implicit none
c*******************************************************************************
        integer ick,i,j
        real(kind) :: Kx,Ky,w,Epsk,r1,r2,epl
        complex(kind) :: csum,z
c*******************************************************************************

        r2=0.5_kind/real(nover,kind)
        csum=(0.0_kind,0.0_kind)
        do i=-nover+1,nover
        do j=-nover+1,nover
          kx=Kc(1,ick) + r2*((real(i,kind)-0.5_kind)*g(1,1) + 
     &                      (real(j,kind)-0.5_kind)*g(2,1)  )
          ky=Kc(2,ick) + r2*((real(i,kind)-0.5_kind)*g(1,2) + 
     &                      (real(j,kind)-0.5_kind)*g(2,2)  )
          epl=-0.5*(cos(kx)+cos(ky)) - 
     &           tprime*(cos(kx)*cos(ky)-1.0)
          csum=csum+
     &        1.0_kind/(w - epl - ed - z)
c         csum=csum+
c     &        1.0_kind/(w - Epsk(kx,ky,tprime) - ed - z)
        end do
        end do    
        func=csum/real((2*nover)**2,kind)
        return
        end

        real(8) function Epsk(kx,ky,tprime)
c*******************************************************************************
        implicit none
c*******************************************************************************
        real(8) kx,ky,tprime
c*******************************************************************************
        Epsk=-0.5*(cos(kx)+cos(ky)) - 
     &           tprime*(cos(kx)*cos(ky)-1.0)
        return
        end
