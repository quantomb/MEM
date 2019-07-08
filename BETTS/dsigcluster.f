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
     &                damping,shift

c       Complex
        complex(kind) :: q,Az,Bz,Cz,z1,z2,z3,znew,denom,P1,P2,P3,
     &                   sigma(nfmax),Grf(nfmax),func,lastsig
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
 
c*******************************************************************************
c       Read in some parameters
c*******************************************************************************
c       Read some parameters from the QMC sigma file
 401    read(92,"(a72)") linemc
        if(index(linemc(1:72),'Results').ne.0 ) then
         read(92,*) nl,nwn,fill,cluster,ntot,ndim
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
        write(6,"('  enter tolerance for the root:  ',$)")
        read(5,*) tolerance
        write(6,"('    enter damping for the root:  ',$)")
        read(5,*) damping

        
        
c       Read in the local spectra A(w) and the frequency grid
        do i=1,1000
          read(9,*,end=15) w(i),Arf(i),dw(i)
        end do
 15     continue
        nw=i-1



c*******************************************************************************
c       Now integrate to find the real part of the Green function
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

c       The Green's fuction energy argument  is  relative to mu, thus 
c       what we denote   G(k,w) below is in fact G(k,w-mu)=G(k,w+ed)
c
c       Now attempt to solve for Sigma(K,w) in
c
c                     ---
c                 Nc  \              1
c       G(K,w) = ----      ------------------------------
c                 N   /    w - ed-Epsk(K+k) - Sigma(K,w)
c                     ---
c                      k
c
c       
c
c       using a Muller complex root finding routine
        write(67,"('#      w   Re[Sigma(K,w)] Im[Sigma(K,w)]  A(K,w)   dw  ')")
        Kx=Kc(1,ick); Ky=Kc(2,ick); if(ndim.eq.3) Kz=Kc(3,ick)
        sigma(:)=(0.0,0.0)
                  
c       First make three guesses for Sigma
        nguess=4
        do i=2,nw-1      ! For each frequency search for a solution
          do iguess=1,nguess
c           For each frequency i allow nguess possible initial conditions
            if(iguess.eq.1) then
c             use the coarse grained result
              z1= w(i)-ed-Epsk(Kx,Ky,Kz)-1.0_kind/Grf(i)
              z2=z1 + delta
              z3=z1 - delta*ii

             else if(iguess.eq.2) then
c             Perhaps  the step in Muller route finder is too large, let delta*0.1, 
c             use the coarse grained result. A corresponding guess would then be:
              z1= w(i)-ed-Epsk(Kx,Ky,Kz)-1.0_kind/Grf(i)
              z2=z1 + delta*0.1
              z3=z1 - delta*0.1*ii
              write(6,"('2nd try, delta=0.1*delta, use coarse grained z1 for i=',i3)") i

            else if(iguess.eq.3) then
c             just use the last values of z1, z2,z3 as IC.
              z1= sigma(i-1)
              z2=z1 + delta
              z3=z1 - delta*ii
            else if(iguess.eq.4) then
c             Perhaps there is a gap so that the real part has abruptly
c             switched sign.  A corresponding guess would then be:
              z1= -conjg(sigma(i-1))
              z2=z1 + delta
              z3=z1 - delta*ii
              write(6,"('4th try, switch ReSig sign i=',i3)") i
            else
c             Apparently everything has failed
              write(6,"('no solution found for i=',i3)") i
              do j=1,nw
               write(400+i,"(e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)")
     &          w(j), real(Grf(j)), aimag(Grf(j)), real(sigma(j)), aimag(sigma(j))
              end do
c              stop   
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
        write(6,"('done first RS 0  Sigma(0)=',2g12.5)") sigma(nw/2)
        do j=1,nw
               write(700,"(e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)")
     &          w(j), real(Grf(j)), aimag(Grf(j)), real(sigma(j)), aimag(sigma(j))
              end do
        lastsig=0.0
c*******************************************************************************
c       Now polish the roots a maximum of npolish times
c*******************************************************************************
        do ipolish=1,npolish

c         Now integrate to find the real part of Sigma again, using 
c         the Kramers-Kronig relations.  Note that these relations dont 
c         include the parts of Sigma which remain finite as |w|-->oo
c
c                 -1 /     Si(w')
c         Sr(w) = ---| dw' ------
c                 pi /      w-w'
c
c
c                  1 /     Sr(w')
c         Si(w) = ---| dw' ------
c                 pi /      w-w'
c
c         Integrate Im[Sigma] to obtain Re[Sigma]
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
          

c         Integrate Re[Sigma] to obtain Im[Sigma]
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
c         Find the average shift of Re[Sigma]
          shift=0.0
          do i=2,nw-1
            shift=shift + (Gr(i)-real(sigma(i)))
          end do
          shift=shift/real(nw-2)
        
          Gr(:)=Gr(:)-shift

c         Find the average shift of Im[Sigma]
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
c         Try to find the roots yet again, starting with Sigma obtained above.
c         Note that the final Sigma(K,w) is obtained from the zeros of 
c         the coarse-graining equation.
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
c        if(abs(sigma(nw/2)-lastsig).lt.0.00001) goto 102
        if(abs(sigma(nw/2)-lastsig).lt.tolerance) goto 102
        lastsig=sigma(nw/2)
        end do  ! ipolish=1,npolish
        
 102    if(Nc.eq.1) then
          Arf(:)=-aimag(
     &            1.0_kind/(w(:)-ed-Epsk(pi,zeror,zeror)-sigma(:))
     &            )/pi
        else
          Arf(:)=-aimag(
     &            1.0_kind/(w(:)-ed-Epsk(Kx,Ky,Kz)-sigma(:))
     &            )/pi

        end if

        do i=1,nw
          write(67,"(f11.6,2x,f11.6,2x,f11.6,2x,f11.6,2x,f11.6)") 
     &      w(i),real(sigma(i)),aimag(sigma(i)),Arf(i),dw(i)
        end do

        stop
        end program dsigcluster
        
        complex(8) function func(w,z,ick)
c*******************************************************************************
        use Global
        implicit none
c*******************************************************************************
        integer ick,i,j
        real(kind) :: Kx,Ky,kz,w,Epsk,r1,r2,epl
        complex(kind) :: csum,z
c*******************************************************************************
        csum=(0.0_kind,0.0_kind)
        if(ndim.eq.2) then
        do i=1,nwsc
            kx=Kc(1,ick) + kt(1,i)
            ky=Kc(2,ick) + kt(2,i)
            epl=-0.5*(cos(kx)+cos(ky)) - tprime*(cos(kx)*cos(ky)-1.0)
     &             - 0.5*tdprime*(cos(2.d0*kx)+cos(2.d0*ky))
     &             - t3prime*(cos(2.d0*kx)*cos(ky)
     &                         +cos(kx)*cos(2.d0*ky))   


            csum=csum + 1.0_kind/(w - epl -ed- z)
        end do
        else
        do i=1,nwsc
            kx=Kc(1,ick) + kt(1,i)
            ky=Kc(2,ick) + kt(2,i)
            kz=Kc(3,ick) + kt(3,i)
            epl=-0.5_kind*(cos(kx)+cos(ky)+cos(kz))
     &         - tprime*(cos(kx)*cos(ky)
     &                  +cos(kx)*cos(kz)
     &                  +cos(ky)*cos(kz))
            csum=csum + 1.0_kind/(w - epl -ed- z)
        end do
        end if

        func=csum/real(nwsc,kind)
        return
        end



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
