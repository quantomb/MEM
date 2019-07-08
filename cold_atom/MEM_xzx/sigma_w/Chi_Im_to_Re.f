        program Chi_Im_to_Re
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
        integer, parameter :: itermax=2000, npolish=10, nfmax=800
        integer :: i,j,ick,run,meas,iter,nw,ised,nguess,iguess,
     &             ipolish,ios
     
c       Reals        
        real(kind), parameter :: delta=0.01_kind
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
        write(6,"('   enter the QMC spectra filename:  ',$)") 
        read(5,"(a72)") linemc
        open(unit=9,file=linemc,status='old')

c       Open a datafile for output spectra
        write(6,"('enter the output spectra filename:  ',$)") 
        read(5,"(a72)") linemc
        open(unit=92,file=linemc,status='new')
c 1000   format(a72)
 
 
c*******************************************************************************
c       Read in the old spectra
c*******************************************************************************
c       In this case, 
c
c       A(w) = -1/pi chi''(w) /(w chi(T)) 
c
c       So that 
c       oo
c       /
c       | dw (A(w) = 1/2           A(-w) = A(w)
c       /
c       0
        
c       Read in the local spectra A(w) and the frequency grid
        do i=1,1000
          read(9,*,end=15) w(i),Arf(i),dw(i)
        end do
 15     continue
        nw=i-1



c*******************************************************************************
c       Now integrate to find the real part of the Green function
c*******************************************************************************
c                oo               oo
c                /     w*A(w)     /     w*A(w)  
c       Gr(w') = | dw -------  -  | dw --------  since A(-w) = A(w)
c                /     w'-w       /      w'+w
c                0                0

        do i=1,nw-1
          Gi(i)=-pi*w(i)*Arf(i)
          Gr(i)=0.0_kind
          do j=1,i-1
            Gr(i)=Gr(i) + dw(j)*w(j)*Arf(j)/(w(i)-w(j))
          end do
          do j=i+1,nw
            Gr(i)=Gr(i) + dw(j)*w(j)*Arf(j)/(w(i)-w(j))
          end do
          do j=1,nw
            Gr(i)=Gr(i) - dw(j)*w(j)*Arf(j)/(w(i)+w(j))
          end do
        end do
        Gi(nw)=-pi*Arf(nw)

        do i=1,nw
          write(92,*) w(i),Gr(i),Gi(i)
        end do

        stop
        end program Chi_Im_to_Re
