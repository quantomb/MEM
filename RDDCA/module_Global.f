        module Global
        implicit none
        save

c       Global declarations for the program analyze_cluster

c       ***********************INTEGER***************************
c       parameters
        integer, parameter :: kind=8
        integer, parameter :: ndim=2
        integer, parameter :: iprint=10       ! set print level 
        integer, parameter :: myrank=0       ! dummy variable
        integer, parameter :: lud=6          ! make so the lun for debugging info
c       variables
        integer :: nwn,nl,Nc,Ncw,nt,ntw,nover,N_sp,ic0,icpi,ngroup,
     &             a(3,3)
c       allocatable variables
        integer, allocatable :: Rc(:,:),ig(:),jg(:),icrequ(:,:),
     &                          icrdiff(:,:),ickdiff(:,:),
     &                          ickplus(:,:),ickequ(:,:),
     &                          ict(:,:),ickdeg(:),ickmap(:),
     &                          ick2map(:,:)
                

c       ************************REAL*****************************
c       parameters
        real(kind), parameter :: pi=3.1415926536,eps=0.0001_kind,
     &                           zeror=0.00_kind,halfr=0.50_kind,
     &                           oner=1.00_kind,twor=2.00_kind,
     &                           threer=3.00_kind

c       variables
        real(kind) :: ed,U,beta,temp,tprime,g(3,3)
c       allocatable variables
        real(kind), allocatable :: Kc(:,:),wn(:),fformfac(:),
     &                             Epsbar(:)

c       ***********************COMPLEX***************************
c       parameters
        complex(kind), parameter :: ii=(0.0_kind,1.0_kind)
c       variables
c       allocatable variables
        



        end module Global

c*************************************************************************

        subroutine allocate_1p

        use Global
        integer :: infot

        infot=0

c       integers
        allocate (Rc(3,Nc),stat=info)
        infot=infot+info
        allocate (icrequ(Nc,8),stat=info)
        infot=infot+info
        allocate (icrdiff(Nc,Nc),stat=info)
        infot=infot+info
        allocate (ickdiff(Nc,Nc),stat=info)
        infot=infot+info
        allocate (ickplus(Nc,Nc),stat=info)
        infot=infot+info
        allocate (ickequ(Nc,8),stat=info)
        infot=infot+info
        allocate (ict(-Nc:Nc,-Nc:Nc),stat=info)
        infot=infot+info
        allocate (ickdeg(Nc),stat=info)
        infot=infot+info
        allocate (ickmap(Nc),stat=info)
        infot=infot+info
        allocate (ick2map(Nc,Nc),stat=info)
        infot=infot+info
        allocate(ig(Nc),stat=info)
        infot=infot+info
        allocate(jg(Nc),stat=info)
        infot=infot+info

c       reals
        allocate (Kc(3,Nc),stat=info)
        infot=infot+info
        allocate (wn(-nwn:nwn),stat=info)
        infot=infot+info
        allocate (fformfac(nwn),stat=info)
        infot=infot+info
        allocate(Epsbar(Nc),stat=info)
        infot=infot+info

        if(infot.ne.0) then
          write(6,*) '********** ERROR in allocate_1p*********'
          write(6,*) 'infot=',infot
          stop
        end if

        return
        end subroutine allocate_1p

c*************************************************************************

        subroutine deallocate_1p

        use Global

c       integers
        deallocate (Rc)
        deallocate (icrequ)
        deallocate (icrdiff)
        deallocate (ickdiff)
        deallocate (ickplus)
        deallocate (ickequ)
        deallocate (ict)
        deallocate (ickdeg)
        deallocate (ickmap)
        deallocate (ick2map)
        deallocate(ig)        
        deallocate(jg)        
        
c       reals
        deallocate (Kc)
        deallocate (wn)
        deallocate (fformfac)
        deallocate (Epsbar)


        return
        end subroutine deallocate_1p

c*************************************************************************
