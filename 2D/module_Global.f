        module Global
        implicit none
        save

c       Global declarations for the program analyze_cluster
        integer :: nwn,nl,Nc,Ncw,nover
        integer :: N_sp
        integer, parameter :: nlmax=200, NCmax=26, nwnmax=800
        integer, parameter :: kind=8
        integer, parameter :: ndim=2
                
c       some global parameters
        real(kind) :: ed,U,beta,temp,tprime
        real(kind), parameter :: pi=3.141592653589793_kind
        real(kind), parameter :: eps=0.000001_kind
        complex(kind), parameter :: ii=(0.0_kind,1.0_kind)
        
        
c       the following are set in tables
        real(kind) :: Kcx(NCmax),Kcy(NCmax),wn(-nwnmax:nwnmax),
     &       g1x,g1y,g2x,g2y,origx,origy,fformfac(nwnmax),
     &       Epsbar(NCmax)
        integer :: ngroup,a1x,a1y,a2x,a2y,Rcx(NCmax),Rcy(NCmax),
     &          icrequ(NCmax,8),icrdiff(NCmax,NCmax),ic00,icpipi,
     &          ickdiff(NCmax,NCmax),ickplus(NCmax,NCmax),
     &          ickequ(NCmax,8),ict(-Ncmax:NCmax,-Ncmax:NCmax),
     &          ickdeg(NCmax),ickmap(NCmax),ick2map(NCmax,NCmax)


        end module Global
