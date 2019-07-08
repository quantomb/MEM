module DQMC_MPI

  implicit none
  include 'mpif.h'

  type MPI_SIMOR
     integer :: level       ! level of parallelization
     integer :: rank        ! global rank
     integer :: size        ! global size  

     integer :: aggr_rank
     integer :: meas_rank
     integer :: gfun_rank

     integer :: aggr_size
     integer :: meas_size
     integer :: gfun_size

     integer :: aggr_root
     integer :: meas_root
     integer :: gfun_root

     integer :: aggr_comm   ! communicator for aggregation
     integer :: meas_comm   ! communicator for measurements
     integer :: gfun_comm   ! communicator for Green's fun
  end type MPI_SIMOR

  integer, parameter :: CHANNEL_AGGR = 1
  integer, parameter :: CHANNEL_MEAS = 2
  integer, parameter :: CHANNEL_GFUN = 3

  
  integer, parameter :: PLEVEL_1 = 1
  integer, parameter :: PLEVEL_2 = 2
  integer, parameter :: PLEVEL_3 = 3
  integer, parameter :: PLEVEL_4 = 4

contains

  subroutine DQMC_MPI_Init(sim, level)
    !
    ! Purpose 
    ! =======
    !   This subroutine initializes MPI procedure.
    !
    ! Argument
    ! ========
    type(MPI_SIMOR)     :: sim
    integer, intent(in) :: level
    
    ! ... Local variables ...
    integer :: ierr, rc

    ! ... Executable ... 
    
    call MPI_INIT(ierr)
    if (ierr .ne. MPI_SUCCESS) then
       print *,'Error starting MPI program. Terminating.'
       call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
    end if

    sim%level = level
    
    if (level .eq. 1) then
       ! Get MPI parameters
       call MPI_COMM_RANK(MPI_COMM_WORLD, sim%rank, ierr)
       call MPI_COMM_SIZE(MPI_COMM_WORLD, sim%size, ierr)

       sim%aggr_rank  = sim%rank
       sim%meas_rank  = sim%rank
       sim%gfun_rank  = sim%rank

       sim%aggr_size  = sim%size
       sim%meas_size  = sim%size
       sim%gfun_size  = sim%size

       sim%aggr_root  = 0
       sim%meas_root  = 0
       sim%gfun_root  = 0
       
       sim%aggr_comm = MPI_COMM_WORLD
       sim%meas_comm = MPI_COMM_WORLD
       sim%gfun_comm = MPI_COMM_WORLD
    end if

  end subroutine DQMC_MPI_Init

  ! =========================================================

  subroutine DQMC_MPI_Final(sim)
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes MPI procedure.
    !
    ! Argument
    ! ========
    type(MPI_SIMOR) :: sim
    
    ! ... Local variables ...
    integer :: ierr, rc

    ! ... Executable ... 
    
    call MPI_FINALIZE(ierr)
    if (ierr .ne. MPI_SUCCESS) then
       print *,'Error in finalize MPI program. Terminating.'
       call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
    end if
  end subroutine DQMC_MPI_Final

  ! =========================================================

  function DQMC_MPI_Is_Root(sim, channel) result(isRoot)
    !
    ! Purpose 
    ! =======
    !   This subroutine finalizes MPI procedure.
    !
    ! Argument
    ! ========
    type(MPI_SIMOR)     :: sim
    integer, intent(in) :: channel
    logical             :: isRoot
    
    ! ... Executable ... 

    select case(channel)
    case(CHANNEL_AGGR)
       isRoot = (sim%aggr_rank .eq. sim%aggr_root)
    end select

  end function DQMC_MPI_Is_Root

end module DQMC_MPI
