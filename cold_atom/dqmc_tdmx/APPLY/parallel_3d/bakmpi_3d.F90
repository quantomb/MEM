program mpi_3d
  !--------------------------------------------------------------------!
  ! This program is a MPI version of DQMC code.
  !
  ! It first performs the warmup sweep on one node and then sends out 
  ! the Hubbard-Statonivich field to others. After the measurement sweep 
  ! (performed in parallel), the root processor gathers all the measurements 
  ! and makes statisitcs.
  !
  ! Related modules
  use DQMC_CFG
  use DQMC_HUBBARD
  use DQMC_3D
  use DQMC_PHY0

  implicit none
  include 'mpif.h'

  ! QMC related arguments
  ! =====================
  !
  !integer, parameter :: OPT = STDOUT

  type(Config)  :: cfg
  type(hubbard) :: Hub
  integer       :: nx, ny, nz
  integer       :: nClass
  integer       :: i, j, k
  integer       :: ntausk
  real          :: s1, s2
  
  ! Parallel related arguments
  ! ===========================
  !
  real(wp)      :: t1, t2
  integer       :: myrank, nproc, ierr, cnt, rc
  integer       :: mat,p0_basic,p0_array,iseed
  integer, parameter:: SOURCE=0
  integer, parameter:: cfg_file=17
  integer, parameter:: out_file=18
  ! Executable

  ! Initialize MPI
  call MPI_INIT(ierr)
  if (ierr .ne. MPI_SUCCESS) then
     print *,'Error starting MPI program. Terminating.'
     call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  end if

  ! start timing
  if(myrank.eq.0) then
     t1 = MPI_WTIME()
  else
     t1 = ZERO
     t2 = ZERO
  endif

  ! Get MPI parameters
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr) !get rank
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)!get processor


  ! Open config file
  open(unit=cfg_file,file='small.in')
  ! Read configuration from standard input
  call DQMC_Read_Config(cfg, cfg_file)
  !call DQMC_Read_Config(cfg, STDIN)

  ! Initialize the Hubbard model
  call CFG_Get(cfg, "nx",nx)
  call CFG_Get(cfg, "ny",ny)
  call CFG_Get(cfg, "nz",nz)
  call CFG_Get(cfg, "tausk",ntausk)
  call CFG_Get(cfg, "seed",iseed)

  
  ! Set new data type
  cnt = Hub%n*Hub%L
  nClass = Hub%P0%nClass
 
  ! Initial 3D Hubbard
  call DQMC_Init_3D(nx,ny,nz,Hub%s) 
  !call DQMC_Init_2DPerl(nx, ny, Hub%S, IMP_TRIANGLE)  
  
  !if (myrank .eq. 0) then
     ! The root processor allocates nTask bins.
  call DQMC_Hub_Config(Hub,cfg)
  iseed=iseed+myrank
  call DQMC_Hub_sprng_init(Hub,iseed)

     ! warm up step on the first node
  do i = 1, Hub%nWarm
     ! The second parameter means no measurement should be made.
     call DQMC_Hub_Sweep(Hub, NO_MEAS0)
     call DQMC_Hub_Sweep2(Hub, Hub%nTry)
  end do
  !else
     ! The other processor allocates 1 bin.
     !call DQMC_Hub_Config(Hub,cfg)
  !end if

  ! Broadcast the HSF to other node
  !call  MPI_BCAST (Hub%HSF, 1, mat, SOURCE, MPI_COMM_WORLD, ierr)
  
  ! Change random seed for each machines
  !Hub%seed = mod(abs(Hub%seed)*(myrank+1), 4095)
  !if (mod(Hub%seed(4),2) .eq. 0) then
  !   Hub%seed(4) = Hub%seed(4) + 1
  !end if
   

  ! Measruement sweep
  ! Number of iteration = total ieration / available processors
  !nBin = nBin/nproc
  !write(*,*) 'nBin=',Hub%nBin,' nIter=',Hub%nIter 
  do i = 1,Hub%nBin
     do j = 1, Hub%nIter
        call DQMC_Hub_Sweep(Hub, Hub%nMeas)
        call DQMC_Hub_Sweep2(Hub, Hub%nTry)
     end do
     call DQMC_Phy0_Avg(Hub%P0)  

     ! Gather data
     call MPI_GATHER(Hub%P0%meas(:,1)  ,  P0_N,MPI_DOUBLE_PRECISION, Hub%P0%meas(:,1:nproc),   &
       P0_N,   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_GATHER(Hub%P0%sign(:,1)  ,  P0_N,MPI_DOUBLE_PRECISION, Hub%P0%sign(:,1:nproc),   &
       P0_N, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_GATHER(Hub%P0%G_fun(:,1) ,nClass,MPI_DOUBLE_PRECISION, Hub%P0%G_fun(:,1:nproc),  &
       nClass, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_GATHER(Hub%P0%SpinXX(:,1),nClass,MPI_DOUBLE_PRECISION, Hub%P0%SpinXX(:,1:nproc), &
       nClass, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_GATHER(Hub%P0%SpinZZ(:,1),nClass,MPI_DOUBLE_PRECISION, Hub%P0%SpinZZ(:,1:nproc), &
       nClass, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_GATHER(Hub%P0%Den0(:,1)  ,nClass,MPI_DOUBLE_PRECISION, Hub%P0%Den0(:,1:nproc),   &
       nClass, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_GATHER(Hub%P0%Den1(:,1)  ,nClass,MPI_DOUBLE_PRECISION, Hub%P0%Den1(:,1:nproc),   &
       nClass, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     write(*,*) 'Done bin',i,' in Processor', myrank
  enddo
  if (myrank .eq. 0) then
     ! Get average result
     call DQMC_Phy0_GetErr(Hub%P0)

     ! print out results
     open(out_file,file='output',access='append')
     call DQMC_Hub_Print(Hub, out_file)
     !call DQMC_Phy0_Print(Hub%P0, Hub%S, out_file)
     !call DQMC_Phy0_Print(Hub%P0, Hub%S, STDOUT)
     write(STDOUT, FMT_DBLINE) 

     ! end timing
     t2 = MPI_WTIME()
     print *, "etime=", t2-t1
  end if

  call DQMC_Hub_Free(Hub)
  call DQMC_Config_Free(cfg)

  call MPI_FINALIZE(ierr)
  
end program mpi_3d
