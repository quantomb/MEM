program mpitest
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
  use DQMC_2DPERL
  use DQMC_PHY0
  use DQMC_TDM

  implicit none
  include 'mpif.h'

  ! QMC related arguments
  ! =====================
  !
  integer, parameter :: OPT = STDOUT

  type(Config)  :: cfg
  type(hubbard) :: Hub
  type(TDM)     :: tm
  integer       :: nx, ny, L
  integer       :: i, j, k 
  integer       :: nBin, nIter, ntausk
  real          :: s1, s2
  
  ! Parallel related arguments
  ! ===========================
  !
  real(wp)      :: t1, t2
  integer       :: myrank, nproc, ierr, cnt, rc
  integer       :: mat, p0_basic, p0_array
  integer       :: n_tdm, tdm_array, iseed
  integer       :: idx_gather, rank_wrt
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
  t1 = MPI_WTIME()

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
  call CFG_Get(cfg, "tausk",ntausk)
  call CFG_Get(cfg, "nbin",nBin)
  call CFG_Get(cfg, "seed",iseed)

  call DQMC_Init_2DPerl(nx, ny, Hub%S, IMP_TRIANGLE)  
  
  if (myrank .eq. 0) then
     ! The root processor allocates nTask bins.
     call DQMC_Hub_Config(Hub,cfg)
     call DQMC_TDM_Init(TDM_ALL, tm, Hub%S, Hub%WS, Hub%B, cfg)
     open(unit=out_file,file='small.in')

     ! warm up step on the first node
     do i = 1, Hub%nWarm
        ! The second parameter means no measurement should be made.
        call DQMC_Hub_Sweep(Hub, NO_MEAS0)
        call DQMC_Hub_Sweep2(Hub, Hub%nTry)
     end do
  else
     ! The other processor allocates 1 bin.
     call DQMC_Hub_Config(Hub,cfg)
     call DQMC_TDM_Init(TDM_ALL, tm, Hub%S, Hub%WS, Hub%B, cfg)
  end if

  iseed=iseed+myrank
  call DQMC_Hub_sprng_init(Hub,iseed)

  cnt = Hub%n*Hub%L
  n_tdm = Hub%S%nClass*L

  ! Set new data typ
  call MPI_TYPE_CONTIGUOUS(cnt, MPI_INTEGER, mat, ierr)
  call MPI_TYPE_COMMIT(mat, ierr)
  call MPI_TYPE_CONTIGUOUS(P0_N, MPI_DOUBLE_PRECISION, p0_basic, ierr)
  call MPI_TYPE_COMMIT(p0_basic, ierr)
  call MPI_TYPE_CONTIGUOUS(Hub%P0%nClass, MPI_DOUBLE_PRECISION, p0_array, ierr)
  call MPI_TYPE_COMMIT(p0_array, ierr)
  call MPI_TYPE_CONTIGUOUS(n_tdm, MPI_DOUBLE_PRECISION, tdm_array,  ierr)
  call MPI_TYPE_COMMIT(tdm_array, ierr)


  ! Broadcast the HSF to other node
  call  MPI_BCAST (Hub%HSF, 1, mat, SOURCE, MPI_COMM_WORLD, ierr)

  ! Change random seed for each machines
  !Hub%seed = mod(abs(Hub%seed)*(myrank+1), 4095)
  !if (mod(Hub%seed(4),2) .eq. 0) then
  !   Hub%seed(4) = Hub%seed(4) + 1
  !end if
 
  ! Measruement sweep
  ! Number of iteration = total ieration / available processors
  nBin = nBin / nproc
  nIter= 500
  do i = 1, nBin
     do j = 1, nIter
        do k = 1, ntausk
           call DQMC_Hub_Sweep(Hub, NO_MEAS0)
           call DQMC_Hub_Sweep2(Hub, Hub%nTry)
        enddo
        call DQMC_Hub_Sweep(Hub, Hub%nMeas)
        call DQMC_TDM_Meas(tm, Hub%G_up, Hub%G_dn) 
     end do
     call DQMC_Phy0_Avg(Hub%P0)
     call DQMC_TDM_Avg(tm)

     ! Gather data
     call MPI_GATHER(Hub%P0%meas(:,1), 1, p0_basic, Hub%P0%meas(:,1:nproc), &
       1, p0_basic, source, MPI_COMM_WORLD, ierr)
     call MPI_GATHER(Hub%P0%sign(:,1), 1, p0_basic, Hub%P0%sign(:,1:nproc), &
       1, p0_basic, source, MPI_COMM_WORLD, ierr)
     call MPI_GATHER(Hub%P0%G_fun(:,1), 1, p0_array, Hub%P0%G_fun(:,1:nproc), &
       1, p0_array, source, MPI_COMM_WORLD, ierr)
     call MPI_GATHER(Hub%P0%SpinXX(:,1), 1, p0_array, Hub%P0%SpinXX(:,1:nproc), &
       1, p0_array, source, MPI_COMM_WORLD, ierr)
     call MPI_GATHER(Hub%P0%SpinZZ(:,1), 1, p0_array, Hub%P0%SpinZZ(:,1:nproc), &
       1, p0_array, source, MPI_COMM_WORLD, ierr)
     call MPI_GATHER(Hub%P0%Den0(:,1), 1, p0_array, Hub%P0%Den0(:,1:nproc), &
       1, p0_array, source, MPI_COMM_WORLD, ierr)
     call MPI_GATHER(Hub%P0%Den1(:,1), 1, p0_array, Hub%P0%Den1(:,1:nproc), &
       1, p0_array, source, MPI_COMM_WORLD, ierr)
     ! Gather data for TDM1
     idx_gather = tm%idx - 1
     rank_wrt   = (i-1)*nproc + myrank
     call MPI_GATHER(tm%T1%gnl(:,:,idx_gather)   , 1, tdm_array,  &
     tm%T1%gnl(:,:,rank_wrt),1, tdm_array, SOURCE, MPI_COMM_WORLD, ierr)
     call MPI_GATHER(tm%T1%chinl(:,:,idx_gather),  1, tdm_array,  &
     tm%T1%chinl(:,:,rank_wrt),1,tdm_array,SOURCE, MPI_COMM_WORLD, ierr)
  enddo

  if (myrank .eq. 0) then
     ! Get average result
     call DQMC_Phy0_GetErr(Hub%P0)

     ! print out results
     call DQMC_Hub_OutputParam(Hub, out_file)
     !call DQMC_Phy0_Print(Hub%P0, Hub%S, STDOUT)
     call DQMC_Phy0_Print(Hub%P0, Hub%S, out_file)
     !call DQMC_Hub_Print(Hub,out_file)
     write(STDOUT, FMT_DBLINE) 
     ! TDM dependent results
     ! Print out result
     call DQMC_TDM_Postprocessing(tm, out_file)

     ! end timing
     t2 = MPI_WTIME()
     print *, "etime=", t2-t1
  end if

  call DQMC_TDM_Free(tm)
  call DQMC_Hub_Free(Hub)
  call DQMC_Config_Free(cfg)

  call MPI_TYPE_FREE(mat, ierr)
  call MPI_TYPE_FREE(p0_basic, ierr)
  call MPI_TYPE_FREE(p0_array, ierr)
  call MPI_TYPE_FREE(tdm_array, ierr)

  call MPI_FINALIZE(ierr)
  
end program mpitest
