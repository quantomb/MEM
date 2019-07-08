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
  use DQMC_2Dperl
!  use DQMC_TDMX
  use DQMC_PHY0

  implicit none
  include 'mpif.h'

  ! QMC related arguments
  ! =====================
  !
  integer, parameter :: OPT = STDOUT

  type(Config)  :: cfg
  type(hubbard) :: Hub
!  type(TDMX)    :: tm
  integer       :: nx, ny, n, L
  integer       :: i, j, k, il, unit_t
  integer       :: ntausk,nClass
  real          :: s1, s2, dtau
  real(wp)      :: alphax
!  real(8),allocatable :: collect_kkt(:,:)!<KK>
!  real(8),allocatable :: collect_kdt(:,:)!<KD+DK>
!  real(8),allocatable :: collect_ddt(:,:)!<DD>
!  real(8),allocatable :: collect_eet(:,:)!<EE>
!  real(8),allocatable :: collect_vcm(:,:)!<K> and <D>
!  real(8),allocatable :: collect_sgn(:)!sign
  ! alphax = (delta U)/(delta J)
  
  ! Parallel related arguments
  ! ===========================
  !
  real(wp)      :: t1, t2
  integer       :: myrank, nproc, ierr, cnt, rc
  integer       :: mat, tdm_array, n_tdm
  integer       :: iseed
  integer       :: rank
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
  !t1 = MPI_WTIME()

  ! Get MPI parameters
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr) !get rank
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc,  ierr) !get processor


  ! Open config file
  open(unit=cfg_file,file='small.in')
  ! Read configuration from config file
  call DQMC_Read_Config(cfg, cfg_file)
  !call DQMC_Read_Config(cfg, STDIN)

  call CFG_Get(cfg, "nx",nx)
  call CFG_Get(cfg, "ny",ny)
  call CFG_Get(cfg, "L", L )
  call CFG_Get(cfg, "tausk",ntausk)
  call CFG_Get(cfg, "alphax",alphax)
  call CFG_Get(cfg, "seed",iseed)


  ! Initialize the 2D or 3D rectangle lattice
  !call DQMC_Init_3D(nx,ny,nz,Hub%S)
  call DQMC_Init_2DPerl(nx, ny, Hub%S, IMP_TRIANGLE)
  
  !if (myrank .eq. 0) then
     ! Initialize the Hubbard model
  
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
     ! Initialize the Hubbard model
  !   call DQMC_Hub_Config(Hub,cfg)
  !end if
  
  ! Initialize TDMX for each processor
!  call DQMC_TDMX_Init(tm, Hub%S, Hub%WS, Hub%B, cfg, alphax)
  
  n=nx*ny
  cnt = n*L
  n_tdm = Hub%S%nClass*L
  nClass= Hub%P0%nClass
  ! Change random seed for each machines
  !Hub%seed = mod(abs(Hub%seed)*(myrank+1), 4095)
  !if (mod(Hub%seed(4),2) .eq. 0) then
  !   Hub%seed(4) = Hub%seed(4) + 1
  !end if

  ! write header of output data file
!  dtau=tm%dtau
  if(myrank.eq.0) then
!     allocate(collect_kkt(L,nproc))
!     allocate(collect_kdt(L,nproc))
!     allocate(collect_ddt(L,nproc))
!     allocate(collect_eet(L,nproc))
!     allocate(collect_vcm(2,nproc))
!     allocate(collect_sgn(nproc))
     open(unit=out_file,file='output',access='append')
!     write(10000,*) 'Results of <K(T)K(0)>'
!     write(10000,*) L,L,dtau*L,0.d0
!     write(20000,*) 'Results of <K(T)D(0)>'
!     write(20000,*) L,L,dtau*L,0.d0
!     write(30000,*) 'Results of <D(T)K(0)>'
!     write(30000,*) L,L,dtau*L,0.d0
!     write(40000,*) 'Results of <D(T)D(0)>'
!     write(40000,*) L,L,dtau*L,0.d0
!     write(5000,*) 'Results of <KK>(w_n)'
!     write(5000,*) L,L,dtau*L,0.d0
!     write(6000,*) 'Results of <KD>(w_n)'
!     write(6000,*) L,L,dtau*L,0.d0
!     write(7000,*) 'Results of <DK>(w_n)'
!     write(7000,*) L,L,dtau*L,0.d0
!     write(8000,*) 'Results of <DD>(w_n)'
!     write(8000,*) L,L,dtau*L,0.d0
!      collect_kkt=0.d0
!      collect_kdt=0.d0
!      collect_ddt=0.d0
!      collect_eet=0.d0
!      collect_vcm=0.d0
!      collect_sgn=0.d0
  endif
 
  ! Broadcast the HSF to other node
  !call MPI_BCAST(Hub%HSF, 1, mat, SOURCE, MPI_COMM_WORLD, ierr)
  
  ! Measruement sweep
  ! Number of data bin = total bin / available processors
 
  do i = 1, Hub%nBin
     ! start timing
     t1 = MPI_WTIME()
     !call DQMC_Hub_Sweep2(Hub, Hub%nTry)
     do j = 1, Hub%nIter
        do k = 1, ntausk
           call DQMC_Hub_Sweep(Hub, NO_MEAS0)
           call DQMC_Hub_Sweep2(Hub, Hub%nTry)
        enddo
        call DQMC_Hub_Sweep(Hub, Hub%nMeas)
!        call DQMC_TDMX_Meas(tm, Hub%G_up, Hub%G_dn)
     enddo

     !Equal time measurements
     call DQMC_Phy0_Avg(Hub%P0)
     !Gathering data of Hub
     call MPI_GATHER(Hub%P0%meas(:,1),  P0_N,  MPI_DOUBLE_PRECISION,Hub%P0%meas(:,1:nproc),   &
       P0_N,   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_GATHER(Hub%P0%sign(:,1),  P0_N,  MPI_DOUBLE_PRECISION,Hub%P0%sign(:,1:nproc),   &
       P0_N, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_GATHER(Hub%P0%G_fun(:,1) ,nClass,MPI_DOUBLE_PRECISION,Hub%P0%G_fun(:,1:nproc),  &
       nClass, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_GATHER(Hub%P0%SpinXX(:,1),nClass,MPI_DOUBLE_PRECISION,Hub%P0%SpinXX(:,1:nproc), &
       nClass, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_GATHER(Hub%P0%SpinZZ(:,1),nClass,MPI_DOUBLE_PRECISION,Hub%P0%SpinZZ(:,1:nproc), &
       nClass, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_GATHER(Hub%P0%Den0(:,1)  ,nClass,MPI_DOUBLE_PRECISION,Hub%P0%Den0(:,1:nproc),   &
       nClass, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_GATHER(Hub%P0%Den1(:,1)  ,nClass,MPI_DOUBLE_PRECISION,Hub%P0%Den1(:,1:nproc),   &
       nClass, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)


     !TDMX dependent measurements
!     call DQMC_TDMX_Avg(tm)
     !!!call DQMC_TDMX_WriteBin(tm,i,myrank)

     !Gathering data to tm
!     call MPI_GATHER(tm%dynmod_tkk(:,tm%avg),L,MPI_DOUBLE_PRECISION,collect_kkt(:,1:nproc),L,&
!     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!     call MPI_GATHER(tm%dynmod_tkd(:,tm%avg),L,MPI_DOUBLE_PRECISION,collect_kdt(:,1:nproc),L,&
!     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!     call MPI_GATHER(tm%dynmod_tdd(:,tm%avg),L,MPI_DOUBLE_PRECISION,collect_ddt(:,1:nproc),L,&
!     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!     call MPI_GATHER(tm%dynmod_t(:,tm%avg),  L,MPI_DOUBLE_PRECISION,collect_eet(:,1:nproc),L,&
!     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!     call MPI_GATHER(tm%dynmod_e(:,tm%avg),  2,MPI_DOUBLE_PRECISION,collect_vcm(:,1:nproc),2,&
!     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!     call MPI_GATHER(tm%sgn,                 1,MPI_DOUBLE_PRECISION,collect_sgn(1:nproc),  1,&
!     MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!     if(myrank.eq.0) then
!       do rank = 1, nproc
!
!          unit_t = 10000
!          write(unit_t,"('#No. bin=',i5,' from processor',i5)") i,rank-1
!          do il = 1, L
!             write(unit_t,"(1X,i4,1x,E16.8)") il, collect_kkt(il,rank)
!          enddo
!          write(unit_t,"('#sgn=',1x,E16.8)") collect_sgn(rank)
!          write(unit_t,"('#<K>=',1x,E16.8)") collect_vcm(1,rank)
!          write(unit_t,"('#<D>=',1x,E16.8)") collect_vcm(2,rank)
!
!
!          unit_t = 20000
!          write(unit_t,"('#No. bin=',i5,' from processor',i5)") i,rank-1
!          do il = 1, L
!             write(unit_t,"(1X,i4,1x,E16.8)") il, collect_kdt(il,rank)
!          enddo
!          write(unit_t,"('#sgn=',1x,E16.8)") collect_sgn(rank)
!          write(unit_t,"('#<K>=',1x,E16.8)") collect_vcm(1,rank)
!          write(unit_t,"('#<D>=',1x,E16.8)") collect_vcm(2,rank)

!          unit_t = 30000
!          write(unit_t,"('#No. bin=',i5,' from processor',i5)") i,rank-1
!          do il = 1, L
!             write(unit_t,"(1X,i4,1x,E16.8)") il, collect_ddt(il,rank)
!          enddo
!          write(unit_t,"('#sgn=',1x,E16.8)") collect_sgn(rank)
!          write(unit_t,"('#<K>=',1x,E16.8)") collect_vcm(1,rank)
!          write(unit_t,"('#<D>=',1x,E16.8)") collect_vcm(2,rank)
!
!          unit_t = 40000
!          write(unit_t,"('#No. bin=',i5,' from processor',i5)") i,rank-1
!          do il = 1, L
!             write(unit_t,"(1X,i4,1x,E16.8)") il, collect_eet(il,rank)
!          enddo
!          write(unit_t,"('#sgn=',1x,E16.8)") collect_sgn(rank)
!          write(unit_t,"('#<K>=',1x,E16.8)") collect_vcm(1,rank)
!          write(unit_t,"('#<D>=',1x,E16.8)") collect_vcm(2,rank)
!
!       enddo
!       collect_kkt=0.d0
!       collect_kdt=0.d0
!       collect_ddt=0.d0
!       collect_eet=0.d0
!       collect_vcm=0.d0
!       collect_sgn=0.d0
!     endif

     ! Gather data for TDM1
     !idx_gather = tm%idx - 1
     !rank_wrt   = (i-1)*nproc + myrank
     !call MPI_GATHER(tm%T1%gnl(:,:,idx_gather)   , 1, tdm_array,  &
     !tm%T1%gnl(:,:,rank_wrt),1, tdm_array, 0, MPI_COMM_WORLD, ierr)
     !call MPI_GATHER(tm%T1%chinl(:,:,idx_gather),  1, tdm_array,  &
     !tm%T1%chinl(:,:,rank_wrt),1,tdm_array,0, MPI_COMM_WORLD, ierr)

     ! end timing
     t2 = MPI_WTIME()
     write(6,*) 'Done bin=',i,' in processor',myrank,' cpu time=',t2-t1
     flush(6)
  enddo

  if (myrank .eq. 0) then
     ! Equal time average results
     call DQMC_Phy0_GetErr(Hub%P0)
     ! print out
     !call DQMC_Phy0_Print(Hub%P0, Hub%S, out_file)
     call DQMC_Hub_Print(Hub,out_file)

     ! TDM dependent results
     ! call DQMC_Hub_OutputParam(Hub, out_file)
     ! Print out result
     ! call DQMC_TDM_Postprocessing(tm, OPT)

!     deallocate(collect_kkt,collect_kdt,collect_ddt,collect_eet)
  end if

  ! end timing
  t2 = MPI_WTIME()
  !print *, "etime=", t2-t1
  
!  call DQMC_TDMX_Free(tm)
  call DQMC_Hub_Free(Hub)
  call DQMC_Config_Free(cfg)
  !write(*,*) 'Done, free arrays!'

  call MPI_FINALIZE(ierr)
  
end program mpitest
