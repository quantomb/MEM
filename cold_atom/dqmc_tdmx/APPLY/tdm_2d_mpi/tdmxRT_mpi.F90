program mpitest
  !--------------------------------------------------------------------!
  ! This program is a MPI version of DQMC code in the level
  ! that each processor run on its own Markov chain.
  !
  ! It performs the warmup sweep on one node and then do measurment. 
  ! After the measurement sweep, the root processor gathers all the 
  ! measurements and makes statisitcs.
  !--------------------------------------------------------------------!
  ! Related modules
  use DQMC_CFG
  use DQMC_HUBBARD
  use DQMC_2Dperl
  use DQMC_TDMX_RT
  use DQMC_PHY0

  implicit none
  include 'mpif.h'

  ! QMC related arguments
  ! =====================
  !
  integer, parameter :: OPT = STDOUT

  type(Config)  :: cfg
  type(hubbard) :: Hub
  type(TDMX_RT) :: tm
  integer       :: nx, ny, n, L, nClass, nl
  integer       :: i, j, k, il, ic, unit_t
  integer       :: ntausk,nphy
  real          :: s1, s2, dtau
  real(wp)      :: alphax,xdd_T,xkk_T
  real(wp),allocatable :: meas_buf(:)
  real(wp),allocatable :: sign_buf(:)
  real(wp),allocatable :: G_fun_buf(:)
  real(wp),allocatable :: SpinXX_buf(:)
  real(wp),allocatable :: SpinZZ_buf(:)
  real(wp),allocatable :: Den0_buf(:)
  real(wp),allocatable :: Den1_buf(:)
  real(8),allocatable :: meas_sum(:)
  real(8),allocatable :: sign_sum(:)
  real(8),allocatable :: G_fun_sum(:)
  real(8),allocatable :: SpinXX_sum(:)
  real(8),allocatable :: SpinZZ_sum(:)
  real(8),allocatable :: Den0_sum(:)
  real(8),allocatable :: Den1_sum(:)
  real(8),allocatable :: collect_kkt(:,:)!<KK>
  real(8),allocatable :: collect_kdt(:,:)!<KD+DK>
  real(8),allocatable :: collect_ddt(:,:)!<DD>
  real(8),allocatable :: collect_eet(:,:)!<EE>
  real(8),allocatable :: collect_x0t(:)!<X(T)>
  real(8),allocatable :: collect_vcm(:)!<K> and <D>
  real(8),allocatable :: collect_sgn(:)!sign
  ! alphax = (delta U)/(delta J)
  
  ! Parallel related arguments
  ! ===========================
  !
  real(wp)      :: t1, t2
  integer       :: myrank, nproc, ierr, cnt, rc
  integer       :: iseed
  integer       :: rank,itag,req(14)
  integer       :: istat(MPI_STATUS_SIZE)
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
  !call DQMC_Init_2DPerl(nx, ny, Hub%S, IMP_TRIANGLE)
  call DQMC_Init_2DPerl(nx, ny, Hub%S, IMP_RECTANGLE)
  
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
  call DQMC_TDMX_Init(tm, Hub%S, Hub%WS, Hub%B, cfg, alphax)
  
  n=nx*ny
  cnt = n*L
  nClass= tm%nClass
  nl = L*nClass
  nphy = P0_N
  allocate(meas_buf(nphy))
  allocate(sign_buf(nphy))
  allocate(G_fun_buf(nClass))
  allocate(SpinXX_buf(nClass))
  allocate(SpinZZ_buf(nClass))
  allocate(Den0_buf(nClass))
  allocate(Den1_buf(nClass))
  allocate(collect_kkt(nClass,L))
  allocate(collect_kdt(nClass,L))
  allocate(collect_ddt(nClass,L))
  allocate(collect_eet(nClass,L))
  allocate(collect_x0t(4))
  allocate(collect_vcm(2))
  allocate(collect_sgn(1))
  meas_buf   = 0.d0
  sign_buf   = 0.d0
  G_fun_buf  = 0.d0
  SpinXX_buf = 0.d0
  SpinZZ_buf = 0.d0
  Den0_buf   = 0.d0
  Den1_buf   = 0.d0
  collect_kkt= 0.d0
  collect_kdt= 0.d0
  collect_ddt= 0.d0
  collect_eet= 0.d0
  collect_x0t= 0.d0
  collect_vcm= 0.d0
  collect_sgn= 0.d0
  xdd_T = 0.d0
  xkk_T = 0.d0

  ! Change random seed for each machines
  !Hub%seed = mod(abs(Hub%seed)*(myrank+1), 4095)
  !if (mod(Hub%seed(4),2) .eq. 0) then
  !   Hub%seed(4) = Hub%seed(4) + 1
  !end if

  ! write header of output data file
  dtau=tm%dtau
  if(myrank.eq.0) then

     allocate(meas_sum(nphy))
     allocate(sign_sum(nphy))
     allocate(G_fun_sum(nClass))
     allocate(SpinXX_sum(nClass))
     allocate(SpinZZ_sum(nClass))
     allocate(Den0_sum(nClass))
     allocate(Den1_sum(nClass))

     open(unit=out_file,file='output',access='append')

     write(10000,*) 'Results of <K(T)K(0)>'
     write(10000,*) L,nClass,dtau*L,0.d0
     write(20000,*) 'Results of <K(T)D(0)>'
     write(20000,*) L,nClass,dtau*L,0.d0
     write(30000,*) 'Results of <D(T)D(0)>'
     write(30000,*) L,nClass,dtau*L,0.d0
     write(40000,*) 'Results of <E(T)E(0)>'
     write(40000,*) L,nClass,dtau*L,0.d0
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
        call DQMC_TDMX_Meas(tm, Hub%G_up, Hub%G_dn)
     enddo

     !Equal time measurements
     call DQMC_Phy0_Avg(Hub%P0)
     !TDMX dependent measurements
     call DQMC_TDMX_Avg(tm)
     !!!call DQMC_TDMX_WriteBin(tm,i,myrank)

     !-----------------------------------------------------------------
     !For MPI data transfer
     !set tag
     itag=1
     ! Buffer data
     ! PHY0 data
     meas_buf(:)   = Hub%P0%meas(:,i)
     sign_buf(:)   = Hub%P0%sign(:,i)
     G_fun_buf(:)  = Hub%P0%G_fun(:,i)
     SpinXX_buf(:) = Hub%P0%SpinXX(:,i)
     SpinZZ_buf(:) = Hub%P0%SpinZZ(:,i)
     Den0_buf(:)   = Hub%P0%Den0(:,i)
     Den1_buf(:)   = Hub%P0%Den1(:,i)
     ! TDMX data
     collect_kkt(:,:) = tm%dynmod_tkk(:,:,tm%avg)
     collect_kdt(:,:) = tm%dynmod_tkd(:,:,tm%avg)
     collect_ddt(:,:) = tm%dynmod_tdd(:,:,tm%avg)
     collect_eet(:,:) = tm%dynmod_t(:,:,tm%avg)
     collect_x0t(:) = tm%dynmod_w0(:,tm%avg)
     collect_vcm(:) = tm%dynmod_e(:,tm%avg)
     collect_sgn(1) = tm%sgn

     tm%dynmod_tkk = ZERO
     tm%dynmod_tkd = ZERO
     tm%dynmod_tdd = ZERO
     tm%dynmod_t = ZERO

     tm%dynmod_w0 = ZERO
     tm%dynmod_e = ZERO
     tm%sgn = ZERO


     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

     !-----------------------------------------------------------------
     ! Transfer data to rank 0, accumulate, write bins
     !
     if(myrank.eq.0) then

       meas_sum   = 0.d0
       sign_sum   = 0.d0
       G_fun_sum  = 0.d0
       SpinXX_sum = 0.d0
       SpinZZ_sum = 0.d0
       Den0_sum   = 0.d0
       Den1_sum   = 0.d0

       do rank = 0, nproc-1

          ! For Phy0 data
          if(rank.ne.0) call MPI_RECV(meas_buf,  nphy,  MPI_REAL8,rank,itag,MPI_COMM_WORLD,istat,ierr)
          if(rank.ne.0) call MPI_RECV(sign_buf,  nphy,  MPI_REAL8,rank,itag,MPI_COMM_WORLD,istat,ierr)
          if(rank.ne.0) call MPI_RECV(G_fun_buf, nClass,MPI_REAL8,rank,itag,MPI_COMM_WORLD,istat,ierr)
          if(rank.ne.0) call MPI_RECV(SpinXX_buf,nClass,MPI_REAL8,rank,itag,MPI_COMM_WORLD,istat,ierr)
          if(rank.ne.0) call MPI_RECV(SpinZZ_buf,nClass,MPI_REAL8,rank,itag,MPI_COMM_WORLD,istat,ierr)
          if(rank.ne.0) call MPI_RECV(Den0_buf,  nClass,MPI_REAL8,rank,itag,MPI_COMM_WORLD,istat,ierr)
          if(rank.ne.0) call MPI_RECV(Den1_buf,  nClass,MPI_REAL8,rank,itag,MPI_COMM_WORLD,istat,ierr)
          meas_sum   = meas_sum   + meas_buf
          sign_sum   = sign_sum   + sign_buf
          G_fun_sum  = G_fun_sum  + G_fun_buf
          SpinXX_sum = SpinXX_sum + SpinXX_buf
          SpinZZ_sum = SpinZZ_sum + SpinZZ_buf
          Den0_sum   = Den0_sum   + Den0_buf
          Den1_sum   = Den1_sum   + Den1_buf
          ! TDMX data and write bins
          if(rank.ne.0) call MPI_RECV(collect_kkt,nl,MPI_REAL8,rank,itag,MPI_COMM_WORLD,istat,ierr)
          if(rank.ne.0) call MPI_RECV(collect_kdt,nl,MPI_REAL8,rank,itag,MPI_COMM_WORLD,istat,ierr)
          if(rank.ne.0) call MPI_RECV(collect_ddt,nl,MPI_REAL8,rank,itag,MPI_COMM_WORLD,istat,ierr)
          if(rank.ne.0) call MPI_RECV(collect_eet,nl,MPI_REAL8,rank,itag,MPI_COMM_WORLD,istat,ierr)
          if(rank.ne.0) call MPI_RECV(collect_x0t,4,MPI_REAL8,rank,itag,MPI_COMM_WORLD,istat,ierr)
          if(rank.ne.0) call MPI_RECV(collect_vcm,2,MPI_REAL8,rank,itag,MPI_COMM_WORLD,istat,ierr)
          if(rank.ne.0) call MPI_RECV(collect_sgn,1,MPI_REAL8,rank,itag,MPI_COMM_WORLD,istat,ierr)

          unit_t = 10000
          write(unit_t,"('#No. bin=',i5,' from processor',i5)") i,rank
          write(unit_t,"(6(1X,f12.8))") ((collect_kkt(ic,j),j=1,L),ic=1,nClass)&
           ,collect_sgn(1),collect_vcm(1)

          unit_t = 20000
          write(unit_t,"('#No. bin=',i5,' from processor',i5)") i,rank
          write(unit_t,"(6(1X,f12.8))") ((collect_kdt(ic,j),j=1,L),ic=1,nClass)&
           ,collect_sgn(1),DSQRT(collect_vcm(1)*collect_vcm(2)*2.d0)


          unit_t = 30000
          write(unit_t,"('#No. bin=',i5,' from processor',i5)") i,rank
          write(unit_t,"(6(1X,f12.8))") ((collect_ddt(ic,j),j=1,L),ic=1,nClass)&
           ,collect_sgn(1),collect_vcm(2)

          unit_t = 40000
          write(unit_t,"('#No. bin=',i5,' from processor',i5)") i,rank
          write(unit_t,"(6(1X,f12.8))") ((collect_eet(ic,j),j=1,L),ic=1,nClass)&
           ,collect_sgn(1),-collect_vcm(1)+collect_vcm(2)*alphax

          !accumulate X(T)
          xdd_T = xdd_T + collect_x0t(3)
          xkk_T = xkk_T + collect_x0t(1)

          ! Clear buffer
          meas_buf   = 0.d0
          sign_buf   = 0.d0
          G_fun_buf  = 0.d0
          SpinXX_buf = 0.d0
          SpinZZ_buf = 0.d0
          Den0_buf   = 0.d0
          Den1_buf   = 0.d0
          collect_kkt= 0.d0
          collect_kdt= 0.d0
          collect_ddt= 0.d0
          collect_eet= 0.d0
          collect_x0t= 0.d0
          collect_vcm= 0.d0
          collect_sgn= 0.d0
       enddo
       Hub%P0%meas(:,i)   = meas_sum(:)  /DBLE(nproc)
       Hub%P0%sign(:,i)   = sign_sum(:)  /DBLE(nproc)
       Hub%P0%G_fun(:,i)  = G_fun_sum(:) /DBLE(nproc)
       Hub%P0%SpinXX(:,i) = SpinXX_sum(:)/DBLE(nproc)
       Hub%P0%SpinZZ(:,i) = SpinZZ_sum(:)/DBLE(nproc)
       Hub%P0%Den0(:,i)   = Den0_sum(:)  /DBLE(nproc)
       Hub%P0%Den1(:,i)   = Den1_sum(:)  /DBLE(nproc)
     else !myrank =x=0 nonblocking send

       call MPI_ISEND(meas_buf,  nphy,  MPI_REAL8,0,itag,MPI_COMM_WORLD,req(1),ierr)
       call MPI_ISEND(sign_buf,  nphy,  MPI_REAL8,0,itag,MPI_COMM_WORLD,req(2),ierr)
       call MPI_ISEND(G_fun_buf, nClass,MPI_REAL8,0,itag,MPI_COMM_WORLD,req(3),ierr)
       call MPI_ISEND(SpinXX_buf,nClass,MPI_REAL8,0,itag,MPI_COMM_WORLD,req(4),ierr)
       call MPI_ISEND(SpinZZ_buf,nClass,MPI_REAL8,0,itag,MPI_COMM_WORLD,req(5),ierr)
       call MPI_ISEND(Den0_buf,  nClass,MPI_REAL8,0,itag,MPI_COMM_WORLD,req(6),ierr)
       call MPI_ISEND(Den1_buf,  nClass,MPI_REAL8,0,itag,MPI_COMM_WORLD,req(7),ierr)

       call MPI_ISEND(collect_kkt,nl,   MPI_REAL8,0,itag,MPI_COMM_WORLD,req(8),ierr)
       call MPI_ISEND(collect_kdt,nl,   MPI_REAL8,0,itag,MPI_COMM_WORLD,req(9),ierr)
       call MPI_ISEND(collect_ddt,nl,   MPI_REAL8,0,itag,MPI_COMM_WORLD,req(10),ierr)
       call MPI_ISEND(collect_eet,nl,   MPI_REAL8,0,itag,MPI_COMM_WORLD,req(11),ierr)
       call MPI_ISEND(collect_x0t,4,    MPI_REAL8,0,itag,MPI_COMM_WORLD,req(12),ierr)
       call MPI_ISEND(collect_vcm,2,    MPI_REAL8,0,itag,MPI_COMM_WORLD,req(13),ierr)
       call MPI_ISEND(collect_sgn,1,    MPI_REAL8,0,itag,MPI_COMM_WORLD,req(14),ierr)

       call MPI_WAITALL(14,req,istat,ierr)
     endif

     !timing one bin
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
     deallocate(meas_sum)
     deallocate(sign_sum)
     deallocate(G_fun_sum)
     deallocate(SpinXX_sum)
     deallocate(SpinZZ_sum)
     deallocate(Den0_sum)
     deallocate(Den1_sum)

     xdd_T = xdd_T/DBLE(Hub%nBin*nproc)
     xkk_T = xkk_T/DBLE(Hub%nBin*nproc)

     write(6,*) "xdd_T:",xdd_T,'xkk_T:',xkk_T
     write(out_file,*) "xdd_T:",xdd_T
     write(out_file,*) "xkk_T:",xkk_T
  end if

  ! end timing
  t2 = MPI_WTIME()
  write(6,*) "Total time=", t2-t1
  
  deallocate(meas_buf)
  deallocate(sign_buf)
  deallocate(G_fun_buf)
  deallocate(SpinXX_buf)
  deallocate(SpinZZ_buf)
  deallocate(Den0_buf)
  deallocate(Den1_buf)
  deallocate(collect_kkt)
  deallocate(collect_kdt)
  deallocate(collect_ddt)
  deallocate(collect_eet)
  deallocate(collect_x0t)
  deallocate(collect_vcm)
  deallocate(collect_sgn)
  call DQMC_TDMX_Free(tm)
  call DQMC_Hub_Free(Hub)
  call DQMC_Config_Free(cfg)
  !write(*,*) 'Done, free arrays!'

  call MPI_FINALIZE(ierr)
  
end program mpitest
