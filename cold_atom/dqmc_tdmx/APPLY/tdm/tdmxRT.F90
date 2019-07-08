program test_tdm

  use DQMC_CFG
  use DQMC_HUBBARD
  use DQMC_2DPERL
  use DQMC_TDMX_RT
  implicit none
  
  ! This program verifies the correctness of time dependent measurement.
  integer, parameter :: OPT = STDOUT

  type(Config)  :: cfg
  type(hubbard) :: Hub
  type(TDMX_RT)    :: tm
  integer       :: nx, ny, L, nClass
  integer       :: i, j, k
  integer       :: ntausk,iseed
  real(wp)      :: s1, s2, dtau
  real(wp)      :: alphax
  integer,parameter:: infile=17
  integer,parameter::outfile=18
  
  ! Execution
  call cpu_time(s1)

  open(UNIT=infile,FILE='small.in',STATUS='OLD')
  call DQMC_Read_Config(cfg,infile)
  !call DQMC_Read_Config(cfg, STDIN)

  call CFG_Get(cfg, "tausk", ntausk)
  call CFG_Get(cfg, "nx", nx)
  call CFG_Get(cfg, "ny", ny)
  call CFG_Get(cfg, "L",  L )
  call CFG_Get(cfg, "dtau",dtau)
  call CFG_Get(cfg, "seed",iseed)
  call CFG_Get(cfg, "alphax",alphax)
  ! Initialize the 2D square lattice
  call DQMC_Init_2DPerl(nx, ny, Hub%S, IMP_RECTANGLE)
  ! Initialize the Hubbard model
  call DQMC_Hub_Config (Hub, cfg)
  call DQMC_Hub_sprng_init(Hub,iseed)

  write(*,*) 'do dynamic modulation'
  ! Warmup sweep
  do i = 1, Hub%nWarm
     ! The second parameter means no measurement should be made.
     call DQMC_Hub_Sweep(Hub, NO_MEAS0)
     ! Do nTry of global update
     !call DQMC_Hub_Sweep2(Hub, Hub%nTry)
  end do
  call DQMC_TDMX_Init(tm, Hub%S, Hub%WS, Hub%B, cfg, alphax )
  nClass=tm%nClass
  write(*,*) 'In tdm, nClass=',nClass
  write(6,*) 'Done warm!'
  write(10000,*) 'Results of <K(T)K(0)>'
  write(10000,*) L,nClass,dtau*L,0.d0
  write(20000,*) 'Results of <K(T)D(0)>'
  write(20000,*) L,nClass,dtau*L,0.d0
  write(30000,*) 'Results of <D(T)D(0)>'
  write(30000,*) L,nClass,dtau*L,0.d0
  write(40000,*) 'Results of <E(T)E(0)>'
  write(40000,*) L,nClass,dtau*L,0.d0
  ! We divide all the measurement into nBin
  do i = 1, Hub%nBin
     do j = 1, Hub%nIter
        do k = 1, ntausk
           call DQMC_Hub_Sweep(Hub, NO_MEAS0)
           !call DQMC_Hub_Sweep2(Hub, Hub%nTry)
        end do
        ! Measuring time-dependent quantities
        call DQMC_Hub_Sweep(Hub, Hub%nMeas)
        call DQMC_TDMX_Meas(tm, Hub%G_up, Hub%G_dn)
     end do
     ! Accumulate results for each bin
     call DQMC_Phy0_Avg(Hub%P0)
     call DQMC_TDMX_Avg(tm)
     call DQMC_TDMX_WriteBin(tm,i,0)
     write(6,*) 'Done bin=',i,'!'
  end do
  
  call DQMC_Phy0_GetErr(Hub%P0)

  ! print out results
  open(UNIT=outfile,FILE='ouput',STATUS='unknown')
  call DQMC_Hub_Print(Hub,outfile)
  !call DQMC_Hub_OutputParam(Hub, STDOUT)
  !call DQMC_Phy0_Print(Hub%P0, Hub%S, STDOUT)

  call cpu_time(s2)
  write(STDOUT,*) "Running time:",  s2-s1, "(second)"

  call DQMC_Hub_Free(Hub)
  call DQMC_TDMX_Free(tm)
  call DQMC_Config_Free(cfg)

end program test_tdm
