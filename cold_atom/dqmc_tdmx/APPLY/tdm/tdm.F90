program test_tdm

  use DQMC_CFG
  use DQMC_HUBBARD
  use DQMC_2DPERL
  use DQMC_TDM
  use DQMC_PHY0
  implicit none
  
  ! This program verifies the correctness of time dependent measurement.
  integer, parameter :: OPT = STDOUT

  type(Config)  :: cfg
  type(hubbard) :: Hub
  type(TDM)     :: tm
  integer       :: nx, ny, i, j, k, nBin, nIter, ntausk 
  real          :: s1, s2
  integer,parameter::out_file=18
  
  ! Execution
  call cpu_time(s1)

  call DQMC_Read_Config(cfg, STDIN)

  ! Initialize the Hubbard model
  call CFG_Get(cfg, "tausk", ntausk)
  call CFG_Get(cfg, "nx", nx)
  call CFG_Get(cfg, "ny", ny)

  call DQMC_Init_2DPerl(nx, ny, Hub%S, IMP_RECTANGLE)
  call DQMC_Hub_Config (Hub, cfg)
  call DQMC_TDM_Init   (TDM_ALL, tm, Hub%S, Hub%WS, Hub%B, cfg)

  ! Warmup sweep
  do i = 1, Hub%nWarm
     ! The second parameter means no measurement should be made.
     call DQMC_Hub_Sweep(Hub, NO_MEAS0)
  end do
 
  ! We divide all the measurement into nBin,
  ! each having nPass/nBin pass.
  nBin   = Hub%P0%nBin 
  nIter  = Hub%nPass/nBin/ntausk
  do i = 1, nBin
     do j = 1, nIter
        do k = 1, ntausk
           call DQMC_Hub_Sweep(Hub, NO_MEAS0)
        end do
        call DQMC_TDM_Meas(tm, Hub%G_up, Hub%G_dn)
     end do
     ! Accumulate results for each bin
     call DQMC_TDM_Avg(tm) 
  end do

  call DQMC_Phy0_GetErr(Hub%P0)
  ! print out
  !call DQMC_Phy0_Print(Hub%P0, Hub%S, out_file)
  open(unit=out_file,file='output',access='append')
  call DQMC_Hub_Print(Hub,out_file)

  ! Get average results
  call DQMC_Hub_OutputParam(Hub, OPT)
  
  ! Print out results
  call DQMC_TDM_Postprocessing(tm, OPT)
  
  call cpu_time(s2)
  write(STDOUT,*) "Running time:",  s2-s1, "(second)"

  call DQMC_Hub_Free(Hub)
  call DQMC_TDM_Free(tm)
  call DQMC_Config_Free(cfg)

end program test_tdm
