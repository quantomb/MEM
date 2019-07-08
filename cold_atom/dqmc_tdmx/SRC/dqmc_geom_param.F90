module DQMC_GEOM_PARAM

  implicit none

  integer   ,   parameter  :: rdim=3
  real*8    ,   parameter  :: toll=1.d-6
  real*8    ,   parameter  :: pi=3.1415926535897932d0
  complex*16,   parameter  :: im=(0.d0,1.d0)
  
  integer,      parameter  :: N_Fields  = 12  ! number of fields
  integer,      parameter  :: NDIM_F    =  1  
  integer,      parameter  :: MU_F      =  2
  integer,      parameter  :: PRIM_F    =  3
  integer,      parameter  :: SUPER_F   =  4
  integer,      parameter  :: ORB_F     =  5
  integer,      parameter  :: HAMILT_F  =  6
  integer,      parameter  :: KPOINT_F  =  7
  integer,      parameter  :: SYMM_F    =  8
  integer,      parameter  :: PHASE_F   =  9
  integer,      parameter  :: BONDS_F   = 10
  integer,      parameter  :: PAIRS_F   = 11
  integer,      parameter  :: DILUT_F   = 12
  
  character*10, parameter  :: INPUT_FIELDS(N_Fields) = (/ &
       '#NDIM','#MU','#PRIM',    &
       '#SUPER','#ORB','#HAMILT','#K-POINT','#SYMM', &
       '#PHASE','#BONDS','#PAIR','#DILUT'/)
  
  integer                  :: inpunit
  logical                  :: Found_Field(N_fields)
  save

contains


  !----------------------------------------------
  ! Determines which Fields are specified
  !----------------------------------------------
  subroutine analyze_input
    integer            :: i,ios
    logical            :: stop_exe
    character(len=100) :: str
    
    !!!roger: can be better 
    do i=1,N_Fields
       rewind(inpunit)
       Found_Field(i)=.false.
       do
          read(inpunit,'(A)',iostat=ios)str
          if(ios/=0)exit
          if(index(str,INPUT_FIELDS(i))>0)Found_Field(i)=.true.
       enddo
    enddo

    stop_exe=.false.
    do i=NDIM_F,HAMILT_F
       if(.not.Found_Field(i))then
          write(*,'(A)')INPUT_FIELDS(i),' is compulsory in input'
          stop_exe=.true.
       endif
    enddo

    if(stop_exe)stop

    if(Found_Field(PAIRS_F).and..not.Found_Field(BONDS_F))then
       write(*,'(A)') '#PAIR requires #BONDS to be specified in input'
       stop
    endif
    ! do i=1,N_Fields
    !   write(*,*)INPUT_FIELDS(i),Found_Field(i)
    ! enddo
  end subroutine analyze_input



  !-----------------------------------------------------------------
  ! Move to record where string is found. If not found return false.
  !-----------------------------------------------------------------
  logical function move_to_record(string,iunit)
    implicit none
    integer iunit,ios,istring
    character*100 line
    character*(*) string
    rewind(iunit)
    do
       read(iunit,'(A)',iostat=ios)line
       if(ios.ne.0)then
          !write(*,*)'Problem finding string',string
          !stop
          move_to_record=.false.
          exit
       endif
       istring=index(line,string)
       if(istring==0)cycle
       move_to_record=.true.
       exit
    enddo
  end function move_to_record

  
  
  !-------------------------------
  !Determinant of a 3x3 matrix
  !-------------------------------
  real*8 function get_det(a)
    real*8, intent(in) :: a(3,3)
    real*8 :: d(3)
    d(1)=a(2,1)*a(3,2)-a(3,1)*a(2,2)
    d(2)=a(3,1)*a(1,2)-a(1,1)*a(3,2)
    d(3)=a(1,1)*a(2,2)-a(2,1)*a(1,2)
    get_det=d(1)*a(1,3)+d(2)*a(2,3)+d(3)*a(3,3)
  end function get_det
  
  


  !-----------------------------
  !Inverse of a 3x3 matrix
  !-----------------------------
  subroutine get_inverse(a,inv)
    real*8 a(3,3),inv(3,3),det
    det=get_det(a)
    inv(1,1)=a(2,2)*a(3,3)-a(2,3)*a(3,2)
    inv(1,2)=a(3,2)*a(1,3)-a(3,3)*a(1,2)
    inv(1,3)=a(1,2)*a(2,3)-a(1,3)*a(2,2)
    inv(2,1)=a(2,3)*a(3,1)-a(2,1)*a(3,3)
    inv(2,2)=a(3,3)*a(1,1)-a(3,1)*a(1,3)
    inv(2,3)=a(1,3)*a(2,1)-a(1,1)*a(2,3)
    inv(3,1)=a(2,1)*a(3,2)-a(2,2)*a(3,1)
    inv(3,2)=a(3,1)*a(1,2)-a(3,2)*a(1,1)
    inv(3,3)=a(1,1)*a(2,2)-a(1,2)*a(2,1)
    inv(:,:)=inv(:,:)/det
  end subroutine get_inverse
  
  
end module DQMC_GEOM_PARAM
