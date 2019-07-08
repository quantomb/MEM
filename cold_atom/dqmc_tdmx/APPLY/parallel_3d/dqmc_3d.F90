module DQMC_3D

  use DQMC_UTIL
  use DQMC_CFG
  use DQMC_STRUCT
  use DQMC_HUBBARD
  implicit none 

  ! 
  ! This module defines subroutines to initialize the data structure of 
  ! a three-dimensional periodic rectangular lattice (3DPerl).
  !     
  !  Subroutine List
  !  ===============
  !    DQMC_Comp_3D(IPT,OPT) : read in data for 3D rectangular 
  !                                    lattice.
  !    DQMC_Init_3D(nx,ny,S) : construct data for a 3D rectangular lattice.
  !

  integer, parameter :: IMP_TRIANGLE  = 1
  integer, parameter :: IMP_RECTANGLE = 2
  
contains
  
  !---------------------------------------------------------------------!

  subroutine DQMC_Comp_3D(IPT, OPT)
    !
    ! Purpose
    ! =======
    !    This subroutine reads in data for a 3D rectangular lattice,
    !    and the calls DQMC_Init_3D to construct the lattice structure.
    !    Basically, the only parameters needed are Nx, Ny, and Nz, which
    !    the size of the lattice along the x, y and z coordinate.
    !
    ! Arguments
    ! =========
    !
    integer, intent(in) :: IPT, OPT ! Input/output handle

    ! ... Local scalar ...
    type(config)  :: cfg
    type(Hubbard) :: Hub          ! Hubbard model
    integer       :: nx, ny, nz

    ! ... Executable ...

    call DQMC_Read_Config(cfg, IPT)

    ! Initialize the geometry
    call CFG_Get(cfg, "nx", nx) 
    call CFG_Get(cfg, "ny", ny)
    call CFG_Get(cfg, "nz", nz)
    call DQMC_Init_3D(nx, ny, nz, Hub%S)

    ! Initialize the rest data
    call DQMC_Hub_Config(Hub, cfg)

    ! Execution MC loop
    call DQMC_Hub_Run(Hub)

    ! Print computed results
    call DQMC_Hub_Print(Hub, OPT)
    
    ! Clean up the used storage
    call DQMC_Hub_Free(Hub)
    call DQMC_Config_Free(cfg)
    
  end subroutine DQMC_Comp_3D
  
  !---------------------------------------------------------------------!

  subroutine DQMC_Init_3D(nx, ny, nz, S)    
    !
    ! Purpose
    ! =======
    !    This subroutine constuctures data structure for a 3D 
    !    periodic rectangular lattice.  
    !
    ! Details
    ! =======
    !    For a Nx*Ny*Nz 3D rectangular lattice.
    !
    !    1. The sites are numbered from 1, which is the site on the 
    !       bottom south-west corner. The numbering is row major in 3D.
    !       It means that the index increases along x direction first 
    !       (west to east) and then y-direction (south to north). After 
    !       sweep through the layer, then increase in z-direction 
    !       (bottom to top).
    !
    !    2. Adjacency (T) has exact 6 elements per site: west
    !       , east, north and south (in same layer) (x-y plane), 
    !       addition with bottom and top (z-direction). Since the 
    !       adjacency is cyclic, the code has had spacial treat for 
    !       boundary sites.
    !       *** The way it computed here is not satisfied the 
    !           'checkboard' order. should change it later.
    !
    !    3. The number of unique distance is computed as follows.
    !       Let long = max(Nx, Ny, Nz) and short = min(Nx, Ny, Nz).
    !       The distinct distance sites form a trapezoid
    !       The bottom is (long/2+1), the top is ((long-short)/2+1)
    !       and the height is (short/2+1). Therefore, 
    !       
    !           nClass = (short/2+1)*(long-short/2+2)/2
    !       
    !    4. The phase (P) is computed in the rules that
    !       (a) adjacent sites have opposite phase.
    !       (b) site 1 is phased +.
    !
    ! Arguments
    ! =========
    !
    integer, intent(in)         :: nx, ny, nz  ! dim of the lattice
    type(Struct), intent(inout) :: S           ! Struct

    ! ... local vars ...
    integer  :: n, nxy, nlink    ! Order of matrix T and D 
    integer  :: i, ix, iy, iz    ! Loop iterator
    integer  :: j, jx, jy, jz    ! Loop iterator
    integer  :: up, dn, lt, rt   ! adj and neighbor
    integer  :: upt, dnt         ! adj and neighbor
    integer  :: dx, dy, dz
    integer  :: tbl(nx,ny,nz)    ! lookup table
    
    integer  :: idx              ! 
    integer  :: tmp(nx*ny*nz, nx*ny*nz)
    real(wp) :: cord(nx*ny*nz,3)
    
    ! ... parameters ...
    integer, parameter :: NADJ  = 6  ! Number of adjacencies
    real(wp), parameter  :: TWOPI     = 6.283185307179586

    ! ... Executable ...

    nxy = nx*ny
    n   = nxy*nz
    S%nSite = n

    write(S%name,'(A,I3,A,I3,A,I3,A,I5)') &
         "3D Lattice (Nx=", nx, ", Ny=", ny, ", Nz=", nz, &
         ") total sites=", S%nSite
    
    ! Compute all distinct sites
    S%nWave = 0

    ! memory allocation
    allocate(S%P(n))
    allocate(S%dim(3))
    allocate(S%map(n))

    S%dim(1) = nx
    S%dim(2) = ny
    S%dim(3) = nz
    S%Map    = 1
    S%nGroup = 1

    ! Build adjacent matrix for checkerboard method
    S%n_t = 1
    tmp = 0
    do j = 1, ny
       do i = 1, nx
          idx = (j-1)*nx+i
          
          ! north link
          up  = idx + nx
          if (up .gt. nxy) then
             up = up - nxy
          end if
          
          ! south link
          dn = idx - nx
          if (dn .le. 0) then
             dn = dn + nxy
          end if
          
          ! west link
          if (i .eq. 1) then
             lt = idx + nx - 1
          else
             lt = idx - 1
          end if
             
          ! east link
          if (i .eq. nx) then
             rt = idx - nx + 1
          else
             rt = idx + 1
          end if

          tmp(idx, up) = 1
          tmp(idx, dn) = 1
          tmp(idx, lt) = 1
          tmp(idx, rt) = 1
       end do
    end do
    
    ! copy the first layer to others
    ! 1st layer (i=0) to nz layer (i=n-nxy=(nz-1)*nxy)
    do i = 0, n-nxy, nxy
       tmp(i+1:i+nxy, i+1:i+nxy) = tmp(1:nxy, 1:nxy)
       ! links between layers
       do j = 1, nxy
          ! down link
          dnt = i+j-nxy
          if(dnt.le.0) then
          ! should go to the last layer.
          ! use 'le' of last site of 1st layer
             dnt=n-nxy+j
          endif
          tmp(i+j,dnt) = 1
          ! up link
          upt = i+j+nxy
          if(upt.gt.n) then
          ! should go to the first layer
          ! use 'gt' for last site of (nz-1) layer
             upt=j
          endif
          tmp(i+j,upt) = 1
          !tmp(i+j-nxy,i+j) = 1
       end do
    end do

    ! link the last layer and the first layer
    !if (nz .gt. 2) then
    !   do j = 1, nxy
    !      tmp(j,n-nxy+j) = 1
    !      tmp(n-nxy+j,j) = 1
    !   end do
    !end if

    nlink = n*NADJ !number of whole adjacencies
    if (nz .le. 2) then
       nlink = nlink - 2*nxy
    end if
    
    !lt = 0
    !do i = 1, n
    !   do j= 1,n
    !      lt = lt+tmp(j,i)
    !   enddo
    !enddo
    !write(*,*) 'Nonzero elements:', lt
    !pause


    call DQMC_CCS_Compress(n, nlink, tmp, S%T)

    ! build up the distance matrix.
    S%nClass = (nx/2+1)*(ny/2+1)*(nz/2+1)
    allocate(S%D(n,n))
    allocate(S%F(S%nClass))
    allocate(S%clabel(S%nClass))
    S%D  = 0
    S%F  = 0

    S%F(1) = n
    write(S%clabel(1),"(I3,I3,I3)") 0, 0, 0
    tbl = 0
    idx = 1  ! the first class is a special one (dx=dy=dz=0)
    cord(1,1) = ZERO
    cord(1,2) = ZERO
    cord(1,3) = ZERO

    do i = 1, n
       !! compute the index of i
       ix = MOD(i-1, nx)+1
       iy = MOD(i-1, nxy)/nx + 1
       iz = INT((i-1)/nxy) + 1

       !! initial the index of j
       do j = i+1, n
          !! compute the index of j
          jx = MOD(j-1, nx)+1
          jy = MOD(j-1, nxy)/nx + 1
          jz = INT((j-1)/nxy) + 1

          !! compute the distance
          dx = ABS(ix-jx)
          dx = MIN(dx,nx-dx)+1
          dy = ABS(iy-jy)
          dy = MIN(dy,ny-dy)+1
          dz = ABS(iz-jz)
          dz = MIN(dz,nz-dz)+1

          ! not found
          if (tbl(dx,dy,dz) .eq. 0) then

             ! Creat a new node
             idx = idx + 1
             tbl(dx, dy, dz) = idx
             S%D(i,j) = idx
             write(S%clabel(idx),"(I3,I3,I3)") dx-1, dy-1, dz-1

             ! Build a new row of COS table
             cord(idx, 1) = dx - 1
             cord(idx, 2) = dy - 1
             cord(idx, 3) = dz - 1
          else ! found

             S%D(i,j) = tbl(dx, dy, dz)

          end if

          ! matrix D is symmetric
          S%D(j,i) = S%D(i,j)

          ! increase count by 2
          S%F(S%D(i,j)) = S%F(S%D(i,j)) + 2

       end do

       ! site i to i
       S%D(i,i) = 1

    end do

    ! Initialize phase matrix
    S%P(1:nx:2) = -1
    S%P(2:nx:2) = 1
    do i = 2, ny, 2
       S%P((i-1)*nx+1:i*nx) =  -S%P(1:nx)
    end do

    do i = 3, ny, 2
       S%P((i-1)*nx+1:i*nx) = S%P(1:nx)
    end do

    do i = 2, nz, 2
       S%P((i-1)*nxy+1:i*nxy) =  -S%P(1:nxy)
    end do

    do i = 3, nz, 2
       S%P((i-1)*nxy+1:i*nxy) = S%P(1:nxy)
    end do

    ! Fourier Transformation
    allocate(S%FT(S%nClass, S%nClass))
    S%FT = ZERO
    do i = 1, S%nClass
       do jx = 0, nx - 1
          do jy = 0, ny - 1
             do jz = 0, nz - 1
                j = S%D(jz*nxy+jy*nx+jx+1, 1)
                S%FT(i,j) = S%FT(i,j) + &
                     cos(TWOPI*(cord(i,1)*jx/nx + &
                     &          cord(i,2)*jy/ny + &
                     &          cord(i,3)*jz/nz ))
             end do
          end do
       end do
    end do
    
    ! Enable the flag
    S%checklist = .true.

    !!! No waves and bonds
    S%checklist(STRUCT_WAVE) = .false.
    S%checklist(STRUCT_BOND) = .false.

  end subroutine DQMC_INIT_3D

end module DQMC_3D
 
