module DQMC_RECLATT

  use DQMC_GEOM_PARAM
  use DQMC_LATT

  implicit none

  type :: recip_lattice_t

     integer                :: ndim
     integer                :: nkpts             !number of k-points (equal to ncell)
     real*8, pointer        :: klist(:,:)        !list of k-points(nkpts,rdim)
     real*8                 :: kpoint(rdim)      !Input k-point 
     real*8                 :: ktwist(rdim)      !Twist vector
     real*8                 :: kcs(rdim,rdim)    !cartesian component of reciprocal superlattice***
     real*8                 :: ks(rdim,rdim)     !fractional components of reciprocal superlattice*** 
     real*8                 :: kc(rdim,rdim)     !cartesian components of reciprocal unit cell***
     !*** rows of these matrices are the vectors

     integer                :: nclass_k          !number of inequivalent k-points
     integer, pointer       :: myclass_k(:)      !class for each k-point (nkpts)
     integer, pointer       :: class_size_k(:)   !number of equivalent k-points in each class (nclass_k)
     integer, pointer       :: class_repr_k(:)   !representative k-point for each class (nclass_k)

     logical                :: initialized
     logical                :: constructed
     logical                :: analyzed
  end type recip_lattice_t

contains

  !---------------------------------------------------------------------
  ! Read and fill most of the variables that define the lattices in
  ! real and reciprocal space.
  !---------------------------------------------------------------------
  subroutine init_recip_latt(lattice,recip_lattice) 
    integer                                     :: ndim,i,j
    real*8                                      :: projk(rdim)
    real*8, pointer                             :: kc(:,:),kcs(:,:),ktwist(:),kpoint(:)
    real*8, pointer                             :: ac(:,:),scc(:,:)
    type(lattice_t),intent(in),target           :: lattice
    type(recip_lattice_t),intent(out),target    :: recip_lattice
    logical                                     :: ldum

    if(.not.lattice%initialized)stop'Need to initialize lattice before recip_lattice'

    !alias arrays
    scc=>lattice%scc
    ac=>lattice%ac
    kc=>recip_lattice%kc
    kcs=>recip_lattice%kcs
    ktwist=>recip_lattice%ktwist
    kpoint=>recip_lattice%kpoint

    ndim=lattice%ndim

    !k-point in fractionary units (e.g. in a square lattice 0.5 0.5 will correspond to pi pi )
    kpoint(:)=0.d0
    ldum=move_to_record(INPUT_FIELDS(KPOINT_F),inpunit)
    if(ldum)then
       read(inpunit,*,iostat=i)(kpoint(i),i=1,ndim)
       if(i/=0)stop 'Problem reading K-Point'
    endif

    !Since QUEST cannot handle complex hoppings we stop if the k-point is /=0.
    if(sum(kpoint**2)>toll)then
       write(*,*)'KPOINTS/=0 are not implemented in DQMC yet'
       stop
    endif

    !Take care of k-space cells
    !Rows of kc are reciprocal lattice basis vectors in cartesian coordinates
    call get_inverse(ac,kc)
    kc(:,:)=kc*2.d0*pi

    !transform k-point in cartesian coordinate
    do j=1,rdim 
       ktwist(j)=sum(kpoint(:)*kc(j,:))
    enddo
    kpoint(:)=ktwist(:)

    !find reciprocal lattice vectors of the supercell
    call get_inverse(scc,kcs)
    kcs(:,:)=kcs(:,:)*2.d0*pi

    !and project it inside the reciprocal supercell
    do i=1,rdim
       projk(i)=sum(kpoint(:)*scc(:,i))/(2.0*pi)
    enddo
    projk(:)=projk(:)-nint(projk(:))
    !then find its cartesian coordinates
    do i=1,rdim
       ktwist(i)=sum(projk(:)*kcs(:,i))
    enddo

    !find components of kcs in units of kc
    do i=1,3 
       do j=1,3
          recip_lattice%ks(i,j)=sum(kcs(i,:)*ac(:,j))/(2.d0*pi)
       enddo
    enddo

    recip_lattice%ndim=ndim
    recip_lattice%nkpts=lattice%ncell

    recip_lattice%initialized=.true.

    !write to stdout
    write(*,*)
    write(*,*)' Reciprocal lattice basis vectors'
    write(*,'(3f14.7)')((kc(i,j),j=1,rdim),i=1,rdim)
    write(*,*)
    write(*,*)' Reciprocal super-lattice basis vectors'
    write(*,'(3f14.7)')((kcs(i,j),j=1,rdim),i=1,rdim)
    write(*,*)
    write(*,*)'Original k-point'
    write(*,'(3f14.7)')(kpoint(i),i=1,rdim)
    write(*,*)
    write(*,*)'Twist vector'
    write(*,'(3f14.7)')(ktwist(i),i=1,rdim)
    write(*,*)'================================================================'

  end subroutine init_recip_latt

  !----------------------------------------------------------------------------
  ! construct the reciprocal lattice (klist).
  !----------------------------------------------------------------------------
  subroutine construct_recip_lattice(recip_lattice)
    integer                :: nkx,nky,nkz,ndim,ikx,iky,ikz,ikv(rdim),ncubex,ncubey,ncubez, &
         & ilist,nedge,nfound,nkpts,i,j
    integer*8, allocatable :: indedge(:)
    real*8                 :: invkc(rdim,rdim)
    real*8, allocatable    :: kset(:,:)
    real*8, pointer        :: klist(:,:)
    type(recip_lattice_t), intent(inout)  :: recip_lattice

    if(.not.recip_lattice%initialized)stop'Need to initialize recip_lattice before construction'

    call get_inverse(recip_lattice%kc,invkc)
    !ktwist is the vector contained in the 1st BZ of the reciprocal 
    !superlattice. It is the true twist vector. To find it:
    !First project vector inside kcs...
    !and finally find its fractional coordinate
    !Define kset as the set of points surrouding the origin.
    ndim=recip_lattice%ndim
    allocate(klist(recip_lattice%nkpts,rdim))
    nkx=0; if(ndim>0)nkx=1 
    nky=0; if(ndim>1)nky=1 
    nkz=0; if(ndim>2)nkz=1
    allocate(kset(3**ndim-1,rdim),indedge(recip_lattice%nkpts))
    i=0
    do ikx=-nkx,nkx; ikv(1)=ikx
       do iky=-nky,nky; ikv(2)=iky
          do ikz=-nkz,nkz; ikv(3)=ikz
             if(sum(ikv**2)>0)then
                i=i+1 
                do j=1,rdim
                   kset(i,j)=sum(ikv(:)*recip_lattice%kc(:,j))
                enddo
             endif
          enddo
       enddo
    enddo
    !write(*,'(3f12.6)')((kset(i,j),j=1,3),i=1,3**ndim-1)
    ncubex=0; ncubey=0; ncubez=0
    ilist=1; klist(ilist,:)=recip_lattice%ktwist(:); nedge=0
    if(ndim>0)then
       do 
          ncubex=ncubex+1
          if(ndim>1)ncubey=ncubey+1
          if(ndim>2)ncubez=ncubez+1
          nfound=0
          do ikz=-ncubez,ncubez
             ikv(3)=ikz
             !constant y edges
             if(ncubey/=0)then
                do iky=-ncubey,ncubey,2*ncubey
                   ikv(2)=iky
                   do ikx=-ncubex,ncubex
                      ikv(1)=ikx
                      call check_and_update
                   enddo
                enddo
             endif
             !write(*,*)'NFOUND',nfound
             !constant x edges
             do ikx=-ncubex,ncubex,2*ncubex
                ikv(1)=ikx
                do iky=min(0,-ncubey+1),max(0,ncubey-1)
                   ikv(2)=iky
                   call check_and_update
                enddo
             enddo
             !write(*,*)'NFOUND',nfound
          enddo
          !constant kz faces
          if(ncubez>0)then
             do ikz=-ncubez,ncubez,2*ncubez
                ikv(3)=ikz
                do iky=-ncubey+1,ncubey-1
                   ikv(2)=iky
                   do ikx=-ncubex+1,ncubex-1
                      ikv(1)=ikx
                      call check_and_update
                   enddo
                enddo
             enddo
          endif
          if(nfound==0)exit
       enddo
    endif
    nkpts=ilist
    call vsort(klist,nkpts)
    recip_lattice%klist=>klist
    write(*,*)'Reciprocal lattice (1st Brillouin zone)'
    write(*,'(/,A,1x,i4)')' Number of k-points found: ',nkpts
    if(nkpts/=recip_lattice%nkpts)then
       write(*,*)'This is different from', recip_lattice%nkpts
       stop 
    endif
    write(*,'(3f12.6)')((klist(i,j),j=1,rdim),i=1,recip_lattice%nkpts)
    write(*,*)'================================================================'

    recip_lattice%constructed=.true.

  contains

    subroutine  check_and_update
      integer i,ie
      real*8 :: kpt(rdim),projk(rdim),diff(rdim)
      logical on_edge,included
      do i=1,3
         kpt(i)=sum(ikv(:)*recip_lattice%kcs(:,i))+recip_lattice%ktwist(i)
      enddo
      if(closer_to_zero(kpt,kset,on_edge,ndim))then
         included=.false.
         if(on_edge)then
            do ie=1,nedge
               diff(:)=klist(indedge(ie),:)-kpt(:)
               do i=1,3
                  projk(i)=sum(diff(:)*invkc(:,i))
               enddo
               projk(:)=projk(:)-nint(projk(:))
               if(sum(projk**2)<1.d-6)then
                  included=.true.
               endif
            enddo
         endif
         if(.not.included)then
            nfound=nfound+1
            ilist=ilist+1
            !write(*,'(3f12.6)')kpt(:)
            klist(ilist,:)=kpt(:)
            if(on_edge)then
               nedge=nedge+1
               indedge(nedge)=ilist
            endif
         endif
      endif
    end subroutine check_and_update

    subroutine vsort(vvec,n)
      integer n,i,j,k
      real*8 vvec(n,rdim),vlen(n),av(rdim),a
      do i=1,n
         vlen(i)=sum(vvec(i,:)**2)
      enddo
      do i=2,n
         a=vlen(i)
         av(:)=vvec(i,:)
         do j=1,i-1
            if(a+1.d-10<vlen(j))exit
         enddo
         do k=i-1,j,-1
            vlen(k+1)=vlen(k)
            vvec(k+1,:)=vvec(k,:)
         enddo
         vlen(j)=a
         vvec(j,:)=av(:)
      enddo
    end subroutine vsort

  end subroutine construct_recip_lattice


  !----------------------------------------------------------
  !Determine whether a k-point (kpt) is closer to k=0 than
  !to any of the points in kset. If half way set to true
  !but also set on_edge to true.
  !----------------------------------------------------------
  logical function closer_to_zero(ktp,kset,on_edge,ndim)
    integer, intent(in) :: ndim
    real*8, intent(in)  :: ktp(rdim),kset(3**ndim-1,rdim)
    real*8              :: diff(rdim),dist0,dist
    logical             :: on_edge
    integer             :: i
    on_edge=.false.
    diff(:)=ktp(:)**2
    dist0=sum(diff)
    do i=1,3**ndim-1
       diff(:)=(ktp(:)-kset(i,:))**2
       dist=sum(diff)
       if(dist+toll<dist0)then
          closer_to_zero=.false.
          return
       else
          if(dist-toll<dist0)on_edge=.true.
       endif
    enddo
    closer_to_zero=.true.
  end function closer_to_zero



end module DQMC_RECLATT

