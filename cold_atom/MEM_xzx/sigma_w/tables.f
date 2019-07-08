      subroutine tables
c****************************************************************************
c     This subroutine makes the lookup tables which are used
c     throughout the rest of the codes.  These tables include:
c     wn(n)
c     
c****************************************************************************
      use Global
      implicit none
c****************************************************************************
      integer :: n
c**************************************************************************
      
      N_sp=sqrt(float(Nc))
      
      
c*****************************************************************************
c     MATSUBARA FREQUENCIES	  
c*****************************************************************************1536

      do n=-nwn,nwn
        wn(n)=(2.0_kind*n+1.0_kind)*pi*temp
      end do

      if(ndim.eq.2) then
         call tables_2D_Betts
      else if(ndim.eq.3) then
         call tables_3D_Betts
      else
         stop 'bad ndim in tables'
      endif

      end subroutine tables

      subroutine tables_2D_Betts      ! TWO DIMENSIONAL BETTS TABLES
c***************************************************************************
c     CHANGES REQUIRED FOR USE IN QMC CODES:
c     1. remove comment(s) in from of
c          use Global
c          define the  fmat matrix
c          call allocate_2p
c     2. change program and end statements to subroutine and return
c     3. delete the "Stuff that belongs in module_Global"...
c     4. remove "write(6,*) 'enter cluster'"...
c     5. set  use_brute_force = .false. and use_lin_alg = .true.
c     7. remove these comments
c     CHANGES REQUIRED FOR USE IN Analysis CODES:
c     1. do the steps above
c     2. remove the delcarations of ict, ig, jg and kg
c***************************************************************************
c     This subroutine makes the lookup tables used throughout the
c     full TWO DIMENSIONAL code.  It employs the definitions and 
c     cluster geometries developed by D.D. Betts, H.Q. Lin, and 
c     J.S. Flynn Canadian Journal of Physics, 77, 353-369.
c
c     Developers include M. Jarrell (symmetry tables), P. Kent (fast 
c     determination of allowed group operations), Th. Maier (WS cells of 
c     the cluster superlattice).  We reserve all rights to this code.
c***************************************************************************
      use Global
      implicit none
      interface transform2D
      function gtransform2D(a,igroup,ixy,isx,isy)
      integer :: a(3,3),gtransform2D(3,3),igroup,ixy(8),isx(8),isy(8)
      end
      
      function Rtransform2D(Rc,igroup,ixy,isx,isy)
      integer :: Rc(3),Rtransform2D(3),igroup,ixy(8),isx(8),isy(8)
      end
      
      function Ktransform2D(Kc,igroup,ixy,isx,isy)
      integer :: igroup, ixy(8),isx(8),isy(8)
      real(8) :: Kc(3),Ktransform2D(3)
      end
      end interface 

c***************************************************************************
c     Locally defined
c***************************************************************************
      logical    :: allowed,denied,denied_lin_alg
      integer, parameter :: nshellmax=50
      integer :: i,i1,i2,i3,ix,iy,it,igroup,ind,ic,ic1,ic2,ick,idim,ik,
     &           itest,j,jmd,jm,j1,j2,jc,jk,jkc,jlow,jhigh,
     &           k,k2,kroniker,m,m1,m2,n,n1,n2,l,p,
     &           ixy(8),isx(8),isy(8),iv1(3),iv2(3),iv3(3),irwgroup(64),
     &           a1(3,3),group(8),
     &           ig(64),jg(64),kg(8),
     &           neighbor(3,3),
     &           nneigh(nshellmax),nneighl(nshellmax),ishellmax,I_F
      real(kind) :: r1,r2,r3,r4,r5,Dkt,c1(3,3),c2(3,3),Dwavebar(64),
     &              rv1(3),rv2(3),rv3(3),Ai(3,3)
c For lin algebra test
      real(kind) :: ainv(3,3),a1ainv(3,3)
      complex(kind) :: cr1,Dwave(64,64)
c Select use of "brute force" vs linear algebra method for 
c checking lattice coincidence. Both methods should agree; lin alg
c method is much faster.
      logical, parameter :: use_brute_force = .false.
      logical, parameter :: use_lin_alg = .true.
c***************************************************************************
c     Local arrays
c     ixy(igroup)      D=2 square point group table(see below)
c     isx(igroup)      D=2 square point group table(see below)
c     isy(igroup)      D=2 square point group table(see below)
c     group(i)         Indices of allowed group operations
c     irwgroup(ic)     the group operation maps K(ic) into IRW
c     ig(ic)           i in the g1 and j in the g2 direction
c     jg(jc)           corresponding to K point ic
c     ict(i,j)         maps i*g1+j*g2 to ic
c     Ai(i,j)          The inverse transpose of a (see below)
c     the rest are just temporary arrays.
c*************************************************************************
      
c     CLUSTER MOMENTUM and DUAL SPACE for the 2D square lattice 
c                                                           _   _ 
c     Different Real-space Cluster Geometries.  The vectors a1, a2
c     define the parallelpipeds which define the periodically replicated 
c     clusters.  Here, we list only the clusters from table 1A with Nc
c     even.
c              _          _         
c  cluster     a1         a2   bipartide  IF  IB  used
c    1A     ( 1, 0)    ( 0, 1)    no               no
c
c    2A     ( 1, 1)    (-1, 1)    yes              no
c
c    4A     ( 2, 0)    ( 0, 2)    yes              no
c
c    8A     ( 2, 2)    (-2, 2)    yes      0   0   yes
c    8C     ( 2, 3)    (-2, 1)    no       0   -   no
c
c   10A     ( 2, 3)    (-2, 2)    no       0   -   no
c   10B     ( 1, 3)    (-3, 1)    yes      1   0   no
c
c   12A     ( 2, 2)    (-4, 2)    yes      2   0   no
c   12B     ( 2, 3)    (-2, 3)    no       0   -   no
c   12D     ( 1, 3)    ( 4, 0)    yes      2   0   no
c
c   14A     ( 3, 2)    (-1, 4)    no       0   -   no
c   14B     ( 4, 2)    (-1, 3)    yes      2   0   no
c
c   16A     ( 4, 2)    ( 0, 4)    yes      1   0   yes
c   16B     ( 4, 0)    ( 0, 4)    yesI     3   1   no
c   16C     ( 3, 2)    (-2, 4)    no       0   -   no
c   16D     ( 4, 1)    ( 0, 4)    no       1   -   no
c   16E     ( 1, 5)    (-3, 1)    yesI     3   1   no
c
c   18A     ( 3, 3)    (-3, 3)    yes      0   0   yes
c   18B     ( 3, 3)    (-2, 4)    yes      0   0   yes
c   18C     ( 4, 2)    (-1, 4)    no       0   -   no
c   18D     ( 5, 3)    (-1, 3)    yesI     4   4   no
c
c   20A     ( 4, 2)    (-2, 4)    yes      1   0   yes
c   20B     ( 3, 4)    (-2, 4)    no       0   -   no
c   20C     ( 1, 4)    (-5, 0)    no       1   -   no
c   20E     ( 5, 1)    ( 0, 4)    yesI     3   1   no
c
c   22A     ( 2, 4)    (-5, 1)    yes      2   0   yes
c   22B     ( 4, 2)    (-3, 4)    no       0   -   no
c   22D     ( 1, 4)    (-5, 2)    no       1   -   no
c
c   24A     ( 4, 4)    (-4, 2)    yes      3   0   yes
c   24B     ( 4, 0)    ( 2, 6)    yesI     5   1   no
c   24C     ( 4, 0)    ( 0, 6)    yesI     7   2   no
c   24D     ( 1, 5)    ( 5, 1)    yes      3   0   yes
c   24E     ( 3, 3)    (-5, 3)    yes      3   0   yes
c   24F     ( 4, 4)    (-3, 3)    yes      3   0   yes
c   24G     ( 4, 3)    (-4, 3)    no       0   -   no
c                                   
c   26A     ( 4, 2)    (-3, 5)    yes      3   0   yes
c   26B     ( 1, 5)    (-5, 1)    yesI     5   1   yes
c   26C     ( 5, 2)    (-3, 4)    no       0   -   no
c
c   28A     ( 2, 4)    (-6, 2)    yesI     4   1   no
c   28B     ( 5, 3)    (-1, 5)    yes      1   0   yes
c   28C     ( 2, 5)    (-4, 4)    no       0   -   no
c   28D     ( 4, 3)    (-4, 4)    no       0   -   no
c
c   30A     ( 5, 3)    (-5, 3)    yes      1   0   yes
c   30B     ( 4, 3)    (-2, 6)    no       0   -   no
c   30C     ( 2, 5)    (-6, 0)    no       1   -   no
c   30D     ( 1, 5)    (-6, 0)    yesI     5   2   no
c   30G     ( 5, 5)    (-2, 4)    yesI     5   2   no
c
c   32A     ( 4, 4)    (-4, 4)    yes      0   0   yes
c   32B     ( 2, 6)    (-4, 4)    yes      0   0   yes
c   32C     ( 3, 5)    (-4, 4)    yes      0   0   no
c   32D     ( 1, 6)    (-5, 2)    no       1   -   no
c   32G     ( 3, 4)    (-5, 4)    no       0   -   no
c 
c   36A     ( 6, 0)    ( 0, 6)    yes
c
c   64A     ( 8, 0)    ( 0, 8)    yes
c
c  used determined form Table 3, bipartide from Ib in table 1A, yesI indicates
c  that the cluster is imperfect
c
c  We need to automate the process of finding which of the 8 point group 
c  operations are retained by the different clusters.  To do this, we must 
c  perform all a cubic group operations on the {a} and see if they form 
c  an equivalent parallelpiped.  
c
c  Group Operations for a square lattice:
c
c  x --> +/- x                  There are 2^D D! operations in D dimensions
c    \/                         or 8 operations in 2D.
c    /\
c  y --> +/- y   
c        .
c  a'=transform2D(a,igroup,ixy,isx,isy)
c  ixy=0(1) dont (do) exchange x and y
c  isx      sign of x 
c  isy      sign of y
c
c  x <-- isx*[(1-ixy)*x + ixy*y] 
c  y <-- isy*[(1-ixy)*y + ixy*x] 
c

      ixy( 1)=0; isx( 1)=+1; isy( 1)=+1 !identity  xyz e
      ixy( 2)=0; isx( 2)=-1; isy( 2)=+1 !
      ixy( 3)=0; isx( 3)=+1; isy( 3)=-1 !
      ixy( 4)=0; isx( 4)=-1; isy( 4)=-1 !
      ixy( 5)=1; isx( 5)=+1; isy( 5)=+1 !
      ixy( 6)=1; isx( 6)=-1; isy( 6)=+1 !
      ixy( 7)=1; isx( 7)=+1; isy( 7)=-1 !
      ixy( 8)=1; isx( 8)=-1; isy( 8)=-1 !
      
c***********************************************************************        
c       determine the principle translation vectors a1, and a2
c***********************************************************************
           If(cluster.eq.'1A') then; a(1,1:2)=(/1,0/); a(2,1:2)=(/ 0,1/)
      else if(cluster.eq.'2A') then; a(1,1:2)=(/1,1/); a(2,1:2)=(/-1,1/)
      else if(cluster.eq.'4A') then; a(1,1:2)=(/2,0/); a(2,1:2)=(/ 0,2/)
      else if(cluster.eq.'8A') then; a(1,1:2)=(/2,2/); a(2,1:2)=(/-2,2/)
      else if(cluster.eq.'8C') then; a(1,1:2)=(/2,3/); a(2,1:2)=(/-2,1/)
      else if(cluster.eq.'10A')then; a(1,1:2)=(/2,3/); a(2,1:2)=(/-2,2/)
      else if(cluster.eq.'10B')then; a(1,1:2)=(/1,3/); a(2,1:2)=(/-3,1/)
      else if(cluster.eq.'12A')then; a(1,1:2)=(/2,2/); a(2,1:2)=(/-4,2/)
      else if(cluster.eq.'12B')then; a(1,1:2)=(/2,3/); a(2,1:2)=(/-2,3/)
      else if(cluster.eq.'12D')then; a(1,1:2)=(/1,3/); a(2,1:2)=(/ 4,0/)
      else if(cluster.eq.'14A')then; a(1,1:2)=(/3,2/); a(2,1:2)=(/-1,4/)
      else if(cluster.eq.'14B')then; a(1,1:2)=(/4,2/); a(2,1:2)=(/-1,3/)
      else if(cluster.eq.'16A')then; a(1,1:2)=(/4,2/); a(2,1:2)=(/ 0,4/)
      else if(cluster.eq.'16B')then; a(1,1:2)=(/4,0/); a(2,1:2)=(/ 0,4/)
      else if(cluster.eq.'16C')then; a(1,1:2)=(/3,2/); a(2,1:2)=(/-2,4/)
      else if(cluster.eq.'16D')then; a(1,1:2)=(/4,1/); a(2,1:2)=(/ 0,4/)
      else if(cluster.eq.'16E')then; a(1,1:2)=(/1,5/); a(2,1:2)=(/-3,1/)
      else if(cluster.eq.'18A')then; a(1,1:2)=(/3,3/); a(2,1:2)=(/-3,3/)
      else if(cluster.eq.'18B')then; a(1,1:2)=(/3,3/); a(2,1:2)=(/-2,4/)
      else if(cluster.eq.'18C')then; a(1,1:2)=(/4,2/); a(2,1:2)=(/-1,4/)
      else if(cluster.eq.'18D')then; a(1,1:2)=(/5,3/); a(2,1:2)=(/-1,3/)
      else if(cluster.eq.'20A')then; a(1,1:2)=(/4,2/); a(2,1:2)=(/-2,4/)
      else if(cluster.eq.'20B')then; a(1,1:2)=(/3,4/); a(2,1:2)=(/-2,4/)
      else if(cluster.eq.'20C')then; a(1,1:2)=(/1,4/); a(2,1:2)=(/-5,0/)
      else if(cluster.eq.'20E')then; a(1,1:2)=(/5,1/); a(2,1:2)=(/ 0,4/)
      else if(cluster.eq.'22A')then; a(1,1:2)=(/2,4/); a(2,1:2)=(/-5,1/)
      else if(cluster.eq.'22B')then; a(1,1:2)=(/4,2/); a(2,1:2)=(/-3,4/)
      else if(cluster.eq.'22D')then; a(1,1:2)=(/1,4/); a(2,1:2)=(/-5,2/)
      else if(cluster.eq.'24A')then; a(1,1:2)=(/4,4/); a(2,1:2)=(/-4,2/)
      else if(cluster.eq.'24B')then; a(1,1:2)=(/4,0/); a(2,1:2)=(/ 2,6/)
      else if(cluster.eq.'24C')then; a(1,1:2)=(/4,0/); a(2,1:2)=(/ 0,6/)
      else if(cluster.eq.'24D')then; a(1,1:2)=(/1,5/); a(2,1:2)=(/ 5,1/)
      else if(cluster.eq.'24E')then; a(1,1:2)=(/3,3/); a(2,1:2)=(/-5,3/)
      else if(cluster.eq.'24F')then; a(1,1:2)=(/4,4/); a(2,1:2)=(/-3,3/)
      else if(cluster.eq.'24G')then; a(1,1:2)=(/4,3/); a(2,1:2)=(/-4,3/)
      else if(cluster.eq.'26A')then; a(1,1:2)=(/4,2/); a(2,1:2)=(/-3,5/)
      else if(cluster.eq.'26B')then; a(1,1:2)=(/1,5/); a(2,1:2)=(/-5,1/)
      else if(cluster.eq.'26C')then; a(1,1:2)=(/5,2/); a(2,1:2)=(/-3,4/)
      else if(cluster.eq.'28A')then; a(1,1:2)=(/2,4/); a(2,1:2)=(/-6,2/)
      else if(cluster.eq.'28B')then; a(1,1:2)=(/5,3/); a(2,1:2)=(/-1,5/)
      else if(cluster.eq.'28C')then; a(1,1:2)=(/2,5/); a(2,1:2)=(/-4,4/)
      else if(cluster.eq.'28D')then; a(1,1:2)=(/4,3/); a(2,1:2)=(/-4,4/)
      else if(cluster.eq.'30A')then; a(1,1:2)=(/5,3/); a(2,1:2)=(/-5,3/)
      else if(cluster.eq.'30B')then; a(1,1:2)=(/4,3/); a(2,1:2)=(/-2,6/)
      else if(cluster.eq.'30C')then; a(1,1:2)=(/2,5/); a(2,1:2)=(/-6,0/)
      else if(cluster.eq.'30D')then; a(1,1:2)=(/1,5/); a(2,1:2)=(/-6,0/)
      else if(cluster.eq.'30G')then; a(1,1:2)=(/5,5/); a(2,1:2)=(/-2,4/)
      else if(cluster.eq.'32A')then; a(1,1:2)=(/4,4/); a(2,1:2)=(/-4,4/)
      else if(cluster.eq.'32B')then; a(1,1:2)=(/2,6/); a(2,1:2)=(/-4,4/)
      else if(cluster.eq.'32C')then; a(1,1:2)=(/3,5/); a(2,1:2)=(/-4,4/)
      else if(cluster.eq.'32D')then; a(1,1:2)=(/1,6/); a(2,1:2)=(/-5,2/)
      else if(cluster.eq.'32G')then; a(1,1:2)=(/3,4/); a(2,1:2)=(/-5,4/)
      else if(cluster.eq.'36A')then; a(1,1:2)=(/6,0/); a(2,1:2)=(/0,6/)
      else if(cluster.eq.'64A')then; a(1,1:2)=(/8,0/); a(2,1:2)=(/0,8/)
      else
         if(myrank.eq.0) write(lud,*) 'Bad cluster specification' !0
         stop
      end if  
      if(iprint.ge.5.and.myrank.eq.0) 
     &   write(lud,"('a=',2('(',i2,1x,i2,')'))") a(1,1:2),a(2,1:2)  

c*********************************************************************************
c     Determine which group operations are allowed    
c*********************************************************************************

      ngroup=0
      if(myrank.eq.0.and.iprint.ge.5) then
        write(lud,*) '----------------------------------'
        write(lud,*) '  trans(a)   igroup L  ixy isx isy'
        write(lud,*) '----------------------------------'
      end if
      do igroup=1,8
         a1=transform2D(a,igroup,ixy,isx,isy)
          
          
         if (use_brute_force) then

c        look at the lattice points generated by a1 and see if they
c        correspond to those produced by a.
c        First, generate the a1 lattice points
         do i=-Nc,Nc
         do j=-Nc,Nc
           denied=.true.
           iv1=(/i,j,0/)
           do l=1,2         ! x,y components of this lattice point
             rv1(l)=sum(iv1(1:2)*a1(1:2,l))
           end do
c          Now generate the a lattice points.  If a and a1 are the 
c          same lattice one of these points MUST correspond to ij
           do n=-6*Nc,6*Nc
           do m=-6*Nc,6*Nc
             iv2=(/n,m,0/)
             do l=1,2         ! x,y components of this lattice point
               rv2(l)=sum(iv2(1:2)*a(1:2,l))
             end do
             if(sum(abs(rv1-rv2)).lt.eps) denied=.false.
           end do
           end do
           if(denied) exit !no corrs. point was found
         end do
         if(denied) exit
         end do
         end if

         if (use_lin_alg) then
c PK If lattice points generated by "a1" vectors can be found in 
c lattice points generated by "a" vectors then (a1)x(a^-1) can only 
c contain integers.
            
            r1 = a(1,1)*a(2,2) - a(2,1)*a(1,2) ! the determinent of a
            ainv=0.
            ainv(1,1) = +a(2,2)/r1
            ainv(1,2) = -a(1,2)/r1       
            ainv(2,1) = -a(2,1)/r1
            ainv(2,2) = +a(1,1)/r1
            a1ainv=matmul(a1,ainv)
            if (all(abs(a1ainv-nint(a1ainv)).le.eps)) then 
               denied_lin_alg=.false. ! matrix only contains integers
            else
               denied_lin_alg=.true.
            end if
         end if

c If both algorithms are enabled, check they agree
         if (use_brute_force.and.use_lin_alg) then
            if (denied.neqv.denied_lin_alg) then
               write (*,*) '*** BUG : Denied algorithm disagree'
               write (*,*) 'a11',a1(:,1)
               write (*,*) 'a12',a1(:,2)
               write (*,*) 'a13',a1(:,3)
               
               write (*,*) 'a1',a(:,1)
               write (*,*) 'a2',a(:,2)
               write (*,*) 'a3',a(:,3)
               write (*,*)
               write (*,*) 'ainv',ainv(:,1)
               write (*,*) 'ainv',ainv(:,2)
               write (*,*) 'ainv',ainv(:,3)
               write (*,*)
               
               write (*,*) 'a1ainv',a1ainv(:,1)
               write (*,*) 'a1ainv',a1ainv(:,2)
               write (*,*) 'a1ainv',a1ainv(:,3)
               write (*,*) 'Denied_lin_alg=',denied_lin_alg
               write (*,*) 'Denied_brute_force=',denied
               write (*,*)
               stop
            end if
         end if
         if (use_lin_alg) denied=denied_lin_alg


         if(denied.eqv..false.) then ! a corresponding point was found
            ngroup=ngroup+1
            group(ngroup)=igroup
         end if    
         if(myrank.eq.0.and.iprint.ge.5) then
           allowed=.true.;if (denied) allowed=.false.
           write(lud,"(2('(',i2,1x,i2,')'),2x,i2,2x,l1,3(2x,i2))")
     &     a1(1,1:2),a1(2,1:2),igroup,allowed,ixy(igroup),
     &     isx(igroup),isy(igroup)
         end if
      end do      
     
      if(myrank.eq.0.and.iprint.ge.5) then
        write(lud,*) '----------------------------------'
        write(lud,"('cluster=',a3,' ngroup=',i2,'  groups=',8(i2,1x))") 
     &  cluster,ngroup,group(1:ngroup)
      end if
     




c******************************************************************************
c       calculate the table of cluster point locations Rc(i,ic), i=x,y
c******************************************************************************
c                                    a1 and a2 need not be orthogonal.    
c  ^ y                               Therefore, if v=A n, where n is a
c  |                                 vector containing the projection of
cc-|-c---c-*-c---c---c---c           v onto a, and                          
c| | |   |   |   |   |   |                                                  
cc-|-c---c---c---c---c---c           A= [ a1x  a1y ]  a= [ a1x  a2x ]       
ca2^ * * * * * * * * |   |              [ a2x  a2y ]     [ a1y  a2y ]       
cc-|-c---c---c---c-*-c---c                                                  
c| | |   P   |   | * |   |           which is the transpose of a, then the  
cc-|-c--/c---c---c-*-c---c           projection of v onto a is              
c| | | /v|   |   | * |   |                  -1                              
cc-|-c/--c---c---c-*-c---c             n = A  v,  where                     
c| | /   |   |   | * |   |                                                  
cc-|/0---c---c---c-*-c---c             -1  [ a2y -a2x ] /                   
c| O--------------->---------> x      A  = [-a1y  a1x ]/ (a1x a2y-a2x a1y)  
cc---c---c---c---ca1-c---c                                                  
c|   |   |   |   |   |   |           v is within the cluster, if the        
cc---c---c---c---c---c---c           elements of n, -eps < n_i < 1-eps      
c
c                                                                          

      Ai=zeror
      r1 = a(1,1)*a(2,2) - a(2,1)*a(1,2) ! the determinent of a
      Ai(1,1) = +a(2,2)/r1
      Ai(1,2) = -a(2,1)/r1       
      Ai(2,1) = -a(1,2)/r1
      Ai(2,2) = +a(1,1)/r1

      if(iprint.ge.5.and.myrank.eq.0) then
        write(lud,*) '  '                   !0
        write(lud,*) 'ic,Rcx(ic),Rcy(ic)'   !0
      end if
      ic=0
      do j=-Nc,Nc  ! Rc(:,1) must be the origin!  Since a2x may be <0
      do i=-Nc,Nc  ! we need to index i (x) faster than j (y).
        rv1=(/i,j,0/)
        rv2=matmul(Ai,rv1)
        if(rv2(1).lt.oner-eps.and.rv2(1).gt.zeror-eps.and.
     &     rv2(2).lt.oner-eps.and.rv2(2).gt.zeror-eps) then
           ic=ic+1
           Rc(:,ic)=(/i,j,0/) ! Rc, like i.e. a(1,1) is integer
            if(iprint.ge.5.and.myrank.eq.0) 
     &           write(lud,"(i2,4x,i3,4x,i3)") ic,Rc(1:2,ic)   !0
           if(sum(abs(Rc(1:2,ic))).eq.0) icR0=ic
           if(sum(Rc(1:2,ic)**2).eq.1) icR1=ic
         end if
      end do
      end do
      if(ic.ne.Nc) then
         if(myrank.eq.0) write(lud,*) 'Bug in Rc(ic) loop'   !0
         stop
      end if
      if(iprint.ge.5.and.myrank.eq.0) write(lud,*) 'icR0=',icR0

c***********************************************************************        
c       create a table mapping points outside the cluster back into it.
c***********************************************************************        
        

      if(iprint.ge.5.and.myrank.eq.0) then
        write(lud,*) '  '                                !0
        write(lud,*) ' Equivalency table r'              !0
        write(lud,*) ' ic  igroup   icrequ(ic,ieqiv)'    !0
      end if
      do ic=1,Nc
        do igroup=1,ngroup
c         Generate the equivalent points using the ngroup operations
          iv1=transform2D(Rc(:,ic),group(igroup),ixy,isx,isy) ! cartesian coordinates
          do n=1,3            
            rv2=matmul(Ai,iv1) ! projection of iv1 onto a1 and a2
            do l=1,2           ! move iv1 and rv2 into the cluster
              if(rv2(l).lt.zeror-eps) then
                 iv1(:)=iv1(:)+a(l,:)
              else if(rv2(l).gt.oner-eps) then
                 iv1(:)=iv1(:)-a(l,:)
              end if
            end do
          end do
c         Figure out which cluster point iv1 is
          itest=0
          do jc=1,Nc
             if(sum(abs(iv1(1:2)-Rc(1:2,jc))).eq.0) then ! iv1=Rc(jc)
                icrequ(ic,igroup)=jc
                itest=1
             end if
          end do
          if(itest.eq.0) then
             write(lud,*) 'No equivalent found in icrequ',ic,igroup
             stop
          end if 
          if(iprint.ge.5.and.myrank.eq.0) 
     &      write(lud,"(2x,i2,5x,i1,8x,i2)") 
     &      ic,igroup,icrequ(ic,igroup) !0
        end do
      end do                               


c***********************************************************************        
c     Now create the tables for the differences of R
c***********************************************************************        
      if(iprint.ge.5.and.myrank.eq.0) then
        write(lud,*) '  '                      !0
        write(lud,*) ' differences in r'       !0
        write(lud,*) ' ic  jc  icrdiff(ic,jc)' !0
      end if
      do i=1,Nc
      do j=1,Nc
        iv1(:)=Rc(:,i)-Rc(:,j) 
        itest=0          ! We dont know where iv1 is.
        do n=-1,1        ! Look in the surrounding 
        do m=-1,1        ! 9 clusters for iv1.
          iv2=(/n,m,0/)
          do l=1,2
            iv3(l)=iv1(l)+sum(a(1:2,l)*iv2(1:2))
          end do
c         Is iv3 in the main cluster?
          rv2=matmul(Ai,iv3)
          if(rv2(1).lt.oner-eps.and.rv2(1).gt.zeror-eps.and.
     &       rv2(2).lt.oner-eps.and.rv2(2).gt.zeror-eps) then
c            The point is located within the cluster.  
c            Now we must find which point it is.
             do ic=1,Nc
               if(sum(abs(iv3(1:2)-Rc(1:2,ic))).eq.0) then
                 icrdiff(i,j)=ic
                 itest=1
               end if
             end do
             if(iprint.ge.5.and.myrank.eq.0) 
     &         write(lud,"(1x,i2,2x,i2,7x,i2)") i,j,icrdiff(i,j) !0
          end if
        end do
        end do
        if(itest.ne.1) then
          write(lud,*) 'no icrdiff(i,j) found',i,j
          stop
        end if
      end do
      end do  
        

c***********************************************************************        
c     Now create the tables for the neighbors of cluster points
c***********************************************************************        
c     Note we use c1 from above, don't overwrite it!!  
c
c     first locate  the neighbors

      neighbor(1,:) = (/ 1, 0, 0/) ! neighbor in the +x direction
      neighbor(2,:) = (/ 0, 1, 0/) ! neighbor in the +y direction
     
      if(iprint.ge.5.and.myrank.eq.0) then
        write(lud,*) '  '                      !0
        write(lud,*) ' nearest neighbors of ic '       !0
        write(lud,*) ' ic  neigh  nneighbor(neigh,ic)' !0
      end if
      do i=1,Nc
      do j=1,ndim
        iv1(:)=Rc(:,i)+ neighbor(j,:)
        itest=0          ! We dont know where iv1 is.
        do n=-1,1        ! Look in the surrounding 
        do m=-1,1        ! 9 clusters for iv1.
          iv2=(/n,m,0/)
          do l=1,2
            iv3(l)=iv1(l)+sum(a(1:2,l)*iv2(1:2))
          end do
c         Is iv3 in the main cluster?
          rv2=matmul(Ai,iv3)
          if(rv2(1).lt.oner-eps.and.rv2(1).gt.zeror-eps.and.
     &       rv2(2).lt.oner-eps.and.rv2(2).gt.zeror-eps) then
c            The point is located within the cluster.  
c            Now we must find which point it is.
             do ic=1,Nc
               if(sum(abs(iv3(1:2)-Rc(1:2,ic))).eq.0) then
                 nneighbor(j,i)=ic
                 itest=1
               end if
             end do
             if(iprint.ge.5.and.myrank.eq.0) 
     &         write(lud,"(1x,i2,2x,i2,7x,i2)") i,j,nneighbor(j,i) !0
          end if
        end do
        end do
        if(itest.ne.1) then
          write(lud,*) 'no nneighbor(j,i) found',j,i
          stop
        end if
      end do
      end do  

c
c     Now count the number of neighbors in each shell  
c


      if(iprint.ge.5.and.myrank.eq.0) then
        write(lud,*) '  ' !0
        write(lud,*) ' neighbors of icR0 ' !0
        write(lud,*) '  ishell  nneigh(ishell)  nneighl(ishell)' !0
      end if
      i1 = 2*(Nc+1)**(oner/ndim); rv1=zeror; nneigh=0; ishellmax=0
      nneighl=0
      
      do i=-i1,i1
      do j=-i1,i1
        iv1=(/ i, j, 0/)
        i2=sum(abs(iv1)) ! the lattice neighbor index (according to Betts)
        if(i2.gt.0) then
          if(i2.le.nshellmax) nneighl(i2)=nneighl(i2)+1
          do n=1,ndim
            rv1(n)=iv1(n) + halfr*sum(a(:,n)) 
          end do
          rv2=matmul(Ai,rv1)
          if(rv2(1).lt.oner-eps.and.rv2(1).gt.zeror-eps.and.
     &       rv2(2).lt.oner-eps.and.rv2(2).gt.zeror-eps) then
c            This point is within the cluster. Now we need to find
c            the shortest distance between this point and icR0, and
c            the corresponding cluster neighbor index i3.
             i3=i2
             do n=-1,1
             do m=-1,1
               iv2(:) = iv1 + n*a(1,:) + m*a(2,:)
               if(sum(abs(iv2)).lt.i2) i3=sum(abs(iv2))
             end do !m
             end do !n
             nneigh(i3)=nneigh(i3)+1
             if(i3.gt.ishellmax) ishellmax=i3
          end if
        end if
      end do !j
      end do !i
      
      if(sum(nneigh)+1.ne.Nc) then
        if(myrank.eq.0) write(lud,*) 'error in nneigh'
        stop
      end if
        
      if(iprint.ge.5.and.myrank.eq.0) then
        do i1=1,ishellmax
          write(lud,"(2x,i4,8x,i4,12x,i4)") i1,nneigh(i1),nneighl(i1)
        end do
        write(lud,"(' biggest cluster shell=',i3)") ishellmax
      end if
c
c     Now determine the perfection.  It seems that Betts defines the
c     ferromagnetic imperfection I_F as the numer of additional spins in
c     the shells beyond j plus the number missing the shells before j,
c     and not counting those in the jth shell.  The integer j seems to
c     be the number which yields the smallest I_F

      do j=1,ishellmax
        i3=0
        do i1=1,j-1
          i3=i3+abs(nneigh(i1)-nneighl(i1))
        end do
        do i1=j+1,ishellmax
          i3=i3+nneigh(i1)
        end do
        if(j.gt.1) then
          if(i3.lt.I_F) I_F=i3 
        else
          i_F=i3
        end if
      end do
      
      if(myrank.eq.0) 
     &  write(lud,"(' ferromagnetic imperfection I_F= ',i3)") I_F
      
c***********************************************************************        
c       Now set up the K-space tables
c***********************************************************************        

c    The reciprocal lattice vectors g1=(g1_x,g1_y) and g2=(g2_x,g2_y) for 
c    a general set of real space lattice vectors a1=(a1_x,a1_y) and 
c    a2=(a2_x,a2_y) is:
c    
c    g1_x =  a2_y * twopi / area
c    g1_y = -a2_x * twopi / area
c    
c    g2_x = -a1_y * twopi / area
c    g2_y =  a1_x * twopi / area
c    
c    with
c    
c    area = |ax*by-ay*bx|

      r1=a(1,1)*a(2,2) - a(2,1)*a(1,2)
      g(1,1) =+2*pi*a(2,2)/r1
      g(1,2) =-2*pi*a(2,1)/r1
      g(2,1) =-2*pi*a(1,2)/r1
      g(2,2) =+2*pi*a(1,1)/r1
      
      if(iprint.ge.5.and.myrank.eq.0) then
       write(lud,*) '  '
       write(lud,"('g=',2(' (',f6.3,1x,f6.3,')  '))") g(1,1:2),g(2,1:2)
      end if



c***********************************************************************        
c       Create the cluster momentum tables Kc(i,ic), i=x,y
c***********************************************************************        

c     First find the points in the irreducible wedge.  They will be
c     indexed from 1 to Ncw
      if(iprint.ge.5.and.myrank.eq.0) then
        write(lud,*) ' Wedge K points '    !0
        write(lud,*) 'ic  i  j    Kcx(ic)    Kcy(ic)' !0
      end if
      ic=0; Ncw=0; icK0=0; icKpi=0
      do i=-Nc,Nc
      do j=-Nc,Nc
        iv1=(/i,j,0/)    ! i in the g1 and j in the g2 direction
        do l=1,2         ! x,y components 
          rv1(l)=sum(iv1(1:2)*g(1:2,l))
        end do
        if(rv1(1).gt.-pi+eps.and.rv1(1).lt.pi+eps.and.
     &     rv1(2).gt.-pi+eps.and.rv1(2).lt.pi+eps) then
c          The point is in the first SC Brillouin zone.  Now see if the 
c          point falls within the IRW.  First see if it is equivalent to 
c          any known point in the IRW (if so, it cannot be in IRW).
           do igroup=2,ngroup ! exclude the identity
              rv2=transform2D(rv1,group(igroup),ixy,isx,isy)
              itest=0
              do jc=1,Ncw
                 if(sum(abs(Kc(1:2,jc)-rv2(1:2))).lt.eps) then ! rv2 and Kc(ic)
                    itest=1                                    ! are the same 
                    exit                                       ! leave the jc loop
                 end if
              end do
              if(itest.eq.1) exit !leave the igroup loop
           end do
           if(itest.eq.0) then ! rv1 is in 1BZ and not eqivalent to any point in 
              Ncw=Ncw+1        ! the known IRW.  It must be in the IRW
              ic=ic+1
              Kc(:,ic)=rv1(:)
              ig(ic)=i
              jg(ic)=j
              if(sum(abs(rv1(1:2)-pi)).lt.eps) icKpi=ic
              if(sum(abs(rv1(1:2)-zeror)).lt.eps) icK0=ic
              if(iprint.ge.5.and.myrank.eq.0) 
     &          write(lud,5) ic,i,j,Kc(1:2,ic)     !0
           end if
        end if 
      end do 
      end do 
      
      if(icK0.eq.0) then  ! we failed to find (0,0)
        if(myrank.eq.0) 
     &    write(lud,*) 'failed to find (0,0) in the IRW'    !0
        stop
      else if(icKpi.eq.0) then  ! we failed to find (pi,pi)
        icKpi=icK0
        if(myrank.eq.0) then
          write(lud,*) '*****************************************'    !0
          write(lud,*) '* failed to find (pi,pi) in the IRW     *'    !0
          write(lud,*) '* apparently the lattice isnt bipartide *'    !0
          write(lud,*) '* icKpi=icK0 so that the code will run  *'    !0
          write(lud,*) '*****************************************'    !0
        end if
      end if



c     Now find the points in the 1BZ but outside the irreducible wedge.  
c     They will be indexed from Ncw+1 to Nc
      if(iprint.ge.5.and.myrank.eq.0) then
        write(lud,*) ' Non Wedge K points '    !0
      end if
      do i=-Nc,Nc              
      do j=-Nc,Nc
         iv1=(/i,j,0/)    ! i in the g1 and j in the g2 direction
         do l=1,2         ! x,y components 
            rv1(l)=sum(iv1(1:2)*g(1:2,l))
         end do
         if(rv1(1).gt.-pi+eps.and.rv1(1).lt.pi+eps.and.
     &      rv1(2).gt.-pi+eps.and.rv1(2).lt.pi+eps) then
c           The point is in the first SC Brillouin zone.  Now see if the 
c           point falls within the IRW.  
            itest=0
            do jc=1,Ncw
               if(sum(abs(Kc(1:2,jc)-rv1(1:2))).lt.eps) then ! rv2 and Kc(ic) 
                  itest=1                                    ! are the same 
                  exit                                       ! leave the jc loop
               end if
            end do
            if(itest.eq.0) then ! rv1 is in 1BZ and not the IRW
               ic=ic+1
               Kc(:,ic)=rv1(:)
               ig(ic)=i
               jg(ic)=j
               if(iprint.ge.5.and.myrank.eq.0) 
     &              write(lud,5) ic,i,j,Kc(1:2,ic)     !0
            end if
         end if 
      end do 
      end do 
      
      
      if(iprint.ge.5.and.myrank.eq.0) 
     &  write(lud,"('icK0=',i2,' icKpi=',i2,' Ncw=',i2,' Nc=',i2)") 
     &  icK0,icKpi,Ncw,ic

 5    format(1x,i2,1x,i2,1x,i2,2x,f9.6,2x,f9.6)
      ntw=nl*Ncw      


c***********************************************************************
c     Create the indirect addressing table ict(i,j,k) which maps a 
c     point indexed by the principle translations in k-space to ic.
c***********************************************************************        
      N_sp=sqrt(float(Nc))
      do i=-N_sp,N_sp
      do j=-N_sp,N_sp
         iv1=(/i,j,0/)
         do l=1,2         !x,y components 
            rv1(l)=sum(iv1(1:2)*g(1:2,l)) ! K point in xyz
         end do
         do it=1,3
           do l=1,2
              if(rv1(l).lt.-pi+eps) rv1(l)=rv1(l)+twor*pi
              if(rv1(l).gt.pi+eps)  rv1(l)=rv1(l)-twor*pi
           end do
         end do
         itest=0
         do jc=1,Nc
            if(sum(abs(rv1(1:2)-Kc(1:2,jc))).lt.eps) then        
               ict(i,j,0)=jc
               itest=1
            end if
         end do
         if(itest.eq.0) then
            write(lud,5) 'ict failed for i,j=',i,j
            stop
         end if
           
      end do
      end do

c***********************************************************************        
c     Create the table ickdiff(ic1,ic2), which takes the difference
c     between two K-points
c***********************************************************************        
      if(iprint.ge.5.and.myrank.eq.0) then
        write(lud,*) '   '                           !0
        write(lud,*) 'ic1   ic2  ickdiff(ic1,ic2)  ' !0
      end if
      do ic1=1,Nc
      do ic2=1,Nc
         rv1=Kc(:,ic1)-Kc(:,ic2)
         do it=1,3           ! move the difference into the 1BZ
           do l=1,2
              if(rv1(l).lt.-pi+eps) rv1(l)=rv1(l)+twor*pi
              if(rv1(l).gt.pi+eps)  rv1(l)=rv1(l)-twor*pi
           end do
         end do
c        Now we need to figure out where this point is!
         itest=0
         do ic=1,Nc
            if(sqrt(sum((rv1(1:2)-Kc(1:2,ic))**2)).lt.eps) then    
               ickdiff(ic1,ic2)=ic
               itest=1
            end if
         end do
         if(itest.eq.0) then
            write(lud,5) 'ickdiff failed for ic1,ic2=',ic1,ic2
            stop
         end if
         if(iprint.ge.5.and.myrank.eq.0) 
     &        write(lud,"(1x,i2,4x,i2,6x,i2)") ic1,ic2,ickdiff(ic1,ic2) !0
      end do
      end do
          

c***********************************************************************        
c     Create the table ickplus(ic1,ic2), which takes the sum
c     between to K-points
c***********************************************************************        
      if(iprint.ge.5.and.myrank.eq.0) then
        write(lud,*) '   '                           !0
        write(lud,*) 'ic1   ic2  ickplus(ic1,ic2)  ' !0
      end if
      do ic1=1,Nc
      do ic2=1,Nc
         rv1=Kc(:,ic1)+Kc(:,ic2)
         do it=1,3           ! move the difference into the 1BZ
           do l=1,2
              if(rv1(l).lt.-pi+eps) rv1(l)=rv1(l)+twor*pi
              if(rv1(l).gt.pi+eps)  rv1(l)=rv1(l)-twor*pi
           end do
         end do
c        Now we need to figure out where this point is!
         itest=0
         do ic=1,nc
            if(sqrt(sum((rv1(1:2)-Kc(1:2,ic))**2)).lt.eps) then    
               ickplus(ic1,ic2)=ic
               itest=1
            end if
         end do
         if(itest.eq.0) then
            write(lud,5) 'ickplus failed for ic1,ic2=',ic1,ic2
            stop
         end if
         if(iprint.ge.5.and.myrank.eq.0) 
     &     write(lud,"(1x,i2,4x,i2,6x,i2)") ic1,ic2,ickplus(ic1,ic2) !0
      end do
      end do

c***********************************************************************
c     Now find the points equivalent to ic, ickequ(ic,igroup) where 
c     igroup is the group operation.
c***********************************************************************        
c
c
      if(iprint.ge.5.and.myrank.eq.0) then
        write(lud,*) '   '                            !0
        write(lud,*) ' equivalency in k'              !0
        write(lud,*) ' ic  igroup  ickequ(ic,igroup)' !0
      end if
      do ic=1,Nc
         i=ig(ic)
         j=jg(ic)
         iv1=(/i,j,0/)
         do l=1,2
            rv2(l)=sum(iv1(1:2)*g(1:2,l))
         end do
c     Kc(ic) falls within the zone, by construction, now look 
c     for equivalent points
c     
         do igroup=1,ngroup
            rv1=transform2D(rv2,group(igroup),ixy,isx,isy) 
c           These are the equivalent points.  Now map to 1st BZ
            do it=1,3           ! move the difference into the 1BZ
              do l=1,2
                 if(rv1(l).lt.-pi+eps) rv1(l)=rv1(l)+twor*pi
                 if(rv1(l).gt.+pi+eps) rv1(l)=rv1(l)-twor*pi
              end do
            end do
c           Now we need to figure out where this point is!
            itest=0
            do jc=1,Nc
               if(sqrt(sum((rv1(1:2)-Kc(1:2,jc))**2)).lt.eps) then
                  ickequ(ic,igroup)=jc
                  itest=1
               end if
            end do
            if(iprint.ge.5.and.myrank.eq.0) 
     &        write(lud,"(2x,i2,4x,i2,6x,i2)") 
     &        ic,igroup,ickequ(ic,igroup) !0
            if(itest.eq.0) then
               if(myrank.eq.0) write(lud,*)'no equivalent K-point found' !0
c               if(myrank.eq.0) write(lud,*)'j1,n2,m2',j1,n2,m2 !0
               if(myrank.eq.0) write(lud,*)'K=',rv2(1:2) !0
               stop
            end if
         end do
      end do
        
c***********************************************************************        
c     Now generate a table which maps any K to an equivalent point
c     in the irreducible wedge (IW), and a table for the degeneracy 
c     of each point in the IW. 
c 
c     The use of the wedge can greatly reduce the storage required.
c     for example for Nc=8, only 4 of the 8 points are in the irreducible 
c     wedge.  Thus, all storage arrays may be reduced in size by a 
c     factor of 2.  As Nc-->oo, the saving factor approaches 8 for
c     square and diamond clusters.  
c
c***********************************************************************        
c       
      
      ickdeg=0; irwgroup=0
      if(iprint.ge.5.and.myrank.eq.0) then
        write(lud,*) '   '                            !0
        write(lud,*) ' degeneracy in k'              !0
        write(lud,*) ' ic  ickmap(ic)  ickdeg(ic) irwgroup(ic)' !0
      end if
c     scan over all of the K points in the 1BZ
      do ic=1,Nc
         i=ig(ic)
         j=jg(ic)
         iv1=(/i,j,0/)
         do l=1,2
            rv2(l)=sum(iv1(1:2)*g(1:2,l))
         end do
c        Scan over all allowed point group operations, to find the
c        one that maps ic into the IRW
         itest=0
         groupops: do igroup=1,ngroup
           rv1=transform2D(rv2,group(igroup),ixy,isx,isy) 
c          These are the equivalent points.  Now map to 1st BZ
           do it=1,3              ! move the difference into the 1BZ
             do l=1,2
                if(rv1(l).lt.-pi+eps) rv1(l)=rv1(l)+twor*pi
                if(rv1(l).gt.+pi+eps) rv1(l)=rv1(l)-twor*pi
             end do
           end do
c          See if this point is in the IRW!
           do jc=1,Ncw
              if(sqrt(sum((rv1(1:2)-Kc(1:2,jc))**2)).le.eps) then
                 ickmap(ic)=jc
                 ickdeg(jc)=ickdeg(jc)+1
                 irwgroup(ic)=igroup
                 itest=1
              end if
           end do
           if(itest.eq.1) exit    ! at first group op. that takes ic into IRW
         end do groupops
         if(itest.eq.0) then
            if(myrank.eq.0) write(lud,*) 'mapping failed for' !0
            if(myrank.eq.0) write(lud,*) 'ic',ic !0
            if(myrank.eq.0) write(lud,*) 'K=',rv2(1:2) !0
            stop
         end if
      end do !ic
      
      if(iprint.ge.5.and.myrank.eq.0) then 
        do ic=1,Nc    
          write(lud,"(2x,i2,6x,i2,9x,i2,9x,i2)") 
     &    ic,ickmap(ic),ickdeg(ic),irwgroup(ic)
        end do
      end if
         


c     Test the degeneracy table.  Since ickdeg holds the degeneracy
c     of points in the IRW, its sum must equal Nc.
      i=0
      do ic=1,Ncw
         i=i+ickdeg(ic)
      end do
      if(i.ne.Nc) then
         if(myrank.eq.0) write(lud,*) 'bug in ickdeg',i
         stop
      end if
        
c     Test the group table.  irwgroup(ic within wedge)=1, the identity
      do ic=1,Ncw
      if(irwgroup(ic).ne.1) then
        if(myrank.eq.0) write(lud,*) 'irwgroup failed for wedge pt.' !0
        if(myrank.eq.0) write(lud,*) 'ic',ic,'irwgroup(ic)',irwgroup(ic) !0
        if(myrank.eq.0) write(lud,"('K=(',f6.3,1x,f6.3,')')") Kc(1:2,ic) !0
      end if
      end do


c***********************************************************************
c     Form Epsbar(K)
c***********************************************************************        

      if(myrank.eq.0.and.iprint.ge.5) write(lud,*) '    '
      if(myrank.eq.0.and.iprint.ge.5) write(lud,*) 
     &               'ic   Epsbar(ic)   Dwavebar(ic)'

c first generate the k-points in the Brillouin-zone of the 
c super-lattice, i.e. in the Wigner-Seitz cell of the reciprocal 
c space defined by the cluster K-points
      
      
      Dkt=oner/real(ntot,kind)
      nwsc=0 ! number of k-points in Wigner-Seitz cell 
      if (ntot.gt.0) then
      do i=-ntot+1,ntot
        do j=1,ntot
          iv1=(/i,j,0/)
          rv1(1)=Dkt*(iv1(1)-halfr)*pi
          rv1(2)=Dkt*(iv1(2)-halfr)*pi
c  Check if this k-point is in Wigner-Seitz cell, i.e. if the closest
c  K-point is K=0 
            r1=rv1(1)*rv1(1)+rv1(2)*rv1(2) ! distance from K=0
            itest=0
            do jc=1,Nc
            if (jc.eq.icK0) cycle
             r2=rv1(1)-Kc(1,jc)
             r3=rv1(2)-Kc(2,jc)
             if((r2-pi).gt.zeror) then; r2=r2-twor*pi; endif
             if((r2+pi).lt.zeror) then; r2=r2+twor*pi; endif
             if((r3-pi).gt.zeror) then; r3=r3-twor*pi; endif
             if((r3+pi).lt.zeror) then; r3=r3+twor*pi; endif
             r2=r2*r2+r3*r3
             if(r2.lt.r1) then
               itest=1
               exit
              endif
            enddo
            if(itest.eq.0) then 
             nwsc=nwsc+1
             kt(1,nwsc) = rv1(1)
             kt(2,nwsc) = rv1(2)
             nwsc=nwsc+1 ! -k must also be part of the WS cell
             kt(1,nwsc) = -rv1(1)
             kt(2,nwsc) = -rv1(2)
             endif
        enddo ! j
       enddo ! i
       else ! needed for FSS calculation
          nwsc=1
          kt(:,1)=zeror
       endif
c
       if(myrank.eq.0.and.iprint.ge.5) then
          write(lud,*) 'Number of points in WS*Nc', nwsc*Nc
          write(lud,*) 'Total number of points in BZ', (2*ntot)**2 
       endif 

      Epsbar=zeror
      Dwavebar=zeror
      do ic=1,Nc
         do i=1,nwsc
            rv1(1) = Kc(1,ic)+kt(1,i)
            rv1(2) = Kc(2,ic)+kt(2,i)
c           Epsbar(ic)=Epsbar(ic) - halfr*sum(cos(rv1(1:2)))
          Epsbar(ic)= Epsbar(ic) -halfr*(cos(rv1(1))+cos(rv1(2))) 
     &                - tprime*(cos(rv1(1))*cos(rv1(2))-oner)
     &              - tdprime*halfr*(cos(twor*rv1(1))+cos(twor*rv1(2))) 
     &           - t3prime*(cos(twor*rv1(1))*cos(rv1(2))
     &                         +cos(rv1(1))*cos(twor*rv1(2))) 


c            Epsbar(ic)= Epsbar(ic) -halfr*(cos(rv1(1))+cos(rv1(2))) 
c     &                - tprime*(cos(rv1(1))*cos(rv1(2))-oner)
            Dwavebar(ic)=Dwavebar(ic) -halfr*(cos(rv1(1))-cos(rv1(2))) 
         end do
         Epsbar(ic)=Epsbar(ic)/real(nwsc,kind)
         Dwavebar(ic)=Dwavebar(ic)/real(nwsc,kind)
      if(myrank.eq.0.and.iprint.ge.5)
     &            write(lud,"(i2,6x,f9.6,2x,f9.6,2x,f9.6,2x,f9.6)") 
     &           ic,Epsbar(ic),Dwavebar(ic),cos(Kc(1,ic))-cos(Kc(2,ic)),
     &           Dwavebar(ic)/(cos(Kc(1,ic))-cos(Kc(2,ic)))
      end do
c     Check 
      r1=sum(Epsbar); r2=sum(Epsbar*ickdeg)
      r3=sum(Dwavebar)
      if(abs(tprime).lt.eps.and.abs(r1).gt.eps
     &   .or.
     &	 abs(tprime).lt.eps.and.abs(r2).gt.eps) then
        if(myrank.eq.0) then
          write(lud,*) 'sum Epsbar wrong.',r1,r2
	end if
c          stop
       endif
       if(abs(r3).gt.eps) then
        if(myrank.eq.0) then
          write(lud,*) 'sum Dwavebar wrong.',r3
	end if
c        stop
      end if

     
c***********************************************************************        
c     Form the Fourier transform tables.  
c     R-->K,  FTCoefs_K_to_R
c     K-->R,  FTCoefs_K_to_R
c
c                  __          the factor of 1/N is required
c               1  \           so that both G(R)~1/w and G(K)~1/w
c     G(R=0) = --- /  G(K)     for large w
c              N   --
c                  K
c***********************************************************************        
      do ic=1,Nc
        do ik=1,Nc
          r1=zeror
          do idim=1,ndim
             r1=r1+Kc(idim,ik)*Rc(idim,ic)
          end do
          FTCoefs_K_to_R(ik,ic)=exp(-ii*r1)/real(Nc,kind)
          FTCoefs_R_to_K(ik,ic)=exp(+ii*r1)
        end do
      end do
      
c     Now test for consistency of Kc and Rc
      do ic=1,Nc
        r1= real(sum(FTCoefs_K_to_R(1:Nc,ic)))-kroniker(ic,icR0)
        r2=aimag(sum(FTCoefs_K_to_R(1:Nc,ic)))
        r3= real(sum(FTCoefs_K_to_R(ic,1:Nc)))-kroniker(ic,icK0)
        r4=aimag(sum(FTCoefs_K_to_R(ic,1:Nc)))
        
        if(abs(r1).gt.eps.or.abs(r2).gt.eps) then
          write(6,*) 'Kc and Rc are inconsistent'
          write(6,*) 'icr,r1,r2=',ic,r1,r2
          stop
        end if
        if(abs(r3).gt.eps.or.abs(r4).gt.eps) then
          write(6,*) 'Kc and Rc are inconsistent, K=0'
          write(6,*) 'ick,r3,r4=',ic,r3,r4
          stop
        end if
      end do

c***********************************************************************        
c       Create the d-wave form factor
c***********************************************************************        
        gdwave=zeror
        gB1g=zeror
	gB2g=zeror
        do ik=1,Nc
          do idim=1,2    ! d_x2-y2 and B1g form factor
            gdwave(ik)=gdwave(ik)+cos(Kc(idim,ik))*(-1)**(idim-1)
          end do
c	  gB2g(ik)=sin(Kc(1,ik))*sin(Kc(2,ik))  ! Incorrect B2g form factor
        end do
c
c Add proper form factor for Raman B1g and B2g that includes the proper
c Q dependence
c
        do ic=1,Nc
          do ik=1,Nc
            gB1g(ik,ic) = cos(Kc(1,ik)+Kc(1,ic)/twor)
     &                   -cos(Kc(2,ik)+Kc(2,ic)/twor)
            gB2g(ik,ic) = sin(Kc(1,ik)+Kc(1,ic)/twor)*
     &                    sin(Kc(2,ik)+Kc(2,ic)/twor)
          end do
        end do

c print out the d-wave order parameter in real space
      if(myrank.eq.0.and.iprint.ge.5) write(lud,*) '     '
      if(myrank.eq.0.and.iprint.ge.5) 
     &        write(lud,*) 'Dwave order parameter'
      if(myrank.eq.0.and.iprint.ge.5) 
     &        write(lud,*) 'ic, jc, Dwave(ic,jc)'
      Dwave(1:32,1:32) = cmplx(zeror,zeror)
      do ic =1,Nc   ! for all R1
       do jc=1,Nc   ! for all R2
        do ik=1,Nc  ! sum over K
          Dwave(ic,jc)=Dwave(ic,jc)+FTCoefs_K_to_R(ik,ic)
     &               *conjg(FTCoefs_K_to_R(ik,jc))*Dwavebar(ik)
        end do
        if (abs(Dwave(ic,jc)).gt.eps) then
         if(myrank.eq.0.and.iprint.ge.5) write(lud,"(i2,3x,i2,3x,f9.6)")
     &      ic,jc,real(Dwave(ic,jc)*Nc)
        end if
      end do
      end do


      return
      end  
      
      
c***********************************************************************        
c***********************************************************************        
c         2D subroutines and functions
c***********************************************************************        
c***********************************************************************        

        
      function gtransform2D(a,igroup,ixy,isx,isy)
c
c  x <-- isx*[(1-ixy)*x + ixy*y] 
c  y <-- isy*[(1-ixy)*y + ixy*x] 
c
      implicit none
      integer :: i,a(3,3),g(3,3),gtransform2D(3,3),igroup,ixy(8),
     &           isx(8),isy(8)
      
      do i=1,3  
         gtransform2D(i,1)=isx(igroup)*((1-ixy(igroup))*a(i,1) 
     &                  + ixy(igroup)*a(i,2))
         gtransform2D(i,2)=isy(igroup)*((1-ixy(igroup))*a(i,2) 
     &                  + ixy(igroup)*a(i,1))
      end do
      
      return
      end
        
      function Rtransform2D(Rc,igroup,ixy,isx,isy)
c
c  x <-- isx*[(1-ixy)*x + ixy*y] 
c  y <-- isy*[(1-ixy)*y + ixy*x] 
c
      implicit none
      integer :: i,Rc(3),g(3),Rtransform2D(3),igroup,ixy(8),
     &           isx(8),isy(8)
      
      Rtransform2D(1)=isx(igroup)*((1-ixy(igroup))*Rc(1)
     &             + ixy(igroup)*Rc(2))
      Rtransform2D(2)=isy(igroup)*((1-ixy(igroup))*Rc(2)
     &             + ixy(igroup)*Rc(1))
                         
      return
      end
      
      function Ktransform2D(Kc,igroup,ixy,isx,isy)
c     
c  x <-- isx*[(1-ixy)*x + ixy*y] 
c  y <-- isy*[(1-ixy)*y + ixy*x] 
c
      implicit none
      integer :: i,igroup, ixy(8),isx(8),isy(8)
      real(8) :: Kc(3),g(3),Ktransform2D(3)
      
      Ktransform2D(1)=isx(igroup)*((1-ixy(igroup))*Kc(1) 
     &             + ixy(igroup)*Kc(2))
      Ktransform2D(2)=isy(igroup)*((1-ixy(igroup))*Kc(2) 
     &             + ixy(igroup)*Kc(1))
      
      return
      end
      
      function kroniker(i,j)
      integer i,j,kroniker
      
      if(i.eq.j) then
        kroniker=1
      else
        kroniker=0   
      end if
      
      return
      end 









      subroutine tables_3D_Betts      ! THREE DIMENSIONAL BETTS TABLES
c***************************************************************************
c     CHANGES REQUIRED FOR USE IN QMC CODES:
c     1. remove comment(s) in front of
c          use Global
c          define the  fmat matrix
c          call allocate_2p
c     2. change program and end statements to subroutine and return
c     3. delete the "Stuff that belongs in module_Global"
c     4. remove "write(6,*) 'enter cluster'"...
c     5. remove any duplicate defintions of kroniker
c     6. remove these comments
c     CHANGES REQUIRED FOR USE IN Analysis CODES:
c     1. do the steps above
c     2. remove the delcarations of ict, ig, jg and kg
c***************************************************************************
c     This subroutine makes the lookup tables used throughout the
c     full THREE DIMENSIONAL code.  It employs the definitions and 
c     cluster geometries developed by D.D. Betts and G.E. Stewart,
c     Canadian Journal of Physics, 75, 47-66.
c
c     Developers include M. Jarrell (symmetry tables), P. Kent (fast 
c     determination of allowed group operations, large clusters nc>26)
c     Th. Maier (WS cells of the cluster superlattice),
c     Th. Pruschke (some initial discussion).  
c     We reserve all rights to this code.
c***************************************************************************
      use Global
      implicit none
c***************************************************************************
      interface transform3D
      function gtransform3D(a,igroup,ixy,ixz,iyz,isx,isy,isz)
      integer :: a(3,3),gtransform3D(3,3),igroup,ixy(48),ixz(48),
     &           iyz(48),isx(48),isy(48),isz(48)
      end
      
      function Rtransform3D(Rc,igroup,ixy,ixz,iyz,isx,isy,isz)
      integer :: Rc(3),Rtransform3D(3),igroup,ixy(48),ixz(48),iyz(48),
     &     isx(48),isy(48),isz(48)
      end
      
      function Ktransform3D(Kc,igroup,ixy,ixz,iyz,isx,isy,isz)
      integer    :: igroup, ixy(48),ixz(48),iyz(48),
     &     isx(48),isy(48),isz(48)
      real(8) :: Kc(3),Ktransform3D(3)
      end
      end interface 
      
      interface Inverse
      function Inverse_3x3(a)
      integer :: a(3,3)
      real(8) :: Inverse_3x3(3,3)
      end function
      end interface
      
      interface determinant
      function Ideterminant(a)
      integer :: a(3,3),Ideterminant
      end function
      function Rdeterminant(a)
      real(8) :: a(3,3),Rdeterminant
      end function
      end interface


c***************************************************************************
c     Locally defined
c***************************************************************************
      logical    :: allowed,denied,denied_lin_alg
      integer, parameter :: nshellmax=50, maxnc=80
      integer :: i,i1,i2,i3,ix,iy,it,igroup,ic,ic1,ic2,ick,idim,ind,
     &           ik,itest,j,jmd,jm,j1,j2,jc,jk,jkc,k,k2,kroniker,
     &           l,m,m1,m2,n,n1,n2,p,Lc,iv1(3),iv2(3),iv3(3),
     &           ig(maxnc),jg(maxnc),kg(maxnc),
     &           a1(3,3),a2(3,3),
     &           ixy(48),ixz(48),iyz(48),isx(48),isy(48),isz(48),
     &           neighbor(3,3),
     &           irwgroup(maxnc),group(48),nneigh(nshellmax),I_F,
     &           nneighl(nshellmax),nneighp(nshellmax),ishellmax
      real(kind) :: r1,r2,r3,r4,r5,r6,r7,r8,r9,Dkt,c1(3,3),c2(3,3),
     &              rv1(3),rv2(3),rv3(3),
     &              corner(3,8),Bettsd(4),Bettsf(6),Bettsl(3),
     &              meand,meanf,meanl,cubicity
c For lin algebra test
      real(kind) :: ainv(3,3),a1ainv(3,3)
c Select use of "brute force" vs linear algebra method for 
c checking lattice coincidence. Both methods should agree; lin alg
c method is much faster.
      logical, parameter :: use_brute_force = .false.
      logical, parameter :: use_lin_alg = .true.

c***************************************************************************
c       Local arrays
c       ixy(igroup)      D=3 cubic point group table(see below)
c       ixz(igroup)      D=3 cubic point group table(see below)
c       iyz(igroup)      D=3 cubic point group table(see below)
c       isx(igroup)      D=3 cubic point group table(see below)
c       isy(igroup)      D=3 cubic point group table(see below)
c       isz(igroup)      D=3 cubic point group table(see below)
c       group(i)         Indices of allowed group operations
c       irwgroup(ic)     the group operation maps K(ic) into IRW
c       ig(ic)           i in the g1 and j in the g2 direction
c       jg(jc)           k in the g3 direction
c       kg(ic)           corresponding to K point ic
c       ict(i,j,k)       maps i*g1+j*g2+k*g3 to ic
c       nneigh(i)        number of neighbors on shell i within the cluster
c       nneighl(i)       number of neighbors on shell i on the lattice
c       nneighp(i)       number of neighbors on shell i on a perfect cluster
c***************************************************************************
        
c       CLUSTER MOMENTUM and DUAL SPACE for the 3D cubic lattice
c                                                             _   _   _
c       Different Real-space Cluster Geometries.  The vectors a1, a2, a3
c       define the parallelpipeds which define the periodically replicated 
c       clusters.  Here, we only use the clusters which Betts considers "good".
c               _             _            _
c  cluster      a1            a2           a3       group  Ng  I_F  bipartide 
c    1A     ( 1, 0, 0)    ( 0, 1, 0)   ( 0, 0, 1)     Oh   48   0     no
c
c    8A     ( 2, 0, 0)    ( 0, 2, 0)   ( 0, 0, 2)     Oh   48   4     yes
c    8B     ( 2, 1, 0)    ( 0, 1, 2)   ( 1,-1, 1)     C2h  4    0     no  (=C2 x Ci)
c
c   10A     ( 1, 0, 2)    ( 1, 2,-1)   ( 1,-2, 0)     Ci   2    0     no
c
c   12A     ( 2, 1, 1)    ( 1,-2, 0)   ( 1, 0,-2)     C2h  4    0     no
c   12B     ( 0, 2, 0)    ( 0, 0, 3)   ( 2, 0, 1)     C2h  4    2     no
c
c   14A     ( 2, 2, 2)    ( 2,-1, 1)   ( 1, 2,-1)     C3i  6    1    yes
c   14B     ( 2, 2, 1)    ( 0,-1, 2)   ( 2,-1, 0)     Ci   2    0     no
c   14C     ( 2, 2, 0)    (-1, 2, 1)   ( 1, 0, 2)     Ci   2    0     no
c
c   16A     ( 2, 2, 0)    ( 2, 0, 2)   ( 0, 2, 2)     Oh  48    2    yes*
c   16B     ( 2, 2, 0)    ( 1,-1,-2)   (-1, 1,-2)     D2h  8    2    yes* (=D2 x Ci)
c   16C     ( 4, 2, 0)    ( 1, 2,-1)   (-1, 2, 1)     D2h  8    2    yes*
c   16D     ( 2,-2, 0)    ( 0, 2,-2)   ( 1, 2, 1)     C2h  4    2    yes
c   16E     ( 2, 2, 0)    ( 2,-2, 0)   ( 1, 0,-2)     C2h  4    0     no
c   16F     ( 1, 2, 0)    (-3, 1,-1)   ( 0, 2, 2)     C2h  4    0     no
c   16Z     ( 2, 0, 0)    ( 0, 2, 0)   ( 0, 0, 4) [Made up BAD cluster]
c
c   18A     ( 2, 2, 2)    ( 2,-1,-1)   ( 1, 1,-2)     D3h 12    3    yes*
c   18B     ( 3, 2, 1)    ( 2,-1,-1)   ( 1, 1,-2)     C3i  6    3    yes*
c   18C     ( 2, 2, 2)    ( 2, 1,-1)   (-1, 2,-1)      Ci  2    3    yes
c   18D     ( 0, 0, 3)    ( 2,-2,-1)   ( 2, 1,-1)     C2h  4    1     no
c   18E     ( 3, 3, 0)    ( 2, 0, 2)   ( 2, 0,-1)     C2h  4    1     no
c   18F     ( 3, 2, 1)    ( 0, 2,-2)   ( 0, 2, 1)     C2h  4    2     no
c
c   20A     ( 3, 2, 1)    ( 2,-1,-1)   (-1, 1,-2)     Ci   2    4     yes
c   20B     ( 4, 1, 1)    ( 0, 2,-2)   ( 1, 2, 1)     Ci   2    4     yes
c   20C     ( 2, 2, 2)    ( 2, 0,-2)   (-1, 1,-2)     Ci   2    4     yes
c   20D     ( 3, 1, 0)    ( 0, 2, 2)   ( 1,-2, 1)     Ci   2    4     yes
c   20E     ( 4, 1, 1)    ( 2,-2, 0)   ( 1,-1,-2)     C2h  4    4     yes
c   20F     ( 2, 2, 1)    ( 2,-2, 1)   ( 1, 0,-2)     C2h  4    0     yes
c   20G     ( 2, 2, 1)    ( 2,-2, 0)   ( 1, 1,-2)     C2h  4    2     no
c
c   22A     ( 4, 2, 0)    ( 0, 3,-1)   (-1, 2, 1)     Ci   2    5     yes
c   22B     ( 2, 2, 2)    ( 3, 0,-1)   (-1, 1,-2)     Ci   2    5     yes
c   22C     ( 3, 2, 0)    ( 0, 2, 2)   ( 1,-1, 2)     Ci   2    4     no
c   22D     ( 1, 2, 0)    ( 2, 1, 2)   ( 2,-2, 1)     Ci   2          no  (BAD)
c   22E     ( 4, 1, 1)    ( 2, 2,-3)   ( 0, 2,-1)     Ci   2    6     no
c
c   24A     ( 3, 1, 0)    ( 0, 2,-2)   ( 0, 2, 2)     C2h  4    6     yes
c   24B     ( 2, 2, 2)    ( 0, 2,-2)   ( 2, 0,-2)     D3h 12    6     yes
c   24C     ( 2, 2, 2)    ( 3, 0,-1)   (-1, 2,-1)     Ci   2    6     yes
c   24D     ( 0, 0, 3)    ( 1, 3, 0)   (-2, 2, 0)     D2h  8    6     no
c   24E     ( 3, 2, 0)    ( 3,-2, 0)   ( 1, 0,-2)     C2h  4    4     no
c                                                           
c   26A     ( 3, 2, 1)    ( 2,-2,-2)   ( 1, 1,-2)     Ci   2    6     yes
c   26B     ( 1, 3, 2)    ( 3,-2, 1)   ( 3, 1, 0)     C3i  6    6     yes
c   26C     ( 2, 1, 3)    ( 2,-1,-2)   ( 2, 2,-1)     Ci   2    2     no
c
c   26Z     ( 1, 2, 3)    ( 3, 3, -2)  (-3,-2,-3)         20          yes
c
c  The bipartide lattices had a eps_0^HA in table 8 (and common sense), 
c  I_F (the feerro imperfection) was calculated below and *not* taken from Betts
c  A * means equivlent to others with same Nc.
c  I.e., lattices 16A=16B=16C, 18A=18B,  20B=20C=20D.
c
c  Larger lattices were calculated with gen_cluster.f Each "A" cluster is the
c  most perfect cluster with the smallest cubicity for a given size. Each "B"
c  cluster is the most perfect cluster for a given size that is also bipartide.
c
c  ----------------------------------------------------------------------------
c                       GROUP OPERATIONS
c  ----------------------------------------------------------------------------
c  Only 16A has all 48 operations and has the point group Oh.  The remaining clusters 
c  have lower symmetry than the cubic lattice, with symmetries designated by the
c  Schoenflies symbols (From Ibach and Luth):
c  Cj (j=2,3,4, 6) j-fold rotation axis 
c  Sj j-fold rotation-inversion axis 
c  Dj j 2-fold rotation axes perpendicular to a j-fold principle rotation axis 
c  T  4 three-and 3 two-fold rotation axes, as in a tetrahedron 
c  O  4 three-and 3 four-fold rotation axes, as in a octahedron 
c  Ci a center of inversion 
c  Cs a mirror plane
c  In addition, their are sufixes for mirror planes 
c     h: horizontal=perpendicular to the c  rotation axis, 
c     v: vertical=parallel to the main rotation axis in the plane, 
c     d: diagonal=parallel to the main rotation axis in the plane  
c        bisecting the two-fold rotation axes.
c 
c  The Clusters identified by Betts, have several different symmetries
c
c  C2h   A two fold axis, plus a single mirror plane that is perpendicular to the axis
c  Ci    Contains only the identity and inversion
c  C3i   Triagonal with a 3-fold axis and inversion (GUESS)??
c          * rotations by 2pi/3 about the axix x=y=z
c          * inversion wrt the origin.
c  Oh    The full symmetry of the cubic lattice 
c  D2h   2 2-fold axes perp. to the 2-fold rotation axis, plus a horizonal mirror plane
c  D3h   3 2-fold axes perp. to the 3-fold rotation axis, plus a horizonal mirror plane
c
c  We need to automate the process of finding which of the 48 point group operations are
c  retained by the different clusters.  To do this, we must perform all 48 cubic group
c  operations on the {a} and see if they form an equivalent parallelpiped.  
c
c  Group Operations for a cubic lattice:
c
c  x --> +/- x                  There are 2^D D! operations in D dimensions
c    \/                         or 48 operations in 3D.
c    /\
c  y --> +/- y   
c    \/
c    /\
c  z --> +/- z   
c        .
c  a'=transform3D(a,igroup,ixy,ixz,iyz,isx,isy,isz)
c  ixy=0(1) dont (do) exchange x and y
c  ixz=0(1) dont (do) exchange x and z
c  iyz=0(1) dont (do) exchange y and z
c  isx      sign of x 
c  isy      sign of y
c  isz      sign of z
c
c  x <-- isx*[(1-ixy)*x + ixy*y]  <-- [(1-ixz)*x + ixz*z]<----------------------
c  y <-- isy*[(1-ixy)*y + ixy*x]  <----------------------<-- [(1-iyz)*y + iyz*z]
c  z <-- isz*<------------------- <-- [(1-ixz)*z + ixz*x]<-- [(1-iyz)*z + iyz*y]
c

      ixy( 1)=0;ixz( 1)=0; iyz( 1)=0; isx( 1)=+1; isy( 1)=+1; isz( 1)=+1 !identity  xyz e
      ixy( 2)=1;ixz( 2)=0; iyz( 2)=0; isx( 2)=+1; isy( 2)=+1; isz( 2)=+1 !ref. x=y  yxz o
      ixy( 3)=0;ixz( 3)=1; iyz( 3)=0; isx( 3)=+1; isy( 3)=+1; isz( 3)=+1 !ref. x=z  zyx o
      ixy( 4)=0;ixz( 4)=0; iyz( 4)=1; isx( 4)=+1; isy( 4)=+1; isz( 4)=+1 !ref. y=z  xzy o
      ixy( 5)=1;ixz( 5)=1; iyz( 5)=0; isx( 5)=+1; isy( 5)=+1; isz( 5)=+1 !          yzx e
      ixy( 6)=1;ixz( 6)=0; iyz( 6)=1; isx( 6)=+1; isy( 6)=+1; isz( 6)=+1 !          zxy e
      
      ixy( 7)=0;ixz( 7)=0; iyz( 7)=0; isx( 7)=-1; isy( 7)=+1; isz( 7)=+1 !
      ixy( 8)=1;ixz( 8)=0; iyz( 8)=0; isx( 8)=-1; isy( 8)=+1; isz( 8)=+1 !
      ixy( 9)=0;ixz( 9)=1; iyz( 9)=0; isx( 9)=-1; isy( 9)=+1; isz( 9)=+1 !
      ixy(10)=0;ixz(10)=0; iyz(10)=1; isx(10)=-1; isy(10)=+1; isz(10)=+1 !
      ixy(11)=1;ixz(11)=1; iyz(11)=0; isx(11)=-1; isy(11)=+1; isz(11)=+1 !
      ixy(12)=1;ixz(12)=0; iyz(12)=1; isx(12)=-1; isy(12)=+1; isz(12)=+1 !
     
      ixy(13)=0;ixz(13)=0; iyz(13)=0; isx(13)=+1; isy(13)=-1; isz(13)=+1 !
      ixy(14)=1;ixz(14)=0; iyz(14)=0; isx(14)=+1; isy(14)=-1; isz(14)=+1 !
      ixy(15)=0;ixz(15)=1; iyz(15)=0; isx(15)=+1; isy(15)=-1; isz(15)=+1 !
      ixy(16)=0;ixz(16)=0; iyz(16)=1; isx(16)=+1; isy(16)=-1; isz(16)=+1 !
      ixy(17)=1;ixz(17)=1; iyz(17)=0; isx(17)=+1; isy(17)=-1; isz(17)=+1 !
      ixy(18)=1;ixz(18)=0; iyz(18)=1; isx(18)=+1; isy(18)=-1; isz(18)=+1 !
      
      ixy(19)=0;ixz(19)=0; iyz(19)=0; isx(19)=-1; isy(19)=-1; isz(19)=+1 !
      ixy(20)=1;ixz(20)=0; iyz(20)=0; isx(20)=-1; isy(20)=-1; isz(20)=+1 !
      ixy(21)=0;ixz(21)=1; iyz(21)=0; isx(21)=-1; isy(21)=-1; isz(21)=+1 !
      ixy(22)=0;ixz(22)=0; iyz(22)=1; isx(22)=-1; isy(22)=-1; isz(22)=+1 !
      ixy(23)=1;ixz(23)=1; iyz(23)=0; isx(23)=-1; isy(23)=-1; isz(23)=+1 !
      ixy(24)=1;ixz(24)=0; iyz(24)=1; isx(24)=-1; isy(24)=-1; isz(24)=+1 !
      
      ixy(25)=0;ixz(25)=0; iyz(25)=0; isx(25)=+1; isy(25)=+1; isz(25)=-1 !
      ixy(26)=1;ixz(26)=0; iyz(26)=0; isx(26)=+1; isy(26)=+1; isz(26)=-1 !
      ixy(27)=0;ixz(27)=1; iyz(27)=0; isx(27)=+1; isy(27)=+1; isz(27)=-1 !
      ixy(28)=0;ixz(28)=0; iyz(28)=1; isx(28)=+1; isy(28)=+1; isz(28)=-1 !
      ixy(29)=1;ixz(29)=1; iyz(29)=0; isx(29)=+1; isy(29)=+1; isz(29)=-1 !
      ixy(30)=1;ixz(30)=0; iyz(30)=1; isx(30)=+1; isy(30)=+1; isz(30)=-1 !
      
      ixy(31)=0;ixz(31)=0; iyz(31)=0; isx(31)=-1; isy(31)=+1; isz(31)=-1 !
      ixy(32)=1;ixz(32)=0; iyz(32)=0; isx(32)=-1; isy(32)=+1; isz(32)=-1 !
      ixy(33)=0;ixz(33)=1; iyz(33)=0; isx(33)=-1; isy(33)=+1; isz(33)=-1 !
      ixy(34)=0;ixz(34)=0; iyz(34)=1; isx(34)=-1; isy(34)=+1; isz(34)=-1 !
      ixy(35)=1;ixz(35)=1; iyz(35)=0; isx(35)=-1; isy(35)=+1; isz(35)=-1 !
      ixy(36)=1;ixz(36)=0; iyz(36)=1; isx(36)=-1; isy(36)=+1; isz(36)=-1 !
      
      ixy(37)=0;ixz(37)=0; iyz(37)=0; isx(37)=+1; isy(37)=-1; isz(37)=-1 !
      ixy(38)=1;ixz(38)=0; iyz(38)=0; isx(38)=+1; isy(38)=-1; isz(38)=-1 !
      ixy(39)=0;ixz(39)=1; iyz(39)=0; isx(39)=+1; isy(39)=-1; isz(39)=-1 !
      ixy(40)=0;ixz(40)=0; iyz(40)=1; isx(40)=+1; isy(40)=-1; isz(40)=-1 !
      ixy(41)=1;ixz(41)=1; iyz(41)=0; isx(41)=+1; isy(41)=-1; isz(41)=-1 !
      ixy(42)=1;ixz(42)=0; iyz(42)=1; isx(42)=+1; isy(42)=-1; isz(42)=-1 !
      
      ixy(43)=0;ixz(43)=0; iyz(43)=0; isx(43)=-1; isy(43)=-1; isz(43)=-1 !
      ixy(44)=1;ixz(44)=0; iyz(44)=0; isx(44)=-1; isy(44)=-1; isz(44)=-1 !
      ixy(45)=0;ixz(45)=1; iyz(45)=0; isx(45)=-1; isy(45)=-1; isz(45)=-1 !
      ixy(46)=0;ixz(46)=0; iyz(46)=1; isx(46)=-1; isy(46)=-1; isz(46)=-1 !
      ixy(47)=1;ixz(47)=1; iyz(47)=0; isx(47)=-1; isy(47)=-1; isz(47)=-1 !
      ixy(48)=1;ixz(48)=0; iyz(48)=1; isx(48)=-1; isy(48)=-1; isz(48)=-1 !

      

c***********************************************************************        
c       determine the principle translation vectors a1,a2, and a3 
c***********************************************************************
      if(iprint.ge.5.and.myrank.eq.0) 
     &     write(lud,*) ' cluster= ',cluster !0

      If(cluster.eq.'1A') then                              
         a(1,:)=(/1, 0, 0/); a(2,:)=(/0, 1, 0/); a(3,:)=(/0, 0, 1/)
         
      else if(cluster.eq.'8A') then                           
         a(1,:)=(/2, 0, 0/); a(2,:)=(/0, 2, 0/); a(3,:)=(/0, 0, 2/)
      else if(cluster.eq.'8B') then                           
         a(1,:)=(/2, 1, 0/); a(2,:)=(/0, 1, 2/); a(3,:)=(/1,-1, 1/)
         
      else if(cluster.eq.'10A') then                           
         a(1,:)=(/1, 0, 2/); a(2,:)=(/1, 2,-1/); a(3,:)=(/1,-2, 0/)
         
      else if(cluster.eq.'12A') then                           
         a(1,:)=(/2, 1, 1/); a(2,:)=(/1,-2, 0/); a(3,:)=(/1, 0,-2/)
      else if(cluster.eq.'12B') then                           
         a(1,:)=(/0, 2, 0/); a(2,:)=(/0, 0, 3/); a(3,:)=(/2, 0, 1/)
         
      else if(cluster.eq.'14A') then                           
         a(1,:)=(/2, 2, 2/); a(2,:)=(/2,-1, 1/); a(3,:)=(/1, 2,-1/)
      else if(cluster.eq.'14B') then                           
         a(1,:)=(/2, 2, 1/); a(2,:)=(/0,-1, 2/); a(3,:)=(/2,-1, 0/)
      else if(cluster.eq.'14C') then                           
         a(1,:)=(/2, 2, 0/); a(2,:)=(/-1, 2, 1/); a(3,:)=(/1, 0, 2/) 
         
      else if(cluster.eq.'16A') then                           
         a(1,:)=(/2, 2, 0/); a(2,:)=(/2, 0, 2/); a(3,:)=(/0, 2, 2/) 
      else if(cluster.eq.'16B') then                           
         a(1,:)=(/2, 2, 0/); a(2,:)=(/1,-1,-2/); a(3,:)=(/-1, 1,-2/) 
      else if(cluster.eq.'16C') then                           
         a(1,:)=(/4, 2, 0/); a(2,:)=(/1, 2,-1/); a(3,:)=(/-1, 2, 1/)
      else if(cluster.eq.'16D') then                           
         a(1,:)=(/2,-2, 0/); a(2,:)=(/0, 2,-2/); a(3,:)=(/ 1, 2, 1/) 
      else if(cluster.eq.'16E') then                           
         a(1,:)=(/2, 2, 0/); a(2,:)=(/2,-2, 0/); a(3,:)=(/ 1, 0,-2/) 
      else if(cluster.eq.'16F') then                           
         a(1,:)=(/1, 2, 0/); a(2,:)=(/-3,1,-1/); a(3,:)=(/ 0, 2, 2/) 
      else if(cluster.eq.'16Z') then                           
         a(1,:)=(/2, 0, 0/); a(2,:)=(/ 0, 2, 0/); a(3,:)=(/ 0, 0, 4/) 
      else if(cluster.eq.'18A') then                           
         a(1,:)=(/ 2, 2, 2/); a(2,:)=(/ 2,-1,-1/); a(3,:)=(/ 1, 1,-2/)
      else if(cluster.eq.'18B') then                           
         a(1,:)=(/ 3, 2, 1/); a(2,:)=(/ 2,-1,-1/); a(3,:)=(/ 1, 1,-2/)
      else if(cluster.eq.'18C') then                           
         a(1,:)=(/ 2, 2, 2/); a(2,:)=(/ 2, 1,-1/); a(3,:)=(/-1, 2,-1/)
      else if(cluster.eq.'18D') then                           
         a(1,:)=(/ 0, 0, 3/); a(2,:)=(/ 2,-2,-1/); a(3,:)=(/ 2, 1,-1/)
      else if(cluster.eq.'18E') then                           
         a(1,:)=(/ 3, 3, 0/); a(2,:)=(/ 2, 0, 2/); a(3,:)=(/ 2, 0,-1/)
      else if(cluster.eq.'18F') then                           
         a(1,:)=(/ 3, 2, 1/); a(2,:)=(/ 0, 2,-2/); a(3,:)=(/ 0, 2, 1/)
         
      else if(cluster.eq.'20A') then                           
         a(1,:)=(/ 3, 2, 1/); a(2,:)=(/ 2,-1,-1/); a(3,:)=(/-1, 1,-2/)
      else if(cluster.eq.'20B') then                           
         a(1,:)=(/ 4, 1, 1/); a(2,:)=(/ 0, 2,-2/); a(3,:)=(/ 1, 2, 1/)
      else if(cluster.eq.'20C') then                           
         a(1,:)=(/ 2, 2, 2/); a(2,:)=(/ 2, 0,-2/); a(3,:)=(/-1, 1,-2/)
      else if(cluster.eq.'20D') then                           
         a(1,:)=(/ 3, 1, 0/); a(2,:)=(/ 0, 2, 2/); a(3,:)=(/ 1,-2, 1/)
      else if(cluster.eq.'20E') then                           
         a(1,:)=(/ 4, 1, 1/); a(2,:)=(/ 2,-2, 0/); a(3,:)=(/ 1,-1,-2/)
      else if(cluster.eq.'20F') then                           
         a(1,:)=(/ 2, 2, 1/); a(2,:)=(/ 2,-2, 1/); a(3,:)=(/ 1, 0,-2/)
      else if(cluster.eq.'20G') then                           
         a(1,:)=(/ 2, 2, 1/); a(2,:)=(/ 2,-2, 0/); a(3,:)=(/ 1, 1,-2/)
      else if(cluster.eq.'22A') then                           
         a(1,:)=(/ 4, 2, 0/); a(2,:)=(/ 0, 3,-1/); a(3,:)=(/-1, 2, 1/)
      else if(cluster.eq.'22B') then                           
         a(1,:)=(/ 2, 2, 2/); a(2,:)=(/ 3, 0,-1/); a(3,:)=(/-1, 1,-2/)
      else if(cluster.eq.'22C') then                           
         a(1,:)=(/ 3, 2, 0/); a(2,:)=(/ 0, 2, 2/); a(3,:)=(/ 1,-1, 2/)
      else if(cluster.eq.'22D') then                           
         a(1,:)=(/ 1, 2, 0/); a(2,:)=(/ 2, 1, 2/); a(3,:)=(/ 2,-2, 1/)
      else if(cluster.eq.'22E') then                           
         a(1,:)=(/ 4, 1, 1/); a(2,:)=(/ 2, 2,-3/); a(3,:)=(/ 0, 2,-1/)
      else if(cluster.eq.'24A') then                           
         a(1,:)=(/ 3, 1, 0/); a(2,:)=(/ 0, 2,-2/); a(3,:)=(/ 0, 2, 2/)
      else if(cluster.eq.'24B') then                           
         a(1,:)=(/ 2, 2, 2/); a(2,:)=(/ 0, 2,-2/); a(3,:)=(/ 2, 0,-2/)
      else if(cluster.eq.'24C') then                           
         a(1,:)=(/ 2, 2, 2/); a(2,:)=(/ 3, 0,-1/); a(3,:)=(/-1, 2,-1/)
      else if(cluster.eq.'24D') then                           
         a(1,:)=(/ 0, 0, 3/); a(2,:)=(/ 1, 3, 0/); a(3,:)=(/-2, 2, 0/)
      else if(cluster.eq.'24E') then                           
         a(1,:)=(/ 3, 2, 0/); a(2,:)=(/ 3,-2, 0/); a(3,:)=(/ 1, 0,-2/)
      else if(cluster.eq.'26A') then                           
         a(1,:)=(/ 3, 2, 1/); a(2,:)=(/ 2,-2,-2/); a(3,:)=(/ 1, 1,-2/)
      else if(cluster.eq.'26B') then                           
         a(1,:)=(/ 1, 3, 2/); a(2,:)=(/ 3,-2, 1/); a(3,:)=(/ 3, 1, 0/)
      else if(cluster.eq.'26C') then                           
         a(1,:)=(/ 2, 1, 3/); a(2,:)=(/ 2,-1,-2/); a(3,:)=(/ 2, 2,-1/)
      else if(cluster.eq.'26Z') then                           
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3, 3, -2/); a(3,:)=(/ -3,-2,-3/)
! Automatically generated code from gen_cluster.f:
      else if(cluster.eq."28A") then
! Cubicity=  1.06296811 Perfection=  0 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 1, 3/); a(2,:)=(/ 3,-1, 1/); a(3,:)=(/ 1, 2,-2/)
      else if(cluster.eq."28B") then
! Cubicity=  1.01773473 Perfection=  5 Bipartide=T ID#=199
         a(1,:)=(/ 1, 1, 2/); a(2,:)=(/ 3, 2,-1/); a(3,:)=(/ 1,-3, 2/)
      else if(cluster.eq."30A") then
! Cubicity=  1.00723632 Perfection=  0 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 2, 2/); a(2,:)=(/ 2, 2,-2/); a(3,:)=(/ 2,-2, 1/)
      else if(cluster.eq."30B") then
! Cubicity=  1.01243421 Perfection=  4 Bipartide=T ID#=124
         a(1,:)=(/ 1, 1, 2/); a(2,:)=(/ 3, 1,-2/); a(3,:)=(/ 3,-2, 1/)
      else if(cluster.eq."32A") then
! Cubicity=  1.02208161 Perfection=  0 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 1, 3/); a(2,:)=(/ 2, 2,-2/); a(3,:)=(/ 2,-2, 1/)
      else if(cluster.eq."32B") then
! Cubicity=  1.02814321 Perfection=  3 Bipartide=T ID#=100
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 2, 0,-2/); a(3,:)=(/ 2,-2, 2/)
      else if(cluster.eq."34A") then
! Cubicity=  1.00907609 Perfection=  0 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 1, 3/); a(2,:)=(/ 3,-2, 0/); a(3,:)=(/ 1, 2,-2/)
      else if(cluster.eq."34B") then
! Cubicity=  1.05741559 Perfection=  2 Bipartide=T ID#= 32
         a(1,:)=(/ 1, 0, 3/); a(2,:)=(/ 2, 2,-2/); a(3,:)=(/ 1,-3,-2/)
      else if(cluster.eq."36A") then
! Cubicity=  1.00419451 Perfection=  0 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 2, 2/); a(2,:)=(/ 3, 0,-2/); a(3,:)=(/ 2,-2, 2/)
      else if(cluster.eq."36B") then
! Cubicity=  1.04017240 Perfection=  3 Bipartide=T ID#=135
         a(1,:)=(/ 1, 0, 3/); a(2,:)=(/ 3, 2,-1/); a(3,:)=(/ 2,-2,-2/)
      else if(cluster.eq."38A") then
! Cubicity=  1.00210705 Perfection=  0 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 1, 3/); a(2,:)=(/ 3, 1,-3/); a(3,:)=(/ 2,-2, 1/)
      else if(cluster.eq."38B") then
! Cubicity=  1.08723039 Perfection=  0 Bipartide=T ID#= 12
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3,-1,-2/); a(3,:)=(/ 2,-2, 2/)
      else if(cluster.eq."40A") then
! Cubicity=  1.00299531 Perfection=  0 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 2, 2/); a(2,:)=(/ 3, 1,-2/); a(3,:)=(/ 1,-3, 2/)
      else if(cluster.eq."40B") then
! Cubicity=  1.04098349 Perfection=  3 Bipartide=T ID#=130
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 2, 2,-2/); a(3,:)=(/ 2,-2, 2/)
      else if(cluster.eq."42A") then
! Cubicity=  1.00530051 Perfection=  0 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 2, 2/); a(2,:)=(/ 3, 0,-2/); a(3,:)=(/ 0, 3,-3/)
      else if(cluster.eq."42B") then
! Cubicity=  1.05586176 Perfection=  2 Bipartide=T ID#= 54
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3,-1, 2/); a(3,:)=(/ 2, 2,-2/)
      else if(cluster.eq."44A") then
! Cubicity=  1.01039979 Perfection=  0 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 2, 2/); a(2,:)=(/ 3, 2,-2/); a(3,:)=(/ 3,-2, 1/)
      else if(cluster.eq."44B") then
! Cubicity=  1.03550915 Perfection=  3 Bipartide=T ID#= 44
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3, 2,-1/); a(3,:)=(/ 2,-2, 2/)
      else if(cluster.eq."46A") then
! Cubicity=  1.01422437 Perfection=  0 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 1, 3/); a(2,:)=(/ 3, 2,-2/); a(3,:)=(/ 3,-2, 0/)
      else if(cluster.eq."46B") then
! Cubicity=  1.01749200 Perfection=  4 Bipartide=T ID#= 34
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3, 1,-2/); a(3,:)=(/ 2,-2, 2/)
      else if(cluster.eq."48A") then
! Cubicity=  1.00860570 Perfection=  0 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 1, 3/); a(2,:)=(/ 3, 2,-2/); a(3,:)=(/ 2,-3,-1/)
      else if(cluster.eq."48B") then
! Cubicity=  1.00194693 Perfection=  5 Bipartide=T ID#=112
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3,-2, 1/); a(3,:)=(/ 2, 2,-2/)
      else if(cluster.eq."50A") then
! Cubicity=  1.00500070 Perfection=  1 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 1, 3/); a(2,:)=(/ 3, 2,-2/); a(3,:)=(/ 2,-3, 1/)
      else if(cluster.eq."50B") then
! Cubicity=  1.01815330 Perfection=  6 Bipartide=T ID#= 71
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3, 2,-1/); a(3,:)=(/ 2,-3, 1/)
      else if(cluster.eq."52A") then
! Cubicity=  1.10921976 Perfection=  1 Bipartide=F ID#=  1
         a(1,:)=(/ 2, 2, 3/); a(2,:)=(/ 3, 2,-2/); a(3,:)=(/ 3,-2,-2/)
      else if(cluster.eq."52B") then
! Cubicity=  1.00288782 Perfection=  7 Bipartide=T ID#= 57
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3, 1,-2/); a(3,:)=(/ 2,-3, 1/)
      else if(cluster.eq."54A") then
! Cubicity=  1.06306393 Perfection=  2 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3, 0,-3/); a(3,:)=(/ 3,-2, 2/)
      else if(cluster.eq."54B") then
! Cubicity=  1.00461045 Perfection=  8 Bipartide=T ID#= 88
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3,-3, 0/); a(3,:)=(/ 2, 2,-2/)
      else if(cluster.eq."56A") then
! Cubicity=  1.00256259 Perfection=  3 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 1, 3/); a(2,:)=(/ 3, 2,-2/); a(3,:)=(/ 3,-3,-1/)
      else if(cluster.eq."56B") then
! Cubicity=  1.02945995 Perfection=  9 Bipartide=T ID#= 34
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3, 2,-3/); a(3,:)=(/ 3,-1, 2/)
      else if(cluster.eq."58A") then
! Cubicity=  1.01420938 Perfection=  3 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 1, 3/); a(2,:)=(/ 3, 2,-2/); a(3,:)=(/ 3,-3, 1/)
      else if(cluster.eq."58B") then
! Cubicity=  1.01101713 Perfection= 10 Bipartide=T ID#= 18
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3,-3, 2/); a(3,:)=(/ 2, 2,-2/)
      else if(cluster.eq."60A") then
! Cubicity=  1.00088509 Perfection=  4 Bipartide=F ID#=  1
         a(1,:)=(/ 2, 0, 3/); a(2,:)=(/ 2, 3,-2/); a(3,:)=(/ 2,-3,-2/)
      else if(cluster.eq."60B") then
! Cubicity=  1.01114956 Perfection= 11 Bipartide=T ID#= 81
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3,-3, 2/); a(3,:)=(/ 2, 1,-3/)
      else if(cluster.eq."62A") then
! Cubicity=  1.08718883 Perfection=  5 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 3, 3/); a(2,:)=(/ 3, 2,-2/); a(3,:)=(/ 3,-3,-1/)
      else if(cluster.eq."62B") then
! Cubicity=  1.00308829 Perfection= 12 Bipartide=T ID#= 12
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3, 2,-1/); a(3,:)=(/ 3,-3, 2/)
      else if(cluster.eq."64A") then
! Cubicity=  1.01251339 Perfection=  6 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3, 1,-3/); a(3,:)=(/ 2,-3, 2/)
      else if(cluster.eq."64B") then
! Cubicity=  1.01039914 Perfection= 12 Bipartide=T ID#= 22
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3, 2,-3/); a(3,:)=(/ 2,-3, 1/)
      else if(cluster.eq."66A") then
! Cubicity=  1.06679699 Perfection=  5 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 3, 3/); a(2,:)=(/ 3, 3,-1/); a(3,:)=(/ 2,-3, 2/)
      else if(cluster.eq."66B") then
! Cubicity=  1.02569249 Perfection= 11 Bipartide=T ID#= 30
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3, 0,-3/); a(3,:)=(/ 3,-3, 2/)
      else if(cluster.eq."68A") then
! Cubicity=  1.05536589 Perfection=  4 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 3, 3/); a(2,:)=(/ 3, 3,-1/); a(3,:)=(/ 3,-2, 2/)
      else if(cluster.eq."68B") then
! Cubicity=  1.05400291 Perfection= 10 Bipartide=T ID#=  9
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3, 3,-2/); a(3,:)=(/ 2,-3, 3/)
      else if(cluster.eq."70A") then
! Cubicity=  1.06296488 Perfection=  3 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 3, 3/); a(2,:)=(/ 3, 3,-2/); a(3,:)=(/-2, 2,-3/)
      else if(cluster.eq."70B") then
! Cubicity=  1.03408521 Perfection=  9 Bipartide=T ID#=  5
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3, 3,-2/); a(3,:)=(/ 3,-2, 3/)
      else if(cluster.eq."72A") then
! Cubicity=  1.01027869 Perfection=  1 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 3, 3/); a(2,:)=(/ 3, 2,-2/); a(3,:)=(/ 3,-2, 2/)
      else if(cluster.eq."72B") then
! Cubicity=  1.13073007 Perfection= 10 Bipartide=T ID#= 58
         a(1,:)=(/ 2, 3, 3/); a(2,:)=(/ 2, 3,-3/); a(3,:)=(/ 2,-3, 3/)
      else if(cluster.eq."74A") then
! Cubicity=  1.06582181 Perfection=  3 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 3, 3/); a(2,:)=(/ 3,-1,-3/); a(3,:)=(/-3, 3,-2/)
      else if(cluster.eq."74B") then
! Cubicity=  1.02702810 Perfection=  9 Bipartide=T ID#=  3
         a(1,:)=(/ 1, 2, 3/); a(2,:)=(/ 3,-3, 2/); a(3,:)=(/ 2, 3,-3/)
      else if(cluster.eq."76A") then
! Cubicity=  1.01047869 Perfection=  0 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 3, 3/); a(2,:)=(/ 3, 2,-2/); a(3,:)=(/ 3,-3, 1/)
      else if(cluster.eq."78A") then
! Cubicity=  1.02391988 Perfection=  1 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 3, 3/); a(2,:)=(/ 3, 2,-2/); a(3,:)=(/ 2,-3, 3/)
      else if(cluster.eq."78B") then
! Cubicity=  1.10279190 Perfection=  5 Bipartide=T ID#= 12
         a(1,:)=(/ 2, 3, 3/); a(2,:)=(/ 3, 3,-2/); a(3,:)=(/ 3,-3,-2/)
      else if(cluster.eq."80A") then
! Cubicity=  1.01381121 Perfection=  1 Bipartide=F ID#=  1
         a(1,:)=(/ 1, 3, 3/); a(2,:)=(/ 3, 1,-3/); a(3,:)=(/ 3,-3, 1/)
      else if(cluster.eq."80B") then
! Cubicity=  1.09411902 Perfection=  4 Bipartide=T ID#=  3
         a(1,:)=(/ 2, 3, 3/); a(2,:)=(/ 3, 3,-2/); a(3,:)=(/-3, 2,-3/)
      else
         if(myrank.eq.0) write(lud,*) 'Bad cluster. Go directly to',
     &        ' jail. Do not pass go, do not collect $200'
         stop
      end if  
      if(iprint.ge.5.and.myrank.eq.0) 
     &   write(lud,"('a=',3('(',i2,1x,i2,1x,i2,') '))") 
     &   a(1,:),a(2,:),a(3,:)  

c     Test the a vectors      
      if(abs(determinant(a)).ne.Nc) then
        if(myrank.eq.0) write(lud,*) 'This cell volume isnt Nc',
     &                                abs(determinant(a))
        stop
      end if
c     Check fixed dimensions (PK: Could be allocated)
      if (nc.gt.maxnc) then
        if(myrank.eq.0) write(lud,*) 'Nc>maxnc parameter',nc,maxnc,
     &        'Need to recompile'
        stop
      end if
c     Calculate the cubicity of the clusters      
c
c       c8-------c7    body diagonals   face diagonals
c      /|       /|     --------------   --------------
c     / |      / |       d1 (c1,c7)       f1 (c1,c3)
c    /  |     /  |       d2 (c2,c8)       f2 (c2,c4)
c   c5-------c6  |       d3 (c3,c5)       f3 (c1,c6)
c   |   c4---|---c3      d4 (c4,c6)       f4 (c2,c5)
c a3|  /     |  /                         f5 (c1,c9)
c   | / a2   | /                          f6 (c4,c5)
c   |/       |/
c   c1-------c2
c        a1
c
      corner(:,1) = zeror
      corner(:,2) = a(1,:)
      corner(:,3) = a(1,:)+a(2,:)
      corner(:,4) = a(2,:)
      corner(:,5) = a(3,:)
      corner(:,6) = a(3,:)+a(1,:)
      corner(:,7) = a(3,:)+a(1,:)+a(2,:)
      corner(:,8) = a(3,:)+a(2,:)
      
      Bettsd(1)=sqrt(sum((corner(:,1)-corner(:,7))**2))
      Bettsd(2)=sqrt(sum((corner(:,2)-corner(:,8))**2))
      Bettsd(3)=sqrt(sum((corner(:,3)-corner(:,5))**2))
      Bettsd(4)=sqrt(sum((corner(:,4)-corner(:,6))**2))
      
      Bettsf(1)=sqrt(sum((corner(:,1)-corner(:,3))**2))
      Bettsf(2)=sqrt(sum((corner(:,2)-corner(:,4))**2))
      Bettsf(3)=sqrt(sum((corner(:,1)-corner(:,6))**2))
      Bettsf(4)=sqrt(sum((corner(:,2)-corner(:,5))**2))
      Bettsf(5)=sqrt(sum((corner(:,1)-corner(:,8))**2))
      Bettsf(6)=sqrt(sum((corner(:,4)-corner(:,5))**2))
      
      Bettsl(1)=sqrt(dot_product(corner(:,2),corner(:,2)))
      Bettsl(2)=sqrt(dot_product(corner(:,4),corner(:,4)))
      Bettsl(3)=sqrt(dot_product(corner(:,5),corner(:,5)))
      
      meand=product(Bettsd)**(0.25000000000_kind)
      meanf=product(Bettsf)**(0.16666666667_kind)
      meanl=product(Bettsl)**(0.33333333333_kind)
      r1=1.732050808_kind*meanl/meand
      r2=1.414213562_kind*meanl/meanf
      cubicity=max(r1,oner/r1)*max(r2,oner/r2)
      if(myrank.eq.0) 
     &  write(lud,*) 'cubicity=',cubicity
      
      

c*********************************************************************************
c     Determine which group operations are allowed    
c*********************************************************************************
      
      ngroup=0; Lc=2*(Nc**.3333+1)
      if(myrank.eq.0.and.iprint.ge.5) then
        write(lud,*) '--------------------------------------------',
     &               '------------------' 
        write(lud,*) '             trans(a)        igroup L  ixy ix',
     &               'z iyz isx isy isz'
        write(lud,*) '--------------------------------------------',
     &               '------------------'
      end if
      do igroup=1,48
         a1=transform3D(a,igroup,ixy,ixz,iyz,isx,isy,isz)

         if (use_brute_force) then

c        look at the lattice points generated by a1 and see if they
c        correspond to those produced by a.
c        First, generate the a1 lattice points
         do i=-Lc,Lc
         do j=-Lc,Lc
         do k=-Lc,Lc
           denied=.true.
           iv1=(/i,j,k/)
           do l=1,3         ! x,ymz components of this lattice point
             rv1(l)=sum(iv1(:)*a1(:,l))
           end do
c          Now generate the a lattice points.  If a and a1 are the 
c          same lattice one of these points MUST correspond to ij
           do n=-4*Lc,4*Lc
           do m=-4*Lc,4*Lc
           do p=-4*Lc,4*Lc
             iv2=(/n,m,p/)
             do l=1,3         ! x,y components of this lattice point
               rv2(l)=sum(iv2(:)*a(:,l))
             end do
             if(sum(abs(rv1-rv2)).lt.eps) denied=.false.
           end do
           end do
           end do
           if(denied) exit !no corrs. point was found
         end do
         if(denied) exit
         end do
         if(denied) exit
         end do
         end if
         if (use_lin_alg) then
c PK If lattice points generated by "a1" vectors can be found in 
c lattice points generated by "a" vectors then (a1)x(a^-1) can only 
c contain integers.
            ainv=Inverse_3x3(a)
            a1ainv=matmul(a1,ainv)
            if (all(abs(a1ainv-nint(a1ainv)).le.eps)) then 
               denied_lin_alg=.false. ! matrix only contains integers
            else
               denied_lin_alg=.true.
            end if
         end if

c If both algorithms are enabled, check they agree
         if (use_brute_force.and.use_lin_alg) then
            if (denied.neqv.denied_lin_alg) then
               write (*,*) '*** BUG : Denied algorithm disagree'
               write (*,*) 'a11',a1(:,1)
               write (*,*) 'a12',a1(:,2)
               write (*,*) 'a13',a1(:,3)
               
               write (*,*) 'a1',a(:,1)
               write (*,*) 'a2',a(:,2)
               write (*,*) 'a3',a(:,3)
               write (*,*)
               write (*,*) 'ainv',ainv(:,1)
               write (*,*) 'ainv',ainv(:,2)
               write (*,*) 'ainv',ainv(:,3)
               write (*,*)
               
               write (*,*) 'a1ainv',a1ainv(:,1)
               write (*,*) 'a1ainv',a1ainv(:,2)
               write (*,*) 'a1ainv',a1ainv(:,3)
               write (*,*) 'Denied_lin_alg=',denied_lin_alg
               write (*,*) 'Denied_brute_force=',denied
               write (*,*)
               stop
            end if
         end if
         if (use_lin_alg) denied=denied_lin_alg

         if(denied.eqv..false.) then ! a corresponding point was found
            ngroup=ngroup+1
            group(ngroup)=igroup
         end if    

         if(myrank.eq.0.and.iprint.ge.5) then
           allowed=.true.;if (denied) allowed=.false.
           write(lud,"(3('(',i2,1x,i2,1x,i2,')'),2x,i2,2x,l,6(2x,i2))")
     &     a1(1,:),a1(2,:),a1(3,:),igroup,allowed,ixy(igroup),
     &     ixz(igroup),iyz(igroup),isx(igroup),isy(igroup),isz(igroup)
         end if
      end do      

     
      if(myrank.eq.0.and.iprint.ge.5) 
     &  write(lud,"('cluster=',a3,' ngroup=',i2,' groups=',48(i2,1x))") 
     &  cluster,ngroup,group(1:ngroup)

      
c*********************************************************************************
c       calculate the table of cluster point locations Rc(i,ic), i=1,2,3 for x,y,z
c*********************************************************************************       
c 
c  ^ y
c  |
cc-|-c---c-*-c---c---c---c               A point P is within the cluster 
c| | |   |   |   |   |   |               if a vector v from the box origin, 
cc-|-c---c---c---c---c---c               O, to the point satisfies
ca2^ * * * * * * * * |   |               
cc-|-c---c---c---c-*-c---c                        
c| | |   P   |   | * |   |                            _   _
cc-|-c--/c---c---c-*-c---c                   -eps  <  v__ai  < 1-eps
c| | | /v|   |   | * |   |                        
cc-|-c/--c---c---c-*-c---c 
c| | /   |   |   | * |   |               for all i=1,3 where v__ai is the
cc-|/0---c---c---c-*-c---c               projection of the vector v along
c| O--------------->---------> x         ai
cc---c---c---c---ca1-c---c                              -1
c|   |   |   |   |   |   |               v = A n;  n = A  v ; where A=transpose(a)
cc---c---c---c---c---c---c
c                                        v = i*a1 + j*a2 + k*a3
c
c                                            [ a1 ]
c                                        a = [ a2 ]
c                                            [ a3 ]
c                                        
      do i=1,3
        a1(:,i)=a(i,:)
      end do
      c1=Inverse(a1)
      if(iprint.ge.5.and.myrank.eq.0) write(lud,*) '  ' !0
      if(iprint.ge.5.and.myrank.eq.0) 
     &     write(lud,*) 'ic,Rcx(ic),Rcy(ic),Rcz(ic)' !0
      ic=0
      do i=-Nc,Nc
      do j=-Nc,Nc
      do k=-Nc,Nc
         rv1=(/i,j,k/)          ! The coordinates of the point.
         rv2=matmul(c1,rv1)
         if(rv2(1).lt.oner-eps.and.rv2(1).gt.zeror-eps.and.
     &      rv2(2).lt.oner-eps.and.rv2(2).gt.zeror-eps.and.
     &      rv2(3).lt.oner-eps.and.rv2(3).gt.zeror-eps ) then
            ic=ic+1
            Rc(:,ic)=(/i,j,k/) ! Rc, like i.e. a(1,1) is integer
            if(iprint.ge.5.and.myrank.eq.0)                             !0
     &           write(lud,"(1x,i2,3x,3(i2,6x))") ic,Rc(1,ic),Rc(2,ic), !0
     &                                            Rc(3,ic)              !0
            if(sum(abs(Rc(:,ic))).eq.0) icR0=ic
         end if
      end do
      end do
      end do
      if(ic.ne.Nc) then
         if(myrank.eq.0) write(lud,*) 'Bug in Rc(ic) loop' !0
         stop
      end if
      if(iprint.ge.5.and.myrank.eq.0) write(lud,*) 'icR0=',icR0
 
c***********************************************************************        
c     create a table mapping points outside the cluster back into it.
c***********************************************************************        
c     Note we use c1 from above!!        

      if(iprint.ge.5.and.myrank.eq.0) then                  !0
         write(lud,*) '  '                                  !0
         write(lud,*) 'Equivalency table r'                 !0
         write(lud,*) 'ic   igroup   icrequ(ic,ieqiv)'      !0
      end if                                                !0
      do ic=1,Nc
        do igroup=1,ngroup
c       Generate the equivalent points using the ngroup operations
          rv3=transform3D(Rc(:,ic),group(igroup),ixy,ixz,iyz,
     &                    isx,isy,isz)

! Map rv3 back into cluster
          rv2=matmul(c1,rv3)    !project rv3 onto a1,a2,a3
          do while (any(rv2.lt.zeror-eps).or.any(rv2.gt.oner-eps))
             do l=1,3
                if(rv2(l).lt.zeror-eps) then
                   rv3(:)=rv3(:)+a(l,:)
                else if(rv2(l).gt.oner-eps) then
                   rv3(:)=rv3(:)-a(l,:)
                end if
             end do
             rv2=matmul(c1,rv3) !project rv3 onto a1,a2,a3
          end do
c         Figure out which point rv3 is
          do jc=1,Nc
            if(sum(abs(rv3(:)-Rc(:,jc))).lt.eps) then ! rv3=Rc(:,ic)
               icrequ(ic,igroup)=jc
               exit
            end if
          end do
          if(jc.gt.Nc) then
            write(lud,*) 'No equivalent found in icrequ',ic,igroup
            stop
          end if 
          if(iprint.ge.5.and.myrank.eq.0) 
     &      write(lud,"(2x,i2,4x,i2,8x,i2)") ic,igroup,icrequ(ic,igroup) !0
        end do
      end do                               

c***********************************************************************        
c     Now create the tables for the differences of R
c***********************************************************************        
c     Note we use c1 from above!!   
     
      if(iprint.ge.5.and.myrank.eq.0) then
        write(lud,*) '  ' !0
        write(lud,*) 'differences in r' !0
        write(lud,*) '  ic   jc  icrdiff(ic,jc)' !0
      end if
      do i=1,Nc
      do j=1,Nc
         iv1(:)=Rc(:,i)-Rc(:,j)
         itest=0
         do n=-1,1        ! Look in the surrounding 
         do m=-1,1        ! 27 clusters for R by adding and
         do p=-1,1        ! subtracting a1, a2 and a3
            iv2=iv1 + n*a(1,:)+m*a(2,:)+p*a(3,:)
            rv2=matmul(c1,iv2)
            if(rv2(1).lt.oner-eps.and.rv2(1).gt.zeror-eps.and.
     &         rv2(2).lt.oner-eps.and.rv2(2).gt.zeror-eps.and.
     &         rv2(3).lt.oner-eps.and.rv2(3).gt.zeror-eps ) then
c              This point is within the cluster
c              Now we must find which point it is.
               do ic=1,Nc
                  if(sum(abs(iv2(:)-Rc(:,ic))).eq.0) then
                     icrdiff(i,j)=ic
                     itest=itest+1   ! complain if multiple solutions are found
                  end if
               end do
               if(iprint.ge.5.and.myrank.eq.0) 
     &              write(lud,"(2x,i2,4x,i2,6x,i2)") i,j,icrdiff(i,j) !0
            end if
         end do
         end do
         end do
         if(itest.ne.1) then
            write(lud,*) 'no icrdiff(i,j) found',i,j
            stop
         end if
      end do
      end do  
        
c***********************************************************************        
c     Now create the tables for the neighbors of cluster points
c***********************************************************************        
c     Note we use c1 from above, don't overwrite it!!  
c
c     first locate  the neighbors

      neighbor(1,:) = (/ 1, 0, 0/) ! neighbor in the +x direction
      neighbor(2,:) = (/ 0, 1, 0/) ! neighbor in the +y direction
      neighbor(3,:) = (/ 0, 1, 0/) ! neighbor in the +y direction
     
      if(iprint.ge.5.and.myrank.eq.0) then
        write(lud,*) '  '                      !0
        write(lud,*) ' nearest neighbors of ic '       !0
        write(lud,*) ' ic  neigh  nneighbor(neigh,ic)' !0
      end if

      do i=1,Nc
      do j=1,ndim
         iv1(:)=Rc(:,i)+ neighbor(j,:)
         itest=0
         do n=-1,1        ! Look in the surrounding 
         do m=-1,1        ! 27 clusters for R by adding and
         do p=-1,1        ! subtracting a1, a2 and a3
            iv2=iv1 + n*a(1,:)+m*a(2,:)+p*a(3,:)
            rv2=matmul(c1,iv2)
            if(rv2(1).lt.oner-eps.and.rv2(1).gt.zeror-eps.and.
     &         rv2(2).lt.oner-eps.and.rv2(2).gt.zeror-eps.and.
     &         rv2(3).lt.oner-eps.and.rv2(3).gt.zeror-eps ) then
c              This point is within the cluster
c              Now we must find which point it is.
               do ic=1,Nc
                  if(sum(abs(iv2(:)-Rc(:,ic))).eq.0) then
                     nneighbor(j,i)=ic
                     itest=itest+1   ! complain if multiple solutions are found
                  end if
               end do
               if(iprint.ge.5.and.myrank.eq.0) 
     &              write(lud,"(2x,i2,4x,i2,6x,i2)") i,j,nneighbor(j,i) !0
            end if
         end do
         end do
         end do
         if(itest.ne.1) then
            write(lud,*) 'no nneighbor(j,i) found',j,i
            stop
         end if
      end do
      end do  


c
c     Now count the number of neighbors in each shell  
c

     
      if(iprint.ge.5.and.myrank.eq.0) then
        write(lud,*) '  ' !0
        write(lud,*) 'neighbors of ic' !0
        write(lud,*) 
     &    'shell  nneigh(shell)  nneighl(shell) nneighp(shell)' !0
      end if
      i1 = 2*(Nc+1)**(oner/ndim); rv1=zeror; nneigh=0; ishellmax=0
      nneighl(:)=0 
      nneighp(:)=0

      if(iprint.ge.5.and.myrank.eq.0) then
         write (*,*) 'i1=',i1
      end if

      do i=-i1,i1
      do j=-i1,i1
      do k=-i1,i1
        iv1=(/ i, j, k/)
        i2=sum(abs(iv1)) ! the lattice neighbor index (according to Betts)
        if(i2.gt.0) then
          if(i2.le.nshellmax) then
            nneighl(i2)=nneighl(i2)+1
          end if
          do n=1,ndim
            rv1(n)=iv1(n) + halfr*sum(a(:,n))
          end do
          rv2=matmul(c1,rv1)
          if(rv2(1).lt.oner-eps.and.rv2(1).gt.zeror-eps.and.
     &       rv2(2).lt.oner-eps.and.rv2(2).gt.zeror-eps.and.
     &       rv2(3).lt.oner-eps.and.rv2(3).gt.zeror-eps ) then
c            This point is within the cluster. Now we need to find
c            the shortest distance between this point and icR0, and
c            the corresponding cluster neighbor index i3.
             i3=i2
             do n=-1,1
             do m=-1,1
             do p=-1,1
               iv2=iv1 + n*a(1,:) + m*a(2,:) + p*a(3,:)
               if(sum(abs(iv2)).lt.i2) i3=sum(abs(iv2))
             end do !p
             end do !m
             end do !n
             nneigh(i3)=nneigh(i3)+1
             if(i3.gt.ishellmax) ishellmax=i3
          end if
        end if
      end do !k
      end do !j
      end do !i
      
      i2=0
      do while (sum(nneighp).lt.Nc-1)
        i2=i2+1
        nneighp(i2)=nneighl(i2)
      end do
      nneighp(i2)=Nc-1-sum(nneighp(1:i2-1))
      
      if(sum(nneigh)+1.ne.Nc.or.sum(nneighp)+1.ne.Nc) then
        if(myrank.eq.0) write(lud,*) 'error in nneigh'
        if(myrank.eq.0) write(lud,*) 'sum(nneigh)=',sum(nneigh)
        if(myrank.eq.0) write(lud,*) 'sum(nneighp)=',sum(nneighp)
        stop
      end if
        
      if(iprint.ge.5.and.myrank.eq.0) then
        do i1=1,ishellmax
          write(lud,"(2x,i3,8x,i3,10x,i3,10x,i3)") 
     &      i1,nneigh(i1),nneighl(i1),nneighp(i1)
        end do
        write(lud,*) 'ishellmax=',ishellmax
      end if
c
c     Now determine the perfection.  It seems that Betts defines the
c     ferromagnetic imperfection I_F as the numer of additional spins in
c     the shells beyond j plus the number missing the shells before j,
c     and not counting those in the jth shell.  The integer j seems to
c     be the number which yeilds the smallest I_F

      do j=1,ishellmax
        i3=0
        do i1=1,j-1
          i3=i3+abs(nneigh(i1)-nneighl(i1))
        end do
        do i1=j+1,ishellmax
          i3=i3+nneigh(i1)
        end do
        if(j.gt.1) then
          if(i3.lt.I_F) I_F=i3 
        else
          i_F=i3
        end if
      end do
      
      if(myrank.eq.0) 
     &  write(lud,*) 'ferromagnetic imperfection I_F=',I_F
     
c     Betts proposes a simpler algorithm to calculate the perfection,
c     involving the difference of two quantities.
c
c     F_f = SUM k nneigh(k)
c            k
c
c     F_p = SUM k nneighp(k)
c            k

      i2=0; i3=0
      do i1=1,ishellmax
        i2=i2+i1*nneigh(i1)
        i3=i3+i1*nneighp(i1)
      end do
      if(myrank.eq.0) 
     &  write(lud,*) 'Betts imperfection I_F=',i2-i3
     
      

c***********************************************************************        
c       Now set up the K-space tables
c***********************************************************************        

c       Now calculate the principle translations in K-space
c   
c  
c                              _    _
c    _     _    _        2 pi aj x ak          
c    a --> g    g_i =  ------------------      plus cyclic permutations
c                      | ai . (aj x ak)|       
c                                           
c       Principle translation vectors for lattice tiling       
c                                        
c
      r1=determinant(a)
      
      g(1,1)=2*pi*( a(2,2)*a(3,3)-a(3,2)*a(2,3) )/r1  !a2y*a3z-a3y*a2z
      g(1,2)=2*pi*( a(2,3)*a(3,1)-a(3,3)*a(2,1) )/r1  !a2z*a3x-a3z*a2x
      g(1,3)=2*pi*( a(2,1)*a(3,2)-a(3,1)*a(2,2) )/r1  !a2x*a3y-a3x*a2y
      
      g(2,1)=2*pi*( a(3,2)*a(1,3)-a(1,2)*a(3,3) )/r1  !a3y*a1z-a1y*a3z
      g(2,2)=2*pi*( a(3,3)*a(1,1)-a(1,3)*a(3,1) )/r1  !a3z*a1x-a1z*a3x
      g(2,3)=2*pi*( a(3,1)*a(1,2)-a(1,1)*a(3,2) )/r1  !a3x*a1y-a1x*a3y
        
      g(3,1)=2*pi*( a(1,2)*a(2,3)-a(2,2)*a(1,3) )/r1  !a1y*a2z-a2y*a1z
      g(3,2)=2*pi*( a(1,3)*a(2,1)-a(2,3)*a(1,1) )/r1  !a1z*a2x-a2z*a1x
      g(3,3)=2*pi*( a(1,1)*a(2,2)-a(2,1)*a(1,2) )/r1  !a1x*a2y-a2x*a1y
      
      r1=abs(determinant(g))
      if(abs((2*pi)**3/r1-Nc).gt.eps) then
        if(myrank.eq.0) write(lud,*) 'This cell volume aint Nc',r1
        stop
      end if
      do i=1,3
      do j=1,3
        r1=sum(a(i,:)*g(j,:))/(2*pi)
        if(i.eq.j.and.abs(r1-oner).gt.eps) then
          if(myrank.eq.0) write(lud,*) 'bug in g',i,j,r1
          stop
        end if
        if(i.ne.j.and.abs(r1).gt.eps) then
          if(myrank.eq.0) write(lud,*) 'bug in g',i,j,r1
          stop
        end if
      end do
      end do
      
      if(iprint.ge.5.and.myrank.eq.0) 
     &  write(lud,"(' g=',3(' (',f5.2,1x,f5.2,1x,f5.2,')'))") 
     &  (g(l,:),l=1,3)
     
     
c***********************************************************************        
c     Create the cluster momentum tables Kc(i,ic), i=x,y,z
c***********************************************************************        

c     First find the points in the irreducible wedge.  They will be
c     indexed from 1 to Ncw
      if(iprint.ge.5.and.myrank.eq.0) then
        write(lud,*) '  Wedge K points  '    !0
        write(lud,*) 'ic  i  j  k  Kcx(ic)  Kcy(ic)  Kcz(ic)  ' !0
      end if
      ic=0; Ncw=0; icK0=0; icKpi=0
      do i=-Nc,Nc    !times a1
      do j=-Nc,Nc    !times a2
      do k=-Nc,Nc    !times a3
        iv1=(/i,j,k/)
        do l=1,3         !x,y,z components 
          rv1(l)=sum(iv1(:)*g(:,l))
        end do
        if(rv1(1).gt.-pi+eps.and.rv1(1).lt.pi+eps.and.
     &     rv1(2).gt.-pi+eps.and.rv1(2).lt.pi+eps.and.
     &     rv1(3).gt.-pi+eps.and.rv1(3).lt.pi+eps) then
c          The point is in the first SC Brillouin zone.  Now see if
c          the point falls within the IRW.  First see if it is 
c          equivalent to any known point in the IRW (if so it cannot be in IRW)
          Kcops: do igroup=2,ngroup ! exclude the identity
            rv2=transform3D(rv1,group(igroup),ixy,ixz,iyz,isx,isy,isz)
            itest=0
            do jc=1,Ncw
              if(sum(abs(Kc(:,jc)-rv2(:))).lt.eps) then ! rv2 = Kc(jc)
                 exit Kcops
              end if
            end do
          end do  Kcops
          if(igroup.gt.ngroup) then ! rv1 is in 1BZ and not equivalent to any point  
             Ncw=Ncw+1        ! in the known IRW.  It must be in the IRW
             ic=ic+1
             Kc(:,ic)=rv1(:)
             ig(ic)=i
             jg(ic)=j
             kg(ic)=k
             if(sum(abs(rv1-pi)).lt.eps) icKpi=ic
             if(sum(abs(rv1-zeror)).lt.eps) icK0=ic
             if(iprint.ge.5.and.myrank.eq.0) 
     &            write(lud,5) ic,i,j,k,Kc(1,ic),Kc(2,ic),Kc(3,ic) !0
          end if
        end if 
      end do 
      end do 
      end do 
        
      if(icK0.eq.0) then  ! we failed to find (0,0)
        if(myrank.eq.0) 
     &    write(lud,*) 'failed to find (0,0) in the IRW'    !0
        stop
      else if(icKpi.eq.0) then  ! we failed to find (pi,pi)
        icKpi=icK0
        if(myrank.eq.0) then
          write(lud,*) '*****************************************'    !0
          write(lud,*) '* failed to find (pi,pi,pi) in the IRW  *'    !0
          write(lud,*) '* apparently the lattice isnt bipartide *'    !0
          write(lud,*) '* icKpi=icK0 so that the code will run  *'    !0
          write(lud,*) '*****************************************'    !0
        end if
      end if

c     Now find the points in the 1BZ but outside the irreducible wedge.  
c     They will be indexed from Ncw+1 to Nc
      if(iprint.ge.5.and.myrank.eq.0) then
        write(lud,*) ' Non Wedge K points '    !0
      end if
      do i=-Nc,Nc              
      do j=-Nc,Nc
      do k=-Nc,Nc
         iv1=(/i,j,k/)
         do l=1,3         !x,y,z components 
            rv1(l)=sum(iv1(:)*g(:,l)) ! K point in xyz
         end do
         if(rv1(1).gt.-pi+eps.and.rv1(1).lt.pi+eps.and.
     &      rv1(2).gt.-pi+eps.and.rv1(2).lt.pi+eps.and.
     &      rv1(3).gt.-pi+eps.and.rv1(3).lt.pi+eps) then
c           The point is in the first SC Brillouin zone.  Now see if
c           the point falls within the IRW.  
            do jc=1,Ncw
               if(sqrt(sum((Kc(:,jc)-rv1(:))**2)).lt.eps) then
c                 rv1 and Kc(ic) are the same point
                  exit
               end if
            end do
            if(jc.gt.Ncw) then !rv1 is in 1BZ and not in the IRW.
               ic=ic+1
               Kc(:,ic)=rv1(:)
               ig(ic)=i
               jg(ic)=j
               kg(ic)=k
               if(iprint.ge.5.and.myrank.eq.0) 
     &              write(lud,5) ic,i,j,k,Kc(1,ic),Kc(2,ic),Kc(3,ic) !0
            end if
         end if 
      end do 
      end do 
      end do 
        
 5    format(1x,i2,1x,i2,1x,i2,1x,i2,2x,f6.3,3x,f6.3,3x,f6.3)
      
      ntw=nl*Ncw      
      
      if(ic.ne.Nc) then
         if(myrank.eq.0) write(lud,*) 'Bug in Kc(ic) loop' !0
         stop
      end if
      if(iprint.ge.5.and.myrank.eq.0) 
     &  write(lud,"('icK0=',i2,' icKpi=',i2,' Ncw=',i2,' Nc=',i2)") 
     &  icK0,icKpi,Ncw,ic


c***********************************************************************
c     Create the indirect addressing table ict(i,j,k) which maps a 
c     point indexed by the principle translations in k-space to ic.
c***********************************************************************        
      N_sp=sqrt(float(Nc))
      do i=-N_sp,N_sp
      do j=-N_sp,N_sp
      do k=-N_sp,N_sp
         iv1=(/i,j,k/)
         do l=1,3         !x,y,z components 
            rv1(l)=sum(iv1(:)*g(:,l)) ! K point in xyz
         end do
! Map rv1 into -pi->+pi BZ
         rv1=rv1-twor*pi*int(rv1/(twor*pi))
         where (rv1.gt.pi+eps) rv1=rv1-twor*pi
         where (rv1.lt.-pi+eps) rv1=rv1+twor*pi

         do jc=1,Nc
            if(sum(abs(rv1(:)-Kc(:,jc))).lt.eps) then    
               ict(i,j,k)=jc
               exit
            end if
         end do
         if(jc.gt.Nc) then
            write(lud,*) 'ict failed for i,j,k=',i,j,k
            write(lud,*) rv1
            stop
         end if
      end do
      end do
      end do

c***********************************************************************        
c     Create the table ickdiff(ic1,ic2), which takes the difference
c     between to K-points
c***********************************************************************        
      if(iprint.ge.5.and.myrank.eq.0) then
        write(lud,*) '   '                           !0
        write(lud,*) 'ic1   ic2  ickdiff(ic1,ic2)  ' !0
      end if
      do ic1=1,Nc
      do ic2=1,Nc
         rv1=Kc(:,ic1)-Kc(:,ic2)
! Map rv1 into -pi->+pi BZ
         rv1=rv1-twor*pi*int(rv1/(twor*pi))
         where (rv1.gt.pi+eps) rv1=rv1-twor*pi
         where (rv1.lt.-pi+eps) rv1=rv1+twor*pi
c        Now we need to figure out where this point is!
         do ic=1,Nc
            if(sum(abs(rv1(:)-Kc(:,ic))).lt.eps) then    
               ickdiff(ic1,ic2)=ic
               exit
            end if
         end do
         if(ic.gt.Nc) then
            write(lud,5) 'ickdiff failed for ic1,ic2=',ic1,ic2
            stop
         end if
         if(iprint.ge.5.and.myrank.eq.0) 
     &        write(lud,"(1x,i2,4x,i2,6x,i2)") ic1,ic2,ickdiff(ic1,ic2) !0
      end do
      end do
          
c***********************************************************************        
c     Create the table ickplus(ic1,ic2), which takes the sum
c     between to K-points
c***********************************************************************        
      if(iprint.ge.5.and.myrank.eq.0) then
        write(lud,*) '   '                           !0
        write(lud,*) 'ic1   ic2  ickplus(ic1,ic2)  ' !0
      end if
      do ic1=1,Nc
      do ic2=1,Nc
         rv1=Kc(:,ic1)+Kc(:,ic2)
! Map rv1 into -pi->+pi BZ
         rv1=rv1-twor*pi*int(rv1/(twor*pi))
         where (rv1.gt.pi+eps) rv1=rv1-twor*pi
         where (rv1.lt.-pi+eps) rv1=rv1+twor*pi
c        Now we need to figure out where this point is!
         do ic=1,Nc
            if(sum(abs(rv1(:)-Kc(:,ic))).lt.eps) then    
               ickplus(ic1,ic2)=ic
               exit
            end if
         end do
         if(ic.gt.Nc) then
            write(lud,5) 'ickplus failed for ic1,ic2=',ic1,ic2
            stop
         end if
         if(iprint.ge.5.and.myrank.eq.0) 
     &        write(lud,"(1x,i2,4x,i2,6x,i2)") ic1,ic2,ickplus(ic1,ic2) !0
      end do
      end do

c***********************************************************************
c     Now find the points equivalent to ic, ickequ(ic,igroup) where 
c     igroup is the group operation.
c***********************************************************************        
c
c
      if(iprint.ge.5.and.myrank.eq.0) then
        write(lud,*) '   '                            !0
        write(lud,*) 'equivalency in k'               !0
        write(lud,*) ' ic  igroup  ickequ(ic,igroup)' !0
      end if
      do ic=1,Nc
         i=ig(ic)
         j=jg(ic)
         k=kg(ic)
         iv1=(/i,j,k/)
         do l=1,3
           rv2(l)=sum(iv1(:)*g(:,l))
         end do
c        Kc(ic) falls within the zone, by construction, now look 
c        for equivalent points
c     
         do igroup=1,ngroup
            rv1=transform3D(rv2,group(igroup),ixy,ixz,iyz,isx,isy,isz) 
c           These are the equivalent points.  Now map to 1st BZ
         rv1=rv1-twor*pi*int(rv1/(twor*pi))
         where (rv1.gt.pi+eps) rv1=rv1-twor*pi
         where (rv1.lt.-pi+eps) rv1=rv1+twor*pi
c           Now we need to figure out where this point is.
            do jc=1,Nc
               if(sqrt(sum((rv1(:)-Kc(:,jc))**2)).lt.eps) then
                  ickequ(ic,igroup)=jc
                  exit
               end if
            end do
            if(iprint.ge.5.and.myrank.eq.0) 
     &        write(lud,"(2x,i2,4x,i2,6x,i2)") 
     &        ic,igroup,ickequ(ic,igroup)                  !0
            if(jc.gt.Nc) then
               if(myrank.eq.0) then
                 write(lud,*)'no equivalent K-point found' !0
                 write(lud,"('ic,igroup,group(igroup)',3(2x,i2))")
     &             ic,igroup,group(igroup)         !0
                 write(lud,"('K =',3(f6.3,1x))") rv2        !0
                 write(lud,"('K1=',3(f6.3,1x))") rv1       !0
               end if
               stop
            end if
         end do
      end do
        



c***********************************************************************        
c       Now generate a table which maps any K to an equivalent point
c       in the irreducible wedge (IW), and a table for the degeneracy 
c       of each point in the IW. 
c 
c       The use of the wedge can greatly reduce the storage required.
c       for example for Nc=8, only 4 of the 8 points are in the irreducible 
c       wedge.  Thus, all storage arrays may be reduced in size by a 
c       factor of 2.  As Nc-->oo, the saving factor approaches 48 for
c       clusters of high symmetry.  
c
c***********************************************************************        
c       
      
      ickdeg=0; irwgroup=0
      if(iprint.ge.5.and.myrank.eq.0) then
        write(lud,*) '   '                            !0
        write(lud,*) ' degeneracy in k'              !0
        write(lud,*) ' ic  ickmap(ic)  ickdeg(ic) irwgroup(ic)' !0
      end if
c     scan over all of the K points in the 1BZ
      do ic=1,Nc
         i=ig(ic)
         j=jg(ic)
         k=kg(ic)
         iv1=(/i,j,k/)
         do l=1,3
            rv2(l)=sum(iv1(:)*g(:,l))
         end do
c        Scan over all allowed point group operations, to find the
c        one that maps ic into the IRW
         groupops: do igroup=1,ngroup
           rv1=transform3D(rv2,group(igroup),ixy,ixz,iyz,isx,isy,isz) 
c          These are the equivalent points.  Now map to 1st BZ
         rv1=rv1-twor*pi*int(rv1/(twor*pi))
         where (rv1.gt.pi+eps) rv1=rv1-twor*pi
         where (rv1.lt.-pi+eps) rv1=rv1+twor*pi
c          See if this point is in the IRW!
           do jc=1,Ncw
              if(sum(abs(rv1(:)-Kc(:,jc))).lt.eps) then
                 ickmap(ic)=jc
                 ickdeg(jc)=ickdeg(jc)+1
                 irwgroup(ic)=igroup
                 exit groupops
              end if
           end do
         end do groupops
         if(igroup.gt.ngroup) then
            if(myrank.eq.0) write(lud,*) 'mapping failed for' !0
            if(myrank.eq.0) write(lud,*) 'ic',ic !0
            if(myrank.eq.0) write(lud,"('K=',3(f6.3,1x))") rv2 !0
            stop
         end if
         
      end do !ic

      if(iprint.ge.5.and.myrank.eq.0) then 
        do ic=1,Nc    
          write(lud,"(2x,i2,6x,i2,9x,i2,9x,i2)") 
     &    ic,ickmap(ic),ickdeg(ic),irwgroup(ic)
        end do
      end if

c     Test the degeneracy table.  Since ickdeg holds the degeneracy
c     of points in the IRW, its sum must equal Nc.
      if(sum(ickdeg(1:Ncw)).ne.Nc) then
         if(myrank.eq.0) then
           write(lud,*) 'bug in ickdeg',sum(ickdeg(1:Ncw))
           write(lud,*) 'ickdeg=',ickdeg(1:Ncw)
         end if
         stop
      end if
        
c     Test the group table.  irwgroup(ic within wedge)=1, the identity
      do ic=1,Ncw
         if(irwgroup(ic).ne.1) then
            if(myrank.eq.0) then
              write(lud,*) 'irwgroup failed for wedge pt.'     !0
              write(lud,*) 'ic',ic,'irwgroup(ic)',irwgroup(ic) !0
              write(lud,"('K=(',3(f6.3,1x),')')") Kc(:,ic)      !0
            end if
            stop
         end if
      end do


c***********************************************************************
c     generate the kt-points in the Brillouin-zone of the 
c     super-lattice, i.e. in the Wigner-Seitz cell of the  
c     reciprocal space defined by the cluster K-points
c***********************************************************************

      Dkt=oner/real(ntot,kind)
      nwsc=0 ! number of k-points in Wigner-Seitz cell 
      if (ntot.gt.0) then
      do i=-ntot+1,ntot
        do j=-ntot+1,ntot
          do k=1,ntot
           iv1=(/i,j,k/)
           rv1(1)=Dkt*(iv1(1)-halfr)*pi
           rv1(2)=Dkt*(iv1(2)-halfr)*pi
           rv1(3)=Dkt*(iv1(3)-halfr)*pi
c  Check if this k-point is in Wigner-Seitz cell, i.e. if the closest
c  K-point is K=0 
            r1=(rv1(1)*rv1(1)+rv1(2)*rv1(2)+rv1(3)*rv1(3)) ! distance from K=0
            itest=0
            do jc=1,Nc
            if (jc.eq.icK0) cycle
             r2=rv1(1)-Kc(1,jc)
             r3=rv1(2)-Kc(2,jc)
             r4=rv1(3)-Kc(3,jc)
             if((r2-pi).gt.zeror) then; r2=r2-twor*pi; endif
             if((r2+pi).lt.zeror) then; r2=r2+twor*pi; endif
             if((r3-pi).gt.zeror) then; r3=r3-twor*pi; endif
             if((r3+pi).lt.zeror) then; r3=r3+twor*pi; endif
             if((r4-pi).gt.zeror) then; r4=r4-twor*pi; endif
             if((r4+pi).lt.zeror) then; r4=r4+twor*pi; endif
             r2=(r2*r2+r3*r3+r4*r4)
             if(r2.lt.r1) then
               itest=1
               exit
             endif
c
            enddo ! jc
            if(itest.eq.0) then 
             nwsc=nwsc+1
             kt(1,nwsc) = rv1(1)
             kt(2,nwsc) = rv1(2)
             kt(3,nwsc) = rv1(3)
             nwsc=nwsc+1 ! -k must also be part of the WS cell
             kt(1,nwsc) = -rv1(1)
             kt(2,nwsc) = -rv1(2)
             kt(3,nwsc) = -rv1(3)
            endif
          enddo ! k
        enddo ! j
       enddo ! i
       else ! needed for FSS calculation
          nwsc=1
          kt(:,1)=zeror
       endif
c
       if(myrank.eq.0.and.iprint.ge.5) then
          write(lud,*) 'Number of points in WS, nwsc', nwsc
          write(lud,*) 'Nc times nwsc', Nc*nwsc
          write(lud,*) 'Total number of points in BZ', (2*ntot)**3 
       endif 

c***********************************************************************
c     Form Epsbar(K)
c***********************************************************************        
      if(myrank.eq.0.and.iprint.ge.5) write(lud,*) '    '
      if(myrank.eq.0.and.iprint.ge.5) write(lud,*) 
     &               'ic   Epsbar(ic)'

      Epsbar=zeror
      do ic=1,Nc
         do i=1,nwsc
            rv1(1) = Kc(1,ic)+kt(1,i)
            rv1(2) = Kc(2,ic)+kt(2,i)
            rv1(3) = Kc(3,ic)+kt(3,i)
            Epsbar(ic)= Epsbar(ic)
     &            - halfr*(cos(rv1(1))+cos(rv1(2))+cos(rv1(3))) 
     &            - tprime*(cos(rv1(1))*cos(rv1(2))
     &                     +cos(rv1(2))*cos(rv1(3))
     &                     +cos(rv1(1))*cos(rv1(3))-oner)
         end do
         Epsbar(ic)=Epsbar(ic)/real(nwsc,kind)
      if(myrank.eq.0.and.iprint.ge.5) write(lud,"(i2,6x,f9.6)") 
     &           ic,Epsbar(ic)
      end do
C     Check      
      r1=sum(Epsbar); r2=sum(Epsbar*ickdeg)
      if(abs(tprime).lt.eps.and.abs(r1).gt.eps
     &   .or.
     &  abs(tprime).lt.eps.and.abs(r2).gt.eps) then
        if(myrank.eq.0) then
          write(lud,*) 'sum Epsbar wrong.',r1,r2
        end if
c        stop
      end if
C      write(lud,*) 'sum Epsbar:    ',r1,r2

        
c***********************************************************************        
c       Form the Fourier transform tables.  
c       R-->K,  FTCoefs_K_to_R
c       K-->R,  FTCoefs_K_to_R
c
c                    __          the factor of 1/N is required
c                 1  \           so that both G(R)~1/w and G(K)~1/w
c       G(R=0) = --- /  G(K)     for large w
c                N   --
c                    K
c
c***********************************************************************        
      do ic=1,Nc
         do ik=1,Nc
            r1=zeror
            do idim=1,ndim
               r1=r1+Kc(idim,ik)*Rc(idim,ic)
            end do
            FTCoefs_K_to_R(ik,ic)=exp(-ii*r1)/real(Nc,kind)
            FTCoefs_R_to_K(ik,ic)=exp(+ii*r1)
         end do
      end do

c     Now test for consistency of Kc and Rc
      do ic=1,Nc
        r1= real(sum(FTCoefs_K_to_R(1:Nc,ic)))-kroniker(ic,icR0)
        r2=aimag(sum(FTCoefs_K_to_R(1:Nc,ic)))
        r3= real(sum(FTCoefs_K_to_R(ic,1:Nc)))-kroniker(ic,icK0)
        r4=aimag(sum(FTCoefs_K_to_R(ic,1:Nc)))
        
        if(abs(r1).gt.eps.or.abs(r2).gt.eps) then
          write(lud,*) 'Kc and Rc are inconsistent'
          write(lud,*) 'icr,r1,r2=',ic,r1,r2
          stop
        end if
        if(abs(r3).gt.eps.or.abs(r4).gt.eps) then
          write(lud,*) 'Kc and Rc are inconsistent, K=0'
          write(lud,*) 'ick,r3,r4=',ic,r3,r4
          stop
        end if
      end do


      return
      end  
      
      
      
c***********************************************************************        
c***********************************************************************        
c         3D subroutines and functions
c***********************************************************************        
c***********************************************************************
      
      
      
        
      function gtransform3D(a,igroup,ixy,ixz,iyz,isx,isy,isz)
c
c  x <-- isx*[(1-ixy)*x + ixy*y]  <-- [(1-ixz)*x + ixz*z]<----------------------
c  y <-- isy*[(1-ixy)*y + ixy*x]  <----------------------<-- [(1-iyz)*y + iyz*z]
c  z <-- isz*<------------------- <-- [(1-ixz)*z + ixz*x]<-- [(1-iyz)*z + iyz*y]
c
      implicit none
      integer :: i,a(3,3),g(3,3),gtransform3D(3,3),igroup,
     &           ixy(48),ixz(48),iyz(48),isx(48),isy(48),isz(48)
      
      do i=1,3  
         gtransform3D(i,1)=a(i,1)
         gtransform3D(i,2)=(1-iyz(igroup))*a(i,2) + iyz(igroup)*a(i,3)
         gtransform3D(i,3)=(1-iyz(igroup))*a(i,3) + iyz(igroup)*a(i,2)
      end do
               
      do i=1,3  
         g(i,1)= (1-ixz(igroup))*gtransform3D(i,1)
     &         + ixz(igroup)*gtransform3D(i,3)
         g(i,2)=                gtransform3D(i,2)
         g(i,3)= (1-ixz(igroup))*gtransform3D(i,3) 
     &         + ixz(igroup)*gtransform3D(i,1)
      end do

      do i=1,3  
         gtransform3D(i,1)=isx(igroup)*((1-ixy(igroup))*g(i,1) 
     &                  + ixy(igroup)*g(i,2))
         gtransform3D(i,2)=isy(igroup)*((1-ixy(igroup))*g(i,2) 
     &                  + ixy(igroup)*g(i,1))
         gtransform3D(i,3)=isz(igroup)*                 g(i,3)         
      end do
      
      return
      end
        
      function Rtransform3D(Rc,igroup,ixy,ixz,iyz,isx,isy,isz)
c
c  x <-- isx*[(1-ixy)*x + ixy*y]  <-- [(1-ixz)*x + ixz*z]<----------------------
c  y <-- isy*[(1-ixy)*y + ixy*x]  <----------------------<-- [(1-iyz)*y + iyz*z]
c  z <-- isz*<------------------- <-- [(1-ixz)*z + ixz*x]<-- [(1-iyz)*z + iyz*y]
c
      implicit none
      integer :: i,Rc(3),g(3),Rtransform3D(3),igroup, 
     &           ixy(48),ixz(48),iyz(48),isx(48),isy(48),isz(48)
      
      Rtransform3D(1)=            Rc(1)
      Rtransform3D(2)=(1-iyz(igroup))*Rc(2) + iyz(igroup)*Rc(3)
      Rtransform3D(3)=(1-iyz(igroup))*Rc(3) + iyz(igroup)*Rc(2)
      
      g(1)=(1-ixz(igroup))*Rtransform3D(1) + ixz(igroup)*Rtransform3D(3)
      g(2)=                Rtransform3D(2)
      g(3)=(1-ixz(igroup))*Rtransform3D(3) + ixz(igroup)*Rtransform3D(1)
      
      Rtransform3D(1)=isx(igroup)*((1-ixy(igroup))*g(1) 
     &             + ixy(igroup)*g(2))
      Rtransform3D(2)=isy(igroup)*((1-ixy(igroup))*g(2) 
     &             + ixy(igroup)*g(1))
      Rtransform3D(3)=isz(igroup)*                 g(3) 
      
      return
      end
      
      function Ktransform3D(Kc,igroup,ixy,ixz,iyz,isx,isy,isz)
c     
c  x <-- isx*[(1-ixy)*x + ixy*y]  <-- [(1-ixz)*x + ixz*z]<----------------------
c  y <-- isy*[(1-ixy)*y + ixy*x]  <----------------------<-- [(1-iyz)*y + iyz*z]
c  z <-- isz*<------------------- <-- [(1-ixz)*z + ixz*x]<-- [(1-iyz)*z + iyz*y]
c
      implicit none
      integer :: i,igroup,
     &           ixy(48),ixz(48),iyz(48),isx(48),isy(48),isz(48)
      real(8) :: Kc(3),g(3),Ktransform3D(3)
      
      Ktransform3D(1)=                Kc(1)
      Ktransform3D(2)=(1-iyz(igroup))*Kc(2) + iyz(igroup)*Kc(3)
      Ktransform3D(3)=(1-iyz(igroup))*Kc(3) + iyz(igroup)*Kc(2)
      
      g(1)=(1-ixz(igroup))*Ktransform3D(1) + ixz(igroup)*Ktransform3D(3)
      g(2)=                Ktransform3D(2)
      g(3)=(1-ixz(igroup))*Ktransform3D(3) + ixz(igroup)*Ktransform3D(1)
      
      Ktransform3D(1)=isx(igroup)*((1-ixy(igroup))*g(1) 
     &             + ixy(igroup)*g(2))
      Ktransform3D(2)=isy(igroup)*((1-ixy(igroup))*g(2) 
     &             + ixy(igroup)*g(1))
      Ktransform3D(3)=isz(igroup)*                 g(3) 
      
      return
      end
      
      function Inverse_3x3(a)   
      integer :: a(3,3)
      real(8) :: Inverse_3x3(3,3),r1
c
c       [ a(1,1) a(1,2) a(1,3) ]
c     A=[ a(2,1) a(2,2) a(2,3) ]
c       [ a(3,1) a(3,2) a(3,3) ]
c
c      -1        i+j 
c     A_ji = (-1)   C_ij  /|A|    C_ij is the ij cofactor of A
c     
      
      r1=   a(1,1)*a(2,2)*a(3,3)  ! the determinent of a
     &    + a(1,2)*a(2,3)*a(3,1)
     &    + a(1,3)*a(2,1)*a(3,2)
     &    - a(3,1)*a(2,2)*a(1,3)
     &    - a(3,2)*a(2,3)*a(1,1)
     &    - a(3,3)*a(2,1)*a(1,2) 
      
      Inverse_3x3(1,1) = +( a(2,2)*a(3,3)-a(3,2)*a(2,3) )/r1
      Inverse_3x3(1,2) = -( a(1,2)*a(3,3)-a(3,2)*a(1,3) )/r1
      Inverse_3x3(1,3) = +( a(1,2)*a(2,3)-a(2,2)*a(1,3) )/r1
      
      Inverse_3x3(2,1) = -( a(2,1)*a(3,3)-a(3,1)*a(2,3) )/r1
      Inverse_3x3(2,2) = +( a(1,1)*a(3,3)-a(3,1)*a(1,3) )/r1
      Inverse_3x3(2,3) = -( a(1,1)*a(2,3)-a(2,1)*a(1,3) )/r1
      
      Inverse_3x3(3,1) = +( a(2,1)*a(3,2)-a(3,1)*a(2,2) )/r1
      Inverse_3x3(3,2) = -( a(1,1)*a(3,2)-a(3,1)*a(1,2) )/r1
      Inverse_3x3(3,3) = +( a(1,1)*a(2,2)-a(2,1)*a(1,2) )/r1
      
      
      return
      end
      
      function Ideterminant(a)
      integer::  a(3,3),Ideterminant

      Ideterminant=
     &         a(1,1)*a(2,2)*a(3,3)
     &       + a(1,2)*a(2,3)*a(3,1)
     &       + a(1,3)*a(2,1)*a(3,2)
     &       - a(3,1)*a(2,2)*a(1,3)
     &       - a(3,2)*a(2,3)*a(1,1)
     &       - a(3,3)*a(2,1)*a(1,2)
      return
      end
      
      function Rdeterminant(a)
      real(8)::  a(3,3),Rdeterminant

      Rdeterminant=
     &         a(1,1)*a(2,2)*a(3,3)
     &       + a(1,2)*a(2,3)*a(3,1)
     &       + a(1,3)*a(2,1)*a(3,2)
     &       - a(3,1)*a(2,2)*a(1,3)
     &       - a(3,2)*a(2,3)*a(1,1)
     &       - a(3,3)*a(2,1)*a(1,2)
      return
      end




