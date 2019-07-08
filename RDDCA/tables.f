         subroutine tables
c*******************************************************************************
c        This subroutine makes the lookup tables which are used
c        throughout the rest of the codes.  These tables include:
c        wn(n)
c        
c*******************************************************************************
         use Global
         implicit none
c*******************************************************************************
        integer i,i1,i2,ix,iy,iodd,ilow,igroup,irwgroup(Nc),ind,
     &          ic,ic1,ic2,ick,idim,ik,itest,iex(8),isx(8),isy(8),
     &          j,jmd,jm,j1,j2,jc,jk,jkc,jhigh,
     &          k,n,n1,n2,m,m1,m2
        real(kind) :: Kx,Ky,r1,r2,Rx,Ry,vx,vy,origx,origy,Dkt
c*******************************************************************************
        
        N_sp=sqrt(float(Nc))

       
c*******************************************************************************
c       MATSUBARA FREQUENCIES       
c*******************************************************************************
        do n=-nwn,nwn
          wn(n)=(2.0_kind*n+1.0_kind)*pi*temp
        end do

c*******************************************************************************
c       Determine the princ. vectors a(3,3) and g(3,3)
c*******************************************************************************
      

c                                                                   side   clust
c               Different Real-space Cluster Geometries             trans  size     
cc---c---c---c---c---c---c---c-*-c---c---c---c---c---c---c---c---c  x  y   Nc
c| ********* |   |   |   |   *   *   |   |   |   |   |   |   |   |  2  0    4
cc-*-c---c-*-c---c---c---c-*-c---c-*-c---c---c---c---c---c-*-c---c  2  2    8  
c| * | 4 | * |   |   |   *   |   |   *   |   |   |   *   |   |   |  3  1   10
cc-*-c---c-*-c---c---c-*-c---c-8-c---c-*-c---c-*-c---c---c---c---c  4  0   16
c| ********* |   *   |   *   |   |   *   |   |   |   |   |   *   |  3  3   18
cc---c---c---c-*-c-*-c---c-*-c---c-*-c---c---c---c---c-10c---c---c  4  2   20
c|   |   |   *   2   *   |   *   *   |   |   |   *   |   |   |   |  5  1   26
cc---c---c---c-*-c-*-c---c---c-*-c---c---c---c---c---c---c---c-*-c  4  4   32
c|   |   |   |   *   |   |   |   |   |   |   |   |   |   *   |   |  5  3   34
cc---c---c---c---c---c---c---c---c---c---c---c---c-*-c---c---c---c  6  0   36
c|   |   |   |   *   |   |   |   |   |   |   |   |   |   |   |   |  6  2   40
cc---c---c---c-*-c-*-c---c---c---c---c---c---c---c---c---c---c---c  5  5   50
c|   |   |   *   |   *   |   |   |   |   | * * * * * * * * * |   |  7  1   50
cc---c---c-*-c---c---c-*-c---c---c---c---c-*-c---c---c---c-*-c---c  6  4   52
c|   |   *   |   |   |   *   |   |   |   | * |   |   |   | * |   |  7  3   58
cc---c-*-c---c---c---c---c-*-c---c---c---c-*-c---c---c---c-*-c---c  8  0   64
c|   *   |   |  18   |   |   *   |   |   | * |   |16 |   | * |   |  8  2   68
cc---c-*-c---c---c---c---c-*-c---c---c---c-*-c---c---c---c-*-c---c  6  6   72
c|   |   *   |   |   |   *   |   |   |   | * |   |   |   | * |   |  7  5   74
cc---c---c-*-c---c---c-*-c---c---c---c---c-*-c---c---c---c-*-c---c  8  4   80
c|   |   |   *   |   *   |   |   |   |   | * * * * * * * * * |   |  9  1   82
cc---c---c---c-*-c-*-c---c---c---c---c---c---c---c---c---c---c---c  9  3   90
c|   |   |   |   *   |   |   |   |   |   |   |   |   |   |   |   |  7  7   98
cc---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c  8  6  100
c                            _
c    _     _    _       2 pi a_i         a1 --> a2 [ x -->  y
c    a --> g    g_i =  -----------                 |
c                      | a1 X a2 |       g1 --> g2 [ y --> -x
c                                           
c       Principle translation                      [ square   (-1/2,-1/2)
c       vectors for lattice tiling        Origins  | diamond  ( 0  ,-1/2)
c                                                  [ odd      (-eps,-eps)
c
c              ^               ^                    ^                  ^
c Nc=1  a1 = 1 x,       a2 = 1 y,         g1 = 2 pi x,       g2 = 2 pi y
c
c       _      ^        _      ^          _       ^          _       ^
c Nc=4  a1 = 2 x,       a2 = 2 y,         g1 = pi x,         g2 = pi y
c
c       _      ^     ^  _       ^     ^   _       ^   ^      _        ^   ^
c Nc=8  a1 = 2 x + 2 y, a2 = -2 x + 2 y,  g1 = pi(x + y)/2,  g2 = pi(-x + y)/2
c
c       _      ^   ^    _     ^     ^     _        ^   ^     _        ^    ^
c Nc=10 a1 = 3 x + y,   a2 = -x + 3 y,    g1 = pi(3x + y)/5, g2 = pi(-x + 3y)/5
c
c       _      ^        _      ^          _       ^          _       ^
c Nc=16 a1 = 4 x,       a2 = 4 y,         g1 = pi x/2,       g2 = pi y/2
c
c       _     ^    ^    _      ^    ^     _       ^   ^      _        ^   ^
c Nc=18 a1 = 3x + 3y,   a2 = -3x + 3y,    g1 = pi(x + y)/3,  g2 = pi(-x + y)/3
c
c  Group Operations for a square lattice:
c
c  x --> +/- x
c    \/
c    /\
c  y --> +/- y       .
c
c  Square and diamond clusters have all 8, odd clusters (Nc=10,20,26) have 
c  only the first 4 below
c                              rot. by pi
c  identity   rot. by pi/2     inversion    rot. by -pi/2
c  x -> x       x ->  y         x -> -x        x -> -y
c  y -> y       y -> -x         y -> -y        y -> x
c  1,0,1,1      2,1,1,-1        3,0,-1,-1      4,1,-1,1
c
c  ref. x       ref. y          ref. x=y    ref. x=y inversion
c  x ->  x      x -> -x         x -> y         x -> -y
c  y -> -y      y ->  y         y -> x         y -> -x
c  5,0,1,-1     6,0,-1,1        7,1,1,1        8,1,-1,-1
c   
c  pgroup(iop,iex,isx,isy) 
c  iex=0(1) dont (do) exchange x and y
c  isx      sign of x 
c  isy      sign of y
c
c  x --> isx*[(1-iex)*x + iex*y]
c  y --> isy*[(1-iex)*y + iex*x]
c
        iex(1)=0; isx(1)= 1; isy(1)= 1
        iex(2)=1; isx(2)= 1; isy(2)=-1
        iex(3)=0; isx(3)=-1; isy(3)=-1
        iex(4)=1; isx(4)=-1; isy(4)= 1
        iex(5)=0; isx(5)= 1; isy(5)=-1
        iex(6)=0; isx(6)=-1; isy(6)= 1
        iex(7)=1; isx(7)= 1; isy(7)= 1
        iex(8)=1; isx(8)=-1; isy(8)=-1

c***********************************************************************        
c       determine the principle translation vectors a1 and a2, and 
c       reciprocal space vectors g1 and g2
c***********************************************************************
        If(Nc.eq.1) then                              
          a(1,1)=1
          a(1,2)=0
          Ncw=1
        else if(Nc.eq.2) then                           
          a(1,1)=1                                       
          a(1,2)=1                                       
          Ncw=2
        else if(Nc.eq.4) then                           
          a(1,1)=2                                       
          a(1,2)=0                                       
          Ncw=3
        else if(Nc.eq.8) then                           
          a(1,1)=2                                 
          a(1,2)=2                                 
       else if(Nc.eq.10) then                          
          a(1,1)=3                                 
          a(1,2)=1                                 
          Ncw=4
        else if(Nc.eq.16) then                          
          a(1,1)=4                                 
          a(1,2)=0
          Ncw=6
        else if(Nc.eq.18) then
          a(1,1)=3
          a(1,2)=3
        else if(Nc.eq.20) then
          a(1,1)=4
          a(1,2)=2
        else if(Nc.eq.26) then
          a(1,1)=5
          a(1,2)=1
        else if(Nc.eq.32) then
          a(1,1)=4
          a(1,2)=4
        else if(Nc.eq.34) then
          a(1,1)=5
          a(1,2)=3
        else if(Nc.eq.36) then
          a(1,1)=6
          a(1,2)=0
        else if(Nc.eq.40) then
          a(1,1)=6
          a(1,2)=2
        else if(Nc.eq.50) then
          a(1,1)=5
          a(1,2)=5
        else if(Nc.eq.52) then
          a(1,1)=6
          a(1,2)=4
        else if(Nc.eq.58) then
          a(1,1)=7
          a(1,2)=3
        else if(Nc.eq.64) then
          a(1,1)=8
          a(1,2)=0
        else if(Nc.eq.72) then
          a(1,1)=6
          a(1,2)=6
        else if(Nc.eq.98) then
          a(1,1)=7
          a(1,2)=7
        else if(Nc.eq.100) then
          a(1,1)=10
          a(1,2)=0
        else if(Nc.eq.128) then
          a(1,1)=8
          a(1,2)=8
        else if(Nc.eq.144) then
          a(1,1)=12
          a(1,2)=0
        else if(Nc.eq.162) then
          a(1,1)=9
          a(1,2)=9
        else if(Nc.eq.196) then
          a(1,1)=14
          a(1,2)=0
        else if(Nc.eq.200) then
          a(1,1)=10
          a(1,2)=10
        else if(Nc.eq.242) then
          a(1,1)=11
          a(1,2)=11
        else if(Nc.eq.256) then
          a(1,1)=16
          a(1,2)=0
        else if(Nc.eq.288) then
          a(1,1)=12
          a(1,2)=12
        else if(Nc.eq.324) then
          a(1,1)=18
          a(1,2)=0
        else if(Nc.eq.338) then
          a(1,1)=13
          a(1,2)=13
        else if(Nc.eq.392) then
          a(1,1)=14
          a(1,2)=14
        else if(Nc.eq.400) then
          a(1,1)=20
          a(1,2)=0
        else if(Nc.eq.450) then
          a(1,1)=15
          a(1,2)=15
        else if(Nc.eq.484) then
          a(1,1)=22
          a(1,2)=0
        else if(Nc.eq.512) then
          a(1,1)=16
          a(1,2)=16
        else if(Nc.eq.576) then
          a(1,1)=24
          a(1,2)=0
        else if(Nc.eq.1024) then
          a(1,1)=32
          a(1,2)=0
        else
          if(myrank.eq.0) write(lud,*) 'Bad Nc, bad, bad'              !0
          stop
        end if    

        a(2,1)=-a(1,2)  ! a2 is perpendicular to a1
        a(2,2)= a(1,1)

c       Now calculate the principle translations in K-space
        g(1,1)=twor*pi*a(1,1)/real(a(1,1)*a(2,2)-a(2,1)*a(1,2),kind)
        g(1,2)=twor*pi*a(1,2)/real(a(1,1)*a(2,2)-a(2,1)*a(1,2),kind)
        g(2,1)=-g(1,2)
        g(2,2)= g(1,1)

        if(iprint.ge.5.and.myrank.eq.0) 
     &    write(lud,*) ' a1x,a1y,g1x,g1y',a(1,1),a(1,2),g(1,1),g(1,2)  !0
        if(iprint.ge.5.and.myrank.eq.0) 
     &    write(lud,*) ' a2x,a2y,g2x,g2y',a(2,1),a(2,2),g(2,1),g(2,2)  !0

c*************************************************************************
c       Set the location of the origin in real space, denote the symmetry
c       (high, iodd=0; low, iodd=1) and the number of group operations.
c*************************************************************************
        if(a(1,1).eq.a(1,2)) then         ! The cell is a diamond (Nc=2,8,18,32,50...)
          origx= zeror
          origy=-halfr
          iodd=0     ! denote a cluster of high symmetry
          ngroup=8   ! with all ngroup=8 symmetry operations
        else if(a(1,1)*a(1,2).eq.0) then  ! The cell is a square (Nc=4,16,36,64...)
          origx=-halfr
          origy=-halfr
          iodd=0     ! denote a cluster of high symmetry
          ngroup=8   ! with all ngroup=8 symmetry operations
        else                              ! The cell is odd (i.e. Nc=10,20,26,34,40,52,58...)
          origx=-eps
          origy=-eps
          iodd=1     ! denote a cluster of lower symmetry
          ngroup=4   ! with only ngroup=4 symmetry operations
        end if


c****************************************************************************       
c       calculate the table of cluster point locations Rc(i,ic), i=1,x or 2,y
c****************************************************************************       
c 
c  ^ y
c  |
cc-|-c---c-*-c---c---c---c               A point P is within the cluster 
c| | |   |   |   |   |   |               if a vector v from the box origin, 
cc-|-c---c---c---c---c---c               O, to the point satisfies
ca2^ * * * * * * * * |   |               
cc-|-c---c---c---c-*-c---c                    _   _
c| | |   P   |   | * |   |                    v . a1
cc-|-c--/c---c---c-*-c---c               0 < -------   < 1
c| | | /v|   |   | * |   |                   |a1.a1|
cc-|-c/--c---c---c-*-c---c 
c| | /   |   |   | * |   |                      and
cc-|/0---c---c---c-*-c---c                    _   _
c| O--------------->---------> x              v . a2
cc---c---c---c---ca1-c---c               0 < -------   < 1
c|   |   |   |   |   |   |                   |a2.a2|
cc---c---c---c---c---c---c

        if(iprint.ge.5.and.myrank.eq.0) write(lud,*) '  '               !0
        if(iprint.ge.5.and.myrank.eq.0) 
     &       write(lud,*) 'ic,Rcx(ic),Rcy(ic)'                          !0
        ic=0
        do j=-Nc,Nc
        do i=-Nc,Nc
          vx=real(i,kind)-origx
          vy=real(j,kind)-origy
          r1=(vx*a(1,1)+vy*a(1,2))/real(a(1,1)**2+a(1,2)**2,kind)
          r2=(vx*a(2,1)+vy*a(2,2))/real(a(2,1)**2+a(2,2)**2,kind)
          if(r1.lt.oner.and.r1.gt.zeror.and.
     &      r2.lt.oner.and.r2.gt.zeror) then
            ic=ic+1
            Rc(1,ic)=i           ! Rc, like i.e. a(1,1) is integer
            Rc(2,ic)=j
            if(iprint.ge.5.and.myrank.eq.0) 
     &           write(lud,*) ic,Rc(1,ic),Rc(2,ic)                      !0
          end if
        end do
        end do
        if(ic.ne.Nc.or.Rc(1,1).ne.0.or.Rc(2,1).ne.0) then
          if(myrank.eq.0) write(lud,*) 'Bug in Rc(ic) loop'             !0
          stop
        end if
        
c***********************************************************************        
c       create a table mapping points outside the cluster back into it.
c***********************************************************************        

        if(iprint.ge.5.and.myrank.eq.0) write(lud,*) '  '               !0
        if(iprint.ge.5.and.myrank.eq.0) 
     &       write(lud,*) 'Equivalency table r'                         !0
        if(iprint.ge.5.and.myrank.eq.0) 
     &       write(lud,*) '   ic     igroup   icrequ(ic,ieqiv)'         !0
        do ic=1,Nc
        do igroup=1,ngroup
c         Generate the equivalent points using the ngroup operations
c         x --> isx*[(1-iex)*x + iex*y]
c         y --> isy*[(1-iex)*y + iex*x]
          n1=isx(igroup)*((1-iex(igroup))*Rc(1,ic)+iex(igroup)*Rc(2,ic))
          m1=isy(igroup)*((1-iex(igroup))*Rc(2,ic)+iex(igroup)*Rc(1,ic))
          do n=1,2  ! map (n1,m1) back into the cluster.
            vx=real(n1,kind)-origx
            vy=real(m1,kind)-origy
            r1=(vx*a(1,1)+vy*a(1,2))/
     &         real(a(1,1)*a(1,1)+a(1,2)*a(1,2),kind)
            r2=(vx*a(2,1)+vy*a(2,2))/
     &         real(a(2,1)*a(2,1)+a(2,2)*a(2,2),kind)
            if(r1.lt.zeror) then
              n1=n1+a(1,1)
              m1=m1+a(1,2)
            else if(r1.gt.oner) then
              n1=n1-a(1,1)
              m1=m1-a(1,2)
            end if
            if(r2.lt.zeror) then
              n1=n1+a(2,1)
              m1=m1+a(2,2)
            else if(r2.gt.oner) then
              n1=n1-a(2,1)
              m1=m1-a(2,2)
            end if
          end do
c         Figure out which point (n1,m1) is
          itest=0
          do jc=1,Nc
            if(n1.eq.Rc(1,jc).and.m1.eq.Rc(2,jc)) then
              icrequ(ic,igroup)=jc
              itest=1
            end if
          end do
          if(itest.ne.1) then
            write(lud,*) 'No equivalent found in icrequ',ic,igroup
            stop
          end if 
          if(iprint.ge.5.and.myrank.eq.0) 
     &        write(lud,*) ic,igroup,icrequ(ic,igroup)                !0
        end do
        end do                               


c***********************************************************************        
c       Now create the tables for the differences of R
c***********************************************************************        
        if(iprint.ge.5.and.myrank.eq.0) write(lud,*) '  '               !0
        if(iprint.ge.5.and.myrank.eq.0) 
     &       write(lud,*) 'differences in r'                            !0
        if(iprint.ge.5.and.myrank.eq.0) 
     &       write(lud,*) '  ic   jc  icrdiff(ic,jc)'                   !0
        do i=1,Nc
        do j=1,Nc
          Rx=real(Rc(1,i)-Rc(1,j),kind)
          Ry=real(Rc(2,i)-Rc(2,j),kind)
          itest=0
          do n=-1,1          ! Look in the surrounding nine
          do m=-1,1          ! clusters for Rx and Ry
            vx=Rx-origx + n*a(1,1) + m*a(2,1)
            vy=Ry-origy + n*a(1,2) + m*a(2,2)
            r1=(vx*a(1,1)+vy*a(1,2))/real(a(1,1)**2+a(1,2)**2,kind)
            r2=(vx*a(2,1)+vy*a(2,2))/real(a(2,1)**2+a(2,2)**2,kind)
            if(r1.lt.oner.and.r1.gt.zeror.and.
     &        r2.lt.oner.and.r2.gt.zeror) then
c             (Rx,Ry) is located within the cluster.  
c             Now we must find which point it is.
              do ic=1,Nc
                r1 = Rx + n*a(1,1) + m*a(2,1)
                r2 = Ry + n*a(1,2) + m*a(2,2)
                if(abs(r1-Rc(1,ic)).lt.eps.and.
     &             abs(r2-Rc(2,ic)).lt.eps) then
                   icrdiff(i,j)=ic
                   itest=1
                end if
              end do
              if(iprint.ge.5.and.myrank.eq.0) 
     &             write(lud,*) i,j,icrdiff(i,j)                        !0
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
c       Now set up the K-space tables
c***********************************************************************        
        
c                          K-space when Nc=8 
c                         -------------------
c
c       a(1,1)=2
c       a(1,2)=2
c
c       g(1,1)=twor*pi*a(1,1)/8  =  pi/2
c       g(1,2)=twor*pi*a(1,2)/8  =  pi/2
c       g(2,1)=-g(1,2)           = -pi/2
c       g(2,2)= g(1,1)           =  pi/2
c
c       _               ^   ^  
c       g1 = (pi/2) * ( x + y )
c
c       _               ^   ^
c       g2 = (pi/2) * (-x + y )
c
c                                ^
c                                | Ky     
cc---c---Q---c---c---c---c---c---c---c---c---c---Q---c---c---c---c   Legend 
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |  --------
cc---c---c---c---c---c---c---c---c---Q---c---c---c---c---c---c---c  K : in BZ
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |  K': in IRW
cc---c---c---c---c---c---Q---c---c---c---c---c---c---c---c---c---Q  Q : not in
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |  w : boundary
cc---c---c---Q*******************K'w w w w w w w w w K'--c---c---c     of IW
c|   |   |   * g |   |   |   |   w1,1|   |   |   | w |2,0|   |   |   
cQ---c---c---*---g---c---c---c---w---c---c---c---w---*---c---c---c 
c|   |   |   *   | g |   |   |   w   |   |   | w |   *   |   |   |  K and K' are
cc---c---c---*---c---g---c---c---w---c---c---w---c---*---c---c---c  indexed i,j
c|   |   |   *   |   | K |   |   w   |   | K'|   |   *   |   |   |  
cc---c---c---*---c---c---g---c---w---c---w-4-c---c---*---Q---c---c  Ncw=4
c|   |   |   *   |   |   | g |   w   | w |   |1,0|   *   |   |   |
cc---Q---c---*---c---c---c---g---w---w---c---c---c---*---c---c---c
c|   |   |   *   |   |   |   | g w w |   |   |   |   *   |   |   |
cc---c---c---*---c---c---c---c---K'--c---c---c---c---K---c---c---c--> Kx
c|   |   |   *   |   |   |   |   |0,0|   |   |   |   2   |   |   |
cc---c---c---*---c---c---c---c---c---c---c---c---c---*---c---Q---c   ilow =0
c|   |   |   *   |   |   |   |   |   |   |   |   |   *   |   |   |   ihigh=Nc
cc---c---Q---*---c---c---c---c---c---c---c---c---c---*---c---c---c   jlow =0
c|   |   |   *   |   | K |   |   |   |   | K |   |   *   |   |   |   jhigh=i
cc---c---c---*---c---c---c---c---c---c---c---c---c---*---c---c---c
c|   |   |   *   |   |   |   |   |   |   |   |   |   *   |   |   |
cc---c---c---*---c---c---c---c---c---c---c---c---c---*---c---c---Q
c|   |   |   *   |   |   |   |   |   |   |   |   |   *   |   |   |
cc---c---c---Q***************************************Q---c---c---c
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
cQ---c---c---c---c---c---c---c---c---c---Q---c---c---c---c---c---c
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
cc---c---c---c---c---c---c---Q---c---c---c---c---c---c---c---c---c
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
cc---c---c---c---Q---c---c---c---c---c---c---c---c---c---Q---c---c
c
cc
c                          K-space when Nc=10 
c                         -------------------
c
c       a(1,1)=3
c       a(1,2)=1
c
c       g(1,1)=twor*pi*a(1,1)/10 = 3 pi/5
c       g(1,2)=twor*pi*a(1,2)/10 =   pi/5
c       g(2,1)=-g(1,2)           =  -pi/5
c       g(2,2)= g(1,1)           = 3 pi/5
c
c       _        ^   ^   
c       g1 = pi(3x + y)/5
c       
c       _        ^    ^
c       g2 = pi(-x + 3y)/5
c                                ^
c                                | Ky     
cc---c---Q---c---c---c---c---c---c---c---c---c---Q---c---c---c---c   Legend 
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |  --------
cc---c---c---c---c---c---c---c---c---Q---c---c---c---c---c---c---c  K : in BZ
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |  K': in IRW
cc---c---c---c---c---c---Q---c---c---c---c---c---c---c---c---c---Q  Q : not in
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
cc---c---c---Q***************************************K'--c---c---c  K and K' are
c|   |   |   *   |   |   |   |   |   |   |   |   |   *2,1|   |   |  indexed i,j
cQ---c---c---*---c---c---c---c---c---c---K'--c---c---*---c---c---c  
c|   |   |   *   |   |   |   |   |   |   |1,1|   |   *   |   |   |  Ncw=4
cc---c---c---*---c---c---c---K---c---c---c---c---c---*---c---c---c
c|   |   |   *   |   |   |   |0,1|   |   |   |   |   *   |   |   |
cc---c---c---*---K---c---c---c---c---c---c---c---c---*---Q---c---c
c|   |   |   * -1,1  |   |   | g2|   |   |   |   |   *   |   |   |
cc---Q---c---*---c---c---c---c---c---c---c---K'--c---*---c---c---c
c|   |   |   *   |   |   |   |   |   | g1|   |1,0|   *   |   |   |
cc---c---c---*---c---c---c---c---K'--c---c---c---c---*---c---c---c--> Kx
c|   |   |   *   |   |   |   |   |0,0|   |   |   |   *   |   |   |
cc---c---c---*---c---K---c---c---c---c---c---c---c---*---c---Q---c   ilow =0 or 1
c|   |   |   *   | -1,0  |   |   |   |   |   |   |   *   |   |   |   ihigh=Nc
cc---c---Q---*---c---c---c---c---c---c---c---c---K---*---c---c---c   jlow =0
c|   |   |   *   |   |   |   |   |   |   |   |  1,-1 *   |   |   |   jhigh=Nc
cc---c---c---*---c---c---c---c---c---K---c---c---c---*---c---c---c
c|   |   |   *   |   |   |   |   |  0,-1 |   |   |   *   |   |   |
cc---c---c---*---c---c---K---c---c---c---c---c---c---*---c---c---Q
c|   |   |   *   |   | -1,-1 |   |   |   |   |   |   *   |   |   |
cc---c---c---Q***************************************Q---c---c---c
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
cQ---c---c---c---c---c---c---c---c---c---Q---c---c---c---c---c---c
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
cc---c---c---c---c---c---c---Q---c---c---c---c---c---c---c---c---c
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
cc---c---c---c---Q---c---c---c---c---c---c---c---c---c---Q---c---c
c
c
c                          K-space when Nc=16
c                         -------------------
c
c       a(1,1)=4
c       a(1,2)=0
c
c       g(1,1)=twor*pi*a(1,1)/16 = pi/2
c       g(1,2)=twor*pi*a(1,2)/16 = 0
c       g(2,1)=-g(1,2)           = 0
c       g(2,2)= g(1,1)           = pi/2
c
c       _             ^ 
c       g1 = (pi/2) * x 
c
c       _             ^
c       g2 = (pi/2) * y 
c
c                                ^
c                                | Ky     
cc---c---Q---c---c---Q---c---c---Q---c---c---Q---c---c---Q---c---c   Legend 
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |  --------
cc---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c  K': in IRW
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |  K : in 1BZ
cc---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c  Q : not in
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |  w : boundary
cc---c---Q*******c***K***c*******K*******c***K***c*******K'--c---c      of IRW
c|   |   *   |   |   |   |   |   |   |   |   |   |   | w w2,2|   | * : 1BZ of 
cc---c---*---c---c---c---c---c---c---c---c---c---c---c---c---c---c     boumdary
c|   |   *   |   |   |   |   |   |   |   |   |   | w |   w   |   |  
cc---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c
c|   |   *   |   |   |   |   |   |   |   |   | w |   |   w   |   |  Ncw=6
cc---c---Q---c---c---K---c---c---K---c---c---K'--c---c---K'--c---c
c|   |   *   |   |   |   |   |   |   |   | w |1,1|   |   w2,1|   |
cc---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c
c|   |   *   |   |   |   |   |   |   | w |   |   |   |   w   |   |
cc---c---*---c---c---c---c---c---c---c---c---c---c---c---c---c---c
c|   |   *   |   |   |   |   |   | w |   |   |   |   |   w   |   |
cc---c---Q---c---c---K---c---c---K'w-c-w-c-w-K'w-c-w-c-w-K'--c---c--> Kx
c|   |   *   |   |   |   |   |   |0,0|   |   |1,0|   |   *2,0|   |
cc---c---*---c---c---c---c---c---c---c---c---c---c---c---*---c---c   ilow =0 
c|   |   *   |   |   |   |   |   |   |   |   |   |   |   *   |   |   ihigh=Nc
cc---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c   jlow =0
c|   |   *   |   |   |   |   |   |   |   |   |   |   |   *   |   |   jhigh=i
cc---c---Q---c---c---K---c---c---K---c---c---K---c---c---K---c---c
c|   |   *   |   |   |   |   |   |   |   |   |   |   |   *   |   |
cc---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c
c|   |   *   |   |   |   |   |   |   |   |   |   |   |   *   |   |
cc---c---*---c---c---c---c---c---c---c---c---c---c---c---*---c---c
c|   |   *   |   |   |   |   |   |   |   |   |   |   |   *   |   |
cc---c---Q*******c***Q***c*******Q*******c***Q***c*******Q---c---c
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
cc---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
cc---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
cc---c---Q---c---c---Q---c---c---Q---c---c---Q---c---c---Q---c---c
c
c
c                          K-space when Nc=18
c                         -------------------
c
c       a(1,1)=3
c       a(1,2)=3
c
c       g(1,1)=twor*pi*a(1,1)/18 = pi/3
c       g(1,2)=twor*pi*a(1,2)/18 = pi/3
c       g(2,1)=-g(1,2)           = -pi/3
c       g(2,2)= g(1,1)           = pi/3
c
c       _               ^   ^  
c       g1 = (pi/3) * ( x + y )
c
c       _               ^   ^
c       g2 = (pi/3) * (-x + y )
c
c                                ^
c                                | Ky     
cc---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c   Legend 
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |  --------
cc---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c  K': in IRW
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |  K : in 1BZ
cc---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c  Q : not in
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |  w : boundary
cc---c---Q*******c*******K*******w w w w K'w w w w w w w K'--c---c      of IRW
c|   |   * g |   |   |   |   |   w   |   |2,1|   |   | w |3,0|   |  * : 1BZ of 
cc---c---*---g---g---c---c---c---w---c---c---c---c---g1--*---c---c      boumdary
c|   |   *   | g |   |   |   |   w   |   |   |   | w |   *   |   |   
cc---c---c---c---K---c---c---c---K'--c---c---c---K'--c---c---c---c
c|   |   *   |   | g |   |   |   w1,1|   |   | w |2,0|   *   |   |  Ncw=6
cc---c---*---c---c---g---c---c---w---c---c---w---c---c---*---c---c
c|   |   *   |   |   | g |   |   w   |   | w |   |   |   *   |   |
cc---c---Q---c---c---c---K---c---w---c---K'--c---c---c---K---c---c
c|   |   *   |   |   |   | g |   w   | w |1,0|   |   |   *   |   |
cc---c---*---c---c---c---c---g---w---w---c---c---c---c---*---c---c
c|   |   *   |   |   |   |   | g w w |   |   |   |   |   *   |   |
cc---c---c---c---K---c---c---c---K'--c---c---c---K---c---c---c---c--> Kx
c|   |   *   |   |   |   |   |   |0,0|   |   |   |   |   *   |   |
cc---c---*---c---c---c---c---c---c---c---c---c---c---c---*---c---c   ilow =0 
c|   |   *   |   |   |   |   |   |   |   |   |   |   |   *   |   |   ihigh=Nc
cc---c---Q---c---c---c---K---c---c---c---K---c---c---c---K---c---c   jlow =0
c|   |   *   |   |   |   |   |   |   |   |   |   |   |   *   |   |   jhigh=i
cc---c---*---c---c---c---c---c---c---c---c---c---c---c---*---c---c
c|   |   *   |   |   |   |   |   |   |   |   |   |   |   *   |   |
cc---c---c---c---K---c---c---c---K---c---c---c---K---c---c---c---c
c|   |   *   |   |   |   |   |   |   |   |   |   |   |   *   |   |
cc---c---*---c---c---c---c---c---c---c---c---c---c---c---*---c---c
c|   |   *   |   |   |   |   |   |   |   |   |   |   |   *   |   |
cc---c---Q*******c*******Q*******c*******Q*******c*******K---c---c
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
cc---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
cc---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
cc---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c
c
c                          K-space when Nc=20 
c                         -------------------
c       a(1,1)=4
c       a(1,2)=2
c
c       g(1,1)=twor*pi*a(1,1)/20 = 2 pi/5
c       g(1,2)=twor*pi*a(1,2)/20 =   pi/5
c       g(2,1)=-g(1,2)           = - pi/5
c       g(2,2)= g(1,1)           = 2 pi/5
c
c
c       _        ^   ^    
c       g1 = pi(2x + y)/5
c
c       _        ^    ^
c       g2 = pi(-x + 2y)/5
c
c
c                                ^
c                                | Ky     
cc---c---Q---c---c---c---c---c---c---c---c---c---Q---c---c---c---c   Legend 
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |  --------
cc---c---c---c---c---c---c---c---c---Q---c---c---c---c---c---c---c  K : in BZ
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |  K': In IRW
cc---c---c---c---c---c---Q---c---c---c---c---c---c---c---c---c---Q  Q : not in BZ
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
cc---c---c---Q*******************K'******************K'--c---c---c  K and K' are
c|   |   |   *   |   |   |   |   |1,2|   |   |   |   *3,1|   |   |  indexed i,j
cQ---c---c---*---c---c---K---c---c---c---c---K'--c---*---c---c---c  
c|   |   |   *   |   |   |0,2|   |   |   |   |2,1|   *   |   |   |  Ncw=7
cc---c---c---*---K---c---c---c---c---K'--c---c---c---*---c---c---c
c|   |   |   * -1,2  |   |   |   |   |1,1|   |   |   *   |   |   |
cc---c---c---*---c---c---c---K---c---c---c---c---K'--*---Q---c---c
c|   |   |   *   |   |   |   |0,1|   |   |   |   |2,0*   |   |   |
cc---Q---c---*---c---K---c---c---c---c---K'--c---c---*---c---c---c
c|   |   |   *   | -1,1  |   |   |   |   |1,0|   |   *   |   |   |
cc---c---c---*---c---c---c---c---K'--c---c---c---c---K'--c---c---c--> Kx
c|   |   |   *   |   |   |   |   |0,0|   |   |   |  2,-1 |   |   |
cc---c---c---*---c---c---K---c---c---c---c---K---c---*---c---Q---c   ilow =0 or 1
c|   |   |   *   |   | -1,0  |   |   |   |  1,-1 |   *   |   |   |   ihigh=Nc
cc---c---Q---*---K---c---c---c---c---K---c---c---c---*---c---c---c   jlow =0
c|   |   |   * -2,0  |   |   |   |  0,-1 |   |   |   *   |   |   |   jhigh= Nc
cc---c---c---*---c---c---c---K---c---c---c---c---K---*---c---c---c
c|   |   |   *   |   |   | -1,-1 |   |   |   |  1,-2 *   |   |   |
cc---c---c---*---c---K---c---c---c---c---K---c---c---*---c---c---Q
c|   |   |   *   | -2,-1 |   |   |   |  0,-2 |   |   *   |   |   |
cc---c---c---Q***************************************Q---c---c---c
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
cQ---c---c---c---c---c---c---c---c---c---Q---c---c---c---c---c---c
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
cc---c---c---c---c---c---c---Q---c---c---c---c---c---c---c---c---c
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
cc---c---c---c---Q---c---c---c---c---c---c---c---c---c---Q---c---c
c
c
c                          K-space when Nc=26
c                         -------------------
c
c       a(1,1)=5
c       a(1,2)=1
c
c       g(1,1)=twor*pi*a(1,1)/26 = 5*pi/13
c       g(1,2)=twor*pi*a(1,2)/26 =   pi/13
c       g(2,1)=-g(1,2)           =  -pi/13
c       g(2,2)= g(1,1)           = 5*pi/13
c
c       _                  ^   ^  
c       g1 = (pi/13) * ( 5 x + y )
c
c       _                ^     ^
c       g2 = (pi/13) * (-x + 5 y )
c
c                                ^
c                                | Ky     
cc---c---c---c---c---c---c---c---c---c---c---c---c---c-Q-c---c---c   Legend 
c|   |   |   |   |   |   |   |   |   |   |   | Q |   |   |   |   |  --------
cc---c---c---c---c---c---c---c---c---Q---c---c---c---c---c---c---c  K': in IRW
c|   |   |   |   |   |   | Q |   |   |   |   |   |   |   |   |   |  K : in 1BZ
cc---c---c---c---Q---c---c---c---c---c---c---c---c---c---c---c---c  Q : not in
c|   | Q*|***|***|***|***|***|***|***|***|***|***|***|***| K'|   |  * : 1BZ of 
cc---c-*-c---c---c---c---c---c---c---c---c---c---K'--c---c-*-c---c      boumdary
c|   | * |   |   |   |   |   |   |   | K'|   |   |2,2|   |3,2|   |  
cc---c-*-c---c---c---c---c---K---c---c---c---c---c---c---c-*-c---c
c|   | * |   |   | K |   |   |   |   |1,2|   |   |   |   | * |   |                 
cc---c-*-K---c---c---c---c---c---c---c---c---c---c---c---c-*-Q---c  K and K' are 
c|   | * |   |   |   |   |   |   |   |   |   |   | K'|   | * |   |  indexed i,j  
cc---c-*-c---c---c---c---c---c---c---c---K'--c---c---c---c-*-c---c  
c|   | * |   |   |   |   |   | K |   |   |1,1|   |2,1|   | * |   |  Ncw=8
cc---c-*-c---c---c---K---c---c---c---c---c---c---c---c---c-*-c---c
c|   | * | K |   |   |   |   |   |   |   |   |   |   |   | * | Q |
cc---c-*-c---c---c---c---c---c---c---c---c---c---c---K'--c-*-c---c
c|   | * |   |   |   |   |   |   |   |   | K'|   |   |2,0| * |   |
cc---c-*-c---c---c---c---c---c---K'--c---c---c---c---c---c-*-c---c--> Kx
c|   | * |   |   |   | K |   |   |0,0|   |1,0|   |   |   | * |   |
cc---c-*-c---K---c---c---c---c---c---c---c---c---c---c---c-*-c---Q   ilow =0 or 1
c|   | * |   |   |   |   |   |   |   |   |   |   |   | K | * |   |   ihigh=Nc
cc---c-*-c---c---c---c---c---c---c---c---c---K---c---c---c-*-c---c   jlow =0
c|   | * |   |   |   |   |   |   | K |   |   |   |   |   | * |   |   jhigh=Nc
cc---c-*-c---c---c---c---K---c---c---c---c---c---c---c---c-*-c---c
c|   | * |   | K |   |   |   |   |   |   |   |   |   |   | * |   |
cc---c-*-c---c---c---c---c---c---c---c---c---c---c---c---K-*-c---c
c|   | * |   |   |   |   |   |   |   |   |   | K |   |   | * |   |
cc---c-*-c---c---c---c---c---c---c---K---c---c---c---c---c-*-c---c
c|   | * |   |   |   |   | K |   |   |   |   |   |   |   | * |   |
cc---c-*-c---c---K---c---c---c---c---c---c---c---c---c---c-*-c---c
c|   | Q*|***|***|***|***|***|***|***|***|***|***|***|***|*Q |   |
cc---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
cc---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
cc---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c
c
c
c***********************************************************************        
c       Create the cluster momentum tables Kc(i,ic), i=1,x 2,y
c***********************************************************************        

c       First find the points in the irreducible wedge.  They will be
c       indexed from 1 to Ncw
        if(iprint.ge.5.and.myrank.eq.0) 
     &       write(lud,*) '  '                                           !0
        if(iprint.ge.5.and.myrank.eq.0) 
     &       write(lud,*) ' ic  i  j   Kcx(ic)   Kcy(ic)  '              !0
        ic=0
        Ncw=0
        ic0=1
        icpi=1

        if(iodd.eq.0) then   ! The cluster is of high symmetry
          do i=0,Nc          ! the wedge is roughly 1/8 of the 1BZ.
          do j=0,i           ! Index j runs only from 0 to i.
            Kx=i*g(1,1) + j*g(2,1)
            Ky=i*g(1,2) + j*g(2,2)
            if(Kx.gt.-pi+eps.and.Kx.lt.pi+eps.and.
     &        Ky.gt.-pi+eps.and.Ky.lt.pi+eps) then
c             Kx and Ky fall within the first Brillouin zone
              ic=ic+1
              Ncw=Ncw+1
              Kc(1,ic)=Kx
              Kc(2,ic)=Ky
              ig(ic)=i
              jg(ic)=j
              if(abs(Kx).lt.eps.and.abs(Ky).lt.eps) ic0=ic
              if(abs(Kx-pi).lt.eps.and.abs(Ky-pi).lt.eps) icpi=ic
              if(iprint.ge.5.and.myrank.eq.0) 
     &            write(lud,5) ic,i,j,Kc(1,ic),Kc(2,ic)                    !0
            end if        
          end do
          end do
        else if(iodd.eq.1) then ! the cluster is of low symmetry, so
          do j=0,Nc             ! the wedge is roughly 1/4 of the 1BZ
          if(j.eq.0) then       ! Include the point at K=(0,0)
            ilow=0
          else                  ! exclude the rest of the points along 
            ilow=1              ! the g2 axis (i=0).
          end if
          do i=ilow,Nc
            Kx=i*g(1,1) + j*g(2,1)
            Ky=i*g(1,2) + j*g(2,2)
            if(Kx.gt.-pi+eps.and.Kx.lt.pi+eps.and.
     &        Ky.gt.-pi+eps.and.Ky.lt.pi+eps) then
c             Kx and Ky fall within the first zone
              ic=ic+1
              Ncw=Ncw+1
              Kc(1,ic)=Kx
              Kc(2,ic)=Ky
              ig(ic)=i
              jg(ic)=j
              if(abs(Kx).lt.eps.and.abs(Ky).lt.eps) ic0=ic
              if(abs(Kx-pi).lt.eps.and.abs(Ky-pi).lt.eps) icpi=ic
              if(iprint.ge.5.and.myrank.eq.0) 
     &            write(lud,5) ic,i,j,Kc(1,ic),Kc(2,ic)                    !0
            end if        
          end do
          end do
        else
          if(myrank.eq.0) write(lud,*) 'woops, iodd=',iodd
          stop
        end if
 5      format(1x,i2,1x,i2,1x,i2,2x,f9.6,2x,f9.6)
        
        ntw=nl*Ncw
        if(iprint.ge.5.and.myrank.eq.0) 
     &       write(lud,*) 'Nc= ',Nc,' Ncw=',Ncw                          !0
        
c       Now find the remaining points within the 1Bz and outside the 
c       wedge.  These points are indexed from Ncw+1 to Nc       
        do i=-Nc,Nc
        do j=-Nc,Nc
          Kx=i*g(1,1) + j*g(2,1)
          Ky=i*g(1,2) + j*g(2,2)
          itest=0
          do ick=1,Ncw
            if(abs(Kx-Kc(1,ick)).lt.eps
     &         .and.
     &         abs(Ky-Kc(2,ick)).lt.eps)then
c              This K is within the wedge angle
               itest=1
            end if
          end do
          if(itest.eq.0.and.Kx.gt.-pi+eps.and.Kx.lt.pi+eps.and.
     &      Ky.gt.-pi+eps.and.Ky.lt.pi+eps) then
c           Kx and Ky fall within the first zone and are not IRW
            ic=ic+1
            Kc(1,ic)=Kx
            Kc(2,ic)=Ky
            ig(ic)=i
            jg(ic)=j
            if(iprint.ge.5.and.myrank.eq.0) 
     &          write(lud,5) ic,i,j,Kc(1,ic),Kc(2,ic)                    !0
          end if
        end do
        end do
        
        if(ic.ne.Nc) then
          if(myrank.eq.0) write(lud,*) 'Bug in Kc(ic) loop'              !0
          stop
        end if
        if(iprint.ge.5.and.myrank.eq.0) 
     &       write(lud,*) 'ic0=',ic0                                     !0
        if(iprint.ge.5.and.myrank.eq.0) 
     &       write(lud,*) 'icpi=',icpi                                   !0
        
        
        
c***********************************************************************        
c       Create the indirect addressing table ict(i,j) which maps a 
c       point indexed by the principle translations in k-space to ic.
c***********************************************************************        
        do i=-N_sp,N_sp
        do j=-N_sp,N_sp
          Kx=i*g(1,1) + j*g(2,1)
          Ky=i*g(1,2) + j*g(2,2)
          if(Kx.lt.-pi+eps) Kx=Kx+twor*pi
          if(Kx.gt.pi+eps) Kx=Kx-twor*pi
          if(Ky.lt.-pi+eps) Ky=Ky+twor*pi
          if(Ky.gt.pi+eps) Ky=Ky-twor*pi
          do jc=1,Nc
            if(abs(Kx-Kc(1,jc)).lt.eps
     &         .and.
     &         abs(Ky-Kc(2,jc)).lt.eps) then    
              ict(i,j)=jc
            end if
          end do
        end do
        end do


c***********************************************************************        
c       Create the table ickdiff(ic1,ic2), which takes the difference
c       between to K-points
c***********************************************************************        
        if(iprint.ge.5.and.myrank.eq.0) write(lud,*) '   '               !0
        if(iprint.ge.5.and.myrank.eq.0) 
     &       write(lud,*) ' ic1    ic2  ickdiff(ic1,ic2)  '              !0
        do ic1=1,Nc
        do ic2=1,Nc
          Kx=Kc(1,ic1)-Kc(1,ic2)
          Ky=Kc(2,ic1)-Kc(2,ic2)
          if(Kx.lt.-pi+eps) Kx=Kx+twor*pi
          if(Kx.gt.pi+eps) Kx=Kx-twor*pi
          if(Ky.lt.-pi+eps) Ky=Ky+twor*pi
          if(Ky.gt.pi+eps) Ky=Ky-twor*pi
c         Kdiffx(ic1,ic2)=Kx
c         Kdiffy(ic1,ic2)=Ky
c         Now we need to figure out where this point is!
          do ic=1,nc
            if(abs(Kx-Kc(1,ic)).lt.eps
     &         .and.
     &         abs(Ky-Kc(2,ic)).lt.eps) then    
              ickdiff(ic1,ic2)=ic
            end if
          end do
          if(iprint.ge.5.and.myrank.eq.0) 
     &       write(lud,*) ic1,ic2,ickdiff(ic1,ic2)                       !0
        end do
        end do

c***********************************************************************        
c       Create the table ickplus(ic1,ic2), which takes the sum
c       between to K-points
c***********************************************************************        
        if(iprint.ge.5.and.myrank.eq.0) write(lud,*) '   '               !0
        if(iprint.ge.5.and.myrank.eq.0) 
     &       write(lud,*) ' ic1    ic2  ickplus(ic1,ic2)  '              !0
        do ic1=1,Nc
        do ic2=1,Nc
          Kx=Kc(1,ic1)+Kc(1,ic2)
          Ky=Kc(2,ic1)+Kc(2,ic2)
          if(Kx.lt.-pi+eps) Kx=Kx+twor*pi
          if(Kx.gt.pi+eps) Kx=Kx-twor*pi
          if(Ky.lt.-pi+eps) Ky=Ky+twor*pi
          if(Ky.gt.pi+eps) Ky=Ky-twor*pi
c         Kplusx(ic1,ic2)=Kx
c         Kplusy(ic1,ic2)=Ky
c         Now we need to figure out where this point is!
          do ic=1,nc
            if(abs(Kx-Kc(1,ic)).lt.eps.and.abs(Ky-Kc(2,ic)).lt.eps) then
              ickplus(ic1,ic2)=ic
            end if
          end do
          if(iprint.ge.5.and.myrank.eq.0) 
     &       write(lud,*) ic1,ic2,ickplus(ic1,ic2)                       !0
        end do
        end do
        
c***********************************************************************                
c       Now find the points equivalent to ic, ickequ(ic,igroup) where 
c       igroup is the group operation.
c***********************************************************************        
c
c
        if(iprint.ge.5.and.myrank.eq.0) write(lud,*) '   '               !0
        if(iprint.ge.5.and.myrank.eq.0) write(lud,*) 'equivalency in k'  !0
        if(iprint.ge.5.and.myrank.eq.0) 
     &              write(lud,*) '   ic     igroup   ickequ(ic,igroup)'  !0
        do ic=1,Nc
          i=ig(ic)
          j=jg(ic)
c         Kx and Ky fall within the zone, by construction, now look 
c         for equivalent points
c
c  n --> isx*[(1-iex)*i + iex*j]
c  m --> isy*[(1-iex)*j + iex*i]
c
          do igroup=1,ngroup
            n=isx(igroup)*((1-iex(igroup))*i + iex(igroup)*j)
            m=isy(igroup)*((1-iex(igroup))*j + iex(igroup)*i)
            Kx=n*g(1,1) + m*g(2,1)
            Ky=n*g(1,2) + m*g(2,2)
c           These are the equivalent points.  Now map to 1st BZ
            if(Kx.lt.-pi+eps) Kx=Kx+twor*pi
            if(Kx.gt.pi+eps) Kx=Kx-twor*pi
            if(Ky.lt.-pi+eps) Ky=Ky+twor*pi
            if(Ky.gt.pi+eps) Ky=Ky-twor*pi
c           Now we need to figure out where this point is!
            itest=0
            do jc=1,Nc
              if(abs(Kx-Kc(1,jc)).lt.eps
     &           .and.
     &           abs(Ky-Kc(2,jc)).lt.eps) then
                ickequ(ic,igroup)=jc
                itest=1
              end if
            end do
            if(iprint.ge.5.and.myrank.eq.0) 
     &           write(lud,*) ic,igroup,ickequ(ic,igroup)                !0
            if(itest.eq.0) then
              if(myrank.eq.0) write(lud,*) 'no equivalent K-point found' !0
              if(myrank.eq.0) write(lud,*) 'j1,n2,m2',j1,n2,m2           !0
              if(myrank.eq.0) write(lud,*) 'Kx,Ky',Kx,Ky                 !0
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
c       for example forNc=8, only 4 of the 8 points are in the irreducible 
c       wedge.  Thus, all storage arrays may be reduced in size by a 
c       factor of 2.  As Nc-->oo, the saving factor approaches 8 for
c       square and diamond clusters.  
c
c       However, for odd clusters, where ngroup=4, the savings is at
c       most 4.  For an example, see the K-space of the Nc=10 and 20
c       clusters pictured above.
c***********************************************************************        
c       

        ickdeg=0
c       scan over all of the K points in the 1BZ
        do ic=1,Nc
          i=ig(ic)
          j=jg(ic)
c         scan over all allowed point group operations, to find the
c         one that maps ic into the IRW
          itest=0
          groupops: do igroup=1,ngroup
            n=isx(igroup)*((1-iex(igroup))*i + iex(igroup)*j)
            m=isy(igroup)*((1-iex(igroup))*j + iex(igroup)*i)
            Kx=n*g(1,1) + m*g(2,1)
            Ky=n*g(1,2) + m*g(2,2)
c           These are the equivalent points.  Now map to 1st BZ
            if(Kx.lt.-pi+eps) Kx=Kx+twor*pi
            if(Kx.gt.pi+eps) Kx=Kx-twor*pi
            if(Ky.lt.-pi+eps) Ky=Ky+twor*pi
            if(Ky.gt.pi+eps) Ky=Ky-twor*pi
c           See if this point is in the IRW!
            do jc=1,Ncw
              if(abs(Kx-Kc(1,jc)).lt.eps
     &           .and.
     &           abs(Ky-Kc(2,jc)).lt.eps) then
                ickmap(ic)=jc
                ickdeg(jc)=ickdeg(jc)+1
                irwgroup(ic)=igroup
                itest=1
              end if
            end do
            if(itest.eq.1) exit ! at first group op. that takes ic into IRW
          end do groupops
          if(itest.eq.0) then
            if(myrank.eq.0) write(lud,*) 'mapping failed for'            !0
            if(myrank.eq.0) write(lud,*) 'ic',ic                         !0
            if(myrank.eq.0) write(lud,*) 'Kx,Ky',Kx,Ky                   !0
            stop
          end if
          
        end do

c       Test the degeneracy table.  Since ickdeg holds the degeneracy
c       of points in the IRW, its sum must equal Nc.
        i=0
        do ic=1,Ncw
          i=i+ickdeg(ic)
        end do
        if(i.ne.Nc) then
          if(myrank.eq.0) write(lud,*) 'bug in ickdeg',i
          stop
        end if
        
c       Test the group table.  irwgroup(ic within wedge)=1, the identity
        do ic=1,Ncw
          if(irwgroup(ic).ne.1) then
            if(myrank.eq.0) write(lud,*) 'irwgroup failed for wedge pt.' !0
            if(myrank.eq.0) write(lud,*) 'ic',ic                         !0
            if(myrank.eq.0) write(lud,*) 'Kx,Ky',Kx,Ky                   !0
          end if
        end do
       
        
c***********************************************************************        
c       Form Epsbar(K)
c***********************************************************************        
        Dkt=halfr/real(nover,kind)
        Epsbar=zeror
        do ic=1,Nc
          do i=-nover+1,nover
          do j=-nover+1,nover
            kx=Kc(1,ic) + Dkt*((real(i,kind)-halfr)*g(1,1) +
     &                       (real(j,kind)-halfr)*g(2,1))
            ky=Kc(2,ic) + Dkt*((real(i,kind)-halfr)*g(1,2) +
     &                       (real(j,kind)-halfr)*g(2,2))
            Epsbar(ic)=Epsbar(ic) -halfr*(cos(kx)+cos(ky)) -
     &                 tprime*(cos(kx)*cos(ky)-oner)
          end do
          end do
          Epsbar(ic)=Epsbar(ic)/real((2*nover)**2,kind)
        end do

        return
        end     
