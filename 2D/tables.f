         subroutine tables
c        This subroutine makes the lookup tables which are used in
c        the apl portion of the program.
c************************************************************
         use Global
         implicit none
c************************************************************
 
        integer :: i,j,k,n,m,m1,n1,n2,m2,jmd,jm,ind,iequiv,i1,j1,
     &          ic,ic1,ic2,jc,ix,jk,ig(Nc),jg(Nc),itest,
     &          isx(8),isy(8),iex(8)
        integer, parameter :: iprint=2
        real(kind) :: etu2,Kx,Ky,r1,r2,Rx,Ry,vx,vy
        complex(kind) :: c1,c2
        
        N_sp=sqrt(float(Nc))
        
 
c       MATSUBARA FREQUENCIES
        do n=-nwn,nwn
          wn(n)=(2.0_kind*n+1.0_kind)*pi*temp
        end do
        

c                                                                   side   clust
c               Different Real-space Cluster Geometries             trans  size     
cc---c---c---c---c---c---c---c-*-c---c---c---c---c---c---c---c---c  x  y   Nc
c| ********* |   |   |   |   *   *   |   |   |   |   |   |   |   |  2  0    4
cc-*-c---c-*-c---c---c---c-*-c---c-*-c---c---c---c---c---c-*-c---c  2  2    8  
c| * | 4 | * |   |   |   *   |   |   *   |   |   |   *   |   |   |  3  1   10
cc-*-c---c-*-c---c---c-*-c---c-8-c---c-*-c---c-*-c---c---c---c---c  4  0   16
c| ********* |   |   |   *   |   |   *   |   |   |   |   |   *   |  3  3   18
cc---c---c---c---c---c---c-*-c---c-*-c---c---c---c---c-10c---c---c  4  2   20
c|   |   |   |   |   |   |   *   *   |   |   |   *   |   |   |   |  5  1   26
cc---c---c---c---c---c---c---c-*-c---c---c---c---c---c---c---c-*-c  4  4   32
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   *   |   |  5  3   34
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
c       vectors for lattice tiling        Origins  | diamond  (0,-1/2)
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
c  Group Operations
c
c  x --> +/- x
c    \/
c    /\
c  y --> +/- y
c
c  square and diamond have all 8, odd clusters (Nc=10,20,26) have 
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

        ngroup=8
        If(Nc.eq.1) then                              
          a1x=1
          a1y=0
          Ncw=1
        else if(Nc.eq.2) then                           
          a1x=1                                       
          a1y=1                                       
          Ncw=2
        else if(Nc.eq.4) then                           
          a1x=2                                       
          a1y=0                                       
          Ncw=3
        else if(Nc.eq.8) then                           
          a1x=2                                 
          a1y=2                                 
       else if(Nc.eq.10) then                          
          a1x=3                                 
          a1y=1                                 
          Ncw=4
          ngroup=4   ! fix ngroup for lower symmetry clusters
        else if(Nc.eq.16) then                          
          a1x=4                                 
          a1y=0
          Ncw=6
        else if(Nc.eq.18) then
          a1x=3
          a1y=3
        else if(Nc.eq.20) then
          a1x=4
          a1y=2
          ngroup=4   ! fix ngroup for lower symmetry clusters
        else if(Nc.eq.26) then
          a1x=5
          a1y=1
          ngroup=4
        else if(Nc.eq.32) then
          a1x=4
          a1y=4
        else if(Nc.eq.34) then
          a1x=5
          a1y=3
          ngroup=4   ! fix ngroup for lower symmetry clusters
        else if(Nc.eq.36) then
          a1x=6
          a1y=0
        else if(Nc.eq.40) then
          a1x=6
          a1y=2
          ngroup=4   ! fix ngroup for lower symmetry clusters
        else if(Nc.eq.50) then
          a1x=5
          a1y=5
        else if(Nc.eq.52) then
          a1x=6
          a1y=4
          ngroup=4   ! fix ngroup for lower symmetry clusters
        else if(Nc.eq.58) then
          a1x=7
          a1y=3
          ngroup=4   ! fix ngroup for lower symmetry clusters
        else if(Nc.eq.64) then
          a1x=8
          a1y=0
        else if(Nc.eq.72) then
          a1x=6
          a1y=6
        else if(Nc.eq.98) then
          a1x=7
          a1y=7
        else if(Nc.eq.100) then
          a1x=10
          a1y=0
        else if(Nc.eq.128) then
          a1x=8
          a1y=8
        else if(Nc.eq.144) then
          a1x=12
          a1y=0
        else if(Nc.eq.162) then
          a1x=9
          a1y=9
        else if(Nc.eq.196) then
          a1x=14
          a1y=0
        else if(Nc.eq.200) then
          a1x=10
          a1y=10
        else if(Nc.eq.242) then
          a1x=11
          a1y=11
        else if(Nc.eq.256) then
          a1x=16
          a1y=0
        else if(Nc.eq.288) then
          a1x=12
          a1y=12
        else if(Nc.eq.324) then
          a1x=18
          a1y=0
        else if(Nc.eq.338) then
          a1x=13
          a1y=13
        else if(Nc.eq.392) then
          a1x=14
          a1y=14
        else if(Nc.eq.400) then
          a1x=20
          a1y=0
        else if(Nc.eq.450) then
          a1x=15
          a1y=15
        else if(Nc.eq.484) then
          a1x=22
          a1y=0
        else if(Nc.eq.512) then
          a1x=16
          a1y=16
        else if(Nc.eq.576) then
          a1x=24
          a1y=0
        else if(Nc.eq.1024) then
          a1x=32
          a1y=0
        else
          write(6,*) 'Bad Nc, bad, bad'                !0
          stop
        end if    

        a2x=-a1y
        a2y= a1x
        g1x=2.0_kind*pi*a1x/real(a1x*a2y-a2x*a1y,kind)
        g1y=2.0_kind*pi*a1y/real(a1x*a2y-a2x*a1y,kind)
        g2x=-g1y
        g2y= g1x
        if(a1x.eq.a1y) then             ! The cell is a diamond (Nc=2,8,18,32,50...)
          origx= 0.0_kind
          origy=-0.5_kind
        else if(a1x*a1y.eq.0) then      ! The cell is a square (Nc=4,16,36,64...)
          origx=-0.5_kind
          origy=-0.5_kind
        else                            ! The cell is odd (i.e. Nc=10,20,26,34,40,52,58...)
          origx=-eps
          origy=-eps
        end if
        
        if(iprint.ge.1) then
          write(6,*) ' a1x,a1y,g1x,g1y',a1x,a1y,g1x,g1y   
          write(6,*) ' a2x,a2y,g2x,g2y',a2x,a2y,g2x,g2y
          write(6,*) ' origx,origy',origx,origy
        end if


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
c| | /   |   |   | * |   |               
cc-|/0---c---c---c-*-c---c                    _   _
c| O--------------->---------> x              v . a2
cc---c---c---c---ca1-c---c               0 < -------   < 1
c|   |   |   |   |   |   |                   |a2.a2|
cc---c---c---c---c---c---c

        if(iprint.gt.1) write(6,*) '  '
        if(iprint.gt.1) write(6,*) 'ic,Rcx(ic),Rcy(ic)'
        ic=0
        do j=-Nc,Nc
        do i=-Nc,Nc
          vx=real(i,kind)-origx
          vy=real(j,kind)-origy
          r1=(vx*a1x+vy*a1y)/real(a1x**2+a1y**2,kind)
          r2=(vx*a2x+vy*a2y)/real(a2x**2+a2y**2,kind)
          if(r1.lt.1.0_kind.and.r1.gt.0.0_kind.and.
     &      r2.lt.1.0_kind.and.r2.gt.0.0_kind) then
            ic=ic+1
            Rcx(ic)=i           ! Rcx, like i.e. a1x is integer
            Rcy(ic)=j
            if(iprint.gt.1) write(6,*) ic,Rcx(ic),Rcy(ic)
          end if
        end do
        end do
        if(ic.ne.Nc.or.Rcx(1).ne.0.or.Rcy(1).ne.0) then
          write(6,*) 'Bug in Rc(ic) loop',ic,Rcx(1),Rcy(1)
          stop
        end if
        
        if(iprint.ge.1) then
          write(6,*) '  '             
          write(6,*) 'Equivalency table r'                       
          write(6,*) '   ic     iequiv   icrequ(ic,ieqiv)'  
        end if     

        do ic=1,Nc
c  x --> isx*[(1-iex)*x + iex*y]
c  y --> isy*[(1-iex)*y + iex*x]
        do iequiv=1,ngroup
          n1=isx(iequiv)*((1-iex(iequiv))*Rcx(ic)+iex(iequiv)*Rcy(ic))  ! Generate the 
          m1=isy(iequiv)*((1-iex(iequiv))*Rcy(ic)+iex(iequiv)*Rcx(ic))  ! equivalent points
          do n=1,2  ! map them back into the cluster.
            vx=real(n1,kind)-origx
            vy=real(m1,kind)-origy
            r1=(vx*a1x+vy*a1y)/
     &         real(a1x*a1x+a1y*a1y,kind)
            r2=(vx*a2x+vy*a2y)/
     &         real(a2x*a2x+a2y*a2y,kind)
            if(r1.lt.0.0) then
              n1=n1+a1x
              m1=m1+a1y
            else if(r1.gt.1.0) then
              n1=n1-a1x
              m1=m1-a1y
            end if
            if(r2.lt.0.0) then
              n1=n1+a2x
              m1=m1+a2y
            else if(r2.gt.1.0) then
              n1=n1-a2x
              m1=m1-a2y
            end if
          end do
c         Figure out which point this is
          itest=0
          do jc=1,Nc
            if(n1.eq.Rcx(jc).and.m1.eq.Rcy(jc)) then
              icrequ(ic,iequiv)=jc
              itest=1
            end if
          end do
          if(itest.ne.1) then
            write(6,*) 'No equivalent found in icrequ',ic,iequiv
            stop
          end if 
          if(iprint.ge.1) 
     &        write(6,*) ic,iequiv,icrequ(ic,iequiv)                
        end do
        end do                               


c       Now create the tables for the differences of R
        if(iprint.ge.1) then
          write(6,*) '  '                
          write(6,*) 'differences in r'            
          write(6,*) '  ic   jc  icrdiff(ic,jc)'
        end if   
        do i=1,Nc
        do j=1,Nc
          Rx=real(Rcx(i)-Rcx(j),kind)
          Ry=real(Rcy(i)-Rcy(j),kind)
          itest=0
          do n=-1,1
          do m=-1,1               
            vx=Rx-origx + n*a1x + m*a2x
            vy=Ry-origy + n*a1y + m*a2y
            r1=(vx*a1x+vy*a1y)/real(a1x**2+a1y**2,kind)
            r2=(vx*a2x+vy*a2y)/real(a2x**2+a2y**2,kind)
            if(r1.lt.1.0_kind.and.r1.gt.0.0_kind.and.
     &        r2.lt.1.0_kind.and.r2.gt.0.0_kind) then
c             (Rx,Ry) is located within the cluster.  
c             Now we must find which point it is.
              do ic=1,Nc
                r1 = Rx + n*a1x + m*a2x
                r2 = Ry + n*a1y + m*a2y
                if(abs(r1-Rcx(ic)).lt.eps.and.
     &             abs(r2-Rcy(ic)).lt.eps) then
                   icrdiff(i,j)=ic
                   itest=1
                end if
              end do
              if(iprint.ge.1) write(6,*) i,j,icrdiff(i,j) 
            end if
          end do
          end do
          if(itest.ne.1) then
            write(6,*) 'no icrdiff(i,j) found',i,j
            stop
          end if
        end do
        end do
        
        
        
c                          K-space when Nc=10 
c                         -------------------
c         _        ^   ^     _        ^    ^
c         g1 = pi(3x + y)/5, g2 = pi(-x + 3y)/5
c
c                                ^
c                                | Ky     
cc---c---Q---c---c---c---c---c---c---c---c---c---Q---c---c---c---c   Legend 
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |  --------
cc---c---c---c---c---c---c---c---c---Q---c---c---c---c---c---c---c  K: in BZ
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |  Q: not in
cc---c---c---c---c---c---Q---c---c---c---c---c---c---c---c---c---Q
c|   |   |   |   |   |   |   |   |   |   |   |   |   | 4 |   |   |      K
cc---c---c---Q***************************************K---c---c---c     i,j
c|   |   |   *   |   |   |   |   |   |   | 3 |   |   *2,1|   |   |
cQ---c---c---*---c---c---c---c---c---c---K---c---c---*---c---c---c
c|   |   |   *   |   |   |   | 9 |   |   |1,1|   |   *   |   |   |
cc---c---c---*---c---c---c---K---c---c---c---c---c---*---c---c---c
c|   |   |   *   | 7 |   |   |0,1|   |   |   |   |   *   |   |   |
cc---c---c---*---K---c---c---c---c---c---c---c---c---*---Q---c---c
c|   |   |   * -1,1  |   |   |   |   |   |   | 2 |   *   |   |   |
cc---Q---c---*---c---c---c---c---c---c---c---K---c---*---c---c---c
c|   |   |   *   |   |   |   |   | 1 |   |   |1,0|   *   |   |   |
cc---c---c---*---c---c---c---c---K---c---c---c---c---*---c---c---c--> Kx
c|   |   |   *   |   | 6 |   |   |0,0|   |   |   |   *   |   |   |
cc---c---c---*---c---K---c---c---c---c---c---c---c---*---c---Q---c
c|   |   |   *   | -1,0  |   |   |   |   |   |   |10 *   |   |   |
cc---c---Q---*---c---c---c---c---c---c---c---c---K---*---c---c---c
c|   |   |   *   |   |   |   |   |   | 8 |   |  1,-1 *   |   |   |
cc---c---c---*---c---c---c---c---c---K---c---c---c---*---c---c---c
c|   |   |   *   |   |   | 5 |   |  0,-1 |   |   |   *   |   |   |
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
c       Create the cluster momentum tables Kx(ic) and Ky(ic)

c       First find the points in the irreducible wedge
        if(iprint.ge.1) write(6,*) '  '                            
        if(iprint.ge.1) write(6,*) ' ic  i  j   Kcx(ic)   Kcy(ic)  '
        ic=0
        ic00=1
        icpipi=1
        do i=0,Nc
        do j=0,i
          Kx=i*g1x + j*g2x
          Ky=i*g1y + j*g2y
          if(Kx.gt.-pi+eps.and.Kx.lt.pi+eps.and.
     &      Ky.gt.-pi+eps.and.Ky.lt.pi+eps) then
c           Kx and Ky fall within the first zone
            ic=ic+1
            Kcx(ic)=Kx
            Kcy(ic)=Ky
            ig(ic)=i
            jg(ic)=j
            if(abs(Kx).lt.eps.and.abs(Ky).lt.eps) ic00=ic
            if(abs(Kx-pi).lt.eps.and.abs(Ky-pi).lt.eps) icpipi=ic
            if(iprint.ge.1) write(6,5) ic,i,j,Kcx(ic),Kcy(ic)                    !0
          end if
        end do
        end do
 5      format(1x,i2,1x,i2,1x,i2,2x,f9.6,2x,f9.6)
        
        Ncw=ic
        if(iprint.ge.1) write(6,*) 'Nc= ',Nc,' Ncw=',Ncw

c       Now find the remaining points within the 1Bz and outside the wedge.     
        do i=-Nc,Nc
        do j=-Nc,Nc
          Kx=i*g1x + j*g2x
          Ky=i*g1y + j*g2y
          if(i.ge.0.and.j.ge.0.and.i.ge.j) then
c           This K is within the wedge angle
            continue
          else if(Kx.gt.-pi+eps.and.Kx.lt.pi+eps.and.
     &      Ky.gt.-pi+eps.and.Ky.lt.pi+eps) then
c           Kx and Ky fall within the first zone
            ic=ic+1
            Kcx(ic)=Kx
            Kcy(ic)=Ky
            ig(ic)=i
            jg(ic)=j
            if(abs(Kx).lt.eps.and.abs(Ky).lt.eps) ic00=ic
            if(abs(Kx-pi).lt.eps.and.abs(Ky-pi).lt.eps) icpipi=ic
            if(iprint.ge.1) write(6,5) ic,i,j,Kcx(ic),Kcy(ic)    
          end if
        end do
        end do
        
        if(ic.ne.Nc) then
          write(6,*) 'Bug in Kc(ic) loop'
          stop
        end if
        if(iprint.ge.1) then
          write(6,*) 'ic00=',ic00       
          write(6,*) 'icpipi=',icpipi
        end if   


c       Test for closure (?)
        do ic1=2,Nc     ! sum on K
          c1=(0.0,0.0)
          do ic2=1,Nc   ! sum on X
            c1=c1+exp(ii*(Kcx(ic1)*Rcx(ic2)+Kcy(ic1)*Rcy(ic2)))
          end do
          c1=c1/Nc
          if(abs(c1).gt.eps) then
            write(6,*) 'trouble with closure'
            stop
          end if
        end do
        
c       Create the ind. add. table ict(i,j) where N_sp=sqrt(float(Nc))

        do i=-N_sp,N_sp
        do j=-N_sp,N_sp
c         i=ig(ic)
c         j=jg(ic)
          Kx=i*g1x + j*g2x
          Ky=i*g1y + j*g2y
          if(Kx.lt.-pi+eps) Kx=Kx+2.0d0*pi
          if(Kx.gt.pi+eps)  Kx=Kx-2.0d0*pi
          if(Ky.lt.-pi+eps) Ky=Ky+2.0d0*pi
          if(Ky.gt.pi+eps)  Ky=Ky-2.0d0*pi
c         Kx=modulo(Kx+pi,2*pi)-pi
c         Ky=modulo(Ky+pi,2*pi)-pi
          if(abs(Ky).gt.pi.or.abs(Ky).gt.pi) then
            write(6,*) 'ict out of range',i,j
            stop
          end if
          itest=0
          do jc=1,nc
            if(abs(Kx-Kcx(jc)).lt.eps.and.abs(Ky-Kcy(jc)).lt.eps) then    
              ict(i,j)=jc
              itest=1
            end if
          end do
          if(itest.ne.1) then
            write(6,*) 'no Kc found ict'
            stop
          end if
        end do
        end do



c       Create the table ickdiff(ic1,ic2)
        if(iprint.ge.1) write(6,*) '   '            
        if(iprint.ge.1) write(6,*) ' ic1    ic2  ickdiff(ic1,ic2)  '     
        do ic1=1,Nc
        do ic2=1,Nc
          Kx=Kcx(ic1)-Kcx(ic2)
          Ky=Kcy(ic1)-Kcy(ic2)
          if(Kx.lt.-pi+eps) Kx=Kx+2.0_kind*pi
          if(Kx.gt.pi+eps) Kx=Kx-2.0_kind*pi
          if(Ky.lt.-pi+eps) Ky=Ky+2.0_kind*pi
          if(Ky.gt.pi+eps) Ky=Ky-2.0_kind*pi
c         Kdiffx(ic1,ic2)=Kx
c         Kdiffy(ic1,ic2)=Ky
c         Now we need to figure out where this point is!
          itest=0
          do ic=1,nc
            if(abs(Kx-Kcx(ic)).lt.eps.and.abs(Ky-Kcy(ic)).lt.eps) then    
              ickdiff(ic1,ic2)=ic
              itest=1
            end if
          end do
          if(itest.ne.1) then
            write(6,*) 'No ickdiff(ic1,ic2) found',ic1,ic2
            stop
          end if
          if(iprint.ge.1) write(6,*) ic1,ic2,ickdiff(ic1,ic2)                       !0
        end do
        end do

c       Create the table ickplus(ic1,ic2)
        if(iprint.ge.1) write(6,*) '   '           
        if(iprint.ge.1) write(6,*) ' ic1    ic2  ickplus(ic1,ic2)  '
        do ic1=1,Nc
        do ic2=1,Nc
          Kx=Kcx(ic1)+Kcx(ic2)
          Ky=Kcy(ic1)+Kcy(ic2)
          if(Kx.lt.-pi+eps) Kx=Kx+2.0_kind*pi
          if(Kx.gt.pi+eps) Kx=Kx-2.0_kind*pi
          if(Ky.lt.-pi+eps) Ky=Ky+2.0_kind*pi
          if(Ky.gt.pi+eps) Ky=Ky-2.0_kind*pi
c         Kplusx(ic1,ic2)=Kx
c         Kplusy(ic1,ic2)=Ky
c         Now we need to figure out where this point is!
          itest=0
          do ic=1,nc
            if(abs(Kx-Kcx(ic)).lt.eps.and.abs(Ky-Kcy(ic)).lt.eps) then    
              ickplus(ic1,ic2)=ic
              itest=1
            end if
          end do
          if(itest.ne.1) then
            write(6,*) 'No ickplus(ic1,ic2) found',ic1,ic2
            stop
          end if
          if(iprint.ge.1) write(6,*) ic1,ic2,ickplus(ic1,ic2)  
        end do
        end do
        
        
c       Now find the points equivalent to (i,j), which should be given by
c       
c       i ->- +/- i            Care must be taken to assure that
c         \ /                  these points fall within the zone.
c          x                   If not, then the point i,j is 
c         / \                  mapped back to itself.
c       j ->- +/- j
        if(iprint.ge.1) write(6,*) '   '               
        if(iprint.ge.1) write(6,*) 'equivalency in k'  
        if(iprint.ge.1) write(6,*) '   ic     iequiv   ickequ(ic,ieqiv)'

        do ic=1,Nc
          i=ig(ic)
          j=jg(ic)
c         Kx and Ky fall within the zone, by construction, now look 
c         for equivalent points
c
c  n --> isx*[(1-iex)*i + iex*j]
c  m --> isy*[(1-iex)*j + iex*i]
c
          do iequiv=1,ngroup
            n=isx(iequiv)*((1-iex(iequiv))*i + iex(iequiv)*j)
            m=isy(iequiv)*((1-iex(iequiv))*j + iex(iequiv)*i)
            Kx=n*g1x + m*g2x
            Ky=n*g1y + m*g2y
c           These are the equivalent points.  Now map to 1st BZ
            if(Kx.lt.-pi+eps) Kx=Kx+2.0*pi
            if(Kx.gt.pi+eps) Kx=Kx-2.0*pi
            if(Ky.lt.-pi+eps) Ky=Ky+2.0*pi
            if(Ky.gt.pi+eps) Ky=Ky-2.0*pi
c           Now we need to figure out where this point is!
            itest=0
            do jc=1,nc
              if(abs(Kx-Kcx(jc)).lt.eps
     &           .and.
     &           abs(Ky-Kcy(jc)).lt.eps) then
                ickequ(ic,iequiv)=jc
                itest=1
              end if
            end do
            if(iprint.ge.1) write(6,*) ic,iequiv,ickequ(ic,iequiv)               
            if(itest.eq.0) then
              write(6,*) 'no equivalent K-point found'
              write(6,*) 'j1,n2,m2',j1,n2,m2          
              write(6,*) 'Kx,Ky',Kx,Ky                
              stop
            end if
          end do
        end do


       
c       Now generate a table which maps any K to an equivalent point
c       in the irreducible wedge (IW), and a table for the degeneracy 
c       of each point in the IW.  For example, for Nc=18
c
c                                ^
c                                | Ky     
cc---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c   Legend 
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |  --------
cc---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c  K: in BZ
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |  c: not in
cc---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c  w: boundary
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |     of IW
cc---c---Q*******c*******c*******K w w w c w w w c w w w K---c---c  
c|   |   * g |   |   |   |   |   4   |   4   |   4   |   1   |   |   
cc---c---*---g---g---c---c---c---w---c---c---c---w---c---*---c---c
c|   |   *   | g |   |   |   |   w   |   |   |   | w |   *   |   |   
cc---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c
c|   |   *   |   | g |   |   |   4   |   8   | w 4   |   *   |   |
cc---c---*---c---c---K---c---c---w---c---c---K---c---c---*---c---c
c|   |   *   |   |   | g |   |   w   |   | w |   |   |   *   |   |
cc---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c
c|   |   *   |   |   |   | g |   4   | w 4   |   |   |   *   |   |
cc---c---*---c---c---c---c---g---w---w---c---c---c---c---*---c---c
c|   |   *   |   |   |   |   | g w w |   |   |   |   |   *   |   |
cc---c---c---c---c---c---c---c---K---c---c---c---c---c---c---c---c--> Kx
c|   |   *   |   |   |   |   |   1   |   |   |   |   |   *   |   |
cc---c---*---c---c---c---c---c---c---c---c---c---c---c---*---c---c
c|   |   *   |   |   |   |   |   |   |   |   |   |   |   *   |   |
cc---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c
c|   |   *   |   |   |   |   |   |   |   |   |   |   |   *   |   |
cc---c---*---c---c---c---c---c---c---c---c---c---c---c---*---c---c
c|   |   *   |   |   |   |   |   |   |   |   |   |   |   *   |   |
cc---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c
c|   |   *   |   |   |   |   |   |   |   |   |   |   |   *   |   |
cc---c---*---c---c---c---c---c---c---c---c---c---c---c---*---c---c
c|   |   *   |   |   |   |   |   |   |   |   |   |   |   *   |   |
cc---c---c*******c*******c*******c*******c*******c*******c---c---c
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
cc---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
cc---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
cc---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c---c
c
c       For Nc=8
c
c                                ^
c                                | Ky     
cc---c---Q---c---c---c---c---c---c---c---c---c---Q---c---c---c---c   Legend 
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |  --------
cc---c---c---c---c---c---c---c---c---Q---c---c---c---c---c---c---c  K: in BZ
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |  Q: not in
cc---c---c---c---c---c---Q---c---c---c---c---c---c---c---c---c---Q  w: boundary
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |     of IW
cc---c---c---Q*******************K w w w w w w w w w K---c---c---c  
c|   |   |   * g |   |   |   |   w   |   |   |   | w 1   |   |   |   
cQ---c---c---*---g---c---c---c---w---c---c---c---w---*---c---c---c
c|   |   |   *   | g |   |   |   w   |   |   | w |   *   |   |   |
cc---c---c---*---c---g---c---c---w---c---c---w---c---*---c---c---c
c|   |   |   *   |   | K |   |   w   |   | K |   |   *   |   |   |
cc---c---c---*---c---c---g---c---w---c---w-4-c---c---*---Q---c---c
c|   |   |   *   |   |   | g |   w   | w |   |   |   *   |   |   |
cc---Q---c---*---c---c---c---g---w---w---c---c---c---*---c---c---c
c|   |   |   *   |   |   |   | g w w |   |   |   |   *   |   |   |
cc---c---c---*---c---c---c---c---K---c---c---c---c---K---c---c---c--> Kx
c|   |   |   *   |   |   |   |   1   |   |   |   |   2   |   |   |
cc---c---c---*---c---c---c---c---c---c---c---c---c---*---c---Q---c
c|   |   |   *   |   |   |   |   |   |   |   |   |   *   |   |   |
cc---c---Q---*---c---c---c---c---c---c---c---c---c---*---c---c---c
c|   |   |   *   |   | K |   |   |   |   | K |   |   *   |   |   |
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
c
c       For Nc=8, only 4 of the 8 points are in the irreducible wedge.
c       Thus, all storage arrays may be reduced in size by a factor 
c       of 2.  As Nc-->oo, the saving factor approaches 8.
c       
c
c       To achieve this mapping to the irreducible wedge, note that 
c       each K can be written as
c             _       _
c       K = n g  +  m g
c              1       2
c
c       The mapping points in the irreducible wedge
c       
c       [n,m] --> [ max(|n|,|m|),min(|n|,|m|) ]
c
c       Create the table ickmap(ic) and ickdeg(ic)
        ickdeg=0
        do ic=1,Nc
          i=ig(ic)
          j=jg(ic)
          n=max(abs(i),abs(j))
          m=min(abs(i),abs(j))
          Kx=n*g1x + m*g2x
          Ky=n*g1y + m*g2y
c         Find the jc that corresponds to (Kx,Ky)     
          itest=0    
          do jc=1,Nc
            if(abs(Kx-Kcx(jc)).lt.eps.and.abs(Ky-Kcy(jc)).lt.eps) then
              ickmap(ic)=jc
              ickdeg(jc)=ickdeg(jc)+1
              itest=1
            end if
          end do
          if(itest.eq.0) then
            write(6,*) 'mapping failed for'
            write(6,*) 'ic',ic
            write(6,*) 'Kx,Ky',Kx,Ky
            stop
          end if
        end do
        
        
c       Form Epsbar(K)
        Epsbar=0.0
        do ic=1,Nc
          do i=-nover+1,nover
          do j=-nover+1,nover
            kx=Kcx(ic) + r1*((real(i,kind)-0.5)*g1x +
     &                       (real(j,kind)-0.5)*g2x)
            ky=Kcy(ic) + r1*((real(i,kind)-0.5)*g1y +
     &                       (real(j,kind)-0.5)*g2y)
            Epsbar(ic)=Epsbar(ic) -0.5*(cos(kx)+cos(ky)) -
     &                 tprime*(cos(kx)*cos(ky)-1.0)
          end do
          end do
          Epsbar(ic)=Epsbar(ic)/real((2*nover)**2,kind)
        end do

        return
        end     
