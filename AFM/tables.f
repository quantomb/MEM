        subroutine tables                    
c********************************************************************
c       This subroutine makes the lookup tables used throughout the
c       code.
c********************************************************************
        use Global
        implicit none
c**************************************************************
        integer :: n,m,k
c**************************************************************
        N_sp=sqrt(float(Nc))
        


c       MATSUBARA FREQUeNCIES
        do n=-nwn,nwn
          wn(n)=(twor*n+oner)*pi*temp
        end do
        

 
        if(ndim.eq.2) then
          call tables_2D
        else
          if(myrank.eq.0) write(66,*) 'bad dimensionality'
        end if
        
        if(myrank.eq.0.and.iprint.ge.9) write(66,*) 'done tables'       !0
        return
        end
        
        
        subroutine tables_2D                  !TWO DIMENSIONAL VERSION
c********************************************************************
c       This subroutine makes the lookup tables used throughout the
c       TWO DIMENSIONAL code.
c********************************************************************
        use Global
        implicit none
c**************************************************************
        integer i,j,k,n,m,m1,n1,m2,n2,jmd,jm,ind,iequiv,i1,j1,
     &          i2,j2,ic,ic1,ic2,jc,ix,iy,jk,iflag,jkc,icR
        real(kind) :: etu2,Kx,Ky,r1,r2,Rx,Ry,vx,vy,origx,origy,
     &                Dkt
c**************************************************************

        
c       CLUSTER MOMENTUM and DUAL SPACE for the 2D square lattice

c        ndim=2
        if(ngroup.ne.8.and.myrank.eq.0) write(66,*) 'ngroup wrong'
 
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
c       Principle translation                      [ square     (1/2,1/2)
c       vectors for lattice tiling        Origins  |
c                                                  [ non-square  (1/2,0)
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



        If(Nc.eq.4) then                           
          a(1,1)=2                                       
          a(1,2)=0                                       
          Ncw=3
        else if(Nc.eq.8) then                           
          a(1,1)=2                                 
          a(1,2)=2    
           Ncw=4                  
       else if(Nc.eq.16) then                          
          a(1,1)=4                                 
          a(1,2)=0
          Ncw=6
        else
          if(myrank.eq.0) write(66,*) 'Bad Nc, bad, bad'                !0
          stop
        end if    
        a(2,1)=-a(1,2)
        a(2,2)= a(1,1)

c       Now calculate the principle translations in K-space
        g(1,1)=twor*pi*a(1,1)/real(a(1,1)*a(2,2)-a(2,1)*a(1,2),kind)
        g(1,2)=twor*pi*a(1,2)/real(a(1,1)*a(2,2)-a(2,1)*a(1,2),kind)
        g(2,1)=-g(1,2)
        g(2,2)= g(1,1)

c       Set the location of some convenient origin in real space
        if(a(1,1).eq.a(1,2)) then             ! The cell is a diamond
          origx= zeror
          origy=-halfr
        else if(a(1,1)*a(1,2).eq.0) then      ! The cell is a square
          origx=-halfr
          origy=-halfr
        end if

        if(iprint.ge.5.and.myrank.eq.0) 
     &    write(66,*) ' a1x,a1y,g1x,g1y',a(1,1),a(1,2),g(1,1),g(1,2)    !0
        if(iprint.ge.5.and.myrank.eq.0) 
     &    write(66,*) ' a2x,a2y,g2x,g2y',a(2,1),a(2,2),g(2,1),g(2,2)    !0

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

        if(iprint.ge.5.and.myrank.eq.0) write(66,*) '  '                !0
        if(iprint.ge.5.and.myrank.eq.0) 
     &       write(66,*) 'ic,Rcx(ic),Rcy(ic),isublat,invicRR,
     &                               icR, icicRR'                   !0
c
c       calculate a table of cluster point locations.
        ic=0
        icr=0
        invicRR=0
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
            if (mod(i+j,2).eq.0) then
               isublat(ic)=0
               icr=icr+1
               icicRR(icr)=ic
               invicRR(ic)=icR
            else 
               isublat(ic)=1
            endif
            if(iprint.ge.5.and.myrank.eq.0) 
     &         write(66,*) ic,Rc(1,ic),Rc(2,ic),isublat(ic),
     &         invicRR(ic),icR,  icicRR(icr)           !0
            end if
        end do
        end do
      
        if(icR.ne.NcR) then
          if(myrank.eq.0) write(66,*) 'Bug in ICRR(ic) loop'              !0
          stop
        end if

c       create a table mapping points outside the cluster back into it.

        if(iprint.ge.5.and.myrank.eq.0) write(66,*) '  '                !0
        if(iprint.ge.5.and.myrank.eq.0) 
     &       write(66,*) 'Equivalency table r'                          !0
        if(iprint.ge.5.and.myrank.eq.0) 
     &      write(66,*) '   ic     iequiv   icrequ(ic,iequiv)'          !0
        do ic=1,Nc  
        iequiv=0
        do j1=0,1                       
          i=j1*Rc(1,ic)+(1-j1)*Rc(2,ic)
          j=j1*Rc(2,ic)+(1-j1)*Rc(1,ic)
          do n2=-1,1,2
          do m2=-1,1,2
            n=n2*i
            m=m2*j
            iequiv=iequiv+1
            n1=n
            m1=m
            vx=real(n,kind)-origx
            vy=real(m,kind)-origy
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
            do jc=1,Nc
              if(n1.eq.Rc(1,jc).and.m1.eq.Rc(2,jc)) 
     &          icrequ(ic,iequiv)=jc
            end do
            if(iprint.ge.5.and.myrank.eq.0) 
     &          write(66,*) ic,iequiv,icrequ(ic,iequiv)                 !0
          end do
          end do
        end do
        end do
                      

c       Now create the tables for the differences of R
        if(iprint.ge.5.and.myrank.eq.0) write(66,*) '  '                !0
        if(iprint.ge.5.and.myrank.eq.0) 
     &       write(66,*) 'differences in r'                             !0
        if(iprint.ge.5.and.myrank.eq.0) 
     &       write(66,*) '  ic   jc  icrdiff(ic,jc)'                    !0
        do i=1,Nc
        do j=1,Nc
          Rx=real(Rc(1,i)-Rc(1,j),kind)
          Ry=real(Rc(2,i)-Rc(2,j),kind)
          do n=-1,1
          do m=-1,1               
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
     &             abs(r2-Rc(2,ic)).lt.eps) icrdiff(i,j)=ic
              end do
              if(iprint.ge.5.and.myrank.eq.0) 
     &             write(66,*) i,j,icrdiff(i,j)                         !0
            end if
          end do
          end do
        end do
        end do
        
c       Now create the tables for the neighbors of each site
        if(iprint.ge.5.and.myrank.eq.0) write(66,*) '  '                !0
        if(iprint.ge.5.and.myrank.eq.0) 
     &       write(66,*) 'neighbors of each site'                       !0
        if(iprint.ge.5.and.myrank.eq.0) 
     &       write(66,*) '  ic   j  neighbor(ic,j)'                     !0
        do i=1,Nc
        j=0
        do i1=-1,1,2    
        do j1=0,1
          ix=j1*i1
          iy=(1-j1)*i1  
          j= j+1  
c         
c         ix  iy  j            |3
c          0  -1  1            |   4
c         -1   0  2        ---------
c          0   1  3        2   |
c          1   0  4            |1
c
          Rx=real(Rc(1,i)+ix,kind)
          Ry=real(Rc(2,i)+iy,kind)
          do n=-1,1
          do m=-1,1               
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
     &             abs(r2-Rc(2,ic)).lt.eps) neighbor(i,j)=ic
              end do
              if(iprint.ge.5.and.myrank.eq.0) 
     &             write(66,*) i,j,neighbor(i,j)                        !0
            end if
          end do
          end do
          end do
        end do
        end do

        do i=1,4
           idwsgn(i)=(-1)**i
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
c|   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |      K
cc---c---c---Q***************************************K---c---c---c     i,j
c|   |   |   *   |   |   |   |   |   |   |   |   |   *2,1|   |   |
cQ---c---c---*---c---c---c---c---c---c---K---c---c---*---c---c---c
c|   |   |   *   |   |   |   |   |   |   |1,1|   |   *   |   |   |
cc---c---c---*---c---c---c---K---c---c---c---c---c---*---c---c---c
c|   |   |   *   |   |   |   |0,1|   |   |   |   |   *   |   |   |
cc---c---c---*---K---c---c---c---c---c---c---c---c---*---Q---c---c
c|   |   |   * -1,1  |   |   |   |   |   |   |   |   *   |   |   |
cc---Q---c---*---c---c---c---c---c---c---c---K---c---*---c---c---c
c|   |   |   *   |   |   |   |   |   |   |   |1,0|   *   |   |   |
cc---c---c---*---c---c---c---c---K---c---c---c---c---*---c---c---c--> Kx
c|   |   |   *   |   |   |   |   |0,0|   |   |   |   *   |   |   |
cc---c---c---*---c---K---c---c---c---c---c---c---c---*---c---Q---c
c|   |   |   *   | -1,0  |   |   |   |   |   |   |   *   |   |   |
cc---c---Q---*---c---c---c---c---c---c---c---c---K---*---c---c---c
c|   |   |   *   |   |   |   |   |   |   |   |  1,-1 *   |   |   |
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




c       Create the cluster momentum tables Kx(ic) and Ky(ic)
c       and reduce zone table KcR 

c       First find the points in the irreducible wedge.  They will be
c       indexed from 1 to Ncw
        if(iprint.ge.5.and.myrank.eq.0) 
     &       write(66,*) '  '                                           !0
        if(iprint.ge.5.and.myrank.eq.0) 
     & write(66,*) ' ic  i  j   Kcx(ic)   Kcy(ic) icR  icicRK invicRK ' !0
        icR=0
        ic=0
        ic0=1
        icpi=1
        invicRK=0
        do i=0,Nc
        do j=0,i
          Kx=i*g(1,1) + j*g(2,1)
          Ky=i*g(1,2) + j*g(2,2)
          if(Kx.gt.-pi+eps.and.Kx.lt.pi+eps.and.
     &      Ky.gt.-pi+eps.and.Ky.lt.pi+eps) then
c           Kx and Ky fall within the first zone
            ic=ic+1
            Kc(1,ic)=Kx
            Kc(2,ic)=Ky
            ig(ic)=i
            jg(ic)=j
            if(abs(Kx).lt.eps.and.abs(Ky).lt.eps) ic0=ic
            if(abs(Kx-pi).lt.eps.and.abs(Ky-pi).lt.eps) icpi=ic
          if(((g(1,1)*g(1,2).eq.0).and.((Kx.gt.0.and.
     &       abs(Kx)+abs(Ky).lt.pi+eps).or.(Kx.le.0.and.
     &       abs(Kx)+abs(Ky).lt.pi-eps))).or.(
     &       (g(1,1).eq.g(1,2)).and.((Ky.gt.0.and.
     &       abs(Kx)+abs(Ky).lt.pi+eps).or.(Ky.le.0.and.
     &       abs(Kx)+abs(Ky).lt.pi-eps)))) then
c            if((Kx.gt.0.and.abs(Kx)+abs(Ky).lt.pi+eps).or.
c     &       (Kx.le.0.and.abs(Kx)+abs(Ky).lt.pi-eps)) then
               icR=icR+1
               icicRK(icR)=ic
               invicRK(ic)=icR
            endif
            if(iprint.ge.5.and.myrank.eq.0) 
     & write(66,5) ic,i,j,Kc(1,ic),Kc(2,ic),icR,icicRK(icR),invicRK(ic)

          end if
        end do
        end do
 5      format(1x,i2,1x,i2,1x,i2,2x,f9.6,2x,f9.6,1x,i2,1x,i2,1x,i2)
        
        Ncw=ic
        NcwR=icR
        ntw=nl*Ncw
        if(iprint.ge.5.and.myrank.eq.0) 
     &       write(66,*) 'Nc= ',Nc,' Ncw=',Ncw,' NcwR=',NcwR               !0
        
c       Now find the remaining points within the 1Bz and outside the 
c       wedge.  These points are indexed from Ncw+1 to Nc       
        do i=-Nc,Nc
        do j=-Nc,Nc
          Kx=i*g(1,1) + j*g(2,1)
          Ky=i*g(1,2) + j*g(2,2)
          if(i.ge.0.and.j.ge.0.and.i.ge.j) then
c           This K is within the wedge angle
            continue
          else if(Kx.gt.-pi+eps.and.Kx.lt.pi+eps.and.
     &      Ky.gt.-pi+eps.and.Ky.lt.pi+eps) then
c           Kx and Ky fall within the first zone
            ic=ic+1
            Kc(1,ic)=Kx
            Kc(2,ic)=Ky
            ig(ic)=i
            jg(ic)=j
            if(abs(Kx).lt.eps.and.abs(Ky).lt.eps) ic0=ic
            if(abs(Kx-pi).lt.eps.and.abs(Ky-pi).lt.eps) icpi=ic
            if(((g(1,1)*g(1,2).eq.0).and.((Kx.gt.0.and.
     &       abs(Kx)+abs(Ky).lt.pi+eps).or.(Kx.le.0.and.
     &       abs(Kx)+abs(Ky).lt.pi-eps))).or.(
     &       (g(1,1).eq.g(1,2)).and.((Ky.gt.0.and.
     &       abs(Kx)+abs(Ky).lt.pi+eps).or.(Ky.le.0.and.
     &       abs(Kx)+abs(Ky).lt.pi-eps)))) then
                icR=icR+1
                icicRK(icR)=ic
                invicRK(ic)=icR
            endif

            if(iprint.ge.5.and.myrank.eq.0) 
     & write(66,5) ic,i,j,Kc(1,ic),Kc(2,ic),icR,icicRK(icR),invicRK(ic)

          end if
        end do
        end do
        
        if(ic.ne.Nc) then
          if(myrank.eq.0) write(66,*) 'Bug in Kc(ic) loop'              !0
          stop
        end if
        if(icR.ne.NcR) then
          if(myrank.eq.0) write(66,*) 'Bug in icR loop'              !0
          stop
        end if

        if(iprint.ge.5.and.myrank.eq.0) 
     &       write(66,*) 'ic0=',ic0                                     !0
        if(iprint.ge.5.and.myrank.eq.0) 
     &       write(66,*) 'icpi=',icpi                                   !0
        


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
        if(iprint.ge.5.and.myrank.eq.0) write(66,*) '   '               !0
        if(iprint.ge.5.and.myrank.eq.0) 
     &       write(66,*) ' ic1    ic2  ickdiff(ic1,ic2)  '              !0
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
     &       write(66,*) ic1,ic2,ickdiff(ic1,ic2)                       !0
        end do
        end do

c       Create the table ickplus(ic1,ic2), which takes the sum
c       between to K-points
        if(iprint.ge.5.and.myrank.eq.0) write(66,*) '   '               !0
        if(iprint.ge.5.and.myrank.eq.0) 
     &       write(66,*) ' ic1    ic2  ickplus(ic1,ic2)  '              !0
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
     &       write(66,*) ic1,ic2,ickplus(ic1,ic2)                       !0
        end do
        end do
        
        
c       Now find the points equivalent to (i,j), which should be given by
c       
c       i ->- +/- i            Care must be taken to assure that
c         \ /                  these points fall within the zone.
c          x                   If not, then the point i,j is 
c         / \                  mapped back to itself.
c       j ->- +/- j
c
        if(iprint.ge.5.and.myrank.eq.0) write(66,*) '   '               !0
        if(iprint.ge.5.and.myrank.eq.0) write(66,*) 'equivalency in k'  !0
        if(iprint.ge.5.and.myrank.eq.0) 
     &               write(66,*) '   ic     iequiv   ickequ(ic,ieqiv)'  !0
        do ic=1,Nc
          i=ig(ic)
          j=jg(ic)
c         Kx and Ky fall within the zone, by construction, now look 
c         for equivalent points
          iequiv=0
          do j1=0,1,1
          do n2=-1,1,2
          do m2=-1,1,2
            n=n2*i
            m=m2*j
            iequiv=iequiv+1
            Kx=j1*(n*g(1,1) + m*g(2,1)) + (1-j1)*(m*g(1,1) + n*g(2,1))
            Ky=j1*(n*g(1,2) + m*g(2,2)) + (1-j1)*(m*g(1,2) + n*g(2,2))
c           Kx=real(j1,kind)*(real(n,kind)*g(1,1) + real(m,kind)*g(2,1)) + 
c     &          real(1-j1,kind)*(real(m,kind)*g(1,1) + real(n,kind)*g(2,1))
c           Ky=real(j1,kind)*(real(n,kind)*g(1,2) + real(m,kind)*g(2,2)) + 
c     &          real(1-j1,kind)*(real(m,kind)*g(1,2) + real(n,kind)*g(2,2))
c           These are the equivalent points.  Now map to 1st BZ
            if(Kx.lt.-pi+eps) Kx=Kx+twor*pi
            if(Kx.gt.pi+eps) Kx=Kx-twor*pi
            if(Ky.lt.-pi+eps) Ky=Ky+twor*pi
            if(Ky.gt.pi+eps) Ky=Ky-twor*pi
c           Now we need to figure out where this point is!
            iflag=0
            do jc=1,nc
              if(abs(Kx-Kc(1,jc)).lt.eps
     &           .and.
     &           abs(Ky-Kc(2,jc)).lt.eps) then
                ickequ(ic,iequiv)=jc
                iflag=1
              end if
            end do
            if(iprint.ge.5.and.myrank.eq.0) 
     &           write(66,*) ic,iequiv,ickequ(ic,iequiv)                !0
            if(iflag.eq.0) then
              if(myrank.eq.0) write(66,*) 'no equivalent K-point found' !0
              if(myrank.eq.0) write(66,*) 'j1,n2,m2',j1,n2,m2           !0
              if(myrank.eq.0) write(66,*) 'Kx,Ky',Kx,Ky                 !0
              stop
            end if
          end do
          end do
          end do
          if(iequiv.ne.8) then
             if(myrank.eq.0) write(66,*) 'bug in equiv loop'            !0
             stop
          end if
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
c       Create the mapping table ickmap(ic) and the degneracy table ickdeg(ic)
        ickdeg=0
        do ic=1,Nc
          i=ig(ic)
          j=jg(ic)
          n=max(abs(i),abs(j))
          m=min(abs(i),abs(j))
          Kx=n*g(1,1) + m*g(2,1)
          Ky=n*g(1,2) + m*g(2,2)
c         Find the jc that corresponds to (Kx,Ky)     
          iflag=0    
          do jc=1,Nc
            if(abs(Kx-Kc(1,jc)).lt.eps
     &        .and.
     &         abs(Ky-Kc(2,jc)).lt.eps) then
              ickmap(ic)=jc
              ickdeg(jc)=ickdeg(jc)+1
              iflag=1
            end if
          end do
          if(iflag.eq.0) then
            if(myrank.eq.0) write(66,*) 'mapping failed for'            !0
            if(myrank.eq.0) write(66,*) 'ic',ic                         !0
            if(myrank.eq.0) write(66,*) 'Kx,Ky',Kx,Ky                   !0
            stop
          end if
        end do

c       Form Epsbar(K)
        Dkt=halfr/real(nover,kind)
        Epsbar1=zeror
        Epsbar2=zeror
        do jc=1,NcR
           ic=icicRK(jc)
          do i=-nover+1,nover
          do j=-nover+1,nover
            kx=Kc(1,ic) + Dkt*((real(i,kind)-halfr)*g(1,1) +
     &                       (real(j,kind)-halfr)*g(2,1))
            ky=Kc(2,ic) + Dkt*((real(i,kind)-halfr)*g(1,2) +
     &                       (real(j,kind)-halfr)*g(2,2))
            Epsbar1(jc) =Epsbar1(jc) -halfr*(cos(kx)+cos(ky)) 
            Epsbar2(jc) =Epsbar2(jc)- tprime*(cos(kx)*cos(ky)-oner)
          end do
          end do
          Epsbar1(jc)=Epsbar1(jc)/real((2*nover)**2,kind)
          Epsbar2(jc)=Epsbar2(jc)/real((2*nover)**2,kind)
        end do




        return
        end


