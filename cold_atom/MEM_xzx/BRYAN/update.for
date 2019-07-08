      subroutine update
c******************************************************************************
c     For a fixed alpha, solves the MaxEnt equation by the
c     quasi-Newton method with a Marquardt cut-off on the step size.
c******************************************************************************
c
c     If Q= alpha S - L, where S is the entropy of the image f and 
c     L is 1/2 times the extensive chi^2., then the condition of an
c     extremum, or a solution for fixed alpha, is
c
c     alpha Del S - Del L =0 -alpha ln(f_i/m_i)= SUM_j T_ji dL/dF_j
c  
c     where F_j is a component of Tf.  The matrix T is replaced by its
c     SVD, U Sigma V[T], where V and U are unitary, and we choose
c
c     f = m exp(Uu)
c
c     so this equation becomes
c
c     -alpha u = Sigma V[T] dL/dF == g
c
c     which defines g
c
c       Bryan uses a Newton search to find this solution, starting from
c     an arbitrary u.  The increment at each step is given by
c
c     J del u = -alpha u - g
c
c     where the Jacobian J
c
c     J= alpha I + dg/du
c
c     dg/du = Sigma V[T] W V Sigma U[T] {f} U == MK
c 
c     where W= {1/sigma_l^2} is the diagonal matrix composed of the 
c     errors and {f} is the diagonal matrix composed of the image.
c
c     M = Sigma V[T] W V Sigma and K= U[T] {f}  are symmetric
c     ns by ns matrices.  Then the Jacobian
c
c     J = alpha I + MK
c
c     and 
c
c     (alpha I + MK) del u = - alpha u - g
c
c     Unfortunately, this algorithm will be numerically unstable, since
c     some values of del u will be large enough that the second order 
c     approximation used in invalid.  In order to define the step size
c     we must define a metric for the space u.  In the image space f, it
c     is just -DelDel S= {1/f}.  This distance may be transformed to u
c    
c     del f[T] {1/f} del f = del u[T] K del u
c
c     so K, defined above, is the metric in u.  The corressponding 
c     step lenght in u is made by augmenting J with mu I
c
c     { (alpha + mu)I + MK} del u = -alpha u - g
c
c     where mu is adjusted so that 
c
c     del u[T] K del u <= SUM_i m_i
c
c       We can make this search even more efficient if the Newton equation 
c     above is diagonalized so that only Order(ns) operations are
c     required for each alpha mu pair.  First, the metric K is diagonalized
c
c     KP = P Theta
c
c     where Theta is the corressponding matrix of eigenvalues {xi} and P is
c     orthogonal.  To diagonalze M, then let
c
c     A == {sqrt(xi)} P[T] M P {sqrt(xi)}
c
c     and we solve the second eigenvalue problem A R =R Lambda, where
c     Lambda is the diagonal matrix of egenvalues of A and R is the 
c     corressponding unitary transformation.  Let 
c
c     Y = P {1/sqrt(xi)} R, so that
c
c     Y[-T] Y[-1] = K.  Then the Newton equation becomes
c
c     {(alpha+mu)I + Lambda} Y[-1] delu = -alpha Y[-1]u - Y[-1] g
c
c     which are ns independent equations for Y[-1] delu and the step
c     length is 
c
c     delu[T] K delu = | Y[-1] delu|^2
c
c*******************************************************************************
      use Global
      implicit none
      integer l,j,i,ier,iupdate, muflag, eql, icontrol
      double precision r1,norm1, norms, norml,
     &                 dnrm2,ddot,ax,bx,tol,cx,xmax,golden,probmax,
     &                 test1,alphanew, mumin, mumax, ls, los

c     declare some automatic arrays on the stack
      double precision G(nf), RHS(nf),YIRHS(nf), ds(nf), dl(nf),
     &                 fv1(nf),fv2(nf)
c     declare some allocatable arrays, to go on the heap (too big for the stack)
      double precision, allocatable :: K(:,:), A(:,:), rmat1(:,:), rmat2(:,:), 
     &                                 YI(:,:), YIT(:,:), MYIT(:,:), work(:)

      external funct

c*******************************************************************************
c     Allocate some arrays that should go on the heap, not the stack.
      allocate (K(nf,nf))
      allocate (A(nf,nf))            
      allocate (rmat1(nf,nf))            
      allocate (rmat2(nf,nf))            
      allocate (YI(nf,nf))            
      allocate (YIT(nf,nf))            
      allocate (MYIT(nf,nf))            
c*******************************************************************************

c     Iterate until a solution for fixed alpha is found.
      do iupdate=1,iupdate_max
     
c       Form the reconstructed data Tf(ndata)=T(ndata,nf)*f(nf)
        call dgemv('N',ndata,nf,1.0d0,T,ndata,f,1,0.0d0,Tf,1)

c       Form the gradient of L
        do l=1,ndata
           dl(l) = (Tf(l)-data(l))/error(l)**2
        end do

c       Form G = sigma*v[T]*dl
c       T*dl = G
        call dgemv('T',ndata,ns,1.0d0,v,ndata,dl,1,0.0d0,G,1)
        G(1:ns) = s(1:ns)*G(1:ns)
        
c       Form RHS = -alpha*u - G
        RHS(1:ns) = - alpha*u(1:ns) - G(1:ns)
        
c       Construct the Jacobian in its diagonal space
c       ref: Golub & vanLoan, Matrix Computations, pg.469; Bryan loc.cit.
c       First form K =U[T]*f*U
        do i=1, nf
           do j = 1, ns
              rmat1(i,j) = f(i)*Umat(i,j)
           end do
        end do
        
c       K(ns,ns) = Umat(ns,nf)*rmat1(nf,ns)
        call dgemm('T','N',ns,ns,nf,1.0d0,Umat,nf,rmat1,nf,
     &             0.0d0,K,nf)
        rmat1=K
      
c       SVD of rmat1 = K = P diag[Z] rmat2
        call dgesvd('A','N',ns,ns,rmat1,nf,Z,P,nf,rmat2,nf,
     &              r1,-1,ier)
        allocate (work(int(r1)+1))
        call dgesvd('A','N',ns,ns,rmat1,nf,Z,P,nf,rmat2,nf,
     &              work,int(r1)+1,ier)
        deallocate (work)
        
c       Form A = sqrt(Z)*P[T]*M*P*sqrt(Z)
        do i = 1, ns
           do j = 1, ns
              if(Z(j).lt.0.0) write(6,*) 'negative z(j), j=',j
              rmat2(i,j) = P(i,j)*sqrt(Z(j))
           end do
        end do
c       M*rmat2 = rmat1     
        call dgemm('N','N',ns,ns,ns,1.0d0,XM,nf,rmat2,nf,0.0d0,
     &             rmat1,nf)
c          
c       rmat2[T] *rmat1 = A
        call dgemm('T','N',ns,ns,ns,1.0d0,rmat2,nf,rmat1,nf,0.0d0,
     &             A,nf)
        rmat1=A
        
c       SVD of A = rmat1 = R diag[LAMBDA] rmat2
        call dgesvd('A','N',ns,ns,rmat1,nf,LAMBDA,R,nf,rmat2,nf,
     &              r1,-1,ier)
        allocate (work(int(r1)+1))
        call dgesvd('A','N',ns,ns,rmat1,nf,LAMBDA,R,nf,rmat2,nf,
     &              work,int(r1)+1,ier)
        deallocate (work)
      

c       Form Y[-1] = R[T]*sqrt(Z)*P[T] and Y[-T] = P*sqrt(Z)*R
        do i = 1, ns
           do j = 1, ns
              rmat1(i,j) = sqrt(Z(i))*P(j,i)
              rmat2(i,j) = sqrt(Z(i))*R(i,j)
           end do
        end do
c       R[T] * rmat1 = YI      
        call dgemm('T','N',ns,ns,ns,1.0d0,R,nf,rmat1,nf,0.0d0,
     &              YI,nf)
c       P * rmat2 = YIT = Y[-T]
        call dgemm('N','N',ns,ns,ns,1.0d0,P,nf,rmat2,nf,0.0d0,
     &             YIT,nf)

c       Form MYIT = M*Y[-T]
        call dgemm('N','N',ns,ns,ns,1.0d0,XM,nf,YIT,nf,0.0d0,
     &             MYIT,nf)
        
c       Form YIRHS = -alpha*Y[-1]*u - Y[-1]*g
        call dgemv('N',ns,ns,1.0d0,YI,nf,RHS,1,0.0d0,YIRHS,1)

     
        
c       Adjust step-size so it is consistent with the linerization.
c       I.e. adjust mu and calculate
c       Y[-1}*du = Y[-1]*(-alpha*u-G)/(alpha*I+LAMBDA+mu)
c       call control
        mu = 0.0d0; mumin=0.0d0; mumax=0.0d0; muflag=0
        
c       Generate los, the maximum step size
        los=step*sum(f(1:nf))
       
        do icontrol=1,icontrol_max
c         Calculate du(i) and ls
          du(1:ns)=YIRHS(1:ns)/(mu+lambda(1:ns)+alpha)
          ls=sum(du(1:ns)**2)

c         Now test if MU-chop is finished, if so, exit with a successful du
          if (ls.le.los.and.muflag.eq.0.or.eql(ls,los,precision).eq.1) exit
          
          if (ls.gt.los) then
c            mu must be increased (mu=0.0, first iterate)
             if (muflag.eq.0) then ! start with mu=10
        	mu=10.0d0; muflag=1
             else
        	call uchop (mu,mumin,mumax)
             endif
          else    ! mu must be decreased
             call dchop (mu,mumin,mumax)
          endif
        end do
        
        if(icontrol.eq.icontrol_max) then
          write(6,*) '*************************************'
          write(6,*) 'WARNING, Control Loop Maximum Reached'
          write(6,*) '*************************************'
        end if
    
        
        
c       Compute new image
c       u = u + du
c       du = (-alpha*u - g - M*Y[-T]*(Y[-1]*du))/(alpha+mu)
c       MYIT du = fv1
        call dgemv('N',ns,ns,1.0d0,MYIT,nf,du,1,0.0d0,fv1,1)
        u(1:ns) = u(1:ns) + (RHS(1:ns) - fv1(1:ns))/(alpha+mu)
        
c       U u = fv1
        call dgemv('N',nf,ns,1.0d0,Umat,nf,u,1,0.0d0,fv1,1)
        do i = 1, nf
           if(fv1(i).lt.-60.0) then
c             write(6,*) 'fv1 woops'
              iuflow=1
              fv1(i)=-60.0d0
           end if
           f(i) = model(i)*exp(fv1(i))
        end do
c       Form Lo.  First form the reconstructed image Tf=T*f
        call dgemv('N',ndata,nf,1.0d0,T,ndata,f,1,0.0d0,Tf,1)
        lo=0.5d0*sum(((Tf(1:ndata)-data(1:ndata))/error(1:ndata))**2)
        aim=2.0d0*lo/float(ndata)
        
        
c       Form So
        so=sum(-f(1:nf)*log(f(1:nf)/model(1:nf))+f(1:nf)-model(1:nf))

c       Measure tests (remember the metric! K=Y[-T]*Y[-1])
        call dgemv('N',ns,ns,1.0d0,YI,nf,u,1,0.0d0,ds,1)
        call dgemv('N',ns,ns,1.0d0,YI,nf,G,1,0.0d0,dl,1)
        norms = dnrm2(ns,ds,1)
        norml = dnrm2(ns,dl,1)
        
        call dgemv('N',ns,ns,1.0d0,YI,nf,RHS,1,0.0d0,fv1,1)
        test1= 2.*ddot(ns,fv1,1,fv1,1)/(alpha*norms+norml)**2
        if (mu.lt.0.001.and.test1.lt.testf) exit
      end do
      if(iupdate.eq.iupdate_max) then
        write(6,*) '*************************************'
        write(6,*) 'WARNING, Update Loop Maximum Reached'
        write(6,*) '*************************************'
      end if

c     Convergence has been obtained, now we must iterate alpha or quit.
      if(cflag.eq.0)then
c       Perform a golden mean search for the maximum of 
c       P(lpha|D,m)~= Prod Sqrt(alpha/(alpha+lambda_i)) exp(Q) P(alpha)
c       in the interval (ax,cx) starting with the value bx.  
c       See NR Sec. 10.1 for details.
c        ax=1.0d-7
c        cx=1.0d6
        bx=alpha
        ax=0.1d0*bx
        cx=10.d0*bx
        tol = 1.0d-05
        probmax = golden(ax,bx,cx,funct,tol,xmax)
        alphanew = xmax
c       The next value of alpha is given by
        alpha=(1.0d0-alphadamp)*alphaold +alphadamp*alphanew
      else
        alphanew=alpha
      end if

c     Now determine the new value of ngood.
      norm1=sum((LAMBDA(1:ns))/(LAMBDA(1:ns)+alpha))
      ngood=norm1

c     For auto error scaling
      sig=sqrt(2.0*lo/(ndata-norm1))

c     Now compute to posterior probability for the present value of alpha
c     and accumulate average values
      post = 1.0d0
      do i = 1, ns
         post = post*alpha/(alpha+LAMBDA(i))
      end do
      if(aflag.eq.0.or.aflag.eq.1) then
c     Implement the Jeffrey prior.
        post = sqrt(post)*exp(alpha*so-lo)/alpha
      else
        post = sqrt(post)*exp(alpha*so-lo)
      end if

c      write(6,"(' P(alpha=',f10.2,'|D,m): ',e10.3)") alpha, post     
      weight = weight + post*abs(alpha-alphaold)
      if(aflag.eq.1.or.aflag.eq.3) then
         padm(iter) = post 
         alpt(iter) = alphaold
         alphabar = alphabar + post*abs(alpha-alphaold)*alpha
         sigbar = sigbar + post*abs(alpha-alphaold)*sig
         ngbar = ngbar + post*abs(alpha-alphaold)*norm1
         fbar(1:nf) = fbar(1:nf) + post*abs(alpha-alphaold)*f(1:nf)
      else
         padm(iter) = post
      endif
c
c     DeAllocate the heap arrays.
      deallocate (K)            
      deallocate (A)            
      deallocate (rmat1)            
      deallocate (rmat2)            
      deallocate (YI)            
      deallocate (YIT)            
      deallocate (MYIT)            

      return
      end


      double precision function funct(x)
c     assumes the Jeffery prior
      use Global
      implicit none
      double precision x
      integer i
      if(x.le.0.0) then
        funct=-1.0d305
        write(6,*) 'woops, too small'
        return
      end if
      if(x.ge.1.0d305) then
        write(6,*) 'woops,too large'
        stop
      end if
      funct = 1.
      do i = 1, ns
         funct = funct/(x+lambda(i))
      end do
c     the mimus sign makes the search look for a maximum
c         Assumes the Jeffery prior.  Its effect is divided out so 
c         alpha can equal zero
      if(aflag.eq.0.or.aflag.eq.1) then
        funct = -sqrt(funct*(x**(ns-2)))*exp(x*so-lo)
      else
        funct = -sqrt(funct*(x**(ns)))*exp(x*so-lo)
      end if

c
      return
      end
      

      double precision FUNCTION GOLDEN(AX,BX,CX,F,TOL,XMIN)
      implicit double precision (a-h,o-z)
      PARAMETER (R=.61803399,C=.38196602)
      X0=AX
      X3=CX
      IF(ABS(CX-BX).GT.ABS(BX-AX))THEN
        X1=BX
        X2=BX+C*(CX-BX)
      ELSE
        X2=BX
        X1=BX-C*(BX-AX)
      ENDIF
      F1=F(X1)
      F2=F(X2)
1     IF(ABS(X3-X0).GT.TOL*(ABS(X1)+ABS(X2)))THEN
        IF(F2.LT.F1)THEN
          X0=X1
          X1=X2
          X2=R*X1+C*X3
          F0=F1
          F1=F2
          F2=F(X2)
        ELSE
          X3=X2
          X2=X1
          X1=R*X2+C*X0
          F3=F2
          F2=F1
          F1=F(X1)
        ENDIF
      GOTO 1
      ENDIF
      IF(F1.LT.F2)THEN
        GOLDEN=F1
        XMIN=X1
      ELSE
        GOLDEN=F2
        XMIN=X2
      ENDIF
      RETURN
      END

      subroutine dchop (x,xmin,xmax)
c     This function chops x downward towards xmin

      double precision x, xmin, xmax
      xmax=x
      x=0.5*(xmax+xmin)

      return
      end
      
      
      subroutine uchop (x,xmin,xmax)
c     This function chops x upward towards xmax; however,
c     if xmax is zero, then it just sets x=2x.
      
      double precision x, xmin, xmax
      xmin=x
      if (xmax.lt.1.0e-14) then
         x=2.0*xmin
      else
         x=0.5*(xmin+xmax)
      endif
      
      return
      end
      
      
      integer function eql (x,y,prec)
c     This integer function determines if x and y are equal
c     to some precision set by prec.
      
      double precision x, y, prec
      eql=0
      if (abs((x-y)/(x+y)).lt.prec) eql=1
      
      return
      end


