      program main
c*******************************************************************************
c     This is the main block of the MaxEnt code.
c
c***************************************************************************
c     For a review of MEM applied to the analytic continuation problem, see
c
c     Bayesian Inference and the Analytic Continuation of Imaginary-Time 
c     Quantum Monte Carlo Data, M. Jarrell, and J.E. Gubernatis, 
c     Physics Reports Vol. 269 #3, pp133-195, (May, 1996).
c
c***************************************************************************
c     Code History                                                         *
c***************************************************************************
c     * Original Gull/Skilling code, ~1987-1990 (Mark Jarrell)
c     * modified to do both automatic error, alpha scaling, and alpha 
c       integration. Implemented by Mark Jarrell and J.E. Gubernatis, 1990
c       following J.Skilling and R.K. Bryan, Mon. Not. R. Ast. Soc., 1984, 
c       211, 111-124 and R.K. Bryan, Eur. Biophys, J, 1990, 18, 165-174.
c     * Intrinsic correlation function were added about 1/5/91 (deleted)
c     * Historic maxent was added in 1991 (MJ)
c     * F90 modernization, linking to BLAS and LAPACK, and greatly enhanced
c       documentation 2005 (MJ).
c     * Stochastic Statistical Inference, allocatable arrays 2005 (MJ).
c***************************************************************************
c     The variables used in this algorithm are explained in module_Global.f
c     The main function of the algorithm is controlled by strongs method 
c     methodssi 
c       
c***************************************************************************
c       linked calls
c       BLAS
c         dnrm2         Compute the Euclidean length (L2 norm) of a vector.
c         ddot          forms the dot product of two vectors.
c         dgemv         matix-vector product
c         dgemm         matrix-matrix product
c       LAPACK
c         dgesvd        Singular Value Decomposition
c***************************************************************************


      use Global
      implicit none
      double precision test
      integer i
      character*30 string


c     First read the data and then initialize derived quantities
      call readin
c           
      call init

     
      if(iprint.gt.0.and.aflag.ne.4) then
        write(6,*) '  '
        write(6,*) ' iteration #     Alpha     P(alpha|D,m,...)'
      else if(iprint.gt.0) then
        write(6,*) '  '
        write(6,*) ' iteration #     Alpha         Chi^2/Ndata'
      end if

      if(aflag.ne.4) then
c*****************************************************************************
c       Search to find the maximum of P(alpha|m,G,...) for all but the
c       historic method (which is treated separately below).  For the
c       Bryan method, we assume that the initial alpha >> alphamode.
c*****************************************************************************
c       find the value of alpha that maximizes the posterior probability
c       for the image
        do i=1,niter
          iter = i
          alphaold=alpha
          call update        
          if(iprint.gt.0) write(6,"(5x,i4,3x,e13.5,3x,e13.5)") i,alpha,padm(iter)
          test = (alpha-alphaold)/(alpha+alphaold)
          if (alpha.ne.alphaold .and. abs(test).le.0.0005) then
            if(iprint.gt.0) write(6,*) 'Alpha converged'
            cflag=1
            exit
          endif
        end do
        alphaold = alpha; alphamode = alpha
        if(cflag.ne.1) then
          write(6,*) '********** WARNING Alpha not converged **************'
          write(6,*) '**** you may want to choose a differnt value for ****'
          write(6,*) '**** the initial alpha or increase niter        *****'
          write(6,*) '*********** WARNING *********************************'
        end if
 


        alphamax = alpha
        if (aflag.eq.1.or.aflag.eq.3) then  ! Use the Bryan Algorithm.
c         Test whether Bryan can be used here
          if(alphainit.lt.alphamode) then
            write(6,*) '*********** ERROR alphainit < alphamode *************'
            write(6,*) '**********Choose a larger initial alpha *************'
            write(6,*) '*********************** ERROR ***********************'
            stop
          else if(10*alphainit.lt.alphamode) then
            write(6,*) '******** WARNING 10*alphainit < alphamode ***********'
            write(6,*) '****** The initial alpha may be too small ***********'
            write(6,*) '****** Inspect posterior.xvgr to be sure  ***********'
            write(6,*) '********************* WARNING ***********************'
          end if
c*****************************************************************************
c         For the Bryan method, we need to marginalize the image over
c         the posterior of alpha, padm=P(alpha|m,G,...).   Alpha is assumed 
c         to be initialized to a large value.  As alpha converges in the
c         loop above, values of padm(alpha) are thus generated for large
c         values of alpha. This section of the code computes padm for 
c         20 small values of alpha.  The 19, 21 and 20.0 are arbitrary.
c*****************************************************************************
c
c         We must first error for the image corresponding to alphamode
          call errcalc ! calculate error with the most likely image

c         Now find the posterior probability of the image for other values 
c         of alpha and integrate the image over it.
          if (niter+19.gt.npmax) then
            write(6,*) 'main: niter+19 > npmax'
            stop
          endif
          do i = 1, 19
            iter = iter + 1
            alphaold = alpha
            alpha = alphamax*float(21-i)/20.0
            call update
            if(iprint.gt.0) write(6,"(5x,i4,3x,e13.5,3x,e13.5)") i,alpha,padm(iter)
          end do
          f(1:nf) = fbar(1:nf)/weight
          padm(1:iter) = padm(1:iter)/weight
          alpha = alphabar/weight; alphamean = alpha
          sig = sigbar/weight
          ngood = ngbar/weight
        endif
      else
c*****************************************************************************
c       Use historic MAXENT
c*****************************************************************************
        cflag=1 !set so that update does not change alpha.
c       now find alphamin and alphamax
        alphamin=0.0  
      
        call update  
        do while(aim.lt.aimaim)
          alphamin=alpha
          alpha=2.0*alpha
          call update
        end do   
        alphamax=alpha

c       Now that we have an acceptable range for alpha, we must iterate
c       to find alpha(aimaim).
        alpha=0.5*(alphamin+alphamax)
    
        do i=1,niter
          write(6,"(5x,i4,3x,e13.5,3x,f13.5)") i,alpha,aim 
          call update
      
          if(abs(aim-aimaim)/(aim+aimaim).lt.0.0001) exit
          if(aim.gt.aimaim)then
            call dchop(alpha,alphamin,alphamax)
          else
            call uchop(alpha,alphamin,alphamax)
          end if
        end do
      end if  ! (aflag.ne.4)

c     Next calculate the error and optionally perform an averaging over
c     the image if the Bryan method is not used
      if(aflag.eq.0.or.aflag.eq.2.or.aflag.eq.4) then
        call errcalc
        if(iasmc.eq.1) then ! f-space MC with constraint and 1/f metric
          call Av_spec_f
        else if(iasmc.eq.2) then ! f-space MC with constraint
          call Av_spec_f_nometric
        else if(iasmc.eq.3) then ! u-space MC with constraint
          call Av_spec_u
        else if(iasmc.eq.4) then ! f-space parallel tempering MC with constraint
          call Av_spec_f_pt
        end if
      end if
      
c     Finally rewrite the image 
      call output
      call deallocate
      
      stop
      end
