        module Global
c*******************************************************************************
c*****************************************************************************
c       This is the global definitions block for the MaxEnt
c       code.  All global variables must be declared here.
c*****************************************************************************
c       SOME GLOBAL VARIABLES
c*****************************************************************************
c       alpha       Lagrange multiplier which relates C and S
c       alphainit   The initial value of alpha
c       alphamode   The value of alpha that maximizes P(alpha|D,m,...)
c       alphamean   The average alpha over P(alpha|D,m,...)
c       alphaold    The value of alpha from the previous iterate
c       alphadamp   alphanew=alphadamp*alphanew + (1-alphadamp)*alphaold
c
c       method      aflag  string that determines how the image is found
c       ----------------------------------------------------------------
c       classicJP     0    Classic with Jeffrey Prior for tighter fit.
c       bryanJP       1    Bryan with Jeffrey Prior for tighter fit.
c       classic       2    Classic without the Jeffrey Prior
c       bryan         3    Bryan without the Jeffrey Prior  
c       historic      4    Historic MEM image fit to chi^2 = ndata * aim
c       
c       methodssi   iasmc if >0 average spectrum MC done if aflag=0,2
c       ----------------------------------------------------------------
c       none          0    no SSI.
c       ssimetric     1    SSI with the 1/f metric.
c       ssi           2    SSI without the 1/f metric
c
c       data        raw data to be deconvolved
c       ndata       number of data points
c       error       standard error bars associated data

c       f           image array (actually f*dw )
c       nf          number of image pixels
c       w           image frequencies
c       dw          frequency step
c       model       default model values of f

c       T           kernel which relates f and data (data=t*f)
c       ns          size of singular space of T 

c       post        current value (unnormalized ) of P(alpha|m,D,..)
c       du          increment in the solution
c       width       The width of an error calc. region
c       lo          current value of L
c       mu          step-size cut-off parameter
c       niter       total number of iterations of alpha
c       precision   assumed numerical precision of the code
c       sig         error scaling parameter
c       so          current value of S
c       step        multiple of sum(f) used to define max step los
c       test1       test constant as defined by alpha*S+G
c       testf       defines an acceptable value for test
c       iuflow      set to 1 if the image falls below e^-60
c       iprint      print level     
c       ixvgr       if 1, then print out xvgr graphics directives
c       shift       the shift in the image needed to make anomalous
c                   green function spectra positive definite
c*****************************************************************************
c*****************************************************************************
        implicit none
        save

c       Global declarations for the program
c       in the order: character, integer, real, then complex

c       *****************************CHARACTER********************************
        character*16 :: method, methodssi

c       *****************************INTEGER**********************************

        integer, parameter :: npmax=3000         ! max num of alpha iterations
        integer, parameter :: iupdate_max=1000   ! max num of Newton steps
        integer, parameter :: icontrol_max=1000  ! max num of Newton steps
        integer :: iter, niter, ndata, nf, ns, aflag, cflag, nregions,  
     &             iuflow, optiono, iasmc, nsweeps, nwarms,
     &             nruns, idummy, iprint, ixvgr
     
     
c       ***********************REAL, DOUBLE PRECISION**************************

        double precision, parameter ::        pi = 3.1415927
        double precision, parameter :: alphadamp = 0.150
        double precision, parameter :: precision = 1.0d-5
        double precision, parameter ::      step = 0.25
        double precision, parameter ::     range = 1.0d-6
        double precision, parameter ::     testf = 1.0d-5
        double precision, parameter ::     zeror = 0.0d0

        double precision, allocatable :: T(:,:), data(:), error(:),
     &     w(:), dw(:), Tf(:), f(:), model(:), u(:), Umat(:,:), s(:),
     &     V(:,:), xm(:,:), mr(:), du(:), lambda(:), fbar(:),
     &     fave(:), faves(:), fsdev(:), P(:,:), Z(:), R(:,:),
     &     Gvalue(:), Dgvalue(:), low(:), high(:), width(:)
     
        double precision ::
     &     padm(npmax), alpt(npmax), mu, alpha, sig, sigbar, 
     &     alphabar, alphaold, alphainit, alphamode, alphamean, 
     &     alphamax, alphamin, alphaave, alphaaves, alphasdev,
     &     weight, lo, so, post, sbar, dssbar, ngood, ngbar,
     &     aim, aimaim, delu, shift,normalor

        end module Global



        subroutine allocate
c*******************************************************************************
c       allocate some arrays
c*******************************************************************************
        use Global
        implicit none
        integer info,infot
c*******************************************************************************

        infot=0
        allocate (T(ndata,nf),stat=info); infot=infot+info            
        allocate (data(ndata),stat=info); infot=infot+info
        allocate (error(ndata),stat=info); infot=infot+info
        allocate (w(nf),stat=info); infot=infot+info
        allocate (dw(nf),stat=info); infot=infot+info
        allocate (Tf(ndata),stat=info); infot=infot+info
        allocate (f(nf),stat=info); infot=infot+info
        allocate (fbar(nf),stat=info); infot=infot+info
        allocate (fave(nf),stat=info); infot=infot+info
        allocate (faves(nf),stat=info); infot=infot+info
        allocate (fsdev(nf),stat=info); infot=infot+info
        allocate (model(nf),stat=info); infot=infot+info
        allocate (u(nf),stat=info); infot=infot+info
        allocate (du(nf),stat=info); infot=infot+info
        allocate (lambda(nf),stat=info); infot=infot+info
        allocate (Umat(nf,nf),stat=info); infot=infot+info
        allocate (P(nf,nf),stat=info); infot=infot+info
        allocate (R(nf,nf),stat=info); infot=infot+info
        allocate (s(nf),stat=info); infot=infot+info
        allocate (Z(nf),stat=info); infot=infot+info
        allocate (V(ndata,nf),stat=info); infot=infot+info
        allocate (xm(nf,nf),stat=info); infot=infot+info
        allocate (mr(nf),stat=info); infot=infot+info
        allocate (Gvalue(nregions),stat=info); infot=infot+info
        allocate (Dgvalue(nregions),stat=info); infot=infot+info
        allocate (low(nregions),stat=info); infot=infot+info
        allocate (high(nregions),stat=info); infot=infot+info
        allocate (width(nregions),stat=info); infot=infot+info
       
        if(infot.ne.0) then
          write(6,*) 'ERROR:  allocate',infot,' failures'
          stop
        end if

        return
        end subroutine allocate
        


        subroutine deallocate
c*******************************************************************************
        use Global
        implicit none
c*******************************************************************************

        deallocate (T)
        deallocate (data)
        deallocate (error)
        deallocate (w)
        deallocate (dw)
        deallocate (Tf)
        deallocate (f)
        deallocate (fbar)
        deallocate (fave)
        deallocate (faves)
        deallocate (fsdev)
        deallocate (model)
        deallocate (u)         
        deallocate (du)         
        deallocate (lambda)         
        deallocate (Umat)         
        deallocate (P)         
        deallocate (R)         
        deallocate (s)         
        deallocate (Z)         
        deallocate (V)         
        deallocate (xm)         
        deallocate (mr)         
        deallocate (Gvalue)         
        deallocate (Dgvalue)         
        deallocate (low)         
        deallocate (high)         
        deallocate (width)         

        return
        end subroutine deallocate
