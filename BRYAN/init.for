      subroutine readin
c*******************************************************************************

c     This block initializes the MaxEnt code by reading in the problem 
c     dependent variables, and generating the problem dependentb arrays.
c*******************************************************************************

      use Global
      implicit none

      character *80 notused
      character*10 astr
      integer i, j, l, nsing, info, ios
      double precision r1, cop, co, normm, normf,r5

c*******************************************************************************


c     input data contains rotated data, rotated kernel (see README for format)
      write(6,"('   enter the input data filename:  ',$)") 
      read(5,"(14a)") notused 
      open (unit=11,file=notused,status='old')

c     read the filename of file containting the default model
      write(6,"('      enter the default filename:  ',$)") 
      read(5,"(14a)") notused 
      open (unit=2,file=notused,status='old')

c     read the filename of file containting the initial image 
c     (redundant with default)
      write(6,"('enter the initial image filename:  ',$)") 
      read(5,"(14a)") notused 
      open (unit=3,file=notused,status='old')

      write(6,"('  enter the final image filename:  ',$)") 
      read(5,"(14a)") notused 
      open (unit=7,file=notused,status='unknown')
      if(index(notused(1:72),'agr' ).ne.0.or.
     &   index(notused(1:72),'xvgr').ne.0.or.
     &   index(notused(1:72),'xmgr').ne.0) ixvgr=1
      write(6,"('                     enter niter:  ',$)") 
c     max number of alphas allowed, typically << 1000 - more than this 
c     likely means that MEM is not working
      read(5,*) niter 
      if (niter.gt.npmax) then
         write(6,*) 'main: niter > npmax'
         stop
      endif

      write(6,"('                  enter nregions:  ',$)") 
c     cannot show error bars on calculated profile (SSI does though)
c     need to divide "image" into regions and code will calc error bars on 
c     total weight in each region error bar on a single point is infinity 
c     (add more discussion on this)
      read(5,*) nregions 

      write(6,"('                     enter alpha:  ',$)") 
c     intial guess at alpha - should start with too large a value for classic 
c     or historic for Bryan MUST take init alpha greater than alpha max so 
c     100 is good start (this will be out on the asymmetric tail of alpha 
c     prob dist
      read(5,*) alphainit 
      alpha=alphainit

c 
c     SET THE MEM METHOD
c
      write(6,"('Set the MEM Method')") 
      write(6,"('classicjp for classic  MEM w/Jeffrey Prior')") 
      write(6,"('bryanjp   for bryan    MEM w/Jeffrey Prior')") 
      write(6,"('classic   for classic  MEM wo/Jeffrey Prior')") 
      write(6,"('bryan     for bryan    MEM wo/Jeffrey Prior')") 
      write(6,"('historic  for historic MEM')") 

      write(6,"('                    enter method:  ',$)") 
      read(5,*) method 
      read(method,"(i3)",iostat=ios) aflag   ! should barf with characters
      write(6,*) 'method=',method
      if (ios.gt.0) then
c       aflag   method          function
c       0       classicjp       Classic with Jeffrey Prior for tighter fit.
c       1       bryanjp         Bryan with Jeffrey Prior for tighter fit.
c       2       classic         Classic
c       3       bryan           Bryan            
c       4       historic        Historic MEM    image fit to chi^2 = ndata
             if(method.eq.'classicjp')then; aflag=0
        else if(method.eq.'bryanjp')  then; aflag=1
        else if(method.eq.'classic')  then; aflag=2
        else if(method.eq.'bryan')    then; aflag=3
        else if(method.eq.'historic') then; aflag=4
        else
        write(6,*) 'invalid method'
        stop
        end if
      end if
      write(6,*) 'aflag=',aflag

c 
c     SET THE SSI METHOD
c
      if(aflag.eq.0.or.aflag.eq.2) then
        write(6,"('Set the SSI Method')") 
        write(6,"('none      for no SSI')") 
        write(6,"('ssimetric for SSI with the 1/f metric')") 
        write(6,"('ssi       for SSI without the 1/f metric')") 
        write(6,"('                    enter method:  ',$)") 
        read(5,*) methodssi 
        read(methodssi,"(i3)",iostat=ios) iasmc   ! should barf with characters
        if (ios.gt.0) then
c         iasmc   methodssi     function
c         0       none          no SSI.
c         1       ssimetric     SSI with the 1/f metric.
c         2       ssi           SSI without the 1/f metric
          if(methodssi.eq.'none')     then; iasmc=0
          else if(methodssi.eq.'ssimetric')  then; iasmc=1
          else if(methodssi.eq.'ssi')  then; iasmc=2
          else
          write(6,*) 'ERROR invalid ssi methodssi'
          stop
          end if
        end if
        if(iasmc.ge.1) then
          write (6,*)'enter nsweeps,nruns,nwarms,idummy,delu'
          read (5,*) nsweeps,nruns,nwarms,idummy,delu
          write(6,"('     enter the SSI data filename:  ',$)") 
          read(5,"(14a)") notused 
          open (unit=88,file=notused,status='unknown')
        end if
      else
        iasmc=0
        if(aflag.eq.4) then
          write(6,"('                       enter aim:  ',$)") 
          read(5,*) aimaim
        end if
      end if

c     print level
      write(6,"('                    enter iprint:  ',$)") 
      read(5,*) iprint 
      write(6,*) 'done1'

c     Read in the data file header 
      read(11,"(72a)") notused
      if(index(notused(1:72),'shift' ).ne.0) then
        read(notused,*) ndata,nf,r1,shift
      else
        read(notused,*) ndata,nf,r1,normalor
        shift=zeror
      end if
c     Now that we know ndata, nf, and nregions
      call allocate

c     Read in the data and error.  Note it is assumed that the covariance
c     of the data is diagonal, or that the covariance has been diagonalized
c     and that the data, covariance, and kernel have been transformed into
c     this space.

      do i=1,ndata
         read(11,*) data(i),error(i)
         if(error(i).eq.0.0d0) then
           write(6,*) '***************************************'
           write(6,*) '***** ERROR an error bar is zero ******'
           write(6,*) '***************************************'
           stop
         else if(abs(error(i)/data(i)).lt.1.0d-6) then
           write(6,*) 'WARNING an error bar is suspiciously small'
         end if      
      end do

c     Read in the kernel  
      do i=1,ndata
         do j=1,nf
            read (11,*) T(i,j)
         end do
      end do

c     Read in the default model
      read(2,*) astr,r5
      do i=1,nf
         read (2,*) r1,model(i),dw(i)
         model(i) = model(i) + shift
         if(abs(model(i)).lt.10d-25) then
            write(6,*) 'WARNING model=zero!'
            model(i)=10e-25
         endif
       end do

c     Read in the frequency grid and present image f.
      read(3,*) astr,r5
      do i=1,nf
         read (3,*) w(i),f(i)
         f(i) = f(i) + shift
         if(abs(f(i)).lt.10d-25) then
            write(6,*) 'WARNING f(i)=zero!'
            f(i)=10e-25
         endif
      end do


c     normalize the image and model so that their sums are integrals.
      normm=sum(model(1:nf)*dw(1:nf))
      normf=sum(    f(1:nf)*dw(1:nf))
      write(6,"('               norm of the model:  ',e9.4)") normm
      write(6,"('       norm of the initial image:  ',e9.4)") normf

c     Now normalize the model and initial image so that their sums are
c     integrals
      do i=1,nf
        mr(i)   = model(i)
        f(i)    = f(i)*dw(i)
        model(i)= mr(i)*dw(i)
      end do

c     Calculate the chi-squared for the default model.
c     First form the reconstructed data Tf(ndata)=T(ndata,nf)*model(nf)
      co=0.0
      cop=0.0
      call dgemv('N',ndata,nf,1.0d0,T,ndata,model,1,0.0d0,Tf,1)
      do l=1,ndata
         co=co+((Tf(l)-data(l))/error(l))**2
         cop=cop+abs(Tf(l)-data(l))
      end do
      co=co/float(ndata)
      cop=cop/float(ndata)
      if(iprint.gt.0) then
        write (6,"('       default model Chisq/ndata:  ',e10.4)")  co
        write (6,"('                 net discrepancy:  ',e10.4)")  cop
      end if

      return
      end subroutine readin
      
      
      subroutine init
c*******************************************************************************

c     This block initializes the MaxEnt code by generating the problem
c     dependentb arrays.
c*******************************************************************************

      use Global
      implicit none

      integer i, j, l, nsing, info, ios
      double precision r1, e(nf)
      double precision, allocatable :: temp1(:,:),temp2(:,:),
     &                                 temp3(:,:),VT(:,:), work(:)

c*******************************************************************************
      allocate (temp1(ndata,nf))
      allocate (temp2(ndata,nf))
      allocate (temp3(nf,ndata))
      allocate (VT(nf,ndata))
c*******************************************************************************
c     Perform singular value decomposition of the kernel T.  First form
c     the transpose of the kernel T.  Here
c
c     T = V Sigma U[T]
c
c     T[T] = U Sigma V[T]  
c
c     where Sigma is a diagonal matrix composed of the singular values of T,
c     which will be stored in the vector s.
c
      do i=1,ndata
      do j=1,nf
        temp3(j,i) = T(i,j)
      end do
      end do
c     temp3(nf,ndata) = T[T]
      call dgesvd('A','A',nf,ndata,temp3,nf,s,Umat,nf,VT,nf,
     &            r1,-1,info)
      allocate (work(int(r1)+1))
      call dgesvd('A','A',nf,ndata,temp3,nf,s,Umat,nf,VT,nf,
     &            work,int(r1)+1,info)
      deallocate (work)
      if(info.ne.0) then
        write(6,*) 'SVD error in init'
        write(6,*) 'info=',info
        stop
      end if
      write(6,*) 'done svd'
      do i=1,nf
      do j=1,ndata
        V(j,i)=VT(i,j)
      end do
      end do
c     Determine the size of the singlular space, ns.  The singular space
c     includes all S(i) > S(1) * range, where S(1) is the largest 
c     singular value. 
      nsing=0
      do i=1,nf
         if (s(i).lt.range*s(1)) then
           nsing=nsing+1
c          s(i)=0.0
         end if
      end do
      ns = nf-nsing
      if(iprint.gt.0) 
     & write (6,"('                              ns:  ',i2)")  ns
c     
c     Form M = sigma*V[T]*W*V*sigma[T]
      do i = 1, ndata
         do j = 1, ns
            temp1(i,j) = V(i,j)*s(j)
            temp2(i,j) = temp1(i,j)/error(i)**2
         end do
      end do
c     temp1[T] * temp2 = xm
      call dgemm('T','N',ns,ns,ndata,1.0d0,temp1,ndata,temp2,ndata,
     &           0.0d0,xm,nf)


c     Generate from the image f the initial guess for u
c     u = U[T]*log(f/model) u(ns)=Umat(ns,nf)e(nf)
      e(1:nf) = dlog(f(1:nf)/model(1:nf))
      call dgemv('T',nf,ns,1.0d0,Umat,nf,e,1,0.0d0,u,1)
    
c     initialize various sums
      weight = 0.0
      alphabar = 0.0
      sigbar = 0.0
      ngbar = 0.0
      iuflow=0
      fbar(:) = 0.0

c     other initializations
      cflag=0
c        
c*******************************************************************************
      deallocate (temp1)
      deallocate (temp2)
      deallocate (temp3)
      deallocate (VT)
c*******************************************************************************
      return
      end subroutine init


