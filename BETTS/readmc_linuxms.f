      program readmc
c*******************************************************************************
c
c     This subroutine reads the monte carlo data and its covariance from a 
c     file.  It then generates the kernel, diagonalizes the covariance and 
c     then rotates the kernel and the data into the diagonal space of the
c     covariance.  The rotated data and kernel are written to a file in a 
c     format appropriate for the bryan MEM code.  
c
c     Depending on the average sign (the decision is made by rddcams, and
c     communicated to this code through isign=0 or 1).  if the average sign
c     is small, then the sign is included in the kernel and the covariance.
c     If the average sign is large, then each bin of data is divided by 
c     the bin average sign.  
c
c     The relationship between the data G, kernel and spectra A is given by
c     
c              /     A(w) exp(-tau w)
c     G(tau) = | dw ------------------
c              /    1 +/- exp(-beta w)    + Fermion / - Boson
c
c     where
c
c            -1/pi Im G(K,w)   Fermions
c     A(w) =
c            -1/pi Im Chi(K,w) Bosons
c
c     Since Im G(K,w) < 0, but Im Chi(K,w)<0 only for w>0.  And also since
c     A(K,w) integrates to one, while -1/pi Im Chi(K,w)/w integrates to 
c     Chi(T), we typically take (for isign=1):
c     
c              /                        D(tau) = G(tau) Fermions
c     D(tau) = | dw f(w) K(tau,w)
c              /                        D(tau) = Chi(tau)/Chi(T) Bosons
c
c     where
c
c            -1/pi Im G(K,w)              
c     f(w) =
c            -1/pi Im Chi(K,w)/(Chi(T) w) 
c
c         /
c     1 = | dw f(w)  for both fermions and Bosons
c         /
c
c     and
c                   exp(-tau w)
c     K(tau,w) = ------------------  Fermions
c                1  + exp(-beta w)   
c
c                 w  exp(-tau w)
c     K(tau,w) = ------------------  Bosons
c                1  - exp(-beta w) 
c
c     If the average sign is small (isign=0), then we must modify this 
c     formalism to include the sign in the spectra, kernel and data. To
c     do this, we expand the data to include the bin average G times the
c     sign <Gs>, and the bin average sign as the last datum <s>.  We 
c     divide all of this by the global average sign <<s>>.  Then, for the
c     first nuse datum, 
c
c     <sD(tau)>   /    
c     --------- = | dw f(w) K(tau,w)  + error
c      <<s>>      /    
c
c     and
c
c       <s>    /    
c     ------ = | dw f(w)  + error
c      <<s>>   /    
c
c     for the last datum (i.e. K=1).
c
c     Now we take into account the summetry of the Boson spectral function
c     chi''(-w) = -chi''(w).  
c
c
c
c*******************************************************************************
c     Declare some parameters
c*******************************************************************************
      implicit none
      
      character*72 linemc

      integer i1,i2,i,j,ier,ikernel,ixvgr,nl,run,nf,nuse,
     &        nneg,l, symques, isign
      integer, parameter :: ndatamax=800,nfmax=1600

      double precision r1,r2,temp,beta,ta,bwidth,range,
     &     gtau(ndatamax),tau(ndatamax),data(ndatamax),
     &     R(ndatamax,nfmax),Rp(ndatamax,nfmax),w(nfmax),dw(nfmax),     
     &     work(5*ndatamax),umat(ndatamax,ndatamax),
     &     cov(ndatamax,ndatamax),eigs(ndatamax),
     &     umatp(ndatamax,ndatamax),eigsp(ndatamax),
     &     Bf(nfmax,nfmax)

      double precision, parameter :: pi=3.1415927
c*******************************************************************************
c     This code links to the blas and lapack calls
c     dgemv, dgemm
c     dgesvd
c*******************************************************************************
c     unit    file
c     5       standard input
c     6       standard output
c     92      QMC data and covariance
c     10      image default
c     9       eigenvalue output
c*******************************************************************************

c*******************************************************************************
c     Ask a series of questions
c*******************************************************************************
      write(6,"('   enter the input data filename:  ',$)") 
      read(5,"(a72)") linemc 
      open(unit=92,file=linemc,status='old')

      write(6,"('  enter the output data filename:  ',$)") 
      read(5,"(a72)") linemc 
      open(unit=96,file=linemc,status='unknown')

      write(6,"('      enter the default filename:  ',$)") 
      read(5,"(a72)") linemc 
      open(unit=10,file=linemc,status='old')

      write(6,"('enter the eig. spectrum filename:  ',$)") 
      read(5,"(a72)") linemc 
      open(unit=9,file=linemc,status='unknown')
      if(index(linemc(1:72),'xvgr').ne.0.or.
     &   index(linemc(1:72),'xmgr').ne.0.or.
     &   index(linemc(1:72),'agr').ne.0) ixvgr=1

      write(6,"('  enter the range of eigenvalues:  ',$)") 
      read(5,*) range
      if(range.le.0.0) write(6,*) 'Only diagonal elements of cov will be used'

      write(6,*) '1   symmetric fermion kernel'
      write(6,*) '2   asymmetric fermion kernel'
      write(6,*) '3   symmetric boson kernel (for chipp(w)/w)'
      write(6,*) '4   symmetric boson structure factor kernel'
      write(6,*) '5   symmetric matsubara frequency kernel'
      write(6,*) '6   asymmetric matsubara frequency kernel'
      
      write(6,"('        enter the type of kernel:  ',$)") 
      read(5,*) ikernel

c*******************************************************************************
c     Read in data from files
c*******************************************************************************
      symques=0
 401  read(92,"(a72)") linemc
      if(index(linemc(1:72),'symmetric').ne.0) symques=1
      if(index(linemc(1:72),'Results').ne.0 .and.
     1   index(linemc(1:72),' nl').ne.0   ) then
         read(92,*) nl,nuse,run,beta,isign
         temp=1.0/beta
      else
         goto 401
      endif
      if(symques.eq.1.and.ikernel.eq.2) ikernel=1
      
      if(nuse.gt.ndatamax) then
        write(6,*) 'nuse> ndatamax, recompile readmc'
        stop
      end if
c
      do i=1,10000
        read(10,*,end=122) w(i),r1,dw(i)
      end do
 122  nf=i-1
c
c     Check the frequency range and number of frequencies to make sure
c     that they make sense
c
      if(nf.gt.nfmax) then
        write(6,*) 'nf> nfmax, recompile readmc'
        stop
      end if
      do i=1,nf
        if(w(i).lt.0.0d0.and.ikernel.ge.3.and.ikernel.lt.6) then
          write(6,*) '**********ERROR**************'
          write(6,*) 'This Boson Kernel assumes w>0'
          write(6,*) '**********ERROR**************'
        end if
      end do

c     check for inconsistencies
      if(w(1).gt.w(nf)) then
        write(6,*) '******ERROR*********'
        write(6,*) 'ERROR, w(1) > w(nf) '
        write(6,*) '******ERROR*********'
        stop
      end if

c
c     Input the raw data
 404  read(92,"(a72)") linemc
      if(index(linemc(1:72),'Gtau(n)').ne.0) then
c         read(92,"(a72)") linemc
         do i=1,nuse
           read(92,*) tau(i),gtau(i)
c          Note, if isign=0, gtau(nuse) is the average sign
         end do
      else
        goto 404
      end if


c*******************************************************************************
c     Now form the kernel of the transform.  Recall that the
c     chi(tau) data is the -+ data, which is 1/2 of the zz data.
c*******************************************************************************

      if(ikernel.eq.1) then
c       The kernel should reflect the symmetry of the data
        if(w(1).lt.0.0d0)then 
          do i=1,nuse-(1-isign)
            ta=tau(i)
            do j=1,nf
              if(w(j).ge.0.0d0) R(i,j)=
     &            0.5d0*(exp(-(beta-ta)*w(j))+exp(-ta*w(j)))
     &           /(1.0d0+exp(-beta*w(j)))
              if(w(j).lt.0.0d0) R(i,j)=
     &            0.5d0*(exp(ta*w(j))+exp((beta-ta)*w(j)))
     &           /(1.0d0+exp(beta*w(j)))
            end do
          end do
        else  ! run w(i) over positive frequencies only
          do i=1,nuse-(1-isign)
            ta=tau(i)
            do j=1,nf
              R(i,j)=(exp(-(beta-ta)*w(j))+exp(-ta*w(j)))
     1           /(1.0d0+exp(-beta*w(j)))
            end do
          end do
        end if
      else if(ikernel.eq.2) then  ! asymmetric fermion kernel
        do i=1,nuse-(1-isign)
          ta=tau(i)
          do j=1,nf
            if(w(j).ge.0.0d0) R(i,j)=exp(-ta*w(j))/
     &                            (1.0+exp(-beta*w(j)))
            if(w(j).lt.0.0d0) R(i,j)=exp((beta-ta)*w(j))/
     &                            (1.0+exp(beta*w(j)))
          end do
        end do
      else if (ikernel.eq.3) then  ! symmetric boson kernel (for chipp(w)/w)
        do i=1,nuse-(1-isign)
          ta=tau(i)
          do j=1,nf
            if(w(j).gt.1.0d-5) then
              R(i,j)=w(j)*(exp(-ta*w(j)) + exp(-(beta-ta)*w(j)))/
     &               (1.0-exp(-beta*w(j)))
            else
              R(i,j)=temp*(exp(-ta*w(j)) + exp(-(beta-ta)*w(j)))
            end if
          end do
        end do
      else if (ikernel.eq.4) then  ! symmetric boson structure factor kernel
        do i=1,nuse-(1-isign)
          ta=tau(i)
          do j=1,nf
            R(i,j) = exp(-ta*w(j)) + exp(-(beta-ta)*w(j))
          end do
        end do
      else if (ikernel.eq.5) then  ! symmetric matsubara freq. kernel (kinda)
        do i=1,nuse-(1-isign)
          ta=tau(i)
          do j=1,nf
            R(i,j) = 1.0/(ta**2+w(j)**2)
          end do
        end do
      else if (ikernel.eq.6) then  ! assymmetric matsubara freq. kernel
        do i=1,nuse
	  ta=tau(i)
	  if(ta.gt.0) then
            do j=1,nf
              R(i,j) = -w(j)/(ta**2+w(j)**2)
            end do
	  else
	    tau(i)=-tau(i)
            do j=1,nf
              R(i,j) = ta/(ta**2+w(j)**2)
            end do
	  end if
        end do
      end if
      if(isign.eq.0.and.ikernel.le.2.and.w(1).lt.0.0d0) then
        R(nuse,:)=1.0
      else if(isign.eq.0.and.ikernel.le.2.and.w(1).ge.0.0d0) then
        R(nuse,:)=2.0
      else if(isign.eq.0.and.ikernel.eq.3) then
        R(nuse,:)=2.0
      end if

c     Input the covariance
 405  read(92,"(a72)") linemc
      if(index(linemc(1:72),'Gij-GiGj').ne.0) then
c         read(92,"(a72)") linemc
         do i=1,nuse
         do j=i,nuse
           read(92,*) i1,i2,cov(i,j)
           cov(j,i)=cov(i,j)
         end do
         end do
      else
        goto 405
      end if
      if(cov(1,1).lt.0.000001) cov(1,1)=cov(2,2)

      if(range.gt.0.0) then
c       Now diagonalize cov with SSVDC
c       SUBROUTINE SSVDC(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,JOB,INFO)
c        call ssvdc(cov,ndatamax,nuse,nuse,eigs,eigsp,umat,ndatamax,
c     1             umatp,ndatamax,work,10,ier)
        call dgesvd('A','N',nuse,nuse,cov,ndatamax,eigs,umat,ndatamax,umatp,
     &             ndatamax,work,5*ndatamax,ier)
        if(ier.ne.0) write(6,*) 'dgesvd, info= ',ier
c       The singular values (eigenvalues) are in eigs arranged in 
c       decending order.  The corresponding eigenvectors are the columns
c       of umat.

c       find out how many eigenvalues to discard.
c       range must be set to the numerical precision of the computer
c       or of the covariance data; whichever is desired.
        nneg=0
        do i=1,nuse
           if(eigs(i).lt.range*eigs(1)) nneg=nneg+1
        end do

c              T                
c       Rp=Umat * R  C=AB
c       R(nuse-nneg,nf)  Rp(nuse-nneg,nf)  Umat(nuse-nneg,nuse-nneg)
c
        do i=1,nuse
        do j=1,nuse
          Umatp(i,j)=Umat(j,i)
        end do
        end do
        call dgemm('N','N',nuse-nneg,nf,nuse-nneg,1.0d0,Umatp,ndatamax,
     &             R,ndatamax,0.0d0,Rp,ndatamax)

c                                          T
c       form the transformed data data=umat *Gtau 
c       data(nuse-nneg)  Umat(nuse-nneg,nuse-nneg)  Gtau(nuse-nneg)
c
        call dgemv('N',nuse-nneg,nuse-nneg,1.0d0,Umatp,ndatamax,
     &             Gtau,1,0.0d0,data,1)

        do i=1,nuse
          eigs(i)=sqrt(abs(eigs(i))/float(run-1))
        end do
      else
c       Use only the diagonal elements of the covariance
        nneg=0
        do i=1,nuse
          data(i)=Gtau(i)
          eigs(i)=sqrt(abs(cov(i,i))/float(run-1))
          do j=1,nf
            Rp(i,j)=R(i,j)
          end do
        end do
      end if
c     
      write(96,*) nuse-nneg,nf,run
      if(ixvgr.eq.1) then
        write(9,*) '@g0 type logy'
        write(9,*) '@    default font 0'
        write(9,*) '@    yaxis  label "eigenvalue"'
        write(9,*) '@    xaxis  label "index"'
        write(9,*) '@    s0 symbol 2'
        write(9,*) '@    s0 linestyle 0'
        write(9,*) '@    s0 symbol size 0.5'
        write(9,*) '@TYPE xy'
      end if
      do i=1,nuse-nneg
         write(96,*) data(i),eigs(i)
         write(9,*) i,eigs(i)
      end do

      do i=1,nuse-nneg
      do j=1,nf
        write(96,*) Rp(i,j)
      end do
      end do

      stop
      end
