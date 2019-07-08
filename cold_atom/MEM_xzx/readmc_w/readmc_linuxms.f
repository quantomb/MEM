      program readmc
!*******************************************************************************
!
!     This subroutine reads the monte carlo data and its covariance from a 
!     file.  It then generates the kernel, diagonalizes the covariance and 
!     then rotates the kernel and the data into the diagonal space of the
!     covariance.  The rotated data and kernel are written to a file in a 
!     format appropriate for the bryan MEM code.  
!
!     Depending on the average sign (the decision is made by rddcams, and
!     communicated to this code through isign=0 or 1).  if the average sign
!     is small, then the sign is included in the kernel and the covariance.
!     If the average sign is large, then each bin of data is divided by 
!     the bin average sign.  
!
!     The relationship between the data G, kernel and spectra A is given by
!     
!              /     A(w) exp(-tau w)
!     G(tau) = | dw ------------------
!              /    1 +/- exp(-beta w)    + Fermion / - Boson
!
!     where
!
!            -1/pi Im G(K,w)   Fermions
!     A(w) =
!            -1/pi Im Chi(K,w) Bosons
!
!     Since Im G(K,w) < 0, but Im Chi(K,w)<0 only for w>0.  And also since
!     A(K,w) integrates to one, while -1/pi Im Chi(K,w)/w integrates to 
!     Chi(T), we typically take (for isign=1):
!     
!              /                        D(tau) = G(tau) Fermions
!     D(tau) = | dw f(w) K(tau,w)
!              /    		        D(tau) = Chi(tau)/Chi(T) Bosons
!
!     where
!
!            -1/pi Im G(K,w)              
!     f(w) =
!            -1/pi Im Chi(K,w)/(Chi(T) w) 
!
!         /
!     1 = | dw f(w)  for both fermions and Bosons
!         /
!
!     and
!                   exp(-tau w)
!     K(tau,w) = ------------------  Fermions
!                1  + exp(-beta w)   
!
!                 w  exp(-tau w)
!     K(tau,w) = ------------------  Bosons
!                1  - exp(-beta w) 
!
!     If the average sign is small (isign=0), then we must modify this 
!     formalism to include the sign in the spectra, kernel and data. To
!     do this, we expand the data to include the bin average G times the
!     sign <Gs>, and the bin average sign as the last datum <s>.  We 
!     divide all of this by the global average sign <<s>>.  Then, for the
!     first nuse datum, 
!
!     <sD(tau)>   /    
!     --------- = | dw f(w) K(tau,w)  + error
!      <<s>>      /    
!
!     and
!
!       <s>    /    
!     ------ = | dw f(w)  + error
!      <<s>>   /    
!
!     for the last datum (i.e. K=1).
!
!     Now we take into account the summetry of the Boson spectral function
!     chi''(-w) = -chi''(w).  
!
!
!
!*******************************************************************************
!     Declare some parameters
!*******************************************************************************
      implicit none
      
      character*72 linemc
      character*10 astr

      integer i1,i2,i,j,ier,ikernel,ixvgr,nl,run,nf,nuse,
     &        nneg,l, symques, isign
      integer, parameter :: ndatamax=400,nfmax=3200

      double precision r1,r2,r3,r4,temp,beta,ta,bwidth,range,
     &     gtau(ndatamax),tau(ndatamax),data(ndatamax),
     &     R(ndatamax,nfmax),Rp(ndatamax,nfmax),w(nfmax),dw(nfmax),
     &     work(5*ndatamax),umat(ndatamax,ndatamax),
     &     cov(ndatamax,ndatamax),eigs(ndatamax),
     &     umatp(ndatamax,ndatamax),eigsp(ndatamax),
     &     Bf(nfmax,nfmax),r5,normalor

      double precision, parameter :: pi=3.1415927
!*******************************************************************************
!     This code links to the blas and lapack calls
!     dgemv, dgemm
!     dgesvd
!*******************************************************************************
!     unit    file
!     5       standard input
!     6       standard output
!     92      QMC data and covariance
!     10      image default
!     9       eigenvalue output
!*******************************************************************************

!*******************************************************************************
!     Ask a series of questions
!*******************************************************************************
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
      if(range.le.0.0) write(6,*) 'Only diagonal elements of cov will 
     &be used'

      write(6,*) '1   symmetric fermion kernel'
      write(6,*) '2   asymmetric fermion kernel'
      write(6,*) '3   symmetric boson kernel (for chipp(w)/w)'
      write(6,*) '4   symmetric boson structure factor kernel'
      write(6,*) '5   symmetric matsubara frequency kernel'
      
      write(6,"('        enter the type of kernel:  ',$)") 
      read(5,*) ikernel

!*******************************************************************************
!     Read in data from files
!*******************************************************************************
      symques=0
 401  read(92,"(a72)") linemc
      if(index(linemc(1:72),'symmetric').ne.0) symques=1
      if(index(linemc(1:72),'Results').ne.0 .and.
     &   index(linemc(1:72),' nl').ne.0   ) then
         read(92,*) nl,nuse,run,beta,normalor
         isign = 0
!         read(92,*) nuse,run,beta,normalor,isign
         temp=1.0/beta
      else
         goto 401
      endif
      if(symques.eq.1.and.ikernel.eq.2) ikernel=1
      
      if(nl.gt.ndatamax) then
        write(6,*) 'nl> ndatamax, recompile readmc'
        stop
      end if
c
      if(symques.eq.1)then
c       the dataset is symmetric and is stored in reduced form
         nuse=nl/2+1
      else
         nuse=nl
      end if
!
      read(10,*) astr,r5
      do i=1,10000
        read(10,*,end=122) w(i),r1,dw(i)
      end do
 122  nf=i-1
!
!     Check the frequency range and number of frequencies to make sure
!     that they make sense
!
      if(nf.gt.nfmax) then
        write(6,*) 'nf> nfmax, recompile readmc'
        stop
      end if
      do i=1,nf
        if(w(i).lt.0.0d0.and.ikernel.ge.3.and.ikernel.ne.5) then
           write(6,*) '**********ERROR**************'
           write(6,*) 'This Boson Kernel assumes w>0'
           write(6,*) '**********ERROR**************'
        end if
      end do

!     check for inconsistencies
      if(w(1).gt.w(nf)) then
        write(6,*) '******ERROR*********'
        write(6,*) 'ERROR, w(1) > w(nf) '
        write(6,*) '******ERROR*********'
        stop
      end if

!
!     Input the raw data
 404  read(92,"(a72)") linemc
      if(index(linemc(1:72),'Gtau(n)').ne.0) then
!         read(92,"(a72)") linemc
         do i=1,nuse
           read(92,*) tau(i),gtau(i),r1,r2,r3,r4
           gtau(i)=abs(gtau(i))
!          Note, if isign=0, gtau(nuse) is the average sign
         end do
      else
        goto 404
      end if


!*******************************************************************************
!     Now form the kernel of the transform.  Recall that the
!     chi(tau) data is the -+ data, which is 1/2 of the zz data.
!*******************************************************************************
      !isign = 1
      if(ikernel.eq.1) then
!       The kernel should reflect the symmetry of the data
        if(w(1).lt.0.0d0)then 
          do i=1,nuse-(1-isign)
            ta=tau(i)
            do j=1,nf
              if(w(j).ge.0.0d0) R(i,j)=
     & 0.5d0*(exp(-(beta-ta)*w(j))+exp(-ta*w(j)))
     & /(1.0d0+exp(-beta*w(j)))
              if(w(j).lt.0.0d0) R(i,j)=
     & 0.5d0*(exp(ta*w(j))+exp((beta-ta)*w(j)))
     & /(1.0d0+exp(beta*w(j)))
            end do
          end do
        else  ! run w(i) over positive frequencies only
          do i=1,nuse-(1-isign)
            ta=tau(i)
            do j=1,nf
              R(i,j)=(exp(-(beta-ta)*w(j))+exp(-ta*w(j)))
     & /(1.0d0+exp(-beta*w(j)))
            end do
          end do
        end if
      else if(ikernel.eq.2) then  ! asymmetric fermion kernel
        do i=1,nuse-(1-isign)
          ta=tau(i)
          do j=1,nf
            if(w(j).ge.0.0d0) R(i,j)=exp(-ta*w(j))/
     & (1.0+exp(-beta*w(j)))
            if(w(j).lt.0.0d0) R(i,j)=exp((beta-ta)*w(j))/
     & (1.0+exp(beta*w(j)))
          end do
        end do
      else if (ikernel.eq.3) then  ! symmetric boson kernel (for chipp(w)/w)
        do i=1,nuse-(1-isign)
          ta=tau(i)
          do j=1,nf
            if(w(j).gt.1.0d-5) then
              R(i,j)=w(j)*(exp(-ta*w(j)) + exp(-(beta-ta)*w(j)))/
     & (1.0-exp(-beta*w(j)))
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
            if(ABS(ta).gt.1.0d-6) then
               R(i,j) = w(j)**2/(ta**2+w(j)**2)
            else
               R(i,j) = 1.d0
            endif
          end do
        end do
      end if
      if(isign.eq.0.and.ikernel.le.2.and.w(1).lt.0.0d0) then
        R(nuse,:)=1.0
      else if(isign.eq.0.and.ikernel.le.2.and.w(1).ge.0.0d0) then
        R(nuse,:)=2.0
      else if(isign.eq.0.and.ikernel.eq.3) then
        R(nuse,:)=2.0
      end if

!     Input the covariance
 405  read(92,"(a72)") linemc
      if(index(linemc(1:72),'Gij-GiGj').ne.0) then
!         read(92,"(a72)") linemc
         do i=1,nuse
         do j=i,nuse
           read(92,*) i1,i2,cov(i,j),r1
           cov(j,i)=cov(i,j)
         end do
         end do
      else
        goto 405
      end if
      if(cov(1,1).lt.0.000001) cov(1,1)=cov(2,2)

      if(range.gt.0.0) then
!       Now diagonalize cov with SSVDC
!       SUBROUTINE SSVDC(X,LDX,N,P,S,E,U,LDU,V,LDV,WORK,JOB,INFO)
!        call ssvdc(cov,ndatamax,nuse,nuse,eigs,eigsp,umat,ndatamax,
!     &             umatp,ndatamax,work,10,ier)
        call dgesvd('A','N',nuse,nuse,cov,ndatamax,eigs,umat,ndatamax,
     &umatp,ndatamax,work,5*ndatamax,ier)
        if(ier.ne.0) write(6,*) 'dgesvd, info= ',ier
!       The singular values (eigenvalues) are in eigs arranged in 
!       decending order.  The corresponding eigenvectors are the columns
!       of umat.

!       find out how many eigenvalues to discard.
!       range must be set to the numerical precision of the computer
!       or of the covariance data; whichever is desired.
        nneg=0
        do i=1,nuse
           if(eigs(i).lt.range*eigs(1)) nneg=nneg+1
        end do

!              T                
!       Rp=Umat * R  C=AB
!       R(nuse-nneg,nf)  Rp(nuse-nneg,nf)  Umat(nuse-nneg,nuse-nneg)
!
        do i=1,nuse
        do j=1,nuse
          Umatp(i,j)=Umat(j,i)
        end do
        end do
        call dgemm('N','N',nuse-nneg,nf,nuse-nneg,1.0d0,Umatp,ndatamax,
     &             R,ndatamax,0.0d0,Rp,ndatamax)

!                                          T
!       form the transformed data data=umat *Gtau 
!       data(nuse-nneg)  Umat(nuse-nneg,nuse-nneg)  Gtau(nuse-nneg)
!
        call dgemv('N',nuse-nneg,nuse-nneg,1.0d0,Umatp,ndatamax,
     &             Gtau,1,0.0d0,data,1)

        do i=1,nuse
          eigs(i)=sqrt(abs(eigs(i))/float(run-1))
        end do
      else
!       Use only the diagonal elements of the covariance
        nneg=0
        do i=1,nuse
          data(i)=Gtau(i)
          eigs(i)=sqrt(abs(cov(i,i))/float(run-1))
          do j=1,nf
            Rp(i,j)=R(i,j)
          end do
        end do
      end if
!     
      write(96,*) nuse-nneg,nf,run,normalor
      if(ixvgr.eq.1) then
        write(9,*) '@g0 type logy'
        write(9,*) '@	 default font 0'
        write(9,*) '@	 yaxis  label "eigenvalue"'
        write(9,*) '@	 xaxis  label "index"'
        write(9,*) '@	 s0 symbol 2'
        write(9,*) '@	 s0 linestyle 0'
        write(9,*) '@	 s0 symbol size 0.5'
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

