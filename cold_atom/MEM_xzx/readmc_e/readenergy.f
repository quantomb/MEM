        program readmc
c
c       This subroutine reads the energy data file specified
c       by file.  The format of the data is assumed to be
c
c            T       E(T)     deltaE(T)
c
c       It also needs a default model for the fermi (and or 
c       bosonic) density of states.  

        implicit none
        integer:: i1,i2,i,j,ier,ikernel,ixvgr,nl,run,nf,nuse
        integer:: nneg,iblur,l
        integer nfmax,ndatamax
        parameter(ndatamax=400,nfmax=1500)
        double precision:: r1,r2,temp,beta,ta,pi,bwidth
        parameter(pi=3.1415927d0)
        double precision:: energy(ndatamax),denergy(ndatamax)
        double precision:: temperature(ndatamax)
        double precision:: kernel(ndatamax,nfmax),w(nfmax),dw(nfmax)
        double precision:: normalor
        character*72 linemc
        character*10 astr

c       unit file
c       5 standard input
c       6 standard output
c       92 QMC data and covariance
c       10 image default
c       9 eigenvalue output

      write(6,11) 
 11   format('   enter the energy data filename:  ',$)
      read(5,10) linemc 
      open(unit=92,file=linemc,status='old')

      write(6,12) 
 12   format('  enter the output data filename:  ',$)
      read(5,10) linemc 
      open(unit=96,file=linemc,status='unknown')

      write(6,13) 
 13   format('      enter the default filename:  ',$)
      read(5,10) linemc 
      open(unit=10,file=linemc,status='old')
 10   format(a72)

c     read in the default model
      read(10,*) astr,normalor
c      write(*,*) astr,normalor
      do i=1,10000
         read(10,*,end=122) w(i),r1,dw(i)
c         write(*,*) w(i),r1
      end do
 122  nf=i-1
      write(*,*) 'nf=',nf
c      pause

c
c     Input the energy data
      read(92,*) linemc
      do i=1,1000
         read(92,*,end=123) temperature(i),energy(i),denergy(i)
c         write(*,*) temperature(i),energy(i),denergy(i)
      end do
  123 nuse=i-1
      write(*,*) 'nuse=',nuse
c      pause

c     Now form the kernel of the transform.  Recall that the
c     chi(tau) data is the -+ data, which is 1/2 of the zz data.

      do i=1,nuse
         ta=temperature(i)
         do j=1,nf
            if(w(j)/ta.lt.30.0d0) then
              kernel(i,j) = w(j)/(DEXP(w(j)/ta)+1.d0)
c              write(*,*) i,j,kernel(i,j)
            else
              kernel(i,j) = 1.0d-13
c              write(*,*) i,j,kernel(i,j)
            end if
         end do
      end do

      run=30
      write(96,*) nuse,nf,run,normalor
      do 35 i=1,nuse
         write(96,*) energy(i),denergy(i),temperature(i)
 35   continue

      do 40 i=1,nuse
      do 40 j=1,nf
         write(96,*) kernel(i,j)
 40   continue

      stop
      end
