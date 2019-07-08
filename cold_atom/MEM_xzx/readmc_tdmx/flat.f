      program default
c code for default model 
      implicit none
      integer*4,parameter :: nfmax=800
      real*8,parameter :: PI=3.1415926535897932354d0
      real*8 :: w(nfmax),f(nfmax),dw(nfmax)
      integer*4 :: i,j,k
      real*8 :: w_tmp, dw_tmp
      
c
      w_tmp = 0.d0
      dw_tmp = 2.d0*PI/DBLE(nfmax)
      open(99,file='deflt.dat',status='unknown')
      open(96,file='f0.dat',status='unknown')
      do i = 1, nfmax
         w(i) = w_tmp
         f(i) = 1.d0
         dw(i) = dw_tmp
         w_tmp = w_tmp + dw_tmp
         write(99,*) w(i),f(i),dw(i)
         write(96,*) w(i),f(i)
      enddo
c
      close(99)
      close(96)
      stop
      end
