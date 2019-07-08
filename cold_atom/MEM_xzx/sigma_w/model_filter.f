      program model_filter


       implicit none
       character*72 linemc
       integer::i, nf
       double precision :: w,r1,dw



      write(6,"('   enter the input data filename:  ',$)") 
      read(5,"(a72)") linemc 
      open(unit=10,file=linemc,status='old')
      open(unit=11,file='model_good',status='unknown')

      do i=1,10000
        read(10,*,end=122) w,r1,dw
        write(11,*)w,max(1.0e-5,r1),dw
      end do
 122  nf=i-1


      close(10)
      close(11)


      end program model_filter
