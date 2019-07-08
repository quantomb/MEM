      subroutine output
c*******************************************************************************
c
c     This program prints out the final MEM spectrum to unit 7
c     (specified in init).  If ixvgr=1, then it also includes
c     AceGR formatting statements (xmgrace uses these formatting 
c     statment to label the graphs).  The residual is printed to
c     file residual.xvgr (unit 29).  If the residual is gaussianly
c     distributed around one, then the information in the data has
c     been exhausted.

c*******************************************************************************
      use Global
      implicit none

      integer i, l
      double precision co, r1

c**********************************************************************
c     Form chi-squared and the residual of the image, and print this to
c     residual.xvgr together with formatting.  The residual is the point
c     by point deviation of the reconstructed data compared to the 
c     orginal data.  You want this quantity to be uncorrelated as a
c     function of the data index l.
      open(unit=29,file='residual.xvgr',status='unknown')
      write(29,*) '@    default font 0'
      write(29,*) '@    yaxis  label "residual"'
      write(29,*) '@    xaxis  label "l"'
      write(29,*) '@    s0 symbol 2'
      write(29,*) '@    s0 linestyle 0'
      write(29,*) '@    s0 symbol size 0.5'
      write(29,*) '@TYPE xy'
      co=0.0
c     Form the reconstructed data Tf(ndata)=T(ndata,nf)*f(nf)
      call dgemv('N',ndata,nf,1.0d0,T,ndata,f,1,0.0d0,Tf,1)
      do l=1,ndata
        co=co+((Tf(l)-data(l))/error(l))**2
        write(29,*) l,(Tf(l)-data(l))/error(l) ! This is the residual 
      end do
      aim=co/float(ndata)
      write(29,*) '&'
      write(29,*) '@TYPE xy'
      write(29,*) '0.0   0.0 '
      write(29,*) ndata,'   0.0 '
      write(29,*) '&'
c**********************************************************************


c**********************************************************************
c     Renormalize the image for printing, and print it to unit 7, 
c     which is formatted for xmgrace if ixvgr=1
c     
c     Now renormalize the model and image so that they can be plotted
      r1=0.0d0
      do i=1,nf
        r1=r1+f(i)
        u(i)=f(i)
        f(i)=u(i)/dw(i)
        model(i)=mr(i)
        u(i)=f(i)/model(i)
      end do

      if(aflag.eq.2.or.aflag.eq.0) weight=2.0*weight
c     Write the present image
      if(ixvgr.eq.1) then
        write(7,*) '#weight: ',weight,' aflag:',aflag
        write(7,*) '# Ngood: ',ngood
        write(7,*) '#alpha: ',alpha
        write(7,*) '#sigma: ',sig,' aim:',aim
        write(7,*) '@    s2 symbol 2'
        write(7,*) '@    s2 linestyle 0'
        write(7,*) '@    s2 symbol size 0.5'
        write(7,*) '@    legend on'
        write(7,*) '@    legend x1 0.95'
        write(7,*) '@    legend y1 0.80'
        write(7,*) '@    legend string 0 "image"'
        write(7,*) '@    legend string 1 "model"'
        write(7,*) '@    legend string 2 "errors"'
        write(7,*) '@TYPE xy'
      end if
      open(unit=99,file='dosQMC',status='unknown')
      write(99,*) '#normalor:',normalor
      do i=1,nf
        write (7,"(f14.8,', ',e16.9,', ',f14.8)") w(i),f(i)-shift,dw(i)
        write (99,"(f14.8,2x,e16.9,2x,f14.8)") w(i),
     &max(1.0e-7,f(i))-shift,dw(i)
c       write (98,"(f14.8,2x,e16.9,2x,f14.8)") w(i),max(1.0e-5,f(i))-shift,dw(i)
      end do
      close(99)
      if(ixvgr.eq.1) then
        write(7,*) '&'
        write(7,*) '@TYPE xy'
      else
        write(7,*) '   '
      end if
         do i=1,nf
          write(7,"(1x,e10.4,', ',e10.4)") w(i),model(i)-shift
        end do
      if(ixvgr.eq.1) then
        write(7,*) '&'
        write(7,*) '@TYPE xydxdy'
        do i=1,nregions
          write(7,"(1x,f7.3,', ',f8.4,', ',f7.3,', ',e10.4,', ',
     & f7.3,', ',f7.3)")
     & 0.5*(low(i)+high(i)),Gvalue(i),width(i),dGvalue(i),low(i),high(i)
        end do        
      end if
      

c*****************************************************************************
c     if the Bryan method was used, print the posterior probability of
c     alpha P(alpha|m,G,...) to the xmgrace formated file posterior.xvgr
  
      if(aflag.eq.1.or.aflag.eq.3) then
         open (unit=8,status='unknown',file='posterior.xvgr')
         write(8,"('@    default font 0')") 
         write(8,90) 
 90      format('@    yaxis  label "P(alpha|D,m...)"')
         write(8,91) 
 91      format('@    xaxis  label "alpha"')
         do i = 1,iter
            if(padm(i).gt.1.0d-60) write (8,"(1x,e10.4,', ',e10.4)") 
     &alpt(i),padm(i)
         end do
         close(unit=8)
      endif

c*****************************************************************************
c     Now print some useful information to standard output
      write(6,*) ' '
      if(iuflow.eq.1) write(6,*) 'An image pixel fell below e^-60 '
      write(6,*) '        Norm: ', r1
      write(6,*) '       Alpha: ', alpha
      write(6,*) '       Sigma: ', sig
      write(6,*) '       Ngood: ', ngood
      write(6,*) '       Chisq: ', co
      write(6,*) ' Chisq/ndata: ', aim
      if(aflag.eq.1.or.aflag.eq.3) then
        write(6,*) '    weight: ',weight
        write(6,*) 'mean alpha: ',alphamean
        write(6,*) 'mode alpha: ',alphamode
      end if
c
      return
      end
