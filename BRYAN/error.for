      subroutine errcalc
c*******************************************************************************
c     This code calculates the integral of the image in a set of nregions
c     regions, together with the associated error bars.  The integrand may 
c     be multiplied by a function g defined at the bottom of this file.
c
c     In order to calculate the error, we need the covariance of the
c     image.  In the gaussian approximation, this would be
c
c           T	      -1		      -1
c     <df df > ~= -ddQ   = (- alpha ddS - ddL)
c
c       							  T T  
c              = {f}/alpha -{f} U Y {lambda/alpha(alpha+lambda)} Y U {f}
c
c       				   T  T
c              = {f}U Y {1/(alpha+lambda)}Y  U  {f}
c
c     The error involved in the integral
c
c             /
c     G   =   | df f(w) g(w)  
c             /
c
c     May be approximated by  < (delta G)^2 >
c
c       		//
c     < (delta G)^2 > = || dx dy g(x) g(y) <delta f(x) delta f(y)>
c       		//
c
c*******************************************************************************
      use Global
      implicit none
      
      double precision g
      double precision, allocatable :: temp(:,:), ddQinv(:,:), Y(:,:)
      integer ilow, ihigh, i, j, k, iregion
c*******************************************************************************
      allocate (temp(nf,nf))
      allocate (ddQinv(nf,nf))
      allocate (Y(nf,nf))
c*******************************************************************************
c
c     Form Y(ns,ns)   Y= P*{1/sqrt(z)}*R
c
      do i=1,ns
      do j=1,ns
        Y(i,j)=0.0
      do k=1,ns
         Y(i,j)=Y(i,j) + P(i,k)*R(k,j)/sqrt(Z(k))
      end do
      end do
      end do
c
c     Form the matrix temp(ns,ns) = Y*{1/(alpha+lambda)}*Y[T]
c      
      do i=1,ns
      do j=1,ns
         temp(i,j)=0.0
      do k=1,ns
         temp(i,j)=temp(i,j) + Y(i,k)*Y(j,k)/(alpha+lambda(k))
      end do
      end do
      end do

c     Now form the product U*temp*U[T]
c
c     first form U*temp(nf,ns) in Y; Y(nf,ns) =  U(nf,ns) temp(ns,ns)
c
      call dgemm('N','N',nf,ns,ns,1.0d0,Umat,nf,temp,nf,
     .  	 0.0d0,Y,nf)
c
c     second, form U*temp*U[T] in ddQinv
c
      call dgemm('N','T',nf,nf,ns,1.0d0,Y,nf,Umat,nf,
     .  	 0.0d0,ddQinv,nf)
c
c     Form {f}*U*Y*{1/(alpha+lambda)}*Y[T]*U[T]*{f} which is a (nf,nf)
c     matrix, or Y*U[T] in ddQinv.  This should be the covariance of 
c     the image <delta f_i delta f_j> = ddQinv(i,j).
      do i=1,nf
        do j=1,nf
          ddQinv(j,i)=f(j)*ddQinv(j,i)*f(i)
        end do
      end do
c
c     Now split the region from w(1) to w(nf) into nregions regions.
c     This will be done by dividing nf by nregions.
c
      do iregion=1,nregions

        low(iregion)=w(1)+(iregion-1)*(w(nf)-w(1))/float(nregions)
        do i=1,nf
         if (w(i).ge.low(iregion)) exit
        end do
        ilow=i
        low(iregion)=w(ilow)-dw(ilow)*0.5d0

        high(iregion)=w(1)+(iregion)*(w(nf)-w(1))/float(nregions)
        do i=1,nf
         if (w(i).gt.high(iregion)) exit
        end do
        ihigh=i-1
	high(iregion)=w(ihigh)+dw(ihigh)*0.5d0

        Gvalue(iregion)=0.0d0
        do i=ilow,ihigh
          Gvalue(iregion)=Gvalue(iregion)+f(i)-(shift*dw(i))*g(w(i))
        end do

        dGvalue(iregion)=0.0d0
        do i=ilow,ihigh
        do j=ilow,ihigh
          dGvalue(iregion)=dGvalue(iregion)+g(w(i))*g(w(j))*ddQinv(i,j)
        end do
        end do
        dGvalue(iregion)=sqrt(dGvalue(iregion))

        width(iregion)=(high(iregion)-low(iregion))*0.5d0 ! the half width of the region

      end do
      open(unit=77,file='error.dat',status='unknown')
      write(77,*) '#            /'
      write(77,*) '#    G   =   | df f(w) g(w)'
      write(77,*) '#            /'
      write(77,*) '# '
      write(77,*) '#  w_ave       G         dG        w_low   w_high'
      do iregion=1,nregions
        write(77,"(1x,f7.3,2x,e10.4,2x,e10.4,2x,f7.3,2x,f7.3,2x)") 
     &             0.5*(low(iregion)+high(iregion)),Gvalue(iregion),
     &             dGvalue(iregion),low(iregion),high(iregion)
      end do
      write(77,*) '#  '
      write(77,*) '#  the image at the mode'
      write(77,*) '#    w         f(w) '
      do i=1,nf
        write(77,"(1x,e10.4,2x,e10.4)") w(i),f(i)
      end do
c
c*******************************************************************************
      deallocate (temp)
      deallocate (ddQinv)
      deallocate (Y)
c*******************************************************************************

      return
      end


      double precision function g(x)
c     Define the function to be masked.  For example, if we
c     wanted to calculate the first moment of f, then g=x
c     or x**2 for the second moment. etc.
      double precision x

      g=1.0d0

      return
      end



      subroutine Av_spec_u
c
c*******************************************************************************
c     This subroutine implements a stochastic sampling of the spectra
c     as described by A. Sandvik.  However, it differs from Anders method
c     in that alpha (his temperature) is taken from classic MEM.  Thus,
c     this method only works with classic mem, aflag=0 or 2.  
c
c     The stochastic search is done in the singular space.  Since ns<<nf,
c     this is very efficient.  The image f is related to u through
c
c     D= Tf        T=V Sigma U^T       f=m exp(Uu)
c
c*******************************************************************************
      use Global
      implicit none
      integer isweep,i,j,iu,l,irun,ialpha,icount
      real ran2
      double precision fprop(nf),fnow(nf),unow(nf),uprop(nf),
     &                 fv1(nf),acrat(nf),ac(nf),normf,
     &                 loprop,lonow,loave,r1,alphanow,normprop,normnow

c*******************************************************************************

      r1=ran2(-idummy) ! see the random number generator
      write(88,*) '#  From  subroutine Av_spec_u'

c     Initiate some accumulators and variables
      unow(1:ns)=u(1:ns); uprop(1:ns)=u(1:ns)
      fnow(1:nf)=f(1:nf); normnow=sum(f(1:nf)); normf=normnow
      lonow=lo

      do ialpha=11,1,-2 ! simulated annealing in alpha
        alphanow=alpha*dble(ialpha)
        write(6,*) '*************************************************'
        write(6,*) ' alphanow= ',ialpha,'*alpha= ',alphanow
        write(6,*) '*************************************************'
        acrat=0.0d0; ac=0.0d0; fave=0.0d0; faves=0.0d0
        loave=0.0d0
        do irun=1,nruns+nwarms
          do isweep=1,nsweeps
            do iu=1,ns
              r1 = u(iu)*(ran2(idummy)-0.5d0)*delu
              fprop(1:nf) = fnow(1:nf)*exp(Umat(1:nf,iu)*r1)
              normprop = sum(fprop(1:nf))
c             Form Lo.  First calculate Tf = T * f
              call dgemv('N',ndata,nf,1.0d0,T,ndata,fprop,1,0.0d0,Tf,1)
              loprop=0.5d0*sum((((normf/normprop)*Tf(1:ndata)-data(1:ndata))/error(1:ndata))**2)

              if(exp(-(loprop-lonow)/alphanow).gt.ran2(idummy)) then
c               Metropolis accept, make changes in the fields
                acrat(iu)  = acrat(iu) + 1.0d0
                lonow      = loprop
                unow(iu)   = unow(iu) + r1
                fnow(1:nf) = fprop(1:nf)
                normnow    = normprop
              end if
            end do  ! iu end one sweep of updates of field u
            
            if(irun.gt.nwarms) then ! start making measurments
              r1=normf/normnow
              fave(1:nf)  = fave(1:nf)  + r1*fnow(1:nf)
              faves(1:nf) = faves(1:nf) + (r1*fnow(1:nf))**2
              loave       = loave+lonow
            end if

          icount=icount+1
          if(icount.gt.100) then
            call dgemv('N',nf,ns,1.0d0,Umat,nf,unow,1,0.0d0,fv1,1)
            fnow(1:nf) = model(1:nf)*exp(fv1(1:nf)) ! proposed image
c           renormalize fnow and recalculate unow
            icount=0
            r1=normf/sum(fnow(1:nf))
            fnow(1:nf)=r1*fnow(1:nf)
            normnow=normf
            fprop(1:nf) = dlog(fnow(1:nf)/model(1:nf))
            call dgemv('T',nf,ns,1.0d0,Umat,nf,fprop,1,0.0d0,unow,1)

          end if
          end do  ! isweep

c         adjust delu at the end of each run so that acrat=1/2.
          ac(:)=ac(:)+acrat(:)
          acrat=acrat/dble(nsweeps)
          delu=0.5d0*delu*(sum(acrat(1:ns))/dble(ns)+1.5d0) 
          delu=min(0.5d0,delu)
c         write(6,*) 'delu=',delu,' acrat=',sum(acrat(1:ns))/dble(ns)
          acrat=0.0d0

        end do  ! irun 

        write(6,*) '  '
        write(6,*) ' i      u(i)       acrat(i) '
        do i=1,ns
           ac(i)=ac(i)/dble((nruns+nwarms)*nsweeps)
           write(6,"(1x,i2,2x,e10.4,2x,e10.4)") i,u(i),ac(i)
        end do
        write(6,*) 'L0_ave=',loave/dble(nruns*nsweeps)
        write(6,*) '  '

        rewind(88)
        rewind(89)
        write(88,*) '#     w      f_ave(w)    sdev_f_ave      f(w)     m(w)'
        do i=1,nf
           fave(i)  = max(fave(i)/(dw(i)*nruns*nsweeps),1.0d-99)
           faves(i) = max(faves(i)/(dw(i)**2*nruns*nsweeps),1.0d-99)
           fsdev(i) = max(sqrt(abs(faves(i)-(fave(i))**2)),1.0d-99)
           write(88,"(f9.4,', ',e11.4,', ',e11.4,', ',e11.4,', ',e11.4)")
     &       w(i),fave(i),fsdev(i),f(i)/dw(i),mr(i)
           write (89,"(f14.8,2x,e16.9,2x,f14.8)") w(i),max(1.0e-7,fave(i)),dw(i)           
        end do

      end do ! ialpha
      
      return
      end

      subroutine Av_spec_f
c*******************************************************************************
c
c     This subroutine implements a stochastic sampling of the spectra
c     as described by A. Sandvik.  However, it differs from Anders method
c     in that alpha (his temperature) is taken from classic MEM.  Thus,
c     this method only works with method classicJP or classic (aflag=0 or 2).  
c
c     The stochastic search is done in the image space. Note that changes 
c     are made so that Int f=1 and proposed changes are proportional to 
c     the model (or the classis MEM image f), so at high alpha (temp), 
c     f=model (or the classic MEM image f).
c
c     The image normalization is allowed to float; i.e. fnow and fprop
c     are not properly normaized.  Rather the normalizations are tracked
c     separately in normnow and normprop.  The properly normalized image
c     is then
c
c      normf                          normf
c     ------- fnow        and        -------- fprop
c     normnow                        normprop
c
c     Because the image, its normalization, product, Tf, etc are updated 
c     continuously, rather than being recalculated at each step, the code
c     is numerically unstable.  The code is stabilized by periodically 
c     recalculating fnow, normnow, Tf, prodnow, etc. from scratch.
c     
c
c*******************************************************************************
      use Global
      implicit none
      integer isweep,i,j,iu,l,irun,ialpha,icount
      real ran2
      double precision fnow(nf),fv1(nf),Tfp(ndata),
     &                 acrat,acratt,loprop,lonow,loave,r1,r2,alphanow,
     &                 normf, normprop,normnow, prodf, prodprop, prodnow


c*******************************************************************************
c     test to be sure the grid is homogeneous
      do i=2,nf
        if(dw(i).ne.dw(i-1)) then
          write(6,*) '*************** WARNING ******************'
          write(6,*) 'isamc=1 only works with a homogeneous grid'
          write(6,*) 'try isamc=2 i.e. no 1/f metric, or dw=constant'
          return
        end if
      end do

      r1=ran2(-idummy) ! see the random number generator
      write(88,*) '#  From  subroutine Av_spec_f'

c     Initiate some accumulators and variables
      fnow(1:nf)=f(1:nf); lonow=lo; icount=0
      call dgemv('N',ndata,nf,1.0d0,T,ndata,fnow,1,0.0d0,Tf,1)
      normf=sum(f(1:nf)); normprop=normf; normnow=normf
      prodf=product(f(1:nf)/dw(1:nf)); prodnow=prodf; prodprop=prodf

      do ialpha=5,1,-1 ! 
        alphanow=alpha*dble(ialpha)
        write(6,*) '*************************************************'
        write(6,*) ' alpha= ',alphanow
        write(6,*) '*************************************************'
        acrat=0.0d0; acratt=0.0d0; fave=0.0d0; faves=0.0d0; loave=0.0d0
        do irun=1,nruns+nwarms
          do isweep=1,nsweeps
            do i=1,nf
              r1=model(i)*(dble(ran2(idummy))-0.5d0)*delu ! the proposed change
              if (fnow(i)+r1.gt.0.0d0) then 
c               find proposed normalization, Tf and L0.
                normprop=normnow+r1
                Tfp(1:ndata) = Tf(1:ndata) + T(1:ndata,i)*r1
                loprop=0.5d0*sum((((normf/normprop)*Tfp(1:ndata)-data(1:ndata))/error(1:ndata))**2)
                prodprop = prodnow*((fnow(i) + r1)/fnow(i))*(normnow/normprop)**nf

                if(sqrt(prodnow/prodprop)*exp(-(loprop-lonow)/alphanow).gt.ran2(idummy)) then  
c                 Metropolis accept, make changes in the fields
                  fnow(i)     = fnow(i) + r1
                  normnow     = normprop
                  prodnow     = prodprop
                  Tf(1:ndata) = Tfp(1:ndata)
                  lonow       = loprop
                  acrat       = acrat + 1.0d0
                  icount      = icount+1
                end if
              end if  
              if(icount.gt.300) then !compensate for roundoff
                icount=0
                fnow(1:nf)=(normf/normnow)*fnow(1:nf)
                normnow=normf
                prodnow=product(fnow(1:nf)/dw(1:nf))
                call dgemv('N',ndata,nf,1.0d0,T,ndata,fnow,1,0.0d0,Tf,1)
              end if  
            end do  ! i end one sweep of updates of field f
            
            if(irun.gt.nwarms) then ! start making measurments
              r2=normf/normnow
              fave(1:nf)  = fave(1:nf)  + r2*fnow(1:nf)
              faves(1:nf) = faves(1:nf) + (r2*fnow(1:nf))**2
              loave       = loave + lonow
            end if

          end do  ! isweep

c         adjust delu at the end of each run so that acrat=0.85.
          acratt=acratt+acrat
          acrat=acrat/dble(nsweeps*nf)
          delu=0.5d0*delu*(acrat+1.5d0) 
          delu=min(2.0d0,delu)
          acrat=0.0d0

        end do  ! irun 

        loave=loave/dble(nruns*nsweeps) 
        acratt=acratt/dble((nruns+nwarms)*nsweeps*nf)
        write(6,*) '  '
        write(6,*) 'L0_ave=',loave
        write(6,*) 'acrat= ',acratt
        write(6,*) 'delu= ',delu

        rewind(88)
        rewind(89)
        write(88,*) '#     w      f_ave(w)   sdev_f_ave      f(w)       m(w)'
        do i=1,nf
           fave(i)  = max(fave(i)/(dw(i)*nruns*nsweeps),1.0d-99)
           faves(i) = max(faves(i)/(dw(i)**2*nruns*nsweeps),1.0d-99)
           fsdev(i) = max(sqrt(abs(faves(i)-(fave(i))**2)),1.0d-99)
           write(88,"(f9.4,', ',e11.4,', ',e11.4,', ',e11.4,', ',e11.4)")
     &       w(i),fave(i),fsdev(i),f(i)/dw(i),mr(i)
           write (89,"(f14.8,2x,e16.9,2x,f14.8)") w(i),max(1.0e-7,fave(i)),dw(i)           
        end do
        write(6,*) 'normave=',sum(fave(1:nf))
        write(6,*) 'normf=',sum(f(1:nf))
      end do ! ialpha
      
      return
      end


      subroutine Av_spec_f_nometric
c
c*******************************************************************************
c     This subroutine implements a stochastic sampling of the spectra
c     as described by A. Sandvik.  However, it differs from Anders method
c     in that alpha (his temperature) is taken from classic MEM.  Thus,
c     this method only works with method classicJP or classic (aflag=0 or 2).  
c
c     The stochastic search is done in the image space. Note that changes 
c     are made so that Int f=1 and proposed changes are proportional to 
c     the model (or the classis MEM image f), so at high alpha (temp), 
c     f=model (or the classic MEM image f).
c
c     The image normalization is allowed to float; i.e. fnow and fprop
c     are not properly normaized.  Rather the normalizations are tracked
c     separately in normnow and normprop.  The properly normalized image
c     is then
c
c      normf                          normf
c     ------- fnow        and        -------- fprop
c     normnow                        normprop
c
c     Because the image, its normalization, etc. are updated 
c     continuously, rather than being recalculated at each step, the code
c     is numerically unstable.  The code is stabilized by periodically 
c     recalculating fnow, normnow, etc from scratch.
c     
c
c*******************************************************************************
      use Global
      implicit none
      integer isweep,i,j,iu,l,irun,ialpha,icount
      real ran2
      double precision fnow(nf),fv1(nf),Tfp(ndata),
     &                 acrat,acratt,loprop,lonow,loave,r1,r2,alphanow,
     &                 normf, normprop,normnow

c*******************************************************************************
      r1=ran2(-idummy) ! see the random number generator
      write(88,*) '#  From  subroutine Av_spec_f_nometric'

c     Initiate some accumulators and variables
      fnow(1:nf)=f(1:nf); lonow=lo; icount=0
      call dgemv('N',ndata,nf,1.0d0,T,ndata,fnow,1,0.0d0,Tf,1)
      normf=sum(f(1:nf)); normprop=normf; normnow=normf

      do ialpha=5,1,-1 ! 
        alphanow=alpha*dble(ialpha)
        write(6,*) '*************************************************'
        write(6,*) ' alpha= ',alphanow
        write(6,*) '*************************************************'
        acrat=0.0d0; acratt=0.0d0; fave=0.0d0; faves=0.0d0; loave=0.0d0
        do irun=1,nruns+nwarms
          do isweep=1,nsweeps
            do i=1,nf
              r1=model(i)*(dble(ran2(idummy))-0.5d0)*delu ! the proposed change
              if (fnow(i)+r1.gt.0.0d0) then 
c               find proposed normalization, Tf and L0.
                normprop=normnow+r1
                Tfp(1:ndata) = Tf(1:ndata) + T(1:ndata,i)*r1
                loprop=0.5d0*sum((((normf/normprop)*Tfp(1:ndata)-data(1:ndata))/error(1:ndata))**2)

                if(exp(-(loprop-lonow)/alphanow).gt.ran2(idummy)) then  
c                 Metropolis accept, make changes in the fields
                  fnow(i)     = fnow(i) + r1
                  normnow     = normprop
                  Tf(1:ndata) = Tfp(1:ndata)
                  lonow       = loprop
                  acrat       = acrat + 1.0d0
                  icount      = icount+1
                end if
              end if  
              if(icount.gt.300) then !compensate for roundoff
                icount=0
                fnow(1:nf)=(normf/normnow)*fnow(1:nf)
                normnow=normf
                call dgemv('N',ndata,nf,1.0d0,T,ndata,fnow,1,0.0d0,Tf,1)
              end if  
            end do  ! i end one sweep of updates of field f
            
            if(irun.gt.nwarms) then ! start making measurments
              r2=normf/normnow
              fave(1:nf)  = fave(1:nf)  + r2*fnow(1:nf)
              faves(1:nf) = faves(1:nf) + (r2*fnow(1:nf))**2
              loave       = loave + lonow
            end if

          end do  ! isweep

c         adjust delu at the end of each run so that acrat=0.85.
          acratt=acratt+acrat
          acrat=acrat/dble(nsweeps*nf)
          delu=0.5d0*delu*(acrat+1.5d0) 
          delu=min(2.0d0,delu)
          acrat=0.0d0

        end do  ! irun 

        loave=loave/dble(nruns*nsweeps) 
        acratt=acratt/dble((nruns+nwarms)*nsweeps*nf)
        write(6,*) '  '
        write(6,*) 'L0_ave=',loave
        write(6,*) 'acrat= ',acratt
        write(6,*) 'delu= ',delu

        fave(1:nf)  = fave(1:nf)/(nruns*nsweeps) 
        faves(1:nf) = faves(1:nf)/(nruns*nsweeps)
        fsdev(1:nf) = sqrt(abs(faves(1:nf)-(fave(1:nf))**2))

       if(ialpha.eq.2) then   
        rewind(88)
        rewind(89)
        write(88,*) '#     w      f_ave(w)   sdev_f_ave      f(w)       m(w)'
        do i=1,nf
           write(88,"(f9.4,', ',e11.4,', ',e11.4,', ',e11.4,', ',e11.4)")
     &       w(i),fave(i)/dw(i),fsdev(i)/dw(i),f(i)/dw(i),mr(i)
           write (89,"(f14.8,2x,e16.9,2x,f14.8)") w(i),max(1.0e-7,fave(i)/dw(i)),dw(i)           
        end do
       end if

        write(6,*) 'normave=',sum(fave(1:nf))
        write(6,*) 'normf=',sum(f(1:nf))
      end do ! ialpha
      
      return
      end


      subroutine Av_spec_f_pt
c
c*******************************************************************************
c     This subroutine implements a stochastic sampling of the spectra as
c     described by A. Sandvik and K. Beach.  However, it differs from their 
c     method in that alpha (his temperature) is taken from classic MEM.  Thus,
c     this method only works with method classicJP or classic (aflag=0 or 2).  
c     The Monte Carlo used to calculate the image also uses parallel tempering
c     with parameters suggested by F. Assaad.
c
c     The stochastic search is done in the image space. Note that changes 
c     are made so that Int f=1 and proposed changes are proportional to 
c     the model (or the classis MEM image f), so at high alpha (temp), 
c     f=model (or the classic MEM image f).
c
c     The image normalization is allowed to float; i.e. fnow and fprop
c     are not properly normaized.  Rather the normalizations are tracked
c     separately in normnow and normprop.  The properly normalized image
c     is then
c
c      normf                          normf
c     ------- fnow        and        -------- fprop
c     normnow                        normprop
c
c     Because the image, its normalization, etc. are updated 
c     continuously, rather than being recalculated at each step, the code
c     is numerically unstable.  The code is stabilized by periodically 
c     recalculating fnow, normnow, etc from scratch.
c     
c*******************************************************************************
      use Global
      implicit none
      character*80 fname
      integer isweep,i,j,k,iu,l,irun,ialpha,jalpha,iswitch,istart,iend
      real ran2
      double precision fnow(nf),fv1(nf),Tfp(ndata),normf, normprop,normnow,
     &                 loprop,lonow,loave,r1,r2,r3,alphanow,alphai,alphaj

c     For parallel temporing
      integer, parameter :: nptt=8, nptskip=30
      double precision, parameter :: Rpt=1.001d0
      integer :: iptskip,icount(nptt),now(3),time1,time2
      double precision :: fpt(nf,nptt),lopt(nptt),acratpt(nptt),acrattpt(nptt),
     &                    TFpt(ndata,nptt), delupt(nptt), loj,loi,betaj,betai,
c     &                    autodate(nruns*nsweeps,nptt),autodat(nruns*nsweeps),
     &                    autotime,temp(nptt),pttries,ptaccepts(nptt)

c*******************************************************************************
      call itime(now) ! get the start time
      time1=now(1)*3600+now(2)*60+now(3)
c      write(6,*)'input data file name :'
c      read(5,"(14a)") fname
c      open(unit=33,file=fname,status='old')
c      do i=1,nf
c         read(33,*) r1,f(i),dw(i)
c         f(i)=f(i)*dw(i)
c      end do
c      write(6,*)'finished reading'
      r1=ran2(-idummy) ! see the random number generator
      write(88,*) '#  From  subroutine Av_spec_f_pt'

c *****************************************************************
c   assignning temperatures for simulation
c******************************************************************
c      do i=1,nptt
c        temp(i)=exp(dble(i-3))
c      end do
      temp(1)=0.6d0
      temp(2)=0.7d0
      temp(3)=1.0d0
      temp(4)=2.0d0
      temp(5)=5.0d0
      temp(6)=10.0d0
      temp(7)=40.0d0
      temp(8)=100.0d0

c     Initiate some accumulators and variables with the MEM results
      delupt=delu; icount=0; pttries=0.0; ptaccepts=0.0
c      call dgemv('N',ndata,nf,1.0d0,T,ndata,f,1,0.0d0,Tf,1)
c      lonow=0.5d0*sum(((Tf(1:ndata)-data(1:ndata))/error(1:ndata))**2)
      lopt=lo
      normf=sum(f(1:nf)); normprop=normf; normnow=normf
      do i=1,nptt
        Tfpt(1:ndata,i) = Tf(1:ndata) 
        fpt(1:nf,i)     = f(1:nf)
      end do
      acratpt=0.0d0; acrattpt=0.0d0; fave=0.0d0; faves=0.0d0; loave=0.0d0
       
      open(unit=21,file='accept',status='unknown')
       do irun=1,nruns+nwarms
c        acratpt=0.0d0
        do isweep=1,nsweeps   
          do ialpha=1,nptt
              alphanow=alpha*temp(ialpha)
              fnow(1:nf)=fpt(1:nf,ialpha); lonow=lopt(ialpha)
              normnow=sum(fnow(1:nf))
           do iptskip=1,nptskip ! do nptskip skips for each temperature
              do i=1,nf
                r1=model(i)*(dble(ran2(idummy))-0.5d0)*delupt(ialpha) ! the proposed change
                if (fnow(i)+r1.gt.0.d0) then 
c                 find proposed normalization, Tf and L0.
                  normprop=normnow+r1
                  Tfp(1:ndata) = Tfpt(1:ndata,ialpha) + T(1:ndata,i)*r1
                  loprop=0.5d0*sum((((normf/normprop)*Tfp(1:ndata)-
     &             data(1:ndata))/error(1:ndata))**2)
c                  if(ialpha.eq.nptt) then
c                    write(21,*)'lonow =',lonow,'loprop =',loprop,exp(-(loprop-lonow)/alphanow)
c                  end if
                  if(exp(-(loprop-lonow)/alphanow).gt.ran2(idummy)) then  
c                   Metropolis accept, make changes in the fields
                    fnow(i)              = fnow(i) + r1
                    normnow              = normprop
                    Tfpt(1:ndata,ialpha) = Tfp(1:ndata)
                    lonow                = loprop
                    acratpt(ialpha)      = acratpt(ialpha) + 1.0d0
                    icount(ialpha)       = icount(ialpha)+1
                  end if
                end if  
                if(icount(ialpha).gt.300) then !compensate for roundoff
                  icount(ialpha)=0
                  fnow(1:nf)=(normf/normnow)*fnow(1:nf)
                  normnow=normf
                  call dgemv('N',ndata,nf,1.0d0,T,ndata,fnow,1,0.0d0,Tfp,1)
                  Tfpt(:,ialpha)=Tfp(:)
                end if  
              end do  ! i end one sweep of updates of field f
           end do  ! iptskip
           fpt(1:nf,ialpha)=fnow(1:nf)
           lopt(ialpha)=lonow
              
           r2=normf/normnow
           if(irun.gt.nwarms.and.ialpha.eq.3) then ! start making measurments
             fave(1:nf)  = fave(1:nf)  + r2*fnow(1:nf)
             faves(1:nf) = faves(1:nf) + (r2*fnow(1:nf))**2
             loave       = loave + lonow
           end if
c****************************************************************************************************************           
c            autodate((irun-1)*nsweeps+isweep,ialpha)=r2*fnow(180) !choose nf=180 point for autotime calculation
c****************************************************************************************************************
          end do  ! ialpha
          
c        if(irun.eq.100) then
c          write(6,*) '**********************************'      
c         after completing nptskip sweeps for each temperature, propose a swap of 
c         configurations with the higher and lower temperatures.  All configurations 
c         are normalized.
          if(iswitch.eq.1) then
            istart=1
            iend=nptt-1
            iswitch=2
          else
            istart=2
            iend=nptt
            iswitch=1
          end if
          do ialpha=istart,iend,2 ! loop over temps and consider the next higher
c            pttries(ialpha)=pttries(ialpha)+1.0
            jalpha=ialpha+1
            if(ialpha.eq.nptt) jalpha=1    ! make a periodic boundary condition of temperature
            betai=1.d0/(alpha*temp(ialpha))  ! InverseTemp.
            betaj=1.d0/(alpha*temp(jalpha))  ! Next higher T.
            loi=lopt(ialpha)
            loj=lopt(jalpha)
            r3=exp((betaj-betai)*(loj-loi))
c            write(21,"(i8,1x,i4,1x,i4,1x,f16.8,1x,f16.8,1x,f16.8,1x,f16.8,1x,f16.8,1x,f16.8)") 
c     &        (irun-1)*nsweeps+isweep,jalpha,ialpha,temp(jalpha),temp(ialpha),loj-loi,betaj-betai,r3
c***********************************************************************************
c           R = exp(-betai*loj-betaj*loi)/exp(-betai*loi-betaj*loj)
c             = exp(betai*loi + betaj*loj - betai*loj - betaj*loi)
c             = exp((betaj-betai)*(loj-loi))
c***********************************************************************************
            if(r3.gt.ran2(idummy)) then  
c             accept the change and interchange the configurations.
              ptaccepts(ialpha)=ptaccepts(ialpha)+1.0
              Tfp(:) = Tfpt(:,ialpha)
              fnow(:) = fpt(:,ialpha)
              lonow   = lopt(ialpha)
              Tfpt(:,ialpha) = Tfpt(:,jalpha)
              fpt(:,ialpha)  = fpt(:,jalpha)
              lopt(ialpha)   = lopt(jalpha)
              Tfpt(:,jalpha) = Tfp(:)
              fpt(:,jalpha) = fnow(:)
              lopt(jalpha) = lonow
            end if        
          end do
c        end if ! irun.eq.100
        end do  ! isweep
          

        write(6,*) 'irun=',irun
c       adjust delupt at the end of each run so that acrat=0.85.
        do ialpha=1,nptt
          write(6,*)'acratpt =',acratpt(ialpha)
          acrattpt(ialpha)=acrattpt(ialpha)+acratpt(ialpha)
          acratpt(ialpha)=acratpt(ialpha)/dble(nsweeps*nf*nptskip)
          delupt(ialpha)=0.5d0*delupt(ialpha)*(acratpt(ialpha)+1.5d0) 
          delupt(ialpha)=min(2.0d0,delupt(ialpha))
          acratpt(ialpha)=0.0d0
        end do
       end do  ! irun 
        
c       do ialpha=1,nptt
c          alphanow=alpha*dble(ialpha)
c          autodat(:)=autodate(:,ialpha)
c          call autot(autodat,nruns,nsweeps,autotime)
c          write(6,*) 'autotime at alpha ',alphanow,' is ',autotime
c       end do
       loave=loave/dble(nruns*nsweeps) 
       acrattpt=acrattpt/dble((nruns+nwarms)*nsweeps*nf*nptskip)

       rewind(88)
       rewind(89)
       write(88,*) '  '
       write(88,*) 'L0_ave= ',loave
       do i=1,nptt
         write(88,*) 'ialpha =',i
         write(88,*) 'acrat  =',acrattpt(i)
         write(88,*) 'delu   =',delupt(i)
         write(88,*) 'ptacrt =',ptaccepts(i),2.d0*ptaccepts(i)/dble((nruns+nwarms)*nsweeps)
       end do  

       fave(1:nf)  = fave(1:nf)/(nruns*nsweeps) 
       faves(1:nf) = faves(1:nf)/(nruns*nsweeps)
       fsdev(1:nf) = sqrt(abs(faves(1:nf)-(fave(1:nf))**2))
         
       write(88,*) '#     w      f_ave(w)   sdev_f_ave      f(w)       m(w)'
        do i=1,nf
          write(88,"(f9.4,', ',e11.4,', ',e11.4,', ',e11.4,', ',e11.4)")
     &       w(i),fave(i)/dw(i),fsdev(i)/dw(i),f(i)/dw(i),mr(i)
          write (89,"(f14.8,2x,e16.9,2x,f14.8)") w(i),max(1.0e-7,fave(i)/dw(i)),dw(i)           
        end do
       write(6,*) 'normave=',sum(fave(1:nf))
       write(6,*) 'normf=',sum(f(1:nf))
       call itime(now) ! get the start time
       time2=now(1)*3600+now(2)*60+now(3)
       write(6,*)'time consumed is ',time2-time1

      
       return
       end


       subroutine autot(autodat,nruns,nsweeps,autotime)

c       Calculate the autocorrelation time of a measurement determined
c       by the value of iauto contained in the arrays autodat(meas*run)
c       and autosgn(meas*run).  This will be done by calculating the
c       variance of the data, vari(l) as a function of bin size.
c       the number of datum is meas*run, we will allow l, the bin 
c       size, to run from 1 to meas.
c
c                    l vari(n)                 1   N             2
c       autolth(n) = ---------     vari(n) = ---- SUM (x_m - <x>)
c                    2 vari(1)               N-1  m=1
c
c       where x_m is the bin averages.  The autocorrelation length is
c       done this way for each Markov process and then the results are
c       reduced to process 0.
c

        integer nsize,nbins,nbin,nruns,nsweeps
        double precision avdat,autodat(nsweeps*nruns),vari(nsweeps),autotime
        double precision, parameter :: zero=0.0d0        

        avdat=sum(autodat(1:nsweeps*nruns))/nsweeps*nruns
        vari(:)=zero
        do nsize=1,nsweeps         ! the bin size
          nbins=(nsweeps*nruns)/nsize  ! the number of bins of this size
          do nbin=1,nbins
c           calculate the bin averages
            r1=zero
            do i=(nbin-1)*nsize+1,nbin*nsize
              r1=r1+autodat(i)
            end do
            vari(nsize)=vari(nsize)+(r1/nsize-avdat)**2
          end do
          if (nbins.gt.1)   vari(nsize)=vari(nsize)/dble(nbins-1)
c          write(6,*) nsize,vari(nsize)
        end do
c       Now put the autocorrelation length in vari
        if (abs(vari(1)).gt.0.0001d0) then
        r1=vari(1)
        do nsize=1,min(nsweeps,40)
          vari(nsize)=0.5d0*nsize*vari(nsize)/r1
        end do
        end if 

        autotime=zero
        do nsize=10,min(40,nsweeps)
          autotime=autotime+vari(nsize)
        end do
        if(min(40,nsweeps).gt.9) 
     &     autotime=autotime/dble(min(40,nsweeps)-9)
        
        return
        end 

        FUNCTION ran2(idum)
        INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
        REAL ran2,AM,EPS,RNMX
        PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     &  IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     &  NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
        INTEGER idum2,j,k,iv(NTAB),iy
        SAVE iv,iy,idum2
        DATA idum2/123456789/, iv/NTAB*0/, iy/0/
        if (idum.le.0) then
          idum=max(-idum,1)
          idum2=idum
          do 11 j=NTAB+8,1,-1
            k=idum/IQ1
            idum=IA1*(idum-k*IQ1)-k*IR1
            if (idum.lt.0) idum=idum+IM1
            if (j.le.NTAB) iv(j)=idum
11        continue
          iy=iv(1)
        endif
        k=idum/IQ1
        idum=IA1*(idum-k*IQ1)-k*IR1
        if (idum.lt.0) idum=idum+IM1
        k=idum2/IQ2
        idum2=IA2*(idum2-k*IQ2)-k*IR2
        if (idum2.lt.0) idum2=idum2+IM2
        j=1+iy/NDIV
        iy=iv(j)-idum2
        iv(j)=idum
        if(iy.lt.1)iy=iy+IMM1
        ran2=min(AM*iy,RNMX)
        return
        END
