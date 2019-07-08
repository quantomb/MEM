        program rddcams
c*******************************************************************************
c       This code is used to process the binned dca data into a form which  
c       canbe used by the preprocessor, readmcms, of the analytic continuation
c       code.  This this code takes the binned data, and produces the
c       averages, and covariances.  The bin size is variable within
c       the code, i.e. bins can be averaged, and multiple data file
c       can be opened and averaged.  In addition it is 
c       capable of qualifying the data by performing a variety of tests.
c       These include tests to measure the kurtosis and skew of the
c       distribution, and a calculation of the probability that the
c       data is gaussian (normally distributed). This last probability
c       should probably be at least 0.005.
c
c
c*******************************************************************************
        use Global
        implicit none
c*******************************************************************************
c       Declare local variables
        character*72 linemc

        integer, parameter :: nbinsmax=8000
        integer :: i,j,k,ir,l,lr,n,nuse,ind,naver,ihisto,
     &             nbins,iopt,isym,nbinsuse,info,
     &             warm,mcskip,meas,ntr,run,niter,nuse1,
     &             nll,measl,mcskipl,nfile,nfiles,ibin,
     &             ick,icr,iequiv,ios,ic1l,ic1r,ic1,ic2,ic2r,
     &             ic2l,ic,ik,jc,idim

        real(kind) :: sus,sus1,sus2,dsus,Pmin,ave,adev,sdev,var,skew,
     &             curt,chsq,prob,error,tau,dtau,mom1,mom2,h,r1,betal,
     &             edl,Ul,tprimel,tdprimel,smeas
        
        real(kind), allocatable :: data(:,:),data2(:),aver(:),
     &                             Atau(:),Aave(:),Aerror(:),
     &                             Acurt(:),Askew(:),Aprob(:),
     &                             datar(:,:),averr(:),
     &                             Ataur(:),Aaver(:),Aerrorr(:),
     &                             Acurtr(:),Askewr(:),Aprobr(:),
     &                             store(:),ds(:,:),
     &                             gsmeas(:,:,:),data1(:,:,:)
     
        
        complex(kind) :: csum

c*******************************************************************************
c       Learn about the input files, and how the data is to be processed
c*******************************************************************************
        
c       Get the output data filename
        write(6,1) 
 1      format('               Enter output data filename:  ',$)
        read(5,201) linemc
        if(index(linemc(1:72),'bins').ne.0 )then
          write(6,*) 'This is an input file'
          stop
        end if
        open(unit=9,file=linemc,status='unknown')


c       Determine which data are to be read, depending on iopt
c       iopt=0, read the cluster single-particle Green function 
c       iopt=1, read the cluster spin susceptibility
c       iopt=2, read the cluster charge susceptibility
        write(6,2) 
 2      format('          Read Gtau (0), CHI (1) CHA (2)?:  ',$)
        read(5,*) iopt
        isym=0

c       Determine how many files are to be read.  If only one, set 
c       the number of bins to be used.  If nbinsuse> then the number
c       of bins in the file, then all of the bins are used from the
c       data file.
        write(6,3) 
 3      format('                    How many input files?:  ',$)
        read(5,*) nfiles
        if(nfiles.eq.1)then
          write(6,9) 
 9        format('            How many bins are to be used?:  ',$)
          read(5,*) nbinsuse
        else
          nbinsuse=nbinsmax
        end if

c       Pmin is used when the data too highly skewed, i.e. as when
c       there is a gap at the Fermi energy.  Pmin is the probability 
c       that the data is guassianly distributed, as determined by the 
c       Ftest (for a discussion of the statistical description of data,
c       see Numerical Recipes in Fortran 77 or C, available online at
c       http://www.nr.com).  The Ftest probability should be greater 
c       than ~10^-4, set by Pmin. If the probablity that the data for 
c       some time slice is gaussian falls below Pmin, then this time 
c       slice and the corresponding data is discarded.
c       
        write(6,4) 
 4      format('                               Enter Pmin:  ',$)
        read(5,*) Pmin

c       naver is used to re-bin the data. If naver=1, the bins read in
c       from the data files are used without re-binning.  If naver>1,
c       then the data is rebinned.  I.e. if naver=2, then data from 
c       two adjacent bins is averaged into one and the number of bins 
c       of data is halved.  If naver=3, then the data from three adjacent
c       bins is averaged and the number of bins is reduced by 1/3, etc.
c

        write(6,5) 
 5      format('Enter course grain size for bin averaging:  ',$)
        read(5,*) naver 

c       if ihisto>0 then the code attempts to generate histograms
c       of the data.  If ihisto=1, a simple fortran formatted 
c       plot is generated of each time slice data.
        ihisto=0
        write(6,6) 
 6      format('                  Print histograms (y/n)?:  ',$)
        read(5,201) linemc
        if(index(linemc(1:72),'y').ne.0 )ihisto=1
        write(6,*) '  '

        sus=0.0
        nbins=0
        ibin=0
        
c*******************************************************************************
c       Loop over the input data files, read in the data
c*******************************************************************************
        do 999 nfile=1,nfiles  
       
          write(6,7) 
 7        format('                Enter input data filename:  ',$)
          read(5,201) linemc
 201      format(a72)
          open(unit=92,file=linemc,status='old')

c         the variables nll, mcskipl, measl... etc are the variables from
c         the previous data file.  They are compared to nl, mcskip, meas...
c         to see if the different data files actually correspond.  If not,
c         the code will stop.
          nll=nl

c         READ IN THE HEADER DATA (NOTE THE ASSUMED FORM)
 400      read(92,201) linemc
          if(index(linemc(1:72),'Results').ne.0 ) then
            read(92,*) nl,nwn,Nc
            nuse=nl
           NcR=Nc/2
           ngroup=8
          else
           goto 400
          endif


          

          if(nfile.eq.1) then
c           Allocate the arrays and form the tables using the inputs from the
c           first input file.
            call allocate_1p
            allocate(data(nl+1,nbinsmax),stat=info)
            allocate(data2(nbinsmax),stat=info)
            allocate(Atau(nl+1),stat=info)
            allocate(Aave(nl+1),stat=info)
            allocate(Aerror(nl+1),stat=info)
            allocate(Acurt(nl+1),stat=info)
            allocate(Askew(nl+1),stat=info)
            allocate(Aprob(nl+1),stat=info)
            allocate(datar(nl+1,nbinsmax),stat=info)
            allocate(averr(nl+1),stat=info)
            allocate(Ataur(nl+1),stat=info)
            allocate(Aaver(nl+1),stat=info)
            allocate(Aerrorr(nl+1),stat=info)
            allocate(Acurtr(nl+1),stat=info)
            allocate(Askewr(nl+1),stat=info)
            allocate(Aprobr(nl+1),stat=info)
            allocate(store(nl+1),stat=info)
            allocate(ds(nl+1,nbinsmax),stat=info)
            allocate(gsmeas(nl,NcR,4),stat=info)
            allocate(data1(nl+1,NcR,4),stat=info)
            allocate(aver(nl+1),stat=info)
            call tables
         write(6,*) ' ick  Kcx(ick)   Kcy(ick), invicRK(ick)'
            do ick=1,Ncw
              write(6,10) ick,Kc(1,ick),Kc(2,ick),invicRK(ick)
            end do
            write(6,*) ' icr  Rcx(ick)   Rcy(ick) isublat(icr)    ',
     &        'invicRR(icr) '
            do icr=1,Nc
              write(6,11) icr,Rc(1,icr),Rc(2,icr), isublat(icr),
     &         invicRR(icr) 
            end do
 10         format(1x,i3,1x,f9.6,2x,f9.6)
 11         format(1x,i3,5x,i2,9x,i2)
            write(6,8) Nc, Ncw
 8          format('   Nc=',i4,'   Ncw=',i4,'        Enter ick:  ',$)
            read(5,*) ick
          end if
 
          mcskipl=mcskip
          measl=meas
 401      read(92,201) linemc
          if(index(linemc(1:72),'warm').ne.0 ) then
            read(92,*) warm,mcskip,meas,ntr,run,niter,nuse1
          else
           goto 401
          endif

          tprimel=tprime
          tdprimel=tdprime
          edl=ed
          Ul=U
          betal=beta
402       read(92,201) linemc
          if(index(linemc(1:72),'tdprime').ne.0 ) then
            read(92,*) ed,tprime,tdprime,U,beta
            temp=1.0/beta
          elseif(index(linemc(1:72),'tprime').ne.0 ) then
            read(92,*) ed,tprime,U,beta
            tprime=0.0
            temp=1.0/beta
          else
           goto 402
          endif
          dtau=beta/float(nl)
          if(abs(ed/U).lt.0.0001.and.abs(tprime).lt.0.001
     &      .and.abs(tdprime).lt.0.001
     &       .or.iopt.eq.1.or.iopt.eq.2)then
            write(6,*) 'Symmetric storage mode used'
            nuse=nuse/2+1
            isym=1
          end if
          if(nfile.gt.1)then
          if(tprimel.ne.tprime.or.edl.ne.ed.or.Ul.ne.U.or.
     1       mcskip.ne.mcskipl.or.tdprimel.ne.tdprime.or.
     2       betal.ne.beta.or.nll.ne.nl.or.measl.ne.meas) then
            write(6,*) 'This data set is inconsistent with prior set'
            stop
          end if
          end if

c         Now read in the binned data.  It is assumed to be stored
c         as Gs(k,tau),s.  Here, Gs is the average product of the
c         Green function and sign and s is the average sign, for a 
c         bin of data.  Note that we want the last data point to 
c         always be s.
          ds=0.0
          gsmeas=0.0
          data1=0.0
          do 15 i=1,nbinsuse
            
c           find a bin header
 408        read(92,201,end=410) linemc
            if(index(linemc(1:72),'bin of  Chi').ne.0.and.iopt.eq.1.or.
     &         index(linemc(1:72),'bin of  Cha').ne.0.and.iopt.eq.2.or.
     &         index(linemc(1:72),'bin of  Gtau').ne.0.and.iopt.eq.0 ) 
     &      then

              
c             readin a bin
              read(92,*,err=66)                     !readin all nl
     &       (((gsmeas(j,icr,ir), j=1,nl),icr=1,NcR),ir=1,4),smeas 
c             Now symmetrize the Green function data

c              gsmeas=gsmeas*smeas ! uncoment this line for old written data

              ibin=ibin+1
              data1=0.0

          do n=1,nl                                                     !0
          do ic=1,NcR    
             ic1=icicRR(ic)
             ic1r=neighbor(icicRR(ic),4)
             ic1l=neighbor(icicRR(ic),2)
            do iequiv=1,ngroup  ! impose some point-group symmetries    !0
               ic2=icrequ(ic1,iequiv)
               ic2r=icrequ(ic1r,iequiv)
               ic2l=icrequ(ic1l,iequiv)
               ic2r=neighbor(ic2r,2)
               ic2l=neighbor(ic2l,4)
               if(invicRR(ic2).eq.0)  stop'bug 1 in put_symm'           !0
               if(invicRR(ic2r).eq.0) stop'bug 2 in put_symm'           !0
               if(invicRR(ic2l).eq.0) stop'bug 3 in put_symm'           !0
                  data1(n,ic,1) = data1(n,ic,1) + 
     &              gsmeas(n,invicRR(ic2),1)
                   
                  data1(n,ic,2) = data1(n,ic,2) + 
     &              gsmeas(n,invicRR(ic2l),2)
             
                  data1(n,ic,3) = data1(n,ic,3) + 
     &              gsmeas(n,invicRR(ic2r),3)

                  data1(n,ic,4) = data1(n,ic,4) + 
     &              gsmeas(n,invicRR(ic2),4) 


            end do                                                      !0
          enddo                                                         !0
          enddo
              data1=data1/float(ngroup) 
              

           if(ick.ne.0) then                           
c               Now fourier transform gsmeas to obtain store (no   
c               factor of 1/Nc is required here)                                  

              if (invicRK(ick).ne.0) then! in the reduced BZ, kk
                 ik=icicRK(ick)
              elseif(invicRK(ick).eq.0)then ! outside the reduced BZ, k+Qk+Q
                 ik=icicRK(ickmap(ickplus(ick,icpi)))
              if(invicRK(ickmap(ickplus(ick,icpi))).eq.0) stop 'BRRRRRR'
              endif

          do n=1,nl
          csum=zeroc
          do jc=1,NcR
              ic=icicRR(jc)
            r1=zeror
            do idim=1,ndim
              r1=r1+Kc(idim,ik)*Rc(idim,ic)
            end do

          if (invicRK(ick).ne.0) then ! in the reduced BZ, kk
            csum=csum+exp(ii*r1)*(data1(n,jc,1)+
     &      data1(n,jc,2)*exp(-ii*Kc(1,ik))+ 
     &      data1(n,jc,3)*exp(ii*Kc(1,ik))+ data1(n,jc,4))

c            csum2=csum2+exp(ii*r1)*(data1(n,jc,1)-
c     &      data1(n,jc,2)*exp(-ii*Kc(1,ik))+ 
c     &      data1(n,jc,3)*exp(ii*Kc(1,ik))-data1(n,jc,4))

c            csum3=csum3+exp(ii*r1)*(data1(n,jc,1)+
c     &      data1(n,jc,2)*exp(-ii*Kc(1,ik))-  
c     &      data1(n,jc,3)*exp(ii*Kc(1,ik))-data1(n,jc,4))

            elseif(invicRK(ick).eq.0)then! outside the reduced BZ, k+Qk+Q
            csum =csum+exp(ii*r1)*(data1(n,jc,1)-
     &      data1(n,jc,2)*exp(-ii*Kc(1,ik))- 
     &      data1(n,jc,3)*exp(ii*Kc(1,ik))+data1(n,jc,4))
            else
               write(6,*) 'something wrong with invicRK(ick)'
               stop
           endif


          end do
          store(n)=halfr*real(csum)                 ! no factor of 1/Nc
 
        end do
        
                store(nl+1)=smeas  ! average sign in a bin
        else if(ick.eq.0) then      
c               calculate the DOS
                do n=1,nl
                  store(n)=halfr*(gsmeas(n,1,1)+gsmeas(n,1,4))  ! i.e. r=0
                end do
                store(nl+1)=smeas  ! average sign in a bin
              else
                write(6,*) 'Bad values for ick'
              end if                                      
             
              do j=1,nl+1
                ds(j,ibin)=store(j)
                data(j,ibin)=0.0
              end do
            else
             goto 408
            endif
 66       continue
 15       continue
 410      continue

 999    continue
        nbins=ibin



c*******************************************************************************
c       Symmetrize the data and calculate the normalization (iopt=2)
c*******************************************************************************

         if(iopt.eq.0) then
             ds(nuse+1,:)=ds(nl+1,:)
             sus=sum(ds(nuse+1,:))/float(nbins)
             ds=ds/sus
          else if(iopt.eq.1.or.iopt.eq.2) then
             sus=0.0
             r1=2.0/3.0
             do l=1,nl
                sus=sus+r1*dtau*sum(ds(l,:))
                r1=2.0-r1
             enddo
             sus=sus/float(nbins) ! average s*chi(T)
             sus1=sum(ds(nl+1,:))/float(nbins) ! average s
             sus2=sum(ds(1,:))/float(nbins)    ! average chi(tau=0)
             ds(nuse+1,:)=ds(nl+1,:)*sus/sus1  ! s*average chi(T)
             ds=ds/sus
          end if 
          nuse=nuse+1


c*******************************************************************************
c       Write some stuff out
c*******************************************************************************
        write(17,*) ick,beta,sus,sus1,sus2
        write(6,*) 'nbins=',nbins
        write(6,*) '   nl=',nl
        write(6,*) '  '
        write(6,*) 'warm,mcskip,meas,ntr,run'
        write(6,*) warm,mcskip,meas,ntr,nbins
        write(6,*) 'tprime,tdprime,ed,U,beta,sus'
        write(6,*) tprime,tdprime,ed,U,beta,sus
        nbins=nbins/naver


c*******************************************************************************
c       Now do bin averaging
c*******************************************************************************
        ind=0
        data=0.0
        do i=1,nbins
          do j=1,naver
            ind=ind+1
              do k=1,nuse
                data(k,i)=data(k,i)+ds(k,ind)/float(naver)
              end do
          end do
        end do
        
       
        

c*******************************************************************************
c       Now analyze the data (the analysis tools are described below)
c******************************************************************************  
        open(333,file='histnewms.dat')
        do l=1,nuse
          data2(:)=data(l,:)
          call moment(data2,nbins,ave,adev,sdev,var,skew,curt)
         
          tau=dtau*float(l-1)
          aver(l)=ave
          if(adev.gt.1.0e-20) then
            call histo(data2,nbinsmax,nbins,ave,var,sdev,chsq,prob,
     &                ihisto)
          end if
          error=sdev/sqrt(float(nbins))
          Atau(l)=tau
          Aerror(l)=error
          Acurt(l)=curt
          Askew(l)=skew
          Aprob(l)=prob
          write(333,*)'****************',l,nuse
        end do
         close(333)

        lr=0
        do l=1,nuse-1
          if(Aprob(l).gt.Pmin.or.l.eq.1.or.l.eq.nuse) then
c           Keep this set of data
            lr=lr+1
            datar(lr,:)=data(l,:)
            Ataur(lr)=Atau(l)
            averr(lr)=aver(l)
            Aerrorr(lr)=Aerror(l)
            Acurtr(lr)=Acurt(l)
            Askewr(lr)=Askew(l)
            Aprobr(lr)=Aprob(l)
          end if
        end do
            lr=lr+1
            datar(lr,:)=data(nuse,:)
            Ataur(lr)=Atau(nuse)
            averr(lr)=aver(nuse)
            Aerrorr(lr)=Aerror(nuse)
            Acurtr(lr)=Acurt(nuse)
            Askewr(lr)=Askew(nuse)
            Aprobr(lr)=Aprob(nuse)



        nuse=lr
        do l=1,nuse
          data(l,:)=datar(l,:)
          Atau(l)=Ataur(l)
          aver(l)=averr(l)
          Aerror(l)=Aerrorr(l)
          Acurt(l)=Acurtr(l)
          Askew(l)=Askewr(l)
          Aprob(l)=Aprobr(l)
        end do

        
        write(9,*) 'Results for nl,nuse,run,beta,Nc,ick='
        write(9,*)  nl,nuse,nbins,beta,Nc,ick
        write(9,*) ' '
        write(9,*) 'warm,mcskip,meas,ntr,run'
        write(9,*) warm,mcskip,meas,ntr,nbins
        write(9,*) ' '
        write(9,*) 'tprime,tdprime,ed,U,beta,sus'
        write(9,203) tprime,tdprime,ed,U,beta,sus
 203    format(2x,f7.3,2x,f6.4,2x,f6.4,2x,f6.3,2x,f6.3,2x,f7.3,2x,f7.4)
        write(9,*) ' '
        write(9,*) ' '
c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
        write(9,*) '   tau      Gtau(n)          +/-       curt  skew ',
     &             ' prob'

        do l=1,nuse
          write(9,202) Atau(l),aver(l),Aerror(l),Acurt(l),
     &                 Askew(l),Aprob(l)
        end do


        
        write(9,*) ' '
        write(9,*) ' Covariance data '
        write(9,*) '   i     j     Gij-GiGj        +/-'
        call covar(data,nl+1,nbinsmax,nuse,nbins,aver,9)
 202    format(1x,f7.4,1x,e14.7,1x,e14.7,1x,f5.2,1x,f5.2,1x,f7.5)



c*******************************************************************************
c       Wrap things up
c*******************************************************************************
        call deallocate_1p
        deallocate(data     )
        deallocate(data2    )
        deallocate(Atau     )
        deallocate(Aave     )
        deallocate(Aerror   )
        deallocate(Acurt    )
        deallocate(Askew    )
        deallocate(Aprob    )
        deallocate(datar    )
        deallocate(averr    )
        deallocate(Ataur    )
        deallocate(Aaver    )
        deallocate(Aerrorr  )
        deallocate(Acurtr   )
        deallocate(Askewr   )
        deallocate(Aprobr   )
        deallocate(store    )
        deallocate(ds       )
        deallocate(gsmeas    )
        deallocate(data1    )
        deallocate(aver     )
                           
        stop               
        end                


c       This series of subroutine calls is used to process
c       binned data stored in an array meas(index,binnumber).
c       
c       Subroutine moment calculated the moments of the data
c       for a particular index.  The data for the particular
c       index must be transferred to a 1-d array: 
c       
c            data(:)=meas(index,:)
c
c       Then execute
c
c       call moment(data,nbins,ave,adev,sdev,var,skew,curt)

c       where nbins is the number of bins of data.  On output
c       ave is the bin average, adev the average deviation,
c       sdev the standard deviation, var the variance, skew
c       the skew and curt the curtosis of the distribution.
c       (for a more detailed discussion of curt and skew, see
c       Numerical Recipes)

c
c
c       One can also generate a histogram of this same data, and
c       compare it to a gaussian distribution.  On input, nbinsmax
c       is the declared dimension of data, nbins the number of bins,
c       ave, var, and sdev that calculated by moments.  On output
c       prob is the probability that the binned data is normal
c       (fora discussion of chsq see Numerical Recipes).  prob
c       should be about 0.5 or bigger.
        
c       call histo(data,nbinsmax,nbins,ave,var,sdev,chsq,prob)

c       If you generate the binned data in two ways, say shuffled and
c       sequential binning, then you can compare the two data sets
c       with an FTEST.  If the binning procedure removed all the 
c       corellations between adjacent bins, then the two data set
c       are identical, and the calculated probability would be 1.

c       Use the ftest on the two data sets of G(tau)
c       write(9,*) ' '
c       write(9,*) ' FTEST DATA '
c       write(9,*) ' '
c       write(9,*) '   l      f -----> prob  '
c       do 35 i=1,nlg
c          data(:)=gsmeas1(i,:)
c          data2(:)=gsmeas2(i,:)
c          call ftest(data,run,data2,run,f,prob)
c          write(9,87) i,f,prob
c 35    continue
c 87      format(1x,i3,3x,f7.4,3x,f8.5)

c       Of course in order to analytically continue numerical data
c       the covariance matrix is necessary.  The following takes
c       the data in meas, and the associated vector of averages
c       ave(index) and calculate and writes to standard output
c       the covariance matrix.

c       Now calculate the covariance using 
c       covar(meas,nm1,n1,n2,ave,nm2)
c
c       Where 
c       meas(nm1,nm2) contains the measurements
c       nm1 declared dimension of meas
c       nm2 declared dimension of meas
c       n1 dimension of meas and ave
c       n2 dimension of meas
c       ave average over runs measurments
c
c       Thus, with partilce-hole symmetry, we could take
c       n1=nlg/2 in order to save space

c       call covar(meas,indexmax,nl,nbins,ave,binsmax)


      SUBROUTINE MOMENT(DATA,N,AVE,ADEV,SDEV,VAR,SKEW,CURT) 
      implicit none
      integer :: N, J
      real(8) DATA(N),AVE,ADEV,SDEV,VAR,SKEW,CURT,S,P
      IF(N.LE.1) then
        write(6,*) 'N must be at least 2' 
        return
      end if
      S=0. 
      DO 11 J=1,N 
        S=S+DATA(J) 
11    CONTINUE 
      AVE=S/N 
      ADEV=0. 
      VAR=0. 
      SKEW=0. 
      CURT=0. 
      DO 12 J=1,N 
        S=DATA(J)-AVE 
        ADEV=ADEV+ABS(S) 
        P=S*S 
        VAR=VAR+P 
        P=P*S 
        SKEW=SKEW+P 
        P=P*S 
        CURT=CURT+P 
12    CONTINUE 
      ADEV=ADEV/N 
      VAR=VAR/(N-1) 
      SDEV=SQRT(VAR) 
      IF(VAR.NE.0.)THEN 
        SKEW=SKEW/(N*SDEV**3) 
        CURT=CURT/(N*VAR**2)-3.   
      ELSE 
        write(6,*) 'no skew or kurtosis when zero variance' 
        return
      ENDIF 

      END SUBROUTINE MOMENT

      SUBROUTINE FTEST(DATA1,N1,DATA2,N2,F,PROB)
      implicit none
      integer :: n1, n2
      real(8) :: DATA1(N1),DATA2(N2), DF1, DF2, PROB, F, 
     &           VAR1, VAR2, AVE1, AVE2, BETAI
      CALL AVEVAR(DATA1,N1,AVE1,VAR1)
      CALL AVEVAR(DATA2,N2,AVE2,VAR2)
      IF(VAR1.GT.VAR2)THEN
         F=VAR1/VAR2
         DF1=N1-1
         DF2=N2-1
      ELSE
         F=VAR2/VAR1
         DF1=N2-1
         DF2=N1-1
      ENDIF
      PROB = BETAI(0.5*DF2,0.5*DF1,DF2/(DF2+DF1*F))
     *    +(1.-BETAI(0.5*DF1,0.5*DF2,DF1/(DF1+DF2/F)))

      END SUBROUTINE FTEST


      SUBROUTINE AVEVAR(DATA,N,AVE,VAR)
      implicit none
      integer :: N, J
      real(8) :: DATA(N), AVE, VAR, S
      AVE=0.0
      VAR=0.0
      DO 11 J=1,N
         AVE=AVE+DATA(J)
 11   CONTINUE
      AVE=AVE/N
      DO 12 J=1,N
         S=DATA(J)-AVE
         VAR=VAR+S*S
 12   CONTINUE
      VAR=VAR/(N-1)

      END SUBROUTINE AVEVAR
      
      real(8) FUNCTION BETAI(A,B,X)
      implicit none
      real(8) :: X, A, B, BT, GAMMLN, BETACF
      IF(X.LT.0..OR.X.GT.1.)then
        write(6,*) 'bad argument X in BETAI'
        return
      end if
      IF(X.EQ.0..OR.X.EQ.1.)THEN
         BT=0.
      ELSE
         BT=EXP(GAMMLN(A+B)-GAMMLN(A)-GAMMLN(B)
     *        +A*LOG(X)+B*LOG(1.-X))
      ENDIF
      IF(X.LT.(A+1.)/(A+B+2.))THEN
         BETAI=BT*BETACF(A,B,X)/A
         RETURN
      ELSE
         BETAI=1.-BT*BETACF(B,A,1.-X)/B
         RETURN
      ENDIF
      END FUNCTION BETAI
      
      
      REAL(8) FUNCTION BETACF(A,B,X)
      implicit none
      real(8), parameter :: EPS=3.0e-7
      integer, parameter :: ITMAX=100
      real(8) :: AM, BM, AZ, BZ, A, B, X, QAB, QAP, BPP, APP,
     &           D, EM, AP, QAM, AOLD, TEM, BP
      integer :: M
      AM=1.
      BM=1.
      AZ=1.
      QAB=A+B
      QAP=A+1.
      QAM=A-1.
      BZ=1.-QAB*X/QAP
      DO 11 M=1,ITMAX
         EM=M
         TEM=EM+EM
         D=EM*(B-M)*X/((QAM+TEM)*(A+TEM))
         AP=AZ+D*AM
         BP=BZ+D*BM
         D=-(A+EM)*(QAB+EM)*X/((A+TEM)*(QAP+TEM))
         APP=AP+D*AZ
         BPP=BP+D*BZ
         AOLD=AZ
         AM=AP/BPP
         BM=BP/BPP
         AZ=APP/BPP
         BZ=1.
         IF(ABS(AZ-AOLD).LT.EPS*ABS(AZ))GO TO 1
 11   CONTINUE
      write(6,*) 'A or B too big, or ITMAX too small'
 1    BETACF=AZ

      END FUNCTION BETACF
      
      

      subroutine histo(data,runmax,run,avg,var,sdev,chsq,prob,iq)
          implicit none
      integer :: run,runmax,iq, i, nx, j
      real(8) :: hist(100), curve(100), data(runmax),width,sdev,
     &           x, xbeg, xend, avg, delx, var, df, chsq, prob
      
      do i = 1, 50
        curve(i) = 0.0
        hist(i) = 0.0
      end do
      
      nx=20
      width=6.0
      xbeg=avg-0.5*width*sdev
      delx=width*sdev/float(nx)
      xend = xbeg + nx*delx
   
      do 40 i = 1, run
         x=data(i)
         if(x.lt.xbeg) goto 40
         if(xend.lt.x) goto 40
         j = (x-xbeg)/delx + 1
         hist(j) = hist(j) + 1.
 40   end do
      
      call gcurve(xbeg,delx,nx,avg,var,run,curve)
      
      call chsone(hist,curve,nx,2,df,chsq,prob)
      if(iq.eq.1)call phist(xbeg,delx,nx,hist,curve)
      
      end subroutine histo
      
      
      subroutine gcurve(xbeg,delx,nx,xm,var,ndata,curve)
          implicit none
      real(8) :: curve(*), xbeg, delx, var, sigma,
     &           x, c, xm
          integer :: nx, ndata, j
      
      sigma = sqrt(var)
      
      do j = 1, nx
         x = xbeg + (j-0.5)*delx
         c = exp(-1.0*(x-xm)**2/(2.*var))/(2.5066*sigma)
         curve(j) = c*ndata*delx
          end do
      
      end subroutine gcurve
      
      SUBROUTINE CHSONE(BINS,EBINS,NBINS,KNSTRN,DF,CHSQ,PROB)
          implicit none
          integer :: j, nbins, KNSTRN
          real(8) :: BINS(NBINS),EBINS(NBINS), df, chsq, prob, gammq
      DF=NBINS-1-KNSTRN
      CHSQ=0.
      DO 11 J=1,NBINS
         IF(EBINS(J).LE.0.) goto 11
         CHSQ=CHSQ+(BINS(J)-EBINS(J))**2/EBINS(J)
 11   CONTINUE
      PROB=GAMMQ(0.5*DF,0.5*CHSQ)

      END SUBROUTINE CHSONE

      real(8) FUNCTION GAMMQ(A,X)
          implicit none
          real(8) :: a, x, gamser, gln, gammcf
      IF(X.LT.0..OR.A.LE.0.)return
      IF(X.LT.A+1.)THEN
        CALL GSER(GAMSER,A,X,GLN)
        GAMMQ=1.-GAMSER
      ELSE
        CALL GCF(GAMMCF,A,X,GLN)
        GAMMQ=GAMMCF
      ENDIF
      END FUNCTION GAMMQ
      
      
      SUBROUTINE GCF(GAMMCF,A,X,GLN)
          implicit none
          integer, parameter :: itmax=100
          real(8), parameter :: EPS=3.E-7
          integer :: N
          real(8) :: GAMMCF, a, x, gln, gold, a0, a1, b0, b1, fac,
     &           GAMMLN, ANA, AN, ANF, G
      GLN=GAMMLN(A)
      GOLD=0.
      A0=1.
      A1=X
      B0=0.
      B1=1.
      FAC=1.
      DO N=1,ITMAX
         AN=FLOAT(N)
         ANA=AN-A
         A0=(A1+A0*ANA)*FAC
         B0=(B1+B0*ANA)*FAC
         ANF=AN*FAC
         A1=X*A0+ANF*A1
         B1=X*B0+ANF*B1
         IF(A1.NE.0.)THEN
            FAC=1./A1
            G=B1*FAC
            IF(ABS((G-GOLD)/G).LT.EPS)GO TO 1
            GOLD=G
         ENDIF
      END DO
      write(6,*) 'A too large, ITMAX too small'
 1    GAMMCF=EXP(-X+A*LOG(X)-GLN)*G

      END SUBROUTINE GCF
      
      
      SUBROUTINE GSER(GAMSER,A,X,GLN)
          implicit none
          integer, parameter :: ITMAX=100
          real(8), parameter :: EPS=3.E-7
          real(8) :: x, gamser, ap, a, sum, del, gammln, gln
          integer :: n
      GLN=GAMMLN(A)
      IF(X.LE.0.)THEN
         IF(X.LT.0.)return
         GAMSER=0.
         RETURN
      ENDIF
      AP=A
      SUM=1./A
      DEL=SUM
      DO 11 N=1,ITMAX
         AP=AP+1.
         DEL=DEL*X/AP
         SUM=SUM+DEL
         IF(ABS(DEL).LT.ABS(SUM)*EPS)GO TO 1
 11   CONTINUE
      write(6,*) 'A too large, ITMAX too small'
 1    GAMSER=SUM*EXP(-X+A*LOG(X)-GLN)

      END SUBROUTINE GSER
      
      
      real(8) FUNCTION GAMMLN(XX)
          implicit none
          integer :: j
      REAL(8) :: COF(6),STP,HALF,ONE,FPF,X,TMP,SER, XX
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *     -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
         X=X+ONE
         SER=SER+COF(J)/X
 11   CONTINUE
      GAMMLN=TMP+LOG(STP*SER)

      END FUNCTION GAMMLN

      subroutine phist(xbeg,delx,nx,hist,curve)
          implicit none
      character*1 :: blank, cross, ast
      character*50 :: zline
          integer :: nx, j, i, k
      real(8) :: hist(*), curve(*), scfact, hmax, x, xbeg, delx
      data blank, cross, ast/' ', 'X', '*'/
      
      scfact = 1.0
      hmax = 0.
      do 2 j = 1, nx
         if (hmax.gt.hist(j)) goto 2
         hmax = hist(j)
 2    continue
      
 3    if (hmax.le.50.) goto 5
      hmax = hmax/2.
      scfact = scfact/2.
      goto 3
      
 5    do 40 j = 1, nx
         x = xbeg + (j-0.5)*delx
         do 10 i = 1, 50
            zline(i:i) = blank
 10      continue
         k = hist(j)*scfact + 0.5
         if (k.le.0) go to 25
         do 20 i = 1, k
            zline(i:i) = cross
 20      continue
 25      k = curve(j)*scfact + 0.5
         if(k.gt.50) goto 40
       if(k.ne.0)  zline(k:k) = ast
         write(6,1000) x, nint(hist(j)), zline
         write(333,*) x, nint(hist(j)), hist(j)*scfact + 0.5,k
 40   continue
      
 1000 format(1h ,f7.3,i5,2h  ,(a))

      end subroutine phist

        subroutine covar(data,i1,i2,n1,n2,ave,lun)
        implicit none
        integer :: n1,n2,i,j,k,i1,i2,lun
        real(8) :: data(i1,i2),ave(i1),cij,cij2,sigg
c       thic code calculates the cavaraince of data(n1,n2) where
c       n1 is the data index and n2 is the bin index.  ave(n1)
c       contains the average values of the measurement.  The covariance,
c       and its error are printed out to logical unit 9.
        do i=1,n1
        do j=i,n1
          cij=0.0
          cij2=0.0
          do k=1,n2
            cij=cij + (data(i,k)-ave(i))*(data(j,k)-ave(j))
            cij2=cij2+((data(i,k)-ave(i))*(data(j,k)-ave(j)))**2
          end do
          cij=cij/float(n2)
          cij2=cij2/float(n2)
          sigg=(cij2-cij**2)/float(n2-1)
c         if(sigg.lt.0.0) write(6,*) ' - variance of cij ',sigg
          sigg=sqrt(abs(sigg))
          write(lun,90) i,j,cij,sigg
        end do  
        end do  

 90     format(1x,i4,2x,i4,2x,e14.7,2x,e14.7)

        end subroutine covar

