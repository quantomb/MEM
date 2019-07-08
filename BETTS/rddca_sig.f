        program rddca_sig
c*******************************************************************************
c       This code is used to process the binned dca data into a form which  
c       can be used by the preprocessor, readmcms, of the analytic 
c       continuation code.  This this code takes the binned data, and 
c       produces the averages, and covariances.  The bin size is variable 
c       within the code, i.e. bins can be averaged, and multiple data files
c       can be opened and averaged.  In addition it is capable of qualifying 
c       the data by performing a variety of tests.  These include tests to 
c       measure the kurtosis and skew of the distribution, and a calculation 
c       of the probability that the data is gaussian (normally distributed). 
c       This last probability should probably be at least 0.005.  This data
c       is written to the out data file.  The data is also compared against
c       a gaussian distribution with the same average and variance. This 
c       comparison is illustrated graphically with ascii art in the file
c       histnewms.dat if you shoose to print the histograms.
c
c*******************************************************************************
        use Global
        implicit none
c*******************************************************************************
c       Declare local variables
        character*72 linemc

        integer, parameter :: nbinsmax=100000
        integer :: i,j,ic,k,l,lr,i1,n,mm,nuse,ind,naver,ihisto,nreal,
     &             nbins,iopt,nbinsuse,nbinsuse_per_proc,
     &             info,ldadata,ldbdata,
     &             warm,mcskip,meas,ntr,run,niter,nuse1,
     &             nll,measl,mcskipl,nfile,nfiles,ibin,ibinn,
     &             ick,icr,iequiv,ios,nproc,proc,iproc

        
        real(kind) :: sus,sus1,sus2,dsus,Pmin,ave,adev,sdev,var,skew,
     &                curt,chsq,prob,error,tau,dtau,mom1,mom2,h,r1,r2,
     &                betal,edl,Ul,tprimel,smeas,nmeas,signa,
     &                fills,chcrt0,sbin(nbinsmax),Epsk,EpskBar

        integer, allocatable :: jc(:)
        
        real(kind), allocatable :: data(:,:),datac(:,:),datad(:,:),data2(:),aver(:),
     &                             Aw(:),Atau(:),Ataur(:),
     &                             Aave(:),Aerror(:),
     &                             Acurt(:),Askew(:),Aprob(:),
     &                             datar(:,:),averr(:),
     &                             Awr(:),Aaver(:),Aerrorr(:),
     &                             Acurtr(:),Askewr(:),Aprobr(:),
     &                             ds(:,:),gsmeas(:,:),
     &                             iwn(:),iwnweight(:),taul(:)
           
        complex(kind) :: csum,rrcmplx,bingckf,sigbin,sigkiwn,Tmtrx,binCumu
        complex(kind), allocatable :: GcKfsc(:,:)

c*******************************************************************************
c       Learn about the input files, and how the data is to be processed
c*******************************************************************************       
c       Get the output data filename
        write(6,"('               Enter output data filename:  ',$)")
        read(5,"(a72)") linemc
        if(index(linemc(1:72),'binsw').ne.0 )then
          write(6,*) 'This is an input file'
          stop
        end if
        open(unit=9,file=linemc,status='unknown')
c*******************************************************************************
c       Determine which data are to be read, depending on iopt
c       iopt=0, read the cluster green function and calculate sig(k,tau)
c       iopt=1, read the cluster green function and calculate sig(k,iwn)
c       iopt=2, read the cluster green function and calculate Imsig(k,iwn)
c       iopt=3, read the cluster green function and calculate Cumulant M(k,iwn)
c*******************************************************************************
        write(6,"('iopt, sig(tau)(0), sig(iwn)(1), Im(sig(iwn))(2), Cumulant(3) M(iwn)=1/(iwn-ed-sig(iwn)):  ',$)")
        read(5,*) iopt

c*******************************************************************************
c       Determine how many files are to be read.  If only one, set 
c       the number of bins to be used.  If nbinsuse> run, the number
c       of bins in the file, then all of the bins are used from the
c       data file.
c*******************************************************************************
        write(6,"('                    How many input files?:  ',$)")
        read(5,*) nfiles
        if(nfiles.eq.1)then
          write(6,"('            How many bins are to be used?:  ',$)")
          read(5,*) nbinsuse
        else
          nbinsuse=nbinsmax
        end if

c*******************************************************************************
c       Pmin is used when the data too highly skewed, i.e. as when
c       there is a gap at the Fermi energy.  Pmin is the probability 
c       that the data is guassianly distributed, as determined by the 
c       Ftest (for a discussion of the statistical description of data,
c       see Numerical Recipes in Fortran 77 or C, available online at
c       http://www.nr.com).  The Ftest probability should be greater 
c       than ~10^-4, set by Pmin. If the probablity that the data for 
c       some time slice is gaussian falls below Pmin, then this time 
c       slice and the corresponding data is discarded.
c*******************************************************************************
c       
        write(6,"('                               Enter Pmin:  ',$)")
        read(5,*) Pmin

c*******************************************************************************
c       naver is used to re-bin the data. If naver=1, the bins read in
c       from the data files are used without re-binning.  If naver>1,
c       then the data is rebinned.  I.e. if naver=2, then data from 
c       two adjacent bins is averaged into one and the number of bins 
c       of data is halved.  If naver=3, then the data from three adjacent
c       bins is averaged and the number of bins is reduced by 1/3, etc.
c*******************************************************************************
c

        write(6,"('Enter course grain size for bin averaging:  ',$)")
        read(5,*) naver 

c       if ihisto>0 then the code attempts to generate histograms
c       of the data.  If ihisto=1, a simple fortran formatted 
c       plot is generated of each time slice data.
        ihisto=0
        write(6,"('                  Print histograms (y/n)?:  ',$)")
        read(5,"(a72)") linemc
        if(index(linemc(1:72),'y').ne.0 )ihisto=1
        write(6,*) '  '
      
          write(6,"('        Enter input binswxxx.dat filename:  ',$)")
          read(5,"(a72)") linemc
          open(unit=92,file=linemc,status='old')

 399      read(92,"(a72)") linemc
          if(index(linemc(1:72),'warm').ne.0 ) then
            read(92,*) warm,mcskip,meas,ntr,run,niter,nuse 
          else
           goto 399
          endif

          write(6,"('                     Enter sigma filename:  ',$)")
          read(5,"(a72)") linemc
          open(unit=73,file=linemc,status='old')

c         READ IN THE HEADER DATA FROM sigmaxxxx.dat, unit=73(NOTE THE ASSUMED FORM)
 400      read(73,"(a72)") linemc
          if(index(linemc(1:72),'Results').ne.0 ) then
            read(73,*) mm,nwn,fills,cluster,ntot,ndim,mm,chcrt0
          else
           goto 400
          endif

          read(cluster,"(i3)",iostat=ios) Nc ! should barf with characters
          if (ios.gt.0) then
            read(cluster,"(i2)",iostat=ios) Nc   ! ignore trailing charac.
            if (ios.gt.0) then
              read(cluster,"(i1)",iostat=ios) Nc ! ignore trailing charac.
              if (ios.gt.0) then
                write(6,*) "Cannot read Nc from string ",cluster
                stop                   
              endif
            endif
          endif

c*******************************************************************************
c         Now read physical parameters from the sigmaxxxx.dat file
c*******************************************************************************
 402      read(73,"(a72)") linemc
          if(index(linemc(1:72),'tprime').ne.0 ) then
             if(index(linemc(1:72),'tdprime').ne.0 ) then
                if(index(linemc(1:72),'t3prime').ne.0 ) then
                   read(73,*) ed,tprime,tdprime,t3prime,U,beta
                else
                   read(73,*) ed,tprime,tdprime,U,beta
                   t3prime=0.d0
                endif 
             else
                read(73,*) ed,U,beta,mm,mm,mm,tprime,r1
                tdprime=0.d0
                t3prime=0.d0
             endif
            temp=1.0/beta
          else
           goto 402
          endif

 403      read(73,"(a72)") linemc
          if(index(linemc(1:72),'nuse').ne.0 ) then
            read(73,*) mm,nuse 
          else
           goto 403
          endif
          nl=nuse
         
c           Allocate the arrays and form the tables using the inputs from the
c           first input file.
            call allocate_1p
            ldadata=2*nuse+2
            ldbdata=nbinsmax
            allocate(data(ldadata,ldbdata),stat=info)            
            allocate(datac(ldadata,ldbdata),stat=info)
            allocate(datad(ldadata,ldbdata),stat=info)
            allocate(datar(ldadata,ldbdata),stat=info)
            allocate(ds(ldadata,ldbdata),stat=info)
            allocate(data2(nbinsmax),stat=info)
            allocate(gsmeas(2*nuse,Nc),stat=info)

            allocate(Aw(2*nuse+1),stat=info)
            allocate(Atau(2*nuse+1),stat=info)
            allocate(Ataur(2*nuse+1),stat=info)
            allocate(Aave(2*nuse+1),stat=info)
            allocate(Aerror(2*nuse+1),stat=info)
            allocate(Acurt(2*nuse+1),stat=info)
            allocate(Askew(2*nuse+1),stat=info)
            allocate(Aprob(2*nuse+1),stat=info)

            allocate(averr(2*nuse+1),stat=info)
            allocate(Awr(2*nuse+1),stat=info)
            allocate(Aaver(2*nuse+1),stat=info)
            allocate(Aerrorr(2*nuse+1),stat=info)
            allocate(Acurtr(2*nuse+1),stat=info)
            allocate(Askewr(2*nuse+1),stat=info)
            allocate(Aprobr(2*nuse+1),stat=info)
            allocate(aver(2*nuse+1),stat=info)
 
            call tables
            allocate(jc(ngroup),stat=info)

            allocate(GcKfsc(-nwn:nwn,Nc),stat=info)
            allocate(iwn(nwn),stat=info)
            allocate(iwnweight(nwn),stat=info)
            allocate(taul(nuse),stat=info)


        sus=0.0
        nbins=0
        ds=0.0
        ibin=0

c*******************************************************************************
c    Now read GcKksc from sigmaxxx.dat
c    in put.f, the GcKfsc is written out as following:
c*******************************************************************************    
c            write(lus,*) 'wn   ic  sigma(n,ic)  sigmapt(n,ic)  GcKfsc(n,ic)  iwn(n)   quadw(n)'
c            do ic=1,Ncw         ! Only write for K in the IW    '     
c            do n=0,nwn                                             
c              z = sigma(n,ic)      !convert to kindg
c              z1= sigmapt(n,ic)     !convert to kindg
c              z2= GcKfsc(n,ic)
c              r7= quadw(n)
c              write(lus,44) n,ic,z,z1,z2,real(2*n+1)*pi/beta,r7 !'
c            end do
c            end do
c44       format(1x,i4,1x,i4,2x,'(',e13.6,',',e13.6,')',x,
c     &     '(',e13.6,',',e13.6,')',x,'(',e13.6,',',e13.6,')',x,f13.6,x,f13.6)                                                                       

c****************************************************************************** 
            write(6,*) 'start read GcKfsc'
 2104       read(73,"(a72)") linemc                                     
            if(index(linemc(1:72),'GcKfsc(').ne.0 ) then 

              do i=1,Ncw              !read data in the IW            
              do n=0,nwn-1                                             
c  only positive frequency is read in. The negative frequency of GcKfsc(-iwn) is just the conjg(GcKfsc(iwn))
                read(73,*,end=1997) mm,mm,sigkiwn,rrcmplx,GcKfsc(n,i),iwn(n+1),iwnweight(n+1)           
                if(n.lt.nwn) GcKfsc(-n-1,i)=conjg(GcKfsc(n,i))          
c  output \sig(K,iwn) from sigmaxxxx.dat to check
c                write(445,*) iwn(n+1), (aimag(sigkiwn))/(chcrt0*U*U),
c     &             (real(sigkiwn)-U*(fills/2.0d0-0.5d0))/(chcrt0*U*U)
c  output \Sig(K,iwn) from sigmaxxxx.dat to check
c                write(556,*) iwn(n+1), aimag(sigkiwn), real(sigkiwn)
c  output global GcKf(K,iwn) from sigmaxxxx.dat to check
c                write(778,*) iwn(n+1), aimag( 1.0d0/(1.0d0/GcKfsc(n,i)-sigkiwn) ), 
c     &           real( 1.0d0/(1.0d0/GcKfsc(n,i)-sigkiwn) )


              end do   
                read(73,*,end=1997) mm,mm
                if( abs(iwn(1)-pi/beta).gt.eps  ) then
                  write(6,*) 'wrong frequency wn, pay attention to extra line '
                end if                                       
              end do   
              if(abs(iwn(1)-pi*temp).gt.eps) then
                write(6,*) '*********************'
                write(6,*) 'read error sigma.dat '
                write(6,*) '*********************'
              end if
            else                                                      
              goto 2104                                               
            endif
1997        continue 
c           Now fill in data outside the IW 
            do i=Ncw+1,Nc                 
            do n=0,nwn-1                                            
              GcKfsc(n,i)=GcKfsc(n,ickmap(i))                         
              GcKfsc(-n-1,i)=conjg(GcKfsc(n,i))                       
            end do                                                  
            end do          
                                        

c*******************************************************************************
c    Print out K information
c    Get ick for the K point we want
c*******************************************************************************
            if(ndim.eq.1) then
              write(6,*) ' ick  Kcx(ick)'
              do i=1,Ncw
                write(6,"(1x,i3,1x,f9.6)") i,Kc(1,i)
              end do
              write(6,*) ' icr  Rcx(icr)'
              do icr=1,Nc
                write(6,"(1x,i3,5x,i2)") icr,Rc(1,icr)
              end do
            else if(ndim.eq.2) then
              write(6,*) ' ick  Kcx(ick)   Kcy(ick)'
              do i=1,Ncw
                write(6,"(1x,i3,1x,f9.6,2x,f9.6)") i,Kc(1,i),Kc(2,i)
              end do
              write(6,*) ' icr  Rcx(icr)   Rcy(icr)'
              do icr=1,Nc
                write(6,"(1x,i3,5x,i2,9x,i2)") icr,Rc(1,icr),Rc(2,icr)
              end do
            else if(ndim.eq.3) then
              write(6,*) ' ick  Kcx(ick)   Kcy(ick)   Kcz(ick)   '
              do i=1,Ncw
                write(6,"(1x,i3,1x,f9.6,2x,f9.6,2x,f9.6)") i,Kc(1,i),Kc(2,i),Kc(3,i)
              end do
              write(6,*) ' icr  Rcx(icr)   Rcy(icr)   Rcz(icr)'
              do icr=1,Nc
                write(6,"(1x,i3,5x,i2,9x,i2,9x,i2)") icr,Rc(1,icr),Rc(2,icr),Rc(3,icr)
              end do
            end if
            write(6,"('      Nc=',i4,' Ncw=',i4,'  Enter ick:  ',$)") Nc, Ncw
            read(5,*) ick


c*******************************************************************************
c     in DCA sumup.f, the G(iwn,K) is written out as following
c*******************************************************************************
c        lubn = 20000+myrank
cc           The bins of Gf  are determined printed onto the data file.
c            write(lubn,*) ' '
c              write(lubn,*) ' '
c              write(lubn,*) '   '
c              if(meastype.eq.1) 
c     &          write(lubn,*) ' bin of  Gw(n) from process  ',myrank
c        write(lubn,402) 
c     &         ((real(measurement(j,ic)), aimag(measurement(j,ic)),
c     &                        j=0,nuse-1),ic=1,Nc),sgnm,n_s
c 402    format(1x,f11.8,1x,f11.8,1x,f11.8,1x,f11.8,1x,f11.8,1x,f11.8)
c*******************************************************************************
c     for homogenous grid, nuse is set as following:
c          if(inhgrid.eq.0) then
c            nuse = max(Int(0.25d0*beta*max(U,1.0d0*real(ndim))),8)
c            if(nuse.ge.nwn) nuse=nwn
c          endif
c****************************************************************************** 
c         Now run through all the data in ifile and read in the measurement of G(iwn,K), <s>, <ns>
           iproc=0
           do i=1,nbinsuse             
c             find a bin header
 408          read(92,"(a72)",end=15) linemc
              if(index(linemc(1:72),'bin of  Gw(').ne.0) then
                read(linemc,"(35x,i20)") proc
                if(iproc.ne.proc) write(6,*)'done proc',iproc
                iproc=proc
c    Here the QMC data start from matsubara frequency wn(0)=pi*temp, and there 
c    are 2*nuse*Nc+1 numbers in one bin
!    PAY ATTENTION -- HERE 2*nuse is used!!! 
                  read(92,*,err=408) ((gsmeas(j,icr), j=1,2*nuse),icr=1,Nc),smeas,nmeas 
                  ibin=ibin+1
                  if(ick.gt.0) then                           
c                   Now symmetrize the mastubara frequency  Greens function G(k,iwn)
                    do n=1,2*nuse                                       
                     do j=1,ngroup
                      ds(n,ibin)=ds(n,ibin) + gsmeas(n,ickequ(ick,j))/float(ngroup)   
                     end do 
                    end do    
                    ds(2*nuse+1,ibin)=smeas  ! average sign in a bin 
                    ds(2*nuse+2,ibin)=nmeas  ! average one spin n times sign in a bin, nfoccm*sign
                  else
                    write(6,*) '*****ERROR********'
                    write(6,*) 'Bad values for ick'
                    write(6,*) '*****ERROR********'
                    stop
                  end if !ick can not be 0 for MemSig                                   
              else
                goto 408
              endif !index
 15        end do 

        nbins=ibin
        if(nbins.gt.nbinsmax) stop 'in rddca_sig.f increase nbinsmax'

        signa=sum(ds(2*nuse+1,1:nbins))/float(nbins)
        fills=sum(ds(2*nuse+2,1:nbins))/(float(nbins)*signa)!This is filling/2
        write(6,*) '     total bins=',nbins
        write(6,*) 'average filling=',fills
        write(6,*) '   average sign=',signa


c*******************************************************************************
c       Now do bin averaging.  This step can help to reduce correlations 
c       between the data at the expense of reducing the number of bins.
c*******************************************************************************
        nbins=nbins/naver
        ind=0; data=0.0; sbin=0.0
        do i=1,nbins
        do j=1,naver
          ind=ind+1
          do k=1,2*nuse+2
            data(k,i)=data(k,i)+ds(k,ind)/float(naver)
          end do
        end do   
        end do

c*****************************************************************
c   recover <gwmeas*sgn>=<(T matrix)*sgn> from GcKf(ick,iwn,ibins) store in Tmtrx Chen and Shiquan believe the sgn
c   in the binsw data actually doesn't correspond to GcKf(ick,iwn,ibins), but to T(ick,iwn,ibins)
c   see from the equation: GcKf=GcKfsc-GcKfsc^2*T , only the T has multiple copies from QMC 
c   associated with sgn, GcKfsc is from DCA myrank=0, only one copy, no corresponding sgn.
c*****************************************************************
c   data(iwn,ibins)=GcKf(ick,iwn,ibins)
        do i=1,nbins
          if(abs(data(2*nuse+1,i)).lt.0.00000001) then
            write(6,*) 'sign 0 in bin=',i
            stop 
          end if
          do j=1,nuse
            Tmtrx=cmplx( data(2*j-1,i),data(2*j,i) ) 
            datad(2*j-1,i)=real( (-Tmtrx+GcKfsc(j-1,ick))/GcKfsc(j-1,ick)**2 )
            datad(2*j ,i)=aimag( (-Tmtrx+GcKfsc(j-1,ick))/GcKfsc(j-1,ick)**2 )
c*******************************************************************************
c      take care of the sign by dividing  <Ts>/<s> Ts from datad(1:2*nuse,:) <s> from data(2*nuse+1,:)
c*******************************************************************************
            do k=1,2*nuse
              datad(k,i)=datad(k,i)/data(2*nuse+1,i)  ! <Ts>/<s>.
            end do
c       calculate GcKf again from sign free T matrix in datad
            Tmtrx=cmplx( datad(2*j-1,i),datad(2*j,i) ) 
            data(2*j-1,i)=real( (-GcKfsc(j-1,ick)**2*Tmtrx)+GcKfsc(j-1,ick) )
            data(2*j ,i)=aimag( (-GcKfsc(j-1,ick)**2*Tmtrx)+GcKfsc(j-1,ick) )                      
          enddo
        enddo 

c*******************************************************************************
c      take care of the sign of each bin average fills by dividing  <ns>/<s>
c*******************************************************************************
          do i=1,nbins
c             if(abs(data(2*nuse+1,i)).lt.0.00000001) then
c               write(6,*) 'sign 0 in bin=',i
c               stop 
c             end if
c             do k=1,2*nuse
c                data(k,i)=data(k,i)/data(2*nuse+1,i)  ! <Gs>/<s>.
c             end do
             data(2*nuse+2,i)=data(2*nuse+2,i)/data(2*nuse+1,i)  ! <ns>/<s>.
          enddo

c*******************************************************************************
c     output global GcKf to check
c*******************************************************************************
       do j=1,nuse
         r1=sum(data(2*j,1:nbins))/real(nbins)    ! Imaginary part  GcKf
         r2=sum(data(2*j-1,1:nbins))/real(nbins)  ! Real part       GcKf
c         write(667,*) iwn(j),r1,r2
       enddo

c********************************************************
c    From here to below, the GcKf(2*nuse,ick) in data(2*nuse,ibin) and nfocc
c    in data(2*nuse+2,ibin) is sign problem free
c********************************************************


c*******************************************************************************
c       Apply Dyson equation on data(2*nuse+2,nbins)=GcKf
c        data(2*nuse,nbins)=\Sig= 1/GcKfsc(2*nuse,ick)-1/data(2*nuse,nbins)
c       after that data(2*nuse+2,nbins)=\Sig
c       nbins*naver= how many raw sets of GcKf in ds
c       only calculate positive frequency \Sigma within low frequency from 0 to nuse-1
c*******************************************************************************
        if(iopt.eq.3) then ! Prepare Cumulant <M> = 1/(iwn-ed-<sig>)


        do i=1,nbins
        do j=1,nuse   ! data(2*nuse,:) stored the sign of the binned GcKf, unchanges here
          bingckf=cmplx( data(2*j-1,i),data(2*j,i) ) 
          sigbin=1.0d0/GcKfsc(j-1,ick)-1.0d0/bingckf
          binCumu=1.0d0/(iwn(j)-ed-sigbin)
          data(2*j-1,i)=real(binCumu)
          data(2*j,i)=aimag(binCumu)
        enddo
        enddo


        else
        do i=1,nbins
        do j=1,nuse   ! data(2*nuse,:) stored the sign of the binned GcKf, unchanges here
          bingckf=cmplx( data(2*j-1,i),data(2*j,i) ) 
          sigbin=1.0d0/GcKfsc(j-1,ick)-1.0d0/bingckf
          data(2*j-1,i)=real(sigbin)
          data(2*j,i)=aimag(sigbin)
        enddo
        enddo
        end if!iopt


c*******************************************************************************
c     output \Sig(iwn,K) to check
c*******************************************************************************
       do j=1,nuse
         r1=sum(data(2*j,1:nbins))/real(nbins)    ! Imaginary part  Sig
         r2=sum(data(2*j-1,1:nbins))/real(nbins)  ! Real part       Sig
c         write(223,*) iwn(j),r1,r2
       enddo

c*******************************************************************************
c       subtract Hatree term
c*******************************************************************************
        do i=1,nbins
        do j=1,nuse   ! data(2*nuse,:) stored the sign of the binned GcKf, unchanges here
          data(2*j-1,i)=data(2*j-1,i)-U*(fills-0.5d0)
c          data(2*j-1,i)=data(2*j-1,i)-U*(data(2*nuse+2,i) - 0.5d0)
        enddo
        enddo

c*******************************************************************************
c      divide U^2*\chi_c(r=0,\tau=0) (old version, use global \chi_c(r=0,\tau=0)
c       is not a good approximation for QCP region
c*******************************************************************************
c       data(1:2*nuse,:)=data(1:2*nuse,:)/(chcrt0*U*U)

c*******************************************************************************
c      divide U^2*nfocc_(averge just in each bin)*(1-nfocc), nfocc is in data(2*nuse+2,:)
c*******************************************************************************
        do i = 1, nbins
          data(1:2*nuse,i)=data(1:2*nuse,i)/(U*U*(1.0d0-fills)*fills)
c          data(1:2*nuse,i)=data(1:2*nuse,i)/(U*U*(1.0d0-data(2*nuse+2,i))*data(2*nuse+2,i))
        enddo

c*******************************************************************************
c     output \sig(iwn,K) to check
c*******************************************************************************
       do j=1,nuse
         r1=sum(data(2*j,1:nbins))/real(nbins)    ! Imaginary part
         r2=sum(data(2*j-1,1:nbins))/real(nbins)  ! Real part
c         write(334,*) iwn(j),r1,r2
       enddo

       if(iopt.eq.0) then
c*******************************************************************************
c        form imaginary time grid taul for n*\delta\tau
c*******************************************************************************   
          do i=1,nuse
            taul(i)=real(i-1)*beta/real(nuse)
          enddo


c*******************************************************************************
c         Fourier transform \sig(iwn,K) to \sig(\tau,K)
c*******************************************************************************
c         similar to the subroutine FT_Kf_to_Kt(GcgKf,GcgKt,itype) in geod.f
          datac=0.0d0
          do i=1,nbins
          do j=1,nuse   !j labels tau in Sig    
             do k=1,nuse
               sigbin=cmplx( data(2*k-1,i), data(2*k,i) )           
               datac(j,i)= datac(j,i)+ iwnweight(k)*
     &                  real( exp(-ii*iwn(k)*taul(j))* (sigbin-1.0d0/(ii*iwn(k)) ) )
             enddo
c  put in a minus sign by hand to consistent with default model
             datac(j,i)=-(datac(j,i)*2.0d0/beta-0.5d0)
          enddo
          enddo
       else if (iopt.eq.1.or.iopt.eq.2) then
         datac=data  ! AC sigma(k,iwn)
c        form matsubara frequency grid taul
          do i=1,nuse
            taul(i)=iwn(i)
          enddo
       else if (iopt.eq.3) then
         datac=data  ! AC Cumulant M(k,iwn)
c        form matsubara frequency grid taul
          do i=1,nuse
            taul(i)=iwn(i)
          enddo
       else

         write(6,*) '*****ERROR********'
         write(6,*) 'Bad value for iopt'
         write(6,*) '*****ERROR********'
         stop
       end if


c*******************************************************************************
c       Write some stuff out
c*******************************************************************************

c        write(6,"(2x,i4,3x,f7.3,3x,f10.4,3x,f10.4)") ick,beta,sus/sus1,sus2/sus1   ! ick, beta, chi(T), chi(tau=0)

        write(6,*) 'nbins=',nbins
        write(6,*) '   nl=',nl
        write(6,*) '  '
        write(6,*) 'warm,mcskip,meas,ntr,run'
        write(6,*) warm,mcskip,meas,ntr,nbins
        write(6,*) 'tprime,ed,U,beta,sus'
        write(6,*) tprime,ed,U,beta,sus
        

c*******************************************************************************
c       Now analyze the data (the analysis tools are described below).
c******************************************************************************  
c       Briefly, we calculate the skew, curtosis and the f-test probability
c       of each time slice l.

c*************************************
c       use iopt to adjust nuse value.  
c       iopt=0 for \sig(K,\tau), nuse=nuse (real)
c       iopt=1 for \sig(K,iwn),nuse=2*nuse (complex)
c       iopt=2 for \Im sig(K,iwn),nuse=nuse (real)
c*************************************
        if (iopt .eq. 0.or.iopt.eq.2) then
          nuse=nuse
        else if (iopt .eq. 1.or.iopt.eq.3) then
          nuse=2*nuse
        else
          write(*,*)'wrong iopt'
          stop
        endif
  
        if(ihisto.eq.1) open(333,file='histnewms.dat')

        do l=1,nuse   ! upper limit is nuse if iopt=0 or 2, is 2*nuse if iopt=1 or 3
          data2(:)=datac(l,:)
          call moment(data2,nbins,ave,adev,sdev,var,skew,curt)

c  tau=(-iopt)**(l-1)*taul(iopt+(l-iopt)/(iopt+1)) ! 12-1 34-2          
          if (iopt .eq. 0) then ! set Atau(l)=\tau or iwn
            tau=taul(l)
          else if (iopt .eq. 1.or.iopt.eq.2.or.iopt.eq.3) then ! set Atau(l)=iwn 1,2->1, 3,4->2, 5,6->3          
            tau=(-1)**(l-1)*taul(1+(l-1)/2)
          endif

          write(*,*)'line540 l,tau',l,tau

          aver(l)=ave
          prob=0.0
          if(adev.gt.1.0e-20) then
            call histo(data2,nbinsmax,nbins,ave,var,sdev,chsq,prob,ihisto)
          end if
          error=sdev/sqrt(float(nbins))
          Atau(l)=tau
          Aerror(l)=error
          Acurt(l)=curt
          Askew(l)=skew
          Aprob(l)=prob
          if(ihisto.eq.1) write(333,*)'****************',l,nuse
        end do
        if(ihisto.eq.1) close(333)

c       Then we reject the data for which the probability falls below Pmin.
        lr=0
        do l=1,nuse

           if(Aprob(l).gt.Pmin.or.l.eq.1.or.l.eq.nuse.and.iopt.ne.7) then  !!origin
c           Keep this set of data
            lr=lr+1
            datar(lr,:)=datac(l,:)
            Ataur(lr)=Atau(l)
            averr(lr)=aver(l)
            Aerrorr(lr)=Aerror(l)
            Acurtr(lr)=Acurt(l)
            Askewr(lr)=Askew(l)
            Aprobr(lr)=Aprob(l)
          end if
        end do

c       Since we may have rejected some of the data, nuse must be adjusted again.
c       Then we reload the data back into the original arrays.
        nuse=lr
        do l=1,nuse
          datac(l,:)=datar(l,:)
          Atau(l)=Ataur(l)
          aver(l)=averr(l)
          Aerror(l)=Aerrorr(l)
          Acurt(l)=Acurtr(l)
          Askew(l)=Askewr(l)
          Aprob(l)=Aprobr(l)
        end do

c       Now write stuff out to the output data file.
c        
        write(9,*) 'Results for nl,nuse,run,beta,isign='
        write(9,*)  nl,nuse,nbins,beta,1 
c   isgin is not relevant in this code, output to be 1 for consistency for the following code.
        write(9,*) ' '
        write(9,*) 'warm,mcskip,meas,ntr,run'
        write(9,*) warm,mcskip,meas,ntr,nbins
        write(9,*) ' '
        write(9,*) 'tprime,ed,U,beta,sus'
        write(9,"(2x,f7.3,2x,f7.4,2x,f6.3,2x,f6.3,2x,f7.3,2x,f7.4)") 
     &       tprime,ed,U,beta,sus
        write(9,*) ' '
        write(9,*) ' '
c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
        write(9,*) '   tau      Gtau(n)          +/-       curt  skew ',
     &             ' prob'
        do l=1,nuse
          write(9,"(1x,f9.4,1x,e14.7,1x,e14.7,1x,f,1x,f,1x,f7.4)") 
     &         Atau(l),aver(l),Aerror(l),Acurt(l),Askew(l),Aprob(l)
        end do


        
        write(9,*) ' '
        write(9,*) ' Covariance data '
        write(9,*) '   i     j     Gij-GiGj        +/- '
        call covar(datac,ldadata,ldbdata,nuse,nbins,aver,9) 



c*******************************************************************************
c       Wrap things up
c*******************************************************************************
        call deallocate_1p
        deallocate(data     )

        deallocate(data2    )
        deallocate(Aw       )
        deallocate(Atau     )
        deallocate(Ataur    )
        deallocate(Aave     )
        deallocate(Aerror   )
        deallocate(Acurt    )
        deallocate(Askew    )
        deallocate(Aprob    )
        deallocate(datar    )
        deallocate(averr    )
        deallocate(Awr      )
        deallocate(Aaver    )
        deallocate(Aerrorr  )
        deallocate(Acurtr   )
        deallocate(Askewr   )
        deallocate(Aprobr   )
        deallocate(ds       )
        deallocate(gsmeas   )
        deallocate(aver     )
        deallocate(jc       )

        deallocate(GcKfsc   )
        deallocate(iwn      )
        deallocate(iwnweight)
        deallocate(taul     )
        deallocate(datac    )
        deallocate(datad    )
                           
        stop              
        end program rddca_sig              



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
c       thic code calculates the covariance of data(n1,n2) where
c       n1 is the data index and n2 is the bin index.  ave(n1)
c       contains the average values of the measurement.  The covariance,
c       and its error are printed out to logical unit 9.
c
c       data    real array of size data(n1,n2) containing the data.
c       i1      the leading index of data as declared in the calling routine
c       i2      the second index of data as declared in the calling routine
c       n1      the data index (i.e. the time slice index l)
c       n2      the number of bins of data
c       ave     on return, array containing the data averages
c       lun     the logical unit to which the covariance data is written
c
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
cccccccccccccccccccccccccccccccccccccccccccccc
        real*8 function Epsk(kx,ky)
        use Global
        implicit none
        real(kind) :: kx, ky    
        
  


         Epsk=-0.5_kind*(cos(kx)+cos(ky)) -
     &       tprime*(cos(kx)*cos(ky)-1.0_kind)
     &             - 0.5*tdprime*(cos(2.d0*kx)+cos(2.d0*ky))
     &             - t3prime*(cos(2.d0*kx)*cos(ky)
     &                         +cos(kx)*cos(2.d0*ky))   


        
        end function Epsk
 

