        program rddcams
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
c       This code also has two ways of treating the minus sign problem.
c
c       
c          If the average sign is LARGE, <s> > signswitch, then it is 
c       assumed that the variance of the sign is may be neglected.  In this 
c       case we set isign=1 and form the average in each bin
c
c       G = <Gs>/<s>  ! average formed in each bin
c
c       We then calculate the covariance of this set of data
c
c       C_ij = << dG_i dG_j >> ! average over the bins
c
c       where dG_i = G_i - <<G_i>> .
c
c          If the average sign is SMALL, <s> <  signswitch, then it is 
c       assumed that the variance of the sign is may not be neglected.
c       In this case we set isign=0 and must not divide by the bin average of 
c       the sign before we calculate the covariance.  It is essential that 
c       the covariance also reflect the sign fluctuations.  In this case the
c       covariance is calculated as 
c
c       C_ij = << dGs_i dGs_j >> ! average over the bins 
c
c       where dGs_i = Gs_i - <<Gs_i>> and i,j < = nl.  There is an 
c       additional row and column in the covariance matrix, due to
c       the sign.  For example,
c
c       C_ij =  << dGs_i s >>  where j = nl+1 or nuse+1
c
c
c
c*******************************************************************************
        use Global
        implicit none
c*******************************************************************************
c       Declare local variables
        character*72 linemc

        integer, parameter :: nbinsmax=8000
        integer :: i,j,k,l,lr,i1,n,nuse,ind,naver,ihisto,
     &             nbins,iopt,isym,nbinsuse,info,isign,
     &             warm,mcskip,meas,ntr,run,niter,nuse1,
     &             nll,measl,mcskipl,nfile,nfiles,ibin,ibinn,
     &             ick,icr,iequiv,ios,nproc,proc,iproc

c        real(kind), parameter :: signswitch=0.1_kind

        real(kind) :: sus,sus1,sus2,dsus,Pmin,ave,adev,sdev,var,skew,
     &                curt,chsq,prob,error,tau,dtau,mom1,mom2,h,r1,
     &                betal,edl,Ul,delta,deltal,smeas,signa,fillsbin(nbinsmax),
     &                sbin(nbinsmax),fills,signswitch

        real(kind), allocatable :: data(:,:),data2(:),aver(:),
     &                             Atau(:),Aave(:),Aerror(:),
     &                             Acurt(:),Askew(:),Aprob(:),
     &                             datar(:,:),averr(:),
     &                             Ataur(:),Aaver(:),Aerrorr(:),
     &                             Acurtr(:),Askewr(:),Aprobr(:),
     &                             ds(:,:),gsmeas(:,:),data1(:,:)


        complex(kind) :: csum

c*******************************************************************************
c       Learn about the input files, and how the data is to be processed
c*******************************************************************************

c       Get the output data filename
        write(6,"('               Enter output data filename:  ',$)")
        read(5,"(a72)") linemc
        if(index(linemc(1:72),'bins').ne.0 )then
          write(6,*) 'This is an input file'
          stop
        end if
        open(unit=9,file=linemc,status='unknown')


c*******************************************************************************
c       Determine which data are to be read, depending on iopt
c       iopt=0, read the cluster single-particle Green function 
c       iopt=1, read the cluster spin susceptibility
c       iopt=2, read the cluster charge susceptibility
c       iopt=3, read the cluster CC/B1G susceptibility
c       iopt=4, read the cluster Sc susceptibility
c       iopt=5, read the cluster B2G susceptibility
c       iopt=6, read the cluster spin susceptibility for O NMR calculation
c*******************************************************************************
        write(6,"('G(0), CHS(1), CHC(2), CC(3), SC(4),B2(5)?:  ',$)")
        read(5,*) iopt
        isym=0
c       set signswitch according to iopt (maybe should be read in??)
        if(iopt.eq.0) then
          signswitch=0.8  ! isign=0 works well for SP spectra
        else if(iopt.eq.1) then
          signswitch=0.5
        else if(iopt.eq.2) then
          signswitch=0.1  ! only use in desperation
        else
          signswitch=0.5  ! dont really know yet
        end if
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


c*******************************************************************************
c       Loop over the input data files, read in the data and learn the
c       the average sign and the number of processors used.
c*******************************************************************************
c*******************************************************************************
c*******************************************************************************
c*******************************************************************************
        sus=0.0
        nbins=0
        ds=0.0
        ibin=0
        do 999 nfile=1,nfiles  

          write(6,"('                Enter input data filename:  ',$)")
          read(5,"(a72)") linemc
          open(unit=92,file=linemc,status='old')

c*******************************************************************************
c         The variables nll, mcskipl, measl... etc are the variables from
c         the previous data file.  They are compared to nl, mcskip, meas...
c         to see if the different data files actually correspond.  If not,
c         the code will stop.
c*******************************************************************************
          nll=nl
c         READ IN THE HEADER DATA (NOTE THE ASSUMED FORM)
 400      read(92,"(a72)") linemc
          if(index(linemc(1:72),'Results').ne.0 ) then
            read(92,*) nl,nwn,cluster,ndim         ! OK
            nuse=nl
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
            allocate(ds(nl+2,nbinsmax),stat=info)
            allocate(gsmeas(nl,Nc),stat=info)
            allocate(data1(nl+1,Nc),stat=info)
            allocate(aver(nl+1),stat=info)
            call tables
            if(ndim.eq.2) then
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
            write(6,"('   Nc=',i4,'   Ncw=',i4,'        Enter ick:  ',$)") Nc, Ncw
            read(5,*) ick
          end if

c*******************************************************************************
c         Read in the first 1024 lines containing Gtau, and try to 
c         determine how many processors were used and the average sign.
c*******************************************************************************
          nproc=0; i1=0; ave=0.0
          do while (i1.lt.1024)            
          read(92,"(a72)",end=444) linemc
          if(index(linemc(1:72),'bin of  Gtau').ne.0) then
            i1=i1+1
c            read(linemc,"(35x,i20)") proc
            read(linemc(index(linemc(1:72),'process')+len('process'):72),*) proc
            if(proc.gt.nproc) nproc=proc
            read(92,*) ((gsmeas(j,icr), j=1,nl),icr=1,Nc),smeas 
            ave=ave+smeas
          end if
          end do
 444      nproc=nproc+1
          rewind(92)
          write(6,*) 'number of processors used=',nproc
          write(6,*) 'average sign', ave/float(i1)
          if(nfile.eq.1) then
            if(ave/float(i1).gt.signswitch) then
              isign=1 ! treat sign in bin averaging <G>=<Gs>/<s>
            else
              isign=0 ! treat sign as part of the covariance
            end if
          end if

c*******************************************************************************
c         Now read physical parameters from the binned data files, and check if
c         they are consistent with previous files (if any).
c*******************************************************************************

          mcskipl=mcskip
          measl=meas
 401      read(92,"(a72)") linemc
          if(index(linemc(1:72),'warm').ne.0 ) then
            read(92,*) warm,mcskip,meas,ntr,run,niter,nuse1
          else
           goto 401
          endif

          deltal=delta; edl=ed; Ul=U; betal=beta
 402      read(92,"(a72)") linemc
          if(index(linemc(1:72),'delta').ne.0 ) then
            read(92,*) ed,delta,U,beta
            temp=1.0/beta
          else
           goto 402
          endif
          dtau=beta/float(nl)
c
c         
c
          isym=1  ! the data is never symmetric
          if(nfile.gt.1)then
          if(deltal.ne.delta.or.edl.ne.ed.or.Ul.ne.U.or.
     1       mcskip.ne.mcskipl.or.
     2       betal.ne.beta.or.nll.ne.nl.or.measl.ne.meas) then
            write(6,*) 'This data set is inconsistent with prior set'
            write (6,*) 'delta,deltal',delta,deltal
            write (6,*) 'edl,ed',edl,ed
            write (6,*) 'Ul,U',Ul,U
            write (6,*) 'mcskipl,mcskip',mcskipl,mcskip
            write (6,*) 'betal,beta',betal,beta
            write (6,*) 'nll,nl',nll,nl
            write (6,*) 'measl,meas',measl,meas
            stop
          end if
          end if

c*******************************************************************************
c         Now read in the binned data.  It is assumed to be stored
c         as Gs(X,tau),s.  Here, Gs is the average product of the
c         Green function and sign and s is the average sign, for a 
c         bin of data.  Note that we want the last datum to always be s.
c*******************************************************************************
c         First, run through all the data in ifile and determine the 
c         sign and filling*sign in each bin; which are placed into 
c         ds(nl+1,ibin) and ds(nl+2,ibin), resp.
          ibinn=ibin
          do iproc=0,nproc-1
            do 16 i=1,nbinsuse             
c             find a bin header
 409          read(92,"(a72)",end=16) linemc
              if(index(linemc(1:72),'bin of  Gtau').ne.0)
     &        then
                read(linemc,"(35x,i20)") proc
                if(proc.eq.iproc) then !readin a bin
                  ibinn=ibinn+1
                  read(92,*) ((gsmeas(j,icr), j=1,nl),icr=1,Nc),smeas 
                  sbin(ibinn)=smeas
                  ds(nl+2,ibinn)=2.0*(smeas-gsmeas(1,1))  ! fills=<ns>
                end if
              else
               goto 409
              endif
 16          continue
            rewind(92)
            write(6,*) 'done first read for proc',iproc
          end do
          rewind(92)

c         Now run through all the data in ifile and read in the measurement.
          do iproc=0,nproc-1
            do 15 i=1,nbinsuse              
c             find a bin header
 408          read(92,"(a72)",end=15) linemc
              if(index(linemc(1:72),'bin of  Chi' ).ne.0.and.iopt.eq.1.or.
     &           index(linemc(1:72),'bin of  Cha' ).ne.0.and.iopt.eq.2.or.
     &           index(linemc(1:72),'bin of  CC'  ).ne.0.and.iopt.eq.3.or.
     &           index(linemc(1:72),'bin of  Sc'  ).ne.0.and.iopt.eq.4.or.
     &           index(linemc(1:72),'bin of  B2'  ).ne.0.and.iopt.eq.5.or.
     &           index(linemc(1:72),'bin of  Chi' ).ne.0.and.iopt.eq.6.or.
     &           index(linemc(1:72),'bin of  Gtau').ne.0.and.iopt.eq.0 )
     &        then
                read(linemc,"(35x,i20)") proc
                if(proc.eq.iproc) then !readin a bin
                  ibin=ibin+1
                  read(92,*) ((gsmeas(j,icr), j=1,nl),icr=1,Nc),smeas 
                  if(abs(smeas-sbin(ibin)).gt.eps) stop 'sign inconsistent'
                  if(ick.ne.0) then                           
c                   Now symmetrize the Green function data
                    data1=0.0
                    do n=1,nl                                           
                      do icr=1,Nc       
                        do iequiv=1,ngroup                               
                          data1(n,icr)=data1(n,icr)+ 
     &                                 gsmeas(n,icrequ(icr,iequiv)) 
                        end do                                  
                      end do                                              
                    end do    
                    data1=data1/float(ngroup) 
c                   Fourier transform G(X,tau) to obtain G(K,tau).      
                    do n=1,nl                                           
                      csum=(0.0,0.0)                                    
                      do icr=1,Nc                                      
                        csum=csum+ FTCoefs_R_to_K(ick,icr)*data1(n,icr)
                      end do                                            
                      ds(n,ibin)=real(csum)                               
                    end do    
                    ds(nl+1,ibin)=smeas  ! average sign in a bin
                  else if(ick.eq.0) then      
c                   calculate the R=0 result (DOS)
                    if(iopt.eq.6) then
                      ds(1:nl,ibin)=gsmeas(1:nl,icR0)+ gsmeas(1:nl,icR1)
		    else
		      ds(1:nl,ibin)=gsmeas(1:nl,icR0)  ! i.e. r=0
		    end if
                    ds(nl+1,ibin)=smeas  ! average sign in a bin
                  else
                    write(6,*) '*****ERROR********'
                    write(6,*) 'Bad values for ick'
                    write(6,*) '*****ERROR********'
                    stop
                  end if                                      
                end if
              else
                goto 408
              endif
 15         continue
            rewind(92)
            write(6,*) 'done second read for proc',iproc
          end do

 999    continue ! nfile
        nbins=ibin
        if(nbins.gt.nbinsmax) stop 'in rddca.f increase nbinsmax'
        signa=sum(ds(nl+1,1:nbins))/float(nbins)
        fills=sum(ds(nl+2,1:nbins))/float(nbins)
        write(6,*) '     total bins=',nbins
        write(6,*) 'average filling=',fills/signa
        write(6,*) '   average sign=',signa

c*******************************************************************************
c       Now do bin averaging.  This step can help to reduce correlations 
c       between the data at the expense of reducing the number of bins.
c*******************************************************************************
        nbins=nbins/naver
        ind=0; data=0.0; fillsbin=0.0; sbin=0.0
        do i=1,nbins
          do j=1,naver
            ind=ind+1
            do k=1,nl+1
              data(k,i)=data(k,i)+ds(k,ind)/float(naver)
            end do
            sbin(i)=sbin(i)+ds(nl+1,ind)/float(naver)
            fillsbin(i)=fillsbin(i)+ds(nl+2,ind)/float(naver)
          end do
        end do

c*******************************************************************************
c       Now, if iopt=2, subtract the vacuum term from the charge 
c       susceptibility.  If the sign is zero in a coarse-grained 
c       bin, then the code will crash.
c*******************************************************************************
        if(iopt.eq.2.and.ick.eq.0) then
        do i=1,nbins
          if(sbin(i).lt.eps) stop 'choose a larger CG size'
          do k=1,nl
            data(k,i)=data(k,i)-fillsbin(i)**2/sbin(i)
c            data(k,i)=data(k,i)-fillsbin(i)*fills/signa
c            data(k,i)=data(k,i)-fills**2/signa
          end do
        end do  
        else if(iopt.eq.2.and.ick.eq.ick0) then
        do i=1,nbins
          if(sbin(i).lt.eps) stop 'choose a larger CG size'
          do k=1,nl
            data(k,i)=data(k,i)-Nc*fillsbin(i)**2/sbin(i)
          end do
        end do  
        end if 
c*******************************************************************************
c       Normalize the data so that the corresponding spectrum integrates to one.
c       For iopt>0 (chi), we need to ensure that the integral of the spectra
c       does not change with temperature so that annealing can be used.
c       In this case, we can use the fact that 
c
c                /
c       chi(T) = | dw chi''(w)/w
c                /
c
c       so, if we divide the data by chi(T) it will be normalized to 1.
c
c         In this case, when the sign is not always one, we are working with
c       bins of chi(tau) times s.  Then
c
c                                   beta
c                /                  /
c       chi(T) = | dw chi''(w)/w  = | dtau <<s chi(tau) >> / <<s>>
c                /                  /
c                                   0
c
c                                 = <<s chi(T)>>/<<s>>
c       so
c
c        /    chi''(w)  w exp(-tau*w)       <<s chi(tau)>>
c        | dw -------- ----------------  =  --------------
c        /       w     1 - exp(-beta*w)         <<s>>
c
c      This may be written in the form of a normalized spectrum:
c
c        /              w exp(-tau*w)       <<s chi(tau)>>
c        | dw sigma(w) ----------------  =  --------------
c        /             1 - exp(-beta*w)      <<s chi(T)>>
c
c                       chi''(w)     <<s>>
c      where sigma(w) = --------- -------------
c                          w       <<s chi(T)>>
c
c      is positive definite and integrates to one.
c
c      When iopt=0 (single-particle G), things are a bit simpler.  Here,  of 
c      course, 
c
c       /
c       | dw A(w) dw = 1
c       /
c
c      So, if we just divide the data <Gs> by the global average sign <<s>>,
c      then it corressponds to a PAD A which is normalized to one.
c
c*******************************************************************************

        if(iopt.eq.0) then
           data(nuse+1,1:nbins)=data(nl+1,1:nbins)    ! load the signs after the last Gs
           sus=sum(data(nuse+1,1:nbins))/float(nbins) ! The global average sign <<s>>.
           data=data/sus                              ! <Gs>/<<s>>  for l <= nuse
c                                                     ! <s>/<<s>>   for l=nuse+1
        else if(iopt.ge.1) then
           sus=0.0
           r1=2.0/3.0
           do l=1,nl
              sus=sus+r1*dtau*sum(data(l,1:nbins))
              r1=2.0-r1
           enddo
           sus=sus/float(nbins)                       ! << s*chi(T) >>
           sus1=sum(data(nl+1,1:nbins))/float(nbins)  ! << s >>
c          
c          chi(T) = sus/sus1
c
           sus2=sum(data(1,1:nbins))/float(nbins)     ! << s*chi(tau=0) >>
           data(nuse+1,1:nbins)=                      ! <s> * << s*chi(T) >>/<< s >>
     &                    data(nl+1,1:nbins)*sus/sus1           
           data=data/sus                              ! < s*chi(tau) >/<< s*chi(T) >>  for l<=nuse
                                                      ! <s>/<<s>> for l=nuse+1
c
c          Now 
c              <<data(nuse+1,:)>> = <s> * << s*chi(T) >>/<< s >>
c          and
c              <<data(l,:)>> = <<s chi(tau)>> / << s*chi(T) >>   for l <= nuse
c
c          so that the integral of the spectrum corresponding to this data
c          is one.
        end if 


c*******************************************************************************
c       Write some stuff out
c*******************************************************************************
        write(17,"(2x,i4,3x,f7.3,3x,f10.4,3x,f10.4)") ick,beta,sus/sus1,sus2/sus1   ! ick, beta, chi(T), chi(tau=0)
        write(6,*) 'nbins=',nbins
        write(6,*) '   nl=',nl
        write(6,*) '  '
        write(6,*) 'warm,mcskip,meas,ntr,run'
        write(6,*) warm,mcskip,meas,ntr,nbins
        write(6,*) 'delta,ed,U,beta,sus'
        write(6,*) delta,ed,U,beta,sus


c*******************************************************************************
c       If isign=1, we dont use the sign data to calculate the covariance.
c       In this case, we must divide the data by the average sign in a bin.
c       Fortunately, both the data for l<=nuse and the sign data have been
c       divided by <<s>>, so the ratio is precisely what we want; e.g.
c
c       <Gs>/<s> when iopt=0.
c       or
c       <s chi(tau)><<s>>/(<<s chi(T)>> <s>) 
c               = <s chi(tau)>/ (<s> chi(T)) when iopt>=1.
c*******************************************************************************
        if(isign.eq.1) then
          do i=1,nbins
             if(abs(data(nuse+1,i)).lt.0.00000001) then
               write(6,*) 'sign 0 in bin=',i
               stop 
             end if
             do k=1,nl
                data(k,i)=data(k,i)/data(nuse+1,i)
             end do
          enddo
        end if

c       If isign=0, then we want to use the sign data in the covariance.
c       In this case, we dont divide by <s>.  In this case, data now contains
c       (to reiterate):
c
c       <Gs> when iopt=1.
c       or
c       <s chi(tau)> /<<s chi(T)>> 
c              = <s chi(tau)>/(<<s>> chi(T)) when iopt>=1.
c


c*******************************************************************************
c       Now analyze the data (the analysis tools are described below).
c******************************************************************************  
c       Briefly, we calculate the skew, curtosis and the f-test probability
c       of each time slice l.  
        if(ihisto.eq.1) open(333,file='histnewms.dat')

        if(isign.eq.0) nuse=nuse+1 ! if <<s>> is small use it in C_ij
        do l=1,nuse ! recall nuse now incudes the sign if isign=0
          data2(:)=data(l,:)
          call moment(data2,nbins,ave,adev,sdev,var,skew,curt)

          tau=dtau*float(l-1)
          aver(l)=ave
          if(adev.gt.1.0e-20) then
            call histo(data2,nbinsmax,nbins,ave,var,sdev,chsq,prob,ihisto)
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

c       Then we reject the data for which the probability falls below Pmin.
        lr=0
        do l=1,nuse-(1-isign)
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
        if(isign.eq.0) then ! treat the sign special
          lr=lr+1 
          datar(lr,:)=data(nuse,:)
          Ataur(lr)=Atau(nuse)
          averr(lr)=aver(nuse)
          Aerrorr(lr)=Aerror(nuse)
          Acurtr(lr)=Acurt(nuse)
          Askewr(lr)=Askew(nuse)
          Aprobr(lr)=Aprob(nuse)
        end if

c       Since we may have rejected some of the data, nuse must be reset.
c       Then we reload the data back into the original arrays.
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

c       Now write stuff out to the output data file.
c        
        write(9,*) 'Results for nl,nuse,run,beta,isign='
        write(9,*)  nl,nuse,nbins,beta,isign
        write(9,*) ' '
        write(9,*) 'warm,mcskip,meas,ntr,run'
        write(9,*) warm,mcskip,meas,ntr,nbins
        write(9,*) ' '
        write(9,*) 'delta,ed,U,beta,sus'
        write(9,"(2x,f7.3,2x,f6.4,2x,f6.3,2x,f6.3,2x,f7.3,2x,f7.4)") 
     &       delta,ed,U,beta,sus
        write(9,*) ' '
        write(9,*) ' '
c        1         2         3         4         5         6         7
c23456789012345678901234567890123456789012345678901234567890123456789012
        write(9,*) '   tau      Gtau(n)          +/-       curt  skew ',
     &             ' prob'

        do l=1,nuse
          write(9,"(1x,f7.4,1x,e14.7,1x,e14.7,1x,f5.2,1x,f5.2,1x,f7.5)") 
     &         Atau(l),aver(l),Aerror(l),Acurt(l),Askew(l),Aprob(l)
        end do



        write(9,*) ' '
        write(9,*) ' Covariance data '
        write(9,*) '   i     j     Gij-GiGj        +/-'
        call covar(data,nl+1,nbinsmax,nuse,nbins,aver,9) 



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

