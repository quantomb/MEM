        program giwbssigiw
c*******************************************************************************
c       This code is used to bootstrap (bs for short) the dca G(iw_n) into Sigma(iw_n) bins data
c       Then Sigma(iw_n) bins will be used in rddcams to obtain covariance matrix.

c*******************************************************************************
        use Global
        implicit none
c*******************************************************************************
c       Declare local variables
        character*72 linemc

        integer, parameter :: nbinsmax=16000
        integer :: i,j,ic,jc,k,l,lr,i1,n,nuse,ind,naver,
     &             nbins,nbinsuse,nbinsuse_per_proc,nbinsbs,
     &             info,isign,iopt,
     &             ranclmn,iseed,
     &             warm,mcskip,meas,ntr,run,niter,nuse1,
     &             nll,measl,mcskipl,nfile,nfiles,ibin,ibinn,
     &             ick,icr,iequiv,ios,nproc,proc,iproc,
     &             numbermeas,measperproc, numbermeasprev,indexofmeas
        integer, allocatable :: indexinproc(:)
        
        real(kind) :: betal,edl,Ul,tprimel,smeas,signa,ave,
     &                fillsbin(nbinsmax),rang
        
        real(kind), allocatable :: rawtobs(:,:),dsbsob(:,:,:),dsbssgn(:),
     &                             G_bs(:,:),Sig_bs(:,:),
     &                             gsmeas(:,:),GcKfsc(:,:)
     
        
        complex(kind) :: csum

c*******************************************************************************
c       Learn about the input files, and how the data is to be processed
c*******************************************************************************
        
c       Get the output data filenames for G(K,iw_n) and Sig(K,iw_n)
        write(6,"('   Enter output data filename for bootstrap G(K,iw_n):  ',$)")
        read(5,"(a72)") linemc
        open(unit=9,file=linemc,status='unknown')
        write(6,"('   Enter output data filename for bootstrap Sig(K,iw_n):  ',$)")
        read(5,"(a72)") linemc
        open(unit=29,file=linemc,status='unknown')


c*******************************************************************************
c       Determine which data are to be read, depending on iopt
c       iopt=0, read the cluster single-particle Green function Gw
c       iopt=1, read the cluster spin susceptibility
c       iopt=2, read the cluster charge susceptibility
c       iopt=3, read the cluster CC/B1G susceptibility
c       iopt=4, read the cluster Sc susceptibility
c       iopt=5, read the cluster B2G susceptibility
c       iopt=6, read the cluster spin susceptibility for O NMR calculation
c       iopt=7, read the cluster anomalous single-particle Green function Fw 
c*******************************************************************************
        write(6,"('G(0), CHS(1), CHC(2), CC(3), SC(4),B2(5)?:  ',$)")
        read(5,*) iopt

c*******************************************************************************
c       Determine how many files are to be read.  If only one, set 
c       the number of bins to be used.  If nbinsuse> run, the number
c       of bins in the file, then all of the bins are used from the
c       data file.
c*******************************************************************************
        write(6,"('                    How many input files?:  ',$)")
        read(5,*) nfiles
        nbinsuse=0
        if(nfiles.eq.1)then
          write(6,"('            How many bins are read in from input file?:  ',$)")
          read(5,*) nbinsuse
        end if
        write(6,"('                    How many bootstrap (bs) bins will be output?:  ',$)")
        read(5,*) nbinsbs
       
c*******************************************************************************
c       Loop over the input data files, read in the data and learn the
c       the average sign and the number of processors used.
c*******************************************************************************

        mcskip=0.0
        meas=0.0
        tprime=0.0
        ed=0.0
        U=0.0
        beta=0.0
        numbermeasprev=0
        do 999 nfile=1,nfiles  
       
          write(6,"('                Enter input data filename:  ',$)")
          read(5,"(a72)") linemc
          open(unit=92,file=linemc,status='old')
          write(6,"(' Enter the associate sigma filename:  ',$)")
          read(5,"(a72)") linemc
          open(unit=82,file=linemc,status='old')

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
            read(92,*) nl,nuse,cluster,ndim
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
            allocate(rawtobs(nbinsbs,nbinsuse),stat=info)
            allocate(dsbsob(2*nuse,Nc,nbinsbs),stat=info)
            allocate(dsbssgn(nbinsbs),stat=info)
            allocate(G_bs(2*nuse,Nc),stat=info)
            allocate(Sig_bs(2*nuse,Nc),stat=info)
            allocate(gsmeas(2*nuse,Nc),stat=info)
            allocate(GcKfsc(2*nuse,Nc),stat=info)
            
            call tables
          end if
 
c*******************************************************************************
c         Read in all the lines containing Gw, and try to 
c         determine how many processors were used and the average sign.
c*******************************************************************************
          nproc=0; i1=0; ave=0.0
          do             
            read(92,"(a72)",end=444) linemc
            if(index(linemc(1:72),'bin of  Gw').ne.0) then
              i1=i1+1
              read(linemc,"(35x,i15)") proc
              if(proc.gt.nproc) nproc=proc
              read(92,*) ((gsmeas(j,icr), j=1,2*nuse),icr=1,Nc),smeas 
              ave=ave+smeas
            end if
          end do
 444      nproc=nproc+1
          numbermeas=i1
          if(nbinsuse.eq.0) nbinsuse=numbermeas  ! if no number was entered, use them all
          nbinsuse=min(nbinsuse,numbermeas)      ! if nbinsuse was entered, its max is numbermeas
          measperproc = numbermeas/nproc   
          allocate(indexinproc(nproc),stat=info) 
          indexinproc(:)=0
          rewind(92)
          write(6,*) 'number of processors used=',nproc
          nbinsuse_per_proc=max(1,nbinsuse/nproc)
          write(6,*) 'bins used per processor=',nbinsuse_per_proc
          write(6,*) 'average sign', ave/float(i1)

c****************************************************************************************************
c         Now read physical parameters from the binned data files, same as the header in binsw file
c****************************************************************************************************

          mcskipl=mcskip
          measl=meas
 401      read(92,"(a72)") linemc
          if(index(linemc(1:72),'warm').ne.0 ) then
            read(92,*) warm,mcskip,meas,ntr,run,niter,nuse1
          else
           goto 401
          endif

          tprimel=tprime; edl=ed; Ul=U; betal=beta
 402      read(92,"(a72)") linemc
          if(index(linemc(1:72),'tprime').ne.0 ) then
             if(index(linemc(1:72),'tdprime').ne.0 ) then
                if(index(linemc(1:72),'t3prime').ne.0 ) then
                   read(92,*) ed,tprime,tdprime,t3prime,U,beta
                else
                   read(92,*) ed,tprime,tdprime,U,beta
                   t3prime=0.d0
                endif 
             else
                read(92,*) ed,tprime,U,beta
                tdprime=0.d0
                t3prime=0.d0
             endif
            temp=1.0/beta
          else
           goto 402
          endif

                   
c***************************************************************************
c        Check if the parameters are consistent with previous and files.          
c***************************************************************************
         if(nfile.gt.1)then
          if(tprimel.ne.tprime.or.edl.ne.ed.or.Ul.ne.U.or.
     1       mcskip.ne.mcskipl.or.
     2       betal.ne.beta.or.nll.ne.nl.or.measl.ne.meas) then
            write(6,*) 'This data set is inconsistent with prior set'
            stop
          end if
          end if

c**************************************
c        Construct bootstrap resample table (rawtobs matrix). Notice that how to bootstrap the raw 
c        data does not depend on the value of raw data, only need the random 
c        number generator to assign the number for rawtobs matrix:
c                           (             )
c                           (   rawtobs   ) bootstrap resample
c                           (   matrix    ) number labeled
c                           (             ) nbinsbs row
c                           (             ) 
c             raw sample number labeled nbinsuse column
c  example: raw sample      (1 1 1 1 1 1 1)
c         bootstrap sample  (1 0 0 3 0 1 2) number of reasamples = 7 = number of raw samples
c**************************************
         rawtobs(:,:)=0
         do i=1, nbinsbs
         do j=1, nbinsuse
           ranclmn=NINT(DBLE(nbinsuse)*rang(iseed)+0.5)
           rawtobs(i,ranclmn)=rawtobs(i,ranclmn)+1
         enddo
         enddo

c*******************************************************************************
c         Now read in the binned data.  It is assumed to be stored
c         as Gs(K,iw_n),s.  Here, Gs is the average product of the
c         Green function and sign and s is the average sign, for a 
c         bin of data.  Note that we want the last datum to always be s.
c         
c         Then construct the bootstrap resamples and output. In the bootstrap, 
c         we use the input Gs, s bin data as the data list, and compute the 
c         object of <G>=<Gs>/<s>=<Gs_bs>/<s_bs>=<G_bs> in bootstrap resamples Gs_bs 
c         and s_bs. Since <G_bs> already reflects the effect of sign, so we fix the
c         sign of <G_bs> as 1 in the output file.
c*******************************************************************************

            indexinproc(:)=0 
            do 16 i=1,nbinsuse                     
c             find a bin header
 409          read(92,"(a72)",end=16) linemc
              if(index(linemc(1:72),'bin of  Gw').ne.0)
     &        then
                read(linemc,"(35x,i20)") proc
                if(proc+1.gt.nproc)then
                  write(6,*)'nproc wrong=',nproc
                  stop 
                endif
                indexinproc(proc+1)=indexinproc(proc+1)+1
                indexofmeas= indexinproc(proc+1)+ proc*measperproc
                ibinn = numbermeasprev + indexofmeas
                 
                write(256,*)'1 indexinproc',indexinproc(proc+1),'indexofmeas',indexofmeas,'ibinn',ibinn
                read(92,*) ((gsmeas(j,icr), j=1,2*nuse),icr=1,Nc),smeas 
              else
               goto 409
              endif
 16         continue
            rewind(92)
            write(6,*) 'done first read for proc',iproc
          
c         Now run through all the data in ifile and read in the measurement.
c         Then use the measurement to construct bootstrap resamples.

            indexinproc(:)=0
            dsbsob(:,:,:)=0
            dsbssgn(:)=0

            do 15 i=1,nbinsuse             
c             find a bin header
 408          read(92,"(a72)",end=15) linemc
              if(index(linemc(1:72),'bin of  Chi' ).ne.0.and.iopt.eq.1.or.
     &           index(linemc(1:72),'bin of  Cha' ).ne.0.and.iopt.eq.2.or.
     &           index(linemc(1:72),'bin of  CC'  ).ne.0.and.iopt.eq.3.or.
     &           index(linemc(1:72),'bin of  Sc'  ).ne.0.and.iopt.eq.4.or.
     &           index(linemc(1:72),'bin of  B2'  ).ne.0.and.iopt.eq.5.or.
     &           index(linemc(1:72),'bin of  Chi' ).ne.0.and.iopt.eq.6.or.
     &           index(linemc(1:72),'bin of  Fw').ne.0.and.iopt.eq.7.or.
     &           index(linemc(1:72),'bin of  Gw').ne.0.and.iopt.eq.0 )
     &        then
                read(linemc,"(35x,i20)") proc	
                   indexinproc(proc+1)=indexinproc(proc+1)+1
                   indexofmeas= indexinproc(proc+1) + proc*measperproc
                   ibin = numbermeasprev + indexofmeas
              write(256,*)'2 indexinproc',indexinproc(proc+1),'indexofmeas',indexofmeas,'ibin',ibin
c                  ibin=ibin+1
                   read(92,*) ((gsmeas(j,icr), j=1,2*nuse),icr=1,Nc),smeas 
 
                   do j=1,nbinsbs
c****************************************************              
c              now construct the bootstrap resamples using the rawtobs matrix
c              and compute the object of (similar to <G>=<Gs>/<s>) <G_bs>=<Gs_bs>/<s_bs> in 
c              bootstrap resamples Gs_bs and s_bs. This step consider the sign effect, after
c              it, no sign anymore, the sign in the resamples will always be assigned as 1.     
c****************************************************
                     dsbsob(:,:,j)=dsbsob(:,:,j)+DBLE(rawtobs(j,i))*gsmeas(:,:) 
                     dsbssgn(j)=dsbssgn(j)+DBLE(rawtobs(j,i))*smeas
                   enddo    

              else
                goto 408
              endif !index
 15         continue
          numbermeasprev = numbermeasprev + numbermeas
 999    continue ! nfile      
        nbins=ibin
        if(nbins.gt.nbinsmax) stop 'in rddca.f increase nbinsmax'

c*******************************************************************************
c       Now output G_bs (No sgn any more, always fix sgn=1 at the end of G_bs(k,iw_n).)
c*******************************************************************************

c       write out fileheader binswbs, the same as that in ctqmc sumup.f
        write(9,*)' '

        do j=1,nbinsbs
c     write out G_bs(k,iw_n), sgn(fix to 1)   , the same as that in ctqmc sumup.f
          write(9,*)dsbsob(:,:,j)/dsbssgn(j),1.0d0
        enddo

c*******************************************************************************
c       Now output Sig_bs (No sgn any more, always fix sgn=1 at the end of Sig_bs(k,iw_n).)
c*******************************************************************************
c       write out fileheader binswbssig, the same as that in ctqmc sumup.f
        write(29,*)' '
     

c************************
c       read in GcKfsc from sigma.dat
c************************
 2042   read(82,1000) linemc
 1000   format(a72)
        if(index(linemc(1:72),'GcKfsc').ne.0 ) then
          do ick=1,Nc
          do n=1,2*nuse
            read(9,*) i,j,GcKfsc(n,ick)           
          end do
          end do
        else
          goto 2042
        endif

        do j=1,nbinsbs
c       Apply Dyson equation: Sig(K,iw_n)=G_0(K,iw_n)^-1-G(K,iw_n)^-1
         Sig_bs(:,:)=1.0d0/GcKfsc(:,:)-1.0d0/(dsbsob(:,:,j)/dsbssgn(j))
c       write out G_bs(k,iw_n), sgn(fix to 1)   , the same as that in ctqmc sumup.f
        write(9,*)Sig_bs(:,:),1.0d0
        enddo


c*******************************************************************************
c       Write some stuff out
c*******************************************************************************

        write(6,*) 'nbins=',nbins
        write(6,*) '   nl=',nl
        write(6,*) '  '
        write(6,*) 'warm,mcskip,meas,ntr,run'
        write(6,*) warm,mcskip,meas,ntr,nbins
        write(6,*) 'tprime,ed,U,beta'
        write(6,*) tprime,ed,U,beta

 

c*******************************************************************************
c       Wrap things up
c*******************************************************************************
        call deallocate_1p      
            deallocate(rawtobs )
            deallocate(dsbsob  )
            deallocate(dsbssgn )
            deallocate(G_bs    )
            deallocate(Sig_bs  )
            deallocate(gsmeas  )                        
        stop               
        end     


!======================================================
      real*8 function  rang(iq)

      integer iq
      real*8 pi, ranmod, theta, ranf

      PI = 3.1415926536D0
      RANMOD = SQRT(-2.D0 * LOG(RANF(iq)))
      THETA  = 2.D0 * PI * RANF(iq)
      rang = RANMOD * COS(THETA)

      return
      end
!======================================================           

