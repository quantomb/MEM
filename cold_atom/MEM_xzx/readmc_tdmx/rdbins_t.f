      program rdbins
c This code is used to process the binned data into data which can be 
c used by readmc, the preprocessor of the analytic continuation code.
c This this code takes the binned data, and produces the averages, 
c and covariances.  The binn size is variable within the code, i.e. 
c bins can be averaged, and multiple data files can be opened and 
c averaged.  In addition, rdbins is capable of qualifying the data 
c by performing a variety of tests.   These include tests to measure 
c the kurtosis and skew of the distribution, and a calculation of 
c the probability that the data is gaussian (normally distributed). 
c This last probability should be at least 0.5.
c
c f77 -r8 -O -o rdbins rdbins.for  

      implicit none
      integer nl,i,j,k,l,skip,nuse,ind,naver,ihisto,isymm
      integer nlmax,nbinsmax,nbins,nbinsuse,ibin,i_tmp
      integer nfiles,nfile,nl_bin,matsu_freq
      integer matsu_flag,sus_flag
      parameter(nlmax=200,nbinsmax=400000)
      real*4 beta,tau,dtau
      real*4 data(nlmax,nbinsmax),data2(nbinsmax),aver(nlmax)
! sgn_avg: averaged sign
! k_avg  : averaged hopping term
! d_avg  : averaged double occupancy
      real*4 sgn_avg, k_avg, d_avg, avg_tmp, sgn_tmp
      real*4 r1
      real*4 ave,adev,sdev,var,skew,curt,chsq,prob,error
      real*4 store(nlmax),ds(nlmax,nbinsmax)
      real*4 PI,normalor,tmp_data(nlmax,6)
      parameter(PI=3.1415926535897932354626d0)
      character*72 linemc
      character*5 chr_tmp

      write(6,11) 
 11   format('               enter output data filename:  ',$)
      read(5,201) linemc
      if(index(linemc(1:72),'bins').ne.0 )then
         write(6,*) 'This is an input file'
         stop
      end if
      open(unit=9,file=linemc,status='unknown')

      write(6,12) 
 12   format('          enter the number of input files:  ',$)
      read(5,*) nfiles

      if(nfiles.eq.1)then
          write(6,13) 
 13       format('      enter the number of bins to be used:  ',$)
          read(5,*) nbinsuse
      else
          nbinsuse=nbinsmax
      end if

c     skip is used when the data is too dependent.  If
c     skip is greater than one, some of the time slices in
c     each bin will be skipped.  Typically (now) skip=1
      write(6,14) 
 14   format('                               enter skip:  ',$)
      read(5,*) skip

      write(6,15) 
 15   format('enter course grain size for bin averaging:  ',$)
      read(5,*) naver

      write(6,61) 
 61   format('Using Matsubara frequencies ?(y/n) :',$)
      read(5,*) linemc
      if(index(linemc(1:72),'y').ne.0) then
         matsu_flag=1
      endif
      if(matsu_flag.eq.1) then 
      write(6,66)
 66   format('enter number of Matsubara frequencies:  ',$)
      read(5,*) matsu_freq
      endif

      write(6,80)
 80   format('Which suscept. ?(1:<KK>, 2:<KD+DK>,3:<DD>,4:<EE>) :',$)
      read(5,*) sus_flag

      ihisto=0
      write(6,16) 
 16   format('                  Print histograms? (y/n):  ',$)
      read(5,201) linemc 
      if(index(linemc(1:72),'y').ne.0) ihisto=1

      isymm=0
      write(6,17) 
 17   format('    Should the data be symmetrised? (y/n):  ',$)
      read(5,201) linemc 
      if(index(linemc(1:72),'y').ne.0) then
         isymm=1
         write(9,*) 'symmetric storage mode'
      end if

      nbins=0
      ibin=0
      do 999 nfile=1,nfiles
c
         write(6,18) 
 18   format('            enter the input data filename:  ',$)
         read(5,201) linemc
         open(unit=92,file=linemc,status='old')
!----------------------------------------------------------------------
! readin file header

 400     read(92,201) linemc
         if(index(linemc(1:72),'Results').ne.0 ) then
            read(92,*) nl,nl_bin,beta
! set number of time slices 
            nuse=(nl_bin-1)/skip+1
            if(matsu_flag.eq.1) then
! for Matsubara frequency case, dtau is w_1
               dtau = 2.d0*PI/beta
            else
! set dtau for imaginary time case
               dtau=beta/float(nl)
            endif
! if symmetrize data, change number of time slice
            if (isymm.eq.1) nuse=nuse/2+1
         else
            goto 400
         endif

!----------------------------------------------------------------------
! readin binned data

         sgn_avg = 0.0
         k_avg   = 0.0
         d_avg   = 0.0

         do i=1,nbinsuse
! find a bin header
 407        read(92,201,end=410) linemc
            if(index(linemc(1:72),'bin').ne.0) then
               ibin=ibin+1
! readin one bin
               read(92,*) (i_tmp,store(j),j=1,nl_bin)
!               read(92,*) (store(j),j=1,nl_bin)
! readin sgn and accumulate
               read(92,*) chr_tmp, sgn_tmp
               sgn_avg = sgn_avg + sgn_tmp
! readin <K>
               read(92,*) chr_tmp, avg_tmp
               k_avg = k_avg + avg_tmp
! readin <D>
               read(92,*) chr_tmp, avg_tmp
               d_avg = d_avg + avg_tmp
! transfer data to array ds and do average <X.s>/<s>
               if(isymm.eq.0) then
                  do j=1,nuse
                     !ds(j,ibin)=store(1+skip*(j-1))
                     ds(j,ibin)=store(1+skip*(j-1))/sgn_tmp
                     data(j,ibin)=0.0
                  enddo
               else
! symmetrize data
                  ds(1,ibin)=store(1)/sgn_tmp
                  do j=2,nuse
!                     ds(j,ibin)=0.5*(store(1+skip*(j-1))+
!     &store(1+skip*(nl+1-j))  )  ! j --> nl+2-j
                     ds(j,ibin)=0.5*(store(1+skip*(j-1))+
     &store(1+skip*(nl+1-j))  )/sgn_tmp
                     data(j,ibin)=0.0
                  enddo
               end if
            else
               goto 407
            endif
         end do

 410     continue ! end of read in one file

 999  continue ! end of read in all file

!----------------------------------------------------------------------
! Normalizing data

! specify how many matsubara frequencies used
      if(matsu_flag.eq.1) then
         nuse = matsu_freq
      endif
! current data bin number
      nbins=ibin
      nbins=nbins/naver
      write(6,19) nbins
 19   format('                                    nbins:  ',i3)
      write(6,20) nl
 20   format('                                       nl:  ',i3)

c      write(9,*) 'Results for nl,nbins,beta'
c      if(matsu_flag.eq.1) then
c         write(9,*) nuse,nbins,beta
c      else
c         write(9,*) (nl_bin-1)/skip+1,nbins,beta
c      endif
c      write(9,*) ' '

! Do average of sign, <K> and <D>, check sign
      sgn_avg = sgn_avg/REAL(ibin)
      k_avg   = k_avg  /REAL(ibin)
      d_avg   = d_avg  /REAL(ibin)
      write(6,*) '<Sign>=',sgn_avg
      write(6,*) '<K>   =',k_avg
      write(6,*) '<D>   =',d_avg
      if(sgn_avg.le.0.01) then
         write(6,*) 'Woops, the sign is too small to be:',sgn_avg
         pause 'Continue or not?'
      endif


! Now do bin average for rebbined data.  
! The data in bins ibin > naver*int(nbins/naver)
! is not used.
      ind=0
      do i=1,nbins
         do j=1,naver
            ind=ind+1
            do k=1,nuse
               data(k,i)=data(k,i)+ds(k,ind)/float(naver)
            end do
         end do
      end do
!
! remove the vacuum contribution <K>**2, <K><D>, or <D>**2, <E>**2
      if(sus_flag.eq.1) then
         avg_tmp = k_avg*k_avg ! <K><K>
      else if(sus_flag.eq.2) then
         avg_tmp = 2.d0*k_avg*d_avg ! <K><D>
      else if(sus_flag.eq.3) then
         avg_tmp = d_avg*d_avg ! <D><D>
      else if(sus_flag.eq.4) then
         avg_tmp = k_avg*k_avg-2.d0*k_avg*d_avg + d_avg*d_avg!<E><E>
      endif
      do i = 1, nbins
      do l = 1, nuse
         data(l,i) = data(l,i) - avg_tmp
      enddo
      enddo
!
! normalize the data !!
! For Matsubara frequency case!
!      if(matsu_flag.eq.1) then
!         normalor = 0.d0
!         do i = 1, nbins
!            normalor = normalor + data(1,i)
!         enddo
!         normalor = normalor/DBLE(nbins)
!         do l = 1, nuse
!            do i = 1, nbins
!               data(l,i) = data(l,i)/normalor
!            enddo
!         enddo
!      else
! For imaginary time case!
         normalor=0.d0
         r1=2.d0/3.d0
         do l=1,nl
            normalor=normalor+r1*dtau*SUM(data(l,1:nbins))
            r1=2.d0-r1
         enddo
         normalor=normalor/float(nbins)              ! << s*chi(T) >>
!         sus1=sum(data(nl+1,1:nbins))/float(nbins) ! << s >>
!          
!         chi(T) = sus/sus1
! As we have obtained sgn_avg = sus1
!         sus2=sum(data(1,1:nbins))/float(nbins)    ! << s*chi(tau=0)>>
!         data(nuse+1,1:nbins)=                     ! <s> * <<s*chi(T) >>/<< s >>
!     &data(nl+1,1:nbins)*sus/sgn_avg
         data=data/normalor                         ! <s*chi(tau)>/<< s*chi(T) >>  for l<=nuse
                                                    ! <s>/<<s>> for l=nuse+1
!      endif
c
c          Now 
c              <<data(nuse+1,:)>> = <s> * << s*chi(T) >>/<< s >>
c          and
c              <<data(l,:)>> = <<s chi(tau)>> / << s*chi(T) >>   for l<= nuse
c
c          so that the integral of the spectrum corresponding to thisdata
c          is one.


!----------------------------------------------------------------------
c Now analyze the data (the analysis tools are described below)
c 
c      write(9,*) '   tau      Gtau(n)         +/-         curt  skew  
c     &prob'
      do l=1,nuse
         do i=1,nbins
            data2(i)=data(l,i)
         end do
         call moment(data2,nbins,ave,adev,sdev,var,skew,curt)
c         tau=dtau*float(l)
         tau=dtau*float(l-1)*float(skip)
         aver(l)=ave
         if(adev.gt.1.0e-20) then
            call histo(data2,nbinsmax,nbins,ave,var,sdev,chsq,
     &prob,ihisto)
         end if
         error=sdev/sqrt(float(nbins))
         tmp_data(l,1) = tau
         tmp_data(l,2) = ave
         tmp_data(l,3) = error
         tmp_data(l,4) = curt
         tmp_data(l,5) = skew
         tmp_data(l,6) = prob
c         write(9,202) tau,ave,error,curt,skew,prob
      end do
c
      write(9,*) 'Results for nl,nbins,beta'
      if(matsu_flag.eq.1) then
         write(9,*) nuse,nbins,beta,normalor
      else
         write(9,*) nuse,nbins,beta,normalor
      endif
      write(9,*) ' '
      write(9,*) '   tau           Gtau(n)            +/-           curt
     &      skew         prob'
      do l = 1, nuse
         write(9,202) (tmp_data(l,i),i=1,6)
      enddo
c
      write(9,*) ' '
      write(9,*) ' Covariance data '
      write(9,*) '   i     j     Gij-GiGj        +/-'
      call covar(data,nlmax,nbinsmax,nuse,nbins,aver,9)
 202  format(x,e14.7,x,e14.7,x,e14.7,x,e12.6,x,e12.6,x,e12.6)
 201  format(a72)

      stop
      end


c This series of subroutine calls is used to process
c binned data stored in an array meas(index,binnumber).
c
c     Subroutine moment calculated the moments of the data
c for a particular index.  The data for the particular
c index must be transferred to a 1-d array: 
c
c data(:)=meas(index,:)
c
c Then execute
c
c call moment(data,nbins,ave,adev,sdev,var,skew,curt)

c where nbins is the number of bins of data.  On output
c ave is the bin average, adev the average deviation,
c sdev the standard deviation, var the variance, skew
c the skew and curt the kurtosis of the distribution.
c (for a more detailed discussion of curt and skew, see
c Numerical Recipes)

c
c
c One can also generate a histogram of this same data, and
c compare it to a gaussian distribution.  On input, nbinsmax
c is the declared dimension of data, nbins the number of bins,
c ave, var, and sdev that calculated by moments.  On output
c prob is the probability that the binned data is normal
c (fora discussion of chsq see Numerical Recipes).  prob
c should be about 0.5 or bigger.
 
c call histo(data,nbinsmax,nbins,ave,var,sdev,chsq,prob)

c If you generate the binned data in two ways, say shuffled and
c sequential binning, then you can compare the two data sets
c with an FTEST.  If the binning procedure removed all the 
c corellations between adjacent bins, then the two data set
c are identical, and the calculated probability would be 1.

c Use the ftest on the two data sets of G(tau)
c      write(6,*) ' '
c      write(6,*) ' FTEST DATA '
c      write(6,*) ' '
c      write(6,*) '   l      f -----> prob  '
c      do 35 i=1,nlg
c        data(:)=gmeas1(i,:)
c        data2(:)=gmeas2(i,:)
c        call ftest(data,run,data2,run,f,prob)
c        write(6,87) i,f,prob
c 35  continue
c 87  format(1x,i3,3x,f7.4,3x,f8.5)

c Of course in order to analytically continue numerical data
c the covariance matrix is necessary.  The following takes
c the data in meas, and the associated vector of averages
c ave(index) and calculate and writes to standard output
c the covariance matrix.

c Now calculate the covariance using 
c covar(meas,nm1,n1,n2,ave,nm2)
c
c Where 
c meas(nm1,nm2) contains the measurements
c nm1 declared dimension of meas
c nm2 declared dimension of meas
c n1 dimension of meas and ave
c n2 dimension of meas
c ave average over runs measurments
c
c Thus, with partilce-hole symmetry, we could take
c n1=nlg/2 in order to save space

c call covar(meas,indexmax,nl,nbins,ave,binsmax)


      SUBROUTINE MOMENT(DATA,N,AVE,ADEV,SDEV,VAR,SKEW,CURT) 
      DIMENSION DATA(N) 
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
      RETURN 
      END 

      SUBROUTINE FTEST(DATA1,N1,DATA2,N2,F,PROB)
      DIMENSION DATA1(N1),DATA2(N2)
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
      RETURN
      END


      SUBROUTINE AVEVAR(DATA,N,AVE,VAR)
      DIMENSION DATA(N)
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
      RETURN
      END
      
      FUNCTION BETAI(A,B,X)
      IF(X.LT.0..OR.X.GT.1.)then
        write(6,*) 'bad argument X in BETAI'
        return
      end if
      IF(X.EQ.0..OR.X.EQ.1.)THEN
         BT=0.
      ELSE
         BT=EXP(GAMMLN(A+B)-GAMMLN(A)-GAMMLN(B)
     *        +A*ALOG(X)+B*ALOG(1.-X))
      ENDIF
      IF(X.LT.(A+1.)/(A+B+2.))THEN
         BETAI=BT*BETACF(A,B,X)/A
         RETURN
      ELSE
         BETAI=1.-BT*BETACF(B,A,1.-X)/B
         RETURN
      ENDIF
      END
      
      
      FUNCTION BETACF(A,B,X)
      PARAMETER (ITMAX=100,EPS=3.E-7)
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
      RETURN
      END
      
      

      subroutine histo(data,runmax,run,avg,var,sdev,chsq,prob,iq)
      integer run,runmax,iq
      real*4 hist(50), curve(50)
      real*4 data(runmax),width,sdev
      
      do 10 i = 1, 50
         curve(i) = 0.0
         hist(i) = 0.0
 10   continue
      
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
 40   continue
      
      call gcurve(xbeg,delx,nx,avg,var,run,curve)
      
      call chsone(hist,curve,nx,2,df,chsq,prob)
      if(iq.eq.1)call phist(xbeg,delx,nx,hist,curve)
      
      return
      end
      
      
      subroutine gcurve(xbeg,delx,nx,xm,var,ndata,curve)
      dimension curve(*)
      
      sigma = sqrt(var)
      
      do 10 j = 1, nx
         x = xbeg + (j-0.5)*delx
         c = exp(-1.0*(x-xm)**2/(2.*var))/(2.5066*sigma)
         curve(j) = c*ndata*delx
 10   continue
      
      return
      end
      
      SUBROUTINE CHSONE(BINS,EBINS,NBINS,KNSTRN,DF,CHSQ,PROB)
      DIMENSION BINS(NBINS),EBINS(NBINS)
      DF=NBINS-1-KNSTRN
      CHSQ=0.
      DO 11 J=1,NBINS
         IF(EBINS(J).LE.0.) goto 11
         CHSQ=CHSQ+(BINS(J)-EBINS(J))**2/EBINS(J)
 11   CONTINUE
      PROB=GAMMQ(0.5*DF,0.5*CHSQ)
      RETURN
      END

      FUNCTION GAMMQ(A,X)
      IF(X.LT.0..OR.A.LE.0.)return
      IF(X.LT.A+1.)THEN
        CALL GSER(GAMSER,A,X,GLN)
        GAMMQ=1.-GAMSER
      ELSE
        CALL GCF(GAMMCF,A,X,GLN)
        GAMMQ=GAMMCF
      ENDIF
      RETURN
      END
      
      
      SUBROUTINE GCF(GAMMCF,A,X,GLN)
      PARAMETER (ITMAX=100,EPS=3.E-7)
      GLN=GAMMLN(A)
      GOLD=0.
      A0=1.
      A1=X
      B0=0.
      B1=1.
      FAC=1.
      DO 11 N=1,ITMAX
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
 11   CONTINUE
      write(6,*) 'A too large, ITMAX too small'
 1    GAMMCF=EXP(-X+A*ALOG(X)-GLN)*G
      RETURN
      END
      
      
      SUBROUTINE GSER(GAMSER,A,X,GLN)
      PARAMETER (ITMAX=100,EPS=3.E-7)
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
      RETURN
      END
      
      
      FUNCTION GAMMLN(XX)
      REAL*4 COF(6),STP,HALF,ONE,FPF,X,TMP,SER
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
      RETURN
      END

      subroutine phist(xbeg,delx,nx,hist,curve)
      character*1 blank, cross, ast
      character*50 zline
      dimension hist(*), curve(*)
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
         zline(k:k) = ast
         write(6,1000) x, nint(hist(j)), zline
 40   continue
      
 1000 format(1h ,f7.3,i5,2h  ,(a))
      return
      end

      subroutine covar(data,i1,i2,n1,n2,ave,lun)
      integer n1,n2,i,j,k,i1,i2,lun
      real*4 data(i1,i2),ave(i1),cij,cij2,sigg
c thic code calculates the cavaraince of data(n1,n2) where
c n1 is the data index and n2 is the bin index.  ave(n1)
c contains the average values of the measurement.  The covariance,
c and its error are printed out to logical unit 9.

      do 50 i=1,n1
      do 50 j=i,n1
         cij=0.0
         cij2=0.0
         do 55 k=1,n2
            cij=cij + (data(i,k)-ave(i))*(data(j,k)-ave(j))
            cij2=cij2+((data(i,k)-ave(i))*(data(j,k)-ave(j)))**2
 55      continue
         cij=cij/float(n2)
         cij2=cij2/float(n2)
         sigg=(cij2-cij**2)/float(n2-1)
c         if(sigg.lt.0.0) write(6,*) ' - variance of cij ',sigg
         sigg=sqrt(abs(sigg))
         write(lun,90) i,j,cij,sigg
 50   continue
 90   format(1x,i4,2x,i4,2x,e14.7,2x,e14.7)

      return
      end

