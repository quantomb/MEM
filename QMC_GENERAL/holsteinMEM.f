
c     MOFIFIED 05/09/05 TO DUM G(q,tau).  ALSO DO X/Y AVERAGING.

c     MODIFIED 05/05/05 FOR DUMPING GREENS FUNCTION FOR ANALYTIC CONTINUATION
c     [x] REMOVED EXTERNAL STATEMENTS
c     [x] FIXED MPHASE DIMENSIONING
c     [x] INSERT nbin AND CHECK GET SAME RESULTS WITH nbin=10
c     [x] CHECK GET SAME RESULTS FOR AVERAGES WITH nbin NOT 10
c     [x] OPEN NEW AC FILE
c     [x] WRITE bdnl AND bdql TO IT
c     [x] CHECK THEY ARE RIGHT BY POST-AVERAGING
c     [x] RUN WITH BOUNDS CHECKING

c     From niyaz@solid.ucdavis.edu Tue Nov 17 15:55:31 1992


c        The main program for the two-dimensional Holstein model with
c        positive U.  Antecedent is Steve's code exact9.for.
c        We have modified spinzz and zz components of af structure factor
c        and susceptibility to measure cdw correlations instead.  Transverse
c        measurements of spin correlations remain.

c****************************************************
        module global
        implicit none
        save
        integer, parameter:: n=4 , l=32
        integer, parameter:: volume=n*n*l,toff=n*n 
        integer, parameter:: nbin=200 !measurements per bin

        integer,  parameter :: tausk = 10
        integer,  parameter :: phonskip = 5
 


        end module global
c****************************************************
        program holstein4

c pars.h -- basic Holstein parameters

        use global
        implicit none

         integer count
         common/acphononi/count
         real*8 adnl2(0:n/2,0:n/2,0:l),asgnph2
         common/acphononr/adnl2,asgnph2

        integer j1,j2,mx,my,lx,ly,ti,mpx,mpy,kx,ky,kmx,kmy
c integers.h -- integer parameters for holstein.
        integer warms,sweeps,msr,nwrap
        integer iran
        common/integers/warms,sweeps,msr,nwrap, iran
c couple.h -- coupling parameters for holstein.
         real*8 move,t,omega,g,mu,dtau,expmu,gam,lambda,tdtau,dens
         real*8 pe,kee
         common/couple/move,t,omega,g,mu,
     1       dtau,expmu,gam,lambda,tdtau,dens,pe,kee
c ivectors.h -- index vectors for indirect addressing and phase vector
        integer xplus(0:toff-1),xminus(0:toff-1)
        integer yplus(0:toff-1),yminus(0:toff-1)

        common/ivectors/xplus,xminus,yplus,yminus
c vectors.h -- vectors for 2d holstein model.
        real*8 hub(0:volume-1)
        real*8 vup(0:volume-1),vdn(0:volume-1)
        real*8 mphase(0:toff)

        common/vectors/hub,vup,vdn, mphase
c m0var.h -- common for measurements. (2dmoo)
        real*8 sgnup,sgndn,sgn,nup,ndn,nud,ke,detup
        real*8 saf,scdw,sferro,sfer2,grfun(0:n/2,0:n/2),safsq,safsq2
        real*8 spinxx(0:n/2,0:n/2),cdw(0:n/2,0:n/2)
        real*8 pairmat(-1:1,-1:1,-1:1,-1:1),errrat,asgnp
        integer nmeas0,redo,noredo,dopair,nmeasp,numpair,dowrite
        real*8 asgnup,asgndn,asgn,anup,andn,anud,ake,asafsq,asafsq2
        real*8 asaf,ascdw,asferro,asfer2,agrfun(0:n/2,0:n/2)
        real*8 aspinxx(0:n/2,0:n/2),acdw(0:n/2,0:n/2)
        real*8 apairmat(-1:1,-1:1,-1:1,-1:1)
        real*8 den(0:n/2,0:n/2,0:1),aden(0:n/2,0:n/2,0:1)
        common/m0var/ sgnup,sgndn,sgn,nup,ndn,nud,ke,safsq,safsq2,
     1       saf,scdw,sferro,sfer2,grfun,
     1       spinxx,cdw,aspinxx,acdw,asafsq,asafsq2,
     2       asgnup,asgndn,asgn,anup,andn,anud,ake,den,aden,
     3       asaf,ascdw,asferro,asfer2,agrfun,errrat,asgnp
     4       ,pairmat,apairmat
        common/m0vari/nmeas0,redo,noredo,dopair,nmeasp,numpair,dowrite

c mtauvar.h -- common for tau-dependent measurements. (2dmoo)
        real*8 gnl(0:n/2,0:n/2,0:l),agnl(0:n/2,0:n/2,0:l)
        real*8 chinl(0:n/2,0:n/2,0:l),achinl(0:n/2,0:n/2,0:l)
        real*8 pairsus(-1:1,-1:1,-1:1,-1:1,0:l)
        real*8 apairsus(-1:1,-1:1,-1:1,-1:1,0:l),asgnt
        integer nmeast
        common/mtauvar/gnl,agnl,chinl,achinl,pairsus,
     1       apairsus,asgnt,nmeast
c
c getgpam.h -- parameters for getg (2dmoo)
c
        integer orthlen,doauto,torth
        real*8 eorth,difflim,errpam
        common/getgpam/eorth,difflim,errpam,orthlen,doauto,
     1                torth
c newm.h -- Common. (2dmoo)
        real*8        peep,apeep,phpe,aphpe,phke,aphke,phnn,aphnn
        common/newm/peep,apeep,phpe,aphpe,phke,aphke,phnn,aphnn
c cdwm.h -- Another silly common. (2dmoo)
        real*8 achicdw,chicdw
        common/cdwm/achicdw,chicdw
c names.h -- common for I/O names.
        character*80 outname, inlat
        common/names/outname, inlat
c     phonon.h -- variables for the phonon propagator.
      real*8 dnl(0:n/2,0:n/2,0:l),adnl(0:n/2,0:n/2,0:l),asgnph
      common/phonon/dnl,adnl,asgnph

c impure.h -- common block for impurity.
      real*8 delta0, delta(0:toff-1), expdel(0:toff-1)

      common/impure/ delta0, delta, expdel

        real*8 twopi
        parameter(twopi=6.283185307179586)

        integer i,j,k,accept,reject,wraps
	     integer accept2,reject2

        real*8 time,aveval,errval,aveval2,errval2
        complex*16 cveval,crrval,cveval2,crrval2

        real*8 bpairv(nbin,9),bpairl(nbin,9,0:l),gpairl(nbin,9,0:l)
        complex*16 bpairw(nbin,9,0:l),gpairw(nbin,9,0:l)
        real*8 waves(-1:1,-1:1,9)
        data waves/0.d0, 0.d0, 0.d0, 0.d0,1.d0,0.d0,0.d0,0.d0 , 0.d0,
     1        0.d0, 0.5d0, 0.d0 , 0.5d0,0.d0,0.5d0,0.d0 , 0.5d0, 0.d0,
     2        0.d0,-0.5d0, 0.d0 , 0.5d0,0.d0,0.5d0,0.d0 ,-0.5d0, 0.d0,
     3        0.5d0, 0.d0, 0.5d0, 0.d0 ,0.d0,0.d0 ,0.5d0, 0.d0 , 0.5d0,
     4       -0.5d0, 0.d0, 0.5d0, 0.d0 ,0.d0,0.d0 ,0.5d0, 0.d0 , -0.5d0,
     5        0.d0, 0.d0 , 0.d0 ,-0.5d0,0.d0,0.5d0,0.d0 , 0.d0 , 0.d0,
     6        0.d0,-0.5d0, 0.d0 , 0.d0 ,0.d0,0.d0 ,0.d0 , 0.5d0, 0.d0,
     7       -0.5d0, 0.d0, 0.d0 , 0.d0 ,0.d0,0.d0 ,0.d0 , 0.d0 , 0.5d0,
     8        0.d0, 0.d0 ,-0.5d0, 0.d0 ,0.d0,0.d0 ,0.5d0, 0.d0 , 0.d0/

c       *****

        real*8 gql(0:n/2,0:n/2,0:l)
        complex*16 gnw(0:n/2,0:n/2,0:l),gqw(0:n/2,0:n/2,0:l)
        complex*16 sigma(0:n/2,0:n/2,0:l)
        real*8 chiql(0:n/2,0:n/2,0:l)
        complex*16 chinw(0:n/2,0:n/2,0:l),chiqw(0:n/2,0:n/2,0:l)
        real*8 dql(0:n/2,0:n/2,0:l)
        complex*16 dnw(0:n/2,0:n/2,0:l),dqw(0:n/2,0:n/2,0:l)

        real*8 bgnl(0:n/2,0:n/2,0:l),bgql(0:n/2,0:n/2,0:l)
        complex*16 bgnw(0:n/2,0:n/2,0:l),bgqw(0:n/2,0:n/2,0:l)
        complex*16 bsigma(0:n/2,0:n/2,0:l)
        real*8 bchinl(0:n/2,0:n/2,0:l),bchiql(0:n/2,0:n/2,0:l)
        complex*16 bchinw(0:n/2,0:n/2,0:l),bchiqw(0:n/2,0:n/2,0:l)
        real*8 bdnl(0:n/2,0:n/2,0:l),bdql(0:n/2,0:n/2,0:l)
        complex*16 bdnw(0:n/2,0:n/2,0:l),bdqw(0:n/2,0:n/2,0:l)
        complex*16 bpi(0:n/2,0:n/2,0:l)

        real*8 sgnl(0:n/2,0:n/2,0:l),sgql(0:n/2,0:n/2,0:l)
        complex*16 sgnw(0:n/2,0:n/2,0:l),sgqw(0:n/2,0:n/2,0:l)
        complex*16 ssigma(0:n/2,0:n/2,0:l)
        real*8 schinl(0:n/2,0:n/2,0:l),schiql(0:n/2,0:n/2,0:l)
        complex*16 schinw(0:n/2,0:n/2,0:l),schiqw(0:n/2,0:n/2,0:l)
        real*8 sdnl(0:n/2,0:n/2,0:l),sdql(0:n/2,0:n/2,0:l)
        complex*16 sdnw(0:n/2,0:n/2,0:l),sdqw(0:n/2,0:n/2,0:l)
        complex*16 spi(0:n/2,0:n/2,0:l)

        real*8 bsgnup(nbin),bsgndn(nbin),bsgn(nbin)
        real*8 bnup(nbin),bndn(nbin)
        real*8 bntot(nbin),bnud(nbin),bke(nbin),benergy(nbin)
        real*8 bsaf(nbin),bscdw(nbin),bsferro(nbin),bsfer2(nbin)
        real*8 bgrfun(nbin,0:n/2,0:n/2),bsafsq(nbin)
        real*8 bden(nbin,0:n/2,0:n/2,0:1)
        real*8 bspinxx(nbin,0:n/2,0:n/2), bcdw(nbin,0:n/2,0:n/2)
        real*8 tcdw(0:n/2,0:n/2),cdwq(0:n/2,0:n/2)
        real*8 bcdwq(nbin,0:n/2,0:n/2)

        real*8 pairgsus(-1:1,-1:1,-1:1,-1:1,0:l)
        real*8 gmatup(0:n*n,0:n*n),gmatdn(0:n*n,0:n*n),epsk,zbrent
        real*8 getden
        integer nmax
        real*8 bpeep(nbin),bphpe(nbin),bphke(nbin),bphnn(nbin)
        real*8 bchicdw(nbin)
        character*80 rstring, gstring,svstring, m0string, mpstring
        character*80 acstring
	integer numtry
	real*8 gsize

        real*8 axeq,bxeq(nbin)
         common/eqpos/axeq

c        set up the phase vectors
        call phaseset()

        call readin()
	write (6,*) 'numtry,gsize'
        read (5,*) numtry,gsize
        count=0
 
c
c open output files
c
        rstring = 'r'//outname
        gstring = 'g'//outname
        svstring = 'sv'//outname
        m0string = 'm0'//outname
        mpstring = 'mp'//outname
        acstring = 'ac'//outname
        open(unit=66, status='new',file=rstring)
        open(unit=60, status='new',file=svstring)
        open(unit=68, status='new',file=gstring)
        open(unit=88, status='new',file=acstring)
        if (dowrite.eq.1) then
          open(unit=67, status='new',file=m0string)
          open(unit=69, status='new',file=mpstring)
        endif

c
c        write header
c
        write(66,*) 'Version holstein4'
        write(66,*) '7/23 -- Writes certain measurements at each sweep'

        write(66,*) '2/21/91 -- Inputs random site impurities'


        write(66,*) ' '
        write(66,*) 'n=',n,'  l=',l
        write(66,*) ' '
        write(66,*) 'afeps = ',0.0

        write(66,*) 'warms=',warms,'  sweeps=',sweeps

        write(66,*) 't=',t,'  omega=',omega,'  g=',g,' move= ',
     1        move,' dens=',dens

        write(66,*)'dtau=',dtau,' nwrap = ',nwrap,
     1         ' difflim= ', difflim,' errrat= ',errrat
         write (66,*) 'doauto has been zeroed'
        write(66,*)' doauto= ',doauto,' orthlen= ',orthlen,
     1   ' eorth= ',eorth,' dopair= ',dopair,' numpair= ',numpair
        write(66,*)' torth= ',torth
        write(66,*)'errpam is ',errpam
        write(66,*) 'iran=',iran
        write(66,*) 'outname=',outname
        write(66,*) 'inlat=',inlat
        write(66,*) 'tausk=',tausk
        write(66,*) 'phonskip=',phonskip
	write (66,*) 'numtry,gsize'
        write (66,*) numtry,gsize

        write (88,*)' number of bins =', sweeps/nbin/phonskip
        write (88,*)' number of measurements per bin =',nbin
        write(88,*) 'Results for l (time points), n (lattice size nxn)'
        write(88,*)  l , n
        write(88,*) 't=',t,'  omega=',omega,'  g=',g,' move= ',
     1        move,' dens=',dens
        write(88,*)'dtau and beta'
        write(88,*) dtau, l*dtau
c        set up the index vectors
        call indexset()

c        set random number generator

        mu = 0.d0
        expmu=dexp(dtau*mu)
        call setvup()

c       call autoset()

        call ranlat()
c        if(dens .eq. 1.d0)then
c          mu = 0.0
c        else if(dens .lt. 0.0) then
c          mu = dens
c        else
c          mu = zbrent(getden,-5.0,5.0,5e-3)
c        endif
c Just get a mu for now... R.N.
        mu = dens
        expmu=dexp(dtau*mu)
        write (66,*) 'Using mu = ',mu
        write (6,*) 'Using mu = ',mu

        write (66,*) 'Impurity values:'
        do 9000 j = 0, toff-4, 4
          write(66,9999) delta(j), delta(j+1), delta(j+2), delta(j+3)
9999      format(4f12.3)
9000    continue


        nmax = l / 2
        accept = 0
        reject = 0
        accept2 = 0
        reject2 = 0
        do 141 ti = 0, l
         do 141 j1 = 0, n/2
          do 141 j2 = 0, n/2
           bgnl(j1,j2,ti) = 0.d0
           bgnw(j1,j2,ti) = 0.d0
           bgql(j1,j2,ti) = 0.d0
           bgqw(j1,j2,ti) = 0.d0
           bsigma(j1,j2,ti) = 0.d0
           bchinl(j1,j2,ti) = 0.d0
           bchinw(j1,j2,ti) = 0.d0
           bchiql(j1,j2,ti) = 0.d0
           bchiqw(j1,j2,ti) = 0.d0
           bdnl(j1,j2,ti) = 0.d0
           bdnw(j1,j2,ti) = 0.d0
           bdql(j1,j2,ti) = 0.d0
           bdqw(j1,j2,ti) = 0.d0
           bpi(j1,j2,ti) = 0.d0
           sgnl(j1,j2,ti) = 0.d0
           sgnw(j1,j2,ti) = 0.d0
           sgql(j1,j2,ti) = 0.d0
           sgqw(j1,j2,ti) = 0.d0
           ssigma(j1,j2,ti) = 0.d0
           schinl(j1,j2,ti) = 0.d0
           schinw(j1,j2,ti) = 0.d0
           schiql(j1,j2,ti) = 0.d0
           schiqw(j1,j2,ti) = 0.d0
           sdnl(j1,j2,ti) = 0.d0
           sdnw(j1,j2,ti) = 0.d0
           sdql(j1,j2,ti) = 0.d0
           sdqw(j1,j2,ti) = 0.d0
           spi(j1,j2,ti) = 0.d0
141     continue
        call setvup()

c        perform warmup sweeps
        wraps = nwrap
        redo = 0
        noredo = 0
        call getgp(vup,0,gmatup,sgnup,detup)
        do 10 i=1,warms
c         write(6,*)'Starting warmup sweep ',i
          call sweep(gmatup,gmatdn,accept,reject,wraps)
          call sweep2(gmatup,gmatdn,accept2,reject2,
     1             wraps,numtry,gsize)
10      continue
        write(66,*)'after warmups, accept ratio is ',
     1                float(accept)/(accept+reject)
	if (numtry.ne.0) then
        write(66,*)'after warmups, accept2 ratio is ',
     1                float(accept2)/(accept2+reject2)
	endif
        write(66,*)'gamma is ',gam
        write(66,*)'redo ratio is ',
     1                float(redo)/(redo+noredo)

c        write(67,*)10000.d0

c        perform measurement sweeps
        call setvup()
        call zeroas()
       
        do 20 i=1,sweeps
c          write(6,*)'Starting measurement sweep ',i
c          if(mod(i,10) .eq. 0) then
c            write(6,*)'accept, redo ratios are ',
c     1                float(accept)/(accept+reject),
c     2                float(redo)/(noredo+redo)
c          endif
          call sweep(gmatup,gmatdn,accept,reject,wraps)
          call sweep2(gmatup,gmatdn,accept2,reject2,
     1           wraps,numtry,gsize)


          if(mod(i,phonskip) .eq. 0) then
            call measphon
          endif
          if(mod(i,tausk) .eq. 0)then           
            call meastau
          endif
          if(mod(i,sweeps/nbin) .eq. 0) then
c            write(66,*)'Finished measurement sweep ',i
            if(dopair .eq. 0)nmeasp = 1
            write(66,9045) asgn/nmeas0,asgnp/nmeasp,
     1            float(accept)/(accept+reject),
     2            float(redo)/(noredo+redo)
9045        format('asgn, asgnp: ',2(f8.3,' '),
     1          ';accept,redo ratios: ',2('  ',f9.4))
            k = (i*nbin)/sweeps
c In the following lines, we have replaced nmeas0 with asgn
c for a different error estimate.
            if(asgn .eq. 0.d0) then
              write(66,*)'adding .01 to asgn'
              asgn = asgn + .01
            endif
            bnup(k)= anup / asgn
            bndn(k)= andn / asgn
            bntot(k)= (anup+andn) / asgn
            bxeq(k)= axeq / asgn
            bsaf(k)= asaf / asgn
            bsafsq(k)= dsqrt(abs(asafsq / asgn))
            bsferro(k)= asferro / asgn
            bsfer2(k)= asfer2 / asgn
            bscdw(k)= ascdw / asgn
            bke(k)= ake / asgn
            benergy(k)= (ake+aphpe+aphke+apeep) / asgn
            bnud(k)= anud / asgn
            bphpe(k)=aphpe/ asgn
            bphke(k)=aphke/ asgn
            bpeep(k)=apeep/ asgn
            bphnn(k)=aphnn/ asgn
c
            bsgnup(k) = asgnup / nmeas0
            bsgndn(k) = asgndn / nmeas0
            bsgn(k) = asgn / nmeas0
            do 41 j1 = 0, n/2
              do 41 j2 = 0, n/2
                bgrfun(k,j1,j2) = agrfun(j1,j2) / asgn
                bden(k,j1,j2,0) = aden(j1,j2,0) / asgn
                bden(k,j1,j2,1) = aden(j1,j2,1) / asgn
                bspinxx(k,j1,j2) = aspinxx(j1,j2) / asgn
                tcdw(j1,j2) = acdw(j1,j2) / asgn
                bcdw(k,j1,j2) = acdw(j1,j2) / asgn
                bchicdw(k)=achicdw/asgnt
c
c     Set gnl and chinl for this tenth of the run.
c

                do 41 ti = 0, l
                  gnl(j1,j2,ti) = agnl(j1,j2,ti)/asgnt
                  gql(j1,j2,ti) = 0.d0
                  chinl(j1,j2,ti) = achinl(j1,j2,ti)/asgnt
                  chiql(j1,j2,ti) = 0.d0
                  dnl(j1,j2,ti) = adnl(j1,j2,ti)/asgnph
                  dql(j1,j2,ti) = 0.d0
41          continue
            call ftntok(gnl,gql,n,l)
            call ftntok(chinl,chiql,n,l)
            call ftltow(gnl,gnw,n,l,dtau,0,nmax)
            call ftltow(gql,gqw,n,l,dtau,0,nmax)
            call ftltow(chinl,chinw,n,l,dtau,1,nmax)
            call ftltow(chiql,chiqw,n,l,dtau,1,nmax)
            call ft2ntok(tcdw,cdwq,n)
            call ftntok(dnl,dql,n,l)
c     The 2nd to last argument to ftltow means use Bose frequencies.
c            
            call ftltow(dnl,dnw,n,l,dtau,1,nmax)
            call ftltow(dql,dqw,n,l,dtau,1,nmax)

            do 586 mx = 0, n/2
             do 586 my = 0, n/2
              bcdwq(k,mx,my) = cdwq(mx,my)
              do 586 ti = 0, l
               bgnl(mx,my,ti) = bgnl(mx,my,ti) + gnl(mx,my,ti)
               bgql(mx,my,ti) = bgql(mx,my,ti) + gql(mx,my,ti)
               bchinl(mx,my,ti) = bchinl(mx,my,ti)+chinl(mx,my,ti)
               bchiql(mx,my,ti) = bchiql(mx,my,ti)+chiql(mx,my,ti)
               bdnl(mx,my,ti) = bdnl(mx,my,ti) + dnl(mx,my,ti)
               bdql(mx,my,ti) = bdql(mx,my,ti) + dql(mx,my,ti)
               sgnl(mx,my,ti) = sgnl(mx,my,ti) + gnl(mx,my,ti)**2
               sgql(mx,my,ti) = sgql(mx,my,ti) + gql(mx,my,ti)**2
               schinl(mx,my,ti)=schinl(mx,my,ti)+chinl(mx,my,ti)**2
               schiql(mx,my,ti)=schiql(mx,my,ti)+chiql(mx,my,ti)**2
               sdnl(mx,my,ti) = sdnl(mx,my,ti) + dnl(mx,my,ti)**2
               sdql(mx,my,ti) = sdql(mx,my,ti) + dql(mx,my,ti)**2

              if(ti .le. nmax)then

               bgnw(mx,my,ti) = bgnw(mx,my,ti) + gnw(mx,my,ti)
               bgqw(mx,my,ti) = bgqw(mx,my,ti) + gqw(mx,my,ti)
               bchinw(mx,my,ti) = bchinw(mx,my,ti)+chinw(mx,my,ti)
               bchiqw(mx,my,ti) = bchiqw(mx,my,ti)+chiqw(mx,my,ti)
               bdnw(mx,my,ti) = bdnw(mx,my,ti) + dnw(mx,my,ti)
               bdqw(mx,my,ti) = bdqw(mx,my,ti) + dqw(mx,my,ti)
               sgnw(mx,my,ti) = sgnw(mx,my,ti) + dcmplx(dreal
     1                 (gnw(mx,my,ti))**2,dimag(gnw(mx,my,ti))**2)
               sgqw(mx,my,ti) = sgqw(mx,my,ti) + dcmplx(dreal
     1                 (gqw(mx,my,ti))**2,dimag(gqw(mx,my,ti))**2)
               schinw(mx,my,ti) = schinw(mx,my,ti) + dcmplx(dreal
     1              (chinw(mx,my,ti))**2,dimag(chinw(mx,my,ti))**2)
               schiqw(mx,my,ti) = schiqw(mx,my,ti) + dcmplx(dreal
     1              (chiqw(mx,my,ti))**2,dimag(chiqw(mx,my,ti))**2)
               epsk = -2.d0*t*(cos(twopi*mx/n)+cos(twopi*my/n))
               sigma(mx,my,ti)=dcmplx(-epsk+mu,(ti+0.5)*twopi/l/dtau)
     1                             + 1.d0/gqw(mx,my,ti)
               bsigma(mx,my,ti) = bsigma(mx,my,ti) + sigma(mx,my,ti)
               ssigma(mx,my,ti) = ssigma(mx,my,ti) + dcmplx(dreal
     1              (sigma(mx,my,ti))**2,dimag(sigma(mx,my,ti))**2)
               sdnw(mx,my,ti) = sdnw(mx,my,ti) + dcmplx(dreal
     1                 (dnw(mx,my,ti))**2,dimag(dnw(mx,my,ti))**2)
               sdqw(mx,my,ti) = sdqw(mx,my,ti) + dcmplx(dreal
     1                 (dqw(mx,my,ti))**2,dimag(dqw(mx,my,ti))**2)
              endif

586         continue

            if(dopair .eq. 1) then

              if(asgnp .eq. 0.d0) then
                write(66,*)'adding .01 to asgnp'
                asgnp = asgnp + .01
              endif
              
            do 8379 mx = -1, 1
            do 8379 my = -1, 1
             do 8379 mpx = -1, 1
             do 8379 mpy = -1, 1
              do 8356 ti = 0,l
8356           pairgsus(mx,my,mpx,mpy,ti) = 0.d0
              do 8379 lx = 0, n-1
              do 8379 ly = 0, n-1
               kx = min(lx,n-lx)
               ky = min(ly,n-ly)
               kmx = mod(n+n+lx-(mx-mpx),n)
               kmx = min(kmx,n-kmx)
               kmy = mod(n+n+ly-(my-mpy),n)
               kmy = min(kmy,n-kmy)
               do 8379 ti = 0,l
              pairgsus(mx,my,mpx,mpy,ti)=pairgsus(mx,my,mpx,mpy,ti)
     1          +gnl(kmx,kmy,ti)*gnl(kx,ky,ti)
8379        continue

            do 2378 j1 = 1,9
             bpairv(k,j1) = 0.d0
             do 2338 ti = 0, l
               gpairl(k,j1,ti) = 0.d0
2338           bpairl(k,j1,ti) = 0.d0
             do 2378 mx = -1, 1
              do 2378 my = -1, 1
               do 2378 mpx = -1, 1
                do 2378 mpy = -1, 1
                  bpairv(k,j1)=bpairv(k,j1)+waves(mx,my,j1)*
     1            apairmat(mx,my,mpx,mpy)*waves(mpx,mpy,j1)/asgnp
                 do 2378 ti = 0, l
                bpairl(k,j1,ti)=bpairl(k,j1,ti)+waves(mx,my,j1)*
     1           apairsus(mx,my,mpx,mpy,ti)*waves(mpx,mpy,j1)/asgnt
                gpairl(k,j1,ti)=gpairl(k,j1,ti)+waves(mx,my,j1)*
     1            pairgsus(mx,my,mpx,mpy,ti)*waves(mpx,mpy,j1)
2378          continue
              call ftltowp(bpairl,bpairw,l,k,dtau,l)
              call ftltowp(gpairl,gpairw,l,k,dtau,l)
            endif
            call zeroas()

          endif
         

20      continue

        write(66,*)'At end, redo ratio is ',float(redo)/(redo+noredo)
        write(66,*)'gamma is ',gam
       write(66,*)'Acceptance ratio = ',float(accept)/(accept+reject)
	if (numtry.ne.0) then
       write(66,*)'Accept2 ratio = ',float(accept2)/(accept2+reject2)
	endif
        call geterr(bsgnup,aveval,errval)
        write(66,*) 'Average up sign =',aveval,' +- ',errval
        call geterr(bsgndn,aveval,errval)
        write(66,*) 'Average dn sign =',aveval,' +- ',errval
        call geterr(bsgn,aveval,errval)
        write(66,*) 'Average total sign =',aveval,' +- ',errval
        call geterr(bntot,aveval,errval)
        write(66,*) 'Average density =',
     1         aveval,' +- ',errval
        call geterr(bnup,aveval,errval)
        write(66,*) 'Average up occupancy =',
     1         aveval,' +- ',errval
        call geterr(bndn,aveval,errval)
        write(66,*) 'Average dn occupancy =',
     1         aveval,' +- ',errval
        call geterr(bxeq,aveval,errval)
        write(66,*) 'Average phonon displacement  =',
     1         aveval,' +- ',errval
        write(88,*) 'Average displacement'
          write(88,*) aveval*sqrt(2.0*omega)
        call geterr(benergy,aveval,errval)
        write(66,*) 'Average Energy =',
     1         aveval,' +- ',errval
        call geterr(bke,aveval,errval)
        write(66,*) 'Average Kinetic Energy =',
     1         aveval,' +- ',errval
        call geterr(bpeep,aveval,errval)
        write(66,*) 'Average EL-PH PE =',
     1         aveval,' +- ',errval
        call geterr(bnud,aveval,errval)
        write(66,*) 'Average Nup*Ndn =',
     1         aveval,' +- ',errval
        call geterr(bphpe,aveval,errval)
        write(66,*) 'Average Phonon PE =',
     1         aveval,' +- ',errval
        call geterr(bphke,aveval,errval)
        write(66,*) 'Average Phonon KE =',
     1         aveval,' +- ',errval
        call geterr(bphnn,aveval,errval)
        write(66,*) 'Average Phonon NN=',
     1         aveval,' +- ',errval
        call geterr(bsaf,aveval,errval)
        write(66,*) 'AF correlation function (xx) = ',
     1         aveval,' +- ',errval
        call geterr(bscdw,aveval,errval)
        write(66,*) ' CDW correlation function= ',
     1         aveval,' +- ',errval
          call geterr(bchicdw,aveval,errval)
          write (66,*) ' CDW susceptibility= ',
     1           aveval,' +/- ',errval
        call geterr(bsferro,aveval,errval)
        write(66,*) 'Ferro correlation function(xx)= ',
     1         aveval,' +- ',errval
        call geterr(bsfer2,aveval,errval)
        write(66,*) 'Ferro correlation function(zz)= ',
     1         aveval,' +- ',errval
        write(66,*)'Green''s function:'
        do 677 j1 = 0, n/2
          do 677 j2 = 0, n/2
            call geterr(bgrfun(1,j1,j2),aveval,errval)
            if(j2 .ge. j1)
     1        write(66,*)j1,j2,aveval,' +- ',errval
            grfun(j1,j2) = aveval
677     continue

        write (66,*) 'SRW CODE'
        write(66,*)'density-density correlation fn: (up-up,up-dn)'
        do 977 j1 = 0, n/2
          do 977 j2 = j1, n/2
            call geterr(bden(1,j1,j2,0),aveval,errval)
            call geterr(bden(1,j1,j2,1),aveval2,errval2)
            write(66,1987)j1,j2,aveval,errval,aveval2,errval2
1987        format(2i4,2('    ',f12.6,' +- ',f12.6))
977     continue

        write(66,*)'density-density correlation function:'
        do 687 j1 = 0, n/2
          do 687 j2 = j1, n/2
            call geterr(bcdw(1,j1,j2),aveval,errval)
            write(66,*)j1,j2,aveval,' +- ',errval
687     continue
c
c     Added 4/27/90 R.N.
c
        write(66,*) 'Scdw(q): '
        do 1576 j1 = 0, n/2
          do 1576 j2 = j1, n/2
            call geterr(bcdwq(1,j1,j2),aveval,errval)
            write(66,*)j1,j2,aveval,' +- ',errval
1576    continue

        write(66,*)'xx Spin correlation function:'
        do 697 j1 = 0, n/2
          do 697 j2 = j1, n/2
            call geterr(bspinxx(1,j1,j2),aveval,errval)
            write(66,*)j1,j2,aveval,' +- ',errval
697     continue
        call geterr(bsafsq,aveval,errval)
        write(66,*) 'RMS AF correlation function (xx) = ',
     1         aveval,' +- ',errval
        write(66,*)' '

        write(68,*)'G(nx,ny,ti):'
1234      format(i5,f14.6,' +- ',f14.6)
        do 987 j1 = 0, n/2
         do 987 j2 = j1, n/2
          write(68,*)'nx = ',j1,' ny = ',j2
          do 987 ti = 0,l
          call geterr2(bgnl(j1,j2,ti),sgnl(j1,j2,ti),aveval,errval)
          write(68,1234)ti,-aveval,errval
987     continue
c
c Write G(q,ti) and <nk> -- <nk> in the results file.
c
        write(68,*)'G(qx,qy,ti):'
        write(66,*)'<nk>:'
        do 9871 j1 = 0, n/2
         do 9871 j2 = j1, n/2
          write(68,*)'qx = ',j1,' qy = ',j2
          do 9871 ti = 0,l
          call geterr2(bgql(j1,j2,ti),sgql(j1,j2,ti),aveval,errval)
          write(68,1234)ti,-aveval,errval
          if (ti.eq.0) then
            write(66,9898) j1, j2, (1.d0 - aveval),errval
          endif
9871     continue
9898     format (2i5,2f14.4)

         write(68,*)'G(nx,ny,omega), omega = (n+.5) 2 pi T :'
         do 9872 j1 = 0, n/2
           do 9872 j2 = j1, n/2
             write(68,*)'nx = ',j1,' ny = ',j2
             do 9872 ti = 0,nmax
               call geterr2(dreal(bgnw(j1,j2,ti)),dreal(sgnw(j1,j2,ti)),
     1              aveval,errval)
               call geterr2(dimag(bgnw(j1,j2,ti)),dimag(sgnw(j1,j2,ti)),
     1              aveval2,errval2)
               write(68,1235)ti,-aveval,errval,-aveval2,errval2
1235           format(i5,'(',f12.6,' +- ',f12.6,') + i * ('
     1              ,f12.6,' +- ',f12.6,')')
9872     continue

        write(68,*)'G(qx,qy,omega):'
        do 9873 j1 = 0, n/2
         do 9873 j2 = j1, n/2
          write(68,*)'qx = ',j1,' qy = ',j2
          do 9873 ti = 0,nmax
          call geterr2(dreal(bgqw(j1,j2,ti)),dreal(sgqw(j1,j2,ti)),
     1                      aveval,errval)
          call geterr2(dimag(bgqw(j1,j2,ti)),dimag(sgqw(j1,j2,ti)),
     1                      aveval2,errval2)
          write(68,1235)ti,-aveval,errval,-aveval2,errval2
9873     continue

        write(68,*)'SIGMA(qx,qy,omega):'
        do 9874 j1 = 0, n/2
         do 9874 j2 = j1, n/2
          write(68,*)'qx = ',j1,' qy = ',j2
          do 9875 ti = nmax,0,-1
           call geterr2(dreal(bsigma(j1,j2,ti)),dreal(ssigma(j1,j2,ti)),
     1                      aveval,errval)
         call geterr2(dimag(bsigma(j1,j2,ti)),dimag(ssigma(j1,j2,ti))
     1                      ,aveval2,errval2)
           write(68,1235)-(2*ti+1), aveval,errval,-aveval2,errval2
9875      continue
          do 9874 ti = 0,nmax
           call geterr2(dreal(bsigma(j1,j2,ti)),dreal(ssigma(j1,j2,ti)),
     1                      aveval,errval)
         call geterr2(dimag(bsigma(j1,j2,ti)),dimag(ssigma(j1,j2,ti))
     1                      ,aveval2,errval2)
           write(68,1235)(2*ti+1), aveval,errval,aveval2,errval2
9874     continue

        write(68,*)'chi(nx,ny,ti):'
        do 417 j1 = 0, n/2
         do 417 j2 = j1, n/2
          write(68,*)'nx = ',j1,' ny = ',j2
          do 417 ti = 0,l
        call geterr2(bchinl(j1,j2,ti),schinl(j1,j2,ti),aveval,errval)
          write(68,1234)ti, aveval,errval
417     continue

        write(68,*)'chi(qx,qy,ti):'
        do 4171 j1 = 0, n/2
         do 4171 j2 = j1, n/2
          write(68,*)'qx = ',j1,' qy = ',j2
          do 4171 ti = 0,l
        call geterr2(bchiql(j1,j2,ti),schiql(j1,j2,ti),aveval,errval)
          write(68,1234)ti, aveval,errval
4171     continue

        write(68,*)'chi(nx,ny,omega), omega = 2 n pi T :'
        do 4172 j1 = 0, n/2
         do 4172 j2 = j1, n/2
          write(68,*)'nx = ',j1,' ny = ',j2
          do 4172 ti = 0,nmax
       call geterr2(dreal(bchinw(j1,j2,ti)),dreal(schinw(j1,j2,ti)),
     1                      aveval,errval)
       call geterr2(dimag(bchinw(j1,j2,ti)),dimag(schinw(j1,j2,ti)),
     1                      aveval2,errval2)
          write(68,1235)ti, aveval,errval,aveval2,errval2
4172     continue

        write(68,*)'chi(qx,qy,omega):'
        do 4173 j1 = 0, n/2
         do 4173 j2 = j1, n/2
          write(68,*)'qx = ',j1,' qy = ',j2
          do 4173 ti = 0,nmax
       call geterr2(dreal(bchiqw(j1,j2,ti)),dreal(schiqw(j1,j2,ti)),
     1                      aveval,errval)
       call geterr2(dimag(bchiqw(j1,j2,ti)),dimag(schiqw(j1,j2,ti)),
     1                      aveval2,errval2)
          write(68,1235)ti, aveval,errval,aveval2,errval2
4173     continue

         write(68,*)'D(qx,qy,ti):'
         do 4609 j1 = 0, n/2
           do 4609 j2 = j1, n/2
             write(68,*)'qx = ',j1,' qy = ',j2
             do 4609 ti = 0,l
               call geterr2(bdql(j1,j2,ti),sdql(j1,j2,ti),
     @              aveval,errval)
               write(68,1234)ti,-aveval,errval
4609      continue

         write(68,*)'D(nx,ny,ti):'
         do 5609 j1 = 0, n/2
           do 5609 j2 = j1, n/2
             write(68,*)'nx = ',j1,' ny = ',j2
             do 5609 ti = 0,l
               call geterr2(bdnl(j1,j2,ti),sdnl(j1,j2,ti),
     @              aveval,errval)
               write(68,1234)ti,-aveval,errval
5609      continue

         write(68,*)'D(nx,ny,omega), omega = n 2 pi T :'
         do 6665 j1 = 0, n/2
           do 6665 j2 = j1, n/2
             write(68,*)'nx = ',j1,' ny = ',j2
             do 6665 ti = 0,nmax
               call geterr2(dreal(bdnw(j1,j2,ti)),dreal(sdnw(j1,j2,ti)),
     1              aveval,errval)
               call geterr2(dimag(bdnw(j1,j2,ti)),dimag(sdnw(j1,j2,ti)),
     1              aveval2,errval2)
               write(68,1235)ti,-aveval,errval,-aveval2,errval2
6665     continue

        write(68,*)'D(qx,qy,omega):'
        do 4921 j1 = 0, n/2
         do 4921 j2 = j1, n/2
          write(68,*)'qx = ',j1,' qy = ',j2
          do 4921 ti = 0,nmax
          call geterr2(dreal(bdqw(j1,j2,ti)),dreal(sdqw(j1,j2,ti)),
     1                      aveval,errval)
          call geterr2(dimag(bdqw(j1,j2,ti)),dimag(sdqw(j1,j2,ti)),
     1                      aveval2,errval2)
          write(68,1235)ti,-aveval,errval,-aveval2,errval2
4921     continue

       if(dopair .eq. 1) then
1236    format(2(f12.5,' +- ',f12.5,'   '))
        write(66,*) 's-wave: (corr. fn, no vertex)'
        write(66,*) '        (suscept., no vertex)'
        call geterr(bpairv(1,1),aveval,errval)
        call geterr(gpairl(1,1,0),aveval2,errval2)
        write(66,1236)aveval,errval,aveval2,errval2
        call geterrc(bpairw(1,1,0),cveval,crrval)
        call geterrc(gpairw(1,1,0),cveval2,crrval2)
        write(66,1236)dreal(cveval),dreal(crrval),dreal(cveval2)
     1               ,dreal(crrval2)



       endif

        saf = 2.d0*grfun(0,0)
        do 3879 lx = 0, n-1
          do 3879 ly = 0, n-1
            kx = min(lx,n-lx)
            ky = min(ly,n-ly)
3879          saf = saf -(-1)**(lx+ly)*2*grfun(kx,ky)*grfun(kx,ky)
        write(66,*)'saf with no vertex is ',saf
          do 777 i=0,volume-1
         write (60,*) i,hub(i)
777        continue


        stop
        end
c
c init.f -- initialization routines for 2d Holstein model.
c
c********* phaseset() - set and load phase vectors **********
        subroutine phaseset()

c pars.h -- basic Holstein parameters
        use global
        implicit none
c vectors.h -- vectors for 2d holstein model.
        real*8 hub(0:volume-1)
        real*8 vup(0:volume-1),vdn(0:volume-1)
        real*8 mphase(0:toff)

        common/vectors/hub,vup,vdn, mphase

        integer ix,iy,addr
        real*8 ph

c        Calculate mphase - the phase vector needed in the measurement process
c        mphase = +1 on odd spatial sites
c        mphase = -1 on even spatial sites

        addr=0
        ph=1.d0
        do 30 iy=0,n-1
        do 25 ix=0,n-1
        mphase(addr)=ph
        addr=addr+1
25        ph=-ph
30        ph=-ph
        mphase(toff)=1.d0

        return
        end

c**************************indexset()************************
c        this subroutine sets and loads the index vectors
        subroutine indexset()

c pars.h -- basic Holstein parameters
        use global
        implicit none
c ivectors.h -- index vectors for indirect addressing and phase vector
        integer xplus(0:toff-1),xminus(0:toff-1)
        integer yplus(0:toff-1),yminus(0:toff-1)

        common/ivectors/xplus,xminus,yplus,yminus

        integer i,neighbor

c        calculate the index vectors
        do 10 i=0,toff-1
            xplus(i) =neighbor(i,1,0,0)
            xminus(i)=neighbor(i,-1,0,0)
            yplus(i) =neighbor(i,0,1,0)
            yminus(i)=neighbor(i,0,-1,0)
10        continue

        return
        end

c*******************************************
        subroutine readin()

c pars.h -- basic Holstein parameters
        use global
        implicit none
c integers.h -- integer parameters for holstein.
        integer warms,sweeps,msr,nwrap
        integer iran
        common/integers/warms,sweeps,msr,nwrap, iran
c
c getgpam.h -- parameters for getg (2dmoo)
c
        integer orthlen,doauto,torth
        real*8 eorth,difflim,errpam
        common/getgpam/eorth,difflim,errpam,orthlen,doauto,
     1                torth
c m0var.h -- common for measurements. (2dmoo)
        real*8 sgnup,sgndn,sgn,nup,ndn,nud,ke
        real*8 saf,scdw,sferro,sfer2,grfun(0:n/2,0:n/2),safsq,safsq2
        real*8 spinxx(0:n/2,0:n/2),cdw(0:n/2,0:n/2)
        real*8 pairmat(-1:1,-1:1,-1:1,-1:1),errrat,asgnp
        integer nmeas0,redo,noredo,dopair,nmeasp,numpair,dowrite
        real*8 asgnup,asgndn,asgn,anup,andn,anud,ake,asafsq,asafsq2
        real*8 asaf,ascdw,asferro,asfer2,agrfun(0:n/2,0:n/2)
        real*8 aspinxx(0:n/2,0:n/2),acdw(0:n/2,0:n/2)
        real*8 apairmat(-1:1,-1:1,-1:1,-1:1)
        real*8 den(0:n/2,0:n/2,0:1),aden(0:n/2,0:n/2,0:1)
        common/m0var/ sgnup,sgndn,sgn,nup,ndn,nud,ke,safsq,safsq2,
     1       saf,scdw,sferro,sfer2,grfun,
     1       spinxx,cdw,aspinxx,acdw,asafsq,asafsq2,
     2       asgnup,asgndn,asgn,anup,andn,anud,ake,den,aden,
     3       asaf,ascdw,asferro,asfer2,agrfun,errrat,asgnp
     4       ,pairmat,apairmat
        common/m0vari/nmeas0,redo,noredo,dopair,nmeasp,numpair,dowrite

c couple.h -- coupling parameters for holstein.
         real*8 move,t,omega,g,mu,dtau,expmu,gam,lambda,tdtau,dens
         real*8 pe,kee
         common/couple/move,t,omega,g,mu,
     1       dtau,expmu,gam,lambda,tdtau,dens,pe,kee
c kinetic.h -- kinetic energy common block. (2dmoo)
        real*8 ch,soc,ch4
        common/kinetic/ch,soc,ch4
c names.h -- common for I/O names.
        character*80 outname, inlat
        common/names/outname, inlat


c impure.h -- common block for impurity.
      real*8 delta0, delta(0:toff-1), expdel(0:toff-1)

      common/impure/ delta0, delta, expdel


        write(6,*) 'enter warms and sweeps'
        read(5,*) warms,sweeps


        write(6,*)'enter t,omega,g,move,dens '
        read(5,*) t,omega,g,move,dens


        write(6,*) 'enter dtau, nwrap, difflim, errrat'
        read(5,*)dtau,nwrap,difflim,errrat

        write(6,*) 'enter doauto,orthlen,eorth,dopair,numpair,dowrite'
        read(5,*)doauto,orthlen,eorth,dopair, numpair, dowrite
         doauto=0

        write(6,*) 'enter torth'
        read(5,*)torth

        gam = 0.5
        errpam = 1.0e13

        tdtau=t*dtau
c
c Note we have defined lambda so that g and omega are as in the
c Hamiltonian with interaction term:
c                        +
c                +g(b + b )(nup + ndn)
c                       
        lambda = g*dsqrt(2.d0*omega)*dtau

        pe=0.5*omega*omega*dtau

        kee=0.5/dtau
        ch = cosh(tdtau)
        soc = sinh(tdtau)/ch
        ch4 = ch**4

        write(6,*)'enter random number seed'
        read(5,*)iran

        write(6,*) 'Enter impurity scale'
        read(5,*) delta0


        write(6,*) 'Enter output file string'
        read(5,2010) outname
        write(6,*) 'Enter input lattice name'
        read(5,2010) inlat
2010    format(A40)
        return
        end

c*****************************ranlat()*****************************

c        This subroutine makes a lattice with spatially
c        random values for the phonon variables.  The phonon field is,
c        however, temporally ordered to start.

        subroutine ranlat()

c pars.h -- basic Holstein parameters
        use global
        implicit none
c integers.h -- integer parameters for holstein.
        integer warms,sweeps,msr,nwrap
        integer iran
        common/integers/warms,sweeps,msr,nwrap, iran
c couple.h -- coupling parameters for holstein.
         real*8 move,t,omega,g,mu,dtau,expmu,gam,lambda,tdtau,dens
         real*8 pe,kee
         common/couple/move,t,omega,g,mu,
     1       dtau,expmu,gam,lambda,tdtau,dens,pe,kee
c vectors.h -- vectors for 2d holstein model.
        real*8 hub(0:volume-1)
        real*8 vup(0:volume-1),vdn(0:volume-1)
        real*8 mphase(0:toff)

        common/vectors/hub,vup,vdn, mphase
c names.h -- common for I/O names.
        character*80 outname, inlat
        common/names/outname, inlat

c impure.h -- common block for impurity.
      real*8 delta0, delta(0:toff-1), expdel(0:toff-1)

      common/impure/ delta0, delta, expdel


        integer istart

        integer i,j,Lp,size
        real*8 scale
        real*8 ran
        real*8 ran2
        
        write (6,*) '1 for random start, 2 for read in start'
        read (5,*) istart
        write (6,*) 'INPUT size of stored lattice in L space'
        write (6,*)  'should be no smaller than  ',L/2
        write (6,*) 'spatial size should be  ',n
        read (5,*) Lp
        if (istart.eq.1) then
          scale=1./omega/(dexp(L*dtau*omega)-1.)
     1         +.5d0/omega
          scale=2.d0*dsqrt(scale)
          write (66,*) 'initial phonon scale is',scale
          do 10 i=0,toff-1

            ran = ran2(iran)
            hub(i) = scale*(ran -0.5d0)

            do 5 j=1,L-1
              hub(i+j*toff)=hub(i)
5           continue
10        continue
        else if (istart.eq.2) then

          open(unit=61,status='old',file=inlat)


          size=Lp*toff
          do 105 i=0,size-1
            read(61,*) j,hub(i)
105       continue
          do 106 i=size,volume-1
            hub(i)=hub(i-size)
106       continue
        else
          write (6,*) 'BAD START INSTRUCTION'
          stop
        endif
        
        call setvup()
        

        do 300 i=0, toff-1

            ran = ran2(iran)
            delta(i) = delta0*(ran -0.5d0)

            expdel(i) = dexp(delta(i)*dtau)
300     continue


        return
        end

c*****************************setvup()*****************************

c        This subroutine sets vup and vdn given hubs.

        subroutine setvup()

c pars.h -- basic Holstein parameters
        use global
        implicit none
c integers.h -- integer parameters for holstein.
        integer warms,sweeps,msr,nwrap
        integer iran
        common/integers/warms,sweeps,msr,nwrap, iran
c couple.h -- coupling parameters for holstein.
         real*8 move,t,omega,g,mu,dtau,expmu,gam,lambda,tdtau,dens
         real*8 pe,kee
         common/couple/move,t,omega,g,mu,
     1       dtau,expmu,gam,lambda,tdtau,dens,pe,kee
c vectors.h -- vectors for 2d holstein model.
        real*8 hub(0:volume-1)
        real*8 vup(0:volume-1),vdn(0:volume-1)
        real*8 mphase(0:toff)

        common/vectors/hub,vup,vdn, mphase

        integer i

        do 20 i=0,volume-1
          vup(i)=dexp(lambda*hub(i))
          vdn(i)=dexp(lambda*hub(i))
20      continue
        call addaf()

        return
        end
c
c util.f -- General Utility routines for holstein model.
c
c****************************siteindx**********************************

c        this function finds the index of a site given its three coordinates

        integer function siteindx(x,y,ti)

c pars.h -- basic Holstein parameters
        use global
        implicit none
        integer x,y,ti

        siteindx = x + n*(y + n*(ti))

        return
        end

c****************************neighbor()*****************************

c        this function finds the index of a neighboring site

        integer function neighbor(site,delx,dely,delt)

c pars.h -- basic Holstein parameters
        use global
        implicit none

        integer site,delx,dely,delt
        integer x,y,ti,siteindx

c        find the coordinates of site
        x = mod(site,n)
        y = mod(site/n,n)
        ti = mod(site/(n*n),l)

c        find the coordinates of the neighbor
        x = mod(x+delx+n,n)
        y = mod(y+dely+n,n)
        ti = mod(ti+delt+l,l)

c        find the index of the neighboring site
        neighbor = siteindx(x,y,ti)

        return
        end

c************************GETERR*****************************
        subroutine geterr(bin,ave,err)
        use global
        real*8 bin(nbin),err,sum,sumsq,ave
        integer i

        sum = 0.d0
        do 10 i = 1, nbin
10        sum = sum + bin(i)
        ave = sum /dfloat(nbin)
        sumsq = 0.d0
        do 20 i = 1, nbin
20          sumsq = sumsq + (bin(i)-ave)**2
        sumsq =  sumsq/dfloat(nbin)
        err = dsqrt(sumsq)/dsqrt(dfloat(nbin-1))
        if(err*1.0e12 .lt. abs(ave)) err = 0.d0
        return
        end

c************************GETERRC*****************************
        subroutine geterrc(bin,ave,err)
        use global
        complex*16 bin(nbin),err,ave
        real*8 rsumsq,isumsq
        integer i

        ave = 0.d0
        do 10 i = 1, nbin
10        ave = ave + bin(i)
        ave = ave /dfloat(nbin)
        rsumsq = 0.d0
        do 20 i = 1, 10
          rsumsq = rsumsq + (dreal(bin(i)-ave))**2
20        isumsq = isumsq + (dimag(bin(i)-ave))**2
        rsumsq =  rsumsq /dfloat(nbin)
        isumsq =  isumsq /dfloat(nbin)
        err = dcmplx( dsqrt(rsumsq)/dsqrt(dfloat(nbin-1)),
     1                dsqrt(isumsq)/dsqrt(dfloat(nbin-1)) )
        return
        end

c************************GETERR2*****************************
        subroutine geterr2(sum,sumsq,ave,err)
        use global

        real*8 err,sum,sumsq,ave,ssq

        ave = sum /dfloat(nbin)
        ssq = sumsq/dfloat(nbin) - ave**2
        err = dsqrt(abs(ssq))/dsqrt(dfloat(nbin-1))
        if(err*1.0e12 .lt. abs(ave)) err = 0.d0
        return
        end

ccccccccccccccccccccc Subroutine addaf() cccccccc
c        This subroutine adds in the small antiferromagnetic field to
c        vup and vdn.

        subroutine addaf()

c pars.h -- basic Holstein parameters
        use global
        implicit none
c couple.h -- coupling parameters for holstein.
         real*8 move,t,omega,g,mu,dtau,expmu,gam,lambda,tdtau,dens
         real*8 pe,kee
         common/couple/move,t,omega,g,mu,
     1       dtau,expmu,gam,lambda,tdtau,dens,pe,kee
c vectors.h -- vectors for 2d holstein model.
        real*8 hub(0:volume-1)
        real*8 vup(0:volume-1),vdn(0:volume-1)
        real*8 mphase(0:toff)

        common/vectors/hub,vup,vdn, mphase

        integer x,y,ti,siteindx,s,i
        real*8 afeps,epsfacup,epsfacdn
        parameter (afeps=0.d0)
c
        if(afeps .ne. 0.d0) then
          do 10 x = 0, n-1
            do 10 y = 0, n-1
c            s = (-1)**(x+y)
              i = siteindx(x,y,0)
              s = mphase(i)
              epsfacup = dexp(-dtau * afeps * s)
              epsfacdn = 1.d0 / epsfacup
              do 10 ti = 0, l-1
                i = siteindx(x,y,ti)
                vup(i)=vup(i) * epsfacup
                vdn(i)=vdn(i) * epsfacdn
10        continue
        endif
        return
        end
c***************************************
        subroutine zeroas()

c pars.h -- basic Holstein parameters
        use global
        implicit none
      
c m0var.h -- common for measurements. (2dmoo)
        real*8 sgnup,sgndn,sgn,nup,ndn,nud,ke
        real*8 saf,scdw,sferro,sfer2,grfun(0:n/2,0:n/2),safsq,safsq2
        real*8 spinxx(0:n/2,0:n/2),cdw(0:n/2,0:n/2)
        real*8 pairmat(-1:1,-1:1,-1:1,-1:1),errrat,asgnp
        integer nmeas0,redo,noredo,dopair,nmeasp,numpair,dowrite
        real*8 asgnup,asgndn,asgn,anup,andn,anud,ake,asafsq,asafsq2
        real*8 asaf,ascdw,asferro,asfer2,agrfun(0:n/2,0:n/2)
        real*8 aspinxx(0:n/2,0:n/2),acdw(0:n/2,0:n/2)
        real*8 apairmat(-1:1,-1:1,-1:1,-1:1)
        real*8 den(0:n/2,0:n/2,0:1),aden(0:n/2,0:n/2,0:1)
        common/m0var/ sgnup,sgndn,sgn,nup,ndn,nud,ke,safsq,safsq2,
     1       saf,scdw,sferro,sfer2,grfun,
     1       spinxx,cdw,aspinxx,acdw,asafsq,asafsq2,
     2       asgnup,asgndn,asgn,anup,andn,anud,ake,den,aden,
     3       asaf,ascdw,asferro,asfer2,agrfun,errrat,asgnp
     4       ,pairmat,apairmat
        common/m0vari/nmeas0,redo,noredo,dopair,nmeasp,numpair,dowrite

c mtauvar.h -- common for tau-dependent measurements. (2dmoo)
        real*8 gnl(0:n/2,0:n/2,0:l),agnl(0:n/2,0:n/2,0:l)
        real*8 chinl(0:n/2,0:n/2,0:l),achinl(0:n/2,0:n/2,0:l)
        real*8 pairsus(-1:1,-1:1,-1:1,-1:1,0:l)
        real*8 apairsus(-1:1,-1:1,-1:1,-1:1,0:l),asgnt
        integer nmeast
        common/mtauvar/gnl,agnl,chinl,achinl,pairsus,
     1       apairsus,asgnt,nmeast
c cdwm.h -- Another silly common. (2dmoo)
        real*8 achicdw,chicdw
        common/cdwm/achicdw,chicdw
c newm.h -- Common. (2dmoo)
        real*8        peep,apeep,phpe,aphpe,phke,aphke,phnn,aphnn
        common/newm/peep,apeep,phpe,aphpe,phke,aphke,phnn,aphnn
c     phonon.h -- variables for the phonon propagator.
      real*8 dnl(0:n/2,0:n/2,0:l),adnl(0:n/2,0:n/2,0:l),asgnph
      common/phonon/dnl,adnl,asgnph
      real*8 axeq
         common/eqpos/axeq

        integer i,j,i2,j2,ti

        nmeas0 = 0
        nmeasp = 0
        nmeast = 0
        call setto0(asgnup,asgndn,asgn,anup,andn)
        call setto0(anud,ake,asaf,ascdw,asferro)
        call setto0(asfer2,asafsq,asafsq2,asgnp,asgnt)
        call setto0(apeep,aphpe,aphke,aphnn,achicdw)
        axeq=0.d0
        asgnph = 0.d0
        do 10 i = 0, n/2
          do 10 j = 0, n/2
            agrfun(i,j) = 0.d0
            aspinxx(i,j) = 0.d0
            acdw(i,j) = 0.d0
            aden(i,j,0) = 0.d0
            aden(i,j,1) = 0.d0
            do 10 ti = 0, l
              agnl(i,j,ti) = 0.d0
              achinl(i,j,ti) = 0.d0
              adnl(i,j,ti) = 0.d0
10      continue
        if(dopair .eq. 1)then
          do 20 i = -1,1
            do 20 j = -1,1
              do 20 i2 = -1,1
                do 20 j2 = -1,1
                  apairmat(i,j,i2,j2) = 0.d0
                  do 20 ti = 0, l
                    apairsus(i,j,i2,j2,ti) = 0.d0
20        continue
        endif

        return
        end
c**************************************************
        subroutine autoset()

c pars.h -- basic Holstein parameters
       use global
        implicit none
c integers.h -- integer parameters for holstein.
        integer warms,sweeps,msr,nwrap
        integer iran
        common/integers/warms,sweeps,msr,nwrap, iran
c
c getgpam.h -- parameters for getg (2dmoo)
c
        integer orthlen,doauto,torth
        real*8 eorth,difflim,errpam
        common/getgpam/eorth,difflim,errpam,orthlen,doauto,
     1                torth
c couple.h -- coupling parameters for holstein.
         real*8 move,t,omega,g,mu,dtau,expmu,gam,lambda,tdtau,dens
         real*8 pe,kee
         common/couple/move,t,omega,g,mu,
     1       dtau,expmu,gam,lambda,tdtau,dens,pe,kee
c vectors.h -- vectors for 2d holstein model.
        real*8 hub(0:volume-1)
        real*8 vup(0:volume-1),vdn(0:volume-1)
        real*8 mphase(0:toff)

        common/vectors/hub,vup,vdn, mphase

        integer start
        real*8 gmatup(0:n*n,0:n*n),diffup,sgnup
        real*8 gmatacc(0:n*n,0:n*n),detup

        write(66,*)'eorth is ',eorth
        start=orthlen
        orthlen = 1
        call ranlat()
        call getgp(vup,0,gmatacc,sgnup,detup)
5       continue
        do 10 orthlen = start, 2, -1
          call getgp(vup,0,gmatup,sgnup,detup)
          call matdif(gmatup,gmatacc,diffup)
          if(doauto .eq. 1)
     1      write(6,*)'diffup for orthlen= ',orthlen,' is ',diffup
          if(doauto .ne. 1 .or. diffup .le. eorth)goto 40
10      continue
        eorth = eorth * 100.d0
        write(66,*)'resetting eorth to ',eorth
        goto 5
40      continue
        write(6,*)'Using orthlen= ',orthlen,' diffup is ',diffup
        write(66,*)'Using orthlen= ',orthlen,' diffup is ',diffup
        return
        end
c
       subroutine setto0(a,b,c,d,e)
       real*8 a,b,c,d,e
       a=0.d0
       b=0.d0
       c=0.d0
       d=0.d0
       e=0.d0
       return
       end
c
       subroutine ftntok(gn,gq,ndim,maxl)
c Fourier transform g(n,l) to get g(q,l).
       implicit none
       real*8 twopi
       parameter(twopi=6.283185307179586)
       integer ndim,maxl
       real*8 gn(0:ndim/2,0:ndim/2,0:maxl)
       real*8 gq(0:ndim/2,0:ndim/2,0:maxl),cfac,cs(0:400)
       integer mx,my,lx,ly,ti,lxp,lyp,first

       if(first .ne. 1273) then
         if(ndim .gt. 20) write(66,*)'Help#1 in ftntok!,',ndim
         first = 1273
         do 100 lx = 0,ndim**2
100        cs(lx) = cos(twopi/ndim*lx)
       endif
       do 10 mx = 0, ndim/2
        do 10 my = 0, ndim/2
         do 15 ti = 0, maxl
15        gq(mx,my,ti) = 0.d0
         do 10 lx = 0, ndim-1
          lxp = min(lx,ndim-lx)
          do 10 ly = 0, ndim-1
           lyp = min(ly,ndim-ly)
           cfac = cs(mx*lx) * cs(my*ly)
           do 10 ti = 0, maxl
            gq(mx,my,ti) = gq(mx,my,ti)+ cfac * gn(lxp,lyp,ti)
10     continue
       return
       end
c
c     Fourier transform f(x) to f(q) with no other garbage.
c
       subroutine ft2ntok(gn,gq,ndim)

       implicit none
       real*8 twopi
       parameter(twopi=6.283185307179586)
       integer ndim
       real*8 gn(0:ndim/2,0:ndim/2)
       real*8 gq(0:ndim/2,0:ndim/2),cfac,cs(0:400)
       integer mx,my,lx,ly,lxp,lyp,first
       
       if(first .ne. 1273) then
         if(ndim .gt. 20) write(66,*)'Help#1 in ftntok!,',ndim
         first = 1273
         do 100 lx = 0,ndim**2
           cs(lx) = cos(twopi/ndim*lx)
100      continue
       endif
       do 10 mx = 0, ndim/2
         do 10 my = 0, ndim/2
           gq(mx,my) = 0.d0
           do 10 lx = 0, ndim-1
             lxp = min(lx,ndim-lx)
             do 10 ly = 0, ndim-1
               lyp = min(ly,ndim-ly)
               cfac = cs(mx*lx) * cs(my*ly)
               gq(mx,my) = gq(mx,my)+ cfac * gn(lxp,lyp)
10     continue
       return
       end
c
       subroutine ftltow(gl,gw,ndim,maxl,dtau,bose,nmax)
c Fourier transform g(n,l) to get g(n,w) (Fermi frequencies).
       implicit none
       real*8 twopi
       parameter(twopi=6.283185307179586)
       integer ndim,maxl,bose,nmax
       real*8 gl(0:ndim/2,0:ndim/2,0:maxl),dtau,omega
       real*8 larray(1000),terpray(1000),tauray(1000)
       real*8 temp,temp2,temp3,rti2
       complex*16 gw(0:ndim/2,0:ndim/2,0:maxl)
       integer mx,my,ti,ti2,intrat
       parameter (intrat=4)
c dtau/intrat is the spacing used for the tau integration.

       do 20 mx = 0, ndim/2
        do 20 my = 0, ndim/2
         do 10 ti = 0, maxl
           larray(ti+1) = gl(mx,my,ti)
10         tauray(ti+1) = ti*dtau
         call spline(tauray,larray,maxl+1,2d30,2d30,terpray)
c Fourier transform over time to get g(n,w).
         do 20 ti = 0, nmax
          if(bose .eq. 0)then
            omega = (ti+0.5d0)*twopi/(maxl)
          else
            omega = ti*twopi/(maxl)
          endif
          call splint(tauray,larray,terpray,maxl+1,0.d0,temp)
          call splint(tauray,larray,terpray,maxl+1,maxl*dtau,temp2)
          rti2 = (intrat*maxl-1.d0)/intrat*dtau
          call splint(tauray,larray,terpray,maxl+1,rti2,temp3)
          gw(mx,my,ti) = temp*dtau / intrat / 3.d0
     1    + temp2*cdexp(dcmplx(0.d0,1.d0)*omega*maxl)*dtau/intrat/3.d0
     2    + temp3*cdexp(dcmplx(0.d0,1.d0)*omega*rti2/dtau)*
     3                                 dtau/intrat*(4.d0/3.d0)
          do 20 ti2 = 1,intrat*maxl-3,2
           call splint(tauray,larray,terpray,maxl+1,
     1                                     ti2*dtau/intrat,temp)
           gw(mx,my,ti) = gw(mx,my,ti) + temp*(4.d0/3.d0) *
     1        cdexp(dcmplx(0.d0,1.d0)*omega*ti2/intrat)*dtau/intrat
           call splint(tauray,larray,terpray,maxl+1,
     1                                (ti2+1)*dtau/intrat,temp)
           gw(mx,my,ti) = gw(mx,my,ti) + temp*(2.d0/3.d0) *
     1       cdexp(dcmplx(0.d0,1.d0)*omega*(ti2+1)/intrat)*dtau/intrat
20      continue
       return
       end
c
       subroutine ftltowp(pairl,pairw,maxl,k,dtau,nmax)
c Fourier transform pair(n,l) to get pair(n,w) (Fermi frequencies).
       use global
       implicit none
       integer maxl,k,j,nmax
       real*8 pairl(nbin,9,0:maxl),dtau,omega
       complex*16 pairw(nbin,9,0:maxl)
       integer ti,ti2,intrat
       real*8 larray(1000),terpray(1000),tauray(1000)
       real*8 temp,temp2,temp3,rti2
       parameter (intrat=4)
       real*8 twopi
        parameter(twopi=6.283185307179586)
c dtau/intrat is the spacing used for the tau integration.


        do 20 j = 1,9
        do 10 ti = 0, maxl
         larray(ti+1) = pairl(k,j,ti)
10       tauray(ti+1) = ti*dtau
        call spline(tauray,larray,maxl+1,2d30,2d30,terpray)
        do 20 ti = 0, nmax
         omega = ti*twopi/(maxl)
         call splint(tauray,larray,terpray,maxl+1,0.d0,temp)
         call splint(tauray,larray,terpray,maxl+1,maxl*dtau,temp2)
         rti2 = (intrat*maxl-1.d0)/intrat*dtau
         call splint(tauray,larray,terpray,maxl+1,rti2,temp3)

          pairw(k,j,ti) = temp*dtau / intrat / 3.d0
     1    + temp2*cdexp(dcmplx(0.d0,1.d0)*omega*maxl)*dtau/intrat/3.d0
     2    + temp3*cdexp(dcmplx(0.d0,1.d0)*omega*rti2/dtau)*
     3                             dtau/intrat*(4.d0/3.d0)


          do 20 ti2 = 1,intrat*maxl-3,2
           call splint(tauray,larray,terpray,maxl+1,
     1                                     ti2*dtau/intrat,temp)

           pairw(k,j,ti) = pairw(k,j,ti) + temp*(4.d0/3.d0) *
     1        cdexp(dcmplx(0.d0,1.d0)*omega*ti2/intrat)*dtau/intrat

           call splint(tauray,larray,terpray,maxl+1,
     1                                (ti2+1)*dtau/intrat,temp)

           pairw(k,j,ti) = pairw(k,j,ti) + temp*(2.d0/3.d0) *
     1       cdexp(dcmplx(0.d0,1.d0)*omega*(ti2+1)/intrat)*dtau/intrat

20      continue
       return
       end
c
       function zbrent(func,x1,x2,tol)
       real*8 func
       parameter (itmax=100,eps=3.e-8)
       a=x1
       b=x2
       fa=func(a)
       fb=func(b)
       if(fb*fa.gt.0.) then
         write(6,*)'root must be bracketed for zbrent.'
         stop
       endif
       fc=fb
       do 10 iter=1,itmax
         if(fb*fc.gt.0.) then
           c=a
           fc=fa
           d=b-a
           e=d
         endif
         if(abs(fc).lt.abs(fb)) then
           a=b
           b=c
           c=a
           fa=fb
           fb=fc
           fc=fa
         endif
         tol1=2.d0*eps*abs(b)+0.5d0*tol
         xm=.5d0*(c-b)
         if(abs(xm).le.tol1 .or. fb.eq.0.)then
           zbrent=b
           return
         endif
         if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
           s=fb/fa
           if(a.eq.c) then
             p=2.d0*xm*s
             q=1.d0-s
           else
             q=fa/fc
             r=fb/fc
             p=s*(2.d0*xm*q*(q-r)-(b-a)*(r-1.d0))
             q=(q-1.d0)*(r-1.d0)*(s-1.d0)
           endif
           if(p.gt.0.) q=-q
           p=abs(p)
           if(2.d0*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
             e=d
             d=p/q
           else
             d=xm
             e=d
           endif
         else
           d=xm
           e=d
         endif
         a=b
         fa=fb
         if(abs(d) .gt. tol1) then
           b=b+d
         else
           b=b+sign(tol1,xm)
         endif
         fb=func(b)
10     continue
       write(6,*)'zbrent exceeding maximum iterations.'
       zbrent=b
       return
       end
c*******************************************
        real*8 function getden(muarg)

c pars.h -- basic Holstein parameters
        use global
        implicit none
c integers.h -- integer parameters for holstein.
        integer warms,sweeps,msr,nwrap
        integer iran
        common/integers/warms,sweeps,msr,nwrap, iran
c
c getgpam.h -- parameters for getg (2dmoo)
c
        integer orthlen,doauto,torth
        real*8 eorth,difflim,errpam
        common/getgpam/eorth,difflim,errpam,orthlen,doauto,
     1                torth
c couple.h -- coupling parameters for holstein.
         real*8 move,t,omega,g,mu,dtau,expmu,gam,lambda,tdtau,dens
         real*8 pe,kee
         common/couple/move,t,omega,g,mu,
     1       dtau,expmu,gam,lambda,tdtau,dens,pe,kee
c kinetic.h -- kinetic energy common block. (2dmoo)
        real*8 ch,soc,ch4
        common/kinetic/ch,soc,ch4
c vectors.h -- vectors for 2d holstein model.
        real*8 hub(0:volume-1)
        real*8 vup(0:volume-1),vdn(0:volume-1)
        real*8 mphase(0:toff)

        common/vectors/hub,vup,vdn, mphase
c m0var.h -- common for measurements. (2dmoo)
        real*8 sgnup,sgndn,sgn,nup,ndn,nud,ke
        real*8 saf,scdw,sferro,sfer2,grfun(0:n/2,0:n/2),safsq,safsq2
        real*8 spinxx(0:n/2,0:n/2),cdw(0:n/2,0:n/2)
        real*8 pairmat(-1:1,-1:1,-1:1,-1:1),errrat,asgnp
        integer nmeas0,redo,noredo,dopair,nmeasp,numpair,dowrite
        real*8 asgnup,asgndn,asgn,anup,andn,anud,ake,asafsq,asafsq2
        real*8 asaf,ascdw,asferro,asfer2,agrfun(0:n/2,0:n/2)
        real*8 aspinxx(0:n/2,0:n/2),acdw(0:n/2,0:n/2)
        real*8 apairmat(-1:1,-1:1,-1:1,-1:1)
        real*8 den(0:n/2,0:n/2,0:1),aden(0:n/2,0:n/2,0:1)
        common/m0var/ sgnup,sgndn,sgn,nup,ndn,nud,ke,safsq,safsq2,
     1       saf,scdw,sferro,sfer2,grfun,
     1       spinxx,cdw,aspinxx,acdw,asafsq,asafsq2,
     2       asgnup,asgndn,asgn,anup,andn,anud,ake,den,aden,
     3       asaf,ascdw,asferro,asfer2,agrfun,errrat,asgnp
     4       ,pairmat,apairmat
        common/m0vari/nmeas0,redo,noredo,dopair,nmeasp,numpair,dowrite

        integer accept,reject,wraps
        integer i
        real*8 gmatup(0:n*n,0:n*n),gmatdn(0:n*n,0:n*n)
        real*8 muarg
        integer denswp
        parameter (denswp=40)

        mu = muarg
        expmu=dexp(dtau*mu)
        call setvup()
        call zeroas()
30      write(6,*)'starting density calc for mu = ',mu
        do 20 i=1,denswp
20        call sweep(gmatup,gmatdn,accept,reject,wraps)
        if(asgn .eq. 0.d0) then
          write(66,*)'asgn = 0, redoing calc.'
          goto 30
        endif
        getden= (anup+andn) / asgn - dens
        write(66,*)'mu,density are ',mu,(anup+andn)/asgn
        write(6,*)'mu,density are ',mu,(anup+andn)/asgn

        return
        end

       subroutine spline(x,y,n,yp1,ypn,y2)
       parameter (nmax=400)
       real*8  x(n),y(n),y2(n),u(nmax)
       real*8 yp1,ypn,sig,p,qn,un
       if (yp1.gt..99e30) then
         y2(1)=0.d0
         u(1)=0.d0
       else
         y2(1)=-0.5d0
         u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
       endif
       do 10 i=2,n-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2.d0
         y2(i)=(sig-1.)/p
         u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     *         /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
10      continue
       if (ypn.gt..99e30) then
         qn=0.d0
         un=0.d0
       else
         qn=0.5d0
         un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
       endif
       y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
       do 20 k=n-1,1,-1
              y2(k)=y2(k)*y2(k+1)+u(k)
20     continue
       return
       end
c
       subroutine splint(xa,ya,y2a,n,x,y)
       real*8 xa(n),ya(n),y2a(n)
       real*8 x,y,h,a,b
       klo=1
       khi=n
1       if (khi-klo.gt.1d0) then
         k=(khi+klo)/2
         if(xa(k).gt.x)then
           khi=k
         else
           klo=k
         endif
       goto 1
       endif
       h=xa(khi)-xa(klo)
       if (h.eq.0.d0) then
         write(6,*) 'bad xa input.'
         stop
       endif
       a=(xa(khi)-x)/h
       b=(x-xa(klo))/h
       y=a*ya(klo)+b*ya(khi)+
     *     ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0
       return
       end
c
c matrix.F -- matrix manipulation routines for 2d holstein model.
c
c*******************************multt(m,lr)********************

c        This subroutine multiplies m by the matrix T on the
c        left or right depending on whether lr = -1 or 1.
c        The factor expmu is also included.
c        This is the real*8 version.

        subroutine multt(mat,lr)

c pars.h -- basic Holstein parameters
        use global
        implicit none
c couple.h -- coupling parameters for holstein.
         real*8 move,t,omega,g,mu,dtau,expmu,gam,lambda,tdtau,dens
         real*8 pe,kee
         common/couple/move,t,omega,g,mu,
     1       dtau,expmu,gam,lambda,tdtau,dens,pe,kee
c kinetic.h -- kinetic energy common block. (2dmoo)
        real*8 ch,soc,ch4
        common/kinetic/ch,soc,ch4

        real*8 rux(0:n*n,0:n*n)
        common/rtmats/rux

        real*8 mat(0:n*n,0:n*n)
        integer i,j,k,lr

        if(lr .eq. -1) then
c Do odd x first
          do 40 i = 1, n*n-3,2
            do 40 j = 0, n*n-1
              rux(i,j) = mat(i,j) + soc*mat(i+1,j)
40            rux(i+1,j) = mat(i+1,j) + soc*mat(i,j)
          do 50 i = n-1, n*n-1,n
            do 50 j = 0, n*n-1

              rux(i,j) = mat(i,j) + soc*mat(i+1-n,j)
50            rux(i+1-n,j) = mat(i+1-n,j) + soc*mat(i,j)

c Do odd y next
          do 60 k = 0, n-1
           do 60 i = n, n*n-3*n,2*n
            do 60 j = 0, n*n-1
              mat(k+i,j) = rux(k+i,j) + soc*rux(k+i+n,j)
60            mat(k+i+n,j) = rux(k+i+n,j) + soc*rux(k+i,j)
          do 70 k = 0, n-1
            do 70 j = 0, n*n-1
              mat(k+n*n-n,j) = rux(k+n*n-n,j) + soc*rux(k,j)
70            mat(k,j) = rux(k,j) + soc*rux(k+n*n-n,j)
c Do even y next
          do 80 k = 0, n-1
           do 80 i = n, n*n-n,2*n
            do 80 j = 0, n*n-1
              rux(k+i,j) = mat(k+i,j) + soc*mat(k+i-n,j)
80            rux(k+i-n,j) = mat(k+i-n,j) + soc*mat(k+i,j)
c Do even x next
          do 90 i = 1, n*n-1,2
            do 90 j = 0, n*n-1
              mat(i,j)=(rux(i,j)+soc*rux(i-1,j)) * expmu*ch4
90            mat(i-1,j)=(rux(i-1,j)+soc*rux(i,j)) * expmu*ch4
        else
c Do even y first, using transposed mat.
          do 180 k = 0, n-1
           do 180 i = n, n*n-n,2*n
            do 180 j = 0, n*n-1
              rux(k+i,j) = mat(j,k+i) + soc*mat(j,k+i-n)
180           rux(k+i-n,j) = mat(j,k+i-n) + soc*mat(j,k+i)
c Do even x next
          do 190 i = 1, n*n-1,2
            do 190 j = 0, n*n-1
              mat(i,j)=(rux(i,j)+soc*rux(i-1,j))*expmu*ch4
190           mat(i-1,j)=(rux(i-1,j)+soc*rux(i,j))*expmu*ch4
c Do odd x next
          do 140 i = 1, n*n-3,2
            do 140 j = 0, n*n-1
              rux(i,j) = mat(i,j) + soc*mat(i+1,j)
140           rux(i+1,j) = mat(i+1,j) + soc*mat(i,j)
          do 150 i = n-1, n*n-1,n
            do 150 j = 0, n*n-1

              rux(i,j) = mat(i,j) + soc*mat(i+1-n,j)
150           rux(i+1-n,j) = mat(i+1-n,j) + soc*mat(i,j)

c Do odd y next
          do 160 k = 0, n-1
           do 160 i = n, n*n-3*n,2*n
            do 160 j = 0, n*n-1
c The following 6 lines includes the transposing of mat(k,j)
              mat(j,k+i) = rux(k+i,j) + soc*rux(k+i+n,j)
160           mat(j,k+i+n) = rux(k+i+n,j) + soc*rux(k+i,j)
          do 170 k = 0, n-1
            do 170 j = 0, n*n-1
              mat(j,k+n*n-n) = rux(k+n*n-n,j) + soc*rux(k,j)
170           mat(j,k) = rux(k,j) + soc*rux(k+n*n-n,j)
        endif

        return
        end
c*******************************multti(m,lr)********************

c        This subroutine multiplies m by the matrix T^-1 on the
c        left or right depending on whether lr = -1 or 1.
c        The factor expmu^-1 is also included.
c        This is the real*8 version.

        subroutine multti(mat,lr)

c pars.h -- basic Holstein parameters
        use global
        implicit none
c couple.h -- coupling parameters for holstein.
         real*8 move,t,omega,g,mu,dtau,expmu,gam,lambda,tdtau,dens
         real*8 pe,kee
         common/couple/move,t,omega,g,mu,
     1       dtau,expmu,gam,lambda,tdtau,dens,pe,kee
c kinetic.h -- kinetic energy common block. (2dmoo)
        real*8 ch,soc,ch4
        common/kinetic/ch,soc,ch4

        real*8 rux(0:n*n,0:n*n)
        common/rtmats/rux

        real*8 mat(0:n*n,0:n*n)
        integer i,j,k,lr

        if(lr .eq. 1) then
c Do odd x first, using transposed mat
          do 40 i = 1, n*n-3,2
            do 40 j = 0, n*n-1
              rux(i,j) = mat(j,i) - soc*mat(j,i+1)
40            rux(i+1,j) = mat(j,i+1) - soc*mat(j,i)
          do 50 i = n-1, n*n-1,n
            do 50 j = 0, n*n-1

              rux(i,j) = mat(j,i) - soc*mat(j,i+1-n)
50            rux(i+1-n,j) = mat(j,i+1-n) - soc*mat(j,i)

c Do odd y next
          do 60 k = 0, n-1
           do 60 i = n, n*n-3*n,2*n
            do 60 j = 0, n*n-1
              mat(k+i,j) = rux(k+i,j) - soc*rux(k+i+n,j)
60            mat(k+i+n,j) = rux(k+i+n,j) - soc*rux(k+i,j)
          do 70 k = 0, n-1
            do 70 j = 0, n*n-1
              mat(k+n*n-n,j) = rux(k+n*n-n,j) - soc*rux(k,j)
70            mat(k,j) = rux(k,j) - soc*rux(k+n*n-n,j)
c Do even y next
          do 80 k = 0, n-1
           do 80 i = n, n*n-n,2*n
            do 80 j = 0, n*n-1
              rux(k+i,j) = mat(k+i,j) - soc*mat(k+i-n,j)
80            rux(k+i-n,j) = mat(k+i-n,j) - soc*mat(k+i,j)
c Do even x next, and transpose back mat.
          do 90 i = 1, n*n-1,2
            do 90 j = 0, n*n-1
              mat(j,i)=(rux(i,j)-soc*rux(i-1,j)) /expmu *ch4
90            mat(j,i-1)=(rux(i-1,j)-soc*rux(i,j)) /expmu *ch4
        else
c for multiplying on left, do even y first.
          do 180 k = 0, n-1
           do 180 i = n, n*n-n,2*n
            do 180 j = 0, n*n-1
              rux(k+i,j) = mat(k+i,j) - soc*mat(k+i-n,j)
180           rux(k+i-n,j) = mat(k+i-n,j) - soc*mat(k+i,j)
c Do even x next
          do 190 i = 1, n*n-1,2
            do 190 j = 0, n*n-1
              mat(i,j)=(rux(i,j)-soc*rux(i-1,j))/expmu* ch4
190           mat(i-1,j)=(rux(i-1,j)-soc*rux(i,j))/expmu* ch4
c Do odd x next
          do 140 i = 1, n*n-3,2
            do 140 j = 0, n*n-1
              rux(i,j) = mat(i,j) - soc*mat(i+1,j)
140           rux(i+1,j) = mat(i+1,j) - soc*mat(i,j)
          do 150 i = n-1, n*n-1,n
            do 150 j = 0, n*n-1

              rux(i,j) = mat(i,j) - soc*mat(i+1-n,j)
150           rux(i+1-n,j) = mat(i+1-n,j) - soc*mat(i,j)

c Do odd y next
          do 160 k = 0, n-1
           do 160 i = n, n*n-3*n,2*n
            do 160 j = 0, n*n-1
              mat(k+i,j) = rux(k+i,j) - soc*rux(k+i+n,j)
160           mat(k+i+n,j) = rux(k+i+n,j) - soc*rux(k+i,j)
          do 170 k = 0, n-1
            do 170 j = 0, n*n-1
              mat(k+n*n-n,j) = rux(k+n*n-n,j) - soc*rux(k,j)
170           mat(k,j) = rux(k,j) - soc*rux(k+n*n-n,j)
        endif

        return
        end
c*******************************multb(m,vvv,ti,lr)********************

c        This subroutine multiplies m by the matrix B(ti) on the
c        left or right depending on whether lr = -1 or 1.
c        This is the real*8 version.

        subroutine multb(mat,vvv,ti,lr)

c pars.h -- basic Holstein parameters
        use global
        implicit none

c impure.h -- common block for impurity.
      real*8 delta0, delta(0:toff-1), expdel(0:toff-1)

      common/impure/ delta0, delta, expdel


        real*8 vvv(0:volume-1)
        real*8 mat(0:n*n,0:n*n)
        integer i,j,ti,lr

        if(lr .eq. -1) then
          call multt(mat,lr)
          do 10 i = 0, n*n-1
            do 10 j = 0, n*n-1
10            mat(i,j) = mat(i,j) * vvv(i+toff*ti)

     @             *expdel(i)

        else
          do 20 j = 0, n*n-1
            do 20 i = 0, n*n-1
20            mat(i,j) = mat(i,j) * vvv(j+toff*ti)

     @             *expdel(j)

          call multt(mat,lr)
        endif

        return
        end
c*******************************multbi(m,vvv,ti,lr)********************

c        This subroutine multiplies m by the matrix B(ti)^-1 on the
c        left or right depending on whether lr = -1 or 1.
c        This is the real*8 version.

        subroutine multbi(mat,vvv,ti,lr)

c pars.h -- basic Holstein parameters
        use global
        implicit none

c impure.h -- common block for impurity.
      real*8 delta0, delta(0:toff-1), expdel(0:toff-1)

      common/impure/ delta0, delta, expdel


        real*8 vvv(0:volume-1)
        real*8 mat(0:n*n,0:n*n)
        integer i,j,ti,lr

        if(lr .eq. 1) then
          call multti(mat,lr)
          do 10 j = 0, n*n-1
            do 10 i = 0, n*n-1
10            mat(i,j) = mat(i,j) / vvv(j+toff*ti)

     @             /expdel(j)

        else
          do 20 i = 0, n*n-1
            do 20 j = 0, n*n-1
20            mat(i,j) = mat(i,j) / vvv(i+toff*ti)

     @             /expdel(i)

          call multti(mat,lr)
        endif

        return
        end
c       ******************unit(mat)********************
c        This subroutine puts the unit matrix in aux, which is real.

        subroutine unit(mat)

c pars.h -- basic Holstein parameters
        use global
        implicit none
        real*8 mat(0:n*n,0:n*n)
        integer i,j
        do 10 i = 0, n*n-1
          do 20 j = 0, n*n-1
20          mat(i,j) = 0.d0
10        mat(i,i) = 1.d0
        return
        end
c       ******************zeromat(mat)********************
c        This subroutine zeroes mat.

        subroutine zeromat(mat)

c pars.h -- basic Holstein parameters
        use global
        implicit none

        real*8 mat(0:n*n,0:n*n)
        integer i,j
        do 20 i = 0, n*n-1
          do 20 j = 0, n*n-1
20          mat(i,j) = 0.d0
        return
        end
c **************************TRANSP(MAT,AUX)**********************
c        This subroutine multiplies mat by the even-x matrix

        subroutine transp(mat,aux)

c pars.h -- basic Holstein parameters
        use global
        implicit none

        integer i,j
        real*8 mat(0:n*n,0:n*n),aux(0:n*n,0:n*n)

        do 20 i = 0, n*n-1, 4
          do 20 j = 0, n*n-1
            aux(i,j)   = mat(j,i)
            aux(i+1,j) = mat(j,i+1)
            aux(i+2,j) = mat(j,i+2)
            aux(i+3,j) = mat(j,i+3)
20      continue
        return
        end


c**********************************
       subroutine matmult(mat1,mat2,mat3)

c pars.h -- basic Holstein parameters
       use global
        implicit none

       integer i,j,k
       real*8 mat1(0:n*n,0:n*n),mat2(0:n*n,0:n*n)
       real*8 mat3(0:n*n,0:n*n)
c
       do 10 i = 0, n*n-1
         do 10 j = 0, n*n-1
10         mat3(i,j) = 0.d0
       do 20 j = 0, n*n-1
         do 20 k = 0, n*n-1
           do 20 i = 0, n*n-1
20           mat3(i,j) = mat3(i,j) + mat1(i,k)*mat2(k,j)
       return
       end
c**********************************
       subroutine matdif(mat1,mat2,diff)

c pars.h -- basic Holstein parameters
       use global
        implicit none

       integer i,j
       real*8 mat1(0:n*n,0:n*n),mat2(0:n*n,0:n*n),diff
c
       diff = 0.d0
       do 10 i = 0, n*n-1
         do 10 j = 0, n*n-1
10         diff = diff + (mat1(i,j) - mat2(i,j))**2
       diff = dsqrt(diff) / (n*n)
       return
       end
c********************************
       subroutine matcop(mat1,mat2)

c pars.h -- basic Holstein parameters
       use global
        implicit none

       integer i,j
       real*8 mat1(0:n*n,0:n*n),mat2(0:n*n,0:n*n)
c
       do 10 i = 0, n*n-1
         do 10 j = 0, n*n-1
10         mat2(i,j) = mat1(i,j)
       return
       end
c**************************************
       subroutine donorm(vec,vnorm)

c pars.h -- basic Holstein parameters
       use global
        implicit none
    
       integer i
       real*8 vec(0:n*n),vnorm,temp
c
       temp = 0.d0
       do 10 i = 0, n*n-1
10       temp = temp + vec(i)**2
       vnorm = dsqrt(temp)
       if(vnorm .ne. 0.d0) then
         temp = 1.d0 / vnorm
         do 20 i = 0, n*n-1
20         vec(i) = vec(i) * temp
       endif
       return
       end

c
c meas.F -- Measurement routines for 2d holstein model.
c
c*******************************meas0(lots of variables)********************

c        This subroutine does the measurements.

        subroutine meas0(gmatup,gmatdn,ti)

c pars.h -- basic Holstein parameters
        use global
        implicit none
    
c integers.h -- integer parameters for holstein.
        integer warms,sweeps,msr,nwrap
        integer iran
        common/integers/warms,sweeps,msr,nwrap, iran
c couple.h -- coupling parameters for holstein.
         real*8 move,t,omega,g,mu,dtau,expmu,gam,lambda,tdtau,dens
         real*8 pe,kee
         common/couple/move,t,omega,g,mu,
     1       dtau,expmu,gam,lambda,tdtau,dens,pe,kee
c ivectors.h -- index vectors for indirect addressing and phase vector
        integer xplus(0:toff-1),xminus(0:toff-1)
        integer yplus(0:toff-1),yminus(0:toff-1)

        common/ivectors/xplus,xminus,yplus,yminus
c vectors.h -- vectors for 2d holstein model.
        real*8 hub(0:volume-1)
        real*8 vup(0:volume-1),vdn(0:volume-1)
        real*8 mphase(0:toff)

        common/vectors/hub,vup,vdn, mphase
c m0var.h -- common for measurements. (2dmoo)
        real*8 sgnup,sgndn,sgn,nup,ndn,nud,ke
        real*8 saf,scdw,sferro,sfer2,grfun(0:n/2,0:n/2),safsq,safsq2
        real*8 spinxx(0:n/2,0:n/2),cdw(0:n/2,0:n/2)
        real*8 pairmat(-1:1,-1:1,-1:1,-1:1),errrat,asgnp
        integer nmeas0,redo,noredo,dopair,nmeasp,numpair,dowrite
        real*8 asgnup,asgndn,asgn,anup,andn,anud,ake,asafsq,asafsq2
        real*8 asaf,ascdw,asferro,asfer2,agrfun(0:n/2,0:n/2)
        real*8 aspinxx(0:n/2,0:n/2),acdw(0:n/2,0:n/2)
        real*8 apairmat(-1:1,-1:1,-1:1,-1:1)
        real*8 den(0:n/2,0:n/2,0:1),aden(0:n/2,0:n/2,0:1)
        common/m0var/ sgnup,sgndn,sgn,nup,ndn,nud,ke,safsq,safsq2,
     1       saf,scdw,sferro,sfer2,grfun,
     1       spinxx,cdw,aspinxx,acdw,asafsq,asafsq2,
     2       asgnup,asgndn,asgn,anup,andn,anud,ake,den,aden,
     3       asaf,ascdw,asferro,asfer2,agrfun,errrat,asgnp
     4       ,pairmat,apairmat
        common/m0vari/nmeas0,redo,noredo,dopair,nmeasp,numpair,dowrite

c newm.h -- Common. (2dmoo)
        real*8        peep,apeep,phpe,aphpe,phke,aphke,phnn,aphnn
        common/newm/peep,apeep,phpe,aphpe,phke,aphke,phnn,aphnn

        real*8 ctemp
        real*8 gmatup(0:n*n,0:n*n),gmatdn(0:n*n,0:n*n)
        real*8 safch,gtemp,xxtemp,zztemp,xfac,yfac

        integer ti,i,j,ix,jx,iy,jy,kx,ky

        saf = 0.d0
        safsq = 0.d0
        sferro = 0.d0
        sfer2 = 0.d0
        scdw = 0.d0
        safsq2 = 0.d0
        safch = 0.d0
        nup = 0.d0
        ndn = 0.d0
        ke = 0.d0
        nud = 0.d0
        scdw=0.d0
        phke=0.d0
        phpe=0.d0
        phnn=0.d0
        peep=0.d0

          do 20 i = 0,n*n-1
           do 20 j = 0,n*n-1
20           gmatdn(i,j) = gmatup(i,j)

        do 30 kx = 0, n/2
          do 30 ky = 0, n/2
            grfun(kx,ky) = 0.d0
            spinxx(kx,ky) = 0.d0
            cdw(kx,ky) = 0.d0
            den(kx,ky,0) = 0.d0
            den(kx,ky,1) = 0.d0
30      continue

        do 200 i = 0,n*n-1
             peep=peep+
     1     hub(i+ti*toff)*(2.d0-gmatup(i,i)-gmatdn(i,i))
          nup = nup + 1.d0 - gmatup(i,i)
          ndn = ndn + 1.d0 - gmatdn(i,i)
          ke=ke+gmatup(i,xplus(i))+gmatup(i,xminus(i))
     1              +gmatup(i,yplus(i))+gmatup(i,yminus(i))
     2              +gmatdn(i,xplus(i))+gmatdn(i,xminus(i))
     3              +gmatdn(i,yplus(i))+gmatdn(i,yminus(i))
          nud=nud+(1.d0-gmatup(i,i))*(1.d0-gmatdn(i,i))
          sferro = sferro + gmatup(i,i)+gmatdn(i,i)
          sfer2 = sfer2 + gmatup(i,i)+gmatdn(i,i)
          saf = saf + gmatup(i,i)+gmatdn(i,i)
200       scdw = scdw + gmatup(i,i)+gmatdn(i,i)
          peep=peep*g
          do 203 i=0,volume-1

              phpe=phpe+pe*hub(i)*hub(i)

203       continue
          phpe=phpe/dtau
          do 204 i=toff,volume-toff-1
              phke=phke+(hub(i)-hub(i-toff))*(hub(i)-hub(i-toff))
              phke=phke+(hub(i)-hub(i+toff))*(hub(i)-hub(i+toff))
204       continue
          do 205 i=0,toff-1
              phke=phke+(hub(i)-hub(i+(L-1)*toff))
     1                    *(hub(i)-hub(i+(L-1)*toff))
              phke=phke+(hub(i)-hub(i+toff))*(hub(i)-hub(i+toff))
205       continue
          do 206 i=volume-toff,volume-1
              phke=phke+(hub(i)-hub(i-toff))*(hub(i)-hub(i-toff))
              phke=phke+(hub(i)-hub(i-(L-1)*toff))
     1                 *(hub(i)-hub(i-(L-1)*toff))
206       continue
          phke=kee*kee*phke
          do 208 ix=0,L-1
          do 207 i=0,toff-1
              iy=ix*toff
              phnn=phnn+hub(i+iy)*(hub(xplus(i)+iy)+hub(xminus(i)+iy)
     1                            +hub(yplus(i)+iy)+hub(yminus(i)+iy))
207       continue
208       continue

        spinxx(0,0) = saf
        cdw(0,0) = saf
        do 210 ix = 0,n-1
         do 210 iy = 0,n-1
          i = ix+n*iy
          do 210 jx = 0,n-1
           do 210 jy = 0,n-1
            j = jx+n*jy
            kx = abs(ix-jx)
            kx = min(kx,n-kx)
            ky = abs(iy-jy)
            ky = min(ky,n-ky)
            gtemp = gmatup(i,j)+ gmatdn(i,j)
            grfun(kx,ky) = grfun(kx,ky) + gtemp
c The following line is not correct for kx=ky=0, but it will be fixed later.
            den(kx,ky,0) = den(kx,ky,0)+
     1             (1.d0-gmatup(i,i))*(1.d0-gmatup(j,j)) +
     2             (1.d0-gmatdn(i,i))*(1.d0-gmatdn(j,j)) -
     3           gmatup(i,j)*gmatup(j,i)-gmatdn(i,j)*gmatdn(j,i)
            den(kx,ky,1) = den(kx,ky,1)+
     1             (1.d0-gmatup(i,i))*(1.d0-gmatdn(j,j))
            xxtemp = -2.d0*gmatup(i,j)*gmatdn(j,i)
            zztemp = (gmatup(i,i)*gmatup(j,j)
     1        +gmatdn(i,i)*gmatdn(j,j)-2.d0*gmatup(i,i)*gmatdn(j,j)
     2        -gmatup(j,i)*gmatup(i,j) - gmatdn(j,i)*gmatdn(i,j))
            spinxx(kx,ky) = spinxx(kx,ky) + xxtemp
            ctemp = (gmatup(i,i)*gmatup(j,j)
     1        +gmatdn(i,i)*gmatdn(j,j)+2.d0*gmatup(i,i)*gmatdn(j,j)
     2        -gmatup(j,i)*gmatup(i,j) - gmatdn(j,i)*gmatdn(i,j))
            cdw(kx,ky) = cdw(kx,ky) + ctemp
            cdw(kx,ky) = cdw(kx,ky) + 4.d0 - 2.d0*( gmatup(i,i)
     1        +gmatup(j,j)+gmatdn(i,i)+gmatdn(j,j)  )
            sferro=sferro + xxtemp
            sfer2=sfer2 + zztemp
            saf=saf + mphase(i)*mphase(j)*xxtemp
210         scdw=scdw + mphase(i)*mphase(j)*ctemp

        saf = saf / toff
        sferro = sferro / toff
        sfer2 = sfer2 / toff
        scdw = scdw / toff
        ke = ke*t / toff
        nud = nud / toff
        nup = nup / toff
        ndn = ndn / toff
        phpe=phpe/volume
        peep=peep/volume
        phke=phke/volume
        phnn=phnn/volume/4.
        safsq = saf**2

        do 41 kx = 0, n/2
          if(kx .eq. 0 .or. kx .eq. n/2)then
            xfac = 1.d0
          else
            xfac = 2.d0
          endif
          do 41 ky = 0, n/2
            if(ky .eq. 0 .or. ky .eq. n/2)then
              yfac = 1.d0
            else
              yfac = 2.d0
            endif
            spinxx(kx,ky) = spinxx(kx,ky) / (n*n*xfac*yfac)
            cdw(kx,ky) = cdw(kx,ky) / (n*n*xfac*yfac)
            grfun(kx,ky) = grfun(kx,ky) / (2*n*n*xfac*yfac)
            den(kx,ky,0) = den(kx,ky,0) / (2*n*n*xfac*yfac)
            den(kx,ky,1) = den(kx,ky,1) / (n*n*xfac*yfac)
41      continue

c See note above:
        den(0,0,0) = (nup + ndn) * 0.5d0

        do 31 kx = 0, n/2
          do 31 ky = 0, n/2
            grfun(kx,ky) = (grfun(kx,ky)+grfun(ky,kx))*0.5d0
            grfun(ky,kx) = grfun(kx,ky)
            den(kx,ky,0) = (den(kx,ky,0)+den(ky,kx,0))*0.5d0
            den(ky,kx,0) = den(kx,ky,0)
            den(kx,ky,1) = (den(kx,ky,1)+den(ky,kx,1))*0.5d0
            den(ky,kx,1) = den(kx,ky,1)
            spinxx(kx,ky) = (spinxx(kx,ky) + spinxx(ky,kx))*0.5d0
            spinxx(ky,kx) = spinxx(kx,ky)
            cdw(kx,ky) = (cdw(kx,ky) + cdw(ky,kx))*0.5d0
            cdw(ky,kx) = cdw(kx,ky)
31      continue

        nmeas0 = nmeas0 + 1
        sgn = sgnup*sgndn
        anup = anup + nup*sgn
        andn = andn + ndn*sgn
        asaf = asaf + saf*sgn
        asafsq = asafsq + safsq*sgn
        asferro = asferro + sferro*sgn
        asfer2 = asfer2 + sfer2*sgn
        ascdw = ascdw + scdw*sgn
        asafsq2 = asafsq2 + safsq2*sgn
        anud = anud + nud*sgn
        ake = ake + ke*sgn
        asgn = asgn + sgn
        asgnup = asgnup + sgnup
        asgndn = asgndn + sgndn
        do 51 kx = 0, n/2
          do 51 ky = 0, n/2
            agrfun(kx,ky) = agrfun(kx,ky) + grfun(kx,ky)*sgn
            aden(kx,ky,0) = aden(kx,ky,0) + den(kx,ky,0)*sgn
            aden(kx,ky,1) = aden(kx,ky,1) + den(kx,ky,1)*sgn
            aspinxx(kx,ky) = aspinxx(kx,ky) + spinxx(kx,ky)*sgn
51          acdw(kx,ky) = acdw(kx,ky) + cdw(kx,ky)*sgn
        aphpe=aphpe+phpe*sgn
        apeep=apeep+peep*sgn
        aphke=aphke+phke*sgn
        aphnn=aphnn+phnn*sgn

        return
        end
c*******************************measpair(gmatup,gmatdn)********************

c        This subroutine does the measurements.

        subroutine measpair(gmatup,gmatdn)

c pars.h -- basic Holstein parameters
        use global
        implicit none
c integers.h -- integer parameters for holstein.
        integer warms,sweeps,msr,nwrap
        integer iran
        common/integers/warms,sweeps,msr,nwrap, iran
c couple.h -- coupling parameters for holstein.
         real*8 move,t,omega,g,mu,dtau,expmu,gam,lambda,tdtau,dens
         real*8 pe,kee
         common/couple/move,t,omega,g,mu,
     1       dtau,expmu,gam,lambda,tdtau,dens,pe,kee
c ivectors.h -- index vectors for indirect addressing and phase vector
        integer xplus(0:toff-1),xminus(0:toff-1)
        integer yplus(0:toff-1),yminus(0:toff-1)

        common/ivectors/xplus,xminus,yplus,yminus
c vectors.h -- vectors for 2d holstein model.
        real*8 hub(0:volume-1)
        real*8 vup(0:volume-1),vdn(0:volume-1)
        real*8 mphase(0:toff)

        common/vectors/hub,vup,vdn, mphase
c m0var.h -- common for measurements. (2dmoo)
        real*8 sgnup,sgndn,sgn,nup,ndn,nud,ke
        real*8 saf,scdw,sferro,sfer2,grfun(0:n/2,0:n/2),safsq,safsq2
        real*8 spinxx(0:n/2,0:n/2),cdw(0:n/2,0:n/2)
        real*8 pairmat(-1:1,-1:1,-1:1,-1:1),errrat,asgnp
        integer nmeas0,redo,noredo,dopair,nmeasp,numpair,dowrite
        real*8 asgnup,asgndn,asgn,anup,andn,anud,ake,asafsq,asafsq2
        real*8 asaf,ascdw,asferro,asfer2,agrfun(0:n/2,0:n/2)
        real*8 aspinxx(0:n/2,0:n/2),acdw(0:n/2,0:n/2)
        real*8 apairmat(-1:1,-1:1,-1:1,-1:1)
        real*8 den(0:n/2,0:n/2,0:1),aden(0:n/2,0:n/2,0:1)
        common/m0var/ sgnup,sgndn,sgn,nup,ndn,nud,ke,safsq,safsq2,
     1       saf,scdw,sferro,sfer2,grfun,
     1       spinxx,cdw,aspinxx,acdw,asafsq,asafsq2,
     2       asgnup,asgndn,asgn,anup,andn,anud,ake,den,aden,
     3       asaf,ascdw,asferro,asfer2,agrfun,errrat,asgnp
     4       ,pairmat,apairmat
        common/m0vari/nmeas0,redo,noredo,dopair,nmeasp,numpair,dowrite

        real*8 gmatup(0:n*n,0:n*n),gmatdn(0:n*n,0:n*n)

        integer i,j,mx,my,mpx,mpy
        integer inn,in0,in1,i0n,i01,i1n,i10,i11,jp(-1:1),jpp

        do 10 i = 0,n*n-1
          do 20 j = 0,n*n-1
            gmatdn(i,j) = gmatup(i,j)
20        continue
10      continue

        do 294 mx = -1,1
         do 294 my = -1,1
          do 294 mpx = -1,1
           do 294 mpy = -1,1
294          pairmat(mx,my,mpx,mpy) = 0.d0
        do 522 i = 0,n*n-1
          inn = xminus(yminus(i))
          in0 = xminus(i)
          in1 = xminus(yplus(i))
          i0n = yminus(i)
          i01 = yplus(i)
          i1n = xplus(yminus(i))
          i10 = xplus(i)
          i11 = xplus(yplus(i))
          do 522 j = 0,n*n-1

            jpp = xplus(j)
            mx = 0

             jpp = xminus(jpp)
             jp(-1) = yplus(jpp)
             jp(0) = jpp
             jp(1) = yminus(jpp)

             my = 0


               pairmat(mx,my,0,0) = pairmat(mx,my,0,0) +
     1             gmatup(jp(my),i)*gmatdn(j,i)

500         continue
522     continue

        nmeasp = nmeasp + 1
        do 494 mx = -1,1
         do 494 my = -1,1
          do 494 mpx = -1,1
           do 494 mpy = -1,1
            pairmat(mx,my,mpx,mpy)=pairmat(mx,my,mpx,mpy)/toff
494         apairmat(mx,my,mpx,mpy)=apairmat(mx,my,mpx,mpy) +
     1                  sgnup*sgndn*pairmat(mx,my,mpx,mpy)
        asgnp = asgnp + sgnup*sgndn

        return
        end
c*******************************meastau********************

c        This subroutine does the measurements.

        subroutine meastau

c pars.h -- basic Holstein parameters
        use global
        implicit none
c integers.h -- integer parameters for holstein.
        integer warms,sweeps,msr,nwrap
        integer iran
        common/integers/warms,sweeps,msr,nwrap, iran
c
c getgpam.h -- parameters for getg (2dmoo)
c
        integer orthlen,doauto,torth
        real*8 eorth,difflim,errpam
        common/getgpam/eorth,difflim,errpam,orthlen,doauto,
     1                torth
c couple.h -- coupling parameters for holstein.
         real*8 move,t,omega,g,mu,dtau,expmu,gam,lambda,tdtau,dens
         real*8 pe,kee
         common/couple/move,t,omega,g,mu,
     1       dtau,expmu,gam,lambda,tdtau,dens,pe,kee
c ivectors.h -- index vectors for indirect addressing and phase vector
        integer xplus(0:toff-1),xminus(0:toff-1)
        integer yplus(0:toff-1),yminus(0:toff-1)

        common/ivectors/xplus,xminus,yplus,yminus
c vectors.h -- vectors for 2d holstein model.
        real*8 hub(0:volume-1)
        real*8 vup(0:volume-1),vdn(0:volume-1)
        real*8 mphase(0:toff)

        common/vectors/hub,vup,vdn, mphase
c m0var.h -- common for measurements. (2dmoo)
        real*8 sgnup,sgndn,sgn,nup,ndn,nud,ke
        real*8 saf,scdw,sferro,sfer2,grfun(0:n/2,0:n/2),safsq,safsq2
        real*8 spinxx(0:n/2,0:n/2),cdw(0:n/2,0:n/2)
        real*8 pairmat(-1:1,-1:1,-1:1,-1:1),errrat,asgnp
        integer nmeas0,redo,noredo,dopair,nmeasp,numpair,dowrite
        real*8 asgnup,asgndn,asgn,anup,andn,anud,ake,asafsq,asafsq2
        real*8 asaf,ascdw,asferro,asfer2,agrfun(0:n/2,0:n/2)
        real*8 aspinxx(0:n/2,0:n/2),acdw(0:n/2,0:n/2)
        real*8 apairmat(-1:1,-1:1,-1:1,-1:1)
        real*8 den(0:n/2,0:n/2,0:1),aden(0:n/2,0:n/2,0:1)
        common/m0var/ sgnup,sgndn,sgn,nup,ndn,nud,ke,safsq,safsq2,
     1       saf,scdw,sferro,sfer2,grfun,
     1       spinxx,cdw,aspinxx,acdw,asafsq,asafsq2,
     2       asgnup,asgndn,asgn,anup,andn,anud,ake,den,aden,
     3       asaf,ascdw,asferro,asfer2,agrfun,errrat,asgnp
     4       ,pairmat,apairmat
        common/m0vari/nmeas0,redo,noredo,dopair,nmeasp,numpair,dowrite

c cdwm.h -- Another silly common. (2dmoo)
        real*8 achicdw,chicdw
        common/cdwm/achicdw,chicdw
c mtauvar.h -- common for tau-dependent measurements. (2dmoo)
        real*8 gnl(0:n/2,0:n/2,0:l),agnl(0:n/2,0:n/2,0:l)
        real*8 chinl(0:n/2,0:n/2,0:l),achinl(0:n/2,0:n/2,0:l)
        real*8 pairsus(-1:1,-1:1,-1:1,-1:1,0:l)
        real*8 apairsus(-1:1,-1:1,-1:1,-1:1,0:l),asgnt
        integer nmeast
        common/mtauvar/gnl,agnl,chinl,achinl,pairsus,
     1       apairsus,asgnt,nmeast

        integer ti,i,j,ix,jx,iy,jy,kx,ky,sx,sy,mx,my,mpx,mpy
        integer inn,in0,in1,i0n,i01,i1n,i10,i11,jp(-1:1),jpp
        real*8 gupt0(0:n*n,0:n*n),gup0t(0:n*n,0:n*n),chitemp
        real*8 gdnt0(0:n*n,0:n*n),gdn0t(0:n*n,0:n*n),gtemp
        real*8 gmatupc(0:n*n,0:n*n),g00(0:n*n,0:n*n)
        real*8 sgnupc,detup
        integer inacc,wraps,ti1

100     do 15 ti = 0, l
         do 15 mx = -1,1
          do 15 my = -1,1
           do 15 mpx = -1,1
            do 15 mpy = -1,1
15           pairsus(mx,my,mpx,mpy,ti) = 0.d0
        chicdw=0.d0
          
      
        do 10 ti = 0, l-1
         call makegt(0,ti,gupt0,gup0t,inacc)
         call matcop(gupt0,gdnt0)
         call matcop(gup0t,gdn0t)
         if(inacc .eq. 1) then
           torth = max(1,torth-1)
           write(6,*)'reducing torth to ',torth
           write(66,*)'reducing torth to ',torth
           goto 100
         endif
         do 40 kx = 0,n/2
          do 40 ky = 0, n/2
           gtemp = 0.d0
           chitemp = 0.d0
           do 30 sx = -1, 1, 2
            do 30 sy = -1, 1, 2
             do 30 ix = 0,n-1
              do 30 iy = 0,n-1
               jx = mod(ix+sx*kx+n,n)
               jy = mod(iy+sy*ky+n,n)
               i = ix+n*iy
               j = jx+n*jy
               gtemp = gtemp + gupt0(i,j)+ gdnt0(i,j)
               chitemp = chitemp - gup0t(j,i)*gdnt0(i,j)
     1                           - gdn0t(j,i)*gupt0(i,j)
30         continue
           gnl(kx,ky,ti) = gtemp / (8*n*n)
           chinl(kx,ky,ti) = chitemp / (4*n*n)
           if(ti .eq. 0) then
             gnl(kx,ky,l) = -gnl(kx,ky,0)
             chinl(kx,ky,l) = chinl(kx,ky,0)
           endif
40       continue
         gnl(0,0,l) = 1.d0 - gnl(0,0,0)

c       BEGIN CDW CALCULATION:

        if (ti.eq.0) then
        do 2020 i=0,n*n-1
             chicdw = chicdw + 2.d0*gupt0(i,i)
2020    continue
        do 2031 i=0,n*n-1
        do 2030 j=0,n*n-1
             chicdw = chicdw + 2.d0*gupt0(j,i)*gup0t(i,j)
     1              *mphase(i)*mphase(j)
2030    continue
2031    continue
        else
          do 2034 i=0,n*n
          do 2033 j=0,n*n
               chicdw=chicdw-2.d0*gupt0(j,i)*gup0t(i,j)
     1            *mphase(i)*mphase(j)
2033      continue
2034      continue
         endif

        if(dopair .eq. 1) then
        do 522 i = 0,n*n-1
          inn = xminus(yminus(i))
          in0 = xminus(i)
          in1 = xminus(yplus(i))
          i0n = yminus(i)
          i01 = yplus(i)
          i1n = xplus(yminus(i))
          i10 = xplus(i)
          i11 = xplus(yplus(i))
          do 522 j = 0,n*n-1

            mx = 0
            jpp = j

             jp(-1) = yplus(jpp)
             jp(0) = jpp
             jp(1) = yminus(jpp)

             my = 0


               pairsus(mx,my,0,0,ti) = pairsus(mx,my,0,0,ti) +
     1             gupt0(jp(my),i)*gdnt0(j,i)

510         continue
            if(ti .eq. 0)then
c g(beta) = 1-g(0); but for ti=0, g0t(0) = g(0)-1, minus signs cancel.

             my =0


               pairsus(mx,my,0,0,l) = pairsus(mx,my,0,0,l) +
     1             gup0t(jp(my),i)*gdn0t(j,i)

515          continue
            endif
500         continue
522     continue
        endif

10     continue

c
        nmeast = nmeast + 1
        asgnt = asgnt + sgnup*sgndn
        do 594 mx = 0, n/2
         do 594 my = 0, n/2
          do 594 ti = 0,l
           agnl(mx,my,ti)=agnl(mx,my,ti)+gnl(mx,my,ti)*sgnup*sgndn
           achinl(mx,my,ti)=achinl(mx,my,ti)+
     1                          chinl(mx,my,ti)*sgnup*sgndn
594     continue
        do 494 mx = -1,1
        do 494 my = -1,1
         do 494 mpx = -1,1
         do 494 mpy = -1,1
          do 494 ti = 0, l
          pairsus(mx,my,mpx,mpy,ti)=pairsus(mx,my,mpx,mpy,ti)/toff
494       apairsus(mx,my,mpx,mpy,ti)=apairsus(mx,my,mpx,mpy,ti) +
     1                  sgnup*sgndn*pairsus(mx,my,mpx,mpy,ti)

c       COMPLETE THE CDW MEASUREMENTS:
c       FIRST GET G(0,0)

        call getgp(vup,0,gmatupc,sgnupc,detup)
        call matcop(gmatupc,g00)

c       THEN GET G(tau,tau) BY WRAPPING (OR RECALCULATING)
c       DO THE MEASUREMENTS

          do 1010 ti = 0, l-1
          call multb(gmatupc,vup,ti,-1)
          call multbi(gmatupc,vup,ti,1)
          do 1005 i=0,toff-1
          do 1005 j=0,toff-1
               chicdw=chicdw+4.d0*(1.d0-gmatupc(i,i))*(1.d0-g00(j,j))
     1               *mphase(i)*mphase(j)
1005      continue
          if(wraps .gt. 0) then
            wraps = wraps - 1
          else
            wraps = nwrap
            ti1 = mod(ti+1,l)
            call getgp(vup,ti1,gmatupc,sgnupc,detup)
          endif
1010        continue

          chicdw=chicdw*dtau/toff
          achicdw=achicdw+chicdw*sgnup*sgndn

        return
        end
c
c*****Measure the phonon green's function. A seperate routine so*******
c     we can measure it more often.
c  
      subroutine measphon

c pars.h -- basic Holstein parameters
      USE GLOBAL
        implicit none
c couple.h -- coupling parameters for holstein.
         real*8 move,t,omega,g,mu,dtau,expmu,gam,lambda,tdtau,dens
         real*8 pe,kee
         common/couple/move,t,omega,g,mu,
     1       dtau,expmu,gam,lambda,tdtau,dens,pe,kee
c vectors.h -- vectors for 2d holstein model.
        real*8 hub(0:volume-1)
        real*8 vup(0:volume-1),vdn(0:volume-1)
        real*8 mphase(0:toff)

        common/vectors/hub,vup,vdn, mphase
c m0var.h -- common for measurements. (2dmoo)
        real*8 sgnup,sgndn,sgn,nup,ndn,nud,ke
        real*8 saf,scdw,sferro,sfer2,grfun(0:n/2,0:n/2),safsq,safsq2
        real*8 spinxx(0:n/2,0:n/2),cdw(0:n/2,0:n/2)
        real*8 pairmat(-1:1,-1:1,-1:1,-1:1),errrat,asgnp
        integer nmeas0,redo,noredo,dopair,nmeasp,numpair,dowrite
        real*8 asgnup,asgndn,asgn,anup,andn,anud,ake,asafsq,asafsq2
        real*8 asaf,ascdw,asferro,asfer2,agrfun(0:n/2,0:n/2)
        real*8 aspinxx(0:n/2,0:n/2),acdw(0:n/2,0:n/2)
        real*8 apairmat(-1:1,-1:1,-1:1,-1:1)
        real*8 den(0:n/2,0:n/2,0:1),aden(0:n/2,0:n/2,0:1)
        common/m0var/ sgnup,sgndn,sgn,nup,ndn,nud,ke,safsq,safsq2,
     1       saf,scdw,sferro,sfer2,grfun,
     1       spinxx,cdw,aspinxx,acdw,asafsq,asafsq2,
     2       asgnup,asgndn,asgn,anup,andn,anud,ake,den,aden,
     3       asaf,ascdw,asferro,asfer2,agrfun,errrat,asgnp
     4       ,pairmat,apairmat
        common/m0vari/nmeas0,redo,noredo,dopair,nmeasp,numpair,dowrite
      real*8 dnl(0:n/2,0:n/2,0:l),adnl(0:n/2,0:n/2,0:l),asgnph
      common/phonon/dnl,adnl,asgnph

c     phonon.h -- variables for the phonon propagator.

      integer ti, kx, ky, sx, sy, ix, iy, jx, jy, i, j, ti1, tii
      real*8 gtemp

      integer count
c        MAKE SURE YOU CHANGE THIS EVERYWHERE IN THE CODE!
c        MAKE SURE YOU CHANGE THIS EVERYWHERE IN THE CODE!
         common/acphononi/count
         real*8 adnl2(0:n/2,0:n/2,0:l),asgnph2
         common/acphononr/adnl2,asgnph2
         real*8 adql2(0:n/2,0:n/2,0:l)
         real*8 atemp


         real*8 xeq,axeq
         common/eqpos/axeq
c
c     Main Loop

      
     
      if (mod(count,nbin).eq.0) then
c         write (88,*) '   '
c         write (88,*) 'setting adnl2;aql2=0 at count=  ',count
c         write (88,*) '   '
         do 2100 ti = 0, l-1
         do 2200 kx = 0,n/2
         do 2200 ky = 0, n/2
              adnl2(kx,ky,ti)=0.d0
              adql2(kx,ky,ti)=0.d0
2200     continue
2100     continue
         asgnph2=0.d0
      endif

      do 100 ti = 0, l-1
        do 200 kx = 0,n/2
          do 200 ky = 0, n/2
            gtemp = 0.d0
            do 300 sx = -1, 1, 2
              do 300 sy = -1, 1, 2
                do 300 ix = 0,n-1
                  do 300 iy = 0,n-1
                    jx = mod(ix+sx*kx+n,n)
                    jy = mod(iy+sy*ky+n,n)
                    i = ix+n*iy
                    j = jx+n*jy
                    do 300 ti1 = 0, l-1
                      tii = mod(ti+ti1,l)
                      gtemp = gtemp + hub(i+tii*toff)*hub(j+ti1*toff)
300         continue
c
c     The factor of omega is needed because xi = sqrt(1/2Mw0)A and
c     D = < A(tau)A+(0) > whereas we measure <xi(tau)xj(0)>.
c     In order to get the right answer for 
c     D0 = -2*w0/(w0**2+wm**2), we need an extra factor of 1/4,  the
c     origin of which eludes me.
c
            dnl(kx,ky,ti) = gtemp * omega / (2.d0*n*n*l)
            if (ti.eq.0) then
              dnl(kx,ky,l) = dnl(kx,ky,0)
            endif
200     continue
100   continue

      do 400 kx = 0, n/2
        do 400 ky = 0, n/2
          do 400 ti = 0,l
            adnl(kx,ky,ti)=adnl(kx,ky,ti)
     @           +dnl(kx,ky,ti)*sgnup*sgndn
            adnl2(kx,ky,ti)=adnl2(kx,ky,ti)
     @           +dnl(kx,ky,ti)*sgnup*sgndn
400   continue
      asgnph = asgnph + sgnup*sgndn
      asgnph2 = asgnph2 + sgnup*sgndn
c      write (88,*) count,sgnup,sgndn,asgnph2

      xeq=0.d0
      do 480 i = 0, toff-1
      do 480 ti = 0,l-1
              xeq=xeq+hub(i+ti*toff)
480   continue
      axeq=axeq+xeq*sgnup*sgndn/dfloat(toff*l)


      count=count+1

      if (mod(count,nbin).eq.0) then
         if (count.ne.0) then
              write (88,*) '   '
              write (88,*) 'dumping adnl2=0 at count=  ',count

c              write (88,*) '   '
c              write (88,*) 'real space bin at ',count
c              write (88,*) '  tau     rx       ry     <AA> '
              do 1100 ti = 0, l-1
              do 1200 kx = 0,n/2
              do 1200 ky = 0,kx
                   atemp=0.5d0*(adnl2(kx,ky,ti)+adnl2(ky,kx,ti))
c                  write (88,1990) ti,kx,ky,atemp/asgnph2
                   adnl2(kx,ky,ti)=atemp
                   adnl2(ky,kx,ti)=atemp
1200          continue
1100          continue
c              write (88,1991) asgnph2
              call ftntok(adnl2,adql2,n,l)
              write (88,*) '   '
              write (88,*) 'momentum space D bin  ',count/nbin
              write (88,*) 'tau      kx       ky     sgn*<Ak Ak> '
              do 1210 kx = 0,n/2
              do 1210 ky = 0,kx
                 do 1110 ti = 0, l-1
                   write (88,1990) ti,kx,ky,adql2(kx,ky,ti)/real(nbin)
1110          continue
1210          continue
              write (88,*) 'sign of bin',count/nbin
              write (88,1991) asgnph2/real(nbin)
1990          format(3i8,f16.8)
1991          format(f12.8)
         endif
         asgnph2=0.d0
      endif




      return
      end
c     
c******************makegt(ti,dti,gtup,gsup,gtdn,gsdn)*************
c This subroutine returns the unequal time Green's functions
c      gtup=Gup(ti+dti,ti)   gsup=Gup(ti,ti+dti)
c      gtdn=Gdn(ti+dti,ti)   gsdn=Gdn(ti,ti+dti)

c In the measurement routines gtup,... should not be changed since unless
c dti=0, it is assumed that they contain Gup(ti+dti-1,ti),...
c Note that for dti=0 gtup and gtdn return the equal time Green's functions,
c while gsup and gsdn return gtup-1 and gtdn-1, which are not needed for
c measurements, but are needed for future calculations of gsup and gsdn.
c Note that dti should range from 0 to l-1.

       subroutine makegt(ti,dti,gtup,gsup,inacc)

c pars.h -- basic Holstein parameters
       use global
        implicit none
c
c getgpam.h -- parameters for getg (2dmoo)
c
        integer orthlen,doauto,torth
        real*8 eorth,difflim,errpam
        common/getgpam/eorth,difflim,errpam,orthlen,doauto,
     1                torth
c vectors.h -- vectors for 2d holstein model.
        real*8 hub(0:volume-1)
        real*8 vup(0:volume-1),vdn(0:volume-1)
        real*8 mphase(0:toff)

        common/vectors/hub,vup,vdn, mphase

       integer ti,dti,i,j,kk,inacc,ii,jj

       real*8 gtup(0:n*n,0:n*n),gsup(0:n*n,0:n*n)
       real*8 gt2(0:n*n,0:n*n),gs2(0:n*n,0:n*n),difft,diffs

c Calculate the Green's functions on the next time slice.
       inacc = 0
       kk=mod(dti,torth)
c Multiply gt on the left by B(ti+dti) and gs on the right by B^-1(ti+dti)
       i=mod(ti+dti-1,l)
c
c Fix so we don't try to multbi random matrices for dti = 0
c
       if (dti.eq.0) then
         do 5600 ii = 0, n*n-1
           do 5610 jj = 0, n*n-1
             gtup(ii,jj) = 0.d0
             gsup(ii,jj) = 0.d0
5610       continue
5600     continue
       else
         call multb(gtup,vup,i,-1)
         call multbi(gsup,vup,i,1)
       endif
c
       if (kk.eq.0) then
c If kk=0 it is time to calculate the Green's functions from scratch.
        call matcop(gtup,gt2)
        call matcop(gsup,gs2)
        call getgtau(vup,ti,dti,gtup,gsup)
c getgtau returns gsup and gsdn with the wrong sign. Correct this.
        do 20 j=0,n*n-1
         do 20 i=0,n*n-1
          gsup(i,j)=-gsup(i,j)
20      continue
        call matdif(gtup,gt2,difft)
        call matdif(gsup,gs2,diffs)
        if(difft .gt. 1.0e-4 .or. diffs .gt. 1.0e-4)then
          if(dti .ne. 0)then
            inacc = 1
            write(6,*)'difft,diffs are ',difft,diffs
            return
          endif
        endif
        if(difft .gt. 1.0e-4 .or. diffs .gt. 1.0e-4)then
          if(dti .ne. 0)then
            inacc = 1
            write(6,*)'difft,diffs are ',difft,diffs
            return
          endif
        endif
       endif

        return
        end
c***********************************
c This getudr includes pivoting.
       subroutine getudr(vvv,ti1,ti2,dodepth,orthlen,omat,bvec,rmat)

c pars.h -- basic Holstein parameters
       use global
        implicit none

       integer i,mx,orthlen,my,kk,ti1,ti2,dodepth
       real*8 vvv(0:volume-1),omat(0:n*n,0:n*n)
       real*8 rmat(0:n*n,0:n*n),rmat3(0:n*n,0:n*n)
       real*8 rmat2(0:n*n,0:n*n),bvec(0:n*n)
       integer depth
c
c Assume ti2 > ti1
        call unit(omat)
        call unit(rmat)
        do 10 mx = 0, n*n-1
10        bvec(mx) = 1.d0
        do 40 i = ti1, ti2
          kk = mod(i+l,l)
          call multb(omat,vvv,kk,-1)
          if(mod(i+1-ti1,orthlen) .eq. 0 .or. i .eq. ti2)then
            do 30 my = 0, n*n-1
              do 30 mx = 0, n*n-1
30              omat(mx,my) = omat(mx,my) * bvec(my)
            if(i .eq. ti2 .or. dodepth .eq. 0)then
              depth = n*n
            else if(mod(i+1-ti1,2*orthlen) .eq. 0)then
              depth = 2*n*n/3
            else
              depth = n*n/3
            endif
            if(n .eq. 2)then
              depth = 4
            endif
            call orthfacp(omat,bvec,rmat2,depth)
            call matmult(rmat2,rmat,rmat3)
            call matcop(rmat3,rmat)
          endif
40      continue
       do 35 my = 0, n*n-1
         do 35 mx = 0, n*n-1
35         omat(mx,my) = omat(mx,my) * bvec(my)
       call orthfacp(omat,bvec,rmat2,n*n)
       call matmult(rmat2,rmat,rmat3)
       call matcop(rmat3,rmat)
       return
       end
c
       subroutine getudri(vvv,ti1,ti2,dodepth,orthlen,omat,bvec,rmat)

c pars.h -- basic Holstein parameters
       use global
        implicit none
      
       integer i,mx,orthlen,my,kk,ti1,ti2,dodepth
       real*8 vvv(0:volume-1),omat(0:n*n,0:n*n)
       real*8 rmat(0:n*n,0:n*n),rmat3(0:n*n,0:n*n)
       real*8 rmat2(0:n*n,0:n*n),bvec(0:n*n)
       integer depth
c
c Assume ti2 > ti1
        call unit(omat)
        call unit(rmat)
        do 10 mx = 0, n*n-1
10        bvec(mx) = 1.d0
        do 40 i = ti2, ti1, -1
          kk = mod(i+l,l)
          call multbi(omat,vvv,kk,-1)
          if(mod(i+1-ti1,orthlen) .eq. 0 .or. i .eq. ti2)then
            do 30 my = 0, n*n-1
              do 30 mx = 0, n*n-1
30              omat(mx,my) = omat(mx,my) * bvec(my)
            if(i .eq. ti2 .or. dodepth .eq. 0)then
              depth = n*n
            else if(mod(i+1-ti1,2*orthlen) .eq. 0)then
              depth = 2*n*n/3
            else
              depth = n*n/3
            endif
            if(n .eq. 2)then
              depth = 4
            endif
            call orthfacp(omat,bvec,rmat2,depth)
            call matmult(rmat2,rmat,rmat3)
            call matcop(rmat3,rmat)
          endif
40      continue
       do 35 my = 0, n*n-1
         do 35 mx = 0, n*n-1
35         omat(mx,my) = omat(mx,my) * bvec(my)
       call orthfacp(omat,bvec,rmat2,n*n)
       call matmult(rmat2,rmat,rmat3)
       call matcop(rmat3,rmat)
       return
       end
c**********************************
c getgp -- routine to include pivoting in the orthogonalization. 1/30/89
       subroutine getgp(vvv,ti,gmat,sgndet,deta)

c pars.h -- basic Holstein parameters
       use global
        implicit none
c
c getgpam.h -- parameters for getg (2dmoo)
c
        integer orthlen,doauto,torth
        real*8 eorth,difflim,errpam
        common/getgpam/eorth,difflim,errpam,orthlen,doauto,
     1                torth

       integer i,j,mx,k,my,kk,ti,depth
       real*8 vvv(0:volume-1)
       real*8 gmat(0:n*n,0:n*n)

       real*8 rmat(0:n*n,0:n*n),rmat2(0:n*n,0:n*n),rmat3(0:n*n,0:n*n)
       real*8 oimat(0:n*n,0:n*n),omat(0:n*n,0:n*n)
       real*8 omat2(0:n*n,0:n*n),omat3(0:n*n,0:n*n)
       common/getgmats/rmat,rmat2,rmat3,oimat,omat,omat2,omat3

       real*8 bvec(0:n*n),sgndet,det,deta
       integer ipvt(toff),nbig

        nbig = 0
        call unit(omat)
        call unit(rmat)
        do 10 mx = 0, toff-1
10        bvec(mx) = 1.d0
        do 40 i = 0, l-1
          kk = mod(i+ti+l,l)
          call multb(omat,vvv,kk,-1)
          if(mod(i+1,orthlen) .eq. 0 .or. i .eq. l-1)then
            do 30 my = 0, toff-1
              do 30 mx = 0, toff-1
30              omat(mx,my) = omat(mx,my) * bvec(my)
            if(i .eq. l-1 )then
              depth = n*n
            else if(mod(i+1,2*orthlen) .eq. 0)then
              depth = 2*n*n/3
            else
              depth = n*n/3
            endif
            if(n .eq. 2)then
              depth = 4
            endif
c
c           depth = toff
            call orthfacp(omat,bvec,rmat2,depth)
            call matmult(rmat2,rmat,rmat3)
            call matcop(rmat3,rmat)
          endif
40      continue
c One last orthogonalization to make sure omat is orthogonal.
        do 35 my = 0, toff-1
          do 35 mx = 0, toff-1
35          omat(mx,my) = omat(mx,my) * bvec(my)
        call orthfacp(omat,bvec,rmat2,depth)
        call matmult(rmat2,rmat,rmat3)
        call matcop(rmat3,rmat)
c
c oimat = omat^-1
        call invertr(omat,oimat,det,toff+1,toff,2,ipvt)
	deta=det
        sgndet = 1.d0
        if(det .lt. 0.d0) sgndet = -1.d0
c rmat2 = rmat^-1
        call invertr(rmat,rmat2,det,toff+1,toff,2,ipvt)
	deta=deta*det

c Calculate omat2 = D + omat^-1 R^-1
c Dimension of omat2 is toff-nbig.
        call zeromat(omat2)
        do 65 i = 0, toff-1
65        omat2(i,i) = bvec(i)
        do 60 k = 0, toff-1
         do 60 j = 0, toff-1
          do 60 i = 0, toff-1
60         omat2(i,j)=omat2(i,j)+oimat(i,k)*rmat2(k,j)
c calculate inverse of omat2, put in omat3
        call invertr(omat2,omat3,det,toff+1,toff-0,2,ipvt)
	deta=deta*det
        if(det .lt. 0.d0) sgndet = -sgndet

        call zeromat(omat)
        do 80 i = 0,toff-1
         do 80 j = 0,toff-1
          do 80 k = 0,toff-1
80         omat(i,k) = omat(i,k) + omat3(i-nbig,j-nbig)*oimat(j,k)

c Multiply to get gmat=rmat2*omat.
c Remember:omat(j,k) = 0 for j < 0
        call zeromat(gmat)
        do 95 i = 0,toff-1
          do 95 j = nbig,toff-1
            do 95 k = 0,toff-1
95            gmat(i,k) = gmat(i,k) + rmat2(i,j)*omat(j,k)

       return
       end

c**********************************
c getgp -- routine to include pivoting in the orthogonalization. 1/30/89
       subroutine debug(vvv,ti,gmat,   det)

c pars.h -- basic Holstein parameters
       use global
        implicit none
c
c getgpam.h -- parameters for getg (2dmoo)
c
        integer orthlen,doauto,torth
        real*8 eorth,difflim,errpam
        common/getgpam/eorth,difflim,errpam,orthlen,doauto,
     1                torth

       integer i,j,mx,k,my,kk,ti,depth
       real*8 vvv(0:volume-1)
       real*8 gmat(0:n*n,0:n*n)

       real*8 rmat(0:n*n,0:n*n),rmat2(0:n*n,0:n*n),rmat3(0:n*n,0:n*n)
       real*8 oimat(0:n*n,0:n*n),omat(0:n*n,0:n*n)
       real*8 omat2(0:n*n,0:n*n),omat3(0:n*n,0:n*n)
       common/getgmats/rmat,rmat2,rmat3,oimat,omat,omat2,omat3

       real*8 bvec(0:n*n),sgndet,det
       integer ipvt(toff),nbig

        nbig = 0
        call unit(omat)
        call unit(rmat)
        do 10 mx = 0, toff-1
10        bvec(mx) = 1.d0
        do 40 i = 0, l-1
          kk = mod(i+ti+l,l)
          call multb(omat,vvv,kk,-1)
          if(mod(i+1,orthlen) .eq. 0 .or. i .eq. l-1)then
            do 30 my = 0, toff-1
              do 30 mx = 0, toff-1
30              omat(mx,my) = omat(mx,my) * bvec(my)
            if(i .eq. l-1 )then
              depth = n*n
            else if(mod(i+1,2*orthlen) .eq. 0)then
              depth = 2*n*n/3
            else
              depth = n*n/3
            endif
            if(n .eq. 2)then
              depth = 4
            endif
c
c           depth = toff
            call orthfacp(omat,bvec,rmat2,depth)
            call matmult(rmat2,rmat,rmat3)
            call matcop(rmat3,rmat)
          endif
40      continue
c One last orthogonalization to make sure omat is orthogonal.
        do 35 my = 0, toff-1
          do 35 mx = 0, toff-1
35          omat(mx,my) = omat(mx,my) * bvec(my)
        call orthfacp(omat,bvec,rmat2,depth)
        call matmult(rmat2,rmat,rmat3)
        call matcop(rmat3,rmat)
c
c oimat = omat^-1
        call invertr(omat,oimat,det,toff+1,toff,2,ipvt)
        sgndet = 1.d0
        if(det .lt. 0.d0) sgndet = -1.d0
c rmat2 = rmat^-1
        call invertr(rmat,rmat2,det,toff+1,toff,2,ipvt)

c Calculate omat2 = D + omat^-1 R^-1
c Dimension of omat2 is toff-nbig.
        call zeromat(omat2)
        do 65 i = 0, toff-1
65        omat2(i,i) = bvec(i)
        do 60 k = 0, toff-1
         do 60 j = 0, toff-1
          do 60 i = 0, toff-1
60         omat2(i,j)=omat2(i,j)+oimat(i,k)*rmat2(k,j)
c calculate inverse of omat2, put in omat3
        call invertr(omat2,omat3,det,toff+1,toff-0,2,ipvt)
        if(det .lt. 0.d0) sgndet = -sgndet

        call zeromat(omat)
        do 80 i = 0,toff-1
         do 80 j = 0,toff-1
          do 80 k = 0,toff-1
80         omat(i,k) = omat(i,k) + omat3(i-nbig,j-nbig)*oimat(j,k)

c Multiply to get gmat=rmat2*omat.
c Remember:omat(j,k) = 0 for j < 0
        call zeromat(gmat)
        do 95 i = 0,toff-1
          do 95 j = nbig,toff-1
            do 95 k = 0,toff-1
95            gmat(i,k) = gmat(i,k) + rmat2(i,j)*omat(j,k)

       return
       end

c****** NEW VERSION ******* (with pivoting)
       subroutine orthogp(mat,rmat,depth)

c pars.h -- basic Holstein parameters
       use global
        implicit none
  
       integer i,j,k,depth, done(0:n*n-1),ib
       real*8 mat(0:n*n,0:n*n),rmat(0:n*n,0:n*n),maxsize,colsize
c
       call unit(rmat)
       do 10 i = 0, n*n-1
10       done(i) = 0
c done(i)=0 means the other columns havent been orthogonalized to col i yet.
       do 20 i = 0, min(depth,n*n)-1
c First find the column which is largest and hasn't been done yet.
c That col. will be "ib".
         maxsize = -10.d0
         ib = -1
         do 30 j = 0, n*n-1
           if(done(j) .eq. 0)then
             colsize = 0.d0
             do 40 k=0, n*n-1
40             colsize = colsize + mat(k,j)**2
             if(colsize .gt. maxsize) then
               maxsize = colsize
               ib = j
             endif
           endif
30       continue
         done(ib) = 1
         call donorm(mat(0,ib),rmat(ib,ib))
c Now orthogonalize all the other columns to ib.
         do 50 j = 0, n*n-1
           if(done(j) .eq. 0) then
             do 60 k=0, n*n-1
60             rmat(ib,j) = rmat(ib,j) + mat(k,ib)*mat(k,j)
             do 70 k=0, n*n-1
70             mat(k,j) = mat(k,j) - rmat(ib,j)*mat(k,ib)
           endif
50       continue

20     continue

       return
       end

c********************************

       subroutine orthfacp(mat,dvec,rmat,depth)

c pars.h -- basic Holstein parameters
       use global
        implicit none

       integer i,j,depth
       real*8 mat(0:n*n,0:n*n),rmat(0:n*n,0:n*n),dvec(0:n*n-1)
       real*8 temp
c
       call orthogp(mat,rmat,depth)
       do 10 i = 0, n*n-1
         dvec(i) = rmat(i,i)
         if(dvec(i) .eq. 0.d0) then
           rmat(i,i) = 1.d0
         else if(dvec(i) .ne. 1.d0) then
             temp = 1.d0 / dvec(i)
             do 20 j = 0, n*n-1
20             rmat(i,j) = rmat(i,j) * temp
         endif
10     continue
       return
       end
c**************************subroutine getgtau*********************
c This subroutine uses Eugene's idea for unequal time Green's functions
c with Steve's factorizing out diagonal matrices trick.
       subroutine getgtau(vvv,ti,dti,gmat,gmat2)

c pars.h -- basic Holstein parameters
       use global
        implicit none
c
c getgpam.h -- parameters for getg (2dmoo)
c
        integer orthlen,doauto,torth
        real*8 eorth,difflim,errpam
        common/getgpam/eorth,difflim,errpam,orthlen,doauto,
     1                torth

       integer i,j,ti,dti,depth
       real*8 vvv(0:volume-1)
       real*8 gmat(0:n*n,0:n*n), gmat2(0:n*n,0:n*n)

       real*8 r1(0:n*n,0:n*n),r1i(0:n*n,0:n*n),r21i(0:n*n,0:n*n)
       real*8 r12i(0:n*n,0:n*n)
       real*8 oimat(0:n*n,0:n*n),u1(0:n*n,0:n*n)
       real*8 u2i1(0:n*n,0:n*n),u2i(0:n*n,0:n*n)
       real*8 u1i2(0:n*n,0:n*n)
       real*8 u2(0:n*n,0:n*n),r2(0:n*n,0:n*n)
       real*8 imat(0:n*n,0:n*n)

       real*8 d1(0:n*n),d2(0:n*n)
       real*8 dbar1i(0:n*n),dbar2i(0:n*n),dtil1(0:n*n),dtil2(0:n*n)
       real*8 det
       integer ipvt(toff)

        depth = n*n
        call getudr(vvv,ti+dti,ti+l-1,0,orthlen,u1,d1,r1)
        call getudri(vvv,ti,ti+dti-1,0,orthlen,u2,d2,r2)
c r1i = r1^-1 ; r21i = r2*r1^-1; r12i = r1*r2^-1
        call matcop(r1,imat)
        call invertr(imat,r1i,det,toff+1,toff,2,ipvt)
        call matmult(r2,r1i,r21i)
        call matcop(r21i,imat)
        call invertr(imat,r12i,det,toff+1,toff,2,ipvt)
c Calculate u2i=u2^-1; u2i1 = u2^-1*u1; u1i2 = u1^-1*u2
        call transp(u2,u2i)
        call matmult(u2i,u1,u2i1)
        call matcop(u2i1,imat)
        call invertr(imat,u1i2,det,toff+1,toff,2,ipvt)
c To make the last inversion stable, factor out some diagonal matrices.
c Define  dbar1i(i) = 1.d0 / max(d1(i),1); dtil1(i)=dbar1i(i)*d1(i)=min(1,d1(i))
c same for dbar2i, etc.
c The formula for G(ti+dti,ti) is:
c  G = r1i dbar1i (dtil2 r21i dbar1i + dbar2i u2i1 dtil1)^-1 dbar2i u2i
c
        do 76 j = 0, n*n-1
          dbar1i(j) = 1.d0 / max(d1(j),1.d0)
          dtil1(j) = d1(j) * dbar1i(j)
          dbar2i(j) = 1.d0 / max(d2(j),1.d0)
76        dtil2(j) = d2(j) * dbar2i(j)
        do 75 j = 0, n*n-1
          do 75 i = 0, n*n-1
            r21i(i,j) = dtil2(i)*r21i(i,j)*dbar1i(j)
75          u2i1(i,j) = dbar2i(i)*u2i1(i,j)*dtil1(j)
c Now add the two matrices and invert
        do 65 i = 0, n*n-1
          do 65 j = 0, n*n-1
            gmat(i,j)=u2i1(i,j)+r21i(i,j)
65      continue
        call invertr(gmat,oimat,det,toff+1,toff,2,ipvt)
        do 85 i = 0, n*n-1
          do 85 j = 0, n*n-1
            oimat(i,j)=oimat(i,j)*dbar1i(i)*dbar2i(j)
85      continue
c  Multiply on the left by r1i and the right by u2i
        call matmult(r1i,oimat,r21i)
        call matmult(r21i,u2i,gmat)
c
c Now do G(ti,ti+dti)
        do 176 j = 0, n*n-1
          d1(j) = 1.d0 / d1(j)
          d2(j) = 1.d0 / d2(j)
          dbar1i(j) = 1.d0 / max(d1(j),1.d0)
          dtil1(j) = d1(j) * dbar1i(j)
          dbar2i(j) = 1.d0 / max(d2(j),1.d0)
176       dtil2(j) = d2(j) * dbar2i(j)
        do 175 j = 0, n*n-1
          do 175 i = 0, n*n-1
            r12i(i,j) = dtil2(j)*r12i(i,j)*dbar1i(i)
175         u1i2(i,j) = dbar2i(j)*u1i2(i,j)*dtil1(i)
        do 165 i = 0, n*n-1
          do 165 j = 0, n*n-1
            gmat2(i,j)=u1i2(i,j)+r12i(i,j)
165      continue
        call invertr(gmat2,oimat,det,toff+1,toff,2,ipvt)
        do 185 i = 0, n*n-1
          do 185 j = 0, n*n-1
            oimat(i,j)=oimat(i,j)*dbar1i(j)*dbar2i(i)
185      continue
c  Multiply on the left by r1i and the right by u2i
        call matmult(u2,oimat,r21i)
        call matmult(r21i,r1,gmat2)

       return
       end

c
c sweep.f -- sweep routines for 2d holstein model.
c
c*******************************************
        subroutine sweep(gmatup,gmatdn,accept,reject,wraps)

c pars.h -- basic Holstein parameters
        use global
        implicit none

c integers.h -- integer parameters for holstein.
        integer warms,sweeps,msr,nwrap
        integer iran
        common/integers/warms,sweeps,msr,nwrap, iran
c
c getgpam.h -- parameters for getg (2dmoo)
c
        integer orthlen,doauto,torth
        real*8 eorth,difflim,errpam
        common/getgpam/eorth,difflim,errpam,orthlen,doauto,
     1                torth
c couple.h -- coupling parameters for holstein.
         real*8 move,t,omega,g,mu,dtau,expmu,gam,lambda,tdtau,dens
         real*8 pe,kee
         common/couple/move,t,omega,g,mu,
     1       dtau,expmu,gam,lambda,tdtau,dens,pe,kee
c kinetic.h -- kinetic energy common block. (2dmoo)
        real*8 ch,soc,ch4
        common/kinetic/ch,soc,ch4
c vectors.h -- vectors for 2d holstein model.
        real*8 hub(0:volume-1)
        real*8 vup(0:volume-1),vdn(0:volume-1)
        real*8 mphase(0:toff)

        common/vectors/hub,vup,vdn, mphase
c m0var.h -- common for measurements. (2dmoo)
        real*8 sgnup,sgndn,sgn,nup,ndn,nud,ke
        real*8 saf,scdw,sferro,sfer2,grfun(0:n/2,0:n/2),safsq,safsq2
        real*8 spinxx(0:n/2,0:n/2),cdw(0:n/2,0:n/2)
        real*8 pairmat(-1:1,-1:1,-1:1,-1:1),errrat,asgnp
        integer nmeas0,redo,noredo,dopair,nmeasp,numpair,dowrite
        real*8 asgnup,asgndn,asgn,anup,andn,anud,ake,asafsq,asafsq2
        real*8 asaf,ascdw,asferro,asfer2,agrfun(0:n/2,0:n/2)
        real*8 aspinxx(0:n/2,0:n/2),acdw(0:n/2,0:n/2)
        real*8 apairmat(-1:1,-1:1,-1:1,-1:1)
        real*8 den(0:n/2,0:n/2,0:1),aden(0:n/2,0:n/2,0:1)
        common/m0var/ sgnup,sgndn,sgn,nup,ndn,nud,ke,safsq,safsq2,
     1       saf,scdw,sferro,sfer2,grfun,
     1       spinxx,cdw,aspinxx,acdw,asafsq,asafsq2,
     2       asgnup,asgndn,asgn,anup,andn,anud,ake,den,aden,
     3       asaf,ascdw,asferro,asfer2,agrfun,errrat,asgnp
     4       ,pairmat,apairmat
        common/m0vari/nmeas0,redo,noredo,dopair,nmeasp,numpair,dowrite

        integer accept,reject,wraps
        integer ti,ti1
        real*8 gmatup(0:n*n,0:n*n),ogmat(0:n*n,0:n*n),diffup
        real*8 gmatdn(0:n*n,0:n*n),diffdn
        real*8 accrat,redorat
        real*8 dett

        integer lastwr,skipair

        do 10 ti = 0, l-1
          call multb(gmatup,vup,ti,-1)
          call multbi(gmatup,vup,ti,1)
          if(wraps .gt. 0) then
            wraps = wraps - 1
          else
            wraps = nwrap
            ti1 = mod(ti+1,l)
            call matcop(gmatup,ogmat)
            call getgp(vup,ti1,gmatup,sgnup,dett)
            call matdif(gmatup,ogmat,diffup)
            diffdn = diffup
            sgndn = sgnup
            if(diffup.gt.difflim .or. diffdn.gt.difflim)then
              redo = redo+1
            else
              noredo=noredo+1
            endif
          endif
          call swpslice(gmatup,gmatdn,sgnup,sgndn,accept,
     1                   reject,ti)
c Now do the tau=0 measurements with the free of charge Green's functions.
c Only do them every third time slice.
          if(mod(ti,12) .eq. 0) then
            call meas0(gmatup,gmatdn,ti)
            if (dowrite.eq.1) then
              write(67,4656) nup+ndn,scdw 
4656          format(f12.3,f12.3)
            endif
          endif
c Measure pair correlations every skipair slices.
c Set skipair so numpair meas. are done every sweep.
          skipair = l/numpair
          if(mod(l,numpair) .ne. 0)skipair = skipair+1
          if(mod(ti,skipair).eq.0) then
            call  measpair(gmatup,gmatdn)
            if (dowrite.eq.1) then
              write(69,5656) pairmat(0,0,0,0)
5656          format(f12.3)
            endif
          endif
10      continue

        if(accept+reject .gt. 10000) then
          accrat = float(accept)/float(accept+reject)
          if(accrat .gt. 0.52d0 .or. accrat .lt. 0.48d0)then
            gam = gam + (accrat - 0.5d0)
            gam = max(0.d0,gam)
            gam = min(1.d0,gam)
c            write(66,*)'accrat is ',accrat,' changing gam to ',gam
            accept = 100*accrat
            reject = 100*(1.d0-accrat)
          endif
        endif

        if(redo .gt. 20) then
          redorat = float(redo)/float(redo+noredo)
          if(redorat.gt. errrat)then
            if (nwrap.gt.0) then
              nwrap = nwrap - 1
            endif
            write(66,*)'reducing nwrap to ',nwrap
            redo = 0
            noredo = 1
          endif
        endif

        if(noredo .gt. 500) then
          redorat = float(redo)/float(redo+noredo)
          if(redorat .lt. 0.2*errrat) then
            if(nwrap .ge. lastwr)then
c Only increase wrap if nwrap has not been reduced since last increase.
              nwrap = nwrap + 2
              write(66,*)'increasing nwrap to ',nwrap
              redo = 0
              noredo = 1
              lastwr = nwrap
            endif
          endif
        endif

        return
        end
c*******************************************
        subroutine swpslice(gmatup,gmatdn,sgnup,sgndn,accept,
     1                   reject,ti)

c pars.h -- basic Holstein parameters
        use global
        implicit none

c integers.h -- integer parameters for holstein.
        integer warms,sweeps,msr,nwrap
        integer iran
        common/integers/warms,sweeps,msr,nwrap, iran
c couple.h -- coupling parameters for holstein.
         real*8 move,t,omega,g,mu,dtau,expmu,gam,lambda,tdtau,dens
         real*8 pe,kee
         common/couple/move,t,omega,g,mu,
     1       dtau,expmu,gam,lambda,tdtau,dens,pe,kee
c vectors.h -- vectors for 2d holstein model.
        real*8 hub(0:volume-1)
        real*8 vup(0:volume-1),vdn(0:volume-1)
        real*8 mphase(0:toff)

        common/vectors/hub,vup,vdn, mphase

        integer accept,reject
        integer ti,nbar,nbari,j1,j2
        real*8 gmatup(0:n*n,0:n*n)
        real*8 gmatdn(0:n*n,0:n*n),sgnup,sgndn
        real*8 bvecup(0:n*n-1),bvecup2(0:n*n-1)
        real*8 vnbarup,vnbardn,ratup,ratdn,ran,rat,p
        real*8 ran2

c        real*8 dbg(0:n*n,0:n*n),detdbg,aaa1,aaa2,dif,maxdif
        real*8 diff1p,diff2p,dhub,hubp,diff1,diff2,bose,de

        do 20 nbar = 0, toff-1
          nbari = nbar + ti * toff

          ran =ran2(iran)
          dhub=move*(ran-0.5d0)

          vnbarup = dexp(lambda*dhub)
          vnbarup = vnbarup - 1.d0
          vnbardn = vnbarup
          ratup = 1.d0 + (1.d0 - gmatup(nbar,nbar))*vnbarup
          ratdn = ratup

          ran = ran2(iran)

          hubp=hub(nbari)+dhub

          de=pe*(hubp*hubp-hub(nbari)*hub(nbari))

          if (ti.eq.0) then
            diff1 =hub(nbari) -hub(nbar+(L-1)*toff)
            diff2 =hub(nbari) -hub(nbari+toff)
            diff1p=hubp-hub(nbar+(L-1)*toff)
            diff2p=hubp-hub(nbari+toff)
          else if (ti.eq.L-1) then
            diff1 =hub(nbari) -hub(nbari-toff)
            diff2 =hub(nbari) -hub(nbar)
            diff1p=hubp-hub(nbari-toff)
            diff2p=hubp-hub(nbar)
          else
            diff1 =hub(nbari) -hub(nbari-toff)
            diff2 =hub(nbari) -hub(nbari+toff)
            diff1p=hubp-hub(nbari-toff)
            diff2p=hubp-hub(nbari+toff)
          endif
          de=de+kee*(diff1p*diff1p+diff2p*diff2p-
     1              diff1 *diff1 -diff2 *diff2  )
          bose=dexp(-de)
          rat = abs(ratup*ratdn)*bose
          if(rat .le. 1.d0) p = rat/(1.d0+gam*rat)
          if(rat .gt. 1.d0) p = rat/(gam+rat)

          if(p .gt. ran) then

            accept = accept + 1
            if(ratup .lt. 0.d0)sgnup = -sgnup
            if(ratdn .lt. 0.d0)sgndn = -sgndn

              do 35 j1 = 0, toff-1
35              bvecup(j1) = -vnbarup*gmatup(nbar,j1)
              bvecup(nbar) = bvecup(nbar) + vnbarup
              do 45 j1 = 0, toff-1
45              bvecup2(j1) = gmatup(j1,nbar)/(1.d0+bvecup(nbar))
              do 55 j1 = 0, toff-1
               do 55 j2 = 0, toff-1
55              gmatup(j1,j2)=gmatup(j1,j2)-bvecup2(j1)*bvecup(j2)

c *******DEBUG********
c                  ti1=mod(ti+1,l)
c                  call debug(vup,ti1,dbg,detdbg)
c                  aaa1=detdbg
c                  call debug(vup,ti1,dbg,detdbg)
c                  aaa2=detdbg
c *******ENDDEBUG*********

            hub(nbari) = hubp
            vup(nbari) = vup(nbari) * (vnbarup+1.d0)
            vdn(nbari) = vdn(nbari) * (vnbardn+1.d0)

c *******DEBUG*********
c                  call debug(vup,ti1,dbg,detdbg)
c                  maxdif=0.d0
c                  do 9000 ixx=0,toff-1
c                  do 9000 iyy=0,toff-1
c                  dif=abs(dbg(ixx,iyy)-gmatup(ixx,iyy))
c9000              if (dif.gt.maxdif) maxdif=dif
c                  write (6,*) maxdif,detdbg/aaa1/ratup
c *******ENDDEBUG*********

          else
            reject = reject + 1
          endif
20      continue

        return
        end

c
c sweep2.f -- sweep2 routines for 2d holstein model.
c
c*******************************************
        subroutine sweep2(gmatup,gmatdn,accept,
     1        reject,wraps,numtry,gsize)

c pars.h -- basic Holstein parameters
        use global
        implicit none

c integers.h -- integer parameters for holstein.
        integer warms,sweeps,msr,nwrap
        integer iran
        integer itry,numtry,isite
	real*8 gsize,msize,de,r
        real*8 ran2
        common/integers/warms,sweeps,msr,nwrap, iran
c
c getgpam.h -- parameters for getg (2dmoo)
c
        integer orthlen,doauto,torth
        real*8 eorth,difflim,errpam
        common/getgpam/eorth,difflim,errpam,orthlen,doauto,
     1                torth
c couple.h -- coupling parameters for holstein.
         real*8 move,t,omega,g,mu,dtau,expmu,gam,lambda,tdtau,dens
         real*8 pe,kee
         common/couple/move,t,omega,g,mu,
     1       dtau,expmu,gam,lambda,tdtau,dens,pe,kee
c kinetic.h -- kinetic energy common block. (2dmoo)
        real*8 ch,soc,ch4
        common/kinetic/ch,soc,ch4
c vectors.h -- vectors for 2d holstein model.
        real*8 hub(0:volume-1)
        real*8 vup(0:volume-1),vdn(0:volume-1)
        real*8 mphase(0:toff)

        common/vectors/hub,vup,vdn, mphase
c m0var.h -- common for measurements. (2dmoo)
        real*8 sgnup,sgndn,sgn,nup,ndn,nud,ke
        real*8 saf,scdw,sferro,sfer2,grfun(0:n/2,0:n/2),safsq,safsq2
        real*8 spinxx(0:n/2,0:n/2),cdw(0:n/2,0:n/2)
        real*8 pairmat(-1:1,-1:1,-1:1,-1:1),errrat,asgnp
        integer nmeas0,redo,noredo,dopair,nmeasp,numpair,dowrite
        real*8 asgnup,asgndn,asgn,anup,andn,anud,ake,asafsq,asafsq2
        real*8 asaf,ascdw,asferro,asfer2,agrfun(0:n/2,0:n/2)
        real*8 aspinxx(0:n/2,0:n/2),acdw(0:n/2,0:n/2)
        real*8 apairmat(-1:1,-1:1,-1:1,-1:1)
        real*8 den(0:n/2,0:n/2,0:1),aden(0:n/2,0:n/2,0:1)
        common/m0var/ sgnup,sgndn,sgn,nup,ndn,nud,ke,safsq,safsq2,
     1       saf,scdw,sferro,sfer2,grfun,
     1       spinxx,cdw,aspinxx,acdw,asafsq,asafsq2,
     2       asgnup,asgndn,asgn,anup,andn,anud,ake,den,aden,
     3       asaf,ascdw,asferro,asfer2,agrfun,errrat,asgnp
     4       ,pairmat,apairmat
        common/m0vari/nmeas0,redo,noredo,dopair,nmeasp,numpair,dowrite

        integer accept,reject,wraps
        integer ti
        real*8 gmatup(0:n*n,0:n*n),gmatdn(0:n*n,0:n*n)
        real*8 gmatupp(0:n*n,0:n*n),detupp,sgnupp,detup
        real*8 detdnp,sgndnp,detdn

c	TRY MOVING A WHOLE LINE OF PHONON COORDINATES (IN TIME)
c	ON numtry RANDOMLY CHOSEN SITE(S)

	if (numtry.ne.0) then
        call getgp(vup,0,gmatup,sgnup,detup)
	sgndn=sgnup
        detdn=detup
	sgn=sgnup*sgndn
	endif

	do 1000 itry=1,numtry
	    isite=int(ran2(iran)*toff)
	    msize=gsize*(ran2(iran)-0.5d0)
	    do 100 ti=0,L-1
	        hub(isite+ti*toff)=hub(isite+ti*toff)+msize
100	    continue
            call setvup()

c	GET THE NEW DETERMINANTS

            call getgp(vup,0,gmatupp,sgnupp,detupp)
	    sgndnp=sgnupp
            detdnp=detupp

c	GET THE CHANGE IN THE PHONON ACTION (NO CHANGE IN KE!)

	    de=0.
	    do 102 ti=0,L-1
                de=de+2.d0*hub(isite+ti*toff)*msize-msize*msize
102	    continue
	    de=de*pe

	    r=ran2(iran)
	    if (r.le.dexp(-de)*abs(detupp*detdnp
     1                            /detup/detdn) ) then
               accept=accept+1
               call matcop(gmatupp,gmatup)
	       sgnup=sgnupp
	       sgndn=sgnupp
	       detup=detupp
	       detdn=detupp
	       sgn=sgnupp*sgnupp
	       else
	       do 200 ti=0,L-1
	           hub(isite+ti*toff)=hub(isite+ti*toff)-msize
200	       continue
	       call setvup()
               reject=reject+1
	    endif
1000    continue

        return
        end



        subroutine invertr(b,c,d,dim,toff,idet,ipvt)   
                   
c       This subroutine finds the inverse of the matrix b
c       using Gaussian elimination with partial pivoting.
c       b is the input matrix.
c       c is the inverse of b.
c       d is the determinent of b.
c       toff is the size of the matrix
c       If idet=0 the inverse of b is returned but not the determinent.
c       If idet=1 the determinent is returned,but not the inverse.
c       If idet=2 both the inverse and the determinent are returned.
c                                     
        implicit none
        integer dim,toff,ipvt(dim) 
        real*8 b(dim,dim),c(dim,dim),d,tb            
        integer iwn,nm1,k,kp1,mb,i,j,idet,kb,km1
 
        ipvt(toff)=1
       iwn=1
       nm1=toff-1
c                                   
c       Finding the pivot
c
       do 35 k=1,nm1
       kp1=k+1
       mb=k
       do 15 i=kp1,toff
       if (abs(b(i,k)).gt.abs(b(mb,k))) mb=i
15       continue
       ipvt(k)=mb                                          
       tb=b(mb,k)
       if (mb.ne.k) then
       ipvt(toff)=-ipvt(toff)
       b(mb,k)=b(k,k)
       b(k,k)=tb
       endif
       if (abs(tb).lt.1.d-10) then
       iwn=0
       go to 35
       endif
c
c       Computing multipliers
c
       do 20 i=kp1,toff
20       b(i,k)=-b(i,k)/tb
c
c       Interchange and elimination by columns
c
       do 30 j=kp1,toff
       tb=b(mb,j)
       if (mb.ne.k) then
       b(mb,j)=b(k,j)
        b(k,j)=tb
       endif
       do 25 i=kp1,toff
25       b(i,j)=b(i,j)+b(i,k)*tb
30       continue
35       continue
       if (iwn.eq.0) then
       write (6,*) 'Warning: Pivot=0'
       d=0.d0    
       go to 150
       endif
       if (idet.gt.0) then
c      
c       Compute the determinent
c
       d=ipvt(toff)
       do 40 i=1,toff
40       d=d*b(i,i)
       endif
       if (idet.ne.1) then
c
c       Compute the inverse
c       Ax=b goes to A~x=Lb~
c
       do 50 i=1,toff
       do 45 j=1,toff
45       c(i,j)=0.d0
50       c(i,i)=1.d0
       do 100 j=1,toff
       do 60 k=1,nm1
       kp1=k+1
       mb=ipvt(k)
       tb=c(mb,j)
       if (mb.ne.k) then
       c(mb,j)=c(k,j)
       c(k,j)=tb
       endif
       do 55 i=kp1,toff
55       c(i,j)=c(i,j)+b(i,k)*tb
60       continue
c
c       Inverting A~
c
       do 70 kb=1,nm1
       km1=toff-kb
       k=km1+1
       c(k,j)=c(k,j)/b(k,k)
       tb=-c(k,j)
       do 65 i=1,km1    
65       c(i,j)=c(i,j)+b(i,k)*tb
70       continue
       c(1,j)=c(1,j)/b(1,1)
100       continue
       endif
150       continue
       return
       end
     
       


      DOUBLE PRECISION FUNCTION RAN2(IDUM)
      implicit none
      double precision rm
      integer ia,ic,idum,iff,ir,iy,j,m
      save 
      PARAMETER (M=714025,IA=1366,IC=150889,RM=1.4005112d-6)
      DIMENSION IR(97)
      DATA IFF /0/
      IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
        IFF=1
        IDUM=MOD(IC-IDUM,M) 
        DO 11 J=1,97
          IDUM=MOD(IA*IDUM+IC,M)
          IR(J)=IDUM
11      CONTINUE
        IDUM=MOD(IA*IDUM+IC,M)
        IY=IDUM
      ENDIF
      J=1+(97*IY)/M
      IF(J.GT.97.OR.J.LT.1)PAUSE
      IY=IR(J)
      RAN2=IY*RM
      IDUM=MOD(IA*IDUM+IC,M)
      IR(J)=IDUM
      RETURN
      END



