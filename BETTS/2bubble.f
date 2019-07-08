        program default
c******************************************************************************
c       This code uses the coarse-grained A(K,w) to generate a two-particle 
c       default model for the spin and charge dynamic susceptibilities.  
c
c	To compile type make 2default (or make all) using the 
c       associated Makefile.  
c
c	This code was designed to read the majority of its inputs from a qmc 
c	sigma file.  
c******************************************************************************
c******************************************************************************
        use Global
        use module_spectra
        implicit none
c******************************************************************************
c       some local declarations 
        integer :: filenum, istat, i, j, n, run, meas, ick, jck, ios,
     &             ik, isk, nkpoints, info, infot, ikx, iky, deg, nfb,
     &             iik, iw, jw, isus
        real(kind) :: r1, r2, r3, r4, r5, ri1, ri2, ri3, ri4, ri5, 
     &                kx, ky, kx1, ky1, kx2, ky2,
     &                kqx, kqy, Exy_k, Exy_kq, tdk,
     &                offset, res, hall, vx, vyy, fermi, delta, y

        
        complex(kind) :: c1, c2,
     &                   zKQwpp, zKwpm, zKwmp, zKQwpm, zKwmm
        
        character*72 linemc,fmt1,string1
c******************************************************************************


c******************************************************************************
c       Read in some parameters
c******************************************************************************
c       Open a datafile for the parameters
        write(6,"('       enter the QMC sigma filenumber:  ',$)") 
        read(5,*) filenum


        if(filenum.lt.10) then  ! P. 402, Chapman
          fmt1="('sigma',i1,'.dat')"
        else if(filenum.lt.100) then
          fmt1="('sigma',i2,'.dat')"
        else if(filenum.lt.1000) then
          fmt1="('sigma',i3,'.dat')"
        else if(filenum.lt.10000) then
          fmt1="('sigma',i4,'.dat')"
        else
          write(6,*) 'filenum to large'
          stop
        end if  
        write(string1,fmt1) filenum
                
c       Now open the QMC self energy sigma.dat
        open(unit=20,file=string1,status='old')

 1000   format(a72)
c       Read some parameters from the QMC sigma file
 401    read(20,1000) linemc
        if(index(linemc(1:72),'Results').ne.0 ) then
          read(20,*) nl,nwn,fill,cluster,ntot,ndim
        else
         goto 401
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



402    read(20,"(a72)") linemc
        if(index(linemc(1:72),'tprime').ne.0 ) then
           if(index(linemc(1:72),'tdprime').ne.0 ) then
              if(index(linemc(1:72),'t3prime').ne.0 ) then
                 read(20,*) ed,U,beta,i,run,meas,tprime,tdprime,t3prime
              else
                 read(20,*) ed,U,beta,i,run,meas,tprime,tdprime
                 t3prime=0.d0
              endif
           else 
              read(20,*) ed,U,beta,i,run,meas,tprime
              tdprime=0.d0
              t3prime=0.d0
           endif
        else
          goto 402
        endif
        temp=1.0/beta   





 

c******************************************************************************
c       This is all the information we need to allocate the arrays
c******************************************************************************
        call allocate_1p
        nf=1000      !set preliminary value just for allocation
        call allocate_Spec0

c******************************************************************************
c       calculate the symmetry tables
c******************************************************************************
        call tables      
        write(6,*)'after tables'  

c*******************************************************************************
c       Write out the K information
c*******************************************************************************       
        if(ndim.eq.2) then
          write(6,*) ' ick  Kcx(ick)   Kcy(ick)'
          do ick=1,Nc
            write(6,"(1x,i3,1x,f9.6,2x,f9.6)") ick,Kc(1,ick),Kc(2,ick)
          end do
        else if(ndim.eq.3) then
          write(6,*) ' ick  Kcx(ick)   Kcy(ick)   Kcz(ick)'
          do ick=1,Nc
             write(6,"(1x,i3,1x,f9.6,2x,f9.6,2x,f9.6)") 
     &		   ick,Kc(1,ick),Kc(2,ick),Kc(3,ick)
          end do
        end if
c******************************************************************************
c       Open the MEM sigma(K,w) files, readin the self-energy, and then spline it.
c******************************************************************************   
        do ick=1,Ncw  ! do all K in the irreducible wedge
c         Open the self energy files.
          if(ick.lt.10) then
            if(filenum.lt.10) then
              fmt1="('Sigma',i1,'_',i1,'.dat')"
            else if(filenum.lt.100) then
              fmt1="('Sigma',i2,'_',i1,'.dat')"
            else if(filenum.lt.1000) then
              fmt1="('Sigma',i3,'_',i1,'.dat')"
            else if(filenum.lt.10000) then
              fmt1="('Sigma',i4,'_',i1,'.dat')"
            end if
          else if(ick.lt.100) then
            if(filenum.lt.10) then
              fmt1="('Sigma',i1,'_',i2,'.dat')"
            else if(filenum.lt.100) then
              fmt1="('Sigma',i2,'_',i2,'.dat')"
            else if(filenum.lt.1000) then
              fmt1="('Sigma',i3,'_',i2,'.dat')"
            else if(filenum.lt.10000) then
              fmt1="('Sigma',i4,'_',i2,'.dat')"
            end if
          else
            write(6,*) 'ick too large'
            stop
          end if
          write(string1,fmt1) filenum,ick
          open(unit=21,file=string1,status='old')
          
c         Read in the real-frequency self energy
          read(21,"(a72)") linemc    
          do i=1,1000
            read(21,*,IOSTAT=istat) w(i),r1,r2,r3,dw(i)
            if(istat.ne.0) goto 998
            Sigma(i,ick)=cmplx(r1,r2)
          end do
 998      continue
          if(ick.eq.1) then
            nf=i-1
          else if(nf.ne.i-1) then       ! trap for errors
            write(6,*) 'Sigma files are not compatible'
            stop
          end if
          
          do i=1,nf
          ReSigma(i,ick)= real(Sigma(i,ick))
          ImSigma(i,ick)=aimag(Sigma(i,ick))
          end do
          call spline(w,ReSigma(:,ick),nf,0.0_8,0.0_8,scReSigma(:,ick))
          call spline(w,ImSigma(:,ick),nf,0.0_8,0.0_8,scImSigma(:,ick))
        end do 
c       Make sure that nf is less than the limit set above        
        if(nf.gt.1000) then
          write(6,*) 'nf in the data file is too large'
          stop
        end if
        
c       Now fill in Sigma(w,k), for k outside the IW.
        do ick=Ncw+1,Nc
        do i=1,nf
          Sigma(i,ick)=Sigma(i,ickmap(ick))
          ReSigma(i,ick)= real(Sigma(i,ick))
          ImSigma(i,ick)=aimag(Sigma(i,ick))
        end do
         call spline(w,ReSigma(:,ick),nf,0.0_8,0.0_8,scReSigma(:,ick))
         call spline(w,ImSigma(:,ick),nf,0.0_8,0.0_8,scImSigma(:,ick))
        end do
 
        write(6,*)'after read in sigma(K,w) and splint it, nf=',nf
 
c******************************************************************************
c       focus on Q=(0,0), ick=icK0=2 for BETTS 2D 8A
        ick=icK0

c******************************************************************************
c       Now generate an inhomogeneous frequency grid (w):
c******************************************************************************
        delta=0.7 
        nfb=300     
        write(6,*)'w(nf)=',w(nf)
        do i=1,nfb
          y=(i-0.5)/float(2*nfb)  !positive freqs. only.
          wb(i)=delta*tan(pi*y)
          if(wb(i).gt.w(nf)) exit
          dwb(i)=pi*(delta**2+wb(i)**2)/(delta*2*nfb)   
        end do 

        write(6,*)'after loop i=1,nfb , i=',i
        nfb=i-1      

c******************************************************************************
c       Also generate an inhomogeneous frequency grid of (w+v): (not code in yet)
c******************************************************************************

c*****************************************************************************
c     and calculate the
c       polarization bubble on this grid (both real and imag. parts).
c*****************************************************************************


c================================================================================================
c \bar(\chi^0) approximation
c================================================================================================
c================================================================================================
c           Now form the inverse unperturbed coarse grained susceptibility
c           matrices
c                       ||
c                       \/
c           ____-1     ____o-1        c-1       co-1
c           CHIM    =  CHIM    +  CHIM   -  CHIM
c               q          q          q         q
c
c           in temat1 (spin, charge, or pairing).  Also calculate the bare
c           reciprocal-space  coarse-grained greens function.
c           as well as the bare extended d-wave bubbles.
c                _______
c           td_k/       \td_k  in Chixd0  where td_k= cos(kx)-cos(ky)
c               \_______/
c
c           Note that for all these susceptibilies, q is defined on the real
c           lattice so that it is (Pi,Pi) if isus=2 even if Nc=1.
c
c=======================================================================
c
c                             ---->---- k    w'-w/2
c       In the ph channel
c                    ---      ----<---- k+q  w'+w/2
c               2 pi \   /
c chiph2(Q,w) = ---- |   | dw' A(K+Q,w'+w/2) A(K,w'-w/2) [f(w'-w/2)-f(w'+w/2)]
c                Nc  /   /
c                    ---
c                     K
c
c In the pp channel: chipp2(Q,w)=0.5*(chi^l+chi^r)
c
c             chi^l                +     chi^r
c          ----<---- -K  -w'+w/2      ---->---- -K  -w'-w/2
c                                  +
c          ----<---- K+Q  w'+w/2      ---->---- K+Q  w'-w/2
c
c                    ---      
c                -pi \    /
c chi^l(Q,w)=   ---- |    | dw' A(K+Q,w'+w/2) A(K,-w'+w/2) [f(-w'+w/2)+f(w'+w/2)-1]
c                Nc  /    /
c                    ---
c                     K
c
c                    ---      
c                pi  \    /
c chi^r(Q,w)=   ---- |    | dw' A(K+Q,w'-w/2) A(K,-w'-w/2) [f(-w'-w/2)+f(w'-w/2)-1]
c                Nc  /    /
c                    ---
c                     K
c==========================================================================
              write(6,*)'before \bar(\chi^0) approach,Nc,ick=',Nc,ick
              write(6,*)'293 iik,ick,ickplus(iik,ick),ickmap(iik)'
              do iik=1,Nc
               write(6,*)iik,ick,ickplus(iik,ick),ickmap(iik)
              enddo

              chipp2=0.0
              chipp2_d=0.0
              chipp2_1=0.0

c             nwsc=1  ! if set nwsc=1, \bar(\chi^0) should have the same result as \chi^0_c
              isus=1 !at the zone corner(pi,pi), isus=2 and at the center(0,0) isus=1

              do iik=1,Nc                   !form pairing CG bare susc d-wave has lower symmetry
                                            !than the lattice so this calc cannot be limited to te IW
              
               isk=ickmap(ickplus(iik,ick)) ! ick=icK0=0, focus on Q=0 case (not the case: ick=icKpi on Q=(pi,pi) isus=2 case.)
               ik=ickmap(iik)               ! map iik to 1BZ IRW for sigma matrix 
               write(6,*)'310 iik,isk,ik=',iik,isk,ik

               do iw=1,nfb !w 
c               write(6,*)'314 iw=',iw             
               do jw=1,nf  !w' sum over frequency
                  

                 call splint(w,ReSigma(:,isk),scReSigma(:,isk),nf,w(jw)+0.5*wb(iw),ri1) !Sig(K+Q,w'+w/2)
                 call splint(w,ImSigma(:,isk),scImSigma(:,isk),nf,w(jw)+0.5*wb(iw),r1)  !Sig(K+Q,w'+w/2)

c                 call splint(w,ReSigma(:,ik),scReSigma(:,ik),nf,w(jw)-0.5*wb(iw),ri2)   !Sig(K,w'-w/2)
c                 call splint(w,ImSigma(:,ik),scImSigma(:,ik),nf,w(jw)-0.5*wb(iw),r2)   !Sig(K,w'-w/2)

                 call splint(w,ReSigma(:,ik),scReSigma(:,ik),nf,-w(jw)+0.5*wb(iw),ri3)  !Sig(K,-w'+w/2)
                 call splint(w,ImSigma(:,ik),scImSigma(:,ik),nf,-w(jw)+0.5*wb(iw),r3)  !Sig(K,-w'+w/2)

c                 for calculation of Im\chi^{0r}(Q,w)
                 call splint(w,ReSigma(:,isk),scReSigma(:,isk),nf,w(jw)-0.5*wb(iw),ri4)   !Sig(K+Q,w'-w/2)
                 call splint(w,ImSigma(:,isk),scImSigma(:,isk),nf,w(jw)-0.5*wb(iw),r4)   !Sig(K+Q,w'-w/2)

                 call splint(w,ReSigma(:,ik),scReSigma(:,ik),nf,-w(jw)-0.5*wb(iw),ri5)  !Sig(K,-w'-w/2)
                 call splint(w,ImSigma(:,ik),scImSigma(:,ik),nf,-w(jw)-0.5*wb(iw),r5)  !Sig(K,-w'-w/2)

c
c Form the w-ed-Sigma(K,w) on the nearest cluster point (Kx,Ky)=Kc(:,ick) to (kx,ky)
c
                  zKQwpp= w(jw)+0.5*wb(iw)-ed-ri1-ii*r1    ! K+Q, w'+w/2, r1
                  zKwpm = w(jw)-0.5*wb(iw)-ed-ri2-ii*r2    ! K,   w'-w/2, r2
                  zKwmp =-w(jw)+0.5*wb(iw)-ed-ri3-ii*r3    ! K,  -w'+w/2, r3
                  zKQwpm= w(jw)-0.5*wb(iw)-ed-ri4-ii*r4    ! K+Q, w'-w/2, r4
                  zKwmm =-w(jw)-0.5*wb(iw)-ed-ri5-ii*r5    ! K,  -w'-w/2, r5
c
                do i=1,nwsc
                  kx=Kc(1,iik)+kt(1,i)
                  ky=Kc(2,iik)+kt(2,i)
                  Exy_k= -0.5d0*(dcos(-kx)+dcos(-ky))
     &                   -tprime*(dcos(-kx)*dcos(-ky)-1.0d0) 
                  
                
                  kqx = kx+pi*dble(isus-1)      ! (k+Q)_x
                  kqy = ky+pi*dble(isus-1)      ! (k+Q)_y
                  Exy_kq = -0.5d0*(dcos(kqx)+dcos(kqy))
     &                    -tprime*(dcos(kqx)*dcos(kqy)-1.0d0)    
c
c                   Particle-Particle
                    tdk = 0.5d0*(dcos(kx)-dcos(ky))   ! lattice k pair form factor
c    need to consider as: chipp=0.5*(chipp^l+ chipp^r)
          chipp2(iw,ick)=chipp2(iw,ick)+dw(jw)*(       
c     &     -r1*r3
     &     -aimag(1.0d0/(zKQwpp-Exy_k))*aimag(1.0d0/(zKwmp-Exy_kq))
     &     *(fermi(-w(jw)+0.5*wb(iw),beta)+fermi(w(jw)+0.5*wb(iw),beta)-1.0)   !chi^l
c     &     +r4*r5
     &     +aimag(1.0d0/(zKQwpm-Exy_k))*aimag(1.0d0/(zKwmm-Exy_kq))
     &     *(fermi(-w(jw)-0.5*wb(iw),beta)+fermi(w(jw)-0.5*wb(iw),beta)-1.0) ) !chi^r

c          chipp2_1(iw,ick)=chipp2_1(iw,ick)+dw(jw)*tdk*(
cc     &     -r1*r3
c     &     -aimag(1.0d0/(zKQwpp-Exy_k))*aimag(1.0d0/(zKwmp-Exy_kq))
c     &     *(fermi(-w(jw)+0.5*wb(iw),beta)+fermi(w(jw)+0.5*wb(iw),beta)-1.0)   !chi^l
cc     &     +r4*r5
c     &     +aimag(1.0d0/(zKQwpm-Exy_k))*aimag(1.0d0/(zKwmm-Exy_kq))
c     &     *(fermi(-w(jw)-0.5*wb(iw),beta)+fermi(w(jw)-0.5*wb(iw),beta)-1.0) ) !chi^r

          chipp2_d(iw,ick)=chipp2_d(iw,ick)+dw(jw)*tdk*tdk*(
c     &     -r1*r3
     &     -aimag(1.0d0/(zKQwpp-Exy_k))*aimag(1.0d0/(zKwmp-Exy_kq))
     &     *(fermi(-w(jw)+0.5*wb(iw),beta)+fermi(w(jw)+0.5*wb(iw),beta)-1.0)   !chi^l
c     &     +r4*r5
     &     +aimag(1.0d0/(zKQwpm-Exy_k))*aimag(1.0d0/(zKwmm-Exy_kq))
     &     *(fermi(-w(jw)-0.5*wb(iw),beta)+fermi(w(jw)-0.5*wb(iw),beta)-1.0) ) !chi^r  
                end do

               end do
               end do
              end do
        chipp2=pi*chipp2/float(Nc*nwsc)/pi/pi           ! normalize       
        chipp2_1=pi*chipp2_1/float(Nc*nwsc)/pi/pi       ! normalize
        chipp2_d=pi*chipp2_d/float(Nc*nwsc)/pi/pi       ! normalize

c       Check that chi2/w is postive (i.e. like the conductivity)
c        do ick=1,Ncw  ! ick fixes to icK0, doesn't loop
          if(chipp2(1,ick).lt.0.0.or.chipp2(nfb/2,ick).lt.0.0) then
            write(6,*) 'sign error in chipp2(i,ick),ick=',ick
            stop
          end if
c        end do ! ick fixes to icK0, doesn't loop

c================================================================================================
c================================================================================================


c******************************************************************************
c   write out dressed-single-particle-green-function-bubble d-wave projected pairing suscetptibility:
c******************************************************************************
c       Form the normalized default model for -chi2(w)/w
        chi2=chipp2_d
c        do ick=1,Ncw  ! ick fixes to icK0, doesn't loop
        r1=0.0  
        do i=1,nfb
          chi2(i,ick)=chi2(i,ick)/(pi*wb(i))
          r1=r1+dwb(i)*chi2(i,ick)
        end do
        chi2(:,ick)=0.5*chi2(:,ick)/r1
c        chi2(:,ick)=0.5*chi2(:,ick)
c        end do  ! ick fixes to icK0, doesn't loop

c        fort.117 contains: T,  chi'_ret^d(w=0,T)
        write(117,*)oner/beta,r1

c        do ick=1,Ncw  ! do all K in the irreducible wedge ! ick fixes to icK0, doesn't loop
c         Open some files.
          if(ick.lt.10) then
            fmt1="('chbrdk',i1)"
          else if(ick.lt.100) then
            fmt1="('chbrdk',i2)"
          else
            write(6,*) 'ick too large'
            stop
          end if
            
          write(string1,fmt1)        ick        !open model_chbk$ick
          open(unit=21,file=string1,status='unknown')

          do i=1,nfb
            write(21,*) wb(i),chi2(i,ick),dwb(i)
          end do
          close(21)
c        end do   ! ick fixes to icK0, doesn't loop
        
c        open(unit=21,file='chbrdr0',status='unknown')
c        do i=1,nfb
c        r1=0.0
c        do ick=1,Ncw
c          r1=r1+ickdeg(ick)*chi2(i,ick)
c        end do
c        write(21,*) wb(i),r1/float(Nc),dwb(i)
c        end do  

c******************************************************************************
c   write out dressed-single-particle-green-function-bubble pairing suscetptibility (not any form factors applied):
c******************************************************************************
c       Form the normalized default model for -chi2(w)/w
        chi2=chipp2
c        do ick=1,Ncw  ! ick fixes to icK0, doesn't loop
        r1=0.0  
        do i=1,nfb
          chi2(i,ick)=chi2(i,ick)/(pi*wb(i))
          r1=r1+dwb(i)*chi2(i,ick)
        end do
        chi2(:,ick)=0.5*chi2(:,ick)/r1
c        chi2(:,ick)=0.5*chi2(:,ick)
c        end do  ! ick fixes to icK0, doesn't loop

c        fort.128 contains: T,  chi'_ret^b(w=0,T)
        write(128,*)oner/beta,r1

c        do ick=1,Ncw  ! do all K in the irreducible wedge  ! ick fixes to icK0, doesn't loop
c         Open some files.
          if(ick.lt.10) then
            fmt1="('chbrbk',i1)"
          else if(ick.lt.100) then
            fmt1="('chbrbk',i2)"
          else
            write(6,*) 'ick too large'
            stop
          end if
            
          write(string1,fmt1)        ick        !open model_chbk$ick
          open(unit=21,file=string1,status='unknown')

          do i=1,nfb
            write(21,*) wb(i),chi2(i,ick),dwb(i)
          end do
          close(21)
c        end do  ! ick fixes to icK0, doesn't loop
        
c        open(unit=21,file='chbrbr0',status='unknown')
c        do i=1,nfb
c        r1=0.0
c        do ick=1,Ncw
c          r1=r1+ickdeg(ick)*chi2(i,ick)
c        end do
c        write(21,*) wb(i),r1/float(Nc),dwb(i)
c        end do  
   
c====================================================================

c       Form the real part, chi1, using Kramers Kronig
c
c                 /    -chi2(x)/pi
c       chi1(w) = | dx -----------
c                 /      w - x
c

c        do ick=1,Ncw   ! ick fixes to icK0, doesn't loop
        do i=1,nfb
          chi1(i,ick)=0.0

c          chi1_B1g(i,ick)=0.0
c          chi1_B1g_1(i,ick)=0.0
c          chi1_B2g(i,ick)=0.0
c          chi1_B2g_1(i,ick)=0.0

          chipp1(i,ick)=0.0
          chipp1_d(i,ick)=0.0
          chipp1_1(i,ick)=0.0
          do j=1,i-1          ! positive x < w
            chi1(i,ick)=chi1(i,ick)+dwb(j)*chi2(j,ick)/(wb(i)-wb(j))
c            chi1_B1g(i,ick)=chi1_B1g(i,ick)+dwb(j)*chi2_B1g(j,ick)/(wb(i)-wb(j))
c            chi1_B1g_1(i,ick)=chi1_B1g_1(i,ick)+dwb(j)*chi2_B1g_1(j,ick)/(wb(i)-wb(j))
c            chi1_B2g(i,ick)=chi1_B2g(i,ick)+dwb(j)*chi2_B2g(j,ick)/(wb(i)-wb(j))
c            chi1_B2g_1(i,ick)=chi1_B2g_1(i,ick)+dwb(j)*chi2_B2g_1(j,ick)/(wb(i)-wb(j))
            chipp1(i,ick)=chipp1(i,ick)+dwb(j)*chipp2(j,ick)/(wb(i)-wb(j))
            chipp1_d(i,ick)=chipp1_d(i,ick)+dwb(j)*chipp2_d(j,ick)/(wb(i)-wb(j))
            chipp1_1(i,ick)=chipp1_1(i,ick)+dwb(j)*chipp2_1(j,ick)/(wb(i)-wb(j))
          end do
          do j=i+1,nfb        ! positive x > w
            chi1(i,ick)=chi1(i,ick)+dwb(j)*chi2(j,ick)/(wb(i)-wb(j))
c            chi1_B1g(i,ick)=chi1_B1g(i,ick)+dwb(j)*chi2_B1g(j,ick)/(wb(i)-wb(j))
c            chi1_B1g_1(i,ick)=chi1_B1g_1(i,ick)+dwb(j)*chi2_B1g_1(j,ick)/(wb(i)-wb(j))
c            chi1_B2g(i,ick)=chi1_B2g(i,ick)+dwb(j)*chi2_B2g(j,ick)/(wb(i)-wb(j))
c            chi1_B2g_1(i,ick)=chi1_B2g_1(i,ick)+dwb(j)*chi2_B2g_1(j,ick)/(wb(i)-wb(j))
            chipp1(i,ick)=chipp1(i,ick)+dwb(j)*chipp2(j,ick)/(wb(i)-wb(j))
            chipp1_d(i,ick)=chipp1_d(i,ick)+dwb(j)*chipp2_d(j,ick)/(wb(i)-wb(j))
            chipp1_1(i,ick)=chipp1_1(i,ick)+dwb(j)*chipp2_1(j,ick)/(wb(i)-wb(j))
          end do
          do j=1,nfb          ! negative x 
            chi1(i,ick)=chi1(i,ick)-dwb(j)*chi2(j,ick)/(wb(i)+wb(j))
c            chi1_B1g(i,ick)=chi1_B1g(i,ick)-dwb(j)*chi2_B1g(j,ick)/(wb(i)+wb(j))
c            chi1_B1g_1(i,ick)=chi1_B1g_1(i,ick)-dwb(j)*chi2_B1g_1(j,ick)/(wb(i)+wb(j))
c            chi1_B2g(i,ick)=chi1_B2g(i,ick)-dwb(j)*chi2_B2g(j,ick)/(wb(i)+wb(j))
c            chi1_B2g_1(i,ick)=chi1_B2g_1(i,ick)-dwb(j)*chi2_B2g_1(j,ick)/(wb(i)+wb(j))
            chipp1(i,ick)=chipp1(i,ick)-dwb(j)*chipp2(j,ick)/(wb(i)+wb(j))
            chipp1_d(i,ick)=chipp1_d(i,ick)-dwb(j)*chipp2_d(j,ick)/(wb(i)+wb(j))
            chipp1_1(i,ick)=chipp1_1(i,ick)-dwb(j)*chipp2_1(j,ick)/(wb(i)+wb(j))
          end do          
        end do
c        end do  ! ick fixes to icK0, doesn't loop
        chi1=-chi1/pi
c        chi1_B1g=-chi1_B1g/pi
c        chi1_B1g_1=-chi1_B1g_1/pi
c        chi1_B2g=-chi1_B2g/pi
c        chi1_B2g_1=-chi1_B2g_1/pi
        chipp1=-chipp1/pi
        chipp1_d=-chipp1_d/pi
        chipp1_1=-chipp1_1/pi

c       Check that chi1(w-->0) is postive

c        do ick=1,Ncw   ! ick fixes to icK0, doesn't loop
          if(chi1(1,ick).lt.0.0) then
            write(6,*) 'sign error in chiph1(1,ick),ick=',ick
            stop
          end if
c          if(chi1_B1g(1,ick).lt.0.0) then
c            write(6,*) 'sign error in chiph1_B1g(1,ick),ick=',ick
c            stop
c          end if
c          if(chi1_B2g(1,ick).lt.0.0) then
c            write(6,*) 'sign error in chiph1_B2g(1,ick),ick=',ick
c            stop
c          end if
 
         if(chipp1(1,ick).lt.0.0) then
            write(6,*) 'sign error in chipp1(1,ick),ick=',ick
            stop
          end if
c        end do   ! ick fixes to icK0, doesn't loop

	


c******************************************************************************
c       Now chiph=chiph1+ii*chiph2 is the ph bare bubble.  To form the PT
c       bubble, we will use
c
c                chiph                    chiph
c       chi_c= ----------        chi_s= ---------   
c               1+U*chiph               1-U*chiph
c
c       in the spin and charge channels.  In the pairing channel, we get
c
c                chipp
c       chi_sc= ---------
c               1+U*chipp
c
c       I think a similar form to that below is required for chi_sc as well.
c       Needs to be checked!!!!!!!!!
c
c       For Raman B1g we get
c
c                             U*chi_B1g_1*chi_1_B1g
c       chi_b1g = chi_B1g  +  ---------------------
c                                   1+U*chiph
c
c       with a similar form for B2g.
c******************************************************************************

        chi_c=2.0*(chi1+ii*chi2)/(1+U*(chi1+ii*chi2))
c        chi_B1g=2.0*((chi1_B1g+ii*chi2_B1g)+
c     &               U*(chi1_B1g_1+ii*chi2_B1g_1)**2/(1+U*(chi1+ii*chi2)))
c        chi_B2g=2.0*((chi1_B2g+ii*chi2_B2g)+
c     &               U*(chi1_B2g_1+ii*chi2_B2g_1)**2/(1+U*(chi1+ii*chi2)))
        chi_p=2.0*((chipp1_d+ii*chipp1_d)+
     &             U*(chipp1_1+ii*chipp2_1)**2/(1+U*(chipp1+ii*chipp2)))
        chi_s=2.0*(chi1+ii*chi2)/(1-U*(chi1+ii*chi2))

        
c******************************************************************************
c       Wrap things up 
c******************************************************************************
        call deallocate_1p
        call deallocate_Spec0
        stop
        end program default
                  
        double precision function fermi(x,beta)
c******************************************************************************
c       A fermi function function
c******************************************************************************
        real(8) :: r1,beta,x
        real(8), parameter :: one=1.0_8, half=0.5_8
        r1=exp(-beta*abs(x))
        fermi=half*((r1+one)+(r1-one)*sign(one,x))/(one+r1)
        return
        end
                
      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2) 
      IMPLICIT real(8) (A - H, O - Z) 
      PARAMETER (NMAX=1000) 
      real(8) X(N),Y(N),Y2(N),U(NMAX) 
      IF (YP1.GT..99E30) THEN 
        Y2(1)=0. 
        U(1)=0. 
      ELSE 
        Y2(1)=-0.5 
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1) 
      ENDIF 
      DO 11 I=2,N-1 
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1)) 
        P=SIG*Y2(I-1)+2. 
        Y2(I)=(SIG-1.)/P 
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1)) 
     *      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P 
11    CONTINUE 
      IF (YPN.GT..99E30) THEN 
        QN=0. 
        UN=0. 
      ELSE 
        QN=0.5 
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1))) 
      ENDIF 
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.) 
      DO 12 K=N-1,1,-1 
        Y2(K)=Y2(K)*Y2(K+1)+U(K) 
12    CONTINUE 
      RETURN 
      END 


      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y) 
      IMPLICIT real(8) (A - H, O - Z) 
      real(8) XA(N),YA(N),Y2A(N) 
      KLO=1 
      KHI=N 
1     IF (KHI-KLO.GT.1) THEN 
        K=(KHI+KLO)/2 
        IF(XA(K).GT.X)THEN 
          KHI=K 
        ELSE 
          KLO=K 
        ENDIF 
      GOTO 1 
      ENDIF 
      H=XA(KHI)-XA(KLO) 
      IF (H.EQ.0.) PAUSE 'Bad XA input.' 
      A=(XA(KHI)-X)/H 
      B=(X-XA(KLO))/H 
      Y=A*YA(KLO)+B*YA(KHI)+ 
     *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6. 
      RETURN 
      END 


