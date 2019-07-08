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
     &             ik, isk, nkpoints, info, infot, ikx, iky, deg, nfb

        real(kind) :: r1, r2, r3, r4, r5, kx, ky, kx1, ky1, kx2, 
     &                ky2, offset, res, hall, vx, vyy, fermi, delta,
     &                y,rsumg
        
        complex(kind) :: c1, c2
        
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
c******************************************************************************
c       Open the MEM dos files, readin the spectra, and then spline it.
c******************************************************************************
        do ick=1,Ncw  ! do all K in the irreducible wedge
c         Open some files.
          if(ick.lt.10) then
            if(filenum.lt.10) then
              fmt1="('dos',i1,'_',i1,'QMC')"
            else if(filenum.lt.100) then
              fmt1="('dos',i2,'_',i1,'QMC')"
            else if(filenum.lt.1000) then
              fmt1="('dos',i3,'_',i1,'QMC')"
            else if(filenum.lt.10000) then
              fmt1="('dos',i4,'_',i1,'QMC')"
            end if
          else if(ick.lt.100) then
            if(filenum.lt.10) then
              fmt1="('dos',i1,'_',i2,'QMC')"
            else if(filenum.lt.100) then
              fmt1="('dos',i2,'_',i2,'QMC')"
            else if(filenum.lt.1000) then
              fmt1="('dos',i3,'_',i2,'QMC')"
            else if(filenum.lt.10000) then
              fmt1="('dos',i4,'_',i2,'QMC')"
            end if
          else
            write(6,*) 'ick too large'
            stop
          end if
            
          

          write(string1,fmt1) filenum,ick        !open dos_ickQMC
          write(6,*)string1

          open(unit=21,file=string1,status='old')
          
c          read(21,1000) linemc    
          do i=1,1000
            read(21,*,IOSTAT=istat) w(i),r1,dw(i)

            if(istat.ne.0) goto 998
            Akw(i,ick)=max(1.0e-4,r1)
          end do
 998      continue

          if(ick.eq.1) then

            nf=i-1
           else if(nf.ne.i-1) then       ! trap for errors
 

            write(6,*) 'dos files are not compatible'
            stop
          end if
          


c          call spline(w,Akw(1,ick),nf,0.0_8,0.0_8,scAkw(1,ick))
          call spline(w,Akw(:,ick),nf,0.0_8,0.0_8,scAkw(:,ick))

c         To evaluate at w(l)
c         call splint(w,Akw(1,ick),scAkw(1,ick),nf,w(l),A_spline(l))    

          
        end do 
        
        write(6,*)'after read in spectra and splint it'
 
c******************************************************************************
c       Now generate an inhomogeneous frequency grid, and calculate the
c       polarization bubble on this grid (both real and imag. parts).
c******************************************************************************

        delta=0.7 
        nfb=1000     
        write(6,*)'w(nf)=',w(nf)
        do i=1,nfb
          y=(i-0.5)/float(2*nfb)  !positive freqs. only.
          wb(i)=delta*tan(pi*y)
          if(wb(i).gt.w(nf)) exit
          dwb(i)=pi*(delta**2+wb(i)**2)/(delta*2*nfb)   
        end do 


        write(6,*)'after loop i=1,nfb , i=',i
        nfb=i-1



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
c
c
c In the pp channel: chipp2(Q,w)=0.5*(chi^l+chi^r)
c
c             chi^l                +     chi^r
c          ----<---- -K  -w'+w/2      ---->---- -K  -w'-w/2
c                                  +
c          ----<---- K+Q  w'+w/2      ---->---- K+Q  w'-w/2
c
c
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
c
        chi2=0.0
        chi2_B1g=0.0
        chi2_B1g_1=0.0
        chi2_B2g=0.0
        chi2_B2g_1=0.0
        chipp2=0.0
        chipp2_d=0.0
        chipp2_1=0.0
        do ick=1,Ncw                   ! Q
        do jck=1,Nc                    ! K sum over full zone
          isk=ickmap(ickplus(ick,jck)) ! K+Q
          ik=ickmap(jck)               ! K mapped into FBZ
        do i=1,nfb                     ! w
        do j=1,nf                      ! w' sum over frequency
          call splint(w,Akw(1,isk),scAkw(1,isk),nf,w(j)+0.5*wb(i),r1) !A(K+Q,w'+w/2)
          call splint(w,Akw(1,ik),scAkw(1,ik),nf,w(j)-0.5*wb(i),r2)   !A(K,w'-w/2)
          call splint(w,Akw(1,ik),scAkw(1,ik),nf,-w(j)+0.5*wb(i),r3)  !A(K,-w'+w/2)

c    for calculation of Im\chi^{0r}(Q,w)
          call splint(w,Akw(1,isk),scAkw(1,isk),nf,w(j)-0.5*wb(i),r4)   !A(K+Q,w'-w/2)
          call splint(w,Akw(1,ik),scAkw(1,ik),nf,-w(j)-0.5*wb(i),r5)  !A(K,-w'-w/2)

          chi2(i,ick)=chi2(i,ick)+dw(j)*r1*r2*
     &         (fermi(w(j)-0.5*wb(i),beta)-fermi(w(j)+0.5*wb(i),beta))
          chi2_B1g(i,ick)=chi2_B1g(i,ick)+dw(j)*gB1g(jck,ick)*gB1g(jck,ick)*r1*r2*
     &         (fermi(w(j)-0.5*wb(i),beta)-fermi(w(j)+0.5*wb(i),beta))
          chi2_B1g_1(i,ick)=chi2_B1g_1(i,ick)+dw(j)*gB1g(jck,ick)*r1*r2*
     &         (fermi(w(j)-0.5*wb(i),beta)-fermi(w(j)+0.5*wb(i),beta))
          chi2_B2g(i,ick)=chi2_B2g(i,ick)+dw(j)*gB2g(jck,ick)*gB2g(jck,ick)*r1*r2*
     &         (fermi(w(j)-0.5*wb(i),beta)-fermi(w(j)+0.5*wb(i),beta))
          chi2_B2g_1(i,ick)=chi2_B2g_1(i,ick)+dw(j)*gB2g(jck,ick)*r1*r2*
     &         (fermi(w(j)-0.5*wb(i),beta)-fermi(w(j)+0.5*wb(i),beta))

c    need to consider as: chipp=0.5*(chipp^l+ chipp^r)
          chipp2(i,ick)=chipp2(i,ick)+dw(j)*(       
     &     -r1*r3*(fermi(-w(j)+0.5*wb(i),beta)+fermi(w(j)+0.5*wb(i),beta)-1.0)   !chi^l
     &     +r4*r5*(fermi(-w(j)-0.5*wb(i),beta)+fermi(w(j)-0.5*wb(i),beta)-1.0) ) !chi^r

          chipp2_1(i,ick)=chipp2_1(i,ick)+dw(j)*gdwave(jck)*(
     &     -r1*r3*(fermi(-w(j)+0.5*wb(i),beta)+fermi(w(j)+0.5*wb(i),beta)-1.0)   !chi^l
     &     +r4*r5*(fermi(-w(j)-0.5*wb(i),beta)+fermi(w(j)-0.5*wb(i),beta)-1.0) ) !chi^r

          chipp2_d(i,ick)=chipp2_d(i,ick)+dw(j)*gdwave(jck)*gdwave(jck)*(
     &     -r1*r3*(fermi(-w(j)+0.5*wb(i),beta)+fermi(w(j)+0.5*wb(i),beta)-1.0)   !chi^l
     &     +r4*r5*(fermi(-w(j)-0.5*wb(i),beta)+fermi(w(j)-0.5*wb(i),beta)-1.0) ) !chi^r

        end do
        end do
        end do
        end do

        chi2=pi*chi2/float(Nc)       ! normalize
        chi2_B1g=pi*chi2_B1g/float(Nc)       ! normalize
        chi2_B1g_1=pi*chi2_B1g_1/float(Nc)       ! normalize
        chi2_B2g=pi*chi2_B2g/float(Nc)       ! normalize
        chi2_B2g_1=pi*chi2_B2g_1/float(Nc)       ! normalize
        chipp2=pi*chipp2/float(Nc)       ! normalize       
        chipp2_1=pi*chipp2_1/float(Nc)       ! normalize
        chipp2_d=pi*chipp2_d/float(Nc)       ! normalize

c       Check that chi2/w is postive (i.e. like the conductivity)
        do ick=1,Ncw
          if(chi2(1,ick).lt.0.0.or.chi2(nfb/2,ick).lt.0.0) then
            write(6,*) 'sign error in chiph2(i,ick),ick=',ick
            stop
          end if
          if(chi2_B1g(1,ick).lt.0.0.or.chi2_B1g(nfb/2,ick).lt.0.0) then
            write(6,*) 'sign error in chiph2_B1g(i,ick),ick=',ick
            stop
          end if
          if(chi2_B2g(1,ick).lt.0.0.or.chi2_B2g(nfb/2,ick).lt.0.0) then
            write(6,*) 'sign error in chiph2_B2g(i,ick),ick=',ick
            stop
          end if
          if(chipp2(1,ick).lt.0.0.or.chipp2(nfb/2,ick).lt.0.0) then
            write(6,*) 'sign error in chipp2(i,ick),ick=',ick
            stop
          end if
        end do
   
c====================================================================

c       Form the real part, chi1, using Kramers Kronig
c
c                 /    -chi2(x)/pi
c       chi1(w) = | dx -----------
c                 /      w - x
c
        do ick=1,Ncw
        do i=1,nfb
          chi1(i,ick)=0.0
          chi1_B1g(i,ick)=0.0
          chi1_B1g_1(i,ick)=0.0
          chi1_B2g(i,ick)=0.0
          chi1_B2g_1(i,ick)=0.0
          chipp1(i,ick)=0.0
          chipp1_d(i,ick)=0.0
          chipp1_1(i,ick)=0.0
          do j=1,i-1          ! positive x < w
            chi1(i,ick)=chi1(i,ick)+dwb(j)*chi2(j,ick)/(wb(i)-wb(j))
            chi1_B1g(i,ick)=chi1_B1g(i,ick)+dwb(j)*chi2_B1g(j,ick)/(wb(i)-wb(j))
            chi1_B1g_1(i,ick)=chi1_B1g_1(i,ick)+dwb(j)*chi2_B1g_1(j,ick)/(wb(i)-wb(j))
            chi1_B2g(i,ick)=chi1_B2g(i,ick)+dwb(j)*chi2_B2g(j,ick)/(wb(i)-wb(j))
            chi1_B2g_1(i,ick)=chi1_B2g_1(i,ick)+dwb(j)*chi2_B2g_1(j,ick)/(wb(i)-wb(j))
            chipp1(i,ick)=chipp1(i,ick)+dwb(j)*chipp2(j,ick)/(wb(i)-wb(j))
            chipp1_d(i,ick)=chipp1_d(i,ick)+dwb(j)*chipp2_d(j,ick)/(wb(i)-wb(j))
            chipp1_1(i,ick)=chipp1_1(i,ick)+dwb(j)*chipp2_1(j,ick)/(wb(i)-wb(j))
          end do
          do j=i+1,nfb        ! positive x > w
            chi1(i,ick)=chi1(i,ick)+dwb(j)*chi2(j,ick)/(wb(i)-wb(j))
            chi1_B1g(i,ick)=chi1_B1g(i,ick)+dwb(j)*chi2_B1g(j,ick)/(wb(i)-wb(j))
            chi1_B1g_1(i,ick)=chi1_B1g_1(i,ick)+dwb(j)*chi2_B1g_1(j,ick)/(wb(i)-wb(j))
            chi1_B2g(i,ick)=chi1_B2g(i,ick)+dwb(j)*chi2_B2g(j,ick)/(wb(i)-wb(j))
            chi1_B2g_1(i,ick)=chi1_B2g_1(i,ick)+dwb(j)*chi2_B2g_1(j,ick)/(wb(i)-wb(j))
            chipp1(i,ick)=chipp1(i,ick)+dwb(j)*chipp2(j,ick)/(wb(i)-wb(j))
            chipp1_d(i,ick)=chipp1_d(i,ick)+dwb(j)*chipp2_d(j,ick)/(wb(i)-wb(j))
            chipp1_1(i,ick)=chipp1_1(i,ick)+dwb(j)*chipp2_1(j,ick)/(wb(i)-wb(j))
          end do
          do j=1,nfb          ! negative x 
            chi1(i,ick)=chi1(i,ick)-dwb(j)*chi2(j,ick)/(wb(i)+wb(j))
            chi1_B1g(i,ick)=chi1_B1g(i,ick)-dwb(j)*chi2_B1g(j,ick)/(wb(i)+wb(j))
            chi1_B1g_1(i,ick)=chi1_B1g_1(i,ick)-dwb(j)*chi2_B1g_1(j,ick)/(wb(i)+wb(j))
            chi1_B2g(i,ick)=chi1_B2g(i,ick)-dwb(j)*chi2_B2g(j,ick)/(wb(i)+wb(j))
            chi1_B2g_1(i,ick)=chi1_B2g_1(i,ick)-dwb(j)*chi2_B2g_1(j,ick)/(wb(i)+wb(j))
            chipp1(i,ick)=chipp1(i,ick)-dwb(j)*chipp2(j,ick)/(wb(i)+wb(j))
            chipp1_d(i,ick)=chipp1_d(i,ick)-dwb(j)*chipp2_d(j,ick)/(wb(i)+wb(j))
            chipp1_1(i,ick)=chipp1_1(i,ick)-dwb(j)*chipp2_1(j,ick)/(wb(i)+wb(j))
          end do          
        end do
        end do
        chi1=-chi1/pi
        chi1_B1g=-chi1_B1g/pi
        chi1_B1g_1=-chi1_B1g_1/pi
        chi1_B2g=-chi1_B2g/pi
        chi1_B2g_1=-chi1_B2g_1/pi
        chipp1=-chipp1/pi
        chipp1_d=-chipp1_d/pi
        chipp1_1=-chipp1_1/pi

c       Check that chi1(w-->0) is postive
        do ick=1,Ncw
          if(chi1(1,ick).lt.0.0) then
            write(6,*) 'sign error in chiph1(1,ick),ick=',ick
            stop
          end if
          if(chi1_B1g(1,ick).lt.0.0) then
            write(6,*) 'sign error in chiph1_B1g(1,ick),ick=',ick
            stop
          end if
          if(chi1_B2g(1,ick).lt.0.0) then
            write(6,*) 'sign error in chiph1_B2g(1,ick),ick=',ick
            stop
          end if
          if(chipp1(1,ick).lt.0.0) then
            write(6,*) 'sign error in chipp1(1,ick),ick=',ick
            stop
          end if
        end do

	


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
        chi_B1g=2.0*((chi1_B1g+ii*chi2_B1g)+
     &               U*(chi1_B1g_1+ii*chi2_B1g_1)**2/(1+U*(chi1+ii*chi2)))
        chi_B2g=2.0*((chi1_B2g+ii*chi2_B2g)+
     &               U*(chi1_B2g_1+ii*chi2_B2g_1)**2/(1+U*(chi1+ii*chi2)))
        chi_p=2.0*((chipp1_d+ii*chipp1_d)+
     &             U*(chipp1_1+ii*chipp2_1)**2/(1+U*(chipp1+ii*chipp2)))
        chi_s=2.0*(chi1+ii*chi2)/(1-U*(chi1+ii*chi2))

c******************************************************************************
c       Write out model for the charge susceptibility
c******************************************************************************

c       Form the normalized default model for -cha2(w)/w
        chi2=aimag(chi_c)
        do ick=1,Ncw
        r1=0.0  
        do i=1,nfb
          chi2(i,ick)=chi2(i,ick)/(pi*wb(i))
          r1=r1+dwb(i)*chi2(i,ick)
        end do
        chi2(:,ick)=0.5*chi2(:,ick)/r1
	write(111,*)'norm for cha, ick',r1,ick
        end do


c       Now that we have the 2-particle charge spectra, open the appropriate files
c       and write it out. 

        do ick=1,Ncw  ! do all K in the irreducible wedge
c         Open some files.
          if(ick.lt.10) then
            fmt1="('model_chak',i1)"
          else if(ick.lt.100) then
            fmt1="('model_chak',i2)"
          else
            write(6,*) 'ick too large'
            stop
          end if
            
          write(string1,fmt1)        ick        !open modelchakick
          open(unit=21,file=string1,status='unknown')
          do i=1,nfb
            write(21,*) wb(i),chi2(i,ick),dwb(i)
          end do
          close(21)
        end do
        
        open(unit=21,file='model_cha',status='unknown')
        do i=1,nfb
        r1=0.0
        do ick=1,Ncw
          r1=r1+ickdeg(ick)*chi2(i,ick)
        end do
        write(21,*) wb(i),r1/float(Nc),dwb(i)
        end do

c******************************************************************************
c       Write out model for the B1g Raman susceptibility
c******************************************************************************

c       Form the normalized default model for -chi2_B1g(w)/w
c        chi2_B1g=aimag(chi_B1g)
c        do ick=1,Ncw
c        r1=0.0  
c        do i=1,nfb
c          chi2_B1g(i,ick)=chi2_B1g(i,ick)/(pi*wb(i))
c          r1=r1+dwb(i)*chi2_B1g(i,ick)
c        end do
c        chi2_B1g(:,ick)=0.5*chi2_B1g(:,ick)/r1
c	write(111,*)'norm for chi_B1g, ick',r1,ick
c        end do


c       Now that we have the 2-particle charge spectra, open the appropriate files
c       and write it out. 

        do ick=1,Ncw  ! do all K in the irreducible wedge
c         Open some files.
          if(ick.lt.10) then
            fmt1="('model_chi_B1gk',i1)"
          else if(ick.lt.100) then
            fmt1="('model_chi_B1gk',i2)"
          else
            write(6,*) 'ick too large'
            stop
          end if
            
          write(string1,fmt1)        ick        !open modelchakick
          open(unit=21,file=string1,status='unknown')
          do i=1,nfb
            write(21,*) wb(i),real(chi_B1g(i,ick)),aimag(chi_B1g(i,ick))
          end do
          close(21)
        end do
        
        open(unit=21,file='model_chi_B1g',status='unknown')
        do i=1,nfb
        r1=0.0
        r2=0.0
        do ick=1,Ncw
          r1=r1+ickdeg(ick)*real(chi_B1g(i,ick))
          r2=r2+ickdeg(ick)*aimag(chi_B1g(i,ick))
        end do
        write(21,*) wb(i),r1/float(Nc),r2/float(Nc)
        end do

c******************************************************************************
c       Write out model for the B2g Raman susceptibility
c******************************************************************************

c       Form the normalized default model for -chi2_B2g(w)/w
c        chi2_B2g=aimag(chi_B2g)
c        do ick=1,Ncw
c        r1=0.0  
c        do i=1,nfb
c          chi2_B2g(i,ick)=chi2_B2g(i,ick)/(pi*wb(i))
c          r1=r1+dwb(i)*chi2_B2g(i,ick)
c        end do
c        chi2_B2g(:,ick)=0.5*chi2_B2g(:,ick)/r1
c	write(111,*)'norm for chi_B2g, ick',r1,ick
c        end do


c       Now that we have the 2-particle charge spectra, open the appropriate files
c       and write it out. 

        do ick=1,Ncw  ! do all K in the irreducible wedge
c         Open some files.
          if(ick.lt.10) then
            fmt1="('model_chi_B2gk',i1)"
          else if(ick.lt.100) then
            fmt1="('model_chi_B2gk',i2)"
          else
            write(6,*) 'ick too large'
            stop
          end if
            
          write(string1,fmt1)        ick        !open modelchakick
          open(unit=21,file=string1,status='unknown')
          do i=1,nfb
            write(21,*) wb(i),real(chi_B2g(i,ick)),aimag(chi_B2g(i,ick))
          end do
          close(21)
        end do
        
        open(unit=21,file='model_chi_B2g',status='unknown')
        do i=1,nfb
        r1=0.0
        r2=0.0
        do ick=1,Ncw
          r1=r1+ickdeg(ick)*real(chi_B2g(i,ick))
          r2=r2+ickdeg(ick)*aimag(chi_B2g(i,ick))
        end do
        write(21,*) wb(i),r1/float(Nc),r2/float(Nc)
        end do

c******************************************************************************
c       Write out model for the spin susceptibility
c******************************************************************************

c       Form the normalized default model for -chi2(w)/w
        chi2=aimag(chi_s)
        do ick=1,Ncw
        r1=0.0  
        do i=1,nfb
          chi2(i,ick)=chi2(i,ick)/(pi*wb(i))
          r1=r1+dwb(i)*chi2(i,ick)
        end do
        chi2(:,ick)=0.5*chi2(:,ick)/r1
	write(111,*)'norm for chi, ick',r1,ick
        end do

        do ick=1,Ncw  ! do all K in the irreducible wedge
c         Open some files.
          if(ick.lt.10) then
            fmt1="('model_chik',i1)"
          else if(ick.lt.100) then
            fmt1="('model_chik',i2)"
          else
            write(6,*) 'ick too large'
            stop
          end if
            
          write(string1,fmt1)        ick        !open modelchakick
          open(unit=21,file=string1,status='unknown')
          do i=1,nfb
            write(21,*) wb(i),chi2(i,ick),dwb(i)
          end do
          close(21)
        end do
        
        open(unit=21,file='model_chi',status='unknown')
        do i=1,nfb
        r1=0.0
        do ick=1,Ncw
          r1=r1+ickdeg(ick)*chi2(i,ick)
        end do
        write(21,*) wb(i),r1/float(Nc),dwb(i)
        end do
        
c******************************************************************************
c       Write out model for the pair susceptibility
c******************************************************************************
c       Form the normalized default model for -chi2(w)/w
        chi2=aimag(chi_p)
        do ick=1,Ncw
        r1=0.0  
        do i=1,nfb
          chi2(i,ick)=chi2(i,ick)/(pi*wb(i))
          r1=r1+dwb(i)*chi2(i,ick)
        end do
        chi2(:,ick)=0.5*chi2(:,ick)/r1
        end do

        do ick=1,Ncw  ! do all K in the irreducible wedge
c         Open some files.
          if(ick.lt.10) then
            fmt1="('model_chpk',i1)"
          else if(ick.lt.100) then
            fmt1="('model_chpk',i2)"
          else
            write(6,*) 'ick too large'
            stop
          end if
            
          write(string1,fmt1)        ick        !open modelchakick
          open(unit=21,file=string1,status='unknown')
          write(21,*)'ncw=',Ncw, 'ick=',ick,'nfb=',nfb
          write(21,*)'wb(i)',wb
          write(21,*)'dwb(i)',dwb
          do i=1,nfb
            write(21,*) wb(i),chi2(i,ick),dwb(i)
          end do
          close(21)
        end do
        
        open(unit=21,file='model_chp',status='unknown')
        do i=1,nfb
        r1=0.0
        do ick=1,Ncw
          r1=r1+ickdeg(ick)*chi2(i,ick)
        end do
        write(21,*) wb(i),r1/float(Nc),dwb(i)
        end do


c******************************************************************************
c   write out dressed-single-particle-green-function-bubble d-wave projected pairing suscetptibility:
c******************************************************************************
c***************************************************
c       calculate the normalize factor for the form factor gdwave as: c_g=\Sig_K g(K)^2=sum_ick^Nc gdwave(ick)^2
c***************************************************
        rsumg=0.0
        do ick=1,Nc
          rsumg=rsumg+gdwave(ick)*gdwave(ick)
        end do

c       Form the normalized default model for -chi2(w)/w
        chi2=chipp2_d
        do ick=1,Ncw
        r1=0.0  
        do i=1,nfb
          chi2(i,ick)=chi2(i,ick)/(rsumg*pi*wb(i))
          r1=r1+dwb(i)*chi2(i,ick)
        end do
c        chi2(:,ick)=0.5*chi2(:,ick)/r1
        chi2(:,ick)=0.5*chi2(:,ick)

        end do
        write(117,*)oner/beta,r1
        do ick=1,Ncw  ! do all K in the irreducible wedge

c         Open some files.
          if(ick.lt.10) then
            fmt1="('model_chdk',i1)"
          else if(ick.lt.100) then
            fmt1="('model_chdk',i2)"
          else
            write(6,*) 'ick too large'
            stop
          end if            
          write(string1,fmt1)        ick        !open model_chdk$ick
          open(unit=21,file=string1,status='unknown')
          do i=1,nfb
            write(21,*) wb(i),chi2(i,ick),dwb(i)
          end do
          close(21)
        end do
        
        open(unit=21,file='model_chd',status='unknown')
        do i=1,nfb
        r1=0.0
        do ick=1,Ncw
          r1=r1+ickdeg(ick)*chi2(i,ick)
        end do
        write(21,*) wb(i),r1/float(Nc),dwb(i)
        end do    
c====================================================================

c******************************************************************************
c   write out dressed-single-particle-green-function-bubble pairing suscetptibility (not any form factors applied):
c******************************************************************************
c       Form the normalized default model for -chi2(w)/w
        chi2=chipp2
        do ick=1,Ncw
        r1=0.0  
        do i=1,nfb
          chi2(i,ick)=chi2(i,ick)/(pi*wb(i))
          r1=r1+dwb(i)*chi2(i,ick)
        end do
c        chi2(:,ick)=0.5*chi2(:,ick)/r1
        chi2(:,ick)=0.5*chi2(:,ick)
        end do
        write(128,*)oner/beta,r1
        do ick=1,Ncw  ! do all K in the irreducible wedge
c         Open some files.
          if(ick.lt.10) then
            fmt1="('model_chbk',i1)"
          else if(ick.lt.100) then
            fmt1="('model_chbk',i2)"
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
        end do
        
        open(unit=21,file='model_chb',status='unknown')
        do i=1,nfb
        r1=0.0
        do ick=1,Ncw
          r1=r1+ickdeg(ick)*chi2(i,ick)
        end do
        write(21,*) wb(i),r1/float(Nc),dwb(i)
        end do 

        
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


