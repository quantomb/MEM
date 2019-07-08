        program sopt_dca
c*******************************************************************************
c       to compile type make sopt_dca (or make all) using the 
c       associated Makefile.  This code was designed to read the
c       magority of its inputs from a qmc sigma file.
c*******************************************************************************
c       This is a simulation of the Hubbard model using second-order
c       conserving pertubation theory within the dynamical cluster 
c       approximation.
c
c                        H=H1+H2
c                              +
c                H1= -t SUM ((C (i,s)C(j,s)) + h.c. )
c                       <ij>s                           ! i,j nn
c                                +
c                     -t' SUM ((C (i,s)C(j,s)) + h.c. )
c                       <<ij>>s                         ! i,j nnn
c
c                H2= U SUM (ndu-1/2)(ndd-1/2) + ed*( ndu + ndd )
c                       i     i        i               i     i
c
c       H is particle-hole symmetric whenever ed=0.  For a two-dimensional
c       system, this means that 
c
c               +             +    i
c       C  --> C exp(iQ.i) = C (-1)
c        i      i             i
c                                                +              +
c       C  = SUM exp(ik.i) C  --> SUM exp(ik.i) C exp(-iQ.i) = C
c        k    i             i                    i              Q-k
c
c       or if ed=0, G (tau) = -G   (-tau) = G   (beta-tau)
c                    k          Q-k          Q-k
c       and
c                   G (iwn) = - G   (-iwn)  real parts opposite sign
c                    k           Q-k        imaginary parts equal sign
c
c*******************************************************************************
        use Global
        implicit none
c*******************************************************************************
        character*72 linemc,string1,string2
        
        integer, parameter :: itermax=2000  ! max number of iterations
        
        integer :: i,j,k,l,n,m,p,ick,icq,ic,icp,iequiv,info,
     &             run,meas,minbound,maxbound,ix,iy,iter,nordpref,
     &             unit,minp,ios

        real(kind) :: test,testl,tolerance,kx,ky,kz,dnded,dndedl, edl,
     &                nfoccl,nfocc,Epsk,r1,r2,damp,damp1,xsc
        
        real(kind), allocatable :: tau(:)

        complex(kind) :: csum        
        complex(kind), allocatable :: GcKf(:,:),GcKfsc(:,:),sigma(:,:),
     &                                sigma1(:,:),chi(:,:),chicc(:)
c*******************************************************************************


c*******************************************************************************
c       Read in some parameters, and open some files.
c*******************************************************************************
        write(6,"('      enter the QMC sigma filename:  ',$)") 
        read(5,"(a72)") linemc
        open(unit=92,file=linemc,status='old')


c       Read some parameters from the QMC sigma file
 401    read(92,"(a72)") linemc
        if(index(linemc(1:72),'Results').ne.0 ) then
          read(92,*) nl,nwn,fill,cluster,ntot,ndim
        else
         goto 401
        endif
 402    read(92,"(a72)") linemc
        if(index(linemc(1:72),'tprime').ne.0 ) then
           if(index(linemc(1:72),'tdprime').ne.0 ) then
              if(index(linemc(1:72),'t3prime').ne.0 ) then
                 read(92,*) ed,U,beta,i,run,meas,tprime,tdprime,t3prime
              else
                 read(92,*) ed,U,beta,i,run,meas,tprime,tdprime
                 t3prime=0.d0
              endif
           else 
              read(92,*) ed,U,beta,i,run,meas,tprime
              tdprime=0.d0
              t3prime=0.d0
           endif
        else
          goto 402
        endif
        temp=1.0/beta   







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

        write(6,"('    QMC ed=',f8.5,'       Enter ed:  ',$)") ed
        read(5,*) ed    
        
c*******************************************************************************
c       Now that we have some parameters, we may allocate the arrays
c       and form the tables.        
c*******************************************************************************
        call allocate_1p
        allocate(tau(nwn),stat=info)
        allocate(GcKf(-nwn:nwn,Nc),stat=info)
        allocate(GcKfsc(-nwn:nwn,Nc),stat=info)
        allocate(sigma(-nwn:nwn,Nc),stat=info)
        allocate(sigma1(-nwn:nwn,Nc),stat=info)
        allocate(chi(-nwn:nwn,Nc),stat=info)
        allocate(chicc(0:nwn),stat=info)
        call tables
        write(6,"('Ncw=',i4)") Ncw
        
c*******************************************************************************
c       Write out the K infomation.
c*******************************************************************************
        open(unit=16,status='unknown',file='Kinfo')
        if(ndim.eq.2) then
          write(16,*) ' ick  Kcx(ick)   Kcy(ick)   degeneracy'
          do ick=1,Ncw
            write(16,"(1x,i3,1x,f9.6,2x,f9.6,4x,i3)") ick,Kc(1,ick),
     &Kc(2,ick),ickdeg(ick)
          end do
          write(16,*) '  '
          do ick=Ncw+1,Nc
            write(16,"(1x,i3,1x,f9.6,2x,f9.6)") ick,Kc(1,ick),Kc(2,ick)
          end do
        else if(ndim.eq.3) then
          write(16,*) ' ick  Kcx(ick)   Kcy(ick)   Kcz(ick)  degeneracy'
          do ick=1,Ncw
            write(16,"(1x,i3,1x,f9.6,2x,f9.6,2x,f9.6,4x,i3)") 
     &         ick,Kc(1,ick),Kc(2,ick),Kc(3,ick),ickdeg(ick)
          end do
          write(16,*) '  '
          do ick=Ncw+1,Nc
            write(16,"(1x,i3,1x,f9.6,2x,f9.6,2x,f9.6)") 
     &         ick,Kc(1,ick),Kc(2,ick),Kc(3,ick)
          end do
        else
          stop 'bad ndim to sopt_dca'
        end if
        close(16)

c*******************************************************************************
c       Read in the rest of the information.
c*******************************************************************************

        write(6,"('                  Enter dn/ded < 0:  ',$)") 
        read(5,*) dnded 
        if(dnded.ge.0.0) stop 'dnded must be <0 '

        write(6,"('       Enter convergence tolerance:  ',$)") 
        read(5,*) tolerance     

        write(6,"('  Enter xsc (0 cons., 1 site-exc.):  ',$)") 
        read(5,*) xsc   

c       calculate the self-energy in a self-consistent loop (999)
c       First initialize the self energy.
        sigma=(0.0_kind,0.0_kind)

c*******************************************************************************
c       SOPT is iterative.  
c*******************************************************************************
        do 999 iter=1,itermax      
c         calculate the coarse-grained green function     
          if(iter.le.5) then    ! use damping to calculate sigma
            damp=0.5_kind       ! fraction of old
          else
            damp=0.2_kind
          end if


c*******************************************************************************
c         Calculate the coarse-grained cluster propagator GcKf and
c         GcKfsc. If xsc=0, then GcKf=GcKfsc, but if xsc=1, GcKfsc
c         is the cluster-excluded propagator.
c*******************************************************************************

          GcKf = (0.0_kind,0.0_kind)
          GcKfsc = (0.0_kind,0.0_kind)
          nfocc=0.0_kind
          do ic = 1,Nc
            do j=1,nwsc
              if(ndim.eq.2) then
                kx=Kc(1,ic) + kt(1,j)
                ky=Kc(2,ic) + kt(2,j)   
            Epsk=-0.5_kind*(cos(kx)+cos(ky)) -
     &       tprime*(cos(kx)*cos(ky)-1.0_kind)
     &             - 0.5*tdprime*(cos(2.d0*kx)+cos(2.d0*ky))
     &             - t3prime*(cos(2.d0*kx)*cos(ky)
     &                         +cos(kx)*cos(2.d0*ky))  
              else if(ndim.eq.3) then
                kx=Kc(1,ic) + kt(1,j)
                ky=Kc(2,ic) + kt(2,j)
                kz=Kc(3,ic) + kt(3,j)
                Epsk=-0.5_kind*(cos(kx)+cos(ky)+cos(ky)) 
     &              - tprime*(cos(kx)*cos(ky)
     &                       +cos(kx)*cos(kz)
     &                       +cos(ky)*cos(kz))
              end if
              
              do n=0,nwn-1
                GcKf(n,ic)=GcKf(n,ic)+
     &             1.0_kind/(ii*wn(n) - Epsk - ed - sigma(n,ic))
              end do
            end do
            do n=0,nwn-1
              GcKf(n,ic)=GcKf(n,ic)/float(nwsc)
              GcKf(-1-n,ic)=conjg(GcKf(n,ic))
              GcKfsc(n,ic)=GcKf(n,ic)/(1.0+xsc*GcKf(n,ic)*sigma(n,ic))
              GcKfsc(-1-n,ic)=conjg(GcKfsc(n,ic))
              nfocc = nfocc + real(1.0_kind/(ii*wn(n))-GcKf(n,ic))
            end do
          end do

          nfocc=0.5_kind - temp*2.0_kind*nfocc/real(Nc)
          damp1=0.1_kind
          dndedl=dnded
          if(abs(nfocc-nfoccl).gt.0.00001_kind.and.iter.gt.3) then
            dnded=min(-0.2_kind,(nfocc-nfoccl)/(ed-edl))
            dnded=max(-4.0_kind,dnded)
            dnded=(1.0-damp1)*dndedl+damp1*dnded
          end if
c         write(6,*) 'nfocc=',nfocc

c*******************************************************************************
c         Symmetrize the GcKfsc and  x <---> +/-x       
c         hence the self energy by      \ /           
c         imposing the point group       X                          
c         symmetries.                   / \                         
c                                    y <---> +/-y 
c*******************************************************************************
          chi=(0.0,0.0)             
          do n=0,nwn-1                                               
          do ic=1,Nc    
            do iequiv=1,ngroup                             
              chi(n,ic)=chi(n,ic) + GcKfsc(n,ickequ(ic,iequiv))
            end do                  
            chi(n,ic)=chi(n,ic)/real(ngroup)
            chi(-n-1,ic)=conjg(chi(n,ic))                   
          end do                                                    
          end do                                                    
                                                                    
          do n=-nwn,nwn                                         
          do ic=1,Nc                                                
            GcKfsc(n,ic)=chi(n,ic)                            
          end do                                                    
          end do  


c*******************************************************************************
c         Now calculate the perturbation theory self-energy sigma(iwn)
c
c
c        ___                 wm,q
c       /   \               ---<---
c      |     |            /         \
c       \___/            |\         /|
c         |              |  --->---  |
c         |              | wm+wp,q+p |      = sigma(wn)
c         |              |           |
c       --<-- wn,k       ------<------  wn,k
c                         wn-wp, k-p
c
c                         2  2 ___                __
c   U(1/2 -gf(tau=0))   -T U   \   Gf(wn-wp,k-p) \   Gf(wm,q) Gf(wm+wp,q+p)
c                      -----   /__               /__
c                      Nc*Nc   wp,p             wm,q
c
c         First we must calculate the bubble in the second-order part.
c         As discussed by Samuel, this is also needed in the FLEX.
c*******************************************************************************

          chi=(0.0_kind,0.0_kind)
          do icp=1,Nc
            do p=-nwn/2,nwn/2
              if(p.ge.0) then
                maxbound=nwn-p-1
                minbound=-nwn+1
              else
                maxbound=nwn-1
                minbound=-nwn-p+1
              end if
              do icq=1,Nc
              do m=minbound,maxbound
                chi(p,icp)=chi(p,icp) - 
     &                       GcKfsc(m,icq)*GcKfsc(m+p,ickplus(icq,icp))
     &                      - 1.0/(wn(m)*wn(m+p))
              end do
              if(p.eq.0) chi(0,icp)=chi(0,icp) + (0.5*beta)**2
              end do
            end do  !p
          end do  !icp

c*******************************************************************************
c         Now that we have chi, complete the calculation of the
c         second order part of sigma.
c*******************************************************************************
          sigma1=(0.0_kind,0.0_kind)

          do ick=1,Nc
          do n=0,nwn-1
            maxbound=nwn-1-abs(n)
            minbound=-nwn+1+abs(n)
            do icp=1,Nc
            do p=minbound,maxbound
              sigma1(n,ick)=sigma1(n,ick)+
     &                      0.5_kind*chi(p,icp)*
     &                      (GcKfsc(n+p,ickplus(ick,icp))+
     &                       GcKfsc(n-p,ickdiff(ick,icp)))
            end do
            end do
c           If damp != 0, then damping is used to improve convergence
c           Here, sigma1 contains the second-order bubble and
c           u*(nfocc-0.5_kind) is the Hartree part.
            sigma(n,ick)=damp*sigma(n,ick)+
     &                  (1.0_kind-damp)*(
     &                                u*(nfocc-0.5_kind)+
     &                  sigma1(n,ick)*(u*temp)**2/real(Nc**2)
     &                               )
            sigma(-1-n,ick)=conjg(sigma(n,ick))
          end do
          end do

          if(iter.gt.2) then
c           test for convergence by inspecting Im[sigma(0,ick) ]
            testl=test
            test=0.0_kind
            do ick=1,Nc
              test = test + abs(sigma(0,ick))
            end do
            r1=abs(testl-test)
            write(6,"(' test= ',e10.4,'  ed=',f8.5,'  nf=',f8.5,
     &                '  dnded=',f10.5)") r1,ed,nfocc,dnded
            if(r1.lt.tolerance) goto 998
          end if
          
c         Now adjust ed to converge to nfocc=fill/2
          if(iter.gt.2) then
            edl=ed
            nfoccl=nfocc
            ed = ed + (0.5_kind*fill-nfocc)/dnded
          end if
 999    continue
 998    write(6,*) 'convergence in ',iter,' iterations'

c*******************************************************************************
c       DONE PT
c       Now do something with the results...
c*******************************************************************************

c       ..like write it out
        open(unit=33,file='Sigma_PT_mats.dat')
        do ic=1,Nc
        do n=0,nwn
          write(33,"(1x,i3,1x,i3,2x,e13.6,1x,e13.6,2x,e13.6,1x,e13.6)") 
     &n,ic,sigma(n,ic),chi(n,ic)
        end do
        end do
c*******************************************************************************
c       ...like a pade analytic continuation (works except when it doesn't)
c*******************************************************************************

c       Get some additional parameters from the call
        if(ndim.eq.2) then
          write(6,*) ' ick  Kcx(ick)   Kcy(ick)'
          do ick=1,Ncw
            write(6,"(1x,i3,1x,f9.6,2x,f9.6)") ick,Kc(1,ick),Kc(2,ick)
          end do
        else if(ndim.eq.3) then
          write(6,*) ' ick  Kcx(ick)   Kcy(ick)   Kcz(ick)'
          do ick=1,Ncw
            write(6,"(1x,i3,1x,f9.6,2x,f9.6,2x,f9.6)") ick,Kc(1,ick),
     &Kc(2,ick),Kc(3,ick)
          end do
        else
          stop 'bad ndim in sopt_dca'
        end if

        unit=19
        nordpref=-1
        minp=min(40,nwn/2)

        if(Nc.gt.1) then
          do ick=1,Ncw  ! do all K in the irreducible wedge
c           Open some files.
            if(ick.lt.10) then
              write(string1,"('model_Akw',i1)") ick
              write(string2,"('Sigma_pt',i1)") ick
            else
              write(string1,"('model_Akw',i2)") ick
              write(string2,"('Sigma_pt',i2)") ick
            end if
            open(unit=59,file=string1,status='unknown') !model_Akw_#
            open(unit=19,file=string2,status='unknown') !Sigma_pt#
c  Write nomalize constant to default model
            write(unit+40,*) " #normalor:   3.14159265358979"

            do i=0,nwn
              chicc(i)=-ii*sigma(i,ick)
            end do

c           if(iter.eq.1) then          ! get Pade order
              call PADE_FC(chicc,minp,0.1_kind,tau,0,0,ick)
              write(6,*) '  '
              write(6,*) ' PADE SELF ENERGY DATA ick=',ick
              write(6,*) ' order    Im[sigma(0)] '
              write(10,*) '  '
              write(10,*) ' PADE SELF ENERGY DATA ick=',ick
              write(10,*) ' order   Im[sigma(0)] '
              do n=1,min0(minp,nwn/2)
                write(6,*) n,tau(n)
                write(10,*) n,tau(n)
              end do
              write(6,*) '  '
              write(10,*) '  '
              write(6,"('ick=',i2,'  Enter the pade order')") ick
              read(5,*) nordpref
c           end if

c           write(6,*) 'Results in units=',unit,unit+40
            call PADE_FC(chicc,nwn/2,0.1_kind,tau,nordpref,unit,ick)
            close(unit) 
            close(unit+40) 
          
          end do
        else
          nordpref=30
        end if
 
        unit=69
        write(6,*) 'DOS in model_dos'
        open(unit=69,file='model_dos',status='unknown') !model_dos
c       Now calculate the DOS
        do i=0,nwn
          csum=(0.0_kind,0.0_kind)
          do ick=1,Nc
            csum=csum+GcKf(i,ick)
          end do
          chicc(i) = -ii*csum/real(Nc)
        end do
        call PADE_FC(chicc,nwn/2,0.1_kind,tau,nordpref,unit,ick)

c*******************************************************************************
c       Wrap things up.
c*******************************************************************************
        call deallocate_1p
        deallocate(tau   )
        deallocate(GcKf  )
        deallocate(GcKfsc)
        deallocate(sigma )
        deallocate(sigma1)
        deallocate(chi   )
        deallocate(chicc )

        stop
        end


        SUBROUTINE PADE_FC(chi,nchi,Escale,tau,nordpref,unit,ick)
c
c       Given the matsubara-frequency susceptibility chi(iwn), which is 
c       assumed to be complex, this code performs a truncated Pade of each order
c       up to nchipp.  For each order it then returns the zero frequency
c       estimate of
c
c                    
c              lim  Im[ sigma(norder)]
c              w->0     
c
c       which is actually evaluated at w = Escale/100. 
c
c       chi      vector with first index of 0 which contains chi(iwn)
c       nchi     the number of frequencies wn to be used
c       Escale   energy scale of the problem
c       tau      vector containing the dc limits of sigma''(w) for each order
c       nordpref is the order at which the ac result is calculated.
c                if this is <=0, then the ac result is not calculated
c       unit     output is written to file fort.unit and fort.(unit+40)
c                if unit=69, PADE_FC assumes that chi is a local
c                green function, otherwise chi is assumed to be 
c                self energy data.
c       Epsk     Epsilon(kx,ky) the bare dispersion of the data (unit!=69).
c*******************************************************************************
        use Global
        implicit none
c*******************************************************************************
        integer, parameter :: nf=1600
        integer :: nchi,n,n1,i,j,nordpref,unit,iKx,iKy,ick

        real(kind), parameter :: wmax=8.0_kind, delta=0.5_kind
        real(kind) :: tau(*),r1,r2,Escale,r3,r4,y,Epsk,kx,ky,kz,
     &                w(nf),dw(nf)
        complex(kind) :: chi(0:*),apade(0:nwn),c1,c2,z,GcKf(nf),
     &                   Sigma(nf)

c*******************************************************************************
c       Calculate the Pade' coefficients 
c*******************************************************************************
c
c                              a(0)
c       chi(wn)=  -------------------------------      =-i Sigma(i wn)
c                   1 + a(1)*(wn-w0)
c                       ---------------------------
c                        1+ a(2)*(wn-w1)
c                           -----------------------
c                           1+ a(3)*(wn-w2)
c                              --------------------
c                               1+ a(4)*(wn-w3)
c                                  ----------------
c                                  1+ a(5)*(wn-w4)
c                                          .
c                                           .
c                                            .
c       
c
c*******************************************************************************
c       a --> apade
c*******************************************************************************
        apade(0)=chi(0)
        do n=1,nchi
           c1=apade(0)/chi(n)
    
           do i=0,n-2
              c2=1/(-1+c1)
              c1=c2*apade(i+1)*(wn(n)-wn(i))
           end do
           
           apade(n)=(-1+c1)/(wn(n)-wn(n-1))
        end do

c*******************************************************************************
c       Now analytically continue [ iwn --> w or wn --> -iw ] and evaluate 
c       the approximant in the dc limit ( w = Escale/100).
c*******************************************************************************
        z=-ii*Escale/100.0_kind    ! essentially the DC limit.

        do n=1,nchi   ! n is the order of the Pade'.
          c1=1.0_kind+apade(n)*(z-wn(n-1))
          if(n.gt.1) then
            do n1=1,n-1
              c1=1.0_kind+apade(n-n1)*(z-wn(n-n1-1))/c1
            end do
          end if
          c1=apade(0)/c1
          tau(n)=dreal(c1)        ! The DC limit, order-by-order n
        end do

        if(unit.ne.69) then
c         we need K(ick) and the dispersion
          kx=Kc(1,ick)
          ky=Kc(2,ick)
           Epsk=-0.5_kind*(cos(kx)+cos(ky)) -
     &       tprime*(cos(kx)*cos(ky)-1.0_kind)
     &             - 0.5*tdprime*(cos(2.d0*kx)+cos(2.d0*ky))
     &             - t3prime*(cos(2.d0*kx)*cos(ky)
     &                         +cos(kx)*cos(2.d0*ky))  

        end if
        
        if (nordpref.gt.0) then
c         for some preferred order, calculate the self energy,
c         and write the result into file unit.

          if(unit.ne.69) write(unit,"('#     w        dw   Re[sigma(w)] 
     &Im[sigma(w)] A(K,w)')") 

          do i=-nf/2+1,nf/2
            n=i+nf/2
            y=(i-0.5_kind)/real(nf)
            w(n)=delta*tan(pi*y)
            dw(n)=pi*(delta**2+(w(n))**2)/(delta*nf)
            z=-ii*w(n)
c           form the approximant.
            c1=1.0_kind+apade(nordpref)*(z-wn(nordpref-1))
            if(nordpref.gt.1) then
              do n1=1,nordpref-1
                c1=1.0_kind+apade(nordpref-n1)*(z-wn(nordpref-n1-1))/c1
              end do
            end if
            c1=apade(0)/c1
            Sigma(n)=ii*c1
            if(abs(w(n)).lt.wmax.and.unit.ne.69) then
              r2=dreal(c1)
              r3=-aimag(c1)
              if(r2.gt.0.0_kind) r2=-0.00001_kind ! causality filter
              c1=dcmplx(r2,-r3)
              r4=-aimag(1.0_kind/(w(n)-Epsk-ii*c1))/pi
              write(unit,"(f10.6,1x,f10.6,1x,f10.6,1x,f10.6,1x,f10.6)") 
     &w(n),dw(n),r3,r2,r4
            end if
          end do

          if(unit.ne.69) then
c           Make a coarse-grained model
            GcKf = (0.0_kind,0.0_kind)
            do i=1,nwsc
              if(ndim.eq.2) then
                kx=Kc(1,ick) + kt(1,i)
                ky=Kc(2,ick) + kt(2,i)
            Epsk=-0.5_kind*(cos(kx)+cos(ky)) -
     &       tprime*(cos(kx)*cos(ky)-1.0_kind)
     &             - 0.5*tdprime*(cos(2.d0*kx)+cos(2.d0*ky))
     &             - t3prime*(cos(2.d0*kx)*cos(ky)
     &                         +cos(kx)*cos(2.d0*ky))  
              else if(ndim.eq.3) then
                kx=Kc(1,ick) + kt(1,i)
                ky=Kc(2,ick) + kt(2,i)
                kz=Kc(3,ick) + kt(2,i)
                Epsk=-0.5_kind*(cos(kx)+cos(ky)+cos(ky)) 
     &              - tprime*(cos(kx)*cos(ky)
     &                       +cos(kx)*cos(kz)
     &                       +cos(ky)*cos(kz))
              end if
              do n=1,nf
                GcKf(n)=GcKf(n)+
     &                  1.0_kind/(w(n) - Epsk - ed - Sigma(n))
              end do
            end do
            do n=1,nf
              GcKf(n)=GcKf(n)/real(nwsc)
              r3=-aimag(GcKf(n))/pi
              if(r3.lt.0.0000001_kind) r3=0.0000001_kind
              if(abs(w(n)).lt.wmax) write(unit+40,
     &"(e12.6,1x,e12.6,1x,e12.6)") w(n),r3,dw(n)
            end do
          else
            write(6,*) 'Write out DOS in unit',unit
            do n=1,nf
              r3=-aimag(Sigma(n))/pi
              if(r3.lt.0.0000001_kind) r3=0.0000001_kind
              if(abs(w(n)).lt.wmax) write(unit,
     &"(e12.6,1x,e12.6,1x,e12.6)") w(n),r3,dw(n)
            end do
          end if
        end if
        
        return
        end
