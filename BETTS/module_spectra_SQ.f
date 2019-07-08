      module module_spectra
        implicit none
        save

        integer, parameter :: ispline=1    ! 0 => bilinear,  1 => bicubic, 2 ==> star
        integer, parameter :: icum=0       ! 0 => interpolate Sigma,  1 =>  interpolate Cumulant
        integer, parameter :: iconductivity=1 !  1 for conductivity if ndim=2
        real(8), parameter :: gamma=0.005d0 ! damping of the ImSigma (cannot be larger than -gamma)

        integer :: nf,nfh
        integer , allocatable:: xo2pw(:,:), xo2mw(:,:)        

c       Allocatable arrays.  First, allocate_Spec0
        complex(8) , allocatable:: Sigma(:,:),chi_c(:,:),chi_s(:,:),chi_p(:,:),
     &                             chi_B1g(:,:),chi_B2g(:,:)

        real(8) , allocatable:: w(:), dw(:), wb(:), dwb(:), 
     &                          wh(:),dwh(:),w3i(:),dw3i(:),
     &                          chi2(:,:),chi2_B1g(:,:),chi2_B2g(:,:),
     &                          chi2_B1g_1(:,:),chi2_B2g_1(:,:), 
     &                          chi1(:,:), chi1_B1g(:,:),chi1_B2g(:,:),
     &                          chi1_B1g_1(:,:),chi1_B2g_1(:,:),
     &                          chipp2(:,:), chipp1(:,:),
     &                          chipp2_d(:,:), chipp1_d(:,:),
     &                          chipp2_1(:,:), chipp1_1(:,:), 
     &                          Akw(:,:),scAkw(:,:),tau(:)
        
c       Now, allocate_Spec1
        real(8), allocatable :: g1_sp(:),g2_sp(:),Sig_sp_r(:,:,:),
     &                          Sig_sp_i(:,:,:),Derivs_sp_r(:,:,:),
     &                          Derivs_sp_i(:,:,:)

        end module module_spectra

c*************************************************************************
        subroutine allocate_Spec0
c*************************************************************************
c       Allocate a few arrays that are used in many of the spectral
c       analysis codes.  We need to know nf and Nc
c*************************************************************************
        use module_spectra
        use Global
        implicit none
        integer :: info,infot
c*************************************************************************

        infot=0

        allocate (Sigma(nf,Nc),stat=info)
        infot=infot+info
        allocate (chi_c(nf,Nc),stat=info)
        infot=infot+info
        allocate (chi_s(nf,Nc),stat=info)
        infot=infot+info
        allocate (chi_p(nf,Nc),stat=info)
        infot=infot+info
        allocate (chi_B1g(nf,Nc),stat=info)
        infot=infot+info
        allocate (chi_B2g(nf,Nc),stat=info)
        infot=infot+info
        allocate (chi1(nf,Nc),stat=info)
        infot=infot+info
        allocate (chi1_B1g(nf,Nc),stat=info)
        infot=infot+info
        allocate (chi1_B2g(nf,Nc),stat=info)
        infot=infot+info
        allocate (chi1_B1g_1(nf,Nc),stat=info)
        infot=infot+info
        allocate (chi1_B2g_1(nf,Nc),stat=info)
        infot=infot+info
        allocate (chi2(nf,Nc),stat=info)
        infot=infot+info
        allocate (chi2_B1g(nf,Nc),stat=info)
        infot=infot+info
        allocate (chi2_B2g(nf,Nc),stat=info)
        infot=infot+info
        allocate (chi2_B1g_1(nf,Nc),stat=info)
        infot=infot+info
        allocate (chi2_B2g_1(nf,Nc),stat=info)
        infot=infot+info
        allocate (chipp1(nf,Nc),stat=info)
        infot=infot+info
        allocate (chipp2(nf,Nc),stat=info)
        infot=infot+info
        allocate (chipp1_d(nf,Nc),stat=info)
        infot=infot+info
        allocate (chipp2_d(nf,Nc),stat=info)
        infot=infot+info
        allocate (chipp1_1(nf,Nc),stat=info)
        infot=infot+info
        allocate (chipp2_1(nf,Nc),stat=info)
        infot=infot+info

        allocate (w(nf),stat=info)
        infot=infot+info
        allocate (dw(nf),stat=info)
        infot=infot+info
        allocate (wb(nf),stat=info)
        infot=infot+info
        allocate (dwb(nf),stat=info)
        infot=infot+info
        
        allocate (Akw(nf,Nc),stat=info)
        infot=infot+info
        allocate (scAkw(nf,Nc),stat=info)
        infot=infot+info
        allocate (tau(nf),stat=info)
        infot=infot+info


        if(infot.ne.0) then
          write(6,*) '********** ERROR in allocate_Spec0*********'
          write(6,*) 'infot=',infot
          stop
        end if

        return
        end subroutine allocate_Spec0

c*************************************************************************

c*************************************************************************
        subroutine deallocate_Spec0
c*************************************************************************
c       Deallocate a few arrays that are used in many of the spectral
c       analysis codes.  
c*************************************************************************
        use module_spectra
c*************************************************************************

        deallocate (Sigma)
        deallocate (chi_c)
        deallocate (chi_s)
        deallocate (chi_p)
        deallocate (chi_B1g)
        deallocate (chi_B2g)
        deallocate (chi1)
        deallocate (chi1_B1g)
        deallocate (chi1_B2g)
        deallocate (chi1_B1g_1)
        deallocate (chi1_B2g_1)
        deallocate (chi2)
        deallocate (chi2_B1g)
        deallocate (chi2_B2g)
        deallocate (chi2_B1g_1)
        deallocate (chi2_B2g_1)
        deallocate (chipp1)
        deallocate (chipp2)
        deallocate (chipp1_d)
        deallocate (chipp2_d)
        deallocate (chipp1_1)
        deallocate (chipp2_1)

        deallocate (w)
        deallocate (dw)
        deallocate (wb)
        deallocate (dwb)

        deallocate (Akw)
        deallocate (scAkw)
        deallocate (tau)

        return
        end subroutine deallocate_Spec0
   
   
c*************************************************************************
        subroutine allocate_Spec1
c*************************************************************************
c       Allocate a few arrays that are used in energy_surface.f and 
c       spectra.f.  We need to know nf and N_sp (from tables)
c*************************************************************************
        use module_spectra
        use Global
        implicit none
        integer :: info,infot
c*************************************************************************

        infot=0

        allocate (g1_sp(-N_sp:N_sp),stat=info)
        infot=infot+info
        allocate (g2_sp(-N_sp:N_sp),stat=info)
        infot=infot+info
        allocate (Sig_sp_r(-N_sp:N_sp,-N_sp:N_sp,nfh),stat=info)
        infot=infot+info
        allocate (Derivs_sp_r(-N_sp:N_sp,-N_sp:N_sp,nfh),stat=info)
        infot=infot+info
        allocate (Sig_sp_i(-N_sp:N_sp,-N_sp:N_sp,nfh),stat=info)
        infot=infot+info
        allocate (Derivs_sp_i(-N_sp:N_sp,-N_sp:N_sp,nfh),stat=info)
        infot=infot+info

        allocate (xo2pw(nf,nf),stat=info)
        infot=infot+info
        allocate (xo2mw(nf,nf),stat=info)
        infot=infot+info

        allocate (wh(nfh),stat=info)
        infot=infot+info
        allocate (dwh(nfh),stat=info)
        infot=infot+info
        allocate (w3i(nf+1),stat=info)
        infot=infot+info
        allocate (dw3i(nf+1),stat=info)
        infot=infot+info


        if(infot.ne.0) then
          write(6,*) '********** ERROR in allocate_Spec1*********'
          write(6,*) 'infot=',infot
          stop
        end if

        return
        end subroutine allocate_Spec1

c*************************************************************************



c*************************************************************************
        subroutine deallocate_Spec1
c*************************************************************************
c       Deallocate a few arrays that are used in many of the spectral
c       analysis codes.  
c*************************************************************************
        use module_spectra
        use Global
c*************************************************************************

        deallocate (g1_sp)
        deallocate (g2_sp)
        deallocate (Sig_sp_r)
        deallocate (Derivs_sp_r)
        deallocate (Sig_sp_i)
        deallocate (Derivs_sp_i)

        deallocate (xo2pw)
        deallocate (xo2mw)

        deallocate (wh)
        deallocate (dwh)
        deallocate (w3i)
        deallocate (dw3i)

        return
        end subroutine deallocate_Spec1

       
