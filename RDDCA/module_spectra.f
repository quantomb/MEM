        module module_spectra
        implicit none
        save

        integer, parameter :: ispline=1     ! 0 <=> bilinear 1 <=> bicubic
        integer :: nf

c       Allocatable arrays.  First, allocate_Spec0
        complex(8) , allocatable:: Sigma(:,:),chi_c(:,:),chi_s(:,:)

        real(8) , allocatable:: w(:), dw(:), wb(:), dwb(:), chi2(:,:), 
     &                          chi1(:,:), Akw(:,:),scAkw(:,:)
        
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
        allocate (chi1(nf,Nc),stat=info)
        infot=infot+info
        allocate (chi2(nf,Nc),stat=info)
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
        deallocate (chi1)
        deallocate (chi2)
        deallocate (w)
        deallocate (dw)
        deallocate (wb)
        deallocate (dwb)
        deallocate (Akw)
        deallocate (scAkw)

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
        allocate (Sig_sp_r(-N_sp:N_sp,-N_sp:N_sp,nf),stat=info)
        infot=infot+info
        allocate (Derivs_sp_r(-N_sp:N_sp,-N_sp:N_sp,nf),stat=info)
        infot=infot+info
        allocate (Sig_sp_i(-N_sp:N_sp,-N_sp:N_sp,nf),stat=info)
        infot=infot+info
        allocate (Derivs_sp_i(-N_sp:N_sp,-N_sp:N_sp,nf),stat=info)
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

        return
        end subroutine deallocate_Spec1

       
