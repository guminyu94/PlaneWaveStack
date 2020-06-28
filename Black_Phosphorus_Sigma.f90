!****************************************************************************
!
!   Module: Black_Phosphorus_Sigma
!
!   PURPOSE: Compute sigma tensor of black phosphorus
!
!****************************************************************************
    Module Black_Phosphorus_Sigma
    use Constants, only : hbar, e, PI
    use Sim_parameters, only : wp
    implicit none
    ! hbar is in eV
    real(wp), parameter :: a_bp = 0.223E-9_wp
    real(wp), parameter :: delta_bp = 2.0_wp * e
    real(wp), parameter :: eta_bp = 10.0E-3_wp * e
    real(wp), parameter :: m_0 = 9.11E-31_wp
    real(wp), parameter :: eta_c_bp = ((hbar*e) ** 2.0_wp) / (0.4_wp * m_0)
    real(wp), parameter :: nu_c_bp = ((hbar*e) ** 2.0_wp) / (1.4_wp * m_0)
    real(wp), parameter :: gamma_bp = 4 * a_bp / PI * e
    real(wp), parameter :: m_arm = ((hbar*e) ** 2.0_wp) / (2.0_wp * (gamma_bp ** 2.0_wp) / delta_bp + eta_c_bp)
    real(wp), parameter :: m_zig = ((hbar*e) ** 2.0_wp) / (2.0_wp * nu_c_bp)
    
    real(wp), parameter :: m_xx = 0.15_wp * m_0
    real(wp), parameter :: m_yy = 0.7 * m_0
    real(wp), parameter :: m_star = (m_xx * m_yy)**0.5_wp
    real(wp), parameter :: xi_bp = eta_bp / (hbar*e)
    
    
    contains 
    
    function bp_sigma(freq,n) result(sigma_mat) 
        real(wp), intent(in) :: freq
        ! carrier density (adjustable)
        real(wp), intent(in) :: n 
        complex(wp) , dimension(2,2) :: sigma_mat
        real(wp) :: omega
        complex(wp) :: sigma_arm 
        complex(wp) :: sigma_zig
        real(wp) :: D_arm, D_zig
        omega = 2 * PI * freq
        D_arm = PI * (e**2.0_wp) * n / m_arm 
        D_zig = PI * (e**2.0_wp) * n / m_zig
        sigma_arm = (0.0_wp,1.0_wp) * D_arm / (PI * (omega + (0.0_wp,1.0_wp) * eta_bp / (hbar*e)))
        sigma_zig = (0.0_wp,1.0_wp) * D_zig / (PI * (omega + (0.0_wp,1.0_wp) * eta_bp / (hbar*e)))
        sigma_mat(1,1) =  sigma_arm
        sigma_mat(1,2) = (0.0_wp,0.0_wp)
        sigma_mat(2,1) = (0.0_wp,0.0_wp)
        sigma_mat(2,2) = sigma_zig
    end function bp_sigma
    
    function bp_sigma_b0(freq,n,b_0) result(sigma_mat)
        real(wp), intent(in) :: freq
        ! carrier density (adjustable)
        real(wp), intent(in) :: n 
        complex(wp) , dimension(2,2) :: sigma_mat
        real(wp), intent(in) :: b_0
        
        real(wp) :: omega
        real(wp) :: omega_c
        omega = 2 * PI * freq
        omega_c = e * b_0 / m_star
        sigma_mat(1,1) = (0.0_wp,1.0_wp) * n * (e**2.0_wp) / m_xx * (omega + (0.0_wp,1.0_wp) * xi_bp) / ( (omega + (0.0_wp,1.0_wp)* xi_bp)**2.0_wp - omega_c**2.0_wp )
        sigma_mat(1,2) = - n * (e**2.0_wp)/m_star * omega_c / ((omega+(0.0_wp,1.0_wp)*xi_bp)**2.0_wp - omega_c**2.0_wp)
        sigma_mat(2,1) = - sigma_mat(1,2)
        sigma_mat(2,2) = (0.0_wp,1.0_wp) * n * (e**2.0_wp) / m_yy * (omega + (0.0_wp,1.0_wp) * xi_bp) / ( (omega + (0.0_wp,1.0_wp)* xi_bp)**2.0_wp - omega_c**2.0_wp )
     end function bp_sigma_b0
    
    !!! plot sigma of bp
    subroutine plot_bp_sigma(freq_start,freq_end,n_p,n,b_0)
        use Plot_Pgplot
        
        real(wp), intent(in) :: freq_start, freq_end
        ! carrier density (adjustable)
        real(wp), intent(in) :: n
        integer, intent(in) :: n_p        
        real(wp), intent(in) :: b_0
        real(wp) :: freq_cur, freq_step
        real, allocatable :: freq_array(:)
        real, allocatable, dimension(:,:) :: sigma_x,sigma_y, sigma_h
        complex(wp), dimension(2,2) :: sigma_bp_cur
        integer :: i
        allocate(freq_array(n_p))
        allocate(sigma_x(2,n_p))
        allocate(sigma_y(2,n_p))
        allocate(sigma_h(2,n_p))

        if (n_p .EQ. 1) then
            freq_step = 0.0_wp
        else
            freq_step = (freq_end-freq_start)/(n_p-1)
        end if
        
        do i = 1,n_p
            freq_cur = freq_start + freq_step * i
            sigma_bp_cur = bp_sigma_b0(freq_cur,n,b_0)
            freq_array(i) = real(freq_cur / 1E12_wp)
            sigma_x(1,i) = real(sigma_bp_cur(1,1))
            sigma_x(2,i) = real(sigma_bp_cur(1,1))
            sigma_y(1,i) = AIMAG(sigma_bp_cur(1,1))
            sigma_y(2,i) = AIMAG(sigma_bp_cur(1,1))
            sigma_h(1,i) = real(sigma_bp_cur(2,2))
            sigma_h(2,i) = real(sigma_bp_cur(2,2))
        end do
        
        call plot_1d(freq_array,sigma_h, x_label = 'Freq(THz)', y_label = 'Z(ohm)', title = 'Sigma Plot')
        
    end subroutine plot_bp_sigma
end module Black_Phosphorus_Sigma