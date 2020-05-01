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
end module Black_Phosphorus_Sigma