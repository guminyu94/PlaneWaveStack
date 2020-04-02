!************************************************************************
!
!   Module: Sim_parameters
!
!   PURPOSE: Paraemters for sim space
!
!****************************************************************************
Module Sim_parameters
    implicit none
    public
    
    ! define the precision and constant
    integer, parameter :: wp = KIND(0.0d0)
    real(wp) :: freq, k_0, lambda_0, omega
    real(wp), parameter :: eps_0 = 8.85E-12_wp, mu_0 = 1.257E-6_wp, c_0 = 2.998E8_wp, eta_0 = 376.730313_wp, PI = 3.14159265358979_wp
    complex(wp) :: k_rho
    
    ! define model varaiables 
    complex(wp), allocatable :: eps_t(:), mu_t(:), sigma_x(:), sigma_y(:), nu_e(:), nu_h(:), sigma_xy(:), sigma_yx(:)
    real(wp), allocatable :: d(:)
    real(wp) :: theta
    real(wp) :: xi
    integer :: n_layers
    
    contains
    subroutine update_freq(freq_in)
        real(wp), intent(in) :: freq_in
        freq = freq_in
        omega = 2.0_wp * PI * freq 
        k_0 = omega / c_0
        lambda_0 = 2.0_wp * PI / k_0   
    end subroutine update_freq
    
end module Sim_Parameters