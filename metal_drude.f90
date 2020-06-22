!****************************************************************************
!
!   Module: Metal_drude
!
!   PURPOSE: Drude model of metal conductivity
!
!****************************************************************************
module Metal_drude
    use Sim_parameters, only : wp
    ! hbar in eV
    use constants, only : hbar, PI, e
    real(wp), parameter :: m = 9.10938356e-31_wp
    
    ! drude relaxation times tau of metals, 273K
    real(wp), parameter :: tau_cu = 2.7e-14_wp
    real(wp), parameter :: tau_ag = 4.0e-14_wp
    real(wp), parameter :: tau_au = 4.0e-14_wp
    real(wp), parameter :: tau_fe = 0.24e-14_wp
    
    ! free electron densities of metals
    real(wp), parameter :: n_cu = 8.47e28_wp
    real(wp), parameter :: n_ag = 5.86e28_wp
    real(wp), parameter :: n_au = 5.9e28_wp
    real(wp), parameter :: n_fe = 17.0e28_wp    
    
    contains   
    ! drude model
    function ac_sig_drude(freq,sigma_0,tau) result(sigma)
        real(wp) ::  omega
        real(wp), intent(in) :: freq, tau, sigma_0
        complex(wp) :: sigma
        omega = 2 * PI * freq
        sigma = sigma_0  / ((1.0_wp,0.0_wp) + (0.0_wp,-1.0_wp)*omega*tau)
    end function ac_sig_drude
    
    function dc_sig_drude(tau,n) result(sigma_0)
        real(wp), intent(in) :: tau, n
        real(wp) :: sigma_0
        sigma_0 = n*(e**2.0_wp)*tau/m
    end function dc_sig_drude
    
    
    function cu_sig(freq) result(sigma)
        real(wp), intent(in) :: freq
        complex(wp) :: sigma
        sigma = ac_sig_drude(freq,dc_sig_drude(tau_cu,n_cu),tau_cu)
    end function cu_sig
    
    function ag_sig(freq) result(sigma)
        real(wp), intent(in) :: freq
        complex(wp) :: sigma
        sigma = ac_sig_drude(freq,dc_sig_drude(tau_au,n_au),tau_au)
    end function ag_sig
    
end module metal_drude