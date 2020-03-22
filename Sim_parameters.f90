Module Sim_parameters
    public
    real :: freq, k_0, lambda_0, omega
    real, parameter :: eps_0 = 8.85E-12, mu_0 = 1.257E-6, c_0 = 2.998E8, PI = 3.1415926E3, eta_0 = 376.730313
    complex :: k_rho
    
    contains
    subroutine update_freq(freq_in)
    freq = freq_in
    omega = 2 * pi * freq 
    k_0 = omega / c_0
    lambda_0 = 2 * PI / k_0   
    end subroutine update_freq
    
end module Sim_Parameters