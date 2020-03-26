!****************************************************************************
!
!   class: Layer_Class
!
!   PURPOSE:  calculate parameters of nth layer
!
!****************************************************************************
module Layer_Class
    use Sim_parameters, only : eta_0, k_0, k_rho
    implicit none    
    
    private
    public :: Layer
    
    type Layer
        complex :: eps_t, mu_t, nu_e, nu_h, Z_e, Z_h, kz_e, kz_h, k_t
        complex, dimension(2,2) :: sigma_n, P_n, Z_n, Y_n  
        real :: d
    end type Layer
    
    interface Layer
        procedure :: initalize_layer
    end interface Layer
    
contains
    
    ! constr function of Layer type
    ! calculate the rest of paramters based on eps, mu, and sigma
    function  initalize_layer(eps_t_in, mu_t_in, sigma_x_in, sigma_y_in, nu_e_in, nu_h_in, d_in) result(new_layer)
        implicit none
        complex, intent(in) :: eps_t_in, mu_t_in, sigma_x_in, sigma_y_in, nu_e_in, nu_h_in
        real, intent(in) :: d_in
        type(Layer) :: new_layer
    
        new_layer%eps_t = eps_t_in
        new_layer%mu_t = eps_t_in
        new_layer%nu_e = nu_e_in
        new_layer%nu_h = nu_h_in
        new_layer%d = d_in
        
        new_layer%sigma_n(1,1) = sigma_x_in
        new_layer%sigma_n(1,2) = (0.0,0.0)
        new_layer%sigma_n(2,1) = (0.0,0.0)
        new_layer%sigma_n(2,2) = sigma_y_in
        
        new_layer%k_t = k_0 * (new_layer%eps_t * new_layer%mu_t) ** 0.5
        new_layer%kz_e = (new_layer%k_t**2-k_rho**2/new_layer%nu_e)**0.5
        new_layer%kz_h = (new_layer%k_t**2-k_rho**2/new_layer%nu_h)**0.5
        new_layer%P_n(1,1) = EXP(-1*(0.0,1.0)*new_layer%kz_e*new_layer%d)
        new_layer%P_n(1,2) = (0.0,0.0)
        new_layer%P_n(2,1) = (0.0,0.0)
        new_layer%P_n(2,2) = EXP(-1*(0.0,1.0)*new_layer%kz_h*new_layer%d)
        
        new_layer%Z_e = eta_0 * new_layer%kz_e / k_0 / new_layer%eps_t
        new_layer%Z_h = eta_0 * k_0 * new_layer%mu_t / new_layer%kz_h
        new_layer%Z_n(1,1) = new_layer%Z_e
        new_layer%Z_n(2,2) = new_layer%Z_h
        new_layer%Z_n(1,2) = (0.0,0.0)
        new_layer%Z_n(2,1) = (0.0,0.0)
        
        new_layer%Y_n(1,1) = new_layer%Z_e**-1
        new_layer%Y_n(2,2) = new_layer%Z_h**-1
        new_layer%Y_n(1,2) = (0.0,0.0)
        new_layer%Y_n(2,1) = (0.0,0.0)
    end function initalize_layer
  
end module Layer_Class