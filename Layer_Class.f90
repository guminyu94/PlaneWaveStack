module Layer_Class
    use Sim_parameters, only : eta_0, k_0
    implicit none
    private
    public :: initalize_layer, Layer
    
    type Layer
        complex :: eps_t, mu_t, Z_e, Z_h, kz_e, kz_h, k_t, nu_e, nu_h
        complex, dimmension(2,2) :: sigma, P_n  
        real :: d
    end type Layer
    
    contains
    
    ! constr function of Layer type
    ! calculate the rest of paramters based on eps, mu, and sigma
    subroutine  initalize_layer(layer,eps_t_in, mu_t_in, sigma_x_in, sigma_y_in, nu_e_in, nu_h_in, d_in)
    implicit none
    type(Layer), intent(inout) :: layer
    complex, intent(in) :: eps_t_in, mu_t_in, sigma_x_in, sigma_y_in, nu_e_in, nu_h_in
    real, intent(in) :: d_in
    
    layer%eps_t = eps_t_in
    layer%mu_t = eps_t_in
    layer%nu_e = nu_e_in
    layer%nu_h = nu_h_in
    
    layer%sigma = (/ sigma_x_in, 0, 0, sigma_y_in/)
    layer%k_t = k_0 * (layer%eps_t * layer%mu_t) ** 0.5
    layer%kz_e = (layer%k_t**2-k_rho**2/layer%nu_e)**0.5
    layer%kz_h = (layer%k_t**2-k_rho**2/layer%nu_h)**0.5
    layer%P_n = (/ EXP(-1*(0.0,1.0)*layer%kz_e*layer%d), 0, 0, EXP(-i*layer%kz_h*layer%d)/)
    layer%Z_e = eta_0 * layer%kz_e / k_0 / layer%eps_t
    layer%Z_h = eta_0 * k_0 * layer%mu_t / layer%kz_h
    end subroutine initalize_layer
  
end module Layer_Class