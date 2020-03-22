module Layer_Class
    use Sim_parameters, only : eta_0, k_0, k_rho
    implicit none
    private
    public :: initalize_layer, Layer
    
    type Layer
        complex :: eps_t, mu_t, Z_e, Z_h, kz_e, kz_h, k_t, nu_e, nu_h
        complex, dimension(2,2) :: sigma, P_n  
        real :: d
    end type Layer
    
    contains
    
    ! constr function of Layer type
    ! calculate the rest of paramters based on eps, mu, and sigma
    subroutine  initalize_layer(this,eps_t_in, mu_t_in, sigma_x_in, sigma_y_in, nu_e_in, nu_h_in, d_in)
    implicit none
    type(Layer), intent(inout) :: this
    complex, intent(in) :: eps_t_in, mu_t_in, sigma_x_in, sigma_y_in, nu_e_in, nu_h_in
    real, intent(in) :: d_in
    
    this%eps_t = eps_t_in
    this%mu_t = eps_t_in
    this%nu_e = nu_e_in
    this%nu_h = nu_h_in
    
    this%sigma(1,1) = sigma_x_in
    this%sigma(1,2) = (0.0,0.0)
    this%sigma(2,1) = (0.0,0.0)
    this%sigma(2,2) = sigma_y_in
    this%k_t = k_0 * (this%eps_t * this%mu_t) ** 0.5
    this%kz_e = (this%k_t**2-k_rho**2/this%nu_e)**0.5
    this%kz_h = (this%k_t**2-k_rho**2/this%nu_h)**0.5
    this%P_n(1,1) = EXP(-1*(0.0,1.0)*this%kz_e*this%d)
    this%P_n(1,2) = (0.0,0.0)
    this%P_n(2,1) = (0.0,0.0)
    this%P_n(2,2) = EXP(-1*(0.0,1.0)*this%kz_h*this%d)
    this%Z_e = eta_0 * this%kz_e / k_0 / this%eps_t
    this%Z_h = eta_0 * k_0 * this%mu_t / this%kz_h
    end subroutine initalize_layer
  
end module Layer_Class