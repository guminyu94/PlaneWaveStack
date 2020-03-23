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
        complex :: eps_t, mu_t, Z_e, Z_h, kz_e, kz_h, k_t, nu_e, nu_h
        complex, dimension(2,2) :: sigma, P_n, Z_n, Y_n  
        real :: d
    end type Layer
    
    interface Layer
        procedure :: initalize_layer
    end interface Layer
    
    contains
    
    ! constr function of Layer type
    ! calculate the rest of paramters based on eps, mu, and sigma
    type(Layer) function  initalize_layer(eps_t_in, mu_t_in, sigma_x_in, sigma_y_in, nu_e_in, nu_h_in, d_in)
        implicit none
        complex, intent(in) :: eps_t_in, mu_t_in, sigma_x_in, sigma_y_in, nu_e_in, nu_h_in
        real, intent(in) :: d_in
    
        initalize_layer%eps_t = eps_t_in
        initalize_layer%mu_t = eps_t_in
        initalize_layer%nu_e = nu_e_in
        initalize_layer%nu_h = nu_h_in
    
        initalize_layer%sigma(1,1) = sigma_x_in
        initalize_layer%sigma(1,2) = (0.0,0.0)
        initalize_layer%sigma(2,1) = (0.0,0.0)
        initalize_layer%sigma(2,2) = sigma_y_in
        initalize_layer%k_t = k_0 * (initalize_layer%eps_t * initalize_layer%mu_t) ** 0.5
        initalize_layer%kz_e = (initalize_layer%k_t**2-k_rho**2/initalize_layer%nu_e)**0.5
        initalize_layer%kz_h = (initalize_layer%k_t**2-k_rho**2/initalize_layer%nu_h)**0.5
        initalize_layer%P_n(1,1) = EXP(-1*(0.0,1.0)*initalize_layer%kz_e*initalize_layer%d)
        initalize_layer%P_n(1,2) = (0.0,0.0)
        initalize_layer%P_n(2,1) = (0.0,0.0)
        initalize_layer%P_n(2,2) = EXP(-1*(0.0,1.0)*initalize_layer%kz_h*initalize_layer%d)
        initalize_layer%Z_e = eta_0 * initalize_layer%kz_e / k_0 / initalize_layer%eps_t
        initalize_layer%Z_h = eta_0 * k_0 * initalize_layer%mu_t / initalize_layer%kz_h
        initalize_layer%Z_n(1,1) = initalize_layer%Z_e
        initalize_layer%Z_n(2,2) = initalize_layer%Z_h
        initalize_layer%Z_n(1,2) = (0.0,0.0)
        initalize_layer%Z_n(2,1) = (0.0,0.0)
        initalize_layer%Y_n(1,1) = initalize_layer%Z_e**-1
        initalize_layer%Y_n(2,2) = initalize_layer%Z_h**-1
        initalize_layer%Y_n(1,2) = (0.0,0.0)
        initalize_layer%Y_n(2,1) = (0.0,0.0)
    end function initalize_layer
  
end module Layer_Class