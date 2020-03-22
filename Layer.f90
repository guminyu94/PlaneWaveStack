module Layer_Class
    use Sim_parameters
    implicit none
    public :: CAL_LAYER_PARAM
    
    type Layer
        complex :: eps, mu, sigma(2,2), Z_e, Z_h, P_n(2,2), kz_e, kz_h
        real :: d
    end type Layer
    
    contains
    subroutine  CAL_LAYER_PARAM(layer)
    implicit none
    ! input
    type(layer), intent(inout)
    layer%sigma = (/ layer%sigma_x_in, 0, 0, layer%sigma_y_in/)
    layer%kz_e = (k_0**2-k_rho**2/eps)**0.5
    layer%kz_h = (k_0**2-k_rho**2/mu)**0.5
    layer%P_n = EXP(/ -i*kz_e*d, 0, 0, -i*kz_h*layer%d)
    layer%Z_e = 
    layer%Z_h =
    end subroutine CAL_LAYER_PARAM
  
end module Layer