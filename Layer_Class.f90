!****************************************************************************
!
!   class: Layer_Class
!
!   PURPOSE:  calculate parameters of nth layer
!
!****************************************************************************
Module Layer_Class
    use Sim_parameters, only : eta_0, k_0, k_rho, wp
    use Math, only : PI
    implicit none    
    
    private
    public :: Layer, print_layer
    
    type Layer
        complex(wp) :: eps_t, mu_t, nu_e, nu_h, Z_e, Z_h, kz_e, kz_h, k_t
        complex(wp), dimension(2,2) :: sigma_n, P_n, Z_n, Y_n  
        real(wp) :: d
    end type Layer
    
    interface Layer
        procedure :: initalize_layer
    end interface Layer
    
    contains
    
    ! constr function of Layer type
    ! calculate the rest of paramters based on eps, mu, and sigma
    function  initalize_layer(eps_t_in, mu_t_in, sigma_x_in, sigma_y_in, nu_e_in, nu_h_in, d_in) result(new_layer)
        implicit none
        complex(wp), intent(in) :: eps_t_in, mu_t_in, sigma_x_in, sigma_y_in, nu_e_in, nu_h_in
        real(wp), intent(in) :: d_in
        type(Layer) :: new_layer
    
        new_layer%eps_t = eps_t_in
        new_layer%mu_t = mu_t_in
        new_layer%nu_e = nu_e_in
        new_layer%nu_h = nu_h_in
        new_layer%d = d_in
        
        new_layer%sigma_n = sigma_n_uv(sigma_x_in, sigma_y_in, (0.0_wp,0.0_wp), (0.0_wp,0.0_wp))
       
        new_layer%k_t = k_0 * ((new_layer%eps_t * new_layer%mu_t)**0.5_wp)
        new_layer%kz_e = (new_layer%k_t**2.0_wp-(k_rho**2.0_wp)/new_layer%nu_e)**0.5_wp
        new_layer%kz_h = (new_layer%k_t**2.0_wp-(k_rho**2.0_wp)/new_layer%nu_h)**0.5_wp
        new_layer%P_n(1,1) = EXP(-1*(0.0_wp,1.0_wp)*new_layer%kz_e*new_layer%d)
        new_layer%P_n(1,2) = (0.0_wp,0.0_wp)
        new_layer%P_n(2,1) = (0.0_wp,0.0_wp)
        new_layer%P_n(2,2) = EXP(-1.0_wp*(0.0_wp,1.0_wp)*new_layer%kz_h*new_layer%d)
        
        new_layer%Z_e = eta_0 * new_layer%kz_e / k_0 / new_layer%eps_t
        new_layer%Z_h = eta_0 * k_0 * new_layer%mu_t / new_layer%kz_h
        new_layer%Z_n(1,1) = new_layer%Z_e
        new_layer%Z_n(2,2) = new_layer%Z_h
        new_layer%Z_n(1,2) = (0.0_wp,0.0_wp)
        new_layer%Z_n(2,1) = (0.0_wp,0.0_wp)
        
        new_layer%Y_n(1,1) = new_layer%Z_e**-1.0_wp
        new_layer%Y_n(2,2) = new_layer%Z_h**-1.0_wp
        new_layer%Y_n(1,2) = (0.0_wp,0.0_wp)
        new_layer%Y_n(2,1) = (0.0_wp,0.0_wp)
    end function initalize_layer
    
    subroutine print_layer(layer_in)
         implicit none
         type(Layer), intent(in) :: layer_in
         print *, "Printing current layer's paramter, ", "thickness: ", layer_in%d
         print *, "eps_t: ", layer_in%eps_t, ", mu_t: ", layer_in%mu_t, ", nu_e: ", layer_in%nu_e, ", nu_h: ", layer_in%nu_h
         print *, "k_t: ", layer_in%k_t, ", Z_e: ", layer_in%Z_e, ", Z_h: ", layer_in%Z_h, ", kz_e: ", layer_in%kz_e, ", kz_h: ", layer_in%kz_h
         print *, "sigma: ", layer_in%sigma_n
         print *, "P_n: ", layer_in%P_n
         print *, "Z_n: ", layer_in%Z_n
         print *, "Y_n: ", layer_in%Y_n  
    end subroutine
    
    function sigma_n_uv(sigma_xx_in, sigma_yy_in, sigma_xy_in, sigma_yx_in) result(sigma_n)
        complex(wp), intent(in) :: sigma_xx_in, sigma_yy_in, sigma_xy_in, sigma_yx_in
        complex(wp), dimension(2,2) :: sigma_n
        real(wp) :: eta, sineta, coseta, sinetasqr, cosetasqr
        if (REAL(k_rho) .EQ. 0.0_wp) then 
            eta = 0.0_wp
        else if (REAL(k_rho) > 0.0_wp) then 
            eta = ATAN(AIMAG(k_rho) / REAL(k_rho) )
        else
            eta = PI - ATAN(AIMAG(k_rho) / REAL(k_rho) )
        end if
        sineta = SIN(eta)
        coseta = COS(eta)
        sinetasqr = sineta**2.0_wp
        cosetasqr = coseta**2.0_wp
        sigma_n(1,1) = sigma_xx_in * cosetasqr + sigma_yy_in * sinetasqr + (sigma_xy_in + sigma_yx_in) * coseta * sineta
        sigma_n(1,2) = sigma_xy_in * cosetasqr - sigma_yx_in * sinetasqr - (sigma_xx_in - sigma_yy_in) * coseta * sineta
        sigma_n(2,1) = sigma_yx_in * cosetasqr - sigma_xy_in * sinetasqr - (sigma_xx_in - sigma_yy_in) * coseta * sineta
        sigma_n(2,2) = sigma_yy_in * cosetasqr + sigma_xx_in * sinetasqr - (sigma_xy_in + sigma_yx_in) * coseta * sineta
    end function sigma_n_uv
    
    
end module Layer_Class