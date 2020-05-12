Module Faraday_Rot
    contains
    subroutine faraday_rot_config(freq_in,layers,inc_field)
        use Sim_parameters
        use Layer_Class
        use GrapheneSig
        use Fields_Class
        use Substrate
        implicit none   
        real(wp), intent(in) :: freq_in
        type(Layer), allocatable, intent(inout) :: layers(:)
        type(Fields), intent(inout) :: inc_field
        
        call update_freq(freq_in)
        if (.NOT. allocated(layers)) then
            n_layers = 2
            allocate(layers(n_layers))
        end if
        
        xi = 0.0_wp
        theta = 0.0_wp
        
        inc_field  = Fields((1.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp))
        
        ! update constant paramters
        call sigmas(real(freq_in),sig_d,sig_h,n_d,n_h)
        sigxx = CMPLX(sig_d,wp)
        sigyy = CMPLX(sig_d,wp)
        sigyx = CMPLX(sig_h,wp)
        sigxy = CMPLX(-sig_h,wp)
        !print*, "simga_d * eta_0"
        !print*, sig_d*eta_0
        !print*, "simga_h * jeta_0"
        !print*, sig_h*eta_0*(0.0,1.0)
        k_rho = k_0*(1.0_wp**0.5_wp) * SIN(theta/180.0_wp*PI)
        
        layers(1)=Layer((1.0_wp,0.0_wp), (1.0_wp,0.0_wp), sigxx, sigyy, sigxy, sigyx, (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp)
        layers(2)=Layer(sic_h6_e, (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), sic_h6_nu, (1.0_wp,0.0_wp), 0.0_wp)
        
        
    end subroutine faraday_rot_config
end module Faraday_Rot