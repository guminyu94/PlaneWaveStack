Module Otto
    contains
    subroutine otto_config(freq_in,layers,inc_field)
        use Sim_parameters
        use Layer_Class
        use GrapheneSig
        use Fields_Class
        implicit none   
        real(wp), intent(in) :: freq_in
        type(Layer), allocatable, intent(inout) :: layers(:)
        type(Fields), intent(inout) :: inc_field
        
        call update_freq(freq_in)
        if (.NOT. allocated(layers)) then
            n_layers = 5
            allocate(layers(n_layers))
        end if
        
        xi = 0.0_wp
        theta = 60.0_wp
        
        inc_field  = Fields((0.5_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp))
        
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
        k_rho = k_0*(12.0_wp**0.5_wp) * SIN(theta/180.0_wp*PI)
        
        layers(1) = Layer((12.0_wp,0.0_wp), (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp)
        layers(2) = Layer((1.0_wp,0.0_wp), (1.0_wp,0.0_wp), sigxx, sigyy, sigxy, sigyx, (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 5E-6_wp)
        layers(3) = Layer((4.0_wp,0.0_wp), (1.0_wp,0.0_wp), sigxx, sigyy, sigxy, sigyx, (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 5E-6_wp)
        layers(4) = Layer((2.25_wp,0.0_wp), (1.0_wp,0.0_wp), sigxx, sigyy, sigxy, sigyx, (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 5E-6_wp)
        layers(5) = Layer((2.0_wp,0.0_wp), (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp)
        
    end subroutine otto_config
end module Otto