Module Otto
    contains
    subroutine otto_config(freq_in,layers,inc_field,theta_in,xi_in,parameters)
        use Sim_parameters
        use Layer_Class
        use GrapheneSig
        use Fields_Class
        implicit none   
        real(wp), intent(in) :: freq_in
        type(Layer), allocatable, intent(inout) :: layers(:)
        type(Fields), intent(inout) :: inc_field
        real(wp), intent(in), optional :: theta_in, xi_in
        real(wp), intent(in), optional, dimension(:) :: parameters
        
        call update_freq(freq_in)
        if (.NOT. allocated(layers)) then
            n_layers = 5
            allocate(layers(n_layers))
        end if
        
        if (present(xi_in)) then
            xi = xi_in
        else
            xi = 45.0_wp
        end if
        
        if (present(theta_in)) then
            ! normal inc plane wave
            theta = theta_in
        else
            theta = 60.0_wp
        end if
        
        inc_field  = Fields((0.5_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp))
        
        ! update constant paramters
        call sigmas(real(freq_in),sig_d,sig_h,n_d,n_h)
        sigxx = CMPLX(sig_d,wp)
        sigyy = CMPLX(sig_d,wp)
        sigyx = CMPLX(sig_h,wp)
        sigxy = CMPLX(-sig_h,wp)
        print*, "simga_d"
        print*, sig_d
        print*, "simga_h"
        print*, sig_h
        k_rho = k_0 * (12.0_wp**0.5_wp) * SIN(theta/180.0_wp*PI)
        
        layers(1) = Layer((12.0_wp,0.0_wp), (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp)
        layers(2) = Layer((1.0_wp,0.0_wp), (1.0_wp,0.0_wp), sigxx, sigyy, sigxy, sigyx, (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 5E-6_wp)
        layers(3) = Layer((4.0_wp,0.0_wp), (1.0_wp,0.0_wp), sigxx, sigyy, sigxy, sigyx, (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 5E-6_wp)
        layers(4) = Layer((2.25_wp,0.0_wp), (1.0_wp,0.0_wp), sigxx, sigyy, sigxy, sigyx, (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 5E-6_wp)
        layers(5) = Layer((2.0_wp,0.0_wp), (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp)
        
    end subroutine otto_config
end module Otto