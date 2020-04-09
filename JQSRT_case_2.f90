module JQSRT_case_2
    contains
    subroutine JQSRT_case_2_config(freq_in,layers_in,inc_field)
        use Sim_parameters
        use Layer_Class
        use GrapheneSig
        use Fields_Class
        implicit none   
        real(wp), intent(in) :: freq_in
        type(Layer), allocatable, intent(inout) :: layers_in(:)
        type(Fields), intent(inout) :: inc_field
        
        call update_freq(freq_in)
        if (.NOT. allocated(layers_in)) then
            n_layers = 6
            allocate(layers_in(n_layers))
        end if
        
        xi = 0.0_wp
        ! normal inc plane wave
        theta = 0.0_wp
        
        ! update constant paramters
        call sigmas(real(freq_in),sig_d,sig_h,n_d,n_h)
        sigxx = CMPLX(sig_d,wp)
        sigyy = CMPLX(sig_d,wp)
        sigyx = CMPLX(sig_h,wp)
        sigxy = CMPLX(-sig_h,wp)
        print*, "simga_d * eta_0"
        print*, sig_d*eta_0
        print*, "simga_h * jeta_0"
        print*, sig_h*eta_0*(0.0,1.0)
        k_rho = k_0*SIN(theta/180.0_wp*PI)
        
        layers_in(1)=Layer((1.0_wp,0.0_wp), (1.0_wp,0.0_wp), sigxx, sigyy, sigxy, sigyx, (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp)
        layers_in(2)=Layer((10.2_wp,0.0_wp), (1.0_wp,0.0_wp), sigxx, sigyy, sigxy, sigyx, (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 10E-6_wp)
        layers_in(3)=Layer((10.2_wp,0.0_wp), (1.0_wp,0.0_wp), sigxx, sigyy, sigxy, sigyx, (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 10E-6_wp)
        layers_in(4)=Layer((10.2_wp,0.0_wp), (1.0_wp,0.0_wp), sigxx, sigyy, sigxy, sigyx, (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 10E-6_wp)
        layers_in(5)=Layer((10.2_wp,0.0_wp), (1.0_wp,0.0_wp), sigxx, sigyy, sigxy, sigyx, (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 10E-6_wp)
        layers_in(6)=Layer((1.0_wp,0.0_wp), (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp)
        
    end subroutine
end module JQSRT_case_2