Module Graphene_PC_OM
    contains 
    subroutine graphene_pc_config(freq_in, layers, inc_field, theta_in, xi_in, parameters)
        use Sim_parameters
        use Layer_Class
        use GrapheneSig
        use Fields_Class
        use Substrate
        
        implicit none   
        real(wp), intent(in) :: freq_in
        type(Layer), allocatable, intent(inout) :: layers(:)
        type(Fields), intent(inout) :: inc_field
        integer :: i 
        real(wp), intent(in), optional :: theta_in, xi_in
        real(wp), intent(in), optional, dimension(:) :: parameters
        
        call update_freq(freq_in)
        if (.NOT. allocated(layers)) then
            n_layers = 15
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
            theta = 0.0_wp
        end if
        
        inc_field  = Fields((1.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp))
        
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
        k_rho = k_0 * (1.0_wp**0.5_wp) * SIN(theta/180.0_wp*PI)
        
        ! vaccum with graphene sheet
        layers(1) = Layer((1.0_wp,0.0_wp), (1.0_wp,0.0_wp), sigxx, sigyy, sigxy, sigyx, (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp)
        ! one thin layer of SiC
        layers(2) = Layer(sic_e, (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.52E-6_wp)
        ! PC
        do i=1,6
            layers(1+i*2) = Layer(sic_e, (1.0_wp,0.0_wp),(0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 4.01E-6_wp)
            layers(2+i*2) = Layer(si_e, (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 6.38E-6_wp)
        end do
        ! last layer is vaccum 
        layers(15) = Layer((1.0_wp,0.0_wp), (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp)
        
    end subroutine graphene_pc_config
    
end module Graphene_PC_OM