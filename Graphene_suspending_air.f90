module Graphene_suspending_air
    use Sim_parameters
    use Substrate
    use Layer_Class
    use Fields_Class
    use GrapheneSig
    use graphene
    use S_Matrix_Class
    implicit none   
    
    contains
    subroutine graphene_air_config(freq, layers, inc_field, theta_in, xi_in, parameters)
        implicit none
        real(wp), intent(in) :: freq
        type(Layer), allocatable, intent(inout) :: layers(:)
        type(Fields), intent(inout) :: inc_field
        integer :: i, j
        real(wp), intent(in), optional :: theta_in, xi_in
        real(wp), intent(in), optional, dimension(:) :: parameters
        complex(wp), dimension(2,2,2) :: txref_coeff_pc
        
        
        ! update paramters relative to freq and allocate array
        call update_freq(freq)
        if (.NOT. allocated(layers)) then
            n_layers =  2
            allocate(layers(n_layers))
        else
            deallocate(layers)
            n_layers =  2
            allocate(layers(n_layers))
        end if
        
        if (present(xi_in)) then
            xi = xi_in
        else
            xi = 0.0_wp
        end if
        
        if (present(theta_In)) then
            ! normal inc plane wave
            theta = theta_in
        else
            theta = 0.0_wp
        end if
        
        k_rho = k_0 * SIN(theta / 180.0_wp * PI)
        inc_field  = Fields((1.0_wp,0.0_wp)*cos(theta/pi*180_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp))
        
        ! these graphene paramters are imported
        b0 = 0.5
        call sigmas(real(freq),sig_d,sig_h,n_d,n_h)
        
        sigxx = CMPLX(sig_d,wp)
        sigyy = CMPLX(sig_d,wp)
        sigyx = (CMPLX(sig_h,wp))
        sigxy = (CMPLX(-sig_h,wp))

        
        layers(1)=Layer((1.0_wp,0.0_wp), (1.0_wp,0.0_wp), sigxx, sigyy, sigxy, sigyx, (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp)
       
        
        layers(2)=Layer(air_e, (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp)

        pec_flag = 0
        
    end subroutine graphene_air_config
    

    
end module Graphene_suspending_air