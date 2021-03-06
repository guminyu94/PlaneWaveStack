module Black_Phosphorus
    use Sim_parameters
    use Substrate, only : sio2_e, si_e
    implicit none   
    complex(wp), parameter :: ds_e= (11.6964_wp,0.0_wp)
    ! defected Si thickness
    real(wp), parameter :: ds_t = 130E-6_wp
    integer :: count = 15, is_bp
        
    contains
    subroutine black_phosphorus_config(freq_in, layers_in, inc_field, theta_in, xi_in, parameters )
        use Layer_Class
        use Black_Phosphorus_Sigma
        use Fields_Class
        implicit none   
        real(wp), intent(in) :: freq_in
        type(Layer), allocatable, intent(inout) :: layers_in(:)
        type(Fields), intent(inout) :: inc_field
        integer :: i
        real(wp), intent(in), optional :: theta_in, xi_in
        complex(wp), dimension(2, 2) :: bp_sigma_mat 
        real(wp), intent(in), optional, dimension(:) :: parameters
        ! defected Si thickness
        
        ! update paramters relative to freq and allocate array
        call update_freq(freq_in)
        if (.NOT. allocated(layers_in)) then
            n_layers =  2 * 15 + 2 + 1
            allocate(layers_in(n_layers))
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
        
        k_rho = k_0 * SIN(theta/180.0_wp*PI)
        inc_field  = Fields((1.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp))
        
        ! update constant paramters, the n can be adjusted
        bp_sigma_mat = bp_sigma(freq_in,6.55E17_wp)

        is_bp = 0
        
        if (is_bp .EQ. 1) then
            ! first layer is bp and incoming media air
            layers_in(1)=Layer((1.0_wp,0.0_wp), (1.0_wp,0.0_wp), bp_sigma_mat(1,1), bp_sigma_mat(2,2), bp_sigma_mat(1,2), bp_sigma_mat(2,1), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp)
        else
            layers_in(1)=Layer((1.0_wp,0.0_wp), (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp)
        end if
        
        ! second layer is defected Si
        layers_in(2) = Layer(ds_e, (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), ds_t)
        ! PC 15 pairs (Si and SiO2)
        
        do i = 1, 15
            ! SiO2
            layers_in(i*2+1) = Layer(sio2_e, (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 39.47E-6_wp)
            ! Si
            layers_in(i*2+2) = Layer(si_e, (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 21.93E-6_wp)
        end do
        
        ! last layer is outcoming media air
        layers_in(33) = Layer(sio2_e, (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp)
        
    end subroutine
end module Black_Phosphorus