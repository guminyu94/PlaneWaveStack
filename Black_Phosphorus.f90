module Black_Phosphorus
    contains
    subroutine black_phosphorus_config(freq_in,layers_in,inc_field)
        use Sim_parameters
        use Layer_Class
        use Black_Phosphorus_Sigma
        use Fields_Class
        implicit none   
        real(wp), intent(in) :: freq_in
        type(Layer), allocatable, intent(inout) :: layers_in(:)
        type(Fields), intent(inout) :: inc_field
        real(wp), parameter :: sio2_e = 1.9_wp**2.0_wp, si_e = 3.42_wp**2.0_wp
        
        complex(wp), dimension(2, 2) :: bp_sigma_mat 
        integer :: i
        ! defected Si thickness
        real(wp) :: ds_t = 130E-6_wp
        ! defected Si thickness
        complex(wp), parameter :: ds_e = 3.42_wp**2.0_wp
        
        ! update paramters relative to freq and allocate array
        call update_freq(freq_in)
        if (.NOT. allocated(layers_in)) then
            n_layers = 15 * 2 + 2 + 1
            allocate(layers_in(n_layers))
        end if
        
        xi = 0.0_wp
        ! normal inc plane wave
        theta = 0.0_wp
        
        ! update constant paramters
        bp_sigma_mat = bp_sigma(freq_in)
        !print*, "simga_bp_arm"
        !print*, bp_sigma_mat(1,1)
        !print*, "simga_bp_zig"
        !print*, bp_sigma_mat(2,2)

        k_rho = k_0*SIN(theta/180.0_wp*PI)

        ! first layer is bp and incoming media air
        layers_in(1)=Layer((1.0_wp,0.0_wp), (1.0_wp,0.0_wp), bp_sigma_mat(1,1), bp_sigma_mat(2,2), bp_sigma_mat(1,2), bp_sigma_mat(2,1), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp)
        ! second layer is detected Si
        layers_in(2)=Layer((1.0_wp,0.0_wp), (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), ds_t)
        ! PC 15 pairs (Si and SiO2)
        
        do i = 1, 15
            ! SiO2
            layers_in(i*2+1)=Layer((sio2_e,0.0_wp), (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 39.47E-6_wp)
            ! Si
            layers_in(i*2+2)=Layer((si_e,0.0_wp), (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 21.93E-6_wp)
        end do
        ! last layer is outcoming media air
        layers_in(33)=Layer((1.0_wp,0.0_wp), (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp)
        
    end subroutine
end module Black_Phosphorus