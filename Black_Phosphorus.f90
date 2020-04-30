module Black_Phosphorus
    use Sim_parameters
    implicit none   
    complex(wp), parameter :: sio2_e = (3.61_wp,0.0_wp), si_e = (11.6964_wp,0.0_wp), ds_e= (11.6964_wp,0.0_wp)
    ! defected Si thickness
    real(wp), parameter :: ds_t = 130E-6_wp
    integer :: count = 0, is_bp

        
    contains
    subroutine black_phosphorus_config(freq_in,layers_in,inc_field)
        use Layer_Class
        use Black_Phosphorus_Sigma
        use Fields_Class
        implicit none   
        real(wp), intent(in) :: freq_in
        type(Layer), allocatable, intent(inout) :: layers_in(:)
        type(Fields), intent(inout) :: inc_field
        integer :: i
        
        complex(wp), dimension(2, 2) :: bp_sigma_mat 
        ! defected Si thickness

        
        ! update paramters relative to freq and allocate array
        call update_freq(freq_in)
        if (.NOT. allocated(layers_in)) then
            n_layers =  2 * 15 + 2 + 1
            allocate(layers_in(n_layers))
        end if
        
        xi = 45.0_wp
        ! normal inc plane wave
        theta = 30.0_wp
        k_rho = k_0*SIN(theta/180.0_wp*PI)
        inc_field  = Fields((1.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp))
        
        ! update constant paramters
        bp_sigma_mat = bp_sigma(freq_in,8.55E17_wp)
        
        is_bp = 1
        
        if (is_bp .EQ. 1) then
            ! first layer is bp and incoming media air
            layers_in(1)=Layer((1.0_wp,0.0_wp), (1.0_wp,0.0_wp), bp_sigma_mat(1,1), bp_sigma_mat(2,2), bp_sigma_mat(1,2), bp_sigma_mat(2,1), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp)
        else
            layers_in(1)=Layer((1.0_wp,0.0_wp), (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp)
        end if
        
        if (count .EQ. 0)   then
            !print*, "simga_bp_arm"
            !print*, bp_sigma_mat(1,1)
            !print*, "simga_bp_zig"
            !print*, bp_sigma_mat(2,2)
            print *, layers_in(1)%sigma_n
            count = 1
        end if 
        
        ! second layer is detected Si
        layers_in(2)=Layer(ds_e, (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), ds_t)
        ! PC 15 pairs (Si and SiO2)
        
        do i = 1, 15
            ! SiO2
            layers_in(i*2+1)=Layer(sio2_e, (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 39.47E-6_wp)
            ! Si
            layers_in(i*2+2)=Layer(si_e, (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 21.93E-6_wp)
        end do
        ! last layer is outcoming media air
        layers_in(33)=Layer(sio2_e, (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp)
        
    end subroutine
end module Black_Phosphorus