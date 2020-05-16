module BP_FR_Sensing
    use Sim_parameters
    use Substrate
    implicit none   

    integer :: p_count = 20, is_bp, is_bp_1
        
    contains
    subroutine bp_fr_sen_config(freq_in, layers_in, inc_field, theta_in, xi_in, parameters )
        use Layer_Class
        use Black_Phosphorus_Sigma
        use Fields_Class
        implicit none   
        real(wp), intent(in) :: freq_in
        type(Layer), allocatable, intent(inout) :: layers_in(:)
        type(Fields), intent(inout) :: inc_field
        integer :: i
        complex(wp), dimension(2, 2) :: bp_sigma_mat 
        real(wp), intent(in), optional :: theta_in, xi_in
        real(wp), intent(in), optional, dimension(:) :: parameters
        ! defected Si thickness
        
        ! update paramters relative to freq and allocate array
        call update_freq(freq_in)
        if (.NOT. allocated(layers_in)) then
            n_layers =  p_count*2 + 3
            allocate(layers_in(n_layers))
        end if
        
        if (present(xi_in)) then
            xi = xi_in
        else
            xi = 0.0_wp
        end if
        
        if (present(theta_in)) then
            ! normal inc plane wave
            theta = theta_in
        else
            theta = 30.0_wp
        end if
        
        k_rho = k_0 * SIN(theta / 180.0_wp * PI)
        inc_field  = Fields((1.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp))
        
        ! update constant paramters, the n can be adjusted
        bp_sigma_mat = bp_sigma_b0(freq_in,6.55E17_wp,5.0_wp)
        
        is_bp_1 = 0
        ! first layer is bp and incoming media air
        if (is_bp_1 .EQ. 1) then 
            layers_in(1)=Layer((1.0_wp,0.0_wp), (1.0_wp,0.0_wp), bp_sigma_mat(1,1), bp_sigma_mat(2,2), bp_sigma_mat(1,2), bp_sigma_mat(2,1), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp)
        else
            layers_in(1)=Layer((1.0_wp,0.0_wp), (1.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp)
        end if
        
        layers_in(2)=Layer((2.8689_wp,0.0_wp), (1.5_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp),  390.0000e-06_wp)
        
        ! add bp 
        is_bp = 0
        
        do i = 1,p_count
            if (is_bp .EQ. 1) then 
                layers_in(i*2+1)=Layer(sio2_e_1t, (1.0_wp,0.0_wp), bp_sigma_mat(1,1), bp_sigma_mat(2,2), bp_sigma_mat(1,2), bp_sigma_mat(2,1), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp),51.3699e-6_wp)
                !layers_in(i*2+1)=Layer(sic_e, (1.0_wp,0.0_wp),bp_sigma_mat(1,1),bp_sigma_mat(2,2), bp_sigma_mat(1,2), bp_sigma_mat(2,1), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 1.1765e-05_wp/4.0_wp)
                layers_in(i*2+2)=Layer(polymethylpentene_e_1t, (1.0_wp,0.0_wp),bp_sigma_mat(1,1),bp_sigma_mat(2,2), bp_sigma_mat(1,2), bp_sigma_mat(2,1), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 39.4737e-6_wp)
            else
                layers_in(i*2+1)=Layer(sio2_e_1t, (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp),51.3699e-6_wp)
               ! layers_in(i*2+1)=Layer(sic_e, (1.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 1.1765e-05_wp/4.0_wp)
                layers_in(i*2+2)=Layer(polymethylpentene_e_1t, (1.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 39.4737e-6_wp)
            end if
        end do

        layers_in(p_count*2+3)=Layer((1.0_wp,0.0_wp), (1.0_wp,0.0_wp),(0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0E-6_wp)
        
    end subroutine bp_fr_sen_config
end module BP_FR_Sensing