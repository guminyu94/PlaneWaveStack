module Black_Phosphorus_FR
    use Sim_parameters
    use Substrate, only : sio2_e, si_e, sic_e
    implicit none   
    complex(wp), parameter :: ds_e= (11.6964_wp,0.0_wp)
    ! defected Si thickness
    real(wp), parameter :: ds_t = 130E-6_wp
    integer :: p_count = 15, is_bp, is_bp_1
        
    contains
    subroutine black_phosphorus_fr_config(freq_in, layers_in, inc_field, theta_in, xi_in, parameters )
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
            n_layers =  p_count*2 + 2
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
            theta = 0.0_wp
        end if
        
        k_rho = k_0 * SIN(theta / 180.0_wp * PI)
        inc_field  = Fields((1.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp))
        
        ! update constant paramters, the n can be adjusted
        bp_sigma_mat = bp_sigma_b0(freq_in,2.55E17_wp,5.0_wp)
        
        is_bp_1 = 1
        ! first layer is bp and incoming media air
        if (is_bp_1 .EQ. 1) then 
            layers_in(1)=Layer((1.0_wp,0.0_wp), (1.0_wp,0.0_wp), bp_sigma_mat(1,1), bp_sigma_mat(2,2), bp_sigma_mat(1,2), bp_sigma_mat(2,1), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp)
        else
            layers_in(1)=Layer((1.0_wp,0.0_wp), (1.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp)
        end if
        
        ! add bp 
        is_bp = 0
        
        do i = 1,p_count
            if (is_bp .EQ. 1) then 
                layers_in(i*2)=Layer(si_e, (1.0_wp,0.0_wp), bp_sigma_mat(1,1), bp_sigma_mat(2,2), bp_sigma_mat(1,2), bp_sigma_mat(2,1), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp),  8.7719e-06_wp/4.0_wp)
                layers_in(i*2+1)=Layer(sic_e, (1.0_wp,0.0_wp),bp_sigma_mat(1,1),bp_sigma_mat(2,2), bp_sigma_mat(1,2), bp_sigma_mat(2,1), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 1.1765e-05_wp/4.0_wp)
                !layers_in(i*2+1)=Layer(sio2_e, (1.0_wp,0.0_wp),bp_sigma_mat(1,1),bp_sigma_mat(2,2), bp_sigma_mat(1,2), bp_sigma_mat(2,1), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 1.5789e-05_wp/4.0_wp)
            else
                layers_in(i*2)=Layer(si_e, (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp),  8.7719e-06_wp/4.0_wp)
                layers_in(i*2+1)=Layer(sic_e, (1.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 1.1765e-05_wp/4.0_wp)
                !layers_in(i*2+1)=Layer(sio2_e, (1.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 1.5789e-05_wp/4.0_wp)
            end if
        end do

        layers_in(p_count*2+2)=Layer((1.0_wp,0.0_wp), (1.0_wp,0.0_wp),(0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0E-6_wp)
        
    end subroutine black_phosphorus_fr_config
end module Black_Phosphorus_FR