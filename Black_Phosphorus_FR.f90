module Black_Phosphorus_FR
    use Sim_parameters
    use Substrate
    use GrapheneSig
    implicit none   

    integer :: p_count = 0, is_bp, is_bp_1
        
    contains
    subroutine bp_fr_config(freq_in, layers, inc_field, theta_in, xi_in, parameters )
        use Layer_Class
        use Black_Phosphorus_Sigma
        use Fields_Class
        implicit none   
        real(wp), intent(in) :: freq_in
        type(Layer), allocatable, intent(inout) :: layers(:)
        type(Fields), intent(inout) :: inc_field
        integer :: i
        complex(wp), dimension(2, 2) :: bp_sigma_mat 
        real(wp), intent(in), optional :: theta_in, xi_in
        real(wp), intent(in), optional, dimension(:) :: parameters
        ! defected Si thickness
        
        ! update paramters relative to freq and allocate array
        call update_freq(freq_in)
        if (.NOT. allocated(layers)) then
            n_layers =  p_count * 3 + 3
            allocate(layers(n_layers))
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
            theta = 15.0_wp
        end if
        
        k_rho = k_0 * SIN(theta / 180.0_wp * PI)
        inc_field  = Fields((1.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp))
        
        ! update constant paramters, the n can be adjusted
        !bp_sigma_mat = bp_sigma_b0(freq_in,2.55E17_wp,3.0001_wp
        call sigmas(real(freq_in),sig_d,sig_h,n_d,n_h)
        sigxx = CMPLX(sig_d,wp)
        sigyy = CMPLX(sig_d,wp)
        sigyx = CMPLX(sig_h,wp)
        sigxy = CMPLX(-sig_h,wp)
        bp_sigma_mat(1,1) = sigxx
        bp_sigma_mat(1,2) = sigxy
        bp_sigma_mat(2,1) = sigyx
        bp_sigma_mat(2,2) = sigyy
        
        is_bp_1 = 0
        ! first layer is bp and incoming media air
        if (is_bp_1 .EQ. 1) then 
            layers(1)=Layer((1.0_wp,0.0_wp), (1.0_wp,0.0_wp), bp_sigma_mat(1,1), bp_sigma_mat(2,2), bp_sigma_mat(1,2), bp_sigma_mat(2,1), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp)
        else
            layers(1)=Layer((1.0_wp,0.0_wp), (1.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp)
        end if
        
        layers(2)=Layer(si_e, (1.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp),  0.0001e-06_wp)
        
        ! add bp 
        is_bp = 1
        
        do i = 1,p_count
            if (is_bp .EQ. 1) then 
                layers(i*3)=Layer(si_e, (1.0_wp,0.0_wp), bp_sigma_mat(1,1), bp_sigma_mat(2,2), bp_sigma_mat(1,2), bp_sigma_mat(2,1), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp),   8.7719e-06_wp)
                layers(i*3+1)=Layer(sic_e, (1.0_wp,0.0_wp),bp_sigma_mat(1,1),bp_sigma_mat(2,2), bp_sigma_mat(1,2), bp_sigma_mat(2,1), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 1.4006e-05_wp)
                layers(i*3+2)=Layer(polymethylpentene_e, (1.0_wp,0.0_wp),bp_sigma_mat(1,1),bp_sigma_mat(2,2), bp_sigma_mat(1,2), bp_sigma_mat(2,1), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 2.0690e-05_wp)
            else
                layers(i*3)=Layer(si_e, (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp),   8.7719e-06_wp)
                layers(i*3+1)=Layer(sic_e, (1.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 1.4006e-05_wp)
                layers(i*3+2)=Layer(polymethylpentene_e,(1.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp),  2.0690e-05_wp)
            end if
        end do

        layers(p_count*3+3)=Layer(si_e, (1.0_wp,0.0_wp),(0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0E-6_wp)
        
    end subroutine bp_fr_config
end module Black_Phosphorus_FR