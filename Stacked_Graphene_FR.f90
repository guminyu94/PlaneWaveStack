module Stacked_Graphene_FR
    use Sim_parameters
    use Substrate
    implicit none   

    integer :: p_count = 6, is_bp, is_bp_1
        
    contains
    subroutine stack_graphene_fr_config(freq_in, layers_in, inc_field, theta_in, xi_in, parameters )
        use Layer_Class
        use Fields_Class
        use GrapheneSig
        use graphene
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
            n_layers =  p_count * 3 + 3
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
            theta = 45.0_wp
        end if
        
        k_rho = k_0 * SIN(theta / 180.0_wp * PI)
        inc_field  = Fields((1.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp))
        
        ! these graphene paramters are imported
        b0 = 5.0
        call sigmas(real(freq_in),sig_d,sig_h,n_d,n_h)
        sigxx = CMPLX(sig_d,wp)
        sigyy = CMPLX(sig_d,wp)
        sigyx = CMPLX(sig_h,wp)
        sigxy = CMPLX(-sig_h,wp)
        
        is_bp_1 = 0
        ! first layer is bp and incoming media air
        if (is_bp_1 .EQ. 1) then 
            layers_in(1)=Layer((1.0_wp,0.0_wp), (1.0_wp,0.0_wp), sigxx, sigyy, sigxy, sigyx, (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp)
        else
            layers_in(1)=Layer((1.0_wp,0.0_wp), (1.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp)
        end if
        
        layers_in(2)=Layer(si_e, (1.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 7.9745e-06_wp)
        
        ! add bp 
        is_bp = 1
        
        do i = 1,p_count
            if (is_bp .EQ. 1) then 
                layers_in(i*3)=Layer(si_e, (1.0_wp,0.0_wp), sigxx, sigyy, sigxy, sigyx, (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 8.7719e-06_wp)
                layers_in(i*3+1)=Layer(sic_e, (1.0_wp,0.0_wp),sigxx, sigyy, sigxy, sigyx, (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 1.4006e-05_wp)
                layers_in(i*3+2)=Layer(polymethylpentene_e, (1.0_wp,0.0_wp),sigxx, sigyy, sigxy, sigyx, (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 2.0690e-05_wp)
            else
                layers_in(i*3)=Layer(si_e, (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp),   8.7719e-06_wp)
                layers_in(i*3+1)=Layer(sic_e, (1.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 1.4006e-05_wp)
                layers_in(i*3+2)=Layer(polymethylpentene_e,(1.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp),  2.0690e-05_wp)
            end if
        end do

        layers_in(p_count*3+3)=Layer(si_e, (1.0_wp,0.0_wp),(0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0E-6_wp)
        
    end subroutine stack_graphene_fr_config
end module Stacked_Graphene_FR