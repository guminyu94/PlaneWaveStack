module Stacked_Graphene_FR
    use Sim_parameters
    use Substrate
    implicit none   
    integer, parameter :: n_mat = 2
    complex(wp), allocatable :: eps_array(:)
    real(wp), allocatable :: thickness_array(:)
    integer :: p_count = 1, is_g, is_g_1
    
    contains
    subroutine stack_graphene_fr_config(freq_in, layers_in, inc_field, theta_in, xi_in, parameters)
        use Layer_Class
        use Fields_Class
        use GrapheneSig
        use graphene
        implicit none
        real(wp), intent(in) :: freq_in
        type(Layer), allocatable, intent(inout) :: layers_in(:)
        type(Fields), intent(inout) :: inc_field
        integer :: i, j
        complex(wp), dimension(2, 2) :: bp_sigma_mat
        real(wp), intent(in), optional :: theta_in, xi_in
        real(wp), intent(in), optional, dimension(:) :: parameters
        
        if (.NOT. allocated(eps_array)) then
            allocate(eps_array(n_mat))
            eps_array(1) = si_e
            eps_array(2) = air
            allocate(thickness_array(n_mat))
            call pc_thickness_config(eps_array,thickness_array,2e12_wp)                        
        end if
        
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
            theta = 15.0_wp
        end if
        
        k_rho = k_0 * SIN(theta / 180.0_wp * PI)
        inc_field  = Fields((1.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp))
        
        ! these graphene paramters are imported
        b0 = 0.1
        call sigmas(real(freq_in),sig_d,sig_h,n_d,n_h)
        
        sigxx = CMPLX(sig_d,wp)
        sigyy = CMPLX(sig_d,wp)
        sigyx = CMPLX(sig_h,wp)
        sigxy = CMPLX(-sig_h,wp)
        
        is_g_1 = 0
        ! first layer is bp and incoming media air
        if (is_g_1 .EQ. 1) then 
            layers_in(1)=Layer((1.0_wp,0.0_wp), (1.0_wp,0.0_wp), sigxx, sigyy, sigxy, sigyx, (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp)
        else
            layers_in(1)=Layer((1.0_wp,0.0_wp), (1.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp)
        end if
        
        layers_in(2)=Layer(si_e, (1.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.9745e-08_wp)
        
        ! if add graphene
        is_g = 1
        if (is_g .EQ. 0) then
            sigxx = (0.0_wp,0.0_wp)
            sigyy = (0.0_wp,0.0_wp)
            sigyx = (0.0_wp,0.0_wp)
            sigxy = (0.0_wp,0.0_wp)
        end if
        
        do i = 1, p_count
            do j = 1, n_mat
                layers_in((i-1)*n_mat+2+j)=Layer(eps_array(j), (1.0_wp,0.0_wp), sigxx, sigyy, sigxy, sigyx, (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), thickness_array(j))
            end do
        end do
        
        layers_in(p_count*n_mat+3)=Layer(si_e, (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0E-6_wp)
        layers_in(p_count*n_mat+3)=Layer(si_e, (1.0_wp,0.0_wp),(0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0E-6_wp)
        
    end subroutine stack_graphene_fr_config
end module Stacked_Graphene_FR