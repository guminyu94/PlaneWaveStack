module Stacked_Graphene_FR
    use Sim_parameters
    use Substrate
    implicit none   
    integer, parameter :: n_mat = 2
    complex(wp), allocatable :: eps_array(:)
    real(wp), allocatable :: thickness_array(:)
    integer, allocatable :: g_array(:)
    integer :: p_count = 0, is_g_1
    
    contains
    subroutine stack_graphene_fr_config(freq, layers, inc_field, theta_in, xi_in, parameters)
        use Layer_Class
        use Fields_Class
        use GrapheneSig
        use graphene
        implicit none
        real(wp), intent(in) :: freq
        type(Layer), allocatable, intent(inout) :: layers(:)
        type(Fields), intent(inout) :: inc_field
        integer :: i, j
        real(wp), intent(in), optional :: theta_in, xi_in
        real(wp), intent(in), optional, dimension(:) :: parameters
        complex(wp), dimension(2,2,2) :: txref_coeff_pc
        
        if (.NOT. allocated(eps_array)) then
            allocate(eps_array(n_mat))
            allocate(g_array(n_mat))
            
            eps_array(1) = si_e
            !eps_array(2) = sic_e
            eps_array(2) = polymethylpentene_e
            
            call pc_thickness_config(eps_array,thickness_array,2.0e12_wp,4.0_wp)   
            
            g_array(1) = 0
            !g_array(2) = 0
            g_array(2) = 0
        end if
        
        ! update paramters relative to freq and allocate array
        call update_freq(freq)
        if (.NOT. allocated(layers)) then
            n_layers =  p_count * n_mat + 3
            allocate(layers(n_layers))
        end if
        
        if (present(xi_in)) then
            xi = xi_In
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
        inc_field  = Fields((1.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp))
        
        ! these graphene paramters are imported
        b0 = 0.5
        call sigmas(real(freq),sig_d,sig_h,n_d,n_h)
        
        sigxx = CMPLX(sig_d,wp)
        sigyy = CMPLX(sig_d,wp)
        sigyx = (CMPLX(sig_h,wp))
        sigxy = (CMPLX(-sig_h,wp))
        
        is_g_1 = 1
        ! first layer is bp and incoming media air
        if (is_g_1 .EQ. 1) then 
            layers(1)=Layer((1.0_wp,0.0_wp), (1.0_wp,0.0_wp), sigxx, sigyy, sigxy, sigyx, (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp)
        else
            layers(1)=Layer((1.0_wp,0.0_wp), (1.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp)
        end if
        
        if ( p_count > 0 ) then
        do i = 1, p_count
            do j = 1, n_mat
                if (g_array(j) .EQ. 1) then
                    layers((i-1)*n_mat+2+j)=Layer(eps_array(j), (1.0_wp,0.0_wp), sigxx, sigyy, sigxy, sigyx, (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), thickness_array(j))
                else
                    layers((i-1)*n_mat+2+j)=Layer(eps_array(j), (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), thickness_array(j))
                end if
            end do
        end do
        end if
        
        layers(p_count*n_mat+3)=Layer((1.0_wp,0.0_wp), (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 1e-06_wp)
        
        ! if pec backed
        pec_flag = 1
        
        txref_coeff_pc = peak_txref_coeff(layers,3,size(layers),pec_flag,0)
        print*, 'PC REF: ',txref_coeff_pc(2,1,1),txref_coeff_pc(2,2,1)
        layers(2)=Layer(si_e, (1.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp),pc_wl(si_e,2.0e12_wp,4.0_wp))
    end subroutine stack_graphene_fr_config
end module Stacked_Graphene_FR