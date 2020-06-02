module Stacked_Graphene_FR
    use Sim_parameters
    use Substrate
    use Layer_Class
    use Fields_Class
    use GrapheneSig
    use graphene
    use S_Matrix_Class
    implicit none   
    integer, parameter :: n_mat = 1
    complex(wp), allocatable :: eps_array(:)
    real(wp), allocatable :: thickness_array(:)
    integer, allocatable :: g_array(:)
    integer :: p_count = 0, is_g_1, is_g_2
    real(wp) :: inter_thickness, bottom_thickness
    complex(wp) :: inter_mat
    
    contains
    subroutine stack_graphene_fr_config(freq, layers, inc_field, theta_in, xi_in, parameters)
        implicit none
        real(wp), intent(in) :: freq
        type(Layer), allocatable, intent(inout) :: layers(:)
        type(Fields), intent(inout) :: inc_field
        integer :: i, j
        real(wp), intent(in), optional :: theta_in, xi_in
        real(wp), intent(in), optional, dimension(:) :: parameters
        complex(wp), dimension(2,2,2) :: txref_coeff_pc
        
        inter_mat = si_e
        bottom_thickness = 12.0e-6_wp
        !bottom_thickness = 9.6e-6_wp
        
        if (.NOT. allocated(eps_array)) then
            allocate(eps_array(n_mat))
            allocate(g_array(n_mat))
            
            eps_array(1) = si_e
            ! eps_array(2) = si_e
            ! eps_array(2) = polymethylpentene_e
            
            call pc_thickness_config(eps_array,thickness_array,2.0e12_wp,4.0_wp*(real(p_count,wp)+1.0_wp)) 
            !print*,'dieletrics thickness: ', thickness_array(1)
            
            g_array(1) = 1
            !g_array(2) = 0
            !g_array(2) = 0
        end if
        
        ! update paramters relative to freq and allocate array
        call update_freq(freq)
        if (.NOT. allocated(layers)) then
            n_layers =  p_count * n_mat + 2
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
        !b0 = 2.5
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
        
        !is_g_2 = 0
        !if (is_g_2 .EQ. 1) then 
        !    layers(2)=Layer(inter_mat, (1.0_wp,0.0_wp),sigxx, sigyy, sigxy, sigyx, (1.0_wp,0.0_wp), (1.0_wp,0.0_wp),inter_thickness)
        !else
        !    layers(2)=Layer(inter_mat, (1.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(1.0_wp,0.0_wp), (1.0_wp,0.0_wp),inter_thickness)
        !end if
        
        if ( p_count > 0 ) then
        do i = 1, p_count
            do j = 1, n_mat
                if (g_array(j) .EQ. 1) then
                    layers((i-1)*n_mat+1+j)=Layer(eps_array(j), (1.0_wp,0.0_wp), sigxx, sigyy, sigxy, sigyx, (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), thickness_array(j))
                else
                    layers((i-1)*n_mat+1+j)=Layer(eps_array(j), (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), thickness_array(j))
                end if
            end do
        end do
        end if
        
        layers(p_count*n_mat+2)=Layer(si_e, (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), bottom_thickness)
        ! layers(p_count*n_mat+3)=Layer((1.0_wp,-1e12_wp), (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp)
        ! if pec backed
        pec_flag = 1
        
        ! txref_coeff_pc = peak_txref_coeff(layers,3,size(layers),pec_flag,0
        
        
    end subroutine stack_graphene_fr_config
    
    ! compute the thickness of the substrate
    subroutine compute_inter_thickness(freq_in)
        type(Layer), allocatable :: layers_temp(:)
        type(Fields) :: inc_field_temp
        real(wp) :: ref_angle
        complex(wp), dimension(2,2,2) :: txref_coeff_pc
        real(wp), intent(in) :: freq_in
        
        call stack_graphene_fr_config(freq_in, layers_temp, inc_field_temp, 0.0_wp, 0.0_wp)
        txref_coeff_pc = peak_txref_coeff(layers_temp,3,size(layers_temp),1,0)
        ref_angle = atan2(aimag(txref_coeff_pc(2,1,1)),real(txref_coeff_pc(2,1,1)))
        
        if (ref_angle < 0.0_wp) then
            ref_angle = ref_angle + 2*PI
        end if
        
        inter_thickness =(2*PI - ref_angle)/2.0_wp/(k_0*(inter_mat**0.5_wp))
        print*, 'dieletrics thickness: ', inter_thickness
        deallocate(layers_temp)
        
    end subroutine compute_inter_thickness
    
end module Stacked_Graphene_FR