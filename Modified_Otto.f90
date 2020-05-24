Module Modified_Otto
    use Sim_parameters
    implicit none
    integer :: if_graphene
    
    contains 
    subroutine mod_otto_config(freq_in,layers,inc_field,theta_in,xi_in,parameters)
        use Layer_Class
        use Graphene_Drude
        !use GrapheneSig
        use Fields_Class
        use Substrate
        implicit none   
        real(wp), intent(in) :: freq_in
        type(Layer), allocatable, intent(inout) :: layers(:)
        type(Fields), intent(inout) :: inc_field
        real(wp), intent(in), optional :: theta_in, xi_in
        real(wp), intent(in), optional, dimension(:) :: parameters
        complex(wp) :: sig_d
        complex(wp) :: field_cos
        
        call update_freq(freq_in)
        if (.NOT. allocated(layers)) then
            n_layers = 4
            allocate(layers(n_layers))
        end if
        
        if( present(xi_in) ) then
            xi = xi_in
        else
            xi = 0.0_wp
        end if
    
        if( present(theta_in) ) then
            theta = theta_in
        else
            theta = 0.0_wp
        end if
        
        field_cos =  CMPLX(COS(theta/180.0_wp * PI),wp)
        ! field_cos = (1.0_wp, 0.0_wp)
        inc_field  = Fields(field_cos,(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp))
        
        ! update constant paramters
        sig_d = graphene_Drude_sig(freq_in)
        !call sigmas(real(freq_in),sig_d,sig_h,n_d,n_h)
        ! print *, sig_d
        
        !print*, "simga_d * eta_0"
        !print*, sig_d * eta_0
        !print*, "simga_h * jeta_0"
        !print*, sig_h*eta_0*(0.0,1.0)
        k_rho = k_0 * (prism_e**0.5_wp) * SIN(theta / 180.0_wp * PI)
        
        if_graphene = 1
        layers(1)=Layer( prism_e, (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp )
        
        if (if_graphene .EQ. 1) then
            layers(2)=Layer((1.0_wp,0.0_wp), (1.0_wp,0.0_wp), cmplx(sig_d,wp), cmplx(sig_d,wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 1E-6_wp)
        else 
            layers(2)=Layer((1.0_wp,0.0_wp), (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 1E-6_wp)
        end if
        
        layers(3)=Layer( polymethylpentene_e, (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 10E-6_wp)
        layers(4)=Layer( prism_e, (1.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (1.0_wp,0.0_wp), (1.0_wp,0.0_wp), 0.0_wp )
        
        
    end subroutine
    
    
end module modified_Otto