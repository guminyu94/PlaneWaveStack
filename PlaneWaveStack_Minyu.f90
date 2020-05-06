!   PlaneWaveStack_Minyu.f90 
!   S Matrix Formation for Stacked Layed Media Exicted by Plane Wave.
!   Minyu Gu, 03/20/2020
!   Intel Fortran Compiler, VS2017
!   Ref[1] Modal transmission line theory of plane wave excited layered media with multiple conductive anisotropic sheets at the interfaces, K.A. Michakski, 2018
!   Ref[2] Electromagnetic field computation in planar multilayers, K.A. Michalski
    
!****************************************************************************
!
!   PROGRAM: PlaneWaveStack_Minyu
!
!   PURPOSE:  Main Program
!   1. Data IO 
!   2. Build layers and S matrix objs 
!   3. Cascade S matrices with star product and compute fields vector
!   
!
!****************************************************************************

Program PlaneWaveStack_Minyu  
    use Sim_parameters
    use Layer_Class
    use GrapheneSig
    use S_Matrix_Class
    use Fields_Class
    use Swapper
    use Modified_Otto
    implicit none    
    
    integer :: use_saved_config
    integer :: isPecBacked
    integer :: isGraphene = 0
    integer :: addSheet = 0
    integer :: i = 0
    
    procedure(fun_temp_angle), pointer :: fun_p
    complex(wp), dimension(2,2,2) :: tx_ref
    
    ! *, "Use Saved configure"
    ! read (*,*) use_saved_config
    use_saved_config = 1
        
    if (use_saved_config .EQ. 1) then
        ! assign config to swapper
        fun_p => mod_otto_config
        ! swap freq
        call theta_swap(fun_p,5.157E12_wp, 0.0_wp, 90.0_wp, 1001, 'data/txangle_bp')
        ! plot field of a single freq
        ! call fields_computation(fun_p,4.0E12_wp, 101, 'data/fields_1t_otto')
        
    else
        print *, "Number of Layers"
        read (*,*) n_layers
        
        allocate(eps_t(n_layers))
        allocate(mu_t(n_layers))
        allocate(sigma_x(n_layers))
        allocate(sigma_y(n_layers))
        allocate(d(n_layers))
        allocate(nu_e(n_layers))
        allocate(nu_h(n_layers))
        allocate(layers(n_layers))
        allocate(S_Matrices(n_layers-1))
        allocate(sigma_xy(n_layers))
        allocate(sigma_yx(n_layers))
    
       ! input frequency   
        print *, "Frequency"
        read (*,*) freq
    
        print *, "theta(degree)"
        read (*,*) theta
        
        print *, "xi(degree)"
        read (*,*) xi

        ! update constant paramters
        call update_freq(freq)
    
        ! Using subroutine to compute graphene sigma tensor ?
        print *, "Input 1 if Using Graphene"
        read(*,*) isGraphene
        if (isGraphene .EQ. 1) then
            call sigmas(real(freq),sig_d,sig_h,n_d,n_h)
            sigxx = CMPLX(sig_d,wp)
            sigyy = CMPLX(sig_d,wp)
            sigyx = CMPLX(sig_h,wp)
            sigxy = CMPLX(-sig_h,wp)
        end if
       
       
       ! input layered parameters
        do i = 1, n_layers
            ! fisrt layer
            if (i .EQ. 1) then
                ! if use graphene
                if (isGraphene .EQ. 1) then
                    print *, i, "th: ", "eps, mu, nu_e, nu_h (anisotropic ratio), first layer"
                    ! for first layer, d is not necessary
                    read(*,*) eps_t(i), mu_t(i), nu_e(i), nu_h(i)
                    ! update k_rho
                    k_rho = k_0 * ((eps_t(i)*mu_t(i))**0.5_wp) * SIN(theta/180.0_wp*PI)
                    print *, "Add graphene ?"
                    read(*,*) addSheet
                    if (addSheet .EQ. 1) then
                        ! add a sheet of graphene
                        layers(i) = Layer(eps_t(i), mu_t(i), sigxx, sigyy, sigxy, sigyx, nu_e(i), nu_h(i), 0.0_wp)
                    else
                       layers(i) = Layer(eps_t(i), mu_t(i), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), nu_e(i), nu_h(i), 0.0_wp)
                    end if
               else
                    ! input sigma 
                    print *, i, "th: ", "eps, mu, sigma_x, sigma_y, sigma_xy, sigma_yx, nu_e, nu_h (anisotropic ratio), first layer"
                    read(*,*) eps_t(i), mu_t(i), sigma_x(i), sigma_y(i), sigma_xy(i), sigma_yx(i), nu_e(i), nu_h(i)
                    ! update k_rho
                    k_rho = k_0 * ((eps_t(i)*mu_t(i))**0.5_wp) * SIN(theta/180.0_wp*PI)
                    ! for first layer, d is not necessarys
                    layers(i) = Layer(eps_t(i), mu_t(i), sigma_x(i), sigma_y(i), sigma_xy(i), sigma_yx(i), nu_e(i), nu_h(i), 0.0_wp)
               end if
            ! last layer
           else if (i .EQ. n_layers) then
               ! is pec backed ?
                print *, "Input 1 if PEC Backed"
                read(*,*) isPecBacked 
                ! if bakced by PEC
                if (isPecBacked .EQ. 1) then
                    print *, i, "th: ", "eps, mu, nu_e, nu_h (anisotropic ratio), thickness"
                    read(*,*) eps_t(i), mu_t(i), nu_e(i), nu_h(i), d(i)
                    ! readin layers' parameters, omit sigma
                    layers(i) = Layer(eps_t(i), mu_t(i), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), nu_e(i), nu_h(i), d(i))
                else
                    print *, i, "th: ", "eps, mu, nu_e, nu_h (anisotropic ratio), last layer"
                    read(*,*) eps_t(i), mu_t(i), nu_e(i), nu_h(i)
                   ! readin layers' parameters, omit sigma and thickness
                    layers(i) = Layer(eps_t(i), mu_t(i), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), nu_e(i), nu_h(i), 0.0_wp)
                end if
            ! layers in the middle
            else
                ! if use graphene
                if (isGraphene .EQ. 1) then
                    print *, i, "th: ", "eps, mu, nu_e, nu_h (anisotropic ratio), thickness"
                    read(*,*) eps_t(i), mu_t(i), nu_e(i), nu_h(i), d(i)
                    print *, "Add graphene ?"
                    read(*,*) addSheet
                    if (addSheet .EQ. 1) then
                        ! add a sheet of graphen
                        layers(i) = Layer(eps_t(i), mu_t(i), sigxx, sigyy, sigxy, sigyx, nu_e(i), nu_h(i), d(i))
                    else
                        layers(i) = Layer(eps_t(i), mu_t(i), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), (0.0_wp,0.0_wp), nu_e(i), nu_h(i), d(i))
                    end if
                else
                    ! input sigma 
                    print *, i, "th: ", "eps, mu, sigma_x, sigma_y, sigma_xy, sigma_yx, nu_e, nu_h (anisotropic ratio), thickness"
                    read(*,*) eps_t(i), mu_t(i), sigma_x(i), sigma_y(i), sigma_xy(i), sigma_yx(i), nu_e(i), nu_h(i), d(i)
                    layers(i) = Layer(eps_t(i), mu_t(i), sigma_x(i), sigma_y(i), sigma_xy(i), sigma_yx(i), nu_e(i), nu_h(i), d(i))
                end if     
            end if
        end do
    end if

    ! end
    print *, 'End of Program, Type in Any Key to Exit'
    read (*,*)

end program PlaneWaveStack_Minyu
    
