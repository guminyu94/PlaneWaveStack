!   PlaneWaveStack_Minyu.f90 
!   S and T Matrix Formation for Stacked Layed Media Exicted by Plane Wave.
!   Minyu Gu, 03/20/2020
!   Intel Fortran Compiler, VS2017
    
    
!****************************************************************************
!
!   PROGRAM: PlaneWaveStack_Minyu
!
!   PURPOSE:  Main Program
!   1. Data IO 
!   2. Call Subroutine - CAL_PARAM, to calculate paramters Z_n, Y_n, P_n used for S matrix 
!   3. 
!   
!****************************************************************************

    program PlaneWaveStack_Minyu

    use Sim_parameters
    use Layer_Class
    implicit none
    
    
    integer :: n_layers, i
    complex, allocatable :: eps_t(:), mu_t(:), sigma_x(:), sigma_y(:), nu_e(:), nu_h(:)
    real, allocatable :: d(:)
    real :: freq_in
    type(Layer), allocatable :: layers(:)
    
    ! input layered parameters        
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
    allocate(S_Matrice(n_layers))
    
    ! input frequency   
    print *, "Frequency"
    read (*,*) freq_in
    print *, "k_rho"
    read (*,*) k_rho
    
    ! update constant paramters
    call update_freq(freq_in)
    print *, lambda_0
    
    do i = 1, n_layers
        print *, i, "th: ", "eps, mu, sigma_x, sigma_y, nu_e, nu(h) (anisotrpic ratio), thickness"
        read(*,*) eps_t(i), mu_t(i), sigma_x(i), sigma_y(i), nu_e(i), nu_h(i), d(i)
        ! readin layers' parameters, and assemble layer obj
        layers(i)=Layer(eps_t(i), mu_t(i), sigma_x(i), sigma_y(i), nu_e(i), nu_h(i), d(i))
    end do
    
    
    do i = 1, n_layers
        ! assemble S Matrix
        S_Matrice(i)=S_Matrix(layers(i),layers(i+1))
    end do
    
    ! calculate the paramters of each layer used for assembling matrix, encapsulated in an array of layer obj
    
    
    ! end
    print *, 'End of Program, Type in Any Key to Exit'
    
    read (*,*)

    end program PlaneWaveStack_Minyu
    

    
