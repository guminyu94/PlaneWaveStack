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
    implicit none
    
    
    integer :: n_layers, i
    complex, allocatable :: eps(:), mu(:), sigma_x(:), sigma_y(:) 
    real, allocatable :: d(:)
    real :: freq_in
    
    
    ! input layered parameters        
    print *, "Number of Layers"
    read (*,*) n_layers
    
    allocate(eps(n_layers))
    allocate(mu(n_layers))
    allocate(sigma_x(n_layers))
    allocate(sigma_y(n_layers))
    allocate(d(n_layers))
    
    ! input frequency   
    print *, "Frequency"
    read (*,*) freq_in
    print *, "k_rho"
    read (*,*) k_rho
    
    ! update constant paramters
    call update_freq(freq_in)
    print *, lambda_0
    
    type(layer), allocatable :: layers(n_layers)
    do i = 1, n_layers
         :: 
        print *, i, "th: ", "eps, mu, sigma_x, sigma_y, thickness"
        read(*,*) eps(i), mu(i), sigma_x(i), sigma_y(i),  d(i)
    end do
    
    ! calculate the paramters of each layer used for assembling matrix, encapsulated in an array of layer obj
    
    
    ! end
    print *, 'End of Program, Type in Any Key to Exit'
    
    read (*,*)

    end program PlaneWaveStack_Minyu
    
!****************************************************************************
!
!   subroutine: CAL_PARAM
!
!   PURPOSE:  calculate paramters Z_n, Y_n, P_n used for S matrix 
!
!****************************************************************************
    
!****************************************************************************
!
!   subroutine: S_Matrix
!
!   PURPOSE:  calculate S matrix of nth layer 
!
!****************************************************************************