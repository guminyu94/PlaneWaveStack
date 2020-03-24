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
!   3. Cascade S matrices with star product and build fields vector
!   
!
!****************************************************************************

program PlaneWaveStack_Minyu

    use Sim_parameters
    use Layer_Class
    use S_Matrix_Class
    use Fields_Class
    implicit none    
    
    integer :: n_layers, i
    complex, allocatable :: eps_t(:), mu_t(:), sigma_x(:), sigma_y(:), nu_e(:), nu_h(:)
    real, allocatable :: d(:)
    real :: freq_in
    type(Layer), allocatable :: layers(:)
    type(S_Matrix), allocatable :: S_Matrice(:)
    type(S_Matrix) :: S_Cascaded
    integer :: isPecBacked
    type(Fields) :: inc_Fields
    
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
    
    ! input frequency   
    print *, "Frequency"
    read (*,*) freq_in
    print *, "k_rho"
    read (*,*) k_rho
    
    ! update constant paramters
    call update_freq(freq_in)
    print *, lambda_0
     
    ! input layered parameters   
    do i = 1, n_layers
        if (i .EQ. 1) then
            print *, i, "th: ", "eps, mu, sigma_x, sigma_y, nu_e, nu_h (anisotrpic ratio), first layer"
            read(*,*) eps_t(i), mu_t(i), sigma_x(i), sigma_y(i), nu_e(i), nu_h(i)
            ! for first and last layer, d is not necessary 
            layers(i)=Layer(eps_t(i), mu_t(i), sigma_x(i), sigma_y(i), nu_e(i), nu_h(i), 0)
        else if (i .EQ. n_layers) then
            print *, i, "th: ", "eps, mu, nu_e, nu_h (anisotrpic ratio), last layer"
            read(*,*) eps_t(i), mu_t(i), nu_e(i), nu_h(i)
            ! readin layers' parameters, and assemble layer obj
            layers(i)=Layer(eps_t(i), mu_t(i), 0, 0, nu_e(i), nu_h(i), 0)
            layers(i)%P_n = unit_matrix
        else
            print *, i, "th: ", "eps, mu, sigma_x, sigma_y, nu_e, nu_h (anisotrpic ratio), thickness"
            read(*,*) eps_t(i), mu_t(i), sigma_x(i), sigma_y(i), nu_e(i), nu_h(i), d(i)
            ! readin layers' parameters, and assemble layer obj
            layers(i)=Layer(eps_t(i), mu_t(i), sigma_x(i), sigma_y(i), nu_e(i), nu_h(i), d(i))            
        end if
    end do
    
    ! asseble the fields vector, pec or free space on the last layer
    
    
   do i = 1, n_layers-1
    ! assemble S Matrix from layer obj
       S_Matrices(i)=S_Matrix(layers(i),layers(i+1))
   end do
   
   ! cascade S Matrices
    S_Cascaded = S_Matrices_Cascaded(S_Matrices)
    
    ! is pec backed ?
    print *, "Input 1 if PEC Backed"
    read(*,*) isPecBacked 
    if (isPecBacked .EQ. 1) then
        
    else
        
        inc_Fields%Field_2 = (/(0.0,0.0),(0.0,0.0)/)
    end if
    
    ! end
    print *, 'End of Program, Type in Any Key to Exit'
    read (*,*)

end program PlaneWaveStack_Minyu
    
