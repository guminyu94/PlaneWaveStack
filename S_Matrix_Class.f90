!****************************************************************************
!
!   class: S_Matrix
!
!   PURPOSE:  calculate S matrix of nth layer, initalizatizate with nth and n+1 th layer as input
!
!****************************************************************************
Module S_Matrix_Class
    use Sim_parameters, only : eta_0, k_0, k_rho, unit_matrix
    use Layer_Class
    implicit none
    
    private
    public :: S_Matrix, S_Matrices_Cascade, star_product
    
    ! define S_Matrix type
    Type S_Matrix
        complex, dimension(2, 2) :: alpha_n, gamma_n, delta_n, beta_n 
    end Type S_Matrix
    
    ! S matrix construction function
    interface S_Matrix
        procedure :: S_Matrix_initalization
    end interface S_Matrix
    
    ! opertor override, star product
    interface operator (*)
        procedure :: star_product
    end interface operator (*)
    
    contains
    ! initalization function of S_Matrix class
    Type(S_Matrix) function S_Matrix_initalization(layer_n_in,layer_n_1_in)
        Type(Layer), intent(in) :: layer_n_in, layer_n_1_in
        S_Matrix_initalization%alpha_n = (-1 * unit_matrix + 2 * ( unit_matrix + layer_n_in%Z_n * (layer_n_1_in%Y_n + layer_n_in%sigma_n))**-1 ) * layer_n_in%P_n
        S_Matrix_initalization%delta_n = 2 * ( unit_matrix + layer_n_in%Z_n * (layer_n_1_in%Y_n + layer_n_in%sigma_n))**-1 * layer_n_in%P_n
        S_Matrix_initalization%gamma_n = (unit_matrix - (unit_matrix + layer_n_in%Z_n * (layer_n_1_in%Y_n + layer_n_in%sigma_n))**-1 * (unit_matrix - layer_n_in%Z_n * (layer_n_1_in%Y_n - layer_n_in%sigma_n)) ) * layer_n_1_in%P_n
        S_Matrix_initalization%beta_n =  -1 * (unit_matrix + layer_n_in%Z_n * (layer_n_1_in%Y_n + layer_n_in%sigma_n))**-1 * (unit_matrix - layer_n_in%Z_n * (layer_n_1_in%Y_n - layer_n_in%sigma_n)) * layer_n_1_in%P_n
    end function S_Matrix_initalization
    
    ! opertor override, star product function for S matrix
    Type(S_Matrix) function star_product(S_Matrix_l, S_Matrix_r)
        Type(S_Matrix), intent(in) :: S_Matrix_l, S_Matrix_r
            star_product%alpha_n =  S_Matrix_l%alpha_n + S_Matrix_l%gamma_n * S_Matrix_r%alpha_n * (unit_matrix - S_Matrix_l%beta_n * S_Matrix_r%alpha_n)**-1 * S_Matrix_l%delta_n
            star_product%gamma_n = S_Matrix_l%gamma_n * (unit_matrix - S_Matrix_r%alpha_n * S_Matrix_l%beta_n)**-1 * S_Matrix_r%gamma_n
            star_product%delta_n = S_Matrix_r%delta_n * (unit_matrix - S_Matrix_l%beta_n * S_Matrix_r%alpha_n)**-1 * S_Matrix_l%delta_n
            star_product%beta_n =  S_Matrix_r%beta_n + S_Matrix_r%delta_n * S_Matrix_l%beta_n * (unit_matrix - S_Matrix_r%alpha_n * S_Matrix_l%beta_n)**-1 * S_Matrix_r%gamma_n
    end function star_product
    
    ! cascade S matrix
    Type(S_Matrix) function S_Matrices_Cascade(S_Matrices, start_index, end_index)
        Type(S_Matrix), allocatable, intent(in)
        integer, intent(in) :: start_index, end_index
        S_Matrices_Cascaded = S_Matrices(start_index)
        do i = start_index+1, end_index
            S_Matrices_Cascaded = S_Matrices_Cascaded * S_Matrices(i)
        end do
    end function
end Module S_Matrix_Class
    
    