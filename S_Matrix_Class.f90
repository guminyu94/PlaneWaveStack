!****************************************************************************
!
!   class: S_Matrix
!
!   PURPOSE:  calculate S matrix of nth layer, initalizatizate with nth and n+1 th layer as input
!
!****************************************************************************
Module S_Matrix_Class
    use Sim_parameters, only : eta_0, k_0, k_rho
    use Layer_Class
    use Math
    implicit none
    private
    public :: S_Matrix, S_Matrices_Cascade, print_S_Matrix
    
    ! define S_Matrix type
    Type S_Matrix
        complex(wp), dimension(2, 2) :: alpha_n, gamma_n, delta_n, beta_n 
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
        implicit none 
        Type(Layer), intent(in) :: layer_n_in, layer_n_1_in
        S_Matrix_initalization%alpha_n = MATMUL( -1.0_wp * unit_matrix + 2.0_wp * ( ( unit_matrix + MATMUL(layer_n_in%Z_n , (layer_n_1_in%Y_n + layer_n_in%sigma_n)) )**-1 )  , layer_n_in%P_n )
        S_Matrix_initalization%delta_n = 2.0_wp * MATMUL(( unit_matrix + MATMUL( layer_n_in%Z_n , (layer_n_1_in%Y_n + layer_n_in%sigma_n)))**-1 , layer_n_in%P_n)
        S_Matrix_initalization%gamma_n = MATMUL( unit_matrix - MATMUL((unit_matrix + MATMUL( layer_n_in%Z_n , layer_n_1_in%Y_n + layer_n_in%sigma_n))**-1 , unit_matrix - MATMUL( layer_n_in%Z_n , layer_n_1_in%Y_n - layer_n_in%sigma_n )), layer_n_1_in%P_n)
        S_Matrix_initalization%beta_n =  -1.0_wp * MATMUL( MATMUL( (unit_matrix + MATMUL( layer_n_in%Z_n , layer_n_1_in%Y_n + layer_n_in%sigma_n ))**-1 , unit_matrix - MATMUl( layer_n_in%Z_n,  layer_n_1_in%Y_n - layer_n_in%sigma_n)) , layer_n_1_in%P_n)
    end function S_Matrix_initalization
    
    ! opertor override, star product function for S matrix
    pure Type(S_Matrix) function star_product(S_Matrix_l, S_Matrix_r)
        implicit none
        Type(S_Matrix), intent(in) :: S_Matrix_l, S_Matrix_r   
        star_product%alpha_n = S_Matrix_l%alpha_n + MATMUL( MATMUL( MATMUL( S_Matrix_l%gamma_n , S_Matrix_r%alpha_n) , (unit_matrix - MATMUL( S_Matrix_l%beta_n, S_Matrix_r%alpha_n))**-1) , S_Matrix_l%delta_n)
        star_product%gamma_n = MATMUL( MATMUL( S_Matrix_l%gamma_n , (unit_matrix - MATMUL( S_Matrix_r%alpha_n , S_Matrix_l%beta_n) )**-1 ) , S_Matrix_r%gamma_n )
        star_product%delta_n = MATMUL( MATMUL( S_Matrix_r%delta_n , (unit_matrix - MATMUL( S_Matrix_l%beta_n , S_Matrix_r%alpha_n) )**-1 )  , S_Matrix_l%delta_n )
        star_product%beta_n = S_Matrix_r%beta_n +  MATMUL( MATMUL( MATMUL( S_Matrix_r%delta_n , S_Matrix_l%beta_n ), (unit_matrix - MATMUL( S_Matrix_r%alpha_n , S_Matrix_l%beta_n))**-1), S_Matrix_r%gamma_n)
    end function star_product
    
    ! cascade S matrix
    pure function S_Matrices_Cascade(S_Matrices, start_index, end_index) result(S_matrices_Cascaded)
        implicit none 
        Type(S_Matrix), allocatable, intent(in) :: S_Matrices(:)
        integer, intent(in) :: start_index, end_index
        Type(S_Matrix) :: S_matrices_Cascaded
        integer :: i
        S_matrices_Cascaded = S_Matrices(start_index)
        do i = start_index+1, end_index
            S_matrices_Cascaded = S_matrices_Cascaded * S_Matrices(i)
        end do
    end function
    
    ! print S matrix
    subroutine print_S_matrix(S_matrix_in)
        implicit none 
        Type(S_Matrix), intent(in) :: S_matrix_in
        print*, "Printing S Matrix: "
        print*, "S Matrix alpha: ", S_matrix_in%alpha_n
        print*, "S Matrix gamma: ", S_matrix_in%gamma_n
        print*, "S Matrix delta: ", S_matrix_in%delta_n
        print*, "S Matrix beta: ", S_matrix_in%beta_n
    end subroutine print_S_matrix
    
end Module S_Matrix_Class
    
    