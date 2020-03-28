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
    public :: S_Matrix, S_Matrices_Cascade, star_product, print_S_Matrix
    
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
    
    ! override the diagonal matirx's inversion operation
    interface operator (**)
        procedure :: diagonal_inverse
    end interface operator (**)
    
    contains
    ! initalization function of S_Matrix class
    Type(S_Matrix) function S_Matrix_initalization(layer_n_in,layer_n_1_in)
        implicit none 
        Type(Layer), intent(in) :: layer_n_in, layer_n_1_in
        S_Matrix_initalization%alpha_n = (-1.0 * unit_matrix + 2.0 * ( unit_matrix + layer_n_in%Z_n * (layer_n_1_in%Y_n + layer_n_in%sigma_n))**-1.0 ) * layer_n_in%P_n
        print*, layer_n_in%Z_n * (layer_n_1_in%Y_n + layer_n_in%sigma_n)
        S_Matrix_initalization%delta_n = 2.0 * ( unit_matrix + layer_n_in%Z_n * (layer_n_1_in%Y_n + layer_n_in%sigma_n))**-1.0 * layer_n_in%P_n
        S_Matrix_initalization%gamma_n = (unit_matrix - (unit_matrix + layer_n_in%Z_n * (layer_n_1_in%Y_n + layer_n_in%sigma_n))**-1 * (unit_matrix - layer_n_in%Z_n * (layer_n_1_in%Y_n - layer_n_in%sigma_n)) ) * layer_n_1_in%P_n
        S_Matrix_initalization%beta_n =  -1.0 * (unit_matrix + layer_n_in%Z_n * (layer_n_1_in%Y_n + layer_n_in%sigma_n))**-1.0 * (unit_matrix - layer_n_in%Z_n * (layer_n_1_in%Y_n - layer_n_in%sigma_n)) * layer_n_1_in%P_n
    end function S_Matrix_initalization
    
    ! opertor override, star product function for S matrix
    Type(S_Matrix) function star_product(S_Matrix_l, S_Matrix_r)
        implicit none
        Type(S_Matrix), intent(in) :: S_Matrix_l, S_Matrix_r   
        star_product%alpha_n =  S_Matrix_l%alpha_n + S_Matrix_l%gamma_n * S_Matrix_r%alpha_n * (unit_matrix - S_Matrix_l%beta_n * S_Matrix_r%alpha_n)**-1 * S_Matrix_l%delta_n
        star_product%gamma_n = S_Matrix_l%gamma_n * (unit_matrix - S_Matrix_r%alpha_n * S_Matrix_l%beta_n)**-1 * S_Matrix_r%gamma_n
        star_product%delta_n = S_Matrix_r%delta_n * (unit_matrix - S_Matrix_l%beta_n * S_Matrix_r%alpha_n)**-1 * S_Matrix_l%delta_n
        star_product%beta_n =  S_Matrix_r%beta_n + S_Matrix_r%delta_n * S_Matrix_l%beta_n * (unit_matrix - S_Matrix_r%alpha_n * S_Matrix_l%beta_n)**-1 * S_Matrix_r%gamma_n
    end function star_product
    
    ! cascade S matrix
    Type(S_Matrix) function S_Matrices_Cascade(S_Matrices, start_index, end_index)
        implicit none 
        Type(S_Matrix), allocatable, intent(in) :: S_Matrices(:)
        integer, intent(in) :: start_index, end_index
        integer :: i
        S_Matrices_Cascade = S_Matrices(start_index)
        do i = start_index+1, end_index
            S_Matrices_Cascade = S_Matrices_Cascade * S_Matrices(i)
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
    
    ! diagonal matrix inversion 
    function diagonal_inverse(matrix_in, power) result(diagonal_inversed_matrix)
        implicit none 
        complex, dimension(2,2), intent(in) :: matrix_in
        complex, dimension(2,2)  :: diagonal_inversed_matrix
        real, intent(in) :: power
        print *, "diagonal matrix inversion"
            ! confirm it's a diagonal matrix
        if ((matrix_in(1,2) .EQ. (0.0,0.0)) .AND. (matrix_in(2,1) .EQ. (0.0,0.0))) then
            diagonal_inversed_matrix(1,1) = matrix_in(1,1)**power
            diagonal_inversed_matrix(1,2) = 0
            diagonal_inversed_matrix(2,1) = 0
            diagonal_inversed_matrix(2,2) = matrix_in(2,2)**power
        else
            diagonal_inversed_matrix(1,1) = matrix_in(1,1)**power
            diagonal_inversed_matrix(1,2) = matrix_in(1,2)**power
            diagonal_inversed_matrix(2,1) = matrix_in(2,1)**power
            diagonal_inversed_matrix(2,2) = matrix_in(2,2)**power
        end if
    end function diagonal_inverse
        
end Module S_Matrix_Class
    
    