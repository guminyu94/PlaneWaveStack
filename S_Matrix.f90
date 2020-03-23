!****************************************************************************
!
!   class: S_Matrix
!
!   PURPOSE:  calculate S matrix of nth layer 
!
!****************************************************************************
Module S_Matrix_Class
    use Sim_parameters, only : eta_0, k_0, k_rho
    use Layer_Class
    implicit none
    
    ! define S_Matrix type, contains all n_layers of S matrix
    Type S_Matrix
        complex, dimension(2, 2) :: matrix 
    end Type S_Matrix
    
    interface S_Matrix
        procedure :: S_Matrix_initalization
    end interface S_Matrix
    
    ! initalization function of S_Matrix class
    Type(S_Matrix) function S_Matrix_initalization(layer_n_in,layers_n_1_in)
        Type(layer), intent(in) :: layer_n_in,layers_n_1_in
        alpha = -1 + 2 * ( 1 + layer_n_in%Z_n*() )
        gamma = 
    end function S_Matrix
end Module S_Matrix_Class