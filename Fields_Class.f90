!****************************************************************************
!
!   class: Fields
!
!   PURPOSE:  contains 4 fields components, forward, backward, e and h
!
!****************************************************************************
Module Fields_Class
    use S_Matrix_Class
    use Math, only : unit_matrix
    use Sim_parameters, only : wp
    implicit none
    private
    public :: Fields, trans_ref_coeff_freespace, trans_ref_coeff_freespace_pec, cal_Fields, S_Matrix_Fields_product
    type Fields 
        complex(wp), dimension(2) :: Field_1, Field_2
    end type Fields
    
    ! Fields construction function override
    interface Fields
        procedure :: Fields_initalization
    end interface Fields
    
    ! opertor override, S_Matrix times Fields
    interface operator (*)
        procedure :: S_Matrix_Fields_product
    end interface operator (*)
    
    contains
    ! constr function of Fields
    pure type(Fields) function  Fields_initalization(Field_1_e, Field_1_h, Field_2_e, Field_2_h)
        implicit none
        complex(wp), intent(in) :: Field_1_e, Field_1_h, Field_2_e, Field_2_h
        Fields_initalization%Field_1(1) = Field_1_e
        Fields_initalization%Field_1(2) = Field_1_h
        Fields_initalization%Field_2(1) = Field_2_e
        Fields_initalization%Field_2(2) = Field_2_h
    end function Fields_initalization
    
    ! Fields and S_Matrix product
    pure type(Fields) function S_Matrix_Fields_product(S_matrix_in, fields_in)
        implicit none
        type(S_Matrix), intent(in) :: S_matrix_in
        type(Fields), intent(in) :: fields_in
        
        S_Matrix_Fields_product%Field_1 = MATMUL(S_matrix_in%alpha_n, fields_in%field_1) + MATMUL(S_matrix_in%gamma_n, fields_in%field_2)
        S_Matrix_Fields_product%Field_2 = MATMUL(S_matrix_in%delta_n, fields_in%field_1) + MATMUL(S_matrix_in%beta_n, fields_in%field_2)
    end function S_Matrix_Fields_product
    
    ! tx and ref coeff of last layer with free space 
    pure function trans_ref_coeff_freespace(S_matrices_in) result(trans_ref_coeff)
        implicit none    
        type(S_Matrix), intent(in), allocatable:: S_matrices_in(:)
        complex(wp), dimension(2,2,2) :: trans_ref_coeff
        type(S_Matrix) :: S_matrices_Cascaded
        integer :: n_S_matrices
        n_S_matrices = SIZEOF(S_matrices_in)/SIZEOF(S_matrices_in(1))
        ! cascade S Matrices
        S_matrices_Cascaded = S_Matrices_Cascade(S_matrices_in, 1, n_S_matrices)
        trans_ref_coeff(1,:,:) = S_matrices_Cascaded%delta_n
        trans_ref_coeff(2,:,:) = S_matrices_Cascaded%alpha_n
    end function trans_ref_coeff_freespace
    
    ! tx and ref coeff of last layer with PEC backed
    pure function trans_ref_coeff_freespace_pec(S_matrices_in, P_n_in) result(trans_ref_coeff)
        implicit none
        type(S_Matrix), intent(in), allocatable :: S_matrices_in(:)
        complex(wp), dimension(2,2), intent(in) :: P_n_in
        complex(wp), dimension(2,2,2) :: trans_ref_coeff
        type(S_Matrix) :: S_matrices_Cascaded
        integer :: n_S_matrices
        n_S_matrices = SIZEOF(S_matrices_in)/SIZEOF(S_matrices_in(1))
        ! cascade S Matrices
        S_matrices_Cascaded = S_Matrices_Cascade(S_matrices_in, 1, n_S_matrices)
        trans_ref_coeff(1,:,:) = MATMUL((unit_matrix + MATMUL(S_matrices_Cascaded%beta_n , P_n_in))**-1 , S_matrices_Cascaded%delta_n)
        trans_ref_coeff(2,:,:) = S_matrices_Cascaded%alpha_n - MATMUL( MATMUL(S_matrices_Cascaded%gamma_n , P_n_in) , trans_ref_coeff(1,:,:))
    end function trans_ref_coeff_freespace_pec
    
    ! calculate Fields on each layer
    pure type(Fields) function cal_Fields(S_matrices_in, fields_in, n)
        type(S_Matrix), intent(in), allocatable :: S_matrices_in(:)
        ! Field_1 for , and Field_2 for
        type(Fields), intent(in) :: fields_in
        integer, intent(in) :: n
        type(S_Matrix) :: S_matrices_Cascaded_1_n, S_matrices_Cascaded_n_N
        integer :: n_S_matrices
        n_S_matrices = SIZEOF(S_matrices_in)/SIZEOF(S_matrices_in(1))
        S_matrices_Cascaded_1_n = S_Matrices_Cascade(S_matrices_in,1,n-1)
        S_matrices_Cascaded_n_N = S_Matrices_Cascade(S_matrices_in,n,n_S_matrices)
        cal_Fields%Field_1 = MATMUL((unit_matrix - MATMUL( S_matrices_Cascaded_1_n%beta_n , S_matrices_Cascaded_n_N%alpha_n))**-1 , MATMUL( S_matrices_Cascaded_1_n%delta_n,fields_in%Field_1) + MATMUL( MATMUL(S_matrices_Cascaded_1_n%beta_n , S_matrices_Cascaded_n_N%gamma_n), fields_in%Field_2))
        cal_Fields%Field_2 = MATMUL((unit_matrix - MATMUL( S_matrices_Cascaded_n_N%alpha_n , S_matrices_Cascaded_1_n%beta_n))**-1, MATMUL( MATMUL( S_matrices_Cascaded_n_N%alpha_n , S_matrices_Cascaded_1_n%delta_n), fields_in%Field_1) + MATMUL(S_matrices_Cascaded_n_N%gamma_n, fields_in%Field_2))
    end function cal_Fields
    
end module Fields_Class