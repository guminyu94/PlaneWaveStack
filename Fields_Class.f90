!****************************************************************************
!
!   class: Fields
!
!   PURPOSE:  contains 4 fields components, forward, backward, e and h
!
!****************************************************************************
    module Fields_Class
    implicit none
    use S_Matrix_Class
    private
    public :: Fields
    type Fields 
        complex, dimension(2) :: Field_1, Field_2
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
    type(Fields) function  Fields_initalization(Field_1_e, Field_1_h, Field_2_e, Field_2_h)
        implicit none
        complex, intent(in) :: Field_1_e, Field_1_h, Field_2_e, Field_2_h
        Fields_initalization%Field_1(1) = Field_1_e
        Fields_initalization%Field_1(2) = Field_1_h
        Fields_initalization%Field_2(1) = Field_2_e
        Fields_initalization%Field_2(2) = Field_2_h
    end function Fields_initalization
    
    ! Fields and S_Matrix product
    type(Fields) function S_Matrix_Fields_product(s_matrix_in, fields_in)
        implicit none
        type(S_Matrix), intent(in) :: s_matrix_in
        type(Fields), intent(in) :: fields_in
        S_Matrix_Fields_product%Field_1 = MATMUL(s_matrix_in%alpha_n, fields_in%field_1) + MATMUL(s_matrix_in%gamma_n, fields_in%field_2)
        S_Matrix_Fields_product%Field_2 = MATMUL(s_matrix_in%delta_n, fields_in%field_1) + MATMUL(s_matrix_in%beta_n, fields_in%field_2)
    end function S_Matrix_Fields_product
    
end module Fields_Class