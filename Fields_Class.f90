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
    public :: Fields
    type Fields 
        complex(wp), dimension(2) :: Field_1, Field_2
    end type Fields
    
    ! Fields construction function override
    interface Fields
        procedure :: Fields_initalization_1, Fields_initalization_1
    end interface Fields
    
    ! opertor override, S_Matrix times Fields
    interface operator (*)
        procedure :: Matrix_Fields_product
    end interface operator (*)
    
    contains
    ! constr function of Fields
    function  Fields_initalization_1(Field_1_e, Field_1_h, Field_2_e, Field_2_h) result (field_out)
        implicit none
        complex(wp), intent(in) :: Field_1_e, Field_1_h, Field_2_e, Field_2_h
        type(Fields):: field_out
        field_out%Field_1(1) = Field_1_e
        field_out%Field_1(2) = Field_1_h
        field_out%Field_2(1) = Field_2_e
        field_out%Field_2(2) = Field_2_h
    end function Fields_initalization_1
    
    pure type(Fields) function  Fields_initalization_2(Field_1, Field_2)
        implicit none
        complex(wp), dimension(2), intent(in) :: Field_1, Field_2
        Fields_initalization_2%Field_1(1) = Field_1(1)
        Fields_initalization_2%Field_1(2) = Field_1(2)
        Fields_initalization_2%Field_2(1) = Field_2(1)
        Fields_initalization_2%Field_2(2) = Field_2(2)
    end function Fields_initalization_2
    
    ! Fields and S_Matrix product
    pure type(Fields) function Matrix_Fields_product(matrix_in, fields_in)
        implicit none
        type(S_Matrix), intent(in) :: matrix_in
        type(Fields), intent(in) :: fields_in
        
        Matrix_Fields_product%Field_1 = MATMUL(matrix_in%alpha_n, fields_in%field_1) + MATMUL(matrix_in%gamma_n, fields_in%field_2)
        Matrix_Fields_product%Field_2 = MATMUL(matrix_in%delta_n, fields_in%field_1) + MATMUL(matrix_in%beta_n, fields_in%field_2)
    end function Matrix_Fields_product
    
    
    
end module Fields_Class