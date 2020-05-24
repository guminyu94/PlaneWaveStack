!****************************************************************************
!
!   Module: Math
!
!   PURPOSE: Math functions and operators
!
!****************************************************************************
Module Math
    use Sim_parameters, only : wp, PI
    implicit none
    public
    complex(wp), dimension(2,2), parameter :: unit_matrix=(/(1.0_wp,0.0_wp),(0.0_wp,0.0_wp),(0.0_wp,0.0_wp),(1.0_wp,0.0_wp)/)

    interface operator(*)
        procedure :: matrix_vector_product
    end interface operator(*)
    
    interface operator(**)
        procedure :: matrix_inverse
    end interface operator(**)
    
    contains
    
    
    pure function matinv2(A) result(B)
        !! Performs a direct calculation of the inverse of a 2×2 matrix.
        complex(wp), intent(in) :: A(2,2)   !! Matrix
        complex(wp)             :: B(2,2)   !! Inverse matrix
        complex(wp)             :: detinv

        ! Calculate the inverse determinant of the matrix
        detinv = 1/(A(1,1)*A(2,2) - A(1,2) * A(2,1))

        ! Calculate the inverse of the matrix
        B(1,1) = +detinv * A(2,2)
        B(2,1) = -detinv * A(2,1)
        B(1,2) = -detinv * A(1,2)
        B(2,2) = +detinv * A(1,1)
    end function

    pure function matinv3(A) result(B)
        !! Performs a direct calculation of the inverse of a 3×3 matrix.
        complex(wp), intent(in) :: A(3,3)   !! Matrix
        complex(wp)             :: B(3,3)   !! Inverse matrix
        complex(wp)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
              - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
              + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
  end function

  pure function matinv4(A) result(B)
    !! Performs a direct calculation of the inverse of a 4×4 matrix.
    complex(wp), intent(in) :: A(4,4)   !! Matrix
    complex(wp)             :: B(4,4)   !! Inverse matrix
    complex(wp)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = &
      1/(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
       - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
       + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
       - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

    ! Calculate the inverse of the matrix
    B(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
    B(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
    B(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
    B(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
    B(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
    B(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
    B(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
    B(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
    B(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
    B(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
    B(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
    B(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
  end function
  
  ! matrix inversion 
   pure function matrix_inverse(matrix_in, power) result(inversed_matrix)
        implicit none 
        complex(wp), dimension(2,2), intent(in) :: matrix_in
        complex(wp), dimension(2,2)  :: inversed_matrix
        integer, intent(in) :: power
            ! check if it's a diagonal matrix
        if (power .EQ. -1) then
            if ((matrix_in(1,2) .EQ. (0.0_wp,0.0_wp)) .AND. (matrix_in(2,1) .EQ. (0.0_wp,0.0_wp))) then
                inversed_matrix(1,1) = matrix_in(1,1)**-1.0_wp
                inversed_matrix(1,2) = 0.0_wp
                inversed_matrix(2,1) = 0.0_wp
                inversed_matrix(2,2) = matrix_in(2,2)**-1.0_wp
            else
                inversed_matrix = matinv2(matrix_in)
            end if
        end if
   end function matrix_inverse
   
   ! matrix_vector_product
   pure function matrix_vector_product(matrix_l_in, vector_r_in) result(vector_out)
        implicit none 
        complex(wp), dimension(2,2), intent(in) :: matrix_l_in
        complex(wp), dimension(2), intent(in) :: vector_r_in
        complex(wp), dimension(2)  :: vector_out
        vector_out = MATMUL(matrix_l_in,vector_r_in)
   end function matrix_vector_product
  
   ! compute the principal axis angle of ellipse
   pure function ellipse_angle(Ex,Ey) result(theta)
        complex(wp), intent(in) :: Ex, Ey
        complex(wp) :: Eref
        real(wp) :: theta, phi
        
        theta = atan2(2.0_wp * real(Ex * CONJG(Ey)) ,(abs(Ex)**2.0_wp - abs(Ey)**2.0_wp)) / 2.0_wp
        
   end function ellipse_angle
   
   subroutine phase_unwrap_1d(xw)
        real, intent(inout), allocatable :: xw(:)
        real, allocatable :: xu(:)
        real :: difference
        integer :: i, n_x
        n_x = size(xw)
        allocate(xu(n_x))
        xu = xw
        ! 1D Phase Unwrap
        do i = 2, n_x
            difference = xw(i)-xw(i-1);
            if (difference > PI) then
                xu(i:n_x) = xu(i:n_x) - 2.0 * PI;
            else if (difference < -PI) then
                xu(i:n_x) = xu(i:n_x) + 2.0 * PI;
            end if
        end do
        xw = xu
        deallocate(xu)
   end subroutine phase_unwrap_1d
   
   
end Module Math