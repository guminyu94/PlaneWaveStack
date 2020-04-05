!****************************************************************************
!
!   Module: Swapper
!
!   PURPOSE: Swap parameters and output to file
!
!****************************************************************************
Module Swapper
    use Sim_parameters, only : wp, n_layers
    use Layer_Class
    use S_Matrix_Class
    use Fields_Class
    implicit none
    type(Layer), allocatable :: layers(:)
    type(S_Matrix), allocatable :: S_Matrices(:)
    
    abstract interface
        subroutine fun_temp(f, larray)
            use Sim_parameters, only : wp
            use Layer_Class
            real(wp), intent(in) :: f
            type(Layer), allocatable, intent(inout) :: larray(:)
        end subroutine fun_temp
    end interface

    contains
    subroutine freq_swap(fun,freq_start,freq_end,n_points,file_name)
        real(wp), intent(in) :: freq_start, freq_end
        integer, intent(in) :: n_points
        procedure (fun_temp), pointer :: fun
        real(wp) :: freq_step
        complex(wp), dimension(2,2,2) :: tx_ref
        character(len=*), intent(in) :: file_name
        integer :: i, j
        
        freq_step = (freq_end - freq_start) / real(n_points-1,wp)
        
        ! output data into a file
        open(1, file = file_name, status = 'OLD')
        do j = 1, n_points
            call fun(real(j-1,wp)*freq_step + freq_start,layers)
            
            ! allocate the S_Matrix 
            if (j .EQ. 1) then
                allocate(S_Matrices(n_layers-1))
            end if
            
            ! compute S matrix
            do i = 1, n_layers-1
                ! assemble S Matrix from layer obj
                S_Matrices(i) = S_Matrix(layers(i),layers(i+1))
            end do
            
            ! obtain reflection coeff
            tx_ref = trans_ref_coeff_freespace(S_Matrices)
            ! write to file
            write(1,*) real(j-1,wp)*freq_step + freq_start, ABS(tx_ref(1,1,1))**2.0_wp
            
            
            print *, "|tx_coeff_ee|^2"
            print *, ABS(tx_ref(1,1,1))**2.0_wp
            print *, "|tx_coeff_eh|^2"
            print *, ABS(tx_ref(1,1,2))**2.0_wp
            print *, "|ref_coeff_ee|^2"
            print *, ABS(tx_ref(2,1,1))**2.0_wp
            print *, "|ref_coeff_hh|^2"
            print *, ABS(tx_ref(2,2,2))**2.0_wp
        end do
    
    close(1) 
    end subroutine freq_swap

end module Swapper