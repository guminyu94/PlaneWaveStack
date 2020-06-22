module Substrate
    use Sim_parameters, only : wp
    implicit none
    
    complex(wp), parameter :: sic_h6_e = (6.70_wp, 0.0_wp)
    complex(wp), parameter :: sic_h6_nu = (0.9731343_wp, 0.0_wp)
    complex(wp), parameter :: sic_e = (6.5025_wp, 0.0_wp)
    complex(wp), parameter :: prism_e = (14.7456_wp, 0.0_wp)
    complex(wp), parameter :: polymethylpentene_e = (2.1025_wp, 0.0_wp)
    complex(wp), parameter :: polymethylpentene_e_1t = (3.6100_wp, 0.0_wp)
    complex(wp), parameter :: sio2_e = (3.61_wp,0.0_wp)
    complex(wp), parameter :: sio2_e_1t = (2.1316_wp, 0.0_wp)
    complex(wp), parameter :: si_e = (11.6964_wp,0.0_wp)
    complex(wp), parameter :: air_e = (1.0_wp,0.0_wp)
    complex(wp), parameter :: si_20t_e = (2.25_wp,0.0_wp)
    contains
    pure function pc_wl(eps,freq,factor) result(thickness)
    ! compute the thickness of dieletric layer to construct 1/4 wavelength photonic crystal based on eps input and freq
        use constants, only : c0
        real(wp) :: thickness
        complex(wp), intent(in) :: eps
        real(wp), intent(in) :: freq
        real(wp), intent(in) :: factor
        thickness = c0 / freq / factor / (real(eps)**0.5_wp)
    end function pc_wl
    
    subroutine pc_thickness_config(dieletric_array,thickness_array,freq,factor)
    ! generate the photonic crystal thickness array based on the dieletrics array input
        complex(wp), allocatable, intent(in) :: dieletric_array(:)
        real(wp), allocatable, intent(inout) :: thickness_array(:)
        real(wp), intent(in) :: freq
        integer :: i
        real(wp), intent(in) :: factor
        
        if (.NOT. allocated(thickness_array)) then
            allocate(thickness_array(size(dieletric_array)))
        end if
        
        do i = 1,size(dieletric_array)
            thickness_array(i) = pc_wl(dieletric_array(i),freq,factor)
        end do
    end subroutine pc_thickness_config
    
    ! peak the tx ref of part of the layers stucture, pass the whole layers, given start and end index
    ! direction 0 stand for left to right, 1 for right to left(higher index to lower index)
    function peak_txref_coeff(layers,start_index,end_index,is_pec_back,direction) result(trans_ref)
        use Layer_Class
        use Swapper, only : trans_ref_coeff_pec, trans_ref_coeff_freespace
        use S_Matrix_Class
        implicit none
        integer, intent(in) :: start_index, end_index
        type(Layer), intent(in), allocatable ::  layers(:)
        integer, intent(in) :: is_pec_back, direction
        type(S_Matrix), allocatable :: S_matrices(:)
        type(Layer), allocatable ::  layers_temp(:)
        complex(wp), dimension(2,2,2) :: trans_ref
        integer :: i

        if (direction .eq. 0) then
            layers_temp = layers(start_index:end_index)
        else if (direction .eq. 1) then
            ! flip the order
            allocate(layers_temp(end_index-start_index+1))
            do i = 1, (end_index-start_index+1)
                layers_temp(i) = layers(end_index-i+1)
            end do
        end if
        
        call compute_S_matrix(layers_temp,S_matrices)
        
        if (is_pec_back .EQ. 1) then
            trans_ref = trans_ref_coeff_pec(S_matrices,layers_temp(size(layers_temp))%P_n)
        else
            trans_ref = trans_ref_coeff_freespace(S_matrices)
        end if
        
        deallocate(S_matrices)
        deallocate(layers_temp)
    end function peak_txref_coeff
    
end module Substrate