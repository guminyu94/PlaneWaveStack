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
    complex(wp), parameter :: air = (1.0_wp,0.0_wp)
    
    contains
    pure function pc_quadwl(eps,freq) result(thickness)
    ! compute the thickness of dieletric layer to construct 1/4 wavelength photonic crystal based on eps input and freq
        use constants, only : c0
        real(wp) :: thickness
        complex(wp), intent(in) :: eps
        real(wp), intent(in) :: freq
        thickness = c0 / freq / 4.0_wp / (real(eps)**0.5_wp)
    end function pc_quadwl
    
    subroutine pc_thickness_config(dieletric_array,thickness_array,freq)
    ! generate the photonic crystal thickness array based on the dieletrics array input
        complex(wp), allocatable, intent(in) :: dieletric_array(:)
        real(wp), allocatable, intent(inout) :: thickness_array(:)
        real(wp), intent(in) :: freq
        integer :: i
        
        if (.NOT. allocated(thickness_array)) then
            allocate(thickness_array(size(dieletric_array)))
        end if
        
        do i = 1,size(dieletric_array)
            thickness_array(i) = pc_quadwl(dieletric_array(i),freq)
        end do
    
        
    end subroutine pc_thickness_config
    
end module Substrate