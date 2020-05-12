module Substrate
    use Sim_parameters, only : wp
    implicit none
    
    complex(wp), parameter :: sic_h6_e = (6.70_wp, 0.0_wp)
    complex(wp), parameter :: sic_h6_nu = (0.9731343_wp, 0.0_wp)
    complex(wp), parameter :: sic_e = (6.5025_wp, 0.0_wp)
    complex(wp), parameter :: prism_e = (14.7456_wp, 0.0_wp)
    complex(wp), parameter :: polymethylpentene_e = (2.1025_wp, 0.0_wp)
    complex(wp), parameter :: sio2_e = (3.61_wp,0.0_wp)
    complex(wp), parameter :: si_e = (11.6964_wp,0.0_wp)
    
    
end module Substrate