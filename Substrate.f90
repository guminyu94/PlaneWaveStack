module Substrate
    use Sim_parameters, only : wp
    implicit none
    
    complex(wp), parameter :: sic_h6_eps_eps_t = (6.70_wp, 0.0_wp)
    complex(wp), parameter :: sic_h6_eps_nu_e = (0.9731343_wp, 0.0_wp)    
    complex(wp), parameter :: prism_e = (14.7456_wp, 0.0_wp)
    complex(wp), parameter :: polymethylpentene_e = (2.1025_wp, 0.0_wp)
end module Substrate