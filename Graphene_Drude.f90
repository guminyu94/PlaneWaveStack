!****************************************************************************
!
!   Module: Graphene_Drude
!
!   PURPOSE: Drude model of graphene sigma
!
!****************************************************************************
    Module Graphene_Drude
    implicit none
    contains
    function graphene_Drude_sig(freq) result(sigma_d)
        use Sim_parameters, only : wp
        ! hbar in eV
        use constants, only : hbar, PI, e
        ! kbt, muc in eV
        use graphene, only : kbt, tau, muc
        
        complex(wp) :: sigma_d
        real(wp) ::  omega
        real(wp), intent(in) :: freq
        omega = 2 * PI * freq
        sigma_d = (0.0_wp,-1.0_wp) * (kbt * e)  / ( PI * (hBar ** 2.0_wp) * (omega + (0.0_wp,-1.0_wp)/tau ) ) * ( muc / kbt  + 2.0_wp * log( exp( - muc / kbt ) + 1.0_wp ))
    end function graphene_Drude_sig
    
end module Graphene_Drude