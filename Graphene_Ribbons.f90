!****************************************************************************
!
!   Module: Graphene_Ribbons
!
!   PURPOSE: Compute effective sigma of graphene ribbons, based on electrostatic approach
!
!****************************************************************************
    Module Graphene_Ribbons
    contains
    function gr_sigma(freq,W,L,G,eps_sub) result(sigma_eff_mat)
        use GrapheneSig
        use Sim_parameters, only : eps_0, wp, PI
        implicit none
        ! W is the width of graphen ribbons, L is the periodic distance between two ribbons, and G is the gap between two ribbons
        real(wp), intent(in) :: W, L, G, freq
        complex(wp), intent(in) :: eps_sub
        complex(wp), dimension(2,2) :: sigma_eff_mat
        complex(wp) :: sigma_c, sigma_d, sigma_h
        sigma_c = (0.0_wp,-2.0_wp) * 2.0_wp * PI * freq * eps_sub * eps_0 * (L/PI) * log(1.0_wp / sin(PI * G / 2 / L))
        ! compute the graphene sigma
        call sigmas(real(freq),sig_d,sig_h,n_d,n_h)
        sigma_d = cmplx(sig_d)
        sigma_h = cmplx(sig_h)
        sigma_eff_mat(1,1) = W * sigma_d * sigma_c / (L * sigma_c + W * sigma_d)
        sigma_eff_mat(1,2) = sigma_eff_mat(1,1) * W / L * -sigma_h / sigma_d
        sigma_eff_mat(2,1) = -sigma_eff_mat(1,2)
        sigma_eff_mat(2,2) = W / L * sigma_d - W / L * -sigma_h * sigma_h / sigma_d / sigma_c + sigma_eff_mat(1,2) * sigma_eff_mat(2,1) / sigma_eff_mat(1,1)
    end function gr_sigma
end module Graphene_Ribbons