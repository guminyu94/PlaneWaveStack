!****************************************************************************
!
!   Module: graphene
!
!   PURPOSE: adjustable parameters for graphene
!
!****************************************************************************
    
    module graphene
    
    use Constants, only : kb, hbar, e
    
    implicit none
    
    real :: muc = 1.2   !1.0 ! 0.2 ! eV

    real, parameter :: tk = 10 ! 300K

    real, parameter :: vf = 1.0e6

    real :: tau = 0.5*1/(20836617636.1327*15) !0.5e-12   !3.0e-12 !  

    real, parameter :: sigmin = e / (4.0 * hbar)

    real, parameter :: kbt = kB * tK
    
    real :: b0 = 1.0
    
end module graphene
