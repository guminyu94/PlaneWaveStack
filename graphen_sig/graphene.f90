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
    
    real :: muc = 0.2   !1.0 ! 0.2 ! eV

    real, parameter :: tk = 10 ! 300K

    real, parameter :: vf = 1.0e6

    real :: tau = 0.5e-12 !0.5e-12   !3.0e-12 !  

    real, parameter :: sigmin = e / (4.0 * hbar)

    real :: kbt = kB * tK
    
    real :: b0 = 0.5
    
end module graphene
