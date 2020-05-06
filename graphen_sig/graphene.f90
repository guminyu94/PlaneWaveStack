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
    
    real, parameter :: muc= 1.2  !1.0 ! 0.2 ! eV

    real, parameter :: tk = 300

    real, parameter :: vf= 1.0e6

    real, parameter :: tau=0.5e-12 !0.5e-12   !3.0e-12 !  s 

    real, parameter :: sigmin=e/(4.0*hbar)

    real, parameter :: kbt=kB*tK
    
    real, parameter :: b0 = 0
    
end module graphene
