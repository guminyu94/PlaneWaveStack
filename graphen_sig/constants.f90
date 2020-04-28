module Constants
use Math, only : PI
implicit none

real, parameter :: kb=8.6173303e-5 ! eV/K

real, parameter :: hbar=6.582119514e-16 !eVs

real, parameter :: twopi=6.283185307179586476925286766559005768394

real, parameter :: e=1.602176634e-19

real, parameter :: c0=299792458

real, parameter :: mu0=4*pi*1.0e-7

real, parameter :: eta0=c0*mu0

complex, parameter :: j=(0.0,1.0)

complex, parameter :: one=(1.0,0.0)

complex, parameter :: zero=(0.0,0.0)

end module Constants