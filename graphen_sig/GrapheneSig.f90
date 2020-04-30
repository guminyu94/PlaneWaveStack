module GrapheneSig
use Sim_parameters, only : wp
implicit none
private
public :: sigmas

! graphene computation
complex, public :: sig_d, sig_h
complex(wp), public :: sigxx, sigyy, sigxy, sigyx
integer, public :: n_d, n_h

contains

subroutine sigmas(freq,sigd,sigh,nd,nh)

use Constants, only : hbar, pi, j, zero

use graphene, only : b0, muc, vf, tau, kbt, sigmin

use fermipass, only : omgel, mucel, kbtel

real, intent(in) :: freq

complex, intent(out) :: sigd, sigh

integer, intent(out) :: nd, nh

real, parameter :: klein=1e-2, eps=1e-8

integer, parameter :: nmax=500000

complex :: s, omh

real :: omega, b, el2, el

integer :: nbeg

omega = 2.0*pi*freq

omh = hbar*(omega-j/tau)

if (abs(b0)<klein)then
    
    sigh=zero

    b=log(2.0*(1.0+cosh(muc/kbt)))
    
    sigd=4.0/(j*pi)*(kbt/omh)*b
    
    sigd=sigd*sigmin
    
else

    el2=2.0*vf**2*abs(b0)*hbar
    
    el=sqrt(el2)
    
    kbtel=kbt/el
    
    mucel=muc/el
    
    omgel=omh/el
    
    nbeg=int((muc/el)**2)

    call sumsig(fund,eps,nbeg,nmax,s,nd)

    sigd=-2.0/(j*pi)*omgel*s
    
    sigd=sigd*sigmin

    call sumsig(funh,eps,nbeg,nmax,s,nh)
    
    sigh=2.0/pi*sign(1.0,b0)*s
    
    sigh=sigh*sigmin

end if

end subroutine sigmas


subroutine sumsig(fun,eps,nbeg,nmax,s,ns)

use Constants, only : zero

interface

    function fun(x,y) result(val)

    real, intent(in) :: x, y
    
    complex :: val

    end function fun
    
end interface

real, intent(in) :: eps

integer, intent(in) :: nbeg, nmax

complex, intent(out) :: s

integer, intent(out) :: ns

real :: em1, em2

complex :: q, sold

integer :: n, n1, n2

em1=sqrt(real(nbeg))

em2=sqrt(real(nbeg+1))
        
sold=fun(em1,em2); ns=1

n1=nbeg

n2=nbeg

do n=1,nmax
    
    n1=n1+1

    em1=sqrt(real(n1))

    em2=sqrt(real(n1+1))
        
    q=fun(em1,em2); ns=ns+1
    
    s=sold+q
    
    n2=n2-1
    
    if(n2>=0)then
        
        em1=sqrt(real(n2))

        em2=sqrt(real(n2+1))
        
        q=fun(em1,em2); ns=ns+1
    
        s=s+q

    end if

    if(abs(s-sold) < eps*abs(s))exit
    
    sold=s

end do

end subroutine sumsig


function nf(x) result(val)

use fermipass, only : mucel, kbtel

real, intent(in) :: x

real :: val, exx

exx=exp((x-mucel)/kbtel)

val=1.0/(1.0+exx)

end function nf


function f1(x,y) result(val)

real, intent(in) :: x, y

real :: val

val=nf(x)-nf(y)

end function f1


function fund(x,y) result(val)

use fermipass, only : omgel

real, intent(in) :: x, y

real :: num1, num2

complex :: val, den1, den2

num1=f1(x,y)-f1(-x,-y)

num2=f1(-x,y)-f1(x,-y)

den1=(y-x)*((y-x)**2-omgel**2)

den2=(y+x)*((y+x)**2-omgel**2)

val=num1/den1+num2/den2

end function fund


function funh(x,y) result(val)

use fermipass, only : omgel

real, intent(in) :: x, y

real :: num

complex :: val, den1, den2

num=f1(x,y)+f1(-x,-y)

den1=(y-x)**2-omgel**2

den2=(y+x)**2-omgel**2

val=num*(1.0/den1+1.0/den2)

end function funh


subroutine SigTensor(xi,sigd,sigh,sigmat)

real, intent(in) :: xi

complex, intent(in) :: sigd, sigh

complex, dimension(2,2), intent(out) :: sigmat

complex :: sigxx, sigxy, sigyx, sigyy, c2, s2, cs

sigxx=sigd

sigyy=sigd

sigyx=sigh

sigxy=-sigh

c2=cos(xi)**2

s2=sin(xi)**2

cs=cos(xi)*sin(xi)

sigmat(1,1)=sigxx*c2+sigyy*s2+(sigxy+sigyx)*cs

sigmat(1,2)=sigxy*c2-sigyx*s2-(sigxx-sigyy)*cs

sigmat(2,1)=sigyx*c2-sigxy*s2-(sigxx-sigyy)*cs

sigmat(2,2)=sigyy*c2+sigxx*s2-(sigxy+sigyx)*cs

end subroutine SigTensor


end module GrapheneSig


