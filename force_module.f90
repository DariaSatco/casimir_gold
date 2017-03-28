module force_module

use dopcasimir
use gold

contains

function Lifshitz(a,x,y)
    ! function needed to integrate in formula for Casimir force
    ! x=dist*r, y=w*dist/c change the variables to dimensioness
    ! the constants in front of the integral: rdr --> 1/dist**2 xdx;  dw --> c/dist
    ! w in eV
    ! c in eV*m (defined as cEv in module parameters)

real(4):: Lifshitz
real(4):: x
real(8):: a,y
real(8):: q,rTMn(2),rTEn(2)

rTMn=rTM(x/a,y)
!calculate the TM Frenel coefficient for the concrete dist as a function of x,y
rTEn=rTE(x/a,y)
!calculate the TE Frenel coefficient for the concrete dist as a function of x,y

q=sqrt(x**2+(a*y/cEv)**2)

Lifshitz=q*x*( ( ( rTMn(1) * rTMn(2) )**(-1)*exp( 2*q ) - 1.0 )**(-1) + &
( ( rTEn(1) * rTEn(2) )**(-1)*exp( 2*q ) - 1.0 )**(-1) )


end function

!------------------------------------------------------------------------

function lif_to_int(arg)
    !integrable Lifshitz function for nonzero temperatures - function of one argument - x <--> k-frequency
    !interval [0,inf]
    !w in eV

real(4):: arg, lif_to_int

    lif_to_int = Lifshitz( dist, arg, w )

end function

!-----------------------------------------------------------------------
function zero_sum_force_Dr(x)
!does not depend on distance if we have done change x --> dist*x
!Drude, Marachevsky, Kramers-Kr

real(4):: zero_sum_force_Dr, x
real(4):: rTm, rTe

rTe=0.0
rTm=1.0

zero_sum_force_Dr = x**2*( ( rTm**(-2)*exp( 2*x ) - 1.0 )**(-1) + &
( rTe**(-2)*exp( 2*x ) - 1.0 )**(-1) )

end function

!----------------------------------------------------------------------
function zero_sum_force_Pl(a,x)
!Generalized plasma
use gold

real(4):: zero_sum_force_Pl, x
real(8):: rTm, rTe, a
real(8):: wp,cEv

cEv=19.746e-8 ! eV*s
wp=parMost(1)

rTm=1._8

rTe=((x/a)-sqrt((x/a)**2+(wp/cEv)**2))/((x/a)+sqrt((x/a)**2+(wp/cEv)**2))

zero_sum_force_Pl = x**2*( ( rTM**(-2)*exp( 2*x ) - 1.0 )**(-1) + &
(rTE**(-2)*exp( 2*x ) - 1.0 )**(-1) )

end function

!-------------------------------------------------------------------------

function zero_sum_int_force_Pl(x)
!Plasma zero summand function to integrate from 0 to 1
use dopcasimir

real(4):: zero_sum_int_force_Pl, x

if (x.eq.(0.0)) then
    zero_sum_int_force_Pl = 0.0
    else
    zero_sum_int_force_Pl = zero_sum_force_Pl(dist,x)
endif

end function

!-------------------------------------------------------------------------


end module force_module
