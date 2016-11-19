module gold

real(8), public:: freqA
!parameters (Drude model)
!vector format: (wp, damp)
!Mostepanenko
real(8), dimension(2), parameter, public:: parMost=(/13.671e15, 0.5317e14/)
!approximation data of the imaginary part (Palik)
real(8), dimension(2), parameter, public:: parAproxIm=(/10.06e15, 1.36e14/)
!approximation data of the real part (Palik)
real(8), dimension(2), parameter, public:: parAproxRe=(/12.06e15, 1.219e14/)

!oscillator model, Mostepanenko parameters
real(8), dimension(6), parameter:: gAu=(/7.091, 41.46, 2.7, 154.7, 44.55, 309.6/)*(1.519e15)**2
real(8), dimension(6), parameter:: gammaAu=(/0.75, 1.85, 1.0, 7.0, 6.0, 9.0/)*1.519e15
real(8), dimension(6), parameter:: wAu=(/3.05, 4.15, 5.4, 8.5, 13.5, 21.5/)*1.519e15


contains

function Drude(x, param)
!Drude-model permettivity function eps(iw)
!param - vector of parameters (wp, damp)
real(8):: x, Drude, param(2)

	Drude=1._8+param(1)**2/(x**2+param(2)*x);

end function

!----------------------------------------------------

function epsMar(x)
!Marachevsky formula fo eps(iw)
real(8):: x,epsMar,num,denum
real(8), dimension(2,4):: paramMar

!Marachevsky model parameters (wl1,wl2//gl1,gl2//gt1,gt2//wt2,0)
paramMar(1,1)=exp(-0.96)*1e16
paramMar(2,1)=exp(0.2866)*1e16
paramMar(1,2)=exp(-2.536)*1e16
paramMar(2,2)=exp(1.255)*1e16
paramMar(1,3)=exp(-4.7922)*1e16
paramMar(2,3)=exp(-0.957)*1e16
paramMar(1,4)=exp(-0.8359)*1e16
paramMar(2,4)=0.0_8

num=(paramMar(1,1)**2+x**2+paramMar(1,2)*x)*(paramMar(2,1)**2+x**2+paramMar(2,2)*x)
denum=(x**2+paramMar(1,3)*x)*(paramMar(1,4)**2+x**2+paramMar(2,3)*x)

epsMar=num/denum

end function

!---------------------------------------------------

function oscil_to_int(x)

real(8):: oscil,x,oscil_to_int

oscil=0.0_8

do i=1,6
oscil=oscil+gAu(i)*gammaAu(i)*x/((wAu(i)**2-x**2)**2+(gammaAu(i)*x)**2)
end do

oscil_to_int=oscil+parAproxRe(1)**2*parAproxRe(2)/((x**2+parAproxRe(2)**2)*(x**2+freqA**2))

end function

!--------------------------------------------------

function drude_to_int(x)

real(8):: x, drude_to_int

drude_to_int=parAproxRe(1)**2*parAproxRe(2)/((x**2+parAproxRe(2)**2)*(x**2+freqA**2))

end function

!-------------------------------------------------
!Generalized plasma-like dielectric permettivity function

function Gen_Plasma(x)

real(8):: x, Gen_Plasma

Gen_Plasma=1._8+(parMost(1)/x)**2

do i=1,6
Gen_Plasma=Gen_Plasma+gAu(i)/(wAu(i)**2+x**2-gammaAu(i)*x)
end do

end function
!-------------------------------------------------

!Want to calculate zero summand

! Drude, Marachevsky, Kramers-Kr
function zero_sum1(a,x)

real(8):: zero_sum1, x
real(8):: rTm, rTe, a

!silicon parameters
!theory - Drude-Loretz model
!parameters from the article

rTe=0._8
rTm=1._8

zero_sum1=x*(log(1._8-rTm**2*exp(-2.*x))+log(1._8-rTe**2*exp(-2.*x)))

end function

!----------------------------------------------------------------------

!Generalized plama-like model

function zero_sum2(a,x)

real(8):: zero_sum2, x
real(8):: rTm, rTe, a, log1, log2
real(8):: wp,c
integer:: k

c=299792458._8
wp=parMost(1)

rTm=1._8

rTe=((x/a)-sqrt((x/a)**2+(wp/c)**2))/((x/a)+sqrt((x/a)**2+(wp/c)**2))

log1 = log(1._8-rTe**2*exp(-2.*x))
log2 = log(1._8-rTm**2*exp(-2.*x))

zero_sum2=x*(log1+log2)

end function
!---------------------------------------------------------------------

function zero_sum_int_Dr(x)
use dopcasimir

real(8):: zero_sum_int_Dr, x

if (x.eq.(0._8)) then
    zero_sum_int_Dr=0._8
    else
zero_sum_int_Dr=zero_sum1(dist,x)
endif

end function
!---------------------------------------------------------------------

function zero_sum_int1_Dr(x)
use dopcasimir

real(8):: zero_sum_int1_Dr, x

if (x.eq.(0._8)) then

    zero_sum_int1_Dr=0._8
    else

    zero_sum_int1_Dr=(1._8/x)**2*zero_sum1(dist,1._8/x)
endif

end function
!---------------------------------------------------------------------


function zero_sum_int_Pl(x)
use dopcasimir

real(8):: zero_sum_int_Pl, x

if (x.eq.(0._8)) then
    zero_sum_int_Pl=0._8
    else
zero_sum_int_Pl=zero_sum2(dist,x)
endif

end function
!---------------------------------------------------------------------

function zero_sum_int1_Pl(x)
use dopcasimir

real(8):: zero_sum_int1_Pl, x

if (x.eq.(0._8)) then

    zero_sum_int1_Pl=0._8
    else

    zero_sum_int1_Pl=(1._8/x)**2*zero_sum2(dist,1._8/x)
endif

end function
!---------------------------------------------------------------------


end module



