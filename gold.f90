module gold

!contains different models for gold dielectric permittivity function
! and calculation of zero summands for Lifshitz formula

interface
    SUBROUTINE WOFZ (XI, YI, U, V, FLAG)
    real(8):: xi, yi, u,v
    logical::flag
    end subroutine
end interface

real(8), public:: freqA
!parameters
!vector format: (wp, damp)
!Mostepanenko
real(8), dimension(2), parameter:: parMost=(/9.0, 0.035/) !eV
!Rakic Lorentz-Drude
real(8), dimension(2), parameter:: parRakic_LD=(/7.87, 0.053/) !eV
!Rakic Brendel-Bormann
real(8), dimension(2), parameter:: parRakic_BB=(/7.92, 0.05/) !eV
!approximation data of the imaginary part (Palik)
real(8), dimension(2), parameter:: parAproxIm=(/9.0, 0.024/) !eV
!approximation data of the real part (Palik)
real(8), dimension(2), parameter:: parAproxRe=(/7.93, 0.08/) !eV

!oscillator model, Mostepanenko parameters
real(8), dimension(6), parameter:: gAu=(/7.091, 41.46, 2.7, 154.7, 44.55, 309.6/) !eV
real(8), dimension(6), parameter:: gammaAu=(/0.75, 1.85, 1.0, 7.0, 6.0, 9.0/) !eV
real(8), dimension(6), parameter:: wAu=(/3.05, 4.15, 5.4, 8.5, 13.5, 21.5/) !eV

!Lorentz-Drude model, Rakic parameters
real(8), dimension(6), parameter:: fR=(/0.024, 0.01, 0.071, 0.601, 4.384, 0.0/) !eV
! gR=fR*wpR**2
real(8), dimension(6), parameter:: gammaR=(/0.241, 0.345, 0.87, 2.494, 2.214, 0.0/) !eV
real(8), dimension(6), parameter:: wR=(/0.415, 0.83, 2.969, 4.304, 13.32, 0.0/) !eV

!Brendel-Bormann model, Rakic parameters
real(8), dimension(6), parameter:: fRb=(/0.054, 0.05, 0.312, 0.719, 1.648, 0.0/)
! gRb=fRb*wpR**2
real(8), dimension(6), parameter:: gammaRb=(/0.074, 0.035, 0.083, 0.125, 0.179, 0.0/)
real(8), dimension(6), parameter:: wRb=(/0.218, 2.885, 4.069, 6.137, 27.97, 0.0/)
real(8), dimension(6), parameter:: sigmaRb=(/0.742, 0.349, 0.830, 1.246, 1.795, 1.0/)

!Racik plasma frequency
real(8), parameter:: wpR=9.03

!Gauss model
!imaginary part
real(8), dimension(6), parameter:: amplitude_im=(/ 6.139, 1.567, 1.524, 0.8617, 1.997, 0.4712 /)
real(8), dimension(6), parameter:: w_band_im=(/ 3.73, 8.153, 5.5, 7.54, 11.3, 20.85 /)
real(8), dimension(6), parameter:: sigma_im=(/ 1.177, 5.202, 0.9631, 1.666, 24.93, 2.232 /)

!real part
real(8), dimension(6), parameter:: amplitude_re=(/ 1.696, 6.351, 1.487, 0.415, 0.2371, 0.0 /)
real(8), dimension(6), parameter:: w_band_re=(/ 1.94, 2.45, 5.5, 7.54, 11.3, 0.0 /)
real(8), dimension(6), parameter:: sigma_re=(/ 2.666, 1.089, 3.633, 0.9662, 1.772, 0.0 /)


contains

function Drude(x, param)
!Drude-model permettivity function eps(iw)
!param - vector of parameters (wp, damp)

real(8):: x, Drude, param(2)

	Drude=1+param(1)**2/(x**2+param(2)*x);

end function
!----------------------------------------------------

function Drude_Re_Im(x,param,marker)
!real/imaginary part of dielectric perm. function, Drude model
!param - vector of parameters (wp, damp)

!marker (integer) allows to fix real or imaginary part
!marker=0 - real part
!marker=1 - imaginary part

real(8):: x, Drude_Re_Im
real(8):: param(2)
complex(8):: Drude_mod, i
integer:: marker

i=(0.0,1.0)

Drude_mod=cmplx(1.0,0.0,8) - cmplx(param(1)**2,0.0,8)/(cmplx(x**2,0.0,8) + i*x*param(2))

if (marker.eq.0) then

Drude_Re_Im=real(Drude_mod)

elseif (marker.eq.1) then
Drude_Re_Im=aimag(Drude_mod)

else
print*, "Incorrect marker. Please, use marker=0 to get real part and marker=1 - imaginary"

endif

end function

!----------------------------------------------------

function epsMar(x)
!Marachevsky formula fo eps(iw)
real(8):: x,epsMar,num,denum, eV
real(8), dimension(2,4):: paramMar

eV=1.519e15
!Marachevsky model parameters (wl1,wl2//gl1,gl2//gt1,gt2//wt2,0)
paramMar(1,1)=exp(-0.96)*1e16
paramMar(2,1)=exp(0.2866)*1e16
paramMar(1,2)=exp(-2.536)*1e16
paramMar(2,2)=exp(1.255)*1e16
paramMar(1,3)=exp(-4.7922)*1e16
paramMar(2,3)=exp(-0.957)*1e16
paramMar(1,4)=exp(-0.8359)*1e16
paramMar(2,4)=0.0_8

paramMar=paramMar/eV

num=(paramMar(1,1)**2+x**2+paramMar(1,2)*x)*(paramMar(2,1)**2+x**2+paramMar(2,2)*x)
denum=(x**2+paramMar(1,3)*x)*(paramMar(1,4)**2+x**2+paramMar(2,3)*x)

epsMar=num/denum

end function

!---------------------------------------------------

function oscil_gold_re_im(x,a,b,c,n,marker)
!real/imaginary part of oscillator approximation (Mostepanenko/Racik/other...) eps(w)

!you can fix your own parameters
!a - vector of oscillator strengths
!b - vector of relaxation parameters
!c - vector of resonant frequancies

!marker (integer) allows to fix real or imaginary part
!marker=0 - real part
!marker=1 - imaginary part

!n - dimension of a,b,c

integer:: n
real(8):: oscil_gold_re_im, x, a(n), b(n), c(n)
integer:: marker
complex(8):: i, oscil_mod
integer:: k

i=(0.0,1.0)
oscil_mod=(0.0,0.0)

do k=1,n
oscil_mod=oscil_mod+a(k)/(c(k)**2-x**2-i*x*b(k))
end do

if (marker.eq.0) then

oscil_gold_re_im=real(oscil_mod)

elseif (marker.eq.1) then
oscil_gold_re_im=aimag(oscil_mod)

else
print*, "Incorrect marker. Please, use marker=0 to get real part and marker=1 - imaginary"

endif


end function
!--------------------------------------------------

function oscil_gold(x,a,b,c)
!Lorentz function (oscillators) eps(iw)

!you can fix your own parameters
!a - vector of oscillator strengths
!b - vector of relaxation parameters
!c - vector of resonant frequancies

real(8):: oscil_gold, x, a(6), b(6), c(6)

oscil_gold=0.0_8

do i=1,size(a)
oscil_gold=oscil_gold+a(i)/(c(i)**2+x**2+b(i)*x)
end do

end function

!--------------------------------------------------

function oscil_to_int(x,a,b,c)
!KK expression under integral for oscillators

real(8):: x, oscil_to_int
real(8):: a(6),b(6),c(6)

oscil_to_int=x/(x**2+freqA**2)*oscil_gold(x,a,b,c)

end function

!--------------------------------------------------

function drude_to_int(x)
!KK expression under integral for Drude model

real(8):: x, drude_to_int, param(2)

param=parAproxIm
drude_to_int=param(1)**2*param(2)/((x**2+param(2)**2)*(x**2+freqA**2))

end function

!-------------------------------------------------

function Plasma(x,param,marker)
!Plasma model eps(w)/eps(iw)
!param - model parameters (wp,damp)

!marker (integer) fix one of two possible cases
!marker=0 <--> eps(w)
!marker=1 <--> eps(iw)

real(8):: Plasma,x,param(2)
integer:: marker

if (marker.eq.0) then

Plasma=1.0_8-param(1)**2/x**2

elseif (marker.eq.1) then
Plasma=1.0_8+param(1)**2/x**2

else
print*, "Incorrect marker. Please, use marker=0 to get eps(w) and marker=1 - eps(iw)"

endif


end function

!-------------------------------------------------

function Gen_Plasma(x,param,a,b,c)
!Generalized plasma-like dielectric permettivity function eps(iw)

!you can fix your own parameters
!a - vector of oscillator strengths
!b - vector of relaxation parameters
!c - vector of resonant frequancies
!param - model parameters (wp,damp)

real(8):: x, Gen_Plasma
real(8):: a(6), b(6), c(6), param(2)

Gen_Plasma=Plasma(x,param,1) + oscil_gold(x,a,b,c)

end function

!-------------------------------------------------

function Lorentz_Drude(x,param,a,b,c)
!Lorentz-Drude model eps(iw)

!you can fix your own parameters
!a - vector of oscillator strengths
!b - vector of relaxation parameters
!c - vector of resonant frequancies

real(8):: x, Lorentz_Drude
real(8):: a(6), b(6), c(6), param(2)

Lorentz_Drude=Drude(x, param) + oscil_gold(x,a,b,c)

end function

!-------------------------------------------------

function Brendel_Bormann_oscil(x,a,b,c,d)
!Brendel-Bormann model, oscillation part in eps(iw)

!you can fix your own parameters
!a - vector of oscillator strengths
!b - vector of relaxation parameters
!c - vector of resonant frequancies
!d - vector of sigma (dispersion)

real, parameter:: pi=3.14159265359

real(8):: x, Brendel_Bormann_oscil, a(6),b(6),c(6),d(6)
real(8):: WvalueRe1, WvalueIm1, WvalueRe2, WvalueIm2
logical:: logvar
real(8):: Re, Re1, Im
integer:: k
complex(8):: i, alpha, const
complex(8):: Brendel_Bormann_oscil_complex

!argW=Re+i*Im

i=(0.0_8, 1.0_8)
Brendel_Bormann_oscil_complex=(0.0, 0.0)

do k=1,size(a)

!calculate argument
alpha=sqrt(cmplx(x*b(k),0.0, 8) - cmplx(x**2,0.0, 8))

Re=real(alpha)-c(k)/(sqrt(2.0)*d(k))
Re1=real(alpha)+c(k)/(sqrt(2.0)*d(k))

Im=aimag(alpha)/(sqrt(2.0)*d(k))

!call Faddeeva function w(z)
call wofz(Re, Im, WvalueRe1, WvalueIm1, logvar) !first summand
call wofz(Re1, Im, WvalueRe2, WvalueIm2, logvar) !second summand


!calculate prefactor
const=i*sqrt(pi)*a(k)/(2*sqrt(2.0)*d(k)*alpha)

if (logvar.eqv. .True.) then
print*, 'Attention, error in wofz subroutine'
else
Brendel_Bormann_oscil_complex=Brendel_Bormann_oscil_complex+(cmplx(WvalueRe1+WvalueRe2,0.0,8)+i*(WvalueIm1+WvalueIm2))*const
end if

end do

Brendel_Bormann_oscil=real(Brendel_Bormann_oscil_complex)

end function

!-------------------------------------------------

function Brendel_Bormann_oscil_Re_Im(x,a,b,c,d,n,marker)
!Brendel-Bormann model oscillations, real/imaginary part of eps(w)

!you can fix your own parameters
!a - vector of oscillator strengths
!b - vector of relaxation parameters
!c - vector of resonant frequancies
!d - vector of sigma (dispersion)

!marker (integer) allows to fix real or imaginary part
!marker=0 - real part
!marker=1 - imaginary part

!n - dimension of a,b,c,d
real, parameter:: pi=3.14159265359

integer:: n
real(8):: x, Brendel_Bormann_oscil_Re_Im, a(n),b(n),c(n),d(n)
real(8):: WvalueRe1, WvalueIm1, WvalueRe2, WvalueIm2
logical:: logvar
real(8):: Re, Re1, Im
integer:: k
integer:: marker
complex(8):: i, alpha, const
real(8):: alpha1,alpha2
complex(8):: Brendel_Bormann_oscil_complex

!argW=Re+i*Im

i=(0.0, 1.0)
Brendel_Bormann_oscil_complex=(0.0, 0.0)

do k=1,n

!calculate argument
!alpha=sqrt(cmplx(x**2,0.0, 8) - i*cmplx(x*b(k),0.0, 8))

alpha1=x/sqrt(2.0)*sqrt(sqrt(1.0+(b(k)/x)**2)+1.0)
alpha2=x/sqrt(2.0)*sqrt(sqrt(1.0+(b(k)/x)**2)-1.0)

alpha=cmplx(alpha1,alpha2,8)

!Re=real(alpha)-c(k)/(sqrt(2.0)*d(k))
!Re1=real(alpha)+c(k)/(sqrt(2.0)*d(k))

Re=(alpha1-c(k))/(sqrt(2.0)*d(k))
Re1=(alpha1+c(k))/(sqrt(2.0)*d(k))

Im=alpha2/(sqrt(2.0)*d(k))

!call Faddeeva function w(z)
call wofz(Re, Im, WvalueRe1, WvalueIm1, logvar) !first summand
call wofz(Re1, Im, WvalueRe2, WvalueIm2, logvar) !second summand


!calculate prefactor
const=i*sqrt(pi)*a(k)/(2*sqrt(2.0)*d(k)*alpha)

if (logvar.eqv. .True.) then
print*, 'Attention, error in wofz subroutine'
else
Brendel_Bormann_oscil_complex=Brendel_Bormann_oscil_complex + (cmplx(WvalueRe1+WvalueRe2,0.0,8) + i*(WvalueIm1+WvalueIm2))*const
end if

end do

if (marker.eq.0) then

Brendel_Bormann_oscil_Re_Im=real(Brendel_Bormann_oscil_complex)

elseif (marker.eq.1) then
Brendel_Bormann_oscil_Re_Im=aimag(Brendel_Bormann_oscil_complex)

else
print*, "Incorrect marker. Please, use marker=0 to get real part and marker=1 - imaginary"

endif

end function
!-------------------------------------------------
function Brendel_Bormann(x,param,a,b,c,d)
!Brendel-Bormann model, eps(iw)

real(8):: Brendel_Bormann, x
real(8):: a(6),b(6),c(6),d(6), param(2)

!you can fix your own parameters
!a - vector of oscillator strengths
!b - vector of relaxation parameters
!c - vector of resonant frequancies
!d - vector of sigma (dispersion)

Brendel_Bormann=Drude(x, param)+Brendel_Bormann_oscil(x, a, b, c, d)

end function

!-------------------------------------------------

!-------------------------------------------------

function Gauss(x,a,w0,sigma)
!Gauss function
!a - amplitude
!w0 - resonans frequency
!sigma^2 - dispersion

real(8):: x, a, w0, sigma
real(8):: Gauss

Gauss=a*exp(-(x-w0)**2/sigma**2)

end function

!-------------------------------------------------

function Gauss_aprox(x,a,b,c,marker)
!Approximation as a sum of 6 Gauss contours

real(8):: x, a(6), b(6), c(6)
real(8):: Gauss_aprox
integer:: marker
integer:: i

Gauss_aprox=0.0D+00

if (marker .eq. 1) then
    do i=1,6
    Gauss_aprox = Gauss_aprox + Gauss(x,a(i),b(i),c(i)) + Gauss(x, a(i),-b(i),c(i))
    end do
        elseif (marker .eq. 2) then
    do i=1,6
    Gauss_aprox = Gauss_aprox + Gauss(x,a(i),b(i),c(i)) + Gauss(x, -a(i),-b(i),c(i))
    end do
        else
    print*, 'marker must be equal 1 for real part, 2 for imaginary. Please, try again!'
endif

end function

!-------------------------------------------------

function gauss_kramers_int(x)
!function to integrate in kramers-kronig to find real part from imaginary

real(4):: gauss_kramers_int
real(4):: x

gauss_kramers_int=x*Gauss_aprox(real(x,8), amplitude_im, w_band_im, sigma_im, 2)/(x+freqA)

end function

!-------------------------------------------------
function gauss_kramers_int1(x)
!function to integrate in kramers-kronig to find imaginary part from real

real(4):: gauss_kramers_int1
real(4):: x

gauss_kramers_int1=real(Gauss_aprox(real(x,8), amplitude_re, w_band_re, sigma_re, 1),4)

end function
!-------------------------------------------------

function eps_gauss_iomega(x,a,b,c,n)
!kramers_kronig from gauss sum

!you can fix your own parameters
!a - vector of amplitudes
!b - vector of resonant frequencies
!c - vector of dispersions

!n - dimension of a,b,c
real, parameter:: pi=3.14159265359

integer:: n
real(8):: x, eps_gauss_iomega, a(n),b(n),c(n)
real(8):: WvalueRe1, WvalueIm1, WvalueRe2, WvalueIm2
logical:: logvar
real(8):: Re, Re1, Im
integer:: k
complex(8):: i
complex(8):: complex_result

!argW=Re+i*Im

i = (0.0, 1.0)
comlex_result = (0.0, 0.0)
eps_gauss_iomega = 0.0D+00

do k=1,n

!calculate argument

Re=b(k)/c(k)
Re1=-b(k)/c(k)

Im=x/c(k)

!call Faddeeva function w(z)
call wofz(Re, Im, WvalueRe1, WvalueIm1, logvar) !first summand
call wofz(Re1, Im, WvalueRe2, WvalueIm2, logvar) !second summand


if (logvar.eqv. .True.) then
print*, 'Attention, error in wofz subroutine'
else
complex_result = complex_result + (cmplx(WvalueRe1-WvalueRe2,0.0,8) + i*(WvalueIm1-WvalueIm2))
end if

eps_gauss_iomega = eps_gauss_iomega + a(k)*aimag(complex_result)

end do

end function
!---------------------------------------------------

!*****************************************
!Want to calculate zero summand
!*****************************************

!-------------------------------------------------
function zero_sum_Dr(x)
!does not depend on distance if we have done change x --> dist*x
!Drude, Marachevsky, Kramers-Kr

real(4):: zero_sum_Dr, x
real(4):: rTm, rTe

rTe=0.0
rTm=1.0

if (x.eq.0.0) then
    zero_sum_Dr = 0.0
    else
    zero_sum_Dr = x*(log(1._8-rTm**2*exp(-2.*x))+log(1._8-rTe**2*exp(-2.*x)))
endif

end function

!----------------------------------------------------------------------

function zero_sum_Pl(a,x)
!Generalized plama-like model

real(4):: zero_sum_Pl, x
real(8):: rTm, rTe, a, log1, log2
real(8):: wp,cEv

cEv=19.746e-8 ! eV*s
wp=parMost(1)

rTm=1._8

rTe=((x/a)-sqrt((x/a)**2+(wp/cEv)**2))/((x/a)+sqrt((x/a)**2+(wp/cEv)**2))

log1 = log(1._8-rTe**2*exp(-2.*x))
log2 = log(1._8-rTm**2*exp(-2.*x))

zero_sum_Pl = x*(log1+log2)

end function

!---------------------------------------------------------------------

function zero_sum_int_Pl(x)
!Plasma zero summand function to integrate from 0 to 1
use dopcasimir

real(4):: zero_sum_int_Pl, x

if (x.eq.(0.0)) then
    zero_sum_int_Pl=0._8
    else
    zero_sum_int_Pl = zero_sum_Pl(dist,x)
endif

end function
!---------------------------------------------------------------------


end module



