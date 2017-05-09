program casimir_gold

use dopcasimir
use silicon
use gold
use force_module
use nickel
use titanium

implicit none

!counters
integer:: i,j,k,l
!---------------------------------------

integer:: N,stepn
!N - number of distance points
!stepn - number of Matsubara frequencies
real(8):: T, wmax
!T - temerature
real(8):: res, rest0, res1
real(8):: radius, coef, radius1, coef1
real(8):: start_point, final_point
real(8):: step

!casimir energy (1st column)/force(2nd column) variables
real(8), allocatable:: casimir_res_kk(:,:)
real(8), allocatable:: casimir_res_dr(:,:)
real(8), allocatable:: casimir_res_pl(:,:)
real(8), allocatable:: casimir_res_ldr(:,:)
real(8), allocatable:: casimir_res_bb(:,:)
real(8), allocatable:: casimir_res_gauss(:,:)
real(8), allocatable:: casimir_res_mar(:,:)
real(8), allocatable:: casimir_res_original(:,:)
real(8), allocatable:: casimir_res_dec1(:), casimir_res_dec2(:)
real(8), allocatable:: casimir_res_dec_dr(:), casimir_res_dec_pl(:)

!experimental data
real(8):: dec_exp_2007(17,2)
real(8):: dec_exp_2016(99,2)
real(8):: lamoreaux_exp(21,2)
integer:: expnumb
real(8), allocatable:: x(:)

!zero summands
real(4):: sume_dr, sume_pl, sumf_dr, sumf_pl
real(4):: sume_dr_dec1, sume_pl_dec1, sume_dr_dec2, sume_pl_dec2

!--------------------------------------------------------------------
!gold permittivity variables
integer,parameter:: ng=663 !rows number in file
!experimental data table Au
real(8):: matrixAu(ng,3)
!permittivity vectors Au
integer,parameter:: nni=353
!experimental data table Ni
real(8):: matrixNi(nni,3)
!permittivity vectors Ni
integer,parameter:: nti=286
!experimental data table Ti
real(8):: matrixTi(nti,3)
!permittivity vectors Ti

real(8), allocatable:: epsAur(:), epsAurDr(:), epsAurMar(:), epsAurGenPl(:), epsAurLD(:), epsAurBB(:), epsAurGauss(:)
real(8), allocatable:: epsNiDr(:), epsTiDr(:)
!integration results
real(8):: integralA1,integralA2
!variables needed for Kr-Kr integral in cspint subroutine
real(8)::funvalAu(ng),SpcoefA(3,ng), exintA(ng), workA(ng)
real(8)::funvalNi(nni),SpcoefNi(3,nni), exintNi(nni), workNi(nni)
real(8)::funvalTi(nti),SpcoefTi(3,nti), exintTi(nti), workTi(nti)
!variables for Kr-Kr integral in qnc79 subroutine
integer(4)::nevalA,ierA
!variables for qagi
real(4):: abserr, result
integer(4):: neval, ier

real(8):: a(2), b(2)


!-----------------------------------------------------------------------
open(unit=11, file='casimirgold.txt', status='replace')
!casimirgold.txt - file with computational results

!title
write(11,150) 'a, micro m', 'Kramers-Kr, nJ', 'Drude', 'Marachevsky', 'Gen_Plasma', &
'Drude-Lor', 'Bren-Borm', 'Gauss', 'T=0 Gen_Pl', 'Casimir'

!-------------------------------------------------------------------------

open(unit=13, file='casimir_force_gold.txt', status='replace')
!casimir_force_gold.txt - file with computational results

write(13,150) 'a, micro m', 'Kramers-Kr, nN', 'Drude', 'Marachevsky', 'Gen_Plasma', &
'Drude-Lor', 'Bren-Borm', 'Gauss', 'Casimir'

!-----------------------------------------------------------------------
open(unit=12, file='eta_to_plot_gold_energy.txt', status='replace')
!file with eta results

!title
write(12,150) 'a, micro m', 'Kramers-Kr', 'Drude', 'Marachevsky', 'Gen_Plasma', &
'Drude-Lor', 'Bren-Borm', 'Gauss', 'T=0 Gen_Pl'

!-----------------------------------------------------------------------
open(unit=14, file='eta_to_plot_gold_force.txt', status='replace')
!file with eta results

!title
write(14,150) 'a, micro m', 'Kramers-Kr', 'Drude', 'Marachevsky', 'Gen_Plasma', &
'Drude-Lor', 'Bren-Borm', 'Gauss'

!-----------------------------------------------------------------------
open(unit=20, file='casimir_lamoreaux.txt', status='replace')
!file with force values in order to compare with Lamoreaux

write(20,150) 'a, nm', 'Kramers-Kr', 'Drude', 'Marachevsky', 'Gen_Plasma', &
'Drude-Lor', 'Bren-Borm', 'Gauss', 'T=0 Gen_Pl', 'Casimir'

!-----------------------------------------------------------------------
open(unit=26, file='casimir_decca_2007.txt', status='replace')
!file with force values in order to compare with Decca 2007

write(26,150) 'a, nm', 'Kramers-Kr, mPa', 'Drude', 'Marachevsky', 'Gen_Plasma', &
'Drude-Lor', 'Bren-Borm', 'Gauss', 'Casimir'

!-----------------------------------------------------------------------
open(unit=25, file='casimir_decca.txt', status='replace')
!file with force values in order to compare with Decca 2016

write(25,150) 'a, nm', 'Drude, fN', 'Gen_Plasma'

!-----------------------------------------------------------------------
open(unit=30, file='decca-2007-paper-data.txt', status='old')
!file with experimental data (Decca 2007)

read(30,*)

do k=1,17
    read(30,*) (dec_exp_2007(k,j),j=1,2)
end do

!-----------------------------------------------------------------------
open(unit=31, file='Decca-2016-06-au-37.csv', status='old')
!file with experimental data (Decca 2016)

read(31,*)

do k=1,99
    read(31,*) (dec_exp_2016(k,j),j=1,2)
end do

!-----------------------------------------------------------------------
open(unit=32, file='lamoreaux-2010-fig2.csv', status='old')
!file with experimental data (Lamoreaux)

read(32,*)

do k=1,21
    read(32,*) (lamoreaux_exp(k,j),j=1,2)
end do


!-----------------------------------------------------------------------
!code for calulation of casimir free energy

open(unit=15, file='gold_eps_im_re_ev+Olmon.txt', status='old')
!read matrix from file
!1st column - frequencies
!2nd column - real part of permittivity
!3rd column - imaginary part of permittivity

read(15,*)

do i=1,ng
read(15,*), (matrixAu(i,j),j=1,3)
end do

open(unit=21, file='tableNi_eV_eps_Re_Im.txt', status='old')
!read matrix from file
!1st column - frequencies
!2nd column - real part of permittivity
!3rd column - imaginary part of permittivity

read(21,*)

do i=1,nni
read(21,*), (matrixNi(nni+1-i,j),j=1,3)
end do

open(unit=22, file='tableTi_eV_eps_Re_Im.txt', status='old')
!read matrix from file
!1st column - frequencies
!2nd column - real part of permittivity
!3rd column - imaginary part of permittivity

read(22,*)

do i=1,nti
read(22,*), (matrixTi(nti+1-i,j),j=1,3)
end do


T=300._8
!T - Temperature of the material

!epxperimental parameters (Decca)
d=1.0e-8
t_Au=37.0e-9

eps(3)=1._8
!fix the dielectric permettivity of the gap (vacuum)

!print*, 'type starting point in nanometers'
!read *, start_point
!
!print*, 'type final point in nanometers'
!read*, final_point
!
!if (final_point .le. start_point) then
!    print*, 'final point must be grater then starting, try again'
!    read*, final_point
!end if
!
!print *, 'type number of points '
!read *, N
!!N - number of points at the plot
!
!step = (final_point - start_point)*1.0e-9/N
!start_point = start_point * 1.0e-9

print*, 'choose the number of experiment you want to compare data with:'
print*, '1 - Decca 2007'
print*, '2 - Decca 2016'
print*, '3 - Lamoreaux'

read*, expnumb

if (expnumb .eq. 1) then
    N=17
    allocate(x(N))
    x=dec_exp_2007(1:17,1)*1e-9
elseif (expnumb .eq. 2) then
    N=99
    allocate(x(N))
    x=dec_exp_2016(1:99,1)*1e-9
elseif (expnumb .eq. 3) then
    N=21
    allocate(x(N))
    x=lamoreaux_exp(1:21,1)*1e-6
else
    print*, 'The number must be equal 1,2 or 3! Please, try again.'
endif


allocate(casimir_res_kk(N,2))
casimir_res_kk = 0.0
allocate(casimir_res_mar(N,2))
casimir_res_mar = 0.0
allocate(casimir_res_dr(N,2))
casimir_res_dr = 0.0
allocate(casimir_res_pl(N,2))
casimir_res_pl = 0.0
allocate(casimir_res_ldr(N,2))
casimir_res_ldr = 0.0
allocate(casimir_res_bb(N,2))
casimir_res_bb = 0.0
allocate(casimir_res_gauss(N,2))
casimir_res_gauss = 0.0
allocate(casimir_res_original(N,2))
casimir_res_original = 0.0
allocate(casimir_res_dec1(N))
casimir_res_dec1 = 0.0
allocate(casimir_res_dec2(N))
casimir_res_dec2 = 0.0
allocate(casimir_res_dec_dr(N))
allocate(casimir_res_dec_pl(N))

model(1)=4
model(2)=4

!    do j=1,N+1
!    !do-cycle for distance
!    dist=start_point+(j-1)*step

    do j=1,N
    dist=x(j)


    !first calculate zero-temperature casimir energy
    a=(/0._8, 0._8/)
    b=(/1._8, 1._8/)
    call monte_carlo_nd (zerotempint, 2, a,b, 1000500, res)
    call trapzd2d(zerotempint,a,b,res1,stepn)
    rest0=res
    call monte_carlo_nd (zerotempint1, 2,  a,b, 1000500, res)
    call trapzd2d(zerotempint1,a,b,res1,stepn)
    rest0=(rest0+res)*c/dist**3

    k=1
    wmax=stepn*cEv/dist
    w=2*pi*kb*T/(h*eV)
    do while (w.le.wmax)
    k=k+1
    w=2*k*pi*kb*T/(h*eV)
    end do

print*, 'k=', k

allocate(epsAur(k))
allocate(epsAurDr(k))
allocate(epsNiDr(k))
allocate(epsTiDr(k))
allocate(epsAurMar(k))
allocate(epsAurGenPl(k))
allocate(epsAurLD(k))
allocate(epsAurBB(k))
allocate(epsAurGauss(k))

        !calculate casimir energy
        do i=1,k
        !do-cycle for frequencies

            w=2*i*pi*kb*T/(h*eV)
            !w - Matsubara frequency
            freqA=w
!GOLD**************************************************************************************
    !integrates the Drude-model approximation
    call qnc79(drude_to_int, 0.0_8, matrixAu(1,1), 1.0e-3_8, integralA1, ierA, nevalA)

    !builds vector to integrate from the Palik data
        do l=1,ng
        funvalAu(l)=matrixAu(l,3)*matrixAu(l,1)/(matrixAu(l,1)**2+freqA**2)
        end do
    !integrates the vector of data
    call cspint(ng, matrixAu(1:ng,1), funvalAu, matrixAu(1,1), matrixAu(ng,1), SpcoefA, exintA, workA, integralA2)

    !Kr-Kr formula for permittivity
    epsAur(i)=(integralA1+integralA2)*2/pi + 1.0

!NICKEL**********************************************************************************
!integrates the Drude-model approximation
    call qnc79(drude_to_int_Ni, 0.0_8, matrixNi(1,1), 1.0e-3_8, integralA1, ierA, nevalA)

    !builds vector to integrate from the Palik data
        do l=1,nni
        funvalNi(l)=matrixNi(l,3)*matrixNi(l,1)/(matrixNi(l,1)**2+freqA**2)
        end do
    !integrates the vector of data
    call cspint(nni, matrixNi(1:nni,1), funvalNi, matrixNi(1,1), matrixNi(nni,1), SpcoefNi, exintNi, workNi, integralA2)

    !Kr-Kr formula for permittivity
    epsNiDr(i)=(integralA1+integralA2)*2/pi + 1.0

!TITANIUM********************************************************************************
!integrates the Drude-model approximation
    call qnc79(drude_to_int_Ti, 0.0_8, matrixTi(1,1), 1.0e-3_8, integralA1, ierA, nevalA)

    !builds vector to integrate from the Palik data
        do l=1,nti
        funvalTi(l)=matrixTi(l,3)*matrixTi(l,1)/(matrixTi(l,1)**2+freqA**2)
        end do
    !integrates the vector of data
    call cspint(nti, matrixTi(1:nti,1), funvalTi, matrixTi(1,1), matrixTi(nti,1), SpcoefTi, exintTi, workTi, integralA2)

    !Kr-Kr formula for permittivity
    epsTiDr(i)=(integralA1+integralA2)*2/pi + 1.0

!****************************************************************************************

    !create vectors of permittivity values within different models
    !Drude model
    epsAurDr(i)=Drude(w,parAproxIm)

    !Marachevsky
    epsAurMar(i)=epsMar(w)

    !Generalized plasma (Mostepanenko)
    epsAurGenPl(i)=Gen_Plasma(w, parMost, gAu, gammaAu, wAu)

    !Lorentz-Drude model (Rakic)
    epsAurLD(i)=Lorentz_Drude(w, parRakic_LD, fR*wpR**2, gammaR, wR)

    !Brenedel-Bormann (Rakic)
    epsAurBB(i)=Brendel_Bormann(w, parRakic_BB, fRb*wpR**2, gammaRb, wRb, sigmaRb)

    !Gauss (mine)
    epsAurGauss(i)=epsAurDr(i)+eps_gauss_iomega(w,amplitude_im,w_band_im,sigma_im,6)
!-----------------------------------------------------------------------

    !calculate casimir (Kramers-Kr)
    !fix dielectric model
    eps(1)=epsAur(i)
    eps(2)=epsAur(i)

    !calculate energy
    call qagi( to_int, 0.0, 1, 1.0e-5, 1.0e-5, result, abserr, neval, ier )
    casimir_res_kk(j,1)=casimir_res_kk(j,1) + result

    !calculate force
    call qagi( lif_to_int, 0.0, 1, 1.0e-10, 1.0e-8, result, abserr, neval, ier )
    casimir_res_kk(j,2)=casimir_res_kk(j,2) + result

!-----------------------------------------------------------------------

    !calculate casimir (Drude)
    !fix dielecric model
    eps(1)=epsAurDr(i)
    eps(2)=epsAurDr(i)

    !calculate energy
    call qagi( to_int, 0.0, 1, 1.0e-5, 1.0e-5, result, abserr, neval, ier )
    casimir_res_dr(j,1)=casimir_res_dr(j,1) + result

    !calculate force
    call qagi( lif_to_int, 0.0, 1, 1.0e-10, 1.0e-8, result, abserr, neval, ier )
    casimir_res_dr(j,2)=casimir_res_dr(j,2) + result

!----------------------------------------------------------------------

    !calculate casimir (Marachevsky)
    !fix dielectric model
    eps(1)=epsAurMar(i)
    eps(2)=epsAurMar(i)

    !calculate energy
    call qagi( to_int, 0.0, 1, 1.0e-5, 1.0e-5, result, abserr, neval, ier )
    casimir_res_mar(j,1)=casimir_res_mar(j,1) + result

    !calculate force
    call qagi( lif_to_int, 0.0, 1, 1.0e-10, 1.0e-8, result, abserr, neval, ier )
    casimir_res_mar(j,2)=casimir_res_mar(j,2) + result

!----------------------------------------------------------------------

    !calculate casimir (Generalized plasma)
    !fix dielectric model
    eps(1)=epsAurGenPl(i)
    eps(2)=epsAurGenPl(i)

    !calculate energy
    call qagi( to_int, 0.0, 1, 1.0e-5, 1.0e-5, result, abserr, neval, ier )
    casimir_res_pl(j,1)=casimir_res_pl(j,1) + result

    !calculate force
    call qagi( lif_to_int, 0.0, 1, 1.0e-10, 1.0e-8, result, abserr, neval, ier )
    casimir_res_pl(j,2)=casimir_res_pl(j,2) + result

!---------------------------------------------------------------------

    !calculate casimir (Lorentz-Drude)
    !fix dielectric model
    eps(1)=epsAurLD(i)
    eps(2)=epsAurLD(i)

    !calculate energy
    call qagi( to_int, 0.0, 1, 1.0e-5, 1.0e-5, result, abserr, neval, ier )
    casimir_res_ldr(j,1)=casimir_res_ldr(j,1) + result

    !calculate force
    call qagi( lif_to_int, 0.0, 1, 1.0e-10, 1.0e-8, result, abserr, neval, ier )
    casimir_res_ldr(j,2)=casimir_res_ldr(j,2) + result

!---------------------------------------------------------------------

    !calculate casimir (Brendel-Bormann)
    !fix dielectric model
    eps(1)=epsAurBB(i)
    eps(2)=epsAurBB(i)

    !calculate energy
    call qagi( to_int, 0.0, 1, 1.0e-5, 1.0e-5, result, abserr, neval, ier )
    casimir_res_bb(j,1)=casimir_res_bb(j,1) + result

    !calculate force
    call qagi( lif_to_int, 0.0, 1, 1.0e-10, 1.0e-8, result, abserr, neval, ier )
    casimir_res_bb(j,2)=casimir_res_bb(j,2) + result

!----------------------------------------------------------------------

!calculate casimir (Gauss)
    !fix dielectric model
    eps(1)=epsAurGauss(i)
    eps(2)=epsAurGauss(i)

    !calculate energy
    call qagi( to_int, 0.0, 1, 1.0e-5, 1.0e-5, result, abserr, neval, ier )
    casimir_res_gauss(j,1)=casimir_res_gauss(j,1) + result

    !calculate force
    call qagi( lif_to_int, 0.0, 1, 1.0e-10, 1.0e-8, result, abserr, neval, ier )
    casimir_res_gauss(j,2)=casimir_res_gauss(j,2) + result

!----------------------------------------------------------------------

    !calculate Casimir for Decca experiment
    epsNi = epsNiDr(i)
    epsTi = epsTiDr(i)
    epsAu = epsAurGauss(i)

call qagi( to_int_Au, 0.0, 1, 1.0e-5, 1.0e-5, result, abserr, neval, ier )
     casimir_res_dec1(j) = casimir_res_dec1(j) + result

call qagi( to_int_Ni, 0.0, 1, 1.0e-5, 1.0e-5, result, abserr, neval, ier )
     casimir_res_dec2(j) = casimir_res_dec2(j) + result

    end do

!****************************************************
!we calculated only 1,...,N summands, let's calculate zero summand

!energy----------------------------------------
!(Drude,Marachevsky,Kramers-Kr)

    call qagi( zero_sum_Dr, 0.0, 1, 1.0e-5, 1.0e-5, sume_dr, abserr, neval, ier )

!(Plasma)

    call qagi( zero_sum_int_Pl, 0.0, 1, 1.0e-5, 1.0e-5, sume_pl, abserr, neval, ier )

!force----------------------------------------
!zero summand for force (Drude)
    call qagi( zero_sum_force_Dr, 0.0, 1, 1.0e-10, 1.0e-8, sumf_dr, abserr, neval, ier )

!zero summand for force (Plasma)
    call qagi( zero_sum_int_force_Pl, 0.0, 1, 1.0e-10, 1.0e-8, sumf_pl, abserr, neval, ier )

!***************************************************
!multiply integrational energy results by necessary constants
! Kramers-Kr
casimir_res_kk(j,1)=1.602e-19*kb*T/(2*pi*dist**2)*(0.5*sume_dr + casimir_res_kk(j,1))
! Drude
casimir_res_dr(j,1)=1.602e-19*kb*T/(2*pi*dist**2)*(0.5*sume_dr + casimir_res_dr(j,1))
! Marachevsky
casimir_res_mar(j,1)=1.602e-19*kb*T/(2*pi*dist**2)*(0.5*sume_dr + casimir_res_mar(j,1))
! Generalized plasma
casimir_res_pl(j,1)=1.602e-19*kb*T/(2*pi*dist**2)*(0.5*sume_pl + casimir_res_pl(j,1))
! Drude-Lorentz
casimir_res_ldr(j,1)=1.602e-19*kb*T/(2*pi*dist**2)*(0.5*sume_dr + casimir_res_ldr(j,1))
! Brendel-Bormann
casimir_res_bb(j,1)=1.602e-19*kb*T/(2*pi*dist**2)*(0.5*sume_dr + casimir_res_bb(j,1))
! Gauss
casimir_res_gauss(j,1)=1.602e-19*kb*T/(2*pi*dist**2)*(0.5*sume_dr + casimir_res_gauss(j,1))

rest0=1.602e-19*h*rest0/(2*pi)**2

casimir_res_original(j,1)=-1.602e-19*(pi)**2*h*c/(720*dist**3)

!write energy data in casimirgold.txt
100 format(10(f15.3))
150 format(10(A15))

write(11,100) dist*1.0e6, abs((casimir_res_kk(j,1))*1.0e9), abs((casimir_res_dr(j,1))*1.0e9), &
abs((casimir_res_mar(j,1))*1.0e9), abs((casimir_res_pl(j,1))*1.0e9), abs((casimir_res_ldr(j,1))*1.0e9), &
abs((casimir_res_bb(j,1))*1.0e9), abs((casimir_res_gauss(j,1))*1.0e9), &
abs(rest0)*1.0e9, abs(casimir_res_original(j,1)*1.0e9)

!want to write eta = free energy / casimir result
write(12,100) dist*1.0e6, casimir_res_kk(j,1)/casimir_res_original(j,1),&
casimir_res_dr(j,1)/casimir_res_original(j,1), &
casimir_res_mar(j,1)/casimir_res_original(j,1), &
casimir_res_pl(j,1)/casimir_res_original(j,1), &
casimir_res_ldr(j,1)/casimir_res_original(j,1), &
casimir_res_bb(j,1)/casimir_res_original(j,1), &
casimir_res_gauss(j,1)/casimir_res_original(j,1), &
rest0/casimir_res_original(j,1)


!write data to compare with Lamoreaux
!radius of Au sphere R=15,6 cm=15.6e-2 m
radius=15.6e-2
coef=2*pi*radius*(dist*1.0e6)**2

!write(20,100) dist*1.0e6, coef*casimir_res_kk(j,1)*1.0e12, coef*casimir_res_dr(j,1)*1.0e12, &
!coef*casimir_res_mar(j,1)*1.0e12, coef*casimir_res_pl(j,1)*1.0e12, coef*casimir_res_ldr(j,1)*1.0e12, &
!coef*casimir_res_bb(j,1)*1.0e12, coef*casimir_res_gauss(j,1)*1.0e12, &
!coef*rest0*1.0e12, coef*casimir_res_original(j,1)*1.0e12

write(20,100) dist*1.0e6, &
coef*casimir_res_kk(j,1)*1.0e12+lamoreaux_exp(j,2), &
coef*casimir_res_dr(j,1)*1.0e12+lamoreaux_exp(j,2), &
coef*casimir_res_mar(j,1)*1.0e12+lamoreaux_exp(j,2), &
coef*casimir_res_pl(j,1)*1.0e12+lamoreaux_exp(j,2), &
coef*casimir_res_ldr(j,1)*1.0e12+lamoreaux_exp(j,2), &
coef*casimir_res_bb(j,1)*1.0e12+lamoreaux_exp(j,2), &
coef*casimir_res_gauss(j,1)*1.0e12+lamoreaux_exp(j,2), &
coef*rest0*1.0e12+lamoreaux_exp(j,2), &
coef*casimir_res_original(j,1)*1.0e12+lamoreaux_exp(j,2)


!***********************************************************************************************
!multiply integrational force results by necessary constants
! Kramers-Kr
casimir_res_kk(j,2)=-1.602e-19*kb*T/(pi*dist**3)*(0.5*sumf_dr + casimir_res_kk(j,2))
! Drude
casimir_res_dr(j,2)=-1.602e-19*kb*T/(pi*dist**3)*(0.5*sumf_dr + casimir_res_dr(j,2))
! Marachevsky
casimir_res_mar(j,2)=-1.602e-19*kb*T/(pi*dist**3)*(0.5*sumf_dr + casimir_res_mar(j,2))
! Generalized plasma
casimir_res_pl(j,2)=-1.602e-19*kb*T/(pi*dist**3)*(0.5*sumf_pl + casimir_res_pl(j,2))
! Drude-Lorentz
casimir_res_ldr(j,2)=-1.602e-19*kb*T/(pi*dist**3)*(0.5*sumf_dr + casimir_res_ldr(j,2))
! Brendel-Bormann
casimir_res_bb(j,2)=-1.602e-19*kb*T/(pi*dist**3)*(0.5*sumf_dr + casimir_res_bb(j,2))
! Gauss
casimir_res_gauss(j,2)=-1.602e-19*kb*T/(pi*dist**3)*(0.5*sumf_dr + casimir_res_gauss(j,2))

!Casimir original result
casimir_res_original(j,2)=-1.602e-19*(pi)**2*h*c/(240*dist**4)

write(13,100) dist*1.0e6, casimir_res_kk(j,2), casimir_res_dr(j,2), &
casimir_res_mar(j,2), casimir_res_pl(j,2), casimir_res_ldr(j,2), &
casimir_res_bb(j,2), casimir_res_gauss(j,2), casimir_res_original(j,2)

!want to write eta = free energy / casimir result
write(14,100) dist*1.0e6, casimir_res_kk(j,2)/casimir_res_original(j,2),&
casimir_res_dr(j,2)/casimir_res_original(j,2), &
casimir_res_mar(j,2)/casimir_res_original(j,2), &
casimir_res_pl(j,2)/casimir_res_original(j,2), &
casimir_res_ldr(j,2)/casimir_res_original(j,2), &
casimir_res_bb(j,2)/casimir_res_original(j,2), &
casimir_res_gauss(j,2)/casimir_res_original(j,2)

!write data to compare with Decca article 2007

!write(26,100) dist*1.0e9, casimir_res_kk(j,2)*1e3, casimir_res_dr(j,2)*1e3, &
!casimir_res_mar(j,2)*1e3, casimir_res_pl(j,2)*1e3, casimir_res_ldr(j,2)*1e3, &
!casimir_res_bb(j,2)*1e3, casimir_res_gauss(j,2)*1e3, casimir_res_original(j,2)*1e3

write(26,100) dist*1.0e9, &
casimir_res_kk(j,2)*1e3+dec_exp_2007(j,2), &
casimir_res_dr(j,2)*1e3+dec_exp_2007(j,2), &
casimir_res_mar(j,2)*1e3+dec_exp_2007(j,2), &
casimir_res_pl(j,2)*1e3+dec_exp_2007(j,2), &
casimir_res_ldr(j,2)*1e3+dec_exp_2007(j,2), &
casimir_res_bb(j,2)*1e3+dec_exp_2007(j,2), &
casimir_res_gauss(j,2)*1e3+dec_exp_2007(j,2), &
casimir_res_original(j,2)*1e3+dec_exp_2007(j,2)

!conduct calculations similar with that form Decca paper---------------------------------

!zero summand
!Drude - Nickel

    call qagi( zero_sum_Dr_Ni, 0.0, 1, 1.0e-5, 1.0e-5, sume_dr_dec2, abserr, neval, ier )

!Plasma - Nickel

    call qagi( zero_sum_Pl_Ni, 0.0, 1, 1.0e-5, 1.0e-5, sume_pl_dec2, abserr, neval, ier )

!Drude - Gold

    call qagi( zero_sum_Dr_AuAu, 0.0, 1, 1.0e-5, 1.0e-5, sume_dr_dec1, abserr, neval, ier )

!Plasma - Gold

    call qagi( zero_sum_Pl_AuAu, 0.0, 1, 1.0e-5, 1.0e-5, sume_pl_dec1, abserr, neval, ier )

!final
radius1=149.3e-6
coef1=2*pi*radius1

casimir_res_dec_dr(j) = (casimir_res_dec1(j) + 0.5*sume_dr_dec1) - (casimir_res_dec2(j) + 0.5*sume_dr_dec2)
casimir_res_dec_dr(j) = 1.602e-19*kb*T/(2*pi*dist**2)*coef1*casimir_res_dec_dr(j)

casimir_res_dec_pl(j) = (casimir_res_dec1(j) + 0.5*sume_pl_dec1) - (casimir_res_dec2(j) + 0.5*sume_pl_dec2)
casimir_res_dec_pl(j) = 1.602e-19*kb*T/(2*pi*dist**2)*coef1*casimir_res_dec_pl(j)

write(25,100) dist*1.0e9, &
casimir_res_dec_dr(j)*1.0e15-dec_exp_2016(j,2), &
casimir_res_dec_pl(j)*1.0e15-dec_exp_2016(j,2)


deallocate(epsAur)
deallocate(epsAurDr)
deallocate(epsNiDr)
deallocate(epsTiDr)
deallocate(epsAurMar)
deallocate(epsAurGenPl)
deallocate(epsAurLD)
deallocate(epsAurBB)
deallocate(epsAurGauss)

    end do


!**********************************************************************************************

deallocate(casimir_res_kk)
deallocate(casimir_res_dr)
deallocate(casimir_res_mar)
deallocate(casimir_res_ldr)
deallocate(casimir_res_bb)
deallocate(casimir_res_gauss)
deallocate(casimir_res_original)
deallocate(casimir_res_dec1)
deallocate(casimir_res_dec2)
deallocate(casimir_res_dec_dr)
deallocate(casimir_res_dec_pl)
deallocate(x)



close(11)
close(10)
close(12)
close(13)
close(14)
close(15)
close(20)
close(21)
close(22)
close(25)
close(26)

print*, 'Done!'

end program


	
