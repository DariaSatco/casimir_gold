program casimir_gold

use dopcasimir
use silicon
use gold
use force_module

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
real(8):: radius, coef

!casimir energy (1st column)/force(2nd column) variables
real(8), allocatable:: casimir_res_kk(:,:)
real(8), allocatable:: casimir_res_dr(:,:)
real(8), allocatable:: casimir_res_pl(:,:)
real(8), allocatable:: casimir_res_ldr(:,:)
real(8), allocatable:: casimir_res_bb(:,:)
real(8), allocatable:: casimir_res_mar(:,:)
real(8), allocatable:: casimir_res_original(:,:)

!zero summands
real(4):: sume_dr, sume_pl, sumf_dr, sumf_pl

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

real(8), allocatable:: epsAur(:), epsAurDr(:), epsAurMar(:), epsAurGenPl(:), epsAurLD(:), epsAurBB(:)
!integration results
real(8):: integralA1,integralA2
!variables needed for Kr-Kr integral in cspint subroutine
real(8)::funvalAu(ng),SpcoefA(3,ng), exintA(ng), workA(ng)
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
'Drude-Lor', 'Bren-Borm', 'T=0 Gen_Pl', 'Casimir'

!-------------------------------------------------------------------------

open(unit=13, file='casimir_force_gold.txt', status='replace')
!casimir_force_gold.txt - file with computational results

write(13,150) 'a, micro m', 'Kramers-Kr, nN', 'Drude', 'Marachevsky', 'Gen_Plasma', &
'Drude-Lor', 'Bren-Borm', 'Casimir'

!-----------------------------------------------------------------------
open(unit=12, file='eta_to_plot_gold_energy.txt', status='replace')
!file with eta results

!title
write(12,150) 'a, micro m', 'Kramers-Kr', 'Drude', 'Marachevsky', 'Gen_Plasma', &
'Drude-Lor', 'Bren-Borm', 'T=0 Gen_Pl'

!-----------------------------------------------------------------------
open(unit=14, file='eta_to_plot_gold_force.txt', status='replace')
!file with eta results

!title
write(14,150) 'a, micro m', 'Kramers-Kr', 'Drude', 'Marachevsky', 'Gen_Plasma', &
'Drude-Lor', 'Bren-Borm'

!-----------------------------------------------------------------------
open(unit=20, file='casimir_decca.txt', status='replace')
!file with force values in order to compare with Decca

write(20,150) 'a, nm', 'Kramers-Kr', 'Drude', 'Marachevsky', 'Gen_Plasma', &
'Drude-Lor', 'Bren-Borm', 'T=0 Gen_Pl', 'Casimir'

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
read(21,*), (matrixNi(i,j),j=1,3)
end do

open(unit=22, file='tableTi_eV_eps_Re_Im.txt', status='old')
!read matrix from file
!1st column - frequencies
!2nd column - real part of permittivity
!3rd column - imaginary part of permittivity

read(22,*)

do i=1,nti
read(22,*), (matrixTi(i,j),j=1,3)
end do


T=300._8
!T - Temperature of the material

eps(3)=1._8
!fix the dielectric permettivity of the gap (vacuum)

print *, 'type number of points '
read *, N
!N - number of points at the plot

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
allocate(casimir_res_original(N,2))
casimir_res_original = 0.0

model(1)=4
model(2)=4

    do j=1,N
    !do-cycle for distance
    dist=1.0e-7*j

    sume_dr=0.0
    sume_pl=0.0
    sumf_dr=0.0
    sumf_pl=0.0

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
allocate(epsAurMar(k))
allocate(epsAurGenPl(k))
allocate(epsAurLD(k))
allocate(epsAurBB(k))

        !calculate casimir energy
        do i=1,k
        !do-cycle for frequencies

            w=2*i*pi*kb*T/(h*eV)
            !w - Matsubara frequency
            freqA=w

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

    end do

!****************************************************
!we calculated only 1,...,N summands, let's calculate zero summand

!energy----------------------------------------
!ssum0 (Drude,Marachevsky,Kramers-Kr)

    call qagi( zero_sum_Dr, 0.0, 1, 1.0e-5, 1.0e-5, sume_dr, abserr, neval, ier )

!ssum01 (Plasma)

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

rest0=1.602e-19*h*rest0/(2*pi)**2

casimir_res_original(j,1)=-1.602e-19*(pi)**2*h*c/(720*dist**3)

!write energy data in casimirgold.txt
100 format(10(f15.3))
150 format(10(A15))

write(11,100) dist*1.0e6, abs((casimir_res_kk(j,1))*1.0e9), abs((casimir_res_dr(j,1))*1.0e9), &
abs((casimir_res_mar(j,1))*1.0e9), abs((casimir_res_pl(j,1))*1.0e9), abs((casimir_res_ldr(j,1))*1.0e9), &
abs((casimir_res_bb(j,1))*1.0e9), abs(rest0)*1.0e9, abs(casimir_res_original(j,1)*1.0e9)

!want to write eta = free energy / casimir result
write(12,100) dist*1.0e6, casimir_res_kk(j,1)/casimir_res_original(j,1),&
casimir_res_dr(j,1)/casimir_res_original(j,1), &
casimir_res_mar(j,1)/casimir_res_original(j,1), &
casimir_res_pl(j,1)/casimir_res_original(j,1), &
casimir_res_ldr(j,1)/casimir_res_original(j,1), &
casimir_res_bb(j,1)/casimir_res_original(j,1), &
rest0/casimir_res_original(j,1)

!write data to compare with Decca
!radius of Au sphere R=149,3 micro m=149.3e-6 m
radius=149.3e-6
coef=2*pi*radius

write(*,100) dist*1.0e9, coef*casimir_res_kk(j,1)*1.0e15, coef*casimir_res_dr(j,1)*1.0e15, &
coef*casimir_res_mar(j,1)*1.0e15, coef*casimir_res_pl(j,1)*1.0e15, coef*casimir_res_ldr(j,1)*1.0e15, &
coef*casimir_res_bb(j,1)*1.0e15, coef*rest0*1.0e15, coef*casimir_res_original(j,1)*1.0e15

!write data to compare with Lamoreaux
!radius of Au sphere R=15,6 cm=15.6e-2 m
radius=15.6e-2
coef=2*pi*radius*(dist*1.0e6)**2

write(20,100) dist*1.0e6, coef*casimir_res_kk(j,1)*1.0e12, coef*casimir_res_dr(j,1)*1.0e12, &
coef*casimir_res_mar(j,1)*1.0e12, coef*casimir_res_pl(j,1)*1.0e12, coef*casimir_res_ldr(j,1)*1.0e12, &
coef*casimir_res_bb(j,1)*1.0e12, coef*rest0*1.0e12, coef*casimir_res_original(j,1)*1.0e12

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

!Casimir original result
casimir_res_original(j,2)=-1.602e-19*(pi)**2*h*c/(240*dist**4)

write(13,100) dist*1.0e6, casimir_res_kk(j,2), casimir_res_dr(j,2), &
casimir_res_mar(j,2), casimir_res_pl(j,2), casimir_res_ldr(j,2), &
casimir_res_bb(j,2), casimir_res_original(j,2)

!want to write eta = free energy / casimir result
write(14,100) dist*1.0e6, casimir_res_kk(j,2)/casimir_res_original(j,2),&
casimir_res_dr(j,2)/casimir_res_original(j,2), &
casimir_res_mar(j,2)/casimir_res_original(j,2), &
casimir_res_pl(j,2)/casimir_res_original(j,2), &
casimir_res_ldr(j,2)/casimir_res_original(j,2), &
casimir_res_bb(j,2)/casimir_res_original(j,2)


deallocate(epsAur)
deallocate(epsAurDr)
deallocate(epsAurMar)
deallocate(epsAurGenPl)
deallocate(epsAurLD)
deallocate(epsAurBB)

    end do


!**********************************************************************************************

deallocate(casimir_res_kk)
deallocate(casimir_res_dr)
deallocate(casimir_res_mar)
deallocate(casimir_res_ldr)
deallocate(casimir_res_bb)
deallocate(casimir_res_original)


close(11)
close(10)
close(12)
close(13)
close(14)
close(15)
close(20)
close(21)
close(22)

print*, 'Done!'

end program


	
