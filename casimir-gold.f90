program casimir_gold

use dopcasimir
use silicon
use gold

implicit none

!counters
integer:: i,j,k,l
!---------------------------------------
!casimir energy variables

integer:: N,stepn
real(8):: res, rest0, res1
real(8), allocatable:: fenergy(:), fenergy1(:), fenergy2(:), fenergy3(:), fenergy4(:), fenergy5(:), ecasimir(:)
real(8), allocatable:: ssum(:), ssum1(:), ssum2(:), ssum3(:), ssum4(:), ssum5(:)
real(8):: ssum0, ssum01
real(8):: T, wmax

!--------------------------------------------------------------------
!gold permettivity variables
integer,parameter:: ng=310 !rows number in file
!Palik data table
real(8):: matrixAu(ng,3)
!permettivity vector
real(8), allocatable:: epsAur(:), epsAurDr(:), epsAurMar(:), epsAurGenPl(:), epsAurLD(:), epsAurBB(:)
!integration results
real(8):: integralA1,integralA2
!variables needed for Kr-Kr integral in cspint subroutine
real(8)::funvalAu(ng),SpcoefA(3,ng), exintA(ng), workA(ng)
!variables for Kr-Kr integral in qnc79 subroutine
integer(4)::nevalA,ierA

real(8):: a(2), b(2)


!-----------------------------------------------------------------------
open(unit=11, file='casimirgold.txt', status='replace')
!casimirgold.txt - file with computational results

!title
write(11,150) 'a, micro m', 'Kramers-Kr', 'Drude', 'Marachevsky', 'Gen_Plasma', &
'Drude-Lor', 'Bren-Borm', 'T=0 Gen_Pl', 'Casimir'

!-----------------------------------------------------------------------
open(unit=12, file='eta_to_plot_gold.txt', status='replace')
!file with eta results

!title
write(12,150) 'a, micro m', 'Kramers-Kr', 'Drude', 'Marachevsky', 'Gen_Plasma', &
'Drude-Lor', 'Bren-Borm', 'T=0 Gen_Pl'

!-----------------------------------------------------------------------
!code for calulation of casimir free energy

open(unit=10, file='resultAu_eV.txt', status='old')
!read matrix from file (gold)
!1st column - frequencies
!2nd column - real part of permittivity
!3rd column - imaginary part of permittivity

do i=1,ng
read(10,*), (matrixAu(i,j),j=1,3)
end do

T=300._8
!T - Temperature of the material

eps(3)=1._8
!fix the dielectric permettivity of the gap (vacuum)

print *, 'type number of points '
read *, N
! N - number of points on the plot

allocate(ssum(N))
allocate(fenergy(N))
allocate(ssum1(N))
allocate(fenergy1(N))
allocate(ssum2(N))
allocate(fenergy2(N))
allocate(ssum3(N))
allocate(fenergy3(N))
allocate(ssum4(N))
allocate(fenergy4(N))
allocate(ssum5(N))
allocate(fenergy5(N))
allocate(ecasimir(N))

model(1)=4
model(2)=4

    do j=1,N
    !do-cycle for distance
    dist=1.0e-7*j
    ssum(j)=0._8
    ssum1(j)=0._8
    ssum2(j)=0._8
    ssum3(j)=0._8
    ssum4(j)=0._8
    ssum5(j)=0._8

    ssum0=0._8
    ssum01=0._8

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

    epsAurDr(i)=Drude(w,parRakic_LD)

    epsAurMar(i)=epsMar(w)

    epsAurGenPl(i)=Gen_Plasma(w, parMost, gAu, gammaAu, wAu)

    epsAurLD(i)=Lorentz_Drude(w, parRakic_LD, fR*wpR**2, gammaR, wR)

    epsAurBB(i)=Brendel_Bormann(w, parRakic_BB, fRb*wpR**2, gammaRb, wRb, sigmaRb)
!-----------------------------------------------------------------------
    !calculate casimir (Kramers-Kr)
    eps(1)=epsAur(i)
    eps(2)=epsAur(i)
    call trapzd(to_int,0._8,1._8,res)

            ssum(j)=ssum(j)+res

    call trapzd(to_int1,0._8,1._8,res)

            ssum(j)=ssum(j)+res


!-----------------------------------------------------------------------

    !calculate casimir (Drude)
    eps(1)=epsAurDr(i)
    eps(2)=epsAurDr(i)
    call trapzd(to_int,0._8,1._8,res)

            ssum1(j)=ssum1(j)+res

    call trapzd(to_int1,0._8,1._8,res)

            ssum1(j)=ssum1(j)+res

!----------------------------------------------------------------------

    !calculate casimir (Marachevsky)
    eps(1)=epsAurMar(i)
    eps(2)=epsAurMar(i)
    call trapzd(to_int,0._8,1._8,res)

            ssum2(j)=ssum2(j)+res

    call trapzd(to_int1,0._8,1._8,res)

            ssum2(j)=ssum2(j)+res

!----------------------------------------------------------------------

    !calculate casimir (Generalized plasma)
    eps(1)=epsAurGenPl(i)
    eps(2)=epsAurGenPl(i)
    call trapzd(to_int,1.0e-8_8,1._8,res)

            ssum3(j)=ssum3(j)+res

    call trapzd(to_int1,0._8,1._8,res)

            ssum3(j)=ssum3(j)+res

!---------------------------------------------------------------------

!calculate casimir (Lorentz-Drude)
    eps(1)=epsAurLD(i)
    eps(2)=epsAurLD(i)
    call trapzd(to_int,1.0e-8_8,1._8,res)

            ssum4(j)=ssum4(j)+res

    call trapzd(to_int1,0._8,1._8,res)

            ssum4(j)=ssum4(j)+res

!---------------------------------------------------------------------

!calculate casimir (Brendel-Bormann)
    eps(1)=epsAurBB(i)
    eps(2)=epsAurBB(i)
    call trapzd(to_int,1.0e-8_8,1._8,res)

            ssum5(j)=ssum5(j)+res

    call trapzd(to_int1,0._8,1._8,res)

            ssum5(j)=ssum5(j)+res

    end do
!*********
!ssum0 (Drude,Marachevsky,Kramers-Kr)
    ssum0=0._8
    !call trapzd(zero_sum_int_Dr, 0._8, 1._8, res)
    call qnc79(zero_sum_int_Dr, 0._8, 1._8, 1.0e-3_8, res, ierA, nevalA)

            ssum0=ssum0+res

    !call trapzd(zero_sum_int1_Dr, 0._8, 1._8, res)
    call qnc79(zero_sum_int1_Dr, 0._8, 1._8, 1.0e-3_8, res, ierA, nevalA)

            ssum0=ssum0+res

!ssum01 (Plasma)
   ssum01=0._8
   !call trapzd(zero_sum_int_Pl, 0._8, 1._8, res1)
   call qnc79(zero_sum_int_Pl, 0._8, 1._8, 1.0e-3_8, res1, ierA, nevalA)

            ssum01=ssum01+res1

   !call trapzd(zero_sum_int1_Pl, 0._8, 1._8, res1)
   call qnc79(zero_sum_int1_Pl, 0._8, 1._8, 1.0e-3_8, res1, ierA, nevalA)

            ssum01=ssum01+res1

!****************
! Kramers-Kr
fenergy(j)=1.602e-19*kb*T/(2*pi*dist**2)*(0.5*ssum0+ssum(j))
! Drude
fenergy1(j)=1.602e-19*kb*T/(2*pi*dist**2)*(0.5*ssum0+ssum1(j))
! Marachevsky
fenergy2(j)=1.602e-19*kb*T/(2*pi*dist**2)*(0.5*ssum0+ssum2(j))
! Generalized plasma
fenergy3(j)=1.602e-19*kb*T/(2*pi*dist**2)*(0.5*ssum01+ssum3(j))
! Drude-Lorentz
fenergy4(j)=1.602e-19*kb*T/(2*pi*dist**2)*(0.5*ssum0+ssum4(j))
! Brendel-Bormann
fenergy5(j)=1.602e-19*kb*T/(2*pi*dist**2)*(0.5*ssum0+ssum5(j))

rest0=1.602e-19*h*rest0/(2*pi)**2

ecasimir(j)=1.602e-19*(pi)**2*h*c/(720*dist**3)

!write data in casimirgold.txt
100 format(10(f15.3))
150 format(10(A15))

!values
write(11,100) dist*1.0e6, abs((fenergy(j))*1.0e9), abs((fenergy1(j))*1.0e9), &
abs((fenergy2(j))*1.0e9), abs((fenergy3(j))*1.0e9), abs((fenergy4(j))*1.0e9), &
abs((fenergy5(j))*1.0e9), abs(rest0)*1.0e9, ecasimir(j)*1.0e9


!want to write eta = free energy / casimir result
write(12,100) dist*1.0e6, abs(fenergy(j))/ecasimir(j), abs(fenergy1(j))/ecasimir(j), &
abs(fenergy2(j))/ecasimir(j), abs(fenergy3(j))/ecasimir(j), abs(fenergy4(j))/ecasimir(j), &
abs(fenergy5(j))/ecasimir(j), abs(rest0)/ecasimir(j)


!***********************************************************************************************

deallocate(epsAur)
deallocate(epsAurDr)
deallocate(epsAurMar)
deallocate(epsAurGenPl)
deallocate(epsAurLD)
deallocate(epsAurBB)

    end do


!**********************************************************************************************

deallocate(ssum)
deallocate(fenergy)
deallocate(ssum1)
deallocate(fenergy1)
deallocate(ssum2)
deallocate(fenergy2)
deallocate(ssum3)
deallocate(fenergy3)
deallocate(ssum4)
deallocate(fenergy4)
deallocate(ssum5)
deallocate(fenergy5)
deallocate(ecasimir)

close(11)
close(10)
close(12)

print*, 'Done!'

end program


	
