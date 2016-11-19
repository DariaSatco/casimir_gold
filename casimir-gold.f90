program casimir_gold

use dopcasimir
use silicon
use gold

implicit none

!counters
integer:: i,j,k,l,s
!---------------------------------------
!casimir energy variables

integer:: N,stepn
real(8):: res, rest0, res1
real(8), allocatable:: fenergy(:), fenergy1(:), fenergy2(:), fenergy3(:), ecasimir(:)
real(8), allocatable:: ssum(:), ssum1(:), ssum2(:), ssum3(:)
real(8):: ssum0, ssum01, ssum02, ssum03
real(8):: T, wmax

!--------------------------------------------------------------------
!gold permettivity variables
integer,parameter:: ng=310 !rows number in file
!Palik data table
real(8):: matrixAu(ng,3)
!permettivity vector
real(8), allocatable:: epsAur(:), epsAurDr(:), epsAurMar(:), epsAurGenPl(:)
!integration results
real(8):: integralA1,integralA2,integralA3
!variables needed fo Kr-Kr integral in cspint subroutine
real(8)::funvalAu(ng),SpcoefA(3,ng), exintA(ng), workA(ng)
!variables for Kr-K integral in qags subroutine
real(4):: abserrA,abserrA2
integer(4)::nevalA,ierA,nevalA2,ierA2, ierr

real(8):: a(2), b(2)


!-----------------------------------------------------------------------
open(unit=11, file='casimirgold.txt', status='replace')
!casimirsilicon.txt, casimirgold.txt - files with computational results

!-----------------------------------------------------------------------
!code for calulation of casimir free energy

open(unit=10, file='resultAu.txt', status='old')
!read matrix from file (gold)
!1st column - frequencies
!2nd column - real part of permittivity
!3rd column - imaginary part of permittivity

do i=1,ng
read(10,*), (matrixAu(ng+1-i,j),j=1,3)
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
    wmax=stepn*c/dist
    w=2*pi*kb*T/h
    do while (w.le.wmax)
    k=k+1
    w=2*k*pi*kb*T/h
    end do

print*, 'k=', k
!k=100

allocate(epsAur(k))
allocate(epsAurDr(k))
allocate(epsAurMar(k))
allocate(epsAurGenPl(k))

        !calculate casimir energy
        do i=1,k
        !do-cycle for frequencies

            w=2*i*pi*kb*T/h
            !w - Matsubara frequency
            freqA=w

    !integrates the Drude-model approximation
    call qnc79(drude_to_int, 0.0_8, matrixAu(1,1), 1.0e-3_8, integralA1, ierA, nevalA)

    !builds vector to integrate from the Palik data
        do k=1,ng
        funvalAu(k)=matrixAu(k,3)*matrixAu(k,1)/(matrixAu(k,1)**2+freqA**2)
        end do
    !integrates the vector of data
    call cspint(ng, matrixAu(1:ng,1), funvalAu, matrixAu(1,1), matrixAu(ng,1), SpcoefA, exintA, workA, integralA2)

    !integrates the oscillator-function
    call monte_carlo(oscil_to_int, matrixAu(ng,1), matrixAu(ng,1)*1.0e4_8, 1000000, integralA3)


    !Kr-Kr formula for permittivity
    epsAur(i)=(integralA1+integralA2+integralA3)*2/pi+1

    epsAurDr(i)=Drude(w,parAproxIm)

    epsAurMar(i)=epsMar(w)

    epsAurGenPl(i)=Gen_Plasma(w)

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

     end do
!*********
!ssum0 (Drude,Marachevsky,Kramers-Kr)
    ssum0=0._8
    call trapzd(zero_sum_int_Dr,0._8,1._8, res)
    !call qnc79(zero_sum_int_Dr, 0._8, 1._8, 1.0e-3_8, res, ierr, s)

            ssum0=ssum0+res

    call trapzd(zero_sum_int1_Dr, 0._8,1._8, res)
    !call qnc79(zero_sum_int1_Dr, 0._8, 1._8, 1.0e-3_8, res, ierr, s)

            ssum0=ssum0+res
            print*, 'ssum0=', ssum0

 !ssum01 (Plasma)
   ssum01=0._8
   call trapzd(zero_sum_int_Pl, 0.1_8, 1._8, res1)
   print*, res1
   !call qnc79(zero_sum_int_Pl, 0._8, 1._8, 1.0e-3_8, res, ierr, s)

            ssum01=ssum01+res1

    call trapzd(zero_sum_int1_Pl, 0.1_8,1._8, res1)
    print*, res1
    !call qnc79(zero_sum_int1_Pl, 0._8, 1._8, 1.0e-3_8, res, ierr, s)

            ssum01=ssum01+res1
            print*, 'ssum01=', ssum01

!****************
! Kramers-Kr
fenergy(j)=1.602e-19*kb*T/(2*pi*dist**2)*(0.5*ssum0+ssum(j))
! Drude
fenergy1(j)=1.602e-19*kb*T/(2*pi*dist**2)*(0.5*ssum0+ssum1(j))
! Marachevsky
fenergy2(j)=1.602e-19*kb*T/(2*pi*dist**2)*(0.5*ssum0+ssum2(j))
! Generalized plasma
fenergy3(j)=1.602e-19*kb*T/(2*pi*dist**2)*(0.5*ssum01+ssum3(j))

rest0=1.602e-19*h*rest0/(2*pi)**2

ecasimir(j)=1.602e-19*(pi)**2*h*c/(720*dist**3)

write(11,*) dist*1.0e7, abs((fenergy(j))*1.0e9), abs((fenergy1(j))*1.0e9), &
abs((fenergy2(j))*1.0e9), abs((fenergy3(j))*1.0e9), abs(rest0)*1.0e9, ecasimir(j)*1.0e9
!*********

deallocate(epsAur)
deallocate(epsAurDr)
deallocate(epsAurMar)
deallocate(epsAurGenPl)

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

close(11)
close(10)

end program


	
