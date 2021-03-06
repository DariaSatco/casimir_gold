module dopcasimir

!contains subroutines for integration and some functions for Lifshitz formula

    real(8), public :: w,dist
	! w - Matsubara frequency
	! dist - the distance between interacting plates
    real(8), public :: eps(3), epsNi, epsTi, epsAu, eps_gap
	! eps - dielectric permittivity vector
    ! eps(1),eps(2) describes materials
	! eps(3) describes gap permettivity
	real(8), public :: d,t_Au
	!d,t - thikness parameter of slabs (see Decca 2016)
    real(8), parameter :: pi=3.1415926
    real(8), parameter :: eV=1.519e15  ! rad/s
    real(8), parameter :: kb=8.617e-5  ! eV/K  Boltzman constant (kb=1.38e-23 J/K)
    real(8), parameter :: h=6.585e-16  ! eV*s  Plank constant (h=1.054e-34J*s)
    real(8), parameter :: c=299792458._8 ! c - velocity of light m/s
    real(8), parameter :: cEv=19.746e-8 ! eV*m

    integer, public :: model(2)

contains

!------------------------------------------------------------------------------------------------

	function eqf(a,x,y,rtm_func,rte_func)
	! function needed to integrate in formula for Casimir free energy 
	! x=dist*r change the variable to dimensioness
    ! the constants in front of the integral: rdr --> 1/dist**2 xdx
    ! w in eV
    ! c in eV*m (defined as cEv in module parameters)

	real(8):: a,x,q,y,rTMn(2),rTEn(2),eqf

    interface

    function rtm_func(kf, freq)
    real(8):: rtm_func(2)
    real(8):: kf, freq
    end function

    function rte_func(kf, freq)
    real(8):: rte_func(2)
    real(8):: kf, freq
    end function

    end interface

		rTMn=rtm_func(x/a,y)
!calculate the TM Frenel coefficient for the concrete dist as a function of x,y
		rTEn=rte_func(x/a,y)
!calculate the TE Frenel coefficient for the concrete dist as a function of x,y
			
		q=sqrt(x**2+(a*y/cEv)**2)
		eqf=x*(log(1._8-rTMn(1)*rTMn(2)*exp(-2.*q))+log(1._8-rTEn(1)*rTEn(2)*exp(-2.*q)))

					
	end function eqf

!-----------------------------------------------------------------------------------------------

    function epsmodel_func(x)
    use gold
    use silicon
    !epsmodel_func choose from the list of models desired one

    real(8):: epsmodel_func(2)
    real(8):: x

    do i=1,2

    select case(model(i))
        case (1)
        epsmodel_func(i)=epsMar(x)
        case (2)
        epsmodel_func(i)=Drude(x, parRakic_LD)
        case (3)
        epsmodel_func(i)=epsSil_art(x)
        case(4)
        epsmodel_func(i)=Gen_Plasma(x, parMost, gAu, gammaAu, wAu)
    end select

    end do

    end function


!-----------------------------------------------------------------------------------------------

    function zerotempint(arg)
	! function - argument of 2-dimendional integral by frequencies and k-frequencies at T=0 
        ! limits [1, Inf) --> (0,1]
        ! change variables xdx -->  1/(arg1**3)darg1;  dy --> 1/(arg2**2)darg2

	real(8):: zerotempint, arg(2)
	real(8):: q, k(2), a
	real(8):: rTMn(2), rTEn(2), epsmodel(2)

    a=dist
	q=sqrt(arg(1)**2+arg(2)**2)

    epsmodel=epsmodel_func(arg(2)*cEv/a)

	!need to write rTE(arg1,arg2), rTM(arg1,arg2)
	do i=1,2
	   k(i)=sqrt((arg(1)/a)**2+epsmodel(i)*(arg(2)/a)**2)
       rTEn(i)=(q-k(i))/(q+k(i))
       rTMn(i)=(epsmodel(i)*q-k(i))/(epsmodel(i)*q+k(i))
    end do

    if (arg(1)*arg(2) .eq. 0.0_8) then

        zerotempint=0.0_8
        else
	    zerotempint=arg(1)*(log(1._8-rTMn(1)*rTMn(2)*exp(-2.*q))+log(1._8-rTEn(1)*rTEn(2)*exp(-2.*q)))
    end if

	end function

!-----------------------------------------------------------------------------------------------
    function zerotempint1(x)

    real(8):: x(2), zerotempint1

    if (x(1)*x(2) .eq. 0._8) then

        zerotempint1=0.0_8
        else
        zerotempint1=(1._8/x(1))**2*(1._8/x(2))**2*zerotempint(x**(-1))

    end if

    end function
!-----------------------------------------------------------------------------------------------
	
	function to_int(arg)
	!integrable function for nonzero temperatures - function of one argument - x <--> k-frequency
    ! w in eV

	real(4):: arg, to_int

	to_int = eqf(dist, real(arg,8), w, rTM, rTE)

	end function


!----------------------------------------------------------------------------------

	function rTM(kf,freq)
	!the Frenel reflection coeffitient for TM-field
	! freq in eV
    ! kf in 1/m

	  real(8):: kf,freq, rTM(2),q,k(2)
	  integer:: i


		q=sqrt(kf**2+(freq/cEv)**2)

		do i=1,2
			k(i)=sqrt(kf**2+eps(i)*(freq/cEv)**2)
			rTM(i)=(eps(i)*q-k(i))/(eps(i)*q+k(i))
		end do

	end function rTM



!--------------------------------------------------------------------------------

	function rTE(kf,freq)
	! the Frenel reflection coeffitient for TE-field
	! freq in eV
	! kf in 1/m

	  real(8):: kf,freq,rTE(2),q,k(2)
	  integer:: i

	  
		q=sqrt(kf**2+(freq/cEv)**2)
		do i=1,2
			k(i)=sqrt(kf**2+eps(i)*(freq/cEv)**2)
			rTE(i)=(q-k(i))/(q+k(i))
		end do
	  
	end function rTE

!---------------------------------------------------------------------------------
!construct several function to calculete Casimir for three slab plate
!for more information see Decca (2016)

function k_prop(eps_prop , kf, freq)

real(8):: k_prop
real(8):: eps_prop, kf, freq

k_prop = sqrt(kf**2 + eps_prop*(freq/cEv)**2)

end function

!---------------------------------------------------------------------------------

function rTE_ab( eps_a, eps_b, kf, freq )

real(8):: rTE_ab
real(8):: eps_a, eps_b, kf, freq
real(8):: k_a, k_b

k_a = k_prop( eps_a, kf, freq )
k_b = k_prop( eps_b, kf, freq )

rTE_ab = ( k_a - k_b )/( k_a + k_b )

end function rTE_ab

!--------------------------------------------------------------------------------

function rTM_ab( eps_a, eps_b, kf, freq )

real(8):: rTM_ab
real(8):: eps_a, eps_b, kf, freq
real(8):: k_a, k_b

k_a = k_prop( eps_a, kf, freq )
k_b = k_prop( eps_b, kf, freq )

rTM_ab = ( eps_b*k_a - eps_a*k_b )/( eps_b*k_a + eps_a*k_b )

end function

!-------------------------------------------------------------------------------

function rTE_TM_123( eps_1, eps_2, eps_3, kf, freq, marker )
!in case Au-Ti-Ni slab - 1 <--> Au, 2 <--> Ti, 3 <--> Ni
!in case Au-Ti-Au slab - 1 <--> Au, 2 <--> Ti, 3 <--> Au
!marker = 1 ==> rTE
!marker = 2 ==> rTM

real(8):: rTE_123
real(8):: kf, freq
real(8):: eps_1, eps_2, eps_3
real(8):: num,denum
integer:: marker

if ( marker == 1 ) then
num = rTE_ab(eps_1, eps_2, kf, freq) + exp( -2 * d * k_prop (eps_2, kf, freq) )  * rTE_ab(eps_2, eps_3, kf, freq)
denum = 1.0D+00 + exp( -2 * d * k_prop (eps_2, kf, freq) ) * rTE_ab(eps_1, eps_2, kf, freq) * rTE_ab(eps_2, eps_3, kf, freq)

rTE_TM_123 = num/denum

elseif ( marker == 2 ) then

num = rTM_ab(eps_1, eps_2, kf, freq) + exp( -2 * d * k_prop (eps_2, kf, freq) )  * rTM_ab(eps_2, eps_3, kf, freq)
denum = 1.0D+00 + exp( -2 * d * k_prop (eps_2, kf, freq) ) * rTM_ab(eps_1, eps_2, kf, freq) * rTM_ab(eps_2, eps_3, kf, freq)

rTE_TM_123 = num/denum
else

print*, 'Incorrect marker, can be equal only to 1 (rTE) or 2 (rTM)'

endif

end function

!-----------------------------------------------------------------------------

function rTE_2_Ni( kf, freq )
!new Frenel TE coefficient for Ni
real(8):: rTE_2_Ni(2)
real(8):: kf, freq
real(8):: num, denum, q, k_Au

k_Au = k_prop(epsAu, kf, freq)
q=sqrt(kf**2+(freq/cEv)**2)

rTE_2_Ni(1)=(q-k_Au)/(q+k_Au)

num = rTE_ab(eps_gap, epsAu, kf, freq) + exp( -2*t_Au*k_Au )*rTE_TM_123(epsAu,epsTi,epsNi,kf,freq,1)
denum = 1.0D+00 + exp( -2*t_Au*k_Au) * rTE_ab(eps_gap, epsAu, kf, freq) * rTE_TM_123(epsAu,epsTi,epsNi,kf,freq,1)

rTE_2_Ni(2) = num/denum

end function

!-----------------------------------------------------------------------------

function rTM_2_Ni( kf, freq )
!new Frenel TM coefficient for Ni
real(8):: rTM_2_Ni(2)
real(8):: kf, freq
real(8):: num, denum, q, k_Au

k_Au = k_prop(epsAu, kf, freq)
q=sqrt(kf**2+(freq/cEv)**2)

rTM_2_Ni(1)=(epsAu*q-k_Au)/(epsAu*q+k_Au)

num = rTM_ab(eps_gap, epsAu, kf, freq) + exp( -2*t_Au*k_Au)*rTE_TM_123(epsAu,epsTi,epsNi,kf,freq,2)
denum = 1.0D+00 + exp( -2*t_Au*k_Au) * rTM_ab(eps_gap, epsAu, kf, freq) * rTE_TM_123(epsAu,epsTi,epsNi,kf,freq,2)

rTM_2_Ni(2) = num/denum

end function

!-----------------------------------------------------------------------------

function rTE_2_Au( kf, freq )
!new Frenel TE coefficient for Au
real(8):: rTE_2_Au(2)
real(8):: kf, freq
real(8):: num, denum, q, k_Au

k_Au = k_prop(epsAu, kf, freq)
q=sqrt(kf**2+(freq/cEv)**2)

rTE_2_Au(1)=(q-k_Au)/(q+k_Au)

num = rTE_ab(eps_gap, epsAu, kf, freq) + exp( -2*t_Au*k_Au)*rTE_TM_123(epsAu,epsTi,epsAu,kf,freq,1)
denum = 1.0D+00 + exp( -2*t_Au*k_Au) * rTE_ab(eps_gap, epsAu, kf, freq) * rTE_TM_123(epsAu,epsTi,epsAu,kf,freq,1)

rTE_2_Au(2) = num/denum

end function

!-----------------------------------------------------------------------------

function rTM_2_Au( kf, freq )
!new Frenel TM coefficient for Au
real(8):: rTM_2_Au(2)
real(8):: kf, freq
real(8):: num, denum, q, k_Au

k_Au = k_prop(epsAu, kf, freq)
q=sqrt(kf**2+(freq/cEv)**2)

rTM_2_Au(1)=(epsAu*q-k_Au)/(epsAu*q+k_Au)

num = rTM_ab(eps_gap, epsAu, kf, freq) + exp( -2*t_Au*k_Au)*rTE_TM_123(epsAu,epsTi,epsAu,kf,freq,2)
denum = 1.0D+00 + exp( -2*t_Au*k_Au) * rTM_ab(eps_gap, epsAu, kf, freq) * rTE_TM_123(epsAu,epsTi,epsAu,kf,freq,2)

rTM_2_Au(2) = num/denum

end function

!-----------------------------------------------------------------------------

function to_int_Au(arg)
    !integrable function for nonzero temperatures - function of one argument - x <--> k-frequency
    ! w in eV

    real(4):: arg, to_int_Au

    to_int_Au = eqf(dist, real(arg,8), w, rTM_2_Au, rTE_2_Au)

end function

!-----------------------------------------------------------------------------

function to_int_Ni(arg)
    !integrable function for nonzero temperatures - function of one argument - x <--> k-frequency
    ! w in eV

    real(4):: arg, to_int_Ni

    to_int_Ni = eqf(dist, real(arg,8), w, rTM_2_Ni, rTE_2_Ni)

end function

!-----------------------------------------------------------------------------
!*****************************
!INTEGRATION SUBROUTINES
!*****************************

!-------------------------------------------------------------------------------------------------	 
	 
	subroutine trapzd(func,lowl,upl,res)
	! approximates the integral with the trapezoidal method
        
    real(8):: func
	real(8):: res
	real(8):: lowl,upl
	real(8):: s,it
	integer:: i,j,n
	
	s=0._8
	res=0.5*abs(upl-lowl)*(func(lowl)+func(upl))
	i=0

		do while (abs(res-s).gt.(1e-3*abs(res)))
			s=res
			res=0._8
			i=i+1
			n=2**i
			it=abs(upl-lowl)/n
				
				do j=0,n-1
					res=res+0.5*it*(func(lowl+j*it)+func(lowl+(j+1)*it))
				enddo

		enddo

	!print *, 'step number=', 2**i
	return

	end subroutine

!------------------------------------------------------------------------------------------------

	subroutine trapzd2d(func,a,b,res,stepn)
	! approximates the integral with the trapezoidal method

	real(8):: func
	real(8):: res
	real(8), allocatable:: q(:,:)
	real(8):: a(2), b(2), x(2)
	real(8):: s,it(2)
	integer:: i,j,n,stepn
	
	s=0.
	res=0.25*abs(b(1)-a(1))*abs(b(2)-a(2))*(func(a)+func(b))
	i=0

		do while (abs(res-s).gt.(1e-2*abs(res)))
			s=res
			res=0._8
			i=i+1
			n=2**i
			it=abs(b-a)/n

			allocate(q(0:n-1,0:n-1))
				
				do j=0,n-1

				do k=0,n-1

				if ((j==0).or.(j==n-1)) then
					if ((k==0).or.(k==n-1)) then
					q(j,k)=0.25
						else
					q(i,j)=0.5
					endif
				else
					if ((k==0).or.(k==n-1)) then
					q(j,k)=0.5
						else
					q(j,k)=1.
					endif
				endif
				
				    x(1)=a(1)+j*it(1)
				    x(2)=a(2)+k*it(2)
					res=res+q(j,k)*(func(x))
				enddo

				enddo
			deallocate(q)

		res=res*it(1)*it(2)

		enddo
	
	stepn=n

	return

	end subroutine
	
!-----------------------------------------------------------------------------------------


end module dopcasimir
