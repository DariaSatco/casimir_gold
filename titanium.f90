module titanium

!Titanium parameters for Drude model
real(8), dimension(2), parameter:: Dr_Ti=(/ 2.42, 0.034 /)

contains

function drude_to_int_Ti(x)
!KK expression under integral for Drude model
use gold

real(8):: x, drude_to_int_Ti, param(2)

param=Dr_Ti
drude_to_int_Ti=param(1)**2*param(2)/((x**2+param(2)**2)*(x**2+freqA**2))

end function

end module
