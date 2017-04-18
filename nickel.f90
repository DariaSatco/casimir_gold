module nickel

!Nickel parameters for Drude model
real(8), dimension(2), parameter:: Dr_Ni=(/ 4.89, 0.0436 /)
real(8), parameter:: miuNi=110.0D+00

contains

function drude_to_int_Ni(x)
!KK expression under integral for Drude model
use gold

real(8):: x, drude_to_int_Ni, param(2)

param=Dr_Ni
drude_to_int_Ni=param(1)**2*param(2)/((x**2+param(2)**2)*(x**2+freqA**2))

end function

!-----------------------------------------------------------------

function zero_sum_Dr_Ni(x)
!does not depend on distance if we have done change x --> dist*x
!Drude model
use dopcasimir

real(4):: zero_sum_Dr_Ni, x
real(4):: rTmn(2), rTen(2)

rTen(1)=0.0
rTen(2)=exp(-2*(t_Au+d)/dist*x)*(miuNi-1.0D+00)/(miuNi+1.0D+00)

rTmn(1)=1.0
rTmn(2)=1.0

if (x.eq.0.0) then
    zero_sum_Dr_Ni = 0.0
    else
    zero_sum_Dr_Ni = x*(log(1._8-rTmn(1)*rTmn(2)*exp(-2.*x))+log(1._8-rTen(1)*rTen(2)*exp(-2.*x)))
endif

end function

!----------------------------------------------------------------------------------

function zero_sum_Pl_Ni(x)
!does not depend on distance if we have done change x --> dist*x
!Drude model
use dopcasimir
use gold
use titanium

real(4):: zero_sum_Pl_Ni, x
real(4):: rTmn(2), rTen(2)
real(8):: wpG, wpNi, wpTi
real(8):: multipG, multipNi, multipTi
real(8):: expAu, expTi
real(8):: rTiNi, rAuTi, rAuTiNi

wpG=parAproxIm(1)
wpNi=Dr_Ni(1)
wpTi=Dr_Ti(1)

multipG=sqrt((x/dist)**2+(wpG/cEv)**2)
multipNi=miuNi*sqrt((x/dist)**2+(wpNi/cEv)**2)
multipTi=sqrt((x/dist)**2+(wpTi/cEv)**2)

rTen(1)=((x/dist)-multipG)/((x/dist)+multipG)

expAu=exp(-2*t_Au*multipG)
expTi=exp(-2*t_Au*multipNi)

rTiNi=(multipTi-multipNi)/(multipTi+multipNi)
rAuTi=(multipG-multipTi)/(multipG+multipTi)

rAuTiNi=(rAuTi+expTi*rTiNi)/(1.0+expTi*rAuTi*rTiNi)

rTen(2)=(rTen(1)+expAu*rAuTiNi)/(1.0+expAu*rTen(1)*rAuTiNi)

rTmn(1)=1.0
rTmn(2)=1.0

if (x.eq.0.0) then
    zero_sum_Pl_Ni = 0.0
    else
    zero_sum_Pl_Ni = x*(log(1._8-rTmn(1)*rTmn(2)*exp(-2.*x))+log(1._8-rTen(1)*rTen(2)*exp(-2.*x)))
endif

end function

!--------------------------------------------------------------------------------

function zero_sum_Dr_AuAu(x)
!does not depend on distance if we have done change x --> dist*x
!Drude model
use dopcasimir

real(4):: zero_sum_Dr_AuAu, x
real(4):: rTmn(2), rTen(2)

rTen(1)=0.0
rTen(2)=0.0

rTmn(1)=1.0
rTmn(2)=1.0

if (x.eq.0.0) then
    zero_sum_Dr_AuAu = 0.0
    else
    zero_sum_Dr_AuAu = x*(log(1._8-rTmn(1)*rTmn(2)*exp(-2.*x))+log(1._8-rTen(1)*rTen(2)*exp(-2.*x)))
endif

end function

!----------------------------------------------------------------------------------

function zero_sum_Pl_AuAu(x)
!does not depend on distance if we have done change x --> dist*x
!Drude model
use dopcasimir
use gold
use titanium

real(4):: zero_sum_Pl_AuAu, x
real(4):: rTmn(2), rTen(2)
real(8):: wpG, wpNi, wpTi
real(8):: multipG, multipNi, multipTi
real(8):: expAu, expTi
real(8):: rAuTi, rAuTiAu

wpG=parAproxIm(1)
wpTi=Dr_Ti(1)

multipG=sqrt((x/dist)**2+(wpG/cEv)**2)
multipTi=sqrt((x/dist)**2+(wpTi/cEv)**2)

rTen(1)=((x/dist)-multipG)/((x/dist)+multipG)

expAu=exp(-2*t_Au*multipG)
expTi=exp(-2*t_Au*multipNi)

rTiNi=(multipTi-multipNi)/(multipTi+multipNi)
rAuTi=(multipG-multipTi)/(multipG+multipTi)

rAuTiAu=(rAuTi+expTi*(-rAuTi))/(1.0-expTi*rAuTi**2)

rTen(2)=(rTen(1)+expAu*rAuTiAu)/(1.0+expAu*rTen(1)*rAuTiAu)

rTmn(1)=1.0
rTmn(2)=1.0

if (x.eq.0.0) then
    zero_sum_Pl_AuAu = 0.0
    else
    zero_sum_Pl_AuAu = x*(log(1._8-rTmn(1)*rTmn(2)*exp(-2.*x))+log(1._8-rTen(1)*rTen(2)*exp(-2.*x)))
endif

end function


end module
