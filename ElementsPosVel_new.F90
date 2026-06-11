
!==================================================================================================
!    Перевод орбитальных элементов в положение и скорость...
!    Все переменные в угловых секундах!!!
!==================================================================================================
      
subroutine ElementsPosVel_new(Q, E, W, OM, I, M, A, AMASS, EPS, XV)

!==================================================================================================
!   Implicit none is your friend!
!==================================================================================================

implicit none

!==================================================================================================
!   Main parameters.
!==================================================================================================

real*8, intent(in) :: Q, E, W, OM, I, M, A
real*8, intent(out), dimension(:) :: XV(6)

real*8 :: AMASS, EPS
real*8 :: EX, U, V, UZ, AI, QQ, P, Z1, R, C0, C1, C2, C3, C4, C5, C6
real*8 :: XV1, XV2, XV3, XV4, XV5, XV6, SE, CE

real*8, parameter :: AK0 = 0.01720209895d0
real*8, parameter :: PI = 4.d0*atan(1.d0)
real*8, parameter :: RAD = 1.d0 !57.2957795130823d0 
real*8, parameter :: NEAR_PARAB = 1.d-10
real*8, parameter :: NEAR_PARAB_MIN = 1.d-14

real(8) :: AK

!==================================================================================================
!   ...
!==================================================================================================

AK = AK0

if(AMASS > 0.d0)then 
    AK = AK * dsqrt(1.d0 / AMASS)
endif

call KeplerEq_new(M / RAD, E, EX, 1.d-13)


if(dabs(1.d0 - E) < NEAR_PARAB .and. dabs(1.d0 - E) > NEAR_PARAB_MIN) then
    V = 2.d0 * datan(EX)
else
    if(E < 1.d0) then
        V = 2.d0 * datan(dsqrt((1.d0 + E) / (1.d0 - E)) * dtan(EX / 2.d0))
    else if(E == 1.d0) then
        V = 2.d0 * datan(EX)
    else
        V = 2.d0 * datan(dsqrt((E + 1.d0) / (E - 1.d0)) * dtanh(EX / 2.d0))  
    endif
endif

U = (W / RAD) + V
UZ = OM / RAD
AI = I / RAD

QQ = 1.d0 + E * dcos(V)
      
if(E == 1.d0) then
    P = 2.d0 * Q
else
    P = A * (1.D0 - E**2.d0)
endif

Z1 = AK / dsqrt(P)
R = P / QQ
C0 = E * dsin(V)

C1 = dcos(U) * dcos(UZ) - dsin(U) * dsin(UZ) * dcos(AI)
C2 = dcos(U) * dsin(UZ) + dsin(U) * dcos(UZ) * dcos(AI)
C3 = dsin(U) * dsin(AI)
C4 = dsin(U) * dcos(UZ) + DCOS(U) * dsin(UZ) * dcos(AI)
C5 = dsin(U) * dsin(UZ) - dcos(U) * dcos(UZ) * dcos(AI)
C6 = dcos(U) * dsin(AI)

XV1 = R * C1
XV2 = R * C2
XV3 = R * C3
XV4 = Z1 * (C0 * C1 - QQ * C4)
XV5 = Z1 * (C0 * C2 - QQ * C5)
XV6 = Z1 * (C0 * C3 + QQ * C6)
SE = dsin(EPS)
CE = dcos(EPS)

XV(1) = XV1
XV(2) = XV2 * CE - XV3 * SE
XV(3) = XV3 * CE + XV2 * SE
XV(4) = XV4
XV(5) = XV5 * CE - XV6 * SE
XV(6) = XV6 * CE + XV5 * SE

return
end

!================================================================================================== 
!    Решение уравнения Кеплера...
!================================================================================================== 

subroutine KeplerEq_new(AM, EXC, EX, PREC)
      
!==================================================================================================
!   Implicit none is your friend!
!==================================================================================================

implicit none

!==================================================================================================
!   Main parameters.
!==================================================================================================

real*8 :: AM, EXC, EX, PREC
real*8 :: EO, SIG, TET, U_TET, V_TET, EPS2

real*8, parameter :: NEAR_PARAB = 1.d-10
real*8, parameter :: NEAR_PARAB_MIN = 1.d-14

integer*4 :: inext
integer, parameter :: inext_lim = 10000
real*8 :: AK = 0.01720209895d0

!==================================================================================================
!   ...
!==================================================================================================

inext = 0

if(dabs(1.d0 - EXC) < NEAR_PARAB .and. dabs(1.d0 - EXC) > NEAR_PARAB_MIN) then

    EX = AM
    
    do
        inext = inext + 1
	    EO = EX
	    EX = EO + (AM - EO - EO**3 / 3.d0) / (1.d0 + EO**2)
	    if(dabs(EX - EO) < PREC) then
	        inext = 0
	        exit
	    endif
	    if(inext > inext_lim) then
	        inext = 0
	        exit
	    endif
    enddo

    SIG = EX
    TET = ((1.d0 - EXC) / 2.d0) * SIG**2
    EO = SIG
    
    inext = 0
    
    do
        U_TET = dsqrt(1.d0 - TET)
        V_TET = (1.d0 / 3.d0  +  (1.d0 / 10.d0) * TET + (3.d0 / 56.d0) * TET**2 + (5.d0 /144.d0) * TET**3 + (35.d0 /1408.d0) * TET**4 + (63.d0 /3328.d0) * TET**5 + (77.d0 /5120.d0) * TET**6)
        
        
        EX = EO + (AM - EO * U_TET - EO**3 * V_TET) * (AK * U_TET / ((1.d0 + EXC * EO**2) * dsqrt(2.d0)))
        
        
        TET = ((1.d0 - EXC) / 2.d0) * EX**2
        if(dabs(EX - EO) < PREC) then
            EPS2 = ((1.d0 - EXC) / 2.d0)
            EX = EX * dsqrt(1.d0 - EPS2) / dsqrt(1.d0 - TET)
            return
        else
            EO = EX
            inext = inext + 1
        endif    
        if(inext > inext_lim) then
            return
        endif
    enddo

else
    
    if(EXC < 1.d0) then

	    EX = AM
	    do
            inext = inext + 1
		    EO = EX
		    EX = EO - (AM + EXC * dsin(EO) - EO) / (EXC * dcos(EO) - 1.d0)
		    if(dabs(EX - EO) < PREC) return
		    if(inext > inext_lim) then
		        return
		    endif
	    enddo

    else if(EXC == 1.d0) then

        EX = AM
	    do
		    inext = inext + 1
		    EO = EX
		    EX = EO + (AM - EO - EO**3 / 3.d0) / (1.d0 + EO**2)
		    if(dabs(EX - EO) < PREC) return
		    if(inext > inext_lim) then
		        return
		    endif
	    enddo

    else

        EX = AM

	    do
		    inext = inext + 1
		    EO = EX
		    EX = EO - (AM - EXC * dsinh(EO) + EO) / (1.d0 - EXC * dcosh(EO))
		    if(dabs(EX - EO) < PREC) return
		    if(inext > inext_lim) then
		        return
		    endif
	    enddo
	    
    endif
    
endif
	
return
end