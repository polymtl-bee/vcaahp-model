!TRNSYS Type3223: Variable capacity heat pump controller.
! ---------------------------------------------------------------------------------------------------------------------
!
! This routine implements a controller for an air-source heat pump with variable speed compressor.
!
!
! Inputs
! ---------------------------------------------------------------------------------------------------------------------
! Nb | Variable     | Description                                                   | Input Units     | Internal Units
! ---|--------------|---------------------------------------------------------------|-----------------|----------------
!  1 | Tset         | Setpoint temperature (command signal)                         | °C              | °C
!  2 | Tr           | Controlled variable                                           | °C              | °C
!  3 | onOff        | ON/OFF signal                                                 | -               | -
!  4 | fmin         | Minimum value for the frequency (control signal)              | -               | -
!  5 | fmax         | Maximum value for the frequency (control signal)              | -               | -
!  6 | Kc           | Gain constant                                                 | any             | any
!  7 | ti           | Integral time constant                                        | h               | h
!  8 | tt           | Tracking time constant                                        | h               | h
!  9 | b            | Proportional setpoint weight                                  | -               | -
!
!
! Outputs
! ---------------------------------------------------------------------------------------------------------------------
! Nb | Variable     | Description                                                   | Output  Units   | Internal Units
! ---|--------------|---------------------------------------------------------------|-----------------|----------------
!  1 | f            | Normalized frequency                                          | -               | -
!  2 | status       | Controller status                                             | -               | -


subroutine Type3223
!export this subroutine for its use in external DLLs
!dec$attributes dllexport :: Type3223

use, intrinsic :: iso_fortran_env, only : wp=>real64    ! Defines a constant "wp" (working precision) that can be used in real numbers, e.g. 1.0_wp, and sets it to real64 (double precision)

use TrnsysConstants
use TrnsysFunctions

implicit none

real(wp) :: time, timeStep  ! TRNSYS time and timestep
real(wp) :: Tset, Tr ! Temperatures
real(wp) :: f, fmin, fmax, fmean  ! Frequencies
real(wp) :: onOff, Kc, ti, tt, b  ! Controller parameters
real(wp) :: status  ! Controller status
real(wp) :: e, es, v, vp, vi  ! Controller signals
real(wp) :: h, N ! timestep, counter
real(wp) :: Tset1, Tr1, vi1, es1, e1  ! Values of the previous timestep

integer :: thisUnit, thisType    ! unit and type numbers

! Get the Global Trnsys Simulation Variables
time = getSimulationTime()
timestep = getSimulationTimeStep()
thisUnit = getCurrentUnit()
thisType = getCurrentType()

!--- Version signing call: set the version number for this Type --------------------------------------------------------

if (GetIsVersionSigningTime()) then

    call SetTypeVersion(18)
    return  ! We are done for this call

endif

! ---- Last call in the simulation (after last timestep has converged) -------------------------------------------------

if (GetIsLastCallofSimulation()) then

    ! This Type should not perform any task during the last call
    return  ! We are done for this call

endif

! --- End of timestep call (after convergence or too many iterations) --------------------------------------------------

if (GetIsEndOfTimestep()) then

    return  ! We are done for this call

endif

! --- Very first call of simulation (initialization call) --------------------------------------------------------------
if(getIsFirstCallofSimulation()) then

	!Tell the TRNSYS Engine How This Type Works
	call SetNumberofParameters(0)
	call SetNumberofInputs(9)
	call SetNumberofDerivatives(0)
	call SetNumberofOutputs(9)
	call SetIterationMode(1)  ! An indicator for the iteration mode (default=1).  Refer to section 8.4.3.5 of the documentation for more details.
	call SetNumberStoredVariables(0,0)  ! The number of static variables that the model wants stored in the global storage array and the number of dynamic variables that the model wants stored in the global storage array
	call SetNumberofDiscreteControls(0)  ! The number of discrete control functions set by this model (a value greater than zero requires the user to use Solver 1: Powell's method)

    h = getSimulationTimeStep()
    status = 2

	return

endIf

! --- Start time call: not a real time step, there are no iterations at the initial time - output initial conditions ---

if (getIsStartTime()) then


    Tset = GetInputValue(1)
    Tr = GetInputValue(2)
    onOff = GetInputValue(3)
    fmin = GetInputValue(4)
    fmax = GetInputValue(5)
    Kc = GetInputValue(6)
    ti = GetInputValue(7)
    tt = GetInputValue(8)
    b = GetInputValue(9)

!Check the Parameters for Problems (#,ErrorType,Text)
!Sample Code: If( PAR1 <= 0.) Call FoundBadParameter(1,'Fatal','The first parameter provided to this model is not acceptable.')

!Set the Initial Values of the Outputs (#,Value)
	call SetOutputValue(1, 0.0_wp) ! Normalized frequency
	call SetOutputValue(2, 0.0_wp) ! Controller status
    call SetOutputValue(3, status)
    call SetOutputValue(4, Tset)
    call SetOutputValue(5, Tr)


!If Needed, Set the Initial Values of the Static Storage Variables (#,Value)
!Sample Code: SetStaticArrayValue(1,0.d0)

!If Needed, Set the Initial Values of the Dynamic Storage Variables (#,Value)
!Sample Code: Call SetDynamicArrayValueThisIteration(1,20.d0)

!If Needed, Set the Initial Values of the Discrete Controllers (#,Value)
!Sample Code for Controller 1 Set to Off: Call SetDesiredDiscreteControlState(1,0)

	return

endIf

! --- TRNSYS has detected that parameters must be re-read - indicates another unit of this Type ------------------------
if(getIsReReadParameters()) then

endif

! --- Read inputs (for all calls except very first call in simulation) -------------------------------------------------

Tset = GetInputValue(1)
Tr = GetInputValue(2)
onOff = GetInputValue(3)
fmin = GetInputValue(4)
fmax = GetInputValue(5)
Kc = GetInputValue(6)
ti = GetInputValue(7)
tt = GetInputValue(8)
b = GetInputValue(9)

e = Tset - Tr  ! Error

if(ErrorFound()) return

! Default values for extra parameters
if (tt < 0.0_wp) then
   tt = ti
endif

if (b < 0.0_wp) then
   b = 1.0_wp
endif

! --- Recall stored values (...1 means at the end of previous time step) ---
status = getOutputValue(3)
Tset1 = getOutputValue(4)
Tr1 = getOutputValue(5)
fmean = getOutputValue(6)
N = getOutputValue(7)
vi1 = getOutputValue(8)
es1 = getOutputValue(9)
e1 = Tset1 - Tr1

if (onOff <= 0) then
    status = 0
    f = 0
elseif (abs(Tset - Tset1) > 5) then  ! Setpoint-change mode if setpoint step > 5�C
    status = 8
    f = fmax
elseif (status == 8) then
    f = fmax
    if (abs(e) < 1) then  ! Tr gets close to the setpoint: exit setpoint-change mode.
        status = 2
        vi = 0  ! Reset the integral
    endif
elseif (status == 2 .and. N*h > 0.2) then  ! After 0.2 h (12 min) in steady-state, keep the frequency constant.
    status = 4                                     !|--This value should be an input
    f = fmean
    N = 0
elseif (status == 4) then  ! steady-state mode: keep the frequency at its mean value.
    f = fmean
    if (abs(e) > 0.5) then  ! Count number of times the error is too high (change this condition, i.e. using variance)
        N = N + 1
    endif
    if (N*h > 0.08) then  ! Too much time with an error too big: switch back to normal mode.
        status = 2
        N = 0
    endif
elseif (status == 2) then  ! Normal mode: PI controller with anti-windup.
    vp = Kc * (b*Tset - Tr)  ! Proportional signal
    if (ti > 0.0_wp) then  ! Integral action
        vi = vi1 + Kc / ti * h * (e + e1) / 2  ! Update the integral (using trapezoidal integration).
        if (tt > 0.0_wp .and. (v < fmin .or. v > fmax)) then
            v = vp + vi  ! Unsaturated signal
            es = v - min(fmax, max(fmin, v))  ! Error with saturated signal
            vi = vi - h * (es + es1) / 2 / tt  ! De-saturate integral signal
        endif
    endif
    f = vp + vi
    if (abs(Tr - Tr1) < 1) then
        fmean = (f + N*fmean) / (N+1)  ! Update the frequency mean
        N = N + 1
    elseif (abs(Tr - Tr1) > 2) then
        fmean = 0
        N = 0
    endif
endif

!Check the Inputs for Problems (#,ErrorType,Text)
!Sample Code: If( IN1 <= 0.) Call FoundBadInput(1,'Fatal','The first input provided to this model is not acceptable.')



!Set the Outputs from this Model (#,Value)
call SetOutputValue(1, f) ! Normalized frequency
call SetOutputValue(2, status) ! Controller status
call SetOutputValue(3, status)
call SetOutputValue(4, Tset1)
call SetOutputValue(5, Tr1)
call SetOutputValue(6, fmean)
call SetOutputValue(7, N)
call SetOutputValue(8, vi1)
call SetOutputValue(9, es1)

!-----------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------
!If Needed, Store the Desired Disceret Control Signal Values for this Iteration (#,State)
!Sample Code:  Call SetDesiredDiscreteControlState(1,1)
!-----------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------
!If Needed, Store the Final value of the Dynamic Variables in the Global Storage Array (#,Value)
!Sample Code:  Call SetDynamicArrayValueThisIteration(1,T_FINAL_1)
!-----------------------------------------------------------------------------------------------------------------------

return
end
!-----------------------------------------------------------------------------------------------------------------------
