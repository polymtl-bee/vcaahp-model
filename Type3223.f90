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

real(wp) :: time, timeStep    ! TRNSYS time and timestep
real(wp) :: Tset, Tr ! Temperatures
real(wp) :: f, fmin, fmax ! Frequencies
real(wp) :: onOff, Kc, ti, tt, b ! Controller parameters
real(wp) :: status ! Controller status

! Local variables
real(wp) :: e, vp  ! Controller signals
real(wp) :: h  ! timestep

integer :: thisUnit, thisType    ! unit and type numbers

! Get the Global Trnsys Simulation Variables
time = getSimulationTime()
timestep = getSimulationTimeStep()
thisUnit = getCurrentUnit()
thisType = getCurrentType()

!--- Version signing call: set the version number for this Type --------------------------------------------------------

if (GetIsVersionSigningTime()) then

    Call SetTypeVersion(18)
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
	call SetNumberofOutputs(2)
	call SetIterationMode(1)  ! An indicator for the iteration mode (default=1).  Refer to section 8.4.3.5 of the documentation for more details.
	call SetNumberStoredVariables(0,0)  ! The number of static variables that the model wants stored in the global storage array and the number of dynamic variables that the model wants stored in the global storage array
	call SetNumberofDiscreteControls(0)  ! The number of discrete control functions set by this model (a value greater than zero requires the user to use Solver 1: Powell's method)
    
    h = getSimulationTimeStep()
    
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
	call SetOutputValue(2, 0.0_w) ! Controller status


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
Tset1 = getOutputValue(8)
Tr1 = getOutputValue(9)
u1 = getOutputValue(10)
v1 = getOutputValue(11)
vp1 = getOutputValue(12)
vi1 = getOutputValue(13)
vd1 = getOutputValue(14)
e1 = Tset1 - Tr1

vp = Kc * (b*Tset - Tr)  ! Proportional signal
if (ti > 0.0_wp) then  ! Integral action
    if (tt > 0.0_wp) then
        vi = vi1 + h/tt * (u1 - v1)
    endif
endif

!Check the Inputs for Problems (#,ErrorType,Text)
!Sample Code: If( IN1 <= 0.) Call FoundBadInput(1,'Fatal','The first input provided to this model is not acceptable.')



!Set the Outputs from this Model (#,Value)
call SetOutputValue(1, 0) ! Normalized frequency
call SetOutputValue(2, 0) ! Controller status

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

