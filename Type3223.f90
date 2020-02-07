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
! 
! Author: Gregor Strugala

subroutine Type3223
!export this subroutine for its use in external DLLs
!dec$attributes dllexport :: Type3223

use, intrinsic :: iso_fortran_env, only : wp=>real64    ! Defines a constant "wp" (working precision) that can be used in real numbers, e.g. 1.0_wp, and sets it to real64 (double precision)

use TrnsysConstants
use TrnsysFunctions

implicit none

real(wp) :: time, timeStep  ! TRNSYS time and timestep
real(wp) :: Tset, Tr ! Temperatures
real(wp) :: fsat, fq, fmin, fmax  ! Frequencies
real(wp) :: onOff, Kc, ti, tt, b  ! Controller parameters
real(wp) :: e, es, f, fp, fi  ! Controller signals
real(wp) :: h ! timestep
real(wp) :: Tset_old, Tr_old, fi_old, es_old, e_old  ! Values of the previous timestep
integer :: N  ! Number of frequency levels
integer :: Nsvar = 4, Noutputs = 1  ! number of of stored variables and outputs returned by the Type

integer :: thisUnit, thisType    ! unit and type numbers

! Set the version number for this Type
if (GetIsVersionSigningTime()) then
    call SetTypeVersion(18)
    return  ! We are done for this call
endif

call GetTRNSYSvariables()
call ExecuteSpecialCases()
call GetInputValues()
if (ErrorFound()) return

e = Tset - Tr  ! Error

! Default values for extra parameters
if (tt < 0.0_wp) then
   tt = ti
endif

if (b < 0.0_wp) then
   b = 1.0_wp
endif

! Recall stored values
call RecallStoredValues()
e_old = Tset_old - Tr_old

if (onOff <= 0) then
    f = 0
else
    fp = Kc * (b*Tset - Tr)  ! Proportional signal
    if (ti > 0.0_wp) then  ! Integral action
        fi = fi_old + Kc / ti * h * (e + e_old) / 2  ! Update the integral (using trapezoidal integration).
    else
        fi = fi_old
    endif
    f = fp + fi  ! Unsaturated signal
    if (tt > 0.0_wp .and. (f < fmin .or. f > fmax)) then
        es = f - min(fmax, max(fmin, f))  ! Error with saturated signal
        fi = fi - h * (es + es_old) / 2 / tt  ! De-saturate integral signal
        f = fp + fi  ! Re-calculate the unsaturated signal
    endif
    fsat = min(fmax, max(fmin, f))  ! Saturated signal
    fq = (1.0_wp * floor(N * fsat)) / (1.0_wp * N)
endif

call StoreValues()
call SetOutputValues()

return

    contains
    
    subroutine StoreValues
    
    call SetOutputValue(Noutputs+1, Tset)
    call SetOutputValue(Noutputs+2, Tr)
    call SetOutputValue(Noutputs+3, fi)
    call SetOutputValue(Noutputs+4, es)
    
    end subroutine StoreValues
    
    
    subroutine RecallStoredValues
    
    Tset_old = GetOutputValue(Noutputs+1)
    Tr_old = GetOutputValue(Noutputs+2)
    fi_old = GetOutputValue(Noutputs+3)
    es_old = GetOutputValue(Noutputs+4)
    
    end subroutine RecallStoredValues
    
    
    subroutine GetInputValues
    
    Tset = GetInputValue(1)
    Tr = GetInputValue(2)
    onOff = GetInputValue(3)
    fmin = GetInputValue(4)
    fmax = GetInputValue(5)
    Kc = GetInputValue(6)
    ti = GetInputValue(7)
    tt = GetInputValue(8)
    b = GetInputValue(9)
    N = GetInputValue(10)
    
    end subroutine GetInputValues
    

    subroutine ExecuteSpecialCases
    
    ! Last call in the simulation (after last timestep has converged)
    if (GetIsLastCallofSimulation()) then
        ! This Type should not perform any task during the last call
        return  ! We are done for this call
    endif

    ! End of timestep call (after convergence or too many iterations)
    if (GetIsEndOfTimestep()) then
        return  ! We are done for this call
    endif

    ! Very first call of simulation (initialization call)
    if(GetIsFirstCallofSimulation()) then

	    !Tell the TRNSYS Engine How This Type Works
	    call SetNumberofParameters(0)
	    call SetNumberofInputs(10)
	    call SetNumberofDerivatives(0)
	    call SetNumberofOutputs(Nsvar + Noutputs)
	    call SetIterationMode(1)  ! An indicator for the iteration mode (default=1).  Refer to section 8.4.3.5 of the documentation for more details.
	    call SetNumberStoredVariables(0,0)  ! The number of static variables that the model wants stored in the global storage array and the number of dynamic variables that the model wants stored in the global storage array
	    call SetNumberofDiscreteControls(0)  ! The number of discrete control functions set by this model (a value greater than zero requires the user to use Solver 1: Powell's method)
        call SetIterationMode(2)
        h = GetSimulationTimeStep()

	    return

    endif

    ! Start time call: not a real time step, there are no iterations at the initial time - output initial conditions
    if (GetIsStartTime()) then

        call GetInputValues()

        !Check the Parameters for Problems (#,ErrorType,Text)
        !Sample Code: If( PAR1 <= 0.) Call FoundBadParameter(1,'Fatal','The first parameter provided to this model is not acceptable.')

        !Set the Initial Values of the Outputs (#,Value)
	    call SetOutputValue(1, 0.0_wp) ! Normalized frequency
        call SetOutputValue(2, Tset)
        call SetOutputValue(3, Tr)

	    return

    endif

    ! TRNSYS has detected that parameters must be re-read - indicates another unit of this Type
    if(GetIsReReadParameters()) then

    endif
    
    end subroutine ExecuteSpecialCases
    
    
    subroutine SetOutputValues
    
    !Set the Outputs from this Model (#,Value)
    call SetOutputValue(1, fq) ! Normalized saturated quantized frequency
    
    return
    
    end subroutine SetOutputValues
    
    subroutine GetTRNSYSvariables
        time = getSimulationTime()
        timestep = getSimulationTimeStep()
        thisUnit = getCurrentUnit()
        thisType = getCurrentType()
    end subroutine GetTRNSYSvariables
    
end subroutine Type3223
