! +---------------------------------------------------------+
! | TRNSYS Type3223: Variable capacity heat pump controller |
! +---------------------------------------------------------+
    
! This routine implements a controller for an air-source heat pump with variable speed compressor.


! Inputs
! ------------------------------------------------------------------------------------------------------
!  # | Variable     | Description                                       | Input Units   | Internal Units
! ------------------------------------------------------------------------------------------------------
!  1 | Tset         | Setpoint temperature (command signal)             | °C            | °C
!  2 | Tr           | Controlled variable                               | °C            | °C
!  3 | onOff        | ON/OFF signal                                     | -             | -
!  4 | fmin         | Minimum value for the frequency (control signal)  | -             | -
!  5 | fmax         | Maximum value for the frequency (control signal)  | -             | -
!  6 | Kc           | Gain constant                                     | any           | any
!  7 | ti           | Integral time constant                            | h             | h
!  8 | tt           | Tracking time constant                            | h             | h
!  9 | b            | Proportional setpoint weight                      | -             | -
! 10 | N            | Number of frequency levels                        | -             | -
! 11 | mode         | 0 = cooling mode                                  | -             | -
!                   | 1 = heating mode                                  |               |
! ------------------------------------------------------------------------------------------------------
    
! Parameters
! --------------------------------------------------------------------------------------------------
!  # | Variable     | Description                                   | Param. Units  | Internal Units
! --------------------------------------------------------------------------------------------------
!  1 | mode_deadband| 2 = Humidity ratio as humidity input          | °C            | °C
!  2 | Nosc_max     | Maximum number of oscillations                | -             | -
! --------------------------------------------------------------------------------------------------

! Outputs
! ------------------------------------------------------------------------------------------------------
!  # | Variable     | Description                                       | Output  Units | Internal Units
! ------------------------------------------------------------------------------------------------------
!  1 | f            | Normalized frequency                              | -             | -
!  2 | mode         | 0 = cooling mode                                  | -             | -
!                   | 1 = heating mode                                  |               |
! ------------------------------------------------------------------------------------------------------
    
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
real(wp) :: mode_deadband  ! Parameters
integer :: Nosc, Nosc_max
integer :: N, mode, prev_mode  ! Number of frequency levels, operating mode
integer :: Nsvar = 3, Noutputs = 2  ! number of of stored variables and outputs returned by the Type

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

e = Tset - Tr
if (mode == -1) then
    prev_mode = GetDynamicArrayValueLastTimestep(1)
    Nosc = GetDynamicArrayValueLastTimestep(2)
    if (e < -mode_deadband / 2.0_wp) then
        mode = 0
    else if (e > mode_deadband / 2.0_wp) then
        mode = 1
    else
        if (mode /= prev_mode) then
            Nosc = Nosc + 1
        else
            Nosc = 0
        endif
        if (Nosc > Nosc_max) then
            mode = 1
        else
            mode = prev_mode
        endif
    end if
    call SetDynamicArrayValueThisIteration(2, Nosc)
end if
call SetDynamicArrayValueThisIteration(1, mode)

e = (Tset - Tr) * (2.0_wp * real(mode, wp) - 1.0_wp)  ! Error

! Default values for extra parameters
if (tt < 0.0_wp) then
   tt = ti
endif

if (b < 0.0_wp) then
   b = 1.0_wp
endif

! Recall stored values
call RecallStoredValues()

if (onOff <= 0) then
    f = 0
else
    fp = Kc * (b*Tset - Tr) * (2.0_wp * real(mode, wp) - 1.0_wp)  ! Proportional signal
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
        call SetOutputValue(Noutputs+1, e)
        call SetOutputValue(Noutputs+2, fi)
        call SetOutputValue(Noutputs+3, es)
    end subroutine StoreValues
    
    
    subroutine RecallStoredValues
        e_old = GetOutputValue(Noutputs+1)
        fi_old = GetOutputValue(Noutputs+2)
        es_old = GetOutputValue(Noutputs+3)
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
        mode = GetInputValue(11)
    end subroutine GetInputValues
    

    subroutine ExecuteSpecialCases
    
    ! Last call in the simulation (after last timestep has converged)
    if (GetIsLastCallofSimulation()) then
        return
    endif

    ! End of timestep call (after convergence or too many iterations)
    if (GetIsEndOfTimestep()) then
        return
    endif

    ! Very first call of simulation (initialization call)
    if(GetIsFirstCallofSimulation()) then
	    call SetNumberofParameters(2)
	    call SetNumberofInputs(11)
	    call SetNumberofDerivatives(0)
	    call SetNumberofOutputs(Nsvar + Noutputs)
	    call SetIterationMode(1)
	    call SetNumberStoredVariables(0, 2)
	    call SetNumberofDiscreteControls(0)
        call SetIterationMode(2)
        h = GetSimulationTimeStep()
	    return
    endif

    if (GetIsStartTime()) then
        call ReadParameters()
        call GetInputValues()
	    call SetOutputValue(1, 0.0_wp)  ! Normalized frequency
        call SetOutputValue(2, 0.0_wp)  ! Operating mode
        
        call SetDynamicArrayInitialValue(1, 0)  ! Operating mode
        call SetDynamicArrayInitialValue(2, 5)  ! maximum oscillations number
	    return
    endif

    if(GetIsReReadParameters()) call ReadParameters()
    
    end subroutine ExecuteSpecialCases
    
    
    subroutine ReadParameters
        mode_deadband = jfix(GetParameterValue(1) + 0.01)
        Nosc_max = jfix(GetParameterValue(2) + 0.01)
    end subroutine ReadParameters
    
    
    subroutine SetOutputValues
        call SetOutputValue(1, fq)  ! Normalized saturated quantized frequency
        call SetOutputValue(2, real(mode, wp))  ! Operating mode
        return
    end subroutine SetOutputValues
    
    
    subroutine GetTRNSYSvariables
        time = getSimulationTime()
        timestep = getSimulationTimeStep()
        thisUnit = getCurrentUnit()
        thisType = getCurrentType()
    end subroutine GetTRNSYSvariables
    
end subroutine Type3223
