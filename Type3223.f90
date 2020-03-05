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

module Type3223Data

use, intrinsic :: iso_fortran_env, only : wp=>real64    ! Defines a constant "wp" (working precision) that can be used in real numbers, e.g. 1.0_wp, and sets it to real64 (double precision)
use TrnsysConstants
use TrnsysFunctions
implicit none

type Type3223DataStruct
    
    ! Parameters
    real(wp) :: e_min(0:1), e_max(0:1)
    integer :: nAFR(0:1), nAFRboost(0:1), nf0(0:1), nZones, nAFR2
    real(wp), allocatable :: AFR(:, :), AFRerror(:, :), AFRdb(:, :)
    real(wp), allocatable :: f0(:, :), Toa0(:, :)
    real(wp), allocatable :: f2(:, :), AFR2(:), Toa2(:), db2(:)
    real(wp) :: f2heat

end type Type3223DataStruct

type(Type3223DataStruct), allocatable, save :: s(:)

end module Type3223Data

subroutine Type3223
!export this subroutine for its use in external DLLs
!dec$attributes dllexport :: Type3223

use, intrinsic :: iso_fortran_env, only : wp=>real64    ! Defines a constant "wp" (working precision) that can be used in real numbers, e.g. 1.0_wp, and sets it to real64 (double precision)

use TrnsysConstants
use TrnsysFunctions
use Type3223Data

implicit none

real(wp) :: time, timeStep  ! TRNSYS time and timestep
real(wp) :: Tset, Tr ! Temperatures
real(wp) :: fsat, fq, fmin, fmax  ! Frequencies
real(wp) :: onOff, Kc, ti, tt, b  ! Controller parameters
real(wp) :: e, es, f, fp, fi  ! Controller signals
real(wp) :: h ! timestep
real(wp) :: Tset_old, Tr_old, fi_old, es_old, e_old  ! Values of the previous timestep
real(wp) :: mode_deadband  ! Parameters
integer :: Nosc, Nosc_max, LUheat, LUcool
integer :: N, mode, prev_mode  ! Number of frequency levels, operating mode
integer :: Nsvar = 3, Noutputs = 2  ! number of of stored variables and outputs returned by the Type
integer :: Ni = 1, Ninstances = 1  ! temporary, should use a kernel function to get the actual instance number.

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
    if (e < 0) then
        mode = 0
    else
        mode = 1
    end if
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
    
    subroutine ReadControlFiles(LUc, LUh)
        character (len=maxPathLength) :: cfCoolPath
        character (len=maxPathLength) :: cfHeatPath
        integer, intent(in) :: LUc, LUh
        integer :: nAFRmax, nf0max, i, j, LUs(2), LUcool(1), LUheat(1)
        LUcool(1) = LUc
        LUheat(1) = LUh
    
        ! Ni = GetCurrentUnit()
        LUs = (/LUc, LUh/)
    
        cfCoolPath = GetLUfileName(LUc)
        cfHeatPath = GetLUfileName(LUh)
        call CheckControlFile(cfCoolPath)
        call CheckControlFile(cfHeatPath)
        
        open(LUh, file=cfHeatPath, status='old')
        open(LUc, file=cfCoolPath, status='old')
    
            call SkipLines(LUs, 6)
        read(LUc, *) s(Ni)%e_min(0)
        read(LUh, *) s(Ni)%e_min(1)
            call SkipLines(LUs, 2)
        read(LUc, *) s(Ni)%e_max(0)
        read(LUh, *) s(Ni)%e_max(1)
            call SkipLines(LUs, 3)
        read(LUc, *) s(Ni)%nAFR(0)
        read(LUh, *) s(Ni)%nAFR(1)
        nAFRmax = maxval(s(Ni)%nAFR)
            call SkipLines(LUs, 1)
        allocate(s(Ni)%AFR(nAFRmax, 0:1))
        read(LUc, *) (s(Ni)%AFR(i, 0), i = 1, s(Ni)%nAFR(0))
        read(LUh, *) (s(Ni)%AFR(i, 1), i = 1, s(Ni)%nAFR(1))
            call SkipLines(LUs, 4)
        allocate(s(Ni)%AFRerror(nAFRmax - 1, 0:1))
        allocate(s(Ni)%AFRdb(nAFRmax - 1, 0:1))
        do i = 1, s(Ni)%nAFR(0) - 1
            read(LUc, *) s(Ni)%AFRerror(i, 0), s(Ni)%AFRdb(i, 0)
            read(LUh, *) s(Ni)%AFRerror(i, 1), s(Ni)%AFRdb(i, 1)
        enddo
            call SkipLines(LUs, 3)
        read(LUc, *) s(Ni)%nf0(0)
        read(LUh, *) s(Ni)%nf0(1)
        nf0max = maxval(s(Ni)%nf0)
            call SkipLines(LUs, 2)
        allocate(s(Ni)%Toa0(nf0max - 1, 0:1))
        allocate(s(Ni)%f0(nf0max, 0:1))
        read(LUc, *) (s(Ni)%Toa0(i, 0), i = 1, s(Ni)%nf0(0) - 1)
        read(LUh, *) (s(Ni)%Toa0(i, 1), i = 1, s(Ni)%nf0(1) - 1)
        read(LUc, *) (s(Ni)%f0(i, 0), i = 1, s(Ni)%nf0(0))
        read(LUh, *) (s(Ni)%f0(i, 1), i = 1, s(Ni)%nf0(1))
            call SkipLines(LUheat, 2)
        read(LUh, *) s(Ni)%f2heat
        close(LUh)
            call SkipLines(LUcool, 3)
        read(LUc, *) s(Ni)%nZones, s(Ni)%nAFR2
            call SkipLines(LUcool, 6)
        allocate(s(Ni)%Toa2(s(Ni)%nZones))
        allocate(s(Ni)%db2(s(Ni)%nZones))
        allocate(s(Ni)%AFR2(s(Ni)%nAFR2))
        allocate(s(Ni)%f2(s(Ni)%nZones + 1, s(Ni)%nAFR2))
        do i = 1, s(Ni)%nZones
            read(LUc, *) s(Ni)%Toa2(i), s(Ni)%db2(i)
        end do
            call SkipLines(LUcool, 1)
        read(LUc, *) (s(Ni)%AFR2(j), j = 1, s(Ni)%nAFR2)
            call SkipLines(LUcool, 1)
        do i = 1, s(Ni)%nZones + 1
            read(LUc, *) (s(Ni)%f2(i, j), j = 1, s(Ni)%nAFR2)
        end do
        
        close(LUc)
        
    end subroutine ReadControlFiles
    
    
    subroutine CheckControlFile(cfPath)
        logical :: ControlFileFound = .false.
        character (len=maxPathLength) :: cfPath
        character (len=maxMessageLength) :: msg
        inquire(file=trim(cfPath), exist=ControlFileFound)
        if ( .not. ControlFileFound ) then
            write(msg,'("""",a,"""")') trim(cfPath)
            msg = "Could not find the specified performance map file. Searched for: " // trim(msg)
            call Messages(-1, msg, 'fatal', thisUnit, thisType)
            return
        end if
    end subroutine CheckControlFile
    
    
    subroutine SkipLines(LUs, N)
        integer, intent(in) :: LUs(:)
        integer :: i, j, N
        do i = 1, size(LUs)
            do j = 1, N
                read(LUs(i), *)
            end do
        end do
    end subroutine SkipLines
    
    
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
    
    ! All the stuff that must be done once at the beginning
    if(GetIsFirstCallofSimulation()) then
	    call SetNumberofParameters(4)
	    call SetNumberofInputs(11)
	    call SetNumberofDerivatives(0)
	    call SetNumberofOutputs(Nsvar + Noutputs)
	    call SetIterationMode(1)
	    call SetNumberStoredVariables(0, 2)
	    call SetNumberofDiscreteControls(0)
        call SetIterationMode(2)
        h = GetSimulationTimeStep()
        
        ! Allocate stored data structure
        if (.not. allocated(s)) then
            allocate(s(Ninstances))
        endif
        
        call ReadParameters()
        call ReadControlFiles(LUcool, LUheat)
        
	    return
    endif
    
    ! Start of the first timestep: no iterations, outputs initial conditions
    if (GetIsStartTime()) then
        call ReadParameters()
        call GetInputValues()
	    call SetOutputValue(1, 0.0_wp)  ! Normalized frequency
        call SetOutputValue(2, 0.0_wp)  ! Operating mode
        
        call SetDynamicArrayInitialValue(1, 0)  ! Operating mode
        call SetDynamicArrayInitialValue(2, 5)  ! maximum oscillations number
	    return
    endif
    
    ! Parameters must be re-read - indicates another unit of this Type
    if(GetIsReReadParameters()) call ReadParameters()
    
    ! End of timestep call (after convergence or too many iterations)
    if (GetIsEndOfTimestep()) then
        return  ! We are done for this call
    endif
    
    if (GetIsLastCallofSimulation()) then
        return  ! We are done for this call
    endif
    
    end subroutine ExecuteSpecialCases
    
    
    subroutine ReadParameters
        mode_deadband = jfix(GetParameterValue(1) + 0.01)
        Nosc_max = jfix(GetParameterValue(2) + 0.01)
        LUcool = GetParameterValue(3)
        LUheat = GetParameterValue(4)
    end subroutine ReadParameters
    
    
    subroutine SetOutputValues
        call SetOutputValue(1, fq)  ! Normalized saturated quantized frequency
        call SetOutputValue(2, real(mode, wp))  ! Operating mode
        return
    end subroutine SetOutputValues
    
    
    subroutine GetTRNSYSvariables
        time = GetSimulationTime()
        timestep = GetSimulationTimeStep()
        thisUnit = GetCurrentUnit()
        thisType = GetCurrentType()
    end subroutine GetTRNSYSvariables
    
end subroutine Type3223
