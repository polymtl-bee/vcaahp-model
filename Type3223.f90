! +---------------------------------------------------------+
! | TRNSYS Type3223: Variable capacity heat pump controller |
! +---------------------------------------------------------+
    
! This routine implements a controller for an air-source heat pump with variable speed compressor.


! Inputs
! --------------------------------------------------------------------------------------------------
!  # | Variable     | Description                                   | Input Units   | Internal Units
! --------------------------------------------------------------------------------------------------
!  1 | Tset         | Setpoint temperature (command signal)         | °C            | °C
!  2 | Tr           | Controlled variable                           | °C            | °C
!  3 | Toa          | Outdoor air temperature                       | °C            | °C
!  4 | onOff        | ON/OFF signal                                 | -             | -
!  5 | fmin         | Minimum value for the output frequency        | -             | -
!  6 | fmax         | Maximum value for the output frequency        | -             | -
!  7 | Kc           | Gain constant                                 | any           | any
!  8 | ti           | Integral time constant                        | h             | h
!  9 | tt           | Tracking time constant                        | h             | h
! 10 | b            | Proportional setpoint weight                  | -             | -
! 11 | N            | Number of frequency levels                    | -             | -
! 12 | mode         | 0 = cooling mode                              | -             | -
!                   | 1 = heating mode                              |               |
! 13 | defrost_mode |-1 = defrost cycles (normal behaviour)         | -             | -
!                   | 0 = defrost (off) mode                        |               |
!                   | 1 = recovery mode (transient)                 |               |
!                   | 2 = steady-state mode                         |               |
! --------------------------------------------------------------------------------------------------
    
! Parameters
! --------------------------------------------------------------------------------------------------
!  # | Variable     | Description                                   | Param. Units  | Internal Units
! --------------------------------------------------------------------------------------------------
!  1 | mode_deadband| Operating mode deadband for automatic mode    | °C            | °C
!  2 | LUcool       | Logical Unit - cooling mode                   | -             | -
!  3 | LUheat       | Logical Unit - heating mode                   | -             | -
! --------------------------------------------------------------------------------------------------

! Outputs
! --------------------------------------------------------------------------------------------------
!  # | Variable     | Description                                   | Output  Units | Internal Units
! --------------------------------------------------------------------------------------------------
!  1 | fq           | Normalized frequency (control signal)         | -             | -
!  2 | AFR          | Normalized air flow rate                      | -             | -
!  3 | mode         | 0 = cooling mode                              | -             | -
!                   | 1 = heating mode                              |               |
!  4 | defrost_mode | 0 = defrost (off) mode                        |               |
!                   | 1 = recovery mode (transient)                 |               |
!                   | 2 = steady-state mode                         |               |
! --------------------------------------------------------------------------------------------------
    
! Author: Gregor Strugala

module Type3223Data

use, intrinsic :: iso_fortran_env, only : wp=>real64    ! Defines a constant "wp" (working precision) that can be used in real numbers, e.g. 1.0_wp, and sets it to real64 (double precision)
implicit none

type Type3223DataStruct
    
    ! Frequency limitation parameters
    real(wp) :: e_min(0:1), e_max(0:1)
    integer :: nAFR(0:1), nAFRboost(0:1), nf0(0:1), nZones, nAFR2
    real(wp), allocatable :: AFR(:, :), AFRerror(:, :), AFRdb(:, :)
    real(wp), allocatable :: f0(:, :), Toa0(:, :)
    real(wp), allocatable :: f2(:, :), AFR2(:), Toa2(:), db2(:)
    real(wp) :: f2heat, t_boost_max, f1f2
    
    ! Defrost parameters
    real(wp) :: Tcutoff, t_cy, t_off, t_rec
    

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

real(wp) :: time, dt  ! TRNSYS time and timestep
real(wp) :: Tset, Tr, Toa ! Temperatures
real(wp) :: fsat, fq, fmin, fmax  ! Frequencies
real(wp) :: onOff, Kc, ti, tt, b  ! Controller parameters
real(wp) :: e, es, f, fp, fi  ! Controller signals
real(wp) :: h ! timestep
real(wp) :: Tset_old, Tr_old, fi_old, es_old, e_old  ! Values of the previous timestep
real(wp) :: mode_deadband  ! Parameters
integer :: LUheat, LUcool
integer :: N, mode, prev_mode  ! Number of frequency levels, operating mode
integer :: Ni = 1, Ninstances = 1  ! temporary, should use a kernel function to get the actual instance number.
integer :: AFRlevel, old_AFRlevel, zone, old_zone, AFR2level, Toalevel
real(wp) :: AFR, t_boost
logical :: modulate, fmaxBoost
integer :: defrost_mode
real(wp) :: t_ld, tau, t_uc, t_oc, recov_penalty

integer :: thisUnit, thisType    ! unit and type numbers

! Set the version number for this Type
if (GetIsVersionSigningTime()) then
    call SetTypeVersion(18)
    return
endif

call GetTRNSYSvariables()

! All the stuff that must be done once at the beginning
if(GetIsFirstCallofSimulation()) then
	call ExecuteFirstCallOfSimulation()
	return
endif

! Parameters must be re-read - indicates another unit of this Type
if(GetIsReReadParameters()) call ReadParameters()

! Start of the first timestep: no iterations, outputs initial conditions
if (GetIsStartTime()) then
    call ExecuteStartTime()
	return
endif

! End of timestep call (after convergence or too many iterations)
if (GetIsEndOfTimestep()) then
    call ExecuteEndOfTimestep()
    return
endif
    
if (GetIsLastCallofSimulation()) then
    call ExecuteLastCallOfSimulation()
    return
endif


call GetInputValues()
if (ErrorFound()) return

e = Tset - Tr
if (mode == -1) then
    prev_mode = int(GetDynamicArrayValueLastTimestep(4))
    if (e < -mode_deadband / 2.0_wp) then
        mode = 0
    else if (e > mode_deadband / 2.0_wp) then
        mode = 1
    else
        mode = prev_mode
    endif
end if
call SetDynamicArrayValueThisIteration(4, real(mode, wp))


! Get air flow rate
old_AFRlevel = int(GetDynamicArrayValueLastTimestep(5))
AFRlevel = GetLevel(s(Ni)%AFRerror(:, mode), s(Ni)%AFRdb(:, mode), e, old_AFRlevel, mode)
AFR = s(Ni)%AFR(AFRlevel, mode)
call SetDynamicArrayValueThisIteration(5, real(AFRlevel, wp))

! Get frequency limits
if (fmin < 0.0_wp) then
    Toalevel = FindLevel(s(Ni)%Toa0(:, mode), Toa, s(Ni)%nf0(mode) - 1)
    fmin = s(Ni)%f0(Toalevel, mode)
end if

if (getSimulationTime() > 3.44) then
    continue
end if

fmaxBoost = .false.
if (fmax < 0.0_wp) then
    if (mode == 0) then
        old_zone = int(GetDynamicArrayValueLastTimestep(6))
        t_boost = int(GetDynamicArrayValueLastTimestep(7))
        zone = GetLevel(s(Ni)%Toa2, s(Ni)%db2, Toa, old_zone, 1)
        call SetDynamicArrayValueThisIteration(6, real(zone, wp))
        AFR2level = FindLevel(s(Ni)%AFR2, AFR, s(Ni)%nAFR2)
        if (AFR2level > s(Ni)%nAFR2) AFR2level = s(Ni)%nAFR2
        if (t_boost < s(Ni)%t_boost_max) then
            fmax = s(Ni)%f2(zone, AFR2level)
        else
            fmax = s(Ni)%f2(zone, AFR2level) * s(Ni)%f1f2
        end if
        fmaxBoost = .true.
    else
        fmax = s(Ni)%f2heat
    end if
end if

! Assign fixed frequency value depending on error signal
modulate = .false.
if (e < s(Ni)%e_min(mode)) then
    fq = (1 - real(mode, wp)) * fmax
else if (e > s(Ni)%e_max(mode)) then
    fq = real(mode, wp) * fmax
else
    modulate = .true.
end if

call RecallStoredPIvalues()
if (modulate) then
    e = (Tset - Tr) * (2.0_wp * real(mode, wp) - 1.0_wp)  ! Error

    ! Default values for extra parameters
    if (tt < 0.0_wp) tt = ti
    if (b < 0.0_wp) b = 1.0_wp

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
        if (f > 0.0_wp) then
            fsat = min(fmax, max(fmin, f))  ! Saturated signal
            fq = (1.0_wp * floor(N * fsat)) / (1.0_wp * N)
        else
            fq = 0.0_wp
        end if
    endif
else if (ti > 0.0_wp) then
    fi = fi_old + Kc / ti * h * (e + e_old) / 2
end if
call StorePIvalues()

if (abs(fq - fmax) < 0.001_wp .and. fmaxBoost) then
    t_boost = t_boost + dt
else
    t_boost = 0.0_wp
endif
call SetDynamicArrayValueThisIteration(7, real(t_boost, wp))

if (defrost_mode == -1) call SetDefrostMode(defrost_mode)
recov_penalty = RecoveryPenalty(defrost_mode, t_ld, tau)

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
        if (ErrorFound()) return
        
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
            call SkipLines(LUheat, 7)
        read(LUh, *) s(Ni)%Tcutoff, s(Ni)%t_cy, s(Ni)%t_off, s(Ni)%t_rec
        close(LUh)
            call SkipLines(LUcool, 3)
        read(LUc, *) s(Ni)%t_boost_max, s(Ni)%f1f2
            call SkipLines(LUcool, 1)
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
    
    
    function GetLevel(centers, deadbands, value, old_level, hyst_dir)
        real(wp), intent(in) :: centers(:), deadbands(:), value
        integer, intent(in) :: old_level, hyst_dir
        real(wp) :: hdb(size(deadbands))
        integer :: level, GetLevel, idx
        hdb = deadbands / 2.0_wp  ! half deadbands
        level = FindLevel(centers, value, size(centers))
        idx = level
        if (level == 1) idx = 2
        if (level == size(centers) + 1) idx = size(centers)
        if (value < centers(idx-1) + hdb(idx-1)) then
            if (old_level < level .and. hyst_dir == 1) level = level - 1
            if (old_level < level .and. hyst_dir == 0) level = level + 1
        else if (value > centers(idx) - hdb(idx)) then
            if (old_level > level .and. hyst_dir == 1) level = level + 1
            if (old_level > level .and. hyst_dir == 0) level = level - 1
        end if
        GetLevel = level
    end function GetLevel
    
    function FindLevel(array, value, extent)
        real(wp), intent(in) :: array(:), value
        integer, intent(in) :: extent
        integer :: FindLevel
        integer :: L, R, mid
        if (value > array(extent)) then
            FindLevel = extent + 1
        else
            L = 1
            R = extent
            do while (L < R)
                mid = (L + R) / 2  ! L & R are integers -> automatic floor
                if (array(mid) < value) then
                    L = mid + 1
                else
                    R = mid
                end if
            end do
            FindLevel = L
        end if
    end function FindLevel
    
    
    subroutine SetDefrostMode(defrost_mode)
        integer :: defrost_mode
        real(wp) :: Tcutoff, t_cy, t_off, t_rec, t_ss
        t_cy = s(Ni)%t_cy
        t_off = s(Ni)%t_off
        t_rec = s(Ni)%t_rec
        t_ss = t_cy - t_off - t_rec
        tau = t_rec - t_off

        call RecallStoredDefrostValues()
        
        if (Toa < Tcutoff) then
            t_uc = t_uc + dt
        else
            t_oc = t_oc + dt
        end if
    
        if (t_ld < t_rec) then
            defrost_mode = 1
        else if (t_ld < t_rec + t_ss) then
            defrost_mode = 2
        else if (t_uc < t_oc) then
            t_uc = 0.0_wp
            t_oc = 0.0_wp
            t_ld = t_rec - dt
            defrost_mode = 2
        else if (t_ld < t_cy) then
            defrost_mode = 0
        else
            t_ld = -dt
            t_uc = 0.0_wp
            t_oc = 0.0_wp
            defrost_mode = 1
        end if
        t_ld = t_ld + dt
        
        call StoreDefrostValues()
    end subroutine SetDefrostMode
    
    
    subroutine StorePIvalues
        call SetDynamicArrayValueThisIteration(1, e)
        call SetDynamicArrayValueThisIteration(2, fi)
        call SetDynamicArrayValueThisIteration(3, es)
    end subroutine StorePIvalues
    
    
    subroutine RecallStoredPIvalues
        e_old = GetDynamicArrayValueLastTimestep(1)
        fi_old = GetDynamicArrayValueLastTimestep(2)
        es_old = GetDynamicArrayValueLastTimestep(3)
    end subroutine RecallStoredPIvalues
    
    
    subroutine StoreDefrostValues
        call SetDynamicArrayValueThisIteration(8, t_ld)
        call SetDynamicArrayValueThisIteration(9, t_uc)
        call SetDynamicArrayValueThisIteration(10, t_oc)
    end subroutine StoreDefrostValues
    
    
    subroutine RecallStoredDefrostValues
        t_ld = GetDynamicArrayValueLastTimestep(8)
        t_uc = GetDynamicArrayValueLastTimestep(9)
        t_oc = GetDynamicArrayValueLastTimestep(10)
    end subroutine RecallStoredDefrostValues
    
    function RecoveryPenalty(defrost_mode, t_ld, tau)
        integer, intent(in) :: defrost_mode
        real(wp), intent(in) :: t_ld, tau
        real(wp) :: recoveryPenalty
        if (defrost_mode == 1) then
            recoveryPenalty = 1.0_wp - exp(-1-t_ld/tau)
        else
            recoveryPenalty = 0.0_wp
        end if
    end function RecoveryPenalty
    
    subroutine GetInputValues
        Tset = GetInputValue(1)
        Tr = GetInputValue(2)
        Toa = GetInputValue(3)
        onOff = GetInputValue(4)
        fmin = GetInputValue(5)
        fmax = GetInputValue(6)
        Kc = GetInputValue(7)
        ti = GetInputValue(8)
        tt = GetInputValue(9)
        b = GetInputValue(10)
        N = GetInputValue(11)
        mode = GetInputValue(12)
        defrost_mode = GetInputValue(13)
    end subroutine GetInputValues
    
    
    subroutine ExecuteFirstCallOfSimulation
        call SetNumberofParameters(3)
	    call SetNumberofInputs(13)
	    call SetNumberofDerivatives(0)
	    call SetNumberofOutputs(5)
	    call SetIterationMode(1)
	    call SetNumberStoredVariables(0, 10)
	    call SetNumberofDiscreteControls(0)
        h = GetSimulationTimeStep()
        
        ! Allocate stored data structure
        if (.not. allocated(s)) then
            allocate(s(Ninstances))
        endif
        
        call ReadParameters()
        call ReadControlFiles(LUcool, LUheat)
    end subroutine ExecuteFirstCallOfSimulation
    
    
    subroutine ExecuteStartTime
        call ReadParameters()
        call GetInputValues()
	    call SetOutputValue(1, 0.0_wp)  ! Normalized frequency
        call SetOutputValue(2, 0.0_wp)  ! Normalized air flow rate
        call SetOutputValue(3, 0.0_wp)  ! Operating mode
        call SetOutputValue(4, 0.0_wp)  ! Defrost mode
        call SetOutputValue(5, 0.0_wp)  ! Defrost recovery penalty
        
        call SetDynamicArrayInitialValue(1, 0.0_wp)  ! error signal
        call SetDynamicArrayInitialValue(2, 0.0_wp)  ! integral value
        call SetDynamicArrayInitialValue(3, 0.0_wp)  ! saturation error
        call SetDynamicArrayInitialValue(4, 0.0_wp)  ! Operating mode
        call SetDynamicArrayInitialValue(5, 1.0_wp)  ! Air flow rate level
        call SetDynamicArrayInitialValue(6, 1.0_wp)  ! Outdoor temperature zone
        call SetDynamicArrayInitialValue(7, 0.0_wp)  ! Boost frequency operation time
        call SetDynamicArrayInitialValue(8, 0.0_wp)  ! Time since last defrost
        call SetDynamicArrayInitialValue(9, 0.0_wp)  ! Time under cutoff frequency
        call SetDynamicArrayInitialValue(10, 0.0_wp)  ! Time over cutoff frequency
    end subroutine ExecuteStartTime

    
    subroutine ExecuteEndOfTimestep
        continue
    end subroutine ExecuteEndOfTimestep
    
    
    subroutine ExecuteLastCallOfSimulation
        continue
    end subroutine ExecuteLastCallOfSimulation    
    
    
    subroutine ReadParameters
        mode_deadband = GetParameterValue(1)
        LUcool = GetParameterValue(2)
        LUheat = GetParameterValue(3)
    end subroutine ReadParameters
    
    
    subroutine SetOutputValues
        call SetOutputValue(1, fq)  ! Normalized saturated quantized frequency
        call SetOutputValue(2, AFR)  ! Normalized air flow rate
        call SetOutputValue(3, real(mode, wp))  ! Operating mode
        call SetOutputValue(4, real(defrost_mode, wp))  ! Defrost mode
        call SetOutputValue(5, recov_penalty)  ! Defrost recovery penalty
        return
    end subroutine SetOutputValues
    
    
    subroutine GetTRNSYSvariables
        time = GetSimulationTime()
        dt = GetSimulationTimeStep()
        thisUnit = GetCurrentUnit()
        thisType = GetCurrentType()
    end subroutine GetTRNSYSvariables
    
end subroutine Type3223
