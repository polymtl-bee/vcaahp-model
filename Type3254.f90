! +-----------------------------------------------------------------------------+
! | TRNSYS Type3254: Variable capacity air-air heat pump with performance files |
! +-----------------------------------------------------------------------------+

! This routine implements an air-air heat pump with variable speed compressor.


! Inputs
! --------------------------------------------------------------------------------------------------
!  # | Variable     | Description                                   | Input  Units  | Internal Units
! --------------------------------------------------------------------------------------------------
!  1 | Tr           | Inlet (return) air temperature                | °C            | °C
!  2 | wr           | Inlet (return) air humidity ratio             | -             | -
!  3 | RHr          | Inlet (return) air relative humidity          | % (base 100)  | -
!  4 | pr           | Inlet (return) air pressure                   | atm           | atm
!  5 | mDot         | Inlet (return) air mass flow rate             | kg/h          | kg/h
!  6 | Toa          | Outdoor air dry bulb temperature              | °C            | °C
!  7 | woa          | Outdoor air humidity ratio                    | -             | -
!  8 | RHoa         | Outdoor air realtive humidity                 | % (base 100)  | -
!  9 | poa          | Outdoor air pressure                          | atm           | atm
! 10 | freq         | Normalized frequency signal                   | -             | -
! 11 | AFR          | Inlet (return) normalized air flow rate       | -             | -
! 12 | mode         | 0 = cooling mode                              | -             | -
!                   | 1 = heating mode                              |               |
! 13 | defrost_mode | 0 = defrost (off) mode                        | -             | -
!                   | 1 = recovery mode (transient)                 |               |
!                   | 2 = steady-state mode                         |               |
! 14 | recov_penalty| Penalty factor for defrost recovery mode      | -             | -
! 15 | PfanI        | Indoor fan power                              | kJ/h          | kJ/h
! 16 | PfanO        | Outdoor fan power                             | kJ/h          | kJ/h
! --------------------------------------------------------------------------------------------------

! Parameters
! --------------------------------------------------------------------------------------------------
!  # | Variable     | Description                                   | Param. Units  | Internal Units
! --------------------------------------------------------------------------------------------------
!  1 | psymode      | 2 = Humidity ratio as humidity input          | -             | -
!                   | 4 = Relative humidity as humidity input       |               |
!  2 | PelcRated    | Rated total cooling power                     | kJ/h          | kJ/h
!  3 | QcRated      | Rated cooling capacity                        | kJ/h          | kJ/h
!  4 | PelhRated    | Rated total heating power                     | kJ/h          | kJ/h
!  5 | QhRated      | Rated heating capacity                        | kJ/h          | kJ/h
!  6 | AFRrated     | Rated inlet air mass flow rate                | kg/h          | kg/h
!  7 | freqRatedHeat| Rated heating frequency                       | 1/s           | 1/s
!  8 | freqRatedCool| Rated cooling frequency                       | 1/s           | 1/s
!  9 | LUcool       | Logical Unit - cooling mode                   | -             | -
! 10 | LUheat       | Logical Unit - heating mode                   | -             | -
! --------------------------------------------------------------------------------------------------

! Outputs
! --------------------------------------------------------------------------------------------------
!  # | Variable     | Description                                   | Output  Units | Internal Units
! --------------------------------------------------------------------------------------------------
!  1 | Ts           | Outlet (supply) air temperature               | °C            | °C
!  2 | ws           | Outlet (supply) air humidity ratio            | -             | -
!  3 | RHs          | Outlet (supply) air % RH                      | % (base 100)  | % (base 100)
!  4 | ps           | Outlet (supply) air pressure                  | atm           | atm
!  5 | mDot         | Outlet (supply) air mass flow rate            | kg/h          | kg/h
!  6 | Qc           | Total cooling rate                            | kJ/h          | kJ/h
!  7 | Qcs          | Sensible cooling rate                         | kJ/h          | kJ/h
!  8 | Qcl          | Latent cooling rate                           | kJ/h          | kJ/h
!  9 | Qrej         | Heat rejection rate                           | kJ/h          | kJ/h
! 10 | Qh           | Total heating rate                            | kJ/h          | kJ/h
! 11 | Qabs         | Heat absorption rate                          | kJ/h          | kJ/h
! 12 | Pel          | Total power consumption                       | kJ/h          | kJ/h
! 13 | COP          | Coefficient of performance                    | -             | -
! 14 | EER          | Energy efficiency rating                      | -             | -
! 15 | PfanI        | Indoor fan power                              | kJ/h          | kJ/h
! 16 | PfanO        | Outdoor fan power                             | kJ/h          | kJ/h
! 17 | Pcomp        | Compressor power                              | kJ/h          | kJ/h
! 18 | fComp        | Compressor frequency                          | 1/s           | 1/s
! 19 | Tc           | Condensate temperature                        | °C            | °C
! 20 | cmfr         | Condensate mass flow rate                     | kg/h          | kg/h
! 21 | defrost_mode | 0 = defrost (off) mode                        | -             | -
!                   | 1 = Recovery mode (transient)                 |               |
!                   | 2 = Steady-state mode                         |               |
! 22 | shutdown     | 0 = No (false)                                | -             | -
!                   | 1 = Yes (true)                                |               |
! --------------------------------------------------------------------------------------------------

module Type3254Data

use, intrinsic :: iso_fortran_env, only : wp=>real64    ! Defines a constant "wp" (working precision) that can be used in real numbers, e.g. 1.0_wp, and sets it to real64 (double precision)
implicit none

type Type3254DataStruct

    ! Parameters
    real(wp), allocatable :: entries(:, :, :), lbounds(:, :), ubounds(:, :)
    integer, allocatable :: extents(:, :)
    integer :: PMClength, PMHlength  ! Length of the flattened performance map

    ! Performance matrices
    real(wp), allocatable :: PelcMap(:, :, :, :, :)
    real(wp), allocatable :: QcsMap(:, :, :, :, :)
    real(wp), allocatable :: QclMap(:, :, :, :, :)
    real(wp), allocatable :: PelhMap(:, :, :, :)
    real(wp), allocatable :: QhMap(:, :, :, :)

end type Type3254DataStruct

type(Type3254DataStruct), allocatable, save :: s(:)

end module Type3254Data


subroutine Type3254
!export this subroutine for its use in external DLLs
!DEC$Attributes DLLexport :: Type3254

use, intrinsic :: iso_fortran_env, only : wp=>real64    ! Defines a constant "wp" (working precision) that can be used in real numbers, e.g. 1.0_wp, and sets it to real64 (double precision)

use TrnsysConstants
use TrnsysFunctions
use Type3254Data

implicit none

integer :: thisUnit, thisType  ! unit and type numbers
real(wp) :: time, dt  ! TRNSYS time and timestep

! Proforma variables
real(wp) :: Tr, wr, RHr, mDot, AFR, pr, Toa, woa, RHoa, poa, freq, PfanI, PfanO  ! Inputs
integer :: mode, defrost_mode = 1  ! assume that the heat pump starts up at the beginning -> start with transient state
integer :: psymode, LUcool, LUheat  ! Parameters
real(wp) :: PelcRated, QcRated, PelhRated, QhRated, AFRrated, freqRatedh, freqRatedc  ! Parameters (rated values)
real(wp) :: Ts, ws, RHs, ps  ! Outputs (supply conditions)
real(wp) :: Pel, Qc, Qcs, Qcl, Qrej, Qh, Qabs, Pcomp  ! Outputs (heat and power)
real(wp) :: fComp, COP, EER, Tc, cmfr, recov_penalty  ! Outputs (misc)


! Local variables
real(wp) :: psydat(9), Twbr, Twboa, hr, hx, hs, dr
integer :: status
integer, parameter :: Ninstances = 1  ! Number of units (should be provided by a function)
integer :: Ni = 1  ! temporary, should use a kernel function to get the actual instance number.
real(wp) :: defrost_corr(2) = 0.0_wp  ! defrost correction factors
logical :: shutdown

! Performance map reading variables
integer, parameter :: Nc = 5, Nh = 4  ! Number of interpolation variables
integer, parameter :: Nmax = max(Nc, Nh)
integer :: i, j, N, Noutc = 3, Nouth = 2, Nout


! Interpolation variables
real(wp), allocatable :: interpolationResults(:), point(:)

! Set the version number for this Type
if (GetIsVersionSigningTime()) then
    call SetTypeVersion(18)
    return
endif

time = GetSimulationTime()
dt = GetSimulationTimeStep()
thisUnit = GetCurrentUnit()
thisType = GetCurrentType()

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
fComp = freq * ((1-mode) * freqRatedc + mode * freqRatedh)
shutdown = .false.
if (ErrorFound()) return

! Ni = GetCurrentUnit()
N = (1-mode) * Nc + mode * Nh  ! number of interpolation variables, depending on the operating mode.
Nout = (1-mode) * Noutc + mode * Nouth  ! number of ouput values for the interpolation

! Determine return air state
psydat(1) = pr
psydat(2) = Tr
psydat(4) = RHr/100.0_wp
psydat(6) = wr
if (mode==0) then
    call MoistAirProperties(thisUnit, thisType, 1, psymode, 1, psydat, 1, status)
    ! (unit, type, si units used, psych inputs, Twb computed, inputs, warning mgmt, warning occurences)
else
    call MoistAirProperties(thisUnit, thisType, 1, psymode, 0, psydat, 1, status)  ! Twb not computed
end if
pr = psydat(1)
Tr = psydat(2)
Twbr = psydat(3)
RHr = psydat(4)  ! RHr between 0 and 1 (not 0 and 100)
wr = psydat(6)
hr = psydat(7)
dr = psydat(9)
if (AFR >= 0) then
    mDot = AFR * AFRrated * dr  ! use normalized AFR as input if it is positive
else
    AFR = mDot / (dr * AFRrated)
endif

if (mode == 1) then ! compute outdoor air wet bulb for the interpolation
    psydat(1) = poa
    psydat(2) = Toa
    psydat(4) = RHoa/100.0_wp
    psydat(6) = woa
    call MoistAirProperties(thisUnit, thisType, 1, psymode, 1, psydat, 1, status)
    Twboa = psydat(3)
    defrost_corr = Correction(defrost_mode, recov_penalty)
end if

if (freq > 0) then
    ! Interpolate using wet bulb in cooling
    allocate(point(N))
    point(1) = Tr
    if (mode == 0) then
        point(2) = Twbr
        point(3) = Toa
        point(4) = AFR
        point(5) = freq
    else
        point(2) = Toa
        point(3) = AFR
        point(4) = freq
    end if
    allocate(interpolationResults(Nout))
    interpolationResults = Interpolate(point, mode, Nout)
    deallocate(point)
    if (mode == 0) Pel = interpolationResults(1) * PelcRated
    if (mode == 1) Pel = interpolationResults(1) * PelhRated * defrost_corr(1)
    Qcs = interpolationResults(2) * QcRated
    Qcl = interpolationResults(3) * QcRated
    Qh = interpolationResults(4) * QhRated * defrost_corr(2)
    Qc = Qcs + Qcl
    deallocate(interpolationResults)
else
    Qcs = 0.0_wp
    Qcl = 0.0_wp
    Qh = 0.0_wp
    Qc = 0.0_wp
    Pel = 0.0_wp
endif


! Determine supply air state

ps = pr  ! Fan pressure drop neglected

! Moist air state
if (Qc < Qcs) then
    Qc = Qcs
    ! Add warning
endif
if (mDot /= 0.0_wp) then
    ws = wr
    if (mode == 0) then
        hs = hr - Qc/mDot
        hx = hr  ! useful when the following if clause is not true
        if (Qcl > 0.0_wp) then  ! compute humidity after condensation
            psydat(1) = pr
            psydat(2) = Tr
            hx = hr - Qcl/mDot
            psydat(7) = hx  ! enthalpy of the state (Tr, ws)
            call MoistAirProperties(thisUnit, thisType, 1, 5, 0, psydat, 1, status)  ! dry-bulb and enthalpy as inputs
            if (ErrorFound()) return
            ws = psydat(6)
        endif
    else
        hs = hr + Qh/mDot
    end if
else
    hs = hr
    ws = wr
    hx = hr
endif

psydat(1) = ps
psydat(6) = ws
psydat(7) = hs
call MoistAirProperties(thisUnit, thisType, 1, 7, 0, psydat, 1, status)  ! humidity ratio and enthalpy as inputs
if (ErrorFound()) return
ps = psydat(1)
Ts = psydat(2)
RHs = psydat(4)
ws = psydat(6)
hs = psydat(7)

! Re-calculate heat transfer whose value is modified if saturation occurs
if (mode == 0 .and. freq > 0.0_wp) then
    Qcs = mDot * (hx - hs)  ! Sensible cooling rate
    Qcl = mDot * (hr - hx)  ! Latent cooling rate
    Qc = Qcs + Qcl  ! Total cooling rate
    Qrej = Qc + Pel  ! Heat rejection in ambient air
    Qabs = 0.0_wp
else if (freq > 0.0_wp) then
    Qrej = 0.0_wp
    Qabs = Qh - Pel  ! Heat absorption in ambient air
end if
if (freq > 0.0_wp) Pcomp = Pel - PfanI - PfanO  ! Compressor power

! COP, EER and condensate
if (Pel /= 0.0_wp) then
    COP = (Qc + Qh) / Pel
else
    COP = 0.0_wp
endif
EER = 3.413_wp * COP
Tc = Ts
cmfr = mDot * (wr - ws)  ! Condensate flow rate - water balance

if (shutdown) fcomp = 0.0_wp  ! force shutdown -> set output frequency to zero

call SetOutputValues()

return

    contains

    subroutine ReadPermap(LUc, LUh)
        integer, intent(in) :: LUc, LUh
        character (len=maxPathLength) :: permapCoolPath
        character (len=maxPathLength) :: permapHeatPath
        integer :: i, j, LUs(2), LUcool(1), LUheat(1)
        integer :: nTr, nTwbr, nToa, nAFR, nfreq  ! number of entries for each variable
        real(wp) :: filler(Nmax)
        LUcool(1) = LUc
        LUheat(1) = LUh

        ! Ni = GetCurrentUnit()
        LUs = (/LUc, LUh/)
        
        permapCoolPath = GetLUfileName(LUc)
        permapHeatPath = GetLUfileName(LUh)
        call CheckPMfile(permapCoolPath)
        call CheckPMfile(permapHeatPath)
        if (ErrorFound()) return
        
        open(LUc, file=permapCoolPath, status='old')
        open(LUh, file=permapHeatPath, status='old')

            call SkipLines(LUs, 6)
        allocate(s(Ni)%extents(Nmax, 0:1))
        allocate(s(Ni)%lbounds(Nmax, 0:1))
        allocate(s(Ni)%ubounds(Nmax, 0:1))
        do i = 1, Nc
            call SkipLines(LUcool, 1)
            read(LUc, *) s(Ni)%extents(i, 0), s(Ni)%lbounds(i, 0), s(Ni)%ubounds(i, 0)
        end do
        do i = 1, Nh
            call SkipLines(LUheat, 1)
            read(LUh, *) s(Ni)%extents(i, 1), s(Ni)%lbounds(i, 1), s(Ni)%ubounds(i, 1)
        end do
        s(Ni)%PMClength = product(s(Ni)%extents(1:Nc, 0))
        s(Ni)%PMHlength = product(s(Ni)%extents(1:Nh, 1))
        allocate(s(Ni)%entries(maxval(s(Ni)%extents), Nmax, 0:1))
        do i = 1, Nc
            call SkipLines(LUcool, 1)
            read(LUc, *) (s(Ni)%entries(j, i, 0), j = 1, s(Ni)%extents(i, 0))
        end do
        do i = 1, Nh
            call SkipLines(LUheat, 1)
            read(LUh, *) (s(Ni)%entries(j, i, 1), j = 1, s(Ni)%extents(i, 1))
        end do
            call SkipLines(LUs, 4)

        nTr = s(Ni)%extents(1, 0)
        nTwbr = s(Ni)%extents(2, 0)
        nToa = s(Ni)%extents(3, 0)
        nAFR = s(Ni)%extents(4, 0)
        nfreq = s(Ni)%extents(5, 0)
        allocate(s(Ni)%PelcMap(nTr, nTwbr, nToa, nAFR, nfreq))
        allocate(s(Ni)%QcsMap(nTr, nTwbr, nToa, nAFR, nfreq))
        allocate(s(Ni)%QclMap(nTr, nTwbr, nToa, nAFR, nfreq))
        do i = 1, s(Ni)%PMClength
            read(LUc, *) (filler(j), j = 1, Nc), Pel, Qcs, Qcl
            call SetPMvalue(s(Ni)%PelcMap, RowToColMajorOrder(i, s(Ni)%extents(1:Nc, 0)), Pel, s(Ni)%PMClength)
            call SetPMvalue(s(Ni)%QcsMap, RowToColMajorOrder(i, s(Ni)%extents(1:Nc, 0)), Qcs, s(Ni)%PMClength)
            call SetPMvalue(s(Ni)%QclMap, RowToColMajorOrder(i, s(Ni)%extents(1:Nc, 0)), Qcl, s(Ni)%PMClength)
        end do

        close(LUc)

        nTr = s(Ni)%extents(1, 1)
        nToa = s(Ni)%extents(2, 1)
        nAFR = s(Ni)%extents(3, 1)
        nfreq = s(Ni)%extents(4, 1)
        allocate(s(Ni)%PelhMap(nTr, nToa, nAFR, nfreq))
        allocate(s(Ni)%QhMap(nTr, nToa, nAFR, nfreq))
        do i = 1, s(Ni)%PMHlength
            read(LUh, *) (filler(j), j = 1, Nh), Pel, Qh
            call SetPMvalue(s(Ni)%PelhMap, RowToColMajorOrder(i, s(Ni)%extents(1:Nh, 1)), Pel, s(Ni)%PMHlength)
            call SetPMvalue(s(Ni)%QhMap, RowToColMajorOrder(i, s(Ni)%extents(1:Nh, 1)), Qh, s(Ni)%PMHlength)
        end do

        close(LUh)

    end subroutine ReadPermap


    subroutine CheckPMfile(permapPath)
        logical :: permapFileFound = .false.
        character (len=maxPathLength) :: permapPath
        character (len=maxMessageLength) :: msg
        inquire(file=trim(permapPath), exist=permapFileFound)
        if ( .not. permapFileFound ) then
            write(msg,'("""",a,"""")') trim(permapPath)
            msg = "Could not find the specified performance map file. Searched for: " // trim(msg)
            call Messages(-1, msg, 'fatal', thisUnit, thisType)
            return
        end if
    end subroutine CheckPMfile
    
    
    subroutine SkipLines(LUs, N)
        integer, intent(in) :: LUs(:)
        integer :: i, j, N
        do i = 1, size(LUs)
            do j = 1, N
                read(LUs(i), *)
            end do
        end do
    end subroutine SkipLines
    
    
    function RowToColMajorOrder(rowIndex, extents) result(colIndex)
    ! RowToColMajorOrder transforms a row-major order index
    ! corresponding to a given array shape into a column-major index.
    !
    ! Inputs
    !   rowIndex (integer) : one-based index of an array in row-major order.
    !   extents (integer array) : shape of the array that rowIndex is indexing.
    !
    ! Output
    !   colIndex (integer) : one-based index corresponding to the array element
    !                        indexed by rowIndex, but with a column-major order.
        integer, intent(in) :: extents(:), rowIndex
        integer :: i, colIndex, j, p
            
        i = rowIndex - 1  ! rowIndex is one-based
        colIndex = 1  ! colIndex is one-based
        p = product(extents)
        do j = size(extents), 1, -1
            p = p / extents(j)  ! ("/" performs integer division)
            colIndex = colIndex + p * modulo(i, extents(j))
            i = i / extents(j)
        end do
    end function RowToColMajorOrder


    function Interpolate(point, mode, Nout) result(interpolation)
    ! Interpolate performs a piecewise multilinear interpolation on the tables in
    ! the structure Type3254DataStruct.
    !
    ! Inputs
    !   point (real(wp) array) : point where the interpolation has to be performed.
    !   mode (integer) : operating mode, cooling (0) or heating(1)
    !   Nout (integer) : number of tables for the interpolation,
    !                    also the number of required output values.
    !
    ! Output
    !   interpolation (real(wp) array) : results of the interpolation.
        real(wp), intent(in) :: point(N)
        real(wp) :: scaled_point(N), left, right, sp
        real(wp), allocatable :: hypercube(:, :)
        integer, intent(in) :: mode
        integer, dimension(N) :: idx, lb_idx, counter_int, ones, zeros
        integer :: i
        logical :: counter_bool(N), increasing(N), in_range
        integer, intent(in) :: Nout
        real(wp) :: interpolation(Noutc + Nouth - 1)
        zeros = 0
        ones = 1
        
        ! Check if point is within the operating range before interpolation
        if (mode == 0) increasing = (/.true., .false., .false., .true., .true./)
        if (mode == 1) increasing = (/.false., .true., .true., .true./)
        in_range = .true.
        shutdown = .false.
        do i = 1, N
            if (increasing(i) .and. (point(i) < s(Ni)%lbounds(i, mode)) .or. &
            (.not. increasing(i)) .and. (point(i) > s(Ni)%ubounds(i, mode))) then
                in_range = .false.
                shutdown = .true.
                interpolation = 0.0_wp
                exit
            end if
        end do
        
        if (in_range) then
            do i = 1, N
                if (increasing(i) .and. (point(i) > s(Ni)%ubounds(i, mode)) .or. &
                (.not. increasing(i)) .and. (point(i) < s(Ni)%lbounds(i, mode))) then
                    in_range = .false.
                    exit
                end if
            end do
        end if
        
        
        if (.not. shutdown) then
            ! Map the point to the unit N-hypercube with the table values bounding the point
            do i = 1, N
                ! Find the index of the lower bound for the ith component of point
                j = Findlb(s(Ni)%entries(:, i, mode), point(i), s(Ni)%extents(i, mode))
                lb_idx(i) = j
                left = s(Ni)%entries(j, i, mode)
                right = s(Ni)%entries(j+1, i, mode)
                ! Compute the scaled point, and keep it between 0 and 1
                scaled_point(i) = min(1.0_wp, max(0.0_wp, (point(i) - left) / (right - left)))
            end do
            allocate(hypercube(Nout, 2**N))
            counter_bool = .true.  ! All N values of the binary counter are initialized at 1
            do i = 1, 2**N
                call Increment(counter_bool)
                counter_int = merge(ones, zeros, counter_bool)  ! convert logical to int
                idx = lb_idx + counter_int  ! index of each of the 2**N vertices
                hypercube(:, i) = Vertex(idx, mode)  ! Get table values associated with the vertices
                if (any(hypercube(:, i) < 0.0_wp)) then
                    deallocate(hypercube)
                    interpolation = 0.0_wp
                    shutdown = .true.
                    exit
                end if
            end do
        end if
        
        if (.not. shutdown) then
            ! Interpolate using the scaled point and the table values
            do i = 1, N
                sp = scaled_point(i)
                j = N - i
                hypercube(:, :2**j) = (1-sp) * hypercube(:, :2**j) &
                                        + sp * hypercube(:, 2**j+1:2**(j+1))
            end do
            if (mode == 0) then  ! cooling -> fill the first Noutc values, the rest are zeros
                do i = 1, Noutc
                    interpolation(i) = hypercube(i, 1)
                end do
                do i = Noutc+1, Noutc+Nouth-1
                    interpolation(i) = 0.0_wp
                end do
            else  ! heating -> first value is still the input power, values 2 through Noutc are zeros
                interpolation(1) = hypercube(1, 1)
                do i = 2, Noutc
                    interpolation(i) = 0.0_wp
                end do
                do i = Noutc+1, Noutc+Nouth-1
                    interpolation(i) = hypercube(i-Noutc+1, 1)
                end do
            end if
            deallocate(hypercube)
        end if
    end function Interpolate


    function Vertex(idx, mode)
    ! Vertex provides the performance map values associated with a certain index,
    ! taking the operating mode into account.
    !
    ! Inputs
    !   idx (integer array) : index of the table values to be retrieved.
    !   mode (integer) : operating mode, cooling (0) or heating (1).
    !
    ! Output
    !   vertex (real(wp) array) : array with the table values, given in the form
    !                             (power, sensible capacity, latent capacity) in
    !                             cooling and (power, total capacity) in heating.
        integer, intent(in) :: idx(:), mode
        real(wp) :: Pel, Qcs, Qcl, Qh, vertex(Nout)
        if (mode == 0) then
            Pel = GetPMvalue(mode, s(Ni)%PelcMap, idx)
            Qcs = GetPMvalue(mode, s(Ni)%QcsMap, idx)
            Qcl = GetPMvalue(mode, s(Ni)%QclMap, idx)
            vertex = (/Pel, Qcs, Qcl/)
        else
            Pel = GetPMvalue(mode, s(Ni)%PelhMap, idx)
            Qh = GetPMvalue(mode, s(Ni)%QhMap, idx)
            vertex = (/Pel, Qh/)
        end if
    end function Vertex


    function GetPMvalue(mode, array, idx) result(value)
    ! GetPMvalue retrieves the value at a given index in a performance map.
    !
    ! Inputs
    !   mode (integer, 0 or 1) : mode (heating or cooling) associated with the performance map.
    !   array : performance map, contained in the structure Type3254DataStruct.
    !   idx (integer array with the same length as the array shape) : index of the element to retrieve.
    !
    ! Output
    !   value : element retireved at index idx.
        integer, intent(in) :: idx(:)
        integer :: mode, i, array_idx
        ! array has a length equal to the number of elements in the appropriate performance map
        real(wp) :: array((1-mode)*s(Ni)%PMClength + mode*s(Ni)%PMHlength)
        real(wp) :: value
        array_idx = idx(1)  ! begin with the leftmost index, since arrays are stored in column-major order
        do i = 2, N
            array_idx = array_idx + product(s(Ni)%extents(1:i-1, mode)) * (idx(i) - 1)
        end do
        value = array(array_idx)
    end function GetPMvalue


    subroutine SetPMvalue(array, idx, value, PMlength)
    ! SetPMvalue sets a value in a performance map at a given index.
    !
    ! Inputs
    !   array : performance map, contained in the structure Type3254DataStruct.
    !   idx (integer) : one-dimensional index of the element to set, in column-major order.
    !   value (real(wp)) : the value to set at index idx.
    !   PMlength (integer) : length of the performance map.
        integer, intent(in) :: idx, PMlength
        real(wp) :: array(PMlength)
        real(wp), intent(in) :: value
        array(idx) = value
    end subroutine SetPMvalue


    function Findlb(array, value, extent) result(lowerBoundIndex)
    ! findlb finds the index of the element of an ordered array which is smaller
    ! and closest to a given value.
    ! Special cases: If the value is exactly equal to an element of the array
    ! with index i, the function returns i-1. If the value is smaller than the first
    ! (and lowest) value of the array, the function returns 1. If the value is higher
    ! than the last (and largest) value of the array, the function returns size(array) - 1.
    !
    ! Inputs
    !   array (1-D real(wp) array) : array to be searched.
    !   value (real(wp)) : value to search.
    !   extent (integer) : size of the array
    !
    ! Output
    !   lowerBoundIndex (integer) : index of the lower bound of the value in the array.
        real(wp), intent(in) :: array(:)
        real(wp), intent(in) :: value
        integer, intent(in) :: extent
        integer :: L, R, mid, lowerBoundIndex
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
        lowerBoundIndex = L - 1
        if (lowerBoundIndex == 0) lowerBoundIndex = 1
    end function Findlb


    function FullAdder(a, b, carry_in) result(results)
    ! FullAdder implements a full adder, to perform addition on binary numbers.
    !
    ! Inputs
    !   a, b (logicals) : operands of the addition
    !   carry_in (logical) : carry from the previous stage.
    !
    ! Outputs (given as a logical array of size 2)
    !   sum (logical) : a + b + carry_in (+ being exclusive or).
    !   carry_out : carry for the next stage.
        implicit none
        logical, intent(in) :: a, b, carry_in
        logical :: sum, carry_out, results(2)
        sum = a .neqv. b .neqv. carry_in
        carry_out = a .and. b .or. carry_in .and. (a .neqv. b)
        results = (/sum, carry_out /)
    end function FullAdder


    subroutine Increment(C)
    ! Increment increments a binary counter.
    !
    ! Inputs
    !   C (logical array) : the counter to increment,
    !                       with the rightmost bit being the lower level.
        implicit none
        logical, intent(inout) :: C(:)
        logical :: sumcarry(2), hasFalse
        integer :: N, k, i
        N = size(C)
        hasFalse = .false.
        do i = 1, N
            if (.not. C(i)) hasFalse = .true.
        end do
        if (.not. hasFalse) then  ! counter has reached max capacity
            C = .false.  ! reset counter
        else  ! increment counter using binary addition
            sumcarry = FullAdder(C(N), .true., .false.)
            C(N) = sumcarry(1)
            k = N - 1
            do while (sumcarry(2))
                sumcarry = FullAdder(C(k), .false., sumcarry(2))
                C(k) = sumcarry(1)
                k = k-1
            end do
        end if
    end subroutine Increment


    function Correction(defrost_mode, recov_penalty)
    ! Correction provides the correction factor to apply to the input power
    ! and heating capacity depending on the defrost mode.
    !
    ! Inputs
    !   defrost_mode (integer) : the defrost mode: defrost(0), recovery (1)
    !                            or steady-state (2)
    ! Output
    !   correction (real(wp) array) : array containing the power correction
    !                                 factor at (1) and the capacity correction
    !                                 factor at (2).
        integer, intent(in) :: defrost_mode
        real(wp), intent(in) :: recov_penalty
        real(wp) :: correction(2)
        if (defrost_mode == 0) then
            correction = (/0.6_wp, 0.0_wp/)
        else if (defrost_mode == 1) then
            correction = (/1.0_wp, recov_penalty/)
        else if (defrost_mode == 2) then
            correction = (/1.0_wp, 1.0_wp/)
        else
            ! add warning
        end if
    end function Correction

        
    subroutine ExecuteFirstCallOfSimulation
  	    call SetNumberofParameters(10)
  	    call SetNumberofInputs(16)
  	    call SetNumberofDerivatives(0)
  	    call SetNumberofOutputs(22)
  	    call SetIterationMode(1)
  	    call SetNumberStoredVariables(0, 0)
  	    call SetNumberofDiscreteControls(0)

        ! Allocate stored data structure
        if (.not. allocated(s)) then
            allocate(s(Ninstances))
        endif

        call ReadParameters()  ! required to get LUcool and LUheat
        call ReadPermap(LUcool, LUheat)

    end subroutine ExecuteFirstCallOfSimulation


    subroutine ExecuteStartTime
        call ReadParameters()
        call GetInputValues()
        call SetOutputValue(1, 0.0_wp)  ! Outlet air temperature
        call SetOutputValue(2, 0.0_wp)  ! Outlet air humidity ratio
        call SetOutputValue(3, 0.0_wp)  ! Outlet air % RH
        call SetOutputValue(4, 0.0_wp)  ! Outlet air pressure
        call SetOutputValue(5, 0.0_wp)  ! Outlet air flow rate
        call SetOutputValue(6, 0.0_wp)  ! Total cooling rate
        call SetOutputValue(7, 0.0_wp)  ! Sensible cooling rate
        call SetOutputValue(8, 0.0_wp)  ! Latent cooling rate
        call SetOutputValue(9, 0.0_wp)  ! Heat rejection rate
        call SetOutputValue(10, 0.0_wp)  ! Total heating rate
        call SetOutputValue(11, 0.0_wp)  ! Heat absorption rate
        call SetOutputValue(12, 0.0_wp)  ! Total power consumption
        call SetOutputValue(13, 0.0_wp)  ! COP
        call SetOutputValue(14, 0.0_wp)  ! EER
        call SetOutputValue(15, 0.0_wp)  ! Indoor fan power
        call SetOutputValue(16, 0.0_wp)  ! Outdoor fan power
        call SetOutputValue(17, 0.0_wp)  ! Compressor power
        call SetOutputValue(18, 0.0_wp)  ! Compressor frequency
        call SetOutputValue(19, 0.0_wp)  ! Condensate temperature
        call SetOutputValue(20, 0.0_wp)  ! Condensate flow rate
        call SetOutputValue(21, 0.0_wp)  ! Defrost mode
        call SetOutputValue(22, 0.0_wp)  ! Force shutdown
    end subroutine ExecuteStartTime


    subroutine ExecuteEndOfTimestep
        continue
    end subroutine ExecuteEndOfTimestep

    subroutine ExecuteLastCallOfSimulation
        deallocate(s(Ni)%extents)
        deallocate(s(Ni)%entries)
        deallocate(s(Ni)%lbounds)
        deallocate(s(Ni)%ubounds)
        deallocate(s(Ni)%PelcMap)
        deallocate(s(Ni)%PelhMap)
        deallocate(s(Ni)%QcsMap)
        deallocate(s(Ni)%QclMap)
        deallocate(s(Ni)%QhMap)
    end subroutine ExecuteLastCallOfSimulation


    subroutine ReadParameters
        psymode = GetParameterValue(1)
        PelcRated = GetParameterValue(2)
        QcRated = GetParameterValue(3)
        PelhRated = GetParameterValue(4)
        QhRated = GetParameterValue(5)
        AFRrated = GetParameterValue(6)
        freqRatedc = GetParameterValue(7)
        freqRatedh = GetParameterValue(8)
        LUcool = GetParameterValue(9)
        LUheat = GetParameterValue(10)
    end subroutine ReadParameters


    subroutine GetInputValues
        Tr = GetInputValue(1)
        wr = GetInputValue(2)
        RHr = GetInputValue(3)
        pr = GetInputValue(4)
        mDot = GetInputValue(5)
        Toa = GetInputValue(6)
        woa = GetInputValue(7)
        RHoa = GetInputValue(8)
        poa = GetInputValue(9)
        freq = GetInputValue(10)
        AFR = GetInputValue(11)
        mode = GetInputValue(12)
        defrost_mode = GetInputValue(13)
        recov_penalty = GetInputValue(14)
        PfanI = GetInputValue(15)
        PfanO = GetInputValue(16)
    end subroutine GetInputValues


    subroutine SetOutputValues
        call SetOutputValue(1, Ts)  ! Outlet air temperature
        call SetOutputValue(2, ws)  ! Outlet air humidity ratio
        call SetOutputValue(3, RHs*100.0_wp)  ! Outlet air % RH
        call SetOutputValue(4, ps)  ! Outlet air pressure
        call SetOutputValue(5, mDot)  ! Outlet air flow rate
        call SetOutputValue(6, Qc)  ! Total cooling rate
        call SetOutputValue(7, Qcs)  ! Sensible cooling rate
        call SetOutputValue(8, Qcl)  ! Latent cooling rate
        call SetOutputValue(9, Qrej)  ! Heat rejection rate
        call SetOutputValue(10, Qh)  ! Total heating rate
        call SetOutputValue(11, Qabs)  ! Heat absorption rate
        call SetOutputValue(12, Pel)  ! Total power consumption
        call SetOutputValue(13, COP)  ! COP
        call SetOutputValue(14, EER)  ! EER
        call SetOutputValue(15, PfanI)  ! Indoor fan power
        call SetOutputValue(16, PfanO)  ! Outdoor fan power
        call SetOutputValue(17, Pcomp)  ! Compressor power
        call SetOutputValue(18, fComp)  ! Compressor frequency
        call SetOutputValue(19, Tc)  ! Condensate temperature
        call SetOutputValue(20, cmfr)  ! Condensate flow rate
        call SetOutputValue(21, real(defrost_mode, wp))  ! Defrost mode
        call SetOutputValue(22, merge(1.0_wp, 0.0_wp, shutdown))  ! Force shutdown
    end subroutine SetOutputValues

end subroutine Type3254
