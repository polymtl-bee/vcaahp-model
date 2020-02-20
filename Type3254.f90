! +-----------------------------------------------------------------------------+
! | TRNSYS Type3254: Variable capacity air-air heat pump with performance files |
! +-----------------------------------------------------------------------------+
    
! This routine implements an air-air heat pump with variable speed compressor.

    
! Inputs
! ----------------------------------------------------------------------------------------------------
! Nb | Variable     | Description                                   | Input  Units    | Internal Units
! ----------------------------------------------------------------------------------------------------
!  1 | Tr           | Inlet (return) air temperature                | °C              | °C
!  2 | wr           | Inlet (return) air humidity ratio             | -               | -
!  3 | RHr          | Inlet (return) air relative humidity          | % (base 100)    | -
!  4 | amfr         | Inlet (return) air mass flow rate             | kg/h            | kg/h
!  5 | pr           | Inlet (return) air pressure                   | atm             | atm
!  6 | mode         | 0 = cooling mode                              | -               | -
!                   | 1 = heating mode                              |                 |
!  7 | Toa          | Outdoor air dry bulb temperature              | °C              | °C
!  8 | RHoa         | Outdoor air realtive humidity                 | % (base 100)    | -
!  9 | freq         | Compressor frequency                          | 1/s             | 1/s
! 10 | PfanI        | Indoor fan power                              | kJ/h            | kJ/h
! 11 | PfanO        | Outdoor fan power                             | kJ/h            | kJ/h
! ----------------------------------------------------------------------------------------------------

! Parameters
! ----------------------------------------------------------------------------------------------------
! Nb | Variable     | Description                                   | Param. Units    | Internal Units
! ----------------------------------------------------------------------------------------------------
!  1 | yHum         | 1 = Humidity ratio as humidity input          | -               | -
!                   | 2 = Relative humidity as humidity input       |                 |
!  2 | PelcRated    | Rated total cooling power                     | kJ/h            | kJ/h
!  3 | QcsRated     | Rated sensible cooling capacity               | kJ/h            | kJ/h
!  4 | QclRated     | Rated latent cooling capacity                 | kJ/h            | kJ/h
!  5 | PelhRated    | Rated total heating power                     | kJ/h            | kJ/h
!  6 | QhRated      | Rated heating capacity                        | kJ/h            | kJ/h
!  7 | amfrRated    | Rated inlet air mass flow rate                | kg/h            | kg/h
!  8 | freqRated    | Rated frequency                               | 1/s             | 1/s
!  9 | t_defrost    | Defrost duration during a cycle               | hr              | hr
! 10 | Tcutoff      | Defrost cutoff temperature                    | °C              | °C
! 11 | LUcool       | Logical Unit - cooling mode                   | -               | -
! 12 | LUheat       | Logical Unit - heating mode                   | -               | -
! ----------------------------------------------------------------------------------------------------

! Outputs
! ----------------------------------------------------------------------------------------------------
! Nb | Variable     | Description                                   | Output  Units   | Internal Units
! ----------------------------------------------------------------------------------------------------
!  1 | Ts           | Outlet (supply) air temperature               | °C              | °C
!  2 | ws           | Outlet (supply) air humidity ratio            | -               | -
!  3 | RHs          | Outlet (supply) air % RH                      | % (base 100)    | % (base 100)
!  4 | amfr         | Outlet (supply) air mass flow rate            | kg/h            | kg/h
!  5 | ps           | Outlet (supply) air pressure                  | atm             | atm
!  6 | Qc           | Total cooling rate                            | kJ/h            | kJ/h
!  7 | Qcs          | Sensible cooling rate                         | kJ/h            | kJ/h
!  8 | Qcl          | Latent cooling rate                           | kJ/h            | kJ/h
!  9 | Qrej         | Heat rejection rate                           | kJ/h            | kJ/h
! 10 | Qh           | Total heating rate                            | kJ/h            | kJ/h
! 11 | Qabs         | Heat absorption rate                          | kJ/h            | kJ/h
! 12 | Pel          | Total power consumption                       | kJ/h            | kJ/h
! 13 | COP          | Coefficient of performance                    | -               | -
! 14 | EER          | Energy efficiency rating                      | -               | -
! 15 | PfanI        | Indoor fan power                              | kJ/h            | kJ/h
! 16 | PfanO        | Outdoor fan power                             | kJ/h            | kJ/h
! 17 | Pcomp        | Compressor power                              | kJ/h            | kJ/h
! 18 | Tc           | Condensate temperature                        | °C              | °C
! 19 | cmfr         | Condensate mass flow rate                     | kg/h            | kg/h
! 20 | defrost_mode | 0 = defrost (off) mode                        | -               | -
!                   | 1 = Recovery mode (transient)                 |                 |
!                   | 2 = Steady-state mode                         |                 |
! ----------------------------------------------------------------------------------------------------

module Type3254Data

use, intrinsic :: iso_fortran_env, only : wp=>real64    ! Defines a constant "wp" (working precision) that can be used in real numbers, e.g. 1.0_wp, and sets it to real64 (double precision)
use TrnsysConstants
use TrnsysFunctions
implicit none

type Type3254DataStruct
    
    ! Parameters
    real(wp), allocatable :: entries(:, :, :)
    integer, allocatable :: extents(:, :)
    integer :: PMClength, PMHlength  ! Length of the flattened performance map

    ! Performance matrices
    real(wp), allocatable :: PelcMap(:, :, :, :, :)
    real(wp), allocatable :: QcsMap(:, :, :, :, :)
    real(wp), allocatable :: QclMap(:, :, :, :, :)
    real(wp), allocatable :: PelhMap(:, :, :, :, :)
    real(wp), allocatable :: QhMap(:, :, :, :, :)

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
real(wp) :: Tr, wr, RHr, amfr, pr, Toa, RHoa, freq, PfanI, PfanO  ! Inputs
integer :: mode  ! operating mode
integer :: yHum, LUcool, LUheat  ! Parameters
real(wp) :: PelcRated, QcsRated, QclRated, PelhRated, QhRated, amfrRated, freqRated  ! Parameters (rated values)
real(wp) :: Tcutoff, t_off  ! defrost parameters
real(wp) :: Ts, ws, RHs, ps  ! Outputs (supply conditions)
real(wp) :: Pel, Qc, Qcs, Qcl, Qrej, Qh, Qabs, Pcomp  ! Outputs (heat and power)
real(wp) :: COP, EER, Tc, cmfr  ! Outputs (misc)


! Local variables
real(wp) :: psydat(9), Twbr, Twboa, hr, hx, hs
integer :: psymode, status
integer, parameter :: Ninstances = 1  ! Number of units
integer :: Ni = 1  ! temporary, should use a kernel function to get the actual instance number.

! Defrost variables
integer :: defrost_mode = 2
real(wp) :: t_uc = 0, t_oc = 0  ! time under/over cutoff temperature
real(wp) :: t_cy  ! duration of a full defrost cycle
real(wp) :: t_rec  ! time between defrost and steady-state (recovery mode)
real(wp) :: t_ss  ! steady-state time in a cycle
real(wp) :: tau  ! time constant for the recovery penalty
real(wp) :: t_ld  ! time since last defrost
real(wp) :: defrost_corr(2) = 0.0_wp  ! correction factors

! Performance map reading variables
integer, parameter :: Nc = 5, Nh = 5  ! Number of interpolation variables
integer, parameter :: Nmax = max(Nc, Nh)
integer :: i, j, N, Noutc = 3, Nouth = 2, Nout


! Interpolation variables
real(wp) :: point(Nmax, 0:1)
real(wp), allocatable :: interpolationResults(:)

! Set the version number for this Type
if (GetIsVersionSigningTime()) then
    call SetTypeVersion(18)
    return
endif

call GetTRNSYSvariables()
call ExecuteSpecialCases()
call GetInputValues()
if (ErrorFound()) return

! Ni = GetCurrentUnit()
N = (1-mode) * Nc + mode * Nh
Nout = Noutc*(1-mode) + Nouth*mode

! Return air state
psydat(1) = pr
psydat(2) = Tr
psydat(4) = RHr/100.0_wp
psydat(6) = wr
if (yHum == 1) then
    psymode = 4
else
    psymode = 2
endif
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

if (mode == 1) then ! compute outdoor air wet bulb
    psydat(1) = pr
    psydat(2) = Toa
    psydat(4) = RHoa/100.0_wp
    call MoistAirProperties(thisUnit, thisType, 1, 2, 1, psydat, 1, status)
    Twboa = psydat(3)
    call SetDefrostMode()
    defrost_corr = Correction(defrost_mode, t_ld, tau)
end if

! Interpolate using wet bulb
point(1, mode) = Tr
if (mode == 0) then
    point(2, mode) = Twbr
    point(3, mode) = Toa
else
    point(2, mode) = Toa
    point(3, mode) = Twboa
end if
point(4, mode) = amfr / amfrRated
point(5, mode) = freq / freqRated
allocate(interpolationResults(Nout))
interpolationResults = interpolate(point(:, mode), mode, Nout)
if (mode == 0) Pel = interpolationResults(1) * PelcRated
if (mode == 1) Pel = interpolationResults(1) * PelhRated * defrost_corr(1)
Qcs = interpolationResults(2) * QcsRated
Qcl = interpolationResults(3) * QclRated
Qh = interpolationResults(4) * QhRated * defrost_corr(2)
Qc = Qcs + Qcl
deallocate(interpolationResults)

! Supply air state
ps = pr  ! Fan pressure drop neglected

! Moist air state
if (Qc < Qcs) then
    Qc = Qcs
    ! Add warning
endif

if (amfr /= 0.0_wp) then
    ws = wr
    if (mode == 0) then
        hs = hr - Qc/amfr
        hx = hr  ! useful when the following if clause is not true
        if (Qcl > 0.0_wp) then  ! compute humidity after condensation
            psydat(1) = pr
            psydat(2) = Tr
            hx = hr - Qcl/amfr
            psydat(7) = hx  ! enthalpy of the state (Tr, ws)
            call MoistAirProperties(thisUnit, thisType, 1, 5, 0, psydat, 1, status)  ! dry-bulb and enthalpy as inputs
            if (ErrorFound()) return
            ws = psydat(6)
        endif
    else
        hs = hr + Qh/amfr
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

if (mode == 0) then
    ! Re-calculate heat transfer whose value is modified if saturation occurs
    Qcs = amfr * (hx - hs)  ! Sensible cooling rate
    Qcl = amfr * (hr - hx)  ! Latent cooling rate
    Qc = Qcs + Qcl  ! Total cooling rate
    Qrej = Qc + Pel  ! Heat rejection
    Qabs = 0.0_wp
else
    Qrej = 0.0_wp
    Qabs = Qh - Pel
end if
Pcomp = Pel - PfanI - PfanO  ! Compressor power
if (Pel /= 0.0_wp) then
    COP = (Qc + Qh) / Pel
else
    COP = 0.0_wp
endif
EER = 3.413_wp * COP
Tc = Ts
cmfr = amfr * (wr - ws)  ! Condensate flow rate - water balance

call SetOutputValues()

return

    contains
    
    subroutine ReadPermap
        character (len=maxPathLength) :: permapCoolPath
        character (len=maxPathLength) :: permapHeatPath
        integer :: nTr, nTwbr, nToa, nTwboa, namfr, nfreq  ! number of entries for each variable
        real(wp) :: filler(Nmax)
    
        ! Ni = GetCurrentUnit()
    
        permapCoolPath = GetLUfileName(LUcool)
        permapHeatPath = GetLUfileName(LUheat)
        call CheckPMfile(permapCoolPath)
        call CheckPMfile(permapHeatPath)
    
        open(LUcool, file=permapCoolPath, status='old')
        open(LUHeat, file=permapHeatPath, status='old')
    
            do i = 1, 6  ! Skip 6 first lines
                read(LUcool, *)
                read(LUHeat, *)
            enddo
        allocate(s(Ni)%extents(Nmax, 0:1))
        do i = 1, Nc
            read(LUcool, *)  ! Skip a line
            read(LUcool, *) s(Ni)%extents(i, 0)
        end do
        do i = 1, Nh
            read(LUheat, *)  ! Skip a line
            read(LUheat, *) s(Ni)%extents(i, 1)
        end do
        s(Ni)%PMClength = product(s(Ni)%extents(:, 0))
        s(Ni)%PMHlength = product(s(Ni)%extents(:, 1))
        allocate(s(Ni)%entries(maxval(s(Ni)%extents), Nmax, 0:1))
        do i = 1, Nc
            read(LUcool, *)  ! Skip a line
            read(LUcool, *) (s(Ni)%entries(j, i, 0), j = 1, s(Ni)%extents(i, 0))
        end do
        do i = 1, Nh
            read(LUheat, *)  ! Skip a line
            read(LUheat, *) (s(Ni)%entries(j, i, 1), j = 1, s(Ni)%extents(i, 1))
        end do
            do i = 1, 4  ! Skip 4 lines
                read(LUcool, *)
                read(LUheat, *)
            end do
            
        nTr = s(Ni)%extents(1, 0)
        nTwbr = s(Ni)%extents(2, 0)
        nToa = s(Ni)%extents(3, 0)
        namfr = s(Ni)%extents(4, 0)
        nfreq = s(Ni)%extents(5, 0)
        allocate(s(Ni)%PelcMap(nTr, nTwbr, nToa, namfr, nfreq))
        allocate(s(Ni)%QcsMap(nTr, nTwbr, nToa, namfr, nfreq))
        allocate(s(Ni)%QclMap(nTr, nTwbr, nToa, namfr, nfreq))
        do i = 1, s(Ni)%PMClength
            read(LUcool, *) (filler(j), j = 1, Nc), Pel, Qcs, Qcl
            call SetPMvalue(s(Ni)%PelcMap, i, Pel, s(Ni)%PMClength)
            call SetPMvalue(s(Ni)%QcsMap, i, Qcs, s(Ni)%PMClength)
            call SetPMvalue(s(Ni)%QclMap, i, Qcl, s(Ni)%PMClength)
        end do
    
        close(LUcool)
        
        nTr = s(Ni)%extents(1, 1)
        nToa = s(Ni)%extents(2, 1)
        nTwboa = s(Ni)%extents(3, 1)
        namfr = s(Ni)%extents(4, 1)
        nfreq = s(Ni)%extents(5, 1)
        allocate(s(Ni)%PelhMap(nTr, nToa, nTwboa, namfr, nfreq))
        allocate(s(Ni)%QhMap(nTr, nToa, nTwboa, namfr, nfreq))
        do i = 1, s(Ni)%PMHlength
            read(LUheat, *) (filler(j), j = 1, Nh), Pel, Qh
            call SetPMvalue(s(Ni)%PelhMap, i, Pel, s(Ni)%PMHlength)
            call SetPMvalue(s(Ni)%QhMap, i, Qh, s(Ni)%PMHlength)
        end do
        
        close(LUheat)
    
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
    
    
    function Interpolate(point, mode, Nout)
        real(wp), intent(in) :: point(N)
        real(wp) :: scaled_point(N), LBvalue, UBvalue, sp
        real(wp), allocatable :: hypercube(:, :)
        integer, intent(in) :: mode
        integer, dimension(N) :: idx, lb_idx, counter_int, ones, zeros
        integer :: i
        logical :: counter_bool(N)
        integer, intent(in) :: Nout
        real(wp) :: interpolate(Noutc + Nouth - 1)
        zeros = 0
        ones = 1
        do i = 1, N
            j = findlb(s(Ni)%entries(:, i, mode), point(i), s(Ni)%extents(i, mode))
            lb_idx(i) = j
            LBvalue = s(Ni)%entries(j, i, mode)
            UBvalue = s(Ni)%entries(j+1, i, mode)
            scaled_point(i) = (point(i) - LBvalue) / (UBvalue - LBvalue)
        end do
        allocate(hypercube(Nout, 2**N))
        counter_bool = .true.
        do i = 1, 2**N
            call increment(counter_bool)
            counter_int = merge(ones, zeros, counter_bool)
            idx = lb_idx + counter_int
            hypercube(:, i) = Vertex(i, idx)
        end do

        do i = 1, N
            sp = scaled_point(i)
            j = N - i
            hypercube(:, :2**j) = (1-sp) * hypercube(:, :2**j) &
                                    + sp * hypercube(:, 2**j+1:2**(j+1))
        end do
        if (mode == 0) then
            do i = 1, Noutc
                interpolate(i) = hypercube(i, 1)
            end do
            do i = Noutc+1, Noutc+Nouth-1
                interpolate(i) = 0.0_wp
            end do
        else
            interpolate(1) = hypercube(1, 1)
            do i = 2, Noutc
                interpolate(i) = 0.0_wp
            end do
            do i = Noutc+1, Noutc+Nouth-1
                interpolate(i) = hypercube(i-Noutc+1, 1)
            end do
        end if
        deallocate(hypercube)
    end function Interpolate
    
    
    function Vertex(i, idx)
        integer, intent(in) :: i, idx(:)
        real(wp) :: Pel, Qcs, Qcl, vertex(Nout)
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
    
    
    function GetPMvalue(mode, array, idx)
        integer, intent(in) :: idx(:)
        integer :: mode, i, array_idx
        real(wp) :: array((1-mode)*s(Ni)%PMClength + mode*s(Ni)%PMHlength)
        real(wp) :: GetPMvalue
        array_idx = idx(N)
        do i = N-1, 1, -1
            array_idx = array_idx + product(s(Ni)%extents(i+1:, mode)) * (idx(i) - 1)
        end do
        GetPMvalue = array(array_idx)
    end function GetPMvalue
    
    
    subroutine SetPMvalue(array, idx, value, PMlength)
        integer, intent(in) :: idx, PMlength
        real(wp) :: array(PMlength)
        real(wp), intent(in) :: value
        array(idx) = value
    end subroutine SetPMvalue
    
    
    function findlb(array, value, extent)
        real(wp), intent(in) :: array(:)
        real(wp), intent(in) :: value
        integer, intent(in) :: extent
        integer :: findlb
        integer :: L, R, mid
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
        findlb = L - 1
        if (findlb == 0) findlb = 1
    end function findlb
    
    
    function full_adder(a, b, carry_in)
        implicit none
        logical, intent(in) :: a, b, carry_in
        logical :: sum, carry_out, full_adder(2)
        sum = a .neqv. b .neqv. carry_in
        carry_out = a .and. b .or. carry_in .and. (a .neqv. b)
        full_adder = (/sum, carry_out /)
    end function full_adder

    
    subroutine increment(C)
        implicit none
        logical, intent(inout) :: C(:)
        logical :: sumcarry(2), hasFalse
        integer :: N, k, i
        N = size(C)
        hasFalse = .false.
        do i = 1, N
            if (C(i) .neqv. .true.) hasFalse = .true.
        end do
        if (hasFalse .neqv. .true.) then  ! reset if all true
            C = .false.
        else
            sumcarry = full_adder(C(N), .true., .false.)
            C(N) = sumcarry(1)
            k = N - 1
            do while (sumcarry(2))
                sumcarry = full_adder(C(k), .false., sumcarry(2))
                C(k) = sumcarry(1)
                k = k-1
            end do
        end if
    end subroutine increment
    
    
    subroutine SetDefrostMode
        t_cy = 156
        t_rec = 30
        t_ss = t_cy - t_off - t_rec
        tau = t_rec - t_off
        t_ld = GetDynamicArrayValueLastTimestep(1)
        t_uc = GetDynamicArrayValueLastTimestep(2)
        t_oc = GetDynamicArrayValueLastTimestep(3)
        
        if (Toa < Tcutoff) then
            t_uc = t_uc + dt
        else
            t_oc = t_oc + dt
        end if
    
        if (t_ld < t_rec) then
            defrost_mode = 1
            t_ld = t_ld + dt
        else if (t_ld < t_rec + t_ss) then
            defrost_mode = 2
            t_ld = t_ld + dt
        else if (t_uc < t_oc) then
            t_uc = 0.0_wp
            t_oc = 0.0_wp
            t_ld = t_rec
            defrost_mode = 2
        else if (t_ld < t_cy) then
            t_ld = t_ld + dt
            defrost_mode = 0
        else
            t_ld = 0.0_wp
            t_uc = 0.0_wp
            t_oc = 0.0_wp
            defrost_mode = 1
        end if
        
        call SetDynamicArrayValueThisIteration(1, t_ld)
        call SetDynamicArrayValueThisIteration(2, t_uc)
        call SetDynamicArrayValueThisIteration(3, t_oc)
    end subroutine SetDefrostMode
    
    
    function Correction(defrost_mode, t_ld, tau)
        integer, intent(in) :: defrost_mode
        real(wp), intent(in) :: t_ld, tau
        real(wp) :: correction(2), recov_penalty
        if (defrost_mode == 0) then
            correction = (/0.6_wp, 0.0_wp/)
        else if (defrost_mode == 1) then
            recov_penalty = 1.0_wp - exp(-1-t_ld/tau)
            correction = (/1.0_wp, recov_penalty/)
        else if (defrost_mode == 2) then
            correction = (/1.0_wp, 1.0_wp/)
        else
            ! add warning
        end if
    end function Correction
    
    
    subroutine ExecuteSpecialCases
    
    ! All the stuff that must be done once at the beginning
    if(getIsFirstCallofSimulation()) then

  	    ! Tell the TRNSYS engine how this Type works
  	    call SetNumberofParameters(12)
  	    call SetNumberofInputs(11)
  	    call SetNumberofDerivatives(0)
  	    call SetNumberofOutputs(20)
  	    call SetIterationMode(1)
  	    call SetNumberStoredVariables(0, 4)
  	    call SetNumberofDiscreteControls(0)
        
        ! Allocate stored data structure
        if (.not. allocated(s)) then
            allocate(s(Ninstances))
        endif
        
        call ReadParameters()
        call ReadPermap()

  	    return

    endif
    
    ! Start of the first timestep: no iterations, outputs initial conditions
    if (getIsStartTime()) then
        
        call ReadParameters()
        call GetInputValues()

        ! Set outputs to zeros at initial time
	    call SetOutputValue(1, 0.0_wp)  ! Supply air temperature
	    call SetOutputValue(2, 0.0_wp)  ! Supply air humidity ratio
	    call SetOutputValue(3, 0.0_wp)  ! Supply air % RH
	    call SetOutputValue(4, 0.0_wp)  ! Supply air flow rate
	    call SetOutputValue(5, 0.0_wp)  ! Supply air pressure
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
	    call SetOutputValue(18, 0.0_wp)  ! Condensate temperature
	    call SetOutputValue(19, 0.0_wp)  ! Condensate flow rate
        call SetOutputValue(20, 0.0_wp)  ! Defrost mode
        
        call SetDynamicArrayInitialValue(1, 0.0_wp)
        call SetDynamicArrayInitialValue(2, 0.0_wp)
        call SetDynamicArrayInitialValue(3, 0.0_wp)

        return

    endif
    
    ! Parameters must be re-read - indicates another unit of this Type
    if(getIsReReadParameters()) call ReadParameters()
    
    ! End of timestep call (after convergence or too many iterations)
    if (GetIsEndOfTimestep()) then
        return  ! We are done for this call
    endif
    
    if (GetIsLastCallofSimulation()) then
        return  ! We are done for this call
    endif
    
    end subroutine ExecuteSpecialCases
    
    
    subroutine ReadParameters
        yHum = GetParameterValue(1)
        PelcRated = GetParameterValue(2)
        QcsRated = GetParameterValue(3)
        QclRated = GetParameterValue(4)
        PelhRated = GetParameterValue(5)
        QhRated = GetParameterValue(6)
        amfrRated = GetParameterValue(7)
        freqRated = GetParameterValue(8)
        t_off = GetParameterValue(9)
        Tcutoff = GetParameterValue(10)
        LUcool = GetParameterValue(11)
        LUheat = GetParameterValue(12)
    end subroutine ReadParameters
    
    
    subroutine GetInputValues
        Tr = GetInputValue(1)
        wr = GetInputValue(2)
        RHr = GetInputValue(3)
        amfr = GetInputValue(4)
        pr = GetInputValue(5)
        mode = GetInputValue(6)
        Toa = GetInputValue(7)
        RHoa = GetInputValue(8)
        freq = GetInputValue(9)
        PfanI = GetInputValue(10)
        PfanO = GetInputValue(11)
    end subroutine GetInputValues
    
    
    subroutine SetOutputValues
        call SetOutputValue(1, Ts)  ! Outlet air temperature
        call SetOutputValue(2, ws)  ! Outlet air humidity ratio
        call SetOutputValue(3, RHs*100.0_wp)  ! Outlet air % RH
        call SetOutputValue(4, amfr)  ! Outlet air flow rate
        call SetOutputValue(5, ps)  ! Outlet air pressure
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
        call SetOutputValue(18, Tc)  ! Condensate temperature
        call SetOutputValue(19, cmfr)  ! Condensate flow rate
        call SetOutputValue(20, real(defrost_mode, wp))  ! Defrost mode
    end subroutine SetOutputValues
    
    
    subroutine GetTRNSYSvariables
        time = getSimulationTime()
        dt = getSimulationTimeStep()
        thisUnit = getCurrentUnit()
        thisType = getCurrentType()
    end subroutine GetTRNSYSvariables

end subroutine Type3254
