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
!  5 | pIn          | Inlet (return) air pressure                   | atm             | atm
!  6 | Toa          | Outdoor air temperature                       | °C              | °C
!  7 | freq         | Compressor frequency                          | 1/s             | 1/s
!  8 | PfanI        | Indoor fan power                              | kJ/h            | kJ/h
!  9 | PfanO        | Outdoor fan power                             | kJ/h            | kJ/h
! ----------------------------------------------------------------------------------------------------

! Parameters
! ----------------------------------------------------------------------------------------------------
! Nb | Variable     | Description                                   | Param. Units    | Internal Units
! ----------------------------------------------------------------------------------------------------
!  1 | yHum         | 1 = Humidity ratio as humidity input          | -               | -
!    |              | 2 = Relative humidity as humidity input       | -               | -
!  2 | PelRated     | Rated total cooling power                     | kJ/h            | kJ/h
!  3 | QcsRated     | Rated sensible cooling capacity               | kJ/h            | kJ/h
!  4 | QclRated     | Rated latent cooling capacity                 | kJ/h            | kJ/h
!  5 | amfrRated    | Rated inlet air mass flow rate                | kg/h            | kg/h
!  6 | freqRated    | Rated frequency                               | 1/s             | 1/s
!  7 | LUcool       | Logical Unit - cooling mode                   | -               | -
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
! 10 | Pel          | Total power consumption                       | kJ/h            | kJ/h
! 11 | COP          | Coefficient of performance                    | -               | -
! 12 | EER          | Energy efficiency rating                      | -               | -
! 13 | PfanI        | Indoor fan power                              | kJ/h            | kJ/h
! 14 | PfanO        | Outdoor fan power                             | kJ/h            | kJ/h
! 15 | Pcomp        | Compressor power                              | kJ/h            | kJ/h
! 16 | Tc           | Condensate temperature                        | °C              | °C
! 17 | cmfr         | Condensate mass flow rate                     | kg/h            | kg/h
! ----------------------------------------------------------------------------------------------------

module Type3254Data

use, intrinsic :: iso_fortran_env, only : wp=>real64    ! Defines a constant "wp" (working precision) that can be used in real numbers, e.g. 1.0_wp, and sets it to real64 (double precision)
use TrnsysConstants
use TrnsysFunctions
implicit none

type Type3254DataStruct
    
    ! Parameters
    real(wp) :: QcRated, QcsRated, PtotRated, mDotInRated, fRated  ! rated values
    real(wp), allocatable :: entries(:, :)
    integer, allocatable :: extents(:)

    ! Performance matrices
    real(wp), allocatable :: PelMap(:, :, :, :, :)
    real(wp), allocatable :: QcsMap(:, :, :, :, :)
    real(wp), allocatable :: QclMap(:, :, :, :, :)

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
real(wp) :: time, timestep  ! TRNSYS time and timestep

! Proforma variables
real(wp) :: Tr, wr, RHr, amfr, pr, Toa, freq, PfanI, PfanO  ! Inputs
integer :: yHum, LUcool  ! Parameters
real(wp) :: PelRated, QcsRated, QclRated, amfrRated, freqRated  ! Parameters (rated values)
real(wp) :: Ts, ws, RHs, ps  ! Outputs (supply conditions)
real(wp) :: Pel, Qc, Qcs, Qcl, Qrej, Pcomp  ! Outputs (heat and power)
real(wp) :: COP, EER, Tc, cmfr  ! Outputs (misc)


! Local variables
real(wp) :: psydat(9), Twbr, hr, hx, hs
integer :: psymode, status
integer, parameter :: Ninstances = 1  ! Number of units
integer :: Ni = 1  ! temporary, should use a kernel function to get the actual instance number.

! Performance map reading variables
integer, parameter :: N = 5  ! Number of interpolation variables
character (len=maxPathLength) :: permapPath
character (len=maxMessageLength) :: msg
logical :: permapFileFound = .false.
integer :: nTr, nRHr, nTo, namfr, nfreq  ! number of entries for each variable
integer :: i, j, line_count = 1, prev_line
real(wp), allocatable, dimension(:) :: TrValues,  RHrValues, ToValues, amfrValues, freqValues
integer :: idx(N)
real(wp) :: filler(N)

! Interpolation variables
real(wp) :: hypercube(2**N, 3), lbvalue, ubvalue
integer :: lb_idx(N)
!real(wp) :: point(N) = (/22.5, 0.4, 24.0, 0.7, 0.2/), scaled_point(N), sp
real(wp) :: point(N), scaled_point(N), sp
integer :: counter_int(N), ones(N) = 1, zeros(N) = 0
logical :: counter_bool(N) = .false.
counter_int = merge(ones, zeros, counter_bool)

! Set the version number for this Type
if (GetIsVersionSigningTime()) then
    call SetTypeVersion(18)
    return
endif

call GetTRNSYSvariables()
call ExecuteSpecialCases()
call GetInputValues()
if (ErrorFound()) return

point(1) = Tr
point(2) = RHr / 100_wp
point(3) = Toa
point(4) = amfr / amfrRated
point(5) = freq / freqRated
do i = 1, N
    j = findlb(s(Ni)%entries(i, :), point(i), s(Ni)%extents(i))
    lb_idx(i) = j
    lbvalue = s(Ni)%entries(i, j)
    ubvalue = s(Ni)%entries(i, j+1)
    scaled_point(i) = (point(i) - lbvalue) / (ubvalue - lbvalue)
end do

Pel = GetPMvalue(s(Ni)%PelMap, lb_idx)
Qcs = GetPMvalue(s(Ni)%QcsMap, lb_idx)
Qcl = GetPMvalue(s(Ni)%QclMap, lb_idx)
hypercube(1, :) = (/Pel, Qcs, Qcl/)
counter_bool = .false.
do i = 2, 2**N
    call increment(counter_bool)
    counter_int = merge(ones, zeros, counter_bool)
    idx = lb_idx + counter_int
    Pel = GetPMvalue(s(Ni)%PelMap, idx)
    Qcs = GetPMvalue(s(Ni)%QcsMap, idx)
    Qcl = GetPMvalue(s(Ni)%QclMap, idx)
    hypercube(i, :) = (/Pel, Qcs, Qcl/)
end do

do i = 1, N
    sp = scaled_point(i)
    j = N - i
    hypercube(:2**j, :) = (1-sp) * hypercube(:2**j, :) &
                            + sp * hypercube(2**j+1:2**(j+1), :)
end do

Pel = hypercube(1, 1) * PelRated
Qcs = hypercube(1, 2) * QcsRated
Qcl = hypercube(1, 3) * QclRated
Qc = Qcs + Qcl

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
call MoistAirProperties(thisUnit, thisType, 1, psymode, 0, psydat, 1, status)
! (unit, type, si units used, psych inputs, Twb not computed, inputs, warning mgmt, warning occurences)
pr = psydat(1)
Tr = psydat(2)
Twbr = psydat(3)
RHr = psydat(4)  ! RHr between 0 and 1 (not 0 and 100)
wr = psydat(6)
hr = psydat(7)

! Supply air state
ps = pr  ! Fan pressure drop neglected

! Moist air state
if (Qc < Qcs) then
    Qc = Qcs
    ! Add warning
endif

if (amfr /= 0.0_wp) then
    hs = hr - Qc/amfr
    ws = wr  ! useful when the following if clause is not true
    hx = hr  ! same
    if (Qcl > 0.0_wp) then
        psydat(1) = pr
        psydat(2) = Tr
        hx = hr - Qcl/amfr
        psydat(7) = hx  ! enthalpy of the state (Tr, ws)
        call MoistAirProperties(thisUnit, thisType, 1, 5, 0, psydat, 1, status)  ! dry-bulb and enthalpy as inputs
        if (ErrorFound()) return
        ws = psydat(6)
    endif
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
Qcs = amfr * (hx - hs)  ! Sensible cooling rate
Qcl = amfr * (hr - hx)  ! Latent cooling rate
Qc = Qcs + Qcl  ! Total cooling rate
Qrej = Qc + Pel  ! Heat rejection
Pcomp = Pel - PfanI - PfanO  ! Compressor power
if (Pel /= 0.) then
    COP = Qc / Pel
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
    
    idx = 1
    !Ni = GetCurrentUnit()
    
    permapPath = GetLUfileName(LUcool)
    inquire(file=trim(permapPath), exist=permapFileFound)
    if ( .not. permapFileFound ) then
        write(msg,'("""",a,"""")') trim(permapPath)
        msg = "Could not find the specified performance map file. Searched for: " // trim(msg)
        call Messages(-1, msg, 'fatal', thisUnit, thisType)
        return
    end if
    
    open(LUcool, file=permapPath, status='old')
    
        do i = 1, 7  ! Skip 7 lines
            read(LUcool, *)
        enddo
    read(LUcool, *) nTr  ! Read number of Tr values
        read(LUcool, *)  ! Skip line
    allocate(TrValues(nTr))  ! Read Tr values
    read(LUcool, *) (TrValues(i), i = 1, nTr)
        read(LUcool, *)  ! Skip line
    read(LUcool, *) nRHr
        read(LUcool, *)
    allocate(RHrValues(nRHr))
    read(LUcool, *) (RHrValues(i), i = 1, nRHr)
        read(LUcool, *)
    read(LUcool, *) nTo
        read(LUcool, *)
    allocate(ToValues(nTo))
    read(LUcool, *) (ToValues(i), i = 1, nTo)
        read(LUcool, *)
    read(LUcool, *) namfr
        read(LUcool, *)
    allocate(amfrValues(namfr))
    read(LUcool, *) (amfrValues(i), i = 1, namfr)
        read(LUcool, *)
    read(LUcool, *) nfreq
        read(LUcool, *)
    allocate(freqValues(nfreq))
    read(LUcool, *) (freqValues(i), i = 1, nfreq)
        do i = 1, 4  ! Skip 4 lines
            read(LUcool, *)
        end do
    allocate(s(Ni)%PelMap(nTr, nRHr, nTo, namfr, nfreq))
    allocate(s(Ni)%QcsMap(nTr, nRHr, nTo, namfr, nfreq))
    allocate(s(Ni)%QclMap(nTr, nRHr, nTo, namfr, nfreq))
    allocate(s(Ni)%extents(N))
    s(Ni)%extents = shape(s(Ni)%PelMap)
    read(LUcool, *) (filler(i), i = 1, N), Pel, Qcs, Qcl
    call SetPMvalue(s(Ni)%PelMap, idx, Pel)
    call SetPMvalue(s(Ni)%QcsMap, idx, Qcs)
    call SetPMvalue(s(Ni)%QclMap, idx, Qcl)
    do while (sum(idx) < sum(s(Ni)%extents))
        prev_line = line_count
        i = N
        do while (prev_line == line_count)
            if (idx(i) == s(Ni)%extents(i)) then
                i = i - 1
                if (i < N) idx(i+1) = 1
            else
                idx(i) = idx(i) + 1
                line_count = line_count + 1
                read(LUcool, *) (filler(j), j = 1, N), Pel, Qcs, Qcl
                call SetPMvalue(s(Ni)%PelMap, idx, Pel)
                call SetPMvalue(s(Ni)%QcsMap, idx, Qcs)
                call SetPMvalue(s(Ni)%QclMap, idx, Qcl)
            end if
        end do
    end do
    
    close(LUcool)
    ! make a check ?
    
    allocate(s(Ni)%entries(N, maxval(s(Ni)%extents)))
    s(Ni)%entries(1, 1:nTr) = TrValues
    s(Ni)%entries(2, 1:nRHr) = RHrValues
    s(Ni)%entries(3, 1:nTo) = ToValues
    s(Ni)%entries(4, 1:namfr) = amfrValues
    s(Ni)%entries(5, 1:nfreq) = freqValues
    
    end subroutine ReadPermap
    
    
    function GetPMvalue(array, idx)
        integer, intent(in) :: idx(:)
        integer :: i, array_idx
        real(wp) :: array(product(s(Ni)%extents))
        real(wp) :: GetPMvalue
        array_idx = idx(1)
        do i = 2, size(idx)
            array_idx = array_idx + product(s(Ni)%extents(:i-1)) * (idx(i) - 1)
        end do
        GetPMvalue = array(array_idx)
    end function GetPMvalue
    
    
    subroutine SetPMvalue(array, idx, value)
        real(wp) :: array(product(s(Ni)%extents))
        real(wp), intent(in) :: value
        integer, intent(in) :: idx(:)
        integer :: i, array_idx
        array_idx = idx(1)
        do i = 2, size(idx)
            array_idx = array_idx + product(s(Ni)%extents(:i-1)) * (idx(i) - 1)
        end do
        array(array_idx) = value
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
        logical :: sumcarry(2)
        integer :: N, k
        N = size(C)
        sumcarry = full_adder(C(N), .true., .false.)
        C(N) = sumcarry(1)
        k = N - 1
        do while (sumcarry(2))
            sumcarry = full_adder(C(k), .false., sumcarry(2))
            C(k) = sumcarry(1)
            k = k-1
        end do
    end subroutine increment
    
    
    subroutine ExecuteSpecialCases
    
    ! All the stuff that must be done once at the beginning
    if(getIsFirstCallofSimulation()) then

  	    ! Tell the TRNSYS engine how this Type works
  	    call SetNumberofParameters(7)
  	    call SetNumberofInputs(9)
  	    call SetNumberofDerivatives(0)
  	    call SetNumberofOutputs(17)
  	    call SetIterationMode(1)
  	    call SetNumberStoredVariables(0,0)
  	    call SetNumberofDiscreteControls(0)
        
        ! Allocate stored data structure
        if (.not. allocated(s)) then
            allocate(s(Ninstances))
        endif
        
        call ReadParameters()
        call ReadPermap()

        ! Set units (optional)
        call SetInputUnits(1, 'TE1')  ! °C
        call SetInputUnits(2, 'DM1')  ! -
        call SetInputUnits(3, 'PC1')  ! %
        call SetInputUnits(4, 'MF1')  ! kg/h
        call SetInputUnits(5, 'PR4')  ! atm
        call SetInputUnits(6, 'TE1')  ! °C
        ! call SetInputUnits(7,)    No frequency units ?
        call SetInputUnits(8, 'PW1')  ! kJ/h
        call SetInputUnits(9, 'PW1')  ! kJ/h

        call SetOutputUnits(1,'TE1')  ! °C
        call SetOutputUnits(2,'DM1')  ! -
        call SetOutputUnits(3,'PC1')  ! %
        call SetOutputUnits(4,'MF1')  ! kg/h
        call SetOutputUnits(5,'PR4')  ! atm
        call SetOutputUnits(6,'PW1')  ! kJ/h
        call SetOutputUnits(7,'PW1')  ! kJ/h
        call SetOutputUnits(8,'PW1')  ! kJ/h
        call SetOutputUnits(9,'PW1')  ! kJ/h
        call SetOutputUnits(10,'PW1')  ! kJ/h
        call SetOutputUnits(11,'DM1')  ! -
        call SetOutputUnits(12,'DM1')  ! -
        call SetOutputUnits(13,'PW1')  ! kJ/h
        call SetOutputUnits(14,'PW1')  ! kJ/h
        call SetOutputUnits(15,'PW1')  ! kJ/h
        call SetOutputUnits(16,'TE1')  ! °C
        call SetOutputUnits(17,'MF1')  ! kg/h

  	    return

    endif
    
    ! Start of the first timestep: no iterations, outputs initial conditions
    if (getIsStartTime()) then
        
        call ReadParameters()

        Tr = GetInputValue(1)
        wr = GetInputValue(2)
        RHr = GetInputValue(3)
        amfr = GetInputValue(4)
        pr = GetInputValue(5)
        Toa = GetInputValue(6)
        freq = GetInputValue(7)
        PfanI = GetInputValue(8)
        PfanO = GetInputValue(9)

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
	    call SetOutputValue(10, 0.0_wp)  ! Total power consumption
	    call SetOutputValue(11, 0.0_wp)  ! COP
	    call SetOutputValue(12, 0.0_wp)  ! EER
	    call SetOutputValue(13, 0.0_wp)  ! Indoor fan power
	    call SetOutputValue(14, 0.0_wp)  ! Outdoor fan power
	    call SetOutputValue(15, 0.0_wp)  ! Compressor power
	    call SetOutputValue(16, 0.0_wp)  ! Condensate temperature
	    call SetOutputValue(17, 0.0_wp)  ! Condensate flow rate

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
        yHum = getParameterValue(1)
        PelRated = getParameterValue(2)
        QcsRated = getParameterValue(3)
        QclRated = getParameterValue(4)
        amfrRated = getParameterValue(5)
        freqRated = getParameterValue(6)
        LUcool = getParameterValue(7)
    end subroutine ReadParameters
    
    
    subroutine GetInputValues
        Tr = GetInputValue(1)
        wr = GetInputValue(2)
        RHr = GetInputValue(3)
        amfr = GetInputValue(4)
        pr = GetInputValue(5)
        Toa = GetInputValue(6)
        freq = GetInputValue(7)
        Pfani = GetInputValue(8)
        Pfano = GetInputValue(9)
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
        call SetOutputValue(10, Pel)  ! Total power consumption
        call SetOutputValue(11, COP)  ! COP
        call SetOutputValue(12, EER)  ! EER
        call SetOutputValue(13, PfanI)  ! Indoor fan power
        call SetOutputValue(14, PfanO)  ! Outdoor fan power
        call SetOutputValue(15, Pcomp)  ! Compressor power
        call SetOutputValue(16, Tc)  ! Condensate temperature
        call SetOutputValue(17, cmfr)  ! Condensate flow rate
    end subroutine SetOutputValues
    
    
    subroutine GetTRNSYSvariables
        time = getSimulationTime()
        timestep = getSimulationTimeStep()
        thisUnit = getCurrentUnit()
        thisType = getCurrentType()
    end subroutine GetTRNSYSvariables

end subroutine Type3254
