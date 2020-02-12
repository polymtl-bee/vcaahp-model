!TRNSYS Type3254: Variable capacity air-source heat pump with performance files
! ---------------------------------------------------------------------------------------------------------------------
!
! This routine implements an air-source heat pump with variable speed compressor.
!
!
! Inputs
! ---------------------------------------------------------------------------------------------------------------------
! Nb | Variable     | Description                                                   | Input  Units    | Internal Units
! ---|--------------|---------------------------------------------------------------|-----------------|----------------
!  1 | Tin          | Inlet air temperature                                         | °C              | °C
!  2 | wIn          | Inlet air humidity ratio                                      | -               | -
!  3 | RHin         | Inlet air relative humidity                                   | % (base 100)    | -
!  4 | mDotIn       | Inlet air flow rate                                           | kg/h            | kg/h
!  5 | pIn          | Inlet air pressure                                            | atm             | atm
!  6 | Toa          | Outside air temperature                                       | °C              | °C
!  7 | f            | Compressor frequency                                          | 1/s             | 1/s
!  8 | Tset         | Indoor air temperature setpoint                               | °C              | °C
!  9 | Tm           | Indoor air temperature measurment                             | °C              | °C
! 10 | Pfani        | Indoor fan power                                              | kJ/h            | kJ/h
! 11 | Pfano        | Outdoor fan power                                             | kJ/h            | kJ/h
!
!
! Parameters
! ---------------------------------------------------------------------------------------------------------------------
! Nb | Variable     | Description                                                   | Param. Units    | Internal Units
! ---|--------------|---------------------------------------------------------------|-----------------|----------------
!  1 | yControl     | 1 = Internal embedded controller                              | -               | -
!    |              | 2 = External user-defined controller                          | -               | -
!  2 | yHum         | 1 = Humidity ratio as humidity input                          | -               | -
!    |              | 2 = Relative humidity as humidity input                       | -               | -
!  3 | LUcool       | Logical Unit - cooling mode                                   | -               | -
!  4 | nDBin        | Number of indoor dry-bulb temperatures                        | -               | -
!  5 | nWBin        | Number of indoor wet-bulb temperatures                        | -               | -
!  6 | nFlowRates   | Number of air flow rates                                      | -               | -
!  7 | nDBoa        | Number of outdoor dry-bulb temperatures                       | -               | -
!  8 | nf           | Number of frequencies                                         | -               | -
!  9 | QcRated      | Rated cooling capacity                                        | kJ/h            | kJ/h
! 10 | QcsRated     | Rated sensible cooling capacity                               | kJ/h            | kJ/h
! 11 | PtotRated    | Rated total cooling power                                     | kJ/h            | kJ/h
! 12 | mDotInRated  | Rated inlet air flow rate                                     | kg/h            | kg/h
! 13 | fRated       | Rated frequency                                               | 1/s             | 1/s
!
!
! Outputs
! ---------------------------------------------------------------------------------------------------------------------
! Nb | Variable     | Description                                                   | Output  Units   | Internal Units
! ---|--------------|---------------------------------------------------------------|-----------------|----------------
!  1 | Tout         | Outlet air temperature                                        | °C              | °C
!  2 | wOut         | Outlet air humidity ratio                                     | -               | -
!  3 | RHout        | Outlet air % RH                                               | % (base 100)    | % (base 100)
!  4 | mDotOut      | Outlet air flow rate                                          | kg/h            | kg/h
!  5 | pOut         | Outlet air pressure                                           | atm             | atm
!  6 | Qc           | Total cooling rate                                            | kJ/h            | kJ/h
!  7 | Qcs          | Sensible cooling rate                                         | kJ/h            | kJ/h
!  8 | Qcl          | Latent cooling rate                                           | kJ/h            | kJ/h
!  9 | Qrej         | Heat rejection rate                                           | kJ/h            | kJ/h
! 10 | Ptot         | Total power consumption                                       | kJ/h            | kJ/h
! 11 | COP          | Coefficient of performance                                    | -               | -
! 12 | EER          | Energy efficiency rating                                      | -               | -
! 13 | Pfani        | Indoor fan power                                              | kJ/h            | kJ/h
! 14 | Pfano        | Outdoor fan power                                             | kJ/h            | kJ/h
! 15 | Pcomp        | Compressor power                                              | kJ/h            | kJ/h
! 16 | Tcond        | Condensate temperature                                        | °C              | °C
! 17 | mDotCond     | Condensate flow rate                                          | kg/h            | kg/h
! 18 | f            | Compressor frequency                                          | 1/s             | 1/s

module Type3254Data

use, intrinsic :: iso_fortran_env, only : wp=>real64    ! Defines a constant "wp" (working precision) that can be used in real numbers, e.g. 1.0_wp, and sets it to real64 (double precision)
use TrnsysConstants
use TrnsysFunctions
implicit none

integer, parameter :: nMaxMap = 20

type Type3254DataStruct
    ! Parameters
    real(wp) :: QcRated, QcsRated, PtotRated, mDotInRated, fRated   ! rated values
    integer :: nTr, nRHr, nTo, nmfr, nfreq   ! number of entries for each variable
    real(wp), allocatable :: entries(:, :)
    real(wp), allocatable :: extents(:)

    ! Performance matrices (Choose better names?)
    real(wp), allocatable :: PelMap(:, :, :, :, :)
    real(wp), allocatable :: QcsMap(:, :, :, :, :)
    real(wp), allocatable :: QclMap(:, :, :, :, :)

end type Type3254DataStruct

type(Type3254DataStruct), allocatable, save :: storedData(:)

end module Type3254Data


subroutine Type3254
!export this subroutine for its use in external DLLs
!DEC$Attributes DLLexport :: Type3254

use, intrinsic :: iso_fortran_env, only : wp=>real64    ! Defines a constant "wp" (working precision) that can be used in real numbers, e.g. 1.0_wp, and sets it to real64 (double precision)

use TrnsysConstants
use TrnsysFunctions
use Type3254Data

implicit none

real(wp) :: time, timeStep    ! TRNSYS time and timestep
real(wp) :: Tin, wIn, RHin, mDotIn, pIn, Toa, f, Tset, Tm, PfanI, PfanO    ! Inputs
integer :: yControl, yHum, LUcool, nDBin, nWBin, nFlowRates, nDBoa, nf    ! Parameters
real(wp) :: QcRated, QcsRated, PtotRated,mDotInRated, fRated    ! Parameters (rated values)
real(wp) :: Tout, wOut, RHout, mDotOut, pOut    ! Outputs (outlet conditions)
real(wp) :: Qc, Qcs, Qcl, Qrej, Pel, Ptot, Pcomp    ! Outputs (heat and power)
!real(wp) :: Pel(nMaxMap, nMaxMap), Qevs(nMaxMap, nMaxMap), Qevl(nMaxMap, nMaxMap)
real(wp) :: COP, EER, Tcond, mDotCond    ! Outputs (misc)
real(wp) :: x, y    ! arrays with inputs and outputs of interpolation
integer :: nx, ny, nval    ! interpolation parameters
dimension :: x(5), y(3), nval(5)
integer :: thisUnit, thisType    ! unit and type numbers
integer, parameter :: N = 5   ! Number of interpolation variables
integer :: extents(N)  ! Number of values for each interpolation variable

! Local variables
!real(wp) :: psyDat(9), TwbIn, hIn, hx, hOut
!integer :: psychMode, status
integer, parameter :: NumOfInstances = 1
integer :: thisInstanceNo = 1   ! temporary, should use a kernel function to get the actual unit number.
real(wp) :: filler(N), PMvalue

! Interpolation variables
character (len=maxPathLength) :: permapPath
character (len=maxMessageLength) :: msg
logical :: permapFileFound = .false.
integer :: nTr, nRHr, nTo, nmfr, nfreq   ! number of entries for each variable
integer :: i, j, line_count = 1, prev_line
real(wp), allocatable, dimension(:) :: TrValues,  RHrValues, ToValues, mfrValues, freqValues
integer, allocatable :: entries(:, :)
real(wp), allocatable :: PelMap(:, :, :, :, :)   ! N dimensions
real(wp), allocatable :: QcsMap(:, :, :, :, :)
real(wp), allocatable :: QclMap(:, :, :, :, :)
integer :: idx(N) = 1

! Set the version number for this Type
if (GetIsVersionSigningTime()) then
    call SetTypeVersion(18)
    return  ! We are done for this call
endif

call GetTRNSYSvariables()
call ExecuteSpecialCases()
call GetInputValues()
if (ErrorFound()) return

nTr = storedData(thisInstanceNo)%nTr
nRHr = storedData(thisInstanceNo)%nRHr
nTo = storedData(thisInstanceNo)%nTo
nmfr = storedData(thisInstanceNo)%nmfr
nfreq = storedData(thisInstanceNo)%nfreq
extents = storedData(thisInstanceNo)%extents
allocate(entries(N, maxval(extents)))
allocate(PelMap(nTr, nRHr, nTo, nmfr, nfreq))
allocate(QcsMap(nTr, nRHr, nTo, nmfr, nfreq))
allocate(QclMap(nTr, nRHr, nTo, nmfr, nfreq))
entries = storedData(thisInstanceNo)%entries
PelMap = storedData(thisInstanceNo)%PelMap
QcsMap = storedData(thisInstanceNo)%QcsMap
QclMap = storedData(thisInstanceNo)%QclMap

!! Inlet air state
!psyDat(1) = pIn
!psyDat(2) = Tin
!psyDat(4) = RHin/100.0_wp
!psyDat(6) = wIn
!if (yHum == 1) then
!    psychMode = 4
!else
!    psychMode = 2
!endif
!call MoistAirProperties(thisUnit, thisType, 1, psychMode, 0, psyDat, 1, status) ! unit, type, si units used, psych inputs, Twb not computed, inputs, warnings mgmt, warning occurrences 
!pIn = psyDat(1)
!Tin = psyDat(2)
!TwbIn = psydat(3)
!RHin = psyDat(4)    ! RHin between 0 and 1 (not 0 and 100)
!wIn = psyDat(6)
!hIn = psyDat(7)
!
!! Find cooling performance from data file
!nx = 5
!nval(5) = nDBin
!nval(4) = nDBoa
!nval(3) = nWBin
!nval(2) = nf
!nval(1) = nFlowRates
!x(1) = mDotIn/mDotInRated
!x(2) = f/fRated
!x(3) = TwbIn
!x(4) = Toa
!x(5) = Tin
!ny = 3
!call InterpolateData(LUcool,nx,nval,ny,x,y)
!if (ErrorFound()) return
!
!Qc = QcRated * y(1)
!Qcs = QcsRated * y(2)
!Ptot = PtotRated * y(3)
!
!! Outlet air state
!mDotOut = mDotIn    ! Dry air mass conservation
!pOut = pIn    ! Fan pressure drop neglected
!
!! Moist air state
!if (Qc < Qcs) then
!    Qc = Qcs
!    ! Add warning
!endif
!
!if (mDotOut /= 0) then
!    hOut = hIn - Qc/mDotOut
!    wOut = wIn    ! useful when the following if clause is not true
!    hx = hIn      ! same
!    if (Qc > Qcs) then    ! nonzero latent heat
!        psyDat(1) = pIn
!        psyDat(2) = Tin
!        hx = hIn - (Qc - Qcs)/mDotIn
!        psyDat(7) = hx    ! enthalpy of the state (Tin, wOut)
!        call MoistAirProperties(thisUnit, thisType, 1, 5, 0, psyDat, 1, status)    ! dry-bulb and enthalpy as inputs
!        if (ErrorFound()) return
!        wOut = psyDat(6)
!    endif
!else
!    hOut = hIn
!    wOut = wIn
!    hx = hIn
!endif
!
!psyDat(1) = pOut
!psyDat(6) = wOut
!psyDat(7) = hOut
!call MoistAirProperties(thisUnit, thisType, 1, 7, 0, psyDat, 1, status)    ! humidity ratio and enthalpy as inputs
!if (ErrorFound()) return
!pOut = psyDat(1)
!Tout = psydat(2)
!RHout = psydat(4)
!wOut = psydat(6)
!hOut = psydat(7)
!
!! Re-calculate heat transfer whose value is modified if saturation occurs
!Qc = mDotOut * (hIn - hOut)    ! Total cooling rate
!Qcs = mDotOut * (hx - hOut)    ! Sensible cooling rate
!Qcl = Qc - Qcs    ! Latent cooling rate
!Qrej = Qc + Ptot    ! Heat rejection
!Pcomp = Ptot - PfanI - PfanO    ! Compressor power
!if (Ptot /= 0.) then
!    COP = Qc/Ptot
!else
!    COP = 0.0_wp
!endif
!EER = 3.413_wp * COP
!Tcond = Tout
!mDotCond = mDotOut * (wIn - wOut)    ! Condensate flow rate - water balance

call SetOutputValues()

return

    contains
    
    subroutine ReadPermap
    
    permapPath = GetLUfileName(LUcool)
    
    inquire(file=trim(permapPath), exist=permapFileFound)
    
    if ( .not. permapFileFound ) then
        write(msg,'("""",a,"""")') trim(permapPath)
        msg = "Could not find the specified performance map file. Searched for: " // trim(msg)
        call Messages(-1, msg, 'fatal', thisUnit, thisType)
        return
    end if
    
    open(LUcool, file=permapPath, status='old')
    
        do i = 1, 7   ! Skip 7 lines
            read(LUcool, *)
        enddo
    read(LUcool, *) nTr   ! Read number of Tr values
        read(LUcool, *)   ! Skip line
    allocate(TrValues(nTr))   ! Read Tr values
    read(LUcool, *) (TrValues(i), i = 1, nTr)
        read(LUcool, *)   ! Skip line
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
    read(LUcool, *) nmfr
        read(LUcool, *)
    allocate(mfrValues(nmfr))
    read(LUcool, *) (mfrValues(i), i = 1, nmfr)
        read(LUcool, *)
    read(LUcool, *) nfreq
        read(LUcool, *)
    allocate(freqValues(nfreq))
    read(LUcool, *) (freqValues(i), i = 1, nfreq)
        do i = 1, 4   ! Skip 4 lines
            read(LUcool, *)
        end do
    allocate(PelMap(nTr, nRHr, nTo, nmfr, nfreq))
    allocate(QcsMap(nTr, nRHr, nTo, nmfr, nfreq))
    allocate(QclMap(nTr, nRHr, nTo, nmfr, nfreq))
    extents = shape(PelMap)
    allocate(entries(N, maxval(extents)))
    
    read(LUcool, *) (filler(i), i = 1, N), Pel, Qcs, Qcl
    call SetPMvalue(PelMap, idx, Pel)
    call SetPMvalue(QcsMap, idx, Qcs)
    call SetPMvalue(QclMap, idx, Qcl)
    do while (sum(idx) < sum(extents))
        prev_line = line_count
        i = 1
        do while (prev_line == line_count)
            if (idx(i) == extents(i)) then
                i = i + 1
                if (i > 1) idx(i-1) = 1
            else
                idx(i) = idx(i) + 1
                line_count = line_count + 1
                read(LUcool, *) (filler(j), j = 1, N), Pel, Qcs, Qcl
                call SetPMvalue(PelMap, idx, Pel)
                call SetPMvalue(QcsMap, idx, Qcs)
                call SetPMvalue(QclMap, idx, Qcl)
            end if
        end do
    end do
    
    close(LUcool)
    
    ! make a check ?
    
    storedData(thisInstanceNo)%nTr = nTr
    storedData(thisInstanceNo)%nRHr = nRHr
    storedData(thisInstanceNo)%nTo = nTo
    storedData(thisInstanceNo)%nmfr = nmfr
    storedData(thisInstanceNo)%nfreq = nfreq
    allocate(storedData(NumOfInstances)%extents(N))
    allocate(storedData(NumOfInstances)%entries(N, maxval(extents)))
    allocate(storedData(NumOfInstances)%PelMap(nTr, nRHr, nTo, nmfr, nfreq))
    allocate(storedData(NumOfInstances)%QcsMap(nTr, nRHr, nTo, nmfr, nfreq))
    allocate(storedData(NumOfInstances)%QclMap(nTr, nRHr, nTo, nmfr, nfreq))
    storedData(thisInstanceNo)%extents = extents
    storedData(thisInstanceNo)%entries = entries
    storedData(thisInstanceNo)%PelMap = PelMap
    storedData(thisInstanceNo)%QcsMap = QcsMap
    storedData(thisInstanceNo)%QclMap = QclMap
    deallocate(entries)
    deallocate(PelMap)
    deallocate(QcsMap)
    deallocate(QclMap)
    
    end subroutine ReadPermap
    
    
    function GetPMvalue(array, idx)
        integer, intent(in) :: idx(:)
        integer :: i, array_idx
        real(wp) :: array(product(extents))
        real(wp) :: GetPMvalue
        array_idx = idx(1)
        do i = 2, size(idx)
            array_idx = array_idx + product(extents(:i-1)) * (idx(i) - 1)
        end do
        GetPMvalue = array(array_idx)
    end function GetPMvalue
    
    
    subroutine SetPMvalue(array, idx, value)
        real(wp) :: array(product(extents))
        real(wp), intent(in) :: value
        integer, intent(in) :: idx(:)
        integer :: i, array_idx
        array_idx = idx(1)
        do i = 2, size(idx)
            array_idx = array_idx + product(extents(:i-1)) * (idx(i) - 1)
        end do
        array(array_idx) = value
    end subroutine SetPMvalue
    
    
    subroutine ExecuteSpecialCases
    
    ! All the stuff that must be done once at the beginning
    if(getIsFirstCallofSimulation()) then

  	    ! Tell the TRNSYS engine how this Type works
  	    call SetNumberofParameters(13)
  	    call SetNumberofInputs(11)
  	    call SetNumberofDerivatives(0)
  	    call SetNumberofOutputs(18)
  	    call SetIterationMode(1)    !An indicator for the iteration mode (default=1).  Refer to section 8.4.3.5 of the documentation for more details.
  	    call SetNumberStoredVariables(0,0)    !The number of static variables that the model wants stored in the global storage array and the number of dynamic variables that the model wants stored in the global storage array
  	    call SetNumberofDiscreteControls(0)   !The number of discrete control functions set by this model (a value greater than zero requires the user to use Solver 1: Powell's method)
        
        ! Allocate stored data structure
        if (.not. allocated(storedData)) then
            allocate(storedData(numOfInstances))
        endif
        
        yControl = getParameterValue(1)
        yHum = getParameterValue(2)
        LUcool = getParameterValue(3)
        nDBin = getParameterValue(4)
        nWBin = getParameterValue(5)
        nFlowRates = getParameterValue(6)
        nDBoa = getParameterValue(7)
        nf = getParameterValue(8)
        QcRated = getParameterValue(9)
        QcsRated = getParameterValue(10)
        PtotRated = getParameterValue(11)
        mDotInRated = getParameterValue(12)
        fRated = getParameterValue(13)
        
        call ReadPermap()

        ! Set units (optional)
        call SetInputUnits(1,'TE1')    ! °C
        call SetInputUnits(2,'DM1')    ! -
        call SetInputUnits(3,'PC1')    ! %
        call SetInputUnits(4,'MF1')    ! kg/h
        call SetInputUnits(5,'PR4')    ! atm
        call SetInputUnits(6,'TE1')    ! °C
        ! call SetInputUnits(7,)    No frequency units ?
        call SetInputUnits(8,'TE1')    ! °C
        call SetInputUnits(9,'TE1')    ! °C
        call SetInputUnits(10,'PW1')    ! kJ/h
        call SetInputUnits(11,'PW1')    ! kJ/h

        call SetOutputUnits(1,'TE1')    ! °C
        call SetOutputUnits(2,'DM1')    ! -
        call SetOutputUnits(3,'PC1')    ! %
        call SetOutputUnits(4,'MF1')    ! kg/h
        call SetOutputUnits(5,'PR4')    ! atm
        call SetOutputUnits(6,'PW1')    ! kJ/h
        call SetOutputUnits(7,'PW1')    ! kJ/h
        call SetOutputUnits(8,'PW1')    ! kJ/h
        call SetOutputUnits(9,'PW1')    ! kJ/h
        call SetOutputUnits(10,'PW1')    ! kJ/h
        call SetOutputUnits(11,'DM1')    ! -
        call SetOutputUnits(12,'DM1')    ! -
        call SetOutputUnits(13,'PW1')    ! kJ/h
        call SetOutputUnits(14,'PW1')    ! kJ/h
        call SetOutputUnits(15,'PW1')    ! kJ/h
        call SetOutputUnits(16,'TE1')    ! °C
        call SetOutputUnits(17,'MF1')    ! kg/h
        ! call SetInputUnits(18,)    No frequency units ?

  	    return

    endif
    
    ! Start of the first timestep: no iterations, outputs initial conditions
    if (getIsStartTime()) then
        
        yControl = getParameterValue(1)
        yHum = getParameterValue(2)
        LUcool = getParameterValue(3)
        nDBin = getParameterValue(4)
        nWBin = getParameterValue(5)
        nFlowRates = getParameterValue(6)
        nDBoa = getParameterValue(7)
        nf = getParameterValue(8)
        QcRated = getParameterValue(9)
        QcsRated = getParameterValue(10)
        PtotRated = getParameterValue(11)
        mDotInRated = getParameterValue(12)
        fRated = getParameterValue(13)


        Tin = GetInputValue(1)
        wIn = GetInputValue(2)
        RHin = GetInputValue(3)
        mDotIn = GetInputValue(4)
        pIn = GetInputValue(5)
        Toa = GetInputValue(6)
        f = GetInputValue(7)
        Tset = GetInputValue(8)
        Tm = GetInputValue(9)
        Pfani = GetInputValue(10)
        Pfano = GetInputValue(11)

       !Check the Parameters for Problems (#,ErrorType,Text)
       !Sample Code: If( PAR1 <= 0.) call FoundBadParameter(1,'Fatal','The first parameter provided to this model is not acceptable.')

        ! Set outputs to zeros at initial time
	    call SetOutputValue(1, 0.0_wp) ! Outlet air temperature
	    call SetOutputValue(2, 0.0_wp) ! Outlet air humidity ratio
	    call SetOutputValue(3, 0.0_wp) ! Outlet air % RH
	    call SetOutputValue(4, 0.0_wp) ! Outlet air flow rate
	    call SetOutputValue(5, 0.0_wp) ! Outlet air pressure
	    call SetOutputValue(6, 0.0_wp) ! Total cooling rate
	    call SetOutputValue(7, 0.0_wp) ! Sensible cooling rate
	    call SetOutputValue(8, 0.0_wp) ! Latent cooling rate
	    call SetOutputValue(9, 0.0_wp) ! Heat rejection rate
	    call SetOutputValue(10, 0.0_wp) ! Total power consumption
	    call SetOutputValue(11, 0.0_wp) ! COP
	    call SetOutputValue(12, 0.0_wp) ! EER
	    call SetOutputValue(13, 0.0_wp) ! Indoor fan power
	    call SetOutputValue(14, 0.0_wp) ! Outdoor fan power
	    call SetOutputValue(15, 0.0_wp) ! Compressor power
	    call SetOutputValue(16, 0.0_wp) ! Condensate temperature
	    call SetOutputValue(17, 0.0_wp) ! Condensate flow rate
	    call SetOutputValue(18, 0.0_wp) ! Compressor frequency

        return

    endif
    
    ! Parameters must be re-read - indicates another unit of this Type
    if(getIsReReadParameters()) then

        yControl = getParameterValue(1)
        yHum = getParameterValue(2)
        LUcool = getParameterValue(3)
        nDBin = getParameterValue(4)
        nWBin = getParameterValue(5)
        nFlowRates = getParameterValue(6)
        nDBoa = getParameterValue(7)
        nf = getParameterValue(8)
        QcRated = getParameterValue(9)
        QcsRated = getParameterValue(10)
        PtotRated = getParameterValue(11)
        mDotInRated = getParameterValue(12)
        fRated = getParameterValue(13)
    
        ! Retrieved Stored data
        !Pel = storedData(thisInstanceNo)%Pel
        !Qevs = storedData(thisInstanceNo)%Qevs
        !Qevl = storedData(thisInstanceNo)%Qevl
    
    endif
    
    ! End of timestep call (after convergence or too many iterations)
    if (GetIsEndOfTimestep()) then
        deallocate(entries)
        deallocate(PelMap)
        deallocate(QcsMap)
        deallocate(QclMap)
        return  ! We are done for this call
    endif
    
    if (GetIsLastCallofSimulation()) then
        ! This Type should not perform any task during the last call
        return  ! We are done for this call
    endif
    
    end subroutine ExecuteSpecialCases
    
    subroutine GetInputValues
        Tin = GetInputValue(1)
        wIn = GetInputValue(2)
        RHin = GetInputValue(3)
        mDotIn = GetInputValue(4)
        pIn = GetInputValue(5)
        Toa = GetInputValue(6)
        f = GetInputValue(7)
        Tset = GetInputValue(8)
        Tm = GetInputValue(9)
        Pfani = GetInputValue(10)
        Pfano = GetInputValue(11)
        !Check the Inputs for Problems (#,ErrorType,Text)
        !Sample Code: If( IN1 <= 0.) call FoundBadInput(1,'Fatal','The first input provided to this model is not acceptable.')
    end subroutine GetInputValues
    
    
    subroutine SetOutputValues
        call SetOutputValue(1, Tout) ! Outlet air temperature
        call SetOutputValue(2, wOut) ! Outlet air humidity ratio
        call SetOutputValue(3, RHout*100.0_wp) ! Outlet air % RH
        call SetOutputValue(4, mDotOut) ! Outlet air flow rate
        call SetOutputValue(5, pOut) ! Outlet air pressure
        call SetOutputValue(6, Qc) ! Total cooling rate
        call SetOutputValue(7, Qcs) ! Sensible cooling rate
        call SetOutputValue(8, Qcl) ! Latent cooling rate
        call SetOutputValue(9, Qrej) ! Heat rejection rate
        call SetOutputValue(10, Ptot) ! Total power consumption
        call SetOutputValue(11, COP) ! COP
        call SetOutputValue(12, EER) ! EER
        call SetOutputValue(13, PfanI) ! Indoor fan power
        call SetOutputValue(14, PfanO) ! Outdoor fan power
        call SetOutputValue(15, Pcomp) ! Compressor power
        call SetOutputValue(16, Tcond) ! Condensate temperature
        call SetOutputValue(17, mDotCond) ! Condensate flow rate
        call SetOutputValue(18, f) ! Compressor frequency
    end subroutine SetOutputValues
    
    subroutine GetTRNSYSvariables
        time = getSimulationTime()
        timestep = getSimulationTimeStep()
        thisUnit = getCurrentUnit()
        thisType = getCurrentType()
    end subroutine GetTRNSYSvariables

end subroutine Type3254
