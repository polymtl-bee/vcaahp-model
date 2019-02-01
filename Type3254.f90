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

subroutine Type3254
!export this subroutine for its use in external DLLs
!DEC$Attributes DLLexport :: Type3254

use, intrinsic :: iso_fortran_env, only : wp=>real64    ! Defines a constant "wp" (working precision) that can be used in real numbers, e.g. 1.0_wp, and sets it to real64 (double precision)

use TrnsysConstants
use TrnsysFunctions

implicit none

real(wp) :: time, timeStep    ! TRNSYS time and timestep
real(wp) :: Tin, wIn, RHin, mDotIn, pIn, Toa, f, Tset, Tm, PfanI, PfanO    ! Inputs
integer :: yControl, yHum, LUcool, nDBin, nWBin, nFlowRates, nDBoa, nf    ! Parameters
real(wp) :: QcRated, QcsRated, PtotRated,mDotInRated, fRated    ! Parameters (rated values)
real(wp) :: Tout, wOut, RHout, mDotOut, pOut    ! Outputs (outlet conditions)
real(wp) :: Qc, Qcs, Qcl, Qrej, Ptot, Pcomp    ! Outputs (heat and power)
real(wp) :: COP, EER, Tcond, mDotCond    ! Outputs (misc)
real(wp) :: x, y    ! arrays with inputs and outputs of interpolation
integer :: nx, ny, nval    ! interpolation parameters
dimension :: x(5), y(3), nval(5)
integer :: thisUnit, thisType    ! unit and type numbers

! Local variables
real(wp) :: psyDat(9), TwbIn, hIn, hOut, Cp, TinK
integer :: psychMode, status

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

  	! Tell the TRNSYS engine how this Type works
  	call SetNumberofParameters(13)
  	call SetNumberofInputs(11)
  	call SetNumberofDerivatives(0)
  	call SetNumberofOutputs(18)
  	call SetIterationMode(1)    !An indicator for the iteration mode (default=1).  Refer to section 8.4.3.5 of the documentation for more details.
  	call SetNumberStoredVariables(0,0)    !The number of static variables that the model wants stored in the global storage array and the number of dynamic variables that the model wants stored in the global storage array
  	call SetNumberofDiscreteControls(0)   !The number of discrete control functions set by this model (a value greater than zero requires the user to use Solver 1: Powell's method)

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

! --- Start time call: not a real time step, there are no iterations at the initial time - output initial conditions ---

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
   !Sample Code: If( PAR1 <= 0.) Call FoundBadParameter(1,'Fatal','The first parameter provided to this model is not acceptable.')

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


   !If Needed, Set the Initial Values of the Static Storage Variables (#,Value)
   !Sample Code: SetStaticArrayValue(1,0.d0)

   !If Needed, Set the Initial Values of the Dynamic Storage Variables (#,Value)
   !Sample Code: Call SetDynamicArrayValueThisIteration(1,20.d0)

   !If Needed, Set the Initial Values of the Discrete Controllers (#,Value)
   !Sample Code for Controller 1 Set to Off: Call SetDesiredDiscreteControlState(1,0)

    return

endif

! --- TRNSYS has detected that parameters must be re-read - indicates another unit of this Type ------------------------

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


endif

! --- Read inputs (for all calls except very first call in simulation) -------------------------------------------------

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
!Sample Code: If( IN1 <= 0.) Call FoundBadInput(1,'Fatal','The first input provided to this model is not acceptable.')

if (ErrorFound()) return
!-----------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------
!    *** PERFORM ALL THE CALCULATION HERE FOR THIS MODEL. ***
!-----------------------------------------------------------------------------------------------------------------------

	!-----------------------------------------------------------------------------------------------------------------------
	!If Needed, Get the Previous Control States if Discrete Controllers are Being Used (#)
	!Sample Code: CONTROL_LAST=getPreviousControlState(1)
	!-----------------------------------------------------------------------------------------------------------------------

	!-----------------------------------------------------------------------------------------------------------------------
	!If Needed, Get the Values from the Global Storage Array for the Static Variables (#)
	!Sample Code: STATIC1=getStaticArrayValue(1)
	!-----------------------------------------------------------------------------------------------------------------------

	!-----------------------------------------------------------------------------------------------------------------------
	!If Needed, Get the Initial Values of the Dynamic Variables from the Global Storage Array (#)
	!Sample Code: T_INITIAL_1=getDynamicArrayValueLastTimestep(1)
	!-----------------------------------------------------------------------------------------------------------------------

! Inlet air state
psyDat(1) = pIn
psyDat(2) = Tin
psyDat(4) = RHin/100.0_wp
psyDat(6) = wIn
if (yHum == 1) then
    psychMode = 4
else
    psychMode = 2
endif
call MoistAirProperties(thisUnit, thisType, 1, psychMode, 0, psyDat, 1, status) ! unit, type, si units used, psych inputs, Twb not computed, inputs, warnings mgmt, warning occurrences 
pIn = psyDat(1)
Tin = psyDat(2)
TwbIn = psydat(3)
RHin = psyDat(4)    ! RHin between 0 and 1 (not 0 and 100)
wIn = psyDat(6)
hIn = psyDat(7)

! Find cooling performance from data file
nx = 5
nval(5) = nDBin
nval(4) = nDBoa
nval(3) = nWBin
nval(2) = nf
nval(1) = nFlowRates
x(1) = mDotIn
x(2) = f
x(3) = TwbIn
x(4) = Toa
x(5) = Tin
ny = 3
call InterpolateData(LUcool,nx,nval,ny,x,y)
if (ErrorFound()) return

Qc = QcRated * y(1)
Qcs = QcsRated * y(2)
Ptot = PtotRated * y(3)

! Outlet air state
mDotOut = mDotIn    ! Dry air mass conservation
pOut = pIn    ! Fan pressure drop neglected

! Moist air state
if (Qc < Qcs) then
    Qc = Qcs
    ! Add warning
endif

if (mDotOut /= 0) then
    hOut = hIn - Qc/mDotOut
    wOut = wIn    ! useful when the following if clause is not true
    if (Qc > Qcs) then    ! nonzero latent heat
        psyDat(1) = pIn
        psyDat(2) = Tin
        psyDat(7) = hIn - (Qc - Qcs)/mDotIn    ! enthalpy of the state (Tin, wOut)
        call MoistAirProperties(thisUnit, thisType, 1, 5, 0, psyDat, 1, status)    ! dry-bulb and enthalpy as inputs
        if (ErrorFound()) return
        wOut = psyDat(6)
    endif
else
    hOut = hIn
    wOut = wIn
endif

psyDat(1) = pOut
psyDat(6) = wOut
psyDat(7) = hOut
call MoistAirProperties(thisUnit, thisType, 1, 7, 0, psyDat, 1, status)    ! humidity ratio and enthalpy as inputs
if (ErrorFound()) return
pOut = psyDat(1)
Tout = psydat(2)
RHout = psydat(4)
wOut = psydat(6)
hOut = psydat(7)

! Re-calculate heat transfer whose value is modified if saturation occurs
Qc = mDotOut * (hIn - hOut)    ! Total cooling rate
Qcs = mDotOut * Cp * (Tin - Tout)    ! Sensible cooling rate
Qcl = Qc - Qcs    ! Latent cooling rate
Qrej = Qc + Ptot    ! Heat rejection
Pcomp = Ptot - PfanI - PfanO    ! Compressor power
if (Ptot /= 0.) then
    COP = Qc/Ptot
else
    COP = 0.0_wp
endif
EER = 3.413_wp * COP
Tcond = Tout
mDotCond = mDotOut * (wIn - wOut)    ! Condensate flow rate - water balance


!Set the Outputs from this Model (#,Value)
Call SetOutputValue(1, Tout) ! Outlet air temperature
Call SetOutputValue(2, wOut) ! Outlet air humidity ratio
Call SetOutputValue(3, RHout*100.0_wp) ! Outlet air % RH
Call SetOutputValue(4, mDotOut) ! Outlet air flow rate
Call SetOutputValue(5, pOut) ! Outlet air pressure
Call SetOutputValue(6, Qc) ! Total cooling rate
Call SetOutputValue(7, Qcs) ! Sensible cooling rate
Call SetOutputValue(8, Qcl) ! Latent cooling rate
Call SetOutputValue(9, Qrej) ! Heat rejection rate
Call SetOutputValue(10, Ptot) ! Total power consumption
Call SetOutputValue(11, COP) ! COP
Call SetOutputValue(12, EER) ! EER
Call SetOutputValue(13, PfanI) ! Indoor fan power
Call SetOutputValue(14, PfanO) ! Outdoor fan power
Call SetOutputValue(15, Pcomp) ! Compressor power
Call SetOutputValue(16, Tcond) ! Condensate temperature
Call SetOutputValue(17, mDotCond) ! Condensate flow rate
Call SetOutputValue(18, f) ! Compressor frequency

return

! contains
! subroutine ...
! ...
! end
! types 155, 28

end subroutine Type3254
