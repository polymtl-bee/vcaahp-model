!TRNSYS Type3254: Variable capacity air-source heat pump with performance files
! ---------------------------------------------------------------------------------------------------------------------
!
! This routine implements an air-source heat pump with variable compressor speed.
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
!
!
! Revision history
! ---------------------------------------------------------------------------------------------------------------------
!
! Created 2018-06-07 - Gregor Strugala, École Polytechnique de Montréal
! ---------------------------------------------------------------------------------------------------------------------

subroutine Type3254
!export this subroutine for its use in external DLLs
!DEC$Attributes DLLexport :: Type3254

use, intrinsic :: iso_fortran_env, only : wp=>real64    ! Defines a constant "wp" (working precision) that can be used in real numbers, e.g. 1.0_wp, and sets it to real64 (double precision)

use TrnsysConstants
use TrnsysFunctions

implicit None

! Local variables
real(wp) :: time, timeStep    ! TRNSYS time and timestep
real(wp) :: Tin, wIn, RHin, mDotIn, pIn, Toa, f, Tset, Tm, Pfani, Pfano    ! Inputs
real(wp) :: yControl, yHum, LUcool, nDBin, nWBin, nFlowRates, nDBoa nf    ! Parameters
real(wp) :: QcRated, QcsRated, PtotRated,mDotInRated, fRated    ! Parameters (rated values)
real(wp) :: Tout, wOut, RHout, mDotOut, pOut    ! Outputs (outlet conditions)
real(wp) :: Qc, Qcs, Qcl, Qrej, Ptot, Pcomp    ! Outputs (heat and power)
real(wp) :: COP, EER, Tcond, mDotCond    ! Outputs (misc)
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

  	! Tell the TRNSYS engine how this Type works
  	call SetNumberofParameters(13)
  	call SetNumberofInputs(11)
  	call SetNumberofDerivatives(0)
  	call SetNumberofOutputs(18)
  	call SetIterationMode(1)    !An indicator for the iteration mode (default=1).  Refer to section 8.4.3.5 of the documentation for more details.
  	call SetNumberStoredVariables(0,0)    !The number of static variables that the model wants stored in the global storage array and the number of dynamic variables that the model wants stored in the global storage array
  	call SetNumberofDiscreteControls(0)   !The number of discrete control functions set by this model (a value greater than zero requires the user to use Solver 1: Powell's method)

  	return

endif

! --- Do All of the First Timestep Manipulations Here - There Are No Iterations at the Intial Time ---------------------

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

   !Set the Initial Values of the Outputs (#,Value)
		call SetOutputValue(1, 0) ! Outlet air temperature
		call SetOutputValue(2, 0) ! Outlet air humidity ratio
		call SetOutputValue(3, 0) ! Outlet air % RH
		call SetOutputValue(4, 0) ! Outlet air flow rate
		call SetOutputValue(5, 0) ! Outlet air pressure
		call SetOutputValue(6, 0) ! Total cooling rate
		call SetOutputValue(7, 0) ! Sensible cooling rate
		call SetOutputValue(8, 0) ! Latent cooling rate
		call SetOutputValue(9, 0) ! Heat rejection rate
		call SetOutputValue(10, 0) ! Total power consumption
		call SetOutputValue(11, 0) ! COP
		call SetOutputValue(12, 0) ! EER
		call SetOutputValue(13, 0) ! Indoor fan power
		call SetOutputValue(14, 0) ! Outdoor fan power
		call SetOutputValue(15, 0) ! Compressor power
		call SetOutputValue(16, 0) ! Condensate temperature
		call SetOutputValue(17, 0) ! Condensate flow rate
		call SetOutputValue(18, 0) ! Compressor frequency


   !If Needed, Set the Initial Values of the Static Storage Variables (#,Value)
   !Sample Code: SetStaticArrayValue(1,0.d0)

   !If Needed, Set the Initial Values of the Dynamic Storage Variables (#,Value)
   !Sample Code: Call SetDynamicArrayValueThisIteration(1,20.d0)

   !If Needed, Set the Initial Values of the Discrete Controllers (#,Value)
   !Sample Code for Controller 1 Set to Off: Call SetDesiredDiscreteControlState(1,0)

		return

endif
!-----------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------
!ReRead the Parameters if Another Unit of This Type Has Been Called Last
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
!-----------------------------------------------------------------------------------------------------------------------

!Read the Inputs
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

if(ErrorFound()) return
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

	!-----------------------------------------------------------------------------------------------------------------------
	!Perform All of the Calculations Here to Set the Outputs from the Model Based on the Inputs

	!Sample Code: OUT1=IN1+PAR1

	!If the model requires the solution of numerical derivatives, set these derivatives and get the current solution
	!Sample Code: T1=getNumericalSolution(1)
	!Sample Code: T2=getNumericalSolution(2)
	!Sample Code: DTDT1=3.*T2+7.*T1-15.
	!Sample Code: DTDT2=-2.*T1+11.*T2+21.
	!Sample Code: Call SetNumericalDerivative(1,DTDT1)
	!Sample Code: Call SetNumericalDerivative(2,DTDT2)

!-----------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------
!Set the Outputs from this Model (#,Value)
		Call SetOutputValue(1, 0) ! Outlet air temperature
		Call SetOutputValue(2, 0) ! Outlet air humidity ratio
		Call SetOutputValue(3, 0) ! Outlet air % RH
		Call SetOutputValue(4, 0) ! Outlet air flow rate
		Call SetOutputValue(5, 0) ! Outlet air pressure
		Call SetOutputValue(6, 0) ! Total cooling rate
		Call SetOutputValue(7, 0) ! Sensible cooling rate
		Call SetOutputValue(8, 0) ! Latent cooling rate
		Call SetOutputValue(9, 0) ! Heat rejection rate
		Call SetOutputValue(10, 0) ! Total power consumption
		Call SetOutputValue(11, 0) ! COP
		Call SetOutputValue(12, 0) ! EER
		Call SetOutputValue(13, 0) ! Indoor fan power
		Call SetOutputValue(14, 0) ! Outdoor fan power
		Call SetOutputValue(15, 0) ! Compressor power
		Call SetOutputValue(16, 0) ! Condensate temperature
		Call SetOutputValue(17, 0) ! Condensate flow rate
		Call SetOutputValue(18, 0) ! Compressor frequency

      return
end subroutine Type3254
