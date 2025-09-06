!**********************************************************************************************************************************
! This file contains modifications to code originally developed by the National Renewable Energy Laboratory (NREL) 
! as part of FAST's Controls and Electrical Drive Module, "ServoDyn".
!**********************************************************************************************************************************
! ORIGINAL WORK
! Copyright (C) 2013-2016  National Renewable Energy Laboratory
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
!**********************************************************************************************************************************
! MODIFICATIONS
! Copyright (C) 2021 NUMATH https://numath.dmae.upm.es
!          
! The modifications in this file are licensed under the MIT License.          
!**********************************************************************************************************************************

MODULE ControllerInterface

   USE SMConstants
   USE MPI_Process_Info
   
   USE, INTRINSIC :: ISO_C_Binding

   IMPLICIT NONE

   ! Module constant parameters
   INTEGER, PARAMETER :: RPC = SELECTED_REAL_KIND(4)
   INTEGER, PARAMETER :: IPC = SELECTED_INT_KIND(9)
   INTEGER, PARAMETER :: CLEN = 1024

   ! DISCON interface constant parameters
   INTEGER, PARAMETER :: R = 165                                   ! Start index of the generator speed look-up table in SWAP array (Bladed version 4.3 and above)
   INTEGER, PARAMETER :: MaxLoggingChannels = 0                    ! Maximum number of logging channels

   INTEGER, PARAMETER :: DISCON_STATUS_FINALISING    = -1          ! Final call at the end of the simulation.
   INTEGER, PARAMETER :: DISCON_STATUS_INITIALISING  =  0          ! First call at time zero.
   INTEGER, PARAMETER :: DISCON_STATUS_DISCRETE_STEP =  1          ! Simulation discrete timestep.
   
   INTEGER, PARAMETER :: DISCON_PITCH_CONTROL_COLLECTIVE = 0       ! Pitch is controlled collectively - use GetCollectivePitchAngle and SetDemandedCollectivePitchAngle.
   INTEGER, PARAMETER :: DISCON_PITCH_CONTROL_INDIVIDUAL = 1       ! Pitch is controlled on each blade individually - use GetPitchAngle and SetDemandedPitchAngle.

   INTEGER, PARAMETER :: ErrID_Info  = 0                           ! Information flag
   INTEGER, PARAMETER :: ErrID_Fatal = 1                           ! Fatal error flag
   INTEGER, PARAMETER :: ErrID_Severe = 2                          ! Severe error flag
   INTEGER, PARAMETER :: ErrID_None  = -1                          ! No error flag

   INTEGER, PARAMETER :: ErrMsgLen = 100                           ! Maximum length of error message
   
   ! Derived type for the external controller (DLL) interface
   TYPE, PUBLIC :: controller_t
      REAL(RPC) , DIMENSION(:), ALLOCATABLE  :: avrSWAP                      ! SWAP array used to pass data to and from the external controller (see Bladed documentation for more info)
      INTEGER(IPC)                           :: avcOUTNAME_LEN = CLEN        ! Length of avcOUTNAME (output channel name)
      REAL(KIND=RP)                          :: azimuth                      ! Rotor azimuth angle (blade 1) (rad)
      REAL(KIND=RPC), DIMENSION(3)           :: blade_pitch                  ! Blade pitch angles (rad)
      REAL(RPC), DIMENSION(3)                :: blade_pitch_com              ! Demanded blade pitch angles (rad)
      INTEGER(IPC)                           :: brake_state                  ! Shaft brake status flag
      REAL(RPC)                              :: controller_dt                ! Time step used for communication between solver and external controller (s)
      REAL(RPC)                              :: controller_t_ini             ! Initial time for the external controller operation (s)
      CHARACTER(len=CLEN)                    :: controller_dll               ! Path to external controller dynamic library file
      INTEGER                                :: controller_id                ! Controller ID (for farm control)
      CHARACTER(CLEN)                        :: controller_inputs            ! Path to the input file used by the external controller
      INTEGER                                :: controller_iter              ! Controller iteration number
      CHARACTER(len=LINE_LENGTH)             :: controller_restart           ! Controller restart file name
      REAL(KIND=RPC)                         :: drivetrain_inertia           ! Moment of intertia of the drivetrain around rotation axis (kg m^2)
      INTEGER                                :: err_status                   ! Error status of the controller
      CHARACTER(len=CLEN)                    :: err_msg                      ! Error message of the controller
      REAL(KIND=RPC)                         :: gear_ratio                   ! LSS to HSS gear ratio
      REAL(KIND=RPC)                         :: gen_eff                      ! Generator efficiency (-)
      REAL(KIND=RPC)                         :: gen_elec_power               ! Generator electrical power (W)
      INTEGER(IPC)                           :: gen_state                    ! Generator contactor state flag
      REAL(KIND=RPC)                         :: gen_trq                      ! Generator torque (Nm)
      REAL(RPC)                              :: gen_trq_com                  ! Demanded generator torque (Nm)
      REAL(KIND=RPC)                         :: HSS_speed                    ! High-speed shaft (HSS) speed (rad/s)
      LOGICAL                                :: initialized                  ! Boolean flag to indicate if the controller has been initialized
      REAL(KIND=RPC)                         :: M_aero                       ! Rotor aerodynamic torque (Nm)
      REAL(KIND=RPC)                         :: M_brake                      ! High-speed shaft (HSS) brake torque (Nm)
      REAL(RPC)                              :: M_brake_com                  ! Demanded brake torque on high-speed shaft (HSS) (Nm)
      REAL(KIND=RPC), DIMENSION(3)           :: My_blade_root                ! Blade root out-of-plane bending moments (Nm)
      INTEGER(IPC)                           :: n_trq_lookup = 0             ! Number of points in the torque-speed lookup table (0 if unused)
      INTEGER(IPC)                           :: pitch_state                  ! Pitch control flag
      CHARACTER(CLEN)                        :: root_name                    ! Root path for writing output files
      REAL(KIND=RPC)                         :: rotor_power                  ! Rotor power (W)
      REAL(KIND=RP)                          :: rotor_speed                  ! Rotor or Lowe-speed shaft (LSS) speed (rad/s)
      INTEGER(IPC)                           :: sim_status                   ! Simulation status flag
      REAL(KIND=RPC)                         :: wind_dir                     ! Wind direction at the wind probe (rad)
      REAL(KIND=RP)                          :: wind_speed                   ! Wind speed at the wind probe (m/s)
      REAL(KIND=RP), DIMENSION(3)            :: wind_probe_coords            ! Coordinates of the wind speed/direction probe for controller (m)
      REAL(RPC)                              :: yaw_rate_com                 ! Demanded yaw rate (rad/s)      
      PROCEDURE(DISCON), POINTER, NOPASS     :: DISCON_ptr                   ! Pointer to the DISCON procedure in the external controller
   END TYPE controller_t

   ! Definition of interface for DISCON subroutine in external controller
   INTERFACE
      SUBROUTINE DISCON (avrSWAP, aviFAIL, accINFILE, avcOUTNAME, avcMSG)  BIND(C, NAME='DISCON')
         USE, INTRINSIC :: ISO_C_Binding   
         REAL(C_FLOAT),          INTENT(INOUT) :: avrSWAP   (*)  ! Array for data exchange between solver and controller
         INTEGER(C_INT),         INTENT(INOUT) :: aviFAIL        ! Status flag set in the controller and returned to solver (0: successful call, >0: successfull call with warning avcMSG, <0: unsuccessful call or stop with error avcMSG)
         CHARACTER(KIND=C_CHAR), INTENT(IN)    :: accINFILE (*)  ! Path to the controller parameter file, commonly named DISCON.IN
         CHARACTER(KIND=C_CHAR), INTENT(IN)    :: avcOUTNAME(*)  ! Path to simulation root for controller oputput logging
         CHARACTER(KIND=C_CHAR), INTENT(INOUT) :: avcMSG    (*)  ! Message returned from controller to solver
      END SUBROUTINE DISCON
   END INTERFACE

   ! Definition of interfaces for dynamic loading of libraries (from the C standard library dlfcn.h)
   INTERFACE
      FUNCTION dlsym(handle, symbol) bind(C, NAME="dlsym")
         USE, INTRINSIC :: ISO_C_Binding
         IMPLICIT NONE 
         TYPE(C_PTR)                                 :: dlsym
         TYPE(C_PTR),                   VALUE        :: handle
         CHARACTER(KIND=C_CHAR, LEN=1), DIMENSION(*) :: symbol
      END FUNCTION dlsym
   END INTERFACE

   INTERFACE
      FUNCTION dlmopen(lmid, filename, mode) bind(c, name="dlmopen")
         USE, INTRINSIC :: ISO_C_Binding
         IMPLICIT NONE
         INTEGER(C_LONG),   VALUE      :: lmid
         CHARACTER(C_CHAR), INTENT(IN) :: filename(*)
         INTEGER(C_INT),    VALUE      :: mode
         TYPE(C_PTR)                   :: dlmopen
      END FUNCTION dlmopen
   END INTERFACE

   INTERFACE
      FUNCTION dlclose(handle) bind(C, NAME="dlclose")
         USE, INTRINSIC :: ISO_C_Binding
         IMPLICIT NONE
         INTEGER(C_INT)     :: dlclose
         TYPE(C_PTR), VALUE :: handle
      END FUNCTION dlclose
   END INTERFACE
   
CONTAINS

   ! This subroutine allocates memory for the avrSWAP array
   SUBROUTINE AllocAvrSWAP (avrSWAP, avrSWAPdim, Descr, ErrStat, ErrMsg)
   
      REAL(RPC),    ALLOCATABLE  :: avrSWAP(:)        !  Array to be allocated
      INTEGER,      INTENT(IN)   :: avrSWAPdim        !  Size of the array
      CHARACTER(*), INTENT(IN)   :: Descr             !  Array description
      INTEGER,      INTENT(OUT)  :: ErrStat           !  Error status
      CHARACTER(*), INTENT(OUT)  :: ErrMsg            !  Error message corresponding to ErrStat
   
      ALLOCATE ( avrSWAP(avrSWAPdim) , STAT=ErrStat )
   
      IF ( ErrStat /= 0 ) THEN
         ErrStat = ErrID_Fatal
         IF ( ALLOCATED(avrSWAP) ) THEN 
            ErrMsg = 'Error allocating memory for the '//TRIM( Descr )//' array; array was already allocated.'
         ELSE
            ErrMsg = 'Error allocating memory for the '//TRIM( Descr )//' array.'
         END IF
      ELSE
         ErrStat = ErrID_None
         ErrMsg  = ''
         avrSWAP = 0.0_RPC
      END IF
   
   END SUBROUTINE AllocAvrSWAP
   
   ! This subroutine calls de external DISCON controller
   SUBROUTINE CallController(controller_data)
      
      TYPE(controller_t), INTENT(INOUT)  :: controller_data                                   ! Data TYPE containing controller data, including avrSWAP, accINFILE, and avcOUTNAME arrays
   
      INTEGER(C_INT)                     :: aviFAIL                                           ! Status flag set in the controller and returned to solver (0: successful call, >0: successfull call with warning avcMSG, <0: unsuccessful call or stop with error avcMSG)
      CHARACTER(KIND=C_CHAR)             :: accINFILE(NINT(controller_data%avrSWAP(50)))      ! Path to the controller parameter file, commonly named DISCON.IN
      CHARACTER(KIND=C_CHAR)             :: avcOUTNAME(NINT(controller_data%avrSWAP(51)))     ! Path to simulation root for controller oputput logging
      CHARACTER(KIND=C_CHAR)             :: avcMSG(NINT(controller_data%avrSWAP(49)))         ! Message returned from controller to solver
      
      ! initialize aviFAIL
      aviFAIL = 0
   
      ! From Fortran CHARACTER to C CHARACTER (trailing C_NULL_CHAR)
      avcOUTNAME  = TRANSFER(TRIM(controller_data%root_name)//C_NULL_CHAR, avcOUTNAME, NINT(controller_data%avrSWAP(51)))
      accINFILE   = TRANSFER(TRIM(controller_data%controller_inputs)//C_NULL_CHAR, accINFILE, NINT(controller_data%avrSWAP(50)))
      avcMSG      = TRANSFER(C_NULL_CHAR, avcMSG, NINT(controller_data%avrSWAP(49))) 
       
      ! Call the DISCON procedure through pointer
      CALL controller_data%DISCON_ptr(controller_data%avrSWAP, aviFAIL, accINFILE, avcOUTNAME, avcMSG)
   
      ! Handle errors from external controller (if any)
      IF (aviFAIL /= 0) THEN
   
         controller_data%err_msg = TRANSFER(avcMSG,controller_data%err_msg) ! From C to Fortran CHARACTER
          
         IF (aviFAIL > 0) THEN
            controller_data%err_status = ErrID_Info
         ELSE
            controller_data%err_status = ErrID_Fatal
         END IF   
      
      ELSE
         controller_data%err_status = ErrID_None
         controller_data%err_msg = ''
      END IF
   
      ! Update simulation status flag
      IF (controller_data%sim_status == DISCON_STATUS_FINALISING) THEN
        controller_data%sim_status = DISCON_STATUS_INITIALISING
      ELSE
         controller_data%sim_status = DISCON_STATUS_DISCRETE_STEP
      END IF
   
   END SUBROUTINE CallController
   
   ! This subroutine initialises all the variables required for the external controller interface
   SUBROUTINE ControllerInterfaceInit(controller_data)
      
      TYPE(controller_t), INTENT(INOUT)  :: controller_data   ! Data TYPE containing controller data, including the avrSWAP
   
      INTEGER                            :: ErrStat2          ! Error status flag
      CHARACTER(20)                      :: ErrMsg2           ! Error message
   
      controller_data%err_status = ErrID_None
      controller_data%err_msg= ''

      ! Load the external controller dynamic library and assign the DISCON procedure pointer
      call LoadController(controller_data)
   
      IF (controller_data%controller_dt < EPSILON(controller_data%controller_dt)) THEN
         ! TO DO: Print error to file (controller time step must be greater than zero)
      END IF
   
      ! Set status flag and initialize avrSWAP
      controller_data%sim_status = DISCON_STATUS_INITIALISING
      
      ! Allocate memory for the avrSWAP array
      CALL AllocAvrSWAP(controller_data%avrSwap,   R+(2*controller_data%n_trq_lookup)-1 + MaxLoggingChannels, 'avrSwap', ErrStat2, ErrMsg2)
   
      ! Initialize dll data stored in OtherState
      controller_data%initialized = .FALSE.
   
   END SUBROUTINE ControllerInterfaceInit
   
   ! This subroutine calls the DLL for the final time (if it was previously called), and frees the corresponding dynamic library.
   SUBROUTINE ControllerInterfaceEnd(controller_data)
   
      TYPE(controller_t), INTENT(INOUT)  :: controller_data  ! Data TYPE containing the avrSWAP and sim_status variables among others
   
      ! Call DLL at final time (only if it was previously called)
      IF (allocated(controller_data%avrSWAP)) THEN
         IF (controller_data%sim_status /= DISCON_STATUS_INITIALISING) THEN
            controller_data%sim_status = DISCON_STATUS_FINALISING
            controller_data%avrSWAP(1) = controller_data%sim_status
            CALL CallController(controller_data)
         END IF
      END IF
   
   END SUBROUTINE ControllerInterfaceEnd
   
   ! This subroutine fills the avrSWAP array with its inputs, as described in Appendices A and B of the Bladed User Manual of Bladed (v4.3 and above).
   SUBROUTINE Fill_avrSWAP(t, controller_data)
      
      REAL,               INTENT(IN)     :: t 
      TYPE(controller_t), INTENT(INOUT)  :: controller_data 
   
      ! Assign the values to the avrSWAP array
      controller_data%avrSWAP( 1) = controller_data%sim_status                     ! Channel  1: Simulation controller call status flag (0: first call, 1: every other time step, -1: final call)
      controller_data%avrSWAP( 2) = REAL(t, RPC)                                   ! Channel  2: Simulation time (s)
      controller_data%avrSWAP( 3) = controller_data%controller_dt                  ! Channel  3: Controller call interval (s) 
      controller_data%avrSWAP( 4) = controller_data%blade_pitch(1)                 ! Channel  4: Pitch angle of blade 1 (rad) 
      
      controller_data%avrSWAP(14) = controller_data%rotor_power                    ! Channel 14: Rotor power (W) 
      controller_data%avrSWAP(15) = controller_data%gen_elec_power                 ! Channel 15: Generator electrical power output (W)
      
      controller_data%avrSWAP(20) = controller_data%HSS_speed                      ! Channel 20: Generator speed (rad/s) 
      controller_data%avrSWAP(21) = REAL(controller_data%rotor_speed, RPC)         ! Channel 21: Rotor speed (rad/s) 
      
      controller_data%avrSWAP(23) = controller_data%gen_trq                        ! Channel 23: Generator torque (Nm)
      controller_data%avrSWAP(24) = 0.0                                            ! Channel 24: Yaw error (rad) 
      
      controller_data%avrSWAP(27) = REAL(controller_data%wind_speed, RPC)          ! Channel 27: Hub wind speed (m/s) 
      
      controller_data%avrSWAP(30) = controller_data%My_blade_root(1)               ! Channel 30: Root out-of-plane bending moment of blade 1 (Nm) 
      controller_data%avrSWAP(31) = controller_data%My_blade_root(2)               ! Channel 31:                 ||                  blade 2 (Nm) 
      controller_data%avrSWAP(32) = controller_data%My_blade_root(3)               ! Channel 32:                 ||                  blade 3 (Nm) 
      controller_data%avrSWAP(33) = controller_data%blade_pitch(2)                 ! Channel 33: Pitch angle of blade 2 (rad) 
      controller_data%avrSWAP(34) = controller_data%blade_pitch(3)                 ! Channel 34:         ||     blade 3 (rad)
      
      controller_data%avrSWAP(49) = LEN(controller_data%err_msg)+1                 ! Channel 49: Maximum no. of characters in avcMSG of DISCON subroutine (C string)
      controller_data%avrSWAP(50) = LEN_TRIM(controller_data%controller_inputs)+1  ! Channel 50: Length of accINFILE in DISCON subroutine (C string)
      controller_data%avrSWAP(51) = controller_data%avcOUTNAME_LEN                 ! Channel 51: Length of avcOUTNAME in DISCON subroutine (C string)
      
      controller_data%avrSWAP(53) = 0.0                                            ! Channel 53: Fore-aft acceleration of tower top (m/s^2) 
      controller_data%avrSWAP(54) = 0.0                                            ! Channel 54: Side-to-side acceleration of tower top (m/s^2) 
      
      controller_data%avrSWAP(60) = REAL(controller_data%azimuth, RPC)             ! Channel 60: Rotor azimuth angle (rad) 
      controller_data%avrSWAP(61) = 3.0                                            ! Channel 61: No. of blades (-)
               
      controller_data%avrSWAP(98) = 0                                              ! Channel 98: set to 0. This channel is used as output in Retrieve_avrSWAP()
      
      controller_data%avrSWAP(117) = 0                                             ! Channel 117: Controller state [always set to 0]
   
      ! Channels 120-129: User-defined variables 1-10
      controller_data%avrSWAP(120) = controller_data%controller_id                 ! Channel 120: User-defined variable 1 - Controller ID (for farm control) (-)
      controller_data%avrSWAP(121) = controller_data%controller_iter               ! Channel 121: User-defined variable 2 - Controller iteration number (-)
      
      controller_data%avrSWAP(129) = size(controller_data%avrSWAP)                 ! Channel 129: Size of the avrSWAP array
   
      ! Channels 130-142, L1 and on are outputs
   
   END SUBROUTINE Fill_avrSWAP
   
   ! This routine retrieves the DLL return VALUEs from the avrSWAP array, as described in Appendices A and B of the Bladed User Manual (v4.3 and above).
   SUBROUTINE Retrieve_avrSWAP(controller_data)
   
      TYPE(controller_t), INTENT(INOUT)  :: controller_data           ! Data for the Bladed DLL
   
      INTEGER                  :: K                                   ! Counter variable
      CHARACTER(*), PARAMETER  :: RoutineName = 'Retrieve_avrSWAP'    ! Name of this routine (for error messages)
   
   
      ! Initialize err_status and err_msg
      controller_data%err_status = ErrID_None
      controller_data%err_msg  = ''
   
      ! Retrieve values from the avrSWAP array
      controller_data%gen_state  = NINT(controller_data%avrSWAP(35))                  ! Channel 35: Generator contactor (-)
      controller_data%brake_state = NINT(controller_data%avrSWAP(36))                 ! Channel 36: Shaft brake status (0 or 1) (-)
   
      IF ( controller_data%pitch_state == DISCON_PITCH_CONTROL_INDIVIDUAL )  THEN
         DO K = 1,3 
            controller_data%blade_pitch_com(K) = controller_data%avrSWAP( 41 + K )    ! Channels 42-44: Demanded individual pitch position of blade K (rad)
         ENDDO
      ELSE
         controller_data%blade_pitch_com(:)   = controller_data%avrSWAP(45)           ! Channel 45: Demanded pitch angle (Collective pitch) (rad)
      ENDIF

      ! Channel 46: Demanded collective pitch rate (rad/s)   (not used)
   
      controller_data%gen_trq_com  = controller_data%avrSWAP(47)                      ! Channel 47: Demanded generator torque (Nm)
      controller_data%yaw_rate_com = controller_data%avrSWAP(48)                      ! Channel 48: Demanded nacelle yaw rate (rad/s)
   
      ! Channel 55: Pitch override flag                      (not used)
      ! Channel 56: Torque override flag                     (not used)
      
      ! Channel 65: Number of variables returned for logging (not used)
      
      ! Channel 102: Yaw control flag                        (not used)
   
      IF (controller_data%brake_state == 16) THEN
         controller_data%M_brake_com = controller_data%avrSWAP(107)                   ! Channel 107: Brake torque demand on high-speed shaft (HSS), used only when avrSWAP(36) is 16 (Nm)
      END IF
   
   END SUBROUTINE Retrieve_avrSWAP
   
   ! This subroutine loads the external controller (dynamic library) and assigns the DISCON procedure pointer to its corresponding DISCON_ptr 
   SUBROUTINE LoadController (controller_data)
      USE, INTRINSIC :: ISO_C_Binding
      TYPE(controller_t), INTENT(INOUT)  :: controller_data   ! Data TYPE containing the controller DLL path and the DISCON procedure pointer
   
      TYPE(C_PTR)                        :: handle            ! Handle to the dynamically loaded library
      TYPE(C_PTR)                        :: sym_ptr           ! Pointer to the DISCON symbol in the library
      TYPE(C_FUNPTR)                     :: func_ptr          ! Function pointer used to transfer DISCON address
      INTEGER(C_INT)                     :: dlerror           ! Error code from dynamic loading
      INTEGER(C_INT),     PARAMETER      :: RTLD_NOW = 2      ! Define RTLD_NOW with platform-specific value 
      INTEGER(C_LONG),    PARAMETER      :: LM_ID_NEWLM = -1  ! Define LM_ID_NEWLM with platform-specific value
      CHARACTER(CLEN)                    :: libname_c         ! Path to controller DLL with C string terminator

      ! Load the shared library
      libname_c = trim(controller_data%controller_dll) // c_null_char
      handle = dlmopen(LM_ID_NEWLM, libname_c, RTLD_NOW)
      IF (.not. c_associated(handle)) THEN
          PRINT *, "Error: Unable to load " // trim(libname_c) // " library."
          STOP
      ELSE
          IF (MPI_Process % isRoot) THEN
            WRITE(STD_OUT,'(30X,A)') "-> INFO: Successfully loaded " // trim(libname_c) // " external controller."
          END IF
      END IF
   
      ! Get the function pointer
      sym_ptr = dlsym(handle, "DISCON" // c_null_char)
   
      func_ptr = transfer(sym_ptr, func_ptr)
      
      CALL C_F_PROCPOINTER(func_ptr, controller_data%DISCON_ptr)
   
      ! Check if the procedure pointer is associated
      IF (ASSOCIATED(controller_data%DISCON_ptr)) THEN
         IF (MPI_Process % isRoot) THEN
            WRITE(STD_OUT,'(30X,A)') "-> INFO: Function pointer to DISCON routine is associated."
         END IF
      ELSE
          PRINT *, "Error: Function pointer to DISCON is not associated."
      END IF
   
   END SUBROUTINE LoadController

END MODULE ControllerInterface
