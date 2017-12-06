!
!////////////////////////////////////////////////////////////////////////
!
!      HORSES3DMain.f90
!      Created: May 21, 2015 at 12:56 PM 
!      By: David Kopriva  
!
!////////////////////////////////////////////////////////////////////////
!
!
!////////////////////////////////////////////////////////////////////////
!
      PROGRAM HORSES3DMain
      
      USE SMConstants
      use FTValueDictionaryClass
      USE PhysicsStorage
      USE SharedBCModule
      USE zoneClass
      USE DGSEMClass
      USE BoundaryConditionFunctions
      USE TimeIntegratorClass
      USE mainKeywordsModule
      USE Headers
      USE pAdaptationClass
      use StopwatchClass
      use MPI_Process_Info
#ifdef _HAS_MPI_
      use mpi
#endif
      
      IMPLICIT NONE
interface
         SUBROUTINE UserDefinedStartup
            IMPLICIT NONE  
         END SUBROUTINE UserDefinedStartup
         SUBROUTINE UserDefinedFinalSetup(mesh , thermodynamics_, &
                                                 dimensionless_, &
                                                     refValues_ )
            use PhysicsStorage
            use HexMeshClass
            IMPLICIT NONE
            CLASS(HexMesh)             :: mesh
            type(Thermodynamics_t),    intent(in)  :: thermodynamics_
            type(Dimensionless_t),     intent(in)  :: dimensionless_
            type(RefValues_t),         intent(in)  :: refValues_
         END SUBROUTINE UserDefinedFinalSetup
         SUBROUTINE UserDefinedFinalize(mesh, time, iter, maxResidual, thermodynamics_, &
                                                    dimensionless_, &
                                                        refValues_, &
                                                          monitors, &
                                                       elapsedTime, &
                                                           CPUTime   )
            use PhysicsStorage
            use HexMeshClass
            use MonitorsClass
            IMPLICIT NONE
            CLASS(HexMesh)                        :: mesh
            REAL(KIND=RP)                         :: time
            integer                               :: iter
            real(kind=RP)                         :: maxResidual
            type(Thermodynamics_t),    intent(in) :: thermodynamics_
            type(Dimensionless_t),     intent(in) :: dimensionless_
            type(RefValues_t),         intent(in) :: refValues_
            type(Monitor_t),          intent(in) :: monitors
            real(kind=RP),             intent(in) :: elapsedTime
            real(kind=RP),             intent(in) :: CPUTime
         END SUBROUTINE UserDefinedFinalize
      SUBROUTINE UserDefinedTermination
         IMPLICIT NONE  
      END SUBROUTINE UserDefinedTermination
      character(len=LINE_LENGTH) function getFileName(inputLine)
         use SMConstants
         implicit none
         character(len=*)    :: inputLine
      end function getFileName
end interface

      TYPE( FTValueDictionary)            :: controlVariables
      TYPE( DGSem )                       :: sem
      TYPE( TimeIntegrator_t )            :: timeIntegrator
      
      LOGICAL                             :: success, saveGradients
      integer                             :: initial_iteration
      INTEGER                             :: ierr
      INTEGER, EXTERNAL                   :: UnusedUnit
      real(kind=RP)                       :: initial_time
      EXTERNAL                            :: externalStateForBoundaryName
      EXTERNAL                            :: ExternalGradientForBoundaryName
      character(len=LINE_LENGTH)          :: solutionFileName
      
      ! For pAdaptation
      INTEGER, ALLOCATABLE                :: Nx(:), Ny(:), Nz(:)
      INTEGER                             :: Nmax
!
!     ---------------
!     Initializations
!     ---------------
!
      call MPI_Process % Init
      CALL Main_Header("HORSES3D High-Order (DG) Spectral Element Solver",__DATE__,__TIME__)

      CALL controlVariables % initWithSize(16)
      CALL UserDefinedStartup
      CALL ConstructSharedBCModule
      
      CALL ReadInputFile( controlVariables )
      CALL CheckInputIntegrity(controlVariables, success)
      IF(.NOT. success)   ERROR STOP "Control file reading error"
      
!
!     ----------------
!     Set up the DGSEM
!     ----------------
!      
      CALL ConstructPhysicsStorage( controlVariables, success )
      IF(.NOT. success)   ERROR STOP "Physics parameters input error"
      
      ! Initialize manufactured solutions if necessary
      sem % ManufacturedSol = controlVariables % containsKey("manufactured solution")
      
      IF (sem % ManufacturedSol) THEN
         CALL InitializeManufacturedSol(controlVariables % StringValueForKey("manufactured solution",LINE_LENGTH))
      END IF
      
      call GetMeshPolynomialOrders(controlVariables,Nx,Ny,Nz,Nmax)
      call InitializeNodalStorage(Nmax)
      
      call sem % construct (  controlVariables  = controlVariables,                                         &
                                 externalState     = externalStateForBoundaryName,                             &
                                 externalGradients = ExternalGradientForBoundaryName,                          &
                                 Nx_ = Nx,     Ny_ = Ny,     Nz_ = Nz,                                                 &
                                 success           = success)
                           
      IF(.NOT. success)   ERROR STOP "Mesh reading error"
      CALL checkBCIntegrity(sem % mesh, success)
      IF(.NOT. success)   ERROR STOP "Boundary condition specification error"
      CALL UserDefinedFinalSetup(sem % mesh, thermodynamics, dimensionless, refValues)
!
!     -------------------------
!     Set the initial condition
!     -------------------------
!
      call sem % SetInitialCondition(controlVariables, initial_iteration, initial_time)
!
!     -----------------------------
!     Construct the time integrator
!     -----------------------------
!
      CALL timeIntegrator % construct (controlVariables, initial_iteration, initial_time)
!
!     -----------------
!     Integrate in time
!     -----------------
!
      CALL timeIntegrator % integrate(sem, controlVariables, sem % monitors)
!
!     --------------------------
!     Show simulation statistics
!     --------------------------
!
      call DisplaySimulationStatistics(sem % numberOftimeSteps, sem % mesh)
!
!     -----------------------------------------------------
!     Let the user perform actions on the computed solution
!     -----------------------------------------------------
!
      CALL UserDefinedFinalize(sem % mesh, timeIntegrator % time, sem % numberOfTimeSteps, &
                              sem % maxResidual, thermodynamics, dimensionless, refValues, &
                              sem % monitors, Stopwatch % ElapsedTime("Solver"), &
                              Stopwatch % CPUTime("Solver"))
#ifdef _HAS_MPI_
      if ( MPI_Process % doMPIAction ) then
         call mpi_barrier(MPI_COMM_WORLD, ierr)
      end if
#endif
!
!     -------------------------------------
!     Save the results to the solution file
!     -------------------------------------
!
      IF(controlVariables % stringValueForKey(solutionFileNameKey,LINE_LENGTH) /= "none")     THEN 
         solutionFileName = trim(getFileName(controlVariables % stringValueForKey(solutionFileNameKey,LINE_LENGTH))) // ".hsol"
         saveGradients    = controlVariables % logicalValueForKey(saveGradientsToSolutionKey)
         CALL sem % mesh % SaveSolution(sem % numberOfTimeSteps, timeIntegrator % time, solutionFileName, saveGradients)
      END IF


      CALL timeIntegrator % destruct()
      CALL sem % destruct()
      CALL destructSharedBCModule
      
      CALL UserDefinedTermination

      call MPI_Process % Close
      
      END PROGRAM HORSES3DMain
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE CheckBCIntegrity(mesh, success)
!
         USE HexMeshClass
         use FTValueDictionaryClass
         USE SharedBCModule
         USE BoundaryConditionFunctions, ONLY:implementedBCNames
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(HexMesh) :: mesh
         LOGICAL       :: success
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER                              :: i, j
         INTEGER                              :: faceID, eId
         CHARACTER(LEN=BC_STRING_LENGTH)      :: bcName, namedBC
         CHARACTER(LEN=BC_STRING_LENGTH)      :: bcType
         TYPE(FTMutableObjectArray), POINTER :: bcObjects
         CLASS(FTValue)             , POINTER :: v
         CLASS(FTObject), POINTER             :: obj
         
         success = .TRUE.
!
!        ----------------------------------------------------------
!        Check to make sure that the boundaries defined in the mesh
!        have an associated name in the control file.
!        ----------------------------------------------------------
         
         DO eID = 1, SIZE( mesh % elements )
            DO faceID = 1, 6
               namedBC = mesh % elements(eId) % boundaryName(faceID)
               IF( namedBC == emptyBCName ) CYCLE
               
               bcName = bcTypeDictionary % stringValueForKey(key             = namedBC,         &
                                                             requestedLength = BC_STRING_LENGTH)
               IF ( LEN_TRIM(bcName) == 0 )     THEN
                  PRINT *, "Control file does not define a boundary condition for boundary name = ", &
                            mesh % elements(eId) % boundaryName(faceID)
                  success = .FALSE.
                  return 
               END IF 
            END DO   
         END DO
!
!        --------------------------------------------------------------------------
!        Check that the boundary conditions to be applied are implemented
!        in the code. Keep those updated in the boundary condition functions module
!        --------------------------------------------------------------------------
!
         bcObjects => bcTypeDictionary % allObjects()
         DO j = 1, bcObjects % COUNT()
            obj => bcObjects % objectAtIndex(j)
            CALL castToValue(obj,v)
            bcType = v % stringValue(requestedLength = BC_STRING_LENGTH)
            DO i = 1, SIZE(implementedBCNames)
               IF ( bcType == implementedBCNames(i) )     THEN
                  success = .TRUE. 
                  EXIT 
               ELSE 
                  success = .FALSE. 
               END IF 
            END DO
            
            IF ( .NOT. success )     THEN
               PRINT *, "Boundary condition ", TRIM(bcType)," not implemented in this code"
               CALL release(bcObjects)
               RETURN 
            END IF  
            
         END DO
         
         CALL release(bcObjects)
         
      END SUBROUTINE checkBCIntegrity
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE CheckInputIntegrity( controlVariables, success )  
         use SMConstants
         USE FTValueDictionaryClass
         USE mainKeywordsModule
         use FTValueClass
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(FTValueDictionary) :: controlVariables
         LOGICAL                 :: success
!
!        ---------------
!        Local variables
!        ---------------
!
         CLASS(FTObject), POINTER :: obj
         INTEGER                  :: i
         character(len=LINE_LENGTH)    :: inviscidDiscretization, discretizationNodes
         interface
            subroutine toLower(str) 
               character(*), intent(in out) :: str
            end subroutine toLower
         end interface
         
         success = .TRUE.
!
!        Control variables with default value
!        ------------------------------------
         obj => controlVariables % objectForKey(saveGradientsToSolutionKey)
         if ( .not. associated(obj) ) then
            call controlVariables % addValueForKey(".false.",saveGradientsToSolutionKey)
         end if

         obj => controlVariables % objectForKey(discretizationNodesKey)
         if ( .not. associated(obj) ) then
            call controlVariables % addValueForKey("Gauss",discretizationNodesKey)
         end if

         obj => controlVariables % objectForKey(inviscidDiscretizationKey)
         if ( .not. associated(obj) ) then
            call controlVariables % addValueForKey("Standard",inviscidDiscretizationKey)
         end if

         obj => controlVariables % objectForKey(splitFormKey)
         if ( .not. associated(obj) ) then
            call controlVariables % addValueForKey("Ducros",splitFormKey)
         end if
!
!        Check for inconsistencies in the input variables
!        ------------------------------------------------
         inviscidDiscretization = trim(controlVariables % stringValueForKey(inviscidDiscretizationKey, LINE_LENGTH))
         discretizationNodes = trim(controlVariables % stringValueForKey(discretizationNodesKey, LINE_LENGTH))

         call toLower(inviscidDiscretization)
         call toLower(discretizationNodes)

         if ( (trim(inviscidDiscretization) .eq. "split-form") .and. (trim(discretizationNodes) .eq. "gauss") ) then
            write(STD_OUT,'(A)') "*** WARNING:    Only Gauss-Lobatto nodes are available for Split-Form discretizations"
            write(STD_OUT,'(A)') "*** WARNING:    Automatically switched to Gauss-Lobatto points"
            call controlVariables % removeObjectForKey(discretizationNodesKey)
            call controlVariables % addValueForKey("Gauss-Lobatto",discretizationNodesKey)
         end if
!
!        Check the controlVariables created
!        ----------------------------------        
         DO i = 1, SIZE(mainKeywords)
            obj => controlVariables % objectForKey(mainKeywords(i))
            IF ( .NOT. ASSOCIATED(obj) )     THEN
               PRINT *, "Input file is missing entry for keyword: ",mainKeywords(i)
               success = .FALSE. 
            END IF  
         END DO  
         
      END SUBROUTINE checkInputIntegrity

      subroutine DisplaySimulationStatistics(iter,mesh)
         use SMConstants
         use HexMeshClass
         use StopwatchClass
         use Headers
         use MPI_Process_Info
         implicit none
         integer,    intent(in)      :: iter
         type(HexMesh),   intent(in) :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                    :: eID
         integer                    :: NDOF
         real(kind=RP)              :: Naverage
         real(kind=RP)              :: t_elaps, t_cpu
   
         if ( .not. MPI_Process % isRoot ) return

         write(STD_OUT,'(/)')
         call Section_Header("Simulation statistics")
         write(STD_OUT,'(/)')
!
!        Get mesh-related quantities
!        ---------------------------
         NDOF = 0
         Naverage = 0
   
         do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )
            NDOF = NDOF + (e % Nxyz(1) + 1)*(e % Nxyz(2) + 1)*(e % Nxyz(3) + 1)      
            Naverage = Naverage + e % Nxyz(1) + e % Nxyz(2) + e % Nxyz(3)
            end associate
         end do

         Naverage = Naverage / (3.0_RP * mesh % no_of_elements)
!
!        Show preprocessing time
!        -----------------------
         t_elaps = Stopwatch % Elapsedtime("Preprocessing")
         t_cpu   = Stopwatch % CPUTime("Preprocessing")

         call Subsection_Header("Preprocessing")

         write(STD_OUT,'(30X,A,I0,A,F5.2,A,I0,A)')      "->   ", mesh % no_of_elements, &
                                                      " elements with polynomial order ",Naverage," (NDOF = ",NDOF,")."
         write(STD_OUT,'(30X,A,A30,ES10.3,A,ES10.3,A)') "->", "Preprocessing time: ",t_elaps," seconds (total CPU time: ",t_cpu,")."

!
!        Show simulation time
!        --------------------
         write(STD_OUT,'(/)')
         call Subsection_Header("Solver")
         if ( iter .le. 0 ) return

         t_elaps = Stopwatch % ElapsedTime("Solver")
         t_cpu   = Stopwatch % CPUTime("Solver")

         write(STD_OUT,'(30X,A,A30,ES10.3,A)') "->", "Simulation elapsed time: ",t_elaps," seconds."
         write(STD_OUT,'(30X,A,A30,ES10.3,A,ES10.3,A)') "->", "Simulation CPU time: ",t_cpu," seconds (ratio is ",t_cpu/t_elaps ,")."
         write(STD_OUT,'(30X,A,A30,ES10.3,A)') "->", "Solver efficiency: " , t_elaps/(NDOF * iter)*1.0e6_RP, " seconds/(1 Million DOF·iter)."

      end subroutine DisplaySimulationStatistics
