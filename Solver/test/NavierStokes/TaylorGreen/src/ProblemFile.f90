!
!////////////////////////////////////////////////////////////////////////
!
!      ProblemFile.f90
!      Created: June 26, 2015 at 8:47 AM 
!      By: David Kopriva  
!
!      The Problem File contains user defined procedures
!      that are used to "personalize" i.e. define a specific
!      problem to be solved. These procedures include initial conditions,
!      exact solutions (e.g. for tests), etc. and allow modifications 
!      without having to modify the main code.
!
!      The procedures, *even if empty* that must be defined are
!
!      UserDefinedStartup
!      UserDefinedInitialCondition(sem)
!      UserDefinedPeriodicOperation(sem)
!      UserDefinedFinalize(sem)
!      UserDefinedTermination
!
!      *** This problem file sets up a Taylor-Green vortex test case (Re=1600 M=0.08) *** 
!
!//////////////////////////////////////////////////////////////////////// 
!
      MODULE UserDefinedDataStorage
         USE SMConstants
         IMPLICIT NONE 
         REAL(KIND=RP) :: rad0, f, h 
      END MODULE UserDefinedDataStorage
!
!//////////////////////////////////////////////////////////////////////// 
! 
      MODULE UserDefinedFunctions

!
!     ========      
      CONTAINS
!     ========
!
         SUBROUTINE UserDefinedStartup  
!
!        --------------------------------
!        Called before any other routines
!        --------------------------------
!
            IMPLICIT NONE  
         END SUBROUTINE UserDefinedStartup
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalSetup(sem, controlVariables)
!
!           ----------------------------------------------------------------------
!           Called after the mesh is read in to allow mesh related initializations
!           or memory allocations.
!           ----------------------------------------------------------------------
!
            USE DGSEMClass
            IMPLICIT NONE
            CLASS(DGSem)             :: sem
            class(FTValueDictionary) :: controlVariables
         END SUBROUTINE UserDefinedFinalSetup
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedInitialCondition(sem, controlVariables)
!
!           ------------------------------------------------
!           Called to set the initial condition for the flow
!           Simulation of the Taylor-Green vortex using
!           high-order flux reconstruction schemes.
!           J. Bull, A. Jameson
!           AIAA Journal 2015, vol 53-9, pages 1-30
!           ------------------------------------------------
!
            USE SMConstants
            USE DGSEMClass
            USE PhysicsStorage
            USE BoundaryConditionFunctions
            IMPLICIT NONE
            
            TYPE(DGSem)      :: sem
            class(FTValueDictionary) :: controlVariables            
            EXTERNAL         :: initialStateSubroutine
            
            REAL(KIND=RP) :: x(3)        
            INTEGER       :: i, j, k, eID
            REAL(KIND=RP) :: rho , u , v , w , p
            REAL(KIND=RP) :: L, u_0, rho_0, p_0
            integer       :: Nx, Ny, Nz
            
            L     = 1.0_RP
            u_0   = 1.0_RP
            rho_0 = 1.0_RP 
            p_0   = 100.0_RP

            DO eID = 1, SIZE(sem % mesh % elements)
               Nx = sem % mesh % elements(eID) % Nxyz(1)
               Ny = sem % mesh % elements(eID) % Nxyz(2)
               Nz = sem % mesh % elements(eID) % Nxyz(3)

               DO k = 0, Nz
                  DO j = 0, Ny
                     DO i = 0, Nx 

                         x = sem % mesh % elements(eID) % geom % x(:,i,j,k)
                       
                         rho = rho_0
                         u   =  u_0 * sin(x(1)/L) * cos(x(2)/L) * cos(x(3)/L) 
                         v   = -u_0 * cos(x(1)/L) * sin(x(2)/L) * cos(x(3)/L)
                         w   =  0.0_RP
                         p   =   p_0 + rho_0 / 16.0_RP * (                          &
                               cos(2.0_RP*x(1)/L)*cos(2.0_RP*x(3)/L) +                  &
                               2.0_RP*cos(2.0_RP*x(2)/L) + 2.0_RP*cos(2.0_RP*x(1)/L) +  &
                               cos(2.0_RP*x(2)/L)*cos(2.0_RP*x(3)/L)                    &
                               )

                         sem % mesh % elements(eID) % Q(i,j,k,1) = rho
                         sem % mesh % elements(eID) % Q(i,j,k,2) = rho*u
                         sem % mesh % elements(eID) % Q(i,j,k,3) = rho*v
                         sem % mesh % elements(eID) % Q(i,j,k,4) = rho*w
                         sem % mesh % elements(eID) % Q(i,j,k,5) = p / (gamma - 1.0_RP) + 0.5_RP * rho * (u*u + v*v + w*w)

                     END DO
                  END DO
               END DO 
               
            END DO 
            
         END SUBROUTINE UserDefinedInitialCondition
!
!//////////////////////////////////////////////////////////////////////// 
! 
!
         SUBROUTINE UserDefinedPeriodicOperation(sem, time)
!
!           ----------------------------------------------------------
!           Called at the output interval to allow periodic operations
!           to be performed
!           ----------------------------------------------------------
!
            USE DGSEMClass
            IMPLICIT NONE
            CLASS(DGSem)        :: sem
            REAL(KIND=RP)       :: time

            LOGICAL, SAVE       :: created = .FALSE.
            REAL(KIND=RP)       :: k_energy = 0.0_RP
            INTEGER             :: fUnit
            INTEGER, SAVE       :: restartfile = 0
            INTEGER             :: restartUnit
            CHARACTER(LEN=100)  :: fileName

            fUnit = 101

            IF (.NOT. created) THEN
                OPEN(fUnit, file = "./kinetic_energy.dat", action = "WRITE" , status = "UNKNOWN" )
                created = .TRUE.
            ELSE
                OPEN(fUnit, file = "./kinetic_energy.dat", action = "WRITE" , status = "OLD" , access = "APPEND")
            END IF
!           -------------------------------
!           Compute the flow kinetic energy
!           -------------------------------
            k_energy = computeKineticEnergy(sem)
            WRITE(fUnit , '(E24.16,E24.16)') time , k_energy

            CLOSE(fUnit)

         END SUBROUTINE UserDefinedPeriodicOperation
!
!//////////////////////////////////////////////////////////////////////// 
! 
         FUNCTION computeKineticEnergy(sem) result(val)
            USE DGSEMClass
            IMPLICIT NONE
            CLASS(DGSem)        :: sem
            REAL(kind=RP)       :: val
            
            REAL(kind=RP)       :: rho, u , v , w
            INTEGER             :: eID , i , j , k
            integer             :: Nx, Ny, Nz

!           Initialization
            val = 0.0_RP

            DO eID = 1 , SIZE(sem % mesh % elements)
               Nx = sem % mesh % elements(eID) % Nxyz(1)
               Ny = sem % mesh % elements(eID) % Nxyz(2)
               Nz = sem % mesh % elements(eID) % Nxyz(3)

                DO k = 0 , Nz
                    DO j = 0 , Ny
                        DO i = 0 , Nx

                            rho = sem % mesh % elements(eID) % Q(i,j,k,1)
                            u = sem % mesh % elements(eID) % Q(i,j,k,2) / rho
                            v = sem % mesh % elements(eID) % Q(i,j,k,3) / rho
                            w = sem % mesh % elements(eID) % Q(i,j,k,4) / rho

                            ASSOCIATE(wx => sem % spA(Nx,Ny,Nz) % wx , &
                                      wy => sem % spA(Nx,Ny,Nz) % wy , &
                                      wz => sem % spA(Nx,Ny,Nz) % wz   )
                              val = val + 0.5_RP*rho*(u*u+v*v+w*w)*wx(i)*wy(j)*wz(k)*      &
                                    sem % mesh % elements (eID) % geom % jacobian (i,j,k)
                            END ASSOCIATE
                             
                        END DO
                    END DO
                END DO
            END DO

            val = val / (2*acos(-1.d0))**3

         END FUNCTION computeKineticEnergy      
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalize(sem, time)
            USE FTAssertions
!
!           --------------------------------------------------------
!           Called after the solution computed to allow, for example
!           error tests to be performed
!           --------------------------------------------------------
!
            USE DGSEMClass
            IMPLICIT NONE
!
!           ---------
!           Arguments
!           ---------
!
            CLASS(DGSem)  :: sem
            REAL(KIND=RP) :: time
!
!           ---------------
!           Local variables
!           ---------------
!
            CHARACTER(LEN=29)                  :: testName           = "Taylor-Green vortex"
            REAL(KIND=RP)                      :: maxError
            REAL(KIND=RP), ALLOCATABLE         :: QExpected(:,:,:,:)
            INTEGER                            :: eID
            INTEGER                            :: i, j, k, N
            TYPE(FTAssertionsManager), POINTER :: sharedManager
            LOGICAL                            :: success

            WRITE(6,*) "This test case has no expected solution yet."

            
         END SUBROUTINE UserDefinedFinalize
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE UserDefinedTermination
!
!        -----------------------------------------------
!        Called at the the end of the main driver after 
!        everything else is done.
!        -----------------------------------------------
!
         IMPLICIT NONE  
      END SUBROUTINE UserDefinedTermination
      
      END MODULE UserDefinedFunctions
!
!=====================================================================================================
!=====================================================================================================
!!
!!
!      SUBROUTINE externalStateForBoundaryName( x, t, nHat, Q, boundaryType )
!!
!!     ----------------------------------------------
!!     Set the boundary conditions for the mesh by
!!     setting the external state for each boundary.
!!     ----------------------------------------------
!!
!      USE BoundaryConditionFunctions
!      USE UserDefinedDataStorage
!      USE MeshTypes
!      
!      IMPLICIT NONE
!!
!!     ---------
!!     Arguments
!!     ---------
!!
!      REAL(KIND=RP)   , INTENT(IN)    :: x(3), t, nHat(3)
!      REAL(KIND=RP)   , INTENT(INOUT) :: Q(N_EQN)
!      CHARACTER(LEN=*), INTENT(IN)    :: boundaryType
!!
!!     ---------------
!!     Local variables
!!     ---------------
!!
!      REAL(KIND=RP)   :: pExt
!      LOGICAL         :: success
!
!      IF ( boundarytype == "freeslipwall" )             THEN
!         CALL FreeSlipWallState( x, t, nHat, Q )
!      ELSE IF ( boundaryType == "noslipadiabaticwall" ) THEN 
!         CALL  NoSlipAdiabaticWallState( x, t, Q)
!      ELSE IF ( boundarytype == "noslipisothermalwall") THEN 
!         CALL NoSlipIsothermalWallState( x, t, Q )
!      ELSE IF ( boundaryType == "outflowspecifyp" )     THEN 
!         pExt =  ExternalPressure()
!         CALL ExternalPressureState ( x, t, nHat, Q, pExt )
!      ELSE 
!         CALL UniformFlowState( x, t, Q ) 
!      END IF
!
!      END SUBROUTINE externalStateForBoundaryName
!!
!////////////////////////////////////////////////////////////////////////
!!
!      SUBROUTINE ExternalGradientForBoundaryName( x, t, nHat, GradU, boundaryType )
!!
!!     ------------------------------------------------
!!     Set the boundary conditions for the mesh by
!!     setting the external gradients on each boundary.
!!     ------------------------------------------------
!!
!      USE BoundaryConditionFunctions
!      USE MeshTypes
!      IMPLICIT NONE
!!
!!     ---------
!!     Arguments
!!     ---------
!!
!      REAL(KIND=RP)   , INTENT(IN)    :: x(3), t, nHat(3)
!      REAL(KIND=RP)   , INTENT(INOUT) :: GradU(3,N_GRAD_EQN)
!      CHARACTER(LEN=*), INTENT(IN)    :: boundaryType
!!
!!     ---------------
!!     Local variables
!!     ---------------
!!
!      REAL(KIND=RP) :: U_x(N_GRAD_EQN), U_y(N_GRAD_EQN), U_z(N_GRAD_EQN)
!
!      U_x(:) = GradU(1,:)
!      U_y(:) = GradU(2,:)
!      U_z(:) = GradU(3,:)
!
!      IF ( boundarytype == "freeslipwall" )                   THEN
!         CALL FreeSlipNeumann( x, t, nHat, U_x, U_y, U_z )
!      ELSE IF ( boundaryType == "noslipadiabaticwall" )       THEN 
!         CALL  NoSlipAdiabaticWallNeumann( x, t, nHat, U_x, U_y, U_z )
!      ELSE IF ( boundarytype == "noslipisothermalwall")       THEN 
!         CALL NoSlipIsothermalWallNeumann( x, t, nHat, U_x, U_y, U_z )
!      ELSE
!         CALL UniformFlowNeumann( x, t, nHat, U_x, U_y, U_z )
!      END IF
!
!      GradU(1,:) = U_x(:)
!      GradU(2,:) = U_y(:)
!      GradU(3,:) = U_z(:)
!
!      END SUBROUTINE ExternalGradientForBoundaryName
!
