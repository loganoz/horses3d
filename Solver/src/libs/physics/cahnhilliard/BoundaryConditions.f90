!
!////////////////////////////////////////////////////////////////////////
!
!      BoundaryConditions.f90
!      Created: 2011-07-21 10:30:14 -0400 
!      By: David Kopriva
!      
!      This file has the BCs for Navier-Stokes computations.
!      The two entry points are:
!
!         SUBROUTINE ExternalStateForBoundaryName( x, y, t, nHat, Q, boundaryName )
!
!            Returns the state Q for the boundary with name boundaryName
!
!         ExternalGradientForBoundaryName( x, y, t, nHat, gradU, boundaryName )
!
!            Returns the gradient gradU for the boundary named boundaryName
!
!      The two entry points call boundary routines depending on the type 
!      of boundary conditions. Currently supported BCs are:
!
!      (1) Free-Stream Inflow-Outflow (constant external state)
!          Dirichlet: Free stream value
!          Neumann  : All gradients = 0
!
!          Set up in:
!
!             SUBROUTINE UniformFlowState( x, t, Q )
!             SUBROUTINE UniformFlowNeumann( x, t, nHat, gradU )
!
!      (2) Free-Slip Wall (symmetry boundary, normal velocity and normal gradients vanish)
!
!             Dirichlet: \vec q^ext \cdot \hat n = - \vec q^int \cdot \hat n
!             Neumann  : \vec \nabla q^ext \cdot \hat n = - \vec \nabla q^int \cdot \hat n
!
!          Set up in:
!
!             SUBROUTINE FreeSlipWallState ( x, t, nHat, Q )
!             SUBROUTINE FreeSlipNeumann ( x, t, nHat, grad )
!
!      (3) Adiabataic Wall (velocity = 0, normal temperature gradient = 0)
!
!             Dirichlet: \vec q^ext = - \vec q^int
!             Neumann  : \vec \nabla q(3)^ext \cdot \hat n = - \vec \nabla q(3)^int \cdot \hat n
!
!          Set up in:
!
!             SUBROUTINE NoSlipAdiabaticWallState ( x, t, Q )
!             SUBROUTINE NoSlipAdiabaticNeumann( x, t, nHat, gradU )
!
!      (4) Isothermal Wall (velocity = 0, temperature set to wall value, gradient from interior)
!
!          Set up in:
!
!             SUBROUTINE NoSlipIsothermalWallState( x, t, Q )
!             SUBROUTINE NoSlipIsothermalWallNeumann( x, t, gradU )
!
!
!////////////////////////////////////////////////////////////////////////
!
      MODULE BoundaryConditionFunctions
         USE SMConstants
         USE Physics
         USE SharedBCModule
         use PhysicsStorage

         private

         public C_BC, MU_BC
         public implementedBCNames

         public NoFluxState, NoFluxNeumann, WallAngleBC
         public UserDefinedState, UserDefinedNeumann

         CHARACTER(LEN=BC_STRING_LENGTH), DIMENSION(5) :: implementedBCNames = &
               ["no-flux             ", &
                "periodic+           ", &
                "periodic-           ", &
                "noslipadiabaticwall ", &
                "user-defined        "   ]
               
         integer, parameter :: C_BC = 1
         integer, parameter :: MU_BC = 2
         integer, parameter :: NO_FLUX_INDEX        = 1
         integer, parameter :: PERIODIC_PLUS_INDEX  = 2
         integer, parameter :: PERIODIC_MINUS_INDEX = 3
         INTEGER, PARAMETER :: USER_DEFINED_INDEX   = 4
!
!     ========         
      CONTAINS
!     ========
! 
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE NoFluxState( x, t, nHat, Q )
!
!     ----------------------------------------------------
!        No flux state: boundary = interior: does nothing
!     ----------------------------------------------------
!
      IMPLICIT NONE 
      REAL(KIND=RP), INTENT(IN)    :: x(3), t
      REAL(KIND=RP), INTENT(IN)    :: nHat(3)
      REAL(KIND=RP), INTENT(INOUT) :: Q(N_EQN)

      END SUBROUTINE NoFluxState
!
!     /////////////////////////////////////////////////////////////////
!
      SUBROUTINE NoFluxNeumann(x, t, nHat, U_x, U_y, U_z )
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), INTENT(IN)    :: x(3), t
      REAL(KIND=RP), INTENT(IN)    :: nHat(3)
      REAL(KIND=RP), INTENT(INOUT) :: U_x(N_GRAD_EQN), U_y(N_GRAD_EQN), U_z(N_GRAD_EQN)

      U_x = 0.0_RP
      U_y = 0.0_RP
      U_z = 0.0_RP

      END SUBROUTINE NoFluxNeumann

      subroutine WallAngleBC(x, t, nHat, U_x, U_y, U_z)
      implicit none
      real(kind=RP), intent(in)  :: x(3), t
      real(kind=RP), intent(in)  :: nHat(3)
      REAL(KIND=RP), INTENT(INOUT) :: U_x(N_GRAD_EQN), U_y(N_GRAD_EQN), U_z(N_GRAD_EQN)

      U_x = thermodynamics % thetaw * nHat(1)
      U_y = thermodynamics % thetaw * nHat(2)
      U_z = thermodynamics % thetaw * nHat(3)

      end subroutine WallAngleBC 
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE UniformFlowState( x, t, Q )
      USE SMConstants
      USE PhysicsStorage
      IMPLICIT NONE
      
      REAL(KIND=RP) :: x(3), t
      REAL(KIND=RP) :: Q(N_EQN)
      
      END SUBROUTINE UniformFlowState
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE UniformFlowNeumann( x, t, nHat, U_x, U_y, U_z  )
!
!     ----------------------------
!     Set normal gradients to zero
!     ----------------------------
!
      USE SMConstants
      USE PhysicsStorage
      IMPLICIT NONE
      
      REAL(KIND=RP) :: x(3), t
      REAL(KIND=RP) :: nHat(3)
      REAL(KIND=RP), INTENT(INOUT) :: U_x(N_GRAD_EQN), U_y(N_GRAD_EQN), U_z(N_GRAD_EQN)
      
      END SUBROUTINE UniformFlowNeumann
!
!////////////////////////////////////////////////////////////////////////
!
!! Perdurbs a state vector for a particular point located in the unit cube
!! using the random number generator and a gaussian function
! 
   SUBROUTINE GaussianPerturbUnitCube(x, Q)
      IMPLICIT NONE
      
      REAL(KIND=RP) :: x(3)
      REAL(KIND=RP) :: Q(N_EQN)
      
      REAL(KIND=RP) :: GaussFac, RandNum
      INTEGER       :: i
      
      IF(ANY(x>1._RP) .OR. ANY(x<0._RP)) THEN
         print*, 'WARNING: Points in geometry outside of unit cube.'
      END IF
      
      GaussFac = exp((-(x(1)-0.5_RP)**2-(x(2)-0.5_RP)**2-(x(3)-0.5_RP)**2)*20)
      
      DO i=1, N_EQN
         CALL RANDOM_NUMBER(RandNum)
         Q(i) = Q(i) + Q(i) * GaussFac * (RandNum-0.5_RP) * 0.5_RP
      END DO
      
   END SUBROUTINE GaussianPerturbUnitCube
!
!////////////////////////////////////////////////////////////////////////
!
!! Perdurbs a state vector for a particular point located in the unit cube
!! using the random number generator and a gaussian function
! 
   SUBROUTINE GaussianPerturbUnitSquare(x, Q)
      IMPLICIT NONE
      
      REAL(KIND=RP) :: x(3)
      REAL(KIND=RP) :: Q(N_EQN)
      
      REAL(KIND=RP) :: GaussFac, RandNum
      INTEGER       :: i
      
      IF(ANY(x(1:2)>1._RP) .OR. ANY(x(1:2)<0._RP)) THEN
         print*, 'WARNING: Points in geometry outside of unit square.'
      END IF
      
      GaussFac = exp((-(x(1)-0.5_RP)**2-(x(2)-0.5_RP)**2)*20)
      
      DO i=1, N_EQN
         CALL RANDOM_NUMBER(RandNum)
         Q(i) = Q(i) + Q(i) * GaussFac * (RandNum-0.5_RP) * 0.5_RP
      END DO
      
   END SUBROUTINE GaussianPerturbUnitSquare

      subroutine UserDefinedState(x, t, nHat, Q)
         implicit none
         real(kind=RP)  :: x(NDIM)
         real(kind=RP)  :: t
         real(kind=RP)  :: nHat(NDIM)
         real(kind=RP)  :: Q(N_EQN)
         interface
            subroutine UserDefinedState1(x, t, nHat, Q, thermodynamics_, dimensionless_, refValues_)
               use SMConstants
               use PhysicsStorage
               implicit none
               real(kind=RP)  :: x(NDIM)
               real(kind=RP)  :: t
               real(kind=RP)  :: nHat(NDIM)
               real(kind=RP)  :: Q(N_EQN)
               type(Thermodynamics_t), intent(in)  :: thermodynamics_
               type(Dimensionless_t),  intent(in)  :: dimensionless_
               type(RefValues_t),      intent(in)  :: refValues_
            end subroutine UserDefinedState1
         end interface

         call UserDefinedState1(x, t, nHat, Q, thermodynamics, dimensionless, refValues)

      end subroutine UserDefinedState
!
!////////////////////////////////////////////////////////////////////////
!
!     ===========
      END MODULE BoundaryConditionFunctions
!     ===========
!
!=====================================================================================================
!=====================================================================================================
!
!
      SUBROUTINE externalStateForBoundaryName( x, t, nHat, Q, boundaryType )
!
!     ----------------------------------------------
!     Set the boundary conditions for the mesh by
!     setting the external state for each boundary.
!     ----------------------------------------------
!
      use SMConstants
      USE BoundaryConditionFunctions
      use PhysicsStorage
      
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP)   , INTENT(IN)    :: x(3), t, nHat(3)
      REAL(KIND=RP)   , INTENT(INOUT) :: Q(N_EQN)
      CHARACTER(LEN=*), INTENT(IN)    :: boundaryType
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP)   :: pExt
      LOGICAL         :: success

      CALL NoFluxState( x, t, nHat, Q )

      END SUBROUTINE externalStateForBoundaryName
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ExternalConcentrationGradientForBoundaryName( x, t, nHat, GradU, boundaryType )
!
!     ------------------------------------------------
!     Set the boundary conditions for the mesh by
!     setting the external gradients on each boundary.
!     ------------------------------------------------
!
      use SMConstants
      USE BoundaryConditionFunctions
      use PhysicsStorage
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP)   , INTENT(IN)    :: x(3), t, nHat(3)
      REAL(KIND=RP)   , INTENT(INOUT) :: GradU(3,N_GRAD_EQN)
      CHARACTER(LEN=*), INTENT(IN)    :: boundaryType
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP) :: U_x(N_GRAD_EQN), U_y(N_GRAD_EQN), U_z(N_GRAD_EQN)

      U_x(:) = GradU(1,:)
      U_y(:) = GradU(2,:)
      U_z(:) = GradU(3,:)

      IF ( boundarytype == "no-flux" )                   THEN
         CALL NoFluxNeumann( x, t, nHat, U_x, U_y, U_z )
      ELSEIF ( boundaryType == "noslipadiabaticwall" ) then
         call WallAngleBC(x, t, nHat, U_x, U_y, U_z)
      END IF

      GradU(1,:) = U_x(:)
      GradU(2,:) = U_y(:)
      GradU(3,:) = U_z(:)

      END SUBROUTINE ExternalConcentrationGradientForBoundaryName

      SUBROUTINE ExternalChemicalPotentialGradientForBoundaryName( x, t, nHat, GradU, boundaryType )
!
!     ------------------------------------------------
!     Set the boundary conditions for the mesh by
!     setting the external gradients on each boundary.
!     ------------------------------------------------
!
      use SMConstants
      USE BoundaryConditionFunctions
      use PhysicsStorage
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP)   , INTENT(IN)    :: x(3), t, nHat(3)
      REAL(KIND=RP)   , INTENT(INOUT) :: GradU(3,N_GRAD_EQN)
      CHARACTER(LEN=*), INTENT(IN)    :: boundaryType
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP) :: U_x(N_GRAD_EQN), U_y(N_GRAD_EQN), U_z(N_GRAD_EQN)

      U_x(:) = GradU(1,:)
      U_y(:) = GradU(2,:)
      U_z(:) = GradU(3,:)

      CALL NoFluxNeumann( x, t, nHat, U_x, U_y, U_z )

      GradU(1,:) = U_x(:)
      GradU(2,:) = U_y(:)
      GradU(3,:) = U_z(:)

      END SUBROUTINE ExternalChemicalPotentialGradientForBoundaryName
