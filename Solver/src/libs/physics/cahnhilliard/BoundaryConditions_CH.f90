!
!//////////////////////////////////////////////////////
!
!   @File:    BoundaryConditions_CH.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Thu Apr 19 17:24:29 2018
!   @Last revision date: Wed May  9 15:26:11 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 030a84de7dedac0cada2e2d9ba22dfd63aa09eb8
!
!//////////////////////////////////////////////////////
!
      MODULE BoundaryConditionFunctions_CH
         USE SMConstants
         USE Physics_CH
         USE SharedBCModule
         use PhysicsStorage_CH
         use FluidData_CH

         private

         public implementedCHBCNames

         public NoFluxState, NoFluxNeumann, WallAngleBC
         public UserDefinedState, UserDefinedNeumann

         CHARACTER(LEN=BC_STRING_LENGTH), DIMENSION(7) :: implementedCHBCNames = &
               ["no-flux             ", &
                "periodic+           ", &
                "periodic-           ", &
                "noslipadiabaticwall ", &
                "inflow              ", &
                "outflowspecifyp     ", &
                "user-defined        "   ]

         enum, bind(C)
            enumerator :: NO_FLUX_INDEX=1, PERIODIC_PLUS_INDEX, PERIODIC_MINUS_INDEX
            enumerator :: INFLOW_INDEX, OUTFLOWSPECIFYP_INDEX
            enumerator :: USER_DEFINED_INDEX
         end enum
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
      REAL(KIND=RP), INTENT(INOUT) :: Q(NCOMP)

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
      REAL(KIND=RP), INTENT(INOUT) :: U_x(NCOMP), U_y(NCOMP), U_z(NCOMP)

      U_x = 0.0_RP
      U_y = 0.0_RP
      U_z = 0.0_RP

      END SUBROUTINE NoFluxNeumann

      subroutine WallAngleBC(x, t, nHat, U_x, U_y, U_z)
      implicit none
      real(kind=RP), intent(in)  :: x(3), t
      real(kind=RP), intent(in)  :: nHat(3)
      REAL(KIND=RP), INTENT(INOUT) :: U_x(NCOMP), U_y(NCOMP), U_z(NCOMP)

      U_x = multiphase % thetaw * nHat(1)
      U_y = multiphase % thetaw * nHat(2)
      U_z = multiphase % thetaw * nHat(3)

      end subroutine WallAngleBC 
!
!////////////////////////////////////////////////////////////////////////
!
!     ===========
      END MODULE BoundaryConditionFunctions_CH
!     ===========
!
!=====================================================================================================
!=====================================================================================================
!
!
      SUBROUTINE externalCHStateForBoundaryName( nEqn, x, t, nHat, Q, boundaryType, boundaryName )
!
!     ----------------------------------------------
!     Set the boundary conditions for the mesh by
!     setting the external state for each boundary.
!     ----------------------------------------------
!
      use SMConstants
      USE BoundaryConditionFunctions_CH
      use PhysicsStorage_CH
      
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      integer         , intent(in)    :: nEqn
      REAL(KIND=RP)   , INTENT(IN)    :: x(3), t, nHat(3)
      REAL(KIND=RP)   , INTENT(INOUT) :: Q(nEqn)
      CHARACTER(LEN=*), INTENT(IN)    :: boundaryType
      CHARACTER(LEN=*), INTENT(IN)    :: boundaryName
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP)   :: pExt
      LOGICAL         :: success

      CALL NoFluxState( x, t, nHat, Q )

      END SUBROUTINE externalCHStateForBoundaryName
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ExternalConcentrationGradientForBoundaryName( nGradEqn, x, t, nHat, GradU, boundaryType, boundaryName )
!
!     ------------------------------------------------
!     Set the boundary conditions for the mesh by
!     setting the external gradients on each boundary.
!     ------------------------------------------------
!
      use SMConstants
      USE BoundaryConditionFunctions_CH
      use PhysicsStorage_CH
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      integer,          intent(in)    :: nGradEqn
      REAL(KIND=RP)   , INTENT(IN)    :: x(3), t, nHat(3)
      REAL(KIND=RP)   , INTENT(INOUT) :: GradU(3,nGradEqn)
      CHARACTER(LEN=*), INTENT(IN)    :: boundaryType
      CHARACTER(LEN=*), INTENT(IN)    :: boundaryName
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP) :: U_x(nGradEqn), U_y(nGradEqn), U_z(nGradEqn)

      U_x(:) = GradU(1,:)
      U_y(:) = GradU(2,:)
      U_z(:) = GradU(3,:)

      IF ( boundarytype == "no-flux" )                   THEN
         CALL NoFluxNeumann( x, t, nHat, U_x, U_y, U_z )
      ELSEIF ( boundaryType == "noslipadiabaticwall" ) then
         call WallAngleBC(x, t, nHat, U_x, U_y, U_z)
      elseif ( boundaryType == "inflow" ) then
         CALL NoFluxNeumann( x, t, nHat, U_x, U_y, U_z )
      elseif ( boundaryType == "outflowspecifyp") then
         CALL NoFluxNeumann( x, t, nHat, U_x, U_y, U_z )
      END IF

      GradU(1,:) = U_x(:)
      GradU(2,:) = U_y(:)
      GradU(3,:) = U_z(:)

      END SUBROUTINE ExternalConcentrationGradientForBoundaryName

      SUBROUTINE ExternalChemicalPotentialGradientForBoundaryName( nGradEqn, x, t, nHat, GradU, boundaryType, boundaryName )
!
!     ------------------------------------------------
!     Set the boundary conditions for the mesh by
!     setting the external gradients on each boundary.
!     ------------------------------------------------
!
      use SMConstants
      USE BoundaryConditionFunctions_CH
      use PhysicsStorage_CH
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      integer,          intent(in)    :: nGradEqn
      REAL(KIND=RP)   , INTENT(IN)    :: x(3), t, nHat(3)
      REAL(KIND=RP)   , INTENT(INOUT) :: GradU(3,nGradEqn)
      CHARACTER(LEN=*), INTENT(IN)    :: boundaryType
      CHARACTER(LEN=*), INTENT(IN)    :: boundaryName
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP) :: U_x(nGradEqn), U_y(nGradEqn), U_z(nGradEqn)

      U_x(:) = GradU(1,:)
      U_y(:) = GradU(2,:)
      U_z(:) = GradU(3,:)

      CALL NoFluxNeumann( x, t, nHat, U_x, U_y, U_z )

      GradU(1,:) = U_x(:)
      GradU(2,:) = U_y(:)
      GradU(3,:) = U_z(:)

      END SUBROUTINE ExternalChemicalPotentialGradientForBoundaryName
