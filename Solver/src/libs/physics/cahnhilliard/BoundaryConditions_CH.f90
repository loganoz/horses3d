!
!//////////////////////////////////////////////////////
!
!   @File:    BoundaryConditions_CH.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Thu Apr 19 17:24:29 2018
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
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

         public C_BC, MU_BC
         public implementedBCNames

         public NoFluxState, NoFluxNeumann, WallAngleBC
         public UserDefinedState, UserDefinedNeumann

         CHARACTER(LEN=BC_STRING_LENGTH), DIMENSION(5) :: implementedCHBCNames = &
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
      subroutine UserDefinedConcentrationState(x, t, nHat, Q)
         implicit none
         real(kind=RP)  :: x(NDIM)
         real(kind=RP)  :: t
         real(kind=RP)  :: nHat(NDIM)
         real(kind=RP)  :: Q(NCOMP)
         interface
            subroutine UserDefinedConcentrationState1(x, t, nHat, Q, multiphase_)
               use SMConstants
               use PhysicsStorage_CH
               use FluidData_CH
               implicit none
               real(kind=RP)      :: x(NDIM)
               real(kind=RP)      :: t
               real(kind=RP)      :: nHat(NDIM)
               real(kind=RP)      :: Q(NCOMP)
               type(Multiphase_t) :: multiphase_
            end subroutine UserDefinedConcentrationState1
         end interface

         call UserDefinedConcentrationState1(x, t, nHat, Q, multiphase)

      end subroutine UserDefinedConcentrationState
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
      SUBROUTINE externalCHStateForBoundaryName( x, t, nHat, Q, boundaryType )
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
      REAL(KIND=RP)   , INTENT(IN)    :: x(3), t, nHat(3)
      REAL(KIND=RP)   , INTENT(INOUT) :: Q(NCOMP)
      CHARACTER(LEN=*), INTENT(IN)    :: boundaryType
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
      SUBROUTINE ExternalConcentrationGradientForBoundaryName( x, t, nHat, GradU, boundaryType )
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
      REAL(KIND=RP)   , INTENT(IN)    :: x(3), t, nHat(3)
      REAL(KIND=RP)   , INTENT(INOUT) :: GradU(3,NCOMP)
      CHARACTER(LEN=*), INTENT(IN)    :: boundaryType
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP) :: U_x(NCOMP), U_y(NCOMP), U_z(NCOMP)

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
      USE BoundaryConditionFunctions_CH
      use PhysicsStorage_CH
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP)   , INTENT(IN)    :: x(3), t, nHat(3)
      REAL(KIND=RP)   , INTENT(INOUT) :: GradU(3,NCOMP)
      CHARACTER(LEN=*), INTENT(IN)    :: boundaryType
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP) :: U_x(NCOMP), U_y(NCOMP), U_z(NCOMP)

      U_x(:) = GradU(1,:)
      U_y(:) = GradU(2,:)
      U_z(:) = GradU(3,:)

      CALL NoFluxNeumann( x, t, nHat, U_x, U_y, U_z )

      GradU(1,:) = U_x(:)
      GradU(2,:) = U_y(:)
      GradU(3,:) = U_z(:)

      END SUBROUTINE ExternalChemicalPotentialGradientForBoundaryName
