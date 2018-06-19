!
!//////////////////////////////////////////////////////
!
!   @File:    BoundaryConditions_iNS.f90
!   @Author:  Juan Manzanero (j.manzanero1992@gmail.com)
!   @Created: Tue Jun 19 17:39:25 2018
!   @Last revision date: Wed Jun 20 18:14:41 2018
!   @Last revision author: Juan Manzanero (j.manzanero1992@gmail.com)
!   @Last revision commit: 9c8ed8b6306ad0912cb55b510aa73d1610bb1cb5
!
!//////////////////////////////////////////////////////
!
MODULE BoundaryConditionFunctions_iNS
   USE SMConstants
   USE Physics_iNS
   USE SharedBCModule
   use PhysicsStorage_iNS
   use FluidData_iNS

   private

   public implementediNSBCNames

   CHARACTER(LEN=BC_STRING_LENGTH), DIMENSION(2) :: implementediNSBCNames = &
         ["periodic-           ", &
          "periodic+           "]

   enum, bind(C)
      enumerator :: PERIODIC_PLUS_INDEX=1, PERIODIC_MINUS_INDEX
   end enum
!
!  ========         
   contains
!  ========
! 
!
!////////////////////////////////////////////////////////////////////////
!
!     ===========
      END MODULE BoundaryConditionFunctions_iNS
!     ===========
!
!=====================================================================================================
!=====================================================================================================
!
!
      SUBROUTINE externalStateForBoundaryName_iNS( nEqn, x, t, nHat, Q, boundaryType, boundaryName )
!
!     ----------------------------------------------
!     Set the boundary conditions for the mesh by
!     setting the external state for each boundary.
!     ----------------------------------------------
!
      use SMConstants
      USE BoundaryConditionFunctions_iNS
      use PhysicsStorage_iNS
      
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      integer         , intent(in)    :: nEqn
      REAL(KIND=RP)   , INTENT(IN)    :: x(3), t, nHat(3)
      REAL(KIND=RP)   , INTENT(INOUT) :: Q(nEqn)
      CHARACTER(LEN=BC_STRING_LENGTH), INTENT(IN)    :: boundaryType
      CHARACTER(LEN=BC_STRING_LENGTH), INTENT(IN)    :: boundaryName
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP)   :: pExt
      LOGICAL         :: success

      END SUBROUTINE externalStateForBoundaryName_iNS
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ExternalGradientForBoundaryName_iNS( nGradEqn, x, t, nHat, GradU, boundaryType, boundaryName )
!
!     ------------------------------------------------
!     Set the boundary conditions for the mesh by
!     setting the external gradients on each boundary.
!     ------------------------------------------------
!
      use SMConstants
      USE BoundaryConditionFunctions_iNS
      use PhysicsStorage_iNS
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      integer,          intent(in)    :: nGradEqn
      REAL(KIND=RP)   , INTENT(IN)    :: x(3), t, nHat(3)
      REAL(KIND=RP)   , INTENT(INOUT) :: GradU(3,nGradEqn)
      CHARACTER(LEN=BC_STRING_LENGTH), INTENT(IN)    :: boundaryType
      CHARACTER(LEN=BC_STRING_LENGTH), INTENT(IN)    :: boundaryName

      END SUBROUTINE ExternalGradientForBoundaryName_iNS
