
!////////////////////////////////////////////////////////////////////////
!
!      DGSEMClass.f95
!      Created: 2008-07-12 13:38:26 -0400 
!      By: David Kopriva
!
!      Basic class for the discontinuous Galerkin spectral element
!      solution of conservation laws.
!
!      Algorithms:
!         Algorithm 136: DGSEM Class
!         Algorithm 129: ConstructDGSem (Constructor)
!         Algorithm 138: ComputeTimeDerivative (TimeDerivative)
!         Algorithm 137: ComputeRiemannFluxes (EdgeFluxes)
!         Algorithm  35: InterpolateToFineMesh (2DCoarseToFineInterpolation)
!
!      Modified for 3D 6/11/15, 11:32 AM by DAK
!
!////////////////////////////////////////////////////////////////////////
!
      Module DGSEMClass
      
      USE NodalStorageClass
      USE HexMeshClass
      USE PhysicsStorage
      
      IMPLICIT NONE 
!
!-------------------------------------------------------------------
! Basic class for the Spectral element solution of conservation laws
!-------------------------------------------------------------------
!
      TYPE DGSem
         TYPE(NodalStorage) :: spA
         TYPE(HexMesh)      :: mesh
!
!        ========         
         CONTAINS
!        ========         
!
         PROCEDURE :: construct => ConstructDGSem
         PROCEDURE :: destruct  => DestructDGSem     
            
      END TYPE DGSem
      
      CONTAINS 
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ConstructDGSem( self, polynomialOrder, meshFileName, success )
!
!     --------------------------
!     Constructor for the class.
!     --------------------------
!
      IMPLICIT NONE
      CLASS(DGSem)     :: self
      INTEGER          :: polynomialOrder
      CHARACTER(LEN=*) :: meshFileName
      LOGICAL          :: success
      
      INTEGER :: k
      
      CALL ConstructNodalStorage( self % spA, polynomialOrder )
      CALL self % mesh % constructFromFile( meshfileName, self % spA, success )
      
      IF(.NOT. success) RETURN 
!
!     ------------------------
!     Allocate and zero memory
!     ------------------------
!
      DO k = 1, SIZE(self % mesh % elements) 
         CALL allocateElementStorage( self % mesh % elements(k), polynomialOrder, &
                                      N_EQN, N_GRAD_EQN, flowIsNavierStokes )
      END DO
      
      END SUBROUTINE ConstructDGSem
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DestructDGSem( self )
      IMPLICIT NONE 
      CLASS(DGSem) :: self
      
      CALL DestructNodalStorage( self % spA )
      CALL DestructMesh( self % mesh )
      
      END SUBROUTINE DestructDGSem
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE SaveSolutionForRestart( self, fUnit ) 
         IMPLICIT NONE
         CLASS(DGSem)      :: self
         INTEGER          :: fUnit
         INTEGER          :: k

         DO k = 1, SIZE(self % mesh % elements) 
            WRITE(fUnit) self % mesh % elements(k) % Q
         END DO

      END SUBROUTINE SaveSolutionForRestart
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE LoadSolutionForRestart( self, fUnit ) 
         IMPLICIT NONE
         CLASS(DGSem)      :: self
         INTEGER          :: fUnit
         INTEGER          :: k

         DO k = 1, SIZE(self % mesh % elements) 
            READ(fUnit) self % mesh % elements(k) % Q
         END DO

      END SUBROUTINE LoadSolutionForRestart
      
   END Module DGSEMClass
