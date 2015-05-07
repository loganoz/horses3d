!
!////////////////////////////////////////////////////////////////////////
!
!      DGSolutionStorage.f95
!      Created: 2008-07-13 15:56:15 -0400 
!      By: David Kopriva  
!
!     Algorithms:
!        Algorithm 110: DGSolutionStorage
!
!////////////////////////////////////////////////////////////////////////
!
      Module DGSolutionStorageClass
      USE SMConstants
      IMPLICIT NONE
!
!     ----------
!     Class type
!     ----------
!
      TYPE DGSolutionStorage
         INTEGER :: nEqn
!
!        ---------------------------------------
!        Interior solutions and time derivatives
!        ---------------------------------------
!
         REAL(KIND=RP), DIMENSION(:,:,:), ALLOCATABLE :: Q, QDot, G
         REAL(KIND=RP), DIMENSION(:,:,:), ALLOCATABLE :: U_x, U_y
!
!        -------------------------------------------------------------
!        Boundary values of: The solution, the inviscid Riemann flux, 
!        the viscous riemann flux
!        -------------------------------------------------------------
!
         REAL(KIND=RP), DIMENSION(:,:,:), ALLOCATABLE :: Qb, Ub, U_xb, U_yb, FStarb
      END TYPE DGSolutionStorage
!
!     ========
      CONTAINS 
!     ========
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ConstructDGSolutionStorage( this, N, nEqn, nGradEqn, flowIsNavierStokes )
      IMPLICIT NONE 
      TYPE(DGSolutionStorage) :: this
      INTEGER                 :: N, nEqn, nGradEqn
      LOGICAL                 :: flowIsNavierStokes
      
      this%nEqn = nEqn
!
!     --------------------------------------
!     Solution and time derivative variables
!     --------------------------------------
!
      ALLOCATE( this%Q(0:N,0:N,nEqn), this%QDot(0:N,0:N,nEqn), this%G(0:N,0:N,nEqn) )
      IF ( flowIsNavierStokes )     THEN
         ALLOCATE( this%U_x(0:N,0:N,nGradEqn), this%U_y(0:N,0:N,nGradEqn) )
      END IF
!
!     ---------------
!     Boundary values
!     ---------------
!
      ALLOCATE( this%Qb(nEqn,0:N,4) )
      IF ( flowIsNavierStokes )     THEN
         ALLOCATE( this%U_xb(nGradEqn,0:N,4), this%U_yb(nGradEqn,0:N,4) )
         ALLOCATE( this%Ub(nGradEqn,0:N,4) )
      END IF

      ALLOCATE( this%FStarb(nEqn,0:N,4) )
!
      
      this%G           = 0.0_RP
      this%Q           = 0.0_RP
      this%QDot        = 0.0_RP
      this%Qb          = 0.0_RP
      
      IF ( flowIsNavierStokes )     THEN
         this%Ub          = 0.0_RP
         this%U_x         = 0.0_RP
         this%U_y         = 0.0_RP
         this%U_xb        = 0.0_RP
         this%U_yb        = 0.0_RP
      END IF
      
      END SUBROUTINE ConstructDGSolutionStorage
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DestructDGSolutionStorage( this )
      IMPLICIT NONE 
      TYPE(DGSolutionStorage) :: this
      
      DEALLOCATE( this%Q, this%QDot, this%G )
      DEALLOCATE( this%Qb, this%FStarb )
      IF ( ALLOCATED(this%Ub) )     THEN
         DEALLOCATE( this%Ub, this%U_x, this%U_y, this%U_xb, this%U_yb )
      END IF

      
      END SUBROUTINE DestructDGSolutionStorage
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE SaveDGSolutionStorageToUnit( this, fUnit )
         IMPLICIT NONE
!
!        -----------------------
!        Save for a restart file
!        -----------------------
!
         TYPE(DGSolutionStorage) :: this
         INTEGER                 :: fUnit
         
         WRITE(funit) this%Q
      
      END SUBROUTINE SaveDGSolutionStorageToUnit
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE LoadDGSolutionStorageFromUnit( this, fUnit )
         IMPLICIT NONE
!
!        -----------------------
!        Save for a restart file
!        -----------------------
!
         TYPE(DGSolutionStorage) :: this
         INTEGER                 :: fUnit
         
         READ(funit) this%Q
      
      END SUBROUTINE LoadDGSolutionStorageFromUnit

      
      END Module DGSolutionStorageClass
