!
!////////////////////////////////////////////////////////////////////////
!
!      ElementClass.f95
!      Created: 2008-06-04 15:34:44 -0400 
!      By: David Kopriva
!
!      Implements Algorithms:
!         Algorithm: 124: ElementClass (QuadElementClass)
!
!       The Quad Element class, Alg. 124. See Sec. 8.2.1.2. This has
!       been modified to add the association of a boundary name to an element
!       edge so that different boundary conditions can be applied to different
!       elements. The names of the boundaries (not necessarily the names of the
!       *boundary conditions* to be applied) are of length BC_STRING_LENGTH.
!       One will associate boundary conditions to boundaries in the routine
!       "ExternalState".
!
!       Modified 2D Code to move solution into element class. 5/14/15, 5:36 PM
!
!////////////////////////////////////////////////////////////////////////
!
      Module ElementClass
      USE SMConstants
      USE PolynomialInterpAndDerivsModule
      USE GaussQuadrature
      USE TransfiniteMapClass
      USE MappedGeometryClass
      USE MeshTypes
      IMPLICIT NONE
      
      
      TYPE Element
          INTEGER                                      :: nodeIDs(4)
          INTEGER                                      :: N
          TYPE(MappedGeometry)                         :: geom
          REAL(KIND=RP), DIMENSION(:,:,:), ALLOCATABLE :: Q, QDot, G
          REAL(KIND=RP), DIMENSION(:,:,:), ALLOCATABLE :: U_x, U_y
!
!         -------------------------------------------------------------
!         Boundary values of: The solution, the inviscid Riemann flux, 
!         the viscous riemann flux
!         -------------------------------------------------------------
!
          REAL(KIND=RP), DIMENSION(:,:,:), ALLOCATABLE :: Qb, Ub, U_xb, U_yb, FStarb
          CHARACTER(LEN=BC_STRING_LENGTH)              :: boundaryName(4)
      END TYPE Element 
      
      INTEGER, PARAMETER :: edgeMap(2,4) = RESHAPE( (/1,2,2,3,4,3,1,4/), (/2,4/) )
      
      CONTAINS 
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ConstructElement( this, ng, nodeIDs, quadMap )
         USE Nodal2DStorageClass
         IMPLICIT NONE
         
         TYPE(Element)            :: this
         TYPE(Nodal2DStorage)     :: ng
         INTEGER                  :: nodeIDs(4)
         TYPE(TransfiniteQuadMap) :: quadMap
         INTEGER                  :: N
         
         this%nodeIDs      = nodeIDs
         this%N            = ng%N
         N                 = ng%N 
         this%boundaryName = "---"
!
!        --------
!        Geometry
!        --------
!
         CALL ConstructMappedGeometry( this%geom, ng, this%N, this%N, quadMap )
!
!        ----------------------------------------
!        Solution Storage is allocated separately
!        ----------------------------------------
!
         
      END SUBROUTINE ConstructElement
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE allocateElementStorage(this, N, nEqn, nGradEqn, flowIsNavierStokes)  
         IMPLICIT NONE
         TYPE(Element) :: this
         INTEGER       :: N, nEqn, nGradEqn
         LOGICAL       :: flowIsNavierStokes
!
!        --------------------------------------
!        Solution and time derivative variables
!        --------------------------------------
!
         ALLOCATE( this%Q(0:N,0:N,nEqn), this%QDot(0:N,0:N,nEqn), this%G(0:N,0:N,nEqn) )
         IF ( flowIsNavierStokes )     THEN
            ALLOCATE( this%U_x(0:N,0:N,nGradEqn), this%U_y(0:N,0:N,nGradEqn) )
         END IF
!
!        ---------------
!        Boundary values
!        ---------------
!
         ALLOCATE( this%Qb(nEqn,0:N,4) )
         IF ( flowIsNavierStokes )     THEN
            ALLOCATE( this%U_xb(nGradEqn,0:N,4), this%U_yb(nGradEqn,0:N,4) )
            ALLOCATE( this%Ub(nGradEqn,0:N,4) )
         END IF
   
         ALLOCATE( this%FStarb(nEqn,0:N,4) )
!
!        -----------------
!        Initialize memory
!        -----------------
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

      END SUBROUTINE allocateElementStorage
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE SetElementBoundaryNames( this, names ) 
         IMPLICIT NONE
         TYPE(Element)                   :: this
         CHARACTER(LEN=BC_STRING_LENGTH) :: names(4)
         this%boundaryName = names
      END SUBROUTINE SetElementBoundaryNames
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DestructElement( this )
         IMPLICIT NONE
         TYPE(Element) :: this
         
         CALL DestructMappedGeometry( this%geom )
         
         DEALLOCATE( this%Q, this%QDot, this%G )
         DEALLOCATE( this%Qb, this%FStarb )
         
         IF ( ALLOCATED(this%Ub) )     THEN
            DEALLOCATE( this%Ub, this%U_x, this%U_y, this%U_xb, this%U_yb )
         END IF

      END SUBROUTINE DestructElement
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE PrintElement( this, id )
         IMPLICIT NONE 
         TYPE(Element) :: this
         INTEGER      :: id
         PRINT *, id, this%nodeIDs, this%boundaryName
      END SUBROUTINE PrintElement
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE SaveSolutionStorageToUnit( this, fUnit )
         IMPLICIT NONE
!
!        -----------------------
!        Save for a restart file
!        -----------------------
!
         TYPE(Element) :: this
         INTEGER       :: fUnit
         
         WRITE(funit) this%Q
      
      END SUBROUTINE SaveSolutionStorageToUnit
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE LoadSolutionFromUnit( this, fUnit )
         IMPLICIT NONE
!
!        -----------------------
!        Save for a restart file
!        -----------------------
!
         TYPE(Element) :: this
         INTEGER       :: fUnit
         
         READ(funit) this%Q
      
      END SUBROUTINE LoadSolutionFromUnit
      
      END Module ElementClass

