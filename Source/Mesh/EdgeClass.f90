!
!////////////////////////////////////////////////////////////////////////
!
!      EdgeClass.f
!      Created: 2008-06-05 14:12:52 -0400 
!      By: David Kopriva  
!
!      Implements Algorithms:
!         Algorithm 125: EdgeClass
!
!////////////////////////////////////////////////////////////////////////
!
      Module EdgeClass
      USE SMConstants
      USE MeshTypes
      IMPLICIT NONE 
!
!     ---------------
!     Edge definition
!     ---------------
!
      TYPE Edge
         INTEGER                         :: edgeType
         INTEGER                         :: bcType          ! QMESH_NEUMANN OR QMESH_DIRICHLET
         INTEGER                         :: nodeIDs(2)
         INTEGER                         :: elementIDs(2)
         INTEGER                         :: elementSide(2)
         INTEGER                         :: nStart, nEnd, nInc
         LOGICAL                         :: maskEdge
         CHARACTER(LEN=BC_STRING_LENGTH) :: boundaryName
         
      END TYPE Edge
!
!     ========
      CONTAINS
!     ========
!
      SUBROUTINE ConstructEdge( this, nodeIDs, elementID, side )
         IMPLICIT NONE 
         TYPE(Edge) :: this
         INTEGER    :: nodeIDs(2), elementID, side
         this%nodeIDS        = nodeIDs
         this%maskEdge       = .TRUE.
         this%bcType         = QMESH_NONE
         this%elementIDs     = QMESH_NONE
         this%elementSide    = QMESH_NONE
         this%elementIDs(1)  = elementID
         this%elementSide(1) = side
         this%nStart         = QMESH_NONE
         this%nEnd           = QMESH_NONE
         this%nInc           = QMESH_NONE
         this%boundaryName   = "none"
      END SUBROUTINE ConstructEdge
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DestructEdge( this )
         IMPLICIT NONE 
         TYPE(Edge) :: this
      END SUBROUTINE DestructEdge
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE PrintEdge( this ) 
      IMPLICIT NONE
      TYPE(Edge) :: this
      PRINT *, "Edge TYPE = "   , this%edgeType
      PRINT *, "Element IDs: "  , this%elementIDs
      PRINT *, "Element Sides: ", this%elementSide
      IF ( this%edgeType == QMESH_INTERIOR )     THEN
         PRINT *, "Neigbor loop data: ", this%nStart, this%nEnd, this%nInc
      ELSE
         PRINT *, "Boundary name = ", this%boundaryName
      END IF
      PRINT *, "-----------------------------------"
      END SUBROUTINE PrintEdge
      
      END Module EdgeClass
