!
!////////////////////////////////////////////////////////////////////////
!
!      NodeClass.f
!      Created: 2008-06-05 14:29:53 -0400 
!      By: David Kopriva  
!
!      Implements Algorithms:
!         Algorithm 123: Node (CornerNode)
!         
!////////////////////////////////////////////////////////////////////////
!
      Module NodeClass 
      USE SMConstants
      USE LinkedListClass
      USE MeshTypes
      IMPLICIT NONE 
!
!     ---------------
!     Node definition
!     ---------------
!
      TYPE Node
         INTEGER                         :: nodeType ! QMESH_INTERIOR or QMESH_BOUNDARY
         INTEGER                         :: bcType   ! QMESH_NEUMANN  or QMESH_DIRICHLET
         REAL(KIND=RP)                   :: x(2)
         TYPE(LinkedList)                :: adjElements
         CHARACTER(LEN=BC_STRING_LENGTH) :: boundaryName
      END TYPE Node
!
!     ========
      CONTAINS
!     ========
!
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ConstructNode( this, x )
         IMPLICIT NONE 
         TYPE(Node)    :: this
         REAL(KIND=RP) :: x(2)
         this%x            = x
         this%nodeType     = QMESH_INTERIOR
         this%bcType       = QMESH_NONE ! QMESH_NEUMANN or QMESH_DIRICHLET
         this%boundaryName = "none"
         CALL ConstructList(this%adjElements)
      END SUBROUTINE ConstructNode
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE PrintNode( this, id )
         IMPLICIT NONE 
         TYPE(Node) :: this
         INTEGER    :: id
         PRINT *, id, this%x, this%nodeType, this%BoundaryName
         CALL PRINT( this%adjElements )
      END SUBROUTINE PrintNode
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DestructNode( this )
         IMPLICIT NONE 
         TYPE(Node) :: this
         CALL Destruct( this%adjElements )
      END SUBROUTINE DestructNode
      
      END Module NodeClass
