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
!      Modified for 3D DG 5/27/15, 12:18 PM
!
!////////////////////////////////////////////////////////////////////////
!
      Module NodeClass 
      USE SMConstants
      USE MeshTypes
      IMPLICIT NONE 
!
!     ---------------
!     Node definition
!     ---------------
!
      TYPE Node
         REAL(KIND=RP)                   :: x(3)
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
         REAL(KIND=RP) :: x(3)
         this%x        = x
      END SUBROUTINE ConstructNode
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE PrintNode( this, id )
         IMPLICIT NONE 
         TYPE(Node) :: this
         INTEGER    :: id
         PRINT *, id, this%x
      END SUBROUTINE PrintNode
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DestructNode( this )
         IMPLICIT NONE 
         TYPE(Node) :: this
      END SUBROUTINE DestructNode
      
      END Module NodeClass
