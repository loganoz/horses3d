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

      private
      public Node, ConstructNode, PrintNode
!
!     ---------------
!     Node definition
!     ---------------
!
      TYPE Node
         integer       :: globID = -1
         REAL(KIND=RP) :: x(3)
         contains
            procedure :: construct => ConstructNode
            procedure :: destruct  => Node_Destruct
      END TYPE Node
!
!     ========
      CONTAINS
!     ========
!
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ConstructNode( this, x, globID )
         IMPLICIT NONE 
         class(Node)   :: this
         REAL(KIND=RP) :: x(3)
         integer       :: globID
         this % x  = x
         this % globID = globID
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
      elemental SUBROUTINE Node_Destruct( this )
         IMPLICIT NONE 
         class(Node), intent(inout) :: this
         
         this % GlobID = -1
      END SUBROUTINE Node_Destruct
      
      END Module NodeClass
