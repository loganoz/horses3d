!
!////////////////////////////////////////////////////////////////////////
!
!      FaceClass.f
!      Created: 2008-06-05 14:12:52 -0400 
!      By: David Kopriva  
!
!      Modified to 3D 5/27/15, 11:13 AM
!
!      Implements Algorithms:
!         Algorithm 125: EdgeClass -> 3D
!
!      A face simply keeps track of which elements share a face and 
!      how they are oriented.
!
!////////////////////////////////////////////////////////////////////////
!
      Module FaceClass
      USE SMConstants
      USE MeshTypes
      IMPLICIT NONE 
!
!     ---------------
!     Face definition
!     ---------------
!
      TYPE Face
         INTEGER                         :: FaceType
         INTEGER                         :: nodeIDs(4)
         INTEGER                         :: elementIDs(2)
         INTEGER                         :: elementSide(2)
         INTEGER                         :: rotation
         CHARACTER(LEN=BC_STRING_LENGTH) :: boundaryName
      END TYPE Face
!
!     ========
      CONTAINS
!     ========
!
      SUBROUTINE ConstructFace( self, nodeIDs, elementID, side )
         IMPLICIT NONE 
         TYPE(Face) :: self
         INTEGER    :: nodeIDs(4), elementID, side
         self % FaceType       = HMESH_UNDEFINED
         self % nodeIDS        = nodeIDs
         self % elementIDs     = HMESH_NONE
         self % elementSide    = HMESH_NONE
         self % elementIDs(1)  = elementID
         self % elementSide(1) = side
         self % boundaryName   = emptyBCName
         self % rotation       = 0
      END SUBROUTINE ConstructFace
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DestructFace( self )
         IMPLICIT NONE 
         TYPE(Face) :: self
      END SUBROUTINE DestructFace
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE PrintFace( self ) 
      IMPLICIT NONE
      TYPE(Face) :: self
      PRINT *, "Face TYPE = "   , self % FaceType
      PRINT *, "Element IDs: "  , self % elementIDs
      PRINT *, "Element Sides: ", self % elementSide
      IF ( self % FaceType == HMESH_INTERIOR )     THEN
         PRINT *, "Neighbor rotation: ", self  %  rotation
      ELSE
         PRINT *, "Boundary name = ", self % boundaryName
      END IF
      PRINT *, "-----------------------------------"
      END SUBROUTINE PrintFace
!
!////////////////////////////////////////////////////////////////////////
!
!  ROUTINE USED TO COMPUTE FACE ROTATION INDEXES
!     TODO: Check if this is enough or if one needs 8 indexes!!
!////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE iijjIndexes(i,j,N,rotation,ii,jj)
      IMPLICIT NONE
      
      INTEGER :: i,j       !<  Input indexes
      INTEGER :: N         !<  Polynomial order
      INTEGER :: rotation  !<  Face rotation
      INTEGER :: ii,jj     !>  Output indexes
      
      SELECT CASE (rotation)
         CASE (0)
            ii = i
            jj = j
         CASE (1) ! is this o.k?
            ii = N - i
            jj = j
         CASE (2)
            ii = N - i
            jj = N - j
         CASE (3)
            ii = j
            jj = N - i
         CASE DEFAULT 
            PRINT *, "ERROR: Unknown rotation in element faces"
      END SELECT
      
   END SUBROUTINE iijjIndexes
!
!////////////////////////////////////////////////////////////////////////
!
      END Module FaceClass
