!
!////////////////////////////////////////////////////////////////////////
!
!      Filename.f
!      Created: 2008-05-29 17:44:55 -0400 
!      By: David Kopriva  
!
!////////////////////////////////////////////////////////////////////////
!
      Module MeshTypes
      IMPLICIT NONE
!
!     ---------
!     Constants
!     ---------
!
      INTEGER, PARAMETER :: QMESH_NONE = 0
      INTEGER, PARAMETER :: QMESH_BOUNDARY = 0, QMESH_INTERIOR = 1
      INTEGER, PARAMETER :: QMESH_NEUMANN = 1 , QMESH_DIRICHLET = 2

      INTEGER, PARAMETER :: BC_STRING_LENGTH = 32

      END Module MeshTypes
