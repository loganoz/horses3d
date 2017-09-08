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
      INTEGER, PARAMETER :: HMESH_NONE                          = 0
      INTEGER, PARAMETER :: HMESH_UNDEFINED                     = -1
      INTEGER, PARAMETER :: HMESH_BOUNDARY = 0, HMESH_INTERIOR  = 1
      INTEGER, PARAMETER :: HMESH_NEUMANN  = 1, HMESH_DIRICHLET = 2
      
      CHARACTER(LEN=3)   :: emptyBCName = "---"

      INTEGER, PARAMETER :: BC_STRING_LENGTH = 32

      END Module MeshTypes
