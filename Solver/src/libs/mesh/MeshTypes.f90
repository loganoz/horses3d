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
      integer, parameter :: HMESH_NONE      = 0       ! Not constructed
      integer, parameter :: HMESH_UNDEFINED = -1      ! Constructed but undefined
      integer, parameter :: HMESH_INTERIOR  = 1       ! Interior face
      integer, parameter :: HMESH_BOUNDARY  = 2       ! Physical boundary face
      integer, parameter :: HMESH_MPI       = 3       ! MPI face
      
      CHARACTER(LEN=3), parameter   :: emptyBCName = "---"


      END Module MeshTypes
