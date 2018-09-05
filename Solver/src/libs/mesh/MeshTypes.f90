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
      use SMConstants
      IMPLICIT NONE
!
!     ---------
!     Constants
!     ---------
!
      private
      public   HMESH_NONE, HMESH_UNDEFINED, HMESH_INTERIOR
      public   HMESH_BOUNDARY, HMESH_MPI, emptyBCName
      public   EFRONT, EBACK, EBOTTOM, ERIGHT, ETOP, ELEFT
      public   SPECMESH, HOPRMESH
   
      public iijjIndexes, coordRotation

      integer, parameter :: EFRONT = 1, EBACK = 2, EBOTTOM = 3
      integer, parameter :: ERIGHT = 4, ETOP = 5, ELEFT = 6

   
      integer, parameter :: HMESH_NONE      = 0       ! Not constructed
      integer, parameter :: HMESH_UNDEFINED = -1      ! Constructed but undefined
      integer, parameter :: HMESH_INTERIOR  = 1       ! Interior face
      integer, parameter :: HMESH_BOUNDARY  = 2       ! Physical boundary face
      integer, parameter :: HMESH_MPI       = 3       ! MPI face
      
      CHARACTER(LEN=3), parameter   :: emptyBCName = "---"
      
      integer, parameter :: SPECMESH = 1
      integer, parameter :: HOPRMESH = 2
      
      contains
!
!////////////////////////////////////////////////////////////////////////
!
!  ROUTINE useD TO COMPUTE FACE ROTATION INDEXES
!     This routine takes indexes on the master Face of a mortar and
!     output the corresponding indexes on the slave Face 
!////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE iijjIndexes(i,j,Nx,Ny,rotation,ii,jj)
      IMPLICIT NONE
      
      integer :: i,j       !<  Input indexes
      integer :: Nx, Ny    !<  Polynomial orders
      integer :: rotation  !<  Face rotation
      integer :: ii,jj     !>  Output indexes
      
      SELECT CASE (rotation)
         CASE (0)
            ii = i
            jj = j
         CASE (1)
            ii = Ny - j
            jj = i
         CASE (2)
            ii = Nx - i
            jj = Ny - j
         CASE (3)
            ii = j
            jj = Nx - i
         CASE (4)
            ii = j
            jj = i
         CASE (5)
            ii = Nx - i
            jj = j
         CASE (6)
            ii = Ny - j
            jj = Nx - i
         CASE (7)
            ii = i
            jj = Ny - j
         CASE DEFAULT 
            PRINT *, "ERROR: Unknown rotation in element faces"
      end SELECT
      
   end SUBROUTINE iijjIndexes

   SUBROUTINE coordRotation(xi,eta,rotation,xiRot, etaRot)
      IMPLICIT NONE
      real(kind=RP), intent(in)   :: xi, eta       ! Master coords 
      integer,       intent(in)   :: rotation      ! Face rotation
      real(kind=RP), intent(out)  :: xiRot, etaRot ! Slave coords
      
      SELECT CASE (rotation)
      CASE (0)
         xiRot  = xi
         etaRot = eta
      CASE (1)
         xiRot  = 1.0_RP - eta
         etaRot = xi
      CASE (2)
         xiRot  = 1.0_RP - xi
         etaRot = 1.0_RP - eta
      CASE (3)
         xiRot  = eta
         etaRot = 1.0_RP - xi
      CASE (4)
         xiRot  = eta
         etaRot = xi
      CASE (5)
         xiRot  = 1.0_RP - xi
         etaRot = eta
      CASE (6)
         xiRot  = 1.0_RP - eta
         etaRot = 1.0_RP - xi
      CASE (7)
         xiRot  = xi
         etaRot = 1.0_RP - eta
      CASE DEFAULT 
         PRINT *, "ERROR: Unknown rotation in element faces"
      end SELECT
      
   end SUBROUTINE coordRotation



      END Module MeshTypes
