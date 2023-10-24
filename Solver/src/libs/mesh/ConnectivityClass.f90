#include "Includes.h"
   module ConnectivityClass
      use SMConstants
      IMPLICIT NONE
!
!     -------------------------------------------------------------------
!     This class stores the information about the connectivities between 
!     the different elements in the mesh      
!     -------------------------------------------------------------------
!
      
      TYPE Connectivity
         integer :: globID       ! Global ID of the element on the other side of the face
         integer :: Nxyz(NDIM)   ! Polynomial orders of the element on the other side of the face
      
         CONTAINS
            PROCEDURE           :: construct => ConstructConnectivity
            PROCEDURE           :: destruct  => DestructConnectivity
      END TYPE Connectivity
          
      CONTAINS 
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ConstructConnectivity( this, globID, Nxyz )
         IMPLICIT NONE
         
         class(Connectivity) :: this
         integer, intent(in) :: globID
         integer, intent(in) :: Nxyz(NDIM)
         
         this % globID = globID
         this % Nxyz   = Nxyz
         
      END SUBROUTINE ConstructConnectivity
!
!////////////////////////////////////////////////////////////////////////
!
      elemental SUBROUTINE DestructConnectivity( this )
         IMPLICIT NONE
         CLASS(Connectivity), intent(inout) :: this
         
         this % globID = 0
         this % Nxyz   = -1
         
      END SUBROUTINE DestructConnectivity

      END Module ConnectivityClass