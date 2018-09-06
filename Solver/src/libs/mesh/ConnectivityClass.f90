!
!////////////////////////////////////////////////////////////////////////
!
!      ConnectivityClass.f90
!      Created: 25 de octubre de 2012 13:27 
!      By: Gonzalo Rubio Calzado  
!
!////////////////////////////////////////////////////////////////////////
!
#include "Includes.h"
      Module ConnectivityClass
      IMPLICIT NONE
!
!     -------------------------------------------------------------------
!     This class stores the information about the connectivities between 
!     the different elements in the mesh      
!     -------------------------------------------------------------------
!
      
      TYPE Connectivity
         INTEGER, ALLOCATABLE   :: ElementIDs(:)
      
         CONTAINS
            PROCEDURE           :: construct => ConstructConnectivity
            PROCEDURE           :: destruct  => DestructConnectivity
      END TYPE Connectivity
          
      CONTAINS 
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ConstructConnectivity( this, NumberOfElements )
         IMPLICIT NONE
         
         CLASS(Connectivity)            :: this
         INTEGER                        :: NumberOfElements
         
         ALLOCATE(this%ElementIDs(NumberOfElements))
         this%ElementIDs = 0
         
      END SUBROUTINE ConstructConnectivity
!
!////////////////////////////////////////////////////////////////////////
!
      elemental subroutine DestructConnectivity( this )
         IMPLICIT NONE
         CLASS(Connectivity), intent(inout) :: this
         
         safedeallocate(this%ElementIDs)
         
      END SUBROUTINE DestructConnectivity

!
!////////////////////////////////////////////////////////////////////////
!

!~      SUBROUTINE SetConnectivities(this,fUnit)
!~         IMPLICIT NONE 
!~         INTEGER               :: fUnit, totalElementIDs
!~         INTEGER, ALLOCATABLE  :: temp(:)
!~         TYPE(Connectivity)     :: this(:)
         
!~         INTEGER            :: i, j, counter
         
!~         totalElementIDs = 0
!~         DO i = 1, 4
!~            totalElementIDs = totalElementIDs + SIZE(this(i)%ElementIDs)
!~         ENDDO
         
!~         ALLOCATE(temp(totalElementIDs))
         
!~         READ(fUnit, * ) temp(:)
!~         PRINT*, temp         
!~         counter = 0
!~         DO i = 1,4
!~            DO j = 1,SIZE(this(i)%ElementIDs)
!~               counter = counter + 1
!~               IF (temp(counter) == 0) THEN 
!~                  !PRINT*, i,j
!~                  !STOP 
!~                  temp(counter) = 6
!~               ENDIF
!~               this(i)%ElementIDs(j) = temp(counter)
!~            ENDDO
!~         ENDDO

!~         DEALLOCATE(temp)
         
!~      END SUBROUTINE SetConnectivities 
!
!////////////////////////////////////////////////////////////////////////
!  
      
      END Module ConnectivityClass
