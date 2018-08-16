!
!//////////////////////////////////////////////////////
!
!   @File:    ConnectivityClass.f90
!   @Author:  Gonzalo Rubio Calzado (g.rubio@upm.es)
!   @Created: Thu Oct  25 13:27:17 2012
!   @Last revision date: Thu Aug 16 16:15:29 2018
!   @Last revision author: Andrés Rueda (am.rueda@upm.es)
!   @Last revision commit: d9871e8d2a08e4b4346bb29d921b80d139c575cd
!
!//////////////////////////////////////////////////////
!
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
      elemental SUBROUTINE DestructConnectivity( this )
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
