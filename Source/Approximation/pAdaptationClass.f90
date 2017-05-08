MODULE pAdaptationClass
   IMPLICIT NONE
!========
 CONTAINS
!========

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE ReadOrderFile(filename, N)
      IMPLICIT NONE
!
!     ----------------------------------------------------------------------
!     Subroutine that reads input file containing polynomial orders for mesh
!     ----------------------------------------------------------------------
!
      CHARACTER(len=*), INTENT(IN) :: filename !<  Name of file containing polynomial orders to initialize
      INTEGER, ALLOCATABLE         :: N(:,:)   !>  Polynomial orders for each element
      !------------------------------------------
      INTEGER                      :: fd       ! File unit
      INTEGER                      :: nelem    ! Number of elements
      INTEGER                      :: i        ! counter
      !------------------------------------------
      
      OPEN(newunit = fd, FILE = filename )   
         READ(fd,*) nelem
         
         ALLOCATE(N(nelem,3))
         
         DO i = 1, nelem
            READ(fd,*) N(i,1), N(i,2), N(i,3)
         ENDDO
      CLOSE(UNIT=fd)
      
   END SUBROUTINE ReadOrderFile
END MODULE pAdaptationClass
