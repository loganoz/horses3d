MODULE pAdaptationClass
   IMPLICIT NONE
!========
 CONTAINS
!========

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE ReadOrderFile(filename, Nx, Ny, Nz)
      IMPLICIT NONE
!
!     ----------------------------------------------------------------------
!     Subroutine that reads input file containing polynomial orders for mesh
!     ----------------------------------------------------------------------
!
      CHARACTER(len=*), INTENT(IN) :: filename          !<  Name of file containing polynomial orders to initialize
      INTEGER, ALLOCATABLE         :: Nx(:),Ny(:),Nz(:) !>  Polynomial orders for each element
      !------------------------------------------
      INTEGER                      :: fd       ! File unit
      INTEGER                      :: nelem    ! Number of elements
      INTEGER                      :: i        ! counter
      !------------------------------------------
      
      OPEN(newunit = fd, FILE = filename )   
         READ(fd,*) nelem
         
         ALLOCATE(Nx(nelem),Ny(nelem),Nz(nelem))
         
         DO i = 1, nelem
            READ(fd,*) Nx(i), Ny(i), Nz(i)
         ENDDO
      CLOSE(UNIT=fd)
      
   END SUBROUTINE ReadOrderFile
END MODULE pAdaptationClass
