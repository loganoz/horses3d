!
!//////////////////////////////////////////////////////
!
!   @File:    pAdaptationClass.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Sun Jan 14 13:23:08 2018
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
MODULE pAdaptationClass
   use SMConstants
   IMPLICIT NONE
!========
 CONTAINS
!========
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
!
   subroutine GetMeshPolynomialOrders(controlVariables,Nx,Ny,Nz,Nmax)
      use FTValueDictionaryClass
      use ReadMeshFile
      implicit none
      !-------------------------------------------------
      type(FTValueDictionary), intent(in)    :: controlVariables
      integer, allocatable   , intent(inout) :: Nx(:), Ny(:), Nz(:)  
      integer                , intent(out)   :: Nmax
      !-------------------------------------------------
      integer                                :: nelem
      !-------------------------------------------------
      
      if (controlVariables % containsKey("polynomial order file")) THEN
         CALL ReadOrderFile( controlVariables % stringValueForKey("polynomial order file", requestedLength = LINE_LENGTH), &
                             Nx, Ny, Nz )
      else
         nelem = NumOfElemsFromMeshFile( controlVariables % stringValueForKey("mesh file name", requestedLength = LINE_LENGTH) )
         allocate( Nx(nelem), Ny(nelem), Nz(nelem) )
         
         if (controlVariables % containsKey("polynomial order")) THEN
            Nx = controlVariables % integerValueForKey("polynomial order")
            Ny = Nx
            Nz = Nx
         else
            if (controlVariables % containsKey("polynomial order i") .AND. &
                controlVariables % containsKey("polynomial order j") .AND. &
                controlVariables % containsKey("polynomial order k") ) THEN
               Nx = controlVariables % integerValueForKey("polynomial order i")
               Ny = controlVariables % integerValueForKey("polynomial order j")
               Nz = controlVariables % integerValueForKey("polynomial order k")
            else
               error stop "The polynomial order(s) must be specified"
            end if
         end if
      end if
      
      Nmax = 0
      if (controlVariables % containsKey("adaptation nmax i")) Nmax = max(Nmax,controlVariables % integerValueForKey("adaptation nmax i"))
      if (controlVariables % containsKey("adaptation nmax j")) Nmax = max(Nmax,controlVariables % integerValueForKey("adaptation nmax j"))
      if (controlVariables % containsKey("adaptation nmax k")) Nmax = max(Nmax,controlVariables % integerValueForKey("adaptation nmax k"))
      
      Nmax = max(Nmax,maxval(Nx),maxval(Ny),maxval(Nz))
      
      end subroutine GetMeshPolynomialOrders
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
