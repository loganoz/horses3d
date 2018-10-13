!
!//////////////////////////////////////////////////////
!
!   @File: ColorsClass.f90
!   @Author:  Carlos Redondo (module for 2D) and Andrés Rueda  (am.rueda@upm.es - 3D implementation and changes) 
!   @Created: Tue Mar 17 17:05:00 2017
!   @Last revision date: Sat Oct 13 21:14:10 2018
!   @Last revision author: Andrés Rueda (am.rueda@upm.es)
!   @Last revision commit: e5b253a108cfb11371774c409a869acee046f770
!
!//////////////////////////////////////////////////////
!
!      Module for computing element colorings in order to generate the numerical Jacobian
!
!////////////////////////////////////////////////////////////////////////
MODULE ColorsClass
   use HexMeshClass
   implicit none
   
   private
   public Colors_t
   
   type Colors_t
         INTEGER                             :: num_of_colors ! number of colors
         INTEGER,DIMENSION(:), ALLOCATABLE   :: elmnts  ! color ordered elements
         INTEGER,DIMENSION(:), ALLOCATABLE   :: bounds  ! idx of the first element on a color
         INTEGER                             :: ntotal  ! total  numer of elements
      CONTAINS
         procedure :: construct  => Colors_Contruct !construct(nbr)
         procedure :: destruct   => Colors_Destruct
         procedure :: info       => getcolorsinfo   ! prints coloring info  
         procedure :: export2Tec => ExportColors2Tec
   END TYPE Colors_t
   
   integer, parameter :: NUM_OF_NEIGHBORS = 6 ! Hardcoded: Hexahedral conforming meshes
   
   contains
      subroutine Colors_Contruct(this, nbr,depth)
         implicit none
         !-arguments------------------------------------------------------
         class(Colors_t), intent(out)        :: this
         type(neighbour), intent(in)         :: nbr(:)
         integer        , intent(in)         :: depth
         !-local-variables------------------------------------------------
         INTEGER                             :: ncolored = 0
         LOGICAL, DIMENSION(:), ALLOCATABLE  :: colored, used
         LOGICAL                             :: allcolored = .FALSE.
         INTEGER                             :: i, j, counter, idx
         INTEGER                             :: ntotal, maxcolor
         INTEGER, DIMENSION(:), ALLOCATABLE  :: colors
         !----------------------------------------------------------------
         
         ntotal = SIZE(nbr)
         this%ntotal = ntotal
         ALLOCATE(used(0:ntotal)) !0 correspond to boundary "neighbour"
         ALLOCATE(colored(ntotal))
         ALLOCATE(colors(ntotal))
         colored(:) = .FALSE.
         maxcolor = 0
         
!        Create colors and assign elements
!        *********************************
         
         do while (.NOT. allcolored)
            used = .FALSE.
            maxcolor = maxcolor + 1
            do i = 1, ntotal
               
!              Check if current element is colored
!              -----------------------------------
               if ( colored(i) ) cycle
               
!              Check if its neighbors(...depth times) were used in this color
!              --------------------------------------------------------------  
               used(0) = .FALSE. !boundary neigbour is empty
               if (neighbors_were_used(used,nbr,i,depth)) cycle
               
!              Mark neighbors as used
!              ----------------------
               call mark_neighbors_as_used(used,nbr,i,depth)
               
!              Color this element
!              ------------------
               colored(i) = .TRUE.
               colors(i) = maxcolor
               ncolored = ncolored + 1
               if (ncolored == ntotal) allcolored = .TRUE.
               
            end do
         end do
         
         
         this % num_of_colors = maxcolor
         ALLOCATE(this%elmnts(ntotal))
         ALLOCATE(this%bounds(this % num_of_colors + 1))
         
!        Order elements acording to colors
!        *********************************
         idx = 1
         DO i = 1, this % num_of_colors
            this%bounds(i)= idx
            counter = 0
            DO j = 1, ntotal
               IF (colors(j) == i) THEN
                  this%elmnts(idx) = j
                  counter = counter + 1
                  idx = idx + 1
               end if      
            end do
         end do
         this%bounds(i)= idx
         
      end subroutine Colors_Contruct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Colors_Destruct(this)
         implicit none
         class(Colors_t), intent(out)        :: this
         
         this % num_of_colors = 0
         this % ntotal        = 0
         
         deallocate (this % elmnts)
         deallocate (this % bounds)
         
      end subroutine Colors_Destruct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE getcolorsinfo(this)
         CLASS(Colors_t)        :: this
         INTEGER              :: i
       
         WRITE(*,'(A13,I2)') "# of colors: ", this % num_of_colors
         WRITE(*,*) "Element list:"         
         WRITE(*,'(*(I5,1X))') (this%elmnts(i), i = 1,this%ntotal)
         WRITE(*,*) "New color indexes"
         WRITE(*,'(*(I5,1X))') (this%bounds(i), i= 1,this % num_of_colors + 1)
      END subroutine getcolorsinfo
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE ExportColors2Tec(this,mesh,filename)
      IMPLICIT NONE
      !-------------------------------------------------------
      CLASS(Colors_t)                     :: this
      TYPE(HexMesh)                       :: mesh
      character(len=*), intent(in)        :: filename
      
      !-------------------------------------------------------
      INTEGER                             :: fd, Nelem, id, N(3), i, j, k
      INTEGER, DIMENSION(:), ALLOCATABLE  :: colors
      !-------------------------------------------------------
      open(newunit = fd, file=trim(filename), action='WRITE')
      
      Nelem = SIZE(mesh % elements)
      
      !! Create colors array
      ALLOCATE(colors(Nelem))
      
      DO i=1, this % num_of_colors
         DO j=this % bounds(i), this % bounds(i+1)-1
            colors(this % elmnts(j)) = i
         end do
      end do
      
      write(fd,*) 'TITLE = "Colors of mesh NSLITE3D"'
      write(fd,*) 'VARIABLES = "X","Y","Z","Color" '
      
      DO id = 1, Nelem
         N = mesh % elements(id) % Nxyz
         WRITE(fd,*) 'ZONE I=', N(1)+1, ",J=",N(2)+1, ",K=",N(3)+1,", F=POINT"
         
         DO k = 0, N(3)
            DO j= 0, N(2)
               DO i = 0, N(1)
                  write(fd,'(3E13.5,x,I0)') mesh % elements(id) % geom % x(1,i,j,k), &
                                            mesh % elements(id) % geom % x(2,i,j,k), &
                                            mesh % elements(id) % geom % x(3,i,j,k), &
                                            colors(id)
               end do
            end do
         end do
      end do
      
      close(fd)
      
   END SUBROUTINE ExportColors2Tec
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  SOME EXTRA UTILITIES
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------------------------
!  Check if any of the neighbors [(depth-1) * "of neighbors"] was used 
!  -------------------------------------------------------------------
   recursive function neighbors_were_used(used,nbr,eID,depth) result (were_used)
      implicit none
      !-arguments-------------------------------------------------
      logical        , intent(in) :: used(0:)
      type(Neighbour), intent(in) :: nbr(:)
      integer        , intent(in) :: eID
      integer        , intent(in) :: depth
      logical                     :: were_used
      !-local-variables-------------------------------------------
      integer :: n_eID
      !-----------------------------------------------------------
      
      if (eID == 0) then
         were_used = .FALSE.
         return
      end if
      
      were_used = any ( used(nbr(eID) % elmnt) )
      if (were_used) return
      
      if (depth > 1) then
         do n_eID = 1, NUM_OF_NEIGHBORS
            were_used = neighbors_were_used(used, nbr, nbr(eID) % elmnt(n_eID), depth-1)
            if (were_used) return
         end do
      end if
      
   end function neighbors_were_used
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------------------------------------------------
!  Mark all the nighbors [(depth-1) * "of neighbors"] as used 
!  ----------------------------------------------------------
   recursive subroutine mark_neighbors_as_used(used,nbr,eID,depth)
      implicit none
      !-arguments-------------------------------------------------
      logical        , intent(inout) :: used(0:)
      type(Neighbour), intent(in)    :: nbr(:)
      integer        , intent(in)    :: eID
      integer        , intent(in)    :: depth
      !-local-variables-------------------------------------------
      integer :: n_eID
      !-----------------------------------------------------------
      
      if (eID == 0) return
      used (nbr(eID)%elmnt) = .TRUE.
      
      if (depth > 1) then
         do n_eID = 1, NUM_OF_NEIGHBORS
            call mark_neighbors_as_used (used, nbr, nbr(eID)%elmnt(n_eID), depth-1)
         end do
      end if
      
   end subroutine mark_neighbors_as_used
   
END MODULE ColorsClass
