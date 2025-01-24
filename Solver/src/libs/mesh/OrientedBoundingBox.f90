#include "Includes.h"
module OrientedBoundingBox

   use SMConstants
   use Utilities
   use TessellationTypes

   implicit none   
   
   integer,       parameter :: TASK_THRESHOLD = 10000, LOCAL = 0, GLOBAL =1, BOXVERTICES = 8
   real(kind=RP), parameter :: SAFETY_FACTOR = 0.001_RP

   public :: OBB

!  **************************************************
!  Main type for the Oriented Bounding Box computations
!  **************************************************
   type OBB_type
      
      real(kind=rp) :: vertices(NDIM,BOXVERTICES), GlobalVertices(NDIM,BOXVERTICES)
      integer       :: maxAxis, minAxis
      
      contains
         procedure :: construct         => OBB_construct
         procedure :: isPointInside     => OBB_isPointInside
         procedure :: plot              => OBB_plot
         procedure :: GetMinAxis        => OBB_GetMinAxis
         procedure :: GetGLobalVertices => OBB_GetGLobalVertices
   end type
   
   type(OBB_type),  allocatable :: OBB(:)

contains
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!   
!  -------------------------------------------------------------------------------------------
!  Subroutine for plotting the OBB
!  -------------------------------------------------------------------------------------------
   subroutine OBB_plot( this, filename, iter )
      use PhysicsStorage
      use MPI_Process_Info
      implicit none
      !-arguments--------------------------
      class(OBB_type),  intent(inout) :: this
      character(len=*), intent(in)    :: filename  
      integer,          intent(in)    :: iter 
      !-local-variables----------------------
      integer                    :: i, funit
      character(len=LINE_LENGTH) :: it 

      if( .not. MPI_Process% isRoot ) return

      funit = UnusedUnit()
      
      write(it,*) iter

      open(funit,file='IBM/OrientedBoundingBox_'//trim(filename)//'_'//trim(adjustl(it))//'.tec', status='unknown')

      write(funit,"(a28)") 'TITLE = "OrientedBoudingBox"'
      write(funit,"(a25)") 'VARIABLES = "x", "y", "z"'
      write(funit,"(a69)") 'ZONE NODES=8, ELEMENTS = 6, DATAPACKING=POINT, ZONETYPE=FETETRAHEDRON'

      do i = 1, 8
         write(funit,'(E17.10,1X,E17.10,1X,E17.10)') Lref*this% vertices(1,i), Lref*this% vertices(2,i), Lref*this% vertices(3,i)
      end do 

      write(funit,'(4i2)') 1, 2, 3, 4
      write(funit,'(4i2)') 1, 5, 8, 4
      write(funit,'(4i2)') 5, 6, 7, 8
      write(funit,'(4i2)') 2, 3, 7, 6
      write(funit,'(4i2)') 4, 8, 7, 3
      write(funit,'(4i2)') 1, 2, 6, 5

      close(unit=funit)

   end subroutine OBB_plot
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine computes the Oriented Bounding Box. 
!  -------------------------------------------------
   subroutine OBB_construct( this, stl, isPlot )
   
      implicit none
      !-arguments-----------------------------------
      class(OBB_type), intent(inout) :: this
      type(STLfile),   intent(in)    :: stl
      logical,         intent(in)    :: isPlot
      !-local-variables-----------------------------
      integer              :: i, j 
      real(kind=RP)        :: max_x, max_y, max_z
      real(kind=RP)        :: min_x, min_y, min_z

      max_x = -huge(1.0_RP); max_y = -huge(1.0_RP); max_z = -huge(1.0_RP)
      min_x =  huge(1.0_RP); min_y =  huge(1.0_RP); min_z =  huge(1.0_RP)

      do i = 1, stl% NumOfObjs
         do j = 1, NumOfVertices
            max_x = max(max_x, stl% ObjectsList(i)% vertices(j)% coords(IX))
            max_y = max(max_y, stl% ObjectsList(i)% vertices(j)% coords(IY))
            max_z = max(max_z, stl% ObjectsList(i)% vertices(j)% coords(IZ))
            min_x = min(min_x, stl% ObjectsList(i)% vertices(j)% coords(IX))
            min_y = min(min_y, stl% ObjectsList(i)% vertices(j)% coords(IY))
            min_z = min(min_z, stl% ObjectsList(i)% vertices(j)% coords(IZ))
         end do 
      end do 

      this% vertices(:,1) = (/min_x, min_y, min_z/)
      this% vertices(:,2) = (/max_x, min_y, min_z/)
      this% vertices(:,3) = (/max_x, max_y, min_z/)
      this% vertices(:,4) = (/min_x, max_y, min_z/)
      this% vertices(:,5) = (/min_x, min_y, max_z/)
      this% vertices(:,6) = (/max_x, min_y, max_z/)
      this% vertices(:,7) = (/max_x, max_y, max_z/)
      this% vertices(:,8) = (/min_x, max_y, max_z/)

      if( isPlot ) call this% plot( stl% NoExtfilename, 0)

   end subroutine OBB_construct
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This function check if a point is inside the OBB or not
!  -------------------------------------------------   
   function OBB_isPointInside( this, Point, coeff ) result( isInsideOBB )
   
      implicit none
      !-arguments------------------------------------
      class(OBB_type),         intent(inout) :: this
      real(kind=rp),           intent(in)    :: Point(NDIM)
      real(kind=rp), optional, intent(in)    :: coeff
      logical                                :: isInsideOBB
      !-local-variables--------------------------------
      real(kind=rp) :: Center(NDIM), Vertices(NDIM,BOXVERTICES)
      real(kind=rp) :: Multcoeff
      integer       :: i      
      
      if( present(coeff) ) then
         Multcoeff = coeff
      else
         Multcoeff = 1.0_RP
      end if
      
      isInsideOBB = .false.

      Center(IX) = sum(this% GlobalVertices(IX,:))/BOXVERTICES
      Center(IY) = sum(this% GlobalVertices(IY,:))/BOXVERTICES
      Center(IZ) = sum(this% GlobalVertices(IZ,:))/BOXVERTICES
      
      do i = 1, BOXVERTICES
         Vertices(IX,i) = Center(IX) + Multcoeff*(this% GlobalVertices(IX,i) - Center(IX))
         Vertices(IY,i) = Center(IY) + Multcoeff*(this% GlobalVertices(IY,i) - Center(IY))
         Vertices(IZ,i) = Center(IZ) + Multcoeff*(this% GlobalVertices(IZ,i) - Center(IZ))
      end do
      !check x y z     
       if( isInsideBox(Point, Vertices) ) isInsideOBB = .true.

   end function OBB_isPointInside

   subroutine OBB_GetGlobalVertices( this )
      use MPI_Utilities
      implicit none 

      class(OBB_type), intent(inout) :: this 
      !-local-variables--------------------------------
      real(kind=rp) :: Vertices(NDIM,BOXVERTICES)
      integer       :: i      
      
      this% GlobalVertices = this% Vertices

      call MPI_MinMax(this% GlobalVertices(IX,1), this% GlobalVertices(IX,7))
      call MPI_MinMax(this% GlobalVertices(IY,1), this% GlobalVertices(IY,7))
      call MPI_MinMax(this% GlobalVertices(IZ,1), this% GlobalVertices(IZ,7))

      this% GlobalVertices(:,2) = (/this% GlobalVertices(IX,7), this% GlobalVertices(IY,1), this% GlobalVertices(IZ,1)/)
      this% GlobalVertices(:,6) = (/this% GlobalVertices(IX,7), this% GlobalVertices(IY,1), this% GlobalVertices(IZ,7)/)

      this% GlobalVertices(:,3) = (/this% GlobalVertices(IX,7), this% GlobalVertices(IY,7), this% GlobalVertices(IZ,1)/)
      this% GlobalVertices(:,5) = (/this% GlobalVertices(IX,1), this% GlobalVertices(IY,1), this% GlobalVertices(IZ,7)/)

      this% GlobalVertices(:,4) = (/this% GlobalVertices(IX,1), this% GlobalVertices(IY,7), this% GlobalVertices(IZ,1)/)
      this% GlobalVertices(:,8) = (/this% GlobalVertices(IX,1), this% GlobalVertices(IY,7), this% GlobalVertices(IZ,7)/)

   end subroutine OBB_GetGlobalVertices

   integer function OBB_GetMinAxis( this, maxAxis, clipAxis )

      implicit none 

      class(OBB_type), intent(inout) :: this 
      integer,         intent(in)    :: maxAxis, clipAxis 

      real(kind=RP) :: dx, dy, dz 
      integer       :: i, vec(NDIM), minAxis

      dx = abs( this% GlobalVertices(IX,7) - this% GlobalVertices(IX,1) )
      dy = abs( this% GlobalVertices(IY,7) - this% GlobalVertices(IY,1) )
      dz = abs( this% GlobalVertices(IZ,7) - this% GlobalVertices(IZ,1) )

      minAxis = minloc((/dx,dy,dz/),dim=1)

      vec = 0
      do i = 1, NDIM
         if( i .eq. ClipAxis ) cycle  
         if( i .eq. maxAxis  ) cycle
         vec(i)  = i 
      end do

      if( vec(minAxis) .eq. 0 ) then 
         do i = 1, NDIM 
            if( vec(i) .ne. 0 ) minAxis = i 
         end do 
      end if 

      OBB_GetMinAxis = minAxis

   end function OBB_GetMinAxis
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This function check if a point is inside a generic box
!  -------------------------------------------------   

   logical function isInsideBox( Point, vertices, equal ) result( isInside )
   
      implicit none

      real(kind=rp), dimension(:),   intent(in) :: Point
      real(kind=rp), dimension(:,:), intent(in) :: vertices
      logical,                       intent(in) :: equal
      
      optional :: equal

      isInside = .false.

      if( present(equal) ) then
         if( .not. equal ) then
            if( (Point(1) > vertices(1,1) .and. Point(1) < vertices(1,7)) .and. &
                (Point(2) > vertices(2,1) .and. Point(2) < vertices(2,7)) .and. &
                (Point(3) > vertices(3,1) .and. Point(3) < vertices(3,7)) ) then
                isInside = .true.
            end if
            return
         end if
      end if

      if( (Point(1) >= vertices(1,1) .and. Point(1) <= vertices(1,7)) .and. &
          (Point(2) >= vertices(2,1) .and. Point(2) <= vertices(2,7)) .and. &
          (Point(3) >= vertices(3,1) .and. Point(3) <= vertices(3,7)) ) isInside = .true.

   end function isInsideBox

   logical function IsPointInside_eBB( corners, coords )

      implicit none

      real(kind=RP), intent(in) :: corners(NDIM,8) 
      real(kind=RP), intent(in) :: coords(NDIM) 

      real(Kind=RP) :: xMAX, xMIN, yMAX, yMIN, zMAX, zMIN

      IsPointInside_eBB = .false.

      xMAX = maxval(corners(IX,:))
      xMIN = minval(corners(IX,:))
      yMAX = maxval(corners(IY,:))
      yMIN = minval(corners(IY,:))
      zMAX = maxval(corners(IZ,:))
      zMIN = minval(corners(IZ,:))

      if( (coords(IX) .ge. xMIN .and. coords(IX) .le. xMAX) .and. & 
          (coords(IY) .ge. yMIN .and. coords(IY) .le. yMAX) .and. & 
          (coords(IZ) .ge. zMIN .and. coords(IZ) .le. zMAX)       ) then 
         IsPointInside_eBB = .true. 
      end if 

   end function IsPointInside_eBB

end module OrientedBoundingBox