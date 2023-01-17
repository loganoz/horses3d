#include "Includes.h"
module MPI_IBMUtilities

   use SMConstants
   use Utilities
   use TessellationTypes
   use OrientedBoundingBox
   use MPI_Process_info
#ifdef _HAS_MPI_
   use mpi
#endif
   
! maximum 9 points after polygon's been split

!~    ----> axis
!~       |       |
!~    plane1   plane2
   
   implicit none
   
   private
   public :: MPI_KDtree_type, KDtree_partition_type
   public :: MPI_KDtree_buildPartition, MPI_KDtree_destroy, recvSTLPartition, SendSTLPartitions
   public :: MPI_M_Points_type, MPI_Pointpartition
   public :: RootSendPointMaskPartition, MaskCandidates, RootRecvrecvPointMaskPartition
   public :: sendPointMaskPartition, RecvPointMaskPartition
   public :: RootRecvPointMask, RootSendPointMask
   public :: recvPointMask, SendPointMask
   public :: RootrecvBandPoint, RootSendBandPoint, MPI_BandPointpartition
   public :: recvBandPointPartition, sendBandPointPartition
   
   integer, parameter :: ON_PLANE = 0, IN_FRONT_PLANE = 1, BEHIND_PLANE = 2
   
   
   type KDtree_partition_type
   
      real(kind=RP), dimension(2)      :: plane
      real(kind=RP), dimension(NDIM,8) :: vertices
      integer                          :: prev, axis, ID
      type(STLfile)                    :: stl

      contains
         
         procedure :: build       => KDtree_partition_build
         procedure :: SetVertices => KDtree_partition_SetVertices
         procedure :: plot        => KDtree_partition_plot
         procedure :: plotObjs    => KDtree_partition_plotObjs

   end type KDtree_partition_type
   
   type MPI_KDtree_type
   
      type(KDtree_partition_type), dimension(:), allocatable :: partition
      real(kind=RP)                                          :: axisLength
      integer                                                :: axis
   
      contains
      
         procedure :: build           => MPI_KDtree_build
         procedure :: SplittingPlanes => MPI_KDtree_SplittingPlanes
         procedure :: checkObj        => MPI_KDtree_checkObj
         procedure :: STL_partition   => MPI_KDtree_STL_partition
   
   end type MPI_KDtree_type
   
   
   type MPI_M_Points_type
   
      type(point_type), dimension(:), allocatable :: x
      integer,          dimension(:), allocatable :: buffer
      integer                                     :: NumOfObjs,    &
                                                     LocNumOfObjs
   
      contains
      
         procedure :: destroy => MPI_M_Points_type_destroy
   
   end type MPI_M_Points_type  
   
   type(MPI_KDtree_type),       public :: MPI_KDtree_all
   type(KDtree_partition_type), public :: MPI_KDtreePartition
   
   type(MPI_M_Points_type), public :: MPI_M_Points_ALL
   type(MPI_M_Points_type), public :: MPI_M_PointsPartition
   
   type(MPI_M_Points_type), public :: BandPoints_ALL
   
contains

   subroutine MPI_BandPointpartition( BandRegion, NumOfObjs, BandPoints)
   
      implicit none
      !-arguments-----------------------------------------
      type(MPI_M_Points_type), intent(inout) :: BandRegion
      integer,                 intent(in)    :: NumOfObjs
      type(PointLinkedList),   intent(inout) :: BandPoints
      !-local-variables-----------------------------------
      type(point_type), pointer :: p
      integer :: i 

      allocate(BandRegion% buffer(MPI_Process% nProcs))
      allocate(BandRegion% x(NumOfObjs))

      BandRegion% buffer    = 0
      BandRegion% NumOfObjs = NumOfObjs

      BandRegion% buffer(1) = BandPoints% NumOfPoints
      
      if( MPI_Process% isRoot ) then
         p => BandPoints% head
         do i = 1, BandPoints% NumOfPoints
            BandRegion% x(i)% index     = p% index
            BandRegion% x(i)% coords(1) = p% coords(1)
            BandRegion% x(i)% coords(2) = p% coords(2)
            BandRegion% x(i)% coords(3) = p% coords(3)
            BandRegion% x(i)% local_Position(1) = p% local_Position(1)
            BandRegion% x(i)% local_Position(2) = p% local_Position(2)
            BandRegion% x(i)% local_Position(3) = p% local_Position(3)
            BandRegion% x(i)% element_index = p% element_index
            BandRegion% x(i)% partition = p% partition
            p => p% next
         end do
      end if

   end subroutine MPI_BandPointpartition
   
   subroutine MPI_Pointpartition(NumOfObjs)
   
      implicit none
      !-arguments-----------------------
      integer, intent(in) :: NumOfObjs
 
      MPI_M_Points_ALL% NumOfObjs = NumOfObjs
       
      allocate(MPI_M_Points_ALL% buffer(MPI_Process% nProcs))
      allocate(MPI_M_Points_ALL% x(NumOfObjs))

      MPI_M_Points_ALL% buffer    = 0
       
      MPI_M_Points_ALL% buffer(1) = MPI_M_PointsPartition% LocNumOfObjs
      MPI_M_Points_ALL% x(:)% isInsideBody = .false.

      if( MPI_Process% isRoot ) then
         MPI_M_Points_ALL% x(1:MPI_M_PointsPartition% LocNumOfObjs) = &
         MPI_M_PointsPartition% x(1:MPI_M_PointsPartition% LocNumOfObjs) 
      end if

   end subroutine MPI_Pointpartition

   subroutine MPI_KDtree_buildPartition( stl )
   
      implicit none
      !-arguments---------------------------
      type(STLfile), intent(inout) :: stl
      !-local-variables---------------------
      integer :: i

      if( MPI_Process% doMPIRootAction ) then
#ifdef _HAS_MPI_
         call MPI_KDtree_all% build( stl% body )
         call MPI_KDtree_all% SplittingPlanes( stl% body )
         call MPI_KDtree_all% STL_partition( stl% ObjectsList )
#endif
      elseif( MPI_Process% nProcs .eq. 1 ) then 
         MPI_KDtreePartition% stl = stl 
         MPI_KDtreePartition% vertices = OBB(stl% body)% LocVertices
      end if

      if( MPI_Process% doMPIRootAction ) then
         do i = 1, MPI_Process% nProcs
            call MPI_KDtree_ALL% partition(i)% stl% DescribePartitions()
         end do
      end if

   end subroutine MPI_KDtree_buildPartition
   
   subroutine MPI_KDtree_build( this, STLNum )
      use MPI_Process_Info
      implicit none
      !-arguments---------------------------------------
      class(MPI_KDtree_type), intent(inout) :: this
      integer,                intent(in)    :: STLNum
      !-local-variables---------------------------------
      real(kind=RP), dimension(NDIM) :: axis
      integer,       dimension(1)    :: maxvec
      integer                        :: i
     
      allocate(this% partition(MPI_Process% nProcs))
       
      axis = (/ OBB(STLNum)% MBR% Length,OBB(STLNum)% MBR% Width,abs(OBB(STLNum)% nMax) + abs(OBB(STLNum)% nMin) /)
      
      maxvec = maxloc( axis )  
      
      this% axis       = maxvec(1)
      this% axisLength = axis(this% axis)
      
      do i = 1, MPI_Process% nProcs
         this% partition(i)% axis = this% axis
      end do
   
   end subroutine MPI_KDtree_build 
   
   subroutine KDtree_partition_build( this, rank )
   
      implicit none
      !-arguments-------------------------------------------
      class(KDtree_partition_type), intent(inout) :: this
      integer,                      intent(in)    :: rank
   
      this% plane    = 0.0_RP
      this% vertices = 0.0_RP
      this% prev     = 0
      this% axis     = 0
      this% ID       = rank
   
      safedeallocate(this% stl% ObjectsList)
   
   end subroutine KDtree_partition_build
   
   subroutine MPI_KDtree_destroy()
      use MPI_Process_Info
      implicit none
      !-local-varaibles----------
#ifdef _HAS_MPI_
      integer :: i
#endif

      if( MPI_Process% doMPIRootAction ) then
#ifdef _HAS_MPI_
         do i = 1, size(MPI_KDtree_all% partition)
            call MPI_KDtree_all% partition(i)% stl% destroy()
         end do
         deallocate(MPI_KDtree_all% partition)
         MPI_KDtree_all% axis       = 0
         MPI_KDtree_all% axisLength = 0.0_RP
#endif
      end if
      
      call MPI_KDtreePartition% stl% destroy()
      MPI_KDtreePartition% plane    = 0.0_RP
      MPI_KDtreePartition% prev     = 0
      MPI_KDtreePartition% axis     = 0
      MPI_KDtreePartition% ID       = 0
  
   end subroutine MPI_KDtree_destroy

   
  subroutine KDtree_partition_plot( this )
      use MPI_Process_Info
      implicit none
      !-arguments-----------------------------------------
      class(KDtree_partition_type), intent(inout) :: this
      !-local-variables----------------------------------- 
      real(kind=RP), dimension(NDIM) :: x_g
      character(len=LINE_LENGTH)     :: filename, myString
      integer                        :: i, funit
      
      funit = UnusedUnit()
      
      write(myString,'(i100)') MPI_Process% rank
      
      filename = 'STL_partition'//trim(adjustl(myString))
      filename = trim(filename)
      
      open(funit,file='IBM/'//trim(filename)//'.tec', status='unknown')
 
      write(funit,"(a25)") 'TITLE = "Partition Block"'
      write(funit,"(a25)") 'VARIABLES = "x", "y", "z"'
      
      write(funit,"(a69)") 'ZONE NODES=8, ELEMENTS = 6, DATAPACKING=POINT, ZONETYPE=FETETRAHEDRON'
      do i = 1, 8
         call OBB(this% stl% body)% ChangeRefFrame( this% vertices(:,i), 'global', x_g )
         write(funit,'(3E13.5)') x_g(1),x_g(2), x_g(3)
      end do 

      write(funit,'(4i2)') 1, 2, 3, 4
      write(funit,'(4i2)') 1, 5, 8, 4
      write(funit,'(4i2)') 5, 6, 7, 8
      write(funit,'(4i2)') 2, 3, 7, 6
      write(funit,'(4i2)') 4, 8, 7, 3
      write(funit,'(4i2)') 1, 2, 6, 5

      close(funit)
      
   end subroutine KDtree_partition_plot
   
   
   subroutine KDtree_partition_plotObjs( this, STLname )
      use MPI_Process_Info
      implicit none
      !-arguments------------------------------------------ 
      class(KDtree_partition_type), intent(inout) :: this
      character(len=*),             intent(in)    :: STLname
      !-local-variables------------------------------------ 
      real(kind=RP), dimension(NDIM) :: x_g
      character(len=LINE_LENGTH)     :: filename, myString
      integer                        :: i, j, funit
      
      funit = UnusedUnit()
      
      write(myString,'(i100)') MPI_Process% rank
      
      if( MPI_Process% nProcs .gt. 1 ) then
         filename = trim(STLname)//'_partition'//trim(adjustl(myString))
         filename = trim(filename)
      else
         filename = trim(STLname)
      end if
      
      open(funit,file='IBM/'//trim(filename)//'.tec', status='unknown')
 
      write(funit,"(a28)") 'TITLE = "Partition objects"'
      write(funit,"(a25)") 'VARIABLES = "x", "y", "z"'
      
      do i = 1, SIZE(this% stl% ObjectsList)
         write(funit,"(a66)") 'ZONE NODES=3, ELEMENTS = 1, DATAPACKING=POINT, ZONETYPE=FETRIANGLE'
         do j = 1, this% stl% ObjectsList(i)% NumOfVertices
            call OBB(this% stl% body)% ChangeRefFrame( this% stl% ObjectsList(i)% vertices(j)% coords, 'global', x_g )
            write(funit,'(3E13.5)') x_g(1), x_g(2), x_g(3)
         end do

         write(funit,'(3i2)') 1, 2, 3 
      end do

      close(funit)
      
   end subroutine KDtree_partition_plotObjs
   
   subroutine MPI_KDtree_STL_partition( this, objs )
      use MPI_Process_Info
      implicit none
      !-arguments--------------------------------------------------------------------
      class(MPI_KDtree_type),          intent(inout) :: this
      type(Object_type), dimension(:), intent(inout) :: objs
      !-local-variables--------------------------------------------------------------
      type(Object_type)                                 :: objTemp, ObjFront, &
                                                           ObjBack
      type(Object_type), dimension(:), allocatable      :: tria
      integer,           dimension(MPI_Process% nProcs) :: index
      integer,           dimension(2)                   :: PARTITION
      real(kind=RP),     dimension(NDIM)                :: plane_normal, plane_point
      integer                                           :: i, j, k

      index = 0

!$omp parallel shared(this,i,objs,index)
!$omp do schedule(runtime) private(j,PARTITION)
      do i = 1, SIZE(objs)
         
         call this% checkObj( maxval(objs(i)% vertices(:)% coords(this% axis)), &
                              minval(objs(i)% vertices(:)% coords(this% axis)), &
                              PARTITION                                         )
          
         objs(i)% partition(1) = PARTITION(1)
         objs(i)% partition(2) = PARTITION(2)
          
         if( PARTITION(2) - PARTITION(1) .eq. 0 ) then
!$omp critical
            index(PARTITION(1)) = index(PARTITION(1)) + 1
!$omp end critical
         elseif( PARTITION(2) - PARTITION(1) .gt. 0 ) then
!$omp critical
            do j = 0, PARTITION(2) - PARTITION(1)
               index(PARTITION(1)+j) = index(PARTITION(1)+j) + 3 !max number of triagnles is 3
            end do
!$omp end critical
         end if
         
      end do
!$omp end do
!$omp end parallel

      do i = 1, MPI_Process% nProcs
         allocate(this% partition(i)% stl% ObjectsList(index(i)))
      end do      
         
      associate( partitions => this% partition )
           
      index = 0
      plane_normal = 0.0_RP
      plane_normal(this% axis) = 1.0_RP 
      
      do i = 1, SIZE(objs)
         if( objs(i)% partition(1) - objs(i)% partition(2) .eq. 0 ) then
            index(objs(i)% partition(1)) = index(objs(i)% partition(1)) + 1
            partitions(objs(i)% partition(1))% stl% ObjectsList(index(objs(i)% partition(1))) = objs(i)
            partitions(objs(i)% partition(1))% stl% ObjectsList(index(objs(i)% partition(1)))% index = index(objs(i)% partition(1))
         else
            objTemp = objs(i)
            plane_point = 0.0_RP
            do j = 0, objs(i)% partition(2) - objs(i)% partition(1) - 1
               plane_point(this% axis) = partitions(objs(i)% partition(1)+j)% plane(2)
               call ClipPloy( objTemp, plane_normal, plane_point, ObjFront, ObjBack )
               call Poly2Triangles( ObjBack, tria )
               index(objs(i)% partition(1)+j) = index(objs(i)% partition(1)+j) + 1
               partitions(objs(i)% partition(1)+j)% stl%  ObjectsList(index(objs(i)% partition(1)+j)) = ObjBack 
               partitions(objs(i)% partition(1)+j)% stl%  ObjectsList(index(objs(i)% partition(1)+j))% index = index(objs(i)% partition(1)+j) 
               if( allocated(tria) ) then
                  do k = 1, size(tria)
                     index(objs(i)% partition(1)+j) = index(objs(i)% partition(1)+j) + 1
                     partitions(objs(i)% partition(1)+j)% stl% ObjectsList(index(objs(i)% partition(1)+j)) = tria(k)
                     partitions(objs(i)% partition(1)+j)% stl% ObjectsList(index(objs(i)% partition(1)+j))% index = index(objs(i)% partition(1)+j)     
                  end do
                  deallocate(tria)
               end if
               objTemp = ObjFront
            end do
            call Poly2Triangles( ObjFront, tria )
            index(objs(i)% partition(2)) = index(objs(i)% partition(2)) + 1
            partitions(objs(i)% partition(2))% stl% ObjectsList(index(objs(i)% partition(2))) = ObjFront
            partitions(objs(i)% partition(2))% stl% ObjectsList(index(objs(i)% partition(2)))% index = index(objs(i)% partition(2))
            if( allocated(tria) ) then
               do k = 1, size(tria)
                  index(objs(i)% partition(2)) = index(objs(i)% partition(2)) + 1
                  partitions(objs(i)% partition(2))% stl%  ObjectsList(index(objs(i)% partition(2))) = tria(k)
                  partitions(objs(i)% partition(2))% stl%  ObjectsList(index(objs(i)% partition(2)))% index = index(objs(i)% partition(2))     
               end do
               deallocate(tria)
            end if            
         end if
      end do
      
      end associate
      
      do i = 1, MPI_Process% nProcs
         this% partition(i)% stl% NumOfObjs = index(i)
         this% partition(i)% stl% partition = i
      end do 
   
   end subroutine MPI_KDtree_STL_partition
   
   
   subroutine MPI_KDtree_SplittingPlanes( this, STLNum )
      use MPI_Process_Info
      implicit none
      !-arguments--------------------------------------------------
      class(MPI_KDtree_type), intent(inout) :: this
      integer,                intent(in)    :: STLNum
      !-local-variables--------------------------------------------
      real(kind=RP) :: partition, v_max, v_min
      integer       :: i, j
    
      do i = 1, MPI_Process% nProcs
         do j = 1, 8
            this% partition(i)% vertices(:,j) = OBB(STLNum)% LocVertices(:,j)
         end do
      end do

      partition = this% axisLength/MPI_Process% nProcs

      v_max = OBB(STLNum)% LocVertices(this% axis,7)
      v_min = OBB(STLNum)% LocVertices(this% axis,1)

      this% partition(1)% plane(1) = v_min
      this% partition(1)% plane(2) = this% partition(1)% plane(1) + partition 

      call this% partition(1)% SetVertices( 1, this% axis, partition )

      if( MPI_Process% nProcs .eq. 1 ) return

      do i = 2, MPI_Process% nProcs-1
         this% partition(i)% plane(1) = this% partition(i-1)% plane(2) 
         this% partition(i)% plane(2) = this% partition(i)% plane(1) + partition 
         call this% partition(i)% SetVertices( i, this% axis, partition )
      end do
      
      this% partition(MPI_Process% nprocs)% plane(1) = this% partition(MPI_Process% nprocs-1)% plane(2)
      this% partition(MPI_Process% nprocs)% plane(2) = v_max
      call this% partition(MPI_Process% nprocs)% SetVertices( MPI_Process% nprocs, this% axis, partition )
 
   end subroutine MPI_KDtree_SplittingPlanes
   
   
   subroutine KDtree_partition_SetVertices( this, nprocs, axis, Locpartition )
   
      implicit none
      !-arguments-----------------------------------------------------------
      class(KDtree_partition_type),   intent(inout) :: this
      real(kind=RP),                  intent(in)    :: Locpartition
      integer,                        intent(in)    :: axis, nprocs
      !-local-variables-----------------------------------------------------
      integer, dimension(4,2) :: v_indeces
   
      if( axis .eq. 1 ) then
         v_indeces(:,1) = (/1,4,8,5/)
         v_indeces(:,2) = (/2,3,7,6/)
      elseif( axis .eq. 2 ) then
         v_indeces(:,1) = (/1,5,6,2/)
         v_indeces(:,2) = (/4,8,7,3/)
      else
         v_indeces(:,1) = (/1,2,3,4/)
         v_indeces(:,2) = (/5,6,7,8/)      
      end if

      this% vertices(axis,v_indeces(:,2)) = this% vertices(axis,v_indeces(:,1)) + &
                                            (nprocs)*Locpartition
      this% vertices(axis,v_indeces(:,1)) = this% vertices(axis,v_indeces(:,1)) + &
                                            (nprocs-1)*Locpartition
   
   end subroutine KDtree_partition_SetVertices
   
   
   subroutine MPI_KDtree_checkObj( this, coord_max, coord_min, PARTITION )
      use MPI_Process_Info
      implicit none
      !-arguments-----------------------------------------------------------
      class(MPI_KDtree_type),       intent(inout) :: this
      real(kind=RP),                intent(in)    :: coord_max, coord_min
      integer,        dimension(2), intent(out)   :: PARTITION
      !-local-variables------------------------------------------------------
      integer :: i
      
      PARTITION = 0
      
      do i = 1, MPI_Process% nProcs
         if( coord_min .gt. this% partition(i)% plane(1) .and. &
             coord_max .lt. this% partition(i)% plane(2) ) then
             PARTITION(1) = i; PARTITION(2) = i
             exit
         elseif( coord_min .lt. this% partition(i)% plane(1) .and. &
                 coord_max .gt. this% partition(i)% plane(1) ) then
!~       an obj can intersect more the one partition, in this case PARTITION(2)-PARTITION(1) =/ 1 and > 0
            if( PARTITION(1) .eq. 0 ) then
               PARTITION(1) = i-1
            end if
            PARTITION(2) = i
         end if   
      end do
   
   end subroutine MPI_KDtree_checkObj
   
   
   subroutine Poly2Triangles( Poly, tria )
   
      implicit none
      !-arguments-----------------------------------------------------------
      type(object_type),                            intent(inout) :: Poly 
      type(object_type), dimension(:), allocatable, intent(out)   :: tria
      !-local-variables-----------------------------------------------------
      real(kind=RP), dimension(NDIM,Poly% NumOfVertices) :: coords
      integer                                            :: j, n_vert
   
      if( Poly% NumOfVertices .eq. 3 ) return
   
      n_vert = Poly% NumOfVertices
      
      do j = 1, n_vert
         coords(:,j) = Poly% vertices(j)% coords
      end do
   
      if( Poly% NumOfVertices .eq. 4 ) then
         allocate(tria(1))
         tria(1) = Poly
         Poly% NumOfVertices = 3
         tria(1)% NumOfVertices = 3
         deallocate(Poly% vertices, tria(1)% vertices)
         allocate(Poly% vertices(3), tria(1)% vertices(3))
      elseif( Poly% NumOfVertices .eq. 5 ) then
         allocate(tria(2))
         tria(1) = Poly
         tria(2) = Poly
         Poly% NumOfVertices = 3
         tria(1)% NumOfVertices = 3
         tria(2)% NumOfVertices = 3
         deallocate(Poly% vertices,tria(1)% vertices,tria(2)% vertices)
         allocate(Poly% vertices(3),tria(1)% vertices(3),tria(2)% vertices(3))
      else
         print *, 'error in MPI_KDtree, NumOfVertices not recognized', Poly% NumOfVertices
      end if
 
      do j = 1, 3
         Poly% vertices(j)% coords = coords(:,j)
         if( j .eq. 3 ) then
            tria(1)% vertices(j)% coords = coords(:,1)
         else
            tria(1)% vertices(j)% coords = coords(:,j+2)
         end if
         if( size(tria) .eq. 2 ) then
            if( j .eq. 3 ) then
               tria(2)% vertices(j)% coords = coords(:,1)
            else
               tria(2)% vertices(j)% coords = coords(:,j+3)
            end if         
         end if
      end do   
      
   end subroutine Poly2Triangles
   
   
   
   subroutine ClipPloy( obj, plane_normal, plane_point, objFront, objBack )
   
      implicit none
      !-arguments-----------------------------------------------------------------
      type(object_type),           intent(in)    :: obj
      real(kind=rp), dimension(:), intent(in)    :: plane_normal, plane_point
      type(object_type),           intent(out)   :: objFront, objBack
      !-local-variables-----------------------------------------------------------
      real(kind=RP), dimension(NDIM,9) :: PointFront, PointBack 
      real(kind=RP), dimension(NDIM)   :: PointA, PointB, Point_inters
      integer                          :: PointA_Is, PointB_Is, n_front, n_back, i
      
      n_front = 0; n_back = 0
      
      pointA = obj% vertices(obj% NumOfVertices)% coords
      
      PointA_Is = Point_wrt_Plane( plane_normal, plane_point, pointA )
   
      do i = 1, obj% NumOfVertices
         PointB    = obj% vertices(i)% coords
         
         PointB_Is = Point_wrt_Plane( plane_normal, plane_point, pointB )
         if( PointB_Is .eq. IN_FRONT_PLANE ) then
            if( PointA_Is .eq. BEHIND_PLANE ) then
               Point_inters = EdgePlaneIntersection( plane_normal, plane_point, PointA, PointB )
               n_front = n_front + 1
               n_back  = n_back + 1
               PointFront(:,n_front) = Point_Inters
               PointBack(:,n_back) = Point_Inters
            end if
            n_front = n_front + 1
            PointFront(:,n_front) = PointB
         elseif( PointB_Is .eq. BEHIND_PLANE ) then
            if( PointA_Is .eq. IN_FRONT_PLANE ) then
               Point_inters = EdgePlaneIntersection( plane_normal, plane_point, PointA, PointB )
               n_front = n_front + 1
               n_back  = n_back + 1
               PointFront(:,n_front) = Point_Inters
               PointBack(:,n_back) = Point_Inters
            elseif( PointA_Is .eq. ON_PLANE ) then
               n_back  = n_back + 1
               PointBack(:,n_back) = PointA
            end if
            n_back  = n_back + 1
            PointBack(:,n_back) = PointB 
         else
            n_front = n_front + 1
            PointFront(:,n_front) = PointB
            if( PointA_Is .eq. BEHIND_PLANE ) then
               n_back  = n_back + 1
               PointBack(:,n_back) = PointB 
            end if
         end if
         PointA = PointB
         PointA_Is = PointB_Is
      end do
      
      call objFront% build( PointFront, obj% normal, n_front, obj% index, obj% computeIntegrals )
      call objBack% build( PointBack, obj% normal, n_back, obj% index, obj% computeIntegrals )
   
   end subroutine ClipPloy
   
   
   integer function Point_wrt_Plane( plane_normal, plane_point, point ) result( PointIs )
      use MappedGeometryClass
      implicit none
      !-arguments-----------------------------------------------------------------
      real(kind=RP), dimension(:), intent(in) :: plane_normal, plane_point, point
      !-local-variables-----------------------------------------------------------
      real(kind=RP) :: d
   
      d = vdot(plane_normal,point) - vdot(plane_normal,plane_point)
      
      if( AlmostEqual(d,0.0_RP) ) then
         PointIs = ON_PLANE
      elseif( d > 0.0d0 ) then
         PointIs = IN_FRONT_PLANE
      elseif( d < 0.0d0 ) then
         PointIs = BEHIND_PLANE
      end if
   
   end function Point_wrt_Plane
  
   function EdgePlaneIntersection( plane_normal, plane_point, PointA, PointB ) result( Point_inters )
      use MappedGeometryClass
      implicit none
      !-arguments-----------------------------------------------------------------
      real(kind=RP), dimension(:), intent(in) :: plane_normal, plane_point, &
                                                 PointA, PointB
      real(kind=RP), dimension(NDIM)          :: Point_inters
      !-local-variables-----------------------------------------------------------
      real(kind=RP), dimension(NDIM) :: n_AB
      real(kind=RP)                  :: d, t
 
      n_AB = PointB - PointA
      
      d = vdot(plane_normal,plane_point)
      
      t = ( d - vdot(plane_normal,PointA) )/vdot(plane_normal,n_AB)
      
      if( almostEqual(t,0.0_RP) ) then
         Point_inters = PointA
      elseif( almostEqual(t,1.0_RP) ) then
         Point_inters = PointB
      elseif( t .gt. 0.0_RP .and. t .lt. 1.0_RP ) then  
         Point_inters = PointA + t*n_AB
      end if
 
   end function EdgePlaneIntersection
   
   
   subroutine recvSTLPartition()
   
      implicit none
      
#ifdef _HAS_MPI_   
      !-local-variables-------------------------------------------------------------------------------------------   
      integer                                  :: ObjsSize, ierr, j, array_of_statuses(MPI_STATUS_SIZE,19), &
                                                  i, recv_req(19)
      real(kind=RP), dimension(:), allocatable :: normal_x, normal_y, normal_z, COORD
      real(kind=RP), dimension(:), allocatable :: COORD_x1, COORD_x2, COORD_x3
      real(kind=RP), dimension(:), allocatable :: COORD_y1, COORD_y2, COORD_y3
      real(kind=RP), dimension(:), allocatable :: COORD_z1, COORD_z2, COORD_z3
      real(kind=RP), dimension(:), allocatable :: vertex1, vertex2, vertex3, vertex4, &
                                                  vertex5, vertex6, vertex7, vertex8
      integer,       dimension(:), allocatable :: indeces
      
      if( MPI_Process% isRoot ) return 
      
      ! receive n of objects
      call mpi_irecv( ObjsSize, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(1), ierr ) 
      
      call mpi_wait(recv_req(1), MPI_STATUS_IGNORE, ierr)
      
      MPI_KDtreePartition% stl% NumOfObjs = ObjsSize
      
      allocate(MPI_KDtreePartition% stl% ObjectsList(ObjsSize))
      
      call mpi_irecv( MPI_KDtreePartition% plane, 2, MPI_DOUBLE_PRECISION, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(2), ierr  )
   
      call mpi_irecv( MPI_KDtreePartition% axis, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(3), ierr  )

      allocate( indeces(ObjsSize),  &
                normal_x(ObjsSize), &
                normal_y(ObjsSize), & 
                normal_z(ObjsSize), &
                COORD_x1(ObjsSize), &
                COORD_x2(ObjsSize), &
                COORD_x3(ObjsSize), &
                COORD_y1(ObjsSize), &
                COORD_y2(ObjsSize), &
                COORD_y3(ObjsSize), &
                COORD_z1(ObjsSize), &
                COORD_z2(ObjsSize), &
                COORD_z3(ObjsSize), &
                vertex1(8),         &
                vertex2(8),         &
                vertex3(8)          )

      call mpi_irecv( indeces, ObjsSize, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(4), ierr )
                
      call mpi_irecv( normal_x, ObjsSize, MPI_DOUBLE_PRECISION, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(5), ierr )  
                      
      call mpi_irecv( normal_y, ObjsSize, MPI_DOUBLE_PRECISION, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(6), ierr )  
                      
      call mpi_irecv( normal_z, ObjsSize, MPI_DOUBLE_PRECISION, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(7), ierr )  
                      
      call mpi_irecv( COORD_x1, ObjsSize, MPI_DOUBLE_PRECISION, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(8), ierr )  
                      
      call mpi_irecv( COORD_x2, ObjsSize, MPI_DOUBLE_PRECISION, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(9), ierr )  
                      
      call mpi_irecv( COORD_x3, ObjsSize, MPI_DOUBLE_PRECISION, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(10), ierr )  
                      
      call mpi_irecv( COORD_y1, ObjsSize, MPI_DOUBLE_PRECISION, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(11), ierr )  
                      
      call mpi_irecv( COORD_y2, ObjsSize, MPI_DOUBLE_PRECISION, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(12), ierr )  
                      
      call mpi_irecv( COORD_y3, ObjsSize, MPI_DOUBLE_PRECISION, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(13), ierr )  
                      
      call mpi_irecv( COORD_z1, ObjsSize, MPI_DOUBLE_PRECISION, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(14), ierr )  
                      
      call mpi_irecv( COORD_z2, ObjsSize, MPI_DOUBLE_PRECISION, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(15), ierr )  
                      
      call mpi_irecv( COORD_z3, ObjsSize, MPI_DOUBLE_PRECISION, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(16), ierr )  
                      
      call mpi_irecv( vertex1, 8, MPI_DOUBLE_PRECISION, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(17), ierr )    
                      
      call mpi_irecv( vertex2, 8, MPI_DOUBLE_PRECISION, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(18), ierr )  
                      
      call mpi_irecv( vertex3, 8, MPI_DOUBLE_PRECISION, 0, MPI_ANY_TAG, MPI_COMM_WORLD, recv_req(19), ierr )    

      call mpi_waitall(19, recv_req, array_of_statuses, ierr)     

      do i = 1, MPI_KDtreePartition% stl% NumOfObjs
         MPI_KDtreePartition% stl% ObjectsList(i)% NumOfVertices = 3
         allocate( MPI_KDtreePartition% stl% ObjectsList(i)% vertices(3) )
         MPI_KDtreePartition% stl% ObjectsList(i)% index     = i
         MPI_KDtreePartition% stl% ObjectsList(i)% normal(1) = normal_x(i)
         MPI_KDtreePartition% stl% ObjectsList(i)% normal(2) = normal_y(i)
         MPI_KDtreePartition% stl% ObjectsList(i)% normal(3) = normal_z(i)
         MPI_KDtreePartition% stl% ObjectsList(i)% vertices(1)% coords(1) = COORD_x1(i)
         MPI_KDtreePartition% stl% ObjectsList(i)% vertices(2)% coords(1) = COORD_x2(i)
         MPI_KDtreePartition% stl% ObjectsList(i)% vertices(3)% coords(1) = COORD_x3(i)
         MPI_KDtreePartition% stl% ObjectsList(i)% vertices(1)% coords(2) = COORD_y1(i)
         MPI_KDtreePartition% stl% ObjectsList(i)% vertices(2)% coords(2) = COORD_y2(i)
         MPI_KDtreePartition% stl% ObjectsList(i)% vertices(3)% coords(2) = COORD_y3(i)
         MPI_KDtreePartition% stl% ObjectsList(i)% vertices(1)% coords(3) = COORD_z1(i)
         MPI_KDtreePartition% stl% ObjectsList(i)% vertices(2)% coords(3) = COORD_z2(i)
         MPI_KDtreePartition% stl% ObjectsList(i)% vertices(3)% coords(3) = COORD_z3(i)
      end do
      
      MPI_KDtreePartition% vertices(1,:) = vertex1
      MPI_KDtreePartition% vertices(2,:) = vertex2
      MPI_KDtreePartition% vertices(3,:) = vertex3
    
      deallocate( indeces,  &
                  normal_x, &
                  normal_y, &
                  normal_z, &
                  COORD_x1, &
                  COORD_x2, &
                  COORD_x3, & 
                  COORD_y1, &
                  COORD_y2, &
                  COORD_y3, &
                  COORD_z1, &
                  COORD_z2, &
                  COORD_z3, &
                  vertex1,  &
                  vertex2,  &
                  vertex3   )
#endif

   end subroutine recvSTLPartition
   
   subroutine SendSTLPartitions()
   
      implicit none
   
#ifdef _HAS_MPI_  
      !-local-variables-------------------------------------------------------------------------
      integer                                    :: ObjsSize, PartSize, nProcs, i, j, ierr,    &
                                                    msg, array_of_statuses(MPI_STATUS_SIZE,19)
      real(kind=RP), dimension(:),   allocatable :: normal_x, normal_y, normal_z
      real(kind=RP), dimension(:),   allocatable :: vertex1, vertex2, vertex3
      real(kind=RP), dimension(:),   allocatable :: COORD_x1, COORD_x2, COORD_x3
      real(kind=RP), dimension(:),   allocatable :: COORD_y1, COORD_y2, COORD_y3
      real(kind=RP), dimension(:),   allocatable :: COORD_z1, COORD_z2, COORD_z3
      integer,       dimension(:),   allocatable :: indeces
      integer,       dimension(:,:), allocatable :: send_req

      allocate(send_req(MPI_Process% nProcs-1,19))

      do nProcs = 2, MPI_Process% nProcs
      
         ObjsSize = MPI_KDtree_all% partition(nProcs)% stl% NumOfObjs
         call mpi_isend( ObjsSize, 1, MPI_INT, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, &
                         send_req(nProcs-1,1), ierr                                   )
         
      end do
      
      do nProcs = 2, MPI_Process% nProcs
         
         partSize = MPI_KDtree_all% partition(nProcs)% stl% NumOfObjs
         
         call mpi_isend( MPI_KDtree_all% partition(nProcs)% plane, 2, MPI_DOUBLE_PRECISION, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, &
                         send_req(nProcs-1,2), ierr)

         call mpi_isend( MPI_KDtree_all% partition(nProcs)% axis, 1, MPI_INT, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, &
                         send_req(nProcs-1,3), ierr)
           
           
         allocate( indeces(partSize),  &
                   normal_x(partSize), &
                   normal_y(partSize), & 
                   normal_z(partSize), &
                   COORD_x1(partSize), &
                   COORD_x2(partSize), &
                   COORD_x3(partSize), &
                   COORD_y1(partSize), &
                   COORD_y2(partSize), &
                   COORD_y3(partSize), &
                   COORD_z1(partSize), &
                   COORD_z2(partSize), &
                   COORD_z3(partSize), &
                   vertex1(8),         &
                   vertex2(8),         &                                           
                   vertex3(8)          )

      do i = 1, partSize
         indeces(i)  = MPI_KDtree_all% partition(nProcs)% stl% ObjectsList(i)% index
         normal_x(i) = MPI_KDtree_all% partition(nProcs)% stl% ObjectsList(i)% normal(1)
         normal_y(i) = MPI_KDtree_all% partition(nProcs)% stl% ObjectsList(i)% normal(2)
         normal_z(i) = MPI_KDtree_all% partition(nProcs)% stl% ObjectsList(i)% normal(3)
         COORD_x1(i) = MPI_KDtree_all% partition(nProcs)% stl% ObjectsList(i)% vertices(1)% coords(1)
         COORD_x2(i) = MPI_KDtree_all% partition(nProcs)% stl% ObjectsList(i)% vertices(2)% coords(1)
         COORD_x3(i) = MPI_KDtree_all% partition(nProcs)% stl% ObjectsList(i)% vertices(3)% coords(1)
         COORD_y1(i) = MPI_KDtree_all% partition(nProcs)% stl% ObjectsList(i)% vertices(1)% coords(2)
         COORD_y2(i) = MPI_KDtree_all% partition(nProcs)% stl% ObjectsList(i)% vertices(2)% coords(2)
         COORD_y3(i) = MPI_KDtree_all% partition(nProcs)% stl% ObjectsList(i)% vertices(3)% coords(2)
         COORD_z1(i) = MPI_KDtree_all% partition(nProcs)% stl% ObjectsList(i)% vertices(1)% coords(3)
         COORD_z2(i) = MPI_KDtree_all% partition(nProcs)% stl% ObjectsList(i)% vertices(2)% coords(3)
         COORD_z3(i) = MPI_KDtree_all% partition(nProcs)% stl% ObjectsList(i)% vertices(3)% coords(3)
      end do

         vertex1 = MPI_KDtree_all% partition(nProcs)% vertices(1,:)
         vertex2 = MPI_KDtree_all% partition(nProcs)% vertices(2,:)
         vertex3 = MPI_KDtree_all% partition(nProcs)% vertices(3,:)

         call mpi_isend( indeces, partSize, MPI_INT, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,4), ierr )

         call mpi_isend( normal_x, partSize, MPI_DOUBLE_PRECISION, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,5), ierr  )         

         call mpi_isend( normal_y, partSize, MPI_DOUBLE_PRECISION, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,6), ierr  )         

         call mpi_isend( normal_z, partSize, MPI_DOUBLE_PRECISION, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,7), ierr  )         

         call mpi_isend( COORD_x1, partSize, MPI_DOUBLE_PRECISION, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,8), ierr  )            

         call mpi_isend( COORD_x2, partSize, MPI_DOUBLE_PRECISION, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,9), ierr  )            

         call mpi_isend( COORD_x3, partSize, MPI_DOUBLE_PRECISION, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,10), ierr )            

         call mpi_isend( COORD_y1, partSize, MPI_DOUBLE_PRECISION, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,11), ierr )            

         call mpi_isend( COORD_y2, partSize, MPI_DOUBLE_PRECISION, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,12), ierr )            

         call mpi_isend( COORD_y3, partSize, MPI_DOUBLE_PRECISION, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,13), ierr )            

         call mpi_isend( COORD_z1, partSize, MPI_DOUBLE_PRECISION, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,14), ierr )            

         call mpi_isend( COORD_z2, partSize, MPI_DOUBLE_PRECISION, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,15), ierr )            

         call mpi_isend( COORD_z3, partSize, MPI_DOUBLE_PRECISION, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,16), ierr )            

         call mpi_isend( vertex1, 8, MPI_DOUBLE_PRECISION, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,17), ierr )

         call mpi_isend( vertex2, 8, MPI_DOUBLE_PRECISION, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,18), ierr )

         call mpi_isend( vertex3, 8, MPI_DOUBLE_PRECISION, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, send_req(nProcs-1,19), ierr )

         call mpi_waitall(19, send_req(nProcs-1,:), array_of_statuses, ierr)
 
         deallocate( indeces,  &
                     normal_x, &
                     normal_y, &
                     normal_z, &
                     COORD_x1, &
                     COORD_x2, &
                     COORD_x3, &
                     COORD_y1, &
                     COORD_y2, &
                     COORD_y3, &
                     COORD_z1, &
                     COORD_z2, &
                     COORD_z3, &
                     vertex1,  &
                     vertex2,  &
                     vertex3   )
         
      end do

      MPI_KDtreePartition = MPI_KDtree_all% partition(1)     
   
      deallocate(send_req)
   
#endif   

   end subroutine SendSTLPartitions
   

! PARTITION FOR MASK POINTS

   subroutine MPI_M_Points_type_Destroy( this )
      
      implicit none
      !-arguments-------------------------------------
      class(MPI_M_Points_type), intent(inout) :: this
   
      if( allocated(this% x) )      deallocate(this% x)
      if( allocated(this% buffer) ) deallocate(this% buffer)
      
      this% NumOfObjs     = 0
      this% LocNumOfObjs  = 0
      
   end subroutine MPI_M_Points_type_Destroy

   subroutine MaskCandidates( elements, no_of_elements, no_of_DoFs, STLNum, NumOfSTL ) 
      use ElementClass
      implicit none
      !-arguments----------------------------------------------------------------
      type(element),  dimension(:), intent(inout) :: elements
      integer,                      intent(in)    :: no_of_elements, no_of_DoFs, &
                                                     STLNum, NumOfSTL
      !-local-variables-----------------------------------------------------------
      integer       :: n, eID, i, j, k
      
      allocate(MPI_M_PointsPartition% x(no_of_DoFs))
      
      n = 0

!$omp parallel shared(n,eID,elements,no_of_elements,MPI_M_PointsPartition,STLNum,NumOfSTL)
!$omp do schedule(runtime) private(i,j,k)
      do eID = 1, no_of_elements
         if( .not. allocated(elements(eID)% isInsideBody) ) then
            call elements(eID)% ConstructIBM(elements(eID)% Nxyz(1), elements(eID)% Nxyz(2), elements(eID)% Nxyz(3), NumOfSTL )
         end if
         
         do k = 0, elements(eID)% Nxyz(3); do j = 0, elements(eID)% Nxyz(2) ; do i = 0, elements(eID)% Nxyz(1)

            if( elements(eID)% isInsideBody(i,j,k) ) cycle

            elements(eID)% isInsideBody(i,j,k) = OBB(STLNum)% isPointInside( coords = elements(eID)% geom% x(:,i,j,k) )

            if( elements(eID)% isInsideBody(i,j,k) ) then
!$omp critical
               n = n + 1
               MPI_M_PointsPartition% x(n)% coords = elements(eID)% geom% x(:,i,j,k)
               MPI_M_PointsPartition% x(n)% local_Position(1) = i
               MPI_M_PointsPartition% x(n)% local_Position(2) = j
               MPI_M_PointsPartition% x(n)% local_Position(3) = k
               MPI_M_PointsPartition% x(n)% element_index = eID
               MPI_M_PointsPartition% x(n)% partition = MPI_Process% rank
               MPI_M_PointsPartition% LocNumOfObjs = n
!$omp end critical
            end if
         end do; end do; end do
          
      end do
!$omp end do
!$omp end parallel 

   end subroutine MaskCandidates
   
   
   subroutine RootRecvrecvPointMaskPartition()
   
      implicit none
   
#ifdef _HAS_MPI_     
      !-local-variables---------------------------------------------------------------- 
      integer                                    :: ObjsSize, nProcs, ierr, i, msg,      &
                                                    array_of_statuses(MPI_STATUS_SIZE,9)
      real(kind=RP), dimension(:),   allocatable :: COORD_x, COORD_y, COORD_z
      integer,       dimension(:),   allocatable :: i_v, j_v, k_v, eID, partition, buffer
      integer,       dimension(:,:), allocatable :: RootMaskrecv_req

      if( MPI_Process% isRoot ) then

         allocate( buffer(MPI_Process% nProcs-1),            &
                   RootMaskrecv_req(MPI_Process% nProcs-1,9) )

         do nProcs = 2, MPI_Process% nProcs
            call mpi_irecv( ObjsSize, 1, MPI_INT, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, RootMaskrecv_req(nProcs-1,1), ierr ) 
            call mpi_wait(RootMaskrecv_req(nProcs-1,1), MPI_STATUS_IGNORE, ierr)
            
            MPI_M_Points_All% buffer(nProcs) = MPI_M_Points_All% buffer(nProcs-1) + ObjsSize
            buffer(nProcs-1) = ObjsSize
            
         end do         
            
         do nProcs = 2, MPI_Process% nProcs
      
             allocate( COORD_x(buffer(nProcs-1)),  &
                       COORD_y(buffer(nProcs-1)),  &
                       COORD_z(buffer(nProcs-1)),  &
                       i_v(buffer(nProcs-1)),      &
                       j_v(buffer(nProcs-1)),      &
                       k_v(buffer(nProcs-1)),      &
                       eID(buffer(nProcs-1)),      &
                       partition(buffer(nProcs-1)) )
                                  
             call mpi_irecv( COORD_x, buffer(nProcs-1), MPI_DOUBLE, nProcs-1, &
                             MPI_ANY_TAG, MPI_COMM_WORLD, RootMaskrecv_req(nProcs-1,2), ierr            )  
          
             call mpi_irecv( COORD_y, buffer(nProcs-1), MPI_DOUBLE, nProcs-1, &
                             MPI_ANY_TAG, MPI_COMM_WORLD, RootMaskrecv_req(nProcs-1,3), ierr             )  

             call mpi_irecv( COORD_z, buffer(nProcs-1), MPI_DOUBLE, nProcs-1, &
                             MPI_ANY_TAG, MPI_COMM_WORLD, RootMaskrecv_req(nProcs-1,4), ierr             )   
      
             call mpi_irecv( i_v, buffer(nProcs-1), MPI_INT, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, RootMaskrecv_req(nProcs-1,5), ierr )  
      
             call mpi_irecv( j_v, buffer(nProcs-1), MPI_INT, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, RootMaskrecv_req(nProcs-1,6), ierr )  
      
             call mpi_irecv( k_v, buffer(nProcs-1), MPI_INT, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, RootMaskrecv_req(nProcs-1,7), ierr )  
      
             call mpi_irecv( eID, buffer(nProcs-1), MPI_INT, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, RootMaskrecv_req(nProcs-1,8), ierr )         
      
             call mpi_irecv( partition, buffer(nProcs-1), MPI_INT, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, RootMaskrecv_req(nProcs-1,9), ierr )         

             call mpi_waitall(9, RootMaskrecv_req(nProcs-1,:), array_of_statuses, ierr)                        
             
             do i = 1, buffer(nProcs-1)             
                MPI_M_Points_All% x(i+MPI_M_Points_All% buffer(nProcs-1))% coords(1) = COORD_x(i)
                MPI_M_Points_All% x(i+MPI_M_Points_All% buffer(nProcs-1))% coords(2) = COORD_y(i)
                MPI_M_Points_All% x(i+MPI_M_Points_All% buffer(nProcs-1))% coords(3) = COORD_z(i)
                MPI_M_Points_All% x(i+MPI_M_Points_All% buffer(nProcs-1))% local_Position(1) = i_v(i)
                MPI_M_Points_All% x(i+MPI_M_Points_All% buffer(nProcs-1))% local_Position(2) = j_v(i)
                MPI_M_Points_All% x(i+MPI_M_Points_All% buffer(nProcs-1))% local_Position(3) = k_v(i)
                MPI_M_Points_All% x(i+MPI_M_Points_All% buffer(nProcs-1))% element_index = eID(i)  
                MPI_M_Points_All% x(i+MPI_M_Points_All% buffer(nProcs-1))% partition = partition(i)  
             end do
                      
             deallocate( COORD_x, COORD_y, COORD_z, i_v, j_v, k_v, eID, partition )

          end do
          
          deallocate(buffer,RootMaskrecv_req)

       end if
#endif          
   
   end subroutine RootRecvrecvPointMaskPartition
   
   subroutine RootSendPointMaskPartition()
   
      implicit none
   
#ifdef _HAS_MPI_  
      !-local-variables----------------------------------------------------------------
      integer                                  :: ObjsSize, i, j, ierr, msg, &
                                                  array_of_statuses(MPI_STATUS_SIZE,9), &
                                                  RootMasksend_req(9)
      real(kind=RP), dimension(:), allocatable :: COORD_x, COORD_y, COORD_z
      integer,       dimension(:), allocatable :: i_v, j_v, k_v, eID, partition
      
      if( MPI_Process% isRoot ) return
      
      ObjsSize = MPI_M_PointsPartition% LocNumOfObjs 
      
      call mpi_isend( ObjsSize, 1, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD, &
                      RootMasksend_req(1), ierr                                 )       
       
      allocate( COORD_x(ObjsSize),  &
                COORD_y(ObjsSize),  &
                COORD_z(ObjsSize),  &
                i_v(ObjsSize),      &
                j_v(ObjsSize),      &
                k_v(ObjsSize),      &
                eID(ObjsSize),      &
                partition(ObjsSize) )

      COORD_x = MPI_M_PointsPartition% x(1:ObjsSize)% coords(1)
      COORD_y = MPI_M_PointsPartition% x(1:ObjsSize)% coords(2)
      COORD_z = MPI_M_PointsPartition% x(1:ObjsSize)% coords(3)
      
      i_v = MPI_M_PointsPartition% x(1:ObjsSize)% local_Position(1)
      j_v = MPI_M_PointsPartition% x(1:ObjsSize)% local_Position(2)
      k_v = MPI_M_PointsPartition% x(1:ObjsSize)% local_Position(3)
      eID = MPI_M_PointsPartition% x(1:ObjsSize)% element_index
      partition = MPI_M_PointsPartition% x(1:ObjsSize)% partition

      call mpi_isend( COORD_x, MPI_M_PointsPartition% LocNumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, &
                      RootMasksend_req(2), ierr                                                                           )
                         
      call mpi_isend( COORD_y, MPI_M_PointsPartition% LocNumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, &
                      RootMasksend_req(3), ierr                                                                           )
                         
      call mpi_isend( COORD_z, MPI_M_PointsPartition% LocNumOfObjs, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, &
                      RootMasksend_req(4), ierr                                                                           )
                         
      call mpi_isend( i_v, MPI_M_PointsPartition% LocNumOfObjs, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD, &
                      RootMasksend_req(5), ierr                                                          )
                         
      call mpi_isend( j_v, MPI_M_PointsPartition% LocNumOfObjs, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD, &
                      RootMasksend_req(6), ierr                                                          )
                         
      call mpi_isend( k_v, MPI_M_PointsPartition% LocNumOfObjs, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD, &
                      RootMasksend_req(7), ierr                                                          )
                         
      call mpi_isend( eID, MPI_M_PointsPartition% LocNumOfObjs, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD, &
                      RootMasksend_req(8), ierr                                                          )
                         
      call mpi_isend( partition, MPI_M_PointsPartition% LocNumOfObjs, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD, &
                      RootMasksend_req(9), ierr                                                                )
           
      call mpi_waitall(9, RootMasksend_req, array_of_statuses, ierr)        

      deallocate( COORD_x, COORD_y, COORD_z, i_v, j_v, k_v, eID, partition )  
#endif  

   end subroutine RootSendPointMaskPartition
 
 
   subroutine recvPointMaskPartition()
    
      implicit none
      
#ifdef _HAS_MPI_    
      !-local-variables-----------------------------------------------------------------  
      integer                                    :: ObjsSize, ierr, i, &
                                                    array_of_statuses(MPI_STATUS_SIZE,8), &
                                                    Maskrecv_req(8)
      real(kind=RP), dimension(:),   allocatable :: COORD_x, COORD_y, COORD_z
      integer,       dimension(:),   allocatable :: i_v, j_v, k_v, eID, partition
      
      if( MPI_Process% isRoot ) return 

      ObjsSize = MPI_M_Points_ALL% NumOfObjs

      allocate( COORD_x(ObjsSize),  &
                COORD_y(ObjsSize),  &
                COORD_z(ObjsSize),  &
                i_v(ObjsSize),      &
                j_v(ObjsSize),      &
                k_v(ObjsSize),      &
                eID(ObjsSize),      &
                partition(ObjsSize) )      
                            
      call mpi_irecv( COORD_x, ObjsSize, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, Maskrecv_req(1), ierr )  

      call mpi_irecv( COORD_y, ObjsSize, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, Maskrecv_req(2), ierr )  

      call mpi_irecv( COORD_z, ObjsSize, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, Maskrecv_req(3), ierr )   

      call mpi_irecv( i_v, ObjsSize, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, Maskrecv_req(4), ierr )  

      call mpi_irecv( j_v, ObjsSize, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, Maskrecv_req(5), ierr )  

      call mpi_irecv( k_v, ObjsSize, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, Maskrecv_req(6), ierr )  

      call mpi_irecv( eID, ObjsSize, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, Maskrecv_req(7), ierr ) 

      call mpi_irecv( partition, ObjsSize, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, Maskrecv_req(8), ierr ) 

      call mpi_waitall(8, Maskrecv_req, array_of_statuses, ierr)     

      do i = 1, ObjsSize
         MPI_M_Points_All% x(i)% coords(1) = COORD_x(i)
         MPI_M_Points_All% x(i)% coords(2) = COORD_y(i)
         MPI_M_Points_All% x(i)% coords(3) = COORD_z(i)
         MPI_M_Points_All% x(i)% local_Position(1) = i_v(i)
         MPI_M_Points_All% x(i)% local_Position(2) = j_v(i)
         MPI_M_Points_All% x(i)% local_Position(3) = k_v(i)
         MPI_M_Points_All% x(i)% element_index = eID(i)  
         MPI_M_Points_All% x(i)% partition = partition(i)  
      end do 
 
      deallocate( COORD_x, COORD_y, COORD_z, i_v, j_v, k_v, eID, partition )
#endif
   
   end subroutine recvPointMaskPartition

   subroutine sendPointMaskPartition()
   
      implicit none
   
#ifdef _HAS_MPI_  
      !-local-variables----------------------------------------------------------------
      integer                                    :: ObjsSize, nProcs, msg, ierr, &
                                                    array_of_statuses(MPI_STATUS_SIZE,9)
      real(kind=RP), dimension(:),   allocatable :: COORD_x, COORD_y, COORD_z
      integer,       dimension(:),   allocatable :: i_v, j_v, k_v, eID, partition
      integer,       dimension(:,:), allocatable :: Masksend_req
   
      ObjsSize = MPI_M_Points_ALL% NumOfObjs
      
      allocate( COORD_x(ObjsSize),                    &
                COORD_y(ObjsSize),                    & 
                COORD_z(ObjsSize),                    &
                i_v(ObjsSize),                        &
                j_v(ObjsSize),                        &
                k_v(ObjsSize),                        &
                eID(ObjsSize),                        &
                partition(ObjsSize),                  &
                Masksend_req(MPI_Process% nProcs-1,8) )
      
      COORD_x = MPI_M_Points_ALL% x(:)% coords(1)
      COORD_y = MPI_M_Points_ALL% x(:)% coords(2)
      COORD_z = MPI_M_Points_ALL% x(:)% coords(3)
      i_v = MPI_M_Points_ALL% x(:)% local_Position(1)
      j_v = MPI_M_Points_ALL% x(:)% local_Position(2)
      k_v = MPI_M_Points_ALL% x(:)% local_Position(3)
      eID = MPI_M_Points_ALL% x(:)% element_index 
      partition = MPI_M_Points_ALL% x(:)% partition
   
      do nProcs = 2, MPI_Process% nProcs        
      
         call mpi_isend( COORD_x, ObjsSize, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, &
                         Masksend_req(nProcs-1,1), ierr                                        )
      
         call mpi_isend( COORD_y, ObjsSize, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, &
                         Masksend_req(nProcs-1,2), ierr                                        )
      
         call mpi_isend( COORD_z, ObjsSize, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, &
                         Masksend_req(nProcs-1,3), ierr                                        )
      
         call mpi_isend( i_v, ObjsSize, MPI_INT, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, &
                         Masksend_req(nProcs-1,4), ierr                                 )
      
         call mpi_isend( j_v, ObjsSize, MPI_INT, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, &
                         Masksend_req(nProcs-1,5), ierr                                 )
      
         call mpi_isend( k_v, ObjsSize, MPI_INT, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, &
                         Masksend_req(nProcs-1,6), ierr                                 )
      
         call mpi_isend( eID, ObjsSize, MPI_INT, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, &
                         Masksend_req(nProcs-1,7), ierr                                  )
      
         call mpi_isend( partition, ObjsSize, MPI_INT, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, &
                         Masksend_req(nProcs-1,8), ierr                                        )
                         
         call mpi_waitall(8, Masksend_req(nProcs-1,:), array_of_statuses, ierr  )                
         
      end do
                         
      deallocate(COORD_x,COORD_y,COORD_z,i_v,j_v,k_v,eID,partition,Masksend_req)
#endif   
   
   end subroutine sendPointMaskPartition
   
      
   subroutine RootRecvPointMask()
   
      implicit none
   
#ifdef _HAS_MPI_      
      !-local-variables---------------------------------------------------------------
      integer                            :: ObjsSize, nProcs, ierr, i, msg,      &
                                            status(MPI_STATUS_SIZE), RootMaskrecv_req
      logical, dimension(:), allocatable :: isInsideBody
      
      if( MPI_Process% isRoot ) then      
      
         allocate( isInsideBody(MPI_M_Points_All% NumOfObjs) )
            
         do nProcs = 2, MPI_Process% nProcs
      

             call mpi_irecv( isInsideBody, MPI_M_Points_All% NumOfObjs, MPI_LOGICAL, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, RootMaskrecv_req, ierr )

             call mpi_wait(RootMaskrecv_req, status, ierr)
                        
             do i = 1, MPI_M_Points_All% NumOfObjs
                if( MPI_M_Points_All% x(i)% isInsideBody ) cycle            
                MPI_M_Points_All% x(i)% isInsideBody = isInsideBody(i)  
             end do
                      

          end do

          deallocate( isInsideBody )

       end if
#endif          
   
   end subroutine RootRecvPointMask
   
   subroutine RootSendPointMask()
   
      implicit none
   
#ifdef _HAS_MPI_  
      !-local-variables-----------------------------------------------------------
      integer                            :: ObjsSize, i, j, ierr, msg,          &
                                            status(MPI_STATUS_SIZE), RootMasksend_req
      logical, dimension(:), allocatable :: isInsideBody
      
      if( MPI_Process% isRoot ) return
      
      ObjsSize = MPI_M_Points_ALL% NumOfObjs       
 
      allocate( isInsideBody(ObjsSize) )

      isInsideBody = MPI_M_Points_ALL% x(:)% isInsideBody

      call mpi_isend( isInsideBody, MPI_M_Points_ALL% NumOfObjs, MPI_LOGICAL, 0, DEFAULT_TAG, MPI_COMM_WORLD, RootMasksend_req, ierr )
      
      call mpi_wait(RootMasksend_req, status, ierr)
  
      deallocate( isInsideBody )
#endif  

   end subroutine RootSendPointMask
   
 
   subroutine recvPointMask()
    
      implicit none
      
#ifdef _HAS_MPI_    
      !-local-variables-----------------------------------------------------------  
      integer                            :: ObjsSize, ierr, i, &
                                            status(MPI_STATUS_SIZE), Maskrecv_req
      logical, dimension(:), allocatable :: isInsideBody
      
      if( MPI_Process% isRoot ) return 

      ObjsSize = MPI_M_Points_ALL% NumOfObjs

      allocate( isInsideBody(ObjsSize) )      

      call mpi_irecv( isInsideBody, ObjsSize, MPI_LOGICAL, 0, MPI_ANY_TAG, MPI_COMM_WORLD, Maskrecv_req, ierr )

      call mpi_wait(Maskrecv_req, status, ierr)

      do i = 1, ObjsSize 
         MPI_M_Points_All% x(i)% isInsideBody = isInsideBody(i)  
      end do 
 
      deallocate( isInsideBody )
#endif
   
   end subroutine recvPointMask

   subroutine sendPointMask()
   
      implicit none
   
#ifdef _HAS_MPI_  
      !-local-variables-----------------------------------------------------------------------------
      integer                              :: ObjsSize, nProcs, msg, ierr, &
                                              array_of_statuses(MPI_STATUS_SIZE,MPI_Process% nProcs-1)
      logical, dimension(:), allocatable :: isInsideBody
      integer, dimension(:), allocatable :: Masksend_req
   
      ObjsSize = MPI_M_Points_ALL% NumOfObjs
      
      allocate( isInsideBody(ObjsSize), Masksend_req(MPI_Process% nProcs-1) )
        
      isInsideBody = MPI_M_Points_ALL% x(:)% isInsideBody
        
      do nProcs = 2, MPI_Process% nProcs          
      
         call mpi_isend( isInsideBody, ObjsSize, MPI_LOGICAL, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, &
                         Masksend_req(nProcs-1), ierr                                          )
         
      end do

      call mpi_waitall(MPI_Process% nProcs-1, Masksend_req, array_of_statuses, ierr  )
      
      deallocate(isInsideBody, Masksend_req)
#endif   
   
   end subroutine sendPointMask  
   
   
   
   
! BAND REGION 
   
    subroutine RootrecvBandPoint( BandRegion )
   
      implicit none
      !-arguments-----------------------------------------------------------------------
      type(MPI_M_Points_type), intent(inout) :: BandRegion
#ifdef _HAS_MPI_      
      !-local-variables-----------------------------------------------------------------
      integer                                    :: ObjsSize, nProcs, ierr, i, msg, &
                                                    array_of_statuses(MPI_STATUS_SIZE,9)
      real(kind=RP), dimension(:),   allocatable :: COORD_x, COORD_y, COORD_z
      integer,       dimension(:),   allocatable :: buffer, i_v, j_v, k_v, eID, partition
      integer,       dimension(:,:), allocatable :: RootBandrecv_req
    
      if( MPI_Process% isRoot ) then

         allocate(buffer(MPI_Process% nProcs-1),RootBandrecv_req(MPI_Process% nProcs-1,9))

         do nProcs = 2, MPI_Process% nProcs
         
            call mpi_irecv( ObjsSize, 1, MPI_INT, nProcs-1, MPI_ANY_TAG, MPI_COMM_WORLD, RootBandrecv_req(nProcs-1,1), ierr ) 
            call mpi_wait(RootBandrecv_req(nProcs-1,1), MPI_STATUS_IGNORE, ierr)
            
!~             BandPoints_All% buffer(nProcs) = BandPoints_All% buffer(nProcs-1) + ObjsSize
            BandRegion% buffer(nProcs) = BandRegion% buffer(nProcs-1) + ObjsSize
            buffer(nProcs-1) = ObjsSize
            
         end do         
                   
         do nProcs = 2, MPI_Process% nProcs
      
             allocate( COORD_x(buffer(nProcs-1)),  &
                       COORD_y(buffer(nProcs-1)),  &
                       COORD_z(buffer(nProcs-1)),  &
                       i_v(buffer(nProcs-1)),      &    
                       j_v(buffer(nProcs-1)),      &
                       k_v(buffer(nProcs-1)),      &
                       eID(buffer(nProcs-1)),      & 
                       partition(buffer(nProcs-1)) ) 
                                  
             call mpi_irecv( COORD_x, buffer(nProcs-1), MPI_DOUBLE, nProcs-1,   &
                             MPI_ANY_TAG, MPI_COMM_WORLD, RootBandrecv_req(nProcs-1,2), ierr )  
          
             call mpi_irecv( COORD_y, buffer(nProcs-1), MPI_DOUBLE, nProcs-1,  &
                             MPI_ANY_TAG, MPI_COMM_WORLD, RootBandrecv_req(nProcs-1,3), ierr ) 

             call mpi_irecv( COORD_z, buffer(nProcs-1), MPI_DOUBLE, nProcs-1,   &
                             MPI_ANY_TAG, MPI_COMM_WORLD, RootBandrecv_req(nProcs-1,4), ierr )   

             call mpi_irecv( i_v, buffer(nProcs-1), MPI_INT, nProcs-1,   &
                             MPI_ANY_TAG, MPI_COMM_WORLD, RootBandrecv_req(nProcs-1,5), ierr )   

             call mpi_irecv( j_v, buffer(nProcs-1), MPI_INT, nProcs-1,   &
                             MPI_ANY_TAG, MPI_COMM_WORLD, RootBandrecv_req(nProcs-1,6), ierr )   

             call mpi_irecv( k_v, buffer(nProcs-1), MPI_INT, nProcs-1,   &
                             MPI_ANY_TAG, MPI_COMM_WORLD, RootBandrecv_req(nProcs-1,7), ierr )   

             call mpi_irecv( eID, buffer(nProcs-1), MPI_INT, nProcs-1,   &
                             MPI_ANY_TAG, MPI_COMM_WORLD, RootBandrecv_req(nProcs-1,8), ierr )   

             call mpi_irecv( partition, buffer(nProcs-1), MPI_INT, nProcs-1,   &
                             MPI_ANY_TAG, MPI_COMM_WORLD, RootBandrecv_req(nProcs-1,9), ierr )   
                             
             call mpi_waitall(9, RootBandrecv_req(nProcs-1,:), array_of_statuses, ierr)                                                
             
             do i = 1, buffer(nProcs-1)               
                BandRegion% x(i+BandRegion% buffer(nProcs-1))% index     = i+BandPoints_All% buffer(nProcs-1)
                BandRegion% x(i+BandRegion% buffer(nProcs-1))% coords(1) = COORD_x(i)
                BandRegion% x(i+BandRegion% buffer(nProcs-1))% coords(2) = COORD_y(i)
                BandRegion% x(i+BandRegion% buffer(nProcs-1))% coords(3) = COORD_z(i)  
                BandRegion% x(i+BandRegion% buffer(nProcs-1))% local_Position(1) = i_v(i)  
                BandRegion% x(i+BandRegion% buffer(nProcs-1))% local_Position(2) = j_v(i)  
                BandRegion% x(i+BandRegion% buffer(nProcs-1))% local_Position(3) = k_v(i)  
                BandRegion% x(i+BandRegion% buffer(nProcs-1))% element_index = eID(i)  
                BandRegion% x(i+BandRegion% buffer(nProcs-1))% partition = partition(i)  
             end do
                      
             deallocate( COORD_x, COORD_y, COORD_z, i_v, j_v, k_v, eID, partition )

          end do
          
          deallocate(buffer,RootBandrecv_req)

       end if
#endif          
   
   end subroutine RootrecvBandPoint
   
   
   subroutine RootSendBandPoint( BandPoints )
   
      implicit none
      !-arguments----------------------------------------------------------------------
      type(PointLinkedList), intent(inout) :: BandPoints
#ifdef _HAS_MPI_  
      !-local-variables----------------------------------------------------------------
      integer                                    :: ObjsSize, i, ierr, &
                                                    array_of_statuses(MPI_STATUS_SIZE,9), &
                                                    RootBandsend_req(9)
      type(point_type), pointer                  :: p
      real(kind=RP), dimension(:),   allocatable :: COORD_x, COORD_y, COORD_z    
      integer,       dimension(:),   allocatable :: i_v, j_v, k_v, eID, partition
 
      if( MPI_Process% isRoot ) return
      
      ObjsSize = BandPoints% NumOfPoints
       
      call mpi_isend( ObjsSize, 1, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD, RootBandsend_req(1), ierr ) !

      allocate( COORD_x(ObjsSize),  &
                COORD_y(ObjsSize),  &
                COORD_z(ObjsSize),  &
                i_v(ObjsSize),      &    
                j_v(ObjsSize),      &
                k_v(ObjsSize),      &
                eID(ObjsSize),      &
                partition(ObjsSize) )

      p => BandPoints% head
      
      do i = 1, ObjsSize 
         COORD_x(i) = p% coords(1)
         COORD_y(i) = p% coords(2)
         COORD_z(i) = p% coords(3)
         i_v(i) = p% local_Position(1)
         j_v(i) = p% local_Position(2)
         k_v(i) = p% local_Position(3)
         eID(i) = p% element_index
         partition(i) = p% partition
         p => p% next
      end do

      call mpi_isend( COORD_x, ObjsSize, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, RootBandsend_req(2), ierr )
  
      call mpi_isend( COORD_y, ObjsSize, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, RootBandsend_req(3), ierr )

      call mpi_isend( COORD_z, ObjsSize, MPI_DOUBLE, 0, DEFAULT_TAG, MPI_COMM_WORLD, RootBandsend_req(4), ierr )
      
      call mpi_isend( i_v, ObjsSize, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD, RootBandsend_req(5), ierr )
                         
      call mpi_isend( j_v, ObjsSize, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD, RootBandsend_req(6), ierr )
                         
      call mpi_isend( k_v, ObjsSize, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD, RootBandsend_req(7), ierr )
                         
      call mpi_isend( eID, ObjsSize, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD, RootBandsend_req(8), ierr  )
      
      call mpi_isend( partition, ObjsSize, MPI_INT, 0, DEFAULT_TAG, MPI_COMM_WORLD, RootBandsend_req(9), ierr  )
           
      call mpi_waitall(9, RootBandsend_req, array_of_statuses, ierr)                     
        
      deallocate( COORD_x,  &
                  COORD_y,  &
                  COORD_z,  &
                  i_v,      &  
                  j_v,      &
                  k_v,      &
                  eID,      &
                  partition )  
#endif  

   end subroutine RootSendBandPoint 
 
   subroutine recvBandPointPartition( BandRegion )
    
      implicit none
      !-arguments-------------------------------------------------------------------------
      type(MPI_M_Points_type), intent(inout) :: BandRegion
#ifdef _HAS_MPI_      
      !-local-variables-------------------------------------------------------------------
      integer                                  :: ObjsSize, ierr, i,                    &
                                                  array_of_statuses(MPI_STATUS_SIZE,9), &
                                                  Bandrecv_req(9)
      real(kind=RP), dimension(:), allocatable :: COORD_x, COORD_y, COORD_z
      integer,       dimension(:), allocatable :: indeces, i_v, j_v, k_v, eID, partition
      
      if( MPI_Process% isRoot ) return 

      ObjsSize = BandRegion% NumOfObjs

      allocate( indeces(ObjsSize),  &
                COORD_x(ObjsSize),  &
                COORD_y(ObjsSize),  &
                COORD_z(ObjsSize),  &
                i_v(ObjsSize),      &
                j_v(ObjsSize),      &
                k_v(ObjsSize),      &
                eID(ObjsSize),      &
                partition(ObjsSize) )      
                            
      call mpi_irecv( indeces, ObjsSize, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, Bandrecv_req(1), ierr )  
      
      call mpi_irecv( COORD_x, ObjsSize, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, Bandrecv_req(2), ierr )  

      call mpi_irecv( COORD_y, ObjsSize, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, Bandrecv_req(3), ierr )  

      call mpi_irecv( COORD_z, ObjsSize, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, Bandrecv_req(4), ierr )    

      call mpi_irecv( i_v, ObjsSize, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, Bandrecv_req(5), ierr ) 
         
      call mpi_irecv( j_v, ObjsSize, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, Bandrecv_req(6), ierr )    
      
      call mpi_irecv( k_v, ObjsSize, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, Bandrecv_req(7), ierr )    
      
      call mpi_irecv( eID, ObjsSize, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, Bandrecv_req(8), ierr )    
      
      call mpi_irecv( partition, ObjsSize, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, Bandrecv_req(9), ierr )    

      call mpi_waitall(9, Bandrecv_req, array_of_statuses, ierr)     

      do i = 1, ObjsSize  
         BandRegion% x(i)% index     = indeces(i)
         BandRegion% x(i)% coords(1) = COORD_x(i)
         BandRegion% x(i)% coords(2) = COORD_y(i)
         BandRegion% x(i)% coords(3) = COORD_z(i)  
         BandRegion% x(i)% local_Position(1) = i_v(i)  
         BandRegion% x(i)% local_Position(2) = j_v(i)  
         BandRegion% x(i)% local_Position(3) = k_v(i)  
         BandRegion% x(i)% element_index = eID(i)  
         BandRegion% x(i)% partition = partition(i)  
      end do 
 
      deallocate( indeces, COORD_x, COORD_y, COORD_z, i_v, j_v, k_v, eID, partition )
#endif
   
   end subroutine recvBandPointPartition

   subroutine sendBandPointPartition( BandRegion )
   
      implicit none
      !-arguments-------------------------------------------------------------------------
      type(MPI_M_Points_type), intent(inout) :: BandRegion
#ifdef _HAS_MPI_  
      !-local-variables--------------------------------------------------------------------------
      integer                                    :: ObjsSize, nProcs, msg, ierr, &
                                                    array_of_statuses(MPI_STATUS_SIZE,9)
      real(kind=RP), dimension(:),   allocatable :: COORD_x, COORD_y, COORD_z
      integer,       dimension(:),   allocatable :: indeces, i_v, j_v, k_v, eID, partition
      integer,       dimension(:,:), allocatable :: Bandsend_req
      integer                                    :: i

      ObjsSize = BandRegion% NumOfObjs
      
      allocate( indeces(ObjsSize),                    &
                COORD_x(ObjsSize),                    & 
                COORD_y(ObjsSize),                    & 
                COORD_z(ObjsSize),                    &
                i_v(ObjsSize),                        &
                j_v(ObjsSize),                        &
                k_v(ObjsSize),                        &
                eID(ObjsSize),                        &
                partition(ObjsSize),                  &
                Bandsend_req(MPI_Process% nProcs-1,9) )

      indeces = BandRegion% x(:)% index
      COORD_x = BandRegion% x(:)% coords(1)
      COORD_y = BandRegion% x(:)% coords(2)
      COORD_z = BandRegion% x(:)% coords(3)
      i_v = BandRegion% x(:)% local_Position(1)
      j_v = BandRegion% x(:)% local_Position(2)
      k_v = BandRegion% x(:)% local_Position(3)
      eID = BandRegion% x(:)% element_index
      partition = BandRegion% x(:)% partition
   
      do nProcs = 2, MPI_Process% nProcs        
      
         call mpi_isend( indeces, ObjsSize, MPI_INT, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, &
                         Bandsend_req(nProcs-1,1), ierr                                     )
      
         call mpi_isend( COORD_x, ObjsSize, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, &
                         Bandsend_req(nProcs-1,2), ierr                                        )
      
         call mpi_isend( COORD_y, ObjsSize, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, &
                         Bandsend_req(nProcs-1,3), ierr                                        )
      
         call mpi_isend( COORD_z, ObjsSize, MPI_DOUBLE, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, &
                         Bandsend_req(nProcs-1,4), ierr                                        )
      
         call mpi_isend( i_v, ObjsSize, MPI_INT, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, &
                         Bandsend_req(nProcs-1,5), ierr                                 )
      
         call mpi_isend( j_v, ObjsSize, MPI_INT, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, &
                         Bandsend_req(nProcs-1,6), ierr                                 )
      
         call mpi_isend( k_v, ObjsSize, MPI_INT, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, &
                         Bandsend_req(nProcs-1,7), ierr                                 )
      
         call mpi_isend( eID, ObjsSize, MPI_INT, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, &
                         Bandsend_req(nProcs-1,8), ierr                                 )
      
         call mpi_isend( partition, ObjsSize, MPI_INT, nProcs-1, DEFAULT_TAG, MPI_COMM_WORLD, &
                         Bandsend_req(nProcs-1,9), ierr                                       )
                         
         call mpi_waitall( 9, Bandsend_req(nProcs-1,:), array_of_statuses, ierr )                
         
      end do
                         
      deallocate( indeces, COORD_x, COORD_y, COORD_z, i_v, j_v, k_v, eID, partition, Bandsend_req )
#endif   
   
   end subroutine sendBandPointPartition

   
end module MPI_IBMUtilities
