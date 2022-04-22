module KDClass

   use SMConstants
   use Utilities
   use TessellationTypes
   use OrientedBoundingBox

   implicit none

!~          5-----------------8
!~         /|                /|
!~        / |               / |
!~       /  |              /  |
!~      6-----------------7   |
!~      |   |             |   |
!~      |   |             |   |
!~      |   |             |   |
!~      |   |             |   |
!~      |   1-------------|---4                  
!~      |  /              |  /                   | z
!~      | /               | /                    |_______  y
!~      |/                |/                    /
!~      2-----------------3                    / x
       
  public :: LastLevel, BoxIndex

  integer :: BoxIndex, LastLevel = -1, depth, NumOfObjsIndx, NumOfKDtreePoints
  
  integer, parameter :: side_L = 4, side_R = 0, ON_PLANE = 0, POINTS_KDTREE = 0, TRIANGLES_KDTREE = 1
  integer, parameter :: START_ = 2, END_ = 0, PLANAR_ = 1, ONLYLEFT = 0, ONLYRIGHT = 1, BOTH = 2, ODD = 1, EVEN = 2
  integer, parameter :: BREADTHFIRST = 1, DEPTHFIRST = 0
  
  real(kind=rp), parameter :: C_TRANSVERSE = 1.0_RP, C_INTERSECT = 1.5_RP, C_1 = 1.2_RP, C_2 = 2.0_RP
  
  integer,  parameter, dimension(8) :: vertices_x = (/ 1,4,5,8,  2,3,6,7 /), &
                                       vertices_y = (/ 1,2,6,5,  3,4,7,8 /), &
                                       vertices_z = (/ 1,2,3,4,  5,6,7,8 /)  
                                       
  type Event
  
     real(kind=rp) :: plane
     integer       :: eType, index, axis
  
  end type
                            
!
!  **************************************************
!  Main type for a kd-tree
!  **************************************************
   type KDtree

      class(KDtree), pointer                         :: child_L, child_R, parent
      type(Object_type), dimension(:),   allocatable :: ObjectsList 
      real(kind=rp),     dimension(3,8)              :: vertices
      integer                                        :: NumOfObjs, level, axis, &
                                                        index, Min_n_of_Objs,   &
                                                        which_KDtree, MaxAxis,  &
                                                        SIDE, N_L, N_R, N_B,    &
                                                        HalfEvents, NumThreads, &
                                                        STLNum
      logical                                        :: isLast, Split
      integer,           dimension(NDIM)             :: NumOfEvents
      real(kind=rp)                                  :: S, SplitCost, SplittingPlane
      type(Event),       dimension(:,:), allocatable :: Events
      integer,           dimension(:),   allocatable :: ObjsIndeces
      
      contains
         procedure :: construct         => KDtree_construct
         procedure :: SetUpRoot         => KDtree_SetUpRoot
         procedure :: FindLeaf          => KDtree_FindLeaf 
         procedure :: plot              => KDtree_plot
         procedure :: plotBlock         => KDtree_plotBlock
         procedure :: Destruct          => KD_treeDestruct
         procedure :: GetArea           => KD_treeGetArea
         procedure :: BuildChild        => KDtree_BuildChild
         procedure :: EvaluateCost      => KDtree_EvaluateCost
         procedure :: SaveObjsIndeces   => KDtree_SaveObjsIndeces
         procedure :: SavePointsIndeces => KDtree_SavePointsIndeces

   end type

   type DepthFirst_type
   
      type(KDtree), pointer                     :: taskTree
      type(Event),  dimension(:,:), allocatable :: taskEvents
      logical                                   :: active = .false.
   
   end type 
   
   type taskPart_dim_type
   
      real(kind=RP) :: splittingPlane
      integer       :: N_L, N_R, N_P, ChunkDim, Starting_index, Values2Check
      
   end type 
   
   type taskSAH_type
    
      real(kind=RP) :: SplitCost, splittingPlane
      integer       :: axis, SIDE
      logical       :: split = .false.
   
   end type 
   
   type taskPart_type
   
      type(taskPart_dim_type), dimension(NDIM) :: taskPartDim

   end type 

   contains 
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  --------------------------------------------------------------
! This function computes the index corresponding to the integer i. 
!  --------------------------------------------------------------
 
   integer function AxisIndex( i ) result( axis )
 
      implicit none
      !-arguments-------------
      integer, intent(in) :: i
       
      if( i .eq. 1 ) axis = 2
      if( i .eq. 2 ) axis = 3
      if( i .eq. 3 ) axis = 1
 
   end function AxisIndex
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This function computes the longest axis of the OBB. 
!  -------------------------------------------------
   integer function FindAxis( this ) result( FirstAxis )
   
      implicit none
      !-arguments
      type(KDtree), intent(in) :: this
      !-local-variables----------------
      integer,  dimension(1) :: maxvec
      
      maxvec = maxloc( (/ OBB(this% STLNum)% MBR% Length, OBB(this% STLNum)% MBR% Width, &
                          abs(OBB(this% STLNum)% nMax) + abs(OBB(this% STLNum)% nMin) /) )
       
      FirstAxis = maxvec(1)
 
   end function FindAxis
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  ----------------------------------------------------------
! This subroutine saves and plots the leaves of the KD tree. 
!  ----------------------------------------------------------
   
   recursive subroutine Plot_leaves( tree, STLNum, funit, which_KDtree )
      use PhysicsStorage
      implicit none
      !-arguments------------------------------
      type(KDtree),                    intent(inout) :: tree
      integer,                         intent(in)    :: STLNum, funit, &
                                                        which_KDtree
      !-local-variables------------------------
      real(kind=rp), dimension(NDIM) :: x_g
      integer                        :: i
      
      if( tree% isLast ) then 
         write(funit,"(a69)") 'ZONE NODES=8, ELEMENTS = 6, DATAPACKING=POINT, ZONETYPE=FETETRAHEDRON'
         select case( which_KDtree )
            case( POINTS_KDTREE )
               do i = 1, 8
                  write(funit,'(3E13.5)') Lref*tree% vertices(1,i),Lref*tree% vertices(2,i), Lref*tree% vertices(3,i)
               end do            
            case( TRIANGLES_KDTREE )
               do i = 1, 8
                  call OBB(STLNum)% ChangeRefFrame( tree% vertices(:,i), 'global', x_g )
                  write(funit,'(3E13.5)') Lref*x_g(1),Lref*x_g(2), Lref*x_g(3)
               end do 
         end select
         
         write(funit,'(4i2)') 1, 2, 3, 4
         write(funit,'(4i2)') 1, 5, 8, 4
         write(funit,'(4i2)') 5, 6, 7, 8
         write(funit,'(4i2)') 2, 3, 7, 6
         write(funit,'(4i2)') 4, 8, 7, 3
         write(funit,'(4i2)') 1, 2, 6, 5
         
      else
         call Plot_leaves( tree% child_L, STLNum, funit, which_KDtree )
         call Plot_leaves( tree% child_R, STLNum, funit, which_KDtree )
      end if
        
   end subroutine Plot_leaves      
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine creates a new file call KDtree in which the leaves will be stored. 
!  ------------------------------------------------  
   subroutine KDtree_plot( this )
      use MPI_Process_Info
      implicit none
      !-arguemnts----------------------------------
      class(KDtree), target, intent(inout) :: this
      !-local-variables----------------------------
      real(kind=rp), dimension(NDIM) :: x_g
      character(len=LINE_LENGTH) :: filename, myString
      integer :: i, funit

      select case( this% which_KDtree )
         case( POINTS_KDTREE )
            if( MPI_Process% isRoot ) then
               ! same kdtree for all the processes
               funit = UnusedUnit()
            
               open(funit,file='IBM/KDTreeBandPoints.tec', status='unknown')
               
               write(funit,"(a28)") 'TITLE = "KD-tree"'
               write(funit,"(a25)") 'VARIABLES = "x", "y", "z"'

               call Plot_leaves( this, this% STLNum, funit, this% which_KDtree )

               close(funit)
               
            end if
         case( TRIANGLES_KDTREE )
         
            write(myString,'(i100)') MPI_Process% rank
      
            if( MPI_Process% nProcs .eq. 1 ) then
               filename = OBB(this% STLNum)% filename
            else
               filename = 'PARTITION'//trim(adjustl(myString))//'_'//OBB(this% STLNum)% filename
               filename = trim(filename)
            end if
     
            funit = UnusedUnit()
     
            open(funit,file='IBM/KDTree_'//trim(filename)//'.tec', status='unknown')
         
            write(funit,"(a28)") 'TITLE = "KD-tree"'
            write(funit,"(a25)") 'VARIABLES = "x", "y", "z"'

            call Plot_leaves( this, this% STLNum, funit, this% which_KDtree )

            close(funit)
         
         case default
            print *, " KD tree type not recognized, stop."
            error stop         
      end select
   
   end subroutine KDtree_plot
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine plot a single block whose index is BlockIndex. if PlotObjs is true
! also the objects belonging to the box are saved and plot 
!  ------------------------------------------------
   subroutine KDtree_plotBlock( this, ObjectsList, PlotObjs )
   
      implicit none
      !-arguments---------------------------------------------
      class(KDtree),  intent(inout) :: this
      type(object_type), intent(in) :: ObjectsList(:)
      logical,        intent(in)    :: PlotObjs
      !-local-variables---------------------------------------
      real(kind=rp), dimension(NDIM) :: x_g
      character(len=LINE_LENGTH)     :: filename, myString
      integer :: i,j, funit, index
      
      funit = UnusedUnit()

      write(myString,'(i100)') this% index
      filename = 'block'//trim(adjustl(myString))//'_'//OBB(this% STLNum)% filename
      filename = trim(filename)
      
      open(funit,file='IBM/'//trim(filename)//'.tec', status='unknown')
 
      write(funit,"(a28)") 'TITLE = "BLOCK"'
      write(funit,"(a25)") 'VARIABLES = "x", "y", "z"'
      
      write(funit,"(a69)") 'ZONE NODES=8, ELEMENTS = 6, DATAPACKING=POINT, ZONETYPE=FETETRAHEDRON'
      do i = 1, 8
         call OBB(this% STLNum)% ChangeRefFrame( this% vertices(:,i), 'global', x_g )
         write(funit,'(3E13.5)') x_g(1),x_g(2), x_g(3)
      end do 

      write(funit,'(4i2)') 1, 2, 3, 4
      write(funit,'(4i2)') 1, 5, 8, 4
      write(funit,'(4i2)') 5, 6, 7, 8
      write(funit,'(4i2)') 2, 3, 7, 6
      write(funit,'(4i2)') 4, 8, 7, 3
      write(funit,'(4i2)') 1, 2, 6, 5
      
      close(funit)
      
      if( PlotObjs .and. this% NumOfObjs .gt. 0 ) then
        filename = 'Objects'//trim(adjustl(myString))//'_'//OBB(this% STLNum)% filename
        filename = trim(filename)
      
         open(funit,file='IBM/'//trim(filename)//'.tec', status='unknown')
 
         write(funit,"(a28)") 'TITLE = "ObjsIndx"'
         write(funit,"(a25)") 'VARIABLES = "x", "y", "z"'
      
         do i = 1, this% NumOfObjs
            write(funit,"(a66)") 'ZONE NODES=3, ELEMENTS = 1, DATAPACKING=POINT, ZONETYPE=FETRIANGLE'
            index = this% ObjsIndeces(i)
      
            do j = 1, ObjectsList(index)% NumOfVertices
               call OBB(this% STLNum)% ChangeRefFrame( ObjectsList(index)% vertices(j)% coords, 'global', x_g )
               write(funit,'(3E13.5)') x_g(1),x_g(2), x_g(3)
            end do

            write(funit,'(3i2)') 1, 2, 3 
         end do
             
         close(funit)
      end if
      
   end subroutine KDtree_plotBlock
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine sets up the starting box of the tree which is called root. 
! It coincides whit the OBB if PointList not present. 
!  ------------------------------------------------   
   subroutine KDtree_SetUpRoot( this, stl, Vertices, PointList )
   
      implicit none
      !-arguments--------------------------------------------
      class(KDtree),                       intent(inout) :: this
      type(STLfile),                       intent(in)    :: stl
      real(kind=RP),    dimension(NDIM,8), intent(in)    :: Vertices
      type(Point_type), dimension(:),      intent(in)    :: PointList
      !-local-variables--------------------------------------
      type(Object_type), pointer     :: vertex => null()
      type(Object_type)              :: Obj, Objvertex  
      type(point_type),  pointer     :: p => null()
      real(kind=rp), dimension(NDIM) :: coords
      real(kind=rp)                  :: MaxX, MaxY, MaxZ, MinX, MinY, MinZ
      integer                        :: i
      
      optional :: PointList
      
      this% NumOfObjs = 0
      
      if( present(PointList) ) then
         MaxX = maxval(PointList(:)% coords(1)); MinX = minval(PointList(:)% coords(1))
         MaxY = maxval(PointList(:)% coords(2)); MinY = minval(PointList(:)% coords(2))
         MaxZ = maxval(PointList(:)% coords(3)); MinZ = minval(PointList(:)% coords(3))
        
         this% NumOfObjs = size(PointList)

         this% vertices(:,1) = (/ MinX, MinY, MinZ /)
         this% vertices(:,2) = (/ MaxX, MinY, MinZ /)
         this% vertices(:,3) = (/ MaxX, MaxY, MinZ /)
         this% vertices(:,4) = (/ MinX, MaxY, MinZ /)  
         this% vertices(:,5) = (/ MinX, MinY, MaxZ /)
         this% vertices(:,6) = (/ MaxX, MinY, MaxZ /)
         this% vertices(:,7) = (/ MaxX, MaxY, MaxZ /)
         this% vertices(:,8) = (/ MinX, MaxY, MaxZ /)

      else      
!
!        Local ref frame
!        ----------------------------------------------------
         do i = 1, 4
            this% vertices(:,i)   = Vertices(:,i)
            this% vertices(:,i+4) = Vertices(:,i+4)
         end do
                  
         if( associated(this% parent) ) nullify(this% parent)
         allocate( this% ObjectsList(stl% NumOfObjs) )
         
         associate( Objs => stl% ObjectsList )

         do i = 1, stl% NumOfObjs
            this% ObjectsList(i) = Objs(i)            
         end do

         end associate

         this% NumOfObjs = stl% NumOfObjs

      end if
   
   end subroutine KDtree_SetUpRoot
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine builds the KD tree. 
!  ------------------------------------------------   
   subroutine KDtree_construct( this, stl, Vertices, isPlot, Min_n_of_Objs, PointList )
      use omp_lib
      implicit none
      !-arguments-----------------------------------------------------------------
      class(KDtree),                       intent(inout) :: this
      type(STLfile),                       intent(in)    :: stl
      real(kind=RP),    dimension(NDIM,8), intent(in)    :: Vertices
      logical,                             intent(in)    :: isPlot
      integer,                             intent(in)    :: Min_n_of_Objs
      type(point_type), dimension(:),      intent(in)    :: PointList
      !-local-varables------------------------------------------------------------
      real(kind=rp)                                      :: NumOfObjs,        &
                                                            NumThreads_RP
      integer,               dimension(:),   allocatable :: ObjsIndx
      type(Event),           dimension(:,:), allocatable :: Events
      type(DepthFirst_type), dimension(:),   allocatable :: Depth_First
      integer                                            :: i, k, NumThreads, &
                                                            DepthFirstLevel,  &
                                                            NumDepthFirst

      optional :: PointList 
      
      if( present( PointList ) ) then
         this% which_KDtree = POINTS_KDTREE
         call this% SetUpRoot( stl, vertices, PointList )   
      else
         this% which_KDtree = TRIANGLES_KDTREE
         call this% SetUpRoot( stl, vertices )
      end if
      
      this% level         = 0
      NumOfObjs           = this% NumOfObjs
      depth               = C_1*log(NumOfObjs)/log(2.0_RP) + C_2!1.2 * log2(N) + 2 --> on improving kd tree for ray shooting Havran
      this% Min_n_of_Objs = Min_n_of_Objs
      this% isLast        = .false.
      BoxIndex            = 1
      this% index         = BoxIndex
      
      select case( this% which_KDtree )
         case( POINTS_KDTREE )
         
            this% axis = FindAxis( this )
         
            allocate(Events(NDIM,this% NumOfObjs))

            call GetPointsEvents( this, Events, PointList, NumThreads )

            this% NumThreads  = NumThreads
            NumThreads_RP     = NumThreads 
            DepthFirstLevel   = floor(log(NumThreads_RP)/log(2.0_RP)) 
            NumDepthFirst     = 2.0_RP**(DepthFirstLevel)
            
            allocate(Depth_First(NumThreads),ObjsIndx(this% NumOfObjs))

            call KDtree_buildPoints_BreadthFirst( this, Events, ObjsIndx, DepthFirstLevel, Depth_First )

            deallocate(ObjsIndx)

!$omp parallel shared(this,Depth_First,k) 
!$omp single 
            do k = 1, count( Depth_First(:)% active )
!$omp task firstprivate(k) private(ObjsIndx) 
               allocate(ObjsIndx(this% NumOfObjs))
               call KDtree_buildPoints_DepthFirst( Depth_First(k)% taskTree, Depth_First(k)% taskEvents, ObjsIndx )
               nullify( Depth_First(k)% taskTree )
               deallocate(ObjsIndx)
!$omp end task
            end do
!$omp end single
!$omp end parallel

            deallocate(Depth_First)

         case( TRIANGLES_KDTREE )

            allocate(Events(NDIM,2*this% NumOfObjs))      
      
            call GetEvents( this, Events, NumThreads )

            this% NumThreads  = NumThreads
            NumThreads_RP     = NumThreads 
            DepthFirstLevel   = floor(log(NumThreads_RP)/log(2.0_RP)) 
            NumDepthFirst     = 2.0_RP**(DepthFirstLevel)
            
            allocate(Depth_First(NumThreads),ObjsIndx(this% NumOfObjs))

            call KDtree_buildSAH_BreadthFirst( this, Events, ObjsIndx, DepthFirstLevel, Depth_First )

            deallocate(ObjsIndx)

!$omp parallel shared(this,Depth_First,k) 
!$omp single 
            do k = 1, count( Depth_First(:)% active )
!$omp task firstprivate(k) private(ObjsIndx) 
               allocate(ObjsIndx(this% NumOfObjs))
               call KDtree_buildSAH_DepthFirst( Depth_First(k)% taskTree, Depth_First(k)% taskEvents, ObjsIndx )
               nullify( Depth_First(k)% taskTree )
               deallocate(ObjsIndx)
!$omp end task
            end do
!$omp end single
!$omp end parallel
            deallocate(Depth_First)

         case default
            print *, 'KD tree type not recognized'
            error stop
      end select  

      if( isPlot ) call this% plot()

   end subroutine KDtree_construct
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine finds the leaf where Point lies
!  ------------------------------------------------   
   subroutine KDtree_FindLeaf( this, Point, tree ) 
   
      implicit none
      !-arguments------------------------------------------
      class(KDtree), target,       intent(inout) :: this
      real(kind=rp), dimension(:), intent(in)    :: Point
      type(KDtree), pointer,       intent(inout) :: tree
      !-local-variables------------------------------------
      real(kind=rp), dimension(NDIM) :: vertices
      integer                        :: level
   
      if( associated(this% parent) ) then
         print *, "The starting point of the KDtree is not the root" 
         error stop 
      end if
   
      tree => this
      
      level = 0
      
      do
         if(  Point(tree% axis) .le. tree% SplittingPlane ) then !lt
            tree => tree% child_L
         elseif( Point(tree% axis) .gt. tree% SplittingPlane ) then
            tree => tree% child_R
         end if             
         if( tree% isLast ) return
      end do
      
   end subroutine KDtree_FindLeaf

   recursive subroutine KD_treeDestruct( this ) 
   
      implicit none
      !-arguments----------------------------------
      class(KDtree), intent(inout) :: this
      !-local-variables----------------------------
      type(KDtree), pointer :: child_L, child_R 
      integer               :: i     
      
      if( associated(this% child_L) .and. .not. this% isLast ) call this% child_L% destruct      
      if( associated(this% child_R) .and. .not. this% isLast ) call this% child_R% destruct      

      if( allocated(this% ObjectsList) ) deallocate(this% ObjectsList)      
      if( allocated(this% ObjsIndeces) ) deallocate(this% ObjsIndeces)
      if( allocated(this% Events) )      deallocate(this% Events)
      
      if( associated(this% child_L) ) deallocate(this% child_L)
      if( associated(this% child_R) ) deallocate(this% child_R)

   end subroutine KD_treeDestruct  
   
   subroutine KD_treeGetArea( tree )
   
      implicit none
      
      class(KDtree), intent(inout) :: tree
      
      real(kind=rp) :: L_x, L_y, L_z
      
      L_x = abs(tree% vertices(1,7)-tree% vertices(1,1))
      L_y = abs(tree% vertices(2,7)-tree% vertices(2,1))
      L_z = abs(tree% vertices(3,7)-tree% vertices(3,1))
      
      tree% S = 2.0_RP*(L_x*L_y + L_x*L_z + L_y*L_z)
      
   end subroutine KD_treeGetArea
   
   subroutine GetEvents( tree, Events, NumThreads )
      use omp_lib
      implicit none
   
      type(KDtree),                intent(inout) :: tree
      type(Event), dimension(:,:), intent(inout) :: Events
      integer,                     intent(out)   :: NumThreads
   
      integer :: k
   
!$omp parallel shared(tree,Events,k)
!$omp single
            NumThreads = omp_get_num_threads()
            do k = 1, NDIM
!$omp task firstprivate(k)
               call Event_Construct( Events(k,:), tree% ObjectsList, tree, k ) 
               call QsortEvents( Events(k,:), 1, tree% NumOfEvents(k) )  
!$omp end task 
            end do           
!$omp end single
!$omp end parallel
   
   end subroutine GetEvents 
   
   subroutine GetPointsEvents( tree, Events, PointList, NumThreads )
      use omp_lib
      implicit none
   
      type(KDtree),                     intent(inout) :: tree
      type(Event),      dimension(:,:), intent(inout) :: Events
      type(point_type), dimension(:),   intent(in)    :: PointList
      integer,                          intent(out)   :: NumThreads
   
      integer :: k
      
!$omp parallel shared(tree,Events,k,PointList)
           NumThreads = omp_get_num_threads()
!$omp single
            do k = 1, NDIM
!$omp task firstprivate(k)
               call Points_Event_Construct( Events(k,:), tree, k, PointList )
               call QsortEvents( Events(k,:), 1, tree% NumOfEvents(k) )
!$omp end task
            end do        
!$omp end single
!$omp end parallel
   
   end subroutine GetPointsEvents 
   
   recursive subroutine KDtree_buildSAH_BreadthFirst( this, Events, ObjsIndx, level, Depth_First )
   
      implicit none
      
      type(KDtree), target,       intent(inout) :: this
      type(Event),  allocatable,  intent(inout) :: Events(:,:)
      integer,                    intent(in)    :: level
      integer,      dimension(:), intent(inout) :: ObjsIndx 
      type(DepthFirst_type),      intent(inout) :: Depth_First(:)
      
      type(KDtree),   pointer :: child_L, child_R 
      type(Event), dimension(:,:), allocatable :: Events_L, Events_R
      integer :: j
   
      call this% GetArea()

      BoxIndex = BoxIndex + 1

      this% index  = BoxIndex

      if( this% level .lt. level .and. this% NumOfObjs .gt. 0 ) then

         call this% EvaluateCost( Events, BREADTHFIRST )

         if( this% split  ) then 

            allocate(this% child_L,this% child_R)
            allocate(this% child_L% parent, this% child_R% parent)
         
            this% child_L% parent => this; this% child_R% parent => this
        
            child_L => this% child_L
            child_R => this% child_R
           
            child_L% level = this% level + 1
            child_R% level = this% level + 1

            this% isLast = .false.   

            call this% BuildChild( child_L, side_L )
            call this% BuildChild( child_R, side_R )
         
            child_L% Min_n_of_Objs = this% Min_n_of_Objs
            child_R% Min_n_of_Objs = this% Min_n_of_Objs
            child_L% NumThreads = this% NumThreads
            child_R% NumThreads = this% NumThreads

            call Event_ClassifyObjs( Events(this% axis,:), this, child_L, child_R, ObjsIndx )

            allocate(Events_L(NDIM,2*child_L% NumOfObjs))
            allocate(Events_R(NDIM,2*child_R% NumOfObjs))

            call Event_BuildLists( Events, this, child_L, child_R, ObjsIndx, Events_L, Events_R, BREADTHFIRST )

            deallocate(Events)
            
            call KDtree_buildSAH_BreadthFirst( child_L, Events_L, ObjsIndx, level, Depth_First )
            call KDtree_buildSAH_BreadthFirst( child_R, Events_R, ObjsIndx, level, Depth_First )
            
         else
            this% isLast = .false.
            if( this% level .lt. level ) this% isLast = .true.
            do j = 1, size(Depth_First)
               if( .not. Depth_First(j)% active ) then 
                  Depth_First(j)% taskTree => this
                  allocate(Depth_First(j)% taskEvents(NDIM,2*this% NumOfObjs))
                  Depth_First(j)% taskEvents = Events
                  Depth_First(j)% active = .true.
                  deallocate(Events)
                  exit
               end if 
            end do  
         end if
      else
         this% isLast = .false.
         if( this% NumOfObjs .eq. 0 ) this% isLast = .true.
         do j = 1, size(Depth_First)
            if( .not. Depth_First(j)% active ) then 
               Depth_First(j)% taskTree => this
               allocate(Depth_First(j)% taskEvents(NDIM,2*this% NumOfObjs))               
               Depth_First(j)% taskEvents = Events
               Depth_First(j)% active = .true.
               deallocate(Events)
               exit
            end if
         end do
      end if
   
   end subroutine KDtree_buildSAH_BreadthFirst
   
   
   recursive subroutine KDtree_buildSAH_DepthFirst( this, Events, ObjsIndx )
   use omp_lib
      implicit none
      
      type(KDtree), target,       intent(inout) :: this
      type(Event),  allocatable,  intent(inout) :: Events(:,:)
      integer,      dimension(:), intent(inout) :: ObjsIndx 
      
      type(KDtree),   pointer :: child_L, child_R 
      type(Event), dimension(:,:), allocatable :: Events_L, Events_R

      call this% GetArea()

!$omp critical
      BoxIndex = BoxIndex + 1
!$omp end critical

      this% index  = BoxIndex

      if( this% level .lt. depth .and. this% NumOfObjs .ge. this% Min_n_of_Objs ) then

         call this% EvaluateCost( Events, DEPTHFIRST )

         if( this% split ) then 

            allocate(this% child_L,this% child_R)
            allocate(this% child_L% parent, this% child_R% parent)
            this% child_L% parent => this; this% child_R% parent => this
        
            child_L => this% child_L
            child_R => this% child_R
           
            child_L% level = this% level + 1
            child_R% level = this% level + 1

            this% isLast = .false.   

            call this% BuildChild( child_L, side_L )
            call this% BuildChild( child_R, side_R )
         
            child_L% Min_n_of_Objs = this% Min_n_of_Objs
            child_R% Min_n_of_Objs = this% Min_n_of_Objs
            child_L% NumThreads = this% NumThreads
            child_R% NumThreads = this% NumThreads

            call Event_ClassifyObjs( Events(this% axis,:), this, child_L, child_R, ObjsIndx )

            allocate(Events_L(NDIM,2*child_L% NumOfObjs))
            allocate(Events_R(NDIM,2*child_R% NumOfObjs))

            call Event_BuildLists( Events, this, child_L, child_R, ObjsIndx, Events_L, Events_R, DEPTHFIRST )

            deallocate(Events)

            call KDtree_buildSAH_DepthFirst( child_L, Events_L, ObjsIndx )
            call KDtree_buildSAH_DepthFirst( child_R, Events_R, ObjsIndx )
 
         else
            this% isLast = .true.
            call this% SaveObjsIndeces( Events )
            deallocate(Events)
         end if
      else
         this% isLast = .true.
         if( this% NumOfObjs .gt. 0 ) then
            call this% SaveObjsIndeces( Events )     
            deallocate(Events)
         end if
      end if

   end subroutine KDtree_buildSAH_DepthFirst
   
   
   
   subroutine ComputeSplitSurface_L(tree, SplittingPlane, axis, S_L)
   
      implicit none
      
      type(KDtree),  intent(in)  :: tree
      real(kind=rp), intent(in)  :: SplittingPlane
      integer,       intent(in)  :: axis
      real(kind=rp), intent(out) :: S_L
   
      real(kind=rp) :: L_x, L_y, L_z
   
      select case( axis )
         case(1)
            L_x = abs(SplittingPlane-tree% vertices(1,1))
            L_y = abs(tree% vertices(2,7)-tree% vertices(2,1))
            L_z = abs(tree% vertices(3,7)-tree% vertices(3,1))    
         case(2)
            L_x = abs(tree% vertices(1,7)-tree% vertices(1,1))
            L_y = abs(SplittingPlane-tree% vertices(2,1))
            L_z = abs(tree% vertices(3,7)-tree% vertices(3,1))
         case(3)
            L_x = abs(tree% vertices(1,7)-tree% vertices(1,1))
            L_y = abs(tree% vertices(2,7)-tree% vertices(2,1))
            L_z = abs(SplittingPlane-tree% vertices(3,1))
      end select

      S_L = 2.0_RP*(L_x*L_y+L_y*L_z+L_x*L_z)
   
   end subroutine ComputeSplitSurface_L
   
   subroutine ComputeSplitSurface_R(tree, SplittingPlane, axis, S_R)
   
      implicit none
      
      type(KDtree),  intent(in)  :: tree
      real(kind=rp), intent(in)  :: SplittingPlane
      integer,       intent(in)  :: axis
      real(kind=rp), intent(out) :: S_R
   
      real(kind=rp) :: L_x, L_y, L_z
   
      select case( axis )
         case(1)
            L_x = abs(tree% vertices(1,7)-SplittingPlane)
            L_y = abs(tree% vertices(2,7)-tree% vertices(2,1))
            L_z = abs(tree% vertices(3,7)-tree% vertices(3,1))    
         case(2)
            L_x = abs(tree% vertices(1,7)-tree% vertices(1,1))
            L_y = abs(tree% vertices(2,7)-SplittingPlane)
            L_z = abs(tree% vertices(3,7)-tree% vertices(3,1))
         case(3)
            L_x = abs(tree% vertices(1,7)-tree% vertices(1,1))
            L_y = abs(tree% vertices(2,7)-tree% vertices(2,1))
            L_z = abs(tree% vertices(3,7)-SplittingPlane)
      end select

      S_R = 2.0_RP*(L_x*L_y+L_y*L_z+L_x*L_z)
   
   end subroutine ComputeSplitSurface_R
   
   real(kind=rp) function computeSplittingCost(S, S_L, S_R, N_L, N_R) result(SplittingCost)
   
      implicit none
      
      real(kind=rp), intent(in) :: S, S_L, S_R
      integer,       intent(in) :: N_L, N_R
     
      SplittingCost = C_TRANSVERSE + C_INTERSECT*(S_L/S*N_L+S_R/S*N_R)
   
   end function computeSplittingCost
   
   subroutine KDtree_EvaluateCost( this, Events, parallelization_type )
      use omp_lib
      implicit none
      
      class(KDtree), intent(inout) :: this
      type(Event),   intent(in)    :: Events(:,:)
      integer,       intent(in)    :: parallelization_type
 
      integer                          :: i, j, k, l, N_L(NDIM), N_R(NDIM), N_P(NDIM), BestSide(1), pminus, pplus, pvert, axis
      real(kind=rp)                    :: SplittingCost(2), S_L, S_R, SplittingPlane
      type(taskPart_type), allocatable :: taskPart(:)
      integer                          :: level, taskNum
      type(KDtree),        allocatable :: taskTree(:)

      this% split = .false.
      this% SplitCost = C_INTERSECT * this% NumOfObjs
      
      if( this% NumOfObjs .eq. 0 ) return

      select case( parallelization_type )
      
      case( BREADTHFIRST ) 
      
      call constrct_taskPart( taskPart, this% NumOfEvents(:), this% NumThreads )
      allocate(taskTree(this% NumThreads))
      taskTree(:)% S         = this% S
      taskTree(:)% split     = .false.
      taskTree(:)% SplitCost = this% SplitCost
      
!$omp parallel shared(Events,taskPart,taskTree)
!$omp single 
      do taskNum = 1, size(taskPart)-1
!$omp task private(i,k,pplus,pminus,pvert,SplittingPlane) firstprivate(taskNum)
         do k = 1, NDIM
             pplus = 0; pminus = 0; pvert = 0
             SplittingPlane = Events(k,taskPart(taskNum)% taskPartDim(k)% Starting_index + &
                              taskPart(taskNum)% taskPartDim(k)% ChunkDim)% plane
             i = 1
             do while(i .le. taskPart(taskNum)% taskPartDim(k)% ChunkDim)
                               
                if( Events(k,taskPart(taskNum)% taskPartDim(k)% Starting_index+i)% eType .lt. 0 ) exit
                               
                do while( i .le. taskPart(taskNum)% taskPartDim(k)% ChunkDim .and.                                         &
                          (Events(k,taskPart(taskNum)% taskPartDim(k)% Starting_index+i)% plane .eq. SplittingPlane  .or.  &
                          Events(k,taskPart(taskNum)% taskPartDim(k)% Starting_index+i)% plane .lt. SplittingPlane ) .and. &
                          Events(k,taskPart(taskNum)% taskPartDim(k)% Starting_index+i)% eType .eq. END_                   )
                   pminus = pminus + 1
                   i = i + 1
                end do
          
                do while( i .le. taskPart(taskNum)% taskPartDim(k)% ChunkDim .and.                                         &
                         (Events(k,taskPart(taskNum)% taskPartDim(k)% Starting_index+i)% plane .eq. SplittingPlane   .or.  &
                          Events(k,taskPart(taskNum)% taskPartDim(k)% Starting_index+i)% plane .lt. SplittingPlane ) .and. &
                          Events(k,taskPart(taskNum)% taskPartDim(k)% Starting_index+i)% eType .eq. PLANAR_                )
                   pvert = pvert + 1
                   if(Events(k,taskPart(taskNum)% taskPartDim(k)% Starting_index+i)% plane .eq. SplittingPlane) &
                   taskPart(taskNum+1)% taskPartDim(k)% N_P = taskPart(taskNum+1)% taskPartDim(k)% N_P + 1
                   i = i + 1
                end do
              
                do while( i .le. taskPart(taskNum)% taskPartDim(k)% ChunkDim .and.                                         &
                          (Events(k,taskPart(taskNum)% taskPartDim(k)% Starting_index+i)% plane .eq. SplittingPlane  .or.  &
                          Events(k,taskPart(taskNum)% taskPartDim(k)% Starting_index+i)% plane .lt. SplittingPlane ) .and. &
                          Events(k,taskPart(taskNum)% taskPartDim(k)% Starting_index+i)% eType .eq. START_                 )                                    
                   pplus = pplus + 1
                   i = i + 1
                end do   
             end do
             taskPart(taskNum+1)% taskPartDim(k)% N_R = taskPart(taskNum+1)% taskPartDim(k)% N_R - pminus - pvert         
             taskPart(taskNum+1)% taskPartDim(k)% N_L = taskPart(taskNum+1)% taskPartDim(k)% N_L + pplus + pvert
             
         end do
!$omp end task
      end do
!$omp taskwait

     do taskNum = 1, size(taskPart)
        do k = 1, NDIM
           if( taskNum .eq. 1 ) then
              taskPart(taskNum)% taskPartDim(k)% N_R = this% NumOfObjs
           else
              taskPart(taskNum)% taskPartDim(k)% N_R = taskPart(taskNum)% taskPartDim(k)% N_R + &
                                                       taskPart(taskNum-1)% taskPartDim(k)% N_R
              taskPart(taskNum)% taskPartDim(k)% N_L = taskPart(taskNum)% taskPartDim(k)% N_L + &
                                                       taskPart(taskNum-1)% taskPartDim(k)% N_L
              taskPart(taskNum)% taskPartDim(k)% N_P = taskPart(taskNum)% taskPartDim(k)% N_P + &
                                                       taskPart(taskNum-1)% taskPartDim(k)% N_P
           end if
        end do
     end do

     do taskNum = 1, size(taskPart)
!$omp task private(i,k,pplus,pminus,pvert,SplittingPlane,N_R,N_L,N_P,SplittingCost,BestSide,S_L,S_R) firstprivate(taskNum)
        do k = 1, NDIM
           N_R(k) = taskPart(taskNum)% taskPartDim(k)% N_R
           N_L(k) = taskPart(taskNum)% taskPartDim(k)% N_L
           N_P(k) = taskPart(taskNum)% taskPartDim(k)% N_P   
           i = 1
           do while ( i .le. taskPart(taskNum)% taskPartDim(k)% ChunkDim )
              pplus = 0; pminus = 0; pvert = 0
              SplittingPlane = Events(k,taskPart(taskNum)% taskPartDim(k)% Starting_index+i)% plane

              if( Events(k,taskPart(taskNum)% taskPartDim(k)% Starting_index+i)% eType .lt. 0 ) exit

              do while( i .le. taskPart(taskNum)% taskPartDim(k)% ChunkDim .and.                                       &
                        Events(k,taskPart(taskNum)% taskPartDim(k)% Starting_index+i)% plane .eq. SplittingPlane .and. &
                        Events(k,taskPart(taskNum)% taskPartDim(k)% Starting_index+i)% eType .eq. END_                 )
                  pminus = pminus + 1
                  i = i + 1
               end do
         
               do while( i .le. taskPart(taskNum)% taskPartDim(k)% ChunkDim .and.                                       &
                         Events(k,taskPart(taskNum)% taskPartDim(k)% Starting_index+i)% plane .eq. SplittingPlane .and. &
                         Events(k,taskPart(taskNum)% taskPartDim(k)% Starting_index+i)% eType .eq. PLANAR_              )
                  pvert = pvert + 1
                  i = i + 1
               end do
         
               do while( i .le. taskPart(taskNum)% taskPartDim(k)% ChunkDim .and.                                       &
                         Events(k,taskPart(taskNum)% taskPartDim(k)% Starting_index+i)% plane .eq. SplittingPlane .and. &
                         Events(k,taskPart(taskNum)% taskPartDim(k)% Starting_index+i)% eType .eq. START_               )                                     
                  pplus = pplus + 1
                  i = i + 1
               end do

               N_R(k) = N_R(k) - pminus - pvert
   
               call ComputeSplitSurface_L(this, SplittingPlane, k, S_L)
               call ComputeSplitSurface_R(this, SplittingPlane, k, S_R)

               SplittingCost(1) = computeSplittingCost( taskTree(taskNum)% S, S_L, S_R, N_L(k)+N_P(k), N_R(k) )
               SplittingCost(2) = computeSplittingCost(taskTree(taskNum)% S, S_L, S_R, N_L(k), N_R(k)+N_P(k) )
               
               if( SplittingCost(1) .lt. SplittingCost(2) ) then
                  BestSide(1) = 1
               else
                  BestSide(1) = 2
               end if

               if( taskTree(taskNum)% SplitCost .gt. SplittingCost(BestSide(1)) ) then
                  taskTree(taskNum)% SplittingPlane = SplittingPlane
                  taskTree(taskNum)% SplitCost      = SplittingCost(BestSide(1))
                  taskTree(taskNum)% split          = .true.
                  taskTree(taskNum)% axis           = k
                  taskTree(taskNum)% SIDE           = BestSide(1) 
               end if

               N_L(k) = N_L(k) + pplus + pvert
               N_P(k) = pvert

            end do
         end do
!$omp end task
      end do
!$omp end single
!$omp end parallel

      if( any( taskTree(:)% split ) ) then
         BestSide = minloc(taskTree(:)% SplitCost)
         this% SplittingPlane = taskTree(BestSide(1))% SplittingPlane
         this% SplitCost      = taskTree(BestSide(1))% SplitCost
         this% split          = .true.
         this% axis           = taskTree(BestSide(1))% axis
         this% SIDE           = taskTree(BestSide(1))% SIDE
      end if
      
      deallocate(taskTree,taskPart)        
 
      case( DEPTHFIRST ) 
      
      N_R = this% NumOfObjs; N_L = 0; N_P = 0
      this% split = .false.
      this% SplitCost = C_INTERSECT * this% NumOfObjs      
      
      do k = 1, NDIM
         i = 1
         do while ( i .le. this% NumOfEvents(k) )
            pplus = 0; pminus = 0; pvert = 0
            SplittingPlane = Events(k,i)% plane

            if( Events(k,i)% eType .lt. 0 ) exit

            do while( i .le. this% NumOfEvents(k)  .and. &
                      Events(k,i)% plane .eq. SplittingPlane .and. Events(k,i)% eType .eq. END_ )
               pminus = pminus + 1
               i = i + 1
            end do
        
            do while( i .le. this% NumOfEvents(k) .and.  &
                     Events(k,i)% plane .eq. SplittingPlane .and. Events(k,i)% eType .eq. PLANAR_ )
               pvert = pvert + 1
               i = i + 1
            end do
         
            do while( i .le. this% NumOfEvents(k) .and. &
                      Events(k,i)% plane .eq. SplittingPlane .and. Events(k,i)% eType .eq. START_ )                                    
               pplus = pplus + 1
               i = i + 1
            end do

            N_R(k) = N_R(k) - pminus - pvert
   
            call ComputeSplitSurface_L(this, SplittingPlane, k, S_L)
            call ComputeSplitSurface_R(this, SplittingPlane, k, S_R)

            SplittingCost(1) = computeSplittingCost( this% S, S_L, S_R, N_L(k)+N_P(k), N_R(k) )
            SplittingCost(2) = computeSplittingCost( this% S, S_L, S_R, N_L(k), N_R(k)+N_P(k) )
               
            BestSide = minloc(SplittingCost) 

            if( this% SplitCost .gt. SplittingCost(BestSide(1)) ) then
               this% SplittingPlane = SplittingPlane
               this% SplitCost = SplittingCost(BestSide(1))
               this% split = .true.
               this% axis  = k
               this% SIDE  = BestSide(1) 
            end if
            N_L(k) = N_L(k) + pplus + pvert
            N_P(k) = pvert

         end do
      end do

      end select

   end subroutine KDtree_EvaluateCost
   
   subroutine Event_ClassifyObjs( this, tree, child_L, child_R, ObjsIndx )

      implicit none
      
      type(Event),  dimension(:), intent(in)    :: this
      type(KDtree),               intent(inout) :: tree
      type(KDtree),               intent(inout) :: child_L, child_R
      integer,      dimension(:), intent(inout) :: ObjsIndx
      
      integer :: i, N_B, N_L, N_R

      ObjsIndx = BOTH
      
      N_B = tree% NumOfObjs; N_L = 0; N_R = 0
           
      do i = 1, tree% NumOfEvents(tree% axis)
         if( this(i)% eType .eq. END_ .and. &
            (this(i)% plane .lt. tree% SplittingPlane .or. this(i)% plane .eq. tree% SplittingPlane) ) then         
            ObjsIndx(this(i)% index) = ONLYLEFT   
            N_L = N_L + 1
         elseif( this(i)% eType .eq. START_ .and. &
               (this(i)% plane .gt. tree% SplittingPlane .or. this(i)% plane .eq. tree% SplittingPlane) ) then 
            ObjsIndx(this(i)% index) = ONLYRIGHT 
            N_R = N_R + 1
         elseif( this(i)% eType .eq. PLANAR_ ) then
            if( this(i)% plane .lt. tree% SplittingPlane .or. &
                (this(i)% plane .eq. tree% SplittingPlane .and. tree% SIDE .eq. 1) ) then
               ObjsIndx(this(i)% index) = ONLYLEFT
               N_L = N_L + 1
            elseif( this(i)% plane .gt. tree% SplittingPlane .or. &
                    (this(i)% plane .eq. tree% SplittingPlane .and. tree% SIDE .eq. 2) ) then 
               ObjsIndx(this(i)% index) = ONLYRIGHT
               N_R = N_R + 1
            end if   
         end if
      end do
      
      N_B = N_B - N_L - N_R
      
      tree% N_L = N_L; tree% N_R = N_R; tree% N_B = N_B
      child_L% NumOfObjs = N_B+N_L 
      child_R% NumOfObjs = N_B+N_R 
      
   end subroutine Event_ClassifyObjs   
   
   subroutine KDtree_BuildChild( this, child, side )
   
      implicit none
   
      class(KDtree), intent(inout) :: this
      type(KDtree),  intent(inout) :: child
      integer,       intent(in)    :: side
      
      integer :: i
   
      child% vertices = this% vertices
      child% isLast = .false.

      select case( this% axis ) 
         case( 1 ) 
            do i = 1, 4
               child% vertices(1,vertices_x(side+i)) = this% SplittingPlane
            end do
         case( 2 ) 
            do i = 1, 4
               child% vertices(2,vertices_y(side+i)) = this% SplittingPlane
            end do
         case( 3 ) 
            do i = 1, 4
               child% vertices(3,vertices_z(side+i)) = this% SplittingPlane
            end do
      end select
   
   end subroutine KDtree_BuildChild

   subroutine Event_BuildLists( this, tree, child_L, child_R, ObjsIndx, Events_L, Events_R, parallelization_type )
      use omp_lib
      implicit none
      
      type(Event),  dimension(:,:), intent(in)    :: this
      type(KDtree),                 intent(in)    :: tree
      type(KDtree),                 intent(inout) :: child_L, child_R
      integer,      dimension(:),   intent(inout) :: ObjsIndx
      type(Event),  dimension(:,:), intent(inout) :: Events_L, Events_R
      integer,                      intent(in)    :: parallelization_type

      type(taskPart_type), allocatable :: taskPart(:)
      integer                          :: i, k, B, L, R, j, level, index, taskNum, N_L_start, N_R_start
      integer                          :: index_L, index_R, N_Events_L(NDIM), N_Events_R(NDIM)
      
      select case( parallelization_type )
    
      case( BREADTHFIRST )

      child_L% NumOfEvents = 0
      child_R% NumOfEvents = 0
      call constrct_taskPart( taskPart, tree% NumOfEvents(:), tree% NumThreads )
      
!$omp parallel shared(N_Events_L,N_Events_R,Events_L,Events_R,this,ObjsIndx,taskPart)
!$omp single
      N_Events_L = 0
      N_Events_R = 0
      do taskNum = 1, size(taskPart)
!$omp task private(i,k,L,R,index) firstprivate(taskNum)
         do k = 1, NDIM
            L = 0; R = 0 
            index = taskPart(taskNum)% taskPartDim(k)% starting_index
            do i = 1, taskPart(taskNum)% taskPartDim(k)% ChunkDim
               if( ObjsIndx(this(k,index+i)% index) .eq. ONLYLEFT ) then 
                  L = L + 1
               elseif( ObjsIndx(this(k,index+i)% index) .eq. ONLYRIGHT ) then
                  R = R + 1
               elseif( ObjsIndx(this(k,index+i)% index) .eq. BOTH ) then
                  L = L + 1
                  R = R + 1
               end if
            end do 
            taskPart(taskNum)% taskPartDim(k)% N_L = L
            taskPart(taskNum)% taskPartDim(k)% N_R = R
         end do
!$omp end task
      end do
!$omp taskwait 

      do k = 1, NDIM
         N_Events_L = sum(taskPart(:)% taskPartDim(k)% N_L)  
         N_Events_R = sum(taskPart(:)% taskPartDim(k)% N_R)
      end do

      do taskNum = size(taskPart),2,-1
        do k = 1, NDIM
           taskPart(taskNum)% taskPartDim(k)% N_R = sum(taskPart(1:taskNum-1)% taskPartDim(k)% N_R)
           taskPart(taskNum)% taskPartDim(k)% N_L = sum(taskPart(1:taskNum-1)% taskPartDim(k)% N_L)
        end do
     end do
     
     taskPart(1)% taskPartDim(:)% N_R = 0
     taskPart(1)% taskPartDim(:)% N_L = 0

     do taskNum = 1, size(taskPart)
!$omp task private(i,k,L,R,index_L,index_R,index) firstprivate(taskNum)
         do k = 1, NDIM
            L = 0; R = 0 
            index_L = taskPart(taskNum)% taskPartDim(k)% N_L
            index_R = taskPart(taskNum)% taskPartDim(k)% N_R
            index = taskPart(taskNum)% taskPartDim(k)% starting_index
            do i = 1, taskPart(taskNum)% taskPartDim(k)% ChunkDim
               if( ObjsIndx(this(k,index+i)% index) .eq. ONLYLEFT ) then 
                  L = L + 1
                  Events_L(k,L+index_L) = this(k,index+i)
               elseif( ObjsIndx(this(k,index+i)% index) .eq. ONLYRIGHT ) then
                  R = R + 1
                  Events_R(k,R+index_R) = this(k,index+i)
               elseif( ObjsIndx(this(k,index+i)% index) .eq. BOTH ) then
                  L = L + 1
                  R = R + 1
                  Events_L(k,L+index_L) = this(k,index+i)
                  Events_R(k,R+index_R) = this(k,index+i)
               end if
            end do 
         end do
!$omp end task
      end do
!$omp end single
!$omp end parallel    
      
      child_L% NumOfEvents = N_Events_L
      child_R% NumOfEvents = N_Events_R

      deallocate(taskPart)
      
      case( DEPTHFIRST )
    
      child_L% NumOfEvents = 0
      child_R% NumOfEvents = 0
      
      do k = 1, NDIM
         L = 0; R = 0
         do i = 1, tree% NumOfEvents(k)
            if( ObjsIndx(this(k,i)% index) .eq. ONLYLEFT ) then 
               L = L + 1
               Events_L(k,L) = this(k,i)
               child_L% NumOfEvents(k) = child_L% NumOfEvents(k) + 1
            elseif( ObjsIndx(this(k,i)% index) .eq. ONLYRIGHT ) then
               R = R + 1
               Events_R(k,R) = this(k,i)
               child_R% NumOfEvents(k) = child_R% NumOfEvents(k) + 1
            elseif( ObjsIndx(this(k,i)% index) .eq. BOTH ) then
               L = L + 1
               R = R + 1
               Events_L(k,L) = this(k,i)
               Events_R(k,R) = this(k,i)
               child_L% NumOfEvents(k) = child_L% NumOfEvents(k) + 1
               child_R% NumOfEvents(k) = child_R% NumOfEvents(k) + 1     
            end if
         end do 
      end do
      
      end select

   end subroutine Event_BuildLists
   
   subroutine KDtree_SaveObjsIndeces( this, Events )
   
      implicit none
      
      class(KDtree), intent(inout) :: this
      type(Event),   intent(in)    :: Events(:,:)
      integer                      :: i, k
      
      allocate(this% ObjsIndeces(this% NumOfObjs))
      
      k = 1
      
      do i = 1, this% NumOfEvents(1)
         if( Events(1,i)% eType .eq. START_ .or. &
             Events(1,i)% eType .eq. PLANAR_      ) then
            this% ObjsIndeces(k) = Events(1,i)% index
            k = k + 1
            if( k .gt. this% NumOfObjs ) exit
         end if
      end do
   
   end subroutine KDtree_SaveObjsIndeces
   
   subroutine Event_Construct( this, ObjectsList, tree, axis )
    
      implicit none
      
      type(Event),       dimension(:), intent(inout) :: this
      type(object_type), dimension(:), intent(in)    :: ObjectsList
      type(KDtree),                    intent(inout) :: tree
      integer,                         intent(in)    :: axis
      
      integer :: i, j, k, shift

      tree% NumOfEvents(axis) = 0

      do i =1, tree% NumOfObjs
         this(i)% index = ObjectsList(i)% index 
         this(i+tree% NumOfObjs)% index = ObjectsList(i)% index 
     
         this(i)% plane = minval(ObjectsList(i)% vertices(:)% coords(axis)) 
         this(i+tree% NumOfObjs)% plane = maxval(ObjectsList(i)% vertices(:)% coords(axis)) 

         if( almostEqual(this(i)% plane,this(i+tree% NumOfObjs)% plane) ) then
            this(i)% eType = PLANAR_
            !pushed to the end
            this(i+tree% NumOfObjs)% plane = huge(1.0_RP)
            this(i+tree% NumOfObjs)% eType = -1
            tree% NumOfEvents(axis) = tree% NumOfEvents(axis) + 2           
         else
            this(i)% eType = START_
            this(i+tree% NumOfObjs)% eType = END_
            tree% NumOfEvents(axis) = tree% NumOfEvents(axis) + 2
         end if
      end do
      
   end subroutine Event_Construct
  
   
!////////////////// SUBROUTINES FOR KD TREE made of points ///////////////////////
   
   recursive subroutine KDtree_buildPoints_BreadthFirst( this, Events, ObjsIndx, level, Depth_First )
   
      implicit none
      
      type(KDtree), target,       intent(inout) :: this
      type(Event),  allocatable,  intent(inout) :: Events(:,:)
      integer,                    intent(in)    :: level
      integer,      dimension(:), intent(inout) :: ObjsIndx 
      type(DepthFirst_type),      intent(inout) :: Depth_First(:)
      
      type(KDtree),   pointer :: child_L, child_R 
      type(Event), dimension(:,:), allocatable :: Events_L, Events_R
      integer :: j

      BoxIndex = BoxIndex + 1

      this% index  = BoxIndex
      this% NumOfEvents(:)  = this% NumOfObjs

      if( this% level .lt. level .and. this% NumOfObjs .gt. 0 ) then

         if( mod(this% NumOfObjs,2) .eq. 0 ) then
            this% HalfEvents = this% NumOfObjs/2
            this% SIDE = EVEN
         else
            this% HalfEvents = this% NumOfObjs/2
            this% HalfEvents = this% HalfEvents+1
            this% SIDE = ODD
         end if

         this% SplittingPlane = Events(this% axis,this% HalfEvents)% plane
         
         allocate(this% child_L,this% child_R)
         allocate(this% child_L% parent, this% child_R% parent)
         
         this% child_L% parent => this; this% child_R% parent => this
        
         child_L => this% child_L
         child_R => this% child_R
           
         child_L% level = this% level + 1
         child_R% level = this% level + 1

         this% isLast = .false.   

         call this% BuildChild( child_L, side_L )
         call this% BuildChild( child_R, side_R )
         
         child_L% Min_n_of_Objs = this% Min_n_of_Objs
         child_R% Min_n_of_Objs = this% Min_n_of_Objs
         child_L% NumThreads = this% NumThreads
         child_R% NumThreads = this% NumThreads

         child_L% axis = AxisIndex( this% axis )
         child_R% axis = AxisIndex( this% axis )

         call Points_Event_ClassifyObjs( Events(this% axis,:), this, child_L, child_R, ObjsIndx )

         allocate(Events_L(NDIM,child_L% NumOfObjs))
         allocate(Events_R(NDIM,child_R% NumOfObjs))

         call Event_BuildLists( Events, this, child_L, child_R, ObjsIndx, Events_L, Events_R, BREADTHFIRST )

         deallocate(Events)
            
         call KDtree_buildPoints_BreadthFirst( child_L, Events_L, ObjsIndx, level, Depth_First )
         call KDtree_buildPoints_BreadthFirst( child_R, Events_R, ObjsIndx, level, Depth_First )

      else
         this% isLast = .false.
         if( this% NumOfObjs .eq. 0  ) this% isLast = .true.
         do j = 1, size(Depth_First)
            if( .not. Depth_First(j)% active ) then 
               Depth_First(j)% taskTree => this
               allocate(Depth_First(j)% taskEvents(NDIM,this% NumOfObjs))
               Depth_First(j)% taskEvents = Events
               Depth_First(j)% active = .true.
               deallocate(Events)
               exit
            end if
         end do

      end if
   
   end subroutine KDtree_buildPoints_BreadthFirst
   
   recursive subroutine KDtree_buildPoints_DepthFirst( this, Events, ObjsIndx )
   use omp_lib
      implicit none
      
      type(KDtree), target,       intent(inout) :: this
      type(Event),  allocatable,  intent(inout) :: Events(:,:)
      integer,      dimension(:), intent(inout) :: ObjsIndx 
      
      type(KDtree),   pointer :: child_L, child_R 
      type(Event), dimension(:,:), allocatable :: Events_L, Events_R

!$omp critical
      BoxIndex = BoxIndex + 1
!$omp end critical

      this% index  = BoxIndex
      this% NumOfEvents(:)  = this% NumOfObjs
      
      if( this% level .lt. depth .and. this% NumOfObjs .ge. this% Min_n_of_Objs ) then

         if( mod(this% NumOfObjs,2) .eq. 0 ) then
            this% HalfEvents = this% NumOfObjs/2
            this% SIDE = EVEN
         else
            this% HalfEvents = this% NumOfObjs/2
            this% HalfEvents = this% HalfEvents+1
            this% SIDE = ODD
         end if

         this% SplittingPlane = Events(this% axis,this% HalfEvents)% plane
         
         allocate(this% child_L,this% child_R)
         allocate(this% child_L% parent, this% child_R% parent)
         this% child_L% parent => this; this% child_R% parent => this
        
         child_L => this% child_L
         child_R => this% child_R
           
         child_L% level = this% level + 1
         child_R% level = this% level + 1

         this% isLast = .false.   

         call this% BuildChild( child_L, side_L )
         call this% BuildChild( child_R, side_R )
         
         child_L% Min_n_of_Objs = this% Min_n_of_Objs
         child_R% Min_n_of_Objs = this% Min_n_of_Objs
         child_L% NumThreads = this% NumThreads
         child_R% NumThreads = this% NumThreads
         
         child_L% axis = AxisIndex( this% axis )
         child_R% axis = AxisIndex( this% axis )

         call Points_Event_ClassifyObjs( Events(this% axis,:), this, child_L, child_R, ObjsIndx )

         allocate(Events_L(NDIM,child_L% NumOfObjs))
         allocate(Events_R(NDIM,child_R% NumOfObjs))

         call Event_BuildLists( Events, this, child_L, child_R, ObjsIndx, Events_L, Events_R, DEPTHFIRST )

         deallocate(Events)

         call KDtree_buildPoints_DepthFirst( child_L, Events_L, ObjsIndx )
         call KDtree_buildPoints_DepthFirst( child_R, Events_R, ObjsIndx )
            
      else
         this% isLast = .true.
         if( this% NumOfObjs .gt. 0 ) then
            call this% SavePointsIndeces(Events)     
            deallocate(Events)
         end if
      end if

   end subroutine KDtree_buildPoints_DepthFirst
   
   
   subroutine KDtree_SavePointsIndeces( this, Events )
   
      implicit none
      
      class(KDtree), intent(inout) :: this
      type(Event),   intent(in)    :: Events(:,:)
      integer                      :: i, k
      
      allocate(this% ObjsIndeces(this% NumOfObjs))
   
      do i = 1, this% NumOfObjs
         this% ObjsIndeces(i) = Events(1,i)% index
      end do   
   
   end subroutine KDtree_SavePointsIndeces
   
   subroutine Points_Event_ClassifyObjs( this, tree, child_L, child_R, ObjsIndx )
    
    use omp_lib
      implicit none
      
      type(Event),  dimension(:), intent(in)    :: this
      type(KDtree),               intent(inout) :: tree
      type(KDtree),               intent(inout) :: child_L, child_R
      integer,      dimension(:), intent(inout) :: ObjsIndx
      
      integer :: i, N_B, N_L, N_R
      
      N_L = 0; N_R = 0; N_B = 0
      
      do i = 1, tree% NumOfObjs
         if( i .lt. tree% HalfEvents ) then        
             ObjsIndx(this(i)% index) = ONLYLEFT   
             N_L = N_L + 1
         elseif( i .gt. tree% HalfEvents ) then 
            ObjsIndx(this(i)% index) = ONLYRIGHT 
            N_R = N_R + 1   
         elseif( i .eq. tree% HalfEvents ) then
            if( tree% SIDE .eq. ODD ) then
               ObjsIndx(this(i)% index) = BOTH 
               N_B = N_B + 1 
            else
               ObjsIndx(this(i)% index) = ONLYLEFT 
               N_L = N_L + 1
            end if
         end if
      end do
            
      tree% N_L = N_L; tree% N_R = N_R; tree% N_B = N_B
      child_L% NumOfObjs = N_L + N_B
      child_R% NumOfObjs = N_R + N_B
      
   end subroutine Points_Event_ClassifyObjs   
   
   subroutine Points_Event_construct(this, tree, axis, PointList)
      use MPI_IBMUtilities
      implicit none
      
      type(Event),      dimension(:), intent(inout) :: this
      type(KDtree),                   intent(inout) :: tree
      integer,                        intent(in)    :: axis
      type(point_type), dimension(:), intent(in)    :: PointList   
      
      integer :: i
      
      tree% NumOfEvents = tree% NumOfObjs
      
      do i = 1, tree% NumOfObjs !no planar evetns thus NumOfObjs = NumOfEvents
         this(i)% index = PointList(i)% index
         this(i)% plane = PointList(i)% coords(axis)
         this(i)% eType = START_
      end  do
   
   end subroutine Points_Event_construct
   
!/////////////////////////////// EVENTS SORTING ///////////////////////

  recursive subroutine QsortEvents( Events, startIdx, endIdx )
  
     implicit none
     !-arguments-------------------------------------------------
     type(Event), dimension(:),   intent(inout) :: Events
     integer,                       intent(in)  :: startIdx, endIdx
     
     type(Event) :: pivot
     integer :: left, right, mid
     
     if( startIdx .gt. endIdx ) return
     
     left = startIdx; right = endIdx
     
     mid = (startIdx + endIdx)/2
     
     pivot = Events(mid)
     
     do while ( left < right )
        do while( EventIsLess(Events(left),pivot) )
           left = left + 1
        end do
        do while( EventIsLess(pivot,Events(right)) )
           right = right - 1
        end do
        if( left .le. right ) then
           call swapEvents(Events(left), Events(right))
           left = left + 1
           right = right - 1
        end if
     end do
     
     if( startIdx .lt. right ) call QsortEvents( Events, startIdx, right )
     if( left .lt. endIdx  ) call QsortEvents( Events, left, endIdx )
 
  end subroutine QsortEvents

  subroutine swapEvents( EventA, EventB )
  
     implicit none
  
     type(Event), intent(inout) :: EventA, EventB
  
     type(Event) :: temp
  
     temp = EventA
     EventA = EventB
     EventB = temp
  
  end subroutine swapEvents
  
  logical function EventIsLess( EventA, EventB ) result( IsLess )
  
     implicit none
  
     type(Event), intent(in) :: EventA, EventB
     
     if( EventA% plane .eq. EventB% plane ) then
        if( EventA% eTYPE .lt. EventB% eTYPE ) then
           IsLess = .true.
        else
           IsLess = .false.
        end if
     elseif( EventA% plane < EventB% plane ) then
        IsLess = .true.
     else
        IsLess = .false.
     end if
  
  end function EventIsLess
   
   subroutine constrct_taskPart( taskPart, TotalChunk, NumOfThreads )
   
      implicit none
      
      type(taskPart_type), allocatable, intent(inout) :: taskPart(:)
      integer,                          intent(in)    :: NumOfThreads
      integer,      dimension(NDIM),    intent(in)    :: TotalChunk
      
      integer, dimension(NDIM,NumOfThreads) :: ChunkDim
      integer :: i, k
   
   
      ChunkDim = ChunkPartition( TotalChunk, NumOfThreads )
      
      if( allocated(taskPart) ) deallocate(taskPart)
      allocate(taskPart(NumOfThreads))
      
      do k = 1, NDIM
         do i = 1, NumOfThreads
            taskPart(i)% taskPartDim(k)% ChunkDim = ChunkDim(k,i)
            taskPart(i)% taskPartDim(k)% Starting_index = (i-1)*ChunkDim(k,1)
         
            taskPart(i)% taskPartDim(k)% N_R = 0
            taskPart(i)% taskPartDim(k)% N_L = 0
            taskPart(i)% taskPartDim(k)% N_P = 0
            taskPart(i)% taskPartDim(k)% Values2Check = 0
         end do
      end do
   
   end subroutine constrct_taskPart
   
   function ChunkPartition( TotalChunk, NumOfThreads ) result(ChunkDim)
   
      implicit none
   
      integer,                    intent(in) :: NumOfThreads
      integer,  dimension(NDIM),  intent(in) :: TotalChunk
      integer,  dimension(NDIM,NumOfThreads) :: ChunkDim
      
      integer :: Part, i, k
   
      do k = 1, NDIM
         Part = TotalChunk(k)/NumOfThreads
      
         do i = 1, NumOfThreads
            ChunkDim(k,i) = Part
            if( i .eq. NumOfThreads ) then
               ChunkDim(k,i) = ChunkDim(k,i) + (TotalChunk(k) - i*Part) 
            end if
         end do
      end do
   
   end function ChunkPartition

end module KDClass
