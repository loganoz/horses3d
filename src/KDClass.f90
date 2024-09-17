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
       
  public :: LastLevel, BoxIndex, POINTS_KDTREE, TRIANGLES_KDTREE_SAH, TRIANGLES_KDTREE_MEDIAN

  integer :: BoxIndex, LastLevel = -1, depth, NumOfObjsIndx, NumOfKDtreePoints
  
  integer, parameter :: side_L = 4, side_R = 0, ON_PLANE = 0, POINTS_KDTREE = 0, TRIANGLES_KDTREE_SAH = 1, TRIANGLES_KDTREE_MEDIAN = 2
  integer, parameter :: START_ = 2, END_ = 0, PLANAR_ = 1, ONLYLEFT = 0, ONLYRIGHT = 1, BOTH = 2, ODD = 1, EVEN = 2
  integer, parameter :: BREADTHFIRST = 1, DEPTHFIRST = 0
  
  real(kind=rp), parameter :: C_TRANSVERSE = 1.0_RP, C_INTERSECT = 1.5_RP, C_1 = 1.2_RP, C_2 = 2.0_RP
  
  integer,  parameter, dimension(8) :: vertices_x = (/ 1,4,5,8,  2,3,6,7 /), &
                                       vertices_y = (/ 1,2,6,5,  3,4,7,8 /), &
                                       vertices_z = (/ 1,2,3,4,  5,6,7,8 /)  
                                       
  type Event
  
     real(kind=rp) :: plane, median
     integer       :: eType, index, axis
  
  end type
                            
!
!  **************************************************
!  Main type for a kd-tree
!  **************************************************
   type KDtree

      class(KDtree), pointer                         :: child_L, child_R, parent
      type(Object_type), dimension(:),   allocatable :: ObjectsList 
      real(kind=rp),     dimension(3,8)              :: vertices   !local
      integer                                        :: NumOfObjs, level, axis, &
                                                        index, Min_n_of_Objs,   &
                                                        which_KDtree, MaxAxis,  &
                                                        SIDE, N_L, N_R, N_B,    &
                                                        HalfEvents, NumThreads, &
                                                        STLNum
      logical                                        :: isLast, Split
      logical                                        :: built_R = .false., built_L = .false.
      integer,           dimension(NDIM)             :: NumOfEvents
      real(kind=rp)                                  :: S, SplitCost, SplittingPlane
      type(Event),       dimension(:,:), allocatable :: Events
      integer,           dimension(:),   allocatable :: ObjsIndeces
      
      contains
         procedure :: construct           => KDtree_construct
         procedure :: SetUpRoot           => KDtree_SetUpRoot
         procedure :: FindLeaf            => KDtree_FindLeaf 
         procedure :: plot                => KDtree_plot
         procedure :: plotBlock           => KDtree_plotBlock
         procedure :: Destruct            => KD_treeDestruct
         procedure :: GetArea             => KD_treeGetArea
         procedure :: BuildChild          => KDtree_BuildChild
         procedure :: EvaluateCostSAH     => KDtree_EvaluateCostSAH
         procedure :: EvaluateCostMEDIAN  => KDtree_EvaluateCostMEDIAN
         procedure :: SaveObjsIndeces     => KDtree_SaveObjsIndeces
         procedure :: SavePointsIndeces   => KDtree_SavePointsIndeces

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
      !-arguments----------------------
      type(KDtree), intent(in) :: this
      !-local-variables----------------
      integer,  dimension(1) :: maxvec
      
      maxvec = maxloc( (/ abs(this% vertices(1,7) - this% vertices(1,1)),   &
                          abs(this% vertices(2,7) - this% vertices(2,1)),   &
                          abs(this% vertices(3,7) - this% vertices(3,1)) /) )
       
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
      !-arguments--------------------------------------
      type(KDtree), intent(inout) :: tree
      integer,      intent(in)    :: STLNum, funit, &
                                    which_KDtree
      !-local-variables--------------------------------
      real(kind=rp) :: x_g(NDIM)
      integer       :: i
      
      if( tree% isLast ) then 
         write(funit,"(a69)") 'ZONE NODES=8, ELEMENTS = 6, DATAPACKING=POINT, ZONETYPE=FETETRAHEDRON'

         do i = 1, 8
            call OBB(STLNum)% ChangeRefFrame( tree% vertices(:,i), GLOBAL, x_g )
            write(funit,'(3E13.5)') Lref*x_g(1),Lref*x_g(2), Lref*x_g(3)
         end do 
         
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
   subroutine KDtree_plot( this, lvl )
      use MPI_Process_Info
      implicit none
      !-arguemnts----------------------------------------
      class(KDtree), target,   intent(inout) :: this
      integer,       optional, intent(in)    :: lvl
      !-local-variables----------------------------------
      real(kind=rp)              :: x_g(NDIM)
      character(len=LINE_LENGTH) :: filename, myString
      integer                    :: i, funit

      select case( this% which_KDtree )
         case( POINTS_KDTREE )
          
            if( .not. MPI_Process% isRoot ) return
                        
            if( present(lvl) ) then
               if( lvl .gt. 0 ) then
                  write(myString,'(i100)') lvl
                  filename = 'KDtreeBandPoints_MGlevel'//trim(adjustl(myString))
               else 
                  filename = 'KDTreeBandPoints'
               end if
            else
               filename = 'KDTreeBandPoints'
            end if

         case( TRIANGLES_KDTREE_SAH )
         
            write(myString,'(i100)') MPI_Process% rank
      
            if( MPI_Process% nProcs .eq. 1 ) then
               filename = 'KDTree_'//OBB(this% STLNum)% filename
            else
               filename = 'KDTree_Partition'//trim(adjustl(myString))//'_'//OBB(this% STLNum)% filename
               filename = trim(filename)
            end if
            
      end select
          
      funit = UnusedUnit()
     
      open(funit,file='IBM/'//trim(filename)//'.tec', status='unknown')
         
      write(funit,"(a28)") 'TITLE = "KD-tree"'
      write(funit,"(a25)") 'VARIABLES = "x", "y", "z"'

      call Plot_leaves( this, this% STLNum, funit, this% which_KDtree )

      close(funit)
   
   end subroutine KDtree_plot
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine plot a single block whose index is BlockIndex. if PlotObjs is true
! also the objects belonging to the box are saved and plot 
!  ------------------------------------------------
   subroutine KDtree_plotBlock( this, ObjectsList, STLNum, PlotObjs )
   
      implicit none
      !-arguments---------------------------------------------
      class(KDtree),     intent(inout) :: this
      type(object_type), intent(in)    :: ObjectsList(:)
      integer,           intent(in)    :: STLNum
      logical,           intent(in)    :: PlotObjs
      !-local-variables---------------------------------------
      real(kind=rp)              :: x_g(NDIM)
      character(len=LINE_LENGTH) :: filename, myString
      integer                    :: i,j, funit, index
      
      funit = UnusedUnit()

      write(myString,'(i100)') this% index
      filename = 'block'//trim(adjustl(myString))//'_'//OBB(STLNum)% filename
      filename = trim(filename)
      
      open(funit,file='IBM/'//trim(filename)//'.tec', status='unknown')
 
      write(funit,"(a28)") 'TITLE = "BLOCK"'
      write(funit,"(a25)") 'VARIABLES = "x", "y", "z"'
      
      write(funit,"(a69)") 'ZONE NODES=8, ELEMENTS = 6, DATAPACKING=POINT, ZONETYPE=FETETRAHEDRON'
      do i = 1, 8
         call OBB(STLNum)% ChangeRefFrame( this% vertices(:,i), GLOBAL, x_g )
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

        filename = 'Objects'//trim(adjustl(myString))//'_'//OBB(STLNum)% filename
        filename = trim(filename)
      
         open(funit,file='IBM/'//trim(filename)//'.tec', status='unknown')
 
         write(funit,"(a28)") 'TITLE = "ObjsIndx"'
         write(funit,"(a25)") 'VARIABLES = "x", "y", "z"'

         do i = 1, this% NumOfObjs
            write(funit,"(a66)") 'ZONE NODES=3, ELEMENTS = 1, DATAPACKING=POINT, ZONETYPE=FETRIANGLE'
            index = this% ObjsIndeces(i)
            do j = 1, 3
               call OBB(STLNum)% ChangeRefFrame( ObjectsList(index)% vertices(j)% coords, GLOBAL, x_g )
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
      class(KDtree),              intent(inout) :: this
      type(STLfile),              intent(in)    :: stl
      real(kind=RP),              intent(in)    :: Vertices(:,:)
      type(Point_type), optional, intent(in)    :: PointList(:)
      !-local-variables--------------------------------------
      type(Object_type) :: Obj 
      integer           :: i
      
      this% NumOfObjs = 0
      
      do i = 1, 4
         this% vertices(:,i)   = Vertices(:,i)
         this% vertices(:,i+4) = Vertices(:,i+4)
      end do
      
      select case( this% which_KDtree )
      
         case( POINTS_KDTREE )
         
            allocate( this% ObjectsList(size(PointList)) )
            this% NumOfObjs = size(PointList)
         
         case( TRIANGLES_KDTREE_MEDIAN,TRIANGLES_KDTREE_SAH )
         
            allocate( this% ObjectsList(stl% NumOfObjs) )
         
            associate( Objs => stl% ObjectsList )

            do i = 1, stl% NumOfObjs
               this% ObjectsList(i) = Objs(i)            
            end do

            end associate

            this% NumOfObjs = stl% NumOfObjs
         
      end select
   
   end subroutine KDtree_SetUpRoot
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine builds the KD tree. 
!  ------------------------------------------------   
   subroutine KDtree_construct( this, stl, Vertices, isPlot, Min_n_of_Objs, PointList, lvl )
      use omp_lib
      implicit none
      !-arguments-------------------------------------------------
      class(KDtree),              intent(inout) :: this
      type(STLfile),              intent(in)    :: stl
      real(kind=RP),              intent(in)    :: Vertices(:,:)
      logical,                    intent(in)    :: isPlot
      integer,                    intent(in)    :: Min_n_of_Objs
      integer,          optional, intent(in)    :: lvl
      type(point_type), optional, intent(in)    :: PointList(:)
      !-local-varables--------------------------------------------
      real(kind=rp)                      :: NumOfObjs,        &
                                            NumThreads_RP
      integer,               allocatable :: ObjsIndx(:)
      type(Event),           allocatable :: Events(:,:)
      type(DepthFirst_type), allocatable :: Depth_First(:)
      integer                            :: i, k, NumThreads, &
                                            DepthFirstLevel,  &
                                            NumDepthFirst
      
      if( this% which_KDtree .eq. POINTS_KDTREE ) then
         if( .not. present(PointList) ) then
            print *, "KDtree_construct:: PointList needed for this type of KDtree"
            error stop
         end if
         call this% SetUpRoot( stl, vertices, PointList )   
      else
         call this% SetUpRoot( stl, vertices )
      end if
      
      this% level         = 0
      NumOfObjs           = this% NumOfObjs
      this% Min_n_of_Objs = Min_n_of_Objs
      this% isLast        = .false.
      BoxIndex            = 1
      this% index         = BoxIndex
      
      select case( this% which_KDtree )
         case( POINTS_KDTREE )

            this% axis = FindAxis( this )

            allocate(Events(NDIM,this% NumOfObjs))

            call GetPointsEvents( this, Events, PointList, NumThreads )
      
            depth             = huge(1)
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

         case( TRIANGLES_KDTREE_SAH )

            allocate(Events(NDIM,2*this% NumOfObjs))      

            call GetEvents( this, Events, NumThreads )
          
            depth             = C_1*log(NumOfObjs)/log(2.0_RP) + C_2!1.2 * log2(N) + 2 --> on improving kd tree for ray shooting Havran
            this% NumThreads  = NumThreads
            NumThreads_RP     = NumThreads 
            DepthFirstLevel   = floor(log(NumThreads_RP)/log(2.0_RP)) 
            NumDepthFirst     = 2.0_RP**(DepthFirstLevel)
            
            allocate(Depth_First(NumThreads),ObjsIndx(this% NumOfObjs))

            call KDtree_buildTRIANGLES_BreadthFirst( this, Events, ObjsIndx, DepthFirstLevel, Depth_First )

            deallocate(ObjsIndx)
!$omp parallel shared(this,Depth_First,k) 
!$omp single 
            do k = 1, count( Depth_First(:)% active )
!$omp task firstprivate(k) private(ObjsIndx) 
               allocate(ObjsIndx(this% NumOfObjs))
               call KDtree_buildTRIANGLES_DepthFirst( Depth_First(k)% taskTree, Depth_First(k)% taskEvents, ObjsIndx )
               nullify( Depth_First(k)% taskTree )
               deallocate(ObjsIndx)
!$omp end task
            end do
!$omp end single
!$omp end parallel
            deallocate(Depth_First)

         case( TRIANGLES_KDTREE_MEDIAN )

            this% axis = FindAxis(this)

            allocate(Events(NDIM,2*this% NumOfObjs))      

            call GetEvents( this, Events, NumThreads )
          
            depth             = log(NumOfObjs)/log(2.0_RP) !--> on improving kd tree for ray shooting Havran
            this% NumThreads  = NumThreads
            NumThreads_RP     = NumThreads 
            DepthFirstLevel   = floor(log(NumThreads_RP)/log(2.0_RP)) 
            NumDepthFirst     = 2.0_RP**(DepthFirstLevel)
            
            allocate(Depth_First(NumThreads),ObjsIndx(this% NumOfObjs))

            call KDtree_buildTRIANGLES_BreadthFirst( this, Events, ObjsIndx, DepthFirstLevel, Depth_First )

            deallocate(ObjsIndx)

!$omp parallel shared(this,Depth_First,k) 
!$omp single 
            do k = 1, count( Depth_First(:)% active )
!$omp task firstprivate(k) private(ObjsIndx) 
               allocate(ObjsIndx(this% NumOfObjs))
               call KDtree_buildTRIANGLES_DepthFirst( Depth_First(k)% taskTree, Depth_First(k)% taskEvents, ObjsIndx )
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

      if( present(lvl) ) then
         if( isPlot ) call this% plot(lvl)
      else
         if( isPlot ) call this% plot()
      end if

   end subroutine KDtree_construct
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine finds the leaf where Point lies
!  ------------------------------------------------   
   subroutine KDtree_FindLeaf( this, Point, tree, RIGHTCHILD ) 
   
      implicit none
      !-arguments----------------------------------------------
      class(KDtree), target, intent(inout) :: this
      real(kind=rp),         intent(in)    :: Point(:)
      type(KDtree), pointer, intent(inout) :: tree
      logical,               intent(in)    :: RIGHTCHILD
      !-local-variables----------------------------------------
      real(kind=rp) :: vertices(NDIM)
      integer       :: level
   
      tree => this
      
      if( tree% isLast ) return
      
      level = 0
      
      do
         if( RIGHTCHILD ) then
            if(  Point(tree% axis) .lt. tree% SplittingPlane ) then !lt
               tree => tree% child_L            
            elseif( Point(tree% axis) .ge. tree% SplittingPlane ) then
               tree => tree% child_R
            end if
         else
            if(  Point(tree% axis) .le. tree% SplittingPlane ) then !lt
               tree => tree% child_L            
            elseif( Point(tree% axis) .gt. tree% SplittingPlane ) then
               tree => tree% child_R
            end if      
         endif       
         if( tree% isLast ) return
      end do
      
   end subroutine KDtree_FindLeaf

   recursive subroutine KD_treeDestruct( this, isChild ) 
   
      implicit none
      !-arguments-----------------------------
      class(KDtree), intent(inout) :: this
      logical,       intent(in)    :: isChild
      !-local-variables-----------------------
      integer :: i, stat       
         
      if( allocated(this% ObjectsList) ) then
         if( this% which_KDtree .ne.  POINTS_KDTREE ) then 
            do i = 1, this% NumOfObjs
               deallocate(this% ObjectsList(i)% vertices)
            end do
         end if
         deallocate(this% ObjectsList)      
      end if

      if( allocated(this% ObjsIndeces) ) deallocate(this% ObjsIndeces)
      if( allocated(this% Events) )      deallocate(this% Events)
      
      if( this% built_L .and. .not. this% isLast ) call this% child_L% destruct( isChild )     
      if( this% built_R .and. .not. this% isLast ) call this% child_R% destruct( isChild ) 

      if( this% built_L ) deallocate(this% child_L) 
      if( this% built_R ) deallocate(this% child_R) 
      nullify( this% parent )

   end subroutine KD_treeDestruct
   
   subroutine KD_treeGetArea( tree )
   
      implicit none
      !-arguments---------------------------
      class(KDtree), intent(inout) :: tree
      !-local-variables---------------------
      real(kind=rp) :: L_x, L_y, L_z
      
      L_x = abs(tree% vertices(1,7)-tree% vertices(1,1))
      L_y = abs(tree% vertices(2,7)-tree% vertices(2,1))
      L_z = abs(tree% vertices(3,7)-tree% vertices(3,1))
      
      tree% S = 2.0_RP*(L_x*L_y + L_x*L_z + L_y*L_z)
      
   end subroutine KD_treeGetArea
   
   subroutine GetEvents( tree, Events, NumThreads )
      use omp_lib
      implicit none
      !-arguments---------------------------------
      type(KDtree), intent(inout) :: tree
      type(Event),  intent(inout) :: Events(:,:)
      integer,      intent(out)   :: NumThreads
      !-local-variables---------------------------
      integer :: k
!$omp parallel shared(tree,Events,k)
!$omp single
#ifdef _OPENMP
            NumThreads = omp_get_num_threads()
#else
            NumThreads = 1
#endif
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
      !-arguments-------------------------------------
      type(KDtree),     intent(inout) :: tree
      type(Event),      intent(inout) :: Events(:,:)
      type(point_type), intent(in)    :: PointList(:)
      integer,          intent(out)   :: NumThreads
      !-local-variables-------------------------------
      integer :: k
!$omp parallel shared(tree,Events,k,PointList)
!$omp single
#ifdef _OPENMP
           NumThreads = omp_get_num_threads()
#else
           NumThreads = 1
#endif
            do k = 1, NDIM
!$omp task firstprivate(k)
               call Points_Event_Construct( Events(k,:), tree, k, PointList )
               call QsortEvents( Events(k,:), 1, tree% NumOfEvents(k) )
!$omp end task
            end do        
!$omp end single
!$omp end parallel
   end subroutine GetPointsEvents
   
   recursive subroutine KDtree_buildTRIANGLES_BreadthFirst( this, Events, ObjsIndx, level, Depth_First )
   
      implicit none
      !-arguments------------------------------------------------
      type(KDtree), target,      intent(inout) :: this
      type(Event),  allocatable, intent(inout) :: Events(:,:)
      integer,                   intent(in)    :: level
      integer,                   intent(inout) :: ObjsIndx(:) 
      type(DepthFirst_type),     intent(inout) :: Depth_First(:)
      !-local-variables------------------------------------------
      type(KDtree), pointer     :: child_L, child_R 
      type(Event),  allocatable :: Events_L(:,:), Events_R(:,:)
      integer                   :: j

      call this% GetArea()

      BoxIndex = BoxIndex + 1

      this% index  = BoxIndex
 
      if( this% level .lt. level .and. this% NumOfObjs .gt. 0 ) then
 
         if( this% which_KDtree .eq. TRIANGLES_KDTREE_SAH ) then
            call this% EvaluateCostSAH( Events )
         else
            call this% EvaluateCostMEDIAN( Events )
         end if

         if( this% split  ) then 

            allocate(this% child_L,this% child_R)
                    
            child_L => this% child_L
            child_R => this% child_R
            
            child_L% parent => this
            child_R% parent => this
           
            child_L% level = this% level + 1
            child_R% level = this% level + 1
            
            this% built_L = .true.
            this% built_R = .true.

            this% isLast = .false.   

            if( this% which_KDtree .eq. TRIANGLES_KDTREE_MEDIAN ) then
               child_L% axis = AxisIndex( this% axis )
               child_R% axis = AxisIndex( this% axis )
            endif

            call this% BuildChild( child_L, side_L )
            call this% BuildChild( child_R, side_R )
         
            child_L% Min_n_of_Objs = this% Min_n_of_Objs
            child_R% Min_n_of_Objs = this% Min_n_of_Objs
            child_L% NumThreads    = this% NumThreads
            child_R% NumThreads    = this% NumThreads
            child_L% which_KDtree  = this% which_KDtree
            child_R% which_KDtree  = this% which_KDtree

            call Event_ClassifyObjs( Events(this% axis,:), this, child_L, child_R, ObjsIndx )

            allocate(Events_L(NDIM,2*child_L% NumOfObjs))
            allocate(Events_R(NDIM,2*child_R% NumOfObjs))

            call Event_BuildLists( Events, this, child_L, child_R, ObjsIndx, Events_L, Events_R )

            deallocate(Events)
            
            call KDtree_buildTRIANGLES_BreadthFirst( child_L, Events_L, ObjsIndx, level, Depth_First )
            call KDtree_buildTRIANGLES_BreadthFirst( child_R, Events_R, ObjsIndx, level, Depth_First )
            
            nullify(child_L,child_R)

         else
            this% isLast = .false.
            if( this% level .lt. level ) this% isLast = .true.
            do j = 1, size(Depth_First)
               if( .not. Depth_First(j)% active ) then 
                  Depth_First(j)% taskTree => this
                  allocate(Depth_First(j)% taskEvents(NDIM,2*this% NumOfObjs+1))
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
     
      if( allocated(Events) ) deallocate( Events )

   end subroutine KDtree_buildTRIANGLES_BreadthFirst
   
   
   recursive subroutine KDtree_buildTRIANGLES_DepthFirst( this, Events, ObjsIndx )
      use omp_lib
      implicit none
      !-arguments----------------------------------------------------------
      type(KDtree), target,      intent(inout) :: this
      type(Event),  allocatable, intent(inout) :: Events(:,:)
      integer,                   intent(inout) :: ObjsIndx(:)
      !-local-variables---------------------------------------------------
      type(KDtree), pointer     :: child_L, child_R 
      type(Event),  allocatable :: Events_L(:,:), Events_R(:,:)

      call this% GetArea()
!$omp critical
      BoxIndex = BoxIndex + 1
!$omp end critical
      this% index  = BoxIndex

      if( this% level .lt. depth .and. this% NumOfObjs .ge. this% Min_n_of_Objs ) then
         
         if( this% which_KDtree .eq. TRIANGLES_KDTREE_SAH ) then
            call this% EvaluateCostSAH( Events )
         else
            call this% EvaluateCostMEDIAN( Events )
         end if

         if( this% split ) then 

            allocate(this% child_L,this% child_R)
        
            child_L => this% child_L
            child_R => this% child_R
            
            child_L% parent => this
            child_R% parent => this
           
            child_L% level = this% level + 1
            child_R% level = this% level + 1

            this% built_L = .true.
            this% built_R = .true.

            this% isLast = .false.   

            if( this% which_KDtree .eq. TRIANGLES_KDTREE_MEDIAN ) then
               child_L% axis = AxisIndex( this% axis )
               child_R% axis = AxisIndex( this% axis )
            endif

            call this% BuildChild( child_L, side_L )
            call this% BuildChild( child_R, side_R )
         
            child_L% Min_n_of_Objs = this% Min_n_of_Objs
            child_R% Min_n_of_Objs = this% Min_n_of_Objs
            child_L% NumThreads    = this% NumThreads
            child_R% NumThreads    = this% NumThreads
            child_L% which_KDtree  = this% which_KDtree
            child_R% which_KDtree  = this% which_KDtree

            call Event_ClassifyObjs( Events(this% axis,:), this, child_L, child_R, ObjsIndx )

            allocate(Events_L(NDIM,2*child_L% NumOfObjs))
            allocate(Events_R(NDIM,2*child_R% NumOfObjs))

            call Event_BuildLists( Events, this, child_L, child_R, ObjsIndx, Events_L, Events_R )

            deallocate(Events)
  
            call KDtree_buildTRIANGLES_DepthFirst( child_L, Events_L, ObjsIndx )
            call KDtree_buildTRIANGLES_DepthFirst( child_R, Events_R, ObjsIndx )
 
            nullify(child_L,child_R)

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

      if( allocated(Events) ) deallocate(Events)

   end subroutine KDtree_buildTRIANGLES_DepthFirst
   
   subroutine ComputeSplitSurface_L(tree, SplittingPlane, axis, S_L)
   
      implicit none
      !-arguments-------------------------------------
      type(KDtree),  intent(in)  :: tree
      real(kind=rp), intent(in)  :: SplittingPlane
      integer,       intent(in)  :: axis
      real(kind=rp), intent(out) :: S_L
      !-local-variables-------------------------------
      real(kind=rp) :: L_x, L_y, L_z
   
      select case( axis )
         case(IX)
            L_x = abs(SplittingPlane-tree% vertices(1,1))
            L_y = abs(tree% vertices(2,7)-tree% vertices(2,1))
            L_z = abs(tree% vertices(3,7)-tree% vertices(3,1))    
         case(IY)
            L_x = abs(tree% vertices(1,7)-tree% vertices(1,1))
            L_y = abs(SplittingPlane-tree% vertices(2,1))
            L_z = abs(tree% vertices(3,7)-tree% vertices(3,1))
         case(IZ)
            L_x = abs(tree% vertices(1,7)-tree% vertices(1,1))
            L_y = abs(tree% vertices(2,7)-tree% vertices(2,1))
            L_z = abs(SplittingPlane-tree% vertices(3,1))
      end select

      S_L = 2.0_RP*(L_x*L_y+L_y*L_z+L_x*L_z)
   
   end subroutine ComputeSplitSurface_L
   
   subroutine ComputeSplitSurface_R(tree, SplittingPlane, axis, S_R)
   
      implicit none
      !-arguments------------------------------------
      type(KDtree),  intent(in)  :: tree
      real(kind=rp), intent(in)  :: SplittingPlane
      integer,       intent(in)  :: axis
      real(kind=rp), intent(out) :: S_R
      !-local-variables------------------------------
      real(kind=rp) :: L_x, L_y, L_z
   
      select case( axis )
         case(IX)
            L_x = abs(tree% vertices(1,7)-SplittingPlane)
            L_y = abs(tree% vertices(2,7)-tree% vertices(2,1))
            L_z = abs(tree% vertices(3,7)-tree% vertices(3,1))    
         case(IY)
            L_x = abs(tree% vertices(1,7)-tree% vertices(1,1))
            L_y = abs(tree% vertices(2,7)-SplittingPlane)
            L_z = abs(tree% vertices(3,7)-tree% vertices(3,1))
         case(IZ)
            L_x = abs(tree% vertices(1,7)-tree% vertices(1,1))
            L_y = abs(tree% vertices(2,7)-tree% vertices(2,1))
            L_z = abs(tree% vertices(3,7)-SplittingPlane)
      end select

      S_R = 2.0_RP*(L_x*L_y+L_y*L_z+L_x*L_z)
   
   end subroutine ComputeSplitSurface_R
   
   real(kind=rp) function computeSplittingCost(S, S_L, S_R, N_L, N_R) result(SplittingCost)
   
      implicit none
      !-arguments------------------------------
      real(kind=rp), intent(in) :: S, S_L, S_R
      integer,       intent(in) :: N_L, N_R
     
      SplittingCost = C_TRANSVERSE + C_INTERSECT*(S_L/S*N_L+S_R/S*N_R)
   
   end function computeSplittingCost
   
   subroutine KDtree_EvaluateCostSAH( this, Events )
      use omp_lib
      implicit none
      !-arguments----------------------------------------------------------------------------
      class(KDtree), intent(inout) :: this
      type(Event),   intent(in)    :: Events(:,:)
      !-local-variables----------------------------------------------------------------------
      integer                          :: i, k, N_L(NDIM), N_R(NDIM), N_P(NDIM), &
                                          BestSide(1), pminus, pplus, pvert, axis
      real(kind=rp)                    :: SplittingCost(2), S_L, S_R, SplittingPlane

      this% split     = .false.
      this% SplitCost = C_INTERSECT * this% NumOfObjs
      
      if( this% NumOfObjs .eq. 0 ) return
      
      N_R = this% NumOfObjs; N_L = 0; N_P = 0
      this% split = .false.
      this% SplitCost = C_INTERSECT * this% NumOfObjs      
      
      do k = 1, NDIM
         i = 1
         do while ( i .le. this% NumOfEvents(k) )
            pplus = 0; pminus = 0; pvert = 0
            SplittingPlane = Events(k,i)% plane

            if( Events(k,i)% eType .lt. 0 ) exit

            if( i .le. this% NumOfEvents(k) ) then
               do while( Events(k,i)% plane .eq. SplittingPlane .and. Events(k,i)% eType .eq. END_ )
                  pminus = pminus + 1
                  i = i + 1
                  if( i .gt. this% NumOfEvents(k) ) exit
               end do
            end if
        
            if( i .le. this% NumOfEvents(k) ) then
               do while( Events(k,i)% plane .eq. SplittingPlane .and. Events(k,i)% eType .eq. PLANAR_ )
                  pvert = pvert + 1
                  i = i + 1
                  if( i .gt. this% NumOfEvents(k) ) exit
               end do
            end if
            
            if( i .le. this% NumOfEvents(k) ) then
               do while( Events(k,i)% plane .eq. SplittingPlane .and. Events(k,i)% eType .eq. START_ )                                    
                  pplus = pplus + 1
                  i = i + 1
                  if( i .gt. this% NumOfEvents(k) ) exit
               end do   
            end if

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
      
   end subroutine KDtree_EvaluateCostSAH

   subroutine KDtree_EvaluateCostMEDIAN( this, Events )
      use omp_lib
      implicit none
      !-arguments----------------------------------
      class(KDtree), intent(inout) :: this
      type(Event),   intent(in)    :: Events(:,:)
      
      this% SplittingPlane = sum(Events(this% axis,:)% median)/this% NumOfEvents(this% axis)
      this% split = .true.

   end subroutine KDtree_EvaluateCostMEDIAN
   
   subroutine Event_ClassifyObjs( this, tree, child_L, child_R, ObjsIndx )

      implicit none
      !-arguments----------------------------------
      type(Event),  intent(in)    :: this(:)
      type(KDtree), intent(inout) :: tree
      type(KDtree), intent(inout) :: child_L, child_R
      integer,      intent(inout) :: ObjsIndx(:)
      !-local-variables----------------------------------
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
      !-arguments----------------------------
      class(KDtree), intent(inout) :: this
      type(KDtree),  intent(inout) :: child
      integer,       intent(in)    :: side
      !-local-variables----------------------
      integer :: i                                        
   
      child% vertices = this% vertices
      child% isLast = .false.

      select case( this% axis ) 
         case( IX ) 
            do i = 1, 4
               child% vertices(1,vertices_x(side+i)) = this% SplittingPlane
            end do
         case( IY ) 
            do i = 1, 4
               child% vertices(2,vertices_y(side+i)) = this% SplittingPlane
            end do
         case( IZ ) 
            do i = 1, 4
               child% vertices(3,vertices_z(side+i)) = this% SplittingPlane
            end do
      end select
   
   end subroutine KDtree_BuildChild

   subroutine Event_BuildLists( this, tree, child_L, child_R, ObjsIndx, Events_L, Events_R )
      use omp_lib
      implicit none
      !-arguments---------------------------------------------------------------------------------
      type(Event),  intent(in)    :: this(:,:)
      type(KDtree), intent(in)    :: tree
      type(KDtree), intent(inout) :: child_L, child_R
      integer,      intent(inout) :: ObjsIndx(:)
      type(Event),  intent(inout) :: Events_L(:,:), Events_R(:,:)
      !-local-variables---------------------------------------------------------------------------
      integer                          :: i, k, B, L, R

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

   end subroutine Event_BuildLists
   
   subroutine KDtree_SaveObjsIndeces( this, Events )
   
      implicit none
      !-arguments----------------------------------
      class(KDtree), intent(inout) :: this
      type(Event),   intent(in)    :: Events(:,:)
      !-local-variables----------------------------
      integer :: i, k
      
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
      !-arguments------------------------------------------
      type(Event),       intent(inout) :: this(:)
      type(object_type), intent(in)    :: ObjectsList(:)
      type(KDtree),      intent(inout) :: tree
      integer,           intent(in)    :: axis
      !-local-variables------------------------------------
      integer :: i, j, k, shift

      tree% NumOfEvents(axis) = 0

      do i =1, tree% NumOfObjs
         this(i)% index = ObjectsList(i)% index 
         this(i+tree% NumOfObjs)% index = ObjectsList(i)% index 
     
         this(i)% plane = minval(ObjectsList(i)% vertices(1:NumOfVertices)% coords(axis)) 
         this(i+tree% NumOfObjs)% plane = maxval(ObjectsList(i)% vertices(1:NumOfVertices)% coords(axis)) 

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
         this(i)% median = sum(ObjectsList(i)% vertices(1:NumOfVertices)% coords(axis))/NumOfVertices
         this(i+tree% NumOfObjs)% median = this(i)% median
      end do
      
   end subroutine Event_Construct
   
!////////////////// SUBROUTINES FOR KD TREE made of points ///////////////////////
   
   recursive subroutine KDtree_buildPoints_BreadthFirst( this, Events, ObjsIndx, level, Depth_First )
   
      implicit none
      !-arguments---------------------------------------------------
      type(KDtree), target,      intent(inout) :: this
      type(Event),  allocatable, intent(inout) :: Events(:,:)
      integer,                   intent(in)    :: level
      integer,                   intent(inout) :: ObjsIndx(:)
      type(DepthFirst_type),     intent(inout) :: Depth_First(:)
      !-local-variables---------------------------------------------
      type(KDtree), pointer     :: child_L, child_R 
      type(Event),  allocatable :: Events_L(:,:), Events_R(:,:)
      integer                   :: j

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
        
         child_L => this% child_L
         child_R => this% child_R
           
         child_L% level = this% level + 1
         child_R% level = this% level + 1

         this% built_L = .true.
         this% built_R = .true.

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

         call Event_BuildLists( Events, this, child_L, child_R, ObjsIndx, Events_L, Events_R )

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
      
      if( allocated( Events ) ) deallocate( Events )
   
   end subroutine KDtree_buildPoints_BreadthFirst
   
   recursive subroutine KDtree_buildPoints_DepthFirst( this, Events, ObjsIndx )
      use omp_lib
      implicit none
      !-arguments------------------------------------------------
      type(KDtree), target,      intent(inout) :: this
      type(Event),  allocatable, intent(inout) :: Events(:,:)
      integer,                   intent(inout) :: ObjsIndx(:)
      !-local-variables------------------------------------------
      type(KDtree), pointer     :: child_L, child_R 
      type(Event),  allocatable :: Events_L(:,:), Events_R(:,:)
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
        
         child_L => this% child_L
         child_R => this% child_R
           
         child_L% level = this% level + 1
         child_R% level = this% level + 1

         this% built_L = .true.
         this% built_R = .true.

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

         call Event_BuildLists( Events, this, child_L, child_R, ObjsIndx, Events_L, Events_R )

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

      if( allocated( Events ) ) deallocate( Events )

   end subroutine KDtree_buildPoints_DepthFirst
   
   subroutine KDtree_SavePointsIndeces( this, Events )
   
      implicit none
      !-arguments----------------------------------
      class(KDtree), intent(inout) :: this
      type(Event),   intent(in)    :: Events(:,:)
      !-local-variables----------------------------
      integer :: i, k
      
      allocate(this% ObjsIndeces(this% NumOfObjs))
   
      do i = 1, this% NumOfObjs
         this% ObjsIndeces(i) = Events(1,i)% index
      end do   
   
   end subroutine KDtree_SavePointsIndeces
   
   subroutine Points_Event_ClassifyObjs( this, tree, child_L, child_R, ObjsIndx )
      use omp_lib
      implicit none
      !-arguments----------------------------------------
      type(Event),  intent(in)    :: this(:)
      type(KDtree), intent(inout) :: tree
      type(KDtree), intent(inout) :: child_L, child_R
      integer,      intent(inout) :: ObjsIndx(:)
      !-local-variables----------------------------------
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
      !-arguments--------------------------------------
      type(Event),      intent(inout) :: this(:)
      type(KDtree),     intent(inout) :: tree
      integer,          intent(in)    :: axis
      type(point_type), intent(in)    :: PointList(:)
      !-local-varibales--------------------------------
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
      !-arguments-------------------------------------
      type(Event), intent(inout) :: Events(:)
      integer,     intent(in)    :: startIdx, endIdx
      !-local-variables-------------------------------
      type(Event) :: pivot
      integer     :: left, right, mid
     
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
      !-arguments------------------------------------
      type(Event), intent(inout) :: EventA, EventB
      !-local-variables------------------------------
      type(Event) :: temp
  
      temp = EventA
      EventA = EventB
      EventB = temp
  
   end subroutine swapEvents
  
   logical function EventIsLess( EventA, EventB ) result( IsLess )
  
     implicit none
     !-arguments----------------------------------
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
   
   logical function TrinagleIntersectBox( vertices, triangle ) result( inside )
      use MappedGeometryClass
      implicit none 

      real(kind=rp),     intent(in) :: vertices(:,:)
      type(object_type), intent(in) :: Triangle

      real(kind=rp) :: center(NDIM), halfL(NDIM), trivert0(NDIM), trivert1(NDIM), trivert2(NDIM), &
                       v0(NDIM), v1(NDIM), v2(NDIM), e0(NDIM), e1(NDIM), e2(NDIM),                &
                       f0(NDIM), f1(NDIM), f2(NDIM), a00(NDIM), a01(NDIM), a02(NDIM),             &
                       a10(NDIM), a11(NDIM), a12(NDIM), a20(NDIM), a21(NDIM), a22(NDIM),          &
                       a23(NDIM), p0, p1, p2, r, plane_normal(NDIM), plane_distance, Bmin(NDIM),  &
                       Bmax(NDIM)
      integer       :: i


      do i = 1, NDIM 
         Bmin(i)   = minval(vertices(1,:))
         Bmax(i)   = maxval(vertices(1,:))
         center(i) = 0.5_RP*(Bmax(i)+Bmin(i))
         halfL(i)  = 0.5_RP*abs(Bmax(i)-Bmin(i))
      end do 

      trivert0 = triangle% vertices(1)% coords
      trivert1 = triangle% vertices(2)% coords
      trivert2 = triangle% vertices(3)% coords

      if( isInsideBox( trivert0, vertices ) .and. &
          isInsideBox( trivert1, vertices ) .and. &
          isInsideBox( trivert2, vertices ) ) then 
         inside = .true.
         return 
      end if 

      v0 = trivert0 - center;
      v1 = trivert1 - center;
      v2 = trivert2 - center;

      e0 = v1 - v0;
      e1 = v2 - v1;
      e2 = v0 - v2;
  
      f0 = trivert1 - trivert0
      f1 = trivert2 - trivert1
      f2 = trivert0 - trivert2
  
      a00 = (/ 0.0_RP, -f0(3), f0(2) /)
      p0 = dot_product(v0, a00)
      p1 = dot_product(v1, a00)
      p2 = dot_product(v2, a00)
      r = halfL(2) * abs(f0(3)) + halfL(3) * abs(f0(2))
      if (max(-max(p0, p1, p2), min(p0, p1, p2)) > r) then 
         inside = .false.
         return 
      end if 
  
      a01 = (/ 0.0_RP, -f1(3), f1(2) /)
      p0 = dot_product(v0, a01)
      p1 = dot_product(v1, a01)
      p2 = dot_product(v2, a01)
      r = halfL(2) * abs(f1(3)) + halfL(3) * abs(f1(2))
      if (max(-max(p0, p1, p2), min(p0, p1, p2)) > r) then 
         inside = .false. 
         return 
      end if 
  
      a02 = (/ 0.0_RP, -f2(3), f2(2) /)
      p0 = dot_product(v0, a02)
      p1 = dot_product(v1, a02)
      p2 = dot_product(v2, a02)
      r = halfL(2) * abs(f2(3)) + halfL(3) * abs(f2(2))
      if ( max(-max(p0, p1, p2), min(p0, p1, p2)) > r ) then 
         inside= .false.
         return 
      end if 
  
      a10 = (/ f0(3), 0.0_RP, -f0(1) /)
      p0 = dot_product(v0, a10)
      p1 = dot_product(v1, a10)
      p2 = dot_product(v2, a10)
      r = halfL(1) * abs(f0(3)) + halfL(3) * abs(f0(1))
      if ( max(-max(p0, p1, p2), min(p0, p1, p2)) > r) then 
         inside = .false.
         return 
      end if 
  
      a11 = (/ f1(3), 0.0_RP, -f1(1) /)
      p0 = dot_product(v0, a11)
      p1 = dot_product(v1, a11)
      p2 = dot_product(v2, a11)
      r = halfL(1) * abs(f1(3)) + halfL(3) * abs(f1(1))
      if (max(-Max(p0, p1, p2), Min(p0, p1, p2)) > r) then 
         inside = .false.
         return 
      end if 
  
      a12 = (/ f2(3), 0.0_RP, -f2(1) /)
      p0 = dot_product(v0, a12)
      p1 = dot_product(v1, a12)
      p2 = dot_product(v2, a12)
      r = halfL(1) * abs(f2(3)) + halfL(3) * abs(f2(1))
      if (max(-Max(p0, p1, p2), Min(p0, p1, p2)) > r) then 
         inside = .false.
         return 
      end if 
  
      a20 = (/ -f0(2), f0(1), 0.0_RP /)
      p0 = dot_product(v0, a20)
      p1 = dot_product(v1, a20)
      p2 = dot_product(v2, a20)
      r = halfL(1) * abs(f0(2)) + halfL(2) * abs(f0(1))
      if (max(-Max(p0, p1, p2), Min(p0, p1, p2)) > r) then 
         inside = .false.
         return 
      end if 
  
      a21 = (/ -f1(2), f1(1), 0.0_RP /)
      p0 = dot_product(v0, a21)
      p1 = dot_product(v1, a21)
      p2 = dot_product(v2, a21)
      r = halfL(1) * abs(f1(2)) + halfL(2) * abs(f1(1))
      if (max(-Max(p0, p1, p2), Min(p0, p1, p2)) > r) then 
         inside = .false.
         return 
      end if 
  
      a22 = (/-f2(2), f2(1), 0.0_RP/)
      p0 = dot_product(v0, a22)
      p1 = dot_product(v1, a22)
      p2 = dot_product(v2, a22)
      r = halfL(1) * abs(f2(2)) + halfL(2) * abs(f2(1))
      if (max(-Max(p0, p1, p2), Min(p0, p1, p2)) > r) then 
         inside = .false.
         return 
      end if 

      if (Max(v0(1), v1(1), v2(1)) < -halfL(1) .or. Min(v0(1), v1(1), v2(1)) > halfL(1)) then 
         inside = .false.
         return 
      end if 

      if (Max(v0(2), v1(2), v2(2)) < -halfL(2) .or. Min(v0(2), v1(2), v2(2)) > halfL(2)) then 
         inside = .false.
         return 
      end if 

      if (Max(v0(3), v1(3), v2(3)) < -halfL(3) .or. Min(v0(3), v1(3), v2(3)) > halfL(3)) then 
         inside = .false.
         return 
      end if 
  
      call vcross(f0, f1, plane_normal)
      plane_distance = abs(dot_product(plane_normal, v0))
  
      r = halfL(1) * abs(plane_normal(1)) + halfL(2) * abs(plane_normal(2)) + halfL(3) * abs(plane_normal(3));
  
      if (plane_distance > r) then 
         inside = .false.
         return 
      end if 
  
      inside = .true.

   end function TrinagleIntersectBox

end module KDClass
