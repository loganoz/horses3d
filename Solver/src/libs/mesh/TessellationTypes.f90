#include "Includes.h"
module TessellationTypes

   use SMConstants
   use Utilities
   use IntegerDataLinkedList
   use ParamfileRegions                , only: readValueInRegion
   use PhysicsStorage

   implicit none

   public :: DescribeSTLPartitions

   integer, parameter :: FORCING_POINT = 1, NOT_FORCING_POINT = 0, ROTATION = 1, LINEAR = 2
   integer, parameter :: POINT_ON_PLANE = 0, POINT_IN_FRONT_PLANE = 1, POINT_BEHIND_PLANE = 2
   integer, parameter :: NumOfVertices = 3, NumOfIntegrationVertices = 3

!
!  **************************************************
!  Main type for a list of points
!  **************************************************
   type PointLinkedList
   
      class(point_type), pointer :: head => null()
      integer                    :: NumOfPoints
      
      contains
         procedure :: add        => PointLinkedList_add
         procedure :: remove     => PointLinkedList_Remove
         procedure :: removeLast => PointLinkedList_RemoveLast
         procedure :: destruct   => PointLinkedList_Destruct
   
   end type PointLinkedList

!
!  **************************************************
!  Main type for a generic point point
!  **************************************************
   type point_type
   
      class(point_type), pointer :: next => null(), prev => null()
      
      real(kind=rp), dimension(NDIM) :: coords, ImagePoint_coords, normal, xi, VectorValue
      real(kind=rp)                  :: theta, dist, Rank, ScalarValue, xiB, xiI
      integer                        :: index, element_index, NumOfIntersections = 0, &
                                        Translate = 0, partition, objIndex, isForcingPoint, &
                                        STLNum, element_in, faceID, domain = 0, N
      integer,       dimension(NDIM) :: local_Position
      logical                        :: delete = .false., isInsideBody = .false., &
                                        forcingPoint = .false., isInsideBox = .false., state = .false.
      real(kind=RP), allocatable     :: invPhi(:,:), b(:), V(:,:,:), bb(:,:)
      real(kind=RP)                  :: Q(NCONS), U_x(NCONS), U_y(NCONS), U_z(NCONS)
      integer,       allocatable     :: domains(:), indeces(:)  !interPoint, index, domain

      contains
         procedure :: copy => point_type_copy
   
   end type   
   
!
!  **************************************************
!  Main type for a list of points
!  **************************************************
   type ObjectLinkedList
   
      class(object_type), pointer :: head => null()
      integer                     :: NumOfObjs
      
      contains
         procedure :: add        => ObjectLinkedList_add
         procedure :: destruct   => ObjectLinkedList_Destruct
   
   end type ObjectLinkedList 
   
!
!  **************************************************
!  Main type for a generic object
!  **************************************************  
   type Object_type

      class(Object_type), pointer :: next => null(), prev => null()
      
      type(point_type), dimension(:),   allocatable :: vertices, IntegrationVertices
      real(kind=rp),    dimension(NDIM)             :: normal, tangent, coords
      integer                                       :: index, NumOfVertices
      integer,          dimension(2)                :: partition

      contains
         procedure :: copy       => object_type_copy
         procedure :: build      => object_type_build
         procedure :: destruct   => object_type_destruct
                  
   end type Object_type
   
!
!  **************************************************
!  Main type for a STL file reader
!  **************************************************
   type STLfile

      type(Object_type), dimension(:), allocatable :: ObjectsList
      integer                                      :: NumOfObjs, partition,      &
                                                      motionAxis, body,          &
                                                      NumOfObjs_OLD, motionType, &
                                                      MaxAxis, bcType
      real(kind=RP)                                :: angularVelocity, ds,       &
                                                      Velocity,                  & 
                                                      rotationMatrix(NDIM,NDIM), &
                                                      rotationCenter(NDIM)
      logical                                      :: move, show, construct = .false., &
                                                      read = .false.,                  &
                                                      BFcorrection = .false.
      character(len=LINE_LENGTH)                   :: filename, maskName
      
       contains
         procedure :: ReadTessellation
         procedure :: Clip                   => STL_Clip
         procedure :: updateNormals          => STL_updateNormals
         procedure :: SetIntegration         => STL_SetIntegration
         procedure :: ComputeVectorIntegral  => STL_ComputeVectorIntegral
         procedure :: ComputeScalarIntegral  => STL_ComputeScalarIntegral
         procedure :: destroy                => STLfile_destroy
         procedure :: Describe               => Describe_STLfile
         procedure :: Copy                   => STLfile_copy
         procedure :: plot                   => STLfile_plot
         procedure :: SetIntegrationPoints   => STL_SetIntegrationPoints
         procedure :: ResetIntegrationPoints => STL_ResetIntegrationPoints
   end type
   
   type ObjsDataLinkedList_t
      type(ObjData_t), pointer     :: head => null()
      integer                      :: no_of_entries = 0
      contains
         procedure   :: Add      => ObjsDataLinkedList_Add
         procedure   :: check    => CheckObj
         procedure   :: Destruct => ObjsDataLinkedList_Destruct
    end type ObjsDataLinkedList_t
   
   type ObjData_t
      integer                   :: value 
      type(ObjData_t), pointer  :: next 
   end type ObjData_t


   type ObjsRealDataLinkedList_t
      class(ObjRealData_t), pointer :: head => NULL()
      integer                       :: no_of_entries = 0
      contains
         procedure   :: Add      => ObjsRealDataLinkedList_Add
         procedure   :: check    => CheckReal
         procedure   :: Destruct => ObjsRealDataLinkedList_Destruct
   end type ObjsRealDataLinkedList_t
   
   type ObjRealData_t
      real(kind=RP)                 :: value
      class(ObjRealData_t), pointer :: next 
   end type ObjRealData_t

   interface ObjsDataLinkedList_t
      module procedure  ConstructObjsDataLinkedList
   end interface 
    
   interface ObjsRealDataLinkedList_t
     module procedure  ConstructObjsRealDataLinkedList
   end interface 
   
   interface PointLinkedList
      module procedure :: PointLinkedList_Construct
   end interface

   interface ObjectLinkedList
      module procedure :: ObjectLinkedList_Construct
   end interface

   contains   

      function ConstructObjsDataLinkedList( )
         implicit none
         type(ObjsDataLinkedList_t) :: ConstructObjsDataLinkedList 

         ConstructObjsDataLinkedList% head => null()
         ConstructObjsDataLinkedList% no_of_entries = 0

      end function ConstructObjsDataLinkedList
      
      function ConstructObjsRealDataLinkedList( )
         implicit none
         type(ObjsRealDataLinkedList_t) :: ConstructObjsRealDataLinkedList 

         ConstructObjsRealDataLinkedList% head => null()
         ConstructObjsRealDataLinkedList% no_of_entries = 0

      end function ConstructObjsRealDataLinkedList

      subroutine ObjsDataLinkedList_Add( this, value ) 
         implicit none
         !-arguments----------------------------------------------
         class(ObjsDataLinkedList_t), intent(inout) :: this
         integer,                     intent(in)    :: value
         !-local-variables----------------------------------------
         type(ObjData_t), pointer :: current
         integer                  :: i

         if ( this% no_of_entries .eq. 0 ) then
            allocate( this% head ) 
            this% head% value = value
            this% no_of_entries = 1
         else
            current => this% head    
            do i = 1, this% no_of_entries-1
               current => current% next
            end do
            allocate(current% next)
            current% next% value = value
            this% no_of_entries = this% no_of_entries + 1 
         end if

      end subroutine ObjsDataLinkedList_Add
      
      subroutine ObjsRealDataLinkedList_Add( this, value ) 
         implicit none
         !-arguments---------------------------------------------
         class(ObjsRealDataLinkedList_t), intent(inout) :: this
         real(kind=RP),                   intent(in)    :: value
         !-local-variables---------------------------------------
         type(ObjRealData_t), pointer  :: current
         integer                       :: i

         if ( this% no_of_entries .eq. 0 ) then
            allocate( this% head ) 
            this% head% value = value
            this% no_of_entries = 1
         else
            current => this% head    
            do i = 1, this% no_of_entries-1
               current => current% next
            end do
            allocate(current% next)
            current% next% value = value
            this% no_of_entries = this% no_of_entries + 1 
         end if

      end subroutine ObjsRealDataLinkedList_Add

      logical function CheckObj( this, value ) result( found )
         implicit none
         !-arguments---------------------------------------------
         class(ObjsDataLinkedList_t), intent(inout) :: this
         integer,                     intent(in)    :: value
         !-local-variables---------------------------------------
         type(ObjData_t), pointer :: current
         integer                  :: i
        
         found = .false.
        
         if( this% no_of_entries .eq.0 ) return
        
         current => this% head
        
         do i = 1, this% no_of_entries
            if( current% value .eq. value ) then
               found = .true.
               exit
            end if
            current => current% next
         end do
     
      end function CheckObj

      logical function CheckReal( this, value ) result( found )
         implicit none
         !-arguments----------------------------------------------
         class(ObjsRealDataLinkedList_t), intent(inout) :: this
         real(kind=RP),                   intent(in)    :: value
         !-local-variables----------------------------------------
         type(ObjRealData_t), pointer :: current
         integer                      :: i
        
         found = .false.
        
         if( this% no_of_entries .eq.0 ) return
        
         current => this% head
        
         do i = 1, this% no_of_entries
            if( almostEqual(current% value,value) ) then
               found = .true.
               exit
            end if
            current => current% next
         end do
     
      end function CheckReal

      elemental subroutine ObjsDataLinkedList_Destruct(this)
         implicit none
         !-arguments---------------------------------------------
         class(ObjsDataLinkedList_t), intent(inout) :: this
         !-local-variables---------------------------------------
         type(ObjData_t), pointer :: data, nextdata
         integer                  :: i

         data => this% head
         do i = 1, this% no_of_entries
            nextdata => data% next

            deallocate(data)
            data => nextdata
         end do
         
         this% no_of_entries = 0

      end subroutine ObjsDataLinkedList_Destruct
      
      elemental subroutine ObjsRealDataLinkedList_Destruct(this)
         implicit none
         !-arguments---------------------------------------------
         class(ObjsRealDataLinkedList_t), intent(inout) :: this
         !-local-variables---------------------------------------
         type(ObjRealData_t), pointer :: data, nextdata
         integer                      :: i

         data => this% head
         do i = 1, this% no_of_entries
            nextdata => data% next

            deallocate(data)
            data => nextdata
         end do
         
         this% no_of_entries = 0

      end subroutine ObjsRealDataLinkedList_Destruct

!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
!  This subroutine initializes a points list
!  ------------------------------------------------ 
   function PointLinkedList_Construct(  )
   
      implicit none
      
      type(PointLinkedList) :: PointLinkedList_Construct
      
      PointLinkedList_Construct% head => null()
   
   end function PointLinkedList_Construct
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
!  This subroutine adds a point to a points list
!  ------------------------------------------------ 
   subroutine PointLinkedList_Add( this, point )
      implicit none
      !-arguments-----------------------------------
      class(PointLinkedList) :: this
      type(point_type)       :: point
      !-local-variables-----------------------------
      type(point_type), pointer :: current,    &
                                   currentNext 
   
      if( .not. associated(this% head) ) then 
         allocate(this% head)
         call this% head% copy(point)
         this% head% next => this% head
         this% head% prev => this% head
      else
         current => this% head% prev
         currentNext => current% next
         allocate(currentNext) 
         call currentNext% copy(point)
         currentNext% next => this% head
         currentNext% prev => current
         nullify(current% next)
         current% next => currentNext
         nullify(this% head% prev)
         this% head% prev => currentNext
      end if
      
      nullify(current, currentNext)
   
   end subroutine PointLinkedList_Add
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
!  This subroutine deletes the point p from a points list
!  ------------------------------------------------  
   subroutine PointLinkedList_Remove( this, p )
      implicit none
      !-arguments--------------------------------------------------------------
      class(PointLinkedList),   intent(inout) :: this
      type(point_type), target, intent(inout) :: p
      !-local-variables--------------------------------------------------------
      type(point_type), pointer :: dataNext => null(), dataPrev => null(), &
                                   dataDel => null()
      
      dataDel => p
      
      dataPrev => p% prev; dataNext => p% next
      
      if( associated(dataDel, this% head) ) then
         error stop ":: Head of list can't be deleted here!"
         return
      end if
      
      dataPrev% next => null()
      dataPrev% next => dataNext
      dataNext% prev => null()
      dataNext% prev => dataPrev
   
      nullify(dataPrev, dataNext)
   
   end subroutine PointLinkedList_Remove
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
!  This subroutine deletes the last point from a points list
!  ------------------------------------------------  
   subroutine PointLinkedList_RemoveLast( this )
      implicit none
      !-arguments--------------------------------------------------------
      class(PointLinkedList), intent(inout) :: this
      !-local-variables--------------------------------------------------
      type(point_type), pointer :: data => null(), dataPrev => null()
      
      data => this% head% prev
   
      if( associated(data, this% head) ) then
         error stop ":: Head of list can't be deleted here!"
         return
      end if
      
      dataPrev => data% prev
      
      deallocate(data)
      
      dataPrev% next => null()
      dataPrev% next => this% head
      this% head% prev => null() 
      this% head% prev => dataPrev 
   
      nullify(dataPrev)
   
   end subroutine PointLinkedList_RemoveLast
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
!  This subroutine deletes a point from a points list
!  ------------------------------------------------ 
   subroutine PointLinkedList_Destruct( this )
      implicit none
      !-arguments--------------------------------------
      class(PointLinkedList), intent(inout) :: this
      !-local-variables--------------------------------
      class(point_type), pointer :: current, next 
      integer                    :: i
      
      if( this% NumOfPoints .eq. 0 ) return
      
      current => this% head
      next    => current% next

      do i = 2, this% NumOfPoints
         deallocate(current)
         current => next
         next    => current% next
      end do
      
      this% NumOfPoints = 0
      nullify(current, next)
      
   end subroutine PointLinkedList_Destruct
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine copies a point
!  ------------------------------------------------  

   subroutine point_type_copy( this, point )
      implicit none
      !-arguments-------------------------------
      class(point_type), intent(inout) :: this
      type(point_type),  intent(in)    :: point
    
      this% coords            = point% coords
      this% ImagePoint_coords = point% ImagePoint_coords
      this% theta             = point% theta
      this% index             = point% index
      this% element_index     = point% element_index
      this% local_position    = point% local_position
      this% Translate         = point% Translate
      this% partition         = point% partition
   
   end subroutine point_type_copy
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine initializes an objects list
!  ------------------------------------------------  
   function ObjectLinkedList_Construct(  )
      implicit none
      !-local-variables-------------------------------------
      type(ObjectLinkedList) :: ObjectLinkedList_Construct
      
      ObjectLinkedList_Construct% head => null()
      ObjectLinkedList_Construct% NumOfObjs = 0
   
   end function ObjectLinkedList_Construct
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine adds an object to an objects list
!  ------------------------------------------------  

   subroutine ObjectLinkedList_Add( this, object )
      implicit none
      !arguemnts--------------------------------------------
      class(ObjectLinkedList) :: this
      type(Object_type)       :: object
      !-local-variables-------------------------------------
      type(object_type), pointer :: current => null(), &
                                    currentNext => null()
   
      if( .not. associated(this% head) ) then 
         allocate(this% head)
         call this% head% copy(object)
         this% head% next => this% head
         this% head% prev => this% head
      else
         current => this% head% prev
         currentNext => current% next
         allocate(currentNext) 
         call currentNext% copy(object)
         currentNext% next => this% head
         currentNext% prev => current
         nullify(current% next)
         current% next => currentNext
         nullify(this% head% prev)
         this% head% prev => currentNext
         nullify(current, currentNext)
      end if
      
      this% NumOfObjs = this% NumOfObjs + 1
   
   end subroutine ObjectLinkedList_Add
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
!  This subroutine destroys an objects list
!  ------------------------------------------------ 
   subroutine ObjectLinkedList_Destruct( this )
      implicit none
      !arguemnts----------------------------------------
      class(ObjectLinkedList), intent(inout) :: this
      !-local-variables---------------------------------
      type(object_type), pointer :: data => null(), &
                                    dataPrev => null()
      
      if( .not. associated(this% head) ) return
      
      data => this% head% prev
      
      if( .not. associated(data, this% head) ) then
         do
            dataPrev => data% prev   
            call data% destruct()
            deallocate(data)
            data => dataPrev
            if( associated(data, this% head) ) exit
         end do
      end if
      
      deallocate(data)
      
      this% NumOfObjs = this% NumOfObjs - 1
      
      nullify(dataPrev)
      
   end subroutine ObjectLinkedList_Destruct
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine builds an object
!  -----------------------------------------------
   subroutine object_type_build( this, Points, normal, NumOfVertices, index )
      implicit none
      !-arguments---------------------------
      class(Object_type), intent(inout) :: this
      real(kind=RP),      intent(in)    :: Points(:,:)
      real(kind=RP),      intent(in)    :: normal(:)
      integer,            intent(in)    :: NumOfVertices, index
      !-local-variables--------------------------------------
      integer :: i
      
      if( allocated(this% vertices) ) deallocate(this% vertices)
      allocate(this% vertices(NumOfVertices))
      
      do i = 1, NumOfVertices
         this% vertices(i)% coords = Points(:,i)
      end do
      
      this% index            = index 
      this% normal           = normal 
      this% NumOfVertices    = NumOfVertices
      
   end subroutine object_type_build
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine copies an object
!  -----------------------------------------------
   subroutine object_type_copy( this, Object )
    
      implicit none
      !-arguments---------------------------
      class(Object_type), intent(inout) :: this
      type(Object_type),  intent(in)    :: Object
      !-local-variables-----------------------
      integer :: i
      
      allocate(this% vertices(Object% NumOfVertices))
      
      do i = 1, Object% NumOfVertices
         call this% vertices(i)% copy(Object% vertices(i))
      end do
      
      this% index            = Object% index 
      this% normal           = Object% normal 
      this% NumOfVertices    = Object% NumOfVertices
      this% partition        = Object% partition
      
   end subroutine object_type_copy
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine copies an object
!  -----------------------------------------------
   subroutine object_type_destruct( this )
    
      implicit none
      !-arguments---------------------------
      class(Object_type), intent(inout) :: this
      
      deallocate(this% vertices)
      
      this% index  = 0
      this% normal = 0.0_RP 
      this% NumOfVertices = 0
      
   end subroutine object_type_destruct
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine adds an object to a list of objects
!  -----------------------------------------------
   subroutine addObj( ObjList, Object )
   
      implicit none
      !-arguments----------------------------------------------------
      type(Object_type), target, intent(inout) :: ObjList
      type(Object_type),         intent(in)    :: Object
      !-local-variables----------------------------------------------
      type(Object_type), pointer :: Obj=>null(), Objprev => null() 
     
      if( .not. associated(ObjList% next)) then
         call ObjList% copy( Object )
         ObjList% next => ObjList
         ObjList% prev => ObjList
         return
      end if
      
      Obj => ObjList% prev
!            
!     previous point
!     ------------
      Objprev => ObjList% prev
      
      Obj => Obj% next
      allocate(Obj)
      
      Obj% next => ObjList
      Obj% prev => Objprev
      
      Objprev% next => Obj
      
      call Obj% copy( Object )
      
      ObjList% prev => null()
      ObjList% prev => Obj
      
      if( associated(Obj% prev, ObjList) ) then
         ObjList% next => Obj 
      end if
      
      nullify(Obj, Objprev)
      
   end subroutine addObj
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine describes the .stl file
!  -----------------------------------------------
   subroutine Describe_STLfile( this, filename )
      use Headers
      use MPI_Process_Info
      use PhysicsStorage
      use BoundaryConditions 
      implicit none
      !-arguments-----------------------------------
      class(STLfile),   intent(inout) :: this
      character(len=*), intent(in)    :: filename

      if( MPI_Process% isRoot) then
       
         write(STD_OUT,'(/)')
         call Section_Header("Reading stl file")
         write(STD_OUT,'(/)')
         call SubSection_Header('Stl file "' // trim(fileName) // '"')

         write(STD_OUT,'(30X,A,A35,I10)') "->" , "Number of objects: " , this% NumOfObjs
         write(STD_OUT,'(30X,A,A35,L10)') "->" , "Moving: " , this% move
         if( .not. this% move ) write(STD_OUT,'(30X,A,A35,A)') "->" , "Boundary conditions: " , trim(implementedBCNames(this% bctype))
         write(STD_OUT,'(30X,A,A35,L10)') "->" , "BF correction: " , this% BFcorrection
         if( this% read ) then 
            write(STD_OUT,'(30X,A,A35,A)') "->" , "Mask file: " , this% maskName
         end if 
         
      end if
   
   end subroutine Describe_STLfile
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine describes the .stl file
!  -----------------------------------------------
   subroutine DescribeSTLPartitions( partition, NumOfObjs )
      use Headers
      use MPI_Process_Info
      implicit none
      !-arguments--------------------------------
      integer, intent(in) :: partition, NumOfObjs
      !-local-variables--------------------------
      character(len=LINE_LENGTH) :: myString
      
      if( .not. MPI_Process% doMPIAction ) return

      write(myString,'(i100)') partition+1
      
      if( partition .eq. 0 ) then
         write(STD_OUT,'(/)')
         call Section_Header("stl partitions")
         write(STD_OUT,'(/)')
      end if
      
      call SubSection_Header('partition ' // trim(adjustl(myString)))
      write(STD_OUT,'(30X,A,A32,I10)') "->" , "Number of objects: " , NumOfObjs
      
   end subroutine DescribeSTLPartitions
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine reads the .stl file
!  -----------------------------------------------
   subroutine  ReadTessellation( this, filename )
      use PhysicsStorage
      use FileReadingUtilities, only: getFileExtension
      use,intrinsic :: iso_c_binding
      implicit none
      !-arguments---------------------------------------
      class(STLfile),             intent(inout) :: this
      character(len=*),           intent(in)    :: filename
      !-local-variables---------------------------------
      integer*4              :: i, j, funit, NumOfTri,      &
                              fileStat, NumOfVertices, ios
      integer*2            :: padding
      real*4, dimension(3) :: norm, vertex
      character(len=80)    :: header
      real(kind=RP)        :: max_x, max_y, max_z
      real(kind=RP)        :: min_x, min_y, min_z

      character(len=7)  :: cCAD , stlAutoCAD = "AutoCAD"
      character(len=21) :: cGMSH, stlGMSH    = "solid Created by Gmsh"
      logical                    :: formatted = .false.

      NumOfVertices = 3 
      
      this% partition = 1

      max_x = -huge(1.0_RP); max_y = -huge(1.0_RP); max_z = -huge(1.0_RP)
      min_x =  huge(1.0_RP); min_y =  huge(1.0_RP); min_z =  huge(1.0_RP)
      
      funit = UnusedUnit()
      open(unit=funit, file=trim(this% filename), status='old', form='unformatted', action='read', access='stream', iostat=fileStat)

      if( fileStat .ne. 0 ) then
         print *, "Read Tessellation: file '",trim(this% filename),"' not found"
         error stop
      end if

      open(unit=funit, file=trim(this% filename), status='old', form='unformatted', action='read', access='stream', iostat=fileStat) 

      read(funit) header
      read(funit) NumOfTri

      this% NumOfObjs = NumOfTri 
      
      allocate(this% ObjectsList(NumOfTri))

      associate( Objs => this% ObjectsList )
      
      do i = 1, NumOfTri
         read(funit) norm(1), norm(2), norm(3)
         Objs(i)% normal = norm
         allocate(Objs(i)% vertices(NumOfVertices))
         do j = 1, NumOfVertices
            read(funit) vertex(1), vertex(2), vertex(3)
            max_x = max(max_x, vertex(1))
            max_y = max(max_y, vertex(2))
            max_z = max(max_z, vertex(3))
            min_x = min(min_x, vertex(1))
            min_y = min(min_y, vertex(2))
            min_z = min(min_z, vertex(3))
            Objs(i)% vertices(j)% coords = vertex !/Lref -> always 1                                     
         end do
         read(funit) padding
         Objs(i)% index = i
         Objs(i)% NumOfVertices = NumOfVertices
         Objs(i)% partition = 1
      end do
      
      end associate  
      
      this% MaxAxis = maxloc((/abs(max_x-min_x),abs(max_y-min_y),abs(max_z-min_z)/),dim=1)

      close(unit=funit)

   end subroutine  ReadTessellation
 
   subroutine STL_SetIntegrationPoints( this )

      implicit none

      class(STLfile), intent(inout) :: this 

      integer       :: i, j, m,                &
                      indecesL(NumOfVertices), &
                      indecesR(NumOfVertices)
!          * 2
!         / \
!      4 *   * 5 
!       /  *   \
!      /   7    \
!   1 *----*-----* 3
!          6
      do i = 1, this% NumOfObjs 
         associate(obj => this% ObjectsList(i))
         if( .not. allocated(obj% IntegrationVertices) ) allocate( obj% IntegrationVertices(NumOfIntegrationVertices) )
         obj% IntegrationVertices(NumOfIntegrationVertices)% coords = 0.0_RP 
         do j = 1, NumOfVertices
            obj% IntegrationVertices(j)% coords                        = obj% vertices(j)% coords
            ! obj% IntegrationVertices(NumOfIntegrationVertices)% coords = obj% IntegrationVertices(NumOfIntegrationVertices)% coords + &
            !                                                              obj% vertices(j)% coords
         end do 
         ! obj% IntegrationVertices(NumOfIntegrationVertices)% coords = obj% IntegrationVertices(NumOfIntegrationVertices)% coords/NumOfVertices
         ! indecesL = (/ 1, 2, 3 /)
         ! indecesR = (/ 2, 3, 1 /)
         ! m = 0
         ! do j = NumOfVertices+1, NumOfIntegrationVertices-1
         !    m = m + 1
         !    obj% IntegrationVertices(j)% coords = 0.5_RP*( obj% IntegrationVertices(indecesL(m))% coords + &
         !                                                   obj% IntegrationVertices(indecesR(m))% coords   )
         ! end do
         end associate 
      end do

   end subroutine STL_SetIntegrationPoints

   subroutine STLfile_plot( this, iter )
      use MPI_Process_Info
      use PhysicsStorage
      use FileReadingUtilities, only :getFileName
      implicit none
      !-arguments----------------------------------------------
      class(STLfile),    intent(inout) :: this
      integer,           intent(in)    :: iter
      !-local-variables----------------------------------------
      character(len=LINE_LENGTH)     :: filename, myString, name 
      integer                        :: i, j, funit, index_, lastIndex_
      integer*2                      :: padding = 0
      real*4, dimension(3)           :: norm, vertex
      character(len=80)              :: header = repeat(' ',80), rank 
            
      if( .not. MPI_Process% isRoot ) return 

      funit = UnusedUnit()
     
      name = getFileName(this% filename)

      myString = trim(name)
      index_   = INDEX(myString, '_') 
      lastIndex_ = 0
      do while (index_ /= 0)
         lastIndex_ = index_
         index_     = INDEX(myString(lastIndex_+1:), '_')
         if (index_ /= 0) then
            index_ = index_ + lastIndex_
         end if
      end do

      if( lastIndex_ .ne. 0 ) then 
         write(filename,'(A,A,I10.10)') trim(name(1:lastIndex_-1)),'_', iter
      else 
         write(filename,'(A,A,I10.10)') trim(name),'_', iter
      end if 

      open(funit,file=trim(filename)//'.stl', status='unknown',access='stream',form='unformatted')
 
      write(funit) header, this% NumOfObjs 
      
      do i = 1, this% NumOfObjs
         norm = this% ObjectsList(i)% normal
         write(funit) norm(IX), norm(IY), norm(IZ)
         do j = 1, NumOfVertices      
            vertex = this% ObjectsList(i)% vertices(j)% coords * Lref
            write(funit) vertex(IX), vertex(IY), vertex(IZ)
         end do
         write(funit) padding      
      end do

      close(funit)
   
   end subroutine STLfile_plot

   subroutine STLfile_copy( this, STL )

      implicit none 

      class(STLfile), intent(inout) :: this
      type(STLfile),  intent(in)    :: STL

      integer :: NumOfObjs, i, j

      this% filename = trim(STL% filename)
      NumOfObjs      = STL% NumOfObjs

      this% NumOfObjs = NumOfObjs

      allocate( this% ObjectsList(NumOfObjs) )
!$omp parallel
!$omp do schedule(runtime) private(j)
      do i = 1, NumOfObjs
         allocate(this% ObjectsList(i)% vertices(NumOfVertices))
         do j = 1, NumOfVertices
            this% ObjectsList(i)% vertices(j)% coords = STL% ObjectsList(i)% vertices(j)% coords
         end do
         this% ObjectsList(i)% normal        = STL% ObjectsList(i)% normal
         this% ObjectsList(i)% index         = STL% ObjectsList(i)% index
         this% ObjectsList(i)% NumOfVertices = NumOfVertices
         this% ObjectsList(i)% partition     = 1
      end do
!$omp end do
!$omp end parallel
   end subroutine STLfile_copy
  
   subroutine STLfile_GetInfo( this, STLfilename )
      use FileReadingUtilities
      use FTValueDictionaryClass
      use PhysicsStorage
      use FileReaders
      use BoundaryConditions
      use GenericBoundaryConditionClass, only: CheckIfStlNameIsContained, CheckIfBoundaryNameIsContained
      implicit none
      !-arguments-------------------------------------------------------
      type(STLfile),    intent(inout) :: this
      character(len=*), intent(inout) :: STLfilename
      !-local-arguments-------------------------------------------------
      integer                    :: fid, i, bctype, io, GetZoneType
      logical                    :: inside, logval
      character(len=LINE_LENGTH) :: currentLine, loweredBname
      character(len=LINE_LENGTH) :: keyword, keyval
      type(FTValueDIctionary)    :: bcdict
      character(LINE_LENGTH)     :: ext, stlname, stlNoPath

      ext = getFileExtension(trim(this% filename))

      if ( trim(ext) .ne. 'stl') then
         print *, "Read Tessellation: extension .",trim(ext),"' not supported"
         error stop
      end if

      stlNoPath = RemovePath( trim(this% filename) )   
      stlname   = getFileName( trim(stlNoPath) )

      STLfilename = trim(stlname)
 
      loweredbName = trim(adjustl(STLfilename))
      call toLower(loweredbName)

      call bcdict % initWithSize(16) 

      open(newunit = fid, file = trim(controlFileName), status = "old", action = "read")

      inside = .false.
      do 
         read(fid, '(A)', iostat=io) currentLine

         if( io .ne. 0 ) exit

         call PreprocessInputLine(currentLine)
         call toLower(currentLine)

         if ( index(trim(currentLine),"#define stl") .ne. 0 ) then
            inside = CheckIfStlNameIsContained(trim(currentLine), trim(loweredbname)) 
         end if

!
!           Get all keywords inside the zone
!           --------------------------------
         if ( inside ) then
            if ( trim(currentLine) .eq. "#end" ) exit

            keyword  = ADJUSTL(GetKeyword(currentLine))
            keyval   = ADJUSTL(GetValueAsString(currentLine))
            call ToLower(keyword)
   
            call bcdict % AddValueForKey(keyval, trim(keyword))
         end if

      end do

      if ( .not. bcdict % ContainsKey("type") ) then
         bctype = 3
      end if

      keyval = bcdict % StringValueForKey("type", LINE_LENGTH)
      call tolower(keyval)
      
      GetZoneType = -1
      do bctype = 1, size(implementedBCnames)
         if ( trim(keyval) .eq. trim(implementedBCnames(bctype)) ) then
            GetZoneType = bctype
         end if
      end do
      
      if ( GetZoneType .eq. -1 ) then
         GetZoneType = 3 
      end if

      if ( .not. bcdict % ContainsKey("moving") ) then
         logval = .false. 
      else 
         logval = bcdict % logicalValueForKey("moving")
      end if 

      this% bctype = GetZoneType
      this% move   = logval

      call bcdict % Destruct
      close(fid)

      ! check if stl is BF correctd

      open(newunit = fid, file = trim(controlFileName), status = "old", action = "read")

      inside = .false.
      do 
         read(fid, '(A)', iostat=io) currentLine

         if( io .ne. 0 ) exit

         call PreprocessInputLine(currentLine)
         call toLower(currentLine)

         if ( index(trim(currentLine),"#define boundary") .ne. 0 ) then
            inside = CheckIfBoundaryNameIsContained(trim(currentLine), trim(loweredbname)) 
         end if

!
!           Get all keywords inside the zone
!           --------------------------------
         if ( inside ) then
            if ( trim(currentLine) .eq. "#end" ) exit

            keyword  = ADJUSTL(GetKeyword(currentLine))
            keyval   = ADJUSTL(GetValueAsString(currentLine))
            call ToLower(keyword)
            
            this% BFcorrection = .true. 
         end if

      end do

      close(fid)

   end subroutine STLfile_GetInfo

   subroutine STL_updateNormals( this, obj )
      use MappedGeometryClass
      implicit none
      !-arguments------------------------------------
      class(STLfile),    intent(inout) :: this
      type(Object_type), intent(inout) :: obj
      !-local-variables------------------------------
      real(kind=rp) :: du(NDIM), dv(NDIM), dw(NDIM)
      integer       :: i 

      du = obj% vertices(2)% coords - obj% vertices(1)% coords
      dw = obj% vertices(3)% coords - obj% vertices(1)% coords
      call vcross(du,dw,dv)
      obj% normal = dv/norm2(dv)      
 
   end subroutine STL_updateNormals
   
   subroutine STLfile_destroy( this )
   
      implicit none
      !-arguments------------------------------
      class(STLfile), intent(inout) :: this
      !-local-variables------------------------
      integer :: i, j
      
      do i = 1, this% NumOfObjs
         deallocate(this% ObjectsList(i)% vertices)
         if( allocated(this% ObjectsList(i)% IntegrationVertices) ) then 
            do j = 1, NumOfIntegrationVertices
               deallocate( this% ObjectsList(i)% IntegrationVertices(j)% domains, &
                           this% ObjectsList(i)% IntegrationVertices(j)% indeces, & 
                           this% ObjectsList(i)% IntegrationVertices(j)% invPhi,  &
                           this% ObjectsList(i)% IntegrationVertices(j)% b        )
            end do 
            deallocate(this% ObjectsList(i)% IntegrationVertices)
         end if 
      end do
      
      deallocate(this% ObjectsList)  

      this% NumOfObjs = 0
   
   end subroutine STLfile_destroy
!//////////////////////////////////////////////
!
!   Procedures form STL triangles splitting   
!   
!//////////////////////////////////////////////
! SPLITTING TRIANGLES 
   subroutine STL_Clip( this, minplane, maxplane, axis, describe )
      implicit none
      !-arguments--------------------------------------------------------------------
      class(STLfile), intent(inout) :: this
      real(kind=RP),  intent(in)    :: minplane, maxplane
      integer,        intent(in)    :: axis
      logical,        intent(in)    :: describe
      !-local-variables--------------------------------------------------------------
      type(ObjectLinkedList)     :: ObjectsLinkedList, ObjectsLinkedListFinal
      type(object_type), pointer :: obj 
      real(kind=RP)              :: minplane_point(NDIM), maxplane_point(NDIM),   & 
                                    minplane_normal(NDIM), maxplane_normal(NDIM), &
                                    Objmax, Objmin, vertices(NDIM,3),             &
                                    AB(NDIM), AC(NDIM), normal(NDIM)
      integer                    :: i, j

      minplane_point  = 0.0_RP
      maxplane_point  = 0.0_RP
      minplane_normal = 0.0_RP
      maxplane_normal = 0.0_RP
      
      minplane_point(axis)  = minplane
      maxplane_point(axis)  = maxplane
      minplane_normal(axis) = -1.0_RP
      maxplane_normal(axis) =  1.0_RP
 
      ObjectsLinkedList      = ObjectLinkedList_Construct()
      ObjectsLinkedListFinal = ObjectLinkedList_Construct()

      do i = 1, this% NumOfObjs
          call ClipPoly( this% ObjectsList(i), maxplane_normal, maxplane_point, ObjectsLinkedList )
      end do

      obj => ObjectsLinkedList% head 

      do i = 1, ObjectsLinkedList% NumOfObjs
         call ClipPoly( obj, minplane_normal, minplane_point, ObjectsLinkedListFinal )
         obj => obj% next 
      end do  
      
      call ObjectsLinkedList% destruct()       

      call this% destroy()

      this% partition = 1
      this% NumOfObjs = ObjectsLinkedListFinal% NumOfObjs
  
      allocate(this% ObjectsList(this% NumOfObjs))
 
      obj => ObjectsLinkedListFinal% head 

      do i = 1, this% NumOfObjs
         allocate(this% ObjectsList(i)% vertices(obj% NumOfVertices))
         do j = 1, obj% NumOfVertices
            this% ObjectsList(i)% vertices(j)% coords = obj% vertices(j)% coords
         end do

         AB = this% ObjectsList(i)% vertices(2)% coords - this% ObjectsList(i)% vertices(1)% coords 
         AC = this% ObjectsList(i)% vertices(3)% coords - this% ObjectsList(i)% vertices(1)% coords

         normal(IX) = AB(2)*AC(3) - AB(3)*AC(2)
         normal(IY) = AB(3)*AC(1) - AB(1)*AC(3)
         normal(IZ) = AB(1)*AC(2) - AB(2)*AC(1)
         
         if( almostEqual(norm2(normal),0.0_RP) ) then 
            normal = 0.0_RP 
         else 
            normal = normal/norm2(normal)
         end if 

         this% ObjectsList(i)% normal        = normal!obj% normal
         this% ObjectsList(i)% index         = i
         this% ObjectsList(i)% NumOfVertices = obj% NumOfVertices
         this% ObjectsList(i)% partition     = 1
         obj => obj% next 
      end do 

      call ObjectsLinkedListFinal% destruct()

      if( describe ) call this% describe(this% filename)

      call this% plot(0)

   end subroutine STL_Clip

   subroutine ClipPoly( obj, plane_normal, plane_point, ObjectsLinkedList )
      use MappedGeometryClass
      implicit none
      !-arguments--------------------------------------------------------------------
      type(object_type),      intent(in)    :: obj
      real(kind=rp),          intent(in)    :: plane_normal(:), plane_point(:)
      type(ObjectLinkedList), intent(inout) :: ObjectsLinkedList
      !-local-variables--------------------------------------------------------------
      real(kind=RP)     :: PointFront(NDIM,4), PointBack(NDIM,4)
      type(object_type) :: objBack
      real(kind=RP)     :: PointA(NDIM), PointB(NDIM), Point_inters(NDIM), v(NDIM),u(NDIM),W(NDIM)
      integer           :: PointA_Is, PointB_Is, n_front, n_back, i 

      n_front = 0; n_back = 0
      
      pointA = obj% vertices(obj% NumOfVertices)% coords
      
      PointA_Is = Point_wrt_Plane( plane_normal, plane_point, pointA )
   
      do i = 1, obj% NumOfVertices
         PointB    = obj% vertices(i)% coords
         PointB_Is = Point_wrt_Plane( plane_normal, plane_point, pointB )
         if( PointB_Is .eq. POINT_IN_FRONT_PLANE ) then
            if( PointA_Is .eq. POINT_BEHIND_PLANE ) then
               Point_inters = EdgePlaneIntersection( plane_normal, plane_point, PointA, PointB )
               n_front = n_front + 1
               n_back  = n_back + 1
               PointFront(:,n_front) = Point_Inters
               PointBack(:,n_back) = Point_Inters
            end if
            n_front = n_front + 1
            PointFront(:,n_front) = PointB
         elseif( PointB_Is .eq. POINT_BEHIND_PLANE ) then
            if( PointA_Is .eq. POINT_IN_FRONT_PLANE ) then
               Point_inters = EdgePlaneIntersection( plane_normal, plane_point, PointA, PointB )
               n_front = n_front + 1
               n_back  = n_back + 1
               PointFront(:,n_front) = Point_Inters
               PointBack(:,n_back) = Point_Inters
            elseif( PointA_Is .eq. POINT_ON_PLANE ) then
               n_back  = n_back + 1
               PointBack(:,n_back) = PointA
            end if
            n_back  = n_back + 1
            PointBack(:,n_back) = PointB 
         else
            n_front = n_front + 1
            PointFront(:,n_front) = PointB
            if( PointA_Is .eq. POINT_BEHIND_PLANE ) then
               n_back  = n_back + 1
               PointBack(:,n_back) = PointB 
            end if
         end if
         PointA    = PointB
         PointA_Is = PointB_Is 
      end do
      
      ! take only back elements !! 
      if( n_back .eq. 3 ) then
         call objBack% build( PointBack(:,1:n_back), obj% normal, obj% NumOfVertices, obj% index )
         call ObjectsLinkedList% add(objBack)
         call objBack% destruct()
      elseif( n_back .eq. 4 ) then
         call objBack% build( PointBack(:,1:n_back-1), obj% normal, obj% NumOfVertices, obj% index )
         call ObjectsLinkedList% add(objBack)
         call objBack% destruct()
         call objBack% build( PointBack(:,(/n_back-1,n_back,1/)), obj% normal, obj% NumOfVertices, obj% index )
         call ObjectsLinkedList% add(objBack)
         call objBack% destruct()
      elseif( n_back .eq. 0 ) then 
      else
         print *, "ClipPloy:: wrong number of vertices: ", n_back
         error stop
      end if 
   
   end subroutine ClipPoly

   subroutine STL_SetIntegration( this, NumOfInterPoints )
      use MPI_Process_Info
      implicit none 

      class(STLfile), intent(inout) :: this 
      integer,        intent(in)    :: NumOfInterPoints

      integer :: i, j
      
      do i = 1, this% NumOfObjs
         do j = 1, NumOfIntegrationVertices
            allocate( this% ObjectsList(i)% IntegrationVertices(j)% domains(NumOfInterPoints),                 &
                      this% ObjectsList(i)% IntegrationVertices(j)% indeces(NumOfInterPoints),                 &
                      this% ObjectsList(i)% IntegrationVertices(j)% invPhi(NumOfInterPoints,NumOfInterPoints), &
                      this% ObjectsList(i)% IntegrationVertices(j)% b(NumOfInterPoints)                        )
         end do 
      end do 
 
   end subroutine STL_SetIntegration

   function STL_ComputeScalarIntegral( this )

      implicit none 

      class(STLfile), intent(inout) :: this 
      real(kind=RP)                 :: STL_ComputeScalarIntegral

      real(kind=RP) :: ScalarVar(NumOfIntegrationVertices)
      integer       :: i, j

      STL_ComputeScalarIntegral = 0.0_RP

      do i = 1, this% NumOfObjs
         associate( obj => this% ObjectsList(i) )
         do j = 1, NumOfIntegrationVertices
            ScalarVar(j) = obj% IntegrationVertices(j)% ScalarValue
         end do
         STL_ComputeScalarIntegral = STL_ComputeScalarIntegral + TriangleScalarIntegral( obj, ScalarVar )
         end associate
      end do

   end function STL_ComputeScalarIntegral

   function STL_ComputeVectorIntegral( this )

      implicit none 

      class(STLfile), intent(inout) :: this 
      real(kind=RP)                 :: STL_ComputeVectorIntegral(NDIM)

      real(kind=RP) :: VectorVar(NDIM,NumOfIntegrationVertices)
      integer       :: i, j

      STL_ComputeVectorIntegral = 0.0_RP

      do i = 1, this% NumOfObjs
         associate( obj => this% ObjectsList(i) )
         do j = 1, NumOfIntegrationVertices
            VectorVar(:,j) = obj% IntegrationVertices(j)% VectorValue
         end do
         STL_ComputeVectorIntegral = STL_ComputeVectorIntegral + TriangleVectorIntegral( obj, VectorVar )
         end associate
      end do

   end function STL_ComputeVectorIntegral

   function TriangleScalarIntegral( obj, ScalarVar ) result( Val )
      use MappedGeometryClass
      implicit none 

      type(object_type), intent(in) :: obj 
      real(kind=RP),     intent(in) :: ScalarVar(NumOfIntegrationVertices)
      real(kind=RP)                 :: Val

      real(kind=RP) :: AB(NDIM), AC(NDIM), S(NDIM), A

      AB = obj% IntegrationVertices(2)% coords - obj% IntegrationVertices(1)% coords 
      AC = obj% IntegrationVertices(3)% coords - obj% IntegrationVertices(1)% coords 

      call vcross(AB,AC,S)
   
      A = 0.5_RP * norm2(S)

      Val = A * (ScalarVar(1) + ScalarVar(2) + ScalarVar(3))/3.0_RP
      ! Val = A/60.0_RP * ( 27.0_RP * ScalarVar(NumOfIntegrationVertices) +                      &
      !                     3.0_RP  * sum(ScalarVar(1:NumOfVertices))     +                      &
      !                     8.0_RP  * sum(ScalarVar(NumOfVertices+1:NumOfIntegrationVertices-1)) )

   end function TriangleScalarIntegral

   function TriangleVectorIntegral( obj, VectorVar ) result( Val )
      use MappedGeometryClass
      implicit none 

      type(object_type), intent(in) :: obj 
      real(kind=RP),     intent(in) :: VectorVar(NDIM,NumOfIntegrationVertices)
      real(kind=RP)                 :: Val(NDIM)
 
      real(kind=RP) :: AB(NDIM), AC(NDIM), S(NDIM), A
      real(kind=RP) :: normal(NDIM)

      AB = obj% IntegrationVertices(2)% coords - obj% IntegrationVertices(1)% coords 
      AC = obj% IntegrationVertices(3)% coords - obj% IntegrationVertices(1)% coords 

      call vcross(AB,AC,S)

      A = 0.5_RP * norm2(S)
      
      Val = A * (VectorVar(:,1) + VectorVar(:,2) + VectorVar(:,3))/3.0_RP
      ! Val = A/60.0_RP * ( 27.0_RP * VectorVar(:,NumOfIntegrationVertices) +                      &
      !                      3.0_RP * sum(VectorVar(:,1:NumOfVertices))     +                      &
      !                      8.0_RP * sum(VectorVar(:,NumOfVertices+1:NumOfIntegrationVertices-1)) )

   end function TriangleVectorIntegral

   subroutine STL_ResetIntegrationPoints( this )
   
      implicit none 
   
      class(STLfile), intent(inout) :: this

      integer :: i, j 
      
      do i = 1, this% NumOfObjs
         do j = 1, NumOfIntegrationVertices
            this% ObjectsList(i)% IntegrationVertices(j)% Q   = 0.0_RP
            this% ObjectsList(i)% IntegrationVertices(j)% U_x = 0.0_RP
            this% ObjectsList(i)% IntegrationVertices(j)% U_y = 0.0_RP
            this% ObjectsList(i)% IntegrationVertices(j)% U_z = 0.0_RP
         end do 
      end do
   
   end subroutine STL_ResetIntegrationPoints
   
   integer function Point_wrt_Plane( plane_normal, plane_point, point ) result( PointIs )
      use MappedGeometryClass
      implicit none
      !-arguments-----------------------------------------------------------------
      real(kind=RP), dimension(:), intent(in) :: plane_normal, plane_point, point
      !-local-variables-----------------------------------------------------------
      real(kind=RP) :: d, eps
   
      d = dot_product(plane_normal,point) - dot_product(plane_normal,plane_point)
      
      if( d > EPSILON(eps) ) then
         PointIs = POINT_IN_FRONT_PLANE
      elseif( d < -EPSILON(eps) ) then
         PointIs = POINT_BEHIND_PLANE
      else 
         PointIs = POINT_ON_PLANE
      end if
   
   end function Point_wrt_Plane
  
   function EdgePlaneIntersection( plane_normal, plane_point, PointA, PointB ) result( Point_inters )
      use MappedGeometryClass
      implicit none
      !-arguments-----------------------------------------------------------------
      real(kind=RP), intent(in) :: plane_normal(:), plane_point(:), &
                                   PointA(:), PointB(:)
      real(kind=RP)             :: Point_inters(NDIM)
      !-local-variables-----------------------------------------------------------
      real(kind=RP) :: B_A(NDIM), d, t
 
      B_A = PointB - PointA
      
      d = dot_product(plane_normal,plane_point)
      
      t = ( d - dot_product(plane_normal,PointA) )/dot_product(plane_normal,B_A)
 
      Point_inters = PointA + t*B_A 
 
   end function EdgePlaneIntersection
!//////////////////////////////////////////////
   
   subroutine TecFileHeader( FileName, Title, NumOfObjs, funit )
   
      implicit none
      !-arguments------------------------------------------------------
      character(len=*), intent(in)  :: FileName
      integer,          intent(in)  :: NumOfObjs
      character(len=*), intent(in)  :: Title
      integer,          intent(out) :: funit
   
      character(LINE_LENGTH) :: itoa 

      funit = UnusedUnit()
   
      open(funit,file=trim(FileName)//'.tec', status='unknown')

      write(itoa, '(I0)') NumOfObjs
      
      write(funit,"(A9,A,A)") 'TITLE = "',trim(Title),'"'
      write(funit,"(A23)") 'VARIABLES = "x","y","z"'
      write(funit, '(A)') 'ZONE T="Points", I=1, J=1, K=' // trim(adjustl(itoa)) // ', F=POINT'
      
   end subroutine TecFileHeader

   subroutine ReadTecFileHeader( FileName, NumOfObjs, funit )
   
      implicit none
      !-arguments------------------------------------------------------
      character(len=*), intent(in)  :: FileName
      integer,          intent(out) :: NumOfObjs
      integer,          intent(out) :: funit 
   
      integer                    :: ios
      character(len=LINE_LENGTH) :: line, zoneinfo
      
      funit = UnusedUnit()
      
      open(funit,file='IBM/'//trim(FileName)//'.tec', status='old', action='read')

      read(funit, '(A)', iostat=ios) line
      if (ios /= 0) stop 'Error reading title line'

      read(funit, '(A)', iostat=ios) line
      if (ios /= 0) stop 'Error reading variables line'

      read(funit, '(A)', iostat=ios) zoneinfo
      if (ios /= 0) stop 'Error reading zone line'

      ! Extract the number of points (NumOfObjs) from the zoneinfo line
      call get_num_of_objs(trim(zoneinfo), NumOfObjs)
      
   end subroutine ReadTecFileHeader

   subroutine get_num_of_objs(zoneinfo, NumOfObjs)

      implicit none 

      character(len=*), intent(in) :: zoneinfo
      integer,          intent(out) :: NumOfObjs

      character(len=LINE_LENGTH) :: temp
      integer                    :: pos, ios

      pos = index(zoneinfo, 'K=')
      if (pos > 0) then
          !read(zoneinfo(pos+2:), '(I)', iostat=ios) NumOfObjs
          read(zoneinfo(pos+2:), *, iostat=ios) NumOfObjs
          if (ios /= 0) stop 'Error reading number of objects'
      else
          stop 'Error: number of objects not found in zone info'
      end if

  end subroutine get_num_of_objs

end module TessellationTypes
