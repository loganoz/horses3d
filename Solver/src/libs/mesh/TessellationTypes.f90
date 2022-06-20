!
!//////////////////////////////////////////////////////
!
!   @File:    TesselationTypes.f90
!
!//////////////////////////////////////////////////////
#include "Includes.h"
module TessellationTypes

   use SMConstants
   use Utilities
   use IntegerDataLinkedList
   use ParamfileRegions                , only: readValueInRegion

   implicit none

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
      
      real(kind=rp), dimension(NDIM) :: coords, ImagePoint_coords, normal
      real(kind=rp)                  :: theta, dist, Rank
      integer                        :: index, element_index, &
                                        Translate = 0, partition, objIndex
      integer,       dimension(NDIM) :: local_Position
      logical                        :: delete = .false., isInsideBody = .false., &
                                        forcingPoint = .false.
!~       real(kind=rp), dimension(NCONS) :: Q_IP
      
      contains
         procedure :: copy => point_type_copy
   
   end type   
   
!
!  **************************************************
!  Main type for a list of points
!  **************************************************
   type ObjectLinkedList
   
      class(object_type), pointer :: head => null()
      
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
      
      type(point_type), dimension(:),   allocatable :: vertices
      real(kind=rp),    dimension(NDIM)             :: normal, tangent, coords
      integer                                       :: index, NumOfVertices
      integer,          dimension(2)                :: partition
      logical                                       :: ComputeIntegrals = .true.

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
      integer                                      :: NumOfObjs, partition, &
                                                      motionAxis, body
      real(kind=RP)                                :: angularVelocity, ds,  &
                                                      Velocity,             & 
                                                      rotationMatrix(NDIM,NDIM)
      logical                                      :: move
      character(len=LINE_LENGTH)                   :: filename, motionType
   
       contains
          procedure :: ReadTesselation
          procedure :: getRotationaMatrix => STLfile_getRotationaMatrix
          procedure :: getDisplacement    => STLfile_getDisplacement
          procedure :: destroy            => STLfile_destroy
          procedure :: Describe           => Describe_STLfile
          procedure :: plot               => STLfile_plot
          procedure :: DescribePartitions => DescribePartitions_STLfile

   end type
   
   
   
   type ObjsDataLinkedList_t
      type(ObjData_t), pointer     :: head => null()
      integer                       :: no_of_entries = 0
      contains
         procedure   :: Add           => ObjsDataLinkedList_Add
         procedure   :: check         => CheckObj
         procedure   :: Destruct      => ObjsDataLinkedList_Destruct
    end type ObjsDataLinkedList_t
   
    type ObjData_t
      integer                     :: value 
      type(ObjData_t), pointer    :: next 
    end type ObjData_t


    type ObjsRealDataLinkedList_t
      class(ObjRealData_t), pointer    :: head => NULL()
      integer                       :: no_of_entries = 0
      contains
         procedure   :: Add      => ObjsRealDataLinkedList_Add
         procedure   :: check    => CheckReal
         procedure   :: Destruct => ObjsRealDataLinkedList_Destruct
    end type ObjsRealDataLinkedList_t
   
   type ObjRealData_t
      real(kind=RP)             :: value
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
   
   character(len=8) :: ROTATION = "rotation"
   character(len=6) :: LINEAR = "linear"

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
         class(ObjsDataLinkedList_t), intent(inout) :: this
         integer,                     intent(in)    :: value
         type(ObjData_t),                 pointer  :: current

         integer :: i

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
         class(ObjsRealDataLinkedList_t), intent(inout) :: this
         real(kind=RP),                   intent(in)    :: value
         type(ObjRealData_t),             pointer  :: current

         integer :: i

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
     
        class(ObjsDataLinkedList_t), intent(inout) :: this
        integer,                     intent(in)    :: value
        
        type(ObjData_t),             pointer       :: current
        integer :: i
        
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
     
        class(ObjsRealDataLinkedList_t), intent(inout) :: this
        real(kind=RP),                     intent(in)  :: value
        
        type(ObjRealData_t),             pointer       :: current
        integer :: i
        
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
         class(ObjsDataLinkedList_t), intent(inout)   :: this
         type(ObjData_t), pointer       :: data, nextdata
         integer     :: i

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
         class(ObjsRealDataLinkedList_t), intent(inout)   :: this
         type(ObjRealData_t), pointer       :: data, nextdata
         integer     :: i

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
   
      class(PointLinkedList) :: this
      type(point_type)       :: point
      
      type(point_type), pointer :: current => null(), currentNext => null()
   
      if( .not. associated(this% head) ) then 
         allocate(this% head)
         call this% head% copy(point)
         this% head% next => this% head
         this% head% prev => this% head
         return
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
         ERROR stop ":: Head of list can't be deleted here!"
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
      
      class(PointLinkedList), intent(inout) :: this
      
      type(point_type), pointer :: data => null(), dataPrev => null()
      
      data => this% head% prev
   
      if( associated(data, this% head) ) then
         ERROR stop ":: Head of list can't be deleted here!"
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
      
      class(PointLinkedList), intent(inout) :: this
      
      class(point_type), pointer :: current, next 
      integer :: i
      
      if( this% NumOfPoints .eq. 0 ) return
      
      current => this% head
      next    => current% next
      
      do i = 1, this% NumOfPoints 
         deallocate(current)
         current => next
         next    => current% next
      end do
      
      this% NumOfPoints = 0
      
   end subroutine PointLinkedList_Destruct
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine copies a point
!  ------------------------------------------------  

   subroutine point_type_copy( this, point )
   
      implicit none
   
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
      
      type(ObjectLinkedList) :: ObjectLinkedList_Construct
      
      ObjectLinkedList_Construct% head => null()
   
   end function ObjectLinkedList_Construct
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine adds an object to an objects list
!  ------------------------------------------------  

   subroutine ObjectLinkedList_Add( this, object )
   
      implicit none
   
      class(ObjectLinkedList) :: this
      type(Object_type)       :: object
      
      type(object_type), pointer :: current => null(), currentNext => null()
   
   
      if( .not. associated(this% head) ) then 
         allocate(this% head)
         call this% head% copy(object)
         this% head% next => this% head
         this% head% prev => this% head
         return
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
      end if
      
      nullify(current, currentNext)
   
   end subroutine ObjectLinkedList_Add     
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
!  This subroutine destroys an objects list
!  ------------------------------------------------ 
   subroutine ObjectLinkedList_Destruct( this )
    
      implicit none
      
      class(ObjectLinkedList), intent(inout) :: this
      
      type(object_type), pointer :: data => null(), dataPrev => null()
      
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
      
      nullify(dataPrev)
      
   end subroutine ObjectLinkedList_Destruct
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine builds an object
!  -----------------------------------------------
   subroutine object_type_build( this, Points, normal, NumOfVertices, index, computeIntegrals )
    
      implicit none
      !-arguments---------------------------
      class(Object_type),            intent(inout) :: this
      real(kind=RP), dimension(:,:), intent(in)    :: Points
      real(kind=RP), dimension(:),   intent(in)    :: normal
      integer,                       intent(in)    :: NumOfVertices, index
      logical,                       intent(in)    :: computeIntegrals
      !-local-variables-----------------------
      integer :: i
      
      if( allocated(this% vertices) ) deallocate(this% vertices)
      allocate(this% vertices(NumOfVertices))
      
      do i = 1, NumOfVertices
         this% vertices(i)% coords = Points(:,i)
      end do
      
      this% index            = index 
      this% normal           = normal 
      this% computeIntegrals = computeIntegrals 
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
      class(Object_type),  intent(inout) :: this
      type(Object_type),   intent(in)    :: Object
      !-local-variables-----------------------
      integer :: i
      
      allocate(this% vertices(Object% NumOfVertices))
      
      do i = 1, Object% NumOfVertices
         call this% vertices(i)% copy(Object% vertices(i))
      end do
      
      this% index            = Object% index 
      this% normal           = Object% normal 
      this% computeIntegrals = Object% computeIntegrals 
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
      class(Object_type),  intent(inout) :: this
      
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
      implicit none
      !-arguments--------------------------------
      class(STLfile),   intent(inout) :: this
      character(len=*), intent(in)    :: filename
      !-local-variables--------------------------

      if( MPI_Process% isRoot) then
       
         write(STD_OUT,'(/)')
         call Section_Header("Reading stl file")
         write(STD_OUT,'(/)')
         call SubSection_Header('Stl file "' // trim(fileName) // '"')
      
         write(STD_OUT,'(30X,A,A32,I10)') "->" , "Number of objects: " , this% NumOfObjs
         write(STD_OUT,'(30X,A,A32,I10)') "->" , "Number of stored points: " , 3*this% NumOfObjs
         if( this% move ) then
            write(STD_OUT,'(30X,A,A32,A10)') "->" , "Motion: " , this% motionType
            if( this% motionType .eq. ROTATION ) then 
               write(STD_OUT,'(30X,A,A32,F10.3,A)') "->" , "Angular Velocity: " , this% angularVelocity/timeref, " rad/s."
            elseif( this% motionType .eq. LINEAR ) then
               write(STD_OUT,'(30X,A,A32,F10.3,A)') "->" , "Translation Velocity: " , this% Velocity*(Lref/timeref), " m/s."
            end if
            write(STD_OUT,'(30X,A,A32,I10)') "->" , "Axis of motion: " , this% motionAxis
         end if
         
      end if
   
   end subroutine Describe_STLfile
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine describes the .stl file
!  -----------------------------------------------
   subroutine DescribePartitions_STLfile( this )
      use Headers
      use MPI_Process_Info
      implicit none
      !-arguments--------------------------------
      class(STLfile), intent(inout) :: this
      !-local-variables--------------------------
      character(len=LINE_LENGTH) :: myString
      
      write(myString,'(i100)') this% partition
      
      if( this% partition .eq. 1 ) then
         write(STD_OUT,'(/)')
         call Section_Header("stl partitions")
         write(STD_OUT,'(/)')
      end if
      
      call SubSection_Header('partition ' // trim(adjustl(myString)))
      write(STD_OUT,'(30X,A,A32,I10)') "->" , "Number of objects: " , this% NumOfObjs
      write(STD_OUT,'(30X,A,A32,I10)') "->" , "Number of stored points: " , 3*this% NumOfObjs
      
   end subroutine DescribePartitions_STLfile
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine reads the .stl file
!  -----------------------------------------------
   subroutine  ReadTesselation( this, filename )
      use PhysicsStorage
      implicit none
      !-arguments--------------------------------
      class(STLfile),   intent(inout) :: this
      character(len=*), intent(in)    :: filename
      !-local-variables---------------------------
      integer              :: i, j, funit, NumOfTri, fileStat, NumOfVertices
      integer*2            :: padding
      real*4, dimension(3) :: norm, vertex
      character(len=80)    :: header
      
      NumOfVertices = 3 
      
      this% partition = 1
      
      funit = UnusedUnit()
      open(unit=funit,file='MESH/'//trim(filename)//'.stl',status='old',access='stream',form='unformatted', iostat=fileStat)
      
      if( fileStat .ne. 0 ) then
         print *, "Read Tesselation: file '",trim(filename),"' not found"
         error stop
      end if
      
      this% filename = filename

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
            Objs(i)% vertices(j)% coords = vertex/Lref
         end do
         read(funit) padding
         Objs(i)% index = i
         Objs(i)% NumOfVertices = NumOfVertices
         Objs(i)% partition = 1
      end do
        
      end associate   
       
      close(unit=funit)
      
      call this% describe( filename )
      
   end subroutine  ReadTesselation

  subroutine STLfile_plot( this, timestep )
      use MPI_Process_Info
      use PhysicsStorage
      implicit none
      !-arguments----------------------------------------------
      class(STLfile), intent(inout) :: this
      integer,        intent(in)    :: timestep
      !-local-variables----------------------------------------
      character(len=LINE_LENGTH)     :: filename, myString
      integer                        :: i, j, funit
      integer*2                      :: padding = 0
      real*4, dimension(3)           :: norm, vertex
      character(len=80)              :: header = repeat(' ',80)
      
      if( .not.  MPI_Process% isRoot ) return
      
      funit = UnusedUnit()
      
      write(myString,'(i100)') timestep
      
      filename = trim(this% filename)//'_'//trim(adjustl(myString))
      
      open(funit,file='MESH/'//trim(filename)//'.stl', status='unknown',access='stream',form='unformatted')
 
      write(funit) header, this% NumOfObjs 
      
      do i = 1, this% NumOfObjs
         norm = this% ObjectsList(i)% normal
         write(funit) norm(1), norm(2), norm(3)
         do j = 1, size(this% ObjectsList(i)% vertices)
            vertex = this% ObjectsList(i)% vertices(j)% coords * Lref
            write(funit) vertex(1), vertex(2), vertex(3)
         end do
         write(funit) padding      
      end do

      close(funit)
   
   end subroutine STLfile_plot 
  
   subroutine STLfile_GetMotionInfo( this, STLfilename, NumOfSTL )
      use FileReadingUtilities
      use FTValueDictionaryClass
      use PhysicsStorage
      implicit none
      !-arguments------------------------------------------
      type(STLfile),           intent(inout) :: this
      character(len=*),        intent(in)    :: STLfilename
      integer,                 intent(in)    :: NumOfSTL
      !-local-arguments------------------------------------
      integer                    :: i
      integer,       allocatable :: motionAxis_STL
      real(kind=RP), allocatable :: angularVelocity_STL, Velocity_STL
      character(len=LINE_LENGTH) :: in_label, paramFile, &
                                    motion_STL, STL_name 
   
      this% move = .false.
      
      do i = 1, NumOfSTL
      
         write(in_label , '(A,I0)') "#define stl motion ", i

         call get_command_argument(1, paramFile)
         call readValueInRegion ( trim ( paramFile )  , "stl name"         , STL_name           , in_label, "#end" ) 
         call readValueInRegion ( trim ( paramFile )  , "type"             , motion_STL         , in_label, "#end" ) 
         call readValueInRegion ( trim ( paramFile )  , "angular velocity" , angularVelocity_STL, in_label, "#end" ) 
         call readValueInRegion ( trim ( paramFile )  , "velocity"         , Velocity_STL       , in_label, "#end" ) 
         call readValueInRegion ( trim ( paramFile )  , "motion axis"      , motionAxis_STL     , in_label, "#end" ) 

         if( trim(STLfilename) .ne. trim(STL_name) ) cycle

         this% move = .true.

         select case( trim(motion_STL) )
           case( "rotation" )
               this% motionType = ROTATION
            case( "linear" )
               this% motionType = LINEAR
            case default
               print *, "STLfile_GetMotionInfo: motion not recognized. Available motions are ", ROTATION,"and",LINEAR,"."
               error stop
         end select
         
         if( allocated(angularVelocity_STL) ) then
            this% angularVelocity = angularVelocity_STL*timeref
         elseif( this% motionType .eq. ROTATION ) then
            print *, "STLfile_GetMotionInfo: 'angular velocity' must be specified for ", ROTATION, " motion."
            error stop            
         end if
         
         if( allocated(Velocity_STL) ) then
            this% Velocity = Velocity_STL/(Lref/timeref)
         elseif( this% motionType .eq. LINEAR ) then
            print *, "STLfile_GetMotionInfo: 'velocity' must be specified for ", LINEAR, " motion."
            error stop            
         end if
         
         if( allocated(motionAxis_STL) ) then
            this% motionAxis = motionAxis_STL
            if( this% motionAxis .gt. 3 .or. this% motionAxis .lt. 1 ) then
               print *, "STLfile_GetMotionInfo: 'motion axis' =", this% motionAxis, " not valid:"
               print *, "select 1 for x-axis, 2 for y-axis or 3 for z-axis."
               error stop                
            end if
         elseif( this% move ) then
            print *, "STLfile_GetMotionInfo: 'motion axis' must be specified."
            error stop            
         end if

         return

      end do
 
   end subroutine STLfile_GetMotionInfo
   
   subroutine STLfile_getRotationaMatrix( this, dt )
   
      implicit none
      !-arguments-----------------------------
      class(STLfile), intent(inout):: this
      real(kind=RP),  intent(in)   :: dt
      !-local-variables-----------------------
      real(kind=RP) :: theta
   
      theta = this% angularVelocity * dt
   
      this% rotationMatrix = 0.0_RP
   
      select case( this% motionAxis )
         case( IX )
            this% rotationMatrix(1,1) = 1.0_RP
            this% rotationMatrix(2,2) = cos(theta)
            this% rotationMatrix(2,3) = -sin(theta)
            this% rotationMatrix(3,2) = sin(theta)
            this% rotationMatrix(3,3) = cos(theta)
         case( IY ) 
            this% rotationMatrix(2,2) = 1.0_RP
            this% rotationMatrix(1,1) = cos(theta)
            this% rotationMatrix(1,3) = sin(theta)
            this% rotationMatrix(3,1) = -sin(theta)
            this% rotationMatrix(3,3) = cos(theta)
         case( IZ )
            this% rotationMatrix(3,3) = 1.0_RP
            this% rotationMatrix(1,1) = cos(theta)
            this% rotationMatrix(1,2) = -sin(theta)
            this% rotationMatrix(2,1) = sin(theta)
            this% rotationMatrix(2,2) = cos(theta)
      end select
        
   end subroutine STLfile_getRotationaMatrix
   
   subroutine STLfile_getDisplacement( this, dt )
   
      implicit none
      !-arguments-----------------------------
      class(STLfile), intent(inout):: this
      real(kind=RP),  intent(in)   :: dt
   
      this% ds = this% Velocity * dt
        
   end subroutine STLfile_getDisplacement
   
   subroutine STLfile_destroy( this )
   
      implicit none
      !-arguments------------------------------
      class(STLfile), intent(inout) :: this
      !-local-variables------------------------
      integer :: i
      
      do i = 1, this% NumOfObjs
         deallocate(this% ObjectsList(i)% vertices)
      end do
      
      deallocate(this% ObjectsList)  
   
   end subroutine STLfile_destroy
   
   
!//////////////////////////////////////////////
   
   subroutine TecFileHeader( FileName, Title, I, J, K, funit, DATAPACKING, ZONETYPE )
   
      implicit none
      
      character(len=*), intent(in)  :: FileName, Title, DATAPACKING
      integer,          intent(in)  :: I, J, K
      character(len=*), optional    :: ZONETYPE
      integer,          intent(out) :: funit
   
      funit = UnusedUnit()
   
      open(funit,file=trim(FileName)//'.tec', status='unknown')
      
      write(funit,"(A9,A,A)") 'TITLE = "',trim(Title),'"'
      write(funit,"(A23)") 'VARIABLES = "x","y","z"'
      if( present(ZONETYPE) ) then
         write(funit,"(A7,I0,A3,I0,A3,I0,A14,A,A11,A)") 'ZONE I=',I,',J=',J,',K=',K,', DATAPACKING=',trim(DATAPACKING),', ZONETYPE=',trim(ZONETYPE)
      else
         write(funit,"(A7,I0,A3,I0,A3,I0,A14,A)") 'ZONE I=',I,',J=',J,',K=',K,', DATAPACKING=',trim(DATAPACKING)
      end if
      
   end subroutine TecFileHeader
    
end module TessellationTypes
