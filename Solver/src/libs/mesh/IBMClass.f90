#include "Includes.h"
module IBMClass

   use SMConstants
   use Utilities
   use FTValueDictionaryClass
   use ElementClass
   use FaceClass
   use TessellationTypes
   use OrientedBoundingBox
   use KDClass
   use IntegerDataLinkedList
   use MPI_IBMUtilities
#ifdef _HAS_MPI_
   use mpi
#endif
   implicit none

   type :: Integral_Obj
   
      integer,       dimension(:,:), allocatable :: PointsIndex
      real(kind=rp), dimension(:,:), allocatable :: x
      real(kind=rp)                              :: Area
   
   end type 

   type :: Integral_t
   
      real(kind=rp),      dimension(:),   allocatable :: xi, eta, weights
      integer                                         :: Order, n_of_Q_points
      type(Integral_Obj), dimension(:),   allocatable :: IntegObjs
      logical                                         :: ListComputed = .false., constructed = .false.
    
      contains
         procedure :: construct => Integral_construct
         procedure :: destruct  => Integral_destruct 
         procedure :: GetCoords => Integral_GetCoords
   
   end type

   type IBM_type

      type(STLfile),              allocatable :: stl(:)
      type(KDtree)                            :: rootPoints
      type(KDtree),               allocatable :: root(:)
      type(PointLinkedList)                   :: BandPoints
      type(MPI_M_Points_type)                 :: BandRegion
      type(Integral_t),           allocatable :: Integral(:)
      character(len=LINE_LENGTH), allocatable :: STLfilename(:)
      character(len=LINE_LENGTH)              :: filename
      logical                                 :: plotOBB, plotKDtree, active = .false., TimePenal, &
                                                 describeIBM, semiImplicit, active_semiImplicit,   &
                                                 symmetry, plotBandPoints, plotMask,               &
                                                 ComputeInterpolation = .false.,                   &
                                                 Wallfunction = .false.
      real(kind=rp)                           :: eta, BandRegionCoeff, IP_Distance = 0.0_RP,       &
                                                 y_plus_target
      real(kind=rp),              allocatable :: symCoords(:), penalization(:)
      integer                                 :: KDtree_Min_n_of_Objs, KDtree_n_of_interPoints,   &
                                                 IntegrationOrder, n_of_INpoints,                 &
                                                 rank, lvl = 0, NumOfSTL, NumOfForcingPoints
      integer,                    allocatable :: symPlanes(:), ImagePoint_NearestPoints(:,:)
      
      contains  
         procedure :: read_info                 => IBM_read_info
         procedure :: construct                 => IBM_construct
         procedure :: constructMask             => IBM_constructMask
         procedure :: constructSTL_KDtree       => IBM_constructSTL_KDtree
         procedure :: CheckPoint                => IBM_checkPoint
         procedure :: constructBandRegion       => IBM_constructBandRegion
         procedure :: build                     => IBM_build
         procedure :: GetMask                   => IBM_GetMask
         procedure :: MPI_send2Root             => IBM_MPI_send2Root
         procedure :: MPI_send2Partitions       => IBM_MPI_send2Partitions
         procedure :: MPI_sendMask2Root         => IBM_MPI_sendMask2Root
         procedure :: MPI_sendMask2Partitions   => IBM_MPI_sendMask2Partitions
         procedure :: MPI_sendBand2Root         => IBM_MPI_sendBand2Root
         procedure :: MPI_sendBand2Partitions   => IBM_MPI_sendBand2Partitions
         procedure :: BandRegionPoints          => IBM_bandRegionPoints
         procedure :: GetForcingPointsGeom      => IBM_GetForcingPointsGeom
         procedure :: GetInfo                   => IBM_GetInfo 
         procedure :: SourceTerm                => IBM_SourceTerm
         procedure :: ComputeIBMWallDistance    => IBM_ComputeIBMWallDistance
         procedure :: SemiImplicitCorrection    => IBM_SemiImplicitCorrection
         procedure :: GetImagePoint_nearest     => IBM_GetImagePoint_nearest
#if defined(NAVIERSTOKES) 
         procedure :: SourceTermTurbulence      => IBM_SourceTermTurbulence
#endif
         procedure :: semiImplicitShiftJacobian => IBM_semiImplicitShiftJacobian
         procedure :: semiImplicitJacobian      => IBM_semiImplicitJacobian
         procedure :: GetSemiImplicitStep       => IBM_GetSemiImplicitStep
         procedure :: upDateNormals             => IBM_upDateNormals
         procedure :: SetIntegration            => IBM_SetIntegration
         procedure :: copyKDtree                => IBM_copyKDtree
         procedure :: MoveBody                  => IBM_MoveBody
         procedure :: CleanMask                 => IBM_CleanMask
         procedure :: BandPoint_state           => IBM_BandPoint_state
         procedure :: Describe                  => IBM_Describe
         procedure :: plot_Mask                 => IBM_plot_Mask
         procedure :: Destruct                  => IBM_Destruct
         procedure :: DestroyKDtree             => IBM_DestroyKDtree

   end type
   
   integer, parameter :: mean_density  = 1
   integer, parameter :: mean_u_II     = 2
   integer, parameter :: mean_mu       = 3

   contains
   
   subroutine Integral_construct( this, Order, NumOfObjs )
   
      implicit none
   
      class(Integral_t), intent(inout) :: this
      integer,           intent(in)    :: Order, NumOfObjs
      
      if( this% constructed ) return
      
      this% Order = Order
      
      allocate(this% IntegObjs(NumOfObjs))
      
      select case( Order )
         case( 2 )
            this% n_of_Q_points = 3
            
            allocate(this% xi(3),this% eta(3), this% weights(3))
         
            this% xi(1) = 1.0_RP/6.0_RP
            this% xi(2) = 2.0_RP/3.0_RP
            this% xi(3) = 1.0_RP/6.0_RP
         
            this% eta(1) = 1.0_RP/6.0_RP
            this% eta(2) = 1.0_RP/6.0_RP
            this% eta(3) = 2.0_RP/3.0_RP
         
            this% weights(1) = 1.0_RP/6.0_RP 
            this% weights(2) = this% weights(1) 
            this% weights(3) = this% weights(1)
      
         case( 3 )
            this% n_of_Q_points = 4
         
            allocate(this% xi(4), this% eta(4), this% weights(4))
            
            this% xi(1) = 1.0_RP/3.0_RP
            this% xi(2) = 1.0_RP/5.0_RP
            this% xi(3) = 1.0_RP/5.0_RP
            this% xi(4) = 3.0_RP/5.0_RP
            
            this% eta(1) = 1.0_RP/3.0_RP
            this% eta(2) = 1.0_RP/5.0_RP
            this% eta(3) = 3.0_RP/5.0_RP
            this% eta(4) = 1.0_RP/5.0_RP
         
            this% weights(1) = -27.0_RP/96.0_RP
            this% weights(2) = 25.0_RP/96.0_RP
            this% weights(3) = this% weights(2)
            this% weights(4) = this% weights(2)
         
         case default 
            write(*,"(a26,I10,a21)")" Intgral_construct:: order", Order ," not implemented yet."
            write(*,"(a25)")" Available orders:: 2 & 3"
            error stop
      end select
   
      this% constructed = .true.
    
   end subroutine Integral_construct
   
   subroutine Integral_destruct( this )
   
      implicit none
      !-arguments------------------------------
      class(Integral_t), intent(inout) :: this
      !-local-variables------------------------
      integer :: i
      
      deallocate(this% xi) 
      deallocate(this% eta)
      deallocate(this% weights)
      
      do i = 1, size(this% IntegObjs)
         if( allocated(this% IntegObjs(i)% PointsIndex) ) then
            deallocate(this% IntegObjs(i)% PointsIndex)
            deallocate(this% IntegObjs(i)% x)
         end if
      end do
      deallocate(this% IntegObjs)
   
   end subroutine Integral_destruct
   
   subroutine Integral_GetCoords( this, vertex1, vertex2, vertex3, index, n_of_Q_points )
   
      implicit none
      !-arguments---------------------------------------------------------------------
      class(Integral_t),                  intent(inout) :: this
      real(kind=rp),     dimension(NDIM), intent(in)    :: vertex1, vertex2, vertex3
      integer,                            intent(in)    :: index, n_of_Q_points
      !-local-variables---------------------------------------------------------------
      integer :: i
   
      do i = 1, n_of_Q_points
         this% IntegObjs(index)% x(:,i) = vertex1 + ( vertex2-vertex1 ) * this% xi(i) + &
                                                    ( vertex3-vertex1 ) * this% eta(i) 
      end do
     
   end subroutine Integral_GetCoords
   
   subroutine IBM_read_info( this, controlVariables )
      use FileReadingUtilities
      implicit none
      
      class(IBM_type), intent(inout) :: this
      class(FTValueDictionary)       :: controlVariables
 
      call this% GetInfo( controlVariables )

   end subroutine IBM_read_info
   
   
   subroutine IBM_construct( this, controlVariables )   
      use mainKeywordsModule
      use FileReadingUtilities
      use MPI_Process_Info
      implicit none
      !-arguments----------------------------------------
      class(IBM_type), intent(inout) :: this  
      class(FTValueDictionary)       :: controlVariables
      !-local-variables----------------------------------
      character(len=LINE_LENGTH) :: filename, MyString
      integer                    :: STLNum
      
      if( this% describeIBM ) call this% describe()         

      allocate( this% stl(this% NumOfSTL),        &
                OBB(this% NumOfSTL),              &
                this% root(this% NumOfSTL),       &
                this% integral(this% NumOfSTL),   &
                this% STLfilename(this% NumOfSTL) )
 
      do STLNum = 1, this% NumOfSTL
         write(MyString, '(i100)') STLNum
         if( STLNum .eq. 1 ) then
            filename = stlFileNameKey
         else
            filename = trim(stlFileNameKey)//trim(adjustl(MyString))
         end if         
         this% STLfilename(STLNum) = controlVariables% stringValueForKey(trim(filename), requestedLength = LINE_LENGTH)
         call STLfile_GetMotionInfo( this% stl(STLNum), this% STLfilename(STLNum), this% NumOfSTL )
         this% stl(STLNum)% body = STLNum       
         call this% stl(STLNum)% ReadTesselation( this% STLfilename(STLNum) ) 
         call OBB(STLNum)% construct( this% stl(STLNum), this% plotOBB )
         this% root(STLNum)% STLNum = STLNum
         call OBB(STLNum)% ChangeObjsRefFrame( this% stl(STLNum)% ObjectsList )
         call this% constructSTL_KDtree( STLNum )    
      end do

   end subroutine IBM_Construct
   
   subroutine IBM_constructSTL_KDtree( this, STLNum )
      use MPI_Process_Info
      implicit none
      
      class(IBM_type), intent(inout) :: this
      integer,         intent(in)    :: STLNum

      call MPI_KDtree_buildPartition( this% stl(STLNum) )

      if ( MPI_Process% doMPIAction ) then   
         call recvSTLPartition()
      end if

      if( MPI_Process% doMPIRootAction ) then
         call SendSTLPartitions()
      end if

      call this% root(STLNum)% construct( stl           = MPI_KDtreePartition% stl,      &
                                          Vertices      = MPI_KDtreePartition% Vertices, &
                                          isPlot        = this% plotKDtree,              & 
                                          Min_n_of_Objs = this% KDtree_Min_n_of_Objs     )

      this% root(STLNum)% MaxAxis = MPI_KDtreePartition% axis

      call this% Integral(STLNum)% construct( this% IntegrationOrder, MPI_KDtreePartition% stl% NumOfObjs )

      call MPI_KDtree_destroy()

   end subroutine IBM_constructSTL_KDtree
   
   subroutine IBM_Destruct( this, isChild )
      use MPI_Process_Info
      implicit none
      !-arguments------------------------------
      class(IBM_type), intent(inout) :: this
      logical,         intent(in)    :: isChild
      !----------------------------------------
      integer :: i

      do i = 1, size(this% root)
         call this% root(i)% destruct( isChild )
         if( .not. isChild ) call this% Integral(i)% destruct()
         call this% BandPoints% destruct()
         !-band region--
         call this% BandRegion% destroy()
         call this% rootPoints% Destruct( .false. )
      end do

      deallocate(this% penalization)

      if( .not. isChild ) then
          if( allocated(this% symCoords) ) deallocate(this% symCoords)
          if( allocated(this% symPlanes) ) deallocate(this% symPlanes)
      end if
   
   end subroutine IBM_Destruct   
   
   subroutine IBM_DestroyKDtree( this, isChild )
      use MPI_Process_Info
      implicit none
      !-arguments------------------------------
      class(IBM_type), intent(inout) :: this
      logical,         intent(in)    :: isChild
      !----------------------------------------
      integer :: i
        
      do i = 1, size(this% root)
         call this% root(i)% destruct( isChild )
      end do
   
   end subroutine IBM_DestroyKDtree   
   
   
   subroutine IBM_build( this, elements, no_of_elements, no_of_DoFs, isChild )
      use MPI_Process_Info
      implicit none
      !-arguments-----------------------------------------------------------------
      class(IBM_type),               intent(inout) :: this
      type(element),   dimension(:), intent(inout) :: elements
      integer,                       intent(in)    :: no_of_elements, no_of_DoFs
      logical,                       intent(in)    :: isChild       
      !-local-variables-----------------------------------------------------------
      integer :: MaskPoints, STLNum
#ifdef _HAS_MPI_
      integer :: localVal, ierr
#endif

      this% n_of_INpoints = 0   

      do STLNum = 1, this% NumOfSTL
         call this% GetMask( elements, no_of_elements, no_of_DoFs, STLNum )
      end do
     
      if( MPI_Process% doMPIAction ) then
#ifdef _HAS_MPI_
         localVal = this% n_of_INpoints
         call mpi_allreduce(localVal, MaskPoints, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
      else
         MaskPoints = this% n_of_INpoints
      end if   
     
      if( MaskPoints .eq. 0 .and. this% lvl .gt. 0 ) then
         print *, "The mask for the multigrid level ", this% lvl, " is made of 0 points."
         print *, "Try to increase the polynomial order or to refine the mesh."
         error stop
      elseif( MaskPoints .eq. 0 ) then
         print *, "The mask is made of 0 points."
         print *, "Try to increase the polynomial order or to refine the mesh."
         error stop         
      end if
      
      allocate( this% penalization(no_of_elements) )
      
      this% penalization = this% eta
      
      if( this% plotMask ) call this% plot_Mask( elements, no_of_elements )
   
      if( isChild .and. .not. this% Wallfunction ) return
      
      call this% constructBandRegion( elements, no_of_elements ) 

      if( .not. isChild ) then
         do STLNum = 1, this% NumOfSTL
            call this% SetIntegration( STLNum )
         end do 
      end if

   end subroutine IBM_build
   
   subroutine IBM_GetMask( this, elements, no_of_elements, no_of_DoFs, STLNum )
      use MPI_Process_Info
      implicit none
      !-arguments------------------------------------------------------
      class(IBM_type),               intent(inout) :: this
      type(element),   dimension(:), intent(inout) :: elements
      integer,                       intent(in)    :: no_of_elements, no_of_DoFs, STLNum

      call this% MPI_send2Root( elements, no_of_elements, no_of_DoFs, STLNum )

      call this% MPI_send2Partitions() 

      call this% constructmask( MPI_M_Points_ALL% x, elements, STLNum )

      call MPI_M_PointsPartition% destroy()

      call MPI_M_Points_ALL% destroy()

   end subroutine IBM_GetMask
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
   subroutine IBM_MPI_send2Root( this, elements, no_of_elements, no_of_DoFs, STLNum )
      use MPI_Process_Info
      implicit none
   
      class(IBM_type),               intent(inout) :: this
      type(element),   dimension(:), intent(inout) :: elements
      integer,                       intent(in)    :: no_of_elements, no_of_DoFs, STLNum
#ifdef _HAS_MPI_      
      integer :: LocalVal, ierr
#endif

      call MaskCandidates( elements, no_of_elements, no_of_DoFs, STLNum, this% NumOfSTL )

      if( MPI_Process% doMPIAction ) then
#ifdef _HAS_MPI_
         localVal = MPI_M_PointsPartition% LocNumOfObjs
         call mpi_allreduce(localVal, MPI_M_Points_All% NumOfObjs, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
      else
         MPI_M_Points_All% NumOfObjs = MPI_M_PointsPartition% LocNumOfObjs
      end if

      call MPI_Pointpartition(MPI_M_Points_All% NumOfObjs)

      if ( MPI_Process% doMPIAction ) then
         call RootRecvrecvPointMaskPartition()
         call RootSendPointMaskPartition()
      end if

   end subroutine IBM_MPI_send2Root
  
  
   subroutine IBM_MPI_send2Partitions( this )
      use MPI_Process_Info
      implicit none
   
      class(IBM_type), intent(inout) :: this

      if ( MPI_Process% doMPIAction ) then
         call recvPointMaskPartition()
      end if 
      
      if( MPI_Process% doMPIRootAction ) then
         call sendPointMaskPartition()
      end if

   end subroutine IBM_MPI_send2Partitions
  
!--------------------------------------
   
   subroutine IBM_MPI_sendMask2Root( this )
      use MPI_Process_Info
      implicit none
   
      class(IBM_type), intent(inout) :: this
   
      if ( MPI_Process% doMPIAction ) then
         call RootRecvPointMask()
         call RootSendPointMask()
      end if
      
   end subroutine IBM_MPI_sendMask2Root
  
  
   subroutine IBM_MPI_sendMask2Partitions( this )
      use MPI_Process_Info
      implicit none
   
      class(IBM_type), intent(inout) :: this
 
      if ( MPI_Process% doMPIAction ) then
         call recvPointMask()
      end if

      if( MPI_Process% doMPIRootAction ) then
         call sendPointMask()
      end if

   end subroutine IBM_MPI_sendMask2Partitions
 
!--------------------------------------
   
   subroutine IBM_MPI_sendBand2Root( this )
      use MPI_Process_Info
      implicit none
   
      class(IBM_type), intent(inout) :: this
   
      if ( MPI_Process% doMPIAction ) then
         call RootrecvBandPoint( this% BandRegion )
         call RootSendBandPoint( this% BandPoints )
      end if
      
   end subroutine IBM_MPI_sendBand2Root
  
  
   subroutine IBM_MPI_sendBand2Partitions( this )
      use MPI_Process_Info
      implicit none
   
      class(IBM_type), intent(inout) :: this
 
      if ( MPI_Process% doMPIAction ) then
         call recvBandPointPartition( this% BandRegion )
      end if 

      if( MPI_Process% doMPIRootAction ) then
        call sendBandPointPartition( this% BandRegion )
      end if

   end subroutine IBM_MPI_sendBand2Partitions  
  
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -----------------------------------------------------------
! This subroutine plots the points outside & inside the body
!  -----------------------------------------------------------

   subroutine IBM_plot_Mask(  this, elements,  no_of_elements, timestep )
      use MPI_Process_Info
      implicit none
      !-arguments------------------------------------------------------
      class(IBM_type),               intent(inout) :: this
      type(element),   dimension(:), intent(inout) :: elements
      integer,                       intent(in)    :: no_of_elements, timestep
      !-local-variables------------------------------------------------
      character(len=LINE_LENGTH) :: PointFile, lvlName, nRank, step
      integer                    :: eID, i, j, k, funit
      logical                    :: add_Point = .false.
      
      optional :: timestep
      
      if( MPI_Process% nProcs .gt. 1 ) then
         write(nRank,"(I100)") MPI_Process% rank
         if( this% lvl .gt. 0 ) then
            write(lvlName,"(I10)") this% lvl
            PointFile = trim(this% filename)//'_MGlevel'//trim(adjustl(lvlName))//'_'//trim(adjustl(nRank))
         else
            PointFile = trim(this% filename)//'_'//trim(adjustl(nRank))
         end if
      else
         if( this% lvl .gt. 0 ) then
            write(lvlName,"(I10)") this% lvl
            PointFile = trim(this% filename)//'_MGlevel'//trim(adjustl(lvlName))
         else
            PointFile = trim(this% filename)
         end if      
      end if
      
      if( this% n_of_INpoints .gt. 0 ) then
         if( present(timestep) ) then
            write(step,"(I10)") timestep
            PointFile = PointFile//trim(step)
         end if
        
         call TecFileHeader( 'IBM/Mask_'//trim(PointFile), 'Mask Points', this% n_of_INpoints/2+mod(this% n_of_INpoints,2),2,1, funit, 'POINT')
     
         if( mod(this% n_of_INpoints,2) .ne. 0 )  add_Point  = .true.

         do eID = 1,  no_of_elements
            associate ( e => elements(eID) )
            do k = 0, e% Nxyz(3); do j = 0, e% Nxyz(2) ; do i = 0, e% Nxyz(1)
               if( e% isInsideBody(i,j,k) ) then
                  write(funit,'(3E13.5)')  e% geom% x(1,i,j,k), e% geom% x(2,i,j,k), e% geom% x(3,i,j,k)
                  if( add_Point ) then
                     write(funit,'(3E13.5)')  e% geom% x(1,i,j,k), e% geom% x(2,i,j,k), e% geom% x(3,i,j,k)
                     add_Point = .false.
                  end if
               end if
            end do; end do; end do
            end associate
         end do
         
         close(funit)
      end if     
      
   end subroutine IBM_plot_Mask    
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine constructs the mask
!  -----------------------------------------------   
   subroutine IBM_constructmask( this, x, elements, STLNum )
      use MPI_Process_Info
      implicit none
      !-arguments---------------------------------------------------- 
      class(IBM_type),                intent(inout) :: this
      type(point_type), dimension(:), intent(inout) :: x
      integer,                        intent(in)    :: STLNum
      type(element),    dimension(:), intent(inout) :: elements
      !-local-variables----------------------------------------------
      real(kind=RP) :: Point(NDIM) 
      integer       :: eID, n, i, j, k

!$omp parallel shared(x,STLNum,n) 
!$omp do schedule(runtime) private(Point)
      do n = 1, size(x)           
         call OBB(STLNum)% ChangeRefFrame( x(n)% coords, 'local', Point )
         if( isInsideBox( Point, this% root(STLNum)% vertices ) ) then
            call this% CheckPoint( x(n)% coords, STLNum, x(n)% isInsideBody )   
         end if
      end do
!$omp end do
!$omp end parallel

       call this% MPI_sendMask2Root()       
       call this% MPI_sendMask2Partitions()

!$omp parallel shared(this,x,elements,STLNum,n)
!$omp do schedule(runtime) private(eID,i,j,k)
       do n = 1, size(x)
          if( x(n)% partition .eq. MPI_Process% rank ) then
             eID = x(n)% element_index
             i = x(n)% local_Position(1)
             j = x(n)% local_Position(2)
             k = x(n)% local_Position(3)
             elements(eID)% isInsideBody(i,j,k) = x(n)% isInsideBody
             if( elements(eID)% isInsideBody(i,j,k) ) then 
                 elements(eID)% STL(STLNum,i,j,k) = STLNum
!$omp critical
                this% n_of_INpoints = this% n_of_INpoints + 1   
!$omp end critical 
             end if 
         end if    
      end do
!$omp end do
!$omp end parallel

   end subroutine IBM_constructmask
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  --------------------------------------------------------------
! This subroutine constructs the band region for the integration
!  --------------------------------------------------------------

   subroutine IBM_constructBandRegion( this, elements, no_of_elements )
      use MPI_Process_Info
      implicit none
      !-arguments------------------------------------------------------
      class(IBM_type),               intent(inout) :: this
      type(element),   dimension(:), intent(inout) :: elements
      integer,                       intent(in)    :: no_of_elements
      !-local-variables------------------------------------------------
#ifdef _HAS_MPI_
      integer localVal, ierr
#endif

      call this% BandRegionPoints( elements, no_of_elements )

      if ( MPI_Process % doMPIAction ) then
#ifdef _HAS_MPI_
         localVal = this% BandPoints% NumOfPoints
         call mpi_allreduce(localVal, this% BandRegion% NumOfObjs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
      else
         this% BandRegion% NumOfObjs = this% BandPoints% NumOfPoints
      end if

      if( this% BandRegion% NumOfObjs .eq. 0 ) then
         print *, "IBM_bandRegionPoints: Number of points in the band region is 0"
         error stop
      end if

      call MPI_BandPointpartition( this% BandRegion, this% BandRegion% NumOfObjs, this% BandPoints )

      call this% MPI_sendBand2Root()      
      call this% MPI_sendBand2Partitions() 

      ! STL & OBB not used for rootpoints
      call this% rootPoints% construct( stl           = this% stl(1),                  &
                                        Vertices      = OBB(1)% LocVertices,           &
                                        isPlot        = this% plotBandPoints,          &
                                        Min_n_of_Objs = this% KDtree_n_of_interPoints, &
                                        PointList     = this% BandRegion% x,           &
                                        lvl           = this% lvl                      )  
 
      call this% BandPoints% destruct()

   end subroutine IBM_constructBandRegion     
   
   subroutine IBM_bandRegionPoints( this, elements, n_of_elements )
      use ElementClass
      use FaceClass
      use MPI_Process_Info
      implicit none
      !-arguments------------------------------------------------------------
      class(IBM_type),               intent(inout) :: this
      type(element),   dimension(:), intent(inout) :: elements
      integer,                       intent(in)    :: n_of_elements
      !-local-variables------------------------------------------------------
      type(point_type)               :: p       
      integer                        :: eID, n, i, j, k, NumOfPoints, STLNum

      this% BandPoints = PointLinkedList()      

      n = 0; NumOfPoints = 0

      do eID = 1, n_of_elements
         associate( e => elements(eID) )
         do i = 0, e% Nxyz(1); do j = 0, e% Nxyz(2); do k = 0, e% Nxyz(3) 
            if( .not. e% isInsideBody(i,j,k) ) then
               do STLNum = 1, this% NumOfSTL
                  if( OBB(STLNum)% isPointInside( e% geom% x(:,i,j,k), this% BandRegionCoeff ) ) then
                     p% coords = e% geom% x(:,i,j,k)
                     p% element_index = eID; p% local_Position = (/ i,j,k /)
                     n = n + 1; NumOfPoints = NumOfPoints + 1 
                     p% index = n; p% partition = MPI_Process% rank
                     call this% BandPoints% add(p)
                     exit
                  end if
               end do
            end if
         end do; end do; end do
         end associate
      end do
      
      this% BandPoints% NumOfPoints = NumOfPoints
      
   end subroutine IBM_bandRegionPoints   
   
   subroutine IBM_copyKDtree( this, KDtree2copy )
   
      implicit none
   
      class(IBM_type),         intent(inout) :: this
      type(KDtree),    target, intent(in)    :: KDtree2copy(:)
   
      integer :: i
   
      allocate(this% root(size(KDtree2copy)))
   
      do i = 1, size(KDtree2copy)
         this% root(i) = KDtree2copy(i)
      end do
   
   end subroutine IBM_copyKDtree
   
   subroutine IBM_upDateNormals( this, STLNum )
      use MappedGeometryClass 
      implicit none
      !-argument--------------------------------
      class(IBM_type), intent(inout) :: this
      integer,         intent(in)    :: STLNum
      !-local-variables-------------------------
      real(kind=rp) :: dv(NDIM),         &
                       Vertices(NDIM,3)
      integer       :: i
      
!$omp parallel shared(this,STLNum,i)
!$omp do schedule(runtime) private(Vertices,dv)
      do i = 1, this% root(STLNum)% NumOfObjs
         call OBB(STLNum)% ChangeRefFrame( this% root(STLNum)% ObjectsList(i)% vertices(1)% coords ,'global', Vertices(:,1))
         call OBB(STLNum)% ChangeRefFrame( this% root(STLNum)% ObjectsList(i)% vertices(2)% coords ,'global', Vertices(:,2))
         call OBB(STLNum)% ChangeRefFrame( this% root(STLNum)% ObjectsList(i)% vertices(3)% coords ,'global', Vertices(:,3))
         call vcross(Vertices(:,2) - Vertices(:,1),Vertices(:,3) - Vertices(:,1),dv)
         dv = dv/norm2(dv)
         this% root(STLNum)% ObjectsList(i)% normal = dv/norm2(dv)      
      end do
!$omp end do
!$omp end parallel   
   end subroutine IBM_upDateNormals  
    
   subroutine IBM_GetForcingPointsGeom( this, elements )
      use MPI_Process_Info
      
      use omp_lib
      
      implicit none
      !-arguments---------------------------------------------
      class(IBM_type),               intent(inout) :: this
      type(element),   dimension(:), intent(inout) :: elements
      !-local-variables--------------------------------------
      real(kind=RP) :: Dist, Point(NDIM), IntersectionPoint(NDIM)
      integer       :: eID, i, j, k  
      
      character(len=LINE_LENGTH) :: MyString
      integer :: funit
      real(kind=rp) :: ImagePoint(NDIM)
      
      logical :: isInsideBody
      
      this% IP_Distance = InitializeDistance( this% y_plus_target )

      this% NumOfForcingPoints = 0
      
!$omp parallel 
!$omp do schedule(runtime) private(i,j,k) 
      do eID = 1, size(elements)
         do k = 0, elements(eID)% Nxyz(3); do j = 0, elements(eID)% Nxyz(2); do i = 0, elements(eID)% Nxyz(1)        
            if( elements(eID)% geom% dWall(i,j,k) .gt. this% IP_Distance .or. elements(eID)% isInsideBody(i,j,k) ) cycle
            elements(eID)% isForcingPoint(i,j,k) = .true. 
!$omp critical
            this% NumOfForcingPoints = this% NumOfForcingPoints + 1 
!$omp end critical
         end do                          ; end do                          ; end do
      end do
!$omp end do
!$omp end parallel

      if( this% NumOfForcingPoints .eq. 0 ) then
         print *, 'Number of forcing point is zero, try to increase y+'
         error stop
      end if


    write(MyString,*) this% lvl

      call TecFileHeader( 'IBM/ForcingPoints'//trim(adjustl(MyString)), 'FP Points', this% NumOfForcingPoints/2+mod(this% NumOfForcingPoints,2),2,1, funit, 'POINT')

         do eID = 1, size(elements)
            associate ( e => elements(eID) )
            do k = 0, e% Nxyz(3); do j = 0, e% Nxyz(2) ; do i = 0, e% Nxyz(1)
               if( e% isForcingPoint(i,j,k) ) then
                  write(funit,'(3E13.5)')  e% geom% x(1,i,j,k), e% geom% x(2,i,j,k), e% geom% x(3,i,j,k)
               end if
            end do; end do; end do
            end associate
         end do
         
         close(funit)
         
         
    call TecFileHeader( 'IBM/ImagePoints'//trim(adjustl(MyString)), 'Image Points', this% NumOfForcingPoints/2+mod(this% NumOfForcingPoints,2),2,1, funit, 'POINT')

      do eID = 1, size(elements)
         associate ( e => elements(eID) )
          do k = 0, e% Nxyz(3); do j = 0, e% Nxyz(2) ; do i = 0, e% Nxyz(1)
            if( e% isForcingPoint(i,j,k) ) then
               Dist = this% IP_distance - elements(eID)% geom% dWall(i,j,k)
               ImagePoint = elements(eID)% geom% x(:,i,j,k) + Dist * elements(eID)% geom% normal(:,i,j,k)
               
!~                call this% CheckPoint( ImagePoint, 1, isInsideBody )
               
!~                if( isInsideBody ) then
!~                   write(*,*) 'imagepoint =', ImagePoint
!~                   write(*,*) 'i,j,k =', i,j,k
!~                   write(*,*) 'eID =', eID
!~                end if
               
               write(funit,'(3E13.5)')  ImagePoint(1), ImagePoint(2), ImagePoint(3)
            end if
         end do; end do; end do
         end associate
      end do
         
      close(funit)

   end subroutine IBM_GetForcingPointsGeom
   
   subroutine IBM_BandPoint_state( this, elements, bpQ, bpU_x, bpU_y, bpU_z )
      use PhysicsStorage
      use MPI_Process_Info
      implicit none
      !-arguments-------------------------------------------------------
      class(IBM_type),                 intent(inout) :: this
      type(element),   dimension(:),   intent(in)    :: elements
      real(kind=rp),   dimension(:,:), intent(inout) :: bpQ, bpU_x,  &
                                                        bpU_y, bpU_z
      !-local-variables-------------------------------------------------
      integer                                  :: i
#ifdef _HAS_MPI_
      real(kind=rp), dimension(:), allocatable :: local_sum
      integer                                  :: ierr
#endif

      optional :: bpU_x, bpU_y, bpU_z
      
      bpQ = 0.0_RP
      if( present(bpU_x) ) bpU_x = 0.0_RP
      if( present(bpU_y) ) bpU_y = 0.0_RP
      if( present(bpU_z) ) bpU_z = 0.0_RP
      
!$omp parallel shared(i,bpQ,bpU_x,bpU_y,bpU_z)
!$omp do schedule(runtime)
      do i = 1, this% BandRegion% NumOfObjs      
         if( this% BandRegion% x(i)% partition .eq. MPI_Process% rank ) then
            bpQ(:,i) = elements(this% BandRegion% x(i)% element_index)% storage%  &
                                   Q(:,this% BandRegion% x(i)% local_Position(1), &
                                       this% BandRegion% x(i)% local_Position(2), &
                                       this% BandRegion% x(i)% local_Position(3)  )
            if( present(bpU_x) ) then
               bpU_x(:,i) = elements(this% BandRegion% x(i)% element_index)% storage%  &
                                      U_x(:,this% BandRegion% x(i)% local_Position(1), &
                                            this% BandRegion% x(i)% local_Position(2), &
                                            this% BandRegion% x(i)% local_Position(3)  )
            end if
            if( present(bpU_y) ) then
               bpU_y(:,i) = elements(this% BandRegion% x(i)% element_index)% storage%  &
                                      U_y(:,this% BandRegion% x(i)% local_Position(1), &
                                            this% BandRegion% x(i)% local_Position(2), &
                                            this% BandRegion% x(i)% local_Position(3)  )
            end if
            if( present(bpU_z) ) then
               bpU_z(:,i) = elements(this% BandRegion% x(i)% element_index)% storage%  &
                                      U_z(:,this% BandRegion% x(i)% local_Position(1), &
                                            this% BandRegion% x(i)% local_Position(2), &
                                            this% BandRegion% x(i)% local_Position(3)  )
            end if
         end if
      end do
!$omp end do
!$omp end parallel

#ifdef _HAS_MPI_ 
      allocate(local_sum(this% BandRegion% NumOfObjs))
      do i = 1, NCONS
         local_sum = bpQ(i,:)
         call mpi_allreduce(local_sum, bpQ(i,:), this% BandRegion% NumOfObjs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
         if( present(bpU_x) ) then
            local_sum = bpU_x(i,:)
            call mpi_allreduce(local_sum, bpU_x(i,:), this% BandRegion% NumOfObjs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
         end if
         if( present(bpU_y) ) then
            local_sum = bpU_y(i,:)
            call mpi_allreduce(local_sum, bpU_y(i,:), this% BandRegion% NumOfObjs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
         end if
         if( present(bpU_z) ) then
            local_sum = bpU_z(i,:)
            call mpi_allreduce(local_sum, bpU_z(i,:), this% BandRegion% NumOfObjs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
         end if
      end do 
      deallocate(local_sum)
#endif
   end subroutine IBM_BandPoint_state
   
   
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine reads the info in the #define region of controlfile
!  -----------------------------------------------   
   subroutine IBM_GetInfo( this, controlVariables )
      use FileReadingUtilities
#if defined(NAVIERSTOKES)
      use WallFunctionDefinitions
#endif    
      implicit none
      !-arguments----------------------------------------------------------------
      class(IBM_type), intent(inout) :: this
      class(FTValueDictionary)       :: controlVariables
      !-local-variables----------------------------------------------------------
      logical,       allocatable :: active_in, describe_in, plotOBB_in,     &
                                    plotKDtree_in, semiImplicit_in,         &
                                    plotBandPoints_in, plotMask_in
      real(kind=rp), allocatable :: penalization_in, BandRegionCoeff_in,    &
                                    y_plus_target_in
      integer,       allocatable :: n_of_Objs_in, n_of_interpoints_in,      &
                                    IntegrationOrder_in, sym_planes(:)
      character(len=LINE_LENGTH) :: in_label, paramFile, name_in, tmp
      integer                    :: i
      
      character(len=LINE_LENGTH), parameter :: NumberOfSTL = "number of stl" 
      
      
      logical :: correct !DELETE
      
!     Read block
!     **********
      write(in_label , '(A)') "#define ibm"
      call get_command_argument(1, paramFile) 
      call readValueInRegion ( trim ( paramFile ), "name", name_in, in_label, "#end" )
      call readValueInRegion ( trim ( paramFile ), "active", active_in, in_label, "#end" )
      call readValueInRegion ( trim ( paramFile ), "penalization", penalization_in, in_label, "#end" )   
      call readValueInRegion ( trim ( paramFile ), "semi implicit", semiImplicit_in, in_label, "#end" )
      call readValueInRegion ( trim ( paramFile ), "target y plus", y_plus_target_in, in_label, "#end" )
      call readValueInRegion ( trim ( paramFile ), "number of objects", n_of_Objs_in, in_label, "#end" )
      call readValueInRegion ( trim ( paramFile ), "number of interpolation points", n_of_interpoints_in, in_label, "#end" )
      call readValueInRegion ( trim ( paramFile ), "band region coefficient", BandRegionCoeff_in, in_label, "#end" )
      call readValueInRegion ( trim ( paramFile ), "integration order", IntegrationOrder_in, in_label, "#end" )     
      call readValueInRegion ( trim ( paramFile ), "describe", describe_in, in_label, "#end" )
      call readValueInRegion ( trim ( paramFile ), "plot obb"  ,plotOBB_in, in_label , "# end" )
      call readValueInRegion ( trim ( paramFile ), "plot kdtree"  ,plotKDtree_in, in_label , "# end" )
      call readValueInRegion ( trim ( paramFile ), "plot mask", plotMask_in, in_label, "#end" )
      call readValueInRegion ( trim ( paramFile ), "plot band points", plotBandPoints_in, in_label, "#end" )

      this% filename = trim(name_in)

      if( allocated(active_in) ) then
         this% active = active_in
      else
         this% active = .FALSE.
      end if
      
      if( .not. this% active) return
      
      this% semiImplicit = .false.
      if( allocated(semiImplicit_in) ) then
         this% active_semiImplicit = semiImplicit_in
      else
         this% active_semiImplicit = .FALSE.
      end if   
        
       if( allocated(describe_in) ) then
         this% describeIBM = describe_in
      else
         this% describeIBM = .FALSE.
      end if        
   
      if( allocated(plotOBB_in) ) then
         this% plotOBB = plotOBB_in
      else
         this% plotOBB = .FALSE.
      end if

      if( allocated(plotKDtree_in) ) then
         this% plotKDtree = plotKDtree_in
      else
         this% plotKDtree = .FALSE.
      end if

      if( allocated(penalization_in) ) then
         this% eta = penalization_in
         this% TimePenal    = .false.
      else
         this% eta = 0.1_RP
         this% TimePenal = .true.
      end if
 
      if( allocated(n_of_Objs_in) ) then
         this% KDtree_Min_n_of_Objs = n_of_Objs_in
      else
         this% KDtree_Min_n_of_Objs = 5
      end if
      
      if( allocated(n_of_interpoints_in) ) then
         this% KDtree_n_of_interPoints = n_of_interpoints_in
      else
         this% KDtree_n_of_interPoints = 15
      end if
    
      if( allocated(IntegrationOrder_in) ) then
         this% IntegrationOrder = IntegrationOrder_in
      else
         this% IntegrationOrder = 2
      end if  
      
      if( allocated(plotBandPoints_in) ) then
         this% plotBandPoints = plotBandPoints_in
      else
         this% plotBandPoints = .false.
      end if
      
      if( allocated(plotMask_in) ) then
         this% plotMask = plotMask_in
      else
         this% plotMask = .false.
      end if
      
      if( allocated(BandRegionCoeff_in) ) then
         this% BandRegionCoeff = BandRegionCoeff_in
      else
         this% BandRegionCoeff = 1.5_RP
      end if

#if defined(NAVIERSTOKES)
     if( controlVariables% containsKey("wall function") ) then
        this% Wallfunction = .true.
        call Initialize_Wall_Fuction(controlVariables, correct)   !TO BE REMOVED
        if( allocated(y_plus_target_in) ) then
           this% y_plus_target = y_plus_target_in
        else
           this% y_plus_target = 100.0_RP
        end if
     else
        this% Wallfunction = .false.
     end if
#endif
   
      if( controlVariables% containsKey(trim(NumberOfSTL)) ) then
         tmp = controlVariables% StringValueForKey(trim(NumberOfSTL),LINE_LENGTH) 
         this% NumOfSTL = GetIntValue(tmp)
      else
         this% NumOfSTL = 1
      end if      

      this% symmetry = .false.
      
      if( controlVariables% containsKey("symmetry planes") ) then
         this% symmetry = .true.
         tmp        = controlVariables% StringValueForKey("symmetry planes",LINE_LENGTH)
         sym_planes = getIntArrayFromString(tmp)
         if( maxval(abs(sym_planes)) .gt. 3 ) then
            error stop " IBM_GetInfo ::  Symmetry planes allowed values -3:3"
         end if 
         allocate( this% symPlanes(size(sym_Planes,1)), this% symCoords(size(sym_Planes,1)) )
         do i = 1, size(sym_Planes,1)
            this% symPlanes(i) = sym_planes(i) 
            if( sym_planes(i) .gt. 0 ) then
                this% symCoords(i) = -huge(1.0_RP)
             else
                this% symCoords(i) = huge(1.0_RP)
             end if
         end do 
       end if     

      
   end subroutine IBM_GetInfo
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine describes the IBM info
!  -----------------------------------------------
   subroutine IBM_Describe( this )
      use Headers
      use mainKeywordsModule
      use MPI_Process_Info
      implicit none
      !-arguments--------------------------------
      class(IBM_type),  intent(inout) :: this
      
      if ( MPI_Process % isRoot ) then
      write(STD_OUT,'(/)')
      call Section_Header("IBM parameters")
      write(STD_OUT,'(/)')
      
      call SubSection_Header('IBM info')

      write(STD_OUT,'(30X,A,A35,L10)') "->" , "Semi implicit treatment: ", this% active_semiImplicit
      if( .not. this% TimePenal ) then
         write(STD_OUT,'(30X,A,A35,1pG10.6)') "->" , "Penalization term: " , this% eta
      else
         write(STD_OUT,'(30X,A,A35,A10)') "->" , "Penalization term: ", " Dt"
      end if

      write(STD_OUT,'(30X,A,A35,I10)') "->" , "Minimum number of objects: ", this% KDtree_Min_n_of_Objs
      write(STD_OUT,'(30X,A,A35,I10)') "->" , "Number of interpolation points: ", this% KDtree_n_of_interPoints
      write(STD_OUT,'(30X,A,A35,I10)') "->" , "Surface integration order: ", this% IntegrationOrder
      if( this% Wallfunction ) then
         write(STD_OUT,'(30X,A,A35,F10.3)') "->" , "Target y+: ", this% y_plus_target
      end if
      if( this% symmetry ) then
         select case( size(this% symPlanes,1) )
         case( 1 )
            write(STD_OUT,'(30X,A,A35,I2,A)') "->" , "Symmetry planes: [", this% symPlanes(1),"]"
         case( 2 )
            write(STD_OUT,'(30X,A,A36,I2,A,I2,A)') "->" , "Symmetry planes: [", this% symPlanes(1),",",this% symPlanes(2),"]"
         case( 3 )
            write(STD_OUT,'(30X,A,A35,I2,A,I2,A,I2,A)') "->" , "Symmetry planes: [", this% symPlanes(1),",",this% symPlanes(2),",",this% symPlanes(3),"]"
         case default
         write(STD_OUT,'(30X,A,A35,I10)') "->" , "Symmetry planes: ", this% symPlanes
         end select
      end if
      end if
   end subroutine IBM_Describe
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine computes the IBM source term 
!  ------------------------------------------------   

   subroutine IBM_SourceTerm( this, eID, Q, Q_target, Source )
      use PhysicsStorage
      implicit none
      !-arguments------------------------------------------------
      class(IBM_type),                 intent(inout) :: this
      integer,                         intent(in)    :: eID
      real(kind=rp), dimension(NCONS), intent(in)    :: Q, Q_target
      real(kind=rp), dimension(NCONS), intent(inout) :: Source
      !-local-variables-----------------------------------
      real(kind=rp) :: rho, rho_s, u, u_s, v, v_s, w, w_s
#if defined(SPALARTALMARAS)
      real(kind=rp) :: theta, theta_s
#endif

      optional :: Q_target
      
      Source = 0.0_RP
      
#if defined(NAVIERSTOKES)
      if( present(Q_target) ) then
         rho_s = Q_target(IRHO)
         u_s   = Q_target(IRHOU)/rho_s
         v_s   = Q_target(IRHOV)/rho_s
         w_s   = Q_target(IRHOW)/rho_s   
#if defined(SPALARTALMARAS)
         theta_s = Q_target(IRHOTHETA)/rho_s  
#endif
      else
         rho_s = Q(IRHO)
         u_s   = 0.0_RP
         v_s   = 0.0_RP
         w_s   = 0.0_RP
#if defined(SPALARTALMARAS)
         theta_s = 0.0_RP  
#endif
      end if   
        
      rho = Q(IRHO)
      u   = Q(IRHOU)/rho
      v   = Q(IRHOV)/rho
      w   = Q(IRHOW)/rho
      
      Source(IRHOU) = rho*(u-u_s)  
      Source(IRHOV) = rho*(v-v_s) 
      Source(IRHOW) = rho*(w-w_s)
      Source(IRHOE) = 0.5_RP * rho*( POW2(u) + POW2(v) + POW2(w) )  & 
                     -0.5_RP * rho_s*( POW2(u_s) + POW2(v_s) + POW2(w_s) ) 
#if defined(SPALARTALMARAS)
      theta = Q(IRHOTHETA)/rho

      Source(IRHOTHETA) = rho*(theta - theta_s)
#endif                     
#endif   

      Source = -1.0_RP/this% penalization(eID) * Source
   
   end subroutine IBM_SourceTerm
   
   subroutine IBM_semiImplicitShiftJacobian( this, eID, Q, dt, invdS_dQ )
      use PhysicsStorage
      implicit none
         !-arguments----------------------------------------------------------------
         class(IBM_type),                       intent(inout) :: this
         integer,                               intent(in)    :: eID
         real(kind=rp), dimension(NCONS),       intent(in)    :: Q
         real(kind=rp),                         intent(in)    :: dt
         real(kind=rp), dimension(NCONS,NCONS), intent(inout) :: invdS_dQ
         !-local-variables----------------------------------------------------------
         real(kind=rp) :: rho, u, v, w

         invdS_dQ = 0.0_RP  
         
         associate( eta => this% penalization(eID) )     

#if defined(NAVIERSTOKES) 
         rho = Q(IRHO)
         u   = Q(IRHOU)/rho 
         v   = Q(IRHOV)/rho 
         w   = Q(IRHOW)/rho 
         
         invdS_dQ(IRHO,IRHO)        = 1.0_RP
         invdS_dQ(IRHOU,IRHOU)      = eta/( dt + eta )
         invdS_dQ(IRHOV,IRHOV)      = eta/( dt + eta )
         invdS_dQ(IRHOW,IRHOW)      = eta/( dt + eta )
         invdS_dQ(IRHOE,IRHO:IRHOE) = (/ 0.5_RP*dt/eta * (POW2(u) + POW2(v) + POW2(w)), &
                                                                      -dt*u/(dt + eta), &
                                                                      -dt*v/(dt + eta), &
                                                                      -dt*w/(dt + eta), &
                                                                                1.0_RP /)
#if defined(SPALARTALMARAS)
         invdS_dQ(IRHOTHETA,IRHOTHETA) = eta/( dt + eta )
#endif                                                                           
#endif
 
       end associate
 
   end subroutine IBM_SemiImplicitShiftJacobian
    
   subroutine IBM_semiImplicitJacobian( this, eID, Q, dS_dQ )
      use PhysicsStorage
      implicit none
      !-arguments----------------------------------------------------
      class(IBM_type),                       intent(inout) :: this
      integer,                               intent(in)    :: eID
      real(kind=rp), dimension(NCONS),       intent(in)    :: Q
      real(kind=rp), dimension(NCONS,NCONS), intent(inout) :: dS_dQ
      !-local-variables----------------------------------------------
      real(kind=rp) :: rho, u, v, w

      dS_dQ = 0.0_RP       

#if defined(NAVIERSTOKES) 
      rho = Q(IRHO)
      u   = Q(IRHOU)/rho 
      v   = Q(IRHOV)/rho 
      w   = Q(IRHOW)/rho 

      dS_dQ(IRHOU,IRHOU)      = 1.0_RP
      dS_dQ(IRHOV,IRHOV)      = 1.0_RP
      dS_dQ(IRHOW,IRHOW)      = 1.0_RP
      dS_dQ(IRHOE,IRHO:IRHOE) = (/ -0.5_RP*( POW2(u) + POW2(v) + POW2(w) ), &
                                                           u, v, w, 0.0_RP /)
#if defined(SPALARTALMARAS)
       dS_dQ(IRHOTHETA,IRHOTHETA) = 1.0_RP
#endif
#endif
       
       dS_dQ = -1.0_RP/this% penalization(eID) * dS_dQ
 
    end subroutine IBM_SemiImplicitJacobian
    
    subroutine IBM_GetSemiImplicitStep( this, eID, dt, Q )
       use PhysicsStorage
       implicit none
       !-arguments----------------------------------------------
       class(IBM_type),                 intent(inout) :: this
       integer,                         intent(in)    :: eID
       real(kind=rp),                   intent(in)    :: dt
       real(kind=rp), dimension(NCONS), intent(inout) :: Q
       !-local-variables----------------------------------------
       real(kind=rp), dimension(NCONS,NCONS) :: dS_dQ, invdS_dQ
       real(kind=rp), dimension(NCONS)       :: IBMSource
       
       call this% semiImplicitJacobian( eID, Q, dS_dQ )
       call this% semiImplicitShiftJacobian( eID, Q, dt, invdS_dQ )
    
       call this% SourceTerm(eID = eID, Q = Q, Source = IBMSource)          
    
       Q = matmul(invdS_dQ, Q + dt*( IBMSource - matmul(dS_dQ,Q) ))
     
    end subroutine IBM_GetSemiImplicitStep   
    
    
    subroutine IBM_SemiImplicitCorrection( this, elements, dt, dt_vec )
    
       implicit none
       !-arguments------------------------------------------------------
       class(IBM_type),               intent(inout) :: this
       type(element),   dimension(:), intent(inout) :: elements
       real(kind=RP),                 intent(in)    :: dt, dt_vec(:)
       !-local-variables------------------------------------------------
       real(kind=RP) :: deltaT
       integer       :: eID, i, j, k
       
       optional :: dt_vec
       
       if( .not. this% semiImplicit ) return
       
!$omp parallel do schedule(runtime) private(i,j,k,deltaT)
       do eID = 1, SIZE( elements )
          if( present(dt_vec) ) then 
             deltaT = dt_vec(eID)
          else
             deltaT = dt
          end if
          
          if( this% TimePenal ) this% penalization(eID) = 0.5_RP*deltaT
          do i = 0, elements(eID)% Nxyz(1); do j = 0, elements(eID)% Nxyz(2); do k = 0, elements(eID)% Nxyz(3)
             if( elements(eID)% isInsideBody(i,j,k) ) then
                associate( Q => elements(eID)% storage% Q(:,i,j,k) )

                call this% GetSemiImplicitStep( eID, 0.5_RP*deltaT, Q ) 
                
                end associate
             end if
          end do; end do; end do
       end do
!$omp end parallel do
    
    end subroutine IBM_SemiImplicitCorrection
    
    
   
   subroutine IBM_SetIntegration( this, STLNum )
      use MPI_Process_Info
      use MappedGeometryClass
      implicit none
      !-arguments---------------------------------------------
      class(IBM_type), intent(inout) :: this
      integer,         intent(in)    :: STLNum
      !-local-variables---------------------------------------
      real(kind=rp)              :: coord, Vertices(NDIM,3), & 
                                    T(NDIM)
      integer                    :: i, j, k, axis,           &
                                    n_of_Q_points,           &
                                    kdtree_n_of_interPoints, &
                                    symPlaneIndex
      
      if( this% symmetry ) then
!$omp parallel shared(this,STLNum,i)
!$omp do schedule(runtime) private(k,j,coord,axis,Vertices)  
         do i = 1, size(this% root(STLNum)% ObjectsList)
            
            this% root(STLNum)% ObjectsList(i)% ComputeIntegrals = .true.
            
            call OBB(STLNum)% ChangeRefFrame( this% root(STLNum)% ObjectsList(i)% vertices(1)% coords, 'global', Vertices(:,1) )
            call OBB(STLNum)% ChangeRefFrame( this% root(STLNum)% ObjectsList(i)% vertices(2)% coords, 'global', Vertices(:,2) )
            call OBB(STLNum)% ChangeRefFrame( this% root(STLNum)% ObjectsList(i)% vertices(3)% coords, 'global', Vertices(:,3) )

            do k = 1, size(this% symPlanes,1)      
               axis = abs(this% symPlanes(k))            
               if( this% symPlanes(k) .lt. 0 .and. this% root(STLNum)% ObjectsList(i)% ComputeIntegrals ) then
                  coord = minval(this% rootPoints% vertices(axis,:))
                  if( Vertices(axis,1) .lt. coord .and. &
                      Vertices(axis,2) .lt. coord .and. &
                      Vertices(axis,3) .lt. coord        ) then
                      this% root(STLNum)% ObjectsList(i)% ComputeIntegrals = .false.
                  else
                     do j = 1, size(this% root(STLNum)% ObjectsList(i)% vertices)
                        if( Vertices(axis,j) .lt. coord ) this% root(STLNum)% ObjectsList(i)% vertices(j)% Translate = k
                     end do
                  end if
               elseif( this% root(STLNum)% ObjectsList(i)% ComputeIntegrals ) then
                  coord = maxval(this% rootPoints% vertices(axis,:))
                  if( Vertices(axis,1) .gt. coord .and. &
                      Vertices(axis,2) .gt. coord .and. &
                      Vertices(axis,3) .gt. coord       ) then
                      this% root(STLNum)% ObjectsList(i)% ComputeIntegrals = .false.
                  else
                     do j = 1, size(this% root(STLNum)% ObjectsList(i)% vertices)
                        if( Vertices(axis,j) .gt. coord ) this% root(STLNum)% ObjectsList(i)% vertices(j)% Translate = k
                     end do       
                  end if
               end if   
            end do            
         end do  
!$omp end do
!$omp end parallel 

         do k = 1, size(this% symPlanes,1)
            axis = abs(this% symPlanes(k))
            if( this% symPlanes(k) .lt. 0 ) then
               coord = minval(this% rootPoints% vertices(axis,:))
            else
               coord = maxval(this% rootPoints% vertices(axis,:))
            end if
            this% symCoords(k) = coord
         end do
      end if      
      
      n_of_Q_points           = this% Integral(STLNum)% n_of_Q_points
      kdtree_n_of_interPoints = this% kdtree_n_of_interPoints

!$omp parallel shared(this,n_of_Q_points,kdtree_n_of_interPoints,i,STLNum)
!$omp do schedule(runtime) private(Vertices,symPlaneIndex,T,j)
      do i = 1, size(this% root(STLNum)% ObjectsList)
         if( this% root(STLNum)% ObjectsList(i)% ComputeIntegrals ) then  
          !if body is moving, we deallocate and reallocate
            if( allocated(this% Integral(STLNum)% IntegObjs(i)% PointsIndex) )      deallocate(this% Integral(STLNum)% IntegObjs(i)% PointsIndex)
            if( allocated(this% Integral(STLNum)% IntegObjs(i)% x          ) )      deallocate(this% Integral(STLNum)% IntegObjs(i)% x)
                            
            allocate( this% Integral(STLNum)% IntegObjs(i)% PointsIndex(kdtree_n_of_interPoints,n_of_Q_points), &
                      this% Integral(STLNum)% IntegObjs(i)% x(NDIM,n_of_Q_points)                               )

            do j = 1, size(this% root(STLNum)% ObjectsList(i)% vertices) 
               Vertices(:,j) = this% root(STLNum)% ObjectsList(i)% vertices(j)% coords
               symPlaneIndex = this% root(STLNum)% ObjectsList(i)% vertices(j)% Translate
               if( symPlaneIndex .gt. 0 .and. this% symmetry ) then
                  call OBB(STLNum)% ChangeRefFrame( Vertices(:,j), 'global', Vertices(:,j) )
                  Vertices(abs(this% symPlanes(symPlaneIndex)),j) = this% symCoords( symPlaneIndex )
                  call OBB(STLNum)% ChangeRefFrame( Vertices(:,j), 'local', Vertices(:,j) )
               end if
            end do
            
            call this% Integral(STLNum)% GetCoords( Vertices(:,1), Vertices(:,2), Vertices(:,3), i, n_of_Q_points )

            call vcross( (Vertices(:,3)-Vertices(:,1)), (Vertices(:,2)-Vertices(:,1)), T )
            
            this% Integral(STLNum)% IntegObjs(i)% Area = norm2(T)
      
         end if
      end do
!$omp end do
!$omp end parallel
      
   end subroutine IBM_SetIntegration
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------------------
! This subroutine returns true if the point is inside the body. 
!  -------------------------------------------------------------
   subroutine IBM_CheckPoint( this, x, STLNum, isInsideBody )
      use MPI_Process_Info
      
      use omp_lib
      
      implicit none
      !-arguments------------------------------------------------------------------
      class(IBM_type),             intent(inout) :: this
      real(kind=rp), dimension(:), intent(in)    :: x
      integer,                     intent(in)    :: STLNum
      logical,                     intent(inout) :: isInsideBody
      !-local-variables------------------------------------------------------------
      type(KDtree),  pointer         :: tree 
      real(kind=rp), dimension(NDIM) :: RayDirection, Point, vecAxis
      integer                        :: Axis, NumOfIntersections, minAxis(1)
      logical                        :: Upward, OnSurface
      type(ObjsDataLinkedList_t)     :: Integer_List  
 
      real(kind=RP) :: EPS = 1.0d-8

      call OBB(STLNum)% ChangeRefFrame(x, 'local', Point)
      
      RayDirection = 0.0_RP
            
      isInsideBody = .false.
      OnSurface    = .false.

      vecAxis = (/OBB(STLNum)% MBR% Length,OBB(STLNum)% MBR% Width,abs(OBB(STLNum)% nMax) + abs(OBB(STLNum)% nMin)/)
      
      minAxis = minloc(vecAxis)
      axis = minAxis(1)
      
      if( axis .eq. this% root(STLNum)% maxAxis ) then
         axis = axis - 1
         if( axis .eq. 0 ) axis = 3
      end if
      
      if( Point(axis) .le. 0.0_RP ) then
         RayDirection(axis) = -1.0_RP
         Upward = .true.
         vecAxis =  (/ -0.5_RP*OBB(STLNum)% MBR% Length, -0.5_RP*OBB(STLNum)% MBR% Width, &
                        OBB(STLNum)% nMin                                                /)
      else
         RayDirection(axis) = 1.0_RP
         Upward = .false.
         vecAxis =  (/ 0.5_RP*OBB(STLNum)% MBR% Length, 0.5_RP*OBB(STLNum)% MBR% Width, &
                       OBB(STLNum)% nMax                                               /)
      end if

      Integer_List = ObjsDataLinkedList_t( )

      NumOfIntersections = 0

      do 

         call this% root(STLNum)% FindLeaf( Point, tree )
         
         call isPointInside( Point, RayDirection, this% root(STLNum)% ObjectsList, tree, Integer_List, NumOfIntersections, OnSurface )

         if( OnSurface ) exit
         
         if( Upward ) then
            Point(axis) = tree% vertices(axis,1) - EPS
            if( Point(axis) .lt. vecAxis(axis) ) exit
         elseif( .not. Upward ) then
            Point(axis) = tree% vertices(axis,7) + EPS
            if( Point(axis) .gt. vecAxis(axis) ) exit
         end if

      end do

      if( mod(NumOfIntersections, 2) .eq. 0 ) then !even
         isInsideBody = .false.
      else
         isInsideBody = .true.
      end if 
      
      if( OnSurface ) isInsideBody = .true.

      call integer_List% Destruct()

   end subroutine IBM_CheckPoint
   
   
   subroutine IBM_ComputeIBMWallDistance( this, elements, faces, elementsList, facesList)

      implicit none
      !-arguments---------------------------------------------------------------
      class(IBM_type),    intent(inout) :: this
      type(element),      intent(inout) :: elements(:)
      type(face),         intent(inout) :: faces(:)
      integer, optional , intent(in)    :: facesList(:)
      integer, optional , intent(in)    :: elementsList(:)
      !-local-variables----------------------------------------------------------
      integer       :: eID, fID, ii, i, j, k, num_of_elems, num_of_faces, STLNum
      real(kind=RP) :: xP(NDIM), IntersectionPoint(NDIM), Dist, Point(NDIM), normal(NDIM)

      num_of_elems = size(elements)   
      num_of_faces = size(faces) 

!$omp parallel
!$omp do schedule(runtime) private(eID,i,j,k,STLNum,Dist,Point,xP,normal)
      do ii = 1, num_of_elems
         if ( present(elementsList) ) then
            eID = elementsList(ii)
         else
            eID = ii
         end if
         allocate(elements(eID)% geom% dWall(0:elements(eID)% Nxyz(1), 0:elements(eID)% Nxyz(2), 0:elements(eID)% Nxyz(3)),      &
                  elements(eID)% geom% normal(NDIM,0:elements(eID)% Nxyz(1), 0:elements(eID)% Nxyz(2), 0:elements(eID)% Nxyz(3)) )

         do k = 0, elements(eID)% Nxyz(3)   ; do j = 0, elements(eID)% Nxyz(2) ; do i = 0, elements(eID)% Nxyz(1)
            xP = elements(eID)% geom% x(:,i,j,k)
            elements(eID)% geom% dWall(i,j,k) = huge(1.0_RP)
            do STLNum = 1, this% NumOfSTL
               if( OBB(STLNum)% isPointInside( xP, this% BandRegionCoeff ) ) then
                  call MinimumDistance( xP, this% root(STLNum), Dist, normal )
                  if( Dist .lt. elements(eID)% geom% dWall(i,j,k) ) then
                     elements(eID)% geom% dWall(i,j,k)    = Dist
                     elements(eID)% geom% normal(:,i,j,k) = normal
                  end if
               end if
            end do            
         end do                  ; end do                ; end do
      end do
!$omp end do
!
!     Get the minimum distance to each face nodal degree of freedom
!     -------------------------------------------------------------

!$omp do schedule(runtime) private(eID,i,j,STLNum,Dist,Point,xP,normal)
      do ii = 1, num_of_faces
         if ( present(facesList) ) then
            eID = facesList(ii)
         else
            eID = ii
         end if

         allocate(faces(eID)% geom% dWall(0:faces(eID)% Nf(1), 0:faces(eID)% Nf(2)))

         do j = 0, faces(eID)% Nf(2) ; do i = 0, faces(eID)% Nf(1)
            xP = faces(eID)% geom% x(:,i,j)
            faces(eID)% geom% dWall(i,j) = huge(1.0_RP)
            do STLNum = 1, this% NumOfSTL
               if( OBB(STLNum)% isPointInside( xP, this% BandRegionCoeff ) ) then
                  call MinimumDistance( xP, this% root(STLNum), Dist, normal )
                  if( Dist .lt. faces(eID)% geom% dWall(i,j) ) then
                     faces(eID)% geom% dWall(i,j)    = Dist
                     faces(eID)% geom% normal(:,i,j) = normal
                  end if
               end if
            end do
         end do                ; end do
      end do
!$omp end do
!$omp end parallel

      if( this% Wallfunction ) then
#if defined(NAVIERSTOKES)
         call this% GetForcingPointsGeom( elements )
         call this% GetImagePoint_nearest( elements )
#endif
      end if    

   end subroutine IBM_ComputeIBMWallDistance
   
   
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!                                  BODY MOTION
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------  
      
   subroutine IBM_MoveBody( this, elements, no_of_elements, no_of_DoFs, isChild, dt, k, autosave )
      use MPI_Process_Info
      implicit none
      !-arguments------------------------------
      class(IBM_type),               intent(inout) :: this
      type(element),   dimension(:), intent(inout) :: elements
      integer,                       intent(in)    :: no_of_elements, no_of_DoFs, k
      logical,                       intent(in)    :: isChild 
      real(kind=RP),                 intent(in)    :: dt
      logical,                       intent(in)    :: autosave
      !-local-variables------------------------
      integer :: STLNum, MaskPoints
#ifdef _HAS_MPI_
      integer :: localVal, ierr
#endif
      
      optional :: k, autosave         

      do STLNum = 1, this% NumOfSTL
         if( this% stl(STLNum)% move .and. .not. isChild ) then
            if( this% stl(STLNum)% motionType .eq. ROTATION ) then
               call this% stl(STLNum)% getRotationaMatrix( dt )
               call OBB(STLNum)% STL_rotate(this% stl(STLNum))
            elseif( this% stl(STLNum)% motionType .eq. LINEAR ) then
               call this% stl(STLNum)% getDisplacement( dt )
               call OBB(STLNum)% STL_translate(this% stl(STLNum))
            end if
            if( .not. isChild .and. present(autosave) .and. autosave ) call this% stl(STLNum)% plot( k )
            call OBB(STLNum)% construct( this% stl(STLNum), this% plotOBB )
            this% root(STLNum)% STLNum = STLNum
            call OBB(STLNum)% ChangeObjsRefFrame( this% stl(STLNum)% ObjectsList )
            call this% root(STLNum)% Destruct( isChild )
            this% plotKDtree = .false.
            call this% constructSTL_KDtree( STLNum ) 
            call this% upDateNormals( STLNum )
            call this% CleanMask( elements, no_of_elements, STLNum )
         elseif( this% stl(STLNum)% move .and. isChild ) then
            call this% CleanMask( elements, no_of_elements, STLNum )
         end if
      end do

      call this% BandRegion% destroy()

      if( .not. isChild ) then
         call this% rootPoints% destruct( isChild )
      end if

      do STLNum = 1, this% NumOfSTL
         if( this% stl(STLNum)% move ) then      
            ! get new mask
            call this% GetMask( elements, no_of_elements, no_of_DoFs, STLNum )
         end if
      end do

      if( MPI_Process% doMPIAction ) then
#ifdef _HAS_MPI_
         localVal = this% n_of_INpoints
         call mpi_allreduce(localVal, MaskPoints, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
      else
         MaskPoints = this% n_of_INpoints
      end if    

      if( MaskPoints .eq. 0 ) then
         print *, "The mask is made of 0 points."
         print *, "Try to increase the polynomial order or refine the mesh."
         error stop
      end if

      if( .not. isChild ) then
         call this% constructBandRegion( elements, no_of_elements )
         do STLNum = 1, this% NumOfSTL
            call this% SetIntegration( STLNum )
         end do
         if( this% Wallfunction ) then
#if defined(NAVIERSTOKES)
            call this% GetForcingPointsGeom( elements )
            call this% GetImagePoint_nearest( elements )
#endif
         end if
      end if

      if( this% plotMask .and. .not. isChild ) call this% plot_Mask( elements, no_of_elements, k )

   end subroutine IBM_MoveBody
   
   
   subroutine IBM_CleanMask( this, elements, no_of_elements, STLNum )
   
      implicit none
   
      class(IBM_type),             intent(inout) :: this
      type(element), dimension(:), intent(inout) :: elements
      integer,                     intent(in)    :: no_of_elements, STLNum
   
      integer :: eID, i, j, k
!$omp parallel shared(this,elements,no_of_elements,STLNum,eID)
!$omp do schedule(runtime) private(i,j,k)
      do eID = 1, no_of_elements
         do i = 0, elements(eID)% Nxyz(1); do j = 0, elements(eID)% Nxyz(2); do k = 0, elements(eID)% Nxyz(3) 
            if( any(elements(eID)% STL(:,i,j,k) .eq. STLNum) .and. elements(eID)% isInsideBody(i,j,k) ) then
               elements(eID)% STL(STLNum,i,j,k) = 0
               if(all(elements(eID)% STL(:,i,j,k) .eq. 0)) then
                  elements(eID)% isInsideBody(i,j,k) = .false.
!$omp critical
                  this% n_of_INpoints = this% n_of_INpoints - 1
!$omp end critical
               end if
            end if
         end do; end do; end do
      end do 
!$omp end do
!$omp end parallel
   end subroutine IBM_CleanMask
   
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine checks if a point intersects one or 
! more of the triangles belonging to a box (tree).
!  ------------------------------------------------
   subroutine isPointInside( Point, RayDirection, ObjectsList, tree, Integer_List, NumOfIntersections, OnSurface )
      use RealDataLinkedList
      use omp_lib
      implicit none
      !-arguments--------------------------------------------------------------
      real(kind=rp), dimension(:),     intent(in)    :: Point, RayDirection
      type(object_type), dimension(:), intent(in)    :: ObjectsList
      type(KDtree),                    intent(inout) :: tree
      type(ObjsDataLinkedList_t),      intent(inout) :: Integer_List
      integer,                         intent(inout) :: NumOfIntersections
      logical,                         intent(inout) :: OnSurface
      !-local-variables--------------------------------------------------------
      logical                       :: Intersect, OnTriBound, found, foundReal
      real(kind=rp)                  :: t
      integer                        :: i
      type(ObjsRealDataLinkedList_t) :: Real_List

      integer :: index

      if( tree% NumOfObjs .eq. 0 ) then
         return
      end if 
      
      Real_List = ConstructObjsRealDataLinkedList()
      
      do i = 1, tree% NumOfObjs
         index = tree% ObjsIndeces(i)
         found = integer_List% Check( index )
         if( .not. found ) then
            call Integer_List% Add( index )
            
            call PointIntersectTriangle( Point,ObjectsList(index)% vertices(1)% coords, &
                                               ObjectsList(index)% vertices(2)% coords, &
                                               ObjectsList(index)% vertices(3)% coords, &
                                         RayDirection, Intersect, OnTriBound, t         ) 
            
            foundReal = Real_List% Check( t )
            foundReal = .false.
            
            if( Intersect .and. .not. foundReal .and. t .ge. 0.0_RP ) then
               call Real_List% add( t )
               NumOfIntersections = NumOfIntersections + 1 
!~                ! If the point is on the triangle
!~                !--------------------------------
               if( almostEqual(t,0.0_RP) ) OnSurface = .true.
            end if
         end if
      end do

      call Real_List% destruct()

   end subroutine isPointInside  
   
   
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  ---------------------------------------------------------------------
! This subroutine checks if a ray (RayDirection) starting from a point (Point) intersects
! a traingle in 3D space. If present, the intersection point Q is Q = P + t*RayDirection.
! If there is more than one intersections, the second one is not counted if it has the same t
! as the one previously found.
! See Fast, Minimum Storage Ray/Trinagle INtersection,  Moller TRumbore
!  ---------------------------------------------------------------------
   
   subroutine PointIntersectTriangle( Point, TriangleVertex1, TriangleVertex2, &
                                      TriangleVertex3, RayDirection,           &
                                      Intersect, OnTriBound, t                 )
      use MappedGeometryClass
      implicit none
      !-arguments---------------------------------------------------------------------------------------
      real(kind=rp), dimension(NDIM), intent(in)  :: Point, RayDirection
      real(kind=rp), dimension(NDIM), intent(in)  :: TriangleVertex1, TriangleVertex2, TriangleVertex3
      logical,                        intent(out) :: Intersect, OnTriBound 
      real(kind=rp),                  intent(out) :: t
      !-local-variables----------------------------------------------------------------------------------
      real(kind=rp), dimension(NDIM) :: E1vec, E2vec, Pvec, Qvec, Tvec, N
      real(kind=rp)                  :: Det, u, v, invDet
      logical                        :: isInside
      
      Intersect  = .false.
      OnTriBound = .false.
      t          = -1.0_RP
      
      E1vec = TriangleVertex2 - TriangleVertex1 
      E2vec = TriangleVertex3 - TriangleVertex1
      Tvec   = Point - TriangleVertex1
      
      call vcross(RayDirection,E2vec,Pvec)
      Det   = dot_product( E1vec, Pvec )

      If( almostEqual(Det,0.0_RP) ) then
      ! If Pvec .ne. (/0,0,0/), then the vector E1vec must lie on the plane made of the 2 vectors RayDirection and E2vec.
      ! Thus we have to check if the ray intersects the triangle or not. In the latter case no intersection is detected. In
      ! the first case, we have to check whether the point is inside the triangle or not. If it is in the triangle, one intersection  
      ! will be detected. 
      !------------------------------------------------------------------------------------------------------------------------------
         call vcross(E1vec, E2vec, N)
         Intersect = .false.
        
         t = -1.0_RP
         ! Check if the ray lies in the same plane of the triangle
         !--------------------------------------------------------
         if( AlmostEqual(dot_product( Tvec, N ), 0.0_RP) ) then 
            isInside = isPointInsideTri( Point, TriangleVertex1, TriangleVertex2, &
                                         TriangleVertex3                          )
                                         
            if( isInside ) t = 0.0_RP
            if( isInside ) Intersect = .true.
         end if
         return
      end if
      
      call vcross(Tvec,E1vec,Qvec)
      
      invDet = 1.0_RP/Det
      
      u = dot_product( Tvec, Pvec )*invDet

      if( u < 0.0_RP .or. u > 1.0_RP ) return
      
      v = dot_product( RayDirection, Qvec )/Det
      

      if( v < 0.0_RP .or. u+v > 1.0_RP ) return

      t  = dot_product( E2vec, Qvec )*invDet
   
      ! Check if the point lies on the boundaries of the triangle
      !----------------------------------------------------------
      if( almostEqual(u,0.0_RP) .or. almostEqual(v,0.0_RP) .or. &
          almostEqual(u+v,1.0_RP) )  OnTriBound = .true.

      Intersect = .true.

   end subroutine PointIntersectTriangle
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This function checks if a point is inside a triangle. 
! The point and the triangle are on the same plane 
!  ------------------------------------------------   
   logical function isPointInsideTri( Point, TriangleVertex1, TriangleVertex2, &
                                      TriangleVertex3 ) result( isInside )
      use MappedGeometryClass
      implicit none
      !-arguments-----------------------------------------------------------------
      real(kind=rp), dimension(:), intent(in) :: Point, TriangleVertex1, &
                                                 TriangleVertex2, TriangleVertex3
      !-local-variables-----------------------------------------------------------
      real(kind=rp), dimension(NDIM) :: bb, E0, E1, dd
      real(kind=rp) :: a, b, c, d, e, f, det, s, t
      integer :: region
      
      isInside = .false.
      
      bb = TriangleVertex1
      E0 = TriangleVertex1 - TriangleVertex2
      E1 = TriangleVertex1 - TriangleVertex3
      dd = bb - Point
   
      a = dot_product(E0,E0)
      b = dot_product(E0,E1)
      c = dot_product(E1,E1)
      d = dot_product(E0,dd)
      e = dot_product(E1,dd)
      f = dot_product(dd,dd)
      
      det = abs( a*c - b*b)
      s    = (b*e - c*d)/det
      t    = (b*d - a*e)/det
      
      region = FindRegion( det, s, t )
      
      if( region .eq. 0 ) isInside = .true.
      
   end function isPointInsideTri
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This function computes the minimum distance from a point to a triangle in 3D. 
! for more ditails see https://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
!  ------------------------------------------------   
   subroutine MinumumPointTriDistance( Point, TriangleVertex1, TriangleVertex2, &
                                       TriangleVertex3, dist, IntersectionPoint ) 
      use MappedGeometryClass
      implicit none
      !-arguments-----------------------------------------------------------------------------------------
      real(kind=rp), dimension(:),    intent(in)  :: Point, TriangleVertex1, &
                                                     TriangleVertex2,        &
                                                     TriangleVertex3
      real(kind=rp), dimension(NDIM), intent(out) :: IntersectionPoint 
      real(kind=rp),                  intent(out) :: dist
      !-local-variables-----------------------------------------------------------------------------------
      real(kind=rp), dimension(NDIM) :: bb, E0, E1, dd 
      real(kind=rp) :: a00, a01, a11, b0, b1, f, det, s, t, tmp1, tmp0, numer, denom
      integer :: region
      
      bb = TriangleVertex1
      E0 = TriangleVertex2 - bb
      E1 = TriangleVertex3 - bb
      dd = bb - Point
   
      a00 = dot_product(E0,E0) 
      a01 = dot_product(E0,E1) 
      a11 = dot_product(E1,E1) 
      b0  = dot_product(dd,E0) 
      b1  = dot_product(dd,E1) 
      f   = dot_product(dd,dd)
      
      det = max(a00*a11 - a01*a01,0.0_RP)  
      s   = a01*b1 - a11*b0      
      t   = a01*b0 - a00*b1
      
      
      if( s + t <= det ) then
         if( s < 0.0_RP ) then
            if( t < 0.0_RP ) then !region 4
               if( b0 < 0.0_RP ) then
                  t = 0.0_RP
                  if( -b0 >= a00 ) then
                     s = 1.0_RP
                  else
                     s = -b0/a00
                  end if
               else
                  s = 0.0_RP
                  if( b1 >= 0.0_RP ) then
                     t = 0.0_RP
                  elseif( -b1 >= a11 ) then
                     t = 1.0_RP
                  else
                     t = -b1/a11
                  end if
               end if
            else !region 3
               s = 0.0_RP
               if( b1 >= 0.0_RP ) then
                  t = 0.0_RP
               elseif( -b1 >= a11 ) then
                  t = 1.0_RP
               else
                  t = -b1/a11
               end if
            end if
         elseif( t < 0.0_RP ) then !region 5
            t = 0.0_RP
            if( b0 >= 0.0_RP ) then
               s = 0.0_RP
            elseif( -b0 >= a00 ) then
               s = 1.0_RP
            else
               s = -b0/a00
            end if
         else !region 0
            s = s/det
            t = t/det
         end if
      else
         if( s < 0.0_RP ) then !region 2
            tmp0 = a01 + b0
            tmp1 = a11 + b1
            if (tmp1 > tmp0) then
               numer = tmp1 - tmp0
               denom = a00 - 2.0_RP * a01 + a11
               if (numer >= denom) then
                  s = 1.0_RP
                  t = 0.0_RP
               else
                  s = numer / denom
                  t = 1.0_RP - s
               end if
            else
               s = 0.0_RP
               if ( tmp1 <= 0.0_RP ) then
                  t = 1.0_RP
               elseif( b1 >= 0.0_RP ) then
                  t = 0.0_RP
               else
                  t = -b1 / a11
               end if
            end if
         elseif( t < 0.0_RP ) then !region 6
            tmp0 = a01 + b1
            tmp1 = a00 + b0
            if( tmp1 > tmp0 ) then
               numer = tmp1 - tmp0
               denom = a00 - 2.0_RP * a01 + a11
               if (numer >= denom ) then
                  t = 1.0_RP
                  s = 0.0_RP
               else
                  t = numer / denom;
                  s = 1.0_RP - t;
               endif
            else
               t = 0.0_RP
               if( tmp1 <= 0.0_RP ) then
                  s = 1.0_RP
               elseif( b0 >= 0.0_RP ) then
                  s = 0.0_RP                 
               else
                  s = -b0 / a00
               end if
            end if
         else  ! region 1
            numer = a11 + b1 - a01 - b0
            if( numer <= 0.0_RP ) then
               s = 0.0_RP
               t = 1.0_RP
            else
               denom = a00 - 2.0_RP * a01 + a11
               if( numer >= denom ) then
                  s = 1.0_RP
                  t = 0.0_RP
               else
                  s = numer / denom
                  t = 1.0_RP - s
               end if
            end if
         end if
      end if

      IntersectionPoint = TriangleVertex1 + s*E0 + t*E1    
        
      dist = norm2(Point - IntersectionPoint)

   end subroutine MinumumPointTriDistance
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  --------------------------------------------------------------------------
! This function detects the region, see below.
!  --------------------------------------------------------------------------
!      \ region2 ^t
!       \        |
!        \       |
!         \      | 
!          \     |
!           \    |
!            \   |
!             \  |
!              \ |
!               \|
!                *
!                |\    P2
!                | \
!                |  \
!                |   \
!                |    \
!                |     \
!                |      \
!                |       \
!                |        \
!       region3  |         \ region1
!                |          \
!                | region0   \
!                |            \
!                |             \ P1
! ---------------*--------------*------->s
!                |P0             \
!  region4       |    region5     \ region6


   integer function FindRegion( det, s, t ) result( region )
      implicit none
      !-arguments-----------------------------------------------------------------------------------------
      real(kind=rp), intent(in) :: det, s, t
      
      if( s+t <= det ) then
         if( s < 0 ) then
            if( t < 0 ) then
               region = 4
            else
               region = 3
            end if
         elseif( t < 0 ) then
            region = 5
         else
            region = 0
         end if
      else
         if( s < 0 ) then
            region = 2
         elseif ( t < 0 ) then
            region = 6
         else
            region = 1
         end if
      end if 
      
   end function FindRegion
   
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine computes the minumum distance from x (global ref. frame) to the body. First,
! the box (tree) where x lies is found, then the min distance between the objects inside the tree and
! the point is computed. If the sphere whose radius is minDist, is enclosed in the box, the subroutine stops.
! If the latter condition is not true, all the boxes that intersects the sphere are checked in the same way as
! the initial one, getting new_minDist. If a lower distance is found, minDist is updated.
!  ------------------------------------------------
   
   subroutine MinimumDistance( xP, root, minDist, normal )
   
      implicit none
      !-arguments---------------------------------------------
      real(kind=rp), dimension(:),    intent(in)    :: xP
      type(KDtree),                   intent(inout) :: root
      real(kind=rp),                  intent(inout) :: minDist
      real(kind=rp), dimension(NDIM), intent(out)   :: normal 
      !-local-variables---------------------------------------
      real(kind=rp), dimension(NDIM) :: IntersPoint, new_IntersPoint, &
                                        dsvec, x, Point, P, IP,       &
                                        IntersectionPoint
      logical                        :: Inside
      type(KDtree), pointer          :: tree, tmpTree
      real(kind=rp)                  :: Dist, New_minDist, Radius, ds
      integer                        :: i, index, LeafIndex
      
      minDist = huge(1.0_RP)

      call OBB(root% STLNum)% ChangeRefFrame( xP, 'local', Point ) 

      call root% FindLeaf( Point, tree )

      LeafIndex = tree% index

      do i = 1, tree% NumOfObjs
         index = tree% ObjsIndeces(i)
         call MinumumPointTriDistance( Point, root% ObjectsList(index)% vertices(1)% coords, &
                                       root% ObjectsList(index)% vertices(2)% coords,        &
                                       root% ObjectsList(index)% vertices(3)% coords, Dist,  &
                                       IntersPoint                                           )
          
         if( Dist .lt. minDist ) then
            minDist           = Dist
            IntersectionPoint = IntersPoint              
         end if
      end do   
 
      if( tree% NumOfObjs .gt. 0 ) then
      ! Check the sphere
      !-----------------
         Radius = sqrt(minDist)
         Inside = CheckHypersphere( tree, Point, Radius )
      else
         dsvec(1) = abs(tree% parent% vertices(1,7)-tree% parent% vertices(1,1))
         dsvec(2) = abs(tree% parent% vertices(2,7)-tree% parent% vertices(2,1))
         dsvec(3) = abs(tree% parent% vertices(3,7)-tree% parent% vertices(3,1))  
         ds = maxval(dsvec)
         Radius = 0.5_RP*ds
         Inside = .true.
      end if

      nullify(tree)

      if( Inside ) then
         New_minDist = huge(1.0_RP)
         call MinimumDistOtherBoxes( Point, root, root, Radius, LeafIndex, New_minDist, New_IntersPoint )
         if( New_minDist .lt. minDist ) then 
            minDist           = New_minDist
            IntersectionPoint = New_IntersPoint
         end if
      end if      
      
      call OBB(root% STLNum)% ChangeRefFrame( Point, 'global', P )
      call OBB(root% STLNum)% ChangeRefFrame( IntersectionPoint, 'global', IP )

      normal = (P - IP)/norm2(P - IP)
   
   end subroutine MinimumDistance
   
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This subroutine computes New_minDist 
!  ------------------------------------------------
   recursive subroutine MinimumDistOtherBoxes( Point, root, tree, Radius, LeafIndex, New_minDist, New_IntersPoint )
                                                                              
      implicit none
      !-arguments----------------------------------------------------------
      real(kind=rp), dimension(:),   intent(in)    :: Point
      type(KDtree),                  intent(in)    :: root
      type(KDtree),                  intent(inout) :: tree
      real(kind=rp),                 intent(in)    :: Radius
      integer,                       intent(in)    :: LeafIndex
      real(kind=rp),                 intent(inout) :: New_minDist
      real(kind=rp), dimension(:),   intent(inout) :: New_IntersPoint 
      !-local-variables---------------------------------------------------
      real(kind=rp)                  :: Dist
      logical                        :: Intersect
      real(kind=rp), dimension(NDIM) :: IntersPoint
      integer                        :: i, index
      
      Intersect = BoxIntersectSphere( Radius, Point, tree% vertices )
      
      if( Intersect ) then
         if( tree% isLast ) then
            if( LeafIndex .ne. tree% index ) then 
               do i = 1, tree% NumOfObjs
                  index = tree% ObjsIndeces(i)
                  call MinumumPointTriDistance( Point, root% ObjectsList(index)% vertices(1)% coords, &
                                                root% ObjectsList(index)% vertices(2)% coords,        &
                                                root% ObjectsList(index)% vertices(3)% coords, Dist,  &
                                                IntersPoint                                           )
                  if( Dist .lt. New_minDist ) then
                     New_minDist     = Dist
                     New_IntersPoint = IntersPoint
                  end if
               end do    
            end if
         else
            call MinimumDistOtherBoxes( Point, root, tree% child_L,     & 
                                        Radius, LeafIndex, New_minDist, &
                                        New_IntersPoint                 )    
            call MinimumDistOtherBoxes( Point, root, tree% child_R,     &
                                        Radius, LeafIndex, New_minDist, &
                                        New_IntersPoint                 )
         end if
      end if
          
   end subroutine MinimumDistOtherBoxes
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  ---------------------------------------------------------------------------------------------------
! This subroutine checks if a the sphere whit radius minDist is inside the box or not. If it's not
! a/ circle/s is/are computed. Each circle is the base of a cylinder used to find possible boxes that 
! can intersect the sphere. 
!  ---------------------------------------------------------------------------------------------------

   logical function CheckHypersphere( tree, Point, minDist) result( Inside )
   
      implicit none
      !-arguments-----------------------------------------------------
      real(kind=rp), dimension(:),       intent(in)    :: Point
      type(KDtree), target,              intent(inout) :: tree
      real(kind=rp),                     intent(inout) :: minDist
      
      Inside = BoxIntersectSphere( minDist, Point, tree% vertices )
   
   end function CheckHypersphere
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  ---------------------------------------------------------
! This function gives true if a circle interscts a rectangle
!  ---------------------------------------------------------

   logical function CircleRectIntersection( RectLength, RectWidth, RectCenter, &
                                            Center, Radius ) result( Intersect )

      implicit none
      !-arguments---------------------------------------------------------------
      real(kind=rp), dimension(:), intent(in) :: RectCenter, Center
      real(kind=rp),               intent(in) :: RectLength, RectWidth, Radius 
      !-local-variables---------------------------------------------------------
      real(kind=rp), dimension(2) :: d

      Intersect = .false.
      
      d = Center - RectCenter

      If( abs(d(1)) .gt. 0.5_RP*RectLength + Radius ) return
      If( abs(d(2)) .gt. 0.5_RP*RectWidth + Radius ) return

      if( abs(d(1)) .le. 0.5_RP*RectLength .and. &
          abs(d(2)) .le. 0.5_RP*RectWidth ) then
         Intersect = .true.
         return
      end if
      
      if( POW2(abs(d(1)) - 0.5_RP*RectLength) + POW2(abs(d(2)) - 0.5_RP*RectWidth) .le. &
          POW2(Radius) ) Intersect = .true.

   end function CircleRectIntersection
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This function computes the distance between a point and a box
!  ------------------------------------------------
   
   real(kind=rp) function PointBoxDistance( Point, vertices ) result( sqrDist )
   
      implicit none
      !-arguments--------------------------------------------
      real(kind=rp), dimension(:),   intent(in) :: Point
      real(kind=rp), dimension(:,:), intent(in) :: vertices
      !-local-variables--------------------------------------
      integer, dimension(NDIM) :: ind
      real(kind=rp)            :: v
      integer                  :: i
   
      sqrDist = 0.0_RP
   
      ind(1) = 2 ! x-axis
      ind(2) = 4 ! y-axis
      ind(3) = 5 ! z-axis
   
      ! if the point's x-coordinate (xp) is less than the x-min coord. of the box, then
      ! the minimum distance on the x-dir is the distance between the xp and x-min, or
      ! vertices(i,1) - x-coordinate. The following loop avoids the computation of the point
      ! with the minimum distance from P. If the point P is inside the box, sqrDist = 0.  
      !--------------------------------------------------------------------------------------      
      do i = 1, NDIM
         v = Point(i)
         if( v .lt. vertices(i,1) )      sqrDist = sqrDist + (vertices(i,1)-v)*(vertices(i,1)-v)
         if( v .gt. vertices(i,ind(i)) ) sqrDist = sqrDist + (v-vertices(i,ind(i)))*(v-vertices(i,ind(i)))
      end do
   
   end function PointBoxDistance
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
! This function checks if a sphere intersects a box 
!  ------------------------------------------------
   
   logical function BoxIntersectSphere( Radius, Center, vertices ) result( Intersect )
   
      implicit none
      !-arguments-------------------------------------------      
      real(kind=rp), dimension(:,:), intent(in) :: vertices
      real(kind=rp), dimension(:),   intent(in) :: Center
      real(kind=rp),                 intent(in) :: Radius
      !-local-variables-------------------------------------
      real(kind=rp) :: sqrDist
      
      Intersect = .false.
      sqrDist = PointBoxDistance( Center, vertices )
      
      if( sqrDist .le. POW2(Radius) ) Intersect = .true.
      
   end function BoxIntersectSphere
   

   subroutine MinimumDistancePoints( Point, root, BandRegion, minDist, LowerBound, actualIndex, &
                                     PointsIndex, forcingPointsMask, elements                   )
      use ElementClass
      implicit none
      !-arguments-------------------------------------------------------
      real(kind=rp),  dimension(:), intent(in)    :: Point
      type(KDtree),                 intent(inout) :: root
      type(MPI_M_Points_type),      intent(in)    :: BandRegion
      real(kind=rp),                intent(inout) :: minDist, LowerBound
      integer,                      intent(in)    :: actualIndex
      integer,        dimension(:), intent(inout) :: PointsIndex
      logical,        optional,     intent(in)    :: forcingPointsMask
      type(element),  optional,     intent(in)    :: elements(:)
      !-local-variables-------------------------------------------------
      real(kind=rp), dimension(NDIM) :: BandPoint, dsvec
      logical                        :: Inside, found, FPmask
      type(KDtree), pointer          :: tree
      real(kind=rp)                  :: sqrDist, New_minsqrDist, Radius, &
                                        ds
      integer                        :: i, LeafIndex, temp_index, eID,   &  
                                        loc_pos(NDIM)
      
      ! ref frame = global
      if( present(forcingPointsMask) ) then
         FPmask = forcingPointsMask
      else
         FPmask = .false.
      end if
      
      minDist = huge(1.0_RP)

      call root% FindLeaf( Point, tree ) 
      
      LeafIndex = tree% index

      do i = 1, tree% NumOfObjs
      
         if( FPMask ) then
            eID     = BandRegion% x(tree% ObjsIndeces(i))% element_index
            loc_pos = BandRegion% x(tree% ObjsIndeces(i))% local_Position 
            if( elements(eID)% isForcingPoint(loc_pos(1),loc_pos(2),loc_pos(3)) ) cycle 
         end if
         
         BandPoint = BandRegion% x(tree% ObjsIndeces(i))% coords
         
         sqrDist = POW2(norm2(Point - BandPoint))

         if( sqrDist .lt. minDist .and. sqrDist .gt. LowerBound .or. &
             AlmostEqual(sqrDist, LowerBound) ) then

             found = .false.
             if( any(PointsIndex .eq. tree% ObjsIndeces(i)) )then
                found = .true.
             end if

            if( .not. found ) then
               minDist = sqrDist
               PointsIndex(actualIndex) = tree% ObjsIndeces(i)
            end if           

         end if
         
      end do

      if( tree% NumOfObjs .gt. 0 ) then
      ! Check the sphere
      !-----------------
         Radius = sqrt(minDist)
         Inside = CheckHypersphere( tree, Point, Radius )
      else 
         dsvec(1) = abs(tree% parent% vertices(1,7)-tree% parent% vertices(1,1))
         dsvec(2) = abs(tree% parent% vertices(2,7)-tree% parent% vertices(2,1))
         dsvec(3) = abs(tree% parent% vertices(3,7)-tree% parent% vertices(3,1))  
         ds = maxval(dsvec)
         Radius = 0.5_RP*ds
         Inside = .true.
      end if
      
      nullify(tree)

      if( Inside ) then
         New_minsqrDist = huge(1.0_RP)
         if( present(elements) ) then
            call MinimumDistOtherBoxesPoints( Point, root, BandRegion, Radius, &
                                              New_minsqrDist, LowerBound,      &
                                              PointsIndex, LeafIndex,          &
                                              temp_index, FPmask, elements     )
         else
            call MinimumDistOtherBoxesPoints( Point, root, BandRegion, Radius, &
                                              New_minsqrDist, LowerBound,      &
                                              PointsIndex, LeafIndex,          &
                                              temp_index, FPmask               )         
         
         end if
         if( New_minsqrDist .le. minDist ) then
            minDist = New_minsqrDist; PointsIndex(actualIndex) = temp_index  
         end if
      end if    

      minDist = sqrt(minDist)
      
   end subroutine MinimumDistancePoints
   
   recursive subroutine MinimumDistOtherBoxesPoints( Point, tree, BandRegion, Radius, &
                                                     New_minsqrDist, LowerBound,      &
                                                     PointsIndex, LeafIndex,          &
                                                     temp_index, forcingPointsMask,   &
                                                     elements                         )
      use elementClass                                                                 
      implicit none
      !-arguments----------------------------------------------------------
      real(kind=rp), dimension(:),   intent(in)    :: Point
      type(KDtree),                  intent(inout) :: tree
      type(MPI_M_Points_type),       intent(in)    :: BandRegion
      real(kind=rp),                 intent(in)    :: Radius, LowerBound
      real(kind=rp),                 intent(inout) :: New_minsqrDist
      integer,                       intent(inout) :: temp_index
      integer,         dimension(:), intent(in)    :: PointsIndex
      integer,                       intent(in)    :: LeafIndex
      logical,                       intent(in)    :: forcingPointsMask
      type(element),   optional,     intent(in)    :: elements(:)
      !-local-variables---------------------------------------------------
      real(kind=rp)                  :: sqrDist, BandPoint(NDIM)
      logical                        :: Intersect, found
      integer                        :: i, eID, loc_pos(NDIM)
      
      Intersect = BoxIntersectSphere( Radius, Point, tree% vertices )
      
      if( Intersect ) then
         if( tree% isLast ) then
            if( LeafIndex .ne. tree% index ) then
               do i = 1, tree% NumOfObjs
               
                  if( forcingPointsMask ) then
                     eID     = BandRegion% x(tree% ObjsIndeces(i))% element_index
                     loc_pos = BandRegion% x(tree% ObjsIndeces(i))% local_Position 
                     if( elements(eID)% isForcingPoint(loc_pos(1),loc_pos(2),loc_pos(3)) ) cycle
                  end if
               
                  BandPoint = BandRegion% x(tree% ObjsIndeces(i))% coords
               
                  sqrDist = POW2(norm2(Point - BandPoint))
                  
                  if( sqrDist .lt. New_minsqrDist .and. sqrDist .gt. LowerBound .or. &
                      AlmostEqual(sqrDist, LowerBound) )then
                      
                      found = .false. 
                      if( any(PointsIndex .eq. tree% ObjsIndeces(i)) ) then
                         found = .true.
                      end if

                     if( .not. found ) then
                        New_minsqrDist = sqrDist
                        temp_index = tree% ObjsIndeces(i) 
                     end if
                     
                  end if
               end do    
            end if
         else
            if( present(elements) ) then
               call MinimumDistOtherBoxesPoints( Point, tree% child_L, BandRegion, &
                                                 Radius, New_minsqrDist,           &
                                                 LowerBound, PointsIndex,          &
                                                 LeafIndex, temp_index,            &
                                                 forcingPointsMask,                &
                                                 elements                          )    
               call MinimumDistOtherBoxesPoints( Point, tree% child_R, BandRegion, &
                                                 Radius, New_minsqrDist,           &
                                                 LowerBound, PointsIndex,          &
                                                 LeafIndex, temp_index,            &
                                                 forcingPointsMask,                &
                                                 elements                          )
            else
               call MinimumDistOtherBoxesPoints( Point, tree% child_L, BandRegion, &
                                                 Radius, New_minsqrDist,           &
                                                 LowerBound, PointsIndex,          &
                                                 LeafIndex, temp_index,            &
                                                 forcingPointsMask                 )
               call MinimumDistOtherBoxesPoints( Point, tree% child_R, BandRegion, &
                                                 Radius, New_minsqrDist,           &
                                                 LowerBound, PointsIndex,          &
                                                 LeafIndex, temp_index,            &
                                                 forcingPointsMask                 )           
            
            end if 
         end if
      end if
          
   end subroutine MinimumDistOtherBoxesPoints

   subroutine GetIDW_value( Point, BandRegion, normal, Q, PointsIndex, value )
      use PhysicsStorage
      use MappedGeometryClass
      implicit none
      !-arguments---------------------------------------------------------
      real(kind=RP),           dimension(:),     intent(in)  :: Point, normal
      type(MPI_M_Points_type),                   intent(in)  :: BandRegion
      real(kind=RP),           dimension(:,:),   intent(in)  :: Q
      integer,                 dimension(:),     intent(in)  :: PointsIndex
      real(kind=RP),           dimension(NCONS), intent(out) :: value
      !-local-variables---------------------------------------------------
      real(kind=RP) :: DistanceNormal(NDIM), d2, d1, distanceSqr, sqrd1, &
                       num(NCONS), den
      real(kind=RP), parameter :: eps = 1.0d-10
      integer       :: i, k
   
      num = 0.0_RP; den = 0.0_RP
   
      do i = 1, size(PointsIndex)
                          
         DistanceNormal = BandRegion% x(PointsIndex(i))% coords - Point
         d2 = dot_product(DistanceNormal, normal)
               
         distanceSqr = 0.0_RP
         do k = 1, NDIM
            distanceSqr = distanceSqr + POW2(BandRegion% x(PointsIndex(i))% coords(k) - Point(k))
         end do

         sqrd1 = distanceSqr - POW2(d2)
          
         d1 = sqrt(abs(sqrd1)) 
               
         num = num + Q(:,i)/(d1+eps)
         den = den + 1.0_RP/(d1+eps)

      end do 
   
      value = num/den
   
   end subroutine GetIDW_value      
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!                                          TURBULENCE 
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!  
!  -------------------------------------------------
!  This subroutine performes iterations so that the 
!  forcing point lies inside the log region 
!  ------------------------------------------------

   subroutine IBM_GetImagePoint_nearest( this, elements )
      use MappedGeometryClass
      use PhysicsStorage
      use MPI_Process_Info
      use ElementClass
      implicit none
      !-arguments---------------------------------------------------------------------------
      class(IBM_type), intent(inout) :: this
      type(element),   intent(inout) :: elements(:)
      !-local-variables----------------------------------------------------------------------
      real(kind=RP) :: Dist, LowerBound, ImagePoint(NDIM)
      integer       :: eID, i, j, k, m, n, Nearest_Points(this% KDtree_n_of_interPoints)
      
      if( allocated(this% ImagePoint_NearestPoints) ) deallocate(this% ImagePoint_NearestPoints)
      allocate(this% ImagePoint_NearestPoints(this% KDtree_n_of_interPoints,this% NumOfForcingPoints))
      
      this% ImagePoint_NearestPoints = 0
      m = 0  

!$omp parallel
!$omp do schedule(runtime) private(i,j,k,Dist,LowerBound,n,ImagePoint,Nearest_Points)
      do eID = 1, size(elements)             
         do k = 0, elements(eID)% Nxyz(3); do j = 0, elements(eID)% Nxyz(2); do i = 0, elements(eID)% Nxyz(1)
            if( elements(eID)% isForcingPoint(i,j,k) ) then
            
               Dist = this% IP_distance - elements(eID)% geom% dWall(i,j,k)
            
               ImagePoint = elements(eID)% geom% x(:,i,j,k) + Dist * elements(eID)% geom% normal(:,i,j,k) 
              
               LowerBound  = -huge(1.0_RP)
               do n = 1, this% KDtree_n_of_interPoints 
                  call MinimumDistancePoints( ImagePoint, this% rootPoints, this% BandRegion, &
                                              Dist, LowerBound, n, Nearest_Points,            &
                                              .true., elements                                ) 
                  LowerBound = POW2(Dist)
               end do
!$omp critical
               m = m + 1
               elements(eID)% IP_index(i,j,k) = m
               this% ImagePoint_NearestPoints(:,m) = Nearest_Points
               if(any(this% ImagePoint_NearestPoints(:,m) .eq. 0) ) then
                  print *, "Can't fine nearest point"
                  error stop
               end if
!$omp end critical
            end if
         end do               ; end do               ; end do                 
      end do
!$omp end do
!$omp end parallel
   end subroutine IBM_GetImagePoint_nearest
   
#if defined(NAVIERSTOKES)  
   subroutine ForcingPointState( Q_IP, y_IP, y_FP, normal, Q_FP )
      use PhysicsStorage
      use WallFunctionDefinitions
      use WallFunctionBC
      use VariableConversion
      use SpallartAlmarasTurbulence
      use FluidData   
      implicit none
      
      real(kind=rp), dimension(:), intent(in)    :: Q_IP, normal
      real(kind=rp),               intent(in)    :: y_IP, Y_FP
      real(kind=rp), dimension(:), intent(inout) :: Q_FP
      !-local-variables--------------------------------------------------------
      real(kind=rp)              :: nu_IP, nu_FP, mu_IP, mu_FP, T_IP, T_FP, &
                                    u_tau, u_IP(NDIM), u_IPt, u_FPt,        &
                                    u_FPn, u_FP(NDIM), u_IP_t(NDIM),        &
                                    tangent(NDIM), yplus_FP, nu_tilde,      &
                                    kappa_IP, kappa_FP

#if defined(SPALARTALMARAS)
      real(kind=rp)              :: Dump, chi, fv1, nu_t
#endif
     
      call get_laminar_mu_kappa(Q_IP,mu_IP,kappa_IP)      
      nu_IP = mu_IP/Q_IP(IRHO)
            
      u_IP   = Q_IP(IRHOU:IRHOW)/Q_IP(IRHO)
      u_IP_t = u_IP - ( dot_product(u_IP,normal) * normal )
   
      tangent = u_IP_t/norm2(u_IP_t)   
   
      u_IPt = dot_product(u_IP,tangent)
            
      if( almostEqual(u_IPt,0.0_RP) ) return

      u_tau = u_tau_f(u_IPt, y_IP, nu_IP, u_tau0=1.0_RP)

      call get_laminar_mu_kappa(Q_FP,mu_FP,kappa_FP)
      nu_FP = mu_FP/Q_FP(IRHO)                 
         
      if( any(isnan(Q_IP)) ) write(*,*) 'Q_IP is nan'
      if( isnan(nu_IP) ) write(*,*) 'nu_IP is nan, rho_ip =', Q_IP(IRHO)
         
      u_FPt = u_plus_f( y_plus_f( y_FP, u_tau, nu_FP) ) * u_tau
      u_FPn = dot_product(u_IP,normal) * y_FP/y_IP
         
      u_FP = u_FPt*tangent + u_FPn*normal 
            
      T_FP = Temperature(Q_IP) + (dimensionless% Pr)**(1._RP/3._RP)/(2.0_RP*thermodynamics% cp) * (POW2(u_IPt) - POW2(u_FPt))
            
#if defined(SPALARTALMARAS)      
!~       yplus_FP = y_plus_f( y_FP, u_tau, nu_FP) 
         
!~       Dump     = 1.0_RP - exp(-yplus_FP/19.0_RP)      
      
!~       chi      = QuarticRealPositiveRoot( 1.0_RP, -kappa*u_tau*y_FP*Dump/nu_FP, 0.0_RP, 0.0_RP, -kappa*u_tau*y_FP*Dump/nu_FP*POW3(SAmodel% cv1) )
      
!~       nu_t = nu_FP * chi
      
      nu_T = kappa * u_tau * y_FP
      
#endif

      T_FP = refValues% T * T_FP

      Q_FP(IRHO)  = Pressure(Q_IP)*refvalues% p/(thermodynamics% R * T_FP)
      Q_FP(IRHO)  = Q_FP(IRHO)/refvalues% rho     
      Q_FP(IRHOU) = Q_FP(IRHO) * u_FP(1)
      Q_FP(IRHOV) = Q_FP(IRHO) * u_FP(2)
      Q_FP(IRHOW) = Q_FP(IRHO) * u_FP(3)
      Q_FP(IRHOE) = pressure(Q_IP)/( thermodynamics% gamma - 1._RP) + 0.5_RP*Q_FP(IRHO)*(POW2(u_FP(1)) + POW2(u_FP(2)) + POW2(u_FP(3)))
#if defined(SPALARTALMARAS)
      Q_FP(IRHOTHETA) = Q_FP(IRHO) * nu_t
#endif              
   end subroutine ForcingPointState
#endif  
   
   
#if defined(NAVIERSTOKES)    
   subroutine IBM_SourceTermTurbulence( this, eID, Q, Qdot, Qbp, normal, x, dWall, IP_index, TurbulenceSource )
      use PhysicsStorage
      implicit none
      !-arguments--------------------------------------------------------------
      class(IBM_type), intent(inout) :: this
      integer,         intent(in)    :: eID, IP_index
      real(kind=rp),   intent(in)    :: Q(:), Qdot(:), Qbp(:,:), x(:), dWall, normal(:)
      real(kind=rp),   intent(inout) :: TurbulenceSource(NCONS)
      !-local-variables--------------------------------------------------------
      real(kind=rp)              :: Q_IP(NCONS), Q_FP(NCONS), Dist, ImagePoint(NDIM)
      
      Dist = this% IP_distance - dWall
            
      ImagePoint = x + Dist * normal
        
      call GetIDW_value( ImagePoint, this% BandRegion, normal,              &
                         Qbp(:,this% ImagePoint_NearestPoints(:,IP_index)), &
                         this% ImagePoint_NearestPoints(:,IP_index),        &
                         Q_IP                                               )
                                 
      Q_FP = Q
                                 
      call ForcingPointState( Q_IP, this% IP_Distance, dWall, normal, Q_FP )

      TurbulenceSource =  1.0_RP/this% penalization(eID) * (Q_FP - Q)
   
   end subroutine IBM_SourceTermTurbulence
#endif 


! estimate the yplus for a flat plate 
   
   real(kind=RP) function InitializeDistance( y_plus ) result( y )
      use FluidData
      implicit none
      !-arguments
      real(kind=rp), intent(in) :: y_plus
      !-local-varirables-------------------------
      real(kind=RP) :: nu, Cf
      
#if defined(NAVIERSTOKES)     
      nu = refValues% mu / refValues% rho
#endif

      Cf = Estimate_Cf()
#if defined(NAVIERSTOKES)      
      y = sqrt(2.0_RP) * y_plus * nu /(refValues% V * sqrt(Cf)) 
#endif  
   end function InitializeDistance
   
   real(kind=RP) function Estimate_Cf() result( Cf )
      use FluidData
      implicit none
      
      Cf = 0.0_RP
      
#if defined(NAVIERSTOKES)   
     ! Schlichting, Hermann (1979), Boundary Layer Theory... ok if Re < 10^9
!~       Cf = (2.0_RP*log10(dimensionless% Re) - 0.65_RP)**(-2.3_RP)     
      Cf = 0.058_RP * dimensionless% Re**(-0.2)     
#endif        
   end function Estimate_Cf
   
   real(kind=RP) function Estimate_u_tau( ) result( u_tau )
      use FluidData
      implicit none
      !-local-variables--------------------------------
      real(kind=RP) :: Cf
      
      u_tau = 0.0_RP
      
      Cf = Estimate_Cf()
#if defined(NAVIERSTOKES) 
      u_tau = sqrt( 0.5_RP * POW2(refValues% V) * Cf ) !sqrt(tau_w/rho)
#endif      
   end function Estimate_u_tau
   
   
   real(kind=RP) function GetEstimated_y( y_plus, nu, u_tau ) result( y )
      use PhysicsStorage 
      implicit none
      !-arguments--------------------------------------
      real(kind=RP), intent(in) :: y_plus, nu, u_tau
   
      y = y_plus * nu / u_tau
      
      y = y/Lref 
   
   end function GetEstimated_y


   real(kind=RP) function QuarticRealPositiveRoot(a0, b0, c0, d0, e0) result( value )
      !https://math.stackexchange.com/questions/115631/quartic-equation-solution-and-conditions-for-real-roots
      implicit none
      !-arguments--------------------------------------
      real(kind=rp), intent(in) :: a0, b0, c0, d0, e0
      !-local-variables--------------------------------
      real(kind=rp) :: a, b, c, d, x, m2, m, n,  &
                       alpha, beta, psi,         & 
                       solution(2)

      a = b0/a0; b = c0/a0; c = d0/a0; d = e0/a0
      
      x = cubic_poly( 1.0_RP, -b, a*c-4.0_RP*d, d*(4.0_RP*b-POW2(a)) - POW2(c) )

      m2 = 0.25_RP*POW2(a)-b+x
      
      if( m2 .gt. 0.0_RP ) then
         m = sqrt(m2)
         n = (a*x-2.0_RP*c)/(4.0_RP*m)
      elseif( almostEqual(m2,0.0_RP) ) then
         m = 0.0_RP
         n = sqrt(0.25_RP*POW2(x) - d)
      else
         print *, 'no real roots found'
         error stop
      end if
      
      alpha = 0.5_RP*POW2(a) - x - b
      
      beta = 4.0_RP * n - a * m
      
      if( (alpha + beta) .ge. 0.0_RP ) then
         psi = sqrt(alpha+beta)
         solution(1) = -0.25_RP*a + 0.5_RP*m + 0.5_RP*psi
         solution(2) = -0.25_RP*a + 0.5_RP*m - 0.5_RP*psi
      elseif( (alpha - beta) .ge. 0.0_RP ) then
         psi = sqrt(alpha-beta)
         solution(1) = -0.25_RP*a - 0.5_RP*m + 0.5_RP*psi
         solution(2) = -0.25_RP*a - 0.5_RP*m - 0.5_RP*psi        
      end if
      !take positive solution
      value = maxval(solution) 

   end function QuarticRealPositiveRoot
   
   real(KIND=rp) function cubic_poly( a0, b0, c0, d0 ) result( x )
   
      implicit none
      !-arguments----------------------------------------------
      real(kind=RP),      intent(in) :: a0, b0, c0, d0
      !-local-variables----------------------------------------
      real(kind=RP) :: a, b, c, Jac, x0, q, r
      real(kind=RP) :: theta, m, s, txx
      integer :: i
      
      
      a = b0/a0; b = c0/a0; c = d0/a0
           
      x0 = 1.0_RP
      
      do  i = 1, 100
      
         Jac = 3.0_RP*a0*POW2(x0) + 2.0_RP*b0*x0 + c0
   
         x = x0 - eval_cubic(a0, b0, c0, d0, x0)/Jac
   
         if( abs(x-x0) .lt. 1.0E-10_RP ) return
   
         x0 = x
   
      end do
    
      print *, "cubic_poly:: Newtown method doesn't coverge"
      error stop
    
   end function cubic_poly

   real(kind=RP) function eval_cubic( a0, b0, c0, d0, x0 ) result( value )
   
      implicit none
      !-arguments----------------------------------
      real(kind=RP), intent(in) :: a0, b0, c0, d0, x0
      
      value = a0*POW3(x0) + b0*POW2(x0) + c0*x0 + d0
   
   end function eval_cubic

end module IBMClass
