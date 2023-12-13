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

   type :: Integral_t

      logical              :: ListComputed = .false., &
                              compute      = .false., &
                              constructed  = .false.
   end type

   type :: bandRegion_t

      type(IBMPoints),  allocatable :: IBMmask(:)

      contains
         procedure :: plot => bandRegion_plot

   end type bandRegion_t

   type :: IBM_type

      type(STLfile),              allocatable :: stl(:)
      type(KDtree),               allocatable :: root(:), rootDistance(:)
      type(bandRegion_t),         allocatable :: BandRegion(:), ImagePoints(:)
      type(Integral_t),           allocatable :: Integral(:)
      character(len=LINE_LENGTH), allocatable :: STLfilename(:)
      character(len=LINE_LENGTH)              :: filename
      logical                                 :: plotOBB              = .false., &
                                                 plotKDtree           = .false., &
                                                 active               = .false., &
                                                 TimePenal            = .false., &
                                                 semiImplicit         = .false., &
                                                 ComputeBandRegion    = .false., &
                                                 plotBandPoints       = .false., &
                                                 plotMask             = .false., &
                                                 ComputeInterpolation = .false., &
                                                 Wallfunction         = .false., &
                                                 ComputeDistance      = .false., &
                                                 AAB                  = .false., &
                                                 HO_IBM               = .false.
      real(kind=rp)                           :: eta, BandRegionCoeff, IP_Distance = 0.0_RP, &
                                                 y_plus_target, minCOORDS, maxCOORDS,        &
                                                 penalCoeff
      real(kind=rp),              allocatable :: penalization(:)
      integer                                 :: KDtree_Min_n_of_Objs, NumOfInterPoints,     &
                                                 n_of_INpoints,  rank, lvl = 0, NumOfSTL,    &
                                                 NumOfForcingPoints, Clipaxis = 0,           &
                                                 Nx, Ny, Nz, LocClipAxis = 0,                &
                                                 InterpolationType, NumOfMaskObjs = 0
      integer,                    allocatable :: ImagePoint_NearestPoints(:,:)
      type(IBMPoints),            allocatable :: IBMmask(:)

      contains
         procedure :: read_info                           => IBM_read_info
         procedure :: construct                           => IBM_construct
         procedure :: constructMask                       => IBM_constructMask
         procedure :: constructSTL_KDtree                 => IBM_constructSTL_KDtree
         procedure :: CheckPoint                          => IBM_checkPoint
         procedure :: constructBandRegion                 => IBM_constructBandRegion
         procedure :: build                               => IBM_build
         procedure :: SetPolynomialOrder                  => IBM_SetPolynomialOrder
         procedure :: GetMask                             => IBM_GetMask
         procedure :: MPI_sendSTLpartitions               => IBM_MPI_sendSTLpartitions
         procedure :: MPI_OBB                             => IBM_MPI_OBB
         procedure :: GetForcingPointsGeom                => IBM_GetForcingPointsGeom
         procedure :: GetInfo                             => IBM_GetInfo
         procedure :: minDistance                         => IBM_minDistance
         procedure :: SourceTerm                          => IBM_SourceTerm
         procedure :: ComputeIBMWallDistance              => IBM_ComputeIBMWallDistance
         procedure :: SemiImplicitCorrection              => IBM_SemiImplicitCorrection
         procedure :: GetImagePoint_nearest               => IBM_GetImagePoint_nearest
         procedure :: GetBandRegionStates                 => IBM_GetBandRegionStates
         procedure :: GetImagePointsStates                => IBM_GetImagePointsStates
         procedure :: GetDomainExtreme                    => IBM_GetDomainExtreme
         procedure :: semiImplicitShiftJacobian           => IBM_semiImplicitShiftJacobian
         procedure :: semiImplicitJacobian                => IBM_semiImplicitJacobian
         procedure :: semiImplicitTurbulenceJacobian      => IBM_semiImplicitTurbulenceJacobian
         procedure :: GetSemiImplicitStep                 => IBM_GetSemiImplicitStep
         procedure :: GetSemiImplicitStepTurbulence       => IBM_GetSemiImplicitStepTurbulence
         procedure :: copy                                => IBM_copy
         procedure :: MoveBody                            => IBM_MoveBody
         procedure :: CleanMask                           => IBM_CleanMask
         procedure :: Describe                            => IBM_Describe
         procedure :: plot_Mask                           => IBM_plot_Mask
         procedure :: Destruct                            => IBM_Destruct
         procedure :: DestroyKDtree                       => IBM_DestroyKDtree
         procedure :: constructDistance_KDtree            => IBM_constructDistance_KDtree
         procedure :: MaskVelocity                        => IBM_MaskVelocity
         procedure :: HOmask                              => IBM_HOmask
         procedure :: GetStencil                          => IBM_GetStencil
         procedure :: HO_IBMstencilState                  => IBM_HO_IBMstencilState
         procedure :: MPI_GatherStancilState              => IBM_MPI_GatherStancilState
         procedure :: buildHOfaces                        => IBM_buildHOfaces
   end type

   public :: expCoeff, EXPONENTIAL

   real(kind=RP)      :: expCoeff
   integer, parameter :: EXPONENTIAL = 1, IDW = 2, POLYHARMONIC_SPLINE = 3, POLYNOMIAL = 4, MLS = 5

   contains
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!

   subroutine bandRegion_plot( this, filename )
      use MPI_Process_Info
      implicit none

      class(bandRegion_t), intent(inout) :: this
      character(len=*),    intent(in)    :: filename

      integer :: domains, i, funit, NumOfObjs
      logical :: add_Point = .false.

      NumOfObjs = sum(this% IBMmask(1:MPI_Process% nProcs)% NumOfObjs)

      if( NumOfObjs .eq. 0 ) return

      call TecFileHeader( 'IBM/'//trim(filename), 'Points', NumOfObjs/2+mod(NumOfObjs,2),2,1, funit, 'POINT' )

      if( mod(NumOfObjs,2) .ne. 0 ) add_Point  = .true.

      do domains = 1, MPI_Process% nProcs
         do i = 1, this% IBMmask(domains)% NumOfObjs
            write(funit,'(3E13.5)')  this% IBMmask(domains)% x(i)% coords(IX), this% IBMmask(domains)% x(i)% coords(IY), this% IBMmask(domains)% x(i)% coords(IZ)
            if( add_Point ) then
               write(funit,'(3E13.5)') this% IBMmask(domains)% x(i)% coords(IX), this% IBMmask(domains)% x(i)% coords(IY), this% IBMmask(domains)% x(i)% coords(IZ)
               add_Point = .false.
            end if
         end do
      end do

      close(funit)

   end subroutine bandRegion_plot

!  -------------------------------------------------
!  Immersed boundary info
!  -------------------------------------------------
   subroutine IBM_read_info( this, controlVariables )
      use FileReadingUtilities
      implicit none

      class(IBM_type), intent(inout) :: this
      class(FTValueDictionary)       :: controlVariables

      call this% GetInfo( controlVariables )

   end subroutine IBM_read_info


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
      logical,       allocatable :: active_in, plotOBB_in,              &
                                    plotKDtree_in, semiImplicit_in,     &
                                    plotBandPoints_in, plotMask_in,     &
                                    BandRegion_in, Distance_in, AAB_in, &
                                    HO_in
      real(kind=rp), allocatable :: penalization_in, y_plus_target_in,  &
                                    BandRegionCoeff_in,                 &
                                    penalization_coeff_in
      integer,       allocatable :: n_of_Objs_in, n_of_interpoints_in,  &
                                    Nx_in, Ny_in, Nz_in, Clipaxis_in
      character(len=LINE_LENGTH) :: in_label, paramFile, name_in, tmp,  &
                                    InterpolationType_in
      real(kind=rp), allocatable :: coords(:)
      logical                    :: correct

      character(len=LINE_LENGTH), parameter :: NumberOfSTL = "number of stl"

!     Read block
!     **********
      write(in_label , '(A)') "#define ibm"
      call get_command_argument(1, paramFile)
      call readValueInRegion( trim( paramFile ), "name",                           name_in,               in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "active",                         active_in,             in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "penalization",                   penalization_in,       in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "penalization coeff",             penalization_coeff_in, in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "semi implicit",                  semiImplicit_in,       in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "nx",                             Nx_in,                 in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "ny",                             Ny_in,                 in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "nz",                             Nz_in,                 in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "target y plus",                  y_plus_target_in,      in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "number of objects",              n_of_Objs_in,          in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "number of interpolation points", n_of_interpoints_in,   in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "clip axis",                      Clipaxis_in,           in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "band region",                    BandRegion_in,         in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "distance",                       Distance_in,           in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "band region coeff",              BandRegionCoeff_in,    in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "aab",                            AAB_in,                in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "intepolation",                   InterpolationType_in,  in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "plot obb",                       plotOBB_in,            in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "plot kdtree" ,                   plotKDtree_in,         in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "plot mask",                      plotMask_in,           in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "plot band points",               plotBandPoints_in,     in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "high order",                     HO_in,                 in_label, "#end" )

      this% filename = trim(name_in)

      if( allocated(active_in) ) then
         this% active = active_in
      else
         this% active = .FALSE.
      end if

      if( .not. this% active) return

      if( allocated(semiImplicit_in) ) then
         this% semiImplicit = semiImplicit_in
      else
         this% semiImplicit = .FALSE.
      end if

      if( allocated(Nx_in) ) then
         this% Nx = Nx_in
      else
         this% Nx = 0
      end if

      if( allocated(Ny_in) ) then
         this% Ny = Ny_in
      else
         this% Ny = 0
      end if

      if( allocated(Nz_in) ) then
         this% Nz = Nz_in
      else
         this% Nz = 0
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
         this% TimePenal = .false.
      else
         this% eta = 0.1_RP
         this% TimePenal = .true.
      end if

      if( allocated(penalization_coeff_in) ) then
         this% penalCoeff = penalization_coeff_in
      else
         this% penalCoeff = 1.0_RP
      end if

      if( allocated(n_of_Objs_in) ) then
         this% KDtree_Min_n_of_Objs = n_of_Objs_in
      else
         this% KDtree_Min_n_of_Objs = 5
      end if

      if( allocated(n_of_interpoints_in) ) then
         this% NumOfInterPoints = n_of_interpoints_in
      else
         this% NumOfInterPoints = 15
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

      if( allocated(BandRegion_in) ) then
         this% ComputeBandRegion = BandRegion_in
      else
         this% ComputeBandRegion = .false.
      end if

      if( allocated(Distance_in) ) then
         this% ComputeDistance = Distance_in
      else
         this% ComputeDistance = .false.
      end if

      if( allocated(Clipaxis_in) ) then
         this% ClipAxis = Clipaxis_in
      else
         this% ClipAxis = 0
      endif

     if( controlVariables% containsKey("wall function") ) then
        this% Wallfunction = .true.
#if defined(NAVIERSTOKES)
        call Initialize_Wall_Function(controlVariables, correct)   !TO BE REMOVED
#endif
        if( allocated(y_plus_target_in) ) then
           this% y_plus_target = y_plus_target_in
        else
           this% y_plus_target = 50.0_RP
        end if

        this% ComputeBandRegion = .true.
        this% ComputeDistance   = .true.

     else
        this% Wallfunction = .false.
     end if

      if( allocated(BandRegionCoeff_in) ) then
         this% BandRegionCoeff = BandRegionCoeff_in
      else
         this% BandRegionCoeff = 2.0_RP
      end if

      if( controlVariables% containsKey(trim(NumberOfSTL)) ) then
         tmp = controlVariables% StringValueForKey(trim(NumberOfSTL),LINE_LENGTH)
         this% NumOfSTL = GetIntValue(tmp)
      else
         this% NumOfSTL = 1
      end if

      if( allocated(AAB_in) ) then
         this% AAB = AAB_in
      else
         this% AAB = .false.
      end if

      if( allocated(HO_in) ) then
         this% HO_IBM = HO_in
         if( this% HO_IBM ) this% ComputeDistance = .true.
      else
         this% HO_IBM = .false.
      end if

      select case(trim(InterpolationType_in))
         case("exp")
            this% InterpolationType =  EXPONENTIAL
         case("idw")
            this% InterpolationType =  IDW
         case("spline")
            this% InterpolationType =  POLYHARMONIC_SPLINE
         case("polynomial")
            this% InterpolationType =  POLYNOMIAL
         case("mls")
            this% InterpolationType =  MLS
         case default
            this% InterpolationType =  IDW
      end select

   end subroutine IBM_GetInfo
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------
!  KDtree procedures
!  -------------------------------------------------
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
      real(kind=RP)              :: axis(NDIM)
      integer                    :: STLNum, j, k, maxAxis, minAxis

      call this% describe()

      allocate( this% stl(this% NumOfSTL),        &
                OBB(this% NumOfSTL),              &
                this% root(this% NumOfSTL),       &
                this% integral(this% NumOfSTL),   &
                this% STLfilename(this% NumOfSTL) )

      if( this% ComputeBandRegion ) allocate( this% BandRegion(this% NumOfSTL)   )
      if( this% ComputeDistance   ) allocate( this% rootDistance(this% NumOfSTL) )
      if( this% Wallfunction      ) allocate( this% ImagePoints(this% NumOfSTL)  )

      do STLNum = 1, this% NumOfSTL
         write(MyString, '(i100)') STLNum
         if( STLNum .eq. 1 ) then
            filename = stlFileNameKey
         else
            filename = trim(stlFileNameKey)//trim(adjustl(MyString))
         end if
         this% STLfilename(STLNum) = controlVariables% stringValueForKey(trim(filename), requestedLength = LINE_LENGTH)
         OBB(STLNum)% filename     = this% STLfilename(STLNum)
         call STLfile_GetMotionInfo( this% stl(STLNum), this% STLfilename(STLNum), this% NumOfSTL )
         if( MPI_Process% isRoot ) then
            this% stl(STLNum)% show = .true.
            call this% stl(STLNum)% ReadTessellation( this% STLfilename(STLNum) )
            if( this% ClipAxis .ne. 0 ) then
               call this% stl(STLNum)% Clip( this% minCOORDS, this% maxCOORDS, this% ClipAxis, .true. )
            else
               call this% stl(STLNum)% describe(this% STLfilename(STLNum))
            end if
            call OBB(STLNum)% construct( this% stl(STLNum), this% plotOBB, this% AAB, .false.)
            OBB(STLNum)% maxAxis = OBB(STLNum)% GetMaxAxis()
            OBB(STLNum)% minAxis = OBB(STLNum)% GetMinAxis( OBB(STLNum)% maxAxis, this% ClipAxis )
         end if
#ifdef _HAS_MPI_
         call this% MPI_sendSTLpartitions( STLNum, OBB(STLNum)% maxAxis )
         call this% MPI_OBB( STLNum )
#endif
         call OBB(STLNum)% ChangeObjsRefFrame( this% stl(STLNum)% ObjectsList, LOCAL )
         call this% constructSTL_KDtree( STLNum )
      end do

   end subroutine IBM_Construct

   subroutine IBM_constructSTL_KDtree( this, STLNum )
      use MPI_Process_Info
      implicit none
      !-arguments-----------------------------------
      class(IBM_type), intent(inout) :: this
      integer,         intent(in)    :: STLNum
      !-local-variables-----------------------------
      real(kind=RP) :: vertices(NDIM,BOXVERTICES)

      this% root(STLNum)% STLNum       = STLNum
      this% root(STLNum)% which_KDtree = TRIANGLES_KDTREE_SAH

      vertices = OBB(STLNum)% LocVertices

      this% root(STLNum)% MaxAxis = OBB(STLNum)% maxAxis

      call this% root(STLNum)% construct( stl           = this% stl(STLNum),         &
                                          vertices      = vertices,                  &
                                          isPlot        = this% plotKDtree,          &
                                          Min_n_of_Objs = this% KDtree_Min_n_of_Objs )

      if( this% ComputeDistance ) call this% constructDistance_KDtree( STLNum )

   end subroutine IBM_constructSTL_KDtree

   subroutine IBM_constructDistance_KDtree( this, STLNum )
      use MPI_Process_Info
      implicit none
      !-arguments----------------------------------------
      class(IBM_type), intent(inout) :: this
      integer,         intent(in)    :: STLNum
      !-local-variables----------------------------------
      real(kind=RP) :: vertices(NDIM,BOXVERTICES)

      this% rootDistance(STLNum)% STLNum       = STLNum
      this% rootDistance(STLNum)% which_KDtree = TRIANGLES_KDTREE_MEDIAN

      vertices = OBB(STLNum)% LocVertices

      call this% rootDistance(STLNum)% construct( stl           = this% stl(STLNum),         &
                                                  Vertices      = vertices,                  &
                                                  isPlot        = .false.,                   &
                                                  Min_n_of_Objs = this% KDtree_Min_n_of_Objs )

   end subroutine IBM_constructDistance_KDtree
!
!  Copy a KD tree
!  ---------------
   subroutine IBM_copy( this, parent, lvl )

      implicit none
      !-arguments----------------------------------------------------------------
      class(IBM_type),         intent(inout) :: this
      type(IBM_type),  target, intent(in)    :: parent
      integer,                 intent(in)    :: lvl
      !-local-variables----------------------------------------------------------
      integer :: STLNum

      allocate( this% root(parent% NumOfSTL),       &
                this% STLfilename(parent% NumOfSTL) )

      if( parent% ComputeBandRegion ) then
         allocate( this% BandRegion(parent% NumOfSTL) )
      end if

      if( parent% ComputeDistance ) allocate(this% rootDistance(parent% NumOfSTL))

      do STLNum = 1, parent% NumOfSTL
         this% STLfilename(STLNum) = parent% STLfilename(STLNum)
         this% root(STLNum)        = parent% root(STLNum)
         if( parent% ComputeDistance ) this% rootDistance(STLNum) = parent% rootDistance(STLNum)
      end do

      this% ClipAxis = parent% ClipAxis

      this% lvl = lvl

   end subroutine IBM_copy
!
!  Destroy the KD tree
!  --------------------
   subroutine IBM_DestroyKDtree( this, isChild, DistanceKDtree )
      use MPI_Process_Info
      implicit none
      !-arguments--------------------------------------------------
      class(IBM_type),           intent(inout) :: this
      logical,                   intent(in)    :: isChild
      logical,         optional, intent(in)    :: DistanceKDtree
      !-local-variables--------------------------------------------
      integer :: STLNum

      if(this% ComputeDistance ) then
         do STLNum = 1, this% NumOfSTL
            call this% rootDistance(STLNum)% destruct( isChild )
         end do
         deallocate(this% rootDistance)
         return
      end if

      do STLNum = 1, this% NumOfSTL
         call this% root(STLNum)% destruct( isChild )
      end do

   end subroutine IBM_DestroyKDtree
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------
!  Mask procedures
!  -------------------------------------------------
  subroutine IBM_GetMask( this, elements, faces, no_of_DoFs, STLNum, iter )
      use MPI_Process_Info
      implicit none
      !-arguments-----------------------------------------------------------
      class(IBM_type), intent(inout) :: this
      type(element),   intent(inout) :: elements(:)
      type(face),      intent(inout) :: faces(:)
      integer,         intent(in)    :: no_of_DoFs, STLNum, iter

      call this% constructmask( elements, STLNum, no_of_DoFs, faces, iter )

   end subroutine IBM_GetMask
!
!
!  mask construction
!  -----------------------------------
   subroutine IBM_constructmask( this, elements, STLNum, no_of_DoFs, faces, iter )
      use MPI_Process_Info
      use ElementConnectivityDefinitions
      implicit none
      !-arguments----------------------------------------------------
      class(IBM_type),  intent(inout) :: this
      integer,          intent(in)    :: STLNum, no_of_DoFs, iter
      type(element),    intent(inout) :: elements(:)
      type(face),       intent(inout) :: faces(:)
      !-local-variables----------------------------------------------
      real(kind=RP) :: Point(NDIM)
      integer       :: eID, n, i, j, k, NumOfObjs, domains, domain, NInters
#ifdef _HAS_MPI_
      integer       :: ierr
#endif
      domain = MPI_Process% rank + 1

      allocate(this% IBMmask(MPI_Process% nProcs))

      call this% IBMmask(domain)% build(no_of_DoFs)

      this% IBMmask(domain)% NumOfObjs = 0

      do eID = 1, size(elements)
         associate( e => elements(eID) )
         call e% ConstructIBM( e% Nxyz(1), e% Nxyz(2), e% Nxyz(3), this% NumOfSTL )
         if( this% HO_IBM ) then
            do k = 1, NODES_PER_ELEMENT
               e% MaskCorners(k) = .false.
               if( OBB(STLNum)% isPointInside( coords = e% SurfInfo% corners(:,k) ) ) then
                  this% IBMmask(domain)% NumOfObjs = this% IBMmask(domain)% NumOfObjs + 1
                  this% IBMmask(domain)% x(this% IBMmask(domain)% NumOfObjs)% coords         = e% SurfInfo% corners(:,k)
                  this% IBMmask(domain)% x(this% IBMmask(domain)% NumOfObjs)% element_index  = eID
                  this% IBMmask(domain)% x(this% IBMmask(domain)% NumOfObjs)% local_Position = (/k,0,0/)
               end if
            end do
         else
            do i = 0, e% Nxyz(1); do j = 0, e% Nxyz(2); do k = 0, e% Nxyz(3)
               e% isInsideBody(i,j,k) = .false.
               if( OBB(STLNum)% isPointInside( coords = e% geom% x(:,i,j,k) ) ) then
                  this% IBMmask(domain)% NumOfObjs = this% IBMmask(domain)% NumOfObjs + 1
                  this% IBMmask(domain)% x(this% IBMmask(domain)% NumOfObjs)% coords         = e% geom% x(:,i,j,k)
                  this% IBMmask(domain)% x(this% IBMmask(domain)% NumOfObjs)% element_index  = eID
                  this% IBMmask(domain)% x(this% IBMmask(domain)% NumOfObjs)% local_Position = (/i,j,k/)
               end if
            end do; end do; end do
         end if
         end associate
      end do
#ifdef _HAS_MPI_
      call castMaskNumOfObjs( this% IBMmask, domain )
      call castMask( this% IBMmask )
#endif
      do domains = 1, MPI_Process% nProcs
         if( this% IBMmask(domains)% NumOfObjs .eq. 0 ) cycle
         do i = 1, this% IBMmask(domains)% NumOfObjs
            call OBB(STLNum)% ChangeRefFrame( this% IBMmask(domains)% x(i)% coords, LOCAL, Point )
            call this% CheckPoint( Point, STLNum, this% IBMmask(domains)% x(i)% NumOfIntersections )
         end do
      end do 
#ifdef _HAS_MPI_
      do domains = 1, MPI_Process% nProcs
         call mpi_allreduce( this% IBMmask(domains)% x(:)% NumOfIntersections, &
                             this% IBMmask(domains)% x(:)% NumOfIntersections, &
                             this% IBMmask(domains)% NumOfObjs, MPI_INT,       &
                             MPI_SUM, MPI_COMM_WORLD, ierr                     )
      end do
#endif
      do i = 1, this% IBMmask(domain)% NumOfObjs
         if( mod(this% IBMmask(domain)% x(i)% NumOfIntersections,2) .ne. 0 ) this% IBMmask(domain)% x(i)% isInsideBody = .true.
      end do

      if( this% HO_IBM) then
         do n = 1, this% IBMmask(domain)% NumOfObjs
            if( this% IBMmask(domain)% x(n)% isInsideBody ) then
               eID = this% IBMmask(domain)% x(n)% element_index
               k   = this% IBMmask(domain)% x(n)% local_Position(1)

               elements(eID)% MaskCorners(k) = .true.
               elements(eID)% STL(1,0,0)     = STLNum

               this% NumOfMaskObjs = this% NumOfMaskObjs + 1
            end if
         end do
#ifdef _HAS_MPI_
         call mpi_allreduce(this% NumOfMaskObjs, NumOfObjs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
#else 
         NumOfObjs = this% NumOfMaskObjs
#endif
         call this% HOmask( elements, faces )
      else
         do n = 1, this% IBMmask(domain)% NumOfObjs
            if( this% IBMmask(domain)% x(n)% isInsideBody ) then
               eID = this% IBMmask(domain)% x(n)% element_index
               i   = this% IBMmask(domain)% x(n)% local_Position(1)
               j   = this% IBMmask(domain)% x(n)% local_Position(2)
               k   = this% IBMmask(domain)% x(n)% local_Position(3)

               elements(eID)% isInsideBody(i,j,k) = .true.
               elements(eID)% STL(i,j,k)          = STLNum

               this% NumOfMaskObjs = this% NumOfMaskObjs + 1
            end if
         end do
#ifdef _HAS_MPI_
         call mpi_allreduce(this% NumOfMaskObjs, NumOfObjs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
#else 
         NumOfObjs = this% NumOfMaskObjs
#endif
      end if
      
      if( NumOfObjs .eq. 0 ) then
         print*, "Mask is made of 0 zero points"
         error stop
      endif

      if( this% plotMask ) call this% plot_Mask( iter, STLNum, NumOfObjs )

      do domains = 1, MPI_Process% nProcs
         if( allocated( this% IBMmask(domains)% x) ) deallocate( this% IBMmask(domains)% x )
      end do

      deallocate(this% IBMmask)

   end subroutine IBM_constructmask

   subroutine IBM_HOmask( this, elements, faces )
      use Meshtypes
      use MPI_Process_Info
      implicit none

      class(IBM_type), intent(inout) :: this
      type(element),   intent(inout) :: elements(:)
      type(face),      intent(inout) :: faces(:)

      real(kind=RP)     :: Point(NDIM)
      integer           :: eID, fID, i, j, STLNum, Nxi, Neta, NumOfPoints
#ifdef _HAS_MPI_
      integer :: ierr
#endif
      this% NumOfMaskObjs = 0
      
      do eID = 1, size(elements)
         associate( e => elements(eID) )
         if(  all(e% MaskCorners((/1,2,6,5/))) .and. .not. ALL(e% MaskCorners((/4,3,7,8/))) ) then
            associate(f => faces(e% faceIDs(EFRONT)))
            if( f% faceType .eq. HMESH_BOUNDARY ) cycle
            Nxi  = f% Nf(1)
            Neta = f% Nf(2)
            if( f% HO_IBM ) then
               f% HO_IBM = .false.
               this% NumOfMaskObjs = this% NumOfMaskObjs - (Nxi+1)*(Neta+1)
               if( f% faceType .eq. HMESH_MPI ) f% fmpi = .false.
            else
               f% HO_IBM = .true.
               allocate(f% stencil(0:Nxi,0:Neta))
               STLNum = e% STL(1,0,0)
               f% HOSIDE = FindHOside( f, e% eID )
               do j = 0, Neta; do i = 0, Nxi
                  f% stencil(i,j)% x   = f% geom% x(:,i,j)
                  this% NumOfMaskObjs = this% NumOfMaskObjs + 1
               end do; end do
               if( f% faceType .eq. HMESH_MPI ) f% fmpi = .true.
            end if
            end associate
         end if
         if( all(e% MaskCorners((/4,3,7,8/))) .and. .not. all( e% MaskCorners((/1,2,6,5/)) ) ) then
            associate(f => faces(e% faceIDs(EBACK)))
            if( f% faceType .eq. HMESH_BOUNDARY ) cycle
            Nxi  = f% Nf(1)
            Neta = f% Nf(2)
            if( f% HO_IBM ) then
               f% HO_IBM = .false.
               this% NumOfMaskObjs = this% NumOfMaskObjs - (Nxi+1)*(Neta+1)
               if( f% faceType .eq. HMESH_MPI ) f% fmpi = .false.
            else
               f% HO_IBM = .true.
               allocate(f% stencil(0:Nxi,0:Neta))
               STLNum = e% STL(1,0,0)
               f% HOSIDE = FindHOside( f, e% eID )
               do j = 0, Neta; do i = 0, Nxi
                  f% stencil(i,j)% x    = f% geom% x(:,i,j)
                  this% NumOfMaskObjs = this% NumOfMaskObjs + 1
               end do; end do
               if( f% faceType .eq. HMESH_MPI ) f% fmpi = .true.
            end if
            end associate
         end if
         if( all(e% MaskCorners((/1,2,3,4/))) .and. .not. all(e% MaskCorners((/2,3,7,6/))) ) then
            associate(f => faces(e% faceIDs(EBOTTOM)))
            if( f% faceType .eq. HMESH_BOUNDARY ) cycle
            Nxi  = f% Nf(1)
            Neta = f% Nf(2)
            if( f% HO_IBM ) then
               f% HO_IBM = .false.
               this% NumOfMaskObjs = this% NumOfMaskObjs - (Nxi+1)*(Neta+1)
               if( f% faceType .eq. HMESH_MPI ) f% fmpi = .false.
            else
               f% HO_IBM = .true.
               allocate(f% stencil(0:Nxi,0:Neta))
               STLNum = e% STL(1,0,0)
               f% HOSIDE = FindHOside( f, e% eID )
               do j = 0, Neta; do i = 0, Nxi
                  f% stencil(i,j)% x   = f% geom% x(:,i,j)
                  this% NumOfMaskObjs = this% NumOfMaskObjs + 1
               end do; end do
               if( f% faceType .eq. HMESH_MPI ) f% fmpi = .true.
            end if
            end associate
         end if
         if( all(e% MaskCorners((/2,3,7,6/))) .and. .not. all(e% MaskCorners((/1,2,3,4/))) ) then
            associate(f => faces(e% faceIDs(ERIGHT)))
            if( f% faceType .eq. HMESH_BOUNDARY ) cycle
            Nxi  = f% Nf(1)
            Neta = f% Nf(2)
            if( f% HO_IBM ) then
               f% HO_IBM = .false.
               this% NumOfMaskObjs = this% NumOfMaskObjs - (Nxi+1)*(Neta+1)
               if( f% faceType .eq. HMESH_MPI ) f% fmpi = .false.
            else
               f% HO_IBM = .true.
               allocate(f% stencil(0:Nxi,0:Neta))
               STLNum = e% STL(1,0,0)
               f% HOSIDE = FindHOside( f, e% eID )
               do j = 0, Neta; do i = 0, Nxi
                  f% stencil(i,j)% x   = f% geom% x(:,i,j)
                  this% NumOfMaskObjs = this% NumOfMaskObjs + 1
               end do; end do
               if( f% faceType .eq. HMESH_MPI ) f% fmpi = .true.
            end if
            end associate
         end if
         if( all(e% MaskCorners((/5,6,7,8/))) .and. .not. all(e% MaskCorners((/1,4,8,5/))) ) then
            associate(f => faces(e% faceIDs(ETOP)))
            if( f% faceType .eq. HMESH_BOUNDARY ) cycle
            Nxi  = f% Nf(1)
            Neta = f% Nf(2)
            if( f% HO_IBM ) then
               f% HO_IBM = .false.
               this% NumOfMaskObjs = this% NumOfMaskObjs - (Nxi+1)*(Neta+1)
               if( f% faceType .eq. HMESH_MPI ) f% fmpi = .false.
            else
               f% HO_IBM = .true.
               allocate(f% stencil(0:Nxi,0:Neta))
               STLNum = e% STL(1,0,0)
               f% HOSIDE = FindHOside( f, e% eID )
               do j = 0, Neta; do i = 0, Nxi
                  f% stencil(i,j)% x   = f% geom% x(:,i,j)
                  this% NumOfMaskObjs = this% NumOfMaskObjs + 1
               end do; end do
               if( f% faceType .eq. HMESH_MPI ) f% fmpi = .true.
            end if
            end associate
         end if
         if( all(e% MaskCorners((/1,4,8,5/))) .and. .not. all(e% MaskCorners((/5,6,7,8/))) ) then
            associate(f => faces(e% faceIDs(ELEFT)))
            if( f% faceType .eq. HMESH_BOUNDARY ) cycle
            Nxi  = f% Nf(1)
            Neta = f% Nf(2)
            if( f% HO_IBM ) then
               f% HO_IBM = .false.
               this% NumOfMaskObjs = this% NumOfMaskObjs - (Nxi+1)*(Neta+1)
               if( f% faceType .eq. HMESH_MPI ) f% fmpi = .false.
            else
               f% HO_IBM = .true.
               allocate(f% stencil(0:Nxi,0:Neta))
               STLNum = e% STL(1,0,0)
               f% HOSIDE = FindHOside( f, e% eID )
               do j = 0, Neta; do i = 0, Nxi
                  f% stencil(i,j)% x   = f% geom% x(:,i,j)
                  this% NumOfMaskObjs = this% NumOfMaskObjs + 1
               end do; end do
               if( f% faceType .eq. HMESH_MPI ) f% fmpi = .true.
            end if
            end associate
         end if
         if(all(e% MaskCorners)) e% HO_IBM = .true.
         end associate
      end do

      do fID = 1, size(faces)
         associate(f => faces(fID))
         if( f% HO_IBM ) then
            elements( f% elementIDs(f% HOSIDE) )% HO_IBM = .true.
         end if
         end associate
      end do

   end subroutine IBM_HOmask

   integer function FindHOside( f, eID ) result( side )

      implicit none

      type(face), intent(inout) :: f
      integer,    intent(in)    :: eID

      ! f% elementSide(1) = 1,2,3,4,5,6 -> side of element on the left where f lies
      ! f% elementIDs(1) = eID -> index of the element on the left of f
      if( f% elementIDs(LEFT) .eq. eID ) then
         side = RIGHT
      elseif( f% elementIDs(RIGHT) .eq. eID ) then
         side = LEFT
      else
         print*, "Can't find face side for HO IBM"
         error stop
      end if

   end function FindHOside
!
!  Mask plot
!  -----------------------------------------------------------
   subroutine IBM_plot_Mask(  this, iter, STLNum, NumOfObjs )
      use MPI_Process_Info
      implicit none
      !-arguments------------------------------------------------------
      class(IBM_type), intent(inout) :: this
      integer,         intent(in)    :: iter, STLNum, NumOfObjs
      !-local-variables------------------------------------------------
      character(len=LINE_LENGTH) :: filename, rank, it, lvl
      integer                    :: i, funit, domains
      logical                    :: add_Point = .false.
#ifdef _HAS_MPI_
      call castMaskPlot(this% IBMmask)
#endif

      if( .not. MPI_Process% isRoot ) return

      write(it,'(I10.10)') iter

      filename = trim(this% STLfilename(STLNum))//'_'//trim(adjustl(it))

      call TecFileHeader( 'IBM/Mask_'//trim(filename), 'Mask Points', NumOfObjs/2+mod(NumOfObjs,2),2,1, funit, 'POINT')
      if( mod(NumOfObjs,2) .ne. 0 )  add_Point  = .true.

      do domains = 1, MPI_Process% nProcs 
         do i = 1, this% IBMmask(domains)% NumOfObjs
            if( this% IBMmask(domains)% x(i)% isInsideBody ) then 
               write(funit,'(3E13.5)') this% IBMmask(domains)% x(i)% coords(IX), this% IBMmask(domains)% x(i)% coords(IY), this% IBMmask(domains)% x(i)% coords(IZ)
               if( add_Point ) then
                  write(funit,'(3E13.5)') this% IBMmask(domains)% x(i)% coords(IX), this% IBMmask(domains)% x(i)% coords(IY), this% IBMmask(domains)% x(i)% coords(IZ)
                  add_Point = .false.
               end if
            end if 
         end do 
      end do   

      close(funit)

   end subroutine IBM_plot_Mask
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------
!  Getting extreme of the domain mesh. These are used to clip the STL file
!  -----------------------------------------------------------------------

   subroutine IBM_GetDomainExtreme( this, elements )
      use MPI_Process_Info
      use ElementConnectivityDefinitions
      implicit none
      !-arguments---------------------------------------------
      class(IBM_type), intent(inout) :: this
      type(element),   intent(in)   :: elements(:)
      !-local-variables---------------------------------------
      real(kind=rp) :: ElemMax, ElemMin
      integer       :: eID, i, axis, j, k
#ifdef _HAS_MPI_
      real(kind=RP) :: localmax, localmin
      integer       :: ierr
#endif
      this% maxCOORDS = -huge(1.0_RP); this% minCOORDS = huge(1.0_RP)
      axis = this% ClipAxis

      if( axis .eq. 0 ) return

      do eID = 1, size(elements)
         associate(e => elements(eID))
         if( this% HO_IBM ) then
            do k = 1, NODES_PER_ELEMENT
               this% maxCOORDS = max(this% maxCOORDS,e% SurfInfo% corners(axis,k)); this% minCOORDS = min(this% minCOORDS,e% SurfInfo% corners(axis,k))
            end do
         else
            do k = 0, elements(eID)% Nxyz(3)   ; do j = 0, elements(eID)% Nxyz(2) ; do i = 0, elements(eID)% Nxyz(1)
               this% maxCOORDS = max(this% maxCOORDS,elements(eID)% geom% x(axis,i,j,k)); this% minCOORDS = min(this% minCOORDS,elements(eID)% geom% x(axis,i,j,k))
            end do; end do; end do
         end if
         end associate
      end do
      this% maxCOORDS = this% maxCOORDS + 1.0e-8
      this% minCOORDS = this% minCOORDS - 1.0e-8
#ifdef _HAS_MPI_
      localmax = this% maxCOORDS; localmin = this% minCOORDS
      call mpi_allreduce(localmax, this% maxCOORDS, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, ierr)
      call mpi_allreduce(localmin, this% minCOORDS, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, ierr)
#endif
   end subroutine IBM_GetDomainExtreme

   subroutine PlotSurface( IBM, faces, STLNum, iter )
      use MPI_Process_Info
      implicit none

      type(IBM_type), intent(in) :: IBM
      type(face),     intent(in) :: faces(:)
      integer,        intent(in) :: STLNum, iter

      character(len=LINE_LENGTH) :: filename, rank
      real(kind=RP)              :: x(NDIM)
      integer                    :: fID, domains, NumOfObjs, funit, i, j
      logical                    :: add_Point = .false.

      if( .not. IBM% HO_IBM ) return

      domains   = MPI_Process% nProcs
      NumOfObjs = 0
      do fID = 1, size(faces)
         associate( f => faces(fID) )
         if( f% HO_IBM ) NumOfObjs = NumOfObjs + (f% Nf(1)+1)*(f% Nf(2)+1)
         end associate
      end do

      if( NumOfObjs .eq. 0 ) return
      
      write(rank,*) MPI_Process% rank

      if( IBM% lvl .gt. 0 ) then
         write(filename,'(A,A,I1)') trim(IBM% STLfilename(STLNum)),'_MGlevel',IBM% lvl
      else
         write(filename,'(A,A,A)') trim(IBM% STLfilename(STLNum)),'_',trim(adjustl(rank))
      end if

      call TecFileHeader( 'IBM/Surface_'//trim(filename), 'Mask Points', NumOfObjs/2+mod(NumOfObjs,2),2,1, funit, 'POINT')

      if( mod(NumOfObjs,2) .ne. 0 )  add_Point  = .true.

      do fID = 1, size(faces)
         associate( f=> faces(fID) )
         if( f% HO_IBM ) then
            do i = 0, f% Nf(1); do j = 0, f% Nf(2)
               x = f% stencil(i,j)% x - f% stencil(i,j)% dist * f% stencil(i,j)% normal
               write(funit,'(3E13.5)')  x(1), x(2), x(3)
               if( add_Point ) then
                  write(funit,'(3E13.5)') x(1), x(2), x(3)
                  add_Point = .false.
               end if
            end do; end do
         end if
         end associate
      end do

   end subroutine PlotSurface

   subroutine PlotStencil( IBM, STLNum, iter )
      use MPI_Process_Info
      implicit none

      type(IBM_type), intent(in) :: IBM
      integer,        intent(in) :: STLNum, iter

      character(len=LINE_LENGTH) :: filename, rank
      real(kind=RP) :: x(NDIM)
      integer       :: i, j, k, fID, NumOfObjs, Nxi, Neta, funit, domains, domain
      logical       :: add_Point = .false.

      if( .not. IBM% HO_IBM  .and. .not. MPI_Process% isRoot ) return!.or. .not. MPI_Process% isRoot ) return

      domains   = MPI_Process% nProcs
      NumOfObjs = sum(IBMStencilPoints(1:domains)% NumOfObjs)
      
      if( IBM% lvl .gt. 0 ) then
         write(filename,'(A,A,I1,A,I10.10)') trim(IBM% STLfilename(STLNum)),'_MGlevel',IBM% lvl,'_',iter
      else
         write(filename,'(A)') trim(IBM% STLfilename(STLNum))
      end if

      call TecFileHeader( 'IBM/StencilPoints_'//trim(filename), 'Mask Points', NumOfObjs/2+mod(NumOfObjs,2),2,1, funit, 'POINT')

      if( mod(NumOfObjs,2) .ne. 0 )  add_Point  = .true.

      do domain = 1, domains
         do fID = 1, IBM_HO_faces(domain)% NumOfFaces
            associate( f => IBM_HO_faces(domain)% faces(fID) )
               Nxi  = f% Nf(1)
               Neta = f% Nf(2)
               do i = 0, Nxi; do j = 0, Neta
                  do k = 0, f% stencil(i,j)% N
                     x = f% stencil(i,j)% x_s(:,k)
                     write(funit,'(3E13.5)')  x(1), x(2), x(3)
                     if( add_Point ) then
                        write(funit,'(3E13.5)') x(1), x(2), x(3)
                        add_Point = .false.
                     end if
                  end do
               end do; end do
            end associate
         end do
      end do

      close(funit)

   end subroutine PlotStencil
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------
!  building the immersed boundary
!  -------------------------------------------------
   subroutine IBM_build( this, elements, faces, no_of_DoFs, isChild, movingSTL, iter )
      use MPI_Process_Info
      use PhysicsStorage
      implicit none
      !-arguments-----------------------------------------------------------------
      class(IBM_type),           intent(inout) :: this
      type(element),             intent(inout) :: elements(:)
      type(face),                intent(inout) :: faces(:)
      integer,                   intent(in)    :: no_of_DoFs
      logical,                   intent(in)    :: isChild
      integer,         optional, intent(in)    :: movingSTL, iter
      !-local-variables-----------------------------------------------------------
      integer :: MaskPoints, STLNum
#ifdef _HAS_MPI_
      integer :: localVal, ierr
#endif
      do STLNum = 1, this% NumOfSTL
         if( .not. present(movingSTL) )then
            call this% GetMask( elements, faces, no_of_DoFs, STLNum, iter )
            if( this% ComputeDistance .and. .not. this% HO_IBM ) call this% ComputeIBMWallDistance( elements=elements, STLNum=STLNum )
         else
            if( STLNum .eq. movingSTL ) then
               call this% GetMask( elements, faces, no_of_DoFs, STLNum, iter )
               if( this% ComputeDistance .and. .not. this% HO_IBM ) call this% ComputeIBMWallDistance( elements=elements, STLNum=STLNum )
            end if
         end if
      end do

      if( .not. present(movingSTL) )then
         allocate( this% penalization(size(elements)) )
         this% penalization = this% eta
      end if

      if( this% ComputeBandRegion ) then
         do STLNum = 1, this% NumOfSTL
            if( .not. present(movingSTL) ) then
               call this% constructBandRegion( elements, no_of_DoFs, STLNum, NCONS )
            else
               if( STLNum .eq. movingSTL ) then
                  call this% constructBandRegion( elements, no_of_DoFs, STLNum, NCONS )
               end if
            end if
         end do
      end if

      if( this% Wallfunction ) then
         call this% GetForcingPointsGeom( elements )
         do STLNum = 1, this% NumOfSTL
            call this% GetImagePoint_nearest( STLNum )
         end do
      end if

   end subroutine IBM_build

   subroutine IBM_buildHOfaces( this, elements, faces )
      use PhysicsStorage
      use MPI_Process_Info
      implicit none

      class(IBM_type), intent(inout) :: this
      type(element),   intent(inout) :: elements(:)
      type(face),      intent(inout) :: faces(:)

      integer :: STLNum
      
      do STLNum = 1, this% NumOfSTL
         call this% ComputeIBMWallDistance( elements=elements, STLNum = STLNum, faces=faces )
      end do
      
      call this% GetStencil( elements, faces )
      call Set_IBM_HO_faces( faces, NCONS )
#ifdef _HAS_MPI_
      call MPIProcedures_IBM_HO_faces( NCONS )
#endif
      do STLNum = 1, this% NumOfSTL
         call PlotSurface( this, faces, STLNum, 0 )
         call PlotStencil( this, STLNum, 0 )
      end do
      
      call IBM_HO_findElements( elements )
      
   end subroutine IBM_buildHOfaces
!
!  Immersed boundary description
!  -----------------------------
   subroutine IBM_Describe( this )
      use Headers
      use mainKeywordsModule
      use MPI_Process_Info
      implicit none
      !-arguments-----------------------------
      class(IBM_type),  intent(inout) :: this

      if ( .not. MPI_Process % isRoot ) return
      write(STD_OUT,'(/)')
      call Section_Header("IBM parameters")
      write(STD_OUT,'(/)')

      call SubSection_Header('IBM info')
      if( .not. this% HO_IBM ) then
         write(STD_OUT,'(30X,A,A35,A20)') "->" ,"Type: ","volume penalization"
      else
         write(STD_OUT,'(30X,A,A35,A10)') "->" ,"Type: ","high order"
      end if
         write(STD_OUT,'(30X,A,A35,L10)') "->" , "Semi implicit treatment: ", this% semiImplicit
      if( .not. this% TimePenal ) then
         write(STD_OUT,'(30X,A,A35,ES14.2)') "->" , "Penalization term: " , this% eta
      else
         write(STD_OUT,'(30X,A,A35,A10)') "->" , "Penalization term: ", " Dt"
      end if

      write(STD_OUT,'(30X,A,A35,I10)') "->" , "Minimum number of objects: ", this% KDtree_Min_n_of_Objs
      write(STD_OUT,'(30X,A,A35,I10)') "->" , "Number of interpolation points: ", this% NumOfInterPoints

   end subroutine IBM_Describe
!
!  Destroy immersed boundary
!  -------------------------
   subroutine IBM_Destruct( this, isChild )
      use MPI_Process_Info
      implicit none
      !-arguments------------------------------
      class(IBM_type), intent(inout) :: this
      logical,         intent(in)    :: isChild
      !----------------------------------------
      integer :: STLNum, i, j, domains

      do STLNum = 1, this% NumOfSTL
         call this% root(STLNum)% destruct( isChild )
         if( this% ComputeDistance ) call this% rootDistance(STLNum)% destruct( isChild )
      end do

      if( this% ComputeBandRegion ) deallocate( this% BandRegion )

      if( this% Wallfunction ) then
         deallocate( this% ImagePoints )
      end if

      deallocate( this% penalization, &
                  this% stl,          &
                  this% Integral,     &
                  this% STLfilename   )

   end subroutine IBM_Destruct
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------------------------
!  Creating the .omesh file for p-adaptation close to the STL
!  ------------------------------------------------------------
   subroutine IBM_SetPolynomialOrder( this, elements, corners )
      implicit none

      !-arguments----------------------------------------------
      class(IBM_type),         intent(inout) :: this
      type(element),           intent(inout) :: elements(:)
      real(kind=RP), optional, intent(in)    :: corners(:,:)
      !-local-variables----------------------------------------
      real(kind=RP) :: corner(NDIM)
      integer       :: STLNum, eID, k

      if( present(corners) ) then
         do eID = 1, size(elements)
            associate(e => elements(eID))
            do k = 1, BOXVERTICES
               if( isInsideBox( e% SurfInfo% corners(:,k), corners ) ) then
                  e% Nxyz(1) = this% Nx
                  e% Nxyz(2) = this% Ny
                  e% Nxyz(3) = this% Nz
                  exit
               end if
            end do
            end associate
         end do
      else
         do eID = 1, size(elements)
            associate(e => elements(eID))
            do STLNum = 1, this% NumOfSTL
               loop: do k = 1, BOXVERTICES
                  call OBB(STLNum)% ChangeRefFrame( e% SurfInfo% corners(:,k), LOCAL, corner )
                  if( isInsideBox( corner, this% BandRegionCoeff*OBB(STLNum)% LocVertices ) ) then
                     if( this% Nx .ne. 0 ) e% Nxyz(1) = this% Nx
                     if( this% Ny .ne. 0 ) e% Nxyz(2) = this% Ny
                     if( this% Nz .ne. 0 ) e% Nxyz(3) = this% Nz
                     exit loop
                  end if
               end do loop
            end do
            end associate
         end do
      end if

   end subroutine IBM_SetPolynomialOrder
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------------------
!  Band region for the computation of the forces and turbulent quantities
!  -----------------------------------------------------------------------

   subroutine IBM_constructBandRegion( this, elements, no_of_DoFs, STLNum, nEqn )
      use MPI_Process_Info
      use PhysicsStorage
      implicit none
      !-arguments------------------------------------------------------
      class(IBM_type), intent(inout) :: this
      type(element),   intent(inout) :: elements(:)
      integer,         intent(in)    :: no_of_DoFs, STLNum, nEqn
      !-local-variables------------------------------------------------
      integer :: eID, i, j, k, domains, domain, NumOfObjs

      allocate(this% BandRegion(STLNum)% IBMmask(MPI_Process% nProcs))

      domain = MPI_Process% rank + 1

      allocate(this% BandRegion(STLNum)% IBMmask(domain)% x(no_of_DoFs))

      this% BandRegion(STLNum)% IBMmask(domain)% NumOfObjs = 0

      do eID = 1, size(elements)
         associate(e => elements(eID))
         do i = 0, e% Nxyz(1); do j = 0, e% Nxyz(2); do k = 0, e% Nxyz(3)
            if( e% isInsideBody(i,j,k) ) cycle
            if( OBB(STLNum)% isPointInside( coords = e% geom% x(:,i,j,k), coeff=this% BandRegionCoeff ) ) then
               this% BandRegion(STLNum)% IBMmask(domain)% NumOfObjs = this% BandRegion(STLNum)% IBMmask(domain)% NumOfObjs + 1
               this% BandRegion(STLNum)% IBMmask(domain)% x(this% BandRegion(STLNum)% IBMmask(domain)% NumOfObjs)% coords         = e% geom% x(:,i,j,k)
               this% BandRegion(STLNum)% IBMmask(domain)% x(this% BandRegion(STLNum)% IBMmask(domain)% NumOfObjs)% element_index  = eID
               this% BandRegion(STLNum)% IBMmask(domain)% x(this% BandRegion(STLNum)% IBMmask(domain)% NumOfObjs)% local_Position = (/i,j,k/)
            end if
         end do; end do; end do
         end associate
      end do
#ifdef _HAS_MPI_
      call castMaskNumOfObjs( this% BandRegion(STLNum)% IBMmask, domain )
      call castMask( this% BandRegion(STLNum)% IBMmask )
#endif
      if( this% plotBandPoints .and. MPI_Process% isRoot ) call this% BandRegion(STLNum)% plot( 'BandPoints' )

   end subroutine IBM_constructBandRegion

   subroutine IBM_minDistance( this, Point, STLNum, indeces, domain )
      use MPI_Process_Info
      implicit none
      !-arguments-----------------------------------------------
      class(IBM_type), intent(inout) :: this
      real(kind=RP),   intent(in)    :: Point(NDIM)
      integer,         intent(in)    :: STLNum
      integer,         intent(inout) :: indeces(:), &
                                        domain(:)
      !-local-variables-----------------------------------------
      real(kind=RP) :: dist, dist_, lastDist
      integer       :: k, domains, i

      lastDist = -huge(1.0_RP)
      do k = 1, this% NumOfInterPoints
         dist = huge(1.0_RP)
         do domains = 1, MPI_Process% nProcs
            do i = 1, this% BandRegion(STLNum)% IBMmask(domains)% NumOfObjs
               dist_ = norm2(Point - this% BandRegion(STLNum)% IBMmask(domains)% x(i)% coords)
               if( dist_ .eq. 0.0_RP ) cycle
               if( dist_ .lt. dist .and. dist_ .gt. lastDist ) then
                  indeces(k) = i
                  domain(k)  = domains
                  dist       = dist_
               end if
            end do
         end do
         lastDist = dist
      end do

   end subroutine IBM_minDistance

   subroutine IBM_GetBandRegionStates( this, elements, nEqn )
      use PhysicsStorage
      implicit none
      !-arguments-----------------------------------------------
      class(IBM_type), intent(inout) :: this
      type(element),   intent(in)    :: elements(:)
      integer,         intent(in)    :: nEqn
      !-local-variables-----------------------------------------
      integer :: STLNum

      do STLNum = 1, this% NumOfSTL
         call BandPointsState( this% BandRegion(STLNum), elements, nEqn )
      end do

   end subroutine IBM_GetBandRegionStates
!
!  Band region points' state is stored
!  -----------------------------------
   subroutine BandPointsState( this, elements, nEqn )
      use PhysicsStorage
      use MPI_Process_Info
      implicit none
      !-arguments------------------------------------------------
      type(bandRegion_t), intent(inout) :: this
      type(element),      intent(in)    :: elements(:)
      integer,            intent(in)    :: nEqn
      !-local-variables------------------------------------------
      integer                    :: n, i, j, k, eID, domains, domain

      domain = MPI_Process% rank + 1

      do n = 1, this% IBMmask(domain)% NumOfObjs
         i   = this% IBMmask(domain)% x(n)% local_Position(1)
         j   = this% IBMmask(domain)% x(n)% local_Position(2)
         k   = this% IBMmask(domain)% x(n)% local_Position(3)
         eID = this% IBMmask(domain)% x(n)% element_index

         this% IBMmask(domain)% x(n)% Q   = elements(eID)% storage% Q(:,i,j,k)
         this% IBMmask(domain)% x(n)% U_x = elements(eID)% storage% U_x(:,i,j,k)
         this% IBMmask(domain)% x(n)% U_y = elements(eID)% storage% U_y(:,i,j,k)
         this% IBMmask(domain)% x(n)% U_z = elements(eID)% storage% U_z(:,i,j,k)
      end do
#ifdef _HAS_MPI_
      call castStateBandRegion( this% IBMmask, nEqn )
      call castGradientsBandRegion( this% IBMmask, nEqn )
#endif
   end subroutine BandPointsState
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------------------------
!  MPI routines
!  ------------------------------------------------------------
   subroutine IBM_MPI_sendSTLpartitions( this, STLNum, axis )
      use MPI_Process_Info
      implicit none
      !-arguments---------------------------
      class(IBM_type), intent(inout) :: this
      integer,         intent(in)    :: STLNum, axis

      if( MPI_Process% isRoot ) then
         call sendSTL2Partitions( this% stl(STLNum), STLNum, axis )
      end if

      if( MPI_Process% doMPIAction ) then
         call receiveSTLpartitions( this% stl(STLNum), STLNum )
      end if

      this% stl(STLNum)% construct = .true.

   end subroutine IBM_MPI_sendSTLpartitions

   subroutine IBM_MPI_OBB( this, STLNum )
      use MPI_Process_Info
      implicit none

      class(IBM_type), intent(inout) :: this
      integer,         intent(in)    :: STLNum

      if( MPI_Process% doMPIAction ) then
         call SendOBB( STLNum )
         call RecvOBB( STLNum )
         call SendAxis( STLNum )
         call RecvAxis( STLNum )
      end if

   end subroutine IBM_MPI_OBB
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------------------------
!  Moving bodies
!  ------------------------------------------------------------
   subroutine IBM_MoveBody( this, elements, faces, no_of_DoFs, isChild, t, iter, autosave )
      use MPI_Process_Info
      implicit none
      !-arguments-----------------------------------------------------------------
      class(IBM_type),           intent(inout) :: this
      type(element),             intent(inout) :: elements(:)
      type(face),                intent(inout) :: faces(:)
      integer,                   intent(in)    :: no_of_DoFs
      logical,                   intent(in)    :: isChild
      real(kind=RP),             intent(in)    :: t
      integer,         optional, intent(in)    :: iter
      logical,         optional, intent(in)    :: autosave
      !-local-variables-----------------------------------------------------------
      integer :: STLNum

      do STLNum = 1, this% NumOfSTL
         if( this% stl(STLNum)% move ) call this% CleanMask( elements, STLNum )
      end do

      do STLNum = 1, this% NumOfSTL
         if( this% stl(STLNum)% move ) then
            call this% root(STLNum)% Destruct( isChild )
            if( this% ComputeDistance ) call this% rootDistance(STLNum)% destruct( isChild )
            call OBB(STLNum)% ChangeObjsRefFrame( this% stl(STLNum)% ObjectsList, GLOBAL )
            select case( this% stl(STLNum)% motionType )
            case( ROTATION )
               call this% stl(STLNum)% getRotationaMatrix( t )
               call this% stl(STLNum)% rotate( .false. )
            case( LINEAR )
               call this% stl(STLNum)% getDisplacement( t )
               call this% stl(STLNum)% translate( .false. )
            end select
            call this% stl(STLNum)% updateNormals()
            this% stl(STLNum)% show = .false.
            this% plotMask          = .false.
            if( present(autosave) .and. autosave ) this% plotMask = .true.
            if( present(autosave) .and. autosave ) call plotSTL( this% stl(STLNum), iter )
            if( MPI_Process% isRoot ) then
               call OBB(STLNum)% construct( this% stl(STLNum), this% plotOBB, this% AAB, .true. )
               OBB(STLNum)% maxAxis = OBB(STLNum)% GetMaxAxis()
               OBB(STLNum)% minAxis = OBB(STLNum)% GetMinAxis( OBB(STLNum)% maxAxis, this% ClipAxis )
            end if
#ifdef _HAS_MPI_
            call this% MPI_OBB(STLNum)
#endif
            call obb(stlnum)% plot(mpi_process% rank)
            call OBB(STLNum)% ChangeObjsRefFrame( this% stl(STLNum)% ObjectsList, LOCAL )
            this% plotKDtree = .false.
            call this% constructSTL_KDtree( STLNum )
         end if
      end do

      do STLNum = 1, this% NumOfSTL
         if( this% stl(STLNum)% move ) call this% build( elements, faces, no_of_DoFs, isChild, STLNum, iter )
      end do

   end subroutine IBM_MoveBody
!
!  Mask is cleaned
!  ---------------
   subroutine IBM_CleanMask( this, elements, STLNum )

      implicit none
      !-arguments------------------------------------------------
      class(IBM_type), intent(inout) :: this
      type(element),   intent(inout) :: elements(:)
      integer,         intent(in)    :: STLNum
      !-local-variables------------------------------------------
      integer :: eID, i, j, k
!$omp parallel
!$omp do schedule(runtime) private(i,j,k)
      do eID = 1, size(elements)
         associate(e => elements(eID))
         do i = 0, e% Nxyz(1); do j = 0, e% Nxyz(2); do k = 0, e% Nxyz(3)
            if( e% STL(i,j,k) .eq. STLNum .and. e% isInsideBody(i,j,k) ) then
               e% STL(i,j,k)          = 0
               e% isInsideBody(i,j,k) = .false.
               this% NumOfMaskObjs    = this% NumOfMaskObjs - 1
            end if
         end do; end do; end do
         end associate
      end do
!$omp end do
!$omp end parallel
   end subroutine IBM_CleanMask
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------------------------------------------
!  Geometrical quantities for the turbulence wall model
!  ----------------------------------------------------
   subroutine IBM_GetForcingPointsGeom( this, elements )
      use MPI_Process_Info
      use MeshTypes
      implicit none
      !-arguments----------------------------------------
      class(IBM_type), intent(inout) :: this
      type(element),   intent(inout) :: elements(:)
      !-local-variables----------------------------------
      real(kind=RP) :: h, d, d_min, Dist
      integer       :: eID, i, j, k, n, Total
#ifdef _HAS_MPI_
      integer       :: ierr
#endif
      this% NumOfForcingPoints = 0
      d = huge(1.0_RP)

      do eID = 1, size(elements)
         h = (elements(eID)% geom% Volume)**(1.0_RP/3.0_RP)
         d = min(d,h)
      end do
      d_min = 1.5_RP*sqrt(3.0_RP)*d
#ifdef _HAS_MPI_
      call mpi_allreduce(d_min, this% IP_Distance, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, ierr)
      Dist = this% IP_Distance
#else
      this% IP_Distance = d_min
      Dist = this% IP_Distance
#endif
      do eID = 1, size(elements)
         associate(e => elements(eID))
         do k = 0, e% Nxyz(3); do j = 0, e% Nxyz(2); do i = 0, e% Nxyz(1)
            if( e% geom% dWall(i,j,k) .gt. Dist ) cycle
            e% isForcingPoint(i,j,k) = .true.
            this% NumOfForcingPoints = this% NumOfForcingPoints + 1
         end do; end do; end do
         end associate
      end do
#ifdef _HAS_MPI_
      call mpi_allreduce(this% NumOfForcingPoints, Total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
      Total = this% NumOfForcingPoints
#endif
      if( Total .eq. 0 .and. MPI_Process% isRoot ) then
         print*, "No forcing points found, change y+ target"
         error stop
      endif

   end subroutine IBM_GetForcingPointsGeom

   subroutine IBM_constructImagePoint( this, elements, STLNum )
      use MPI_Process_Info
      implicit none

      class(IBM_type), intent(inout) :: this
      type(element),   intent(inout) :: elements(:)
      integer,         intent(in)    :: STLNum

      real(kind=RP) :: dist
      integer       :: domain, eID, i, j, k

      domain = MPI_Process% rank + 1

      this% ImagePoints(STLNum)% IBMmask(domain)% NumOfObjs = 0

      do eID = 1, size(elements)
         associate(e => elements(eID))
         do k = 0, e% Nxyz(3); do j = 0, e% Nxyz(2); do i = 0, e% Nxyz(1)
            if( e% isForcingPoint(i,j,k) .and. e% STL(i,j,k) .eq. STLNum ) then
               this% ImagePoints(STLNum)% IBMmask(domain)% NumOfObjs = this% ImagePoints(STLNum)% IBMmask(domain)% NumOfObjs + 1
               e% IP_index = this% ImagePoints(STLNum)% IBMmask(domain)% NumOfObjs
               dist        = this% IP_distance - e% geom% dWall(i,j,k)

               this% ImagePoints(STLNum)% IBMmask(domain)% x(this% ImagePoints(STLNum)% IBMmask(domain)% NumOfObjs)% coords         = e% geom% x(:,i,j,k) + dist * e% geom% normal(:,i,j,k)
               this% ImagePoints(STLNum)% IBMmask(domain)% x(this% ImagePoints(STLNum)% IBMmask(domain)% NumOfObjs)% element_index  = eID
               this% ImagePoints(STLNum)% IBMmask(domain)% x(this% ImagePoints(STLNum)% IBMmask(domain)% NumOfObjs)% local_Position = (/i,j,k/)
            end if
         end do; end do; end do
         end associate
      end do
#ifdef _HAS_MPI_
      call castMaskNumOfObjs( this% ImagePoints(STLNum)% IBMmask, domain )
      call castMask(this% ImagePoints(STLNum)% IBMmask )
#endif
      if(MPI_Process% isRoot) call this% ImagePoints(STLNum)% plot( 'ImagePoints' )

   end subroutine IBM_constructImagePoint

   subroutine Plot_Forcing_Imagepoints( IBM, elements )
      use MPI_Process_Info
      implicit none
      !-arguments----------------------------------------
      type(IBM_type), intent(in) :: IBM
      type(element),  intent(in) :: elements(:)
      !-local-variables----------------------------------
      real(kind=RP)              :: ImagePoint(NDIM), &
                                    Dist
      character(len=LINE_LENGTH) :: filenameFP,       &
                                    filenameIP,       &
                                    lvl, rank
      integer                    :: eID, i, j, k,     &
                                    funit
      logical                    :: add_Point = .false.

      if( .not. MPI_Process% isRoot ) return

      call TecFileHeader( 'IBM/'//trim(filenameFP), 'Forcing Points', IBM% NumOfForcingPoints/2+mod(IBM% NumOfForcingPoints,2),2,1, funit, 'POINT' )

      if( mod(IBM% NumOfForcingPoints,2) .ne. 0 ) add_Point  = .true.

      do eID = 1, size(elements)
         associate ( e => elements(eID) )
         do k = 0, e% Nxyz(3); do j = 0, e% Nxyz(2) ; do i = 0, e% Nxyz(1)
            if( e% isForcingPoint(i,j,k) ) then
               write(funit,'(3E13.5)')  e% geom% x(1,i,j,k), e% geom% x(2,i,j,k), e% geom% x(3,i,j,k)
               if( add_Point ) then
                  write(funit,'(3E13.5)') e% geom% x(1,i,j,k), e% geom% x(2,i,j,k), e% geom% x(3,i,j,k)
                  add_Point = .false.
               end if
            end if
         end do; end do; end do
         end associate
      end do

      close(funit)

      if( IBM% NumOfForcingPoints .eq. 0 ) return

      if( MPI_Process% nProcs .gt. 1 ) then
         write(rank,*) MPI_Process% rank
         if( IBM% lvl .gt. 0 ) then
            write(lvl,*)  IBM% lvl
            filenameFP = 'ForcingPoints_'//trim(IBM% filename)//'_MGlevel'//trim(adjustl(lvl))//'_Process'//trim(adjustl(rank))
         else
            filenameFP = 'ForcingPoints_'//trim(IBM% filename)//'_Process'//trim(adjustl(rank))
         end if
      else
         if( IBM% lvl .gt. 0 ) then
            write(lvl,*)  IBM% lvl
            filenameFP = 'ForcingPoints_'//trim(IBM% filename)//'_MGlevel'//trim(adjustl(lvl))
         else
            filenameFP = 'ForcingPoints_'//trim(IBM% filename)
         end if
      endif

      call TecFileHeader( 'IBM/'//trim(filenameFP), 'Forcing Points', IBM% NumOfForcingPoints/2+mod(IBM% NumOfForcingPoints,2),2,1, funit, 'POINT' )

      if( mod(IBM% NumOfForcingPoints,2) .ne. 0 ) add_Point  = .true.

      do eID = 1, size(elements)
         associate ( e => elements(eID) )
         do k = 0, e% Nxyz(3); do j = 0, e% Nxyz(2) ; do i = 0, e% Nxyz(1)
            if( e% isForcingPoint(i,j,k) ) then
               write(funit,'(3E13.5)')  e% geom% x(1,i,j,k), e% geom% x(2,i,j,k), e% geom% x(3,i,j,k)
               if( add_Point ) then
                  write(funit,'(3E13.5)') e% geom% x(1,i,j,k), e% geom% x(2,i,j,k), e% geom% x(3,i,j,k)
                  add_Point = .false.
               end if
            end if
         end do; end do; end do
         end associate
      end do

      close(funit)

   end subroutine Plot_Forcing_Imagepoints

!
!  DoFs closer to each image point
!  -------------------------------
   subroutine IBM_GetImagePoint_nearest( this, STLNum )
      use MappedGeometryClass
      use PhysicsStorage
      use MPI_Process_Info
      use ElementClass
      implicit none
      !-arguments-------------------------------------------------------------------------------------
      class(IBM_type), intent(inout) :: this
      integer,         intent(in)    :: STLNum
      !-local-variables-------------------------------------------------------------------------------
      type(point_type) :: x(this% NumOfInterPoints)
      real(kind=RP)    :: Point(NDIM)
      integer          :: domains, i, k, index, domain

      do domains = 1, MPI_Process% nProcs
         do i = 1, this% ImagePoints(STLNum)% IBMmask(domains)% NumOfObjs
            Point = this% ImagePoints(STLNum)% IBMmask(domains)% x(i)% coords
            call this% minDistance( Point, STLNum, this% ImagePoints(STLNum)% IBMmask(domains)% x(i)% indeces, &
                                                   this% ImagePoints(STLNum)% IBMmask(domains)% x(i)% domains  )

            do k = 1, this% NumOfInterPoints
               domain = this% ImagePoints(STLNum)% IBMmask(domains)% x(i)% domains(k)
               index  = this% ImagePoints(STLNum)% IBMmask(domains)% x(i)% indeces(k)
               x(k)% coords = this% bandRegion(STLNum)% IBMmask(domain)% x(index)% coords
            end do

            call GetMatrixInterpolationSystem( this% ImagePoints(STLNum)% IBMmask(domains)% x(i)% coords, &
                                               x,                                                         &
                                               this% ImagePoints(STLNum)% IBMmask(domains)% x(i)% invPhi, &
                                               this% ImagePoints(STLNum)% IBMmask(domains)% x(i)% b,      &
                                               this% InterpolationType                                    )
         end do
      end do

   end subroutine IBM_GetImagePoint_nearest


   subroutine IBM_GetImagePointsStates( this, nEqn )
      use PhysicsStorage
      implicit none
      !-arguments-----------------------------------------------
      class(IBM_type), intent(inout) :: this
      integer,         intent(in)    :: nEqn
      !-local-variables-----------------------------------------
      integer :: STLNum

      do STLNum = 1, this% NumOfSTL
         call ImagePointsState( this% ImagePoints(STLNum), this% bandRegion(STLNum), nEqn, this% NumOfInterPoints, this% InterpolationType )
      end do

   end subroutine IBM_GetImagePointsStates

    subroutine ImagePointsState( this, bandRegion, nEqn, NumOfInterPoints, InterpolationType )
      use PhysicsStorage
      use MPI_Process_Info
      implicit none
      !-arguments------------------------------------------------
      type(bandRegion_t), intent(inout) :: this
      type(bandRegion_t), intent(in)    :: bandRegion
      integer,            intent(in)    :: nEqn, NumOfInterPoints, &
                                           InterpolationType

      real(kind=RP) ::  Q(nEqn,NumOfInterPoints)
      integer       :: domain, i, k, domains, index

      domain = MPI_Process% rank + 1

      do i = 1, this% IBMmask(domain)% NumOfObjs
         associate( Point => this% IBMmask(domain)% x(i) )
         do k = 1, NumOfInterPoints
            domains = Point% domains(k)
            index   = Point% indeces(k)
            Q(:,k)  = bandRegion% IBMmask(domains)% x(index)% Q
         end do
         do k = 1, NCONS
            Point% Q(k) = GetInterpolatedValue( Q(k,:), Point% invPhi, Point% b, InterpolationType )
         end do
         end associate
      end do

    end subroutine ImagePointsState
!
!  Distance from the wall needed
!  -----------------------------
   subroutine IBM_ComputeIBMWallDistance( this, elements, STLNum, movingSTL, faces )
      use MPI_Process_Info
      implicit none
      !-arguments--------------------------------------
      class(IBM_type),      intent(inout) :: this
      type(element),        intent(inout) :: elements(:)
      integer,              intent(in)    :: STLNum
      integer, optional ,   intent(in)    :: movingSTL
      type(face), optional, intent(inout) :: faces(:)
      !-local-variables--------------------------------
      real(kind=RP) :: Point(NDIM), normal(NDIM), dist
      integer       :: domain, domains, i, j, k, n, fID, eID
#ifdef _HAS_MPI_
      integer       :: ierr
#endif
      domain = MPI_Process% rank + 1
      allocate(this% IBMmask(MPI_Process% nProcs))
      allocate(this% IBMmask(domain)% x(this% NumOfMaskObjs))

      this% IBMmask(domain)% NumOfObjs = 0

      if( this% HO_IBM ) then
         do fID = 1, size(faces)
            associate(f => faces(fID))
            if( f% HO_IBM ) then
               do i = 0, f% Nf(1); do j = 0, f% Nf(2)
                  this% IBMmask(domain)% NumOfObjs = this% IBMmask(domain)% NumOfObjs + 1
                  this% IBMmask(domain)% x(this% IBMmask(domain)% NumOfObjs)% coords         = f% geom% x(:,i,j)
                  this% IBMmask(domain)% x(this% IBMmask(domain)% NumOfObjs)% element_index  = fID
                  this% IBMmask(domain)% x(this% IBMmask(domain)% NumOfObjs)% local_Position = (/i,j,0/)
               end do; end do
            end if
            end associate
         end do
      else
         do eID = 1, size(elements)
            associate(e => elements(eID))
            do i = 0, e% Nxyz(1); do j = 0, e% Nxyz(2); do k = 0, e% Nxyz(3)
               if( e% isInsideBody(i,j,k) ) cycle
               if( OBB(STLNum)% isPointInside( coords = e% geom% x(:,i,j,k), coeff=this% BandRegionCoeff ) ) then
                  e% STL(i,j,k) = STLNum
                  this% IBMmask(domain)% NumOfObjs = this% IBMmask(domain)% NumOfObjs + 1
                  this% IBMmask(domain)% x(this% IBMmask(domain)% NumOfObjs)% coords         = e% geom% x(:,i,j,k)
                  this% IBMmask(domain)% x(this% IBMmask(domain)% NumOfObjs)% element_index  = eID
                  this% IBMmask(domain)% x(this% IBMmask(domain)% NumOfObjs)% local_Position = (/i,j,k/)
               end if
            end do; end do; end do
            end associate
         end do
      end if
#ifdef _HAS_MPI_
      call castMaskNumOfObjs( this% IBMmask, domain )
      call castMask( this% IBMmask )
#endif
      do domains = 1, MPI_Process% nProcs
         do n = 1, this% IBMmask(domains)% NumOfObjs
            this% IBMmask(domains)% x(n)% Dist = HUGE(1.0_RP)
            call OBB(STLNum)% ChangeRefFrame(this% IBMmask(domains)% x(n)% coords,LOCAL,Point)
            call MinimumDistance( Point, this% rootDistance(STLNum), dist, normal )
            if( Dist .lt. this% IBMmask(domains)% x(n)% dist ) then
               this% IBMmask(domains)% x(n)% dist   = dist
               this% IBMmask(domains)% x(n)% normal = normal
            end if
         end do
      end do
#ifdef _HAS_MPI_
      call gatherMaskGeom( this% IBMmask )
#endif
      if( this% HO_IBM ) then
         do n = 1, this% IBMmask(domain)% NumOfObjs
            fID = this% IBMmask(domain)% x(n)% element_index
            i   = this% IBMmask(domain)% x(n)% local_Position(1)
            j   = this% IBMmask(domain)% x(n)% local_Position(2)

            faces(fID)% stencil(i,j)% dist   = this% IBMmask(domain)% x(n)% dist
            faces(fID)% stencil(i,j)% normal = this% IBMmask(domain)% x(n)% normal
         end do
      else
         do n = 1, this% IBMmask(domain)% NumOfObjs
            eID = this% IBMmask(domain)% x(n)% element_index
            i   = this% IBMmask(domain)% x(n)% local_Position(1)
            j   = this% IBMmask(domain)% x(n)% local_Position(2)
            k   = this% IBMmask(domain)% x(n)% local_Position(3)

            elements(eID)% geom% dwall(i,j,k)    = this% IBMmask(domain)% x(n)% dist
            elements(eID)% geom% normal(:,i,j,k) = this% IBMmask(domain)% x(n)% normal
         end do
      end if
      
      do domains = 1, MPI_Process% nProcs
         if( allocated(this% IBMmask(domains)% x) ) deallocate( this% IBMmask(domains)% x )
      end do

      deallocate(this% IBMmask)

   end subroutine IBM_ComputeIBMWallDistance

   subroutine IBM_GetStencil( this, elements, faces )
      use NodalStorageClass, only: NodalStorage
      implicit none

      class(IBM_type), intent(inout) :: this
      type(element),   intent(in)    :: elements(:)
      type(face),      intent(inout) :: faces(:)

      real(kind=RP) :: L, dist, dl, x0(NDIM), h, normal(NDIM), xi(NDIM)
      integer       :: fID, eID, Nxi, Neta, N, i, j, k

      do fID = 1, size(faces)
         associate( f=> faces(fID) )
         if( f% HO_IBM ) then
            Nxi  = f% Nf(1)
            Neta = f% Nf(2)
            N    = max(Nxi,Neta)

            eID = f% elementIDs(1)
            if( eID .eq. 0 ) eID = f% elementIDs(2)

            h = (elements(eID)% geom% Volume)**(1.0_RP/3.0_RP)
            L = 1._RP/(N+1) * ABS(NodalStorage(N)% x(0) - NodalStorage(N)% x(1)) * h

            do i = 0, Nxi; do j = 0, Neta
               dist   = f% stencil(i,j)% dist
               normal = -f% stencil(i,j)% normal
               dl     = dist + sqrt(3.0_RP) * h   
               x0     = f% stencil(i,j)% x + dl * normal

               f% stencil(i,j)% N   = N
               f% stencil(i,j)% xiI = -1.0_RP - 2.0_RP * dl/L
               f% stencil(i,j)% xiB = -1.0_RP - 2.0_RP * (-dist+dl)/L

               call f% stencil(i,j)% build( x0, normal, L, N )
      
            end do; end do
         end if
         end associate
      end do

   end subroutine IBM_GetStencil

   subroutine IBM_HO_IBMstencilState( this, nEqn, elements )
      use PhysicsStorage
      implicit none

      class(IBM_type), intent(inout) :: this
      integer,         intent(in)    :: nEqn
      type(element),   intent(inout) :: elements(:)

      real(kind=RP) :: xi(NDIM), Qs(nEqn)
      integer       :: fID, i, j, k, eID, domain, ierr, n, m
      
      call IBM_HO_GetState( elements, nEqn )
#ifdef _HAS_MPI_
      call GatherHOfacesState( nEqn )
#else
      domain = MPI_Process% rank + 1
      do n = 1, IBMStencilPoints(domain)% NumOfObjs
         m = IBMStencilPoints(domain)% x(n)% faceID
         i = IBMStencilPoints(domain)% x(n)% local_position(IX)
         j = IBMStencilPoints(domain)% x(n)% local_position(IY)
         k = IBMStencilPoints(domain)% x(n)% local_position(IZ)

         IBM_HO_faces(domain)% faces(m)% stencil(i,j)% Q(:,k) = IBMStencilPoints(domain)% x(n)% Q
      end do 
#endif
   end subroutine IBM_HO_IBMstencilState

   subroutine IBM_MPI_GatherStancilState( this, nEqn )
      use MPI_Process_info

      use Meshtypes
      implicit none

      class(IBM_type), intent(inout) :: this
      integer,         intent(in)    :: nEqn

      integer :: domains, fID, i, j, k, domain

      domain = MPI_Process% rank + 1
      
      do fID = 1, IBM_HO_faces(domain)% NumOfFaces
         associate( f => IBM_HO_faces(domain)% faces(fID) )
         do i = 0, f% Nf(1); do j = 0, f% Nf(2)
            call f% stencil(i,j)% ComputeState()
         end do; end do
         end associate
      end do

   end subroutine IBM_MPI_GatherStancilState
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------
!  Source terms for the immersed boundary
!  ------------------------------------------------
   subroutine IBM_SourceTerm( this, eID, Q, Q_target, Source, Wallfunction )
      use PhysicsStorage
      implicit none
      !-arguments--------------------------------------------
      class(IBM_type),           intent(inout) :: this
      integer,                   intent(in)    :: eID
      real(kind=rp),             intent(in)    :: Q(NCONS)
      real(kind=rp),   optional, intent(in)    :: Q_target(NCONS)
      real(kind=rp),             intent(inout) :: Source(NCONS)
      logical,                   intent(in)    :: Wallfunction
      !-local-variables--------------------------------------
      real(kind=rp) :: rho, rho_s, u, u_s, v, v_s, w, w_s
#if defined(SPALARTALMARAS)
      real(kind=rp) :: theta, theta_s
#endif
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

      Source(IRHOU) = rho*u-rho_s*u_s
      Source(IRHOV) = rho*v-rho_s*v_s
      Source(IRHOW) = rho*w-rho_s*w_s
      Source(IRHOE) = 0.5_RP * rho*( POW2(u) + POW2(v) + POW2(w) )  &
                     -0.5_RP * rho_s*( POW2(u_s) + POW2(v_s) + POW2(w_s) )
#if defined(SPALARTALMARAS)
      theta = Q(IRHOTHETA)/rho
      Source(IRHOTHETA) = rho*theta - rho_s*theta_s
#endif
#endif
      if( Wallfunction ) then
         Source = -1.0_RP/(this% penalCoeff * this% penalization(eID)) * Source
      else
         Source = -1.0_RP/this% penalization(eID) * Source
      end if
   end subroutine IBM_SourceTerm

   function IBM_MaskVelocity( this, Q, nEqn, STLNum, x, t ) result( Q_target )
      use PhysicsStorage
      use FluidData
      use TessellationTypes
      use OrientedBoundingBox
      implicit none

      class(IBM_type), intent(inout) :: this
      integer,         intent(in)    :: nEqn, STLNum
      real(kind=RP),   intent(in)    :: Q(nEqn), x(NDIM), t
      real(kind=RP)                  :: Q_target(nEqn)

      real(kind=RP) :: R, v_plane, time, theta, &
                       rho, u_s, v_s, w_s

      Q_target = Q

      u_s = 0.0_RP; v_s = 0.0_RP; w_s = 0.0_RP
#if defined(NAVIERSTOKES)
      if( this% stl(STLNum)% motionType .eq. LINEAR ) then
         select case( this% stl(STLNum)% motionAxis )
         case( IX )
            u_s = this% stl(STLNum)% Velocity/refValues% V
         case( IY )
            v_s = this% stl(STLNum)% Velocity/refValues% V
         case( IZ )
            w_s = this% stl(STLNum)% Velocity/refValues% V
         end select
      elseif( this% stl(STLNum)% motionType .eq. ROTATION ) then
         R = norm2( x - this% stl(STLNum)% rotationCenter )

         v_plane = this% stl(STLNum)% angularVelocity * R * Lref
         v_plane = v_plane/refValues% V

         time    = t * Lref/refValues% V
         theta   = this% stl(STLNum)% angularVelocity * time

         select case( this% stl(STLNum)% motionAxis )
         case( IX )
            v_s = -sin(theta) * v_plane
            w_s =  cos(theta) * v_plane
         case( IY )
            u_s = -cos(theta) * v_plane
            w_s = -sin(theta) * v_plane
         case( IZ )
            u_s = -sin(theta) * v_plane
            v_s =  cos(theta) * v_plane
         end select
      end if

      rho = Q_target(IRHO)

      Q_target(IRHOU) = rho * u_s
      Q_target(IRHOV) = rho * v_s
      Q_target(IRHOW) = rho * w_s
#endif
   end function IBM_MaskVelocity
!
!  Analytical jacobian: dS/dQ
!  ---------------------------
   subroutine IBM_semiImplicitJacobian( this, eID, Q, dS_dQ )
      use PhysicsStorage
      implicit none
      !-arguments----------------------------------------------------
      class(IBM_type), intent(inout) :: this
      integer,         intent(in)    :: eID
      real(kind=rp),   intent(in)    :: Q(NCONS)
      real(kind=rp),   intent(inout) :: dS_dQ(NCONS,NCONS)
      !-local-variables----------------------------------------------
      real(kind=rp) :: rho, u, v, w

      dS_dQ = 0.0_RP
#if defined(NAVIERSTOKES)
      rho = Q(IRHO)
      u   = Q(IRHOU)/rho
      v   = Q(IRHOV)/rho
      w   = Q(IRHOW)/rho

      dS_dQ(IRHOU,IRHOU) = 1.0_RP
      dS_dQ(IRHOV,IRHOV) = 1.0_RP
      dS_dQ(IRHOW,IRHOW) = 1.0_RP
      dS_dQ(IRHOE,IRHO)  = -0.5_RP*( POW2(u) + POW2(v) + POW2(w) )
      dS_dQ(IRHOE,IRHOU) = u
      dS_dQ(IRHOE,IRHOV) = v
      dS_dQ(IRHOE,IRHOW) = w
#if defined(SPALARTALMARAS)
      dS_dQ(IRHOTHETA,IRHOTHETA) = 1.0_RP
#endif
#endif
       dS_dQ = -1.0_RP/this% penalization(eID) * dS_dQ

    end subroutine IBM_semiImplicitJacobian

   subroutine IBM_semiImplicitTurbulenceJacobian( this, e, eID, Q, Q_target, i, j, k, dS_dQt )
      use PhysicsStorage
      implicit none
      !-arguments-----------------------------------------------------------------
      class(IBM_type),  intent(inout) :: this
      type(element),    intent(in)    :: e
      integer,          intent(in)    :: eID
      real(kind=rp),    intent(in)    :: Q(NCONS), Q_target(NCONS)
      integer,          intent(in)    :: i, j, k
      real(kind=rp),    intent(inout) :: dS_dQt(NCONS,NCONS)
      !-local-variables-----------------------------------------------------------
      real(kind=rp) :: Q_IP(NCONS), Q_F(NCONS), &
                       dS_dQ_F(NCONS,NCONS),    &
                       dS_dQ(NCONS,NCONS), eps
      integer       :: n

      Q_IP = Q
      Q_F  = Q_target

      eps    = sqrt(EPSILON(eps))
#if defined(NAVIERSTOKES)
      do n = 1, NCONS
         Q_IP(n) = Q_IP(n) + eps
         Q_F     = ForcingPointState( Q,                                        &
                                      this% IP_Distance, e% geom% dWall(i,j,k), &
                                      e% geom% normal(:,i,j,k), NCONS           )
         dS_dQ_F(:,n) = (Q_F - Q)/eps
         Q_IP = Q
         Q_F  = Q_target
      end do
#endif
      call this% semiImplicitJacobian( eID, Q, dS_dQ )

      dS_dQt = -1.0_RP/(this% penalCoeff*this% penalization(eID)) * (dS_dQ - dS_dQ_F)

   end subroutine IBM_semiImplicitTurbulenceJacobian
!
!  Analytical inverse jacobian: (1 - dt*dS/dQ)^(-1)
!  --------------------------------------------------
   subroutine IBM_semiImplicitShiftJacobian( this, eID, Q, dt, invdS_dQ )
      use PhysicsStorage
      implicit none
         !-arguments----------------------------------------------------------------
         class(IBM_type), intent(inout) :: this
         integer,         intent(in)    :: eID
         real(kind=rp),   intent(in)    :: Q(NCONS)
         real(kind=rp),   intent(in)    :: dt
         real(kind=rp),   intent(inout) :: invdS_dQ(NCONS,NCONS)
         !-local-variables----------------------------------------------------------
         real(kind=rp) :: rho, u, v, w

         invdS_dQ = 0.0_RP

         associate( eta => this% penalization(eID) )
#if defined(NAVIERSTOKES)
         rho = Q(IRHO)
         u   = Q(IRHOU)/rho
         v   = Q(IRHOV)/rho
         w   = Q(IRHOW)/rho

         invdS_dQ(IRHO,IRHO)   = 1.0_RP
         invdS_dQ(IRHOU,IRHOU) = eta/( dt + eta )
         invdS_dQ(IRHOV,IRHOV) = eta/( dt + eta )
         invdS_dQ(IRHOW,IRHOW) = eta/( dt + eta )
         invdS_dQ(IRHOE,IRHO)  = 0.5_RP*dt/eta * (POW2(u) + POW2(v) + POW2(w))
         invdS_dQ(IRHOE,IRHOU) = -dt*u/( dt + eta )
         invdS_dQ(IRHOE,IRHOV) = -dt*v/( dt + eta )
         invdS_dQ(IRHOE,IRHOW) = -dt*w/( dt + eta )
         invdS_dQ(IRHOE,IRHOE) =  1.0_RP
#if defined(SPALARTALMARAS)
         invdS_dQ(IRHOTHETA,IRHOTHETA) = eta/( dt + eta )
#endif
#endif
       end associate

   end subroutine IBM_SemiImplicitShiftJacobian
!
!  (1 - dt*dS_dQt)
!
   subroutine SemiImplicitTurbulenceShiftJacobian( dt, dS_dQt )
      use PhysicsStorage
      implicit none
      !-arguments-----------------------------------
      real(kind=RP),   intent(in)    :: dt
      real(kind=RP),   intent(inout) :: dS_dQt(NCONS,NCONS)
      !-local-variables-----------------------------
      integer :: i

      dS_dQt = -dt*dS_dQt

      do i = 1, NCONS
         dS_dQt(i,i) = 1.0_RP + dS_dQt(i,i)
      end do

   end subroutine SemiImplicitTurbulenceShiftJacobian
!
!  Second order Strang splitting correction Q^*. (1/dt - dS/dQ)^(-1)*Q^* = Q + dt*(S - dS/dQ*Q^*)
!  ------------------------------------------------------------------------------------------------
   subroutine IBM_GetSemiImplicitStep( this, eID, dt, Q, Q_target )
      use PhysicsStorage
      implicit none
      !-arguments-----------------------------------------------------
      class(IBM_type),           intent(inout) :: this
      integer,                   intent(in)    :: eID
      real(kind=rp),             intent(in)    :: dt
      real(kind=rp),             intent(inout) :: Q(NCONS)
      real(kind=rp),   optional, intent(in)    :: Q_target(NCONS)
      !-local-variables-----------------------------------------------
      real(kind=rp) :: dS_dQ(NCONS,NCONS), invdS_dQ(NCONS,NCONS), &
                       IBMSource(NCONS)

      call this% semiImplicitJacobian( eID, Q, dS_dQ )

      call this% semiImplicitShiftJacobian( eID, Q, dt, invdS_dQ )
      if( present(Q_target) ) then
         call this% SourceTerm(eID = eID, Q = Q, Q_target = Q_target, Source = IBMSource, Wallfunction  =.false.)
      else
         call this% SourceTerm(eID = eID, Q = Q, Source = IBMSource, Wallfunction  =.false.)
      end if

      Q = matmul(invdS_dQ, Q + dt*( IBMSource - matmul(dS_dQ,Q) ))

   end subroutine IBM_GetSemiImplicitStep
!
!  Second order Strang splitting correction Q^*. (1/dt - dS/dQ)^(-1)*Q^* = Q + dt*(S - dS/dQ*Q^*)
!  ------------------------------------------------------------------------------------------------
   subroutine IBM_GetSemiImplicitStepTurbulence( this, e, eID, Q, Q_target, i, j, k, dt )
      use PhysicsStorage
      use DenseMatUtilities
      implicit none
      !-arguments-------------------------------------------------------
      class(IBM_type), intent(inout) :: this
      type(element),   intent(in)    :: e
      integer,         intent(in)    :: eID, i, j, k
      real(kind=rp),   intent(in)    :: dt
      real(kind=rp),   intent(inout) :: Q(NCONS), Q_target(NCONS)
      !-local-variables-------------------------------------------------
      real(kind=rp) :: dS_dQ(NCONS,NCONS), invdS_dQ(NCONS,NCONS), &
                       dS_dQt(NCONS,NCONS), Source(NCONS)

      call this% semiImplicitTurbulenceJacobian( e, eID, Q, Q_target, i, j, k, dS_dQt )

      dS_dQ = dS_dQt

      call SemiImplicitTurbulenceShiftJacobian( dt, dS_dQt )

      invdS_dQ = inverse(dS_dQt)

      call this% SourceTerm( eID, Q, Q_target, Source, this% wallfunction )

      Q = matmul(invdS_dQ, Q + dt*( Source - matmul(dS_dQ,Q) ))

    end subroutine IBM_GetSemiImplicitStepTurbulence
!
!   Splitting correction is applied. In the splitting, dt = dt/2
!   ------------------------------------------------------------
   subroutine IBM_SemiImplicitCorrection( this, elements, t, dt )
      use PhysicsStorage
      implicit none
      !-arguments-----------------------------------
      class(IBM_type), intent(inout) :: this
      type(element),   intent(inout) :: elements(:)
      real(kind=RP),   intent(in)    :: t, dt
      !-local-variables-----------------------------
#if defined(NAVIERSTOKES)
      real(kind=RP) :: Q_target(NCONS)
#endif
      integer       :: eID, i, j, k, domain, STLNum

      if( .not. this% semiImplicit ) return
#if defined(NAVIERSTOKES)
      if( this% Wallfunction ) then
         call this% GetBandRegionStates( elements, NCONS )
         call this% GetImagePointsStates( NCONS )
      end if
!$omp parallel
!$omp do schedule(runtime) private(i,j,k,Q_target)
      do eID = 1, SIZE( elements )
         associate(e => elements(eID))
         do i = 0, e% Nxyz(1); do j = 0, e% Nxyz(2); do k = 0, e% Nxyz(3)
            if( e% isInsideBody(i,j,k) ) then
               if( this% stl(e% STL(i,j,k))% move ) then
                  Q_target = this% MaskVelocity( e% storage% Q(:,i,j,k), NCONS, e% STL(i,j,k), e% geom% x(:,i,j,k), t )
                  call this% GetSemiImplicitStep( eID, 0.5_RP*dt, e% storage% Q(:,i,j,k), Q_target )
               else
                  call this% GetSemiImplicitStep( eID, 0.5_RP*dt, e% storage% Q(:,i,j,k) )
               end if
            end if
         end do; end do; end do
         end associate
      end do
!$omp end do
      if( this% Wallfunction ) then
!$omp do schedule(runtime) private(i,j,k,domain,STLNum,Q_target)
         do eID = 1, SIZE( elements )
            associate(e => elements(eID))
            do i = 0, e% Nxyz(1); do j = 0, e% Nxyz(2); do k = 0, e% Nxyz(3)
               if( e% isForcingPoint(i,j,k) ) then
                  STLNum   = e% STL(i,j,k)
                  domain   =  MPI_Process% rank+1
                  Q_target = ForcingPointState( this% ImagePoints(STLNum)% IBMmask(domain)% x(e% IP_index)% Q, &
                                                this% IP_Distance, e% geom% dWall(i,j,k),                      &
                                                e% geom% normal(:,i,j,k), NCONS                                )

                  call this% GetSemiImplicitStepTurbulence( e, eID, this% ImagePoints(STLNum)% IBMmask(domain)% x(e% IP_index)% Q, &
                                                            Q_target, i, j, k, 0.5_RP*dt                                           )
               end if
            end do; end do; end do
            end associate
         end do
!$omp end do
      end if
!$omp end parallel
#endif
    end subroutine IBM_SemiImplicitCorrection
!
!   State to be imposed on the forcing points due to the wall model
!   ---------------------------------------------------------------
#if defined(NAVIERSTOKES)
   function ForcingPointState( Q_IP, y_IP, y_FP, normal, nEqn ) result( Q_FP )
      use PhysicsStorage
      use WallFunctionDefinitions
      use WallFunctionBC
      use VariableConversion
      use FluidData
#if defined(SPALARTALMARAS)
      use SpallartAlmarasTurbulence
#endif
      implicit none
      !-arguments--------------------------------------------------------------
      real(kind=rp), intent(in) :: Q_IP(nEqn), normal(NDIM)
      real(kind=rp), intent(in) :: y_IP, y_FP
      integer,       intent(in) :: nEqn
      real(kind=rp)             :: Q_FP(nEqn)
      !-local-variables--------------------------------------------------------
      real(kind=rp)            :: nu_IP, nu_FP, mu_IP, mu_FP, T_IP, T_FP, &
                                  u_tau, u_IP(NDIM), u_IPt, u_FPt,        &
                                  u_FPn, u_FP(NDIM), u_IP_t(NDIM),        &
                                  tangent(NDIM), yplus_FP, nu_tilde,      &
                                  kappa_IP, kappa_FP
#if defined(SPALARTALMARAS)
      real(kind=rp)            :: Dump, chi, fv1, nu_t
#endif
      call get_laminar_mu_kappa(Q_IP,mu_IP,kappa_IP)
      nu_IP = mu_IP/Q_IP(IRHO)

      u_IP   = Q_IP(IRHOU:IRHOW)/Q_IP(IRHO)

      u_IP_t = u_IP - (dot_product(u_IP,normal) * normal)

      if( almostEqual(norm2(u_IP_t),0.0_RP) ) return

      tangent = u_IP_t/norm2(u_IP_t)

      u_IPt = dot_product(u_IP,tangent)

      u_tau = u_tau_f(u_IPt, y_IP, nu_IP, u_tau0=.1_RP)

      call get_laminar_mu_kappa(Q_FP,mu_FP,kappa_FP)
      nu_FP = mu_FP/Q_FP(IRHO)

      u_FPt = u_plus_f(y_plus_f(y_FP, u_tau, nu_FP)) * u_tau

      u_FPn = dot_product(u_IP,normal) * y_FP/y_IP

      u_FP = u_FPt*tangent + u_FPn*normal

      T_FP = Temperature(Q_IP) + (dimensionless% Pr)**(1._RP/3._RP)/(2.0_RP*thermodynamics% cp) * (POW2(u_IPt) - POW2(u_FPt))
#if defined(SPALARTALMARAS)
      !yplus_FP = y_plus_f( y_FP, u_tau, nu_FP)
      !Dump     = POW2(1.0_RP - exp(-yplus_FP/19.0_RP))

      !chi  = QuarticRealPositiveRoot( 1.0_RP, -kappa*u_tau*y_FP*Dump/nu_FP, 0.0_RP, 0.0_RP, -kappa*u_tau*y_FP*Dump/nu_FP*POW3(SAmodel% cv1) )

      !nu_t = nu_FP * chi
      nu_t = kappa * u_tau * y_FP
#endif
     ! Q_FP(IRHO)  = Pressure(Q_IP)*refvalues% p/(thermodynamics% R * refValues% T * T_FP)
     ! Q_FP(IRHO)  = Q_FP(IRHO)/refvalues% rho

      Q_FP(IRHOU:IRHOW) = Q_FP(IRHO) * u_FP

     ! Q_FP(IRHOE) = pressure(Q_IP)/( thermodynamics% gamma - 1._RP) + 0.5_RP*Q_FP(IRHO)*(POW2(u_FP(1)) + POW2(u_FP(2)) + POW2(u_FP(3)))
      Q_FP(IRHOE) = 0.5_RP*Q_FP(IRHO)*(POW2(u_FP(1)) + POW2(u_FP(2)) + POW2(u_FP(3)))
#if defined(SPALARTALMARAS)
      Q_FP(IRHOTHETA) = Q_FP(IRHO) * nu_t
#endif
   end function ForcingPointState
#endif
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------------------
!  Ray tracing
!  -------------------------------------------------------------
   subroutine IBM_CheckPoint( this, Point, STLNum, NumOfIntersections )
      use MPI_Process_Info
      implicit none
      !-arguments------------------------------------------------------------------
      class(IBM_type), intent(inout) :: this
      real(kind=rp),   intent(inout) :: Point(:)
      integer,         intent(in)    :: STLNum
      integer,         intent(inout) :: NumOfIntersections
      !-local-variables------------------------------------------------------------
      type(KDtree),  pointer     :: tree
      real(kind=rp)              :: RayDirection(NDIM), vecAxis(NDIM)
      integer                    :: Axis, lastIndex
      logical                    :: Upward, OnSurface, OnTriBound, delete
      type(ObjsDataLinkedList_t) :: Integer_List

      RayDirection = 0.0_RP

      OnSurface = .false.; delete = .false.

      axis = OBB(STLNum)% minAxis

      RayDirection(axis) = 1.0_RP

      Upward  = .false.
      vecAxis = OBB(STLNum)% LocVertices(:,7)

      Integer_List = ObjsDataLinkedList_t()

      NumOfIntersections = 0; lastIndex = -1

      do

         call this% root(STLNum)% FindLeaf( Point, tree, .false. )

         if( tree% index .eq. lastIndex ) then
            call this% root(STLNum)% FindLeaf( Point, tree, .true. )
         endif

         lastIndex = tree% index

         call isPointInside( Point, RayDirection, this% root(STLNum)% ObjectsList, tree, Integer_List, NumOfIntersections, OnTriBound )

         if( OnTriBound ) delete = .true.

         if( Upward ) then
            Point(axis) = tree% vertices(axis,1)
            if( Point(axis) .le. vecAxis(axis) ) exit
         elseif( .not. Upward ) then
            Point(axis) = tree% vertices(axis,7)
            if( Point(axis) .ge. vecAxis(axis) ) exit
         end if

      end do

      if( delete ) NumOfIntersections = NumOfIntersections - 1

      call integer_List% Destruct()

   end subroutine IBM_CheckPoint
!
!  Intersection between a ray an a set of triangles
!  ------------------------------------------------
   subroutine isPointInside( Point, RayDirection, ObjectsList, tree, Integer_List, NumOfIntersections, OnTriBound )
      use RealDataLinkedList
      use omp_lib
      implicit none
      !-arguments----------------------------------------------------------------
      real(kind=rp),              intent(in)    :: Point(:), RayDirection(:)
      type(object_type),          intent(in)    :: ObjectsList(:)
      type(KDtree),               intent(inout) :: tree
      type(ObjsDataLinkedList_t), intent(inout) :: Integer_List
      integer,                    intent(inout) :: NumOfIntersections
      logical,                    intent(inout) :: OnTriBound
      !-local-variables----------------------------------------------------------
      logical                        :: Intersect, found
      integer                        :: i, index

      OnTriBound = .false.

      if( tree% NumOfObjs .eq. 0 ) return

      do i = 1, tree% NumOfObjs
         index = tree% ObjsIndeces(i)
         found = integer_List% Check( index )
         
         if( .not. found ) then
            call Integer_List% Add( index )

            call PointIntersectTriangle( Point,ObjectsList(index)% vertices(1)% coords, &
                                               ObjectsList(index)% vertices(2)% coords, &
                                               ObjectsList(index)% vertices(3)% coords, &
                                               RayDirection, Intersect, OnTriBound      )
                                               
            if( Intersect ) then
               NumOfIntersections = NumOfIntersections + 1
            end if
         end if
      end do

   end subroutine isPointInside
!
! This subroutine checks if a ray (RayDirection) starting from a point (Point) intersects
! a triangle in 3D space. If present, the intersection point Q is Q = P + t*RayDirection.
! If there is more than one intersections, the second one is not counted if it has the same t
! as the one previously found.
! See Fast, Minimum Storage Ray/Trinagle INtersection,  Moller TRumbore
!  --------------------------------------------------------------------------------------------
   subroutine PointIntersectTriangle( Point, TriangleVertex1, TriangleVertex2, &
                                      TriangleVertex3, RayDirection,           &
                                      Intersect, OnTriBound                    )
      use MappedGeometryClass
      implicit none
      !-arguments---------------------------------------------------------------------------------------
      real(kind=rp), intent(in)    :: Point(NDIM), RayDirection(NDIM)
      real(kind=rp), intent(in)    :: TriangleVertex1(NDIM), TriangleVertex2(NDIM), TriangleVertex3(NDIM)
      logical,       intent(inout) :: Intersect, OnTriBound
      !-local-variables----------------------------------------------------------------------------------
      real(kind=rp) :: E1vec(NDIM), E2vec(NDIM), Pvec(NDIM), &
                       Qvec(NDIM), Tvec(NDIM), N(NDIM),      &
                       Det, u, v, invDet, t
      logical       :: isInside

      Intersect  = .false.
      OnTriBound = .false.
      t          = -1.0_RP

      E1vec = TriangleVertex2 - TriangleVertex1
      E2vec = TriangleVertex3 - TriangleVertex1
      Tvec   = Point - TriangleVertex1

      call vcross(RayDirection,E2vec,Pvec)
      Det = dot_product( E1vec, Pvec )

      If( almostEqual(Det,0.0_RP) ) return

      call vcross(Tvec,E1vec,Qvec)

      invDet = 1.0_RP/Det

      u = dot_product( Tvec, Pvec )*invDet

      if( u < 0.0_RP .or. u > 1.0_RP ) return

      v = dot_product( RayDirection, Qvec )*invDet

      if( v < 0.0_RP .or. u+v > 1.0_RP ) return

      t  = dot_product( E2vec, Qvec )*invDet
      ! Check if the point lies on the boundaries of the triangle
      !----------------------------------------------------------
      if( almostEqual(u,0.0_RP) .and. ( v .ge. 0.0_RP .and. v .le. 1.0_RP ) ) OnTriBound = .true.
      if( almostEqual(v,0.0_RP) .and. ( u .ge. 0.0_RP .and. u .le. 1.0_RP ) ) OnTriBound = .true.
      if( almostEqual(u+v,1.0_RP) )  OnTriBound = .true.

      if( t .gt. 0.0_RP ) Intersect = .true.

   end subroutine PointIntersectTriangle
!
! Point-triangle intersection (point and triangle are on the same plane)
!  ------------------------------------------------
   logical function isPointInsideTri( Point, TriangleVertex1, TriangleVertex2, &
                                      TriangleVertex3 ) result( isInside )
      use MappedGeometryClass
      implicit none
      !-arguments-----------------------------------------------------------------
      real(kind=rp), intent(in) :: Point(:), TriangleVertex1(:),          &
                                   TriangleVertex2(:), TriangleVertex3(:)
      !-local-variables-----------------------------------------------------------
      real(kind=rp) :: bb(NDIM), E0(NDIM), E1(NDIM), dd(NDIM), &
                       a, b, c, d, e, f, det, s, t
      integer       :: region

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

   subroutine RayAxis( maxAxis, ClipAxis, axis )
      implicit none
      !-arguments----------------------------------
      integer, intent(in)    :: maxAxis, ClipAxis
      integer, intent(inout) :: axis
      !-local-variables----------------------------
      integer :: vecAxis(NDIM), i

      vecAxis = (/1,2,3/)

      vecAxis(maxAxis)  = 0
      if( ClipAxis .ne. 0 ) vecAxis(ClipAxis) = 0

      if( axis .eq. maxAxis .or. axis .eq. ClipAxis ) then
         do i = 1, NDIM
            if(vecAxis(i) .eq. 0) cycle
            axis = i; exit
         end do
      end if

   end subroutine RayAxis
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------
!  Procedures for computing the point-STL distance
!  -------------------------------------------------
!
! This function computes the minimum distance from a point to a triangle in 3D.
! for more ditails see https://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
!  ------------------------------------------------
   subroutine MinimumPointTriDistance( Point, TriangleVertex1, TriangleVertex2, &
                                       TriangleVertex3, dist, IntersectionPoint )
      use MappedGeometryClass
      implicit none
      !-arguments--------------------------------------------------------------------
      real(kind=rp), intent(in)  :: Point(:), TriangleVertex1(:), &
                                    TriangleVertex2(:),        &
                                    TriangleVertex3(:)
      real(kind=rp), intent(out) :: IntersectionPoint(NDIM)
      real(kind=rp), intent(out) :: dist
      !-local-variables--------------------------------------------------------------
      real(kind=rp) :: bb(NDIM), E0(NDIM), E1(NDIM), dd(NDIM), &
                       a, b, c, d, e, f, det, s, t,    &
                       tmp1, tmp0, numer, denom
      integer       :: region

      bb = TriangleVertex1
      E0 = TriangleVertex2 - bb
      E1 = TriangleVertex3 - bb
      dd = bb - Point

      a = dot_product(E0,E0)
      b = dot_product(E0,E1)
      c = dot_product(E1,E1)
      d = dot_product(E0,dd)
      e = dot_product(E1,dd)
      f = dot_product(dd,dd)

      det = a*c - b*b
      s   = b*e - c*d
      t   = b*d - a*e

      if( (s + t) <= det ) then
         if( s < 0.0_RP ) then
            if( t < 0.0_RP ) then
               region = 4
            else
               region = 3
            end if
         else if ( t < 0.0_RP ) then
            region = 5
         else
            region = 0
         end if
      else
         if( s < 0.0_RP ) then
            region = 2
         elseif( t < 0.0_RP ) then
            region = 6
         else
            region = 1
         end if
      end if

      select case( region )
      case( 0 )
         s = s/det
         t = t/det
      case( 1 )
         numer = (c + e) - (b + d)
         if( numer <= 0.0_RP ) then
            s = 0.0_RP
         else
            denom = a - 2.0_RP * b + c
            if( numer >= denom ) then
               s = 1.0_RP
            else
               s = numer/denom
            end if
         end if
         t = 1.0_RP - s
      case( 2 )
         tmp0 = b + d
         tmp1 = c + e
         if( tmp1 > tmp0 ) then
            numer = tmp1 - tmp0
            denom = a - 2.0_RP*b + c
            if( numer >= denom ) then
               s = 1.0_RP
            else
               s = numer/denom
            end if
            t = 1.0_RP - s
         else
            s = 0.0_RP
            if ( tmp1 <= 0.0_RP ) then
               t = 1.0_RP
            elseif( e >= 0.0_RP ) then
               t = 0.0_RP
            else
               t = -e/c
            end if
         end if
      case( 3 )
         s = 0.0_RP
         if( e >= 0.0_RP ) then
            t = 0.0_RP
         elseif( -e >= c ) then
            t = 1.0_RP
         else
            t = -e/c
         end if
      case( 4 )
         if( d < 0.0_RP ) then
            t = 0.0_RP
            if( -d >= a ) then
               s = 1.0_RP
            else
               s = -d/a
            end if
         else
            s = 0.0_RP
            if( e >= 0.0_RP ) then
               t = 0.0_RP
            elseif( -e >= c ) then
               t = 1.0_RP
            else
               t = -e/c
            end if
         end if
      case( 5 )
         t = 0.0_RP
         if( d >= 0.0_RP ) then
            s = 0.0_RP
         elseif( -d >= a ) then
            s = 1.0_RP
         else
            s = -d/a
         end if
      case( 6 )
         tmp0 = b + e
         tmp1 = a + d
         if( tmp1 > tmp0 ) then
            numer = tmp1 - tmp0
            denom = a - 2.0_RP*b + c
            if( numer >= denom ) then
               t = 1.0_RP
            else
               t = numer/denom
            end if
            s = 1.0_RP - t
         else
            t = 0.0_RP
            if( tmp1 <= 0.0_RP ) then
               s = 1.0_RP
            elseif( d >= 0.0_RP ) then
               s = 0.0_RP
            else
               s = -d/a
            end if
         end if
      end select

      IntersectionPoint = TriangleVertex1 + s*E0 + t*E1

      dist = norm2(Point - IntersectionPoint)

   end subroutine MinimumPointTriDistance
!
!  Identification of the section, see below.
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
      !-arguments-----------------------------
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
! This subroutine computes the minumum distance from x (global ref. frame) to the body. First,
! the box (tree) where x lies is found, then the min distance between the objects inside the tree and
! the point is computed. If the sphere whose radius is minDist, is enclosed in the box, the subroutine stops.
! If the latter condition is not true, all the boxes that intersects the sphere are checked in the same way as
! the initial one, getting new_minDist. If a lower distance is found, minDist is updated.
! minDist is the SQUARE of the actual distance.
!  -------------------------------------------------------------------------------------------------------------
   subroutine MinimumDistance( Point, root, minDist, normal )
      implicit none
      !-arguments----------------------------------------------------------
      real(kind=rp), intent(in)    :: Point(:)
      type(KDtree),  intent(inout) :: root
      real(kind=rp), intent(inout) :: minDist
      real(kind=rp), intent(out)   :: normal(NDIM)
      !-local-variables-----------------------------------------------------
      real(kind=rp)         :: IntersPoint(NDIM), new_IntersPoint(NDIM), &
                               dsvec(NDIM), x(NDIM), P(NDIM), IP(NDIM),  &
                               IntersectionPoint(NDIM), Dist,            &
                               New_minDist, Radius, ds
      logical               :: Intersect
      type(KDtree), pointer :: tree, newGuess
      integer               :: i, index, LeafIndex,             &
                               TriangleIndex, New_TriangleIndex

      minDist = huge(1.0_RP)

      call root% FindLeaf( Point, tree, .false. )

      LeafIndex = tree% index

      do i = 1, tree% NumOfObjs
         index = tree% ObjsIndeces(i)

         call MinimumPointTriDistance( Point, root% ObjectsList(index)% vertices(1)% coords, &
                                       root% ObjectsList(index)% vertices(2)% coords,        &
                                       root% ObjectsList(index)% vertices(3)% coords, Dist,  &
                                       IntersPoint                                           )

         if( Dist .lt. minDist ) then
            minDist           = Dist
            IntersectionPoint = IntersPoint
            TriangleIndex     = index
         end if
      end do

      if( tree% NumOfObjs .gt. 0 ) then
      ! Check the sphere
      !-----------------
         Radius    = POW2(minDist)
         Intersect = CheckHypersphere( tree, Point, Radius )
      else
         Intersect = .true.
      end if

      nullify(tree)

      if( Intersect ) then
         New_minDist = huge(1.0_RP)
         call MinimumDistOtherBoxes( Point, root, root, Radius, LeafIndex, New_minDist, New_IntersPoint, New_TriangleIndex )
         if( New_minDist .lt. minDist ) then
            minDist           = New_minDist
            IntersectionPoint = New_IntersPoint
            TriangleIndex     = New_TriangleIndex
         end if
      end if

      call OBB(root% STLNum)% ChangeRefFrame( Point, GLOBAL, P )
      call OBB(root% STLNum)% ChangeRefFrame( IntersectionPoint, GLOBAL, IP )

      normal = (P - IP)/norm2(P - IP)

   end subroutine MinimumDistance
!
!  Distance between the point and triangles in other boxes, i.e. not the one containing it
!  -----------------------------------------------------------------------------------------
   recursive subroutine MinimumDistOtherBoxes( Point, root, tree, Radius, LeafIndex,           &
                                               New_minDist, New_IntersPoint, New_TriangleIndex )

      implicit none
      !-arguments-----------------------------------------------
      real(kind=rp), intent(in)    :: Point(:)
      type(KDtree),  intent(in)    :: root
      type(KDtree),  intent(inout) :: tree
      real(kind=rp), intent(in)    :: Radius
      integer,       intent(in)    :: LeafIndex
      real(kind=rp), intent(inout) :: New_minDist
      real(kind=rp), intent(inout) :: New_IntersPoint(:)
      integer,       intent(inout) :: New_TriangleIndex
      !-local-variables-----------------------------------------
      real(kind=rp)                  :: IntersPoint(NDIM), Dist
      logical                        :: Intersect
      integer                        :: i, index

      Intersect = BoxIntersectSphere( Radius, Point, tree% vertices )

      if( Intersect ) then
         if( tree% isLast ) then
            if( LeafIndex .ne. tree% index ) then
               do i = 1, tree% NumOfObjs
                  index = tree% ObjsIndeces(i)
                  call MinimumPointTriDistance( Point, root% ObjectsList(index)% vertices(1)% coords, &
                                                root% ObjectsList(index)% vertices(2)% coords,        &
                                                root% ObjectsList(index)% vertices(3)% coords, Dist,  &
                                                IntersPoint                                           )
                  if( Dist .lt. New_minDist ) then
                     New_minDist       = Dist
                     New_IntersPoint   = IntersPoint
                     New_TriangleIndex = index
                  end if
               end do
            end if
         else
            call MinimumDistOtherBoxes( Point, root, tree% child_L,        &
                                        Radius, LeafIndex, New_minDist,    &
                                        New_IntersPoint, New_TriangleIndex )
            call MinimumDistOtherBoxes( Point, root, tree% child_R,        &
                                        Radius, LeafIndex, New_minDist,    &
                                        New_IntersPoint, New_TriangleIndex )
         end if
      end if

   end subroutine MinimumDistOtherBoxes
!
!  Once the first distance is found, a sphere whose radius is equal to the distance is computed. If
!  it lies inside the starting box, nothing is done; otherwise the boxes intersecting the sphere are
!  selected
!  ---------------------------------------------------------------------------------------------------
   logical function CheckHypersphere( tree, Point, minDist ) result( Intersect )

      implicit none
      !-arguments-----------------------------------------------------
      real(kind=rp),         intent(in)    :: Point(:)
      type(KDtree),  target, intent(inout) :: tree
      real(kind=rp),         intent(inout) :: minDist

      Intersect = BoxIntersectSphere( minDist, Point, tree% vertices )

   end function CheckHypersphere

!  Point-box distance
!  ------------------------------------------------
   real(kind=rp) function PointBoxDistance( Point, vertices ) result( sqrDist )

      implicit none
      !-arguments--------------------------------
      real(kind=rp), intent(in) :: Point(:)
      real(kind=rp), intent(in) :: vertices(:,:)
      !-local-variables--------------------------
      integer :: i

      sqrDist = 0.0_RP

      ! if the point's x-coordinate (xp) is less than the x-min coord. of the box, then
      ! the minimum distance on the x-dir is the distance between the xp and x-min, or
      ! vertices(i,1) - x-coordinate. The following loop avoids the computation of the point
      ! with the minimum distance from P. If the point P is inside the box, sqrDist = 0.
      !--------------------------------------------------------------------------------------
      do i = 1, NDIM
         if( Point(i) .lt. minval(vertices(i,:)) ) sqrDist = sqrDist + POW2(minval(vertices(i,:))-Point(i))
         if( Point(i) .gt. maxval(vertices(i,:)) ) sqrDist = sqrDist + POW2(Point(i)-maxval(vertices(i,:)))
      end do

   end function PointBoxDistance
!
!  Sphere-box intersection, the radius must be SQUARED
!  ------------------------------------------------
   logical function BoxIntersectSphere( Radius, Center, vertices ) result( Intersect )

      implicit none
      !-arguments-------------------------------------------
      real(kind=rp), intent(in) :: vertices(:,:)
      real(kind=rp), intent(in) :: Center(:)
      real(kind=rp), intent(in) :: Radius
      !-local-variables-------------------------------------
      real(kind=rp) :: sqrDist

      Intersect = .false.
      sqrDist = PointBoxDistance( Center, vertices )

      if( sqrDist .le. Radius ) Intersect = .true.

   end function BoxIntersectSphere
!
!  Nearest-neighbor algorithm for selecting the DoFs closer to Point.
!  --------------------------------------------------------------------
   subroutine MinimumDistancePoints( Point, root, BandRegion, minDist, actualIndex, PointsIndex )
      use ElementClass
      implicit none
      !-arguments-------------------------------------------------------
      real(kind=rp),   intent(in)    :: Point(:)
      type(KDtree),    intent(inout) :: root
      type(IBMpoints), intent(in)    :: BandRegion
      real(kind=rp),   intent(inout) :: minDist
      integer,         intent(in)    :: actualIndex
      integer,         intent(inout) :: PointsIndex(:)
      !-local-variables-------------------------------------------------
      real(kind=rp)         :: BandPoint(NDIM), sqrDist, &
                               New_sqrDist, Radius
      logical               :: Intersect
      type(KDtree), pointer :: tree
      integer               :: i, LeafIndex, new_index, k

      minDist = huge(1.0_RP)

      call root% FindLeaf( Point, tree, .false. )

      LeafIndex = tree% index

      do i = 1, tree% NumOfObjs

         if( any(PointsIndex .eq. tree% ObjsIndeces(i)) ) cycle

         BandPoint = BandRegion% x(tree% ObjsIndeces(i))% coords

         sqrDist = 0.0_RP
         do k = 1, NDIM
            sqrDist = sqrDist + POW2(Point(k) - BandPoint(k))
         end do

         if( sqrDist .lt. minDist .or. almostEqual(sqrDist,minDist) ) then
               minDist                  = sqrDist
               PointsIndex(actualIndex) = tree% ObjsIndeces(i)
         end if

      end do

      if( tree% NumOfObjs .gt. 0 ) then
      ! Check the sphere
      !-----------------
         if( PointsIndex(actualIndex) .eq. 0 ) then
            BandPoint = BandRegion% x(PointsIndex(actualIndex-1))% coords
         else
            BandPoint = BandRegion% x(PointsIndex(actualIndex))% coords
         end if

         sqrDist = 0.0_RP
         do k = 1, NDIM
            sqrDist = sqrDist + POW2(Point(k) - BandPoint(k))
         end do

         Radius = sqrDist  !minDist
         Intersect = CheckHypersphere( tree, Point, Radius )
      else
         print *, "IBM:: MinimumDistance: "
         print *, "Can't find triangles in leaf ", LeafIndex
         error stop
      end if

      nullify(tree)

      if( Intersect ) then
         New_sqrDist = huge(1.0_RP)
         call MinimumDistOtherBoxesPoints( Point, root, BandRegion, Radius, &
                                           New_sqrDist, PointsIndex,        &
                                           LeafIndex, new_index             )
         if( New_sqrDist .le. minDist ) then
            minDist = New_sqrDist; PointsIndex(actualIndex) = new_index
         end if
      end if

      minDist = sqrt(minDist)

   end subroutine MinimumDistancePoints
!
!  Distance between the point and DoFs in other boxes, i.e. not the one containing it
!  -----------------------------------------------------------------------------------------
   recursive subroutine MinimumDistOtherBoxesPoints( Point, tree, BandRegion, Radius, &
                                                     New_sqrDist,  PointsIndex,       &
                                                     LeafIndex,  new_index            )
      use elementClass
      implicit none
      !-arguments---------------------------------------------------
      real(kind=rp),   intent(in)    :: Point(:)
      type(KDtree),    intent(inout) :: tree
      type(IBMpoints), intent(in)    :: BandRegion
      real(kind=rp),   intent(in)    :: Radius
      real(kind=rp),   intent(inout) :: New_sqrDist
      integer,         intent(inout) :: new_index
      integer,         intent(in)    :: PointsIndex(:)
      integer,         intent(in)    :: LeafIndex
      !-local-variables---------------------------------------------
      real(kind=rp) :: sqrDist, BandPoint(NDIM)
      logical       :: Intersect
      integer       :: i, k

      Intersect = BoxIntersectSphere( Radius, Point, tree% vertices )

      if( Intersect ) then
         if( tree% isLast ) then
            if( LeafIndex .ne. tree% index ) then
               do i = 1, tree% NumOfObjs

                  if( any(PointsIndex .eq. tree% ObjsIndeces(i)) ) cycle

                  BandPoint = BandRegion% x(tree% ObjsIndeces(i))% coords

                  sqrDist = 0.0_RP
                  do k = 1, NDIM
                     sqrDist = sqrDist + POW2(Point(k) - BandPoint(k))
                  end do

                  if( sqrDist .lt. New_sqrDist .or. almostEqual(sqrDist,New_sqrDist)  )then
                     New_sqrDist = sqrDist
                     new_index   = tree% ObjsIndeces(i)
                  end if

               end do
            end if
         else
            call MinimumDistOtherBoxesPoints( Point, tree% child_L, BandRegion, &
                                              Radius, New_sqrDist,PointsIndex,  &
                                              LeafIndex, new_index              )
            call MinimumDistOtherBoxesPoints( Point, tree% child_R, BandRegion, &
                                              Radius, New_sqrDist, PointsIndex, &
                                              LeafIndex, new_index              )
         end if
      end if

   end subroutine MinimumDistOtherBoxesPoints


   subroutine GetMatrixInterpolationSystem( Point, x, invPhi, b, INTERPOLATION )
      use DenseMatUtilities
      implicit none
      !-arguments------------------------------------
      real(kind=RP),     intent(in)    :: Point(:)
      type(point_type),  intent(in)    :: x(:)
      real(kind=RP),     intent(inout) :: invPhi(:,:), b(:)
      integer,           intent(in)    :: INTERPOLATION
      !-local-variables------------------------------
      real(kind=RP) :: Phi(size(x),size(x)),   &
                       dist(size(x),size(x)),  &
                       L(size(x)), d,          &
                       P(size(x),size(x)),     &
                       P_T(size(x),size(x)),   &
                       V(2*size(x),2*size(x)), &
                       lambda(size(x)),        &
                       pp(size(x)), w(size(b))
      integer       :: i, j, k

      select case( INTERPOLATION )
         case( EXPONENTIAL )

            do i = 1, size(x)
               do j = i, size(x)
                  dist(i,j) = norm2(x(i)% coords - x(j)% coords)
                  dist(j,i) = dist(i,j)
               end do
            end do

            expCoeff = 0.001_RP!norm2(Point - x(1)% coords)

            do i = 1, size(x)
               do j = i, size(x)
                  Phi(i,j) = interpolationfunction(dist(i,j), EXPONENTIAL )
                  Phi(j,i) = Phi(i,j)
               end do
               d    = norm2(Point - x(i)% coords)
               b(i) = interpolationfunction(d, EXPONENTIAL)
            end do

            invPhi = inverse(Phi)

         case( IDW )

            b = 0.0_RP; invPhi = 0.0_RP
            do i = 1, size(x)
               d           = norm2(Point - x(i)% coords)
               invPhi(i,i) = 1.0_RP/d
               b           = b + invPhi(i,i)
            end do
            do i = 1, size(x)
               b(i) = 1.0_RP/b(i)
            end do

         case( POLYHARMONIC_SPLINE )

            b = 0.0_RP; invPhi = 0.0_RP; Phi = 0.0_RP

            do i = 1, size(x)
               do j = i, size(x)
                  dist(i,j) = norm2(x(i)% coords - x(j)% coords)
                  dist(j,i) = dist(i,j)
               end do
            end do

            do i = 1, size(x)
               do j = i, size(x)
                  Phi(i,j) = interpolationfunction(dist(i,j), POLYHARMONIC_SPLINE )
                  Phi(j,i) = Phi(i,j)
               end do
               d    = norm2(Point - x(i)% coords)
               b(i) = interpolationfunction(d, POLYHARMONIC_SPLINE)
               call PolynomialVector( x(i)% coords(1), x(i)% coords(2), x(i)% coords(3), P(i,:) )
            end do

            P_T = transpose(P)

            call buildMatrixPolySpline( Phi, P, P_T, V )

            invPhi = inverse(V)

         case( POLYNOMIAL )

            Phi = 0.0_RP; invPhi = 0.0_RP; b = 0.0_RP
            do i = 1, size(x)
               call PolynomialVector( x(i)% coords(1), x(i)% coords(2), x(i)% coords(3), Phi(i,:) )
            end do

            call PolynomialVector( Point(1), Point(2), Point(3), b )

            invPhi = inverse(Phi)

         case( MLS )

            Phi = 0.0_RP; invPhi = 0.0_RP; b = 0.0_RP
            expCoeff = norm2(Point - x(1)% coords)
            do i = 1, size(x)
               call PolynomialVector( x(i)% coords(1)-Point(1), x(i)% coords(2)-Point(2), x(i)% coords(3)-Point(3), P(i,:) )
               d    = norm2(Point - x(i)% coords)
               W(i) = interpolationfunction( d, MLS )
            end do

            call PolynomialVector( 0.0_RP, 0.0_RP, 0.0_RP, pp )

            do j = 1, size(p,2)
               do k = 1, size(p,2)
                  do i = 1, size(x)
                     Phi(j,k) = Phi(j,k) + P(i,j)*P(i,k)*W(i)
                  end do
               end do
            end do

            lambda = matmul(inverse(Phi),pp)

            do i = 1, size(p,2)
               do j = 1, size(p,2)
                  b(i) = b(i) + lambda(j) * p(i,j)
               end do
               b(i) = W(i)*b(i)
            end do

      end select

   end subroutine GetMatrixInterpolationSystem

   real(kind=RP) function interpolationfunction( x, INTERPOLATION, k ) result(f)
      implicit none
      real(kind=RP),           intent(in) :: x
      integer,                 intent(in) :: INTERPOLATION
      real(kind=rp), optional, intent(in) :: k

      select case( INTERPOLATION )
         case( EXPONENTIAL )

            f = exp(-POW2(x/expCoeff))

         case( POLYHARMONIC_SPLINE )

            if( abs(x) .lt. 1.0d-12 ) then
               f = 0.0_RP
               return
            end if
            if( present(k) ) then
               f = x**k * log(x)
            else
               f = POW2(x) * log(x)
            end if

         case( MLS )

            f = 1.0_RP/(POW2(x) + 1.0d-12)

      end select
   end function interpolationfunction

   real(kind=RP) function GetInterpolatedValue( forcing, invPhi, b, INTERPOLATION, x, y, z ) result( value )
   implicit none
   !-arguments------------------------------------------------------------
   real(kind=RP),           intent(in) :: forcing(:), invPhi(:,:), b(:)
   integer,                 intent(in) :: INTERPOLATION
   real(kind=RP), optional, intent(in) :: x, y, z
   !-local-variables------------------------------------------------------
   real(kind=RP), allocatable :: weights(:), f(:)

   select case( INTERPOLATION )
      case( POLYHARMONIC_SPLINE )

         allocate( weights(2*size(b)), &
                   f(2*size(b))        )
         f            = 0.0_RP
         f(1:size(b)) = forcing
         weights = matmul(invPhi,f)
         value = dot_product(weights(1:size(b)),b) +                   &
                 EvaluateInterp( weights(size(b)+1:2*size(b)), x, y, z )
         deallocate(weights,f)
      case( MLS )

         value = dot_product(forcing,b)

      case default
         allocate(weights(size(b)))
         weights = matmul(invPhi,forcing)
         value = dot_product(weights,b)
         deallocate(weights)
   end select

   end function GetInterpolatedValue

! estimate the yplus for a flat plate

   real(kind=RP) function InitializeDistance( y_plus ) result( y )
      use FluidData
      implicit none
      !-arguments--------------------------
      real(kind=rp), intent(in) :: y_plus
      !-local-varirables-------------------
      real(kind=RP) :: nu, Cf

      y = 0.0_RP
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
      Cf = 0.058_RP * dimensionless% Re**(-0.2_RP)
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
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
!  ---------------------------------------------------
!  Subourinte for generating a .tec file with forces
!  ---------------------------------------------------
   subroutine STLScalarTEC( x, y, z, scalar, STLNum, FileName, title, variables )

      implicit none
      !-arguments--------------------------------------------------
      real(kind=RP),     intent(in) :: x(:), y(:), z(:), scalar(:)
      integer,           intent(in) :: STLNum
      character(len=*),  intent(in) :: FileName, title, variables
      !-local-variables--------------------------------------------
      integer :: NumOfObjs, funit

      NumOfObjs = size(x)

      call TecHeader( funit, FileName, title, variables )
      call WriteScalarQuatity( funit, NumOfObjs, x, y, z, scalar )
      call TECtriangle_3points( funit, NumOfObjs )

      close(funit)

   end subroutine STLScalarTEC

   subroutine STLvectorTEC( x, y, z, vector_x, vector_y, vector_z, STLNum, FileName, title, variables )

      implicit none
      !-arguments---------------------------------------------------------
      real(kind=RP),     intent(in) :: x(:), y(:), z(:), vector_x(:), &
                                       vector_y(:), vector_z(:)
      integer,           intent(in) :: STLNum
      character(len=*),  intent(in) :: FileName, title, variables
      !-local-variables--------------------------------------------------
      integer                    :: NumOfObjs, funit

      NumOfObjs = size(x)

      call TecHeader( funit, FileName, title, variables )
      call WriteVectorQuatity( funit, NumOfObjs, x, y, z, vector_x, vector_y, vector_z )
      call TECtriangle_3points( funit, NumOfObjs )

      close(funit)

   end subroutine STLvectorTEC

   subroutine TecHeader( funit, FileName, title, variables )

      implicit none
      !-arguments-------------------
      integer,          intent(inout) :: funit
      character(len=*), intent(in)    :: FileName, title, &
                                         variables

      funit = UnusedUnit()

      open(funit,file='RESULTS/'//trim(FileName), status='unknown')

      write(funit,'(A)') 'Title="' //trim(title)// '"'
      write(funit,'(A)') 'Variables='//trim(variables)

   end subroutine TecHeader

   subroutine WriteScalarQuatity( funit, NumOfObjs, x, y, z, scalar )

      implicit none
      !-arguments-------------------------------------------------------
      integer,       intent(in) :: funit, NumOfObjs
      real(kind=RP), intent(in) :: x(:), y(:), z(:), scalar(:)
      !-local-variables-------------------------------------------------
      integer :: i, j

      write(funit,'(a,i6,a,i6,a)') 'ZONE NODES =', NumOfObjs, ',ELEMENTS =', NumOfObjs/3,', DATAPACKING=POINT, ZONETYPE=FETRIANGLE'

      do i = 1, NumOfObjs
         write(funit,'(4g15.6)') x(i), y(i), z(i), scalar(i)
      end do

   end subroutine WriteScalarQuatity

   subroutine WriteVectorQuatity( funit, NumOfObjs, x, y, z, vector_x, vector_y, vector_z )

      implicit none
      !-arguments---------------------------------------------------------
      integer,       intent(in) :: funit, NumOfObjs
      real(kind=RP), intent(in) :: x(:), y(:), z(:), vector_x(:), &
                                   vector_y(:), vector_z(:)
      !-local-variables---------------------------------------------------
      integer :: i, j

      write(funit,'(a,i6,a,i6,a)') 'ZONE NODES=', NumOfObjs, ',ELEMENTS =', NumOfObjs/3,', DATAPACKING=POINT, ZONETYPE=FETRIANGLE'

      do i = 1, NumOfObjs
         write (funit,'(6g15.6)') x(i), y(i), z(i), vector_x(i), vector_y(i), vector_z(i)
      end do

   end subroutine WriteVectorQuatity

   subroutine WriteTimeFile( t, STLNum )
      use MPI_Process_Info
      implicit none
      !-arguments---------------------------------------------------------
      real(kind=RP), intent(in) :: t
      integer ,      intent(in) :: STLNum
      !-local-variables---------------------------------------------------
      integer :: funit

      if( .not. MPI_Process% isRoot ) return

      funit = UnusedUnit()

      open(funit,file='RESULTS/'//trim(OBB(STLNum)% FileName)//'_time.dat', action = "write" , access = "append" , status = "unknown")

      write(funit,'(1g15.6)') t

      close(funit)

   end subroutine WriteTimeFile

   subroutine TECtriangle_3points( funit, NumOfObjs )

      implicit none
      !-arguments------------------------------------
      integer, intent(in) :: funit, NumOfObjs
      !-local-variables------------------------------
      integer :: i, index

      index = 0

      do i = 1, NumOfObjs/3
         write(funit,'(i0,a,i0,a,i0)') index + 1,' ', index + 2,' ', index + 3
         index = index + 3
      end do

   end subroutine TECtriangle_3points

end module IBMClass