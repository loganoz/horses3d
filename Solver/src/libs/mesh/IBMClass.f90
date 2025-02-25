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
   use BoundaryConditions
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
      type(KDtree),               allocatable :: root(:)
      type(bandRegion_t),         allocatable :: BandRegion(:), ImagePoints(:)
      type(Integral_t),           allocatable :: Integral(:)
      character(len=LINE_LENGTH), allocatable :: STLfilename(:)
      logical,                    allocatable :: zoneMask(:)
      character(len=LINE_LENGTH)              :: filename
      logical                                 :: plotOBB              = .false., &
                                                 plotKDtree           = .false., &
                                                 active               = .false., &
                                                 TimePenal            = .false., &
                                                 semiImplicit         = .false., &
                                                 ComputeBandRegion    = .false., &
                                                 plotBandPoints       = .false., &
                                                 plot                 = .false., &
                                                 ComputeInterpolation = .false., &
                                                 Wallfunction         = .false., &
                                                 ComputeDistance      = .false., &
                                                 AAB                  = .false., &
                                                 HO_IBM               = .false., &
                                                 Implicit             = .false., &
                                                 autosave             = .false., &
                                                 twoD                 = .false., &
                                                 MPIfixed             = .false.
      real(kind=rp)                           :: eta, BandRegionCoeff, IP_dWall = 0.0_RP, &
                                                 y_plus_target, minCOORDS, maxCOORDS,        &
                                                 penalCoeff, t = 0.0_RP, dt, dl, d, L
      real(kind=rp),              allocatable :: penalization(:)
      integer                                 :: KDtree_Min_n_of_Objs, NumOfInterPoints,     &
                                                 n_of_INpoints,  rank, lvl = 0, NumOfSTL,    &
                                                 NumOfForcingPoints, Clipaxis = 0,           &
                                                 Nx, Ny, Nz, LocClipAxis = 0,                &
                                                 InterpolationType, NumOfMaskObjs = 0,       &
                                                 iter = 0, twoD_axis = 0, N = 0
      integer,                    allocatable :: bctype(:)
      type(IBMPoints),            allocatable :: IBMmask(:), IBMStencilPoints(:), IBM_HOIntegrationPoints(:)

      contains
         procedure :: read_info                           => IBM_read_info
         procedure :: construct                           => IBM_construct
         procedure :: constructMask                       => IBM_constructMask
         procedure :: constructzonemask                   => IBM_constructzonemask
         procedure :: constructSTL_KDtree                 => IBM_constructSTL_KDtree
         procedure :: CheckPoint                          => IBM_checkPoint
         procedure :: constructBandRegion                 => IBM_constructBandRegion
         procedure :: build                               => IBM_build
         procedure :: SetPolynomialOrder                  => IBM_SetPolynomialOrder
         procedure :: MPI_sendSTLpartitions               => IBM_MPI_sendSTLpartitions
         procedure :: GetForcingPointsGeom                => IBM_GetForcingPointsGeom
         procedure :: GetInfo                             => IBM_GetInfo
         procedure :: minDistance                         => IBM_minDistance
         procedure :: GetPointInterpolation               => IBM_GetPointInterpolation
         procedure :: SourceTerm                          => IBM_SourceTerm
         procedure :: TurbulentSourceTerm                 => IBM_TurbulentSourceTerm
         procedure :: ComputeIBMWallDistance              => IBM_ComputeIBMWallDistance
         procedure :: SemiImplicitCorrection              => IBM_SemiImplicitCorrection
         procedure :: GetDomainExtreme                    => IBM_GetDomainExtreme
         procedure :: semiImplicitShiftJacobian           => IBM_semiImplicitShiftJacobian
         procedure :: semiImplicitJacobian                => IBM_semiImplicitJacobian
         procedure :: semiImplicitTurbulenceJacobian      => IBM_semiImplicitTurbulenceJacobian
         procedure :: GetSemiImplicitStep                 => IBM_GetSemiImplicitStep
         procedure :: GetImplicitStep                     => IBM_GetImplicitStep
         procedure :: GetSemiImplicitStepTurbulence       => IBM_GetSemiImplicitStepTurbulence
         procedure :: copy                                => IBM_copy
         procedure :: MoveBody                            => IBM_MoveBody
         procedure :: Describe                            => IBM_Describe
         procedure :: plotMask                            => IBM_plotMask
         procedure :: Destruct                            => IBM_Destruct
         procedure :: DestroyKDtree                       => IBM_DestroyKDtree
         procedure :: DestroyBandRegion                   => IBM_DestroyBandRegion
         procedure :: HOmask                              => IBM_HOmask
         procedure :: GetStencil                          => IBM_GetStencil
         procedure :: HO_IBMstencilState                  => IBM_HO_IBMstencilState
         procedure :: HO_IBMstencilGradient               => IBM_HO_IBMstencilGradient
         procedure :: MPI_GatherStancilState              => IBM_MPI_GatherStancilState
         procedure :: buildHOfaces                        => IBM_buildHOfaces
         procedure :: buildHOIntegrationPoints            => IBM_buildHOIntegrationPoints
         procedure :: constructImagePoint                 => IBM_constructImagePoint
         procedure :: updateImagePoint                    => IBM_updateImagePoint
         procedure :: RelaxingSourceTerm                  => IBM_RelaxingSourceTerm
   end type

   public :: expCoeff, EXPONENTIAL, GetPointState, GetPointGrads

   real(kind=RP)      :: expCoeff
   integer, parameter :: EXPONENTIAL = 1, IDW = 2

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

      NumOfObjs = sum(this% IBMmask(1:MPI_Process% nProcs)% NumOfObjs)

      if( NumOfObjs .eq. 0 ) return

      call TecFileHeader( 'IBM/'//trim(filename), 'Points', NumOfObjs, funit )

      do domains = 1, MPI_Process% nProcs
         do i = 1, this% IBMmask(domains)% NumOfObjs
            write(funit,'(3E13.5)')  this% IBMmask(domains)% coords(i,IX), this% IBMmask(domains)% coords(i,IY), this% IBMmask(domains)% coords(i,IZ)
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
                                    HO_in, Implicit_in
      real(kind=rp), allocatable :: penalization_in, y_plus_target_in,  &
                                    BandRegionCoeff_in,                 &
                                    penalization_coeff_in
      integer,       allocatable :: n_of_Objs_in, n_of_interpoints_in,  &
                                    Nx_in, Ny_in, Nz_in, Clipaxis_in,   &
                                    interval_in, twoD_axis_in
      character(len=LINE_LENGTH) :: in_label, paramFile, name_in, tmp,  &
                                    InterpolationType_in
      real(kind=rp), allocatable :: coords(:)
      logical                    :: correct

      character(len=LINE_LENGTH), parameter :: NumberOfSTL = "number of stl"

!     Read block
!     **********
      write(in_label , '(A)') "#define ibm"
      call get_command_argument(1, paramFile)
      call readValueInRegion( trim( paramFile ), "active",                         active_in,             in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "penalization",                   penalization_in,       in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "penalization coeff",             penalization_coeff_in, in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "semi implicit",                  semiImplicit_in,       in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "implicit",                       Implicit_in,           in_label, "#end" )
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
      call readValueInRegion( trim( paramFile ), "plot band points",               plotBandPoints_in,     in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "high order",                     HO_in,                 in_label, "#end" )
      call readValueInRegion( trim( paramFile ), "2d axis",                        twoD_axis_in,          in_label, "#end" )

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

      if( allocated(Implicit_in) ) then
         this% Implicit = Implicit_in
      else
         this% Implicit = .FALSE.
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
         this% KDtree_Min_n_of_Objs = 10
      end if

      if( allocated(n_of_interpoints_in) ) then
         this% NumOfInterPoints = n_of_interpoints_in
      else
         this% NumOfInterPoints = 5
      end if

      if( allocated(plotBandPoints_in) ) then
         this% plotBandPoints = plotBandPoints_in
      else
         this% plotBandPoints = .false.
      end if

      this% plot = .true.

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

      if( allocated(twoD_axis_in) ) then 
         this% twoD     = .true. 
         this% twoD_axis = twoD_axis_in  
      else 
         this% twoD = .false. 
      end if 

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
         if( .not. this% HO_IBM )  then
            this% ComputeBandRegion = .true.
            this% ComputeDistance   = .true.
         end if
      else
         this% Wallfunction = .false.
      end if

      select case(trim(InterpolationType_in))
         case("exp")
            this% InterpolationType =  EXPONENTIAL
         case("idw")
            this% InterpolationType =  IDW
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
      real(kind=RP)              :: axis(NDIM), Area, localVal
      integer                    :: STLNum, j, k, maxAxis, minAxis, ierr
      logical                    :: describe = .true.

      if ( .not. allocated(OBB) ) then 
         call this% describe()
         allocate(OBB(this% NumOfSTL))
         describe = .true.
      else 
         describe = .false.
      end if

      allocate( this% stl(this% NumOfSTL),         &
                this% root(this% NumOfSTL),        &
                this% integral(this% NumOfSTL),    &
                this% STLfilename(this% NumOfSTL), &
                this% bctype(this% NumOfSTL),      &
                this% zoneMask(this% NumOfSTL)     )

      if( this% ComputeBandRegion ) allocate( this% BandRegion(this% NumOfSTL)   )

      this% zoneMask = .false. 

      if( this% Wallfunction .and. .not. this% HO_IBM ) allocate( this% ImagePoints(this% NumOfSTL)  )
      
      do STLNum = 1, this% NumOfSTL
         write(MyString, '(i100)') STLNum
         if( STLNum .eq. 1 ) then
            filename = stlFileNameKey
         else
            filename = trim(stlFileNameKey)//trim(adjustl(MyString))
         end if

         this% stl(STLNum)% filename = controlVariables% stringValueForKey(trim(filename), requestedLength = LINE_LENGTH)
         
         call STLfile_GetInfo( this% stl(STLNum), this% STLfilename(STLNum) )
         
         if( .not. this% HO_IBM ) this% stl(STLNum)% BFcorrection = .false. 

         this% bctype(STLNum) = this% stl(STLNum)% bctype
         
         if( MPI_Process% isRoot ) then
            this% stl(STLNum)% show = .true.
            call this% stl(STLNum)% ReadTessellation( this% STLfilename(STLNum) )
            if( this% ClipAxis .ne. 0 ) then
               call this% stl(STLNum)% Clip( this% minCOORDS, this% maxCOORDS, this% ClipAxis, describe )
            else
               if( describe ) call this% stl(STLNum)% describe(this% STLfilename(STLNum))
            end if
         end if
#ifdef _HAS_MPI_
         call this% MPI_sendSTLpartitions( STLNum, this% stl(STLNum)% maxAxis ) 
#endif
         Area = this% stl(STLNum)% ComputeArea()
#ifdef _HAS_MPI_
         localVal = Area
         call mpi_allreduce(localVal, Area, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
         this% stl(STLNum)% A = Area
         
         if( this% ComputeBandRegion .and. .not. this% HO_IBM ) then 
            call this% stl(STLNum)% SetIntegrationPoints()
            call this% stl(STLNum)% SetIntegration( this% NumOfInterPoints ) 
         end if 

         call OBB(STLNum)% construct( this% stl(STLNum), this% plotOBB )

         OBB(STLNum)% maxAxis = this% stl(STLNum)% maxAxis
         ! for  mpi 
         call OBB(STLNum)% GetGLobalVertices()

         OBB(STLNum)% minAxis = OBB(STLNum)% GetMinAxis( OBB(STLNum)% maxAxis, this% ClipAxis )

         if( .not. this% stl(STLNum)% BFcorrection ) call this% constructSTL_KDtree( STLNum )

         call this% stl(STLNum)% StoreInitialPosition( )

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

      vertices = OBB(STLNum)% Vertices

      this% root(STLNum)% MaxAxis = OBB(STLNum)% maxAxis
      
      call this% root(STLNum)% construct( stl           = this% stl(STLNum),         &
                                          vertices      = vertices,                  &
                                          isPlot        = this% plotKDtree,          &
                                          Min_n_of_Objs = this% KDtree_Min_n_of_Objs )
 
   end subroutine IBM_constructSTL_KDtree
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
 return
      allocate( this% root(parent% NumOfSTL),        &
                this% STLfilename(parent% NumOfSTL), &
                this% bctype(this% NumOfSTL),        &
                this% zoneMask(this% NumOfSTL),      &
                this% stl(this% NumOfSTL)            )

      this% ComputeBandRegion = parent% ComputeBandRegion
      this% HO_IBM            = parent% HO_IBM
      this% wallfunction      = parent% wallfunction
      this% ComputeDistance   = parent% ComputeDistance

      if( parent% ComputeBandRegion ) then
         allocate( this% BandRegion(parent% NumOfSTL) )
      end if

      do STLNum = 1, parent% NumOfSTL
         this% stl(STLNum)% BFcorrection = parent% stl(STLNum)% BFcorrection
         this% STLfilename(STLNum)       = parent% STLfilename(STLNum)
         this% bcType(STLNum)            = parent% bcType(STLNum)
         this% zoneMask(STLNum)          = parent% zoneMask(STLNum)
         this% root(STLNum)              = parent% root(STLNum)
      end do

      this% ClipAxis = parent% ClipAxis

      this% lvl = lvl

   end subroutine IBM_copy
!
!  Destroy the KD tree
!  --------------------
   subroutine IBM_DestroyKDtree( this, STLNum, isChild )
      use MPI_Process_Info
      implicit none
      !-arguments--------------------------------------------------
      class(IBM_type),           intent(inout) :: this
      integer,                   intent(in)    :: STLNum
      logical,                   intent(in)    :: isChild

      call this% root(STLNum)% destruct( isChild )

   end subroutine IBM_DestroyKDtree

   subroutine IBM_DestroyBandRegion( this, STLNum )
      use MPI_Process_Info
      implicit none
      !-arguments--------------------------------------------------
      class(IBM_type), intent(inout) :: this
      integer,         intent(in)    :: STLNum 

      integer :: domains 
      
      do domains = 1, MPI_Process% nProcs
         call this% BandRegion(STLNum)% IBMmask(domains)% destroy()
      end do

   end subroutine IBM_DestroyBandRegion


   subroutine IBM_constructzonemask( this, faces, elements, STLNum )
      use Meshtypes
      use ElementConnectivityDefinitions
      implicit none 

      class(IBM_type), intent(inout) :: this 
      type(face),      intent(inout) :: faces(:)
      type(element),   intent(inout)    :: elements(:)
      integer,         intent(in)    :: STLNum

      integer :: fID, i, j, N, eID

      if( .not. this% stl(STLNum)% BFcorrection ) return  
      
      this% NumOfMaskObjs = 0

      do fID = 1, size(faces)
         associate( f => faces(fID) )
         if( f% FaceType .eq. HMESH_INTERIOR .or. f% FaceType .eq. HMESH_MPI ) cycle
         if( trim(this% STLfilename(STLNum)) .eq. trim(f% boundaryName) ) then 
            if( f% HO_IBM ) cycle  
            f% HO_IBM  = .true.
            f% monitor = .true.
            allocate(f% stencil(0:f% Nf(1),0:f% Nf(2)))
            f% HOSIDE = minloc(f% elementIDs,dim=1)
            f% STLNum = STLNum
            do j = 0, f% Nf(2); do i = 0, f% Nf(1)
               f% stencil(i,j)% x  = f% geom% x(:,i,j)
               this% NumOfMaskObjs = this% NumOfMaskObjs + 1
            end do; end do
            eID = f% elementIDs(maxloc(f% elementIDs,dim=1))
            elements(eID)% HOcorrection = .true. 
         end if   
         end associate
      end do 

      this% zoneMask(STLNum) = .true.   

   end subroutine IBM_constructzonemask
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
      type(IBMPoints)      :: tocopy
      integer              :: eID, n, i, j, k, NumOfObjs, domains, domain
      integer, allocatable :: NumOfIntersections(:)
#ifdef _HAS_MPI_
      integer :: ierr
#endif
      if( this% zoneMask(STLNum) ) return 
      
      domain = MPI_Process% rank + 1

      allocate(this% IBMmask(MPI_Process% nProcs))

      call tocopy% build(no_of_DoFs)
      tocopy% NumOfObjs = 0

      do eID = 1, size(elements)
         associate( e => elements(eID) )
         call e% ConstructIBM( e% Nxyz(1), e% Nxyz(2), e% Nxyz(3), this% NumOfSTL, this% wallfunction )
         if( this% HO_IBM ) then
            do k = 1, NODES_PER_ELEMENT
               e% MaskCorners(k) = .false.
               if( OBB(STLNum)% isPointInside( Point = e% SurfInfo% corners(:,k) ) ) then
                  tocopy% NumOfObjs                           = tocopy% NumOfObjs + 1
                  tocopy% coords(tocopy% NumOfObjs,:)         = e% SurfInfo% corners(:,k)
                  tocopy% element_index(tocopy% NumOfObjs)    = eID
                  tocopy% local_position(tocopy% NumOfObjs,:) = (/k,0,0/)
               end if
            end do
         else
            do i = 0, e% Nxyz(1); do j = 0, e% Nxyz(2); do k = 0, e% Nxyz(3)
               if( e% STL(i,j,k) .eq. STLNum ) then 
                  e% STL(i,j,k)          = 0 
                  e% isInsideBody(i,j,k) = .false.
               end if
               if( OBB(STLNum)% isPointInside( Point = e% geom% x(:,i,j,k) ) ) then
                  tocopy% NumOfObjs                           = tocopy% NumOfObjs + 1
                  tocopy% coords(tocopy% NumOfObjs,:)         = e% geom% x(:,i,j,k)
                  tocopy% element_index(tocopy% NumOfObjs)    = eID
                  tocopy% local_position(tocopy% NumOfObjs,:) = (/i,j,k/)
                  tocopy% STLNum(tocopy% NumOfObjs)           = STLNum
               end if
            end do; end do; end do
         end if
         end associate
      end do

      call this% IBMmask(domain)% build( tocopy% NumOfObjs )
      call this% IBMmask(domain)% copy( tocopy             )
      call tocopy% destroy()
#ifdef _HAS_MPI_
      if( MPI_Process% doMPIAction ) then
         call castMaskNumOfObjs( this% IBMmask, domain )
         call castMask( this% IBMmask, domain )
      end if
#endif
      do domains = 1, MPI_Process% nProcs
         if( this% IBMmask(domains)% NumOfObjs .eq. 0 ) cycle
         do i = 1, this% IBMmask(domains)% NumOfObjs
            call this% CheckPoint( this% IBMmask(domains)% coords(i,:), STLNum, this% IBMmask(domains)% NumOfIntersections(i) )
         end do
      end do
#ifdef _HAS_MPI_
      do domains = 1, MPI_Process% nProcs
         NumOfObjs = this% IBMmask(domains)% NumOfObjs
         allocate(NumOfIntersections(NumOfObjs))
         NumOfIntersections = this% IBMmask(domains)% NumOfIntersections(1:NumOfObjs)
         call mpi_allreduce( NumOfIntersections, &
                           this% IBMmask(domains)% NumOfIntersections, &
                           this% IBMmask(domains)% NumOfObjs, MPI_INT, &
                           MPI_SUM, MPI_COMM_WORLD, ierr               )
         deallocate(NumOfIntersections)
      end do
#endif
      do i = 1, this% IBMmask(domain)% NumOfObjs
         if( mod(this% IBMmask(domain)% NumOfIntersections(i),2) .ne. 0 ) this% IBMmask(domain)% isInsideBody(i) = .true.
      end do
      
      this% NumOfMaskObjs = 0

      if( this% HO_IBM ) then
         do n = 1, this% IBMmask(domain)% NumOfObjs
            if( this% IBMmask(domain)% isInsideBody(n) ) then
               eID = this% IBMmask(domain)% element_index(n)
               k   = this% IBMmask(domain)% local_Position(n,1) 
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
         call this% HOmask( elements, faces, STLNum )
      else
         do n = 1, this% IBMmask(domain)% NumOfObjs
            if( this% IBMmask(domain)% isInsideBody(n) ) then
               eID = this% IBMmask(domain)% element_index(n)
               i   = this% IBMmask(domain)% local_Position(n,1)
               j   = this% IBMmask(domain)% local_Position(n,2)
               k   = this% IBMmask(domain)% local_Position(n,3)

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

      if( this% plot ) call this% plotMask( iter, STLNum, NumOfObjs )  

      do domains = 1, MPI_Process% nProcs
         call this% IBMmask(domains)% destroy()
      end do

      deallocate(this% IBMmask)
      
   end subroutine IBM_constructmask

   subroutine IBM_HOmask( this, elements, faces, STLNum )
      use Meshtypes
      use MPI_Process_Info
      implicit none

      class(IBM_type), intent(inout) :: this
      type(element),   intent(inout) :: elements(:)
      type(face),      intent(inout) :: faces(:)
      integer,         intent(in)    :: STLNum

      real(kind=RP)     :: Point(NDIM)
      integer           :: eID, fID, i, j, Nxi, Neta, NumOfPoints
#ifdef _HAS_MPI_
      integer :: ierr
#endif

      character(len=LINE_LENGTH) :: filename, it
      integer :: NumOfObjs, funit

      this% NumOfMaskObjs = 0

      do eID = 1, size(elements)
         associate( e => elements(eID) )
         if(all(e% MaskCorners)) then 
               e% HO_IBM = .true.
               cycle 
         end if 
         if(  all(e% MaskCorners((/4,3,7,8/))) .and. .not. ALL(e% MaskCorners((/1,2,6,5/))) ) then
            associate(f => faces(e% faceIDs(EBACK)))             
            call SetHOface( f, eID, STLNum, this% NumOfMaskObjs )      
            end associate
            e% HOcorrection = .true.
         end if
          if( all(e% MaskCorners((/1,2,6,5/))) .and. .not. all( e% MaskCorners((/4,3,7,8/)) ) ) then
            associate(f => faces(e% faceIDs(EFRONT)))                 
            call SetHOface( f, eID, STLNum, this% NumOfMaskObjs )      
            end associate
            e% HOcorrection = .true.
         end if
         if( all(e% MaskCorners((/5,6,7,8/))) .and. .not. all(e% MaskCorners((/1,2,3,4/))) ) then
            associate(f => faces(e% faceIDs(ETOP)))
            call SetHOface( f, eID, STLNum, this% NumOfMaskObjs )  
            end associate
            e% HOcorrection = .true.
         end if
         if( all(e% MaskCorners((/1,2,3,4/))) .and. .not. all(e% MaskCorners((/5,6,7,8/))) ) then
            associate(f => faces(e% faceIDs(EBOTTOM)))
            call SetHOface( f, eID, STLNum, this% NumOfMaskObjs )  
            end associate
            e% HOcorrection = .true.
         end if
         if( all(e% MaskCorners((/1,4,8,5/))) .and. .not. all(e% MaskCorners((/2,3,7,6/))) ) then
            associate(f => faces(e% faceIDs(ELEFT)))
            call SetHOface( f, eID, STLNum, this% NumOfMaskObjs )   
            end associate
            e% HOcorrection = .true.
         end if
         if( all(e% MaskCorners((/2,3,7,6/))) .and. .not. all(e% MaskCorners((/1,4,8,5/))) ) then
            associate(f => faces(e% faceIDs(ERIGHT)))
            call SetHOface( f, eID, STLNum, this% NumOfMaskObjs )  
            end associate
            e% HOcorrection = .true.
         end if

         end associate
      end do


            
      write(it,'(I10.10)') 0

      filename = trim(this% STLfilename(STLNum))//'_'//trim(adjustl(it))

      NumOfObjs = 0 
      do fID = 1, size(faces)
         if( .not. faces(fID)% HO_IBM ) NumOfObjs = NumOfObjs + (faces(fID)% Nf(1) + 1) + (faces(fID)% Nf(2) + 1) 
      end do

      call TecFileHeader( 'IBM/shiftedface_'//trim(filename), 'Mask Points', NumOfObjs, funit ) 

      do fID = 1, size(faces)
         if( .not. faces(fID)% HO_IBM ) cycle 
         do i = 0, faces(fID)% Nf(1); do j = 0, faces(fID)% Nf(2)
            write(funit,'(1X,1E15.8,1X,1E15.8,1X,1E15.8)') faces(fID)% geom% x(IX,i,j), faces(fID)% geom% x(IY,i,j), faces(fID)% geom% x(IZ,i,j)
         end do; end do 
      end do

      close(funit)

   end subroutine IBM_HOmask

   subroutine SetHOface( f, eID, STLNum, NumOfMaskObjs )
      use Meshtypes
      implicit none 

      type(face), intent(inout) :: f 
      integer,    intent(in)    :: eID, STLNum
      integer,    intent(inout) :: NumOfMaskObjs

      integer :: i, j 

      if( f% HO_IBM ) then
        f% HO_IBM = .false.
        NumOfMaskObjs = NumOfMaskObjs - (f% Nf(1)+1)*(f% Nf(2)+1)
        f% HOSIDE   = 0
        f% STLNum   = 0
      else
         f% HO_IBM  = .true.
         f% monitor = .true.
         allocate(f% stencil(0:f% Nf(1),0:f% Nf(2)))
         f% HOSIDE = FindHOside( f, eID )
         f% STLNum = STLNum
         do j = 0, f% Nf(2); do i = 0, f% Nf(1)
            f% stencil(i,j)% x  = f% geom% x(:,i,j)
            NumOfMaskObjs       = NumOfMaskObjs + 1
         end do; end do
         if( f% faceType .eq. HMESH_MPI ) f% HOSIDE = minloc(f% elementIDs, dim=1)
      end if

   end subroutine SetHOface

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
   subroutine IBM_plotMask(  this, iter, STLNum, NumOfObjs )
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
      if( MPI_Process% doMPIAction) then 
         call MaskLogical2Root(this% IBMmask, MPI_Process% rank + 1)
      end if
#endif
      if( .not. MPI_Process% isRoot ) return

      write(it,'(I10.10)') iter

      filename = trim(this% STLfilename(STLNum))//'_'//trim(adjustl(it))

      call TecFileHeader( 'IBM/Mask_'//trim(filename), 'Mask Points', NumOfObjs, funit ) 

      do domains = 1, MPI_Process% nProcs
         do i = 1, this% IBMmask(domains)% NumOfObjs
            if( this% IBMmask(domains)% isInsideBody(i) ) then
               write(funit,'(1X,1E15.8,1X,1E15.8,1X,1E15.8)') this% IBMmask(domains)% coords(i,IX), this% IBMmask(domains)% coords(i,IY), this% IBMmask(domains)% coords(i,IZ)
            end if
         end do
      end do

      close(funit)

   end subroutine IBM_plotMask
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
      real(kind=rp) :: epsilon = 1.0d-12
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
         do k = 1, NODES_PER_ELEMENT
            this% maxCOORDS = max(this% maxCOORDS,e% SurfInfo% corners(axis,k)); this% minCOORDS = min(this% minCOORDS,e% SurfInfo% corners(axis,k))
         end do
         end associate
      end do
      this% maxCOORDS = this% maxCOORDS + epsilon; this% minCOORDS = this% minCOORDS - epsilon 
#ifdef _HAS_MPI_
      localmax = this% maxCOORDS; localmin = this% minCOORDS
      call mpi_allreduce(localmax, this% maxCOORDS, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, ierr)
      call mpi_allreduce(localmin, this% minCOORDS, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, ierr)
#endif
   end subroutine IBM_GetDomainExtreme

   subroutine PlotSurface( IBM, STLNum, iter )
      use MPI_Process_Info
      implicit none

      type(IBM_type), intent(in) :: IBM
      integer,        intent(in) :: STLNum, iter

      character(len=LINE_LENGTH) :: filename
      integer                    :: domains, domain, NumOfObjs, funit, i, counter, N 

      if( .not. IBM% HO_IBM  .and. .not. MPI_Process% isRoot ) return

      domains   = MPI_Process% nProcs
      NumOfObjs = 0 

      do domain = 1, domains
         do i = 1, IBM% IBMStencilPoints(domain)% NumOfObjs
            if( IBM% IBMStencilPoints(domain)% STLNum(i) .eq. STLNum  .and. IBM% IBMStencilPoints(domain)% local_position(i,IZ) .eq. 0 ) NumOfObjs = NumOfObjs + 1 
         end do 
      end do 

      if( NumOfObjs .eq. 0 ) return

      if( IBM% lvl .gt. 0 ) then
         write(filename,'(A,A,I1,A,I10.10)') trim(IBM% STLfilename(STLNum)),'_MGlevel',IBM% lvl,'_',iter
      else
         write(filename,'(A)') trim(IBM% STLfilename(STLNum))
      end if

      call TecFileHeader( 'IBM/SurfacePoints_'//trim(filename), 'Mask Points', NumOfObjs, funit )

      do domain = 1, domains
         do i = 1, IBM% IBMStencilPoints(domain)% NumOfObjs
            if( IBM% IBMStencilPoints(domain)% STLNum(i) .eq. STLNum  .and. IBM% IBMStencilPoints(domain)% local_position(i,IZ) .eq. 0 ) &  
            write(funit,'(E17.10,1X,E17.10,1X,E17.10)')  IBM% IBMStencilPoints(domain)% coords(i,IX), IBM% IBMStencilPoints(domain)% coords(i,IY), IBM% IBMStencilPoints(domain)% coords(i,IZ)
         end do
      end do

      close(funit)

   end subroutine PlotSurface

   subroutine PlotStencil( IBM, STLNum, iter )
      use MPI_Process_Info
      implicit none

      type(IBM_type), intent(in) :: IBM
      integer,        intent(in) :: STLNum, iter

      character(len=LINE_LENGTH) :: filename
      integer                    :: i, NumOfObjs, funit, domains, domain

      if( .not. IBM% HO_IBM  .and. .not. MPI_Process% isRoot ) return

      domains   = MPI_Process% nProcs
      NumOfObjs = 0
      
      do domain = 1, domains
         do i = 1, IBM% IBMStencilPoints(domain)% NumOfObjs
            if( IBM% IBMStencilPoints(domain)% STLNum(i) .eq. STLNum ) NumOfObjs = NumOfObjs + 1
         end do 
      end do 
   
      if( IBM% lvl .gt. 0 ) then
         write(filename,'(A,A,I1,A,I10.10)') trim(IBM% STLfilename(STLNum)),'_MGlevel',IBM% lvl,'_',iter
      else
         write(filename,'(A)') trim(IBM% STLfilename(STLNum))
      end if

      call TecFileHeader( 'IBM/StencilPoints_'//trim(filename), 'Mask Points', NumOfObjs, funit )
      
      do domain = 1, domains
         do i = 1, IBM% IBMStencilPoints(domain)% NumOfObjs
            if( IBM% IBMStencilPoints(domain)% STLNum(i) .ne. STLNum ) cycle 
            write(funit,'(E17.10,1X,E17.10,1X,E17.10)')  IBM% IBMStencilPoints(domain)% coords(i,IX), IBM% IBMStencilPoints(domain)% coords(i,IY), IBM% IBMStencilPoints(domain)% coords(i,IZ)
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
   subroutine IBM_build( this, elements, faces, MPIfaces, no_of_DoFs, STLNum, isChild, movingSTL, iter )
      use MPI_Process_Info
      use PhysicsStorage
      use MPI_Face_Class
      implicit none
      !-arguments-----------------------------------------------------------------
      class(IBM_type),      intent(inout) :: this
      type(element),        intent(inout) :: elements(:)
      type(face),           intent(inout) :: faces(:)
      type(MPI_FacesSet_t), intent(inout) :: MPIfaces
      integer,              intent(in)    :: no_of_DoFs, STLNum 
      logical,              intent(in)    :: isChild, movingSTL
      integer,              intent(in)    :: iter

      call this% constructzonemask( faces, elements, STLNum )

      call this% constructmask( elements, STLNum, no_of_DoFs, faces, iter )
#ifdef _HAS_MPI_
      if( this% HO_IBM ) call FixingmpiFaces( faces, MPIfaces, this% NumOfMaskObjs, this% MPIfixed )
#endif 
      if( this% ComputeDistance ) call this% ComputeIBMWallDistance( elements=elements, STLNum=STLNum, faces=faces )

      if( .not. movingSTL ) then 
         if( this% ComputeBandRegion .and. .not. this% HO_IBM ) call this% constructBandRegion( elements, faces, no_of_DoFs, STLNum, NCONS )
         if( this% wallfunction  .and. .not. this% HO_IBM ) then 
            call this% GetForcingPointsGeom( elements )
            call this% constructImagePoint( elements, no_of_DoFs, STLNum, NCONS )
         end if 
      end if 

      call this% DestroyKDtree( STLNum, isChild )

   end subroutine IBM_build

   subroutine IBM_buildHOfaces( this, elements, faces )
      use PhysicsStorage
      use MPI_Process_Info
      implicit none

      class(IBM_type),              intent(inout) :: this
      type(element),                intent(inout) :: elements(:)
      type(face),                   intent(inout) :: faces(:)

      integer :: STLNum
      
      call this% GetStencil( elements, faces )
      
      allocate(this% IBMStencilPoints(MPI_Process% nProcs))
      
      call Set_IBM_HO_faces( this% IBMStencilPoints, faces, NCONS )
#ifdef _HAS_MPI_
      call MPIProcedures_IBM_HO_faces( this% IBMStencilPoints, NCONS )
#endif
      do STLNum = 1, this% NumOfSTL
         call PlotSurface( this, STLNum, this% iter )
         call PlotStencil( this, STLNum, this% iter )
      end do

      call IBM_HO_findElements( this% IBMStencilPoints, elements, this% NumOfSTL, this% clipAxis )

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

      integer :: STLNum 

      if ( .not. MPI_Process % isRoot ) return
      write(STD_OUT,'(/)')
      call Section_Header("IBM parameters")
      write(STD_OUT,'(/)')

      call SubSection_Header('IBM info')
      if( .not. this% HO_IBM ) then
         write(STD_OUT,'(30X,A,A35,A20)') "->" ,"Type: ","volume penalization"
         write(STD_OUT,'(30X,A,A35,L10)') "->" , "Semi-implicit treatment: ", this% semiImplicit
         if( .not. this% TimePenal ) then
            write(STD_OUT,'(30X,A,A35,ES14.2)') "->" , "Penalization term: " , this% eta
         else
            write(STD_OUT,'(30X,A,A35,A10)') "->" , "Penalization term: ", " Dt"
         end if
      else
         write(STD_OUT,'(30X,A,A35,A10)') "->" ,"Type: ","high-order"
      end if

      if( this% ComputeBandRegion .and. .not. this% HO_IBM ) then 
         write(STD_OUT,'(30X,A,A35,I10)') "->" , "Minimum number of objects: ", this% KDtree_Min_n_of_Objs
         write(STD_OUT,'(30X,A,A35,I10)') "->" , "Number of interpolation points: ", this% NumOfInterPoints
      end if 
      if( this% Wallfunction ) write(STD_OUT,'(30X,A,A35,A10)') "->" , "Wall function: "," active"

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
      end do

      if( this% ComputeBandRegion ) deallocate( this% BandRegion )

      if( this% Wallfunction .and. .not. this% HO_IBM ) then
         deallocate( this% ImagePoints )
      end if

      deallocate( this% penalization, &
                  this% stl,          &
                  this% Integral,     &
                  this% STLfilename,  &
                  this% bctype,       &
                  this% zoneMask      )

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
            do STLNum = 1, this% NumOfSTL
               do k = 1, BOXVERTICES
                  if( OBB(STLNum)% isPointInside( e% SurfInfo% corners(:,k), 1.5_RP ) ) then
                     e% Nxyz(1) = this% Nx
                     e% Nxyz(2) = this% Ny
                     e% Nxyz(3) = this% Nz
                     exit
                  end if
               end do
            end do 
            end associate
         end do
      else
         do eID = 1, size(elements)
            associate(e => elements(eID))
            do STLNum = 1, this% NumOfSTL
               loop: do k = 1, BOXVERTICES
                  if( OBB(STLNum)% isPointInside( e% SurfInfo% corners(:,k), 1.5_RP ) ) then
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

   subroutine IBM_constructBandRegion( this, elements, faces, no_of_DoFs, STLNum, nEqn )
      use MPI_Process_Info
      use PhysicsStorage
      use NodalStorageClass
      use PolynomialInterpAndDerivsModule
      implicit none
      !-arguments------------------------------------------------------
      class(IBM_type), intent(inout) :: this
      type(element),   intent(inout) :: elements(:)
      type(face),      intent(in)    :: faces(:)
      integer,         intent(in)    :: no_of_DoFs, STLNum, nEqn
      !-local-variables------------------------------------------------
      type(IBMPoints)               :: tocopy
      integer                       :: eID, i, j, k, n, domains, domain, NumOfObjs, indeces(4)
      type(NodalStorage_t), pointer :: spA 
      real(kind=RP)                 :: ds, h, L, coords(NDIM)
      real(kind=RP), allocatable    :: nodes(:)
      
      if( this% HO_IBM ) return 

      if( .not. allocated(this% BandRegion(STLNum)% IBMmask) ) allocate(this% BandRegion(STLNum)% IBMmask(MPI_Process% nProcs))

      domain = MPI_Process% rank + 1

      call tocopy% build(no_of_DoFs)
      tocopy% NumOfObjs = 0
      
      do eID = 1, size(elements)
         associate(e => elements(eID))
         do i = 0, e% Nxyz(1); do j = 0, e% Nxyz(2); do k = 0, e% Nxyz(3)
            if( e% isInsideBody(i,j,k) ) cycle 
            if( OBB(STLNum)% isPointInside( Point = e% geom% x(:,i,j,k), coeff = this% BandRegionCoeff ) ) then
               tocopy% NumOfObjs                           = tocopy% NumOfObjs + 1
               tocopy% coords(tocopy% NumOfObjs,:)         = e% geom% x(:,i,j,k)
               tocopy% element_index(tocopy% NumOfObjs)    = eID
               tocopy% local_position(tocopy% NumOfObjs,:) = (/i,j,k/)
            end if
         end do; end do; end do
         end associate
      end do
      
      call this% BandRegion(STLNum)% IBMmask(domain)% build( tocopy% NumOfObjs )
      call this% BandRegion(STLNum)% IBMmask(domain)% copy( tocopy )
      call tocopy% destroy()
#ifdef _HAS_MPI_
      if( MPI_Process% doMPIAction ) then
         call castMaskNumOfObjs( this% BandRegion(STLNum)% IBMmask, domain )
         call castMask( this% BandRegion(STLNum)% IBMmask, domain )
      end if
#endif
      do domains = 1, MPI_Process% nProcs
         call this% BandRegion(STLNum)% IBMmask(domains)% buildState( this% BandRegion(STLNum)% IBMmask(domains)% NumOfObjs, nEqn )
      end do

      if( this% plotBandPoints .and. MPI_Process% isRoot ) call this% BandRegion(STLNum)% plot( 'BandPoints' )
      
   end subroutine IBM_constructBandRegion 

   subroutine IBM_minDistance( this, Point, IBMmask, indeces, domain )
      use MPI_Process_Info
      implicit none
      !-arguments-----------------------------------------------
      class(IBM_type), intent(inout) :: this
      type(IBMPoints), intent(in)    :: IBMmask(:) 
      real(kind=RP),   intent(in)    :: Point(NDIM)
      integer,         intent(inout) :: indeces(this% NumOfInterPoints), &
                                        domain(this% NumOfInterPoints)
      !-local-variables-----------------------------------------
      real(kind=RP) :: dist, dist_, lastDist
      integer       :: k, domains, i

      lastDist = -huge(1.0_RP)
      do k = 1, this% NumOfInterPoints
         dist = huge(1.0_RP)
         do domains = 1, MPI_Process% nProcs
            do i = 1, IBMmask(domains)% NumOfObjs
               dist_ = norm2(Point - IBMmask(domains)% coords(i,:))
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

   subroutine IBM_GetPointInterpolation( this, point, IBMmask )

      implicit none 

      class(IBM_type),  intent(inout) :: this 
      type(point_type), intent(inout) :: point 
      type(IBMPoints),  intent(in)    :: IBMmask(:)

      call this% minDistance( point% coords, IBMmask, point% indeces, point% domains )
      
      call GetMatrixInterpolationSystem( point, this% NumOfInterPoints, IBMmask, this% InterpolationType )

   end subroutine IBM_GetPointInterpolation
!
!  Band region points' state is stored
!  -----------------------------------
   subroutine BandPointsState( this, elements, nEqn, HO_IBM, zoneMask, stl )
      use PhysicsStorage
      use MPI_Process_Info
      implicit none
      !-arguments------------------------------------------------
      type(bandRegion_t), intent(inout) :: this
      type(element),      intent(inout) :: elements(:)
      integer,            intent(in)    :: nEqn
      logical,            intent(in)    :: HO_IBM, zoneMask
      type(STLfile),      intent(inout) :: stl
      !-local-variables------------------------------------------
      integer                    :: n, i, j, k, eID, domain

      domain = MPI_Process% rank + 1

      do n = 1, this% IBMmask(domain)% NumOfObjs
         i   = this% IBMmask(domain)% local_Position(n,1)
         j   = this% IBMmask(domain)% local_Position(n,2)
         k   = this% IBMmask(domain)% local_Position(n,3)
         eID = this% IBMmask(domain)% element_index(n)

         this% IBMmask(domain)% Q(n,:)   = elements(eID)% storage% Q(:,i,j,k)
         this% IBMmask(domain)% U_x(n,:) = elements(eID)% storage% U_x(:,i,j,k)
         this% IBMmask(domain)% U_y(n,:) = elements(eID)% storage% U_y(:,i,j,k)
         this% IBMmask(domain)% U_z(n,:) = elements(eID)% storage% U_z(:,i,j,k)
      end do
#ifdef _HAS_MPI_
      if( MPI_Process% doMPIAction ) then
         call castStateBandRegion( this% IBMmask, nEqn )
         call castGradientsBandRegion( this% IBMmask, nEqn )
      endif
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

      call SendAxis( this% stl(STLNum) ) 
      call RecvAxis( this% stl(STLNum) )

   end subroutine IBM_MPI_sendSTLpartitions
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------------------------
!  Moving bodies
!  ------------------------------------------------------------
   subroutine IBM_MoveBody( this, elements, faces, MPIfaces, no_of_DoFs, isChild, t, iter, autosave )
      use MPI_Process_Info
      use FluidData
      use MPI_Face_Class
      implicit none
      !-arguments-----------------------------------------------------------------
      class(IBM_type),           intent(inout) :: this
      type(element),             intent(inout) :: elements(:)
      type(face),                intent(inout) :: faces(:)
      type(MPI_FacesSet_t),      intent(inout) :: MPIfaces
      integer,                   intent(in)    :: no_of_DoFs
      logical,                   intent(in)    :: isChild
      real(kind=RP),             intent(in)    :: t
      integer,         optional, intent(in)    :: iter
      logical,         optional, intent(in)    :: autosave
      !-local-variables-----------------------------------------------------------
      integer       :: i, domains, NumOfObjs, j, STLNum
      real(kind=RP) :: cL, cD 

      this% dt = t - this% t 
 
      if( any(this% stl(:)% move) ) this% MPIfixed = .false.
      
      do STLNum = 1, this% NumOfSTL

         if( .not. this% stl(STLNum)% move ) cycle 
         
         do i = 1, this% stl(STLNum)% NumOfObjs
            do j = 1, NumOfVertices
#if defined(NAVIERSTOKES)
               call BCsIBM(STLNum)% bc% PositionMoving_IBM( this% stl(STLNum)% ObjectsList(i)% vertices(j)% coords, t, this% dt, cL, cD )
#endif 
            end do 
            call this% stl(STLNum)% updateNormals( this% stl(STLNum)% ObjectsList(i) )
         end do 

         if( this% ComputeBandRegion ) call this% stl(STLNum)% SetIntegrationPoints()

         call OBB(STLNum)% construct( this% stl(STLNum), this% plotOBB )

         if( autosave ) call plotSTL( this% stl(STLNum), iter )

         call OBB(STLNum)% GetGLobalVertices()
         
         this% plot              = autosave
         this% plotKDtree        = .false.
         this% stl(STLNum)% show = .false.
         
         call this% constructSTL_KDtree( STLNum )
         
         call this% build( elements, faces, MPIfaces, no_of_DoFs, STLNum, isChild, .true., iter )

         call this% stl(STLNum)% ResetInitialPosition() 
         
      end do 

      this% t = t 

   end subroutine IBM_MoveBody

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
      real(kind=RP) :: h, d, d_min
      integer       :: eID, i, j, k, NumOfObjs
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
      call mpi_allreduce(d_min, this% IP_dWall, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, ierr)
#else
      this% IP_dWall = d_min
#endif
   end subroutine IBM_GetForcingPointsGeom

   subroutine IBM_constructImagePoint( this, elements, no_of_DoFs, STLNum, nEqn )
      use MPI_Process_Info
      implicit none

      class(IBM_type), intent(inout) :: this
      type(element),   intent(inout) :: elements(:)
      integer,         intent(in)    :: no_of_DoFs, STLNum, nEqn
      !-local-variables------------------------------------------------
      type(IBMPoints) :: tocopy
      integer         :: eID, i, j, k, n

      if( .not. allocated(this% ImagePoints(STLNum)% IBMmask) ) allocate(this% ImagePoints(STLNum)% IBMmask(1))

      call tocopy% buildPoints(no_of_DoFs)
      tocopy% NumOfObjs = 0

      do eID = 1, size(elements)
         associate(e => elements(eID))
         if( any(e% isInsideBody) ) cycle
         do i = 0, e% Nxyz(1); do j = 0, e% Nxyz(2); do k = 0, e% Nxyz(3)
            if( e% geom% dWall(i,j,k) .lt. this% IP_dWall ) then 
               tocopy% NumOfObjs                            = tocopy% NumOfObjs + 1 
               tocopy% x(tocopy% NumOfObjs)% element_index  = eID  
               tocopy% x(tocopy% NumOfObjs)% local_position = (/i,j,k/)  
               tocopy% x(tocopy% NumOfObjs)% STLNum         = STLNum  
               e% isForcingPoint(i,j,k)                     = .true.
            end if
         end do; end do; end do
         end associate
      end do
      
      call this% ImagePoints(STLNum)% IBMmask(1)% buildPoints(tocopy% NumOfObjs)
      call this% ImagePoints(STLNum)% IBMmask(1)% copyPoints(tocopy)
      call tocopy% destroyPoints()

      do n = 1, this% ImagePoints(STLNum)% IBMmask(1)% NumOfObjs
         eID = this% ImagePoints(STLNum)% IBMmask(1)% x(n)% element_index
         i   = this% ImagePoints(STLNum)% IBMmask(1)% x(n)% local_position(IX)
         j   = this% ImagePoints(STLNum)% IBMmask(1)% x(n)% local_position(IY)
         k   = this% ImagePoints(STLNum)% IBMmask(1)% x(n)% local_position(IZ)

         this% ImagePoints(STLNum)% IBMmask(1)% x(n)% coords = elements(eID)% geom% x(:,i,j,k) + &
                                                               this% IP_dWall * elements(eID)% geom% normal(:,i,j,k)

         elements(eID)% forcingPointIndex(i,j,k) = n
         allocate( this% ImagePoints(STLNum)% IBMmask(1)% x(n)% invPhi(this% NumOfInterPoints,this% NumOfInterPoints), &
                   this% ImagePoints(STLNum)% IBMmask(1)% x(n)% domains(this% NumOfInterPoints),                       &
                   this% ImagePoints(STLNum)% IBMmask(1)% x(n)% b(this% NumOfInterPoints)                              )
      end do 

   end subroutine IBM_constructImagePoint

   subroutine IBM_updateImagePoint( this, nEqn, elements, no_of_DoFs )

      implicit none 

      class(IBM_type), intent(inout) :: this 
      integer,         intent(in)    :: nEqn, no_of_DoFs
      type(element),   intent(inout) :: elements(:)

      integer :: STLNum , i

      if( .not. this% wallfunction .or. this% HO_IBM ) return 

      do STLNum = 1, this% NumOfSTL  
         call this% DestroyBandRegion( STLNum )
#ifdef FLOW
    !     call this% constructBandRegion( elements, no_of_DoFs, STLNum, nEqn )
#endif
         do i = 1, this% ImagePoints(STLNum)% IBMmask(1)% NumOfObjs
            call this% GetPointInterpolation( this% ImagePoints(STLNum)% IBMmask(1)% x(i), this% BandRegion(STLnum)% IBMmask )
            call GetPointState( nEqn, this% ImagePoints(STLNum)% IBMmask(1)% x(i), this% BandRegion(STLnum)% IBMmask, this% NumOfInterPoints, &
                                this% InterpolationType, this% ImagePoints(STLNum)% IBMmask(1)% x(i)% Q                                       )
         end do 
      end do 

   end subroutine IBM_updateImagePoint

   subroutine GetPointState( nEqn, Point, IBMmask, NumOfInterPoints, InterpolationType, Q )
      use PhysicsStorage
      use MPI_Process_Info
      implicit none
      !-arguments------------------------------------------------
      type(point_type), intent(in)    :: Point
      type(IBMpoints),  intent(in)    :: IBMmask(:)
      integer,          intent(in)    :: nEqn, NumOfInterPoints, &
                                         InterpolationType
      real(kind=RP),    intent(inout) :: Q(nEqn)

      real(kind=RP) :: Q_(nEqn,NumOfInterPoints)
      integer       :: k, n

      do k = 1, NumOfInterPoints
         Q_(:,k) = IBMmask(point% domains(k))% Q(point% indeces(k),:)
      end do 
      do n = 1, nEqn 
         Q(n) = GetInterpolatedValue( Q_(n,:), Point% invPhi, Point% b, InterpolationType )
      end do 
      
    end subroutine GetPointState

    subroutine GetPointGrads( nEqn, Point, IBMmask, NumOfInterPoints, InterpolationType, U_x, U_y, U_z )
      use PhysicsStorage
      use MPI_Process_Info
      implicit none
      !-arguments------------------------------------------------
      type(point_type), intent(in)    :: Point
      type(IBMpoints),  intent(in)    :: IBMmask(:)
      integer,          intent(in)    :: nEqn, NumOfInterPoints, &
                                         InterpolationType
      real(kind=RP),    intent(inout) :: U_x(nEqn), U_y(nEqn), U_z(nEqn)

      real(kind=RP) :: U_x_(nEqn,NumOfInterPoints), U_y_(nEqn,NumOfInterPoints), U_z_(nEqn,NumOfInterPoints)
      integer       :: k, n

      do k = 1, NumOfInterPoints
         U_x_(:,k) = IBMmask(Point% domains(k))% U_x(Point% indeces(k),:)
         U_y_(:,k) = IBMmask(Point% domains(k))% U_y(Point% indeces(k),:)
         U_z_(:,k) = IBMmask(Point% domains(k))% U_z(Point% indeces(k),:)
      end do 

      do n = 1, nEqn 
         U_x(n) = GetInterpolatedValue( U_x_(n,:), Point% invPhi, Point% b, InterpolationType )
         U_y(n) = GetInterpolatedValue( U_y_(n,:), Point% invPhi, Point% b, InterpolationType )
         U_z(n) = GetInterpolatedValue( U_z_(n,:), Point% invPhi, Point% b, InterpolationType )
      end do 

    end subroutine GetPointGrads
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
      real(kind=RP) :: Point(NDIM), normal(NDIM), dist, s, t, s_, t_ 
      real(kind=RP) :: IntersectionPoint(NDIM), P(NDIM), xP(NDIM)
      integer       :: domain, domains, i, j, k, n, fID, eID, index
#ifdef _HAS_MPI_
      integer       :: ierr
#endif
      domain = MPI_Process% rank + 1
      allocate(this% IBMmask(MPI_Process% nProcs))
      call this% IBMmask(domain)% build(this% NumOfMaskObjs)

      this% IBMmask(domain)% NumOfObjs = 0

      if( this% HO_IBM ) then
         do fID = 1, size(faces)
            associate(f => faces(fID))
            if( f% HO_IBM ) then
               if( f% STLNum .ne. STLNum ) cycle 
               do i = 0, f% Nf(1); do j = 0, f% Nf(2)
                  this% IBMmask(domain)% NumOfObjs = this% IBMmask(domain)% NumOfObjs + 1
                  this% IBMmask(domain)% coords(this% IBMmask(domain)% NumOfObjs,:)         = f% geom% x(:,i,j)
                  this% IBMmask(domain)% element_index(this% IBMmask(domain)% NumOfObjs)    = fID
                  this% IBMmask(domain)% local_position(this% IBMmask(domain)% NumOfObjs,:) = (/i,j,0/)
               end do; end do
            end if
            end associate
         end do
      else
         do eID = 1, size(elements)
            associate(e => elements(eID))
            do i = 0, e% Nxyz(1); do j = 0, e% Nxyz(2); do k = 0, e% Nxyz(3)
               if( e% isInsideBody(i,j,k) ) cycle
               if( OBB(STLNum)% isPointInside( Point = e% geom% x(:,i,j,k), coeff=this% BandRegionCoeff ) ) then
                  e% STL(i,j,k) = STLNum
                  this% IBMmask(domain)% NumOfObjs = this% IBMmask(domain)% NumOfObjs + 1
                  this% IBMmask(domain)% coords(this% IBMmask(domain)% NumOfObjs,:)         = e% geom% x(:,i,j,k)
                  this% IBMmask(domain)% element_index(this% IBMmask(domain)% NumOfObjs)    = eID
                  this% IBMmask(domain)% local_position(this% IBMmask(domain)% NumOfObjs,:) = (/i,j,k/)
               end if
            end do; end do; end do
            end associate
         end do
      end if
#ifdef _HAS_MPI_
      if( MPI_Process% doMPIAction ) then
         call castMaskNumOfObjs( this% IBMmask, domain )
         call castMask( this% IBMmask, domain )
      end if
#endif
      do domains = 1, MPI_Process% nProcs
         allocate( this% IBMmask(domains)% dist(this% IBMmask(domains)% NumOfObjs),       &
                   this% IBMmask(domains)% normal(this% IBMmask(domains)% NumOfObjs,NDIM) )
      end do 

      do domains = 1, MPI_Process% nProcs
         do n = 1, this% IBMmask(domains)% NumOfObjs
            this% IBMmask(domains)% dist(n) = HUGE(1.0_RP)
            Point = this% IBMmask(domains)% coords(n,:)
            do i = 1, this% stl(STLNum)% NumOfObjs
               call MinimumPointTriDistance( Point, this% stl(STLNum)% ObjectsList(i)% vertices(1) % coords, &
                                                    this% stl(STLNum)% ObjectsList(i)% vertices(2) % coords, &
                                                    this% stl(STLNum)% ObjectsList(i)% vertices(3) % coords, &
                                                    dist, IntersectionPoint , s, t                           )
               if( Dist .lt. this% IBMmask(domains)% dist(n) ) then
                  this% IBMmask(domains)% dist(n) = dist
                  normal = (Point - IntersectionPoint)/norm2( (Point - IntersectionPoint) )
                  index  = i 
               end if 
            end do
            this% IBMmask(domains)% normal(n,:) = normal 
            ! if(  Point(IX) .gt. 0.0136_RP .and. Point(IX) .lt. 0.0138_RP .and. Point(IY) .gt. 0.0154_RP .and. Point(IY) .lt. 0.0156_RP ) then 
            !    write(*,*) 'coords =', Point 
            !    write(*,*) 'dot product =',sign(1.0_RP,dot_product(normal, this% stl(STLNum)% ObjectsList(index)% normal))
            !    this% IBMmask(domains)% dist(n)     = -this% IBMmask(domains)% dist(n)
            ! end if
            !GENERAL
            if( this% zoneMask(STLNum) ) then 
                if( sign(1.0_RP,dot_product(normal, this% stl(STLNum)% ObjectsList(index)% normal)) .gt. 0.0_RP ) then 
                   this% IBMmask(domains)% normal(n,:) = -normal
                   !this% IBMmask(domains)% normal(n,:) = normal
                elseif( sign(1.0_RP,dot_product(normal, this% stl(STLNum)% ObjectsList(index)% normal)) .lt. 0.0_RP ) then 
                  this% IBMmask(domains)% dist(n)     = -this% IBMmask(domains)% dist(n)
                  !this% IBMmask(domains)% dist(n)     = this% IBMmask(domains)% dist(n)
               end if
            else 
               this% IBMmask(domains)% normal(n,:) = -normal
            end if 
            if( AlmostEqual(this% IBMmask(domains)% dist(n),0.0_RP) ) this% IBMmask(domains)% normal(n,:) = this% stl(STLNum)% ObjectsList(index)% normal
         end do
      end do
#ifdef _HAS_MPI_
      call gatherMaskGeom( this% IBMmask )
#endif
      if( this% HO_IBM ) then
         do n = 1, this% IBMmask(domain)% NumOfObjs
            fID = this% IBMmask(domain)% element_index(n)
            i   = this% IBMmask(domain)% local_Position(n,1)
            j   = this% IBMmask(domain)% local_Position(n,2)
            
            faces(fID)% stencil(i,j)% dist   = this% IBMmask(domain)% dist(n) 
            faces(fID)% stencil(i,j)% normal = this% IBMmask(domain)% normal(n,:)
         end do
      else
         do n = 1, this% IBMmask(domain)% NumOfObjs
            eID = this% IBMmask(domain)% element_index(n)
            i   = this% IBMmask(domain)% local_Position(n,1)
            j   = this% IBMmask(domain)% local_Position(n,1)
            k   = this% IBMmask(domain)% local_Position(n,1)

            elements(eID)% geom% dwall(i,j,k)    = this% IBMmask(domain)% dist(n)
            elements(eID)% geom% normal(:,i,j,k) = this% IBMmask(domain)% normal(n,:)
         end do
      end if

      do domains = 1, MPI_Process% nProcs
         call this% IBMmask(domains)% destroy()
         if( allocated(this% IBMmask(domains)% dist  ) ) deallocate(this% IBMmask(domains)% dist  )
         if( allocated(this% IBMmask(domains)% normal) ) deallocate(this% IBMmask(domains)% normal)
      end do
      
      deallocate(this% IBMmask)

   end subroutine IBM_ComputeIBMWallDistance

   subroutine IBM_GetStencil( this, elements, faces )
      use NodalStorageClass
      use MPI_Process_Info
      implicit none

      class(IBM_type), intent(inout) :: this
      type(element),   intent(in)    :: elements(:)
      type(face),      intent(inout) :: faces(:)

      real(kind=RP)                 :: alpha, L, dist, dl, h, d_, dl_, L_, sign_, x_mid(NDIM)
      real(kind=RP)                 :: x0(NDIM), normal(NDIM), xi(NDIM), Center(NDIM), R, x, y 
      integer                       :: fID, eID, N, M, i, j, k, NumOfObjs_, N_, NumOfObjs, info
      type(NodalStorage_t), pointer :: spA
      
      do fID = 1, size(faces)
         associate( f=> faces(fID) )
         if( f% HO_IBM ) then
            N = max(f% Nf(1),f% Nf(2)) 
            M = N - 1

            spA => NodalStorage(N)

            call spA_s% construct(GAUSS, M) 

            eID = f% elementIDs(maxloc(f% elementIDs, dim=1))
            h   = hGeom(150._RP) 
            L   = hgeom(10._RP)

            do i = 0, f% Nf(1); do j = 0, f% Nf(2)
               
               dl                      = abs(f% stencil(i,j)% dist) + h
               f% stencil(i,j)% N      = N
               f% stencil(i,j)% d      = f% stencil(i,j)% dist
               f% stencil(i,j)% L      = L
               f% stencil(i,j)% dl     = dl
               f% stencil(i,j)% eID    = f% elementIDs(f% HOSIDE)
               f% stencil(i,j)% domain = MPI_Process% rank + 1

               call f% stencil(i,j)% build( M )
               f% stencil(i,j)% x_s(:,0) = f% stencil(i,j)% x + f% stencil(i,j)% dist * f% stencil(i,j)% normal

               if( this% Wallfunction ) f% stencil(i,j)% wallfunction = .true. 

            end do; end do
         end if
         end associate
      end do

   end subroutine IBM_GetStencil

   real(kind=RP) function hGeom( y_plus )
      use FluidData
      implicit none 

      real(kind=RP), intent(in) :: y_plus

      real(kind=RP) :: Cf, Tau_w, u_star
#if defined(NAVIERSTOKES)
      associate( Re => Dimensionless% Re ) 
 
      Cf     = (2.0_RP * log10(Re) - 0.65_RP)**(-2.3_RP)
      Tau_w  = 0.5_RP * Cf * refValues% V**2
      u_star = sqrt(Tau_w)

      hGeom = y_plus * refValues% mu/(refValues% rho * u_star)
      
      end associate
#else 
      hGeom = 0.0_RP 
#endif 
   end function hGeom

   subroutine ShiftingFirstPointWF( faces )

      implicit none 

      type(face), intent(inout) :: faces(:)

      integer :: fID, i, j 

      do fID = 1, size(faces)
         associate( f=> faces(fID) )
         if( f% HO_IBM ) then
            do i = 0, f% Nf(1); do j = 0, f% Nf(2)
               f% stencil(i,j)% x_s(:,0) = 0.5_RP*(f% stencil(i,j)% x_s(:,0)+f% stencil(i,j)% x_s(:,1))
               f% stencil(i,j)% dWall(0) = 0.5_RP * f% stencil(i,j)% dWall(1)
            end do; end do
         end if
         end associate
      end do

   end subroutine ShiftingFirstPointWF


   subroutine IBM_HO_IBMstencilState( this, nEqn, elements, faces )
      use PhysicsStorage
      implicit none

      class(IBM_type), intent(inout) :: this
      integer,         intent(in)    :: nEqn
      type(element),   intent(inout) :: elements(:)
      type(face),      intent(inout) :: faces(:)

      real(kind=RP) :: xi(NDIM), Qs(nEqn)
      integer       :: fID, i, j, k, eID, domain, ierr, n, m

      call IBM_HO_GetState( this% IBMStencilPoints, elements, nEqn )
#ifdef _HAS_MPI_
      call GatherHOfacesState( this% IBMStencilPoints, nEqn, faces )
#else
      domain = MPI_Process% rank + 1
      do n = 1, this% IBMStencilPoints(domain)% NumOfObjs
         fID = this% IBMStencilPoints(domain)% fIDs(n)
         i   = this% IBMStencilPoints(domain)% local_position(n,IX)
         j   = this% IBMStencilPoints(domain)% local_position(n,IY)
         k   = this% IBMStencilPoints(domain)% local_position(n,IZ)

         faces(fID)% stencil(i,j)% Q(:,k) = this% IBMStencilPoints(domain)% Q(n,:)
      end do
#endif
   end subroutine IBM_HO_IBMstencilState

   subroutine IBM_HO_IBMstencilGradient( this, nEqn, elements, faces )
      use PhysicsStorage
      implicit none

      class(IBM_type), intent(inout) :: this
      integer,         intent(in)    :: nEqn
      type(element),   intent(inout) :: elements(:)
      type(face),      intent(inout) :: faces(:)

      real(kind=RP) :: xi(NDIM), Qs(nEqn)
      integer       :: fID, i, j, k, eID, domain, ierr, n, m

      call IBM_HO_GetGradient( this% IBMStencilPoints, elements, nEqn )
#ifdef _HAS_MPI_
      call GatherHOfacesGradient_x( this% IBMStencilPoints, nEqn, faces )
      call GatherHOfacesGradient_y( this% IBMStencilPoints, nEqn, faces )
      call GatherHOfacesGradient_z( this% IBMStencilPoints, nEqn, faces )
#else
      domain = MPI_Process% rank + 1
      do n = 1, this% IBMStencilPoints(domain)% NumOfObjs
         fID = this% IBMStencilPoints(domain)% fIDs(n)
         i   = this% IBMStencilPoints(domain)% local_position(n,IX)
         j   = this% IBMStencilPoints(domain)% local_position(n,IY)
         k   = this% IBMStencilPoints(domain)% local_position(n,IZ)

         faces(fID)% stencil(i,j)% U_x(:,k) = this% IBMStencilPoints(domain)% U_x(n,:)
         faces(fID)% stencil(i,j)% U_y(:,k) = this% IBMStencilPoints(domain)% U_y(n,:)
         faces(fID)% stencil(i,j)% U_z(:,k) = this% IBMStencilPoints(domain)% U_z(n,:)
      end do
#endif
   end subroutine IBM_HO_IBMstencilGradient

   subroutine IBM_MPI_GatherStancilState( this, nEqn, faces, time )
      use MPI_Process_info
      use BoundaryConditions
      use Meshtypes
      implicit none

      class(IBM_type), intent(inout) :: this
      type(face),      intent(inout) :: faces(:)
      integer,         intent(in)    :: nEqn
      real(kind=RP),   intent(in)    :: time 

      real(kind=RP) :: dt 
      integer       :: domains, fID, i, j, k, domain, ID

      do fID = 1, size(faces)
         associate( f => faces(fID) )
         if( .not. f% HO_IBM ) cycle
         do i = 0, f% Nf(1); do j = 0, f% Nf(2)
            dt = time - f% stencil(i,j)% time
            f% stencil(i,j)% time = time 
            call f% stencil(i,j)% ComputeState( f% geom% normal(:,i,j), f% geom% t1(:,i,j), f% geom% t2(:,i,j), f% STLNum, dt )
         end do; end do
         end associate
      end do

   end subroutine IBM_MPI_GatherStancilState
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------
!  HO Integration points procedure
!  ------------------------------------------------
   subroutine IBM_buildHOIntegrationPoints( this, elements, STLnum )
      use MPI_Process_Info
      implicit none

      class(IBM_type), intent(inout) :: this
      type(element),   intent(inout) :: elements(:)
      integer,         intent(in)    :: STLNum

      real(kind=RP) :: Point(NDIM)
      integer       :: i, j, domain, counter

      domain = MPI_Process% rank + 1

      allocate(this% IBM_HOIntegrationPoints(MPI_Process% nProcs))

      call this% IBM_HOIntegrationPoints(domain)% build(this% stl(STLNum)% NumOfObjs*NumOfIntegrationVertices)

      counter = 1

      do i = 1, this% stl(STLNum)% NumOfObjs
         do j = 1, NumOfIntegrationVertices
            this% IBM_HOIntegrationPoints(domain)% x(counter)% coords         = this% stl(STLNum)% ObjectsList(i)% IntegrationVertices(j)% coords
            this% IBM_HOIntegrationPoints(domain)% x(counter)% local_position = (/i,j,0/)
            counter = counter + 1
         end do
      end do
#ifdef _HAS_MPI_
      call MPIProcedures_IBM_HOIntegrationPoints( this% IBM_HOIntegrationPoints )
#endif
      call IBM_HO_findElements( this% IBM_HOIntegrationPoints, elements, this% NumOfSTL, this% clipAxis )

   end subroutine IBM_buildHOIntegrationPoints

   subroutine IBM_RelaxingSourceTerm( this, nEqn, x, dt, Q, h, eID, STLNum, Source )
      use MPI_Process_Info
      implicit none 

      class(IBM_type),   intent(inout) :: this
      integer,           intent(in)    :: STLNum, nEqn, eID 
      real(kind=RP),     intent(in)    :: dt, h, Q(nEqn) 
      real(kind=RP),     intent(inout) :: Source(nEqn), x(NDIM)

      integer       :: i, j, ierr 
      real(kind=RP) :: dist, dist_, epsilon, xP(NDIM)

      Source = 0.0_RP 

      if( .not. this% stl(STLNum)% move ) return 
      if( .not.  OBB(STLNum)% isPointInside( x, 1.5_RP) ) return 

      dist_ = huge(1.0_RP)

      do i = 1, this% stl(STLNum)% NumOfObjs 
         do j = 1, NumOfVertices
            xP     = this% stl(STLNum)% ObjectsList(i)% vertices(j)% coords
            dist_  = min(dist_,norm2(xP - x))
         end do 
      end do
#ifdef _HAS_MPI_
      call mpi_allreduce(dist_, dist, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, ierr)
#else 
      dist = dist_ 
#endif 
      epsilon = h**(1._RP/3._RP) 
      
      !call this% SourceTerm( nEqn, eID, x, dt, Q, STLNum, Source )

      Source = 0.5_RP*(1.0_RP + tanh(dist/epsilon)) * Source 

   end subroutine IBM_RelaxingSourceTerm
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------
!  Source terms for the immersed boundary
!  ------------------------------------------------
   subroutine IBM_SourceTerm( this, nEqn, eID, x, t, dt, Q, STLNum, Source )
      use PhysicsStorage
      implicit none
      !-arguments--------------------------------------------
      class(IBM_type),           intent(inout) :: this
      integer,                   intent(in)    :: eID, STLNum, nEqn
      real(kind=rp),             intent(inout) :: x(NDIM)
      real(kind=rp),             intent(in)    :: t, dt
      real(kind=rp),             intent(in)    :: Q(nEqn)
      real(kind=rp),             intent(inout) :: Source(nEqn)

      real(kind=RP) :: Qsb(nEqn), cL, cD

      Source = 0.0_RP
#if defined(NAVIERSTOKES)
      if( this% stl(STLNum)% move ) then 
         call BCsIBM(STLNum)% bc% FlowStateMoving_IBM( Q, x, t, dt, cL, cD, Qsb )
      else
         call BCsIBM(STLNum)% bc% FlowState_VPIBM( Q, x, Qsb ) 
      end if 
#endif
      Source = Q - Qsb 
      
      if( this% wallfunction ) then
         Source = -1.0_RP/(this% penalCoeff * this% penalization(eID)) * Source
      else
         Source = -1.0_RP/this% penalization(eID) * Source
      end if
   end subroutine IBM_SourceTerm

   subroutine IBM_TurbulentSourceTerm( this, nEqn, eID, x, nHat, Q, STLNum, Source )
      use PhysicsStorage
      implicit none
      !-arguments--------------------------------------------
      class(IBM_type),  intent(inout) :: this
      integer,          intent(in)    :: eID, STLNum, nEqn
      type(point_type), intent(inout) :: x
      real(kind=rp),    intent(in)    :: Q(nEqn), nHat(NDIM)
      real(kind=rp),    intent(inout) :: Source(nEqn)

      real(kind=RP) :: Q_ref(nEqn), dWall_ref, dWall, Qsb(nEqn)

      Source = 0.0_RP

      Q_ref     = x% Q
      dWall_ref = 2.0_RP * this% IP_dWall 
      dWall     = this% IP_dWall 
      Qsb       = Q 

      call ImagePointState( nEqn, Q_ref, dWall_ref, nHat, dWall, Qsb )

      call this% GetPointInterpolation( x, this% BandRegion(STLnum)% IBMmask ) 

      Source = Q - Qsb 
      
      if( this% wallfunction ) then
         Source = -1.0_RP/(this% penalCoeff * this% penalization(eID)) * Source
      else
         Source = -1.0_RP/this% penalization(eID) * Source
      end if

   end subroutine IBM_TurbulentSourceTerm
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
         !Q_F     = ForcingPointState( Q,                                        &
         !                             this% IP_Distance, e% geom% dWall(i,j,k), &
         !                             e% geom% normal(:,i,j,k), NCONS           )
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

   subroutine IBM_GetImplicitStep( this, element_, eID, dt )

      implicit none

      class(IBM_type), intent(inout) :: this
      type(element),   intent(inout) :: element_
      integer,         intent(in)    :: eID
      real(kind=RP),   intent(in)    :: dt

      real(kind=RP) :: invdS_dQ(NCONS,NCONS)
      integer       :: i, j, k, n

      do i = 0, element_% Nxyz(1); do j = 0, element_% Nxyz(2); do k = 0, element_% Nxyz(3)
         if( element_% isInsideBody(i,j,k) ) then
            call this% semiImplicitShiftJacobian( eID, element_% storage% Q(:,i,j,k), dt, invdS_dQ )
         else
            invdS_dQ = 0.0_RP
            do n = 1, NCONS
               invdS_dQ(n,n) = 1.0_RP
            end do
         end if
         element_% storage% Q(:,i,j,k) = element_% storage% Q(:,i,j,k) + matmul(invdS_dQ, dt * element_% storage% Qdot(:,i,j,k) )
      end do; end do; end do

   end subroutine IBM_GetImplicitStep

!
!  Second order Strang splitting correction Q^*. (1/dt - dS/dQ)^(-1)*Q^* = Q + dt*(S - dS/dQ*Q^*)
!  ------------------------------------------------------------------------------------------------
   subroutine IBM_GetSemiImplicitStep( this, nEqn, eID, x, t, dt, Q, STLNum )
      use PhysicsStorage
      implicit none
      !-arguments-----------------------------------------------------
      class(IBM_type),           intent(inout) :: this
      integer,                   intent(in)    :: nEqn, eID, STLNum
      real(kind=rp),             intent(inout) :: x(NDIM)
      real(kind=rp),             intent(in)    :: t, dt
      real(kind=rp),             intent(inout) :: Q(nEqn)
      !-local-variables-----------------------------------------------
      real(kind=rp) :: dS_dQ(nEqn,nEqn), invdS_dQ(nEqn,nEqn), &
                       Source(nEqn)

      call this% semiImplicitJacobian( eID, Q, dS_dQ )

      call this% semiImplicitShiftJacobian( eID, Q, dt, invdS_dQ )

      call this% SourceTerm( nEqn, eID, x, t, dt, Q, STLNum, Source )

      Q = matmul(invdS_dQ, Q + dt*( Source - matmul(dS_dQ,Q) ))

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

      !call this% SourceTerm( eID, Q, Q_target, Source, this% wallfunction )

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
      integer       :: eID, i, j, k

      if( .not. this% semiImplicit ) return
!$omp parallel
!$omp do schedule(runtime) private(i,j,k)
      do eID = 1, SIZE( elements )
         associate(e => elements(eID))
         do i = 0, e% Nxyz(1); do j = 0, e% Nxyz(2); do k = 0, e% Nxyz(3)
            if( e% isInsideBody(i,j,k) ) then
#if defined(NAVIERSTOKES) 
               call this% GetSemiImplicitStep( NCONS, eID, elements(eID)% geom% x(:,i,j,k), t, dt, e% storage% Q(:,i,j,k), e% STL(i,j,k) )
#endif
            end if
         end do; end do; end do
         end associate
      end do
!$omp end do
!$omp end parallel
    end subroutine IBM_SemiImplicitCorrection
!
!   State to be imposed on the forcing points due to the wall model
!   ---------------------------------------------------------------

   subroutine ImagePointState( nEqn, Q_ref, dWall_ref, nHat, dWall, Q ) 
      use PhysicsStorage
      use FluidData
      use VariableConversion
#if defined(NAVIERSTOKES)
      use WallFunctionDefinitions
      use WallFunctionBC
#endif
      implicit none
      !-arguments--------------------------------------------------------------
      integer,       intent(in)    :: nEqn
      real(kind=RP), intent(in)    :: Q_ref(nEqn)
      real(kind=RP), intent(in)    :: dWall_ref, dWall, nHat(NDIM) 
      real(kind=RP), intent(inout) :: Q(nEqn) 
      !-local-variables--------------------------------------------------------
      real(kind=RP) :: u_parallel(NDIM), x_II(NDIM), u_II, v_II, U_ref(NDIM), U(NDIM)
      real(kind=RP) :: nu, mu, kappa_ref, nu_ref, mu_ref, u_tau, y_plus, theta, P
#if defined (SPALARTALMARAS)
      real(kind=RP), parameter :: A = 17.0_RP 
      real(kind=RP)            :: DD, D2 
#endif 
#if defined(NAVIERSTOKES)
      U_ref      = Q_ref(IRHOU:IRHOW)/Q_ref(IRHO)
      u_parallel = U_ref - dot_product(U_ref, nHat) * nHat
      x_II       = u_parallel / norm2(u_parallel)
      u_II       = dot_product(U_ref, x_II)

      call get_laminar_mu_kappa(Q_ref, mu_ref, kappa_ref)
      nu_ref = mu_ref/Q_ref(IRHO)
 
      u_tau = u_tau_f(u_II, dWall_ref, nu_ref, u_tau0=.1_RP)

      call get_laminar_mu_kappa(Q,mu,kappa_ref)
      nu = mu/Q(IRHO)

      y_plus = y_plus_f(dWall, u_tau, nu)
      v_II   = u_plus_f(y_plus) * u_tau

      U = dot_product(U_ref,nHat) * nHat + v_II*x_II 
      P = pressure(Q)

      Q(IRHOU:IRHOW) = Q(IRHO) * U
      Q(IRHOE)       = P/thermodynamics % gammaMinus1 + 0.5_RP*Q(IRHO)*sum(U*U)
#if defined(SPALARTALMARAS)
      DD    = (1.0_RP - exp(-y_plus/A))
      D2    = DD*DD 
      theta = kappa * u_tau * dWall * D2
      Q(IRHOTHETA) = Q(IRHO) * theta
#endif
#endif
   end subroutine ImagePointState
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
      real(kind=rp),   intent(in)    :: Point(NDIM)
      integer,         intent(in)    :: STLNum
      integer,         intent(inout) :: NumOfIntersections
      !-local-variables------------------------------------------------------------
      type(KDtree),  pointer     :: tree
      real(kind=rp)              :: RayDirection(NDIM), vecAxis(NDIM), direction(NDIM)
      integer                    :: axis, lastIndex, OnBoundary_1, OnBoundary_2
      logical                    :: Upward
      type(ObjsDataLinkedList_t) :: Integer_List

      RayDirection = 0.0_RP

      axis = OBB(STLNum)% minAxis

      RayDirection(axis) = 1.0_RP

      vecAxis   = OBB(STLNum)% Vertices(:,7)
      direction = Point 

      NumOfIntersections = 0; OnBoundary_1 = 0; OnBoundary_2 = 0; lastIndex = -1
      Integer_List = ObjsDataLinkedList_t()
      do

         call this% root(STLNum)% FindLeaf( direction, tree, .false. )
         ! this check is needed when the ray 
         if( tree% index .eq. lastIndex ) call this% root(STLNum)% FindLeaf( direction, tree, .true. )

         lastIndex = tree% index

         if( tree% NumOfObjs .gt. 0 )                                                                &
         call isPointInside( Point, RayDirection, this% stl(STLNum)% ObjectsList, tree% ObjsIndeces, &
                             Integer_List, NumOfIntersections, OnBoundary_1, OnBoundary_2            )


         direction(axis) = tree% vertices(axis,7)
         if( direction(axis) .ge. vecAxis(axis) ) exit

      end do 

      NumOfIntersections = NumOfIntersections + abs(OnBoundary_1 - OnBoundary_2) 

      call integer_List% Destruct()

   end subroutine IBM_CheckPoint
!
!  Intersection between a ray an a set of triangles
!  ------------------------------------------------
   subroutine isPointInside( Point, RayDirection, ObjectsList, LeafIndeces, Integer_List, NumOfIntersections, OnBoundary_1, OnBoundary_2 )
      use RealDataLinkedList
      use omp_lib
      implicit none
      !-arguments----------------------------------------------------------------
      real(kind=rp),              intent(in)    :: Point(NDIM), RayDirection(NDIM)
      type(object_type),          intent(in)    :: ObjectsList(:)
      integer,                    intent(in)    :: LeafIndeces(:)
      type(ObjsDataLinkedList_t), intent(inout) :: Integer_List
      integer,                    intent(inout) :: NumOfIntersections, OnBoundary_1, OnBoundary_2
      !-local-variables----------------------------------------------------------
      real(kind=RP) :: n_sign
      logical :: Intersect, OnTriBound, found
      integer :: i, index

      do i = 1, size(LeafIndeces)
         index = LeafIndeces(i)
         found = integer_List% Check( index )
         if( index .eq. 0 ) cycle 
         if( .not. found ) then
            call Integer_List% Add( index )
            call PointIntersectTriangle( Point, ObjectsList(index)% vertices(1)% coords, &
                                                ObjectsList(index)% vertices(2)% coords, &
                                                ObjectsList(index)% vertices(3)% coords, &
                                                RayDirection, Intersect, OnTriBound      )
            if( Intersect )  NumOfIntersections = NumOfIntersections + 1

            n_sign = sign(1.0_RP,dot_product(RayDirection,ObjectsList(index)% normal))
            if( n_sign > 0.0_RP .and. OnTriBound ) then  
               OnBoundary_1 = 1
            elseif( n_sign < 0.0_RP .and. OnTriBound ) then 
               OnBoundary_2 = 1
            endif
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
                       Det, invDet, u, v, t
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

      if( u < 0.0_RP .or. (u-1.0_RP) > 0.0_RP ) return

      v = dot_product( RayDirection, Qvec )*invDet

      if( v < 0.0_RP .or. (u+v-1.0_RP) > 0 ) return

      t  = dot_product( E2vec, Qvec )*invDet 

      if( t > 0.0_RP .and. .not. almostEqual(t,0.0_RP) ) Intersect = .true.

      if( almostEqual(v,0.0_RP)            .or. &
          almostEqual((u+v-1.0_RP),0.0_RP) .or. &
          almostEqual((u-1.0_RP),0.0_RP)   .or. &
          almostEqual(u,0.0_RP)             ) then
         if(Intersect) then  
            OnTriBound = .true.
            intersect = .false.
         end if 
         return 
      end if  

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
                                       TriangleVertex3, dist, IntersectionPoint,s, t )
      use MappedGeometryClass
      implicit none
      !-arguments--------------------------------------------------------------------
      real(kind=rp), intent(in)  :: Point(:), TriangleVertex1(:), &
                                    TriangleVertex2(:),        &
                                    TriangleVertex3(:)
      real(kind=rp), intent(out) :: IntersectionPoint(NDIM)
      real(kind=rp), intent(out) :: dist, s, t
      !-local-variables--------------------------------------------------------------
      real(kind=rp) :: bb(NDIM), E0(NDIM), E1(NDIM), dd(NDIM), &
                       a, b, c, d, e, f, det, &! s, t,    &
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

   subroutine GetMatrixInterpolationSystem( point, NumOfInterPoints, IBMmask, INTERPOLATION )
      use DenseMatUtilities
      implicit none
      !-arguments------------------------------------
      type(point_type),  intent(inout) :: point
      integer,           intent(in)    :: NumOfInterPoints
      type(IBMpoints),   intent(in)    :: IBMmask(:) 
      integer,           intent(in)    :: INTERPOLATION
      !-local-variables------------------------------
      real(kind=RP) :: Phi(NumOfInterPoints,NumOfInterPoints), dist(NumOfInterPoints,NumOfInterPoints), d
      integer       :: i, j

      select case( INTERPOLATION )
      case( EXPONENTIAL )
         do i = 1, NumOfInterPoints
            do j = i, NumOfInterPoints
               dist(i,j) = norm2(IBMmask(point% domains(i))% coords(point% indeces(i),:) - IBMmask(point% domains(j))% coords(point% indeces(j),:))
               dist(j,i) = dist(i,j)
            end do
         end do
         expCoeff = 0.001_RP
         do i = 1, NumOfInterPoints
            do j = i, NumOfInterPoints
               Phi(i,j) = interpolationfunction(dist(i,j), EXPONENTIAL )
               Phi(j,i) = Phi(i,j)
            end do
            d    = norm2(point% coords - IBMmask(point% domains(i))% coords(point% indeces(i),:))
            point% b(i) = interpolationfunction(d, EXPONENTIAL)
         end do
         point% invPhi = inverse(Phi)

      case( IDW )
         point% b = 0.0_RP; point% invPhi = 0.0_RP
         do i = 1, NumOfInterPoints
            d           = norm2(point% coords - IBMmask(point% domains(i))% coords(point% indeces(i),:))
            point% invPhi(i,i) = 1.0_RP/d
            point% b           = point% b + point% invPhi(i,i)
         end do
         do i = 1, NumOfInterPoints
            point% b(i) = 1.0_RP/point% b(i)
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

   subroutine WriteTimeFile( t, filename )
      use MPI_Process_Info
      implicit none
      !-arguments---------------------------------------------------------
      real(kind=RP)    , intent(in) :: t
      character(len=*) , intent(in) :: filename
      !-local-variables---------------------------------------------------
      integer :: funit

      if( .not. MPI_Process% isRoot ) return

      funit = UnusedUnit()

      open(funit,file='RESULTS/'//trim(filename)//'_time.dat', action = "write" , access = "append" , status = "unknown")

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
