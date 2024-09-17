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
      type(KDtree),               allocatable :: root(:), rootDistance(:)
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
                                                 iter = 0, twoD_axis = 0, N
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
         procedure :: constructDistance_KDtree            => IBM_constructDistance_KDtree
         procedure :: HOmask                              => IBM_HOmask
         procedure :: GetStencil                          => IBM_GetStencil
         procedure :: HO_IBMstencilState                  => IBM_HO_IBMstencilState
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
      real(kind=RP)              :: axis(NDIM)
      integer                    :: STLNum, j, k, maxAxis, minAxis

      call this% describe()

      allocate( this% stl(this% NumOfSTL),         &
                OBB(this% NumOfSTL),               &
                this% root(this% NumOfSTL),        &
                this% integral(this% NumOfSTL),    &
                this% STLfilename(this% NumOfSTL), &
                this% bctype(this% NumOfSTL),      &
                this% zoneMask(this% NumOfSTL)     )

      if( this% ComputeBandRegion ) allocate( this% BandRegion(this% NumOfSTL)   )
      if( this% ComputeDistance   ) allocate( this% rootDistance(this% NumOfSTL) )

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

         OBB(STLNum)% filename = this% STLfilename(STLNum)
         this% bctype(STLNum)  = this% stl(STLNum)% bctype
         
         if( MPI_Process% isRoot ) then
            this% stl(STLNum)% show = .true.
            call this% stl(STLNum)% ReadTessellation( this% STLfilename(STLNum) )
            if( this% ClipAxis .ne. 0 ) then
               call this% stl(STLNum)% Clip( this% minCOORDS, this% maxCOORDS, this% ClipAxis, .true. )
            else
               call this% stl(STLNum)% describe(this% STLfilename(STLNum))
            end if
         end if
         this% AAB = .true.
#ifdef _HAS_MPI_
         call this% MPI_sendSTLpartitions( STLNum, this% stl(STLNum)% maxAxis ) 
#endif
         if( this% ComputeBandRegion ) then 
            call this% stl(STLNum)% SetIntegrationPoints()
            if( .not. this% HO_IBM ) call this% stl(STLNum)% SetIntegration( this% NumOfInterPoints ) 
         end if 

         call OBB(STLNum)% construct( this% stl(STLNum), this% plotOBB, this% AAB, .false.)

         OBB(STLNum)% maxAxis = this% stl(STLNum)% maxAxis

         call OBB(STLNum)% GetGLobalVertices()

         OBB(STLNum)% minAxis = OBB(STLNum)% GetMinAxis( OBB(STLNum)% maxAxis, this% ClipAxis )

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
   subroutine IBM_DestroyKDtree( this, STLNum, isChild )
      use MPI_Process_Info
      implicit none
      !-arguments--------------------------------------------------
      class(IBM_type),           intent(inout) :: this
      integer,                   intent(in)    :: STLNum
      logical,                   intent(in)    :: isChild

      if( this% ComputeDistance ) call this% rootDistance(STLNum)% destruct( isChild )

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


   subroutine IBM_constructzonemask( this, faces, STLNum )
      use Meshtypes
      implicit none 

      class(IBM_type), intent(inout) :: this 
      type(face),      intent(inout) :: faces(:)
      integer,         intent(in)    :: STLNum

      integer :: fID, i, j 

      if( .not. this% stl(STLNum)% BFcorrection ) return  
      
      this% NumOfMaskObjs = 0

      do fID = 1, size(faces)
         associate( f => faces(fID) )
         if (f% FaceType .eq. HMESH_INTERIOR .or. f% FaceType .eq. HMESH_MPI) cycle
         if( trim(this% STLfilename(STLNum)) .eq. trim(f% boundaryName) ) then 
            f% HO_IBM = .true.
            allocate(f% stencil(0:f% Nf(1),0:f% Nf(2)))
            f% HOSIDE = minloc(f% elementIDs,dim=1)
            f% STLNum = STLNum
            do j = 0, f% Nf(2); do i = 0, f% Nf(1)
               f% stencil(i,j)% x  = f% geom% x(:,i,j)
               this% NumOfMaskObjs = this% NumOfMaskObjs + 1
            end do; end do
         end if 
         end associate
      end do 
 
      if( this% NumOfMaskObjs .gt. 0 ) this% zoneMask(STLNum) = .true.   

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
      type(IBMPoints) :: tocopy
      real(kind=RP)   :: Point(NDIM)
      integer         :: eID, n, i, j, k, NumOfObjs, domains, domain, NInters
#ifdef _HAS_MPI_
      integer         :: ierr
#endif
      integer, allocatable :: NumOfIntersections(:)

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
               if( OBB(STLNum)% isPointInsideAAB( coords = e% SurfInfo% corners(:,k) ) ) then
                  tocopy% NumOfObjs                            = tocopy% NumOfObjs + 1
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
               if( OBB(STLNum)% isPointInsideAAB( coords = e% geom% x(:,i,j,k) ) ) then
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
            call OBB(STLNum)% ChangeRefFrame( this% IBMmask(domains)% coords(i,:), LOCAL, Point )
            call this% CheckPoint( Point, STLNum, this% IBMmask(domains)% NumOfIntersections(i) )
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
      this% NumOfMaskObjs = 0

      do eID = 1, size(elements)
         associate( e => elements(eID) )
         if(  all(e% MaskCorners((/4,3,7,8/))) .and. .not. ALL(e% MaskCorners((/1,2,6,5/))) ) then
            associate(f => faces(e% faceIDs(EBACK)))             
            call SetHOface( f, eID, STLNum, this% NumOfMaskObjs )      
            end associate
         end if
          if( all(e% MaskCorners((/1,2,6,5/))) .and. .not. all( e% MaskCorners((/4,3,7,8/)) ) ) then
            associate(f => faces(e% faceIDs(EFRONT)))                 
            call SetHOface( f, eID, STLNum, this% NumOfMaskObjs )      
            end associate
         end if
         if( all(e% MaskCorners((/5,6,7,8/))) .and. .not. all(e% MaskCorners((/1,2,3,4/))) ) then
            associate(f => faces(e% faceIDs(ETOP)))
            call SetHOface( f, eID, STLNum, this% NumOfMaskObjs )  
            end associate
         end if
         if( all(e% MaskCorners((/1,2,3,4/))) .and. .not. all(e% MaskCorners((/5,6,7,8/))) ) then
            associate(f => faces(e% faceIDs(EBOTTOM)))
            call SetHOface( f, eID, STLNum, this% NumOfMaskObjs )  
            end associate
         end if
         if( all(e% MaskCorners((/1,4,8,5/))) .and. .not. all(e% MaskCorners((/2,3,7,6/))) ) then
            associate(f => faces(e% faceIDs(ELEFT)))
            call SetHOface( f, eID, STLNum, this% NumOfMaskObjs )   
            end associate
         end if
         if( all(e% MaskCorners((/2,3,7,6/))) .and. .not. all(e% MaskCorners((/1,4,8,5/))) ) then
            associate(f => faces(e% faceIDs(ERIGHT)))
            call SetHOface( f, eID, STLNum, this% NumOfMaskObjs )  
            end associate
         end if
         if(all(e% MaskCorners)) e% HO_IBM = .true.
         end associate
      end do

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
        f% HOeID    = 0
        f% HOdomain = 0
        f% STLNum   = 0
      else
         f% HO_IBM = .true.
         allocate(f% stencil(0:f% Nf(1),0:f% Nf(2)))
         f% HOSIDE = FindHOside( f, eID )
         f% STLNum = STLNum
         do j = 0, f% Nf(2); do i = 0, f% Nf(1)
            f% stencil(i,j)% x  = f% geom% x(:,i,j)
            NumOfMaskObjs       = NumOfMaskObjs + 1
         end do; end do
         f% HOeID = eID
         f% HOdomain = MPI_Process% rank + 1
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
         do k = 1, NODES_PER_ELEMENT
            this% maxCOORDS = max(this% maxCOORDS,e% SurfInfo% corners(axis,k)); this% minCOORDS = min(this% minCOORDS,e% SurfInfo% corners(axis,k))
         end do
         end associate
      end do
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
            write(funit,'(3E13.5)')  IBM% IBMStencilPoints(domain)% coords(i,IX), IBM% IBMStencilPoints(domain)% coords(i,IY), IBM% IBMStencilPoints(domain)% coords(i,IZ)
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
            write(funit,'(3E13.5)')  IBM% IBMStencilPoints(domain)% coords(i,IX), IBM% IBMStencilPoints(domain)% coords(i,IY), IBM% IBMStencilPoints(domain)% coords(i,IZ)
         end do
      end do

      close(funit)

   end subroutine PlotStencil

   subroutine PlotMPIFaces( IBM, faces, STLNum )
      use MPI_Process_Info
      use Meshtypes
      implicit none

      type(IBM_type), intent(in) :: IBM
      type(face),     intent(in) :: faces(:)
      integer,        intent(in) :: STLNum

      integer :: NumOfObjs, i, j, fID, funit
      character(len=LINE_LENGTH) :: rank, filename

      NumOfObjs = 0 

      write(rank,'(I0)') MPI_Process% rank + 1

      do fID = 1, size(faces)
         associate( f => faces(fID) )
         if( f% HO_IBM .and. f% STLNum .eq. STLNum ) then !.and. f% faceType .eq. HMESH_MPI ) then 
            NumOfObjs = NumOfObjs + (f% Nf(1) +1) * (f% Nf(2) + 1)
         end if 
         end associate
      end do 

      if( NumOfObjs .eq. 0 ) return 

      write(filename,'(A)') trim(IBM% STLfilename(STLNum))//'_'//trim(rank)

      call TecFileHeader( 'IBM/MPIFacePoints_'//trim(filename), 'Mask Points', NumOfObjs, funit )

      do fID = 1, size(faces)
         associate( f => faces(fID) )
         if( f% HO_IBM .and. f% STLNum .eq. STLNum  ) then !.and. f% faceType .eq. HMESH_MPI ) then 
            do i = 0, f% Nf(1); do j = 0, f% Nf(2)
               write(funit,'(3E13.5)')  f% geom% x(IX,i,j), f% geom% x(IY,i,j), f% geom% x(IZ,i,j)
            end do; end do  
         end if
         end associate
      end do

   end subroutine PlotMPIFaces
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

      call this% constructzonemask( faces, STLNum )
      call this% constructmask( elements, STLNum, no_of_DoFs, faces, iter )
#ifdef _HAS_MPI_
      if( this% HO_IBM ) call FixingmpiFaces( faces, MPIfaces, this% NumOfMaskObjs, this% MPIfixed )
#endif 
      if( this% ComputeDistance ) call this% ComputeIBMWallDistance( elements=elements, STLNum=STLNum, faces=faces )

      if( .not. movingSTL ) then 
         if( this% ComputeBandRegion ) call this% constructBandRegion( elements, no_of_DoFs, STLNum, NCONS )
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

      class(IBM_type), intent(inout) :: this
      type(element),   intent(inout) :: elements(:)
      type(face),      intent(inout) :: faces(:)

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
         call PlotMPIFaces( this, faces, STLNum )
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
         if( this% ComputeDistance ) call this% rootDistance(STLNum)% destruct( isChild )
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
      use NodalStorageClass
      use PolynomialInterpAndDerivsModule
      implicit none
      !-arguments------------------------------------------------------
      class(IBM_type), intent(inout) :: this
      type(element),   intent(inout) :: elements(:)
      integer,         intent(in)    :: no_of_DoFs, STLNum, nEqn
      !-local-variables------------------------------------------------
      type(IBMPoints)               :: tocopy
      integer                       :: eID, i, j, k, n, domains, domain, NumOfObjs
      type(NodalStorage_t), pointer :: spA 
      real(kind=RP)                 :: ds
      real(kind=RP), allocatable    :: nodes(:)
      
      if( .not. allocated(this% BandRegion(STLNum)% IBMmask) ) allocate(this% BandRegion(STLNum)% IBMmask(MPI_Process% nProcs))

      domain = MPI_Process% rank + 1

      call tocopy% build(no_of_DoFs)
      tocopy% NumOfObjs = 0
      
      if( this% HO_IBM ) then
         N = this% N 
         associate( spA => NodalStorage(N) )
         do i = 1, this% stl(STLNum)% NumOfObjs
            do j = 1, NumOfIntegrationVertices
               if( this% zoneMask(STLNum) ) then 
                  do k = 0, N 
                     tocopy% NumOfObjs                           = tocopy% NumOfObjs + 1
                     ds                                          = this% dl + 0.5_RP*(1._RP + spA% x(k))*this% L
                     tocopy% coords(tocopy% NumOfObjs,:)         = this% stl(STLNum)% ObjectsList(i)% IntegrationVertices(j)% coords + this% stl(STLNum)% ObjectsList(i)% normal * ds
                     tocopy% element_index(tocopy% NumOfObjs)    = 0
                     tocopy% local_position(tocopy% NumOfObjs,:) = (/i,j,k/)
                  end do
              else
                  tocopy% NumOfObjs                           = tocopy% NumOfObjs + 1
                  tocopy% coords(tocopy% NumOfObjs,:)         = this% stl(STLNum)% ObjectsList(i)% IntegrationVertices(j)% coords
                  tocopy% element_index(tocopy% NumOfObjs)    = 0
                  tocopy% local_position(tocopy% NumOfObjs,:) = (/i,j,0/)
               end if 
            end do
         end do
        end associate 
      else
         do eID = 1, size(elements)
            associate(e => elements(eID))
            if( any(e% isInsideBody) ) cycle
            do i = 0, e% Nxyz(1); do j = 0, e% Nxyz(2); do k = 0, e% Nxyz(3)
               if( OBB(STLNum)% isPointInside( coords = e% geom% x(:,i,j,k), coeff=this% BandRegionCoeff ) ) then
                  tocopy% NumOfObjs                           = tocopy% NumOfObjs + 1
                  tocopy% coords(tocopy% NumOfObjs,:)         = e% geom% x(:,i,j,k)
                  tocopy% element_index(tocopy% NumOfObjs)    = eID
                  tocopy% local_position(tocopy% NumOfObjs,:) = (/i,j,k/)
               end if
            end do; end do; end do
            end associate
         end do
      end if
      
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
      
      if( this% HO_IBM ) then 
         do domains = 1, MPI_Process% nProcs
            call this% BandRegion(STLNum)% IBMmask(domains)% buildBandRegion( this% BandRegion(STLNum)% IBMmask(domains)% NumOfObjs )
         end do  
         if( this% zoneMask(STLNum) ) then
            associate( spA => NodalStorage(N) )
            if( .not. allocated(this% BandRegion(STLNum)% IBMmask(domain)% lj) ) allocate( this% BandRegion(STLNum)% IBMmask(domain)% lj(0:N) )
            allocate(nodes(0:N))
            do i = 0, N 
               ds       = this% dl + 0.5_RP * (1._RP + spA% x(i))*this% L
               nodes(i) = 2.0_RP * ds/(this% dl + this% L) -1.0_RP
            end do 

            do i = 0, N
               this% BandRegion(STLNum)% IBMmask(domain)% lj(i) = LagrangeInterpolatingPolynomial( i, -1._RP, N, nodes )
            end do
            deallocate(nodes)
            end associate
         end if
         call IBM_HO_findElements( this% BandRegion(STLNum)% IBMmask, elements, this% NumOfSTL, this% clipAxis )
      end if 
      
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

      if( HO_IBM ) then
         call IBM_HO_GetState   ( this% IBMmask, elements, nEqn )
         call IBM_HO_GetGradient( this% IBMmask, elements, nEqn )
         if( zoneMask ) call stl% ResetIntegrationPoints()
         if( MPI_Process% doMPIAction ) then 
#ifdef _HAS_MPI_
            call GatherHOIntegrationPointsState( this% IBMmask, stl% ObjectsList, zoneMask, nEqn )
#endif
         else 
            if( zoneMask ) then 
               do n = 1, this% IBMmask(domain)% NumOfObjs
                  i = this% IBMmask(domain)% local_position(n,IX)
                  j = this% IBMmask(domain)% local_position(n,IY)
                  k = this% IBMmask(domain)% local_position(n,IZ)

                  stl% ObjectsList(i)% IntegrationVertices(j)% Q   = stl% ObjectsList(i)% IntegrationVertices(j)% Q + &
                                                                     this% IBMmask(domain)% lj(k) * this% IBMmask(domain)% Q(n,:)
                  stl% ObjectsList(i)% IntegrationVertices(j)% U_x = stl% ObjectsList(i)% IntegrationVertices(j)% U_x + &
                                                                     this% IBMmask(domain)% lj(k) * this% IBMmask(domain)% U_x(n,:)
                  stl% ObjectsList(i)% IntegrationVertices(j)% U_y = stl% ObjectsList(i)% IntegrationVertices(j)% U_y + &
                                                                     this% IBMmask(domain)% lj(k) * this% IBMmask(domain)% U_y(n,:)
                  stl% ObjectsList(i)% IntegrationVertices(j)% U_z = stl% ObjectsList(i)% IntegrationVertices(j)% U_z + &
                                                                     this% IBMmask(domain)% lj(k) * this% IBMmask(domain)% U_z(n,:)
               end do
            else
               do n = 1, this% IBMmask(domain)% NumOfObjs
                  i = this% IBMmask(domain)% local_position(n,IX)
                  j = this% IBMmask(domain)% local_position(n,IY)

                  stl% ObjectsList(i)% IntegrationVertices(j)% Q   = this% IBMmask(domain)% Q(n,:)
                  stl% ObjectsList(i)% IntegrationVertices(j)% U_x = this% IBMmask(domain)% U_x(n,:)
                  stl% ObjectsList(i)% IntegrationVertices(j)% U_y = this% IBMmask(domain)% U_y(n,:)
                  stl% ObjectsList(i)% IntegrationVertices(j)% U_z = this% IBMmask(domain)% U_z(n,:)
               end do 
            end if 
         end if 
      else
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
      end if

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

      if( any(this% stl(:)% move) )  this% MPIfixed = .false.
      
      do STLNum = 1, this% NumOfSTL

         if( .not. this% stl(STLNum)% move ) cycle 
         
         call OBB(STLNum)% ChangeObjsRefFrame( this% stl(STLNum)% ObjectsList, GLOBAL )
         do i = 1, this% stl(STLNum)% NumOfObjs
            do j = 1, NumOfVertices
#if defined(NAVIERSTOKES)
               call BCsIBM(STLNum)% bc% PositionMoving_IBM( this% stl(STLNum)% ObjectsList(i)% vertices(j)% coords, this% dt, cL, cD )
#endif 
            end do 
            call this% stl(STLNum)% updateNormals( this% stl(STLNum)% ObjectsList(i) )
            if( this% ComputeBandRegion ) then 
               do j = 1, NumOfIntegrationVertices
#if defined(NAVIERSTOKES)
                  call BCsIBM(STLNum)% bc% PositionMoving_IBM( this% stl(STLNum)% ObjectsList(i)% IntegrationVertices(j)% coords, this% dt, cL, cD )
#endif
               end do 
            end if 
         end do 
         
         call OBB(STLNum)% construct( this% stl(STLNum), this% plotOBB, this% AAB, .false.)

         if( autosave ) call plotSTL( this% stl(STLNum), iter )

         call OBB(STLNum)% GetGLobalVertices()
         call OBB(STLNum)% ChangeObjsRefFrame( this% stl(STLNum)% ObjectsList, LOCAL )
         
         this% plot              = autosave
         this% plotKDtree        = .false.
         this% stl(STLNum)% show = .false.
         
         call this% constructSTL_KDtree( STLNum )
         
         call this% build( elements, faces, MPIfaces, no_of_DoFs, STLNum, isChild, .true., iter )
         
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
         call this% constructBandRegion( elements, no_of_DoFs, STLNum, nEqn )
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
      real(kind=RP) :: Point(NDIM), normal(NDIM), dist, sign_
      integer       :: domain, domains, i, j, k, n, fID, eID
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
               if( OBB(STLNum)% isPointInside( coords = e% geom% x(:,i,j,k), coeff=this% BandRegionCoeff ) ) then
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
            call OBB(STLNum)% ChangeRefFrame(this% IBMmask(domains)% coords(n,:), LOCAL, Point)
            call MinimumDistance( Point, this% rootDistance(STLNum), dist, normal )
            if( Dist .lt. this% IBMmask(domains)% dist(n) ) then
               this% IBMmask(domains)% dist(n)     = dist
               this% IBMmask(domains)% normal(n,:) = normal
            end if
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

            sign_ = sign( 1.0_RP,dot_product(faces(fID)% geom% normal(:,i,j),this% IBMmask(domain)% normal(n,:)))

            faces(fID)% stencil(i,j)% dist = this% IBMmask(domain)% dist(n)
            if( this% zoneMask(STLNum) ) then 
               faces(fID)% stencil(i,j)% normal = -sign_*this% IBMmask(domain)% normal(n,:)
            else 
               faces(fID)% stencil(i,j)% normal = -this% IBMmask(domain)% normal(n,:)
            end if 
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
      implicit none

      class(IBM_type), intent(inout) :: this
      type(element),   intent(in)    :: elements(:)
      type(face),      intent(inout) :: faces(:)

      real(kind=RP)                 :: alpha, L, dist, dl, h, d_, dl_, L_
      real(kind=RP)                 :: x0(NDIM), normal(NDIM), xi(NDIM)
      integer                       :: fID, eID, N, M, i, j, k, NumOfObjs_, N_, NumOfObjs
      type(NodalStorage_t), pointer :: spA
#if _HAS_MPI_
      integer :: ierr 
#endif 

      spA_s% Constructed = .false. 

      d_         = 0.0_RP
      L_         = 0.0_RP
      dl_        = 0.0_RP
      NumOfObjs_ = 0 
      N_         = 0
      
      do fID = 1, size(faces)
         associate( f=> faces(fID) )
         if( f% HO_IBM ) then
            N = max(f% Nf(1),f% Nf(2)) 
            M = N - 1

            spA => NodalStorage(N)

            call spA_s% construct(GAUSS, M) 

            eID = f% elementIDs(maxloc(f% elementIDs, dim=1))
            h   = 0.0025_RP*sqrt(2.0_RP) * (elements(eID)% geom% Volume)**(1.0_RP/3.0_RP)

            L   = ABS(spA% x(0) - spA% x(1)) * h

            do i = 0, f% Nf(1); do j = 0, f% Nf(2)

                allocate(f% stencil(i,j)% nodes(0:N),f% stencil(i,j)% dWall(0:N))

               dist   = f% stencil(i,j)% dist
               normal = f% stencil(i,j)% normal 
               
               dl                 = dist + h
               x0                 = f% stencil(i,j)% x + dl * normal
               f% stencil(i,j)% N = N

               call f% stencil(i,j)% build( x0, normal, L, M )

               f% stencil(i,j)% x_s(:,0) = f% stencil(i,j)% x + dist * normal

               f% stencil(i,j)% d  = dist
               f% stencil(i,j)% L  = L
               f% stencil(i,j)% dl = dl

               d_         =  d_        + dist
               L_         =  L_        + L
               dl_        =  dl_       + dl
               NumOfObjs_ = NumOfObjs_ + 1 
               N_         = max(N_,N)
               if( this% Wallfunction ) f% stencil(i,j)% wallfunction = .true. 

            end do; end do
         end if
         end associate
      end do
#ifdef _HAS_MPI_
         call mpi_allreduce(NumOfObjs_, NumOfObjs, 1, MPI_INT   , MPI_SUM, MPI_COMM_WORLD, ierr)
         call mpi_allreduce(d_        , this% d  , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
         call mpi_allreduce(L_        , this% L  , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
         call mpi_allreduce(dl_       , this% dl , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
         call mpi_allreduce(N_        , this% N  , 1, MPI_INT   , MPI_MAX, MPI_COMM_WORLD, ierr)
#else
         this% d  = d_
         this% L  = L_
         this% dl = dl_
         this% N  = N_
#endif
      this% d  = this% d/NumOfObjs 
      this% L  = this% L/NumOfObjs 
      this% dl = this% dl/NumOfObjs 

   end subroutine IBM_GetStencil

   real(kind=RP) function hGeom( e, twoD_axis )

      implicit none 

      type(element), intent(in) :: e 
      integer,       intent(in) :: twoD_axis 

      real(kind=RP) :: L1, L2, L3
      integer       :: i, j 

      L1 = huge(1.0_RP); L2 = huge(1.0_RP); L3 = huge(1.0_RP)

      if( twoD_axis .ne. 1 ) then 
         L1 = norm2( e% SurfInfo% corners(:,4) - e% SurfInfo% corners(:,1) )   
      endif 
      if( twoD_axis .ne. 2 ) then
         L2 = norm2( e% SurfInfo% corners(:,2) - e% SurfInfo% corners(:,1) )
      end if 
      if( twoD_axis .ne. 3 ) then
         L3 = norm2( e% SurfInfo% corners(:,5) - e% SurfInfo% corners(:,1) )
      end if 

      hGeom = min(L1, L2, L3)  

   end function hGeom


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

   subroutine IBM_MPI_GatherStancilState( this, nEqn, faces, time )
      use MPI_Process_info

      use Meshtypes
      implicit none

      class(IBM_type), intent(inout) :: this
      type(face),      intent(inout) :: faces(:)
      integer,         intent(in)    :: nEqn
      real(kind=RP),   intent(in)    :: time 

      integer :: domains, fID, i, j, k, domain, ID

      do fID = 1, size(faces)
         associate( f => faces(fID) )
         if( .not. f% HO_IBM ) cycle
         do i = 0, f% Nf(1); do j = 0, f% Nf(2)
            f% stencil(i,j)% time = time 
            call f% stencil(i,j)% ComputeState( f% geom% normal(:,i,j), f% geom% t1(:,i,j), f% geom% t2(:,i,j), f% STLNum )
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
            call OBB(STLNum)% ChangeRefFrame(this% stl(STLNum)% ObjectsList(i)% IntegrationVertices(j)% coords, GLOBAL, Point)
            this% IBM_HOIntegrationPoints(domain)% x(counter)% coords         = Point
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

      integer       :: i, j
      real(kind=RP) :: dist, dist_, epsilon, xP(NDIM)

      Source = 0.0_RP 

      if( .not. this% stl(STLNum)% moving ) return 
      if( .not.  OBB(STLNum)% isPointInside(x, 1.5_RP) ) return 

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
      
      call this% SourceTerm( nEqn, eID, x, dt, Q, STLNum, Source )

      Source = 0.5_RP*(1.0_RP + tanh(dist/epsilon)) * Source 

   end subroutine IBM_RelaxingSourceTerm
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------------------
!  Source terms for the immersed boundary
!  ------------------------------------------------
   subroutine IBM_SourceTerm( this, nEqn, eID, x, dt, Q, STLNum, Source )
      use PhysicsStorage
      implicit none
      !-arguments--------------------------------------------
      class(IBM_type),           intent(inout) :: this
      integer,                   intent(in)    :: eID, STLNum, nEqn
      real(kind=rp),             intent(inout) :: x(NDIM)
      real(kind=rp),             intent(in)    :: dt
      real(kind=rp),             intent(in)    :: Q(nEqn)
      real(kind=rp),             intent(inout) :: Source(nEqn)

      real(kind=RP) :: Qsb(nEqn), cL, cD

      Source = 0.0_RP
#if defined(NAVIERSTOKES)
      if( this% stl(STLNum)% move ) then 
         call BCsIBM(STLNum)% bc% FlowStateMoving_IBM( Q, x, dt, cL, cD, Qsb )
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

      call this% SourceTerm( nEqn, eID, x, dt, Q, STLNum, Source )

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
