!
!//////////////////////////////////////////////////////
!
!   @File:    StorageClass.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Thu Oct  5 09:17:17 2017
!   @Last revision date: Sat May 12 21:51:10 2018
!   @Last revision author: Juan (juan.manzanero@upm.es)
!   @Last revision commit: 0a98ff59a5332051367a2a5c89543fa1ed797190
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module StorageClass
   use, intrinsic :: iso_c_binding
   use SMConstants
   use PhysicsStorage
   implicit none

   private
   public   ElementStorage_t, FaceStorage_t, Storage_t
   public   GetStorageEquations
  
   enum, bind(C)
      enumerator :: OFF = 0, NS, C, MU
   end enum 

   type Statistics_t
      real(kind=RP), dimension(:,:,:,:),  allocatable    :: data
      contains
         procedure   :: Construct => Statistics_Construct
         procedure   :: Destruct  => Statistics_Destruct
   end type Statistics_t
!  
!  Class for storing variables in the whole domain
!  ***********************************************
   type Storage_t
      real(kind=RP),                 pointer     :: Q(:)
      real(kind=RP),                 pointer     :: QDot(:)
      real(kind=RP),                 pointer     :: PrevQ(:,:)
#if defined(NAVIERSTOKES)
      real(kind=RP), dimension(:)  , allocatable :: QdotNS
      real(kind=RP), dimension(:)  , allocatable :: QNS
      real(kind=RP), dimension(:,:), allocatable :: PrevQNS ! Previous solution(s) in the whole domain
#endif
#if defined(CAHNHILLIARD)
      real(kind=RP), dimension(:)  , allocatable :: cDot
      real(kind=RP), dimension(:)  , allocatable :: c
      real(kind=RP), dimension(:,:), allocatable :: Prevc(:,:)
#endif      
      contains
         procedure :: Construct => Storage_Construct
         procedure :: Destruct  => Storage_Destruct
   end type Storage_t
!  
!  Class for pointing to previous solutions in an element
!  ******************************************************
   type ElementPrevSol_t
#if defined(NAVIERSTOKES)
      real(kind=RP), dimension(:,:,:,:),  pointer     :: QNS
#endif
#if defined(CAHNHILLIARD)
      real(kind=RP), dimension(:,:,:,:),  pointer     :: c
#endif
   end type ElementPrevSol_t
!  
!  Class for storing variables element-wise
!     (Q and Qdot are not owned by ElementStorage_t) 
!  ****************************************
   type ElementStorage_t
      integer                                         :: currentlyLoaded
      integer                                         :: NDOF              ! Number of degrees of freedom of element
      integer                                         :: firstIdx          ! Position in the global solution array
      real(kind=RP), dimension(:,:,:,:),  pointer, contiguous     :: Q           ! Pointers to the appropriate storage (NS or CH)
      real(kind=RP), dimension(:,:,:,:),  pointer, contiguous     :: QDot        !
      real(kind=RP), dimension(:,:,:,:),  pointer, contiguous     :: U_x         !
      real(kind=RP), dimension(:,:,:,:),  pointer, contiguous     :: U_y         !
      real(kind=RP), dimension(:,:,:,:),  pointer, contiguous     :: U_z         !
      type(ElementPrevSol_t),  allocatable :: PrevQ(:)           ! Previous solution
#if defined(NAVIERSTOKES)
      real(kind=RP),           pointer    , contiguous :: QNS(:,:,:,:)         ! NSE State vector
      real(kind=RP), private,  pointer    , contiguous :: QDotNS(:,:,:,:)      ! NSE State vector time derivative
      real(kind=RP), private,  allocatable :: U_xNS(:,:,:,:)       ! NSE x-gradients
      real(kind=RP), private,  allocatable :: U_yNS(:,:,:,:)       ! NSE y-gradients
      real(kind=RP), private,  allocatable :: U_zNS(:,:,:,:)       ! NSE z-gradients
      real(kind=RP),           allocatable :: gradRho(:,:,:,:)
      real(kind=RP),           allocatable :: G_NS(:,:,:,:)        ! NSE auxiliar storage
      real(kind=RP),           allocatable :: S_NS(:,:,:,:)        ! NSE source term
      type(Statistics_t)                   :: stats                ! NSE statistics
!
!     Pointers to face Jacobians
!     -> Currently only for df/dq⁺, For off-diagonal blocks, add df/dq⁻
!     ----------------------------------------------------------------
      real(kind=RP), dimension(:,:,:,:), pointer      :: dfdq_fr  ! FRONT
      real(kind=RP), dimension(:,:,:,:), pointer      :: dfdq_ba  ! BACK
      real(kind=RP), dimension(:,:,:,:), pointer      :: dfdq_bo  ! BOTTOM
      real(kind=RP), dimension(:,:,:,:), pointer      :: dfdq_to  ! TOP
      real(kind=RP), dimension(:,:,:,:), pointer      :: dfdq_ri  ! RIGHT
      real(kind=RP), dimension(:,:,:,:), pointer      :: dfdq_le  ! LEFT
      
      real(kind=RP), dimension(:,:,:,:,:,:), pointer    :: dfdGradQ_fr  ! FRONT
      real(kind=RP), dimension(:,:,:,:,:,:), pointer    :: dfdGradQ_ba  ! BACK
      real(kind=RP), dimension(:,:,:,:,:,:), pointer    :: dfdGradQ_bo  ! BOTTOM
      real(kind=RP), dimension(:,:,:,:,:,:), pointer    :: dfdGradQ_to  ! TOP
      real(kind=RP), dimension(:,:,:,:,:,:), pointer    :: dfdGradQ_ri  ! RIGHT
      real(kind=RP), dimension(:,:,:,:,:,:), pointer    :: dfdGradQ_le  ! LEFT
#endif
#if defined(CAHNHILLIARD)
      real(kind=RP), dimension(:,:,:,:),   pointer    , contiguous :: c     ! CHE concentration
      real(kind=RP), dimension(:,:,:,:),   pointer    , contiguous :: cDot  ! CHE concentration time derivative
      real(kind=RP), dimension(:,:,:,:),   allocatable :: c_x   ! CHE concentration x-gradient
      real(kind=RP), dimension(:,:,:,:),   allocatable :: c_y   ! CHE concentration y-gradient
      real(kind=RP), dimension(:,:,:,:),   allocatable :: c_z   ! CHE concentration z-gradient
      real(kind=RP), dimension(:,:,:,:),   allocatable :: mu    ! CHE chemical potential
      real(kind=RP), dimension(:,:,:,:),   allocatable :: mu_x  ! CHE chemical potential x-gradient
      real(kind=RP), dimension(:,:,:,:),   allocatable :: mu_y  ! CHE chemical potential y-gradient
      real(kind=RP), dimension(:,:,:,:),   allocatable :: mu_z  ! CHE chemical potential z-gradient
      real(kind=RP), dimension(:,:,:,:),   allocatable :: v     ! CHE flow field velocity
      real(kind=RP), dimension(:,:,:,:),   allocatable :: G_CH  ! CHE auxiliar storage   
#endif
      contains
         procedure   :: Assign            => ElementStorage_Assign
         generic     :: assignment(=)     => Assign
         procedure   :: Construct         => ElementStorage_Construct
         procedure   :: Destruct          => ElementStorage_Destruct
#if defined(NAVIERSTOKES)
         procedure   :: SetStorageToNS    => ElementStorage_SetStorageToNS
#endif
#if defined(CAHNHILLIARD)
         procedure   :: SetStorageToCH_c  => ElementStorage_SetStorageToCH_c
         procedure   :: SetStorageToCH_mu => ElementStorage_SetStorageToCH_mu
#endif
   end type ElementStorage_t
!  
!  Class for storing variables in the faces
!  ****************************************
   type FaceStorage_t
      integer                                          :: currentlyLoaded
      real(kind=RP), dimension(:,:,:),     pointer     :: Q
      real(kind=RP), dimension(:,:,:),     pointer     :: U_x, U_y, U_z
      real(kind=RP), dimension(:,:,:),     pointer     :: FStar
      real(kind=RP), dimension(:,:,:,:),   pointer     :: unStar
      real(kind=RP), dimension(:),         allocatable :: genericInterfaceFluxMemory ! unStar and fStar point to this memory simultaneously. This seems safe.
#if defined(NAVIERSTOKES)
      real(kind=RP), dimension(:,:,:),     allocatable :: QNS
      real(kind=RP), dimension(:,:,:),     allocatable :: U_xNS, U_yNS, U_zNS
      real(kind=RP), dimension(:,:,:),     allocatable :: gradRho          ! Gradient of density
      ! Inviscid Jacobians
      real(kind=RP), dimension(:,:,:,:)  , allocatable :: dFStar_dqF   ! In storage(1), it stores dFStar/dqL, and in storage(2), it stores dFStar/dqR on the mortar points
      real(kind=RP), dimension(:,:,:,:,:), allocatable :: dFStar_dqEl  ! Stores both dFStar/dqL and dFStar/dqR on the face-element points of the corresponding side
      ! Viscous Jacobians
      real(kind=RP), dimension(:,:,:,:,:,:), allocatable :: dFv_dGradQF  ! In storage(1), it stores dFv*/d∇qL, and in storage(2), it stores dFv*/d∇qR on the mortar points
      real(kind=RP), dimension(:,:,:,:,:,:), allocatable :: dFv_dGradQEl ! In storage(1), it stores dFv*/d∇qL, and in storage(2), it stores dFv*/d∇qR on the face-element points ... NOTE: this is enough for the diagonal blocks of the Jacobian, for off-diagonal blocks the crossed quantities must be computed and stored
#endif
#if defined(CAHNHILLIARD)
      real(kind=RP), dimension(:,:,:),   allocatable :: c
      real(kind=RP), dimension(:,:,:),   allocatable :: c_x
      real(kind=RP), dimension(:,:,:),   allocatable :: c_y
      real(kind=RP), dimension(:,:,:),   allocatable :: c_z
      real(kind=RP), dimension(:,:,:),   allocatable :: mu
      real(kind=RP), dimension(:,:,:),   allocatable :: mu_x
      real(kind=RP), dimension(:,:,:),   allocatable :: mu_y
      real(kind=RP), dimension(:,:,:),   allocatable :: mu_z
      real(kind=RP), dimension(:,:,:),   allocatable :: v
#endif
      contains
         procedure   :: Construct => FaceStorage_Construct
         procedure   :: Destruct  => FaceStorage_Destruct
#if defined(NAVIERSTOKES)
         procedure   :: SetStorageToNS => FaceStorage_SetStorageToNS
#endif
#if defined(CAHNHILLIARD)
         procedure   :: SetStorageToCH_c  => FaceStorage_SetStorageToCH_c
         procedure   :: SetStorageToCH_mu => FaceStorage_SetStorageToCH_mu
#endif
   end type FaceStorage_t
!
!  ========
   contains
!  ========
!
!
!///////////////////////////////////////////////////////////////////////////////////////////
!
!           Global Storage procedures
!           --------------------------
!
!///////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Storage_Construct(self, NDOF, bdf_order)
         implicit none
         !----------------------------------------------
         class(Storage_t)    :: self
         integer, intent(in) :: NDOF
         integer, intent(in) :: bdf_order
         !----------------------------------------------
         
#if defined(NAVIERSTOKES)
         allocate ( self % QNS   (NCONS*NDOF) )
         allocate ( self % QdotNS(NCONS*NDOF) )
         
         if (bdf_order /= 0) then
            allocate ( self % PrevQNS (NCONS*NDOF,bdf_order) )
         end if
#endif
#if defined(CAHNHILLIARD)
         allocate ( self % c(NCOMP*NDOF) )
         allocate ( self % cDot(NCOMP*NDOF) )
         
         if (bdf_order /= 0) then
            allocate ( self % PrevC (NCOMP*NDOF,bdf_order) )
         end if
#endif
      end subroutine Storage_Construct
!
!/////////////////////////////////////////////////
!
      subroutine Storage_Destruct(self)
         implicit none
         class(Storage_t)    :: self

#if defined(NAVIERSTOKES)         
         safedeallocate(self % QNS)
         safedeallocate(self % QdotNS)
         safedeallocate(self % PrevQNS)
#endif
#if defined(CAHNHILLIARD)
         safedeallocate(self % c)
         safedeallocate(self % cDot)
         safedeallocate(self % PrevC)
#endif
      end subroutine Storage_Destruct
!
!///////////////////////////////////////////////////////////////////////////////////////////
!
!           Element Storage procedures
!           --------------------------
!
!///////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine ElementStorage_Construct(self, Nx, Ny, Nz, computeGradients,globalStorage,firstIdx)
         implicit none
         !------------------------------------------------------------
         class(ElementStorage_t)     :: self                               !<> Storage to be constructed
         integer        , intent(in) :: Nx, Ny, Nz                         !<  Polynomial orders in every direction
         logical        , intent(in) :: computeGradients                   !<  Compute gradients?
         type(Storage_t), intent(in), target :: globalStorage    !<? Global storage to point to
         integer        , intent(in)         :: firstIdx         !<? Position of the solution of the element in the global array
         !------------------------------------------------------------
         integer :: k, num_prevSol
         integer :: bounds(2)
         !------------------------------------------------------------
!
!        --------------------------------
!        Get number of degrees of freedom
!        --------------------------------
!
         self % NDOF = (Nx + 1) * (Ny + 1) * (Nz + 1)
!
!        ----------------
!        Volume variables
!        ----------------
!
#if defined(NAVIERSTOKES)
         bounds(1) = (firstIdx-1)*NCONS + 1
         bounds(2) = bounds(1) + NCONS * self % NDOF - 1
         self % QNS   (1:NCONS,0:Nx,0:Ny,0:Nz) => globalStorage % QNS   (bounds(1) : bounds(2))
         self % QdotNS(1:NCONS,0:Nx,0:Ny,0:Nz) => globalStorage % QdotNS(bounds(1) : bounds(2))
         ! Previous solution
         num_prevSol = size(globalStorage % PrevQ,2)
         allocate ( self % PrevQ(num_prevSol) )
         do k=1, num_prevSol
            self % PrevQ(k) % QNS(1:NCONS,0:Nx,0:Ny,0:Nz) => globalStorage % PrevQNS(bounds(1):bounds(2),k)
         end do

         ALLOCATE( self % G_NS   (NCONS,0:Nx,0:Ny,0:Nz) )
         ALLOCATE( self % S_NS   (NCONS,0:Nx,0:Ny,0:Nz) )
         
         ALLOCATE( self % U_xNS (NGRAD,0:Nx,0:Ny,0:Nz) )
         ALLOCATE( self % U_yNS (NGRAD,0:Nx,0:Ny,0:Nz) )
         ALLOCATE( self % U_zNS (NGRAD,0:Nx,0:Ny,0:Nz) )
         allocate( self % gradRho(NDIM,0:Nx,0:Ny,0:Nz) )
!
!        Point to NS by default
!        ----------------------
         call self % SetStorageToNS
#endif

#if defined(CAHNHILLIARD)
         bounds(1) = (firstIdx-1)*NCOMP + 1
         bounds(2) = bounds(1) + NCOMP * self % NDOF - 1
         self % c   (1:NCOMP,0:Nx,0:Ny,0:Nz) => globalStorage % c   (bounds(1) : bounds(2))
         self % cDot(1:NCOMP,0:Nx,0:Ny,0:Nz) => globalStorage % cDot(bounds(1) : bounds(2))
         ! Previous solution
         num_prevSol = size(globalStorage % PrevC,2)
         allocate ( self % PrevQ(num_prevSol) )
         do k=1, num_prevSol
            self % PrevQ(k) % c(1:NCOMP,0:Nx,0:Ny,0:Nz) => globalStorage % PrevC(bounds(1):bounds(2),k)
         end do

         allocate(self % c_x (NCOMP, 0:Nx, 0:Ny, 0:Nz))
         allocate(self % c_y (NCOMP, 0:Nx, 0:Ny, 0:Nz))
         allocate(self % c_z (NCOMP, 0:Nx, 0:Ny, 0:Nz))
         allocate(self % mu  (NCOMP, 0:Nx, 0:Ny, 0:Nz))
         allocate(self % mu_x(NCOMP, 0:Nx, 0:Ny, 0:Nz))
         allocate(self % mu_y(NCOMP, 0:Nx, 0:Ny, 0:Nz))
         allocate(self % mu_z(NCOMP, 0:Nx, 0:Ny, 0:Nz))
         ALLOCATE(self % G_CH(NCOMP,0:Nx,0:Ny,0:Nz) )
         allocate(self % v   (1:NDIM, 0:Nx, 0:Ny, 0:Nz))
#endif
!         
!        -----------------
!        Initialize memory
!        -----------------
!
#if defined(NAVIERSTOKES)
         self % G_NS   = 0.0_RP
         self % S_NS   = 0.0_RP
         self % QNS    = 0.0_RP
         self % QDotNS = 0.0_RP
         
         self % U_xNS = 0.0_RP
         self % U_yNS = 0.0_RP
         self % U_zNS = 0.0_RP
#endif

#if defined(CAHNHILLIARD)
         self % c     = 0.0_RP
         self % c_x   = 0.0_RP
         self % c_y   = 0.0_RP
         self % c_z   = 0.0_RP
         self % mu    = 0.0_RP
         self % mu_x  = 0.0_RP
         self % mu_y  = 0.0_RP
         self % mu_z  = 0.0_RP
         self % G_CH  = 0.0_RP
         self % v     = 0.0_RP
#endif
      
      end subroutine ElementStorage_Construct

      subroutine ElementStorage_Assign(to, from)
!
!        ***********************************
!        We need an special assign procedure
!        ***********************************
!
         implicit none
         class(ElementStorage_t), intent(out) :: to
         type(ElementStorage_t),  intent(in)  :: from
!
!        Copy the storage
!        ----------------
#if defined(NAVIERSTOKES)
         to % QNS    = from % QNS
         to % U_xNS  = from % U_xNS
         to % U_yNS  = from % U_yNS
         to % U_zNS  = from % U_zNS
         to % QDotNS = from % QDotNS
         to % G_NS   = from % G_NS
#endif
#if defined(CAHNHILLIARD)
         to % c    = from % c
         to % c_x  = from % c_x
         to % c_y  = from % c_y
         to % c_z  = from % c_z
         to % mu   = from % mu
         to % mu_x = from % mu_x
         to % mu_y = from % mu_y
         to % mu_z = from % mu_z
         to % v    = from % v
         to % cDot = from % cDot
         to % G_CH = from % G_CH
#endif

         select case ( to % currentlyLoaded ) 
         case (OFF)
            to % Q    => NULL()
            to % U_x  => NULL()
            to % U_y  => NULL()
            to % U_z  => NULL()
            to % QDot => NULL()

#if defined(NAVIERSTOKES)
         case (NS)
            call to % SetStorageToNS   
#endif
#if defined(CAHNHILLIARD)
         case (C)
            call to % SetStorageToCH_c

         case (MU)
            call to % SetStorageToCH_mu
#endif
         end select

      end subroutine ElementStorage_Assign
!
!///////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine ElementStorage_Destruct(self)
         implicit none
         class(ElementStorage_t) :: self
         integer                 :: num_prevSol, k

         self % currentlyLoaded = OFF
         self % NDOF = 0
         self % firstIdx = -1
         
         self % Q => NULL()
         self % QDot => NULL()
         self % U_x => NULL()
         self % U_y => NULL()
         self % U_z => NULL()

#if defined(NAVIERSTOKES)
         self % QNS => NULL()
         self % QDotNS => NULL()

         num_prevSol = size(self % PrevQ)
         do k=1, num_prevSol
            nullify( self % PrevQ(k) % QNS )
         end do

         safedeallocate(self % G_NS)
         safedeallocate(self % S_NS)
         safedeallocate(self % U_xNS)
         safedeallocate(self % U_yNS)
         safedeallocate(self % U_zNS)
         safedeallocate(self % gradRho)
         
         nullify ( self % dfdq_fr )
         nullify ( self % dfdq_ba )
         nullify ( self % dfdq_bo )
         nullify ( self % dfdq_to )
         nullify ( self % dfdq_ri )
         nullify ( self % dfdq_le )
         
         nullify ( self % dfdGradQ_fr )
         nullify ( self % dfdGradQ_ba )
         nullify ( self % dfdGradQ_bo )
         nullify ( self % dfdGradQ_to )
         nullify ( self % dfdGradQ_ri )
         nullify ( self % dfdGradQ_le )
#endif
#if defined(CAHNHILLIARD)
         self % c => NULL()
         self % cDot => NULL()

         num_prevSol = size(self % PrevQ)
         do k=1, num_prevSol
            nullify( self % PrevQ(k) % c )
         end do

         safedeallocate(self % c_x)
         safedeallocate(self % c_y)
         safedeallocate(self % c_z)
         safedeallocate(self % mu)
         safedeallocate(self % mu_x)
         safedeallocate(self % mu_y)
         safedeallocate(self % mu_z)
         safedeallocate(self % G_CH)
         safedeallocate(self % v)
#endif
         safedeallocate(self % PrevQ)

      end subroutine ElementStorage_Destruct
#if defined(NAVIERSTOKES)
      subroutine ElementStorage_SetStorageToNS(self)
!
!        *****************************************
!        This subroutine selects the Navier-Stokes
!        state vector as current storage.
!        *****************************************
!
         implicit none
         class(ElementStorage_t), target   :: self

         self % currentlyLoaded = NS
         self % Q   (1:,0:,0:,0:) => self % QNS
         self % U_x (1:,0:,0:,0:) => self % U_xNS
         self % U_y (1:,0:,0:,0:) => self % U_yNS
         self % U_z (1:,0:,0:,0:) => self % U_zNS
         self % QDot(1:,0:,0:,0:) => self % QDotNS

      end subroutine ElementStorage_SetStorageToNS
#endif
#if defined(CAHNHILLIARD)
      subroutine ElementStorage_SetStorageToCH_c(self)
!
!        *********************************************
!        This subroutine selects the concentration as
!        current storage.
!        *********************************************
!
         implicit none
         class(ElementStorage_t), target   :: self
      
         self % currentlyLoaded = C
!
!        Point to the one dimensional pointers with generic arrays
!        ---------------------------------------------------------
         self % Q   (1:,0:,0:,0:) => self % c
         self % U_x (1:,0:,0:,0:) => self % c_x
         self % U_y (1:,0:,0:,0:) => self % c_y
         self % U_z (1:,0:,0:,0:) => self % c_z
         self % QDot(1:,0:,0:,0:) => self % cDot
   
      end subroutine ElementStorage_SetStorageToCH_c

      subroutine ElementStorage_SetStorageToCH_mu(self)
!
!        *************************************************
!        This subroutine selects the chemical potential as
!        current storage, with the particularity that
!        selects also cDot as QDot.
!        *************************************************
!
         implicit none
         class(ElementStorage_t), target   :: self

         self % currentlyLoaded = MU

         self % Q   (1:,0:,0:,0:) => self % mu
         self % U_x (1:,0:,0:,0:) => self % mu_x
         self % U_y (1:,0:,0:,0:) => self % mu_y
         self % U_z (1:,0:,0:,0:) => self % mu_z
         self % QDot(1:,0:,0:,0:) => self % cDot
   
      end subroutine ElementStorage_SetStorageToCH_mu
#endif
!
!////////////////////////////////////////////////////////////////////////////////////////////
!
!        Face storage procedures
!        -----------------------
!
!////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine FaceStorage_Construct(self, NDIM, Nf, Nel, computeGradients)
         implicit none
         class(FaceStorage_t)     :: self
         integer, intent(in)  :: NDIM
         integer, intent(in)  :: Nf(2)              ! Face polynomial order
         integer, intent(in)  :: Nel(2)             ! Element face polynomial order
         logical              :: computeGradients 
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: interfaceFluxMemorySize

         interfaceFluxMemorySize = 0

#if defined(NAVIERSTOKES)

         ALLOCATE( self % QNS   (NCONS,0:Nf(1),0:Nf(2)) )
         ALLOCATE( self % U_xNS(NGRAD,0:Nf(1),0:Nf(2)) )
         ALLOCATE( self % U_yNS(NGRAD,0:Nf(1),0:Nf(2)) )
         ALLOCATE( self % U_zNS(NGRAD,0:Nf(1),0:Nf(2)) )
!
!        Biggest Interface flux memory size is u\vec{n}
!        ----------------------------------------------
         interfaceFluxMemorySize = NGRAD * nDIM * product(Nf + 1)
!
!        TODO: JMT, if (implicit..?)
         allocate( self % dFStar_dqF (NCONS,NCONS, 0: Nf(1), 0: Nf(2)) )
         allocate( self % dFStar_dqEl(NCONS,NCONS, 0:Nel(1), 0:Nel(2),2) )
         
         allocate( self % dFv_dGradQF (NCONS,NCONS,NDIM,2,0: Nf(1),0: Nf(2)) )
         allocate( self % dFv_dGradQEl(NCONS,NCONS,NDIM,2,0:Nel(1),0:Nel(2)) )
         
         allocate( self % gradRho   (NDIM,0:Nf(1),0:Nf(2)) )
#endif
#if defined(CAHNHILLIARD)
         allocate(self % c   (NCOMP , 0:Nf(1), 0:Nf(2)))
         allocate(self % c_x (NCOMP , 0:Nf(1), 0:Nf(2)))
         allocate(self % c_y (NCOMP , 0:Nf(1), 0:Nf(2)))
         allocate(self % c_z (NCOMP , 0:Nf(1), 0:Nf(2)))
         allocate(self % mu  (NCOMP , 0:Nf(1), 0:Nf(2)))
         allocate(self % mu_x(NCOMP , 0:Nf(1), 0:Nf(2)))
         allocate(self % mu_y(NCOMP , 0:Nf(1), 0:Nf(2)))
         allocate(self % mu_z(NCOMP , 0:Nf(1), 0:Nf(2)))
         allocate(self % v   (1:NDIM, 0:Nf(1), 0:Nf(2)))
!
!        CH will never be the biggest memory requirement unless NSE are disabled
!        -----------------------------------------------------------------------
         interfaceFluxMemorySize = max(interfaceFluxMemorySize, NCOMP*nDIM*product(Nf+1))
#endif
!
!        Reserve memory for the interface fluxes
!        ---------------------------------------
         allocate(self % genericInterfaceFluxMemory(interfaceFluxMemorySize))

#if defined(NAVIERSTOKES)
!
!        Point to NS by default
!        ----------------------
         call self % SetStorageToNS
#endif
!
!        -----------------
!        Initialize memory
!        -----------------
!
#if defined(NAVIERSTOKES)
         self % QNS    = 0.0_RP
         
         self % U_xNS = 0.0_RP
         self % U_yNS = 0.0_RP
         self % U_zNS = 0.0_RP

         self % dFStar_dqF  = 0.0_RP
         self % dFStar_dqEl = 0.0_RP
#endif

#if defined(CAHNHILLIARD)
         self % c     = 0.0_RP
         self % c_x   = 0.0_RP
         self % c_y   = 0.0_RP
         self % c_z   = 0.0_RP
         self % mu    = 0.0_RP
         self % mu_x  = 0.0_RP
         self % mu_y  = 0.0_RP
         self % mu_z  = 0.0_RP
         self % v     = 0.0_RP
#endif

      end subroutine FaceStorage_Construct

      subroutine FaceStorage_Destruct(self)
         implicit none
         class(FaceStorage_t)     :: self
   
         self % currentlyLoaded = OFF

#if defined(NAVIERSTOKES)
         safedeallocate(self % QNS)
         safedeallocate(self % U_xNS)
         safedeallocate(self % U_yNS)
         safedeallocate(self % U_zNS)
         safedeallocate(self % dFStar_dqF)
         safedeallocate(self % dFStar_dqEl)
         safedeallocate(self % dFv_dGradQF)
         safedeallocate(self % dFv_dGradQEl)
         safedeallocate(self % gradRho)
#endif
#if defined(CAHNHILLIARD)
         safedeallocate(self % c)
         safedeallocate(self % c_x)
         safedeallocate(self % c_y)
         safedeallocate(self % c_z)
         safedeallocate(self % mu)
         safedeallocate(self % mu_x)
         safedeallocate(self % mu_y)
         safedeallocate(self % mu_z)
         safedeallocate(self % v)
#endif
         safedeallocate(self % genericInterfaceFluxMemory)

         self % Q      => NULL()
         self % U_x    => NULL() ; self % U_y => NULL() ; self % U_z => NULL()
         self % unStar => NULL()
         self % fStar  => NULL()

      end subroutine FaceStorage_Destruct
#if defined(NAVIERSTOKES)
      subroutine FaceStorage_SetStorageToNS(self)
         implicit none
         class(FaceStorage_t), target    :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: Nx, Ny

         self % currentlyLoaded = NS
!
!        Get sizes
!        ---------
         Nx   = size(self % QNS,2) - 1 
         Ny   = size(self % QNS,3) - 1 

         self % Q   (1:,0:,0:)            => self % QNS
         self % fStar(1:NCONS, 0:Nx, 0:Ny) => self % genericInterfaceFluxMemory

         self % genericInterfaceFluxMemory = 0.0_RP

         if ( allocated(self % U_xNS) ) then
            self % U_x (1:,0:,0:) => self % U_xNS
            self % U_y (1:,0:,0:) => self % U_yNS
            self % U_z (1:,0:,0:) => self % U_zNS
            self % unStar(1:NGRAD, 1:NDIM, 0:Nx, 0:Ny) => self % genericInterfaceFluxMemory
         end if

      end subroutine FaceStorage_SetStorageToNS
#endif
#if defined(CAHNHILLIARD)
      subroutine FaceStorage_SetStorageToCH_c(self)
         implicit none
         class(FaceStorage_t), target  :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: Nx, Ny

         self % currentlyLoaded = C
!
!        Get sizes
!        ---------
         Nx   = size(self % c,2) - 1 
         Ny   = size(self % c,3) - 1 

         self % Q(1:,0:,0:)   => self % c
         self % U_x(1:,0:,0:) => self % c_x
         self % U_y(1:,0:,0:) => self % c_y
         self % U_z(1:,0:,0:) => self % c_z

         self % fStar(1:NCOMP,0:Nx,0:Ny)            => self % genericInterfaceFluxMemory
         self % unStar(1:NCOMP, 1:NDIM, 0:Nx, 0:Ny) => self % genericInterfaceFluxMemory

         self % genericInterfaceFluxMemory = 0.0_RP

      end subroutine FaceStorage_SetStorageToCH_c

      subroutine FaceStorage_SetStorageToCH_mu(self)
         implicit none
         class(FaceStorage_t), target  :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: Nx, Ny

         self % currentlyLoaded = MU
!
!        Get sizes
!        ---------
         Nx   = size(self % mu,2) - 1 
         Ny   = size(self % mu,3) - 1 

         self % Q(1:,0:,0:)   => self % mu
         self % U_x(1:,0:,0:) => self % mu_x
         self % U_y(1:,0:,0:) => self % mu_y
         self % U_z(1:,0:,0:) => self % mu_z

         self % fStar(1:NCOMP,0:Nx,0:Ny)            => self % genericInterfaceFluxMemory
         self % unStar(1:NCOMP, 1:NDIM, 0:Nx, 0:Ny) => self % genericInterfaceFluxMemory

         self % genericInterfaceFluxMemory = 0.0_RP

      end subroutine FaceStorage_SetStorageToCH_mu
#endif
!
!/////////////////////////////////////////////////////////////////////////////////////
!
!           Statistics procedures
!           ---------------------
!
!/////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Statistics_Construct(self, no_of_variables, N)
         implicit none
         class(Statistics_t)           :: self
         integer,          intent(in)  :: no_of_variables
         integer,          intent(in)  :: N(3)
!
!        Allocate and initialize
!        -----------------------
         allocate( self % data(no_of_variables, 0:N(1), 0:N(2), 0:N(3) ) ) 
         self % data = 0.0_RP

      end subroutine Statistics_Construct
   
      subroutine Statistics_Destruct(self)
         implicit none
         class(Statistics_t)     :: self

         safedeallocate( self % data )

      end subroutine Statistics_Destruct

      subroutine GetStorageEquations(off_, ns_, c_, mu_)
         implicit none  
         integer, intent(out) :: off_, ns_, c_, mu_

         off_ = OFF
         ns_  = NS
         c_   = C
         mu_  = MU

      end subroutine GetStorageEquations

end module StorageClass
