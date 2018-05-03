!
!//////////////////////////////////////////////////////
!
!   @File:    StorageClass.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Thu Oct  5 09:17:17 2017
!   @Last revision date: Thu May  3 16:26:16 2018
!   @Last revision author: Juan (juan.manzanero@upm.es)
!   @Last revision commit: 5a86eb6fbfa5f685edfa7826a0b6714de7b3cf7c
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
   public   Storage_t, FaceStorage_t

   type Statistics_t
      real(kind=RP), dimension(:,:,:,:),  allocatable    :: data
      contains
         procedure   :: Construct => Statistics_Construct
         procedure   :: Destruct  => Statistics_Destruct
   end type Statistics_t


   enum, bind(C)
      enumerator :: OFF = 0, NS, C, MU
   end enum 

   type Storage_t
!
!     *************************
!     Generic storage variables
!     *************************
!
      integer                                         :: currentlyLoaded
      real(kind=RP), dimension(:,:,:,:),  pointer     :: Q           ! Pointers to the appropriate storage (NS or CH)
      real(kind=RP), dimension(:,:,:,:),  pointer     :: QDot        !
      real(kind=RP), dimension(:,:,:,:),  pointer     :: U_x         !
      real(kind=RP), dimension(:,:,:,:),  pointer     :: U_y         !
      real(kind=RP), dimension(:,:,:,:),  pointer     :: U_z         !
#if defined(NAVIERSTOKES)
      real(kind=RP),           allocatable :: QNS(:,:,:,:)         ! NSE State vector
      real(kind=RP), private,  allocatable :: U_xNS(:,:,:,:)       ! NSE x-gradients
      real(kind=RP), private,  allocatable :: U_yNS(:,:,:,:)       ! NSE y-gradients
      real(kind=RP), private,  allocatable :: U_zNS(:,:,:,:)       ! NSE z-gradients
      real(kind=RP), private,  allocatable :: QDotNS(:,:,:,:)      ! NSE State vector time derivative
      real(kind=RP),           allocatable :: G_NS(:,:,:,:)        ! NSE auxiliar storage
      real(kind=RP),           allocatable :: S(:,:,:,:)           ! NSE source term
      type(Statistics_t)                   :: stats                ! NSE statistics
#endif
#if defined(CAHNHILLIARD)
      real(kind=RP), dimension(:,:,:,:),   allocatable :: c     ! CHE concentration
      real(kind=RP), dimension(:,:,:,:),   allocatable :: c_x   ! CHE concentration x-gradient
      real(kind=RP), dimension(:,:,:,:),   allocatable :: c_y   ! CHE concentration y-gradient
      real(kind=RP), dimension(:,:,:,:),   allocatable :: c_z   ! CHE concentration z-gradient
      real(kind=RP), dimension(:,:,:,:),   allocatable :: cDot  ! CHE concentration time derivative
      real(kind=RP), dimension(:,:,:,:),   allocatable :: mu    ! CHE chemical potential
      real(kind=RP), dimension(:,:,:,:),   allocatable :: mu_x  ! CHE chemical potential x-gradient
      real(kind=RP), dimension(:,:,:,:),   allocatable :: mu_y  ! CHE chemical potential y-gradient
      real(kind=RP), dimension(:,:,:,:),   allocatable :: mu_z  ! CHE chemical potential z-gradient
      real(kind=RP), dimension(:,:,:,:),   allocatable :: v     ! CHE flow field velocity
      real(kind=RP), dimension(:,:,:,:),   allocatable :: G_CH  ! CHE auxiliar storage   
#endif
      contains
         procedure   :: Construct         => Storage_Construct
         procedure   :: Destruct          => Storage_Destruct
#if defined(NAVIERSTOKES)
         procedure   :: SetStorageToNS    => Storage_SetStorageToNS
#endif
#if defined(CAHNHILLIARD)
         procedure   :: SetStorageToCH_c  => Storage_SetStorageToCH_c
         procedure   :: SetStorageToCH_mu => Storage_SetStorageToCH_mu
#endif
   end type Storage_t

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
      real(kind=RP), dimension(:,:,:,:),   allocatable :: dFStar_dqF   ! In storage(1), it stores dFStar/dqL, and in storage(2), it stores dFStar/dqR on the mortar points
      real(kind=RP), dimension(:,:,:,:,:), allocatable :: dFStar_dqEl  ! Stores both dFStar/dqL and dFStar/dqR on the face-element points of the corresponding side
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
!///////////////////////////////////////////////////////////////////////////////////////////
!
!           Storage procedures
!           ------------------
!
!///////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Storage_Construct(self, Nx, Ny, Nz, computeGradients)
         implicit none
         class(Storage_t)     :: self
         integer, intent(in)  :: Nx, Ny, Nz
         logical              :: computeGradients
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: Nmax
!
!        --------------
!        Get dimensions         
!        --------------
!
!        ----------------
!        Volume variables
!        ----------------
!
#if defined(NAVIERSTOKES)
         ALLOCATE( self % QNS    (NCONS,0:Nx,0:Ny,0:Nz) )
         ALLOCATE( self % QDotNS (NCONS,0:Nx,0:Ny,0:Nz) )
         ALLOCATE( self % G_NS   (NCONS,0:Nx,0:Ny,0:Nz) )
         ALLOCATE( self % S      (NCONS,0:Nx,0:Ny,0:Nz) )
         
         IF ( computeGradients )     THEN
            ALLOCATE( self % U_xNS (NGRAD,0:Nx,0:Ny,0:Nz) )
            ALLOCATE( self % U_yNS (NGRAD,0:Nx,0:Ny,0:Nz) )
            ALLOCATE( self % U_zNS (NGRAD,0:Nx,0:Ny,0:Nz) )
         END IF
!
!        Point to NS by default
!        ----------------------
         call self % SetStorageToNS
#endif

#if defined(CAHNHILLIARD)
         allocate(self % c   (NCOMP, 0:Nx, 0:Ny, 0:Nz))
         allocate(self % c_x (NCOMP, 0:Nx, 0:Ny, 0:Nz))
         allocate(self % c_y (NCOMP, 0:Nx, 0:Ny, 0:Nz))
         allocate(self % c_z (NCOMP, 0:Nx, 0:Ny, 0:Nz))
         allocate(self % mu  (NCOMP, 0:Nx, 0:Ny, 0:Nz))
         allocate(self % mu_x(NCOMP, 0:Nx, 0:Ny, 0:Nz))
         allocate(self % mu_y(NCOMP, 0:Nx, 0:Ny, 0:Nz))
         allocate(self % mu_z(NCOMP, 0:Nx, 0:Ny, 0:Nz))
         allocate(self % cDot(NCOMP, 0:Nx, 0:Ny, 0:Nz))
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
         self % S      = 0.0_RP
         self % QNS    = 0.0_RP
         self % QDotNS = 0.0_RP
         
         IF ( computeGradients )     THEN
            self % U_xNS = 0.0_RP
            self % U_yNS = 0.0_RP
            self % U_zNS = 0.0_RP
         END IF
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
      
      end subroutine Storage_Construct

      subroutine Storage_Destruct(self)
         implicit none
         class(Storage_t)     :: self
   
         self % currentlyLoaded = OFF

#if defined(NAVIERSTOKES)
         safedeallocate(self % QNS)
         safedeallocate(self % QDotNS)
         safedeallocate(self % G_NS)
         safedeallocate(self % S)
         safedeallocate(self % U_xNS)
         safedeallocate(self % U_yNS)
         safedeallocate(self % U_zNS)

         call self % stats % Destruct()
#endif

#if defined(CAHNHILLIARD)
         safedeallocate(self % c)
         safedeallocate(self % c_x)
         safedeallocate(self % c_y)
         safedeallocate(self % c_z)
         safedeallocate(self % cDot)
         safedeallocate(self % mu)
         safedeallocate(self % mu_x)
         safedeallocate(self % mu_y)
         safedeallocate(self % mu_z)
         safedeallocate(self % G_CH)
         safedeallocate(self % v)
#endif

         self % Q    => NULL()
         self % U_x  => NULL()
         self % U_y  => NULL()
         self % U_z  => NULL()
         self % QDot => NULL()
         


      end subroutine Storage_Destruct
#if defined(NAVIERSTOKES)
      subroutine Storage_SetStorageToNS(self)
!
!        *****************************************
!        This subroutine selects the Navier-Stokes
!        state vector as current storage.
!        *****************************************
!
         implicit none
         class(Storage_t), target   :: self

         self % currentlyLoaded = NS
         self % Q   (1:,0:,0:,0:) => self % QNS
         self % U_x (1:,0:,0:,0:) => self % U_xNS
         self % U_y (1:,0:,0:,0:) => self % U_yNS
         self % U_z (1:,0:,0:,0:) => self % U_zNS
         self % QDot(1:,0:,0:,0:) => self % QDotNS

      end subroutine Storage_SetStorageToNS
#endif
#if defined(CAHNHILLIARD)
      subroutine Storage_SetStorageToCH_c(self)
!
!        *********************************************
!        This subroutine selects the concentration as
!        current storage.
!        *********************************************
!
         implicit none
         class(Storage_t), target   :: self
      
         self % currentlyLoaded = C
!
!        Point to the one dimensional pointers with generic arrays
!        ---------------------------------------------------------
         self % Q   (1:,0:,0:,0:) => self % c
         self % U_x (1:,0:,0:,0:) => self % c_x
         self % U_y (1:,0:,0:,0:) => self % c_y
         self % U_z (1:,0:,0:,0:) => self % c_z
         self % QDot(1:,0:,0:,0:) => self % cDot
   
      end subroutine Storage_SetStorageToCH_c

      subroutine Storage_SetStorageToCH_mu(self)
!
!        *************************************************
!        This subroutine selects the chemical potential as
!        current storage, with the particularity that
!        selects also cDot as QDot.
!        *************************************************
!
         implicit none
         class(Storage_t), target   :: self

         self % currentlyLoaded = MU

         self % Q   (1:,0:,0:,0:) => self % mu
         self % U_x (1:,0:,0:,0:) => self % mu_x
         self % U_y (1:,0:,0:,0:) => self % mu_y
         self % U_z (1:,0:,0:,0:) => self % mu_z
         self % QDot(1:,0:,0:,0:) => self % cDot
   
      end subroutine Storage_SetStorageToCH_mu
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
         if (computeGradients) then
            ALLOCATE( self % U_xNS(NGRAD,0:Nf(1),0:Nf(2)) )
            ALLOCATE( self % U_yNS(NGRAD,0:Nf(1),0:Nf(2)) )
            ALLOCATE( self % U_zNS(NGRAD,0:Nf(1),0:Nf(2)) )
!
!           Biggest Interface flux memory size is u\vec{n}
!           ----------------------------------------------
            interfaceFluxMemorySize = NGRAD * nDIM * product(Nf + 1)

         else
!
!           Biggers Interface flux memory is fStar
!           --------------------------------------
            interfaceFluxMemorySize = NCONS * product(Nf + 1)

         end if
!
!        TODO: JMT, if (implicit..?)
         allocate( self % dFStar_dqF (NCONS,NCONS, 0: Nf(1), 0: Nf(2)) )
         allocate( self % dFStar_dqEl(NCONS,NCONS, 0:Nel(1), 0:Nel(2),2) )

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
         
         IF ( computeGradients )     THEN
            self % U_xNS = 0.0_RP
            self % U_yNS = 0.0_RP
            self % U_zNS = 0.0_RP
         END IF

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

end module StorageClass
