!
!//////////////////////////////////////////////////////
!
!   @File:    StorageClass.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Thu Oct  5 09:17:17 2017
!   @Last revision date: Wed Apr 11 16:48:12 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: bbf58e3f84d7b02811cdad12966f27e084fd8a36
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module StorageClass
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

   type Storage_t
      real(kind=RP), dimension(:,:,:,:),  allocatable :: Q
      real(kind=RP), dimension(:,:,:,:),  allocatable :: QDot
      real(kind=RP), dimension(:,:,:,:),  allocatable :: G
      real(kind=RP), dimension(:,:,:,:),  allocatable :: S
      real(kind=RP), dimension(:,:,:,:),  allocatable :: U_x
      real(kind=RP), dimension(:,:,:,:),  allocatable :: U_y
      real(kind=RP), dimension(:,:,:,:),  allocatable :: U_z
      type(Statistics_t)                              :: stats
#if defined(CAHNHILLIARD)
      real(kind=RP), dimension(:,:,:),   allocatable :: c     ! Cahn-Hilliard concentration
      real(kind=RP), dimension(:,:,:,:), allocatable :: gradC ! Cahn-Hilliard concentration gradient
      real(kind=RP), dimension(:,:,:),   allocatable :: mu    ! Cahn-Hilliard chemical pot.
      real(kind=RP), dimension(:,:,:,:), allocatable :: v   ! External flow field
#endif
      contains
         procedure   :: Construct => Storage_Construct
         procedure   :: Destruct  => Storage_Destruct
   end type Storage_t

   type FaceStorage_t
      real(kind=RP), dimension(:,:,:),     allocatable :: Q
      real(kind=RP), dimension(:,:,:),     allocatable :: U_x, U_y, U_z
      real(kind=RP), dimension(:,:,:),     allocatable :: FStar
      real(kind=RP), dimension(:,:,:,:)  , allocatable :: unStar
      real(kind=RP), dimension(:,:,:,:)  , allocatable :: dFStar_dqF   ! In storage(1), it stores dFStar/dqL, and in storage(2), it stores dFStar/dqR on the mortar points
      real(kind=RP), dimension(:,:,:,:,:), allocatable :: dFStar_dqEl  ! Stores both dFStar/dqL and dFStar/dqR on the face-element points of the corresponding side
#if defined(CAHNHILLIARD)
      real(kind=RP), dimension(:,:), allocatable :: c 
      real(kind=RP), dimension(:,:), allocatable :: mu 
      real(kind=RP), dimension(:,:,:), allocatable :: v 
#endif
      contains
         procedure   :: Construct => FaceStorage_Construct
         procedure   :: Destruct => FaceStorage_Destruct
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
      subroutine Storage_Construct(self, Nx, Ny, Nz, nEqn, nGradEqn, computeGradients)
         implicit none
         class(Storage_t)     :: self
         integer, intent(in)  :: Nx, Ny, Nz
         integer, intent(in)  :: nEqn, nGradEqn
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
         ALLOCATE( self % Q   (nEqn,0:Nx,0:Ny,0:Nz) )
         ALLOCATE( self % QDot(nEqn,0:Nx,0:Ny,0:Nz) )
         ALLOCATE( self % G   (nEqn,0:Nx,0:Ny,0:Nz) )
         ALLOCATE( self % S   (nEqn,0:Nx,0:Ny,0:Nz) )
         
         IF ( computeGradients )     THEN
            ALLOCATE( self % U_x(nGradEqn,0:Nx,0:Ny,0:Nz) )
            ALLOCATE( self % U_y(nGradEqn,0:Nx,0:Ny,0:Nz) )
            ALLOCATE( self % U_z(nGradEqn,0:Nx,0:Ny,0:Nz) )
         END IF

#if defined(CAHNHILLIARD)
         allocate( self % mu(0:Nx, 0:Ny, 0:Nz) )
         allocate( self % c (0:Nx, 0:Ny, 0:Nz) )
         allocate( self % gradC (1:NDIM,0:Nx, 0:Ny, 0:Nz) )
         allocate( self % v (1:NDIM,0:Nx, 0:Ny, 0:Nz) )
#endif
!         
!        -----------------
!        Initialize memory
!        -----------------
!
         self % G           = 0.0_RP
         self % S           = 0.0_RP
         self % Q           = 0.0_RP
         self % QDot        = 0.0_RP

#if defined(CAHNHILLIARD)
         self % mu    = 0.0_RP
         self % c     = 0.0_RP
         self % gradC = 0.0_RP
         self % v     = 0.0_RP
#endif
      
         IF ( computeGradients )     THEN
            self % U_x         = 0.0_RP
            self % U_y         = 0.0_RP
            self % U_z         = 0.0_RP
         END IF

      end subroutine Storage_Construct

      subroutine Storage_Destruct(self)
         implicit none
         class(Storage_t)     :: self
   
         safedeallocate(self % Q)
         safedeallocate(self % QDot)
         safedeallocate(self % G)
         safedeallocate(self % S)
         safedeallocate(self % U_x)
         safedeallocate(self % U_y)
         safedeallocate(self % U_z)

#if defined(CAHNHILLIARD)
         safedeallocate(self % mu)
         safedeallocate(self % c)
         safedeallocate(self % gradC)
         safedeallocate(self % v)
#endif

         call self % stats % Destruct()

      end subroutine Storage_Destruct
!
!////////////////////////////////////////////////////////////////////////////////////////////
!
!        Face storage procedures
!        -----------------------
!
!////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine FaceStorage_Construct(self, NDIM, Nf, Nel, nEqn, nGradEqn, computeGradients)
         implicit none
         class(FaceStorage_t)     :: self
         integer, intent(in)  :: NDIM
         integer, intent(in)  :: Nf(2)              ! Face polynomial order
         integer, intent(in)  :: Nel(2)             ! Element face polynomial order
         integer, intent(in)  :: nEqn, nGradEqn     ! Equations
         logical              :: computeGradients 
!
!        ---------------
!        Local variables
!        ---------------
!
!        ----------------
!        Volume variables
!        ----------------
!
         ALLOCATE( self % Q   (nEqn,0:Nf(1),0:Nf(2)) )
         allocate( self % fStar(nEqn, 0:Nel(1), 0:Nel(2)) )
         
         allocate( self % dFStar_dqF (nEqn,nEqn, 0: Nf(1), 0: Nf(2)) )
         allocate( self % dFStar_dqEl(nEqn,nEqn, 0:Nel(1), 0:Nel(2),2) )
         
         ALLOCATE( self % U_x(nGradEqn,0:Nf(1),0:Nf(2)) )
         ALLOCATE( self % U_y(nGradEqn,0:Nf(1),0:Nf(2)) )
         ALLOCATE( self % U_z(nGradEqn,0:Nf(1),0:Nf(2)) )
         ALLOCATE( self % unStar(nGradEqn,NDIM,0:Nel(1),0:Nel(2)) )

#if defined(CAHNHILLIARD)
         allocate( self % mu(0:Nf(1),0:Nf(2)) )
         allocate( self % c (0:Nf(1),0:Nf(2)) )
         allocate( self % v (1:NDIM,0:Nf(1),0:Nf(2)) )
#endif
!
!        -----------------
!        Initialize memory
!        -----------------
!
         self % Q           = 0.0_RP
         self % fStar       = 0.0_RP
      
         self % U_x         = 0.0_RP
         self % U_y         = 0.0_RP
         self % U_z         = 0.0_RP
         self % unStar      = 0.0_RP

#if defined(CAHNHILLIARD)
         self % mu = 0.0_RP
         self % c  = 0.0_RP
         self % v  = 0.0_RP
#endif

      end subroutine FaceStorage_Construct

      subroutine FaceStorage_Destruct(self)
         implicit none
         class(FaceStorage_t)     :: self
   
         safedeallocate(self % Q)
         safedeallocate(self % fStar)
         safedeallocate(self % U_x)
         safedeallocate(self % U_y)
         safedeallocate(self % unStar)
         safedeallocate(self % U_z)
         safedeallocate(self % dFStar_dqF)
         safedeallocate(self % dFStar_dqEl)
#if defined(CAHNHILLIARD)
         safedeallocate(self % mu)
         safedeallocate(self % c)
         safedeallocate(self % v)
#endif

      end subroutine FaceStorage_Destruct
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
