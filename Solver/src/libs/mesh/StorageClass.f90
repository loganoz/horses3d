!
!//////////////////////////////////////////////////////
!
!   @File:    StorageClass.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Thu Oct  5 09:17:17 2017
!   @Last revision date: Tue Jan 16 13:25:53 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 5143fa03eb24e8282a2043aa22cb178df572b474
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module StorageClass
   use SMConstants
   use PhysicsStorage
   implicit none

   private
   public   ElementStorage_t, FaceStorage_t, Storage_t

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
      real(kind=RP), dimension(:)  , allocatable :: Qdot
      real(kind=RP), dimension(:)  , allocatable :: Q
      real(kind=RP), dimension(:,:), allocatable :: PrevQ ! Previous solution(s) in the whole domain
      
      contains
         procedure :: Construct => Storage_Construct
         procedure :: Destruct  => Storage_Destruct
   end type Storage_t
   
   type ElementPrevSol_t
      real(kind=RP), dimension(:,:,:,:),  pointer     :: Q
   end type ElementPrevSol_t
!  
!  Class for storing variables element-wise
!     (Q and Qdot are not owned by ElementStorage_t) 
!  ****************************************
   type ElementStorage_t
      real(kind=RP), dimension(:,:,:,:), pointer, contiguous :: Q
      real(kind=RP), dimension(:,:,:,:), pointer, contiguous :: QDot
      type(ElementPrevSol_t)           ,  allocatable :: PrevQ(:)
      real(kind=RP), dimension(:,:,:,:),  allocatable :: G
      real(kind=RP), dimension(:,:,:,:),  allocatable :: S
      real(kind=RP), dimension(:,:,:,:),  allocatable :: U_x
      real(kind=RP), dimension(:,:,:,:),  allocatable :: U_y
      real(kind=RP), dimension(:,:,:,:),  allocatable :: U_z
      integer                                         :: NDOF              ! Number of degrees of freedom of element
      integer                                         :: firstIdx          ! Position in the global solution array
      logical                                         :: pointed = .TRUE.  ! .TRUE. (default) if Q and Qdot are pointed instead of allocated (needed for destruction since there's no other way to check this)
      type(Statistics_t)                              :: stats
#if defined(CAHNHILLIARD)
      real(kind=RP), dimension(:,:,:),   allocatable :: c   ! Cahn-Hilliard concentration
      real(kind=RP), dimension(:,:,:,:), allocatable :: gradC
      real(kind=RP), dimension(:,:,:),   allocatable :: mu  ! Cahn-Hilliard chemical pot.
#endif
      contains
         procedure   :: Construct => ElementStorage_Construct
         procedure   :: Destruct  => ElementStorage_Destruct
   end type ElementStorage_t

!  
!  Class for storing variables in the faces
!  ****************************************
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
         
         allocate ( self % Q   (NDOF) )
         allocate ( self % Qdot(NDOF) )
         
         if (bdf_order /= 0) then
            allocate ( self % PrevQ (NDOF,bdf_order) )
         end if
      end subroutine Storage_Construct
!
!/////////////////////////////////////////////////
!
      subroutine Storage_Destruct(self)
         implicit none
         !----------------------------------------------
         class(Storage_t)    :: self
         !----------------------------------------------
         
         safedeallocate(self % Q)
         safedeallocate(self % Qdot)
         safedeallocate(self % PrevQ)
      end subroutine Storage_Destruct
!
!///////////////////////////////////////////////////////////////////////////////////////////
!
!           Element Storage procedures
!           --------------------------
!
!///////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine ElementStorage_Construct(self, Nx, Ny, Nz, nEqn, nGradEqn, computeGradients,globalStorage,firstIdx)
         implicit none
         !------------------------------------------------------------
         class(ElementStorage_t)     :: self                               !<> Storage to be constructed
         integer        , intent(in) :: Nx, Ny, Nz                         !<  Polynomial orders in every direction
         integer        , intent(in) :: nEqn, nGradEqn                     !<  Number of equations and gradient equations
         logical        , intent(in) :: computeGradients                   !<  Compute gradients?
         type(Storage_t), intent(in), target, optional :: globalStorage    !<? Global storage to point to
         integer        , intent(in)        , optional :: firstIdx         !<? Position of the solution of the element in the global array
         !------------------------------------------------------------
         integer :: k, num_prevSol
         !------------------------------------------------------------
!
!        --------------------------------
!        Get number of degrees of freedom
!        --------------------------------
!
         self % NDOF = (Nx + 1) * (Ny + 1) * (Nz + 1) * nEqn
!
!        ----------------
!        Volume variables
!        ----------------
!
         if ( present(globalStorage) .and. present(firstIdx) ) then
            ! Solution and its derivative:
            self % Q   (1:nEqn,0:Nx,0:Ny,0:Nz) => globalStorage % Q   (firstIdx : firstIdx + self % NDOF-1)
            self % Qdot(1:nEqn,0:Nx,0:Ny,0:Nz) => globalStorage % Qdot(firstIdx : firstIdx + self % NDOF-1)
            
            ! Previous solution
            num_prevSol = size(globalStorage % PrevQ,2)
            allocate ( self % PrevQ(num_prevSol) )
            do k=1, num_prevSol
               self % PrevQ(k) % Q(1:nEqn,0:Nx,0:Ny,0:Nz) => globalStorage % PrevQ(firstIdx : firstIdx + self % NDOF-1,k)
            end do
            
            self % pointed = .TRUE.
         else
            ALLOCATE( self % Q   (nEqn,0:Nx,0:Ny,0:Nz) )
            ALLOCATE( self % QDot(nEqn,0:Nx,0:Ny,0:Nz) )
            self % pointed = .FALSE.
         end if
         
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
         self % mu = 0.0_RP
         self % c  = 0.0_RP
         self % gradC = 0.0_RP
#endif
      
         IF ( computeGradients )     THEN
            self % U_x         = 0.0_RP
            self % U_y         = 0.0_RP
            self % U_z         = 0.0_RP
         END IF

      end subroutine ElementStorage_Construct
!
!///////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine ElementStorage_Destruct(self)
         implicit none
         class(ElementStorage_t) :: self
         integer                 :: num_prevSol, k
         
         if (self % pointed) then
            nullify(self % Q)
            nullify(self % Qdot)
            num_prevSol = size(self % PrevQ)
            do k=1, num_prevSol
               nullify( self % PrevQ(k) % Q )
            end do
            safedeallocate(self % PrevQ)
         else
            deallocate(self % Q)
            deallocate(self % QDot)
         end if
         safedeallocate(self % G)
         safedeallocate(self % S)
         safedeallocate(self % U_x)
         safedeallocate(self % U_y)
         safedeallocate(self % U_z)

#if defined(CAHNHILLIARD)
         safedeallocate(self % mu)
         safedeallocate(self % c)
         safedeallocate(self % gradC)
#endif

         call self % stats % Destruct()

      end subroutine ElementStorage_Destruct
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
