!
!//////////////////////////////////////////////////////
!
!   @File:    StorageClass.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Thu Oct  5 09:17:17 2017
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module StorageClass
   use SMConstants
   implicit none

   private
   public   Storage_t

   type Statistics_t
      real(kind=RP), dimension(:,:,:,:),  allocatable    :: data
      contains
         procedure   :: Construct => Statistics_Construct
         procedure   :: Destruct  => Statistics_Destruct
   end type Statistics_t

   type Storage_t
      integer                                            :: N(3)
      integer                                            :: nEqn, nGradEqn
!
!     Element storage
!     ---------------
      real(kind=RP), dimension(:,:,:,:),  allocatable    :: Q
      real(kind=RP), dimension(:,:,:,:),  allocatable    :: QDot
      real(kind=RP), dimension(:,:,:,:),  allocatable    :: G
      real(kind=RP), dimension(:,:,:,:),  allocatable    :: S
      real(kind=RP), dimension(:,:,:,:),  allocatable    :: U_x
      real(kind=RP), dimension(:,:,:,:),  allocatable    :: U_y
      real(kind=RP), dimension(:,:,:,:),  allocatable    :: U_z
!
!     Face storage
!     ------------
      real(kind=RP), dimension(:,:,:,:), allocatable    :: Qb, Ub
      real(kind=RP), dimension(:,:,:,:), allocatable    :: U_xb, U_yb, U_zb
      real(kind=RP), dimension(:,:,:,:), allocatable    :: FStarb
!
!     Statistics
!     ----------
      type(Statistics_t)                               :: stats
      contains
         procedure   :: Construct => Storage_Construct
         procedure   :: Destruct  => Storage_Destruct
      
   end type Storage_t
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
      subroutine Storage_Construct(self, Nx, Ny, Nz, nEqn, nGradEqn, flowIsNavierStokes)
         implicit none
         class(Storage_t)     :: self
         integer, intent(in)  :: Nx, Ny, Nz
         integer, intent(in)  :: nEqn, nGradEqn
         logical              :: flowIsNavierStokes
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
         self % N = (/Nx, Ny, Nz/)
         self % nEqn = nEqn
         self % nGradEqn = nGradEqn
!
!        ----------------
!        Volume variables
!        ----------------
!
         ALLOCATE( self % Q   (nEqn,0:Nx,0:Ny,0:Nz) )
         ALLOCATE( self % QDot(nEqn,0:Nx,0:Ny,0:Nz) )
         ALLOCATE( self % G   (nEqn,0:Nx,0:Ny,0:Nz) )
         ALLOCATE( self % S   (nEqn,0:Nx,0:Ny,0:Nz) )
         
         IF ( flowIsNavierStokes )     THEN
            ALLOCATE( self % U_x(nGradEqn,0:Nx,0:Ny,0:Nz) )
            ALLOCATE( self % U_y(nGradEqn,0:Nx,0:Ny,0:Nz) )
            ALLOCATE( self % U_z(nGradEqn,0:Nx,0:Ny,0:Nz) )
         END IF
!
!        ---------------
!        Boundary values
!        ---------------
!
         ! Temporarily allocating with maximum (TODO: this is not very efficient and has to be changed) DGBoundaryStorage TYPE!!
         Nmax = MAX(Nx,Ny,Nz)
         ALLOCATE( self % Qb    (nEqn,0:Nmax,0:Nmax,6) )
         ALLOCATE( self % FStarb(nEqn,0:Nmax,0:Nmax,6) )
         
         IF ( flowIsNavierStokes )     THEN
            ALLOCATE( self % U_xb(nGradEqn,0:Nmax,0:Nmax,6) )
            ALLOCATE( self % U_yb(nGradEqn,0:Nmax,0:Nmax,6) )
            ALLOCATE( self % U_zb(nGradEqn,0:Nmax,0:Nmax,6) )
            ALLOCATE( self % Ub  (nGradEqn,0:Nmax,0:Nmax,6) )
         END IF
!
!        -----------------
!        Initialize memory
!        -----------------
!
         self % G           = 0.0_RP
         self % S           = 0.0_RP
         self % Q           = 0.0_RP
         self % QDot        = 0.0_RP
         self % Qb          = 0.0_RP
         self % FStarb      = 0.0_RP
      
         IF ( flowIsNavierStokes )     THEN
            self % Ub          = 0.0_RP
            self % U_x         = 0.0_RP
            self % U_y         = 0.0_RP
            self % U_z         = 0.0_RP
            self % U_xb        = 0.0_RP
            self % U_yb        = 0.0_RP
            self % U_zb        = 0.0_RP
         END IF

      end subroutine Storage_Construct

      subroutine Storage_Destruct(self)
         implicit none
         class(Storage_t)     :: self
   
         self % N = 0
         self % nEqn = 0
         self % nGradEqn = 0
         
         safedeallocate(self % Q)
         safedeallocate(self % QDot)
         safedeallocate(self % G)
         safedeallocate(self % S)
         safedeallocate(self % U_x)
         safedeallocate(self % U_y)
         safedeallocate(self % U_z)
         safedeallocate(self % Qb)
         safedeallocate(self % Ub)
         safedeallocate(self % U_xb)
         safedeallocate(self % U_yb)
         safedeallocate(self % U_zb)
         safedeallocate(self % FStarb)

         call self % stats % Destruct()

      end subroutine Storage_Destruct
!
!/////////////////////////////////////////////////////////////////////////////////////
!
!           Statistics procedures
!           ---------------------
!
!/////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Statistics_Construct(self, no_of_variables, storage)
         implicit none
         class(Statistics_t)           :: self
         integer,          intent(in)  :: no_of_variables
         class(Storage_t), intent(in)  :: storage
!
!        Allocate and initialize
!        -----------------------
         allocate( self % data(no_of_variables, 0:storage % N(1), 0:storage % N(2), 0:storage % N(3) ) ) 
         self % data = 0.0_RP

      end subroutine Statistics_Construct
   
      subroutine Statistics_Destruct(self)
         implicit none
         class(Statistics_t)     :: self

         safedeallocate( self % data )

      end subroutine Statistics_Destruct

end module StorageClass
