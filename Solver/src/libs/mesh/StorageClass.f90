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
      contains
         procedure   :: Construct => Storage_Construct
         procedure   :: Destruct  => Storage_Destruct
      
   end type Storage_t


   contains

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
         ALLOCATE( self % Q   (0:Nx,0:Ny,0:Nz,nEqn) )
         ALLOCATE( self % QDot(0:Nx,0:Ny,0:Nz,nEqn) )
         ALLOCATE( self % G   (0:Nx,0:Ny,0:Nz,nEqn) )
         ALLOCATE( self % S   (0:Nx,0:Ny,0:Nz,nEqn) )
         
         IF ( flowIsNavierStokes )     THEN
            ALLOCATE( self % U_x(0:Nx,0:Ny,0:Nz,nGradEqn) )
            ALLOCATE( self % U_y(0:Nx,0:Ny,0:Nz,nGradEqn) )
            ALLOCATE( self % U_z(0:Nx,0:Ny,0:Nz,nGradEqn) )
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
         
         checkdealloc(self % Q)
         checkdealloc(self % QDot)
         checkdealloc(self % G)
         checkdealloc(self % S)
         checkdealloc(self % U_x)
         checkdealloc(self % U_y)
         checkdealloc(self % U_z)
         checkdealloc(self % Qb)
         checkdealloc(self % Ub)
         checkdealloc(self % U_xb)
         checkdealloc(self % U_yb)
         checkdealloc(self % U_zb)
         checkdealloc(self % FStarb)

      end subroutine Storage_Destruct

end module StorageClass
