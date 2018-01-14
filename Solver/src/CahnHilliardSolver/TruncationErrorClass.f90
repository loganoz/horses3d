!
!//////////////////////////////////////////////////////
!
!   @File:    TruncationErrorClass.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Sun Jan 14 17:14:45 2018
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
module TruncationErrorClass
   use SMConstants
   use MultigridTypes
   use DGSEMClass
   
   implicit none
   
   private
   public TruncationError_t, EstimateTruncationError, InitializeForTauEstimation
   public NON_ISOLATED_TE, ISOLATED_TE
   
   !---------------------------------------------------------------------------------------------------------
   ! Type for storing the truncation error for one element in one direction according to the polynomial order
   !  -> Currently, we are only storing the infinity norm. This needs less storage and is indeed  a conservative 
   !     creterion |\tau| <= |\tau_1| + |\tau_2| (see paper "Truncation Error Estimation in the p-Anisotropic Discontinuous Galerkin Spectral Element Method")
   !---------------------------------------------------------------------------------------------------------
   type :: TruncErrorPol_t
      real(kind=RP), allocatable  :: maxTE(:)   ! |\tau|∞
      integer                     :: P          ! Polynomial order in this direction
   end type TruncErrorPol_t
   
   !-----------------------------------------------------------------------------
   ! Type for storing the truncation error of one element in the three directions 
   !-----------------------------------------------------------------------------
   type :: TruncationError_t
      type(TruncErrorPol_t)       :: Dir(3)
      
      contains
         procedure :: construct => ConstructTruncationError
         procedure :: destruct  => DestructTruncationError
   end type TruncationError_t
   
!
!  ----------------
!  Module variables
!  ----------------
!
   !
   !  Interface for the time derivative
   !  ---------------------------------
   abstract interface
      subroutine TimeDerivative_t(sem, t)
         use DGSEMClass
         type(DGSem)   :: sem
         real(kind=RP) :: t
      end subroutine TimeDerivative_t
   end interface
   
   procedure(TimeDerivative_t), pointer :: TimeDerivative
   
   !! Parameters
   integer, parameter :: ISOLATED_TE = 0
   integer, parameter :: NON_ISOLATED_TE = 1
   
!========
 contains
!========
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ROUTINES FOR TRUNCATION ERROR
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------------------------------------------------------------
!  Subroutine that constructs the variable that stores tau in one element
!  ----------------------------------------------------------------------
   subroutine ConstructTruncationError(TE,Nmin,Nmax)
      implicit none
      !------------------------------------------
      class(TruncationError_t), intent(inout) :: TE         !<> Variable that stores the truncation error
      integer                 , intent(in)    :: Nmin       !<  Minimum polynomial order for estimation
      integer                 , intent(in)    :: Nmax(3)    !<  Maximum polynomial order for estimation
      !------------------------------------------
      
      !------------------------------------------
      
      ! Allocate maxTE
      allocate (TE % Dir(1) % maxTE(NMin:Nmax(1)))
      allocate (TE % Dir(2) % maxTE(NMin:Nmax(2)))
      allocate (TE % Dir(3) % maxTE(NMin:Nmax(3)))
!
!     --------------------------
!     Initialize maxTE to zero 
!        Remarks:
!           a. This value will be overwritten when the TE is estimated.
!           b. If the polynomial order specified by the user is NMin, no TE estimation can be done and therefore the value TE = 0 is assumed.
!     --------------------------
!
      TE % Dir(1) % maxTE = 0._RP
      TE % Dir(2) % maxTE = 0._RP
      TE % Dir(3) % maxTE = 0._RP
   end subroutine ConstructTruncationError
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------------------------------------------------------------
!  Subroutine that destructs the variable that stores tau in one element
!  ----------------------------------------------------------------------
   subroutine DestructTruncationError(TE)
      implicit none
      !------------------------------------------
      class(TruncationError_t) :: TE     !<> Variable that stores the truncation error
      !------------------------------------------
      
      deallocate (TE % Dir(1) % maxTE)
      deallocate (TE % Dir(2) % maxTE)
      deallocate (TE % Dir(3) % maxTE)
      
   end subroutine DestructTruncationError
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------------------------------------------------------------
!  Subroutine that sets P for all elements in TE
!  ----------------------------------------------------------------------
   subroutine InitializeForTauEstimation(TE,sem,TruncErrorType)
      implicit none
      !------------------------------------------
      type(TruncationError_t) :: TE(:)
      type(DGSem), intent(in) :: sem
      integer    , intent(in) :: TruncErrorType !<  Either NON_ISOLATED_TE or ISOLATED_TE
      !------------------------------------------
      integer                 :: eID
      !------------------------------------------
      
      do eID=1, size(sem % mesh % elements)
         TE(eID) % Dir(1) % P = sem % Nx(eID)
         TE(eID) % Dir(2) % P = sem % Ny(eID)
         TE(eID) % Dir(3) % P = sem % Nz(eID)
      end do
      
      if (TruncErrorType == ISOLATED_TE) then
         error stop ':: Isolated truncation error is not yet implemented'
!~         TimeDerivative => ComputeTimeDerivativeIsol 
      else
         TimeDerivative => ComputeTimeDerivative
      end if
      
   end subroutine InitializeForTauEstimation
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ------------------------------------------------------------------------
!  Routine for computing the truncation error associated with every element
!  in a specified direction for a specific evaluation of the fine solution
!     (a-posteriori or quasio a-priori without correction)
!  ------------------------------------------------------------------------
   subroutine EstimateTruncationError(TE,sem,t,Var,Dir)
      use NodalStorageClass
      implicit none
      !-----------------------------
      type(TruncationError_t)  :: TE(:)       !<> Type containing the truncation error estimation
      type(DGSem)              :: sem           !<  sem  containing the Qdot after evaluating the fine solution in a given coarser mesh
      real(kind=RP)            :: t
      type(MGSolStorage_t)     :: Var(:)        !<  Type containing the source term in the mesh where the truncation error is being estimated
      integer                  :: Dir           !<  Direction in which the truncation error is being estimated
      !-----------------------------
      integer                  :: N(3)          !   Polynomial orders in element
      integer                  :: iEl           !   Element counter
      integer                  :: iEQ           !   Equation counter
      integer                  :: i,j,k         !   Coordinate index counters
      real(kind=RP)            :: wx, wy, wz    ! 
      real(kind=RP)            :: Jac
      real(kind=RP)            :: maxTE
      !-----------------------------
      
      call TimeDerivative(sem,t)
      
      do iEl = 1, size(sem % mesh % elements)
         N = sem % mesh % elements(iEl) % Nxyz
         
         if (N(Dir) >= TE(iEl) % Dir(Dir) % P) cycle ! it is actually never going to be greater than... but for security
         
         maxTE = 0._RP
         
         ! loop over all the degrees of freedom of the element
         do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1) ; do iEQ = 1, N_EQN
            wx  = NodalStorage(N(1)) % w (i)
            wy  = NodalStorage(N(2)) % w (j)
            wz  = NodalStorage(N(3)) % w (k)
            Jac = sem % mesh % elements(iEl) % geom % jacobian(i,j,k)
            
            maxTE =  MAX(maxTE , wx * wy * wz * Jac * ABS  (sem % mesh % elements(iEl) % storage % Qdot (iEQ,i,j,k) + &
                                                                                        Var(iEl) % Scase(iEQ,i,j,k))  )
            
         end do         ; end do         ; end do         ; end do
         
         TE(iEl) % Dir(Dir) % maxTE(N(Dir)) = maxTE
         
      end do
   end subroutine EstimateTruncationError
end module TruncationErrorClass
