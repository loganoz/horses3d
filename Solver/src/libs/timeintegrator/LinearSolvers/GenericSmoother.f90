!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!      GenericSmoother.f90
!      Created: Wed Nov 6 17:45:26 CET 2019
!      Version: 1.0 (Wed Nov 6 17:45:26 CET 2019)
!      Author: Wojciech Laskowski (wj.laskowski@upm.es).  
!
!      Generalized classs for handling all the different types of iterative solvers that can efficiently smooth out high-frequency errors. Routine has been initially implemented for p-Multigrid (MultigridSolverClass) to act as a smoother there. 
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
module GenericSmoother
   use SMConstants
   use MatrixClass
   ! use PhysicsStorage
   ! use PolynomialInterpAndDerivsModule
   ! use GaussQuadrature
   ! use DGSEMClass
   ! use TimeIntegratorDefinitions
   ! use NumericalJacobian      , only: NumJacobian_t
   ! use AnalyticalJacobian     , only: AnJacobian_t
   ! use JacobianComputerClass  , only: JacobianComputer_t, GetJacobianFlag
   implicit none
   
   private 
   public WeightedPointJacobi, ElementBlockJacobi, Smooth
   public S_NOTDEF, S_POINTJAC, S_BLOCKJAC, S_ILU, S_BILU
   public BJSmooth_t

   integer, parameter :: S_NOTDEF     = 0
   integer, parameter :: S_POINTJAC   = 1
   integer, parameter :: S_BLOCKJAC   = 2
   integer, parameter :: S_ILU        = 3
   integer, parameter :: S_BILU       = 4

   type :: BlockPrec_t
      real(KIND=RP), dimension(:,:), allocatable :: PLU        ! LU factorization of elemental preconditioner matrix
      integer      , dimension(:)  , allocatable :: LUpivots   ! LU pivots
   end type BlockPrec_t

   type :: BJSmooth_t
      class(DenseBlockDiagMatrix_t), pointer :: A_p           ! pointer to Block-Jacobian matrix
      type(BlockPrec_t), allocatable         :: BlockPrec(:)

   contains
      procedure                                  :: Construct => ConstructBlockJacobi
      procedure                                  :: Destruct  => DestructBlockJacobi
   end type BJSmooth_t

contains
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Smooth(A,x,b,n,niter,s_type,BJsmoother)
      implicit none
!-----Arguments---------------------------------------------------
      class(Matrix_t),  intent(inout)                           :: A                  ! Matrix to solve
      real(kind=RP),    intent(inout), dimension(:)             :: x                  ! Initial solution vector
      real(kind=RP),    intent(in),    dimension(:)             :: b                  ! Right hand side
      integer,          intent(in)                              :: n                  ! System size
      integer,          intent(in)                              :: niter              ! # of iterations
      integer,          intent(in)                              :: s_type             ! dummy smoother type
      type(BJSmooth_t), intent(in)                              :: BJSmoother         ! BJ smoother
!-----Local-Variables---------------------------------------------

      select case (s_type)
         case (S_NOTDEF)
            ERROR stop 'GenericSmoother :: Wrong smoother type'
         case (S_POINTJAC)
            select type (Acsr => A)
            type is (csrMat_t)
               call WeightedPointJacobi(Acsr,x,b,n,niter)
            end select 
         case (S_BLOCKJAC)
            select type (Acsr => A)
            type is (csrMat_t)
               call ElementBlockJacobi(Acsr,x,b,n,niter,BJSmoother)
            end select 
      end select 

   end subroutine Smooth
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine WeightedPointJacobi(A,x,b,n,SmoothIters)
      implicit none
!-----Arguments---------------------------------------------------
      class(csrMat_t), intent(inout)                           :: A                     ! Matrix to solve
      real(kind=rp),   intent(inout), dimension(:)             :: x                     ! Initial solution vector
      real(kind=rp),   intent(in),    dimension(:)             :: b                     ! Right hand side
      integer,         intent(in)                              :: n                     ! System siz
      integer,         intent(in)                              :: SmoothIters           ! # of iterations
!-----Local-Variables---------------------------------------------
      real(kind=rp)                           :: r(n) ! Residual
      real(kind=rp), parameter                :: w = 2._RP/3._RP  ! Weight (optimal for FD laplacian... but DGSEM?)
      integer                                 :: i,j              ! Counters
      ! real(kind=rp)                           :: rnorm            ! Residual norm
      !--------------------------------------------

      print *, " we're doing ", SmoothIters, " smoothings!"

      ! print *, "--------------"
      do i=1,SmoothIters

         r = CSR_MatVecMul( A, x ) ! CSR matrix product
         ! print *, "--------------"
         do j=1,n
            r(j) = b(j) - r(j)
            x(j) = x(j) + w * r(j) / A % Values(A % Diag(j))
            ! print *, j, A % Values(A % Diag(j))
            ! print *, j, x(j)
         end do
         ! print *, "--------------"

         ! rnorm = norm2(r)       ! Saves relative tolerance (one iteration behind)
         print*, "Iteration: ", i, " ; Norm: ", norm2(r) 
         
      end do
      ! print *, "--------------"

   end subroutine WeightedPointJacobi
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine ElementBlockJacobi(A,x,b,n,SmoothIters,BJSmoother)
      use DenseMatUtilities
      implicit none
!-----Arguments---------------------------------------------------
      class(csrMat_t),  intent(inout)                           :: A                     ! Matrix to solve
      real(kind=rp),    intent(inout), dimension(:)             :: x                     ! Initial solution vector
      real(kind=rp),    intent(in),    dimension(:)             :: b                     ! Right hand side
      integer,          intent(in)                              :: n                     ! System siz
      integer,          intent(in)                              :: SmoothIters           ! # of iterations
      type(BJSmooth_t), intent(in)                              :: BJSmoother            ! BJ smoother
!-----Local-Variables---------------------------------------------
      real(kind=rp)                           :: r(n) ! Residual
      integer                                 :: i,j, idx1, idx2 ! Counters
      real(kind=rp),    dimension(n)          :: xtmp            ! tmp solution vector
      ! real(kind=rp)                           :: rnorm            ! Residual norm
      !--------------------------------------------

      ! print *, " we're doing ", SmoothIters, " smoothings!"

      ! print *, "------- XXX --------"
      ! print *, x
      ! print *, "------- XXX --------"
      ! print *, "------- BBB --------"
      ! print *, b
      ! print *, "------- BBB --------"
      xtmp = 0.0_RP
      do i = 1, SmoothIters

         ! r = A*x
         r = CSR_MatVecMul( A, x ) ! CSR matrix product
         ! xtmp = x

         ! x = x + r/B 
! !$omp parallel do private(idx1,idx2) schedule(runtime)
         do j=1, BJSmoother % A_p % num_of_Blocks
            idx1 = BJSmoother % A_p % BlockIdx(j)
            idx2 = BJSmoother % A_p % BlockIdx(j+1)-1

            r(idx1:idx2) = b(idx1:idx2) - r(idx1:idx2)

            ! print *, "Here we go again"
            ! print *, BJSmoother % BlockPrec(j) % PLU
            ! print *, BJSmoother % BlockPrec(j) % LUpivots
            ! error stop "ELO"

            call SolveLU(ALU      = BJSmoother % BlockPrec(j) % PLU, &
                         LUpivots = BJSmoother % BlockPrec(j) % LUpivots, &
                         x = xtmp(idx1:idx2), &
                         b = r   (idx1:idx2))


            ! xtmp(idx1:idx2) =  matmul(BJSmoother % A_p % Blocks(j) % Matrix , r(idx1:idx2))
            
            x(idx1:idx2) = x(idx1:idx2) + xtmp(idx1:idx2)
         end do
! !$omp end parallel do
         ! print *, x


         ! rnorm = norm2(r)       ! Saves relative tolerance (one iteration behind)
         ! print*, "Iteration: ", i, " ; Norm: ", norm2(r) 
         ! r = CSR_MatVecMul( A, x ) ! CSR matrix product
         ! r = b - r
         ! print*, "Iteration: ", i, " ; Norm after: ", norm2(r) 
         ! error stop "HHHHHHHHHHHH"
      !    
      end do
      ! ! print *, "--------------"

   end subroutine ElementBlockJacobi
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine ConstructBlockJacobi(this,Amg)
      implicit none
!-----Arguments---------------------------------------------------
      class(BJSmooth_t), intent(inout)                  :: this                     ! Matrix to solve
      class(DenseBlockDiagMatrix_t), intent(in), target :: Amg                      ! Block-Jacobian matrix from MG classs
!-----Local-Variables---------------------------------------------
      integer :: k              ! Counters
      integer :: ndofelm        ! Dummies
      !--------------------------------------------

      this % A_p => Amg 
!
!     ------------------------------------------------
!     Allocate important variables for preconditioners
!     ------------------------------------------------
!
      allocate (this % BlockPrec(this % A_p % num_of_blocks))
      DO k = 1, this % A_p % num_of_blocks
         ndofelm = this % A_p % BlockSizes(k)
         print *, "counter is ", k, ". NDOF is ", ndofelm
         allocate (this % BlockPrec(k) % PLU(ndofelm,ndofelm) )
         allocate (this % BlockPrec(k) % LUpivots   (ndofelm) )
      end do
   end subroutine ConstructBlockJacobi
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine DestructBlockJacobi(this)
      implicit none
!-----Arguments---------------------------------------------------
      class(BJSmooth_t), intent(inout)             :: this                     ! Matrix to solve
!-----Local-Variables---------------------------------------------
      integer                                 :: i,j              ! Counters
      !--------------------------------------------
   end subroutine DestructBlockJacobi
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
end module GenericSmoother