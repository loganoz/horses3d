#include "Includes.h"
MODULE CuSparseSolverClass
   USE GenericLinSolverClass
   USE CSRMatrixClass            , only: csrMat_t
   use PETScMatrixClass
   USE SMConstants
   use DGSEMClass
   use TimeIntegratorDefinitions
   use Utilities                 , only: AlmostEqual
   use StopwatchClass            , only: Stopwatch
   use MPI_Process_Info          , only: MPI_Process

#ifdef HAS_CUDA
   ! use cublas
   ! use cusolverDn
   ! use cudafor
   use cusparse
#endif
   use iso_c_binding

   implicit none
   
   TYPE, EXTENDS(GenericLinSolver_t) :: CuSparseSolver_t
      type(csrMat_t)                             :: A                                  ! Jacobian matrix
      type(csrMat_t), pointer                    :: ALU                                ! LU-Factorized Jacobian matrix
      real(kind=RP), DIMENSION(:), ALLOCATABLE   :: x                                  ! Solution vector
      real(kind=RP), DIMENSION(:), ALLOCATABLE   :: b                                  ! Right hand side
      real(kind=RP)                              :: Ashift
      LOGICAL                                    :: AIsPrealloc
      logical                                    :: Variable_dt                        ! Is the time-step variable?
   CONTAINS
      !Subroutines:
      PROCEDURE :: construct                    => ConstructCuSparseSolver
      procedure :: ComputeAndFactorizeJacobian  => CS_ComputeAndFactorizeJacobian
      procedure :: ReFactorizeJacobian          => CS_ReFactorizeJacobian
      PROCEDURE :: solve
      procedure :: SolveLUDirect                => CS_SolveLUDirect
      procedure :: SetRHSValue                  => CS_SetRHSValue
      procedure :: SetRHS                       => CS_SetRHS
      PROCEDURE :: GetXValue                    => CS_GetXValue
      PROCEDURE :: GetX                         => CS_GetX
      PROCEDURE :: destroy                      => CS_destroy
      PROCEDURE :: SetOperatorDt
      PROCEDURE :: ReSetOperatorDt
      procedure :: CS_ComputeJacobian
      procedure :: FactorizeJacobian            => CS_FactorizeJacobian
      procedure :: SetJacobian                  => CS_SetJacobian
      !Functions:
      PROCEDURE :: Getxnorm                     => CS_GetXnorm    !Get solution norm
      PROCEDURE :: Getrnorm    !Get residual norm
   END TYPE CuSparseSolver_t
   
   private
   public CuSparseSolver_t, GenericLinSolver_t
   
!
!  Useful interfaces
!  -----------------

   ! interface

   !    ! cudaMalloc
   !    integer (c_int) function cudaMalloc ( buffer, size ) bind (C, name="cudaMalloc" ) 
   !       use iso_c_binding
   !       implicit none
   !       type (c_ptr)  :: buffer
   !       integer (c_size_t), value :: size
   !    end function cudaMalloc

   ! end interface

      ! ! cudaMemcpy 
      ! integer (c_int) function cudaMemcpy ( dst, src, count, kind ) bind (C, name="cudaMemcpy" )
      !    ! note: cudaMemcpyHostToDevice = 1
      !    ! note: cudaMemcpyDeviceToHost = 2
      !    use iso_c_binding
      !    type (C_PTR), value :: dst, src
      !    integer (c_size_t), value :: count, kind
      ! end function cudaMemcpy

      ! ! cudaFree
      ! integer (c_int) function cudaFree(buffer)  bind(C, name="cudaFree")
      !    use iso_c_binding
      !    implicit none
      !    type (C_PTR), value :: buffer
      ! end function cudaFree

      ! integer (c_int) function cudaMemGetInfo(fre, tot) bind(C, name="cudaMemGetInfo")
      !    use iso_c_binding
      !    implicit none
      !    type(c_ptr),value :: fre
      !    type(c_ptr),value :: tot
      ! end function cudaMemGetInfo

      ! integer(c_int) function cusolverDnCreate(cusolver_Hndl) bind(C,name="cusolverDnCreate")      
      !    use iso_c_binding
      !    implicit none
      !    type(c_ptr)::cusolver_Hndl
         
      ! end function
         
      ! integer(c_int) function cusolverDnDestroy(cusolver_Hndl) bind(C,name="cusolverDnDestroy")
      !    use iso_c_binding
      !    implicit none
      !    type(c_ptr),value::cusolver_Hndl
      
      ! end function

      ! integer(c_int) function cusolverDnSgetrf_bufferSize(cusolver_Hndl,m,n,d_A,lda,Lwork) bind(C,name="cusolverDnSgetrf_bufferSize") 
      !    use iso_c_binding
      !    implicit none
      
      !    type(c_ptr),value::cusolver_Hndl
      !    integer(c_int),value::m
      !    integer(c_int),value::n
      !    type(c_ptr),value::d_A
      !    integer(c_int),value::lda
      !    type(c_ptr),value::Lwork
      ! end function

      ! integer(c_int) function cusolverDnSgetrf(cusolver_Hndl,m,n,d_A,lda,d_WS,d_Ipiv,d_devInfo) bind(C, name="cusolverDnSgetrf")
      !    use iso_c_binding
      !    implicit none

      !    type(c_ptr),value::cusolver_Hndl
      !    integer(c_int),value::m
      !    integer(c_int),value::n
      !    type(c_ptr),value::d_A
      !    integer(c_int),value::lda
      !    type(c_ptr),value::d_WS
      !    type(c_ptr),value::d_Ipiv
      !    type(c_ptr),value::d_devInfo
      
      ! end function 

      ! integer (c_int) function cusolverDnSgetrs(cusolver_Hndl,trans,n,nrhs,d_A,lda,d_Ipiv,d_B,ldb,d_devInfo) bind(C, name="cusolverDnSgetrs")
      !    use iso_c_binding
      !    implicit none

      !    type(c_ptr),value::cusolver_Hndl
      !    integer(c_int), value::trans
      !    integer(c_int), value::n
      !    integer(c_int), value::nrhs
      !    type(c_ptr),value::d_A
      !    integer(c_int), value::lda    
      !    type(c_ptr),value::d_Ipiv
      !    type(c_ptr),value::d_B
      !    integer(c_int),value::ldb
      !    type(c_ptr),value::d_devInfo
         
      ! end function
   
!========
 CONTAINS
!========
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------
!  Direct solver constructor
!  -----------------------
   subroutine ConstructCuSparseSolver(this,DimPrb, globalDimPrb, nEqn,controlVariables,sem,MatrixShiftFunc)
      implicit none
      !-----------------------------------------------------------
      class(CuSparseSolver_t), intent(inout), TARGET :: this
      integer                  , intent(in)            :: DimPrb
      integer                  , intent(in)            :: globalDimPrb
      integer                  , intent(in)            :: nEqn
      TYPE(FTValueDictionary)  , intent(in), OPTIONAL  :: controlVariables
      TYPE(DGSem), TARGET                  , OPTIONAL  :: sem
      procedure(MatrixShift_FCN)                       :: MatrixShiftFunc
      !-----------------------------------------------------------
      integer :: A_mem_stat
      type(c_ptr) :: d_A, d_B
      type(c_ptr) :: d_C
      integer :: A_size, status
      integer, target :: tmp_i_1

! #ifdef HAS_CUDA
!       type(cusparseHandle) :: h
! #endif
      !-----------------------------------------------------------

      tmp_i_1 = 0
      d_C = c_loc(tmp_i_1)

! #ifdef HAS_CUDA
!       status = cusparseCreate(h)
! #endif

      call this % GenericLinSolver_t % construct(DimPrb, globalDimPrb, nEqn,controlVariables,sem,MatrixShiftFunc)
      
      if (MPI_Process % doMPIRootAction) then
         ERROR stop 'CuSparseSolver_t cannot be used as a distributed solver'
      end if
      
      if ( present(sem) ) then
         this % p_sem => sem
      end if
      
      MatrixShift => MatrixShiftFunc
      
      this % DimPrb = DimPrb
      
      allocate(this % x(DimPrb))
      allocate(this % b(DimPrb))      

      call this % A % construct(num_of_Rows = DimPrb, withMPI = .false.)

!
!     Configure the Jacobian (this includes matrix preallocation -if needed)
!     ----------------------------------------------------------------------
      if ( present(sem) ) then
         call this % Jacobian % Configure (sem % mesh, nEqn, this % A)
      end if
      
!
!     Set variables from input file
!     -----------------------------
!
      if ( present(controlVariables) ) then
!
!        See if the time-step is constant
!           If it's not, the factorized matrix cannot be stored in the matrix structure
!        ------------------------------------------------------------------------------
         if ( controlVariables % containsKey("dt") ) then
            this % Variable_dt = .FALSE.
            this % ALU => this % A
         else
            this % Variable_dt = .TRUE.
            allocate (this % ALU)
         end if
         
      else
         this % Variable_dt = .TRUE.
         allocate (this % ALU)
      end if
      
      call Stopwatch % CreateNewEvent("Sparse LU-Factorization")
      
   end subroutine ConstructCuSparseSolver
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine CS_SetRHSValue(this, irow, value)
      implicit none
      !-----------------------------------------------------------
      class(CuSparseSolver_t), intent(inout) :: this
      INTEGER                  , intent(in)    :: irow
      real(kind=RP)            , intent(in)    :: value
      !-----------------------------------------------------------
      
      this % b (irow) = value
      
   end subroutine CS_SetRHSValue
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine CS_SetRHSValues(this, nvalues, irow, values)
      class(CuSparseSolver_t)   , intent(inout)     :: this
      INTEGER                     , intent(in)        :: nvalues
      INTEGER      , DIMENSION(1:), intent(in)        :: irow
      real(kind=RP), DIMENSION(1:), intent(in)        :: values
      !------------------------------------------------------
      integer                                        :: i
      
      do i=1, nvalues
         if (irow(i)<0) cycle
         this % b(irow(i)) = values(i)
      end do
      
   end subroutine CS_SetRHSValues
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine CS_SetRHS(this, RHS)
      implicit none
      !-arguments-----------------------------------------------------------
      class(CuSparseSolver_t), intent(inout) :: this
      real(kind=RP)            , intent(in)    :: RHS(this % DimPrb)
      !---------------------------------------------------------------------
      
      this % b = RHS
      
   end subroutine CS_SetRHS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine solve(this, nEqn, nGradEqn, ComputeTimeDerivative,tol,maxiter,time,dt,ComputeA) 
      implicit none
!
!     ----------------------------------------------------
!     Main subroutine for solving system using direct solver
!     ----------------------------------------------------
!
      !-----------------------------------------------------------
      class(CuSparseSolver_t), target, intent(inout) :: this
      integer,       intent(in)                :: nEqn
      integer,       intent(in)                :: nGradEqn
      procedure(ComputeTimeDerivative_f)               :: ComputeTimeDerivative
      real(kind=RP), OPTIONAL                  :: tol
      INTEGER      , OPTIONAL                  :: maxiter
      real(kind=RP), OPTIONAL                  :: time
      real(kind=RP), OPTIONAL                  :: dt
      logical      , optional      , intent(inout) :: ComputeA
      !-----------------------------------------------------------
      integer                                  :: error
      !-----------------------------------------------------------
      
!
!     Compute Jacobian matrix if needed
!     -----------------------------------------------------
      
      if ( present(ComputeA)) then
         if (ComputeA) then
            call this % CS_ComputeJacobian(dt,time,nEqn,nGradEqn,ComputeTimeDerivative)
            ComputeA = .FALSE.
         end if
      else
         call this % CS_ComputeJacobian(dt,time,nEqn,nGradEqn,ComputeTimeDerivative)
      end if
      
      call this % SolveLUDirect(error)

      if (error .NE. 0) THEN
         WRITE(*,*) 'Direct Solver ERROR:', error
         stop
      else
         this % converged = .TRUE.
      end if
      
   end subroutine solve
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine CS_ComputeJacobian(this,dt,time,nEqn,nGradEqn,ComputeTimeDerivative)
      use DenseBlockDiagonalMatrixClass
      implicit none
      !-----------------------------------------------------------
      class(CuSparseSolver_t), intent(inout) :: this
      real(kind=RP), intent(in)                :: dt
      real(kind=RP), intent(in)                :: time
      integer,       intent(in)                :: nEqn
      integer,       intent(in)                :: nGradEqn
      procedure(ComputeTimeDerivative_f)       :: ComputeTimeDerivative
      !-----------------------------------------------------------
      type(csrMat_t) :: B, Cmat !debug
      

      call this % Jacobian % Compute (this % p_sem, nEqn, time, this % A, ComputeTimeDerivative, ComputeTimeDerivative)
      call this % SetOperatorDt(dt)
      call this % FactorizeJacobian

!~         !<debug
         
!~         call this % A % Visualize('AnJac_visu.dat')
         
!~         !------------
!~         this % JacobianComputation = NUMERICAL_JACOBIAN
!~         call B % construct(num_of_Rows = this % DimPrb, withMPI = .false.)
         
         
!~         call this % ComputeJacobian(B,time,nEqn,nGradEqn,ComputeTimeDerivative)
         
!~         call B % Visualize('NumJac_visu.dat')
         
!~         !------------
!~         call  this % A % MatAdd(B, Cmat,-1._RP)
         
!~         print*, 'Error(L2)  = ',  norm2( (Cmat % Values) )
!~         print*, 'Error(inf) = ',  maxval( abs(Cmat % Values) )
!~         print*, '       pos = ',  maxloc( abs(Cmat % Values) ), Cmat % Cols(maxloc( abs(Cmat % Values) ))
         
!~         call Cmat % Visualize('C_visu.dat')
         
!~         stop
         
         
!~         !debug>
      
   end subroutine CS_ComputeJacobian
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine CS_GetXValue(this,irow,x_i)       
      implicit none
      !-----------------------------------------------------------
      class(CuSparseSolver_t), intent(inout) :: this
      INTEGER                  , intent(in)    :: irow
      real(kind=RP)            , INTENT(OUT)   :: x_i
      !-----------------------------------------------------------
      
      x_i = this % x(irow)
      
   end subroutine CS_GetXValue
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function CS_GetX(this) result(x)
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(CuSparseSolver_t), INTENT(INOUT) :: this
      REAL(KIND=RP)                            :: x(this % DimPrb)
      !-----------------------------------------------------------
      
      x = this % x
      
   end function CS_GetX
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine CS_destroy(this)       
      implicit none
      !-----------------------------------------------------------
      class(CuSparseSolver_t), intent(inout) :: this
      !-----------------------------------------------------------
      
      call this % A % destruct
      
      DEALLOCATE(this % b)
      DEALLOCATE(this % x)
      
      if (this % Variable_dt) then
         call this % ALU % destruct
         deallocate (this % ALU)
      else
         nullify (this % ALU)
      end if
      
      this % AIsPrealloc = .FALSE.
      
    end subroutine CS_destroy
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SetOperatorDt(this,dt)       
      implicit none
      !-----------------------------------------------------------
      class(CuSparseSolver_t), intent(inout) :: this
      real(kind=RP)            , intent(in)    :: dt
      !-----------------------------------------------------------
      
      this % Ashift = MatrixShift(dt)

      call this % A % Shift(this % Ashift)
      
    end subroutine SetOperatorDt
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine ReSetOperatorDt(this,dt)       
      implicit none
      !-----------------------------------------------------------
      class(CuSparseSolver_t), intent(inout) :: this
      real(kind=RP)            , intent(in)    :: dt
      !-----------------------------------------------------------
      real(kind=RP)                            :: shift
      !-----------------------------------------------------------
      
      shift = MatrixShift(dt)
      if ( AlmostEqual(shift,this % Ashift) ) return
      call this % A % Shift(-this % Ashift)
      call this % A % Shift(shift)
      this % Ashift = shift
      call this % FactorizeJacobian
      
    end subroutine ReSetOperatorDt
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function CS_GetXnorm(this,TypeOfNorm) result(xnorm)
      implicit none
      !-----------------------------------------------------------
      class(CuSparseSolver_t), intent(inout) :: this
      CHARACTER(len=*)                         :: TypeOfNorm
      real(kind=RP)                            :: xnorm
      !-----------------------------------------------------------
      
      select case (TypeOfNorm)
         CASE ('infinity')
            xnorm = MAXVAL(ABS(this % x))
         CASE ('l2')
            xnorm = NORM2(this % x)
         CASE DEFAULT
            stop 'MKLPardisoSolverClass ERROR: Norm not implemented yet'
      end select
   end function CS_GetXnorm
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function Getrnorm(this) result(rnorm)
      implicit none
!
!     ----------------------------------------
!     Currently implemented with infinity norm
!     ----------------------------------------
!
      !-----------------------------------------------------------
      class(CuSparseSolver_t), intent(inout) :: this
      real(kind=RP)                            :: rnorm
      !-----------------------------------------------------------
      real(kind=RP)                            :: residual(this % DimPrb)
      !-----------------------------------------------------------
      
      residual = this % b - this % A % MatVecMul (this % x)
      rnorm = MAXVAL(ABS(residual))
      
      
      !rnorm = NORM2(this % x)
      
   end function Getrnorm
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!

   subroutine CS_ComputeAndFactorizeJacobian(self,nEqn, nGradEqn, F_J, dt, eps, mode_in)
!
!     *************************************************************************************
!     This subroutine performs the following:
!           -> Construct the numerical Jacobian of the function F_J, with perturbation eps.
!           -> Shifts the Jacobian -dt * I.
!           -> Converts the Jacobian to CSR.
!           -> Factorizes the Jacobian (jobs 12 in Pardiso).
!
!     *************************************************************************************
!
      implicit none
      class(CuSparseSolver_t), intent(inout) :: self
      integer,                   intent(in)    :: nEqn, nGradEqn
      procedure(ComputeTimeDerivative_f)       :: F_J
      real(kind=RP), intent(in)                :: dt
      real(kind=RP), intent(in)                :: eps
      integer,       intent(in)                :: mode_in



      call self % Jacobian % Compute (self % p_sem, nEqn, 0._RP, self % A, F_J, F_J, eps, .false., mode_in)
      call self % SetOperatorDt(dt)
      self % A % values = -dt * self % A % values

      call self % FactorizeJacobian
      
   end subroutine CS_ComputeAndFactorizeJacobian

   subroutine CS_ReFactorizeJacobian(self) 
! 
!     ************************************************************************************* 
!        This subroutine changes the gamma0 coefficient from 1.0 to 1.5 by 
!        adding 0.5 to the diagonal. 
!     ************************************************************************************* 
! 
      implicit none 
      class(CuSparseSolver_t), intent(inout) :: self 
! 
!     --------------- 
!     Local variables 
!     --------------- 
! 
      integer     :: error 
! 
!     Shift the Jacobian 
!     ------------------ 
      self % A % values(self % A % diag) = self % A % values(self % A % diag) + 0.5_RP 

! 
!     Perform the factorization 
!     ------------------------- 
    ! TBD
      ! #ifdef HAS_MKL 
!       call pardiso(self % Pardiso_pt, 1, 1, self % mtype, 12, self % ALU % num_of_Rows, self % ALU % values, &
!                    self % ALU % rows, self % ALU % cols, self % perm, 1, self % Pardiso_iparm, 0, &
!                    self % b, self % x, error)
! #else 
!       STOP 'MKL not linked correctly' 
! #endif 
 
   end subroutine CS_ReFactorizeJacobian 

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine CS_FactorizeJacobian(self)
      implicit none
      !-arguments---------------------------------------------------
      class(CuSparseSolver_t), intent(inout)  :: self
      !-local-variables---------------------------------------------
      integer     :: error
      !-------------------------------------------------------------
#ifdef HAS_CUDA
      ! type(cusparseHandle) :: h
#endif
      !-----------------------------------------------------------

#ifdef HAS_CUDA
      ! status = cusparseCreate(h)
      print *, "Factorize later ---"
#else
      error stop "Only CUDA routine"
#endif

      call Stopwatch % Start("Sparse LU-Factorization")

      ! TBD
      
!       if (self % Variable_dt) then
!          call self % ALU % destruct
!          call self % ALU % constructWithCSRArrays  (self % A % Rows, &
!                                                  self % A % Cols, &
!                                                  self % A % Values)
!       end if
      
! !
! !     Perform the factorization
! !     -------------------------
! #ifdef HAS_MKL
!       call pardiso(self % Pardiso_pt, 1, 1, self % mtype, 12, self % ALU % num_of_Rows, self % ALU % values, &
!                    self % ALU % rows, self % ALU % cols, self % perm, 1, self % Pardiso_iparm, 0, &
!                    self % b, self % x, error)
! #else
!       stop 'MKL not linked correctly'
! #endif
      
      call Stopwatch % Pause("Sparse LU-Factorization")
   end subroutine CS_FactorizeJacobian
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine CS_SolveLUDirect(self,error_out)
      implicit none
      class(CuSparseSolver_t), intent(inout)  :: self
      integer, optional        , intent(out)    :: error_out
!
!     ---------------
!     Local variables
!     ---------------
!
      integer     :: error
      integer :: status, cusparse_create 

#ifdef HAS_CUDA
      type(cusparseHandle) :: handle
#endif

#ifdef HAS_CUDA
!       cusolverSpDcsrlsvlu

!       cusolverStatus_t 
! cusolverSpDcsrlsvlu[Host](cusolverSpHandle_t handle,
!                  int n,
!                  int nnzA,
!                  const cusparseMatDescr_t descrA,
!                  const double *csrValA,
!                  const int *csrRowPtrA,
!                  const int *csrColIndA,
!                  const double *b,
!                  double tol,
!                  int reorder,
!                  double *x,
!                  int *singularity);]

      status = cusparseCreate(handle)
#else
      error stop "Only CUDA routine"
#endif
error stop "TBD"

      ! TBD
! #ifdef HAS_MKL
!       call pardiso(self % Pardiso_pt, 1, 1, self % mtype, 33, self % ALU % num_of_Rows, self % ALU % values, &
!                    self % ALU % rows, self % ALU % cols, self % perm, 1, self % Pardiso_iparm, 0, &
!                    self % b, self % x, error)
      
!       if ( present(error_out) ) then
!          error_out = error
!       end if
      
! #else
!       stop 'MKL not linked correctly'
! #endif

   end subroutine CS_SolveLUDirect
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine CS_SetJacobian(this,Matrix)
      implicit none
      !-arguments-----------------------------------------------------------
      class(CuSparseSolver_t), intent(inout)  :: this
      class(Matrix_t)          , intent(in)     :: Matrix
      !---------------------------------------------------------------------
      
      select type(Matrix)
      class is(csrMat_t)
         this % A = Matrix
      class default
         ERROR stop 'CS_SetJacobian :: Wrong matrix type'
      end select
      
   end subroutine CS_SetJacobian
END MODULE CuSparseSolverClass