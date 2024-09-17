!//////////////////////////////////////////////////////
!
!     Class for solving linear systems using MKL version of Pardiso
!     -> It is possible to construct the matrix using PETSc. This option is currently deactivated... deprecate??
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "Includes.h"
#ifdef HAS_PETSC
#include "petsc/finclude/petsc.h"
#endif
MODULE MKLPardisoSolverClass
   USE GenericLinSolverClass
   USE CSRMatrixClass            , only: csrMat_t
   use PETScMatrixClass
   USE SMConstants
   use DGSEMClass
   use TimeIntegratorDefinitions
   use Utilities                 , only: AlmostEqual
   use StopwatchClass            , only: Stopwatch
   use MPI_Process_Info          , only: MPI_Process
#ifdef HAS_PETSC
   use petsc
#endif
   implicit none
   
   TYPE, EXTENDS(GenericLinSolver_t) :: MKLPardisoSolver_t
      type(csrMat_t)                             :: A                                  ! Jacobian matrix
      type(csrMat_t), pointer                    :: ALU                                ! LU-Factorized Jacobian matrix
      type(PETSCMatrix_t)                        :: PETScA
      real(kind=RP), DIMENSION(:), ALLOCATABLE   :: x                                  ! Solution vector
      real(kind=RP), DIMENSION(:), ALLOCATABLE   :: b                                  ! Right hand side
      real(kind=RP)                              :: Ashift
      LOGICAL                                    :: AIsPrealloc
      logical                                    :: Variable_dt                        ! Is the time-step variable?
      
      !Variables for creating Jacobian in PETSc context:
      LOGICAL                                    :: AIsPetsc = .false.
      
      !Variables directly related with mkl pardiso solver
      INTEGER                                    :: mtype                              ! Matrix type. See construct
      INTEGER, ALLOCATABLE                       :: perm(:)
      INTEGER, POINTER                           :: Pardiso_iparm(:) => NULL()         ! Parameters for mkl version of pardiso
      INTEGER(KIND=AddrInt), POINTER             :: Pardiso_pt(:)    => NULL()  
   CONTAINS
      !Subroutines:
      PROCEDURE :: construct                    => ConstructMKLContext
      procedure :: ComputeAndFactorizeJacobian  => MKL_ComputeAndFactorizeJacobian
      procedure :: ReFactorizeJacobian          => MKL_ReFactorizeJacobian
      PROCEDURE :: solve
      procedure :: SolveLUDirect                => MKL_SolveLUDirect
      procedure :: SetRHSValue                  => MKL_SetRHSValue
      procedure :: SetRHS                       => MKL_SetRHS
      PROCEDURE :: GetXValue                    => MKL_GetXValue
      PROCEDURE :: GetX                         => MKL_GetX
      PROCEDURE :: destroy                      => MKL_destroy
      PROCEDURE :: SetOperatorDt
      PROCEDURE :: ReSetOperatorDt
      procedure :: ComputeJacobianMKL
      procedure :: FactorizeJacobian            => MKL_FactorizeJacobian
      procedure :: SetJacobian                  => MKL_SetJacobian
      !Functions:
      PROCEDURE :: Getxnorm                     => MKL_GetXnorm    !Get solution norm
      PROCEDURE :: Getrnorm    !Get residual norm
   END TYPE MKLPardisoSolver_t
   
   private
   public MKLPardisoSolver_t, GenericLinSolver_t
   
!
!  Useful interfaces
!  -----------------
   interface
      subroutine pardisoinit(pt, mtype, iparm)
         use SMConstants
         implicit none
         integer(kind=AddrInt) :: pt(*)
         integer               :: mtype
         integer               :: iparm(*)
      end subroutine pardisoinit
      
      subroutine pardiso(pt, maxfct, mnum, mtype, phase, n, &
                        values, rows, cols, perm, nrhs, iparm, msglvl, b, x, ierror)
         use SMConstants
         real(kind=RP)           :: values(*), b(*), x(*)
         integer(kind=AddrInt)   :: pt(*)
         integer                 :: perm(*), nrhs, iparm(*), msglvl, ierror
         integer                 :: maxfct, mnum, mtype, phase, n, rows(*), cols(*)
      end subroutine pardiso
   end interface
   
!========
 CONTAINS
!========
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------
!  MKL pardiso constructor
!  -----------------------
   subroutine ConstructMKLContext(this,DimPrb, globalDimPrb, nEqn,controlVariables,sem,MatrixShiftFunc)
      implicit none
      !-----------------------------------------------------------
      class(MKLPardisoSolver_t), intent(inout), TARGET :: this
      integer                  , intent(in)            :: DimPrb
      integer                  , intent(in)            :: globalDimPrb
      integer                  , intent(in)            :: nEqn
      TYPE(FTValueDictionary)  , intent(in), OPTIONAL  :: controlVariables
      TYPE(DGSem), TARGET                  , OPTIONAL  :: sem
      procedure(MatrixShift_FCN)                       :: MatrixShiftFunc
      !-----------------------------------------------------------
#ifdef HAS_PETSC
      PetscErrorCode :: ierr
#endif
      !-----------------------------------------------------------
      
      call this % GenericLinSolver_t % construct(DimPrb, globalDimPrb, nEqn,controlVariables,sem,MatrixShiftFunc)
      
      if (MPI_Process % doMPIRootAction) then
         error stop 'MKLPardisoSolver_t cannot be used as a distributed solver'
         !TODO: Implement cluster_sparse_solver (MKL) or use the actual pardiso solver (http://www.pardiso-project.org/)
      end if
      
      if ( present(sem) ) then
         this % p_sem => sem
      end if
      
      MatrixShift => MatrixShiftFunc
      
      this % DimPrb = DimPrb
      
      allocate(this % x(DimPrb))
      allocate(this % b(DimPrb))
      
      this % mtype = 11 !Set matrix type to real unsymmetric (change?)
    
      ALLOCATE(this % Pardiso_pt(64))
      ALLOCATE(this % Pardiso_iparm(64))
      
      ALLOCATE(this % perm(DimPrb))
      this % perm = 0
      
      if(this % AIsPetsc) then
#ifdef HAS_PETSC
         call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
#else
         error stop "MKL-Pardiso needs PETSc for constructung the Jacobian Matrix"
#endif
         call this % PETScA % construct (num_of_Rows = DimPrb, withMPI = .FALSE.)
      else
         call this % A % construct(num_of_Rows = DimPrb, withMPI = .false.)
         
      end if
!
!     Configure the Jacobian (this includes matrix preallocation -if needed)
!     ----------------------------------------------------------------------
      if ( present(sem) ) then
         call this % Jacobian % Configure (sem % mesh, nEqn, this % A)
      end if

#ifdef HAS_MKL
      call pardisoinit(this % Pardiso_pt, this % mtype, this % Pardiso_iparm)
#else
      error stop 'MKL not linked correctly'
#endif
      
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
      
   end subroutine ConstructMKLContext
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MKL_SetRHSValue(this, irow, value)
      implicit none
      !-----------------------------------------------------------
      class(MKLPardisoSolver_t), intent(inout) :: this
      INTEGER                  , intent(in)    :: irow
      real(kind=RP)            , intent(in)    :: value
      !-----------------------------------------------------------
      
      this % b (irow) = value
      
   end subroutine MKL_SetRHSValue
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MKL_SetRHSValues(this, nvalues, irow, values)
      class(MKLPardisoSolver_t)   , intent(inout)     :: this
      INTEGER                     , intent(in)        :: nvalues
      INTEGER      , DIMENSION(1:), intent(in)        :: irow
      real(kind=RP), DIMENSION(1:), intent(in)        :: values
      !------------------------------------------------------
      integer                                        :: i
      
      do i=1, nvalues
         if (irow(i)<0) cycle
         this % b(irow(i)) = values(i)
      end do
      
   end subroutine MKL_SetRHSValues
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MKL_SetRHS(this, RHS)
      implicit none
      !-arguments-----------------------------------------------------------
      class(MKLPardisoSolver_t), intent(inout) :: this
      real(kind=RP)            , intent(in)    :: RHS(this % DimPrb)
      !---------------------------------------------------------------------
      
      this % b = RHS
      
   end subroutine MKL_SetRHS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine solve(this, nEqn, nGradEqn, ComputeTimeDerivative,tol,maxiter,time,dt,ComputeA) 
      implicit none
!
!     ----------------------------------------------------
!     Main subroutine for solving system using mkl pardiso
!     ----------------------------------------------------
!
      !-----------------------------------------------------------
      class(MKLPardisoSolver_t), target, intent(inout) :: this
      integer,       intent(in)                :: nEqn
      integer,       intent(in)                :: nGradEqn
      procedure(ComputeTimeDerivative_f)               :: ComputeTimeDerivative
      real(kind=RP), OPTIONAL                  :: tol
      INTEGER      , OPTIONAL                  :: maxiter
      real(kind=RP), OPTIONAL                  :: time
      real(kind=RP), OPTIONAL                  :: dt
      logical      , optional      , intent(inout) :: ComputeA
      !-----------------------------------------------------------
#ifdef HAS_MKL
      integer                                  :: error
      !-----------------------------------------------------------
      
!
!     Compute Jacobian matrix if needed
!        (done in petsc format and then transformed to CSR since the CSR cannot be filled by the Jacobian calculators)
!     -----------------------------------------------------
      
      if ( present(ComputeA)) then
         if (ComputeA) then
            call this % ComputeJacobianMKL(dt,time,nEqn,nGradEqn,ComputeTimeDerivative)
            ComputeA = .FALSE.
         end if
      else
         call this % ComputeJacobianMKL(dt,time,nEqn,nGradEqn,ComputeTimeDerivative)
      end if
      
      call this % SolveLUDirect(error)

      if (error .NE. 0) THEN
         WRITE(*,*) 'MKL Pardiso ERROR:', error
         error stop
      else
         this % converged = .TRUE.
      end if
    
#else
      error stop 'MKL is not linked properly'
#endif
      
   end subroutine solve
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine ComputeJacobianMKL(this,dt,time,nEqn,nGradEqn,ComputeTimeDerivative)
      use DenseBlockDiagonalMatrixClass
      implicit none
      !-----------------------------------------------------------
      class(MKLPardisoSolver_t), intent(inout) :: this
      real(kind=RP), intent(in)                :: dt
      real(kind=RP), intent(in)                :: time
      integer,       intent(in)                :: nEqn
      integer,       intent(in)                :: nGradEqn
      procedure(ComputeTimeDerivative_f)       :: ComputeTimeDerivative
      !-----------------------------------------------------------
      type(csrMat_t) :: B, Cmat !debug
      
      if (this % AIsPetsc) then
         call this % Jacobian % Compute (this % p_sem, nEqn, time, this % PETScA, ComputeTimeDerivative, ComputeTimeDerivative)
         
         call this % PETScA % GetCSRMatrix(this % A)
         call this % SetOperatorDt(dt)
         this % AIsPetsc = .FALSE.
         call this % PETScA % destruct
      else
         call this % Jacobian % Compute (this % p_sem, nEqn, time, this % A, ComputeTimeDerivative, ComputeTimeDerivative)
         
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
         
         
         call this % SetOperatorDt(dt)
      end if
      
      call this % FactorizeJacobian
      
   end subroutine ComputeJacobianMKL
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MKL_GetXValue(this,irow,x_i)       
      implicit none
      !-----------------------------------------------------------
      class(MKLPardisoSolver_t), intent(inout) :: this
      INTEGER                  , intent(in)    :: irow
      real(kind=RP)            , INTENT(OUT)   :: x_i
      !-----------------------------------------------------------
      
      x_i = this % x(irow)
      
   end subroutine MKL_GetXValue
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function MKL_GetX(this) result(x)
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(MKLPardisoSolver_t), INTENT(INOUT) :: this
      REAL(KIND=RP)                            :: x(this % DimPrb)
      !-----------------------------------------------------------
      
      x = this % x
      
   end function MKL_GetX
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MKL_destroy(this)       
      implicit none
      !-----------------------------------------------------------
      class(MKLPardisoSolver_t), intent(inout) :: this
      !-----------------------------------------------------------
      
      call this % A % destruct
      
      DEALLOCATE(this % b)
      DEALLOCATE(this % x)
      DEALLOCATE(this % Pardiso_pt)
      DEALLOCATE(this % Pardiso_iparm)
      DEALLOCATE(this % perm)
      
      if (this % Variable_dt) then
         call this % ALU % destruct
         deallocate (this % ALU)
      else
         nullify (this % ALU)
      end if
      
      this % AIsPrealloc = .FALSE.
      
    end subroutine MKL_destroy
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SetOperatorDt(this,dt)       
      implicit none
      !-----------------------------------------------------------
      class(MKLPardisoSolver_t), intent(inout) :: this
      real(kind=RP)            , intent(in)    :: dt
      !-----------------------------------------------------------
      
      this % Ashift = MatrixShift(dt)
      if (this % AIsPetsc) THEN
         call this % PETScA % shift(this % Ashift)
      else
         call this % A % Shift(this % Ashift)
      end if
      
    end subroutine SetOperatorDt
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine ReSetOperatorDt(this,dt)       
      implicit none
      !-----------------------------------------------------------
      class(MKLPardisoSolver_t), intent(inout) :: this
      real(kind=RP)            , intent(in)    :: dt
      !-----------------------------------------------------------
      real(kind=RP)                            :: shift
      !-----------------------------------------------------------
      
      shift = MatrixShift(dt)
      if ( AlmostEqual(shift,this % Ashift) ) return
      
      if (this % AIsPetsc) THEN
         call this % PETScA % shift(shift)
      else
         call this % A % Shift(-this % Ashift)
         call this % A % Shift(shift)
      end if
      this % Ashift = shift
      
      call this % FactorizeJacobian
      
    end subroutine ReSetOperatorDt
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function MKL_GetXnorm(this,TypeOfNorm) result(xnorm)
      implicit none
      !-----------------------------------------------------------
      class(MKLPardisoSolver_t), intent(inout) :: this
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
   end function MKL_GetXnorm
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
      class(MKLPardisoSolver_t), intent(inout) :: this
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

   subroutine MKL_ComputeAndFactorizeJacobian(self,nEqn, nGradEqn, F_J, dt, eps, mode_in)
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
      class(MKLPardisoSolver_t), intent(inout) :: self
      integer,                   intent(in)    :: nEqn, nGradEqn
      procedure(ComputeTimeDerivative_f)       :: F_J
      real(kind=RP), intent(in)                :: dt
      real(kind=RP), intent(in)                :: eps
      integer,       intent(in)                :: mode_in


!
!     Compute numerical Jacobian in the PETSc matrix
!     ----------------------------------------------
      if ( self % AIsPetsc) then
         call self % Jacobian % Compute (self % p_sem, nEqn, 0._RP, self % PETScA, F_J, eps_in = eps)
!
!        Transform the Jacobian to CSRMatrix
!        -----------------------------------
         call self % PETScA % GetCSRMatrix(self % A)
         call self % SetOperatorDt(dt)
!
!        Correct the shifted Jacobian values
!        -----------------------------------
         self % A % values = -dt * self % A % values
      
      else
         call self % Jacobian % Compute (self % p_sem, nEqn, 0._RP, self % A, F_J, F_J, eps, .false., mode_in)
         call self % SetOperatorDt(dt)
         
         self % A % values = -dt * self % A % values
      end if
      
      call self % FactorizeJacobian
      
   end subroutine MKL_ComputeAndFactorizeJacobian

   subroutine MKL_ReFactorizeJacobian(self) 
! 
!     ************************************************************************************* 
!        This subroutine changes the gamma0 coefficient from 1.0 to 1.5 by 
!        adding 0.5 to the diagonal. 
!     ************************************************************************************* 
! 
      implicit none 
      class(MKLPardisoSolver_t), intent(inout) :: self 
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
#ifdef HAS_MKL 
      call pardiso(self % Pardiso_pt, 1, 1, self % mtype, 12, self % ALU % num_of_Rows, self % ALU % values, &
                   self % ALU % rows, self % ALU % cols, self % perm, 1, self % Pardiso_iparm, 0, &
                   self % b, self % x, error)
#else 
      error stop 'MKL not linked correctly' 
#endif 
 
   end subroutine MKL_ReFactorizeJacobian 

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MKL_FactorizeJacobian(self)
      implicit none
      !-arguments---------------------------------------------------
      class(MKLPardisoSolver_t), intent(inout)  :: self
      !-local-variables---------------------------------------------
      integer     :: error
      !-------------------------------------------------------------
      
      call Stopwatch % Start("Sparse LU-Factorization")
      
      if (self % Variable_dt) then
         call self % ALU % destruct
         call self % ALU % constructWithCSRArrays  (self % A % Rows, &
                                                 self % A % Cols, &
                                                 self % A % Values)
      end if
      
!
!     Perform the factorization
!     -------------------------
#ifdef HAS_MKL
      call pardiso(self % Pardiso_pt, 1, 1, self % mtype, 12, self % ALU % num_of_Rows, self % ALU % values, &
                   self % ALU % rows, self % ALU % cols, self % perm, 1, self % Pardiso_iparm, 0, &
                   self % b, self % x, error)
#else
      error stop 'MKL not linked correctly'
#endif
      
      call Stopwatch % Pause("Sparse LU-Factorization")
   end subroutine MKL_FactorizeJacobian
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MKL_SolveLUDirect(self,error_out)
      implicit none
      class(MKLPardisoSolver_t), intent(inout)  :: self
      integer, optional        , intent(out)    :: error_out
!
!     ---------------
!     Local variables
!     ---------------
!
      integer     :: error

#ifdef HAS_MKL
      call pardiso(self % Pardiso_pt, 1, 1, self % mtype, 33, self % ALU % num_of_Rows, self % ALU % values, &
                   self % ALU % rows, self % ALU % cols, self % perm, 1, self % Pardiso_iparm, 0, &
                   self % b, self % x, error)
      
      if ( present(error_out) ) then
         error_out = error
      end if
      
#else
      error stop 'MKL not linked correctly'
#endif
   end subroutine MKL_SolveLUDirect
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine MKL_SetJacobian(this,Matrix)
      implicit none
      !-arguments-----------------------------------------------------------
      class(MKLPardisoSolver_t), intent(inout)  :: this
      class(Matrix_t)          , intent(in)     :: Matrix
      !---------------------------------------------------------------------
      
      select type(Matrix)
      class is(csrMat_t)
         this % A = Matrix
      class default
         error stop 'MKL_SetJacobian :: Wrong matrix type'
      end select
      
   end subroutine MKL_SetJacobian
END MODULE MKLPardisoSolverClass