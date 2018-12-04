!
!//////////////////////////////////////////////////////
!
!   @File:    MKLPardisoSolverClass.f90
!   @Author:  Andrés Rueda (am.rueda@upm.es)
!   @Created: 2017-04-10 10:006:00 +0100
!   @Last revision date: Tue Dec  4 21:53:46 2018
!   @Last revision author: Andrés Rueda (am.rueda@upm.es)
!   @Last revision commit: 9b3844379fde2350e64816efcdf3bf724c8b3041
!
!//////////////////////////////////////////////////////
!
!     Class for solving linear systems using MKL version of Pardiso
!     -> It is possible to construct the matrix using PETSc. Thhis option is currently deactivated... deprecate??
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "Includes.h"
#ifdef HAS_PETSC
#include "petsc/finclude/petsc.h"
#endif
MODULE MKLPardisoSolverClass
   USE GenericLinSolverClass
   USE CSRMatrixClass
   use PETScMatrixClass
   USE SMConstants
   use DGSEMClass
   use TimeIntegratorDefinitions
#ifdef HAS_PETSC
   use petsc
#endif
   implicit none
#ifdef HAS_PETSC
#include <petsc.h>
#endif
   TYPE, EXTENDS(GenericLinSolver_t) :: MKLPardisoSolver_t
      TYPE(csrMat_t)                             :: A                                  ! Jacobian matrix
      type(PETSCMatrix_t)                        :: PETScA
      real(kind=RP), DIMENSION(:), ALLOCATABLE   :: x                                  ! Solution vector
      real(kind=RP), DIMENSION(:), ALLOCATABLE   :: b                                  ! Right hand side
      real(kind=RP)                              :: Ashift
      LOGICAL                                    :: AIsPrealloc
      
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
      PROCEDURE :: solve
      procedure :: SolveLUDirect                => MKL_SolveLUDirect
      procedure :: SetRHSValue                  => MKL_SetRHSValue
      procedure :: SetRHS                       => MKL_SetRHS
      PROCEDURE :: GetXValue
      PROCEDURE :: destroy
      PROCEDURE :: SetOperatorDt
      PROCEDURE :: ReSetOperatorDt
      procedure :: ComputeJacobianMKL
      !Functions:
      PROCEDURE :: Getxnorm    !Get solution norm
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
   
   subroutine ConstructMKLContext(this,DimPrb,controlVariables,sem,MatrixShiftFunc)
      implicit none
      !-----------------------------------------------------------
      class(MKLPardisoSolver_t), intent(inout), TARGET :: this
      integer                  , intent(in)            :: DimPrb
      TYPE(FTValueDictionary)  , intent(in), OPTIONAL  :: controlVariables
      TYPE(DGSem), TARGET                  , OPTIONAL  :: sem
      procedure(MatrixShift_FCN)                       :: MatrixShiftFunc
      !-----------------------------------------------------------
#ifdef HAS_PETSC
      PetscErrorCode :: ierr
#endif
      !-----------------------------------------------------------
      
      if ( controlVariables % containsKey("jacobian flag") ) then
         this % JacobianComputation = controlVariables % integerValueForKey("jacobian flag")
      end if
      
      MatrixShift => MatrixShiftFunc
      
      this % p_sem => sem
      
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
         ERROR stop "MKL-Pardiso needs PETSc for constructung the Jacobian Matrix"
#endif
         call this % PETScA % construct (num_of_Rows = DimPrb, withMPI = .FALSE.)
      else
         call this % A % construct(num_of_Rows = DimPrb, withMPI = .false.)
         
      end if

#ifdef HAS_MKL
      call pardisoinit(this % Pardiso_pt, this % mtype, this % Pardiso_iparm)
#else
      stop 'MKL not linked correctly'
#endif
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
      class(MKLPardisoSolver_t), intent(inout) :: this
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
      INTEGER                                  :: error
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
      
!~    	call mkl_set_num_threads( 4 )
      
      !-----------------------
      ! Solve the system using MKL - Pardiso!!
!~    phase = 33            ! Solve, iterative refinement
      call pardiso(  pt      = this % Pardiso_pt    ,     &
                     maxfct  = 1                    ,     &     ! Set up space for 1 matrix at most
                     mnum    = 1                    ,     &     ! Matrix to use in the solution phase (1st and only one)
                     mtype   = this % mtype         ,     &
                     phase   = 13                   ,     &     !  
                     n       = this % DimPrb        ,     &     ! Number of equations
                     values  = this % A % Values    ,     & 
                     rows    = this % A % Rows      ,     &
                     cols    = this % A % Cols      ,     &
                     perm    = this % perm          ,     &     ! ...
                     nrhs    = 1                    ,     &     ! Only one right hand side 
                     iparm   = this % Pardiso_iparm ,     &
                     msglvl  = 0                    ,     &     ! 1: verbose... Too much printing
                     b       = this % b             ,     &
                     x       = this % x             ,     &
                     ierror  = error              )

      if (error .NE. 0) THEN
         WRITE(*,*) 'MKL Pardiso ERROR:', error
         stop
      else
         this % converged = .TRUE.
      end if
    
#else
      stop 'MKL is not linked properly'
#endif
      
   end subroutine solve
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine ComputeJacobianMKL(this,dt,time,nEqn,nGradEqn,ComputeTimeDerivative)
      implicit none
      !-----------------------------------------------------------
      class(MKLPardisoSolver_t), intent(inout) :: this
      real(kind=RP), intent(in)                :: dt
      real(kind=RP), intent(in)                :: time
      integer,       intent(in)                :: nEqn
      integer,       intent(in)                :: nGradEqn
      procedure(ComputeTimeDerivative_f)       :: ComputeTimeDerivative
      !-----------------------------------------------------------
      
      if (this % AIsPetsc) then
         call this % ComputeJacobian(this % PETScA,dt,time,nEqn,nGradEqn,ComputeTimeDerivative)
         
         call this % PETScA % GetCSRMatrix(this % A)
         this % AIsPetsc = .FALSE.
      else
         call this % ComputeJacobian(this % A,dt,time,nEqn,nGradEqn,ComputeTimeDerivative)
      end if
      
   end subroutine ComputeJacobianMKL
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine GetXValue(this,irow,x_i)       
      implicit none
      !-----------------------------------------------------------
      class(MKLPardisoSolver_t), intent(inout) :: this
      INTEGER                  , intent(in)    :: irow
      real(kind=RP)            , INTENT(OUT)   :: x_i
      !-----------------------------------------------------------
      
      x_i = this % x(irow)
      
   end subroutine GetXValue
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function GetX(this) result(x)
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(MKLPardisoSolver_t), INTENT(INOUT) :: this
      REAL(KIND=RP)                            :: x(this % DimPrb)
      !-----------------------------------------------------------
      
      x = this % x
      
   end function GetX
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine destroy(this)       
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
      
      call this % PETScA % destruct
      this % AIsPetsc    = .TRUE.
      this % AIsPrealloc = .FALSE.
      
    end subroutine destroy
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
      if (this % AIsPetsc) THEN
         call this % PETScA % shift(shift)
      else
         call this % A % Shift(-this % Ashift)
         call this % A % Shift(shift)
      end if
      this % Ashift = shift
      
    end subroutine ReSetOperatorDt
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function Getxnorm(this,TypeOfNorm) result(xnorm)
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
   end function Getxnorm
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
      
      residual = this % b - CSR_MatVecMul(this % A, this % x)
      rnorm = MAXVAL(ABS(residual))
      
      
      !rnorm = NORM2(this % x)
      
   end function Getrnorm
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!

   subroutine MKL_ComputeAndFactorizeJacobian(self,nEqn, nGradEqn, F_J, dt, eps)
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
!
!     ---------------
!     Local variables
!     ---------------
!
      integer     :: error

!
!     Compute numerical Jacobian in the PETSc matrix
!     ----------------------------------------------
      if ( self % AIsPetsc) then
         call self % ComputeJacobian(self  % PETScA,dt,0._RP,nEqn,nGradEqn,F_J,eps)
!
!        Transform the Jacobian to CSRMatrix
!        -----------------------------------
         call self % PETScA % GetCSRMatrix(self % A)
!
!        Correct the shifted Jacobian values
!        -----------------------------------
         self % A % values = -dt * self % A % values
      
      else
         call self % ComputeJacobian(self % A,dt,0._RP,nEqn,nGradEqn,F_J,eps)
         
         self % A % values = -dt * self % A % values
      end if
!
!     Perform the factorization
!     -------------------------
#ifdef HAS_MKL
      call pardiso(self % Pardiso_pt, 1, 1, self % mtype, 12, self % A % num_of_Rows, self % A % values, &
                   self % A % rows, self % A % cols, self % perm, 1, self % Pardiso_iparm, 0, &
                   self % b, self % x, error)
#else
      stop 'MKL not linked correctly'
#endif

   end subroutine MKL_ComputeAndFactorizeJacobian

   subroutine MKL_SolveLUDirect(self)
      implicit none
      class(MKLPardisoSolver_t), intent(inout)  :: self
!
!     ---------------
!     Local variables
!     ---------------
!
      integer     :: error

#ifdef HAS_MKL
      call pardiso(self % Pardiso_pt, 1, 1, self % mtype, 33, self % A % num_of_Rows, self % A % values, &
                   self % A % rows, self % A % cols, self % perm, 1, self % Pardiso_iparm, 0, &
                   self % b, self % x, error)
#else
      stop 'MKL not linked correctly'
#endif
   end subroutine MKL_SolveLUDirect

END MODULE MKLPardisoSolverClass
