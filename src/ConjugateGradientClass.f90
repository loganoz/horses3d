!//////////////////////////////////////////////////////
!
!      Class for solving linear systems using CG
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
module ConjugateGradientClass

   use GenericLinSolverClass  , only: GenericLinSolver_t, MatrixShift_FCN, MatrixShift
   use SMConstants            , only: RP, STD_OUT, LINE_LENGTH
   use DGSEMClass             , only: DGSem, ComputeTimeDerivative_f
   use FTValueDictionaryClass , only: FTValueDictionary
   use MPI_Utilities          , only: infNorm, L2Norm, MPI_SumAll
   use PhysicsStorage        , only:  CTD_IGNORE_MODE

   USE CSRMatrixClass

   implicit none

   private
   public   :: ConjugateGradient_t

   abstract interface
      function MatrixAction_t(x)
         use SMConstants
         real(kind=RP), intent(in) :: x(:)
         real(kind=RP)             :: MatrixAction_t(size(x))
      end function
   end interface

!  ***********************
!  Conjugate Gradient class
!  ***********************

   type, extends(GenericLinSolver_t) :: ConjugateGradient_t

      !integer                              :: maxiter = 400
      !real(kind=RP)                        :: res = -1._RP
      !real(kind=RP)                        :: tol = 1e-15_RP

      TYPE(csrMat_t)                       :: A                                  ! Jacobian matrix
      real(kind=RP), allocatable           :: RHS(:)
      real(kind=RP), allocatable           :: x(:)
      real(kind=RP), allocatable           :: res(:)


      real(kind=RP)                        :: Ashift           
      real(kind=RP)                        :: dtsolve          ! dt for the solution
      real(kind=RP)                        :: timesolve        ! Time at the solution


!     Variables for matrix-free matrix action
!---------------------------------------
      real(kind=RP), allocatable           :: F_Ur(:)          ! Qdot at the beginning of solving procedure
      real(kind=RP), allocatable           :: Ur(:)            ! Q at the beginning of solving procedure



   contains
      ! Overridden procedures:
      procedure                           :: Construct         => CG_Construct
      procedure                           :: Destroy           => CG_Destruct
      procedure                           :: Solve             => CG_Solve
      procedure                           :: SetRHS            
      procedure                           :: SetRHSValues
      procedure                           :: SetRHSValue
      procedure                           :: GetX
      procedure                           :: GetXValue
      procedure                           :: Getxnorm
      procedure                           :: Getrnorm
      procedure                           :: SetOperatorDt
      procedure                           :: ReSetOperatorDt

      ! Own procedures

      ! Matrixfree 
      procedure                           :: MatrixAction  => CG_MatrixAction  
      procedure                           :: p_F   ! Get the time derivative for a specific global Q

      ! Preconditioners

   end type ConjugateGradient_t

contains

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!

   subroutine CG_Construct(this, DimPrb, globalDimPrb, nEqn, controlVariables, sem, MatrixShiftFunc)
      implicit none
      !------------------------------------------------
      class(ConjugateGradient_t)     , intent(inout), target :: this
      integer                , intent(in)            :: DimPrb
      integer                , intent(in)            :: globalDimPrb
      integer                , intent(in)            :: nEqn
      TYPE(FTValueDictionary), intent(in), optional  :: controlVariables
      TYPE(DGSem), target                , optional  :: sem
      procedure(MatrixShift_FCN)                     :: MatrixShiftFunc


      call this % GenericLinSolver_t % construct(DimPrb, globalDimPrb, nEqn,controlVariables,sem,MatrixShiftFunc)

      IF (.NOT. PRESENT(sem)) error stop 'Fatal error: IterativeSolver needs sem.'

      MatrixShift => MatrixShiftFunc

      this % p_sem => sem
      this % DimPrb = DimPrb

      WRITE(*,*), " this % DimPrb",  this % DimPrb
      allocate(this % x   (DimPrb))
      allocate(this % RHS (DimPrb))
      allocate(this % res (DimPrb))
      allocate(this % F_Ur(DimPrb))
      allocate(this % Ur  (DimPrb))

      CALL this % A % construct (num_of_Rows = DimPrb)


   end subroutine CG_Construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine CG_Destruct(this)
      implicit none
      !-arguments--------------------------------------
      class(ConjugateGradient_t), intent(inout)          :: this
      !------------------------------------------------

      deallocate(this % x)
      deallocate(this % RHS )
      deallocate(this % res )
      deallocate(this % F_Ur)
      deallocate(this % Ur  )

   end subroutine CG_Destruct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!

   subroutine CG_solve(this,nEqn, nGradEqn, ComputeTimeDerivative,tol,maxiter,time,dt,computeA)
      implicit none
      !----------------------------------------------------
      class(ConjugateGradient_t), target, intent(inout)      :: this
      integer,       intent(in)                 :: nEqn
      integer,       intent(in)                 :: nGradEqn
      procedure(ComputeTimeDerivative_f)        :: ComputeTimeDerivative
      real(kind=RP), optional                   :: tol
      integer      , optional                   :: maxiter
      real(kind=RP), optional                   :: time
      real(kind=RP), optional                   :: dt
      logical      , optional  , intent(inout)  :: computeA

      ! CG ddedicated vars 
      real(kind=RP), allocatable           :: P(:)
      real(kind=RP)                        :: alpha, beta  
      ! Local vars
      INTEGER                                 :: i

      ! if matrix based
      if ( present(ComputeA)) then
         if (ComputeA) then
            call this % Jacobian % Compute (this % p_sem, nEqn, time, this % A, ComputeTimeDerivative)
            call this % SetOperatorDt(dt)
            ComputeA = .FALSE.
         end if
      else 
         call this % Jacobian % Compute (this % p_sem, nEqn, time, this % A, ComputeTimeDerivative)
         call this % SetOperatorDt(dt)
      end if

      ! if matrix free
      if ( present(time) ) then
         this % timesolve = time
      else
         error stop ':: MatFreeCG needs the solution time'
      end if
      if ( present(dt) ) then
         this % dtsolve = dt
      else
         error stop ':: MatFreeCG needs the dt'
      end if

      if (.not. associated(this % p_sem) ) error stop 'MatFreeGMRES needs sem or MatrixAction'
      call this % p_sem % mesh % storage % local2GlobalQ(this % p_sem % NDOF)
      this % Ur   = this % p_sem % mesh % storage % Q
      this % F_Ur = this % p_F (this % Ur, ComputeTimeDerivative)    !

      !assert that A is symmtric before continuing
      WRITE(*,*) "A si symmetric ?",  this % A % IsSymmetric() 

      this % x = 0.0 
      this % niter = 0 

      !Initialize x Is it better than 0?
      ! DO i=1, this % DimPrb
      !    this % x(i) = this % RHS(i) / this % A % Values(this%A%Diag(i))
      ! END DO

   !!!!!!!!!!! matrix based 
   ! !    this % res =  this % RHS - CSR_MatVecMul(this % A ,this % x)
 
   ! !    allocate( P ( this % DimPrb))
   ! !    P = this % res

   ! !   ! Write(*,*)  "this % x", L2Norm( this % x), "this % RHS", L2Norm( this % RHS), "this % res", L2Norm( this % res)
      
   ! !    do while (L2Norm(this % res) .GE. tol .AND. this%niter .LE. maxiter) 
   ! !       alpha =  sum(this % res * this % res)/sum(P * (CSR_MatVecMul(this % A ,P)))
   ! !       beta = sum( this % res* this % res) 
   ! !       this % x =  this % x + alpha*P
   ! !       this % res = this % res  - alpha*CSR_MatVecMul(this % A ,P)

   ! !       beta = sum( this % res* this % res) /beta
   ! !       P = this % res + beta * P
   ! !       this % niter = this % niter + 1 
   ! !    end do

   ! !    !Write(*,*)  "this % x", L2Norm( this % x), "this % RHS", L2Norm( this % RHS), "this % res", L2Norm( this % res), "this % A", this % A % FrobeniusNorm()

     
   ! !    !call this % A % Visualize("mat.txt")

   ! !    deallocate(P)

      
          this % res =  this % RHS - this % MatrixAction(this % x, ComputeTimeDerivative)
 
      allocate( P ( this % DimPrb))
      P = this % res

     ! Write(*,*)  "this % x", L2Norm( this % x), "this % RHS", L2Norm( this % RHS), "this % res", L2Norm( this % res)
      
      do while (L2Norm(this % res) .GE. tol .AND. this%niter .LE. maxiter) 
         alpha =  sum(this % res * this % res)/sum(P * (this % MatrixAction(P,ComputeTimeDerivative)))
         beta = sum( this % res* this % res) 
         this % x =  this % x + alpha*P
         this % res = this % res  - alpha*this % MatrixAction(P,ComputeTimeDerivative)

         beta = sum( this % res* this % res) /beta
         P = this % res + beta * P
         this % niter = this % niter + 1 
      end do

      !Write(*,*)  "this % x", L2Norm( this % x), "this % RHS", L2Norm( this % RHS), "this % res", L2Norm( this % res), "this % A", this % A % FrobeniusNorm()

     
      !call this % A % Visualize("mat.txt")

      deallocate(P)

   end subroutine CG_solve

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   function GetX(this) result(x)
      implicit none
      !-----------------------------------------------------------
      class(ConjugateGradient_t), intent(inout)    :: this
      real(KIND=RP)                        :: x(this % DimPrb)
      !-----------------------------------------------------------

      x = this % x

   end function GetX
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE GetXValue(this,irow,x_i)
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(ConjugateGradient_t), INTENT(INOUT) :: this
      INTEGER              , INTENT(IN)    :: irow
      REAL(KIND=RP)        , INTENT(OUT)   :: x_i
      !-----------------------------------------------------------

      x_i = this % x(irow)

   END SUBROUTINE GetXValue
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SetRHS(this, RHS)
      implicit none
      !-arguments-----------------------------------------------------------
      class(ConjugateGradient_t)       , intent(inout)  :: this
      real(kind=RP)            , intent(in) :: RHS(this % DimPrb)
      !---------------------------------------------------------------------

      this % RHS = RHS
   end subroutine SetRHS

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SetRHSValue(this, irow, value)
      implicit none
      !-----------------------------------------------------------
      CLASS(ConjugateGradient_t), INTENT(INOUT) :: this
      INTEGER                  , INTENT(IN)    :: irow
      REAL(KIND=RP)            , INTENT(IN)    :: value
      !-----------------------------------------------------------

      this % RHS (irow) = value

   end subroutine SetRHSValue
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SetRHSValues(this, nvalues, irow, values)
      IMPLICIT NONE
      CLASS(ConjugateGradient_t)   , INTENT(INOUT)     :: this
      INTEGER                     , INTENT(IN)        :: nvalues
      INTEGER      , DIMENSION(1:), INTENT(IN)        :: irow
      REAL(KIND=RP), DIMENSION(1:), INTENT(IN)        :: values
      !------------------------------------------------------
      INTEGER                                        :: i

      do i=1, nvalues
         IF (irow(i)<0) CYCLE
         this % RHS(irow(i)) = values(i)
      end do

   end subroutine SetRHSValues
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   FUNCTION Getxnorm(this,TypeOfNorm) RESULT(xnorm)
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(ConjugateGradient_t), INTENT(INOUT) :: this
      CHARACTER(len=*)                         :: TypeOfNorm
      REAL(KIND=RP)                            :: xnorm
      !-----------------------------------------------------------

      SELECT CASE (TypeOfNorm)
       CASE ('infinity')
         xnorm = infNorm(this % x)
       CASE ('l2')
         xnorm = L2Norm(this % x)
       CASE DEFAULT
         error stop 'ConjugateGradientClass ERROR: Norm not implemented yet'
      END SELECT
   END FUNCTION Getxnorm
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   FUNCTION Getrnorm(this) RESULT(rnorm)
      IMPLICIT NONE
      !
      !     ----------------------------------------
      !     Currently implemented with infinity norm
      !     ----------------------------------------
      !
      !-----------------------------------------------------------
      CLASS(ConjugateGradient_t), INTENT(INOUT) :: this
      REAL(KIND=RP)                            :: rnorm
      !-----------------------------------------------------------
      REAL(KIND=RP)                            :: residual(this % DimPrb)
      !-----------------------------------------------------------

      rnorm = L2Norm( this % res)

   END FUNCTION Getrnorm
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE SetOperatorDt(this,dt)       
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(ConjugateGradient_t), INTENT(INOUT) :: this
      REAL(KIND=RP)           , INTENT(IN)    :: dt
      !-----------------------------------------------------------
      
      this % dtsolve = dt

      this % Ashift = MatrixShift(dt)
      CALL this % A % Shift(this % Ashift)
      
    END SUBROUTINE SetOperatorDt
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
    SUBROUTINE ReSetOperatorDt(this,dt)       
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(ConjugateGradient_t), INTENT(INOUT) :: this
      REAL(KIND=RP)           , INTENT(IN)    :: dt
      !-----------------------------------------------------------
      REAL(KIND=RP)                            :: shift
      !-----------------------------------------------------------
      
      shift = MatrixShift(dt)
      CALL this % A % Shift(-this % Ashift)
      CALL this % A % Shift(shift)
      
      this % Ashift = shift
      
    END SUBROUTINE ReSetOperatorDt
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ---------------------------------------------------
!     Returns the (matrix-free) matrix vector product 
!     ---------------------------------------------------
    function CG_MatrixAction(this,x, ComputeTimeDerivative) result(Ax)
      implicit none
      !-arguments-----------------------------------------------
      class(ConjugateGradient_t), intent(inout) :: this
      real(kind=RP), intent(in)            :: x (this % DimPrb)
      procedure(ComputeTimeDerivative_f)   :: ComputeTimeDerivative
      real(kind=RP)                        :: Ax(this % DimPrb)

      !-local-variables-----------------------------------------
      real(kind=RP) :: eps, shift
      !---------------------------------------------------------
      

      eps = 1e-8_RP * (1._RP + L2Norm(x) )
      shift = MatrixShift(this % dtsolve)
         
      Ax = ( this % p_F(this % Ur + x * eps, ComputeTimeDerivative) - this % F_Ur)/ eps  + shift * x     !First order 
      !Ax = ( this % p_F(this % Ur + x * eps, ComputeTimeDerivative) - this % p_F(Ur - x * eps, ComputeTimeDerivative))  /(2._RP * eps)  - x / Comp_Dt   !Second order

      
   end function CG_MatrixAction

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ---------------------------------------------------
!     Function to return the time derivative using a specific solution vector
!     ---------------------------------------------------
   function p_F(this,u, computeTimeDerivative) result(F)
      implicit none
      !-arguments-----------------------------------------------
      class(ConjugateGradient_t), INTENT(INOUT) :: this
      REAL(KIND=RP)        , INTENT(IN)    :: u(this % DimPrb)
      procedure(ComputeTimeDerivative_f)           :: ComputeTimeDerivative
      REAL(KIND=RP)                        :: F(this % DimPrb)
      !-local-variables-----------------------------------------
      REAL(KIND=RP)                        :: u_p(this % DimPrb)
      !---------------------------------------------------------
      
       ! Save original Q
       u_p = this % p_sem % mesh % storage % Q
      
       ! Obtain derivative with new Q
       this % p_sem % mesh % storage % Q = u
       call this % p_sem % mesh % storage % global2LocalQ
       CALL ComputeTimeDerivative(this % p_sem % mesh, this % p_sem % particles, this % timesolve + this % dtsolve, CTD_IGNORE_MODE)
       call this % p_sem % mesh % storage % local2GlobalQdot(this % p_sem % NDOF)
      
       F = this % p_sem % mesh % storage % Qdot

       ! Restore original Q
       this % p_sem % mesh % storage % Q = u_p   ! TODO: this step can be avoided if Ur is not read in the "child" GMRES (the preconditioner)
       call this % p_sem % mesh % storage % global2LocalQ

   end function p_F
   
end module ConjugateGradientClass
