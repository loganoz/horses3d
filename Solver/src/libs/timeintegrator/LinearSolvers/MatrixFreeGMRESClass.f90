!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!      MatrixFreeGMRESClass.f90
!      Created: 2017-04-10 10:006:00 +0100 
!      By: Carlos Redondo
!          Andrés Rueda (adapted to HORSES3D)
!
!      Class for solving linear systems using a matrix free GMRES
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
module MatrixFreeGMRESClass
   use GenericLinSolverClass
   use SMConstants
   use DGSEMClass
   use FTValueDictionaryClass
   use BDFFunctions
   implicit none
   
   private
   public   :: MatFreeGMRES_t 
   
   type, extends(GenericLinSolver_t) :: MatFreeGMRES_t
      integer                              :: m = 60           ! Number of GMRES iterations before restart -- Default petsc value m=30 
      integer                              :: maxiter = 500
      real(kind=RP)                        :: tol = 1e-15_RP
      real(kind=RP)                        :: res = -1._RP
      integer                              :: ERROR_CODE = 0
      real(kind=RP)                        :: norm0 = -1._RP
      real(kind=RP), allocatable           :: RHS(:)
      real(kind=RP), allocatable           :: x(:)
      real(kind=RP), allocatable           :: x0(:)
      real(kind=RP)                        :: rnorm
      integer                              :: Preconditioner   ! Which preconditioner is being used (PC_NONE, PC_GMRES, PC_BlockJacobi -last one to be developed)
      
      type(DGSem), pointer                 :: p_sem            ! Associated sem
      real(kind=RP), allocatable           :: F_Ur(:)          ! Qdot at the beginning of solving procedure
      real(kind=RP), allocatable           :: Ur(:)            ! Q at the beginning of solving procedure
      
      real(kind=RP)                        :: timesolve        ! Time at the solution
      real(kind=RP)                        :: dtsolve          ! dt for the solution
      
      ! Krylov subspace variables:
      real(kind=RP), allocatable, private  :: H(:,:)
      real(kind=RP), allocatable, private  :: W(:)
      real(kind=RP), allocatable, private  :: V(:,:)           ! Orthogonal vectors of Krylov subspace (Arnoldi)
      real(kind=RP), allocatable, private  :: Z(:,:)
      real(kind=RP), allocatable, private  :: Y(:)
      real(kind=RP), allocatable, private  :: cc(:)
      real(kind=RP), allocatable, private  :: ss(:)
      real(kind=RP), allocatable, private  :: g(:)
      
      TYPE(MatFreeGMRES_t), pointer, private :: PCsolver ! Inner GMRES solver for preconditioning
      
      contains
         ! Overriden procedures:
         procedure                           :: Construct => ConstructSolver
         procedure                           :: Destroy   => DestructSolver
         procedure                           :: SetRHSValue
         procedure                           :: SetRHSValues
         procedure                           :: SetOperatorDt
         procedure                           :: ReSetOperatorDt
         procedure                           :: GetXValue
         procedure                           :: SetRHS
         procedure                           :: Getxnorm    !Get solution norm
         procedure                           :: Getrnorm    !Get residual norm
         procedure                           :: Solve     => SolveGMRES
         ! Own procedures
         procedure                           :: SetTol
         procedure                           :: SetMaxIter
         procedure                           :: SetMaxInnerIter
         procedure                           :: SetInitialGuess
         
         ! Internal procedures:
         procedure :: p_F         ! Get the time derivative for a specific global Q
         procedure :: MatFreeAx   ! Matrix action
         procedure :: PC_GMRES_Ax ! Matrix free multiplication used for GMRES recursive preconditioning
   end type MatFreeGMRES_t

   abstract interface
      subroutine matmultsub(v,x, ComputeTimeDerivative)
         use SMConstants
         use DGSEMClass, only: ComputeQDot_FCN
         real(kind = RP), intent(in)         :: v(:)
         real(kind = RP), intent(out)        :: x(:)
         procedure(ComputeQDot_FCN)          :: ComputeTimeDerivative
      end subroutine
   end interface

!
!  Module variables
!  ----------------
   integer, parameter :: MAX_OUTER_ITER = 5000 ! this is a "security" limit, should never be reached....
   ! Preconditioners
   integer, parameter :: PC_NONE        = 0
   integer, parameter :: PC_GMRES       = 1
   integer, parameter :: PC_BlockJacobi = 2
   
contains
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine ConstructSolver(this,DimPrb,controlVariables, sem)
         implicit none
         !------------------------------------------------
         class(MatFreeGMRES_t)  , intent(inout), target :: this
         integer                , intent(in)            :: DimPrb
!~         integer, intent(in),optional                   :: m          ! to be deprecated
         TYPE(FTValueDictionary), intent(in), optional  :: controlVariables
         TYPE(DGSem), target                , optional  :: sem
         !------------------------------------------------
         character(len=LINE_LENGTH)                 :: pc
         !------------------------------------------------
         
!~         if(PRESENT(m)) this%m = m
         if (.not. present(sem)) ERROR stop ':: Matrix free GMRES needs sem'
         
         this % DimPrb = DimPrb
         allocate(this % RHS (DimPrb))
         allocate(this % x0  (DimPrb))
         allocate(this % F_Ur(DimPrb))
         allocate(this % Ur  (DimPrb))
         this%x0 = 0.0_RP
         allocate(this%V (DimPrb,this%m+1))
         allocate(this%H (this%m+1,this%m))
         allocate(this%W (DimPrb))
         allocate(this%Y (this%m))
         allocate(this%cc(this%m+1))
         allocate(this%ss(this%m+1))
         allocate(this%g (this%m+1))
         
         this % p_sem => sem
         
!        Set preconditioner
!        ------------------
         if ( present(controlVariables) ) then
            if ( controlVariables % containsKey("preconditioner") ) THEN
               pc = controlVariables % StringValueForKey("preconditioner",LINE_LENGTH)
               select case(pc)
!                 ********************
!                 GMRES preconditioner
!                 ********************
                  case('GMRES')
                     ! Allocate extra storage
                     allocate(this%Z(this%DimPrb,this%m+1))
                     ! Construct inner GMRES solver
                     allocate (this % PCsolver)
                     call this % PCsolver % Construct(dimprb,sem = sem)
                     call this % PCsolver % SetMaxInnerIter(15)      ! Hardcoded to 15
                     call this % PCsolver % SetMaxIter(30)           ! Hardcoded to 30... old: 15
                     ! Change this solver's definitions
                     this % maxiter = 60                        ! Hardcoded to 60... old: 30
                     this % Preconditioner = PC_GMRES
                     
!                 *****************
!                 No preconditioner
!                 *****************
                  case default
                     write(STD_OUT,*) pc, ' preconditioner not found. No preconditioner will be used for GMRES'
                     this % Preconditioner = PC_NONE
               end select
            else
               this % Preconditioner = PC_NONE
            end if
         end if
      end subroutine ConstructSolver
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE SetRHSValue(this, irow, value)
         IMPLICIT NONE
         !-----------------------------------------------------------
         CLASS(MatFreeGMRES_t), INTENT(INOUT) :: this
         INTEGER              , INTENT(IN)    :: irow
         REAL(KIND=RP)        , INTENT(IN)    :: value
         !-----------------------------------------------------------
         
         this % RHS (irow+1) = value
         
      END SUBROUTINE SetRHSValue
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE SetRHSValues(this, nvalues, irow, values)
         IMPLICIT NONE
         !------------------------------------------------------
         CLASS(MatFreeGMRES_t)       , INTENT(INOUT)     :: this
         INTEGER                     , INTENT(IN)        :: nvalues
         INTEGER      , DIMENSION(1:), INTENT(IN)        :: irow
         REAL(KIND=RP), DIMENSION(1:), INTENT(IN)        :: values
         !------------------------------------------------------
         INTEGER                                        :: i
         
         DO i=1, nvalues
            IF (irow(i)<0) CYCLE
            this % RHS(irow(i)+1) = values(i)
         END DO
         
      END SUBROUTINE SetRHSValues
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE GetXValue(this,irow,x_i)       
         IMPLICIT NONE
         !-----------------------------------------------------------
         CLASS(MatFreeGMRES_t), INTENT(INOUT) :: this
         INTEGER              , INTENT(IN)    :: irow
         REAL(KIND=RP)        , INTENT(OUT)   :: x_i
         !-----------------------------------------------------------
         
         x_i = this % x(irow+1)
         
      END SUBROUTINE GetXValue
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   FUNCTION Getxnorm(this,TypeOfNorm) RESULT(xnorm)
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(MatFreeGMRES_t), INTENT(INOUT) :: this
      CHARACTER(len=*)                         :: TypeOfNorm
      REAL(KIND=RP)                            :: xnorm
      !-----------------------------------------------------------
      
      SELECT CASE (TypeOfNorm)
         CASE ('infinity')
            xnorm = MAXVAL(ABS(this % x))
         CASE ('l2')
            xnorm = NORM2(this % x)
         CASE DEFAULT
            STOP 'MatFreeSmoothClass ERROR: Norm not implemented yet'
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
      CLASS(MatFreeGMRES_t), INTENT(INOUT) :: this
      REAL(KIND=RP)                            :: rnorm
      !-----------------------------------------------------------
      REAL(KIND=RP)                            :: residual(this % DimPrb)
      !-----------------------------------------------------------
      
      rnorm = this % rnorm
      
   END FUNCTION Getrnorm
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine SetMaxInnerIter(this,m)
         implicit none
         !------------------------------------------------
         class(MatFreeGMRES_t), intent(inout) :: this 
         integer                              :: m
         !------------------------------------------------
         
         this % m = m
         
      end subroutine SetMaxInnerIter
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine DestructSolver(this)
         implicit none
         class(MatFreeGMRES_t), intent(inout)          :: this
         
         deallocate(this%RHS)
         deallocate(this%x0)
         deallocate(this%V)
         deallocate(this%H)
         deallocate(this%W)
         deallocate(this%Y)
         deallocate(this%cc)
         deallocate(this%ss)
         deallocate(this%g)
         
         if (this % Preconditioner == PC_GMRES) then
            deallocate(this%Z)
            call this % PCsolver % destroy
         end if
      end subroutine DestructSolver
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE SetOperatorDt(this, dt)  ! Modify for matrix peconditioners...
         IMPLICIT NONE
         !------------------------------------------------
         CLASS(MatFreeGMRES_t), INTENT(INOUT) :: this
         REAL(KIND=RP)        , INTENT(IN)    :: dt
         !------------------------------------------------
         
         this % dtsolve = dt
         
      END SUBROUTINE SetOperatorDt
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ReSetOperatorDt(this, dt)  ! Modify for matrix peconditioners...
         IMPLICIT NONE
         !------------------------------------------------
         CLASS(MatFreeGMRES_t), INTENT(INOUT)  :: this
         REAL(KIND=RP)        , INTENT(IN)     :: dt
         !------------------------------------------------
         
         this % dtsolve = dt
         
      END SUBROUTINE ReSetOperatorDt
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine SetTol(this,tol)
         implicit none
         class(MatFreeGMRES_t), intent(inout)      :: this
         real(kind = RP)                        :: tol
         this%tol = MAX(tol, EPSILON(1._RP)) ! This limits the solver tolerance to the used real precicision
         
      end subroutine SetTol
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!    
      subroutine SetMaxIter(this,maxiter)
         implicit none
         class(MatFreeGMRES_t), intent(inout)          :: this
         integer                                :: maxiter
         
         this%maxiter = maxiter
      end subroutine SetMaxIter
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
! 
      subroutine SetRHS(this,RHS)
         implicit none
         class(MatFreeGMRES_t), intent(inout)   :: this
         real(kind = RP)      , intent(in)      :: RHS(:)
         
         this%RHS = RHS
      end subroutine  SetRHS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine SetInitialGuess(this,x0)
         implicit none
         class(MatFreeGMRES_t), intent(inout)         :: this
         real(kind = RP)                           :: x0(:)
         
         this%x0 = x0
      end subroutine
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      recursive subroutine innerGMRES(this, ComputeTimeDerivative)
         implicit none
         class(MatFreeGMRES_t), intent(inout)      :: this
         procedure(ComputeQDot_FCN)             :: ComputeTimeDerivative
         integer                                :: i,j,k, l, ii,kk, m
         real(kind = RP)                        :: tmp1, tmp2
         
         !Compute first krylov vector
         call this % MatFreeAx(this%x0,this%V(:,1), ComputeTimeDerivative)
         
         this%V(:,1) = this%RHS - this%V(:,1)   
         this%g(1) = NORM2(this%V(:,1))
         this%V(:,1) = this%V(:,1) / this%g(1)
         
         this%H = 0.0_RP
         m = this%m

         do j = 1,m ! Krylov loop
            select case (this % Preconditioner)
               case (PC_GMRES)
                  call this%PC_GMRES_Ax(this%V(:,j),this%Z(:,j), ComputeTimeDerivative)
                  call this % MatFreeAx(this%Z(:,j),this%W, ComputeTimeDerivative)
               case default ! PC_NONE
                  call this % MatFreeAx(this%V(:,j),this%W, ComputeTimeDerivative)
            end select
            
            do i = 1,j
               this%H(i,j) = dot_product(this%W,this%V(:,i))
               this%W = this%W - this%H(i,j) * this%V(:,i)
            end do
            this%H(j+1,j) = NORM2(this%W)
            if ((ABS(this%H(j+1,j)) .LT. this%tol)) then
               this%CONVERGED = .TRUE.
               this%res = this%tol
               this%niter = this%niter + 1                 
               m = j
               exit
            end if
            this%V(:,j+1) =  this%W / this%H(j+1,j) 
            do i = 1, j-1
               tmp1 = this%H(i,j)
               tmp2 = this%H(i+1,j)
               this%H(i,j) = this%cc(i) * tmp1 + this%ss(i) * tmp2
               this%H(i+1,j) = this%cc(i) * tmp2 - this%ss(i) * tmp1
            end do 
            tmp1 = SQRT(this%H(j,j)*this%H(j,j) + this%H(j+1,j)*this%H(j+1,j) )    
            if (ABS(tmp1) .LT. 1e-15_RP) then
               this%ERROR_CODE = -1
               RETURN
            end if
            this%cc(j) = this%H(j,j) / tmp1
            this%ss(j) = this%H(j+1,j) / tmp1
            this%g(j+1) = -this%ss(j) * this%g(j)
            this%g(j) = this%cc(j) * this%g(j)
            this%H(j,j) = this%cc(j) * this%H(j,j) + this%ss(j) * this%H(j+1,j) 
            this%res = ABS(this%g(j+1))
            this%niter = this%niter + 1
            if (this%res .LT. this%tol) then
               this%CONVERGED = .TRUE.
               m = j
               exit
            end if
            if (this%niter .GE. this%maxiter) then
               m = j
               this%CONVERGED = .FALSE.
               exit
            end if
         end do ! End of Krylov loop
         
         if (m > 0) then
            this%y(m) = this%g(m) / this%H(m,m)
            do ii = 1, m-1
               kk = m - ii
               tmp1 = this%g(kk)
               l = kk+1
               do while (l .LE. m)
                  tmp1 = tmp1 - this%H(kk,l) * this%y(l)
                  this%y(kk) = tmp1 / this%H(kk,kk)
                  l = l+1
               end do 
            end do
         end if ! m > 0
         
         select case (this % Preconditioner)
            case (PC_GMRES)
               this%x = this%x0 + MATMUL(this%Z(:,1:m),this%y(1:m))
            case default !PC_NONE
               this%x = this%x0 + MATMUL(this%V(:,1:m),this%y(1:m))
         end select
       end subroutine innerGMRES
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
! 
      recursive subroutine SolveGMRES(this, ComputeTimeDerivative, tol, maxiter,time,dt,computeA)
         implicit none
         !----------------------------------------------------
         class(MatFreeGMRES_t), intent(inout)      :: this
         procedure(ComputeQDot_FCN)                :: ComputeTimeDerivative
         real(kind=RP), optional                   :: tol
         integer      , optional                   :: maxiter
         real(kind=RP), optional                   :: time
         real(kind=RP), optional                   :: dt
         logical      , optional  , intent(inout)  :: computeA                !<> In case of block preconditioning, this tells the solver if the block preconditioner should be calculated
         !----------------------------------------------------
         integer                                :: i
         !----------------------------------------------------
         
!        Set the optional values
!        -----------------------
         
         if (present(tol)) this % tol = max(tol, epsilon(1._RP))
         if (present(maxiter)) this % maxiter = maxiter
         
         if ( present(time) ) then
            this % timesolve = time
         else
            ERROR stop ':: MatFreeGMRES needs the solution time'
         end if
         if ( present(dt) ) then
            this % dtsolve = dt
         else
            ERROR stop ':: MatFreeGMRES needs the dt'
         end if
         
         this % Ur   = this % p_sem % mesh % storage % Q
         this % F_Ur = this % p_F (this % Ur, ComputeTimeDerivative)    ! need to compute the time derivative?
         
         this%niter = 0
         this%CONVERGED = .FALSE.
         this%x0 = 0._RP
         
         if (this % Preconditioner == PC_GMRES) call this % PCsolver % SetTol(this % tol) ! Set the same tolerance for the first iteration of preconditioner... This is a first approximation and can be changed...
!
!        Outer (restart) loop
!        --------------------
         do i = 1, MAX_OUTER_ITER
            call innerGMRES(this, ComputeTimeDerivative)
            if (this%ERROR_CODE .NE. 0) then
               PRINT*, 'ERROR IN GMRES, ERROR CODE: ', this%ERROR_CODE
               call this%destroy()
               STOP
            end if
            if(this%CONVERGED .OR. (this%niter .GE. this%maxiter)) then
               call this % MatFreeAx(this % x, this % x0, ComputeTimeDerivative) ! using x0 as a temporary storing variable
               this % rnorm = norm2(this % RHS - this % x0)
               RETURN
            else
               this%x0 = this%x
            end if
         end do
         write(*,*) " ** WARNING!!: Reached max outer iter limit in GMRES ** "
          
!        Restore default values
!        ----------------------

         
      end subroutine SolveGMRES
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ---------------------------------------------------
!     Returns the preconditioning product   Pv = P⁻¹ * v 
!     ---------------------------------------------------
      subroutine PC_GMRES_Ax(this,v, Pv, ComputeTimeDerivative)
         implicit none
         !---------------------------------------------------------
         class(MatFreeGMRES_t), intent(inout) :: this
         real(kind=RP)        , intent(in)    :: v(:)
         real(kind=RP)        , intent(out)   :: Pv(:)
         procedure(ComputeQDot_FCN)           :: ComputeTimeDerivative
         !---------------------------------------------------------
          
         CALL this % PCSolver % SetRHS(v)
         CALL this % PCSolver % Solve(ComputeTimeDerivative,time = this % timesolve, dt = this % dtsolve)
         CALL this % PCSolver % Settol(1e-1_RP)                 ! TODO: aquí habrá que poner algo más elaborado
         Pv = this % PCSolver % x
!         n_preco_iter = n_preco_iter + PCSolver%niter     ! Not really needed
                        
      END subroutine PC_GMRES_Ax
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ---------------------------------------------------
!     Returns the (matrix-free) matrix vector product 
!     ---------------------------------------------------
      subroutine MatFreeAx(this,x, Ax, ComputeTimeDerivative)
         implicit none
         !---------------------------------------------------------
         class(MatFreeGMRES_t), intent(inout) :: this
         real(kind=RP), intent(in)            :: x (this % DimPrb)
         real(kind=RP), intent(out)           :: Ax(this % DimPrb)
         procedure(ComputeQDot_FCN)           :: ComputeTimeDerivative
         !---------------------------------------------------------
         real(kind=RP) :: eps, shift
         !---------------------------------------------------------
          
         eps = 1e-8_RP * (1._RP + norm2(x) )
         shift = BDF_MatrixShift(this % dtsolve)
         
         Ax = ( this % p_F(this % Ur + x * eps, ComputeTimeDerivative) - this % F_Ur)/ eps  + shift * x     !First order 
         !Ax = ( this % p_F(this % Ur + x * eps, ComputeTimeDerivative) - this % p_F(Ur - x * eps, ComputeTimeDerivative))  /(2._RP * eps)  - x / Comp_Dt   !Second order
      
      end subroutine MatFreeAx
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ---------------------------------------------------
!     Function to return the time derivative using a specific solution vector
!     ---------------------------------------------------
      FUNCTION p_F(this,u, computeTimeDerivative) RESULT(F)
         IMPLICIT NONE
         !---------------------------------------------------------
         CLASS(MatFreeGMRES_t), INTENT(INOUT) :: this
         REAL(KIND=RP)        , INTENT(IN)    :: u(this % DimPrb)
         procedure(ComputeQDot_FCN)           :: ComputeTimeDerivative
         REAL(KIND=RP)                        :: F(this % DimPrb)
         !---------------------------------------------------------
         REAL(KIND=RP)                        :: u_p(this % DimPrb)
         !---------------------------------------------------------
         
         ! Save original Q
         u_p = this % p_sem % mesh % storage % Q
         
         ! Obtain derivative with new Q
         this % p_sem % mesh % storage % Q = u
         CALL ComputeTimeDerivative(this % p_sem % mesh, this % timesolve + this % dtsolve, this % p_sem % externalState, this % p_sem % externalGradients)
         F = this % p_sem % mesh % storage % Qdot

         ! Restore original Q
         this % p_sem % mesh % storage % Q = u_p   ! TODO: this step can be avoided if Ur is not read in the "child" GMRES (the preconditioner)
      END FUNCTION p_F
end module MatrixFreeGMRESClass
