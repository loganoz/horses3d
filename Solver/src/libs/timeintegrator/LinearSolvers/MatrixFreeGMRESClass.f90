!//////////////////////////////////////////////////////
!
!      Class for solving linear systems using a matrix free GMRES
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
module MatrixFreeGMRESClass
   use GenericLinSolverClass  , only: GenericLinSolver_t, MatrixShift_FCN, MatrixShift
   use SMConstants            , only: RP, STD_OUT, LINE_LENGTH
   use DGSEMClass             , only: DGSem, ComputeTimeDerivative_f
   use FTValueDictionaryClass , only: FTValueDictionary
   use MatrixClass            , only: Matrix_t, DenseBlockDiagMatrix_t, SparseBlockDiagMatrix_t
   use PhysicsStorage        , only: NCONS, NGRAD, CTD_IGNORE_MODE
   use MPI_Utilities          , only: infNorm, L2Norm, MPI_SumAll
   implicit none
   
   private
   public   :: MatFreeGMRES_t
   
   abstract interface
      function MatrixAction_t(x)
         use SMConstants
         real(kind=RP), intent(in) :: x(:)
         real(kind=RP)             :: MatrixAction_t(size(x))
      end function
   end interface
   
!  ***********************
!  Matrix-Free GMRES class
!  ***********************
   type, extends(GenericLinSolver_t) :: MatFreeGMRES_t
      integer                              :: m = 500           ! Number of GMRES iterations before restart -- Default petsc value m=30 
      integer                              :: maxiter = 500
      real(kind=RP)                        :: tol = 1e-15_RP
      real(kind=RP)                        :: res = -1._RP
      integer                              :: ERROR_CODE = 0
      real(kind=RP)                        :: norm0 = -1._RP
      real(kind=RP), allocatable           :: RHS(:)
      real(kind=RP), allocatable           :: x(:)
      real(kind=RP), allocatable           :: x0(:)
      real(kind=RP)                        :: rnorm
      
      real(kind=RP)                        :: timesolve        ! Time at the solution
      real(kind=RP)                        :: dtsolve          ! dt for the solution
      
!     Variables for matrix-free matrix action
!     ---------------------------------------
      real(kind=RP), allocatable           :: F_Ur(:)          ! Qdot at the beginning of solving procedure
      real(kind=RP), allocatable           :: Ur(:)            ! Q at the beginning of solving procedure
      
!     Variables for user-defined matrix action
!     ----------------------------------------
      logical                                    :: UserDef_Ax = .FALSE.
      procedure(MatrixAction_t), nopass, pointer :: User_MatrixAction => null()
      
!     Variables for preconditioners
!     -----------------------------
      integer                                :: Preconditioner    ! Which preconditioner is being used (PC_NONE, PC_GMRES, PC_BlockJacobi)
      ! Block-Jacobi:
      class(Matrix_t), allocatable           :: BlockA            ! Block-diagonal Jacobian matrix for BlockJacobi preconditioner
      class(Matrix_t), allocatable           :: BlockPreco        ! LU factorized Block-diagonal Jacobian matrix for BlockJacobi preconditioner (in each block Matrix = L+U, Indexes = LU-pivots)
      ! GMRES:
      type(MatFreeGMRES_t), pointer, private :: PCsolver          ! Inner GMRES solver for preconditioning
      
!     Krylov subspace variables
!     -------------------------
      real(kind=RP), allocatable, private  :: H(:,:)
      real(kind=RP), allocatable, private  :: W(:)
      real(kind=RP), allocatable, private  :: V(:,:)           ! Orthogonal vectors of Krylov subspace (Arnoldi)
      real(kind=RP), allocatable, private  :: Z(:,:)
      real(kind=RP), allocatable, private  :: Y(:)
      real(kind=RP), allocatable, private  :: cc(:)
      real(kind=RP), allocatable, private  :: ss(:)
      real(kind=RP), allocatable, private  :: g(:)
      
      
      contains
         ! Overridden procedures:
         procedure                           :: Construct         => GMRES_Construct
         procedure                           :: Destroy           => GMRES_Destruct
         procedure                           :: SetRHSValue
         procedure                           :: SetRHSValues
         procedure                           :: SetOperatorDt
         procedure                           :: ReSetOperatorDt
         procedure                           :: GetXValue
         procedure                           :: GetX
         procedure                           :: SetRHS            => GMRES_SetRHS
         procedure                           :: Getxnorm
         procedure                           :: Getrnorm
         procedure                           :: Solve             => GMRES_Solve
         ! Own procedures
         procedure                           :: SetTol
         procedure                           :: SetMaxIter
         procedure                           :: SetMaxInnerIter
         procedure                           :: SetInitialGuess
         procedure                           :: SetMatrixAction   => GMRES_SetMatrixAction
         
         ! Internal procedures:
         procedure :: p_F                                   ! Get the time derivative for a specific global Q
         procedure :: MatrixAction  => GMRES_MatrixAction   ! Matrix action
         
         ! Preconditioner action procedures
         procedure :: PC_GMRES_Ax         ! P⁻¹x for GMRES recursive preconditioning
         procedure :: PC_BlockJacobi_Ax   ! P⁻¹x for Block-Jacobi preconditioning
   end type MatFreeGMRES_t

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
      recursive subroutine GMRES_Construct(this,DimPrb, globalDimPrb, nEqn,controlVariables, sem,MatrixShiftFunc)
         implicit none
         !------------------------------------------------
         class(MatFreeGMRES_t)  , intent(inout), target :: this
         integer                , intent(in)            :: DimPrb
         integer                , intent(in)            :: globalDimPrb
         integer                , intent(in)            :: nEqn
         TYPE(FTValueDictionary), intent(in), optional  :: controlVariables
         TYPE(DGSem), target                , optional  :: sem
         procedure(MatrixShift_FCN)                     :: MatrixShiftFunc
         !------------------------------------------------
         character(len=LINE_LENGTH)                 :: pc
         integer                                    :: nelem
         integer                                    :: k, Nx, Ny, Nz
         integer, allocatable                       :: ndofelm(:)
         !------------------------------------------------
         
         call this % GenericLinSolver_t % construct(DimPrb, globalDimPrb, nEqn,controlVariables,sem,MatrixShiftFunc)
         
         if ( present(sem) ) this % p_sem => sem
         
         MatrixShift => MatrixShiftFunc
         
         this % DimPrb = DimPrb
         
!        ***********************************         
!        Set variables from controlVariables
!        ***********************************
         if ( present(controlVariables) ) then
            
!           Inner iterations (before restart)
!           *********************************
            if ( controlVariables % containsKey("gmres inner iterations") ) then
               this % m = controlVariables % integerValueForKey("gmres inner iterations")
            end if
            
!           Preconditioner
!           **************
            pc = controlVariables % StringValueForKey("preconditioner",LINE_LENGTH)
            select case(pc)
!              
!              GMRES preconditioner
!              --------------------
               case('GMRES')
                  if (.not. associated(this % p_sem) ) error stop 'MatFreeGMRES needs sem for "GMRES" preconditioner'
               
                  ! Allocate extra storage
                  allocate(this%Z(this%DimPrb,this%m+1))
                  ! Construct inner GMRES solver
                  allocate (this % PCsolver)
                  call this % PCsolver % Construct(dimprb, globalDimPrb, nEqn,sem = sem, MatrixShiftFunc = MatrixShiftFunc)
                  call this % PCsolver % SetMaxInnerIter(15)      ! Hardcoded to 15
                  call this % PCsolver % SetMaxIter(30)           ! Hardcoded to 30... old: 15
                  ! Change this solver's definitions
                  this % maxiter = 60                             ! Hardcoded to 60... old: 30
                  this % Preconditioner = PC_GMRES
!              
!              Block-Jacobi preconditioner
!              ---------------------------
               case('Block-Jacobi')
                  if (.not. associated(this % p_sem) ) error stop 'MatFreeGMRES needs sem for "BlockJacobi" preconditioner'
                  
                  allocate(this%Z(this%DimPrb,this%m+1))
                  
                  if ( controlVariables % containsKey("matrix type") ) then
                     
                     select case ( trim( controlVariables % StringValueForKey("matrix type",LINE_LENGTH) ) )
                        case ('dense block-diagonal')
                           allocate(DenseBlockDiagMatrix_t :: this % BlockA)
                           allocate(DenseBlockDiagMatrix_t :: this % BlockPreco)
                        case ('sparse block-diagonal')
                           allocate(DenseBlockDiagMatrix_t :: this % BlockA)
                           allocate(DenseBlockDiagMatrix_t :: this % BlockPreco)
                        case default
                           write(STD_OUT,'(3A)') 'Not recognized matrix type "', trim( controlVariables % StringValueForKey("matrix type",LINE_LENGTH) ), '"'
                           write(STD_OUT,'(A)' ) '-> Defaulting to "dense block-diagonal"'
                           allocate(SparseBlockDiagMatrix_t :: this % BlockA)
                           allocate(SparseBlockDiagMatrix_t :: this % BlockPreco)
                     end select
                  else
                     allocate(DenseBlockDiagMatrix_t :: this % BlockA)
                     allocate(DenseBlockDiagMatrix_t :: this % BlockPreco)
                  end if
                  
                  nelem = sem % mesh % no_of_elements
                  call this % BlockA     % construct (num_of_Blocks = nelem)
                  call this % BlockPreco % construct (num_of_Blocks = nelem)
                  
                  allocate ( ndofelm(nelem) )
                  do k = 1, nelem
                     Nx = sem % mesh % elements(k) % Nxyz(1)
                     Ny = sem % mesh % elements(k) % Nxyz(2)
                     Nz = sem % mesh % elements(k) % Nxyz(3)
                     ndofelm(k) = NCONS*(Nx+1)*(Ny+1)*(Nz+1)
                  end do
                  call this % BlockPreco % PreAllocate(nnzs = ndofelm)
                  
                  deallocate ( ndofelm )
                  
                  this % Preconditioner = PC_BlockJacobi
!                 
!              No preconditioner
!              -----------------
               case default
                  write(STD_OUT,*) trim(pc), ' preconditioner not found. No preconditioner will be used for GMRES'
                  this % Preconditioner = PC_NONE
            end select
         else
            this % Preconditioner = PC_NONE
         end if
         
         if ( present(sem) .and. this % Preconditioner == PC_BlockJacobi ) then
            call this % Jacobian % Configure (sem % mesh, nEqn, this % BlockA)
         end if
         
!        ******************
!        Allocate variables
!        ******************
         
         allocate(this % RHS (DimPrb))
         allocate(this % x0  (DimPrb))
         allocate(this % x   (DimPrb))
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
         
      end subroutine GMRES_Construct
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
         
         this % RHS (irow) = value
         
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
            this % RHS(irow(i)) = values(i)
         END DO
         
      END SUBROUTINE SetRHSValues
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine GMRES_SetRHS(this, RHS)
      implicit none
      !-arguments-----------------------------------------------------------
      class(MatFreeGMRES_t), intent(inout)  :: this
      real(kind=RP)            , intent(in) :: RHS(this % DimPrb)
      !---------------------------------------------------------------------
      
      this % RHS = RHS
      
   end subroutine GMRES_SetRHS
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
         
         x_i = this % x(irow)
         
      END SUBROUTINE GetXValue
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      function GetX(this) result(x)
         IMPLICIT NONE
         !-----------------------------------------------------------
         CLASS(MatFreeGMRES_t), INTENT(INOUT) :: this
         REAL(KIND=RP)                        :: x(this % DimPrb)
         !-----------------------------------------------------------
         
         x = this % x
         
      end function GetX
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
            xnorm = infNorm(this % x)
         CASE ('l2')
            xnorm = L2Norm(this % x)
         CASE DEFAULT
            error stop 'MatFreeSmoothClass ERROR: Norm not implemented yet'
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
      
!~      rnorm = this % rnorm
      rnorm = this % res
      
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
      recursive subroutine GMRES_Destruct(this)
         implicit none
         !-arguments--------------------------------------
         class(MatFreeGMRES_t), intent(inout)          :: this
         !------------------------------------------------
         
         deallocate(this % RHS )
         deallocate(this % x0  )
         deallocate(this % F_Ur)
         deallocate(this % Ur  )
         
         deallocate(this%V)
         deallocate(this%H)
         deallocate(this%W)
         deallocate(this%Y)
         deallocate(this%cc)
         deallocate(this%ss)
         deallocate(this%g)
         
         select case (this % Preconditioner)
            case (PC_GMRES)
               deallocate(this%Z)
               call this % PCsolver % destroy
            case (PC_BlockJacobi)
               deallocate(this%Z)
               call this % BlockA % destruct
               call this % BlockPreco % destruct
         end select
         
         this % UserDef_Ax        = .FALSE.
         this % User_MatrixAction => null()
         
      end subroutine GMRES_Destruct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE SetOperatorDt(this, dt)  ! Modify for matrix peconditioners...
         IMPLICIT NONE
         !-arguments--------------------------------------
         CLASS(MatFreeGMRES_t), INTENT(INOUT) :: this
         REAL(KIND=RP)        , INTENT(IN)    :: dt
         !------------------------------------------------
         REAL(KIND=RP) :: shift

         this % dtsolve = dt
         
         if (this % Preconditioner == PC_BlockJacobi) then
            shift = MatrixShift(dt)
            call this % BlockA % shift( shift )
            call this % BlockA % FactorizeBlocks_LU(this % BlockPreco)
         end if
         
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
         REAL(KIND=RP) :: shift
         
         this % dtsolve = dt
         
         if (this % Preconditioner == PC_BlockJacobi) then
            shift = MatrixShift(dt)
            call this % BlockA % shift( shift )
            call this % BlockA % FactorizeBlocks_LU(this % BlockPreco)
         end if
         
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
         class(MatFreeGMRES_t), intent(inout)   :: this
         procedure(ComputeTimeDerivative_f)     :: ComputeTimeDerivative
         integer                                :: i,j,k, l, ii,kk, m
         real(kind = RP)                        :: tmp1, tmp2
         
         !Compute first krylov vector
         call this % MatrixAction(this%x0,this%V(:,1), ComputeTimeDerivative)
         
         this%V(:,1) = this%RHS - this%V(:,1)   
         this%g(1) = L2Norm(this%V(:,1))
         this%V(:,1) = this%V(:,1) / this%g(1)
         
         this%H = 0.0_RP
         m = this%m
         do j = 1,m ! Krylov loop
            select case (this % Preconditioner)
               case (PC_GMRES)
                  call this % PC_GMRES_Ax(this%V(:,j),this%Z(:,j), ComputeTimeDerivative)
                  call this % MatrixAction(this%Z(:,j),this%W, ComputeTimeDerivative)
               case (PC_BlockJacobi)
                  call this % PC_BlockJacobi_Ax(this%V(:,j),this%Z(:,j))
                  call this % MatrixAction(this%Z(:,j),this%W, ComputeTimeDerivative)
               case default ! PC_NONE
                  call this % MatrixAction(this%V(:,j),this%W, ComputeTimeDerivative)
            end select
            
            do i = 1,j
               this%H(i,j) = dot_product(this%W,this%V(:,i))
               call MPI_SumAll(this%H(i,j))
               
               this%W = this%W - this%H(i,j) * this%V(:,i)
            end do
            
            this%H(j+1,j) = L2Norm(this%W)
            
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
            case (PC_NONE)
               this%x = this%x0 + MATMUL(this%V(:,1:m),this%y(1:m))
            case default !PC_GMRES, PC_BlockJacobi
               this%x = this%x0 + MATMUL(this%Z(:,1:m),this%y(1:m))
         end select
       end subroutine innerGMRES
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
! 
      recursive subroutine GMRES_Solve(this, nEqn, nGradEqn, ComputeTimeDerivative, tol, maxiter,time,dt,computeA)
         implicit none
         !----------------------------------------------------
         class(MatFreeGMRES_t), target, intent(inout)      :: this
         integer,       intent(in)                :: nEqn
         integer,       intent(in)                :: nGradEqn
         procedure(ComputeTimeDerivative_f)                :: ComputeTimeDerivative
         real(kind=RP), optional                   :: tol
         integer      , optional                   :: maxiter
         real(kind=RP), optional                   :: time
         real(kind=RP), optional                   :: dt
         logical      , optional  , intent(inout)  :: computeA                !<> In case of block preconditioning, this tells the solver if the block preconditioner should be calculated
         !----------------------------------------------------
         integer                                :: i
         real(kind=RP) :: xnorm
         !----------------------------------------------------
         
!        Set the optional values
!        -----------------------
         
         if (present(tol)    ) this % tol = max(tol, epsilon(1._RP))
         if (present(maxiter)) this % maxiter = maxiter
         
         if ( present(time) ) then
            this % timesolve = time
         else
            error stop ':: MatFreeGMRES needs the solution time'
         end if
         if ( present(dt) ) then
            this % dtsolve = dt
         else
            error stop ':: MatFreeGMRES needs the dt'
         end if
         
         if (.not. this % UserDef_Ax) then
            if (.not. associated(this % p_sem) ) error stop 'MatFreeGMRES needs sem or MatrixAction'
            call this % p_sem % mesh % storage % local2GlobalQ(this % p_sem % NDOF)
            this % Ur   = this % p_sem % mesh % storage % Q
            this % F_Ur = this % p_F (this % Ur, ComputeTimeDerivative)    ! need to compute the time derivative?
         end if
         
!        Preconditioner initializations
!        ------------------------------
         select case (this % Preconditioner)
            case (PC_BlockJacobi)
               if ( present(ComputeA)) then
                  if (ComputeA) then
                     call this % Jacobian % Compute (this % p_sem, nEqn, time, this % BlockA, ComputeTimeDerivative, BlockDiagonalized = .TRUE.)
                     call this % SetOperatorDt(dt) ! the matrix is factorized inside
                     
                     ComputeA = .FALSE.
                  end if
               else
                  call this % Jacobian % Compute (this % p_sem, nEqn, time, this % BlockA, ComputeTimeDerivative, BlockDiagonalized = .TRUE.)
                  call this % SetOperatorDt(dt) ! the matrix is factorized inside
                  
               end if
         end select
         
!        Initializations
!        ---------------
         
         this%niter = 0
         this%CONVERGED = .FALSE.
         
         ! Set initial guess to Knoll: P⁻¹b
!~         select case (this % Preconditioner)
!~            case (PC_GMRES)
!~               call this % PC_GMRES_Ax(this%RHS,this%x0, ComputeTimeDerivative)
!~            case (PC_BlockJacobi)
!~               call this % PC_BlockJacobi_Ax(this%RHS,this%x0)
!~            case default ! PC_NONE
               this % x0 = 0._RP ! The actual "preconditioned" guess would be this%RHS, but 0 is a good initial guess for the Newton linearized system 
!~         end select
         
         if (this % Preconditioner == PC_GMRES) call this % PCsolver % SetTol(this % tol) ! Set the same tolerance for the first iteration of preconditioner... This is a first approximation and can be changed...
!
!        Outer (restart) loop
!        --------------------
         do i = 1, MAX_OUTER_ITER
            call innerGMRES(this, ComputeTimeDerivative)
            if (this%ERROR_CODE .NE. 0) then
               PRINT*, 'ERROR IN GMRES, ERROR CODE: ', this%ERROR_CODE
               call this%destroy()
               error stop
            end if
            if(this%CONVERGED .OR. (this%niter .GE. this%maxiter)) then
               call this % MatrixAction(this % x, this % x0, ComputeTimeDerivative) ! using x0 as a temporary storing variable
               this % rnorm = norm2(this % RHS - this % x0)
               RETURN
            else
               this%x0 = this%x
            end if
         end do
         write(*,*) " ** WARNING!!: Reached max outer iter limit in GMRES ** "
          
!        Restore default values
!        ----------------------

         
      end subroutine GMRES_Solve
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     -----------------------------------------------------------------------------------------
!     Returns the preconditioning product   Pv = P⁻¹ * v for the GMRES recursive preconditioner
!     -----------------------------------------------------------------------------------------
      subroutine PC_GMRES_Ax(this,v, Pv, ComputeTimeDerivative)
         implicit none
         !---------------------------------------------------------
         class(MatFreeGMRES_t), intent(inout) :: this
         real(kind=RP)        , intent(in)    :: v(:)
         real(kind=RP)        , intent(out)   :: Pv(:)
         procedure(ComputeTimeDerivative_f)           :: ComputeTimeDerivative
         !---------------------------------------------------------
          
         CALL this % PCSolver % SetRHS(v)
         CALL this % PCSolver % Solve(NCONS, NGRAD, ComputeTimeDerivative,time = this % timesolve, dt = this % dtsolve)
         CALL this % PCSolver % Settol(1e-1_RP)                 ! TODO: aquí habrá que poner algo más elaborado
         Pv = this % PCSolver % x
!         n_preco_iter = n_preco_iter + PCSolver%niter     ! Not really needed
                        
      END subroutine PC_GMRES_Ax
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     --------------------------------------------------------------------------------------
!     Returns the preconditioning product   Pv = P⁻¹ * v for the Block-Jacobi preconditioner
!     --------------------------------------------------------------------------------------
      subroutine PC_BlockJacobi_Ax(this, v, Pv)
         implicit none
         !---------------------------------------------------------
         class(MatFreeGMRES_t), intent(inout) :: this
         real(kind=RP)        , intent(in)    :: v(:)
         real(kind=RP)        , intent(out)   :: Pv(:)
         !---------------------------------------------------------
         integer :: eID
         real(kind=RP), allocatable :: x(:)
         !---------------------------------------------------------
         
         call this % BlockPreco % SolveBlocks_LU(Pv,v)
         
      end subroutine PC_BlockJacobi_Ax
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ---------------------------------------------------
!     Returns the (matrix-free) matrix vector product 
!     ---------------------------------------------------
      subroutine GMRES_MatrixAction(this,x, Ax, ComputeTimeDerivative)
         implicit none
         !-arguments-----------------------------------------------
         class(MatFreeGMRES_t), intent(inout) :: this
         real(kind=RP), intent(in)            :: x (this % DimPrb)
         real(kind=RP), intent(out)           :: Ax(this % DimPrb)
         procedure(ComputeTimeDerivative_f)   :: ComputeTimeDerivative
         !-local-variables-----------------------------------------
         real(kind=RP) :: eps, shift
         !---------------------------------------------------------
         
         if (this % UserDef_Ax) then
            Ax = this % User_MatrixAction(x)
         else
            eps = 1e-8_RP * (1._RP + L2Norm(x) )
            shift = MatrixShift(this % dtsolve)
            
            Ax = ( this % p_F(this % Ur + x * eps, ComputeTimeDerivative) - this % F_Ur)/ eps  + shift * x     !First order 
            !Ax = ( this % p_F(this % Ur + x * eps, ComputeTimeDerivative) - this % p_F(Ur - x * eps, ComputeTimeDerivative))  /(2._RP * eps)  - x / Comp_Dt   !Second order
         end if
         
      end subroutine GMRES_MatrixAction
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ---------------------------------------------------
!     Function to return the time derivative using a specific solution vector
!     ---------------------------------------------------
      FUNCTION p_F(this,u, computeTimeDerivative) RESULT(F)
         IMPLICIT NONE
         !-arguments-----------------------------------------------
         class(MatFreeGMRES_t), INTENT(INOUT) :: this
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
      END FUNCTION p_F
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ------------------------------------------------
!     GMRES_SetMatrixAction:
!     Routine to set a user-defined Matrix-Action (Ax)
!     ------------------------------------------------
      subroutine GMRES_SetMatrixAction(this,MatrixAx)
         implicit none
         !-arguments-----------------------------------------------
         class(MatFreeGMRES_t), intent(inout) :: this
         procedure(MatrixAction_t)            :: MatrixAx
         !---------------------------------------------------------
         
         this % UserDef_Ax        = .TRUE.
         this % User_MatrixAction => MatrixAx
         
      end subroutine GMRES_SetMatrixAction

end module MatrixFreeGMRESClass