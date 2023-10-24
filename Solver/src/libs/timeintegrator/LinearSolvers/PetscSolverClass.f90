!//////////////////////////////////////////////////////
!
!   Class for solving linear systems using the Krylov Subspace Methods of PETSc library
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef HAS_PETSC
#include "petsc/finclude/petsc.h"
#endif
module PetscSolverClass
   use GenericLinSolverClass
   use MatrixClass   
   use SMConstants
   use DGSEMClass             , only: DGSem, computetimederivative_f
   use MPI_Process_Info       , only: MPI_Process
#ifdef HAS_PETSC
   use petsc
#endif
#ifdef _HAS_MPI_
   use mpi
#endif 
   implicit none
   
   type, extends(GenericLinSolver_t) :: PetscKspLinearSolver_t
      type(PETSCMatrix_t)                           :: A
      character(len=LINE_LENGTH)                    :: preconditioner
#ifdef HAS_PETSC
      Vec                                           :: x                                  ! Solution vector
      Vec                                           :: b                                  ! Right hand side
      KSP                                           :: ksp                                ! 
      PC                                            :: pc
      PetscReal                                     :: tol = 1d-15
      PetscInt                                      :: maxiter = 50
      PetscInt                                      :: nz = 0
      PetscScalar                                   :: Ashift                              ! Stores the shift to the Jacobian due to time integration
      PetscBool                                     :: init_context = PETSC_FALSE
#endif
      CONTAINS
         !Subroutines
         procedure :: construct           => PETSc_construct
         procedure :: SetRHSValues        => PETSc_SetRHSValues
         procedure :: SetRHSValue         => PETSc_SetRHSValue
         procedure :: GetXValues          => PETSc_GetXValues
         procedure :: GetXValue           => PETSc_GetXValue
         procedure :: GetX                => PETSc_GetX
         procedure :: SetOperatorDt       => PETSc_SetOperatorDt
         procedure :: ReSetOperatorDt     => PETSc_ReSetOperatorDt
         procedure :: AssemblyRHS         => PETSc_AssemblyRHS
         procedure :: SaveMat
         procedure :: solve               => PETSc_solve
         procedure :: destroy             => PETSc_Destroy
         procedure :: SetRHS              => PETSc_SetRHS
         procedure :: SetPreconditioner   => PETSc_SetPreconditioner
         !Functions
         procedure :: GetXnorm         => PETSc_GetXnorm
         procedure :: GetRnorm         => PETSc_GetRnorm
   end type PetscKspLinearSolver_t
   
  
   private                                          
   public                                           :: PetscKspLinearSolver_t
   
   

!========
 CONTAINS
!========
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine CheckPetscErr(ierr,msg)
      implicit none
      !-arguments-----------------------------------------------------------
      character(LEN=*), optional                   :: msg
#ifdef HAS_PETSC
      PetscErrorCode, intent(in)                   :: ierr
      !---------------------------------------------------------------------
      
      if(ierr .EQ. 0) then
         return
      else
         if (.NOT. PRESENT(msg)) msg = 'error in petsc'
         write(*,*) msg,' **** Petsc call returned an error. Code: ' ,ierr
         error stop
      end if
#else
      integer                                      :: ierr
      error stop ':: PETSc is not linked correctly'
#endif
   end subroutine CheckPetscErr
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine PETSc_construct(this, DimPrb, globalDimPrb, nEqn, controlVariables,sem,MatrixShiftFunc)
      implicit none
      !-arguments-----------------------------------------------------------
      class(PetscKspLinearSolver_t), intent(inout), TARGET :: this
      integer                      , intent(in)            :: nEqn
      type(FTValueDictionary)      , intent(in), optional  :: controlVariables
      type(DGSem), TARGET                      , optional  :: sem
      procedure(MatrixShift_FCN)                           :: MatrixShiftFunc
#ifdef HAS_PETSC
      PetscInt, intent(in)                                 :: DimPrb
      PetscInt, intent(in)                                 :: globalDimPrb
      !-local-variables-----------------------------------------------------
      PetscErrorCode                                       :: ierr
      !---------------------------------------------------------------------
      
      call this % GenericLinSolver_t % construct(DimPrb,globalDimPrb, nEqn,controlVariables,sem,MatrixShiftFunc)
      
      if ( this % withMPI .and. (globalDimPrb /= MPI_Process % nprocs * DimPrb)) then
         print*, "IMPORTANT WARNING (PetscSolverClass)"
         print*,  "-> There is a problem (likely a BUG) when the MPI partitions don't \n",     &
                 "    have the exact same number of degrees of freedom \n",                    &
                 " -> This simulation will probably crash \n",                                 &
                 " -> To make it work, make sure the number of partitions is a divisor \n",    &
                 "    of the number of elements (of course, if all elements have the same DOFs) \n", &
                 "            ... or fix the bug"
      end if
      
      MatrixShift => MatrixShiftFunc
      
      !Initialisation of the PETSc variables
      ! call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-ksp_gmres_restart","100",ierr)
      ! call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-ksp_gmres_modifiedgramschmidt","true",ierr)
      ! call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-pc_bjacobi_blocks",size(sem % mesh % elements),ierr)
      call PetscInitialize(PETSC_NULL_character,ierr)

!     PETSc matrix A 
      call this % A % construct(num_of_Rows = DimPrb, num_of_TotalRows = globalDimPrb)
      
      if ( present(sem) ) then
         this % p_sem => sem
         call this % Jacobian % Configure (sem % mesh, nEqn, this % A)
      end if
      
!     Petsc vectors x and b (of A x = b)
      
      if (this % withMPI) then ! Only possible if this % A is preallocated
         call MatCreateVecs(this % A % A, this % x, this % b,ierr) ; call CheckPetscErr(ierr,'error creating MPI Petsc vector')
         
!~         call VecCreateMPI(PETSC_COMM_WORLD,dimPrb,globalDimPrb,this % x,ierr)
!~         call CheckPetscErr(ierr,'error creating MPI Petsc vector')
      else
         call VecCreate  (PETSC_COMM_WORLD,this % x,ierr)          ; call CheckPetscErr(ierr,'error creating Petsc vector')
         call VecSetSizes(this % x,dimPrb,globalDimPrb,ierr)       ; call CheckPetscErr(ierr,'error setting Petsc vector options')
         call VecSetFromOptions(this % x,ierr)                     ; call CheckPetscErr(ierr,'error setting Petsc vector options')
         call VecDuplicate(this % x,this % b,ierr)                 ; call CheckPetscErr(ierr,'error creating Petsc vector')
      end if

!     Petsc ksp solver context      
      call KSPCreate(PETSC_COMM_WORLD,this%ksp,ierr)                    ; call CheckPetscErr(ierr,'error in KSPCreate')
      call KSPSetFromOptions(this%ksp,ierr) ! debug
!     Petsc preconditioner 
      call KSPGetPC(this%ksp,this%pc,ierr)                              ; call CheckPetscErr(ierr,'error in KSPGetPC')
      
      if ( controlVariables % containsKey("preconditioner") ) then
         this % preconditioner = controlVariables % stringValueForKey("preconditioner",LINE_LENGTH)
      else
         this % preconditioner = 'Block-Jacobi'
      end if

!     Finish up
!     ---------
      
      this%init_context = PETSC_TRUE
      this%dimprb = DimPrb
      
#else
      integer, intent(in)                       :: DimPrb
      integer, intent(in)                       :: globalDimPrb
      error stop ':: PETSc is not linked correctly'
#endif
   end subroutine PETSc_construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine PETSc_SetPreconditioner(this)
      implicit none
      !-arguments-----------------------------------------------------------
      class(PetscKspLinearSolver_t), intent(inout)      :: this
#ifdef HAS_PETSC
      !-local-variables-----------------------------------------------------
      PetscErrorCode                                  :: ierr
      !---------------------------------------------------------------------
!      
!     Set preconditioner settings in KSP (this only has to be done once in theory, but it's needed when the matrix is reconstructed (like in the static-condensation solver)
!     ----------------------------------
      select case ( trim(this % preconditioner) )
         case ('Block-Jacobi')
            
            call MatSetVariableBlockSizes (this % A % A, size(this % Jacobian % ndofelm_l), this % Jacobian % ndofelm_l(1), ierr)  ; call CheckPetscErr(ierr, 'error in MatSetVariableBlockSizes')     ! PCVPBJACOBI
            call PCSetType(this%pc,PCVPBJACOBI,ierr)                 ; call CheckPetscErr(ierr, 'error in PCSetType')
            
         case ('Jacobi')
            
            call PCSetType(this%pc,PCJACOBI,ierr)                 ; call CheckPetscErr(ierr, 'error in PCSetType')
         case ('ILU')
            
            call PCSetType(this%pc,PCILU,ierr)                 ; call CheckPetscErr(ierr, 'error in PCSetType') 
         case default
         
            error stop 'PETSc_SetPreconditioner: Not recognized preconditioner'
      end select
!      
!     Set operators for KSP
!     ---------------------
      call KSPSetOperators(this%ksp, this%A%A, this%A%A, ierr)     ; call CheckPetscErr(ierr, 'error in KSPSetOperators')
      
!~      call PCSetType(this%pc,PCILU,ierr)       ;call CheckPetscErr(ierr, 'error in PCSetType')
#else
      error stop ':: PETSc is not linked correctly'
#endif
   end subroutine PETSc_SetPreconditioner 
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine PETSc_solve(this, nEqn, nGradEqn, ComputeTimeDerivative, tol, maxiter, time,dt, ComputeA)
      implicit none
      !-arguments-----------------------------------------------------------
      class(PetscKspLinearSolver_t), target, intent(inout)  :: this
      integer,       intent(in)                             :: nEqn, nGradEqn
      procedure(ComputeTimeDerivative_f)                    :: ComputeTimeDerivative
      real(kind=RP), optional                               :: time
      real(kind=RP), optional                               :: dt
      logical      , optional               , intent(inout) :: ComputeA
#ifdef HAS_PETSC
      PetscReal    , optional                               :: tol
      PetscInt     , optional                               :: maxiter
      !-local-variables-----------------------------------------------------
      PetscErrorCode                                  :: ierr
      PetscInt                                        :: nresvec=500
      PetscReal ,dimension(500)                       :: resvec
      type(csrMat_t) :: Afull ! to save & visualize the matrix if needed
      !---------------------------------------------------------------------
      
      if ( present(ComputeA)) then
         if (ComputeA) then
            call this % Jacobian % Compute (this % p_sem, nEqn, time, this % A, ComputeTimeDerivative)
            ! call this % A % GetCSRMatrix (Afull) ! FINDME!
            ! call Afull % Visualize('Afull_f.txt') ! FINDME1
            call this % SetOperatorDt(dt)
            ComputeA = .FALSE.
            
            call this % SetPreconditioner
         end if
      else 
         call this % Jacobian % Compute (this % p_sem, nEqn, time, this % A, ComputeTimeDerivative)
         ! call this % A % GetCSRMatrix (Afull) ! FINDME1
         ! call Afull % Visualize('Afull_f.txt') ! FINDME1
         call this % SetOperatorDt(dt)
         
         call this % SetPreconditioner
      end if
      
      ! call this % A % GetCSRMatrix (Afull)
      ! call Afull % Visualize('Afull_f.txt') ! visualize
      
      ! Set , if given, solver tolerance and max number of iterations
      if (PRESENT(tol)) then
         this%tol = tol
      else
         this%tol = PETSC_DEFAULT_real
      end if
      
      if (PRESENT(maxiter)) then
         this%maxiter = maxiter
      else
         this%maxiter = PETSC_DEFAULT_integer
      end if
      
      call KSPSetTolerances(this%ksp,this % tol,this%tol,PETSC_DEFAULT_REAL,this%maxiter,ierr)
      call CheckPetscErr(ierr, 'error in KSPSetTolerances')
      
      ! Set initial guess to P⁻¹b
      call KSPSetInitialGuessKnoll(this%ksp,PETSC_TRUE,ierr)
      call CheckPetscErr(ierr, 'error in KSPSetInitialGuessKnoll')

      ! set vector for residual history
      call KSPSetResidualHistory(this%ksp,resvec,nresvec,PETSC_TRUE,ierr)
      call CheckPetscErr(ierr, 'error in KSPSetResidualHistory')

      ! set type of solver
      call KSPSetType(this % ksp,KSPGMRES,ierr) ; call CheckPetscErr(ierr, 'error in KSetType')
      ! call KSPSetType(this % ksp,KSPRICHARDSON,ierr) ; call CheckPetscErr(ierr, 'error in KSetType')
      
      call KSPSolve(this%ksp,this%b,this%x,ierr)               ; call CheckPetscErr(ierr, 'error in KSPSolve')

      ! get residual vector
      call KSPGetResidualHistory(this%ksp,resvec,nresvec,ierr)
      call CheckPetscErr(ierr, 'error in KSPGetResidualHistory')
      
      call KSPGetIterationNumber(this%ksp,this%niter,ierr)     ; call CheckPetscErr(ierr,'error in KSPGetIterationNumber')
!~       call KSPGetResidualNorm(this%ksp, this%residual, ierr)   ; call CheckPetscErr(ierr,'error in KSPGetResidualNorm')
!~       call VecNorm(this%x,NORM_INFINITY,this%xnorm,ierr)       ; call CheckPetscErr(ierr,'error in VecNorm')

!~      ! print stuff 
!~      print *, "No iterations: ", this % niter
!~      print *, "Residual history: "
!~      write(*,"(ES14.7)") resvec(1:this%niter)

      if (this%niter < maxiter) then
         this%converged = .TRUE.
      else
         this%converged = .FALSE.
      end if
      
#else
      real*8 , optional                     :: tol
      integer, optional                     :: maxiter
      error stop ':: PETSc is not linked correctly'
#endif
   end subroutine PETSc_solve
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine PETSc_SetOperatorDt(this, dt)
      implicit none
      !-arguments-----------------------------------------------------------
      class(PetscKspLinearSolver_t),     intent(inout)     :: this
#ifdef HAS_PETSC
      PetscScalar,                     intent(in)        :: dt
      !-local-variables-----------------------------------------------------
      PetscScalar                                        :: shift
      PetscScalar                                        :: eps = 1e-10
      !---------------------------------------------------------------------
      
      shift = MatrixShift(dt) !
      if (ABS(shift) .GT. eps) then                  
         call this % A % shift(shift) ! A = A + shift * I
         this % Ashift = shift
      end if
#else
      real*8,                     intent(in)        :: dt
      error stop ':: PETSc is not linked correctly'
#endif
   end subroutine PETSc_SetOperatorDt
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------------------
!  Removes previous shift in order to insert new one 
!              (important when Jacobian is reused)
!  --------------------------------------------------
   subroutine PETSc_ReSetOperatorDt(this, dt)
      implicit none
      !-arguments-----------------------------------------------------------
      class(PetscKspLinearSolver_t),     intent(inout)     :: this
#ifdef HAS_PETSC
      PetscScalar,                     intent(in)        :: dt
      !-local-variables-----------------------------------------------------
      PetscScalar                                        :: shift
      PetscScalar                                        :: eps = 1e-10
      !---------------------------------------------------------------------
      
      shift = MatrixShift(dt) !
      if (ABS(shift) .GT. eps) then
         call this % A % Reshift (shift) ! A = A + shift * I
         this % Ashift = shift
      end if
#else
      real*8,                     intent(in)        :: dt
      error stop ':: PETSc is not linked correctly'
#endif
   end subroutine PETSc_ReSetOperatorDt
!
!/////////////////////////////////////////////////////////////////////////////////////////////////   
!
   subroutine PETSc_SetRHSValues(this, nvalues, irow, values)
      implicit none
      !-arguments-----------------------------------------------------------
      class(PetscKspLinearSolver_t),     intent(inout)     :: this
#ifdef HAS_PETSC
      PetscInt,                        intent(in)        :: nvalues
      PetscInt, dimension(:),          intent(in)        :: irow
      PetscScalar, dimension(:),       intent(in)        :: values
      !-local-variables-----------------------------------------------------
      PetscErrorCode                                     :: ierr
      !---------------------------------------------------------------------
      
      call VecSetValues(this%b,nvalues, irow-1,values,INSERT_VALUES, ierr)
      call CheckPetscErr(ierr, 'error in VecSetValues')
#else
      integer,                        intent(in)        :: nvalues
      integer, dimension(:),          intent(in)        :: irow
      real*8    , dimension(:),       intent(in)        :: values
      error stop ':: PETSc is not linked correctly'
#endif
   end subroutine PETSc_SetRHSValues
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine PETSc_SetRHSValue(this, irow, value)
      implicit none
      !-arguments-----------------------------------------------------------
      class(PetscKspLinearSolver_t),     intent(inout)     :: this
#ifdef HAS_PETSC
      PetscInt,                        intent(in)        :: irow
      PetscScalar,                     intent(in)        :: value
      !-local-variables-----------------------------------------------------
      PetscErrorCode                                     :: ierr
      !---------------------------------------------------------------------
      
      call VecSetValue(this%b, irow-1,value,INSERT_VALUES, ierr)
      call CheckPetscErr(ierr, 'error in VecSetValues')
#else
      integer,           intent(in)        :: irow
      real*8    ,        intent(in)        :: value
      error stop ':: PETSc is not linked correctly'
#endif
   end subroutine PETSc_SetRHSValue
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine PETSc_SetRHS(this, RHS)
      implicit none
      !-arguments-----------------------------------------------------------
      class(PetscKspLinearSolver_t), intent(inout) :: this
#ifdef HAS_PETSC
      PetscScalar                  , intent(in)    :: RHS(this % DimPrb)
      !-local-variables-----------------------------------------------------
      integer              :: i, counter, ndof
      integer, allocatable :: ind(:)
      integer              :: eID, globID
      PetscErrorCode       :: ierr
      PetscInt :: ranges(this % DimPrb + 1)
      !---------------------------------------------------------------------
      
      if (this % withMPI) then ! Assuming the Jacobian was constructed
         counter = 1
         do eID = 1, size(this % Jacobian % globIDs_l)
            globID = this % Jacobian % globIDs_l (eID)
            ndof   = this % Jacobian % ndofelm(globID)
            
            allocate ( ind (ndof) )
            ind = [(i, i=this % Jacobian % firstIdx(globID)      - 1 , &
                         this % Jacobian % firstIdx(globID + 1 ) - 2) ]
            call VecSetValues  (this%b, ndof, ind, RHS(counter:counter+ndof-1), INSERT_VALUES, ierr)
            counter = counter + ndof
            deallocate(ind)
         end do
      else
         call VecSetValues  (this%b, this % DimPrb, [(i, i=0, this % DimPrb-1)] , RHS, INSERT_VALUES, ierr)
         call CheckPetscErr(ierr, 'error in VecSetValues')
      end if
      
#else
      real(kind=RP)                , intent(in)    :: RHS(this % DimPrb)
      error stop ':: PETSc is not linked correctly'
#endif
   end subroutine PETSc_SetRHS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine PETSc_AssemblyRHS(this)
      implicit none
      !-arguments-----------------------------------------------------------
      class(PetscKspLinearSolver_t),     intent(inout)   :: this
      !-local-variables-----------------------------------------------------
#ifdef HAS_PETSC
      PetscErrorCode                                     :: ierr
      !---------------------------------------------------------------------
      
      call VecAssemblyBegin(this%b, ierr);  call CheckPetscErr(ierr," Assembly B in PETSc Begin")      
      call VecAssemblyEnd(this%b, ierr)  ;  call CheckPetscErr(ierr," Assembly B in PETSc End")  
#else
      error stop ':: PETSc is not linked correctly'
#endif
   end subroutine PETSc_AssemblyRHS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine PETSc_GetXValues(this, nvalues, irow, values)
      implicit none
      !-arguments-----------------------------------------------------------
      class(PetscKspLinearSolver_t),     intent(inout)      :: this
#ifdef HAS_PETSC
      PetscInt,                        intent(in)         :: nvalues
      PetscInt, dimension(:),          intent(in)         :: irow
      PetscScalar, dimension(:),       intent(out)        :: values
      !-local-variables-----------------------------------------------------
      PetscErrorCode                                      :: ierr
      !---------------------------------------------------------------------
      
      call VecGetValues(this%x,nvalues,irow-1,values, ierr)
      call CheckPetscErr(ierr, 'error in VecGetValues')
#else
      integer,                        intent(in)        :: nvalues
      integer, dimension(:),          intent(in)        :: irow
      real*8    , dimension(:),       intent(in)        :: values
      error stop ':: PETSc is not linked correctly'
#endif
   end subroutine PETSc_GetXValues
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine PETSc_GetXValue(this, irow, x_i)
      implicit none
      !-arguments-----------------------------------------------------------
      class(PetscKspLinearSolver_t),     intent(inout)      :: this
#ifdef HAS_PETSC
      PetscInt,                        intent(in)         :: irow
      PetscScalar,                     intent(out)        :: x_i
      !-local-variables-----------------------------------------------------
      PetscErrorCode                                      :: ierr
      PetscInt    :: i(1)
      PetscScalar :: x(1)
      !---------------------------------------------------------------------
      
      i = irow-1
      
      call VecGetValues(this%x, 1, i, x, ierr)  ! TODO: Fix problem here?
      call CheckPetscErr(ierr, 'error in VecGetValue')
      
      x_i = x(1)
#else
      integer,           intent(in)        :: irow
      real*8    ,        intent(out)       :: x_i
      error stop ':: PETSc is not linked correctly'
#endif
   end subroutine PETSc_GetXValue
!
!//////////////////////////////////////////////////////////////////////////////////////////////////
!   
   function PETSc_GetX(this) result(x)
      implicit none
      !-arguments-----------------------------------------------------------
      class(PetscKspLinearSolver_t),     intent(inout)      :: this
#ifdef HAS_PETSC
      PetscScalar                                           :: x(this % DimPrb)
      !-local-variables-----------------------------------------------------
      PetscInt             :: irow(this % DimPrb)
      PetscErrorCode       :: ierr
      PetscScalar, pointer :: xout(:)
!~      integer :: ndof, counter, eID, i, globID
!~      integer, allocatable :: ind(:)
      !------------------------------------------
      
!~      if (this % withMPI) then ! Assuming the Jacobian was constructed
!~         counter = 1
!~         do eID = 1, size(this % Jacobian % globIDs_l)
!~            globID = this % Jacobian % globIDs_l (eID)
!~            ndof   = this % Jacobian % ndofelm(globID)
            
!~            allocate ( ind (ndof) )
!~            ind = [(i, i=this % Jacobian % firstIdx(globID)      - 1 , &
!~                         this % Jacobian % firstIdx(globID + 1 ) - 2) ]
!~            call VecGetValues  (this%x, ndof, ind, x(counter:counter+ndof-1), ierr)
!~            call CheckPetscErr(ierr, 'error in VecGetValue')
!~            counter = counter + ndof
!~            deallocate(ind)
!~         end do
!~      else
!~         irow = (/ (i, i=0, this % DimPrb-1) /)
!~         call VecGetValues(this%x,this % DimPrb ,irow,x, ierr) ; call CheckPetscErr(ierr, 'error in VecGetValue')
!~      end if
      
      ! TODO: Check if this works for non-consecutive (in the numbering) elements in a single partition
      call VecGetArrayReadF90(this%x,xout,ierr)
      x = xout
      call VecRestoreArrayReadF90(this%x,xout,ierr)
      
#else
      integer         :: irow
      real(kind=RP)                        :: x(this % DimPrb)
      error stop ':: PETSc is not linked correctly'
#endif
   end function PETSc_GetX
!
!//////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SaveMat(this,filename)
      implicit none
      !-arguments-----------------------------------------------------------
      class(PetscKspLinearSolver_t), intent(inout)         :: this
      character(LEN=*), optional                         :: filename
#ifdef HAS_PETSC
      !-local-variables-----------------------------------------------------
      PetscViewer                                        :: viewer
      PetscErrorCode                                     :: ierr
      !---------------------------------------------------------------------
      
      !call MatView(this % A % A,PETSC_VIEWER_DRAW_SELF)
!~      read(*,*)
!~       if (.NOT. PRESENT(filename)) filename = &
!~                             '/home/andresrueda/Dropbox/PhD/03_Initial_Codes/3D/Implicit/nslite3d/Tests/Euler/NumJac/MatMatlab.dat'
!~       call PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename , viewer, ierr)    ; call CheckPetscErr(ierr)
!~       call PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB , ierr)      ; call CheckPetscErr(ierr)
!~       call MatView(this%A, viewer, ierr)                                      ; call CheckPetscErr(ierr)
!~       call PetscViewerDestroy(viewer, ierr)                                   ; call CheckPetscErr(ierr)
#else
      error stop ':: PETSc is not linked correctly'
#endif
   end subroutine SaveMat
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine PETSc_Destroy(this)
      implicit none
      !-arguments-----------------------------------------------------------
      class(PetscKspLinearSolver_t), intent(inout)       :: this
#ifdef HAS_PETSC
      !-local-variables-----------------------------------------------------
      PetscErrorCode                                   :: ierr1, ierr2, ierr3, ierr4
      !---------------------------------------------------------------------
      call VecDestroy(this%x,ierr1)
      call VecDestroy(this%b,ierr2)
      
      call this % A % destruct
      call KSPDestroy(this%ksp,ierr3)  ! this % pc is destructed inside
      call PetscFinalize(ierr4)
      
      call CheckPetscErr(ierr1,'error in VecDestroy x')
      call CheckPetscErr(ierr2,'error in VecDestroy b')
      call CheckPetscErr(ierr3,'error in KSPDestroy')
      call CheckPetscErr(ierr4,'error in PetscFinalize')
#else
      error stop ':: PETSc is not linked correctly'
#endif
   end subroutine PETSc_Destroy
!
!////////////////////////////////////////////////////////////////////////////////////////////////// 
!
   function PETSc_GetXnorm(this,TypeOfNorm) RESULT(xnorm)
      implicit none
      !-arguments-----------------------------------------------------------
      class(PetscKspLinearSolver_t), intent(inout) :: this
      character(len=*)                             :: TypeOfNorm
      
      !-local-variables-----------------------------------------------------
#ifdef HAS_PETSC
      PetscScalar :: xnorm
      PetscErrorCode                               :: ierr
      !--------------------------------------------------------------
      
      select case(TypeOfNorm)
         case('infinity')
            call VecNorm(this%x,NORM_INFINITY,xnorm,ierr)       ; call CheckPetscErr(ierr,'error in VecNorm')
         case default
            error stop 'PetscSolverClass error: Type of Norm not defined'
      end select
#else
      real(kind=RP)                                :: xnorm
      error stop ':: PETSc is not linked correctly'
#endif
   end function PETSc_GetXnorm
!
!////////////////////////////////////////////////////////////////////////////////////////////////// 
!
   function PETSc_GetRnorm(this) RESULT(rnorm)
      implicit none
      !-arguments-----------------------------------------------------------
      class(PetscKspLinearSolver_t), intent(inout) :: this
      real(kind=RP)                                :: rnorm
      !-local-variables-----------------------------------------------------
#ifdef HAS_PETSC
      PetscErrorCode                               :: ierr
      !--------------------------------------------------------------
      
      ! I don't know which type of norm PETSc computes!
      call KSPGetResidualNorm(this%ksp, rnorm, ierr)   ; call CheckPetscErr(ierr,'error in KSPGetResidualNorm')
#else
      error stop ':: PETSc is not linked correctly'
#endif
   end function PETSc_GetRnorm
end module PetscSolverClass