!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!      PetscSolverClass.f90
!      Created: 2017-04-10 10:006:00 +0100 
!      By: Andrés Rueda
!
!      Class for solving linear systems using MKL version of Pardiso
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
MODULE MKLPardisoSolverClass
   USE GenericLinSolverClass
   USE CSR_Matrices
   USE SMConstants
   USE PetscSolverClass   ! For allocating Jacobian matrix
   IMPLICIT NONE
   
   TYPE, EXTENDS(GenericLinSolver_t) :: MKLPardisoSolver_t
      TYPE(csrMat_t)                             :: A                                  ! Jacobian matrix
      REAL(KIND=RP), DIMENSION(:), ALLOCATABLE   :: x                                  ! Solution vector
      REAL(KIND=RP), DIMENSION(:), ALLOCATABLE   :: b                                  ! Right hand side
      REAL(KIND=RP)                              :: Ashift
      LOGICAL                                    :: AIsPrealloc   
      
      !Variables for creating Jacobian in PETSc context:
      TYPE(PetscKspLinearSolver_t)               :: PetscSolver                        ! PETSc solver (created only for allocating the matrix -to be deprecated?)
      LOGICAL                                    :: AIsPetsc = .TRUE.
      
      !Variables directly related with mkl pardiso solver
      INTEGER                                    :: mtype                              ! Matrix type. See construct
      INTEGER, ALLOCATABLE                       :: perm(:)
      INTEGER, POINTER                           :: Pardiso_iparm(:) => NULL()         ! Parameters for mkl version of pardiso
      INTEGER(KIND=AddrInt), POINTER             :: Pardiso_pt(:)    => NULL()  
   CONTAINS
      !Subroutines:
      PROCEDURE                                  :: construct => ConstructMKLContext
      PROCEDURE                                  :: PreallocateA
      PROCEDURE                                  :: ResetA
      PROCEDURE                                  :: SetAColumn 
      PROCEDURE                                  :: AssemblyA
      PROCEDURE                                  :: SetBValue
      PROCEDURE                                  :: solve
      PROCEDURE                                  :: GetCSRMatrix
      PROCEDURE                                  :: GetXValue
      PROCEDURE                                  :: destroy
      PROCEDURE                                  :: SetOperatorDt
      PROCEDURE                                  :: ReSetOperatorDt
      !Functions:
      PROCEDURE                                  :: Getxnorm    !Get solution norm
      PROCEDURE                                  :: Getrnorm    !Get residual norm
   END TYPE MKLPardisoSolver_t
   
   PRIVATE
   PUBLIC MKLPardisoSolver_t, GenericLinSolver_t
   
   
!========
 CONTAINS
!========
   
   SUBROUTINE ConstructMKLContext(this,DimPrb)
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(MKLPardisoSolver_t), INTENT(INOUT) :: this
      INTEGER                  , INTENT(IN)    :: DimPrb
      !-----------------------------------------------------------
      INTERFACE
         SUBROUTINE pardisoinit(pt, mtype, iparm)
            USE SMConstants
            IMPLICIT NONE
            INTEGER(KIND=AddrInt) :: pt(*)
            INTEGER :: mtype
            INTEGER :: iparm(*)
         END SUBROUTINE pardisoinit
      END INTERFACE
      !-----------------------------------------------------------
      
      this % DimPrb = DimPrb
      
      ALLOCATE(this % x(DimPrb))
      ALLOCATE(this % b(DimPrb))
      
      this % mtype = 11 !Set matrix type to real unsymmetric (change?)
    
      ALLOCATE(this % Pardiso_pt(64))
      ALLOCATE(this % Pardiso_iparm(64))
      
      ALLOCATE(this % perm(DimPrb))
      this % perm = 0
      
      IF(this % AIsPetsc) CALL this % PetscSolver % construct (DimPrb)
#ifdef HAS_MKL
      CALL pardisoinit(this % Pardiso_pt, this % mtype, this % Pardiso_iparm)
#else
      STOP 'MKL not linked correctly'
#endif
   END SUBROUTINE ConstructMKLContext
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE PreallocateA(this,nnz)
      IMPLICIT NONE
!
!  --------------------------------------------
!  Subroutine for preallocating Jacobian matrix
!     IMPORTANT:
!        - To date, it's using the PETSc library for preallocating due to its efficiency
!        - A proposed alternative to emulate PETSc allocation is to use 
!          a linked list matrix and then transform it to CSR format
!  --------------------------------------------
!
      !-----------------------------------------------------------
      CLASS(MKLPardisoSolver_t), INTENT(INOUT) :: this
      INTEGER                                  :: nnz
      !-----------------------------------------------------------
      LOGICAL, SAVE                            :: isfirst = .TRUE.
      !-----------------------------------------------------------
      
      
      IF (this % AIsPetsc) THEN
         CALL this % PetscSolver % PreallocateA(nnz)
      ELSE
         IF (.NOT. this % AIsPrealloc) THEN
            stop 'routines for preallocating matrix not implemented'
         ELSE
            print*, 'Matrix already preallocated'
         END IF
      END IF
      this % AIsPrealloc = .TRUE.
      
   END SUBROUTINE PreallocateA
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE ResetA(this)         
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(MKLPardisoSolver_t), INTENT(INOUT) :: this
      !-----------------------------------------------------------
      
      IF (this % AIsPetsc) THEN
         CALL this % PetscSolver % ResetA
      ELSE
         this % A % Values = 0._RP
      END IF
      
   END SUBROUTINE ResetA
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE SetAColumn(this,nvalues,irow,icol,values)         
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(MKLPardisoSolver_t), INTENT(INOUT)  :: this
      INTEGER       , INTENT(IN)                :: nvalues
      INTEGER       , INTENT(IN), DIMENSION(:)  :: irow
      INTEGER       , INTENT(IN)                :: icol
      REAL(KIND=RP) , INTENT(IN), DIMENSION(:)  :: values
      !-----------------------------------------------------------
      
      IF (this % AIsPetsc) THEN
         CALL this % PetscSolver % SetAColumn(nvalues,irow,icol,values)
      ELSE
         CALL this % A % SetColumn(irow+1,icol+1,values)
      END IF
      
   END SUBROUTINE SetAColumn
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE AssemblyA(this)         
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(MKLPardisoSolver_t), INTENT(INOUT) :: this
      !-----------------------------------------------------------
      
      IF (this % AIsPetsc) THEN
         CALL this % PetscSolver % AssemblyA
         CALL this % PetscSolver % GetCSRMatrix(this % A)
!~         CALL this % PetscSolver % destroy
!~         this % AIsPetsc = .FALSE.
      ELSE
         print*, 'A is assembled'
      END IF
      
   end SUBROUTINE AssemblyA
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE SetBValue(this, irow, value)
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(MKLPardisoSolver_t), INTENT(INOUT) :: this
      INTEGER                  , INTENT(IN)    :: irow
      REAL(KIND=RP)            , INTENT(IN)    :: value
      !-----------------------------------------------------------
      
      this % b (irow+1) = value
      
   END SUBROUTINE SetBValue
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE solve(this,tol,maxiter)         
      IMPLICIT NONE
!
!     ----------------------------------------------------
!     Main subroutine for solving system using mkl pardiso
!     ----------------------------------------------------
!
      !-----------------------------------------------------------
      CLASS(MKLPardisoSolver_t), INTENT(INOUT) :: this
      REAL(KIND=RP), OPTIONAL                  :: tol
      INTEGER      , OPTIONAL                  :: maxiter
      !-----------------------------------------------------------
#ifdef HAS_MKL
      INTEGER                                  :: error
      !-----------------------------------------------------------
      INTERFACE
         SUBROUTINE pardiso(pt, maxfct, mnum, mtype, phase, n, &
                           values, rows, cols, perm, nrhs, iparm, msglvl, b, x, ierror)
            USE SMConstants
            IMPLICIT NONE
            REAL(KIND=RP) :: values(*), b(*), x(*)
            INTEGER(KIND=AddrInt) :: pt(*)
            INTEGER :: perm(*), nrhs, iparm(*), msglvl, ierror
            INTEGER :: maxfct, mnum, mtype, phase, n, rows(*), cols(*)
         END SUBROUTINE pardiso
      END INTERFACE
      !-----------------------------------------------------------
      
!~    	call mkl_set_num_threads( 4 )
      
      !-----------------------
      ! Solve the system using MKL - Pardiso!!
!~    phase = 33            ! Solve, iterative refinement
      CALL pardiso(  pt      = this % Pardiso_pt    ,     &
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
                     msglvl  = 1                    ,     &     ! 1: verbose... Too much printing
                     b       = this % b             ,     &
                     x       = this % x             ,     &
                     ierror  = error              )
                
!~    msglvl    = 0 ! Do not write out any info
!~    nrhs      = 1 ! Use only one RHS

      IF (error .NE. 0) THEN
         WRITE(*,*) 'MKL Pardiso ERROR:', error
         STOP
      ELSE
         this % converged = .TRUE.
      END IF
    
#else
      STOP 'MKL is not linked properly'
#endif
      
   END SUBROUTINE solve
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE GetCSRMatrix(this,Acsr)         
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(MKLPardisoSolver_t), INTENT(IN)  :: this
      TYPE(csrMat_t)             , INTENT(OUT) :: Acsr 
      !-----------------------------------------------------------
      
      Acsr = this % A
   END SUBROUTINE GetCSRMatrix
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE GetXValue(this,irow,x_i)       
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(MKLPardisoSolver_t), INTENT(INOUT) :: this
      INTEGER                  , INTENT(IN)    :: irow
      REAL(KIND=RP)            , INTENT(OUT)   :: x_i
      !-----------------------------------------------------------
      
      x_i = this % x(irow+1)
      
   END SUBROUTINE GetXValue
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE destroy(this)       
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(MKLPardisoSolver_t), INTENT(INOUT) :: this
      !-----------------------------------------------------------
      
      CALL this % A % destroy()
      
      DEALLOCATE(this % b)
      DEALLOCATE(this % x)
      DEALLOCATE(this % Pardiso_pt)
      DEALLOCATE(this % Pardiso_iparm)
      DEALLOCATE(this % perm)
      
      IF (this % AIsPetsc) CALL this % PetscSolver % destroy()
      this % AIsPetsc    = .TRUE.
      this % AIsPrealloc = .FALSE.
      
    END SUBROUTINE destroy
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE SetOperatorDt(this,dt)       
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(MKLPardisoSolver_t), INTENT(INOUT) :: this
      REAL(KIND=RP)            , INTENT(IN)    :: dt
      !-----------------------------------------------------------
      
      this % Ashift = -1._RP/dt
      IF (this % AIsPetsc) THEN
         CALL this % PetscSolver % SetOperatorDt(dt)
      ELSE
         CALL this % A % SetMatShift(this % Ashift)
      END IF
      
    END SUBROUTINE SetOperatorDt
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE ReSetOperatorDt(this,dt)       
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(MKLPardisoSolver_t), INTENT(INOUT) :: this
      REAL(KIND=RP)            , INTENT(IN)    :: dt
      !-----------------------------------------------------------
      REAL(KIND=RP)                            :: shift
      !-----------------------------------------------------------
      
      shift = -1._RP/dt
      IF (this % AIsPetsc) THEN
         CALL this % PetscSolver % ReSetOperatorDt(dt)
      ELSE
         CALL this % A % SetMatShift(-this % Ashift)
         CALL this % A % SetMatShift(shift)
      END IF
      this % Ashift = shift
      
    END SUBROUTINE ReSetOperatorDt
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   FUNCTION Getxnorm(this,TypeOfNorm) RESULT(xnorm)
      IMPLICIT NONE
      !-----------------------------------------------------------
      CLASS(MKLPardisoSolver_t), INTENT(INOUT) :: this
      CHARACTER(len=*)                         :: TypeOfNorm
      REAL(KIND=RP)                            :: xnorm
      !-----------------------------------------------------------
      
      SELECT CASE (TypeOfNorm)
         CASE ('infinity')
            xnorm = MAXVAL(ABS(this % x))
         CASE ('l2')
            xnorm = NORM2(this % x)
         CASE DEFAULT
            STOP 'MKLPardisoSolverClass ERROR: Norm not implemented yet'
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
      CLASS(MKLPardisoSolver_t), INTENT(INOUT) :: this
      REAL(KIND=RP)                            :: rnorm
      !-----------------------------------------------------------
      REAL(KIND=RP)                            :: residual(this % DimPrb)
      !-----------------------------------------------------------
      
      residual = this % b - CSR_MatVecMul(this % A, this % x)
      rnorm = MAXVAL(ABS(residual))
      
      
      !rnorm = NORM2(this % x)
      
   END FUNCTION Getrnorm
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
END MODULE MKLPardisoSolverClass
