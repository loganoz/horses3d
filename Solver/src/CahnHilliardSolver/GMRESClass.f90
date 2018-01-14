!
!//////////////////////////////////////////////////////
!
!   @File:    GMRESClass.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Sun Jan 14 13:23:02 2018
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
MODULE GMRESClass
   USE SMConstants,                 ONLY: RP
   IMPLICIT NONE
   
   TYPE GmresSolver
      INTEGER                                :: m = 60                 ! arueda: Number of GMRES iterations before restart (?) -- Default petsc value m=30 
      INTEGER                                :: psize
      INTEGER                                :: maxiter = 500
      INTEGER                                :: niter 
      REAL(KIND = RP)                        :: tol = 1e-15_RP
      REAL(KIND = RP)                        :: res = -1._RP
      LOGICAL                                :: CONVERGED = .FALSE.
      LOGICAL                                :: FLEXIBLE = .FALSE.
      INTEGER                                :: ERROR_CODE = 0
      PROCEDURE(matmultsub), POINTER, NOPASS :: MatMult => NULL()
      PROCEDURE(matmultsub), POINTER, NOPASS :: PCMult => NULL()
      REAL(KIND = RP)                        :: norm0 = -1._RP
      REAL(KIND = RP), ALLOCATABLE           :: RHS(:)
      REAL(KIND = RP), ALLOCATABLE           :: x(:)
      REAL(KIND = RP), ALLOCATABLE           :: x0(:)
      
      REAL(KIND = RP), ALLOCATABLE, PRIVATE  :: H(:,:)
      REAL(KIND = RP), ALLOCATABLE, PRIVATE  :: W(:)
      REAL(KIND = RP), ALLOCATABLE, PRIVATE  :: V(:,:)                ! arueda: Orthogonal vectors of Krylov subspace (Arnoldi)
      REAL(KIND = RP), ALLOCATABLE, PRIVATE  :: Z(:,:)
      REAL(KIND = RP), ALLOCATABLE, PRIVATE  :: Y(:)
      REAL(KIND = RP), ALLOCATABLE, PRIVATE  :: cc(:)
      REAL(KIND = RP), ALLOCATABLE, PRIVATE  :: ss(:)
      REAL(KIND = RP), ALLOCATABLE, PRIVATE  :: g(:)
      CONTAINS
         PROCEDURE                           :: Construct => ConstructSolver
         PROCEDURE                           :: Destruct  => DestructSolver
         PROCEDURE                           :: SetRHS
         PROCEDURE                           :: SetMatMult
         PROCEDURE                           :: SetPCMult
         PROCEDURE                           :: SetTol
         PROCEDURE                           :: SetMaxIter
         PROCEDURE                           :: Solve => SolveGMRES
         PROCEDURE                           :: SetInitialGuess
   END TYPE GmresSolver

   ABSTRACT INTERFACE
      SUBROUTINE matmultsub(v,x)
         IMPORT RP
         REAL(KIND = RP), INTENT(IN)         :: v(:)
         REAL(KIND = RP), INTENT(OUT)        :: x(:)
      END SUBROUTINE
   END INTERFACE
   


   PRIVATE
   PUBLIC   :: GmresSolver 
   
   CONTAINS
 !/////////////////////////////////////////////////////////////////////////
      SUBROUTINE ConstructSolver(this,psize,m)
         CLASS(GmresSolver), INTENT(INOUT)       :: this
         INTEGER, INTENT(IN)                     :: psize
         INTEGER, INTENT(IN),OPTIONAL            :: m
         
         IF(PRESENT(m)) this%m = m
         this%psize = psize
         ALLOCATE(this%RHS(psize))
         ALLOCATE(this%x0(psize))
         this%x0 = 0.0_RP
         ALLOCATE(this%V(psize,this%m+1))
         ALLOCATE(this%H(this%m+1,this%m))
         ALLOCATE(this%W(psize))
         ALLOCATE(this%Y(this%m))
         ALLOCATE(this%cc(this%m+1))
         ALLOCATE(this%ss(this%m+1))
         ALLOCATE(this%g(this%m+1))
      END SUBROUTINE ConstructSolver
 !/////////////////////////////////////////////////////////////////////////      
      SUBROUTINE DestructSolver(this)
         CLASS(GmresSolver), INTENT(INOUT)          :: this
         
         DEALLOCATE(this%RHS)
         DEALLOCATE(this%x0)
         DEALLOCATE(this%V)
         DEALLOCATE(this%H)
         DEALLOCATE(this%W)
         DEALLOCATE(this%Y)
         DEALLOCATE(this%cc)
         DEALLOCATE(this%ss)
         DEALLOCATE(this%g)
         this%MatMult => NULL()
         IF (ALLOCATED(this%Z)) DEALLOCATE(this%Z)
      END SUBROUTINE DestructSolver
  !/////////////////////////////////////////////////////////////////////////     
      SUBROUTINE SetTol(this,tol)
         CLASS(GmresSolver), INTENT(INOUT)      :: this
         REAL(KIND = RP)                        :: tol
         this%tol = MAX(tol, EPSILON(1._RP)) ! This limits the solver tolerance to the used real precicision
         
      END SUBROUTINE SetTol
 !/////////////////////////////////////////////////////////////////////////      
      SUBROUTINE SetMaxIter(this,maxiter)
         CLASS(GmresSolver), INTENT(INOUT)          :: this
         INTEGER                                :: maxiter
         
         this%maxiter = maxiter
      END SUBROUTINE SetMaxIter
 !/////////////////////////////////////////////////////////////////////////      
      SUBROUTINE SetMatMult(this,MatMult)
         CLASS(GmresSolver), INTENT(INOUT)      :: this
         PROCEDURE(matmultsub)                  :: MatMult
         !EXTERNAL                               :: MatMult
         
         this%MatMult => MatMult
      END SUBROUTINE  SetMatMult
 !/////////////////////////////////////////////////////////////////////////      
      SUBROUTINE SetPCMult(this,PCMult)
         CLASS(GmresSolver), INTENT(INOUT)      :: this
         PROCEDURE(matmultsub)                  :: PCMult
         !EXTERNAL                               :: PCMult
         
         ! Sets preconditioning Matrix-vector multiplication 
         ! and allotates extra array
         
         IF (.NOT. ALLOCATED(this%Z)) ALLOCATE(this%Z(this%psize,this%m+1))
         this%PCMult => PCMult
         this%FLEXIBLE = .TRUE.
         
      END SUBROUTINE  SetPCMult
!/////////////////////////////////////////////////////////////////////////      
      SUBROUTINE SetRHS(this,RHS) ! TODO: by doing things in this way vector RHS is duplicated, this can be done better  
         CLASS(GmresSolver), INTENT(INOUT)         :: this
         REAL(KIND = RP)                           :: RHS(:)
         
         this%RHS = RHS
      END SUBROUTINE  SetRHS
!/////////////////////////////////////////////////////////////////////////   
      SUBROUTINE SetInitialGuess(this,x0)
         CLASS(GmresSolver), INTENT(INOUT)         :: this
         REAL(KIND = RP)                           :: x0(:)
         
         this%x0 = x0
      END SUBROUTINE
!///////////////////////////////////////////////////////////////////////// 
      RECURSIVE SUBROUTINE innerGMRES(this)   !arueda: This subroutine has been made recursive since GMRES can be preconditioned using another GMRES
         CLASS(GmresSolver), INTENT(INOUT)      :: this
         INTEGER                                :: i,j,k, l, ii,kk, m
         REAL(KIND = RP)                        :: tmp1, tmp2
         
         !Compute first krylov vector
         CALL this%MatMult(this%x0,this%V(:,1))
         this%V(:,1) = this%RHS - this%V(:,1)   
         this%g(1) = NORM2(this%V(:,1))
         this%V(:,1) = this%V(:,1) / this%g(1)
         
         this%H = 0.0_RP
         m = this%m

         DO j = 1,m ! Krylov loop
            IF (this%FLEXIBLE) THEN
               CALL this%PCMult(this%V(:,j),this%Z(:,j))
               CALL this%MatMult(this%Z(:,j),this%W)
            ELSE
               CALL this%MatMult(this%V(:,j),this%W)
            ENDIF
            
            DO i = 1,j
               this%H(i,j) = DOT_PRODUCT(this%W,this%V(:,i))
               this%W = this%W - this%H(i,j) * this%V(:,i)
            END DO
            this%H(j+1,j) = NORM2(this%W)
            IF ((ABS(this%H(j+1,j)) .LT. this%tol)) THEN
               this%CONVERGED = .TRUE.
               this%res = this%tol
               this%niter = this%niter + 1                 
               m = j
               EXIT
            ENDIF
            this%V(:,j+1) =  this%W / this%H(j+1,j) 
            DO i = 1, j-1
               tmp1 = this%H(i,j)
               tmp2 = this%H(i+1,j)
               this%H(i,j) = this%cc(i) * tmp1 + this%ss(i) * tmp2
               this%H(i+1,j) = this%cc(i) * tmp2 - this%ss(i) * tmp1
            END DO 
            tmp1 = SQRT(this%H(j,j)*this%H(j,j) + this%H(j+1,j)*this%H(j+1,j) )    
            IF (ABS(tmp1) .LT. 1e-15_RP) THEN
               this%ERROR_CODE = -1
               RETURN
            END IF
            this%cc(j) = this%H(j,j) / tmp1
            this%ss(j) = this%H(j+1,j) / tmp1
            this%g(j+1) = -this%ss(j) * this%g(j)
            this%g(j) = this%cc(j) * this%g(j)
            this%H(j,j) = this%cc(j) * this%H(j,j) + this%ss(j) * this%H(j+1,j) 
            this%res = ABS(this%g(j+1))
            this%niter = this%niter + 1
            IF (this%res .LT. this%tol) THEN
               this%CONVERGED = .TRUE.
               m = j
               EXIT
            ENDIF
            IF (this%niter .GE. this%maxiter) THEN
               m = j
               this%CONVERGED = .FALSE.
               EXIT
            END IF
         END DO ! End of Krylov loop
         
         IF (m > 0) THEN
            this%y(m) = this%g(m) / this%H(m,m)
            DO ii = 1, m-1
               kk = m - ii
               tmp1 = this%g(kk)
               l = kk+1
               DO WHILE (l .LE. m)
                  tmp1 = tmp1 - this%H(kk,l) * this%y(l)
                  this%y(kk) = tmp1 / this%H(kk,kk)
                  l = l+1
               END DO 
            END DO
         END IF ! m > 0
         
         IF (this%FLEXIBLE) THEN
            this%x = this%x0 + MATMUL(this%Z(:,1:m),this%y(1:m))
         ELSE
            this%x = this%x0 + MATMUL(this%V(:,1:m),this%y(1:m))
         ENDIF
       END SUBROUTINE innerGMRES
 !/////////////////////////////////////////////////////////////////////////      
      RECURSIVE SUBROUTINE SolveGMRES(this)    !arueda: This subroutine has been made recursive since GMRES can be preconditioned using another GMRES
         CLASS(GmresSolver), INTENT(INOUT)      :: this
         
         INTEGER                                :: i
         INTEGER                                :: MAX_OUTER_ITER = 5000 ! this is a "security" limit, should never be reached....
       
         this%niter = 0
         this%CONVERGED = .FALSE.
         this%x0 = 0._RP
         DO i = 1, MAX_OUTER_ITER
            CALL innerGMRES(this)
            IF (this%ERROR_CODE .NE. 0) THEN ! Some day, and that day may never come, this will be useful
               PRINT*, 'ERROR IN GMRES, ERROR CODE: ', this%ERROR_CODE
               CALL this%Destruct()
               STOP
            ENDIF
            IF(this%CONVERGED .OR. (this%niter .GE. this%maxiter)) THEN
               RETURN
            ELSE
               this%x0 = this%x
            ENDIF
          END DO
          WRITE(*,*) " ** Reached max outer iter limit in GMRES ** "
      END SUBROUTINE SolveGMRES
      
!/////////////////////////////////////////////////////////////////////////      
!      SUBROUTINE MPRINT(A)
!         REAL(KIND = RP), INTENT(IN)         :: A(:,:)
!         CHARACTER(LEN=40)                   :: fmt
!         INTEGER                             :: i, mn(2)
!         
!         mn = SHAPE(A)
!         WRITE(fmt,*) "(1x,",mn(2),"f7.3)"
!         WRITE(*,fmt) (A(i,:), i = 1,mn(1))
!      END SUBROUTINE
   
END MODULE GMRESClass
!/////////////////////////////////////////////////////////////////////////
!/////////////////////////////////////////////////////////////////////////
