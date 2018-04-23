!
!////////////////////////////////////////////////////////////////////////
!
!      Implicit_NJ.f90
!      Created: 2017-04-09 16:30:00 +0100 
!      By:  Carlos Redondo (module for 2D) 
!           Andrés Rueda   (3D implementation and changes) 
!      Implicit module using BDF1 and numerical Jacobian computed using colorings
!
!////////////////////////////////////////////////////////////////////////
MODULE Implicit_NJ
   use AnalyticalJacobian
   USE SMConstants                  
   USE DGSEMClass,                  ONLY: DGSem
   USE ElementClass,                ONLY: Element, allocateElementStorage    !arueda: No DGSolutionStorage implemented in nslite3d... Using whole element definitions
   USE PhysicsStorage
   use HexMeshClass
   USE LinearSolverClass
   USE CSRMatrixClass
   USE FTValueDictionaryClass
   use TimeIntegratorDefinitions
   use MatrixClass
   use DGSEMClass, only: ComputeQDot_FCN
   implicit none
   
   PRIVATE                          
   PUBLIC TakeBDFStep_NJ
   
   real(kind=RP) :: time               ! Time at the beginning of each inner(!) time step
   logical       :: computeA = .TRUE.  ! Compute Jacobian? (only valid if it is meant to be computed according to the convergence)
   
CONTAINS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!   
   SUBROUTINE TakeBDFStep_NJ (sem, t , dt , controlVariables, ComputeTimeDerivative)

      IMPLICIT NONE
      TYPE(DGSem),                  INTENT(INOUT)           :: sem                  !<>DGSem class with solution storage 
      REAL(KIND=RP),                INTENT(IN)              :: t                    !< Time at the beginning of time step
      REAL(KIND=RP),                INTENT(IN)              :: dt                   !< Initial (outer) time step (can internally, the subroutine can use a smaller one depending on convergence)
      TYPE(FTValueDictionary),      INTENT(IN)              :: controlVariables     !< Input file variables
      procedure(ComputeQDot_FCN)                            :: ComputeTimeDerivative
      !--------------------------------------------------------
      CHARACTER(len=LINE_LENGTH)                            :: LinearSolver
        
      CLASS(GenericLinSolver_t), POINTER                    :: linsolver           ! Linear solver (as an abstract type, it must be declared as CLASS)
      INTEGER                                               :: cli, clf, clrate
      INTEGER                                               :: k, nelm, DimPrb, newtonit
      INTEGER                                               :: ninner = 1
      LOGICAL                                               :: isfirst = .TRUE.
      REAL(KIND=RP)                                         :: ConvRate
      REAL(KIND=RP)                                         :: inner_dt
      REAL(KIND=RP), DIMENSION(:), ALLOCATABLE              :: U_n                                   !Solution at the beginning of time step (even for inner time steps)
      LOGICAL                                               :: PRINT_NEWTON_INFO, CONVERGED
      LOGICAL                                               :: JacByConv                               ! 
      LOGICAL                                               :: TimeAccurate = .FALSE., UserNewtonTol = .FALSE.
      ! Not used variables?
      CHARACTER(LEN=15)                                     :: filename
      REAL(KIND=RP)                                         :: ctime
      integer                                               :: nEqnJac, nGradJac
      
      !NEWTON PARAMETERS, TODO: this should be defined in a better place
      REAL(KIND=RP)  :: minrate = 0.5_RP              ! If newton convergence rate lower this value,  newton loop stops and inner_dt is reduced  
      REAL(KIND=RP)  :: maxrate = 1.7_RP              ! If newton loop convergence rate passes this value, inner_dt is increased
      REAL(KIND=RP)  :: NEWTON_TOLERANCE =1.e-6_RP    ! newton iter tolerance relative to the first iter norm !=1.e-6_RP
      INTEGER        :: MAX_NEWTON_ITER = 30          ! If newton iter reachs this limit, this iteration is marked as  not converged 
      INTEGER        :: LIM_NEWTON_ITER = 12          ! If Newton converges but this limit is reached, jacobian matrix will be recomputed
      
      
      SAVE isfirst, ninner, JacByConv, PRINT_NEWTON_INFO
      SAVE u_N, DimPrb, nelm, linsolver, TimeAccurate, UserNewtonTol

#if defined(NAVIERSTOKES)
      nEqnJac = NCONS
      nGradJac = NGRAD
#endif
      
      IF (isfirst) THEN           
         isfirst = .FALSE.
         
         !Which linear solver?
         LinearSolver = controlVariables % StringValueForKey("linear solver",LINE_LENGTH)
         SELECT CASE (LinearSolver)
            CASE('petsc')
               ALLOCATE (PetscKspLinearSolver_t :: linsolver)
            CASE('pardiso')
               ALLOCATE (MKLPardisoSolver_t     :: linsolver)
            CASE('smooth')
               ALLOCATE (IterativeSolver_t      :: linsolver)
            CASE('multigrid')
               ALLOCATE (MultigridSolver_t      :: linsolver)
            CASE DEFAULT
               print*, "Keyword 'linear solver' missing... Using PETSc as default"
               ALLOCATE (PetscKspLinearSolver_t :: linsolver)
         END SELECT
         
         PRINT_NEWTON_INFO = controlVariables % logicalValueForKey("print newton info")
         
         nelm = SIZE(sem%mesh%elements)
         
         DimPrb = sem % NDOF
         
         
         ALLOCATE(U_n(0:Dimprb-1))

         CALL linsolver%construct(DimPrb,controlVariables,sem)             !Constructs linear solver 
         JacByConv = controlVariables % LogicalValueForKey("jacobian by convergence")
         
         IF (controlVariables % StringValueForKey("time integration",LINE_LENGTH) == 'time-accurate') TimeAccurate = .TRUE.
         
         IF (controlVariables % containsKey("newton tolerance")) THEN
            UserNewtonTol = .TRUE.
            NEWTON_TOLERANCE = controlVariables % doublePrecisionValueForKey("newton tolerance")
         END IF
      ENDIF
      
      IF (.NOT. TimeAccurate .AND. .NOT. UserNewtonTol) THEN
         NEWTON_TOLERANCE = sem % MaxResidual* 1e-3_RP
      END IF
      
      inner_dt = dt            ! first inner_dt is the outer step dt     !arueda: out of isfirst for solving time-accurate problems
      time = t
      
      !**************************
      ! If the Jacobian must only be computed sometimes
       IF (JacByConv) THEN
         IF (.not. computeA) THEN
            CALL linsolver%ReSetOperatorDt(inner_dt)
         END IF
       ENDIF
      ! 
      !**************************
      
      CALL sem % GetQ(U_n)      !stores sem%mesh%elements(:)% storage % Q in Vector U_n
      
      DO                                                 
         CALL NewtonSolve(sem, nEqnJac, nGradJac, time+inner_dt, inner_dt, linsolver, nelm, U_n, MAX_NEWTON_ITER, NEWTON_TOLERANCE, &
                          PRINT_NEWTON_INFO, minrate,JacByConv,ConvRate, newtonit,CONVERGED, ComputeTimeDerivative)
         
         IF (CONVERGED) THEN
            time = time + inner_dt
            CALL sem % GetQ(U_n) 
            
            !*************************************************
            !
            IF (JacByConv .AND. newtonit .GT. LIM_NEWTON_ITER) THEN   !Recomputes jacobian Matrix if convergence rate is poor
               IF (PRINT_NEWTON_INFO) THEN
                  WRITE(*,*) "Convergence rate is poor, Jacobian matrix will be computed in next iteration..."
               ENDIF
               computeA = .TRUE.                                        
            ENDIF
            !
            !*************************************************
            
            
            !*************************************************
            !
            IF (ABS((time)-(t+dt)) < 10 * EPSILON(1._RP)) THEN       ! If outer t+dt is reached, the time integration is done
               EXIT                                            
            ENDIF
            
            !Increase Comp_Dt if good convergence in previous step ¿¿IF (newtonit .LE. LIM_NEWTON_ITER) THEN??
            IF (ConvRate > maxrate) THEN
               inner_dt = inner_dt * 2.0_RP
               IF (JacByConv)  CALL linsolver%ReSetOperatorDt(inner_dt)    ! Resets the operator with the new dt
               
               IF (PRINT_NEWTON_INFO) WRITE(*,*) "Increasing  dt  = ", inner_dt
            ENDIF
            
            ! Adjust inner_dt to prevent "inner_t" be greater than outer Dt 
            IF ( time+inner_dt > t + dt) THEN  ! Adjusts inner dt to achieve exact outer Dt in the last substep
               inner_dt = t + dt - time
               IF (JacByConv)  CALL linsolver%ReSetOperatorDt(inner_dt)    ! Resets the operator with the new dt
               
               IF (PRINT_NEWTON_INFO) WRITE(*,*) "Adjusting dt = ", inner_dt
            ENDIF
            
            !
            !*************************************************
            
         ELSE  ! Reduce dt
            inner_dt = inner_dt / 2._RP
            IF (JacByConv)  CALL linsolver%ReSetOperatorDt(inner_dt)    ! Resets the operator with the new dt
            
            CALL sem % SetQ(U_n)          ! restores Q in sem to begin a new newton iteration       
            IF (PRINT_NEWTON_INFO) WRITE(*,*) "Newton loop did not converge, trying a smaller dt = ", inner_dt
         END IF
      
      END DO
 
      IF (PRINT_NEWTON_INFO) WRITE(*,'(A10,f5.2)') "ConvRate: ", ConvRate
      !WRITE(*,'(A11,1p,e8.2,A11,1p,e8.2)') "Outer DT = ", dt, "  Inner DT = ", inner_dt
      
      !**************************
      ! for computing sometimes
      IF (JacByConv .AND. ConvRate <0.65_RP .AND. newtonit .LT. LIM_NEWTON_ITER) THEN
         computeA = .TRUE.
         !computeA = linsolver % ComputeANextStep()
      END IF
      ! for computing sometimes
      !**************************
      
!~       IF (MAXVAL(maxResidual) > sem % maxResidual) computeA = .TRUE.
      
   END SUBROUTINE TakeBDFStep_NJ
!
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE NewtonSolve(sem, nEqn, nGradEqn, t, dt, linsolver, nelm, U_n, MAX_NEWTON_ITER, NEWTON_TOLERANCE, &
                          INFO, minrate,JacByConv,ConvRate, niter,CONVERGED, ComputeTimeDerivative)
!     
!     ----------------------
!     Input-Output arguments
!     ----------------------
!      
      TYPE(DGSem),                  INTENT(INOUT)           :: sem
      integer,                      intent(in)              :: nEqn
      integer,                      intent(in)              :: nGradEqn
      REAL(KIND=RP),                INTENT(IN)              :: t
      REAL(KIND=RP),                INTENT(IN)              :: dt              !< Inner dt
      CLASS(GenericLinSolver_t),    INTENT(INOUT)           :: linsolver       !Linear operator is calculate outside this subroutine
      INTEGER,                      INTENT(IN)              :: nelm
      REAL(KIND=RP), DIMENSION(0:), INTENT(IN)              :: u_N
      INTEGER,                      INTENT(IN)              :: MAX_NEWTON_ITER
      REAL(KIND=RP),                INTENT(IN)              :: NEWTON_TOLERANCE
      LOGICAL,                      INTENT(IN)              :: INFO
      REAL(KIND=RP),                INTENT(IN)              :: minrate
      LOGICAL,                      INTENT(IN)              :: JacByConv         !< Must the Jacobian be computed for bad convergence? if .false., the Jacobian is computed at the beginning of every newton it
      REAL(KIND=RP),                INTENT(OUT)             :: ConvRate
      INTEGER,                      INTENT(OUT)             :: niter
      LOGICAL,                      INTENT(OUT)             :: CONVERGED   
      procedure(ComputeQDot_FCN)                            :: ComputeTimeDerivative
!~       TYPE(csrMat_t) :: Acsr        !   CSR matrix 
!     
!     ------------------
!     Internal variables
!     ------------------
! 
      INTEGER(8)                                               :: cli, clf, clrate           
      INTEGER                                               :: newtonit
      REAL(KIND=RP)                                         :: norm, norm_old, rel_tol, norm1
!~       LOGICAL, SAVE :: isfirst = .TRUE.
!~       SAVE norm1
      
!~       IF (isfirst) THEN
         norm = 1.0_RP
!~          isfirst = .FALSE.
!~       ELSE
!~          norm = norm1
!~       END IF
      norm_old = -1.0_RP  !Must be initialized to -1 to avoid bad things in the first newton iter
      ConvRate = 1.0_RP
   
      IF (INFO) THEN
         PRINT*, "Newton it     Newton abs_err   Newton rel_err   LinSolverErr   # ksp iter   Iter wall time (s)"
      END IF
      
      CALL SYSTEM_CLOCK(COUNT_RATE=clrate)
      
      DO newtonit = 1, MAX_NEWTON_ITER                                 !NEWTON LOOP
         if (.not. JacByConv) computeA = .TRUE.
         
         CALL ComputeTimeDerivative( sem % mesh, sem % particles, t, sem % BCFunctions)
         CALL ComputeRHS(sem, dt, U_n, nelm, linsolver )               ! Computes b (RHS) and stores it into linsolver
         
         CALL SYSTEM_CLOCK(COUNT=cli)
         CALL linsolver%solve(nEqn, nGradEqn, tol=norm*1.e-3_RP, maxiter=500, time= t, dt=dt, &
                              ComputeTimeDerivative = ComputeTimeDerivative, computeA = computeA)        ! Solve (J-I/dt)·x = (Q_r- U_n)/dt - Qdot_r
         CALL SYSTEM_CLOCK(COUNT=clf)
         IF (.NOT. linsolver%converged) THEN                           ! If linsolver did not converge, return converged=false
            converged = .FALSE.
            RETURN
         ENDIF
         CALL UpdateNewtonSol(sem, nelm, linsolver)                    ! Q_r+1 = Q_r + x
         
         norm = linsolver%Getxnorm('infinity')

         IF (norm_old .NE. -1.0_RP) THEN
            ConvRate = ConvRate + (LOG10(norm_old/norm)-ConvRate)/newtonit 
         ENDIF
         norm_old = norm
         niter = newtonit
         IF (newtonit == 1) THEN
            norm1 = norm
            rel_tol = norm1 * NEWTON_TOLERANCE
         ENDIF
         IF (INFO) THEN
            WRITE(*, "(I8,1p,E18.3,E18.3,E15.3,I10,F18.5)")newtonit, norm, norm/norm1, linsolver%Getrnorm(),&
                                                      linsolver%niter,0.1_RP*(clf-cli)/real(clrate,RP)  !!!! I have NO IDEA why I have to multiply by 0.1!!!
         ENDIF
         
         IF (ConvRate < minrate .OR. newtonit == MAX_NEWTON_ITER .OR. ISNAN(norm)) THEN
            IF (INFO) print*, 'ConvRate: ', ConvRate
            converged = .FALSE.
            RETURN
         ENDIF
        
         IF (norm < MAX(rel_tol,NEWTON_TOLERANCE)) THEN
            converged = .TRUE.
            RETURN
         ENDIF
         
      ENDDO
   
   END SUBROUTINE
!  
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE ComputeRHS(sem, dt, U_n, nelm, linsolver )

      TYPE(DGSem),                INTENT(IN)       :: sem
      REAL(KIND=RP),              INTENT(IN)       :: dt
      REAL(KIND=RP),              INTENT(IN)       :: U_n(0:)
      INTEGER,                    INTENT(IN)       :: nelm
      CLASS(GenericLinSolver_t),  INTENT (INOUT)   :: linsolver

      INTEGER                                      :: Nx, Ny, Nz, l, i, j, k, elmnt, counter   
      REAL(KIND=RP)                                :: value
#if defined(NAVIERSTOKES)

!     Right-Hand side BDF1 (stored into the linsolver vector b):
!     b = [(Q(n+1) - Q(n))/dt] - Qdot  == [(U_r - U_n ) / dt] - F(U_r)
!     U_r is stored in sem%dgs% storage % Q

!      CALL ComputeTimeDerivative( sem ) ! computes Qdot
      
      counter = 0
      DO elmnt = 1, nelm
         Nx = sem%mesh%elements(elmnt)%Nxyz(1)
         Ny = sem%mesh%elements(elmnt)%Nxyz(2)
         Nz = sem%mesh%elements(elmnt)%Nxyz(3)
         DO k = 0, Nz
            DO j = 0, Ny
               DO i = 0, Nx
                  DO l = 1,NCONS
                     value = (sem%mesh%elements(elmnt)% storage % Q(l,i,j,k) - U_n(counter))/dt - &
                              sem%mesh%elements(elmnt)% storage % QDot(l,i,j,k)
                     CALL linsolver%SetBValue(counter, value)
                     counter =  counter + 1
                  END DO
               END DO
            END DO
         END DO
      END DO
      CALL linsolver%AssemblyB     ! b must be assembled before using
#endif
   END SUBROUTINE ComputeRHS
!//////////////////////////////////////////////////////////////////////////////////////////////
   SUBROUTINE UpdateNewtonSol(sem, nelm, linsolver)

      TYPE(DGSem),                     INTENT(INOUT)    :: sem
      INTEGER,                         INTENT(IN)       :: nelm
      CLASS(GenericLinSolver_t),       INTENT(INOUT)    :: linsolver

      REAL(KIND=RP)                                     :: value
      INTEGER                                           :: Nx, Ny, Nz, l, i, j, k, counter, elm
#if defined(NAVIERSTOKES)

      counter = 0
      DO elm = 1, nelm
         Nx = sem%mesh%elements(elm)%Nxyz(1)
         Ny = sem%mesh%elements(elm)%Nxyz(2)
         Nz = sem%mesh%elements(elm)%Nxyz(3)
         DO k = 0, Nz
            DO j = 0, Ny
               DO i = 0, Nx
                  DO l = 1, NCONS
                     CALL linsolver%GetXValue(counter,value)
                     sem%mesh%elements(elm)% storage % Q(l,i,j,k) = sem%mesh%elements(elm)% storage % Q(l,i,j,k) + value
                     counter =  counter + 1
                  END DO    
               END DO
            END DO
         END DO
      END DO

#endif
   END SUBROUTINE UpdateNewtonSol
!
!////////////////////////////////////////////////////////////////////////////////////////////
!  
   SUBROUTINE WriteEigenFiles(Mat,sem,nEqn,FileName)
      IMPLICIT NONE
!
!     -----------------------------------------------------------
!     Writes files for performing eigenvalue analysis using TAUev
!        This only works for isotropic order meshes.........................TODO: Change that
!     -----------------------------------------------------------
!
      TYPE(csrMat_t)      :: Mat      !< Jacobian matrix
      TYPE(DGSem)         :: sem      !< DGSem class containing mesh
      integer, intent(in) :: nEqn
      CHARACTER(len=*)    :: FileName !< ...
!     -----------------------------------------------------------
      INTEGER           :: fd
!     -----------------------------------------------------------
#if defined(NAVIERSTOKES)
      
      ! .frm file
      OPEN(newunit=fd, file=TRIM(FileName)//'.frm', action='WRITE')
         WRITE(fd,*)
         WRITE(fd,*) SIZE(Mat % Values), SIZE(Mat % Rows)-1, 1, NCONS, 1
         WRITE(fd,*) sem % mesh % elements(1) % Nxyz(1), SIZE(sem % mesh % elements)
      CLOSE (fd)
      
      ! .amg file
      CALL Mat % Visualize(TRIM(FileName)//'.amg',FirstRow=.FALSE.)
      
      ! .coo file
      CALL sem % mesh % WriteCoordFile(nEqn,TRIM(FileName)//'.coo')
      
#endif
      
   END SUBROUTINE WriteEigenFiles
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
END MODULE Implicit_NJ
