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

   USE SMConstants                  
   USE DGSEMClass,                  ONLY: DGSem, ComputeTimeDerivative
   USE ElementClass,                ONLY: Element, allocateElementStorage    !arueda: No DGSolutionStorage implemented in nslite3d... Using whole element definitions
   USE PhysicsStorage,              ONLY: N_EQN, N_GRAD_EQN, flowIsNavierStokes
   USE Jacobian,                    ONLY: Neighbour, Look_for_neighbour           
   USE ColorsClass,                 ONLY: Colors
   USE LinearSolverClass
!~   USE PetscSolverClass
!~   USE MKLPardisoSolverClass
   
   USE CSR_Matrices
   USE FTValueDictionaryClass
   
   IMPLICIT NONE
   
   TYPE(Neighbour),ALLOCATABLE, SAVE                     :: nbr(:)
   TYPE(Colors), SAVE                                    :: ecolors
   
   REAL(KIND = RP)                                       :: time         ! Time at the beginning of each inner(!) time step

   PRIVATE                          
   PUBLIC                           :: TakeBDFStep_NJ

   CONTAINS
   !/////////////////////////////////////////////////////////////////////////////////////////////////
   SUBROUTINE TakeBDFStep_NJ (sem, t , dt , maxResidual, controlVariables)

      IMPLICIT NONE
      TYPE(DGSem),                  INTENT(INOUT)           :: sem                  !<>DGSem class with solution storage 
      REAL(KIND=RP),                INTENT(IN)              :: t                    !< Time at the beginning of time step
      REAL(KIND=RP),                INTENT(IN)              :: dt                   !< Initial (outer) time step (can internally, the subroutine can use a smaller one depending on convergence)
      REAL(KIND=RP)                                         :: maxResidual(N_EQN)   !> Residuals obtained after time integration process
      TYPE(FTValueDictionary),      INTENT(IN)              :: controlVariables     !< Input file variables
      !--------------------------------------------------------
      !Variables used for residual computation
      REAL(KIND=RP)                                         :: localMaxResidual(N_EQN)
      INTEGER                                               :: id, eq
      !--------------------------------------------------------
      CHARACTER(len=LINE_LENGTH)                            :: LinearSolver
        
      CLASS(GenericLinSolver_t), POINTER                    :: linsolver           ! Linear solver (as an abstract type, it must be declared as CLASS)
      INTEGER                                               :: cli, clf, clrate
      INTEGER                                               :: k, nelm, DimPrb, newtonit
      INTEGER                                               :: ninner = 1
      INTEGER, DIMENSION(:), ALLOCATABLE                    :: Nx, Ny, Nz
      LOGICAL                                               :: isfirst = .TRUE., computeA = .TRUE.
      REAL(KIND=RP)                                         :: ConvRate
      REAL(KIND=RP)                                         :: inner_dt
      REAL(KIND=RP), DIMENSION(:), ALLOCATABLE              :: U_n                                   !Solution at the beginning of time step (even for inner time steps)
      LOGICAL                                               :: PRINT_NEWTON_INFO = .TRUE., CONVERGED
      LOGICAL                                               :: JacByConv                               ! 
      ! Not used variables?
      CHARACTER(LEN=15)                                     :: filename
      REAL(KIND=RP)                                         :: ctime
      
      !NEWTON PARAMETERS, TODO: this should be defined in a better place
      REAL(KIND=RP)  :: minrate = 0.5_RP              ! If newton convergence rate lower this value,  newton loop stops and inner_dt is reduced  
      REAL(KIND=RP)  :: maxrate = 1.7_RP              ! If newton loop convergence rate passes this value, inner_dt is increased
      REAL(KIND=RP)  :: NEWTON_TOLERANCE = 1e-6_RP    ! newton iter tolerance relative to the first iter norm 
      INTEGER        :: MAX_NEWTON_ITER = 30          ! If newton iter reachs this limit, this iteration is marked as  not converged 
      INTEGER        :: LIM_NEWTON_ITER = 12          ! If Newton converges but this limit is reached, jacobian matrix will be recomputed
      
      
      
      SAVE isfirst, computeA, ninner, JacByConv
      SAVE u_N, DimPrb, nelm, linsolver
      

      IF (isfirst) THEN           
         isfirst = .FALSE.
         
         !Which linear solver?
         LinearSolver = controlVariables % StringValueForKey("linear solver",LINE_LENGTH)
         SELECT CASE (LinearSolver)
            CASE('petsc')
               ALLOCATE (PetscKspLinearSolver_t :: linsolver)
            CASE('pardiso')
               ALLOCATE (MKLPardisoSolver_t :: linsolver)
            CASE DEFAULT
               print*, "Keyword 'linear solver' missing... Using PETSc as default"
               ALLOCATE (PetscKspLinearSolver_t :: linsolver)
         END SELECT
         
         nelm = SIZE(sem%mesh%elements)
         ALLOCATE(Nx(nelm))
         ALLOCATE(Ny(nelm))
         ALLOCATE(Nz(nelm))
         DO k = 1, nelm
            Nx(k) = sem%mesh%elements(k)%N ! arueda: the routines were originally developed for a code that allows different polynomial orders in different directions. Notation conserved just for the sake of generality (future improvement -?)
            Ny(k) = sem%mesh%elements(k)%N
            Nz(k) = sem%mesh%elements(k)%N 
         END DO
         DimPrb = 0
         DO k = 1, nelm
            DimPrb = DimPrb + N_EQN*(Nx(k)+1)*(Ny(k)+1)*(Nz(k)+1)
         END DO

         ALLOCATE(nbr(nelm))
         CALL Look_for_neighbour(nbr, sem)    
         ALLOCATE(U_n(0:Dimprb-1))
         CALL ecolors%construct(nbr,flowIsNavierStokes)       
         !CALL ecolors%info
         CALL linsolver%construct(DimPrb)             !Constructs linear solver 
         JacByConv = controlVariables % LogicalValueForKey("jacobian by convergence")
      ENDIF
      
      inner_dt = dt            ! first inner_dt is the outer step dt     !arueda: out of isfirst for solving time-accurate problems
      time = t
      
      !**************************
      ! If the Jacobian must only be computed sometimes
       IF (JacByConv) THEN
         IF (computeA) THEN
            CALL ComputeNumJac(sem, time+inner_dt, ecolors, nbr, linsolver, nelm, .TRUE.,.FALSE.) 
            CALL linsolver%SetOperatorDt(inner_dt)
            computeA = .FALSE. 
         ELSE
            CALL linsolver%ReSetOperatorDt(inner_dt)
         END IF
       ENDIF
      ! 
      !**************************
      
      CALL StoreSolution( sem, nelm, U_n )             !stores sem%dgS(elmnt)%Q in Vector u_N
      
      
      DO                                                 
         CALL NewtonSolve(sem, time+inner_dt, inner_dt, ecolors, nbr, linsolver, nelm, U_n, MAX_NEWTON_ITER, NEWTON_TOLERANCE, &
                          PRINT_NEWTON_INFO, minrate,JacByConv,ConvRate, newtonit,CONVERGED)
         
         IF (CONVERGED) THEN
            time = time + inner_dt
            CALL StoreSolution( sem, nelm, U_n )
            
            !*************************************************
            !
            IF (JacByConv .AND. newtonit .GT. LIM_NEWTON_ITER) THEN   !Recomputes jacobian Matrix if convergence rate is poor          !arueda: outcommented cause it's only for steady-state... Add flag and proof if it's better this way
               IF (PRINT_NEWTON_INFO) THEN
                  WRITE(*,*) "Convergence rate is poor,  recomputing jacobian matrix..."
               ENDIF
               CALL ComputeNumJac(sem, time+inner_dt, ecolors, nbr, linsolver, nelm, .TRUE., .FALSE.)
               CALL linsolver%SetOperatorDt(inner_dt)                                         
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
            
            CALL RestoreSolution( sem, nelm, U_n )       ! restores Q to begin a new newton iteration       
             IF (PRINT_NEWTON_INFO) WRITE(*,*) "Newton loop did not converge, trying a smaller dt = ", inner_dt
         END IF
      
      END DO
 
      IF (PRINT_NEWTON_INFO) WRITE(*,'(A10,f5.2)') "ConvRate: ", ConvRate
      WRITE(*,'(A11,1p,e8.2,A11,1p,e8.2)') "Outer DT = ", dt, "  Inner DT = ", inner_dt
      
      !**************************
      ! for computing sometimes
      IF (JacByConv .AND. ConvRate <0.65_RP .AND. newtonit .LT. LIM_NEWTON_ITER) computeA = .TRUE.  !last condition cause' if that's true, the Jacobian was just computed
      ! for computing sometimes
      !**************************
!
!     ----------------
!     Compute residual
!     ----------------
!
      maxResidual = 0.0_RP
      DO id = 1, SIZE( sem % mesh % elements )
         DO eq = 1 , N_EQN
            localMaxResidual(eq) = MAXVAL(ABS(sem % mesh % elements(id) % QDot(:,:,:,eq)))
            maxResidual(eq) = MAX(maxResidual(eq),localMaxResidual(eq))
         END DO
      END DO
      
   END SUBROUTINE TakeBDFStep_NJ
!
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE NewtonSolve(sem, t, dt, ecolors, nbr, linsolver, nelm, U_n, MAX_NEWTON_ITER, NEWTON_TOLERANCE, &
                          INFO, minrate,JacByConv,ConvRate, niter,CONVERGED)
!     
!     ----------------------
!     Input-Output arguments
!     ----------------------
!      
      TYPE(DGSem),                  INTENT(INOUT)           :: sem
      REAL(KIND=RP),                INTENT(IN)              :: t
      REAL(KIND=RP),                INTENT(IN)              :: dt              !< Inner dt
      TYPE(Colors),                 INTENT(IN)              :: ecolors
      TYPE(Neighbour),DIMENSION(:), INTENT(IN)              :: nbr
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
!     
!     ------------------
!     Internal variables
!     ------------------
! 
      INTEGER                                               :: cli, clf, clrate           
      INTEGER                                               :: newtonit
      REAL(KIND=RP)                                         :: norm, norm_old, rel_tol, norm1
      LOGICAL                                               :: STRICT_NEWTON = .TRUE.
      
      ! Temp
!~      TYPE(csrMat)                                       :: Matcsr
        
      
      norm = 1.0_RP
      norm_old = -1.0_RP  !Must be initialized to -1 to avoid bad things in the first newton iter
      ConvRate = 1.0_RP
   
      IF (INFO) THEN
         PRINT*, "Newton it     Newton abs_err   Newton rel_err   LinSolverErr   # ksp iter   Iter wall time (s)"
      END IF
      
      DO newtonit = 1, MAX_NEWTON_ITER                                 !NEWTON LOOP
         
         IF (.NOT. JacByConv) THEN !If Jacobian must be computed always
            CALL ComputeNumJac(sem, t, ecolors, nbr, linsolver, nelm, .TRUE., .FALSE.) 
            CALL linsolver%SetOperatorDt(dt)
         ELSE
            CALL ComputeTimeDerivative( sem, t )
         END IF
!~          CALL linsolver%GetCSRMatrix(Matcsr)
!~          print*, 'after getting csr'
!~          CALL WriteEigenFiles(Matcsr,sem,'SmallMat')
!~          STOP
         
         CALL SYSTEM_CLOCK(COUNT=cli, COUNT_RATE=clrate)         
         CALL ComputeRHS(sem, dt, u_N, nelm, linsolver )               ! Computes b (RHS) and stores it into linsolver
         CALL linsolver%solve(tol=norm*1.e-12_RP, maxiter=1000)        ! Solve (J-I/dt)·x = (Q_r- U_n)/dt - Qdot_r
         IF (.NOT. linsolver%converged) THEN                           ! If linsolver did not converge, return converged=false
            converged = .FALSE.
            RETURN
         ENDIF
         CALL UpdateNewtonSol(sem, nelm, linsolver)                    ! Q_r+1 = Q_r + x
         CALL SYSTEM_CLOCK(COUNT=clf)
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
            WRITE(*, "(I8,1p,E18.3,E18.3,E15.3,I10,f18.5)"),newtonit, norm, norm/norm1, linsolver%Getrnorm(),&
                                                      linsolver%niter,(clf-cli)/real(clrate)   
         ENDIF
         
         IF (ConvRate < minrate .OR. newtonit == MAX_NEWTON_ITER) THEN
            converged = .FALSE.
            RETURN
         ENDIF
        
         IF (norm < MAX(rel_tol,NEWTON_TOLERANCE)) THEN
            converged = .TRUE.
            RETURN
         ENDIF
         
      ENDDO
   
   END SUBROUTINE

!/////////////////////////////////////////////////////////////////////////////////////////////////
   SUBROUTINE StoreSolution( sem, nelm, u_store )

      TYPE(DGSem),   INTENT(IN)                    :: sem
      INTEGER,       INTENT(IN)                    :: nelm
      REAL(KIND=RP), INTENT(INOUT)                 :: u_store(0:)

      INTEGER                                      :: counter, i, j, k, l, elm, Nx, Ny, Nz

      !U_n vector <-- sem%dgS(elmnt)%Q
      counter = 0
      DO elm = 1, nelm
         Nx = sem%mesh%elements(elm)%N ! arueda: the routines were originally developed for a code that allows different polynomial orders in different directions. Notation conserved just for the sake of generality (future improvement -?)
         Ny = sem%mesh%elements(elm)%N
         Nz = sem%mesh%elements(elm)%N
         DO k = 0, Nz
            DO j = 0, Ny
               DO i = 0, Nx
                  DO l = 1, N_EQN
                     u_store(counter) = sem%mesh%elements(elm)%Q(i,j,k,l)
                     counter =  counter + 1
                  END DO 
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE StoreSolution
!/////////////////////////////////////////////////////////////////////////////////////////////////
   SUBROUTINE RestoreSolution( sem, nelm, u_store )

      TYPE(DGSem),   INTENT(INOUT)                 :: sem
      INTEGER,       INTENT(IN)                    :: nelm
      REAL(KIND=RP), INTENT(IN)                    :: u_store(0:)

      INTEGER                                      :: counter, i, j, k, l, elm, Nx, Ny, Nz

      !sem%dgS(elmnt)%Q <-- U_n vector
      counter = 0
      DO elm = 1, nelm
         Nx = sem%mesh%elements(elm)%N ! arueda: the routines were originally developed for a code that allows different polynomial orders in different directions. Notation conserved just for the sake of generality (future improvement -?)
         Ny = sem%mesh%elements(elm)%N
         Nz = sem%mesh%elements(elm)%N
         DO k = 0, Nz
            DO j = 0, Ny
               DO i = 0, Nx
                  DO l = 1, N_EQN
                     sem%mesh%elements(elm)%Q(i,j,k,l) = u_store(counter)
                     counter =  counter + 1
                  END DO 
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE RestoreSolution
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

!     Right-Hand side BDF1 (stored into the linsolver vector b):
!     b = [(Q(n+1) - Q(n))/dt] - Qdot  == [(U_r - U_n ) / dt] - F(U_r)
!     U_r is stored in sem%dgs%Q

!      CALL ComputeTimeDerivative( sem ) ! computes Qdot
      
      counter = 0
      DO elmnt = 1, nelm
         Nx = sem%mesh%elements(elmnt)%N ! arueda: the routines were originally developed for a code that allows different polynomial orders in different directions. Notation conserved just for the sake of generality (future improvement -?)
         Ny = sem%mesh%elements(elmnt)%N
         Nz = sem%mesh%elements(elmnt)%N
         DO k = 0, Nz
            DO j = 0, Ny
               DO i = 0, Nx
                  DO l = 1,N_EQN
                     value = (sem%mesh%elements(elmnt)%Q(i,j,k,l) - U_n(counter))/dt - sem%mesh%elements(elmnt)%QDot(i,j,k,l)
                     CALL linsolver%SetBValue(counter, value)
                     counter =  counter + 1
                  END DO
               END DO
            END DO
         END DO
      END DO
      CALL linsolver%AssemblyB     ! b must be assembled before using them
   END SUBROUTINE ComputeRHS
!//////////////////////////////////////////////////////////////////////////////////////////////
   SUBROUTINE UpdateNewtonSol(sem, nelm, linsolver)

      TYPE(DGSem),                     INTENT(INOUT)    :: sem
      INTEGER,                         INTENT(IN)       :: nelm
      CLASS(GenericLinSolver_t),       INTENT(INOUT)    :: linsolver

      REAL(KIND=RP)                                     :: value
      INTEGER                                           :: Nx, Ny, Nz, l, i, j, k, counter, elm

      counter = 0
      DO elm = 1, nelm
         Nx = sem%mesh%elements(elm)%N ! arueda: the routines were originally developed for a code that allows different polynomial orders in different directions. Notation conserved just for the sake of generality (future improvement -?)
         Ny = sem%mesh%elements(elm)%N
         Nz = sem%mesh%elements(elm)%N
         DO k = 0, Nz
            DO j = 0, Ny
               DO i = 0, Nx
                  DO l = 1,N_EQN
                     CALL linsolver%GetXValue(counter,value)
                     sem%mesh%elements(elm)%Q(i,j,k,l) = sem%mesh%elements(elm)%Q(i,j,k,l) + value
                     counter =  counter + 1
                  END DO    
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE UpdateNewtonSol
!
!////////////////////////////////////////////////////////////////////////////////////////////
!  
   SUBROUTINE WriteEigenFiles(Mat,sem,FileName)
      IMPLICIT NONE
!
!     -----------------------------------------------------------
!     Writes files for performing eigenvalue analysis using TAUev
!     -----------------------------------------------------------
!
      TYPE(csrMat)      :: Mat      !< Jacobian matrix
      TYPE(DGSem)       :: sem      !< DGSem class containing mesh
      CHARACTER(len=*)  :: FileName !< ...
!     -----------------------------------------------------------
      INTEGER           :: fd
!     -----------------------------------------------------------
      
      ! .frm file
      OPEN(newunit=fd, file=TRIM(FileName)//'.frm', action='WRITE')
         WRITE(fd,*)
         WRITE(fd,*) SIZE(Mat % Values), SIZE(Mat % Rows)-1, 1, N_EQN, 1
         WRITE(fd,*) sem % mesh % elements(1) % N, SIZE(sem % mesh % elements)
      CLOSE (fd)
      
      ! .amg file
      CALL Mat % Visualize(TRIM(FileName)//'.amg',FirstRow=.FALSE.)
      
      ! .coo file
      CALL sem % mesh % WriteCoordFile(TRIM(FileName)//'.coo')
      
      
   END SUBROUTINE WriteEigenFiles
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!
!     Routines for computing the numerical Jacobian
!
!
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE ComputeNumJac(sem, t, ecolors, nbr, linsolver, nelm, PINFO ,PRINT_JAC)
!
!     ---------------
!     Input arguments
!     ---------------
!
      TYPE(DGSem),                INTENT(INOUT), TARGET  :: sem
      REAL(KIND=RP),              INTENT(IN)             :: t
      TYPE(Colors),               INTENT(IN)             :: ecolors
      TYPE(Neighbour),            INTENT(IN)             :: nbr(:)
      CLASS(GenericLinSolver_t),  INTENT(INOUT)          :: linsolver
      INTEGER,                    INTENT(IN)             :: nelm     
      LOGICAL,                    OPTIONAL               :: PINFO
      LOGICAL,                    OPTIONAL               :: PRINT_JAC
!
!     ------------------
!     Internal variables
!     ------------------
!
      INTEGER                                            :: i, j                                   ! generic counters 
      INTEGER                                            :: thiscolor, thiselmidx, thiselm         ! specific counters
      INTEGER                                            :: thisdof, elmnbr, nbrnbr                ! specific counters
      INTEGER, ALLOCATABLE, DIMENSION(:), SAVE           :: used                                   ! array containing index of elements whose contributions to Jacobian has already been considered
      INTEGER                                            :: usedctr                                ! counter to fill positions of used
      
      TYPE(Element), ALLOCATABLE, SAVE                   :: dgs_clean(:)                           ! arueda: not only defining storage, but all element definitions (change it?)
      INTEGER                                            :: ielm, felm, ndofelm, nnz
      INTEGER                                            :: Nx, Ny, Nz                             ! Polynomial order
      
      INTEGER                                            :: icol, cli, clf, clrate
      INTEGER, DIMENSION(4)                              :: ijkl                                   ! Indexes to locate certain degree of freedom i,j,k...l:equation number
      INTEGER, ALLOCATABLE, DIMENSION(:), SAVE           :: irow_0, irow
!~      REAL(KIND=RP), POINTER, DIMENSION(:)               :: pbuffer                              ! Outcommented, new definition below (arueda: previous was more efficient.. change after defining storage)
      REAL(KIND=RP), ALLOCATABLE, DIMENSION(:), SAVE     :: pbuffer                                ! arueda: changed to allocatable variable because of i,j,k,l distribution in nslite3d
      REAL(KIND=RP)                                      :: ctime, maxQ
      REAL(KIND=RP), SAVE                                :: jaceps=1e-8_RP, eps
      CHARACTER(LEN=15)                                  :: filename 
      
            
      CALL SYSTEM_CLOCK(cli,clrate)
      
      !This suposes same polinomial order in all elements  TODO:Con la nueva estructura extender esto
      Nx = sem%mesh%elements(1)%N 
      Ny = sem%mesh%elements(1)%N
      Nz = sem%mesh%elements(1)%N
      ndofelm = N_EQN * (Nx+1) * (Ny+1) * (Nz+1)
!
!     ---------------------------------------------------------------
!     Allocate arrays that will be used throughout all the simulation
!     ---------------------------------------------------------------
!
      IF (.NOT. ALLOCATED(pbuffer)) ALLOCATE(pbuffer(ndofelm))  ! arueda: this can be done here and once only because all elements have same order. TODO: generalize!
      
      IF (.NOT. ALLOCATED(irow)) THEN
         ALLOCATE(irow  (0:ndofelm-1))
         ALLOCATE(irow_0(0:ndofelm-1))
         DO i = 0, ndofelm-1
            irow_0(i) = i
         ENDDO
      ENDIF
      
      IF (.NOT. ALLOCATED(dgs_clean)) THEN
         ALLOCATE(dgs_clean(nelm))
         
         
         DO i = 1, nelm
            CALL allocateElementStorage( dgs_clean(i), Nx, N_EQN, N_GRAD_EQN, flowIsNavierStokes ) !(Only using Nx since 3D code doesn't support something else)
!~             maxQ=MAX(maxQ,MAXVAL(ABS(sem%mesh%elements(i)%Q(:,:,:,:))))
         ENDDO
         
      ENDIF
      
      IF (.NOT. ALLOCATED(used)) THEN
         IF (flowIsNavierStokes) THEN ! .AND. BR1 (only implementation to date)
            ALLOCATE(used(26))   ! 25 neighbors (including itself) and a last entry that will be 0 always (boundary index)
         ELSE
            ALLOCATE(used(8))    ! 8 neighbors (including itself) and a last entry that will be 0 always (boundary index)
         END IF
      END IF
!
!     ---------------------------------------------
!     Set value of eps
!        See:
!           > Knoll, Dana A., and David E. Keyes. "Jacobian-free Newton–Krylov methods: a survey of approaches and applications." Journal of Computational Physics 193.2 (2004): 357-397.
!           > Mettot, Clément, Florent Renac, and Denis Sipp. "Computation of eigenvalue sensitivity to base flow modifications in a discrete framework: Application to open-loop control." Journal of Computational Physics 269 (2014): 234-258.
!     --------------------------------------------
!
      maxQ=0.0_RP
      DO i = 1, nelm
         maxQ=MAX(maxQ,MAXVAL(ABS(sem%mesh%elements(i)%Q(:,:,:,:))))
      ENDDO
      eps=SQRT(EPSILON(eps))*(maxQ+1._RP)
!
!     -------------------------------------------------------------------------
!     Set max number of nonzero values expected in a row of the Jacobian matrix 
!           Assumes Legendre-Gauss quadrature and neglects zero values in each 
!              block (taken into account later when assembling)
!           For Legendre-Gauss-Lobatto much less entries are expected (a node on the
!              interface has more cols than an interior node)
!           IMPORTANT: These numbers assume conforming meshes!
!     -------------------------------------------------------------------------
!
      IF (flowIsNavierStokes) THEN ! .AND. BR1 (only implementation to date)
         nnz = ndofelm * 25
      ELSE
         nnz = ndofelm * 7
      END IF
      
      CALL linsolver%PreallocateA(nnz)
      CALL linsolver%ResetA
      
      CALL ComputeTimeDerivative( sem, t )
!
!     ------------------------------------------
!     Compute numerical Jacobian using colorings
!     ------------------------------------------
!
      dgs_clean = sem%mesh%elements
      DO thiscolor = 1 , ecolors%ncolors
         ielm = ecolors%bounds(thiscolor)             
         felm = ecolors%bounds(thiscolor+1)
         DO thisdof = 1, ndofelm                      ! computes one column for each dof within an elment (all elements must have the same number of DOFs - change?)
            ijkl = local2ijk(thisdof,N_EQN,Nx,Ny,Nz)
            DO thiselmidx = ielm, felm-1              !perturbs a dof in all elements within current color
               thiselm = ecolors%elmnts(thiselmidx)
               sem%mesh%elements(thiselm)%Q(ijkl(1),ijkl(2),ijkl(3),ijkl(4)) = &
                                                   sem%mesh%elements(thiselm)%Q(ijkl(1),ijkl(2),ijkl(3),ijkl(4)) + eps 
            ENDDO
            
            CALL ComputeTimeDerivative( sem, t )            
            
            DO thiselmidx = ielm, felm-1
               thiselm = ecolors%elmnts(thiselmidx)
               ! Redifine used array and counter
               used    = 0
               usedctr = 1
               
               DO i = 1,SIZE(nbr(thiselm)%elmnt)
                  elmnbr = nbr(thiselm)%elmnt(i) 
                  
                  IF (.NOT. ANY(used == elmnbr)) THEN  !(elmnbr .NE. 0)
                     sem%mesh%elements(elmnbr)%QDot = (sem%mesh%elements(elmnbr)%QDot - dgs_clean(elmnbr)%QDot) / eps                      
!~                     pbuffer(1:ndofelm) => sem%mesh%elements(elmnbr)%QDot                     !maps Qdot array into a 1D pointer 
                     CALL GetElemQdot(sem%mesh%elements(elmnbr),pbuffer)
                     irow = irow_0 + ndofelm * (elmnbr - 1)                         !generates the row indices vector
                     WHERE (ABS(pbuffer(1:ndofelm)) .LT. jaceps) irow = -1          !MatSetvalues will ignore entries with irow=-1
                     icol = (thiselm - 1) * ndofelm  + thisdof - 1
                     CALL linsolver%SetAColumn(ndofelm, irow, icol, pbuffer )
                     
                     used(usedctr) = elmnbr
                     usedctr = usedctr + 1
                  ENDIF
                  
                  ! If we are using BR1, we also have to get the contributions of the neighbors of neighbors
                  IF(flowIsNavierStokes) THEN ! .AND. BR1 (only implementation to date)
                     IF (elmnbr .NE. 0) THEN
                        DO j=1, SIZE(nbr(elmnbr)%elmnt)
                           nbrnbr = nbr(elmnbr)%elmnt(j)                          
                           IF (.NOT. ANY(used == nbrnbr)) THEN  !(elmnbr .NE. 0)
                                             
                              sem%mesh%elements(nbrnbr)%QDot = (sem%mesh%elements(nbrnbr)%QDot - dgs_clean(nbrnbr)%QDot) / eps                      
         !~                     pbuffer(1:ndofelm) => sem%mesh%elements(elmnbr)%QDot                     !maps Qdot array into a 1D pointer 
                              CALL GetElemQdot(sem%mesh%elements(nbrnbr),pbuffer)
                              irow = irow_0 + ndofelm * (nbrnbr - 1)                         !generates the row indices vector
                              WHERE (ABS(pbuffer(1:ndofelm)) .LT. jaceps) irow = -1          !MatSetvalues will ignore entries with irow=-1
                              icol = (thiselm - 1) * ndofelm  + thisdof - 1
                              CALL linsolver%SetAColumn(ndofelm, irow, icol, pbuffer )
                              
                              used(usedctr) = nbrnbr
                              usedctr = usedctr + 1                        
                           ENDIF
                        END DO
                     END IF
                  END IF
                  
               ENDDO
            END DO           
            DO thiselmidx = ielm, felm-1                              !Cleans modified Qs
               thiselm = ecolors%elmnts(thiselmidx)
               sem%mesh%elements(thiselm)%Q = dgs_clean(thiselm)%Q           
            END DO                                                
         ENDDO
      ENDDO
      sem%mesh%elements = dgs_clean ! Cleans sem % mesh % elements completely
      
      CALL linsolver%AssemblyA                                 ! Matrix A needs to be assembled before being used in PETSc (at least)
      
      CALL SYSTEM_CLOCK(clf)
      ctime = (clf - cli) / REAL(clrate)
      IF (PRESENT(PINFO)) THEN
         IF (PINFO) PRINT*, "Implicit operator computing and assembly time: ", ctime, "s"
      ENDIF
      IF (PRESENT(PRINT_JAC)) THEN
         IF(PRINT_JAC) THEN
            WRITE(filename,"(A2,f6.4,A4)") "A_",t,".dat"
!~             CALL linsolver%SaveMat(filename)
            CALL WriteJacInfo(N_EQN,nelm,linsolver%dimprb,ctime,'Jac.info') 
         ENDIF
      ENDIF            
   END SUBROUTINE ComputeNumJac
!
!/////////////////////////////////////////////////////////////////////////////////////////////////     
!
   SUBROUTINE WriteJacInfo(NEq,NoE,DimPrb,t,filename)

      INTEGER,          INTENT(IN)     :: NEq,NoE,DimPrb
      REAL(KIND=RP),    INTENT(IN)     :: t
      CHARACTER(LEN=*), OPTIONAL       :: filename
      
      IF (.NOT. PRESENT(filename)) filename = 'Jac.info'
      OPEN(44, FILE=filename)
      WRITE(44,*) " # eqs  # elements      # DOFs    Jac comp time (s)"
      WRITE(44,"(I7,I12,I12,f6.0)") NEq,NoE,DimPrb,t
      CLOSE(44) 
   END SUBROUTINE WriteJacInfo
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!  
!  Returns the local index relative to an element from the local coordinates: i(lagrange node x), j(lagrange node y), 
!  k(lagrange node z), l(equation number)
!  N are the polinomial orders in x, y and z directions, N_EQN is the number of equations
   
   FUNCTION ijk2local(i,j,k,l,N_EQN,Nx,Ny,Nz) RESULT(idx)
      IMPLICIT NONE
      
      INTEGER, INTENT(IN)   :: i, j, k, l, Nx, Ny, Nz, N_EQN
      INTEGER               :: idx
      
      IF (l < 1 .OR. l > N_EQN)  STOP 'error in ijk2local, l has wrong value'
      IF (i < 0 .OR. i > Nx)     STOP 'error in ijk2local, i has wrong value'
      IF (j < 0 .OR. j > Ny)     STOP 'error in ijk2local, j has wrong value'
      IF (k < 0 .OR. k > Nz)     STOP 'error in ijk2local, k has wrong value'
      
      idx = k*(Nx+1)*(Ny+1)*N_EQN + j*(Nx+1)*N_EQN + i*N_EQN + l
   END FUNCTION
   
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  Returns the coordinates relative to an element: i(lagrange node x), j(lagrange node y), k(lagrange node z), l(equation number)
!  from the local index  
!  N are the polinomial orders in x, y and z directions, N_EQN is the number of equations
   
   FUNCTION local2ijk(idx,N_EQN,Nx,Ny,Nz) RESULT (indices)
   
      INTEGER, INTENT(IN)   :: idx, Nx, Ny, Nz, N_EQN
      INTEGER               :: indices(4)
      INTEGER               :: tmp1, tmp2
      
      IF (idx < 1 .OR. idx > (Nx+1)*(Ny+1)*(Nz+1)*N_EQN) STOP 'error in local2ijk, idx has wrong value'
      
      indices(3) = (idx-1) / ((Nx+1)*(Ny+1) * N_EQN)
      tmp1       = MOD((idx-1),((Nx+1)*(Ny+1) * N_EQN) )
      indices(2) = tmp1 / ((Nx+1)*N_EQN)
      tmp2       = MOD(tmp1,((Nx+1)*N_EQN))
      indices(1) = tmp2 / (N_EQN)
      indices(4) = MOD(tmp2, N_EQN) + 1
   END FUNCTION
   
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

!////////////////////////////////////////////////////////////////////////////////////////      
!  Subroutine for extracting Qdot of a single element as a 1 dimensional array
!  arueda: Originally, this process was done by a pbuffer 1D pointer (better performance)... However, copying procedure had to be introduced since 
!  NSLITE3D organizes the element information in a different manner than BNLITE2D...
!   
   SUBROUTINE GetElemQdot(CurrEl,Qdot) !arueda: check ordering of variables in solution vector
      TYPE(Element)                                 :: CurrEl
      REAL(KIND = RP),     INTENT(OUT)              :: Qdot(:)
      
      INTEGER                                       :: Nx, Ny, Nz, l, i, j, k, counter
      
      counter = 1
      
      Nx = CurrEl%N ! arueda: the routines were originally developed for a code that allows different polynomial orders in different directions. Notation conserved just for the sake of generality (future improvement -?)
      Ny = CurrEl%N
      Nz = CurrEl%N
      DO k = 0, Nz
         DO j = 0, Ny
            DO i = 0, Nx
               DO l = 1,N_EQN
                  Qdot(counter)  = CurrEl%Qdot(i, j, k, l) 
                  counter =  counter + 1
               END DO
            END DO
         END DO
      END DO
      
   END SUBROUTINE GetElemQdot
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
END MODULE Implicit_NJ
