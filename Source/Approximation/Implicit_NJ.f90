!
!////////////////////////////////////////////////////////////////////////
!
!      Implicit_NJ.f90
!      Created: 2017-03-XX XX:XX:XX -XXXX 
!      By:  Carlos Redondo (module for 2D) 
!           Andrés Rueda   (3D implementation and changes) 
!      Implicit module using BDF1 and numerical Jacobian computed using colorings
!
!////////////////////////////////////////////////////////////////////////
MODULE Implicit_NJ

   USE SMConstants,                 ONLY: RP                  
   USE DGSEMClass,                  ONLY: DGSem, ComputeTimeDerivative
   USE ElementClass,                ONLY: Element, allocateElementStorage    !arueda: No DGSolutionStorage implemented in nslite3d... Using whole element definitions
   USE PhysicsStorage,              ONLY: N_EQN, N_GRAD_EQN, flowIsNavierStokes
   USE Jacobian,                    ONLY: Neighbour, Look_for_neighbour                    !arueda: not ready (Neighbor is actually defined in DGSEMClass.)
   USE ColorsClass,                 ONLY: Colors
   USE PetscSolverClass,            ONLY: PetscKspLinearSolver
   USE CSR_Matrices
   
   IMPLICIT NONE
   
   TYPE(Neighbour),ALLOCATABLE, SAVE                     :: nbr(:)
   TYPE(Colors), SAVE                                    :: ecolors
   
   REAL(KIND = RP)                                       :: time         ! Time at the beginning of each inner time step

   PRIVATE                          
   PUBLIC                           :: TakeBDFStep_NJ

   CONTAINS
   !/////////////////////////////////////////////////////////////////////////////////////////////////
   SUBROUTINE TakeBDFStep_NJ (sem, t , dt , maxResidual)

      IMPLICIT NONE
      TYPE(DGSem),                  INTENT(INOUT)           :: sem
      REAL(KIND=RP),                INTENT(IN)              :: t
      REAL(KIND=RP),                INTENT(IN)              :: dt
      REAL(KIND=RP)                                         :: maxResidual(N_EQN)
      
      !--------------------------------------------------------
      REAL(KIND=RP)                                         :: localMaxResidual(N_EQN)
      INTEGER                                               :: id, eq
      !--------------------------------------------------------
      
        
      TYPE(PetscKspLinearSolver)                            :: linsolver
      INTEGER                                               :: cli, clf, clrate
      INTEGER                                               :: k, nelm, DimPrb, newtonit
      INTEGER                                               :: ninner = 1
      INTEGER, DIMENSION(:), ALLOCATABLE                    :: Nx, Ny, Nz
      LOGICAL                                               :: isfirst = .TRUE., computeA = .TRUE.   !arueda: variable INNERIT deleted (not used)
      REAL(KIND=RP)                                         :: ConvRate, ctime
      CHARACTER(LEN=15)                                     :: filename
      REAL(KIND=RP)                                         :: coeff 
      REAL(KIND=RP)                                         :: inner_dt
      REAL(KIND=RP), DIMENSION(:), ALLOCATABLE              :: U_n                                   !Solution at the beginning of time step (even for inner time steps)
      LOGICAL                                               :: PRINT_NEWTON_INFO = .TRUE., PRINT_A = .FALSE., CONVERGED
      
      
      !NEWTON PARAMETERS, TODO: this should be defined in a better place
      REAL(KIND=RP)  :: minrate = 0.5_RP              ! If newton convergence rate lower this value,  newton loop stops and inner_dt is reduced  
      REAL(KIND=RP)  :: maxrate = 1.7_RP              ! If newton loop convergence rate passes this value, inner_dt is increased
      REAL(KIND=RP)  :: NEWTON_TOLERANCE = 1e-6_RP   ! newton iter tolerance relative to the first iter norm 
      INTEGER        :: MAX_NEWTON_ITER = 30          ! If newton iter reachs this limit, this iteration is marked as  not converged 
      INTEGER        :: LIM_NEWTON_ITER = 12          ! If Newton converges but this limit is reached, jacobian matrix will be recomputed
      LOGICAL        :: FORCE_CALCJAC = .TRUE.        ! If true, jacobian matrix is calculated every newton iteration and previous paremters are ignored
      
      SAVE isfirst, computeA, ninner
      SAVE u_N, DimPrb, nelm, linsolver

      IF (isfirst) THEN                         
         isfirst = .FALSE.
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
         CALL linsolver%construct(DimPrb)             !Constructs Petsc linear solver 
      ENDIF
      
      inner_dt = dt            ! first inner_dt is the outer step dt     !arueda: out of isfirst for solving time-accurate problems
      time = t
     
!~       IF (computeA) THEN
!~          CALL ComputeJacMatrix(sem, time+inner_dt, ecolors, nbr, linsolver, nelm, .TRUE., .FALSE.) !arueda: hereIam   ! J     !arueda: removed "ExternalState, ExternalGradients" (not listed arguments in routine definition)
!~          CALL linsolver%SetOperatorDt(inner_dt)                                         ! A = J - I/dt
!~!          computeA = .FALSE.     !arueda: this line is commented out so that the Jacobian is computed in every time step
!~       ENDIF

      CALL StoreSolution( sem, nelm, U_n )             !stores sem%dgS(elmnt)%Q in Vector u_N
      
      
      DO                                                 
         CALL NewtonSolve(sem, time+inner_dt, inner_dt, ecolors, nbr, linsolver, nelm, U_n, MAX_NEWTON_ITER, NEWTON_TOLERANCE, &
                          PRINT_NEWTON_INFO, minrate,ConvRate, newtonit,CONVERGED)                              
                          !old arguments ( sem, linsolver, inner_dt, nelm, U_n, MAX_NEWTON_ITER, NEWTON_TOLERANCE, &
                                          !PRINT_NEWTON_INFO, minrate, ConvRate, newtonit, CONVERGED)                     !arueda: removed "ExternalState, ExternalGradients"!
         
         IF (CONVERGED) THEN
!~             IF (newtonit .GT. LIM_NEWTON_ITER) THEN   !Recomputes jacobian Matrix if convergence rate is poor          !arueda: outcommented cause it's only for steady-state... Add flag and proof if it's better this way
!~                IF (PRINT_NEWTON_INFO) THEN
!~                   WRITE(*,*) "Convergence rate is poor,  recomputing jacobian matrix..."
!~                ENDIF
!~                CALL ComputeJacMatrix(sem, t, ecolors, nbr, linsolver, nelm, .TRUE.,.FALSE.,ExternalState, ExternalGradients)
!~                CALL linsolver%SetOperatorDt(inner_dt)                                         
!~             ENDIF           
            
            !**************************************************
            ! This is new
            time = time + inner_dt
            CALL StoreSolution( sem, nelm, U_n )
            
            IF (ABS((time)-(t+dt)) < 10 * EPSILON(1._RP)) THEN       ! If outer t+dt is reached, the time integration is done
               
               EXIT                                            
            ENDIF
            
            !Increase Comp_Dt if good convergence in previous step ¿¿IF (newtonit .LE. LIM_NEWTON_ITER) THEN??
            IF (ConvRate > maxrate) THEN
               inner_dt = inner_dt * 2.0_RP
!~                CALL linsolver%SetOperatorDt(inner_dt)    ! Sets the operator with the new dt
               IF (PRINT_NEWTON_INFO) WRITE(*,*) "Increasing  dt  = ", inner_dt
            ENDIF
            
            ! Adjust inner_dt to prevent "inner_t" be greater than outer Dt 
            IF ( time+inner_dt > t + dt) THEN  ! Adjusts inner dt to achieve exact outer Dt in the last substep
               inner_dt = t + dt - time
!~                CALL linsolver%SetOperatorDt(inner_dt)    ! Sets the operator with the new dt
            ENDIF
            
            !
            !*************************************************
            
            
!~             IF (time2dt - inner_dt < 0.99*inner_dt) THEN         !if time2dt is 0 (smaller than inner dt) bdftimestep is done
!~                time = time + inner_dt
!~                EXIT
!~             ELSE
!~                time2dt = time2dt - inner_dt
!~                CALL StoreSolution( sem, nelm, U_n )      !stores sem%dgS(elmnt)%Q in Vector u_N
!~                print*, newtonit, LIM_NEWTON_ITER
!~             ENDIF
         ELSE  ! Reduce dt
            inner_dt = inner_dt / 2._RP
!~             CALL linsolver%SetOperatorDt(inner_dt)    ! Sets the operator with the new dt
            CALL RestoreSolution( sem, nelm, U_n )    ! restores Q to begin a new newton iteration       
             IF (PRINT_NEWTON_INFO) WRITE(*,*) "Newton loop did not converge, trying a smaller dt = ", inner_dt
         END IF
      
      END DO
 
      IF (PRINT_NEWTON_INFO) WRITE(*,'(A10,f5.2)') "ConvRate: ", ConvRate
      WRITE(*,'(A11,1p,e8.2,A11,1p,e8.2)') "Outer DT = ", dt, "  Inner DT = ", inner_dt
      
      !IF (ConvRate <1.3_RP) computeA = .TRUE.
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
   
!/////////////////////////////////////////////////////////////////////////////////////////////////
   SUBROUTINE NewtonSolve(sem, t, dt, ecolors, nbr, linsolver, nelm, U_n, MAX_NEWTON_ITER, NEWTON_TOLERANCE, &
                          INFO, minrate,ConvRate, niter,converged)   ! arueda: ", ExternalState, ExternalGradients" deleted since they are not used
!     
!     ----------------------
!     Input-Output arguments
!     ----------------------
!      
      TYPE(DGSem),                  INTENT(INOUT)           :: sem
      REAL(KIND=RP),                INTENT(IN)              :: t
      REAL(KIND=RP),                INTENT(IN)              :: dt              !arueda: why inner_dt and dt??? (inner_dt is not used: deleted!!)
      TYPE(Colors),                 INTENT(IN)              :: ecolors
      TYPE(Neighbour),DIMENSION(:), INTENT(IN)              :: nbr
      TYPE(PetscKspLinearSolver),   INTENT(INOUT)           :: linsolver   !Linear operator is calculate outside this subroutine
      INTEGER,                      INTENT(IN)              :: nelm
      REAL(KIND=RP), DIMENSION(0:), INTENT(IN)              :: u_N
      INTEGER,                      INTENT(IN)              :: MAX_NEWTON_ITER
      REAL(KIND=RP),                INTENT(IN)              :: NEWTON_TOLERANCE
      LOGICAL,                      INTENT(IN)              :: INFO
      REAL(KIND=RP),                INTENT(IN)              :: minrate
      REAL(KIND=RP),                INTENT(OUT)             :: ConvRate
      INTEGER,                      INTENT(OUT)             :: niter
      LOGICAL,                      INTENT(OUT)             :: converged    
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
      TYPE(csrMat)                                       :: Matcsr
        
      
      norm = 1.0_RP
      norm_old = -1.0_RP  !Must be initialized to -1 to avoid bad things in the first newton iter
      newtonit = 1
      ConvRate = 1.0_RP
   
      IF (INFO) THEN
         PRINT*, "Newton it     Newton abs_err   Newton rel_err   LinSolverErr   # ksp iter   Iter wall time (s)"
      END IF
      
      DO newtonit = 1, MAX_NEWTON_ITER                                 !NEWTON LOOP
         CALL ComputeJacMatrix(sem, t, ecolors, nbr, linsolver, nelm, .FALSE., .FALSE.) !arueda: hereIam   ! J     !arueda: removed "ExternalState, ExternalGradients" (not listed arguments in routine definition)
!~         CALL linsolver%SaveMat &
!~                   ('/home/andresrueda/Dropbox/PhD/03_Initial_Codes/3D/Implicit/nslite3d/Tests/Euler/NumJac/MatMatlab.dat')
         CALL linsolver%SetOperatorDt(dt)    
         
!~          CALL linsolver%GetCSRMatrix(Matcsr)
!~          print*, 'after getting csr'
!~          CALL Matcsr%Visualize('Mat.dat')
!~          STOP
         
         CALL SYSTEM_CLOCK(COUNT=cli, COUNT_RATE=clrate)         
         CALL ComputeRHS(sem, dt, u_N, nelm, linsolver )               ! Computes b (RHS) and stores it into PetscVec
         CALL linsolver%solve(norm*1e-3_RP, 500)                       ! Solve (J-I/dt)·x = (Q_r- U_n)/dt - Qdot_r
         IF (.NOT. linsolver%converged) THEN                           ! If linsolver did not converge, return converged=false
            converged = .FALSE.
            RETURN
         ENDIF
         CALL UpdateNewtonSol(sem, nelm, linsolver)                    ! Q_r+1 = Q_r + x
         CALL SYSTEM_CLOCK(COUNT=clf)
         norm = linsolver%xnorm

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
            WRITE(*, "(I8,1p,E18.3,E18.3,E15.3,I10,f18.5)"),newtonit, norm, norm/norm1, linsolver%residual,&
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
         
!~          IF (STRICT_NEWTON) THEN  ! If true, jacobian matrix is recomputed for each newton iteration
!~             CALL ComputeJacMatrix(sem, t, ecolors, nbr, linsolver, nelm, .FALSE.,.FALSE.)  !arueda: second argument: t... not 0._RP!!
!~             CALL linsolver%SetOperatorDt(dt)
!~          ENDIF
         
      ENDDO
   
   END SUBROUTINE

!/////////////////////////////////////////////////////////////////////////////////////////////////
   SUBROUTINE StoreSolution( sem, nelm, u_store )

      TYPE(DGSem),   INTENT(IN)                    :: sem
      INTEGER,       INTENT(IN)                    :: nelm
      REAL(KIND=RP), INTENT(INOUT)                 :: u_store(0:)

      INTEGER                                      :: counter, i, j, k, l, elm, Nx, Ny, Nz

      !sem%dgS(elmnt)%Q --> U_n vector
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

      !sem%dgS(elmnt)%Q --> U_n vector
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
!/////////////////////////////////////////////////////////////////////////////////////////////////     
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
!/////////////////////////////////////////////////////////////////////////////////////////////////
   SUBROUTINE ComputeRHS(sem, dt, U_n, nelm, linsolver )

      TYPE(DGSem),                INTENT(IN)       :: sem
      REAL(KIND=RP),              INTENT(IN)       :: dt
      REAL(KIND=RP),              INTENT(IN)       :: U_n(0:)
      INTEGER,                    INTENT(IN)       :: nelm
      TYPE(PetscKspLinearSolver), INTENT (INOUT)   :: linsolver

      INTEGER                                      :: Nx, Ny, Nz, l, i, j, k, elmnt, counter   
      REAL(KIND=RP)                                :: value

!     Right-Hand side BDF1 (stored into the PETSc vector b):
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
      TYPE(PetscKspLinearSolver),      INTENT(INOUT)    :: linsolver

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
!~                CALL ConservativeToPrimitive2P(sem%dgS(elm)%Q(1:N_EQN,i,j), sem%dgS(elm)%Qp(1:N_PRIM,i,j))       ! arueda: outcommented for being of multiphase flow solver     
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE UpdateNewtonSol
!////////////////////////////////////////////////////////////////////////////////////////////
   SUBROUTINE ComputeJacMatrix(sem, t, ecolors, nbr, linsolver, nelm, PINFO ,PRINT_JAC)
!
!     ---------------
!     Input arguments
!     ---------------
!
      TYPE(DGSem),                INTENT(INOUT), TARGET  :: sem
      REAL(KIND=RP),              INTENT(IN)             :: t
      TYPE(Colors),               INTENT(IN)             :: ecolors
      TYPE(Neighbour),            INTENT(IN)             :: nbr(:)
      TYPE(PetscKspLinearSolver), INTENT(INOUT)          :: linsolver
      INTEGER,                    INTENT(IN)             :: nelm     
      LOGICAL,                    OPTIONAL               :: PINFO
      LOGICAL,                    OPTIONAL               :: PRINT_JAC
!
!     ------------------
!     Internal variables
!     ------------------
!
      REAL(KIND=RP)                                      :: coeffBDF, coeff
      INTEGER                                            :: i, j, kk 
      TYPE(Element), ALLOCATABLE, SAVE                   :: dgs_clean(:)                           !arueda: not only defining storage, but all element definitions
      INTEGER                                            :: ielm, felm, ndofelm, nnz
      INTEGER                                            :: Nx, Ny, Nz
      INTEGER                                            :: thiscolor, thiselmidx, thiselm, thisdof, elmnbr
      INTEGER                                            :: icol, cli, clf, clrate
      INTEGER, DIMENSION(4)                              :: ijkl                                   ! Indexes to locate certain degree of freedom i,j,k...l:equation number
      INTEGER, ALLOCATABLE, DIMENSION(:), SAVE           :: irow_0, irow
!~      REAL(KIND=RP), POINTER, DIMENSION(:)               :: pbuffer                              ! Outcommented, new definition below
      REAL(KIND=RP), ALLOCATABLE, DIMENSION(:)           :: pbuffer                                ! arueda: changed to allocatable variable because of i,j,k,l distribution in nslite3d
      REAL(KIND=RP)                                      :: ctime, maxQ
      REAL(KIND=RP), SAVE                                :: jaceps=1e-8_RP, eps
      CHARACTER(LEN=15)                                  :: filename 
            
      CALL SYSTEM_CLOCK(cli,clrate)
      
      !This suposes same polinomial order in all elements  TODO:Con la nueva estructura extender esto
      Nx = sem%mesh%elements(1)%N 
      Ny = sem%mesh%elements(1)%N
      Nz = sem%mesh%elements(1)%N
      ndofelm = N_EQN * (Nx+1) * (Ny+1) * (Nz+1)
      
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
         maxQ=0.0_RP
         DO i = 1, nelm
            CALL allocateElementStorage( dgs_clean(i), Nx, N_EQN, N_GRAD_EQN, flowIsNavierStokes ) !arueda flowIsNavierStokes defined in Physics (Only using Nx since 3D code doesn't support something else)
            maxQ=MAX(maxQ,MAXVAL(ABS(sem%mesh%elements(i)%Q(:,:,:,:))))
         ENDDO
         eps=1e-8_RP*maxQ
!         write(*,*) maxQ
!         read(*,*)
      ENDIF
    
      nnz = ndofelm * 7                 !max nonzero values expected in a row in the jacobian matrix. arueda: "for 2D GL: 5*ndofelm, for 3D GL: 7*ndofelm, for LGL depends on row!!! (a node on the interface has more cols than an interior node)"

      CALL linsolver%PreallocateA(nnz)
      CALL linsolver%ResetA                                                 

      coeffBDF = 1.0_RP
!      eps = 1e-8_RP

!      CALL ComputeTimeDerivative( sem )
      
      CALL ComputeTimeDerivative( sem, t )
      
      dgs_clean = sem%mesh%elements
      DO thiscolor = 1 , ecolors%ncolors
         ielm = ecolors%bounds(thiscolor)             
         felm = ecolors%bounds(thiscolor+1)
         DO thisdof = 1, ndofelm                      ! computes one column for each dof within an elment
            ijkl = local2ijk(thisdof,N_EQN,Nx,Ny,Nz)
            DO thiselmidx = ielm, felm-1              !perturbs a dof in all elements within current color
               thiselm = ecolors%elmnts(thiselmidx)
               sem%mesh%elements(thiselm)%Q(ijkl(1),ijkl(2),ijkl(3),ijkl(4)) = &
                                                   sem%mesh%elements(thiselm)%Q(ijkl(1),ijkl(2),ijkl(3),ijkl(4)) + eps 
            ENDDO
!            CALL ComputeTimeDerivative( sem )
            CALL ComputeTimeDerivative( sem, t )            
            DO thiselmidx = ielm, felm-1
               thiselm = ecolors%elmnts(thiselmidx)
               DO i = 1,7 !hardcoded: 7 neighbours (only valid for hexahedrals in conforming meshes)
                  elmnbr = nbr(thiselm)%elmnt(i)                 
                  IF (elmnbr .NE. 0) THEN                                                       
                     sem%mesh%elements(elmnbr)%QDot = (sem%mesh%elements(elmnbr)%QDot - dgs_clean(elmnbr)%QDot) / eps                      
!~                     pbuffer(1:ndofelm) => sem%mesh%elements(elmnbr)%QDot                     !maps Qdot array into a 1D pointer 
                     CALL GetElemQdot(sem%mesh%elements(elmnbr),pbuffer)
                     irow = irow_0 + ndofelm * (elmnbr - 1)                         !generates the row indices vector
                     WHERE (ABS(pbuffer(1:ndofelm)) .LT. jaceps) irow = -1          !MatSetvalues will ignore entries with irow=-1
                     icol = (thiselm - 1) * ndofelm  + thisdof - 1
                     CALL linsolver%SetAColumn(ndofelm, irow, icol, pbuffer )
                  ENDIF
               ENDDO
            END DO           
            DO i = 1,nelm                                      !Cleans dgs
               sem%mesh%elements(i)%Q = dgs_clean(i)%Q           
            END DO                                                
         ENDDO
      ENDDO
            
      
      CALL linsolver%AssemblyA                                 ! Petsc Matrix A needs to be assembled before being used
      linsolver%Ashift = 0.0_RP                                ! Shift must be set to zero when a new jac matrix is calculated
      
      CALL SYSTEM_CLOCK(clf)
      ctime = (clf - cli) / REAL(clrate)
      IF (PRESENT(PINFO)) THEN
         IF (PINFO) PRINT*, "Implicit operator computing and assembly time: ", ctime, "s"
      ENDIF
      IF (PRESENT(PRINT_JAC)) THEN
         IF(PRINT_JAC) THEN
            WRITE(filename,"(A2,f6.4,A4)") "A_",t,".dat"
            CALL linsolver%SaveMat(filename)
            CALL WriteJacInfo(N_EQN,nelm,linsolver%dimprb,ctime,'Jac.info') 
         ENDIF
      ENDIF            
   END SUBROUTINE ComputeJacMatrix

!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!  
!  Returns the local index relative to an element from the local coordinates: i(lagrange node x), j(lagrange node y), 
!  k(lagrange node z), l(equation number)
!  M,N are the polinomial orders in x, y directions, N_EQN is the number of equations
   
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
!  Returns the coordinates relative to an elements i(lagrange node x), j(lagrange node y), k(variable) from the local index   
!  M,N are the polinomial orders in x, y directions, N_EQN is the number of equations
   
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
   SUBROUTINE WriteEigenFiles(Mat,sem,FileName)
      IMPLICIT NONE
!
!     -----------------------------------------------------------
!     Writes files for performing eigenvalue analysis using TAUev
!     -----------------------------------------------------------
!
      TYPE(csrMat)      :: Mat      !< Jacobian matrix
      TYPE(DGSem)       :: sem      !< DG class with mesh inside
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

END MODULE Implicit_NJ
