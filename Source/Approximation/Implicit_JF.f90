!
!////////////////////////////////////////////////////////////////////////
!
!      Implicit_JF.f90
!      Created: 2017-03-XX XX:XX:XX -XXXX 
!      By:  Carlos Redondo (module for 2D) 
!           Andrés Rueda   (3D implementation and changes) -- to be listed here
!      Implicit module using BDF1 and Jacobian Free Newton-Krylov
!
!////////////////////////////////////////////////////////////////////////
MODULE Implicit_JF
   USE SMConstants,                 ONLY: RP                  
   USE DGSEMClass,                  ONLY: DGSem, ComputeTimeDerivative
   USE PhysicsStorage,              ONLY: N_EQN, N_GRAD_EQN
   USE JFNKClass,                   ONLY: JFNK_Integrator
   
   IMPLICIT NONE
   
   INTEGER                                            :: nelm, DimPrb
   INTEGER, DIMENSION(:), ALLOCATABLE                 :: Nx, Ny, Nz
   REAL(KIND = RP), DIMENSION(:), ALLOCATABLE         :: U_n
   TYPE(JFNK_Integrator)                              :: integrator
   TYPE(DGSem), POINTER                               :: p_sem
   
   PRIVATE
   PUBLIC      :: TakeBDFStep_JF
 
!========
 CONTAINS
!========
   
   SUBROUTINE TakeBDFStep_JF(sem, t , dt , maxResidual)
      TYPE(DGSem),TARGET,  INTENT(INOUT)           :: sem
      REAL(KIND=RP),       INTENT(IN)              :: t
      REAL(KIND=RP),       INTENT(IN)              :: dt
      !--------------------------------------------------------
      LOGICAL                                      :: isfirst = .TRUE.
      INTEGER                                      :: k
      
      REAL(KIND=RP) :: localMaxResidual(N_EQN), maxResidual(N_EQN)
      INTEGER       :: id, eq
      !--------------------------------------------------------
      
      IF (isfirst) THEN                                      
         isfirst = .FALSE.
         nelm = SIZE(sem%mesh%elements)
         ALLOCATE(Nx(nelm))
         ALLOCATE(Ny(nelm))
         ALLOCATE(Nz(nelm))
         DO k = 1, nelm
            Nx(k) = sem%mesh%elements(k)%N       ! arueda: the routines were originally developed for a code that allows different polynomial orders in different directions. Notation conserved just for the sake of generality (future improvement -?)
            Ny(k) = sem%mesh%elements(k)%N
            Nz(k) = sem%mesh%elements(k)%N
         END DO
         DimPrb = 0
         DO k = 1, nelm
            DimPrb = DimPrb + N_EQN*(Nx(k)+1)*(Ny(k)+1)*(Nz(k)+1)   ! arueda: Currently, it doesn't allow for p-adaptation (an adaptation flag should be added and this computation moved outside of isfirst)
         END DO
         ALLOCATE(U_n(1:Dimprb))
         !CALL integrator%Construct(DimPrb)                  !Constructs jfnk solver
         CALL integrator%Construct(DimPrb, pc = 'GMRES')    !Constructs jfnk solver using preconditioning (hardcoded: GMRES preco)
         p_sem => sem                                       !Makes sem visible by other routines whitin the module
      ENDIF
      
      CALL integrator%SetTime(t)                            ! arueda: sets time at the beginning of iteration
      CALL integrator%SetDt(dt)                             ! Sets Dt in the time integrator     !arueda: sets it always that a time step is "taken"
      
      CALL GetSemQ(p_sem,U_n)                               ! Stores sem%Q into U_n    ! TODO: esto no es necesario hacerlo siempre                      
      CALL integrator%SetUn(U_n)                            ! Sets U_n  in the time integrator
      CALL integrator%SetTimeDerivative(DGTimeDerivative)   ! Sets DGTime derivative as Time derivative function in the integrtor
      CALL integrator%Integrate                             ! Performs the time step implicit integrator
      CALL integrator%GetUnp1(U_n)                          ! Gets U_n+1 and stores it in U n
      CALL SetSemQ(p_sem,U_n)                               ! Stores U_n (with the updated U_n+1) into sem%Q  ! TODO: esto tampoco es necesario hacerlo siempre 
         
!
!     ----------------
!     Compute residual
!     ----------------
!
      maxResidual = 0.0_RP
      DO id = 1, SIZE( p_sem % mesh % elements )
         DO eq = 1 , N_EQN
            localMaxResidual(eq) = MAXVAL(ABS(sem % mesh % elements(id) % QDot(:,:,:,eq)))
            maxResidual(eq) = MAX(maxResidual(eq),localMaxResidual(eq))
         END DO
      END DO
         
   END SUBROUTINE TakeBDFStep_JF
!
!//////////////////////////////////////////////////////////////////////////////////////// 
!
   FUNCTION DGTimeDerivative(u,time) RESULT(F)
      REAL(KIND = RP), INTENT(IN)                  :: u(:)
      REAL(KIND = RP)                              :: F(size(u)) 
      REAL(KIND = RP)                              :: time
  
      CALL SetSemQ(p_sem, u)
      CALL ComputeTimeDerivative(p_sem,time)
      CALL GetSemQdot(p_sem, F)
    
   END FUNCTION DGTimeDerivative
!////////////////////////////////////////////////////////////////////////////////////////       
      SUBROUTINE SetSemQ(sem,Q)
         TYPE(DGSem),         INTENT(INOUT)           :: sem 
         REAL(KIND = RP),     INTENT(IN)              :: Q(:)   
         
         INTEGER                                       :: Nx, Ny, Nz, l, i, j, k, counter, elm
         
         counter = 1
         DO elm = 1, nelm
            Nx = sem%mesh%elements(elm)%N ! arueda: the routines were originally developed for a code that allows different polynomial orders in different directions. Notation conserved just for the sake of generality (future improvement -?)
            Ny = sem%mesh%elements(elm)%N
            Nz = sem%mesh%elements(elm)%N
            DO k = 0, Nz
               DO j = 0, Ny
                  DO i = 0, Nx
                     DO l = 1,N_EQN
                         sem%mesh%elements(elm)%Q(i, j, k, l) = Q(counter) 
                        counter =  counter + 1
                     END DO
                  END DO
               END DO
            END DO
         END DO 
         
      END SUBROUTINE SetSemQ
!////////////////////////////////////////////////////////////////////////////////////////       
      SUBROUTINE GetSemQ(sem,Q)
         TYPE(DGSem),         INTENT(INOUT)            :: sem
         REAL(KIND = RP),     INTENT(OUT)              :: Q(:)
         
         INTEGER                                       :: Nx, Ny, Nz, l, i, j, k, counter, elm
         
         counter = 1
         DO elm = 1, nelm
            Nx = sem%mesh%elements(elm)%N ! arueda: the routines were originally developed for a code that allows different polynomial orders in different directions. Notation conserved just for the sake of generality (future improvement -?)
            Ny = sem%mesh%elements(elm)%N
            Nz = sem%mesh%elements(elm)%N
            DO k = 0, Nz
               DO j = 0, Ny
                   DO i = 0, Nx
                     DO l = 1,N_EQN
                        Q(counter)  = sem%mesh%elements(elm)%Q(i, j, k, l)
                        counter =  counter + 1
                     END DO
                   END DO
               END DO
            END DO
         END DO
         
      END SUBROUTINE GetSemQ
 !////////////////////////////////////////////////////////////////////////////////////////      
      SUBROUTINE GetSemQdot(sem,Qdot)
         TYPE(DGSem),         INTENT(INOUT)            :: sem
         REAL(KIND = RP),     INTENT(OUT)              :: Qdot(:)
         
         INTEGER                                       :: Nx, Ny, Nz, l, i, j, k, counter, elm
         
         counter = 1
         DO elm = 1, nelm
            Nx = sem%mesh%elements(elm)%N ! arueda: the routines were originally developed for a code that allows different polynomial orders in different directions. Notation conserved just for the sake of generality (future improvement -?)
            Ny = sem%mesh%elements(elm)%N
            Nz = sem%mesh%elements(elm)%N
            DO k = 0, Nz
               DO j = 0, Ny
                  DO i = 0, Nx
                     DO l = 1,N_EQN
                        Qdot(counter)  = sem%mesh%elements(elm)%Qdot(i, j, k, l) 
                        counter =  counter + 1
                     END DO
                  END DO
               END DO
            END DO
         END DO
         
      END SUBROUTINE GetSemQdot

END MODULE Implicit_JF
