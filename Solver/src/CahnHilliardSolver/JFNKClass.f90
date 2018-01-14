!
!//////////////////////////////////////////////////////
!
!   @File:    JFNKClass.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Sun Jan 14 13:23:05 2018
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
!
!////////////////////////////////////////////////////////////////////////
!
!      JFNKClass.f90
!      Created: 2017-03-XX XX:XX:XX -XXXX 
!      By:  Carlos Redondo (module for 2D) 
!           Andrés Rueda   (3D implementation and changes) -- to be listed here
!      Jacobian-Free Newton Krylov method using GMRES
!           Implemented using BDF1 and variable time step   !arueda: change this to allow BDF of higher order 
!
!////////////////////////////////////////////////////////////////////////
MODULE JFNKClass
   USE SMConstants,                 ONLY: RP
   USE GMRESClass,                  ONLY: GmresSolver
   IMPLICIT NONE
      
   !Class public variables
   TYPE JFNK_Integrator
      INTEGER                                         :: Dimprb
      REAL(KIND = RP)                                 :: Dt
      REAL(KIND = RP)                                 :: last_Dt
      REAL(KIND = RP)                                 :: inner_DT                       ! inner_Dt
      REAL(KIND = RP), ALLOCATABLE                    :: Un(:)                          ! Solution in last step (for inner iterations): Updated every successful Newton solver call
      REAL(KIND = RP)                                 :: norm_x0                        ! Norm of solution at the beginning of inner loop
      REAL(KIND = RP)                                 :: norm_x_abs             
      REAL(KIND = RP)                                 :: norm_x_rel
      REAL(KIND = RP)                                 :: norm_G                         ! Norm of F(Un) at the beginning of inner loop
      LOGICAL                                         :: F_SET = .FALSE.
      LOGICAL                                         :: CONVERGED = .FALSE.
      INTEGER                                         :: MAX_NEWTON_ITER = 50           !old: 15
      LOGICAL                                         :: MAX_ITER_REACHED = .FALSE.     ! !arueda: allow parameters to be read from input file (and set default values if not in input)
      INTEGER                                         :: niter = 0
      INTEGER                                         :: n_linsolver_iter
      INTEGER                                         :: nsteps
      INTEGER                                         :: substep_red_ratio = 4
      REAL(KIND = RP)                                 :: substep_inc_ratio = 1.5_RP
      REAL(KIND = RP)                                 :: reltol = 1.e-12_RP
      REAL(KIND = RP)                                 :: abstol = 1.e-12_RP
      LOGICAL                                         :: PRINT_NEWTON_INFO = .TRUE.
      TYPE(GmresSolver)                               :: linsolver
      LOGICAL                                         :: PC_GMRES = .FALSE.
      
     
      CONTAINS
         PROCEDURE                                    :: Construct => Construct_JFNK_Integrator
         PROCEDURE                                    :: Desctruct => Destruct_JFNK_Integrator      !arueda: where does "this" get destructed?
         PROCEDURE                                    :: Integrate => JFNK_Integrate
         PROCEDURE                                    :: SetUn
         PROCEDURE                                    :: GetUnp1
         PROCEDURE                                    :: SetTime                                    !arueda: New procedure to set initial time for step.
         PROCEDURE                                    :: SetDt
         PROCEDURE                                    :: SetTimeDerivative
   END TYPE JFNK_Integrator
   

   !Module Interfaces
   ABSTRACT INTERFACE
      FUNCTION TimeDerivative(u,tt) RESULT(F)
         IMPORT RP
         REAL(KIND = RP), INTENT(IN)  :: u(:)
         REAL(KIND = RP)              :: F(size(u))
         REAL(KIND = RP)              :: tt
      END FUNCTION
   END INTERFACE
   
   !Module private variables                      
   PROCEDURE(TimeDerivative), POINTER              :: p_F => NULL()
   REAL(KIND = RP), ALLOCATABLE                    :: F_Ur(:)
   REAL(KIND = RP), ALLOCATABLE                    :: Ur(:)                         ! arueda: Solution updated at the end of each Newton iteration
   REAL(KIND = RP), ALLOCATABLE                    :: Un_sub(:)
   REAL(KIND = RP), ALLOCATABLE                    :: RHS(:)                        !
   REAL(KIND = RP)                                 :: eps =1e-8_RP
   REAL(KIND = RP)                                 :: gmres_tol
   REAL(KIND = RP)                                 :: inner_t
   REAL(KIND = RP)                                 :: inner_Dt
   REAL(KIND = RP)                                 :: Comp_Dt
   REAL(KIND = RP)                                 :: time                          ! Time at the beginning of each inner time step
   INTEGER                                         :: n_preco_iter
   
   TYPE(GmresSolver)                               :: PCsolver


   PRIVATE
   PUBLIC         :: JFNK_Integrator

   CONTAINS
      SUBROUTINE Construct_JFNK_Integrator(this,dimprb,pc)
         CLASS(JFNK_Integrator), INTENT(INOUT)     :: this
         INTEGER, INTENT(IN)                       :: dimprb
         CHARACTER(LEN=*), OPTIONAL                :: pc
         !Used module variables: F_Ur, Ur, RHS
         
         this%dimprb = dimprb
         ALLOCATE(this%Un(dimprb))
         ALLOCATE(F_Ur(dimprb))
         ALLOCATE(Ur(dimprb))
         ALLOCATE(RHS(dimprb))
         CALL this%linsolver%Construct(dimprb)
         CALL this%linsolver%SetMatMult(JacFreeAx)
         
         
         IF (PRESENT(pc)) THEN
            SELECT CASE(pc)
               CASE('GMRES')
                  CALL this%linsolver%SetPCMult(GMRES_Preco)
                  CALL PCsolver%Construct(dimprb,15)
                  CALL PCsolver%SetMatMult(JacFreeAx_PRECO)
                  CALL PCsolver%SetMaxIter(30)           !old: 15
                  CALL this%linsolver%SetMaxiter(60)     !old: 30
                  this%PC_GMRES = .TRUE.
               CASE DEFAULT
                  WRITE(*,*) pc, ' preconditioner not found. No preconditioner will be used for GMRES'
             END SELECT
          END IF
       
      END SUBROUTINE
!/////////////////////////////////////////////////////////////////////////           
      SUBROUTINE Destruct_JFNK_Integrator(this)
         CLASS(JFNK_Integrator), INTENT(INOUT)     :: this
         !Used module variables: p_F, RHS, Ur
         
         DEALLOCATE(this%Un)
         DEALLOCATE(F_Ur)
         DEALLOCATE(Ur)
         DEALLOCATE(RHS)
         p_F => NULL()
         CALL this%linsolver%Destruct()
         ! arueda: destruct PCSolver here??
      
      END SUBROUTINE
!/////////////////////////////////////////////////////////////////////////     
      SUBROUTINE JFNK_Integrate(this)
         CLASS(JFNK_Integrator), INTENT(INOUT)     :: this
         LOGICAL                                   :: ISFIRST = .TRUE.
         !Used module variables: F_Ur, Ur                        
         
         IF(ISFIRST) THEN
            this%last_Dt = this%Dt
            ISFIRST = .FALSE.
         END IF
         
         this%nsteps = 0                  !Reset iteration and step counters

         inner_t = 0._RP                  ! Reset iteration time
         
         Comp_Dt = this%Dt                ! arueda: Here, I don't use lastDt, because I want the code to run time-accurate cases
         
         DO ! While inner_t < Dt
            Ur = this%Un
            F_Ur = p_F(Ur,time+Comp_Dt)
            this%norm_G = UNORM(F_Ur,'L2')
            this%norm_x0 = UNORM(Ur,'L2')       
            this%norm_x_abs = this%norm_x0
            this%norm_x_rel = 1._RP
            this%n_linsolver_iter = 0
            n_preco_iter = 0
            IF (this%PRINT_NEWTON_INFO) THEN
               CALL Newton_Output(this,'HEADER')
               CALL Newton_Output(this, 'FIRST')
            ENDIF
            
            ! Newton solver will use value of Comp_Dt at this point to perform computations
            CALL NewtonSolver(this)
            IF((this%linsolver%CONVERGED) .AND. (this%CONVERGED)) THEN  
               this%nsteps = this%nsteps + 1
               inner_t = inner_t + Comp_Dt                                  ! Increases integration time
               time = time + Comp_Dt
               this%Un = Ur                                                 ! Updates Un after each successful Newton Solver call
             
!WRITE(*,*)   'linsolver iter: ', this%n_linsolver_iter      
!WRITE(*,'(A)', ADVANCE='NO') 'Inner_t / out_Dt '
!WRITE(*,'(2(1p,e10.2))') inner_t,   this%Dt
!WRITE(*,'(A)', ADVANCE='NO') '  Comp Dt = '
!WRITE(*,'(1p,e8.2)') Comp_Dt 
             
              ! Check if Outer Dt has been reached
              IF (ABS((this%Dt - inner_t)) < 10 * EPSILON(1._RP)) THEN       ! If outer Dt is reached, the time integration is done
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  ! arueda: uncomment but as messages
!~                   WRITE(*,*)
!~                   WRITE(*,*) 'TIME STEP DONE'
!~                   WRITE(*,'(A)', ADVANCE='NO') ' # sub-steps:  '
!~                   WRITE(*, '(I4)'), this%nsteps
!~                   WRITE(*,'(A)', ADVANCE='NO') ' # ComputeTimeDerivative calls : '
!~                   WRITE(*, '(I6)'), this%n_linsolver_iter + n_preco_iter
                  !
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                  WRITE(*,'(A)', ADVANCE='NO') ' inner time   /  DT  '
!                  WRITE(*, '(2(1p,e21.10))'), inner_t, this%Dt 
!                  WRITE(*,*)
!~                   WRITE(*,*) !arueda: uncomment as message
                  RETURN                                             
               ENDIF
               
              !Increase Comp_Dt if good convergence in previous step
              IF ((this%n_linsolver_iter < 40) .AND. (Comp_Dt < this%Dt)) THEN
                  Comp_Dt = Comp_Dt * REAL(this%substep_inc_ratio, KIND=RP) 
                  this%last_Dt = Comp_Dt
                  WRITE(*,'(A)', ADVANCE='NO') 'Increasing Substep Delta T... Control_Dt = '
                  WRITE(*,'(1p,e8.2)', ADVANCE='NO') this%Dt
                  WRITE(*,'(A)', ADVANCE='NO') '  Comp Dt = '
                  WRITE(*,'(1p,e8.2)') Comp_Dt 
              ENDIF
              
              ! Adjust Comp_Dt to prevent inner_t be greater than outer Dt 
              IF ( (Comp_DT > (this%Dt - inner_t ))) THEN  ! Adjusts inner dt to achieve exact outer Dt in the last substep
                !  WRITE(*,*) 'Esta debería ser la última sub-iteración....'
                  Comp_Dt = this%Dt - inner_t                                ! Do not update last_Dt value
                !  WRITE(*,*) 'Comp_DT =', Comp_Dt 
              ENDIF
            
            ELSE ! Reduce Comp Dt if newton or linear solver did not converge
               print*, 'arueda: comp_dt: ', Comp_Dt
               Comp_Dt = Comp_Dt /REAL(this%substep_red_ratio, KIND=RP) 
               print*, 'arueda: comp_dt: ', Comp_Dt
               this%last_Dt = Comp_Dt
               this%nsteps = 0                                 
               WRITE(*,'(A)', ADVANCE='NO') 'Reducing Substep Delta T... Control_Dt = '
               WRITE(*,'(1p,e8.2)', ADVANCE='NO') this%Dt
               WRITE(*,'(A)', ADVANCE='NO') ' Comp Dt = '
               WRITE(*,'(1p,e8.2)') Comp_Dt 
            ENDIF
         END DO
     
      END SUBROUTINE JFNK_Integrate
!/////////////////////////////////////////////////////////////////////////           
      SUBROUTINE NewtonSolver(this)
         CLASS(JFNK_Integrator), INTENT(INOUT)     :: this
         INTEGER                                   :: i
                    
         DO i = 1, this%MAX_NEWTON_ITER
            this%niter = i
            
            CALL NewtonInnerIt(this)
            
            IF(.NOT. this%linsolver%CONVERGED) THEN                  
               RETURN
            ENDIF
            
            F_Ur = p_F(Ur,time+Comp_Dt)                              ! Update F(Ur)

            this%norm_x_abs = UNORM(this%linsolver%x,'L2')           ! Compute L_2 norm of x (U_r correction)
            this%norm_x_rel = this%norm_x_abs / this%norm_x0         ! Adimensionalize norm with the firsrt norm
            this%norm_G = UNORM(RHS, 'L2')/ this%DimPrb           
            IF (this%PRINT_NEWTON_INFO) THEN
               CALL Newton_Output(this)
            ENDIF
            IF ((this%norm_x_rel .LT. this%reltol) &
                  .OR. (this%norm_x_abs .LT. this%abstol)) THEN      !Check for convergence in DeltaU
               this%CONVERGED = .TRUE.
               RETURN
            ENDIF
            IF (this%norm_G .LT. this%abstol) THEN                   !Check for convergence in F
               this%CONVERGED = .TRUE.
               RETURN
            ENDIF
         END DO   
         this%CONVERGED = .FALSE.                                    !If MAX_NEWTON_ITER reached iteration is set as not converged       
         
      END SUBROUTINE NewtonSolver
!/////////////////////////////////////////////////////////////////////////     
      SUBROUTINE NewtonInnerIt(this)
         CLASS(JFNK_Integrator), INTENT(INOUT)     :: this
         ! Used module variables: Ur, RHS
         
         IF (this%niter == 1) THEN ! Use a tolerance based on F(U_r) and Dt for the first iteration         
            gmres_tol = MIN(this%norm_x_abs* 1.e-2_RP, this%norm_G * inner_Dt * 1.e-1_RP )
         ELSE ! Use a tolerance based on the previous solution and the convergence ratio (2) of Newton solver
            gmres_tol = this%norm_x_abs* 1.e-3_RP
         END IF
         
         IF (this%PC_GMRES) THEN  ! Here the tolerance of GMRES used for precondiotioning is set TODO: aquí hay que mirar un poco más
            CALL PCsolver%SetTol(gmres_tol)
            PCSolver%x0 = 0._RP
         END IF
         
         eps = SQRT((1._RP + this%norm_x_abs) * this%abstol)    
            
         CALL ComputeRHS(this)                     
         CALL this%linsolver%SetRHS(RHS)
         CALL this%linsolver%SetTol(gmres_tol)
         CALL this%linsolver%Solve
         this%n_linsolver_iter = this%n_linsolver_iter + this%linsolver%niter
         IF (this%linsolver%ERROR_CODE  .NE. 0) THEN
            print*, 'error in linear solver: ',  this%linsolver%ERROR_CODE
         ENDIF
         IF (.NOT. this%linsolver%CONVERGED) THEN
            PRINT*, '*** WARNING: linsolver max iter reached'
            RETURN
         ENDIF
         
         Ur = Ur + this%linsolver%x
          
      END SUBROUTINE NewtonInnerIt
!/////////////////////////////////////////////////////////////////////////          
      SUBROUTINE SetTimeDerivative(this, F)
         CLASS(JFNK_Integrator), INTENT(INOUT)     :: this
         INTERFACE
            FUNCTION F(u, time) RESULT(Func)
               IMPORT RP
               REAL(KIND = RP), INTENT(IN)  :: u(:)
               REAL(KIND = RP)              :: Func(size(u))
               REAL(KIND = RP)              :: time
            END FUNCTION
         END INTERFACE
         !Used module variables: p_F
         
         p_F    => F ! Shared Pointer to TimeDerivative in this module
         this%F_SET = .TRUE.
      
      END SUBROUTINE SetTimeDerivative
!/////////////////////////////////////////////////////////////////////////          
      SUBROUTINE SetDt(this, Dt_in)
         CLASS(JFNK_Integrator), INTENT(INOUT)     :: this
         REAL(KIND = RP), INTENT(IN)               :: Dt_in
         !Used module variables: Dt  
         
         this%Dt = Dt_in
         inner_Dt = Dt_in   ! TODO: This reset inner_Dt if a new outter Dt is set. Is this a good idea??????
         this%inner_Dt = Dt_in
      
      END SUBROUTINE SetDt
!/////////////////////////////////////////////////////////////////////////    
      !---------------------------------------------------------------- arueda
      SUBROUTINE SetTime(this, t)
      !>   Sets initial time for time-stepping 
      !---------------------------------------------------------------- 
         CLASS(JFNK_Integrator), INTENT(INOUT)     :: this
         REAL(KIND = RP), INTENT(IN)               :: t
         !Used module variables: time  
         time = t
      !---------------------------------------------------------------- 
      END SUBROUTINE SetTime
      !---------------------------------------------------------------- 
!/////////////////////////////////////////////////////////////////////////       
      SUBROUTINE SetUn(this, Un)
         CLASS(JFNK_Integrator), INTENT(INOUT)     :: this
         REAL(KIND = RP), INTENT(IN)               :: Un(:)
        
         this%Un = Un
      END SUBROUTINE SetUn
!/////////////////////////////////////////////////////////////////////////          
      SUBROUTINE GetUnp1(this, Unp1)
         CLASS(JFNK_Integrator), INTENT(INOUT)     :: this
         REAL(KIND = RP), INTENT(OUT)              :: Unp1(:)
        
         Unp1 = this%Un
      END SUBROUTINE GetUnp1
!/////////////////////////////////////////////////////////////////////////          
      SUBROUTINE ComputeRHS(this)
         CLASS(JFNK_Integrator), INTENT(INOUT)     :: this
         !Module variables: RHS, Ur, Dt_Comp, F_Ur
         
         RHS = (Ur - this%Un) / Comp_Dt - F_Ur   ! arueda: this is also defined only for BDF1... 
      
      END SUBROUTINE ComputeRHS
!/////////////////////////////////////////////////////////////////////////  
      SUBROUTINE JacFreeAx(x, Ax)
         REAL(KIND = RP), INTENT(IN)         :: x(:)
         REAL(KIND = RP), INTENT(OUT)        :: Ax(:)
         !Module variables: Dt_Comp, p_F, eps, F_Ur
          
         eps = 1e-8_RP * (1._RP + UNORM(x,'L2'))
         Ax = ( p_F(Ur + x * eps,time+Comp_Dt) - F_Ur)  / eps  - x / Comp_Dt                          !First order   ! arueda: this is also defined only for BDF1
         !Ax = ( p_F(Ur + x * eps) - p_F(Ur - x * eps))  /(2._RP * eps)  - x / Comp_Dt   !Second order
      
      END SUBROUTINE JacFreeAx
!/////////////////////////////////////////////////////////////////////////          
      SUBROUTINE JacFreeAx_PRECO(x, Ax)
         REAL(KIND = RP), INTENT(IN)         :: x(:)
         REAL(KIND = RP), INTENT(OUT)        :: Ax(:)
         !Module variables: Dt_Comp, p_F, eps, F_Ur
          
         eps = 1e-8_RP * (1._RP + UNORM(x,'L2'))
         Ax = ( p_F(Ur + x * eps,time+Comp_Dt) - F_Ur)  / eps  - x / (Comp_Dt)                          !First order   ! arueda: this is also defined only for BDF1
      
      END SUBROUTINE JacFreeAx_PRECO
!/////////////////////////////////////////////////////////////////////////          
      SUBROUTINE GMRES_Preco(v, Pv)
         ! Returns  the preconditioning product   Pv = P⁻¹ * v 
         REAL(KIND = RP), INTENT(IN)         :: v(:)
         REAL(KIND = RP), INTENT(OUT)        :: Pv(:)
         
          
        CALL PCSolver%SetRHS(v)
        CALL PCSolver%Solve
        CALL PCSolver%Settol(1e-1_RP)                 ! TODO: aquí habrá que poner algo más elaborado
        Pv = PCSolver%x
        n_preco_iter = n_preco_iter + PCSolver%niter
                        
      END SUBROUTINE GMRES_Preco
!/////////////////////////////////////////////////////////////////////////          
      SUBROUTINE Newton_Output(this, opt)
         CLASS(JFNK_Integrator), INTENT(INOUT)     :: this
         CHARACTER(LEN = *), OPTIONAL              :: opt
         CHARACTER(LEN = 50)                       :: fmt1
         
         fmt1 = '(I8,1p,E14.2,E11.2,E13.2,I10,E17.2,E10.2)'
         IF(.NOT. PRESENT(opt)) THEN
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! arueda: uncomment but with message
!~             WRITE(*,fmt1)  this%niter, this%norm_x_abs,  this%norm_x_rel, this%norm_G, &
!~                            this%linsolver%niter, this%linsolver%res, this%linsolver%tol
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ELSE
            SELECT CASE (opt)
               CASE ('HEADER')
                !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                ! arueda: uncomment but with message
!~                   WRITE(*,*)
!~                   WRITE(*,*) ' NewtonIt | Newton: AbsErr  RelErr | Nrm(F_Ur) | #GMRESiter | linsolver Res  Tol'
                !
                !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               CASE ('FIRST')
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  ! arueda: : uncomment but with message
!~                   WRITE(*,fmt1)  0, this%norm_x_abs,  this%norm_x_rel, this%norm_G, &
!~                            -1, -1., -1. 
                  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               CASE DEFAULT
            END SELECT
         ENDIF
      
      END SUBROUTINE Newton_Output
!/////////////////////////////////////////////////////////////////////////          
      FUNCTION UNORM(x, opt) RESULT(norm)
         REAL(KIND = RP)                           :: x(:)
         REAL(KIND = RP)                           :: norm
         CHARACTER(LEN = *), OPTIONAL              :: opt
         
         IF(.NOT. PRESENT(opt)) THEN  !Default norm L2-norm
            norm = NORM2(x)
         ELSE
            SELECT CASE (opt)
               CASE('L2')
                  norm = NORM2(x)
               CASE('Linf')
                  norm = MAXVAL(ABS(x))
               CASE DEFAULT
                  PRINT*, 'ERROR: selected norm', opt ,' not implemented'
                  STOP
            END SELECT 
         END IF
      END FUNCTION

END MODULE JFNKClass
