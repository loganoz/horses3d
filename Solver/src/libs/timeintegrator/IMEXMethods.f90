!
!//////////////////////////////////////////////////////
!
!   @File:    IMEXMethods.f90
!   @Author:  Juan (juan.manzanero@upm.es)
!   @Created: Tue Apr 17 16:55:49 2018
!   @Last revision date: Sun May 19 16:54:10 2019
!   @Last revision author: Andrés Rueda (am.rueda@upm.es)
!   @Last revision commit: 8958d076d5d206d1aa118cdd3b9adf6d8de60aa3
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
MODULE IMEXMethods
   USE SMConstants                  
   USE DGSEMClass,                  ONLY: DGSem
   USE ElementClass,                ONLY: Element
   USE PhysicsStorage
   use FluidData
   use HexMeshClass
   use MKLPardisoSolverClass
   USE CSRMatrixClass
   USE FTValueDictionaryClass
   use TimeIntegratorDefinitions
   use MatrixClass
   use DGSEMClass, only: ComputeTimeDerivative_f
   use BoundaryConditions, only: C_BC, NS_BC
   implicit none
   
   PRIVATE                          
   PUBLIC TakeIMEXStep
   
   real(kind=RP) :: time               ! Time at the beginning of each inner(!) time step
   logical       :: computeA = .TRUE.  ! Compute Jacobian? (only valid if it is meant to be computed according to the convergence)
!
!  ========
   contains
!  ========
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!   
   subroutine TakeIMEXStep(sem, t, dt, controlVariables, ComputeTimeDerivative)
      implicit none
      TYPE(DGSem),                  intent(inout)           :: sem                  !<>DGSem class with solution storage 
      REAL(KIND=RP),                intent(in)              :: t                    !< Time at the beginning of time step
      REAL(KIND=RP),                intent(in)              :: dt                   !< Initial (outer) time step (can internally, the subroutine can use a smaller one depending on convergence)
      TYPE(FTValueDictionary),      intent(in)              :: controlVariables     !< Input file variables
      procedure(ComputeTimeDerivative_f)                            :: ComputeTimeDerivative

      select case(solver)
      case(NAVIERSTOKES_SOLVER, INCNS_SOLVER)
         print*, "IMEX solver not implemented for monophase Navier-Stokes equations"
         stop

      case(CAHNHILLIARD_SOLVER)
         call TakeIMEXEulerStep_CH(sem, t , dt , controlVariables, ComputeTimeDerivative)

      case(MULTIPHASE_SOLVER)
         call TakeIMEXBDF2Step_MU(sem, t , dt , controlVariables, ComputeTimeDerivative)

      case default
         print*, "Solver not recognized"
         stop
      end select

   end subroutine TakeIMEXStep

   SUBROUTINE TakeIMEXEulerStep_CH (sem, t , dt , controlVariables, ComputeTimeDerivative)
      IMPLICIT NONE
      TYPE(DGSem),                  intent(inout) :: sem                  !<>DGSem class with solution storage
      REAL(KIND=RP),                intent(in)    :: t                    !< Time at the beginning of time step
      REAL(KIND=RP),                intent(in)    :: dt                   !< Initial (outer) time step (can internally, the subroutine can use a smaller one depending on convergence)
      TYPE(FTValueDictionary),      intent(in)    :: controlVariables     !< Input file variables
      procedure(ComputeTimeDerivative_f)          :: ComputeTimeDerivative
!
!     ---------------
!     Local variables
!     ---------------
!
      type(MKLPardisoSolver_t), save           :: linsolver
      integer                                  :: nelm, DimPrb, globalDimPrb
      logical, save                            :: isfirst = .TRUE.
      REAL(KIND=RP), DIMENSION(:), allocatable :: U_n                                   !Solution at the beginning of time step (even for inner time steps)
      logical                                  :: TimeAccurate = .true.
      integer                                  :: nEqnJac, nGradJac
      REAL(KIND=RP), DIMENSION(3)              :: a = (/0.0_RP       , -5.0_RP /9.0_RP , -153.0_RP/128.0_RP/)
      REAL(KIND=RP), DIMENSION(3)              :: b = (/0.0_RP       ,  1.0_RP /3.0_RP ,    3.0_RP/4.0_RP  /)
      REAL(KIND=RP), DIMENSION(3)              :: c = (/1.0_RP/3.0_RP,  15.0_RP/16.0_RP,    8.0_RP/15.0_RP /)
      real(kind=RP)                            :: tk
      integer                                  :: k, id
      SAVE DimPrb, nelm, TimeAccurate

#if (!defined(FLOW)) && defined(CAHNHILLIARD)

      nEqnJac = NCOMP
      nGradJac = NCOMP
      
      IF (isfirst) THEN           
!
!        ***********************************************************************
!           Construct the Jacobian, and perform the factorization in the first
!           call.
!        ***********************************************************************
!
         isfirst = .FALSE.
         nelm = SIZE(sem%mesh%elements)
         DimPrb = sem % NDOF * NCOMP
         globalDimPrb = sem % totalNDOF * NCOMP
         
         ALLOCATE(U_n(0:Dimprb-1))

         CALL linsolver%construct(DimPrb,globalDimPrb, nEqnJac,controlVariables,sem, IMEX_MatrixShift) 

         call linsolver%ComputeAndFactorizeJacobian(nEqnJac,nGradJac, ComputeTimeDerivative, dt, 1.0_RP)
         
      ENDIF
      
      time = t

!     TODO do i need this?      
      !CALL sem % GetQ(U_n)      !stores sem%mesh%elements(:)% storage % Q in Vector U_n
!
!     Compute the non linear time derivative
!     -------------------------------------- 
      call ComputeTimeDerivative(sem % mesh, sem % particles, time, CTD_IMEX_EXPLICIT)
!
!     Compute the RHS
!     ---------------
      call ComputeRHS(sem, dt, nelm, linsolver)
!
!     Solve the linear system
!     -----------------------
      call linsolver % SolveLUDirect
!
!     Return the computed state vector to storage
!     -------------------------------------------
      call sem % mesh % storage % local2GlobalQ (sem % NDOF)
      sem % mesh % storage % Q = linsolver % x
      call sem % mesh % storage % global2LocalQ
!
!     Compute the standard time derivative to get residuals
!     -----------------------------------------------------
      call ComputeTimeDerivative(sem % mesh, sem % particles, time, CTD_IGNORE_MODE)
      
#endif
   END SUBROUTINE TakeIMEXEulerStep_CH

   SUBROUTINE TakeIMEXRKStep_MU (sem, t , dt , controlVariables, ComputeTimeDerivative)
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!           This subroutine implements an IMplicit-EXplicit third order RK scheme
!        See the details in: https://doi.org/10.1016/S0168-9274(97)00056-1
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      implicit none
      type(DGSem),                  intent(inout) :: sem                  !<>DGSem class with solution storage
      real(kind=RP),                intent(in)    :: t                    !< Time at the beginning of time step
      real(kind=RP),                intent(in)    :: dt                   !< Initial (outer) time step (can internally, the subroutine can use a smaller one depending on convergence)
      TYPE(FTValueDictionary),      intent(in)    :: controlVariables     !< Input file variables
      procedure(ComputeTimeDerivative_f)          :: ComputeTimeDerivative
!
!     ---------------
!     Local variables
!     ---------------
!
      type(MKLPardisoSolver_t), save           :: linsolver
      integer                                  :: nelm, DimPrb, globalDimPrb
      logical, save                            :: isfirst = .TRUE.
      logical                                  :: TimeAccurate = .true.
      integer                                  :: nEqnJac, nGradJac
      integer,       parameter                 :: NSTAGES = 2
      real(kind=RP), parameter                 :: gamma = (2.0_RP - sqrt(2.0_RP))/(2.0_RP)
      real(kind=RP), parameter                 :: delta = -2.0_RP * sqrt(2.0_RP) / 3.0_RP
      real(kind=RP), parameter                 :: a(2,2) = RESHAPE((/gamma,1.0_RP-gamma,0.0_RP,gamma/),(/2,2/))
      real(kind=RP), parameter                 :: b(2) = [1.0_RP-gamma,gamma]
      real(kind=RP), parameter                 :: hatA(3,3) = RESHAPE((/0.0_RP,gamma,delta,0.0_RP,0.0_RP,1-delta,0.0_RP,0.0_RP,0.0_RP/),(/3,3/))
      real(Kind=RP), parameter                 :: hatB(3) = [0.0_RP, 1-gamma, gamma]
!      real(kind=RP), parameter                 :: gamma = (3.0_RP + sqrt(3.0_RP))/(6.0_RP)
!      real(kind=RP), parameter                 :: a(2,2) = RESHAPE((/gamma,1.0_RP-2.0_RP*gamma,0.0_RP,gamma/),(/2,2/))
!      real(kind=RP), parameter                 :: b(2) = [0.5_RP,0.5_RP]
!      real(kind=RP), parameter                 :: c(2) = [gamma, 1.0_RP-gamma]
!      real(kind=RP), parameter                 :: hatA(3,3) = RESHAPE((/0.0_RP,gamma,gamma-1.0_RP,0.0_RP,0.0_RP,2.0_RP*(1.0_RP-gamma),0.0_RP,0.0_RP,0.0_RP/),(/3,3/))
!      real(Kind=RP), parameter                 :: hatB(3) = [0.0_RP, 0.5_RP, 0.5_RP]
!      real(Kind=RP), parameter                 :: hatC(3) = [0.0_RP,gamma,1.0_RP-gamma]
      real(kind=RP)                            :: tk
      integer                                  :: k, id, s
      SAVE DimPrb, nelm, TimeAccurate

#if defined(MULTIPHASE) && defined(CAHNHILLIARD)
      nEqnJac = NCOMP
      nGradJac = NCOMP

      IF (isfirst) THEN           
!
!        ***********************************************************************
!           Construct the Jacobian, and perform the factorization in the first
!           call.
!        ***********************************************************************
!
         isfirst = .FALSE.
         nelm = SIZE(sem%mesh%elements)
         DimPrb = sem % NDOF * NCOMP
         globalDimPrb = sem % totalNDOF * NCOMP
         
         CALL linsolver%construct(DimPrb,globalDimPrb, nEqnJac,controlVariables,sem, IMEX_MatrixShift) 
         call linsolver%ComputeAndFactorizeJacobian(nEqnJac,nGradJac, ComputeTimeDerivative, gamma*dt, 1.0_RP)
         
      ENDIF
!
!     First stage
!     -----------
      call ComputeTimeDerivative(sem % mesh, sem % particles, time, CTD_IMEX_EXPLICIT)

!$omp parallel do schedule(runtime)
      do id = 1, sem % mesh % no_of_elements
         sem % mesh % elements(id) % storage % RKSteps(1) % hatK = sem % mesh % elements(id) % storage % QDot 
      end do
!$omp end parallel do
!
!     Rest of the stages            
!     ------------------
      do s = 1, NSTAGES
!
!        First step: get the implicit Runge-Kutta coefficient for stage-s
!        ----------------------------------------------------------------
!
!        Load Q
!        ------
!$omp parallel do schedule(runtime) private(k)
         do id = 1, sem % mesh % no_of_elements
!
!         - Explicit coefficients
            do k = 1, s
               sem % mesh % elements(id) % storage % Q = sem % mesh % elements(id) % storage % Q + dt*hatA(s+1,k)*sem % mesh % elements(id) % storage % RKSteps(k) % hatK
            end do
!
!         - Implicit coefficients
            do k = 1,s-1
               sem % mesh % elements(id) % storage % Q(IMC,:,:,:) = sem % mesh % elements(id) % storage % Q(IMC,:,:,:) &
                                                                    + dt*A(s,k)*sem % mesh % elements(id) % storage % RKSteps(k) % K(1,:,:,:)
            end do
         end do
!$omp end parallel do
!
!        Get the Right Hand Side for the linear solver
!        ---------------------------------------------
         call ComputeTimeDerivative(sem % mesh, sem % particles, time, CTD_IMEX_IMPLICIT)
!
!        Solve the linear system
!        -----------------------
         call sem % mesh % SetStorageToEqn(C_BC)
         call ComputeRHS_MURK(sem, dt, nelm, linsolver)
         call linsolver % SolveLUDirect
         call sem % mesh % storage % local2GlobalQ (sem % NDOF)
         sem % mesh % storage % Q = linsolver % x
         call sem % mesh % storage % global2LocalQ
         call sem % mesh % SetStorageToEqn(NS_BC)
!
!        Set the implicit Runge-Kutta coefficient (stored in c) and add it to the solution
!        ---------------------------------------------------------------------------------
!$omp parallel do schedule(runtime)
         do id = 1, sem % mesh % no_of_elements
            sem % mesh % elements(id) % storage % RKSteps(s) % K(1,:,:,:) = sem % mesh % elements(id) % storage % c(1,:,:,:)
            sem % mesh % elements(id) % storage % Q(IMC,:,:,:) = sem % mesh % elements(id) % storage % Q(IMC,:,:,:) &
                                                                 + dt*A(s,s)*sem % mesh % elements(id) % storage % RKSteps(s) % K(1,:,:,:)
         end do
!$omp end parallel do
!
!        Compute QDot -> RKStep(s+1) % hatK
!        ----------------------------------
         call ComputeTimeDerivative(sem % mesh, sem % particles, time, CTD_IMEX_EXPLICIT)
!
!        Set new explicit Runge-Kutta coefficient + recover the original solution
!        ------------------------------------------------------------------------
!$omp parallel do schedule(runtime) private(k)
         do id = 1, sem % mesh % no_of_elements
            sem % mesh % elements(id) % storage % RKSteps(s+1) % hatK = sem % mesh % elements(id) % storage % QDot
!
!         - Explicit coefficients
            do k = 1, s
               sem % mesh % elements(id) % storage % Q = sem % mesh % elements(id) % storage % Q - dt*hatA(s+1,k)*sem % mesh % elements(id) % storage % RKSteps(k) % hatK
            end do
!
!         - Implicit coefficients
            do k = 1,s
               sem % mesh % elements(id) % storage % Q(IMC,:,:,:) = sem % mesh % elements(id) % storage % Q(IMC,:,:,:) &
                                                                    - dt*A(s,k)*sem % mesh % elements(id) % storage % RKSteps(k) % K(1,:,:,:)
            end do

         end do
!$omp end parallel do
      end do
!
!     Perform the time step
!     ---------------------
!$omp parallel do schedule(runtime) private(s)
      do id = 1, sem % mesh % no_of_elements
!
!      - Explicit coefficients
         do s = 1, NSTAGES + 1
            sem % mesh % elements(id) % storage % Q = sem % mesh % elements(id) % storage % Q + dt*hatB(s)*sem % mesh % elements(id) % storage % RKSteps(s) % hatK
         end do
!
!      - Implicit coefficients
         do s = 1, NSTAGES
            sem % mesh % elements(id) % storage % Q(IMC,:,:,:) = sem % mesh % elements(id) % storage % Q(IMC,:,:,:) &
                                                                 + dt*B(s)*sem % mesh % elements(id) % storage % RKSteps(s) % K(1,:,:,:)
   
         end do
      end do
!$omp end parallel do
!
!     The "good" residuals CTD call
!     -----------------------------
      call ComputeTimeDerivative(sem % mesh, sem % particles, time, CTD_IGNORE_MODE)

#endif
   end subroutine TakeIMEXRKStep_MU

   SUBROUTINE TakeIMEXBDF2Step_MU (sem, t , dt , controlVariables, ComputeTimeDerivative)
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!           This subroutine implements an IMplicit-EXplicit third order RK scheme
!        See the details in: https://doi.org/10.1016/S0168-9274(97)00056-1
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      implicit none
      type(DGSem),                  intent(inout) :: sem                  !<>DGSem class with solution storage
      real(kind=RP),                intent(in)    :: t                    !< Time at the beginning of time step
      real(kind=RP),                intent(in)    :: dt                   !< Initial (outer) time step (can internally, the subroutine can use a smaller one depending on convergence)
      TYPE(FTValueDictionary),      intent(in)    :: controlVariables     !< Input file variables
      procedure(ComputeTimeDerivative_f)          :: ComputeTimeDerivative
!
!     ---------------
!     Local variables
!     ---------------
!
      type(MKLPardisoSolver_t), save           :: linsolver
      integer                                  :: nelm, DimPrb, globalDimPrb
      logical, save                            :: isfirst = .TRUE.
      logical, save                            :: isSecond = .false.
      logical                                  :: TimeAccurate = .true.
      integer                                  :: nEqnJac, nGradJac
      integer,       parameter                 :: NSTAGES = 2
      real(kind=RP)                            :: tk
      integer                                  :: k, id, s
      SAVE DimPrb, nelm, TimeAccurate

#if defined(MULTIPHASE) && defined(CAHNHILLIARD)
      nEqnJac = NCOMP
      nGradJac = NCOMP

      if (isSecond) then
!
!        Store the explicit approximation of Q
!        -------------------------------------
!$omp parallel do schedule(runtime)
         do id = 1, size(sem % mesh % elements)
!
!           Compute Q^{*,n+1}
!           -----------------
            sem % mesh % elements(id) % storage % Q = 2.0_RP * sem % mesh % elements(id) % storage % Q - sem % mesh % elements(id) % storage % PrevQ(1) % Q

         end do
!$omp end parallel do
!
!        Compute the time derivative (Explicit)
!        --------------------------------------
         call ComputeTimeDerivative(sem % mesh, sem % particles, time, CTD_IMEX_EXPLICIT)
!
!        Change the y^{*,n+1} in Q to \hat{y}
!        ------------------------------------
!$omp parallel do schedule(runtime)
         do id = 1, size(sem % mesh % elements)
            sem % mesh % elements(id) % storage % Q = sem % mesh % elements(id) % storage % Q + 0.5_RP * sem % mesh % elements(id) % storage % prevQ(1) % Q
            sem % mesh % elements(id) % storage % PrevQ(1) % Q = 0.5_RP * (sem % mesh % elements(id) % storage % Q + 0.5_RP * sem % mesh % elements(id) % storage % PrevQ(1) % Q)
         end do
!$omp end parallel do
!
!        Perform the implicit time-step
!        ------------------------------
         call ComputeRHS_MUBDF2(sem, dt, nelm, linsolver)
         call linsolver % SolveLUDirect
         call sem % mesh % SetStorageToEqn(C_BC)
         call sem % mesh % storage % local2GlobalQ (sem % NDOF)
         sem % mesh % storage % Q = linsolver % x
         call sem % mesh % storage % global2LocalQ
         call sem % mesh % SetStorageToEqn(NS_BC)
!
!$omp parallel do schedule(runtime)
         do id = 1, sem % mesh % no_of_elements
!
!           Set the concentration from the linear solver solution
!           -----------------------------------------------------
            sem % mesh % elements(id) % storage % Q(IMC,:,:,:) = sem % mesh % elements(id) % storage % c(1,:,:,:)
!
!           Perform a time-step on the rest of the variables
!           ------------------------------------------------
            sem % mesh % elements(id) % storage % Q(IMC+1:,:,:,:) = (2.0_RP / 3.0_RP ) * (sem % mesh % elements(id) % storage % Q(IMC+1:,:,:,:) & 
                                                                    + dt*sem % mesh % elements(id) % storage % QDot(IMC+1:,:,:,:))
         end do 
!$omp end parallel do

      end if

      IF (isfirst) THEN           
!
!        ***********************************************************************
!           Construct the Jacobian, and perform the factorization in the first
!           call.
!        ***********************************************************************
!
         isfirst = .false.
         isSecond = .true.
         nelm = SIZE(sem%mesh%elements)
         DimPrb = sem % NDOF * NCOMP
         globalDimPrb = sem % totalNDOF * NCOMP
         
         CALL linsolver%construct(DimPrb,globalDimPrb, nEqnJac,controlVariables,sem, IMEX_MatrixShift) 
         call linsolver%ComputeAndFactorizeJacobian(nEqnJac,nGradJac, ComputeTimeDerivative, dt, 1.0_RP)

         call ComputeTimeDerivative(sem % mesh, sem % particles, time, CTD_IMEX_EXPLICIT)
   
         call ComputeRHS_MUBDF2(sem, dt, nelm, linsolver)
         call linsolver % SolveLUDirect
         call sem % mesh % SetStorageToEqn(C_BC)
         call sem % mesh % storage % local2GlobalQ (sem % NDOF)
         sem % mesh % storage % Q = linsolver % x
         call sem % mesh % storage % global2LocalQ
         call sem % mesh % SetStorageToEqn(NS_BC)
!
!$omp parallel do schedule(runtime)
         do id = 1, sem % mesh % no_of_elements
!
!           Set old solution
!           ----------------
            sem % mesh % elements(id) % storage % PrevQ(1) % Q = sem % mesh % elements(id) % storage % Q
!
!           Set the concentration from the linear solver solution
!           -----------------------------------------------------
            sem % mesh % elements(id) % storage % Q(IMC,:,:,:) = sem % mesh % elements(id) % storage % c(1,:,:,:)
!
!           Perform a time-step on the rest of the variables
!           ------------------------------------------------
            sem % mesh % elements(id) % storage % Q(IMC+1:,:,:,:) = sem % mesh % elements(id) % storage % Q(IMC+1:,:,:,:) & 
                                                                    + dt*sem % mesh % elements(id) % storage % QDot(IMC+1:,:,:,:)
         end do 
!$omp end parallel do
!
!        Change the coefficient in the Jacobian for the BDF2
!        ---------------------------------------------------
         call linsolver % ReFactorizeJacobian
         print*, "Jacobian re-factorized"
         
      ENDIF
!
!     The "good" residuals CTD call
!     -----------------------------
      call ComputeTimeDerivative(sem % mesh, sem % particles, time, CTD_IGNORE_MODE)

#endif
   end subroutine TakeIMEXBDF2Step_MU
!  
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE ComputeRHS(sem, dt, nelm, linsolver )
      implicit none
      TYPE(DGSem),                intent(in)       :: sem
      REAL(KIND=RP),              intent(in)       :: dt
      integer,                    intent(in)       :: nelm
      CLASS(GenericLinSolver_t),  intent(inout)   :: linsolver

      integer                                      :: Nx, Ny, Nz, l, i, j, k, elmnt, counter   
      REAL(KIND=RP)                                :: value

#if defined(CAHNHILLIARD)
      counter = 0
      DO elmnt = 1, nelm
         Nx = sem%mesh%elements(elmnt)%Nxyz(1)
         Ny = sem%mesh%elements(elmnt)%Nxyz(2)
         Nz = sem%mesh%elements(elmnt)%Nxyz(3)
         DO k = 0, Nz
            DO j = 0, Ny
               DO i = 0, Nx
                  DO l = 1,NCOMP
                     value = sem%mesh%elements(elmnt)%storage%c(l,i,j,k) + &
                          dt*sem%mesh%elements(elmnt)%storage%cDot(l,i,j,k)
                     counter =  counter + 1
                     CALL linsolver%SetRHSValue(counter, value)
                     
                  END DO
               END DO
            END DO
         END DO
      END DO

!      CALL linsolver%AssemblyB     ! b must be assembled before using
#endif

   END SUBROUTINE ComputeRHS

   SUBROUTINE ComputeRHS_MURK(sem, dt, nelm, linsolver )
      implicit none
      TYPE(DGSem),                intent(in)       :: sem
      REAL(KIND=RP),              intent(in)       :: dt
      integer,                    intent(in)       :: nelm
      CLASS(GenericLinSolver_t),  intent(inout)   :: linsolver

      integer                                      :: Nx, Ny, Nz, l, i, j, k, elmnt, counter   
      REAL(KIND=RP)                                :: value

#if defined(CAHNHILLIARD)
      counter = 0
      DO elmnt = 1, nelm
         Nx = sem%mesh%elements(elmnt)%Nxyz(1)
         Ny = sem%mesh%elements(elmnt)%Nxyz(2)
         Nz = sem%mesh%elements(elmnt)%Nxyz(3)
         DO k = 0, Nz
            DO j = 0, Ny
               DO i = 0, Nx
                  value = sem%mesh%elements(elmnt)%storage%cDot(1,i,j,k)
                  counter =  counter + 1
                  CALL linsolver%SetRHSValue(counter, value)
               END DO
            END DO
         END DO
      END DO
#endif

   END SUBROUTINE ComputeRHS_MURK

   SUBROUTINE ComputeRHS_MUBDF2(sem, dt, nelm, linsolver )
      implicit none
      TYPE(DGSem),                intent(in)       :: sem
      REAL(KIND=RP),              intent(in)       :: dt
      integer,                    intent(in)       :: nelm
      CLASS(GenericLinSolver_t),  intent(inout)   :: linsolver

      integer                                      :: Nx, Ny, Nz, l, i, j, k, elmnt, counter   
      REAL(KIND=RP)                                :: value

#if defined(CAHNHILLIARD)
      counter = 0
      DO elmnt = 1, nelm
         Nx = sem%mesh%elements(elmnt)%Nxyz(1)
         Ny = sem%mesh%elements(elmnt)%Nxyz(2)
         Nz = sem%mesh%elements(elmnt)%Nxyz(3)
         DO k = 0, Nz
            DO j = 0, Ny
               DO i = 0, Nx
                     value = sem%mesh%elements(elmnt)%storage%Q(1,i,j,k) + &
                          dt*sem%mesh%elements(elmnt)%storage%QDot(1,i,j,k)
                     counter =  counter + 1
                     CALL linsolver%SetRHSValue(counter, value)
                     
               END DO
            END DO
         END DO
      END DO

!      CALL linsolver%AssemblyB     ! b must be assembled before using
#endif

   END SUBROUTINE ComputeRHS_MUBDF2

   function IMEX_MatrixShift(dt) result(Ashift)
      implicit none
      !------------------------------------------------------
      real(kind=RP), intent(in) :: dt
      real(kind=RP)             :: Ashift
      !------------------------------------------------------
      
      Ashift = -1.0_RP/dt
      
   end function IMEX_MatrixShift
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
END MODULE IMEXMethods
