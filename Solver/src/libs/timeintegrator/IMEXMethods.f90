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
   PUBLIC TakeIMEXStep, Enable_CTD_AFTER_STEPS_IMEX
   
   real(kind=RP) :: time               ! Time at the beginning of each inner(!) time step
   logical       :: computeA = .TRUE.  ! Compute Jacobian? (only valid if it is meant to be computed according to the convergence)

   logical, protected :: CTD_AFTER_STEPS = .false.
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
         error stop

      case(CAHNHILLIARD_SOLVER)
         call TakeIMEXEulerStep_CH(sem, t , dt , controlVariables, ComputeTimeDerivative)

      case(MULTIPHASE_SOLVER)
         call TakeIMEXBDF2Step_MU(sem, t , dt , controlVariables, ComputeTimeDerivative)

      case default
         print*, "Solver not recognized"
         error stop
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

         call linsolver%ComputeAndFactorizeJacobian(nEqnJac,nGradJac, ComputeTimeDerivative, dt, 1.0_RP, CTD_IMEX_IMPLICIT)
         
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
      type(MKLPardisoSolver_t), save           :: linsolver_Lap
      integer                                  :: nelm, DimPrb, globalDimPrb
      logical, save                            :: isfirst = .TRUE.
      logical, save                            :: isSecond = .false.
      logical, save                            :: isSecondOrder = .true.
      logical                                  :: TimeAccurate = .true.
      integer                                  :: nEqnJac, nGradJac
      real(kind=RP)                            :: tk
      integer                                  :: k, id, s, i, j, counter
      integer, save                            :: IMEX_order = 2
      logical, save                            :: ACM2_MODEL = .false.
      SAVE DimPrb, nelm, TimeAccurate

      tk = t + dt

#if defined(MULTIPHASE) && defined(CAHNHILLIARD)
      nEqnJac = NCOMP
      nGradJac = NCOMP

      if (isSecond .and. isSecondOrder) then
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
!        Compute the t derivative (Explicit)
!        --------------------------------------
         call ComputeTimeDerivative(sem % mesh, sem % particles, tk, CTD_IMEX_EXPLICIT)
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
!        Perform the implicit t-step
!        ------------------------------
         call ComputeRHS_MUBDF2(sem, dt, nelm, linsolver)
         call linsolver % SolveLUDirect
         call sem % mesh % SetStorageToEqn(C_BC)
         sem % mesh % storage % Q = linsolver % x
         call sem % mesh % storage % global2LocalQ
         call sem % mesh % SetStorageToEqn(NS_BC)

         if (ACM2_MODEL) then
!
!           solve the pressure laplacian
!           ----------------------------
            call ComputeRHS_Lap(sem, dt, nelm, linsolver_Lap)
            call linsolver_Lap % solveLUdirect
            counter = 0
            do id = 1, sem % mesh % no_of_elements
               associate(e => sem % mesh % elements(id))
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  counter = counter + 1
                  e % storage % QDot(IMP,i,j,k) = linsolver_Lap % x(counter)
               end do                ; end do                ; end do
               end associate
            end do
         end if

!
!$omp parallel do schedule(runtime)
         do id = 1, sem % mesh % no_of_elements
!
!           Set the concentration from the linear solver solution
!           -----------------------------------------------------
            sem % mesh % elements(id) % storage % Q(IMC,:,:,:) = sem % mesh % elements(id) % storage % c(1,:,:,:)
!
!           Perform a t-step on the rest of the variables
!           ------------------------------------------------
            sem % mesh % elements(id) % storage % Q(IMC+1:,:,:,:) = (2.0_RP / 3.0_RP ) * (sem % mesh % elements(id) % storage % Q(IMC+1:,:,:,:) & 
                                                                    + dt*sem % mesh % elements(id) % storage % QDot(IMC+1:,:,:,:))
         end do 
!$omp end parallel do

      end if

      IF (isfirst) THEN           
         if ( controlVariables % containsKey("imex order") ) then
            IMEX_order = controlVariables % IntegerValueForKey("imex order")
         else
            IMEX_order = 2
         end if
         print*, "IMEX order = ", IMEX_order

         select case (IMEX_order)
         case(1)
            isSecondOrder = .false.
         case(2)
            isSecondOrder = .true.
         case default
            print*, "IMEX order should be 1 or 2"
            errorMessage(STD_OUT)
            error stop
         end select

!
!        ***********************************************************************
!           Construct the Jacobian, and perform the factorization in the first
!           call.
!        ***********************************************************************


         isSecond = .true.
         nelm = SIZE(sem%mesh%elements)
         DimPrb = sem % NDOF * NCOMP
         globalDimPrb = sem % totalNDOF * NCOMP
!
!        Construct the CahnHilliard Jacobian         
         CALL linsolver%construct(DimPrb,globalDimPrb, nEqnJac,controlVariables,sem, IMEX_MatrixShift) 
         call linsolver%ComputeAndFactorizeJacobian(nEqnJac,nGradJac, ComputeTimeDerivative, dt, 1.0_RP, CTD_IMEX_IMPLICIT)

         if (ACM2_MODEL) then
!
!           Construct the pressure laplacian Jacobian
            CALL linsolver_Lap%construct(DimPrb,globalDimPrb, nEqnJac,controlVariables,sem, Laplacian_MatrixShift) 
            call linsolver_Lap%ComputeAndFactorizeJacobian(nEqnJac,nGradJac, ComputeTimeDerivative, -1.0_RP, 1.0_RP, CTD_LAPLACIAN)
         end if
      end if

      if ( isFirst .or. (.not. isSecondOrder) ) then
         ACM2_MODEL = controlvariables % LogicalValueForKey("enable acm model 2")
         isfirst = .false.
         call ComputeTimeDerivative(sem % mesh, sem % particles, tk, CTD_IMEX_EXPLICIT)
!
!        Solve the implicit part of the CHE   
!        ----------------------------------
         call ComputeRHS_MUBDF2(sem, dt, nelm, linsolver)
         call linsolver % SolveLUDirect
         call sem % mesh % SetStorageToEqn(C_BC)
         sem % mesh % storage % Q = linsolver % x
         call sem % mesh % storage % global2LocalQ
         call sem % mesh % SetStorageToEqn(NS_BC)

         if (ACM2_MODEL) then
!
!           solve the pressure laplacian
!           ----------------------------
            call ComputeRHS_Lap(sem, dt, nelm, linsolver_Lap)
            call linsolver_Lap % solveLUdirect
            counter = 0
            do id = 1, sem % mesh % no_of_elements
               associate(e => sem % mesh % elements(id))
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  counter = counter + 1
                  e % storage % QDot(IMP,i,j,k) = linsolver_Lap % x(counter)
               end do                ; end do                ; end do
               end associate
            end do
         end if

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
!           Perform a t-step on the rest of the variables
!           ------------------------------------------------
            sem % mesh % elements(id) % storage % Q(IMC+1:,:,:,:) = sem % mesh % elements(id) % storage % Q(IMC+1:,:,:,:) & 
                                                                    + dt*sem % mesh % elements(id) % storage % QDot(IMC+1:,:,:,:)

         end do 
!$omp end parallel do
!
!        Change the coefficient in the Jacobian for the BDF2
!        ---------------------------------------------------
         if ( isSecondOrder ) then
            call linsolver % ReFactorizeJacobian
            print*, "Jacobian re-factorized"
         end if
         
      ENDIF
!
!     The "good" residuals CTD call
!     -----------------------------
      if ( CTD_AFTER_STEPS) call ComputeTimeDerivative(sem % mesh, sem % particles, tk, CTD_IGNORE_MODE)

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

   SUBROUTINE ComputeRHS_Lap(sem, dt, nelm, linsolver )
      implicit none
      TYPE(DGSem),                intent(in)       :: sem
      REAL(KIND=RP),              intent(in)       :: dt
      integer,                    intent(in)       :: nelm
      CLASS(GenericLinSolver_t),  intent(inout)   :: linsolver

      integer                                      :: Nx, Ny, Nz, l, i, j, k, elmnt, counter   
      REAL(KIND=RP)                                :: value

#if defined(CAHNHILLIARD) && defined(MULTIPHASE)
      counter = 0
      DO elmnt = 1, nelm
         Nx = sem%mesh%elements(elmnt)%Nxyz(1)
         Ny = sem%mesh%elements(elmnt)%Nxyz(2)
         Nz = sem%mesh%elements(elmnt)%Nxyz(3)
         DO k = 0, Nz
            DO j = 0, Ny
               DO i = 0, Nx
                     value = -sem%mesh%elements(elmnt)%storage%QDot(IMP,i,j,k) 
                     counter =  counter + 1
                     CALL linsolver%SetRHSValue(counter, value)
               END DO
            END DO
         END DO
      END DO

!      CALL linsolver%AssemblyB     ! b must be assembled before using
#endif

   END SUBROUTINE ComputeRHS_Lap

   function IMEX_MatrixShift(dt) result(Ashift)
      implicit none
      !------------------------------------------------------
      real(kind=RP), intent(in) :: dt
      real(kind=RP)             :: Ashift
      !------------------------------------------------------
      
      Ashift = -1.0_RP/dt
      
   end function IMEX_MatrixShift

   function Laplacian_MatrixShift(dt) result(Ashift)
      implicit none
      !------------------------------------------------------
      real(kind=RP), intent(in) :: dt
      real(kind=RP)             :: Ashift
      !------------------------------------------------------
      
      Ashift = 0.0_RP
      
   end function Laplacian_MatrixShift

   subroutine Enable_CTD_AFTER_STEPS_IMEX()
      implicit none
      CTD_AFTER_STEPS = .true.
   end subroutine Enable_CTD_AFTER_STEPS_IMEX


!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
END MODULE IMEXMethods