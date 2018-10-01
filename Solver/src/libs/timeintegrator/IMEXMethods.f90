!
!//////////////////////////////////////////////////////
!
!   @File:    IMEXMethods.f90
!   @Author:  Juan (juan.manzanero@upm.es)
!   @Created: Tue Apr 17 16:55:49 2018
!   @Last revision date: Mon Aug 20 17:10:12 2018
!   @Last revision author: Andrés Rueda (am.rueda@upm.es)
!   @Last revision commit: 9fb80d209ec1b9ae1b044040a2af4e790b2ecd64
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
      TYPE(DGSem),                  INTENT(INOUT)           :: sem                  !<>DGSem class with solution storage 
      REAL(KIND=RP),                INTENT(IN)              :: t                    !< Time at the beginning of time step
      REAL(KIND=RP),                INTENT(IN)              :: dt                   !< Initial (outer) time step (can internally, the subroutine can use a smaller one depending on convergence)
      TYPE(FTValueDictionary),      INTENT(IN)              :: controlVariables     !< Input file variables
      procedure(ComputeTimeDerivative_f)                            :: ComputeTimeDerivative

      select case(solver)
      case(NAVIERSTOKES_SOLVER, INCNS_SOLVER)
         print*, "IMEX solver not implemented for monophase Navier-Stokes equations"
         stop

      case(CAHNHILLIARD_SOLVER, NSCH_SOLVER)
         call TakeIMEXEulerStep_NSCH (sem, t , dt , controlVariables, ComputeTimeDerivative)

      case(INSCH_SOLVER)
         call TakeIMEXEulerStep_MU (sem, t , dt , controlVariables, ComputeTimeDerivative)

      case default
         print*, "Solver not recognized"
         stop
      end select

   end subroutine TakeIMEXStep

   SUBROUTINE TakeIMEXEulerStep_NSCH (sem, t , dt , controlVariables, ComputeTimeDerivative)

      IMPLICIT NONE
      TYPE(DGSem),                  INTENT(INOUT)           :: sem                  !<>DGSem class with solution storage 
      REAL(KIND=RP),                INTENT(IN)              :: t                    !< Time at the beginning of time step
      REAL(KIND=RP),                INTENT(IN)              :: dt                   !< Initial (outer) time step (can internally, the subroutine can use a smaller one depending on convergence)
      TYPE(FTValueDictionary),      INTENT(IN)              :: controlVariables     !< Input file variables
      procedure(ComputeTimeDerivative_f)                            :: ComputeTimeDerivative
!
!     ---------------
!     Local variables
!     ---------------
!
      type(MKLPardisoSolver_t), save           :: linsolver
      INTEGER                                  :: nelm, DimPrb
      LOGICAL, save                            :: isfirst = .TRUE.
      REAL(KIND=RP), DIMENSION(:), ALLOCATABLE :: U_n                                   !Solution at the beginning of time step (even for inner time steps)
      LOGICAL                                  :: TimeAccurate = .true.
      integer                                  :: nEqnJac, nGradJac
      REAL(KIND=RP), DIMENSION(3)              :: a = (/0.0_RP       , -5.0_RP /9.0_RP , -153.0_RP/128.0_RP/)
      REAL(KIND=RP), DIMENSION(3)              :: b = (/0.0_RP       ,  1.0_RP /3.0_RP ,    3.0_RP/4.0_RP  /)
      REAL(KIND=RP), DIMENSION(3)              :: c = (/1.0_RP/3.0_RP,  15.0_RP/16.0_RP,    8.0_RP/15.0_RP /)
      real(kind=RP)                            :: tk
      INTEGER                                  :: k, id
      SAVE DimPrb, nelm, TimeAccurate

#if defined(CAHNHILLIARD)
      nEqnJac = NCOMP
      nGradJac = NCOMP
#endif
      
      IF (isfirst) THEN           
!
!        ***********************************************************************
!           Construct the Jacobian, and perform the factorization in the first
!           call.
!        ***********************************************************************
!
         isfirst = .FALSE.
         nelm = SIZE(sem%mesh%elements)
#if (!defined(CAHNHILLIARD))
         print*, "IMEX Methods only configured to solve Cahn-Hilliard"
         stop
#else
         DimPrb = sem % NDOF * NCOMP
#endif
         
         ALLOCATE(U_n(0:Dimprb-1))

         CALL linsolver%construct(DimPrb,controlVariables,sem, IMEXEuler_MatrixShift) 

         call linsolver%ComputeAndFactorizeJacobian(nEqnJac,nGradJac, ComputeTimeDerivative, dt, 1.0_RP)
         
      ENDIF
      
      time = t

!     TODO do i need this?      
      !CALL sem % GetQ(U_n)      !stores sem%mesh%elements(:)% storage % Q in Vector U_n
!
!     Compute the non linear time derivative
!     -------------------------------------- 
#if defined(CAHNHILLIARD)
      call ComputeTimeDerivative(sem % mesh, sem % particles, time, CTD_ONLY_CH_NONLIN)
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
#endif

#if (!defined(NAVIERSTOKES))
!
!     Compute the standard time derivative to get residuals
!     -----------------------------------------------------
      call ComputeTimeDerivative(sem % mesh, sem % particles, time, CTD_IGNORE_MODE)

#else
!
!     *****************************
!     Perform a RK3 time step in NS
!     *****************************
!
!     Compute the new chemical potential
!     ----------------------------------
#if defined(CAHNHILLIARD)
      CALL ComputeTimeDerivative( sem % mesh, sem % particles, tk, CTD_ONLY_CH)
#endif

      DO k = 1,3
         tk = t + b(k)*dt
         CALL ComputeTimeDerivative( sem % mesh, sem % particles, tk, CTD_ONLY_NS)
         
!$omp parallel do schedule(runtime)
         DO id = 1, SIZE( sem % mesh % elements )
             sem % mesh % elements(id) % storage % G_NS = a(k)* sem % mesh % elements(id) % storage % G_NS +          sem % mesh % elements(id) % storage % QDot
!             sem % mesh % elements(id) % storage % Q    =       sem % mesh % elements(id) % storage % Q    + c(k)*dt* sem % mesh % elements(id) % storage % G_NS
         END DO
!$omp end parallel do
      END DO
#endif
      
   END SUBROUTINE TakeIMEXEulerStep_NSCH

   SUBROUTINE TakeIMEXEulerStep_MU (sem, t , dt , controlVariables, ComputeTimeDerivative)

      IMPLICIT NONE
      TYPE(DGSem),                  INTENT(INOUT)           :: sem                  !<>DGSem class with solution storage 
      REAL(KIND=RP),                INTENT(IN)              :: t                    !< Time at the beginning of time step
      REAL(KIND=RP),                INTENT(IN)              :: dt                   !< Initial (outer) time step (can internally, the subroutine can use a smaller one depending on convergence)
      TYPE(FTValueDictionary),      INTENT(IN)              :: controlVariables     !< Input file variables
      procedure(ComputeTimeDerivative_f)                            :: ComputeTimeDerivative
!
!     ---------------
!     Local variables
!     ---------------
!
      type(MKLPardisoSolver_t), save           :: linsolver
      INTEGER                                  :: nelm, DimPrb
      LOGICAL, save                            :: isfirst = .TRUE.
      LOGICAL                                  :: TimeAccurate = .true.
      integer                                  :: nEqnJac, nGradJac
      REAL(KIND=RP), DIMENSION(3)              :: a = (/0.0_RP       , -5.0_RP /9.0_RP , -153.0_RP/128.0_RP/)
      REAL(KIND=RP), DIMENSION(3)              :: b = (/0.0_RP       ,  1.0_RP /3.0_RP ,    3.0_RP/4.0_RP  /)
      REAL(KIND=RP), DIMENSION(3)              :: c = (/1.0_RP/3.0_RP,  15.0_RP/16.0_RP,    8.0_RP/15.0_RP /)
      real(kind=RP)                            :: tk
      INTEGER                                  :: k, id
      SAVE DimPrb, nelm, TimeAccurate

#if (defined(INCNS) && defined(CAHNHILLIARD))

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
         
         CALL linsolver%construct(DimPrb,controlVariables,sem, IMEXEuler_MatrixShift) 

         call linsolver%ComputeAndFactorizeJacobian(nEqnJac,nGradJac, ComputeTimeDerivative, dt, 1.0_RP)
         
      ENDIF
      
      time = t
!
!     *************************************************
!     1) Compute cDot (full) to obtain \nabla c, and mu
!     *************************************************
!
      call ComputeTimeDerivative(sem % mesh, sem % particles, t, CTD_ONLY_CH)
!
!     ********************************
!     2) Perform a RK3 time step in NS
!     ********************************
!
      DO k = 1,3
         tk = t + b(k)*dt
         CALL ComputeTimeDerivative( sem % mesh, sem % particles, tk, CTD_ONLY_NS)
         
!$omp parallel do schedule(runtime)
         DO id = 1, SIZE( sem % mesh % elements )
             sem % mesh % elements(id) % storage % G_NS = a(k)* sem % mesh % elements(id) % storage % G_NS +          sem % mesh % elements(id) % storage % QDot
             sem % mesh % elements(id) % storage % QNS  =       sem % mesh % elements(id) % storage % QNS + c(k)*dt* sem % mesh % elements(id) % storage % G_NS
         END DO
!$omp end parallel do
      END DO
!
!     ******************************************************************
!     3) Compatibilize the concentration with the new (advected) density
!     ******************************************************************
!
      call sem % mesh % ConvertDensityToPhaseField
!
!     *****************************************
!     4) Compute Cahn-Hilliard non-linear terms
!     *****************************************
!
      call sem % mesh % SetStorageToEqn(2)
      call ComputeTimeDerivative(sem % mesh, sem % particles, time, CTD_ONLY_CH_NONLIN)
!
!     Compute the RHS
!     ---------------
      call ComputeRHS(sem, dt, nelm, linsolver)
!
!     ************************************
!     5) Solve Cahn-Hilliard linear system
!     ************************************
!
      call linsolver % SolveLUDirect
!
!     Return the computed state vector to storage
!     -------------------------------------------
      sem % mesh % storage % Q = linsolver % x
      call sem % mesh % storage % global2LocalQ
!
!     Return NS as main storage
!     -------------------------
      call sem % mesh % SetStorageToEqn(1)
!
!     ******************************************************
!     6) Update the density with the new concentration value
!     ******************************************************
!
      call sem % mesh % ConvertPhaseFieldToDensity

#else
      print*, "Multiphase IMEX solver only works with both Navier-Stokes and Cahn-Hilliard equations."
      errorMessage(STD_OUT)
      stop

#endif
      
   END SUBROUTINE TakeIMEXEulerStep_MU

!  
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE ComputeRHS(sem, dt, nelm, linsolver )
      implicit none
      TYPE(DGSem),                INTENT(IN)       :: sem
      REAL(KIND=RP),              INTENT(IN)       :: dt
      INTEGER,                    INTENT(IN)       :: nelm
      CLASS(GenericLinSolver_t),  INTENT (INOUT)   :: linsolver

      INTEGER                                      :: Nx, Ny, Nz, l, i, j, k, elmnt, counter   
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

                     CALL linsolver%SetRHSValue(counter, value)
                     counter =  counter + 1
                  END DO
               END DO
            END DO
         END DO
      END DO

!      CALL linsolver%AssemblyB     ! b must be assembled before using
#endif

   END SUBROUTINE ComputeRHS

   function IMEXEuler_MatrixShift(dt) result(Ashift)
      implicit none
      !------------------------------------------------------
      real(kind=RP), intent(in) :: dt
      real(kind=RP)             :: Ashift
      !------------------------------------------------------
      
      Ashift = -1.0_RP/dt
      
   end function IMEXEuler_MatrixShift
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
END MODULE IMEXMethods
