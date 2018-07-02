!
!//////////////////////////////////////////////////////
!
!   @File:    IMEXMethods.f90
!   @Author:  Juan (juan.manzanero@upm.es)
!   @Created: Tue Apr 17 16:55:49 2018
!   @Last revision date: Sat Jun 23 10:20:38 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: fce351220409e80ce5df1949249c2b870dd847aa
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
   use DGSEMClass, only: ComputeQDot_FCN
   implicit none
   
   PRIVATE                          
   PUBLIC TakeIMEXStep, SetIMEXComputeQDotProcedures
   
   real(kind=RP) :: time               ! Time at the beginning of each inner(!) time step
   logical       :: computeA = .TRUE.  ! Compute Jacobian? (only valid if it is meant to be computed according to the convergence)

   procedure(ComputeQDot_FCN), protected, pointer    :: CTD_onlyRK3       => NULL()
   procedure(ComputeQDot_FCN), protected, pointer    :: CTD_onlyLinear    => NULL()
   procedure(ComputeQDot_FCN), protected, pointer    :: CTD_onlyNonLinear => NULL()
   
CONTAINS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!   
   subroutine SetIMEXComputeQDotProcedures(CTD_onlyLinear_, CTD_onlyNonLinear_, CTD_onlyRK3_)
      implicit none
      procedure(ComputeQDot_FCN)           :: CTD_onlyLinear_, CTD_onlyNonLinear_
      procedure(ComputeQDot_FCN), optional :: CTD_onlyRK3_

      CTD_onlyLinear    => CTD_onlyLinear_
      CTD_onlyNonLinear => CTD_onlyNonLinear_

      if ( present(CTD_onlyRK3_) ) then
         CTD_onlyRK3       => CTD_onlyRK3_
      end if

   end subroutine SetIMEXComputeQDotProcedures

   subroutine TakeIMEXStep(sem, t, dt, controlVariables, ComputeTimeDerivative)
      implicit none
      TYPE(DGSem),                  INTENT(INOUT)           :: sem                  !<>DGSem class with solution storage 
      REAL(KIND=RP),                INTENT(IN)              :: t                    !< Time at the beginning of time step
      REAL(KIND=RP),                INTENT(IN)              :: dt                   !< Initial (outer) time step (can internally, the subroutine can use a smaller one depending on convergence)
      TYPE(FTValueDictionary),      INTENT(IN)              :: controlVariables     !< Input file variables
      procedure(ComputeQDot_FCN)                            :: ComputeTimeDerivative

      select case(trim(solver))
      case("cahn-hilliard", "nsch")
         call TakeIMEXEulerStep_NSCH (sem, t , dt , controlVariables, ComputeTimeDerivative)

      case("multiphase")
         call TakeIMEXEulerStep_MU (sem, t , dt , controlVariables, ComputeTimeDerivative)

      case("navier-stokes")
         print*, "IMEX solver not implemented for Navier-Stokes"

      end select

   end subroutine TakeIMEXStep

   SUBROUTINE TakeIMEXEulerStep_NSCH (sem, t , dt , controlVariables, ComputeTimeDerivative)

      IMPLICIT NONE
      TYPE(DGSem),                  INTENT(INOUT)           :: sem                  !<>DGSem class with solution storage 
      REAL(KIND=RP),                INTENT(IN)              :: t                    !< Time at the beginning of time step
      REAL(KIND=RP),                INTENT(IN)              :: dt                   !< Initial (outer) time step (can internally, the subroutine can use a smaller one depending on convergence)
      TYPE(FTValueDictionary),      INTENT(IN)              :: controlVariables     !< Input file variables
      procedure(ComputeQDot_FCN)                            :: ComputeTimeDerivative
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

         call linsolver%ComputeAndFactorizeJacobian(nEqnJac,nGradJac, CTD_onlyLinear, dt, 1.0_RP)
         
      ENDIF
      
      time = t

!     TODO do i need this?      
      !CALL sem % GetQ(U_n)      !stores sem%mesh%elements(:)% storage % Q in Vector U_n
!
!     Compute the non linear time derivative
!     -------------------------------------- 
#if defined(CAHNHILLIARD)
      call CTD_onlyNonLinear(sem % mesh, sem % particles, time, sem % BCFunctions)
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
      call sem % SetQ(linsolver % x, NCOMP)
#endif

#if (!defined(NAVIERSTOKES))
!
!     Compute the standard time derivative to get residuals
!     -----------------------------------------------------
      call ComputeTimeDerivative(sem % mesh, sem % particles, time, sem % BCFunctions)

#else
!
!     Perform a RK3 time step in NS
!     -----------------------------
      DO k = 1,3
         tk = t + b(k)*dt
         CALL ComputeTimeDerivative( sem % mesh, sem % particles, tk, sem % BCFunctions)
         
!$omp parallel do schedule(runtime)
         DO id = 1, SIZE( sem % mesh % elements )
             sem % mesh % elements(id) % storage % G_NS = a(k)* sem % mesh % elements(id) % storage % G_NS +          sem % mesh % elements(id) % storage % QDot
             sem % mesh % elements(id) % storage % Q    =       sem % mesh % elements(id) % storage % Q    + c(k)*dt* sem % mesh % elements(id) % storage % G_NS
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
      procedure(ComputeQDot_FCN)                            :: ComputeTimeDerivative
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

#if (defined(NAVIERSTOKES) && defined(CAHNHILLIARD))

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

         call linsolver%ComputeAndFactorizeJacobian(nEqnJac,nGradJac, CTD_onlyLinear, dt, 1.0_RP)
         
      ENDIF
      
      time = t
   
      call ComputeTimeDerivative(sem % mesh, sem % particles, t, sem % BCFunctions)
!
!     Perform a RK3 time step in NS
!     -----------------------------
      DO k = 1,3
         tk = t + b(k)*dt
         CALL CTD_onlyRK3( sem % mesh, sem % particles, tk, sem % BCFunctions)
         
!$omp parallel do schedule(runtime)
         DO id = 1, SIZE( sem % mesh % elements )
             sem % mesh % elements(id) % storage % G_NS = a(k)* sem % mesh % elements(id) % storage % G_NS +          sem % mesh % elements(id) % storage % QDot
             sem % mesh % elements(id) % storage % Q    =       sem % mesh % elements(id) % storage % Q    + c(k)*dt* sem % mesh % elements(id) % storage % G_NS
         END DO
!$omp end parallel do
      END DO
!
!     Compute the non linear time derivative
!     -------------------------------------- 
      call CTD_onlyNonLinear(sem % mesh, sem % particles, time, sem % BCFunctions)
!
!     Change the concentration to its updated value from the density
!     --------------------------------------------------------------
!$omp parallel do schedule(runtime)
      do id = 1, size(sem % mesh % elements)
         sem % mesh % elements(id) % storage % c(1,:,:,:) = (sem % mesh % elements(id) % storage % QNS(IRHO,:,:,:) - multiphase % barRho)/multiphase % tildeRho
      end do
!$omp end parallel do
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
      call sem % SetQ(linsolver % x, NCOMP)
!
!     Return NS as main storage
!     -------------------------
      call sem % mesh % SetStorageToEqn(1)
!
!     Update the density with the new concentration value
!     ---------------------------------------------------
!$omp parallel do schedule(runtime)
      do id = 1, size(sem % mesh % elements)
         sem % mesh % elements(id) % storage % Q(IRHO,:,:,:) = multiphase % tildeRho * sem % mesh % elements(id) % storage % c(1,:,:,:) + multiphase % barRho
      end do
!$omp end parallel do

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
