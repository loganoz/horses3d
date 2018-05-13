!
!//////////////////////////////////////////////////////
!
!   @File:    IMEXMethods.f90
!   @Author:  Juan (juan.manzanero@upm.es)
!   @Created: Tue Apr 17 16:55:49 2018
!   @Last revision date: Sun May 13 11:22:07 2018
!   @Last revision author: Juan (juan.manzanero@upm.es)
!   @Last revision commit: 664796b96ada01ab3f21660a398ffe36d0c767ef
!
!//////////////////////////////////////////////////////
!
MODULE IMEXMethods
   USE SMConstants                  
   USE DGSEMClass,                  ONLY: DGSem
   USE ElementClass,                ONLY: Element, allocateElementStorage    !arueda: No DGSolutionStorage implemented in nslite3d... Using whole element definitions
   USE PhysicsStorage
   use HexMeshClass
   use MKLPardisoSolverClass
   USE CSRMatrixClass
   USE FTValueDictionaryClass
   use TimeIntegratorDefinitions
   use MatrixClass
   use DGSEMClass, only: ComputeQDot_FCN
   implicit none
   
   PRIVATE                          
   PUBLIC TakeIMEXEulerStep
   
   real(kind=RP) :: time               ! Time at the beginning of each inner(!) time step
   logical       :: computeA = .TRUE.  ! Compute Jacobian? (only valid if it is meant to be computed according to the convergence)
   
CONTAINS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!   
   SUBROUTINE TakeIMEXEulerStep (sem, t , dt , controlVariables, ComputeTimeDerivative, CTD_onlyLinear, CTD_onlyNonLinear)

      IMPLICIT NONE
      TYPE(DGSem),                  INTENT(INOUT)           :: sem                  !<>DGSem class with solution storage 
      REAL(KIND=RP),                INTENT(IN)              :: t                    !< Time at the beginning of time step
      REAL(KIND=RP),                INTENT(IN)              :: dt                   !< Initial (outer) time step (can internally, the subroutine can use a smaller one depending on convergence)
      TYPE(FTValueDictionary),      INTENT(IN)              :: controlVariables     !< Input file variables
      procedure(ComputeQDot_FCN)                            :: ComputeTimeDerivative
      procedure(ComputeQDot_FCN)                            :: CTD_onlyLinear
      procedure(ComputeQDot_FCN)                            :: CTD_onlyNonLinear
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
         DimPrb = sem % NDOF
#else
         DimPrb = sem % NDOF / NTOTALVARS * NCOMP
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
      
   END SUBROUTINE TakeIMEXEulerStep
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
