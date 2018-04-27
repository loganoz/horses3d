!
!//////////////////////////////////////////////////////
!
!   @File:    IMEXMethods.f90
!   @Author:  Juan (juan.manzanero@upm.es)
!   @Created: Tue Apr 17 16:55:49 2018
!   @Last revision date: Fri Apr 27 12:22:05 2018
!   @Last revision author: Juan (juan.manzanero@upm.es)
!   @Last revision commit: c3532365f3cc0c1e6e95281cbe9836354994daea
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
         DimPrb = sem % NDOF
         
         ALLOCATE(U_n(0:Dimprb-1))

         CALL linsolver%construct(DimPrb,controlVariables,sem) 

         call linsolver%ComputeAndFactorizeJacobian(nEqnJac,nGradJac, CTD_onlyLinear, dt, 1.0_RP)
         
      ENDIF
      
      time = t

!     TODO do i need this?      
      !CALL sem % GetQ(U_n)      !stores sem%mesh%elements(:)% storage % Q in Vector U_n
!
!     Compute the non linear time derivative
!     -------------------------------------- 
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
      call sem % SetQ(linsolver % x)
!
!     Compute the standard time derivative to get residuals
!     -----------------------------------------------------
      call ComputeTimeDerivative(sem % mesh, sem % particles, time, sem % BCFunctions)
      
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

                     CALL linsolver%SetBValue(counter, value)
                     counter =  counter + 1
                  END DO
               END DO
            END DO
         END DO
      END DO

!      CALL linsolver%AssemblyB     ! b must be assembled before using
#endif

   END SUBROUTINE ComputeRHS
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
END MODULE IMEXMethods
