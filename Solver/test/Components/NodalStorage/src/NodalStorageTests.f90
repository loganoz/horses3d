!
!////////////////////////////////////////////////////////////////////////
!
!      NodalStorageTests.f90
!      Created: May 22, 2015 at 1:03 PM 
!      By: David Kopriva  
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE testNodalStorage  
         USE FTAssertions
         USE NodalStorageClass
         use SMConstants
         use PolynomialInterpAndDerivsModule
         IMPLICIT NONE 
         INTEGER                    :: N
         TYPE(NodalStorage_t)         :: spA
         REAL(KIND=RP), ALLOCATABLE :: f(:), d(:), dExact(:)
         INTEGER                    :: j
         REAL(KIND=RP)              :: s, e
!
!        -----
!        Setup
!        -----
!
         N = 8
         CALL spA % construct(GAUSS,N)
         ALLOCATE(f(0:N))
         ALLOCATE(d(0:N))
         ALLOCATE(dExact(0:N))
         
         DO j = 0, N
            f(j)      = sPa % x(j)**2 
            dExact(j) = 2.0_RP*spA % x(j)
         END DO
!
!        ----------------
!        Integration test
!        ----------------
!
         s = 0.0_RP
         DO j = 0, N
            s = s + f(j)*spA % w(j) 
         END DO
         
         CALL FTAssertEqual(expectedValue = 2.0_RP/3.0_RP,&
                            actualValue = s,              &
                            tol = 1.d-12,                 &
                            msg = "Quadrature of quadratic function")
!
!        --------------------
!        Differentiation test
!        --------------------
!
         CALL MatrixMultiplyDeriv( f, d, spa % DT, N, transp = MXV_DIRECT)
         e = MAXVAL(ABS(d-dExact))
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-12,              &
                            msg           = "Differentiation of quadratic function")
         
      END SUBROUTINE testNodalStorage