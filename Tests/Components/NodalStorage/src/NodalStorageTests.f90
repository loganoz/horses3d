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
         IMPLICIT NONE 
         INTEGER                    :: N
         TYPE(NodalStorage)         :: spA
         REAL(KIND=RP), ALLOCATABLE :: f(:), d(:), dExact(:)
         INTEGER                    :: j
         REAL(KIND=RP)              :: s, e
!
!        -----
!        Setup
!        -----
!
         N = 8
         CALL spA % construct(N)
         ALLOCATE(f(0:N))
         ALLOCATE(d(0:N))
         ALLOCATE(dExact(0:N))
         
         DO j = 0, N
            f(j)      = sPa % xi(j)**2 
            dExact(j) = 2.0_RP*spA % xi(j)
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
         CALL MatrixMultiplyDeriv( f, d, spa % standardDerivativeMatrix, N, transp = MXV_DIRECT)
         e = MAXVAL(ABS(d-dExact))
         
         CALL FTAssertEqual(expectedValue = 0.0_RP,              &
                            actualValue   = e,                   &
                            tol           = 1.d-12,              &
                            msg           = "Differentiation of quadratic function")
         
      END SUBROUTINE testNodalStorage
