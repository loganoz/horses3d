!
!////////////////////////////////////////////////////////////////////////
!
!   @File:    LegendreAlgorithms.f90
!   @Author:  David Kopriva
!   @Created: Created: 2009-12-08 13:47:42 -0500 
!   @Last revision date: Wed Jun  2 18:14:23 2021
!   @Last revision author: Wojciech Laskowski (wj.laskowski@upm.es)
!   @Last revision commit: 28e347f1c2e07ee09992b923a6ed72dacfe9ffb9
!
!      Contains:
!            ALGORITHM 20: FUNCTION LegendrePolynomial( N, x ) RESULT(L_N)
!            ALGORITHM 22: SUBROUTINE LegendrePolyAndDerivative( N, x, L_N, LPrime_N )
!            ALGORITHM 23: SUBROUTINE GaussLegendreNodesAndWeights( N, x, w )
!            ALGORITHM 24: SUBROUTINE qAndLEvaluation( N, x, Q, Q_prime, L_N )
!            ALGORITHM 25: SUBROUTINE LegendreLobattoNodesAndWeights( N, x, w )
!
!////////////////////////////////////////////////////////////////////////
!
      FUNCTION LegendrePolynomial( N, x ) RESULT(L_N)
!
!     ----------------------------------------------------------------------
!     Compute the  Legendre Polynomial by the three point recursion
!
!           L (x) = 2k-1 xL     - k-1 L
!            k      ----   k-1    ---  k-2
!                     k            k
!
!     ------------------------------------------------------------------------
!
      USE SMConstants
      IMPLICIT NONE
!
!     Compute the Legendre Polynomial of degree k and its derivative
!
!     -----------------
!     Input parameters:
!     -----------------
!
      INTEGER      , INTENT(IN) :: N
      REAL(KIND=RP), INTENT(IN) :: x
!
!     ------------------
!     Output parameters:
!     ------------------
!
      REAL(KIND=RP) :: L_N
!
!     ----------------
!     Local Variables:
!     ----------------
!
      INTEGER       :: k
      REAL(KIND=RP) :: L_NM1, L_NM2
      
      IF( N == 0 )     THEN
         L_N      = 1.0_rp
      ELSE IF ( N == 1 )     THEN
         L_N      = x
      ELSE
         L_NM2 = 1.0_rp
         L_NM1 = x
         DO k = 2, N
            L_N        = ((2*k-1)*x*L_NM1 - (k-1)*L_NM2)/k
            L_NM2      = L_NM1
            L_NM1      = L_N
         END DO
      END IF
      
      END FUNCTION LegendrePolynomial
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      MODULE GaussQuadrature
      USE SMConstants
      IMPLICIT NONE
      
      PUBLIC  :: GaussLegendreNodesAndWeights, LegendreLobattoNodesAndWeights
      PUBLIC  :: LegendrePolyAndDerivative
      PRIVATE :: qAndLEvaluation
!
!     ========
      CONTAINS
!     ========
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE LegendrePolyAndDerivative( N, x, L_N, LPrime_N )
!
!     Compute the Legendre Polynomial of degree k and its derivative
!
!     -----------------
!     Input parameters:
!     -----------------
!
      INTEGER      , INTENT(IN) :: N
      REAL(KIND=RP), INTENT(IN) :: x
!
!     ------------------
!     Output parameters:
!     ------------------
!
      REAL(KIND=RP), INTENT(OUT) :: L_N, LPrime_N
!
!     ----------------
!     Local Variables:
!     ----------------
!
      INTEGER       :: k
      REAL(KIND=RP) :: L_NM1, L_NM2, LPrime_NM2, LPrime_NM1
      
      IF( N == 0 )     THEN
         L_N      = 1.0_rp
         LPrime_N = 0.0_rp
      ELSE IF ( N == 1 )     THEN
         L_N      = x
         LPrime_N = 1.0_rp
      ELSE
         L_NM2 = 1.0_rp
         LPrime_NM2 = 0.0_rp
         L_NM1 = x
         LPrime_NM1 = 1.0_rp
         DO k = 2, N
            L_N        = ((2*k-1)*x*L_NM1 - (k-1)*L_NM2)/k
            LPrime_N   = LPrime_NM2 + (2*k-1)*L_NM1
            L_NM2      = L_NM1
            L_NM1      = L_N
            LPrime_NM2 = LPrime_NM1
            LPrime_NM1 = LPrime_N
         END DO
      END IF
      
      END SUBROUTINE LegendrePolyAndDerivative
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE GaussLegendreNodesAndWeights( N, x, w )
!
!     Compute the Gauss-legendre quadrature nodes and
!     weights
!
!     -----------------
!     Input parameters:
!     -----------------
!
      INTEGER, INTENT(IN) :: N
!
!     ------------------
!     Output parameters:
!     ------------------
!
      REAL(KIND=RP), DIMENSION(0:N), INTENT(OUT) :: x, w
!
!     ----------------
!     Local Variables:
!     ----------------
!
      REAL(KIND=RP) :: xj, L_NP1, LPrime_NP1, delta, tolerance
      INTEGER       :: j, k, NDiv2
!
!     ----------
!     Constants:
!     ----------
!
      INTEGER, PARAMETER       :: noNewtonIterations = 10
      REAL(KIND=RP), PARAMETER :: toleranceFactor    = 4.0_RP
      
      tolerance = toleranceFactor*EPSILON(L_NP1)
      IF( N == 0 )     THEN
         x(0) = 0.0_RP
         w(0) = 2.0_RP
         RETURN
      ELSE IF( N == 1 )     THEN
         x(0) = -SQRT(1.0_RP/3.0_RP)
         w(0) =  1.0_RP
         x(1) = -x(0)
         w(1) =  w(0)
      ELSE
!
!        ----------------------------------
!        Iterate on half the interior nodes
!        ----------------------------------
!
         NDiv2 = (N+1)/2
         DO j = 0, NDiv2-1
            xj = -COS( (2*j+1)*PI/(2*N+2) )
            DO k = 0, noNewtonIterations
               CALL LegendrePolyAndDerivative( N+1, xj, L_NP1, LPrime_NP1 )
               delta = -L_NP1/LPrime_NP1
               xj = xj + delta
               IF( ABS(delta) <=  tolerance*ABS(xj) )     EXIT
            END DO
            CALL LegendrePolyAndDerivative( N+1, xj, L_NP1, LPrime_NP1 )
            x(j)   = xj
            w(j)   = 2.0_RP/( (1.0_RP - xj**2)*LPrime_NP1**2 )
            x(N-j) = -xj
            w(N-j) = w(j)
         END DO
      END IF
!
!     ---------------------------
!     Fill in middle if necessary
!     ---------------------------
!
      IF( MOD(N,2) == 0 )     THEN
         CALL LegendrePolyAndDerivative( N+1, 0.0_RP, L_NP1, LPrime_NP1 )
         x(N/2) = 0.0_RP
         w(N/2) = 2.0_RP/LPrime_NP1**2
      END IF
      
      END SUBROUTINE GaussLegendreNodesAndWeights
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE qAndLEvaluation( N, x, Q, Q_prime, L_N )
!
!     Compute the function Q_N = L_{N+1} - L_{N-1} and
!     its derivative to find the Gauss-Lobatto points and
!     weights.
!
!     -----------------
!     Input parameters:
!     -----------------
!
      INTEGER      , INTENT(IN) :: N
      REAL(KIND=RP), INTENT(IN) :: x
!
!     ------------------
!     Output parameters:
!     ------------------
!
      REAL(KIND=RP), INTENT(OUT) :: Q, Q_prime, L_N
!
!     ----------------
!     Local Variables:
!     ----------------
!
      INTEGER       :: k
      REAL(KIND=RP) :: L_kM1, L_kM2, LPrime_kM2, LPrime_kM1, L_k, LPrime_k
      
      IF( N == 0 )     THEN       !Should never be called
         L_N      = 1.0_rp
         Q = HUGE(Q)
         Q_Prime = Huge(Q_Prime)
      ELSE IF ( N == 1 )     THEN ! Should never be called
         L_N      = x
         Q = HUGE(Q)
         Q_Prime = Huge(Q_Prime)
      ELSE
         L_kM2 = 1.0_rp
         LPrime_kM2 = 0.0_rp
         L_kM1 = x
         LPrime_kM1 = 1.0_rp
         DO k = 2, N
            L_k        = ((2*k-1)*x*L_kM1 - (k-1)*L_kM2)/k
            LPrime_k   = LPrime_kM2 + (2*k-1)*L_kM1
            L_kM2      = L_kM1
            L_kM1      = L_k
            LPrime_kM2 = LPrime_kM1
            LPrime_kM1 = LPrime_k
         END DO
         k = N+1
         L_k        = ((2*k-1)*x*L_kM1 - (k-1)*L_kM2)/k
         LPrime_k   = LPrime_kM2 + (2*k-1)*L_kM1
         
         Q = L_k - L_kM2
         Q_prime = LPrime_k - LPrime_kM2
         L_N = L_kM1
      END IF
      
      END SUBROUTINE qAndLEvaluation
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE LegendreLobattoNodesAndWeights( N, x, w )
!
!     Compute the Gauss-Lobatto-legendre quadrature nodes and
!     weights
!
!     -----------------
!     Input parameters:
!     -----------------
!
      INTEGER, INTENT(IN) :: N
!
!     ------------------
!     Output parameters:
!     ------------------
!
      REAL(KIND=RP), DIMENSION(0:N), INTENT(OUT) :: x, w
!
!     ----------------
!     Local Variables:
!     ----------------
!
      REAL(KIND=RP) :: xj, Q, Q_prime, L_N, delta, tolerance, LPrime_NP1
      INTEGER       :: j, k, NDiv2
!
!     ----------
!     Constants:
!     ----------
!
      INTEGER, PARAMETER       :: noNewtonIterations = 10
      REAL(KIND=RP), PARAMETER :: toleranceFactor    = 4.0_RP
      
      tolerance = toleranceFactor*EPSILON(L_N)
      IF( N == 0 )     THEN ! Error - Must have at least two points, +/-1
         x(0) = 0.0_RP
         w(0) = 0.0_RP
         RETURN
      ELSE IF( N == 1 )     THEN
         x(0) = -1.0_RP
         w(0) =  1.0_RP
         x(1) =  1.0_RP
         w(1) =  w(0)
      ELSE
!
!        ---------
!        Endpoints
!        ---------
!
         x(0) = -1.0_RP
         w(0) =  2.0_RP/(N*(N+1))
         x(N) =  1.0_RP
         w(N) =  w(0)
!
!        ----------------------------------
!        Iterate on half the interior nodes
!        ----------------------------------
!
         NDiv2 = (N+1)/2
         DO j = 1, NDiv2-1
            xj = -COS( (j+0.25_RP)*PI/N - 3.0_RP/(8*N*PI*(j+0.25_RP)) )
            DO k = 0, noNewtonIterations
               CALL qAndLEvaluation( N, xj, Q, Q_prime, L_N )
               delta = -Q/Q_prime
               xj = xj + delta
               IF( ABS(delta) <=  tolerance*ABS(xj) )     EXIT
            END DO
            CALL qAndLEvaluation( N, xj, Q, Q_prime, L_N )
            x(j)   = xj
            w(j)   = 2.0_RP/(N*(N+1)*L_N**2)
            x(N-j) = -xj
            w(N-j) = w(j)
         END DO
      END IF
!
!     ---------------------------
!     Fill in middle if necessary
!     ---------------------------
!
      IF( MOD(N,2) == 0 )     THEN
         CALL LegendrePolyAndDerivative( N, 0.0_RP, L_N, LPrime_NP1 )
         x(N/2) = 0.0_RP
         w(N/2) = 2.0_RP/(N*(N+1)*L_N**2)
      END IF
      
      END SUBROUTINE LegendreLobattoNodesAndWeights
!
!     ==========      
      END MODULE GaussQuadrature
!     ==========      