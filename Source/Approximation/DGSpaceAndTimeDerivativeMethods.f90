!
!////////////////////////////////////////////////////////////////////////
!
!      DGTimeDerivativeRoutines.f95
!      Created: 2008-07-13 16:13:12 -0400 
!      By: David Kopriva  
!
!      3D version by D.A. Kopriva 6/17/15, 12:35 PM
!
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      MODULE DGTimeDerivativeMethods
      USE SMConstants
!
!     ========      
      CONTAINS 
!     ========      
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE LocalTimeDerivative( e, spA, t )
      USE ElementClass
      USE NodalStorageClass
      USE PhysicsStorage
      IMPLICIT NONE
!
!     -----------------
!     Input parameters:
!     -----------------
!
      TYPE(Element)      :: e
      TYPE(NodalStorage) :: spA
      REAL(KIND=RP)      :: t
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP), DIMENSION( 0:spA % N, &
                                0:spA % N, &
                                0:spA % N, &
                                N_EQN, 3 )  :: contravariantFlux
      
      CALL ComputeContravariantFlux( e, contravariantFlux )
      
      IF ( flowIsNavierStokes )     THEN
!         CALL AddViscousContravariantFluxes(  dgsem % mesh % elements(eID), contravariantFlux )
      END IF
      
      CALL ComputeDGDivergence( contravariantFlux, e, spA, e % Qdot ) !QDot saves the divergence
!
!     --------------------------------------------------------
!     Finish up - move divergence to left side of the equation
!     --------------------------------------------------------
!
      e % QDot = -e % QDot
    
      END SUBROUTINE LocalTimeDerivative
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ComputeContravariantFlux( e, contravariantFlux )
!
!     --------------------------------------
!     As described, compute
      
!     \[
!        \tilde f^i = J\vec a^i \cdot \vec F
!     \]
!     --------------------------------------
!      
      USE ElementClass
      USE PhysicsStorage
      USE Physics
      IMPLICIT NONE
!
!     -----------------
!     Input parameters:
!     -----------------
!
      TYPE(Element)                          :: e
      REAL(KIND=RP), dimension(0:,0:,0:,:,:) :: contravariantFlux 
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER                         :: n, m, l, nv
      REAL(KIND=RP), DIMENSION(N_EQN) :: ff, gg, hh
      
      DO l = 0, e % N
         DO m = 0, e % N
            DO n = 0, e % N

               CALL xflux( e % Q(n,m,l,:), ff )
               CALL yflux( e % Q(n,m,l,:), gg )
               CALL zflux( e % Q(n,m,l,:), hh )

               DO nv = 1, N_EQN
                  contravariantFlux(n,m,l,nv,1) = e % geom % jGradXi(1,n,m,l)  *ff(nv) +   &
                                                  e % geom % jGradXi(2,n,m,l)  *gg(nv) +   &
                                                  e % geom % jGradXi(3,n,m,l)  *hh(nv)
                  contravariantFlux(n,m,l,nv,2) = e % geom % jGradEta(1,n,m,l) *ff(nv) +   &
                                                  e % geom % jGradEta(2,n,m,l) *gg(nv) +   &
                                                  e % geom % jGradEta(3,n,m,l) *hh(nv)
                  contravariantFlux(n,m,l,nv,3) = e % geom % jGradZeta(1,n,m,l)*ff(nv) +   &
                                                  e % geom % jGradZeta(2,n,m,l)*gg(nv) +   &
                                                  e % geom % jGradZeta(3,n,m,l)*hh(nv)
               END DO
               
            END DO
         END DO
      END DO
    
      END SUBROUTINE ComputeContravariantFlux
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ComputeDGDivergence( contravariantFlux, e, spA, divFlux )
      USE ElementClass
      USE NodalStorageClass
      USE PhysicsStorage
!
!     --------------------------------------------------------
!     Compute the divergence of the input contravariant fluxes
!     --------------------------------------------------------
!
      IMPLICIT NONE 
!     ---------
!     Arguments
!     ---------
!
      TYPE(Element)                                               :: e
      REAL(KIND=RP), dimension(0:e % N,0:e % N,0:e % N,1:N_EQN,3) :: contravariantFlux 
      TYPE(NodalStorage)                                          :: spA
      REAL(KIND=RP), DIMENSION(0:e % N,0:e % N,0:e % N,1:N_EQN)   :: divFlux
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP), DIMENSION(0:e % N,1:N_EQN) :: fx
      
      INTEGER :: N, j, i, k, nv
      
      N = e % N
!
!     ---------------------------
!     Xi-derivative contributions
!     ---------------------------
!
      DO k = 0, N
         DO j = 0, N
            CALL DGSpaceDerivative( contravariantFlux(:,j,k,:,1), &
                                    e % FStarb(:,j,k,ELEFT), e % FStarb(:,j,k,ERIGHT), &
                                    N, spA % D, spA % b, fx )
            divFlux(:,j,k,:) = fx
         END DO
      END DO 
!
!     ----------------------------
!     Eta-Derivative Contributions
!     ----------------------------
!
      DO k = 0, N
         DO i = 0, N
            CALL DGSpaceDerivative( contravariantFlux(i,:,k,:,2), &
                                    e % FStarb(:,i,k,EFRONT), e % FStarb(:,i,k,EBACK), &
                                    N, spA % D, spA % b, fx )
            divFlux(i,:,k,:) = divFlux(i,:,k,:) + fx
         END DO
      END DO 
!
!     -----------------------------
!     Zeta-Derivative Contributions
!     -----------------------------
!
      DO j = 0, N
         DO i = 0, N
            CALL DGSpaceDerivative( contravariantFlux(i,j,:,:,3), &
                                    e % FStarb(:,i,j,EBOTTOM), e % FStarb(:,i,j,ETOP), &
                                    N, spA % D, spA % b, fx )
            divFlux(i,j,:,:) = divFlux(i,j,:,:) + fx
         END DO
      END DO 
!
!     ---------
!     Finish up
!     ---------
!
      DO k = 0, N
         DO j = 0, N 
            DO i = 0, N 
               DO nv = 1, N_EQN
                  divFlux(i,j,k,nv) = divFlux(i,j,k,nv)/e % geom % jacobian(i,j,k)
               END DO
            END DO 
         END DO
      END DO
      

      END SUBROUTINE ComputeDGDivergence
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE AddViscousContravariantFluxes( e, geom, fFlux, gFlux, N )
      USE PhysicsStorage
      USE MappedGeometryClass
      USE ElementClass
      IMPLICIT NONE
!
!     -----------------
!     Input parameters:
!     -----------------
!
      INTEGER                                         :: N
      TYPE(Element)                                   :: e
      REAL(KIND=RP)       , DIMENSION(0:N,0:N, N_EQN) :: fFlux
      REAL(KIND=RP)       , DIMENSION(0:N,0:N, N_EQN) :: gFlux
      TYPE(MappedGeometry)                            :: geom
!
!     ---------------
!     Local variables
!     ---------------
!
!      INTEGER                                :: i, j, M, k, indx
!      REAL(KIND=RP), DIMENSION(N_EQN)        :: ff, gg
!      REAL(KIND=RP), DIMENSION(2,N_GRAD_EQN) :: grad
!      INTEGER      , DIMENSION(2)            :: horizontalIndices = (/EBOTTOM, ETOP /)
!      INTEGER      , DIMENSION(2)            :: verticalIndices   = (/ELEFT, ERIGHT /)
!      
!      M = N
!
!     ----------------------
!     Interior contributions
!     ----------------------
!
!      DO j = 0, M
!         DO i = 0, N 
!            grad(1,:) = e % U_x(i,j,:)
!            grad(2,:) = e % U_y(i,j,:)
!            
!            CALL xDiffusiveFlux( e % Q(i,j,:), grad, ff )
!            CALL yDiffusiveFlux( e % Q(i,j,:), grad, gg )
!            
!            fFlux(i,j,:) =  fFlux(i,j,:) - ( geom % Y_eta(i,j)*ff - geom % X_eta(i,j)*gg )
!            gFlux(i,j,:) =  gFlux(i,j,:) - ( -geom % Y_xi(i,j)*ff + geom % X_xi(i,j) *gg )
!         END DO
!      END DO
!!
!!     ----------------------------------------------------------------------
!!     Boundary contributions
!!     At this point the boundary values of Qb (Storing Ub) and U_xb and U_yb
!!     have already been averaged.
!!     ----------------------------------------------------------------------
!!
!       DO k = 1, 2
!          indx = horizontalIndices(k)
!          DO i = 0, N
!            grad(1,:) = e % U_xb(:,i,indx)
!            grad(2,:) = e % U_yb(:,i,indx)
!            
!            CALL xDiffusiveFlux( e % Qb(:,i,indx), grad, ff )
!            CALL yDiffusiveFlux( e % Qb(:,i,indx), grad, gg )
!            
!            e % FStarb(:,i,indx) = e % FStarb(:,i,indx) - &
!                                   (ff*geom % normal(i,1,indx) + gg*geom % normal(i,2,indx))*geom % scal(i,indx)
!          END DO
!       END DO
!       
!       DO k = 1, 2
!          indx = verticalIndices(k)
!          DO j = 0, N
!            grad(1,:) = e % U_xb(:,j,indx)
!            grad(2,:) = e % U_yb(:,j,indx)
!            
!            CALL xDiffusiveFlux( e % Qb(:,j,indx), grad, ff )
!            CALL yDiffusiveFlux( e % Qb(:,j,indx), grad, gg )
!            
!            e % FStarb(:,j,indx) = e % FStarb(:,j,indx) - &
!                                   (ff*geom % normal(j,1,indx) + gg*geom % normal(j,2,indx))*geom % scal(j,indx)
!          END DO
!       END DO
    
      END SUBROUTINE AddViscousContravariantFluxes
!
!////////////////////////////////////////////////////////////////////////
!
!@mark -
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ProlongToFaces( e, spA )
!
!     -----------------------------------------------------------
!     For Gauss point approximations, we interpolate to each face
!     of the element and store the result in the face solution 
!     array, Qb
!     -----------------------------------------------------------
!
         USE PhysicsStorage
         USE NodalStorageClass
         USE ElementClass
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(NodalStorage) :: spA
         TYPE(Element)      :: e
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER :: N, i, j, k, nv
         
         N = e % N
!
!        --------------
!        Left and right
!        --------------
!
         DO k = 0, N
            DO j = 0, N 
               DO nv = 1, N_EQN 
                  CALL InterpolateToBoundary( e % Q(:,j,k,nv), spA % v(:,LEFT) , N, e % Qb(nv,j,k,ELEFT) )
                  CALL InterpolateToBoundary( e % Q(:,j,k,nv), spA % v(:,RIGHT), N, e % Qb(nv,j,k,ERIGHT))
               END DO
            END DO
         END DO 
!
!        --------------
!        Front and back
!        --------------
!
         DO k = 0, N
            DO i = 0, N 
               DO nv = 1, N_EQN 
                  CALL InterpolateToBoundary( e % Q(i,:,k,nv), spA % v(:,FRONT), N, e % Qb(nv,i,k,EFRONT) )
                  CALL InterpolateToBoundary( e % Q(i,:,k,nv), spA % v(:,BACK) , N, e % Qb(nv,i,k,EBACK) )
               END DO
            END DO
         END DO 
!
!        --------------
!        Bottom and Top
!        --------------
!
         DO j = 0, N
            DO i = 0, N 
               DO nv = 1, N_EQN 
                  CALL InterpolateToBoundary( e % Q(i,j,:,nv), spA % v(:,BOTTOM), N, e % Qb(nv,i,j,EBOTTOM) )
                  CALL InterpolateToBoundary( e % Q(i,j,:,nv), spA % v(:,TOP)   , N, e % Qb(nv,i,j,ETOP)    )
               END DO
            END DO
         END DO 

      END SUBROUTINE ProlongToFaces
!
!////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE DGSpaceDerivative( u, uL, uR,  N, D, b, deriv ) 
!
!  ----------------------------------------------------------
!  The one dimensional space derivative for the DG method for
!  a vector state. This is Algorithm 92 in the book.
!  ----------------------------------------------------------
!
      USE PhysicsStorage
      USE PolynomialInterpAndDerivsModule
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      INTEGER                              :: N
      REAL(KIND=RP), DIMENSION(0:N, N_EQN) :: u
      REAL(KIND=RP), DIMENSION(N_EQN)      :: uL, uR
      REAL(KIND=RP), DIMENSION(0:N, 0:N)   :: D
      REAL(KIND=RP), DIMENSION(0:N, 2)     :: b
!
!     -----------------
!     Output variables:
!     -----------------
!
      REAL(KIND=RP), DIMENSION(0:N,N_EQN), INTENT(OUT) :: deriv
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER            :: j, k
      INTEGER, PARAMETER :: LEFT = 1, RIGHT = 2
!
!     ----------------------------
!     Internal points contribution
!     ----------------------------
!
      DO k = 1, N_EQN
         CALL PolyDirectMatrixMultiplyDeriv( u(:,k), deriv(:,k), D, N )
      END DO
!
!     ----------------------------
!     Boundary points contribution
!     ----------------------------
!
      DO j = 0,N  
         deriv(j,:) = deriv(j,:) + uR*b(j,RIGHT) + uL*b(j,LEFT)
      END DO
     
   END SUBROUTINE DGSpaceDerivative
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE InterpolateToBoundary( u, v, N, bValue )
!
!     -------------------------------------------------------------
!     Interpolation to the boundary is a dot product for each row or
!     column. Using here the intrinsic Fortran function, without
!     having tested that it is faster or slower than a direct
!     computation for the values of N that we used in the DGSEM.
!     -------------------------------------------------------------
!
         USE SMConstants
         IMPLICIT NONE
         INTEGER                      , INTENT(IN)  :: N
         REAL(KIND=RP), DIMENSION(0:N), INTENT(IN)  :: u, v
         REAL(KIND=RP)                , INTENT(OUT) :: bValue
!
         bValue = DOT_PRODUCT( u, v )
         
      END SUBROUTINE InterpolateToBoundary

   END MODULE DGTimeDerivativeMethods
