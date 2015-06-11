!
!////////////////////////////////////////////////////////////////////////
!
!      SpatialApproximations.f90
!      Created: 2011-07-28 15:44:27 -0400 
!      By: David Kopriva  
!
!      Algorithms:
!         Algorithm 112: ProlongToFaces (DG2DProlongToFaces)
!         Algorithm  61: InterpolateToBoundary
!         Algorithm  92: DGSpaceDerivative (SystemDGDerivative)
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ProlongToFaces( e, spA )
!
!     -------------------------------------
!     Prolong the solution Q onto the faces
!     -------------------------------------
!
         USE Nodal2DStorageClass
         USE MappedGeometryClass
         USE ElementClass
         IMPLICIT NONE
         
         TYPE(Nodal2DStorage)    :: spA
         TYPE(Element)           :: e
         
         INTEGER :: nEqn, N, M, i, j, k
         
         nEqn = SIZE(e % Q,3)
         N    = spA % N
         M    = N
!
!        --------------
!        Left and right
!        --------------
!
         DO j = 0, M 
            DO k = 1, nEqn 
               CALL InterpolateToBoundary( e % Q(:,j,k), spA % v(:,LEFT) , N, e % Qb(k,j,ELEFT) )
               CALL InterpolateToBoundary( e % Q(:,j,k), spA % v(:,RIGHT), N, e % Qb(k,j,ERIGHT))
            END DO
         END DO
!
!        --------------
!        Top and bottom
!        --------------
!
         DO i = 0, N 
            DO k = 1, nEqn 
               CALL InterpolateToBoundary( e % Q(i,:,k), spA % v(:,BOTTOM), M, e % Qb(k,i,EBOTTOM) )
               CALL InterpolateToBoundary( e % Q(i,:,k), spA % v(:,TOP)   , M, e % Qb(k,i,ETOP) )
            END DO
         END DO
         
      END SUBROUTINE ProlongToFaces
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ProlongGradientToFaces( e, spA )
!
!     -----------------------------------------
!     Prolong the graident gradU onto the faces
!     Store the results in Qb.
!     -----------------------------------------
!
         USE Nodal2DStorageClass
         USE MappedGeometryClass
         USE ElementClass
         USE PhysicsStorage, ONLY: N_GRAD_EQN
         IMPLICIT NONE
         
         TYPE(Nodal2DStorage)    :: spA
         TYPE(MappedGeometry)    :: geom
         TYPE(Element)           :: e
         
         INTEGER :: nEqn, N, M, i, j, k
         
         nEqn = SIZE(e % Q,3)
         N    = spA % N
         M    = N
!
!        --------------
!        Left and right
!        --------------
!
         DO j = 0, M 
            DO k = 1, N_GRAD_EQN 
               CALL InterpolateToBoundary( e % U_x(:,j,k), spA % v(:,LEFT) , N, e % U_xb(k,j,ELEFT) )
               CALL InterpolateToBoundary( e % U_x(:,j,k), spA % v(:,RIGHT), N, e % U_xb(k,j,ERIGHT))
               CALL InterpolateToBoundary( e % U_y(:,j,k), spA % v(:,LEFT) , N, e % U_yb(k,j,ELEFT) )
               CALL InterpolateToBoundary( e % U_y(:,j,k), spA % v(:,RIGHT), N, e % U_yb(k,j,ERIGHT))
            END DO
         END DO
!
!        --------------
!        Top and bottom
!        --------------
!
         DO i = 0, N 
            DO k = 1, N_GRAD_EQN
               CALL InterpolateToBoundary( e % U_x(i,:,k), spA % v(:,BOTTOM), M, e % U_xb(k,i,EBOTTOM) )
               CALL InterpolateToBoundary( e % U_x(i,:,k), spA % v(:,TOP)   , M, e % U_xb(k,i,ETOP) )
               CALL InterpolateToBoundary( e % U_y(i,:,k), spA % v(:,BOTTOM), M, e % U_yb(k,i,EBOTTOM) )
               CALL InterpolateToBoundary( e % U_y(i,:,k), spA % v(:,TOP)   , M, e % U_yb(k,i,ETOP) )
            END DO
         END DO
         
      END SUBROUTINE ProlongGradientToFaces
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ComputeDGDivergence( fFlux, gFlux, FStarB, jacobian, N, nEqn, D, b, divFlux )
!
!     -------------------------------------------------------------------
!     Compute the divergence of the input contravariant fluxes
!     -------------------------------------------------------------------
!
      USE SMConstants
      USE MappedGeometryClass, ONLY: ELEFT, ERIGHT, ETOP, EBOTTOM
      IMPLICIT NONE 
!     ---------
!     Arguments
!     ---------
!
      INTEGER                                      :: N, nEqn
      REAL(KIND=RP), DIMENSION(0:N,0:N,1:nEqn)     :: fFlux
      REAL(KIND=RP), DIMENSION(0:N,0:N,1:nEqn)     :: gFlux
      REAL(KIND=RP), DIMENSION(0:N,0:N)            :: D
      REAL(KIND=RP), DIMENSION(1:nEqn,0:N,4)       :: FStarb
      REAL(KIND=RP), DIMENSION(0:N,2)              :: b
      REAL(KIND=RP), DIMENSION(0:N,0:N,1:nEqn)     :: divFlux
      REAL(KIND=RP), DIMENSION(0:N,0:N)            :: jacobian
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP), DIMENSION(0:N,nEqn) :: fx
      REAL(KIND=RP), DIMENSION(0:N,nEqn) :: gy
      
      INTEGER :: M, j, i, k
      
      M = N
!
!     --------------------------
!     X-derivative contributions
!     --------------------------
!
      DO j = 0, M
         CALL DGSpaceDerivative( fFlux(:,j,:), FStarb(:,j,ELEFT), FStarb(:,j,ERIGHT), N, nEqn, D, b, fx )
         divFlux(:,j,:) = fx
      END DO
!
!     ----------------------------
!     Eta-Derivative Contributions
!     ----------------------------
!
      DO i = 0, N
         CALL DGSpaceDerivative( gFlux(i,:,:), FStarb(:,i,EBOTTOM), FStarb(:,i,ETOP), M, nEqn, D, b, gy )
         divFlux(i,:,:) = divFlux(i,:,:) + gy
      END DO
!
!     ---------
!     Finish up
!     ---------
!
      DO j = 0, M 
         DO i = 0, N 
            DO k = 1, nEqn
               divFlux(i,j,k) = divFlux(i,j,k)/jacobian(i,j)
            END DO
         END DO
      END DO
      

      END SUBROUTINE ComputeDGDivergence
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ComputeDGGradient( Q, Ub, geom, N, D, b, U_x, U_y )
!
!     ---------------------------------
!     Compute the gradient of the input
!     ---------------------------------
!
      USE SMConstants
      USE MappedGeometryClass
      USE PDEModule
      IMPLICIT NONE 
!     ---------
!     Arguments
!     ---------
!
      INTEGER                                    :: N
      
      REAL(KIND=RP), DIMENSION(0:N,0:N,1:N_EQN)  :: Q
      REAL(KIND=RP), DIMENSION(N_GRAD_EQN,0:N,4) :: Ub
      
      TYPE(MappedGeometry)                       :: geom

      REAL(KIND=RP), DIMENSION(0:N,0:N)          :: D
      REAL(KIND=RP), DIMENSION(0:N,2)            :: b
      
      REAL(KIND=RP), DIMENSION(0:N,0:N,1:N_GRAD_EQN)  :: U_x, U_y
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP), DIMENSION(N_GRAD_EQN,0:N,4)     :: fStarB
      REAL(KIND=RP), DIMENSION(0:N,0:N,1:N_GRAD_EQN) :: f, g
      real(KIND=RP), dimension(N_GRAD_EQN)           :: U
      
      INTEGER            :: M, j, i
      
      M = N
!
!     ------------------------------------------------------------------
!                                 \hat x component
!     ------------------------------------------------------------------
!
!     ---------------
!     Interior fluxes
!     ---------------
!
      DO j = 0, M 
         DO i = 0, N
            CALL GradientValuesForQ( Q(i,j,:), U )
            f(i,j,:) =  U*geom % Y_eta(i,j)
            g(i,j,:) = -U*geom % Y_xi(i,j)
         END DO
      END DO
!
!     ---------------
!     Boundary fluxes
!     ---------------
!
      DO j = 0, M 
         fStarB(:,j,ELEFT)  = geom % normal(j,1,ELEFT)*geom % scal(j, ELEFT)*Ub(1:N_GRAD_EQN,j,ELEFT)
         fStarB(:,j,ERIGHT) = geom % normal(j,1,ERIGHT)*geom % scal(j,ERIGHT)*Ub(1:N_GRAD_EQN,j,ERIGHT)
      END DO
      
      DO i = 0, N 
         fStarB(:,i,EBOTTOM)  =  geom % normal(i,1, EBOTTOM)*geom % scal(i,EBOTTOM)*Ub(1:N_GRAD_EQN,i,EBOTTOM)
         fStarB(:,i,ETOP)     =  geom % normal(i,1, ETOP)*geom % scal(i,ETOP)      *Ub(1:N_GRAD_EQN,i,ETOP)
      END DO
!
!     -------------------------------------------
!     Computing the divergence gives the gradient
!     -------------------------------------------
!
      CALL ComputeDGDivergence( f, g, FStarB, geom % jacobian, N, N_GRAD_EQN, D, b, U_x )
!
!     ------------------------------------------------------------------
!                                 \hat y component
!     ------------------------------------------------------------------
!
!     ---------------
!     Interior fluxes
!     ---------------
!
      DO j = 0, N 
         DO i = 0, N 
            CALL GradientValuesForQ( Q(i,j,:), U )
            f(i,j,:) = -U*geom % X_eta(i,j)
            g(i,j,:) =  U*geom % X_xi(i,j)
         END DO
      END DO
      
      DO j = 0, M 
         fStarB(:,j,ELEFT)  =  geom % normal(j,2, ELEFT)*geom % scal(j, ELEFT)*Ub(1:N_GRAD_EQN,j,ELEFT)
         fStarB(:,j,ERIGHT) =  geom % normal(j,2,ERIGHT)*geom % scal(j,ERIGHT)*Ub(1:N_GRAD_EQN,j,ERIGHT)
      END DO
      
      DO i = 0, N 
         fStarB(:,i,EBOTTOM)  =  geom % normal(i,2, EBOTTOM)*geom % scal(i,EBOTTOM)*Ub(1:N_GRAD_EQN,i,EBOTTOM)
         fStarB(:,i,ETOP)     =  geom % normal(i,2, ETOP)   *geom % scal(i,ETOP)   *Ub(1:N_GRAD_EQN,i,ETOP)
      END DO
!
!     -------------------------------------------
!     Computing the divergence gives the gradient
!     -------------------------------------------
!
      CALL ComputeDGDivergence( f, g, FStarB, geom % jacobian, N, N_GRAD_EQN, D, b, U_y )

      END SUBROUTINE ComputeDGGradient
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE InterpolateToBoundary( u, v, N, bValue )
!
!     -------------------------------------------------------------
!     Interpolation to the boundary is matrix-vector multiplication
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
!
!////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE DGSpaceDerivative( u, uL, uR,  N, nEqn, D, b, deriv ) 
      USE PolynomialInterpAndDerivsModule
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      INTEGER                             :: N, nEqn
      REAL(KIND=RP), DIMENSION(0:N, nEqn) :: u
      REAL(KIND=RP), DIMENSION(nEqn)      :: uL, uR
      REAL(KIND=RP), DIMENSION(0:N, 0:N)  :: D
      REAL(KIND=RP), DIMENSION(0:N, 2)    :: b
!
!     -----------------
!     Output variables:
!     -----------------
!
      REAL(KIND=RP), DIMENSION(0:N,nEqn), INTENT(OUT) :: deriv
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER :: j, k
      INTEGER, PARAMETER :: LEFT = 1, RIGHT = 2
!
!     ----------------------------
!     Internal points contribution
!     ----------------------------
!
      DO k = 1, nEqn
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

