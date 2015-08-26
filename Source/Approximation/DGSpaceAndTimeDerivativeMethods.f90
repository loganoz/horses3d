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
         CALL AddViscousContravariantFluxes(  e, contravariantFlux )
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
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ComputeDGGradient( e, spA, t )
      USE ElementClass
      USE NodalStorageClass
      USE PhysicsStorage
      USE Physics, ONLY : GradientValuesForQ
!
!     ---------------------------------
!     Compute the gradient of the input
!     ---------------------------------
!

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
      REAL(KIND=RP), DIMENSION(0:e % N, 0:e % N, 0:e % N, N_GRAD_EQN) :: f
      REAL(KIND=RP), DIMENSION(0:e % N, 0:e % N, 0:e % N, N_GRAD_EQN) :: g
      REAL(KIND=RP), DIMENSION(0:e % N, 0:e % N, 0:e % N, N_GRAD_EQN) :: h
      REAL(KIND=RP), DIMENSION(N_GRAD_EQN) :: U
      REAL(KIND=RP), DIMENSION(N_EQN) :: Q
      REAL(KIND=RP), DIMENSION(0:e % N, N_GRAD_EQN) :: fx
      
      INTEGER :: N, j, i, k, nv

      
      N = e % N

      DO k = 0, e % N
         DO j = 0, e % N
            DO i = 0, e % N
               Q = e % Q(i,j,k,:)

               CALL GradientValuesForQ( Q = Q, U = U )
               
               !CALL xflux( e % U(n,m,l,:), ff )
               !CALL yflux( e % U(n,m,l,:), gg )
               !CALL zflux( e % U(n,m,l,:), hh )

               DO nv = 1, N_GRAD_EQN
               f(i,j,k,nv) = U(nv) * e % geom % jGradXi(1,i,j,k)
               g(i,j,k,nv) = U(nv) * e % geom % jGradEta(1,i,j,k)
               h(i,j,k,nv) = U(nv) * e % geom % jGradZeta(1,i,j,k) 
               
!                  contravariantFlux(n,m,l,:,1) = e % geom % jGradXi(1,n,m,l)  *ff(:) +   &
!                                                  e % geom % jGradXi(2,n,m,l)  *gg(:) +   &
!                                                  e % geom % jGradXi(3,n,m,l)  *hh(:)
!                  contravariantFlux(n,m,l,:,2) = e % geom % jGradEta(1,n,m,l) *ff(:) +   &
!                                                  e % geom % jGradEta(2,n,m,l) *gg(:) +   &
!                                                  e % geom % jGradEta(3,n,m,l) *hh(:)
!                  contravariantFlux(n,m,l,:,3) = e % geom % jGradZeta(1,n,m,l)*ff(:) +   &
!                                                  e % geom % jGradZeta(2,n,m,l)*gg(:) +   &
!                                                  e % geom % jGradZeta(3,n,m,l)*hh(:)
               END DO
               
            END DO
         END DO
      END DO
!
!     ---------------------------
!     Xi-derivative contributions
!     ---------------------------
!
      DO k = 0, N
         DO j = 0, N
            CALL DGGradSpaceDerivative( f(:,j,k,:), &
                                    e % geom%normal(1,j,k,ELEFT) * e % geom%scal(j,k, ELEFT) * e % Ub(:,j,k,ELEFT), &
                                    e % geom%normal(1,j,k,ERIGHT) * e % geom%scal(j,k, ERIGHT) * e % Ub(:,j,k,ERIGHT), &
                                    N, spA % D, spA % b, fx )
            e % U_x(:,j,k,:) = fx
         END DO
      END DO 
!
!     ----------------------------
!     Eta-Derivative Contributions
!     ----------------------------
!
      DO k = 0, N
         DO i = 0, N
            CALL DGGradSpaceDerivative( g(i,:,k,:), &
                                    e % geom%normal(1,i,k,EFRONT) * e % geom%scal(i,k, EFRONT) * e % Ub(:,i,k,EFRONT), &
                                    e % geom%normal(1,i,k,EBACK) * e % geom%scal(i,k, EBACK) * e % Ub(:,i,k,EBACK), &
                                    N, spA % D, spA % b, fx )
            e % U_x(i,:,k,:) = e % U_x(i,:,k,:) + fx
         END DO
      END DO 
!
!     -----------------------------
!     Zeta-Derivative Contributions
!     -----------------------------
!
      DO j = 0, N
         DO i = 0, N
            CALL DGGradSpaceDerivative( h(i,j,:,:), &
                                    e % geom%normal(1,i,j,EBOTTOM) * e % geom%scal(i,j, EBOTTOM) * e % Ub(:,i,j,EBOTTOM), &
                                    e % geom%normal(1,i,j,ETOP) * e % geom%scal(i,j, ETOP) * e % Ub(:,i,j,ETOP), &
                                    N, spA % D, spA % b, fx )
            e % U_x(i,j,:,:) = e % U_x(i,j,:,:) + fx
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
               DO nv = 1, N_GRAD_EQN
                  e % U_x(i,j,k,nv) = e % U_x(i,j,k,nv) / e % geom % jacobian(i,j,k) 
               END DO
            END DO 
         END DO
      END DO
      
      DO k = 0, e % N
         DO j = 0, e % N
            DO i = 0, e % N

               CALL GradientValuesForQ( e % Q(i,j,k,:), U )

               !CALL xflux( e % U(n,m,l,:), ff )
               !CALL yflux( e % U(n,m,l,:), gg )
               !CALL zflux( e % U(n,m,l,:), hh )

               DO nv = 1, N_GRAD_EQN
               f(i,j,k,nv) = U(nv) * e % geom % jGradXi(2,i,j,k)
               g(i,j,k,nv) = U(nv) * e % geom % jGradEta(2,i,j,k)
               h(i,j,k,nv) = U(nv) * e % geom % jGradZeta(2,i,j,k) 
               
!                  contravariantFlux(n,m,l,:,1) = e % geom % jGradXi(1,n,m,l)  *ff(:) +   &
!                                                  e % geom % jGradXi(2,n,m,l)  *gg(:) +   &
!                                                  e % geom % jGradXi(3,n,m,l)  *hh(:)
!                  contravariantFlux(n,m,l,:,2) = e % geom % jGradEta(1,n,m,l) *ff(:) +   &
!                                                  e % geom % jGradEta(2,n,m,l) *gg(:) +   &
!                                                  e % geom % jGradEta(3,n,m,l) *hh(:)
!                  contravariantFlux(n,m,l,:,3) = e % geom % jGradZeta(1,n,m,l)*ff(:) +   &
!                                                  e % geom % jGradZeta(2,n,m,l)*gg(:) +   &
!                                                  e % geom % jGradZeta(3,n,m,l)*hh(:)
               END DO
               
            END DO
         END DO
      END DO
!
!     ---------------------------
!     Xi-derivative contributions
!     ---------------------------
!
      DO k = 0, N
         DO j = 0, N
            CALL DGGradSpaceDerivative( f(:,j,k,:), &
                                    e % geom%normal(2,j,k,ELEFT) * e % geom%scal(j,k, ELEFT) * e % Ub(:,j,k,ELEFT), &
                                    e % geom%normal(2,j,k,ERIGHT) * e % geom%scal(j,k, ERIGHT) * e % Ub(:,j,k,ERIGHT), &
                                    N, spA % D, spA % b, fx )
            e % U_y(:,j,k,:) = fx
         END DO
      END DO 
!
!     ----------------------------
!     Eta-Derivative Contributions
!     ----------------------------
!
      DO k = 0, N
         DO i = 0, N
            CALL DGGradSpaceDerivative( g(i,:,k,:), &
                                    e % geom%normal(2,i,k,EFRONT) * e % geom%scal(i,k, EFRONT) * e % Ub(:,i,k,EFRONT), &
                                    e % geom%normal(2,i,k,EBACK) * e % geom%scal(i,k, EBACK) * e % Ub(:,i,k,EBACK), &
                                    N, spA % D, spA % b, fx )
            e % U_y(i,:,k,:) = e % U_y(i,:,k,:) + fx
         END DO
      END DO 
!
!     -----------------------------
!     Zeta-Derivative Contributions
!     -----------------------------
!
      DO j = 0, N
         DO i = 0, N
            CALL DGGradSpaceDerivative( h(i,j,:,:), &
                                    e % geom%normal(2,i,j,EBOTTOM) * e % geom%scal(i,j, EBOTTOM) * e % Ub(:,i,j,EBOTTOM), &
                                    e % geom%normal(2,i,j,ETOP) * e % geom%scal(i,j, ETOP) * e % Ub(:,i,j,ETOP), &
                                    N, spA % D, spA % b, fx )
            e % U_y(i,j,:,:) = e % U_y(i,j,:,:) + fx
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
               DO nv = 1, N_GRAD_EQN
                  e % U_y(i,j,k,nv) = e % U_y(i,j,k,nv) / e % geom % jacobian(i,j,k) 
               END DO
            END DO 
         END DO
      END DO

      DO k = 0, e % N
         DO j = 0, e % N
            DO i = 0, e % N

               CALL GradientValuesForQ( e % Q(i,j,k,:), U )

               !CALL xflux( e % U(n,m,l,:), ff )
               !CALL yflux( e % U(n,m,l,:), gg )
               !CALL zflux( e % U(n,m,l,:), hh )

               DO nv = 1, N_GRAD_EQN
               f(i,j,k,nv) = U(nv) * e % geom % jGradXi(3,i,j,k)
               g(i,j,k,nv) = U(nv) * e % geom % jGradEta(3,i,j,k)
               h(i,j,k,nv) = U(nv) * e % geom % jGradZeta(3,i,j,k) 
               
!                  contravariantFlux(n,m,l,:,1) = e % geom % jGradXi(1,n,m,l)  *ff(:) +   &
!                                                  e % geom % jGradXi(2,n,m,l)  *gg(:) +   &
!                                                  e % geom % jGradXi(3,n,m,l)  *hh(:)
!                  contravariantFlux(n,m,l,:,2) = e % geom % jGradEta(1,n,m,l) *ff(:) +   &
!                                                  e % geom % jGradEta(2,n,m,l) *gg(:) +   &
!                                                  e % geom % jGradEta(3,n,m,l) *hh(:)
!                  contravariantFlux(n,m,l,:,3) = e % geom % jGradZeta(1,n,m,l)*ff(:) +   &
!                                                  e % geom % jGradZeta(2,n,m,l)*gg(:) +   &
!                                                  e % geom % jGradZeta(3,n,m,l)*hh(:)
               END DO
               
            END DO
         END DO
      END DO
!
!     ---------------------------
!     Xi-derivative contributions
!     ---------------------------
!
      DO k = 0, N
         DO j = 0, N
            CALL DGGradSpaceDerivative( f(:,j,k,:), &
                                    e % geom%normal(3,j,k,ELEFT) * e % geom%scal(j,k, ELEFT) * e % Ub(:,j,k,ELEFT), &
                                    e % geom%normal(3,j,k,ERIGHT) * e % geom%scal(j,k, ERIGHT) * e % Ub(:,j,k,ERIGHT), &
                                    N, spA % D, spA % b, fx )
            e % U_z(:,j,k,:) = fx
         END DO
      END DO 
!
!     ----------------------------
!     Eta-Derivative Contributions
!     ----------------------------
!
      DO k = 0, N
         DO i = 0, N
            CALL DGGradSpaceDerivative( g(i,:,k,:), &
                                    e % geom%normal(3,i,k,EFRONT) * e % geom%scal(i,k, EFRONT) * e % Ub(:,i,k,EFRONT), &
                                    e % geom%normal(3,i,k,EBACK) * e % geom%scal(i,k, EBACK) * e % Ub(:,i,k,EBACK), &
                                    N, spA % D, spA % b, fx )
            e % U_z(i,:,k,:) = e % U_z(i,:,k,:) + fx
         END DO
      END DO 
!
!     -----------------------------
!     Zeta-Derivative Contributions
!     -----------------------------
!
      DO j = 0, N
         DO i = 0, N
            CALL DGGradSpaceDerivative( h(i,j,:,:), &
                                    e % geom%normal(3,i,j,EBOTTOM) * e % geom%scal(i,j, EBOTTOM) * e % Ub(:,i,j,EBOTTOM), &
                                    e % geom%normal(3,i,j,ETOP) * e % geom%scal(i,j, ETOP) * e % Ub(:,i,j,ETOP), &
                                    N, spA % D, spA % b, fx )
            e % U_z(i,j,:,:) = e % U_z(i,j,:,:) + fx
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
               DO nv = 1, N_GRAD_EQN
                  e % U_z(i,j,k,nv) = e % U_z(i,j,k,nv) / e % geom % jacobian(i,j,k) 
               END DO
            END DO 
         END DO
      END DO
      



    
      END SUBROUTINE ComputeDGGradient
      
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE AddViscousContravariantFluxes( e, contravariantFlux )
      USE PhysicsStorage
      USE MappedGeometryClass
      USE ElementClass
      USE Physics
      IMPLICIT NONE
!
!     -----------------
!     Input parameters:
!     -----------------
!
      TYPE(Element)                                   :: e
      REAL(KIND=RP), DIMENSION( 0:e % N, &
                                0:e % N, &
                                0:e % N, &
                                N_EQN, 3 )  :: contravariantFlux
!
!     -----------------
!     Local variables:
!     -----------------
!                                
      REAL(KIND=RP)   :: grad(3, N_GRAD_EQN)
      REAL(KIND=RP)   :: ff(N_EQN), gg(N_EQN), hh(N_EQN)
      INTEGER         :: l,m,n,nv,i,j,k
!
!     ----------------------
!     Interior contributions
!     ----------------------
!
      DO l = 0, e % N
         DO m = 0, e % N
            DO n = 0, e % N

               grad(1,:) = e % U_x(n,m,l,:)
               grad(2,:) = e % U_y(n,m,l,:)
               grad(3,:) = e % U_z(n,m,l,:)

               CALL xDiffusiveFlux( e % Q(n,m,l,:), grad, ff )
               CALL yDiffusiveFlux( e % Q(n,m,l,:), grad, gg )
               CALL zDiffusiveFlux( e % Q(n,m,l,:), grad, hh )

               DO nv = 1, N_EQN
                  
                  contravariantFlux(n,m,l,nv,1) =   contravariantFlux(n,m,l,nv,1) -          &
                                                  ( e % geom % jGradXi(1,n,m,l)  *ff(nv) +   &
                                                    e % geom % jGradXi(2,n,m,l)  *gg(nv) +   &
                                                    e % geom % jGradXi(3,n,m,l)  *hh(nv) )
                  contravariantFlux(n,m,l,nv,2) = contravariantFlux(n,m,l,nv,2)   -          & 
                                                  ( e % geom % jGradEta(1,n,m,l) *ff(nv) +   &
                                                    e % geom % jGradEta(2,n,m,l) *gg(nv) +   &
                                                    e % geom % jGradEta(3,n,m,l) *hh(nv) )
                  contravariantFlux(n,m,l,nv,3) = contravariantFlux(n,m,l,nv,3)   -          & 
                                                  ( e % geom % jGradZeta(1,n,m,l)*ff(nv) +   &
                                                    e % geom % jGradZeta(2,n,m,l)*gg(nv) +   &
                                                    e % geom % jGradZeta(3,n,m,l)*hh(nv) )   
!                  IF (ABS( ( e % geom % jGradXi(1,n,m,l)  *ff(nv) +   &
!                                                    e % geom % jGradXi(2,n,m,l)  *gg(nv) +   &
!                                                    e % geom % jGradXi(3,n,m,l)  *hh(nv) ) )  > 1.d-13 ) THEN
                                                    
!                  PRINT*, "n,m,l"
!                  PRINT*, n,m,l
!                  PRINT*, contravariantFlux(n,m,l,nv,:)
!                  PRINT*, ( e % geom % jGradXi(1,n,m,l)  *ff(nv) +   &
!                                                    e % geom % jGradXi(2,n,m,l)  *gg(nv) +   &
!                                                    e % geom % jGradXi(3,n,m,l)  *hh(nv) ) 
!                  PRINT*, "----------------------------"
!                  
!                  ENDIF 
!                  
!                  IF (ABS( ( e % geom % jGradEta(1,n,m,l) *ff(nv) +   &
!                                                    e % geom % jGradEta(2,n,m,l) *gg(nv) +   &
!                                                    e % geom % jGradEta(3,n,m,l) *hh(nv) ) )  > 1.d-13 ) THEN
!                                                    
!                  PRINT*, "n,m,l"
!                  PRINT*, n,m,l
!                  PRINT*, contravariantFlux(n,m,l,nv,:)
!                  PRINT*,  ( e % geom % jGradEta(1,n,m,l) *ff(nv) +   &
!                                                    e % geom % jGradEta(2,n,m,l) *gg(nv) +   &
!                                                    e % geom % jGradEta(3,n,m,l) *hh(nv) )
!                  PRINT*, "----------------------------"
!                  
!                  ENDIF 
!                  
!                  IF (ABS( ( e % geom % jGradZeta(1,n,m,l)*ff(nv) +   &
!                                                    e % geom % jGradZeta(2,n,m,l)*gg(nv) +   &
!                                                    e % geom % jGradZeta(3,n,m,l)*hh(nv) ) )  > 1.d-13 ) THEN
!                                                    
!                  PRINT*, "n,m,l"
!                  PRINT*, n,m,l
!                  PRINT*, contravariantFlux(n,m,l,nv,:)
!                  PRINT*, ( e % geom % jGradZeta(1,n,m,l)*ff(nv) +   &
!                                                    e % geom % jGradZeta(2,n,m,l)*gg(nv) +   &
!                                                    e % geom % jGradZeta(3,n,m,l)*hh(nv) ) 
!                  PRINT*, "----------------------------"
!                  
!                  ENDIF 
                                                      
               END DO
               
            END DO
         END DO
      END DO
!
!     ----------------------------------------------------------------------
!     Boundary contributions
!     At this point the boundary values of Qb (Storing Ub)? and U_xb and U_yb
!     have already been averaged.
!     ----------------------------------------------------------------------
!      
      DO k = 1, 6
          DO j = 0, e % N
             DO i = 0, e % N
               grad(1,:) = e % U_xb(:,i,j,k)
               grad(2,:) = e % U_yb(:,i,j,k)
               grad(3,:) = e % U_zb(:,i,j,k)
            
               CALL xDiffusiveFlux( e % Qb(:,i,j,k), grad, ff )
               CALL yDiffusiveFlux( e % Qb(:,i,j,k), grad, gg )
               CALL zDiffusiveFlux( e % Qb(:,i,j,k), grad, hh )
            
               e % FStarb(:,i,j,k) = e % FStarb(:,i,j,k) - &
                                   ( ff * e % geom % normal(1,i,j,k) + &
                                   & gg * e % geom % normal(2,i,j,k) + &
                                   & hh * e % geom % normal(3,i,j,k)) &
                                   * e % geom % scal(i,j,k)
!               IF (MAXVAL (ABS( ff * e % geom % normal(1,i,j,k) + &
!                                   & gg * e % geom % normal(2,i,j,k) + &
!                                   & hh * e % geom % normal(3,i,j,k)) &
!                                   * e % geom % scal(i,j,k))>1.d-13) THEN 
!                                   PRINT*, "error"
!               ENDIF 
             END DO 
          END DO
      ENDDO 
          
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
      SUBROUTINE ProlongGradientToFaces( e, spA )
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
               DO nv = 1, N_GRAD_EQN 
                  CALL InterpolateToBoundary( e % U_x(:,j,k,nv), spA % v(:,LEFT) , N, e % U_xb(nv,j,k,ELEFT) )
                  CALL InterpolateToBoundary( e % U_x(:,j,k,nv), spA % v(:,RIGHT), N, e % U_xb(nv,j,k,ERIGHT))
                  CALL InterpolateToBoundary( e % U_y(:,j,k,nv), spA % v(:,LEFT) , N, e % U_yb(nv,j,k,ELEFT) )
                  CALL InterpolateToBoundary( e % U_y(:,j,k,nv), spA % v(:,RIGHT), N, e % U_yb(nv,j,k,ERIGHT))
                  CALL InterpolateToBoundary( e % U_z(:,j,k,nv), spA % v(:,LEFT) , N, e % U_zb(nv,j,k,ELEFT) )
                  CALL InterpolateToBoundary( e % U_z(:,j,k,nv), spA % v(:,RIGHT), N, e % U_zb(nv,j,k,ERIGHT))                  
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
               DO nv = 1, N_GRAD_EQN 
                  CALL InterpolateToBoundary( e % U_x(i,:,k,nv), spA % v(:,FRONT), N, e % U_xb(nv,i,k,EFRONT) )
                  CALL InterpolateToBoundary( e % U_x(i,:,k,nv), spA % v(:,BACK) , N, e % U_xb(nv,i,k,EBACK) )
                  CALL InterpolateToBoundary( e % U_y(i,:,k,nv), spA % v(:,FRONT), N, e % U_yb(nv,i,k,EFRONT) )
                  CALL InterpolateToBoundary( e % U_y(i,:,k,nv), spA % v(:,BACK) , N, e % U_yb(nv,i,k,EBACK) )
                  CALL InterpolateToBoundary( e % U_z(i,:,k,nv), spA % v(:,FRONT), N, e % U_zb(nv,i,k,EFRONT) )
                  CALL InterpolateToBoundary( e % U_z(i,:,k,nv), spA % v(:,BACK) , N, e % U_zb(nv,i,k,EBACK) )                                    
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
               DO nv = 1, N_GRAD_EQN 
                  CALL InterpolateToBoundary( e % U_x(i,j,:,nv), spA % v(:,BOTTOM), N, e % U_xb(nv,i,j,EBOTTOM) )
                  CALL InterpolateToBoundary( e % U_x(i,j,:,nv), spA % v(:,TOP)   , N, e % U_xb(nv,i,j,ETOP)    )
                  CALL InterpolateToBoundary( e % U_y(i,j,:,nv), spA % v(:,BOTTOM), N, e % U_yb(nv,i,j,EBOTTOM) )
                  CALL InterpolateToBoundary( e % U_y(i,j,:,nv), spA % v(:,TOP)   , N, e % U_yb(nv,i,j,ETOP)    )
                  CALL InterpolateToBoundary( e % U_z(i,j,:,nv), spA % v(:,BOTTOM), N, e % U_zb(nv,i,j,EBOTTOM) )
                  CALL InterpolateToBoundary( e % U_z(i,j,:,nv), spA % v(:,TOP)   , N, e % U_zb(nv,i,j,ETOP)    )                  
               END DO
            END DO
         END DO 

      END SUBROUTINE ProlongGradientToFaces      
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
   SUBROUTINE DGGradSpaceDerivative( u, uL, uR,  N, D, b, deriv ) 
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
      REAL(KIND=RP), DIMENSION(0:N, N_GRAD_EQN) :: u
      REAL(KIND=RP), DIMENSION(N_GRAD_EQN)      :: uL, uR
      REAL(KIND=RP), DIMENSION(0:N, 0:N)   :: D
      REAL(KIND=RP), DIMENSION(0:N, 2)     :: b
!
!     -----------------
!     Output variables:
!     -----------------
!
      REAL(KIND=RP), DIMENSION(0:N,N_GRAD_EQN), INTENT(OUT) :: deriv
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
      DO k = 1, N_GRAD_EQN
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
     
   END SUBROUTINE DGGradSpaceDerivative   
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
