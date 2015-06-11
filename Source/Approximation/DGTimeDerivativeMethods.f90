!
!////////////////////////////////////////////////////////////////////////
!
!      DGTimeDerivativeRoutines.f95
!      Created: 2008-07-13 16:13:12 -0400 
!      By: David Kopriva  
!
!      Algorithms:
!         Algorithm 114: LocalTimeDerivative (MappedDG2DTimeDerivative)
!
!     Plus:
!         MaximumEigenvalue
!         ComputeContravariantFluxes
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE LocalTimeDerivative( t, spA, e )
      
      USE Nodal2DStorageClass
      USE MappedGeometryClass
      USE PDEModule
      USE ElementClass
      IMPLICIT NONE
!
!     -----------------
!     Input parameters:
!     -----------------
!
      REAL(KIND=RP)           :: t
      TYPE(Nodal2DStorage)    :: spA
      TYPE(Element)           :: e
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER                                             :: N, M
      REAL(KIND=RP), DIMENSION(0:spA % N,0:spA % N,N_EQN) :: fFlux
      REAL(KIND=RP), DIMENSION(0:spA % N,0:spA % N,N_EQN) :: gFlux
      TYPE(MappedGeometry)    :: geom
      
      N = spA % N
      M = spA % N
      
      geom = e % geom
      
      CALL ComputeContravariantFluxes( e % Q, geom, fFlux, gFlux, N )
      IF ( flowIsNavierStokes )     THEN
         CALL AddViscousContravariantFluxes( e, geom, fFlux, gFlux, N )
      END IF
      CALL ComputeDGDivergence( fFlux, gFlux, e % FStarb, geom % jacobian, &
                              N, nEqn, spA % D, spA % b, e % Qdot) !QDot saves the divergence
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
      SUBROUTINE ComputeContravariantFluxes( Q, geom, fFlux, gFlux, N )
      
      USE MappedGeometryClass
      USE PDEModule
      IMPLICIT NONE
!
!     -----------------
!     Input parameters:
!     -----------------
!
      INTEGER                                         :: N
      REAL(KIND=RP)       , DIMENSION(0:N,0:N,N_EQN)  :: Q
      REAL(KIND=RP)       , DIMENSION(0:N,0:N, N_EQN) :: fFlux
      REAL(KIND=RP)       , DIMENSION(0:N,0:N, N_EQN) :: gFlux
      TYPE(MappedGeometry)                            :: geom
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER                        :: i, j, M
      REAL(KIND=RP), DIMENSION(nEqn) :: ff, gg
      
      M = N
!
!     --------------------------
!     X-derivative contributions
!     --------------------------
!
      DO j = 0, M
         DO i = 0, N 
            CALL xFlux( Q(i,j,:), ff )
            CALL yFlux( Q(i,j,:), gg )
            fFlux(i,j,:) = geom % Y_eta(i,j)*ff - geom % X_eta(i,j)*gg
         END DO
      END DO
!
!     ----------------------------
!     Eta-Derivative Contributions
!     ----------------------------
!
      DO j = 0, N
         DO i = 0, M 
            CALL xFlux( Q(i,j,:), ff )
            CALL yFlux( Q(i,j,:), gg )
            gFlux(i,j,:) = -geom % Y_xi(i,j)*ff + geom % X_xi(i,j)*gg
         END DO
      END DO
    
      END SUBROUTINE ComputeContravariantFluxes
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE AddViscousContravariantFluxes( e, geom, fFlux, gFlux, N )
      
      USE MappedGeometryClass
      USE PDEModule
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
      INTEGER                                :: i, j, M, k, indx
      REAL(KIND=RP), DIMENSION(N_EQN)        :: ff, gg
      REAL(KIND=RP), DIMENSION(2,N_GRAD_EQN) :: grad
      INTEGER      , DIMENSION(2)            :: horizontalIndices = (/EBOTTOM, ETOP /)
      INTEGER      , DIMENSION(2)            :: verticalIndices   = (/ELEFT, ERIGHT /)
      
      M = N
!
!     ----------------------
!     Interior contributions
!     ----------------------
!
      DO j = 0, M
         DO i = 0, N 
            grad(1,:) = e % U_x(i,j,:)
            grad(2,:) = e % U_y(i,j,:)
            
            CALL xDiffusiveFlux( e % Q(i,j,:), grad, ff )
            CALL yDiffusiveFlux( e % Q(i,j,:), grad, gg )
            
            fFlux(i,j,:) =  fFlux(i,j,:) - ( geom % Y_eta(i,j)*ff - geom % X_eta(i,j)*gg )
            gFlux(i,j,:) =  gFlux(i,j,:) - ( -geom % Y_xi(i,j)*ff + geom % X_xi(i,j) *gg )
         END DO
      END DO
!
!     ----------------------------------------------------------------------
!     Boundary contributions
!     At this point the boundary values of Qb (Storing Ub) and U_xb and U_yb
!     have already been averaged.
!     ----------------------------------------------------------------------
!
       DO k = 1, 2
          indx = horizontalIndices(k)
          DO i = 0, N
            grad(1,:) = e % U_xb(:,i,indx)
            grad(2,:) = e % U_yb(:,i,indx)
            
            CALL xDiffusiveFlux( e % Qb(:,i,indx), grad, ff )
            CALL yDiffusiveFlux( e % Qb(:,i,indx), grad, gg )
            
            e % FStarb(:,i,indx) = e % FStarb(:,i,indx) - &
                                   (ff*geom % normal(i,1,indx) + gg*geom % normal(i,2,indx))*geom % scal(i,indx)
          END DO
       END DO
       
       DO k = 1, 2
          indx = verticalIndices(k)
          DO j = 0, N
            grad(1,:) = e % U_xb(:,j,indx)
            grad(2,:) = e % U_yb(:,j,indx)
            
            CALL xDiffusiveFlux( e % Qb(:,j,indx), grad, ff )
            CALL yDiffusiveFlux( e % Qb(:,j,indx), grad, gg )
            
            e % FStarb(:,j,indx) = e % FStarb(:,j,indx) - &
                                   (ff*geom % normal(j,1,indx) + gg*geom % normal(j,2,indx))*geom % scal(j,indx)
          END DO
       END DO
    
      END SUBROUTINE AddViscousContravariantFluxes
!
!////////////////////////////////////////////////////////////////////////
!
   FUNCTION MaximumEigenvalue( sem )
!
!  -------------------------------------------------------------------
!  Estimate the maximum eigenvalue of the system. This
!  routine computes a heuristic based on the smallest mesh
!  spacing (which goes as 1/N^2) AND the eigenvalues of the
!  particular hyperbolic system being solved. This is not
!  the only way to estimate the eigenvalues, but it works in practice.
!  Other people use variations on this and we make no assertions that
!  it is the only or best way. Other variations look only at the smallest
!  mesh values, other account for differences across the element.
!  -------------------------------------------------------------------
!
      USE DGSEMClass
      
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      TYPE(DGSem)    :: sem
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP)         :: lambdaX, lambdaY, dXi, dEta, dx, dy
      REAL(KIND=RP)         :: eValues(2)
      INTEGER               :: i, j, k
      REAL(KIND=RP)         :: MaximumEigenvalue
!            
      MaximumEigenvalue = 0.0_RP
      
      DO k = 1, SIZE(sem % mesh % elements)
         DO j = 1, sem % spA % N
            dEta = sem % spA % eta(j) - sem % spA % eta(j-1)
            DO i = 1, sem % spA % N 
               dXi = sem % spA % xi(i) - sem % spA % xi(i-1)
               
               CALL ComputeEigenvalues( sem % mesh % elements(k) % Q(i,j,:), eValues )
               
               dx = (dXi*sem % mesh % elements(k) % geom % X_xi(i,j))**2 + (dEta*sem % mesh % elements(k) % geom % X_eta(i,j))**2
               dx = SQRT(dx)
               dy = (dXi*sem % mesh % elements(k) % geom % Y_xi(i,j))**2 + (dEta*sem % mesh % elements(k) % geom % Y_eta(i,j))**2
               dy = SQRT(dy)
               
               lambdaX = ABS(eValues(1)/dx)
               lambdaY = ABS(evalues(2)/dy)
   
               MaximumEigenvalue = MAX( MaximumEigenvalue, lambdaX + lambdaY )
            END DO
         END DO
      END DO 
      
   END FUNCTION MaximumEigenvalue
!
!////////////////////////////////////////////////////////////////////////
!
   FUNCTION EstimateMaximumEigenvalue( sem, ExternalState, ExternalGradients ) RESULT( ev )
      USE DGSEMClass
      IMPLICIT NONE 
!
!-------------------------------------------------------------------
! Use the Power Method to approximate the magnitude of the largest
! eigenvalue of the spatial discretization. Seems rather rough for
! wall BC applied
!-------------------------------------------------------------------
!
!
!     -----------------
!     Input parameters:
!     -----------------
!
      TYPE(DGSem) :: sem
      EXTERNAL    :: ExternalState, ExternalGradients
!
!     ------------------
!     Output parameters:
!     ------------------
!
      REAL(KIND=RP) :: ev
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER :: k, kMax = 30
      INTEGER       :: eID
      REAL(KIND=RP) :: norm
      
      DO k = 1, kMax
         CALL ComputeTimeDerivative( sem, 0.0_RP, ExternalState, ExternalGradients )
         ev   = 0.0_RP
         norm = 0.0_RP
         DO eID = 1, SIZE(sem % mesh % elements)
            ev   = MAX(ev  ,MAXVAL(ABS(sem % mesh % elements(eID) % QDot)))
            norm = MAX(norm,MAXVAL(ABS(sem % mesh % elements(eID) % Q)))
         END DO

         ev = ev/norm
         DO eID = 1, SIZE(sem % mesh % elements)
            sem % mesh % elements(eID) % Q = sem % mesh % elements(eID) % QDot/norm
         END DO 
         PRINT *, k, ev
      END DO

      

   END FUNCTION EstimateMaximumEigenvalue
