
!////////////////////////////////////////////////////////////////////////
!
!      DGSEMClass.f95
!      Created: 2008-07-12 13:38:26 -0400 
!      By: David Kopriva
!
!      Basic class for the discontinuous Galerkin spectral element
!      solution of conservation laws.
!
!      Algorithms:
!         Algorithm 136: DGSEM Class
!         Algorithm 129: ConstructDGSem (Constructor)
!         Algorithm 138: ComputeTimeDerivative (TimeDerivative)
!         Algorithm 137: ComputeRiemannFluxes (EdgeFluxes)
!         Algorithm  35: InterpolateToFineMesh (2DCoarseToFineInterpolation)
!
!      Modified for 3D 6/11/15, 11:32 AM by DAK
!
!////////////////////////////////////////////////////////////////////////
!
      Module DGSEMClass
      
      USE NodalStorageClass
      USE HexMeshClass
      USE PhysicsStorage
      
      IMPLICIT NONE
      
      ABSTRACT INTERFACE
         SUBROUTINE externalStateSubroutine(x,t,nHat,Q,boundaryName)
            USE SMConstants
            REAL(KIND=RP)   , INTENT(IN)    :: x(3), t, nHat(3)
            REAL(KIND=RP)   , INTENT(INOUT) :: Q(:)
            CHARACTER(LEN=*), INTENT(IN)    :: boundaryName
         END SUBROUTINE externalStateSubroutine
         
         SUBROUTINE externalGradientsSubroutine(x,t,nHat,gradU,boundaryName)
            USE SMConstants
            REAL(KIND=RP)   , INTENT(IN)    :: x(3), t, nHat(3)
            REAL(KIND=RP)   , INTENT(INOUT) :: gradU(:,:)
            CHARACTER(LEN=*), INTENT(IN)    :: boundaryName
         END SUBROUTINE externalGradientsSubroutine
      END INTERFACE

      TYPE DGSem
         REAL(KIND=RP)                                         :: maxResidual
         INTEGER                                               :: numberOfTimeSteps
         TYPE(NodalStorage)                                    :: spA
         TYPE(HexMesh)                                         :: mesh
         PROCEDURE(externalStateSubroutine)    , NOPASS, POINTER :: externalState => NULL()
         PROCEDURE(externalGradientsSubroutine), NOPASS, POINTER :: externalGradients => NULL()
!
!        ========         
         CONTAINS
!        ========         
!
         PROCEDURE :: construct => ConstructDGSem
         PROCEDURE :: destruct  => DestructDGSem     
            
      END TYPE DGSem
      
      CONTAINS 
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ConstructDGSem( self, polynomialOrder, meshFileName, &
                                 externalState, externalGradients, success )
      IMPLICIT NONE
!
!     --------------------------
!     Constructor for the class.
!     --------------------------
!
      CLASS(DGSem)     :: self
      INTEGER          :: polynomialOrder
      CHARACTER(LEN=*) :: meshFileName
      LOGICAL          :: success
      EXTERNAL         :: externalState, externalGradients
      
      INTEGER :: k
      INTERFACE
         SUBROUTINE externalState(x,t,nHat,Q,boundaryName)
            USE SMConstants
            REAL(KIND=RP)   , INTENT(IN)    :: x(3), t, nHat(3)
            REAL(KIND=RP)   , INTENT(INOUT) :: Q(:)
            CHARACTER(LEN=*), INTENT(IN)    :: boundaryName
         END SUBROUTINE externalState
         
         SUBROUTINE externalGradients(x,t,nHat,gradU,boundaryName)
            USE SMConstants
            REAL(KIND=RP)   , INTENT(IN)    :: x(3), t, nHat(3)
            REAL(KIND=RP)   , INTENT(INOUT) :: gradU(:,:)
            CHARACTER(LEN=*), INTENT(IN)    :: boundaryName
         END SUBROUTINE externalGradients
      END INTERFACE
      
      CALL self % spA % construct( polynomialOrder )
      CALL self % mesh % constructFromFile( meshfileName, self % spA, success )
      
      IF(.NOT. success) RETURN 
!
!     ------------------------
!     Allocate and zero memory
!     ------------------------
!
      DO k = 1, SIZE(self % mesh % elements) 
         CALL allocateElementStorage( self % mesh % elements(k), polynomialOrder, &
                                      N_EQN, N_GRAD_EQN, flowIsNavierStokes )
      END DO
!
!     -----------------------
!     Set boundary conditions
!     -----------------------
!
      self % externalState     => externalState
      self % externalGradients => externalGradients
      
      call assignBoundaryConditions(self)
      
      END SUBROUTINE ConstructDGSem
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DestructDGSem( self )
      IMPLICIT NONE 
      CLASS(DGSem) :: self
      
      CALL self % spA % destruct()
      CALL DestructMesh( self % mesh )
      self % externalState     => NULL()
      self % externalGradients => NULL()
      
      END SUBROUTINE DestructDGSem
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE SaveSolutionForRestart( self, fUnit ) 
         IMPLICIT NONE
         CLASS(DGSem)     :: self
         INTEGER          :: fUnit
         INTEGER          :: k

         DO k = 1, SIZE(self % mesh % elements) 
            WRITE(fUnit) self % mesh % elements(k) % Q
         END DO

      END SUBROUTINE SaveSolutionForRestart
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE LoadSolutionForRestart( self, fUnit ) 
         IMPLICIT NONE
         CLASS(DGSem)     :: self
         INTEGER          :: fUnit
         INTEGER          :: k

         DO k = 1, SIZE(self % mesh % elements) 
            READ(fUnit) self % mesh % elements(k) % Q
         END DO

      END SUBROUTINE LoadSolutionForRestart
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE assignBoundaryConditions(self)
!
!        ------------------------------------------------------------
!        Assign the boundary condition type to the boundaries through
!        their boundary names
!        ------------------------------------------------------------
!
         USE SharedBCModule
         IMPLICIT NONE  
!
!        ---------
!        Arguments
!        ---------
!
         CLASS(DGSem)     :: self
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER                         :: eID, k
         CHARACTER(LEN=BC_STRING_LENGTH) :: boundaryType, boundaryName
         
         DO eID = 1, SIZE( self % mesh % elements)
            DO k = 1,6
               boundaryName = self % mesh % elements(eID) % boundaryName(k)
               IF ( boundaryName /= emptyBCName )     THEN
                  boundaryType = bcTypeDictionary % stringValueForKey(key             = boundaryName, &
                                                                      requestedLength = BC_STRING_LENGTH)
                  IF( LEN_TRIM(boundaryType) > 0) self % mesh % elements(eID) % boundaryType(k) = boundaryType
               END IF 
                                                                   
            END DO  
         END DO
         
      END SUBROUTINE assignBoundaryConditions
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ComputeTimeDerivative( self, time )
         USE DGTimeDerivativeMethods
         IMPLICIT NONE 
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(DGSem)   :: self
         REAL(KIND=RP) :: time
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER :: k, N
!
         N = self % spA % N
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
!$omp parallel
!$omp do
         DO k = 1, SIZE(self % mesh % elements) 
            CALL ProlongToFaces( self % mesh % elements(k), self % spA )
         END DO
!$omp end do
!
!        -------------------------------------------------------
!        Inviscid Riemann fluxes from the solutions on the faces
!        -------------------------------------------------------
!
!$omp do
         DO k = 1, self % mesh % numberOfFaces
            CALL ComputeRiemannFluxes( self, time )
         END DO
!$omp end do
         
!         IF ( flowIsNavierStokes )     THEN
!!
!!           --------------------------------------
!!           Set up the face Values on each element
!!           --------------------------------------
!!
!            DO k = 1, self%mesh%numberOfEdges 
!               CALL ComputeSolutionRiemannFluxes( self%mesh%edges(k), self%mesh%elements, self%dgS, N, t, ExternalState  )
!            END DO
!!
!!           -----------------------------------
!!           Compute the gradients over the mesh
!!           -----------------------------------
!!
!            DO k = 1, SIZE(self%mesh%elements) 
!               CALL ComputeDGGradient( self%dgS(k)%Q, self%dgS(k)%Ub, self%mesh%elements(k)%geom, N, &
!                                       self%spA%D, self%spA%b, self%dgS(k)%U_x, self%dgS(k)%U_y )
!            END DO
!!
!!           ----------------------------------
!!           Prolong the gradients to the faces
!!           ----------------------------------
!!
!            DO k = 1, SIZE(self%mesh%elements) 
!               CALL ProlongGradientToFaces( self%spA, self%mesh%elements(k)%geom, self%dgS(k) )
!            END DO
!!
!!           -------------------------
!!           Compute gradient averages
!!           -------------------------
!!
!            DO k = 1, self%mesh%numberOfEdges 
!               CALL ComputeGradientAverages( self%mesh%edges(k), self%mesh%elements, self%dgS, N, t, ExternalGradients  )
!            END DO         
!         END IF
!
!        ------------------------
!        Compute time derivatives
!        ------------------------
!
!$omp do
         DO k = 1, SIZE(self % mesh % elements) 
            CALL LocalTimeDerivative( self % mesh % elements(k), self % spA, time )
         END DO
!$omp end do
!$omp end parallel

      END SUBROUTINE ComputeTimeDerivative
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE computeRiemannFluxes( self, time )
         USE Physics
         USE BoundaryConditionFunctions
         IMPLICIT NONE 
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(DGSem)   :: self
         REAL(KIND=RP) :: time
!
!        ---------------
!        Local Variables
!        ---------------
!
         INTEGER       :: faceID
         INTEGER       :: eIDLeft, eIDRight
         INTEGER       :: fIDLeft, fIDright
         INTEGER       :: N
        
         N = self % spA % N
         
         DO faceID = 1, SIZE( self % mesh % faces)
            eIDLeft  = self % mesh % faces(faceID) % elementIDs(1) 
            eIDRight = self % mesh % faces(faceID) % elementIDs(2)
            fIDLeft  = self % mesh % faces(faceID) % elementSide(1)
            
            IF ( eIDRight == HMESH_NONE )     THEN
!
!              -------------
!              Boundary face
!              -------------
!
               CALL computeBoundaryFlux(self % mesh % elements(eIDLeft), fIDLeft, time, self % externalState)
               
            ELSE 
!
!              -------------
!              Interior face
!              -------------
!
               fIDRight =  self % mesh % faces(faceID) % elementSide(2)
               
               CALL computeElementInterfaceFlux(eL = self % mesh % elements(eIDLeft) ,fIDLeft  = fIDLeft, &
                                                eR = self % mesh % elements(eIDRight),fIDRight = fIDright,&
                                                N  = N,                                                   &
                                                rotation = self % mesh % faces(faceID) % rotation)
            END IF 

         END DO           
         
         
      END SUBROUTINE computeRiemannFluxes
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE computeElementInterfaceFlux( eL, fIDLeft, eR, fIDRight, N, rotation)
         USE Physics  
         IMPLICIT NONE  
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(Element) :: eL, eR
         INTEGER       :: fIDLeft, fIdright
         INTEGER       :: rotation
         INTEGER       :: N
!
!        ---------------
!        Local variables
!        ---------------
!
         REAL(KIND=RP) :: flux(N_EQN)
         INTEGER       :: i,j,ii,jj
                  
         SELECT CASE ( rotation )
            CASE( 0 ) 
               DO j = 0, N
                  DO i = 0, N
                     CALL RiemannSolver(QLeft  = eL % QB(:,i,j,fIDLeft), &
                                        QRight = eR % QB(:,i,j,fIDright),&
                                        nHat   = eL % geom % normal(:,i,j,fIDLeft), &
                                        flux   = flux)
                     eL % FStarb(:,i,j,fIDLeft)  =  flux*eL % geom % scal(i,j,fIDLeft)
                     eR % FStarb(:,i,j,fIdright) = -flux*eR % geom % scal(i,j,fIdright)
                  END DO   
               END DO   
            CASE( 1 )
                DO j = 0, N
                  jj = j
                  DO i = 0, N
                     ii = N - i
                     CALL RiemannSolver(QLeft  = eL % QB(:,i ,j,fIDLeft), &
                                        QRight = eR % QB(:,ii,jj,fIDright),&
                                        nHat   = eL % geom % normal(:,i,j,fIDLeft), &
                                        flux   = flux)
                     eL % FStarb(:,i ,j,fIDLeft)  =   flux*eL % geom % scal(i ,j,fIDLeft)
                     eR % FStarb(:,ii,jj,fIdright) = -flux*eR % geom % scal(ii,jj,fIdright)
                  END DO   
               END DO   
           CASE( 2 )
               DO j = 0, N
                  jj = N - j
                  DO i = 0, N
                     ii = N - i
                     CALL RiemannSolver(QLeft  = eL % QB(:,i,j,fIDLeft), &
                                        QRight = eR % QB(:,ii,jj,fIDright),&
                                        nHat   = eL % geom % normal(:,i,j,fIDLeft), &
                                        flux   = flux)
                     eL % FStarb(:,i,j,fIDLeft)    =  flux*eL % geom % scal(i,j,fIDLeft)
                     eR % FStarb(:,ii,jj,fIdright) = -flux*eR % geom % scal(ii,jj,fIdright)
                  END DO   
               END DO   
           CASE( 3 )
                DO j = 0, N
                  ii = j
                  DO i = 0, N
                     jj = N - i
                     CALL RiemannSolver(QLeft  = eL % QB(:,i,j,fIDLeft), &
                                        QRight = eR % QB(:,ii,jj,fIDright),&
                                        nHat   = eL % geom % normal(:,i,j,fIDLeft), &
                                        flux   = flux)
                     eL % FStarb(:,i,j,fIDLeft)  =   flux*eL % geom % scal(i,j,fIDLeft)
                     eR % FStarb(:,ii,jj,fIdright) = -flux*eR % geom % scal(ii,jj,fIdright)
                  END DO   
               END DO   
           CASE DEFAULT 
           PRINT *, "Unknown rotation in element faces"
         END SELECT 
         
      END SUBROUTINE computeElementInterfaceFlux
!
!////////////////////////////////////////////////////////////////////////
!
      REAL(KIND=RP) FUNCTION MaxTimeStep( self, cfl ) 
         IMPLICIT NONE
         TYPE(DGSem)    :: self
         REAL(KIND=RP)  :: cfl
         
         MaxTimeStep  = cfl/MaximumEigenvalue( self )
      
      END FUNCTION MaxTimeStep
!
!////////////////////////////////////////////////////////////////////////
!
   FUNCTION MaximumEigenvalue( self )
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
      USE Physics
      
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      TYPE(DGSem)    :: self
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP)         :: eValues(3)
      REAL(KIND=RP)         :: dcsi, deta, dzet, jac, lamcsi, lamzet, lameta
      INTEGER               :: i, j, k, id
      INTEGER               :: N
      REAL(KIND=RP)         :: MaximumEigenvalue
      EXTERNAL              :: ComputeEigenvaluesForState
!            
      MaximumEigenvalue = 0.0_RP
      
      N = self % spA % N
      
      dcsi = 1.0_RP / self % spA % xi(1)
      deta = 1.0_RP / self % spA % eta(1)
      dzet = 1.0_RP / self % spA % zeta(1)

      DO id = 1, SIZE(self % mesh % elements) 
         DO k = 1, N
            DO j = 1, N
               DO i = 1, N
!
!                 ------------------------------------------------------------
!                 The maximum eigenvalues for a particular state is determined
!                 by the physics.
!                 ------------------------------------------------------------
!
                  CALL ComputeEigenvaluesForState( self % mesh % elements(id) % Q(i,j,k,:), eValues )
!
!                 ----------------------------
!                 Compute contravariant values
!                 ----------------------------
!              
                  lamcsi =  ( self % mesh % elements(id) % geom % jGradXi(1,i,j,k)   * eValues(1) + &
        &                     self % mesh % elements(id) % geom % jGradXi(2,i,j,k)   * eValues(2) + &
        &                     self % mesh % elements(id) % geom % jGradXi(3,i,j,k)   * eValues(3) ) * dcsi
        
                  lameta =  ( self % mesh % elements(id) % geom % jGradEta(1,i,j,k)  * eValues(1) + &
        &                     self % mesh % elements(id) % geom % jGradEta(2,i,j,k)  * eValues(2) + &
        &                     self % mesh % elements(id) % geom % jGradEta(3,i,j,k)  * eValues(3) ) * deta
        
                  lamzet =  ( self % mesh % elements(id) % geom % jGradZeta(1,i,j,k) * eValues(1) + &
        &                     self % mesh % elements(id) % geom % jGradZeta(2,i,j,k) * eValues(2) + &
        &                     self % mesh % elements(id) % geom % jGradZeta(3,i,j,k) * eValues(3) ) * dzet
        
                  jac               = self % mesh % elements(id) % geom % jacobian(i,j,k)
                  MaximumEigenvalue = MAX( MaximumEigenvalue, ABS(lamcsi/jac) + ABS(lameta/jac) + ABS(lamzet/jac) )
               END DO
            END DO
         END DO
      END DO 
      
   END FUNCTION MaximumEigenvalue
      
   END Module DGSEMClass
