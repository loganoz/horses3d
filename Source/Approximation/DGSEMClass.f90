
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
      USE DGTimeDerivativeMethods
      
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

         PROCEDURE :: SaveSolutionForRestart
         PROCEDURE :: LoadSolutionForRestart  
            
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
!
!     -------------------------
!     Build the different zones
!     -------------------------
!
      call self % mesh % ConstructZones()
!
!     -----------------------------------------
!     Initialize Spatial discretization methods
!     -----------------------------------------
!
      call Initialize_SpaceAndTimeMethods
      
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
      IF ( ALLOCATED(InviscidMethod) ) DEALLOCATE( InviscidMethod )
      IF ( ALLOCATED(ViscousMethod ) ) DEALLOCATE( ViscousMethod ) 
      
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
         USE ProlongToFacesProcedures
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
!$omp do schedule(runtime)
         DO k = 1, SIZE(self % mesh % elements) 
            CALL ProlongToFaces( self % mesh % elements(k), self % spA )
         END DO
!$omp end do
!
!        -----------------
!        Compute gradients
!        -----------------
!
         if ( flowIsNavierStokes ) then
            CALL DGSpatial_ComputeGradient( self % mesh , self % spA , time , self % externalState , self % externalGradients )
         end if
!
!        -------------------------------------------------------
!        Inviscid Riemann fluxes from the solutions on the faces
!        -------------------------------------------------------
!
         CALL ComputeRiemannFluxes( self, time )
!
!        -----------------------
!        Compute time derivative
!        -----------------------
!
         call TimeDerivative_ComputeQDot( self % mesh , self % spA , time )
!$omp end parallel
!
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
!$omp barrier
!$omp do schedule(runtime)
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
               CALL computeBoundaryFlux(self % mesh % elements(eIDLeft), fIDLeft, time, self % externalState , self % externalGradients)
               
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
!$omp end do
         
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
         REAL(KIND=RP) :: inv_flux(N_EQN)
         REAL(KIND=RP) :: visc_flux(N_EQN)
         INTEGER       :: i,j,ii,jj
                  
         DO j = 0, N
            DO i = 0, N
               CALL iijjIndexes(i,j,N,rotation,ii,jj)                              ! This turns according to the rotation of the elements
               CALL RiemannSolver(QLeft  = eL % QB(:,i ,j ,fIDLeft ), &
                                  QRight = eR % QB(:,ii,jj,fIDright), &
                                  nHat   = eL % geom % normal(:,i,j,fIDLeft), &
                                  flux   = inv_flux)
               CALL ViscousMethod % RiemannSolver( QLeft = eL % QB(:,i,j,fIDLeft) , &
                                                  QRight = eR % QB(:,ii,jj,fIDright) , &
                                                  U_xLeft = eL % U_xb(:,i,j,fIDLeft) , &
                                                  U_yLeft = eL % U_yb(:,i,j,fIDLeft) , &
                                                  U_zLeft = eL % U_zb(:,i,j,fIDLeft) , &
                                                  U_xRight = eR % U_xb(:,ii,jj,fIDRight) , &
                                                  U_yRight = eR % U_yb(:,ii,jj,fIDRight) , &
                                                  U_zRight = eR % U_zb(:,ii,jj,fIDRight) , &
                                                    nHat = eL % geom % normal(:,i,j,fIDLeft) , &
                                                   flux  = visc_flux )
               eL % FStarb(:,i ,j,fIDLeft)  =   (inv_flux - visc_flux ) * eL % geom % scal(i ,j,fIDLeft)
               eR % FStarb(:,ii,jj,fIDright) = -(inv_flux - visc_flux ) * eR % geom % scal(ii,jj,fIdright)
            END DO   
         END DO  
         
      END SUBROUTINE computeElementInterfaceFlux

      SUBROUTINE computeBoundaryFlux(elementOnLeft, faceID, time, externalStateProcedure , externalGradientsProcedure )
      USE ElementClass
      USE DGViscousDiscretization
      USE Physics
      USE BoundaryConditionFunctions
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      TYPE(Element)           :: elementOnLeft
      INTEGER                 :: faceID
      REAL(KIND=RP)           :: time
      EXTERNAL                :: externalStateProcedure
      EXTERNAL                :: externalGradientsProcedure
      
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER                         :: i, j
      INTEGER                         :: N
      REAL(KIND=RP)                   :: bvExt(N_EQN), inv_flux(N_EQN)
      REAL(KIND=RP)                   :: UGradExt(NDIM , N_GRAD_EQN) , visc_flux(N_EQN)
      CHARACTER(LEN=BC_STRING_LENGTH) :: boundaryType
            
      N            = elementOnLeft % N
      boundaryType = elementOnLeft % boundaryType(faceID)
      
      DO j = 0, N
         DO i = 0, N
!
!           Inviscid part
!           -------------
            bvExt = elementOnLeft % Qb(:,i,j,faceID)
            CALL externalStateProcedure( elementOnLeft % geom % xb(:,i,j,faceID), &
                                         time, &
                                         elementOnLeft % geom % normal(:,i,j,faceID), &
                                         bvExt,&
                                         boundaryType )
            CALL RiemannSolver(QLeft  = elementOnLeft % Qb(:,i,j,faceID), &
                               QRight = bvExt, &
                               nHat   = elementOnLeft % geom % normal(:,i,j,faceID), &
                               flux   = inv_flux)
!
!           ViscousPart
!           -----------
            if ( flowIsNavierStokes ) then

            UGradExt(IX,:) = elementOnLeft % U_xb(:,i,j,faceID)
            UGradExt(IY,:) = elementOnLeft % U_yb(:,i,j,faceID)
            UGradExt(IZ,:) = elementOnLeft % U_zb(:,i,j,faceID)

            CALL externalGradientsProcedure(  elementOnLeft % geom % xb(:,i,j,faceID), &
                                              time, &
                                              elementOnLeft % geom % normal(:,i,j,faceID), &
                                              UGradExt,&
                                              boundaryType )

            CALL ViscousMethod % RiemannSolver( QLeft = elementOnLeft % Qb(:,i,j,faceID) , &
                                                QRight = bvExt , &
                                                U_xLeft = elementOnLeft % U_xb(:,i,j,faceID) , &
                                                U_yLeft = elementOnLeft % U_yb(:,i,j,faceID) , &
                                                U_zLeft = elementOnLeft % U_zb(:,i,j,faceID) , &
                                                U_xRight = UGradExt(IX,:) , &
                                                U_yRight = UGradExt(IY,:) , &
                                                U_zRight = UGradExt(IZ,:) , &
                                                   nHat = elementOnLeft % geom % normal(:,i,j,faceID) , &
                                                flux = visc_flux )
            else
               visc_flux = 0.0_RP
            end if

            elementOnLeft % FStarb(:,i,j,faceID) = (inv_flux - visc_flux)*elementOnLeft % geom % scal(i,j,faceID)
         END DO   
      END DO   

      END SUBROUTINE computeBoundaryFlux
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE computeElementInterfaceAverage( eL, fIDLeft, eR, fIDRight, N, rotation)
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
         REAL(KIND=RP) :: UL(N_GRAD_EQN), UR(N_GRAD_EQN)
         REAL(KIND=RP) :: d(N_GRAD_EQN)
         INTEGER       :: i,j,ii,jj
         
         DO j = 0, N
            DO i = 0, N
               CALL iijjIndexes(i,j,N,rotation,ii,jj)
!
!                 --------------
!                 u,v,T averages
!                 --------------
!
               CALL GradientValuesForQ( Q  = eL % QB(:,i,j,fIDLeft), U = UL )
               CALL GradientValuesForQ( Q  = eR % QB(:,ii,jj,fIDright), U = UR )

               d = 0.5_RP*(UL + UR)
               
               eL % Ub ( : , i  , j  , fIDLeft  ) = d
               eR % Ub ( : , ii , jj , fIDright ) = d
                                  
               eL % QB(:,i,j,fIDLeft)    = 0.5_RP * ( eL % QB(:,i,j,fIDLeft) + eR % QB(:,ii,jj,fIDright) )
               eR % QB(:,ii,jj,fIDright) = eL % QB(:,i,j,fIDLeft)
               
            END DO   
         END DO
         
      END SUBROUTINE computeElementInterfaceAverage   
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE computeElementInterfaceGradientAverage( eL, fIDLeft, eR, fIDRight, N, rotation)
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
         REAL(KIND=RP) :: UL(N_GRAD_EQN), UR(N_GRAD_EQN)
         REAL(KIND=RP) :: d(N_GRAD_EQN)
         INTEGER       :: i,j,ii,jj
         
         
         DO j = 0, N
            DO i = 0, N
               CALL iijjIndexes(i,j,N,rotation,ii,jj)                    ! This turns according to the rotation of the elements
!
!                 --------
!                 x values
!                 --------
!
               UL = eL % U_xb(:,i,j,fIDLeft)
               UR = eR % U_xb(:,ii,jj,fIDright)

               d = 0.5_RP*(UL + UR)
               
               eL % U_xb(:,i,j,fIDLeft) = d
               eR % U_xb(:,ii,jj,fIDright) = d
!
!                 --------
!                 y values
!                 --------
!
               UL = eL % U_yb(:,i,j,fIDLeft)
               UR = eR % U_yb(:,ii,jj,fIDright)

               d = 0.5_RP*(UL + UR)
               
               eL % U_yb(:,i,j,fIDLeft) = d
               eR % U_yb(:,ii,jj,fIDright) = d
!
!                 --------
!                 z values
!                 --------
!
               UL = eL % U_zb(:,i,j,fIDLeft)
               UR = eR % U_zb(:,ii,jj,fIDright)

               d = 0.5_RP*(UL + UR)
               
               eL % U_zb(:,i,j,fIDLeft) = d
               eR % U_zb(:,ii,jj,fIDright) = d
               
            END DO   
         END DO
         
      END SUBROUTINE computeElementInterfaceGradientAverage            
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
      TYPE(DGSem), TARGET    :: self
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
      REAL(KIND=RP)         :: Q(N_EQN)
!            
      MaximumEigenvalue = 0.0_RP
      
      N = self % spA % N
!
!     -----------------------------------------------------------
!     TODO:
!     dcsi, deta and dzet have been modified so they work for N=2
!     However, I'm not sure if this modification is OK.
!     Besides, problems are expected for N=0.
!     -----------------------------------------------------------
!
      IF ( N==0 ) THEN 
         PRINT*, "Error in MaximumEigenvalue function (N<1)"    
      ENDIF
        
      dcsi = 1.0_RP / abs( self % spA % xi(1) - self % spA % xi(0) )   
      deta = 1.0_RP / abs( self % spA % eta(1) - self % spA % eta(0) )
      dzet = 1.0_RP / abs( self % spA % zeta(1) - self % spA % zeta(0) )

      DO id = 1, SIZE(self % mesh % elements) 
         DO k = 0, N
            DO j = 0, N
               DO i = 0, N
!
!                 ------------------------------------------------------------
!                 The maximum eigenvalues for a particular state is determined
!                 by the physics.
!                 ------------------------------------------------------------
!
                  Q(1:N_EQN) = self % mesh % elements(id) % Q(i,j,k,1:N_EQN)
                  CALL ComputeEigenvaluesForState( Q , eValues )
!
!                 ----------------------------
!                 Compute contravariant values
!                 ----------------------------
!              
                  lamcsi =  ( self % mesh % elements(id) % geom % jGradXi(IX,i,j,k)   * eValues(IX) + &
        &                     self % mesh % elements(id) % geom % jGradXi(IY,i,j,k)   * eValues(IY) + &
        &                     self % mesh % elements(id) % geom % jGradXi(IZ,i,j,k)   * eValues(IZ) ) * dcsi
        
                  lameta =  ( self % mesh % elements(id) % geom % jGradEta(IX,i,j,k)  * eValues(IX) + &
        &                     self % mesh % elements(id) % geom % jGradEta(IY,i,j,k)  * eValues(IY) + &
        &                     self % mesh % elements(id) % geom % jGradEta(IZ,i,j,k)  * eValues(IZ) ) * deta
        
                  lamzet =  ( self % mesh % elements(id) % geom % jGradZeta(IX,i,j,k) * eValues(IX) + &
        &                     self % mesh % elements(id) % geom % jGradZeta(IY,i,j,k) * eValues(IY) + &
        &                     self % mesh % elements(id) % geom % jGradZeta(IZ,i,j,k) * eValues(IZ) ) * dzet
        
                  jac               = self % mesh % elements(id) % geom % jacobian(i,j,k)
                  MaximumEigenvalue = MAX( MaximumEigenvalue, ABS(lamcsi/jac) + ABS(lameta/jac) + ABS(lamzet/jac) )
               END DO
            END DO
         END DO
      END DO 
      
   END FUNCTION MaximumEigenvalue
      
   END Module DGSEMClass
