
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
!      Modified for 3D       6/11/15, 11:32 AM by DAK
!      Modified for mortars 25/04/17, 18:00    by arueda
!
!////////////////////////////////////////////////////////////////////////
!
      Module DGSEMClass
      
      USE NodalStorageClass
      USE HexMeshClass
      USE PhysicsStorage
      USE SpatialDiscretization
      USE ManufacturedSolutions
      
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
         REAL(KIND=RP)                                           :: maxResidual
         INTEGER                                                 :: numberOfTimeSteps
         INTEGER                                                 :: NDOF                         ! Number of degrees of freedom
         TYPE(NodalStorage), ALLOCATABLE                         :: spA(:,:,:)
         INTEGER           , ALLOCATABLE                         :: Nx(:), Ny(:), Nz(:)
         TYPE(HexMesh)                                           :: mesh
         PROCEDURE(externalStateSubroutine)    , NOPASS, POINTER :: externalState => NULL()
         PROCEDURE(externalGradientsSubroutine), NOPASS, POINTER :: externalGradients => NULL()
         LOGICAL                                                 :: ManufacturedSol = .FALSE.   ! Use manifactured solutions? default .FALSE.
!
!        ========         
         CONTAINS
!        ========         
!
         PROCEDURE :: construct => ConstructDGSem
         PROCEDURE :: destruct  => DestructDGSem   
         
         PROCEDURE :: GetQ
         PROCEDURE :: SetQ
         PROCEDURE :: GetQdot
         
         PROCEDURE :: SaveSolutionForRestart
         PROCEDURE :: LoadSolutionForRestart
   
         procedure :: SetInitialCondition => DGSEM_SetInitialCondition
            
      END TYPE DGSem
      
      
      CONTAINS 
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ConstructDGSem( self, meshFileName, &
                                 externalState, externalGradients, polynomialOrder, Nx_, Ny_, Nz_, success )
      IMPLICIT NONE
!
!     --------------------------
!     Constructor for the class.
!     --------------------------
!
      !-----------------------------------------------------------------
      CLASS(DGSem)                :: self                               !<> Class to be constructed
      CHARACTER(LEN=*)            :: meshFileName                       !<  Name of mesh file
      EXTERNAL                    :: externalState, externalGradients   !<  External procedures that define the BCs
      INTEGER, OPTIONAL           :: polynomialOrder(3)                 !<  Uniform polynomial order
      INTEGER, OPTIONAL, TARGET   :: Nx_(:), Ny_(:), Nz_(:)             !<  Non-uniform polynomial order
      LOGICAL, OPTIONAL           :: success                            !>  Construction finalized correctly?
      !-----------------------------------------------------------------
      INTEGER                     :: i,j,k,el                           ! Counters
      INTEGER, POINTER            :: Nx(:), Ny(:), Nz(:)                ! Orders of every element in mesh (used as pointer to use less space)
      INTEGER                     :: nelem                              ! Number of elements in mesh
      INTEGER                     :: fUnit
      !-----------------------------------------------------------------
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
      !-----------------------------------------------------------------
      
!
!     ---------------------------------------
!     Get polynomial orders for every element
!     ---------------------------------------
!
      IF (PRESENT(Nx_) .AND. PRESENT(Ny_) .AND. PRESENT(Nz_)) THEN
         Nx => Nx_
         Ny => Ny_
         Nz => Nz_
         nelem = SIZE(Nx)
      ELSEIF (PRESENT(polynomialOrder)) THEN
         OPEN(newunit = fUnit, FILE = meshFileName )  
            READ(fUnit,*) k, nelem, k                    ! Here k is used as default reader since this variables are not important now
         CLOSE(fUnit)
         
         ALLOCATE (Nx(nelem),Ny(nelem),Nz(nelem))
         Nx = polynomialOrder(1)
         Ny = polynomialOrder(2)
         Nz = polynomialOrder(3)
      ELSE
         ERROR STOP 'ConstructDGSEM: Polynomial order not specified'
      END IF
      
      ! Now store everything in sem
      IF (ALLOCATED(self % Nx)) DEALLOCATE (self % Nx)
      IF (ALLOCATED(self % Ny)) DEALLOCATE (self % Ny)
      IF (ALLOCATED(self % Nz)) DEALLOCATE (self % Nz)
      ALLOCATE (self % Nx(nelem),self % Ny(nelem),self % Nz(nelem))
      self % Nx = Nx
      self % Ny = Ny
      self % Nz = Nz
!
!     -------------------------------------------------------------
!     Construct the polynomial storage for the elements in the mesh
!     -------------------------------------------------------------
!
      IF (ALLOCATED(self % spa)) DEALLOCATE(self % spa)
      ALLOCATE(self % spa(0:MAXVAL(Nx),0:MAXVAL(Ny),0:MAXVAL(Nz)))
      
      self % NDOF = 0
      DO k=1, nelem
         self % NDOF = self % NDOF + N_EQN * (Nx(k) + 1) * (Ny(k) + 1) * (Nz(k) + 1)
         IF (self % spA(Nx(k),Ny(k),Nz(k)) % Constructed) CYCLE
         CALL self % spA(Nx(k),Ny(k),Nz(k)) % construct( Nx(k),Ny(k),Nz(k) )
      END DO
!
!     ------------------
!     Construct the mesh
!     ------------------
!
      CALL self % mesh % constructFromFile( meshfileName, self % spA, Nx, Ny, Nz,  success )
      IF(.NOT. success) RETURN 
!
!     ------------------------
!     Allocate and zero memory
!     ------------------------
!
      DO k = 1, SIZE(self % mesh % elements) 
         CALL allocateElementStorage( self % mesh % elements(k), Nx(k),Ny(k),Nz(k), &
                                      N_EQN, N_GRAD_EQN, flowIsNavierStokes )
      END DO
!
!     ----------------------------------------------------
!     Get manufactured solution source term (if requested)
!     ----------------------------------------------------
!
      IF (self % ManufacturedSol) THEN
         DO el = 1, SIZE(self % mesh % elements) 
            DO k=0, Nz(el)
               DO j=0, Ny(el)
                  DO i=0, Nx(el)
                     IF (flowIsNavierStokes) THEN
                        CALL ManufacturedSolutionSourceNS(self % mesh % elements(el) % geom % x(:,i,j,k), &
                                                          0._RP, &
                                                          self % mesh % elements(el) % S (i,j,k,:)  )
                     ELSE
                        CALL ManufacturedSolutionSourceEuler(self % mesh % elements(el) % geom % x(:,i,j,k), &
                                                             0._RP, &
                                                             self % mesh % elements(el) % S (i,j,k,:)  )
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END IF
!
!     ------------------------
!     Construct the mortar
!     ------------------------
!
      DO k=1, SIZE(self % mesh % faces)
         IF (self % mesh % faces(k) % FaceType == HMESH_INTERIOR) THEN !The mortar is only needed for the interior edges
            CALL ConstructMortarStorage( self % mesh % faces(k), N_EQN, N_GRAD_EQN, self % mesh % elements )
         END IF
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
      
      NULLIFY(Nx,Ny,Nz)
      
      END SUBROUTINE ConstructDGSem
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DestructDGSem( self )
      IMPLICIT NONE 
      CLASS(DGSem) :: self
      INTEGER      :: i,j,k      !Counter
      
      DO k=0, UBOUND(self % spA,3)
         DO j=0, UBOUND(self % spA,2)
            DO i=0, UBOUND(self % spA,1)
               IF (.NOT. self % spA(i,j,k) % Constructed) CYCLE
               CALL self % spA(i,j,k) % destruct()
            END DO
         END DO
      END DO
      
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
!  Routine to set the solution in each element with a global solution vector
!
   SUBROUTINE SetQ(self,Q)
      IMPLICIT NONE
      CLASS(DGSem)   ,     INTENT(INOUT)           :: self 
      REAL(KIND = RP),     INTENT(IN)              :: Q(:)   
      
      INTEGER                                      :: Nx, Ny, Nz, l, i, j, k, counter, elm
      
      IF (SIZE(Q) /= self % NDOF) ERROR STOP 'Size mismatch in DGSEM:SetQ'
      
      counter = 1
      DO elm = 1, size(self%mesh%elements)
         Nx = self%mesh%elements(elm)%Nxyz(1)
         Ny = self%mesh%elements(elm)%Nxyz(2)
         Nz = self%mesh%elements(elm)%Nxyz(3)
         DO k = 0, Nz
            DO j = 0, Ny
               DO i = 0, Nx
                  DO l = 1,N_EQN
                     self%mesh%elements(elm)%Q(i, j, k, l) = Q(counter) ! This creates a temporary array: storage must be modified to avoid that
                     counter =  counter + 1
                  END DO
               END DO
            END DO
         END DO
      END DO 
         
   END SUBROUTINE SetQ
!
!////////////////////////////////////////////////////////////////////////////////////////       
!
!  Routine to get the solution in each element as a global solution vector
!
   SUBROUTINE GetQ(self,Q)
      IMPLICIT NONE
      CLASS(DGSem),        INTENT(INOUT)            :: self
      REAL(KIND = RP),     INTENT(OUT)              :: Q(:)
      
      INTEGER                                       :: Nx, Ny, Nz, l, i, j, k, counter, elm
      
      IF (SIZE(Q) /= self % NDOF) ERROR STOP 'Size mismatch in DGSEM:GetQ'
      counter = 1
      DO elm = 1, size(self%mesh%elements)
         Nx = self%mesh%elements(elm)%Nxyz(1)
         Ny = self%mesh%elements(elm)%Nxyz(2)
         Nz = self%mesh%elements(elm)%Nxyz(3)
         DO k = 0, Nz
            DO j = 0, Ny
                DO i = 0, Nx
                  DO l = 1,N_EQN
                     Q(counter)  = self%mesh%elements(elm)%Q(i, j, k, l) ! This creates a temporary array: storage must be modified to avoid that
                     counter =  counter + 1
                  END DO
                END DO
            END DO
         END DO
      END DO
      
   END SUBROUTINE GetQ
!
!////////////////////////////////////////////////////////////////////////////////////////      
!
!  Routine to get the solution's time derivative in each element as a global solution vector
!
   SUBROUTINE GetQdot(self,Qdot)
      IMPLICIT NONE
      CLASS(DGSem),        INTENT(INOUT)            :: self
      REAL(KIND = RP),     INTENT(OUT)              :: Qdot(:)
      
      INTEGER                                       :: Nx, Ny, Nz, l, i, j, k, counter, elm
      
      IF (SIZE(Qdot) /= self % NDOF) ERROR STOP 'Size mismatch in DGSEM:GetQdot'
      counter = 1
      DO elm = 1, size(self%mesh%elements)
         Nx = self%mesh%elements(elm)%Nxyz(1)
         Ny = self%mesh%elements(elm)%Nxyz(2)
         Nz = self%mesh%elements(elm)%Nxyz(3)
         DO k = 0, Nz
            DO j = 0, Ny
               DO i = 0, Nx
                  DO l = 1,N_EQN
                     Qdot(counter)  = self%mesh%elements(elm)%Qdot(i, j, k, l) ! This creates a temporary array: storage must be modified to avoid that
                     counter =  counter + 1
                  END DO
               END DO
            END DO
         END DO
      END DO
      
   END SUBROUTINE GetQdot
!
!////////////////////////////////////////////////////////////////////////////////////////      
!
!  -----------------------------------
!  Compute maximum residual L_inf norm
!  -----------------------------------
   FUNCTION ComputeMaxResidual(self) RESULT(maxResidual)
      IMPLICIT NONE
      !----------------------------------------------
      CLASS(DGSem)  :: self
      REAL(KIND=RP) :: maxResidual(N_EQN)
      !----------------------------------------------
      INTEGER       :: id , eq
      REAL(KIND=RP) :: localMaxResidual(N_EQN)
      !----------------------------------------------
      
      maxResidual = 0.0_RP
      DO id = 1, SIZE( self % mesh % elements )
         DO eq = 1 , N_EQN
            localMaxResidual(eq) = MAXVAL(ABS(self % mesh % elements(id) % QDot(:,:,:,eq)))
            maxResidual(eq) = MAX(maxResidual(eq),localMaxResidual(eq))
         END DO
      END DO
   END FUNCTION ComputeMaxResidual
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
         USE SpatialDiscretization
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
         INTEGER :: k, Nx, Ny, Nz
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
!$omp parallel
!$omp do private(Nx,Ny,Nz) schedule(runtime)
         DO k = 1, SIZE(self % mesh % elements) 
            Nx = self%mesh%elements(k)%Nxyz(1)
            Ny = self%mesh%elements(k)%Nxyz(2)
            Nz = self%mesh%elements(k)%Nxyz(3)
            CALL ProlongToFaces( self % mesh % elements(k), self % spA(Nx,Ny,Nz) )
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
         INTEGER       :: fIDLeft

!$omp barrier
!$omp do private(eIDLeft,eIDRight,fIDLeft) schedule(runtime)
         DO faceID = 1, SIZE( self % mesh % faces)
            eIDLeft  = self % mesh % faces(faceID) % elementIDs(1) 
            eIDRight = self % mesh % faces(faceID) % elementIDs(2)
            IF ( eIDRight == HMESH_NONE )     THEN
!
!              -------------
!              Boundary face
!              -------------
!
               fIDLeft  = self % mesh % faces(faceID) % elementSide(1)
               CALL computeBoundaryFlux(self % mesh % elements(eIDLeft), fIDLeft, time, &
                                        self % externalState , self % externalGradients)
               
            ELSE 
!
!              -------------
!              Interior face
!              -------------
!
               CALL computeElementInterfaceFlux ( eL       = self % mesh % elements(eIDLeft)  , &
                                                  eR       = self % mesh % elements(eIDRight) , &
                                                  thisface = self % mesh % faces(faceID)      )
            END IF 

         END DO           
!$omp end do
         
      END SUBROUTINE computeRiemannFluxes
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ComputeGradientAverages( self, time, externalGradientsProcedure )
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
         
         EXTERNAL      :: externalGradientsProcedure
!
!        ---------------
!        Local Variables
!        ---------------
!
         INTEGER       :: faceID
         INTEGER       :: eIDLeft, eIDRight
         INTEGER       :: fIDLeft
         INTEGER       :: N(2)

         REAL(KIND=RP) :: UGradExt(3,N_GRAD_EQN)
         REAL(KIND=RP) :: UL(N_GRAD_EQN), UR(N_GRAD_EQN), d(N_GRAD_EQN)    
         
         INTEGER       :: i, j
         
!$omp do private(eIDLeft,eIDRight,fIDLeft,N,i,j,UGradExt,UL,UR,d) 
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
               N = self % mesh % elements(eIDLeft) % Nxyz (axisMap(:,fIDLeft))
               DO j = 0, N(2)
                  DO i = 0, N(1)
                  
                     UGradExt(1,:) = self % mesh % elements(eIDLeft) % U_xb(:,i,j,fIDLeft)
                     UGradExt(2,:) = self % mesh % elements(eIDLeft) % U_yb(:,i,j,fIDLeft)
                     UGradExt(3,:) = self % mesh % elements(eIDLeft) % U_zb(:,i,j,fIDLeft)
                     
                     CALL externalGradientsProcedure  (self % mesh % elements(eIDLeft) % geom % xb(:,i,j,fIDLeft), &
                                                       time, &
                                                       self % mesh % elements(eIDLeft) % geom % normal(:,i,j,fIDLeft), &
                                                       UGradExt,&
                                                       self % mesh % elements(eIDLeft) % boundaryType(fIDLeft) )
!
!                 --------
!                 x values
!                 --------
!
                     UL = self % mesh % elements(eIDLeft) % U_xb(:,i,j,fIDLeft)
                     UR = UGradExt(1,:)

                     d = 0.5_RP*(UL + UR)

                     self % mesh % elements(eIDLeft) % U_xb(:,i,j,fIDLeft) = d
!
!                 --------
!                 y values
!                 --------
!
                     UL = self % mesh % elements(eIDLeft) % U_yb(:,i,j,fIDLeft)
                     UR = UGradExt(2,:)

                     d = 0.5_RP*(UL + UR)

                     self % mesh % elements(eIDLeft) % U_yb(:,i,j,fIDLeft) = d
!
!                 --------
!                 z values
!                 --------
!
                     UL = self % mesh % elements(eIDLeft) % U_zb(:,i,j,fIDLeft)
                     UR = UGradExt(3,:)

                     d = 0.5_RP*(UL + UR)

                     self % mesh % elements(eIDLeft) % U_zb(:,i,j,fIDLeft) = d

                  END DO   
               END DO   
            
            ELSE 
!
!              -------------
!              Interior face
!              -------------
!
               CALL computeElementInterfaceGradientAverage  ( eL       = self % mesh % elements(eIDLeft)  , &
                                                              eR       = self % mesh % elements(eIDRight) , &
                                                              thisface = self % mesh % faces(faceID)      )
               
            END IF 

         END DO           
!$omp enddo         
         
      END SUBROUTINE computeGradientAverages
!
!//////////////////////////////////////////////////////////////////////// 
! 
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE computeElementInterfaceFlux( eL, eR, thisface)
         IMPLICIT NONE
         !-----------------------------------------
         TYPE(Element), INTENT(INOUT) :: eL, eR       !<> elements
         TYPE(Face)   , INTENT(INOUT) :: thisface     !<> Mortar
         !-----------------------------------------
         INTEGER       :: fIDLeft, fIDRight
         REAL(KIND=RP) :: norm(3), scal
         INTEGER       :: i,j,ii,jj
         INTEGER       :: Nxy(2)       ! Polynomial orders on the interface
         INTEGER       :: NL(2), NR(2)
         INTEGER       :: rotation
!
!        ------------------------------------------------------------------------
!        The following are auxiliar variables required until mortars modification 
!        ------------------------------------------------------------------------
!
         real(kind=RP) :: QL      (1:N_EQN,0:thisface % NPhi(1),0:thisface % NPhi(2))
         real(kind=RP) :: QR      (1:N_EQN,0:thisface % NPhiR(1),0:thisface % NPhiR(2))
         real(kind=RP) :: U_xLeft (1:N_EQN,0:thisface % NPhi(1),0:thisface % NPhi(2))
         real(kind=RP) :: U_xRight(1:N_EQN,0:thisface % NPhiR(1),0:thisface % NPhiR(2))
         real(kind=RP) :: U_yLeft (1:N_EQN,0:thisface % NPhi(1),0:thisface % NPhi(2))
         real(kind=RP) :: U_yRight(1:N_EQN,0:thisface % NPhiR(1),0:thisface % NPhiR(2))
         real(kind=RP) :: U_zLeft (1:N_EQN,0:thisface % NPhi(1),0:thisface % NPhi(2))
         real(kind=RP) :: U_zRight(1:N_EQN,0:thisface % NPhiR(1),0:thisface % NPhiR(2))
         real(kind=RP) :: inv_flux(1:N_EQN,0:thisface % NPhi(1),0:thisface % NPhi(2))
         real(kind=RP) :: visc_flux(1:N_EQN,0:thisface % NPhi(1),0:thisface % NPhi(2))
         
         fIDLeft  = thisface % elementSide(1)
         fIDRight = thisface % elementSide(2)
         Nxy      = thisface % NPhi
         NL       = thisface % NL
         NR       = thisface % NR
         rotation = thisface % rotation         
!
!        ---------------------
!        Projection to mortars
!        ---------------------
!
         call ProjectToMortar(thisface, eL % Qb(:,0:NL(1),0:NL(2),fIDLeft), eR % Qb(:,0:NR(1),0:NR(2),fIDright), N_EQN)
         QL = thisface % Phi % L
         QR = thisface % Phi % R

         if ( flowIsNavierStokes ) then
            call ProjectToMortar(thisface, eL % U_xb(:,0:NL(1),0:NL(2),fIDLeft), eR % U_xb(:,0:NR(1),0:NR(2),fIDRight), N_GRAD_EQN)
            U_xLeft = thisface % Phi % L
            U_xRight = thisface % Phi % R

            call ProjectToMortar(thisface, eL % U_yb(:,0:NL(1),0:NL(2),fIDLeft), eR % U_yb(:,0:NR(1),0:NR(2),fIDRight), N_GRAD_EQN)
            U_yLeft = thisface % Phi % L
            U_yRight = thisface % Phi % R

            call ProjectToMortar(thisface, eL % U_zb(:,0:NL(1),0:NL(2),fIDLeft), eR % U_zb(:,0:NR(1),0:NR(2),fIDRight), N_GRAD_EQN)
            U_zLeft = thisface % Phi % L
            U_zRight = thisface % Phi % R

         else
            U_xLeft = 0.0_RP
            U_yLeft = 0.0_RP
            U_zLeft = 0.0_RP
            U_xRight = 0.0_RP
            U_yRight = 0.0_RP
            U_zRight = 0.0_RP

         end if
!
!        --------------
!        Invscid fluxes: Rotation is not accounted in the Mortar projection
!        --------------
!
         norm = eL % geom % normal(:,1,1,fIDLeft)
         DO j = 0, Nxy(2)
            DO i = 0, Nxy(1)
               CALL iijjIndexes(i,j,Nxy(1),Nxy(2),rotation,ii,jj)
               CALL RiemannSolver(QLeft  = QL(:,i,j), &
                                  QRight = QR(:,ii,jj), &
                                  nHat   = norm, &    ! This works only for flat faces
                                  flux   = inv_flux(:,i,j) )

               CALL ViscousMethod % RiemannSolver( QLeft = QL(:,i,j), &
                                                  QRight = QR(:,ii,jj) , &
                                                  U_xLeft = U_xLeft(:,i,j) , &
                                                  U_yLeft = U_yLeft(:,i,j) , &
                                                  U_zLeft = U_zLeft(:,i,j) , &
                                                  U_xRight = U_xRight(:,ii,jj) , &
                                                  U_yRight = U_yRight(:,ii,jj) , &
                                                  U_zRight = U_zRight(:,ii,jj) , &
                                                   nHat = norm , &
                                                   flux  = visc_flux(:,i,j) )
            END DO   
         END DO  

         thisface % Phi % C = (inv_flux - visc_flux) 
!
!        ---------------------------
!        Return the flux to elements: The sign in eR % FstarB has already been accouted.
!        ---------------------------
!
         call ProjectFluxToElement( thisface , &
                                eL % FStarb(:,0:NL(1),0:NL(2),fIDLeft), & 
                                eR % FStarb(:,0:NR(1),0:NR(2),fIDRight), & 
                                N_EQN ) 
!
!        ------------------------
!        Multiply by the Jacobian: TODO this has to be performed before projection
!        ------------------------
!
         do j = 0 , NL(2)  ;  do i = 0 , NL(1)
            eL % FstarB(:,i,j,fIDLeft) = eL % FstarB(:,i,j,fIDLeft) * eL % geom % scal(i,j,fIDLeft)
         end do            ;  end do

         do j = 0 , NR(2)  ;  do i = 0 , NR(1)
            eR % FstarB(:,i,j,fIDRight) = eR % FstarB(:,i,j,fIDRight) * eR % geom % scal(i,j,fIDRight)
         end do            ;  end do



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
      INTEGER, DIMENSION(2)           :: N
      REAL(KIND=RP)                   :: bvExt(N_EQN), inv_flux(N_EQN)
      REAL(KIND=RP)                   :: UGradExt(NDIM , N_GRAD_EQN) , visc_flux(N_EQN)
      CHARACTER(LEN=BC_STRING_LENGTH) :: boundaryType
            
      N            = elementOnLeft % Nxyz (axisMap(:,faceID))
      boundaryType = elementOnLeft % boundaryType(faceID)
      
      DO j = 0, N(2)
         DO i = 0, N(1)
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
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE computeElementInterfaceGradientAverage( eL, eR, thisface)
         USE Physics  
         IMPLICIT NONE  
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(Element) :: eL, eR        !<> Left and right elements on interface
         TYPE(Face)    :: thisface      !<> Face inbetween
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER       :: i,j,ii,jj
!TODO<<<<<<< HEAD
!         
!         
!         DO j = 0, N
!            DO i = 0, N
!               CALL iijjIndexes(i,j,N,rotation,ii,jj)                    ! This turns according to the rotation of the elements
!=======
         INTEGER       :: fIDLeft, fIdright
         INTEGER       :: rotation
         INTEGER       :: Nxy(2)
         INTEGER       :: NL(2), NR(2)
         
         fIDLeft  = thisface % elementSide(1)
         fIDRight = thisface % elementSide(2)
         Nxy      = thisface % NPhi
         NL       = thisface % NL
         NR       = thisface % NR
         rotation = thisface % rotation
         
!
!        ----------------------
!        Compute interface flux
!        Using BR1 (averages)
!        ----------------------
!
!
!              --------
!              x values
!              --------
         CALL ProjectToMortar(thisface, eL % U_xb(:,0:NL(1),0:NL(2),fIDLeft), eR % U_xb(:,0:NR(1),0:NR(2),fIDright), N_GRAD_EQN)
         
         DO j = 0, Nxy(2)
            DO i = 0, Nxy(1)
               CALL iijjIndexes(i,j,Nxy(1),Nxy(2),rotation,ii,jj)                    ! This turns according to the rotation of the elements
               
               thisface % Phi % Caux(:,i,j) = 0.5_RP* (thisface % Phi % L(1:N_GRAD_EQN,i ,j ) + &
                                                       thisface % Phi % R(1:N_GRAD_EQN,ii,jj) )
               
            END DO   
         END DO 
         
         CALL ProjectToElement(thisface                             , &
                               thisface % Phi % Caux                , &
                               eL % U_xb(:,0:NL(1),0:NL(2),fIDLeft) , &
                               eR % U_xb(:,0:NR(1),0:NR(2),fIDright), &
                               N_GRAD_EQN)
!
!              --------
!              y values
!              --------
!        
         CALL ProjectToMortar(thisface, eL % U_yb(:,0:NL(1),0:NL(2),fIDLeft), eR % U_yb(:,0:NR(1),0:NR(2),fIDright), N_GRAD_EQN) 
         
         DO j = 0, Nxy(2)
            DO i = 0, Nxy(1)
               CALL iijjIndexes(i,j,Nxy(1),Nxy(2),rotation,ii,jj)                    ! This turns according to the rotation of the elements

               thisface % Phi % Caux(:,i,j) = 0.5_RP* (thisface % Phi % L(1:N_GRAD_EQN,i , j) + &
                                                       thisface % Phi % R(1:N_GRAD_EQN,ii,jj) )

            END DO   
         END DO 
         
         CALL ProjectToElement(thisface                             , &
                               thisface % Phi % Caux                , &
                               eL % U_yb(:,0:NL(1),0:NL(2),fIDLeft) , &
                               eR % U_yb(:,0:NR(1),0:NR(2),fIDright), &
                               N_GRAD_EQN)
!
!              --------
!              z values
!              --------
!         
         CALL ProjectToMortar(thisface, eL % U_zb(:,0:NL(1),0:NL(2),fIDLeft), eR % U_zb(:,0:NR(1),0:NR(2),fIDright), N_GRAD_EQN) 
         
         DO j = 0, Nxy(2)
            DO i = 0, Nxy(1)
               CALL iijjIndexes(i,j,Nxy(1),Nxy(2),rotation,ii,jj)                    ! This turns according to the rotation of the elements
               
               thisface % Phi % Caux(:,i,j) = 0.5_RP* (thisface % Phi % L(1:N_GRAD_EQN, i, j) + &
                                                       thisface % Phi % R(1:N_GRAD_EQN,ii,jj) )
               
            END DO   
         END DO   
         
         CALL ProjectToElement(thisface                             , &
                               thisface % Phi % Caux                , &
                               eL % U_zb(:,0:NL(1),0:NL(2),fIDLeft) , &
                               eR % U_zb(:,0:NR(1),0:NR(2),fIDright), &
                               N_GRAD_EQN)

         
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
      REAL(KIND=RP)               :: eValues(3)
      REAL(KIND=RP)               :: dcsi, deta, dzet, jac, lamcsi, lamzet, lameta
      INTEGER                     :: i, j, k, id
      INTEGER                     :: N(3)
      REAL(KIND=RP)               :: MaximumEigenvalue
      EXTERNAL                    :: ComputeEigenvaluesForState
      REAL(KIND=RP)               :: Q(N_EQN)
      TYPE(NodalStorage), POINTER :: spA_p
!            
      MaximumEigenvalue = 0.0_RP
      
!
!     -----------------------------------------------------------
!     TODO:
!     dcsi, deta and dzet have been modified so they work for N=2
!     However, I'm not sure if this modification is OK.
!     Besides, problems are expected for N=0.
!     -----------------------------------------------------------
!
      DO id = 1, SIZE(self % mesh % elements) 
         N = self % mesh % elements(id) % Nxyz
         spA_p => self % spA(N(1),N(2),N(3))
         IF ( ANY(N<1) ) THEN 
            PRINT*, "Error in MaximumEigenvalue function (N<1)"    
         ENDIF         
         
         dcsi = 1.0_RP / abs( spA_p % xi(1)   - spA_p % xi  (0) )   
         deta = 1.0_RP / abs( spA_p % eta(1)  - spA_p % eta (0) )
         dzet = 1.0_RP / abs( spA_p % zeta(1) - spA_p % zeta(0) )
         DO k = 0, N(3)
            DO j = 0, N(2)
               DO i = 0, N(1)
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

   subroutine DGSEM_SetInitialCondition( self, controlVariables ) 
      use FTValueDictionaryClass
      implicit none
      class(DGSEM)   :: self
      class(FTValueDictionary), intent(in)   :: controlVariables
interface
         SUBROUTINE UserDefinedInitialCondition(sem, thermodynamics_, &
                                                        dimensionless_,&
                                                        refValues_)
            USE SMConstants
            use PhysicsStorage
            import DGSEM
            implicit none
            class(DGSEM)                  :: sem
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
         END SUBROUTINE UserDefinedInitialCondition
end interface

      call UserDefinedInitialCondition(self, thermodynamics, &
                                              dimensionless, &
                                                  refValues )

   end subroutine DGSEM_SetInitialCondition
      
   END Module DGSEMClass