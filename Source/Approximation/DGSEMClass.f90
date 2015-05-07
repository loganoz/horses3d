!
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
!     Plus:
!         ExportToTecplot : Output the results in one of tecplot's formats. Can
!                           also be read by the free VisIt program from LLNL.
!         ExportToTecplotI : Interpolate the solution to a fine mesh before 
!                            exporting.
!
!////////////////////////////////////////////////////////////////////////
!
      Module DGSEMClass
      USE Nodal2DStorageClass
      USE QuadMeshClass
      USE DGSolutionStorageClass
      USE PDEModule
      USE BoundaryConditionFunctions
      IMPLICIT NONE 
!
!-------------------------------------------------------------------
! Basic class for the Spectral element solution of conservation laws
!-------------------------------------------------------------------
!
      TYPE DGSem
         TYPE(Nodal2DStorage)                 :: spA
         TYPE(QuadMesh)                       :: mesh
         TYPE(DGSolutionStorage), ALLOCATABLE :: dgS(:)
      END TYPE DGSem
!
!     --------
!     Generics
!     --------
!
      INTERFACE Construct
         MODULE PROCEDURE ConstructDGSem
      END INTERFACE Construct
      
      INTERFACE Destruct
         MODULE PROCEDURE DestructDGSem
      END INTERFACE Destruct
      
      CONTAINS 
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ConstructDGSem( this, N, meshFile )
!
!     --------------------------
!     Constructor for the class.
!     --------------------------
!
      IMPLICIT NONE
      TYPE(DGSem)      :: this
      INTEGER          :: N
      CHARACTER(LEN=*) :: meshFile
      
      INTEGER :: k
      CALL ConstructNodal2DStorage( this%spA, N )
      CALL ConstructMesh_FromFile_( this%mesh, this%spA, meshFile, N )
!
!     ------------------------
!     Allocate and zero memory
!     ------------------------
!
      ALLOCATE( this%dgS(SIZE(this%mesh%elements)) )
      DO k = 1, SIZE(this%mesh%elements) 
         CALL ConstructDGSolutionStorage( this%dgS(k), N, nEqn, N_GRAD_EQN, flowIsNavierStokes )
      END DO
      
      END SUBROUTINE ConstructDGSem
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DestructDGSem( this )
      IMPLICIT NONE 
      TYPE(DGSem)      :: this
      
      INTEGER :: k
      
      CALL DestructNodal2DStorage( this%spA )
      DO k = 1, SIZE(this%mesh%elements) 
         CALL DestructDGSolutionStorage( this%dgS(k) )
      END DO
      DEALLOCATE( this%dgS )
      CALL DestructMesh( this%mesh )
      
      END SUBROUTINE DestructDGSem
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE SaveSolutionForRestart( this, fUnit ) 
         IMPLICIT NONE
         TYPE(DGSem)      :: this
         INTEGER          :: fUnit
         INTEGER          :: k

         DO k = 1, SIZE(this%mesh%elements) 
            CALL SaveDGSolutionStorageToUnit( this%dgS(k), fUnit )
         END DO

      END SUBROUTINE SaveSolutionForRestart
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE LoadSolutionForRestart( this, fUnit ) 
         IMPLICIT NONE
         TYPE(DGSem)      :: this
         INTEGER          :: fUnit
         INTEGER          :: k

         DO k = 1, SIZE(this%mesh%elements) 
            CALL LoadDGSolutionStorageFromUnit( this%dgS(k), fUnit )
         END DO

      END SUBROUTINE LoadSolutionForRestart
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ComputeTimeDerivative( this, t, ExternalState, ExternalGradients )
         IMPLICIT NONE 
         TYPE(DGSem)   :: this
         REAL(KIND=RP) :: t
         EXTERNAL      :: ExternalState, ExternalGradients
         
         INTEGER :: k, N
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
!$omp parallel
!$omp do
         DO k = 1, SIZE(this%mesh%elements) 
            CALL ProlongToFaces( this%spA, this%mesh%elements(k)%geom, this%dgS(k) )
         END DO
!$omp end do
!
!        -------------------------------------------------------
!        Inviscid Riemann fluxes from the solutions on the faces
!        -------------------------------------------------------
!
         N = this%spA%N
!$omp do
         DO k = 1, this%mesh%numberOfEdges
            CALL ComputeRiemannFluxes( this%mesh%edges(k), this%mesh%elements, this%dgS, N, t, ExternalState )
         END DO
!$omp end do
         
         IF ( flowIsNavierStokes )     THEN
!
!           --------------------------------------
!           Set up the face Values on each element
!           --------------------------------------
!
!$omp do
            DO k = 1, this%mesh%numberOfEdges 
               CALL ComputeSolutionRiemannFluxes( this%mesh%edges(k), this%mesh%elements, this%dgS, N, t, ExternalState  )
            END DO
!$omp end do
!
!           -----------------------------------
!           Compute the gradients over the mesh
!           -----------------------------------
!
!$omp do
            DO k = 1, SIZE(this%mesh%elements) 
               CALL ComputeDGGradient( this%dgS(k)%Q, this%dgS(k)%Ub, this%mesh%elements(k)%geom, N, &
                                       this%spA%D, this%spA%b, this%dgS(k)%U_x, this%dgS(k)%U_y )
            END DO
!$omp end do
!
!           ----------------------------------
!           Prolong the gradients to the faces
!           ----------------------------------
!
!$omp do
            DO k = 1, SIZE(this%mesh%elements) 
               CALL ProlongGradientToFaces( this%spA, this%mesh%elements(k)%geom, this%dgS(k) )
            END DO
!$omp end do
!
!           -------------------------
!           Compute gradient averages
!           -------------------------
!
!$omp do
            DO k = 1, this%mesh%numberOfEdges 
               CALL ComputeGradientAverages( this%mesh%edges(k), this%mesh%elements, this%dgS, N, t, ExternalGradients  )
            END DO         
!$omp end do
         END IF
!
!        ------------------------
!        Compute time derivatives
!        ------------------------
!
!$omp do
         DO k = 1, SIZE(this%mesh%elements) 
            CALL LocalTimeDerivative( t, this%spA, this%mesh%elements(k)%geom, this%dgS(k), t )
         END DO
!$omp end do
!$omp end parallel

      END SUBROUTINE ComputeTimeDerivative
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ComputeRiemannFluxes( edg, elements, dgS, NN, t, ExternalState ) 
         IMPLICIT NONE 
         TYPE(Edge)              :: edg
         TYPE(Element)           :: elements(:)
         TYPE(DGSolutionStorage) :: dgS(:)
         INTEGER                 :: NN
         REAL(KIND=RP)           :: t
         EXTERNAL                :: ExternalState
         
         INTEGER                        :: j, n, e1, s1, e2, s2
         REAL(KIND=RP), DIMENSION(nEqn) :: flux     
         
         IF( edg%edgeType == QMESH_INTERIOR )     THEN
!
!           ---------------
!           Interior fluxes
!           ---------------
!
            n = edg%nStart - edg%nInc
            DO j = 0, NN
               e1  = edg%elementIDs(1)
               s1  = ABS(edg%elementSide(1))
               
               e2  = edg%elementIDs(2)
               s2  = ABS(edg%elementSide(2))
               
               CALL RiemannSolver( dgS(e1)%Qb(:,j,s1), dgS(e2)%Qb(:,n,s2), elements(e1)%geom%normal(j,:,s1), flux )
               dgS(e1)%FStarb(:,j,s1) =  flux*elements(e1)%geom%scal(j,s1)
               dgS(e2)%FStarb(:,n,s2) = -flux*elements(e2)%geom%scal(n,s2)
               n = n + edg%nInc
            END DO
            
         ELSE
!
!           ---------------
!           Boundary fluxes
!           ---------------
!
            DO j = 0, NN
               e1 = edg%elementIDs(1)
               s1  = ABS(edg%elementSide(1))
               CALL ComputeBoundaryFlux( t, elements(e1), s1, j, dgS(e1), ExternalState )
            END DO
         END IF 
         
      END SUBROUTINE ComputeRiemannFluxes
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ComputeSolutionRiemannFluxes( edg, elements, dgS, NN, t, ExternalState )
!
!     ------------------------------------------------------------
!     Bassi-Rebay average of the solution values at the boundaries
!     Store average of u,v,T in boundary solution vector, Qb.
!     ------------------------------------------------------------
!
         IMPLICIT NONE
         
         TYPE(Edge)              :: edg
         TYPE(Element)           :: elements(:)
         TYPE(DGSolutionStorage) :: dgS(:)
         INTEGER                 :: NN
         REAL(KIND=RP)           :: t
         EXTERNAL                :: ExternalState
         
         INTEGER                              :: j, n, e1, s1, e2, s2
         REAL(KIND=RP), DIMENSION(N_EQN)      :: bvExt
         REAL(KIND=RP), DIMENSION(N_GRAD_EQN) :: UL, UR
         REAL(KIND=RP)                        :: d(N_GRAD_EQN)
         
         IF( edg%edgeType == QMESH_INTERIOR )     THEN
!
!           --------------------------
!           Interior solution averages
!           --------------------------
!
            n = edg%nStart - edg%nInc
            DO j = 0, NN
               e1 = edg%elementIDs(1)
               s1  = ABS(edg%elementSide(1))
               
               e2 = edg%elementIDs(2)
               s2  = ABS(edg%elementSide(2))
!
!              --------------
!              u,v,T averages
!              --------------
!
               CALL GradientValuesForQ( dgS(e1)%Qb(:,j,s1), UL )
               CALL GradientValuesForQ( dgS(e2)%Qb(:,n,s2), UR )
               
               d = 0.5*(UL + UR)
               
               dgS(e1)%Ub(:,j,s1) = d
               dgS(e2)%Ub(:,n,s2) = d
!
!              -----------------
!              Solution averages
!              -----------------
!
               dgS(e1)%Qb(:,j,s1) = 0.5_RP*( dgS(e1)%Qb(:,j,s1) + dgS(e2)%Qb(:,n,s2) )
               dgS(e2)%Qb(:,n,s2) = dgS(e1)%Qb(:,j,s1)
               
               n = n + edg%nInc
            END DO
            
         ELSE
!
!           -----------------
!           Boundary averages
!           -----------------
!
            DO j = 0, NN
               e1 = edg%elementIDs(1)
               s1  = ABS(edg%elementSide(1))
               
               bvExt  = dgS(e1)%Qb(:,j,s1)
               CALL ExternalState( elements(e1)%geom%xb(1,j,s1), elements(e1)%geom%xb(2,j,s1), t, &
                                    elements(e1)%geom%normal(j,:,s1), bvExt , &
                                    elements(e1)%boundaryName(s1) )
!
!              ---------------
!              u,v, T averages
!              ---------------
!
               CALL GradientValuesForQ( dgS(e1)%Qb(:,j,s1), UL )
               CALL GradientValuesForQ( bvExt, UR )
               
               d = 0.5*(UL + UR)
               
               dgS(e1)%Ub(:,j,s1) = d
!
!              -----------------
!              Solution averages
!              -----------------
!
               dgS(e1)%Qb(:,j,s1) = 0.5_RP*( dgS(e1)%Qb(:,j,s1) + bvExt )
               
            END DO
         END IF 
         
      END SUBROUTINE ComputeSolutionRiemannFluxes
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ComputeGradientAverages( edg, elements, dgS, NN, t, ExternalState )
!
!     --------------------------------------------------------------
!     Bassi-Rebay average of the gradient values at the boundaries
!     Store average of u,v,T in boundary solution vector, U_xb, U_yb.
!     --------------------------------------------------------------
!
         IMPLICIT NONE
         
         TYPE(Edge)              :: edg
         TYPE(Element)           :: elements(:)
         TYPE(DGSolutionStorage) :: dgS(:)
         INTEGER                 :: NN
         REAL(KIND=RP)           :: t
         EXTERNAL                :: ExternalState
         
         INTEGER                              :: j, n, e1, s1, e2, s2
         REAL(KIND=RP), DIMENSION(N_GRAD_EQN) :: UxExt, UyExt
         REAL(KIND=RP), DIMENSION(N_GRAD_EQN) :: UL, UR
         REAL(KIND=RP)                        :: d(N_GRAD_EQN)
         
         IF( edg%edgeType == QMESH_INTERIOR )     THEN
!
!           --------------------------
!           Interior gradient averages
!           --------------------------
!
            n = edg%nStart - edg%nInc
            DO j = 0, NN
               e1 = edg%elementIDs(1)
               s1  = ABS(edg%elementSide(1))
               
               e2 = edg%elementIDs(2)
               s2  = ABS(edg%elementSide(2))
!
!              --------
!              x values
!              --------
!
               UL = dgS(e1)%U_xb(:,j,s1)
               UR = dgS(e2)%U_xb(:,n,s2)
               
               d = 0.5*(UL + UR)
               
               dgS(e1)%U_xb(:,j,s1) = d
               dgS(e2)%U_xb(:,n,s2) = d
!
!              --------
!              y values
!              --------
!
               UL = dgS(e1)%U_yb(:,j,s1)
               UR = dgS(e2)%U_yb(:,n,s2)
               
               d = 0.5*(UL + UR)
               
               dgS(e1)%U_yb(:,j,s1) = d
               dgS(e2)%U_yb(:,n,s2) = d
               
               n = n + edg%nInc
            END DO
            
         ELSE
!
!           -----------------
!           Boundary averages
!           -----------------
!
            DO j = 0, NN
               e1 = edg%elementIDs(1)
               s1  = ABS(edg%elementSide(1))
               
               UxExt  = dgS(e1)%U_xb(:,j,s1)
               UyExt  = dgS(e1)%U_yb(:,j,s1)
               CALL ExternalState( elements(e1)%geom%xb(1,j,s1), elements(e1)%geom%xb(2,j,s1), t, &
                                   elements(e1)%geom%normal(j,:,s1), UxExt, UyExt , &
                                    elements(e1)%boundaryName(s1) )
!
!              --------
!              x values
!              --------
!
               UL = dgS(e1)%U_xb(:,j,s1)
               UR = UxExt
               
               d = 0.5*(UL + UR)
               
               dgS(e1)%U_xb(:,j,s1) = d
!
!              --------
!              y values
!              --------
!
               UL = dgS(e1)%U_yb(:,j,s1)
               UR = UyExt
               
               d = 0.5*(UL + UR)
               
               dgS(e1)%U_yb(:,j,s1) = d
               
            END DO
         END IF 
         
      END SUBROUTINE ComputeGradientAverages
      
   END Module DGSEMClass
