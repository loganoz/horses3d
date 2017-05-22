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
      use DGInviscidDiscretization
      use DGViscousDiscretization
      use DGWeakIntegrals
!
!

      class(InviscidMethod_t), allocatable         :: InviscidMethod
      class(ViscousMethod_t), allocatable          :: ViscousMethod
!
!     ========      
      CONTAINS 
!     ========      
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine Initialize_SpaceAndTimeMethods()
         use PhysicsStorage
         implicit none
!
!        TODO: this should be selected from the paramfile options
         allocate( StandardDG_t  :: InviscidMethod )

         if ( flowIsNavierStokes ) then
            allocate( BassiRebay1_t :: ViscousMethod  ) 
   
         else
            allocate( ViscousMethod_t  :: ViscousMethod )

         end if

      end subroutine Initialize_SpaceAndTimeMethods
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine TimeDerivative_ComputeQDot( mesh , spA , t )
         use HexMeshClass
         use ElementClass
         use NodalStorageClass
         use PhysicsStorage
         implicit none
         type(HexMesh)              :: mesh
         type(NodalStorage)         :: spA
         real(kind=RP)              :: t
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: eID , iVar 

         do eID = 1 , size(mesh % elements)
!
!           Perform volume integrals
!           ------------------------            
            call TimeDerivative_VolumetricContribution( mesh % elements(eID) , spA , t)
!
!           Perform surface integrals
!           -------------------------
            call TimeDerivative_FacesContribution( mesh % elements(eID) , spA , t)
!
!           Scale with the Jacobian
!           -----------------------
            do iVar = 1 , N_EQN
               mesh % elements(eID) % QDot(:,:,:,iVar) = mesh % elements(eID) % QDot(:,:,:,iVar) &
                                          / mesh % elements(eID) % geom % jacobian
            end do

         end do

      end subroutine TimeDerivative_ComputeQDot

      subroutine TimeDerivative_VolumetricContribution( e , spA , t )
         use HexMeshClass
         use ElementClass
         use NodalStorageClass
         use PhysicsStorage
         implicit none
         type(Element)      :: e
         type(NodalStorage) :: spA
         real(kind=RP)      :: t
!
!        ---------------
!        Local variables
!        ---------------
!
         real ( kind=RP ) :: inviscidContravariantFlux ( 0:spA % N , 0:spA % N , 0:spA % N , 1 : N_EQN , 1:NDIM ) 
         real ( kind=RP ) :: viscousContravariantFlux  ( 0:spA % N , 0:spA % N , 0:spA % N , 1 : N_EQN , 1:NDIM ) 
         real ( kind=RP ) :: contravariantFlux         ( 0:spA % N , 0:spA % N , 0:spA % N , 1 : N_EQN , 1:NDIM ) 
         integer          :: eID

!
!        Compute inviscid and viscous contravariant fluxes
!        -------------------------------------------------
         call InviscidMethod % ComputeInnerFluxes ( e , spA , inviscidContravariantFlux ) 
         call ViscousMethod  % ComputeInnerFluxes ( e , spA , viscousContravariantFlux  ) 
!
!        Compute the total Navier-Stokes flux
!        ------------------------------------
         if ( flowIsNavierStokes ) then
            contravariantFlux = inviscidContravariantFlux - viscousContravariantFlux
         else
            contravariantFlux = inviscidContravariantFlux
         end if
!
!        Perform the Weak Volume Green integral
!        --------------------------------------
         e % QDot = ScalarWeakIntegrals % StdVolumeGreen ( N_EQN , e , spA , contravariantFlux ) 

      end subroutine TimeDerivative_VolumetricContribution

      subroutine TimeDerivative_FacesContribution( e , spA , t )
         use HexMeshClass
         use NodalStorageClass
         use PhysicsStorage
         implicit none
         type(Element)           :: e
         type(NodalStorage)      :: spA
         real(kind=RP)           :: t
!
!        LEFT face
!        ---------
         e % QDot = e % QDot - ScalarWeakIntegrals % StdFace( e , spA , ELEFT , e % Fstarb(:,:,:,ELEFT) ) 
!
!        RIGHT face
!        ---------
         e % QDot = e % QDot - ScalarWeakIntegrals % StdFace( e , spA , ERIGHT , e % Fstarb(:,:,:,ERIGHT) )
!
!        BOTTOM face
!        ---------
         e % QDot = e % QDot - ScalarWeakIntegrals % StdFace( e , spA , EBOTTOM , e % Fstarb(:,:,:,EBOTTOM) )
!
!        TOP face
!        ---------
         e % QDot = e % QDot - ScalarWeakIntegrals % StdFace( e , spA , ETOP , e % Fstarb(:,:,:,ETOP) )
!
!        BACK face
!        ---------
         e % QDot = e % QDot - ScalarWeakIntegrals % StdFace( e , spA , EBACK , e % Fstarb(:,:,:,EBACK) )
!
!        FRONT face
!        ---------
         e % QDot = e % QDot - ScalarWeakIntegrals % StdFace( e , spA , EFRONT , e % Fstarb(:,:,:,EFRONT) )

      end subroutine TimeDerivative_FacesContribution
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!              GRADIENT PROCEDURES
!              -------------------
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine DGSpatial_ComputeGradient( mesh , spA , time , externalStateProcedure )
         use HexMeshClass
         use NodalStorageClass
         use PhysicsStorage
         implicit none
         type(HexMesh)                  :: mesh
         type(NodalStorage), intent(in) :: spA
         real(kind=RP),      intent(in) :: time
         EXTERNAL                       :: externalStateProcedure

         call ViscousMethod % ComputeGradient( mesh , spA , time , externalStateProcedure)

      end subroutine DGSpatial_ComputeGradient

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
      DO nv = 1, N_EQN
         divFlux(:,:,:,nv) = divFlux(:,:,:,nv)/e % geom % jacobian
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
      REAL(KIND=RP), DIMENSION(0:e % N, 0:e % N, 0:e % N, N_GRAD_EQN) :: U
      REAL(KIND=RP), DIMENSION(0:e % N, N_GRAD_EQN) :: fx
      
      INTEGER :: N, j, i, k, nv

      
      N = e % N

      CALL GradientValuesForQ( Q = e % Q, U = U )
      
      DO nv = 1, N_GRAD_EQN

         f(0:N,0:N,0:N,nv) = U(0:N,0:N,0:N,nv) * e % geom % jGradXi(IX,0:N,0:N,0:N)
         g(0:N,0:N,0:N,nv) = U(0:N,0:N,0:N,nv) * e % geom % jGradEta(IX,0:N,0:N,0:N)
         h(0:N,0:N,0:N,nv) = U(0:N,0:N,0:N,nv) * e % geom % jGradZeta(IX,0:N,0:N,0:N) 

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

      DO nv = 1, N_GRAD_EQN
         e % U_x(0:N,0:N,0:N,nv) = e % U_x(0:N,0:N,0:N,nv) / e % geom % jacobian 
      END DO
      
      DO nv = 1, N_GRAD_EQN

         f(0:N,0:N,0:N,nv) = U(0:N,0:N,0:N,nv) * e % geom % jGradXi(IY,0:N,0:N,0:N)
         g(0:N,0:N,0:N,nv) = U(0:N,0:N,0:N,nv) * e % geom % jGradEta(IY,0:N,0:N,0:N)
         h(0:N,0:N,0:N,nv) = U(0:N,0:N,0:N,nv) * e % geom % jGradZeta(IY,0:N,0:N,0:N) 

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
      DO nv = 1, N_GRAD_EQN
         e % U_y(0:N,0:N,0:N,nv) = e % U_y(0:N,0:N,0:N,nv) / e % geom % jacobian 
      END DO
      
      DO nv = 1, N_GRAD_EQN

         f(0:N,0:N,0:N,nv) = U(0:N,0:N,0:N,nv) * e % geom % jGradXi(IZ,0:N,0:N,0:N)
         g(0:N,0:N,0:N,nv) = U(0:N,0:N,0:N,nv) * e % geom % jGradEta(IZ,0:N,0:N,0:N)
         h(0:N,0:N,0:N,nv) = U(0:N,0:N,0:N,nv) * e % geom % jGradZeta(IZ,0:N,0:N,0:N) 

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

      DO nv = 1, N_GRAD_EQN
         e % U_z(0:N,0:N,0:N,nv) = e % U_z(0:N,0:N,0:N,nv) / e % geom % jacobian 
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
!        Initialization
!        --------------
!
         e % Qb = 0.0_RP
!
!        --------------
!        Left and right
!        --------------
!
         CALL InterpolateToBoundary( e % Q, spA % v(:,LEFT) , N, IX, e % Qb(:,:,:,ELEFT) , N_EQN)
         CALL InterpolateToBoundary( e % Q, spA % v(:,RIGHT), N, IX, e % Qb(:,:,:,ERIGHT), N_EQN)
!
!        --------------
!        Front and back
!        --------------
!
         CALL InterpolateToBoundary( e % Q, spA % v(:,FRONT), N, IY, e % Qb(:,:,:,EFRONT) , N_EQN)
         CALL InterpolateToBoundary( e % Q, spA % v(:,BACK) , N, IY, e % Qb(:,:,:,EBACK) , N_EQN)
!
!        --------------
!        Bottom and Top
!        --------------
!
         CALL InterpolateToBoundary( e % Q, spA % v(:,BOTTOM), N, IZ, e % Qb(:,:,:,EBOTTOM) , N_EQN)
         CALL InterpolateToBoundary( e % Q, spA % v(:,TOP)   , N, IZ, e % Qb(:,:,:,ETOP)    , N_EQN)

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
!        Initialization
!        --------------
!  
         e % U_xb = 0.0_RP
         e % U_yb = 0.0_RP
         e % U_zb = 0.0_RP

!
!        --------------
!        Left and Right
!        --------------
!
         call InterpolateToBoundary( e % U_x , spA % v(:,LEFT ) , N , IX , e % U_xb(:,:,:,ELEFT  ) , N_GRAD_EQN) 
         call InterpolateToBoundary( e % U_x , spA % v(:,RIGHT) , N , IX , e % U_xb(:,:,:,ERIGHT ) , N_GRAD_EQN) 
         call InterpolateToBoundary( e % U_y , spA % v(:,LEFT ) , N , IX , e % U_yb(:,:,:,ELEFT  ) , N_GRAD_EQN) 
         call InterpolateToBoundary( e % U_y , spA % v(:,RIGHT) , N , IX , e % U_yb(:,:,:,ERIGHT ) , N_GRAD_EQN) 
         call InterpolateToBoundary( e % U_z , spA % v(:,LEFT ) , N , IX , e % U_zb(:,:,:,ELEFT  ) , N_GRAD_EQN) 
         call InterpolateToBoundary( e % U_z , spA % v(:,RIGHT) , N , IX , e % U_zb(:,:,:,ERIGHT ) , N_GRAD_EQN) 
!
!        --------------
!        Front and back
!        --------------
!
         CALL InterpolateToBoundary( e % U_x , spA % v(:,FRONT) , N , IY , e % U_xb(:,:,:,EFRONT ) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_x , spA % v(:,BACK)  , N , IY , e % U_xb(:,:,:,EBACK  ) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_y , spA % v(:,FRONT) , N , IY , e % U_yb(:,:,:,EFRONT ) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_y , spA % v(:,BACK)  , N , IY , e % U_yb(:,:,:,EBACK  ) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_z , spA % v(:,FRONT) , N , IY , e % U_zb(:,:,:,EFRONT ) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_z , spA % v(:,BACK)  , N , IY , e % U_zb(:,:,:,EBACK  ) , N_GRAD_EQN )
!
!        --------------
!        Bottom and Top
!        --------------
!
         CALL InterpolateToBoundary( e % U_x, spA % v(:,BOTTOM), N, IZ , e % U_xb(:,:,:,EBOTTOM) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_x, spA % v(:,TOP)   , N, IZ , e % U_xb(:,:,:,ETOP)    , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_y, spA % v(:,BOTTOM), N, IZ , e % U_yb(:,:,:,EBOTTOM) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_y, spA % v(:,TOP)   , N, IZ , e % U_yb(:,:,:,ETOP)    , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_z, spA % v(:,BOTTOM), N, IZ , e % U_zb(:,:,:,EBOTTOM) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_z, spA % v(:,TOP)   , N, IZ , e % U_zb(:,:,:,ETOP)    , N_GRAD_EQN )                  

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
      SUBROUTINE InterpolateToBoundary( u, v, N, which_dim , bValue , NEQ)
!
!     -------------------------------------------------------------
!     Interpolation to the boundary is a dot product for each row or
!     column. Using here the intrinsic Fortran function, without
!     having tested that it is faster or slower than a direct
!     computation for the values of N that we used in the DGSEM.
!     -------------------------------------------------------------
!
         USE SMConstants
         USE Physics
         IMPLICIT NONE
         INTEGER                      , INTENT(IN)  :: N
         real(kind=RP)                , intent(in)  :: u(0:,0:,0:,1:) , v(0:)
         integer                      , intent(in)  :: which_dim
         REAL(KIND=RP)                , INTENT(INOUT) :: bValue(1:,0:,0:)
         integer                      , intent(in)  :: NEQ
!        --------------------------------------------------------------------------------
         integer                                    :: i , j , k , eq

         select case (which_dim)
            case (IX)

               do eq=1,NEQ ; do j = 0,N ; do i = 0,N ; do k=0,N
                  bValue(eq,i,j) = bValue(eq,i,j) + u(k,i,j,eq) * v(k)
               end do ;      end do     ; end do     ; end do

            case (IY)

               do eq=1,NEQ ; do j = 0,N ; do k = 0,N ; do i=0,N
                  bValue(eq,i,j) = bValue(eq,i,j) + u(i,k,j,eq) * v(k)
               end do ;      end do     ; end do     ; end do

            case (IZ)

               do eq=1,NEQ ; do k = 0,N ; do j = 0,N ; do i=0,N
                  bValue(eq,i,j) = bValue(eq,i,j) + u(i,j,k,eq) * v(k)
               end do ;      end do     ; end do     ; end do

         end select
         
      END SUBROUTINE InterpolateToBoundary

   END MODULE DGTimeDerivativeMethods
