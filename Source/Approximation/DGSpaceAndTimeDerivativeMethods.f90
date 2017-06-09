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
!     ========      
      CONTAINS 
!     ========      
!
!////////////////////////////////////////////////////////////////////////
!
<<<<<<< HEAD
      subroutine Initialize_SpaceAndTimeMethods()
         use PhysicsStorage
         implicit none
!
!        TODO: this should be selected from the paramfile options
         allocate( StandardDG_t  :: InviscidMethod )
=======
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
      REAL(KIND=RP), DIMENSION( 0:spA % Nx, &
                                0:spA % Ny, &
                                0:spA % Nz, &
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
      
!
!     ---------------
!     Add source term
!     ---------------
!
      e % Qdot = e % Qdot + e % S
      
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
      
      DO l = 0, e % Nxyz(3)
         DO m = 0, e % Nxyz(2)
            DO n = 0, e % Nxyz(1)
>>>>>>> aruedaDevelopment

         if ( flowIsNavierStokes ) then
            allocate( BassiRebay1_t :: ViscousMethod  ) 
   
         else
            allocate( ViscousMethod_t  :: ViscousMethod )

<<<<<<< HEAD
         end if
=======
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
      TYPE(Element)                                                          :: e
      TYPE(NodalStorage)                                                     :: spA
      REAL(KIND=RP), DIMENSION(0:spA % Nx, 0:spA % Ny, 0:spA % Nz,1:N_EQN,3) :: contravariantFlux
      REAL(KIND=RP), DIMENSION(0:spA % Nx, 0:spA % Ny, 0:spA % Nz,1:N_EQN)   :: divFlux
!
!     ---------------
!     Local variables
!     ---------------
!
      REAL(KIND=RP), DIMENSION(0:spA % Nx, N_EQN) :: fx
      REAL(KIND=RP), DIMENSION(0:spA % Ny, N_EQN) :: gy
      REAL(KIND=RP), DIMENSION(0:spA % Nz, N_EQN) :: hz
      
      
      INTEGER :: Nx, Ny, Nz
      INTEGER :: j, i, k, nv
      
      Nx = spA % Nx
      Ny = spA % Ny
      Nz = spA % Nz
!
!     ---------------------------
!     Xi-derivative contributions
!     ---------------------------
!
      DO k = 0, Nz
         DO j = 0, Ny
            CALL DGSpaceDerivative( contravariantFlux(:,j,k,:,1), &
                                    e % FStarb(:,j,k,ELEFT), e % FStarb(:,j,k,ERIGHT), &
                                    Nx, spA % Dx, spA % bx, fx )
            divFlux(:,j,k,:) = fx
         END DO
      END DO 
!
!     ----------------------------
!     Eta-Derivative Contributions
!     ----------------------------
!
      DO k = 0, Nz
         DO i = 0, Nx
            CALL DGSpaceDerivative( contravariantFlux(i,:,k,:,2), &
                                    e % FStarb(:,i,k,EFRONT), e % FStarb(:,i,k,EBACK), &
                                    Ny, spA % Dy, spA % by, gy )
            divFlux(i,:,k,:) = divFlux(i,:,k,:) + gy
         END DO
      END DO 
!
!     -----------------------------
!     Zeta-Derivative Contributions
!     -----------------------------
!
      DO j = 0, Ny
         DO i = 0, Nx
            CALL DGSpaceDerivative( contravariantFlux(i,j,:,:,3), &
                                    e % FStarb(:,i,j,EBOTTOM), e % FStarb(:,i,j,ETOP), &
                                    Nz, spA % Dz, spA % bz, hz )
            divFlux(i,j,:,:) = divFlux(i,j,:,:) + hz
         END DO
      END DO 
!
!     ---------
!     Finish up
!     ---------
!
      DO k = 0, Nz
         DO j = 0, Ny
            DO i = 0, Nx
               DO nv = 1, N_EQN
                  divFlux(i,j,k,nv) = divFlux(i,j,k,nv)/e % geom % jacobian(i,j,k)
               END DO
            END DO 
         END DO
      END DO
      
>>>>>>> aruedaDevelopment

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
<<<<<<< HEAD
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: eID , iVar 
!
!$omp barrier
!$omp do schedule(runtime)
         do eID = 1 , size(mesh % elements)
!
!           Perform volume integrals
!           ------------------------            
            call TimeDerivative_VolumetricContribution( mesh % elements(eID) , spA , t)
=======

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
      REAL(KIND=RP), DIMENSION(0:spA % Nx, 0:spA % Ny, 0:spA % Nz, N_GRAD_EQN) :: f
      REAL(KIND=RP), DIMENSION(0:spA % Nx, 0:spA % Ny, 0:spA % Nz, N_GRAD_EQN) :: g
      REAL(KIND=RP), DIMENSION(0:spA % Nx, 0:spA % Ny, 0:spA % Nz, N_GRAD_EQN) :: h
      REAL(KIND=RP), DIMENSION(0:spA % Nx, 0:spA % Ny, 0:spA % Nz, N_GRAD_EQN) :: U
      REAL(KIND=RP), DIMENSION(0:spA % Nx, N_GRAD_EQN) :: fx
      REAL(KIND=RP), DIMENSION(0:spA % Ny, N_GRAD_EQN) :: gy
      REAL(KIND=RP), DIMENSION(0:spA % Nz, N_GRAD_EQN) :: hz
      
      INTEGER :: Nx, Ny, Nz
      INTEGER :: j, i, k, nv

      Nx = spA % Nx
      Ny = spA % Ny
      Nz = spA % Nz

      CALL GradientValuesForQ( Q = e % Q, U = U )
!!
!!    ╔═════════╗
!!    ║ Get U_x ║
!!    ╚═════════╝
!!
      DO nv = 1, N_GRAD_EQN

         f(0:Nx,0:Ny,0:Nz,nv) = U(0:Nx,0:Ny,0:Nz,nv) * e % geom % jGradXi  (IX,0:Nx,0:Ny,0:Nz)
         g(0:Nx,0:Ny,0:Nz,nv) = U(0:Nx,0:Ny,0:Nz,nv) * e % geom % jGradEta (IX,0:Nx,0:Ny,0:Nz)
         h(0:Nx,0:Ny,0:Nz,nv) = U(0:Nx,0:Ny,0:Nz,nv) * e % geom % jGradZeta(IX,0:Nx,0:Ny,0:Nz) 

      END DO
!
!     ---------------------------
!     Xi-derivative contributions
!     ---------------------------
!
      DO k = 0, Nz
         DO j = 0, Ny
            CALL DGGradSpaceDerivative( f(:,j,k,:), &
                                    e % geom%normal(1,j,k,ELEFT) * e % geom%scal(j,k, ELEFT) * e % Ub(:,j,k,ELEFT) , &
                                    e % geom%normal(1,j,k,ERIGHT)* e % geom%scal(j,k, ERIGHT)* e % Ub(:,j,k,ERIGHT), &
                                    Nx, spA % Dx, spA % bx, fx )
            e % U_x(:,j,k,:) = fx
         END DO
      END DO 
!
!     ----------------------------
!     Eta-Derivative Contributions
!     ----------------------------
!
      DO k = 0, Nz
         DO i = 0, Nx
            CALL DGGradSpaceDerivative( g(i,:,k,:), &
                                    e % geom%normal(1,i,k,EFRONT)* e % geom%scal(i,k, EFRONT)* e % Ub(:,i,k,EFRONT), &
                                    e % geom%normal(1,i,k,EBACK) * e % geom%scal(i,k, EBACK) * e % Ub(:,i,k,EBACK) , &
                                    Ny, spA % Dy, spA % by, gy )
            e % U_x(i,:,k,:) = e % U_x(i,:,k,:) + gy
         END DO
      END DO 
!
!     -----------------------------
!     Zeta-Derivative Contributions
!     -----------------------------
!
      DO j = 0, Ny
         DO i = 0, Nx
            CALL DGGradSpaceDerivative( h(i,j,:,:), &
                                    e % geom%normal(1,i,j,EBOTTOM)* e % geom%scal(i,j, EBOTTOM)* e % Ub(:,i,j,EBOTTOM), &
                                    e % geom%normal(1,i,j,ETOP)   * e % geom%scal(i,j, ETOP)   * e % Ub(:,i,j,ETOP)   , &
                                    Nz, spA % Dz, spA % bz, hz )
            e % U_x(i,j,:,:) = e % U_x(i,j,:,:) + hz
         END DO
      END DO 
!
!     ---------
!     Finish up
!     ---------
!
      DO nv = 1, N_GRAD_EQN
         e % U_x(0:Nx,0:Ny,0:Nz,nv) = e % U_x(0:Nx,0:Ny,0:Nz,nv) / e % geom % jacobian 
      END DO
!!
!!    ╔═════════╗
!!    ║ Get U_y ║
!!    ╚═════════╝
!!
      DO nv = 1, N_GRAD_EQN

         f(0:Nx,0:Ny,0:Nz,nv) = U(0:Nx,0:Ny,0:Nz,nv) * e % geom % jGradXi  (IY,0:Nx,0:Ny,0:Nz)
         g(0:Nx,0:Ny,0:Nz,nv) = U(0:Nx,0:Ny,0:Nz,nv) * e % geom % jGradEta (IY,0:Nx,0:Ny,0:Nz)
         h(0:Nx,0:Ny,0:Nz,nv) = U(0:Nx,0:Ny,0:Nz,nv) * e % geom % jGradZeta(IY,0:Nx,0:Ny,0:Nz) 

      END DO
!
!     ---------------------------
!     Xi-derivative contributions
!     ---------------------------
!
      DO k = 0, Nz
         DO j = 0, Ny
            CALL DGGradSpaceDerivative( f(:,j,k,:), &
                                    e % geom%normal(2,j,k,ELEFT) * e % geom%scal(j,k, ELEFT) * e % Ub(:,j,k,ELEFT) , &
                                    e % geom%normal(2,j,k,ERIGHT)* e % geom%scal(j,k, ERIGHT)* e % Ub(:,j,k,ERIGHT), &
                                    Nx, spA % Dx, spA % bx, fx )
            e % U_y(:,j,k,:) = fx
         END DO
      END DO 
!
!     ----------------------------
!     Eta-Derivative Contributions
!     ----------------------------
!
      DO k = 0, Nz
         DO i = 0, Nx
            CALL DGGradSpaceDerivative( g(i,:,k,:), &
                                    e % geom%normal(2,i,k,EFRONT)* e % geom%scal(i,k, EFRONT)* e % Ub(:,i,k,EFRONT), &
                                    e % geom%normal(2,i,k,EBACK) * e % geom%scal(i,k, EBACK) * e % Ub(:,i,k,EBACK) , &
                                    Ny, spA % Dy, spA % by, gy )
            e % U_y(i,:,k,:) = e % U_y(i,:,k,:) + gy
         END DO
      END DO 
!
!     -----------------------------
!     Zeta-Derivative Contributions
!     -----------------------------
!
      DO j = 0, Ny
         DO i = 0, Nx
            CALL DGGradSpaceDerivative( h(i,j,:,:), &
                                    e % geom%normal(2,i,j,EBOTTOM)* e % geom%scal(i,j, EBOTTOM)* e % Ub(:,i,j,EBOTTOM), &
                                    e % geom%normal(2,i,j,ETOP)   * e % geom%scal(i,j, ETOP)   * e % Ub(:,i,j,ETOP)   , &
                                    Nz, spA % Dz, spA % bz, hz )
            e % U_y(i,j,:,:) = e % U_y(i,j,:,:) + hz
         END DO
      END DO
!
!     ---------
!     Finish up
!     ---------
!
      DO nv = 1, N_GRAD_EQN
         e % U_y(0:Nx,0:Ny,0:Nz,nv) = e % U_y(0:Nx,0:Ny,0:Nz,nv) / e % geom % jacobian 
      END DO
!!
!!    ╔═════════╗
!!    ║ Get U_z ║
!!    ╚═════════╝
!!
      DO nv = 1, N_GRAD_EQN

         f(0:Nx,0:Ny,0:Nz,nv) = U(0:Nx,0:Ny,0:Nz,nv) * e % geom % jGradXi  (IZ,0:Nx,0:Ny,0:Nz)
         g(0:Nx,0:Ny,0:Nz,nv) = U(0:Nx,0:Ny,0:Nz,nv) * e % geom % jGradEta (IZ,0:Nx,0:Ny,0:Nz)
         h(0:Nx,0:Ny,0:Nz,nv) = U(0:Nx,0:Ny,0:Nz,nv) * e % geom % jGradZeta(IZ,0:Nx,0:Ny,0:Nz) 

      END DO
!
!     ---------------------------
!     Xi-derivative contributions
!     ---------------------------
!
      DO k = 0, Nz
         DO j = 0, Ny
            CALL DGGradSpaceDerivative( f(:,j,k,:), &
                                    e % geom%normal(3,j,k,ELEFT) * e % geom%scal(j,k, ELEFT) * e % Ub(:,j,k,ELEFT), &
                                    e % geom%normal(3,j,k,ERIGHT)* e % geom%scal(j,k, ERIGHT)* e % Ub(:,j,k,ERIGHT), &
                                    Nx, spA % Dx, spA % bx, fx )
            e % U_z(:,j,k,:) = fx
         END DO
      END DO 
!
!     ----------------------------
!     Eta-Derivative Contributions
!     ----------------------------
!
      DO k = 0, Nz
         DO i = 0, Nx
            CALL DGGradSpaceDerivative( g(i,:,k,:), &
                                    e % geom%normal(3,i,k,EFRONT)* e % geom%scal(i,k, EFRONT)* e % Ub(:,i,k,EFRONT), &
                                    e % geom%normal(3,i,k,EBACK) * e % geom%scal(i,k, EBACK) * e % Ub(:,i,k,EBACK), &
                                    Ny, spA % Dy, spA % by, gy )
            e % U_z(i,:,k,:) = e % U_z(i,:,k,:) + gy
         END DO
      END DO 
!
!     -----------------------------
!     Zeta-Derivative Contributions
!     -----------------------------
!
      DO j = 0, Ny
         DO i = 0, Nx
            CALL DGGradSpaceDerivative( h(i,j,:,:), &
                                    e % geom%normal(3,i,j,EBOTTOM)* e % geom%scal(i,j, EBOTTOM)* e % Ub(:,i,j,EBOTTOM), &
                                    e % geom%normal(3,i,j,ETOP)   * e % geom%scal(i,j, ETOP)   * e % Ub(:,i,j,ETOP), &
                                    Nz, spA % Dz, spA % bz, hz )
            e % U_z(i,j,:,:) = e % U_z(i,j,:,:) + hz
         END DO
      END DO 
!
!     ---------
!     Finish up
!     ---------
!
      DO nv = 1, N_GRAD_EQN
         e % U_z(0:Nx,0:Ny,0:Nz,nv) = e % U_z(0:Nx,0:Ny,0:Nz,nv) / e % geom % jacobian 
      END DO
    
      END SUBROUTINE ComputeDGGradient
      
>>>>>>> aruedaDevelopment
!
!           Perform surface integrals
!           -------------------------
            call TimeDerivative_FacesContribution( mesh % elements(eID) , spA , t)
!
<<<<<<< HEAD
!           Scale with the Jacobian
!           -----------------------
            do iVar = 1 , N_EQN
               mesh % elements(eID) % QDot(:,:,:,iVar) = mesh % elements(eID) % QDot(:,:,:,iVar) &
                                          / mesh % elements(eID) % geom % jacobian
            end do
=======
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
      REAL(KIND=RP), DIMENSION( 0:e % Nxyz(1), &
                                0:e % Nxyz(2), &
                                0:e % Nxyz(3), &
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
      DO l = 0, e % Nxyz(3)
         DO m = 0, e % Nxyz(2)
            DO n = 0, e % Nxyz(1)
>>>>>>> aruedaDevelopment

         end do
!$omp end do

      end subroutine TimeDerivative_ComputeQDot

<<<<<<< HEAD
      subroutine TimeDerivative_VolumetricContribution( e , spA , t )
         use HexMeshClass
         use ElementClass
         use NodalStorageClass
         use PhysicsStorage
         implicit none
         type(Element)      :: e
         type(NodalStorage) :: spA
         real(kind=RP)      :: t
=======
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
          DO j = 0, e % Nxyz (axisMap(2,k))
             DO i = 0, e % Nxyz (axisMap(1,k))
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
>>>>>>> aruedaDevelopment
!
!        ---------------
!        Local variables
!        ---------------
!
<<<<<<< HEAD
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
=======
         INTEGER :: Nx, Ny, Nz, i, j, k, nv
         
         Nx = e % Nxyz(1)
         Ny = e % Nxyz(2)
         Nz = e % Nxyz(3)
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
         CALL InterpolateToBoundary( e % Q, spA % vx(:,LEFT) , Nx, Ny, Nz, IX, e % Qb(:,:,:,ELEFT) , N_EQN)
         CALL InterpolateToBoundary( e % Q, spA % vx(:,RIGHT), Nx, Ny, Nz, IX, e % Qb(:,:,:,ERIGHT), N_EQN)
>>>>>>> aruedaDevelopment
!
!        BOTTOM face
!        ---------
         e % QDot = e % QDot - ScalarWeakIntegrals % StdFace( e , spA , EBOTTOM , e % Fstarb(:,:,:,EBOTTOM) )
!
<<<<<<< HEAD
!        TOP face
!        ---------
         e % QDot = e % QDot - ScalarWeakIntegrals % StdFace( e , spA , ETOP , e % Fstarb(:,:,:,ETOP) )
=======
         CALL InterpolateToBoundary( e % Q, spA % vy(:,FRONT), Nx, Ny, Nz, IY, e % Qb(:,:,:,EFRONT) , N_EQN)
         CALL InterpolateToBoundary( e % Q, spA % vy(:,BACK) , Nx, Ny, Nz, IY, e % Qb(:,:,:,EBACK)  , N_EQN)
>>>>>>> aruedaDevelopment
!
!        BACK face
!        ---------
         e % QDot = e % QDot - ScalarWeakIntegrals % StdFace( e , spA , EBACK , e % Fstarb(:,:,:,EBACK) )
!
<<<<<<< HEAD
!        FRONT face
!        ---------
         e % QDot = e % QDot - ScalarWeakIntegrals % StdFace( e , spA , EFRONT , e % Fstarb(:,:,:,EFRONT) )
=======
         CALL InterpolateToBoundary( e % Q, spA % vz(:,BOTTOM), Nx, Ny, Nz, IZ, e % Qb(:,:,:,EBOTTOM) , N_EQN)
         CALL InterpolateToBoundary( e % Q, spA % vz(:,TOP)   , Nx, Ny, Nz, IZ, e % Qb(:,:,:,ETOP)    , N_EQN)
>>>>>>> aruedaDevelopment

      end subroutine TimeDerivative_FacesContribution
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!              GRADIENT PROCEDURES
!              -------------------
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
<<<<<<< HEAD
      subroutine DGSpatial_ComputeGradient( mesh , spA , time , externalStateProcedure , externalGradientsProcedure )
         use HexMeshClass
         use NodalStorageClass
         use PhysicsStorage
         implicit none
         type(HexMesh)                  :: mesh
         type(NodalStorage), intent(in) :: spA
         real(kind=RP),      intent(in) :: time
         EXTERNAL                       :: externalStateProcedure
         EXTERNAL                       :: externalGradientsProcedure

         call ViscousMethod % ComputeGradient( mesh , spA , time , externalStateProcedure , externalGradientsProcedure )
=======
         INTEGER :: Nx, Ny, Nz
         INTEGER :: i, j, k, nv
         
         Nx = e % Nxyz(1)
         Ny = e % Nxyz(2)
         Nz = e % Nxyz(3)
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
         call InterpolateToBoundary( e % U_x , spA % vx(:,LEFT ) , Nx, Ny, Nz , IX , e % U_xb(:,:,:,ELEFT  ) , N_GRAD_EQN) 
         call InterpolateToBoundary( e % U_x , spA % vx(:,RIGHT) , Nx, Ny, Nz , IX , e % U_xb(:,:,:,ERIGHT ) , N_GRAD_EQN) 
         call InterpolateToBoundary( e % U_y , spA % vx(:,LEFT ) , Nx, Ny, Nz , IX , e % U_yb(:,:,:,ELEFT  ) , N_GRAD_EQN) 
         call InterpolateToBoundary( e % U_y , spA % vx(:,RIGHT) , Nx, Ny, Nz , IX , e % U_yb(:,:,:,ERIGHT ) , N_GRAD_EQN) 
         call InterpolateToBoundary( e % U_z , spA % vx(:,LEFT ) , Nx, Ny, Nz , IX , e % U_zb(:,:,:,ELEFT  ) , N_GRAD_EQN) 
         call InterpolateToBoundary( e % U_z , spA % vx(:,RIGHT) , Nx, Ny, Nz , IX , e % U_zb(:,:,:,ERIGHT ) , N_GRAD_EQN) 
!
!        --------------
!        Front and back
!        --------------
!
         CALL InterpolateToBoundary( e % U_x , spA % vy(:,FRONT) , Nx, Ny, Nz , IY , e % U_xb(:,:,:,EFRONT ) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_x , spA % vy(:,BACK)  , Nx, Ny, Nz , IY , e % U_xb(:,:,:,EBACK  ) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_y , spA % vy(:,FRONT) , Nx, Ny, Nz , IY , e % U_yb(:,:,:,EFRONT ) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_y , spA % vy(:,BACK)  , Nx, Ny, Nz , IY , e % U_yb(:,:,:,EBACK  ) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_z , spA % vy(:,FRONT) , Nx, Ny, Nz , IY , e % U_zb(:,:,:,EFRONT ) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_z , spA % vy(:,BACK)  , Nx, Ny, Nz , IY , e % U_zb(:,:,:,EBACK  ) , N_GRAD_EQN )
!
!        --------------
!        Bottom and Top
!        --------------
!
         CALL InterpolateToBoundary( e % U_x, spA % vz(:,BOTTOM), Nx, Ny, Nz , IZ , e % U_xb(:,:,:,EBOTTOM) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_x, spA % vz(:,TOP)   , Nx, Ny, Nz , IZ , e % U_xb(:,:,:,ETOP)    , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_y, spA % vz(:,BOTTOM), Nx, Ny, Nz , IZ , e % U_yb(:,:,:,EBOTTOM) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_y, spA % vz(:,TOP)   , Nx, Ny, Nz , IZ , e % U_yb(:,:,:,ETOP)    , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_z, spA % vz(:,BOTTOM), Nx, Ny, Nz , IZ , e % U_zb(:,:,:,EBOTTOM) , N_GRAD_EQN )
         CALL InterpolateToBoundary( e % U_z, spA % vz(:,TOP)   , Nx, Ny, Nz , IZ , e % U_zb(:,:,:,ETOP)    , N_GRAD_EQN )                  
>>>>>>> aruedaDevelopment

      end subroutine DGSpatial_ComputeGradient
!
<<<<<<< HEAD
!////////////////////////////////////////////////////////////////////////////////////////
!
=======
      USE PhysicsStorage
      USE PolynomialInterpAndDerivsModule
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      INTEGER                                   :: N
      REAL(KIND=RP), DIMENSION(0:N, N_GRAD_EQN) :: u
      REAL(KIND=RP), DIMENSION(N_GRAD_EQN)      :: uL, uR
      REAL(KIND=RP), DIMENSION(0:N, 0:N)        :: D
      REAL(KIND=RP), DIMENSION(0:N, 2)          :: b
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
      SUBROUTINE InterpolateToBoundary( u, v, Nx, Ny, Nz, which_dim , bValue , NEQ)
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
         INTEGER                      , INTENT(IN)    :: Nx, Ny, Nz
         real(kind=RP)                , intent(in)    :: u(0:Nx,0:Ny,0:Nz,1:NEQ) , v(0:)
         integer                      , intent(in)    :: which_dim
         REAL(KIND=RP)                , INTENT(INOUT) :: bValue(1:,0:,0:)
         integer                      , intent(in)    :: NEQ
!        --------------------------------------------------------------------------------
         integer                                    :: i , j , k , eq

         select case (which_dim)
            case (IX)

               do eq = 1 , NEQ
                  do j = 0 , Nz
                     do i = 0 , Ny
                        do k = 0 , Nx

                           bValue(eq,i,j) = bValue(eq,i,j) + u(k,i,j,eq) * v(k)
   
                        end do
                     end do
                  end do
               end do

            case (IY)

                do eq = 1 , NEQ
                  do j = 0 , Nz
                     do k = 0 , Ny
                        do i = 0 , Nx

                           bValue(eq,i,j) = bValue(eq,i,j) + u(i,k,j,eq) * v(k)
   
                        end do
                     end do
                  end do
               end do

            case (IZ)

               do eq = 1 , NEQ
                  do k = 0 , Nz
                     do j = 0 , Ny
                        do i = 0 , Nx

                           bValue(eq,i,j) = bValue(eq,i,j) + u(i,j,k,eq) * v(k)
   
                        end do
                     end do
                  end do
               end do

         end select
         
      END SUBROUTINE InterpolateToBoundary

>>>>>>> aruedaDevelopment
   END MODULE DGTimeDerivativeMethods
