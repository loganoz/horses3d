module DGViscousDiscretization
   implicit none
!
!
   private
   public   ViscousMethod_t , BassiRebay1_t , ViscousMethod
!
!
!  *****************************
!  Generic viscous methods class
!  *****************************
!
   type ViscousMethod_t
      contains
         procedure      :: ComputeGradient    => BaseClass_ComputeGradient
         procedure      :: ComputeInnerFluxes => BaseClass_ComputeInnerFluxes
         procedure      :: RiemannSolver      => BaseClass_RiemannSolver
   end type ViscousMethod_t

   type, extends(ViscousMethod_t)   :: BassiRebay1_t
      contains
         procedure      :: ComputeGradient    => BR1_ComputeGradient
         procedure      :: ComputeInnerFluxes => BR1_ComputeInnerFluxes
         procedure      :: RiemannSolver      => BR1_RiemannSolver
   end type BassiRebay1_t

   type, extends(ViscousMethod_t)   :: InteriorPenalty_t
      contains
!         procedure      :: ComputeGradient     => IP_ComputeGradient
!         procedure      :: ComputeInnerFluxes   => IP_ComputeInnerFluxes
   end type InteriorPenalty_t
!
!
   class(ViscousMethod_t), allocatable          :: ViscousMethod

!
!  ========
   contains
!  ========
!
!
!///////////////////////////////////////////////////////////////////////////////////
!
!           BaseClass Procedures
!           --------------------
!///////////////////////////////////////////////////////////////////////////////////
!
      subroutine BaseClass_ComputeGradient( self , mesh , spA , time , externalStateProcedure , externalGradientsProcedure )
         use HexMeshClass
         use NodalStorageClass
         use PhysicsStorage
         use Physics
         implicit none
         class(ViscousMethod_t),    intent(in) :: self
         class(HexMesh)                        :: mesh
         class(NodalStorage),       intent(in) :: spA
         real(kind=RP),             intent(in) :: time
         external                              :: externalStateProcedure
         external                              :: externalGradientsProcedure
!
!        ---------------------------
!        The base class does nothing
!        ---------------------------
!
      end subroutine BaseClass_ComputeGradient

      subroutine BaseClass_ComputeInnerFluxes( self , e , spA , contravariantFlux )
         use ElementClass
         use NodalStorageClass
         use PhysicsStorage
         implicit none
         class(ViscousMethod_t) ,  intent (in)   :: self
         type(Element)                           :: e
         type(NodalStorage)      ,  intent (in)  :: spA
         real(kind=RP)           ,  intent (out) :: contravariantFlux(0:spA % Nx , 0:spA % Ny , 0:spA % Nz , 1:N_EQN, 1:NDIM)
!
!        ---------------------------
!        The base class does nothing
!        ---------------------------
!
         contravariantFlux = 0.0_RP

      end subroutine BaseClass_ComputeInnerFluxes

      subroutine BaseClass_RiemannSolver ( self , QLeft , QRight , U_xLeft , U_yLeft , U_zLeft , U_xRight , U_yRight , U_zRight , nHat , flux )
         use SMConstants
         use PhysicsStorage
         implicit none
         CLASS(ViscousMethod_t)               :: self
         REAL(KIND=RP), DIMENSION(N_EQN)      :: QLeft
         REAL(KIND=RP), DIMENSION(N_EQN)      :: QRight
         REAL(KIND=RP), DIMENSION(N_GRAD_EQN) :: U_xLeft
         REAL(KIND=RP), DIMENSION(N_GRAD_EQN) :: U_yLeft
         REAL(KIND=RP), DIMENSION(N_GRAD_EQN) :: U_zLeft
         REAL(KIND=RP), DIMENSION(N_GRAD_EQN) :: U_xRight
         REAL(KIND=RP), DIMENSION(N_GRAD_EQN) :: U_yRight
         REAL(KIND=RP), DIMENSION(N_GRAD_EQN) :: U_zRight
         REAL(KIND=RP), DIMENSION(NDIM)       :: nHat
         REAL(KIND=RP), DIMENSION(N_EQN)      :: flux
!
!        ---------------------------
!        The base class does nothing
!        ---------------------------
!
         flux = 0.0_RP

      end subroutine BaseClass_RiemannSolver

!
!///////////////////////////////////////////////////////////////////////////////////
!
!           Bassi-Rebay 1 Procedures
!           ------------------------
!///////////////////////////////////////////////////////////////////////////////////
!
      subroutine BR1_ComputeGradient( self , mesh , spA , time , externalStateProcedure , externalGradientsProcedure)
         use HexMeshClass
         use NodalStorageClass
         use PhysicsStorage
         use Physics
         use ProlongToFacesProcedures
         implicit none
         class(BassiRebay1_t), intent(in) :: self
         class(HexMesh)                   :: mesh
         class(NodalStorage),  intent(in) :: spA(0:,0:,0:)
         real(kind=RP),        intent(in) :: time
         external                         :: externalStateProcedure
         external                         :: externalGradientsProcedure
         INTEGER                          :: Nx, Ny, Nz
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                :: eID , fID , dimID , eqID
!
!        Compute the averaged states
!        ---------------------------
         call BR1_ComputeSolutionRiemannSolver( self , mesh , spA , time, externalStateProcedure )
!
!        Perform volume loops
!        --------------------
!$omp barrier
!$omp do schedule(runtime)
         do eID = 1 , size(mesh % elements)
            Nx = mesh % elements(iID) % Nxyz(1)
            Ny = mesh % elements(iID) % Nxyz(2)
            Nz = mesh % elements(iID) % Nxyz(3)
!
!           Add the volumetric integrals
!           ----------------------------
            call BR1_GradientVolumeLoop( self , mesh % elements(eID) , spA(Nx,Ny,Nz) ) 
!
!           Add the surface integrals
!           -------------------------
            call BR1_GradientFaceLoop( self , mesh % elements(eID) , spA(Nx,Ny,Nz) )
!
!           Perform the scaling
!           -------------------               
            do eqID = 1 , N_GRAD_EQN
               mesh % elements(eID) % U_x(:,:,:,eqID) = &
                           mesh % elements(eID) % U_x(:,:,:,eqID) / mesh % elements(eID) % geom % jacobian
               mesh % elements(eID) % U_y(:,:,:,eqID) = &
                           mesh % elements(eID) % U_y(:,:,:,eqID) / mesh % elements(eID) % geom % jacobian
               mesh % elements(eID) % U_z(:,:,:,eqID) = &
                           mesh % elements(eID) % U_z(:,:,:,eqID) / mesh % elements(eID) % geom % jacobian
            end do

            CALL ProlongGradientToFaces( mesh % elements(eID), spA(Nx,Ny,Nz) )

         end do
!$omp end do

      end subroutine BR1_ComputeGradient
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine BR1_GradientVolumeLoop( self , e , spA )
         use ElementClass
         use NodalStorageClass
         use PhysicsStorage
         use Physics
         use DGWeakIntegrals
         implicit none
         class(BassiRebay1_t),   intent(in)     :: self
         class(Element)                         :: e
         class(NodalStorage),    intent(in)     :: spA
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)          :: U(0:spA%Nx , 0:spA%Ny , 0:spA%Nz,1:N_GRAD_EQN)
!
!        Compute gradient variables
!        --------------------------
         CALL GradientValuesForQ( Q = e % Q, U = U )
!
!        Perform the weak integral
!        -------------------------
         call VectorWeakIntegrals % StdVolumeGreen ( N_GRAD_EQN , e , spA , U , e % U_x , e % U_y , e % U_z )

         e % U_x = -e % U_x
         e % U_y = -e % U_y
         e % U_z = -e % U_z

      end subroutine BR1_GradientVolumeLoop
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine BR1_GradientFaceLoop( self, e , spA )
         use ElementClass
         use NodalStorageClass
         use PhysicsStorage
         use Physics
         use DGWeakIntegrals
         implicit none
         class(BassiRebay1_t),   intent(in)  :: self
         class(Element)                      :: e
         class(NodalStorage),    intent(in)  :: spA
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)        :: faceInt_x(0:spA % Nx , 0:spA % Ny , 0:spA % Nz , N_GRAD_EQN )
         real(kind=RP)        :: faceInt_y(0:spA % Nx , 0:spA % Ny , 0:spA % Nz , N_GRAD_EQN )
         real(kind=RP)        :: faceInt_z(0:spA % Nx , 0:spA % Ny , 0:spA % Nz , N_GRAD_EQN )
!
!        LEFT face
!        ---------
         call VectorWeakIntegrals % StdFace( N_GRAD_EQN , e , spA , ELEFT , e % Ub , faceInt_x , faceInt_y , faceInt_z )
         e % U_x = e % U_x + faceInt_x
         e % U_y = e % U_y + faceInt_y
         e % U_z = e % U_z + faceInt_z
!
!        RIGHT face
!        ----------
         call VectorWeakIntegrals % StdFace( N_GRAD_EQN , e , spA , ERIGHT , e % Ub , faceInt_x , faceInt_y , faceInt_z )
         e % U_x = e % U_x + faceInt_x
         e % U_y = e % U_y + faceInt_y
         e % U_z = e % U_z + faceInt_z
!
!        TOP face
!        --------
         call VectorWeakIntegrals % StdFace( N_GRAD_EQN , e , spA , ETOP , e % Ub , faceInt_x , faceInt_y , faceInt_z )
         e % U_x = e % U_x + faceInt_x
         e % U_y = e % U_y + faceInt_y
         e % U_z = e % U_z + faceInt_z
!
!        BOTTOM face
!        -----------
         call VectorWeakIntegrals % StdFace( N_GRAD_EQN , e , spA , EBOTTOM , e % Ub , faceInt_x , faceInt_y , faceInt_z )
         e % U_x = e % U_x + faceInt_x
         e % U_y = e % U_y + faceInt_y
         e % U_z = e % U_z + faceInt_z
!
!        BACK face
!        ---------
         call VectorWeakIntegrals % StdFace( N_GRAD_EQN , e , spA , EBACK , e % Ub , faceInt_x , faceInt_y , faceInt_z )
         e % U_x = e % U_x + faceInt_x
         e % U_y = e % U_y + faceInt_y
         e % U_z = e % U_z + faceInt_z
!
!        FRONT face
!        ----------
         call VectorWeakIntegrals % StdFace( N_GRAD_EQN , e , spA , EFRONT , e % Ub , faceInt_x , faceInt_y , faceInt_z )
         e % U_x = e % U_x + faceInt_x
         e % U_y = e % U_y + faceInt_y
         e % U_z = e % U_z + faceInt_z

      end subroutine BR1_GradientFaceLoop
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine BR1_ComputeSolutionRiemannSolver( self , mesh , spA , time, externalStateProcedure )
         use HexMeshClass
         use NodalStorageClass
         USE Physics
         USE BoundaryConditionFunctions
         IMPLICIT NONE 
!
!        ---------
!        Arguments
!        ---------
!
         type(BassiRebay1_t) :: self
         type(HexMesh)       :: mesh
         type(NodalStorage)  :: spA
         REAL(KIND=RP)       :: time
         EXTERNAL            :: externalStateProcedure
!
!        ---------------
!        Local Variables
!        ---------------
!
         INTEGER       :: faceID
         INTEGER       :: eIDLeft, eIDRight
         INTEGER       :: fIDLeft, fIDright
         INTEGER       :: N(2)
         
         REAL(KIND=RP) :: bvExt(N_EQN), UL(N_GRAD_EQN), UR(N_GRAD_EQN), d(N_GRAD_EQN)     
         
         INTEGER       :: i, j
         
!$omp barrier
!$omp do private(eIDLeft,eIDRight,fIDLEft,N,i,j,bvExt,d) schedule(runtime)
         DO faceID = 1, SIZE(  mesh % faces)
            eIDLeft  =  mesh % faces(faceID) % elementIDs(1) 
            eIDRight =  mesh % faces(faceID) % elementIDs(2)
            
            IF ( eIDRight == HMESH_NONE )     THEN
               fIDLeft  =  mesh % faces(faceID) % elementSide(1)
!
!              -------------
!              Boundary face
!              -------------
!
               N = mesh % elements(eIDLeft) % Nxyz(axisMap(:,fIDLeft) )
               DO j = 0, N(2)
                  DO i = 0, N(1)

                     bvExt =  mesh % elements(eIDLeft) % Qb(:,i,j,fIDLeft)

                     CALL externalStateProcedure( mesh % elements(eIDLeft) % geom % xb(:,i,j,fIDLeft)    , &
                                                  time                                                   , &
                                                  mesh % elements(eIDLeft) % geom % normal(:,i,j,fIDLeft), &
                                                  bvExt                                                  , &
                                                  mesh % elements(eIDLeft) % boundaryType(fIDLeft) )                                                  
!
!                    ---------------
!                    u, v, w, T averages
!                    ---------------
!
                     CALL GradientValuesForQ(  mesh % elements(eIDLeft) % Qb(:,i,j,fIDLeft), UL )
                     CALL GradientValuesForQ( bvExt, UR )

                     d = 0.5_RP*(UL + UR)
               
                     mesh % elements(eIDLeft) % Ub (:,i,j,fIDLeft) = d

                  END DO   
               END DO   
            
            ELSE 
!
!              -------------
!              Interior face
!              -------------
!
               CALL BR1_ComputeElementInterfaceAverage(eL =  mesh % elements(eIDLeft) ,&
                                                       eR =  mesh % elements(eIDRight),&
                                                       thisface = self % mesh % faces(faceID))

            END IF 

         END DO           
!$omp end do
         
      end subroutine BR1_ComputeSolutionRiemannSolver
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE BR1_ComputeElementInterfaceAverage( eL, eR , thisFace )
         USE Physics  
         use ElementClass
         use FaceClass
         IMPLICIT NONE  
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(Element) :: eL, eR
         TYPE(Face)    :: thisFace
!
!        ---------------
!        Local variables
!        ---------------
!
         REAL(KIND=RP) :: UL(N_GRAD_EQN), UR(N_GRAD_EQN)
         INTEGER       :: i,j,ii,jj
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
!        -----------------
!        Project to mortar
!        -----------------
!
         CALL ProjectToMortar(thisface, eL % QB(:,0:NL(1),0:NL(2),fIDLeft), eR % QB(:,0:NR(1),0:NR(2),fIDright), N_EQN)
!
!        ----------------------
!        Compute interface flux
!        Using BR1 (averages)
!        ----------------------
!
         DO j = 0, Nxy(2)
            DO i = 0, Nxy(1)
               CALL iijjIndexes(i,j,Nxy(1),Nxy(2),rotation,ii,jj)                              ! This turns according to the rotation of the elements
!
!              ----------------
!              u,v,w,T averages
!              ----------------
!
               CALL GradientValuesForQ( Q  = thisface % Phi % L(:,i ,j ), U = UL )
               CALL GradientValuesForQ( Q  = thisface % Phi % R(:,ii,jj), U = UR )
               
               thisface % Phi % Caux(:,i,j) = 0.5_RP*(UL + UR)
               
            END DO   
         END DO 
!
!        ------------------------
!        Project back to elements
!        ------------------------
!
         CALL ProjectToElement(thisface                           , &
                               thisface % Phi % Caux              , &
                               eL % Ub(:,0:NL(1),0:NL(2),fIDLeft) , &
                               eR % Ub(:,0:NR(1),0:NR(2),fIDright), &
                               N_GRAD_EQN)
         
      END SUBROUTINE BR1_ComputeElementInterfaceAverage   
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine BR1_ComputeInnerFluxes( self , e , spA , contravariantFlux )
         use ElementClass
         use NodalStorageClass
         use PhysicsStorage
         use Physics
         implicit none
         class(BassiRebay1_t) ,     intent (in) :: self
         type(Element)                          :: e
         type(NodalStorage)      , intent (in)  :: spA
         real(kind=RP)           , intent (out) :: contravariantFlux(0:spA % Nx , 0:spA % Ny , 0:spA % Nz , 1:N_EQN , 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)       :: cartesianFlux(0:spA % Nx , 0:spA % Ny , 0:spA % Nz , 1:N_EQN , 1:NDIM)
         integer             :: nv

         cartesianFlux = ViscousFlux( spA % Nx , spA % Ny , spA % Nz  , e % Q , e % U_x , e % U_y , e % U_z )

         do nv = 1 , N_EQN
         
            contravariantFlux(:,:,:,nv,IX) =    cartesianFlux(:,:,:,nv,IX) * e % geom % jGradXi(IX,:,:,:)  &
                                             +  cartesianFlux(:,:,:,nv,IY) * e % geom % jGradXi(IY,:,:,:)  &
                                             +  cartesianFlux(:,:,:,nv,IZ) * e % geom % jGradXi(IZ,:,:,:)


            contravariantFlux(:,:,:,nv,IY) =    cartesianFlux(:,:,:,nv,IX) * e % geom % jGradEta(IX,:,:,:)  &
                                             +  cartesianFlux(:,:,:,nv,IY) * e % geom % jGradEta(IY,:,:,:)  &
                                             +  cartesianFlux(:,:,:,nv,IZ) * e % geom % jGradEta(IZ,:,:,:)


            contravariantFlux(:,:,:,nv,IZ) =    cartesianFlux(:,:,:,nv,IX) * e % geom % jGradZeta(IX,:,:,:)  &
                                             +  cartesianFlux(:,:,:,nv,IY) * e % geom % jGradZeta(IY,:,:,:)  &
                                             +  cartesianFlux(:,:,:,nv,IZ) * e % geom % jGradZeta(IZ,:,:,:)

         end do

      end subroutine BR1_ComputeInnerFluxes

      function BR1_RiemannSolver ( self , NEQ , QLeft , QRight ) result ( QStar )
         use SMConstants
         use PhysicsStorage
         use Physics
         implicit none
         CLASS(BassiRebay1_t)          :: self
         REAL(KIND=RP), DIMENSION(NEQ) :: QLeft
         REAL(KIND=RP), DIMENSION(NEQ) :: QRight
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: Q(NEQ) 
!
!>       Old implementation: 1st average, then compute
!        ------------------
         QStar   = 0.5_RP * ( QLeft + QRight)

      end subroutine BR1_RiemannSolver

end module DGViscousDiscretization
