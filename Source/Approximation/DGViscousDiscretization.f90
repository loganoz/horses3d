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
         real(kind=RP)           ,  intent (out) :: contravariantFlux(0:spA % N , 0:spA % N , 0:spA % N , 1:N_EQN, 1:NDIM)
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
         class(NodalStorage),  intent(in) :: spA
         real(kind=RP),        intent(in) :: time
         external                         :: externalStateProcedure
         external                         :: externalGradientsProcedure
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
!
!           Add the volumetric integrals
!           ----------------------------
            call BR1_GradientVolumeLoop( self , mesh % elements(eID) , spA ) 
!
!           Add the surface integrals
!           -------------------------
            call BR1_GradientFaceLoop( self , mesh % elements(eID) , spA )
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

            CALL ProlongGradientToFaces( mesh % elements(eID), spA )

         end do
!$omp end do

      end subroutine BR1_ComputeGradient

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
         real(kind=RP)          :: U(0:spA%N , 0:spA%N , 0:spA%N,1:N_GRAD_EQN)
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
         real(kind=RP)        :: faceInt_x(0:spA % N , 0:spA % N , 0:spA % N , N_GRAD_EQN )
         real(kind=RP)        :: faceInt_y(0:spA % N , 0:spA % N , 0:spA % N , N_GRAD_EQN )
         real(kind=RP)        :: faceInt_z(0:spA % N , 0:spA % N , 0:spA % N , N_GRAD_EQN )
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
         INTEGER       :: N
         
         REAL(KIND=RP) :: bvExt(N_EQN), UL(N_GRAD_EQN), UR(N_GRAD_EQN), d(N_GRAD_EQN)     
         
         INTEGER       :: i, j
        
         N = spA % N
!$omp barrier
!$omp do schedule(runtime)
         DO faceID = 1, SIZE(  mesh % faces)
            eIDLeft  =  mesh % faces(faceID) % elementIDs(1) 
            eIDRight =  mesh % faces(faceID) % elementIDs(2)
            fIDLeft  =  mesh % faces(faceID) % elementSide(1)
            
            IF ( eIDRight == HMESH_NONE )     THEN
!
!              -------------
!              Boundary face
!              -------------
!
               DO j = 0, N
                  DO i = 0, N

                     bvExt =  mesh % elements(eIDLeft) % Qb(:,i,j,fIDLeft)

                     CALL externalStateProcedure(  mesh % elements(eIDLeft) % geom % xb(:,i,j,fIDLeft), &
                                                  time, &
                                                   mesh % elements(eIDLeft) % geom % normal(:,i,j,fIDLeft), &
                                                  bvExt,&
                                                   mesh % elements(eIDLeft) % boundaryType(fIDLeft) )                                                  
!
!              ---------------
!              u,v, T averages
!              ---------------
!
                     CALL GradientValuesForQ(  mesh % elements(eIDLeft) % Qb(:,i,j,fIDLeft), UL )
                     CALL GradientValuesForQ( bvExt, UR )

                     d = 0.5_RP*(UL + UR)
               
                      mesh % elements(eIDLeft) % Ub (:,i,j,fIDLeft) = d
!
!              -----------------
!              Solution averages
!              -----------------
!                             Qb + bvExt
!TODO: this makes that Qb <= ------------       (Remove?)
!                                2
!
!                     CALL DiffusionRiemannSolution( self % mesh % elements(eIDLeft) % geom % normal(:,i,j,fIDLeft), &
!                                                    self % mesh % elements(eIDLeft) % Qb(:,i,j,fIDLeft), &
!                                                    bvExt, &
!                                                    self % mesh % elements(eIDLeft) % Qb(:,i,j,fIDLeft) )
                                                    
                     !self % mesh % elements(eIDLeft) % Qb(:,i,j,fIDLeft)    = &
                     !& 0.5_RP*( self % mesh % elements(eIDLeft) % Qb(:,i,j,fIDLeft) + bvExt )

                  END DO   
               END DO   
            
            ELSE 
!
!              -------------
!              Interior face
!              -------------
!
               fIDRight =   mesh % faces(faceID) % elementSide(2)
               
               CALL BR1_ComputeElementInterfaceAverage(eL =  mesh % elements(eIDLeft) ,fIDLeft  = fIDLeft, &
                                                       eR =  mesh % elements(eIDRight),fIDRight = fIDright,&
                                                       N  = N,                                             &
                                                 rotation =  mesh % faces(faceID) % rotation)

            END IF 

         END DO           
!$omp end do
         
      end subroutine BR1_ComputeSolutionRiemannSolver

      SUBROUTINE BR1_ComputeElementInterfaceAverage( eL, fIDLeft, eR, fIDRight, N, rotation)
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
               
            END DO   
         END DO
         
      END SUBROUTINE BR1_ComputeElementInterfaceAverage   

      subroutine BR1_ComputeInnerFluxes( self , e , spA , contravariantFlux )
         use ElementClass
         use NodalStorageClass
         use PhysicsStorage
         use Physics
         implicit none
         class(BassiRebay1_t) ,     intent (in) :: self
         type(Element)                          :: e
         type(NodalStorage)      , intent (in)  :: spA
         real(kind=RP)           , intent (out) :: contravariantFlux(0:spA % N , 0:spA % N , 0:spA % N , 1:N_EQN , 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)       :: cartesianFlux(0:spA % N , 0:spA % N , 0:spA % N , 1:N_EQN , 1:NDIM)
         integer             :: nv

         cartesianFlux = ViscousFlux( spA % N , e % Q , e % U_x , e % U_y , e % U_z )

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

      subroutine BR1_RiemannSolver ( self , QLeft , QRight , U_xLeft , U_yLeft , U_zLeft , U_xRight , U_yRight , U_zRight , nHat , flux )
         use SMConstants
         use PhysicsStorage
         use Physics
         implicit none
         CLASS(BassiRebay1_t)                 :: self
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
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: Q(NCONS) , U_x(N_GRAD_EQN) , U_y(N_GRAD_EQN) , U_z(N_GRAD_EQN)
         real(kind=RP)     :: flux_vec(NCONS,NDIM)

!
!>       Old implementation: 1st average, then compute
!        ------------------
         Q   = 0.5_RP * ( QLeft + QRight)
         U_x = 0.5_RP * ( U_xLeft + U_xRight) 
         U_y = 0.5_RP * ( U_yLeft + U_yRight) 
         U_z = 0.5_RP * ( U_zLeft + U_zRight) 

         flux_vec = ViscousFlux(Q,U_x,U_y,U_z)

         flux = flux_vec(:,IX) * nHat(IX) + flux_vec(:,IY) * nHat(IY) + flux_vec(:,IZ) * nHat(IZ)

      end subroutine BR1_RiemannSolver

end module DGViscousDiscretization
