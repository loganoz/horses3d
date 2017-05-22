module DGViscousDiscretization
   implicit none
!
!
   private
   public   ViscousMethod_t , BassiRebay1_t
!
!
!  *****************************
!  Generic viscous methods class
!  *****************************
!
   type ViscousMethod_t
      contains
         procedure      :: ComputeGradient      => BaseClass_ComputeGradient
         procedure      :: ComputeInnerFluxes    => BaseClass_ComputeInnerFluxes
   end type ViscousMethod_t

   type, extends(ViscousMethod_t)   :: BassiRebay1_t
      contains
         procedure      :: ComputeGradient      => BR1_ComputeGradient
         procedure      :: ComputeInnerFluxes    => BR1_ComputeInnerFluxes

   end type BassiRebay1_t

   type, extends(ViscousMethod_t)   :: InteriorPenalty_t
      contains
!         procedure      :: ComputeGradient     => IP_ComputeGradient
!         procedure      :: ComputeInnerFluxes   => IP_ComputeInnerFluxes
   end type InteriorPenalty_t
!
!

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
      subroutine BaseClass_ComputeGradient( self , mesh , spA , time , externalStateProcedure )
         use HexMeshClass
         use NodalStorageClass
         use PhysicsStorage
         use Physics
         implicit none
         class(ViscousMethod_t),    intent(in)     :: self
         class(HexMesh)                            :: mesh
         class(NodalStorage),       intent(in)     :: spA
         real(kind=RP),             intent(in)     :: time
         external                                  :: externalStateProcedure
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
      end subroutine BaseClass_ComputeInnerFluxes
!
!///////////////////////////////////////////////////////////////////////////////////
!
!           Bassi-Rebay 1 Procedures
!           ------------------------
!///////////////////////////////////////////////////////////////////////////////////
!
      subroutine BR1_ComputeGradient( self , mesh , spA , time , externalStateProcedure )
         use HexMeshClass
         use NodalStorageClass
         use PhysicsStorage
         use Physics
         implicit none
         class(BassiRebay1_t), intent(in) :: self
         class(HexMesh)                   :: mesh
         class(NodalStorage),  intent(in) :: spA
         real(kind=RP),        intent(in) :: time
         external                         :: externalStateProcedure
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                :: eID , fID , dimID , eqID
!
!        Perform volume loops
!        --------------------
         do eID = 1 , size(mesh % elements)
            call BR1_GradientVolumeLoop( self , mesh % elements(eID) , spA ) 
         end do
!
!        Perform face loops
!        ------------------
         call BR1_GradientFaceLoop( self , mesh , spA , time, externalStateProcedure )


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
         call VectorWeakIntegrals % StdVolumeGreen ( N_GRAD_EQN , e , spA , -U , e % U_x , e % U_y , e % U_z )

      end subroutine BR1_GradientVolumeLoop

      subroutine BR1_GradientFaceLoop( self , mesh , spA , time, externalStateProcedure )
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
!TODO: omp loop ($omp do)        
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
!TODO ($omp enddo)       
         
      end subroutine BR1_GradientFaceLoop

      SUBROUTINE BR1_ComputeElementInterfaceAverage( eL, fIDLeft, eR, fIDRight, N, rotation)
         USE Physics  
         use ElementClass
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
                                  
! TODO: Okay, so I have just removed the solution average
!               eL % QB(:,i,j,fIDLeft)    = 0.5_RP * ( eL % QB(:,i,j,fIDLeft) + eR % QB(:,ii,jj,fIDright) )
!               eR % QB(:,ii,jj,fIDright) = eL % QB(:,i,j,fIDLeft)
               
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
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////
!
!           AUXILIAR PROCEDURES
!           -------------------
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////
!

end module DGViscousDiscretization
