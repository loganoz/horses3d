module DGViscousDiscretization
   use SMConstants
   use MeshTypes
   use Physics
   use PhysicsStorage
   use MPI_Face_Class
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
      subroutine BaseClass_ComputeGradient( self , mesh , spA, time , externalStateProcedure , externalGradientsProcedure)
         use HexMeshClass
         use PhysicsStorage
         use Physics
         implicit none
         class(ViscousMethod_t), intent(in) :: self
         class(HexMesh)                   :: mesh
         class(NodalStorage)              :: spA(0:)
         real(kind=RP),        intent(in) :: time
         external                         :: externalStateProcedure
         external                         :: externalGradientsProcedure
!
!        ---------------------------
!        The base class does nothing
!        ---------------------------
!
      end subroutine BaseClass_ComputeGradient

      subroutine BaseClass_ComputeInnerFluxes( self , e , contravariantFlux )
         use ElementClass
         use PhysicsStorage
         implicit none
         class(ViscousMethod_t) ,  intent (in)   :: self
         type(Element)                           :: e
         real(kind=RP)           ,  intent (out) :: contravariantFlux(1:N_EQN, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM)
!
!        ---------------------------
!        The base class does nothing
!        ---------------------------
!
         contravariantFlux = 0.0_RP

      end subroutine BaseClass_ComputeInnerFluxes

      subroutine BaseClass_RiemannSolver ( self , QLeft , QRight , U_xLeft , U_yLeft , U_zLeft , U_xRight , U_yRight , U_zRight , &
                                           nHat , flux )
         use SMConstants
         use PhysicsStorage
         implicit none
         class(ViscousMethod_t)               :: self
         real(kind=RP), dimension(N_EQN)      :: QLeft
         real(kind=RP), dimension(N_EQN)      :: QRight
         real(kind=RP), dimension(N_GRAD_EQN) :: U_xLeft
         real(kind=RP), dimension(N_GRAD_EQN) :: U_yLeft
         real(kind=RP), dimension(N_GRAD_EQN) :: U_zLeft
         real(kind=RP), dimension(N_GRAD_EQN) :: U_xRight
         real(kind=RP), dimension(N_GRAD_EQN) :: U_yRight
         real(kind=RP), dimension(N_GRAD_EQN) :: U_zRight
         real(kind=RP), dimension(NDIM)       :: nHat
         real(kind=RP), dimension(N_EQN)      :: flux
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
      subroutine BR1_ComputeGradient( self , mesh , spA, time , externalStateProcedure , externalGradientsProcedure)
         use HexMeshClass
         use PhysicsStorage
         use Physics
         use MPI_Process_Info
         implicit none
         class(BassiRebay1_t), intent(in) :: self
         class(HexMesh)                   :: mesh
         class(NodalStorage)              :: spA(0:)
         real(kind=RP),        intent(in) :: time
         external                         :: externalStateProcedure
         external                         :: externalGradientsProcedure
         integer                          :: Nx, Ny, Nz
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                :: i, j, k
         integer                :: eID , fID , dimID , eqID, fIDs(6)
!
!        ************
!        Volume loops
!        ************
!
!$omp do schedule(runtime)
         do eID = 1 , size(mesh % elements)
            call BR1_GradientVolumeLoop( self , mesh % elements(eID) ) 
         end do
!$omp end do nowait
!
!        *******************************************
!        Compute Riemann solvers of non-shared faces
!        *******************************************
!
!$omp do schedule(runtime) 
         do fID = 1, SIZE(mesh % faces) 
            associate(f => mesh % faces(fID)) 
            select case (f % faceType) 
            case (HMESH_INTERIOR) 
               call BR1_ComputeElementInterfaceAverage(f) 
            
            case (HMESH_BOUNDARY) 
               call BR1_ComputeBoundaryFlux(f, time, externalStateProcedure) 
 
            end select 
            end associate 
         end do            
!$omp end do 
!
!$omp do schedule(runtime) 
         do eID = 1, size(mesh % elements) 
            associate(e => mesh % elements(eID))
            if ( e % hasSharedFaces ) cycle
!
!           Add the surface integrals
!           -------------------------
            call BR1_GradientFaceLoop( self , e, mesh)
!
!           Perform the scaling
!           -------------------               
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1)
               e % storage % U_x(:,i,j,k) = &
                           e % storage % U_x(:,i,j,k) / e % geom % jacobian(i,j,k)
               e % storage % U_y(:,i,j,k) = &
                           e % storage % U_y(:,i,j,k) / e % geom % jacobian(i,j,k)
               e % storage % U_z(:,i,j,k) = &
                           e % storage % U_z(:,i,j,k) / e % geom % jacobian(i,j,k)
            end do         ; end do          ; end do
!
!           Prolong gradients
!           -----------------
            fIDs = e % faceIDs
            call e % ProlongGradientsToFaces(mesh % faces(fIDs(1)),&
                                             mesh % faces(fIDs(2)),&
                                             mesh % faces(fIDs(3)),&
                                             mesh % faces(fIDs(4)),&
                                             mesh % faces(fIDs(5)),&
                                             mesh % faces(fIDs(6)) )

            end associate
         end do
!$omp end do

!$omp single
         if ( MPI_Process % doMPIAction ) then 
            call WaitUntilSolutionIsReady(mpi_faces)
         end if
!$omp end single


!$omp do schedule(runtime) 
         do fID = 1, SIZE(mesh % faces) 
            associate(f => mesh % faces(fID)) 
            select case (f % faceType) 
            case (HMESH_MPI) 
               call BR1_ComputeMPIFaceAverage(f) 
 
            end select 
            end associate 
         end do            
!$omp end do 
!
!$omp do schedule(runtime) 
         do eID = 1, size(mesh % elements) 
            associate(e => mesh % elements(eID))
            if ( .not. e % hasSharedFaces ) cycle
!
!           Add the surface integrals
!           -------------------------
            call BR1_GradientFaceLoop( self , e, mesh)
!
!           Perform the scaling
!           -------------------               
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)  ; do i = 0, e % Nxyz(1)
               e % storage % U_x(:,i,j,k) = &
                           e % storage % U_x(:,i,j,k) / e % geom % jacobian(i,j,k)
               e % storage % U_y(:,i,j,k) = &
                           e % storage % U_y(:,i,j,k) / e % geom % jacobian(i,j,k)
               e % storage % U_z(:,i,j,k) = &
                           e % storage % U_z(:,i,j,k) / e % geom % jacobian(i,j,k)
            end do         ; end do          ; end do
!
!           Prolong gradients
!           -----------------
            fIDs = e % faceIDs
            call e % ProlongGradientsToFaces(mesh % faces(fIDs(1)),&
                                             mesh % faces(fIDs(2)),&
                                             mesh % faces(fIDs(3)),&
                                             mesh % faces(fIDs(4)),&
                                             mesh % faces(fIDs(5)),&
                                             mesh % faces(fIDs(6)) )

            end associate
         end do
!$omp end do

      end subroutine BR1_ComputeGradient
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine BR1_GradientVolumeLoop( self , e )
         use ElementClass
         use PhysicsStorage
         use Physics
         use DGWeakIntegrals
         implicit none
         class(BassiRebay1_t),   intent(in)     :: self
         class(Element)                         :: e
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)          :: U(0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3),1:N_GRAD_EQN)
!
!        Compute gradient variables
!        --------------------------
         call GradientValuesForQ(e%Nxyz(1), e%Nxyz(2), e%Nxyz(3), Q = e % storage % Q, U = U )
!
!        Perform the weak integral
!        -------------------------
         call VectorWeakIntegrals % StdVolumeGreen (e , U , e % storage % U_x , e % storage % U_y , e % storage % U_z )

         e % storage % U_x = -e % storage % U_x
         e % storage % U_y = -e % storage % U_y
         e % storage % U_z = -e % storage % U_z

      end subroutine BR1_GradientVolumeLoop
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine BR1_GradientFaceLoop( self, e, mesh )
         use ElementClass
         use HexMeshClass
         use PhysicsStorage
         use Physics
         use DGWeakIntegrals
         implicit none
         class(BassiRebay1_t),   intent(in)  :: self
         class(Element)                      :: e
         class(HexMesh)                      :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)        :: faceInt_x(N_GRAD_EQN, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3) )
         real(kind=RP)        :: faceInt_y(N_GRAD_EQN, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3) )
         real(kind=RP)        :: faceInt_z(N_GRAD_EQN, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3) )

         call VectorWeakIntegrals % StdFace(e, &
               mesh % faces(e % faceIDs(EFRONT))  % storage(e % faceSide(EFRONT))  % unStar, &
               mesh % faces(e % faceIDs(EBACK))   % storage(e % faceSide(EBACK))   % unStar, &
               mesh % faces(e % faceIDs(EBOTTOM)) % storage(e % faceSide(EBOTTOM)) % unStar, &
               mesh % faces(e % faceIDs(ERIGHT))  % storage(e % faceSide(ERIGHT))  % unStar, &
               mesh % faces(e % faceIDs(ETOP))    % storage(e % faceSide(ETOP))    % unStar, &
               mesh % faces(e % faceIDs(ELEFT))   % storage(e % faceSide(ELEFT))   % unStar, &
               faceInt_x, faceInt_y, faceInt_z )

         e % storage % U_x = e % storage % U_x + faceInt_x
         e % storage % U_y = e % storage % U_y + faceInt_y
         e % storage % U_z = e % storage % U_z + faceInt_z

      end subroutine BR1_GradientFaceLoop
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine BR1_ComputeElementInterfaceAverage(f)
         use Physics  
         use ElementClass
         use FaceClass
         implicit none  
!
!        ---------
!        Arguments
!        ---------
!
         type(Face)    :: f
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: UL(N_GRAD_EQN), UR(N_GRAD_EQN)
         real(kind=RP) :: Uhat(N_GRAD_EQN)
         real(kind=RP) :: Hflux(N_GRAD_EQN,NDIM,0:f % Nf(1), 0:f % Nf(2))

         integer       :: i,j
         
         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            call GradientValuesForQ(Q = f % storage(1) % Q(:,i,j), U = UL)
            call GradientValuesForQ(Q = f % storage(2) % Q(:,i,j), U = UR)
   
            Uhat = 0.5_RP * (UL + UR) * f % geom % scal(i,j)
            Hflux(:,IX,i,j) = Uhat * f % geom % normal(IX,i,j)
            Hflux(:,IY,i,j) = Uhat * f % geom % normal(IY,i,j)
            Hflux(:,IZ,i,j) = Uhat * f % geom % normal(IZ,i,j)
         end do               ; end do

         call f % ProjectGradientFluxToElements(HFlux,(/1,2/))
         
      end subroutine BR1_ComputeElementInterfaceAverage   

      subroutine BR1_ComputeMPIFaceAverage(f)
         use Physics  
         use ElementClass
         use FaceClass
         implicit none  
!
!        ---------
!        Arguments
!        ---------
!
         type(Face)    :: f
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: UL(N_GRAD_EQN), UR(N_GRAD_EQN)
         real(kind=RP) :: Uhat(N_GRAD_EQN)
         real(kind=RP) :: Hflux(N_GRAD_EQN,NDIM,0:f % Nf(1), 0:f % Nf(2))
         integer       :: i,j, thisSide
         
         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            call GradientValuesForQ(Q = f % storage(1) % Q(:,i,j), U = UL)
            call GradientValuesForQ(Q = f % storage(2) % Q(:,i,j), U = UR)
   
            Uhat = 0.5_RP * (UL + UR) * f % geom % scal(i,j)
            Hflux(:,IX,i,j) = Uhat * f % geom % normal(IX,i,j)
            Hflux(:,IY,i,j) = Uhat * f % geom % normal(IY,i,j)
            Hflux(:,IZ,i,j) = Uhat * f % geom % normal(IZ,i,j)
         end do               ; end do

         thisSide = maxloc(f % elementIDs, dim = 1)
         call f % ProjectGradientFluxToElements(HFlux,(/thisSide, HMESH_NONE/))
         
      end subroutine BR1_ComputeMPIFaceAverage   

      subroutine BR1_ComputeBoundaryFlux(f, time, externalState)
         use Physics
         use FaceClass
         implicit none
         type(Face)    :: f
         real(kind=RP) :: time
         external      :: externalState
         integer       :: i, j
         real(kind=RP) :: Uhat(N_GRAD_EQN), UL(N_GRAD_EQN), UR(N_GRAD_EQN)
         real(kind=RP) :: bvExt(N_EQN)

         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)

            bvExt =  f % storage(1) % Q(:,i,j)
   
            call externalState( f % geom % x(:,i,j), &
                                time               , &
                                f % geom % normal(:,i,j)      , &
                                bvExt              , &
                                f % boundaryType )  
!   
!           -------------------
!           u, v, w, T averages
!           -------------------
!   
            call GradientValuesForQ( f % storage(1) % Q(:,i,j), UL )
            call GradientValuesForQ( bvExt, UR )
   
            Uhat = 0.5_RP * (UL + UR) * f % geom % scal(i,j)
            
            f % storage(1) % unStar(:,1,i,j) = Uhat * f % geom % normal(1,i,j)
            f % storage(1) % unStar(:,2,i,j) = Uhat * f % geom % normal(2,i,j)
            f % storage(1) % unStar(:,3,i,j) = Uhat * f % geom % normal(3,i,j)

         end do ; end do   
         
      end subroutine BR1_ComputeBoundaryFlux
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine BR1_ComputeInnerFluxes( self , e , contravariantFlux )
         use ElementClass
         use PhysicsStorage
         use Physics
         implicit none
         class(BassiRebay1_t) ,     intent (in) :: self
         type(Element)                          :: e
         real(kind=RP)           , intent (out) :: contravariantFlux(1:NCONS, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3), 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)       :: cartesianFlux(1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM)
         integer             :: i, j, k

         cartesianFlux = ViscousFlux( e%Nxyz(1) , e%Nxyz(2) , e%Nxyz(3)  , e % storage % Q , e % storage % U_x , e % storage % U_y , e % storage % U_z )

         do k = 0, e%Nxyz(3)   ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            contravariantFlux(:,i,j,k,IX) =     cartesianFlux(:,i,j,k,IX) * e % geom % jGradXi(IX,i,j,k)  &
                                             +  cartesianFlux(:,i,j,k,IY) * e % geom % jGradXi(IY,i,j,k)  &
                                             +  cartesianFlux(:,i,j,k,IZ) * e % geom % jGradXi(IZ,i,j,k)


            contravariantFlux(:,i,j,k,IY) =     cartesianFlux(:,i,j,k,IX) * e % geom % jGradEta(IX,i,j,k)  &
                                             +  cartesianFlux(:,i,j,k,IY) * e % geom % jGradEta(IY,i,j,k)  &
                                             +  cartesianFlux(:,i,j,k,IZ) * e % geom % jGradEta(IZ,i,j,k)


            contravariantFlux(:,i,j,k,IZ) =     cartesianFlux(:,i,j,k,IX) * e % geom % jGradZeta(IX,i,j,k)  &
                                             +  cartesianFlux(:,i,j,k,IY) * e % geom % jGradZeta(IY,i,j,k)  &
                                             +  cartesianFlux(:,i,j,k,IZ) * e % geom % jGradZeta(IZ,i,j,k)

         end do               ; end do            ; end do

      end subroutine BR1_ComputeInnerFluxes

      subroutine BR1_RiemannSolver ( self , QLeft , QRight , U_xLeft , U_yLeft , U_zLeft , U_xRight , U_yRight , U_zRight , &
                                            nHat , flux )
         use SMConstants
         use PhysicsStorage
         use Physics
         implicit none
         class(BassiRebay1_t)                 :: self
         real(kind=RP), dimension(N_EQN)      :: QLeft
         real(kind=RP), dimension(N_EQN)      :: QRight
         real(kind=RP), dimension(N_GRAD_EQN) :: U_xLeft
         real(kind=RP), dimension(N_GRAD_EQN) :: U_yLeft
         real(kind=RP), dimension(N_GRAD_EQN) :: U_zLeft
         real(kind=RP), dimension(N_GRAD_EQN) :: U_xRight
         real(kind=RP), dimension(N_GRAD_EQN) :: U_yRight
         real(kind=RP), dimension(N_GRAD_EQN) :: U_zRight
         real(kind=RP), dimension(NDIM)       :: nHat
         real(kind=RP), dimension(N_EQN)      :: flux
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
