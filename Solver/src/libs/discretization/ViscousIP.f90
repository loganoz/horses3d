!
!//////////////////////////////////////////////////////
!
!   @File:    ViscousIP.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Dec 12 13:32:09 2017
!   @Last revision date: Wed Dec 13 20:00:30 2017
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: dd41d1867321d548cfb82fb2c4ae6d579d7ba5bb
!
!//////////////////////////////////////////////////////
!
module ViscousIP
   use SMConstants
   use MeshTypes
   use ElementClass
   use HexMeshClass
   use Physics
   use PhysicsStorage
   use MPI_Face_Class
   use ViscousMethodClass
   implicit none
!
!
   private
   public   InteriorPenalty_t, SIPG, IIPG, NIPG

   integer, parameter   :: SIPG = -1
   integer, parameter   :: IIPG = 0
   integer, parameter   :: NIPG = 1

   type, extends(ViscousMethod_t)   :: InteriorPenalty_t
      real(kind=RP)        :: sigma = 1.0_RP
      integer              :: IPmethod = SIPG
      contains
         procedure      :: ComputeGradient    => IP_ComputeGradient
         procedure      :: ComputeInnerFluxes => IP_ComputeInnerFluxes
         procedure      :: RiemannSolver      => IP_RiemannSolver
   end type InteriorPenalty_t
!
!  ========
   contains
!  ========
!
      subroutine IP_ComputeGradient( self , mesh , time , externalStateProcedure , externalGradientsProcedure)
         use HexMeshClass
         use PhysicsStorage
         use Physics
         use MPI_Process_Info
         implicit none
         class(InteriorPenalty_t), intent(in) :: self
         class(HexMesh)                   :: mesh
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
!        *********************************
!        Volume loops and prolong to faces
!        *********************************
!
!$omp do schedule(runtime)
         do eID = 1, size(mesh % elements)
            associate( e => mesh % elements(eID) )
            call e % ComputeLocalGradient
!
!           Prolong to faces
!           ----------------
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
!
!        **********************************************
!        Compute interface solution of non-shared faces
!        **********************************************
!
!$omp do schedule(runtime) 
         do fID = 1, SIZE(mesh % faces) 
            associate(f => mesh % faces(fID)) 
            select case (f % faceType) 
            case (HMESH_INTERIOR) 
               call IP_GradientInterfaceSolution(f) 
            
            case (HMESH_BOUNDARY) 
               call IP_GradientInterfaceSolutionBoundary(f, time, externalStateProcedure) 
 
            end select 
            end associate 
         end do            
!$omp end do 
!
!        **********************
!        Compute face integrals
!        **********************
!
!$omp do schedule(runtime) 
         do eID = 1, size(mesh % elements) 
            associate(e => mesh % elements(eID))
            if ( e % hasSharedFaces ) cycle
            call IP_ComputeGradientFaceIntegrals(self, e, mesh)
            end associate
         end do
!$omp end do
!
!        ******************
!        Wait for MPI faces
!        ******************
!
!$omp single
         if ( MPI_Process % doMPIAction ) then 
            call mesh % GatherMPIFacesSolution
         end if
!$omp end single
!
!        *******************************
!        Compute MPI interface solutions
!        *******************************
!
!$omp do schedule(runtime) 
         do fID = 1, SIZE(mesh % faces) 
            associate(f => mesh % faces(fID)) 
            select case (f % faceType) 
            case (HMESH_MPI) 
               call IP_GradientInterfaceSolutionMPI(f) 
 
            end select 
            end associate 
         end do            
!$omp end do 
!
!        **************************************************
!        Compute face integrals for elements with MPI faces
!        **************************************************
!
!$omp do schedule(runtime) 
         do eID = 1, size(mesh % elements) 
            associate(e => mesh % elements(eID))
            if ( .not. e % hasSharedFaces ) cycle
            call IP_ComputeGradientFaceIntegrals(self, e, mesh)
            end associate
         end do
!$omp end do

      end subroutine IP_ComputeGradient
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine IP_ComputeGradientFaceIntegrals( self, e, mesh)
         use ElementClass
         use HexMeshClass
         use PhysicsStorage
         use Physics
         use DGWeakIntegrals
         implicit none
         class(InteriorPenalty_t),   intent(in) :: self
         class(Element)                         :: e
         class(HexMesh)                         :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
         integer              :: i, j, k
         real(kind=RP)        :: invjac
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
!
!        Add the integrals weighted with the Jacobian
!        --------------------------------------------
         do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2)    ; do i = 0, e % Nxyz(1)
            invjac = (- self % IPmethod) / e % geom % jacobian(i,j,k)
            e % storage % U_x(:,i,j,k) = e % storage % U_x(:,i,j,k) + faceInt_x(:,i,j,k) * invjac
            e % storage % U_y(:,i,j,k) = e % storage % U_y(:,i,j,k) + faceInt_y(:,i,j,k) * invjac
            e % storage % U_z(:,i,j,k) = e % storage % U_z(:,i,j,k) + faceInt_z(:,i,j,k) * invjac
         end do                  ; end do                   ; end do
!
      end subroutine IP_ComputeGradientFaceIntegrals
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine IP_GradientInterfaceSolution(f)
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
   
            Uhat = 0.5_RP * (UR - UL) * f % geom % jacobian(i,j)
            Hflux(:,IX,i,j) = Uhat * f % geom % normal(IX,i,j)
            Hflux(:,IY,i,j) = Uhat * f % geom % normal(IY,i,j)
            Hflux(:,IZ,i,j) = Uhat * f % geom % normal(IZ,i,j)
         end do               ; end do

         call f % ProjectGradientFluxToElements(HFlux,(/1,2/))
         
      end subroutine IP_GradientInterfaceSolution   

      subroutine IP_GradientInterfaceSolutionMPI(f)
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
   
            Uhat = 0.5_RP * (UR - UL) * f % geom % jacobian(i,j)
            Hflux(:,IX,i,j) = Uhat * f % geom % normal(IX,i,j)
            Hflux(:,IY,i,j) = Uhat * f % geom % normal(IY,i,j)
            Hflux(:,IZ,i,j) = Uhat * f % geom % normal(IZ,i,j)
         end do               ; end do

         thisSide = maxloc(f % elementIDs, dim = 1)
         call f % ProjectGradientFluxToElements(HFlux,(/thisSide, HMESH_NONE/))
         
      end subroutine IP_GradientInterfaceSolutionMPI   

      subroutine IP_GradientInterfaceSolutionBoundary(f, time, externalState)
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
   
            Uhat = 0.5_RP * (UR - UL) * f % geom % jacobian(i,j)
            
            f % storage(1) % unStar(:,1,i,j) = Uhat * f % geom % normal(1,i,j)
            f % storage(1) % unStar(:,2,i,j) = Uhat * f % geom % normal(2,i,j)
            f % storage(1) % unStar(:,3,i,j) = Uhat * f % geom % normal(3,i,j)

         end do ; end do   
         
      end subroutine IP_GradientInterfaceSolutionBoundary
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine IP_ComputeInnerFluxes( self , e , contravariantFlux )
         use ElementClass
         use PhysicsStorage
         use Physics
         implicit none
         class(InteriorPenalty_t) ,     intent (in) :: self
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

      end subroutine IP_ComputeInnerFluxes

      subroutine IP_RiemannSolver(self , f, QLeft , QRight , U_xLeft , U_yLeft , U_zLeft , U_xRight , U_yRight , U_zRight , &
                                            nHat , flux )
         use SMConstants
         use PhysicsStorage
         use Physics
         implicit none
         class(InteriorPenalty_t)                 :: self
         class(Face), intent(in)                  :: f
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
         real(kind=RP)     :: sigma

!
!>       Old implementation: 1st average, then compute
!        ------------------
         Q   = 0.5_RP * ( QLeft + QRight)
         U_x = 0.5_RP * ( U_xLeft + U_xRight)
         U_y = 0.5_RP * ( U_yLeft + U_yRight)
         U_z = 0.5_RP * ( U_zLeft + U_zRight)

         flux_vec = ViscousFlux(Q,U_x,U_y,U_z)
!
!        Shahbazi estimate
!        -----------------
         self % sigma = 1.0_RP
         sigma = 0.5_RP * self % sigma * (maxval(f % Nf)+1)*(maxval(f % Nf)+2) / f % geom % h * dimensionless % mu

         flux = flux_vec(:,IX) * nHat(IX) + flux_vec(:,IY) * nHat(IY) + flux_vec(:,IZ) * nHat(IZ) + sigma * (QRight - QLeft)

      end subroutine IP_RiemannSolver

end module ViscousIP
