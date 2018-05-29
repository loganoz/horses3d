#include "Includes.h"
module EllipticBR1
   use SMConstants
   use Headers
   use MeshTypes
   use Physics
   use VariableConversion
   use PhysicsStorage
   use MPI_Process_Info
   use MPI_Face_Class
   use EllipticDiscretizationClass
   use DGSEMClass, only: BCState_FCN
   use FluidData
   implicit none
!
!
   private
   public   BassiRebay1_t

   type, extends(EllipticDiscretization_t)   :: BassiRebay1_t
      contains
         procedure      :: ComputeGradient           => BR1_ComputeGradient
         procedure      :: ComputeInnerFluxes        => BR1_ComputeInnerFluxes
         procedure      :: RiemannSolver             => BR1_RiemannSolver
#if defined(NAVIERSTOKES)
         procedure      :: ComputeInnerFluxesWithSGS => BR1_ComputeInnerFluxesWithSGS
         procedure      :: RiemannSolverWithSGS      => BR1_RiemannSolverWithSGS
#endif
         procedure      :: Describe => BR1_Describe
   end type BassiRebay1_t
!
!  ========
   contains
!  ========
!
      subroutine BR1_Describe(self)
         implicit none
         class(BassiRebay1_t),   intent(in)  :: self  
!
!        Display the configuration
!        -------------------------
         if (MPI_Process % isRoot) write(STD_OUT,'(/)')
         call Subsection_Header("Viscous discretization")

         if (.not. MPI_Process % isRoot ) return

         write(STD_OUT,'(30X,A,A30,A)') "->","Numerical scheme: ","BR1"

      end subroutine BR1_Describe

      subroutine BR1_ComputeGradient(self, nEqn, nGradEqn, mesh, time, externalStateProcedure, GetGradients0D, GetGradients3D )
         use HexMeshClass
         use PhysicsStorage
         use Physics
         implicit none
         class(BassiRebay1_t), intent(in) :: self
         integer,              intent(in) :: nEqn, nGradEqn
         class(HexMesh)                   :: mesh
         real(kind=RP),        intent(in) :: time
         procedure(BCState_FCN)           :: externalStateProcedure
         procedure(GetGradientValues0D_f) :: GetGradients0D
         procedure(GetGradientValues3D_f) :: GetGradients3D
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
            call BR1_GradientVolumeLoop( self , nEqn, nGradEqn, mesh % elements(eID), GetGradients3D ) 
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
               call BR1_ComputeElementInterfaceAverage(f, nEqn, nGradEqn, GetGradients0D) 
            
            case (HMESH_BOUNDARY) 
               call BR1_ComputeBoundaryFlux(f, nEqn, nGradEqn, time, GetGradients0D, externalStateProcedure) 
 
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
            call BR1_GradientFaceLoop( self , nGradEqn, e, mesh)
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
            call e % ProlongGradientsToFaces(nGradEqn, mesh % faces(fIDs(1)),&
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
            call mesh % GatherMPIFacesSolution(nEqn)
         end if
!$omp end single

!$omp do schedule(runtime) 
         do fID = 1, SIZE(mesh % faces) 
            associate(f => mesh % faces(fID)) 
            select case (f % faceType) 
            case (HMESH_MPI) 
               call BR1_ComputeMPIFaceAverage(f, nEqn, nGradEqn, GetGradients0D) 
 
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
            call BR1_GradientFaceLoop(self, nGradEqn, e, mesh)
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
            call e % ProlongGradientsToFaces(nGradEqn, mesh % faces(fIDs(1)),&
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
      subroutine BR1_GradientVolumeLoop( self, nEqn, nGradEqn, e, GetGradients )
         use ElementClass
         use PhysicsStorage
         use Physics
         use DGWeakIntegrals
         implicit none
         class(BassiRebay1_t),   intent(in)     :: self
         integer,                intent(in)     :: nEqn, nGradEqn
         class(Element)                         :: e
         procedure(GetGradientValues3D_f)       :: GetGradients
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)          :: U(0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3),1:nGradEqn)
!
!        Compute gradient variables
!        --------------------------
         call GetGradients(nEqn, nGradEqn, e%Nxyz(1), e%Nxyz(2), e%Nxyz(3), Q = e % storage % Q, U = U )
!
!        Perform the weak integral
!        -------------------------
         call VectorWeakIntegrals % StdVolumeGreen (e , nGradEqn, U , e % storage % U_x , e % storage % U_y , e % storage % U_z )

         e % storage % U_x = -e % storage % U_x
         e % storage % U_y = -e % storage % U_y
         e % storage % U_z = -e % storage % U_z

      end subroutine BR1_GradientVolumeLoop
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine BR1_GradientFaceLoop( self, nGradEqn, e, mesh )
         use ElementClass
         use HexMeshClass
         use PhysicsStorage
         use Physics
         use DGWeakIntegrals
         implicit none
         class(BassiRebay1_t),   intent(in)  :: self
         integer,                intent(in)  :: nGradEqn
         class(Element)                      :: e
         class(HexMesh)                      :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)        :: faceInt_x(nGradEqn, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3) )
         real(kind=RP)        :: faceInt_y(nGradEqn, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3) )
         real(kind=RP)        :: faceInt_z(nGradEqn, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3) )

         call VectorWeakIntegrals % StdFace(e, nGradEqn, &
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
      subroutine BR1_ComputeElementInterfaceAverage(f, nEqn, nGradEqn, GetGradients)
         use Physics  
         use ElementClass
         use FaceClass
         implicit none  
!
!        ---------
!        Arguments
!        ---------
!
         type(Face)                       :: f
         integer,    intent(in)           :: nEqn, nGradEqn
         procedure(GetGradientValues0D_f) :: GetGradients
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: UL(nGradEqn), UR(nGradEqn)
         real(kind=RP) :: Uhat(nGradEqn)
         real(kind=RP) :: Hflux(nGradEqn,NDIM,0:f % Nf(1), 0:f % Nf(2))

         integer       :: i,j
         
         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            call GetGradients(nEqn, nGradEqn, Q = f % storage(1) % Q(:,i,j), U = UL)
            call GetGradients(nEqn, nGradEqn, Q = f % storage(2) % Q(:,i,j), U = UR)
   
            Uhat = 0.5_RP * (UL + UR) * f % geom % jacobian(i,j)
            Hflux(:,IX,i,j) = Uhat * f % geom % normal(IX,i,j)
            Hflux(:,IY,i,j) = Uhat * f % geom % normal(IY,i,j)
            Hflux(:,IZ,i,j) = Uhat * f % geom % normal(IZ,i,j)
         end do               ; end do

         call f % ProjectGradientFluxToElements(nGradEqn, HFlux,(/1,2/),-1)
         
      end subroutine BR1_ComputeElementInterfaceAverage   

      subroutine BR1_ComputeMPIFaceAverage(f, nEqn, nGradEqn, GetGradients)
         use Physics  
         use ElementClass
         use FaceClass
         implicit none  
!
!        ---------
!        Arguments
!        ---------
!
         type(Face)                       :: f
         integer, intent(in)              :: nEqn, nGradEqn
         procedure(GetGradientValues0D_f) :: GetGradients
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: UL(nGradEqn), UR(nGradEqn)
         real(kind=RP) :: Uhat(nGradEqn)
         real(kind=RP) :: Hflux(nGradEqn,NDIM,0:f % Nf(1), 0:f % Nf(2))
         integer       :: i,j, thisSide
         
         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            call GetGradients(nEqn, nGradEqn, Q = f % storage(1) % Q(:,i,j), U = UL)
            call GetGradients(nEqn, nGradEqn, Q = f % storage(2) % Q(:,i,j), U = UR)
   
            Uhat = 0.5_RP * (UL + UR) * f % geom % jacobian(i,j)
            Hflux(:,IX,i,j) = Uhat * f % geom % normal(IX,i,j)
            Hflux(:,IY,i,j) = Uhat * f % geom % normal(IY,i,j)
            Hflux(:,IZ,i,j) = Uhat * f % geom % normal(IZ,i,j)
         end do               ; end do

         thisSide = maxloc(f % elementIDs, dim = 1)
         call f % ProjectGradientFluxToElements(nGradEqn, HFlux,(/thisSide, HMESH_NONE/),-1)
         
      end subroutine BR1_ComputeMPIFaceAverage   

      subroutine BR1_ComputeBoundaryFlux(f, nEqn, nGradEqn, time, GetGradients, externalState)
         use Physics
         use FaceClass
         implicit none
         type(Face)                       :: f
         integer, intent(in)              :: nEqn, nGradEqn
         real(kind=RP), intent(in)        :: time
         procedure(GetGradientValues0D_f) :: GetGradients
         procedure(BCState_FCN)           :: externalState
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: i, j
         real(kind=RP) :: Uhat(nGradEqn), UL(nGradEqn), UR(nGradEqn)
         real(kind=RP) :: bvExt(nEqn)

         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)

            bvExt =  f % storage(1) % Q(:,i,j)
   
            call externalState( nEqn, &
                                f % geom % x(:,i,j), &
                                time               , &
                                f % geom % normal(:,i,j)      , &
                                bvExt              , &
                                f % boundaryType, f % boundaryName )  
!   
!           -------------------
!           u, v, w, T averages
!           -------------------
!   
            call GetGradients(nEqn, nGradEqn, f % storage(1) % Q(:,i,j), UL )
            call GetGradients(nEqn, nGradEqn, bvExt, UR )
   
            Uhat = 0.5_RP * (UL + UR) * f % geom % jacobian(i,j)
            
            f % storage(1) % unStar(:,1,i,j) = Uhat * f % geom % normal(1,i,j)
            f % storage(1) % unStar(:,2,i,j) = Uhat * f % geom % normal(2,i,j)
            f % storage(1) % unStar(:,3,i,j) = Uhat * f % geom % normal(3,i,j)

         end do ; end do   
         
      end subroutine BR1_ComputeBoundaryFlux
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine BR1_ComputeInnerFluxes( self , nEqn, nGradEqn, e , contravariantFlux )
         use ElementClass
         use PhysicsStorage
         use Physics
         implicit none
         class(BassiRebay1_t) ,     intent (in) :: self
         integer,                   intent(in)  :: nEqn
         integer,                   intent(in)  :: nGradEqn
         type(Element)                          :: e
         real(kind=RP)           , intent (out) :: contravariantFlux(1:nEqn, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3), 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)       :: delta
         real(kind=RP)       :: cartesianFlux(1:nEqn, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM)
         real(kind=RP)       :: mu(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: kappa(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         integer             :: i, j, k

#if defined(CAHNHILLIARD)
         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
            call self % GetViscosity(e % storage % c(1,i,j,k), mu(i,j,k))      
         end do                ; end do                ; end do
#else
         mu = dimensionless % mu

#endif

#if defined(NAVIERSTOKES)
         kappa = 1.0_RP / ( thermodynamics % gammaMinus1 * &
                               POW2( dimensionless % Mach) * dimensionless % Pr ) * mu
#else
         kappa = 0.0_RP
#endif

         call self % EllipticFlux3D( nEqn, nGradEqn, e%Nxyz, e % storage % Q , e % storage % U_x , e % storage % U_y , e % storage % U_z, mu, kappa, cartesianFlux )

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
#if defined(NAVIERSTOKES)
      subroutine BR1_ComputeInnerFluxesWithSGS( self , e , contravariantFlux )
         use ElementClass
         use PhysicsStorage
         use Physics
         use LESModels
         implicit none
         class(BassiRebay1_t) ,     intent (in) :: self
         type(Element)                          :: e
         real(kind=RP)           , intent (out) :: contravariantFlux(1:NCONS, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3), 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: delta
         real(kind=RP) :: cartesianFlux(1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM)
         real(kind=RP) :: mu(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP) :: kappa(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP) :: tauSGS(1:NDIM,1:NDIM, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP) :: qSGS(1:NDIM, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         integer       :: i, j, k

         mu = dimensionless % mu
         kappa = dimensionless % kappa
!
!        Compute subgrid-scale modelling tensor   
!        --------------------------------------
         delta = (e % geom % Volume / product(e % Nxyz + 1)) ** (1.0_RP / 3.0_RP)
         call LESModel % ComputeSGSTensor(delta, e % Nxyz, e % geom % dWall, &
                                                           e % storage % U_x, &
                                                           e % storage % U_y, &
                                                           e % storage % U_z, &
                                                                tauSGS, qSGS    )

         call ViscousFlux(NCONS, NGRAD, e%Nxyz, e % storage % Q , e % storage % U_x , e % storage % U_y , e % storage % U_z, mu, kappa, tauSGS, qSGS, cartesianFlux )

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

      end subroutine BR1_ComputeInnerFluxesWithSGS
#endif
      subroutine BR1_RiemannSolver ( self , nEqn, nGradEqn, f, QLeft , QRight , U_xLeft , U_yLeft , U_zLeft , U_xRight , U_yRight , U_zRight , &
                                            mu, nHat , dWall, flux )
         use SMConstants
         use PhysicsStorage
         use Physics
         use FaceClass
         implicit none
         class(BassiRebay1_t)            :: self
         integer,       intent(in)       :: nEqn
         integer,       intent(in)       :: nGradEqn
         class(Face),   intent(in)       :: f
         real(kind=RP), intent(in)       :: QLeft(nEqn)
         real(kind=RP), intent(in)       :: QRight(nEqn)
         real(kind=RP), intent(in)       :: U_xLeft(nGradEqn)
         real(kind=RP), intent(in)       :: U_yLeft(nGradEqn)
         real(kind=RP), intent(in)       :: U_zLeft(nGradEqn)
         real(kind=RP), intent(in)       :: U_xRight(nGradEqn)
         real(kind=RP), intent(in)       :: U_yRight(nGradEqn)
         real(kind=RP), intent(in)       :: U_zRight(nGradEqn)
         real(kind=RP), intent(in)       :: mu
         real(kind=RP), intent(in)       :: nHat(NDIM)
         real(kind=RP), intent(in)       :: dWall
         real(kind=RP), intent(out)      :: flux(nEqn)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: Q(nEqn) , U_x(nGradEqn) , U_y(nGradEqn) , U_z(nGradEqn)
         real(kind=RP)     :: flux_vec(nEqn,NDIM)
         real(kind=RP)     :: kappa
!
!>       Old implementation: 1st average, then compute
!        ------------------
         Q   = 0.5_RP * ( QLeft + QRight)
         U_x = 0.5_RP * ( U_xLeft + U_xRight)
         U_y = 0.5_RP * ( U_yLeft + U_yRight)
         U_z = 0.5_RP * ( U_zLeft + U_zRight)

#if defined(NAVIERSTOKES)
         kappa = 1.0_RP / ( thermodynamics % gammaMinus1 * &
                            POW2( dimensionless % Mach) * dimensionless % Pr ) * mu
#else
         kappa = 0.0_RP
#endif


         call self % EllipticFlux0D(nEqn, nGradEqn, Q,U_x,U_y,U_z, mu, kappa, flux_vec)

         flux = flux_vec(:,IX) * nHat(IX) + flux_vec(:,IY) * nHat(IY) + flux_vec(:,IZ) * nHat(IZ)

      end subroutine BR1_RiemannSolver
#if defined(NAVIERSTOKES)
      subroutine BR1_RiemannSolverWithSGS ( self , f, QLeft , QRight , U_xLeft , U_yLeft , U_zLeft , U_xRight , U_yRight , U_zRight , &
                                            nHat , dWall, flux )
         use SMConstants
         use PhysicsStorage
         use Physics
         use FaceClass
         use LESModels
         implicit none
         class(BassiRebay1_t)            :: self
         class(Face),   intent(in)       :: f
         real(kind=RP), dimension(NCONS) :: QLeft
         real(kind=RP), dimension(NCONS) :: QRight
         real(kind=RP), dimension(NGRAD) :: U_xLeft
         real(kind=RP), dimension(NGRAD) :: U_yLeft
         real(kind=RP), dimension(NGRAD) :: U_zLeft
         real(kind=RP), dimension(NGRAD) :: U_xRight
         real(kind=RP), dimension(NGRAD) :: U_yRight
         real(kind=RP), dimension(NGRAD) :: U_zRight
         real(kind=RP), dimension(NDIM)  :: nHat
         real(kind=RP)                   :: dWall
         real(kind=RP), dimension(NCONS) :: flux
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: Q(NCONS) , U_x(NGRAD) , U_y(NGRAD) , U_z(NGRAD)
         real(kind=RP)     :: flux_vec(NCONS,NDIM)
         real(kind=RP)     :: mu, kappa, tauSGS(NDIM, NDIM), qSGS(NDIM), delta

!
!>       Old implementation: 1st average, then compute
!        ------------------
         Q   = 0.5_RP * ( QLeft + QRight)
         U_x = 0.5_RP * ( U_xLeft + U_xRight)
         U_y = 0.5_RP * ( U_yLeft + U_yRight)
         U_z = 0.5_RP * ( U_zLeft + U_zRight)
!
!        Compute subgrid-scale modelling tensor   
!        --------------------------------------
         delta = sqrt(f % geom % surface / product(f % Nf + 1))
         call LESModel % ComputeSGSTensor(delta, dWall, U_x, U_y, U_z, tauSGS, qSGS) 

         mu    = dimensionless % mu
         kappa = dimensionless % kappa

         call ViscousFlux(NCONS, NGRAD, Q,U_x,U_y,U_z, mu, kappa, tauSGS, qSGS, flux_vec)

         flux = flux_vec(:,IX) * nHat(IX) + flux_vec(:,IY) * nHat(IY) + flux_vec(:,IZ) * nHat(IZ)

      end subroutine BR1_RiemannSolverWithSGS
#endif
end module EllipticBR1
