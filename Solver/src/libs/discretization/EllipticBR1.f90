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
   use FluidData
   use BoundaryConditions, only: BCs
   implicit none
!
!
   private
   public   BassiRebay1_t

   type, extends(EllipticDiscretization_t)   :: BassiRebay1_t
      contains
         procedure      :: ComputeGradient           => BR1_ComputeGradient
         procedure      :: LiftGradients             => BR1_LiftGradients
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

         select case (self % eqName)
         case (ELLIPTIC_NS,ELLIPTIC_iNS,ELLIPTIC_MU)
            call Subsection_Header("Viscous discretization")
      
         case (ELLIPTIC_CH)
            call Subsection_Header("Cahn--Hilliard discretization")

         end select
   

         if (.not. MPI_Process % isRoot ) return


         write(STD_OUT,'(30X,A,A30,A)') "->","Numerical scheme: ","BR1"
         write(STD_OUT,'(30X,A,A30,F6.3)') "->","Penalty parameter: ",self % sigma

      end subroutine BR1_Describe

      subroutine BR1_ComputeGradient(self, nEqn, nGradEqn, mesh, time, GetGradients)
         use HexMeshClass
         use PhysicsStorage
         use Physics
         implicit none
         class(BassiRebay1_t), intent(in) :: self
         integer,              intent(in) :: nEqn, nGradEqn
         class(HexMesh)                   :: mesh
         real(kind=RP),        intent(in) :: time
         procedure(GetGradientValues_f)   :: GetGradients
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
            call BR1_GradientVolumeLoop( self , nEqn, nGradEqn, mesh % elements(eID), GetGradients) 
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
               call BR1_ComputeElementInterfaceAverage(self, f, nEqn, nGradEqn, GetGradients) 
            
            case (HMESH_BOUNDARY) 
               call BR1_ComputeBoundaryFlux(self, f, nEqn, nGradEqn, time, GetGradients) 
 
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

#ifdef _HAS_MPI_
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
               call BR1_ComputeMPIFaceAverage(self, f, nEqn, nGradEqn, GetGradients) 
 
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
#endif

      end subroutine BR1_ComputeGradient

      subroutine BR1_LiftGradients(self, nEqn, nGradEqn, mesh, time, GetGradients)
!
!        **********************************************************************
!        This subroutine performs the gradient lifting in the same way the IP
!        does. Thus, I recycle the subroutines from the IP method to avoid
!        duplicates.
!        **********************************************************************
!
         use HexMeshClass
         use PhysicsStorage
         use Physics
         use EllipticIP, only: IP_GradientInterfaceSolution, IP_GradientInterfaceSolutionBoundary, &
                               IP_GradientInterfaceSolutionMPI
         implicit none
         class(BassiRebay1_t), intent(in) :: self
         integer,              intent(in) :: nEqn, nGradEqn
         class(HexMesh)                   :: mesh
         real(kind=RP),        intent(in) :: time
         procedure(GetGradientValues_f)   :: GetGradients
         integer                          :: Nx, Ny, Nz
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                :: i, j, k
         integer                :: eID , fID , dimID , eqID, fIDs(6)
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
               call IP_GradientInterfaceSolution(f, nEqn, nGradEqn, GetGradients) 
            
            case (HMESH_BOUNDARY) 
               call IP_GradientInterfaceSolutionBoundary(f, nEqn, nGradEqn, time, GetGradients) 
 
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
            call ComputeLiftGradientFaceIntegrals(nGradEqn, e, mesh)
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

      end subroutine BR1_LiftGradients

      subroutine ComputeLiftGradientFaceIntegrals(nGradEqn, e, mesh)
         use ElementClass
         use HexMeshClass
         use PhysicsStorage
         use Physics
         use DGIntegrals
         implicit none
         integer,                    intent(in) :: nGradEqn
         class(Element)                         :: e
         class(HexMesh)                         :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
         integer              :: i, j, k
         real(kind=RP)        :: invjac
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
!
!        Add the integrals weighted with the Jacobian
!        --------------------------------------------
         do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2)    ; do i = 0, e % Nxyz(1)
            e % storage % U_x(:,i,j,k) = e % storage % U_x(:,i,j,k) - faceInt_x(:,i,j,k) * e % geom % InvJacobian(i,j,k)
            e % storage % U_y(:,i,j,k) = e % storage % U_y(:,i,j,k) - faceInt_y(:,i,j,k) * e % geom % InvJacobian(i,j,k)
            e % storage % U_z(:,i,j,k) = e % storage % U_z(:,i,j,k) - faceInt_z(:,i,j,k) * e % geom % InvJacobian(i,j,k)
         end do                  ; end do                   ; end do
!
      end subroutine ComputeLiftGradientFaceIntegrals
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine BR1_GradientVolumeLoop( self, nEqn, nGradEqn, e, GetGradients )
         use ElementClass
         use PhysicsStorage
         use Physics
         use DGIntegrals
         implicit none
         class(BassiRebay1_t),   intent(in)     :: self
         integer,                intent(in)     :: nEqn, nGradEqn
         class(Element)                         :: e
         procedure(GetGradientValues_f)         :: GetGradients
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                :: i,j,k
         real(kind=RP)          :: U(1:nGradEqn,0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3))
!
!        Compute gradient variables
!        --------------------------
#ifdef MULTIPHASE
         do k = 0, e%Nxyz(3) ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            call GetGradients(nEqn, nGradEqn, e % storage % Q(:,i,j,k), U(:,i,j,k), e % storage % rho(i,j,k) )
         end do         ; end do         ; end do
#else
         do k = 0, e%Nxyz(3) ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            call GetGradients(nEqn, nGradEqn, e % storage % Q(:,i,j,k), U(:,i,j,k))
         end do         ; end do         ; end do
#endif

#ifdef MULTIPHASE
         select case (self % eqName)
         case(ELLIPTIC_MU)
!
!           The multiphase solver needs the Chemical potential as first entropy variable
!           ----------------------------------------------------------------------------
            U(IGMU,:,:,:) = e % storage % mu(1,:,:,:)
         end select
#endif
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
         use DGIntegrals
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
      subroutine BR1_ComputeElementInterfaceAverage(self, f, nEqn, nGradEqn, GetGradients)
         use Physics  
         use ElementClass
         use FaceClass
         implicit none  
!
!        ---------
!        Arguments
!        ---------
!
         class(BassiRebay1_t),   intent(in)  :: self
         type(Face)                       :: f
         integer,    intent(in)           :: nEqn, nGradEqn
         procedure(GetGradientValues_f)   :: GetGradients
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: UL(nGradEqn), UR(nGradEqn)
         real(kind=RP) :: Uhat(nGradEqn)
         real(kind=RP) :: Hflux(nGradEqn,NDIM,0:f % Nf(1), 0:f % Nf(2))

         integer       :: i,j
         integer       :: Sidearray(2)
         
         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
#ifdef MULTIPHASE
            call GetGradients(nEqn, nGradEqn, Q = f % storage(1) % Q(:,i,j), U = UL, rho_ = f % storage(1) % rho(i,j))
            call GetGradients(nEqn, nGradEqn, Q = f % storage(2) % Q(:,i,j), U = UR, rho_ = f % storage(2) % rho(i,j))
#else
            call GetGradients(nEqn, nGradEqn, Q = f % storage(1) % Q(:,i,j), U = UL)
            call GetGradients(nEqn, nGradEqn, Q = f % storage(2) % Q(:,i,j), U = UR)
#endif

#ifdef MULTIPHASE
            select case (self % eqName)
            case (ELLIPTIC_MU)
!
!              The multiphase solver needs the Chemical potential as first entropy variable
!              ----------------------------------------------------------------------------
               UL(IGMU) = f % storage(1) % mu(1,i,j)
               UR(IGMU) = f % storage(2) % mu(1,i,j)
            end select
#endif

   
            Uhat = 0.5_RP * (UL + UR) * f % geom % jacobian(i,j)
            Hflux(:,IX,i,j) = Uhat * f % geom % normal(IX,i,j)
            Hflux(:,IY,i,j) = Uhat * f % geom % normal(IY,i,j)
            Hflux(:,IZ,i,j) = Uhat * f % geom % normal(IZ,i,j)
         end do               ; end do
         
         Sidearray = (/1,2/)
         call f % ProjectGradientFluxToElements(nGradEqn, HFlux,Sidearray,-1)
         
      end subroutine BR1_ComputeElementInterfaceAverage   

      subroutine BR1_ComputeMPIFaceAverage(self, f, nEqn, nGradEqn, GetGradients)
         use Physics  
         use ElementClass
         use FaceClass
         implicit none  
!
!        ---------
!        Arguments
!        ---------
!
         class(BassiRebay1_t),   intent(in)  :: self
         type(Face)                       :: f
         integer, intent(in)              :: nEqn, nGradEqn
         procedure(GetGradientValues_f)   :: GetGradients
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: UL(nGradEqn), UR(nGradEqn)
         real(kind=RP) :: Uhat(nGradEqn)
         real(kind=RP) :: Hflux(nGradEqn,NDIM,0:f % Nf(1), 0:f % Nf(2))
         integer       :: i,j
         integer       :: Sidearray(2)
         
         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
#ifdef MULTIPHASE
            call GetGradients(nEqn, nGradEqn, Q = f % storage(1) % Q(:,i,j), U = UL, rho_ = f % storage(1) % rho(i,j))
            call GetGradients(nEqn, nGradEqn, Q = f % storage(2) % Q(:,i,j), U = UR, rho_ = f % storage(2) % rho(i,j))
#else
            call GetGradients(nEqn, nGradEqn, Q = f % storage(1) % Q(:,i,j), U = UL)
            call GetGradients(nEqn, nGradEqn, Q = f % storage(2) % Q(:,i,j), U = UR)
#endif

#ifdef MULTIPHASE
            select case (self % eqName)
            case (ELLIPTIC_MU)
!
!              The multiphase solver needs the Chemical potential as first entropy variable
!              ----------------------------------------------------------------------------
               UL(IGMU) = f % storage(1) % mu(1,i,j)
               UR(IGMU) = f % storage(2) % mu(1,i,j)
            end select
#endif

   
            Uhat = 0.5_RP * (UL + UR) * f % geom % jacobian(i,j)
            Hflux(:,IX,i,j) = Uhat * f % geom % normal(IX,i,j)
            Hflux(:,IY,i,j) = Uhat * f % geom % normal(IY,i,j)
            Hflux(:,IZ,i,j) = Uhat * f % geom % normal(IZ,i,j)
         end do               ; end do

         Sidearray = (/maxloc(f % elementIDs, dim = 1), HMESH_NONE/)
         call f % ProjectGradientFluxToElements(nGradEqn, HFlux,Sidearray,-1)
         
      end subroutine BR1_ComputeMPIFaceAverage   

      subroutine BR1_ComputeBoundaryFlux(self, f, nEqn, nGradEqn, time, GetGradients)
         use Physics
         use FaceClass
         implicit none
         class(BassiRebay1_t),   intent(in)  :: self
         type(Face)                       :: f
         integer, intent(in)              :: nEqn, nGradEqn
         real(kind=RP), intent(in)        :: time
         procedure(GetGradientValues_f)   :: GetGradients
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: i, j, bcID
         real(kind=RP) :: Uhat(nGradEqn)
         real(kind=RP) :: bvExt(nEqn)

         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)

#ifdef MULTIPHASE
            call GetGradients(nEqn, nGradEqn, f % storage(1) % Q(:,i,j), uHat, f % storage(1) % rho(i,j))
#else
            call GetGradients(nEqn, nGradEqn, f % storage(1) % Q(:,i,j), uHat)
#endif

#ifdef MULTIPHASE
            select case (self % eqName)
            case (ELLIPTIC_MU)
!
!              The multiphase solver needs the Chemical potential as first entropy variable
!              ----------------------------------------------------------------------------
               uHat(IGMU) = f % storage(1) % mu(1,i,j)
            end select
#endif
   
            call BCs(f % zone) % bc % GradVarsForEqn( nEqn,&
                                nGradEqn,                  &
                                f % geom % x(:,i,j),       &
                                time               ,       &
                                f % geom % normal(:,i,j),  &
                                f % storage(1) % Q(:,i,j), &
                                Uhat                         )
   
            f % storage(1) % unStar(:,1,i,j) = Uhat * f % geom % normal(1,i,j) * f % geom % jacobian(i,j)
            f % storage(1) % unStar(:,2,i,j) = Uhat * f % geom % normal(2,i,j) * f % geom % jacobian(i,j)    
            f % storage(1) % unStar(:,3,i,j) = Uhat * f % geom % normal(3,i,j) * f % geom % jacobian(i,j)

         end do ; end do   

      end subroutine BR1_ComputeBoundaryFlux
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine BR1_ComputeInnerFluxes( self , nEqn, nGradEqn, EllipticFlux, GetViscosity, e , contravariantFlux )
         use ElementClass
         use PhysicsStorage
         use Physics
         implicit none
         class(BassiRebay1_t) ,     intent (in) :: self
         integer,                   intent(in)  :: nEqn
         integer,                   intent(in)  :: nGradEqn
         procedure(EllipticFlux_f)              :: EllipticFlux
         procedure(GetViscosity_f)              :: GetViscosity
         type(Element)                          :: e
         real(kind=RP)           , intent (out) :: contravariantFlux(1:nEqn, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3), 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)       :: delta
         real(kind=RP)       :: cartesianFlux(1:nEqn, 1:NDIM)
         real(kind=RP)       :: mu(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: kappa(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: beta(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         integer             :: i, j, k

#if defined(NAVIERSTOKES)
         mu = dimensionless % mu + e % storage % mu_art(1,:,:,:)
         kappa = 1.0_RP / ( thermodynamics % gammaMinus1 * &
                               POW2( dimensionless % Mach) * dimensionless % Pr ) * dimensionless % mu + e % storage % mu_art(3,:,:,:)
         beta  = e % storage % mu_art(2,:,:,:)

#elif defined(INCNS)
         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
            call GetViscosity(e % storage % Q(INSRHO,i,j,k), mu(i,j,k))      
         end do                ; end do                ; end do

         kappa = 0.0_RP
         beta  = 0.0_RP

#elif defined(MULTIPHASE)
         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
            call GetViscosity(e % storage % Q(IMC,i,j,k), mu(i,j,k))      
         end do                ; end do                ; end do

         kappa = 0.0_RP
         beta  = multiphase % M0_star

#endif


         do k = 0, e%Nxyz(3)   ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            call EllipticFlux( nEqn, nGradEqn, e % storage % Q(:,i,j,k) , e % storage % U_x(:,i,j,k) , & 
                               e % storage % U_y(:,i,j,k) , e % storage % U_z(:,i,j,k), mu(i,j,k), beta(i,j,k), kappa(i,j,k), cartesianFlux)
            contravariantFlux(:,i,j,k,IX) =     cartesianFlux(:,IX) * e % geom % jGradXi(IX,i,j,k)  &
                                             +  cartesianFlux(:,IY) * e % geom % jGradXi(IY,i,j,k)  &
                                             +  cartesianFlux(:,IZ) * e % geom % jGradXi(IZ,i,j,k)


            contravariantFlux(:,i,j,k,IY) =     cartesianFlux(:,IX) * e % geom % jGradEta(IX,i,j,k)  &
                                             +  cartesianFlux(:,IY) * e % geom % jGradEta(IY,i,j,k)  &
                                             +  cartesianFlux(:,IZ) * e % geom % jGradEta(IZ,i,j,k)


            contravariantFlux(:,i,j,k,IZ) =     cartesianFlux(:,IX) * e % geom % jGradZeta(IX,i,j,k)  &
                                             +  cartesianFlux(:,IY) * e % geom % jGradZeta(IY,i,j,k)  &
                                             +  cartesianFlux(:,IZ) * e % geom % jGradZeta(IZ,i,j,k)

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
         real(kind=RP) :: cartesianFlux(1:NCONS, 1:NDIM)
         real(kind=RP) :: mu(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP) :: kappa(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP) :: tauSGS(1:NDIM,1:NDIM)
         real(kind=RP) :: qSGS(1:NDIM)
         integer       :: i, j, k

         mu = dimensionless % mu
         kappa = dimensionless % kappa
!
!        Compute subgrid-scale modelling tensor   
!        --------------------------------------
         delta = (e % geom % Volume / product(e % Nxyz + 1)) ** (1.0_RP / 3.0_RP)

         do k = 0, e%Nxyz(3)   ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            call LESModel % ComputeSGSTensor(delta, e % geom % dWall(i,j,k), &
                                                           e % storage % Q(:,i,j,k), &
                                                           e % storage % U_x(:,i,j,k), &
                                                           e % storage % U_y(:,i,j,k), &
                                                           e % storage % U_z(:,i,j,k), &
                                                                tauSGS, qSGS    )

            call ViscousFlux_withSGS(NCONS, NGRAD, e % storage % Q(:,i,j,k) , e % storage % U_x(:,i,j,k) , &
                                     e % storage % U_y(:,i,j,k), e % storage % U_z(:,i,j,k), mu(i,j,k), kappa(i,j,k), tauSGS, qSGS, cartesianFlux)


            contravariantFlux(:,i,j,k,IX) =     cartesianFlux(:,IX) * e % geom % jGradXi(IX,i,j,k)  &
                                             +  cartesianFlux(:,IY) * e % geom % jGradXi(IY,i,j,k)  &
                                             +  cartesianFlux(:,IZ) * e % geom % jGradXi(IZ,i,j,k)


            contravariantFlux(:,i,j,k,IY) =     cartesianFlux(:,IX) * e % geom % jGradEta(IX,i,j,k)  &
                                             +  cartesianFlux(:,IY) * e % geom % jGradEta(IY,i,j,k)  &
                                             +  cartesianFlux(:,IZ) * e % geom % jGradEta(IZ,i,j,k)


            contravariantFlux(:,i,j,k,IZ) =     cartesianFlux(:,IX) * e % geom % jGradZeta(IX,i,j,k)  &
                                             +  cartesianFlux(:,IY) * e % geom % jGradZeta(IY,i,j,k)  &
                                             +  cartesianFlux(:,IZ) * e % geom % jGradZeta(IZ,i,j,k)

         end do               ; end do            ; end do

      end subroutine BR1_ComputeInnerFluxesWithSGS
#endif
      subroutine BR1_RiemannSolver ( self , nEqn, nGradEqn, EllipticFlux, f, QLeft , QRight , U_xLeft , U_yLeft , U_zLeft , U_xRight , U_yRight , U_zRight , &
                                           mu, beta, kappa, nHat , dWall, &
#ifdef MULTIPHASE
sigma, & 
#endif
flux )
         use SMConstants
         use PhysicsStorage
         use Physics
         use FaceClass
         implicit none
         class(BassiRebay1_t)            :: self
         integer,       intent(in)       :: nEqn
         integer,       intent(in)       :: nGradEqn
         procedure(EllipticFlux_f)       :: EllipticFlux
         class(Face),   intent(in)       :: f
         real(kind=RP), intent(in)       :: QLeft(nEqn)
         real(kind=RP), intent(in)       :: QRight(nEqn)
         real(kind=RP), intent(in)       :: U_xLeft(nGradEqn)
         real(kind=RP), intent(in)       :: U_yLeft(nGradEqn)
         real(kind=RP), intent(in)       :: U_zLeft(nGradEqn)
         real(kind=RP), intent(in)       :: U_xRight(nGradEqn)
         real(kind=RP), intent(in)       :: U_yRight(nGradEqn)
         real(kind=RP), intent(in)       :: U_zRight(nGradEqn)
         real(kind=RP), intent(in)       :: mu, beta, kappa
         real(kind=RP), intent(in)       :: nHat(NDIM)
         real(kind=RP), intent(in)       :: dWall
#ifdef MULTIPHASE
         real(kind=RP), intent(in)       :: sigma(nEqn)
#endif
         real(kind=RP), intent(out)      :: flux(nEqn)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)     :: Q(nEqn) , U_x(nGradEqn) , U_y(nGradEqn) , U_z(nGradEqn)
         real(kind=RP)     :: flux_vec(nEqn,NDIM), fL(nEqn,NDIM), fR(nEqn,NDIM)
         real(kind=RP)     :: sigma0

         call EllipticFlux(nEqn, nGradEqn, QLeft, U_xLeft, U_yLeft, U_zLeft, mu, beta, kappa, fL)
         call EllipticFlux(nEqn, nGradEqn, QRight, U_xRight, U_yRight, U_zRight, mu, beta, kappa, fR)

         flux_vec = 0.5_RP * (fL + fR)

         flux = flux_vec(:,IX) * nHat(IX) + flux_vec(:,IY) * nHat(IY) + flux_vec(:,IZ) * nHat(IZ) 

#ifdef MULTIPHASE
         sigma0 = 0.5_RP * self % sigma * (maxval(f % Nf))*(maxval(f % Nf)+1) / f % geom % h
         flux = flux - sigma0 * sigma * (QLeft-QRight)
#endif

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
         call LESModel % ComputeSGSTensor(delta, dWall, Q, U_x, U_y, U_z, tauSGS, qSGS) 

         mu    = dimensionless % mu
         kappa = dimensionless % kappa

         call ViscousFlux_withSGS(NCONS, NGRAD, Q,U_x,U_y,U_z, mu, kappa, tauSGS, qSGS, flux_vec)

         flux = flux_vec(:,IX) * nHat(IX) + flux_vec(:,IY) * nHat(IY) + flux_vec(:,IZ) * nHat(IZ)

      end subroutine BR1_RiemannSolverWithSGS
#endif
end module EllipticBR1
