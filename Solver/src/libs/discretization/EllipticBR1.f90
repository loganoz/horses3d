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
         case (ELLIPTIC_NS,ELLIPTIC_NSSA,ELLIPTIC_iNS,ELLIPTIC_MU)
            call Subsection_Header("Viscous discretization")
      
         case (ELLIPTIC_CH)
            call Subsection_Header("Cahn--Hilliard discretization")

         end select
   

         if (.not. MPI_Process % isRoot ) return


         write(STD_OUT,'(30X,A,A30,A)') "->","Numerical scheme: ","BR1"
         write(STD_OUT,'(30X,A,A30,F6.3)') "->","Penalty parameter: ",self % sigma

#ifdef NAVIERSTOKES
         select case (self % eqName)
         case (ELLIPTIC_NS,ELLIPTIC_NSSA)
            select case (grad_vars)
            case(GRADVARS_STATE);   write(STD_OUT,'(30X,A,A30,A)') "->","Gradient variables: ","State"
            case(GRADVARS_ENTROPY); write(STD_OUT,'(30X,A,A30,A)') "->","Gradient variables: ","Entropy"
            case(GRADVARS_ENERGY);  write(STD_OUT,'(30X,A,A30,A)') "->","Gradient variables: ","Energy"
            end select
         end select
#endif

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
         integer                :: eID , fID , dimID , eqID, fIDs(6), iFace, iEl
         logical                :: set_mu
!
!        ************
!        Volume loops
!        ************
!
#ifdef MULTIPHASE
         select case (self % eqName)
         case(ELLIPTIC_MU)
            set_mu = .true.
         case default
            set_mu = .false.
         end select
#else
         set_mu = .false.
#endif
!$omp do schedule(runtime)
         do eID = 1 , size(mesh % elements)
            call mesh % elements(eID) % ComputeLocalGradient(nEqn, nGradEqn, GetGradients, set_mu)
         end do
!$omp end do nowait

         call self % LiftGradients(nEqn, nGradEqn, mesh, time, GetGradients)
   
      end subroutine BR1_ComputeGradient
!
!///////////////////////////////////////////////////////////////////////////////////////////////
!
!        Subroutines for a strong computation of the gradient
!        ----------------------------------------------------
!
!///////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine BR1_LiftGradients(self, nEqn, nGradEqn, mesh, time, GetGradients)
!
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
         integer                :: eID , fID , dimID , eqID, fIDs(6), iFace, iEl
!
!        *******************************************
!        Compute Riemann solvers of non-shared faces
!        *******************************************
!
!$omp do schedule(runtime) private(fID)
         do iFace = 1, size(mesh % faces_interior)
            fID = mesh % faces_interior(iFace)
            call BR1_ComputeElementInterfaceAverage(self, mesh % faces(fID), nEqn, nGradEqn, GetGradients)
         end do
!$omp end do nowait

!$omp do schedule(runtime) private(fID)
         do iFace = 1, size(mesh % faces_boundary)
            fID = mesh % faces_boundary(iFace)
            call BR1_ComputeBoundaryFlux(self, mesh % faces(fID), nEqn, nGradEqn, time, GetGradients)
         end do
!$omp end do 
!
!$omp do schedule(runtime) private(eID)
         do iEl = 1, size(mesh % elements_sequential)
            eID = mesh % elements_sequential(iEl)
            associate(e => mesh % elements(eID))
!
!           Add the surface integrals
!           -------------------------
            call BR1_GradientFaceLoop( self , nGradEqn, e, mesh)
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

!$omp do schedule(runtime) private(fID)
         do iFace = 1, size(mesh % faces_mpi)
            fID = mesh % faces_mpi(iFace)
            call BR1_ComputeMPIFaceAverage(self, mesh % faces(fID), nEqn, nGradEqn, GetGradients)
         end do
!$omp end do 
!
!$omp do schedule(runtime) private(eID) 
         do iEl = 1, size(mesh % elements_mpi)
            eID = mesh % elements_mpi(iEl)
            associate(e => mesh % elements(eID))
!
!           Add the surface integrals
!           -------------------------
            call BR1_GradientFaceLoop(self, nGradEqn, e, mesh)
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

      end subroutine BR1_LiftGradients
!
!//////////////////////////////////////////////////////////////////////////////////////////
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
         integer  :: i,j,k

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
            e % storage % U_x(:,i,j,k) = e % storage % U_x(:,i,j,k) + faceInt_x(:,i,j,k) * e % geom % InvJacobian(i,j,k)
            e % storage % U_y(:,i,j,k) = e % storage % U_y(:,i,j,k) + faceInt_y(:,i,j,k) * e % geom % InvJacobian(i,j,k)
            e % storage % U_z(:,i,j,k) = e % storage % U_z(:,i,j,k) + faceInt_z(:,i,j,k) * e % geom % InvJacobian(i,j,k)
         end do                  ; end do                   ; end do

      end subroutine BR1_GradientFaceLoop
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
         real(kind=RP) :: uStar(nGradEqn)
         real(kind=RP) :: uStar_n(nGradEqn,NDIM,0:f % Nf(1), 0:f % Nf(2))

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

   
            uStar = 0.5_RP * (UR - UL) * f % geom % jacobian(i,j)
            uStar_n(:,IX,i,j) = uStar * f % geom % normal(IX,i,j)
            uStar_n(:,IY,i,j) = uStar * f % geom % normal(IY,i,j)
            uStar_n(:,IZ,i,j) = uStar * f % geom % normal(IZ,i,j)
         end do               ; end do
         
         Sidearray = (/1,2/)
         call f % ProjectGradientFluxToElements(nGradEqn, uStar_n,Sidearray,1)
         
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
         real(kind=RP) :: uStar(nGradEqn)
         real(kind=RP) :: uStar_n(nGradEqn,NDIM,0:f % Nf(1), 0:f % Nf(2))
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

   
            uStar = 0.5_RP * (UR - UL) * f % geom % jacobian(i,j)
            uStar_n(:,IX,i,j) = uStar * f % geom % normal(IX,i,j)
            uStar_n(:,IY,i,j) = uStar * f % geom % normal(IY,i,j)
            uStar_n(:,IZ,i,j) = uStar * f % geom % normal(IZ,i,j)
         end do               ; end do

         Sidearray = (/maxloc(f % elementIDs, dim = 1), HMESH_NONE/)
         call f % ProjectGradientFluxToElements(nGradEqn, uStar_n,Sidearray,1)
         
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
         real(kind=RP) :: u_int(nGradEqn), u_star(nGradEqn)

         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)

#ifdef MULTIPHASE
            call GetGradients(nEqn, nGradEqn, f % storage(1) % Q(:,i,j), u_int, f % storage(1) % rho(i,j))
#else
            call GetGradients(nEqn, nGradEqn, f % storage(1) % Q(:,i,j), u_int)
#endif

#ifdef MULTIPHASE
            select case (self % eqName)
            case (ELLIPTIC_MU)
!
!              The multiphase solver needs the Chemical potential as first entropy variable
!              ----------------------------------------------------------------------------
               u_int(IGMU) = f % storage(1) % mu(1,i,j)
            end select
#endif
   
            u_star = u_int
            call BCs(f % zone) % bc % GradVarsForEqn( nEqn,&
                                nGradEqn,                  &
                                f % geom % x(:,i,j),       &
                                time               ,       &
                                f % geom % normal(:,i,j),  &
                                f % storage(1) % Q(:,i,j), &
                                u_star, GetGradients           )

            f % storage(1) % unStar(:,1,i,j) = (u_star-u_int) * f % geom % normal(1,i,j) * f % geom % jacobian(i,j)
            f % storage(1) % unStar(:,2,i,j) = (u_star-u_int) * f % geom % normal(2,i,j) * f % geom % jacobian(i,j)    
            f % storage(1) % unStar(:,3,i,j) = (u_star-u_int) * f % geom % normal(3,i,j) * f % geom % jacobian(i,j)

         end do ; end do   

      end subroutine BR1_ComputeBoundaryFlux
!
!//////////////////////////////////////////////////////////////////////////////////////////////////
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

#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
         mu = e % storage % mu_ns(1,:,:,:)
         kappa = e % storage % mu_ns(2,:,:,:)
         beta  = 0.0_RP

#elif defined(NAVIERSTOKES) && (SPALARTALMARAS)
         mu    = e % storage % mu_ns(1,:,:,:)
         kappa = e % storage % mu_ns(2,:,:,:)
         beta  = e % storage % mu_ns(3,:,:,:)

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

      subroutine BR1_RiemannSolver ( self , nEqn, nGradEqn, EllipticFlux, f, QLeft , QRight , U_xLeft , U_yLeft , U_zLeft , U_xRight , U_yRight , U_zRight , &
                                           mu_left, mu_right, nHat , dWall, &
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
         real(kind=RP), intent(in)       :: mu_left(3), mu_right(3)
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

         call EllipticFlux(nEqn, nGradEqn, QLeft, U_xLeft, U_yLeft, U_zLeft, mu_left(1), mu_left(2), mu_left(3), fL)
         call EllipticFlux(nEqn, nGradEqn, QRight, U_xRight, U_yRight, U_zRight, mu_right(1), mu_right(2), mu_right(3), fR)

         flux_vec = 0.5_RP * (fL + fR)

         flux = flux_vec(:,IX) * nHat(IX) + flux_vec(:,IY) * nHat(IY) + flux_vec(:,IZ) * nHat(IZ) 

#ifdef MULTIPHASE
         sigma0 = 0.5_RP * self % sigma * (maxval(f % Nf))*(maxval(f % Nf)+1) / f % geom % h
         flux = flux - sigma0 * sigma * (QLeft-QRight)
#endif

      end subroutine BR1_RiemannSolver
end module EllipticBR1