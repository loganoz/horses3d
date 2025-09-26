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
   public   BassiRebay1_t, BR1_RiemannSolver_acc

   type, extends(EllipticDiscretization_t)   :: BassiRebay1_t
      contains
         procedure      :: ComputeGradient           => BR1_ComputeGradient
         procedure      :: LiftGradients             => BR1_LiftGradients
         procedure      :: LiftGradientsHO           => BR1_LiftGradientsHO
         ! procedure      :: ComputeInnerFluxes        => BR1_ComputeInnerFluxes
         procedure      :: RiemannSolver             => BR1_RiemannSolver
         procedure      :: Describe                  => BR1_Describe
         procedure      :: CreateDeviceData          => BR1_CreateDeviceData
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

      subroutine BR1_CreateDeviceData(self)
         implicit none
         class(BassiRebay1_t), intent(in)  :: self

         !$acc update device(self)
         !$acc update device(self % eqName)

      end subroutine BR1_CreateDeviceData

      subroutine BR1_ComputeGradient(self, nEqn, nGradEqn, mesh, time, GetGradients, HO_Elements)
         use HexMeshClass
         use PhysicsStorage
         use Physics
         use ElementClass
         implicit none
         class(BassiRebay1_t), intent(in) :: self
         integer,              intent(in) :: nEqn, nGradEqn
         type(HexMesh),        intent(inout) :: mesh
         real(kind=RP),        intent(in) :: time
         procedure(GetGradientValues_f)   :: GetGradients
         logical, intent(in), optional    :: HO_Elements
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                :: i, j, k
         integer                :: eID , fID , dimID , eqID, fIDs(6), iFace, iEl
         logical                :: set_mu
         logical                :: HOElements


         if (present(HO_Elements)) then
            HOElements = HO_Elements
         else
            HOElements = .false.
         end if
!
!        ************
!        Volume loops
!        ************
!
      !$acc data copyin(self)

      if (HOElements) then

!!$omp do schedule(runtime) private(eID)
!         do i = 1 , size(mesh % HO_Elements)
!            eID = mesh % HO_Elements(i)
!         !   call mesh % elements(eID) % ComputeLocalGradient(nEqn, nGradEqn, GetGradients, set_mu)
!         end do
!!$omp end do nowait
!         call self % LiftGradientsHO(nEqn, nGradEqn, mesh, time, GetGradients)
      else
!!$omp do schedule(runtime)
         !!$acc parallel loop gang present(mesh, mesh % elements) 
         !do eID = 1 , size(mesh % elements)
!        !    call HexElement_ComputeLocalGradient(mesh % elements(eID), NCONS, NGRAD)
         !end do
         !!$acc end parallel loop
!!$omp end do nowait
         call self % LiftGradients(nEqn, nGradEqn, mesh, time, GetGradients)
      end if

      !$acc end data

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
         use ElementClass
         implicit none
         class(BassiRebay1_t), intent(in) :: self
         integer,              intent(in) :: nEqn, nGradEqn
         type(HexMesh),        intent(inout) :: mesh
         real(kind=RP),        intent(in) :: time
         procedure(GetGradientValues_f)   :: GetGradients
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                :: i, j, k
         integer                :: eID , fID , dimID , eqID, fIDs(6), iFace, iEl, nZones, zoneID

!
!        *******************************************
!        Compute Riemann solvers of non-shared faces
!        *******************************************
!
!$omp do schedule(runtime) private(fID)
         !$acc parallel loop gang present(mesh, self) async(1)
         do iFace = 1, size(mesh % faces_interior)
            fID = mesh % faces_interior(iFace)
            call BR1_ComputeElementInterfaceAverage(self, mesh % faces(fID), nEqn, nGradEqn)
         end do
         !$acc end parallel loop
!$omp end do nowait

         nZones = size(mesh % zones)
!$omp do schedule(runtime) private(zoneID)
         do zoneID=1, nZones
            CALL BCs(zoneID) % bc % GradVarsForEqn(mesh, mesh % zones(zoneID)) 
         enddo
!$omp end do 

!$omp do schedule(runtime) private(eID)
         !$acc parallel loop gang present(mesh) async(1) 
         do iEl = 1, size(mesh % elements_sequential)
            eID = mesh % elements_sequential(iEl)
!
!           Add the surface integrals
!           -------------------------
            call BR1_GradientFaceLoop(nGradEqn, mesh % elements(eID), mesh)
         end do
         !$acc end parallel loop
!$omp end do

         call HexMesh_ProlongGradientsToFaces(mesh, size(mesh % elements_sequential), mesh % elements_sequential, nGradEqn)

#ifdef _HAS_MPI_
!$omp single
         if ( MPI_Process % doMPIAction ) then 
            call HexMesh_GatherMPIFacesSolution(mesh, nEqn)
         end if
!$omp end single

!$omp do schedule(runtime) private(fID)
         !$acc parallel loop gang present(mesh, self) async(1) 
         do iFace = 1, size(mesh % faces_mpi)
            fID = mesh % faces_mpi(iFace)
            call BR1_ComputeMPIFaceAverage(self, mesh % faces(fID), nEqn, nGradEqn)
         end do
         !$acc end parallel loop
!$omp end do 
!

!$omp do schedule(runtime) private(eID)
!$acc parallel loop gang vector_length(128) present(mesh, self) async(1)
         do iEl = 1, size(mesh % elements_mpi)
            eID = mesh % elements_mpi(iEl)
!
!           Add the surface integrals
!           -------------------------
            call BR1_GradientFaceLoop(nGradEqn, mesh % elements(eID), mesh)
         end do
!$omp end do
!$acc end parallel loop
!
!           Prolong gradients
!           -----------------
!
         call HexMesh_ProlongGradientsToFaces(mesh, size(mesh % elements_mpi), mesh % elements_mpi, nGradEqn)

#endif

      end subroutine BR1_LiftGradients

      subroutine BR1_LiftGradientsHO(self, nEqn, nGradEqn, mesh, time, GetGradients)
         use HexMeshClass
         use PhysicsStorage
         use Physics
         implicit none
         class(BassiRebay1_t), intent(in) :: self
         integer,              intent(in) :: nEqn, nGradEqn
         type(HexMesh),        intent(inout) :: mesh
         real(kind=RP),        intent(in) :: time
         procedure(GetGradientValues_f)   :: GetGradients
      end subroutine BR1_LiftGradientsHO
!
!//////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine BR1_GradientFaceLoop(nGradEqn, e, mesh )
         !$acc routine vector
         use ElementClass
         use HexMeshClass
         use PhysicsStorage
         use Physics
         use DGIntegrals
         implicit none
         integer,               intent(in)     :: nGradEqn
         type(Element),         intent(inout)  :: e
         type(HexMesh),         intent(inout)  :: mesh

         call  VectorWeakIntegrals_StdFace(e, nGradEqn, &
               mesh % faces(e % faceIDs(EFRONT))  % storage(e % faceSide(EFRONT))  % unStar, &
               mesh % faces(e % faceIDs(EBACK))   % storage(e % faceSide(EBACK))   % unStar, &
               mesh % faces(e % faceIDs(EBOTTOM)) % storage(e % faceSide(EBOTTOM)) % unStar, &
               mesh % faces(e % faceIDs(ERIGHT))  % storage(e % faceSide(ERIGHT))  % unStar, &
               mesh % faces(e % faceIDs(ETOP))    % storage(e % faceSide(ETOP))    % unStar, &
               mesh % faces(e % faceIDs(ELEFT))   % storage(e % faceSide(ELEFT))   % unStar, &
               e % storage % U_x, e % storage % U_y, e % storage % U_z )

      end subroutine BR1_GradientFaceLoop
!
      subroutine BR1_ComputeElementInterfaceAverage(self, f, nEqn, nGradEqn)
         !$acc routine vector
         use Physics  
         use ElementClass
         use FaceClass
         implicit none  
!
!        ---------
!        Arguments
!        ---------
!
         type(BassiRebay1_t),   intent(in)  :: self
         type(Face)                       :: f
         integer,    intent(in)           :: nEqn, nGradEqn
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: UL(NCONS), UR(NCONS)
         real(kind=RP) :: uStar, jacobian

         integer       :: i,j,eq

         !$acc loop vector collapse(2) private(UL, UR, ustar, jacobian)
         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
#ifdef MULTIPHASE
            select case (self % eqName)
               case (ELLIPTIC_MU)
                  call mGradientVariables(nEqn, nGradEqn, f % storage(1) % Q(:,i,j), UL(1:nGradEqn), f % storage(1) % rho(i,j))
                  call mGradientVariables(nEqn, nGradEqn, f % storage(2) % Q(:,i,j), UR(1:nGradEqn), f % storage(2) % rho(i,j))
               !
               !              The multiphase solver needs the Chemical potential as first entropy variable
               !              ----------------------------------------------------------------------------
                  UL(IGMU) = f % storage(1) % mu(1,i,j)
                  UR(IGMU) = f % storage(2) % mu(1,i,j)

               case(ELLIPTIC_CH)
                  call chGradientVariables(nEqn, nGradEqn, f % storage(1) % Q(1:nEqn,i,j), UL(1:nGradEqn))
                  call chGradientVariables(nEqn, nGradEqn, f % storage(2) % Q(1:nEqn,i,j), UR(1:nGradEqn))
            end select
            
#else
            call NSGradientVariables_STATE(nEqn, nGradEqn, f % storage(1) % Q(:,i,j), UL)
            call NSGradientVariables_STATE(nEqn, nGradEqn, f % storage(2) % Q(:,i,j), UR)
#endif

            jacobian = f % geom % jacobian(i,j)

            !$acc loop seq
            do eq =1, nEqn
               uStar = 0.5_RP * (UR(eq) - UL(eq)) * jacobian

               f % storage(1) % unStar(eq,IX,i,j) = uStar * f % geom % normal(IX,i,j)
               f % storage(1) % unStar(eq,IY,i,j) = uStar * f % geom % normal(IY,i,j)
               f % storage(1) % unStar(eq,IZ,i,j) = uStar * f % geom % normal(IZ,i,j)
            enddo
         end do               ; end do
         
         call Face_ProjectGradientFluxToElements(f, nGradEqn, f % storage(1) % unStar, 1,1)
         call Face_ProjectGradientFluxToElements(f, nGradEqn, f % storage(1) % unStar, 2,1)

      end subroutine BR1_ComputeElementInterfaceAverage   

      subroutine BR1_ComputeMPIFaceAverage(self, f, nEqn, nGradEqn)
         !$acc routine vector
         use Physics  
         use ElementClass
         use FaceClass
         implicit none  
!
!        ---------
!        Arguments
!        ---------
!
         type(BassiRebay1_t),   intent(in)  :: self
         type(Face)                       :: f
         integer, intent(in)              :: nEqn, nGradEqn
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: uStar
         integer       :: i,j,eq, Sidearray, maxId
         real(kind=RP) :: UL(NCONS), UR(NCONS)

         !$acc loop vector collapse(2) private(UL,UR)
         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
#ifdef MULTIPHASE
            select case (self % eqName)
               case (ELLIPTIC_MU)
                  call mGradientVariables(nEqn, nGradEqn, f % storage(1) % Q(:,i,j), UL, f % storage(1) % rho(i,j))
                  call mGradientVariables(nEqn, nGradEqn, f % storage(2) % Q(:,i,j), UR, f % storage(2) % rho(i,j))
               !
               !              The multiphase solver needs the Chemical potential as first entropy variable
               !              ----------------------------------------------------------------------------
                  UL(IGMU) = f % storage(1) % mu(1,i,j)
                  UR(IGMU) = f % storage(2) % mu(1,i,j)

               case(ELLIPTIC_CH)
                  call chGradientVariables(nEqn, nGradEqn, f % storage(1) % Q(:,i,j), UL)
                  call chGradientVariables(nEqn, nGradEqn, f % storage(2) % Q(:,i,j), UR)
            end select
            
#else
            call NSGradientVariables_STATE(nEqn, nGradEqn, f % storage(1) % Q(:,i,j), UL)
            call NSGradientVariables_STATE(nEqn, nGradEqn, f % storage(2) % Q(:,i,j), UR)
#endif

            !$acc loop seq
            do eq =1, nEqn
               ! uStar = 0.5_RP * (f % storage(2) % Q(eq,i,j) - f % storage(1) % Q(eq,i,j)) * f % geom % jacobian(i,j)
               
               uStar = 0.5_RP * (UR(eq) - UL(eq)) * f % geom % jacobian(i,j)
               
               f % storage(1) % unStar(eq,IX,i,j) = uStar * f % geom % normal(IX,i,j)
               f % storage(1) % unStar(eq,IY,i,j) = uStar * f % geom % normal(IY,i,j)
               f % storage(1) % unStar(eq,IZ,i,j) = uStar * f % geom % normal(IZ,i,j)
            enddo
         end do               ; end do

         !Code to replace maxloc that is not supported 
         maxId=MAXVAL(f % elementIDs)
         do i=1,SIZE(f % elementIDs)
             if(f % elementIDs(i)==maxId)THEN
               Sidearray=i
                 exit
             endif
         end do
         !Sidearray = MAXLOC(f % elementIDs,1)
         call Face_ProjectGradientFluxToElements(f, nGradEqn, f % storage(1) % unStar, Sidearray,1)
         
      end subroutine BR1_ComputeMPIFaceAverage   
!
!//////////////////////////////////////////////////////////////////////////////////////////////////
!
!
! This subroutine is divided in the new version.
!
!       subroutine BR1_ComputeInnerFluxes( self , nEqn, nGradEqn, EllipticFlux, GetViscosity, e , contravariantFlux )
!          use ElementClass
!          use PhysicsStorage
!          use Physics
!          implicit none
!          class(BassiRebay1_t) ,     intent (in) :: self
!          integer,                   intent(in)  :: nEqn
!          integer,                   intent(in)  :: nGradEqn
!          procedure(EllipticFlux_f)              :: EllipticFlux
!          procedure(GetViscosity_f)              :: GetViscosity
!          type(Element)                          :: e
!          real(kind=RP)           , intent (out) :: contravariantFlux(1:nEqn, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3), 1:NDIM)
! !
! !        ---------------
! !        Local variables
! !        ---------------
! !
!          real(kind=RP)       :: delta
!          real(kind=RP)       :: cartesianFlux(1:nEqn, 1:NDIM)
!          real(kind=RP)       :: mu(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
!          real(kind=RP)       :: kappa(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
!          real(kind=RP)       :: beta(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
!          integer             :: i, j, k

! #if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
!          mu = e % storage % mu_ns(1,:,:,:)
!          kappa = e % storage % mu_ns(2,:,:,:)
!          beta  = 0.0_RP

! #elif defined(NAVIERSTOKES) && (SPALARTALMARAS)
!          mu    = e % storage % mu_ns(1,:,:,:)
!          kappa = e % storage % mu_ns(2,:,:,:)
!          beta  = e % storage % mu_ns(3,:,:,:)

! #elif defined(INCNS)
!          do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
!             call GetViscosity(e % storage % Q(INSRHO,i,j,k), mu(i,j,k))      
!          end do                ; end do                ; end do

!          kappa = 0.0_RP
!          beta  = 0.0_RP

! #elif defined(MULTIPHASE)
!          do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
!             call GetViscosity(e % storage % Q(IMC,i,j,k), mu(i,j,k))      
!          end do                ; end do                ; end do

!          kappa = 0.0_RP
!          beta  = multiphase % M0_star

! #endif


!          do k = 0, e%Nxyz(3)   ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
!             call EllipticFlux( nEqn, nGradEqn, e % storage % Q(:,i,j,k) , e % storage % U_x(:,i,j,k) , & 
!                                e % storage % U_y(:,i,j,k) , e % storage % U_z(:,i,j,k), mu(i,j,k), beta(i,j,k), kappa(i,j,k), cartesianFlux)
            
!             contravariantFlux(:,i,j,k,IX) =     cartesianFlux(:,IX) * e % geom % jGradXi(IX,i,j,k)  &
!                                              +  cartesianFlux(:,IY) * e % geom % jGradXi(IY,i,j,k)  &
!                                              +  cartesianFlux(:,IZ) * e % geom % jGradXi(IZ,i,j,k)


!             contravariantFlux(:,i,j,k,IY) =     cartesianFlux(:,IX) * e % geom % jGradEta(IX,i,j,k)  &
!                                              +  cartesianFlux(:,IY) * e % geom % jGradEta(IY,i,j,k)  &
!                                              +  cartesianFlux(:,IZ) * e % geom % jGradEta(IZ,i,j,k)


!             contravariantFlux(:,i,j,k,IZ) =     cartesianFlux(:,IX) * e % geom % jGradZeta(IX,i,j,k)  &
!                                              +  cartesianFlux(:,IY) * e % geom % jGradZeta(IY,i,j,k)  &
!                                              +  cartesianFlux(:,IZ) * e % geom % jGradZeta(IZ,i,j,k)

!          end do               ; end do            ; end do

!       end subroutine BR1_ComputeInnerFluxes

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

      subroutine BR1_RiemannSolver_acc ( f, nEqn, nGradEqn,&
#ifdef MULTIPHASE
   sigma, BR1_sigma, & 
#endif
   flux ) 
         !$acc routine vector
         use SMConstants
         use PhysicsStorage
         use Physics
         use FaceClass
         implicit none
         type(Face),    intent(in)      :: f
         integer,       intent(in)      :: nEqn
         integer,       intent(in)      :: nGradEqn
         real(kind=RP), intent(out)     :: flux(1:NCONS,0:f% Nf(1),0:f% Nf(2))
#ifdef MULTIPHASE
         real(kind=RP), intent(in)       :: sigma(1:nEqn)
         real(kind=RP), intent(in)       :: BR1_sigma
#endif
         !
         !        ---------------
         !        Local variables
         !        ---------------
         !
         integer :: fID, i, j, eq
         real(kind=RP)     :: sigma0

#ifdef MULTIPHASE
         sigma0 = 0.5_RP * BR1_sigma * (maxval(f % Nf))*(maxval(f % Nf)+1) / f % geom % h
#endif
          
         !$acc loop vector collapse(3)
         do j = 0, f % Nf(2) ;  do i = 0, f % Nf(1) ; do eq = 1, nEqn
               flux (eq,i,j) = 0.5_RP * (f % storage(1) % unStar(eq,IX,i,j) + f % storage(2) % unStar(eq,IX,i,j)) * f % geom % normal(IX,i,j) + &
                               0.5_RP * (f % storage(1) % unStar(eq,IY,i,j) + f % storage(2) % unStar(eq,IY,i,j)) * f % geom % normal(IY,i,j) + &
                               0.5_RP * (f % storage(1) % unStar(eq,IZ,i,j) + f % storage(2) % unStar(eq,IZ,i,j)) * f % geom % normal(IZ,i,j) 
#ifdef MULTIPHASE
               flux (eq,i,j) = flux (eq,i,j) - sigma0 * sigma(eq) * (f % storage(1) % Q(eq,i,j) - f % storage(2) % Q(eq,i,j))
#endif
         enddo ; enddo ; enddo
         
      end subroutine BR1_RiemannSolver_acc

end module EllipticBR1
