!
!//////////////////////////////////////////////////////
!
!   @File:    EllipticIP.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Dec 12 13:32:09 2017
!   @Last revision date: Tue May 29 17:43:58 2018
!   @Last revision author: Juan Manzanero (j.manzanero1992@gmail.com)
!   @Last revision commit: 3c1e755ecd17ea60f252dec3daa7823c04603dcd
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module EllipticIP
   use SMConstants
   use Headers
   use MeshTypes
   use ElementClass
   use HexMeshClass
   use Physics
   use PhysicsStorage
   use VariableConversion
   use MPI_Process_Info
   use MPI_Face_Class
   use EllipticDiscretizationClass
   use DGSEMClass
   use FluidData
   implicit none
!
!
   private
   public   InteriorPenalty_t, SIPG, IIPG, NIPG

   integer, parameter   :: SIPG = -1
   integer, parameter   :: IIPG = 0
   integer, parameter   :: NIPG = 1

   type, extends(EllipticDiscretization_t)   :: InteriorPenalty_t
      procedure(PenaltyParameter_f), pointer   :: PenaltyParameter
      real(kind=RP)        :: sigma = 1.0_RP
      integer              :: IPmethod = SIPG
      contains
         procedure      :: Construct              => IP_Construct
         procedure      :: ComputeGradient         => IP_ComputeGradient
         procedure      :: ComputeInnerFluxes      => IP_ComputeInnerFluxes
         procedure      :: RiemannSolver           => IP_RiemannSolver
#if defined(NAVIERSTOKES)
         procedure      :: ComputeInnerFluxesWithSGS => IP_ComputeInnerFluxesWithSGS
         procedure      :: RiemannSolverWithSGS      => IP_RiemannSolverWithSGS
         procedure      :: RiemannSolver_Jacobians => IP_RiemannSolver_Jacobians
#endif
         procedure      :: Describe                => IP_Describe
   end type InteriorPenalty_t

   abstract interface
      function PenaltyParameter_f(self, f)
         use SMConstants
         use FaceClass
         import InteriorPenalty_t
         implicit none
         class(InteriorPenalty_t)   :: self
         class(Face), intent(in)    :: f
         real(kind=RP)              :: PenaltyParameter_f
      end function PenaltyParameter_f
   end interface
!
!  ========
   contains
!  ========
!
      subroutine IP_Construct(self, controlVariables, EllipticFlux0D, EllipticFlux2D, EllipticFlux3D, GetViscosity, eqname)
         use FTValueDictionaryClass
         use Utilities, only: toLower
         use mainKeywordsModule
         use MPI_Process_Info
         use PhysicsStorage
         implicit none
         class(InteriorPenalty_t)              :: self
         class(FTValueDictionary), intent(in)  :: controlVariables
         procedure(EllipticFlux0D_f)           :: EllipticFlux0D
         procedure(EllipticFlux2D_f)           :: EllipticFlux2D
         procedure(EllipticFlux3D_f)           :: EllipticFlux3D
         procedure(GetViscosity_f)             :: GetViscosity
         character(len=*),         intent(in)  :: eqname
!
!        ---------------
!        Local variables
!        ---------------
!
         character(len=LINE_LENGTH)            :: eqnameaux
         character(len=LINE_LENGTH)            :: IPvariant
!
!        ----------------------------------------------------------
!        Set the particular procedures to compute the elliptic flux
!        ----------------------------------------------------------
!
         self % EllipticFlux0D => EllipticFlux0D
         self % EllipticFlux2D => EllipticFlux2D
         self % EllipticFlux3D => EllipticFlux3D
         self % GetViscosity   => GetViscosity

         eqnameaux = eqname
         call tolower(eqnameaux)
         select case(trim(eqnameaux))
         case("ns")
            self % PenaltyParameter => PenaltyParameterNS

         case("ch")
            self % PenaltyParameter => PenaltyParameterCH

         case default
            print*, "Unrecognized equation name."
            print*, "Options available are:"
            print*, "   * NS"
            print*, "   * CH"    
            errorMessage(STD_OUT)
            stop

         end select
!
!        Request the penalty parameter
!        -----------------------------
         if ( controlVariables % containsKey("penalty parameter") ) then
            self % sigma = controlVariables % doublePrecisionValueForKey("penalty parameter")

         else
!            
!           Set 1.0 by default
!           ------------------
            self % sigma = 1.0_RP

         end if
!
!        Request the interior penalty variant
!        ------------------------------------
         if ( controlVariables % containsKey("interior penalty variant") ) then
            IPvariant = controlVariables % stringValueForKey("interior penalty variant", LINE_LENGTH)
            call toLower(IPVariant)
      
         else
!
!           Select SIPG by default
!           ----------------------
            IPvariant = "sipg"
   
         end if

         select case (trim(IPvariant))
         case("sipg")
            self % IPmethod = SIPG

         case("iipg")
            self % IPmethod = IIPG

         case("nipg")
            self % IPmethod = NIPG

         case default
            if ( MPI_Process % isRoot ) then
            print*, "Unknown selected IP variant", trim(IPvariant), "."
            print*, "Available options are:"
            print*, "   * SIPG"
            print*, "   * IIPG"
            print*, "   * NIPG"
            errorMessage(STD_OUT)
            stop
            end if
         end select

      end subroutine IP_Construct

      subroutine IP_Describe(self)
         implicit none
         class(InteriorPenalty_t),  intent(in)  :: self
!
!        Display the configuration
!        -------------------------
         if (MPI_Process % isRoot) write(STD_OUT,'(/)')
         call Subsection_Header("Viscous discretization")

         if (.not. MPI_Process % isRoot ) return

         write(STD_OUT,'(30X,A,A30,A)') "->","Numerical scheme: ","IP"

         select case(self % IPmethod)
         case(SIPG)
            write(STD_OUT,'(30X,A,A30,A)') "->","Interior penalty variant: ","SIPG"

         case(NIPG)
            write(STD_OUT,'(30X,A,A30,A)') "->","Interior penalty variant: ","NIPG"

         case(IIPG)
            write(STD_OUT,'(30X,A,A30,A)') "->","Interior penalty variant: ","IIPG"
            
         end select

         write(STD_OUT,'(30X,A,A30,F10.3)') "->","Penalty parameter: ", self % sigma
            
      end subroutine IP_Describe

      subroutine IP_ComputeGradient(self, nEqn, nGradEqn, mesh, time, externalStateProcedure, GetGradients0D, GetGradients3D)
         use HexMeshClass
         use PhysicsStorage
         use Physics
         use MPI_Process_Info
         implicit none
         class(InteriorPenalty_t), intent(in) :: self
         integer,                  intent(in) :: nEqn, nGradEqn
         class(HexMesh)                       :: mesh
         real(kind=RP),        intent(in)     :: time
         procedure(BCState_FCN)               :: externalStateProcedure
         procedure(GetGradientValues0D_f)     :: GetGradients0D
         procedure(GetGradientValues3D_f)     :: GetGradients3D
!
!        ---------------
!        Local variables
!        ---------------
!
         integer :: Nx, Ny, Nz
         integer :: i, j, k
         integer :: eID , fID , dimID , eqID, fIDs(6)
!
!        *********************************
!        Volume loops and prolong to faces
!        *********************************
!
!$omp do schedule(runtime)
         do eID = 1, size(mesh % elements)
            associate( e => mesh % elements(eID) )
            call e % ComputeLocalGradient(nEqn, nGradEqn, GetGradients3D)
!
!           Prolong to faces
!           ----------------
            fIDs = e % faceIDs
            call e % ProlongGradientsToFaces(nGradEqn, &
                                             mesh % faces(fIDs(1)),&
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
               call IP_GradientInterfaceSolution(f, nEqn, nGradEqn, GetGradients0D) 
            
            case (HMESH_BOUNDARY) 
               call IP_GradientInterfaceSolutionBoundary(f, nEqn, nGradEqn, time, GetGradients0D, externalStateProcedure) 
 
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
            call IP_ComputeGradientFaceIntegrals(self, nGradEqn, e, mesh)
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
            call mesh % GatherMPIFacesSolution(nEqn)
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
               call IP_GradientInterfaceSolutionMPI(f, nEqn, nGradEqn, GetGradients0D) 
 
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
            call IP_ComputeGradientFaceIntegrals(self, nGradEqn, e, mesh)
            end associate
         end do
!$omp end do

      end subroutine IP_ComputeGradient
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine IP_ComputeGradientFaceIntegrals( self, nGradEqn, e, mesh)
         use ElementClass
         use HexMeshClass
         use PhysicsStorage
         use Physics
         use DGIntegrals
         implicit none
         class(InteriorPenalty_t),   intent(in) :: self
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
            invjac = self % IPmethod / e % geom % jacobian(i,j,k)
            e % storage % U_x(:,i,j,k) = e % storage % U_x(:,i,j,k) + faceInt_x(:,i,j,k) * invjac
            e % storage % U_y(:,i,j,k) = e % storage % U_y(:,i,j,k) + faceInt_y(:,i,j,k) * invjac
            e % storage % U_z(:,i,j,k) = e % storage % U_z(:,i,j,k) + faceInt_z(:,i,j,k) * invjac
         end do                  ; end do                   ; end do
!
      end subroutine IP_ComputeGradientFaceIntegrals
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine IP_GradientInterfaceSolution(f, nEqn, nGradEqn, GetGradients)
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

         integer       :: i,j
         
         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            call GetGradients(nEqn, nGradEqn, Q = f % storage(1) % Q(:,i,j), U = UL)
            call GetGradients(nEqn, nGradEqn, Q = f % storage(2) % Q(:,i,j), U = UR)

            Uhat = 0.5_RP * (UL - UR) * f % geom % jacobian(i,j)
            Hflux(:,IX,i,j) = Uhat * f % geom % normal(IX,i,j)
            Hflux(:,IY,i,j) = Uhat * f % geom % normal(IY,i,j)
            Hflux(:,IZ,i,j) = Uhat * f % geom % normal(IZ,i,j)
         end do               ; end do

         call f % ProjectGradientFluxToElements(nGradEqn, HFlux,(/1,2/),1)
         
      end subroutine IP_GradientInterfaceSolution   

      subroutine IP_GradientInterfaceSolutionMPI(f, nEqn, nGradEqn, GetGradients)
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
         integer       :: i,j, thisSide
         
         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            call GetGradients(nEqn, nGradEqn, Q = f % storage(1) % Q(:,i,j), U = UL)
            call GetGradients(nEqn, nGradEqn, Q = f % storage(2) % Q(:,i,j), U = UR)
   
            Uhat = 0.5_RP * (UL - UR) * f % geom % jacobian(i,j)
            Hflux(:,IX,i,j) = Uhat * f % geom % normal(IX,i,j)
            Hflux(:,IY,i,j) = Uhat * f % geom % normal(IY,i,j)
            Hflux(:,IZ,i,j) = Uhat * f % geom % normal(IZ,i,j)
         end do               ; end do

         thisSide = maxloc(f % elementIDs, dim = 1)
         call f % ProjectGradientFluxToElements(nGradEqn, HFlux,(/thisSide, HMESH_NONE/),1)
         
      end subroutine IP_GradientInterfaceSolutionMPI   

      subroutine IP_GradientInterfaceSolutionBoundary(f, nEqn, nGradEqn, time, GetGradients, externalState)
         use Physics
         use FaceClass
         implicit none
         type(Face)                       :: f
         integer,    intent(in)           :: nEqn
         integer,    intent(in)           :: nGradEqn
         real(kind=RP)                    :: time
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
   
            call externalState( nEqn, f % geom % x(:,i,j), &
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
   
            Uhat = 0.5_RP * (UL - UR) * f % geom % jacobian(i,j)
            
            f % storage(1) % unStar(:,1,i,j) = Uhat * f % geom % normal(1,i,j)
            f % storage(1) % unStar(:,2,i,j) = Uhat * f % geom % normal(2,i,j)
            f % storage(1) % unStar(:,3,i,j) = Uhat * f % geom % normal(3,i,j)

         end do ; end do   
         
      end subroutine IP_GradientInterfaceSolutionBoundary
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine IP_ComputeInnerFluxes( self , nEqn, nGradEqn, e , contravariantFlux )
         use ElementClass
         use PhysicsStorage
         use Physics
         implicit none
         class(InteriorPenalty_t) ,     intent(in)  :: self
         integer,                       intent(in)  :: nEqn, nGradEqn
         type(Element)                              :: e
         real(kind=RP)           , intent (out)     :: contravariantFlux(1:nEqn, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3), 1:NDIM)
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

      end subroutine IP_ComputeInnerFluxes
#if defined(NAVIERSTOKES)
      subroutine IP_ComputeInnerFluxesWithSGS( self , e , contravariantFlux )
         use ElementClass
         use PhysicsStorage
         use Physics
         use LESModels
         implicit none
         class(InteriorPenalty_t) ,     intent (in) :: self
         type(Element)                          :: e
         real(kind=RP)           , intent (out) :: contravariantFlux(1:NCONS, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3), 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)       :: delta
         real(kind=RP)       :: cartesianFlux(1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM)
         real(kind=RP)       :: mu(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: kappa(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: tauSGS(1:NDIM,1:NDIM, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: qSGS(1:NDIM, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         integer             :: i, j, k

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

      end subroutine IP_ComputeInnerFluxesWithSGS
#endif
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ----------------------------------------
!     Function to get the IP penalty parameter
!     ----------------------------------------
      function PenaltyParameterNS(self, f)
         use FaceClass
         implicit none
         class(InteriorPenalty_t)   :: self
         class(Face), intent(in)    :: f
         real(kind=RP)              :: PenaltyParameterNS

         PenaltyParameterNS = 0.5_RP*self % sigma * (maxval(f % Nf)+1)*(maxval(f % Nf)+2) / f % geom % h 

      end function PenaltyParameterNS

      function PenaltyParameterCH(self, f)
         use FaceClass
         implicit none
         class(InteriorPenalty_t)   :: self
         class(Face), intent(in)    :: f
         real(kind=RP)              :: PenaltyParameterCH

         PenaltyParameterCH = 0.5_RP*self % sigma * (maxval(f % Nf)+1)*(maxval(f % Nf)+2) / f % geom % h 
         !PenaltyParameterCH = 0.25_RP * self % sigma * (maxval(f % Nf))  *(maxval(f % Nf)+1) / f % geom % h 

      end function PenaltyParameterCH

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine IP_RiemannSolver ( self , nEqn, nGradEqn, f, QLeft , QRight , U_xLeft , U_yLeft , U_zLeft , U_xRight , U_yRight , U_zRight , &
                                           mu, nHat , dWall, flux )
         use SMConstants
         use PhysicsStorage
         use Physics
         use FaceClass
         implicit none
         class(InteriorPenalty_t)        :: self
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
         real(kind=RP)     :: flux_vecL(nEqn,NDIM)
         real(kind=RP)     :: flux_vecR(nEqn,NDIM)
         real(kind=RP)     :: kappa, delta, sigma
         
         sigma = self % PenaltyParameter(f)

#if defined(NAVIERSTOKES)
         kappa = 1.0_RP / ( thermodynamics % gammaMinus1 * &
                            POW2( dimensionless % Mach) * dimensionless % Pr ) * mu
#else
         kappa = 0.0_RP
#endif

         call self % EllipticFlux0D(nEqn, nGradEqn, QLeft , U_xLeft , U_yLeft , U_zLeft, mu, kappa, flux_vecL )
         call self % EllipticFlux0D(nEqn, nGradEqn, QRight , U_xRight , U_yRight , U_zRight, mu, kappa, flux_vecR )

         flux_vec = 0.5_RP * (flux_vecL + flux_vecR)

         flux = flux_vec(:,IX) * nHat(IX) + flux_vec(:,IY) * nHat(IY) + flux_vec(:,IZ) * nHat(IZ) - sigma * mu * (QLeft - QRight)

      end subroutine IP_RiemannSolver
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     -----------------------------------------------------------------------------
!     Subroutine to get the Jacobian (with respect to ∇q⁺ and ∇q⁻) of the numerical
!     contravariant fluxes on a face. Stored in:
!     -> f % storage(side) % dFv_dGradQF(:,:,:,:,i,j)
!                    |                   |_| | | |_|
!                    |                    |  | |  | 
!                    |                    |  | |  |__Coordinate indexes in face 
!                    |                    |  | |_____1 for inner term, 2 for outer term
!                    |                    |  |_______∇q component: 1, 2, 3
!                    |                    |__________Jacobian for this component
!                    |_______________________________1 for ∇q⁺ and 2 for ∇q⁻
!     -----------------------------------------------------------------------------
#if defined(NAVIERSTOKES)
      subroutine IP_RiemannSolver_Jacobians( self, f) 
         use FaceClass
         use Physics
         use PhysicsStorage
         implicit none
         !--------------------------------------------
         class(InteriorPenalty_t), intent(in)    :: self
         type(Face)              , intent(inout) :: f
         !--------------------------------------------
         real(kind=RP), DIMENSION(NCONS,NCONS,NDIM,NDIM) :: df_dgradq   ! Cartesian Jacobian tensor
         real(kind=RP), DIMENSION(NCONS,NCONS,NDIM)      :: dfdq_
         real(kind=RP), parameter :: SideSign(2) = (/ 1._RP, -1._RP /)
         real(kind=RP) :: mu, sigma
         integer :: i,j    ! Face coordinate counters
         integer :: n, m ! Index of G_xx
         integer :: side
         !--------------------------------------------
!
!        Initializations
!        ---------------
         mu    = dimensionless % mu             ! TODO: change for Cahn-Hilliard
         sigma = self % PenaltyParameter(f)
         sigma = sigma * mu
         
         do side = 1, 2
            do j = 0, f % Nf(2) ; do i = 0, f % Nf(1)
!
!           ************************************************
!           Jacobian with respect to ∇q: dF/d∇q⁺ and dF/d∇q⁻
!           ************************************************
!
!            
!              Definitions
!              -----------
               associate( Q             => f % storage(side) % Q  (:,i,j)                , &
                          U_x           => f % storage(side) % U_x(:,i,j)                , &
                          U_y           => f % storage(side) % U_y(:,i,j)                , &
                          U_z           => f % storage(side) % U_z(:,i,j)                , &
                          nHat          => f % geom % normal(:,i,j)                    , &
                          dFStar_dq     => f % storage(side) % dFStar_dqF(:,:,i,j)     , &
                          dF_dGradQ_in  => f % storage(side) % dFv_dGradQF(:,:,:,1,i,j), & 
                          dF_dGradQ_out => f % storage(side) % dFv_dGradQF(:,:,:,2,i,j) )
               
               call ViscousJacobian(Q, U_x, U_y, U_z, (/1._RP, 1._RP, 1._RP/), df_dgradq, dfdq_)
!
!            For the inner surface integral
!            ******************************
               
!              Construct face point Jacobians
!              ------------------------------
               dF_dGradQ_in = 0._RP
               do n = 1, NDIM ; do m = 1, NDIM
                  dF_dGradQ_in(:,:,1) = dF_dGradQ_in(:,:,1) + df_dgradq(:,:,m,n) * f % geom % GradXi  (n,i,j) * nHat(m)
                  dF_dGradQ_in(:,:,2) = dF_dGradQ_in(:,:,2) + df_dgradq(:,:,m,n) * f % geom % GradEta (n,i,j) * nHat(m)
                  dF_dGradQ_in(:,:,3) = dF_dGradQ_in(:,:,3) + df_dgradq(:,:,m,n) * f % geom % GradZeta(n,i,j) * nHat(m)
               end do          ; end do
               
!              Scale according to scheme and multipĺy by the jacobian (surface integral) 
!              -------------------------------------------------------------------------
               dF_dGradQ_in = dF_dGradQ_in * (0.5_RP) * f % geom % jacobian(i,j) ! TODO: Should the constant be -0.5???
               
!               
!            For the outer surface integral
!            ******************************
!            
!              Construct face point Jacobians
!              ------------------------------
               dF_dGradQ_out = 0._RP
               do m = 1, NDIM ; do n = 1, NDIM
                  dF_dGradQ_out(:,:,1) = dF_dGradQ_out(:,:,1) + df_dgradq(:,:,n,m) * f % geom % GradXi  (n,i,j) * nHat(m)
                  dF_dGradQ_out(:,:,2) = dF_dGradQ_out(:,:,2) + df_dgradq(:,:,n,m) * f % geom % GradEta (n,i,j) * nHat(m)
                  dF_dGradQ_out(:,:,3) = dF_dGradQ_out(:,:,3) + df_dgradq(:,:,n,m) * f % geom % GradZeta(n,i,j) * nHat(m)
               end do          ; end do
               
!              Multiply by 1/2 (IP scheme) and the jacobian (surface integral) 
!              ---------------------------------------------------------------
               dF_dGradQ_out = dF_dGradQ_out * 0.5_RP * f % geom % jacobian(i,j)
               
!
!           *********************************************
!           Jacobian with respect to q: dF/dq⁺ and dF/dq⁻
!           *********************************************
!
!
!              Penalty contribution (shifts dFStar_dq matrix)
!              ----------------------------------------------
            
               do n = 1, NCONS
                  dFStar_dq(n,n) = dFStar_dq(n,n) + SideSign(side) * sigma * f % geom % jacobian(i,j)
               end do
               
               end associate
               
            end do              ; end do
            
         end do
         
      end subroutine IP_RiemannSolver_Jacobians
#endif
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
#if defined(NAVIERSTOKES)
      subroutine IP_RiemannSolverWithSGS ( self , f, QLeft , QRight , U_xLeft , U_yLeft , U_zLeft , U_xRight , U_yRight , U_zRight , &
                                            nHat , dWall, flux )
         use SMConstants
         use PhysicsStorage
         use Physics
         use FaceClass
         use LESModels
         implicit none
         class(InteriorPenalty_t)        :: self
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
         real(kind=RP)     :: flux_vecL(NCONS,NDIM), tauSGS_L(NDIM, NDIM), qSGS_L(NDIM)
         real(kind=RP)     :: flux_vecR(NCONS,NDIM), tauSGS_R(NDIM, NDIM), qSGS_R(NDIM)
         real(kind=RP)     :: mu, kappa, delta, sigma
!
!        Compute subgrid-scale modelling tensor   
!        --------------------------------------
         delta = sqrt(f % geom % surface / product(f % Nf + 1))
         
         call LESModel % ComputeSGSTensor(delta, dWall, U_xLeft , U_yLeft , U_zLeft , tauSGS_L, qSGS_L)
         call LESModel % ComputeSGSTensor(delta, dWall, U_xRight, U_yRight, U_zRight, tauSGS_R, qSGS_R) 

         mu    = dimensionless % mu
         kappa = dimensionless % kappa
         
         call ViscousFlux(NCONS, NGRAD, QLeft , U_xLeft , U_yLeft , U_zLeft , mu, kappa, tauSGS_L, qSGS_L, flux_vecL)
         call ViscousFlux(NCONS, NGRAD, QRight, U_xRight, U_yRight, U_zRight, mu, kappa, tauSGS_R, qSGS_R, flux_vecR)
         
         flux_vec = 0.5_RP * (flux_vecL + flux_vecR)
!
!        Shahbazi estimate
!        -----------------
         sigma = 0.5_RP * self % sigma * mu * (maxval(f % Nf)+1)*(maxval(f % Nf)+2) / f % geom % h 

         flux = flux_vec(:,IX) * nHat(IX) + flux_vec(:,IY) * nHat(IY) + flux_vec(:,IZ) * nHat(IZ) - sigma * (QLeft - QRight)

      end subroutine IP_RiemannSolverWithSGS
#endif
end module EllipticIP
