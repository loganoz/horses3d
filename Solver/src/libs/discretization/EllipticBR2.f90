!
!//////////////////////////////////////////////////////
!
!   @File:    EllipticBR2.f90
!   @Author:  Juan (juan.manzanero@upm.es)
!   @Created: Fri Dec 15 10:18:31 2017
!   @Last revision date: Sat Aug  3 23:57:14 2019
!   @Last revision author: Andrés Rueda (am.rueda@upm.es)
!   @Last revision commit: 3919d52a3f75c1991f290d63ceec488de9bdd35a
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module EllipticBR2
   use SMConstants
   use Headers
   use MeshTypes
   use ElementClass
   use HexMeshClass
   use PhysicsStorage
   use Physics
   use MPI_Process_Info
   use MPI_Face_Class
   use EllipticDiscretizationClass
   use VariableConversion
   use FluidData
   use BoundaryConditions, only: BCs
   implicit none
!
!
   private
   public   BassiRebay2_t

   type, extends(EllipticDiscretization_t)   :: BassiRebay2_t
      real(kind=RP)        :: eta = 1.0_RP
      contains
         procedure      :: Construct         => BR2_Construct
         procedure      :: ComputeGradient    => BR2_ComputeGradient
         procedure      :: ComputeInnerFluxes => BR2_ComputeInnerFluxes
         procedure      :: RiemannSolver      => BR2_RiemannSolver
#if defined(NAVIERSTOKES)
         procedure      :: ComputeInnerFluxesWithSGS => BR2_ComputeInnerFluxesWithSGS
         procedure      :: RiemannSolverWithSGS      => BR2_RiemannSolverWithSGS
#endif
         procedure      :: Describe           => BR2_Describe
   end type BassiRebay2_t
!
!  ========
   contains
!  ========
!
      subroutine BR2_Construct(self, controlVariables, EllipticFlux0D, EllipticFlux2D, EllipticFlux3D, GetViscosity, eqname)
         use FTValueDictionaryClass
         use mainKeywordsModule
         use MPI_Process_Info
         use PhysicsStorage
         implicit none
         class(BassiRebay2_t)                  :: self
         class(FTValueDictionary), intent(in)  :: controlVariables
         procedure(EllipticFlux0D_f)           :: EllipticFlux0D
         procedure(EllipticFlux2D_f)           :: EllipticFlux2D
         procedure(EllipticFlux3D_f)           :: EllipticFlux3D
         procedure(GetViscosity_f)             :: GetViscosity
         integer,          intent(in)          :: eqname
!
!        ---------------
!        Local variables
!        ---------------
!
         character(len=LINE_LENGTH)            :: BR2variant
!
!        ----------------------------------------------------------
!        Set the particular procedures to compute the elliptic flux
!        ----------------------------------------------------------
!
         self % EllipticFlux0D => EllipticFlux0D
         self % EllipticFlux2D => EllipticFlux2D
         self % EllipticFlux3D => EllipticFlux3D
         self % GetViscosity   => GetViscosity

         select case (eqName)
         case(ELLIPTIC_NS)
            self % eqName = ELLIPTIC_NS

         case(ELLIPTIC_iNS)
            self % eqName = ELLIPTIC_iNS

         case(ELLIPTIC_CH)
            self % eqName = ELLIPTIC_CH

         case(ELLIPTIC_MU)
            self % eqName = ELLIPTIC_MU

         case default
            print*, "Unrecognized equation"
            errorMessage(STD_OUT)
            stop

         end select

!
!        Request the penalty parameter
!        -----------------------------
         if ( controlVariables % containsKey("penalty parameter") ) then
            self % eta = controlVariables % doublePrecisionValueForKey("penalty parameter")

         else
!            
!           Set 3.0 by default
!           ------------------
            self % eta = 2.0_RP

         end if
            
      end subroutine BR2_Construct
   
      subroutine BR2_Describe(self)
         implicit none
         class(BassiRebay2_t),   intent(in)  :: self
!
!        Display the configuration
!        -------------------------
         if (MPI_Process % isRoot) write(STD_OUT,'(/)')
         call Subsection_Header("Viscous discretization")

         if (.not. MPI_Process % isRoot ) return

         write(STD_OUT,'(30X,A,A30,A)') "->","Numerical scheme: ","BR2"

         write(STD_OUT,'(30X,A,A30,F10.3)') "->","Penalty parameter: ", self % eta

      end subroutine BR2_Describe

      subroutine BR2_ComputeGradient( self , nEqn, nGradEqn, mesh , time , GetGradients0D, GetGradients3D)
         use HexMeshClass
         use PhysicsStorage
         use Physics
         use MPI_Process_Info
         implicit none
         class(BassiRebay2_t), intent(in) :: self
         integer,              intent(in) :: nEqn
         integer,              intent(in) :: nGradEqn
         class(HexMesh)                   :: mesh
         real(kind=RP),        intent(in) :: time
         procedure(GetGradientValues0D_f) :: GetGradients0D
         procedure(GetGradientValues3D_f) :: GetGradients3D
!
!        ---------------
!        Local variables
!        ---------------
!
         integer :: Nx, Ny, Nz
         integer :: i, j, k
         integer :: eID , fID , dimID , eqID, fIDs(6)
!
!        ***********************
!        Compute local gradients
!        ***********************
!
!$omp do schedule(runtime)
         do eID = 1, size(mesh % elements)
            associate( e => mesh % elements(eID) )
            call e % ComputeLocalGradient(nEqn, nGradEqn, GetGradients3D)
!
!           Prolong to faces
!           ----------------
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
               call BR2_GradientInterfaceSolution(f, nEqn, nGradEqn, GetGradients0D) 
            
            case (HMESH_BOUNDARY) 
               call BR2_GradientInterfaceSolutionBoundary(f, nEqn, nGradEqn, time, GetGradients0D) 
 
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
            call BR2_ComputeGradientFaceIntegrals(self, nGradEqn, e, mesh)
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
               call BR2_GradientInterfaceSolutionMPI(f, nEqn, nGradEqn, GetGradients0D) 
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
            call BR2_ComputeGradientFaceIntegrals(self, nGradEqn, e, mesh)
            end associate
         end do
!$omp end do

      end subroutine BR2_ComputeGradient
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine BR2_ComputeGradientFaceIntegrals( self, nGradEqn, e, mesh)
!
!        *******************************************************
!              The surface integrals in the BR2 method considers
!           also the correction of the surface gradients with
!           the stabilizing term:
!              \nabla u^* = \nabla u + \eta r_f([[u]])
!
!           where
!              \int_{e} r_f([[u]])\tau = -0.5\int_{f} [[u]] \tau ds
!
!        *******************************************************
!
         use ElementClass
         use HexMeshClass
         use PhysicsStorage
         use Physics
         use DGIntegrals
         use NodalStorageClass, only: NodalStorage
         implicit none
         class(BassiRebay2_t),   intent(in) :: self
         integer,                intent(in) :: nGradEqn
         class(Element)                     :: e
         class(HexMesh)                     :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: i, j, k
         real(kind=RP) :: invjac(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP) :: faceInt_x(nGradEqn, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3) )
         real(kind=RP) :: faceInt_y(nGradEqn, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3) )
         real(kind=RP) :: faceInt_z(nGradEqn, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3) )
         real(kind=RP) :: bv_x(0:e % Nxyz(1),2)
         real(kind=RP) :: bv_y(0:e % Nxyz(2),2)
         real(kind=RP) :: bv_z(0:e % Nxyz(3),2)

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
            invjac(i,j,k) = 1.0_RP / e % geom % jacobian(i,j,k)
            e % storage % U_x(:,i,j,k) = e % storage % U_x(:,i,j,k) - faceInt_x(:,i,j,k) * invjac(i,j,k)
            e % storage % U_y(:,i,j,k) = e % storage % U_y(:,i,j,k) - faceInt_y(:,i,j,k) * invjac(i,j,k)
            e % storage % U_z(:,i,j,k) = e % storage % U_z(:,i,j,k) - faceInt_z(:,i,j,k) * invjac(i,j,k)
         end do                  ; end do                   ; end do
!
!        ******************************************
!        Perform the interface gradients correction
!        ******************************************
!
         associate(spAxi   => NodalStorage(e % Nxyz(1)), &
                   spAeta  => NodalStorage(e % Nxyz(2)), &
                   spAzeta => NodalStorage(e % Nxyz(3)) )
         bv_x = spAxi % b * spAxi % v
         bv_y = spAeta % b * spAeta % v
         bv_z = spAzeta % b * spAzeta % v
         end associate
!
!        ----------------
!>       Xi-contributions
!        ----------------
!

         associate(U_x => mesh % faces(e % faceIDs(ELEFT)) % storage(e % faceSide(ELEFT)) % U_x, &
                   U_y => mesh % faces(e % faceIDs(ELEFT)) % storage(e % faceSide(ELEFT)) % U_y, &
                   U_z => mesh % faces(e % faceIDs(ELEFT)) % storage(e % faceSide(ELEFT)) % U_z, &
                   unStar => mesh % faces(e % faceIDs(ELEFT)) % storage(e % faceSide(ELEFT)) % unStar )

         do k = 0, e%Nxyz(3) ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            U_x(:,j,k) = U_x(:,j,k) - self % eta * unStar(:,1,j,k) * bv_x(i,LEFT) * invjac(i,j,k)
            U_y(:,j,k) = U_y(:,j,k) - self % eta * unStar(:,2,j,k) * bv_x(i,LEFT) * invjac(i,j,k)
            U_z(:,j,k) = U_z(:,j,k) - self % eta * unStar(:,3,j,k) * bv_x(i,LEFT) * invjac(i,j,k)
         end do                 ; end do                ; end do
         end associate

         associate(U_x => mesh % faces(e % faceIDs(ERIGHT)) % storage(e % faceSide(ERIGHT)) % U_x, &
                   U_y => mesh % faces(e % faceIDs(ERIGHT)) % storage(e % faceSide(ERIGHT)) % U_y, &
                   U_z => mesh % faces(e % faceIDs(ERIGHT)) % storage(e % faceSide(ERIGHT)) % U_z, &
                   unStar => mesh % faces(e % faceIDs(ERIGHT)) % storage(e % faceSide(ERIGHT)) % unStar )

         do k = 0, e%Nxyz(3) ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            U_x(:,j,k) = U_x(:,j,k) - self % eta * unStar(:,1,j,k) * bv_x(i,RIGHT) * invjac(i,j,k)
            U_y(:,j,k) = U_y(:,j,k) - self % eta * unStar(:,2,j,k) * bv_x(i,RIGHT) * invjac(i,j,k)
            U_z(:,j,k) = U_z(:,j,k) - self % eta * unStar(:,3,j,k) * bv_x(i,RIGHT) * invjac(i,j,k)
         end do                 ; end do                ; end do
         end associate
!
!        -----------------
!>       Eta-contributions
!        -----------------
!
         associate(U_x => mesh % faces(e % faceIDs(EFRONT)) % storage(e % faceSide(EFRONT)) % U_x, &
                   U_y => mesh % faces(e % faceIDs(EFRONT)) % storage(e % faceSide(EFRONT)) % U_y, &
                   U_z => mesh % faces(e % faceIDs(EFRONT)) % storage(e % faceSide(EFRONT)) % U_z, &
                   unStar => mesh % faces(e % faceIDs(EFRONT)) % storage(e % faceSide(EFRONT)) % unStar )

         do k = 0, e%Nxyz(3) ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            U_x(:,i,k) = U_x(:,i,k) - self % eta * unStar(:,1,i,k) * bv_y(j,LEFT) * invjac(i,j,k)
            U_y(:,i,k) = U_y(:,i,k) - self % eta * unStar(:,2,i,k) * bv_y(j,LEFT) * invjac(i,j,k)
            U_z(:,i,k) = U_z(:,i,k) - self % eta * unStar(:,3,i,k) * bv_y(j,LEFT) * invjac(i,j,k)
         end do                 ; end do                ; end do
         end associate

         associate(U_x => mesh % faces(e % faceIDs(EBACK)) % storage(e % faceSide(EBACK)) % U_x, &
                   U_y => mesh % faces(e % faceIDs(EBACK)) % storage(e % faceSide(EBACK)) % U_y, &
                   U_z => mesh % faces(e % faceIDs(EBACK)) % storage(e % faceSide(EBACK)) % U_z, &
                   unStar => mesh % faces(e % faceIDs(EBACK)) % storage(e % faceSide(EBACK)) % unStar )

         do k = 0, e%Nxyz(3) ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            U_x(:,i,k) = U_x(:,i,k) - self % eta * unStar(:,1,i,k) * bv_y(j,RIGHT) * invjac(i,j,k)
            U_y(:,i,k) = U_y(:,i,k) - self % eta * unStar(:,2,i,k) * bv_y(j,RIGHT) * invjac(i,j,k)
            U_z(:,i,k) = U_z(:,i,k) - self % eta * unStar(:,3,i,k) * bv_y(j,RIGHT) * invjac(i,j,k)
         end do                 ; end do                ; end do
         end associate
!
!        ------------------
!>       Zeta-contributions
!        ------------------
!
         associate(U_x => mesh % faces(e % faceIDs(EBOTTOM)) % storage(e % faceSide(EBOTTOM)) % U_x, &
                   U_y => mesh % faces(e % faceIDs(EBOTTOM)) % storage(e % faceSide(EBOTTOM)) % U_y, &
                   U_z => mesh % faces(e % faceIDs(EBOTTOM)) % storage(e % faceSide(EBOTTOM)) % U_z, &
                   unStar => mesh % faces(e % faceIDs(EBOTTOM)) % storage(e % faceSide(EBOTTOM)) % unStar )

         do k = 0, e%Nxyz(3) ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            U_x(:,i,j) = U_x(:,i,j) - self % eta * unStar(:,1,i,j) * bv_z(k,LEFT) * invjac(i,i,j)
            U_y(:,i,j) = U_y(:,i,j) - self % eta * unStar(:,2,i,j) * bv_z(k,LEFT) * invjac(i,i,j)
            U_z(:,i,j) = U_z(:,i,j) - self % eta * unStar(:,3,i,j) * bv_z(k,LEFT) * invjac(i,i,j)
         end do                 ; end do                ; end do
         end associate

         associate(U_x => mesh % faces(e % faceIDs(ETOP)) % storage(e % faceSide(ETOP)) % U_x, &
                   U_y => mesh % faces(e % faceIDs(ETOP)) % storage(e % faceSide(ETOP)) % U_y, &
                   U_z => mesh % faces(e % faceIDs(ETOP)) % storage(e % faceSide(ETOP)) % U_z, &
                   unStar => mesh % faces(e % faceIDs(ETOP)) % storage(e % faceSide(ETOP)) % unStar )

         do k = 0, e%Nxyz(3) ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            U_x(:,i,j) = U_x(:,i,j) - self % eta * unStar(:,1,i,j) * bv_z(k,RIGHT) * invjac(i,j,k)
            U_y(:,i,j) = U_y(:,i,j) - self % eta * unStar(:,2,i,j) * bv_z(k,RIGHT) * invjac(i,j,k)
            U_z(:,i,j) = U_z(:,i,j) - self % eta * unStar(:,3,i,j) * bv_z(k,RIGHT) * invjac(i,j,k)
         end do                 ; end do                ; end do
         end associate

      end subroutine BR2_ComputeGradientFaceIntegrals
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine BR2_GradientInterfaceSolution(f, nEqn, nGradEqn, GetGradients)
!
!        ************************************************
!           The BR2 is written in strong form, since it
!           is more efficient for the interpolation to
!           boundaries of the local gradients. Hence,
!           the numerical flux is compensated with the
!           interior solution to yield the interface
!           jumps:
!              U_x += -0.5\int_{S} [[u]]\tau ds
!
!           We compute here the interface fluxes:
!              unStar = -0.5*[[u]]dS
!        ************************************************
!
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
         
      end subroutine BR2_GradientInterfaceSolution   

      subroutine BR2_GradientInterfaceSolutionMPI(f, nEqn, nGradEqn, GetGradients)
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
   
            Uhat = 0.5_RP * (UL - UR) * f % geom % jacobian(i,j)
            Hflux(:,IX,i,j) = Uhat * f % geom % normal(IX,i,j)
            Hflux(:,IY,i,j) = Uhat * f % geom % normal(IY,i,j)
            Hflux(:,IZ,i,j) = Uhat * f % geom % normal(IZ,i,j)
         end do               ; end do

         thisSide = maxloc(f % elementIDs, dim = 1)
         call f % ProjectGradientFluxToElements(nGradEqn, HFlux,(/thisSide, HMESH_NONE/),1)
         
      end subroutine BR2_GradientInterfaceSolutionMPI   

      subroutine BR2_GradientInterfaceSolutionBoundary(f, nEqn, nGradEqn, time, GetGradients)
         use Physics
         use FaceClass
         implicit none
         type(Face)                       :: f
         integer,    intent(in)           :: nEqn
         integer,    intent(in)           :: nGradEqn
         real(kind=RP)                    :: time
         procedure(GetGradientValues0D_f) :: GetGradients
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
   
            call BCs(f % zone) % bc % StateForEqn( nEqn, f % geom % x(:,i,j), &
                                time               , &
                                f % geom % normal(:,i,j)      , &
                                bvExt              )
!   
!           -------------------
!           u, v, w, T averages
!           -------------------
!   
            call GetGradients( nEqn, nGradEqn, f % storage(1) % Q(:,i,j), UL )
            call GetGradients( nEqn, nGradEqn, bvExt, UR )
   
            Uhat = 0.5_RP * (UL - UR) * f % geom % jacobian(i,j)
            
            f % storage(1) % unStar(:,1,i,j) = Uhat * f % geom % normal(1,i,j)
            f % storage(1) % unStar(:,2,i,j) = Uhat * f % geom % normal(2,i,j)
            f % storage(1) % unStar(:,3,i,j) = Uhat * f % geom % normal(3,i,j)

         end do ; end do   
         
      end subroutine BR2_GradientInterfaceSolutionBoundary
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine BR2_ComputeInnerFluxes(self, nEqn, nGradEqn, e, contravariantFlux )
         use ElementClass
         use PhysicsStorage
         use Physics
         implicit none
         class(BassiRebay2_t) ,     intent (in) :: self
         integer,                   intent(in)  :: nEqn, nGradEqn
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
         real(kind=RP)       :: beta(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: kappa(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         integer             :: i, j, k

#if (!defined(CAHNHILLIARD))

#if defined(NAVIERSTOKES)
         mu = dimensionless % mu + e % storage % mu_art(1,:,:,:)
         kappa = 1.0_RP / ( thermodynamics % gammaMinus1 * &
                               POW2( dimensionless % Mach) * dimensionless % Pr ) * dimensionless % mu + e % storage % mu_art(3,:,:,:)
         beta  = e % storage % mu_art(2,:,:,:)
#elif defined(INCNS)
         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
            call self % GetViscosity(e % storage % Q(INSRHO,i,j,k), mu(i,j,k))      
         end do                ; end do                ; end do

         kappa = 0.0_RP
         beta  = 0.0_RP
#endif

#else /* !(defined(CAHNHILLIARD) */ 

#if defined(NAVIERSTOKES)
         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
            call self % GetViscosity(e % storage % c(1,i,j,k), mu(i,j,k))      
         end do                ; end do                ; end do
         kappa = 1.0_RP / ( thermodynamics % gammaMinus1 * &
                               POW2( dimensionless % Mach) * dimensionless % Pr ) * mu

         beta  = 0.0_RP
#elif defined(INCNS)
         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
            call self % GetViscosity(e % storage % c(1,i,j,k), mu(i,j,k))      
         end do                ; end do                ; end do

         kappa = 0.0_RP
         beta  = 0.0_RP
#endif
#endif

         call self % EllipticFlux3D(nEqn, nGradEqn, e%Nxyz, e % storage % Q , e % storage % U_x , e % storage % U_y , e % storage % U_z, mu, beta, kappa, cartesianFlux )

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

      end subroutine BR2_ComputeInnerFluxes
#if defined(NAVIERSTOKES)
      subroutine BR2_ComputeInnerFluxesWithSGS( self , e , contravariantFlux )
         use ElementClass
         use PhysicsStorage
         use Physics
         use LESModels
         implicit none
         class(BassiRebay2_t) ,     intent (in) :: self
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
                                                           e % storage % Q, &
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

      end subroutine BR2_ComputeInnerFluxesWithSGS
#endif
      subroutine BR2_RiemannSolver ( self , nEqn, nGradEqn, f, QLeft , QRight , U_xLeft , U_yLeft , U_zLeft , U_xRight , U_yRight , U_zRight , &
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
         class(BassiRebay2_t)            :: self
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
         real(kind=RP)     :: flux_vec(nEqn,NDIM)

!
!>       Old implementation: 1st average, then compute
!        ------------------
         Q   = 0.5_RP * ( QLeft + QRight)
         U_x = 0.5_RP * ( U_xLeft + U_xRight)
         U_y = 0.5_RP * ( U_yLeft + U_yRight)
         U_z = 0.5_RP * ( U_zLeft + U_zRight)

         call self % EllipticFlux0D(nEqn, nGradEqn, Q,U_x,U_y,U_z, mu, beta, kappa, flux_vec)

         flux = flux_vec(:,IX) * nHat(IX) + flux_vec(:,IY) * nHat(IY) + flux_vec(:,IZ) * nHat(IZ)

      end subroutine BR2_RiemannSolver
#if defined(NAVIERSTOKES)
      subroutine BR2_RiemannSolverWithSGS ( self , f, QLeft , QRight , U_xLeft , U_yLeft , U_zLeft , U_xRight , U_yRight , U_zRight , &
                                            nHat , dWall, flux )
         use SMConstants
         use PhysicsStorage
         use Physics
         use FaceClass
         use LESModels
         implicit none
         class(BassiRebay2_t)                 :: self
         class(Face),   intent(in)            :: f
         real(kind=RP), dimension(NCONS)      :: QLeft
         real(kind=RP), dimension(NCONS)      :: QRight
         real(kind=RP), dimension(NGRAD) :: U_xLeft
         real(kind=RP), dimension(NGRAD) :: U_yLeft
         real(kind=RP), dimension(NGRAD) :: U_zLeft
         real(kind=RP), dimension(NGRAD) :: U_xRight
         real(kind=RP), dimension(NGRAD) :: U_yRight
         real(kind=RP), dimension(NGRAD) :: U_zRight
         real(kind=RP), dimension(NDIM)       :: nHat
         real(kind=RP)                        :: dWall
         real(kind=RP), dimension(NCONS)      :: flux
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

         call ViscousFlux(NCONS, NGRAD, Q,U_x,U_y,U_z, mu, kappa, tauSGS, qSGS, flux_vec)

         flux = flux_vec(:,IX) * nHat(IX) + flux_vec(:,IY) * nHat(IY) + flux_vec(:,IZ) * nHat(IZ)

      end subroutine BR2_RiemannSolverWithSGS
#endif
end module EllipticBR2
