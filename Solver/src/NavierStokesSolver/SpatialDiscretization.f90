!
!////////////////////////////////////////////////////////////////////////
!
!      DGTimeDerivativeRoutines.f95
!      Created: 2008-07-13 16:13:12 -0400 
!      By: David Kopriva  
!
!      3D version by D.A. Kopriva 6/17/15, 12:35 PM
!
!
!////////////////////////////////////////////////////////////////////////////////////////
!
#include "Includes.h"
module SpatialDiscretization
      use SMConstants
      use HyperbolicDiscretizations
      use EllipticDiscretizations
      use LESModels
      use SpectralVanishingViscosity
      use DGIntegrals
      use MeshTypes
      use HexMeshClass
      use ElementClass
      use PhysicsStorage
      use Physics
      use MPI_Face_Class
      use MPI_Process_Info
      use DGSEMClass
      use ParticlesClass
      use FluidData
      use VariableConversion, only: NSGradientValuesForQ_0D, NSGradientValuesForQ_3D, GetNSViscosity
      use ProblemFileFunctions, only: UserDefinedSourceTermNS_f
      use BoundaryConditions
#ifdef _HAS_MPI_
      use mpi
#endif

      private
      public   ComputeTimeDerivative, ComputeTimeDerivativeIsolated, viscousDiscretizationKey
      public   Initialize_SpaceAndTimeMethods, Finalize_SpaceAndTimeMethods


      abstract interface
         SUBROUTINE computeElementInterfaceFluxF(f)
            use FaceClass
            IMPLICIT NONE
            TYPE(Face)   , INTENT(inout) :: f   
         end subroutine computeElementInterfaceFluxF

         SUBROUTINE computeMPIFaceFluxF(f)
            use FaceClass
            IMPLICIT NONE
            TYPE(Face)   , INTENT(inout) :: f   
         end subroutine computeMPIFaceFluxF

         SUBROUTINE computeBoundaryFluxF(f, time)
            use SMConstants
            use FaceClass,  only: Face
            IMPLICIT NONE
            type(Face),    intent(inout) :: f
            REAL(KIND=RP)                :: time
         end subroutine computeBoundaryFluxF
      end interface
      
      procedure(computeElementInterfaceFluxF), pointer :: computeElementInterfaceFlux
      procedure(computeMPIFaceFluxF),          pointer :: computeMPIFaceFlux
      procedure(computeBoundaryFluxF),         pointer :: computeBoundaryFlux
   
      character(len=LINE_LENGTH), parameter  :: viscousDiscretizationKey = "viscous discretization"
!
!     ========      
      CONTAINS 
!     ========      
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Initialize_SpaceAndTimeMethods(controlVariables, mesh)
         use FTValueDictionaryClass
         use Utilities, only: toLower
         use mainKeywordsModule
         use Headers
         use MPI_Process_Info
         implicit none
         class(FTValueDictionary),  intent(in)  :: controlVariables
         class(HexMesh)                         :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
         character(len=LINE_LENGTH)       :: inviscidDiscretizationName
         character(len=LINE_LENGTH)       :: viscousDiscretizationName
         
         if (.not. mesh % child) then ! If this is a child mesh, all these constructs were already initialized for the parent mesh
         
            if ( MPI_Process % isRoot ) then
               write(STD_OUT,'(/)')
               call Section_Header("Spatial discretization scheme")
               write(STD_OUT,'(/)')
            end if
   !
   !        Initialize inviscid discretization
   !        ----------------------------------
            inviscidDiscretizationName = controlVariables % stringValueForKey(inviscidDiscretizationKey,requestedLength = LINE_LENGTH)

            call toLower(inviscidDiscretizationName)
         
            select case ( trim(inviscidDiscretizationName) )
            case ( "standard" )
               if (.not. allocated(HyperbolicDiscretization)) allocate( StandardDG_t  :: HyperbolicDiscretization )

            case ( "split-form")
               if (.not. allocated(HyperbolicDiscretization)) allocate( SplitDG_t     :: HyperbolicDiscretization)

            case default
               write(STD_OUT,'(A,A,A)') 'Requested inviscid discretization "',trim(inviscidDiscretizationName),'" is not implemented.'
               write(STD_OUT,'(A)') "Implemented discretizations are:"
               write(STD_OUT,'(A)') "  * Standard"
               write(STD_OUT,'(A)') "  * Split-Form"
               errorMessage(STD_OUT)
               stop 

            end select
               
            call HyperbolicDiscretization % Initialize(controlVariables)
   !
   !        Initialize viscous discretization
   !        ---------------------------------         
            if ( flowIsNavierStokes ) then
               if ( .not. controlVariables % ContainsKey(viscousDiscretizationKey) ) then
                  print*, "Input file is missing entry for keyword: viscous discretization"
                  errorMessage(STD_OUT)
                  stop
               end if

               viscousDiscretizationName = controlVariables % stringValueForKey(viscousDiscretizationKey, requestedLength = LINE_LENGTH)
               call toLower(viscousDiscretizationName)
               
               select case ( trim(viscousDiscretizationName) )
               case("br1")
                  allocate(BassiRebay1_t     :: ViscousDiscretization)

               case("br2")
                  allocate(BassiRebay2_t     :: ViscousDiscretization)

               case("ip")
                  allocate(InteriorPenalty_t :: ViscousDiscretization)

               case default
                  write(STD_OUT,'(A,A,A)') 'Requested viscous discretization "',trim(viscousDiscretizationName),'" is not implemented.'
                  write(STD_OUT,'(A)') "Implemented discretizations are:"
                  write(STD_OUT,'(A)') "  * BR1"
                  write(STD_OUT,'(A)') "  * BR2"
                  write(STD_OUT,'(A)') "  * IP"
                  errorMessage(STD_OUT)
                  stop 

               end select

               call ViscousDiscretization % Construct(controlVariables, ViscousFlux0D, ViscousFlux2D, ViscousFlux3D, GetNSViscosity, "NS")
               call ViscousDiscretization % Describe
      
            else
               if (.not. allocated(ViscousDiscretization)) allocate(EllipticDiscretization_t :: ViscousDiscretization)
               call ViscousDiscretization % Construct(controlVariables, ViscousFlux0D, ViscousFlux2D, ViscousFlux3D, GetNSViscosity, "NS")
               
            end if

   !
   !        Initialize models
   !        -----------------
            call InitializeLESModel(LESModel, controlVariables)
         
         end if
!
!        Initialize SVV
!        --------------
         call InitializeSVV(SVV, controlVariables, mesh)
         
         if (.not. mesh % child) then
            if ( SVV % enabled ) then
               computeElementInterfaceFlux => computeElementInterfaceFlux_SVV
               computeMPIFaceFlux          => computeMPIFaceFlux_SVV
               computeBoundaryFlux         => computeBoundaryFlux_SVV
            else
               computeElementInterfaceFlux => computeElementInterfaceFlux_NS
               computeMPIFaceFlux          => computeMPIFaceFlux_NS
               computeBoundaryFlux         => computeBoundaryFlux_NS
            end if
         end if
         
      end subroutine Initialize_SpaceAndTimeMethods
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine Finalize_SpaceAndTimeMethods
         implicit none
         IF ( ALLOCATED(HyperbolicDiscretization) ) DEALLOCATE( HyperbolicDiscretization )
         IF ( ALLOCATED(LESModel) )       DEALLOCATE( LESModel )
         
         call SVV % destruct()
      end subroutine Finalize_SpaceAndTimeMethods
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ComputeTimeDerivative( mesh, particles, time, mode)
         IMPLICIT NONE 
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(HexMesh), target           :: mesh
         type(Particles_t)               :: particles
         REAL(KIND=RP)                   :: time
         integer, intent(in)             :: mode
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER :: k

         call SetBoundaryConditionsEqn(NS_BC)
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
!$omp parallel shared(mesh, time)
         call mesh % ProlongSolutionToFaces(NCONS)

!        ----------------
!        Update MPI Faces
!        ----------------
!
#ifdef _HAS_MPI_
!$omp single
         call mesh % UpdateMPIFacesSolution(NCONS)
!$omp end single
#endif
!
!        -----------------
!        Compute gradients
!        -----------------
!
         if ( computeGradients ) then
            CALL DGSpatial_ComputeGradient(mesh , time)
         end if

#ifdef _HAS_MPI_
!$omp single
         if ( flowIsNavierStokes ) then
            call mesh % UpdateMPIFacesGradients(NGRAD)
         end if
!$omp end single
#endif
!         call ComputeArtificialViscosity(mesh)
!
!        -----------------------
!        Compute time derivative
!        -----------------------
!
         call TimeDerivative_ComputeQDot(mesh = mesh , &
                                         particles = particles, &
                                         t    = time)
!$omp end parallel
!
      END SUBROUTINE ComputeTimeDerivative
!
!////////////////////////////////////////////////////////////////////////
!
!     This routine computes the time derivative element by element, without considering the Riemann Solvers
!     This is useful for estimating the isolated truncation error
!
      SUBROUTINE ComputeTimeDerivativeIsolated( mesh, particles, time, mode)
         use EllipticDiscretizationClass
         IMPLICIT NONE 
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(HexMesh), target      :: mesh
         type(Particles_t)          :: particles
         REAL(KIND=RP)              :: time
         integer,             intent(in)  :: mode
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER :: k
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
!$omp parallel shared(mesh, time)
         call mesh % ProlongSolutionToFaces(NCONS)
!
!        -----------------------------------------------------
!        Compute LOCAL gradients and prolong them to the faces
!        -----------------------------------------------------
!
         if ( computeGradients ) then
            CALL BaseClass_ComputeGradient( ViscousDiscretization, NCONS, NGRAD, mesh , time , NSGradientValuesForQ_0D, NSGradientValuesForQ_3D )
!
!           The prolongation is usually done in the viscous methods, but not in the BaseClass
!           ---------------------------------------------------------------------------------
            call mesh % ProlongGradientsToFaces(NGRAD)
         end if

!
!        -----------------------
!        Compute time derivative
!        -----------------------
!
         call TimeDerivative_ComputeQDotIsolated(mesh = mesh , &
                                                 t    = time )
!$omp end parallel
!
      END SUBROUTINE ComputeTimeDerivativeIsolated

      subroutine TimeDerivative_ComputeQDot( mesh , particles, t)
         implicit none
         type(HexMesh)              :: mesh
         type(Particles_t)          :: particles
         real(kind=RP)              :: t
         procedure(UserDefinedSourceTermNS_f) :: UserDefinedSourceTermNS
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: eID , i, j, k, ierr, fID
!
!        ****************
!        Volume integrals
!        ****************
!
!$omp do schedule(runtime) 
         do eID = 1 , size(mesh % elements)
            call TimeDerivative_VolumetricContribution( mesh % elements(eID) , t)
         end do
!$omp end do nowait
!
!        ******************************************
!        Compute Riemann solver of non-shared faces
!        ******************************************
!
!$omp do schedule(runtime) 
         do fID = 1, size(mesh % faces) 
            associate( f => mesh % faces(fID)) 
            select case (f % faceType) 
            case (HMESH_INTERIOR) 
               CALL computeElementInterfaceFlux( f ) 
 
            case (HMESH_BOUNDARY) 
               CALL computeBoundaryFlux(f, t) 
 
            end select 
            end associate 
         end do 
!$omp end do 
!
!        ***************************************************************
!        Surface integrals and scaling of elements with non-shared faces
!        ***************************************************************
! 
!$omp do schedule(runtime) private(i,j,k)
         do eID = 1, size(mesh % elements) 
            associate(e => mesh % elements(eID)) 
            if ( e % hasSharedFaces ) cycle
            call TimeDerivative_FacesContribution(e, t, mesh) 
 
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1) 
               e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) / e % geom % jacobian(i,j,k) 
            end do         ; end do          ; end do 
            end associate 
         end do
!$omp end do
!
!        ****************************
!        Wait until messages are sent
!        ****************************
!
#ifdef _HAS_MPI_
         if ( MPI_Process % doMPIAction ) then
!$omp single
            if ( flowIsNavierStokes ) then 
               call mesh % GatherMPIFacesGradients(NGRAD)
            else  
               call mesh % GatherMPIFacesSolution(NCONS)
            end if          
!$omp end single
!
!           **************************************
!           Compute Riemann solver of shared faces
!           **************************************
!
!$omp do schedule(runtime) 
            do fID = 1, size(mesh % faces) 
               associate( f => mesh % faces(fID)) 
               select case (f % faceType) 
               case (HMESH_MPI) 
                  CALL computeMPIFaceFlux( f ) 
               end select 
               end associate 
            end do 
!$omp end do 
!
!           ***********************************************************
!           Surface integrals and scaling of elements with shared faces
!           ***********************************************************
! 
!$omp do schedule(runtime) private(i,j,k)
            do eID = 1, size(mesh % elements) 
               associate(e => mesh % elements(eID)) 
               if ( .not. e % hasSharedFaces ) cycle
               call TimeDerivative_FacesContribution(e, t, mesh) 
   
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1) 
                  e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) / e % geom % jacobian(i,j,k) 
               end do         ; end do          ; end do 
               end associate 
            end do
!$omp end do
!
!           Add an MPI Barrier
!           ------------------
!$omp single
            call mpi_barrier(MPI_COMM_WORLD, ierr)
!$omp end single
         end if
#endif
!
!        *****************************************************************************************************************************
!        Compute contributions to source term
!        ATTENTION: This is deactivated for child multigrid meshes since they have specially tailored source terms (already computed).
!                   If you are going to add contributions to the source term, do it adding to e % storage % S_NS inside the condition!
!        *****************************************************************************************************************************
         if (.not. mesh % child) then
!
!           Add physical source term
!           ************************
!$omp do schedule(runtime) private(i,j,k)
            do eID = 1, mesh % no_of_elements
               associate ( e => mesh % elements(eID) )
               do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  call UserDefinedSourceTermNS(e % geom % x(:,i,j,k), e % storage % Q(:,i,j,k), t, e % storage % S_NS(:,i,j,k), thermodynamics, dimensionless, refValues)
               end do                  ; end do                ; end do
               end associate
            end do
!$omp end do
! May be in the future
! !
! !        Add gravity
! !        ***********
! !
! !$omp do schedule(runtime) private(i,j,k)
!          do eID = 1, mesh % no_of_elements
!             associate(e => mesh % elements(eID))
!             do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
!                e % storage % S_NS(IRHOU:IRHOW,i,j,k) = e % storage % S_NS(IRHOU:IRHOW,i,j,k) + &
!                                                    e % storage % Q(IRHO,i,j,k) * &
!                                  dimensionless % invFr2 * dimensionless % gravity_dir
!             end do                ; end do                ; end do
!             end associate
!          end do
! !$omp end do
!
!           Add Particles source
!           ********************
            if ( particles % active ) then            
!$omp do schedule(runtime)            
               do eID = 1, mesh % no_of_elements
                  associate ( e => mesh % elements(eID) )            
                     e % storage % S_NSP = 0.0_RP
                  end associate
               enddo 
!$omp end do             
               
!$omp do schedule(runtime) 
               do i = 1, particles % injection % injected + 1
                  if (particles % particle(i) % active) then 

                     associate ( eID => particles % particle(i) % eID )
                        
                     call particles % AddSource(i, mesh % elements( eID ), &
                        t, thermodynamics, dimensionless, refValues)

                     ! If this is uncommented, private(j) should be added to openmp.
                        !this commented section permits the computation of source term in neighbor elements
                     !do j = 1, 6
                     !   if (particles % particle(i) % mesh%elements( eID )%NumberOfConnections(j) > 0) then  
                     !      call particles % AddSource(i, &
                     !      mesh % elements( particles % particle(i) % mesh%elements( eID )%Connection(j)%ElementIDs(1)  ), &
                     !      t, thermodynamics, dimensionless, refValues)
                     !   else 
                     !      !
                     !   end if
                     !end do   
                     end associate   
                  endif 
               end do
!$omp end do
!$omp do schedule(runtime)  
               do eID = 1, mesh % no_of_elements
                  associate ( e => mesh % elements(eID) )            
                     e % storage % S_NS = e % storage % S_NS + e % storage % S_NSP
                  end associate
               enddo
!$omp end do
            end if
         end if !(.not. mesh % child)
!
!        ***********************
!        Now add the source term
!        ***********************
!$omp do schedule(runtime) private(i,j,k)
         do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )
            do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) + e % storage % S_NS(:,i,j,k)
            end do                  ; end do                ; end do
            end associate
         end do
!$omp end do

      end subroutine TimeDerivative_ComputeQDot
!
!////////////////////////////////////////////////////////////////////////
!
!     -------------------------------------------------------------------------------
!     This routine computes Qdot neglecting the interaction with neighboring elements
!     and boundaries. Therefore, the external states are not needed.
!     -------------------------------------------------------------------------------
      subroutine TimeDerivative_ComputeQDotIsolated( mesh , t )
         implicit none
         type(HexMesh)              :: mesh
         real(kind=RP)              :: t
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: eID , i, j, k, fID
         procedure(UserDefinedSourceTermNS_f) :: UserDefinedSourceTermNS
!
!        ****************
!        Volume integrals
!        ****************
!
!$omp do schedule(runtime) 
         do eID = 1 , size(mesh % elements)
            call TimeDerivative_StrongVolumetricContribution( mesh % elements(eID) , t)
         end do
!$omp end do
!
!        *******************
!        Scaling of elements
!        *******************
! 
!$omp do schedule(runtime) private(i,j,k)
         do eID = 1, size(mesh % elements) 
            associate(e => mesh % elements(eID))

            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1) 
               e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) / e % geom % jacobian(i,j,k) 
            end do         ; end do          ; end do 
            end associate 
         end do
!$omp end do
!
!        Add a source term
!        -----------------
         if (.not. mesh % child) then
!$omp do schedule(runtime) private(i,j,k)
            do eID = 1, mesh % no_of_elements
               associate ( e => mesh % elements(eID) )
               do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  call UserDefinedSourceTermNS(e % geom % x(:,i,j,k), e % storage % Q(:,i,j,k), t, e % storage % S_NS(:,i,j,k), thermodynamics, dimensionless, refValues)
               end do                  ; end do                ; end do
               end associate
            end do
!$omp end do
         end if
         
!$omp do schedule(runtime) private(i,j,k)
         do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )
            do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) + e % storage % S_NS(:,i,j,k)
            end do                  ; end do                ; end do
            end associate
         end do
!$omp end do
         
      end subroutine TimeDerivative_ComputeQDotIsolated
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine TimeDerivative_StrongVolumetricContribution( e , t )
         use HexMeshClass
         use ElementClass
         implicit none
         type(Element)      :: e
         real(kind=RP)      :: t

!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: inviscidContravariantFlux ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         real(kind=RP) :: fSharp(1:NCONS, 0:e%Nxyz(1), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP) :: gSharp(1:NCONS, 0:e%Nxyz(2), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP) :: hSharp(1:NCONS, 0:e%Nxyz(3), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP) :: viscousContravariantFlux  ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         real(kind=RP) :: SVVContravariantFlux  ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         real(kind=RP) :: contravariantFlux         ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         integer       :: eID
!
!        *************************************
!        Compute interior contravariant fluxes
!        *************************************
!
!        Compute inviscid contravariant flux
!        -----------------------------------
         call HyperbolicDiscretization % ComputeInnerFluxes ( e , EulerFlux3D, inviscidContravariantFlux ) 
!
!        Compute viscous contravariant flux
!        ----------------------------------
         if ( .not. LESModel % active ) then
!
!           Without LES model
!           -----------------
            call ViscousDiscretization  % ComputeInnerFluxes ( NCONS, NGRAD, e , viscousContravariantFlux) 

         else
!
!           With LES model
!           --------------
            call ViscousDiscretization  % ComputeInnerFluxesWithSGS ( e , viscousContravariantFlux  ) 

         end if
!
!        Compute the SVV dissipation
!        ---------------------------
         if ( .not. SVV % enabled ) then
            SVVcontravariantFlux = 0.0_RP
         else
            call SVV % ComputeInnerFluxes(e, ViscousFlux3D, SVVContravariantFlux)
         end if
!
!        ************************
!        Perform volume integrals
!        ************************
!
         select type ( HyperbolicDiscretization )
         type is (StandardDG_t)
!
!           Compute the total Navier-Stokes flux
!           ------------------------------------
            contravariantFlux = inviscidContravariantFlux - viscousContravariantFlux - SVVContravariantFlux
!
!           Perform the Weak Volume Green integral
!           --------------------------------------
            e % storage % QDot = ScalarStrongIntegrals % StdVolumeGreen ( e , NCONS, contravariantFlux ) 

         type is (SplitDG_t)
            ERROR stop ':: TimeDerivative_StrongVolumetricContribution not implemented for split form'
!~ !
!~ !           Compute sharp fluxes for skew-symmetric approximations
!~ !           ------------------------------------------------------
!~             call HyperbolicDiscretization % ComputeSplitFormFluxes(e, inviscidContravariantFlux, fSharp, gSharp, hSharp)
!~ !
!~ !           Peform the Weak volume green integral
!~ !           -------------------------------------
!~             viscousContravariantFlux = viscousContravariantFlux + SVVContravariantFlux

!~             e % storage % QDot = -ScalarWeakIntegrals % SplitVolumeDivergence( e, fSharp, gSharp, hSharp, viscousContravariantFlux)

         end select

      end subroutine TimeDerivative_StrongVolumetricContribution
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine TimeDerivative_VolumetricContribution( e , t )
         use HexMeshClass
         use ElementClass
         implicit none
         type(Element)      :: e
         real(kind=RP)      :: t

!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: inviscidContravariantFlux ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         real(kind=RP) :: fSharp(1:NCONS, 0:e%Nxyz(1), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP) :: gSharp(1:NCONS, 0:e%Nxyz(2), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP) :: hSharp(1:NCONS, 0:e%Nxyz(3), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP) :: viscousContravariantFlux  ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         real(kind=RP) :: SVVContravariantFlux  ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         real(kind=RP) :: contravariantFlux         ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         integer       :: eID
!
!        *************************************
!        Compute interior contravariant fluxes
!        *************************************
!
!        Compute inviscid contravariant flux
!        -----------------------------------
         call HyperbolicDiscretization % ComputeInnerFluxes ( e , EulerFlux3D, inviscidContravariantFlux ) 
!
!        Compute viscous contravariant flux
!        ----------------------------------
         if ( .not. LESModel % active ) then
!
!           Without LES model
!           -----------------
            call ViscousDiscretization  % ComputeInnerFluxes ( NCONS, NGRAD, e , viscousContravariantFlux) 

         else
!
!           With LES model
!           --------------
            call ViscousDiscretization  % ComputeInnerFluxesWithSGS ( e , viscousContravariantFlux  ) 

         end if
!
!        Compute the SVV dissipation
!        ---------------------------
         if ( .not. SVV % enabled ) then
            SVVcontravariantFlux = 0.0_RP
         else
            call SVV % ComputeInnerFluxes(e, ViscousFlux3D, SVVContravariantFlux)
         end if
!
!        ************************
!        Perform volume integrals
!        ************************
!
         select type ( HyperbolicDiscretization )
         type is (StandardDG_t)
!
!           Compute the total Navier-Stokes flux
!           ------------------------------------
            contravariantFlux = inviscidContravariantFlux - viscousContravariantFlux - SVVContravariantFlux
!
!           Perform the Weak Volume Green integral
!           --------------------------------------
            e % storage % QDot = ScalarWeakIntegrals % StdVolumeGreen ( e, NCONS, contravariantFlux ) 

         type is (SplitDG_t)
!
!           Compute sharp fluxes for skew-symmetric approximations
!           ------------------------------------------------------
            call HyperbolicDiscretization % ComputeSplitFormFluxes(e, inviscidContravariantFlux, fSharp, gSharp, hSharp)
!
!           Peform the Weak volume green integral
!           -------------------------------------
            viscousContravariantFlux = viscousContravariantFlux + SVVContravariantFlux

            e % storage % QDot = -ScalarWeakIntegrals % SplitVolumeDivergence( e, fSharp, gSharp, hSharp, viscousContravariantFlux)

         end select

      end subroutine TimeDerivative_VolumetricContribution
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine TimeDerivative_FacesContribution( e , t , mesh)
         use HexMeshClass
         implicit none
         type(Element)           :: e
         real(kind=RP)           :: t
         type(HexMesh)           :: mesh

         e % storage % QDot = e % storage % QDot - ScalarWeakIntegrals % StdFace( e, NCONS, &
                      mesh % faces(e % faceIDs(EFRONT))  % storage(e % faceSide(EFRONT))  % fStar, &
                      mesh % faces(e % faceIDs(EBACK))   % storage(e % faceSide(EBACK))   % fStar, &
                      mesh % faces(e % faceIDs(EBOTTOM)) % storage(e % faceSide(EBOTTOM)) % fStar, &
                      mesh % faces(e % faceIDs(ERIGHT))  % storage(e % faceSide(ERIGHT))  % fStar, &
                      mesh % faces(e % faceIDs(ETOP))    % storage(e % faceSide(ETOP))    % fStar, &
                      mesh % faces(e % faceIDs(ELEFT))   % storage(e % faceSide(ELEFT))   % fStar )

      end subroutine TimeDerivative_FacesContribution
!
!///////////////////////////////////////////////////////////////////////////////////////////// 
! 
!        Riemann solver drivers 
!        ---------------------- 
! 
!///////////////////////////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE computeElementInterfaceFlux_NS(f)
         use FaceClass
         use RiemannSolvers_NS
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j
         real(kind=RP) :: inv_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: mu, beta, kappa
         integer       :: Sidearray(2)

         if ( .not. LESModel % active ) then
         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)

               mu    = dimensionless % mu + f % storage(1) % mu_art(1,i,j)
               beta  = f % storage(1) % mu_art(2,i,j)
               kappa = dimensionless % kappa + f % storage(1) % mu_art(3,i,j)
!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               CALL ViscousDiscretization % RiemannSolver(nEqn = NCONS, nGradEqn = NGRAD, &
                                                  f = f, &
                                                  QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                  U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                  U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                  U_xRight = f % storage(2) % U_x(:,i,j), &
                                                  U_yRight = f % storage(2) % U_y(:,i,j), &
                                                  U_zRight = f % storage(2) % U_z(:,i,j), &
                                                  mu   = mu, beta = beta, kappa = kappa, &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  dWall = f % geom % dWall(i,j), &
                                                  flux  = visc_flux(:,i,j) )

            end do
         end do

         else
         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)

!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               CALL ViscousDiscretization % RiemannSolverWithSGS(f = f, &
                                                  QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                  U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                  U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                  U_xRight = f % storage(2) % U_x(:,i,j), &
                                                  U_yRight = f % storage(2) % U_y(:,i,j), &
                                                  U_zRight = f % storage(2) % U_z(:,i,j), &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  dWall = f % geom % dWall(i,j), &
                                                  flux  = visc_flux(:,i,j) )

            end do
         end do
         end if

         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)
!      
!              --------------
!              Invscid fluxes
!              --------------
!      
               CALL RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                                  QRight = f % storage(2) % Q(:,i,j), &
                                  nHat   = f % geom % normal(:,i,j), &
                                  t1     = f % geom % t1(:,i,j), &
                                  t2     = f % geom % t2(:,i,j), &
                                  flux   = inv_flux(:,i,j) )

               
!
!              Multiply by the Jacobian
!              ------------------------
               flux(:,i,j) = ( inv_flux(:,i,j) - visc_flux(:,i,j)) * f % geom % jacobian(i,j)
               
            END DO   
         END DO  
!
!        ---------------------------
!        Return the flux to elements
!        ---------------------------
!
         Sidearray = (/1,2/)
         call f % ProjectFluxToElements(NCONS, flux, Sidearray)

      END SUBROUTINE computeElementInterfaceFlux_NS

      SUBROUTINE computeMPIFaceFlux_NS(f)
         use FaceClass
         use RiemannSolvers_NS
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j
         integer       :: thisSide
         real(kind=RP) :: inv_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: mu, kappa, beta
         integer       :: Sidearray(2)

         if ( .not. LESModel % active ) then
         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)

               mu    = dimensionless % mu + f % storage(1) % mu_art(1,i,j)
               beta  = f % storage(1) % mu_art(2,i,j)
               kappa = dimensionless % kappa + f % storage(1) % mu_art(3,i,j)
!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               CALL ViscousDiscretization % RiemannSolver(nEqn = NCONS, nGradEqn = NGRAD, &
                                                  f = f, &
                                                  QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                  U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                  U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                  U_xRight = f % storage(2) % U_x(:,i,j), &
                                                  U_yRight = f % storage(2) % U_y(:,i,j), &
                                                  U_zRight = f % storage(2) % U_z(:,i,j), &
                                                  mu   = mu, &
                                                  beta = beta, &
                                                  kappa= kappa, &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  dWall = f % geom % dWall(i,j), &
                                                  flux  = visc_flux(:,i,j) )

            end do
         end do

         else
         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)

!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               CALL ViscousDiscretization % RiemannSolverWithSGS(f = f, &
                                                  QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                  U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                  U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                  U_xRight = f % storage(2) % U_x(:,i,j), &
                                                  U_yRight = f % storage(2) % U_y(:,i,j), &
                                                  U_zRight = f % storage(2) % U_z(:,i,j), &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  dWall = f % geom % dWall(i,j), &
                                                  flux  = visc_flux(:,i,j) )

            end do
         end do
         end if

         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)
!      
!              --------------
!              Invscid fluxes
!              --------------
!      
               CALL RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                                  QRight = f % storage(2) % Q(:,i,j), &
                                  nHat   = f % geom % normal(:,i,j), &
                                  t1     = f % geom % t1(:,i,j), &
                                  t2     = f % geom % t2(:,i,j), &
                                  flux   = inv_flux(:,i,j) )

               
!
!              Multiply by the Jacobian
!              ------------------------
               flux(:,i,j) = ( inv_flux(:,i,j) - visc_flux(:,i,j)) * f % geom % jacobian(i,j)
               
            END DO   
         END DO  
!
!        ---------------------------
!        Return the flux to elements: The sign in eR % storage % FstarB has already been accouted.
!        ---------------------------
!
         thisSide = maxloc(f % elementIDs, dim = 1)
         
         Sidearray = (/thisSide, HMESH_NONE/)
         call f % ProjectFluxToElements(NCONS, flux, Sidearray )

      end subroutine ComputeMPIFaceFlux_NS

      SUBROUTINE computeBoundaryFlux_NS(f, time)
      USE ElementClass
      use FaceClass
      USE RiemannSolvers_NS
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      type(Face),    intent(inout) :: f
      REAL(KIND=RP)                :: time
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER                         :: i, j
      INTEGER, DIMENSION(2)           :: N
      REAL(KIND=RP)                   :: inv_flux(NCONS)
      real(kind=RP)                   :: visc_flux(NCONS, 0:f % Nf(1), 0:f % Nf(2))
      real(kind=RP)                   :: fStar(NCONS, 0:f % Nf(1), 0: f % Nf(2))
      real(kind=RP)                   :: mu, beta, kappa
      integer       :: Sidearray(2)
!
!     -------------------
!     Get external states
!     -------------------
!
      do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
         f % storage(2) % Q(:,i,j) = f % storage(1) % Q(:,i,j)
         CALL BCs(f % zone) % bc % FlowState(f % geom % x(:,i,j), &
                                      time, &
                                      f % geom % normal(:,i,j), &
                                      f % storage(2) % Q(:,i,j))

      end do               ; end do


      if ( flowIsNavierStokes ) then
         if ( .not. LESModel % active ) then
         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)
            f % storage(2) % U_x(:,i,j) = f % storage(1) % U_x(:,i,j)
            f % storage(2) % U_y(:,i,j) = f % storage(1) % U_y(:,i,j)
            f % storage(2) % U_z(:,i,j) = f % storage(1) % U_z(:,i,j)

            CALL BCs(f % zone) % bc % FlowNeumann(&
                                              f % geom % x(:,i,j), &
                                              time, &
                                              f % geom % normal(:,i,j), &
                                              f % storage(2) % Q(:,i,j), &
                                              f % storage(2) % U_x(:,i,j), &
                                              f % storage(2) % U_y(:,i,j), &
                                              f % storage(2) % U_z(:,i,j))

               mu    = dimensionless % mu + f % storage(1) % mu_art(1,i,j)
               beta  = f % storage(1) % mu_art(2,i,j)
               kappa = dimensionless % kappa + f % storage(1) % mu_art(3,i,j)
!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               CALL ViscousDiscretization % RiemannSolver(nEqn = NCONS, nGradEqn = NGRAD, &
                                                  f = f, &
                                                  QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                  U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                  U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                  U_xRight = f % storage(2) % U_x(:,i,j), &
                                                  U_yRight = f % storage(2) % U_y(:,i,j), &
                                                  U_zRight = f % storage(2) % U_z(:,i,j), &
                                                  mu   = mu, beta = beta, kappa = kappa, &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  dWall = f % geom % dWall(i,j), &
                                                  flux  = visc_flux(:,i,j) )

            end do
         end do

         else
         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)
               f % storage(2) % U_x(:,i,j) = f % storage(1) % U_x(:,i,j)
               f % storage(2) % U_y(:,i,j) = f % storage(1) % U_y(:,i,j)
               f % storage(2) % U_z(:,i,j) = f % storage(1) % U_z(:,i,j)
   
               CALL BCs(f % zone) % bc % FlowNeumann(&
                                                 f % geom % x(:,i,j), &
                                                 time, &
                                                 f % geom % normal(:,i,j), &
                                                 f % storage(2) % Q(:,i,j), &
                                                 f % storage(2) % U_x(:,i,j), &
                                                 f % storage(2) % U_y(:,i,j), &
                                                 f % storage(2) % U_z(:,i,j))

!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               CALL ViscousDiscretization % RiemannSolverWithSGS(f = f, &
                                                  QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                  U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                  U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                  U_xRight = f % storage(2) % U_x(:,i,j), &
                                                  U_yRight = f % storage(2) % U_y(:,i,j), &
                                                  U_zRight = f % storage(2) % U_z(:,i,j), &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  dWall = f % geom % dWall(i,j), &
                                                  flux  = visc_flux(:,i,j) )

            end do
         end do
         end if
      else
         visc_flux = 0.0_RP

      end if

      DO j = 0, f % Nf(2)
         DO i = 0, f % Nf(1)
!
!           Hyperbolic part
!           -------------
            CALL RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                               QRight = f % storage(2) % Q(:,i,j), &   
                               nHat   = f % geom % normal(:,i,j), &
                               t1     = f % geom % t1(:,i,j), &
                               t2     = f % geom % t2(:,i,j), &
                               flux   = inv_flux)

            fStar(:,i,j) = (inv_flux - visc_flux(:,i,j) ) * f % geom % jacobian(i,j)
         END DO   
      END DO   

      Sidearray = (/1, HMESH_NONE/)
      call f % ProjectFluxToElements(NCONS, fStar, Sidearray)

      END SUBROUTINE computeBoundaryFlux_NS

      SUBROUTINE computeElementInterfaceFlux_SVV(f)
         use FaceClass
         use RiemannSolvers_NS
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j
         real(kind=RP) :: inv_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: SVV_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: mu, beta, kappa
         integer       :: Sidearray(2)
!
!        ----------
!        SVV fluxes
!        ----------
!
         if ( SVV % enabled ) then 
         CALL SVV % RiemannSolver(f = f, &
                        EllipticFlux = ViscousFlux2D, &
                              QLeft = f % storage(1) % Q, &
                             QRight = f % storage(2) % Q, &
                            U_xLeft = f % storage(1) % U_x, &
                            U_yLeft = f % storage(1) % U_y, &
                            U_zLeft = f % storage(1) % U_z, &
                           U_xRight = f % storage(2) % U_x, &
                           U_yRight = f % storage(2) % U_y, &
                           U_zRight = f % storage(2) % U_z, &
                              flux  = SVV_flux               )
         else
            SVV_Flux = 0.0_RP

         end if
!
         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)

!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               mu    = dimensionless % mu + f % storage(1) % mu_art(1,i,j)
               beta  = f % storage(1) % mu_art(2,i,j)
               kappa = dimensionless % kappa + f % storage(1) % mu_art(3,i,j)

               CALL ViscousDiscretization % RiemannSolver(nEqn = NCONS, nGradEqn = NGRAD, &
                                                  f = f, &
                                                  QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                  U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                  U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                  U_xRight = f % storage(2) % U_x(:,i,j), &
                                                  U_yRight = f % storage(2) % U_y(:,i,j), &
                                                  U_zRight = f % storage(2) % U_z(:,i,j), &
                                                  mu = mu, beta = beta, kappa = kappa, &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  dWall = f % geom % dWall(i,j), &
                                                  flux  = visc_flux(:,i,j) )

!      
!              --------------
!              Invscid fluxes
!              --------------
!      
               CALL RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                                  QRight = f % storage(2) % Q(:,i,j), &
                                  nHat   = f % geom % normal(:,i,j), &
                                  t1     = f % geom % t1(:,i,j), &
                                  t2     = f % geom % t2(:,i,j), &
                                  flux   = inv_flux(:,i,j) )

               
!
!              Multiply by the Jacobian
!              ------------------------
               flux(:,i,j) = ( inv_flux(:,i,j) - visc_flux(:,i,j) - SVV_flux(:,i,j) ) * f % geom % jacobian(i,j)
               
            END DO   
         END DO  
!
!        ---------------------------
!        Return the flux to elements
!        ---------------------------
!
         Sidearray = (/1,2/)
         call f % ProjectFluxToElements(NCONS, flux, Sidearray )

      END SUBROUTINE computeElementInterfaceFlux_SVV

      SUBROUTINE computeMPIFaceFlux_SVV(f)
         use FaceClass
         use RiemannSolvers_NS
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j
         integer       :: thisSide
         real(kind=RP) :: inv_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: SVV_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: mu, beta, kappa
         integer       :: Sidearray(2)
!
!        ----------
!        SVV fluxes
!        ----------
!
         CALL SVV % RiemannSolver(f = f, &
                            EllipticFlux = ViscousFlux2D, &
                              QLeft = f % storage(1) % Q, &
                             QRight = f % storage(2) % Q, &
                            U_xLeft = f % storage(1) % U_x, &
                            U_yLeft = f % storage(1) % U_y, &
                            U_zLeft = f % storage(1) % U_z, &
                           U_xRight = f % storage(2) % U_x, &
                           U_yRight = f % storage(2) % U_y, &
                           U_zRight = f % storage(2) % U_z, &
                              flux  = SVV_flux               )

!
!        --------------
!        Invscid fluxes
!        --------------
!
         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)
!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               mu    = dimensionless % mu + f % storage(1) % mu_art(1,i,j)
               beta  = f % storage(1) % mu_art(2,i,j)
               kappa = dimensionless % kappa + f % storage(1) % mu_art(3,i,j)

               CALL ViscousDiscretization % RiemannSolver(nEqn = NCONS, nGradEqn = NGRAD, &
                                                  f = f, &
                                                  QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                  U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                  U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                  U_xRight = f % storage(2) % U_x(:,i,j), &
                                                  U_yRight = f % storage(2) % U_y(:,i,j), &
                                                  U_zRight = f % storage(2) % U_z(:,i,j), &
                                                  mu = mu, beta = beta, kappa = kappa, &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  dWall = f % geom % dWall(i,j), &
                                                  flux  = visc_flux(:,i,j) )

!
               CALL RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                                  QRight = f % storage(2) % Q(:,i,j), &
                                  nHat   = f % geom % normal(:,i,j), &
                                  t1     = f % geom % t1(:,i,j), &
                                  t2     = f % geom % t2(:,i,j), &
                                  flux   = inv_flux(:,i,j) )
!
!              Multiply by the Jacobian
!              ------------------------
               flux(:,i,j) = ( inv_flux(:,i,j) - visc_flux(:,i,j) - SVV_flux(:,i,j) ) * f % geom % jacobian(i,j)
               
            END DO   
         END DO  
!
!        ---------------------------
!        Return the flux to elements: The sign in eR % storage % FstarB has already been accouted.
!        ---------------------------
!
         thisSide = maxloc(f % elementIDs, dim = 1)
         
         Sidearray = (/thisSide, HMESH_NONE/)
         call f % ProjectFluxToElements(NCONS, flux, Sidearray)


      end subroutine ComputeMPIFaceFlux_SVV

      SUBROUTINE computeBoundaryFlux_SVV(f, time)
      USE ElementClass
      use FaceClass
      USE RiemannSolvers_NS
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      type(Face),    intent(inout) :: f
      REAL(KIND=RP)                :: time
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER                         :: i, j
      INTEGER, DIMENSION(2)           :: N
      REAL(KIND=RP)                   :: inv_flux(NCONS)
      real(kind=RP)                   :: visc_flux(NCONS, 0:f % Nf(1), 0:f % Nf(2))
      real(kind=RP)                   :: SVV_flux(NCONS, 0:f % Nf(1), 0:f % Nf(2))
      real(kind=RP)                   :: fStar(NCONS, 0:f % Nf(1), 0: f % Nf(2))
      real(kind=RP)                   :: mu, beta, kappa
      integer       :: Sidearray(2)
!
!     -------------------
!     Get external states
!     -------------------
!
      do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
         f % storage(2) % Q(:,i,j) = f % storage(1) % Q(:,i,j)
         CALL BCs(f % zone) % bc % FlowState( &
                                      f % geom % x(:,i,j), &
                                      time, &
                                      f % geom % normal(:,i,j), &
                                      f % storage(2) % Q(:,i,j))

      end do               ; end do

      if ( flowIsNavierStokes ) then
         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            f % storage(2) % U_x(:,i,j) = f % storage(1) % U_x(:,i,j)
            f % storage(2) % U_y(:,i,j) = f % storage(1) % U_y(:,i,j)
            f % storage(2) % U_z(:,i,j) = f % storage(1) % U_z(:,i,j)

            CALL BCs(f % zone) % bc % FlowNeumann(&
                                              f % geom % x(:,i,j), &
                                              time, &
                                              f % geom % normal(:,i,j), &
                                              f % storage(2) % Q(:,i,j), &
                                              f % storage(2) % U_x(:,i,j), &
                                              f % storage(2) % U_y(:,i,j), &
                                              f % storage(2) % U_z(:,i,j))

!   
!           --------------
!           Viscous fluxes
!           --------------
!   
            mu    = dimensionless % mu + f % storage(1) % mu_art(1,i,j)
            beta  = f % storage(1) % mu_art(2,i,j)
            kappa = dimensionless % kappa + f % storage(1) % mu_art(3,i,j)

            CALL ViscousDiscretization % RiemannSolver(nEqn = NCONS, nGradEqn = NGRAD, &
                                               f = f, &
                                               QLeft = f % storage(1) % Q(:,i,j), &
                                               QRight = f % storage(2) % Q(:,i,j), &
                                               U_xLeft = f % storage(1) % U_x(:,i,j), &
                                               U_yLeft = f % storage(1) % U_y(:,i,j), &
                                               U_zLeft = f % storage(1) % U_z(:,i,j), &
                                               U_xRight = f % storage(2) % U_x(:,i,j), &
                                               U_yRight = f % storage(2) % U_y(:,i,j), &
                                               U_zRight = f % storage(2) % U_z(:,i,j), &
                                               mu = mu, beta = beta, kappa = kappa, &
                                               nHat = f % geom % normal(:,i,j) , &
                                               dWall = f % geom % dWall(i,j), &
                                               flux  = visc_flux(:,i,j) )

         end do               ; end do
      else
         visc_flux = 0.0_RP

      end if
!
!     ----------
!     SVV fluxes
!     ----------
!
      CALL SVV % RiemannSolver(f = f, &
                           EllipticFlux = ViscousFlux2D, &
                           QLeft = f % storage(1) % Q, &
                          QRight = f % storage(2) % Q, &
                         U_xLeft = f % storage(1) % U_x, &
                         U_yLeft = f % storage(1) % U_y, &
                         U_zLeft = f % storage(1) % U_z, &
                        U_xRight = f % storage(2) % U_x, &
                        U_yRight = f % storage(2) % U_y, &
                        U_zRight = f % storage(2) % U_z, &
                           flux  = SVV_flux               )

      DO j = 0, f % Nf(2)
         DO i = 0, f % Nf(1)
!
!           Hyperbolic part
!           -------------
            CALL RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                               QRight = f % storage(2) % Q(:,i,j), &   
                               nHat   = f % geom % normal(:,i,j), &
                               t1     = f % geom % t1(:,i,j), &
                               t2     = f % geom % t2(:,i,j), &
                               flux   = inv_flux)

            fStar(:,i,j) = (inv_flux - visc_flux(:,i,j) - SVV_flux(:,i,j) ) * f % geom % jacobian(i,j)
         END DO   
      END DO   

      Sidearray = (/1, HMESH_NONE/)
      call f % ProjectFluxToElements(NCONS, fStar, Sidearray)

      END SUBROUTINE computeBoundaryFlux_SVV
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!              GRADIENT PROCEDURES
!              -------------------
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine DGSpatial_ComputeGradient( mesh , time)
         use HexMeshClass
         implicit none
         type(HexMesh)                  :: mesh
         real(kind=RP),      intent(in) :: time

         call ViscousDiscretization % ComputeGradient( NCONS, NGRAD, mesh , time, NSGradientValuesForQ_0D, NSGradientValuesForQ_3D)

      end subroutine DGSpatial_ComputeGradient

      subroutine ComputeArtificialViscosity(mesh)
         use ArtificialViscosity
         implicit none
         class(HexMesh)    :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: eID, i, j, k, fID
         real(kind=RP) :: h, sBeta, muAver(3)
         real(kind=RP), parameter   :: kBeta = 1.5_RP

!$omp do private(h,i,j,k,fID,muAver,sBeta)
         do eID = 1, mesh % no_of_elements
            associate(e => mesh % elements(eID))
            h = (e % geom % volume/product(e % Nxyz+1))**(1.0_RP/3.0_RP) 
!
!           Compute Pointwise viscosity
!           ---------------------------
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               call ComputeShockSensor(e % storage % Q(:,i,j,k), e % storage % U_x(:,i,j,k), &
                                       e % storage % U_y(:,i,j,k) , e % storage % U_z(:,i,j,k), h, sBeta)

               e % storage % mu_art(1,i,j,k) = e % storage % Q(IRHO,i,j,k) * h * kBeta * sBeta
               e % storage % mu_art(2,i,j,k) = e % storage % Q(IRHO,i,j,k) * h * kBeta * sBeta * 0.1_RP
               e % storage % mu_art(3,i,j,k) = e % storage % mu_art(1,i,j,k) / ( dimensionless % Mach**2 * thermodynamics % gammaMinus1 * 2.0_RP ) 
            end do                ; end do                ; end do
!
!           Compute the average in each cell
!           --------------------------------
            muAver = 0.0_RP
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               muAver(:) = muAver(:) + e % storage % mu_art(:,i,j,k) * e % geom % Jacobian(i,j,k)
            end do                ; end do                ; end do
            muAver = muAver / e % geom % volume

            muAver = min(0.1_RP, muAver)

            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               e % storage % mu_art(:,i,j,k) = muAver
            end do                ; end do                ; end do
!
!           Assign to each of the neighbouring faces
!           ----------------------------------------
            do fID = 1, 6
               do j = 0, mesh % faces(e % faceIDs(fID)) % Nf(2) 
                  do i = 0, mesh % faces(e % faceIDs(fID)) % Nf(1) 
                     mesh % faces(e % faceIDs(fID)) % storage(e % faceSide(fID)) % mu_art(:,i,j) = muAver
                  end do
               end do
            end do

            end associate
         end do
!$omp end do
!
!        Average each face values
!        ------------------------
!$omp do private(muAver,i,j)
         do fID = 1, size(mesh % faces)
            if ( .not. (mesh % faces(fID) % faceType .eq. HMESH_INTERIOR) ) cycle
            muAver = max(mesh % faces(fID) % storage(1) % mu_art(:,0,0), mesh % faces(fID) % storage(2) % mu_art(:,0,0))
            do j = 0, mesh % faces(fID) % Nf(2) 
               do i = 0, mesh % faces(fID) % Nf(1) 
                  mesh % faces(fID) % storage(1) % mu_art(:,i,j) = muAver
                  mesh % faces(fID) % storage(2) % mu_art(:,i,j) = muAver
               end do
            end do
         end do
!$omp end do

!!$omp do private(h,i,j,k,fID,muAver,sBeta)
!         do fID = 1, size(mesh % faces)
!            associate(f => mesh % faces(fID))
!            h = (f % geom % surface/product(f % Nf+1))**(1.0_RP/2.0_RP) 
!!
!!           Compute Pointwise viscosity
!!           ---------------------------
!            do j = 0, f % Nf(2) ; do i = 0, f % Nf(1)
!               call ComputeShockSensor(f % storage(1) % Q(:,i,j), f % storage(1) % U_x(:,i,j), &
!                                       f % storage(1) % U_y(:,i,j) , f % storage(1) % U_z(:,i,j), h, sBeta)
!
!               f % storage(1) % mu_art(1,i,j) = f % storage(1) % Q(IRHO,i,j) * h * kBeta * sBeta
!               f % storage(1) % mu_art(3,i,j) = f % storage(1) % mu_art(1,i,j) / ( dimensionless % Mach**2 * thermodynamics % gammaMinus1 * 0.75_RP )
!            end do                ; end do
!
!            if (f % faceType .eq. HMESH_INTERIOR) then
!               do j = 0, f % Nf(2) ; do i = 0, f % Nf(1)
!                  call ComputeShockSensor(f % storage(2) % Q(:,i,j), f % storage(2) % U_x(:,i,j), &
!                                          f % storage(2) % U_y(:,i,j) , f % storage(2) % U_z(:,i,j), h, sBeta)
!
!                  f % storage(2) % mu_art(1,i,j) = f % storage(2) % Q(IRHO,i,j) * h * kBeta * sBeta
!                  f % storage(2) % mu_art(3,i,j) = f % storage(2) % mu_art(1,i,j) / ( dimensionless % Mach**2 * thermodynamics % gammaMinus1 * 0.75_RP )
!               end do                ; end do
!            else
!               f % storage(2) % mu_art = f % storage(1) % mu_art
!            end if
!
!!           Compute the average in each cell
!!           --------------------------------
!            muAver = 0.0_RP
!            do j = 0, f % Nf(2) ; do i = 0, f % Nf(1)
!               muAver(:) = muAver(:) + f % storage(1) % mu_art(:,i,j) * f % geom % jacobian(i,j)
!               muAver(:) = muAver(:) + f % storage(2) % mu_art(:,i,j) * f % geom % jacobian(i,j)
!            end do                ; end do
!            muAver = muAver / (2.0_RP * f % geom % surface)
!
!            do j = 0, f % Nf(2) ; do i = 0, f % Nf(1)
!               f % storage(1) % mu_art(:,i,j) = muAver
!               f % storage(2) % mu_art(:,i,j) = muAver
!            end do                ; end do
!            end associate
!         end do
!!$omp end do


!
!        Gather the face values onto the closest nodal value
!        ---------------------------------------------------
!!$omp do private(i,j,k)
!         do eID = 1, mesh % no_of_elements
!            associate(e => mesh % elements(eID))
!            do k = 1, e % Nxyz(3)-1 ; do j = 1, e % Nxyz(2)-1
!               e % storage % mu_art(:,0,j,k) = mesh % faces(e % faceIDs(6)) % storage(e % faceSide(6)) % mu_art(:,0,0)
!               e % storage % mu_art(:,e % Nxyz(1),j,k) = mesh % faces(e % faceIDs(4)) % storage(e % faceSide(4)) % mu_art(:,0,0)
!            end do                ; end do
!
!            do k = 1, e % Nxyz(3)-1 ; do i = 1, e % Nxyz(1)-1
!               e % storage % mu_art(:,i,0,k) = mesh % faces(e % faceIDs(1)) % storage(e % faceSide(1)) % mu_art(:,0,0)
!               e % storage % mu_art(:,i,e % Nxyz(2),k) = mesh % faces(e % faceIDs(2)) % storage(e % faceSide(2)) % mu_art(:,0,0)
!            end do                ; end do
!
!            do j = 1, e % Nxyz(2)-1 ; do i = 1, e % Nxyz(1)-1
!               e % storage % mu_art(:,i,j,0) = mesh % faces(e % faceIDs(3)) % storage(e % faceSide(3)) % mu_art(:,0,0)
!               e % storage % mu_art(:,i,j,e % Nxyz(3)) = mesh % faces(e % faceIDs(5)) % storage(e % faceSide(5)) % mu_art(:,0,0)
!            end do                ; end do
!            end associate
!         end do
!!$omp end do

      end subroutine ComputeArtificialViscosity
!
!////////////////////////////////////////////////////////////////////////////////////////
!
end module SpatialDiscretization
