#include "Includes.h"
module SpatialDiscretization
      use SMConstants
      use HyperbolicDiscretizations
      use EllipticDiscretizations
      use LESModels
      use ShockCapturing
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
      use VariableConversion, only: NSGradientVariables_STATE, GetNSViscosity, NSGradientVariables_ENTROPY, &
                                    GetGradientValues_f, NSGradientVariables_ENERGY, get_laminar_mu_kappa, &
                                    set_getVelocityGradients
      use ProblemFileFunctions, only: UserDefinedSourceTermNS_f
      use BoundaryConditions
      use IBMClass
#ifdef _HAS_MPI_
      use mpi
#endif

      private
      public   ComputeTimeDerivative, ComputeTimeDerivativeIsolated, viscousDiscretizationKey
      public   Initialize_SpaceAndTimeMethods, Finalize_SpaceAndTimeMethods

      procedure(GetGradientValues_f), pointer :: GetGradients
      procedure(EllipticFlux_f),      pointer :: ViscousFlux

      character(len=LINE_LENGTH), parameter  :: viscousDiscretizationKey = "viscous discretization"
!
!     ========
      CONTAINS
!     ========
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Initialize_SpaceAndTimeMethods(controlVariables, sem)
         use FTValueDictionaryClass
         use Utilities, only: toLower
         use mainKeywordsModule
         use Headers
         use MPI_Process_Info
         use WallFunctionConnectivity
         implicit none
         class(FTValueDictionary),  intent(in)  :: controlVariables
         class(DGSem)                           :: sem
!
!        ---------------
!        Local variables
!        ---------------
!
         character(len=LINE_LENGTH)       :: inviscidDiscretizationName
         character(len=LINE_LENGTH)       :: viscousDiscretizationName
         character(len=*), parameter      :: gradient_variables_key = "gradient variables"
         character(len=LINE_LENGTH)       :: gradient_variables
         real(RP)                         :: hnmin, hnmax

         if (.not. sem % mesh % child) then ! If this is a child mesh, all these constructs were already initialized for the parent mesh

            call hnRange(sem % mesh, hnmin, hnmax)

            if ( MPI_Process % isRoot ) then
               write(STD_OUT,'(/)')
               call Section_Header("Spatial discretization scheme")
               write(STD_OUT,'(/)')

               write(STD_OUT,'(30X,A,A30,1pG10.3)') "->", "Minimum h/N: ", hnmin
               write(STD_OUT,'(30X,A,A30,1pG10.3)') "->", "Maximum h/N: ", hnmax
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
               error stop

            end select

            call HyperbolicDiscretization % Initialize(controlVariables)
   !
   !        Initialize viscous discretization
   !        ---------------------------------
            if ( flowIsNavierStokes ) then
               if ( controlVariables % ContainsKey(gradient_variables_key) ) then
                  gradient_variables = controlVariables % StringValueForKey(gradient_variables_key, LINE_LENGTH)
                  call toLower(gradient_variables)

                  select case (trim(gradient_variables))
                  case ("state")
                     call SetGradientVariables(GRADVARS_STATE)
                     GetGradients => NSGradientVariables_STATE
                     ViscousFlux  => ViscousFlux_STATE
                     call set_getVelocityGradients(GRADVARS_STATE)

                  case ("entropy")
                     call SetGradientVariables(GRADVARS_ENTROPY)
                     GetGradients => NSGradientVariables_ENTROPY
                     ViscousFlux  => ViscousFlux_ENTROPY
                     call set_getVelocityGradients(GRADVARS_ENTROPY)

                  case ("energy")
                     call SetGradientVariables(GRADVARS_ENERGY)
                     GetGradients => NSGradientVariables_ENERGY
                     ViscousFlux  => ViscousFlux_ENERGY
                     call set_getVelocityGradients(GRADVARS_ENERGY)

                  case default
                     print*, 'Entropy variables "',trim(gradient_variables),'" are not currently implemented'
                     write(STD_OUT,'(A)') "Implemented options are:"
                     write(STD_OUT,'(A)') "  * State"
                     write(STD_OUT,'(A)') "  * Entropy"
                     write(STD_OUT,'(A)') "  * Energy"
                     errorMessage(STD_OUT)
                     error stop
                  end select

               else
!
!                 Set state as default option
!                 ---------------------------
                  call SetGradientVariables(GRADVARS_STATE)
                  GetGradients => NSGradientVariables_STATE
                  ViscousFlux  => ViscousFlux_STATE
                  call set_getVelocityGradients(GRADVARS_STATE)

               end if

               if ( .not. controlVariables % ContainsKey(viscousDiscretizationKey) ) then
                  print*, "Input file is missing entry for keyword: viscous discretization"
                  errorMessage(STD_OUT)
                  error stop
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
                  error stop

               end select

               call ViscousDiscretization % Construct(controlVariables, ELLIPTIC_NS)
               call ViscousDiscretization % Describe

            else
               if (.not. allocated(ViscousDiscretization)) allocate(EllipticDiscretization_t :: ViscousDiscretization)
               call ViscousDiscretization % Construct(controlVariables, ELLIPTIC_NS)
!
!              Set state as default option
!              ---------------------------
               call SetGradientVariables(GRADVARS_STATE)
               GetGradients => NSGradientVariables_STATE
               ViscousFlux  => ViscousFlux_STATE
               call set_getVelocityGradients(GRADVARS_STATE)

            end if

!
!           Initialize models
!           -----------------
            call InitializeLESModel(LESModel, controlVariables)
!
!           Initialize Shock-Capturing
!           --------------------------
            call Initialize_ShockCapturing(ShockCapturingDriver, controlVariables, sem, &
                                           ComputeTimeDerivative, ComputeTimeDerivativeIsolated)
            call ShockCapturingDriver % Describe

         end if

      end subroutine Initialize_SpaceAndTimeMethods
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine Finalize_SpaceAndTimeMethods
         implicit none


         if ( allocated(HyperbolicDiscretization) ) deallocate( HyperbolicDiscretization )
         if ( allocated(LESModel) )                 deallocate( LESModel )
         if ( allocated(ShockCapturingDriver) )     deallocate( ShockCapturingDriver )


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
         type(HexMesh), target           :: mesh 
         type(Particles_t)               :: particles
         REAL(KIND=RP)                   :: time
         integer, intent(in)             :: mode
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER :: k, nZones, zoneID, eID

         call SetBoundaryConditionsEqn(NS_BC)

         print*, "I am in Spatial Discretization line 247"
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
!$omp parallel shared(mesh, time)
         call HexMesh_ProlongSolToFaces(mesh, NCONS)

         print*, "I am in Spatial Discretization line 255"

!        ----------------
!        Update MPI Faces
!        ----------------
!
#ifdef _HAS_MPI_
!$omp single
         call mesh % UpdateMPIFacesSolution(NCONS)
!$omp end single
#endif

!        ------------------------------------------
!        Apply the Boundary conditions to the state
!        ------------------------------------------
!        This was done in the compute boundary flux before 
!        but it was called twice because we call it once in this file
!        and one in the Elliptic discretisation. So now we compute it
!        only once at the begining of time derivative and store it
! 
         print*, "I am in Spatial Discretization line 277"

         nZones = size(mesh % zones)
         do zoneID=1, nZones
            CALL BCs(zoneID) % bc % FlowState(mesh, zoneID)  
         enddo
!  
!        -----------------
!        Compute gradients
!        -----------------
!
         print*, "I am in Spatial Discretization line 289"
   
!$omp do schedule(runtime)
         !$acc parallel loop gang num_gangs(size(mesh % faces_interior)) vector_length(64) present(mesh, mesh % elements) async
         do eID = 1 , size(mesh % elements)
            call HexElement_ComputeLocalGradient(mesh % elements(eID))
         end do
         !$acc end parallel loop
!$omp end do nowait

         !$acc wait

         if ( computeGradients ) then
            call ViscousDiscretization % ComputeGradient( NCONS, NGRAD, mesh, time, GetGradients)
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
         print*, "I am in Spatial Discretization line 307"

         call TimeDerivative_ComputeQDot(mesh = mesh , &
                                         particles = particles, &
                                         t    = time)
         
         print*, "I am in Spatial Discretization line 313"

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
            CALL BaseClass_ComputeGradient( ViscousDiscretization, NCONS, NGRAD, mesh , time , GetGradients)
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
      use WallFunctionConnectivity
         use TripForceClass, only: randomTrip
         use ActuatorLine, only: farm
         use SpongeClass, only: sponge
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
         integer     :: eID , i, j, k, ierr, fID, iFace, iEl, iP, STLNum, n 
         real(kind=RP)  :: mu_smag, delta, Source(NCONS), TurbulentSource(NCONS), Q_target(NCONS)
         real(kind=RP), allocatable :: Source_HO(:,:,:,:)
         integer,       allocatable :: i_(:), j_(:), k_(:)
!
!        ***********************************************
!        Compute the viscosity at the elements and faces
!        ***********************************************
!
         if (flowIsNavierStokes) then
!$omp do schedule(runtime) private(i,j,k)
            !$acc parallel loop gang present(mesh)
            do eID = 1, size(mesh % elements)
               !$acc loop vector collapse(3)
               do k = 0, mesh % elements(eID) % Nxyz(3) ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1)
                  call get_laminar_mu_kappa(mesh % elements(eID) % storage % Q(:,i,j,k), &
                                            mesh % elements(eID) % storage % mu_NS(1,i,j,k),& 
                                            mesh % elements(eID) % storage % mu_NS(2,i,j,k))
               end do                ; end do                ; end do
            end do
            !$acc end parallel loop
!$omp end do
         end if

         print*, "I am in Spatial Discretization line 421"


         if ( LESModel % active) then
!$omp do schedule(runtime) private(i,j,k,delta,mu_smag)
            !$acc parallel loop gang present(mesh, LESModel)
            do eID = 1, size(mesh % elements)
               delta = (mesh % elements(eID) % geom % Volume / product(mesh % elements(eID) % Nxyz + 1)) ** (1.0_RP / 3.0_RP)
               !$acc loop vector collapse(3)
               do k = 0, mesh % elements(eID) % Nxyz(3) ; do j = 0, mesh % elements(eID) % Nxyz(2) ; do i = 0, mesh % elements(eID) % Nxyz(1)
                  call LESModel_Selector(LESModel, delta, mesh % elements(eID) % geom % dWall(i,j,k), &
                                                          mesh % elements(eID) % storage % Q(:,i,j,k),   &
                                                          mesh % elements(eID) % storage % U_x(:,i,j,k), &
                                                          mesh % elements(eID) % storage % U_y(:,i,j,k), &
                                                          mesh % elements(eID) % storage % U_z(:,i,j,k), &
                                                          mesh % elements(eID) % storage % mu_turb_NS(i,j,k) )
                                                                                   ! mu_smag)
                  ! e % storage % mu_NS(1,i,j,k) = e % storage % mu_NS(1,i,j,k) + mu_smag
                  ! e % storage % mu_NS(2,i,j,k) = e % storage % mu_NS(2,i,j,k) + mu_smag * dimensionless % mu_to_kappa
                  mesh % elements(eID) % storage % mu_NS(1,i,j,k) = mesh % elements(eID) % storage % mu_NS(1,i,j,k) + &
                                                                    mesh % elements(eID) % storage % mu_turb_NS(i,j,k)
                  mesh % elements(eID) % storage % mu_NS(2,i,j,k) = mesh % elements(eID) % storage % mu_NS(2,i,j,k) + &
                                                                    mesh % elements(eID) % storage % mu_turb_NS(i,j,k) * dimensionless % mu_to_kappa
               end do                ; end do                ; end do
            end do
            !$acc end parallel loop
!$omp end do
      end if
!
!        Compute viscosity at interior and boundary faces
!        ------------------------------------------------
         call compute_viscosity_at_faces(size(mesh % faces_interior), 2, mesh % faces_interior, mesh)
         call compute_viscosity_at_faces(size(mesh % faces_boundary), 1, mesh % faces_boundary, mesh)
         print*, "I am in Spatial Discretization line 450"

!
!        ****************
!        Volume integrals
!        ****************
!
!$omp do schedule(runtime)
         !$acc parallel loop gang present(mesh, mesh % elements, mesh % elements % geom)
         do eID = 1 , size(mesh % elements)
            call TimeDerivative_VolumetricContribution(mesh % elements(eID))         
         end do
         !$acc end parallel loop 
!$omp end do

         print*, "I am in Spatial Discretization line 465"

#if defined(_HAS_MPI_)
!$omp single
         if (ShockCapturingDriver % isActive) then
            call mesh % UpdateMPIFacesAviscflux(NCONS)
         end if
!$omp end single
#endif
!
!        ******************************************
!        Compute Riemann solver of non-shared faces
!        ******************************************
!
      print*, "I am in Spatial Discretization line 479"

!$omp do schedule(runtime) private(fID)
         !$acc parallel loop gang num_gangs(size(mesh % faces_interior)) present(mesh, mesh % faces, mesh % faces_interior) private(fID)
         do iFace = 1, size(mesh % faces_interior)
            fID = mesh % faces_interior(iFace)
            call computeElementInterfaceFlux(mesh % faces(fID))
         end do
         !$acc end parallel loop
!$omp end do nowait

         print*, "I am in Spatial Discretization line 490"

         call computeBoundaryFlux(mesh, t)

         print*, "I am in Spatial Discretization line 494"

!
!        ***************************************************************
!        Surface integrals and scaling of elements with non-shared faces
!        ***************************************************************
!
!$omp do schedule(runtime) private(i,j,k,eID)
         !$acc parallel loop gang present(mesh, mesh % elements, mesh % elements_sequential) copyin(t)
         do iEl = 1, size(mesh % elements_sequential)
            eID = mesh % elements_sequential(iEl)
            call TimeDerivative_FacesContribution(mesh % elements(eID), t, mesh)
         end do
         !$acc end parallel loop 
!$omp end do
!
         print*, "I am in Spatial Discretization line 510"

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
!           Compute viscosity at MPI faces
!           ------------------------------
            call compute_viscosity_at_faces(size(mesh % faces_mpi), 2, mesh % faces_mpi, mesh)

!$omp single
            if ( flowIsNavierStokes ) then
               if ( ShockCapturingDriver % isActive ) then
                  call mpi_barrier(MPI_COMM_WORLD, ierr)     ! TODO: This can't be the best way :(
                  call mesh % GatherMPIFacesAviscflux(NCONS)
               end if
            end if
!$omp end single
!
!           **************************************
!           Compute Riemann solver of shared faces
!           **************************************
!
!$omp do schedule(runtime) private(fID)
            do iFace = 1, size(mesh % faces_mpi)
               fID = mesh % faces_mpi(iFace)
               call computeMPIFaceFlux(mesh % faces(fID))
            end do
!$omp end do
!
!           ***********************************************************
!           Surface integrals and scaling of elements with shared faces
!           ***********************************************************
!
!$omp do schedule(runtime) private(i,j,k,eID)
            do iEl = 1, size(mesh % elements_mpi)
               eID = mesh % elements_mpi(iEl)
               associate(e => mesh % elements(eID))
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
               ! the source term is reset to 0 each time Qdot is calculated to enable the possibility to add source terms to
               ! different contributions and not accumulate each call
               e % storage % S_NS = 0.0_RP
               do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  call UserDefinedSourceTermNS(e % geom % x(:,i,j,k), e % storage % Q(:,i,j,k), t, e % storage % S_NS(:,i,j,k), thermodynamics, dimensionless, refValues)
                  call randomTrip % getTripSource( e % geom % x(:,i,j,k), e % storage % S_NS(:,i,j,k) )
                  call farm % ForcesFarm(e % geom % x(:,i,j,k), e % storage % Q(:,i,j,k), e % storage % S_NS(:,i,j,k), t)
               end do                  ; end do                ; end do
               end associate
            end do
!$omp end do
            ! for the sponge, loops are in the internal subroutine as values are precalculated
            call sponge % addSource(mesh)
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
!
!        *********************
!        Add IBM source term
!        *********************
         if( mesh% IBM% active ) then
            if( .not. mesh% IBM% semiImplicit ) then 
!$omp do schedule(runtime) private(i,j,k,Source,Q_target)
                  do eID = 1, mesh % no_of_elements  
                     associate ( e => mesh % elements(eID) ) 
                     do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                        if( e% isInsideBody(i,j,k) ) then
                           if( mesh% IBM% stl(e% STL(i,j,k))% move ) then 
                              Q_target = mesh% IBM% MaskVelocity( e% storage% Q(:,i,j,k), NCONS, e% STL(i,j,k), e% geom% x(:,i,j,k), t )
                              call mesh% IBM% SourceTerm( eID = eID, Q = e % storage % Q(:,i,j,k), Q_target = Q_target, Source = Source, wallfunction = .false. )
                           else 
                              call mesh% IBM% SourceTerm( eID = eID, Q = e % storage % Q(:,i,j,k), Source = Source, wallfunction = .false. )
                           end if 
                           e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) + Source
                        end if
                     end do                  ; end do                ; end do
                     end associate
                  end do
!$omp end do       
               if( mesh% IBM% Wallfunction ) then
!$omp single
                  call mesh% IBM% GetBandRegionStates( mesh% elements )
!$omp end single
!$omp do schedule(runtime) private(i,j,k,TurbulentSource)
                  do iP = 1, mesh% IBM% NumOfForcingPoints
                     associate( e    => mesh% elements(mesh% IBM% ImagePoints(iP)% element_index), &
                                e_in => mesh% elements(mesh% IBM% ImagePoints(iP)% element_in)     )
                     i = mesh% IBM% ImagePoints(iP)% local_position(1)
                     j = mesh% IBM% ImagePoints(iP)% local_position(2)
                     k = mesh% IBM% ImagePoints(iP)% local_position(3)
                     call mesh % IBM % SourceTermTurbulence( mesh% IBM% ImagePoints(iP), e% storage% Q(:,i,j,k), &
                                                             e% geom% normal(:,i,j,k), e% geom% dWall(i,j,k),    &
                                                             e% STL(i,j,k), TurbulentSource                      )             
                     e% storage% QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) + TurbulentSource  
                     end associate
                  end do
!$omp end do
               end if 
            end if 
         end if


      end subroutine TimeDerivative_ComputeQDot

      subroutine compute_viscosity_at_faces(no_of_faces, no_of_sides, face_ids, mesh)
         implicit none
         integer, intent(in)           :: no_of_faces
         integer, intent(in)           :: no_of_sides
         integer, intent(in)           :: face_ids(no_of_faces)
         class(HexMesh), intent(inout) :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: iFace, i, j, side, fID
         real(kind=RP) :: delta, mu_smag

         if (flowIsNavierStokes) then
!$omp do schedule(runtime) private(i,j)
            !$acc parallel loop gang present(mesh, face_ids)
            do iFace = 1, no_of_faces
               fID = face_ids(iFace)
               !$acc loop vector collapse(3)
               do j = 0, mesh % faces(fID) % Nf(2) ; do i = 0, mesh % faces(fID) % Nf(1)
                  do side = 1, no_of_sides
                      call get_laminar_mu_kappa(mesh % faces(fID) % storage(side) % Q(:,i,j), &
                                                mesh % faces(fID) % storage(side) % mu_NS(1,i,j), &
                                                mesh % faces(fID) % storage(side) % mu_NS(2,i,j))
                  end do
               end do              ; end do
            end do
            !$acc end parallel loop
!$omp end do
         end if

         if ( LESModel % Active ) then
!$omp do schedule(runtime) private(i,j,delta,mu_smag)
            !$acc parallel loop gang present(mesh, LESModel)
            do iFace = 1, no_of_faces
               delta = sqrt(mesh % faces(face_ids(iFace)) % geom % surface / product(mesh % faces(face_ids(iFace)) % Nf + 1))
               !$acc loop vector collapse(3)
               do j = 0, mesh % faces(face_ids(iFace)) % Nf(2) ; do i = 0, mesh % faces(face_ids(iFace)) % Nf(1)
                  do side = 1, no_of_sides
                     call LESModel_Selector(LESModel, delta, mesh % faces(face_ids(iFace)) % geom % dWall(i,j), &
                                                             mesh % faces(face_ids(iFace)) % storage(side) % Q(:,i,j),   &
                                                             mesh % faces(face_ids(iFace)) % storage(side) % U_x(:,i,j), &
                                                             mesh % faces(face_ids(iFace)) % storage(side) % U_y(:,i,j), &
                                                             mesh % faces(face_ids(iFace)) % storage(side) % U_z(:,i,j), &
                                                                                    mu_smag)
                     mesh % faces(face_ids(iFace)) % storage(side) % mu_NS(1,i,j) = mesh % faces(face_ids(iFace)) % storage(side) % mu_NS(1,i,j) + mu_smag
                     mesh % faces(face_ids(iFace)) % storage(side) % mu_NS(2,i,j) = mesh % faces(face_ids(iFace)) % storage(side) % mu_NS(2,i,j) + mu_smag * dimensionless % mu_to_kappa
                  end do
               end do              ; end do
            end do
            !$acc end parallel loop
!$omp end do
         end if




      end subroutine compute_viscosity_at_faces
!
!////////////////////////////////////////////////////////////////////////
!
!     -------------------------------------------------------------------------------
!     This routine computes Qdot neglecting the interaction with neighboring elements
!     and boundaries. Therefore, the external states are not needed.
!     -------------------------------------------------------------------------------
      subroutine TimeDerivative_ComputeQDotIsolated( mesh , t )
         use TripForceClass, only: randomTrip
         use ActuatorLine, only: farm
         implicit none
         type(HexMesh)              :: mesh
         real(kind=RP)              :: t
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: eID , i, j, k, fID, iP
         procedure(UserDefinedSourceTermNS_f) :: UserDefinedSourceTermNS
         real(kind=rp) :: Source(NCONS), TurbulentSource(NCONS)
!
!        ****************
!        Volume integrals
!        ****************
!
!$omp do schedule(runtime)
         do eID = 1 , size(mesh % elements)
            call TimeDerivative_StrongVolumetricContribution(mesh, mesh % elements(eID) , t)
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
               e % storage % S_NS = 0.0_RP
               do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  call UserDefinedSourceTermNS(e % geom % x(:,i,j,k), e % storage % Q(:,i,j,k), t, e % storage % S_NS(:,i,j,k), thermodynamics, dimensionless, refValues)
                  call randomTrip % getTripSource( e % geom % x(:,i,j,k), e % storage % S_NS(:,i,j,k) )
                  call farm % ForcesFarm(e % geom % x(:,i,j,k), e % storage % Q(:,i,j,k), e % storage % S_NS(:,i,j,k), t)
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

!
!        *********************
!        Add IBM source term
!        *********************

         if( mesh% IBM% active ) then
            if( .not. mesh% IBM% semiImplicit ) then 
!$omp do schedule(runtime) private(i,j,k,Source)
               do eID = 1, mesh % no_of_elements
                  associate ( e => mesh % elements(eID) )
                  do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                     if( e% isInsideBody(i,j,k) ) then
                        call mesh% IBM% SourceTerm( eID = eID, Q = e % storage % Q(:,i,j,k), Source = Source, wallfunction = .false. )
                        e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) + Source
                     end if
                  end do                  ; end do                ; end do
                  end associate
               end do
!$omp end do      
               if( mesh% IBM% Wallfunction ) then
!$omp single
                  call mesh% IBM% GetBandRegionStates( mesh% elements )
!$omp end single 
!$omp do schedule(runtime) private(i,j,k,TurbulentSource)
                  do iP = 1, mesh% IBM% NumOfForcingPoints
                     associate( e    => mesh% elements(mesh% IBM% ImagePoints(iP)% element_index) )
                     i = mesh% IBM% ImagePoints(iP)% local_position(1)
                     j = mesh% IBM% ImagePoints(iP)% local_position(2)
                     k = mesh% IBM% ImagePoints(iP)% local_position(3)
                     call mesh % IBM % SourceTermTurbulence( mesh% IBM% ImagePoints(iP), e% storage% Q(:,i,j,k), &
                                                             e% geom% normal(:,i,j,k), e% geom% dWall(i,j,k),    &
                                                             e% STL(i,j,k), TurbulentSource                      )              
                     e% storage% QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) + TurbulentSource  
                     end associate
                  end do
!$omp end do
               end if 
            end if
         endif

      end subroutine TimeDerivative_ComputeQDotIsolated
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine TimeDerivative_StrongVolumetricContribution(mesh, e, t)
         use HexMeshClass
         use ElementClass
         implicit none
         type(HexMesh)             :: mesh
         type(Element)             :: e
         real(kind=RP), intent(in) :: t

!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: inviscidContravariantFlux ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM )
         real(kind=RP) :: fSharp(1:NCONS, 0:e%Nxyz(1), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP) :: gSharp(1:NCONS, 0:e%Nxyz(2), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP) :: hSharp(1:NCONS, 0:e%Nxyz(3), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP) :: viscousContravariantFlux ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM )
         real(kind=RP) :: AviscContravariantFlux   ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM )
         real(kind=RP) :: contravariantFlux        ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM )

!
!        *************************************
!        Compute interior contravariant fluxes
!        *************************************
!
!        Compute inviscid contravariant flux
!        -----------------------------------
         call HyperbolicDiscretization % ComputeInnerFluxes ( e , EulerFlux, inviscidContravariantFlux )
!
!        Compute viscous contravariant flux
!        ----------------------------------
         if (flowIsNavierStokes) then

            call ViscousDiscretization  % ComputeInnerFluxes ( NCONS, NGRAD, ViscousFlux, GetNSViscosity, e, viscousContravariantFlux)
!
!           Compute the artificial dissipation
!           ----------------------------------
            if (ShockCapturingDriver % isActive) then
               call ShockCapturingDriver % ComputeViscosity(mesh, e, AviscContravariantFlux)
            else
               AviscContravariantFlux = 0.0_RP
            end if

         else

            viscousContravariantFlux = 0.0_RP
            AviscContravariantFlux   = 0.0_RP

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
            contravariantFlux = inviscidContravariantFlux - viscousContravariantFlux - AviscContravariantFlux
!
!           Perform the Weak Volume Green integral
!           --------------------------------------
            e % storage % QDot = ScalarStrongIntegrals % StdVolumeGreen ( e , NCONS, contravariantFlux )

         type is (SplitDG_t)
!
!           Compute sharp fluxes for skew-symmetric approximations
!           ------------------------------------------------------
            call HyperbolicDiscretization % ComputeSplitFormFluxes(e, inviscidContravariantFlux, fSharp, gSharp, hSharp)
!
!           Perform the Weak volume green integral
!           --------------------------------------
            viscousContravariantFlux = viscousContravariantFlux + AviscContravariantFlux

            e % storage % QDot = -ScalarStrongIntegrals % SplitVolumeDivergence( e, fSharp, gSharp, hSharp, viscousContravariantFlux)

         end select

      end subroutine TimeDerivative_StrongVolumetricContribution
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine TimeDerivative_VolumetricContribution(element)
         !$acc routine vector
         use HexMeshClass
         use ElementClass
         use DGIntegrals
         implicit none
         type(Element), intent (inout)           :: element
!
!        ---------------
!        Local variables
!        ---------------
!
         integer            :: i, j, k, l, eID
         real(kind=RP)      :: mu, kappa, beta

         real(kind=RP) :: inviscidFlux(1:NCONS, 1:NDIM)
         real(kind=RP) :: viscousFlux(1:NCONS, 1:NDIM)

!
!        *************************************
!        Compute interior contravariant fluxes
!        *************************************
!
!        Compute inviscid contravariant flux
!        -----------------------------------
         
         !$acc loop vector collapse(3) private(inviscidFlux, viscousFlux)
         do k = 0, element % Nxyz(3) ; do j = 0, element % Nxyz(2) ; do i = 0, element % Nxyz(1)
                  
            call EulerFlux(element % storage % Q(:,i,j,k), inviscidFlux, element % storage % rho(i,j,k))

            mu = element % storage % mu_ns(1,i,j,k)
            kappa = element % storage % mu_ns(2,i,j,k)
            beta  = 0.0_RP

            call ViscousFlux_STATE( NCONS, NGRAD, element % storage % Q(:,i,j,k) , element % storage % U_x(:,i,j,k) , & 
                                    element % storage % U_y(:,i,j,k) , element % storage % U_z(:,i,j,k), mu, beta, kappa, viscousFlux)
                  
            inviscidFlux = inviscidFlux - viscousFlux
                  
            element % storage % FluxF(:,i,j,k) = inviscidFlux(:,IX) * element % geom % jGradXi(IX,i,j,k)  &
                                               + inviscidFlux(:,IY) * element % geom % jGradXi(IY,i,j,k)  &
                                               + inviscidFlux(:,IZ) * element % geom % jGradXi(IZ,i,j,k)

            element % storage % FluxG(:,i,j,k) = inviscidFlux(:,IX) * element % geom % jGradEta(IX,i,j,k)  &
                                               + inviscidFlux(:,IY) * element % geom % jGradEta(IY,i,j,k)  &
                                               + inviscidFlux(:,IZ) * element % geom % jGradEta(IZ,i,j,k)
                  
            element % storage % FluxH(:,i,j,k) = inviscidFlux(:,IX) * element % geom % jGradZeta(IX,i,j,k)  &
                                               + inviscidFlux(:,IY) * element % geom % jGradZeta(IY,i,j,k)  &
                                               + inviscidFlux(:,IZ) * element % geom % jGradZeta(IZ,i,j,k)
         end do               ; end do                ; end do

         call ScalarWeakIntegrals_StdVolumeGreen( element % Nxyz, NCONS,&
                                                  element % storage % FluxF, &
                                                  element % storage % FluxG, & 
                                                  element % storage % FluxH, & 
                                                  element % storage % QDot)

      end subroutine TimeDerivative_VolumetricContribution

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine TimeDerivative_FacesContribution( element , t , mesh)
         !$acc routine vector
         use HexMeshClass
         implicit none
         type(Element)           :: element
         real(kind=RP)           :: t
         type(HexMesh)           :: mesh

         integer                 :: i,j,k,eID,eq

         call ScalarWeakIntegrals_StdFace( NCONS, element % Nxyz, &
                      mesh % faces(element % faceIDs(EFRONT))  % storage(element % faceSide(EFRONT))  % fStar, &
                      mesh % faces(element % faceIDs(EBACK))   % storage(element % faceSide(EBACK))   % fStar, &
                      mesh % faces(element % faceIDs(EBOTTOM)) % storage(element % faceSide(EBOTTOM)) % fStar, &
                      mesh % faces(element % faceIDs(ERIGHT))  % storage(element % faceSide(ERIGHT))  % fStar, &
                      mesh % faces(element % faceIDs(ETOP))    % storage(element % faceSide(ETOP))    % fStar, &
                      mesh % faces(element % faceIDs(ELEFT))   % storage(element % faceSide(ELEFT))   % fStar, &
                      element % storage % QDot )

         !$acc loop vector collapse(3)
         do k = 0, element % Nxyz(3) ; do j = 0, element % Nxyz(2) ; do i = 0, element % Nxyz(1)
            !$acc loop seq
            do eq = 1,NCONS
            element % storage % QDot(eq,i,j,k) = element % storage % QDot(eq,i,j,k)  / element % geom % jacobian(i,j,k)
            enddo
         end do         ; end do          ; end do
         
      end subroutine TimeDerivative_FacesContribution
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
!        Riemann solver drivers
!        ----------------------
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine computeElementInterfaceFlux(fc)
         !$acc routine vector
         use FaceClass
         use RiemannSolvers_NS
         use EllipticBR1
         implicit none
         type(Face)   , intent(inout) :: fc
         
         integer       :: i, j, eq

         call ViscousFlux_selector(NCONS, NGRAD, fc % Nf(1), &
                                   fc % Nf(2), 0, &
                                   fc % storage(1) % Q , &
                                   fc % storage(1) % U_x, &
                                   fc % storage(1) % U_y, &
                                   fc % storage(1) % U_z, &
                                   fc % storage(1) % mu_NS, &
                                   fc % storage(1) % unStar)
           
         call ViscousFlux_selector(NCONS, NGRAD, fc % Nf(1), & 
                                   fc % Nf(2), 0, &
                                   fc % storage(2) % Q , &
                                   fc % storage(2) % U_x, &
                                   fc % storage(2) % U_y, &
                                   fc % storage(2) % U_z, &
                                   fc % storage(2) % mu_NS, &
                                   fc % storage(2) % unStar)


         call BR1_RiemannSolver_acc(fc, NCONS, NGRAD, fc % storage(2) % FStar)

         call RiemannSolver_Selector(fc % Nf(1), &
                                     fc % Nf(2), &
                                     fc % storage(1) % Q, &
                                     fc % storage(2) % Q, &
                                     fc % geom % normal, &
                                     fc % geom % t1, &
                                     fc % geom % t2, &
                                     fc % storage(1) % FStar )
        
!        ------------------------
!        Multiply by the Jacobian
!        ------------------------
         !$acc loop vector collapse(2)
         do j = 0, fc % Nf(2) ; do i = 0, fc % Nf(1)
            !$acc loop seq
            do eq = 1, NCONS
               fc % storage(1) % FStar(eq,i,j) = (fc % storage(1) % FStar(eq,i,j) - fc % storage(2) % FStar(eq,i,j)) * fc % geom % jacobian(i,j)
            enddo
         end do ;  end do
!
!        ---------------------------
!        Return the flux to elements
!        ---------------------------
!         
        call Face_ProjectFluxToElements(fc, NCONS, fc % storage(1) % FStar, 1)
        call Face_ProjectFluxToElements(fc, NCONS, fc % storage(1) % FStar, 2)

      end subroutine computeElementInterfaceFlux
   
      subroutine computeMPIFaceFlux(f)
         use FaceClass
         use RiemannSolvers_NS
         implicit none
         type(Face)   , intent(inout) :: f
         integer       :: i, j
         integer       :: thisSide
         real(kind=RP) :: inv_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: Avisc_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: mu_left(3), mu_right(3)
         integer       :: Sidearray(2)
!
!        ---------------------------
!        Artificial viscosity fluxes
!        ---------------------------
!
         if ( ShockCapturingDriver % isActive ) then
            Avisc_flux = 0.5_RP * (f % storage(1) % AviscFlux + f % storage(2) % AviscFlux)
         else
            Avisc_flux = 0.0_RP
         end if
!
!        --------------
!        Viscous fluxes
!        --------------
!
         if (flowIsNavierStokes) then
            do j = 0, f % Nf(2)
               do i = 0, f % Nf(1)

                  mu_left(1) = f % storage(1) % mu_NS(1,i,j)
                  mu_left(2) = 0.0_RP
                  mu_left(3) = f % storage(1) % mu_NS(2,i,j)

                  mu_right(1) = f % storage(2) % mu_NS(1,i,j)
                  mu_right(2) = 0.0_RP
                  mu_right(3) = f % storage(2) % mu_NS(2,i,j)

                  call ViscousDiscretization % RiemannSolver(nEqn = NCONS, nGradEqn = NGRAD, &
                                                     EllipticFlux = ViscousFlux, &
                                                     f = f, &
                                                     QLeft = f % storage(1) % Q(:,i,j), &
                                                     QRight = f % storage(2) % Q(:,i,j), &
                                                     U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                     U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                     U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                     U_xRight = f % storage(2) % U_x(:,i,j), &
                                                     U_yRight = f % storage(2) % U_y(:,i,j), &
                                                     U_zRight = f % storage(2) % U_z(:,i,j), &
                                                     mu_left  = mu_left, &
                                                     mu_right = mu_right, &
                                                     nHat = f % geom % normal(:,i,j) , &
                                                     dWall = f % geom % dWall(i,j), &
                                                     flux  = visc_flux(:,i,j) )

               end do
            end do
         else
            visc_flux = 0.0_RP
         end if

         do j = 0, f % Nf(2)
            do i = 0, f % Nf(1)
!
!              --------------
!              Invscid fluxes
!              --------------
!
               !call RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
               !                   QRight = f % storage(2) % Q(:,i,j), &
               !                   nHat   = f % geom % normal(:,i,j), &
               !                   t1     = f % geom % t1(:,i,j), &
               !                   t2     = f % geom % t2(:,i,j), &
               !                   flux   = inv_flux(:,i,j) )
!
!              Multiply by the Jacobian
!              ------------------------
               flux(:,i,j) = ( inv_flux(:,i,j) - visc_flux(:,i,j)) * f % geom % jacobian(i,j) - Avisc_flux(:,i,j)

            end do
         end do
!
!        ---------------------------
!        Return the flux to elements: The sign in eR % storage % FstarB has already been accouted.
!        ---------------------------
!
         thisSide = maxloc(f % elementIDs, dim = 1)

         Sidearray = (/thisSide, HMESH_NONE/)
         !call f % ProjectFluxToElements(NCONS, flux, Sidearray )

      end subroutine ComputeMPIFaceFlux

      SUBROUTINE computeBoundaryFlux(mesh, time)
      USE ElementClass
      use FaceClass
      USE RiemannSolvers_NS
      use WallFunctionBC
      use WallFunctionConnectivity
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP)                :: time
      type(HexMesh), intent(in)    :: mesh
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER                         :: i, j, eq
      INTEGER                         :: nZones, zoneID, zonefID, fID
      REAL(KIND=RP)                   :: inv_flux(NCONS)
      integer                         :: Sidearray(2)
!
!     -------------------
!     Get external states
!     -------------------
!
      nZones = size(mesh % zones)
      do zoneID=1, nZones

         !CALL BCs(zoneID) % bc % FlowState(mesh, zoneID)      
         
         !$acc parallel loop gang present(mesh) async(zoneID)
         do zonefID = 1, mesh % zones(zoneID) % no_of_faces
            fID =  mesh % zones(zoneID) % faces(zonefID)

            call ViscousFlux_selector( NCONS, NGRAD, mesh % faces(fID) % Nf(1), &
                                         mesh % faces(fID) % Nf(2), 0, &
                                         mesh % faces(fID) % storage(1) % Q , &
                                         mesh % faces(fID) % storage(1) % U_x, &
                                         mesh % faces(fID) % storage(1) % U_y, &
                                         mesh % faces(fID) % storage(1) % U_z, &
                                         mesh % faces(fID) % storage(1) % mu_NS, &
                                         mesh % faces(fID) % storage(1) % unStar)
            
            !$acc loop vector collapse(2)
            do j = 0, mesh % faces(fID) % Nf(2) ;  do i = 0, mesh % faces(fID) % Nf(1)
               !$acc loop seq
               do eq = 1, NCONS
                  mesh % faces(fID) % storage(2) % FStar(eq,i,j) = mesh % faces(fID) % storage(1) % unStar(eq,IX,i,j)* mesh % faces(fID) % geom % normal(IX,i,j) &
                                                                 + mesh % faces(fID) % storage(1) % unStar(eq,IY,i,j)* mesh % faces(fID) % geom % normal(IY,i,j) &
                                                                 + mesh % faces(fID) % storage(1) % unStar(eq,IZ,i,j)* mesh % faces(fID) % geom % normal(IZ,i,j)
               enddo
            enddo ; enddo
         enddo
         !$acc end parallel loop 

         !$acc wait

         CALL BCs(zoneID) % bc % FlowNeumann(mesh, zoneID)                             
         
         !$acc wait

         !$acc parallel loop gang present(mesh) private(inv_flux) async(zoneID)
         do zonefID = 1, mesh % zones(zoneID) % no_of_faces
            fID =  mesh % zones(zoneID) % faces(zonefID)

            call RiemannSolver_Selector(Nx = mesh % faces(fID) % Nf(1), &
                                        Ny = mesh % faces(fID) % Nf(2), &
                                        QLeft  = mesh % faces(fID) % storage(1) % Q, &
                                        QRight = mesh % faces(fID) % storage(2) % Q, &
                                        nHat   = mesh % faces(fID) % geom % normal, &
                                        t1     = mesh % faces(fID) % geom % t1, &
                                        t2     = mesh % faces(fID) % geom % t2, &
                                        flux   = mesh % faces(fID) % storage(1) % FStar )
!           ------------------------
!           Multiply by the Jacobian
!           ------------------------
            !$acc loop vector collapse(2)
            do j = 0, mesh % faces(fID) % Nf(2) ; do i = 0, mesh % faces(fID) % Nf(1)
               !$acc loop seq
               do eq = 1, NCONS
                  mesh % faces(fID) % storage(1) % FStar(eq,i,j) = (mesh % faces(fID) % storage(1) % FStar(eq,i,j)  - &
                                                                    mesh % faces(fID) % storage(2) % FStar(eq,i,j)) * &
                                                                    mesh % faces(fID) % geom % jacobian(i,j)
               enddo
            end do ;  end do
            
!
!           ---------------------------
!           Return the flux to elements
!           ---------------------------
!
            call Face_ProjectFluxToElements(mesh % faces(fID), NCONS, mesh % faces(fID) % storage(1) % FStar, 1)
            
         enddo
         !$acc end parallel loop 
      enddo

      !$acc wait
      end subroutine computeBoundaryFlux

      subroutine ViscousFlux_selector(nEqn, nGradEqn, Nx, Ny, Nz, Q, U_x, U_y, U_z, mu, flux_cart)
         !$acc routine vector
         use SMConstants
         use HexMeshClass
         implicit none

         integer,       intent(in)       :: nEqn
         integer,       intent(in)       :: nGradEqn
         integer,       intent(in)       :: Nx
         integer,       intent(in)       :: Ny
         integer,       intent(in)       :: Nz
         real(kind=RP), intent(in)       :: Q(1:nEqn, 0:Nx, 0:Ny, 0:Nz)
         real(kind=RP), intent(in)       :: U_x(1:nGradEqn, 0:Nx, 0:Ny, 0:Nz)
         real(kind=RP), intent(in)       :: U_y(1:nGradEqn, 0:Nx, 0:Ny, 0:Nz)
         real(kind=RP), intent(in)       :: U_z(1:nGradEqn, 0:Nx, 0:Ny, 0:Nz)
         real(kind=RP), intent(in)       :: mu(1:3, 0:Nx, 0:Ny, 0:Nz)
         real(kind=RP), intent(out)    :: flux_cart(1:nEqn, 1:NDIM, 0:Nx, 0:Ny, 0:Nz)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer :: i,j,k
         
         !select case(which_viscousflux)
         !case(VSCFlux_STATE)
            !$acc loop vector collapse(3)
            do k = 0, Nz
               do j = 0, Ny
                  do i = 0, Nx
                     call ViscousFlux_STATE(nEqn, nGradEqn, Q(:,i,j,k),  U_x(:,i,j,k), U_y(:,i,j,k), U_z(:,i,j,k), &
                                            mu(1,i,j,k), 0.0_RP, mu(3,i,j,k), flux_cart(:,:,i,j,k))
                  enddo
               enddo
            enddo
         !end select
      
      end subroutine ViscousFlux_selector

end module SpatialDiscretization
