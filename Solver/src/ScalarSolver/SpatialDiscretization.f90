#include "Includes.h"
module SpatialDiscretization
      use SMConstants
      !use HyperbolicDiscretizations
      use EllipticDiscretizations
      use DGIntegrals
      use MeshTypes
      use HexMeshClass
      use ElementClass
      use PhysicsStorage
      use Physics_SLR
      use MPI_Face_Class
      use MPI_Process_Info
      use DGSEMClass
      use ParticlesClass
      use FluidData
      use VariableConversion_SLR, only: SLRGradientVariables, GetViscosity 
      use ProblemFileFunctions
      use BoundaryConditions, only: BCs
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
         !character(len=LINE_LENGTH)       :: inviscidDiscretizationName
         character(len=LINE_LENGTH)       :: viscousDiscretizationName

         if (.not. mesh % child) then ! If this is a child mesh, all these constructs were already initialized for the parent mesh
         
            if ( MPI_Process % isRoot ) then
               write(STD_OUT,'(/)')
               call Section_Header("Spatial discretization scheme")
               write(STD_OUT,'(/)')
            end if

   !
   !        Initialize viscous discretization
   !        ---------------------------------         
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

               call ViscousDiscretization % Construct(controlVariables, ELLIPTIC_SLR)
               call ViscousDiscretization % Describe


!
!              Compute wall distances
!              ----------------------
               !call mesh % ComputeWallDistances

          end if

      end subroutine Initialize_SpaceAndTimeMethods
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine Finalize_SpaceAndTimeMethods
         implicit none
        ! IF ( ALLOCATED(HyperbolicDiscretization) ) DEALLOCATE( HyperbolicDiscretization )
      end subroutine Finalize_SpaceAndTimeMethods
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ComputeTimeDerivative( mesh, particles, time, mode, HO_Elements)
         IMPLICIT NONE 
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(HexMesh), target           :: mesh
         type(Particles_t)               :: particles
         REAL(KIND=RP)                   :: time
         integer,             intent(in) :: mode
         logical, intent(in), optional   :: HO_Elements
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER :: k, eID, fID


!$omp parallel shared(mesh, time) private(k, eID)

!
!        ------------------------------
!        Change memory to scalar
!        ------------------------------
!
!$omp do schedule(runtime)
         do eID = 1, mesh % no_of_elements
            call mesh % elements(eID) % storage % SetStorageToSLR_slr
         end do
!$omp end do

!$omp do schedule(runtime)
         do fID = 1, size(mesh % faces)
            call mesh % faces(fID) % storage(1) % SetStorageToSLR_slr
            call mesh % faces(fID) % storage(2) % SetStorageToSLR_slr
         end do
!$omp end do

! !$omp single
!          call SetBoundaryConditionsEqn(C_BC)
! !$omp end single

!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
         call mesh % ProlongSolutionToFaces(NCONS)
!
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
         call mesh % UpdateMPIFacesGradients(NCONS)
!$omp end single
#endif

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
      SUBROUTINE ComputeTimeDerivativeIsolated( mesh, particles, time, mode, HO_Elements)
         use EllipticDiscretizationClass
         IMPLICIT NONE 
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(HexMesh), target           :: mesh
         type(Particles_t)               :: particles
         REAL(KIND=RP)                   :: time
         integer,             intent(in) :: mode
         logical, intent(in), optional   :: HO_Elements
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER :: k

!
      END SUBROUTINE ComputeTimeDerivativeIsolated
!
!////////////////////////////////////////////////////////////////////////////////////
!
!           Navier--Stokes procedures
!           -------------------------
!
!////////////////////////////////////////////////////////////////////////////////////
!
      subroutine TimeDerivative_ComputeQDot( mesh , particles, t )
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
               CALL computeElementInterfaceFlux_SLR( f ) 
 
            case (HMESH_BOUNDARY) 
               CALL computeBoundaryFlux_SLR(f, t) 
 
            end select 
            end associate 
         end do 
!$omp end do 
!
!        **************************************************************
!        Surface integrals and scaling of elements without shared faces
!        **************************************************************
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
            call mesh % GatherMPIFacesGradients(NCONS)
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
                  CALL computeMPIFaceFlux_SLR( f ) 
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
!           Add a MPI Barrier
!           -----------------
!$omp single
            call mpi_barrier(MPI_COMM_WORLD, ierr)
!$omp end single
         end if
#endif
!
!           ***************
!           Add source term
!           ***************
!$omp do schedule(runtime) private(i,j,k)
            do eID = 1, mesh % no_of_elements
               associate ( e => mesh % elements(eID) )
               do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  call UserDefinedSourceTermNS(e % geom % x(:,i,j,k), e % storage % Q(:,i,j,k), t, e % storage % S_SLR(:,i,j,k), thermodynamics, dimensionless, refValues)
               end do                  ; end do                ; end do
               end associate
            end do
!$omp end do

!$omp do schedule(runtime) private(i,j,k)
            do eID = 1, mesh % no_of_elements
               associate ( e => mesh % elements(eID) )
               do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) + e % storage % S_SLR(:,i,j,k)
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

         error stop 'ComputeTimeDerivativeIsolated not implemented'
         
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
         real(kind=RP) :: viscousContravariantFlux  ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         real(kind=RP) :: contravariantFlux         ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
!
!        Compute viscous contravariant flux
         call ViscousDiscretization  % ComputeInnerFluxes ( NCONS, NCONS, slrViscousFlux, GetViscosity , e , viscousContravariantFlux) 
!        
!        ************************
!        Perform volume integrals
!        ************************

         contravariantFlux = - viscousContravariantFlux

         e % storage % QDot = ScalarWeakIntegrals % StdVolumeGreen ( e, NCONS, contravariantFlux ) 

         
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
      SUBROUTINE computeElementInterfaceFlux_SLR(f)
         use FaceClass
         use RiemannSolvers_SLR
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j
         !real(kind=RP) :: inv_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: muL, muR, mu

         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)

               call GetViscosity(f % storage(1) % Q(1,i,j), muL)
               call GetViscosity(f % storage(2) % Q(1,i,j), muR)

               mu = 0.5_RP * (muL + muR)
!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               CALL ViscousDiscretization % RiemannSolver(nEqn = NCONS, nGradEqn = NCONS, &
                                                  EllipticFlux = slrViscousFlux, &
                                                  f = f, &
                                                  QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = f % storage(1) % slr_x(:,i,j), &
                                                  U_yLeft = f % storage(1) % slr_y(:,i,j), &
                                                  U_zLeft = f % storage(1) % slr_z(:,i,j), &
                                                  U_xRight = f % storage(2) % slr_x(:,i,j), &
                                                  U_yRight = f % storage(2) % slr_y(:,i,j), &
                                                  U_zRight = f % storage(2) % slr_z(:,i,j), &
                                                  mu_left  = [mu, 0.0_RP, 0.0_RP], &
                                                  mu_right = [mu, 0.0_RP, 0.0_RP], &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  dWall = f % geom % dWall(i,j), &
                                                  flux  = visc_flux(:,i,j) )

            end do
         end do

         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)
!      
!              --------------
!              Invscid fluxes
!              --------------
!      
               ! CALL RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
               !                    QRight = f % storage(2) % Q(:,i,j), &
               !                    nHat   = f % geom % normal(:,i,j), &
               !                    t1     = f % geom % t1(:,i,j), &
               !                    t2     = f % geom % t2(:,i,j), &
               !                    flux   = inv_flux(:,i,j) )

               
!
!              Multiply by the Jacobian
!              ------------------------
               !flux(:,i,j) = ( inv_flux(:,i,j) - visc_flux(:,i,j)) * f % geom % jacobian(i,j)
               flux(:,i,j) = (- visc_flux(:,i,j)) * f % geom % jacobian(i,j)
               
            END DO   
         END DO  
!
!        ---------------------------
!        Return the flux to elements
!        ---------------------------
!
         call f % ProjectFluxToElements(NCONS, flux, (/1,2/))

      END SUBROUTINE computeElementInterfaceFlux_SLR

      SUBROUTINE computeMPIFaceFlux_SLR(f)
         use FaceClass
         use RiemannSolvers_SLR
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j
         integer       :: thisSide
         !real(kind=RP) :: inv_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: mu
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
               call GetViscosity(f % storage(1) % Q(1,i,j), mu)

               CALL ViscousDiscretization % RiemannSolver(nEqn = NCONS, nGradEqn = NCONS, &
                                                  EllipticFlux = slrViscousFlux, &
                                                  f = f, &
                                                  QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                  U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                  U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                  U_xRight = f % storage(2) % U_x(:,i,j), &
                                                  U_yRight = f % storage(2) % U_y(:,i,j), &
                                                  U_zRight = f % storage(2) % U_z(:,i,j), &
                                                  mu_left  = [mu, 0.0_RP, 0.0_RP], &
                                                  mu_right = [mu, 0.0_RP, 0.0_RP], &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  dWall = f % geom % dWall(i,j), &
                                                  flux  = visc_flux(:,i,j) )

!
               ! CALL RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
               !                    QRight = f % storage(2) % Q(:,i,j), &
               !                    nHat   = f % geom % normal(:,i,j), &
               !                    t1     = f % geom % t1(:,i,j), &
               !                    t2     = f % geom % t2(:,i,j), &
               !                    flux   = inv_flux(:,i,j) )
!
!              Multiply by the Jacobian
!              ------------------------
               !flux(:,i,j) = ( inv_flux(:,i,j) - visc_flux(:,i,j)) * f % geom % jacobian(i,j)
               flux(:,i,j) = (- visc_flux(:,i,j)) * f % geom % jacobian(i,j)
               
            END DO   
         END DO  
!
!        ---------------------------
!        Return the flux to elements: The sign in eR % storage % FstarB has already been accouted.
!        ---------------------------
!
         thisSide = maxloc(f % elementIDs, dim = 1)
         call f % ProjectFluxToElements(NCONS, flux, (/thisSide, HMESH_NONE/))

      end subroutine ComputeMPIFaceFlux_SLR

      SUBROUTINE computeBoundaryFlux_SLR(f, time)
         USE ElementClass
         use FaceClass
         USE RiemannSolvers_SLR
         IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      type(Face),    intent(inout) :: f
      REAL(KIND=RP)                :: time
!


      END SUBROUTINE computeBoundaryFlux_SLR
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
         use PhysicsStorage, only: NCONS
         implicit none
         type(HexMesh)                  :: mesh
         real(kind=RP),      intent(in) :: time

         call ViscousDiscretization % ComputeGradient( NCONS, NCONS, mesh , time , SLRGradientVariables)

      end subroutine DGSpatial_ComputeGradient

!
!////////////////////////////////////////////////////////////////////////////////////////
!
end module SpatialDiscretization