!
!//////////////////////////////////////////////////////
!
!   @File:    SpatialDiscretization.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Sun Jan 14 17:14:44 2018
!   @Last revision date: Wed Apr 11 13:13:18 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 354405a2601df9bc6ed4885b661cc83e9e92439b
!
!//////////////////////////////////////////////////////
!
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
      use EllipticDiscretizations
      use DGWeakIntegrals
      use MeshTypes
      use HexMeshClass
      use ElementClass
      use Physics
      use PhysicsStorage
      use MPI_Face_Class
      use MPI_Process_Info
      use DGSEMClass
      use BoundaryConditionFunctions, only: C_BC, MU_BC
#ifdef _HAS_MPI_
      use mpi
#endif

      private
      public  ComputeLaplacian, DGSpatial_ComputeGradient
      public  Initialize_SpaceAndTimeMethods, ComputeTimeDerivative, ComputeTimeDerivativeIsolated
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
         character(len=LINE_LENGTH)       :: inviscidDiscretization
         character(len=LINE_LENGTH)       :: viscousDiscretization

         if ( MPI_Process % isRoot ) then
            write(STD_OUT,'(/)')
            call Section_Header("Spatial discretization scheme")
            write(STD_OUT,'(/)')
         end if
!
!        Initialize viscous discretization
!        ---------------------------------         
         call BassiRebay1     % Initialize(controlVariables)
         call BassiRebay2     % Initialize(controlVariables)
         call InteriorPenalty % Initialize(controlVariables)

         viscousDiscretization = controlVariables % stringValueForKey(viscousDiscretizationKey, requestedLength = LINE_LENGTH)
         call toLower(viscousDiscretization)
         
         select case ( trim(viscousDiscretization) )
         case("br1")
            EllipticDiscretization => BassiRebay1

         case("br2")
            EllipticDiscretization => BassiRebay2

         case("ip")
            EllipticDiscretization => InteriorPenalty

         case default
            write(STD_OUT,'(A,A,A)') 'Requested viscous discretization "',trim(viscousDiscretization),'" is not implemented.'
            write(STD_OUT,'(A)') "Implemented discretizations are:"
            write(STD_OUT,'(A)') "  * BR1"
            write(STD_OUT,'(A)') "  * BR2"
            write(STD_OUT,'(A)') "  * IP"
            errorMessage(STD_OUT)
            stop 

         end select

         call EllipticDiscretization % Describe
!
!        Compute wall distances
!        ----------------------
         call mesh % ComputeWallDistances
      
      end subroutine Initialize_SpaceAndTimeMethods
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ComputeTimeDerivative( mesh, time, BCFunctions)
         IMPLICIT NONE 
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(HexMesh), target      :: mesh
         REAL(KIND=RP)              :: time
         type(BCFunctions_t), intent(in)  :: BCFunctions(no_of_BCsets)
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER :: k, eID
         class(Element), pointer  :: e
         interface
            subroutine UserDefinedSourceTerm(mesh, time, thermodynamics_, dimensionless_, refValues_)
               use SMConstants
               USE HexMeshClass
               use PhysicsStorage
               IMPLICIT NONE
               CLASS(HexMesh)                        :: mesh
               REAL(KIND=RP)                         :: time
               type(Thermodynamics_t),    intent(in) :: thermodynamics_
               type(Dimensionless_t),     intent(in) :: dimensionless_
               type(RefValues_t),         intent(in) :: refValues_
            end subroutine UserDefinedSourceTerm
         end interface

!
!        **************************************
!        Compute chemical potential: Q stores c
!        **************************************
!
!$omp parallel shared(mesh, time) private(e)
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
         call mesh % ProlongSolutionToFaces()
!
!        ----------------
!        Update MPI Faces
!        ----------------
!
!$omp single
         call mesh % UpdateMPIFacesSolution
!$omp end single
!
!        -----------------
!        Compute gradients
!        -----------------
!
         CALL DGSpatial_ComputeGradient( mesh , time , BCFunctions(C_BC) % externalState)

!$omp single
         call mesh % UpdateMPIFacesGradients
!$omp end single
!
!        ------------------------------
!        Compute the chemical potential
!        ------------------------------
!
         call ComputeLaplacian(mesh = mesh , &
                               t    = time, &
                  externalState     = BCFunctions(C_BC) % externalState, &
                  externalGradients = BCFunctions(C_BC) % externalGradients )

         associate(c_alpha => thermodynamics % c_alpha, &
                   c_beta  => thermodynamics % c_beta    ) 
!$omp do schedule(runtime)
         do eID = 1, mesh % no_of_elements
            e => mesh % elements(eID)
            e % storage % mu = - POW2(dimensionless % eps) * e % storage % QDot(1,:,:,:)
            e % storage % c  = e % storage % Q(1,:,:,:)
            e % storage % gradC(1,:,:,:) = e % storage % U_x(1,:,:,:)
            e % storage % gradC(2,:,:,:) = e % storage % U_y(1,:,:,:)
            e % storage % gradC(3,:,:,:) = e % storage % U_z(1,:,:,:)
            call QuarticDWPDerivative(e % Nxyz, e % storage % c, c_alpha, c_beta, e % storage % mu)
            e % storage % Q(1,:,:,:) = e % storage % mu
         end do
!$omp end do
         end associate
!
!        *************************
!        Compute cDot: Q stores mu
!        *************************
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
         call mesh % ProlongSolutionToFaces()
!
!        ----------------
!        Update MPI Faces
!        ----------------
!
!$omp single
         call mesh % UpdateMPIFacesSolution
!$omp end single
!
!        -----------------
!        Compute gradients
!        -----------------
!
         CALL DGSpatial_ComputeGradient( mesh , time , BCFunctions(MU_BC) % externalState)
!
!$omp single
         call mesh % UpdateMPIFacesGradients
!$omp end single
!
!        ------------------------------
!        Compute the chemical potential
!        ------------------------------
!
         call ComputeLaplacian(mesh = mesh , &
                               t    = time, &
                  externalState     = BCFunctions(MU_BC) % externalState, &
                  externalGradients = BCFunctions(MU_BC) % externalGradients )
!
!        Scale QDot with the Peclet number
!        ---------------------------------
!$omp do schedule(runtime)
         do eID = 1, mesh % no_of_elements
            e => mesh % elements(eID)
            e % storage % QDot = (1.0_RP / Pe) * e % storage % QDot
         end do
!$omp end do

!
!        Add a source term
!        -----------------
!$omp single
         call UserDefinedSourceTerm(mesh, time, thermodynamics, dimensionless, refValues)
!$omp end single
!
!        *****************************
!        Return the concentration to Q
!        *****************************
!
!$omp do schedule(runtime)
         do eID = 1, mesh % no_of_elements
            e => mesh % elements(eID)
!            associate(e => mesh % elements(eID))
            e % storage % Q(1,:,:,:) = e % storage % c
!            end associate
         end do
!$omp end do
!$omp end parallel

      END SUBROUTINE ComputeTimeDerivative
      
      SUBROUTINE ComputeTimeDerivativeIsolated( mesh, time, BCFunctions)
         IMPLICIT NONE 
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(HexMesh), target      :: mesh
         REAL(KIND=RP)              :: time
         type(BCFunctions_t), intent(in)  :: BCFunctions(no_of_BCsets)
         
         ERROR stop 'ComputeTimeDerivativeIsolated not implemented for Cahn-Hilliard'
      END SUBROUTINE ComputeTimeDerivativeIsolated

      subroutine ComputeLaplacian( mesh , t, externalState, externalGradients )
         implicit none
         type(HexMesh)              :: mesh
         real(kind=RP)              :: t
         external                   :: externalState, externalGradients
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
               CALL computeBoundaryFlux(f, t, externalState, externalGradients) 
            
            case (HMESH_MPI)

            case default
               print*, "Unrecognized face type"
               errorMessage(STD_OUT)
               stop
                
            end select 
            end associate 
         end do 
!$omp end do 
!
!        ***************************************************************
!        Surface integrals and scaling of elements with non-shared faces
!        ***************************************************************
! 
!$omp do schedule(runtime) private(i, j, k)
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
            call mesh % GatherMPIFacesGradients
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
!$omp do schedule(runtime) 
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

      end subroutine ComputeLaplacian
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
         real(kind=RP) :: contravariantFlux  ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         integer       :: eID
!
!        *************************************
!        Compute interior contravariant fluxes
!        *************************************
!
!        Compute contravariant flux
!        --------------------------
         call EllipticDiscretization  % ComputeInnerFluxes ( e , ViscousFlux3D, contravariantFlux  ) 
!
!        ************************
!        Perform volume integrals
!        ************************
!
         e % storage % QDot = - ScalarWeakIntegrals % StdVolumeGreen ( e , contravariantFlux ) 

      end subroutine TimeDerivative_VolumetricContribution
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine TimeDerivative_FacesContribution( e , t , mesh)
         use HexMeshClass
         use PhysicsStorage
         implicit none
         type(Element)           :: e
         real(kind=RP)           :: t
         type(HexMesh)           :: mesh

         e % storage % QDot = e % storage % QDot + ScalarWeakIntegrals % StdFace( e, &
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
      SUBROUTINE computeElementInterfaceFlux(f)
         use FaceClass
         use Physics
         use PhysicsStorage
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j
         real(kind=RP) :: flux(1:N_EQN,0:f % Nf(1),0:f % Nf(2))

         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)

!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               CALL EllipticDiscretization % RiemannSolver(f = f, &
                                                  EllipticFlux = ViscousFlux0D, &
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
                                                  flux  = flux(:,i,j) )

               flux(:,i,j) = flux(:,i,j) * f % geom % jacobian(i,j)

            end do
         end do
!
!        ---------------------------
!        Return the flux to elements
!        ---------------------------
!
         call f % ProjectFluxToElements(flux, (/1,2/))

      END SUBROUTINE computeElementInterfaceFlux

      SUBROUTINE computeMPIFaceFlux(f)
         use FaceClass
         use Physics
         use PhysicsStorage
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j
         integer       :: thisSide
         real(kind=RP) :: flux(1:N_EQN,0:f % Nf(1),0:f % Nf(2))

         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)
!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               CALL EllipticDiscretization % RiemannSolver(f = f, &
                                                  EllipticFlux = ViscousFlux0D, &
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
                                                  flux  = flux(:,i,j) )

               flux(:,i,j) = flux(:,i,j) * f % geom % jacobian(i,j)

            END DO   
         END DO  
!
!        ---------------------------
!        Return the flux to elements: The sign in eR % storage % FstarB has already been accouted.
!        ---------------------------
!
         thisSide = maxloc(f % elementIDs, dim = 1)
         call f % ProjectFluxToElements(flux, (/thisSide, HMESH_NONE/))

      end subroutine ComputeMPIFaceFlux

      SUBROUTINE computeBoundaryFlux(f, time, externalStateProcedure , externalGradientsProcedure)
      USE ElementClass
      use FaceClass
      USE EllipticDiscretizations
      USE Physics
      use PhysicsStorage
      USE BoundaryConditionFunctions
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      type(Face),    intent(inout) :: f
      REAL(KIND=RP)                :: time
      procedure(BCState_FCN)     :: externalState
      procedure(BCGradients_FCN) :: externalGradients
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER                         :: i, j
      INTEGER, DIMENSION(2)           :: N
      REAL(KIND=RP)                   :: UGradExt(NDIM , N_GRAD_EQN) 
      real(kind=RP)                   :: flux(N_EQN, 0:f % Nf(1), 0:f % Nf(2))
      
      CHARACTER(LEN=BC_STRING_LENGTH) :: boundaryType
            
      boundaryType = f % boundaryType
!
!     -------------------
!     Get external states
!     -------------------
!
      do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
         f % storage(2) % Q(:,i,j) = f % storage(1) % Q(:,i,j)
         CALL externalStateProcedure( f % geom % x(:,i,j), &
                                      time, &
                                      f % geom % normal(:,i,j), &
                                      f % storage(2) % Q(:,i,j),&
                                      boundaryType )

      end do               ; end do

      do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
         UGradExt(IX,:) = f % storage(1) % U_x(:,i,j)
         UGradExt(IY,:) = f % storage(1) % U_y(:,i,j)
         UGradExt(IZ,:) = f % storage(1) % U_z(:,i,j)

         CALL externalGradientsProcedure(  f % geom % x(:,i,j), &
                                           time, &
                                           f % geom % normal(:,i,j), &
                                           UGradExt,&
                                           boundaryType )

         f % storage(1) % U_x(:,i,j) = UGradExt(IX,:)
         f % storage(1) % U_y(:,i,j) = UGradExt(IY,:)
         f % storage(1) % U_z(:,i,j) = UGradExt(IZ,:)

         f % storage(2) % U_x(:,i,j) = UGradExt(IX,:)
         f % storage(2) % U_y(:,i,j) = UGradExt(IY,:)
         f % storage(2) % U_z(:,i,j) = UGradExt(IZ,:)

!   
!           --------------
!           Viscous fluxes
!           --------------
!   
         CALL EllipticDiscretization % RiemannSolver(f = f, &
                                            EllipticFlux = ViscousFlux0D, &
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
                                            flux  = flux(:,i,j) )

         flux(:,i,j) = flux(:,i,j) * f % geom % jacobian(i,j)

      end do               ; end do

      call f % ProjectFluxToElements(flux, (/1, HMESH_NONE/))

      END SUBROUTINE computeBoundaryFlux
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!              GRADIENT PROCEDURES
!              -------------------
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine DGSpatial_ComputeGradient( mesh , time , externalStateProcedure)
         use HexMeshClass
         use PhysicsStorage
         implicit none
         type(HexMesh)                  :: mesh
         real(kind=RP),      intent(in) :: time
         procedure(BCState_FCN)         :: externalStateProcedure

         call EllipticDiscretization % ComputeGradient( mesh , time , externalStateProcedure)

      end subroutine DGSpatial_ComputeGradient
!
!////////////////////////////////////////////////////////////////////////////////////////
!
end module SpatialDiscretization
