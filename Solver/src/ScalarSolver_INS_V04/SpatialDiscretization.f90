#include "Includes.h"
module SpatialDiscretization
      use SMConstants
      ! use HyperbolicDiscretizations
      use EllipticDiscretizations
      use DGIntegrals
      use MeshTypes
      use HexMeshClass
      use ElementClass
      use PhysicsStorage
      use Physics_SLR_INS_V04
      use MPI_Face_Class
      use MPI_Process_Info
      use DGSEMClass
      use ParticlesClass
      use FluidData
      use VariableConversion_SLR_INS_V04, only: SLR_INS_V04GradientVariables, GetViscosity 
      use ProblemFileFunctions
      use BoundaryConditions, only: BCs
      use NodalStorageClass
      use ElementConnectivityDefinitions, only: axisMap, normalAxis
      use FaceClass, only: Face
#ifdef _HAS_MPI_
      use mpi
#endif

      private
      public   ComputeTimeDerivative, ComputeTimeDerivativeIsolated, viscousDiscretizationKey
      public   ComputeNonlinearStep1, TimeDerivative_ComputeQDotStep1
      public   ComputeTimeDerivative_Second
      public   ComputeTimeDerivative_Third
      public   ComputeRHSPoisson
      public   Initialize_SpaceAndTimeMethods, Finalize_SpaceAndTimeMethods
      public   computeL2Error, SCALAR_INS_V04_ComputeInnerFluxes
      public   Custom_var_ComputeGradient


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

            ! if (.not. allocated(HyperbolicDiscretization)) allocate( StandardDG_t  :: HyperbolicDiscretization )


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

               call ViscousDiscretization % Construct(controlVariables, ELLIPTIC_iNS)
               call ViscousDiscretization % Describe


!
!              Compute wall distances
!              ----------------------
               call mesh % ComputeWallDistances

          end if

      end subroutine Initialize_SpaceAndTimeMethods
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine Finalize_SpaceAndTimeMethods
         implicit none
      !   IF ( ALLOCATED(HyperbolicDiscretization) ) DEALLOCATE( HyperbolicDiscretization )
      end subroutine Finalize_SpaceAndTimeMethods


      SUBROUTINE ComputeNonlinearStep1( mesh, particles, time, dt, mode, gamma, alpha0, alpha1, beta0, beta1, HO_Elements)
         IMPLICIT NONE 
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(HexMesh), target           :: mesh
         type(Particles_t)               :: particles
         REAL(KIND=RP)                   :: time
         REAL(KIND=RP)                   :: dt
         integer, optional,   intent(in) :: mode

         real(kind=RP), optional,  intent(in)       :: gamma, alpha0, alpha1, beta0, beta1         

         logical, intent(in), optional   :: HO_Elements
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER :: i, j, k, eID
         INTEGER :: ierr, fID, iFace, iEl, iP, STLNum, n 
! 

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




!$omp do schedule(runtime) private(i,j,k)
         do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )
               IF (ANY(ISNAN(e % storage % slr2))) STOP "Error: slr2 contains NaN"
               IF (ANY(ISNAN(e % storage % Nterm2))) STOP "Error: Nterm2 contains NaN"
               IF (ANY(ISNAN(e % storage % Q))) STOP "Error: Q contains NaN"
               IF (ANY(ISNAN(e % storage % QDot))) STOP "Error: QDot contains NaN"
      
      
            do k = 0, e % Nxyz(3)   ; 
               do j = 0, e % Nxyz(2) ; 
                  do i = 0, e % Nxyz(1)


                     e % storage % Q   (1:1+N_INS-1,i,j,k) = e % storage % Q   (5:5+N_INS-1,i,j,k)
                     e % storage % slr2(1:1+N_INS-1,i,j,k) = e % storage % slr1(1:1+N_INS-1,i,j,k)
                     e % storage % slr1(1:1+N_INS-1,i,j,k) = e % storage % Q   (1:1+N_INS-1,i,j,k)

                  end do                  ; 
               end do                ; 
            end do
            end associate
         end do
!$omp end do
         


!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
         call mesh % ProlongSolutionToFaces(NCONS)

!        ****************
!        Volume integrals
!        ****************
!
!$omp do schedule(runtime) 

         do eID = 1 , size(mesh % elements)
            call TimeDerivative_VolumetricContribution( mesh % elements(eID) , time)
         end do
!$omp end do nowait

!
!        -----------------------
!        Compute time derivative
!        -----------------------
    
!$omp do schedule(runtime) private(i,j,k)
         do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )
            do k = 0, e % Nxyz(3)   ; 
               do j = 0, e % Nxyz(2) ; 
                  do i = 0, e % Nxyz(1)
                     
                     e % storage % Nterm2(1:N_INS,i,j,k) = e % storage % QDot  (1:N_INS,i,j,k)
                     e % storage % QDot  (1:N_INS,i,j,k) = e % storage % Nterm1(1:N_INS,i,j,k) 

                  end do                  ; 
               end do                ; 
            end do
            end associate
         end do
!$omp end do
!


         !
!        ******************************************
!        Compute Riemann solver of non-shared faces
!        ******************************************
!
!$omp do schedule(runtime) private(fID)
         do iFace = 1, size(mesh % faces_interior)
            fID = mesh % faces_interior(iFace)
            call computeElementInterfaceFlux_SLR(mesh % faces(fID))
         end do
!$omp end do nowait

!$omp do schedule(runtime) private(fID)
         do iFace = 1, size(mesh % faces_boundary)
            fID = mesh % faces_boundary(iFace)
            call computeBoundaryFlux_SLR(mesh % faces(fID), time, mesh)
         end do
!$omp end do
! !



!$omp do schedule(runtime) private(i,j,k)
         do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )
            if ( e % hasSharedFaces ) cycle
            ! associate ( e => mesh % elements(eID)) 
            ! if ( .not. e % hasSharedFaces ) cycle
   
            call TimeDerivative_FacesContribution(e, time, mesh)

            write(*,*) "gamma, alpha0, alpha1, beta0, beta1 ComputeNonlinearStep1 = ", gamma, alpha0, alpha1, beta0, beta1

            do k = 0, e % Nxyz(3)   ; 
               do j = 0, e % Nxyz(2) ; 
                  do i = 0, e % Nxyz(1)

                     ! write(*,*) "e % storage % slr1(1:N_INS,i,j,k) ====",e % storage % slr1(1:N_INS,i,j,k)
                     ! write(*,*) "e % storage % slr2(1:N_INS,i,j,k) ====",e % storage % slr2(1:N_INS,i,j,k)
                     ! write(*,*) "e % storage % QDot  (1:N_INS,i,j,k) ====",e % storage % QDot  (1:N_INS,i,j,k)
                     ! write(*,*) "231 e % storage % Nterm2(1:N_INS,i,j,k) ====",e % storage % Nterm2(1:N_INS,i,j,k)
                     ! write(*,*) "e % geom % jacobian(i,j,k) ====",e % geom % jacobian(i,j,k)
                     ! write(*,*) "time ====",time
                     ! write(*,*) "dt ====",dt

                     e % storage % Q(1:N_INS,i,j,k) =   1.0_RP / gamma * &
                                                      ( alpha0 * e % storage % Q(1:N_INS,i,j,k) &
                                                      + alpha1 * e % storage % slr2(1:N_INS,i,j,k) &
                                                      + dt * beta0 * e % storage % QDot  (1:N_INS,i,j,k) / e % geom % jacobian(i,j,k)  &
                                                      + dt * beta1 * e % storage % Nterm2(1:N_INS,i,j,k) / e % geom % jacobian(i,j,k)  &
                                                      )
 
                  end do             
               end do                
            end do
            end associate
         end do
!$omp end do
                                       
! =================================
!$omp end parallel
!     


      END SUBROUTINE ComputeNonlinearStep1
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
!        -----------------
!        Compute gradients
!        -----------------
!
         if ( computeGradients ) then
            CALL DGSpatial_ComputeGradient(mesh , time)
         end if   


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
      SUBROUTINE ComputeTimeDerivative_Third( mesh, particles, time, dt, mode, gamma, alpha0, alpha1, beta0, beta1, HO_Elements)
!
!        ---------
!        Arguments
!        ---------
!
         use SMConstants
         use HexMeshClass
         use ParticlesClass
         IMPLICIT NONE
         type(HexMesh), target           :: mesh
 
         type(Particles_t)               :: particles
 
         REAL(KIND=RP)                   :: time
         REAL(KIND=RP)                   :: dt
         integer,   optional,     intent(in) :: mode
         real(kind=RP), optional,  intent(in)       :: gamma, alpha0, alpha1, beta0, beta1         
         logical, intent(in), optional   :: HO_Elements
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER :: k, eID, fID


!$omp parallel shared(mesh, time) private(k, eID)

!        ------------------------------
!        Change memory to scalar
!        ------------------------------
!
!$omp do schedule(runtime)
         do eID = 1, mesh % no_of_elements
            call mesh % elements(eID) % storage % SetStorageToSLR_slr
         end do
!$omp end do

         call Custom_var_ComputeGradient(1,1,4,mesh , time )

         ! write (*,*) Custom_var_ComputeGradient, 
         write(*,*) "gamma, alpha0, alpha1, beta0, beta1 ComputeTimeDerivative_Third = ", gamma, alpha0, alpha1, beta0, beta1


         call TimeDerivative_ComputeQDot_Third(mesh = mesh , &
                                       particles = particles, &
                                       t    = dt,  &
                                       gamma=gamma, alpha0=alpha0, alpha1= alpha1, beta0 = beta0, beta1 = beta1  &
                                       )


!$omp end parallel
!     


      END SUBROUTINE ComputeTimeDerivative_Third


!       SUBROUTINE ComputeTimeDerivative_Second( mesh, particles, time, mode, HO_Elements)
!          IMPLICIT NONE 
! !
! !        ---------
! !        Arguments
! !        ---------
! !
!          TYPE(HexMesh), target           :: mesh
!          type(Particles_t)               :: particles
!          REAL(KIND=RP)                   :: time
!          integer,             intent(in) :: mode
!          logical, intent(in), optional   :: HO_Elements
      SUBROUTINE ComputeTimeDerivative_Second( mesh, particles, time, dt, mode, gamma, alpha0, alpha1, beta0, beta1, HO_Elements)

   !        ---------
   !        Arguments
   !        ---------
   !
            use SMConstants
            use HexMeshClass
            use ParticlesClass
            IMPLICIT NONE
            type(HexMesh), target           :: mesh
    
            type(Particles_t)               :: particles
    
            REAL(KIND=RP)                   :: time
            REAL(KIND=RP)                   :: dt
            integer,  optional,    intent(in) :: mode
            real(kind=RP), optional,  intent(in)       :: gamma, alpha0, alpha1, beta0, beta1         
            logical, intent(in), optional   :: HO_Elements

         !
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER :: k, eID, fID


!$omp parallel shared(mesh, time) private(k, eID)

!        ------------------------------
!        Change memory to scalar
!        ------------------------------
!
!$omp do schedule(runtime)
         do eID = 1, mesh % no_of_elements
            call mesh % elements(eID) % storage % SetStorageToSLR_slr
         end do
!$omp end do


         ! call ViscousDiscretization % ComputeGradient( NCONS, NCONS, mesh , time , SLR_INS_V04GradientVariables)
         call Custom_var_ComputeGradient(N_INS,N_INS,1,mesh , time )
         call Custom_var_ComputeGradient(N_INS,N_INS,5,mesh , time )

         write(*,*) "gamma, alpha0, alpha1, beta0, beta1 ComputeTimeDerivative_Second = ", gamma, alpha0, alpha1, beta0, beta1

         call TimeDerivative_ComputeQDot_second(mesh = mesh , &
                                       particles = particles, &
                                       t    = dt,  &
                                       gamma=gamma, alpha0=alpha0, alpha1= alpha1, beta0 = beta0, beta1 = beta1  &
                                       )


!$omp end parallel
!     


      END SUBROUTINE ComputeTimeDerivative_Second
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ComputeRHSPoisson( mesh, particles, time, mode, HO_Elements)
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
!        -----------------
!        Compute gradients
!        -----------------
!
         if ( computeGradients ) then
            CALL DGSpatial_ComputeGradient(mesh , time)
         end if   


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


      END SUBROUTINE ComputeRHSPoisson
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
      subroutine TimeDerivative_ComputeQDotStep1( mesh , particles, t, dt )
         use HexMeshClass
         use ElementClass

         implicit none
         ! type(HexMesh), target, intent(inout) :: mesh
         type(HexMesh), intent(inout)         :: mesh
         type(Particles_t)          :: particles
         real(kind=RP)              :: t
         real(kind=RP)              :: dt
         procedure(UserDefinedSourceTermNS_f) :: UserDefinedSourceTermNS
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: eID , i, j, k, ierr, fID

         type(Element), pointer :: e
         type(NodalStorage_t), pointer :: spAxi, spAeta, spAzeta
         real(kind=RP)              :: rTerm(NCONS)
         integer ::  elSide, side
         ! type(Face), target, intent(in)    :: f       !<  Face connecting elements
         type(Face)   , pointer :: f

!        **************************************************************
!        Surface integrals and scaling of elements without shared faces
!        **************************************************************
! 
!$omp do schedule(runtime) private(i,j,k)
         do eID = 1, size(mesh % elements) 
            associate(e => mesh % elements(eID)) 
            if ( e % hasSharedFaces ) cycle
            ! call TimeDerivative_FacesContribution(e, t, mesh) 


            ! mesh % elements(eID) % storage % Q(:,:,:,:) = mesh % elements(eID) % storage % Q(:,:,:,:) + 0.1

 

            do k = 0, e % Nxyz(3)
               do j = 0, e % Nxyz(2)
                  do i = 0, e % Nxyz(1) 
                     write(*,*)  " == ****-e % storage % Q(1:N_INS,i,j,k)              ---",e % storage % Q   (:,i,j,k) 
                     write(*,*)  " == ****-e % storage % slr(1:N_INS,i,j,k)            ---",e % storage % slr (:,i,j,k) 
                     
                     e % storage % QDot(1:N_INS,i,j,k) = e % storage % Nterm2(1:N_INS,i,j,k) 
                     ! e % storage % QDot(1:N_INS,i,j,k) = e % storage % Nterm2(1:N_INS,i,j,k) / e % geom % jacobian(i,j,k)
                     ! e % storage % Q(1:N_INS,i,j,k) = e % storage % Q(1:N_INS,i,j,k) + dt * e % storage % QDot(1:N_INS,i,j,k)
                     e % storage % Q(1,i,j,k) = e % storage % Q(1,i,j,k) + dt * e % storage % QDot(1,i,j,k) / e % geom % jacobian(i,j,k)
                     ! mesh % elements(eID) % storage % slr(3,i,j,k) = mesh % elements(eID) % storage % slr(3,i,j,k) + 0.1

                     
                     write(*,*)  " == ****---TimeDerivative_ComputeQDotStep1-----N_INS----",e % storage % QDot(:,i,j,k) 
                     write(*,*)  " == ****-e % storage % Q(1:N_INS,i,j,k)              ---",e % storage % Q   (:,i,j,k) 
                     write(*,*)  " == ****-e % storage % slr(1:N_INS,i,j,k)            ---",e % storage % slr (:,i,j,k) 
                  end do
                  write(*,*)
               end do       
               write(*,*)
            end do
            



            end associate 

         end do
!$omp end do

      end subroutine TimeDerivative_ComputeQDotStep1
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine TimeDerivative_ComputeQDot( mesh , particles, t )
         implicit none
         type(HexMesh), target    , intent(inout) :: mesh
      !   type(HexMesh)              :: mesh
         type(Particles_t)          :: particles
         real(kind=RP)              :: t
         procedure(UserDefinedSourceTermNS_f) :: UserDefinedSourceTermNS
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: eID , i, j, k, ierr, fID

         type(Element), pointer :: e
         type(NodalStorage_t), pointer :: spAxi, spAeta, spAzeta
         real(kind=RP)              :: rTerm(NCONS)
         integer ::  elSide, side
         ! type(Face), target, intent(in)    :: f       !<  Face connecting elements
         type(Face)   , pointer :: f
 
!        **************************************************************
!        Surface integrals and scaling of elements without shared faces
!        **************************************************************
! 
!$omp do schedule(runtime) private(i,j,k)
         do eID = 1, size(mesh % elements) 
            associate(e => mesh % elements(eID)) 
            if ( e % hasSharedFaces ) cycle
 
            spAxi => NodalStorage(e%Nxyz(1))
            spAeta => NodalStorage(e%Nxyz(2))
            spAzeta => NodalStorage(e%Nxyz(3))

            do k = 0, e % Nxyz(3)
               do j = 0, e % Nxyz(2)
                  do i = 0, e % Nxyz(1) 
                     call UserDefinedSourceTermNS(e % geom % x(:,i,j,k), e % storage % Q(:,i,j,k),        &
                                                      t, e % storage % S_SLR(:,i,j,k), thermodynamics,    &
                                                      dimensionless, refValues)
                     
                     ! write(*,*)  " == ****12345678**** rTerm==== ",e % storage % QDot(:,i,j,k) 
                                                      ! write (*,*) ", i,j,k, e % storage % QDot(:,i,j,k)", i,j,k, e % storage % QDot(:,i,j,k)
                     ! write(*,'(A1 1F12.8 A2)', advance="no") " ", e % storage % QDot(:,i,j,k), ", "
                     rTerm = e % storage % S_SLR(:,i,j,k)
                     rTerm = rTerm  
                     ! write(*,*)  " == ****12345678***-----------------* rTerm==== ",rTerm
                     ! rTerm = rTerm * e % geom % jacobian(i,j,k) * spAxi%w(i)* spAeta%w(j)* spAzeta%w(k)
                     e % storage % QDot(N_INS:NCONS,i,j,k) = rTerm(N_INS:NCONS)
                     ! write(*,*)  " == e % storage % QDot(N_INS:NCONS,i,j,k)0=================== ", e % storage % QDot(:,i,j,k)
                     e % storage % pre_source(1,i,j,k) = rTerm(4)
                     e % storage % vel_source(1:N_INS,i,j,k) = rTerm(5:7)
                     ! write(*,*)  " == ****12345678***--------e % storage % pre_source(4:4,i,j,k)= ", e % storage % pre_source(:,i,j,k)
                     ! e % storage % QDot(:,i,j,k) = rTerm
                     ! write(*,*)  " e % storage % S_NS(:,i,j,k)",e % storage % S_SLR(:,i,j,k)
                     ! write(*,*)  " == rTerm==== ",e % storage % QDot(:,i,j,k) 
                     ! write(*,'(A1 1F26.13 A2)', advance="no") " ", e % storage % QDot(:,i,j,k), ", "
                     ! write (*,*) ", i,j,k, e % storage % QDot(:,i,j,k)", i,j,k, e % storage % QDot(:,i,j,k)
                  end do
                  write(*,*)
               end do       
               write(*,*)
            end do
   !        One block for every neighbor element
    !        ------------------------------------
         !    do elSide = 1, 6
         !       ! write (*,*) "eID, fe_plus % NumberOfConnections(elSide) ======= ", eID, e_plus % NumberOfConnections(elSide)
         !       if (e % NumberOfConnections(elSide) == 0) then 

         !           fID  = e % faceIDs(elSide)
         !           side = e % faceSide(elSide)

         !           write (*,*) "fID, side ====BCBCBCBC RHSRHSRHS====**********==**********========= ", fID, side

         !           f => mesh % faces(fID)


         !           if(f%FaceType == HMESH_BOUNDARY) then
         !               call Local_Get_BC_Dir_RHS(f,e,t,side)
         !           endif
           
         !       endif
          
         !   end do
            
            end associate 

         end do
!$omp end do

      end subroutine TimeDerivative_ComputeQDot
!////////////////////////////////////////////////////////////////////////
!
      subroutine TimeDerivative_ComputeQDot_Third( mesh , particles, t,gamma, alpha0, alpha1, beta0, beta1)
         implicit none
         type(HexMesh), target    , intent(inout) :: mesh
      !   type(HexMesh)              :: mesh
         type(Particles_t)          :: particles
         real(kind=RP)              :: t

         real(kind=RP), optional,  intent(in)       :: gamma, alpha0, alpha1, beta0, beta1         

         procedure(UserDefinedSourceTermNS_f) :: UserDefinedSourceTermNS
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: eID , i, j, k, ierr, fID

         type(Element), pointer :: e
         type(NodalStorage_t), pointer :: spAxi, spAeta, spAzeta
         real(kind=RP)              :: rTerm(NCONS)
         integer ::  elSide, side
         type(Face)   , pointer :: f
 
!        **************************************************************
!        Surface integrals and scaling of elements without shared faces
!        **************************************************************
         write(*,*) "gamma, alpha0, alpha1, beta0, beta1 TimeDerivative_ComputeQDot_Third = ", gamma, alpha0, alpha1, beta0, beta1
! 
!$omp do schedule(runtime) private(i,j,k)
         do eID = 1, size(mesh % elements) 
            associate(e => mesh % elements(eID)) 
            if ( e % hasSharedFaces ) cycle
 
            spAxi => NodalStorage(e%Nxyz(1))
            spAeta => NodalStorage(e%Nxyz(2))
            spAzeta => NodalStorage(e%Nxyz(3))

            ! e % storage % slr1(1:N_INS,:,:,:) = e % storage % Q(5:7,:,:,:)

            do k = 0, e % Nxyz(3)
               do j = 0, e % Nxyz(2)
                  do i = 0, e % Nxyz(1) 
                     ! call UserDefinedSourceTermNS(e % geom % x(:,i,j,k), e % storage % Q(:,i,j,k),        &
                     !                                  t, e % storage % S_SLR(:,i,j,k), thermodynamics,    &
                     !                                  dimensionless, refValues)
                     
                     ! e % storage % vel_source(1,i,j,k) =   1.0_RP * alpha0 /t *e % storage % Q(5,i,j,k) + 1.0_RP * alpha1 /t *e % storage % slr2(1,i,j,k) ! - 1.0_RP * e % storage % U_x(4,i,j,k) 
                     ! e % storage % vel_source(2,i,j,k) =   1.0_RP * alpha0 /t *e % storage % Q(6,i,j,k) + 1.0_RP * alpha1 /t *e % storage % slr2(2,i,j,k) ! - 1.0_RP * e % storage % U_y(4,i,j,k) 
                     ! e % storage % vel_source(3,i,j,k) =   1.0_RP * alpha0 /t *e % storage % Q(7,i,j,k) + 1.0_RP * alpha1 /t *e % storage % slr2(3,i,j,k) ! - 1.0_RP * e % storage % U_z(4,i,j,k) 
                     ! e % storage % vel_source(1,i,j,k) =  gamma*1.0_RP/t *e % storage % Q(5,i,j,k) - 1.0_RP * e % storage % U_x(4,i,j,k) 
                     e % storage % vel_source(1,i,j,k) =   1.0_RP * gamma /t *e % storage % Q(1,i,j,k) - 1.0_RP * e % storage % U_x(4,i,j,k) 
                     e % storage % vel_source(2,i,j,k) =   1.0_RP * gamma /t *e % storage % Q(2,i,j,k) - 1.0_RP * e % storage % U_y(4,i,j,k) 
                     e % storage % vel_source(3,i,j,k) =   1.0_RP * gamma /t *e % storage % Q(3,i,j,k) - 1.0_RP * e % storage % U_z(4,i,j,k) 
                     ! ! e % storage % vel_source(1,i,j,k) =  gamma*1.0_RP/t *e % storage % Q(5,i,j,k) - 1.0_RP * e % storage % U_x(4,i,j,k) 
                     ! e % storage % vel_source(2,i,j,k) =  gamma*1.0_RP/t *e % storage % Q(6,i,j,k) - 1.0_RP * e % storage % U_y(4,i,j,k) 
                     ! e % storage % vel_source(3,i,j,k) =  gamma*1.0_RP/t *e % storage % Q(7,i,j,k) - 1.0_RP * e % storage % U_z(4,i,j,k) 
                     ! e % storage % vel_source(1,i,j,k) =  1.0_RP/t *e % storage % Q(5,i,j,k) - 1.0_RP * e % storage % U_x(4,i,j,k) 

                     ! e % storage % vel_source(1,i,j,k) =  gamma*1.0_RP/t *e % storage % Q(1,i,j,k) - 1.0_RP * e % storage % U_x(4,i,j,k) 
                     ! e % storage % vel_source(2,i,j,k) =  gamma*1.0_RP/t *e % storage % Q(2,i,j,k) - 1.0_RP * e % storage % U_y(4,i,j,k) 
                     ! e % storage % vel_source(3,i,j,k) =  gamma*1.0_RP/t *e % storage % Q(3,i,j,k) - 1.0_RP * e % storage % U_z(4,i,j,k) 
                     ! ! e % storage % vel_source(1,i,j,k) =  1.0_RP/t *e % storage % Q(5,i,j,k) - 1.0_RP * e % storage % U_x(4,i,j,k) 
                     ! e % storage % vel_source(2,i,j,k) =  1.0_RP/t *e % storage % Q(6,i,j,k) - 1.0_RP * e % storage % U_y(4,i,j,k) 
                     ! e % storage % vel_source(3,i,j,k) =  1.0_RP/t *e % storage % Q(7,i,j,k) - 1.0_RP * e % storage % U_z(4,i,j,k) 

                  end do
                  ! write(*,*)
               end do       
               ! write(*,*)
            end do
            
            end associate 

         end do
!$omp end do

      end subroutine TimeDerivative_ComputeQDot_Third


      subroutine TimeDerivative_ComputeQDot_second( mesh , particles, t,gamma, alpha0, alpha1, beta0, beta1)
         implicit none
         type(HexMesh), target    , intent(inout) :: mesh
      !   type(HexMesh)              :: mesh
         type(Particles_t)          :: particles
         real(kind=RP)              :: t

         real(kind=RP), optional,  intent(in)       :: gamma, alpha0, alpha1, beta0, beta1         


         procedure(UserDefinedSourceTermNS_f) :: UserDefinedSourceTermNS
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: eID , i, j, k, ierr, fID

         type(Element), pointer :: e
         type(NodalStorage_t), pointer :: spAxi, spAeta, spAzeta
         real(kind=RP)              :: rTerm(NCONS)
         integer ::  elSide, side
         ! type(Face), target, intent(in)    :: f       !<  Face connecting elements
         type(Face)   , pointer :: f

         ! real
         real(kind=RP)              :: div_u0
         real(kind=RP)              :: div_ur
 
!        **************************************************************
!        Surface integrals and scaling of elements without shared faces
!        **************************************************************
         write(*,*) "gamma, alpha0, alpha1, beta0, beta1 TimeDerivative_ComputeQDot_second = ", gamma, alpha0, alpha1, beta0, beta1
! 
!$omp do schedule(runtime) private(i,j,k)
         do eID = 1, size(mesh % elements) 
            associate(e => mesh % elements(eID)) 
            if ( e % hasSharedFaces ) cycle
 
            spAxi => NodalStorage(e%Nxyz(1))
            spAeta => NodalStorage(e%Nxyz(2))
            spAzeta => NodalStorage(e%Nxyz(3))

            do k = 0, e % Nxyz(3)
               do j = 0, e % Nxyz(2)
                  do i = 0, e % Nxyz(1) 
                     call UserDefinedSourceTermNS(e % geom % x(:,i,j,k), e % storage % Q(:,i,j,k),        &
                                                      t, e % storage % S_SLR(:,i,j,k), thermodynamics,    &
                                                      dimensionless, refValues)

                     div_u0 = -1.0_RP * gamma / t *(                                         &
                                      e % storage % U_x(1,i,j,k)           &
                                    + e % storage % U_y(2,i,j,k)           &
                                    + e % storage % U_z(3,i,j,k)            &
                                    )

                     
                     div_ur = -1.0 * gamma  / t *(                                         &
                                    e % storage % U_x(5,i,j,k)           &
                                  + e % storage % U_y(6,i,j,k)           &
                                  + e % storage % U_z(7,i,j,k)            &
                                  )
                     
                     ! rTerm = e % storage % S_SLR(:,i,j,k)
                     write (*,*) "t, div_u0, div_u1 ======= ", t, div_u0, div_ur
                     ! write (*,*) "t, gamma, div_u0, div_u1 ======= ", t,gamma, div_u0, div_ur

                     e % storage % pre_source(1,i,j,k) = div_u0

                  end do
                  ! write(*,*)
               end do       
               ! write(*,*)
            end do
            
            end associate 

         end do
!$omp end do

      end subroutine TimeDerivative_ComputeQDot_second
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
         real(kind=RP) :: inviscidContravariantFlux ( 1:N_INS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         real(kind=RP) :: contravariantFlux         ( 1:N_INS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         integer       :: eID
!
!
!

         call SCALAR_INS_V04_ComputeInnerFluxes( e, inviscidContravariantFlux ) 

       
!        ************************
!        Perform volume integrals
!        ************************

         contravariantFlux = inviscidContravariantFlux
         ! contravariantFlux = - viscousContravariantFlux



         e % storage % Nterm1 = ScalarWeakIntegrals % StdVolumeGreen ( e, N_INS, contravariantFlux ) 

         
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

         e % storage % QDot(1:N_INS,:,:,:) = e % storage % QDot(1:N_INS,:,:,:) - ScalarWeakIntegrals % StdFace( e, N_INS, &
                                             mesh % faces(e % faceIDs(EFRONT))  % storage(e % faceSide(EFRONT))  % fStar( 1:N_INS,:,: ), &
                                             mesh % faces(e % faceIDs(EBACK))   % storage(e % faceSide(EBACK))   % fStar( 1:N_INS,:,: ), &
                                             mesh % faces(e % faceIDs(EBOTTOM)) % storage(e % faceSide(EBOTTOM)) % fStar( 1:N_INS,:,: ), &
                                             mesh % faces(e % faceIDs(ERIGHT))  % storage(e % faceSide(ERIGHT))  % fStar( 1:N_INS,:,: ), &
                                             mesh % faces(e % faceIDs(ETOP))    % storage(e % faceSide(ETOP))    % fStar( 1:N_INS,:,: ), &
                                             mesh % faces(e % faceIDs(ELEFT))   % storage(e % faceSide(ELEFT))   % fStar( 1:N_INS,:,: ) )

         ! write (*,*) "===============TimeDerivative_FacesContribution( e , t , mesh)================="

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
         use RiemannSolvers_SLR_INS_V04
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f  
         integer       :: i, j
         
         real(kind=RP) :: inv_flux(1:N_INS,0:f % Nf(1),0:f % Nf(2)) 
         real(kind=RP) :: fStar(N_INS, 0:f % Nf(1), 0: f % Nf(2))

         integer       :: Sidearray(2)

!              --------------
!              Invscid fluxes
!              --------------

         do j = 0, f % Nf(2)
            do i = 0, f % Nf(1)
!
!              --------------
!              Invscid fluxes
!              --------------
! ! !
               call LxFRiemannSolver_burges( QLeft  = f % storage(1) % Q(1:N_INS,i,j), &
                                             QRight = f % storage(2) % Q(1:N_INS,i,j), &
                                             nHat   = f % geom % normal(:,i,j), &
                                             t1     = f % geom % t1(:,i,j), &
                                             t2     = f % geom % t2(:,i,j), &
                                             flux   = inv_flux(1:N_INS,i,j)       &
                                             )
             
! !
!              Multiply by the Jacobian
               ! write(*,*) "i,j, flux   = inv_flux ========== Inner ",i, j, inv_flux(:,i,j)
!              ------------------------
               ! flux(:,i,j) = ( inv_flux(:,i,j) ) * f % geom % jacobian(i,j)
               fStar(1:N_INS,i,j) = ( inv_flux(1:N_INS,i,j) ) * f % geom % jacobian(i,j)

            end do
         end do
!
!        ---------------------------
!        Return the flux to elements
!        ---------------------------
!
         Sidearray = (/1,2/)
         call f % ProjectFluxToElements(N_INS, fstar, Sidearray)


      END SUBROUTINE computeElementInterfaceFlux_SLR

      SUBROUTINE computeMPIFaceFlux_SLR(f)
         use FaceClass
         use RiemannSolvers_iNS
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   

      end subroutine ComputeMPIFaceFlux_SLR

      SUBROUTINE computeBoundaryFlux_SLR(f, time, mesh)
      ! SUBROUTINE computeBoundaryFlux    (f, time, mesh)
         USE ElementClass
         use FaceClass
         USE RiemannSolvers_SLR_INS_V04
         IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      type(Face),    intent(inout) :: f
      REAL(KIND=RP)                :: time
      type(HexMesh), intent(in)    :: mesh
      !
!     ---------------
!     Local variables
!     ---------------
!
         INTEGER                         :: i, j
         INTEGER, DIMENSION(2)           :: N
         ! REAL(KIND=RP)                   :: inv_flux(N_INS)
         real(kind=RP)                   :: visc_flux(N_INS, 0:f % Nf(1), 0:f % Nf(2))
         real(kind=RP)                   :: fStar(N_INS, 0:f % Nf(1), 0: f % Nf(2))
         real(kind=RP)                   :: inv_flux(1:N_INS,0:f % Nf(1),0:f % Nf(2)) 

         integer       :: Sidearray(2)
         !     -------------------
!     Get external states
!     -------------------
!
      do j = 0, f % Nf(2)  ; 
         do i = 0, f % Nf(1)
            f % storage(2) % Q(1:N_INS,i,j) = f % storage(1) % Q(1:N_INS,i,j)
            ! write (*,*) "=== Q2  ======bc value ======",  f % storage(2) % Q(1:N_INS,i,j)
            ! write (*,*) "=== Q1  ======bc value ======",  f % storage(1) % Q(1:N_INS,i,j) 
            ! write (*,*) "=== Q1  =====time===",  time 

            CALL BCs(f % zone) % bc % ScalarState(f % geom % x(:,i,j), &
                                         time, &
                                         f % geom % normal(:,i,j), &
                                         f % storage(2) % Q(:,i,j))

            ! write (*,*) "=======bc value ======= f % storage(2) % Q(1:N_INS,i,j) ", f % storage(2) % Q(1:N_INS,i,j)
            ! write (*,*) "=======bc value ======= f % storage(1) % Q(1:N_INS,i,j) ", f % storage(1) % Q(1:N_INS,i,j)


            ! CALL BCs(f % zone) % bc % ScalarNeumann(&
            !                              f % geom % x(:,i,j), &
            !                              time, &
            !                              f % geom % normal(:,i,j), &
            !                              f % storage(1) % Q(:,i,j), &
            !                              f % storage(1) % U_x(:,i,j), &
            !                              f % storage(1) % U_y(:,i,j), &
            !                              f % storage(1) % U_z(:,i,j), &
            !                              visc_flux(:,i,j))

         end do               ; 
      end do


      do j = 0, f % Nf(2)
         do i = 0, f % Nf(1)
!
!           Hyperbolic part
!           -------------
            ! call LxFRiemannSolver( QLeft  = f % storage(1) % Q(:,i,j), &
            ! call LxFRiemannSolver_BC_Test( QLeft  = f % storage(1) % Q(:,i,j), &
            !                            QRight = f % storage(2) % Q(:,i,j), &
            !                            nHat   = f % geom % normal(:,i,j), &
            !                            t1__     = f % geom % t1(:,i,j), &
            !                            t2__     = f % geom % t2(:,i,j), &
            !                            flux   = inv_flux)
            call LxFRiemannSolver_burges( QLeft  = f % storage(1) % Q(1:N_INS,i,j), &
                                          QRight = f % storage(2) % Q(1:N_INS,i,j), &
                                          nHat   = f % geom % normal(:,i,j), &
                                          t1     = f % geom % t1(:,i,j), &
                                          t2     = f % geom % t2(:,i,j), &
                                          flux   = inv_flux(1:N_INS,i,j))
                                          
            ! call LxFRiemannSolverTest( QLeft  = f % storage(1) % Q(:,i,j), &
            !                            QRight = f % storage(2) % Q(:,i,j), &
            !                            nHat   = f % geom % normal(:,i,j), &
            !                            t1     = f % geom % t1(:,i,j), &
            !                            t2     = f % geom % t2(:,i,j), &
            !                            flux   = inv_flux(:,i,j))
            ! write(*,*) "i, j,  flux   = inv_flux ========== BC ", i, j, inv_flux(:,i,j)

            fStar(1:N_INS,i,j) = ( inv_flux(:,i,j) ) * f % geom % jacobian(i,j)
            ! fStar(:,i,j) = 0
         end do
      end do

      Sidearray = (/1, HMESH_NONE/)
      call f % ProjectFluxToElements(N_INS, fStar, Sidearray)


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

         call ViscousDiscretization % ComputeGradient( NCONS, NCONS, mesh , time , SLR_INS_V04GradientVariables)

      end subroutine DGSpatial_ComputeGradient





    subroutine Local_Get_BC_Dir_RHS(f, e_plus,time, side)
      implicit none
      !-arguments----------------------------------------------------------------------
      type(Face), target, intent(in)    :: f       !<  Face connecting elements
      type(Element), intent(in)    :: e_plus  !<  The off-diagonal block is the contribution to this element's equations
      integer, intent(in)    :: side    !<  side of face where e_plus is
      REAL(KIND=RP)                :: time
      !-local-variables----------------------------------------------------------------
      integer :: i1, j1, k1, eq1          ! variable counters
      integer :: nXi1, nEta1              ! Number of nodes in every direction
      integer :: EtaSpa1, ZetaSpa1        ! Spacing for these two coordinate directions
      integer :: nXi2, nEta2              ! Number of nodes in every direction
      integer :: EtaSpa2, ZetaSpa2        ! Spacing for these two coordinate directions
      integer :: elSide_plus              ! Element side where f is on e⁻
      integer :: elSide_minus             ! Element side where f is on e⁺
      integer :: elInd_plus(3)            ! Element indexes on e⁺
      integer :: elInd_minus(3)           ! Element indexes on e⁻
      integer :: normAx_plus              ! Normal axis to f on e⁺
      integer :: normAx_minus             ! Normal axis to f on e⁻
      integer :: tanAx_plus(2)           ! Tangent axes to f on e⁺
      integer :: tanAx_minus(2)           ! Tangent axes to f on e⁻
      integer :: normAxSide_plus          ! Side of the normal axis that is in contact with f on e⁺
      integer :: normAxSide_minus         ! Side of the normal axis that is in contact with f on e⁻
      integer :: faceInd_plus(2)         ! Face indexes on e⁺
      integer :: faceInd_minus(2)         ! Face indexes on e⁻
      integer :: NxyFace_plus(2)         ! Polynomial orders of face on element e⁺
      integer :: NxyFace_minus(2)         ! Polynomial orders of face on element e⁻ (only needed for viscous fluxes)
      integer :: faceInd_plus2minus(2)    ! Face indexes on e⁺ passed to the reference frame of e⁻
      integer :: faceInd_minus2plus(2)    ! Face indexes on e⁻ passed to the reference frame of e⁺ (only needed for viscous fluxes)
      integer :: Deltas                   ! Number of Kronecker deltas /= 0
      integer :: elInd_plus_norm(3)       ! Element indexes on e⁺, but the index in the normal direction has been replaced by a specified value "r"
      integer :: elInd_plus_tan1(3)       ! Element indexes on e⁺, but the index in the first tangent direction has been replaced by the index on e⁻ (in the reference frame of e⁺)
      integer :: elInd_plus_tan2(3)       ! Element indexes on e⁺, but the index in the second tangent direction has been replaced by the index on e⁻ (in the reference frame of e⁺)
      type(NodalStorage_t), target  :: spA_plus(3)       ! Nodal storage in the different directions for e_plus  - local copy
      type(NodalStorage_t), target  :: spA_minus(3)      ! Nodal storage in the different directions for e_minus - local copy
      type(NodalStorage_t), pointer :: spAnorm_plus      ! Nodal storage in the direction that is normal to the face for e⁺
      type(NodalStorage_t), pointer :: spAtan1_plus      ! Nodal storage in the tangent direction "1" to the face for e⁺ (only needed for viscous fluxes)
      type(NodalStorage_t), pointer :: spAtan2_plus      ! Nodal storage in the tangent direction "2" to the face for e⁺ (only needed for viscous fluxes)
      type(NodalStorage_t), pointer :: spAnorm_minus     ! Nodal storage in the direction that is normal to the face for e⁻
      type(NodalStorage_t), pointer :: spAtan1_minus     ! Nodal storage in the tangent direction "1" to the face for e⁻ (only needed for viscous fluxes)
      type(NodalStorage_t), pointer :: spAtan2_minus     ! Nodal storage in the tangent direction "2" to the face for e⁻ (only needed for viscous fluxes)
      real(kind=RP), pointer :: dfdq(:, :, :, :)     !
      real(kind=RP), pointer :: df_dGradQ_f(:, :, :, :, :)       ! Pointer to the Jacobian with respect to gradQ on face
      real(kind=RP), pointer :: df_dGradQ_e(:, :, :, :, :, :, :) ! Pointer to the Jacobian with respect to gradQ on face
      real(kind=RP) :: dtan1, dtan2             ! Kronecker deltas in the tangent directions
      real(kind=RP) :: dtan1_minus, dtan2_minus ! Kronecker deltas in the tangent directions in the reference frame of e⁻ (only needed for viscous fluxes)
      real(kind=RP) :: a_minus
      real(kind=RP) :: normAux(NCONS, NCONS)
      real(kind=RP), pointer :: Gvec_norm(:, :, :)       ! Auxiliary vector containing values of dFv_dgradQ in the direction normal to the face
      real(kind=RP), pointer :: Gvec_tan1(:, :, :)       ! Auxiliary vector containing values of dFv_dgradQ in the first tangent direction to the face
      real(kind=RP), pointer :: Gvec_tan2(:, :, :)       ! Auxiliary vector containing values of dFv_dgradQ in the second tangent direction to the face
      real(kind=RP), allocatable :: nHat(:, :, :)
      !--------------------------------------------------------------------------------

      ! *********************************************************************************************
      ! integer :: elInd_minus(3)            ! Element indexes on e⁺
      ! integer :: elInd_minus2(3)            ! Element indexes on e⁻
      integer :: i12, j12, k12              ! i1=i2, j1=j2, k1=k2

      real(kind=RP) :: epsilon, sigma      

      integer :: elInd_plus2(3)            ! Element indexes on e⁺
      integer :: elInd_minus2(3)           ! Element indexes on e⁻
      real(kind=RP) :: JGradXi0(3), JGradXi02(3)       ! Auxiliary vector containing values of dFv_dgradQ in the second tangent direction to the face
      real(kind=RP) :: normCartesian(3)     ! Auxiliary vector containing values of dFv_dgradQ in the second tangent direction to the face
      real(kind=RP) :: normLR    ! Auxiliary vector containing values of dFv_dgradQ in the second tangent direction to the face
      real(kind=RP) :: tmp01,tmp02    ! the variables in the face after interpolation
      integer :: faceInd_plus2(2)         ! Face indexes on e⁺
      real(kind=RP) :: BcDirichlet(NCONS)     ! Boundary condition values
      real(kind=RP) :: MatEntries_RHS(NCONS)           ! Values of the matrix entries

      ! *********************************************************************************************

      !
  !     ***********
  !     Definitions
  !     ***********
  !

      epsilon = -0.0_RP
      sigma = 10_RP

      ! sigma = 2 * 3.0_RP * sigma * ( maxval(f % Nf)+1 )*( maxval(f % Nf)+2 ) / f % geom % h
      sigma =  ( maxval(f % Nf) )**2_RP *10_RP / (f % geom % h)**1
        ! sigma = 0.5_RP* sigma * ( maxval(f % Nf)+1 )*( maxval(f % Nf)+2 ) / f % geom % h
      write (*,*) "sigma ===BCBC==============", sigma
      BcDirichlet = 3.0

      ! Entry spacing for element e⁺
      nXi1 = e_plus%Nxyz(1) + 1
      nEta1 = e_plus%Nxyz(2) + 1
      EtaSpa1 = NCONS*nXi1
      ZetaSpa1 = NCONS*nXi1*nEta1
      ! write (*,*) "nXi1, nEta1, EtaSpa1, ZetaSpa1 === ", nXi1, nEta1, EtaSpa1, ZetaSpa1

      ! Entry spacing for element e⁻
      nXi2 = e_plus%Nxyz(1) + 1
      nEta2 = e_plus%Nxyz(2) + 1
      EtaSpa2 = NCONS*nXi2
      ZetaSpa2 = NCONS*nXi2*nEta2
      ! write (*,*) "nXi2, nEta2, EtaSpa2, ZetaSpa2 === ", nXi2, nEta2, EtaSpa2, ZetaSpa2

      ! Element sides
      elSide_plus = f%elementSide(side)


      ! Normal and tangent axes
      normAx_plus = normalAxis(elSide_plus)
      tanAx_plus = axisMap(:, f%elementSide(side))

      ! Side of axis where f is
      if (normAx_plus < 0) then
          normAxSide_plus = LEFT
          normAxSide_minus = RIGHT
      else
          normAxSide_plus = RIGHT
          normAxSide_minus = LEFT
         end if
      normAx_plus = abs(normAx_plus)

      ! write(*,*) "((fsideBC---normAx_plus))normAxSide_plus))", normAx_plus,normAxSide_plus
  

      ! Nodal storage
      ! --------------------------------------
      ! TODO: Why this doesn't work since ifort ver. 19.1?
      ! --------------------------------------
      ! spA_plus  = NodalStorage(e_plus  % Nxyz)
      ! spA_minus = NodalStorage(e_minus % Nxyz)

      spA_plus(1) = NodalStorage(e_plus%Nxyz(1))
      spA_plus(2) = NodalStorage(e_plus%Nxyz(2))
      spA_plus(3) = NodalStorage(e_plus%Nxyz(3))


      spAnorm_plus => spA_plus(normAx_plus)
      spAtan1_plus => spA_plus(tanAx_plus(1))
      spAtan2_plus => spA_plus(tanAx_plus(2))



      ! if (elInd_plus(tanAx_plus(1))== elInd_plus(tanAx_plus(2) ) ) then
      if (e_plus%Nxyz(normAx_plus) /= 0 ) then
      
          ! Polynomial orders
          NxyFace_plus = e_plus%Nxyz(tanAx_plus)
      
          ! write (*,*) "NxyFace_plus, NxyFace_minus=== ", NxyFace_plus, NxyFace_minus
          ! write (*,*) "elInd_plus(normAx_plus) --==--------------=== ",elInd_plus(normAx_plus)
          ! write (*,*) "Hello BC ----------------=== "
      

          if (side == LEFT) then
              normLR =  1.0_RP
          else
              normLR =    1.0_RP
          end if
          ! write (*,*) "normLR ================, ", normLR
          ! write (*,*) "f%geom%jacobian BC ================, ", f%geom%jacobian
          ! write (*,*) "f%geom%normal BC ================, ", f%geom%normal
          ! write (*,*) "normLR BC ================, ", normLR
      
          do k1 = 0, e_plus % Nxyz(3)
              do j1 = 0, e_plus % Nxyz(2)
                 do i1 = 0, e_plus % Nxyz(1) 
                    ! write (*,*) ", i,j,k, e % storage % QDot(:,i,j,k)", i,j,k, e % storage % QDot(:,i,j,k)
                    ! write(*,'(A1 1F12.8 A2)', advance="no") " ", e % storage % QDot(:,i,j,k), ", "
              
                     elInd_plus = [i1, j1, k1]
                     elInd_plus2 = [i1, j1, k1]
                     faceInd_plus  = elInd_plus ( tanAx_plus  )
                     faceInd_plus2 = elInd_plus ( tanAx_plus )

        
                     MatEntries_RHS = 0

                     CALL BCs(f % zone) % bc % ScalarState(                                              &
                                 f % geom % x(:, elInd_plus(tanAx_plus(1)),elInd_plus(tanAx_plus(2)) ), &
                                 time,                                                                     &
                                 f % geom % normal(:, elInd_plus(tanAx_plus(1)),elInd_plus(tanAx_plus(2)) ), &
                                 f % storage(2) % Q(:, elInd_plus(tanAx_plus(1)),elInd_plus(tanAx_plus(2)) ) &
                              )

                     BcDirichlet = f % storage(2) % Q(:, elInd_plus(tanAx_plus(1)),elInd_plus(tanAx_plus(2)) )

                     ! write(*,*) "BcDirichlet ==========", BcDirichlet


              
                      if(normAx_plus == 1) then 
                          JGradXi0 = e_plus%geom%jGradXi  (:,i1, j1, k1)
                          JGradXi02 = e_plus%geom%jGradXi  (:,i1, j1, k1)
                      elseIF(normAx_plus == 2) then
                          JGradXi0 = e_plus%geom%jGradEta  (:,i1, j1, k1)
                          JGradXi02 = e_plus%geom%jGradEta  (:,i1, j1, k1)
                      else 
                          JGradXi0 = e_plus%geom%jGradZeta  (:,i1, j1, k1)
                          JGradXi02 = e_plus%geom%jGradZeta  (:,i1, j1, k1)
                      endif

                  
                     MatEntries_RHS = MatEntries_RHS                                          &
                                   + epsilon                                                   &
                                   * spAtan1_plus %w (elInd_plus(tanAx_plus(1)))                            &
                                   * spAtan2_plus %w (elInd_plus(tanAx_plus(2)))   &
                                   * spAnorm_plus %vd(elInd_plus2(normAx_plus), normAxSide_plus) &
                                   * f%geom%jacobian(elInd_plus(tanAx_plus(1)), elInd_plus(tanAx_plus(2))) &
                                   * normLR                                                      &  
                                   * BcDirichlet                                                &
                                   * DOT_PRODUCT(JGradXi02                                         &
                                               / e_plus%geom%jacobian(i1, j1, k1),                &
                                               f%geom%normal(:, elInd_plus(tanAx_plus(1)), elInd_plus(tanAx_plus(2)) )  &
                                   )

                     MatEntries_RHS = MatEntries_RHS                   &
                                    + 2_RP * sigma                           &
                                    * (f%geom%jacobian(elInd_plus(tanAx_plus(1)), elInd_plus(tanAx_plus(2)))) &
                                    * spAtan1_plus  %w(elInd_plus2(tanAx_plus(1)))  &
                                    * spAtan2_plus  %w(elInd_plus2(tanAx_plus(2)))  &
                                    * spAnorm_plus  %v (elInd_plus  (normAx_plus), normAxSide_plus) &
                                    * BcDirichlet

                     MatEntries_RHS = MatEntries_RHS &
                                    / spAtan1_plus  %w(elInd_plus2(tanAx_plus(1)))  &
                                    / spAtan2_plus  %w(elInd_plus2(tanAx_plus(2)))  &
                                    / spAnorm_plus  %w(elInd_plus (normAx_plus)  )  &
                                    / e_plus%geom%jacobian(i1, j1, k1)   

                  
                     e_plus % storage % QDot(:,i1,j1,k1) = e_plus % storage % QDot(:,i1,j1,k1) + MatEntries_RHS

                 end do
              !    write(*,*)
              end do       
              ! write(*,*)
          end do 
      end if

      !     *********
      !     Finish up
      !     *********
      !
      nullify (spAnorm_plus, spAtan1_plus, spAtan2_plus, spAnorm_minus, spAtan1_minus, spAtan2_minus)

  end subroutine Local_Get_BC_Dir_RHS !





  subroutine SCALAR_INS_V04_ComputeInnerFluxes(e , contravariantFlux )
      use ElementClass
      use Physics
      use PhysicsStorage
      implicit none
      type(Element),           intent(in)  :: e
      real(kind=RP),           intent(out) :: contravariantFlux(1:N_INS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM)
   !
   !        ---------------
   !        Local variables
   !        ---------------
   !
      integer            :: i, j, k
      real(kind=RP)      :: cartesianFlux(1:N_INS, 1:NDIM)


      do k = 0, e%Nxyz(3)   ; 
         do j = 0, e%Nxyz(2)    ; 
            do i = 0, e%Nxyz(1)


               call iEulerFlux_INS_burges( e % storage % Q(1:N_INS,i,j,k), cartesianFlux(:,:) )
               ! call iEulerFlux_INSTest( e % storage % Q(:,i,j,k), cartesianFlux(:,:) )
               ! call iEulerFlux_INS( e % storage % Q(:,i,j,k), cartesianFlux(:,:) )

               ! write (*,*) "a bug here ************* ", 1.0_RP/(i-1)



               ! write (*,*) "iEulerFlux_INS ======= ", cartesianFlux

               contravariantFlux(:,i,j,k,IX) =    cartesianFlux(:,IX) * e % geom % jGradXi(IX,i,j,k)  &
                                                + cartesianFlux(:,IY) * e % geom % jGradXi(IY,i,j,k)  &
                                                + cartesianFlux(:,IZ) * e % geom % jGradXi(IZ,i,j,k)


               contravariantFlux(:,i,j,k,IY) =   cartesianFlux(:,IX) * e % geom % jGradEta(IX,i,j,k)  &
                                                + cartesianFlux(:,IY) * e % geom % jGradEta(IY,i,j,k)  &
                                                + cartesianFlux(:,IZ) * e % geom % jGradEta(IZ,i,j,k)


               contravariantFlux(:,i,j,k,IZ) =   cartesianFlux(:,IX) * e % geom % jGradZeta(IX,i,j,k)  &
                                                + cartesianFlux(:,IY) * e % geom % jGradZeta(IY,i,j,k)  &
                                                + cartesianFlux(:,IZ) * e % geom % jGradZeta(IZ,i,j,k)

            end do               ; 
         end do                ; 
      end do

   end subroutine SCALAR_INS_V04_ComputeInnerFluxes


   subroutine Custom_var_ComputeGradient(nEqn, nGradEqn,startVarNum, mesh, time, HO_Elements)
         use HexMeshClass
         use PhysicsStorage
         use Physics
         implicit none
         integer,              intent(in) :: nEqn, nGradEqn
         integer,              intent(in) :: startVarNum
         class(HexMesh)                   :: mesh
         real(kind=RP),        intent(in) :: time
         logical, intent(in), optional    :: HO_Elements
!
!        ---------------
!        Local variables
!        ---------------
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
      set_mu = .false.

         ! subroutine HexElement_ComputeLocalGradient_SLR(self, nEqn, nGradEqn, startVarNum, set_mu)
      if (HOElements) then
!$omp do schedule(runtime) private(eID)
         do i = 1 , size(mesh % HO_Elements)
            eID = mesh % HO_Elements(i)
            call mesh % elements(eID) % HexElement_ComputeLocalGradient_SLR(nEqn, nGradEqn, startVarNum, set_mu)
         end do
!$omp end do nowait
      else
         ! write (*,*) "no HOElements ********* in SpatialDiscretization"
!$omp do schedule(runtime)
         do eID = 1 , size(mesh % elements)
            call mesh % elements(eID) % HexElement_ComputeLocalGradient_SLR(nEqn, nGradEqn, startVarNum, set_mu)
         end do
!$omp end do nowait
      end if
   
   end subroutine Custom_var_ComputeGradient

  !
!////////////////////////////////////////////////////////////////////////////////////
!
!          computeL2__Error
!           -------------------------
!
!////////////////////////////////////////////////////////////////////////////////////
!
  subroutine computeL2Error( mesh )
   implicit none
   type(HexMesh), target    , intent(inout) :: mesh
!   type(HexMesh)              :: mesh
   ! type(Particles_t)          :: particles
   ! real(kind=RP)              :: t
!        ---------------
!        Local variables
!        ---------------
!
   integer     :: eID , i, j, k, ierr, fID

   type(Element), pointer :: e
   type(NodalStorage_t), pointer :: spAxi, spAeta, spAzeta
   real(kind=RP)              :: rTerm, BcDir, BcNum
   integer ::  elSide, side
   ! type(Face), target, intent(in)    :: f       !<  Face connecting elements
   type(Face)   , pointer :: f

! 
!$omp do schedule(runtime) private(i,j,k)
   do eID = 1, size(mesh % elements) 
      associate(e => mesh % elements(eID)) 
      if ( e % hasSharedFaces ) cycle
 
      ! write (*,*) " ==========TimeDerivative_ComputeQDot=============  eID, 1/(3.0_RP - eID) ===", eID, 1_RP/(3 - eID)
      write(*,*) "eID, e % Nxyz(3), e % Nxyz(2), e % Nxyz(1)  ======",&
                  eID, e % Nxyz(3), e % Nxyz(2), e % Nxyz(1)

      spAxi => NodalStorage(e%Nxyz(1))
      spAeta => NodalStorage(e%Nxyz(2))
      spAzeta => NodalStorage(e%Nxyz(3))

      do k = 0, e % Nxyz(3)
         do j = 0, e % Nxyz(2)
            do i = 0, e % Nxyz(1) 
               rTerm = spAxi%x(i)
               e % storage % slr(:,i,j,k) = rTerm
               write(*,*)  " == rTerm___RES___==== ",rTerm
               ! write(*,'(A1 1F26.13 A2)', advance="no") " ", e % storage % QDot(:,i,j,k), ", "
            end do
            write(*,*)
         end do       
         write(*,*)
      end do
    
      
      end associate 

   end do
!$omp end do

end subroutine computeL2Error




!
!////////////////////////////////////////////////////////////////////////////////////////
!
end module SpatialDiscretization