#include "Includes.h"
#if defined(NAVIERSTOKES) || (defined(INCNS))
module SurfaceIntegrals
   use SMConstants
   use PhysicsStorage
   use Physics
   use FaceClass
   use ElementClass
   use HexMeshClass
   use VariableConversion, only: Pressure
   use NodalStorageClass
#ifdef _HAS_MPI_
   use mpi
#endif
   implicit none

   private
   public   SURFACE, TOTAL_FORCE, PRESSURE_FORCE, VISCOUS_FORCE, MASS_FLOW, FLOW_RATE, PRESSURE_DISTRIBUTION
   public   ScalarSurfaceIntegral, VectorSurfaceIntegral, ScalarDataReconstruction, VectorDataReconstruction

   integer, parameter   :: SURFACE = 1
   integer, parameter   :: TOTAL_FORCE = 2
   integer, parameter   :: PRESSURE_FORCE = 3
   integer, parameter   :: VISCOUS_FORCE = 4
   integer, parameter   :: MASS_FLOW = 5
   integer, parameter   :: FLOW_RATE = 6
   integer, parameter   :: PRESSURE_DISTRIBUTION = 7
   integer, parameter   :: USER_DEFINED = 99
   !$acc declare copyin(SURFACE, TOTAL_FORCE, PRESSURE_FORCE, VISCOUS_FORCE, MASS_FLOW, FLOW_RATE, PRESSURE_DISTRIBUTION, USER_DEFINED)
!
!  ========
   contains
!  ========
!
!////////////////////////////////////////////////////////////////////////////////////////
!
!           SCALAR INTEGRALS PROCEDURES
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      function ScalarSurfaceIntegral(mesh, zoneID, integralType, iter) result(val)
!
!        -----------------------------------------------------------
!           This function computes scalar integrals, that is, those
!           in the form:
!                 val = \int \vec{v}·\vec{n}dS
!           Implemented integrals are:
!              * Surface: computes the zone surface.
!              * Mass flow: computes the mass flow across the zone.
!              * Flow: computes the volumetric flow across the zone.
!        -----------------------------------------------------------
!
         implicit none
         class(HexMesh),      intent(inout), target  :: mesh
         integer,             intent(in)    :: zoneID
         integer,             intent(in)    :: integralType, iter
         real(kind=RP)                      :: val, localval
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: zonefID, fID, eID, fIDs(6), ierr
!
!        Initialization
!        --------------
         val = 0.0_RP

         call ProlongToFaces(mesh, zoneID)
!
!        Loop the zone to get faces and elements
!        ---------------------------------------
!$acc parallel loop gang reduction(+:val) present(mesh)
!$omp do private(fID) reduction(+:val) schedule(runtime)
         do zonefID = 1, mesh % zones(zoneID) % no_of_faces
!
!           Face global ID
!           --------------
            fID = mesh % zones(zoneID) % faces(zonefID)
!
!           Compute the integral
!           --------------------
            val = val + ScalarSurfaceIntegral_Face(mesh % faces(fID), integralType)

         end do
!$omp end do
!$omp end parallel
!$acc end parallel loop

#ifdef _HAS_MPI_
         localval = val
         call mpi_allreduce(localval, val, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

      end function ScalarSurfaceIntegral

      function ScalarSurfaceIntegral_Face(f, integralType) result(val)
         !$acc routine vector
         implicit none
         class(Face),                 intent(in)     :: f
         integer,                     intent(in)     :: integralType
         real(kind=RP)                               :: val
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                       :: i, j      ! Face indices
         real(kind=RP)                 :: p
!
!        Initialization
!        --------------
         val = 0.0_RP
!
!        Perform the numerical integration
!        ---------------------------------
         select case ( integralType )
         case ( SURFACE )
!
!           **********************************
!           Computes the surface integral
!              val = \int dS
!           **********************************
!
            !$acc loop vector collapse(2) reduction(+:val)
            do j = 0, f % Nf(2) ;    do i = 0, f % Nf(1)
               val = val + NodalStorage(f % Nf(1)) % w(i) &
                         * NodalStorage(f % Nf(2)) % w(j) & 
                         * f % geom % jacobian(i,j)
            end do          ;    end do

         case ( MASS_FLOW )
!
!           ***********************************
!           Computes the mass-flow integral
!              I = \int rho \vec{v}·\vec{n}dS
!           ***********************************
!
#if defined(NAVIERSTOKES)
            !$acc loop vector collapse(2) reduction(+:val)
            do j = 0, f % Nf(2) ;    do i = 0, f % Nf(1)
!
!              Compute the integral
!              --------------------
               val = val   +(f % storage(1) % Q(IRHOU,i,j) * f % geom % normal(1,i,j)  &
                           + f % storage(1) % Q(IRHOV,i,j) * f % geom % normal(2,i,j)  &
                           + f % storage(1) % Q(IRHOW,i,j) * f % geom % normal(3,i,j) ) &
                           * NodalStorage(f % Nf(1)) % w(i) &
                           * NodalStorage(f % Nf(2)) % w(j) &
                           * f % geom % jacobian(i,j)

            end do          ;    end do
#endif
#if defined(INCNS)
            !$acc loop vector collapse(2) reduction(+:val)
            do j = 0, f % Nf(2) ;    do i = 0, f % Nf(1)
!
!              Compute the integral
!              --------------------
               val = val   +(f % storage(1) % Q(INSRHOU,i,j) * f % geom % normal(1,i,j)  &
                           + f % storage(1) % Q(INSRHOV,i,j) * f % geom % normal(2,i,j)  &
                           + f % storage(1) % Q(INSRHOW,i,j) * f % geom % normal(3,i,j) ) &
                           * NodalStorage(f % Nf(1)) % w(i) &
                           * NodalStorage(f % Nf(2)) % w(j) &
                           * f % geom % jacobian(i,j)

            end do          ;    end do
#endif
         case ( FLOW_RATE )
!
!           ***********************************
!           Computes the flow integral
!              val = \int \vec{v}·\vec{n}dS
!           ***********************************
!
#if defined(NAVIERSTOKES)
            !$acc loop vector collapse(2) reduction(+:val)
            do j = 0, f % Nf(2) ;    do i = 0, f % Nf(1)
!
!              Compute the integral
!              --------------------
               val = val + (1.0_RP / f % storage(1) % Q(IRHO,i,j))*(f % storage(1) % Q(IRHOU,i,j) * f % geom % normal(1,i,j)  &
                                          + f % storage(1) % Q(IRHOV,i,j) * f % geom % normal(2,i,j)  &
                                          + f % storage(1) % Q(IRHOW,i,j) * f % geom % normal(3,i,j) ) &
                                          * NodalStorage(f % Nf(1)) % w(i) &
                                          * NodalStorage(f % Nf(2)) % w(j) &
                                          * f % geom % jacobian(i,j)
            end do          ;    end do
#endif
#if defined(INCNS)
            !$acc loop vector collapse(2) reduction(+:val)
            do j = 0, f % Nf(2) ;    do i = 0, f % Nf(1)
!
!              Compute the integral
!              --------------------
               val = val + (1.0_RP / f % storage(1) % Q(INSRHO,i,j))*(f % storage(1) % Q(INSRHOU,i,j) * f % geom % normal(1,i,j)  &
                                             + f % storage(1) % Q(INSRHOV,i,j) * f % geom % normal(2,i,j)  &
                                             + f % storage(1) % Q(INSRHOW,i,j) * f % geom % normal(3,i,j) ) &
                                             * NodalStorage(f % Nf(1)) % w(i) &
                                             * NodalStorage(f % Nf(2)) % w(j) &
                                             * f % geom % jacobian(i,j)
            end do          ;    end do
#endif
         case ( PRESSURE_FORCE )
!
!           ***********************************
!           Computes the pressure integral
!              val = \int pdS
!           ***********************************
!
            !$acc loop vector collapse(2) reduction(+:val)
            do j = 0, f % Nf(2) ;    do i = 0, f % Nf(1)
!
!              Compute the integral
!              --------------------
               p = Pressure(f % storage(1) % Q(:,i,j))
               val = val + p * NodalStorage(f % Nf(1)) % w(i) &
                             * NodalStorage(f % Nf(2)) % w(j) &
                             * f % geom % jacobian(i,j)
            end do          ;    end do


         case ( USER_DEFINED )   ! TODO
         end select

      end function ScalarSurfaceIntegral_Face
!
!////////////////////////////////////////////////////////////////////////////////////////
!
!           VECTOR INTEGRALS PROCEDURES
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      function VectorSurfaceIntegral(mesh, zoneID, integralType, iter) result(val)
!
!        -----------------------------------------------------------
!           This function computes scalar integrals, that is, those
!           in the form:
!                 val = \int \vec{v}·\vec{n}dS
!           Implemented integrals are:
!              * Surface: computes the zone surface.
!              * Mass flow: computes the mass flow across the zone.
!              * Flow: computes the volumetric flow across the zone.
!        -----------------------------------------------------------
!
#ifdef _HAS_MPI_
         use mpi
#endif
         implicit none
         type(HexMesh),       intent(inout) :: mesh 
         integer,             intent(in)    :: zoneID
         integer,             intent(in)    :: integralType, iter
         real(kind=RP)                      :: val(NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: zonefID, fID, eID, i, j, ierr
         real(kind=RP)  :: localval(NDIM)
         real(kind=RP)  :: valx, valy, valz
         real(kind=RP)  :: localx, localy, localz
         real(kind=RP)  :: p, tau(1:NDIM, 1:NDIM)
!
!        Initialization
!        --------------
         val = 0.0_RP
         localval = 0.0_RP
         valx = 0.0_RP
         valy = 0.0_RP
         valz = 0.0_RP

         call ProlongToFaces(mesh, zoneID)

         select case ( integralType )
         case ( SURFACE )

!$acc parallel loop gang present(mesh) num_gangs(mesh % zones(zoneID) % no_of_faces) reduction(+:valx, valy, valz)
!$omp do private(fID,localVal, valx, valy, valz) reduction(+:valx, valy, valz) schedule(runtime)
         do zonefID = 1, mesh % zones(zoneID) % no_of_faces
            !
            !           Face global ID
            !           --------------
                        fID = mesh % zones(zoneID) % faces(zonefID)
            
                        localx = 0.0_RP
                        localy = 0.0_RP
                        localz = 0.0_RP
            !
            !           **********************************
            !           Computes the surface integral
            !              val = \int \vec{n} dS
            !           **********************************
            !
                        !$acc loop vector collapse(2) reduction(+:localx, localy, localz) private(val)
                        do j = 0, mesh % faces(fid) % Nf(2) ;    do i = 0, mesh % faces(fid) % Nf(1)
                           val = NodalStorage(mesh % faces(fid) % Nf(1)) % w(i) * NodalStorage(mesh % faces(fid) % Nf(2)) % w(j) &
                                     * mesh % faces(fid) % geom % jacobian(i,j) * mesh % faces(fid) % geom % normal(:,i,j)
            
                           localx = localx + val(1)
                           localy = localy + val(2)
                           localz = localz + val(3)
               
                        end do          ;    end do
            
                        valx = valx + localx
                        valy = valy + localy
                        valz = valz + localz
            
                     end do
!$omp end do
!$omp end parallel
!$acc end parallel loop
                     
         case ( TOTAL_FORCE )

!$acc parallel loop gang present(mesh) num_gangs(mesh % zones(zoneID) % no_of_faces) reduction(+:valx, valy, valz)
!$omp do private(fID,localVal, valx, valy, valz) reduction(+:valx, valy, valz) schedule(runtime)
         do zonefID = 1, mesh % zones(zoneID) % no_of_faces
            !
            !           Face global ID
            !           --------------
                        fID = mesh % zones(zoneID) % faces(zonefID)
            
                        localx = 0.0_RP
                        localy = 0.0_RP
                        localz = 0.0_RP
            !
            !           ************************************************
            !           Computes the total force experienced by the zone
            !              F = \int p \vec{n}ds - \int tau'·\vec{n}ds
            !           ************************************************
            !
                        !$acc loop vector collapse(2) reduction(+:localx, localy, localz) private(tau, localval)
                        do j = 0, mesh % faces(fid) % Nf(2) ;    do i = 0, mesh % faces(fid) % Nf(1)
            !
            !              Compute the integral
            !              --------------------
                           p = Pressure(mesh % faces(fid) % storage(1) % Q(:,i,j))
                           call getStressTensor(mesh % faces(fid) % storage(1) % Q(:,i,j),  mesh % faces(fid) % storage(1) % U_x(:,i,j),&
                                                mesh % faces(fid) % storage(1) % U_y(:,i,j),mesh % faces(fid) % storage(1) % U_z(:,i,j), tau)
            
                           localval  = ( p * mesh % faces(fid) % geom % normal(:,i,j) - matmul(tau,mesh % faces(fid) % geom % normal(:,i,j)) ) &
                                     * mesh % faces(fid) % geom % jacobian(i,j) &
                                     * NodalStorage(mesh % faces(fid) % Nf(1)) % w(i) * NodalStorage(mesh % faces(fid) % Nf(2)) % w(j)
                           
                           localx = localx + localval(1)
                           localy = localy + localval(2)
                           localz = localz + localval(3)
                        end do          ;    end do
                        
                        valx = valx + localx
                        valy = valy + localy
                        valz = valz + localz
            
         end do
!$omp end do
!$omp end parallel
!$acc end parallel loop
         
         case ( PRESSURE_FORCE )

!$acc parallel loop gang present(mesh) num_gangs(mesh % zones(zoneID) % no_of_faces) reduction(+:valx, valy, valz) 
!$omp do private(fID,localVal, valx, valy, valz) reduction(+:valx, valy, valz) schedule(runtime)
         do zonefID = 1, mesh % zones(zoneID) % no_of_faces
            !
            !           Face global ID
            !           --------------
                        fID = mesh % zones(zoneID) % faces(zonefID)
            
                        localx = 0.0_RP
                        localy = 0.0_RP
                        localz = 0.0_RP
            !
            !           ****************************************************
            !           Computes the pressure forces experienced by the zone
            !              F = \int p \vec{n}ds
            !           ****************************************************
            !           
                        !$acc loop vector collapse(2) reduction(+:localx, localy, localz) private(val)
                        do j = 0, mesh % faces(fid) % Nf(2) ;    do i = 0, mesh % faces(fid) % Nf(1)
            !
            !              Compute the integral
            !              --------------------
                           p = Pressure(mesh % faces(fid) % storage(1) % Q(:,i,j))
            
                           val = ( p * mesh % faces(fid) % geom % normal(:,i,j) ) * mesh % faces(fid) % geom % jacobian(i,j) &
                                     * NodalStorage(mesh % faces(fid) % Nf(1)) % w(i) &
                                     * NodalStorage(mesh % faces(fid) % Nf(2)) % w(j)
            
                           localx = localx + val(1)
                           localy = localy + val(2)
                           localz = localz + val(3)
                        end do          ;    end do
                        
                        valx = valx + localx
                        valy = valy + localy
                        valz = valz + localz
            
         end do
!$omp end do
!$omp end parallel
!$acc end parallel loop

         case ( VISCOUS_FORCE )

!$acc parallel loop gang present(mesh) num_gangs(mesh % zones(zoneID) % no_of_faces) reduction(+:valx, valy, valz)
!$omp do private(fID,localVal, valx, valy, valz) reduction(+:valx, valy, valz) schedule(runtime)
         do zonefID = 1, mesh % zones(zoneID) % no_of_faces
            !
            !           Face global ID
            !           --------------
                        fID = mesh % zones(zoneID) % faces(zonefID)
            
                        localx = 0.0_RP
                        localy = 0.0_RP
                        localz = 0.0_RP
            !
            !           ************************************************
            !           Computes the total force experienced by the zone
            !              F =  - \int tau'·\vec{n}ds
            !           ************************************************
            !
                        !$acc loop vector collapse(2) reduction(+:localx, localy, localz) private(tau, val)
                        do j = 0, mesh % faces(fid) % Nf(2) ;    do i = 0, mesh % faces(fid) % Nf(1)
            !
            !              Compute the integral
            !              --------------------
                           call getStressTensor(mesh % faces(fid) % storage(1) % Q(:,i,j),  mesh % faces(fid) % storage(1) % U_x(:,i,j),&
                                                mesh % faces(fid) % storage(1) % U_y(:,i,j),mesh % faces(fid) % storage(1) % U_z(:,i,j), tau)
                           
                           val = - matmul(tau,mesh % faces(fid) % geom % normal(:,i,j)) * mesh % faces(fid) % geom % jacobian(i,j) &
                                       * NodalStorage(mesh % faces(fid) % Nf(1)) % w(i) &
                                       * NodalStorage(mesh % faces(fid) % Nf(2)) % w(j)
            
                           localx = localx + val(1)
                           localy = localy + val(2)
                           localz = localz + val(3)
                        end do          ;    end do
            
                        valx = valx + localx
                        valy = valy + localy
                        valz = valz + localz
            
         end do
!$omp end do
!$omp end parallel
!$acc end parallel loop
           
         end select
            
         val(1:3) = [valx, valy, valz]
         
#ifdef _HAS_MPI_
         localVal = val
         call mpi_allreduce(localVal, val, NDIM, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

      end function VectorSurfaceIntegral
!
!////////////////////////////////////////////////////////////////////////////////////////
!
!           INTEGRALS PROCEDURES FOR IBM DATA RECONSTRUCTION
!
!                          SURFACE INTEGRALS
!
!////////////////////////////////////////////////////////////////////////////////////////
   subroutine ScalarDataReconstruction( IBM, elements, STLNum, integralType, iter, autosave, dt ) 
      use TessellationTypes
      use MappedGeometryClass
      use IBMClass
      use OrientedBoundingBox
      use KDClass
      use MPI_Process_Info
      use MPI_IBMUtilities
#ifdef _HAS_MPI_
      use mpi
#endif
!
!        -----------------------------------------------------------------------------------------
!           This function computes Scalar integrals, that is, those
!           in the form:
!                 val = \int \vec{v}·\vec{n}dS
!           The data at the boundary point (BP) is computed through a Inverse Distance Weight
!           procedure. 
!        -----------------------------------------------------------------------------------------
      implicit none
      !-arguments--------------------------------------------------------
      type(IBM_type), intent(inout) :: IBM
      type(element),  intent(inout) :: elements(:)
      integer,        intent(in)    :: integralType, STLNum, iter
      real(kind=RP),  intent(in)    :: dt
      !-local-variables-------------------------------------------------
      real(kind=rp), allocatable  :: Qsurf(:,:), U_xsurf(:,:), U_ysurf(:,:), U_zsurf(:,:)
      integer                     :: i, j
      logical                     :: found, autosave
       
      if( .not. IBM% Integral(STLNum)% compute ) return
      
      allocate( Qsurf(NCONS,IBM% NumOfInterPoints) )
      call IBM% BandPoint_state( elements, STLNum, .true. )
 
      if( IBM% stlSurfaceIntegrals(STLNum)% move ) then
         if( IBM% stlSurfaceIntegrals(STLNum)% motionType .eq. ROTATION ) then
            call IBM% stlSurfaceIntegrals(STLNum)% getRotationaMatrix( dt )
            call OBB(STLNum)% STL_rotate( IBM% stlSurfaceIntegrals(STLNum), .true. )
         elseif( IBM% stlSurfaceIntegrals(STLNum)% motionType .eq. LINEAR ) then
            call IBM% stlSurfaceIntegrals(STLNum)% getDisplacement( dt )
            call OBB(STLNum)% STL_translate( IBM% stlSurfaceIntegrals(STLNum), .true. )
         end if
      end if

      if( .not. MPI_Process% isRoot ) return 
!$omp parallel
!$omp do schedule(runtime) private(j,found)
      do i = 1, IBM% stlSurfaceIntegrals(STLNum)% NumOfObjs

         do j = 1, NumOfVertices + 4 
            call GetSurfaceState( IBM, IBM% stlSurfaceIntegrals(STLNum)% ObjectsList(i), IBM% stlSurfaceIntegrals(STLNum)% ObjectsList(i)% vertices(j), STLNum ) 
            Qsurf = IBM% BandRegion(STLNum)% Q(:,IBM% stlSurfaceIntegrals(STLNum)% ObjectsList(i)% vertices(j)% nearestPoints)   
         end do

         do j = 1, NumOfVertices + 4
            IBM% stlSurfaceIntegrals(STLNum)% ObjectsList(i)% vertices(j)% ScalarValue = IntegratedScalarValue( Q                 = Qsurf,                                                         &
                                                                                                                vertex            = IBM% stlSurfaceIntegrals(STLNum)% ObjectsList(i)% vertices(j), &
                                                                                                                normal            = IBM% stlSurfaceIntegrals(STLNum)% ObjectsList(i)% normal,      &
                                                                                                                integralType      = integralType,                                                  &
                                                                                                                InterpolationType = IBM% InterpolationType                                         )   
         end do
      end do
!$omp end do
!$omp end parallel
      if( IBM% stl(STLNum)% move ) then
         IBM% Integral(STLNum)% ListComputed = .false.
      else 
         IBM% Integral(STLNum)% ListComputed = .true.
      end if 
      
      if( autosave ) call GenerateScalarmonitorTECfile( IBM, STLNum, integralType, iter )

   end subroutine ScalarDataReconstruction
!
!////////////////////////////////////////////////////////////////////////////////////////
!
!                          VECTOR INTEGRALS
!
!////////////////////////////////////////////////////////////////////////////////////////
   subroutine VectorDataReconstruction( IBM, elements, STLNum, integralType, iter, autosave, dt )
      use TessellationTypes
      use MappedGeometryClass
      use IBMClass
      use OrientedBoundingBox
      use KDClass
      use MPI_Process_Info
      use MPI_IBMUtilities
      use omp_lib
#ifdef _HAS_MPI_
      use mpi
#endif
!
!        -----------------------------------------------------------------------------------------
!           This function computes Vector integrals, that is, those
!           in the form:
!                 val = \int \vec{v}·\vec{n}dS
!           The data at the boundary point (BP) is computed through a Inverse Distance Weight
!           procedure. 
!        -----------------------------------------------------------------------------------------
      implicit none
      !-arguments---------------------------------------------------------------------------------
      type(IBM_type), intent(inout) :: IBM
      type(element),  intent(inout) :: elements(:)
      integer,        intent(in)    :: integralType, STLNum, iter
      real(kind=RP),  intent(in)    :: dt
      !-local-variables---------------------------------------------------------------------------
      real(kind=rp), allocatable  :: Qsurf(:,:,:), U_xsurf(:,:,:), U_ysurf(:,:,:), U_zsurf(:,:,:)
      integer                     :: i, j
      logical                     :: found, autosave

      if( .not. IBM% Integral(STLNum)% compute ) return 

      allocate( Qsurf(NCONS,IBM% NumOfInterPoints,NumOfVertices + 4),   &
                U_xsurf(NCONS,IBM% NumOfInterPoints,NumOfVertices + 4), &
                U_ysurf(NCONS,IBM% NumOfInterPoints,NumOfVertices + 4), &
                U_zsurf(NCONS,IBM% NumOfInterPoints,NumOfVertices + 4)  )
      call IBM% BandPoint_state( elements, STLNum, .true. )
 
      if( .not. MPI_Process% isRoot ) return 

      if( IBM% stlSurfaceIntegrals(STLNum)% move ) then
         if( IBM% stlSurfaceIntegrals(STLNum)% motionType .eq. ROTATION ) then
            call IBM% stlSurfaceIntegrals(STLNum)% getRotationaMatrix( dt )
            call OBB(STLNum)% STL_rotate( IBM% stlSurfaceIntegrals(STLNum), .true. )
         elseif( IBM% stlSurfaceIntegrals(STLNum)% motionType .eq. LINEAR ) then
            call IBM% stlSurfaceIntegrals(STLNum)% getDisplacement( dt )
            call OBB(STLNum)% STL_translate( IBM% stlSurfaceIntegrals(STLNum), .true. )
         end if
      end if
!$omp parallel
!$omp do schedule(runtime) private(j,found,Qsurf,U_xsurf,U_ysurf,U_zsurf)
      do i = 1, IBM% stlSurfaceIntegrals(STLNum)% NumOfObjs

         do j = 1, NumOfVertices + 4 
            call GetSurfaceState( IBM, IBM% stlSurfaceIntegrals(STLNum)% ObjectsList(i), IBM% stlSurfaceIntegrals(STLNum)% ObjectsList(i)% vertices(j), STLNum ) 

            Qsurf(:,:,j)   = IBM% BandRegion(STLNum)% Q  (:,IBM% stlSurfaceIntegrals(STLNum)% ObjectsList(i)% vertices(j)% nearestPoints)   
            U_xsurf(:,:,j) = IBM% BandRegion(STLNum)% U_x(:,IBM% stlSurfaceIntegrals(STLNum)% ObjectsList(i)% vertices(j)% nearestPoints) 
            U_ysurf(:,:,j) = IBM% BandRegion(STLNum)% U_y(:,IBM% stlSurfaceIntegrals(STLNum)% ObjectsList(i)% vertices(j)% nearestPoints) 
            U_zsurf(:,:,j) = IBM% BandRegion(STLNum)% U_z(:,IBM% stlSurfaceIntegrals(STLNum)% ObjectsList(i)% vertices(j)% nearestPoints) 
         end do
         
         do j = 1, NumOfVertices + 4
            IBM% stlSurfaceIntegrals(STLNum)% ObjectsList(i)% vertices(j)% VectorValue = IntegratedVectorValue( Q                 = Qsurf(:,:,j),                                                  &
                                                                                                                U_x               = U_xsurf(:,:,j),                                                &  
                                                                                                                U_y               = U_ysurf(:,:,j),                                                &
                                                                                                                U_z               = U_zsurf(:,:,j),                                                &
                                                                                                                vertex            = IBM% stlSurfaceIntegrals(STLNum)% ObjectsList(i)% vertices(j), &
                                                                                                                normal            = IBM% stlSurfaceIntegrals(STLNum)% ObjectsList(i)% normal,      &
                                                                                                                y                 = IBM% IP_Distance,                                              &
                                                                                                                Wallfunction      = IBM% Wallfunction,                                             &
                                                                                                                integralType      = integralType,                                                  &
                                                                                                                InterpolationType = IBM% InterpolationType                                         )   
         end do
      end do 
!$omp end do
!$omp end parallel
      deallocate( Qsurf, U_xsurf, U_ysurf, U_zsurf )

      if( IBM% stl(STLNum)% move ) then
         IBM% Integral(STLNum)% ListComputed = .false.
      else
         IBM% Integral(STLNum)% ListComputed = .true.
      end if

      if( autosave ) call GenerateVectormonitorTECfile( IBM, STLNum, integralType, iter )

   end subroutine VectorDataReconstruction

   subroutine GetSurfaceState( IBM, obj, vertex, STLNum )
      use TessellationTypes
      use IBMClass
      use MPI_Process_Info
      use omp_lib
#ifdef _HAS_MPI_
      use mpi
#endif
      implicit none

      class(IBM_type),   intent(inout) :: IBM 
      type(object_type), intent(inout) :: obj
      type(point_type),  intent(inout) :: vertex
      integer,           intent(in)    :: STLNum

      real(kind=RP) :: Dist
      integer       :: i, j, k 

      if( IBM% Integral(STLNum)% ListComputed ) return 

      vertex% nearestPoints = 0
      do k = 1, IBM% NumOfInterPoints
         if( IBM% Wallfunction ) then 
            call MinimumDistancePoints( vertex% coords + IBM% IP_Distance * obj% Normal,  &
                                        IBM% rootPoints(STLNum), IBM% BandRegion(STLNum), &
                                        Dist, k, vertex% nearestPoints                    )
         else           
            call MinimumDistancePoints( vertex% coords, IBM% rootPoints(STLNum), &
                                        IBM% BandRegion(STLNum), Dist, k,        &
                                        vertex% nearestPoints                    )   
         end if 
      end do
      call GetMatrixInterpolationSystem( vertex% coords,                                    &
                                         IBM% BandRegion(STLNum)% x(vertex% nearestPoints), &
                                         vertex% invPhi,                                    &
                                         vertex% b, IBM% InterpolationType                  )

   end subroutine GetSurfaceState 

   subroutine GetSurfaceState_HO( IBM, obj, vertex, STLNum, elements, Qs, U_xs, U_ys, U_zs, gradients, found )
      use TessellationTypes
      use IBMClass
      use MPI_Process_Info
      use omp_lib
#ifdef _HAS_MPI_
      use mpi
#endif
      implicit none

      class(IBM_type),             intent(in)    :: IBM 
      type(object_type),           intent(inout) :: obj
      type(point_type),            intent(inout) :: vertex
      integer,                     intent(in)    :: STLNum
      type(element),               intent(inout) :: elements(:)
      real(kind=RP),               intent(inout) :: Qs(NCONS,1)
      real(kind=RP),      intent(inout) :: U_xs(NCONS,1), U_ys(NCONS,1), U_zs(NCONS,1)
      logical,                     intent(in)    :: gradients
      logical,                     intent(out)   :: found

      real(kind=RP)                 :: xi(NDIM)
      integer                       :: eID, i, j, k 

      Qs = 0.0_RP
      if( gradients ) then 
         U_xs = 0.0_RP; U_ys = 0.0_RP; U_zs = 0.0_RP
      end if 

      if( IBM% Integral(STLNum)% ListComputed ) then 
         if( vertex% partition .eq. MPI_Process% rank ) then 
            eID = vertex% element_index
            xi  = vertex% xi

            associate( e => elements(eID) )
   
            Qs(:,1) = elements(eID)% EvaluateSolutionAtPoint(NCONS, xi) 
         
            if( gradients ) then 
               U_xs(:,1) = elements(eID)% EvaluateGradientAtPoint(NCONS, xi, IX)
               U_ys(:,1) = elements(eID)% EvaluateGradientAtPoint(NCONS, xi, IY)
               U_zs(:,1) = elements(eID)% EvaluateGradientAtPoint(NCONS, xi, IZ)
            end if

            end associate 
            found = .true.
         else 
            found = .false. 
         end if
         return 
      end if 

      do eID = 1, size(elements)  
         associate(e => elements(eID) )
         found = e% FindPointWithCoords( vertex% coords, 0, xi )
         if( found ) then 
            vertex% element_index = eID
            vertex% partition     = MPI_Process% rank 
            vertex% xi            = xi

            Qs(:,1) = elements(eID)% EvaluateSolutionAtPoint(NCONS, xi) 
         
            if( gradients ) then 
               U_xs(:,1) = elements(eID)% EvaluateGradientAtPoint(NCONS, xi, IX)
               U_ys(:,1) = elements(eID)% EvaluateGradientAtPoint(NCONS, xi, IY)
               U_zs(:,1) = elements(eID)% EvaluateGradientAtPoint(NCONS, xi, IZ)
            end if
            exit
         end if 
         end associate
      end do

   end subroutine GetSurfaceState_HO
!
!////////////////////////////////////////////////////////////////////////////////////////
!
!           INVERSE DISTANCE WEIGHTED INTERPOLATION PROCEDURES FOR IBM DATA RECONSTRUCTION
!
!                                   SCALAR INTERPOLATION
!
!//////////////////////////////////////////////////////////////////////////////////////// 
   function IntegratedScalarValue( Q, vertex, normal, integralType, InterpolationType ) result( outvalue )
      use IBMClass
      implicit none
!
!        -----------------------------------------------------------
!           This function computes the IDW interpolat for a scalar
!           quantity in the point "Point".
!           Available scalars are:
!           Mass flow
!           Flow rate
!           Pressure
!        -----------------------------------------------------------
      !-arguments--------------------------------------------------------------
      real(kind=rp),           intent(in)    :: Q(:,:), normal(:)
      type(point_type),        intent(inout) :: vertex
      integer,                 intent(in)    :: integralType, InterpolationType
      real(kind=rp)                          :: outvalue
      !-local-variables--------------------------------------------------------
      real(kind=rp) :: Qi(NCONS), P
      integer       :: i

      outvalue = 0.0_RP
      
      select case( integralType )

         case( MASS_FLOW )
               
            do i = 1, NCONS 
               Qi(i) = GetInterpolatedValue( Q(i,:), vertex% invPhi, vertex% b, InterpolationType )
            end do 

#if defined(NAVIERSTOKES)
            outvalue = - (1.0_RP / Qi(IRHO))*(Qi(IRHOU)*normal(1) + Qi(IRHOV)*normal(2) + Qi(IRHOW)*normal(3))       
#endif
#if defined(INCNS)
            outvalue = - (1.0_RP / Qi(INSRHO))*(Qi(INSRHOU)*normal(1) + Qi(INSRHOV)*normal(2) + Qi(INSRHOW)*normal(3))       
#endif    
         case ( FLOW_RATE )

            do i = 1, NCONS 
               Qi(i) = GetInterpolatedValue( Q(i,:), vertex% invPhi, vertex% b, InterpolationType )
            end do 

#if defined(NAVIERSTOKES)
            outvalue = - (Qi(IRHOU)*normal(1) + Qi(IRHOV)*normal(2) + Qi(IRHOW)*normal(3)) 
#endif
#if defined(INCNS)
            outvalue = - (Qi(INSRHOU)*normal(1) + Qi(INSRHOV)*normal(2) + Qi(INSRHOW)*normal(3)) 
#endif
         case( PRESSURE_DISTRIBUTION )

            do i = 1, NCONS 
               Qi(i) = GetInterpolatedValue( Q(i,:), vertex% invPhi, vertex% b, InterpolationType )
            end do  

            outvalue = pressure(Qi)
         case ( USER_DEFINED )   ! TODO  

      end select 

   end function IntegratedScalarValue
!
!////////////////////////////////////////////////////////////////////////////////////////
!
!                          VECTOR INTERPOLATION
!
!////////////////////////////////////////////////////////////////////////////////////////         
   function IntegratedVectorValue( Q, U_x, U_y, U_z, vertex, normal,  &
                                   y, Wallfunction, integralType,     &
                                   InterpolationType                  ) result( outvalue )
      use IBMClass
      use VariableConversion
      use FluidData
#if defined(NAVIERSTOKES)
      use WallFunctionBC
#endif
      implicit none
!
!        -----------------------------------------------------------
!           This function computes the IDW interpolat for a vector
!           quantity in the point "Point".
!           Available scalars are:
!           Total force
!           Pressure force
!           Viscous force
!        -----------------------------------------------------------
      !-arguments-----------------------------------------------------------------
      real(kind=rp),           intent(in)    :: Q(:,:), U_x(:,:), U_y(:,:),       &
                                                U_z(:,:), normal(NDIM)
      type(point_type),        intent(inout) :: vertex
      real(kind=rp),           intent(in)    :: y
      logical,                 intent(in)    :: Wallfunction
      integer,                 intent(in)    :: integralType, InterpolationType
      real(kind=rp)                          :: outvalue(NDIM)
      !-local-variables-----------------------------------------------------------
      integer       :: i
      real(kind=rp) :: viscStress(NDIM), U(NDIM), U_t(NDIM), tangent(NDIM),   &
                       Qi(NCONS), U_xi(NCONS), U_yi(NCONS), U_zi(NCONS),      & 
                       tau(NDIM,NDIM), P, T, T_w, rho_w, mu, nu, u_II, u_tau, &
                       tau_w, kappa_                                        
      
      outvalue = 0.0_RP

      select case( integralType )

         case ( TOTAL_FORCE )

            do i = 1, NCONS 
               Qi(i) = GetInterpolatedValue( Q(i,:), vertex% invPhi, vertex% b, InterpolationType )
            end do

            P = pressure(Qi)

            if( Wallfunction ) then
#if defined(NAVIERSTOKES) 
               T  = Temperature(Qi)
               call get_laminar_mu_kappa(Qi,mu,kappa_) 
               nu = mu/Qi(IRHO)
                
               U   = Qi(IRHOU:IRHOW)/Qi(IRHO)
               U_t = U - ( dot_product(U,normal) * normal )
 
               tangent = U_t/norm2(U_t)

               u_II  = dot_product(U,tangent)
               
               u_tau = u_tau_f( u_II, y, nu, u_tau0=0.1_RP )
            
               T_w = T + (dimensionless% Pr)**(1._RP/3._RP)/(2.0_RP*thermodynamics% cp) * POW2(u_II)
               T_w = T_w * refvalues% T
               rho_w = P*refvalues% p/(thermodynamics% R * T_w)
               rho_w = rho_w/refvalues% rho
#endif
               tau_w = rho_w*POW2(u_tau)
               
               viscStress = tau_w*tangent
            else

               do i = 1, NCONS 
                  U_xi(i) = GetInterpolatedValue( U_x(i,:), vertex% invPhi, vertex% b, InterpolationType )
                  U_yi(i) = GetInterpolatedValue( U_y(i,:), vertex% invPhi, vertex% b, InterpolationType )
                  U_zi(i) = GetInterpolatedValue( U_z(i,:), vertex% invPhi, vertex% b, InterpolationType )
               end do 

               call getStressTensor(Qi, U_xi, U_yi, U_zi, tau)
               
               viscStress = matmul(tau,normal)
            end if
            
            outvalue = -P * normal + viscStress   
                  
         case( PRESSURE_FORCE )

            do i = 1, NCONS 
               Qi(i) = GetInterpolatedValue( Q(i,:), vertex% invPhi, vertex% b, InterpolationType )
            end do

            P = pressure(Qi)
            
            outvalue = -P * normal
            
         case( VISCOUS_FORCE )

            if( Wallfunction ) then
#if defined(NAVIERSTOKES) 
               T  = Temperature(Qi)
               call get_laminar_mu_kappa(Qi,mu,kappa_) 
               nu = mu/Qi(IRHO)
                
               U   = Qi(IRHOU:IRHOW)/Qi(IRHO)
               U_t = U - ( dot_product(U,normal) * normal )
 
               tangent = U_t/norm2(U_t)

               u_II  = dot_product(U,tangent)
               
               u_tau = u_tau_f( u_II, y, nu, u_tau0=0.1_RP )
            
               T_w = T + (dimensionless% Pr)**(1._RP/3._RP)/(2.0_RP*thermodynamics% cp) * POW2(u_II)
               T_w = T_w * refvalues% T
               rho_w = P*refvalues% p/(thermodynamics% R * T_w)
               rho_w = rho_w/refvalues% rho
#endif
               tau_w = rho_w*POW2(u_tau)
               
               viscStress = tau_w*tangent
            else

               do i = 1, NCONS 
                  U_xi(i) = GetInterpolatedValue( U_x(i,:), vertex% invPhi, vertex% b, InterpolationType )
                  U_yi(i) = GetInterpolatedValue( U_y(i,:), vertex% invPhi, vertex% b, InterpolationType )
                  U_zi(i) = GetInterpolatedValue( U_z(i,:), vertex% invPhi, vertex% b, InterpolationType )
               end do 
               
               call getStressTensor(Qi, U_xi, U_yi, U_zi, tau)
               
               viscStress = matmul(tau,normal)
            end if 
            
            outvalue = viscStress
            
         case ( USER_DEFINED )   ! TODO  

      end select 

   end function IntegratedVectorValue
   
   subroutine GenerateScalarmonitorTECfile( IBM, STLNum, integralType, iter )
      use MPI_Process_Info
      use TessellationTypes
      use MPI_IBMUtilities
      use IBMClass
      implicit none
      !-arguments-------------------------------------------------------
      type(IBM_type), intent(in) :: IBM
      integer,        intent(in) :: STLNum, integralType, iter 
      !-local-variables-------------------------------------------------
      real(kind=RP), allocatable :: x(:), y(:), z(:), scalar(:), &
                                    local_sum(:), global_sum(:)
      integer                    :: i, j, index, NumOfObjs
      character(len=LINE_LENGTH) :: FileName, FinalName
#ifdef _HAS_MPI_
      integer :: ierr 
#endif     
      NumOfObjs = NumOfVertices * IBM% stlSurfaceIntegrals(STLNum)% NumOfObjs

      allocate( x(NumOfObjs),     &
                y(NumOfObjs),     &
                z(NumOfObjs),     &
                scalar(NumOfObjs) )

      index = 0

      do i = 1, IBM% stlSurfaceIntegrals(STLNum)% NumOfObjs
         associate( obj => IBM% stlSurfaceIntegrals(STLNum)% ObjectsList(i) )
         do j = 1, NumOfVertices
            index = index + 1
            x(index)      = obj% vertices(j)% coords(IX)
            y(index)      = obj% vertices(j)% coords(IY)
            z(index)      = obj% vertices(j)% coords(IZ)
            scalar(index) = obj% vertices(j)% ScalarValue
         end do
         end associate
      end do
#ifdef _HAS_MPI_
      if( MPI_Process% doMPIAction ) then     
         allocate(local_sum(NumOfObjs),global_sum(NumOfObjs))
         local_sum = scalar
         call mpi_allreduce(local_sum, global_sum, NumOfObjs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr) 
         scalar = global_sum  
         deallocate(local_sum,global_sum)
      end if
#endif
      if( .not. MPI_Process% isRoot ) then 
         deallocate(x, y, z, scalar)
         return 
      end if      
      
      if( .not. MPI_Process% isRoot ) return

      select case(integralType)
         case( MASS_FLOW )
            FileName = 'MASS_FLOW_'
            write(FinalName,'(A,A,I10.10,A)')  trim(FileName),trim(OBB(STLNum)% FileName)//'_',iter,'.tec'
            call STLScalarTEC( x, y, z, scalar, STLNum, FinalName, 'MASS FLOW', '"x","y","z","MassFlow"' )
         case( FLOW_RATE )
            FileName = 'FLOW_RATE_FORCE_'
            write(FinalName,'(A,A,I10.10,A)')  trim(FileName),trim(OBB(STLNum)% FileName)//'_',iter,'.tec'
            call STLScalarTEC( x, y, z, scalar, STLNum, FinalName, 'FLOW RATE', '"x","y","z","FlowRate"' )
         case( PRESSURE_DISTRIBUTION )
            FileName = 'PRESSURE_'
            write(FinalName,'(A,A,I10.10,A)')  trim(FileName),trim(OBB(STLNum)% FileName)//'_',iter,'.tec'
            call STLScalarTEC( x, y, z, scalar, STLNum, FinalName, 'PRESSURE DISTRIBUTION', '"x","y","z","Pressure"' )
      end select

      deallocate(x, y, z, scalar)

  end subroutine GenerateScalarmonitorTECfile
  
  subroutine GenerateVectormonitorTECfile( IBM, STLNum, integralType, iter )
      use MPI_Process_Info
      use TessellationTypes
      use MPI_IBMUtilities
      use IBMClass
      implicit none
      !-arguments---------------------------------------------------------
      type(IBM_type), intent(in) :: IBM
      integer,        intent(in) :: STLNum, integralType, iter 
      !-local-variables---------------------------------------------------
      real(kind=RP), allocatable :: x(:), y(:), z(:), vector_x(:),   &
                                    vector_y(:), vector_z(:),        &
                                    local_sum(:), global_sum(:)
      character(len=LINE_LENGTH) :: FileName, FinalName
      integer                    :: index, NumOfObjs, i, j
#ifdef _HAS_MPI_
      integer :: ierr 
#endif
      NumOfObjs = NumOfVertices * IBM% stlSurfaceIntegrals(STLNum)% NumOfObjs

      allocate( x(NumOfObjs),        &
                y(NumOfObjs),        &
                z(NumOfObjs),        &
                vector_x(NumOfObjs), &
                vector_y(NumOfObjs), &
                vector_z(NumOfObjs)  )

      index = 0

      do i = 1, IBM% stlSurfaceIntegrals(STLNum)% NumOfObjs
         associate( obj => IBM% stlSurfaceIntegrals(STLNum)% ObjectsList(i) )
         do j = 1, NumOfVertices
            index = index + 1
            x(index)        = obj% vertices(j)% coords(IX)
            y(index)        = obj% vertices(j)% coords(IY)
            z(index)        = obj% vertices(j)% coords(IZ)
            vector_x(index) = obj% vertices(j)% VectorValue(IX)
            vector_y(index) = obj% vertices(j)% VectorValue(IY)
            vector_z(index) = obj% vertices(j)% VectorValue(IZ)
         end do
         end associate
      end do
#ifdef _HAS_MPI_
      if( MPI_Process% doMPIAction ) then     
         allocate(local_sum(NumOfObjs),global_sum(NumOfObjs))
         local_sum = vector_x
         call mpi_allreduce(local_sum, global_sum, NumOfObjs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr) 
         vector_x = global_sum  
         local_sum = vector_y
         call mpi_allreduce(local_sum, global_sum, NumOfObjs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr) 
         vector_y = global_sum
         local_sum = vector_z
         call mpi_allreduce(local_sum, global_sum, NumOfObjs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr) 
         vector_z = global_sum
         deallocate(local_sum,global_sum)
      end if
#endif
      if( .not. MPI_Process% isRoot ) then 
         deallocate(x, y, z, vector_x, vector_y, vector_z)
         return 
      end if

      select case(integralType)
         case( TOTAL_FORCE )
            FileName = 'TOTAL_FORCE_'
            write(FinalName,'(A,A,I10.10,A)')  trim(FileName),trim(OBB(STLNum)% FileName)//'_',iter,'.tec'            
            call STLvectorTEC( x, y, z, vector_x, vector_y, vector_z, STLNum, FinalName, 'TOTAL FORCE', '"x","y","z","Ftot_x","Ftot_y","Ftot_z"' )
         case( PRESSURE_FORCE )
            FileName = 'PRESSURE_FORCE_'
            write(FinalName,'(A,A,I10.10,A)')  trim(FileName),trim(OBB(STLNum)% FileName)//'_',iter,'.tec'
            call STLvectorTEC( x, y, z, vector_x, vector_y, vector_z, STLNum, FinalName, 'PRESSURE FORCE', '"x","y","z","Fpres_x","Fpres_y","Fpres_z"' )
         case( VISCOUS_FORCE )
            FileName = 'VISCOUS_FORCE_'
            write(FinalName,'(A,A,I10.10,A)')  trim(FileName),trim(OBB(STLNum)% FileName)//'_',iter,'.tec'
            call STLvectorTEC( x, y, z, vector_x, vector_y, vector_z, STLNum, FinalName, 'VISCOUS FORCE', '"x","y","z","Fvisc_x","Fvisc_y","Fvisc_z"' )
      end select
   
      deallocate(x, y, z, vector_x, vector_y, vector_z)

   end subroutine GenerateVectormonitorTECfile

   subroutine ProlongToFaces(mesh, zoneID)
      implicit none
      type(HexMesh),       intent(inout) :: mesh 
      integer,             intent(in)    :: zoneID

      !-------- local variables ----------
      integer :: zonefID, fID, eID

         !$acc wait 

         select case ( mesh %nodeType )
         case(1) !Gauss-Legendre
         !$acc parallel loop gang present(mesh) 
         do zonefID = 1, mesh % zones(zoneID) % no_of_faces
            fID =  mesh % zones(zoneID) % faces(zonefID)
            eID = mesh % faces(fID) % elementIDs(1)
            call HexElement_ProlongSolToFaces(mesh % elements(eID), NCONS, mesh % faces(fID), mesh % faces(fID) % elementSide(1))                        
         end do 
         !$acc end parallel loop
         
         case(2)
         !$acc parallel loop gang present(mesh)
         do zonefID = 1, mesh % zones(zoneID) % no_of_faces
            fID =  mesh % zones(zoneID) % faces(zonefID)
            eID = mesh % faces(fID) % elementIDs(1)
            call HexElement_ProlongSolToFaces_GL(mesh % elements(eID), NCONS, mesh % faces(fID), mesh % faces(fID) % elementSide(1))                        
         end do 
         !$acc end parallel loop
         end select

         if ( computeGradients ) then
            select case ( mesh %nodeType )
            case(1) !Gauss

            !$omp do schedule(runtime) private(eID)
            !$acc parallel loop gang present(mesh) 
            do zonefID = 1, mesh % zones(zoneID) % no_of_faces
               fID =  mesh % zones(zoneID) % faces(zonefID)
               eID = mesh % faces(fID) % elementIDs(1)
               call HexElement_ProlongGradientsToFaces(mesh % elements(eID), NGRAD, &
                                                       mesh % faces(fID), &
                                                       mesh % elements(eID) % storage % U_x, mesh % faces(fID) % elementSide(1), 1)

               call HexElement_ProlongGradientsToFaces(mesh % elements(eID), NGRAD, &
                                                       mesh % faces(fID), &
                                                       mesh % elements(eID) % storage % U_y, mesh % faces(fID) % elementSide(1), 2)

               call HexElement_ProlongGradientsToFaces(mesh % elements(eID), NGRAD, &
                                                       mesh % faces(fID), &
                                                       mesh % elements(eID) % storage % U_z, mesh % faces(fID) % elementSide(1), 3)
            end do
            !$acc end parallel loop
            !$omp end do

            case(2) !Gauss-Lobatto

            !$omp do schedule(runtime)
            !$acc parallel loop gang present(mesh) private(fID)
            do zonefID = 1, mesh % zones(zoneID) % no_of_faces
               fID =  mesh % zones(zoneID) % faces(zonefID)
               eID = mesh % faces(fID) % elementIDs(1)
               call HexElement_ProlongGradientsToFaces_GL(mesh % elements(eID), NGRAD, &
                                                          mesh % faces(fID), &
                                                          mesh % elements(eID) % storage % U_x, mesh % faces(fID) % elementSide(1),1)                        
            end do
            !$acc end parallel loop
            !$omp end do

            !$omp do schedule(runtime)
            !$acc parallel loop gang present(mesh) private(fID)
            do zonefID = 1, mesh % zones(zoneID) % no_of_faces
               fID =  mesh % zones(zoneID) % faces(zonefID)
               eID = mesh % faces(fID) % elementIDs(1)
               call HexElement_ProlongGradientsToFaces_GL(mesh % elements(eID), NGRAD, &
                                                          mesh % faces(fID), &
                                                          mesh % elements(eID) % storage % U_y, mesh % faces(fID) % elementSide(1),2)                        
            end do
            !$acc end parallel loop
            !$omp end do
         
            !$omp do schedule(runtime)
            !$acc parallel loop gang present(mesh) private(fID)
            do zonefID = 1, mesh % zones(zoneID) % no_of_faces
               fID =  mesh % zones(zoneID) % faces(zonefID)
               eID = mesh % faces(fID) % elementIDs(1)
               call HexElement_ProlongGradientsToFaces_GL(mesh % elements(eID), NGRAD, &
                                                          mesh % faces(fID), &
                                                          mesh % elements(eID) % storage % U_z, mesh % faces(fID) % elementSide(1),3)                        
            end do
            !$acc end parallel loop
            !$omp end do
            end select
         end if
!
!        Loop the zone to get faces and elements
!        ---------------------------------------
         !$acc wait

   end subroutine ProlongToFaces

end module SurfaceIntegrals
#endif
