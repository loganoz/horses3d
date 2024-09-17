#include "Includes.h"
#if defined(NAVIERSTOKES)
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
   public   ScalarSurfaceIntegral, VectorSurfaceIntegral
   public   ScalarDataReconstruction, VectorDataReconstruction, IBMSurfaceIntegral, IBMVectorIntegral

   integer, parameter   :: SURFACE = 1
   integer, parameter   :: TOTAL_FORCE = 2
   integer, parameter   :: PRESSURE_FORCE = 3
   integer, parameter   :: VISCOUS_FORCE = 4
   integer, parameter   :: MASS_FLOW = 5
   integer, parameter   :: FLOW_RATE = 6
   integer, parameter   :: PRESSURE_DISTRIBUTION = 7
   integer, parameter   :: USER_DEFINED = 99
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
         class(Element), pointer    :: elements(:)
!
!        Initialization
!        --------------
         val = 0.0_RP
!
!        Loop the zone to get faces and elements
!        ---------------------------------------
         elements => mesh % elements
!$omp parallel private(fID, eID, fIDs) shared(elements,mesh,NodalStorage,zoneID,integralType,val,&
!$omp&                                          computeGradients)
!$omp single
         do zonefID = 1, mesh % zones(zoneID) % no_of_faces
            fID = mesh % zones(zoneID) % faces(zonefID)

            eID = mesh % faces(fID) % elementIDs(1)
            fIDs = mesh % elements(eID) % faceIDs

!$omp task depend(inout:elements(eID))
            call elements(eID) % ProlongSolutionToFaces(NCONS, mesh % faces(fIDs(1)),&
                                            mesh % faces(fIDs(2)),&
                                            mesh % faces(fIDs(3)),&
                                            mesh % faces(fIDs(4)),&
                                            mesh % faces(fIDs(5)),&
                                            mesh % faces(fIDs(6)) )
            if ( computeGradients ) then
               call elements(eID) % ProlongGradientsToFaces(NGRAD, mesh % faces(fIDs(1)),&
                                                mesh % faces(fIDs(2)),&
                                                mesh % faces(fIDs(3)),&
                                                mesh % faces(fIDs(4)),&
                                                mesh % faces(fIDs(5)),&
                                                mesh % faces(fIDs(6)) )
            end if
!$omp end task
         end do
!$omp end single
!
!        Loop the zone to get faces and elements
!        ---------------------------------------
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

#ifdef _HAS_MPI_
         localval = val
         call mpi_allreduce(localval, val, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

      end function ScalarSurfaceIntegral

      function ScalarSurfaceIntegral_Face(f, integralType) result(val)
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
         type(NodalStorage_t), pointer :: spAxi, spAeta
!
!        Initialization
!        --------------
         val = 0.0_RP
         spAxi  => NodalStorage(f % Nf(1))
         spAeta => NodalStorage(f % Nf(2))
!
!        Perform the numerical integration
!        ---------------------------------
         associate( Q => f % storage(1) % Q )
         select case ( integralType )
         case ( SURFACE )
!
!           **********************************
!           Computes the surface integral
!              val = \int dS
!           **********************************
!
            do j = 0, f % Nf(2) ;    do i = 0, f % Nf(1)
               val = val + spAxi % w(i) * spAeta % w(j) * f % geom % jacobian(i,j)
            end do          ;    end do

         case ( MASS_FLOW )
!
!           ***********************************
!           Computes the mass-flow integral
!              I = \int rho \vec{v}·\vec{n}dS
!           ***********************************
!
            do j = 0, f % Nf(2) ;    do i = 0, f % Nf(1)
!
!              Compute the integral
!              --------------------
               val = val +  (Q(IRHOU,i,j) * f % geom % normal(1,i,j)  &
                          + Q(IRHOV,i,j) * f % geom % normal(2,i,j)  &
                          + Q(IRHOW,i,j) * f % geom % normal(3,i,j) ) &
                       * spAxi % w(i) * spAeta % w(j) * f % geom % jacobian(i,j)

            end do          ;    end do

         case ( FLOW_RATE )
!
!           ***********************************
!           Computes the flow integral
!              val = \int \vec{v}·\vec{n}dS
!           ***********************************
!
            do j = 0, f % Nf(2) ;    do i = 0, f % Nf(1)
!
!              Compute the integral
!              --------------------
               val = val + (1.0_RP / Q(IRHO,i,j))*(Q(IRHOU,i,j) * f % geom % normal(1,i,j)  &
                                             + Q(IRHOV,i,j) * f % geom % normal(2,i,j)  &
                                             + Q(IRHOW,i,j) * f % geom % normal(3,i,j) ) &
                                          * spAxi % w(i) * spAeta % w(j) * f % geom % jacobian(i,j)
            end do          ;    end do

         case ( PRESSURE_FORCE )
!
!           ***********************************
!           Computes the pressure integral
!              val = \int pdS
!           ***********************************
!
            do j = 0, f % Nf(2) ;    do i = 0, f % Nf(1)
!
!              Compute the integral
!              --------------------
               p = Pressure(Q(:,i,j))
               val = val + p * spAxi % w(i) * spAeta % w(j) * f % geom % jacobian(i,j)
            end do          ;    end do


         case ( USER_DEFINED )   ! TODO
         end select
         end associate

         nullify (spAxi, spAeta)
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
         class(HexMesh),      intent(inout), target  :: mesh 
         integer,             intent(in)    :: zoneID
         integer,             intent(in)    :: integralType, iter
         real(kind=RP)                      :: val(NDIM)
         real(kind=RP)                      :: localVal(NDIM)
         real(kind=RP)                      :: valx, valy, valz
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: zonefID, fID, eID, fIDs(6), ierr
         class(Element), pointer  :: elements(:)
!
!        Initialization
!        --------------
         val = 0.0_RP
         valx = 0.0_RP
         valy = 0.0_RP
         valz = 0.0_RP
!
!        *************************
!        Perform the interpolation
!        *************************
!
         elements => mesh % elements
!$omp parallel private(fID, eID, fIDs, localVal) shared(elements,mesh,NodalStorage,zoneID,integralType,val,&
!$omp&                                        valx,valy,valz,computeGradients)
!$omp single
         do zonefID = 1, mesh % zones(zoneID) % no_of_faces
            fID = mesh % zones(zoneID) % faces(zonefID)

            eID = mesh % faces(fID) % elementIDs(1)
            fIDs = mesh % elements(eID) % faceIDs

!$omp task depend(inout:elements(eID))
            call elements(eID) % ProlongSolutionToFaces(NCONS, mesh % faces(fIDs(1)),&
                                            mesh % faces(fIDs(2)),&
                                            mesh % faces(fIDs(3)),&
                                            mesh % faces(fIDs(4)),&
                                            mesh % faces(fIDs(5)),&
                                            mesh % faces(fIDs(6)) )
            if ( computeGradients ) then
               call elements(eID) % ProlongGradientsToFaces(NGRAD, mesh % faces(fIDs(1)),&
                                                mesh % faces(fIDs(2)),&
                                                mesh % faces(fIDs(3)),&
                                                mesh % faces(fIDs(4)),&
                                                mesh % faces(fIDs(5)),&
                                                mesh % faces(fIDs(6)) )
            end if
!$omp end task
         end do
!$omp end single
!
!        Loop the zone to get faces and elements
!        ---------------------------------------
!$omp do private(fID,localVal) reduction(+:valx,valy,valz) schedule(runtime)
         do zonefID = 1, mesh % zones(zoneID) % no_of_faces
!
!           Face global ID
!           --------------
            fID = mesh % zones(zoneID) % faces(zonefID)
!
!           Compute the integral
!           --------------------
            localVal = VectorSurfaceIntegral_Face(mesh % faces(fID), integralType)
            valx = valx + localVal(1)
            valy = valy + localVal(2)
            valz = valz + localVal(3)

         end do
!$omp end do
!$omp end parallel

         val = (/valx, valy, valz/)

#ifdef _HAS_MPI_
         localVal = val
         call mpi_allreduce(localVal, val, NDIM, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

      end function VectorSurfaceIntegral

      function VectorSurfaceIntegral_Face(f, integralType) result(val)
         implicit none
         class(Face),                 intent(in)     :: f
         integer,                     intent(in)     :: integralType
         real(kind=RP)                               :: val(NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                       :: i, j      ! Face indices
         real(kind=RP)                 :: p, tau(NDIM,NDIM)
         type(NodalStorage_t), pointer :: spAxi, spAeta
!
!        Initialization
!        --------------
         val = 0.0_RP
         spAxi  => NodalStorage(f % Nf(1))
         spAeta => NodalStorage(f % Nf(2))
!
!        Perform the numerical integration
!        ---------------------------------
         associate( Q => f % storage(1) % Q, &
                  U_x => f % storage(1) % U_x, &
                  U_y => f % storage(1) % U_y, &
                  U_z => f % storage(1) % U_z   )
         select case ( integralType )
         case ( SURFACE )
!
!           **********************************
!           Computes the surface integral
!              val = \int \vec{n} dS
!           **********************************
!
            do j = 0, f % Nf(2) ;    do i = 0, f % Nf(1)
               val = val + spAxi % w(i) * spAeta % w(j) * f % geom % jacobian(i,j) &
                         * f % geom % normal(:,i,j)
            end do          ;    end do

         case ( TOTAL_FORCE )
!
!           ************************************************
!           Computes the total force experienced by the zone
!              F = \int p \vec{n}ds - \int tau'·\vec{n}ds
!           ************************************************
!
            do j = 0, f % Nf(2) ;    do i = 0, f % Nf(1)
!
!              Compute the integral
!              --------------------
               p = Pressure(Q(:,i,j))
               call getStressTensor(Q(:,i,j),U_x(:,i,j),U_y(:,i,j),U_z(:,i,j), tau)

               val = val + ( p * f % geom % normal(:,i,j) - matmul(tau,f % geom % normal(:,i,j)) ) &
                           * f % geom % jacobian(i,j) * spAxi % w(i) * spAeta % w(j)

            end do          ;    end do

         case ( PRESSURE_FORCE )
!
!           ****************************************************
!           Computes the pressure forces experienced by the zone
!              F = \int p \vec{n}ds
!           ****************************************************
!
            do j = 0, f % Nf(2) ;    do i = 0, f % Nf(1)
!
!              Compute the integral
!              --------------------
               p = Pressure(Q(:,i,j))

               val = val + ( p * f % geom % normal(:,i,j) ) * f % geom % jacobian(i,j) &
                         * spAxi % w(i) * spAeta % w(j)

            end do          ;    end do

         case ( VISCOUS_FORCE )
!
!           ************************************************
!           Computes the total force experienced by the zone
!              F =  - \int tau'·\vec{n}ds
!           ************************************************
!
            do j = 0, f % Nf(2) ;    do i = 0, f % Nf(1)
!
!              Compute the integral
!              --------------------
               call getStressTensor(Q(:,i,j),U_x(:,i,j),U_y(:,i,j),U_z(:,i,j), tau)
               val = val - matmul(tau,f % geom % normal(:,i,j)) * f % geom % jacobian(i,j) &
                           * spAxi % w(i) * spAeta % w(j)

            end do          ;    end do

         case ( USER_DEFINED )   ! TODO

         end select
         end associate
         nullify (spAxi, spAeta)
      end function VectorSurfaceIntegral_Face

!
!////////////////////////////////////////////////////////////////////////////////////////
!
!           INTEGRALS PROCEDURES FOR IBM DATA RECONSTRUCTION
!
!                          SURFACE INTEGRALS
!
!////////////////////////////////////////////////////////////////////////////////////////
   subroutine ScalarDataReconstruction( IBM, elements, STLNum, integralType, iter, autosave ) 
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
      !-local-variables-------------------------------------------------
      real(kind=rp) :: Qsurf(NCONS,IBM% NumOfInterPoints)
      real(kind=RP) :: Q(NCONS)
      integer       :: i, j, k, n
      logical       :: autosave
       
      if( .not. IBM% Integral(STLNum)% compute ) return

      call BandPointsState( IBM% BandRegion(STLNum), elements, NCONS, IBM% HO_IBM, IBM% zoneMask(STLNum), IBM% stl(STLNum) )

      if( IBM% stl(STLNum)% move ) IBM% Integral(STLNum)% ListComputed = .false.

      do i = 1, IBM% stl(STLNum)% NumOfObjs
         associate(obj => IBM% stl(STLNum)% ObjectsList(i))
         do j = 1, NumOfIntegrationVertices
            if( IBM% HO_IBM ) then 
               Q = obj% IntegrationVertices(j)% Q
            else
               if( .not. IBM% Integral(STLNum)% ListComputed ) call IBM% GetPointInterpolation( obj% IntegrationVertices(j), IBM% BandRegion(STLnum)% IBMmask ) 
               
               call GetPointState( NCONS, obj% IntegrationVertices(j), IBM% BandRegion(STLnum)% IBMmask, IBM% NumOfInterPoints, IBM% InterpolationType, Q )
               
            end if 

            obj% IntegrationVertices(j)% ScalarValue = IntegratedScalarValue( Q, obj% normal, integralType )   

         end do
         end associate
      end do

      IBM% Integral(STLNum)% ListComputed = .true.
      
      if( autosave ) call GenerateScalarmonitorTECfile( IBM, STLNum, integralType, iter )

   end subroutine ScalarDataReconstruction

   function IBMSurfaceIntegral( IBM, STLNum )
      use IBMClass 
      use MPI_Process_Info
      implicit none 

      type(IBM_type), intent(inout) :: IBM
      integer,        intent(in)    :: STLNum
      real(kind=RP)                 :: IBMSurfaceIntegral

      real(kind=RP) :: s 
#ifdef _HAS_MPI_
      real(kind=RP) :: localval
      integer       :: ierr 
#endif 
      s = IBM% stl(STLNum)% ComputeScalarIntegral()
#ifdef _HAS_MPI_
      localval = s
      call mpi_allreduce(localval, s, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
      IBMSurfaceIntegral = s

   end function IBMSurfaceIntegral
!
!////////////////////////////////////////////////////////////////////////////////////////
!
!                          VECTOR INTEGRALS
!
!////////////////////////////////////////////////////////////////////////////////////////
   subroutine VectorDataReconstruction( IBM, elements, STLNum, integralType, iter, autosave )
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
      !-local-variables---------------------------------------------------------------------------
      real(kind=rp) :: Qsurf(NCONS,IBM% NumOfInterPoints),   &
                       U_xsurf(NCONS,IBM% NumOfInterPoints), &
                       U_ysurf(NCONS,IBM% NumOfInterPoints), &
                       U_zsurf(NCONS,IBM% NumOfInterPoints)
      real(kind=RP) :: Q(NCONS), U_x(NCONS), U_y(NCONS), U_z(NCONS), sum_
      integer       :: i, j, k, n
      logical       :: found, autosave
      
      if( .not. IBM% Integral(STLNum)% compute ) return 
 
      call BandPointsState( IBM% BandRegion(STLNum), elements, NCONS, IBM% HO_IBM, IBM% zoneMask(STLNum), IBM% STL(STLNum) ) 

      if( IBM% stl(STLNum)% move ) IBM% Integral(STLNum)% ListComputed = .false.
      
      do i = 1, IBM% stl(STLNum)% NumOfObjs
         associate( obj =>  IBM% stl(STLNum)% ObjectsList(i) )
         do j = 1, NumOfIntegrationVertices
            if( IBM% HO_IBM ) then
               Q   = obj% IntegrationVertices(j)% Q 
               U_x = obj% IntegrationVertices(j)% U_x
               U_y = obj% IntegrationVertices(j)% U_y
               U_z = obj% IntegrationVertices(j)% U_z
            else
               if( .not. IBM% Integral(STLNum)% ListComputed ) call IBM% GetPointInterpolation( obj% IntegrationVertices(j), IBM% BandRegion(STLnum)% IBMmask ) !call GetSurfaceState( IBM, obj% IntegrationVertices(j), STLNum ) 
               
               call GetPointState( NCONS, obj% IntegrationVertices(j), IBM% BandRegion(STLnum)% IBMmask, IBM% NumOfInterPoints, IBM% InterpolationType, Q )

               if( any(isnan(Q))) then 
                  write(*,*) 'indeces =', obj% IntegrationVertices(j)% indeces
                  write(*,*) 'domains =', obj% IntegrationVertices(j)% domains
                  write(*,*) 'b =', obj% IntegrationVertices(j)% b 
                  error stop
               end if 

               call GetPointGrads( NCONS, obj% IntegrationVertices(j), IBM% BandRegion(STLnum)% IBMmask, IBM% NumOfInterPoints, IBM% InterpolationType, U_x, U_y, U_z )
               if( any(isnan(U_x))) then 
                  write(*,*) 'indeces =', obj% IntegrationVertices(j)% indeces
                  write(*,*) 'domains =', obj% IntegrationVertices(j)% domains
                  write(*,*) 'b =', obj% IntegrationVertices(j)% b 
                  error stop
               end if
               if( any(isnan(U_y))) then 
                  write(*,*) 'indeces =', obj% IntegrationVertices(j)% indeces
                  write(*,*) 'domains =', obj% IntegrationVertices(j)% domains
                  write(*,*) 'b =', obj% IntegrationVertices(j)% b 
                  error stop
               end if
               if( any(isnan(U_z))) then 
                  write(*,*) 'indeces =', obj% IntegrationVertices(j)% indeces
                  write(*,*) 'domains =', obj% IntegrationVertices(j)% domains
                  write(*,*) 'b =', obj% IntegrationVertices(j)% b 
                  error stop
               end if
            end if 
            
            obj% IntegrationVertices(j)% VectorValue = IntegratedVectorValue( Q, U_x, U_y, U_z, obj% normal, IBM% IP_dWall, IBM% Wallfunction, integralType )

            if( any(isnan(obj% IntegrationVertices(j)% VectorValue))) then 
               write(*,*) 'Q =', q
               write(*,*) 'U_X =', U_X
               write(*,*) 'U_y =', U_y
               write(*,*) 'U_z =', U_z
               write(*,*) 'normal =', obj% normal
               error stop
            end if

         end do 
         end associate
      end do
      
      IBM% Integral(STLNum)% ListComputed = .true.

      if( autosave ) call GenerateVectormonitorTECfile( IBM, STLNum, integralType, iter )

   end subroutine VectorDataReconstruction

   function IBMVectorIntegral( IBM, STLNum ) 
      use IBMClass
      use MPI_Process_Info
      implicit none 

      type(IBM_type), intent(inout) :: IBM
      integer,        intent(in)    :: STLNum
      real(kind=RP)                 :: IBMVectorIntegral(NDIM)

      real(kind=RP) :: F(NDIM)
#ifdef _HAS_MPI_
      real(kind=RP) :: localval(NDIM)
      integer       :: ierr 
#endif 
      F = IBM% stl(STLNum)% ComputeVectorIntegral()
#ifdef _HAS_MPI_
      localval = F
      call mpi_allreduce(localval, F, NDIM, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
      IBMVectorIntegral = F 

   end function IBMVectorIntegral


   subroutine GetSurfaceState( IBM, IntegrationVertex, STLNum )
      use TessellationTypes
      use IBMClass
      use MPI_Process_Info
      use omp_lib
#ifdef _HAS_MPI_
      use mpi
#endif
      implicit none

      class(IBM_type),   intent(inout) :: IBM 
      type(point_type),  intent(inout) :: IntegrationVertex
      integer,           intent(in)    :: STLNum

      type(point_type) :: x(IBM% NumOfInterPoints)
      integer          :: k 

      if( IBM% Integral(STLNum)% ListComputed ) return 
      
      !call IBM% minDistance( IntegrationVertex% coords, STLNum, IntegrationVertex% indeces, IntegrationVertex% domains )
      
      do k = 1, IBM% NumOfInterPoints
         x(k)% coords = IBM% bandRegion(STLNum)% IBMmask(IntegrationVertex% domains(k))% coords(IntegrationVertex% indeces(k),:) 
      end do 
      
      ! call GetMatrixInterpolationSystem( IntegrationVertex% coords, &
      !                                    x,                         &
      !                                    IntegrationVertex% invPhi, &
      !                                    IntegrationVertex% b,      &
      !                                    IBM% InterpolationType     )
                                         
   end subroutine GetSurfaceState 
!
!////////////////////////////////////////////////////////////////////////////////////////
!
!           INVERSE DISTANCE WEIGHTED INTERPOLATION PROCEDURES FOR IBM DATA RECONSTRUCTION
!
!                                   SCALAR INTERPOLATION
!
!//////////////////////////////////////////////////////////////////////////////////////// 
   function IntegratedScalarValue( Q, normal, integralType ) result( outvalue )
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
      real(kind=rp),           intent(in)    :: Q(:), normal(:)
      integer,                 intent(in)    :: integralType
      real(kind=rp)                          :: outvalue

      outvalue = 0.0_RP
      
      select case( integralType )

         case( MASS_FLOW )
               
            outvalue = (1.0_RP / Q(IRHO))*(Q(IRHOU)*normal(IX) + Q(IRHOV)*normal(IY) + Q(IRHOW)*normal(IZ))       
            
         case ( FLOW_RATE )

            outvalue = (Q(IRHOU)*normal(IX) + Q(IRHOV)*normal(IY) + Q(IRHOW)*normal(IZ)) 
               
         case( PRESSURE_DISTRIBUTION )

            outvalue = pressure( Q )

         case ( USER_DEFINED )   ! TODO  

      end select 

   end function IntegratedScalarValue
!
!////////////////////////////////////////////////////////////////////////////////////////
!
!                          VECTOR INTERPOLATION
!
!////////////////////////////////////////////////////////////////////////////////////////         
   function IntegratedVectorValue( Q, U_x, U_y, U_z, normal, y, Wallfunction, integralType ) result( outvalue )
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
      real(kind=rp), intent(in) :: Q(NCONS), U_x(NCONS), U_y(NCONS), &
                                   U_z(NCONS), normal(NDIM)
      real(kind=rp), intent(in) :: y
      logical,       intent(in) :: Wallfunction
      integer,       intent(in) :: integralType
      real(kind=rp)             :: outvalue(NDIM)
      !-local-variables-----------------------------------------------------------
      integer       :: i
      real(kind=rp) :: viscStress(NDIM), U(NDIM), U_t(NDIM), tangent(NDIM),   &
                       Qi(NCONS), U_xi(NCONS), U_yi(NCONS), U_zi(NCONS),      & 
                       tau(NDIM,NDIM), P, T, T_w, rho_w, mu, nu, u_II, u_tau, &
                       tau_w, kappa_                                        
      
      outvalue = 0.0_RP

      select case( integralType )

         case ( TOTAL_FORCE )

            P = pressure(Q)

            if( Wallfunction ) then
#if defined(NAVIERSTOKES) 
               T  = Temperature(Q)

               call get_laminar_mu_kappa(Q,mu,kappa_) 

               nu = mu/Q(IRHO)
                
               U   = Q(IRHOU:IRHOW)/Q(IRHO)
               U_t = U - ( dot_product(U,normal) * normal )
 
               tangent = U_t/norm2(U_t)

               u_II  = dot_product(U,tangent)
               
               u_tau = u_tau_f( u_II, y, nu, u_tau0=0.1_RP )
            
               T_w   = T + (dimensionless% Pr)**(1._RP/3._RP)/(2.0_RP*thermodynamics% cp) * POW2(u_II)
               T_w   = T_w * refvalues% T
               rho_w = P*refvalues% p/(thermodynamics% R * T_w)
               rho_w = rho_w/refvalues% rho
#endif
               tau_w = rho_w*POW2(u_tau)
               
               viscStress = tau_w*tangent
            else

               call getStressTensor(Q, U_x, U_y, U_z, tau)
               
               viscStress = matmul(tau,normal)
            end if
            
            outvalue = -P * normal + viscStress   
                  
         case( PRESSURE_FORCE )

            P = pressure(Q)
            
            outvalue = -P * normal
            
         case( VISCOUS_FORCE )

            if( Wallfunction ) then
#if defined(NAVIERSTOKES) 
               T  = Temperature(Q)
               call get_laminar_mu_kappa(Q,mu,kappa_) 
               nu = mu/Q(IRHO)
                
               U   = Q(IRHOU:IRHOW)/Q(IRHO)
               U_t = U - ( dot_product(U,normal) * normal )
 
               tangent = U_t/norm2(U_t)

               u_II  = dot_product(U,tangent)
               
               u_tau = u_tau_f( u_II, y, nu, u_tau0=0.1_RP )
            
               T_w   = T + (dimensionless% Pr)**(1._RP/3._RP)/(2.0_RP*thermodynamics% cp) * POW2(u_II)
               T_w   = T_w * refvalues% T
               rho_w = P*refvalues% p/(thermodynamics% R * T_w)
               rho_w = rho_w/refvalues% rho
#endif
               tau_w = rho_w*POW2(u_tau)
               
               viscStress = tau_w*tangent
            else
               
               call getStressTensor(Q, U_x, U_y, U_z, tau)
               
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
      if( MPI_Process% doMPIAction ) then 
         call sendScalarPlotRoot( IBM% stl(STLNum), STLNum )
      end if 

      if( MPI_Process% isRoot ) then 
         call recvScalarPlotRoot( IBM% stl(STLNum), STLNum, x, y, z, scalar )
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
      end if 

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
                                    vector_y(:), vector_z(:)
      character(len=LINE_LENGTH) :: FileName, FinalName
#ifdef _HAS_MPI_
      integer :: ierr 
#endif
      if( MPI_Process% doMPIAction ) then 
         call sendVectorPlotRoot( IBM% stl(STLNum), STLNum )
      end if 
      if( MPI_Process% isRoot ) then 
         call recvVectorPlotRoot(  IBM% stl(STLNum), STLNum, x, y, z, vector_x, vector_y, vector_z )
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
      end if

   end subroutine GenerateVectormonitorTECfile

end module SurfaceIntegrals
#endif