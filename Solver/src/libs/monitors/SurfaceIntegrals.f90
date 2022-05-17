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
   public   SURFACE, TOTAL_FORCE, PRESSURE_FORCE, VISCOUS_FORCE, MASS_FLOW, FLOW_RATE
   public   ScalarSurfaceIntegral, VectorSurfaceIntegral

   integer, parameter   :: SURFACE = 1
   integer, parameter   :: TOTAL_FORCE = 2
   integer, parameter   :: PRESSURE_FORCE = 3
   integer, parameter   :: VISCOUS_FORCE = 4
   integer, parameter   :: MASS_FLOW = 5
   integer, parameter   :: FLOW_RATE = 6
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
      function ScalarSurfaceIntegral(mesh, zoneID, integralType) result(val)
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
         integer,             intent(in)    :: integralType
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
         
         
         if( mesh% IBM% active ) then 
            val = ScalarDataReconstruction( mesh% IBM, mesh% elements, integralType, zoneID  )
            return
         end if
         
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
      function VectorSurfaceIntegral(mesh, zoneID, integralType) result(val)
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
         integer,             intent(in)    :: integralType
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

         if( mesh% IBM% active ) then 
            val = VectorDataReconstruction( mesh% IBM, mesh% elements, integralType, zoneID )
            return
         end if
         
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
   function ScalarDataReconstruction( IBM, elements, integralType, STLNum ) result(val)
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
!           procedure. To do that the closest n points to the BP are selected using a Nearest
!           Neighbor Algorithm and their distance is stored in the vector Dist.
!           The surface integal is computed on the triangles the body is made of. The BP point
!           are computed using a linear map for each triangle. 
!        -----------------------------------------------------------------------------------------
      implicit none
      !-arguments---------------------------------------------------------------------------------
      type(IBM_type),                intent(inout) :: IBM
      type(element), dimension(:),   intent(in)    :: elements
      integer,                       intent(in)    :: integralType, STLNum
      real(kind=rp)                                :: val
      !-local-variables---------------------------------------------------------------------------
      real(kind=rp), dimension(NDIM)              :: Point
      real(kind=rp)                               :: LowerBound, IDW_Value,  &
                                                     ObjIntegral,            &
                                                     Dist, LocalVal
      real(kind=rp), dimension(:,:), allocatable  :: bpQ                                                                                                                           
      integer                                     :: i, j, k, ierr
      
      val = 0.0_RP

      allocate( bpQ(NCONS,IBM% BandRegion% NumOfObjs) )
                
      call IBM% BandPoint_state(elements, bpQ)

!$omp parallel shared(IBM,val,bpQ,integralType,STLNum,i)
!$omp do schedule(runtime) private(j,k,Point,IDW_Value,ObjIntegral,LowerBound,Dist) 
      do i = 1, size(IBM% root(STLNum)% ObjectsList)
     
         if( IBM% root(STLNum)% ObjectsList(i)% ComputeIntegrals ) then   

            ObjIntegral = 0.0_RP    
 
            if( integralType .eq. SURFACE ) then
               ObjIntegral = 1.0_RP
            else
               do j = 1, IBM% Integral(STLNum)% n_of_Q_points  
               
                  if( .not. IBM% Integral(STLNum)% ListComputed ) then
                     
                     call OBB(STLNum)% ChangeRefFrame(IBM% Integral(STLNum)% IntegObjs(i)% x(:,j), 'global', Point)
                     LowerBound  = -huge(1.0_RP)
                     
                     do k = 1, IBM% kdtree_n_of_interPoints
                        call MinimumDistancePoints( Point, IBM% rootPoints, IBM% BandRegion, Dist, LowerBound, & 
                                                    k, IBM% Integral(STLNum)% IntegObjs(i)% PointsIndex(:,j)   ) 
                        LowerBound = POW2(Dist)
                     end do 
                     
                  end if
                                                                                              
                  IDW_Value = IDWScalarValue( Point, IBM% BandRegion,                                       &
                                              IBM% Integral(STLNum)% IntegObjs(i)% PointsIndex(:,j),        &
                                              IBM% root(STLNum)% ObjectsList(i)% normal,                    &
                                              bpQ(:,IBM% Integral(STLNum)% IntegObjs(i)% PointsIndex(:,j)), &
                                              integralType                                                  ) 
                 
                  ObjIntegral = ObjIntegral + IBM% Integral(STLNum)% weights(j) * IDW_Value
               
               end do
            end if
         
!$omp critical
            val = val + IBM% Integral(STLNum)% IntegObjs(i)% Area * ObjIntegral
!$omp end critical
         end if
     
      end do  
!$omp end do   
!$omp end parallel   
   
      if( IBM% stl(STLNum)% move ) then
         IBM% Integral(STLNum)% ListComputed = .false.
      else
         IBM% Integral(STLNum)% ListComputed = .true.
      end if
      
#ifdef _HAS_MPI_
      localVal = val
      call mpi_allreduce(localVal, val, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
   
      deallocate(bpQ)
   
   end function ScalarDataReconstruction
!
!////////////////////////////////////////////////////////////////////////////////////////
!
!                          VECTOR INTEGRALS
!
!////////////////////////////////////////////////////////////////////////////////////////        
   function VectorDataReconstruction( IBM, elements, integralType, STLNum ) result(val)
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
!           procedure. To do that the closest n points to the BP are selected using a Nearest
!           Neighbor Algorithm and their distance is stored in the vector Dist.
!           The surface integal is computed on the triangles the body is made of. The BP point
!           are computed using a linear map for each triangle. 
!        -----------------------------------------------------------------------------------------
      implicit none
      !-arguments---------------------------------------------------------------------------------
      type(IBM_type),                intent(inout) :: IBM
      type(element), dimension(:),   intent(in)    :: elements
      integer,                       intent(in)    :: integralType, STLNum
      real(kind=rp), dimension(NDIM)               :: val
      !-local-variables---------------------------------------------------------------------------
      real(kind=rp), dimension(NDIM)              :: IDW_Value, ObjIntegral,  &
                                                     localVal, Point
      real(kind=rp)                               :: LowerBound, Dist
      integer                                     :: i, j, k, symPlaneIndex,  &
                                                     ierr
      real(kind=RP), dimension(:,:), allocatable  :: bpQ, bpU_x, bpU_y, bpU_z

      val = 0.0_RP
      
      allocate( bpQ(NCONS,IBM% BandRegion% NumOfObjs),   &
                bpU_x(NCONS,IBM% BandRegion% NumOfObjs), &
                bpU_y(NCONS,IBM% BandRegion% NumOfObjs), &
                bpU_z(NCONS,IBM% BandRegion% NumOfObjs)  )
      
      call IBM% BandPoint_state(elements, bpQ, bpU_x, bpU_y, bpU_z)
 
!$omp parallel shared(IBM,val,bpQ,bpU_x,bpU_y,bpU_z,integralType,STLNum,i)
!$omp do schedule(runtime) private(j,k,Point,IDW_Value,ObjIntegral,LowerBound,Dist)
      do i = 1, size(IBM% root(STLNum)% ObjectsList)

         if( IBM% root(STLNum)% ObjectsList(i)% ComputeIntegrals ) then

            ObjIntegral = 0.0_RP

            do j = 1, IBM% Integral(STLNum)% n_of_Q_points  

               if( .not. IBM% Integral(STLNum)% ListComputed ) then
               
                  call OBB(STLNum)% ChangeRefFrame(IBM% Integral(STLNum)% IntegObjs(i)% x(:,j), 'global', Point)
       
                  if( IBM% Wallfunction ) then
                     Point = Point + IBM% IP_Distance*IBM% root(STLNum)% ObjectsList(i)% normal
                  end if
                  
                  LowerBound = -huge(1.0_RP)
                                  
                  IBM% Integral(STLNum)% IntegObjs(i)% PointsIndex(:,j) = 0
             
                  do k = 1, IBM% kdtree_n_of_interPoints               
                     call MinimumDistancePoints( Point, IBM% rootPoints, IBM% BandRegion, Dist, LowerBound, &
                                                 k, IBM% Integral(STLNum)% IntegObjs(i)% PointsIndex(:,j)   )   
                     LowerBound = POW2(Dist)
                  end do

               end if
               
               IDW_Value  = IDWVectorValue( Point = Point, BandRegion = IBM% BandRegion,                          &
                                            PointsIndex = IBM% Integral(STLNum)% IntegObjs(i)% PointsIndex(:,j),  &
                                            normal = IBM% root(STLNum)% ObjectsList(i)% normal,                   &
                                            Q = bpQ(:,IBM% Integral(STLNum)% IntegObjs(i)% PointsIndex(:,j)),     &
                                            U_x = bpU_x(:,IBM% Integral(STLNum)% IntegObjs(i)% PointsIndex(:,j)), & 
                                            U_y = bpU_y(:,IBM% Integral(STLNum)% IntegObjs(i)% PointsIndex(:,j)), &
                                            U_z = bpU_z(:,IBM% Integral(STLNum)% IntegObjs(i)% PointsIndex(:,j)), &
                                            y   = IBM% IP_Distance,                                               &
                                            Wallfunction = IBM% Wallfunction, integralType = integralType         )
               
               ObjIntegral = ObjIntegral + IBM% Integral(STLNum)% weights(j) * IDW_Value
             
            end do
         
!$omp critical
            val = val + IBM% Integral(STLNum)% IntegObjs(i)% Area * ObjIntegral
!$omp end critical

         end if    
      end do  
!$omp end do   
!$omp end parallel   

      if( IBM% stl(STLNum)% move ) then
         IBM% Integral(STLNum)% ListComputed = .false.
      else
         IBM% Integral(STLNum)% ListComputed = .true.
      end if
      
#ifdef _HAS_MPI_
      localVal = val
      call mpi_allreduce(localVal, val, NDIM, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
      
      deallocate(bpQ, bpU_x, bpU_y, bpU_z)

   end function VectorDataReconstruction
!
!////////////////////////////////////////////////////////////////////////////////////////
!
!           INVERSE DISTANCE WEIGHTED INTERPOLATION PROCEDURES FOR IBM DATA RECONSTRUCTION
!
!                                   SCALAR INTERPOLATION
!
!//////////////////////////////////////////////////////////////////////////////////////// 
   function IDWScalarValue( Point, BandRegion, PointsIndex, normal, Q, integralType ) result( outvalue )
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
      integer,                 dimension(:),    intent(in) :: PointsIndex
      type(MPI_M_Points_type),                  intent(in) :: BandRegion
      real(kind=rp),           dimension(NDIM), intent(in) :: Point, normal
      real(kind=rp),           dimension(:,:),  intent(in) :: Q 
      integer,                                  intent(in) :: integralType
      real(kind=rp)                                        :: outvalue
      !-local-variables--------------------------------------------------------
      real(kind=rp), dimension(NCONS) :: idwQ
      real(kind=rp)                   :: P

      outvalue = 0.0_RP
      
      select case( integralType )
      
         case( MASS_FLOW )
         
            call GetIDW_value( Point, BandRegion, normal, Q, PointsIndex, idwQ ) 
            
            outvalue = - (1.0_RP / idwQ(IRHO))*(idwQ(IRHOU)*normal(1) + idwQ(IRHOV)*normal(2) + idwQ(IRHOW)*normal(3))       
            
         case ( FLOW_RATE )
         
            call GetIDW_value( Point, BandRegion, normal, Q, PointsIndex, idwQ ) 
            
            outvalue = - (idwQ(IRHOU)*normal(1) + idwQ(IRHOV)*normal(2) + idwQ(IRHOW)*normal(3)) 
         
                  
         case( PRESSURE_FORCE )
         
            call GetIDW_value( Point, BandRegion, normal, Q, PointsIndex, idwQ ) 
            
            P = pressure(idwQ)
            
            outvalue = - P
 
         case ( USER_DEFINED )   ! TODO

      end select 
   
   end function IDWScalarValue
!
!////////////////////////////////////////////////////////////////////////////////////////
!
!                          VECTOR INTERPOLATION
!
!////////////////////////////////////////////////////////////////////////////////////////         
   function IDWVectorValue( Point, BandRegion, PointsIndex, normal, &
                            Q, U_x, U_y, U_z, y, Wallfunction,      &
                            integralType ) result( outvalue )
      use IBMClass
      use WallFunctionBC
      use VariableConversion
      use FluidData
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
      !-arguments--------------------------------------------------------------
      integer,                 dimension(:),    intent(in) :: PointsIndex
      type(MPI_M_Points_type),                  intent(in) :: BandRegion
      real(kind=rp),           dimension(NDIM), intent(in) :: Point, normal
      real(kind=rp),           dimension(:,:),  intent(in) :: Q, U_x, U_y, U_z
      integer,                                  intent(in) :: integralType
      logical,                                  intent(in) :: Wallfunction
      real(kind=rp),                            intent(in) :: y
      real(kind=rp),           dimension(NDIM)             :: outvalue
      !-local-variables--------------------------------------------------------
      real(kind=rp), dimension(NDIM)      :: viscStress, U, U_t, tangent
      real(kind=rp), dimension(NCONS)     :: idwQ, idwU_x, idwU_y, idwU_z 
      real(kind=rp), dimension(NDIM,NDIM) :: tau
      real(kind=rp)                       :: P, T, T_w, rho_w, mu, nu,  &
                                             u_II, u_tau, tau_w                                        
    
      optional :: y
      
      outvalue = 0.0_RP
      
      select case( integralType )
      
         case ( TOTAL_FORCE )
         
            call GetIDW_value( Point, BandRegion, normal, Q, PointsIndex, idwQ )     
            
            P = pressure(idwQ)
           
            if( Wallfunction ) then
#if defined(NAVIERSTOKES) 
               T  = Temperature(idwQ)
               mu = dimensionless% mu * SutherlandsLaw(T)
               nu = mu/idwQ(IRHO)
                
               U   = idwQ(IRHOU:IRHOW)
               U_t = U - ( dot_product(U,normal) * normal )
               
               tangent = U_t/norm2(U_t)

               u_II  = dot_product(U,tangent)
               u_tau = u_tau_f( u_II, y, nu )
            
               T_w = T + (dimensionless% Pr)**(1/3)/(2.0_RP*thermodynamics% cp) * POW2(u_II)
               rho_w = P/(thermodynamics% R * T_w)
#endif
               tau_w = rho_w*POW2(u_tau)
               
               viscStress = tau_w*tangent
            else
               call GetIDW_value( Point, BandRegion, normal, U_x, PointsIndex, idwU_x )   
               call GetIDW_value( Point, BandRegion, normal, U_y, PointsIndex, idwU_y )   
               call GetIDW_value( Point, BandRegion, normal, U_z, PointsIndex, idwU_z ) 
               
               call getStressTensor(idwQ, idwU_x, idwU_y, idwU_z, tau)
               
               viscStress = matmul(tau,normal)
            end if
            
            outvalue = -P * normal + viscStress   
                  
         case( PRESSURE_FORCE )
         
            call GetIDW_value( Point, BandRegion, normal, Q, PointsIndex, idwQ )
            
            P = pressure(idwQ)
            
            outvalue = -P * normal
            
         case( VISCOUS_FORCE )
               
           if( Wallfunction ) then
#if defined(NAVIERSTOKES) 
               call GetIDW_value( Point, BandRegion, normal, Q, PointsIndex, idwQ )  
               T  = Temperature(idwQ)
               mu = dimensionless% mu * SutherlandsLaw(T)
               nu = mu/idwQ(IRHO)
                
               U   = idwQ(IRHOU:IRHOW)
               U_t = U - ( dot_product(U,normal) * normal )
               
               tangent = U_t/norm2(U_t)

               u_II  = dot_product(U,tangent)
               u_tau = u_tau_f( u_II, y, nu )
            
               T_w = T + (dimensionless% Pr)**(1/3)/(2.0_RP*thermodynamics% cp) * POW2(u_II)
               rho_w = P/(thermodynamics% R * T_w)
#endif
               tau_w = rho_w*POW2(u_tau)
               
               viscStress = tau_w*tangent
            else
               call GetIDW_value( Point, BandRegion, normal, U_x, PointsIndex, idwU_x )   
               call GetIDW_value( Point, BandRegion, normal, U_y, PointsIndex, idwU_y )   
               call GetIDW_value( Point, BandRegion, normal, U_z, PointsIndex, idwU_z )
               call getStressTensor(idwQ, idwU_x, idwU_y, idwU_z, tau)
               viscStress = matmul(tau,normal)
            end if  
            
            outvalue = viscStress
            
         case ( USER_DEFINED )   ! TODO

      end select 
      
!~       outvalue = num/den
   
   end function IDWVectorValue

end module SurfaceIntegrals
#endif
