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

      allocate( bpQ(NCONS,BandPoints_ALL% NumOfObjs) )
                
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
                        call MinimumDistancePoints( Point, IBM% rootPoints, Dist, LowerBound, k,          &
                                                    IBM% Integral(STLNum)% IntegObjs(i)% PointsIndex(:,j) ) 
                        LowerBound = POW2(Dist)
                     end do 
                     
                  end if
                                                                                              
                  IDW_Value = IDWScalarValue( Point, IBM% Integral(STLNum)% IntegObjs(i)% PointsIndex(:,j), &
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
                                                     localVal, Point,         &
                                                     tau_w
      real(kind=rp)                               :: LowerBound, Dist
      integer                                     :: i, j, k, symPlaneIndex,  &
                                                     ierr
      real(kind=RP), dimension(:,:), allocatable  :: bpQ, bpU_x, bpU_y, bpU_z


   real(kind=rp) :: newdist
   integer :: l

      val = 0.0_RP
      
      allocate( bpQ(NCONS,BandPoints_ALL% NumOfObjs),   &
                bpU_x(NCONS,BandPoints_ALL% NumOfObjs), &
                bpU_y(NCONS,BandPoints_ALL% NumOfObjs), &
                bpU_z(NCONS,BandPoints_ALL% NumOfObjs)  )
      
      call IBM% BandPoint_state(elements, bpQ, bpU_x, bpU_y, bpU_z)
 
!$omp parallel shared(IBM,val,bpQ,bpU_x,bpU_y,bpU_z,integralType,STLNum,i)
!$omp do schedule(runtime) private(j,k,Point,IDW_Value,ObjIntegral,LowerBound,Dist,tau_w)
      do i = 1, size(IBM% root(STLNum)% ObjectsList)

         if( IBM% root(STLNum)% ObjectsList(i)% ComputeIntegrals ) then

            ObjIntegral = 0.0_RP

            do j = 1, IBM% Integral(STLNum)% n_of_Q_points  

               if( .not. IBM% Integral(STLNum)% ListComputed ) then
                          
                  call OBB(STLNum)% ChangeRefFrame(IBM% Integral(STLNum)% IntegObjs(i)% x(:,j), 'global', Point)
                  LowerBound    = -huge(1.0_RP)
                                  
                  IBM% Integral(STLNum)% IntegObjs(i)% PointsIndex(:,j) = 0
                  
                  do k = 1, IBM% kdtree_n_of_interPoints               
                     call MinimumDistancePoints( Point, IBM% rootPoints, Dist, LowerBound, k,          &
                                                 IBM% Integral(STLNum)% IntegObjs(i)% PointsIndex(:,j) )   
                     LowerBound = POW2(Dist)
                  end do

               end if
               
               !should be performed at each iteration? i.e. does the yplus change at each iteration?
!~                if( IBM% Wallfunction ) then
!~                   IBM% Integral(STLNum)% NearestPointsTurbulence(IBM% root(STLNum)% ObjectsList(i)% index)% PointsIndex(:,l) = 0 
!~                   call GetIDW_value( Point, normal, Qtmp(:,IBM% Integral(STLNum)% NearestPoints(IBM% root(STLNum)% ObjectsList(i)% index)% PointsIndex(:,l)), &
!~                                      IBM% Integral(STLNum)% NearestPoints(IBM% root(STLNum)% ObjectsList(i)% index)% PointsIndex(:,l), Q_Point                )
!~                   call GetIBM_WallShearStress( Point, Q_Point, IBM% root(STLNum)% ObjectsList(i), IBM% rootPoints, STLNum,                                 &
!~                                                IBM% Integral(STLNum)% NearestPointsTurbulence(IBM% root(STLNum)% ObjectsList(i)% index)% PointsIndex(:,l), &
!~                                                Qtmp, tau_w                                                                                                 )
!~                else
                  tau_w = 0.0_RP
!~                end if

               IDW_Value  = IDWVectorValue( Point, IBM% Integral(STLNum)% IntegObjs(i)% PointsIndex(:,j),   &
                                            IBM% root(STLNum)% ObjectsList(i)% normal,                      &
                                            bpQ(:,IBM% Integral(STLNum)% IntegObjs(i)% PointsIndex(:,j)),   &
                                            bpU_x(:,IBM% Integral(STLNum)% IntegObjs(i)% PointsIndex(:,j)), & 
                                            bpU_y(:,IBM% Integral(STLNum)% IntegObjs(i)% PointsIndex(:,j)), &
                                            bpU_z(:,IBM% Integral(STLNum)% IntegObjs(i)% PointsIndex(:,j)), &
                                            tau_w, IBM% Wallfunction, integralType                          ) 

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
   function IDWScalarValue( Point, PointsIndex, normal, Q, integralType ) result( outvalue )
      use TessellationTypes
      use MappedGeometryClass
      use KDClass
      use MPI_IBMUtilities
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
      integer,       dimension(:),    intent(in) :: PointsIndex
      real(kind=rp), dimension(NDIM), intent(in) :: Point, normal
      real(kind=rp), dimension(:,:),  intent(in) :: Q 
      integer,                        intent(in) :: integralType
      real(kind=rp)                              :: outvalue
      !-local-variables--------------------------------------------------------
      real(kind=rp), dimension(NDIM) :: DistanceNormal
      real(kind=rp)                  :: sqrd1, d1, d2, num, den, P, &
                                        tau(NDIM,NDIM), distanceSqr 
      integer                        :: i, k

      num = 0.0_RP; den = 0.0_RP
      
      select case( integralType )
      
         case( MASS_FLOW )
         
            do i = 1, size(PointsIndex)
                          
               DistanceNormal = BandPoints_ALL% x(PointsIndex(i))% coords - Point
               d2 = vDot(DistanceNormal, normal)
               
               distanceSqr = 0.0_RP
               do k = 1, NDIM
                  distanceSqr = distanceSqr + POW2(BandPoints_ALL% x(PointsIndex(i))% coords(k) - Point(k))
               end do

               sqrd1 = distanceSqr - POW2(d2)
               if( AlmostEqual(sqrd1,0.0_RP) ) sqrd1 = 0.0_RP
               d1 = sqrt(sqrd1) 
               
               num = num - (1.0_RP / Q(1,i))*(Q(2,i)*normal(1) + Q(3,i)*normal(2) + Q(4,i)*normal(3))/d1
               den = den + 1.0_RP/d1

            end do         
            
         case ( FLOW_RATE )
         
            do i = 1, size(PointsIndex)
                          
               DistanceNormal = BandPoints_ALL% x(PointsIndex(i))% coords - Point
               d2    = vDot(DistanceNormal, normal)

               distanceSqr = 0.0_RP
               do k = 1, NDIM
                  distanceSqr = distanceSqr + POW2(BandPoints_ALL% x(PointsIndex(i))% coords(k) - Point(k))
               end do

               sqrd1 = distanceSqr - POW2(d2)               
               if( AlmostEqual(sqrd1,0.0_RP) ) sqrd1 = 0.0_RP
               d1 = sqrt(sqrd1) 
               
               num = num - (Q(2,i)*normal(1) + Q(3,i)*normal(2) + Q(4,i)*normal(3))/d1
               den = den + 1.0_RP/d1

            end do         
                  
         case( PRESSURE_FORCE )
         
            do i = 1, size(PointsIndex)
     
               DistanceNormal = BandPoints_ALL% x(PointsIndex(i))% coords - Point
               d2 = vDot(DistanceNormal, normal)
               
               distanceSqr = 0.0_RP
               do k = 1, NDIM
                  distanceSqr = distanceSqr + POW2(BandPoints_ALL% x(PointsIndex(i))% coords(k) - Point(k))
               end do

               sqrd1 = distanceSqr - POW2(d2)
               if( AlmostEqual(sqrd1,0.0_RP) ) sqrd1 = 0.0_RP
               d1 = sqrt(sqrd1)
               
               P = pressure( Q(:,i) )
               
               num = num - P/d1
               den = den + 1.0_RP/d1

            end do
 
         case ( USER_DEFINED )   ! TODO

      end select 
                
      outvalue = num/den
   
   end function IDWScalarValue
!
!////////////////////////////////////////////////////////////////////////////////////////
!
!                          VECTOR INTERPOLATION
!
!////////////////////////////////////////////////////////////////////////////////////////         
   function IDWVectorValue( Point, PointsIndex, normal, Q, U_x, U_y, U_z, tau_w, Wallfunction, integralType ) result( outvalue )
      use TessellationTypes
      use MappedGeometryClass
      use KDClass
      use MPI_IBMUtilities
      use IBMClass
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
      integer,       dimension(:),    intent(in) :: PointsIndex
      real(kind=rp), dimension(NDIM), intent(in) :: Point, normal, tau_w
      real(kind=rp), dimension(:,:),  intent(in) ::  Q, U_x, U_y, U_z
      integer,                        intent(in) :: integralType
      logical,                        intent(in) :: Wallfunction
      real(kind=rp), dimension(NDIM)             :: outvalue
      !-local-variables--------------------------------------------------------
      real(kind=rp), dimension(NDIM) :: DistanceNormal, num
      real(kind=rp)                  :: sqrd1, d1, d2, den, P, tau(NDIM,NDIM), &
                                        num1, distanceSqr
      integer                        :: i, k
    
      num = 0.0_RP; den = 0.0_RP; outvalue = 0.0_RP
      
      select case( integralType )
      
         case ( TOTAL_FORCE )

            do i = 1, size(PointsIndex)        
                   
               DistanceNormal = BandPoints_ALL% x(PointsIndex(i))% coords - Point
               d2 = vDot(DistanceNormal, normal)               

               distanceSqr = 0.0_RP
               do k = 1, NDIM
                  distanceSqr = distanceSqr + POW2(BandPoints_ALL% x(PointsIndex(i))% coords(k) - Point(k))
               end do

               sqrd1 = distanceSqr - POW2(d2)
               if( AlmostEqual(sqrd1,0.0_RP) ) sqrd1 = 0.0_RP
               d1 = sqrt(sqrd1)

               P = pressure( Q(:,i) )
               
               if( Wallfunction ) then
                  num = num + (-P * normal)/d1
               else
                  call getStressTensor(Q(:,i), U_x(:,i), U_y(:,i), U_z(:,i), tau)
                  num = num + (-P * normal + matmul(tau,normal))/d1 
               end if
               
               den = den + 1.0_RP/d1

            end do     
            
            if( Wallfunction ) then
               outvalue = num/den + tau_w
            else
               outvalue = num/den
            end if
                  
         case( PRESSURE_FORCE )
         
            do i = 1, size(PointsIndex)
            
               DistanceNormal = BandPoints_ALL% x(PointsIndex(i))% coords - Point
               d2 = vDot(DistanceNormal, normal)

               distanceSqr = 0.0_RP
               do k = 1, NDIM
                  distanceSqr = distanceSqr + POW2(BandPoints_ALL% x(PointsIndex(i))% coords(k) - Point(k))
               end do

               sqrd1 = distanceSqr - POW2(d2)
               
               if( AlmostEqual(sqrd1,0.0_RP) ) sqrd1 = 0.0_RP
               d1 = sqrt(sqrd1)
               
               P = pressure( Q(:,i) )

               num = num - P*normal/d1
               den = den + 1.0_RP/d1
               
            end do
            
            outvalue = num/den
            
         case( VISCOUS_FORCE )
         
            if( Wallfunction ) then
               outvalue = tau_w
               return
            end if
         
            do i = 1, size(PointsIndex)   
                          
               DistanceNormal = BandPoints_ALL% x(PointsIndex(i))% coords - Point
               d2    = vDot(DistanceNormal, normal)
               
               distanceSqr = 0.0_RP
               do k = 1, NDIM
                  distanceSqr = distanceSqr + POW2(BandPoints_ALL% x(PointsIndex(i))% coords(k) - Point(k))
               end do

               sqrd1 = distanceSqr - POW2(d2)
               if( AlmostEqual(sqrd1,0.0_RP) ) sqrd1 = 0.0_RP
               d1 = sqrt(sqrd1)
               
               call getStressTensor(Q(:,i), U_x(:,i), U_y(:,i), U_z(:,i), tau)
               
               num = num + matmul(tau,normal)/d1
               den = den + 1.0_RP/d1

            end do    
            
            outvalue = num/den     
            
         case ( USER_DEFINED )   ! TODO

      end select 
   
   end function IDWVectorValue

end module SurfaceIntegrals
#endif
