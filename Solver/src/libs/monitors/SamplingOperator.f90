#include "Includes.h"

module SamplingOperator
   use SMConstants
   use PhysicsStorage
   use Physics
   use FaceClass
   use ElementClass
   use HexMeshClass
   use FluidData
   use VariableConversion
   use NodalStorageClass
#ifdef _HAS_MPI_
   use mpi
#endif
   implicit none
   
   private
   public   SHEAR_STRESS_TANGENT, SHEAR_STRESS_X, SHEAR_STRESS_Y, SHEAR_STRESS_Z, PRESSURE_SURF, Q1, Q2, Q3, Q4, Q5
   public   VectorSurfaceSampling

   integer, parameter 	:: SHEAR_STRESS_TANGENT = 1 
   integer, parameter 	:: SHEAR_STRESS_X = 2
   integer, parameter 	:: SHEAR_STRESS_Y = 3
   integer, parameter 	:: SHEAR_STRESS_Z = 4
   integer, parameter 	:: PRESSURE_SURF = 5
   integer, parameter 	:: Q1 = 6
   integer, parameter 	:: Q2 = 7
   integer, parameter 	:: Q3 = 8
   integer, parameter 	:: Q4 = 9
   integer, parameter 	:: Q5 = 10
   integer, parameter   :: USER_DEFINED = 99
!
!  ========
   contains
!  ========
!
!           SUBROUTINE TO GET CONSTRUCT THE DATA 
!////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine VectorSurfaceSampling(mesh, zoneID, integralType, monitorName, data_out)
		 use MonitorDefinitions 
#ifdef _HAS_MPI_
         use mpi
#endif
         implicit none
         class(HexMesh),      intent(inout), target :: mesh
         integer,             intent(in)    		:: zoneID
         integer,             intent(in)    		:: integralType
		 real(kind=RP), allocatable, intent(out)  	:: data_out(:)

!
!        ---------------
!        Local variables
!        ---------------
!
         integer  							:: zonefID, fID, eID, fIDs(6), ierr, Nx, Ny, fsID
         class(Element), pointer  			:: elements(:)
		 real(kind=RP) , allocatable        :: data_proc(:)
		 logical 				  			:: file_exists
		 character(len=STR_LEN_MONITORS) 	:: fileName
		 character(len=STR_LEN_MONITORS) 	:: monitorName
		 
		 !
!        Get the number of Order and the sampling interval
!        -------------------------------------------------
		 if (mesh % zones(zoneID) % no_of_faces .gt.0) then
			Nx = mesh % faces (mesh % zones(zoneID) % faces(1))%Nf(1)+1
			Ny = mesh % faces (mesh % zones(zoneID) % faces(1))%Nf(2)+1
		 else
			Nx=0
			Ny=0
		 end if 
		 
		 ALLOCATE (data_out(Nx*Ny*mesh % zones(zoneID) % no_of_faces), data_proc(Nx*Ny))
		 
!
!        *************************
!        Perform the interpolation
!        *************************
!
#if defined(NAVIERSTOKES) && (!(INCNS))
         elements => mesh % elements
!$omp parallel private(fID, eID, fIDs,data_proc) shared(elements,mesh,NodalStorage,zoneID,integralType,&
!$omp&                                        computeGradients,data_out,Nx)
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
!$omp do  private(fID,data_proc) schedule(runtime)
         do zonefID = 1, mesh % zones(zoneID) % no_of_faces         
!
!           Face global ID
!           --------------
            fID = mesh % zones(zoneID) % faces(zonefID)
!
!           Compute the integral
!           --------------------
			
            CALL VectorSurfaceSampling_Face(mesh % faces(fID), integralType, fID, data_proc)

			data_out((zonefID-1)*Nx*Ny+1:zonefID*Nx*Ny)=data_proc

         end do
		 
!$omp end do 
!$omp end parallel

#endif

      end subroutine VectorSurfaceSampling
!
!           SUBROUTINE TO GET CONSTRUCT THE DATA in a FACE
!////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine VectorSurfaceSampling_Face(f, integralType, fID, data_out) 
         implicit none
         class(Face),                 intent(in)     :: f
         integer,                     intent(in)     :: integralType, fID
		 real(kind=RP)				  ,intent(out)	  :: data_out(1:((f%Nf(1)+1)*(f%Nf(2)+1)))
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                       :: i, j      ! Face indices
		 integer					   :: k
         real(kind=RP)                 :: p, tau(NDIM,NDIM)
         type(NodalStorage_t), pointer :: spAxi, spAeta
		 real(kind=RP)				   :: shear(NDIM)
		 real(kind=RP)				   :: shearTangent
		 
!
!        Initialization
!        --------------
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
#if defined(NAVIERSTOKES) && (!(INCNS))
!
!           *************************************************
!           Computes the shear stress, tangent to the surface
!           *************************************************
!
		 case ( SHEAR_STRESS_TANGENT ) 
		    k=0
            do j = 0, f % Nf(2);    do i = 0, f % Nf(1)
			   k=k+1
!
!              Get the stress tensor and the tangent component
!              ----------------------------------------------
               call getStressTensor(Q(:,i,j),U_x(:,i,j),U_y(:,i,j),U_z(:,i,j), tau)
	           shear=- matmul(tau,f % geom % normal(:,i,j))* POW2(refValues % V) *refValues % rho
			   shearTangent=dot_product(shear,f % geom % t1(:,i,j))
			   data_out(k)=shearTangent
            end do          ;    end do
!
!           *********************************************
!           Computes the shear stress, in the x direction
!           *********************************************
!
		 case ( SHEAR_STRESS_X ) 
		    k=0
            do j = 0, f % Nf(2);    do i = 0, f % Nf(1)
			   k=k+1
!
!              Get the stress tensor and the tangent component
!              ----------------------------------------------
               call getStressTensor(Q(:,i,j),U_x(:,i,j),U_y(:,i,j),U_z(:,i,j), tau)
	           shear=- matmul(tau,f % geom % normal(:,i,j))* POW2(refValues % V) *refValues % rho
			   data_out(k)=shear(1)
            end do          ;    end do
!
!           *********************************************
!           Computes the shear stress, in the x direction
!           *********************************************
!
		 case ( SHEAR_STRESS_Y ) 
		    k=0
            do j = 0, f % Nf(2);    do i = 0, f % Nf(1)
			   k=k+1
!
!              Get the stress tensor
!              ---------------------
               call getStressTensor(Q(:,i,j),U_x(:,i,j),U_y(:,i,j),U_z(:,i,j), tau)
	           shear=- matmul(tau,f % geom % normal(:,i,j))* POW2(refValues % V) *refValues % rho
			   data_out(k)=shear(2)
            end do          ;    end do
!
!           *********************************************
!           Computes the shear stress, in the x direction
!           *********************************************
!
		 case ( SHEAR_STRESS_Z ) 
		    k=0
            do j = 0, f % Nf(2);    do i = 0, f % Nf(1)
			   k=k+1
!
!              Get the stress tensor
!              ---------------------
               call getStressTensor(Q(:,i,j),U_x(:,i,j),U_y(:,i,j),U_z(:,i,j), tau)
	           shear=- matmul(tau,f % geom % normal(:,i,j))* POW2(refValues % V) *refValues % rho
			   data_out(k)=shear(3)
            end do          ;    end do
!
!           *********************************************
!           Computes the shear stress, in the x direction
!           *********************************************
!
		 case ( PRESSURE_SURF ) 
		    k=0
            do j = 0, f % Nf(2);    do i = 0, f % Nf(1)
			   k=k+1
!
!              Get the pressure
!              ---------------------
			   data_out(k)=Pressure(Q(:,i,j))* POW2(refValues % V) *refValues % rho
            end do          ;    end do	
			
		 case ( Q1 ) 
		    k=0
            do j = 0, f % Nf(2);    do i = 0, f % Nf(1)
			   k=k+1
			   data_out(k)=Q(1,i,j)
            end do          ;    end do	
			
		 case ( Q2 ) 
		    k=0
            do j = 0, f % Nf(2);    do i = 0, f % Nf(1)
			   k=k+1
			   data_out(k)=Q(2,i,j)
            end do          ;    end do	
			
		 case ( Q3 ) 
		    k=0
            do j = 0, f % Nf(2);    do i = 0, f % Nf(1)
			   k=k+1
			   data_out(k)=Q(3,i,j)
            end do          ;    end do	
			
		 case ( Q4 ) 
		    k=0
            do j = 0, f % Nf(2);    do i = 0, f % Nf(1)
			   k=k+1
			   data_out(k)=Q(4,i,j)
            end do          ;    end do	
			
		 case ( Q5 ) 
		    k=0
            do j = 0, f % Nf(2);    do i = 0, f % Nf(1)
			   k=k+1
			   data_out(k)=Q(5,i,j)
            end do          ;    end do	
#endif
         end select
         end associate
         nullify (spAxi, spAeta)
      end subroutine VectorSurfaceSampling_Face

	  
end module SamplingOperator


