!
! ////////////////////////////////////////////////////////////////////
! HORSES3D Converter
!     Main program horsesConverter
!
!      This module convert VTK results from OpenFOAM to horses3d .hsol (identical nodes)
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!
#include "Includes.h"
MODULE convertVTK2Horses
    USE SMConstants
    USE InterpolationMatrices
    USE SharedSpectralBasis
    USE foamCreateMeshFileConverter
    USE Storage
    USE NodalStorageClass
	use SolutionFile
	use PhysicsStorage
	use foamResultVTKStorageConverter
	use PolynomialInterpAndDerivsModule
	use convertSolution	
!
!     ========
      CONTAINS
!     ========
!
!
!////////////////////////////////////////////////////////////////////////
!
        SUBROUTINE convertOFVTK2Horses (meshName, boundaryFile, Nout,  VTKfile, Ref)
            IMPLICIT NONE
            CHARACTER(LEN=LINE_LENGTH), INTENT(IN)     :: meshName
			CHARACTER(LEN=LINE_LENGTH), INTENT(IN)     :: boundaryFile
			INTEGER, INTENT(IN)     				   :: Nout(NDIM)
			CHARACTER(LEN=LINE_LENGTH), INTENT(IN)     :: VTKfile
			REAL (KIND=RP)            , INTENT(IN)     :: Ref(4)

!
!        ---------------
!        Local variables
!        ---------------
!
            type(Mesh_t)                               :: mesh
			type(VTKResult_t)                          :: vtkResult
            integer                                    :: eID, pointID
            real(kind=RP)                              :: x(NDIM)
			real(kind=RP)                              :: xi(0:Nout(1)), eta(0:Nout(2)), zeta(0:Nout(3))
            integer                                    :: i, j, k, ii, fid, iSol, pIDstart, pIDstartGlobal, counter
			integer                       			   :: pos, pos2
			character(len=LINE_LENGTH) 				   :: dir, time
			real(kind=RP), parameter   				   :: TOL = 1.0e-4_RP
			
!
!  		Write Header Log
!  		------------------------	
			write(STD_OUT,'(A,A)') "/*-------------------------------------- HORSES3D - FOAM FILE --------------------------------------*\"
			write(STD_OUT,'(A,A)') "####################################################################################################"
			write(STD_OUT,'(A,A)') "#                                                                                                  #"
			write(STD_OUT,'(A,A)') "#            HORSES3D High-Order (DG) Spectral Element Sequential Navier-Stokes Solver             #"
			write(STD_OUT,'(A,A)') "#                            Convert Foam VTK result to HORSES3D                                   #"
			write(STD_OUT,'(A,A)') "#                                                                                                  #"
			write(STD_OUT,'(A,A)') "####################################################################################################"
			write(STD_OUT,'(A,A)') "\*--------------------------------------------------------------------------------------------------*/"
			write(STD_OUT,'(A,A)') "\*     Note: -VTK file must containts Points location, variables (rho p U) as Points Data (ASCII)   */"
			write(STD_OUT,'(A,A)') "\*           -Foam mesh (polyMesh) must originated from horsesMesh2OF convertion (Horses Mesh)      */"
			write(STD_OUT,'(A,A)') "\*           -Mesh and polynomial must be identical with horsesMesh2OF                              */"
			write(STD_OUT,'(A,A)') "\*--------------------------------------------------------------------------------------------------*/"
			write(STD_OUT,'(A,A)') "\*Require Input: Mesh Filename 1, Boundary Filename 1, Polynomial Order, VTK file                   */"
			write(STD_OUT,'(A,A)') "\*               Reynolds Number, Mach Number, Reference pressure (Pa), Reference temperature (K)   */"
			write(STD_OUT,'(A,A)') "\*--------------------------------------------------------------------------------------------------*/"
!
!        	Describe Input
!        	--------------
			write(STD_OUT,'(/)')
			write(STD_OUT,'(10X,A,A)') "Input Control File:"
			write(STD_OUT,'(10X,A,A)') "-------------------"
			write(STD_OUT,'(30X,A,A30,A30)') "->","Task: ", "OF2Horses"
			write(STD_OUT,'(30X,A,A30,A30)') "->","Mesh Filename 1: ", trim(meshName)
			write(STD_OUT,'(30X,A,A30,A30)') "->","Boundary Filename 1: ", trim(boundaryFile)
			write(STD_OUT,'(30X,A,A30,I5,I5,I5)') "->","Polynomial Mesh: ", Nout(1), Nout(2), Nout(3)	
			write(STD_OUT,'(30X,A,A30,A30)') "->","Discretization nodes: ", "Gauss-Lobatto"		
			write(STD_OUT,'(30X,A,A30,A30)') "->","VTK File: ", trim(VTKfile)
			write(STD_OUT,'(30X,A,A30,F9.1)') "->","Reynolds Number: ", Ref(1)
			write(STD_OUT,'(30X,A,A30,F6.4)') "->","Mach Number: ", Ref(2)
			write(STD_OUT,'(30X,A,A30,F8.1)') "->","Reference pressure (Pa): ", Ref(3)
			write(STD_OUT,'(30X,A,A30,F6.2)') "->","Reference temperature (K): ", Ref(4)
			
			boundaryFileName=boundaryFile
			hasBoundaries=.true.

!
!        Read the mesh and solution data
!        -------------------------------
			call VTKresult % Construct(VTKfile, Ref)
!
!        Read the mesh and solution data
!        -------------------------------
			call mesh % ReadMesh(meshName)
!
!        Create GaussLobatto Nodes-Nout order
!        ------------------------------------	
			 call addNewSpectralBasis(spA, Nout, GAUSSLOBATTO) ! must be Gauss Lobatto
			 xi   = spA(Nout(1))% x
			 eta  = spA(Nout(2)) % x
			 zeta = spA(Nout(3)) % x
!
!        Write each element zone
!        -----------------------
         do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )
			
            e % Nout = Nout
!
!           Construct spectral basis for both mesh and solution
!           ---------------------------------------------------
            call addNewSpectralBasis(spA, e % Nmesh, mesh % nodeType)
!
!           Construct interpolation matrices for the mesh
!           ---------------------------------------------
            call addNewInterpolationMatrix(Tset, e % Nmesh(1), spA(e % Nmesh(1)), e % Nout(1), xi)
            call addNewInterpolationMatrix(Tset, e % Nmesh(2), spA(e % Nmesh(2)), e % Nout(2), eta)     
            call addNewInterpolationMatrix(Tset, e % Nmesh(3), spA(e % Nmesh(3)), e % Nout(3), zeta)       
!
!           Perform interpolation
!           ---------------------

            call ProjectStoragePoints(e, Tset(e % Nout(1), e % Nmesh(1)) % T, &
                                                    Tset(e % Nout(2), e % Nmesh(2)) % T, &
                                                    Tset(e % Nout(3), e % Nmesh(3)) % T)
            end associate
         end do
			
!
!        Fill Data
!        -------------------------------
!
!        Get the solution file type
!        --------------------------
		 mesh % nodeType     = GAUSSLOBATTO 
         mesh % hasTimeDeriv = .false.
         mesh % isStatistics = .false.
         mesh % hasGradients = .false.
         mesh % hasSensor    = .false.
		 mesh % time         = 0_RP
		 
		 mesh % refs (GAMMA_REF) = 1.4_RP
		 mesh % refs (RGAS_REF) = 287.15_RP
		 mesh % refs (V_REF) = VTKresult % VRef
		 mesh % refs (RHO_REF) = VTKresult % rhoRef
		 mesh % refs (T_REF) = VTKresult % TRef
		 mesh % refs (MACH_REF) = VTKresult % Mach
!
!        Transfer Data
!        -------------
		 write(STD_OUT,'(/)') 
		 write(STD_OUT,'(10X,A,A)') "Transfering VTK Result to .hsol:"
		 write(STD_OUT,'(10X,A,A)') "-------------------------------"
         write(STD_OUT,'(30X,A,A30,ES10.3)') "->","Time: ", mesh % time
         write(STD_OUT,'(30X,A,A30,F7.3)') "->","Reference velocity: ", mesh % refs(V_REF)
         write(STD_OUT,'(30X,A,A30,F7.3)') "->","Reference density: ", mesh % refs(RHO_REF)
         write(STD_OUT,'(30X,A,A30,F7.3)') "->","Reference Temperature: ", mesh % refs(T_REF)
         write(STD_OUT,'(30X,A,A30,F7.3)') "->","Reference Mach number: ", mesh % refs(MACH_REF)

		 write(STD_OUT,'(30X,A,A30)') "->","Looking for element points: "
		 pIDstartGlobal=0
		 counter=0
!$omp parallel shared(mesh, VTKresult,pIDstartGlobal, counter)
!$omp do schedule(runtime) private(i,j,k,x,ii,pointID,pIDstart)			 
		 DO eID=1, mesh % no_of_elements
			pIDstart=pIDstartGlobal
			associate ( e => mesh % elements(eID) )
			
!$omp critical
		    counter = counter + 1
			IF (mod(counter,int(mesh % no_of_elements/10)).eq.0)then
				write(STD_OUT,'(25X,A,A,I10,A,I10,A)') "->  ","Looping Elements: ", counter," of ", mesh % no_of_elements
			END IF			
!$omp end critical
			
			allocate( e % Qout(1:5,0:e % Nout(1),0:e % Nout(2),0:e % Nout(3)) )
			e % Qout=0.0_RP
			DO k = 0, e % Nout(3) ; DO j = 0, e % Nout(2) ; DO i = 0, e % Nout(1)
				x= e % xOut (:,i,j,k)
				DO ii=1, VTKresult % nPoints
					pointID = pIDstart+INT((-1_RP)**(ii+1_RP)*CEILING(real(ii)/2_RP))
					if (pointID.le.0) pointID=pointID+VTKresult % nPoints
				    if (pointID.gt.VTKresult % nPoints) pointID=pointID-VTKresult % nPoints
					if ( maxval(abs(VTKresult % x % data(1:3,pointID)-x)).lt.TOL) then
						e % Qout(1:5,i,j,k)=VTKresult % Q(1:5,pointID)
						pIDstart=pointID
						exit
					end if
					if (ii.eq.VTKresult % nPoints) then 
						write(STD_OUT,'(10X,A)') "ERROR-Node is not in VTK file"
						write(STD_OUT,'(10X,A)') "CHECK-Mesh and polynomial used in horsesMesh2OF and OF2Horses - must identical"
						CALL EXIT(0)
				    end if 
				END DO 
			end do               ; end do                ; end do
			end associate
			pIDstartGlobal=pIDstart
		 END DO 
!$omp end do
!$omp end parallel
!
!        	Write Solution of VTK result to .hsol
!        	-------------------------------------	
			call saveSolution(mesh, 0, mesh % time, "Result_OF.hsol", .false.)		
			
			write(STD_OUT,'(/)')
			write(STD_OUT,'(10X,A,A)') "Finish - OF2Horses"
			write(STD_OUT,'(10X,A,A)') "------------------"

        END SUBROUTINE convertOFVTK2Horses

END MODULE convertVTK2Horses
