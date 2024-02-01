!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!   HORSES3D - convertSolution
!
!      This module converts result of one particular .mesh file to another .mesh file with identical geometry
!
!/////////////////////////////////////////////////////////////////////////////////////////////////////////
!
#include "Includes.h"
MODULE convertSolution
    USE SMConstants
    USE InterpolationMatrices
    USE SharedSpectralBasis
    USE Storage
    USE NodalStorageClass
	use SolutionFile
	use PhysicsStorage
	use TransfiniteMapClass
	use ElementConnectivityDefinitions, only : NODES_PER_ELEMENT
	
	public ProjectStoragePoints
   
!
!     ========
      CONTAINS
!     ========
!
!
!////////////////////////////////////////////////////////////////////////
!
        SUBROUTINE convertSolutionToDiffMesh (meshFile1, boundaryFile1, resultFile1, meshFile2, boundaryFile2, polyOrder)
            IMPLICIT NONE
            CHARACTER(LEN=LINE_LENGTH), INTENT(IN)     :: meshFile1, boundaryFile1, resultFile1
			CHARACTER(LEN=LINE_LENGTH), INTENT(IN)     :: meshFile2, boundaryFile2
			INTEGER                   , INTENT(IN)     :: polyOrder(3)
!
!        ---------------
!        Local variables
!        ---------------
!
            type(Mesh_t)                               :: mesh
			type(Mesh_t)							   :: mesh2
            integer                                    :: eID
            real(kind=RP)                              :: xi(0:polyOrder(1)), eta(0:polyOrder(2)), zeta(0:polyOrder(3))
            integer                                    :: i,fid, iSol, Nout(3)
			integer                       			   :: pos, pos2
			character(len=LINE_LENGTH) 				   :: dir, time
!
!  Write Header Log
!  ------------------------	
			write(STD_OUT,'(A,A)') "/*-------------------------------------- HORSES3D - FOAM FILE --------------------------------------*\"
			write(STD_OUT,'(A,A)') "####################################################################################################"
			write(STD_OUT,'(A,A)') "#                                                                                                  #"
			write(STD_OUT,'(A,A)') "#            HORSES3D High-Order (DG) Spectral Element Sequential Navier-Stokes Solver             #"
			write(STD_OUT,'(A,A)') "#                             Convert .hsol Solution to Different Mesh                             #"
			write(STD_OUT,'(A,A)') "#                                                                                                  #"
			write(STD_OUT,'(A,A)') "####################################################################################################"
			write(STD_OUT,'(A,A)') "\*Require Input: Mesh Filename 1, Boundary Filename 1, Result 1                                     */"
			write(STD_OUT,'(A,A)') "\*               Mesh Filename 2, Boundary Filename 2, Polynomial Order                             */"
			write(STD_OUT,'(A,A)') "\*--------------------------------------------------------------------------------------------------*/"
!
!        	Describe Input
!        	--------------
			write(STD_OUT,'(/)')
			write(STD_OUT,'(10X,A,A)') "Input Control File:"
			write(STD_OUT,'(10X,A,A)') "-------------------"
			write(STD_OUT,'(30X,A,A25,A30)') "->","Task: ", "meshInterpolation"
			write(STD_OUT,'(30X,A,A25,A30)') "->","Mesh Filename 1: ", trim(meshFile1)
			write(STD_OUT,'(30X,A,A25,A30)') "->","Boundary Filename 1: ", trim(boundaryFile1)
			write(STD_OUT,'(30X,A,A25,A30)') "->","Result 1: ", trim(resultFile1)
			write(STD_OUT,'(30X,A,A25,A30)') "->","Mesh Filename 2: ", trim(meshFile2)
			write(STD_OUT,'(30X,A,A25,A30)') "->","Boundary Filename 2: ", trim(boundaryFile2)
			write(STD_OUT,'(30X,A,A25,I5,I5,I5)') "->","Polynomial Result 2: ", polyOrder(1), polyOrder(2), polyOrder(3)
!
!        	Read the meshes and solution data
!        	---------------------------------
            call mesh % ReadMesh(meshFile1)
			call mesh2 % ReadMesh(meshFile2)			
   
			call mesh % ReadSolution(resultFile1)

			mesh2 % refs = mesh % refs
			mesh2 % time = mesh % time

			
			Nout=polyOrder
			
			 call addNewSpectralBasis(spA, Nout, mesh2 % nodeType)
			 xi   = spA(Nout(1))% x
			 eta  = spA(Nout(2)) % x
			 zeta = spA(Nout(3)) % x
!
!        	Write each element zone
!        	-----------------------
			do eID = 1, mesh2 % no_of_elements
				associate ( e => mesh2 % elements(eID) )
				e % Nout = Nout
!
!           	Construct spectral basis mesh2
!           	------------------------------
				call addNewSpectralBasis(spA, e % Nmesh, mesh % nodeType)
!
!           	Construct interpolation matrices for the mesh2
!           	----------------------------------------------
				call addNewInterpolationMatrix(Tset, e % Nmesh(1), spA(e % Nmesh(1)), e % Nout(1), xi)
				call addNewInterpolationMatrix(Tset, e % Nmesh(2), spA(e % Nmesh(2)), e % Nout(2), eta)
				call addNewInterpolationMatrix(Tset, e % Nmesh(3), spA(e % Nmesh(3)), e % Nout(3), zeta)
!
!           	Perform mesh2 interpolation
!           	---------------------------
				call ProjectStoragePoints(e, Tset(e % Nout(1), e % Nmesh(1)) % T, &
														Tset(e % Nout(2), e % Nmesh(2)) % T, &
														Tset(e % Nout(3), e % Nmesh(3)) % T)
				end associate
			end do
!
!        	Perform Interpolation
!        	---------------------
			call interpolation (mesh, mesh2)
!
!        	Write Solution of Mesh 2
!        	------------------------		 
			call saveSolution(mesh2, 0, mesh2 % time, "Result_interpolation.hsol", .false.)

			 write(STD_OUT,'(/)')
			 write(STD_OUT,'(10X,A,A)') "Finish - meshInterpolation"
			 write(STD_OUT,'(10X,A,A)') "--------------------------"
			
        END SUBROUTINE convertSolutionToDiffMesh
		
	  SUBROUTINE interpolation (mesh1, mesh2)
         implicit none
         class(Mesh_t),      intent(in)  :: mesh1
		 class(Mesh_t),      intent(inout)  :: mesh2 
	  
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                       :: i, j, k, l, m, n, counter
		 integer					   :: eID1, eID2, eIDOut, eIDf, eID
		 real(kind=RP)  			   :: xi(NDIM), rhoEClosest, rhoClosest
		 logical                       :: success = .false.
		 logical                       :: success2 = .false.
		 
         real(kind=RP) , allocatable   :: lxi(:), leta(:), lzeta(:)
         type(NodalStorage_t), pointer :: spAxi, spAeta, spAzeta
		 
		 eIDOut=1
		 counter=1

		 call addNewSpectralBasis(spA, mesh1 % elements(1) % Nmesh, mesh1 % nodeType)
!
!        Interpolation Process
!        ---------------------
		 write(STD_OUT,'(/)')
         write(STD_OUT,'(10X,A,A)') "Interpolate Result:"
         write(STD_OUT,'(10X,A,A)') "------------------"

	    rhoEClosest=1.0_RP
        rhoClosest =1.0_RP
!$omp parallel shared(mesh1, mesh2, counter, eIDOut)
!$omp do schedule(runtime) private(i,j,k,l,m,n,eID1,xi,lxi,leta,lzeta, success, eIDf, eID,rhoEClosest)
		 
		 do eID2=1, mesh2 % no_of_elements
			eIDf = eIDOut
			associate( e => mesh2 % elements(eID2) )

				if (eID2.eq.counter*int(mesh2 % no_of_elements/10))then
					write(STD_OUT,'(30X,A,A15,I6,A4,I7)') "->","Element: ", eID2," of ",  mesh2 % no_of_elements
					counter=counter+1
				end if 
						  
!
!              	 Allocate memory for the coordinates
!              	 -----------------------------------            
				 allocate( e % Qout(1:NVARS,0:e % Nout(1),0:e % Nout(2),0:e % Nout(3)) )
				 
				 e % Qout (:,:,:,:)=0.0_RP
			do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
				  success = .false.
				  
!    			  Search in high order element 
					do eID=1, mesh1 % no_of_elements
						eID1=eIDf-1+INT((-1)**(eID+1)*CEILING(REAL(eID/2_RP)))
						if (eID1.le.0) eID1=eID1+mesh1 % no_of_elements
						if (eID1.gt.mesh1 % no_of_elements) eID1=eID1-mesh1 % no_of_elements
						success = FindPointWithCoords(mesh1 % elements(eID1), e % xOut(:,i, j, k), 0.01_RP, xi)
						if (success) then 
							eIDf=eID1
							exit
						end if 
					 end do 
					 
					 if (.not.success) then
						write(STD_OUT,'(10X,A,I6)') "WARNING-Element Outside Range, eID ", eID2
						write(STD_OUT,'(10X,A )') "CHECK-Meshes geometry, must be similar - Increasing Tolerance"
						success2=.false.
						do eID=1, mesh1 % no_of_elements
							eID1=eIDf-1+INT((-1)**(eID+1)*CEILING(REAL(eID/2_RP)))
							if (eID1.le.0) eID1=eID1+mesh1 % no_of_elements
							if (eID1.gt.mesh1 % no_of_elements) eID1=eID1-mesh1 % no_of_elements
							success2 = FindPointWithCoords(mesh1 % elements(eID1), e % xOut(:,i, j, k), 0.2_RP, xi)
							if (success2) then 
								exit
							end if 
						end do 
						if (.not.success2) then
							write(STD_OUT,'(20X,A,I6)') "WARNING-Element Outside Range - assigned zero velocity (wall), eID ", eID2
							e % Qout(1,i,j,k)    = rhoClosest
                            e % Qout(2:4,i,j,k)  = 0.0_RP
                            e % Qout(5,i,j,k)    = rhoEClosest
							CYCLE 
						end if 
!						CALL EXIT(0)
					 end if 
					 success=success2					  

!
!        		Get the Lagrange interpolants
!        		-----------------------------
					  associate(e1 => mesh1 % elements(eID1))
					  call addNewSpectralBasis(spA, e1 % Nmesh, mesh1 % nodeType)
					  associate( spAxi   => spA(e1 % Nmesh(1)), &
								 spAeta  => spA(e1 % Nmesh(2)), &
								 spAzeta => spA(e1 % Nmesh(3)) )
								 
					  if (allocated(lxi)) deallocate(lxi  )
					  if (allocated(leta)) deallocate(leta  )
					  if (allocated(lzeta)) deallocate(lzeta  )
					  
					  allocate( lxi(0 : e1 % Nmesh(1)) )
					  allocate( leta(0 : e1 % Nmesh(2)) )
					  allocate( lzeta(0 : e1 % Nmesh(3)) )
					  lxi = spAxi % lj(xi(1))
					  leta = spAeta % lj(xi(2))
					  lzeta = spAzeta % lj(xi(3))
	  
				  do n = 0, e1 % Nmesh(3)    ; do m = 0, e1 % Nmesh(2)  ; do l = 0, e1 % Nmesh(1)
						e % Qout(:,i,j,k) = e % Qout(:,i,j,k) + e1 % Q (:,l,m,n) * lxi(l) * leta(m) *  lzeta(n)
				  end do               ; end do             ; end do
				  
                  rhoClosest = e % Qout(1,i,j,k)
                  rhoEClosest = e % Qout(5,i,j,k)	 
				  
				  deallocate( lxi, leta, lzeta)
				  end associate
				  end associate
            end do         ; end do         ; end do
			end associate
			eIDOut = eIDf
		 end do
!$omp end do
!$omp end parallel			 
	  
	  END SUBROUTINE interpolation
!
!////////////////////////////////////////////////////////////////////////
!
      logical function FindPointWithCoords(self, x, INSIDE_TOL, xi)
!
!        **********************************************************
!
!           This function finds whether a point is inside or not
!           of the element. This is done solving
!           the mapping non-linear system
!
!        **********************************************************
!
!
         use Utilities, only: SolveThreeEquationLinearSystem
         implicit none
         class(Element_t),      intent(in)  :: self
         real(kind=RP)	,       intent(in)  :: x(NDIM)
		 real(kind=RP)  ,       intent(in)  :: INSIDE_TOL												   
         real(kind=RP)  ,       intent(out) :: xi(NDIM)
!
!        ----------------------------------
!        Newton iterative solver parameters
!        ----------------------------------
!
         integer,       parameter   :: N_MAX_ITER = 50
         real(kind=RP), parameter   :: TOL = 1.0e-6_RP
         integer,       parameter   :: STEP = 1.0_RP
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                       :: i, j, k, iter
         real(kind=RP)                 :: lxi   (0:self % Nmesh(1))
         real(kind=RP)                 :: leta  (0:self % Nmesh(2))
         real(kind=RP)                 :: lzeta (0:self % Nmesh(3))
         real(kind=RP)                 :: dlxi   (0:self % Nmesh(1))
         real(kind=RP)                 :: dleta  (0:self % Nmesh(2))
         real(kind=RP)                 :: dlzeta (0:self % Nmesh(3))
         real(kind=RP)                 :: F(NDIM)
         real(kind=RP)                 :: Jac(NDIM,NDIM)
         real(kind=RP)                 :: dx(NDIM)
         type(NodalStorage_t), Pointer :: spAxi, spAeta, spAzeta

         associate( spAxi   => spA(self % Nmesh(1)), &
         spAeta  => spA(self % Nmesh(2)), &
         spAzeta => spA(self % Nmesh(3)))
		

!
!        Initial seed
!        ------------
         xi = 0.0_RP

         do iter = 1 , N_MAX_ITER
!
!           Get Lagrange polynomials and derivatives
!           ----------------------------------------
            lxi     = spAxi % lj   (xi(1))
            leta    = spAeta % lj  (xi(2))
            lzeta   = spAzeta % lj (xi(3))

            F = 0.0_RP
            do k = 0, spAzeta % N   ; do j = 0, spAeta % N ; do i = 0, spAxi % N
               F = F + self % x(:,i,j,k) * lxi(i) * leta(j) * lzeta(k)
            end do               ; end do             ; end do


            F = F - x
!
!           Stopping criteria: there are several
!           ------------------------------------
            if ( maxval(abs(F)) .lt. TOL ) exit            
			if ( abs(xi(1)) .ge. 2.5_RP ) exit
            if ( abs(xi(3)) .ge. 2.5_RP ) exit
            if ( abs(xi(2)) .ge. 2.5_RP ) exit
			
!
!           Perform a step
!           --------------
            dlxi    = spAxi % dlj  (xi(1))
            dleta   = spAeta % dlj (xi(2))
            dlzeta  = spAzeta % dlj(xi(3))

            Jac = 0.0_RP
            do k = 0, spAzeta % N   ; do j = 0, spAeta % N ; do i = 0, spAxi % N
               Jac(:,1) = Jac(:,1) + self % x(:,i,j,k) * dlxi(i) * leta(j) * lzeta(k)
               Jac(:,2) = Jac(:,2) + self % x(:,i,j,k) * lxi(i) * dleta(j) * lzeta(k)
               Jac(:,3) = Jac(:,3) + self % x(:,i,j,k) * lxi(i) * leta(j) * dlzeta(k)
            end do               ; end do             ; end do

            dx = solveThreeEquationLinearSystem( Jac , -F )
            xi = xi + STEP * dx

         end do

!         nullify (spAxi, spAeta, spAzeta)
		 
		 end associate

         if ( (abs(xi(1)) .lt. 1.0_RP + INSIDE_TOL) .and. &
              (abs(xi(2)) .lt. 1.0_RP + INSIDE_TOL) .and. &
              (abs(xi(3)) .lt. 1.0_RP + INSIDE_TOL)          ) then
!
!           Solution is valid
!           -----------------
            FindPointWithCoords = .true.

         else
!
!           Solution is not valid
!           ---------------------
            FindPointWithCoords = .false.

         end if

      end function FindPointWithCoords
!
!////////////////////////////////////////////////////////////////////////
!
         subroutine ProjectStoragePoints(e, TxMesh, TyMesh, TzMesh)
         use NodalStorageClass
         implicit none
         type(Element_t)     :: e
         real(kind=RP),       intent(in)  :: TxMesh(0:e % Nout(1), 0:e % Nmesh(1))
         real(kind=RP),       intent(in)  :: TyMesh(0:e % Nout(2), 0:e % Nmesh(2))
         real(kind=RP),       intent(in)  :: TzMesh(0:e % Nout(3), 0:e % Nmesh(3))
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: i, j, k, l, m, n
!
!        Project mesh
!        ------------
         allocate( e % xOut(1:3,0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)) )
         e % xOut = 0.0_RP

         do n = 0, e % Nmesh(3) ; do m = 0, e % Nmesh(2) ; do l = 0, e % Nmesh(1)
            do k = 0, e % Nout(3) ; do j = 0, e % Nout(2) ; do i = 0, e % Nout(1)
               e % xOut(:,i,j,k) = e % xOut(:,i,j,k) + e % x(:,l,m,n) * TxMesh(i,l) * TyMesh(j,m) * TzMesh(k,n)
            end do            ; end do            ; end do
         end do            ; end do            ; end do

      end subroutine ProjectStoragePoints
!
!////////////////////////////////////////////////////////////////////////
!
     subroutine saveSolution(self, iter, time, name, saveGradients)
         use SolutionFile
         use MPI_Process_Info
         implicit none
         class(Mesh_t)                          :: self
         integer,             intent(in)        :: iter
         real(kind=RP),       intent(in)        :: time
         character(len=*),    intent(in)        :: name
         logical,             intent(in)        :: saveGradients
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                          :: fid, eID, padding, offsetIO
         integer(kind=AddrInt)            :: pos
         real(kind=RP), allocatable       :: Q(:,:,:,:)
         logical                          :: computeGradients = .true.
!
!        Saving file 
!        -----------		 
		 write(STD_OUT,'(/)')
         write(STD_OUT,'(10X,A,A)') "Saving Result:"
         write(STD_OUT,'(10X,A,A)') "--------------"
!
!        Create new file
!        ---------------
		 call CreateNewSolutionFile(trim(name),SOLUTION_FILE, self % nodeType, &
								   self % no_of_elements, iter, time, self % refs)
		 padding = NCONS
!
!        Write arrays
!        ------------
         fID = putSolutionFileInWriteDataMode(trim(name))
		 offsetIO = 0
         do eID = 1, self % no_of_elements
            associate( e => self % elements(eID) )

            allocate(Q(NCONS, 0:e % Nout(1), 0:e % Nout(2), 0:e % Nout(3)))
            Q(1:NCONS,:,:,:)  = e % Qout

            pos = POS_INIT_DATA + (eID-1)*5_AddrInt*SIZEOF_INT + padding* offsetIO * SIZEOF_RP
            call writeArray(fid, Q, position=pos)

            deallocate(Q)
			offsetIO=offsetIO + product(e % Nout +1)
            end associate
         end do
         close(fid)
!
!        Close the file
!        --------------
         call SealSolutionFile(trim(name))
		 write(STD_OUT,'(30X,A,A30,A,A4,I7)') "->","Result Filename: ", trim(name)

      end subroutine saveSolution


END MODULE convertSolution
