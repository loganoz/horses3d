!////////////////////////////////////////////////////////////////////////
!
!      NumericalJacobian.f90
!      Created: 2017-03-21 17:07:00 +0100 
!      By: Andrés Rueda (based on Carlos Redondo's implementation for 2D code)
!
!      Routines for computing the Jacobian matrix numerically using the colorings technique
!  ! TODO: Implement as a class with a destructor to prevent memory leaking
!////////////////////////////////////////////////////////////////////////
module NumericalJacobian
   use SMConstants
   use MatrixClass
   use ColorsClass
   use HexMeshClass
   use DGSEMClass
   use ElementClass
   use Jacobian,  only: JACEPS, local2ijk, Look_for_neighbour
   use PhysicsStorage
   use Utilities, only: Qsort
   implicit none
   
   type(Neighbour),allocatable :: nbr(:)  ! Neighbors information
   type(Colors)                :: ecolors
   
contains
   subroutine NumericalJacobian_Compute(sem, nEqn, nGradEqn, t, Matrix, ComputeTimeDerivative, PINFO, eps_in )
      use StopwatchClass
      !-------------------------------------------------------------------
      type(DGSem),                intent(inout), target  :: sem
      integer,                    intent(in)             :: nEqn, nGradEqn
      real(kind=RP),              intent(IN)             :: t
      class(Matrix_t)          ,  intent(inout)          :: Matrix
      procedure(ComputeQDot_FCN)                         :: ComputeTimeDerivative      !   
      logical,                    OPTIONAL               :: PINFO                      !<? Print information?
      real(kind=RP),              optional               :: eps_in
      !-------------------------------------------------------------------
      integer                                            :: nelm
      integer                                            :: thiscolor, thiselmidx, thiselm         ! specific counters
      integer                                            :: thisdof, elmnbr, nbrnbr                ! specific counters
      integer, allocatable, dimension(:), save           :: used                                   ! array containing index of elements whose contributions to Jacobian has already been considered
      integer                                            :: usedctr                                ! counter to fill positions of used
      integer                                            :: ielm, felm, ndof                      
      integer, save                                      :: nnz, totalnnz
      integer                           , save           :: maxndofel
      integer, allocatable, dimension(:), save           :: ndofelm, firstIdx                      ! Number of degrees of freedom and relative position in Jacobian for each element 
      integer, allocatable, dimension(:), save           :: Nx, Ny, Nz                             ! Polynomial orders
      integer, allocatable, dimension(:), save           :: ndofcol                                ! Maximum number of degrees of freedom in each color        
      integer, allocatable                               :: cols(:)
      integer, allocatable                               :: rows(:)
      integer, allocatable                               :: diag(:)
      real(kind=RP), allocatable, save                   :: Q0(:), QDot0(:)
      
      integer :: i, j ! General counters
      integer                                            :: icol
      integer, dimension(4)                              :: ijkl                                   ! Indexes to locate certain degree of freedom i,j,k...l:equation number
      integer, allocatable, dimension(:), save           :: irow_0, irow
      real(kind=RP), POINTER, dimension(:)               :: pbuffer                                ! Buffer to point to an element's Qdot
      real(kind=RP), save                                :: eps                                    ! Perturbation magnitude
      
      logical, save                                      :: isfirst = .TRUE.
#if (!defined(NAVIERSTOKES))
      logical                                            :: computeGradients = .true.
#endif
      !-------------------------------------------------------------------
      
!
!     --------------------------------------------------------------------
!     Initialize variables that will be used throughout all the simulation
!     --------------------------------------------------------------------
!
      IF (isfirst) call Stopwatch % CreateNewEvent("Numerical Jacobian construction")
      call Stopwatch % Start("Numerical Jacobian construction")
      
      IF (isfirst) THEN   
         nelm = size(sem % mesh % elements)
         
!
!        Initialize the colorings structure
!        ----------------------------------
         allocate(nbr(nelm))
         CALL Look_for_neighbour(nbr, sem % mesh)
#if (!defined(CAHNHILLIARD))
         CALL ecolors%construct(nbr,flowIsNavierStokes)
#else
         CALL ecolors%construct(nbr, .true. )
#endif
         
         allocate(ndofelm(nelm), firstIdx(nelm+1))
         allocate(Nx(nelm), Ny(nelm), Nz(nelm))
         firstIdx(1) = 1
         
         do i = 1, nelm
            Nx(i) = sem % mesh % elements(i) % Nxyz(1)
            Ny(i) = sem % mesh % elements(i) % Nxyz(2)
            Nz(i) = sem % mesh % elements(i) % Nxyz(3)

            ndofelm(i) = nEqn * (Nx(i)+1)*(Ny(i)+1)*(Nz(i)+1)
            firstIdx(i+1) = firstIdx(i) + ndofelm(i)
         end do         

         maxndofel = MAXVAL(ndofelm)                                             ! TODO: if there's p-adaptation, this value has to be recomputed
!
!        -------------------
!        Row position arrays
!        -------------------
!
         allocate(irow  (maxndofel))
         allocate(irow_0(maxndofel))
         
         irow_0(1:maxndofel) = (/ (i, i=0,maxndofel-1) /)
!
!        ---------------------------------------------------------------
!        Allocate the used array that will contain the information about
!        which neighbor elements were already used in the numerical
!        computation of the Jacobian matrix entries
!        ---------------------------------------------------------------
!
#if (!defined(CAHNHILLIARD))
         IF (flowIsNavierStokes) THEN ! .AND. BR1 (only implementation to date)
            allocate(used(26))   ! 25 neighbors (including itself) and a last entry that will be 0 always (boundary index)
         ELSE
            allocate(used(8))    ! 7 neighbors (including itself) and a last entry that will be 0 always (boundary index)
         END IF
#else
         allocate(used(26))
#endif
!
!        -------------------------------------------------------------------------
!        Set max number of nonzero values expected in a row of the Jacobian matrix    TODO: if there's p-adaptation, this has to be recomputed
!              Assumes Legendre-Gauss quadrature and neglects zero values in each 
!                 block (taken into account later when assembling)
!              For Legendre-Gauss-Lobatto much less entries are expected (a node on the
!                 interface has more cols than an interior node)
!              IMPORTANT: These numbers assume conforming meshes!
!        -------------------------------------------------------------------------
!
#if (!defined(CAHNHILLIARD))
         IF (flowIsNavierStokes) THEN ! .AND. BR1 (only implementation to date)
            nnz = maxndofel * 25
         ELSE
            nnz = maxndofel * 7
         END IF
#else
         nnz = maxndofel * 25
#endif
!
!        --------------------------------------------------------------
!        Compute the maximum number of degrees of freedom in each color               TODO: if there's p-adaptation, this has to be recomputed
!        --------------------------------------------------------------
!
         allocate(ndofcol(ecolors % ncolors))
         ndofcol = 0
         DO thiscolor = 1 , ecolors%ncolors
            ielm = ecolors%bounds(thiscolor)             
            felm = ecolors%bounds(thiscolor+1)
            DO thiselmidx = ielm, felm-1              !perturbs a dof in all elements within current color
               thiselm = ecolors%elmnts(thiselmidx)
               ndofcol(thiscolor) = MAX(ndofcol(thiscolor),ndofelm(thiselm))
            END DO
         END DO
         
         allocate(Q0(size(sem % mesh % storage % Q)))
         allocate(QDot0(size(sem % mesh % storage % QDot)))
         
         ! All initializations done!
         isfirst = .FALSE.
      END IF
!
!     ---------------------------------------------
!     Set value of eps (currently using Mettot et al. approach with L2 norm because it seems to work)
!        See:
!           > Mettot, Clément, Florent Renac, and Denis Sipp. "Computation of eigenvalue sensitivity to base flow modifications in a discrete framework: Application to open-loop control." Journal of Computational Physics 269 (2014): 234-258.
!           > Knoll, Dana A., and David E. Keyes. "Jacobian-free Newton–Krylov methods: a survey of approaches and applications." Journal of Computational Physics 193.2 (2004): 357-397.
!     --------------------------------------------
!
      if (present(eps_in)) then
         eps = eps_in
      else
         associate (Q => sem % mesh % storage % Q)
         eps = SQRT(EPSILON(eps))*(NORM2(Q)+1._RP)
         end associate
      end if
!
!     ---------------------------
!     Preallocate Jacobian matrix
!     ---------------------------
!
      select type(Matrix_p => Matrix)
         type is(DenseBlockDiagMatrix_t)
            call Matrix_p % Preallocate(nnzs=ndofelm) ! Constructing with block size
            CALL Matrix % Reset
         type is(CSRMat_t)
!~             call GetRowsAndColsVector(sem, nEqn, Matrix_p % numRows, totalnnz, firstIdx, rows, cols, diag)
!~             call Matrix_p % PreAllocateWithStructure(totalnnz, rows, cols, diag) 
            call Matrix_p % Preallocate
         class default ! Construct with nonzeros in each row
            call Matrix % Preallocate(nnz)
            CALL Matrix % Reset
      end select
      
      CALL ComputeTimeDerivative( sem % mesh, sem % particles, t, sem % BCFunctions )
!
!     Save base state in Q0 and QDot0
!     -------------------------------

#if defined(CAHNHILLIARD)
      call sem % mesh % SetStorageToEqn(2)
#endif

      Q0    = sem % mesh % storage % Q
      QDot0 = sem % mesh % storage % QDot
!
!     ------------------------------------------
!     Compute numerical Jacobian using colorings
!     ------------------------------------------
!
      DO thiscolor = 1 , ecolors%ncolors
         ielm = ecolors%bounds(thiscolor)             
         felm = ecolors%bounds(thiscolor+1)
         DO thisdof = 1, ndofcol(thiscolor)           ! Computes one column for each dof within an elment (iterates to the maximum DOF of all elements in thiscolor) 
            
            DO thiselmidx = ielm, felm-1              ! Perturbs a dof in all elements within current color
               thiselm = ecolors%elmnts(thiselmidx)
               IF (ndofelm(thiselm)<thisdof) CYCLE    ! Do nothing if the DOF exceeds the NDOF of thiselm
               
               ijkl = local2ijk(thisdof,nEqn,Nx(thiselm),Ny(thiselm),Nz(thiselm))
               
               sem%mesh%elements(thiselm)% storage % Q(ijkl(1),ijkl(2),ijkl(3),ijkl(4)) = &
                                                   sem%mesh%elements(thiselm)% storage % Q(ijkl(1),ijkl(2),ijkl(3),ijkl(4)) + eps 
            ENDDO
            
            CALL ComputeTimeDerivative( sem % mesh, sem % particles, t, sem % BCFunctions )  

            sem % mesh % storage % QDot = (sem % mesh % storage % QDot - QDot0) / eps
            
            DO thiselmidx = ielm, felm-1
               thiselm = ecolors%elmnts(thiselmidx)
               IF (ndofelm(thiselm)<thisdof) CYCLE
               ! Redifine used array and counter
               used    = 0
               usedctr = 1
               
               DO i = 1,SIZE(nbr(thiselm)%elmnt)
                  elmnbr = nbr(thiselm)%elmnt(i) 
               
                  IF (.NOT. ANY(used == elmnbr)) THEN
                     ndof   = ndofelm(elmnbr)
                     pbuffer(1:ndof) => sem%mesh%elements(elmnbr)% storage % QDot                     !maps Qdot array into a 1D pointer
                     irow = irow_0 + firstIdx(elmnbr)                                                 !generates the row indices vector
                     WHERE (ABS(pbuffer(1:ndof)) .LT. jaceps) irow = -1                               !MatSetvalues will ignore entries with irow=-1
                     icol = firstIdx(thiselm) + thisdof - 1  
                     CALL Matrix % SetColumn(ndof, irow(1:ndof), icol, pbuffer(1:ndof) )
                     
                     used(usedctr) = elmnbr
                     usedctr = usedctr + 1
                  ENDIF
                  
                  ! If we are using BR1, we also have to get the contributions of the neighbors of neighbors
#if (!defined(CAHNHILLIARD))
                  IF(flowIsNavierStokes) THEN ! .AND. BR1 (only implementation to date)
#else
                  if ( .true. ) then
#endif
                     IF (elmnbr .NE. 0) THEN
                        DO j=1, SIZE(nbr(elmnbr)%elmnt)
                           nbrnbr = nbr(elmnbr)%elmnt(j)                          
                           
                           IF (.NOT. ANY(used == nbrnbr)) THEN
                              ndof   = ndofelm(nbrnbr)
                                             
                              pbuffer(1:ndof) => sem%mesh%elements(nbrnbr)% storage % QDot       !maps Qdot array into a 1D pointer
                              irow = irow_0 + firstIdx(nbrnbr)                                   !generates the row indices vector
                              WHERE (ABS(pbuffer(1:ndof)) .LT. jaceps) irow = -1                 !SetColumn will ignore entries with irow=-1
                              icol = firstIdx(thiselm) + thisdof - 1 
                              CALL Matrix % SetColumn(ndof, irow(1:ndof), icol, pbuffer(1:ndof) )
                              used(usedctr) = nbrnbr
                              usedctr = usedctr + 1                        
                           ENDIF
                        END DO
                     END IF
                  END IF
                  
               ENDDO
            END DO           
!
!           Restore original values for Q
!           -----------------------------
            sem % mesh % storage % Q = Q0
         ENDDO
      ENDDO
      
      CALL Matrix % Assembly(firstIdx,ndofelm)                             ! Matrix A needs to be assembled before being used
      
      call Stopwatch % Pause("Numerical Jacobian construction")
      IF (PRESENT(PINFO)) THEN
         IF (PINFO) PRINT*, "Numerical Jacobian construction: ", Stopwatch % ElapsedTime("Numerical Jacobian construction"), "seconds"
      ENDIF
      call Stopwatch % Reset("Numerical Jacobian construction")
                
   END subroutine NumericalJacobian_Compute

   subroutine GetRowsAndColsVector(sem, nEqn, nRows, nnz, firstIdx, rows, cols, diag)
      implicit none
      class(DGSEM)         :: sem
      integer, intent(in)  :: nEqn, nRows
      integer, intent(in)  :: firstIdx(1:sem % mesh % no_of_elements)
      integer, intent(out)  :: nnz
      integer, allocatable, intent(out) :: rows(:)
      integer, allocatable, intent(out) :: cols(:)
      integer, allocatable, intent(out) :: diag(:)
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                 :: eID, i, j, k, counter, csr_pos, ieq, nID
      integer                 :: ii, jj, kk, iieq
      integer                 :: pos_i, pos_j
      integer                 :: lb, ub
      integer, dimension(26)  :: neighbours
  
!
!     *********************
!     First loop to get nnz: TODO could I already set rows here?
!     *********************
!
!!$omp parallel private(counter, neighbours, i, j)
!!$omp do reduction(+:nnz)
      nnz = 0
      do eID = 1, sem % mesh % no_of_elements
         associate(e => sem % mesh % elements(eID))
!
!     1/ For each element, get its neighbours
!        ------------------------------------
         counter = 0
         neighbours = -1

         do i = 1, size(nbr(eID) % elmnt)
            if ( nbr(eID) % elmnt(i) .le. 0 ) cycle
            do j = 1, size(nbr(nbr(eID) % elmnt(i)) % elmnt)
               if (nbr(nbr(eID) % elmnt(i)) % elmnt(j) .le. 0) cycle

               if ( .not. any(neighbours .eq. nbr(nbr(eID) % elmnt(i)) % elmnt(j)) ) then
                  counter = counter + 1 
                  neighbours(counter) = nbr(nbr(eID) % elmnt(i)) % elmnt(j)
               end if
            end do
         end do

         do i = 1, counter
            associate(eL => sem % mesh % elements(neighbours(i)))
            nnz = nnz + nEqn*(eL % Nxyz(1)+1)*(eL % Nxyz(2)+1)*(eL % Nxyz(3)+1)*nEqn*(e % Nxyz(1)+1)*(e % Nxyz(2)+1)*(e % Nxyz(3)+1)
            end associate
         end do
         end associate
      end do
!!$omp end do
!!$omp end parallel

      allocate(rows(1:nRows+1))
      allocate(cols(1:nnz))
      allocate(diag(1:nRows))
!
!     ****************************
!     We need to set rows and cols
!     ****************************
!
      csr_pos = 1
      do eID = 1, sem % mesh % no_of_elements
         associate(e => sem % mesh % elements(eID))
!
!     1/ For each element, get its neighbours
!        ------------------------------------
         counter = 0
         neighbours = -1

         do i = 1, size(nbr(eID) % elmnt)
            if ( nbr(eID) % elmnt(i) .le. 0 ) cycle
            do j = 1, size(nbr(nbr(eID) % elmnt(i)) % elmnt)
               if (nbr(nbr(eID) % elmnt(i)) % elmnt(j) .le. 0) cycle

               if ( .not. any(neighbours .eq. nbr(nbr(eID) % elmnt(i)) % elmnt(j)) ) then
                  counter = counter + 1 
                  neighbours(counter) = nbr(nbr(eID) % elmnt(i)) % elmnt(j)
               end if
            end do
         end do
         call Qsort(neighbours(1:counter))

         pos_i = firstIdx(eID)

         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1) ; do ieq = 1, nEqn
            rows(pos_i) = csr_pos
            pos_i = pos_i + 1 
            do nID = 1, counter
               associate(neigh => sem % mesh % elements(neighbours(nID)))

               pos_j = firstIdx(neighbours(nID))

               do kk = 0, neigh % Nxyz(3) ; do jj = 0, neigh % Nxyz(2) ; do ii = 0, neigh % Nxyz(1) ; do iieq = 1, nEqn
                    
                  cols(csr_pos) = pos_j

                  pos_j = pos_j + 1
                  csr_pos = csr_pos + 1 
               end do                  ; end do                ; end do                 ; end do
               end associate
            end do   
         end do                ; end do                ; end do                ; end do

         end associate
      end do

      rows(nRows+1) = nnz+1
!
!     **************************
!     Get the diagonal positions
!     **************************
!
      do i = 1, nRows
         lb = rows(i)
         ub = rows(i+1)-1

         pos_i = minloc(abs(cols(lb:ub)-i),dim=1)
         diag(i) = lb + pos_i - 1
      end do

   end subroutine GetRowsAndColsVector
end module NumericalJacobian
