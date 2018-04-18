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
   use Jacobian
   use PhysicsStorage
   implicit none
   
   real(kind=RP), parameter :: jaceps=1e-8_RP ! Minimum value of a Jacobian entry (smaller values are considered as 0._RP)
   
   
   type(Neighbour),allocatable :: nbr(:)  ! Neighbors information
   type(Colors)                :: ecolors
   
contains
   subroutine NumericalJacobian_Compute(sem, t, Matrix, ComputeTimeDerivative, PINFO )
      use StopwatchClass
      !-------------------------------------------------------------------
      type(DGSem),                intent(inout), target  :: sem
      real(kind=RP),              intent(IN)             :: t
      class(Matrix_t)          ,  intent(inout)          :: Matrix
      procedure(ComputeQDot_FCN)                         :: ComputeTimeDerivative      !   
      logical,                    OPTIONAL               :: PINFO                      !<? Print information?
      !-------------------------------------------------------------------
      integer                                            :: nelm
      integer                                            :: thiscolor, thiselmidx, thiselm         ! specific counters
      integer                                            :: thisdof, elmnbr, nbrnbr                ! specific counters
      integer, allocatable, dimension(:), save           :: used                                   ! array containing index of elements whose contributions to Jacobian has already been considered
      integer                                            :: usedctr                                ! counter to fill positions of used
      integer                                            :: ielm, felm, ndof                      
      integer, save                                      :: nnz
      integer                           , save           :: maxndofel
      integer, allocatable, dimension(:), save           :: ndofelm, firstIdx                      ! Number of degrees of freedom and relative position in Jacobian for each element 
      integer, allocatable, dimension(:), save           :: Nx, Ny, Nz                             ! Polynomial orders
      type(Element), allocatable        , save           :: dgs_clean(:)                           ! Clean elements array used for computations
      integer, allocatable, dimension(:), save           :: ndofcol                                ! Maximum number of degrees of freedom in each color        
      
      integer :: i, j ! General counters
      integer                                            :: icol
      integer, dimension(4)                              :: ijkl                                   ! Indexes to locate certain degree of freedom i,j,k...l:equation number
      integer, allocatable, dimension(:), save           :: irow_0, irow
      real(kind=RP), POINTER, dimension(:)               :: pbuffer                                ! Buffer to point to an element's Qdot
      real(kind=RP), save                                :: eps                                    ! Perturbation magnitude
      
      logical, save                                      :: isfirst = .TRUE.
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
#if defined(NAVIERSTOKES)
         CALL ecolors%construct(nbr,flowIsNavierStokes)
#elif defined(CAHNHILLIARD)
         CALL ecolors%construct(nbr, .true. )
#endif
         
         allocate(ndofelm(nelm), firstIdx(nelm+1))
         allocate(Nx(nelm), Ny(nelm), Nz(nelm))
         allocate(dgs_clean(nelm))
         firstIdx = 1
         DO i=1, nelm
            Nx(i) = sem%mesh%elements(i)%Nxyz(1)
            Ny(i) = sem%mesh%elements(i)%Nxyz(2)
            Nz(i) = sem%mesh%elements(i)%Nxyz(3)
!
!           --------------------------------------
!           Get block sizes and position in matrix
!           --------------------------------------
! 
            ndofelm(i)  = N_EQN * (Nx(i)+1) * (Ny(i)+1) * (Nz(i)+1)              ! TODO: if there's p-adaptation, this value has to be recomputed
            IF (i>1) firstIdx(i) = firstIdx(i-1) + ndofelm(i-1)
!
!           -------------------------------------------------------
!           Allocate the element storage of the clean element array
!           -------------------------------------------------------
!
            CALL allocateElementStorage( dgs_clean(i), N_EQN, N_GRAD_EQN, computeGradients, Nx(i), Ny(i), Nz(i) )
         END DO
         firstIdx(nelm+1) = firstIdx(nelm) + ndofelm(nelm)
         
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
#if defined(NAVIERSTOKES)
         IF (flowIsNavierStokes) THEN ! .AND. BR1 (only implementation to date)
            allocate(used(26))   ! 25 neighbors (including itself) and a last entry that will be 0 always (boundary index)
         ELSE
            allocate(used(8))    ! 7 neighbors (including itself) and a last entry that will be 0 always (boundary index)
         END IF
#elif defined(CAHNHILLIARD)
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
#if defined(NAVIERSTOKES)
         IF (flowIsNavierStokes) THEN ! .AND. BR1 (only implementation to date)
            nnz = maxndofel * 25
         ELSE
            nnz = maxndofel * 7
         END IF
#elif defined(CAHNHILLIARD)
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
         
         ! All initializarions done!
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
      associate (Q => sem % mesh % storage % Q)
      eps = SQRT(EPSILON(eps))*(NORM2(Q)+1._RP)
      end associate
!
!     ---------------------------
!     Preallocate Jacobian matrix
!     ---------------------------
!
      select type(Matrix_p => Matrix)
         type is(DenseBlockDiagMatrix_t)
            call Matrix_p % Preallocate(nnzs=ndofelm) ! Constructing with block size
         class default ! Construct with nonzeros in each row
            call Matrix % Preallocate(nnz)
      end select
      CALL Matrix % Reset
      
      CALL ComputeTimeDerivative( sem % mesh, t, sem % externalState, sem % externalGradients)
!
!     Save base state in dgs_clean
!     ----------------------------
!$omp parallel do schedule(runtime)
      do i=1, nelm
         dgs_clean(i) % storage % Q    = sem%mesh%elements(i) % storage % Q
         dgs_clean(i) % storage % Qdot = sem%mesh%elements(i) % storage % Qdot
      end do
!$omp end parallel do
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
               
               ijkl = local2ijk(thisdof,N_EQN,Nx(thiselm),Ny(thiselm),Nz(thiselm))
               
               sem%mesh%elements(thiselm)% storage % Q(ijkl(1),ijkl(2),ijkl(3),ijkl(4)) = &
                                                   sem%mesh%elements(thiselm)% storage % Q(ijkl(1),ijkl(2),ijkl(3),ijkl(4)) + eps 
            ENDDO
            
            CALL ComputeTimeDerivative( sem % mesh, t, sem % externalState, sem % externalGradients )            
            
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
                     
                     sem%mesh%elements(elmnbr)% storage % QDot = (sem%mesh%elements(elmnbr)% storage % QDot - dgs_clean(elmnbr)% storage % QDot) / eps                      
                     pbuffer(1:ndof) => sem%mesh%elements(elmnbr)% storage % QDot                     !maps Qdot array into a 1D pointer
                     irow = irow_0 + firstIdx(elmnbr)                                                 !generates the row indices vector
                     WHERE (ABS(pbuffer(1:ndof)) .LT. jaceps) irow = -1                               !MatSetvalues will ignore entries with irow=-1
                     icol = firstIdx(thiselm) + thisdof - 1  
                     CALL Matrix % SetColumn(ndof, irow(1:ndof), icol, pbuffer(1:ndof) )
                     
                     used(usedctr) = elmnbr
                     usedctr = usedctr + 1
                  ENDIF
                  
                  ! If we are using BR1, we also have to get the contributions of the neighbors of neighbors
#if defined(NAVIERSTOKES)
                  IF(flowIsNavierStokes) THEN ! .AND. BR1 (only implementation to date)
#elif defined(CAHNHILLIARD)
                  if ( .true. ) then
#endif
                     IF (elmnbr .NE. 0) THEN
                        DO j=1, SIZE(nbr(elmnbr)%elmnt)
                           nbrnbr = nbr(elmnbr)%elmnt(j)                          
                           
                           IF (.NOT. ANY(used == nbrnbr)) THEN
                              ndof   = ndofelm(nbrnbr)
                                             
                              sem%mesh%elements(nbrnbr)% storage % QDot = (sem%mesh%elements(nbrnbr)% storage % QDot - dgs_clean(nbrnbr)% storage % QDot) / eps                      
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
            DO thiselmidx = ielm, felm-1                              !Cleans modified Qs
               thiselm = ecolors%elmnts(thiselmidx)
               sem%mesh%elements(thiselm)% storage % Q = dgs_clean(thiselm)% storage % Q           
            END DO                                                
         ENDDO
      ENDDO
      
      CALL Matrix % Assembly(firstIdx,ndofelm)                             ! Matrix A needs to be assembled before being used
      
      call Stopwatch % Pause("Numerical Jacobian construction")
      IF (PRESENT(PINFO)) THEN
         IF (PINFO) PRINT*, "Numerical Jacobian construction: ", Stopwatch % ElapsedTime("Numerical Jacobian construction"), "seconds"
      ENDIF
      call Stopwatch % Reset("Numerical Jacobian construction")
                
   END subroutine NumericalJacobian_Compute
end module NumericalJacobian
