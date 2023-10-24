!
!//////////////////////////////////////////////////////
!
!      Linked list matrix
!
!////////////////////////////////////////////////////////////////////////
module LinkedListMatrixClass
   use SMConstants
   use GenericMatrixClass
   use JacobianDefinitions, only: JACEPS
#include "Includes.h"
   implicit none
   
   private
   public LinkedListMatrix_t
   
   !-----------------
   ! Type for entries
   !-----------------
   type Entry_t
      real(kind=RP)  :: value
      integer        :: col
      class(Entry_t), pointer    :: next => NULL()
   end type Entry_t
   
   !--------------
   ! Type for rows
   !--------------
   type Row_t
      type(Entry_t), pointer :: head => null()
      integer :: num_of_entries
   end type Row_t
   
   !------------------------------
   ! Main type: Linked list matrix
   !------------------------------
   type, extends(Matrix_t) :: LinkedListMatrix_t
      type(Row_t), allocatable :: rows(:)
      integer                  :: num_of_entries
      contains
      procedure :: construct
      procedure :: setEntry
      procedure :: Reset
      procedure :: ResetEntry
      procedure :: AddToEntry
      procedure :: ForceAddToEntry
      procedure :: setColumn
      procedure :: AddToColumn
      procedure :: destruct
      procedure :: PointToEntry
      procedure :: getCSRarrays
      procedure :: Shift
   end type LinkedListMatrix_t
   
contains
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------------
!  Constructor
!  -------------------------------------
   subroutine construct(this,num_of_Rows,num_of_Cols,num_of_Blocks,num_of_TotalRows,withMPI)
      implicit none
      !-arguments-----------------------------------
      class(LinkedListMatrix_t)     :: this     !<> This matrix
      integer, optional, intent(in) :: num_of_Rows
      integer, optional, intent(in) :: num_of_Cols
      integer, optional, intent(in) :: num_of_Blocks
      integer, optional, intent(in) :: num_of_TotalRows
      logical, optional, intent(in) :: WithMPI
      !-local-variables-----------------------------
      integer :: i
      !---------------------------------------------
      
      if ( .not. present(num_of_Rows) ) then
         error stop 'LinkedListMatrix_t needs num_of_Rows'
      end if
      
      this % num_of_Rows = num_of_Rows
      
      allocate ( this % rows(num_of_Rows) )
      this % num_of_entries = 0
      do i=1, num_of_Rows
         this % rows(i) % num_of_entries = 0
      end do
      
   end subroutine construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Reset(this, ForceDiagonal)
      implicit none
      !---------------------------------------------
      class(LinkedListMatrix_t),     intent(inout)     :: this
      logical, optional            , intent(in)        :: ForceDiagonal
      !-local-variables------------------------------------------------
      type(Entry_t), pointer :: CEntry, Prev
      integer                :: i
      logical                :: mustForceDiagonal
      !---------------------------------------------
      
      do i=1, this % num_of_Rows
         
         CEntry => this % rows(i) % head
         do while( associated(CEntry) )
            CEntry % value = 0._RP
            CEntry => CEntry % next
         end do
         
      end do
      
      if ( present(ForceDiagonal) ) then
         mustForceDiagonal = ForceDiagonal
      else
         mustForceDiagonal = .FALSE.
      end if
      
      if (mustForceDiagonal) then
         do i = 1, this % num_of_Rows
            call this % ForceAddToEntry(i,i,0._RP)
         end do
      end if
      
   end subroutine Reset
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine ResetEntry(this,iRow,iCol)
      implicit none
      !---------------------------------------------
      class(LinkedListMatrix_t),     intent(inout)     :: this
      integer, intent(in) :: iCol, iRow
      !-local-variables------------------------------------------------
      type(Entry_t), pointer :: CEntry, Prev
      integer                :: i
      !----------------------------------------------------------------
      
      CEntry => this % rows(iRow) % head
      
      do while( associated(CEntry) )
         if (CEntry % col == iCol) CEntry % value = 0._RP
         CEntry => CEntry % next
      end do
      
   end subroutine ResetEntry
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -------------------------------
!  Set entry of linked-list matrix
!  ------------------------------- 
   subroutine SetEntry(this, row, col, value )
      implicit none
      !-arguments-----------------------------------
      class(LinkedListMatrix_t), intent(inout) :: this
      integer        , intent(in)    :: row
      integer        , intent(in)    :: col
      real(kind=RP)  , intent(in)    :: value
      !-local-variables------------------------------
      type(Entry_t), pointer :: Entry
      !----------------------------------------------
      
      if (abs(value) < JACEPS) return
      
      if (row<1) return
      if (col<1) return
!$omp critical
      Entry => this % PointToEntry(row,col)
!$omp end critical
      Entry % value = value
      
   end subroutine SetEntry
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ---------------------------------------------
!  Add a value to an entry of linked-list matrix
!  --------------------------------------------- 
   subroutine AddToEntry(this, row, col, value )
      implicit none
      !-arguments-----------------------------------
      class(LinkedListMatrix_t), intent(inout) :: this
      integer        , intent(in)    :: row
      integer        , intent(in)    :: col
      real(kind=RP)  , intent(in)    :: value
      !-local-variables------------------------------
      type(Entry_t), pointer :: Entry
      !----------------------------------------------
      
      if (abs(value) < JACEPS) return
      
      if (row<1) return
      if (col<1) return
      
!$omp critical
      Entry => this % PointToEntry(row,col)
      Entry % value = Entry % value + value
!$omp end critical
   end subroutine AddToEntry
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ---------------------------------------------
!  Same as AddToEntry, but without checking JACEPS
!  --------------------------------------------- 
   subroutine ForceAddToEntry(this, row, col, value )
      implicit none
      !-arguments-----------------------------------
      class(LinkedListMatrix_t), intent(inout) :: this
      integer        , intent(in)    :: row
      integer        , intent(in)    :: col
      real(kind=RP)  , intent(in)    :: value
      !-local-variables------------------------------
      type(Entry_t), pointer :: Entry
      !----------------------------------------------
      
      if (row<1) return
      if (col<1) return
      
!$omp critical
      Entry => this % PointToEntry(row,col)
      Entry % value = Entry % value + value
!$omp end critical
   end subroutine ForceAddToEntry
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SetColumn(this,nvalues, irow, icol, values )
      implicit none
      !-arguments-----------------------------------
      class(LinkedListMatrix_t)  , intent(inout) :: this
      integer                    , intent(in)    :: nvalues
      integer, dimension(:)      , intent(in)    :: irow
      integer                    , intent(in)    :: icol
      real(kind=RP), dimension(:), intent(in)    :: values
      !-local-variables------------------------------
      integer :: i
      !----------------------------------------------
      
      if (nvalues .NE. size(Values)) then
         write (*,*) 'CSR_AddToCol: Dimension error (Values-RowIndexes)'
         error stop
      end if
      
      if ( icol <= 0 ) then
         write (*,*) 'CSR_AddToCol: icol error'
         error stop
      end if
      
      do i=1, nvalues
         call this % SetEntry(irow(i),icol,values(i))
      end do
      
   end subroutine SetColumn
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine AddToColumn(this,nvalues, irow, icol, values )
      implicit none
      !-arguments-----------------------------------
      class(LinkedListMatrix_t)  , intent(inout) :: this
      integer                    , intent(in)    :: nvalues
      integer, dimension(:)      , intent(in)    :: irow
      integer                    , intent(in)    :: icol
      real(kind=RP), dimension(:), intent(in)    :: values
      !-local-variables------------------------------
      integer :: i
      !----------------------------------------------
      
      if (nvalues .NE. size(Values)) then
         write (*,*) 'CSR_AddToCol: Dimension error (Values-RowIndexes)'
         error stop
      end if
      
      if ( icol <= 0 ) then
         write (*,*) 'CSR_AddToCol: icol error'
         error stop
      end if
      
      do i=1, nvalues
         call this % AddToEntry(irow(i),icol,values(i))
      end do
      
   end subroutine AddToColumn
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------
!  Destructor
!  ----------
   subroutine destruct(this)
      implicit none
      !-arguments-----------------------------------
      class(LinkedListMatrix_t), intent(inout) :: this
      !-local-variables------------------------------
      type(Entry_t), pointer :: Centry, next
      integer       :: i
      !----------------------------------------------
      
      do i=1, this % num_of_Rows
         CEntry => this % rows(i) % head
         
         do while ( associated(CEntry) )
            next => CEntry % next
            deallocate(CEntry)
            Centry => next
         end do
         this % rows(i) % head => null()
      end do
      
      deallocate ( this % rows )
      
   end subroutine destruct
   
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  -----------------------------------------------------------
!  Function to point to a specific entry of the matrix. If the
!  entry has not been created, then it creates it
!  -----------------------------------------------------------
   function PointToEntry(Matrix,i,j ) result(Entry)
      implicit none
      !-arguments------------------------------------------------------
      class(LinkedListMatrix_t), intent(inout) :: Matrix
      integer                 , intent(in)     :: i,j
      type(Entry_t)           , pointer        :: Entry
      !-local-variables------------------------------------------------
      type(Entry_t), pointer :: CEntry, Prev
      !----------------------------------------------------------------
     
      CEntry => Matrix % rows(i) % head

      nullify(Prev)
      
      do while( associated(CEntry) )
         if ( CEntry % col >= j ) exit
         Prev   => CEntry
         CEntry => CEntry % next
      end do

      if ( associated(CEntry) ) then
         if ( CEntry % col == j ) then
            Entry => CEntry
            return
         end if
      end if

      allocate( Entry )
      Entry % Value = 0._RP
      Entry % col = j
      Entry % next => CEntry
      
      if ( associated(Prev) ) then
         Prev % next => Entry
      else
         Matrix % rows(i) % head => Entry
      end if
      
      Matrix % rows(i) % num_of_entries = Matrix % rows(i) % num_of_entries + 1
      Matrix % num_of_entries = Matrix % num_of_entries + 1
      
   end function PointToEntry
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  ----------------------------------------
!  Convert linked list matrix to CSR matrix
!  ----------------------------------------
   subroutine getCSRarrays(this,Values,Cols,Rows)
      implicit none
      !-arguments------------------------------------------------------
      class(LinkedListMatrix_t), intent(inout)  :: this
      real(kind=RP), dimension(:), allocatable  :: Values
      integer      , dimension(:), allocatable  :: Cols, Rows
      !-local-variables------------------------------------------------
      integer :: i, j
      type(Entry_t), pointer :: CEntry
      !----------------------------------------------------------------
      
      safedeallocate(Rows)   ; allocate ( Rows  (this % num_of_Rows + 1) )
      safedeallocate(Cols)   ; allocate ( Cols  (this % num_of_entries) )
      safedeallocate(Values) ; allocate ( Values(this % num_of_entries) )
      
!     Set first row
!     -------------
      
      Rows(1) = 1
      CEntry => this % rows(1) % head
      do j=0, this % rows(1) % num_of_entries - 1
         Values(1+j) = CEntry % value
         Cols  (1+j) = CEntry % col
         
         CEntry => CEntry % next
      end do
      
!     Set the rest
!     ------------
      
      do i=2, this % num_of_Rows
         Rows(i) = Rows(i-1) + this % rows(i-1) % num_of_entries
         
         CEntry => this % rows(i) % head
         do j=0, this % rows(i) % num_of_entries - 1
            Values(Rows(i)+j) = CEntry % value
            Cols  (Rows(i)+j) = CEntry % col
            
            CEntry => CEntry % next
         end do
      end do
      
      Rows(this % num_of_Rows + 1) = Rows(this % num_of_Rows) + this % rows(this % num_of_Rows) % num_of_entries
      
   end subroutine getCSRarrays
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Shift(this, shiftval)
      implicit none
      !-arguments-----------------------------------
      class(LinkedListMatrix_t), intent(inout) :: this
      real(kind=RP),             intent(in)    :: shiftval
      !-local-variables-----------------------------
      type(Entry_t), pointer :: Entry
      integer :: i
      !---------------------------------------------
      
      if (abs(shiftval) < JACEPS) return
      
      do i=1, this % num_of_Rows
         Entry => this % PointToEntry(i,i)
         Entry % value = Entry % value + shiftval
      end do
      
   end subroutine Shift
end module LinkedListMatrixClass