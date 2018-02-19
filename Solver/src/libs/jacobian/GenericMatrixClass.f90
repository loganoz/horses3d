module GenericMatrixClass
   use SMConstants
   implicit none
   
   private
   public Matrix_t
   
   type Matrix_t
      integer :: NumRows
      contains
         procedure :: construct
         procedure :: Preallocate
         procedure :: Reset
         procedure :: SetColumn
         procedure :: AddToColumn
         procedure :: Shift
         procedure :: ReShift
         procedure :: PreAssembly
         procedure :: Assembly
         procedure :: destruct
   end type Matrix_t
contains
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine construct(this,dimPrb,withMPI)
      implicit none
      !---------------------------------------------
      class(Matrix_t) :: this
      integer, intent(in)           :: dimPrb
      logical, optional, intent(in) :: WithMPI
      
      ERROR stop ' :: construct not implemented for current matrix type'
   end subroutine construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE Preallocate(this, nnz, nnzs)
      IMPLICIT NONE
      !---------------------------------------------
      CLASS(Matrix_t), INTENT(INOUT) :: this
      INTEGER, optional, intent(in)  :: nnz
      INTEGER, optional, intent(in)  :: nnzs(:)
      
      ERROR stop ' :: Preallocate not implemented for current matrix type'
   END SUBROUTINE Preallocate
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE Reset(this)
      IMPLICIT NONE
      !---------------------------------------------
      CLASS(Matrix_t),     INTENT(INOUT)     :: this
      
      ERROR stop ' :: Reset not implemented for current matrix type'
   END SUBROUTINE Reset
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE SetColumn(this,nvalues, irow, icol, values )
      IMPLICIT NONE
      !---------------------------------------------
      CLASS(Matrix_t), INTENT(INOUT)         :: this
      INTEGER, INTENT(IN)                               :: nvalues
      INTEGER, DIMENSION(:), INTENT(IN)                 :: irow
      INTEGER, INTENT(IN)                               :: icol
      real(kind=RP) , DIMENSION(:), INTENT(IN)                 :: values
      
      ERROR stop ' :: SetColumn not implemented for current matrix type'
   END SUBROUTINE SetColumn
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
! 
   SUBROUTINE AddToColumn(this,nvalues, irow, icol, values )
      IMPLICIT NONE
      CLASS(Matrix_t), INTENT(INOUT)         :: this
      INTEGER, INTENT(IN)                               :: nvalues
      INTEGER, DIMENSION(:), INTENT(IN)                 :: irow
      INTEGER, INTENT(IN)                               :: icol
      real(kind=RP) , DIMENSION(:), INTENT(IN)                 :: values
      
      ERROR stop ' :: AddToColumn not implemented for current matrix type'
   END SUBROUTINE AddToColumn
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE Shift(this, shiftval)
      IMPLICIT NONE
      !---------------------------------------------
      CLASS(Matrix_t), INTENT(INOUT)     :: this
      real(kind=RP),                     INTENT(IN)        :: shiftval
      
      ERROR stop ' :: Shift not implemented for current matrix type'
   END SUBROUTINE Shift
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------------------
!  Removes previous shift in order to insert new one 
!              (important when Jacobian is reused)
!  --------------------------------------------------
   SUBROUTINE ReShift(this, shiftval)
      IMPLICIT NONE
      !---------------------------------------------
      CLASS(Matrix_t), INTENT(INOUT)     :: this
      real(kind=RP),               INTENT(IN)        :: shiftval
      
      ERROR stop ' :: ReShift not implemented for current matrix type'
   END SUBROUTINE ReShift
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE PreAssembly(this)
      IMPLICIT NONE
      !---------------------------------------------
      CLASS(Matrix_t),     INTENT(INOUT)   :: this
      
      ERROR stop ' :: PreAssembly not implemented for current matrix type'
   END SUBROUTINE PreAssembly
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE Assembly(this,BlockIdx,BlockSize)
      IMPLICIT NONE
      !---------------------------------------------
      CLASS(Matrix_t),     INTENT(INOUT)   :: this
      INTEGER, TARGET, OPTIONAL    ,     INTENT(IN)      :: BlockIdx(:)
      INTEGER, TARGET, OPTIONAL, INTENT(IN)    :: BlockSize(:)
      !---------------------------------------------
      
      ERROR stop ' :: Assembly not implemented for current matrix type'
   END SUBROUTINE Assembly
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE destruct(this)
      IMPLICIT NONE
      !---------------------------------------------
      CLASS(Matrix_t),     INTENT(INOUT)     :: this
      
      ERROR stop ' :: destruct not implemented for current matrix type'
   END SUBROUTINE destruct
   
end module GenericMatrixClass
