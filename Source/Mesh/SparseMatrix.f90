!
!////////////////////////////////////////////////////////////////////////
!
!      SparseMatrix.f
!      Created: 2008-06-24 10:21:11 -0400 
!      By: David Kopriva  
!
!      Implements Algorithms:
!         Algorithm 146: SparseMatrix
!         Algorithm 147: SparseMatrix:Procedures
!
!////////////////////////////////////////////////////////////////////////
!
      Module SparseMatrixClass
      USE LinkedListClass
      IMPLICIT NONE 
      
      TYPE SparseMatrix
         TYPE(LinkedList), DIMENSION(:), POINTER :: table
      END TYPE SparseMatrix
!
!     --------
!     Generics
!     --------
!
      INTERFACE Construct
         MODULE PROCEDURE ConstructMatrix
      END INTERFACE Construct
      INTERFACE Destruct
         MODULE PROCEDURE DestructMatrix
      END INTERFACE Destruct
!
!     ========
      CONTAINS 
!     ========
!
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ConstructMatrix( this, N )
         IMPLICIT NONE 
         TYPE(SparseMatrix) :: this
         INTEGER            :: N, k
         ALLOCATE( this%table(N) ) 
         DO k = 1, N 
            CALL Construct( this%table(k) )
         END DO
        
      END SUBROUTINE ConstructMatrix
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DestructMatrix( this )
         IMPLICIT NONE 
         TYPE(SparseMatrix)        :: this
         INTEGER                   :: k
         DO k = 1, SIZE(this%table) 
            CALL Destruct( this%table(k) )
         END DO
         DEALLOCATE( this%table )
      END SUBROUTINE DestructMatrix
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE AddDataForKeys( this, i, j, d )
         IMPLICIT NONE 
         TYPE(SparseMatrix) :: this
         TYPE(ListData)     :: d
         INTEGER            :: i,j
         IF ( .NOT.ContainsKeys(this,i,j) )     THEN
            d%key = j
            CALL Add( this%table(i), d )
         END IF
      END SUBROUTINE AddDataForKeys
!
!////////////////////////////////////////////////////////////////////////
!
      LOGICAL FUNCTION ContainsKeys( this, i,j ) 
      IMPLICIT NONE 
         TYPE(SparseMatrix) :: this
         INTEGER            :: i, j
         TYPE(ListData)     :: d
         ContainsKeys = .false.
         IF( .NOT.ASSOCIATED( this%table(i)%head ) ) RETURN
         this%table(i)%current => this%table(i)%head
         DO WHILE ( ASSOCIATED(this%table(i)%current) )
            CALL GetCurrentData( this%table(i), d )
            IF ( d%key == j )     THEN
               ContainsKeys = .true.
               EXIT
            END IF
            CALL MoveToNext( this%table(i) )
         END DO 
      END FUNCTION ContainsKeys
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DataForKeys( this, i, j, d )
         IMPLICIT NONE 
         TYPE(SparseMatrix) :: this
         TYPE(ListData)     :: d
         INTEGER            :: i,j
         IF( .NOT.ASSOCIATED( this%table(i)%head ) ) RETURN
         this%table(i)%current => this%table(i)%head
         DO WHILE ( ASSOCIATED(this%table(i)%current) )
            CALL GetCurrentData( this%table(i), d )
            IF ( d%key == j )     EXIT 
            CALL MoveToNext( this%table(i) )
         END DO 
      END SUBROUTINE DataForKeys
      
      END Module SparseMatrixClass
      