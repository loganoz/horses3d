!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!      DenseMatUtilities.f90
!      Created: 2017-05-23 10:06:00 +0100 
!      By: Andrés Rueda
!
!      Module containing useful routines for dense matrices
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
MODULE DenseMatUtilities
   USE SMConstants
   IMPLICIT NONE

   private
   public inverse, Mat2File, ComputeLU, SolveLU
!========
 CONTAINS
!========

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine ComputeLU(A,ALU,LUpivots)
      implicit none
      !------------------------------------------------------------
      real(kind=RP), dimension(:,:)                , intent(in)  :: A          !<  Matrix to factorize
      real(kind=RP), dimension(size(A,1),size(A,1)), intent(out) :: ALU        !>  Factorized LU matrix
      integer      , dimension(size(A,1))          , intent(out) :: LUpivots   !>  LU pivots for factorization
      !------------------------------------------------------------
      integer :: n      ! Matrix size
      integer :: info   ! Lapack error code
      !------------------------------------------------------------
#ifdef HAS_LAPACK

      ! Store A in Ainv to prevent it from being overwritten by LAPACK
      ALU = A
      n = size(A,1)
      
      ! DGETRF computes an LU factorization of a general M-by-N matrix A
      ! using partial pivoting with row interchanges.
      call DGETRF(n, n, ALU, n, LUpivots, info)
      
      if (info /= 0) then
         print*, 'Matrix is numerically singular. ERROR:', info
         stop
      end if
#else
      STOP ':: Matrix LU factorization routine needs LAPACK.'
#endif
   end subroutine
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SolveLU(ALU,LUpivots,x,b)
      implicit none
      !------------------------------------------------------------
      real(kind=RP), dimension(:,:)        , intent(in)  :: ALU              !<  Factorized LU matrix
      integer      , dimension(size(ALU,1)), intent(in)  :: LUpivots   !<  LU pivots for factorization
      real(kind=RP), dimension(size(ALU,1)), intent(in)  :: b
      real(kind=RP), dimension(size(ALU,1)), intent(out) :: x
      !------------------------------------------------------------
      integer :: n      ! Matrix size
      integer :: info   ! Lapack error code
      !------------------------------------------------------------
#ifdef HAS_LAPACK
      
      n = size(ALU,1)
      x = b
      
      call dgetrs ( 'N' , n, 1   , ALU, n  , LUpivots, x ,n  , info )
      if (info /= 0) then
         print*,  '*** System could not be solved. ERROR:', info
         stop
      end if
#else
      STOP ':: System solving routine needs LAPACK.'
#endif
   end subroutine
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   !------------------------------------------------------------------------------
   ! Returns the inverse of a matrix calculated by LU decomposition.  Depends on LAPACK.
   !   Modified from: http://fortranwiki.org/fortran/show/Matrix+inversion
   function inverse(A) result(Ainv)
   !------------------------------------------------------------------------------
      IMPLICIT NONE
      !------------------------------------------------------------
      real(KIND=RP), dimension(:,:), intent(in) :: A
      real(KIND=RP), dimension(size(A,1),size(A,2)) :: Ainv
      
      real(KIND=RP), dimension(size(A,1)) :: work  ! work array for LAPACK
      integer, dimension(size(A,1)) :: ipiv   ! pivot indices
      integer :: n, info
      !------------------------------------------------------------
      
#ifdef HAS_LAPACK
      ! External procedures defined in LAPACK
      external DGETRF
      external DGETRI
      
      ! Store A in Ainv to prevent it from being overwritten by LAPACK
      Ainv = A
      n = size(A,1)
      
      ! DGETRF computes an LU factorization of a general M-by-N matrix A
      ! using partial pivoting with row interchanges.
      call DGETRF(n, n, Ainv, n, ipiv, info)
      
      if (info /= 0) then
         stop 'Matrix is numerically singular!'
      end if
      
      ! DGETRI computes the inverse of a matrix using the LU factorization
      ! computed by DGETRF.
      call DGETRI(n, Ainv, n, ipiv, work, n, info)
      
      if (info /= 0) then
         stop 'Lapack matrix inversion failed!'
      end if
#else
      STOP ':: Matrix inversion routine needs LAPACK.'
#endif
   !------------------------------------------------------------------------------
   end function
   !------------------------------------------------------------------------------
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   !---------------------------------------------------------------
   SUBROUTINE Mat2File(Mat,FileName,Formato)
   !       Writes a dense matrix to a File
   !---------------------------------------------------------------
      IMPLICIT NONE
      !------------------------------------------
      REAL(KIND=RP)               :: Mat(:,:)   !< Matrix
      CHARACTER(len=*)            :: FileName
      CHARACTER(len=*),OPTIONAL   :: Formato
      !------------------------------------------
      INTEGER                     :: i, j    ! Counters
      CHARACTER(len=50)           :: Forma   ! Final format
      INTEGER                     :: fd      ! File unit for writing
      !------------------------------------------
      
      IF (PRESENT(Formato)) THEN
         Forma = Formato
      ELSE
         Forma = '(300000F14.8)' !Case default
      END IF
      
      OPEN(newunit=fd,file=FileName)
      DO i=1, SIZE(Mat,1)
         WRITE(fd,Forma) (Mat(i,j),j=1,SIZE(Mat,2))
      END DO
      ClOSE(fd)
      
   !---------------------------------------------------------------
   END SUBROUTINE Mat2File
   !---------------------------------------------------------------
END MODULE DenseMatUtilities
