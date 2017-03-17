module MatrixOperations
   use SMConstants   

   interface InnerProduct
      module procedure InnerProduct2D
   end interface

   contains
      function vectorOuterProduct(a,b) result(val)
         implicit none  
         real(kind=RP)              :: a(:)
         real(kind=RP)              :: b(:)
         real(kind=RP), allocatable :: val(:,:)
         integer                    :: N , M

         N = size(a,1)
         M = size(b,1)

         allocate(val(1:N , 1:M) )

         val = spread(a(1:N),dim=2,ncopies=M)*spread(b(1:M),dim=1,ncopies=N)

      end function vectorOuterProduct

      pure subroutine MatrixTimesVector( A , X , Y , trA , reset )
!
!     -------------------------------------------------------------
!        Computes the product
!           Y = A * X (+ Y)? , to compute X*A just set trA to .true.
!     -------------------------------------------------------------
!
         implicit none
         real(kind=RP), intent(in)           :: A(:,:)
         real(kind=RP), intent(in)           :: X(:)
         logical      , intent(in), optional :: trA
         logical      , intent(in), optional :: reset
         real(kind=RP), intent(inout)        :: Y(:)
         logical                             :: tA , rst
         integer                 :: N , M , K 
!
!        ----------------------------------------
!           Variables for lapack dgemv
!        ----------------------------------------
!
#ifdef _USE_LAPACK  
         integer                 :: LDA , INCX , INCY
         real(kind=RP)           :: beta
#endif 
!

         if (present(trA)) then
            tA = trA
         else
            tA = .false.
         end if

         if ( present(reset) ) then
            rst = reset
         else
            rst = .true.
         end if 

!           Set dimensions
            M = size(A,1)
            N = size(A,2)

            if (tA) then
               K = size(A,2)
            else
               K = size(A,1)
            end if

            if ( size(Y) .ne. K ) then
!               print*, "Matrices sizes are not consistent"
!               stop "Stopped."
            end if


#ifdef _USE_LAPACK

            LDA = M
            INCX = 1
            INCY = 1

            if ( rst ) then
               beta = 0.0_RP
            else
               beta = 1.0_RP
            end if

            if (.not. tA) then
               call dgemv( "N" , M , N , 1.0_RP , A , LDA , X , INCX , beta , Y , INCY )
            elseif ( tA ) then
               call dgemv( "T" , M , N , 1.0_RP , A , LDA , X , INCX , beta , Y , INCY )
            end if
#else

            if ((.not. tA) .and. rst) then
               Y = matmul(A,X)
            elseif ((.not. tA) .and. (.not. rst) ) then
               Y = matmul(A,X) + Y
            elseif ((tA) .and. rst) then
               Y = matmul(X,A)
            elseif (tA .and. (.not. rst)) then
               Y = matmul(X,A) + Y
            end if
#endif

      end subroutine MatrixTimesVector

      pure function MatrixTimesVector_F( A , X , Nout , trA ) result( Y )
!
!     -------------------------------------------------------------
!        Computes the product
!           Y = A * X, to compute X*A just set trA to .true.
!           Nout -> Size of the output: Nrows(A) if trA = .false.
!                                       Ncols(A) if trA = .true.
!     -------------------------------------------------------------
!
         implicit none
         real(kind=RP), intent(in)           :: A(:,:)
         real(kind=RP), intent(in)           :: X(:)
         integer, intent(in)                 :: Nout
         logical      , intent(in), optional :: trA
         logical                             :: tA
         real(kind=RP)                       :: Y(1:Nout)

         if (present(trA)) then
            tA = trA
         else
            tA = .false.
         end if

         call MatrixTimesVector( A=A , X=X , Y=Y , trA = trA , reset = .true. )

      end function MatrixTimesVector_F

      function MatrixByVectorInIndex_F( A , X , M , N , index ) result( B ) 
         implicit none
         integer, intent(in)                 :: M , N
         real(kind=RP), intent(in)           :: A(M,N)
         real(kind=RP), intent(in)           :: X(:)
         integer                             :: index
         real(kind=RP)                       :: B(M,N)
         integer                             :: i , j 

         if ( index .eq. 1 ) then
           do j = 1 , N
               do i = 1 , M
                  B(i,j) = A(i,j) * X(i) 
               end do
           end do
         elseif ( index .eq. 2 ) then
           do j = 1 , N
               do i = 1 , M
                  B(i,j) = A(i,j) * X(j) 
               end do
           end do
         end if
         

      end function MatrixByVectorInIndex_F

      pure subroutine BilinearForm( A , X , Y , B , trA )
!     -----------------------------
!        Computes the product
!           B = X^T A Y
!     -----------------------------
         implicit none
         real(kind=RP), intent(in)           :: A(:,:)
         real(kind=RP), intent(in)           :: X(:)
         real(kind=RP), intent(in)           :: Y(:)
         real(kind=RP), intent(out)          :: B
         logical,       intent(in), optional :: trA
         logical                             :: tA

         if (present(trA)) then
            tA = trA
         else
            tA = .false.
         end if 
          
         if (.not. tA) then
            B = dot_product(X , MatrixTimesVector_F(A,Y,size(A,1)) ) 
         else
            B = dot_product(Y , MatrixTimesVector_F(A,X,size(A,2)) )
         end if
      end subroutine BilinearForm

      pure function BilinearForm_F( A , X , Y , trA ) result( B )
!     -----------------------------
!        Computes the product
!           B = X^T A Y
!     -----------------------------
         implicit none
         real(kind=RP), intent(in)           :: A(:,:)
         real(kind=RP), intent(in)           :: X(:)
         real(kind=RP), intent(in)           :: Y(:)
         real(kind=RP)                       :: B
         logical,       intent(in), optional :: trA
         logical                             :: tA

         if (present(trA) ) then
            call BilinearForm(A , X , Y , B , trA )
         else
            call BilinearForm(A , X , Y , B )
         end if

      end function BilinearForm_F
       
      pure subroutine Mat_x_Mat( A , B , C , trA , trB , reset )
!     -------------------------------------------------------
!        Computes the product
!           C = op(A) * op(B), in which op(A) = A or A^T
!     -------------------------------------------------------
         implicit none
         real(kind=RP), intent(in)           :: A(:,:)
         real(kind=RP), intent(in)           :: B(:,:)
         real(kind=RP), intent(out)          :: C(:,:)
         logical      , intent(in), optional :: trA
         logical      , intent(in), optional :: trB
         logical      , intent(in), optional :: reset
!        -----------------------------------------------
         logical                          :: rst
         logical                          :: tA
         logical                          :: tB
!
!        ----------------------------------------
!           Variables for lapack dgemm
!        ----------------------------------------
!
#ifdef _USE_LAPACK  
         integer                 :: N , M , K 
         integer                 :: LDA , LDB , LDC
         real(kind=RP)           :: beta
#endif 
!
         if (present(reset)) then
            rst = reset
         else
            rst = .true.
         end if

         if (present(trA)) then
            tA = trA
         else
            tA = .false.
         end if

         if (present(trB)) then
            tB = trB
         else
            tB = .false.
         end if

#ifdef _USE_LAPACK  

            if (rst) then
               beta = 0.0_RP        ! Do not add its current value
            else
               beta = 1.0_RP        ! Add its current value
            end if
!           Set dimensions
            if (tA) then
               M = size(A,2)
            else 
               M = size(A,1)
            end if

            if (tB) then
               N = size(B,1)
            else
               N = size(B,2)
            end if

            if (tA) then
               K = size(A,1)
            else
               K = size(A,2)
            end if

            if (tA) then
               LDA = K
            else
               LDA = M
            end if

            if (tB) then
               LDB = N
            else
               LDB = K
            end if

            LDC = M

            if ( (.not. tA) .and. (.not. tB) ) then
               call dgemm( "N" , "N" , M , N , K , 1.0_RP , A , LDA , B , LDB , beta , C , LDC )
            elseif ( (.not. tA) .and. ( tB ) ) then
               call dgemm( "N" , "T" , M , N , K , 1.0_RP , A , LDA , B , LDB , beta , C , LDC )
            elseif ( ( tA ) .and. (.not. tB) ) then
               call dgemm( "T" , "N" , M , N , K , 1.0_RP , A , LDA , B , LDB , beta , C , LDC )
            elseif ( ( tA ) .and. ( tB ) ) then
               call dgemm( "T" , "T" , M , N , K , 1.0_RP , A , LDA , B , LDB , beta , C , LDC )
            end if
#else 
            if ( (.not. tA) .and. (.not. tB) ) then
               if (rst) then
                  C = matmul(A,B)
               else
                  C = C + matmul(A,B)
               end if
            elseif ( (.not. tA) .and. ( tB ) ) then
               if (rst) then
                  C = matmul(A,transpose(B))
               else
                  C = C + matmul(A,transpose(B))
               end if
            elseif ( ( tA ) .and. (.not. tB) ) then
               if (rst) then
                  C = matmul(transpose(A),B)
               else
                  C = C + matmul(transpose(A),B)
               end if
            elseif ( ( tA ) .and. ( tB ) ) then
               if (rst) then
                  C = matmul(transpose(A),transpose(B))
               else
                  C = C + matmul(transpose(A),transpose(B))
               end if   
            end if
#endif 

      end subroutine Mat_x_Mat

      pure function Mat_x_Mat_F( A , B , rowC , colC ,  trA , trB) result ( C )
!     -----------------------------
!        Computes the product
!           C = op(A) * op(B)
!     -----------------------------
         implicit none
         logical      , intent (in), optional :: trA
         logical      , intent (in), optional :: trB
         integer,       intent (in)           :: rowC
         integer,       intent (in)           :: colC
         real(kind=RP), intent (in)           :: A(:,:)
         real(kind=RP), intent (in)           :: B(:,:)
         real(kind=RP)                        :: C(rowC,colC)
!        ---------------------------------------------------------------
         logical                 :: tA , tB
         integer                 :: N , M 

         if (present(trA)) then
            tA = trA
         else
            tA = .false.
         end if

         if (present(trB)) then
            tB = trB
         else
            tB = .false.
         end if


!        Set dimensions
         if (tA) then
            M = size(A,2)
         else 
            M = size(A,1)
         end if

         if (tB) then
            N = size(B,1)
         else
            N = size(B,2)
         end if

         call Mat_x_Mat( A=A , B=B , C=C , trA=trA , trB=trB , reset=.true. )

      end function Mat_x_Mat_F

      subroutine TripleMatrixProduct( A , B , C , val , trA , trB , trC )
!
!        ***********************************
!           Computes the product 
!              val = A B C
!        ***********************************
!
         implicit none
         real(kind=RP), intent(in)    :: A(:,:)
         real(kind=RP), intent(in)    :: B(:,:)
         real(kind=RP), intent(in)    :: C(:,:)
         real(kind=RP), intent(inout) :: val(:,:)
         logical,       optional      :: trA , trB , trC
         logical                      :: tA , tB , tC
         integer                      :: rowAB , colAB , colABC

         if ( present ( trA ) ) then
            tA = trA
   
         else
            tA = .false.

         end if

         if ( present ( trB ) ) then
            tB = trB
   
         else
            tB = .false.

         end if

         if ( present ( trC ) ) then
            tC = trC
   
         else
            tC = .false.

         end if

         if ( .not. tA ) then
            rowAB = size(A,1)
         else
            rowAB = size(A,2)
         end if

         if ( .not. tB ) then
            colAB = size(B,2)
         else
            colAB = size(B,1)
         end if

         if ( .not. tC ) then
            colABC = size(C,2)
         else
            colABC = size(C,1)
         end if

         val = Mat_x_Mat_F( Mat_X_Mat_F( A , B , rowAB , colAB , trA = tA , trB = tB ) , C , rowAB , colABC , trB = tC )

      end subroutine TripleMatrixProduct

      function MatrixMultiplyInIndex_F( A , B , d1C , d2C , d3C , d4C , index) result( C )
         use, intrinsic    :: iso_c_binding
!
!     ----------------------------------------------------
!        Computes the product
!           C = A(:,...i,...,:) * B(i,:)
!        for a chosen index position within 1 and 3
!     ----------------------------------------------------
!
         implicit none
         integer, intent(in)                :: d1C , d2C , d3C , d4C
         real(kind=RP), target, intent(in)  :: A(:,:,:,:)
         real(kind=RP), target, intent(in)  :: B(:,:)
         integer                            :: index
         real(kind=RP), target              :: C(d1C,d2C,d3C,d4C)
         real(kind=RP), pointer             :: PC(:,:)
         real(kind=RP), pointer             :: P1C(:)
         real(kind=RP), pointer             :: PA(:,:)
         real(kind=RP), pointer             :: P1A(:)
         integer                            :: I1 , I2 , I3 , I4
         integer                            :: i , j , k
         
         if (size(A,index) .ne. size(B,1)) then
            print*, 'Inconsistent dimensions in "MatrixMultiplyInIndex_F"'
            stop "Stopped."
         end if

         I1 = size(A,1)
         I2 = size(A,2)
         I3 = size(A,3)
         I4 = size(A,4)
         
         if (index .eq. 1) then
            I1 = size(B,2)
         elseif ( index .eq. 2) then
            I2 = size(B,2)
         elseif ( index .eq. 3) then
            I3 = size(B,2)
         end if

         if (index .eq. 1) then
            do i = 1 , I2  ; do j = 1 , I3 ; do k = 1 , I4
               C(1:I1,i,j,k) = MatrixTimesVector_F ( A=B , X=A(:,i,j,k) , Nout=size(B,2), trA=.true. ) 
            end do         ; end do        ; end do
         
         elseif (index .eq. 2) then
            do i = 1 , I3  ; do j = 1 , I4
               C(1:I1,1:I2,i,j) = Mat_X_Mat_F( A=A(:,:,i,j) , B=B , rowC = I1 , colC = I2 ) 
            end do         ; end do

         elseif (index .eq. 3) then
            do i = 1 , I4
               call c_f_pointer ( c_loc( A(:,:,:,i) ) , P1A , [size(A(:,:,:,i))] )
               call c_f_pointer ( c_loc( C(:,:,:,i) ) , P1C , [size(C(:,:,:,i))] )
               PA(1:size(A,1)*size(A,2),1:size(A,3)) => P1A(1:)
               PC(1:I1*I2,1:I3)  => P1C(1:)
               PC = Mat_X_Mat_F( A=PA , B=B , rowC = I1*I2 , colC = I3 ) 
            end do
         end if

      end function MatrixMultiplyInIndex_F

      subroutine innerProduct2D( A , M , rowA , colA , val )
!     
!     ******************************************************
!        Computes the product
!           val = tr(A) M A
!           
!           rowA and colA are the dimensions of A
!           M(rowA , rowA)
!           val(colA , colA)
!     ******************************************************
!
         implicit none
         integer                       :: rowA , colA
         real(kind=RP), intent(in)     :: A(rowA,colA)
         real(kind=RP), intent(in)     :: M(rowA,rowA)
         real(kind=RP), intent(out)     :: val(colA,colA)

      
         val = Mat_x_Mat_F ( trA = .true. , trB = .false. , A=A , B=M , rowC = colA , colC = rowA )
         val = Mat_x_Mat_F ( val , A , rowC = colA , colC = colA )

      end subroutine innerProduct2D

#if _USE_LAPACK
      function inv(A) result(Ainv)
!
!--------------------------------------------------------------
! Returns the inverse of a matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.
!     Author: Gonzalo Rubio (g.rubio@upm.es)
!--------------------------------------------------------------

        real(KIND = RP), dimension(:,:), intent(in) :: A
        real(KIND = RP), dimension(size(A,1),size(A,2)) :: Ainv
      
        real(KIND = RP), dimension(size(A,1)) :: work  ! work array for LAPACK
        integer, dimension(size(A,1)) :: ipiv   ! pivot indices
        integer :: n, info
      
        ! External procedures defined in LAPACK
        external DGETRF
        external DGETRI
      
        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv = A
        n = size(A,1)
      
        ! DGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        !PRINT*, "Call DGETRF"
        call DGETRF(n, n, Ainv, n, ipiv, info)
        !PRINT*, "info", info
        if (info /= 0) then
           stop 'Matrix is numerically singular!'
        end if
      
        ! DGETRI computes the inverse of a matrix using the LU factorization
        ! computed by DGETRF.
        call DGETRI(n, Ainv, n, ipiv, work, n, info)
      
        if (info /= 0) then
           stop 'Matrix inversion failed!'
        end if
end function inv
#endif
!
!////////////////////////////////////////////////////////////////////////
!


end module MatrixOperations
