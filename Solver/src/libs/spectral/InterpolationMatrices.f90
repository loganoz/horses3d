!
!////////////////////////////////////////////////////////////////////////
!
!      InterpolationAndDerivatives.f90
!      Created: 2009-12-15 15:36:24 -0500 
!      By: David Kopriva  
!
!
!////////////////////////////////////////////////////////////////////////
!
module InterpolationMatrices
   use PolynomialInterpAndDerivsModule
   use NodalStorageClass
   use SMConstants
   implicit none
   
   private
   public   Tset, InterpolationMatrix_t, ConstructInterpolationMatrices, Interp3DArrays, Interp3DArraysOneDir

!
!  ---------------------------------
!  Interpolation matrix derived type
!  ---------------------------------
!
   type InterpolationMatrix_t
      logical                :: Constructed = .false.
      real(kind=RP), pointer :: T(:,:)
   end type InterpolationMatrix_t
!
!  **************************************************
!  The set of interpolation matrices are stored here.
!  >> The following criteria is adopted:
!        * Tset(N,M) with N=M is simply a pointer to the identity matrix.
!        * Tset(N,M) with N<M is a forward matrix.
!        * Tset(N,M) with N>M is a backwards matrix.
!  **************************************************
   integer                    , parameter :: NMAX = 40
   type(InterpolationMatrix_t), target    :: Tset(0:NMAX,0:NMAX)
   type(InterpolationMatrix_t), target    :: Identity
   
!========
 contains
!========

subroutine ConstructInterpolationMatrices(Norigin, Ndest)
!
!     ****************************************************************
!        This subroutine computes the interpolation matrices
!        from Norigin to Ndest (forward) and backwards
!     TODO: If other spA's are needed, add spA as an optional argument 
!     ****************************************************************
!
      implicit none
      integer, intent(in)  :: Norigin
      integer, intent(in)  :: Ndest
!
!     ---------------
!     Local variables
!     ---------------
!
      integer  :: i, j
      
!
!     Evaluate if it's necessary to construct the matrix
!     --------------------------------------------------
      
      if ( Tset(Norigin, Ndest) % Constructed ) return
      
      if ( Norigin .eq. Ndest ) then
         if (.not. Identity % Constructed) call ConstructIdentity
         Tset(Ndest,Norigin) % T(0:,0:) => Identity % T(0:Ndest,0:Norigin)
         Tset(Ndest,Norigin) % Constructed = .TRUE.
         return
      end if
!
!     Matrix construction
!     -------------------      
      
      associate ( spAo => NodalStorage(Norigin), &
                  spAd => NodalStorage(Ndest)      )
      
!
!        Allocate memory
!        ---------------
      allocate( Tset(Norigin, Ndest) % T(0:Ndest, 0:Norigin) )
      allocate( Tset(Ndest, Norigin) % T(0:Norigin, 0:Ndest) )
!
!        Construct the forward interpolation matrix
!        ------------------------------------------
      if (Norigin < Ndest) then
         call PolynomialInterpolationMatrix(Norigin, Ndest, spAo % x, spAo % wb, spAd % x, &
                                            Tset(Norigin, Ndest) % T)
      else
         call PolynomialInterpolationMatrix(Ndest, Norigin, spAd % x, spAd % wb, spAo % x, &
                                            Tset(Ndest, Norigin) % T)
      end if
!
!        Construct the backwards interpolation matrix
!        --------------------------------------------
      if (Norigin < Ndest) then
         do j = 0, Ndest   ; do i = 0, Norigin
            Tset(Ndest,Norigin) % T(i,j) = Tset(Norigin,Ndest) % T(j,i) * spAd % w(j) / spAo % w(i)
         end do            ; end do
      else
         do j = 0, Norigin ; do i = 0, Ndest
            Tset(Norigin,Ndest) % T(i,j) = Tset(Ndest,Norigin) % T(j,i) * spAo % w(j) / spAd % w(i)
         end do            ; end do
      end if
!
!        Set constructed flag
!        --------------------          
      Tset(Norigin,Ndest) % Constructed = .true.
      Tset(Ndest,Norigin) % Constructed = .true.

      end associate
      
   end subroutine ConstructInterpolationMatrices
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------------------
!  Subroutine for interpolating an array of 3D arrays
!  --------------------------------------------------
   subroutine Interp3DArrays(Nvars, Nin, inArray, Nout, outArray, interpXi, interpEta, interpZeta )
      implicit none
      !-------------------------------------------------------
      integer                                                        , intent(in)  :: Nvars
      integer      , dimension(3)                                    , intent(in)  :: Nin
      integer      , dimension(3)                                    , intent(in)  :: Nout
      real(kind=RP), dimension(Nvars,0:Nin (1), 0:Nin (2), 0:Nin (3)), intent(in)  :: inArray
      real(kind=RP), dimension(Nvars,0:Nout(1), 0:Nout(2), 0:Nout(3)), intent(out) :: outArray
      real(kind=RP), dimension(0:Nout(1), 0:Nin(1)), optional, target, intent(in)  :: interpXi
      real(kind=RP), dimension(0:Nout(2), 0:Nin(2)), optional, target, intent(in)  :: interpEta
      real(kind=RP), dimension(0:Nout(3), 0:Nin(3)), optional, target, intent(in)  :: interpZeta
      !-------------------------------------------------------
      real(kind=RP), dimension(:,:), pointer :: Txi
      real(kind=RP), dimension(:,:), pointer :: Teta
      real(kind=RP), dimension(:,:), pointer :: Tzeta
      integer :: i,j,k,l,m,n
      !-------------------------------------------------------
      
      if (present(interpXi) .and. present(interpEta) .and. present(interpZeta)) then
         Txi   => interpXi
         Teta  => interpEta
         Tzeta => interpZeta
      else
         Txi   => Tset(Nin(1),Nout(1)) % T
         Teta  => Tset(Nin(2),Nout(2)) % T
         Tzeta => Tset(Nin(3),Nout(3)) % T
      end if
      
      outArray = 0.0_RP
      
      do n = 0, Nin(3)  ; do k = 0, Nout(3)
         do m = 0, Nin(2)  ; do j = 0, Nout(2)   
            do l = 0, Nin(1)  ; do i = 0, Nout(1)
               outArray(:,i,j,k) = outArray(:,i,j,k) +   Txi   (i,l) &
                                                       * Teta  (j,m) &
                                                       * Tzeta (k,n) &
                                                       * inArray(:,l,m,n)
            end do             ; end do
         end do             ; end do
      end do             ; end do
      
   end subroutine Interp3DArrays
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------------------
!  Subroutine for interpolating an array of 3D arrays
!  --------------------------------------------------
   subroutine Interp3DArraysOneDir(Nvars, Nin, inArray, Nout, outArray, Dir, interpMat )
      implicit none
      !-------------------------------------------------------
      integer                                                            , intent(in)  :: Nvars
      integer      , dimension(3)                                        , intent(in)  :: Nin
      integer      , dimension(3)                                        , intent(in)  :: Nout
      integer                                                            , intent(in)  :: Dir
      real(kind=RP), dimension(Nvars,0:Nin (1), 0:Nin (2), 0:Nin (3))    , intent(in)  :: inArray
      real(kind=RP), dimension(Nvars,0:Nout(1), 0:Nout(2), 0:Nout(3))    , intent(out) :: outArray
      real(kind=RP), dimension(0:Nout(Dir), 0:Nin(Dir)), optional, target, intent(in)  :: interpMat
      !-------------------------------------------------------
      real(kind=RP), dimension(:,:), pointer :: T     ! Interpolation matrix
      integer :: i,j,k,l
      !-------------------------------------------------------
      
      if ( present(interpMat) ) then
         T   => interpMat
      else
         T   => Tset(Nin(Dir),Nout(Dir)) % T
      end if
      
      outArray = 0.0_RP
      
      select case(Dir)
         case (1)
            
            do k = 0, Nout(3)
               do j = 0, Nout(2)   
                  do l = 0, Nin(1)  ; do i = 0, Nout(1)
                     outArray(:,i,j,k) = outArray(:,i,j,k) + T(i,l) * inArray(:,l,j,k)
                  end do             ; end do
               end do
            end do
            
         case (2)
            
            do k = 0, Nout(3)
               do l = 0, Nin(2)  ; do j = 0, Nout(2)   
                  do i = 0, Nout(1)
                     outArray(:,i,j,k) = outArray(:,i,j,k) + T(j,l) * inArray(:,i,l,k)
                  end do
               end do             ; end do
            end do
            
         case (3)
            
            do l = 0, Nin(3)  ; do k = 0, Nout(3)
               do j = 0, Nout(2)   
                  do i = 0, Nout(1)
                     outArray(:,i,j,k) = outArray(:,i,j,k) + T (k,l) * inArray(:,i,j,l)
                  end do
               end do
            end do             ; end do
            
      end select
      
   end subroutine Interp3DArraysOneDir
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine ConstructIdentity
      implicit none
      !----------------------------
      integer :: i
      !----------------------------
      
      allocate (Identity % T(0:NMAX,0:NMAX))
      
      Identity % T = 0._RP
      do i=0, NMAX
         Identity % T (i,i) = 1._RP
      end do
      
      Identity % Constructed = .TRUE.
   end subroutine ConstructIdentity
   
end module InterpolationMatrices
