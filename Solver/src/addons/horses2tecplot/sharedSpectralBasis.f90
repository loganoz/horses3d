!
!//////////////////////////////////////////////////////
!
!   @File:    sharedSpectralBasis.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Sun Oct 15 13:07:03 2017
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
module SharedSpectralBasis
   use SMConstants
   use NodalStorageClass

   private
   public   spA, Tset
   public   ConstructSpectralBasis, addNewSpectralBasis, addNewInterpolationMatrix


   integer,   parameter    :: NMAX = 40


   type  InterpolationMatrices_t
      logical                          :: Constructed = .false.
      real(kind=RP),    allocatable    :: T(:,:)
   end type InterpolationMatrices_t

   type(NodalStorage),            allocatable  :: spA(:,:,:)
   type(InterpolationMatrices_t), allocatable  :: Tset(:,:)


   contains

      subroutine ConstructSpectralBasis()
         implicit none
!
!        Allocate nodal storage and interpolation matrices
!        -------------------------------------------------
         allocate( spA(0:NMAX, 0:NMAX, 0:NMAX) )
         allocate( Tset(0:NMAX, 0:NMAX) )

      end subroutine ConstructSpectralBasis
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
!     Handling spectral approximations and transformation matrices
!     ------------------------------------------------------------
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine addNewSpectralBasis( spA, N)
         use NodalStorageClass
         implicit none
         class(NodalStorage)     :: spA(0:,0:,0:)
         integer                 :: N(3)

         if ( .not. spA(N(1),N(2), N(3)) % Constructed ) then
            call spA(N(1), N(2), N(3) ) % Construct( N(1), N(2), N(3) )
         end if

      end subroutine addNewSpectralBasis

      subroutine addNewInterpolationMatrix( Tset, Nold, spAold, Nnew, xiNew)
         use NodalStorageClass
         implicit none
         class(InterpolationMatrices_t)     :: Tset(0:,0:)
         integer, intent(in)              :: Nold, Nnew
         class(NodalStorage), intent(in)  :: spAold
         real(kind=RP),       intent(in)  :: xiNew(0:Nnew)

         if ( .not. Tset(Nnew,Nold) % Constructed ) then
            allocate ( Tset(Nnew,Nold) % T(0:Nnew,0:Nold) )
            call PolynomialInterpolationMatrix( Nold, Nnew, spAold % xi, spAold % wbx, xiNew, &
                                                Tset(Nnew,Nold) % T)
            Tset(Nnew,Nold) % Constructed = .true.
         end if

      end subroutine addNewInterpolationMatrix

end module SharedSpectralBasis
