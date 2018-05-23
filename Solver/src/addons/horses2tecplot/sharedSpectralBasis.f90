!
!//////////////////////////////////////////////////////
!
!   @File:    sharedSpectralBasis.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Sun Oct 15 13:07:03 2017
!   @Last revision date: Wed May 23 12:57:22 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 7fde177b098184b58177a3a163cefdfebe7af55f
!
!//////////////////////////////////////////////////////
!
module SharedSpectralBasis
   use SMConstants
   use NodalStorageClass
   use InterpolationMatrices
   use PolynomialInterpAndDerivsModule

   private
   public   spA
   public   ConstructSpectralBasis, addNewSpectralBasis, addNewInterpolationMatrix


   integer,   parameter    :: NMAX = 40
   type(NodalStorage_t),            allocatable  :: spA(:)


   contains

      subroutine ConstructSpectralBasis()
         implicit none
!
!        Allocate nodal storage
!        ----------------------
         allocate( spA(0:NMAX) )

      end subroutine ConstructSpectralBasis
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
!     Handling spectral approximations and transformation matrices
!     ------------------------------------------------------------
!
!/////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine addNewSpectralBasis( spA, N, nodeType)
         use NodalStorageClass
         implicit none
         class(NodalStorage_t)     :: spA(0:)
         integer                 :: N(3)
         integer, intent(in)     :: nodeType
         
         call spA(N(1)) % Construct( nodeType, N(1) )
         call spA(N(2)) % Construct( nodeType, N(2) )
         call spA(N(3)) % Construct( nodeType, N(3) )
         

      end subroutine addNewSpectralBasis

      subroutine addNewInterpolationMatrix( Tset_, Nold, spAold, Nnew, xiNew)
         use NodalStorageClass
         implicit none
         class(InterpolationMatrix_t)     :: Tset_(0:,0:)
         integer, intent(in)              :: Nold, Nnew
         class(NodalStorage_t), intent(in)  :: spAold
         real(kind=RP),       intent(in)  :: xiNew(0:Nnew)

         if ( .not. Tset_(Nnew,Nold) % Constructed ) then
            allocate ( Tset_(Nnew,Nold) % T(0:Nnew,0:Nold) )
            call PolynomialInterpolationMatrix( Nold, Nnew, spAold % x, spAold % wb, xiNew, &
                                                Tset_(Nnew,Nold) % T)
            Tset_(Nnew,Nold) % Constructed = .true.
         end if

      end subroutine addNewInterpolationMatrix

end module SharedSpectralBasis
