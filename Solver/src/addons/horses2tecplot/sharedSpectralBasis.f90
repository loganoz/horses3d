!
!//////////////////////////////////////////////////////
!
!   @File:    sharedSpectralBasis.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Sun Oct 15 13:07:03 2017
!   @Last revision date: Wed Oct 25 18:53:00 2017
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 5edaf46ab67ee96cdf80ff143c0ab65970c05b73
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
   type(NodalStorage),            allocatable  :: spA(:)


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
         class(NodalStorage)     :: spA(0:)
         integer                 :: N(3)
         integer, intent(in)     :: nodeType
         
         call spA(N(1)) % Construct( nodeType, N(1) )
         call spA(N(2)) % Construct( nodeType, N(2) )
         call spA(N(3)) % Construct( nodeType, N(3) )
         

      end subroutine addNewSpectralBasis

      subroutine addNewInterpolationMatrix( Tset, Nold, spAold, Nnew, xiNew)
         use NodalStorageClass
         implicit none
         class(InterpolationMatrix_t)     :: Tset(0:,0:)
         integer, intent(in)              :: Nold, Nnew
         class(NodalStorage), intent(in)  :: spAold
         real(kind=RP),       intent(in)  :: xiNew(0:Nnew)

         if ( .not. Tset(Nnew,Nold) % Constructed ) then
            allocate ( Tset(Nnew,Nold) % T(0:Nnew,0:Nold) )
            call PolynomialInterpolationMatrix( Nold, Nnew, spAold % x, spAold % wb, xiNew, &
                                                Tset(Nnew,Nold) % T)
            Tset(Nnew,Nold) % Constructed = .true.
         end if

      end subroutine addNewInterpolationMatrix

end module SharedSpectralBasis
