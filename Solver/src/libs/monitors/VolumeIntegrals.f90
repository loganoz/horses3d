module VolumeIntegrals
   use SMConstants
   use Physics
   use NodalStorageClass
   use HexMeshClass
#include "Includes.h"
   
   private
   public   VOLUME, KINETIC_ENERGY, KINETIC_ENERGY_RATE, ENSTROPHY
   public   ScalarVolumeIntegral, VectorVolumeIntegral

   integer, parameter      :: VOLUME = 1
   integer, parameter      :: KINETIC_ENERGY = 2
   integer, parameter      :: KINETIC_ENERGY_RATE = 3
   integer, parameter      :: ENSTROPHY = 4
!
!  ========
   contains
!  ========
!
!////////////////////////////////////////////////////////////////////////////////////////
!
!           SCALAR INTEGRALS PROCEDURES
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      function ScalarVolumeIntegral(mesh, spA, integralType) result(val)
!
!        -----------------------------------------------------------
!           This function computes scalar integrals, that is, those
!           in the form:
!                 val = \int \vec{v}·\vec{n}dS
!           Implemented integrals are:
!              * Volume: computes the zone surface.
!              * Mass flow: computes the mass flow across the zone.
!              * Flow: computes the volumetric flow across the zone.
!        -----------------------------------------------------------
!
         implicit none
         class(HexMesh),      intent(in)  :: mesh
         class(NodalStorage), intent(in)  :: spA(0:,0:,0:)
         integer,             intent(in)  :: integralType
         real(kind=RP)                    :: val
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: eID
!
!        Initialization
!        --------------            
         val = 0.0_RP
!
!        Loop the mesh
!        -------------
!$omp parallel do reduction(+:val) private(eID) schedule(guided)
         do eID = 1, mesh % no_of_elements
!
!           Compute the integral
!           --------------------
            val = val + ScalarVolumeIntegral_Local(mesh % elements(eID), &
                                                                    spA, &
                                                           integralType    )

         end do
!$omp end parallel do

      end function ScalarVolumeIntegral

      function ScalarVolumeIntegral_Local(e, spA, integralType) result(val)
         implicit none
         class(Element),      target, intent(in)     :: e
         class(NodalStorage), target, intent(in)     :: spA(0:,0:,0:)
         integer,                     intent(in)     :: integralType
         real(kind=RP)                               :: val
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: Nel(3)    ! Element polynomial order
         integer     :: i, j, k
         real(kind=RP)           :: KinEn(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)           :: uvw(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP), pointer  :: Qb(:)

         Nel = e % Nxyz

         associate ( wx => spA(Nel(1),Nel(2),Nel(3)) % wx, &
                     wy => spA(Nel(1),Nel(2),Nel(3)) % wy, &
                     wz => spA(Nel(1),Nel(2),Nel(3)) % wz    )
!
!        Initialization
!        --------------
         val = 0.0_RP
!
!        Perform the numerical integration
!        ---------------------------------
         select case ( integralType )
   
         case ( VOLUME )
!
!           **********************************
!           Computes the volume integral
!              val = \int dV
!           **********************************
!
            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               val = val + wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k)
            end do            ; end do           ; end do

         case ( KINETIC_ENERGY )
!
!           ***********************************
!              Computes the kinetic energy
!              integral:
!              K = \int \rho V^2 dV
!           ***********************************
!
            
            KinEn =         POW2(e % storage % Q(:,:,:,IRHOU)) 
            KinEn = KinEn + POW2( e % storage % Q(:,:,:,IRHOV) )
            KinEn = KinEn + POW2( e % storage % Q(:,:,:,IRHOW) )
            KinEn = 0.5_RP * KinEn / e % storage % Q(:,:,:,IRHO)

            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               val = val +   wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * kinEn(i,j,k)
            end do            ; end do           ; end do

         case ( KINETIC_ENERGY_RATE )
!
!           ***********************************
!              Computes the kinetic energy
!              time derivative:
!              K_t = (d/dt)\int \rho V^2 dV
!           ***********************************
!
            uvw = e % storage % Q(:,:,:,IRHOU) / e % storage % Q(:,:,:,IRHO)
            KinEn = uvw * e % storage % QDot(:,:,:,IRHOU) - 0.5_RP * POW2(uvw) * e % storage % QDot(:,:,:,IRHO)

            uvw = e % storage % Q(:,:,:,IRHOV) / e % storage % Q(:,:,:,IRHO)
            KinEn = KinEn + uvw * e % storage % QDot(:,:,:,IRHOV) - 0.5_RP * POW2(uvw) * e % storage % QDot(:,:,:,IRHO)

            uvw = e % storage % Q(:,:,:,IRHOW) / e % storage % Q(:,:,:,IRHO)
            KinEn = KinEn + uvw * e % storage % QDot(:,:,:,IRHOW) - 0.5_RP * POW2(uvw) * e % storage % QDot(:,:,:,IRHO)

            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               val = val +   wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * kinEn(i,j,k)
            end do            ; end do           ; end do


         case ( ENSTROPHY )
!
!           ***************************
!           Computes the flow enstrophy
!           ***************************
!
            KinEn =   POW2( e % storage % U_y(:,:,:,IRHOW) - e % storage % U_z(:,:,:,IRHOV) ) &
                    + POW2( e % storage % U_z(:,:,:,IRHOU) - e % storage % U_x(:,:,:,IRHOW) ) &
                    + POW2( e % storage % U_x(:,:,:,IRHOV) - e % storage % U_y(:,:,:,IRHOU) )

            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               val = val +   wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * kinEn(i,j,k)
            end do            ; end do           ; end do

         end select

         end associate

      end function ScalarVolumeIntegral_Local

end module VolumeIntegrals
