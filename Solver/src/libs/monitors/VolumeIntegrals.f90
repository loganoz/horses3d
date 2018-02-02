module VolumeIntegrals
   use SMConstants
   use Physics
   use VariableConversion
   use HexMeshClass
#ifdef _HAS_MPI_
   use mpi
#endif
#include "Includes.h"
   
   private
   public   VOLUME

#if defined(NAVIERSTOKES)
   public KINETIC_ENERGY, KINETIC_ENERGY_RATE, ENSTROPHY, VELOCITY
   public ENTROPY, ENTROPY_RATE

#elif defined(CAHNHILLIARD)
   public FREE_ENERGY

#endif

   public   ScalarVolumeIntegral, VectorVolumeIntegral


   integer, parameter      :: VOLUME              = 1

#if defined(NAVIERSTOKES)
   integer, parameter      :: KINETIC_ENERGY      = 2
   integer, parameter      :: KINETIC_ENERGY_RATE = 3
   integer, parameter      :: ENSTROPHY           = 4
   integer, parameter      :: VELOCITY            = 5
   integer, parameter      :: ENTROPY             = 6
   integer, parameter      :: ENTROPY_RATE        = 7

#elif defined(CAHNHILLIARD)   
   integer, parameter      :: FREE_ENERGY         = 2

#endif
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
      function ScalarVolumeIntegral(mesh, integralType) result(val)
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
         integer,             intent(in)  :: integralType
         real(kind=RP)                    :: val, localVal
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: eID, ierr
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
                                                           integralType    )

         end do
!$omp end parallel do

#ifdef _HAS_MPI_
            localVal = val
            call mpi_allreduce(localVal, val, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

      end function ScalarVolumeIntegral

      function ScalarVolumeIntegral_Local(e, integralType) result(val)
         implicit none
         class(Element),      target, intent(in)     :: e
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
         real(kind=RP)           :: p, s, dtP
         real(kind=RP), pointer  :: Qb(:)
         real(kind=RP)           :: free_en, fchem

         Nel = e % Nxyz

         associate ( wx => e % spAxi % w, &
                     wy => e % spAeta % w, &
                     wz => e % spAzeta % w    )
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

#if defined(NAVIERSTOKES)
         case ( KINETIC_ENERGY )
!
!           ***********************************
!              Computes the kinetic energy
!              integral:
!              K = \int \rho V^2 dV
!           ***********************************
!
            
            KinEn =         POW2(e % storage % Q(IRHOU,:,:,:)) 
            KinEn = KinEn + POW2( e % storage % Q(IRHOV,:,:,:) )
            KinEn = KinEn + POW2( e % storage % Q(IRHOW,:,:,:) )
            KinEn = 0.5_RP * KinEn / e % storage % Q(IRHO,:,:,:)

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
            uvw = e % storage % Q(IRHOU,:,:,:) / e % storage % Q(IRHO,:,:,:)
            KinEn = uvw * e % storage % QDot(IRHOU,:,:,:) - 0.5_RP * POW2(uvw) * e % storage % QDot(IRHO,:,:,:)

            uvw = e % storage % Q(IRHOV,:,:,:) / e % storage % Q(IRHO,:,:,:)
            KinEn = KinEn + uvw * e % storage % QDot(IRHOV,:,:,:) - 0.5_RP * POW2(uvw) * e % storage % QDot(IRHO,:,:,:)

            uvw = e % storage % Q(IRHOW,:,:,:) / e % storage % Q(IRHO,:,:,:)
            KinEn = KinEn + uvw * e % storage % QDot(IRHOW,:,:,:) - 0.5_RP * POW2(uvw) * e % storage % QDot(IRHO,:,:,:)

            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               val = val +   wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * kinEn(i,j,k)
            end do            ; end do           ; end do


         case ( ENSTROPHY )
!
!           ***************************
!           Computes the flow enstrophy
!           ***************************
!
            KinEn =   POW2( e % storage % U_y(IGW,:,:,:) - e % storage % U_z(IGV,:,:,:) ) &
                    + POW2( e % storage % U_z(IGU,:,:,:) - e % storage % U_x(IGW,:,:,:) ) &
                    + POW2( e % storage % U_x(IGV,:,:,:) - e % storage % U_y(IGU,:,:,:) )

            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               val = val +   wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * kinEn(i,j,k)
            end do            ; end do           ; end do

         case ( VELOCITY )

            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               val = val +   wx(i) * wy(j) * wz(k) * sqrt(   POW2(e % storage % Q(IRHOU,i,j,k)) &
                                                           + POW2(e % storage % Q(IRHOV,i,j,k)) & 
                                                           + POW2(e % storage % Q(IRHOW,i,j,k))  ) &
                                        / e % storage % Q(IRHO,i,j,k) * e % geom % jacobian(i,j,k)
            end do            ; end do           ; end do

         case ( ENTROPY )
!
!           ********************************************
!              Computes the specific entropy integral
!           ********************************************
!
            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               p = Pressure( e % storage % Q )
               s = ( log(p) - thermodynamics % gamma * log(e % storage % Q(IRHO,i,j,k)) )
               val = val +   wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * s
            end do            ; end do           ; end do

          case ( ENTROPY_RATE )
!
!           ************************************************************
!              Computes the specific entropy integral time derivative
!           ************************************************************
!
            uvw = e % storage % Q(IRHOU,:,:,:) / e % storage % Q(IRHO,:,:,:)
            KinEn = uvw * e % storage % QDot(IRHOU,:,:,:) - 0.5_RP * POW2(uvw) * e % storage % QDot(IRHO,:,:,:)

            uvw = e % storage % Q(IRHOV,:,:,:) / e % storage % Q(IRHO,:,:,:)
            KinEn = KinEn + uvw * e % storage % QDot(IRHOV,:,:,:) - 0.5_RP * POW2(uvw) * e % storage % QDot(IRHO,:,:,:)

            uvw = e % storage % Q(IRHOW,:,:,:) / e % storage % Q(IRHO,:,:,:)
            KinEn = KinEn + uvw * e % storage % QDot(IRHOW,:,:,:) - 0.5_RP * POW2(uvw) * e % storage % QDot(IRHO,:,:,:)

            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               p = Pressure( e % storage % Q )
               dtP = thermodynamics % gammaMinus1 * ( e % storage % QDot(IRHOE,i,j,k) - KinEn(i,j,k) )
               s   =  dtP/p - thermodynamics % gamma * e % storage % QDot(IRHO,i,j,k) / e % storage % Q(IRHO,i,j,k) 
               val = val +   wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * s
            end do            ; end do           ; end do
            
#elif defined(CAHNHILLIARD)
         case (FREE_ENERGY)

            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               call QuarticDWP(e % storage % c(i,j,k), thermodynamics % c_alpha, thermodynamics % c_beta, fchem)
               free_en = fchem + 0.5_RP * POW2(dimensionless % eps) *(   POW2(e % storage % gradC(1,i,j,k)) &
                                                                       + POW2(e % storage % gradC(2,i,j,k)) &
                                                                       + POW2(e % storage % gradC(3,i,j,k)) )
               val = val + wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * free_en
            end do            ; end do           ; end do
           
#endif
         end select

         end associate

      end function ScalarVolumeIntegral_Local

end module VolumeIntegrals
