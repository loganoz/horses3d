module VolumeIntegrals
   use SMConstants
   use PhysicsStorage
   use Physics
   use VariableConversion
   use ElementClass
   use HexMeshClass
   use FluidData
#ifdef _HAS_MPI_
   use mpi
#endif
#include "Includes.h"
   
   private
   public   VOLUME

#if defined(NAVIERSTOKES)
   public KINETIC_ENERGY, KINETIC_ENERGY_RATE, ENSTROPHY, VELOCITY
   public ENTROPY, ENTROPY_RATE, INTERNAL_ENERGY, MOMENTUM, SOURCE, PSOURCE
#endif

#if defined(INCNS)
   public MASS, ENTROPY, KINETIC_ENERGY_RATE, ENTROPY_RATE
#endif

#if defined(MULTIPHASE)
   public ENTROPY_RATE
#endif

#if defined(CAHNHILLIARD)
   public FREE_ENERGY
#endif

   public   ScalarVolumeIntegral, VectorVolumeIntegral



   enum, bind(C)
      enumerator :: VOLUME  
#if defined(NAVIERSTOKES)
      enumerator :: KINETIC_ENERGY, KINETIC_ENERGY_RATE
      enumerator :: ENSTROPHY, VELOCITY, ENTROPY, ENTROPY_RATE, INTERNAL_ENERGY, MOMENTUM, SOURCE, PSOURCE
#endif
#if defined(INCNS)
      enumerator :: MASS, ENTROPY, KINETIC_ENERGY_RATE, ENTROPY_RATE
#endif
#if defined(MULTIPHASE)
      enumerator :: ENTROPY_RATE
#endif
#if defined(CAHNHILLIARD)
      enumerator :: FREE_ENERGY
#endif
   end enum
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
!                 val = \int v dx
!           Implemented integrals are:
!              * VOLUME
!              * KINETIC_ENERGY
!              * KINETIC_ENERGY_RATE
!              * ENSTROPHY
!              * VELOCITY
!              * ENTROPY
!              * ENTROPY_RATE
!              * INTERNAL ENERGY
!              * FREE_ENERGY            
!        -----------------------------------------------------------
!

         implicit none
         class(HexMesh),      intent(in)  :: mesh
         integer,             intent(in)  :: integralType
         real(kind=RP)                    :: val
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: localVal
         integer       :: eID, ierr
         
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
         real(kind=RP)           :: U_x(NDIM,0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)           :: U_y(NDIM,0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)           :: U_z(NDIM,0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)           :: uvw(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)           :: p, s, dtP
         real(kind=RP), pointer  :: Qb(:)
         real(kind=RP)           :: free_en, fchem, entr
         real(kind=RP)           :: Strain(NDIM,NDIM)
         real(kind=RP)           :: mu

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
            call getVelocityGradients(e % Nxyz, e % storage % Q,e % storage % U_x,e % storage % U_y,e % storage % U_z, U_x, U_y, U_z)
            
            KinEn =   POW2( U_y(IZ,:,:,:) - U_z(IY,:,:,:) ) &
                    + POW2( U_z(IX,:,:,:) - U_x(IZ,:,:,:) ) &
                    + POW2( U_x(IY,:,:,:) - U_y(IX,:,:,:) )

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
         case ( INTERNAL_ENERGY )
            !
            !           ***********************************
            !              Computes the internal energy
            !              integral:
            !              \rho e = \int \rho e dV
            !           ***********************************
            !
            
            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               val = val + wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * e % storage % Q(IRHOE,i,j,k)
            end do            ; end do           ; end do
            
#endif
#if defined(INCNS)
         case (MASS)
            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               val = val +   wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * e % storage % Q(INSRHO,i,j,k)
            end do            ; end do           ; end do

         case (ENTROPY)
            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               entr =   0.5_RP*(sum(POW2(e % storage % Q(INSRHOU:INSRHOW,i,j,k))))/e % storage % Q(INSRHO,i,j,k) &
                         + 0.5_RP * POW2(e % storage % Q(INSP,i,j,k)) / thermodynamics % rho0c02
               val = val +   wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * entr
            end do            ; end do           ; end do

         case (KINETIC_ENERGY_RATE)
!
!           ***********************************
!              Computes the kinetic energy
!              time derivative:
!              K_t = (d/dt)\int \rho V^2 dV
!           ***********************************
!
            uvw = e % storage % Q(INSRHOU,:,:,:) / e % storage % Q(INSRHO,:,:,:)
            KinEn = uvw * e % storage % QDot(INSRHOU,:,:,:) - 0.5_RP * POW2(uvw) * e % storage % QDot(INSRHO,:,:,:)

            uvw = e % storage % Q(INSRHOV,:,:,:) / e % storage % Q(INSRHO,:,:,:)
            KinEn = KinEn + uvw * e % storage % QDot(INSRHOV,:,:,:) - 0.5_RP * POW2(uvw) * e % storage % QDot(INSRHO,:,:,:)

            uvw = e % storage % Q(INSRHOW,:,:,:) / e % storage % Q(INSRHO,:,:,:)
            KinEn = KinEn + uvw * e % storage % QDot(INSRHOW,:,:,:) - 0.5_RP * POW2(uvw) * e % storage % QDot(INSRHO,:,:,:)

            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               val = val +   wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * kinEn(i,j,k)
            end do            ; end do           ; end do

         case (ENTROPY_RATE)
!
!           ***********************************
!              Computes the kinetic energy
!              time derivative:
!              K_t = (d/dt)\int \rho V^2 dV
!           ***********************************
!
            uvw = e % storage % Q(INSRHOU,:,:,:) / e % storage % Q(INSRHO,:,:,:)
            KinEn = uvw * e % storage % QDot(INSRHOU,:,:,:) - 0.5_RP * POW2(uvw) * e % storage % QDot(INSRHO,:,:,:)

            uvw = e % storage % Q(INSRHOV,:,:,:) / e % storage % Q(INSRHO,:,:,:)
            KinEn = KinEn + uvw * e % storage % QDot(INSRHOV,:,:,:) - 0.5_RP * POW2(uvw) * e % storage % QDot(INSRHO,:,:,:)

            uvw = e % storage % Q(INSRHOW,:,:,:) / e % storage % Q(INSRHO,:,:,:)
            KinEn = KinEn + uvw * e % storage % QDot(INSRHOW,:,:,:) - 0.5_RP * POW2(uvw) * e % storage % QDot(INSRHO,:,:,:)

            KinEn = KinEn + (1.0_RP/thermodynamics % rho0c02)*e % storage % Q(INSP,:,:,:)*e % storage % QDot(INSP,:,:,:)

            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               val = val +   wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * kinEn(i,j,k)
            end do            ; end do           ; end do


#endif
#ifdef MULTIPHASE
         case (ENTROPY_RATE)
!
!           ***********************************
!              Computes the kinetic energy
!              time derivative:
!              K_t = (d/dt)\int \rho V^2 dV
!           ***********************************
!
            KinEn = e % storage % QDot(IMC,:,:,:) * e % storage % mu(1,:,:,:)
            KinEn = KinEn + e % storage % Q(IMSQRHOU,:,:,:)*e % storage % QDot(IMSQRHOU,:,:,:)
            KinEn = KinEn + e % storage % Q(IMSQRHOV,:,:,:)*e % storage % QDot(IMSQRHOV,:,:,:)
            KinEn = KinEn + e % storage % Q(IMSQRHOW,:,:,:)*e % storage % QDot(IMSQRHOW,:,:,:)
            KinEn = KinEn + dimensionless % Ma2*e % storage % Q(IMP,:,:,:)*e % storage % QDot(IMP,:,:,:)

            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)

!               mu = dimensionless % mu(2) + (dimensionless % mu(1) - dimensionless % mu(2))*e % storage % Q(IMC,i,j,k)
!               Strain(1,1) = e % storage % U_x(IGU,i,j,k)
!               Strain(2,2) = e % storage % U_y(IGV,i,j,k)
!               Strain(3,3) = e % storage % U_z(IGW,i,j,k)
!
!               Strain(1,2) = 0.5_RP * (e % storage % U_x(IGV,i,j,k) + e % storage % U_y(IGU,i,j,k))
!               Strain(1,3) = 0.5_RP * (e % storage % U_x(IGW,i,j,k) + e % storage % U_z(IGU,i,j,k))
!               Strain(2,3) = 0.5_RP * (e % storage % U_y(IGW,i,j,k) + e % storage % U_z(IGV,i,j,k))
!         
!               Strain(2,1) = Strain(1,2)
!               Strain(3,1) = Strain(1,3)
!               Strain(3,2) = Strain(2,3)
               Strain = 0.0_RP

               val = val +   wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * (kinEn(i,j,k) + 2.0_RP*mu*sum(Strain*Strain))
            end do            ; end do           ; end do

#endif
#if defined(CAHNHILLIARD)            
         case (FREE_ENERGY)

            do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
               call QuarticDWP(e % storage % c(1,i,j,k), fchem)
               free_en = fchem + 0.5_RP * POW2(multiphase % eps) *(   POW2(e % storage % c_x(1,i,j,k)) &
                                                                       + POW2(e % storage % c_y(1,i,j,k)) &
                                                                       + POW2(e % storage % c_z(1,i,j,k)) )
               val = val + wx(i) * wy(j) * wz(k) * e % geom % jacobian(i,j,k) * free_en
            end do            ; end do           ; end do
           
#endif
         end select

         end associate

      end function ScalarVolumeIntegral_Local
!
!////////////////////////////////////////////////////////////////////////////////////////
!
!           VECTOR INTEGRALS PROCEDURES
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      function VectorVolumeIntegral(mesh, integralType, num_of_vars) result(val)
!
!        -----------------------------------------------------------
!           This function computes vector integrals, that is, those
!           in the form:
!                 val = \int \vec{v} dx
!           Implemented integrals are:
!              * VELOCITY
!              * MOMENTUM
!        -----------------------------------------------------------
!
         implicit none
         class(HexMesh),      intent(in)  :: mesh
         integer,             intent(in)  :: integralType
         integer,             intent(in)  :: num_of_vars
         real(kind=RP)                    :: val(num_of_vars)
!
!        ---------------
!        Local variables
!        ---------------
!
         logical       :: fiveVars
         integer       :: eID, ierr
         real(kind=RP) :: localVal(num_of_vars)
         real(kind=RP) :: valAux(num_of_vars)
         real(kind=RP) :: val1, val2, val3, val4, val5
         
         if (num_of_vars == 5) then    ! Ugly hack.. But only way to make it work with ifort....
            fiveVars = .TRUE.
         else
            fiveVars = .FALSE.
         end if
         
!
!        Initialization
!        --------------            
         val1 = 0.0_RP
         val2 = 0.0_RP
         val3 = 0.0_RP
         val4 = 0.0_RP
         val5 = 0.0_RP
!
!        Loop the mesh
!        -------------
!$omp parallel do reduction(+:val1,val2,val3,val4,val5) private(valAux) schedule(guided)
         do eID = 1, mesh % no_of_elements
!
!           Compute the integral
!           --------------------
            valAux = VectorVolumeIntegral_Local(mesh % elements(eID), integralType  , num_of_vars  )
            
            val1 = val1 + valAux(1)
            val2 = val2 + valAux(2)
            val3 = val3 + valAux(3)
            
            if (fiveVars) then
               val4 = val4 + valAux(4)
               val5 = val5 + valAux(5)
            end if
         end do
!$omp end parallel do

         val(1:3) = [val1, val2, val3]
         
         if (fiveVars) val(4:5) = [val4, val5] 
         
#ifdef _HAS_MPI_
         localVal = val
         call mpi_allreduce(localVal, val, num_of_vars, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

      end function VectorVolumeIntegral
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      function VectorVolumeIntegral_Local(e, integralType, num_of_vars) result(val)
         implicit none
         !-arguments---------------------------------------------------
         class(Element),      target, intent(in)     :: e
         integer,                     intent(in)     :: integralType
         integer                    , intent(in)     :: num_of_vars
         real(kind=RP)                               :: val(num_of_vars)
         !-local-variables---------------------------------------------
         integer     :: Nel(3)    ! Element polynomial order
         integer     :: i, j, k
         !-------------------------------------------------------------
         
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

#if defined(NAVIERSTOKES)
            case ( VELOCITY )
               
               do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
                  val = val +   wx(i) * wy(j) * wz(k) * e % storage % Q(IRHOU:IRHOW,i,j,k) / e % storage % Q(IRHO,i,j,k) * e % geom % jacobian(i,j,k)
               end do            ; end do           ; end do
               
            case ( MOMENTUM )
               
               do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
                  val = val +   wx(i) * wy(j) * wz(k) * e % storage % Q(IRHOU:IRHOW,i,j,k) * e % geom % jacobian(i,j,k)
               end do            ; end do           ; end do
               
            case ( SOURCE )
               
               do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
                  val = val +   wx(i) * wy(j) * wz(k) * e % storage % S_NS(1:num_of_vars,i,j,k) * e % geom % jacobian(i,j,k)
               end do            ; end do           ; end do

            case ( PSOURCE )
               
               do k = 0, Nel(3)  ; do j = 0, Nel(2) ; do i = 0, Nel(1)
                  val = val +   wx(i) * wy(j) * wz(k) * e % storage % S_NSP(1:num_of_vars,i,j,k) * e % geom % jacobian(i,j,k)
               end do            ; end do           ; end do               
#endif
            case default
               
               write(STD_OUT,'(A,A)') 'VectorVolumeIntegral :: ERROR: Not defined integral type'
               stop 99
               
         end select
         
         end associate
      end function VectorVolumeIntegral_Local
end module VolumeIntegrals
