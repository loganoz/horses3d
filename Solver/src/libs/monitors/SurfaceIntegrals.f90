module SurfaceIntegrals
   use SMConstants
   use Physics
   use NodalStorageClass
   use HexMeshClass
   use ProlongToFacesProcedures
   
   private
   public   SURFACE, TOTAL_FORCE, PRESSURE_FORCE, VISCOUS_FORCE, MASS_FLOW, FLOW
   public   ScalarSurfaceIntegral, VectorSurfaceIntegral

   integer, parameter   :: SURFACE = 1
   integer, parameter   :: TOTAL_FORCE = 2
   integer, parameter   :: PRESSURE_FORCE = 3
   integer, parameter   :: VISCOUS_FORCE = 4
   integer, parameter   :: MASS_FLOW = 5
   integer, parameter   :: FLOW = 6
   integer, parameter   :: USER_DEFINED = 99
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
      function ScalarSurfaceIntegral(mesh, spA, zoneID, integralType) result(val)
!
!        -----------------------------------------------------------
!           This function computes scalar integrals, that is, those
!           in the form:
!                 val = \int \vec{v}·\vec{n}dS
!           Implemented integrals are:
!              * Surface: computes the zone surface.
!              * Mass flow: computes the mass flow across the zone.
!              * Flow: computes the volumetric flow across the zone.
!        -----------------------------------------------------------
!
         implicit none
         class(HexMesh),      intent(in)  :: mesh
         class(NodalStorage), intent(in)  :: spA(0:,0:,0:)
         integer,             intent(in)  :: zoneID
         integer,             intent(in)  :: integralType
         real(kind=RP)                    :: val
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: zonefID, fID, eID 
!
!        Initialization
!        --------------            
         val = 0.0_RP
!
!        Loop the zone to get faces and elements
!        ---------------------------------------
         do zonefID = 1, mesh % zones(zoneID) % no_of_faces         
!
!           Face global ID
!           --------------
            fID = mesh % zones(zoneID) % faces(zonefID)
!
!           Associated element (only boundary faces are expected)
!           -----------------------------------------------------
            eID = mesh % faces(fID) % elementIDs(1)
!
!           Compute the integral
!           --------------------
            val = val + ScalarSurfaceIntegral_Face(mesh % elements(eID), &
                                                                    spA, &
                                     mesh % faces(fID) % elementSide(1), &
                                                           integralType    )

         end do

      end function ScalarSurfaceIntegral

      function ScalarSurfaceIntegral_Face(e, spA, elSide, integralType) result(val)
         implicit none
         class(Element),      target, intent(in)     :: e
         class(NodalStorage), target, intent(in)     :: spA(0:,0:,0:)
         integer,                     intent(in)     :: elSide
         integer,                     intent(in)     :: integralType
         real(kind=RP)                               :: val
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: Nf(2)     ! Face polynomial order
         integer     :: Nel(3)    ! Element polynomial order
         integer     :: i, j      ! Face indices
         real(kind=RP), pointer  :: wx(:), wy(:)   ! Face quadrature weights
         real(kind=RP), pointer  :: ds(:,:)
         real(kind=RP), pointer  :: Qb(:)
         real(kind=RP), pointer  :: n(:,:,:)

         Nel = e % Nxyz
         Nf  = e % Nxyz(axisMap(:,elSide))
!
!        Get the weights
!        ---------------
         select case ( elSide )
            case(1,2)
               wx => spA(Nel(1),Nel(2),Nel(3)) % wx
               wy => spA(Nel(1),Nel(2),Nel(3)) % wz
            case(3,5)
               wx => spA(Nel(1),Nel(2),Nel(3)) % wx
               wy => spA(Nel(1),Nel(2),Nel(3)) % wy
            case(4,6)
               wx => spA(Nel(1),Nel(2),Nel(3)) % wy
               wy => spA(Nel(1),Nel(2),Nel(3)) % wz
         end select
!
!        Get the surface Jacobian and normal vector
!        ------------------------------------------
         ds(0:,0:)   => e % geom % scal(:,:,elSide)
         n(1:,0:,0:) => e % geom % normal(:,:,:,elSide)
!
!        Prolong variables to faces
!        --------------------------            
         call ProlongToFaces(e, spA(Nel(1),Nel(2),Nel(3)) )
         if (flowIsNavierStokes) call ProlongGradientToFaces(e, spA(Nel(1),Nel(2),Nel(3)))
!
!        Initialization
!        --------------
         val = 0.0_RP
!
!        Perform the numerical integration
!        ---------------------------------
         select case ( integralType )
   
         case ( SURFACE )
!
!           **********************************
!           Computes the surface integral
!              val = \int dS
!           **********************************
!
            do j = 0, Nf(2) ;    do i = 0, Nf(1)
               val = val + wx(i) * wy(j) * ds(i,j)
            end do          ;    end do

         case ( MASS_FLOW )
!
!           ***********************************
!           Computes the mass-flow integral
!              I = \int rho \vec{v}·\vec{n}dS         
!           ***********************************
!
            do j = 0, Nf(2) ;    do i = 0, Nf(1)
!
!              Get the state vector
!              --------------------
               Qb => e % Qb(:,i,j,elSide)
!
!              Compute the integral
!              --------------------
               val = val +  (   Qb(IRHOU) * n(1,i,j)  &
                          + Qb(IRHOV) * n(2,i,j)  &
                          + Qb(IRHOW) * n(3,i,j) ) &
                       * wx(i) * wy(j) * ds(i,j)

            end do          ;    end do

         case ( FLOW ) 
!
!           ***********************************
!           Computes the flow integral
!              val = \int \vec{v}·\vec{n}dS         
!           ***********************************
!
            do j = 0, Nf(2) ;    do i = 0, Nf(1)
!
!              Get the state vector
!              --------------------
               Qb => e % Qb(:,i,j,elSide)
!
!              Compute the integral
!              --------------------
               val = val + (1.0_RP / Qb(IRHO))*(   Qb(IRHOU) * n(1,i,j)  &
                                             + Qb(IRHOV) * n(2,i,j)  &
                                             + Qb(IRHOW) * n(3,i,j) ) &
                                          * wx(i) * wy(j) * ds(i,j) 
            end do          ;    end do

         case ( USER_DEFINED )   ! TODO
         end select

      end function ScalarSurfaceIntegral_Face
!
!////////////////////////////////////////////////////////////////////////////////////////
!
!           VECTOR INTEGRALS PROCEDURES
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      function VectorSurfaceIntegral(mesh, spA, zoneID, integralType) result(val)
!
!        -----------------------------------------------------------
!           This function computes scalar integrals, that is, those
!           in the form:
!                 val = \int \vec{v}·\vec{n}dS
!           Implemented integrals are:
!              * Surface: computes the zone surface.
!              * Mass flow: computes the mass flow across the zone.
!              * Flow: computes the volumetric flow across the zone.
!        -----------------------------------------------------------
!
         implicit none
         class(HexMesh),      intent(in)  :: mesh
         class(NodalStorage), intent(in)  :: spA(0:,0:,0:)
         integer,             intent(in)  :: zoneID
         integer,             intent(in)  :: integralType
         real(kind=RP)                    :: val(NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: zonefID, fID, eID 
!
!        Initialization
!        --------------            
         val = 0.0_RP
!
!        Loop the zone to get faces and elements
!        ---------------------------------------
         do zonefID = 1, mesh % zones(zoneID) % no_of_faces         
!
!           Face global ID
!           --------------
            fID = mesh % zones(zoneID) % faces(zonefID)
!
!           Associated element (only boundary faces are expected)
!           -----------------------------------------------------
            eID = mesh % faces(fID) % elementIDs(1)
!
!           Compute the integral
!           --------------------
            val = val + VectorSurfaceIntegral_Face(mesh % elements(eID), &
                                                                    spA, &
                                     mesh % faces(fID) % elementSide(1), &
                                                           integralType    )

         end do

      end function VectorSurfaceIntegral

      function VectorSurfaceIntegral_Face(e, spA, elSide, integralType) result(val)
         implicit none
         class(Element),      target, intent(in)     :: e
         class(NodalStorage), target, intent(in)     :: spA(0:,0:,0:)
         integer,                     intent(in)     :: elSide
         integer,                     intent(in)     :: integralType
         real(kind=RP)                               :: val(NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: Nf(2)     ! Face polynomial order
         integer     :: Nel(3)    ! Element polynomial order
         integer     :: i, j      ! Face indices
         real(kind=RP), pointer  :: wx(:), wy(:)   ! Face quadrature weights
         real(kind=RP), pointer  :: ds(:,:)
         real(kind=RP), pointer  :: Qb(:)
         real(kind=RP), pointer  :: U_xb(:), U_yb(:), U_zb(:)
         real(kind=RP), pointer  :: n(:,:,:)
         real(kind=RP)           :: p, tau(NDIM,NDIM)

         Nel = e % Nxyz
         Nf  = e % Nxyz(axisMap(:,elSide))
!
!        Get the weights
!        ---------------
         select case ( elSide )
            case(1,2)
               wx => spA(Nel(1),Nel(2),Nel(3)) % wx
               wy => spA(Nel(1),Nel(2),Nel(3)) % wz
            case(3,5)
               wx => spA(Nel(1),Nel(2),Nel(3)) % wx
               wy => spA(Nel(1),Nel(2),Nel(3)) % wy
            case(4,6)
               wx => spA(Nel(1),Nel(2),Nel(3)) % wy
               wy => spA(Nel(1),Nel(2),Nel(3)) % wz
         end select
!
!        Get the surface Jacobian and normal vector
!        ------------------------------------------
         ds(0:,0:) => e % geom % scal(:,:,elSide)
         n(1:,0:,0:)  => e % geom % normal(:,:,:,elSide)
!
!        Prolong variables to faces
!        --------------------------            
         call ProlongToFaces(e, spA(Nel(1),Nel(2),Nel(3)) )
         if (flowIsNavierStokes) call ProlongGradientToFaces(e, spA(Nel(1),Nel(2),Nel(3)))
!
!        Initialization
!        --------------
         val = 0.0_RP
!
!        Perform the numerical integration
!        ---------------------------------
         select case ( integralType )
   
         case ( SURFACE )
!
!           **********************************
!           Computes the surface integral
!              val = \int \vec{n} dS
!           **********************************
!
            do j = 0, Nf(2) ;    do i = 0, Nf(1)
               val = val + wx(i) * wy(j) * ds(i,j) * n(:,i,j)
            end do          ;    end do

         case ( TOTAL_FORCE )
!
!           ************************************************
!           Computes the total force experienced by the zone 
!              F = \int p \vec{n}ds - \int tau'·\vec{n}ds 
!           ************************************************
!
            do j = 0, Nf(2) ;    do i = 0, Nf(1)
!
!              Get the state vector
!              --------------------
               Qb => e % Qb(:,i,j,elSide)
               U_xb => e % U_xb(:,i,j,elSide)
               U_yb => e % U_yb(:,i,j,elSide)
               U_zb => e % U_zb(:,i,j,elSide)
               
!
!              Compute the integral
!              --------------------
               p = Pressure(Qb)
               tau = getStressTensor(Qb,U_xb,U_yb,U_zb)

               val = val + ( p * n(:,i,j) - matmul(tau,n(:,i,j)) ) * ds(i,j) * wx(i) * wx(j)

            end do          ;    end do

         case ( PRESSURE_FORCE ) 
!
!           ****************************************************
!           Computes the pressure forces experienced by the zone 
!              F = \int p \vec{n}ds 
!           ****************************************************
!
            do j = 0, Nf(2) ;    do i = 0, Nf(1)
!
!              Get the state vector
!              --------------------
               Qb => e % Qb(:,i,j,elSide)
!
!              Compute the integral
!              --------------------
               p = Pressure(Qb)

               val = val + ( p * n(:,i,j) ) * ds(i,j) * wx(i) * wx(j)

            end do          ;    end do

         case ( VISCOUS_FORCE ) 
!
!           ************************************************
!           Computes the total force experienced by the zone 
!              F = \int p \vec{n}ds - \int tau'·\vec{n}ds 
!           ************************************************
!
            do j = 0, Nf(2) ;    do i = 0, Nf(1)
!
!              Get the state vector
!              --------------------
               Qb => e % Qb(:,i,j,elSide)
               U_xb => e % U_xb(:,i,j,elSide)
               U_yb => e % U_yb(:,i,j,elSide)
               U_zb => e % U_zb(:,i,j,elSide)
               
!
!              Compute the integral
!              --------------------
               tau = getStressTensor(Qb,U_xb,U_yb,U_zb)
               val = val - matmul(tau,n(:,i,j))  * ds(i,j) * wx(i) * wx(j)

            end do          ;    end do

         case ( USER_DEFINED )   ! TODO

         end select

      end function VectorSurfaceIntegral_Face

end module SurfaceIntegrals
