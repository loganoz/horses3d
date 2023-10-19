!
!////////////////////////////////////////////////////////////////////////
!
!      The Problem File contains user defined procedures
!      that are used to "personalize" i.e. define a specific
!      problem to be solved. These procedures include initial conditions,
!      exact solutions (e.g. for tests), etc. and allow modifications 
!      without having to modify the main code.
!
!      The procedures, *even if empty* that must be defined are
!
!      UserDefinedSetUp
!      UserDefinedInitialCondition(mesh)
!      UserDefinedPeriodicOperation(mesh)
!      UserDefinedFinalize(mesh)
!      UserDefinedTermination
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedStartup
!
!        --------------------------------
!        Called before any other routines
!        --------------------------------
!
            IMPLICIT NONE  
         END SUBROUTINE UserDefinedStartup
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalSetup(mesh &
#if defined(NAVIERSTOKES)
                                        , thermodynamics_   &
                                        , dimensionless_    &
                                        , refValues_        &
#endif
#if defined(CAHNHILLIARD)
                                        , multiphase_       &
#endif
)
!
!           ----------------------------------------------------------------------
!           Called after the mesh is read in to allow mesh related initializations
!           or memory allocations.
!           ----------------------------------------------------------------------
!
            use SMConstants
            USE HexMeshClass
            use FluidData
            IMPLICIT NONE
            class(HexMesh)                      :: mesh
#if defined(NAVIERSTOKES)
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
#endif
#if defined(CAHNHILLIARD)
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
         END SUBROUTINE UserDefinedFinalSetup
!
!//////////////////////////////////////////////////////////////////////// 
! 
         subroutine UserDefinedInitialCondition(mesh              &
#if defined(NAVIERSTOKES)
                                              , thermodynamics_   &
                                              , dimensionless_    &
                                              , refvalues_        &
#endif
#if defined(CAHNHILLIARD)
                                              , multiphase_       &
#endif

)
!
!           ------------------------------------------------
!           called to set the initial condition for the flow
!              - by default it sets an uniform initial
!                 condition.
!           ------------------------------------------------
!
            use SMConstants
            use PhysicsStorage
            use HexMeshClass
            use FluidData
            implicit none
            class(HexMesh)                      :: mesh
#if defined(NAVIERSTOKES)
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refvalues_
#endif
#if defined(CAHNHILLIARD)
            type(Multiphase_t),     intent(in)  :: multiphase_
#endif
!
!           ---------------
!           local variables
!           ---------------
!
            integer        :: eid, i, j, k
            real(kind=RP)  :: qq, u, v, w, p, x(NDIM)
#if defined(NAVIERSTOKES)
            real(kind=RP)  :: Q(NCONS), phi, theta
#endif

#if defined(NAVIERSTOKES)
            associate ( gammaM2 => dimensionless_ % gammaM2, &
                        gamma   => thermodynamics_ % gamma )
            theta = refvalues_ % aoatheta*(pi/180.0_RP)
            phi   = refvalues_ % aoaphi*(pi/180.0_RP)
      
            do eid = 1, mesh % no_of_elements
               associate( nx => mesh % elements(eid) % nxyz(1), &
                          ny => mesh % elements(eid) % nxyz(2), &
                          nz => mesh % elements(eid) % nxyz(3) )
               do k = 0, nz;  do j = 0, ny;  do i = 0, nx 
                  qq = 1.0_RP
                  u  = qq*cos(theta)*cos(phi)
                  v  = qq*sin(theta)*cos(phi)
                  w  = qq*sin(phi)
      
                  q(1) = 1.0_RP
                  p    = 1.0_RP/(gammam2)
                  q(2) = q(1)*u
                  q(3) = q(1)*v
                  q(4) = q(1)*w
                  q(5) = p/(gamma - 1._RP) + 0.5_RP*q(1)*(u**2 + v**2 + w**2)

                  mesh % elements(eid) % storage % q(:,i,j,k) = q 

               end do;        end do;        end do
               end associate

            end do
            end associate
#endif

#if defined(CAHNHILLIARD)
            do eID = 1, mesh % no_of_elements
               associate(e => mesh % elements(eID))
               do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2)    ; do i = 0, e % Nxyz(1)
                  x = e % geom % x(:,i,j,k)

                  if ( (abs(x(1)-20.0_RP) .le. 7.0_RP) .and. (abs(x(3)-20.0_RP) .le. 7.0_RP)) then
                     e % storage % c(1,i,j,k) = -1.0_RP
                  else
                     e % storage % c(1,i,j,k) =  1.0_RP
                  end if

               end do                  ; end do                   ; end do
               end associate
            end do
#endif
            
         end subroutine UserDefinedInitialCondition
#if defined(NAVIERSTOKES)
         subroutine UserDefinedState1(x, t, nHat, Q, thermodynamics_, dimensionless_, refValues_)
!
!           -------------------------------------------------
!           Used to define an user defined boundary condition
!           -------------------------------------------------
!
            use SMConstants
            use PhysicsStorage
            use FluidData
            implicit none
            real(kind=RP), intent(in)     :: x(NDIM)
            real(kind=RP), intent(in)     :: t
            real(kind=RP), intent(in)     :: nHat(NDIM)
            real(kind=RP), intent(inout)  :: Q(NCONS)
            type(Thermodynamics_t),    intent(in)  :: thermodynamics_
            type(Dimensionless_t),     intent(in)  :: dimensionless_
            type(RefValues_t),         intent(in)  :: refValues_
         end subroutine UserDefinedState1

         subroutine UserDefinedGradVars1(x, t, nHat, Q, U, GetGradients, thermodynamics_, dimensionless_, refValues_)
            use SMConstants
            use PhysicsStorage
            use FluidData
            use VariableConversion, only: GetGradientValues_f
            implicit none
            real(kind=RP), intent(in)          :: x(NDIM)
            real(kind=RP), intent(in)          :: t
            real(kind=RP), intent(in)          :: nHat(NDIM)
            real(kind=RP), intent(in)          :: Q(NCONS)
            real(kind=RP), intent(inout)       :: U(NGRAD)
            procedure(GetGradientValues_f)     :: GetGradients
            type(Thermodynamics_t), intent(in) :: thermodynamics_
            type(Dimensionless_t),  intent(in) :: dimensionless_
            type(RefValues_t),      intent(in) :: refValues_
         end subroutine UserDefinedGradVars1

         subroutine UserDefinedNeumann1(x, t, nHat, U_x, U_y, U_z)
!
!           --------------------------------------------------------
!           Used to define a Neumann user defined boundary condition
!           --------------------------------------------------------
!
            use SMConstants
            use PhysicsStorage
            use FluidData
            implicit none
            real(kind=RP), intent(in)     :: x(NDIM)
            real(kind=RP), intent(in)     :: t
            real(kind=RP), intent(in)     :: nHat(NDIM)
            real(kind=RP), intent(inout)  :: U_x(NGRAD)
            real(kind=RP), intent(inout)  :: U_y(NGRAD)
            real(kind=RP), intent(inout)  :: U_z(NGRAD)
         end subroutine UserDefinedNeumann1
#endif
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedPeriodicOperation(mesh, time, dt, Monitors)
!
!           ----------------------------------------------------------
!           Called at the output interval to allow periodic operations
!           to be performed
!           ----------------------------------------------------------
!
            USE HexMeshClass
            use SMConstants
            use MonitorsClass
            IMPLICIT NONE
            class(HexMesh)               :: mesh
            REAL(KIND=RP)                :: time
            type(Monitor_t), intent(in) :: monitors
            integer     :: eID

#if defined(CAHNHILLIARD)
            do eID = 1, mesh % no_of_elements
               if ( any(isnan(mesh % elements(eID) % storage % c) )) then
                  print*, "NAN!!!!"
                  error stop
               end if
            end do
#endif
            
         END SUBROUTINE UserDefinedPeriodicOperation
!
!//////////////////////////////////////////////////////////////////////// 
! 
<<<<<<< HEAD
         subroutine UserDefinedSourceTerm(mesh, time &
#if defined(NAVIERSTOKES)
                                        , thermodynamics_, dimensionless_, refValues_ &
#endif
#if defined(CAHNHILLIARD)
                                        , multiphase_ &
#endif
)
=======
         subroutine UserDefinedSourceTerm(x, time, S, thermodynamics_, dimensionless_, refValues_)
>>>>>>> master
!
!           --------------------------------------------
!           Called to apply source terms to the equation
!           --------------------------------------------
!
            use SMConstants
            USE HexMeshClass
            use PhysicsStorage
            use FluidData
            IMPLICIT NONE
<<<<<<< HEAD
            class(HexMesh)                        :: mesh
            REAL(KIND=RP)                         :: time
#if defined(NAVIERSTOKES)
            type(Thermodynamics_t),    intent(in) :: thermodynamics_
            type(Dimensionless_t),     intent(in) :: dimensionless_
            type(RefValues_t),         intent(in) :: refValues_
#endif
#if defined(CAHNHILLIARD)
            type(Multiphase_t),        intent(in) :: multiphase_
#endif
!
!           ---------------
!           Local variables
!           ---------------
!
            integer  :: i, j, k, eID
=======
            real(kind=RP),             intent(in)  :: x(NDIM)
            real(kind=RP),             intent(in)  :: time
            real(kind=RP),             intent(inout) :: S(NCONS)
            type(Thermodynamics_t),    intent(in)  :: thermodynamics_
            type(Dimensionless_t),     intent(in)  :: dimensionless_
            type(RefValues_t),         intent(in)  :: refValues_
>>>>>>> master
!
!           Usage example
!           -------------
!           S(:) = x(1) + x(2) + x(3) + time
   
         end subroutine UserDefinedSourceTerm
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalize(mesh, time, iter, maxResidual  &
#if defined(NAVIERSTOKES)
                                        , thermodynamics_, dimensionless_, refValues_ &
#endif
#if defined(CAHNHILLIARD)
                                        , multiphase_ &
#endif
                                        , monitors, elapsedTime, CPUTime )
!
!           --------------------------------------------------------
!           Called after the solution computed to allow, for example
!           error tests to be performed
!           --------------------------------------------------------
!
            use SMConstants
            USE HexMeshClass
            use PhysicsStorage
            use MonitorsClass
            use FluidData
            IMPLICIT NONE
            class(HexMesh)                        :: mesh
            REAL(KIND=RP)                         :: time
            integer                               :: iter
            real(kind=RP)                         :: maxResidual
#if defined(NAVIERSTOKES)
            type(Thermodynamics_t),    intent(in) :: thermodynamics_
            type(Dimensionless_t),     intent(in) :: dimensionless_
            type(RefValues_t),         intent(in) :: refValues_
#endif
#if defined(CAHNHILLIARD)
            type(Multiphase_t),        intent(in) :: multiphase_
#endif
            type(Monitor_t),           intent(in) :: monitors
            real(kind=RP),             intent(in) :: elapsedTime
            real(kind=RP),             intent(in) :: CPUTime

         END SUBROUTINE UserDefinedFinalize
!
!//////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE UserDefinedTermination
!
!        -----------------------------------------------
!        Called at the the end of the main driver after 
!        everything else is done.
!        -----------------------------------------------
!
         IMPLICIT NONE  
      END SUBROUTINE UserDefinedTermination
      