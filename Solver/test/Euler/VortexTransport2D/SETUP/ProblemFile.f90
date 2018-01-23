!
!//////////////////////////////////////////////////////
!
!   @File:    ProblemFile.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Jan 23 16:27:56 2018
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
!
!////////////////////////////////////////////////////////////////////////
!
!      ProblemFile.f90
!      Created: June 26, 2015 at 8:47 AM 
!      By: David Kopriva  
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
         SUBROUTINE UserDefinedFinalSetup(mesh , thermodynamics_, &
                                                 dimensionless_, &
                                                     refValues_ )
!
!           ----------------------------------------------------------------------
!           Called after the mesh is read in to allow mesh related initializations
!           or memory allocations.
!           ----------------------------------------------------------------------
!
            USE HexMeshClass
            use PhysicsStorage
            IMPLICIT NONE
            CLASS(HexMesh)                      :: mesh
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_

         END SUBROUTINE UserDefinedFinalSetup
!
!//////////////////////////////////////////////////////////////////////// 
! 
#if defined(NAVIERSTOKES)
         subroutine userdefinedinitialcondition(mesh, thermodynamics_, &
                                                      dimensionless_, &
                                                          refvalues_  )
!
!           ------------------------------------------------
!           called to set the initial condition for the flow
!              - by default it sets an uniform initial
!                 condition.
!           ------------------------------------------------
!
            use smconstants
            use physicsstorage
            use hexmeshclass
            implicit none
            class(hexmesh)                      :: mesh
            type(thermodynamics_t), intent(in)  :: thermodynamics_
            type(dimensionless_t),  intent(in)  :: dimensionless_
            type(refvalues_t),      intent(in)  :: refvalues_
!
!           ---------------
!           local variables
!           ---------------
!
            integer        :: eid, i, j, k
            real(kind=rp)  :: q(n_eqn)

            real(kind=RP), parameter :: R = 0.1_RP                  ! This is the "vortex radius" (dimensionless)
            real(kind=RP), parameter :: Beta = 0.05_RP              ! This is the "vortex strength"
            real(kind=RP), parameter :: XC = 0.0_RP                 ! Vortex X position (in dimensionless coordinates)
            real(kind=RP), parameter :: ZC = 0.0_RP                 ! Vortex Y position (in dimensionless coordinates)
            real(kind=RP), parameter :: AngleOfAttack = 0.0_RP
            real(kind=RP)            :: r2 , rho , u , w , T, p
            real(kind=RP)            :: x(NDIM)
         
            associate ( gamma => Thermodynamics_ % Gamma , Mach => Dimensionless_ % Mach , cv => Dimensionless_ % cv )

            do eid = 1, mesh % no_of_elements
               associate( nx => mesh % elements(eid) % nxyz(1), &
                          ny => mesh % elements(eid) % nxyz(2), &
                          nz => mesh % elements(eid) % nxyz(3) )
               do k = 0, nz;  do j = 0, ny;  do i = 0, nx 
                  x = mesh % elements(eID) % geom % x(:,i,j,k)
            r2 = ((x(1) - XC)*(x(1) - XC) + (x(3) - ZC)*(x(3) - ZC)) / (R*R)
         
            u =  (cos(AngleOfAttack) - Beta * (x(3) - ZC) / R * exp(-0.5_RP * r2))
            w =  (sin(AngleOfAttack) + Beta * (x(1) - XC) / R * exp(-0.5_RP * r2))
             
            T = (1.0_RP - gamma * Mach * Mach * beta * beta / (2.0_RP * Dimensionless_ % cp) * exp(-r2) )
            rho = (T)**( dimensionless_ % cv)
            p  = rho * T / dimensionless_ % gammaM2
         
            q(1) = rho
            q(2) = rho * u
            q(3) = 0.0_RP
            q(4) = rho * w
            q(5) = p/(gamma - 1.0_RP) + 0.5_RP * q(1)*(u**2 + w**2)
      
                  mesh % elements(eid) % storage % q(:,i,j,k) = q 
               end do;        end do;        end do
               end associate
            end do
            end associate

         end subroutine userdefinedinitialcondition
#elif defined(CAHNHILLIARD)
         subroutine userdefinedinitialcondition(mesh, thermodynamics_, &
                                                      dimensionless_, &
                                                          refvalues_  )
            use smconstants
            use physicsstorage
            use hexmeshclass
            implicit none
            class(hexmesh)                      :: mesh
            type(thermodynamics_t), intent(in)  :: thermodynamics_
            type(dimensionless_t),  intent(in)  :: dimensionless_
            type(refvalues_t),      intent(in)  :: refvalues_
            
         end subroutine userdefinedinitialcondition
#endif
         subroutine UserDefinedState1(x, t, nHat, Q, thermodynamics_, dimensionless_, refValues_)
!
!           -------------------------------------------------
!           Used to define an user defined boundary condition
!           -------------------------------------------------
!
            use SMConstants
            use PhysicsStorage
            implicit none
            real(kind=RP), intent(in)     :: x(NDIM)
            real(kind=RP), intent(in)     :: t
            real(kind=RP), intent(in)     :: nHat(NDIM)
            real(kind=RP), intent(inout)  :: Q(N_EQN)
            type(Thermodynamics_t),    intent(in)  :: thermodynamics_
            type(Dimensionless_t),     intent(in)  :: dimensionless_
            type(RefValues_t),         intent(in)  :: refValues_
         end subroutine UserDefinedState1

         subroutine UserDefinedNeumann(x, t, nHat, U_x, U_y, U_z)
!
!           --------------------------------------------------------
!           Used to define a Neumann user defined boundary condition
!           --------------------------------------------------------
!
            use SMConstants
            use PhysicsStorage
            implicit none
            real(kind=RP), intent(in)     :: x(NDIM)
            real(kind=RP), intent(in)     :: t
            real(kind=RP), intent(in)     :: nHat(NDIM)
            real(kind=RP), intent(inout)  :: U_x(N_GRAD_EQN)
            real(kind=RP), intent(inout)  :: U_y(N_GRAD_EQN)
            real(kind=RP), intent(inout)  :: U_z(N_GRAD_EQN)
         end subroutine UserDefinedNeumann

!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedPeriodicOperation(mesh, time, Monitors)
!
!           ----------------------------------------------------------
!           Called at the output interval to allow periodic operations
!           to be performed
!           ----------------------------------------------------------
!
            USE HexMeshClass
#if defined(NAVIERSTOKES)
            use MonitorsClass
#endif
            IMPLICIT NONE
            CLASS(HexMesh)               :: mesh
            REAL(KIND=RP)                :: time
#if defined(NAVIERSTOKES)
            type(Monitor_t), intent(in) :: monitors
#else
            logical, intent(in) :: monitors
#endif
            
         END SUBROUTINE UserDefinedPeriodicOperation
!
!//////////////////////////////////////////////////////////////////////// 
! 
         subroutine UserDefinedSourceTerm(mesh, time, thermodynamics_, dimensionless_, refValues_)
!
!           --------------------------------------------
!           Called to apply source terms to the equation
!           --------------------------------------------
!
            USE HexMeshClass
            use PhysicsStorage
            IMPLICIT NONE
            CLASS(HexMesh)                        :: mesh
            REAL(KIND=RP)                         :: time
            type(Thermodynamics_t),    intent(in) :: thermodynamics_
            type(Dimensionless_t),     intent(in) :: dimensionless_
            type(RefValues_t),         intent(in) :: refValues_
!
!           ---------------
!           Local variables
!           ---------------
!
            integer  :: i, j, k, eID
!
!           Usage example (by default no source terms are added)
!           ----------------------------------------------------
!           do eID = 1, mesh % no_of_elements
!              associate ( e => mesh % elements(eID) )
!              do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
!                 e % QDot(:,i,j,k) = e % QDot(:,i,j,k) + Source(:)
!              end do                  ; end do                ; end do
!           end do
   
         end subroutine UserDefinedSourceTerm
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalize(mesh, time, iter, maxResidual, thermodynamics_, &
                                                    dimensionless_, &
                                                        refValues_, &  
                                                          monitors, &
                                                       elapsedTime, &
                                                           CPUTime   )
!
!           --------------------------------------------------------
!           Called after the solution computed to allow, for example
!           error tests to be performed
!           --------------------------------------------------------
!
            USE HexMeshClass
            use PhysicsStorage
#if defined(NAVIERSTOKES)
            use MonitorsClass
#endif
            IMPLICIT NONE
            CLASS(HexMesh)                        :: mesh
            REAL(KIND=RP)                         :: time
            integer                               :: iter
            real(kind=RP)                         :: maxResidual
            type(Thermodynamics_t),    intent(in) :: thermodynamics_
            type(Dimensionless_t),     intent(in) :: dimensionless_
            type(RefValues_t),         intent(in) :: refValues_
#if defined(NAVIERSTOKES)
            type(Monitor_t),          intent(in) :: monitors
#else
            logical, intent(in)  :: monitors
#endif
            real(kind=RP),             intent(in)  :: elapsedTime
            real(kind=RP),             intent(in)  :: CPUTime

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
      
