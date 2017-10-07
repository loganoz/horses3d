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
            CLASS(HexMesh)             :: mesh
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
         END SUBROUTINE UserDefinedFinalSetup
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedInitialCondition(mesh, thermodynamics_, &
                                                      dimensionless_, &
                                                          refValues_  )
!
!           ------------------------------------------------
!           Called to set the initial condition for the flow
!              - By default it sets an uniform initial
!                 condition.
!           ------------------------------------------------
!
            USE SMConstants
            use PhysicsStorage
            use HexMeshClass
            implicit none
            class(HexMesh)                        :: mesh
            type(Thermodynamics_t), intent(in)  :: thermodynamics_
            type(Dimensionless_t),  intent(in)  :: dimensionless_
            type(RefValues_t),      intent(in)  :: refValues_
!
!           ---------------
!           Local variables
!           ---------------
!
            integer        :: eID, i, j, k
            real(kind=RP)  :: qq, u, v, w, p
            real(kind=RP)  :: Q(N_EQN), phi, theta

            associate ( gammaM2 => dimensionless_ % gammaM2, &
                        gamma => thermodynamics_ % gamma )
            theta = refValues_ % AOATheta*(PI/180.0_RP)
            phi   = refValues_ % AOAPhi*(PI/180.0_RP)
      
            do eID = 1, mesh % no_of_elements
               associate( Nx => mesh % elements(eID) % Nxyz(1), &
                          Ny => mesh % elements(eID) % Nxyz(2), &
                          Nz => mesh % elements(eID) % Nxyz(3) )
               do k = 0, Nz;  do j = 0, Ny;  do i = 0, Nx 
                  qq = 1.0_RP
                  u  = qq*cos(theta)*COS(phi)
                  v  = qq*sin(theta)*COS(phi)
                  w  = qq*SIN(phi)
      
                  Q(1) = 1.0_RP
                  p    = 1.0_RP/(gammaM2)
                  Q(2) = Q(1)*u
                  Q(3) = Q(1)*v
                  Q(4) = Q(1)*w
                  Q(5) = p/(gamma - 1._RP) + 0.5_RP*Q(1)*(u**2 + v**2 + w**2)

                  mesh % elements(eID) % storage % Q(i,j,k,:) = Q 
               end do;        end do;        end do
               end associate
            end do

            end associate
            
         END SUBROUTINE UserDefinedInitialCondition

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
         SUBROUTINE UserDefinedPeriodicOperation(mesh, time)
!
!           ----------------------------------------------------------
!           Called at the output interval to allow periodic operations
!           to be performed
!           ----------------------------------------------------------
!
            USE HexMeshClass
            IMPLICIT NONE
            CLASS(HexMesh)  :: mesh
            REAL(KIND=RP) :: time
            
         END SUBROUTINE UserDefinedPeriodicOperation
!
!//////////////////////////////////////////////////////////////////////// 
! 
         SUBROUTINE UserDefinedFinalize(mesh, time, iter, maxResidual, thermodynamics_, &
                                                    dimensionless_, &
                                                        refValues_   )
!
!           --------------------------------------------------------
!           Called after the solution computed to allow, for example
!           error tests to be performed
!           --------------------------------------------------------
!
            USE HexMeshClass
            use PhysicsStorage
            IMPLICIT NONE
            CLASS(HexMesh)                        :: mesh
            REAL(KIND=RP)                         :: time
            integer                               :: iter
            real(kind=RP)                         :: maxResidual
            type(Thermodynamics_t),    intent(in) :: thermodynamics_
            type(Dimensionless_t),     intent(in) :: dimensionless_
            type(RefValues_t),         intent(in) :: refValues_

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
      
