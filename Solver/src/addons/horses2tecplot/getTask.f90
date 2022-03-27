!
!//////////////////////////////////////////////////////////////////////////////////////////////////
!
!     This module explores the command line flags introduced
!  by the user. Software usage is:
!
!     horses2tecplot *.hmesh *.hsol --output-order=N --output-basis=Gauss/Homogeneous/Chebyshev
!
!  Flags are optional. 
!     * The flag --output-order forces the tecplot file to a fixed polynomial order. It can
!        be specified as:
!           --output-order=N: for isotropic outputs
!           --output-order=Nx,Ny,Nz: for anisotropic outputs
!
!     * The flag --output-basis selects the output basis:
!           --output-basis=Gauss: Gauss-Legendre points.
!           --output-basis=Homogeneous: Equally spaced points.
!
!/////////////////////////////////////////////////////////////////////////////////////////////////
!
#include "Includes.h"
module getTask
   use SMConstants
   use OutputVariables
   implicit none
   
   private
   public   MESH_2_PLT, SOLUTION_2_PLT, UNKNOWN_JOB
   public   EXPORT_GAUSS, EXPORT_HOMOGENEOUS, OUTPUT_VARIABLES_FLAG, MODE_MULTIZONE, MODE_FINITEELM

   public   getTaskType

   integer, parameter   :: UNKNOWN_JOB = 0
   integer, parameter   :: MESH_2_PLT = 1
   integer, parameter   :: SOLUTION_2_PLT = 2
   
   integer, parameter   :: EXPORT_GAUSS = 0
   integer, parameter   :: EXPORT_HOMOGENEOUS = 1
   
   integer, parameter   :: MODE_MULTIZONE = 0
   integer, parameter   :: MODE_FINITEELM = 1

   character(len=*), parameter   :: OUTPUT_ORDER_FLAG="--output-order="
   character(len=*), parameter   :: OUTPUT_MODE_FLAG ="--output-mode="
   character(len=*), parameter   :: OUTPUT_BASIS_FLAG="--output-basis="
   character(len=*), parameter   :: OUTPUT_VARIABLES_FLAG="--output-variables="

   contains

      integer function getTaskType(meshName, no_of_solutions, solutionNames, solutionTypes, fixedOrder, Nout, basis,mode)
         use SolutionFile
         implicit none
         character(len=*),                        intent(out) :: meshName
         integer,                                 intent(out) :: no_of_solutions
         character(len=LINE_LENGTH), allocatable, intent(out) :: solutionNames(:)
         integer, allocatable,                    intent(out) :: solutionTypes(:)
         logical,                                 intent(out) :: fixedOrder
         integer,                                 intent(out) :: Nout(3)
         integer,                                 intent(out) :: basis
         integer,                                 intent(out) :: mode
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: no_of_arguments
         integer     :: i, sol, pos, pos2
         integer     :: io, fid
         logical     :: meshFilePresent
         logical     :: basisPresent
         character(len=LINE_LENGTH) :: auxiliarName, basisName, modeName
!
!        Get number of command arguments         
!        -------------------------------
         no_of_arguments = command_argument_count()
!
!        Exit if no input arguments are specified
!        ----------------------------------------
         if ( no_of_arguments .eq. 0 ) then
            write(STD_OUT,'(A)') "No mesh file and/or solution file specified"
            getTaskType = UNKNOWN_JOB
            return
         end if
!
!        Check if the mesh file is present
!        ---------------------------------
         meshFilePresent = .false.
         do i = 1, no_of_arguments
            call get_command_argument(i, meshName)
            open(newunit=fid, file=trim(meshName), action="read", form="unformatted", access="stream", iostat=io)
            close(fid)
            if ( io .ne. 0 ) cycle

            if ( getSolutionFileType(meshName) .eq. MESH_FILE .or. getSolutionFileType(meshName) .eq. ZONE_MESH_FILE ) then
               meshFilePresent = .true.
               exit
            end if 
         end do
!
!        Exit if the mesh file is not present: always required.
!        ------------------------------------
         if ( .not. meshFilePresent ) then
            write(STD_OUT,'(A)') "Mesh file was not specified"
            errorMessage(STD_OUT)
            getTaskType = UNKNOWN_JOB
            return
         end if
!
!        Check if the solution file is present
!        -------------------------------------
!
!        Loop to get number of files
!        ---------------------------
         no_of_solutions = 0
         do i = 1, no_of_arguments
            call get_command_argument(i, auxiliarName)
            open(newunit=fid, file=trim(auxiliarName), action="read", form="unformatted", access="stream", iostat=io)
            close(fid)
            if ( io .ne. 0 ) cycle

            if ( getSolutionFileType(auxiliarName) .ne. MESH_FILE .and. getSolutionFileType(auxiliarName) .ne. ZONE_MESH_FILE ) then
               no_of_solutions = no_of_solutions + 1 
            end if
         end do
!
!        Loop to get solution file names
!        -------------------------------
         if ( no_of_solutions .ne. 0 ) then
            allocate( solutionNames(no_of_solutions) )
            allocate( solutionTypes(no_of_solutions) )
            
            sol = 0
            do i = 1, no_of_arguments
               call get_command_argument(i, auxiliarName)
               open(newunit=fid, file=trim(auxiliarName), action="read", form="unformatted", access="stream", iostat=io)
               close(fid)
               if ( io .ne. 0 ) cycle

               if ( getSolutionFileType(auxiliarName) .ne. MESH_FILE ) then   
                  sol = sol + 1 
                  solutionNames(sol) = trim(auxiliarName)
                  solutionTypes(sol) = getSolutionFileType(auxiliarName)
               end if
            end do

         end if
!
!        Select the job type: Mesh if no solutions files are present, solution otherwise.
!        -------------------             
         if ( no_of_solutions .ne. 0 ) then
            getTaskType = SOLUTION_2_PLT
   
         else
            getTaskType = MESH_2_PLT

         end if
!
!        Check if fixed polynomials are used: the flag is --output-order=N (or Nx,Ny,Nz)
!        -----------------------------------
         fixedOrder = .false.
         Nout = 0
         do i = 1, no_of_arguments
            call get_command_argument(i, auxiliarName)
            pos = index(trim(auxiliarName), OUTPUT_ORDER_FLAG)
            
            if ( pos .ne. 0 ) then
!
!              Set job as fixed order
!              ----------------------
               fixedOrder = .true.

               if ( index(trim(auxiliarName(pos+1:len_trim(auxiliarName))),",") .eq. 0 ) then
!
!                 Constant polynomial order: --output-order=N
!                 -------------------------
                  read(auxiliarName(pos+len_trim(OUTPUT_ORDER_FLAG):len_trim(auxiliarName)),*) Nout(1)
                  Nout(2:3) = Nout(1)
               else
!
!                 Anisotropic polynomial order: --output-order=Nx,Ny,Nz
!                 ----------------------------
                  pos2 = index(trim(auxiliarName),",")                  ! Read Nx
                  read(auxiliarName(len_trim(OUTPUT_ORDER_FLAG)+1:pos2),*) Nout(1)

                  pos = index(trim(auxiliarName),",",BACK=.true.)
!                                                  Read Ny
                  read(auxiliarName(pos2+1:pos),*) Nout(2)         ! Read Nz  
                  read(auxiliarName(pos+1:len_trim(auxiliarName)),*) Nout(3)
                  
               end if
            end if
         end do
!
!        Get the basis
!        -------------
         basisPresent = .false.
         do i = 1, no_of_arguments
            call get_command_argument(i, auxiliarName)
            pos = index(trim(auxiliarName), OUTPUT_BASIS_FLAG)
            
            if ( pos .ne. 0 ) then
               read(auxiliarName(pos+len_trim(OUTPUT_BASIS_FLAG):len_trim(auxiliarName)),*) basisName
               basisPresent = .true.
               exit
            end if
         end do
!
!        Get the mode
!        ------------
         mode = MODE_MULTIZONE ! Default
         do i = 1, no_of_arguments
            call get_command_argument(i, auxiliarName)
            pos = index(trim(auxiliarName), OUTPUT_MODE_FLAG)
            
            if ( pos .ne. 0 ) then
               read(auxiliarName(pos+len_trim(OUTPUT_MODE_FLAG):len_trim(auxiliarName)),*) modeName
               select case(modeName)
               case ("multizone")
                  mode = MODE_MULTIZONE
               case ("FE")
                  mode = MODE_FINITEELM
               end select
               exit
            end if
         end do
!
!        *******************************************************
!        If the basis is not present, the choice depends
!        on whether a new polynomial order is given.
!           * default polynomial: default basis is Gauss.
!           * given polynoamial: default basis is Homogeneous.
!        *******************************************************
!
         if ( basisPresent ) then
      
            select case (trim(basisName))
            case ("Homogeneous")
               basis = EXPORT_HOMOGENEOUS
   
            case ("Gauss")
               basis = EXPORT_GAUSS
   
            case default
               print*, 'Basis "',trim(basisName),'" not recognized.'
               print*, "Options:"
               print*, "   * Homogeneous"
               print*, "   * Gauss"
               stop
      
            end select

         else

            if ( fixedOrder ) then
!
!              Given polynomial: homogeneous
!              -----------------------------
               basis = EXPORT_HOMOGENEOUS

            else
!
!              Default polynomial: Gauss
!              -------------------------
               basis = EXPORT_GAUSS

            end if

         end if

      end function getTaskType
end module getTask
