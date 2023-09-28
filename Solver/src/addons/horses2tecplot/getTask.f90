!
!//////////////////////////////////////////////////////////////////////////////////////////////////
!
!     This module explores the command line flags introduced by the user. Software usage is:
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
   character(len=*), parameter   :: BOUNDARY_FILE_FLAG="--boundary-file="
   character(len=*), parameter   :: PARTITION_FILE_FLAG="--partition-file="

   contains

      subroutine getTaskType(taskType, meshName, no_of_solutions, solutionNames, solutionTypes, fixedOrder, Nout, basis,mode, useCommandArgs,oldStats)
         use Storage               , only: isOldStats
         implicit none
         integer,                                 intent(out) :: taskType
         character(len=*),                        intent(out) :: meshName
         integer,                                 intent(out) :: no_of_solutions
         character(len=LINE_LENGTH), allocatable, intent(out) :: solutionNames(:)
         integer, allocatable,                    intent(out) :: solutionTypes(:)
         logical,                                 intent(out) :: fixedOrder
         integer,                                 intent(out) :: Nout(3)
         integer,                                 intent(out) :: basis
         integer,                                 intent(out) :: mode
         logical,                                 intent(out) :: useCommandArgs
         logical,                                 intent(out) :: oldStats
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: no_of_arguments
         integer     :: i, sol, pos, pos2

!
!        Get number of command arguments         
!        -------------------------------
         no_of_arguments = command_argument_count()
!
!        If the control file is used, there is only 1 argument
!        -----------------------------------------------------
         if (no_of_arguments .eq. 1) then
             !todo: check that it is actually a control file
             useCommandArgs = .false.
             call getTaskTypeControl(taskType, meshName, no_of_solutions, solutionNames, solutionTypes, fixedOrder, Nout, basis, mode, oldStats)
         else
             useCommandArgs = .true.
             call getTaskTypeCommand(taskType, meshName, no_of_solutions, solutionNames, solutionTypes, fixedOrder, Nout, basis, mode, oldStats)
         end if 

         isOldStats = oldStats

      End Subroutine getTaskType

      subroutine getTaskTypeControl(taskType, meshName, no_of_solutions, solutionNames, solutionTypes, fixedOrder, Nout, basis,mode, oldStats)
         use FTValueDictionaryClass, only: FTValueDictionary
         use FileReaders           , only: ReadControlFile 
         use FileReadingUtilities, only: getCharArrayFromString
         use SolutionFile
         use Storage
         use OutputVariables       , only: outScale, hasVariablesFlag, askedVariables, Lreference
         use Utilities, only: toLower
         implicit none
         integer,                                 intent(out) :: taskType
         character(len=*),                        intent(out) :: meshName
         integer,                                 intent(out) :: no_of_solutions
         character(len=LINE_LENGTH), allocatable, intent(out) :: solutionNames(:)
         integer, allocatable,                    intent(out) :: solutionTypes(:)
         logical,                                 intent(out) :: fixedOrder
         integer,                                 intent(out) :: Nout(3)
         integer,                                 intent(out) :: basis
         integer,                                 intent(out) :: mode
         logical,                                 intent(out) :: oldStats

!
!        ---------------
!        Local variables
!        ---------------
!
         type( FTValueDictionary)                               :: controlVariables
         logical                                                :: meshFilePresent, basisPresent, modePresent, solutionPresent, patternPresent
         character(len=LINE_LENGTH)                             :: basisName, modeName, auxiliarName
         character(len=LINE_LENGTH)                             :: solutionsPattern, fullExpression
         real(kind=RP)                                           :: r
         integer                                                :: pos, pos2
         character(len=LINE_LENGTH)                             :: additionalVariablesStr, addVar
         character(len=LINE_LENGTH), dimension(:), allocatable  :: additionalVariablesArr
         integer                                                :: i, fID, reason
!
         call controlVariables % initWithSize(16)
         call ReadControlFile( controlVariables )
!
         Nout = 0
         no_of_solutions = 0
!
!        Check if the mesh file is present
!        ---------------------------------
         meshFilePresent = controlVariables % containsKey("hmesh file")
         if ( .not. meshFilePresent ) then
            write(STD_OUT,'(A)') "Mesh file was not specified"
            errorMessage(STD_OUT)
            taskType = UNKNOWN_JOB
            return
         end if
         meshName = controlVariables % stringValueForKey("hmesh file", LINE_LENGTH)

         solutionPresent = controlVariables % containsKey("hsol file")
         patternPresent = controlVariables % containsKey("hsol files pattern")
         if (solutionPresent) then
             no_of_solutions = 1
             allocate(solutionNames(no_of_solutions), solutionTypes(no_of_solutions))
             solutionNames(1) = controlVariables % stringValueForKey("hsol file", LINE_LENGTH)
             solutionTypes(1) = getSolutionFileType(solutionNames(1))

         ! use pattern match to get an array of soultions
         elseif (patternPresent) then
            solutionsPattern = controlVariables % stringValueForKey("hsol files pattern", LINE_LENGTH)
            ! get files in temporal txt
            write(fullExpression,'(A,A,A)') "ls ", trim(solutionsPattern), " > horses_temporal_file.txt"
            call system(trim(fullExpression))

            open ( newunit = fID , file = "horses_temporal_file.txt", status = "old" , action = "read" ) 
            i = 0
            do
                read(fID,fmt='(a)',iostat=reason) r
                if (reason/=0) exit
                i = i+1
            end do
            no_of_solutions = i

            allocate(solutionNames(no_of_solutions), solutionTypes(no_of_solutions))
            rewind(fID)
            do i = 1, no_of_solutions
                read(fID, '(A)') solutionNames(i)
                solutionTypes(i) = getSolutionFileType(solutionNames(i))
            end do 
         end if
!
!        Select the job type: Mesh if no solutions files are present, solution otherwise.
!        -------------------             
         if ( no_of_solutions .ne. 0 ) then
            taskType = SOLUTION_2_PLT
         else
            taskType = MESH_2_PLT
        end if
!
!        Get the order
!        -------------
         fixedOrder = controlVariables % containsKey("output order")
         if (fixedOrder) then
             auxiliarName = controlVariables % stringValueForKey("output order", LINE_LENGTH)
             if ( index(trim(auxiliarName(pos+1:len_trim(auxiliarName))),",") .eq. 0 ) then
!                Constant polynomial order: output order=N
!                -------------------------
                 Nout(1) = controlVariables % integerValueForKey("output order")
                 Nout(2:3) = Nout(1)
             else
!                Anisotropic polynomial order: output order=Nx,Ny,Nz
!                ----------------------------
                 pos2 = index(trim(auxiliarName),",")                  ! Read Nx
                 read(auxiliarName(1:pos2),*) Nout(1)
                 pos = index(trim(auxiliarName),",",BACK=.true.)       !Read Ny, and Nz
                 read(auxiliarName(pos2+1:pos),*) Nout(2)
                 read(auxiliarName(pos+1:len_trim(auxiliarName)),*) Nout(3)
             end if
         end if
!
!        Get the basis
!        -------------
         basisPresent = controlVariables % containsKey("output basis")
         if (basisPresent) basisName = controlVariables % stringValueForKey("output basis", LINE_LENGTH)
         call getBasis(basisPresent, basisName, fixedOrder, basis)
!
!        Get the mode
!        ------------
         mode = MODE_MULTIZONE ! Default
         modePresent = controlVariables % containsKey("output mode")
         if (modePresent) then
             modeName = controlVariables % stringValueForKey("output mode", LINE_LENGTH)
             select case(modeName)
             case ("multizone")
                 mode = MODE_MULTIZONE
             case ("FE")
                 mode = MODE_FINITEELM
             end select
         end if
!
!        Get the boundary and partition file
!        ------------------------------------
         hasBoundaries = controlVariables % containsKey("boundary file")
         if (hasBoundaries) boundaryFileName = controlVariables % stringValueForKey("boundary file", LINE_LENGTH)
         hasMPIranks = controlVariables % containsKey("partition file")
         if (hasMPIranks) partitionFileName = controlVariables % stringValueForKey("partition file", LINE_LENGTH)
         hasVariablesFlag = controlVariables % containsKey("output variables")
         if (hasVariablesFlag) askedVariables = controlVariables % stringValueForKey("output variables", LINE_LENGTH)
         outScale = .not. controlVariables % logicalValueForKey("dimensionless")
!
!        Get misc values
!        ------------------------------------
         oldStats = controlVariables % logicalValueForKey("legacy stats")
         Lreference = controlVariables % getValueOrDefault("reference length (m)", 1.0_RP)
         hasExtraGradients = controlVariables % logicalValueForKey("has gradients")
         if (controlVariables % containsKey("flow equations")) then
            flowEq = controlVariables%stringValueForKey("flow equations", LINE_LENGTH)
            call toLower(flowEq)
         else
             flowEq = "ns"
         end if

         if (controlVariables % containsKey("additional variables")) then
             additionalVariablesStr = controlVariables % stringValueForKey("additional variables", LINE_LENGTH)
             call getCharArrayFromString(trim(additionalVariablesStr), LINE_LENGTH, additionalVariablesArr)

             do i = 1, size(additionalVariablesArr)
                addVar = additionalVariablesArr(i)
                call toLower(addVar)
                select case (trim(addVar))
                    case ("u_tau")
                        hasUt_NS = .true.
                    case ("turb")
                        hasMu_NS = .true.
                        hasWallY     = .true.
                    case ("les")
                        hasMu_NS = .true.
                        hasMu_sgs     = .true.
                    case default
                        write(STD_OUT,'(A,A,A)') "The variable asked, ", trim(addVar), " is not implemented"
                end select
             end do
         end if

         isNewVariable = controlVariables % logicalValueForKey("is new variable")

      End Subroutine getTaskTypeControl

      subroutine getTaskTypeCommand(taskType, meshName, no_of_solutions, solutionNames, solutionTypes, fixedOrder, Nout, basis,mode,oldStats)
         use SolutionFile
         use Storage               , only: hasMPIranks, hasBoundaries, partitionFileName, boundaryFileName, flowEq
         use OutputVariables       , only: outScale, hasVariablesFlag, askedVariables, Lreference
         implicit none
         integer,                                 intent(out) :: taskType
         character(len=*),                        intent(out) :: meshName
         integer,                                 intent(out) :: no_of_solutions
         character(len=LINE_LENGTH), allocatable, intent(out) :: solutionNames(:)
         integer, allocatable,                    intent(out) :: solutionTypes(:)
         logical,                                 intent(out) :: fixedOrder
         integer,                                 intent(out) :: Nout(3)
         integer,                                 intent(out) :: basis
         integer,                                 intent(out) :: mode
         logical,                                 intent(out) :: oldStats
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                                              :: no_of_arguments
         integer                                              :: i, sol, pos, pos2
         integer                                              :: io, fid
         logical                                              :: meshFilePresent
         logical                                              :: basisPresent
         character(len=LINE_LENGTH)                           :: auxiliarName, basisName, modeName
!
!        Get number of command arguments         
!        -------------------------------
         no_of_arguments = command_argument_count()
!
!        Exit if no input arguments are specified
!        ----------------------------------------
         if ( no_of_arguments .eq. 0 ) then
            write(STD_OUT,'(A)') "No mesh file and/or solution file specified"
            taskType = UNKNOWN_JOB
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
            taskType = UNKNOWN_JOB
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

               if ( getSolutionFileType(auxiliarName) .ne. MESH_FILE .and. getSolutionFileType(auxiliarName) .ne. ZONE_MESH_FILE ) then
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
            taskType = SOLUTION_2_PLT
   
         else
            taskType = MESH_2_PLT

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
         call getBasis(basisPresent, basisName, fixedOrder, basis)
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
!        Get the boundary and partition file
!        Also, check if the dimensionless version is requested
!        ------------------------------------
         hasBoundaries    = .false.
         hasMPIranks      = .false.
         outScale         = .true.
         hasVariablesFlag = .false.
         do i = 1, no_of_arguments
            call get_command_argument(i, auxiliarName)
            pos = index(trim(auxiliarName),BOUNDARY_FILE_FLAG)
            if ( pos .ne. 0 ) then
               hasBoundaries = .TRUE.
               write(boundaryFileName,'(A)') auxiliarName(pos+len_trim(BOUNDARY_FILE_FLAG):len_trim(auxiliarName))
               cycle
            end if
            pos = index(trim(auxiliarName),PARTITION_FILE_FLAG)
            if ( pos .ne. 0 ) then
               hasMPIranks = .TRUE.
               read(auxiliarName(pos+len_trim(PARTITION_FILE_FLAG):len_trim(auxiliarName)),*) partitionFileName
               cycle
            end if
            pos = index(trim(auxiliarName),"--dimensionless")
            if ( pos .ne. 0 ) then
               outScale = .false.
               cycle
            end if
            pos = index(trim(auxiliarName),OUTPUT_VARIABLES_FLAG)
            if ( pos .ne. 0 ) then
               hasVariablesFlag = .true.
               write(askedVariables,'(A)') auxiliarName(pos+len_trim(OUTPUT_VARIABLES_FLAG):len_trim(auxiliarName))
               cycle
            end if
         end do
         oldStats = .false.
         Lreference = 1.0_RP

      end subroutine getTaskTypeCommand

      Subroutine getBasis(basisPresent, basisName, fixedOrder, basis)
         logical,                                 intent(in) :: basisPresent
         character(len=LINE_LENGTH),              intent(in) :: basisName
         logical,                                 intent(in) :: fixedOrder
         integer,                                 intent(out) :: basis
!
!        *******************************************************
!        If the basis is not present, the choice depends
!        on whether a new polynomial order is given.
!           * default polynomial: default basis is Gauss.
!           * given polynomial: default basis is Homogeneous.
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
!
      End Subroutine getBasis
!
end module getTask
