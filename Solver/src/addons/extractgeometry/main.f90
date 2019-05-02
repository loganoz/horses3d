!
!//////////////////////////////////////////////////////
!
!   @File:    main.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Wed Nov  1 19:56:54 2017
!   @Last revision date: Mon Feb 25 16:07:49 2019
!   @Last revision author: Andrés Rueda (am.rueda@upm.es)
!   @Last revision commit: 17d60e4e57235a57aa406023ebe4c26157bc211a
!
!////////////////////////////////////////////////////
!
program ExtractGeometry
   use SMConstants
   use Headers
   use NodalStorageClass
   use HexMeshClass
   use ConstructMeshAndSpectralBasis_MOD
   use FTValueDictionaryClass
   use getInputData_MOD
   use GeometryClass
   use MPI_Process_Info
   implicit none
   character(len=LINE_LENGTH)      :: meshFile
   character(len=LINE_LENGTH)      :: solutionFile
   type(HexMesh)                   :: mesh
   type(FTValueDictionary)         :: controlVariables
   class(Geometry_t), pointer      :: geometry
!
!  Init MPI
!  --------
   call MPI_Process % Init
!
!  Display header
!  --------------
   call Main_Header("Extract geometry utility",__DATE__,__TIME__)
!
!  Initialize control variables
!  ----------------------------
   call controlVariables % InitWithSize(16)
   call getInputData(controlVariables)
!
!  Construct mesh and spectral basis
!  ---------------------------------
   call ConstructMeshAndSpectralBasis( controlVariables % StringValueForKey(meshFileKey, LINE_LENGTH), &
                                       controlVariables % StringValueForKey(solutionFileKey, LINE_LENGTH), &
                                       controlVariables % StringValueForKey(ControlFileKey, LINE_LENGTH), &
                                       mesh, controlVariables)
!
!  Create geometry
!  ---------------
   geometry => Construct(mesh, NodalStorage, controlVariables)
!
!  Get new points coordinates and elements
!  ---------------------------------------
   call geometry % Compute(mesh, NodalStorage)
!
!  Close MPI
!  ---------
   call MPI_Process % Close


end program ExtractGeometry
