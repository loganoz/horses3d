!
!//////////////////////////////////////////////////////
!
!   @File:    main.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Wed Nov  1 19:56:54 2017
!   @Last revision date: Sun Nov  5 19:15:13 2017
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 431d0b8be7da5b914d9e787fc6ac9e78aceca4ef
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
   implicit none
   character(len=LINE_LENGTH)      :: meshFile
   character(len=LINE_LENGTH)      :: solutionFile
   type(NodalStorage), allocatable :: spA(:)
   type(HexMesh)                   :: mesh
   type(FTValueDictionary)         :: controlVariables
   class(Geometry_t), pointer      :: geometry
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
                                       mesh, spA)
!
!  Create geometry
!  ---------------
   geometry => Construct(mesh, spA, controlVariables)
!
!  Get new points coordinates and elements
!  ---------------------------------------
   call geometry % Compute(mesh, spA)


end program ExtractGeometry
