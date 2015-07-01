# Modify these variables for a given installation
#
F90 = /usr/local/bin/gfortran
FTOLibPath = /Users/kopriva/Documents/Research/FortranCode/LibrarySource/FTObjectLibrary/Source
# end
FFLAGS = -O  -fopenmp
##########################
# Object Files for build #
##########################

OBJS = \
Assert.o \
BoundaryConditions.o \
Comparisons.o \
DGSEMClass.o \
DGSpaceAndTimeDerivativeMethods.o \
FaceClass.o \
FacePatchClass.o \
FileReading.o \
FTDictionaryClass.o \
FTLinkedListClass.o \
FTMultiIndexTable.o \
FTObjectArrayClass.o \
FTObjectClass.o \
FTTimerClass.o \
FTValueClass.o \
FTValueDictionaryClass.o \
Hash.o \
HexElementClass.o \
HexElementConnectivityDefinitions.o \
HexMesh.o \
InterpolationAndDerivatives.o \
LegendreAlgorithms.o \
MappedGeometry.o \
MeshTypes.o \
NodalStorageClass.o \
NodeClass.o \
NSLite3DMain.o \
Physics.o \
Plotter.o \
PlotterDataSource.o \
ProblemFile.o \
ReadInputFile.o \
SharedBCModule.o \
SMConstants.o \
TestSuiteManagerClass.o \
TimeIntegrator.o \
TransfiniteMaps3D.o \
Utilities.o \

NSLite3D : $(OBJS)
	 ${F90} -fopenmp -o $@ $(OBJS)

#######################################
# Object dependencies and compilation #
#######################################
Assert.o : /Users/kopriva/Documents/Research/FortranCode/LibrarySource/FTObjectLibrary/Source/FTTesting/Assert.f90 \
Comparisons.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/LibrarySource/FTObjectLibrary/Source/FTTesting/Assert.f90

BoundaryConditions.o : /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Tests/Euler/diffuser/BoundaryConditions.f90 \
Physics.o \
MeshTypes.o \
SMConstants.o \
HexElementClass.o \
SharedBCModule.o \
Physics.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Tests/Euler/diffuser/BoundaryConditions.f90

Comparisons.o : /Users/kopriva/Documents/Research/FortranCode/LibrarySource/FTObjectLibrary/Source/FTTesting/Comparisons.f90
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/LibrarySource/FTObjectLibrary/Source/FTTesting/Comparisons.f90

DGSEMClass.o : /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Approximation/DGSEMClass.f90 \
DGSpaceAndTimeDerivativeMethods.o \
Physics.o \
NodalStorageClass.o \
SharedBCModule.o \
BoundaryConditions.o \
HexMesh.o \
Physics.o \
SMConstants.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Approximation/DGSEMClass.f90

DGSpaceAndTimeDerivativeMethods.o : /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Approximation/DGSpaceAndTimeDerivativeMethods.f90 \
HexElementClass.o \
Physics.o \
NodalStorageClass.o \
MappedGeometry.o \
InterpolationAndDerivatives.o \
Physics.o \
SMConstants.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Approximation/DGSpaceAndTimeDerivativeMethods.f90

FaceClass.o : /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Mesh/FaceClass.f90 \
SMConstants.o \
MeshTypes.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Mesh/FaceClass.f90

FacePatchClass.o : /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Geometry/FacePatchClass.f90 \
SMConstants.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Geometry/FacePatchClass.f90

FileReading.o : /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Foundation/FileReading.f90 \
SMConstants.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Foundation/FileReading.f90

FTDictionaryClass.o : /Users/kopriva/Documents/Research/FortranCode/LibrarySource/FTObjectLibrary/Source/FTObjects/FTDictionaryClass.f90 \
FTObjectArrayClass.o \
FTLinkedListClass.o \
FTObjectClass.o \
FTLinkedListClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/LibrarySource/FTObjectLibrary/Source/FTObjects/FTDictionaryClass.f90

FTLinkedListClass.o : /Users/kopriva/Documents/Research/FortranCode/LibrarySource/FTObjectLibrary/Source/FTObjects/FTLinkedListClass.f90 \
FTObjectArrayClass.o \
FTObjectClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/LibrarySource/FTObjectLibrary/Source/FTObjects/FTLinkedListClass.f90

FTMultiIndexTable.o : /Users/kopriva/Documents/Research/FortranCode/LibrarySource/FTObjectLibrary/Source/FTObjects/FTMultiIndexTable.f90 \
FTObjectClass.o \
FTLinkedListClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/LibrarySource/FTObjectLibrary/Source/FTObjects/FTMultiIndexTable.f90

FTObjectArrayClass.o : /Users/kopriva/Documents/Research/FortranCode/LibrarySource/FTObjectLibrary/Source/FTObjects/FTObjectArrayClass.f90 \
FTObjectClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/LibrarySource/FTObjectLibrary/Source/FTObjects/FTObjectArrayClass.f90

FTObjectClass.o : /Users/kopriva/Documents/Research/FortranCode/LibrarySource/FTObjectLibrary/Source/FTObjects/FTObjectClass.f90
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/LibrarySource/FTObjectLibrary/Source/FTObjects/FTObjectClass.f90

FTTimerClass.o : /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Foundation/FTTimerClass.f90
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Foundation/FTTimerClass.f90

FTValueClass.o : /Users/kopriva/Documents/Research/FortranCode/LibrarySource/FTObjectLibrary/Source/FTObjects/FTValueClass.f90 \
FTObjectClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/LibrarySource/FTObjectLibrary/Source/FTObjects/FTValueClass.f90

FTValueDictionaryClass.o : /Users/kopriva/Documents/Research/FortranCode/LibrarySource/FTObjectLibrary/Source/FTObjects/FTValueDictionaryClass.f90 \
FTValueClass.o \
FTDictionaryClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/LibrarySource/FTObjectLibrary/Source/FTObjects/FTValueDictionaryClass.f90

Hash.o : /Users/kopriva/Documents/Research/FortranCode/LibrarySource/FTObjectLibrary/Source/FTObjects/Hash.f90
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/LibrarySource/FTObjectLibrary/Source/FTObjects/Hash.f90

HexElementClass.o : /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Mesh/HexElementClass.f90 \
TransfiniteMaps3D.o \
HexElementConnectivityDefinitions.o \
NodalStorageClass.o \
MeshTypes.o \
LegendreAlgorithms.o \
InterpolationAndDerivatives.o \
MappedGeometry.o \
SMConstants.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Mesh/HexElementClass.f90

HexElementConnectivityDefinitions.o : /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Mesh/HexElementConnectivityDefinitions.f90 \
SMConstants.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Mesh/HexElementConnectivityDefinitions.f90

HexMesh.o : /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Mesh/HexMesh.f90 \
FTValueClass.o \
FTMultiIndexTable.o \
TransfiniteMaps3D.o \
FaceClass.o \
MeshTypes.o \
NodeClass.o \
HexElementClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Mesh/HexMesh.f90

InterpolationAndDerivatives.o : /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Spectral/InterpolationAndDerivatives.f90 \
SMConstants.o \
LegendreAlgorithms.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Spectral/InterpolationAndDerivatives.f90

LegendreAlgorithms.o : /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Spectral/LegendreAlgorithms.f90 \
SMConstants.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Spectral/LegendreAlgorithms.f90

MappedGeometry.o : /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Geometry/MappedGeometry.f90 \
TransfiniteMaps3D.o \
SMConstants.o \
NodalStorageClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Geometry/MappedGeometry.f90

MeshTypes.o : /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Mesh/MeshTypes.f90
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Mesh/MeshTypes.f90

NodalStorageClass.o : /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Approximation/NodalStorageClass.f90 \
LegendreAlgorithms.o \
SMConstants.o \
InterpolationAndDerivatives.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Approximation/NodalStorageClass.f90

NodeClass.o : /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Mesh/NodeClass.f90 \
SMConstants.o \
MeshTypes.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Mesh/NodeClass.f90

NSLite3DMain.o : /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/NSLite3DMain.f90 \
Physics.o \
TimeIntegrator.o \
ProblemFile.o \
HexMesh.o \
Plotter.o \
FTTimerClass.o \
SharedBCModule.o \
BoundaryConditions.o \
ReadInputFile.o \
SMConstants.o \
DGSEMClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/NSLite3DMain.f90

Physics.o : /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Physics/Physics.f90 \
SMConstants.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Physics/Physics.f90

Plotter.o : /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Foundation/Plotter.f90 \
NodalStorageClass.o \
HexElementClass.o \
PlotterDataSource.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Foundation/Plotter.f90

PlotterDataSource.o : /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Tests/Euler/diffuser/PlotterDataSource.f90 \
SMConstants.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Tests/Euler/diffuser/PlotterDataSource.f90

ProblemFile.o : /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Tests/Euler/diffuser/ProblemFile.f90 \
Physics.o \
Physics.o \
BoundaryConditions.o \
MeshTypes.o \
Assert.o \
SMConstants.o \
DGSEMClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Tests/Euler/diffuser/ProblemFile.f90

ReadInputFile.o : /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/IO/ReadInputFile.f90 \
SharedBCModule.o \
SMConstants.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/IO/ReadInputFile.f90

SharedBCModule.o : /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Physics/SharedBCModule.f90 \
FTValueDictionaryClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Physics/SharedBCModule.f90

SMConstants.o : /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Foundation/SMConstants.f90
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Foundation/SMConstants.f90

TestSuiteManagerClass.o : /Users/kopriva/Documents/Research/FortranCode/LibrarySource/FTObjectLibrary/Source/FTTesting/TestSuiteManagerClass.f90 \
Assert.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/LibrarySource/FTObjectLibrary/Source/FTTesting/TestSuiteManagerClass.f90

TimeIntegrator.o : /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Approximation/TimeIntegrator.f90 \
Physics.o \
Plotter.o \
ProblemFile.o \
SMConstants.o \
DGSEMClass.o \
InterpolationAndDerivatives.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Approximation/TimeIntegrator.f90

TransfiniteMaps3D.o : /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Geometry/TransfiniteMaps3D.f90 \
SMConstants.o \
FacePatchClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Geometry/TransfiniteMaps3D.f90

Utilities.o : /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Foundation/Utilities.f90 \
SMConstants.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ /Users/kopriva/Documents/Research/FortranCode/NSLite3D/Source/Foundation/Utilities.f90

