# Modify these variables for a given installation
#
F90        = /usr/local/bin/gfortran
FTOLibPath = /Users/kopriva/Documents/Research/FortranCode/LibrarySource/FTObjectLibrary/Source
NSLitePath = ../../..
# end
#
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
Assert.o : $(FTOLibPath)/FTTesting/Assert.f90 \
Comparisons.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(FTOLibPath)/FTTesting/Assert.f90

BoundaryConditions.o : $(NSLitePath)/Source/Physics/BoundaryConditions.f90 \
Physics.o \
Physics.o \
MeshTypes.o \
SharedBCModule.o \
SMConstants.o \
HexElementClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(NSLitePath)/Source/Physics/BoundaryConditions.f90

Comparisons.o : $(FTOLibPath)/FTTesting/Comparisons.f90
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(FTOLibPath)/FTTesting/Comparisons.f90

DGSEMClass.o : $(NSLitePath)/Source/Approximation/DGSEMClass.f90 \
DGSpaceAndTimeDerivativeMethods.o \
Physics.o \
NodalStorageClass.o \
SharedBCModule.o \
BoundaryConditions.o \
HexMesh.o \
Physics.o \
SMConstants.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(NSLitePath)/Source/Approximation/DGSEMClass.f90

DGSpaceAndTimeDerivativeMethods.o : $(NSLitePath)/Source/Approximation/DGSpaceAndTimeDerivativeMethods.f90 \
HexElementClass.o \
Physics.o \
NodalStorageClass.o \
MappedGeometry.o \
InterpolationAndDerivatives.o \
Physics.o \
SMConstants.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(NSLitePath)/Source/Approximation/DGSpaceAndTimeDerivativeMethods.f90

FaceClass.o : $(NSLitePath)/Source/Mesh/FaceClass.f90 \
SMConstants.o \
MeshTypes.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(NSLitePath)/Source/Mesh/FaceClass.f90

FacePatchClass.o : $(NSLitePath)/Source/Geometry/FacePatchClass.f90 \
SMConstants.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(NSLitePath)/Source/Geometry/FacePatchClass.f90

FileReading.o : $(NSLitePath)/Source/Foundation/FileReading.f90 \
SMConstants.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(NSLitePath)/Source/Foundation/FileReading.f90

FTDictionaryClass.o : $(FTOLibPath)/FTObjects/FTDictionaryClass.f90 \
FTObjectArrayClass.o \
FTLinkedListClass.o \
FTObjectClass.o \
FTLinkedListClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(FTOLibPath)/FTObjects/FTDictionaryClass.f90

FTLinkedListClass.o : $(FTOLibPath)/FTObjects/FTLinkedListClass.f90 \
FTObjectArrayClass.o \
FTObjectClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(FTOLibPath)/FTObjects/FTLinkedListClass.f90

FTMultiIndexTable.o : $(FTOLibPath)/FTObjects/FTMultiIndexTable.f90 \
FTObjectClass.o \
FTLinkedListClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(FTOLibPath)/FTObjects/FTMultiIndexTable.f90

FTObjectArrayClass.o : $(FTOLibPath)/FTObjects/FTObjectArrayClass.f90 \
FTObjectClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(FTOLibPath)/FTObjects/FTObjectArrayClass.f90

FTObjectClass.o : $(FTOLibPath)/FTObjects/FTObjectClass.f90
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(FTOLibPath)/FTObjects/FTObjectClass.f90

FTTimerClass.o : $(NSLitePath)/Source/Foundation/FTTimerClass.f90
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(NSLitePath)/Source/Foundation/FTTimerClass.f90

FTValueClass.o : $(FTOLibPath)/FTObjects/FTValueClass.f90 \
FTObjectClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(FTOLibPath)/FTObjects/FTValueClass.f90

FTValueDictionaryClass.o : $(FTOLibPath)/FTObjects/FTValueDictionaryClass.f90 \
FTValueClass.o \
FTDictionaryClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(FTOLibPath)/FTObjects/FTValueDictionaryClass.f90

Hash.o : $(FTOLibPath)/FTObjects/Hash.f90
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(FTOLibPath)/FTObjects/Hash.f90

HexElementClass.o : $(NSLitePath)/Source/Mesh/HexElementClass.f90 \
TransfiniteMaps3D.o \
HexElementConnectivityDefinitions.o \
NodalStorageClass.o \
MeshTypes.o \
LegendreAlgorithms.o \
InterpolationAndDerivatives.o \
MappedGeometry.o \
SMConstants.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(NSLitePath)/Source/Mesh/HexElementClass.f90

HexElementConnectivityDefinitions.o : $(NSLitePath)/Source/Mesh/HexElementConnectivityDefinitions.f90 \
SMConstants.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(NSLitePath)/Source/Mesh/HexElementConnectivityDefinitions.f90

HexMesh.o : $(NSLitePath)/Source/Mesh/HexMesh.f90 \
FTValueClass.o \
FTMultiIndexTable.o \
TransfiniteMaps3D.o \
FaceClass.o \
MeshTypes.o \
NodeClass.o \
HexElementClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(NSLitePath)/Source/Mesh/HexMesh.f90

InterpolationAndDerivatives.o : $(NSLitePath)/Source/Spectral/InterpolationAndDerivatives.f90 \
SMConstants.o \
LegendreAlgorithms.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(NSLitePath)/Source/Spectral/InterpolationAndDerivatives.f90

LegendreAlgorithms.o : $(NSLitePath)/Source/Spectral/LegendreAlgorithms.f90 \
SMConstants.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(NSLitePath)/Source/Spectral/LegendreAlgorithms.f90

MappedGeometry.o : $(NSLitePath)/Source/Geometry/MappedGeometry.f90 \
TransfiniteMaps3D.o \
SMConstants.o \
NodalStorageClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(NSLitePath)/Source/Geometry/MappedGeometry.f90

MeshTypes.o : $(NSLitePath)/Source/Mesh/MeshTypes.f90
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(NSLitePath)/Source/Mesh/MeshTypes.f90

NodalStorageClass.o : $(NSLitePath)/Source/Approximation/NodalStorageClass.f90 \
LegendreAlgorithms.o \
SMConstants.o \
InterpolationAndDerivatives.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(NSLitePath)/Source/Approximation/NodalStorageClass.f90

NodeClass.o : $(NSLitePath)/Source/Mesh/NodeClass.f90 \
SMConstants.o \
MeshTypes.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(NSLitePath)/Source/Mesh/NodeClass.f90

NSLite3DMain.o : $(NSLitePath)/Source/NSLite3DMain.f90 \
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
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(NSLitePath)/Source/NSLite3DMain.f90

Physics.o : $(NSLitePath)/Source/Physics/Physics.f90 \
SMConstants.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(NSLitePath)/Source/Physics/Physics.f90

Plotter.o : $(NSLitePath)/Source/Foundation/Plotter.f90 \
NodalStorageClass.o \
HexElementClass.o \
PlotterDataSource.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(NSLitePath)/Source/Foundation/Plotter.f90

PlotterDataSource.o : $(NSLitePath)/Source/Plotting/PlotterDataSource.f90 \
SMConstants.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(NSLitePath)/Source/Plotting/PlotterDataSource.f90

ProblemFile.o : $(NSLitePath)/Tests/Euler/UniformFlow/ProblemFile.f90 \
BoundaryConditions.o \
SMConstants.o \
Assert.o \
DGSEMClass.o \
Physics.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(NSLitePath)/Tests/Euler/UniformFlow/ProblemFile.f90

ReadInputFile.o : $(NSLitePath)/Source/IO/ReadInputFile.f90 \
SharedBCModule.o \
SMConstants.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(NSLitePath)/Source/IO/ReadInputFile.f90

SharedBCModule.o : $(NSLitePath)/Source/Physics/SharedBCModule.f90 \
FTValueDictionaryClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(NSLitePath)/Source/Physics/SharedBCModule.f90

SMConstants.o : $(NSLitePath)/Source/Foundation/SMConstants.f90
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(NSLitePath)/Source/Foundation/SMConstants.f90

TestSuiteManagerClass.o : $(FTOLibPath)/FTTesting/TestSuiteManagerClass.f90 \
Assert.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(FTOLibPath)/FTTesting/TestSuiteManagerClass.f90

TimeIntegrator.o : $(NSLitePath)/Source/Approximation/TimeIntegrator.f90 \
Physics.o \
Plotter.o \
ProblemFile.o \
SMConstants.o \
DGSEMClass.o \
InterpolationAndDerivatives.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(NSLitePath)/Source/Approximation/TimeIntegrator.f90

TransfiniteMaps3D.o : $(NSLitePath)/Source/Geometry/TransfiniteMaps3D.f90 \
SMConstants.o \
FacePatchClass.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(NSLitePath)/Source/Geometry/TransfiniteMaps3D.f90

Utilities.o : $(NSLitePath)/Source/Foundation/Utilities.f90 \
SMConstants.o
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ $(NSLitePath)/Source/Foundation/Utilities.f90

