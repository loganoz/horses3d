mkdir -p include
mkdir -p build
mkdir -p ./Tests/Euler/BoxAroundCircle/RestartFiles
mkdir -p ./Tests/Euler/BoxAroundCircle/PlotFiles
mkdir -p ./Tests/Euler/diffuser/RestartFiles
mkdir -p ./Tests/Euler/diffuser/PlotFiles
mkdir -p ./Tests/Euler/PeriodicFlow/RestartFiles
mkdir -p ./Tests/Euler/PeriodicFlow/PlotFiles
mkdir -p ./Tests/Euler/UniformFlow/RestartFiles
mkdir -p ./Tests/Euler/UniformFlow/PlotFiles
mkdir -p ./Tests/NavierStokes/Cylinder/RestartFiles
mkdir -p ./Tests/NavierStokes/Cylinder/PlotFiles
mkdir -p ./Tests/NavierStokes/FlatPlate/RestartFiles
mkdir -p ./Tests/NavierStokes/FlatPlate/PlotFiles
cp --verbose ./Tests/make.inc ./Tests/Components/FacePatches/make.inc
cp --verbose ./Tests/make.inc ./Tests/Components/Gradients/make.inc
cp --verbose ./Tests/make.inc ./Tests/Components/HexMappings/make.inc
cp --verbose ./Tests/make.inc ./Tests/Components/HexMesh/make.inc
cp --verbose ./Tests/make.inc ./Tests/Components/MappedGeometry/make.inc
cp --verbose ./Tests/make.inc ./Tests/Components/NodalStorage/make.inc
cp --verbose ./Tests/make.inc ./Tests/Euler/BoxAroundCircle/make.inc
cp --verbose ./Tests/make.inc ./Tests/Euler/diffuser/make.inc
cp --verbose ./Tests/make.inc ./Tests/Euler/PeriodicFlow/make.inc
cp --verbose ./Tests/make.inc ./Tests/Euler/UniformFlow/make.inc
cp --verbose ./Tests/make.inc ./Tests/NavierStokes/Cylinder/make.inc
cp --verbose ./Tests/make.inc ./Tests/NavierStokes/FlatPlate/make.inc


