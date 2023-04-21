#!/bin/bash

renumberMesh -overwrite > log.renumberMesh &
wait
checkMesh > log.checkMesh &
wait
decomposePar -latestTime > log.decomposePar &
wait
mpirun -np 8 rhoPimpleFoam -parallel > log.rhoPimpleFoam_1 &
wait
reconstructPar -latestTime > log.reconstructPar &
wait
foamToVTK -fields "(U p T rho)" -ascii -latestTime > log.foamToVTK_1 &
wait





