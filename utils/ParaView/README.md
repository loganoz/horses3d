# Postprocessing of HORSES3D results in ParaView with MPI in HPC platforms

This directory contains a set of `pvbatch` scripts that can serve as a guideline on how to postprocess large HORSES3D result files in an HPC platform using ParaView and MPI. Some of the most common postprocessing tasks are addressed:

- Distributing large solution files into a partitioned format and writing them to disk for efficient parallel handling using MPI.
- Computing additional fields using gradient and divergence operators among others.
- Generating 2D slices from 3D results and writing them to disk in *.csv format for easy plotting.
- Rendering 3D views of CFD results using thresholds and opacity schemes.

Additionally, `python` scripts to perform other pre- and post-proecssing tasks are provided

>[!IMPORTANT]
> These `pvbatch` scripts have been tested on ParaView versions 5.11.0-MPI, 5.12.0-MPI and 5.13.2-MPI. Using different versions might lead to unexpected issues.

>[!TIP]
> The batch python script tracing tool included in ParaView's GUI provides an easy way to quickly adapt these batch scripts to perform any given task under different ParaView versions.

## Conversion to partitioned format (*.pvtu)

The `pvbatch` script [pv_distribute](./pv_distribute.py) redistributes the input solution file (in *.hdf, *.vtk, *.vtkhdf, or *.tec format) data to `np` MPI processes and writes the resulting partitioned data to disk in *.pvtu format for its subsequent efficient use in ParaView with MPI. Performing this preprocessing of the data is particularly adventageous when multiple ParaView postprocessing tasks are going to be executed since redistribution operations can be time consuming. Below is an example of how to use this script:

To print script `help` with a brief description of its usage and arguments definition, use

```
<path_to_pvbatch_mpi> pv_distribute.py -h
```

## Extracting 2D slices

The `pvbatch` script [pv_slice_mpi.py](./pv_slice_mpi.py) extracts a set of 2D slices from an input 3D solution file according to the specified positional arguments (plane origin coordinates, normal vector, output fields...), having previously computed some additional fields if requested by the user. In this particular case, the script will compute the Lamb vector and its divergence if either `Lamb` or `divLamb` are specified as output fields. The auxiliary vorticity field `omega` will also be computed unless already present in the input 3D solution file (and specified via script flag). The resulting 2D slices are writen to files in ASCII *.csv format for easy postprocessing.

To print script `help` with a brief description of usage and arguments definition, use 

```
<path_to_pvbatch_mpi> pv_slice_mpi.py -h
```

### Plotting extracted 2D slices

The python script [plot_pv_slices.py](./postprocessing/plot_pv_slices.py) uses matplotlib to create plots of 2D slices written as *.csv files in the same format as the outputs of [pv_slice_mpi.py](./pv_slice_mpi.py). It can be used simply as

```
python plot_pv_slices.py <infile_path>
```

where `<infile_path>` is the path to the input *.json driver file, which contains all the information necessary to generate the desired plots. The file [plot_slices_inputs_example.json](./postprocessing/plot_slices_inputs_example.json) can serve as a format example.

Additionally, the following table gathers the definition of every input parameter

<div style="margin-left: auto;
            margin-right: auto;
            width: 90%">

| Parameter name      | Type      |  Definition     |
| ------------- | ------------- | ------------- |
| dataPathList | list | List of paths to CSV files containing the data to be plotted |
| nx | int | Resolution of the x-axis data in the output plots |
| ny | int | Resolution of the y-axis data in the output plots |
| xlabel | str | Label of the x-axis |
| ylabel | str or list(str) | Label(s) of the y-axis (or axes) |
| xticks | list(float) | Values of x-axis ticks |
| yticks | list(float) or list(list(float)) | Values of y-axis (or axes) ticks  |
| xtickLabels | list(str) | Values of x-axis tick labels |
| ytickLabels | list(str) or list(list(str)) | Values of y-axis (or axes) ticks  |
| xCol | int or list(int) | Index or list of indexes of the column in the input CSV file(s) representing the x-axis values |
| yCol | int or list(int) | Index or list of indexes of the column in the input CSV file(s) representing the y-axis values |
| dataCol | int or list(int) | Index or list of indexes of the column in the input CSV file(s) representing the values of the field that wants to be plotted |
| dataLabel | str | Label of the field that wants to be plotted |
| figwidth | float | Total width of the figure (in inches) |
| outputPath | str | Path of the output figure, including file format extension. DPI for non-vector formats will be 300 by default |
| plotArray | list(int), optional | List defining the subplot array when including multiple plots in a single figure, [Nrows, Ncols]. Default value of [Nplots, 1] |
| yxRatio | float, optional | Aspect ratio to apply to y-axis data with respect to x-axis data when plotting. Default value of 1.0 |
| invertXaxis | bool, optional | Whether to invert the x-axis for representation (True) or not (False). Default value is False |
| colorMap | str, optional | String identifier of the matplotlib colorbar to use. Default value is "coolwarm" |
| colorMapMin | float, optional | Color map minimum value for field representation. If unspecified, it will automatically fit to the minimum field value in the figure |
| colorMapMax | float, optional | Color map maximum value for field representation. If unspecified, it will automatically fit to the maximum field value in the figure |
| colorBarTicks | list(float), optional | Colorbar tick values. If unspecified, they will automatically adjust to accomodate the field data range in the figure |
| colorBarTicks | list(float), optional | Colorbar tick values. If unspecified, they will automatically adjust to accomodate the field data range in the figure | 
| colorMapNorm | str, optional | Colobar norm for field reopresentation, either "linear" or "log". Default value is "linear" |
| fontsize | float, optional | Font size in pt. Default value is 9 |
| figLabelLoc | list(float), optional | Position of each subplot label ['(a)', '(b)', '(c)', ...] with respect to the bottom left corner of the plot area, [x, y] (expressed as fraction of the axis length). Labels are only included inthe figure if this argument is defined. Default is set to None |
| xLim | list(float) or list(list(float)), optional | X-axis (or axes) limits for field representation, [xmin, xmax]. Default values are the minimum and maximum x-axis (or axes) values in the input file(s) |
| yLim | list(float) or list(list(float)), optional | Y-axis (or axes) limits for field representation, [ymin, ymax]. Default values are the minimum and maximum y-axis (or axes) values in the input file(s) |
| xMask | list(str), optional | List of masks to apply to the input data fields, defined as a string function of "y". If x_point < x(y), field value for that point is set to 'NaN'. Numpy functions can be used under "np.*". Default value is None.
| yMask | list(str), optional | List of masks to apply to the input data fields, defined as a string function of "x". If y_point < y(x), field value for that point is set to 'NaN'. Numpy functions can be used under "np.*". Default value is None. 

</div>

>[!IMPORTANT]
> This python-based tool relies not only on the matplotlib package, but also on numpy and scipy. Additionally, the Qt5Agg backend to matplotlib must be available, as well as a latex compiler for text rendering (this last option can be easily modified by setting "text.usetex" to False in the script).

## Rendering 3D views

The `pvbatch` script [pv_render_mpi.py](pv_render_mpi.py) serves as an example on how to obtain 3D renders of HORSES3D solutions using ParaView and MPI. Particularly, it generates a 3D representation of the vorticity of the input solution file. If the vorticity field 'omega_abs' is already included in the input file, the script flag `--omega-in-file` can be used to avoid recomputing this quantity. If not specified, ParaView's gradient operator will be used to compute all the velocity component derivatives and ultimately the vorticity magnitude. Threshold values are specified in the script (not as arguments) but can be easily identified and modified. A custom opacity map is defined via a hyperbolic tangent function for coloring the resulting threshold view. Finally, some external STL files are loaded into the view for the representation of solid closed surfaces.

To print script help with a brief description of usage and arguments definition, use

```
<path_to_pvbatch_mpi> pv_render_mpi.py -h
```

## Usage example

Below there are two examples of how to use these `pvbatch` scripts in different HPC platforms where ParaView MPI is available. 

The [first example](#flexo) depicts how to convert an input VTK file (generated from a HORSES3D HDF file generated using the `horses2plt` tool and then converted to VTK PolyData using [vtkhdf_to_vtk.py](./preprocessing/vtkhdf_to_vtk.py)) as `ntasks` partitioned VTU files (PVTU). Then, the PVTU files are used to extract three 2D slices defined by three points (first positional argument) and three normal vectors (second positional argument) and the list of specified fields (third positional argument) will be output to the specified file paths (fourth and last positional argument). The `--omega-in-file` flag specified that the input PVTU file already contains the vorticity field, hence there is no need to recompute it. Finally, a 3D render of a vorticity field is generated and writen to the specified output PNG file.

The [second example](#marenostrum5) depicts how to extract the same three 2D slices directly from an HDF file generated using the `horses2plt` tool.

>[!NOTE]
>Make sure the output directories exist before executing the `pvbatch` scripts to avoid i/o errors.

### Flexo
```
#!/bin/bash
#SBATCH --job-name=pvbatch
#SBATCH --nodelist=n005
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=04:00:00
#SBATCH -e pvbatch-err-%j.log
#SBATCH -o pvbatch-out-%j.log
 
export OMP_SCHEDULE="guided"
export OMP_NUM_THREADS=1
 
module purge
module load paraview/5.11

export LD_LIBRARY_PATH=/opt/ohpc/pub/apps/miniconda/miniconda3/pkgs/xorg-libxcursor-1.2.3-hb9d3cd8_0/lib:$LD_LIBRARY_PATH
 
echo "##########################################################################"
echo "#"
echo "# Running with $SLURM_NTASKS tasks and $SLURM_CPUS_PER_TASK cpus/task"
echo "# On nodes $SLURM_JOB_NODELIST"
echo "#"
echo "##########################################################################"
 
EXEC=pvbatch

SCRIPT_DIST=scripts/pv_distribute.py
SCRIPT_SLICE=scripts/pv_slice_mpi.py
SCRIPT_RENDER=scripts/pv_render_mpi.py

VTK_FILE=files/horses3d_solution.vtk
PVTU_FILE=outputs/horses3d_solution.pvtu

OUT_RENDER=outputs/horses3d_solution_voricity.png

base_name=$(basename $VTK_FILE .vtk)

mpiexec -np $SLURM_NTASKS --hosts $SLURM_NODELIST $EXEC $SCRIPT_DIST $VTK_FILE $PVTU_FILE

wait

mpiexec -np $SLURM_NTASKS --hosts $SLURM_NODELIST $EXEC $SCRIPT_SLICE $PVTU_FILE [[0,0.001,0],[0,0,119.0],[-7.07,0,0]] [[0,1,0],[0,0,1],[1,0,0]] "['u','v','w','p','rho','Nav','Lamb','divLamb']" [\"paraview/${base_name}_xz_0D_slice.csv\",\"paraview/${base_name}_xy_zHub_slice.csv\",\"paraview/${base_name}_yz_rotPlane_slice.csv\"] --omega_in_file

wait

mpiexec -np $SLURM_NTASKS --hosts $SLURM_NODELIST $EXEC $SCRIPT_RENDER $PVTU_FILE $OUT_RENDER --omega_in_file
```

### Marenostrum5
```
#!/bin/bash
#SBATCH --job-name=pv_slice
#SBATCH --qos=gp_resa
#SBATCH --chdir=.
#SBATCH --ntasks=336
#SBATCH --cpus-per-task=4
#SBATCH --ntasks-per-node=28
#SBATCH --time=04:00:00
#SBATCH -e pv_slice-err-%j.log
#SBATCH -o pv_slice-out-%j.log

export OMP_SCHEDULE="guided"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Load modules for pvbatch
module purge
module load paraview/5.12.0-osmesa-MPI-Linux-Python3.10

EXEC=pvbatch
SCRIPT=scripts/pv_slice_mpi.py
HDF_FILE=files/horses3d_solution.hdf
base_name=$(basename $HDF_FILE .hdf)


# Run the slicing command
srun $EXEC $SCRIPT $HDF_FILE [[0,0.001,0],[0,0,119.0],[-7.07,0,0]] [[0,1,0],[0,0,1],[1,0,0]] "['u','v','w','p','rho','Nav','Lamb','divLamb']" [\"paraview/${base_name}_xz_0D_slice.csv\",\"paraview/${base_name}_xy_zHub_slice.csv\",\"paraview/${base_name}_yz_rotPlane_slice.csv\"] --omega_in_file

```


## Most common issues and recommendations

* Segmentation faults when performing different tasks (slicing, applying thresholds, ...) on redistributed data may arise depdending on the ParaView version and the input file format. The safest workflow to avoid these issues consists in writing the HORSES3D solution file in `vtkhdf` format when using `horses2plt` and then converting it to a VTK PolyData format file. This format typically requires more disk space, but its handling in ParaView when using MPI is more robust. The preprocessing python script [vtkhdf_to_vtk.py](./preprocessing/vtkhdf_to_vtk.py) can be used to do so. The resulting output VTK file can then be re-written in a partitioned format using [pv_distribute.py](pv_distribute.py) or directly loaded ParaView and internally redistributed before its processing.

* Depending on the ParaView-MPI version and the input file format, the element limit per process after redistributing data may be exceeded, leading to errors when attempting to perform any other tasks. Once more, the recommended workflow is to use VTK PolyData for the input files (see previous bullet point).

* Depending on the HPC platform you are using, running `pvbatch` with `mpiexec` or `srun` might be the right option. When encountering errors, try both options before considering other possible issues.

* Some ParaView installations and HPC setups might lead to OOM errors when executing multiple `pvbatch` scripts under a single `slurm` job, even when done sequentially, due to memory deallocation issues. If this happens, consider submitting each `pvbatch` execution as a single job.
