from paraview.simple import *
from paraview import servermanager

from mpi4py import MPI

import os
import gc
import ast
import logging
import argparse
import numpy as np

# Initialise MPI 
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Functions
def parse_list(arg):
    try:
        return ast.literal_eval(arg)
    except ValueError:
        raise argparse.ArgumentTypeError("Argument must be a list")
    
def parse_list_or_string(arg):
    try:
        return ast.literal_eval(arg)
    except (ValueError, SyntaxError):
        return arg
    
def slice_domain(domain, name, origin, normal, outfields, outfile, logger):
    
    logger.info(f"Creating slice {name} with origin {origin} and normal {normal}")
    
    # Create a new 'Slice'
    slice = Slice(registrationName=name, Input=domain)

    # Set slice origin and normal vector
    slice.SliceType.Origin = origin
    slice.SliceType.Normal = normal

    slice.UpdatePipeline()

    logger.info(f"Saving slice to file {os.path.abspath(outfile)} with fields {outfields}")
    # Save point data to*.csv output file
    SaveData(os.path.abspath(outfile), proxy=slice, ChooseArraysToWrite=1,
             PointDataArrays=outfields,
             Precision=6)
    
    logger.info(f"Deleting slice")
    Delete(slice)
    del slice
    
def get_input_list(input_list, length, idx):

    assert(isinstance(input_list, list))

    if all(isinstance(elem, list) for elem in input_list) and len(input_list) == length:
        output = input_list[idx]
    else:
        output = input_list

    return output

def compute_vorticity(reader, logger): 
  # Compute gradients of each velocity component
  logger.info(f"Computing gradient of velocity component 'u' using {size} processes...")
  gradient_u = Gradient(registrationName='Gradient_u', Input=reader)
  gradient_u.ScalarArray = ['POINTS', 'u']
  gradient_u.ResultArrayName = 'Gradient_u'
  # gradient_u.UpdatePipeline()

  logger.info(f"Computing gradient of velocity component 'v' using {size} processes...")
  gradient_v = Gradient(registrationName='Gradient_v', Input=gradient_u)
  gradient_v.ScalarArray = ['POINTS', 'v']
  gradient_v.ResultArrayName = 'Gradient_v'
  # gradient_v.UpdatePipeline()

  logger.info(f"Computing gradient of velocity component 'w' using {size} processes...")
  gradient_w = Gradient(registrationName='Gradient_w', Input=gradient_v)
  gradient_w.ScalarArray = ['POINTS', 'w']
  gradient_w.ResultArrayName = 'Gradient_w'
  # gradient_w.UpdatePipeline()

  # Compute vorticity components
  logger.info(f"Computing voriticty component 'omega_x' using {size} processes...")
  calc_omega_x = Calculator(registrationName='Calc_omega_x', Input=gradient_w)
  calc_omega_x.ResultArrayName = 'omega_x'
  calc_omega_x.Function = 'Gradient_w_Y-Gradient_v_Z'
  # calc_omega_x.UpdatePipeline()

  logger.info(f"Computing voriticty component 'omega_y' using {size} processes...")
  calc_omega_y = Calculator(registrationName='Calc_omega_y', Input=calc_omega_x)
  calc_omega_y.ResultArrayName = 'omega_y'
  calc_omega_y.Function = 'Gradient_u_Z-Gradient_w_X'
  # calc_omega_y.UpdatePipeline()

  logger.info(f"Computing voriticty component 'omega_z' using {size} processes...")
  calc_omega_z = Calculator(registrationName='Calc_omega_z', Input=calc_omega_y)
  calc_omega_z.ResultArrayName = 'omega_z'
  calc_omega_z.Function = 'Gradient_v_X-Gradient_u_Y'
  # calc_omega_z.UpdatePipeline()

  logger.info(f"Computing voriticty magnitude 'omega_abs' using {size} processes...")
  calc_omega_abs = Calculator(registrationName='Calc_omega_abs', Input=calc_omega_z)
  calc_omega_abs.ResultArrayName = 'omega_abs'
  calc_omega_abs.Function = 'sqrt(omega_x^2+omega_y^2+omega_z^2)'

  calc_omega_abs.UpdatePipeline()

  return calc_omega_abs

def compute_div_lamb(reader, logger, with_omega=False):

  # Compute vorticity if necessary
  if not with_omega:

    last_item = compute_vorticity(reader, logger)

  else: 
    last_item = reader

  # Compute Lamb vector from velocity and vorticity fields
  logger.info(f"Computing Lamb vector 'Lamb' using {size} processes...")
  calc_lamb = Calculator(registrationName='Calc_lamb', Input=last_item)
  calc_lamb.ResultArrayName = 'Lamb'
  calc_lamb.Function = '(v*omega_z-w*omega_y)*iHat + (w*omega_x - u*omega_z)*jHat + (u*omega_y - v*omega_x)*kHat'
  calc_lamb.UpdatePipeline()

  # Compute divergence of Lamb vector
  logger.info(f"Computing Lamb vector divergence 'divLamb' using {size} processes...")
  calc_divLamb = Gradient(registrationName='Gradient_divLamb', Input=calc_lamb)
  calc_divLamb.ComputeGradient = 0
  calc_divLamb.ComputeDivergence = 1
  calc_divLamb.DivergenceArrayName = 'divLamb'
  calc_divLamb.UpdatePipeline()

  return calc_divLamb

def main():
    
  parser = argparse.ArgumentParser(description="Slices a 3D domain dataset in *.hdf, *.vtk, *.vtkhdf or *.pvtu format and writes the resulting slice point data to an ASCII file in *.csv format.")
  
  # Define postional arguments
  parser.add_argument('infile_path', type=str, help='Path to the input *.hdf, *.vtk, *.vtkhdf or *.pvtu file.')
  parser.add_argument('outfile_path', type=str, help='Path to the output *.pdf containing the 3D rendering.')
  parser.add_argument('--omega_in_file', action='store_true', help='Use this flag if the omega field is already in the input file. Default is False, which means that the omega field will be computed from the velocity field.')

  # Parse arguments
  args = parser.parse_args()

  # Configure logging
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
  logger = logging.getLogger()

  _, ext = os.path.splitext(args.infile_path)

  if ext in [".hdf", ".vtkhdf"]:
    logger.info(f"Loading input {args.infile_path} into master node...")
    reader = VTKHDFReader(registrationName=os.path.basename(args.infile_path), FileName=[args.infile_path])
    #reader.PointArrayStatus = args.outfile_fields
    reader.UpdatePipeline()

  elif ext == ".vtk":
    logger.info(f"Loading input {args.infile_path} into master node...")
    reader = LegacyVTKReader(registrationName=os.path.basename(args.infile_path), FileNames=[args.infile_path])
    reader.UpdatePipeline()
    
  elif ext == ".tec":
    logger.info(f"Loading input {args.infile_path} into master node...")
    reader = TecplotReader(registrationName=os.path.basename(args.infile_path), FileNames=[args.infile_path])
    reader.UpdatePipeline()

  elif ext == ".pvtu":
    logger.info(f"Reading partitioned input {args.infile_path} into {size} processes...")
    reader = XMLPartitionedUnstructuredGridReader(registrationName=os.path.basename(args.infile_path), FileName=[args.infile_path])
    reader.TimeArray = 'None'
    redistribute = reader

  if ext in [".hdf", ".vtk", ".vtkhdf"]:
    # Redistribute data
    logger.info(f"Redistributing data into {size} processes...")
    redistribute = RedistributeDataSet(Input=reader)
    redistribute.UpdatePipeline()
    redistributed_data = servermanager.Fetch(redistribute)
    logger.info(f"Data correctly distributed and balanced. Process {rank}: Points = {redistributed_data.GetNumberOfPoints()}")

  logger.info(f"Rendering initial view using {size} processes...")
  renderView1 = GetActiveViewOrCreate('RenderView')
  renderView1.UseLight = 0
  display = Show(redistribute, renderView1, 'UnstructuredGridRepresentation')
  display.Representation = 'Surface'
  renderView1.Update()
  Hide(redistribute, renderView1)

  ######################### INSERT HERE CUSTOM RENDER CODE HERE  #########################

  # Compute vorticity if necessary
  if not args.omega_in_file:
    calc_vorticity = compute_vorticity(redistribute, logger)  
  else:
    calc_vorticity = redistribute

  # create a new 'Threshold'
  logger.info(f"Applying threshold to vorticity field 'omega_abs' using {size} processes...")
  threshold1 = Threshold(registrationName='Threshold1', Input=calc_vorticity)

  # Properties modified on threshold1
  threshold1.Scalars = ['POINTS', 'omega_abs']
  threshold1.LowerThreshold = 0.5
  threshold1.UpperThreshold = 8.5

  # show data in view
  threshold1Display = Show(threshold1, renderView1, 'UnstructuredGridRepresentation')

  # trace defaults for the display properties.
  threshold1Display.Representation = 'Surface'

  # hide data in view
  Hide(calc_vorticity, renderView1)

  # update the view to ensure updated data information
  renderView1.Update()

  # create a new 'Clip'
  logger.info(f"Clipping thresholded data using {size} processes...")
  clip1 = Clip(registrationName='Clip1', Input=threshold1)

  # Properties modified on clip1.ClipType
  clip1.ClipType.Origin = [713.1999969482422, -14.44927978515625, 15.0]
  clip1.ClipType.Normal = [0.0, 0.0, -1.0]

  # show data in view
  clip1Display = Show(clip1, renderView1, 'UnstructuredGridRepresentation')

  # trace defaults for the display properties.
  clip1Display.Representation = 'Surface'

  # hide data in view
  Hide(threshold1, renderView1)

  # set scalar coloring
  ColorBy(clip1Display, ('POINTS', 'omega_abs'))

  # show color bar/color legend
  clip1Display.SetScalarBarVisibility(renderView1, True)

  # get color transfer function/color map for 'omega_abs'
  uLUT = GetColorTransferFunction('omega_abs')
  uLUT.ApplyPreset('Rainbow Desaturated', True)
  uLUT.EnableOpacityMapping = 1
  uLUT.UseOpacityControlPointsFreehandDrawing = 1
  uPWF = GetOpacityTransferFunction('omega_abs')

  # Properties modified on uPWF
  n_points = 100
  x_center = 0.20
  steepness = 20
  x = np.linspace(0,1,n_points)
  opacity = (1+np.tanh(steepness*(x-x_center)))/2

  pwf_points = [0.0]*x.shape[0]*4
  for i in range(x.shape[0]):
     pwf_points[i*4] = threshold1.LowerThreshold + x[i]*(threshold1.UpperThreshold -threshold1.LowerThreshold)
     pwf_points[i*4 + 1] = opacity[i]
     pwf_points[i*4 + 2] = 0.5
     pwf_points[i*4 + 3] = 0.0

  uPWF.Points = pwf_points

  # update the view to ensure updated data information
  renderView1.Update()

  # toggle interactive widget visibility (only when running from the GUI)
  HideInteractiveWidgets(proxy=clip1.ClipType)


  ########################################### LOAD STLs ###########################################
  
  logger.info(f"Loading and transforming STLs...")

  # create STL Readers
  spinner_stl = STLReader(registrationName='spinner.stl', FileNames=['paraview/stl/DTU_10MW_hub_spinner_scaled.stl'])
  tower_stl = STLReader(registrationName='tower.stl', FileNames=['paraview/stl/DTU_10MW_tower_offshore_scaled.stl'])
  blade_1_stl = STLReader(registrationName='blade_1.stl', FileNames=['paraview/stl/DTU10MW_blade_1.stl'])

  # show data in view
  spinner_stl_Display = Show(spinner_stl, renderView1, 'GeometryRepresentation')
  spinner_stl_Display.Representation = 'Surface'

  tower_stl_Display = Show(tower_stl, renderView1, 'GeometryRepresentation')
  tower_stl_Display.Representation = 'Surface'

  blade_1_stl_Display = Show(blade_1_stl, renderView1, 'GeometryRepresentation')
  blade_1_stl_Display.Representation = 'Surface'

  # update the view to ensure updated data information
  renderView1.Update()

  # Color STLs black
  for display in [spinner_stl_Display, tower_stl_Display, blade_1_stl_Display]:
      display.AmbientColor = [0.0, 0.0, 0.0]
      display.DiffuseColor = [0.0, 0.0, 0.0]

  # Create all 3 blade STL objects with the corresponding angle
  azimuth_blade_1 = 42 # in deg
  rotor_center = [-7.07, 0.0, 119.0] # in m

  Hide(blade_1_stl, renderView1)

  for i in range(3):

    if i == 0:
      transform_input = blade_1_stl
      azimuth_offset = azimuth_blade_1
    else:
      transform_input = transform_blade
      azimuth_offset = 0.0
     
    transform_blade = Transform(registrationName=f'Transform_blade_{i+1}', Input=transform_input)

    # Properties modified on transform1.Transform
    transform_blade.Transform = 'RotateAroundOriginTransform'
    transform_blade.Transform.Originofrotation = rotor_center
    transform_blade.Transform.Rotate = [120.0+azimuth_offset, 0.0, 0.0]

    # show data in view
    transform_blade_Display = Show(transform_blade, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    transform_blade_Display.Representation = 'Surface'

    # change solid color
    transform_blade_Display.AmbientColor = [0.0, 0.0, 0.0]
    transform_blade_Display.DiffuseColor = [0.0, 0.0, 0.0]

    # update the view to ensure updated data information
    renderView1.Update()
  
  ##################################################################################################

  logger.info(f"Getting layout and saving screenshot...")

  # get layout
  layout1 = GetLayout()

  # layout/tab size in pixels
  layout1.SetSize(1457, 739)

  # Camera placement
  renderView1.CameraPosition = [-404.9666414577699, -455.0907779130123, 252.58339882694676]
  renderView1.CameraFocalPoint = [542.0855373817037, 158.04690416302589, -34.06012187234646]
  renderView1.CameraViewUp = [0.1686295180815408, 0.1915570846957671, 0.9668867404895188]
  renderView1.CameraParallelScale = 1144.1029297516338

  # Save screenshot
  SaveScreenshot(os.path.abspath(args.outfile_path), viewOrLayout=renderView1, location=16, ImageResolution=[1457, 739],
                 OverrideColorPalette='WhiteBackground')

  #####################################################################################################
  
  logger.info("Cleaning up ParaView pipeline and memory...")
 
  Delete(calc_vorticity)
  del calc_vorticity
  
  Delete(redistribute)
  del redistribute

  Delete(reader)
  del reader

  # Disconnect the ParaView session
  Disconnect()
  
  # Force Python garbage collection
  gc.collect()
  
  logger.info("Cleanup complete. Exiting.")

if __name__ == "__main__":
    main()