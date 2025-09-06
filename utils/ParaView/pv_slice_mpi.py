from paraview.simple import *
from paraview import servermanager

from mpi4py import MPI

import os
import gc
import ast
import logging
import argparse

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

def compute_div_lamb(reader, logger, with_omega=False):

  # Compute vorticity if necessary
  if not with_omega:
      
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

    last_item = calc_omega_z

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
  parser.add_argument('slice_origin', type=parse_list, help='Slice plane origin, [x,y,z]. List must be provided either between quotation marks or with no blank spaces.')
  parser.add_argument('slice_normal', type=parse_list, help='Slice plane normal vector, [x,y,z]. List must be provided either between quotation marks or with no blank spaces.')
  parser.add_argument('outfile_fields', type=parse_list, help='List of fields to write in output file: ["field_name_1","field_name_2",...,"field_name_N"]. List must be provided either between quotation marks or with no blank spaces.')
  parser.add_argument('outfile_path', type=parse_list_or_string, help='Path or list of paths to the output *.csv file(s) containing the slice data. List must be provided either between quotation marks or with no blank spaces.')
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
    reader.PointArrayStatus = args.outfile_fields
    reader.UpdatePipeline()

  elif ext == ".vtk":
    logger.info(f"Loading input {args.infile_path} into master node...")
    reader = LegacyVTKReader(registrationName=os.path.basename(args.infile_path), FileNames=[args.infile_path])
    
  elif ext == ".tec":
    logger.info(f"Loading input {args.infile_path} into master node...")
    reader = TecplotReader(registrationName=os.path.basename(args.infile_path), FileNames=[args.infile_path])

  elif ext == ".pvtu":
    logger.info(f"Reading partitioned input {args.infile_path} into {size} processes...")
    reader = XMLPartitionedUnstructuredGridReader(registrationName=os.path.basename(args.infile_path), FileName=[args.infile_path])
    reader.TimeArray = 'None'
    redistribute = reader

  reader.UpdatePipeline()

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

  ######################### INSERT HERE CUSTOM CODE TO COMPUTE ADDITIONAL FILEDS #########################

  # Example of computation of additional field: Lamb vector divergence
  if 'Lamb' in args.outfile_fields or 'divLamb' in args.outfile_fields:
    calc_divLamb = compute_div_lamb(redistribute, logger, with_omega=args.omega_in_file)
  else:
    calc_divLamb = redistribute

  ########################################################################################################

  # create a new 'Slice'
  logger.info(f"Slicing domain")
  if isinstance(args.slice_origin[0], list):
    for i, origin in enumerate(args.slice_origin):
      name = 'Slice'+str(i+1)
      normal = get_input_list(args.slice_normal, len(args.slice_origin), i)
      outfields = get_input_list(args.outfile_fields, len(args.slice_origin), i)
      
      if isinstance(args.outfile_path, list):
        outfile = args.outfile_path[i]
      else:
        outfile = os.path.splitext(args.outfile_path)[0] + "_" + str(i) + os.path.splitext(args.outfile_path)[1]

      slice_domain(calc_divLamb, name, origin, normal, outfields, outfile, logger)

  else:
    slice_domain(calc_divLamb, 'Silce1', args.slice_origin, args.slice_normal, args.outfile_fields, args.outfile_path, logger)

  
  logger.info("Cleaning up ParaView pipeline and memory...")

  Delete(calc_divLamb)
  del calc_divLamb

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