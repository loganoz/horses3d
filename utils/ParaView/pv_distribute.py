from paraview.simple import *
from mpi4py import MPI

import os
import gc
import logging
import argparse

# Initialise MPI 
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def main():
    
  parser = argparse.ArgumentParser(description="Redistributes a dataset in *.hdf, *.vtk, *.vtkhdf or *.tec format and writes it to *.pvtu format for quick parallel loading using MPI.")
  
  # Define postional arguments
  parser.add_argument('infile_path', type=str, help='Path to the input *.hdf, *.vtk, *.vtkhdf or *.tec file.')
  parser.add_argument('outfile_path', type=str, help='Path to the output *.pvtu file. A directory with the same base name will be created containing all the partition files.')
  
  # Parse arguments
  args = parser.parse_args()

  # Check formats of input and output files
  _, ext_in = os.path.splitext(args.infile_path)
  _, ext_out = os.path.splitext(args.outfile_path)
  if ext_in not in [".hdf", ".vtk", ".vtkhdf", ".tec"]:
    raise ValueError(f"Unsupported input file format: {ext_in}. Supported formats are: .hdf, .vtk, .vtkhdf, .tec")
  if ext_out != ".pvtu":
    raise ValueError(f"Unsupported output file format: {ext_out}. Supported format is: .pvtu")
  
  # Configure logging 
  logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
  logger = logging.getLogger()

  # Load your data
  _, ext = os.path.splitext(args.infile_path)

  if ext in [".hdf", ".vtkhdf"]:
    logger.info(f"Loading input {args.infile_path} into master node...")
    reader = VTKHDFReader(registrationName=os.path.basename(args.infile_path), FileName=[args.infile_path])
    
  elif ext == ".tec":
    logger.info(f"Loading input {args.infile_path} into master node...")
    reader = TecplotReader(registrationName=os.path.basename(args.infile_path), FileNames=[args.infile_path])

  elif ext == ".vtk":
    logger.info(f"Loading input {args.infile_path} into master node...")
    reader = LegacyVTKReader(registrationName=os.path.basename(args.infile_path), FileNames=[args.infile_path])
    
  reader.UpdatePipeline()

  # Redistribute the data
  logger.info(f"Redistributing data into {size} partitions...")
  redistribute = RedistributeDataSet(Input=reader)
  redistribute.UpdatePipeline()

  # Save the data in parallel format
  logger.info(f"Saving partitioned data to {args.outfile_path}...")
  SaveData(args.outfile_path, proxy=redistribute)

  logger.info(f"Done!")

  # Force Python garbage collection
  gc.collect()
  
  logger.info("Cleanup complete. Exiting.")

if __name__ == "__main__":
    main()
