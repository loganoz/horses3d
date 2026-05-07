#!/usr/bin/env bash

sim_name="$1"

# rm -r RESULTS
# mkdir RESULTS

if [ -z "$sim_name" ]; then
  echo "Usage: ./run.sh sim_name"
  exit 1
fi

sbatch --job-name="$sim_name" --export=ALL,SIM_NAME="$sim_name" script.slurm