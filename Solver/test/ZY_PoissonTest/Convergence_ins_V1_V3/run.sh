#!/bin/bash

input_file="ConvergenceBDF.control"
output_file="convergenceTest1103__GaussPoint_order10.txt"
program="  /home/yu/AbbasHorses3d/horses3d/Solver/bin/horses3d.slr"  
parallelRun= "OMP_NUM_THREADS=1"

length=1

# 
for i in {1..10}; do
    new_line=$(printf "                               Mesh file name = MESH/meshTest1024_Per_length1_%dx%dx%d_STR.msh" $length $length $length)

    sed -i "7s|.*|$new_line|" "$input_file"

    result=$($parallelRun $program "$input_file")

    last_line=$(echo "$result" | tail -n 1)

    echo "$last_line" >> "$output_file"

    length=$((length + 1))
done
