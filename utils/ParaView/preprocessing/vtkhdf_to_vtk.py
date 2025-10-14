import vtk

input_path = "<input_file_path>.hdf"
output_path = "<output_file_path>.vtk"

polyData_flag = True

# Read the VTKHDF file
print(f"Reading {input_path} VTKHDF file...")
reader = vtk.vtkHDFReader()
reader.SetFileName(input_path)
reader.Update()

# Get the output data
print(f"Getting output data...")
data = reader.GetOutput()

# Convert vtkUnstructuredGrid to vtkPolyData
if polyData_flag:
  print(f"Converting from vtkUnstructuredGrid to vtkPolyData")
  geometryFilter = vtk.vtkGeometryFilter()
  geometryFilter.SetInputData(data)
  geometryFilter.Update()

  polyData = geometryFilter.GetOutput()

# Write to VTK legacy format
print(f"Writing VTK output to {output_path}...")
writer = vtk.vtkDataSetWriter()
writer.SetFileName(output_path)

if polyData_flag:
  writer.SetInputData(polyData)
else:
  writer.SetInputData(data)
writer.Write()

print(f"Done!")
