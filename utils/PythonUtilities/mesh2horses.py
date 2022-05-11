#!/usr/bin/env python3

import meshio
from sys import argv, exit


def check_args(args):

    # We only read the first argument
    if len(args) < 3 or args[1] == "--help" or args[1] == "-h":

        # Output the help
        print(("Description: Utility to convert meshio-supported "
               "meshes to HORSES3D format."))
        print("Usage:       mesh2horses [-h|--help] input output")
        print("  -h  --help  print this help and exit")

        exit(0)


def face_indices(ID):
    r"""
           v
    3----------2
    |\     ^   |\
    | \    |   | \
    |  \   |   |  \
    |   7------+---6
    |   |  +-- |-- | -> u
    0---+---\--1   |
     \  |    \  \  |
      \ |     \  \ |
       \|      w  \|
        4----------5
    """

    if ID.lower() == "front":       # v-
        indices = [0, 1, 5, 4]

    elif ID.lower() == "back":      # v+
        indices = [3, 2, 6, 7]

    elif ID.lower() == "bottom":    # w-
        indices = [0, 1, 2, 3]

    elif ID.lower() == "right":     # u+
        indices = [1, 2, 6, 5]

    elif ID.lower() == "top":       # w+
        indices = [4, 5, 6, 7]

    elif ID.lower() == "left":      # u-
        indices = [0, 3, 7, 4]

    else:
        indices = []

    return indices


# Future work
def get_connectivities(mesh):

    pass


def write_mesh(mesh, file_name):

    with open(file_name, 'w') as f:

        n_points = len(mesh.points)
        n_elements = len(mesh.cells_dict['hexahedron'])
        # The last number is the order of the curved faces -> 0 so far
        f.write(f"{n_points:12d}{n_elements:12d}{0:12d}\n")

        for point in mesh.points:
            f.write(f"{point[0]: .16f}  {point[1]: .16f}  {point[2]: .16f}\n")

        # boundaries = {"BC_name": array([face IDs])
        boundaries = {}
        for field in mesh.field_data:
            try:
                boundaries[field] = mesh.cell_sets_dict[field]['quad']
            except Exception:
                pass

        # boundary_elems = array([[face1_node1, face1_node2...]])
        boundary_elems = mesh.get_cells_type("quad")

        # faces = {"BC_name": array([[face1_node1, face1_node2...]])}
        faces = {}
        for field in boundaries:
            for face in boundaries[field]:
                if field in faces:
                    faces[field].append(sorted(boundary_elems[face]))
                else:
                    faces[field] = [sorted(boundary_elems[face])]

        # Loop over the elements and write their nodes and BCs
        elements = mesh.get_cells_type("hexahedron")
        for elem in elements:

            # Nodes and element order (only 0 so far)
            f.write(("{:12d}"*8 + '\n').format(*(elem+1)))
            f.write("           0"*8 + '\n')

            # Now get the BCs that the faces belong to
            for pos in ("front", "back", "bottom", "right", "top", "left"):
                ind = face_indices(pos)
                nodes = sorted(elem[ind])
                face = "---"
                for BC in faces:
                    if nodes in faces[BC]:
                        face = BC
                        break
                # Write one BC at a time
                f.write(f" {face:15s}")
            # Next line for the next element
            f.write('\n')


if __name__ == "__main__":

    # Read the mesh
    check_args(argv)
    mesh = meshio.read(argv[1])

    # And save to the output file
    write_mesh(mesh, argv[2])
