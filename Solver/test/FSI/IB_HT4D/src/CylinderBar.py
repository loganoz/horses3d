import numpy as np
from stl import mesh
from math import pi

def body_coordinates():
    n_points_cylinder = 61
    n_points_bar = 200

    radius = 0.5
    x_cylinder_center = 0
    y_cylinder_center = 0

    alpha_cylinder_noend = np.linspace(pi,0.21,n_points_cylinder//2)
    x_cylinder_noend = np.cos(alpha_cylinder_noend)*radius
    y_cylinder_noend = np.sin(alpha_cylinder_noend)*radius
    cylinder_surface_upper = np.vstack((x_cylinder_noend + x_cylinder_center, y_cylinder_noend + y_cylinder_center)).T
    cylinder_surface_lower = np.vstack((x_cylinder_noend + x_cylinder_center, -y_cylinder_noend + y_cylinder_center)).T

    x_bar_upper = np.linspace(0.489898, 4.5, n_points_bar//2)
    y_bar_upper = np.full_like(x_bar_upper, 0.1)
    x_bar_lower = x_bar_upper
    y_bar_lower = np.full_like(x_bar_upper, -0.1)
    bar_surface_upper = np.vstack((x_bar_upper, y_bar_upper)).T
    bar_surface_lower = np.vstack((x_bar_lower, y_bar_lower)).T

    # coordinates = np.vstack((cylinder_surface_upper, bar_surface_upper, bar_surface_lower[::-1], cylinder_surface_lower[::-1]))
    coordinates = np.vstack((cylinder_surface_lower, bar_surface_lower, bar_surface_upper[::-1], cylinder_surface_upper[::-1]))

    return coordinates

def create_body_stl(filename, length_z=2.0):
    coordinates = body_coordinates()
    
    vertices = []
    faces = []
    
    # Main surface of the airfoil
    for i in range(len(coordinates) - 1):
        p1 = [coordinates[i][0], coordinates[i][1], -length_z]
        p2 = [coordinates[i + 1][0], coordinates[i + 1][1], -length_z]
        p3 = [coordinates[i][0], coordinates[i][1], length_z]
        p4 = [coordinates[i + 1][0], coordinates[i + 1][1], length_z]
        
        vertices.append(p1)
        vertices.append(p2)
        vertices.append(p3)
        vertices.append(p4)
        
        faces.append([len(vertices) - 4, len(vertices) - 3, len(vertices) - 2])
        faces.append([len(vertices) - 3, len(vertices) - 1, len(vertices) - 2])
    
    # Front face
    n_points = len(coordinates) // 2
    for i in range(n_points - 1):
        p1 = [coordinates[i][0], coordinates[i][1], -length_z]
        p2 = [coordinates[i + 1][0], coordinates[i + 1][1], -length_z]
        p3 = [coordinates[-(i + 1)][0], coordinates[-(i + 1)][1], -length_z]
        p4 = [coordinates[-(i + 2)][0], coordinates[-(i + 2)][1], -length_z]
        
        vertices.extend([p1, p2, p3, p4])
        faces.append([len(vertices) - 4, len(vertices) - 3, len(vertices) - 2])
        faces.append([len(vertices) - 3, len(vertices) - 1, len(vertices) - 2])
    
    # Back face
    for i in range(n_points - 1):
        p1 = [coordinates[i][0], coordinates[i][1], length_z]
        p2 = [coordinates[i + 1][0], coordinates[i + 1][1], length_z]
        p3 = [coordinates[-(i + 1)][0], coordinates[-(i + 1)][1], length_z]
        p4 = [coordinates[-(i + 2)][0], coordinates[-(i + 2)][1], length_z]
        
        vertices.extend([p1, p2, p3, p4])
        faces.append([len(vertices) - 4, len(vertices) - 3, len(vertices) - 2])
        faces.append([len(vertices) - 3, len(vertices) - 1, len(vertices) - 2])
        
    vertices = np.array(vertices)
    faces = np.array(faces)

    body_mesh = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
    for i, f in enumerate(faces):
        for j in range(3):
            body_mesh.vectors[i][j] = vertices[f[j],:]

    body_mesh.save(filename)

# Create and save the STL file
create_body_stl('body.stl')
