// Gmsh project created on Mon Sep 27 13:13:45 2021

//Import the geometry module
SetFactory("OpenCASCADE");
//Define your domain
xmax =  20;
xmin = -20;
ymax =  10;
ymin = -10;
zmax =  1.0;
zmin =  0.0;
shift = 0.0;
R = 0.5;

//Define some points to divide the domain in regions with 4 sides
Point(1) = {xmin, ymin, zmin, 3.0}; //{x_pos, y_pos, z_pos, mesh_tolerance}
Point(2) = { xmax, ymin, zmin, 3.0};
Point(3) = { xmax,  ymax, zmin, 3.0};
Point(4) = {xmin,  ymax, zmin, 3.0};
Point(5) = {-Cos(45*Pi/180)*R+shift,  Sin(45*Pi/180)*R+shift, zmin, 0.01};
Point(6) = {-Cos(45*Pi/180)*R+shift, -Sin(45*Pi/180)*R+shift, zmin, 0.01};
Point(7) = { Cos(45*Pi/180)*R+shift, -Sin(45*Pi/180)*R+shift, zmin, 0.01};
Point(8) = { Cos(45*Pi/180)*R+shift,  Sin(45*Pi/180)*R+shift, zmin, 0.01};
Point(9) = {0+shift, 0+shift, zmin, 1.0};

//Define the lines that connect the points
Line(1) = {1, 2}; //Line from Point1 to Point2 
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {4, 5};
Line(6) = {1, 6};
Line(7) = {7, 2};
Line(8) = {8, 3};

//Create the cylinder (a circle in 2D) with 4 arcs using Point9 as the middle point
Circle(9)  = {5, 9, 8}; //Circle arc from Point5 to Point8 and Point9 as the center
Circle(10) = {8, 9, 7};
Circle(11) = {7, 9, 6};
Circle(12) = {6, 9, 5};

//Create the Curve Loops (groups of 4 lines to create a closed region)
//In this step, it is very important to consider the orientation of each line.
//The 4 lines must create a circular loop.
Curve Loop(1) = {6, 12, -5, 4}; //Loop: Line6, Circle12, Line5(reversed), Line4
Curve Loop(2) = {1, -7, 11, -6};
Curve Loop(3) = {2, -8, 10, 7};
Curve Loop(4) = {3, 5, 9, 8};

//Create one surface for each Curve Loop
Plane Surface(1) = {1}; //Surface 1 is inside Curve Loop 1
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};

//To create a high-quality mesh, it’s important to select a seed in each line.
//This way, you can control the accuracy in each region.
//1) Create a seed of 50 points in Line5(reversed), Line6(reversed), Line7 and Line8
//   Then, apply a progression of 1.1 instead of a uniform seed. The progression
//   allows you to gather the points near one end of each line.

Transfinite Curve{-5, -6, 7, 8} = 50 Using Progression 1.1;

//2) Create a uniform seed of 15 points in the external boundaries of your domain.
//Also, convert each surface to transfinite to create hexahedral elements.

//Bottom
Transfinite Curve{11, 1} = 15;
Transfinite Surface{2};

//Left
Transfinite Curve{4, 12} = 15;
Transfinite Surface{1};

//Up
Transfinite Curve{3, 9} = 15;
Transfinite Surface{4};

//Right
Transfinite Curve{2, 10} = 15;
Transfinite Surface{3};


//Recombine the surfaces to create the desired pattern
Recombine Surface {1, 2, 3, 4};

//Extrude the 2D mesh into a 3D mesh with 1 layer of height=zmax in the z-direction
Extrude {0, 0, zmax} {
  Surface{1,2,3,4};
  Recombine;
  Layers{1}; 
}

//Finally, assign labels to each surface and volume. These names will be used to define
// the boundary conditions in the control file used by HORSES3D.
//Before this step, load the *.geo file in the Gmsh GUI and check manually the number
// of each surface and volume.
//First, let’s create a new Physical Surface with the number 33 (use an empty number
// for the new surface) based on the original Surface8.
Physical Surface("inlet", 33) = {8};
//Repeat the operation for every boundary
Physical Surface("outlet", 34) = {14};
Physical Surface("up", 35) = {18};
Physical Surface("down", 36) = {10};
Physical Surface("wall1", 37) = {13, 17, 9, 20};
Physical Surface("wall2", 38) = {2, 3, 4, 1};
Physical Surface("cylinder", 39) = {19, 6, 16, 12};
//The whole volume is limited by Surfaces 1, 2, 3 and 4.
Physical Volume("fluid", 40) = {4, 3, 2, 1};
// uncoment next line to assing a order 3
// Mesh.ElementOrder = 3; 
