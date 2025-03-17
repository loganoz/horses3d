// Inputs

height = 150.0;
width = 125.0;
thickness = 25.0;

step_height = 0.2;
step_depth = 0.6;
gridsize = 5.0;

// Point coordinates
Point(1) = {-(width-thickness)/2,height-thickness,0,gridsize};
Point(2) = {0,height-thickness,0,gridsize};
Point(3) = {0,0,0,gridsize};
Point(4) = {thickness,0,0,gridsize};
Point(5) = {thickness,height-thickness,0,gridsize};
Point(6) = {(width+thickness)/2,height-thickness,0,gridsize};
Point(7) = {(width+thickness)/2,height,0,gridsize};
Point(8) = {-(width-thickness)/2,height,0,gridsize};

Line(9) = {1,2};
Line(10) = {2,3};
Line(11) = {3,4};
Line(12) = {4,5};
Line(13) = {5,6};
Line(14) = {6,7};
Line(15) = {7,8};
Line(16) = {8,1};

Line Loop(17) = {9,10,11,12,13,14,15,16};

Plane Surface(19) = 17;

#Transfinite Line{11} = 5;
#Transfinite Line{10,12} = 25;
#Transfinite Line{15} = 30;
#Transfinite Line{16,14} = 5;
#Transfinite Line{9,13} = 10;
#Recombine Surface{19};

Physical Line("boundary") = {9,10,11,12,13,14,15,16};
Physical Surface("interior") = {19};
