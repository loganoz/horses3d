// rectangle_refined.geo

lc_global  = 1;
// lc_refined = 0.002;

lg_x1 = -10;
lg_x2 = 25;
lg_y1 = -10;
lg_y2 = 10;
lg_z  = 0;
lg_zex= 1;

lr_x1 = -1.5;
lr_x2 = 5.5;
lr_y1 = -1.5;
lr_y2 = 1.5;



Point(1) = {lg_x1, lg_y1, lg_z, lc_global};
Point(2) = {lg_x1, lg_y2, lg_z, lc_global};
Point(3) = {lg_x2, lg_y2, lg_z, lc_global};
Point(4) = {lg_x2, lg_y1, lg_z, lc_global};

Point(5)  = {lg_x1, lr_y1, lg_z, lc_global};
Point(6)  = {lg_x1, lr_y2, lg_z, lc_global};
Point(7)  = {lr_x1, lg_y2, lg_z, lc_global};
Point(8)  = {lr_x2, lg_y2, lg_z, lc_global};
Point(9)  = {lg_x2, lr_y2, lg_z, lc_global};
Point(10) = {lg_x2, lr_y1, lg_z, lc_global};
Point(11) = {lr_x2, lg_y1, lg_z, lc_global};
Point(12) = {lr_x1, lg_y1, lg_z, lc_global};
Point(13) = {lr_x1, lr_y1, lg_z, lc_global};
Point(14) = {lr_x1, lr_y2, lg_z, lc_global};
Point(15) = {lr_x2, lr_y2, lg_z, lc_global};
Point(16) = {lr_x2, lr_y1, lg_z, lc_global};


Line(1) = {1, 5};
Line(2) = {5, 6};
Line(3) = {6, 2};
Line(4) = {2, 7};
Line(5) = {7, 8};
Line(6) = {8, 3};
Line(7) = {3, 9};
Line(8) = {9, 10};
Line(9) = {10, 4};
Line(10) = {4, 11};
Line(11) = {11, 12};
Line(12) = {12, 1};
Line(13) = {13, 14};
Line(14) = {14, 15};
Line(15) = {15, 16};
Line(16) = {16, 13};
Line(17) = {13, 5};
Line(18) = {14, 6};
Line(19) = {14, 7};
Line(20) = {15, 8};
Line(21) = {15, 9};
Line(22) = {16, 10};
Line(23) = {16, 11};
Line(24) = {13, 12};



Curve Loop(1) = {13, 14, 15, 16};
Plane Surface(1) = {1};

Curve Loop(2) = {1, -17, 24, 12};
Plane Surface(2) = {2};

Curve Loop(3) = {2, -18, -13, 17};
Plane Surface(3) = {3};

Curve Loop(4) = {3, 4, -19, 18};
Plane Surface(4) = {4};

Curve Loop(5) = {19, 5, -20, -14};
Plane Surface(5) = {5};

Curve Loop(6) = {20, 6, 7, -21};
Plane Surface(6) = {6};

Curve Loop(7) = {-15, 21, 8, -22};
Plane Surface(7) = {7};

Curve Loop(8) = {-23, 22, 9, 10};
Plane Surface(8) = {8};

Curve Loop(9) = {-24, -16, 23, 11};
Plane Surface(9) = {9};

Transfinite Surface {1,2,3,4,5,6,7,8,9};


/*mid*/
Transfinite Line {13,15,2,8} = 55 Using Progression 1;
Transfinite Line {-14,16,-5,11} =105 Using Progression 1;
Transfinite Line {-12,-17,-18,4} = 16 Using Progression 0.8;
Transfinite Line {10,-21,-22,-6} = 28 Using Progression 0.865;
Transfinite Line {1,-24,-23,-9,-3,-19,-20,7} = 16 Using Progression 0.8;

Recombine Surface {1,2,3,4,5,6,7,8,9};

Extrude {0, 0, lg_zex} {
  Surface{1,2,3,4,5,6,7,8,9};Layers{1}; Recombine;
}

Physical Surface("inlet") = {55,77,99};
Physical Surface("top") = {103,125,147};
Physical Surface("bottom") = {67,221,199};
Physical Surface("outlet") = {151,173,195};
Physical Surface("back") = {1,2,3,4,5,6,7,8,9};
Physical Surface("front")={68,90,112,222,46,134,156,178,200};
Physical Volume("fluid") = {1,2,3,4,5,6,7,8,9};

// Point(100) = {0, 0.03, 0, lc_global};
// Point(101) = {0, -0.03, 0, lc_global};
