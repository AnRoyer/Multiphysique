// Gmsh project created on Wed Apr 27 15:48:08 2016
nbelem = 2;

Point(1) = {0, 0, 0, 0.1};
Point(2) = {2, 0, 0, 1.0};
Point(3) = {2, 1, 0, 1.0};
Point(4) = {1, 1, 0, 1.0};
Point(5) = {1, 2, 0, 1.0};
Point(6) = {0, 2, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
Line Loop(7) = {6, 1, 2, 3, 4, 5};
Plane Surface(8) = {7};
Physical Line("t1") = {5};
Physical Line("t2") = {2};
Physical Line("border") = {6, 4, 3, 1};
Physical Surface("plane") = {8};
Transfinite Line {6, 1} = 2*nbelem+1 Using Progression 1;
Transfinite Line {5, 2, 4, 3} = nbelem+1 Using Progression 1;
