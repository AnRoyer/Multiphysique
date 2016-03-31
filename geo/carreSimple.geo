cl__1 = 1;

NN=10;

Point(1) = {0, 0, 0, 1};
Point(2) = {1, 0, 0, 1};
Point(3) = {1, 1, 0, 1};
Point(4) = {0, 1, 0, 1};
Line(1) = {1, 2};
Transfinite Line {1} = NN+1 Using Progression 1;
Line(2) = {2, 3};
Transfinite Line {2} = NN+1 Using Progression 1;
Line(3) = {3, 4};
Transfinite Line {3} = NN+1 Using Progression 1;
Line(4) = {4, 1};
Transfinite Line {4} = NN+1 Using Progression 1;
Line Loop(6) = {3, 4, 1, 2};
Plane Surface(6) = {6};
Physical Line("left") = {4};
Physical Line("right") = {2};
Physical Line("top") = {3};
Physical Line("bottom") = {1};
Physical Surface("plane") = {6};
