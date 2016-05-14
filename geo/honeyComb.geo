// Gmsh project to generate a honeycomb geometry : 

pc_rect = 1;
//points_largeur = 21;
gc_rect = 2;
//points_longueur = 41;
e = 0.11*pc_rect; // Espace entre les hexagones

// Calcul des longueurs et largeurs :

l = pc_rect/2 - e/2;
L = (gc_rect + e)/3;
aspect_ratio = gc_rect/pc_rect;

// Définition du rectangle de base

Point(1) = {0, 0, 0, 1.0};
Point(2) = {pc_rect, 0, 0, 1.0};
Point(3) = {pc_rect, gc_rect, 0, 1.0};
Point(4) = {0, gc_rect, 0, 1.0};

// Définition des points 

Point(5) = {0, L, 0, 1.0};
Point(6) = {l, L/2, 0, 1.0};
Point(7) = {l, 0, 0, 1.0};
Point(8) = {pc_rect - l, 0, 0, 1.0};
Point(9) = {pc_rect - l, L/2, 0, 1.0};
Point(10) = {pc_rect, L, 0, 1.0};
Point(11) = {pc_rect, gc_rect - L, 0, 1.0};
Point(12) = {pc_rect - l, gc_rect - L/2, 0, 1.0};
Point(13) = {pc_rect - l, gc_rect, 0, 1.0};
Point(14) = {l, gc_rect, 0, 1.0};
Point(15) = {l, gc_rect - L/2, 0, 1.0};
Point(16) = {0, gc_rect - L, 0, 1.0};
Point(17) = {l + e/2, L/2 + e/2, 0, 1.0};
Point(18) = {pc_rect - e/2, L + e/2, 0, 1.0};
Point(19) = {pc_rect - e/2, gc_rect - L - e/2, 0, 1.0};
Point(20) = {l + e/2, gc_rect - L/2 - e/2, 0, 1.0};
Point(21) = {e/2, gc_rect - L - e/2, 0, 1.0};
Point(22) = {e/2,  L + e/2, 0, 1.0};

// Définition des lignes:

Line(5) = {14, 15};
Line(6) = {15, 16};
Line(7) = {5, 6};
Line(8) = {6, 7};
Line(9) = {9, 8};
Line(10) = {9, 10};
Line(11) = {11, 12};
Line(12) = {12, 13};
Line(13) = {20, 21};
Line(14) = {21, 22};
Line(15) = {22, 17};
Line(16) = {17, 18};
Line(17) = {18, 19};
Line(18) = {19, 20};
Line(19) = {4, 16};
Line(20) = {13, 3};
Line(21) = {3, 11};
Line(22) = {11, 10};
Line(23) = {10, 2};
Line(24) = {2, 8};
Line(25) = {8, 7};
Line(26) = {7, 1};
Line(27) = {1, 5};
Line(28) = {5, 16};
Line(29) = {4, 14};
Line(30) = {14, 13};

// Définition des surfaces:

Line Loop(31) = {19, -6, -5, -29};
Plane Surface(32) = {31};
Line Loop(33) = {20, 21, 11, 12};
Plane Surface(34) = -{33};
Line Loop(35) = {23, 24, -9, 10};
Plane Surface(36) = -{35};
Line Loop(37) = {26, 27, 7, 8};
Plane Surface(38) = -{37};
Line Loop(39) = {15, 16, 17, 18, 13, 14};
Plane Surface(40) = {39};
Line Loop(41) = {30, -12, -11, 22, -10, 9, 25, -8, -7, 28, -6, -5};
Plane Surface(42) = {39, 41};

// Définition des physicals:

Physical Line("top") = {29, 30, 20};
Physical Line("bottom") = {26, 25, 24};
Physical Line("left") = {27, 28, 19};
Physical Line("right") = {23, 22, 21};
Physical Surface("ext") = {32, 34, 36, 38, 40};
Physical Surface("int") = {42};

// Maillage

Transfinite Line {5, 12, 8, 9} = 5 Using Progression 1;
Transfinite Line {6, 13, 28, 14, 15, 7, 16, 10, 17, 22, 11, 18} = 8 Using Progression 1;
Transfinite Line {30, 25} = 2 Using Progression 1;
Transfinite Line {19, 27, 23, 21} = 4 Using Progression 1;
Transfinite Line {29, 20, 24, 26} = 3 Using Progression 1;
