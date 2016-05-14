// Gmsh project to generate a honeycomb geometry : 

pc_rect = 1;
//points_largeur = 21;
gc_rect = 2;
//points_longueur = 41;
e = 0.1*pc_rect; // Espace entre les hexagones

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

// Définition des surfaces

Line Loop(31) = {5, 6, -19, 29};
Plane Surface(32) = {31};
Line Loop(33) = {20, 21, 11, 12};
Plane Surface(34) = {33};
Line Loop(35) = {10, 23, 24, -9};
Plane Surface(36) = {35};
Line Loop(37) = {7, 8, 26, 27};
Plane Surface(38) = {37};
Line Loop(39) = {15, 16, 17, 18, 13, 14};
Plane Surface(40) = {39};
Line Loop(46) = {30, -12, -11, 22, -10, 9, 25, -8, -7, 28, -6, -5};
Plane Surface(47) = {39, 46};

// Définition des Physical group : 

Physical Line(48) = {29, 30, 20, 21, 22, 23, 24, 25, 26, 27, 28, 19}; // Contour
Physical Line(49) = {5, 6, 13, 14, 15, 7, 8, 9, 10, 16, 17, 18, 11, 12}; // Chemin
Physical Surface("Resist") = {32, 34, 40, 38, 36}; // Cellule résistance
Physical Surface("Conductive") = {46}; // Cellule conductrice

// Définition maillage 


Transfinite Line {19, 21, 23, 27} = 8 Using Progression 1;
Transfinite Line {29, 20, 26, 24} = 5 Using Progression 1;
Transfinite Line {25, 12, 30} = 3 Using Progression 1;
Transfinite Line {28, 14, 6, 13, 18, 11, 17, 22, 16, 10, 15, 7} = 15 Using Progression 1;
Transfinite Line {5, 12, 8, 9} = 10 Using Progression 1;
