// Gmsh project created on Sat Feb 13 16:51:54 2016

alpha = 0.8;
beta = 1-alpha;
period = 3;
point = 21;

// Carré de base
Point(1) = {0, 0, 0, 1.0};
Point(2) = {0, 1, 0, 1.0};

Line(1) = {1, 2};

Transfinite Line {1} = point Using Progression 1;
Physical Line("Left") = {1};

Physical Line("Top") = {};
Physical Line("Bottom") = {};
Physical Line("Inner") = {};
Physical Surface("Alpha") = {};
Physical Surface("Beta") = {};

p = 3;
l = 2;
i = 1;
// Différentes régions
For(1 : period)
	//Zone alpha
	Point(p) = {alpha*1/period + (i-1)/period, 0, 0, 1.0};
	Point(p+1) = {alpha*1/period + (i-1)/period, 1, 0, 1.0};

	Line(l) = {p, p+1};
	Line(l+1) = {p-2,p};
	Line(l+2) = {p-1,p+1};

	//Zone beta
	Point(p+2) = {i/period, 0, 0, 1.0};
	Point(p+3) = {i/period, 1, 0, 1.0};

	Line(l+3) = {p,p+2};
	Line(l+4) = {p+1,p+3};
	Line(l+5) = {p+2, p+3};

	//Surface
	If(i == 1)
		Line Loop(l+6) = {-(l-1), (l+1), (l), -(l+2)};
	Else
		Line Loop(l+6) = {-(l-3), (l+1), (l), -(l+2)};
	EndIf
	Plane Surface(l+6) = {l+6};
	Physical Surface("Alpha") += {l+6};

	Line Loop(l+7) = {-(l), (l+3), (l+5), -(l+4)};
	Plane Surface(l+7) = {l+7};
	Physical Surface("Beta") += {l+7};

	//Transfinite et physical
	Transfinite Line {l} = point Using Progression 1;
	Transfinite Line {l+1} = alpha*point/period Using Progression 1;
	Transfinite Line {l+2} = alpha*point/period Using Progression 1;
	
	Transfinite Line {l+3} = beta*point/period Using Progression 1;
	Transfinite Line {l+4} = beta*point/period Using Progression 1;
	Transfinite Line {l+5} = point Using Progression 1;

	Physical Line("Top") += {l+2, l+4};
	Physical Line("Bottom") += {l+1, l+3};
	Physical Line("Inner") += {l};
	If(i == period)
		Physical Line("Right") = {l+5};
	Else
		Physical Line("Inner") += {l+5};
	EndIf

	p = p + 4;
	l = l + 8;
	i = i + 1;
EndFor
