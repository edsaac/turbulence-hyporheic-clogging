/*
Purpose: Mesh to be used in a 3D LES/RANS flow simulation
	of surface water over dunes following dimensions in
	Fox et al. 2018 (WRR). It takes number of elements in each
    direction as an "argument" rather than the size of the element.
*/

// Options for export
Mesh.SaveAll = 0;
Mesh.MshFileVersion = 2.2;

convertToMeter = 1.0/100;

//#### Dimensons ####################

Hd = 1.50 * convertToMeter;
Yd = 9.75 * convertToMeter;
Lee = 5.00 * convertToMeter;
Stoss = 10.00 * convertToMeter;
Width = 29.00 * convertToMeter;

//### Mesh Size #####################
nElementsX = 50;
nElementsY = 50;
nElementsZ = 40;

meshSizeX = Abs((Lee+Stoss)*convertToMeter/nElementsX);
meshSizeY = Abs(Yd*convertToMeter/nElementsY);
meshSizeZ = Abs(Width*convertToMeter/nElementsZ);

Printf("median mesh size:");
Printf("  X = %.2E",meshSizeX);
Printf("  Y = %.2E",meshSizeY);
Printf("  Z = %.2E",meshSizeZ);

// This growth ratio gives a del0 = 0.01mm
progY = 1.09252660067602;
progZ = 0.1;

// Dummy value
meshSize  = 1.0;

Point(1) = {0,0,0,meshSize};
Point(2) = {Stoss/2,Hd/2,0,meshSize};
Point(3) = {Stoss/2+Lee,-Hd/2,0,meshSize};
Point(4) = {Stoss+Lee,0,0,meshSize};

Point(5) = {Stoss+Lee,Yd,0,meshSize};
Point(6) = {Stoss/2+Lee,Yd,0,meshSize};
Point(7) = {Stoss/2,Yd,0,meshSize};
Point(8) = {0,Yd,0,meshSize};

Spline(1) = {1,2,3,4};
Line(2) = {4,5};
Line(3) = {5,6,7,8};
Line(4) = {8,1};

Transfinite Line{1} = nElementsX;
Transfinite Line{3} = nElementsX;
Transfinite Line{2} = nElementsY Using Progression progY;
Transfinite Line{-4} = nElementsY Using Progression progY;

Point(9) = {0,0,Width,meshSize};
Point(10) = {Stoss/2,Hd/2,Width,meshSize};
Point(11) = {Stoss/2+Lee,-Hd/2,Width,meshSize};
Point(12) = {Stoss+Lee,0,Width,meshSize};

Point(13) = {Stoss+Lee,Yd,Width,meshSize};
Point(14) = {Stoss/2+Lee,Yd,Width,meshSize};
Point(15) = {Stoss/2,Yd,Width,meshSize};
Point(16) = {0,Yd,Width,meshSize};

Spline(5) = {9,10,11,12};
Line(6) = {12,13};
Line(7) = {13,14,15,16};
Line(8) = {16,9};

Transfinite Line{5} = nElementsX;
Transfinite Line{7} = nElementsX;
Transfinite Line{6} = nElementsY Using Progression progY;
Transfinite Line{-8} = nElementsY Using Progression progY;

Line Loop(1) = {1 ... 4};
Plane Surface(1) = {1};
Transfinite Surface {1} = {1,4,5,8};
Recombine Surface {1};

Line Loop(2) = {5 ... 8};
Plane Surface(2) = {2};
Transfinite Surface {2} = {9,12,13,16};
Recombine Surface {2};

Line(9)  = {1,9};
Line(10) = {4,12};
Line(11) = {5,13};
Line(12) = {8,16};

Transfinite Line{9}  = nElementsZ Using Bump progZ;
Transfinite Line{10} = nElementsZ Using Bump progZ;
Transfinite Line{11} = nElementsZ Using Bump progZ;
Transfinite Line{12} = nElementsZ Using Bump progZ;

Line Loop(3) = {-1,9,5,-10};
Surface(3) = {3};
Transfinite Surface {3} = {1,4,12,9};
Recombine Surface {3};

Line Loop(4) = {-2,10,6,-11};
Plane Surface(4) = {4};
Transfinite Surface {4} = {4,12,5,13};
Recombine Surface {4};

Line Loop(5) = {-3,11,7,-12};
Plane Surface(5) = {5};
Transfinite Surface {5} = {8,5,13,16};
Recombine Surface {5};

Line Loop(6) = {-4,12,8,-9};
Plane Surface(6) = {6};
Transfinite Surface {6} = {1,9,16,8};
Recombine Surface {6};

Physical Surface("back") = {1};
Physical Surface("front") = {2};
Physical Surface("bottom") = {3};
Physical Surface("right") = {4};
Physical Surface("top") = {5};
Physical Surface("left") = {6};

Surface Loop(1) = {1 ... 6};
Volume(1) = {1};

Physical Volume("bed") = {1};
Transfinite Volume "*";

