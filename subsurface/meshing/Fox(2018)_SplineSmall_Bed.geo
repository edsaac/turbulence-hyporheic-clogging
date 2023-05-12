/*
Purpose: Mesh to be used in a filtration model of subsurface water under sand dunes following dimensions in	Fox et al. 2018 (WRR)
*/

// Options for export
Mesh.SaveAll = 0;
Mesh.MshFileVersion = 2.2;

convertToMeter = 1.0/100.0;

//#### Dimensions ####################

Hd = 1.50 * convertToMeter;
Yd = -20.0 * convertToMeter;
Lee = 5.00 * convertToMeter;
Stoss = 10.00 * convertToMeter;
meshSize = 0.30 * convertToMeter;
Width = 29.00 * convertToMeter;

// To get a cell of 0.20cm in the interface:
K_Y = 0.99456572240492;

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

nElementsX = Ceil((Lee+Stoss)/meshSize);
Printf("%f",nElementsX);
nElementsY = Ceil(Fabs(Yd)/meshSize);

Transfinite Line{1} = nElementsX;
Transfinite Line{3} = nElementsX;

Transfinite Line{-2} = nElementsY Using Progression K_Y;
Transfinite Line{4} = nElementsY Using Progression K_Y;

Line Loop(1) = {1 ... 4};
Plane Surface(1) = {1};
Transfinite Surface {1} = {1,4,5,8};
Recombine Surface {1};
extrusionName[] = Extrude {0, 0, Width} {
    Surface{1};
    Layers{1};
    Recombine;
  };

Physical Surface("back") = {1};
Physical Surface("front") = extrusionName[0];
Physical Volume("bed") = extrusionName[1];
Physical Surface("top") = extrusionName[2];
Physical Surface("right") = extrusionName[3];
Physical Surface("bottom") = extrusionName[4];
Physical Surface("left") = extrusionName[5];
