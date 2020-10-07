//--------------------------------------------------------------------
// Geometry parameters - inputs
//--------------------------------------------------------------------
hp  =   1;              // parent channel height
hd1 = 0.5;              // 1st daughter channel height == 0.5
hd2 = 0.5;              // 2nd daughter channel height == 0.5
L   =  2;               // channel lengths
L2  =   1;              // Daughter channel contraction lengths == 1
nCells =  30;           // number of cells in transverse direction
xCells = 60;            // number of cells in streamwise direction
//--------------------------------------------------------------------
// Geometry parameters - calculated
//--------------------------------------------------------------------
midp = hp*Sqrt(3)/2;
lp = hp/nCells;
ld1 = hd1/nCells;
ld2 = hd2/nCells;
dm = 2*midp/hp;
dx = (L2^2/(1+dm^2))^(1/2);
dy = dm*dx;
dx1 = 0.5*((hp-hd1)^2/(1+(1/dm)^2))^(1/2);
dy1 = dx1/dm;
dx2 = 0.5*((hp-hd2)^2/(1+(1/dm)^2))^(1/2);
dy2 = dx2/dm;

//--------------------------------------------------------------------
// Points
//--------------------------------------------------------------------
// Junction - triangle
//Point(1) = {0, hp/2, 0, lp};
//Point(2) = {0, -hp/2, 0, lp};
//Point(3) = {midp, 0, 0, lp};
// Junction - contractions
//Point(4) = {dx+dx1, hp/2+dy-dy1, 0, ld1};
//Point(5) = {midp+dx-dx1, dy+dy1, 0, ld1};
//Point(6) = {dx+dx2, -(hp/2+dy-dy2), 0, ld2};
//Point(7) = {midp+dx-dx2, -(dy+dy2), 0, ld2};

// Junction - triangle
Point(1) = {0, hp/2, 0, lp};
Point(2) = {0, -hp/2, 0, lp};
Point(3) = {-midp, 0, 0, lp};
// Junction - contractions
Point(4) = {-(dx+dx1), hp/2+dy-dy1, 0, ld1};
Point(5) = {-(midp+dx-dx1), dy+dy1, 0, ld1};
Point(6) = {-(dx+dx2), -(hp/2+dy-dy2), 0, ld2};
Point(7) = {-(midp+dx-dx2), -(dy+dy2), 0, ld2};


//--------------------------------------------------------------------
// Lines
//--------------------------------------------------------------------
// Junction - triangle
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 1};
// Junction - contractions
Line(4) = {1, 4};
Line(5) = {4, 5};
Line(6) = {5, 3};
Line(7) = {3, 7};
Line(8) = {7, 6};
Line(9) = {6, 2};



//--------------------------------------------------------------------
// Surfaces
//--------------------------------------------------------------------
Line Loop(1) = {1, 2, 3};
Plane Surface(1) = {1};
Line Loop(2) = {4, 5, 6, 3};
Plane Surface(2) = {2};
Line Loop(3) = {7, 8, 9, 2};
Plane Surface(3) = {3};


//--------------------------------------------------------------------
// Channels
//--------------------------------------------------------------------
///pS[] = Extrude {-L, 0, 0} {
 // Line{1};
//  Layers{xCells};
//  Recombine;
//};
//d1S[] = Extrude {dx*L/L2, dy*L/L2, 0} {
//  Line{5};
//  Layers{xCells};
//  Recombine;
//};
//
//d2S[] = Extrude {dx*L/L2, -dy*L/L2, 0} {
//  Line{8};
//  Layers{xCells};
//  Recombine;
//};

pS[] = Extrude {L, 0, 0} {
  Line{1};
  Layers{xCells};
  Recombine;
};
d1S[] = Extrude {-(dx*L/L2), dy*L/L2, 0} {
  Line{5};
  Layers{xCells};
  Recombine;
};

d2S[] = Extrude {-(dx*L/L2), -dy*L/L2, 0} {
  Line{8};
  Layers{xCells};
  Recombine;
};


//--------------------------------------------------------------------
// Unit depth
//--------------------------------------------------------------------
zV[] = Extrude {0, 0, -1} {
  Surface{1,2,3,pS[1],d1S[1],d2S[1]};
  Layers{1};
  Recombine;
};


//--------------------------------------------------------------------
// Physical surfaces
//--------------------------------------------------------------------
Physical Surface("outlet") = {zV[21]};
Physical Surface("inlet1") = {zV[27]};
Physical Surface("inlet2") = {zV[33]};
Physical Surface("walls") = {zV[7],zV[9],zV[13],zV[15],zV[20],zV[22],zV[26],zV[28],zV[32],zV[34]};
Physical Surface("frontAndBack") = {1,2,3,pS[1],d1S[1],d2S[1],zV[0],zV[5],zV[11],zV[17],zV[23],zV[29]};
Physical Volume("internalMesh") = {zV[1],zV[6],zV[12],zV[18],zV[24],zV[30]};


