SetFactory("OpenCASCADE");

L  = 4.0;
H  = 1.0;
lc = 0.1;

Point(1) = {0, 0, 0, lc};
Point(2) = {L, 0, 0, lc};
Point(3) = {L, H, 0, lc};
Point(4) = {0, H, 0, lc};

Line(1) = {1, 2}; // bottom
Line(2) = {2, 3}; // kathode
Line(3) = {3, 4}; // top
Line(4) = {4, 1}; // anode

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Curve {1, 3} = 101;
Transfinite Curve {2, 4} = 11;
Transfinite Surface {1};

// No Recombine Surface: keep triangular mesh

Physical Surface("Si", 1) = {1};

Physical Curve("kathode", 2) = {2};
Physical Curve("anode", 3) = {4};
Physical Curve("insulation", 4) = {1, 3};