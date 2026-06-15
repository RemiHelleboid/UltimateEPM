SetFactory("OpenCASCADE");

Rectangle(1) = {0, 0, 0, 1, 1, 0};

Physical Surface("bulk", 13) = {1};

Physical Curve("edge_0", 9)  = {1};
Physical Curve("edge_1", 10) = {2};
Physical Curve("edge_2", 11) = {3};
Physical Curve("edge_3", 12) = {4};