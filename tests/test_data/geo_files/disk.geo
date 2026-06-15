//+
SetFactory("OpenCASCADE");
//+
SetFactory("OpenCASCADE");
//+
Disk(1) = {0, 0, 0, 1, 1};
//+
Physical Surface("Circle", 2) = {1};
//+
Physical Curve("edge", 3) = {1};
