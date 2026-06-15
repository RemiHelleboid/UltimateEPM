SetFactory("OpenCASCADE");

Box(1) = {0, 0, 0, 1, 1, 1};

Physical Surface("x_min", 2) = {1};
Physical Surface("x_max", 3) = {2};
Physical Surface("y_min", 4) = {3};
Physical Surface("y_max", 5) = {4};
Physical Surface("z_min", 6) = {5};
Physical Surface("z_max", 7) = {6};

Physical Volume("bulk", 1) = {1};