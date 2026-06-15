SetFactory("OpenCASCADE");

// Unit sphere centered at the origin
Sphere(1) = {0, 0, 0, 1};

// Physical groups
Physical Surface("boundary", 2) = {1};
Physical Volume("bulk", 1) = {1};
