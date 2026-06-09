= Device Mesh Library

*Source:* `src/DEVICE_MESH` \
*CMake target:* `libmesh` \
*Alias:* `uepm::mesh` \
*Kind:* static library

== Responsibility

The Device Mesh library represents device geometry and discretization. It
defines vertices, one-, two-, and three-dimensional elements, regions,
contacts, datasets, spatial functions, grids, and device-level mesh data. It
also provides quadtree and octree spatial indexing and mesh output.

== Dependencies

The target publicly depends on Gmsh, `plog`, `rapidcsv`, and `uepm::common`.
OpenMP support is conditional in the current build definition.

== Geometry and Topology

_To document:_ entity ownership, element orientation, region and boundary
identifiers, adjacency, and dimensional assumptions.

== Spatial Indexing

_To document:_ quadtree and octree construction, lookup complexity, stopping
criteria, and behavior at cell boundaries.

== Data and I/O

_To document:_ datasets attached to meshes, interpolation rules, supported
writer formats, and Gmsh conventions.

== Validation

_To document:_ topology invariants, interpolation reference cases, and mesh
round-trip tests.
