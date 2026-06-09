= Brillouin-Zone Mesh Library

*Source:* `src/BZ_MESH` \
*CMake target:* `lib_bzmesh` \
*Kind:* static library

== Responsibility

The Brillouin-zone mesh library represents reciprocal-space vertices,
tetrahedra, states, and spatial indices. It computes quantities including Fermi
levels, dielectric data, electron-phonon interactions, and impact-ionization
data on a Brillouin-zone discretization.

== Dependencies

The target publicly depends on `libepp`, `uepm::common`, Eigen, `rapidcsv`,
`fmt`, Gmsh, and OpenMP.

== Mesh Representation

_To document:_ reciprocal-coordinate conventions, tetrahedron orientation,
band and state indexing, symmetry operations, and octree acceleration.

== Physical Calculations

_To document:_ Fermi-level solution, dielectric mesh construction,
electron-phonon coupling models, deformation potentials, overlap integrals,
and impact-ionization rates.

== Numerical Algorithms

_To document:_ tetrahedral interpolation and integration, iso-energy surface
construction, convergence criteria, and parallel decomposition.

== Validation

Geometry primitives and dielectric-function behavior have tests under
`tests/bz_mesh`. Reference calculations and mesh-refinement studies remain to
be documented.
