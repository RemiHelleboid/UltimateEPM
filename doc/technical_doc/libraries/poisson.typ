= Poisson Library

*Source:* `src/Poisson` \
*CMake target:* `libpoisson` \
*Alias:* `uepm::poisson` \
*Kind:* static library

== Responsibility

The Poisson library implements finite-element infrastructure and electrostatic
Poisson solvers for two- and three-dimensional device meshes.

== Dependencies

The target privately depends on `libmesh` and publicly depends on `plog`,
POSIX threads, Eigen, and `libphysics`.

== Governing Equation

_To document:_ the exact Poisson equation and sign convention, permittivity
model, charge-density definition, boundary conditions, and unit system.

== Finite-Element Formulation

_To document:_ weak formulation, element basis functions, quadrature, global
assembly, boundary-condition treatment, and linear solver configuration.

== Solver Interfaces

_To document:_ inputs and outputs of the 2D and 3D solvers, mesh requirements,
convergence reporting, and failure modes.

== Validation

Finite elements and both solver dimensions have direct tests under
`tests/poisson`. The analytical cases and accepted tolerances should be listed
here.
