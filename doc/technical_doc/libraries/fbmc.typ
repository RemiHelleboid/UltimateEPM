= Full-Band Monte Carlo Library

*Source:* `src/FBMC` \
*CMake target:* `lib_fbmc` \
*Kind:* static library

== Responsibility

The FBMC library provides particle representation and single-particle
full-band Monte Carlo transport based on electronic structure sampled over the
Brillouin zone.

== Dependencies

The target publicly depends on `libepp`, `lib_bzmesh`, Eigen, `rapidcsv`,
`fmt`, Gmsh, and OpenMP.

== Particle State

_To document:_ position, wave vector, band index, energy, velocity, time, and
the units and invariants of each state component.

== Transport Algorithm

_To document:_ free-flight sampling, field-driven state evolution, scattering
selection, Brillouin-zone boundary handling, and random-number generation.

== Validation

_To document:_ deterministic seeded cases, equilibrium distributions,
transport coefficients, and comparisons against reference results.
