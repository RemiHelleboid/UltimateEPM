= Empirical Pseudopotential Library

*Source:* `src/EPP` \
*CMake target:* `libepp` \
*Kind:* static library

== Responsibility

The EPP library implements empirical pseudopotential electronic-structure
calculations. It contains material and symmetry-point definitions, Hamiltonian
construction, local and nonlocal pseudopotential terms, spin-orbit coupling,
band-structure calculations, dielectric-function calculations, and
Brillouin-zone mesh-file support.

== Dependencies

The target publicly depends on `uepm::common`, Eigen, `yaml-cpp`, `rapidcsv`,
`fmt`, OpenMP, and Gmsh.

== Mathematical Model

_To document:_ the plane-wave basis, Hamiltonian matrix elements,
pseudopotential form factors, nonlocal corrections, spin-orbit terms, and
eigenvalue conventions.

== Main Components

_To document:_ `Material`, `Hamiltonian`, `Pseudopotential`,
`BandStructure`, `DielectricFunction`, and the symmetry-point representation.

== Inputs and Outputs

_To document:_ material parameter files, reciprocal-space paths, basis-size
selection, energy units, band indexing, and exported results.

== Validation

The Bessel-function implementation has direct coverage under `tests/epp`.
Additional reference-band and dielectric-function validation remains to be
documented.
