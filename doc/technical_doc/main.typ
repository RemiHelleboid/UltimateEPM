#set document(
  title: "UltimateEPM Technical Documentation",
  author: "UltimateEPM contributors",
)

#set page(
  paper: "a4",
  margin: (x: 2.2cm, y: 2cm),
  numbering: "1",
)
#set text(size: 10.5pt)
#set par(justify: true, leading: 0.65em)
#set heading(numbering: "1.1")
#set outline(indent: auto)

#align(center)[
  #v(3cm)
  #text(size: 24pt, weight: "bold")[UltimateEPM]
  #v(0.4cm)
  #text(size: 16pt)[Technical Documentation]
  #v(1.2cm)
  #text(size: 11pt)[Architecture, numerical methods, and library reference]
  #v(1fr)
  #text(size: 9pt)[UltimateEPM contributors]
  #v(1cm)
]

#pagebreak()
#outline(title: [Contents])
#pagebreak()

= Introduction

UltimateEPM is a C++20 scientific computing project for semiconductor
simulation. Its capabilities include empirical pseudopotential band-structure
calculations, Brillouin-zone processing, device meshing, Poisson solvers, and
full-band and analytical Monte Carlo transport.

This document describes the internal architecture and numerical
responsibilities of the libraries under `src/`. User-facing command-line
instructions remain in the project `README.md` and generated API details remain
the responsibility of Doxygen.

== Document organization

Each source library has a dedicated chapter. A chapter should explain:

- the responsibility and boundaries of the library;
- its public types and principal entry points;
- the mathematical and numerical models it implements;
- its inputs, outputs, units, and invariants;
- its dependencies and interactions with other libraries;
- relevant validation and test coverage.

= Current Scope

The present revision focuses exclusively on the analytical Monte Carlo (AMC)
carrier-transport model. Device electrostatics, circuit coupling, avalanche
detection, and the other UltimateEPM libraries will be documented in later
revisions.

#include "libraries/amc.typ"

= Documentation Conventions

Physical quantities must be documented with their units. When a quantity is
stored in a unit different from the one used by an equation or input file, both
the storage unit and the conversion point must be stated.

Public type and function names are written as code, for example
`uepm::poisson`. File paths are relative to the repository root. Mathematical
symbols should be defined on first use and used consistently throughout a
chapter.

Statements about implemented behavior should be traceable to source code or a
test. Planned behavior must be explicitly labeled as such.
