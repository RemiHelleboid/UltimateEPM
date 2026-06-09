= Parser Library

*Source:* `src/Parser` \
*CMake target:* `libfile` \
*Kind:* compiled library

== Responsibility

The Parser library handles project file abstractions, Gmsh mesh input, and VTK
output through `file`, `msh_file`, and `vtkWriter`.

== Dependencies

The target privately depends on Gmsh and `libmesh`. Optional STF/HDF5 support
is enabled by the `ENABLE_ST_VERSION` build setting.

== Supported Formats

_To document:_ accepted Gmsh versions and entities, VTK output structure, STF
support, required fields, and unsupported constructs.

== Error Handling

_To document:_ malformed-input behavior, diagnostics, partial-read semantics,
and exceptions or error return values.

== Validation

_To document:_ parser fixtures, round-trip tests, and compatibility samples.
