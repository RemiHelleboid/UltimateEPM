= Physics Library

*Source:* `src/Physics` \
*CMake target:* `libphysics` \
*Kind:* static library

== Responsibility

The Physics library defines material properties, physical constants, and
statistical relations shared by device and transport simulations. Material
data are loaded from YAML parameter files.

== Dependencies

The target publicly depends on `yaml-cpp` and `uepm::common`.

== Material Model

_To document:_ the material schema, supported parameter sets, defaults,
validation rules, and units used by `materials.hpp`.

== Physical Statistics

_To document:_ distribution functions, approximations, valid temperature and
energy ranges, and conventions for carrier populations.

== Validation

_To document:_ parameter-file validation and numerical reference cases.
