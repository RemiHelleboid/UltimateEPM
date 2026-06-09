= Common Library

*Source:* `src/Common` \
*CMake target:* `uepm_common` \
*Alias:* `uepm::common` \
*Kind:* header-only interface library

== Responsibility

The Common library provides reusable numerical and infrastructure utilities.
Its headers cover physical constants, unit conversion, vector operations,
numerical integration, statistical functions, floating-point comparison, CSV
helpers, and generic container utilities.

== Dependencies

The interface target propagates `rapidcsv`, `fmt::fmt`, and the C++20 language
requirement to consumers.

== Public API

_To document:_ group the public headers by numerical, physical, data-export,
and general-purpose responsibilities. State units and numerical assumptions for
all physical and integration helpers.

== Numerical Considerations

_To document:_ tolerances, convergence behavior, invalid-input handling, and
the precision expected by the integration and comparison routines.

== Validation

_To document:_ direct unit tests and the higher-level libraries that exercise
these utilities indirectly.
