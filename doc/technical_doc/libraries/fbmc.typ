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

The bulk solver uses null-collision free flights. During each flight the
electric field advances the wave vector, the state is folded back into the
stored Brillouin-zone domain when necessary, and energy and group velocity are
interpolated from the full-band mesh. A real event is selected from the eight
electron--phonon channels and, for electrons by default, the Keldysh
impact-ionization channel.

`--maxenergy` is a strict validity boundary. Carrier energy is not clamped and
phonon rates are not extrapolated beyond it. A carrier that crosses the
boundary is discarded while the remaining carriers continue. The run emits a
warning and records the discarded count and an incomplete-run flag.

== Impact Ionization

Electron impact ionization is enabled by default in `fbmc.epm`; use
`--disable-impact-ionization` to remove the channel. Hole impact ionization is
not implemented and hole runs must disable it.

The Keldysh event rate is

$ Gamma_"II"(E) =
  cases(
    0 & E < E_"th",
    P_0 (E - E_"th")^gamma & E >= E_"th".
  ) $

On an accepted event, the tracked electron loses exactly $E_"th"$ and is
redrawn on the remaining-energy surface. The generated electron--hole pair is
not added to the FBMC ensemble.

== Observables

Steady-state observables are accumulated after the configured warm-up
fraction. In addition to velocity, mobility, energy, and impact-ionization
event rate, FBMC exports two ionization-coefficient estimators.

The existing speed-based estimator is

$ alpha_"speed" =
  frac(N_"II", integral abs(v dot hat(E)) dif t). $

The endpoint-displacement estimator is

$ alpha_"endpoint" =
  frac(sum_i N_{"II",i}, sum_i norm(X_{f,i} - X_{0,i})). $

$X_0$ is the carrier position at the end of warm-up and $X_f$ is its final
position, or its position when discarded at the energy boundary. Both
coefficients are exported in inverse centimetres. The endpoint denominator is
also exported in metres.

The output additionally includes `discarded_carriers_over_max_energy` and
`run_complete`. Results with `run_complete = 0` are diagnostic because
discarding high-energy carriers biases ensemble averages.

== Validation

Deterministic seeded runs are used for parameter comparisons. Production
validation should use several seeds, larger ensembles, longer trajectories,
fields outside the fitting set, and an energy window for which no carriers are
discarded.
