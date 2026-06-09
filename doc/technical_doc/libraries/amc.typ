= Analytical Monte Carlo Transport

*Source:* `src/AMC` \
*CMake target:* `lib_amc` \
*Alias:* `uepm::amc` \
*Language:* C++20

== Scope and Status

The analytical Monte Carlo (AMC) model simulates semiclassical electron and
hole transport in silicon. A numerical particle represents one or more
physical carriers and evolves through alternating deterministic drift and
stochastic scattering.

This chapter documents the transport model implemented by:

- `valley_model.hpp` and `valley_model.cpp`;
- `amc_material_model.hpp` and `amc_material_model.cpp`;
- `amc_scattering_model.hpp` and `amc_scattering_model.cpp`;
- `intervalley_phonon.hpp` and `intervalley_phonon.cpp`;
- `amc_transport_kernel.hpp` and `amc_transport_kernel.cpp`;
- `particle_amc.hpp` and `particle_amc.cpp`;
- the transport loops in `bulk_amc_simulation.cpp` and
  `device_amc_simulation.cpp`.

The current model is silicon-specific. Electrons use six anisotropic,
nonparabolic Delta valleys. Holes use isotropic parabolic heavy-hole and
light-hole bands. The transport kernel does not read these parameters from the
project YAML files: they are currently compiled into
`amc_material_model.cpp` and `intervalley_phonon.cpp`.

Self-consistent Poisson coupling, induced-current computation, external circuit
coupling, avalanche detection, and quench detection are outside the present
scope.

== Semiclassical State

Each `particle_amc` owns a `particle_state`. The transport-relevant state is

$ S = (t, bold(r), bold(k)_v, bold(v), E, gamma, nu), $

where:

- $t$ is the particle time in seconds;
- $bold(r)$ is the position in device-mesh coordinates;
- $bold(k)_v$ is the wave vector in the principal frame of the active valley
  or band, in $"m"^(-1)$;
- $bold(v)$ is the group velocity in the global device frame, in $"m"/"s"$;
- $E$ is kinetic energy measured from the active valley or band edge, in eV;
- $gamma$ is the transformed kinetic energy, in eV;
- $nu$ is the active electron-valley or hole-band index.

The particle also stores the local electric field, doping and impurity
concentrations, its containing mesh element, its previous position, and a flag
indicating collection by a contact.

=== Carrier charge and numerical weight

The signed carrier charge is

$ q_c = cases(
  -q_e & "for electrons",
  +q_e & "for holes",
) $

with $q_e = 1.602176634 times 10^(-19) "C"$.

The numerical particle weight $w > 0$ is independent of the microscopic
transport trajectory. It scales charge and observables when one simulated
particle represents several physical carriers. Scattering probabilities are
not weight-dependent.

=== Coordinate and unit conventions

Three coordinate conventions meet in the transport code:

- the device mesh stores position in micrometres;
- the transport kernel stores velocity in $"m"/"s"$ and electric field in
  $"V"/"m"$;
- `particle_state::local_k` is expressed in the active valley's local
  principal-axis frame.

Consequently, the position update explicitly multiplies displacement in metres
by $10^6$ before adding it to the mesh position.

In a device simulation, the interpolated electric field is treated as
$"V"/"cm"$
and multiplied by $10^2$ before being passed to the transport kernel.

The impurity fields are stored in $"cm"^(-3)$. For device transport, the positive
scattering-centre density is set to the absolute value of the interpolated net
doping concentration:

$ N_I = abs(N_"dop"). $

== Band and Valley Model

=== General ellipsoidal dispersion

Each `valley_model` defines a principal frame with transverse masses $m_t$ and
longitudinal mass $m_l$. In that frame,

$ bold(k)_v = (k_"t1", k_"t2", k_l). $

The transformed energy is the ellipsoidal quadratic form

$ gamma(bold(k)_v)
  = frac(planck.reduce^2, 2 q_e)
    (frac(k_"t1"^2, m_t)
     + frac(k_"t2"^2, m_t)
     + frac(k_l^2, m_l)), $

where division by $q_e$, numerically equal to the joule-per-electronvolt
conversion, expresses $gamma$ in eV.

For a parabolic band,

$ gamma(E) = E. $

For a Kane band,

$ gamma(E) = E (1 + alpha E), $

where $alpha$ is the nonparabolicity in $"eV"^(-1)$. The inverse relation used by
the code is

$ E(gamma) = frac(sqrt(1 + 4 alpha gamma) - 1, 2 alpha). $

The parabolic relation is used when the band is explicitly parabolic or when
$alpha = 0$.

=== Group velocity

The group velocity in the valley frame follows from
$bold(v) = planck.reduce^(-1) nabla_k E$:

$ v_i = frac(planck.reduce k_i, m_i (1 + 2 alpha E)). $

For parabolic bands the denominator correction is one. The local velocity is
rotated into the global device frame before it is stored in the particle.

Each valley owns an orthonormal rotation matrix $R$ whose columns orient the
local frame in global coordinates:

$ bold(k)_v = R^T bold(k)_g, quad
   bold(v)_g = R bold(v)_v. $

The constructor verifies $R^T R = I$ to a tolerance of $10^(-10)$.

=== Constant-energy state generation

After a scattering event, the outgoing direction is sampled isotropically in
the scaled ellipsoidal coordinate system. A unit vector $bold(u)$ is drawn
uniformly on the sphere and converted to a valley-frame wave vector:

$ bold(k)_v =
  frac(sqrt(2 q_e gamma(E)), planck.reduce)
  (sqrt(m_t) u_x, sqrt(m_t) u_y, sqrt(m_l) u_z). $

This construction preserves the requested kinetic energy. It assumes
isotropic angular redistribution in transformed momentum space. In
particular, the screened-Coulomb rate is a momentum-relaxation rate, but the
outgoing direction is still sampled isotropically rather than from the
differential screened-Coulomb cross section.

=== Silicon electron valleys

Electrons occupy six equivalent Delta valleys:

#table(
  columns: (1.4fr, 1fr, 1fr, 1fr, 1fr),
  align: (left, center, center, center, center),
  table.header(
    [Valleys], [$m_t / m_0$], [$m_l / m_0$], [$alpha$], [Axis],
  ),
  [$Delta_x^+, Delta_x^-$], [0.19], [0.916], [$0.5 "eV"^(-1)$], [$x$],
  [$Delta_y^+, Delta_y^-$], [0.19], [0.916], [$0.5 "eV"^(-1)$], [$y$],
  [$Delta_z^+, Delta_z^-$], [0.19], [0.916], [$0.5 "eV"^(-1)$], [$z$],
)

All six minima have zero energy offset in the implemented model. The sign in
the valley name identifies the opposite member of an axis pair, but both
members currently use the same mass tensor and rotation matrix. The valley
centre in the full Brillouin zone is not added to `local_k`; transport evolves
only the wave vector relative to the selected minimum.

The conductivity and density-of-states masses exposed by `valley_model` are

$ m_c = frac(3, 1 / m_l + 2 / m_t), quad
  m_d = (m_l m_t^2)^(1/3). $

=== Silicon hole bands

Holes use two isotropic, parabolic bands:

#table(
  columns: (1.4fr, 1fr, 1fr, 1fr),
  align: (left, center, center, center),
  table.header([Band], [$m_t / m_0$], [$m_l / m_0$], [$alpha$]),
  [Heavy hole], [0.87], [0.87], [0],
  [Light hole], [0.24], [0.24], [0],
)

Both band edges have zero offset. Warping, split-off holes, anisotropy, and
nonparabolicity are therefore omitted from the current hole model.

== Initial Distribution

Particle initialization first selects a valley or band externally. Bulk and
device simulations distribute initial particles cyclically over the available
indices.

The initial kinetic energy is sampled from a gamma distribution with shape
$3/2$ and scale $k_B T$:

$ p(E) prop E^(1/2) exp(-E / (k_B T)). $

This is the three-dimensional Maxwell-Boltzmann energy distribution for a
parabolic density of states. The same sampler is also used for Kane electron
valleys; the nonparabolic density-of-states Jacobian is not included in the
initial-energy distribution.

An isotropic transformed-space direction is then drawn and used to construct
`local_k`, $gamma$, kinetic energy, and global velocity consistently.

== Deterministic Free Flight

Between scattering events, the semiclassical acceleration theorem is

$ planck.reduce frac(d bold(k)_g, d t) = q_c bold(F)_g. $

The kernel rotates the electric field into the active valley frame and applies
the exact constant-field wave-vector increment over a step $Delta t$:

$ bold(k)_v^(n+1)
  = bold(k)_v^n + frac(q_c, planck.reduce) R^T bold(F)_g Delta t. $

It then recomputes $gamma$, $E$, and velocity from the updated wave vector.
Position is advanced with the trapezoidal velocity:

$ bold(r)^(n+1)
  = bold(r)^n
  + frac(bold(v)^n + bold(v)^(n+1), 2) Delta t. $

The displacement is converted from metres to micrometres at this point.
Particle time is incremented by $Delta t$.

The field is assumed constant during one call to `drift_particle`. In device
transport it is interpolated once at the beginning of the global time step.
The implementation does not locate an exact mesh-face or contact crossing
during the drift.

== Scattering Framework

At a particle state $(nu, E)$, the kernel constructs a finite list of allowed
channels $j$. Each channel contains:

- a physical mechanism;
- a rate $Gamma_j$ in $"s"^(-1)$;
- the prescribed final kinetic energy;
- a destination band when applicable;
- an intervalley phonon branch and absorption/emission tag when applicable.

The total physical rate is

$ Gamma_"tot"(nu, E) = sum_j Gamma_j(nu, E). $

Given that a real event occurs, channel selection uses the normalized rates:

$ P(j | "event") = frac(Gamma_j, Gamma_"tot"). $

The channel container has a compile-time capacity of 16. The current electron
model can populate acoustic, ten intervalley absorption/emission channels,
impurity, and impact-ionization channels, which fits this capacity.

=== Acoustic deformation-potential scattering

Acoustic scattering is treated as elastic and equiprobable in angle. For
electrons, the implemented rate can be written

$ Gamma_"ac"(E) =
  frac(sqrt(2) k_B T m_d^(3/2) D_"ac"^2,
       pi planck.reduce^4 rho u_"avg"^2)
  sqrt(gamma_J(E)) (1 + 2 alpha_J E_J), $

with:

- silicon mass density $rho = 2329 "kg"/"m"^3$;
- longitudinal sound speed $u_l = 9.0 times 10^3 "m"/"s"$;
- transverse sound speed $u_t = 5.4 times 10^3 "m"/"s"$;
- $u_"avg" = (u_l + 2 u_t) / 3$;
- acoustic deformation potential $D_"ac" = 6.6 "eV"$;
- $gamma_J(E) = E_J (1 + alpha_J E_J)$.

For holes, the same functional dependence is used with:

- sound speed $u_s = 6.6 times 10^3 "m"/"s"$;
- deformation potential $D_"ac" = 5.5 "eV"$;
- overlap factor $1/2$.

The rate evaluation clamps energy to at least $10^(-9) "eV"$ before evaluating
the square-root factor. Applying the event preserves kinetic energy, band or
valley, and randomizes direction.

=== Electron intervalley phonon scattering

The electron model contains five silicon intervalley phonon branches:

#table(
  columns: (1.1fr, 0.65fr, 0.8fr, 1fr, 1.2fr, 0.8fr),
  align: (left, center, center, center, center, center),
  table.header(
    [Branch], [Family], [Order], [$planck.reduce omega$], [Deformation potential],
    [$Z_f$],
  ),
  [`g1_TA`], [$g$], [First], [11.4 "meV"], [$D_1 = 3.0 "eV"$], [1],
  [`g2_LA`], [$g$], [First], [18.8 "meV"], [$D_1 = 3.0 "eV"$], [1],
  [`g3_LO`], [$g$], [Zeroth], [63.2 "meV"], [$D_0 = 3.4 times 10^10 "eV"/"m"$], [1],
  [`f1_TA`], [$f$], [First], [21.9 "meV"], [$D_1 = 3.0 "eV"$], [4],
  [`f2_LA`], [$f$], [Zeroth], [46.3 "meV"], [$D_0 = 3.4 times 10^10 "eV"/"m"$], [4],
)

For a phonon energy $planck.reduce omega$, the equilibrium Bose-Einstein
occupation is

$ N_"op" = frac(1, exp(planck.reduce omega / (k_B T)) - 1). $

Absorption uses the factor $N_"op"$ and

$ E_f = E_i + planck.reduce omega. $

Emission uses $N_"op" + 1$ and

$ E_f = E_i - planck.reduce omega, $

and is disabled when $E_f < 0$.

For a zeroth-order branch, the implemented rate is

$ Gamma_0 =
  frac(sqrt(2) Z_f m_d^(3/2) D_0^2,
       pi rho planck.reduce^2 (planck.reduce omega))
  N_"proc"
  sqrt(gamma_J(E_f))
  (1 + 2 alpha_J E_"f,J"), $

where $N_"proc"$ is the absorption or emission phonon factor.

For a first-order branch,

$ Gamma_1 =
  frac(sqrt(2) Z_f m_d^(5/2) D_1^2,
       pi rho planck.reduce^4 (planck.reduce omega))
  N_"proc"
  sqrt(gamma_J(E_f))
  (1 + 2 alpha_J E_"f,J")
  (gamma_J(E_f) + gamma_J(E_i)). $

The branch multiplicity $Z_f$ is included in the rate. At event application,
a $g$ process switches to the opposite valley on the same axis. An $f$
process selects uniformly among the four valleys on the other two axes. The
outgoing direction is randomized in the destination valley.

=== Hole optical phonon scattering

The hole model uses a single optical phonon energy of $63 "meV"$ and four
transitions:

#table(
  columns: (1fr, 1fr, 1fr),
  align: (left, center, center),
  table.header([Transition], [Destination], [Overlap]),
  [heavy to heavy], [heavy], [0.5],
  [heavy to light], [light], [1.0],
  [light to heavy], [heavy], [1.0],
  [light to light], [light], [0.5],
)

All four use a deformation potential of $8.0 times 10^10 "eV"/"m"$. The rate has
the same density-of-final-states form as the electron zeroth-order intervalley
rate, multiplied by the transition overlap factor. Both absorption and
emission are included, and the outgoing direction is randomized in the
destination band.

=== Ionized-impurity scattering

Impurity scattering is optional. Bulk transport uses a configured uniform
background density. Device transport uses the particle-local density
interpolated from the mesh.

Two rate models are available.

==== Empirical mobility model

The Caughey-Thomas mobility is

$ mu(N_I) =
  mu_"min"
  + frac(mu_0 - mu_"min",
         1 + (N_I / N_"ref")^beta). $

The implementation removes the lattice-limited contribution using Matthiessen's
rule:

$ frac(1, mu_I) = frac(1, mu(N_I)) - frac(1, mu_0), $

then converts mobility to a scalar momentum-relaxation rate:

$ Gamma_I = frac(q_e, m_d mu_I). $

The compiled silicon parameters are:

#table(
  columns: (1fr, 1fr, 1fr, 1fr, 0.8fr),
  align: (left, center, center, center, center),
  table.header(
    [Carrier], [$mu_0$], [$mu_"min"$], [$N_"ref"$], [$beta$],
  ),
  [Electron], [$1417 "cm"^2/("V" "s")$], [$52.2 "cm"^2/("V" "s")$],
  [$9.68 times 10^16 "cm"^(-3)$], [0.68],
  [Hole], [$470.5 "cm"^2/("V" "s")$], [$44.9 "cm"^2/("V" "s")$],
  [$2.23 times 10^17 "cm"^(-3)$], [0.70],
)

This model produces an energy-independent rate for a given carrier type and
impurity density.

==== Screened-Coulomb model

The alternative model uses a Brooks-Herring-type screened-Coulomb
momentum-relaxation rate. Screening is Debye-like:

$ q_s^2 = frac(q_e^2 n_s, epsilon_s k_B T), quad
  epsilon_s = 11.7 epsilon_0. $

With

$ k^2 = frac(2 m_d gamma_J(E), planck.reduce^2), quad
  b = frac(4 k^2, q_s^2), $

the angular momentum-relaxation factor is

$ A(b) = ln(1 + b) - frac(b, 1 + b). $

The implementation forms the rate from $A(b)/(4 k^4)$, impurity density, the
silicon dielectric constant, and the Kane density-of-states correction.

In this function, electron nonparabolicity is fixed internally to
$0.5 "eV"^(-1)$ and hole nonparabolicity to zero, rather than read from the
active `valley_model`. The screening density is currently set equal to the
impurity density.

Both impurity models are applied as elastic, isotropic direction-randomizing
events in the same band or valley.

=== Impact ionization

Impact ionization is optional and uses a threshold power law:

$ Gamma_"II"(E) =
  cases(
    0 & E <= E_"th",
    A ((E - E_"th") / E_"th")^p & E > E_"th".
  ) $

The compiled parameters are:

#table(
  columns: (1fr, 1fr, 1fr, 1fr),
  align: (left, center, center, center),
  table.header([Carrier], [$E_"th"$], [$A$], [$p$]),
  [Electron], [1.12 "eV"], [$1.1 times 10^14 "s"^(-1)$], [2.5],
  [Hole], [1.49 "eV"], [$1.4 times 10^12 "s"^(-1)$], [3.4],
)

When the event is applied, the primary carrier remains in its current band or
valley, its kinetic energy is reduced by the threshold, and its direction is
randomized:

$ E_f = max(0, E_i - E_"th"). $

In device transport, the event may additionally enqueue one electron and one
hole at the event position, with the same numerical weight as the initiating
particle. Those secondary particles receive independent thermal initial
states. Energy and momentum are not partitioned among the three outgoing
carriers, so this is a carrier-multiplication model rather than a
microscopically energy-conserving final-state model.

== Monte Carlo Time Integration

The code provides two distinct transport algorithms. They should not be
described as numerically identical.

=== Fixed-time-step algorithm

`bulk_amc_simulation::run()` and device transport use a global step
$Delta t$:

1. Interpolate or assign the electric field.
2. Drift the particle for the full $Delta t$.
3. Evaluate all rates at the post-drift state.
4. Compute the probability of at least one event,
   $P_"sc" = 1 - exp(-Gamma_"tot" Delta t)$.
5. If an event is accepted, select and apply exactly one channel.

This method reproduces the no-event probability for a constant rate over the
step, but permits at most one event. It is therefore a time-discretized
approximation when $Gamma_"tot" Delta t$ is not small or when the field changes
the rate significantly during the step. A practical accuracy condition is

$ max(Gamma_"tot" Delta t) lt.double 1. $

The fixed-step path does not use the null-collision bound
`gamma_max`.

=== Event-driven null-collision algorithm

`bulk_amc_simulation::run_self_scattering_emc()` implements a
self-scattering, or null-collision, ensemble Monte Carlo algorithm.

For each band or valley $nu$, initialization samples the total rate on a
uniform energy grid from zero to `m_max_energy_eV`. The majorant rate is

$ Gamma_"max,nu"
  = s max_(0 <= E <= E_"max") Gamma_"tot"(nu, E), $

where $s$ is the configured safety factor. The default values are
$E_"max" = 2 "eV"$, $s = 1.2$, and 1000 energy samples.

The candidate free-flight time is exponentially distributed:

$ tau = -frac(ln U, Gamma_"max,nu"), quad U in (0, 1). $

After drifting for $tau$, the real-event acceptance probability is

$ P_"real" = frac(Gamma_"tot"(nu, E), Gamma_"max,nu"). $

If the candidate is rejected, a self-scattering event is recorded and the
physical state is unchanged. If accepted, a real channel is selected from its
normalized rate.

If an observed total rate exceeds its valley majorant, the thread-local kernel
raises that majorant to the observed rate times the safety factor before
acceptance. This prevents an acceptance probability above one for that
candidate, but the already sampled flight time came from the old majorant.
For unbiased null-collision transport, configuration should therefore ensure
the precomputed majorant covers the physically reached energy range.

Because an intervalley event can change $nu$, the next free-flight sample uses
the destination valley's majorant.

=== Random-number streams and parallel execution

The kernel uses `std::mt19937_64`. Constructors may seed from
`std::random_device` or from an explicit 64-bit seed.

Parallel bulk and device runs create one transport kernel per OpenMP thread.
Thread seeds are offset by 7919; electron and hole device streams receive
adjacent base seeds. Results are reproducible for a fixed seed, thread count,
particle ordering, and scheduling behavior, but changing the number of
threads changes the random streams assigned to particles.

== Device Transport Boundary Handling

After all particles drift for one device time step, the simulation updates
mesh ownership and boundaries before scattering:

- a particle entering a contact is flagged and removed;
- a particle remaining inside its previous element keeps that element;
- otherwise, the device searches for a new containing element;
- if no element is found, or the new element is not silicon, the particle is
  returned to its pre-step position;
- reflection negates the complete velocity and local wave vector.

This is a full back-reflection, not specular reflection about a computed
surface normal. The crossing time is not resolved, and the remaining fraction
of the time step is not transported after reflection or contact collection.

The two-dimensional self-consistent implementation may additionally apply a
periodic condition in the out-of-plane direction. The base three-dimensional
device transport leaves this operation as a no-op.

== Observables and Histories

A particle history records snapshots of time, position, local wave vector,
velocity, kinetic energy, transformed energy, and band or valley index. It
also counts:

- acoustic events;
- intervalley absorption;
- intervalley emission;
- impurity events;
- impact-ionization events;
- self-scattering events.

Bulk observables are accumulated after a configurable warm-up fraction. The
code time-weights velocity and kinetic energy. The reported impact-ionization
coefficient is based on

$ alpha_"II" =
  frac("event rate per carrier", "average drift speed"), $

with conversion from inverse metres to inverse centimetres.

The event-driven bulk path also exposes a per-particle raw estimate based on
event count divided by displacement along $x$. That estimate assumes $x$ is
the drift direction and should not be used for arbitrary field orientation.

== Configuration Parameters

The transport kernel is controlled by `amc_transport_config`:

#table(
  columns: (1.5fr, 0.7fr, 2fr),
  align: (left, center, left),
  table.header([Parameter], [Default], [Role]),
  [`m_carrier_type`], [electron], [Select electron valleys or hole bands.],
  [`m_lattice_temperature`], [300 K], [Phonon occupation and initial thermal energy.],
  [`m_max_energy_eV`], [2 eV], [Upper energy used to construct null-collision bounds.],
  [`m_self_scattering_safety_factor`], [1.2], [Multiplier on sampled maximum rates.],
  [`m_gamma_max_energy_samples`], [1000], [Uniform samples used for each majorant.],
  [`m_background_impurity_density_cm_3`], [0], [Uniform bulk impurity density.],
  [`m_enable_impact_ionization`], [false], [Enable threshold power-law channel.],
  [`m_enable_impurity_scattering`], [false], [Enable an impurity channel.],
  [`m_impurity_scattering_model`], [mobility], [Select empirical or screened-Coulomb rate.],
  [`m_impurity_density_source`], [background], [Use configured or particle-local density.],
)

`m_gamma_max_energy_samples` must be at least two in command-line option
validation. The kernel's `initialize()` method itself assumes this condition
when forming the uniform energy grid.

== Implemented Invariants

The transport classes enforce several local invariants:

- particle weight is strictly positive;
- kinetic energy and $gamma$ cannot be set negative through their setters;
- drift and scattering time increments cannot be negative;
- effective masses are positive;
- nonparabolicity is non-negative;
- valley degeneracy is at least one;
- valley rotations are orthonormal;
- final phonon-emission energy cannot be negative;
- channel rates passed to event application cannot be negative.

Direct public access through `particle_amc::state()` can bypass some particle
setter checks. Simulation code must preserve state consistency:

$ gamma = gamma(bold(k)_v), quad
  E = E(gamma), quad
  bold(v) = bold(v)(bold(k)_v, E). $

The transport kernel restores these relations after initialization, drift, and
every real scattering event.

== Validation Coverage

`tests/amc/test_amc_transport_kernel.cpp` currently verifies:

- the sum of constructed channel rates equals the reported total rate;
- a per-valley null-collision majorant is positive and no larger than the
  global maximum;
- a zero-duration fixed-step scattering call leaves the particle unchanged.

These tests exercise basic bookkeeping but do not yet validate the physical
rate equations or transport observables.

Recommended model-validation cases are:

- energy and velocity consistency for every valley orientation;
- Maxwell-Boltzmann equilibrium at zero field with impact ionization disabled;
- detailed-balance ratios for phonon absorption and emission;
- low-field electron and hole mobility versus temperature and doping;
- velocity-field curves for bulk silicon;
- valley populations under fields along the principal crystal axes;
- convergence with fixed time step;
- convergence with null-collision energy-grid resolution and safety factor;
- impact-ionization coefficient versus electric field;
- statistical reproducibility and confidence intervals across independent
  random seeds.

== Model Limitations

The following limitations are part of the current implementation and should be
kept visible when interpreting results:

- only silicon is supported;
- material and scattering parameters are compiled into C++ sources;
- electron transport uses six analytical Delta valleys with no explicit
  Brillouin-zone valley centres;
- hole transport uses isotropic parabolic heavy- and light-hole bands;
- the initial Kane-electron energy sampler uses the parabolic
  Maxwell-Boltzmann density of states;
- all real scattering events use isotropic outgoing directions in transformed
  momentum space;
- acoustic scattering is elastic and uses an equipartition approximation;
- no carrier-carrier scattering or Pauli blocking is included;
- the screened-Coulomb differential angular distribution is not sampled;
- impact ionization does not conserve the complete outgoing carrier energy and
  momentum microscopically;
- device fields are constant during a global transport step;
- device boundary crossing times and surface normals are not resolved;
- the fixed-step algorithm permits at most one event per step;
- null-collision accuracy depends on the configured energy range containing
  all reached carrier energies.

These limitations do not make the model unusable, but they define the level of
physical fidelity that validation must establish for each intended operating
regime.
