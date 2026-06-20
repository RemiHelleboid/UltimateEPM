# Electron–phonon and FBMC runtime guide

This guide describes the practical workflow for computing electron–phonon
scattering data and using it in full-band Monte Carlo simulations.

It is intentionally organized around *why* each runtime step exists. For the
precise list of command-line options, use:

```bash
./build/apps/elph.epm --help
./build/apps/fbmc.epm --help
```

## The three runtime artifacts

The workflow produces three conceptually different kinds of data:

1. **Rate kernels** contain the expensive band-structure and phase-space part
   of electron–phonon scattering.
2. **Scattering rates** are obtained by applying a selected deformation
   potential to those kernels.
3. **FBMC results** are transport trajectories and averages generated using
   the scattering rates.

The intended flow is:

```text
BZ mesh and phonon physics
          |
          v
 reusable rate-kernel CSV       expensive, normally computed once
          |
          | apply A, B and energy threshold from a YAML profile
          v
 phonon-rate CSV                cheap to regenerate while fitting
          |
          v
 bulk FBMC simulation
```

This separation matters because the tetrahedron search, final-state density of
states and overlap calculations are expensive. The deformation potential is
only an initial-energy-dependent multiplier in the current model, so changing
it does not require repeating those calculations.

## Preparing the mesh

Both programs expect a Gmsh BZ mesh containing the band energies and
gradients. In the examples below, that file is named `bz.msh`.

The carrier and band counts must be chosen consistently:

```text
electrons: --carrier electron -c <conduction bands> -v 0
holes:     --carrier hole     -c 0 -v <valence bands>
```

The code shifts conduction energies to the conduction-band minimum. For holes,
valence energies are converted into positive hole kinetic energies measured
from the valence-band maximum.

The initial carrier energy used by the deformation-potential model follows
this same convention.

## Generating a reusable kernel

Kernel generation is the expensive step. A typical electron command is:

```bash
./build/apps/elph.epm \
  -f bz.msh \
  -m Si \
  --carrier electron \
  --phonon-params remi-2026 \
  -c 2 -v 0 \
  -T 300 \
  -E 0.3 \
  -j 16 \
  --export-kernels \
  --kernels-out electron_kernels_300K.csv \
  --rates-only
```

For holes:

```bash
./build/apps/elph.epm \
  -f bz.msh \
  -m Si \
  --carrier hole \
  --phonon-params remi-2026 \
  -c 0 -v 2 \
  -T 300 \
  -E 0.3 \
  -j 16 \
  --export-kernels \
  --kernels-out hole_kernels_300K.csv \
  --rates-only
```

The kernel includes the effects of:

- the BZ mesh and electronic bands;
- final-state tetrahedron DOS;
- phonon dispersion and Bose occupation;
- electron or hole overlap;
- momentum transfer, including the acoustic `q²` factor;
- material density;
- temperature;
- absorption/emission and longitudinal/transverse branches.

It excludes the fitted deformation-potential values `A`, `B`, and the energy
threshold.

Consequently, changing only these deformation-potential values does not
invalidate the kernel. Changing the mesh, temperature, dispersion, carrier,
band selection, density, overlap model or BZ-domain convention does.

The CSV is sparse: states with zero kernels are omitted. A missing row therefore
means a zero kernel, not missing data.

Newly generated kernel CSVs also receive a `<kernel>.meta.yaml` sidecar. It
records the temperature, mesh size, carrier, band counts, energy window,
BZ-domain convention, Wigner–Seitz radius, and phonon dispersion. The fitting
script validates this metadata before starting. Older kernels without a
sidecar remain usable, but their compatibility cannot be checked automatically
and the fitter prints a warning.

### Energy window

`-E` selects the carrier-energy range for which kernels are calculated.

Use a window large enough for both:

- the measurements or reference data being fitted;
- the energies later reached by FBMC.

A small window is useful for quick validation, but a production FBMC run can
leave that range under a strong field.

### Temperature

Temperature is part of the kernel because phonon occupation depends on it.
A 300 K kernel must not be reused for a 200 K or 400 K rate calculation.

For a temperature sweep, generate one kernel per temperature.

## Producing rates from a kernel

The deformation-potential profile is stored under:

```text
data/materials/<material>/electron_phonon/<parameter-set>.yaml
```

For example:

```text
data/materials/Si/electron_phonon/remi-2026.yaml
```

For fitting runs or temporary profiles, an arbitrary YAML file can instead be
selected directly:

```bash
--phonon-params-file /full/path/to/trial-parameters.yaml
```

`--phonon-params-file` and `--phonon-params` are mutually exclusive. The
external file must still declare the selected material (for example
`material: Si`), `model: electron_phonon`, and a `parameter_set` label. Both
`elph.epm` and `fbmc.epm` support this option.

After editing `A`, `B`, or `energy-threshold`, reconstruct the rates with:

```bash
./build/apps/elph.epm \
  -f bz.msh \
  -m Si \
  --carrier electron \
  --phonon-params remi-2026 \
  -c 2 -v 0 \
  -T 300 \
  -E 0.3 \
  --kernel-file electron_kernels_300K.csv \
  --export-rates \
  --rates-out electron_rates.csv \
  --rates-only
```

This operation loads the mesh and kernel, applies the current YAML deformation
potential and writes the ordinary eight-channel rate table. It does not repeat
the final-state DOS calculation.

`--rates-only` exits after writing the requested files. It is useful for a
fitting loop where mobility and plotting are not needed.

Without `--rates-only`, export options are non-terminal. For example, this
command writes the reconstructed rate CSV and then continues to the electron
diagnostic, MRTA mobility and `--plot` workflow:

```bash
./build/apps/elph.epm \
  -f bz.msh \
  -m Si \
  --carrier electron \
  --phonon-params remi-2026 \
  -c 2 -v 0 \
  -T 300 \
  -E 0.3 \
  --kernel-file electron_kernels_300K.csv \
  --export-rates \
  --rates-out electron_rates.csv \
  --plot
```

The rate applied to each channel is:

```text
rate(Ei) = kernel(Ei) × [A + B min(Ei, energy-threshold)]
```

The acoustic kernel already contains `q²`. The optical kernel does not.

This command is the natural inner operation for fitting experiments: modify or
generate a parameter profile, reconstruct rates, evaluate the objective, and
repeat.

### Fitting low-field mobility to Arora

`python/FBMC/fit_elph_arora.py` performs the deterministic first fitting stage.
It reuses one electron kernel per temperature, generates temporary parameter
profiles, runs the MRTA mobility calculation, and compares it with the
zero-doping Arora mobility from `data/materials/Si/admc/arora-canali.yaml`.

For example:

```bash
python3 python/FBMC/fit_elph_arora.py \
  --mesh bz.msh \
  --kernel 200=electron_kernels_200K.csv \
  --kernel 250=electron_kernels_250K.csv \
  --kernel 300=electron_kernels_300K.csv \
  --kernel 350=electron_kernels_350K.csv \
  --kernel 400=electron_kernels_400K.csv \
  --ncbands 2 \
  --nvbands 4 \
  --output-dir fit_elph_arora
```

By default, only the zero-energy acoustic and optical strengths are fitted.
The high-energy slopes and energy threshold remain at their `remi-2026`
values because low-field mobility does not constrain them reliably. Use
`--evaluate-only` to evaluate the starting profile without optimization.

The fitting directory contains `best.yaml`, a per-evaluation `history.csv`,
the complete `elph.epm` logs, and `run_manifest.json` with SHA-256 hashes of
the mesh, profiles, and kernels. High-field FBMC data should be used in a
separate second stage to fit the energy dependence and threshold.

The fitter loads four valence bands by default because `elph.epm` solves the
intrinsic Fermi level before evaluating MRTA mobility. Kernel generation,
fitting, and manual validation should use the same band counts. The kernel
energy window may be larger than the MRTA fitting window; it only needs to
cover it. A single large-window kernel can therefore be reused for low-field
MRTA at that temperature, although smaller kernels are cheaper to generate.

### Fitting the high-field Canali curve

After the low-field fit, generate a 300 K kernel with an energy window large
enough for the high-field trajectories. For example:

```bash
./build/apps/elph.epm \
  -f bz.msh -m Si \
  --carrier electron \
  --phonon-params remi-2026 \
  -c 2 -v 4 \
  -T 300 -E 2.5 -j 32 -w \
  --export-kernels \
  --kernels-out electron_kernels_300K_2p5eV.csv \
  --rates-only
```

Then use the low-field result as the starting profile:

```bash
python3 python/FBMC/fit_elph_canali.py \
  --mesh bz.msh \
  --kernel electron_kernels_300K_2p5eV.csv \
  --base-params fit_elph_arora/best.yaml \
  --temperature 300 \
  --fields 10000,30000,100000,300000 \
  --max-energy 2.5 \
  --npart 300 \
  --time 20e-12 \
  --warmup 0.3 \
  --nthreads 32 \
  --seed 1234 \
  --output-dir fit_elph_canali
```

For every trial, the high-field fitter reconstructs the rate CSV, checks the
MRTA mobility, runs fixed-seed FBMC simulations, and compares the projected
drift speed with the zero-doping Canali curve. By default it fits the acoustic
and optical strengths at the energy threshold plus the threshold itself,
while preserving the low-energy `A` values obtained in the Arora stage.

Start with `--evaluate-only` to estimate runtime and inspect the baseline
curve. The main outputs are:

```text
best.yaml
best_rates.csv
history.csv
run_manifest.json
runs/eval_*/
```

The same random seed is reused for every parameter trial, providing common
Monte Carlo random numbers and reducing noise in parameter comparisons.
Final validation should use more particles, longer trajectories, multiple
seeds, and fields not used by the fit.

The two energy requirements are intentionally different:

- low-field MRTA only needs coverage of the thermally occupied states;
- high-field FBMC needs coverage of every carrier energy reached during the
  trajectory.

The high-field fitter rejects a requested FBMC `--max-energy` larger than the
kernel metadata window. It does not require the low-field and high-field
kernels to have the same window.

### Why the same phonon profile still matters

The profile selected by `--phonon-params` contains both phonon dispersion and
deformation-potential parameters.

When reusing a kernel, its dispersion section must remain the same as when the
kernel was generated. Only the deformation-potential section may safely vary.
Using a profile with a different dispersion can produce a numerically valid
rate file that is physically inconsistent with the kernel.

## Direct rate generation

Rates can still be computed directly without keeping a kernel:

```bash
./build/apps/elph.epm \
  -f bz.msh \
  -m Si \
  --carrier electron \
  --phonon-params remi-2026 \
  -c 2 -v 0 \
  -T 300 \
  -E 0.3 \
  -j 16 \
  --export-rates \
  --rates-out electron_rates.csv \
  --rates-only
```

This repeats the expensive calculation. It remains useful as a reference check
for the kernel workflow.

Kernels and rates can also be written during one traversal:

```bash
./build/apps/elph.epm \
  -f bz.msh \
  -m Si \
  --carrier electron \
  --phonon-params remi-2026 \
  -c 2 -v 0 \
  -T 300 \
  -E 0.3 \
  -j 16 \
  --export-kernels \
  --kernels-out electron_kernels_300K.csv \
  --export-rates \
  --rates-out electron_rates.csv \
  --rates-only
```

This is useful when validating a new mesh or model. Directly generated and
kernel-reconstructed rate files should agree.

## Running electron FBMC

FBMC consumes the ordinary rate CSV, not the kernel CSV:

```bash
./build/apps/fbmc.epm \
  -f bz.msh \
  -p electron_rates.csv \
  -m Si \
  --carrier electron \
  --phonon-params remi-2026 \
  -c 2 -v 0 \
  -N 1000 \
  -j 16 \
  -T 300 \
  -e 0.3 \
  -t 50e-12 \
  --warmup 0.2 \
  --Ex 1.0e4 \
  --seed 1234 \
  -d fbmc_electron
```

The important runtime choices are:

- `-N`: number of independent particles;
- `-t`: simulated time in seconds;
- `--warmup`: initial fraction discarded from steady-state averages;
- `--Ex`, `--Ey`, `--Ez`: electric field components in V/cm;
- `-e`: maximum modeled carrier energy in eV;
- `--seed`: reproducible random streams;
- `-j`: number of OpenMP threads.

The same carrier, bands, temperature and phonon profile used to create the rate
file should be used in FBMC.

FBMC still loads the phonon profile because it needs the dispersion and overlap
model when selecting a final state after a scattering channel has been chosen.
The deformation-potential factor does not affect the conditional final-state
choice: for a fixed initial state and channel, it is common to all candidates
and cancels during normalization.

### Maximum energy and null collisions

FBMC uses `-e` when constructing its maximum scattering-rate bound. It should
not be smaller than the energy range represented by the rate data or the
energies expected under the applied field.

If FBMC reports a self-scattering bound violation:

1. check that the rate file covers the required energy range;
2. increase the kernel/rate energy window if necessary;
3. increase FBMC `-e`;
4. use a slightly larger `--gamma-safety` if the bound only needs numerical
   margin.

Increasing `--gamma-safety` is not a substitute for missing high-energy rate
data.

### Warmup and reproducibility

`--warmup 0.2` ignores the first 20% of each trajectory when accumulating
steady-state observables. This reduces sensitivity to the thermal initial
state.

Set `--seed` when comparing parameter sets. Without it, FBMC generates a new
seed for every run, which adds Monte Carlo variation to the comparison.

## Running hole FBMC

Generate a hole rate file from the corresponding hole kernel, then run:

```bash
./build/apps/fbmc.epm \
  -f bz.msh \
  -p hole_rates.csv \
  -m Si \
  --carrier hole \
  --phonon-params remi-2026 \
  -c 0 -v 2 \
  -N 1000 \
  -j 16 \
  -T 300 \
  -e 0.3 \
  -t 50e-12 \
  --warmup 0.2 \
  --Ex 1.0e4 \
  --seed 1234 \
  -d fbmc_hole
```

Hole impact ionization is not currently implemented.

## Electron FBMC with impact ionization

Impact ionization is optional and configured independently from
electron–phonon scattering:

```bash
./build/apps/fbmc.epm \
  -f bz.msh \
  -p electron_rates.csv \
  -m Si \
  --carrier electron \
  --phonon-params remi-2026 \
  --enable-impact-ionization \
  --impact-ionization-params keldysh \
  -c 2 -v 0 \
  -N 1000 \
  -j 16 \
  -T 300 \
  -e 2.5 \
  -t 50e-12 \
  --Ex 5.0e5 \
  --seed 1234 \
  -d fbmc_impact
```

Impact-ionization profiles live under:

```text
data/materials/<material>/impact_ionization/
```

Bulk FBMC records the event and reduces the incident electron energy. It does
not currently add the generated electron–hole pair to the simulated ensemble.

## Useful optional outputs

`fbmc.epm --export-history` writes per-particle histories. This is valuable for
debugging individual trajectories but can create many large CSV files, so it is
disabled by default.

`fbmc.epm --test-elph` runs an electron–phonon diagnostic before transport.
It is useful during model development but adds work and is unnecessary for
routine sweeps.

The FBMC output directory also contains `run_info.txt`, which records the mesh,
profiles, field, temperature, random seed and impact-ionization configuration.
Keep this file with simulation results; it is the simplest runtime provenance
record.

For a nonzero electric field, each bulk FBMC run also prints and exports a
rough single-run mobility:

```text
mobility = |mean velocity projected along the field| / |field|
```

The `observables.csv` columns include the projected drift velocity, signed
mobility in m²/(V·s), and positive mobility magnitude in cm²/(V·s). This is a
quick diagnostic only; use `python/FBMC/run_field_sweep.py` with several
positive and negative low fields for a reliable zero-intercept mobility fit.

## Consistency checklist

Before trusting a reconstructed rate file or FBMC result, verify:

- the same BZ mesh is used;
- the same carrier and band counts are used;
- kernel and rates use the same temperature;
- the phonon dispersion has not changed since kernel generation;
- the same BZ-domain convention is used;
- the rate energy window covers the FBMC trajectory energies;
- FBMC uses the profile whose deformation potential generated the rate CSV;
- comparisons use a fixed random seed and adequate particle statistics.

The kernel CSV currently carries vertex indices and energies but not a complete
fingerprint of every physical setting. Compatibility therefore remains a
researcher-controlled runtime responsibility.

## Suggested fitting workflow

A practical fitting loop is:

1. Generate one kernel for each carrier and temperature of interest.
2. Keep the mesh and dispersion fixed.
3. Update `A`, `B`, and thresholds in a dedicated YAML profile.
4. Reconstruct a rate CSV from each kernel.
5. Compare rates or FBMC observables against the reference data.
6. Repeat with the same FBMC seed while tuning parameters.
7. Confirm the final profile with larger particle counts and several seeds.

The kernel makes steps 3–5 cheap. It does not remove Monte Carlo noise from
FBMC, so rate-level fitting is generally the faster first stage, followed by
transport-level validation.
