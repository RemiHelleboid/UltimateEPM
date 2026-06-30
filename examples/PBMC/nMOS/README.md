# PBMC nMOS examples

This directory contains 2D self-consistent PBMC examples for an nMOS device.
The mesh must provide the contact regions `source`, `drain`, `body`, and `gate`.
Input paths such as `nmos.msh` are resolved relative to the YAML file.

## Warnigs

The simulations are rather long to run ~1h each. 
You can plot the result in the notebook during the simulation, without waiting for them to finish.
ANd you can run them in parrallel, they are independant of each others.



## Cases

### `config_nmos_close.yaml`

Baseline nMOS run with:

- `Vsource = 0 V`
- `Vdrain = 0.1 V`
- `Vbody = 0 V`
- `Vgate = 0 V`
- Poisson mixing disabled
- `transport.max_energy_eV = 0.4`
- output directory: `nMOS_close`

This is the nominal closed/off-gate case. It is useful to check leakage-like
currents, contact injection, and whether the channel remains depleted.

Run:

```bash
./build/apps/device_pbmc.epm \
  --config examples/PBMC/nMOS/config_nmos_close.yaml
```

### `config_nmos_open.yaml`

Template for an open/on-gate nMOS run.

At the moment this file is identical to `config_nmos_close.yaml` except for the
output directory (`nMOS_open`). To make it an actual open-channel case, set the
gate voltage either in the YAML:

```yaml
contacts:
  voltages_V:
    gate: 1.0
```

or from the command line:

```bash
./build/apps/device_pbmc.epm \
  --config examples/PBMC/nMOS/config_nmos_open.yaml \
  --set contacts.voltages_V.gate=1.0
```

Use this case to compare the channel current with the gate biased above the
threshold-like condition of the simulated structure.

### `config_nmos_close_relaxed.yaml`

Relaxed off-gate case. It starts from the same electrical bias as
`config_nmos_close.yaml`, but changes numerical stabilization settings:

- output directory: `nMOS_close_relaxed`
- Poisson mixing enabled
- `poisson_mixing.old_solution_fraction = 0.8`
Run:

```bash
./build/apps/device_pbmc.epm \
  --config examples/PBMC/nMOS/config_nmos_close_relaxed.yaml
```
This eneabme more stable results, at the cost of loosing the temporal resolution. Plot of Y vs time are irrelevant physically with this config, it is only there to extract steady-state result (e.g. the drain current here).

## Outputs

Each run writes into the YAML `run.output_directory`.

Important files:

- `device_history.csv`: final full history exported at the end.
- `nmos_history.csv`: live history updated during the run.
- `simulation_manifest.txt`: configuration and output summary.
- `final_particle_state.csv`: final particle state, usable as a future initial state.
- `trajectory/mesh/mesh.pvd`: ParaView mesh time series.
- `trajectory/particles/particles.pvd`: ParaView particle time series.
- `trajectory/open_scene.py`: ParaView helper scene.

Open the ParaView scene with:

```bash
paraview --script nMOS_close/trajectory/open_scene.py
```

Change `nMOS_close` to the output directory of the case you ran.

## Current Columns

The history CSV uses PBMC-compatible column names. The most useful current
columns are:

- `ramo_current`: total Ramo current over the full device.
- `probe_ramo_current`: Ramo current restricted to the configured
  `current_probe` box.

For nMOS IV extraction, `probe_ramo_current` is usually preferable because the
box can be placed around the channel region.

The current probe is configured in each YAML:

```yaml
current_probe:
  enabled: true
  x_min_um: 0.1
  x_max_um: 0.25
  y_min_um: 0.12
  y_max_um: 0.20
```

This means that we only extract current from particle in this box, it helps reducing noise and background current due to the contact acceses. Because there is muuuuch more particle in the acceses so they have quite a huge impact.

## Gate IV Sweep

To sweep the gate voltage and extract an IV curve from `probe_ramo_current`:

```bash
python3 python/PBMC/amc_IV_runner.py \
  --exe ./build/apps/device_pbmc.epm \
  --config examples/PBMC/nMOS/config_nmos_open.yaml \
  --outdir IV_gate_sweep \
  --vmin 0.0 \
  --vmax 1.0 \
  --vstep 0.1 \
  --swept-contact gate \
  --sweep-voltage swept \
  --voltage-axis swept \
  --current-column probe_ramo_current \
  --jobs 1
```

This keeps the drain/source/body voltages from the YAML and only varies
`contacts.voltages_V.gate`.

For a drain-bias sweep at fixed gate voltage:

```bash
python3 python/PBMC/amc_IV_runner.py \
  --exe ./build/apps/device_pbmc.epm \
  --config examples/PBMC/nMOS/config_nmos_open.yaml \
  --outdir IV_drain_sweep \
  --vmin 0.0 \
  --vmax 0.5 \
  --vstep 0.05 \
  --swept-contact drain \
  --reference-contact source \
  --sweep-voltage bias \
  --voltage-axis bias \
  --current-column probe_ramo_current \
  --set contacts.voltages_V.gate=1.0 \
  --jobs 1
```

