#!/usr/bin/env bash
set -euo pipefail

# Electron-phonon fitting workflow for silicon.
#
# Run from an empty working directory. Existing meshes and kernels are reused,
# so the script can be restarted after an interruption.

PROJECT_ROOT=${PROJECT_ROOT:-"$HOME/UltimateEPM"}
PYTHON_DIR="$PROJECT_ROOT/python/FBMC"
NTHREADS=${NTHREADS:-8}
NCBANDS=5
NVBANDS=4

# The low-field fit needs at least as many temperatures as fitted parameters.
# Five temperatures constrain the temperature dependence more reliably than
# the mathematical minimum of two.
LOW_FIELD_TEMPERATURES=(200 250 300 350 400)
LOW_FIELD_WINDOW_EV=1.0
HIGH_FIELD_WINDOW_EV=10.0

# Generate the Brillouin-zone mesh and compute its band structure.
if [[ ! -f bz.msh ]]; then
  bz_meshing.epm --level 0 --mesh 0.1 -o bz.msh --nogui
  BandsOnBZ \
    -f bz.msh \
    -m Si \
    -v "$NVBANDS" \
    -c "$NCBANDS" \
    -o bz.msh \
    -n 10 \
    -j "$NTHREADS" \
    -w \
    -d local-remi-2026
fi

generate_kernel() {
  local temperature=$1
  local energy_window=$2
  local output=$3

  if [[ -f "$output" && -f "${output}.meta.yaml" ]]; then
    echo "Reusing $output"
    return
  fi

  elph.epm \
    -f bz.msh \
    -m Si \
    --carrier electron \
    --phonon-params remi-2026 \
    -c "$NCBANDS" \
    -v "$NVBANDS" \
    -T "$temperature" \
    -E "$energy_window" \
    -j "$NTHREADS" \
    --export-kernels \
    --kernels-out "$output" \
    --rates-only \
    --skip-mesh-vtk \
    -w
}

# Keep the 300 K kernel large enough for the later high-field fit. The other
# temperatures only need the smaller low-field window.
for temperature in "${LOW_FIELD_TEMPERATURES[@]}"; do
  energy_window=$LOW_FIELD_WINDOW_EV
  if [[ "$temperature" == 300 ]]; then
    energy_window=$HIGH_FIELD_WINDOW_EV
  fi
  generate_kernel \
    "$temperature" \
    "$energy_window" \
    "electron_kernels_${temperature}K.csv"
done

# Fit the zero-energy acoustic and optical strengths to Arora mobility.
python3 "$PYTHON_DIR/fit_elph_arora.py" \
  --mesh bz.msh \
  --kernel 200=electron_kernels_200K.csv \
  --kernel 250=electron_kernels_250K.csv \
  --kernel 300=electron_kernels_300K.csv \
  --kernel 350=electron_kernels_350K.csv \
  --kernel 400=electron_kernels_400K.csv \
  --ncbands "$NCBANDS" \
  --nvbands "$NVBANDS" \
  --nthreads "$NTHREADS" \
  --energy-window "$LOW_FIELD_WINDOW_EV" \
  --output-dir fit_elph_arora

# Fit the high-energy strengths and threshold to the 300 K Canali curve.
# Start with --evaluate-only if you only want to validate the setup/runtime.
python3 "$PYTHON_DIR/fit_elph_canali.py" \
  --mesh bz.msh \
  --kernel electron_kernels_300K.csv \
  --base-params fit_elph_arora/best.yaml \
  --temperature 300 \
  --fields 10000,30000,100000,300000 \
  --ncbands "$NCBANDS" \
  --nvbands "$NVBANDS" \
  --nthreads "$NTHREADS" \
  --max-energy 8.0 \
  --npart 300 \
  --time 20e-12 \
  --warmup 0.3 \
  --seed 1234 \
  --output-dir fit_elph_canali

# Validate the fitted profile with a field sweep. Impact ionization is enabled
# by default; add --disable-impact-ionization for a phonon-only sweep.
python3 "$PYTHON_DIR/run_field_sweep.py" \
  --exe "$(command -v fbmc.epm)" \
  --mesh bz.msh \
  --phonon-params-file fit_elph_canali/best.yaml \
  --phonon-rates fit_elph_canali/best_rates.csv \
  --outdir fbmc_best_sweep \
  --fields=-1000,1000,10000,30000,100000,300000 \
  --ncbands "$NCBANDS" \
  --nvbands "$NVBANDS" \
  --nbthreads "$NTHREADS" \
  --npart 1000 \
  --time 20e-12 \
  --warmup 0.3 \
  --temperature 300 \
  --max-energy 8.0 \
  --gamma-safety 1.3 \
  --seed 1234
