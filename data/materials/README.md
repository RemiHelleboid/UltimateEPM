# Material repository

Each material owns one directory named by its chemical symbol:

```text
Si/
  material.yaml
  epm/
    local-cohen.yaml
    potz-vogl.yaml
  electron_phonon/
    kamakura.yaml
  impact_ionization/
    keldysh.yaml
  admc/
    arora-canali.yaml
  pbmc/
    default.yaml
```

`material.yaml` is the single source of truth for common identity and physical
values. Module profiles contain only model-specific parameters and identify
their material, model, and parameter set in metadata.

Code should access these files through `uepm::physics::material_repository`
instead of constructing paths directly.
