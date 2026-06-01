# Git Diff Summary and Rationale

Date: 2026-06-01

This document summarizes the current unstaged git diff and explains the likely rationale for each change.

## High-level Theme

The changes migrate the wrapper from SciPy-bundled OpenBLAS symbol conventions to direct BLAS/LAPACK symbol names, and switch Meson linkage from `scipy_openblas` to `flexiblas`.

## File-by-file Changes

### asv.conf.json

- Changed benchmark branch from `tracker` to `main`.

Rationale:
- Aligns Airspeed Velocity benchmarks with the default branch so benchmark runs target the active development line.

### meson.build

- Replaced:
  - `dependency('scipy_openblas', method: 'pkg-config', required: true)`
  - `openblas_dep = declare_dependency(...)`
- With:
  - `dependency('flexiblas', method: 'pkg-config', required: true)`
  - `flexiblas_dep = declare_dependency(...)`

Rationale:
- Updates the project to link against FlexiBLAS instead of SciPy's packaged OpenBLAS.
- Makes the BLAS provider explicit at system level through pkg-config.

### openblas_wrap/__init__.py

- Removed eager import of `scipy_openblas32`.
- Changed `PREFIX` from `scipy_` to empty string.

Rationale:
- No longer depends on SciPy's symbol preloading path.
- Generated/runtime symbol lookup now expects unprefixed BLAS/LAPACK names (for example `dgemm` instead of `scipy_dgemm`).

### openblas_wrap/_distributor_init.py

- Removed `import scipy_openblas32`.

Rationale:
- Eliminates package-level preload hook that was specific to SciPy OpenBLAS distribution mechanics.
- Keeps initialization consistent with new external BLAS linking strategy.

### openblas_wrap/blas_lapack.pyf.src

- Updated all template prefixes:
  - `scipy_s/scipy_d/scipy_c/scipy_z` -> `s/d/c/z`
  - Similar removals in `prefix2`, `prefix2c`, `prefix3`, `prefix4`, `prefix6`.
- Renamed explicit function wrapper `scipy_ddot` -> `ddot` and adjusted:
  - return variable name
  - `intent(c)` name
  - `fortranname F_FUNC(...)` mapping

Rationale:
- Matches the generated f2py interface to canonical BLAS symbol names without SciPy-specific prefixes.
- Ensures C/Fortran binding declarations remain internally consistent after renaming.

### openblas_wrap/meson.build

- Updated `_flapack` extension dependencies from `[openblas_dep, fortranobject_dep]` to `[flexiblas_dep, fortranobject_dep]`.

Rationale:
- Completes build-system dependency rename so the extension links to FlexiBLAS consistently.

### pyproject.toml

- Removed `scipy_openblas32` from `build-system.requires`.
- Removed `scipy_openblas32` from `project.dependencies`.

Rationale:
- Drops Python package dependency that is no longer required when linking against system BLAS via Meson/pkg-config.
- Reduces wheel/runtime dependency surface.

## Additional Working Tree Note

- Untracked directory detected: `.asv/`.

Rationale:
- Likely local benchmark cache/artifacts; typically not committed unless intentionally versioned.

## Risk/Validation Checklist

- Confirm `flexiblas` and its pkg-config metadata are available in all target environments.
- Rebuild extension and verify symbol resolution for representative functions (`ddot`, `dgemm`, `dnrm2`).
- Run benchmark and import smoke tests to verify behavior parity after prefix migration.