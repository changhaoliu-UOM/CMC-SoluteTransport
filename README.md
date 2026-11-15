# Transient solute transport in a reconstructed pore structure with combinatorial mesh calculus

MATLAB code for the transport calculations in

> C. Liu, P. Boom, M. Azeem, K. Berbatov, M. Sedighi, A. Jivkov. *Predicting transient solute transport from statistically reconstructed pore structures with combinatorial mesh calculus: arrival, dispersion and the role of the inlet concentration history.*

The code is based on the combinatorial mesh calculus (CMC) implementation of advection–diffusion by Dr Pieter Boom, [CMC-AdvDiff](https://github.com/boompiet/CMC-AdvDiff).

The repository contains one reconstruction of the sintered glass sample (SG1). Voluminous and sheet-like pores are assigned to the cells of a Laguerre tessellation, M. Steady flow and transient advection–diffusion are then solved on its Forman subdivision, K. Two inlet concentration histories with the same 10–90% rise time are considered, an exponential and a sigmoidal (erf) one, each at two flow rates.

## Requirements

MATLAB R2019a or later. No toolboxes are needed.

## Usage

```matlab
run_SG1        % builds the model and solves both flow rates -> results/SG1_low.mat, SG1_high.mat
analyse_SG1    % arrival velocities and filling durations    -> results/SG1_summary.csv, SG1_filling_duration.csv
```

## Files

| File | Content |
|---|---|
| `run_SG1.m` | Settings and calculation |
| `analyse_SG1.m` | Arrival velocity and filling duration |
| `src/build_complexes.m` | Physical complex M and computational complex K |
| `src/assign_pores.m` | Pore assignment and hydraulic and diffusion coefficients |
| `src/pore_storage.m` | Pore storage volume of the nodes of K |
| `src/solve_transport.m` | Steady flow and transient transport |
| `src/setup_complex.m`, `read_tess.m`, `vgeo.m`, `vcon.m`, `v2f.m`, `fcon.m`, `fgeo2026.m`, `build_BO.m` | Cell complexes, geometry and boundary operators, from [CMC-AdvDiff](https://github.com/boompiet/CMC-AdvDiff) |
| `data/SG1.tess` | Laguerre tessellation generated with Neper |
| `data/SG_voluminous_pores.csv` | Volumes of the voluminous pores (µm³) |
| `data/SG_sheet_pores_cdf.csv` | Volume (µm³) and cumulative probability of the sheet-like pores |

Lengths are in µm and times in s.

## License

BSD 2-Clause License, as for [CMC-AdvDiff](https://github.com/boompiet/CMC-AdvDiff). See `LICENSE`.
