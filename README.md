# Transient solute transport in statistically reconstructed pore structures with combinatorial mesh calculus

MATLAB code for

> C. Liu, P. Boom, M. Azeem, K. Berbatov, M. Sedighi, A. Jivkov. *Predicting transient solute transport from statistically reconstructed pore structures with combinatorial mesh calculus: arrival, dispersion and the role of the inlet concentration history.*

Pores are assigned to the cells of a Laguerre tessellation, M, from measured pore-volume distributions. Steady flow and transient advection–diffusion are solved on its Forman subdivision, K, with combinatorial mesh calculus, for a step inlet and an inlet concentration history at two flow rates. The cell-complex functions in `src/` are from [CMC-AdvDiff](https://github.com/boompiet/CMC-AdvDiff) by Dr Pieter Boom.

## Usage

MATLAB R2020a or later (`prctile` needs the Statistics and Machine Learning Toolbox in older releases). Run from this folder:

| Script | Content | Output |
|---|---|---|
| `run_transport.m` | Pore assignment, bulk permeability and transport for every tessellation in `data/` | `results/runs/*.mat`, `results/summary.csv` |
| `run_inlet_shape.m` | Exponential and sigmoidal inlet histories, SG1 and BS1 (Appendix A) | `results/inlet_shape/inlet_shape.csv` |

Saved runs are read, not recomputed. Settings (time step 0.1 s, end time, flow rates, diffusion coefficient, porosities, inlet rise times) are in `src/config_material.m`; each saved run records them in `R.opt`.

Expected output of `run_transport.m`:

| Realisation | Flow rate (µL/s) | k (m²) | v_p (µm/s) | v_arrival, step (µm/s) | v_arrival, history (µm/s) | T_f, history, first → last interval (s) |
|---|---|---|---|---|---|---|
| SG1 | 0.142 | 1.885 × 10⁻¹¹ | 31.90 | 45.90 | 37.01 | 182.5 → 212.1 |
| SG1 | 0.284 | | 63.80 | 91.08 | 74.21 | 92.2 → 107.4 |
| BS1 | 0.142 | 6.122 × 10⁻¹² | 45.94 | 68.43 | 54.77 | 187.2 → 213.3 |
| BS1 | 0.284 | | 91.88 | 135.00 | 109.22 | 102.1 → 115.0 |
| SGL1 | 0.142 | 1.890 × 10⁻¹¹ | 31.90 | 41.65 | 36.19 | 183.2 → 237.1 |
| SGL1 | 0.284 | | 63.80 | 83.31 | 72.72 | 92.5 → 120.6 |

## Data

| File | Content |
|---|---|
| `data/SG1.tess`, `data/BS1.tess` | Laguerre tessellations (Neper) of one realisation of the sintered glass (SG) and the Bentheimer sandstone (BS), 4 × 4 × 5 mm |
| `data/SGL1.tess` | Extended SG model, 4 × 4 × 10 mm |
| `data/<SG or BS>_voluminous_pores.csv` | Volumes of the voluminous pores (µm³) |
| `data/<SG or BS>_sheet_pores_cdf.csv` | Volume (µm³) and cumulative probability of the sheet-like pores |

The pore statistics come from the micro-CT images of Van Offenwert, Bultreys and Cnudde (2019), Digital Rocks Portal, https://doi.org/10.17612/5WW0-GH34. Further tessellations placed in `data/` as `SG<id>.tess` or `BS<id>.tess` are included by `run_transport.m`. Units: µm, s, Pa.

## License

BSD 2-Clause License, as for [CMC-AdvDiff](https://github.com/boompiet/CMC-AdvDiff). See `LICENSE`.
