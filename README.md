# Transient solute transport in statistically reconstructed pore structures with combinatorial mesh calculus

MATLAB code for the calculations in

> C. Liu, P. Boom, M. Azeem, K. Berbatov, M. Sedighi, A. Jivkov. *Predicting transient solute transport from statistically reconstructed pore structures with combinatorial mesh calculus: arrival, dispersion and the role of the inlet concentration history.*

Voluminous and sheet-like pores are assigned to the cells of a Laguerre tessellation, M, from measured pore-volume distributions. Steady flow and transient advection–diffusion are solved on its Forman subdivision, K, with combinatorial mesh calculus (CMC). Conductances of the 1-cells of K follow from the pore composition of their cross-sections (diagonal Hodge star); advection uses exponential fitting (Scharfetter–Gummel) and time integration backward Euler. Two inlet conditions are solved together: an instantaneous step and an inlet concentration history with the measured 10–90% rise time.

The cell-complex and geometry functions in `src/` are from [CMC-AdvDiff](https://github.com/boompiet/CMC-AdvDiff) by Dr Pieter Boom.

## Requirements

MATLAB R2020a or later; no toolbox is needed in recent releases. In older releases, `prctile` requires the Statistics and Machine Learning Toolbox.

## Usage

Run the scripts from MATLAB in this folder (paths are relative to the scripts):

| Script | Content | Output |
|---|---|---|
| `run_transport.m` | Pore assignment, bulk permeability, step inlet and inlet history at both flow rates, for every realisation in `data/` | `results/runs/*.mat`, `results/summary.csv` |
| `run_inlet_shape.m` | Exponential and sigmoidal inlet histories (Appendix A) | `results/inlet_shape/inlet_shape.csv` |

Saved runs are read, not recomputed, so a script can be stopped and started again.

## Settings

All scripts take their settings from `src/config_material.m`:

| Setting | Value |
|---|---|
| Time step `dt` | 0.1 s |
| Simulated time `t_end` | 650 s (SG, BS), 1300 s (extended SG model) |
| Flow rates | `Q_low` = 1.42 × 10⁸ µm³/s (0.142 µL/s) and 2 `Q_low` |
| Molecular diffusion coefficient `Dm` | 1770 µm²/s |
| Dynamic viscosity `mu` | 10⁻³ Pa s |
| Monitoring planes `z_fraction`, slab half width `slab_halfwidth` | z/L = 0.1, 0.2, …, 0.9; 0.02 L |
| Per material: porosity `vf`, inlet rise time `Tf_in`, measured values `v_meas`, `dTf_meas` | SG, BS, SGL |

The end time of the inlet-shape study (800 s, to allow for the delay of the sigmoidal history) is set in `run_inlet_shape.m`. Every saved run records the settings it used in `R.opt` (`Q_in`, `Tf_in`, `dt`, `t_end`, `z_fraction`, `slab_halfwidth`, `art_scale`, `history`).

Expected output of `run_transport.m` for the two realisations included:

| Realisation | Flow rate (µL/s) | k (m²) | v_p (µm/s) | v_arrival, step (µm/s) | v_arrival, history (µm/s) | T_f, history, first → last interval (s) | T_f, step, outlet (s) | D_L,eq / D_m |
|---|---|---|---|---|---|---|---|---|
| SG1 | 0.142 | 1.885 × 10⁻¹¹ | 31.90 | 45.90 | 37.01 | 182.5 → 212.1 | 96.4 | 2.78 |
| SG1 | 0.284 | | 63.80 | 91.08 | 74.21 | 92.2 → 107.4 | 48.6 | 5.67 |
| BS1 | 0.142 | 6.122 × 10⁻¹² | 45.94 | 68.43 | 54.77 | 187.2 → 213.3 | 82.1 | 6.21 |
| BS1 | 0.284 | | 91.88 | 135.00 | 109.22 | 102.1 → 115.0 | 41.6 | 12.82 |

The flow rates are those of the experiments scaled to the cross-section of the models (Section 3.3).

## Files

| File | Content |
|---|---|
| `src/config_material.m` | Material data (porosity, inlet rise times, measured values) and run settings |
| `src/build_model.m` | Complexes, pore assignment, conductances and pore storage of one realisation |
| `src/build_complexes.m` | Physical complex M from a Neper tessellation and its Forman subdivision K |
| `src/assign_pores.m` | Voluminous pores on 3-cells and sheet-like pores on 2-cells of M |
| `src/pore_conductances.m` | Hydraulic and diffusive conductances of the 1-cells of K (composition of the cross-sections) |
| `src/pore_storage.m` | Pore volume stored at the vertices of K |
| `src/steady_flow.m` | Steady flow with a unit pressure difference |
| `src/bulk_permeability.m` | Bulk permeability |
| `src/solve_transport.m` | Steady flow scaled to the flow rate and transient advection–diffusion |
| `src/analyse_transport.m` | Arrival velocity, filling duration and equivalent dispersion coefficient (Ogata–Banks) |
| `src/threshold_metrics.m` | Arrival velocity and change of filling duration over a chosen set of nodes (Appendix A) |
| `src/find_realisations.m` | Realisations with a tessellation in `data/` |
| `src/setup_complex.m`, `read_tess.m`, `vgeo.m`, `vcon.m`, `v2f.m`, `fcon.m`, `fgeo2026.m`, `build_BO.m` | Cell complexes, geometry and boundary operators, from CMC-AdvDiff |
| `data/SG1.tess`, `data/BS1.tess` | Laguerre tessellations of one realisation of the sintered glass (SG) and the Bentheimer sandstone (BS), generated with Neper |
| `data/<MAT>_voluminous_pores.csv` | Volumes of the voluminous pores (µm³) |
| `data/<MAT>_sheet_pores_cdf.csv` | Volume (µm³) and cumulative probability of the sheet-like pores |

The pore statistics were obtained from the micro-CT images of Van Offenwert, Bultreys and Cnudde (2019), *Saturated solute transport micro-CT dataset in sintered glass and Bentheimer sandstone*, Digital Rocks Portal, https://doi.org/10.17612/5WW0-GH34.

The paper uses 20 realisations of each material and an extended SG model (4 × 4 × 10 mm). Only SG1 and BS1 are included here; with further tessellations `data/SG<id>.tess`, `data/BS<id>.tess` or `data/SGL1.tess` in `data/`, `run_transport.m` includes them.

Lengths are in µm, times in s and pressures in Pa.

## License

BSD 2-Clause License, as for [CMC-AdvDiff](https://github.com/boompiet/CMC-AdvDiff). See `LICENSE`.
