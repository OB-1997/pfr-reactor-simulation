# Dynamic Simulation of a Tubular Chemical Reactor

**Fortran 2018 finite-volume solver + Python experiment pipeline for a stiff system of 2,000 coupled differential equations**

[![CI](https://github.com/OB-1997/pfr-reactor-simulation/actions/workflows/ci.yml/badge.svg)](https://github.com/OB-1997/pfr-reactor-simulation/actions/workflows/ci.yml)
![Fortran](https://img.shields.io/badge/Fortran-2018-734f96?logo=fortran)
![Python](https://img.shields.io/badge/Python-3.10+-3776ab?logo=python&logoColor=white)
![License: MIT](https://img.shields.io/badge/License-MIT-green)

Bachelor thesis project — *Dynamic Modelling of Tubular Reactor for Conversion of Pyrolytic Gas*, University of Chemistry and Technology, Prague (2026). 📄 [Read the full thesis (PDF)](thesis/main.pdf)

---

## What this is

A chemical reactor is, mathematically, a **partial differential equation**: chemical concentrations change along the reactor's length (flow, diffusion) and over time (reactions). This project builds a simulator for one such reactor — a 20-metre tube that converts plastic-waste pyrolysis gas into hydrogen-rich synthesis gas — entirely from scratch: no simulation packages, just a numerical method implemented directly.

The pipeline from physics to result:

1. **Model** — 10 chemical species, 8 reactions with exponential (Arrhenius) temperature dependence, and 9 side-injection points along the tube. Governing equations: 1-D convection–diffusion–reaction PDEs.
2. **Discretise** — the finite volume method (first-order upwind convection, central diffusion) turns the PDE into a system of **10 species × 200 grid cells = 2,000 coupled ordinary differential equations**. The system is *stiff*: reaction time-scales span many orders of magnitude, which is exactly what makes it numerically interesting.
3. **Integrate** — three solvers, implemented and benchmarked head-to-head:
   - **forward Euler** (explicit, simple, stability-limited),
   - **backward Euler** with a Newton iteration solving a banded Jacobian system each step (implicit, unconditionally stable),
   - **DLSODE** from ODEPACK (adaptive variable-order BDF — the production-grade reference).
4. **Verify & explore** — convergence studies against a closed-form analytical solution, wall-clock benchmarking, and parametric sweeps over temperature, feed velocity, and injection flow rate.

Simulated steady-state concentration profiles along the reactor (the steps are the side-injection points):

![Species concentration profiles along the reactor](thesis/figures/profiles_base.png)

## Highlights

**Numerical methods**
- Finite-volume discretisation of a convection–diffusion–reaction PDE (method of lines)
- Explicit vs implicit time integration on a genuinely stiff system, including a hand-written Newton solver with banded LU factorisation
- **Verification against an analytical solution**: a companion project swaps the 8-reaction chemistry for a single reaction with a known closed-form solution, isolating pure discretisation error. All three solvers reproduce the theoretical first-order convergence rate until they hit the accuracy floor:

![Grid convergence: error falls as 1/N for all three solvers](thesis/figures/test2_error_vs_ncells.png)

**Software engineering**
- ~1,700 lines of modular Fortran 2018 (parameters / model / solvers / I/O / CLI as separate modules), compiled with `-Wall -Wextra -fcheck=all`
- Runtime configuration via namelist files — every experiment is a config file, not a code edit
- **Full reproducibility**: each thesis figure maps to one Python orchestrator (runs the simulations) + one plot script; every figure in the results chapter can be regenerated with two commands
- Programmatic validation suite (`validate.py`): mass conservation, steady-state detection, and error bounds against two independent analytical references
- CI: both solvers build and run a smoke simulation on every push

## Repository map

| Path | What it is |
|---|---|
| [`code/isothermal_model/`](code/isothermal_model/) | The main simulator: Fortran solver + Python experiment pipeline. Its README is a complete user manual. |
| [`code/simple_pfr_analytical_test/`](code/simple_pfr_analytical_test/) | Verification benchmark against a closed-form analytical solution (same solver code, simplified chemistry). |
| [`code/fvm_python_prototype/`](code/fvm_python_prototype/) | The finite-volume method in ~200 lines of NumPy — the accessible entry point to the numerics. |
| [`thesis/`](thesis/) | LaTeX source and the compiled [thesis PDF](thesis/main.pdf). |

## Quickstart

Requires `gfortran`, `make`, and Python 3 with NumPy + Matplotlib (per-platform install commands in the [model README](code/isothermal_model/README.md)).

```bash
cd code/isothermal_model
make                            # build the Fortran binary
python3 runs/base_case.py       # simulate the design operating point
python3 plots/base_case.py      # -> results/[base_case]/profiles.png
```

Or the 200-line Python version, no compiler needed:

```bash
python3 code/fvm_python_prototype/main.py
```

## Selected results

- **Solver choice is problem-dependent, and measurably so.** The adaptive BDF integrator (DLSODE) needs ~100× fewer time steps than forward Euler, yet forward Euler wins on wall-clock time at engineering-relevant grid resolutions — its steps are so cheap that stability limits don't hurt until the grid gets fine. The crossover is mapped empirically (thesis §4.3–4.5).
- **All solvers converge at the theoretical first-order rate** on the analytical benchmark, confirming the implementation is correct and the dominant error is spatial discretisation, not time integration (§4.4).
- **The reactor has a narrow optimal operating window**: parametric sweeps over temperature, feed velocity, and side-injection flow rate locate a single optimum in syngas yield (§4.2).

## Why this project (a note for the ML/DS-minded reader)

This is the same toolbox that underlies modern scientific machine learning: systems of ODEs and their integrators (neural ODEs), discretising PDEs on grids (physics-informed networks, simulation surrogates), stiffness and conditioning, convergence analysis, and the discipline of validating numerical code against known ground truth. The project demonstrates that toolbox end-to-end — model derivation, implementation in a compiled language, verification, benchmarking, and automated, reproducible experiments.

## Author

**Ivan Hromakov** — ivan.gromakov@gmail.com

Thesis supervised by doc. Ing. Alexandr Zubov, Ph.D., Department of Chemical Engineering, UCT Prague.

Code is MIT-licensed (see [LICENSE](LICENSE)); the bundled [ODEPACK](https://computing.llnl.gov/projects/odepack) solvers (`code/*/lib/`) are public-domain LLNL software.
