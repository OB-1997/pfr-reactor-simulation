# Simple A → B PFR — analytical accuracy benchmark

Sibling project to [`../isothermal_model/`](../isothermal_model/). Same FVM
machinery (first-order upwind convection, central diffusion,
method-of-lines discretisation) and **byte-identical** `solvers.f90`,
`output.f90`, `main.f90`, but with a single first-order reaction `A → B`
in place of the 8-reaction kinetics of the real model.

A closed-form analytical steady-state solution exists for the simplified
PDE, so the spatial-discretisation error of each solver can be measured
directly against the analytical reference. The figure produced by this
project — `plots/error_vs_ncells.png` — is §4.4 of the bachelor thesis
*"Dynamic Modelling of Tubular Reactor for Conversion of Pyrolytic Gas"*
(VŠCHT Prague, 2026).

This README is the user manual: it tells you how to build the code, how to
reproduce the §4.4 figure, what the analytical reference is, and how to
extend the benchmark to other test problems.

---

## 1 · Why a separate project

The 8-reaction model has no analytical solution, so the only way to
measure its spatial-discretisation error in isolation is to substitute the
chemistry with something that *does* have one. Replacing the eight
Arrhenius reactions with a single linear `A → B` step gives a **direct**
error measurement while leaving the discretisation, boundary treatment,
and time-integration code untouched. The accuracy ranking observed here
therefore transfers — qualitatively — to the real model: any solver that
reaches the spatial-discretisation floor on this test does so on the
real model too, because the time integrators see the same FVM
right-hand side.

---

## 2 · What this code produces

The §4.4 figure of the thesis: `plots/error_vs_ncells.png`. It plots the
steady-state max-norm relative error of the simulated `C_A` profile
against the analytical reference, for each of the three solvers (DLSODE,
forward Euler, backward Euler) across six grid resolutions
`N_cells ∈ {50, 100, 200, 400, 800, 1000}`. The accompanying
`plots/error_vs_ncells.csv` contains the underlying numbers.

---

## 3 · Toolchain

Identical to `../isothermal_model/`: gfortran 13+ (Fortran 2018), GNU make,
Python 3 with NumPy + Matplotlib. See that project's README §2 for the
per-platform install commands. If `make` works there, it works here.

---

## 4 · Build

From the project root (`code/simple_pfr_analytical_test/`):

```bash
make
```

This produces the binary `./pfr_simple` (or `pfr_simple.exe` on Windows).
The orchestrator scripts rebuild automatically when sources change, so
you rarely need to run `make` by hand. `make clean` wipes build artifacts.

---

## 5 · Reproduce the §4.4 figure

One command:

```bash
python3 sweep.py
```

`sweep.py` builds the binary if needed, runs `./pfr_simple` 18 times
(3 solvers × 6 grids), computes the max-norm relative error of each run
against the analytical reference, writes `plots/error_vs_ncells.csv`, and
plots `plots/error_vs_ncells.png`.

If you only want to re-plot from the existing CSV (e.g. after tweaking the
plotting style), skip the simulations:

```bash
python3 sweep.py --plot-only
```

---

## 6 · Running the Fortran binary directly

```bash
./pfr_simple --solver=<name> [--config=<path>]
```

| `--solver=…`       | Algorithm                                              |
|--------------------|--------------------------------------------------------|
| `odepack`          | DLSODE (BDF + banded FD Jacobian, `mf = 25`)           |
| `explicit_euler`   | Forward Euler, fixed Δt = 1 ms                         |
| `implicit_euler`   | Backward Euler with banded-Jacobian Newton, Δt = 10 ms |

Each invocation writes an auto-incremented folder
`output/<solver>_<N>/` with the same `parameters.txt` / `output.txt` /
`performance.txt` triple as the real-model project. Output format:
`output.txt` is `(t, C_A[1], C_B[1], C_A[2], C_B[2], …, C_A[N], C_B[N])`
on each row, species varying fastest.

A pre-shipped namelist `runs/standardized.nml` reproduces the thesis test
conditions explicitly:

```bash
./pfr_simple --solver=odepack --config=runs/standardized.nml
```

The `runs/sweep_N*.nml` files are generated automatically by `sweep.py`;
they are kept in the repo as a reference of what was actually run.

---

## 7 · Programmatic sanity check

`validate.py` runs five checks against an existing pair of output folders:

```bash
python3 validate.py output/odepack_1 output/odepack_6
```

The first argument is any smoke run; the second is the reference
high-resolution DLSODE run used for the analytical comparison. Checks:

1. Build artefacts and output files are present
2. **Mass conservation** at steady state: `max |C_A + C_B − C_A,in| < 1e-4`
3. **Steady state reached**: last vs second-to-last sample differ by less than `1e-4`
4. **Analytical match** vs the Dirichlet/Neumann solution at N = 1000: max-norm relative error ≤ `5e-3` (typically ~`1e-3`)
5. **Cross-check** vs the Wehner–Wilhelm closed-closed Danckwerts solution
The shipped `validation_log.txt` is a saved transcript of a passing run.

---

## 8 · Model

Two species, A and B. RHS for cell *j*:

```
dC_A/dt = -u (C_A,j  - C_A,j-1) / dx + D_ax (C_A,j+1 - 2 C_A,j + C_A,j-1) / dx² - k C_A,j
dC_B/dt = -u (C_B,j  - C_B,j-1) / dx + D_ax (C_B,j+1 - 2 C_B,j + C_B,j-1) / dx² + k C_A,j
```

Boundary conditions: Dirichlet ghost at the inlet (`C_A = C_A,in`,
`C_B = 0`), Neumann zero-gradient at the outlet. Initial condition: vacuum
(`C(:, 0) = 0`), matching the real-model convention.

**Parameters**

| symbol      | value      | unit       | source                                       |
|-------------|-----------:|------------|----------------------------------------------|
| L           | 20.0       | m          | inherited from real model (Table 3.1)        |
| u           | 0.516      | m/s        | inherited                                    |
| D_ax        | 1×10⁻²     | m²/s       | inherited (Levenspiel turbulent band)        |
| k           | 0.026      | 1/s        | sized for Da ≈ 1 (residence ≈ reaction time) |
| C_A,in      | 10.4       | mol/m³     | ≈ P/(RT) of real model (10.39)               |

**Derived**

| quantity               | value    |
|------------------------|---------:|
| Pe = u·L / D_ax        | 1032     |
| Da = k·L / u           | 1.008    |
| τ  = L / u             | 38.76 s  |

---

## 9 · Analytical reference

The steady-state ODE is `D·d²C/dx² − u·dC/dx − k·C = 0`. In dimensionless
form with `z = x/L`, `Pe = u·L/D`, `Da = k·L/u`, `β = √(1 + 4·Da/Pe)`, the
characteristic roots are `m± = Pe(1 ± β)/2`.

**Dirichlet/Neumann reference (used for the §4.4 figure).** With the BCs
the discrete FVM enforces — Dirichlet at the inlet, zero-gradient at the
outlet — the solution is

```
C(z) / C_A,in  =  X / (1 + X) · exp(m+ · z)  +  1 / (1 + X) · exp(m- · z)
                  with  X = −m- · exp(m- − m+) / m+
```

At our `Pe ≈ 10³` the factor `exp(m- − m+) = exp(−β·Pe)` underflows to
zero, so `X → 0` and the solution collapses to the clean single-exponential
form

```
C(z) / C_A,in  ≈  exp(m- · z)  ≈  exp(−Da · z)
```

This is what the FVM scheme converges to in the `dx → 0` limit — so the
gap between simulation and reference is *purely* spatial-discretisation
error, with no boundary-condition mismatch. This is the right reference
for the §4.4 figure.

**Wehner–Wilhelm sanity check.** The same ODE under closed-closed
Danckwerts BCs (Wehner & Wilhelm, *Chem. Eng. Sci.* 6 (1956) 89–93;
Levenspiel ch. 13; Fogler ch. 14) is the textbook reference. At our
Pe ≈ 10³ it differs from the Dirichlet/Neumann solution by O(1/Pe) ~ 10⁻³
across the whole reactor, comparable to the upwind discretisation error.
`validate.py` reports both reference errors so the gap is documented; only
the Dirichlet/Neumann column is plotted in §4.4.

Both formulae are implemented in `validate.py` (and re-imported by
`sweep.py`), with numerically stable branches for the small-`X` regime.

---

## 10 · Drift management with the real model

`solvers.f90`, `main.f90`, and `output.f90` are intended to be
**byte-identical** to `../isothermal_model/src/`. If a bug surfaces in any
of those files, fix it in the real model first and re-copy here — drift
would invalidate the benchmark, since the whole point is that the time
integrators see the same code on both projects.

Only `params.f90` and `model.f90` are project-specific. To make the
verbatim files compile against `N_species = 2`, the legacy species-index
symbols (`iCH4`, `iC9H10`, `iO2`, `iH2`, `iCO`, `iCO2`, `iH2O`) are
retained as aliases mapped onto species 1 (A) or 2 (B); the legacy
Arrhenius arrays (`A1..A8`, `E1..E8`) are zero-valued placeholders. The
`parameters.txt` echo is therefore misleading by design — the simulation
data is correct.

---

## 11 · Project layout

```
simple_pfr_analytical_test/
├── src/                 Fortran 2018 sources (5 modules)
│   ├── params.f90       reduced-model parameters; legacy aliases for src-shared code
│   ├── model.f90        single-reaction RHS, FVM stencil
│   ├── solvers.f90      forward / backward Euler + DLSODE wrapper (verbatim copy)
│   ├── output.f90       per-run folder + writers (verbatim copy)
│   └── main.f90         CLI driver (verbatim copy)
├── lib/                 ODEPACK F77 sources (DLSODE + LINPACK)
├── runs/                shipped namelists; sweep_N*.nml are regenerated by sweep.py
├── output/              per-run Fortran outputs (auto-incremented)
├── plots/               error_vs_ncells.{png,csv} + profiles_overlay.png
├── build/               gfortran .o + .mod (regenerated by make)
├── Makefile             build / clean targets
├── sweep.py             reproduces the §4.4 figure end-to-end
├── validate.py          programmatic sanity checks against analytical references
├── validation_log.txt   saved transcript of a passing validate.py run
└── README.md            this file
```

---
