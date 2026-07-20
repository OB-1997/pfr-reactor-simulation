# Isothermal Tubular Reformer — Fortran model
This README is the user manual: it tells you how to build the code, how to
reproduce every figure in chapter 4, and where to look if you want to extend
the model.

---

## 1 · What this code produces

Every figure in §4.2 / §4.3 / §4.5 of the thesis comes out of this folder.
Each thesis subsection is backed by one **experiment**, and each experiment
is a paired *orchestrator* + *plotter*:

| Thesis section                    | Experiment              | Orchestrator                          | Plotter                              | Output figure                                          |
|-----------------------------------|-------------------------|---------------------------------------|--------------------------------------|--------------------------------------------------------|
| §4.2.1 design operating point     | `base_case`             | `runs/base_case.py`                   | `plots/base_case.py`                 | `results/[base_case]/profiles.png`                     |
| §4.2.2 + §4.2.3 (T, u) sweep      | `parametric_sweep`      | `runs/parametric_sweep.py`            | `plots/parametric_sweep.py`          | `results/[parametric_sweep]/profiles_{T,u}_sweep.png`  |
| §4.2.4 side-injection F sweep     | `parametric_sweep_F`    | `runs/parametric_sweep_F.py`          | `plots/parametric_sweep_F.py`        | `results/[parametric_sweep_F]/profiles_F_sweep.png`    |
| §4.3 wall-clock vs grid           | `wall_time_data`        | `runs/wall_time_data.py`              | `plots/wall_time_data.py`            | `results/[wall_time_data]/walltime.png` (+ .csv)       |
| §4.5 explicit vs implicit Euler   | `euler_comparison`      | `runs/euler_comparison.py`            | `plots/euler_comparison.py`          | `results/[euler_comparison]/cpu_vs_n.png` (+ .csv)     |

The chapter-4 accuracy figure (§4.4) lives in the sibling project
[`../simple_pfr_analytical_test/`](../simple_pfr_analytical_test/).

---

## 2 · Toolchain

You need a Fortran 2018 compiler, GNU make, and Python 3 with NumPy +
Matplotlib. Tested with gfortran 13.3, Python 3.12, Ubuntu 24.04.

**Linux (Ubuntu/Debian):**

```bash
sudo apt install gfortran make python3 python3-numpy python3-matplotlib
```

**macOS (Homebrew):**

```bash
brew install gcc make python
pip3 install numpy matplotlib
```

**Windows (MSYS2):** install [MSYS2](https://www.msys2.org), open the
*MSYS2 MINGW64* shell, then

```bash
pacman -Syu
pacman -S --needed mingw-w64-x86_64-gcc-fortran mingw-w64-x86_64-make \
                   mingw-w64-x86_64-python mingw-w64-x86_64-python-numpy \
                   mingw-w64-x86_64-python-matplotlib
```

Use `mingw32-make` in place of `make`. The produced binary is
`prototype.exe`.

---

## 3 · Build

From the project root (`code/isothermal_model/`):

```bash
make
```

This produces the binary `./prototype` (or `prototype.exe` on Windows).
The orchestrator scripts in `runs/` will rebuild automatically when any
`src/*.f90` file is newer than the binary, so in practice you rarely need
to call `make` by hand. To wipe build artifacts:

```bash
make clean
```

Compiled-output directories under `output/` are preserved by `make clean`;
delete them manually if you want a fully fresh state.

---

## 4 · Reproduce a thesis figure

Each experiment is a single Python command that (a) ensures the binary is
up to date, (b) runs the Fortran simulation across the experiment's
parameter set, and (c) leaves results in `output/[<experiment>]/`. A second
Python command reads those results and writes the PNG into
`results/[<experiment>]/`.

### §4.2.1 — design operating point

```bash
python3 runs/base_case.py
python3 plots/base_case.py
```

One ODEPACK run at the canonical operating point (T = 1173.15 K, u = 0.516
m/s, N = 200). Produces `results/[base_case]/profiles.png` — the
two-panel figure with all ten species along the reactor.

### §4.2.2 + §4.2.3 — temperature and feed-velocity sweep

```bash
python3 runs/parametric_sweep.py        # 9 runs (3 T × 3 u)
python3 plots/parametric_sweep.py
```

Produces both `profiles_T_sweep.png` and `profiles_u_sweep.png` from the
same 3 × 3 grid.

### §4.2.4 — side-injection flow-rate sweep

```bash
python3 runs/parametric_sweep_F.py      # 3 runs (F_scale = 0.5, 1.0, 1.5)
python3 plots/parametric_sweep_F.py
```

Produces `profiles_F_sweep.png`.

### §4.3 — wall-clock vs grid resolution

```bash
python3 runs/wall_time_data.py          # 3 solvers × 6 grids = 18 runs
python3 plots/wall_time_data.py
```

Produces `walltime.png` and the underlying `walltime.csv`. This experiment
takes the longest — the DLSODE run at N = 800 alone is about one minute on
the reference workstation (see §4.1 of the thesis).

### §4.5 — explicit vs implicit Euler crossover

```bash
python3 runs/euler_comparison.py        # 4 explicit + 20 implicit = 24 runs
python3 plots/euler_comparison.py
```

Produces `cpu_vs_n.png` and `cpu_vs_n.csv`. The implicit-Euler step is
swept over Δt ∈ {5, 10, 25, 50, 100} ms across four grids.

### All figures in one shot

```bash
for exp in base_case parametric_sweep parametric_sweep_F wall_time_data euler_comparison; do
    python3 "runs/$exp.py" && python3 "plots/$exp.py"
done
```

---

## 5 · Running the Fortran binary directly

`./prototype` is a self-contained CLI. The orchestrators above call it for
you; this section is for ad-hoc runs.

```bash
./prototype --solver=<name> [--config=<path>] [--dt-step=<seconds>]
```

| `--solver=…`       | Algorithm                                                              |
|--------------------|------------------------------------------------------------------------|
| `odepack`          | DLSODE (variable-order BDF + banded FD Jacobian, `mf = 25`)            |
| `explicit_euler`   | Forward Euler, fixed Δt = 1 ms (hardcoded — CFL bound on this problem) |
| `implicit_euler`   | Backward Euler with banded-Jacobian Newton, default Δt = 10 ms         |

`--config=<path>` points at a Fortran NAMELIST file (see §6). `--dt-step=<seconds>`
overrides the implicit-Euler step size (used by `runs/euler_comparison.py`).

Each invocation writes an auto-incremented folder
`output/<solver>_<N>/` containing:

- `parameters.txt` — full echo of geometry, kinetics, inlet, injectors,
  time window, solver settings (open this first to verify the run setup)
- `output.txt` — concentration snapshots. Row format: `(t, C[1,1], C[2,1],
  …, C[N_species, N_cells])`. Species index varies fastest, cell index
  slowest.
- `performance.txt` — wall-clock time, CPU time, RHS evaluations, Jacobian
  evaluations, internal step count (NST), final BDF order

The orchestrators in `runs/` then move these `<solver>_<N>` folders into
the per-experiment buckets `output/[<experiment>]/`.

---

## 6 · Configuration via NAMELIST

Geometry, operating point, dispersion, grid resolution, side-injection
multiplier, and time-integration window are all runtime-configurable via a
Fortran NAMELIST file. Three groups are recognised — any of them can be
omitted, in which case its variables fall back to the compiled-in defaults
declared in `src/params.f90`.

```
&reactor
    L_reactor = 20.0          ! reactor length, m
    d_reactor = 0.60          ! inner diameter, m
    u_vel     = 0.516         ! superficial feed velocity, m/s
    T_react   = 1173.15       ! isothermal wall temperature, K
    P_react   = 101325.0      ! operating pressure, Pa
    D_ax      = 1.0e-2        ! axial dispersion coefficient, m²/s
    N_cells   = 200           ! finite-volume grid resolution
    F_scale   = 1.0           ! global side-injection mass-flow multiplier
/

&time_window
    t_start = 0.0             ! s
    t_end   = 50.0            ! s
    dt_out  = 0.5             ! output cadence, s
/

&solver
    dt_step_impl = 1.0e-2     ! implicit-Euler step in seconds (overridable from CLI too)
/
```

Per-experiment namelists live alongside the orchestrators in
`runs/[<experiment>]/`. The shipped files only override the variables that
distinguish each run (e.g. `T_react`, `u_vel`, `N_cells`, `F_scale`);
everything else falls through to the defaults.

**Compile-time vs runtime.** The kinetic model itself — species count,
stoichiometry matrix, Arrhenius pre-exponentials and activation energies,
the injector schedule (which species enters at which axial position) — is
**not** namelist-configurable. To change the reaction set you must edit
`src/params.f90` (kinetics, stoichiometry, injector schedule) and
`src/model.f90` (RHS and Jacobian-pattern matching the new stoichiometry),
then rebuild.

---

## 7 · Project layout

```
isothermal_model/
├── src/                 Fortran 2018 sources (5 modules)
│   ├── params.f90       reactor + kinetic parameters, namelist reader
│   ├── model.f90        8-reaction kinetics, FVM RHS, injector setup
│   ├── solvers.f90      forward Euler, backward Euler (Newton), DLSODE wrapper
│   ├── output.f90       per-run output folders + parameters/output/performance writers
│   └── main.f90         CLI driver
├── lib/                 ODEPACK F77 sources (DLSODE + LINPACK), legacy-flagged in Makefile
├── runs/                Python orchestrators + per-experiment namelist buckets
│   ├── base_case.py
│   ├── parametric_sweep.py
│   ├── parametric_sweep_F.py
│   ├── wall_time_data.py
│   ├── euler_comparison.py
│   ├── [base_case]/canonical.nml
│   ├── [parametric_sweep]/T{1073,1173,1273}_u{400,516,700}.nml      (9 files)
│   ├── [parametric_sweep_F]/F{050,100,150}.nml
│   ├── [wall_time_data]/sweep_N{50,100,200,400,800,1000}.nml
│   └── [euler_comparison]/N{100,200,400,800}.nml
├── plots/               Python plot scripts (one per experiment) + _common.py
├── output/              per-run Fortran outputs; bucketed by experiment after orchestration
├── results/             figures + CSVs that go into the thesis, one folder per experiment
├── build/               gfortran .o + .mod (regenerated by make)
├── Makefile             build / clean targets
└── README.md            this file
```

---

