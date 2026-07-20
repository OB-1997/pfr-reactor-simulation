"""
validate.py - programmatic checks for the simple A->B PFR build.

Five checks (matches the validation list in the implementation plan):
  1. build              (assumed already done; this script just verifies output exists)
  2. smoke run          (assumed already done; this script reads the run dirs)
  3. mass conservation  (max |C_A(j) + C_B(j) - C_A,in| at steady state)
  4. steady state       (last vs second-to-last sample agree)
  5. analytical match   (max-norm rel-err vs analytical reference at N=1000)

Two analytical references are provided:
  - Dirichlet/Neumann (D/N): the exact steady-state solution of the same ODE
    with the same BCs the discrete FVM enforces. This is what the simulation
    converges to in the dx -> 0 limit, so it is the "right" reference for
    measuring spatial discretisation error.
  - Wehner-Wilhelm (W-W) closed-closed Danckwerts (1956): the academic
    closed-form often quoted in textbooks. At large Pe (our regime, Pe~10^3)
    it differs from D/N by O(1/Pe) ~ 10^-3 — at the same level as the
    upwind discretisation error. Useful as an independent sanity check.

Run from the project root:

    python3 validate.py output/odepack_1 output/odepack_2

The first run dir is the smoke check (any solver / N), the second is the N=1000
DLSODE reference run that gets compared to the analytical solutions.
"""

import sys
import math
from pathlib import Path

import numpy as np


# ----------------------------------------------------------- analytical


def analytical_dirichlet_neumann(x_star, pe, da):
    """
    Steady-state C(x)/C_in for the dimensionless ODE
        d^2C/dz^2 - Pe dC/dz - Pe Da C = 0    (z = x/L)
    with the BCs the discrete FVM enforces:
        C(0)        = C_in        (Dirichlet inlet)
        dC/dz |_1   = 0           (Neumann zero-gradient outlet)

    Characteristic roots: m+- = Pe (1 +- beta) / 2, beta = sqrt(1 + 4 Da/Pe).
    For Pe >> 1 we have m+ very large positive, m- O(Da) negative; the m+
    branch is essentially absent because A = -B m_- exp(m_- - m_+) / m_+
    underflows to zero. The implementation handles both regimes.
    """
    beta = math.sqrt(1.0 + 4.0 * da / pe)
    m_plus  = pe * (1.0 + beta) / 2.0
    m_minus = pe * (1.0 - beta) / 2.0    # <= 0 for Da >= 0

    # Solve [1, 1; m+ exp(m+), m- exp(m-)] [A; B] = [1; 0].
    # Equivalent stable form: factor out exp(m+) so coefficients stay finite.
    # X := -m- exp(m- - m+) / m+ ;   A = X/(1+X) ;   B = 1/(1+X).
    log_factor = m_minus - m_plus        # = -beta * Pe (very negative for Pe >> 1)
    x = np.asarray(x_star, dtype=float)

    if log_factor < -700.0:
        # m_+ branch underflows; A is exactly zero in double precision.
        return np.exp(m_minus * x)

    X = -m_minus * math.exp(log_factor) / m_plus
    one_over = 1.0 / (1.0 + X)

    # B exp(m- z) — well-scaled directly.
    bexp = one_over * np.exp(m_minus * x)

    # A exp(m+ z) = (X / (1+X)) exp(m+ z) but X already contains exp(m- - m+),
    # so the net exponent is m_+ (z - 1) + m_- (always <= 0 for z in [0, 1]).
    # Rewrite:   A exp(m+ z) = -m-/m+ * exp(m- + m+ (z - 1)) / (1 + X).
    aexp = (-m_minus / m_plus) * np.exp(m_minus + m_plus * (x - 1.0)) * one_over

    return aexp + bexp


def analytical_wehner_wilhelm(x_star, pe, da):
    """
    Closed-closed Danckwerts BCs (Wehner & Wilhelm 1956):
        z = 0:  u C_in = u C(0) - D dC/dz |_0      (continuity of total flux)
        z = 1:  dC/dz |_1 = 0
    For the same ODE as above. Different inlet treatment than the discrete
    model (small back-diffusion drop at z = 0+).

    Same characteristic roots m+- = Pe (1 +- beta)/2.
    BC system:
        A (1 - m+/Pe) + B (1 - m-/Pe) = 1
        A m+ exp(m+)  + B m- exp(m-)  = 0
    Note 1 - m_+/Pe = (1 - beta)/2 = m_-/Pe and 1 - m_-/Pe = m_+/Pe.
    """
    beta = math.sqrt(1.0 + 4.0 * da / pe)
    m_plus  = pe * (1.0 + beta) / 2.0
    m_minus = pe * (1.0 - beta) / 2.0

    log_factor = m_minus - m_plus

    if log_factor < -700.0:
        # A -> 0; from the inlet BC,   B (1 - m-/Pe) = 1   =>   B = 2/(1+beta).
        B = 2.0 / (1.0 + beta)
        x = np.asarray(x_star, dtype=float)
        return B * np.exp(m_minus * x)

    # General case (kept for low-Pe sanity).
    expm_minus = math.exp(m_minus)
    expm_plus  = math.exp(m_plus)
    M = np.array([[1.0 - m_plus / pe, 1.0 - m_minus / pe],
                  [m_plus * expm_plus, m_minus * expm_minus]])
    rhs = np.array([1.0, 0.0])
    A, B = np.linalg.solve(M, rhs)
    x = np.asarray(x_star, dtype=float)
    return A * np.exp(m_plus * x) + B * np.exp(m_minus * x)


# ----------------------------------------------------------- run-dir IO


def _read_n_cells(run_dir: Path) -> int:
    text = (run_dir / "parameters.txt").read_text()
    for line in text.splitlines():
        if "N_cells" in line and "=" in line:
            return int(line.split("=", 1)[1].strip())
    raise RuntimeError(f"N_cells not found in {run_dir}/parameters.txt")


def _read_run(run_dir: Path):
    """
    Parse output/<solver>_<N>/output.txt into (t, C) where t.shape = (Nstep,)
    and C.shape = (Nstep, N_cells, N_species). Fortran y((j-1)*N_species + i)
    means species varies fastest, cell varies slowest.
    """
    n_cells = _read_n_cells(run_dir)
    raw = np.loadtxt(run_dir / "output.txt")
    t = raw[:, 0]
    y = raw[:, 1:]
    n_species = y.shape[1] // n_cells
    assert n_species == 2, f"expected N_species=2, got {n_species}"
    C = y.reshape(y.shape[0], n_cells, n_species)
    return t, C, n_cells


# ----------------------------------------------------------------- checks


def check_mass_conservation(C, c_a_in, label, tol=1.0e-4):
    """
    Final-time mass conservation: per-cell |C_A + C_B - C_A_in|.
    Tolerance is set well above DLSODE's default rtol*C_in ~ 1e-4 so we catch
    real bugs but tolerate ODE-tolerance-limited approach to steady state.
    """
    final = C[-1]
    s = final[:, 0] + final[:, 1]
    err = np.max(np.abs(s - c_a_in))
    ok = err < tol
    print(f"  [{label}] max |C_A + C_B - C_A_in| at t_end = {err:.3e}  -> {'OK' if ok else 'FAIL'} (tol {tol:.0e})")
    return ok


def check_steady_state(t, C, label, tol=1.0e-4):
    """
    Last vs second-to-last sample relative diff.
    """
    diff = np.max(np.abs(C[-1] - C[-2]))
    ok = diff < tol
    print(f"  [{label}] max |C(t_end) - C(t_end-dt)| = {diff:.3e}     "
          f"({t[-1]-t[-2]:.2f} s apart) -> {'OK' if ok else 'FAIL'} (tol {tol:.0e})")
    return ok


def check_analytical(C, n_cells, label, l_reactor=20.0, u_vel=0.516, d_ax=1e-2,
                     k_rxn=0.026, c_a_in=10.4, tol=5.0e-3):
    pe = u_vel * l_reactor / d_ax
    da = k_rxn * l_reactor / u_vel
    print(f"  [{label}] Pe = {pe:.2f}, Da = {da:.4f}")

    j = np.arange(1, n_cells + 1)
    x_star = (j - 0.5) / n_cells

    c_a_dn = c_a_in * analytical_dirichlet_neumann(x_star, pe, da)
    c_a_ww = c_a_in * analytical_wehner_wilhelm(x_star, pe, da)
    c_a_pf = c_a_in * np.exp(-da * x_star)
    c_a_sim = C[-1, :, 0]

    mid = n_cells // 2
    print(f"  [{label}] sample (mid cell, x* = {x_star[mid]:.4f}):")
    print(f"           sim                = {c_a_sim[mid]:.4f}")
    print(f"           Dirichlet/Neumann  = {c_a_dn[mid]:.4f}     (matches discrete BCs)")
    print(f"           Wehner-Wilhelm     = {c_a_ww[mid]:.4f}     (closed-closed Danckwerts)")
    print(f"           plug-flow limit    = {c_a_pf[mid]:.4f}     (Pe -> infinity)")

    rel_dn = np.abs(c_a_sim - c_a_dn) / c_a_dn
    rel_ww = np.abs(c_a_sim - c_a_ww) / c_a_ww
    print(f"  [{label}] max-norm rel error vs Dirichlet/Neumann = {np.max(rel_dn):.3e}")
    print(f"  [{label}] max-norm rel error vs Wehner-Wilhelm    = {np.max(rel_ww):.3e}")
    ok = np.max(rel_dn) < tol
    print(f"  -> {'OK' if ok else 'FAIL'} (tol {tol:.0e}; plan expects ~1e-3 at N=1000)")
    return ok


# ------------------------------------------------------------------ main


def main(argv):
    if len(argv) < 2:
        print(__doc__.strip())
        return 1

    smoke_dir = Path(argv[1])
    print(f"\n=== Smoke run ({smoke_dir}) ===")
    t, C, n_cells = _read_run(smoke_dir)
    print(f"  Nstep = {len(t)}, N_cells = {n_cells}, t in [{t[0]:.2f}, {t[-1]:.2f}] s")
    ok1 = check_mass_conservation(C, c_a_in=10.4, label="smoke")
    ok2 = check_steady_state(t, C, label="smoke")

    if len(argv) > 2:
        ref_dir = Path(argv[2])
        print(f"\n=== Reference run ({ref_dir}) ===")
        t, C, n_cells = _read_run(ref_dir)
        print(f"  Nstep = {len(t)}, N_cells = {n_cells}, t in [{t[0]:.2f}, {t[-1]:.2f}] s")
        ok3 = check_mass_conservation(C, c_a_in=10.4, label="ref")
        ok4 = check_steady_state(t, C, label="ref")
        ok5 = check_analytical(C, n_cells, label="ref")
        all_ok = ok1 and ok2 and ok3 and ok4 and ok5
    else:
        all_ok = ok1 and ok2

    print()
    print(f"OVERALL: {'PASS' if all_ok else 'FAIL'}")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
