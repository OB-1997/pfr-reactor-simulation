"""
FVM solver for steady-state convection-diffusion in a tube. 
Homework: process modelling - finite volume method

Governing equation: steady state, 1D, constant v and D:
    0 = -v dc/dz + D * (d2c/dz2)
Boundary conditions:
    c(0) = c0,  c(L) = cL
Analytical solution:
    (c(z) - c0) / (cL - c0) = (exp(vz/D) - 1) / (exp(vL/D) - 1)

"""

import numpy as np
import matplotlib.pyplot as plt

# physical parameters:

c0 = 1.0  # concentration at z = 0 [mol/m3]
cL = 0.0  # concentration at z = 0 [mol/m3]
L  = 1.0  # tube length [m]
D  = 0.1  # diffusion coefficient [m2/s]

# test cases as a list of dictionaries:

cases = [
    {"label": "v=0.1, N=5",  "v": 0.1, "N": 5},
    {"label": "v=0.5, N=5",  "v": 0.5, "N": 5},
    {"label": "v=2.5, N=5",  "v": 2.5, "N": 5},
    {"label": "v=2.5, N=20", "v": 2.5, "N": 20}
]

# analytical solution:

def analytical_solution(z,v):
    ratio = (np.exp(v * z / D) - 1.0) / (np.exp(v * L / D) - 1.0)
    return c0 + (cL - c0) * ratio

# FVM solver - central and upwind schemes:

def solve_fvm(v, N, scheme="central"):
    dz = L / N
    z = np.linspace(dz/2, L - dz/2, N)

    # check cell Peclet number:
    Pe_cell = v * dz / D

    # build the linear system A @ c = b
    A = np.zeros((N,N))
    b = np.zeros(N)

    for i in range(N):
        #  Left face of cell i:
        if i == 0:
            # boundary face: distance from c0 to c[0] is dz/2
            jd_in = D * (c0 - 0) / (dz / 2)
            # coefficient: D/(dz/2) * c0 goes to b, -D(dz/2) * c[0] goes to A
            A[i, i] += -D / (dz/2)       # from c[i]
            b[i]    += -D / (dz/2) * c0  # from c[0] (known, to RHS)
        else:
            # interior face: distance between centers is dz
            A[i, i]   += -D / dz  # from c[i]
            A[i, i-1] +=  D / dz  # from c[i-1]

        # Right face of cell i:
        if i == N - 1:
            # Boundary face: distance from c[N-1] to cL is dz/2
            A[i, i] += -D / (dz / 2)
            b[i]    += -D / (dz / 2) * cL
        else:
            # interior face 
            A[i, i]     += -D / dz
            A[i, i + 1] +=  D / dz
        
        # convective flux - depends on scheme choice
        if scheme == "central":
            # Left face (flux IN, positive = entering cell)
            if i == 0: 
                b[i] += -v * c0
            else:
                A[i, i - 1] += v * 0.5   
                A[i, i]     += v * 0.5    
        
            # Right face (flux OUT, negative = leaving cell)
            if i == N - 1:
                # jc_out = v * cL 
                b[i] += v * cL 
            else:
                # jc_out = v * 0.5*(c[i] + c[i+1])
                A[i, i]     += -v*0.5 
                A[i, i + 1] += -v*0.5
        
        elif scheme == "upwind":
            # Left face: fluid comes from the left, so value = left neighbor
            if i == 0: 
                b[i] += -v * c0
            else:
                A[i, i - 1] += v    
            
            # Right face: fluid goes right, so value = current cell 
            if i == N - 1:
                # jc_out = v * c[N-1]
                A[i, i] += -v
            else:
                A[i, i] += -v
    c = np.linalg.solve(A,b)

    return z, c

# solution and graphical interpretation:

results = []

# === thesis-wide unified plot style (matches Fig 2.1, rtd_cascade) ===
plt.rcParams.update({
    "font.size":        11,
    "axes.titlesize":   11,
    "axes.labelsize":   11,
    "xtick.labelsize":  9,
    "ytick.labelsize":  9,
    "legend.fontsize":  9,
    "axes.linewidth":   0.8,
    "lines.linewidth":  1.4,
    "lines.markersize": 5,
    "grid.linestyle":   ":",
    "grid.linewidth":   0.5,
    "grid.alpha":       0.4,
    "savefig.bbox":     "tight",
})

# Sequential viridis sample for central / upwind curves; analytical stays black.
_cmap = plt.get_cmap("viridis")
COLOUR_CENTRAL = _cmap(0.30)
COLOUR_UPWIND  = _cmap(0.70)
# === end unified style ===

fig, axes = plt.subplots(2, 2, figsize=(10.0, 6.0))
axes_flat = axes.flatten()

for idx, case in enumerate(cases):
    v  = case["v"]
    N  = case["N"]
    ax = axes_flat[idx]

    # analytical solution:
    z_fine  = np.linspace(0,L,500)
    c_exact = analytical_solution(z_fine, v)

    # FVM with central scheme:
    z_c, c_central = solve_fvm(v, N, scheme="central")

    # FVM with upwind scheme:
    z_u, c_upwind = solve_fvm(v, N, scheme="upwind")

    # Analytical at cell centres (for error calculation)
    c_exact_at_cells = analytical_solution(z_c, v)
    err_central      = np.max(np.abs(c_central - c_exact_at_cells))
    err_upwind       = np.max(np.abs(c_upwind  - c_exact_at_cells))

    Pe_cell = v * (L / N) / D
    results.append({
        "label":       case["label"],
        "Pe_cell":     Pe_cell,
        "err_central": err_central,
        "err_upwind":  err_upwind,
    })

    # plot
    ax.plot(z_fine, c_exact,
            color="black", linestyle="-", linewidth=1.0,
            label="analytical")
    ax.plot(z_c, c_central,
            color=COLOUR_CENTRAL, linestyle="--", linewidth=1.2,
            marker="o", markersize=5,
            markerfacecolor="white", markeredgewidth=1.2,
            label="central")
    ax.plot(z_u, c_upwind,
            color=COLOUR_UPWIND, linestyle="--", linewidth=1.2,
            marker="s", markersize=5,
            markerfacecolor="white", markeredgewidth=1.2,
            label="upwind")

    ax.set_xlabel(r"axial position $z$ / m")
    ax.set_ylabel(r"concentration / (mol m$^{-3}$)")
    ax.set_title(case["label"])
    ax.legend(frameon=False)
    ax.set_ylim(-0.2, 2.5)
    ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.4)

plt.tight_layout()
plt.savefig("fvm_comparison.png", dpi=200)

print()
print("=" * 70)
print(f"{'Case':<18} {'Pe_cell':>8} {'Err(central)':>14} {'Err(upwind)':>14}")
print("-" * 70)
for r in results:
    flag = " !!!" if r["Pe_cell"] >= 2 else ""
    print(f"{r['label']:<18} {r['Pe_cell']:>8.2f} {r['err_central']:>14.4e} {r['err_upwind']:>14.4e}{flag}")
print("-" * 70)
print("!!!  = Pe_cell >= 2, central scheme unstable")
print("=" * 70)

plt.show()
print(f"\nSaved: fvm_comparison.png")












