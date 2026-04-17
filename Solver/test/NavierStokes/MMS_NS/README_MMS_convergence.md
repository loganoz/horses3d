# MMS Convergence Verification — HORSES3D

This directory contains an automated Method of Manufactured Solutions (MMS) convergence study for HORSES3D 

---

## Directory layout

```
<solver>/
├── SETUP/
│   ├── generate_problemfile_<solver>.py   # symbolic source-term generator  ← edit here
│   └── Makefile                           # compiles ProblemFile.f90
├── MESH/
│   └── mesh<N>.h5                         # structured hex meshes (N³ elements)
├── CONTROL/
│   ├── mms_N<N>_P<P>.control             # auto-generated per run
│   └── mms_N<N>_P<P>.log                 # solver log per run
├── RESULTS/
│   └── MMS_N<N>_P<P>.hsol                # solver output
├── control_template.control               # solver settings  ← edit here
├── run_convergence_<solver>.py            # study configuration and orchestration  ← edit here
├── plot_convergence_<solver>.py           # plotting (optional edits)
└── errors.csv                             # accumulated results (auto-created)
```

Files under `CONTROL/` and `RESULTS/` are generated automatically — there is no reason to edit them by hand.

---

## Running a study

### Step 1 — Generate the problem file

```bash
cd <solver>/SETUP
python3 generate_problemfile_<solver>.py
make
```

This derives the MMS source terms symbolically (via SymPy), writes them into `ProblemFile.f90`, and compiles that file into the solver binary. You need to redo this step whenever you change the manufactured solution or any physical parameter.

### Step 2 — Prepare your meshes

Structured hexahedral meshes go in `MESH/`, named `mesh<N>.h5` where `N` is the number of elements per side — for example `mesh3.h5`, `mesh4.h5`, `mesh5.h5`, `mesh6.h5`.

### Step 3 — Configure the run

Open `run_convergence_<solver>.py` and set the paths and study parameters at the top of the file:

| Parameter | What it controls |
|---|---|
| `HORSES_BINARY` | Path to the compiled solver binary |
| `MESH_SIZES` | List of `N` values (elements per side) to sweep |
| `P_VALUES` | Polynomial orders to test |
| `DT` | Fixed time step (see the note on time-integration floors below) |
| `T_FINAL` | End time |

Other solver settings — Riemann solver choice, viscosity model, boundary conditions — live in `control_template.control`. At runtime the script fills in the following placeholders automatically:

| Token | Replaced with |
|---|---|
| `{MESH_FILE}` | Path to `mesh<N>.h5` |
| `{SOLUTION_FILE}` | Path to output `.hsol` |
| `{P}` | Polynomial order |
| `{DT}` | Time step |
| `{N_STEPS}` | `round(T_FINAL / DT)` |

### Step 4 — Run

```bash
python3 run_convergence_<solver>.py

# Preview what would run without executing anything:
python3 run_convergence_<solver>.py --dry-run
```

Each `(N, P)` pair appends one row to `errors.csv`. If the script is interrupted and restarted, it reads the CSV and skips cases that already finished — so you will not lose work.

### Step 5 — Plot

```bash
python3 plot_convergence_<solver>.py

# Override paths:
python3 plot_convergence_<solver>.py --csv path/to/errors.csv --outdir path/to/out
```

This produces `convergence.pdf` and `convergence.png`, each with three panels:

- **Left** — L² error vs. √NDOF, giving a combined h/p overview
- **Centre** — h-convergence per polynomial order in log-log, with best-fit slopes annotated against the theoretical p+1
- **Right** — p-convergence per mesh, with empirical reduction factors between successive orders

---

## What to expect from the results

For a smooth manufactured solution, a DG scheme of polynomial order p should give:

- **h-convergence** (refine the mesh at fixed p): slope ≈ p+1 in log(error) vs. log(h)
- **p-convergence** (increase p at fixed mesh): exponential error reduction

**Even-odd oscillations in h-slopes** are a known feature of the DGSEM on Gauss-Lobatto nodes. Super-convergence alternates between p+1 and p+2 depending on the parity of p, and its magnitude depends on the symmetry of your manufactured solution. This is not a defect.

**Time-integration floor:** if `DT` is too large, temporal error starts to dominate and the error curves flatten out. The plot marks these cases with faded markers (flagged when consecutive p-errors differ by less than a factor of 2). Fix: reduce `DT` by 10× and re-run.

---

## Changing the manufactured solution

Edit Section 2 of `generate_problemfile_<solver>.py`. Any smooth SymPy expression that is periodic on [−1, 1]³ will work. There are a few constraints to keep in mind per solver:

### Compressible NS

Pressure and density must both be strictly positive everywhere. The simplest way to enforce this is to add a constant offset large enough to dominate the oscillation amplitude. A non-solenoidal velocity field is allowed, but it introduces a large density source term that can make the system stiff and may require a smaller `DT`.

```python
rho = RHO0 * (1 + 0.1*sp.sin(pi*x)*sp.sin(pi*y)*sp.sin(pi*z)*sp.cos(t))
u   =  sp.sin(pi*x)*sp.cos(pi*y)*sp.cos(pi*z)*sp.sin(t)
v   =  sp.cos(pi*x)*sp.sin(pi*y)*sp.cos(pi*z)*sp.sin(t)
w   = -2*sp.cos(pi*x)*sp.cos(pi*y)*sp.sin(pi*z)*sp.sin(t)
p   =  1 + 0.1*sp.sin(pi*x)*sp.sin(pi*y)*sp.sin(pi*z)*sp.cos(t)
```

### Incompressible NS (iNS)

The velocity field must be exactly divergence-free: ∂u/∂x + ∂v/∂y + ∂w/∂z = 0. The easiest construction is to derive `u = curl(A)` for some vector potential `A`. The default field satisfies this exactly; verify any new field before running.

```python
# d/dx(sin·cos·cos) + d/dy(cos·sin·cos) + d/dz(-2cos·cos·sin) = 0 ✓
u =  sp.sin(pi*x)*sp.cos(pi*y)*sp.cos(pi*z)*sp.sin(t)
v =  sp.cos(pi*x)*sp.sin(pi*y)*sp.cos(pi*z)*sp.sin(t)
w = -2*sp.cos(pi*x)*sp.cos(pi*y)*sp.sin(pi*z)*sp.sin(t)
```

### Multiphase (MU)

The phase field must satisfy c ∈ [0, 1] pointwise. Using an amplitude |A| ≤ 1 around a mean of 0.5 guarantees this:

```python
c = 0.5*(1 + A*sp.cos(pi*x)*sp.cos(pi*y)*sp.cos(pi*z)*sp.sin(t))  # |A| <= 1
```

If you want clean p+1 slopes rather than the uniform super-convergence offset that the default fully-separable form produces, choose a `cMMS` that breaks the spatial symmetry.

> **Remember to regenerate and recompile after every change:** `python3 generate_problemfile_<solver>.py && make`

---

## Adding a body force

To verify a module that adds gravity or another body force, you only need to touch the generator. Find the source-term derivation block (look for the `# Source term derivation` comment) and add the force contribution to the appropriate equation. For example, to add gravitational acceleration in the z-direction to the NS momentum equations:

```python
g_val = 9.81   # m/s² or your non-dimensional equivalent
S_mom_z = S_mom_z - rho * g_val
```

SymPy will produce the correct modified source automatically. No other file needs changing.

---

## Verifying a new physics module

The general procedure is the same regardless of the module:

1. Identify which equations are affected and which terms change.
2. Add the corresponding SymPy expressions to the source-term derivation block in the generator.
3. Regenerate `ProblemFile.f90` and recompile.
4. Run the sweep and check the slopes.

Reading the results:

| What you see | Likely explanation |
|---|---|
| All slopes consistent with p+1 | Module is correctly implemented |
| Slopes systematically low across all (N, p) | Bug in the new module, or a mismatch between the generator parameters and the control file |
| Only the finest-mesh / highest-p points look low (faded markers) | Time-integration floor — reduce `DT` |
| Slopes erratic with no clear trend | Check that physical parameters are consistent between the generator and `control_template.control` |

---

## Keeping parameters consistent

The physical parameters in `generate_problemfile_<solver>.py` (things like `RHO0`, `GAMMA`, `RE`) **must match** what is set in `control_template.control`. A mismatch will not cause a crash — the solver will simply run with incorrect source terms, and the slopes will be wrong in ways that can look confusing.

---

## SymPy simplification

Section 3 of the generator exposes two settings:

| `SIMPLIFY_STRATEGY` | Behaviour |
|---|---|
| `"auto"` (default) | Tries `simplify()`, falls back to `trigsimp` after a timeout |
| `"simplify"` | Always uses `simplify()` — thorough but can be slow |
| `"trigsimp"` | Always uses `trigsimp()` — fast; works well for pure trigonometric solutions |

`SIMPLIFY_TIMEOUT` (default 60 s) controls how long `simplify()` is allowed to run before the fallback kicks in. Seeing `simplify timed out` on stderr is normal for complex expressions — set `SIMPLIFY_STRATEGY = "trigsimp"` to skip the wait entirely.

---

## Quick reference — what to edit

| File | Edit? | Purpose |
|---|---|---|
| `SETUP/generate_problemfile_*.py` | **Yes** — Sections 1–2 | Physical parameters, manufactured solution, new physics terms |
| `run_convergence_*.py` | **Yes** — Sections 1–2 | Paths, mesh sizes, polynomial orders, `DT`, `T_FINAL` |
| `control_template.control` | **Yes** | Solver-level settings (Riemann solver, BCs, viscosity model, etc.) |
| `plot_convergence_*.py` | Optional — Section 1 | Colour scheme, reference slopes, floor detection threshold |
| `SETUP/Makefile` | Rarely | Compilation flags |
| `errors.csv` | **Never** | Auto-generated; delete it if you want to force a clean sweep |
| `CONTROL/*.control` | **Never** | Auto-generated from the template |
| `CONTROL/*.log` | **Never** | Auto-generated solver logs |

---

## Troubleshooting

**`ERROR: template not found`**
Run the script from the directory that contains `control_template.control`, or update `TEMPLATE_FILE` in Section 1.

**`ERROR: could not read mms_l2_error.dat`**
The solver ran but `UserDefinedFinalize` did not write output. Check `CONTROL/mms_N<N>_P<P>.log`. The most common cause is that the binary was compiled without the correct preprocessor flag (`NAVIERSTOKES`, `CAHNHILLIARD`, etc.).

**Slopes consistently below p+1 for every case**
First, reduce `DT` by 10× and re-run. If slopes stay low, verify that the physical parameters in the generator match `control_template.control` exactly.

**Negative pressure or NaN during an NS run**
The manufactured pressure is not strictly positive somewhere. Increase the constant offset in the pressure expression or reduce the oscillation amplitude.

**`simplify timed out`** (printed to stderr during generation)
This is normal for complex expressions. Set `SIMPLIFY_STRATEGY = "trigsimp"` in Section 3 of the generator.

---

## Dependencies

| Package | Purpose | Install |
|---|---|---|
| Python ≥ 3.9 | Runner and plotter | — |
| SymPy | Symbolic source-term derivation | `pip install sympy` |
| NumPy | Quadrature weight computation | `pip install numpy` |
| Matplotlib | Convergence plots | `pip install matplotlib` |
| HORSES3D binary | Solver execution | Build from source |