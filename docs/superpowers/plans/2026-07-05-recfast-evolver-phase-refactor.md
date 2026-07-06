# RECFAST evolver phase-mode refactor — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Collapse PR #362's two RECFAST derivative-mode enums onto master's phase-free RHS, carry a single `RecfastPhase` in the workspace behind one derivs shim and one output shim, and pick `tol_thermo_integration` by a convergence study.

**Architecture:** The physics RHS returns to master's `thermodynamics_derivs_with_recfast_member` verbatim (its `x_H > trigger` and `x_He ∈ [5e-9, trigger2]` state tests already produce the correct frozen/evolved behavior in every phase). The evolver-facing plumbing dispatches on `ws->recfast_phase ∈ {helium, full}`; the analytic regime keeps calling the RHS directly for the `cb2` column. Then `tol_thermo_integration` is set from a sweep against a `1e-9` reference.

**Tech Stack:** C++17 (CLASSpp), the shared `evolver_ndf15`/`evolver_rkdp45`/`evolver_rk` interface, `make` build producing `./class`, Python+NumPy for output comparison.

## Global Constraints

- **C++-only, no C-compat guards** — plain C++ headers (no `extern "C"`, no `#ifdef __cplusplus`).
- **Never `git add -A`** — the working tree has many untracked build artifacts (`build_prof/`, `generated/`, `*.pdf`, `/private/tmp/*.ini`). Stage explicit paths only.
- **`cclassy.pxd` is auto-generated** — do not hand-edit; not touched by this plan anyway.
- **Default evolver is `ndf15`** (`evolver_type{rk=0, ndf15=1, rkdp45=2}`; `precision.h` default `evolver_type::ndf15`).
- **Verification altitude:** the *structural refactor* (Task 2) must be inert — agreement with the pre-refactor branch tip to machine precision (`max|rel| < 1e-10`). The *tol change* (Task 3) is a deliberate output change, judged by the convergence criterion below and a ~0.1% CMB-TT sanity check vs master.
- **Scratchpad root** (temp inis/outputs/scripts, never committed):
  `/private/tmp/claude-497981304/-Users-au192734-Projects-class-claude/09d6c8ce-241a-4b7e-8990-17a794a6aa03/scratchpad`

---

## File Structure

- `source/thermodynamics_module.h` — replace the two mode enums + the phase-specific method decls with `RecfastPhase` + two shim decls; add `recfast_phase` to the workspace struct.
- `source/thermodynamics_module.cpp` — revert the RHS to master; replace three static derivs shims + three output callbacks with one of each; rewire the driver; drop the dead zero-sample path.
- `include/precision.h` — set `tol_thermo_integration` to the swept value (Task 3 only).
- Scratchpad (test harness, not committed): comparison `.ini`s and `compare_thermo.py`.

---

### Task 1: A/B baseline + comparison harness

Capture the current branch tip's thermodynamics output as the golden for the inert-refactor check, and build the numeric comparison tool. **No source changes in this task.**

**Files:**
- Create: `<scratch>/ab/base_ndf15.ini`, `<scratch>/ab/base_rkdp45.ini`
- Create: `<scratch>/compare_thermo.py`

**Interfaces:**
- Produces: `<scratch>/ab/ndf15_baseline.dat`, `<scratch>/ab/rkdp45_baseline.dat` (RECFAST thermodynamics tables from the pre-refactor code), consumed by Task 2.
- Produces: `compare_thermo.py <fileA> <fileB>` → prints per-column `max|rel|`, an `OVERALL max|rel|`, and windowed (`800<z<3000`) `x_e`/`Tb` deviations; exits non-zero on shape mismatch.

- [ ] **Step 1: Confirm you are on the branch tip, not detached**

Run: `git -C /Users/au192734/Projects/class_claude status -sb | head -1`
Expected: `## recfast-evolver-thermodynamics` (or ahead of it by the spec commits). If not, stop.

- [ ] **Step 2: Write the two baseline inis**

```bash
SCRATCH=/private/tmp/claude-497981304/-Users-au192734-Projects-class-claude/09d6c8ce-241a-4b7e-8990-17a794a6aa03/scratchpad
mkdir -p "$SCRATCH/ab"
cat > "$SCRATCH/ab/base_ndf15.ini" <<EOF
recombination = RECFAST
write thermodynamics = yes
root = $SCRATCH/ab/ndf15_
EOF
cat > "$SCRATCH/ab/base_rkdp45.ini" <<EOF
recombination = RECFAST
write thermodynamics = yes
evolver = 2
root = $SCRATCH/ab/rkdp45_
EOF
```

- [ ] **Step 3: Write the comparison script**

```python
# <scratch>/compare_thermo.py
import sys
import numpy as np

a = np.loadtxt(sys.argv[1])
b = np.loadtxt(sys.argv[2])
if a.shape != b.shape:
    print(f"SHAPE MISMATCH {a.shape} vs {b.shape}")
    sys.exit(1)

cols = ["z", "conf.time", "x_e", "kappa'", "exp(-kappa)", "g", "Tb", "w_b", "c_b^2", "tau_d"]
den = np.abs(b) + 1e-300
rel = np.abs(a - b) / den
for j in range(a.shape[1]):
    name = cols[j] if j < len(cols) else f"col{j}"
    print(f"{name:12s} max|rel|={rel[:, j].max():.3e}")
overall = rel.max()
print(f"OVERALL      max|rel|={overall:.3e}")

# windowed recombination check on x_e (col 2) and Tb (col 6)
z = b[:, 0]
m = (z > 800.0) & (z < 3000.0)
if m.any():
    print(f"[800<z<3000] x_e max|rel|={rel[m, 2].max():.3e}  "
          f"Tb max|rel|={rel[m, 6].max():.3e}")
```

- [ ] **Step 4: Build the current branch tip**

Run: `cd /Users/au192734/Projects/class_claude && make -j`
Expected: build succeeds, `./class` present. (If the build fails on the untouched tip, stop — that is a pre-existing problem, not this work.)

- [ ] **Step 5: Generate the baselines and snapshot them**

```bash
cd /Users/au192734/Projects/class_claude
SCRATCH=/private/tmp/claude-497981304/-Users-au192734-Projects-class-claude/09d6c8ce-241a-4b7e-8990-17a794a6aa03/scratchpad
./class "$SCRATCH/ab/base_ndf15.ini"
./class "$SCRATCH/ab/base_rkdp45.ini"
cp "$SCRATCH/ab/ndf15_thermodynamics.dat"  "$SCRATCH/ab/ndf15_baseline.dat"
cp "$SCRATCH/ab/rkdp45_thermodynamics.dat" "$SCRATCH/ab/rkdp45_baseline.dat"
wc -l "$SCRATCH/ab/ndf15_baseline.dat" "$SCRATCH/ab/rkdp45_baseline.dat"
```
Expected: both runs print "Computing thermodynamics"; both `*_baseline.dat` exist with the same (nonzero) line count.

- [ ] **Step 6: Self-check the harness (compare a baseline to itself)**

Run: `python3 "$SCRATCH/compare_thermo.py" "$SCRATCH/ab/ndf15_baseline.dat" "$SCRATCH/ab/ndf15_baseline.dat"`
Expected: every line `max|rel|=0.000e+00`, `OVERALL max|rel|=0.000e+00`.

No commit — this task produces only scratchpad artifacts.

---

### Task 2: Revert RHS to master, consolidate onto `RecfastPhase`, rewire driver

> **Implemented differently — see commit `01348d56` and the revised design spec §1.**
> The "revert the RHS to master verbatim" approach in the steps below was found
> **non-inert**: master's `x_He > x_He0_trigger2` Heflag gate chatters at the helium
> seed and moves ndf15 `c_b^2` by ~3e-3. The shipped design instead moves the two
> double-encoded regime switches (H-freeze, He-corrections) onto the phase and lets
> the RHS read them as flags (`hydrogen_frozen = phase!=full`,
> `helium_corrections = phase!=analytic`), plus the `He_Boltz` `std::min` tidy. That
> is bit-identically inert (thermo A/B `max|rel|=0`, both evolvers). The steps below
> are retained for provenance; Steps 1–4, 8 (shim/driver plumbing) landed as written.

One coherent change: the header, the RHS, the shims, and the driver must compile together. After it, the RECFAST output must match Task 1's baselines to machine precision.

**Files:**
- Modify: `source/thermodynamics_module.h` (method-decl block ~124-166; workspace struct ~201-210)
- Modify: `source/thermodynamics_module.cpp` (static shims ~42-99; driver helium/full section ~3518-3580; RHS members ~3605-3894)

**Interfaces:**
- Produces: `enum class RecfastPhase { helium, full };` (file scope in the header, before the workspace struct).
- Produces: workspace field `RecfastPhase recfast_phase = RecfastPhase::full;`.
- Produces: `static void thermodynamics_recfast_derivs(double minus_z, double* y, double* dy, void*)` and `static void thermodynamics_recfast_output(double minus_z, double y[], double dy[], int index_x, void*)`.
- Produces: `void thermodynamics_derivs_with_recfast_member(double z, double* y, double* dy, void*)` with master's body (no mode parameters).
- Removed: `RecfastHydrogenMode`, `RecfastHeliumMode`, `thermodynamics_recfast_derivs_member`, `thermodynamics_derivs_with_recfast_{full,helium}_member`, `thermodynamics_derivs_with_recfast_{minus_z,helium}`, `thermodynamics_recfast_output_{full,helium,none}`.

- [ ] **Step 1: Header — replace the enums + phase-method decls**

In `source/thermodynamics_module.h`, replace this block (the two enums are the first two lines after `thermodynamics_recombination_with_recfast(...)`, through the last `..._output_full` decl):

```cpp
  enum class RecfastHydrogenMode { legacy_trigger, frozen, evolved };
  enum class RecfastHeliumMode { legacy_trigger, evolved };
```
Delete both of the above lines.

Then, further down, replace the derivs/shim/output decl group. Delete these:

```cpp
  void thermodynamics_recfast_derivs_member(double z,
                                            double* y,
                                            double* dy,
                                            void* fixed_parameters,
                                            RecfastHydrogenMode hydrogen_mode,
                                            RecfastHeliumMode helium_mode);
```
```cpp
  void thermodynamics_derivs_with_recfast_full_member(double minus_z,
                                                      double* y,
                                                      double* dy,
                                                      void* fixed_parameters);
  void thermodynamics_derivs_with_recfast_helium_member(double minus_z,
                                                        double* y,
                                                        double* dy,
                                                        void* fixed_parameters);
  static void thermodynamics_derivs_with_recfast_minus_z(double minus_z,
                                                         double* y,
                                                         double* dy,
                                                         void* fixed_parameters);
  static void thermodynamics_derivs_with_recfast_helium(double minus_z,
                                                        double* y,
                                                        double* dy,
                                                        void* fixed_parameters);
```
```cpp
  static void thermodynamics_recfast_output_helium(
      double minus_z, double y[], double dy[], int index_x, void* parameters_and_workspace);
  static void thermodynamics_recfast_output_none(
      double minus_z, double y[], double dy[], int index_x, void* parameters_and_workspace);
  static void thermodynamics_recfast_output_full(
      double minus_z, double y[], double dy[], int index_x, void* parameters_and_workspace);
```

Keep `thermodynamics_derivs_with_recfast_member(...)` (member) and the static `thermodynamics_derivs_with_recfast(...)` and `thermodynamics_recfast_timescale(...)`. Immediately after the static `thermodynamics_derivs_with_recfast(...)` decl, add the two new shims:

```cpp
  static void thermodynamics_recfast_derivs(double minus_z,
                                            double* y,
                                            double* dy,
                                            void* fixed_parameters);
  static void thermodynamics_recfast_output(
      double minus_z, double y[], double dy[], int index_x, void* parameters_and_workspace);
```

- [ ] **Step 2: Header — add the file-scope enum and the workspace field**

In `source/thermodynamics_module.h`, between the class's closing `};` and the `struct thermodynamics_parameters_and_workspace` doc comment, add:

```cpp
enum class RecfastPhase { helium, full };
```

Then add the field to the struct (after `recfast_output_index_offset`):

```cpp
  int recfast_output_index_offset = 0;
  RecfastPhase recfast_phase       = RecfastPhase::full;
```

- [ ] **Step 3: cpp — replace the three static derivs shims + three output callbacks with one of each**

In `source/thermodynamics_module.cpp`, delete `thermodynamics_derivs_with_recfast_minus_z`, `thermodynamics_derivs_with_recfast_helium`, `thermodynamics_recfast_output_helium`, `thermodynamics_recfast_output_none`, and `thermodynamics_recfast_output_full` (keep `thermodynamics_derivs_with_recfast` and `thermodynamics_recfast_timescale`). In their place add:

```cpp
// Single evolver derivative shim: integrates in minus_z = -z, dispatches on phase.
void ThermodynamicsModule::thermodynamics_recfast_derivs(double minus_z,
                                                         double* y,
                                                         double* dy,
                                                         void* parameters_and_workspace) {
  auto ws        = static_cast<thermodynamics_parameters_and_workspace*>(parameters_and_workspace);
  auto* module   = ws->thermodynamics_module;
  const double z = -minus_z;

  if (ws->recfast_phase == RecfastPhase::helium) {
    // State is {x_He, Tmat}; hydrogen follows its Saha branch.
    double y_full[3] = {module->thermodynamics_recfast_hydrogen_saha_xH(ws->preco, z), y[0], y[1]};
    double dy_full_dz[3];
    module->thermodynamics_derivs_with_recfast_member(z,
                                                      y_full,
                                                      dy_full_dz,
                                                      parameters_and_workspace);
    dy[0] = -dy_full_dz[1];
    dy[1] = -dy_full_dz[2];
  }
  else {
    // State is {x_H, x_He, Tmat}.
    double dy_dz[_RECFAST_INTEG_SIZE_];
    module->thermodynamics_derivs_with_recfast_member(z, y, dy_dz, parameters_and_workspace);
    for (int index_y = 0; index_y < _RECFAST_INTEG_SIZE_; index_y++)
      dy[index_y] = -dy_dz[index_y];
  }
}

// Single dense-output shim: stores one recombination-table row, dispatching on phase.
void ThermodynamicsModule::thermodynamics_recfast_output(
    double minus_z, double y[], double dy[], int index_x, void* parameters_and_workspace) {
  auto ws                = static_cast<thermodynamics_parameters_and_workspace*>(parameters_and_workspace);
  auto* module           = ws->thermodynamics_module;
  recombination* preco   = ws->preco;
  const double z         = -minus_z;
  const int sample_index = ws->recfast_output_index_offset + index_x;

  if (ws->recfast_phase == RecfastPhase::helium) {
    const double y_full[3] = {module->thermodynamics_recfast_hydrogen_saha_xH(preco, z), y[0], y[1]};
    const double xe        = module->thermodynamics_recfast_xe_after_helium_ode(preco, z, y_full);
    module->thermodynamics_recfast_store_row(preco, sample_index, z, xe, y[1], -dy[1]);
  }
  else {
    const double xe = module->thermodynamics_recfast_xe_after_full_ode(preco, z, y);
    module->thermodynamics_recfast_store_row(preco, sample_index, z, xe, y[2], -dy[2]);
  }
}
```

- [ ] **Step 4: cpp — delete the two member adapters and the thin member wrapper**

Delete these three functions entirely (they sit just above `thermodynamics_recfast_derivs_member`):
- `void ThermodynamicsModule::thermodynamics_derivs_with_recfast_helium_member(...)`
- `void ThermodynamicsModule::thermodynamics_derivs_with_recfast_member(double z, double* y, double* dy, void* parameters_and_workspace) { thermodynamics_recfast_derivs_member(..., legacy_trigger, legacy_trigger); }` (the thin one)
- `void ThermodynamicsModule::thermodynamics_derivs_with_recfast_full_member(...)`

- [ ] **Step 5: cpp — rename the core RHS and drop its mode parameters**

Change the signature of `thermodynamics_recfast_derivs_member` to master's name and parameter list:

```cpp
void ThermodynamicsModule::thermodynamics_recfast_derivs_member(double z,
                                                                double* y,
                                                                double* dy,
                                                                void* parameters_and_workspace,
                                                                RecfastHydrogenMode hydrogen_mode,
                                                                RecfastHeliumMode helium_mode) {
```
becomes
```cpp
void ThermodynamicsModule::thermodynamics_derivs_with_recfast_member(double z,
                                                                     double* y,
                                                                     double* dy,
                                                                     void* parameters_and_workspace) {
```

- [ ] **Step 6: cpp — restore master's Heflag gate (remove the helium mode branch)**

Inside that function, replace:

```cpp
  int Heflag = 0;
  if (helium_mode == RecfastHeliumMode::evolved) {
    if (x_He >= 5.e-9)
      Heflag = ppr->recfast_Heswitch;
  }
  else {
    if ((x_He >= 5.e-9) && (x_He <= ppr->recfast_x_He0_trigger2))
      Heflag = ppr->recfast_Heswitch;
  }
```
with master's state-only gate:
```cpp
  int Heflag = 0;
  if ((x_He < 5.e-9) || (x_He > ppr->recfast_x_He0_trigger2))
    Heflag = 0;
  else
    Heflag = ppr->recfast_Heswitch;
```

- [ ] **Step 7: cpp — restore master's hydrogen freeze test (remove the hydrogen mode branch)**

Replace:
```cpp
  if ((hydrogen_mode == RecfastHydrogenMode::frozen) ||
      ((hydrogen_mode == RecfastHydrogenMode::legacy_trigger) && (x_H > ppr->recfast_x_H0_trigger)))
    dy[0] = 0.;
  else {
```
with:
```cpp
  if (x_H > ppr->recfast_x_H0_trigger)
    dy[0] = 0.;
  else {
```

- [ ] **Step 8: cpp — rewire the driver helium + full phases**

In `thermodynamics_recombination_with_recfast`, replace the helium+full section (from the `y[0] = thermodynamics_recfast_hydrogen_saha_xH(preco, z_helium_ode_start);` seed through the end of the `if (first_full_sample < Nz) { ... }` block) with:

```cpp
  /** - Evolve Helium and baryon temperature while Hydrogen follows its Saha branch. */

  y[0] = thermodynamics_recfast_hydrogen_saha_xH(preco, z_helium_ode_start);
  y[1] = ppr->recfast_x_He0_trigger;
  y[2] = preco->Tnow * (1. + z_helium_ode_start);

  const int first_helium_sample = first_sample_below(i, z_helium_ode_start);
  const int helium_sample_count = sample_count_down_to(first_helium_sample, z_hydrogen_ode_start);

  class_test(helium_sample_count <= 0,
             "RECFAST helium phase contains no sampled redshifts; increase recfast_Nz0");

  double y_helium[2]               = {y[1], y[2]};
  tpaw.recfast_phase               = RecfastPhase::helium;
  tpaw.recfast_output_index_offset = first_helium_sample;

  generic_evolver(thermodynamics_recfast_derivs,
                  -z_helium_ode_start,
                  -z_hydrogen_ode_start,
                  y_helium,
                  used_in_output_helium.data(),
                  2,
                  &tpaw,
                  ppr->tol_thermo_integration,
                  ppr->smallest_allowed_variation,
                  thermodynamics_recfast_timescale,
                  z_helium_ode_start - z_hydrogen_ode_start,
                  minus_z_sampling.data() + first_helium_sample,
                  helium_sample_count,
                  thermodynamics_recfast_output,
                  nullptr);

  y[0] = ppr->recfast_x_H0_trigger;
  y[1] = y_helium[0];
  y[2] = y_helium[1];

  /** - Evolve Hydrogen, Helium, and baryon temperature with dense-output sampling. */

  const int first_full_sample = first_sample_below(first_helium_sample + helium_sample_count,
                                                   z_hydrogen_ode_start);

  if (first_full_sample < Nz) {
    tpaw.recfast_phase               = RecfastPhase::full;
    tpaw.recfast_output_index_offset = first_full_sample;

    generic_evolver(thermodynamics_recfast_derivs,
                    -z_hydrogen_ode_start,
                    minus_z_sampling[Nz - 1],
                    y,
                    used_in_output_full.data(),
                    _RECFAST_INTEG_SIZE_,
                    &tpaw,
                    ppr->tol_thermo_integration,
                    ppr->smallest_allowed_variation,
                    thermodynamics_recfast_timescale,
                    z_hydrogen_ode_start,
                    minus_z_sampling.data() + first_full_sample,
                    Nz - first_full_sample,
                    thermodynamics_recfast_output,
                    nullptr);
  }
```

(This deletes the `helium_boundary_sample`, `helium_sampling`, `helium_sampling_count`, and `helium_output` locals.)

- [ ] **Step 9: Build**

Run: `cd /Users/au192734/Projects/class_claude && make -j`
Expected: compiles clean. If the compiler reports an undeclared `RecfastPhase` in the workspace struct, the file-scope enum from Step 2 is misplaced (it must precede the struct). If it reports leftover references to a deleted symbol, grep for it: `grep -n "RecfastHydrogenMode\|RecfastHeliumMode\|_derivs_member\|_output_none\|_output_full\|_output_helium\|_minus_z\|_full_member\|_helium_member\|_with_recfast_helium\b" source/thermodynamics_module.*` should return nothing in the cpp/h except the decls you intend.

- [ ] **Step 10: Verify the refactor is inert vs the Task 1 baselines**

```bash
cd /Users/au192734/Projects/class_claude
SCRATCH=/private/tmp/claude-497981304/-Users-au192734-Projects-class-claude/09d6c8ce-241a-4b7e-8990-17a794a6aa03/scratchpad
./class "$SCRATCH/ab/base_ndf15.ini"
./class "$SCRATCH/ab/base_rkdp45.ini"
python3 "$SCRATCH/compare_thermo.py" "$SCRATCH/ab/ndf15_thermodynamics.dat"  "$SCRATCH/ab/ndf15_baseline.dat"
python3 "$SCRATCH/compare_thermo.py" "$SCRATCH/ab/rkdp45_thermodynamics.dat" "$SCRATCH/ab/rkdp45_baseline.dat"
```
Expected: both comparisons print `OVERALL max|rel|` below `1e-10` (ideally `0.000e+00`). If any column shows a larger deviation, a delta was transcribed wrong — do not proceed; re-check Steps 5-8 against master (`git show master:source/thermodynamics_module.cpp`).

- [ ] **Step 11: Run the evolver scenario suite**

Run: `cd /Users/au192734/Projects/class_claude && test/scenarios/compare_evolver.sh ./class`
Expected: same result as PR #362 — only the pre-existing `RUNFAIL(ndf15) type3_scf_veta` ("Failure in sp_ludcmp. Possibly singular matrix!"), which also fails on master. No new failures.

- [ ] **Step 12: Commit**

```bash
cd /Users/au192734/Projects/class_claude
git add source/thermodynamics_module.h source/thermodynamics_module.cpp
git commit -m "$(cat <<'EOF'
Collapse RECFAST derivative modes onto master's phase-free RHS

Revert thermodynamics_derivs_with_recfast_member to master's state-based body
(the x_H>trigger and x_He<=trigger2 tests already yield the correct
frozen/evolved behavior in every phase), delete RecfastHydrogenMode /
RecfastHeliumMode, and route both evolver-driven phases through a single
RecfastPhase-dispatched derivs shim and output shim carried in the workspace.
Drop the dead zero-sample helium path. Output is unchanged at default precision.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01Uzbwf8RtwEDNs9vJNE4dDM
EOF
)"
```

---

### Task 3: Tolerance convergence study and value selection

Pick `tol_thermo_integration` empirically against a `1e-9` reference, set it in `precision.h`, and validate end-to-end.

**Files:**
- Create: `<scratch>/tol/tol_<val>.ini` for `val ∈ {1e-4,1e-5,1e-6,1e-7,1e-8,1e-9}`
- Create: `<scratch>/tol/tt_chosen.ini`, `<scratch>/tol/tt_master.ini` (CMB-TT sanity)
- Modify: `include/precision.h` (the `tol_thermo_integration` default)

**Interfaces:**
- Consumes: `compare_thermo.py` from Task 1; the Task 2 binary.
- Produces: the chosen `tol_thermo_integration` value and a sweep table for the PR body.

- [ ] **Step 1: Write the sweep inis (isolate thermo cost: no Cl output)**

```bash
SCRATCH=/private/tmp/claude-497981304/-Users-au192734-Projects-class-claude/09d6c8ce-241a-4b7e-8990-17a794a6aa03/scratchpad
mkdir -p "$SCRATCH/tol"
for v in 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9; do
  tag=${v/./p}
  cat > "$SCRATCH/tol/tol_$tag.ini" <<EOF
recombination = RECFAST
write thermodynamics = yes
tol_thermo_integration = $v
root = $SCRATCH/tol/tol_${tag}_
EOF
done
```

- [ ] **Step 2: Run the sweep and record wall time**

```bash
cd /Users/au192734/Projects/class_claude
SCRATCH=/private/tmp/claude-497981304/-Users-au192734-Projects-class-claude/09d6c8ce-241a-4b7e-8990-17a794a6aa03/scratchpad
for v in 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9; do
  tag=${v/./p}
  echo "== tol=$v =="; ( time ./class "$SCRATCH/tol/tol_$tag.ini" ) 2>&1 | grep -E "real|Computing thermodynamics"
done
```
Expected: six runs complete; note the `real` time for each (thermo dominates since there is no perturbation output).

- [ ] **Step 3: Compare each candidate to the `1e-9` reference**

```bash
SCRATCH=/private/tmp/claude-497981304/-Users-au192734-Projects-class-claude/09d6c8ce-241a-4b7e-8990-17a794a6aa03/scratchpad
REF="$SCRATCH/tol/tol_1pe-9_thermodynamics.dat"
for v in 1e-4 1e-5 1e-6 1e-7 1e-8; do
  tag=${v/./p}
  echo "== tol=$v vs 1e-9 =="
  python3 "$SCRATCH/compare_thermo.py" "$SCRATCH/tol/tol_${tag}_thermodynamics.dat" "$REF" | grep -E "\[800<z<3000\]|OVERALL"
done
```
Expected: the `[800<z<3000] x_e`/`Tb` deviations shrink monotonically as tol tightens. Record them in a table.

- [ ] **Step 4: Apply the decision criterion and pick the value**

Choose the **loosest** tol whose windowed `x_e` and `Tb` deviation from the `1e-9` reference is `<= 1e-5`. (Confirm with a TT check in Step 6 that this also gives `C_ell^TT` within `1e-4`.) Prefer the looser of two converged values unless the runtime difference is negligible. Note the decision explicitly, e.g. "1e-6 gives x_e 4e-7 / Tb 9e-8 → keep 1e-6" or "1e-5 already meets 1e-5 → loosen to 1e-5".

- [ ] **Step 5: Set the chosen value in precision.h**

If the chosen value differs from the current `1e-6`, edit `include/precision.h`:
```cpp
  double tol_thermo_integration = <chosen value>;
```
Leave the explanatory comment (already updated by PR #362) intact. If the chosen value **is** `1e-6`, make no code change here and record in the PR that the sweep confirms `1e-6`.

- [ ] **Step 6: CMB-TT observable check at the chosen tol (on-branch, isolates the tol effect)**

The decision criterion's observable half is `C_ell^TT` deviation `<= 1e-4` versus the `1e-9`
reference — measured on this branch so it isolates the tol's impact rather than conflating it
with PR #362's structural change. (The refactor's fidelity to master is already covered by
Task 2 Step 10 + PR #362's prior master validation, so no master rebuild is needed here.)

```bash
cd /Users/au192734/Projects/class_claude
SCRATCH=/private/tmp/claude-497981304/-Users-au192734-Projects-class-claude/09d6c8ce-241a-4b7e-8990-17a794a6aa03/scratchpad
CH=<chosen tol>            # e.g. 1e-6
for pair in "chosen:$CH" "ref:1e-9"; do
  name=${pair%%:*}; v=${pair##*:}
  cat > "$SCRATCH/tol/tt_$name.ini" <<EOF
recombination = RECFAST
output = tCl
tol_thermo_integration = $v
root = $SCRATCH/tol/tt_${name}_
EOF
  ./class "$SCRATCH/tol/tt_$name.ini"
done
python3 - <<'PY'
import numpy as np, glob
g = "/private/tmp/claude-497981304/-Users-au192734-Projects-class-claude/09d6c8ce-241a-4b7e-8990-17a794a6aa03/scratchpad/tol"
a = np.loadtxt(glob.glob(f"{g}/tt_chosen_*cl.dat")[0])
b = np.loadtxt(glob.glob(f"{g}/tt_ref_*cl.dat")[0])
n = min(len(a), len(b))
rel = np.abs(a[:n,1]-b[:n,1])/(np.abs(b[:n,1])+1e-300)   # col 1 = TT, no zero crossings
print(f"TT max|rel| (chosen vs 1e-9) = {rel.max():.3e}")
PY
```
Expected: `TT max|rel|` `<= 1e-4`. If larger, the chosen tol is too loose — tighten and repeat Steps 4-6.
Optional end-to-end confirmation (only if you want a master cross-check): build master in a
throwaway worktree (`git worktree add "$SCRATCH/master_wt" master` — never switch this working
tree) and compare TT; expect agreement within ~0.1%. Remove the worktree afterwards
(`git worktree remove "$SCRATCH/master_wt" --force`).

- [ ] **Step 7: Confirm the chosen value also converges under rkdp45**

```bash
SCRATCH=/private/tmp/claude-497981304/-Users-au192734-Projects-class-claude/09d6c8ce-241a-4b7e-8990-17a794a6aa03/scratchpad
CH=<chosen tol>            # e.g. 1e-6
cat > "$SCRATCH/tol/rkdp_chosen.ini" <<EOF
recombination = RECFAST
write thermodynamics = yes
evolver = 2
tol_thermo_integration = $CH
root = $SCRATCH/tol/rkdp_chosen_
EOF
./class "$SCRATCH/tol/rkdp_chosen.ini"
python3 "$SCRATCH/compare_thermo.py" "$SCRATCH/tol/rkdp_chosen_thermodynamics.dat" "$SCRATCH/tol/tol_1pe-9_thermodynamics.dat" | grep -E "\[800<z<3000\]"
```
Expected: rkdp45 at the chosen tol also meets the `<= 1e-5` windowed criterion vs the (ndf15) `1e-9` reference.

- [ ] **Step 8: Commit (only if precision.h changed)**

```bash
cd /Users/au192734/Projects/class_claude
git add include/precision.h
git commit -m "$(cat <<'EOF'
Set tol_thermo_integration from RECFAST convergence sweep

Sweep vs a 1e-9 reference: chose the loosest tol whose recombination-window
x_e/Tb deviation is <= 1e-5 (and CMB-TT within ~0.1% of master). Sweep table
recorded in the PR description.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01Uzbwf8RtwEDNs9vJNE4dDM
EOF
)"
```
If the sweep confirmed `1e-6`, skip the commit and instead paste the sweep table into the PR body.

- [ ] **Step 9: Update the PR description**

Add the sweep table (tol, `x_e`/`Tb` windowed deviation, wall time) and the TT sanity result to PR #362's body, and note that the RHS now matches master verbatim (restoring the `x_He0_trigger2` Heflag gate).
Run: `gh pr view 362 --json body -q .body` to read the current body before editing with `gh pr edit 362 --body-file <file>`.

---

## Self-Review

**1. Spec coverage:**
- Revert RHS to master verbatim → Task 2 Steps 5-7. ✓
- Delete `RecfastHydrogenMode`/`RecfastHeliumMode`/`legacy_trigger` → Task 2 Steps 1, 5-7. ✓
- `RecfastPhase` in workspace → Task 2 Step 2. ✓
- One derivs shim + one output shim → Task 2 Step 3. ✓
- Analytic phase keeps direct RHS call → unchanged (static `thermodynamics_derivs_with_recfast` at cpp:3514 retained; Task 2 Step 1 keeps its decl). ✓
- Drop dead zero-sample path + `class_test` guard → Task 2 Step 8. ✓
- Refactor inert vs branch tip → Task 1 + Task 2 Step 10. ✓
- Restore `x_He0_trigger2` Heflag gate (intended change) → Task 2 Step 6, validated inert at defaults by Step 10. ✓
- Tolerance sweep (ref 1e-9, candidates 1e-4..1e-8), criterion, set value → Task 3 Steps 1-5. ✓
- Sequenced verification (inert first, then tol) → Task 2 before Task 3. ✓
- ndf15 + rkdp45 both exercised → Task 1/Task 2 (both), Task 3 Step 7. ✓
- compare_evolver.sh unchanged → Task 2 Step 11. ✓
- CMB-TT ~0.1% sanity → Task 3 Step 6. ✓

**2. Placeholder scan:** `<chosen tol>` / `<val>` in Task 3 are deliberate — the value is the empirical output of Steps 3-4; every command that needs it either derives it (`CH=...`) or is gated on the decision. No `TBD`/`TODO`/"handle edge cases". ✓

**3. Type consistency:** `RecfastPhase{helium, full}` and `recfast_phase` used identically in header (Step 2), derivs shim, output shim, and driver (Steps 3, 8). Shim signatures match the header decls (Step 1) and the `generic_evolver` function-pointer types. `thermodynamics_derivs_with_recfast_member(double,double*,double*,void*)` is the single name used by the static forwarder, both shims, and the analytic call. ✓
