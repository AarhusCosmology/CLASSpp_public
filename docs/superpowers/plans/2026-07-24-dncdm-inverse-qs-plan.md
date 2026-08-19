# DNCDM Inverse Decays + Quantum Statistics Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend the decay-only `ncdm_decay_dr` model (massive WDM parent → bulk DR) to the full ν_H ↔ ν_l + φ system of arXiv:2011.01502 — momentum-resolved fermion + boson daughters, inverse decays, and (linear-order) quantum statistics — at background and synchronous-gauge perturbation level, without changing behaviour when the new flags are off.

**Architecture:** A pure, unit-tested conservative transition-network kernel (`tools/decay_transition_kernel`) discretizes the joint band integral once per RHS and derives all three species' derivatives (and the per-ℓ perturbation collision operator) from the SAME discrete transitions, so number/energy conservation and detailed balance hold on the grid to machine precision. One massless momentum-resolved daughter class (`DrPsdSpecies`, statistics=fermion|boson) plays both daughters. A composite (`DNCDMInvSpecies`) owns the kernel, a collision-owned parent `DNCDMSpecies`, and the two daughters, coupling them via `BackgroundDerivs` + `AddCouplingDerivs`. The existing `DNCDM_DR_Species` decay-only path is untouched.

**Tech Stack:** C++17, CLASSpp species-as-plugins architecture, ndf15 (default) + rkdp45 (`evolver=2`) evolvers, CMake/Makefile build, per-instance dot-syntax `.ini` input, classyref A/B harness for master comparison.

**Design doc (FIXED — do not redesign):** `docs/superpowers/specs/2026-07-24-dncdm-inverse-decay-qs-design.md`
**Research (read for physics/lessons):**
- `<scratchpad>/research-plugins.md` — architecture map with file:line refs (primary integration source)
- `<scratchpad>/research-thesis.md` — equations (eqs. 6.2–6.9), numerical lessons, validation ladder
- `<scratchpad>/research-branch.md` — prior-code lessons (moving-endpoint discontinuity, QS q/ε bug, dlnf/dlnq saga)
- `<scratchpad>/formalism-notes.md` — Λ kernel, kinematic bands, conservation identities
- `<scratchpad>` = `/private/tmp/claude-497981304/-Users-au192734-Projects-class-claude/e3a19e79-f3f2-4c25-a9cd-26f7a482f9e9/scratchpad`

## Global Constraints

Every task's requirements implicitly include this section. Values copied verbatim from the design doc and repo memory.

- **C++ only.** No `extern "C"`, no `#ifdef __cplusplus`. Only `hyrec/` is C.
- **Never hand-edit `cclassy.pxd`.** It is regenerated from headers by `generate_wrapper.py` at build. A new species needs zero `classy.pyx` edits unless it exposes a new named derived parameter (use `GetSpeciesParam`/`BaseSpecies::GetParam`).
- **Error-severity conventions.** `class_test` (→ `runtime_error` → CosmoComputationError → reject point) for any NUMERIC-value check a sampler could vary. `class_test_severe` (→ `invalid_argument` → CosmoSevereError → abort) ONLY for structural impossibilities (conflicting/missing keys, unsupported gauge/mode). A severe check must never gate on a sampler-varyable value.
- **Verify to ~0.1% tolerance, never bit-identical.** Handle Cl^TE zero-crossings (scale-relative metric, never blind max-rel-diff). The exception: the decay-only regression's background column layout may be asserted bit-identical (index counts), but physics values compare at ~0.1%.
- **Explicit-path git staging only.** Never `git add -A` (in-source CMake/Xcode artifacts get swept in). `git add -f` only for gitignored `.ini` test files if needed.
- **Comments state constraints, not narration.** Follow the doc-comment style already in `wdm_decay_product.h` (explain the load-bearing invariant, cite the spec/PR).
- **New test executables MUST be registered in BOTH** `CMakeLists.txt` (`add_executable` + the `foreach(_t IN ITEMS ...)` list at line 237) AND `Makefile` `TEST_TARGETS` (lines 4–24). CI's `make test` only builds what is in `TEST_TARGETS` (recurring trap).
- **Synchronous gauge only; no tensors with inverse decays; no ncdm fluid approximation for the coupled trio** — all `class_test`/`class_test_severe` at input time (type3 / dcdm_wdm precedent).
- **Never materialize `inf`/`inf()` under `-ffast-math`.** Compare in log form and return caps before exponentiating (see `WdmDecayProductSpecies::SigmaAt`).
- **Branch:** `dncdm-inverse-qs` (already checked out). NEVER `git checkout`/`switch`/`reset` — the working tree is shared with other agents.

---

## The Resolved Design Question: parent-PSD ownership when `inverse_decays=on`

**Question:** With inverse decays off, `DNCDMSpecies` evolves `ln f` per momentum bin with a diagonal decay RHS (`dncdm_species.cpp:400-414`) and a separately-integrated `dlnfdlnq` variable. With inverse decays on the parent needs `f`-variables and a kernel-supplied (non-diagonal, gain-carrying) RHS, and decay-only mode must remain **byte-for-byte the existing code path**.

**Resolution: a mode flag on `DNCDMSpecies` (`collision_owned_`), not composite-owned indices.** Rationale: the parent owns its PSD representation (grid, `ComputeMomenta`, `GetDlnf0Dlnq`, its perturbation hierarchy, the `RescaledPerturbations`/`GetRescalingFactor` underflow machinery). The composite cannot cleanly take over `ComputeBackground`/`ComputeMomenta`. So the parent keeps ownership of its integration + background slots in BOTH modes; only the *integration variable* (ln f → f) and the *RHS source* change. The composite supplies the RHS by writing into the parent's `f`-slots in its `BackgroundDerivs` — exactly the pattern `DCDM_WDM_Species::BackgroundDerivs` already uses to write `dy[wdm_->bi_f_index()+i] = J_i` (`dcdm_wdm_species.cpp:49-64`).

Concretely, add `bool collision_owned_ = false;` + `void SetCollisionOwned(bool v) { collision_owned_ = v; }` to `DNCDMSpecies`, set by the `DNCDMInvSpecies` factory immediately after constructing the parent child (before any index registration). Every change below is gated on `collision_owned_`; when it is false, the method body is the current code verbatim.

| Method | `collision_owned_ == false` (decay-only, UNCHANGED) | `collision_owned_ == true` (inverse) |
|---|---|---|
| `RegisterBackgroundIndices` | number,rho,p,pseudo_p + `lnf`,`dlnfdlnq`,`dlnfdlnq_sep` columns (unchanged) | **identical** (keeps the same bg columns so perturbations read them unchanged; `dlnfdlnq_sep` column is written 0) |
| `RegisterIntegrationIndices` | `index_bi_lnf_decay_dr1_` (q) + `index_bi_dlnfdlnq_separate_decay_` (q) | register `index_bi_f_parent_` (q) ONLY (drops the ln f + separate-dlnfdlnq variables) |
| `SetBackgroundInitialConditions` | seed `lnf = ln(w_/dq_)`, seed separate `dlnfdlnq = -q e^q/(e^q+1)` | seed `f = w_[i]/dq_[i]` into `index_bi_f_parent_` |
| `ComputeBackground` | spline ln f → dlnf/dlnq; fill bg columns (unchanged) | read `f` from `index_bi_f_parent_`, floor at `kFParentFloor`, take `ln f` → spline vs ln q → dlnf/dlnq; write `index_bg_lnf_decay_dr1_ = ln f`, `index_bg_dlnfdlnq_decay_`, `index_bg_dlnfdlnq_sep_ = 0`; then `ComputeMomenta` (unchanged) |
| `BackgroundDerivs` | write diagonal decay term (unchanged) | **no-op** (the composite assigns `dy[index_bi_f_parent_+i]` from the kernel) |
| perturbations (`PerturbDerivs`, `StressEnergy`, `GetDlnf0Dlnq`, `RescaledPerturbations`) | UNCHANGED in both modes — they read `index_bg_lnf_decay_dr1_` / `index_bg_dlnfdlnq_decay_`, which `ComputeBackground` fills identically in both modes |

New accessor for the composite: `int bi_f_parent_index() const { return index_bi_f_parent_; }`.

This keeps decay-only byte-for-byte (no `collision_owned_` branch is ever taken in a `DNCDM_DR_Species` run), gives the parent `f`-variables + kernel RHS when inverse is on, routes the parent's decay LOSS through the kernel (never double-counted with the diagonal ln f term, which is not registered in collision-owned mode), and leaves the entire perturbation surface untouched because both modes publish the same `index_bg_lnf_decay_dr1_`/`index_bg_dlnfdlnq_decay_` columns.

---

## File Structure

**Create:**
- `tools/decay_transition_kernel.h` / `.cpp` — the pure conservative transition-network kernel (grids, node placement, gather/scatter weights, Λ evaluation, per-ℓ angular operator). No CLASS module plumbing.
- `tools/decay_kernel_test.cpp` — kernel unit tests (Phases 1 + 4 operator identities).
- `species/dr_psd_species.h` / `.cpp` — `DrPsdSpecies` (kTypeName `dr_psd`): one massless momentum-resolved PSD, statistics=fermion|boson.
- `species/dr_psd_test.cpp` — standalone daughter construction/grid/free-streaming/switch-copy tests.
- `species/dncdm_inv_species.h` / `.cpp` — `DNCDMInvSpecies` composite (parent + fermion + boson, kernel-owning).
- `species/dncdm_inv_test.cpp` — composite background plumbing + conservation + guards + switch-copy tests.
- `python/tests/dncdm_inv_ab.py` — classyref A/B comparison driver (Phases 3–5; not a make target).

**Modify:**
- `species/dncdm_species.h` / `.cpp` — add `collision_owned_` mode (see design-question table).
- `species/all_species.h` — `#include "dr_psd_species.h"`, `#include "dncdm_inv_species.h"`; add one factory row for `DrPsdSpecies`. The `ncdm_decay_dr` row's function is retargeted to a dispatcher (below).
- `species/dncdm_dr_species.cpp` — `CreateAll` becomes the per-instance dispatcher: `inverse_decays=no` → build `DNCDM_DR_Species` (existing code); `inverse_decays=yes` → hand off to `DNCDMInvSpecies::Create`.
- `species/ncdm_family.cpp` — leave `kFluidApproximationConsumers[]` as is (DNCDM stays; `dr_psd` is deliberately excluded — no fluid approximation).
- `species/species_collection.cpp` — extend the `has_ncdm_` detection (`:54-56`) to also catch `DNCDMInvSpecies` (or make it inherit a type already caught).
- `include/precision.h` — add `int l_max_dncdm_col` near `l_max_dr_col` (`:319`).
- `source/input_module.cpp` — read `l_max_dncdm_col` next to `l_max_dr_col` (`:2542`).
- `CMakeLists.txt` — add `tools/decay_transition_kernel.cpp` to `CLASS_TOOLS_FILES` (`:108`); add the three species `.cpp` to `CLASS_SPECIES_FILES` (`:73`); add three `add_executable` + `foreach` entries (`:217-237`).
- `Makefile` — add three `TEST_TARGETS` (`:4-24`).
- `species/species_type_name_test.cpp` — pin `DrPsdSpecies::kTypeName` (and confirm `ncdm_decay_dr` unchanged).

---

# PHASE 1 — The transition kernel (`tools/decay_transition_kernel`)

Pure component, no CLASS plumbing. TDD: tests first. This is the riskiest phase (weight algebra); conservation is proven by machine-precision unit tests before any species touches it.

## Task 1.1: Kernel API header + skeleton

**Files:**
- Create: `tools/decay_transition_kernel.h`
- Create: `tools/decay_transition_kernel.cpp` (skeleton that links)
- Modify: `CMakeLists.txt:108` (add `tools/decay_transition_kernel.cpp` to `CLASS_TOOLS_FILES`)

**Interfaces:**
- Produces (consumed by Tasks 1.2+, Phase 3, Phase 4): the `DecayTransitionKernel` class, `GridView`, `Statistics`, `DecayTransitionKernel::Config`, `DecayTransitionKernel::Moments`.

- [ ] **Step 1: Write the header.** Create `tools/decay_transition_kernel.h`:

```cpp
#pragma once
#include <vector>

// Non-owning view of a STATIC momentum grid. The kernel never owns these; the
// species owns the buffers (parent reuses its DNCDM grid, daughters own their
// log-trapezoid grids). q ascending, size n; dq = trapezoid cell widths.
struct GridView {
  const double* q  = nullptr;
  const double* dq = nullptr;
  int n            = 0;
};

// Statistics sign for the (1±f) factors: Fermion -> Pauli blocking (1-f);
// Boson -> Bose enhancement (1+f).
enum class Statistics { Fermion, Boson };

/**
 * Conservative discrete transition network for ν_H ↔ ν_l + φ (massless
 * daughters). Discretizes the joint band integral ONCE per RHS evaluation and
 * derives all three background derivatives AND the per-ℓ perturbation collision
 * operator from the SAME discrete transitions, so number+energy conservation and
 * detailed balance hold on the grid to machine precision for any resolution.
 * Physics: research-thesis.md §1.2 (eqs. 6.2–6.9), formalism-notes.md, design §4.
 *
 * Continuum band factor (design §4.1), with the cubic f_H f_l f_φ cancelling:
 *   Λ(q1,q2,q3) = f_l(q2) f_φ(q3) − f_H(q1) + f_H(q1)[f_l(q2) − f_φ(q3)]
 *                 └── inv ──┘   └dec┘   └──────── qs (linear order) ────────┘
 * Kinematics (massless daughters): ε1 = √(q1²+a²m²); daughter band
 * q2 ∈ [(ε1−q1)/2, (ε1+q1)/2]; partner q3 = ε1 − q2 (so q2*+q3* = ε1 exactly,
 * making the two-bin linear deposit conserve number AND energy per transition).
 */
class DecayTransitionKernel {
 public:
  struct Config {
    int n_gauss             = 4;      // Gauss–Legendre nodes per parent band (S ≈ 4–8)
    bool inverse_decays     = true;   // include the l+φ→H repopulation (f_l f_φ) term
    bool quantum_statistics = false;  // include the (1±f) linear-order term f_H(f_l−f_φ)
  };
  struct Moments { double N_H, N_l, N_phi, E_H, E_l, E_phi; };

  DecayTransitionKernel(GridView parent, GridView fermion, GridView boson,
                        Statistics fermion_stat, Statistics boson_stat, Config cfg);

  // Background RHS. Zeroes df_* then accumulates. K = a² m Γ (design §4.1).
  // f_*/df_* sized to the matching grid's n. Pure; hot-path scratch preallocated.
  void ComputeBackgroundDerivs(double a, double m, double Gamma,
                               const double* f_H, const double* f_l, const double* f_phi,
                               double* df_H, double* df_l, double* df_phi) const;

  // Diagnostics: N_i = factor-free Σ dq q² df_i, E_i = Σ dq q² ε_i df_i.
  // (factor_/deg spin weights are applied by the caller; the kernel returns the
  //  bare grid moments so the conservation test can weight them.)
  Moments ComputeMoments(double a, double m,
                         const double* df_H, const double* df_l, const double* df_phi) const;

  // Precompute the transition geometry (node momenta q2*/q3*, deposit bins+split
  // weights on both daughter grids, gathered background f, Λ_s, and the angle
  // cosines cosα*/cosβ*/cosγ*) for the current (a,m,Γ,f) into internal scratch.
  // Background derivs call this internally; the perturbation module calls it once
  // per RHS so ApplyPerturbationOperator can run for every ℓ without recomputing
  // kinematics.
  void PrepareTransitions(double a, double m, double Gamma,
                          const double* f_H, const double* f_l, const double* f_phi) const;

  // Per-ℓ perturbation collision operator on the SAME network. Gathers Ψ at the
  // nodes with P_ℓ(cos·*) (recurrence at exact node kinematics — never interpolate
  // P_ℓ), scatters the collision rate into dΨ_* with the transpose weights.
  // Requires a prior PrepareTransitions for the current (a,m). ℓ=0 → number/energy
  // identity of the perturbed operator; ℓ=1 → momentum identity.
  void ApplyPerturbationOperator(int l,
                                 const double* Psi_H, const double* Psi_l, const double* Psi_phi,
                                 double* dPsi_H, double* dPsi_l, double* dPsi_phi) const;

  // ── Test introspection (do not use on the hot path) ──────────────────────
  int n_nodes_per_parent() const { return cfg_.n_gauss; }
  // Node kinematics for parent bin i, node s: fills q2*, q3*, cosα*, cosβ*, cosγ*.
  void NodeKinematics(int i, int s, double a, double m,
                      double& q2, double& q3, double& cos_a, double& cos_b, double& cos_g) const;

 private:
  GridView parent_, fermion_, boson_;
  Statistics fermion_stat_, boson_stat_;
  Config cfg_;
  std::vector<double> gl_nodes_, gl_weights_;  // Gauss–Legendre on [-1,1], size n_gauss
  // mutable scratch buffers filled by PrepareTransitions (sized in the ctor).
  mutable std::vector<double> q2_star_, q3_star_, lambda_gather_ /*...*/;
  // Helper: two-bin linear split of q* on a GridView -> (j, weight λ into j, 1-λ into j+1).
  static void TwoBinSplit(const GridView& g, double q_star, int& j, double& lambda);
};
```

- [ ] **Step 2: Write the `.cpp` skeleton** with the constructor (compute Gauss–Legendre nodes/weights on `[-1,1]` for `n_gauss`; size scratch to `parent_.n * n_gauss`), `TwoBinSplit` (binary search for the bracketing pair `j,j+1` with `g.q[j] ≤ q_star ≤ g.q[j+1]`; `lambda = (g.q[j+1]-q_star)/(g.q[j+1]-g.q[j])`; clamp to the edge bins when `q_star` is off-grid, depositing fully into the nearest edge bin — mirrors `GaussWeights`' edge-clamp so nothing leaks), and empty `ComputeBackgroundDerivs`/`PrepareTransitions`/`ApplyPerturbationOperator`/`ComputeMoments`/`NodeKinematics` bodies that compile.

- [ ] **Step 3: Register in CMake.** In `CMakeLists.txt:108-121` add `tools/decay_transition_kernel.cpp` to `CLASS_TOOLS_FILES` (alphabetical: after `tools/common.cpp`).

- [ ] **Step 4: Verify it builds.**
Run: `cmake --build build/cmake --target classpp --parallel 2>&1 | tail -5`
Expected: builds with no errors (skeleton links; unused-parameter warnings acceptable).

- [ ] **Step 5: Commit.**
```bash
git add tools/decay_transition_kernel.h tools/decay_transition_kernel.cpp CMakeLists.txt
git commit -m "kernel: DecayTransitionKernel API skeleton + CMake registration"
```

## Task 1.2: Background conservation + detailed balance (TDD)

**Files:**
- Create: `tools/decay_kernel_test.cpp`
- Modify: `CMakeLists.txt:217-244`, `Makefile:4-24` (register `test-decay-kernel`)
- Modify: `tools/decay_transition_kernel.cpp` (implement `ComputeBackgroundDerivs`, `PrepareTransitions`, `ComputeMoments`, `NodeKinematics`)

**Interfaces:**
- Consumes: the Task 1.1 API.
- Produces: a machine-precision-conserving background network trusted by Phase 3.

- [ ] **Step 1: Write the failing test file** `tools/decay_kernel_test.cpp` (single `assert`-based `main`, mirroring `dcdm_wdm_test.cpp` structure). Build three log-spaced grids (`q_min=1e-3,q_max=1e2,N=64` parent; daughters same), fill random positive `f_H,f_l,f_phi` (deterministic seed) and random `(a,m,Gamma)`. Assertions:

```cpp
// MOMENT CONVENTION (Fable review — load-bearing): ComputeMoments returns BARE
// per-dof grid moments N_i = Σ dq q² df_i, E_i = Σ dq q² ε_i df_i. The ×2 boson
// dof ratio (g_H·g_l/g_φ) is folded into df_phi (paper eq. 4.14 has a 2π
// prefactor where the others have 4π). With g_H = g_l = 2, g_phi = 1, one decay
// H → l + φ therefore gives bare moments N_l = −N_H and N_phi = −2·N_H, and the
// paper's App. C identities reduce to EXACTLY:
//   N_H + N_l == 0            (fermion-leg number)
//   2·N_H + N_phi == 0        (boson-leg number)
//   2·E_H + 2·E_l + E_phi == 0  (g-weighted energy; per transition
//                                2ε1 − 2q2* − 2q3* == 0 since q2*+q3* = ε1)
{
  DecayTransitionKernel::Config cfg; cfg.inverse_decays = true; cfg.quantum_statistics = true;
  DecayTransitionKernel K(parent_view, fermion_view, boson_view,
                          Statistics::Fermion, Statistics::Boson, cfg);
  K.ComputeBackgroundDerivs(a, m, Gamma, fH.data(), fl.data(), fphi.data(),
                            dfH.data(), dfl.data(), dfphi.data());
  auto M = K.ComputeMoments(a, m, dfH.data(), dfl.data(), dfphi.data());
  assert(std::fabs(M.N_H + M.N_l) < 1e-12 * Scale(M.N_H, M.N_l));       // number, fermion leg
  assert(std::fabs(2.*M.N_H + M.N_phi) < 1e-12 * Scale(2.*M.N_H, M.N_phi)); // number, boson leg
  assert(std::fabs(2.*M.E_H + 2.*M.E_l + M.E_phi)
         < 1e-12 * Scale3(2.*M.E_H, 2.*M.E_l, M.E_phi));                 // energy, g-weighted
  assert(std::fabs(M.N_H) > 0.);  // guard against a trivially-zero kernel "passing"
}
// Equilibrium residual (Fable review — replaces the original 1e-11 detailed-balance
// assertion, which is UNACHIEVABLE with linear gather interpolation): for FD
// parent+fermion, BE boson with common T and μ_H = μ_l + μ_φ, the continuum Λ ≡ 0
// on the manifold ε1 = q2* + q3*; on the grid the gathered f_l, f_φ are linear
// interpolants, so Λ_s = O(Δq²) — a SMALL residual that is tangent to the
// conservation manifold (N/E identities above must STILL hold to 1e-12 here).
// Test it as a convergence ratio: doubling the daughter grids must shrink the
// equilibrium residual by ~4 (second order). Accept ratio > 2.5.
{
  // f_H(q)=1/(exp((ε1(q)-muH)/T)+1) with ε1=√(q²+a²m²); f_l=1/(exp((q-mul)/T)+1);
  // f_phi=1/(exp((q-muphi)/T)-1) with muphi < q_min (keep BE regular); muH = mul + muphi.
  double r64  = MaxAbsDf(K_64bins,  ...);   // ‖df‖_∞ / typical |f·K-rate| at N_daughter=64
  double r128 = MaxAbsDf(K_128bins, ...);   // same at N_daughter=128
  assert(r128 < r64 / 2.5);                 // O(h²) shrinkage
  // and N/E conservation asserts (as above) at BOTH resolutions — exact regardless.
}
// Mode reduction: cfg.inverse_decays=false, cfg.quantum_statistics=false ->
// df_H bin i == -K/ε1 · f_H(q1_i) (pure diagonal decay, matches dncdm_species.cpp:410
// in the continuum: d f_H/dτ = -a²mΓ/ε1 · f_H); daughters gain the decay integral only.
```

Provide `Scale`/`Scale3` = `max(|terms|, 1e-300)` to avoid 0/0. Add a `main` that runs all blocks and prints `"decay kernel test passed"`.

- [ ] **Step 2: Register the test executable in BOTH build files.**
  - `CMakeLists.txt:236` add `add_executable(test-decay-kernel tools/decay_kernel_test.cpp)` and append `test-decay-kernel` to the `foreach(_t IN ITEMS ...)` list at `:237`.
  - `Makefile:4-24` append `\` + `\ttest-decay-kernel` to `TEST_TARGETS`.

- [ ] **Step 3: Run to verify it fails.**
Run: `cmake --build build/cmake --target test-decay-kernel --parallel && ./build/cmake/test-decay-kernel`
Expected: FAIL (empty kernel bodies → assertions fire / all-zero moments trivially "pass" only for the decay-off block; ensure the conservation block fails because df are all zero and N_H≠0). If a trivial all-zero passes, tighten by also asserting `M.N_H != 0` when `Gamma>0` and inverse/decay on.

- [ ] **Step 4: Implement `PrepareTransitions` + `ComputeBackgroundDerivs`.** Per design §4.2 / thesis §1.2:
  - For each parent bin `i`: `ε1 = √(q1_i²+a²m²)`; band `[qlo,qhi] = [(ε1−q1_i)/2,(ε1+q1_i)/2]`; place `S` Gauss–Legendre nodes `q2*_s = 0.5(qhi−qlo)·(gl+1)+qlo`, node weight `wn_s = gl_weights[s]·0.5(qhi−qlo)`; partner `q3*_s = ε1 − q2*_s`.
  - Gather `f_l(q2*_s)`, `f_phi(q3*_s)` by the two-bin linear interpolation (`TwoBinSplit` on the fermion/boson grids).
  - `Λ_s = f_l·f_phi − f_H(q1_i)` `+ (qs ? f_H(q1_i)·(f_l − f_phi) : 0)` `− (inverse ? 0 : f_l·f_phi)` (i.e. drop the `f_l f_phi` inverse term when inverse off; drop qs term when qs off).
  - **Exact number-space bookkeeping (Fable review — this REPLACES the earlier
    hand-wavy `parent_measure` wording; get this exactly right or conservation
    fails):** define per transition (parent bin i, node s) the signed number rate
    `δN_s ≡ K · dq1[i] · (q1_i/ε1) · wn_s · Λ_s`   with `K = a²·m·Gamma`
    (δN_s < 0 for net decay). Then:
      - parent:  `df_H[i] += δN_s / (q1_i² · dq1[i])`
                 (algebraically = `(K/(ε1·q1_i))·wn_s·Λ_s`, the direct continuum form)
      - fermion: `df_l[j]   += −λ_s     · δN_s / (q2_j²   · dq2[j])`
                 `df_l[j+1] += −(1−λ_s) · δN_s / (q2_{j+1}² · dq2[j+1])`
      - boson:   `df_phi[k]   += −2·μ_s     · δN_s / (q3_k²   · dq3[k])`
                 `df_phi[k+1] += −2·(1−μ_s) · δN_s / (q3_{k+1}² · dq3[k+1])`
        (×2 = g_H·g_l/g_φ dof ratio, folded into df_phi as in paper eq. 4.14).
    The `1/(q² dq)` factors convert a number deposit back into a PSD derivative at
    the deposit bin; with `λ q2_j + (1−λ) q2_{j+1} = q2*` the deposits conserve
    number AND energy per transition exactly (q2* + q3* = ε1). Verify against the
    continuum: summing fermion deposits over (i,s) reproduces
    `−K/q2² ∫dq1 (q1/ε1) Λ` as the grid refines.
  - **Off-grid deposits:** when `q*` falls outside a daughter grid, clamp the full
    deposit into the edge bin (number-preserving) and accumulate the resulting
    energy mismatch `(q_edge − q*)·δN` into a public diagnostic
    (`double clamped_energy_residual()` on the kernel, reset by PrepareTransitions).
    The conservation unit test must choose grids that cover the bands (no clamping,
    residual == 0); Phase 3 wires the diagnostic into a background column and the
    daughter grid defaults must cover the band support of the run
    (q_max_daughter ≳ 0.6·a_end·m in T_dr units for late-decay configs — document
    as an input-guard `class_test` warning path, not a hard failure).
  - **Transpose consistency (the load-bearing invariant):** the same `TwoBinSplit` weight `λ` used to GATHER `f_l(q2*)` must be used to SCATTER into `df_l`. Store `λ` per node. This is what makes conservation exact and replaces the A+B·q calibration entirely (design §2 item 3, §4.2). Store all node data (q2*,q3*,λ_fermion,j_fermion,λ_boson,j_boson,f_l,f_phi,Λ_s,cosα*,cosβ*,cosγ*) in scratch for the perturbation operator.
  - `NodeKinematics`: cos angles per formalism-notes: `cosα* = (2ε1·q2 − a²m²)/(2 q1 q2)`, `cosβ* = (2ε1·q3 − a²m²)/(2 q1 q3)`, `cosγ* = 1 − a²m²/(2 q2 q3)`.
  - `ComputeMoments`: `N_i = Σ_bins dq[b]·q[b]²·df_i[b]`, `E_i = Σ dq[b]·q[b]²·ε_i[b]·df_i[b]` (ε_daughter = q; ε_H = √(q²+a²m²)).

- [ ] **Step 5: Run to verify it passes.**
Run: `./build/cmake/test-decay-kernel`
Expected: PASS, prints `decay kernel test passed`. If the N/E conservation identities fail at `1e-12` (they must hold for ARBITRARY f arrays, equilibrium or not), the gather/scatter weights are NOT transpose-consistent or a `1/(q² dq)` conversion is wrong — fix the algebra, never loosen the tolerance (per research-thesis lesson 5 and design §4.2). The equilibrium-residual block is the only one with a soft (convergence-ratio) criterion.

- [ ] **Step 6: Commit.**
```bash
git add tools/decay_kernel_test.cpp tools/decay_transition_kernel.cpp CMakeLists.txt Makefile
git commit -m "kernel: background transition network — machine-precision N/E conservation + detailed balance"
```

## Task 1.3: RHS continuity across band edges (TDD)

**Files:** Modify `tools/decay_kernel_test.cpp`.

- [ ] **Step 1: Write the failing continuity test.** Sweep `a` finely so a band edge `(ε1±q1)/2` crosses a daughter grid node; assert the max jump in any `df` between adjacent `a` samples is `O(Δa)` (C⁰), NOT a step:
```cpp
// Sweep a across a node crossing; the RHS must be continuous (kinks only C0 as a
// Gauss node crosses a grid point — design §4.2; the 2020–2024 branch's #1 bug was
// a DISCONTINUOUS jump here, research-branch.md §8/§9).
double prev = kernel_df_at(a0), max_ratio = 0;
for (double a = a0; a < a1; a += da) {
  double cur = kernel_df_at(a);
  max_ratio = std::max(max_ratio, std::fabs(cur - prev) / (da * scale));
  prev = cur;
}
assert(max_ratio < 50.);  // bounded slope; a discontinuity would blow this up
```

- [ ] **Step 2: Run to verify.** `./build/cmake/test-decay-kernel` — the Gauss-node placement moves smoothly inside the band, so this should PASS immediately if Task 1.2 is correct. If it fails, the node placement is being spliced instead of scaled — revisit Step 1.2.4.

- [ ] **Step 3: Commit.**
```bash
git add tools/decay_kernel_test.cpp
git commit -m "kernel: RHS-continuity test across band-edge node crossings"
```

---

# PHASE 2 — The daughter species (`DrPsdSpecies`)

ONE massless momentum-resolved PSD class, statistics=fermion|boson parameter (thesis "sibling" lesson: code one daughter equation). Standalone it is a free-streaming massless PSD species; the composite injects the source.

## Task 2.1: `DrPsdSpecies` construction, grid, background plumbing (TDD)

**Files:**
- Create: `species/dr_psd_species.h` / `.cpp`, `species/dr_psd_test.cpp`
- Modify: `species/all_species.h`, `CMakeLists.txt`, `Makefile`, `species/species_type_name_test.cpp`

**Interfaces:**
- Produces (Phase 3/4 consume): `DrPsdSpecies` (inherits `NCDMBaseSpecies`), accessors `bi_f_index()`, `bg_f_index()`, `bg_dlnfdlnq_index()`, `Statistics statistics()`, `const std::vector<double>& q_bg()`, and the perturbation `PerturbLayout` (NCDM-family shape).

- [ ] **Step 1: Write the header** `species/dr_psd_species.h`. Inherit `NCDMBaseSpecies` with `DeferInit{}` (build the grid ourselves, like `WdmDecayProductSpecies`). kTypeName `dr_psd`. Members: `Statistics stat_`; log-trapezoid grid params `q_min_,q_max_,N_q_` (perturbation grid) + `q_min_bg_,q_max_bg_,N_q_bg_` (background grid, `N_bg ≥ N_pt`, design §6); `index_bi_f_` (bg-grid f, per bin), `index_bg_f_`, `index_bg_dlnfdlnq_`, `index_bg_number_`, `index_bg_rho_`, `index_bg_p_`, `index_bg_pseudo_p_`. Set `M_ = 0` (massless) and `factor_` via `SetDegAndFactor` with `deg = 1` (boson) or `2` (fermion) per the g-factor convention — but note the ×2 boson dof is folded into the kernel deposit, so keep `deg_` = the species' spin dof and be consistent with `ComputeMomenta`. Public overrides mirror `WdmDecayProductSpecies` (massless: `M_=0`): `RegisterBackgroundIndices`, `RegisterIntegrationIndices`, `SetBackgroundInitialConditions`, `ComputeBackground`, `Rho/P/PPrime`, `RhoDotOverRho`, `RegisterPerturbationIndices`, `PerturbDerivs`, `StressEnergy`, `ApplyInitialConditions`, `GetOmega0()==0` (default; the composite reserves the sector), decay-product neutrality overrides (`NeutrinoOmega0`/`NeffContribution`/`BackgroundAIni`/`CheckUltraRelativisticAtIc`/`IsUltraRelativisticAtIc` as in `wdm_decay_product.h:148-167`), `IsFreestreaming()==true`. `GetDlnf0Dlnq(iq,pvecback)` returns `pvecback[index_bg_dlnfdlnq_+iq]` (live spline value — the DNCDM seam).

```cpp
class DrPsdSpecies : public NCDMBaseSpecies {
 public:
  static constexpr std::string_view kTypeName = "dr_psd";
  static constexpr double kFFloor = 1e-100;   // positivity floor before ln f (gains dominate)
  static constexpr double kFSeed  = 1e-10;    // evolver-friendly nonzero start (wdm_decay_product.h:247)
  DrPsdSpecies(FileContent* pfc, const std::string& instance_name,
               const NcdmSettings& settings, const background* pba,
               const BackgroundModule* bgm, Statistics stat);
  Statistics statistics() const { return stat_; }
  int bi_f_index() const { return index_bi_f_; }
  int bg_f_index() const { return index_bg_f_; }
  int bg_dlnfdlnq_index() const { return index_bg_dlnfdlnq_; }
  const std::vector<double>& q_bg() const { return q_bg_; }
  const std::vector<double>& dq_bg() const { return dq_bg_; }
  struct Named { std::string key; std::unique_ptr<DrPsdSpecies> species; };
  static std::vector<::Named> CreateAll(const SpeciesBuildContext& ctx); // standalone dr_psd
  // ... overrides ...
 private:
  void BuildLogTrapezoidGrids();      // fills q_/dq_ (pt) and q_bg_/dq_bg_ (bg)
  std::vector<double> dq_bg_;
  Statistics stat_ = Statistics::Fermion;
  double q_min_=1e-2,q_max_=25.,q_min_bg_=1e-3,q_max_bg_=1e2;
  int N_q_=15,N_q_bg_=100;
  int index_bi_f_=-1,index_bg_f_=-1,index_bg_dlnfdlnq_=-1;
  int index_bg_number_=-1,index_bg_rho_bg_=-1 /*use base index_bg_rho_*/,index_bg_p_=-1,index_bg_pseudo_p_=-1;
};
```

Constructor reads dot-keys via `SpeciesInput`: `statistics` (`fermion`/`boson`, or passed by the composite), `dr_q_min/dr_q_max/dr_N_q` (+ `_pt` variants; design §5), optional initial abundance (`Omega_ini`/`T`/`ksi`, default zero). Build grids with `BuildLogTrapezoidGrids` (log-spaced trapezoid with explicit `q_min`, thesis §3 lesson 8: `x[i]=exp(log(qmin)+i·h)`, half-weight endpoints). Enforce `class_test(N_q_bg_ < N_q_, ...)` → `N_bg ≥ N_pt`.

- [ ] **Step 2: Write the failing test** `species/dr_psd_test.cpp`:
```cpp
// Construction + grid: sizes, log spacing, N_bg >= N_pt, massless M()==0.
// Background plumbing: RegisterBackgroundIndices / RegisterIntegrationIndices index
//   counts; SetBackgroundInitialConditions seeds kFSeed; populate one bg bin by hand,
//   ComputeBackground, assert Rho matches closed-form Σ factor/a⁴ dq q² ε f with ε=q.
// dlnf/dlnq: seed a pure power-law f ∝ q^p on the grid -> ComputeBackground -> assert
//   pvecback[bg_dlnfdlnq_ + i] ≈ p at interior bins (spline correctness).
// Standalone free-streaming: build a thermal-BE dr_psd with an initial abundance,
//   run a background-only Cosmology, assert rho redshifts ∝ a^-4 (radiation).
// Statistics guard: statistics key must be fermion|boson (class_test_severe otherwise).
```

- [ ] **Step 3: Register** `test-dr-psd` in `CMakeLists.txt` (`add_executable(test-dr-psd species/dr_psd_test.cpp)` + `foreach`) and `Makefile` `TEST_TARGETS`. Add `species/dr_psd_species.cpp` to `CLASS_SPECIES_FILES`. In `all_species.h` add `#include "dr_psd_species.h"` and a factory row `SpeciesFactoryEntry{DrPsdSpecies::kTypeName, &DrPsdSpecies::CreateAll}`. Pin the kTypeName in `species/species_type_name_test.cpp`.

- [ ] **Step 4: Run to verify it fails.**
Run: `cmake --build build/cmake --target test-dr-psd --parallel && ./build/cmake/test-dr-psd`
Expected: FAIL (unimplemented bodies).

- [ ] **Step 5: Implement `dr_psd_species.cpp`.** Background: `RegisterBackgroundIndices` claims number,rho,p,pseudo_p + `f`(N_bg),`dlnfdlnq`(N_bg); `RegisterIntegrationIndices` claims `index_bi_f_`(N_bg). `SetBackgroundInitialConditions` seeds `f = f0_thermal(q_bg) + kFSeed` (thermal FD/BE if an initial abundance is given, else `kFSeed`). `ComputeBackground`: subtract seed, floor at `kFFloor`, `w_bg_[i]=f·dq_bg_[i]`, spline `ln f` vs `ln q` → `index_bg_dlnfdlnq_`, `ComputeMomenta` (M_=0 ⇒ ε=q). Perturbations: NORMALIZED Ψ hierarchy identical in FORM to `DNCDMSpecies::PerturbDerivs` with `M_=0` (ε=q, so `qk_div_epsilon = k`), driven by `dlnf0_dlnq = pvecback[index_bg_dlnfdlnq_+iq]`; empty bins (`f<kFFloor`) start Ψ=0. `StressEnergy`: `w0 = f·dq_bg_` weights on the pt grid via the bridge (`InterpolateBgToPt`), or — for v1 — set `N_pt=N_bg` (single grid) so no bridge is needed; the plan permits `N_pt == N_bg` initially and flags the separate-grid bridge as a Phase-5 accuracy refinement (design §6: "the grid only affects resolution, not budgets"). `ApplyInitialConditions`: seed Ψ from the adiabatic contrast on populated bins, 0 on empty. `CopyPerturbationsAcrossSwitch`: inherit the base `NCDMBaseSpecies` shape-preserving copy (no fluid approximation).

- [ ] **Step 6: Run to verify it passes.**
Run: `./build/cmake/test-dr-psd && cmake --build build/cmake --target test-species-types --parallel && ./build/cmake/test-species-types`
Expected: both PASS.

- [ ] **Step 7: Commit.**
```bash
git add species/dr_psd_species.h species/dr_psd_species.cpp species/dr_psd_test.cpp species/all_species.h species/species_type_name_test.cpp CMakeLists.txt Makefile
git commit -m "dr_psd: massless momentum-resolved daughter (fermion|boson), standalone free-streaming"
```

## Task 2.2: `DrPsdSpecies` switch-copy regression (TDD)

**Files:** Modify `species/dr_psd_test.cpp`.

- [ ] **Step 1: Add a switch-copy test block** mirroring `dncdm_switch_copy_test.cpp:87-122`: build a standalone `DrPsdSpecies`, fabricate two `NCDMBaseSpecies::PerturbLayout` at different offsets, fill markers, `CopyPerturbationsAcrossSwitch`, assert every `(iq,l)` slot migrated and out-of-layout slots untouched. This guards #372 for the new per-bin layout.

- [ ] **Step 2: Run.** `./build/cmake/test-dr-psd` — Expected PASS (inherits the base copy).

- [ ] **Step 3: Commit.**
```bash
git add species/dr_psd_test.cpp
git commit -m "dr_psd: CopyPerturbationsAcrossSwitch regression (#372 class)"
```

---

# PHASE 3 — The composite (`DNCDMInvSpecies`), background

Parent `DNCDMSpecies` (collision-owned) + fermion `DrPsdSpecies` + boson `DrPsdSpecies`, kernel-driven `BackgroundDerivs`. Gates every validation step. See the resolved-design-question table above for the parent-ownership mechanism.

## Task 3.1: `collision_owned_` mode on `DNCDMSpecies` (TDD)

**Files:** Modify `species/dncdm_species.h`, `species/dncdm_species.cpp`; extend `species/dncdm_switch_copy_test.cpp` OR add to `species/dncdm_inv_test.cpp` (created in 3.2). Use the existing `test-dncdm-switch-copy` executable for the index-count assertions to avoid a premature new target.

**Interfaces:**
- Produces: `SetCollisionOwned(bool)`, `bi_f_parent_index()`, and the mode-gated background behaviour in the design-question table.

- [ ] **Step 1: Write a failing unit test** in `species/dncdm_switch_copy_test.cpp` (a new `{ }` block): construct a `DNCDMSpecies`, call `SetCollisionOwned(true)`, then `RegisterIntegrationIndices` and assert `index_bi == q_size()` (only the f-variable, not `2*q_size()`). Seed ICs, `ComputeBackground` on a synthetic `f`-state, assert `Rho` matches the closed-form quadrature and `pvecback[bg_lnf_index()+i] == ln(f_i)`. Assert a decay-only `DNCDMSpecies` (flag false) still registers `2*q_size()` (unchanged).

- [ ] **Step 2: Run to verify it fails.**
Run: `cmake --build build/cmake --target test-dncdm-switch-copy --parallel && ./build/cmake/test-dncdm-switch-copy`
Expected: FAIL (no `SetCollisionOwned`).

- [ ] **Step 3: Implement the mode.** Add `bool collision_owned_=false;` + `SetCollisionOwned`/`bi_f_parent_index` + `int index_bi_f_parent_=-1;` + `static constexpr double kFParentFloor=1e-100;` to the header. In `.cpp`, gate `RegisterIntegrationIndices`, `SetBackgroundInitialConditions`, `ComputeBackground`, `BackgroundDerivs` exactly per the design-question table. Every `else`/false branch is the CURRENT code verbatim (do not reformat it — keep the diff to the added `if (collision_owned_) {...} else {<existing>}`).

- [ ] **Step 4: Run to verify it passes.**
Run: `./build/cmake/test-dncdm-switch-copy`
Expected: PASS.

- [ ] **Step 5: Decay-only byte-for-byte check.** Build `class`, run a decay-only `ncdm_decay_dr` `.ini` (no `inverse_decays`) before and after this commit and confirm identical background output.
Run: `cmake --build build/cmake --target class --parallel && ./build/cmake/class <a decay-only dncdm .ini writing background>` then `diff` the background table against a pre-change run.
Expected: identical (no `collision_owned_` branch taken). If a `.ini` is not handy, use `python/tests/dncdm_inv_ab.py --decay-only` (Task 3.4) against `classyref`.

- [ ] **Step 6: Commit.**
```bash
git add species/dncdm_species.h species/dncdm_species.cpp species/dncdm_switch_copy_test.cpp
git commit -m "dncdm: collision-owned mode (f-variables + kernel RHS) gated; decay-only path unchanged"
```

## Task 3.2: `DNCDMInvSpecies` construction, factory dispatch, guards (TDD)

**Files:**
- Create: `species/dncdm_inv_species.h` / `.cpp`, `species/dncdm_inv_test.cpp`
- Modify: `species/dncdm_dr_species.cpp` (factory dispatch), `species/all_species.h`, `species/species_collection.cpp`, `CMakeLists.txt`, `Makefile`

**Interfaces:**
- Consumes: `DNCDMSpecies` (collision-owned), two `DrPsdSpecies`, `DecayTransitionKernel`.
- Produces: `DNCDMInvSpecies::Create(std::unique_ptr<DNCDMSpecies>, const SpeciesBuildContext&)` returning `Named`; the composite class following the `child_layouts` contract (children order: `kParent=0, kFermion=1, kBoson=2`).

- [ ] **Step 1: Write the header** `species/dncdm_inv_species.h`. `class DNCDMInvSpecies : public CompositeSpecies`. Members: `DNCDMSpecies* parent_`, `DrPsdSpecies* fermion_`, `DrPsdSpecies* boson_`, `std::unique_ptr<DecayTransitionKernel> kernel_`, scratch `df_H_,df_l_,df_phi_`. Enum `ChildIndex {kParent=0,kFermion=1,kBoson=2}` and typed layout views (mirror `dncdm_dr_species.h:26-34`). Override `BackgroundDerivs`, `AddCouplingDerivs`, `RegisterPerturbationIndices`, `GetShootingTargets`/`ComputeShootingGuess`/`ComputeShootingResidual`/`GetOmega0` (reuse the DNCDM combined-reserve logic), `WriteBackgroundColumnTitles`/`WriteBackgroundData` (parent + both daughters + conservation columns), `FillSources`, `WriteOutputColumns`. Declare `static Named Create(std::unique_ptr<DNCDMSpecies> parent, const SpeciesBuildContext& ctx);`.

- [ ] **Step 2: Implement `Create` + ctor.** `Create` reads (via `SpeciesInput` on the parent's instance name) `inverse_decays`, `quantum_statistics`, `dr_*` grid keys, initial daughter abundances. Guards (input-time):
```cpp
// yes/no read (repo convention, input_module.cpp:1829): first char 'y'/'Y'/'1'/'t'.
const bool inv = ParseYesNo(in.get<std::string>("inverse_decays"), false);
const bool qs  = ParseYesNo(in.get<std::string>("quantum_statistics"), false);
class_test(qs && !inv,
           "species '%s': quantum_statistics requires inverse_decays = yes "
           "(the (1±f) terms are only meaningful with both dec+inv present)",
           name.c_str());
```
Gauge guard (mirror `dcdm_wdm_species.cpp:144-150`): `class_test_severe` if `gauge` starts with "new". Tensor guard: `class_test_severe` on `modes` containing tensors when inverse on. Fluid-approx guard: `class_test` if `<inst>.fluid_approximation` or the global `ncdm_fluid_approximation` is set with these coupled species present (follow type3 precedent). Ctor: `SetCollisionOwned(true)` on the parent; build the two `DrPsdSpecies` (fermion, boson) from the instance's `dr_*` keys; construct `kernel_` with the three species' **background** grid views and `Config{n_gauss, inv, qs}`; push children in order `[parent, fermion, boson]`; size scratch.

- [ ] **Step 3: Implement `BackgroundDerivs`.** Mirror `dcdm_wdm_species.cpp:49-64`:
```cpp
void DNCDMInvSpecies::BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) {
  CompositeSpecies::BackgroundDerivs(tau, y, dy, pvecback);  // daughters' dilution; parent no-op (collision-owned)
  const double a = pvecback[bgm_->index_bg_a_];
  const double m = parent_->GetMass();
  const double Gamma = parent_->Gamma();
  // Gather the three background PSDs (parent f from its integration slots; daughters from theirs).
  kernel_->ComputeBackgroundDerivs(a, m, Gamma,
      &y[parent_->bi_f_parent_index()], &y[fermion_->bi_f_index()], &y[boson_->bi_f_index()],
      df_H_.data(), df_l_.data(), df_phi_.data());
  for (int i=0;i<parent_->q_size();++i)  dy[parent_->bi_f_parent_index()+i] = df_H_[i];
  for (int i=0;i<fermion_->q_bg().size();++i) dy[fermion_->bi_f_index()+i] += df_l_[i];   // += daughter also dilutes? massless: no dilution term (ε=q), so DrPsd BackgroundDerivs is a no-op -> use =
  for (int i=0;i<boson_->q_bg().size();++i)   dy[boson_->bi_f_index()+i]   += df_phi_[i];
}
```
Note: massless daughters have no background dilution term in comoving `f` (the redshift is in `ε=q` comoving), so `DrPsdSpecies::BackgroundDerivs` writes 0 for its f-slots and the composite `+=` the kernel source (or assign `=`; be consistent — assign `=` matching dcdm_wdm, and make `DrPsdSpecies` not register a BackgroundDerivs write).

- [ ] **Step 4: Factory dispatch.** In `species/dncdm_dr_species.cpp` `DNCDM_DR_Species::CreateAll`, after `auto dncdm_vec = DNCDMSpecies::CreateAll(ctx);`, branch per instance:
```cpp
for (auto& e : dncdm_vec) {
  SpeciesInput in(ctx.pfc, e.key);
  const bool inv = ParseYesNo(in.get<std::string>("inverse_decays"), false);
  if (!inv) {
    result.push_back({e.key, std::make_unique<DNCDM_DR_Species>(std::move(e.species), ctx.pba, ctx.bgm)}); // UNCHANGED
  } else {
    result.push_back(DNCDMInvSpecies::Create(std::move(e.species), ctx));
  }
}
```
Add `ParseYesNo` as a small shared helper (e.g. in `species/species_input.h`). Add `#include "dncdm_inv_species.h"` where needed. In `all_species.h`, add `#include "dncdm_inv_species.h"` (the factory row for `ncdm_decay_dr` stays pointed at `DNCDM_DR_Species::CreateAll`, which now dispatches). In `species_collection.cpp:54-56`, extend the `has_ncdm_` detection to also catch `DNCDMInvSpecies` (add `|| dynamic_cast<DNCDMInvSpecies*>(...)`).

- [ ] **Step 5: Write the failing test** `species/dncdm_inv_test.cpp` (register `test-dncdm-inv` in BOTH build files, add `species/dncdm_inv_species.cpp` to `CLASS_SPECIES_FILES`):
```cpp
// Guards: qs without inverse -> ThrowsComputation (class_test/runtime_error).
//         newtonian gauge with inverse -> Throws (class_test_severe/invalid_argument).
// Construction: DNCDMInvSpecies::Create builds parent(collision-owned)+2 daughters;
//   children order [parent, fermion, boson]; kernel n_gauss from precision/default.
// Background plumbing: RegisterBackgroundIndices/RegisterIntegrationIndices counts =
//   parent(q + bg cols) + 2×daughter(N_bg); IC seeding; ComputeBackground on a
//   synthetic state -> Rho sums children.
// Conservation columns on a synthetic state: assemble a random f-state, call
//   BackgroundDerivs, read the kernel Moments back -> g_H N_H + ½(g_l N_l+g_φ N_φ)=0
//   and energy identity to 1e-11 (the SAME invariant as the kernel test, now through
//   the composite plumbing).
// Switch-copy: composite CopyPerturbationsAcrossSwitch migrates all 3 children
//   (mirror dncdm_switch_copy_test.cpp:124-168, three child_layouts).
```

- [ ] **Step 6: Run to verify it fails, implement remaining bodies, run to pass.**
Run: `cmake --build build/cmake --target test-dncdm-inv --parallel && ./build/cmake/test-dncdm-inv`
Expected: FAIL → implement → PASS (`dncdm inv test passed`). Also run `test-decay-kernel`, `test-dr-psd`, `test-composite-layout`, `test-dncdm-switch-copy` to confirm no regressions.

- [ ] **Step 7: Commit.**
```bash
git add species/dncdm_inv_species.h species/dncdm_inv_species.cpp species/dncdm_inv_test.cpp species/dncdm_dr_species.cpp species/all_species.h species/species_collection.cpp species/species_input.h species/species_input.cpp CMakeLists.txt Makefile
git commit -m "dncdm_inv: kernel-owning composite (parent+fermion+boson), factory dispatch, guards, conservation columns"
```

## Task 3.3: Shooting closure for the inverse composite

**Files:** Modify `species/dncdm_inv_species.cpp` (shooter hooks), `species/dncdm_inv_test.cpp`.

- [ ] **Step 1: Implement the shooter hooks** by reusing the DNCDM combined-reserve logic (`dncdm_dr_species.cpp:295-364`): `GetOmega0()` returns the pinned `Omega_dncdmdr`; `GetShootingTargets`/`ComputeShootingGuess`/`ComputeShootingResidual` drive the combined `(rho_parent + rho_fermion + rho_boson)/H0²` to the target (initial-abundance fixed-point OR combined mode). The residual sums all three children instead of parent+DR.

- [ ] **Step 2: Add a full-run smoke test** in `species/dncdm_inv_test.cpp`: a background-only `Cosmology` with `dncdm1.inverse_decays=yes`, `dncdm1.quantum_statistics=no`, small `Gamma`; assert `GetOmega0Species("dncdm1")` is pinned within range and the run completes without throwing.

- [ ] **Step 3: Run.** `./build/cmake/test-dncdm-inv` — Expected PASS. If shooting fails to converge, seed the guess from `DegGuessFromOmegaToday` (already available on the parent).

- [ ] **Step 4: Commit.**
```bash
git add species/dncdm_inv_species.cpp species/dncdm_inv_test.cpp
git commit -m "dncdm_inv: combined-sector shooting closure + background smoke test"
```

## Task 3.4: Background validation ladder (A/B vs master + equilibrium physics)

**Files:** Create `python/tests/dncdm_inv_ab.py` (not a make target — a documented driver).

- [ ] **Step 1: Write the classyref A/B driver.** Per `reference_classyref_testing`: build `classyref` (fresh master) + this branch's `classy`, run the SAME `.ini` through both, compare with a scale-relative metric. Modes:
  - `--decay-only`: `inverse_decays=no` → background AND Cl/P(k) must match master `ncdm_decay_dr` to ~0.1% (this is the byte-for-byte-in-spirit regression anchor; TE zero-crossing handled by the scale-relative metric).
  - `--tiny-gamma`: `inverse_decays=yes` with `Gamma→0` → must reduce to the ΛCDM+ν limit (matches a no-DDM reference to ~0.1%).

- [ ] **Step 2: Run the decay-only regression.**
Run: `python python/tests/dncdm_inv_ab.py --decay-only`
Expected: max scale-relative deviation < 1e-3 on background columns and TT/EE/P(k). (Do NOT gate on `COMPARE_OUTPUT_REF` goldens for the new physics — they will be stale; regenerate only after intentional output changes.)

- [ ] **Step 3: Add background-physics checks** (thesis §7 ladder, design §7.3) to the driver: dec+inv reaches quasi-equilibrium (kink in parent density); final PSDs fit `1/(exp((q−μ)/T)±1)` with `μ_H = μ_l + μ_φ`; the non-relativistic limit (`α≫1`, heavy `m` / long `τ`) shuts inverse+qs off automatically and reproduces `ncdm_decay_dr`; the comoving sector energy `a⁴ρ` redshifts as radiation during equilibrium (conservation columns stay at machine precision). Compare morphology against Barenboim et al. fig. 2 (dec / dec+inv / dec+inv+qs `a⁴ρ` panels) qualitatively.

- [ ] **Step 4: Run.**
Run: `python python/tests/dncdm_inv_ab.py --equilibrium`
Expected: equilibrium reached; FD/BE fits within tolerance; NR limit matches `ncdm_decay_dr` to ~0.1%; `a⁴ρ_total` flat during equilibrium.

- [ ] **Step 5: Commit.**
```bash
git add python/tests/dncdm_inv_ab.py
git commit -m "dncdm_inv: classyref A/B driver + background validation ladder (decay-only regression, equilibrium fits, NR shutoff)"
```

---

# PHASE 4 — Perturbations

Daughter per-q hierarchies already exist (`DrPsdSpecies::PerturbDerivs`, Phase 2). Add the composite collision coupling via the kernel's per-ℓ operator. Decay-only perturbations must reproduce master exactly (the parent has NO pure-decay collision term — it cancels identically; inverse decays are the first thing to source it, design §4.3).

## Task 4.1: Precision parameter `l_max_dncdm_col`

**Files:** Modify `include/precision.h:319`, `source/input_module.cpp:2542`.

- [ ] **Step 1: Add the field.** In `include/precision.h` after `l_max_dr_col` (`:320`): `int l_max_dncdm_col = 17; /**< max multipole receiving the DNCDM inverse-decay collision term; must be <= l_max_ncdm */`.

- [ ] **Step 2: Wire the reader.** In `source/input_module.cpp` after `:2542` (`read(fc, "l_max_dr_col", l_max_dr_col);`) add `read(fc, "l_max_dncdm_col", l_max_dncdm_col);`. Add a `class_test(ppr->l_max_dncdm_col > ppr->l_max_ncdm, ...)` next to the existing `l_max_dr_col` checks in `perturbations_module.cpp:542-546`.

- [ ] **Step 3: Verify build.**
Run: `cmake --build build/cmake --target class --parallel 2>&1 | tail -3`
Expected: builds clean.

- [ ] **Step 4: Commit.**
```bash
git add include/precision.h source/input_module.cpp source/perturbations_module.cpp
git commit -m "precision: l_max_dncdm_col collision-truncation parameter"
```

## Task 4.2: Perturbation-operator ℓ=0/ℓ=1 identities (kernel TDD)

**Files:** Modify `tools/decay_kernel_test.cpp`, `tools/decay_transition_kernel.cpp` (implement `ApplyPerturbationOperator`).

**Interfaces:**
- Consumes: `PrepareTransitions` scratch from Task 1.2.
- Produces: `ApplyPerturbationOperator` trusted by the composite `AddCouplingDerivs`.

- [ ] **Step 1: Write the failing operator-identity tests.** After `PrepareTransitions(a,m,Γ,f...)`, for random Ψ arrays:

**OPERATOR SEMANTICS (Fable review — load-bearing):** `ApplyPerturbationOperator`
returns `dΨ_i = (1/f̄_i)(dF_i/dτ)^(1)_{C,ℓ}` — the pure F-space collision rate
divided by f̄. It must NOT include the `−(f̄̇_i/f̄_i)Ψ_{i,ℓ}` background-evolution
subtraction; that term is added by the composite in Task 4.3 (where, for the
parent in decay-only mode, it cancels the operator's decay piece EXACTLY — the
paper §4.3 cancellation — so getting the split right is what makes the decay-only
regression pass). `PrepareTransitions` already computes the background `df̄` arrays;
cache them and expose `df_bg_H()/df_bg_l()/df_bg_phi()` accessors so the composite
can form `f̄̇/f̄` from the SAME discrete rates (never from a separately-interpolated
background column — mismatched rates would break the cancellation).

```cpp
// ℓ=0: the perturbed operator conserves number & energy (design §4.3, f̄→f̄Ψ linearization).
// Identities use the SAME bare-moment/g-factor convention as Task 1.2 (×2 boson dof
// folded into the boson scatter): with num_i = Σ dq q² f̄_i dΨ_i,0, e_i likewise with ε:
K.ApplyPerturbationOperator(0, PsiH.data(), Psil.data(), Psiphi.data(),
                            dPsiH.data(), dPsil.data(), dPsiphi.data());
assert(std::fabs(numH + numl) < 1e-11*scale);        // fermion-leg number
assert(std::fabs(2.*numH + numphi) < 1e-11*scale);   // boson-leg number
assert(std::fabs(2.*eH + 2.*el + ephi) < 1e-11*scale); // g-weighted energy
// ℓ=1: momentum identity — q⃗1 = q⃗2 + q⃗3 projects to q1 = q2*cosα* + q3*cosβ* at each node
//   (formalism-notes ℓ=1). With P_1(cosθ)=cosθ the first-moment scatter preserves
//   Σ_i g_i ∫dq q³ Ċ_{i,1}=0, which in bare moments (mom_i = Σ dq q³ f̄_i dΨ_i,1) reads:
K.ApplyPerturbationOperator(1, ...);
assert(std::fabs(2.*momH + 2.*moml + momphi) < 1e-10*scale);
```
Also assert `ApplyPerturbationOperator` builds `P_ℓ` by recurrence (spot-check `P_2 = 1.5cos²−0.5` at a node via `NodeKinematics`), never from a table (research-branch lesson 5).

- [ ] **Step 2: Run to verify it fails.** `./build/cmake/test-decay-kernel` — Expected FAIL (empty operator).

- [ ] **Step 3: Implement `ApplyPerturbationOperator`.** Per thesis eqs. 6.4–6.8: for each parent bin/node, gather `Ψ_H(q1_i)`, `Ψ_l(q2*)`, `Ψ_φ(q3*)` (Ψ_daughter via the SAME two-bin split used for `f`), form the collision rate for multipole ℓ using the cached `Λ_s`/`f_l`/`f_φ` and the angular factors `P_ℓ(cosα*)`, `P_ℓ(cosβ*)`, `P_ℓ(cosγ*)` (recurrence at the exact node cosines), and scatter into `dPsi_H/dPsi_l/dPsi_phi` with the transpose weights. The `1/f_i` normalization (thesis `k_i = .../(q² f_i)`) is applied per daughter bin with a guard `f_i < kFFloor → skip` (empty bins have no collision). Fold the fermion/boson `∓` QS sign and the boson ×2 dof exactly as in the background scatter.

- [ ] **Step 4: Run to verify it passes.** `./build/cmake/test-decay-kernel` — Expected PASS.

- [ ] **Step 5: Commit.**
```bash
git add tools/decay_transition_kernel.cpp tools/decay_kernel_test.cpp
git commit -m "kernel: per-ℓ perturbation collision operator — discrete ℓ=0 (N/E) and ℓ=1 (momentum) identities"
```

## Task 4.3: Composite `AddCouplingDerivs` via the operator

**Files:** Modify `species/dncdm_inv_species.cpp`, `species/dncdm_inv_test.cpp`.

- [ ] **Step 1: Implement `AddCouplingDerivs`.** In the composite (`composite_species.h:206` seam; called after each child's free-streaming `PerturbDerivs`): read the current `(a,m,Γ)` and the three background PSDs from `ppw->pvecback` (the interpolated background columns — `index_bg_lnf_decay_dr1_`/daughters' `index_bg_f_`); `kernel_->PrepareTransitions(a,m,Γ, f_H, f_l, f_phi)` ONCE; then for `l = 0 .. min(l_max_dncdm_col, layout.l_max)` gather the per-ℓ Ψ slices from `y` (parent `dncdm_layout(...).index_per_q[iq]+l`, daughters analogously) and `kernel_->ApplyPerturbationOperator(l, ...)`, accumulating into the matching `dy` slots. **Then add the background-evolution term (Fable review):** for every species, bin iq, and ℓ ≤ layout.l_max (NOT truncated at l_max_dncdm_col — f̄̇/f̄ multiplies every multipole):
```cpp
dy[Psi_i(iq,l)] += kernel_dPsi_i[iq]           // (1/f̄)(dF/dτ)_C,ℓ   (l ≤ l_max_dncdm_col)
                 − (kernel_->df_bg_i()[iq] / f_i[iq]) * y[Psi_i(iq,l)];  // −(f̄̇/f̄)Ψ, all ℓ
```
using the kernel's own cached `df_bg` (same discrete rates ⇒ the parent's decay piece cancels to machine precision in decay-only mode, reproducing the paper §4.3 cancellation; a separately-interpolated f̄̇ would leave a spurious residual collision term). Guard `f_i[iq] < kFFloor → skip both terms` (empty bins carry no collision). Guards: `ppw->approx[index_ap_ncdmfa]` must be off for the coupled trio (already forbidden at input); synchronous only (guarded at Create). Use `dlnf/dlnq` live from the background column (never a frozen table) — it is already published by both `DNCDMSpecies::ComputeBackground` (collision-owned) and `DrPsdSpecies::ComputeBackground`.

- [ ] **Step 2: Write the composite-level discrete-identity test** in `species/dncdm_inv_test.cpp`: assemble a synthetic perturbed state (background f-columns + Ψ slots in a hand-built `y`), call the composite's `PerturbDerivs`, and assert the ℓ=0/ℓ=1 conservation of the collision contribution (subtract the free-streaming part by running with `Gamma=0` as the baseline). This is the perturbation analogue of the Phase-3 conservation-column check.

- [ ] **Step 3: Run.** `cmake --build build/cmake --target test-dncdm-inv --parallel && ./build/cmake/test-dncdm-inv` — Expected PASS.

- [ ] **Step 4: Commit.**
```bash
git add species/dncdm_inv_species.cpp species/dncdm_inv_test.cpp
git commit -m "dncdm_inv: AddCouplingDerivs collision coupling via kernel per-ℓ operator (synchronous, l_max_dncdm_col)"
```

## Task 4.4: Perturbation validation (regression + convergence)

**Files:** Modify `python/tests/dncdm_inv_ab.py`.

- [ ] **Step 1: Decay-only perturbation regression.** Extend the driver: `inverse_decays=no` full Cl/P(k) must be IDENTICAL to master `ncdm_decay_dr` to ~0.1% (the parent has no pure-decay collision term, so this is a strong anchor). Prefer `evolver=2` (rkdp45) for the coupled runs per the dcdm_wdm 14× lesson, but verify BOTH evolvers complete (the nonlinear coupling can defeat ndf15's numjac sparsity lock — research-plugins §3.2/§8 checklist item 6).

- [ ] **Step 2: Parent→daughter convergence in equilibrium.** `inverse_decays=yes` (qs off): assert the parent perturbations `δ,θ,σ` converge to the daughters' during the equilibrium phase (thesis §7 check 5, Fig. 6.9) — they never converge without inverse decays.

- [ ] **Step 3: Run.**
Run: `python python/tests/dncdm_inv_ab.py --pert-regression --pert-convergence`
Expected: decay-only Cl/P(k) match to ~0.1% on both evolvers; parent/daughter perturbations converge in equilibrium.

- [ ] **Step 4: Commit.**
```bash
git add python/tests/dncdm_inv_ab.py
git commit -m "dncdm_inv: perturbation validation (decay-only regression both evolvers, equilibrium convergence)"
```

---

# PHASE 5 — Lasing, diagnostics, docs, PR prep

## Task 5.1: Lasing run + q_min convergence of observables

**Files:** Modify `python/tests/dncdm_inv_ab.py`, `species/dncdm_inv_test.cpp`.

- [ ] **Step 1: Lasing configuration.** Add a `--lasing` mode: `m=0.5 eV`, large `Γ`, `quantum_statistics=yes` (thesis §5 / §7 check 4). Assert observables (`ρ_i`, `N_i`, the fermion-energy spike) are BOUNDED and q_min-convergent as `dr_q_min` is lowered (1e-2 → 1e-3 → 1e-4), even while `f_φ(q_min)` grows large. Assert the conservation columns stay at machine precision throughout (this is what exact conservation buys — the runaway saturates when the inversion `∫(f_H−f_l)` is consumed, design §4.2).

- [ ] **Step 2: `bose_occupancy` diagnostic column** (design §4.2): add a background column flagging condensate-like behaviour (`max_q f_φ`), written by `DNCDMInvSpecies::WriteBackgroundData`.

- [ ] **Step 3: Run.**
Run: `python python/tests/dncdm_inv_ab.py --lasing`
Expected: `ρ_φ`, `N_φ`, fermion-energy spike converge in `q_min` to a few %; conservation residual < 1e-10; `bose_occupancy` column present. (If observables do NOT converge, that is a finding, not a plan failure — record it; the design's central claim is that exact conservation makes the observables q_min-convergent even when `f(q_min)` diverges.)

- [ ] **Step 4: Commit.**
```bash
git add python/tests/dncdm_inv_ab.py species/dncdm_inv_species.cpp
git commit -m "dncdm_inv: lasing run + q_min convergence of observables + bose_occupancy diagnostic"
```

## Task 5.2: Anisotropic-stress diagnostic (the physics deliverable)

**Files:** Modify `species/dncdm_inv_species.cpp` (perturbation output), `python/tests/dncdm_inv_ab.py`.

- [ ] **Step 1: Emit `a⁴Π_νφ`** (thesis eq. 6.19): sum `C Σ_i g_i ∫ q² dq (q²/ε_i)(f_i + Ψ_{2,i})` over the three species into a perturbation-output column via `WriteOutputColumns`/`FillSources`.

- [ ] **Step 2: Measure suppression with qs OFF** (design §7.6, Emil's explicit recommendation): run dec+inv (qs off) in the relativistic regime and report whether `a⁴Π_νφ` is suppressed during equilibrium and with which effective γ-power (`Γ_T ≃ Γγ³` vs `Γγ⁵`). Output a plot/table; this is deliverable-grade output for the γ³/γ⁵ dispute, not a pass/fail gate.

- [ ] **Step 3: Run.**
Run: `python python/tests/dncdm_inv_ab.py --anisotropic-stress`
Expected: `a⁴Π_νφ(τ)` emitted for `k ∈ {0.08, 0.15, 0.5} Mpc⁻¹`; suppression (or its absence) recorded.

- [ ] **Step 4: Commit.**
```bash
git add species/dncdm_inv_species.cpp python/tests/dncdm_inv_ab.py
git commit -m "dncdm_inv: a^4 Pi_nuphi anisotropic-stress diagnostic (gamma^3 vs gamma^5 deliverable)"
```

## Task 5.3: Docs, cleanup, PR prep

**Files:** Modify docs/notebooks as appropriate; run `simplify`/`/code-review`.

- [ ] **Step 1: Full test sweep.**
Run: `make test`
Expected: all `TEST_TARGETS` (including `test-decay-kernel`, `test-dr-psd`, `test-dncdm-inv`) PASS via ctest.

- [ ] **Step 2: Cleanup.** Remove any dead debug scaffolding (no `int wef=0;`-style breakpoints, research-branch lesson 10). Confirm comments state constraints not narration. Confirm no `cclassy.pxd` hand-edits, no `git add -A` in history.

- [ ] **Step 3: Documentation.** Add a short tutorial notebook or doc section for `dncdm1.inverse_decays`/`quantum_statistics`/`dr_*`/`l_max_dncdm_col`, referencing the design doc. Add the new species to the docs site CURATED list only if executed cleanly (see `project_docs_site_pr381`).

- [ ] **Step 4: Review + PR.** Run `/code-review` on the diff; address findings. Open the PR referencing arXiv:2011.01502 and the design doc, with the validation-ladder results (decay-only regression, equilibrium fits, lasing convergence, anisotropic-stress table) in the body.

- [ ] **Step 5: Commit.**
```bash
git add <explicit doc/notebook paths>
git commit -m "dncdm_inv: docs, cleanup, PR prep"
```

---

## Execution notes for the orchestrator

- **Phase ordering is strict** — each phase gates the next (design §7 validation ladder). No cross-phase parallelism. Within a phase, TDD lets test-writing precede implementation (the test-first steps are explicit).
- **Riskiest steps, in order:**
  1. **Task 1.2 kernel weight algebra** (Phase 1) — the transpose gather/scatter consistency is what makes conservation + detailed balance exact and replaces the failed A+B·q calibration. If the `1e-11`/`1e-12` tolerances can't be met, the weights are not transpose-consistent — FIX, never loosen (research-thesis lesson 5, design §4.2). Everything downstream trusts this.
  2. **Task 3.2 composite index registration** (Phase 3) — the `child_layouts` contract (#358/#359) and the collision-owned parent's `bi_f_parent_index` wiring. A mis-registered slot segfaults the hot path or silently corrupts the RHS. The switch-copy test (Task 3.2 Step 5) and `test-composite-layout` are the guards.
  3. **Task 4.3 nonlinear-coupling vs ndf15 numjac** — the kernel RHS is nonlinear in the daughter occupations, which defeats the `wdm_decay_product` static-sparsity floor trick (research-plugins §3.2 caveat). Default the coupled runs to `evolver=2` (rkdp45); verify ndf15 completes but treat it as secondary.
  4. **`DrPsdSpecies` dlnf/dlnq live spline** through `f→0` (Task 2.1/4.3) — the branch's most turbulent area (research-branch §3/§7 dlnf/dlnq saga). Positivity floor + `kFSeed` + Ψ=0 on empty bins; the moment-identity/convergence tests are the guards.
- **A/B against master via classyref** (`reference_classyref_testing`): the ONLY use of the classyref harness is confirming the CLASSICAL LIMIT (`inverse_decays=no`, or `Gamma→0`) reproduces master to ~0.1% with a scale-relative metric (TE zero-crossings). The full inverse+qs physics will NOT match master — do not gate it on `COMPARE_OUTPUT_REF` goldens (they go stale; regenerate after intentional changes).
- **Build/run cheat sheet:** `make <target>` builds via the CMake shim into `build/cmake/`; run `./build/cmake/<target>` (working dir = repo root) or `ctest --test-dir build/cmake -R <target> --output-on-failure`. `make test` builds + runs the whole `TEST_TARGETS` set. Every new test MUST be in BOTH `CMakeLists.txt` and `Makefile TEST_TARGETS` or CI silently skips it.
- **Scope discipline (design §8, out of scope for v1):** massive daughters (keep the kernel kinematics general but implement/test massless only), Newtonian gauge, tensors with inverse on, fluid approximation for the trio, the DCDM-approximation speed switch, 2↔2 scattering, multiple simultaneous inverse-decay families. Do not gold-plate these.

## Self-Review

- **Spec coverage:** Phase 1 ↔ design §4 + §7.1 (kernel, conservation, detailed balance, ℓ=1, continuity). Phase 2 ↔ §5 `dr_psd`. Phase 3 ↔ §5 composite + §7.2/§7.3 (decay-only regression, equilibrium, NR shutoff, conservation columns) + the resolved parent-ownership question. Phase 4 ↔ §4.3 + §7.5 (operator identities, decay-only regression, convergence, `l_max_dncdm_col`). Phase 5 ↔ §7.4/§7.6 (lasing, anisotropic stress). Guards (§5), grids (§6), evolver (§5) all placed.
- **Design question:** resolved explicitly (mode flag on `DNCDMSpecies`, table of gated behaviours, decay-only byte-for-byte).
- **Type consistency:** `DecayTransitionKernel`, `GridView`, `Statistics`, `Config`, `Moments`, `ApplyPerturbationOperator`, `DrPsdSpecies` accessors (`bi_f_index`, `bg_dlnfdlnq_index`, `q_bg`), `DNCDMSpecies::SetCollisionOwned`/`bi_f_parent_index`, `DNCDMInvSpecies::Create`, `ChildIndex{kParent,kFermion,kBoson}` are used consistently across tasks.
- **Both build files** updated for each of the three new test executables (`test-decay-kernel`, `test-dr-psd`, `test-dncdm-inv`).
