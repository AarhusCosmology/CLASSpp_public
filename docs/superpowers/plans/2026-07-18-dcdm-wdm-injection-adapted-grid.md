# dcdm_wdm Injection-Adapted Momentum Grid — Implementation Plan (v3)

> **Executed 2026-07-23; Task-2 gate and Task-3 comparison exposed two field bugs fixed in
> `5addd728` and `1586d250` — see ledger below.** All four tasks are done or superseded by an
> amendment recorded in the task reports (`.superpowers/sdd/task-{1,2,3,4}-report.md`).

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the landed order-1 CIC injection deposit with a smooth, conservative **cell-integrated Gaussian (erf) energy deposit** on the (kept) injection-adapted quantile grid, restore `ndf15` compatibility, prove the hang is gone, and deliver the old-vs-new convergence comparison that previous sessions never completed.

**Architecture:** The quantile grid (equal decayed fraction `F = 1−e^{−Γt}` per bin, fiducial `H(t)` map, grid = f(Γ) only) is committed and correct — untouched. The deposit becomes: normal-CDF differences at the `N+1` cell edges with off-grid tails clamped into the edge bins (`Σ W = 1` identically ⇒ energy sum-rule exact, no D-renormalization), width `σ = clamp(Δu_loc, 0.03, 0.25)` where `Δu_loc` is the *analytic* local quantile-bin width — smooth (`C^∞`) in time and always resolvable by the background table (`back_integration_stepsize = 7e−3` in ln a). `∂J/∂ln q` is the cell-averaged profile derivative via normal-pdf differences at the edges (no edge special cases). Both evolvers work; the `ndf15` constructor guard is deleted.

**Tech Stack:** C++17, CLASSpp species framework, CMake/CTest; Python (classy + classyref) for validation.

**History:** v1 grid + v2 CIC landed as commits `93b9bf00`, `5539964c`, `2340f3f3`. Spec v3 (`d61e8f8f`) withdrew the CIC: it banned ndf15 AND effectively hung RKDP45 (8h 100%-CPU sweep, root-caused to injection-column aliasing in `AddCouplingDerivs`). This plan is the v3 delta.

## Global Constraints

- **Spec (authoritative):** `docs/superpowers/specs/2026-07-18-dcdm-wdm-adaptive-momentum-grid-design.md`, §4.3 v3.
- **C++-only** — plain C++ headers, no `extern "C"` / `#ifdef __cplusplus`.
- **Never hand-edit `cclassy.pxd`** — auto-generated at build.
- **Exact energy sum-rule** — `Σᵢ dqᵢ qᵢ² εᵢ Jᵢ = A` must hold to ~1e−12 for every `u_cut` (on-grid and off). The erf weights guarantee it via `Σ Wᵢ = 1` (telescoping + edge clamps). Never break this.
- **Grid untouched** — `BuildInjectionAdaptedGrid`, `q_edge_tol`, `q_min_ratio` floor, fallback grid: no changes.
- **`kFSeed = 1e-10` untouched.**
- **`-ffast-math` discipline** — clamp every `exp` argument (master-style `x² / (2σ²) ≥ 60 → 0`), clamp `erf` arguments (`|x| > 8 → ±1`), never materialize `inf` (evaluate `Δu_loc` in log form and compare before exponentiating).
- **Numeric-range input checks use `class_test`**, structural ones `class_test_severe`.
- **No new test executable** — `test-dcdm-wdm` is registered in both `CMakeLists.txt` and Makefile `TEST_TARGETS`.
- **Build + run (repo root):** `make test-dcdm-wdm && ./build/cmake/test-dcdm-wdm` → ends `dcdm_wdm daughter test passed`.
- **Do NOT run `git checkout`/`switch`/`reset`/`branch` anywhere, including in `.worktrees/master-oldgrid`** — stay on `dcdm-wdm-injection-adapted-grid`; `git add` only the explicit paths named per step.
- **Verification tolerance philosophy:** ~0.1% scale-relative physics, Cl^TE zero-crossing aware; never blind max-rel-diff, never bit-identical.
- **Long runs must be hang-proof:** subprocess per run + hard timeout + incremental result saving. Never launch an unbounded `c.compute()` loop in the foreground.

**Key files:** `species/wdm_decay_product.{h,cpp}`, `species/dcdm_wdm_test.cpp`, `notebooks/dcdm_wdm_convergence.py`, `notebooks/dcdm_wdm_convergence_compare.py`, `notebooks/dcdm-wdm-massive-decay-products.ipynb`.

---

## Task 1: Erf deposit + ndf15 restoration (C++, TDD)

**Files:**
- Modify: `species/wdm_decay_product.h`
- Modify: `species/wdm_decay_product.cpp`
- Test: `species/dcdm_wdm_test.cpp`

**Interfaces:**
- Consumes (unchanged, already in tree): `BuildInjectionAdaptedGrid()` and members `q_`, `u_`, `dq_`, `u_edge_`, `scratch_w_`; `FiducialCosmicTime(a)`; `Gamma_`, `M_`, `factor_`, `kQKick`, `n_bins_`.
- Produces: `double SigmaAt(double u_cut) const` (public), `void GaussWeights(double u_cut, double sigma, double* w) const` (private), reworked `FillInjection`, new members `F_lo_`, `F_hi_`, new public accessor `u_edge()`. Removed: `CicWeights`, the ndf15 constructor guard.

- [x] **Step 1: Adapt the tests (RED first)** — edits to `species/dcdm_wdm_test.cpp`:

  (a) In `BaseFc()`, **delete** the comment block and the line `fc.set("evolver", "2");` (near line 38–42). Also delete every other `fc.set("evolver", ...)` line in the guard-test blocks (near lines 150, 169, 182) — the default evolver is valid again.

  (b) **Delete** the whole evolver-guard test block (the `{ ... }` around lines 190–200 containing `fc.set("evolver", "1");  // ndf15 -> reject` and `assert(Throws(fc))` / `assert(!Throws(fc))`).

  (c) **Replace** the CIC centroid block (the `{ ... }` titled `// ── Order-1 CIC deposit: energy centroid (first moment) is exactly u_cut ────`, near lines 257–287) with:

```cpp
  // ── Erf deposit: interior energy centroid sits at u_cut to within the kernel
  // width. The cell-integrated Gaussian is symmetric in u, so the deposited
  // energy centroid tracks the kick location; sub-cell effects bound the error
  // by a fraction of sigma. (The old CIC's 1e-9 exactness was a stencil
  // property, not a physical requirement — the over-time golden is the
  // physical requirement.)
  {
    FileContent fc = BaseFc();
    WdmDecayProductSpecies sp(&fc, "ddm", TestSettings(), &pba, nullptr);
    const int N = sp.q_size();
    std::vector<double> J(N);
    for (int idx : {N / 4, N / 2, 3 * N / 4}) {
      const double u_cut = 0.5 * (sp.u()[idx] + sp.u()[idx + 1]);  // between centres
      const double a     = std::exp(u_cut) / 10.0;                 // u_cut = ln(a*kQKick)
      const double sigma = sp.SigmaAt(u_cut);
      sp.FillInjection(a, 0.7, J.data(), nullptr);
      double A = 0., Au = 0.;
      bool negative = false;
      for (int i = 0; i < N; ++i) {
        const double q  = sp.q()[i];
        const double ep = std::sqrt(q * q + a * a * sp.M() * sp.M());
        const double e  = sp.dq()[i] * q * q * ep * J[i];  // injected energy in bin i
        A  += e;
        Au += e * sp.u()[i];
        if (J[i] < 0.)
          negative = true;
      }
      assert(!negative);                             // weights in [0,1]
      assert(A > 0.);
      assert(std::fabs(Au / A - u_cut) < 0.5 * sigma);  // centroid within the kernel
      assert(sigma >= 0.03 - 1e-12 && sigma <= 0.25 + 1e-12);  // clamp respected
    }
  }
```

  (d) In the bottom-edge clamp block (the `{ ... }` with `auto gate_ratio = [&](double u_cut)`, near lines 225–255): keep the three `r_at/r_mid/r_far == 1` asserts (the erf clamp conserves energy identically). **Replace** the trailing `dJ` asserts (the lines from `assert(dJ[0] < 0.);` through the loop asserting `dJ[i] == 0.`) with:

```cpp
      // g-source structure below the grid: the deposited profile inside the
      // grid is the falling upper tail of the (off-grid) Gaussian, so every
      // bin's dJ/dlnq <= 0, strictly negative in the bottom bin half a cell
      // below the edge (where the pdf is still finite).
      sp.FillInjection(std::exp(u0 - 0.5 * sp01) / 10.0, 0.7, J.data(), dJ.data());
      assert(dJ[0] < 0.);
      for (int i = 0; i < N; ++i)
        assert(dJ[i] <= 0.);
```

  (e) In the high-side clamp block (near lines 289–322): keep `sum_w == 1` and `J[N-1] > 0.` asserts. **Replace** the trailing `dJ` asserts (`assert(dJ[N - 1] > 0.);` plus the `for (i < N-1) assert(dJ[i] == 0.)` loop) with:

```cpp
    // Above the grid the in-grid profile is the rising lower tail of the
    // off-grid Gaussian: dJ/dlnq >= 0 everywhere, strictly positive at the top.
    assert(dJ[N - 1] > 0.);
    for (int i = 0; i < N; ++i)
      assert(dJ[i] >= 0.);
```

  (f) In the `dJ/dlnq structure` block (near lines 323–347): the sign asserts are correct for the erf deposit *except* possibly in the single cell containing the cutoff (bin centres are F-medians, not u-midpoints). **Replace** the loop guard so the straddling cell is skipped precisely:

```cpp
    for (int i = 1; i < N - 1; ++i) {
      if (J[i] < 1e-30)
        continue;
      if (sp.u_edge()[i + 1] < u_cut)   // cell entirely below the cutoff
        assert(dJ[i] >= 0.);
      if (sp.u_edge()[i] > u_cut)       // cell entirely above the cutoff
        assert(dJ[i] <= 0.);
    }
```

  If the `J` peak-location assert in the same block (`|u_peak − u_cut| ≤ 1.5 spacings`) fails from the kernel's `q³ε`-measure tilt, loosen 1.5 → 2.5 with a one-line comment; do not loosen further without investigation.

  (g) In the golden-run comment (near line 454): update the sentence "The CIC deposit clamps below-grid injection..." to "The deposit clamps below-grid injection into the bottom bin (Σ W = 1 always), so energy is conserved from the earliest decays". The golden now runs under the **default evolver (ndf15)** because BaseFc no longer forces `evolver=2` — that *is* the hang-regression test; do not re-add an evolver key.

- [x] **Step 2: Run to verify RED**

Run: `make test-dcdm-wdm 2>&1 | tail -20`
Expected: compile FAILURE — `SigmaAt` and `u_edge()` don't exist yet. (If it compiles, the test edits weren't applied.)

- [x] **Step 3: Header changes** — in `species/wdm_decay_product.h`:

  (a) Public section, next to the `u()` accessor, add:

```cpp
  const std::vector<double>& u_edge() const {
    return u_edge_;
  }
  /** Injection kernel width at u_cut: the analytic local quantile-bin width
   *  Δu_loc = [(F_hi−F_lo)/N] · H_fid(a) e^{+Γ t_fid(a)} / Γ  (a = e^{u_cut}/kQKick),
   *  clamped to [kSigmaMin, kSigmaMax]. Smooth in u_cut; evaluated in log form so
   *  no overflow/inf is ever materialized (-ffast-math discipline). Floor = the
   *  background-table resolvability bound; cap = the placement-bias budget. */
  double SigmaAt(double u_cut) const;
```

  (b) Replace the `CicWeights` declaration with:

```cpp
  /** Cell-integrated Gaussian (erf) energy-deposit weights for a kick at u_cut:
   *  w[i] = Φ((u_edge[i+1]−u_cut)/σ) − Φ((u_edge[i]−u_cut)/σ), plus the below-/
   *  above-grid tail mass clamped into the corresponding edge bin, so Σ w = 1
   *  identically (energy conserved for every u_cut, no D-renormalization) and
   *  the source is C-infinity in time. */
  void GaussWeights(double u_cut, double sigma, double* w) const;
```

  (c) Private members: after `q_min_ratio_floor_`, add:

```cpp
  static constexpr double kSigmaMin = 0.03;  // ≈ 4 background-table samples (7e-3 in ln a)
  static constexpr double kSigmaMax = 0.25;  // bounds the wide-bin placement bias
  double F_lo_ = 0.;  // realized decayed-fraction span of the grid (0,0 = fallback grid)
  double F_hi_ = 0.;
```

  (d) Update the `FillInjection` doc comment (the `/** Injection source Jᵢ ... */` block): replace the CIC description with:

```cpp
  /**
   * Injection source Jᵢ = dfᵢ/dτ (and optionally dJᵢ/dln q) for the current
   * (a, ρ_dcdm). Conservative cell-integrated Gaussian (erf) energy deposit in
   * u = ln q at u_cut = ln(a·kQKick): the kick's energy A = aΓρ_dcdm·a⁴/factor_
   * is split by GaussWeights over the cells (σ = SigmaAt(u_cut)), off-grid tails
   * clamped into the edge bins, so Σ Wᵢ = 1 identically and the exact sum rule
   * Σ dqᵢ qᵢ² εᵢ Jᵢ = A holds for every u_cut. dJ/dln q is the cell-averaged
   * profile derivative (normal-pdf differences at the cell edges) — smooth, no
   * edge special cases. dJdlnq may be nullptr. Background-thread only
   * (uses mutable scratch).
   */
```

  (e) Update the class-level doc comment: replace the sentence about the "conservative order-1 CIC energy deposit" with "injection via a conservative cell-integrated Gaussian (erf) energy deposit in ln q — the energy-injection sum rule Σᵢ dqᵢ qᵢ² εᵢ Jᵢ · factor/a⁴ = aΓρ_dcdm holds exactly for every cutoff position, and the source is smooth in time (both evolvers integrate it)".

- [x] **Step 4: Constructor + grid bookkeeping** — in `species/wdm_decay_product.cpp`:

  (a) **Delete** the whole ndf15-guard block in the constructor (the comment starting `// The order-1 CIC injection deposit is only C0 ...` through the closing `}` of `if (pfc != nullptr) { ... }`, lines ~71–84).

  (b) In `BuildInjectionAdaptedGrid()`: after `F_lo`/`F_hi` are computed, store them; in the fallback branch, zero them. I.e. immediately before `q_.resize(N);` insert:

```cpp
  F_lo_ = F_lo;
  F_hi_ = F_hi;
```

  and inside the fallback branch (`if (!(F_hi > F_lo + q_edge_tol_)) {`), first lines:

```cpp
    F_lo_ = 0.;
    F_hi_ = 0.;  // signals SigmaAt to use the uniform fallback spacing
```

- [x] **Step 5: Implement `SigmaAt` + `GaussWeights`** — in `species/wdm_decay_product.cpp`, **replace the entire `CicWeights` definition** with:

```cpp
double WdmDecayProductSpecies::SigmaAt(double u_cut) const {
  // Fallback (negligible-decay) grid: uniform spacing, constant width.
  if (!(F_hi_ > F_lo_)) {
    const double du = u_edge_[1] - u_edge_[0];
    return std::min(kSigmaMax, std::max(kSigmaMin, du));
  }
  // Analytic local quantile-bin width (spec §4.3 v3):
  //   dF/du = Γ e^{−Γ t_fid(a)} / H_fid(a),  a = e^{u_cut}/kQKick
  //   Δu_loc = [(F_hi−F_lo)/N] · H_fid(a) e^{+Γ t_fid(a)} / Γ
  // Γ t_fid can exceed 700 (fast decays, late times): compare in log form and
  // return the cap BEFORE exponentiating — never materialize inf (-ffast-math).
  const double a  = std::exp(u_cut) / kQKick;
  const double H0 = kFidH * 1.0e5 / _c_;
  const double H  = H0 * std::sqrt(kFidOmegaM / (a * a * a) + kFidOmegaR / (a * a * a * a));
  const double gt = Gamma_ * FiducialCosmicTime(a);
  const double ln_du =
      std::log((F_hi_ - F_lo_) / n_bins_ * H / Gamma_) + gt;
  if (ln_du >= std::log(kSigmaMax))
    return kSigmaMax;
  return std::max(kSigmaMin, std::exp(ln_du));
}

void WdmDecayProductSpecies::GaussWeights(double u_cut, double sigma, double* w) const {
  const int N        = n_bins_;
  const double inv_s = 1.0 / (sigma * std::sqrt(2.0));
  // Normal CDF at a cell edge; erf saturates by |x| ~ 6, clamp at 8 so no
  // huge-argument libm call happens under -ffast-math.
  auto cdf = [&](double u_e) {
    const double x = (u_e - u_cut) * inv_s;
    if (x >= 8.)
      return 1.0;
    if (x <= -8.)
      return 0.0;
    return 0.5 * (1.0 + std::erf(x));
  };
  double c_lo = cdf(u_edge_[0]);
  w[0]        = c_lo;  // below-grid tail clamped into the bottom bin
  for (int i = 0; i < N; ++i) {
    const double c_hi = cdf(u_edge_[i + 1]);
    if (i > 0)
      w[i] = 0.;
    w[i] += c_hi - c_lo;
    c_lo = c_hi;
  }
  w[N - 1] += 1.0 - c_lo;  // above-grid tail clamped into the top bin
  // Σ w = 1 identically for every u_cut (telescoping) — the exact energy
  // sum-rule needs no renormalization.
}
```

- [x] **Step 6: Replace `FillInjection`** — in `species/wdm_decay_product.cpp`, replace the whole function body with:

```cpp
void WdmDecayProductSpecies::FillInjection(double a,
                                           double rho_dcdm,
                                           double* J,
                                           double* dJdlnq) const {
  const int N        = n_bins_;
  const double u_cut = std::log(a * kQKick);

  auto zero_all = [&]() {
    for (int i = 0; i < N; ++i) {
      J[i] = 0.;
      if (dJdlnq != nullptr)
        dJdlnq[i] = 0.;
    }
  };
  if (rho_dcdm <= 0. || Gamma_ == 0.) {
    zero_all();
    return;
  }

  // Conservative cell-integrated Gaussian (erf) energy deposit (spec §4.3 v3):
  // smooth in time (both evolvers; tabulated injection columns stay
  // representable), Σ w = 1 identically (exact sum rule, edges clamped).
  const double a2    = a * a;
  const double A     = a * Gamma_ * rho_dcdm * a2 * a2 / factor_;
  const double sigma = SigmaAt(u_cut);
  GaussWeights(u_cut, sigma, scratch_w_.data());

  for (int i = 0; i < N; ++i) {
    const double w = scratch_w_[i];
    if (w == 0.) {
      J[i] = 0.;
      continue;
    }
    const double epsilon = std::sqrt(q_[i] * q_[i] + a2 * M_ * M_);
    J[i]                 = A * w / (dq_[i] * q_[i] * q_[i] * epsilon);
  }

  if (dJdlnq == nullptr)
    return;
  // g-source: cell-averaged d/dlnq of the deposited profile, via the normal
  // pdf at the cell edges. Positive below the cutoff (profile rising),
  // negative above, ~0 deep inside a wide cell or off-grid — no edge special
  // cases, and Σ dq q² ε dJdlnq telescopes to ~0.
  const double inv_s = 1.0 / sigma;
  const double norm  = inv_s / std::sqrt(2.0 * M_PI);
  auto pdf = [&](double u_e) {
    const double x  = (u_e - u_cut) * inv_s;
    const double x2 = 0.5 * x * x;
    return (x2 >= 60.) ? 0.0 : norm * std::exp(-x2);  // master-style exp clamp
  };
  double p_lo = pdf(u_edge_[0]);
  for (int i = 0; i < N; ++i) {
    const double p_hi = pdf(u_edge_[i + 1]);
    const double diff = p_hi - p_lo;
    p_lo              = p_hi;
    if (diff == 0.) {
      dJdlnq[i] = 0.;
      continue;
    }
    const double epsilon = std::sqrt(q_[i] * q_[i] + a2 * M_ * M_);
    dJdlnq[i]            = A * diff / (dq_[i] * q_[i] * q_[i] * epsilon);
  }
}
```

- [x] **Step 7: Run to verify GREEN**

Run: `make test-dcdm-wdm && ./build/cmake/test-dcdm-wdm`
Expected: PASS — ends `dcdm_wdm daughter test passed`. Watch for: every `energy conservation at z=...: rel. dev. ...` line `< 5e-3` (this now runs under the **default ndf15** — completing at all is the hang-regression check), and `full perturbed dcdm_wdm pipeline ran`.

Timing note: if the perturbed smoke block takes > ~3 min under ndf15, add `fc.set("evolver", "2");` to that block only, with the comment `// speed only — both evolvers work (see golden above, which runs on ndf15)`, and note the measured times in the report.

If a golden line is `>= 5e-3`: STOP, use superpowers:systematic-debugging. Prime suspects: (i) `SigmaAt` formula error (check against a finite-difference of the quantile map at three interior points), (ii) a sign slip in the pdf differences. Do not loosen the tolerance.

- [x] **Step 8: Full unit suite**

Run: `make test 2>&1 | tail -15`
Expected: `100% tests passed`. A failure outside dcdm_wdm means unexpected coupling — STOP and debug systematically.

- [x] **Step 9: Commit**

```bash
git add species/wdm_decay_product.h species/wdm_decay_product.cpp species/dcdm_wdm_test.cpp
git commit -m "dcdm_wdm: smooth erf energy deposit replaces CIC; ndf15 restored"
```

---

## Task 2: Hang-regression gate (the exact sweep that spun 8h)

**Files:**
- Create: `/private/tmp/claude-497981304/-Users-au192734-Projects-class-claude/36871189-16c5-48c8-9727-5a7ec35e3a61/scratchpad/hang_gate.py` (scratch — not committed)

**Interfaces:**
- Consumes: the branch build of `classy` (Task 1 tree). Install with `pip install . 2>&1 | tail -3` from the repo root (takes ~5–10 min).
- Produces: a timing table proving every previously-hanging configuration completes; gates Task 3.

- [x] **Step 1: Install the branch classy**

Run (repo root): `pip install . 2>&1 | tail -3`
Expected: `Successfully installed classy-community-...`

- [x] **Step 2: Write the gate script** — the exact sweep that hung on 2026-07-23, made hang-proof (subprocess + timeout, incremental output):

```python
"""Hang-regression gate: the exact (model x bins) sweep that spun 8h on the CIC
build must now complete. Each run in a subprocess with a hard timeout."""
import json, subprocess, sys, time

RUN_SRC = r'''
import sys, time
from classy import Class
ig, vk, bins, evolver = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
GYR = 977.792
c = Class()
c.set({"h": 0.6736, "omega_b": 0.02237, "output": "tCl,pCl,lCl,mPk", "lensing": "yes",
       "l_max_scalars": 2500, "P_k_max_h/Mpc": 1.0, "gauge": "synchronous",
       "evolver": evolver, "Omega_cdm": 0.26 * 0.9, "ddm.type": "dcdm_wdm",
       "ddm.Gamma": GYR / ig, "ddm.vkick": vk, "ddm.Omega_ini": 0.26 * 0.1,
       "ddm.momenta_bins": bins})
t0 = time.time()
c.compute()
print(f"OK {time.time() - t0:.1f}")
'''

TIMEOUT = 900  # 15 min per run — the CIC build blew this by >30x
results = {}
for ig, vk, tag in [(1.0, 0.011, "slow"), (0.1, 0.03, "early")]:
    for bins in [32, 64, 96, 128]:
        for evolver in ["2", "1"]:
            if evolver == "1" and bins > 64:
                continue  # ndf15 is ~14x slower; two bin counts suffice as proof
            key = f"{tag}_b{bins}_ev{evolver}"
            t0 = time.time()
            try:
                out = subprocess.run([sys.executable, "-c", RUN_SRC,
                                      str(ig), str(vk), str(bins), evolver],
                                     capture_output=True, text=True, timeout=TIMEOUT)
                ok = out.returncode == 0 and "OK" in out.stdout
                results[key] = {"ok": ok, "wall_s": round(time.time() - t0, 1),
                                "tail": (out.stdout + out.stderr)[-300:]}
            except subprocess.TimeoutExpired:
                results[key] = {"ok": False, "wall_s": TIMEOUT, "tail": "TIMEOUT"}
            print(key, results[key]["ok"], results[key]["wall_s"], flush=True)
            with open("hang_gate_results.json", "w") as f:
                json.dump(results, f, indent=1)
bad = [k for k, v in results.items() if not v["ok"]]
print("FAILED:" if bad else "ALL PASSED", bad)
```

- [x] **Step 3: Run it** (from the scratchpad directory, backgrounded; total expected ≲ 60 min)

Run: `cd <scratchpad> && python3 hang_gate.py 2>&1 | tee hang_gate.log`
Expected: every line `... True <seconds>`, final `ALL PASSED []`. The `evolver=2` runs should be ~1–5 min each; `evolver=1` runs slower but finite.

Any `TIMEOUT` or failure: STOP, use superpowers:systematic-debugging — sample the stuck subprocess (`sample <pid> 5`) before killing it, and report the stack. Do not proceed to Task 3.

- [x] **Step 4: Record** — paste the timing table into the task report (it is evidence for the PR description). Nothing to commit.

---

## Task 3: Old-vs-new convergence comparison (the failed deliverable)

**Files:**
- Modify: `notebooks/dcdm_wdm_convergence.py` (hang-proof rewrite)
- Modify: `notebooks/dcdm_wdm_convergence_compare.py` (REF=192 + robustness)
- Create: `notebooks/dcdm_wdm_convergence.png` / `.pdf` (figure)
- Create (not committed): `notebooks/conv_new.npz`, `notebooks/conv_old.npz`

**Interfaces:**
- Consumes: branch `classy` (installed in Task 2); a master build installed as module **`classyref`**.
- Produces: `conv_{new,old}.npz` with keys `KK`, `ELL`, `BINS`, `{model}_b{bins}_{pk,cl,t}`; the convergence figure; the bins-for-target table.

- [x] **Step 1: Build classyref from the existing master worktree.** The worktree `.worktrees/master-oldgrid` is already checked out at master (`e7e522a7`) — verify, do NOT checkout anything:

Run: `git -C .worktrees/master-oldgrid status --porcelain --branch | head -3`
Expected: `## master...` and no modified files. If dirty or not on master: STOP and report (do not fix with checkout/reset).

Find the reference-build recipe: `grep -rn "classyref\|CLASS_PYTHON_MODULE_NAME" Makefile .github/workflows/ | head`. Use the recipe found there; the known-good form is:

Run (inside `.worktrees/master-oldgrid`): `pip install . --config-settings=cmake.define.CLASS_PYTHON_MODULE_NAME=classyref 2>&1 | tail -3`
(plus whatever project-name switch the recipe specifies so the wheel installs alongside `classy`).
Verify: `python3 -c "import classyref, classy; print('both import OK')"`
Expected: `both import OK`. Also verify the module identities: `python3 -c "import classyref; c=classyref.Class(); c.set({'ddm.type':'dcdm_wdm','ddm.Gamma':97.8,'ddm.vkick':0.01,'ddm.Omega_ini':0.03,'ddm.kernel_width':'1.0'}); print('classyref accepts kernel_width => master build')"` — master accepts `kernel_width`; the branch build does not.

- [x] **Step 2: Rewrite the driver hang-proof** — replace `notebooks/dcdm_wdm_convergence.py` with:

```python
"""Momentum-grid convergence sweep for dcdm_wdm: run once per classy module.

Usage:
    python dcdm_wdm_convergence.py <module> <tag> [out.npz]

<module> is the python module to drive ("classy" = this branch's
injection-adapted grid, "classyref" = master's uniform-ln-q grid). Each
(model x bins) run executes in a SUBPROCESS with a hard timeout, and the
output .npz is rewritten after every run — a hang or crash costs one point,
never the sweep. Identical inputs on both modules (same evolver=2, same
k / l grids): an apples-to-apples accuracy-vs-bins test.
"""
import os
import subprocess
import sys
import time

import numpy as np

MODULE = sys.argv[1] if len(sys.argv) > 1 else "classy"
TAG = sys.argv[2] if len(sys.argv) > 2 else ("new" if MODULE == "classy" else "old")
OUT = sys.argv[3] if len(sys.argv) > 3 else f"conv_{TAG}.npz"
HERE = os.path.dirname(os.path.abspath(__file__))
TIMEOUT = 2400  # s per run; the 192-bin reference is the slowest

KK = np.logspace(-3, 0, 300)  # h/Mpc
ELL = np.arange(2, 2001)
MODELS = {
    "g0p1_v0p03": (0.1, 0.03),   # Gamma^-1 = 0.1 Gyr, v = 0.03  (early, hard)
    "g1_v0p011":  (1.0, 0.011),  # Gamma^-1 = 1.0 Gyr, v = 0.011 (notebook default)
}
BINS = [16, 24, 32, 48, 96, 192]  # 192 = reference ("truth")

RUN_SRC = r'''
import importlib, sys, time
import numpy as np
module, ig, vk, bins, out = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
cls = importlib.import_module(module).Class
GYR = 977.792
c = cls()
c.set({"h": 0.6736, "omega_b": 0.02237, "output": "tCl,pCl,lCl,mPk", "lensing": "yes",
       "l_max_scalars": 2500, "P_k_max_h/Mpc": 1.0, "gauge": "synchronous", "evolver": "2",
       "Omega_cdm": 0.26 * 0.9, "ddm.type": "dcdm_wdm", "ddm.Gamma": GYR / ig,
       "ddm.vkick": vk, "ddm.Omega_ini": 0.26 * 0.1, "ddm.momenta_bins": bins})
t0 = time.time()
c.compute()
t = time.time() - t0
h = c.h()
kk = np.logspace(-3, 0, 300)
pk = np.array([c.pk(k * h, 0.0) for k in kk])
cl = np.array(c.lensed_cl(2000)["tt"][2:2001])
np.savez(out, pk=pk, cl=cl, t=t)
print(f"OK {t:.1f}")
'''


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] [{TAG}] {msg}", flush=True)


data = {"KK": KK, "ELL": ELL, "BINS": np.array(BINS), "TAG": np.array(TAG)}
out_path = os.path.join(HERE, OUT)
for name, (ig, vk) in MODELS.items():
    for b in BINS:
        shard = os.path.join(HERE, f"_shard_{TAG}_{name}_{b}.npz")
        t0 = time.time()
        try:
            r = subprocess.run([sys.executable, "-c", RUN_SRC, MODULE,
                                str(ig), str(vk), str(b), shard],
                               capture_output=True, text=True, timeout=TIMEOUT)
            ok = r.returncode == 0 and os.path.exists(shard)
        except subprocess.TimeoutExpired:
            ok = False
            log(f"{name} bins={b}: TIMEOUT after {TIMEOUT}s")
        if ok:
            s = np.load(shard)
            data[f"{name}_b{b}_pk"] = s["pk"]
            data[f"{name}_b{b}_cl"] = s["cl"]
            data[f"{name}_b{b}_t"] = float(s["t"])
            os.remove(shard)
            log(f"{name} bins={b:4d}: {data[f'{name}_b{b}_t']:6.1f} s compute, "
                f"{time.time() - t0:6.1f} s wall")
        else:
            data[f"{name}_b{b}_pk"] = np.full_like(KK, np.nan)
            data[f"{name}_b{b}_cl"] = np.full(len(ELL), np.nan)
            data[f"{name}_b{b}_t"] = np.nan
            if 'r' in dir():
                log(f"{name} bins={b}: FAILED -- {(r.stdout + r.stderr)[-300:]}")
        np.savez(out_path, **data)  # incremental: saved after EVERY run
log(f"saved -> {out_path}")
```

- [x] **Step 3: Run both sweeps** (backgrounded, sequential; expect ~1–3 h total):

Run (in `notebooks/`): `python3 dcdm_wdm_convergence.py classy new 2>&1 | tee conv_new.log`
then: `python3 dcdm_wdm_convergence.py classyref old 2>&1 | tee conv_old.log`
Expected: 12 timestamped `OK` lines each, `conv_new.npz` / `conv_old.npz` written. A TIMEOUT on any single point: note it, continue (NaN is handled), but a TIMEOUT on a *new*-build point reopens the Task-2 gate — STOP and debug.

- [x] **Step 4: Update the compare script** — in `notebooks/dcdm_wdm_convergence_compare.py`: change `REF = 256` to `REF = 192`. Before plotting, load the `dataviz` skill and restyle the figure to its standards. Keep: the old-192 vs new-192 consistency print (this validates the shared continuum limit — expect ≲ 1e-2; if larger, STOP and investigate before trusting any conclusion), the per-model error curves, and the bins-for-target table for targets `1% / 0.3% / 0.1%`.

- [x] **Step 5: Run the comparison**

Run: `python3 dcdm_wdm_convergence_compare.py`
Expected: consistency `≲ 1e-2` per model; figure saved; a table like

```
model            target  old bins  new bins  reduction
g0p1_v0p03         1.0%       ...       ...      ...x
```

Record the table verbatim in the task report. The headline claim (injection-adapted grid reaches fixed accuracy at materially fewer bins) must be supported by the numbers — if it is NOT (reduction ≲ 1.5x), report that honestly; do not massage.

- [x] **Step 6: Commit** (scripts + figure; npz/logs stay untracked)

```bash
git add notebooks/dcdm_wdm_convergence.py notebooks/dcdm_wdm_convergence_compare.py notebooks/dcdm_wdm_convergence.png notebooks/dcdm_wdm_convergence.pdf
git commit -m "dcdm_wdm: hang-proof convergence harness + old-vs-new grid comparison"
```

---

## Task 4: Notebook, Fig. 3, docs refresh + final green

**Files:**
- Modify: `notebooks/dcdm-wdm-massive-decay-products.ipynb` (§5 convergence, §6 limitations)
- Regenerate if shifted: `notebooks/dcdm_wdm_fig3_cache.npz`, `notebooks/dcdm_wdm_fig3_9models.{png,pdf}`
- Modify: `docs/superpowers/specs/2026-07-18-dcdm-wdm-adaptive-momentum-grid-design.md` (status header), this plan (checkboxes)
- Modify: `species/dcdm_wdm_test.cpp` (stale comment only, if not already clean)

**Interfaces:** verification + documentation only; no new code interfaces.

- [x] **Step 1: Scrub stale wording** — grep `species/wdm_decay_product.{h,cpp}` and `species/dcdm_wdm_test.cpp` for `CIC`, `CicWeights`, `kernel_width`, `LocalDu`, `du_bin`, `evolver=2 (RKDP45) required`, `hangs the stiff`: no stale references may remain in code or comments (historical mentions in the spec/plan docs are fine). In the test file, replace the long `DEVIATION from the brief: momenta_bins=32 ...` comment (perturbed smoke block) with:

```cpp
  // Perturbed run smoke test: the full pipeline runs and Cls are finite.
  // (Uses momenta_bins=48 to stay fast; the historical momenta_bins=32
  // -ffast-math background failure no longer reproduces on master.)
```

- [x] **Step 2: Notebook §5** — in `notebooks/dcdm-wdm-massive-decay-products.ipynb`: update the convergence cell to sweep `(16, 24, 32, 48, 96)` against a 192 reference, re-execute, and rewrite the §5 markdown: the grid is injection-adapted (equal-decay F-quantiles via a fiducial H(t), grid = f(Γ) only), the deposit is a conservative cell-integrated Gaussian (erf) with σ = clamp(local bin width, 0.03, 0.25), `q_edge_tol` is the grid-span knob (`kernel_width` is gone), and both evolvers work (`evolver=2` recommended for speed, ~14×). Quote the Task-3 bins-for-target numbers.

- [x] **Step 3: Notebook §6** — remove/rewrite the stale bullets: the `momenta_bins=32 -ffast-math` bullet (resolved on master) and the "very short lifetimes need lower q_min_ratio" bullet (the window auto-follows the decay now).

- [x] **Step 4: Fig. 3 driver** — run `python3 notebooks/dcdm_wdm_fig3_9models.py`. The converged physics must be unchanged: compare regenerated curves against the committed cache with a scale-relative metric; shifts must be ≲ 0.3% (discretization improvement, not new physics). If larger: STOP, systematic-debugging. Commit regenerated cache/figures only if they actually shifted beyond noise.

- [x] **Step 5: Full suite + docs**

Run: `make test 2>&1 | tail -10` → `100% tests passed`.
Update the spec's top Status block: v3 erf deposit BUILT (commit hash), hang-regression + comparison results one-liner. Tick this plan's checkboxes.

- [x] **Step 6: Commit**

```bash
git add notebooks/dcdm-wdm-massive-decay-products.ipynb docs/superpowers/specs/2026-07-18-dcdm-wdm-adaptive-momentum-grid-design.md docs/superpowers/plans/2026-07-18-dcdm-wdm-injection-adapted-grid.md species/dcdm_wdm_test.cpp
# plus notebooks/dcdm_wdm_fig3_* only if regenerated beyond noise
git commit -m "dcdm_wdm: notebook + docs refresh for the erf-deposit quantile grid"
```

---

## Self-Review Notes (author)

- **Spec coverage:** §4.3 v3 width/deposit/dJdlnq → Task 1 Steps 3–6 (exact code); §4.4 evolver-guard removal → Task 1 Steps 1(b)/4(a); §6.2 golden under ndf15 → Task 1 Step 7 (BaseFc no longer forces an evolver); §6.5 deposit tests → Task 1 Step 1(c–f); §6.6 hang-regression → Task 1 Step 7 + Task 2 (field-scale); §6 comparison harness → Task 3 (subprocess+timeout+incremental, classyref, consistency gate, bins-for-target); §6 notebook → Task 4.
- **Type consistency:** `SigmaAt(double) const` public (test calls it); `GaussWeights(double, double, double*) const` private; `u_edge()` accessor added (test 1(f) uses it); `F_lo_`/`F_hi_` set in both grid branches; `scratch_w_` reused as weight buffer.
- **Load-bearing numerics:** exact sum rule via telescoping Σw=1; exp clamp 60 (pdf) and erf clamp 8 (cdf); log-form Δu_loc (no inf); `kFSeed` untouched.
- **Placeholder scan:** none — all code steps complete.
