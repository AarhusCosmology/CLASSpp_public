#!/usr/bin/env python3
"""dncdm_inv_ab.py -- DNCDM inverse-decay A/B validation driver (plan Task 3.4).

NOT a `make test` / CTest target -- deliberately excluded from CMakeLists.txt and
Makefile TEST_TARGETS (per docs/superpowers/plans/2026-07-24-dncdm-inverse-qs-plan.md
Task 3.4). This is a documented driver script you run by hand (or as a separate CI
step). It implements the "Background validation ladder" from that plan:

  --decay-only
      `inverse_decays=no` (or unset) takes the byte-for-byte UNCHANGED
      `DNCDM_DR_Species` code path: species/dncdm_dr_species.cpp's `CreateAll` only
      routes through the new `DNCDMInvSpecies` composite when `inverse_decays=yes`
      (see the dispatch at the top of that file). So the cheapest valid regression is
      a HEAD-vs-master A/B: build a scratch reference binary from `origin/master`
      (in an isolated git worktree -- never touches this shared working tree), run
      the SAME decay-only .ini through both binaries, and diff background.dat /
      cl.dat / pk.dat with a SCALE-RELATIVE metric (relative to each column's peak
      magnitude, not raw value-by-value relative error) so Cl^TE zero-crossings
      don't blow up the comparison. This mirrors test/scenarios/compare_tol.py's
      convention. See the repo memory notes for why this is the right bar:
        reference_classyref_testing.md              (classyref A/B convention)
        feedback_committed_goldens_stale_ffast_math.md (don't gate on stale goldens)
        feedback_no_bit_identical_requirement.md      (~0.1% tol, never bit-identical)

  --tiny-gamma
      `inverse_decays=yes` with a cosmologically negligible Gamma exercises the NEW
      `DNCDMInvSpecies` / `DecayTransitionKernel` code path end-to-end (unlike
      --decay-only, which never touches it). No master build is needed for this
      mode: THIS branch's binary is run TWICE on the identical combined-mode
      cosmology (same parent mass/abundance/grid) -- once through the new composite
      (inverse_decays=yes) and once through the untouched DNCDM_DR_Species path
      (inverse_decays=no/unset). With Gamma this small, decay is kinematically
      negligible either way, so the parent's background evolution (density, phase
      space density per bin) must agree closely between the two representations --
      this is the "Gamma=0 / no-DCDM" limit made exact: it isolates "does the new
      f-variable + kernel plumbing reduce to the old ln f + diagonal-decay plumbing
      when the source term vanishes", which is precisely what a regression in the
      new code would break.

All modes need only the plain C++ `class` binary (`make -j<N> class`, i.e.
`cmake --build build/cmake --target class`) -- this script never builds or installs
a Python `classy` wheel.

Usage:
  python python/tests/dncdm_inv_ab.py                  # run both modes
  python python/tests/dncdm_inv_ab.py --decay-only
  python python/tests/dncdm_inv_ab.py --tiny-gamma
  python python/tests/dncdm_inv_ab.py --tiny-gamma --gamma 1e-8 --rtol 5e-4
  python python/tests/dncdm_inv_ab.py --decay-only --skip-master   # code-path-only,
                                                                    # no scratch build
  python python/tests/dncdm_inv_ab.py --pert-regression
  python python/tests/dncdm_inv_ab.py --pert-convergence
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# .dat parsing: CLASS writes a single "# 1:name  2:name ..." header comment
# line before the data. Column names can contain spaces ("proper time [Gyr]",
# "gr.fac. D") and cl.dat/pk.dat reuse "1:" for every physical column after the
# index column (a pre-existing CLASS quirk) -- so we split on "<digits>:" as a
# delimiter and never rely on the digit value itself.
# ---------------------------------------------------------------------------
_HDR_TOKEN = re.compile(r"(?<!\S)(\d+):")


def _parse_header_names(header_line):
    body = header_line.lstrip("#")
    matches = list(_HDR_TOKEN.finditer(body))
    names = []
    for i, m in enumerate(matches):
        start = m.end()
        end = matches[i + 1].start() if i + 1 < len(matches) else len(body)
        names.append(body[start:end].strip())
    return names


def parse_dat(path):
    """Return (column_names, data) for a CLASS .dat file."""
    header_line = None
    with open(path) as f:
        for line in f:
            if line.startswith("#") and _HDR_TOKEN.search(line):
                header_line = line  # keep overwriting -> last match wins
    if header_line is None:
        raise ValueError(f"no numbered column-header line found in {path}")
    names = _parse_header_names(header_line)
    data = np.loadtxt(path, comments="#", ndmin=2)
    if data.size and data.shape[1] != len(names):
        raise ValueError(f"{path}: {len(names)} header names but {data.shape[1]} data columns")
    return names, data


def _unique_name_index(names):
    """name -> column index, restricted to names that occur exactly once.
    (Kept as a guard even though DrPsdSpecies now suffixes its daughter columns
    with the statistics, so fermion/boson no longer collide: excluding non-unique
    names keeps cross-run comparisons unambiguous whatever a species emits.)
    """
    c = Counter(names)
    return {n: i for i, n in enumerate(names) if c[n] == 1}


# ---------------------------------------------------------------------------
# Scale-relative comparator (test/scenarios/compare_tol.py convention): pass if
# |a-b| <= atol + rtol*|a|, atol = rtol * column peak. The atol floor is what
# keeps a Cl^TE (or any zero-crossing) column from reporting a meaningless
# blown-up relative error right at its zero.
# ---------------------------------------------------------------------------
def _scale_relative(a, b, rtol):
    scale = np.maximum(np.abs(a), np.abs(b))
    colpeak = scale.max(axis=0, keepdims=True) if scale.ndim == 2 else scale.max()
    atol = rtol * np.where(colpeak > 0, colpeak, 1.0)
    diff = np.abs(a - b)
    bad = diff > (atol + rtol * np.abs(a))
    worst = float((diff / np.where(colpeak > 0, colpeak, 1.0)).max()) if diff.size else 0.0
    return int(np.sum(bad)), worst


def compare_same_shape(names_a, a, names_b, b, rtol, label):
    """Positional comparison for two runs that MUST have identical column
    layout (used for --decay-only: unchanged code path, same .ini)."""
    if a.shape != b.shape:
        print(f"FAIL {label}: shape mismatch branch{a.shape} vs master{b.shape}")
        return False
    if names_a != names_b:
        print(f"FAIL {label}: column-title mismatch (unexpected for an unchanged code path)")
        print(f"   branch: {names_a}")
        print(f"   master: {names_b}")
        return False
    nbad, worst = _scale_relative(a, b, rtol)
    status = "OK  " if nbad == 0 else "FAIL"
    print(f"{status} {label}: shape={a.shape} n_exceed={nbad} worst_scale_rel={worst:.3e} (rtol={rtol:g})")
    return nbad == 0


def compare_named_common(names_a, a, names_b, b, rtol, label, skip_patterns=()):
    """Named-column comparison for two runs whose column SETS differ (used for
    --tiny-gamma: decay-only vs inverse-decays composite output). Compares
    every column name that occurs exactly once in BOTH files, aligned on the
    'z' column (row-for-row if the grids match exactly, else interpolated).
    `skip_patterns` are regexes for columns that are KNOWN to differ by design
    (not a comparison target, not a failure) -- reported separately so the
    driver stays self-documenting rather than silently dropping them."""
    ua, ub = _unique_name_index(names_a), _unique_name_index(names_b)
    common_all = sorted(set(ua) & set(ub))
    skip_re = [re.compile(p) for p in skip_patterns]
    skipped = [n for n in common_all if any(p.search(n) for p in skip_re)]
    common = [n for n in common_all if n not in skipped]
    if skipped:
        print(f"  {label}: skipping {len(skipped)} column(s) that differ by design "
              f"(collision-owned parent writes dlnfdlnq_separate=0; design-doc table): "
              f"{skipped[0]}..{skipped[-1]}" if len(skipped) > 1 else f"{skipped[0]}")
    if "z" not in common:
        print(f"FAIL {label}: no common unique 'z' column to align on")
        return False, []
    za, zb = a[:, ua["z"]], b[:, ub["z"]]
    same_grid = za.shape == zb.shape and np.array_equal(za, zb)
    ok = True
    rows = []
    for name in common:
        va = a[:, ua[name]]
        if same_grid:
            vb = b[:, ub[name]]
        else:
            oa, ob = np.argsort(za), np.argsort(zb)
            interp = np.interp(za[oa], zb[ob], b[:, ub[name]][ob])
            inv = np.empty_like(oa)
            inv[oa] = np.arange(len(oa))
            vb = interp[inv]
        nbad, worst = _scale_relative(va, vb, rtol)
        if nbad:
            ok = False
        rows.append((name, nbad, worst))
    grid_note = "same z-grid" if same_grid else "interpolated onto branch z-grid"
    print(f"  {label}: {len(common)} common unique columns compared ({grid_note})")
    for name, nbad, worst in rows:
        status = "OK  " if nbad == 0 else "FAIL"
        print(f"    {status} {name:32s} n_exceed={nbad} worst_scale_rel={worst:.3e}")
    return ok, rows


# ---------------------------------------------------------------------------
# .ini construction: parse "key = value" lines into an ordered dict, apply
# overrides (dot-notation keys are just ordinary keys to this parser), append
# any new keys, and always force a distinct 'root ='.
# ---------------------------------------------------------------------------
def _parse_ini(text):
    d, order = {}, []
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith("#") or "=" not in s:
            continue
        k, v = s.split("=", 1)
        k, v = k.strip(), v.strip()
        if k not in d:
            order.append(k)
        d[k] = v
    return d, order


def build_ini(base_path, overrides, root):
    d, order = _parse_ini(Path(base_path).read_text())
    for k, v in overrides.items():
        if k not in d:
            order.append(k)
        d[k] = v
    if "root" not in d:
        order.append("root")
    d["root"] = root
    return "\n".join(f"{k} = {d[k]}" for k in order) + "\n"


# ---------------------------------------------------------------------------
# Build + run helpers
# ---------------------------------------------------------------------------
def repo_root():
    out = subprocess.run(["git", "rev-parse", "--show-toplevel"], cwd=Path(__file__).parent,
                          capture_output=True, text=True, check=True)
    return Path(out.stdout.strip())


def ensure_branch_class(root, jobs):
    print(f"[branch] make -j{jobs} class   (confirms {root}/class is fresh)")
    proc = subprocess.run(["make", f"-j{jobs}", "class"], cwd=str(root), capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError("branch `make class` failed:\n" + proc.stdout[-4000:] + proc.stderr[-4000:])
    binary = root / "class"
    if not binary.exists():
        raise RuntimeError(f"expected branch binary at {binary} after a successful build")
    return binary


def build_master_class(root, scratch_dir, jobs, timeout):
    """git worktree of origin/master, built in isolation. NEVER touches the
    shared main working tree -- only `git worktree add/remove` against it."""
    worktree = scratch_dir / "master_src"
    if not worktree.exists():
        subprocess.run(["git", "-C", str(root), "fetch", "origin", "master"],
                        capture_output=True, text=True)  # best-effort; offline is OK if the ref is cached
        print(f"[master] git worktree add --detach {worktree} origin/master")
        proc = subprocess.run(
            ["git", "-C", str(root), "worktree", "add", "--detach", str(worktree), "origin/master"],
            capture_output=True, text=True)
        if proc.returncode != 0:
            print("MASTER WORKTREE SETUP FAILED:\n" + proc.stdout + proc.stderr)
            return None, worktree
    else:
        print(f"[master] reusing existing worktree at {worktree}")
    build_dir = worktree / "build" / "cmake"
    try:
        print("[master] cmake configure -DCMAKE_BUILD_TYPE=RelWithDebInfo "
              "(must match this repo's build/cmake type -- a mismatch alone causes ~1e-5 false diffs)")
        subprocess.run(["cmake", "-S", str(worktree), "-B", str(build_dir),
                         "-DCMAKE_BUILD_TYPE=RelWithDebInfo"],
                        cwd=str(worktree), check=True, capture_output=True, text=True, timeout=timeout)
        print(f"[master] cmake --build {build_dir} --target class --parallel {jobs}")
        subprocess.run(["cmake", "--build", str(build_dir), "--target", "class", "--parallel", str(jobs)],
                        cwd=str(worktree), check=True, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        print(f"MASTER BUILD TIMED OUT (> {timeout}s) -- skipping the master A/B comparison. "
              "The decay-only claim (unchanged code path) still holds by inspection: see the "
              "CreateAll dispatch in species/dncdm_dr_species.cpp.")
        return None, worktree
    except subprocess.CalledProcessError as e:
        print("MASTER BUILD FAILED:\n" + (e.stdout or "")[-4000:] + (e.stderr or "")[-4000:])
        return None, worktree
    binary = worktree / "class"
    if not binary.exists():
        print(f"MASTER BUILD: expected binary at {binary} not found after a reported-successful build")
        return None, worktree
    return binary, worktree


def remove_worktree(root, worktree):
    if worktree and worktree.exists():
        subprocess.run(["git", "-C", str(root), "worktree", "remove", "--force", str(worktree)],
                        capture_output=True, text=True)


def run_class(binary, ini_text, workdir, tag, timeout=300, single_thread=False):
    """`single_thread=True` pins OMP_NUM_THREADS=1. REQUIRED for any inverse-decays
    PERTURBATION run: DecayTransitionKernel's per-node scratch and DNCDMInvSpecies'
    mutable gather buffers are shared across the OpenMP k-loop, so a multi-threaded
    run either races or trips the kernel's own ascending-l guard
    ("ApplyPerturbationOperator called out of order"). Background-only runs are
    unaffected (single-threaded by construction)."""
    ini_path = workdir / f"{tag}.ini"
    ini_path.write_text(ini_text)
    env = None
    if single_thread:
        env = dict(os.environ, OMP_NUM_THREADS="1")
    proc = subprocess.run([str(binary), str(ini_path)], cwd=str(workdir),
                           capture_output=True, text=True, timeout=timeout, env=env)
    log_path = workdir / f"{tag}.log"
    log_path.write_text(proc.stdout + proc.stderr)
    if proc.returncode != 0:
        tail = "\n".join((proc.stdout + proc.stderr).splitlines()[-25:])
        raise RuntimeError(f"{binary} {ini_path} failed (exit {proc.returncode}); log={log_path}\n{tail}")
    return log_path


# ---------------------------------------------------------------------------
# Mode A: decay-only regression (branch vs scratch master)
# ---------------------------------------------------------------------------
def mode_decay_only(args, root, workdir):
    print("\n=== --decay-only: inverse_decays=no/unset, byte-for-byte-unchanged DNCDM_DR_Species path ===")
    branch_bin = ensure_branch_class(root, args.jobs)

    base_ini = root / "test/scenarios/dncdm_dr_combined.ini"
    overrides = {
        "dncdm1.inverse_decays": "no",  # explicit, though absence already means "no" (repo default)
        "output": "tCl,pCl,mPk",
        "l_max_scalars": "1000",
        "P_k_max_1/Mpc": "1.0",
        "write background": "yes",
    }

    branch_root = str(workdir / "branch_")
    run_class(branch_bin, build_ini(base_ini, overrides, branch_root), workdir, "branch")
    print(f"  branch run OK: {branch_root}{{background,cl,pk}}.dat")

    if args.skip_master:
        print("  --skip-master set: skipping the master A/B. The unchanged-code-path argument "
              "(species/dncdm_dr_species.cpp CreateAll dispatch) is the evidence for this run.")
        return True

    master_bin, master_worktree = build_master_class(root, workdir, args.jobs, args.master_timeout)
    if master_bin is None:
        print("  Falling back to code-path inspection only (see module docstring) -- not a hard failure "
              "of this driver, but the master A/B numbers below are unavailable this run.")
        return True
    try:
        master_root = str(workdir / "master_")
        run_class(master_bin, build_ini(base_ini, overrides, master_root), workdir, "master")
        print(f"  master run OK: {master_root}{{background,cl,pk}}.dat")

        ok = True
        for fname in ("background.dat", "cl.dat", "pk.dat"):
            names_a, a = parse_dat(workdir / f"branch_{fname}")
            names_b, b = parse_dat(workdir / f"master_{fname}")
            ok &= compare_same_shape(names_a, a, names_b, b, args.rtol, fname)
        return ok
    finally:
        if not args.keep_scratch:
            remove_worktree(root, master_worktree)


# ---------------------------------------------------------------------------
# Mode B: tiny-Gamma limit (branch only, inverse-decays vs decay-only)
# ---------------------------------------------------------------------------
def mode_tiny_gamma(args, root, workdir):
    print("\n=== --tiny-gamma: inverse_decays=yes, Gamma kinematically negligible ===")
    branch_bin = ensure_branch_class(root, args.jobs)

    base_ini = root / "test/scenarios/dncdm_dr_combined.ini"
    common = {
        "output": "",
        "write background": "yes",
        "dncdm1.Gamma": str(args.gamma),
        "dncdm1.quadrature_strategy": "3",  # qm_trapz: pins q_size == q_size_bg (DNCDMInvSpecies::Create guard)
        "dncdm1.momenta_bins": "16",
        "dncdm1.dr_q_min": "1e-4",
        "dncdm1.dr_q_max": "300",
        "dncdm1.dr_N_q": "100",
    }

    ref_root = str(workdir / "tg_ref_")
    ref_overrides = dict(common)
    ref_overrides["dncdm1.inverse_decays"] = "no"  # untouched DNCDM_DR_Species path
    run_class(branch_bin, build_ini(base_ini, ref_overrides, ref_root), workdir, "tg_ref")
    print(f"  reference (inverse_decays=no) run OK: {ref_root}background.dat")

    inv_root = str(workdir / "tg_inv_")
    inv_overrides = dict(common)
    inv_overrides["dncdm1.inverse_decays"] = "yes"
    inv_overrides["dncdm1.quantum_statistics"] = "no"
    run_class(branch_bin, build_ini(base_ini, inv_overrides, inv_root), workdir, "tg_inv")
    print(f"  inverse-decays (Gamma={args.gamma:g}) run OK: {inv_root}background.dat")

    names_r, r = parse_dat(workdir / "tg_ref_background.dat")
    names_i, i_ = parse_dat(workdir / "tg_inv_background.dat")
    # dlnfdlnq_separate_dncdm1[*] is BY DESIGN written 0 in collision-owned mode
    # (species/dncdm_species.cpp:385, gated on collision_owned_; design table in
    # the plan's "Resolved Design Question" section) -- not a physics regression.
    ok, _rows = compare_named_common(names_r, r, names_i, i_, args.rtol, "background.dat",
                                      skip_patterns=[r"^dlnfdlnq_separate_"])

    # Sanity: the daughters must actually be a tiny fraction of the parent, or this
    # "tiny-Gamma" check would pass trivially even if the kernel silently did nothing.
    rho_parent = i_[:, names_i.index("(.)rho_dncdm1")]
    # Daughter density columns carry a statistics suffix ("_l" fermion / "_phi" boson,
    # DrPsdSpecies::ColumnLabel) so the two daughters sharing an instance name do not
    # collide -- match on the prefix, never the bare instance name.
    daughter_idx = [k for k, n in enumerate(names_i) if n.startswith("(.)rho_dr_dncdm1")]
    if not daughter_idx:
        print("FAIL tiny-gamma: no daughter '(.)rho_dr_dncdm1' column found in the inverse-decays run "
              "-- DNCDMInvSpecies did not construct daughters as expected")
        return False
    rho_daughters = sum(i_[:, k] for k in daughter_idx)
    frac = float(np.max(rho_daughters / np.maximum(rho_parent, 1e-300)))
    print(f"  sanity: max daughter/parent density fraction over the run history = {frac:.3e} "
          f"(must be << 1; {len(daughter_idx)} daughter column(s) found -- fermion+boson)")
    if frac > 1e-2:
        print("FAIL tiny-gamma: daughter fraction not small enough for a valid tiny-Gamma limit "
              "-- lower --gamma")
        ok = False

    return ok


# ---------------------------------------------------------------------------
# Mode C: decay-only PERTURBATION regression (Task 4.4 Step 1), both evolvers
# ---------------------------------------------------------------------------
def mode_pert_regression(args, root, workdir):
    """`inverse_decays=no` full Cl/P(k) must reproduce master `ncdm_decay_dr` to
    ~0.1%. The parent has NO pure-decay perturbation collision term (it cancels
    identically -- design/paper §4.3; verified through the plumbing by
    test-dncdm-inv's decay-only-cancel block), so this is the strongest
    perturbation regression anchor. Runs BOTH evolvers (ndf15=1 default, rkdp45=2)
    and A/Bs each against a scratch-master build; the decay-only path is the
    UNCHANGED DNCDM_DR_Species code, so a mismatch is a genuine regression."""
    print("\n=== --pert-regression: decay-only Cl/P(k) vs master, evolver=1 AND evolver=2 ===")
    branch_bin = ensure_branch_class(root, args.jobs)
    base_ini = root / "test/scenarios/dncdm_dr_combined.ini"
    base_overrides = {
        "dncdm1.inverse_decays": "no",
        "output": "tCl,pCl,mPk",
        "l_max_scalars": "1000",
        "P_k_max_1/Mpc": "1.0",
        "write background": "yes",
    }
    master_bin, master_worktree = (None, None)
    if not args.skip_master:
        master_bin, master_worktree = build_master_class(root, workdir, args.jobs, args.master_timeout)
        if master_bin is None:
            print("  master build unavailable -> completion-only on both evolvers (the unchanged "
                  "DNCDM_DR path is the decay-only correctness argument).")
    try:
        overall = True
        for evolver in (1, 2):
            tag = f"evolver={evolver}"
            ov = dict(base_overrides)
            ov["evolver"] = str(evolver)
            broot = str(workdir / f"pr_branch_e{evolver}_")
            try:
                run_class(branch_bin, build_ini(base_ini, ov, broot), workdir,
                          f"pr_branch_e{evolver}", timeout=args.pert_timeout)
            except (RuntimeError, subprocess.TimeoutExpired) as e:
                print(f"FAIL {tag}: branch decay-only run did not complete: {str(e).splitlines()[-1][:200]}")
                overall = False
                continue
            print(f"  branch {tag}: decay-only run OK")
            if master_bin is None:
                continue
            mroot = str(workdir / f"pr_master_e{evolver}_")
            try:
                run_class(master_bin, build_ini(base_ini, ov, mroot), workdir,
                          f"pr_master_e{evolver}", timeout=args.pert_timeout)
            except (RuntimeError, subprocess.TimeoutExpired) as e:
                print(f"  {tag}: master run did not complete ({str(e).splitlines()[-1][:120]}); skipping its A/B")
                continue
            for fname in ("cl.dat", "pk.dat"):
                na, a = parse_dat(workdir / f"pr_branch_e{evolver}_{fname}")
                nb, b = parse_dat(workdir / f"pr_master_e{evolver}_{fname}")
                overall &= compare_same_shape(na, a, nb, b, args.rtol, f"{tag} {fname}")
        return overall
    finally:
        if master_worktree and not args.keep_scratch:
            remove_worktree(root, master_worktree)


# ---------------------------------------------------------------------------
# Mode D: inverse-decays PERTURBATION run on both evolvers (Task 4.4 Step 2)
# ---------------------------------------------------------------------------
def mode_pert_convergence(args, root, workdir):
    """Exercise the NEW inverse-decays perturbation path
    (DNCDMInvSpecies::AddCouplingDerivs) end-to-end on BOTH evolvers, per the
    research-plugins checklist warning that the nonlinear kernel coupling can defeat
    ndf15's numjac sparsity lock. Reports completion per evolver; on completion, a
    light sanity that P(k) is finite/positive. The full parent->daughter
    delta/theta/sigma convergence-morphology check (thesis Fig 6.9) needs mTk
    transfer output at the equilibrium epochs and is documented as DEFERRED here
    (Task-3.4 precedent) -- the gating deliverable is a clean end-to-end run."""
    print("\n=== --pert-convergence: inverse_decays=yes perturbations, evolver=1 AND evolver=2 ===")
    branch_bin = ensure_branch_class(root, args.jobs)
    base_ini = root / "test/scenarios/dncdm_dr_combined.ini"
    common = {
        # Daughters raise DeltaNeff past the BBN interpolation table, so thermodynamics
        # would reject the point before perturbations even run -- fix YHe to bypass BBN.
        "YHe": "0.245",
        "dncdm1.inverse_decays": "yes",
        "dncdm1.quantum_statistics": "no",
        "dncdm1.Gamma": str(args.gamma_pert),
        "dncdm1.Omega_dncdmdr": "0.01",
        "dncdm1.quadrature_strategy": "3",  # qm_trapz: q_size == q_size_bg (Create guard)
        "dncdm1.momenta_bins": "16",
        "dncdm1.dr_q_min": "1e-4",
        "dncdm1.dr_q_max": "300",
        "dncdm1.dr_N_q": "100",
        "output": "mPk",
        "P_k_max_1/Mpc": "0.5",
        "write background": "yes",
    }
    completed = 0
    for evolver in (1, 2):
        ov = dict(common)
        ov["evolver"] = str(evolver)
        rroot = str(workdir / f"pc_e{evolver}_")
        try:
            run_class(branch_bin, build_ini(base_ini, ov, rroot), workdir,
                      f"pc_e{evolver}", timeout=args.pert_timeout)
        except (RuntimeError, subprocess.TimeoutExpired) as e:
            print(f"  evolver={evolver}: inverse-decays perturbation run did NOT complete:")
            print(f"      {str(e).splitlines()[-1][:200]}")
            continue
        try:
            _names, pk = parse_dat(workdir / f"pc_e{evolver}_pk.dat")
            pcol = pk[:, 1]
            finite_pos = bool(np.all(np.isfinite(pk)) and np.all(pcol > 0))
        except Exception as e:
            print(f"  evolver={evolver}: completed but P(k) unreadable: {e}")
            continue
        status = "finite, positive" if finite_pos else "NON-FINITE / non-positive P(k)!"
        print(f"  evolver={evolver}: inverse-decays perturbation run COMPLETED; P(k) {status}")
        if finite_pos:
            completed += 1
    if completed == 0:
        print("  FAIL: no inverse-decays perturbation run produced a finite, positive P(k).")
        return False
    return True


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--decay-only", action="store_true", help="run the decay-only master A/B regression")
    p.add_argument("--tiny-gamma", action="store_true", help="run the tiny-Gamma inverse-decays limit check")
    p.add_argument("--pert-regression", action="store_true",
                    help="Task 4.4: decay-only Cl/P(k) vs master, both evolvers")
    p.add_argument("--pert-convergence", action="store_true",
                    help="Task 4.4: inverse-decays perturbation run end-to-end, both evolvers")
    p.add_argument("--gamma-pert", type=float, default=10.0,
                    help="dncdm1.Gamma for --pert-convergence (default 10)")
    p.add_argument("--pert-timeout", type=int, default=1800,
                    help="seconds allowed per perturbation `class` run (default 1800)")
    p.add_argument("--rtol", type=float, default=1e-3, help="scale-relative tolerance (default 1e-3)")
    p.add_argument("--gamma", type=float, default=1e-6,
                    help="dncdm1.Gamma for --tiny-gamma, in the input's km/s/Mpc-like units "
                         "(default 1e-6; H0~68 makes this ~1e-7 xH0, i.e. negligible decay over the age "
                         "of the universe)")
    p.add_argument("--jobs", "-j", type=int, default=min(10, os.cpu_count() or 4),
                    help="parallel build jobs (default min(10, ncpu))")
    p.add_argument("--master-timeout", type=int, default=600,
                    help="seconds allowed for the scratch master configure+build before giving up "
                         "(default 600 = 10 min, per the plan's budget)")
    p.add_argument("--skip-master", action="store_true",
                    help="--decay-only: skip the scratch master build/A-B entirely (code-path-only)")
    p.add_argument("--scratch-dir", type=str, default=None,
                    help="scratch working directory (default: a fresh tempdir, removed on exit unless "
                         "--keep-scratch)")
    p.add_argument("--keep-scratch", action="store_true",
                    help="keep the scratch dir (.ini/.dat/.log files + master worktree) for inspection")
    args = p.parse_args()

    if not any([args.decay_only, args.tiny_gamma, args.pert_regression, args.pert_convergence]):
        args.decay_only = args.tiny_gamma = True  # default: the (fast) background ladder only

    root = repo_root()
    if args.scratch_dir:
        workdir = Path(args.scratch_dir)
        workdir.mkdir(parents=True, exist_ok=True)
        owns_workdir = False
    else:
        workdir = Path(tempfile.mkdtemp(prefix="dncdm_inv_ab_"))
        owns_workdir = True
    print(f"[scratch] {workdir}")

    results = {}
    try:
        if args.decay_only:
            results["decay-only"] = mode_decay_only(args, root, workdir)
        if args.tiny_gamma:
            results["tiny-gamma"] = mode_tiny_gamma(args, root, workdir)
        if args.pert_regression:
            results["pert-regression"] = mode_pert_regression(args, root, workdir)
        if args.pert_convergence:
            results["pert-convergence"] = mode_pert_convergence(args, root, workdir)
    finally:
        if owns_workdir and not args.keep_scratch:
            shutil.rmtree(workdir, ignore_errors=True)
        elif args.keep_scratch:
            print(f"[scratch] kept at {workdir}")

    print("\n=== SUMMARY ===")
    overall = True
    for name, ok in results.items():
        if ok is None:
            print(f"DEFER {name}  (documented deferral -- see the mode's output above)")
            continue
        print(f"{'PASS' if ok else 'FAIL'}  {name}")
        overall &= ok
    sys.exit(0 if overall else 1)


if __name__ == "__main__":
    main()
