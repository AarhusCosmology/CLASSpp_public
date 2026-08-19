#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

#include "background.h"
#include "composite_species.h"
#include "dark_radiation_species.h"
#include "dncdm_dr_species.h"
#include "dncdm_species.h"
#include "ncdm_base_species.h"
#include "parser.h"
#include "perturb_source_context.h"
#include "species/background_ic_context.h"

// Issue #372: DNCDMSpecies must migrate its Boltzmann hierarchy when the
// perturbation vector is reallocated at an approximation switch (TCA/RSA of
// OTHER species — DNCDM itself has no approximations, so the layout shape is
// always preserved). The new vector arrives zero-initialized; a species that
// does not copy restarts its hierarchy from zero mid-integration.

namespace {

// Fill an NCDM-family layout: q_size bins of (l_max+1) contiguous multipoles,
// starting at base offset `idx0` (different offsets for old/new simulate the
// slot shift when another species' hierarchy appears/disappears).
void FillNcdmLayout(NCDMBaseSpecies::PerturbLayout& l, int q_size, int l_max, int idx0) {
  l.q_size = q_size;
  l.l_max  = l_max;
  l.index_per_q.clear();
  for (int iq = 0; iq < q_size; ++iq)
    l.index_per_q.push_back(idx0 + iq * (l_max + 1));
}

// Distinct nonzero marker per (iq, l) slot.
double Marker(int iq, int l) {
  return 1. + 100. * iq + l;
}

int CheckNcdmSlotsCopied(const NCDMBaseSpecies::PerturbLayout& old_l,
                         const NCDMBaseSpecies::PerturbLayout& new_l,
                         const std::vector<double>& old_y,
                         const std::vector<double>& new_y,
                         const char* what) {
  int rc = 0;
  for (int iq = 0; iq < new_l.q_size; ++iq) {
    for (int l = 0; l <= new_l.l_max; ++l) {
      const double expect = old_y[old_l.index_per_q[iq] + l];
      const double got    = new_y[new_l.index_per_q[iq] + l];
      if (got != expect) {
        std::fprintf(stderr,
                     "FAIL: %s slot (iq=%d, l=%d) not migrated: got %g, expected %g\n",
                     what,
                     iq,
                     l,
                     got,
                     expect);
        rc = 1;
      }
    }
  }
  return rc;
}

std::unique_ptr<DNCDMSpecies> MakeDncdm(const background* pba) {
  FileContent fc;
  fc.set("dncdm1.type", "ncdm_decay_dr");
  fc.set("dncdm1.m", "1.0");
  fc.set("dncdm1.Gamma", "1e3");
  fc.set("dncdm1.Omega_dncdmdr", "0.001");
  NcdmSettings settings{};
  settings.h           = 0.67556;
  settings.T_cmb       = 2.7255;
  settings.tol_ncdm    = 1.e-3;
  settings.tol_ncdm_bg = 1.e-5;
  settings.tol_M_ncdm  = 1.e-5;
  return std::make_unique<DNCDMSpecies>(&fc, "dncdm1", settings, pba, /*bgm=*/nullptr);
}

// Single-grid DNCDM (q_bg == q): a manual quadrature strategy with equal bins.
// DNCDM's ComputeBackground evolves its background PSD on the perturbation grid
// (ctor comment), so both decay-only and collision-owned modes require q_size ==
// q_size_bg; the default qm_auto with differing tol_ncdm/tol_ncdm_bg does not
// guarantee that. The inverse composite enforces this at Create (Task 3.2).
std::unique_ptr<DNCDMSpecies> MakeDncdmSingleGrid(const background* pba) {
  FileContent fc;
  fc.set("dncdm1.type", "ncdm_decay_dr");
  fc.set("dncdm1.m", "1.0");
  fc.set("dncdm1.Gamma", "1e3");
  fc.set("dncdm1.Omega_dncdmdr", "0.001");
  fc.set("dncdm1.quadrature_strategy", "3");  // qm_trapz (manual)
  fc.set("dncdm1.momenta_bins", "20");
  NcdmSettings settings{};
  settings.h           = 0.67556;
  settings.T_cmb       = 2.7255;
  settings.tol_ncdm    = 1.e-3;
  settings.tol_ncdm_bg = 1.e-5;
  settings.tol_M_ncdm  = 1.e-5;
  return std::make_unique<DNCDMSpecies>(&fc, "dncdm1", settings, pba, /*bgm=*/nullptr);
}

}  // namespace

int main() {
  background pba{};
  pba.H0 = 1e-4;

  int rc = 0;
  const PerturbSwitchContext ctx{};  // plain copies ignore the context body

  // ── 1. Standalone DNCDMSpecies: hierarchy must survive a vector switch ────
  {
    auto sp          = MakeDncdm(&pba);
    const int q_size = sp->q_size();
    const int l_max  = 5;
    if (q_size <= 0) {
      std::fprintf(stderr, "FAIL: DNCDMSpecies quadrature is empty (q_size=%d)\n", q_size);
      return 1;
    }

    NCDMBaseSpecies::PerturbLayout old_l, new_l;
    FillNcdmLayout(old_l, q_size, l_max, /*idx0=*/7);
    FillNcdmLayout(new_l, q_size, l_max, /*idx0=*/3);

    std::vector<double> old_y(7 + q_size * (l_max + 1), 0.);
    std::vector<double> new_y(7 + q_size * (l_max + 1), 0.);
    for (int iq = 0; iq < q_size; ++iq)
      for (int l = 0; l <= l_max; ++l)
        old_y[old_l.index_per_q[iq] + l] = Marker(iq, l);

    sp->CopyPerturbationsAcrossSwitch(old_l, new_l, old_y.data(), new_y.data(), ctx);
    rc |= CheckNcdmSlotsCopied(old_l, new_l, old_y, new_y, "DNCDMSpecies");

    // Slots outside the species' layout must be left alone (still zero).
    for (int i = 0; i < 3; ++i) {
      if (new_y[i] != 0.) {
        std::fprintf(stderr, "FAIL: DNCDMSpecies wrote outside its layout (slot %d)\n", i);
        rc = 1;
      }
    }

    // Uninitialized layouts (shape -1, e.g. mode without this species) must be
    // a safe no-op.
    NCDMBaseSpecies::PerturbLayout empty_old, empty_new;
    sp->CopyPerturbationsAcrossSwitch(empty_old, empty_new, old_y.data(), new_y.data(), ctx);
  }

  // ── 2. DNCDM_DR composite: BOTH children must migrate ─────────────────────
  {
    DNCDM_DR_Species dncdm_dr(MakeDncdm(&pba), &pba, /*bgm=*/nullptr);
    const int q_size  = dncdm_dr.dncdm().q_size();
    const int l_max   = 5;
    const int dr_lmax = 3;

    auto old_base = dncdm_dr.CreatePerturbLayout();
    auto new_base = dncdm_dr.CreatePerturbLayout();
    auto& old_l   = static_cast<CompositeSpecies::PerturbLayout&>(*old_base);
    auto& new_l   = static_cast<CompositeSpecies::PerturbLayout&>(*new_base);

    auto& old_nl = static_cast<NCDMBaseSpecies::PerturbLayout&>(*old_l.child_layouts[0]);
    auto& new_nl = static_cast<NCDMBaseSpecies::PerturbLayout&>(*new_l.child_layouts[0]);
    auto& old_dl = static_cast<DarkRadiationSpecies::PerturbLayout&>(*old_l.child_layouts[1]);
    auto& new_dl = static_cast<DarkRadiationSpecies::PerturbLayout&>(*new_l.child_layouts[1]);

    // Old vector: [0..6 other species][dr hierarchy][dncdm hierarchy]
    // New vector: [0..2 other species][dr hierarchy][dncdm hierarchy]
    old_dl.l_max  = dr_lmax;
    old_dl.idx_F0 = 7;
    new_dl.l_max  = dr_lmax;
    new_dl.idx_F0 = 3;
    FillNcdmLayout(old_nl, q_size, l_max, old_dl.idx_F0 + dr_lmax + 1);
    FillNcdmLayout(new_nl, q_size, l_max, new_dl.idx_F0 + dr_lmax + 1);

    const int old_size = old_nl.index_per_q[0] + q_size * (l_max + 1);
    std::vector<double> old_y(old_size, 0.);
    std::vector<double> new_y(old_size, 0.);
    for (int l = 0; l <= dr_lmax; ++l)
      old_y[old_dl.idx_F0 + l] = 1000. + l;
    for (int iq = 0; iq < q_size; ++iq)
      for (int l = 0; l <= l_max; ++l)
        old_y[old_nl.index_per_q[iq] + l] = Marker(iq, l);

    dncdm_dr.CopyPerturbationsAcrossSwitch(old_l, new_l, old_y.data(), new_y.data(), ctx);

    for (int l = 0; l <= dr_lmax; ++l) {
      if (new_y[new_dl.idx_F0 + l] != old_y[old_dl.idx_F0 + l]) {
        std::fprintf(stderr, "FAIL: DNCDM_DR dr child slot l=%d not migrated\n", l);
        rc = 1;
      }
    }
    rc |= CheckNcdmSlotsCopied(old_nl, new_nl, old_y, new_y, "DNCDM_DR dncdm child");
  }

  // ── 3. collision-owned mode: f-variables + kernel-supplied RHS (Task 3.1) ──
  // With SetCollisionOwned(true) the parent drops the ln f + separate-dlnfdlnq
  // integration variables for a single f-per-bin variable; the composite writes
  // the RHS. Background columns and the perturbation surface are unchanged.
  {
    auto sp          = MakeDncdmSingleGrid(&pba);
    const int q_size = sp->q_size();
    // The parent reuses its DNCDM grid as a single grid (design §5); the
    // f-per-bin variable and ComputeMomenta integration must share one grid.
    if (sp->q_size_bg() != q_size) {
      std::fprintf(stderr,
                   "FAIL: DNCDM single-grid invariant broken: q_size=%d q_size_bg=%d\n",
                   q_size,
                   sp->q_size_bg());
      rc = 1;
    }

    sp->SetCollisionOwned(true);
    int index_bg = 0, index_bi = 0;
    sp->RegisterBackgroundIndices(index_bg);
    sp->RegisterIntegrationIndices(index_bi);
    // Collision-owned integration registers ONLY the f-variable (q_size), not the
    // decay-only ln f + separate-dlnfdlnq pair (2*q_size).
    if (index_bi != q_size) {
      std::fprintf(stderr, "FAIL: collision-owned index_bi=%d, expected %d\n", index_bi, q_size);
      rc = 1;
    }

    std::vector<double> bi(index_bi, 0.), bg(index_bg, 0.);
    BackgroundICContext icc{};
    icc.a_ini                = 1e-14;
    icc.pvecback_integration = bi.data();
    sp->SetBackgroundInitialConditions(icc);
    // IC seeds f = w_/dq_ (positive) into the single f-variable.
    if (!(bi[sp->bi_f_parent_index()] > 0.)) {
      std::fprintf(stderr, "FAIL: collision-owned IC did not seed a positive f\n");
      rc = 1;
    }

    // Overwrite with a synthetic positive f-state and check the background readout.
    const auto& q  = sp->GetQ();
    const auto& dq = sp->dq();
    const double a = 0.3;
    std::vector<double> f(q_size);
    for (int i = 0; i < q_size; ++i) {
      f[i] = 0.5 / (1. + q[i] * q[i]);
      // Integration units: the f-slot carries DNCDMSpecies::kFScale, the published
      // ln f column does not (ComputeBackground divides it back out).
      bi[sp->bi_f_parent_index() + i] = DNCDMSpecies::kFScale * f[i];
    }
    sp->ComputeBackground(a, bi.data(), bg.data());

    // Rho matches the closed-form massive quadrature ρ = factor/a⁴ Σ dq q² ε f.
    double rho_exp = 0.;
    for (int i = 0; i < q_size; ++i) {
      const double eps  = std::sqrt(q[i] * q[i] + sp->M() * sp->M() * a * a);
      rho_exp          += q[i] * q[i] * eps * f[i] * dq[i];
    }
    rho_exp *= sp->factor() / std::pow(a, 4);
    if (std::fabs(sp->Rho(bg.data()) - rho_exp) > 1e-10 * rho_exp) {
      std::fprintf(stderr,
                   "FAIL: collision-owned Rho=%g, closed-form=%g\n",
                   sp->Rho(bg.data()),
                   rho_exp);
      rc = 1;
    }
    // The published lnf column is ln f of the synthetic state (perturbations read it).
    for (int i = 0; i < q_size; ++i) {
      const double want = std::log(f[i]);
      if (std::fabs(bg[sp->bg_lnf_index() + i] - want) > 1e-10 * std::fabs(want)) {
        std::fprintf(stderr,
                     "FAIL: collision-owned lnf[%d]=%g, expected %g\n",
                     i,
                     bg[sp->bg_lnf_index() + i],
                     want);
        rc = 1;
      }
    }
    // dlnfdlnq_separate column is written 0 in collision-owned mode.
    if (bg[sp->bg_dlnfdlnq_sep_index()] != 0.) {
      std::fprintf(stderr, "FAIL: collision-owned dlnfdlnq_separate not zeroed\n");
      rc = 1;
    }

    // ── Soft floor arithmetic ───────────────────────────────────────────────
    // log(max(f,0) + eps), not log(max(f, eps)). The soft part keeps the column
    // free of the kink that a hard floor puts into a curve the perturbation module
    // cubic-splines; the max(f,0) keeps the logarithm total when the integrator
    // overshoots to a negative occupation, which f + eps alone does not (it
    // returned NaN at Gamma=1e9).
    //
    // FloorReport is NOT exercised here: it is gated on
    // BackgroundModule::StoringBackgroundTable() so that it sees stored rows and
    // not the evolver's throwaway trial states, and this driver has no module
    // (bgm_ == nullptr). The gate staying shut is itself the thing worth pinning --
    // an ungated version of this guard reported the lasing regime as broken when
    // the accepted solution never came within 75 decades of the floor.
    for (int i = 0; i < q_size; ++i)
      bi[sp->bi_f_parent_index() + i] = DNCDMSpecies::kFScale * 1e-120;  // under the floor
    sp->ComputeBackground(a, bi.data(), bg.data());
    const double pinned = bg[sp->bg_lnf_index()];
    if (!(std::fabs(pinned - std::log(1e-100)) < 1e-6)) {
      std::fprintf(stderr,
                   "FAIL: pinned lnf = %g, expected ln(1e-100) = %g\n",
                   pinned,
                   std::log(1e-100));
      rc = 1;
    }
    // Negative occupation must clamp, not produce NaN.
    for (int i = 0; i < q_size; ++i)
      bi[sp->bi_f_parent_index() + i] = -DNCDMSpecies::kFScale * 1e-40;
    sp->ComputeBackground(a, bi.data(), bg.data());
    const double clamped = bg[sp->bg_lnf_index()];
    if (!std::isfinite(clamped) || std::fabs(clamped - std::log(1e-100)) > 1e-6) {
      std::fprintf(stderr, "FAIL: negative f gave lnf = %g, expected finite ln(1e-100)\n", clamped);
      rc = 1;
    }
    if (sp->BackgroundFloorReport().touched || sp->BackgroundFloorReport().negative) {
      std::fprintf(stderr, "FAIL: FloorReport fired without a BackgroundModule\n");
      rc = 1;
    }
  }

  // ── 4. decay-only mode is unchanged: still registers 2*q_size integration vars ──
  {
    auto sp      = MakeDncdm(&pba);  // collision_owned_ defaults to false
    int index_bi = 0;
    sp->RegisterIntegrationIndices(index_bi);
    if (index_bi != 2 * sp->q_size()) {
      std::fprintf(stderr,
                   "FAIL: decay-only index_bi=%d, expected %d\n",
                   index_bi,
                   2 * sp->q_size());
      rc = 1;
    }
  }

  // ── 5. decay-only mode: a deeply decayed parent must STAY decayed ──────────
  // The state here is ln f = -1000 in every bin, i.e. a parent that has decayed
  // away entirely -- reached by Gamma * t_0 = 1000, so by Gamma ~ 7e4 km/s/Mpc in
  // a standard cosmology, and by a factor 20 more at the Gamma ~ 1e6..1e8 this
  // sector is used at.
  //
  // ComputeBackground used to answer such a state by substituting the PRISTINE
  // Fermi-Dirac f0 in every bin (guard: `ln f <= -460 in any bin`), which put the
  // parent's full undecayed density back into the budget and let the composite's
  // a*Gamma*M*n term go on sourcing decay radiation forever. Under evolver=rkdp45
  // that aborted the run (rho_crit <= 0); under the default ndf15 it did not abort
  // at all and silently returned a universe that was 94% decay radiation, 5.5 Gyr
  // old, with P(k) wrong by a factor 29.
  //
  // Pinned here rather than through a full solve because the substitution is a
  // pure function of one background row: ln f in, ln f out.
  {
    auto sp          = MakeDncdmSingleGrid(&pba);  // collision_owned_ defaults to false
    const int q_size = sp->q_size();
    int index_bi = 0, index_bg = 0;
    sp->RegisterIntegrationIndices(index_bi);
    sp->RegisterBackgroundIndices(index_bg);

    const double kLnf = -1000.;
    std::vector<double> bi(index_bi, 0.), bg(index_bg, 0.);
    for (int i = 0; i < q_size; ++i) {
      bi[sp->bi_lnf_index() + i]          = kLnf;
      bi[sp->bi_dlnfdlnq_sep_index() + i] = -0.5 * (i + 1);  // any distinctive value
    }
    sp->ComputeBackground(/*a=*/0.5, bi.data(), bg.data());

    for (int i = 0; i < q_size; ++i) {
      if (bg[sp->bg_lnf_index() + i] != kLnf) {
        std::fprintf(stderr,
                     "FAIL: decayed parent resurrected: published lnf[%d]=%g, state was %g\n",
                     i,
                     bg[sp->bg_lnf_index() + i],
                     kLnf);
        rc = 1;
        break;
      }
    }
    // exp(-1000) underflows to 0, so the density must be exactly zero -- the
    // statement that there is no parent left to source decay radiation.
    if (sp->Rho(bg.data()) != 0. || bg[sp->bg_number_index()] != 0.) {
      std::fprintf(stderr,
                   "FAIL: decayed parent still carries rho=%g, n=%g (expected 0)\n",
                   sp->Rho(bg.data()),
                   bg[sp->bg_number_index()]);
      rc = 1;
    }
    // The separately integrated dlnf/dlnq column must be published, not left stale
    // (the old fallback branch skipped it).
    for (int i = 0; i < q_size; ++i) {
      if (bg[sp->bg_dlnfdlnq_sep_index() + i] != -0.5 * (i + 1)) {
        std::fprintf(stderr,
                     "FAIL: dlnfdlnq_separate[%d]=%g not published (expected %g)\n",
                     i,
                     bg[sp->bg_dlnfdlnq_sep_index() + i],
                     -0.5 * (i + 1));
        rc = 1;
        break;
      }
    }
  }

  if (rc != 0)
    return 1;

  std::printf("dncdm switch-copy test passed\n");
  return 0;
}
