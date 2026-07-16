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

  if (rc != 0)
    return 1;

  std::printf("dncdm switch-copy test passed\n");
  return 0;
}
