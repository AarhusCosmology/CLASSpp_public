// Unit test for DNCDMInvSpecies (ν_H ↔ ν_l + φ composite, Phase 3): input-time
// guards, construction + child order, background index plumbing, the kernel
// conservation identities through the composite plumbing, the combined-sector
// shooting closure, and the #372 switch-copy regression across three children.
#include <cassert>
#include <cmath>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "background.h"
#include "decay_transition_kernel.h"
#include "parser.h"
#include "species/background_ic_context.h"
#include "species/dncdm_inv_species.h"
#include "species/dncdm_species.h"
#include "species/dr_psd_species.h"
#include "species/ncdm_base_species.h"
#include "species/perturb_source_context.h"
#include "species/shooting_target.h"
#include "species/species_build_context.h"

static NcdmSettings TestSettings() {
  NcdmSettings s{};
  s.h           = 0.67;
  s.T_cmb       = 2.7255;
  s.tol_ncdm    = 1e-3;
  s.tol_ncdm_bg = 1e-5;
  s.tol_M_ncdm  = 1e-8;
  return s;
}

// Common inverse-decay instance: single-grid parent (manual quadrature, equal
// bins) + a daughter grid wide enough to cover the decay bands (no clamping).
static void SetInvBase(FileContent& fc) {
  fc.set("dncdm1.type", "ncdm_decay_dr");
  fc.set("dncdm1.m", "1.0");
  fc.set("dncdm1.Gamma", "10");               // Gamma/H0 < 1: finite ctor deg guess
  fc.set("dncdm1.Omega_dncdmdr", "0.001");    // combined mode (daughters start empty)
  fc.set("dncdm1.quadrature_strategy", "3");  // qm_trapz (single grid: q == q_bg)
  fc.set("dncdm1.momenta_bins", "16");
  fc.set("dncdm1.dr_q_min", "1e-4");
  fc.set("dncdm1.dr_q_max", "300");
  fc.set("dncdm1.dr_N_q", "100");
}

static background MakeBackground() {
  background pba{};
  pba.H0       = 2.2e-4;
  pba.Omega0_g = 5.4e-5;  // finite radiation so the ctor's k_rad guess is regular
  return pba;
}

static std::unique_ptr<DNCDMInvSpecies> BuildComposite(FileContent& fc,
                                                       background& pba,
                                                       NcdmSettings& settings) {
  auto parent = std::make_unique<DNCDMSpecies>(&fc, "dncdm1", settings, &pba, nullptr);
  SpeciesBuildContext ctx{};
  ctx.pfc           = &fc;
  ctx.pba           = &pba;
  ctx.ncdm_settings = &settings;
  ctx.bgm           = nullptr;
  Named n           = DNCDMInvSpecies::Create(std::move(parent), ctx);
  auto* raw         = dynamic_cast<DNCDMInvSpecies*>(n.species.get());
  assert(raw != nullptr);
  n.species.release();
  return std::unique_ptr<DNCDMInvSpecies>(raw);
}

// Build the composite DIRECTLY with an explicit kernel Config (bypasses the Create
// factory + guards). Only way to get a decay-only (inverse=off) composite for the
// §4.3 cancellation test, since Create routes inverse=off to DNCDM_DR instead.
static std::unique_ptr<DNCDMInvSpecies> BuildCompositeDirect(
    FileContent& fc, background& pba, NcdmSettings& settings, bool inv, bool qs) {
  auto parent = std::make_unique<DNCDMSpecies>(&fc, "dncdm1", settings, &pba, nullptr);
  auto fermion =
      std::make_unique<DrPsdSpecies>(&fc, "dncdm1", settings, &pba, nullptr, Statistics::Fermion);
  auto boson =
      std::make_unique<DrPsdSpecies>(&fc, "dncdm1", settings, &pba, nullptr, Statistics::Boson);
  DecayTransitionKernel::Config cfg;
  cfg.inverse_decays     = inv;
  cfg.quantum_statistics = qs;
  return std::make_unique<DNCDMInvSpecies>(std::move(parent),
                                           std::move(fermion),
                                           std::move(boson),
                                           cfg,
                                           &pba,
                                           nullptr);
}

static double Scale(double a, double b) {
  const double m = std::max(std::fabs(a), std::fabs(b));
  return m > 1e-300 ? m : 1e-300;
}
static double Scale3(double a, double b, double c) {
  const double m = std::max(std::fabs(a), std::max(std::fabs(b), std::fabs(c)));
  return m > 1e-300 ? m : 1e-300;
}

static void FillNcdmLayout(NCDMBaseSpecies::PerturbLayout& l, int q_size, int l_max, int idx0) {
  l.q_size = q_size;
  l.l_max  = l_max;
  l.index_per_q.clear();
  for (int iq = 0; iq < q_size; ++iq)
    l.index_per_q.push_back(idx0 + iq * (l_max + 1));
}
static double Marker(int child, int iq, int l) {
  return 1. + 1e4 * child + 100. * iq + l;
}

int main() {
  background pba        = MakeBackground();
  NcdmSettings settings = TestSettings();

  // ── Guards ────────────────────────────────────────────────────────────────
  {
    // quantum_statistics without inverse_decays -> class_test_severe
    // (invalid_argument): a structural flag-combination error, decidable from
    // which keys are set, not from any sampler-varyable numeric value.
    FileContent fc;
    SetInvBase(fc);
    fc.set("dncdm1.inverse_decays", "no");
    fc.set("dncdm1.quantum_statistics", "yes");
    bool threw = false;
    try {
      BuildComposite(fc, pba, settings);
    }
    catch (const std::invalid_argument&) {
      threw = true;
    }
    assert(threw);
  }
  {
    // newtonian gauge with inverse decays -> class_test_severe (invalid_argument).
    FileContent fc;
    SetInvBase(fc);
    fc.set("dncdm1.inverse_decays", "yes");
    fc.set("gauge", "newtonian");
    bool threw = false;
    try {
      BuildComposite(fc, pba, settings);
    }
    catch (const std::invalid_argument&) {
      threw = true;
    }
    assert(threw);
  }

  // ── Construction + child order + collision-owned parent ─────────────────────
  {
    FileContent fc;
    SetInvBase(fc);
    fc.set("dncdm1.inverse_decays", "yes");
    fc.set("dncdm1.quantum_statistics", "no");
    auto comp = BuildComposite(fc, pba, settings);

    assert(comp->fermion().statistics() == Statistics::Fermion);
    assert(comp->boson().statistics() == Statistics::Boson);
    // Parent is collision-owned -> single f-per-bin integration variable.
    int index_bi = 0;
    comp->RegisterIntegrationIndices(index_bi);
    const int expect = comp->parent().q_size() + static_cast<int>(comp->fermion().q_bg().size()) +
                       static_cast<int>(comp->boson().q_bg().size());
    assert(index_bi == expect);
    // Child f-slots are contiguous and correctly ordered (parent, fermion, boson).
    assert(comp->parent().bi_f_parent_index() == 0);
    assert(comp->fermion().bi_f_index() == comp->parent().q_size());
  }

  // ── Background plumbing: Rho sums children on a synthetic state ──────────────
  {
    FileContent fc;
    SetInvBase(fc);
    fc.set("dncdm1.inverse_decays", "yes");
    auto comp = BuildComposite(fc, pba, settings);

    int index_bg = 0, index_bi = 0;
    comp->RegisterBackgroundIndices(index_bg);
    comp->RegisterIntegrationIndices(index_bi);
    std::vector<double> bi(index_bi, 0.), bg(index_bg, 0.);

    BackgroundICContext icc{};
    icc.a_ini                = 1e-14;
    icc.pvecback_integration = bi.data();
    comp->SetBackgroundInitialConditions(icc);
    // Seed a synthetic positive f-state on all three species.
    auto& parent = comp->parent();
    auto& ferm   = comp->fermion();
    auto& bos    = comp->boson();
    for (int i = 0; i < parent.q_size(); ++i)
      // The parent's f-slot is an INTEGRATION variable and carries kFScale (see
      // DNCDMSpecies::kFScale); the published ln f column does not.
      bi[parent.bi_f_parent_index() + i] = DNCDMSpecies::kFScale * 0.5 /
                                           (1. + parent.GetQ()[i] * parent.GetQ()[i]);
    for (int i = 0; i < (int) ferm.q_bg().size(); ++i)
      bi[ferm.bi_f_index() + i] = 0.3 / (1. + ferm.q_bg()[i]);
    for (int i = 0; i < (int) bos.q_bg().size(); ++i)
      bi[bos.bi_f_index() + i] = 0.2 / (1. + bos.q_bg()[i]);

    const double a = 1e-3;
    comp->ComputeBackground(a, bi.data(), bg.data());
    const double rho_sum = parent.Rho(bg.data()) + ferm.Rho(bg.data()) + bos.Rho(bg.data());
    assert(rho_sum > 0.);
    assert(std::fabs(comp->Rho(bg.data()) - rho_sum) <= 1e-12 * rho_sum);
  }

  // ── Conservation identities through the composite plumbing (1e-11) ──────────
  // Assemble a random f-state, run the kernel via BackgroundDerivs' worker, read
  // the bare kernel moments back. The Fable-review identities (×2 boson dof folded
  // into df_phi) are exact regardless of the f-state: this proves the composite
  // gathers f and scatters df into the correct slots, not just the kernel algebra.
  {
    FileContent fc;
    SetInvBase(fc);
    fc.set("dncdm1.inverse_decays", "yes");
    fc.set("dncdm1.quantum_statistics", "yes");
    auto comp = BuildComposite(fc, pba, settings);

    int index_bg = 0, index_bi = 0;
    comp->RegisterBackgroundIndices(index_bg);
    comp->RegisterIntegrationIndices(index_bi);
    std::vector<double> y(index_bi, 0.), dy(index_bi, 0.);

    auto& parent = comp->parent();
    auto& ferm   = comp->fermion();
    auto& bos    = comp->boson();
    unsigned s   = 0x1234567u;
    auto rnd     = [&]() {
      s = s * 1664525u + 1013904223u;
      return 0.1 + 0.9 * (((s >> 8) & 0xffffu) / 65535.0);
    };
    for (int i = 0; i < parent.q_size(); ++i)
      y[parent.bi_f_parent_index() + i] = DNCDMSpecies::kFScale * rnd();
    for (int i = 0; i < (int) ferm.q_bg().size(); ++i)
      y[ferm.bi_f_index() + i] = rnd();
    for (int i = 0; i < (int) bos.q_bg().size(); ++i)
      y[bos.bi_f_index() + i] = rnd();

    const double a = 1e-3;
    comp->ApplyKernelBackgroundDerivs(a, y.data(), dy.data());
    const double clamp = comp->ClampedEnergyResidual();
    const auto M       = comp->ConservationMoments(a, dy.data());

    const double id_Nf = std::fabs(M.N_H + M.N_l);
    const double id_Nb = std::fabs(2. * M.N_H + M.N_phi);
    // The energy identity is exact up to the energy the LUMPED daughter loss
    // misplaces to buy positivity, and that is booked exactly -- so this still pins
    // a machine-precision statement (nothing goes missing unaccounted), not a
    // loosened tolerance. The f here is uniform random per bin, i.e. the roughest
    // input the split can be handed; on the smooth PSDs of a real run the residual
    // is orders smaller (see test_daughter_positivity for its convergence).
    const double split = comp->SplitEnergyResidual();
    const double id_E  = std::fabs(2. * M.E_H + 2. * M.E_l + M.E_phi - split);
    std::printf(
        "  conservation: clamp=%.3e  N_H+N_l=%.3e  2N_H+N_phi=%.3e  |2E_H+2E_l+E_phi - booked|=%.3e"
        "  (booked %.3e)  N_H=%.3e\n",
        clamp,
        id_Nf,
        id_Nb,
        id_E,
        split,
        M.N_H);

    assert(clamp == 0.);  // daughter grids cover the bands (no off-grid deposits)
    assert(id_Nf < 1e-11 * Scale(M.N_H, M.N_l));
    assert(id_Nb < 1e-11 * Scale(2. * M.N_H, M.N_phi));
    assert(id_E < 1e-11 * Scale3(2. * M.E_H, 2. * M.E_l, M.E_phi));
    assert(split != 0.);            // else the assertion above is vacuous
    assert(std::fabs(M.N_H) > 0.);  // guards against a trivially-zero "pass"
  }

  // ── Combined-sector shooting closure (Task 3.3) ─────────────────────────────
  // Combined mode (Omega_dncdmdr): one target driving deg, GetOmega0 pinned to the
  // reserve, and the residual summing all THREE children's today density.
  {
    FileContent fc;
    SetInvBase(fc);
    fc.set("dncdm1.inverse_decays", "yes");
    auto comp = BuildComposite(fc, pba, settings);

    const auto targets = comp->GetShootingTargets();
    assert(targets.size() == 1);
    assert(targets[0].target_name == "dncdm1.Omega_dncdmdr");
    assert(targets[0].unknown_param == "dncdm1.deg");
    assert(comp->GetOmega0() == 0.001);  // pinned combined reserve

    std::vector<double> guess, dxdy;
    SpeciesBuildContext gctx{};
    gctx.pfc           = &fc;
    gctx.pba           = &pba;
    gctx.ncdm_settings = &settings;
    comp->ComputeShootingGuess(gctx, guess, dxdy);
    assert(guess.size() == 1 && guess[0] > 0.);  // finite deg guess

    int index_bg = 0, index_bi = 0;
    comp->RegisterBackgroundIndices(index_bg);
    comp->RegisterIntegrationIndices(index_bi);
    std::vector<double> bi(index_bi, 0.), bg(index_bg, 0.);
    BackgroundICContext icc{};
    icc.pvecback_integration = bi.data();
    comp->SetBackgroundInitialConditions(icc);
    auto& parent = comp->parent();
    auto& ferm   = comp->fermion();
    auto& bos    = comp->boson();
    for (int i = 0; i < parent.q_size(); ++i)
      // The parent's f-slot is an INTEGRATION variable and carries kFScale (see
      // DNCDMSpecies::kFScale); the published ln f column does not.
      bi[parent.bi_f_parent_index() + i] = DNCDMSpecies::kFScale * 0.5 /
                                           (1. + parent.GetQ()[i] * parent.GetQ()[i]);
    for (int i = 0; i < (int) ferm.q_bg().size(); ++i)
      bi[ferm.bi_f_index() + i] = 0.3 / (1. + ferm.q_bg()[i]);
    for (int i = 0; i < (int) bos.q_bg().size(); ++i)
      bi[bos.bi_f_index() + i] = 0.2 / (1. + bos.q_bg()[i]);
    comp->ComputeBackground(1.0, bi.data(), bg.data());  // a = 1 (today)

    const double combined = (parent.Rho(bg.data()) + ferm.Rho(bg.data()) + bos.Rho(bg.data())) /
                            (pba.H0 * pba.H0);
    ShootingResidualContext rctx{};
    rctx.pba       = &pba;
    rctx.bg_today  = bg.data();
    const double r = comp->ComputeShootingResidual(rctx, targets[0]);
    assert(std::fabs(r - (combined - 0.001)) <= 1e-10 * std::max(std::fabs(combined), 1e-30));
  }

  // ── Switch-copy across all three children (#372 class) ──────────────────────
  {
    FileContent fc;
    SetInvBase(fc);
    fc.set("dncdm1.inverse_decays", "yes");
    auto comp = BuildComposite(fc, pba, settings);

    const int qp    = comp->parent().q_size();
    const int qf    = static_cast<int>(comp->fermion().q_bg().size());
    const int qb    = static_cast<int>(comp->boson().q_bg().size());
    const int l_max = 4;

    auto old_base = comp->CreatePerturbLayout();
    auto new_base = comp->CreatePerturbLayout();
    auto& old_l   = static_cast<CompositeSpecies::PerturbLayout&>(*old_base);
    auto& new_l   = static_cast<CompositeSpecies::PerturbLayout&>(*new_base);
    assert(old_l.child_layouts.size() == 3);

    auto& oP = static_cast<NCDMBaseSpecies::PerturbLayout&>(*old_l.child_layouts[0]);
    auto& oF = static_cast<NCDMBaseSpecies::PerturbLayout&>(*old_l.child_layouts[1]);
    auto& oB = static_cast<NCDMBaseSpecies::PerturbLayout&>(*old_l.child_layouts[2]);
    auto& nP = static_cast<NCDMBaseSpecies::PerturbLayout&>(*new_l.child_layouts[0]);
    auto& nF = static_cast<NCDMBaseSpecies::PerturbLayout&>(*new_l.child_layouts[1]);
    auto& nB = static_cast<NCDMBaseSpecies::PerturbLayout&>(*new_l.child_layouts[2]);

    // Old vector has a 7-slot lead-in; new vector a 3-slot lead-in (slot shift).
    int o = 7;
    FillNcdmLayout(oP, qp, l_max, o);
    o += qp * (l_max + 1);
    FillNcdmLayout(oF, qf, l_max, o);
    o += qf * (l_max + 1);
    FillNcdmLayout(oB, qb, l_max, o);
    const int total = o + qb * (l_max + 1);

    int noff = 3;
    FillNcdmLayout(nP, qp, l_max, noff);
    noff += qp * (l_max + 1);
    FillNcdmLayout(nF, qf, l_max, noff);
    noff += qf * (l_max + 1);
    FillNcdmLayout(nB, qb, l_max, noff);

    std::vector<double> old_y(total, 0.), new_y(total, 0.);
    const NCDMBaseSpecies::PerturbLayout* olay[3] = {&oP, &oF, &oB};
    const int qn[3]                               = {qp, qf, qb};
    for (int c = 0; c < 3; ++c)
      for (int iq = 0; iq < qn[c]; ++iq)
        for (int l = 0; l <= l_max; ++l)
          old_y[olay[c]->index_per_q[iq] + l] = Marker(c, iq, l);

    const PerturbSwitchContext sctx{};
    comp->CopyPerturbationsAcrossSwitch(old_l, new_l, old_y.data(), new_y.data(), sctx);

    const NCDMBaseSpecies::PerturbLayout* nlay[3] = {&nP, &nF, &nB};
    for (int c = 0; c < 3; ++c)
      for (int iq = 0; iq < qn[c]; ++iq)
        for (int l = 0; l <= l_max; ++l)
          assert(new_y[nlay[c]->index_per_q[iq] + l] == Marker(c, iq, l));
    for (int i = 0; i < 3; ++i)
      assert(new_y[i] == 0.);  // lead-in untouched
  }

  // ── Perturbation coupling: decay-only cancellation (design §4.3) ─────────────
  // Run at dr_bg_refine = 1 (PT grid == BG grid) AND 9 (PT grid a 12-point exact
  // subsample of the 100-point BG grid): the composite must drive the collision on
  // the PT-resolution kernel while the background stays on the fine grid, so the
  // daughters' gather/scatter must be PT-sized and their background columns read at
  // bg index j*refine. An off-by-one there shows up here as a broken cancellation.
  for (int refine : {1, 9}) {
    // MIXED REPRESENTATION (load-bearing): the PARENT evolves normalized Ψ_H, the two
    // daughters evolve UNNORMALIZED F = f̄Ψ = δf (a daughter fills from zero abundance,
    // where f̄̇/f̄ is smooth but ~1e9 at early times, so 1/f̄ must not appear in its
    // equation at all). So AddCouplingDerivs gives the parent (1/f̄_H)(dF_H/dτ)_C,ℓ PLUS
    // −(f̄̇_H/f̄_H)Ψ_H, and the daughters the bare (dF/dτ)_C,ℓ with NO dilution term (it
    // cancels structurally in F-space). In decay-only mode the parent's operator
    // (−K/ε1 Ψ_H) and background term (+K/ε1 Ψ_H) cancel EXACTLY, so the parent adds
    // nothing to dy (reproducing master); the daughters still receive the parent's
    // perturbation through the collision deposit. This drives the whole regression anchor.
    {
      FileContent fc;
      SetInvBase(fc);
      fc.set("dncdm1.dr_bg_refine", std::to_string(refine));
      auto comp    = BuildCompositeDirect(fc, pba, settings, /*inv=*/false, /*qs=*/false);
      auto& parent = comp->parent();
      auto& ferm   = comp->fermion();
      auto& bos    = comp->boson();
      assert(ferm.refine() == refine);
      assert(ferm.q_size() == (static_cast<int>(ferm.q_bg().size()) - 1) / refine + 1);

      int index_bg = 0;
      comp->RegisterBackgroundIndices(index_bg);
      std::vector<double> bg(index_bg, 0.);
      const int qp = parent.q_size();
      // Layouts are PT-sized; the background columns stay BG-sized (filled in full).
      const int qf = ferm.q_size();
      const int qb = bos.q_size();
      for (int i = 0; i < qp; ++i)
        bg[parent.bg_lnf_index() + i] = std::log(0.5 / (1. + parent.GetQ()[i] * parent.GetQ()[i]));
      for (int j = 0; j < (int) ferm.q_bg().size(); ++j)
        bg[ferm.bg_f_index() + j] = 0.3 / (1. + ferm.q_bg()[j]);
      for (int k = 0; k < (int) bos.q_bg().size(); ++k)
        bg[bos.bg_f_index() + k] = 0.2 / (1. + bos.q_bg()[k]);

      const int l_max = 4;
      auto base       = comp->CreatePerturbLayout();
      auto& cl        = static_cast<CompositeSpecies::PerturbLayout&>(*base);
      auto& pL        = static_cast<NCDMBaseSpecies::PerturbLayout&>(*cl.child_layouts[0]);
      auto& fL        = static_cast<NCDMBaseSpecies::PerturbLayout&>(*cl.child_layouts[1]);
      auto& bL        = static_cast<NCDMBaseSpecies::PerturbLayout&>(*cl.child_layouts[2]);
      int off         = 0;
      FillNcdmLayout(pL, qp, l_max, off);
      off += qp * (l_max + 1);
      FillNcdmLayout(fL, qf, l_max, off);
      off += qf * (l_max + 1);
      FillNcdmLayout(bL, qb, l_max, off);
      off += qb * (l_max + 1);

      std::vector<double> y(off, 0.), dy(off, 0.);
      unsigned s = 0xBEEF01u;
      auto rnd   = [&]() {
        s = s * 1664525u + 1013904223u;
        return -1.0 + 2.0 * (((s >> 9) & 0xffffu) / 65535.0);
      };
      for (int iq = 0; iq < qp; ++iq)
        for (int l = 0; l <= l_max; ++l)
          y[pL.index_per_q[iq] + l] = rnd();
      for (int iq = 0; iq < qf; ++iq)
        for (int l = 0; l <= l_max; ++l)
          y[fL.index_per_q[iq] + l] = rnd();
      for (int iq = 0; iq < qb; ++iq)
        for (int l = 0; l <= l_max; ++l)
          y[bL.index_per_q[iq] + l] = rnd();

      const double a = 1e-3;
      // Per-k scratch: the composite's perturbation kernel and gather buffers now
      // live on the workspace, not on the species (thread safety), so the driver
      // supplies one exactly as perturb_workspace_init does.
      auto scratch_owner = comp->CreatePerturbScratch();
      auto& scratch      = static_cast<DNCDMInvSpecies::Scratch&>(*scratch_owner);
      comp->ApplyKernelPerturbDerivs(a, l_max, bg.data(), pL, fL, bL, y.data(), dy.data(), scratch);
      assert(comp->ClampedEnergyResidualPert(scratch) == 0.);

      // Decay-only, stated so it does not depend on the parent's representation
      // (#386): pure decay must leave the parent's NORMALIZED perturbation
      // Ψ_H = F_H/f̄_H untouched, because decay removes δf and f̄ at the same per-bin
      // rate. In Ψ-space that read "dy ≈ 0" — a cancellation of two large terms; in
      // F-space the same statement is that the collision rate matches the background
      // rate bin by bin, dF_H/F_H == f̄̇_H/f̄_H, which is what is checked here. Both
      // sides are ratios, so the kappa unit boundary (#385) drops out.
      const std::vector<double>& dfH = comp->KernelDfBgH(scratch);
      const double kappa             = comp->KappaStoredToBare();
      double max_rate_mismatch = 0., max_daughter = 0.;
      for (int iq = 0; iq < qp; ++iq) {
        const double fbar    = kappa * std::exp(bg[parent.bg_lnf_index() + iq]);
        const double bg_rate = dfH[iq] / fbar;
        for (int l = 0; l <= l_max; ++l) {
          const double psi     = y[pL.index_per_q[iq] + l];
          const double pt_rate = dy[pL.index_per_q[iq] + l] / psi;
          max_rate_mismatch    = std::max(max_rate_mismatch,
                                          std::fabs(pt_rate - bg_rate) / Scale(pt_rate, bg_rate));
        }
      }
      for (int iq = 0; iq < qf; ++iq)
        for (int l = 0; l <= l_max; ++l)
          max_daughter = std::max(max_daughter, std::fabs(dy[fL.index_per_q[iq] + l]));
      for (int iq = 0; iq < qb; ++iq)
        for (int l = 0; l <= l_max; ++l)
          max_daughter = std::max(max_daughter, std::fabs(dy[bL.index_per_q[iq] + l]));
      std::printf(
          "  decay-only, parent Psi invariant (refine=%d): max|dF_H/F_H - fdot/f|=%.3e "
          "max|daughter dy|=%.3e\n",
          refine,
          max_rate_mismatch,
          max_daughter);
      assert(max_daughter > 0.);          // daughters DO get the collision deposit
      assert(max_rate_mismatch < 1e-12);  // parent decay leaves Psi_H exactly invariant
    }

    // ── Perturbation coupling: F-space collision conserves through the plumbing ──
    // Full mode (inv+qs). The kernel's F-space rate obeys the ℓ=0 number/energy and
    // ℓ=1 momentum identities. Reconstructing it from the composite's total coupling is
    // representation-dependent (see the mixed-representation note above):
    //   parent   (Ψ-space): F-rate = f̄_H·dy + f̄̇_H·Ψ_H  (the −(f̄̇/f̄)Ψ term added back)
    //   daughters (F-space): F-rate = dy      (nothing to undo — dy IS the collision)
    // The daughter moment weights are the PT-grid measures q()/dq() -- the grid the
    // perturbation kernel actually deposits on. At refine>1 these differ from
    // q_bg()/dq_bg(), and the identities still hold to machine precision because the
    // coarse grid is an exact subsample (same endpoints, same trapezoid convention).
    {
      FileContent fc;
      SetInvBase(fc);
      fc.set("dncdm1.dr_bg_refine", std::to_string(refine));
      fc.set("dncdm1.quantum_statistics", "yes");
      auto comp    = BuildCompositeDirect(fc, pba, settings, /*inv=*/true, /*qs=*/true);
      auto& parent = comp->parent();
      auto& ferm   = comp->fermion();
      auto& bos    = comp->boson();

      int index_bg = 0;
      comp->RegisterBackgroundIndices(index_bg);
      std::vector<double> bg(index_bg, 0.);
      const int qp = parent.q_size();
      const int qf = ferm.q_size();
      const int qb = bos.q_size();
      for (int i = 0; i < qp; ++i)
        bg[parent.bg_lnf_index() + i] = std::log(0.5 / (1. + parent.GetQ()[i] * parent.GetQ()[i]));
      for (int j = 0; j < (int) ferm.q_bg().size(); ++j)
        bg[ferm.bg_f_index() + j] = 0.3 / (1. + ferm.q_bg()[j]);
      for (int k = 0; k < (int) bos.q_bg().size(); ++k)
        bg[bos.bg_f_index() + k] = 0.2 / (1. + bos.q_bg()[k]);

      const int l_max = 4;
      auto base       = comp->CreatePerturbLayout();
      auto& cl        = static_cast<CompositeSpecies::PerturbLayout&>(*base);
      auto& pL        = static_cast<NCDMBaseSpecies::PerturbLayout&>(*cl.child_layouts[0]);
      auto& fL        = static_cast<NCDMBaseSpecies::PerturbLayout&>(*cl.child_layouts[1]);
      auto& bL        = static_cast<NCDMBaseSpecies::PerturbLayout&>(*cl.child_layouts[2]);
      int off         = 0;
      FillNcdmLayout(pL, qp, l_max, off);
      off += qp * (l_max + 1);
      FillNcdmLayout(fL, qf, l_max, off);
      off += qf * (l_max + 1);
      FillNcdmLayout(bL, qb, l_max, off);
      off += qb * (l_max + 1);

      std::vector<double> y(off, 0.), dy(off, 0.);
      unsigned s = 0x51A7Eu;
      auto rnd   = [&]() {
        s = s * 1664525u + 1013904223u;
        return -1.0 + 2.0 * (((s >> 9) & 0xffffu) / 65535.0);
      };
      for (int iq = 0; iq < qp; ++iq)
        for (int l = 0; l <= l_max; ++l)
          y[pL.index_per_q[iq] + l] = rnd();
      for (int iq = 0; iq < qf; ++iq)
        for (int l = 0; l <= l_max; ++l)
          y[fL.index_per_q[iq] + l] = rnd();
      for (int iq = 0; iq < qb; ++iq)
        for (int l = 0; l <= l_max; ++l)
          y[bL.index_per_q[iq] + l] = rnd();

      const double a = 1e-3;
      // Per-k scratch: the composite's perturbation kernel and gather buffers now
      // live on the workspace, not on the species (thread safety), so the driver
      // supplies one exactly as perturb_workspace_init does.
      auto scratch_owner = comp->CreatePerturbScratch();
      auto& scratch      = static_cast<DNCDMInvSpecies::Scratch&>(*scratch_owner);
      comp->ApplyKernelPerturbDerivs(a, l_max, bg.data(), pL, fL, bL, y.data(), dy.data(), scratch);
      assert(comp->ClampedEnergyResidualPert(scratch) == 0.);
      const double A = a * a * parent.GetMass() * parent.GetMass();

      // The parent's dy IS its F-space rate now (#386): all three species carry
      // δf, so nothing has to be un-normalized here. What does remain is the #385
      // unit boundary — the parent stores kappa-suppressed occupation while the
      // kernel's identities are stated on bare occupation — so its leg, and only
      // its leg, is scaled by kappa, exactly as ApplyKernelPerturbDerivs does.
      const double kappa = comp->KappaStoredToBare();
      double numH = 0, numl = 0, numphi = 0, eH = 0, el = 0, ephi = 0;
      double momH = 0, moml = 0, momphi = 0;
      for (int i = 0; i < qp; ++i) {
        const double q = parent.GetQ()[i], w = parent.dq()[i] * q * q, eps = std::sqrt(q * q + A);
        const double Fr0  = kappa * dy[pL.index_per_q[i] + 0];
        const double Fr1  = kappa * dy[pL.index_per_q[i] + 1];
        numH             += w * Fr0;
        eH               += w * eps * Fr0;
        momH             += parent.dq()[i] * q * q * q * Fr1;
      }
      for (int j = 0; j < qf; ++j) {
        const double q = ferm.q()[j], w = ferm.dq()[j] * q * q;  // PT grid: where the op deposits
        const double Fr0  = dy[fL.index_per_q[j] + 0];
        const double Fr1  = dy[fL.index_per_q[j] + 1];
        numl             += w * Fr0;
        el               += w * q * Fr0;
        moml             += ferm.dq()[j] * q * q * q * Fr1;
      }
      for (int k = 0; k < qb; ++k) {
        const double q = bos.q()[k], w = bos.dq()[k] * q * q;
        const double Fr0  = dy[bL.index_per_q[k] + 0];
        const double Fr1  = dy[bL.index_per_q[k] + 1];
        numphi           += w * Fr0;
        ephi             += w * q * Fr0;
        momphi           += bos.dq()[k] * q * q * q * Fr1;
      }
      const double nf = std::fabs(numH + numl) / Scale(numH, numl);
      const double nb = std::fabs(2. * numH + numphi) / Scale(2. * numH, numphi);
      // The q-weighted identities carry the lumped split's O(dq) misplacement, booked
      // per multipole by the kernel; the NUMBER ones above are still exactly zero. Both
      // legs are already in the kernel's bare units here (the parent leg was scaled by
      // kappa above), which is the same convention the residuals are booked in.
      const double bkE = comp->SplitEnergyResidualPert(scratch);
      const double bkM = comp->SplitMomentumResidualPert(scratch);
      const double en  = std::fabs(2. * eH + 2. * el + ephi - bkE) / Scale3(2. * eH, 2. * el, ephi);
      const double mm  = std::fabs(2. * momH + 2. * moml + momphi - bkM) /
                         Scale3(2. * momH, 2. * moml, momphi);
      std::printf(
          "  pert conservation (refine=%d, N_pt=%d/%d): numH+numl=%.2e 2numH+numphi=%.2e "
          "|2eH+2el+ephi - booked|=%.2e |2momH+2moml+momphi - booked|=%.2e\n",
          refine,
          qf,
          (int) ferm.q_bg().size(),
          nf,
          nb,
          en,
          mm);
      assert(nf < 1e-10);
      assert(nb < 1e-10);
      assert(en < 1e-10);
      assert(mm < 1e-9);
      assert(std::fabs(numH) > 0.);
    }

    // ── Perturbation Jacobian diagonal (evolver_etd) ─────────────────────────────
    // ApplyKernelPerturbDiagonal is the object evolver_etd integrates exactly. That
    // the KERNEL's diagonal is the right quantity is pinned elsewhere — in
    // decay_kernel_test, against finite differences of ComputeBackgroundDerivs and
    // against ApplyPerturbationOperator's own per-ℓ diagonal. What is new here, and
    // what this checks, is the SCATTER: one kernel pass, written into every populated
    // ℓ slot of its bin, accumulated rather than assigned, and stopping at l_max_col.
    //
    // Run with l_max_col = 2 against layouts at l_max = 4, so the truncation is
    // exercised rather than assumed. The "accumulate, do not zero" contract is checked
    // by calling the worker TWICE on one buffer and requiring exactly 2x, not by
    // pre-filling a sentinel: the diagonal here is O(1e-9) and an O(1) sentinel would
    // lose seven digits to cancellation on the way back out, which measures the
    // subtraction rather than the code. x + x is exact; that route stays at zero ULP.
    {
      FileContent fc;
      SetInvBase(fc);
      fc.set("dncdm1.dr_bg_refine", std::to_string(refine));
      fc.set("dncdm1.quantum_statistics", "yes");
      auto comp    = BuildCompositeDirect(fc, pba, settings, /*inv=*/true, /*qs=*/true);
      auto& parent = comp->parent();
      auto& ferm   = comp->fermion();
      auto& bos    = comp->boson();

      int index_bg = 0;
      comp->RegisterBackgroundIndices(index_bg);
      std::vector<double> bg(index_bg, 0.);
      const int qp = parent.q_size();
      const int qf = ferm.q_size();
      const int qb = bos.q_size();
      for (int i = 0; i < qp; ++i)
        bg[parent.bg_lnf_index() + i] = std::log(0.5 / (1. + parent.GetQ()[i] * parent.GetQ()[i]));
      for (int j = 0; j < (int) ferm.q_bg().size(); ++j)
        bg[ferm.bg_f_index() + j] = 0.3 / (1. + ferm.q_bg()[j]);
      for (int k = 0; k < (int) bos.q_bg().size(); ++k)
        bg[bos.bg_f_index() + k] = 0.2 / (1. + bos.q_bg()[k]);

      const int l_max     = 4;
      const int l_max_col = 2;
      auto base           = comp->CreatePerturbLayout();
      auto& cl            = static_cast<CompositeSpecies::PerturbLayout&>(*base);
      auto& pL            = static_cast<NCDMBaseSpecies::PerturbLayout&>(*cl.child_layouts[0]);
      auto& fL            = static_cast<NCDMBaseSpecies::PerturbLayout&>(*cl.child_layouts[1]);
      auto& bL            = static_cast<NCDMBaseSpecies::PerturbLayout&>(*cl.child_layouts[2]);
      int off             = 0;
      FillNcdmLayout(pL, qp, l_max, off);
      off += qp * (l_max + 1);
      FillNcdmLayout(fL, qf, l_max, off);
      off += qf * (l_max + 1);
      FillNcdmLayout(bL, qb, l_max, off);
      off += qb * (l_max + 1);

      std::vector<double> diag(off, 0.);

      const double a     = 1e-3;
      auto scratch_owner = comp->CreatePerturbScratch();
      auto& scratch      = static_cast<DNCDMInvSpecies::Scratch&>(*scratch_owner);
      comp->ApplyKernelPerturbDiagonal(a, l_max_col, bg.data(), pL, fL, bL, diag.data(), scratch);

      // Independent reference: repeat the worker's gather, then take ONE
      // CollisionDiagonal straight off the perturbation kernel. bgm_ is null under
      // BuildCompositeDirect, so the worker's a'/a cap is 0 here — match it.
      std::vector<double> refH(qp, 0.), refl(qf, 0.), refb(qb, 0.);
      {
        const double kappa = comp->KappaStoredToBare();
        std::vector<double> fH(qp);
        for (int i = 0; i < qp; ++i)
          fH[i] = kappa * std::exp(bg[parent.bg_lnf_index() + i]);
        scratch.kernel_pt->PrepareTransitions(a,
                                              parent.GetMass(),
                                              parent.Gamma(),
                                              fH.data(),
                                              &bg[ferm.bg_f_index()],
                                              &bg[bos.bg_f_index()],
                                              0.);
        scratch.kernel_pt->CollisionDiagonal(refH.data(), refl.data(), refb.data());
      }

      // Magnitude before agreement: with inverse decays and quantum statistics on,
      // all three legs must carry a nonzero diagonal, or every assertion below is
      // vacuous. (Decay-only would zero g_l and g_phi and make the daughters trivial.)
      double magH = 0., magl = 0., magb = 0.;
      for (int i = 0; i < qp; ++i)
        magH = std::max(magH, std::fabs(refH[i]));
      for (int j = 0; j < qf; ++j)
        magl = std::max(magl, std::fabs(refl[j]));
      for (int k = 0; k < qb; ++k)
        magb = std::max(magb, std::fabs(refb[k]));
      assert(magH > 0.);
      assert(magl > 0.);
      assert(magb > 0.);

      const NCDMBaseSpecies::PerturbLayout* lay[3] = {&pL, &fL, &bL};
      const std::vector<double>* ref[3]            = {&refH, &refl, &refb};
      const int qn[3]                              = {qp, qf, qb};
      double worst_val                             = 0.;
      for (int c = 0; c < 3; ++c) {
        for (int iq = 0; iq < qn[c]; ++iq) {
          const double expect = (*ref[c])[iq];
          for (int l = 0; l <= l_max; ++l) {
            const int idx = lay[c]->index_per_q[iq] + l;
            if (l > l_max_col) {
              assert(diag[idx] == 0.);  // above l_max_col: never written
              continue;
            }
            worst_val = std::max(worst_val,
                                 std::fabs(diag[idx] - expect) / Scale(diag[idx], expect));
            // ℓ-independence AS WRITTEN: every populated slot of a bin holds the one
            // value, bit for bit — the whole reason one kernel pass serves the
            // hierarchy. Not a tolerance; the same double is stored to each slot.
            assert(diag[idx] == diag[lay[c]->index_per_q[iq]]);
          }
        }
      }

      std::printf(
          "  pert diagonal (refine=%d, l_max_col=%d of l_max=%d): max dev from "
          "CollisionDiagonal=%.2e; |diag| peaks H/l/phi = %.2e/%.2e/%.2e\n",
          refine,
          l_max_col,
          l_max,
          worst_val,
          magH,
          magl,
          magb);
      // Agreement, not reproducibility. The reference above re-evaluates the gather's
      // std::exp in a loop of its own, and -ffast-math lets the vectoriser swap either
      // loop for glibc's vector libm, which is a few ULP from the scalar call; the
      // kernel then amplifies that by O(1). Demanding bit-for-bit here measures the
      // libm the build happened to pick. The two exactness claims -- ell-independence
      // above and accumulation below -- compare the buffer against itself, so they
      // stay at zero ULP everywhere.
      assert(worst_val < 1e-12);

      // "Accumulate into diag; do not zero it" (base_species.h). A second pass over
      // the same buffer must double every written slot exactly and still leave the
      // truncated ones alone -- measured against pass one's own output, since that is
      // the only doubling the contract actually promises.
      const std::vector<double> once = diag;
      comp->ApplyKernelPerturbDiagonal(a, l_max_col, bg.data(), pL, fL, bL, diag.data(), scratch);
      for (int c = 0; c < 3; ++c)
        for (int iq = 0; iq < qn[c]; ++iq)
          for (int l = 0; l <= l_max; ++l) {
            const int idx = lay[c]->index_per_q[iq] + l;
            assert(diag[idx] == (l <= l_max_col ? 2. * once[idx] : 0.));
          }
    }
  }  // refine loop

  // ── Physical conservation (#385) ─────────────────────────────────────────────
  // The section above pins the kernel's BARE per-dof identities (exact regardless
  // of the f-state, by construction of the kernel). This section instead checks
  // the PHYSICAL (factor-weighted) identities using only public per-species
  // accessors (factor(), GetT()) -- no kappa, no (2*pi)^3 anywhere here -- at a
  // shooting-realistic deg_H = 0.015 (an under-occupied FD amplitude, the regime
  // bug B/#385 broke: physical conservation held only in bare-kernel units, not
  // in the factor-weighted units the rest of the code actually uses for rho/p).
  {
    // deg drives the parent's closure directly here (NOT SetInvBase's Omega_dncdmdr
    // combined mode -- DNCDMSpecies rejects specifying both, dncdm_species.cpp:176).
    FileContent fc;
    fc.set("dncdm1.type", "ncdm_decay_dr");
    fc.set("dncdm1.m", "1.0");
    fc.set("dncdm1.Gamma", "10");
    fc.set("dncdm1.quadrature_strategy", "3");  // qm_trapz (single grid: q == q_bg)
    fc.set("dncdm1.momenta_bins", "16");
    fc.set("dncdm1.dr_q_min", "1e-4");
    fc.set("dncdm1.dr_q_max", "300");
    fc.set("dncdm1.dr_N_q", "100");
    fc.set("dncdm1.inverse_decays", "yes");
    fc.set("dncdm1.quantum_statistics", "yes");
    fc.set("dncdm1.deg", "0.015");
    auto comp = BuildComposite(fc, pba, settings);

    int index_bg = 0, index_bi = 0;
    comp->RegisterBackgroundIndices(index_bg);
    comp->RegisterIntegrationIndices(index_bi);
    std::vector<double> y(index_bi, 0.), dy(index_bi, 0.);

    auto& parent = comp->parent();
    auto& ferm   = comp->fermion();
    auto& bos    = comp->boson();
    unsigned s   = 0x9E3779B9u;
    auto rnd     = [&]() {
      s = s * 1664525u + 1013904223u;
      return 0.1 + 0.9 * (((s >> 8) & 0xffffu) / 65535.0);
    };
    for (int i = 0; i < parent.q_size(); ++i)
      y[parent.bi_f_parent_index() + i] = DNCDMSpecies::kFScale * rnd();
    for (int i = 0; i < (int) ferm.q_bg().size(); ++i)
      y[ferm.bi_f_index() + i] = rnd();
    for (int i = 0; i < (int) bos.q_bg().size(); ++i)
      y[bos.bi_f_index() + i] = rnd();

    const double a = 1e-3;
    comp->ApplyKernelBackgroundDerivs(a, y.data(), dy.data());

    // Per-species raw grid moments of dy (parent: bg grid GetQ()/dq(), eps from
    // GetMass(); daughters: q_bg()/dq_bg(), eps == q). Public accessors only.
    double dN_H = 0., dE_H = 0.;
    const double a2m2 = a * a * parent.GetMass() * parent.GetMass();
    for (int i = 0; i < parent.q_size(); ++i) {
      const double q    = parent.GetQ()[i];
      const double eps  = std::sqrt(q * q + a2m2);
      const double w    = parent.dq()[i] * q * q;
      const double d    = dy[parent.bi_f_parent_index() + i] / DNCDMSpecies::kFScale;
      dN_H             += w * d;
      dE_H             += w * eps * d;
    }
    double dN_l = 0., dE_l = 0.;
    for (int i = 0; i < (int) ferm.q_bg().size(); ++i) {
      const double q  = ferm.q_bg()[i];
      const double w  = ferm.dq_bg()[i] * q * q;
      const double d  = dy[ferm.bi_f_index() + i];
      dN_l           += w * d;
      dE_l           += w * q * d;  // massless daughter: eps == q
    }
    double dN_phi = 0., dE_phi = 0.;
    for (int i = 0; i < (int) bos.q_bg().size(); ++i) {
      const double q  = bos.q_bg()[i];
      const double w  = bos.dq_bg()[i] * q * q;
      const double d  = dy[bos.bi_f_index() + i];
      dN_phi         += w * d;
      dE_phi         += w * q * d;
    }

    // PHYSICAL identities (factor_i in CLASS's 8piG/3 units, ~ deg_i*T_i^4):
    const double e_lhs = parent.factor() * dE_H + ferm.factor() * dE_l + bos.factor() * dE_phi;
    const double e_scale =
        Scale3(parent.factor() * dE_H, ferm.factor() * dE_l, bos.factor() * dE_phi);
    const double nf_lhs   = (parent.factor() / parent.GetT()) * dN_H +
                            (ferm.factor() / ferm.GetT()) * dN_l;
    const double nf_scale = Scale((parent.factor() / parent.GetT()) * dN_H,
                                  (ferm.factor() / ferm.GetT()) * dN_l);
    const double nb_lhs   = (parent.factor() / parent.GetT()) * dN_H +
                            (bos.factor() / bos.GetT()) * dN_phi;
    const double nb_scale = Scale((parent.factor() / parent.GetT()) * dN_H,
                                  (bos.factor() / bos.GetT()) * dN_phi);
    // The energy leg carries the lumped-loss residual (the number legs do NOT -- the
    // lumped split removes exactly the same total number, which is why every N
    // identity below is still asserted at zero). Converting it into these physical
    // units needs one constant, and the boson fixes it: the daughters' dy is already
    // bare, so bos.factor() multiplies the kernel's bare E_phi directly, and the
    // g_l/g_phi = 2 of the #385 convention carries the other two legs.
    const double split_phys = bos.factor() * comp->SplitEnergyResidual();
    std::printf(
        "  physical conservation (#385): E_rel=%.3e  N_fermion_rel=%.3e  N_boson_rel=%.3e  "
        "dN_H=%.3e  (booked split %.3e)\n",
        std::fabs(e_lhs - split_phys) / e_scale,
        std::fabs(nf_lhs) / nf_scale,
        std::fabs(nb_lhs) / nb_scale,
        dN_H,
        split_phys / e_scale);
    assert(split_phys != 0.);  // else the energy assertion below is vacuous
    assert(std::fabs(e_lhs - split_phys) < 1e-10 * e_scale);
    assert(std::fabs(nf_lhs) < 1e-10 * nf_scale);
    assert(std::fabs(nb_lhs) < 1e-10 * nb_scale);
    // Non-triviality: dN_H must be a non-negligible FRACTION of the raw moment
    // scale, not merely nonzero -- a bare "> 0." check is satisfied by a
    // uniformly-tiny (near-denormal) dy too, since Scale()'s 1e-300 floor then
    // makes every relative check above pass right alongside it.
    const double n_moment_scale = Scale3(dN_H, dN_l, dN_phi);
    assert(std::fabs(dN_H) > 1e-6 * n_moment_scale);
  }

  std::printf("dncdm inv test passed\n");
  return 0;
}
