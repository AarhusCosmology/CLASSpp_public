#include "dncdm_inv_species.h"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <future>
#include <memory>
#include <set>
#include <string>
#include <tuple>

#include "background_module.h"
#include "errors.h"
#include "perturbations_module.h"
#include "reduced_collision_operator.h"
#include "species/shooting_target.h"
#include "species/species_build_context.h"
#include "species/species_input.h"
#include "thread_pool.h"  // Tools::TaskSystem — the precompute is row-parallel

namespace {

/** The reduced scheme is on and `dr_table_max_dlna` was left at a default that has only
 *  ever been validated to Gamma = 1e9. Say so, once per instance.
 *
 *  This warns; it does not size anything. Deriving the spacing from Gamma here would be
 *  one line and is deliberately not done — see kTableMaxDlnaDefault. What the caller
 *  cannot be allowed to do is walk into it unaware, because the failure is silent: the
 *  identities the kernel asserts are number identities, and a mis-interpolated operator
 *  satisfies them all the way into a 4e10% runaway. */
void WarnTableSpacingUnset(const std::string& name, double gamma_kms_mpc) {
  if (!(gamma_kms_mpc > DNCDMInvSpecies::kTableDefaultValidatedTo))
    return;
  // Keyed on the instance and Gamma, so a re-Compute from the python wrapper with a new
  // cosmology warns again. Species construction is single-threaded, so no lock is needed.
  static std::set<std::pair<std::string, double>> warned;
  if (!warned.emplace(name, gamma_kms_mpc).second)
    return;
  fprintf(stderr,
          "WARNING: species '%s': dr_reduced_moments is set and dr_table_max_dlna is not, so "
          "the reduced collision table is interpolated at the default %g. That default is "
          "validated only to Gamma = %g and this run is at Gamma = %g; the table's accuracy "
          "requirement grows with Gamma and a fixed spacing does not follow it. Measured "
          "C_2^TT error at this default, m = 0.3 eV: -2.4%% at Gamma = 1e10, -30%% at 1e11, "
          "a full runaway by 1e12 -- and nothing else in the run reports it. Set "
          "dr_table_max_dlna explicitly; the sizing rule is in dncdm_inv_species.h.\n",
          name.c_str(),
          DNCDMInvSpecies::kTableMaxDlnaDefault,
          DNCDMInvSpecies::kTableDefaultValidatedTo,
          gamma_kms_mpc);
}

}  // namespace

// ─────────────────────────────────────────────────────────────────────────────
// Factory
// ─────────────────────────────────────────────────────────────────────────────

Named DNCDMInvSpecies::Create(std::unique_ptr<DNCDMSpecies> parent,
                              const SpeciesBuildContext& ctx) {
  const std::string name = parent->name();
  SpeciesInput in(ctx.pfc, name);

  // Same default as the dispatch in DNCDM_DR_Species::CreateAll, which is what
  // routed the instance here. `no` is the decay-only rung: the same pure-decay
  // physics as the `integrated` representation but through the resolved-PSD kernel,
  // which is the only independent check this sector has.
  const bool inv = in.get_flag("inverse_decays", false);
  // Defaults to whatever inverse_decays is, so one key selects the full model or the
  // decay-only rung and the two cannot silently disagree.
  const bool qs = in.get_flag("quantum_statistics", inv);

  // quantum_statistics only makes sense with inverse decays present: the linear
  // (1±f) term is f_H(f_l−f_φ), meaningful only alongside the f_l f_φ inverse
  // term. Structural flag-combination check (mirrors the T_opt/gstar_opt XOR in
  // axion_ncdm_species.cpp:22-25) — decidable from which keys are set, not from
  // any numeric value a sampler could vary — so severe, not rejectable.
  class_test_severe(qs && !inv,
                    "species '%s': quantum_statistics requires inverse_decays = yes (the (1±f) "
                    "terms are only meaningful with both decay and inverse present)",
                    name.c_str());

  // Synchronous gauge only (mirror dcdm_wdm): any gauge starting with "new".
  if (auto gauge = ctx.pfc->get<std::string>("gauge")) {
    std::string g = *gauge;
    std::transform(g.begin(), g.end(), g.begin(), [](unsigned char c) { return std::tolower(c); });
    class_test_severe(g.rfind("new", 0) == 0,
                      "species '%s': inverse decays support the synchronous gauge only "
                      "(arXiv:2011.01502)",
                      name.c_str());
  }

  // No tensor modes with inverse decays (the tensor collision operator is out of
  // scope for v1); "modes" containing t/T requests tensors.
  if (auto modes = ctx.pfc->get<std::string>("modes")) {
    const bool has_tensor = modes->find('t') != std::string::npos ||
                            modes->find('T') != std::string::npos;
    class_test_severe(has_tensor,
                      "species '%s': inverse decays do not support tensor modes (scalars only)",
                      name.c_str());
  }

  // No NCDM fluid approximation for the coupled trio: the kernel collision
  // operator needs the full Boltzmann hierarchy on every species. Unsupported-
  // feature-combination check (mirrors the use_psd_file guard in
  // axion_ncdm_species.cpp:44) — severe, not a numeric-value rejection.
  class_test_severe(in.get<std::string>("fluid_approximation").has_value(),
                    "species '%s': inverse decays require the full Boltzmann hierarchy; remove "
                    "'%s.fluid_approximation'",
                    name.c_str(),
                    name.c_str());

  // The parent must be a single momentum grid: its f-per-bin variable and its
  // ComputeMomenta integration share one grid (design §5). qm_auto with differing
  // tol_ncdm/tol_ncdm_bg breaks this — pin momenta_bins_bg == momenta_bins.
  class_test_severe(parent->q_size() != parent->q_size_bg(),
                    "species '%s': inverse decays require a single parent momentum grid "
                    "(q_size=%d != q_size_bg=%d); set momenta_bins_bg == momenta_bins (or a "
                    "manual quadrature_strategy, or tol_ncdm == tol_ncdm_bg)",
                    name.c_str(),
                    parent->q_size(),
                    parent->q_size_bg());

  // Kernel boundary (#385): KappaStoredToBare assumes the parent's stored f is a
  // single deg_H-suppressed Fermi-Dirac; a file-based PSD (arbitrary shape) has no
  // such reading. Structural (key-presence) check, mirrors the fluid_approximation
  // guard above.
  class_test_severe(parent->UsesPsdFile(),
                    "species '%s': inverse decays require an analytic (non-file) parent PSD; "
                    "remove '%s.use_psd_file' (the kernel boundary's bare-occupation reading "
                    "assumes a single Fermi-Dirac shape, #385)",
                    name.c_str(),
                    name.c_str());

  // Nonzero ksi: the stored f_H = FD(q-ksi)+FD(q+ksi) sum has no single bare-
  // occupation reading (the kappa conversion assumes ONE FD shape). Numeric
  // (sampler-varyable) check, not structural.
  class_test(parent->GetKsi() != 0.,
             "species '%s': inverse decays require the parent's chemical potential ksi = 0 "
             "(FD(q-ksi)+FD(q+ksi) has no single bare-occupation reading at the kernel "
             "boundary, #385)",
             name.c_str());

  // dr_T is pinned to the parent's T for kernel-unit consistency (#385): every
  // species must report occupation at the SAME dimensionless q = p/T at the kernel
  // boundary (see the DNCDMInvSpecies ctor's PinTemperature calls). Numeric
  // (sampler-varyable) check: an explicit mismatched dr_T rejects the point rather
  // than aborting the run. Note the legacy side effect (dr_psd_species.cpp's f_ini_
  // fallback): setting dr_T at all — even to the value this guard requires — seeds
  // BOTH daughters at full equilibrium, not just a pinned temperature. To seed the
  // daughters ASYMMETRICALLY (the published A = [0, 1, 1]: phi empty, nu_l full
  // Fermi-Dirac) use the per-daughter amplitudes dr_f_ini_phi / dr_f_ini_l, which
  // override that fallback; the combined-mode default (no key at all) starts both
  // daughters empty.
  if (auto dr_T = in.get<double>("dr_T")) {
    class_test(std::fabs(*dr_T - parent->GetT()) > 1e-12 * std::fabs(parent->GetT()),
               "species '%s': dr_T (=%g) must equal the parent's T (=%g); the daughters' "
               "temperature is pinned to the parent's for kernel-unit consistency (#385). Note: "
               "setting dr_T at all also seeds BOTH daughters with a thermal initial abundance "
               "-- omit it unless that is intended, and use dr_f_ini_l / dr_f_ini_phi to set the "
               "two daughters' initial occupations independently",
               name.c_str(),
               *dr_T,
               parent->GetT());
  }

  DecayTransitionKernel::Config cfg;
  cfg.inverse_decays     = inv;
  cfg.quantum_statistics = qs;
  // Well-balanced band factor (see Config::balanced_gather). Only meaningful with
  // inverse decays: it exists to stop the linear gather's convexity bias reading as a
  // spurious f_l f_phi repopulation, and on the decay-only rung there is no such term,
  // so the gather never enters Lambda and the flag is a no-op. Defaulting it to `inv`
  // would silently change the reference method's published numbers, so it defaults off
  // and is opted into.
  cfg.balanced_gather = in.get_flag("balanced_gather", false);
  class_test_severe(cfg.balanced_gather && !inv,
                    "species '%s': balanced_gather requires inverse_decays = yes (it corrects the "
                    "f_l f_phi band factor, which the decay-only rung does not have)",
                    name.c_str());
  // The lumped daughter loss buys positivity from the LINEAR gather's overestimate at
  // an injection front. The balanced gather removes that overestimate at the root -- a
  // geometric mean vanishes with any factor, so an empty bin is billed exactly nothing
  // -- which makes the exact split Metzler again and lets the lumping's O(dq) energy
  // misplacement go. That misplacement is the dominant term in the daughter-grid
  // requirement at high Gamma, so the two flags are one scheme and this is where they
  // are coupled. Overridable: `lumped_loss` set explicitly wins, which is what makes
  // the 2x2 measurable from an ini.
  cfg.lumped_loss = in.get_flag("lumped_loss", !cfg.balanced_gather);
  class_test_severe(!cfg.lumped_loss && !cfg.balanced_gather,
                    "species '%s': lumped_loss = no requires balanced_gather = yes. Without it the "
                    "exact two-bin split has a negative off-diagonal, so an empty daughter bin is "
                    "billed for its neighbour's particles and the right-hand side drives it "
                    "negative -- no integrator can repair that",
                    name.c_str());
  // The balanced gather's CHAIN RULE in the perturbation operator (see
  // Config::balanced_pert). Without it the operator is a consistent linearisation of the
  // LINEAR-gather band factor while the background integrates the balanced one -- a
  // 1.1e-1 drift from a finite-difference Jacobian, against 8.8e-12 on the linear path.
  // Off by default because it changes measured perturbation output; requires the
  // background flag, since there is no chain factor without it.
  cfg.balanced_pert = in.get_flag("balanced_pert", false);
  class_test_severe(cfg.balanced_pert && !cfg.balanced_gather,
                    "species '%s': balanced_pert = yes requires balanced_gather = yes (it is the "
                    "chain rule OF that gather; with a linear band factor there is nothing to "
                    "differentiate)",
                    name.c_str());
  cfg.chain_cap = in.get_or("chain_cap", cfg.chain_cap);
  class_test_severe(!(cfg.chain_cap > 1.),
                    "species '%s': chain_cap must exceed 1 (got %g) -- 1 is the linear "
                    "gather's own sensitivity",
                    name.c_str(),
                    cfg.chain_cap);
  // Emission-quadrature order; see Config::emission_gauss. The default reproduces the
  // historical midpoint rule exactly.
  cfg.emission_gauss = in.get_or("emission_gauss", cfg.emission_gauss);
  class_test_severe(cfg.emission_gauss < 1 || cfg.emission_gauss > 4,
                    "species '%s': emission_gauss must be 1..4 (got %d)",
                    name.c_str(),
                    cfg.emission_gauss);
  // No stiffness cap: this is the reference method, and capping it is an
  // approximation of unquantified accuracy. DNCDMProxySpecies, whose whole point is
  // cheapness, does cap -- see its dr_rate_cap.

  // Daughter grid sizing. The emission band's top is closed-form (see
  // MaxDaughterMomentum) and SATURATES at a = 1, so the daughters can be sized once,
  // here, from the parent that will feed them. This is the composite's job and not
  // the daughter's: DrPsdSpecies is also a standalone species, and a daughter with no
  // decay channel feeding it has no kinematic requirement to meet.
  const double q_star_max = DecayTransitionKernel::MaxDaughterMomentum(parent->GetQ().back(),
                                                                       parent->GetMass());

  // ONE class plays both daughters — statistics passed explicitly (the composite
  // owns the pairing; it does not parse a statistics key off the instance).
  auto fermion = std::make_unique<DrPsdSpecies>(ctx.pfc,
                                                name,
                                                *ctx.ncdm_settings,
                                                ctx.pba,
                                                ctx.bgm,
                                                Statistics::Fermion,
                                                q_star_max);
  auto boson   = std::make_unique<DrPsdSpecies>(ctx.pfc,
                                                name,
                                                *ctx.ncdm_settings,
                                                ctx.pba,
                                                ctx.bgm,
                                                Statistics::Boson,
                                                q_star_max);

  // An explicit dr_q_max is always honoured — a narrow grid is sometimes a deliberate
  // cost choice, and refusing it would make published configurations unreproducible.
  // But it is worth exactly one warning, because the failure mode is silent: the
  // deposit is CLAMPED into the top bin rather than dropped, so number conservation
  // (which is what the kernel's exact identities test) stays perfect while the sector
  // loses a factor q_max/q* of every decay's energy. Emitted here rather than in
  // DrPsdSpecies so it fires once for the pair, not once per daughter.
  //
  // Latched to once per RUN, not once per Create: the shooter rebuilds the whole
  // species stack at every Newton iteration, so an unlatched warning printed three
  // times for a single run. Keyed on the instance AND the two momenta, so a genuinely
  // different configuration (a second instance, or a re-Compute from the python
  // wrapper with a new cosmology) still gets its own warning. Species construction is
  // single-threaded (input + background setup), so the bare static needs no lock.
  const double q_max_used = fermion->q_bg().back();
  static std::set<std::tuple<std::string, double, double>> warned;
  if (q_max_used < q_star_max && warned.emplace(name, q_max_used, q_star_max).second) {
    fprintf(stderr,
            "WARNING: species '%s': dr_q_max (=%g) is below the highest momentum this decay can "
            "emit (q* = %g at a = 1, from the parent's q_max = %g and M = %g). Off-band deposits "
            "are CLAMPED into the top daughter bin, which conserves number but files the energy "
            "at q_max instead of at q*, so the decay sector silently loses roughly a factor "
            "q_max/q* of the energy released once a*M/2 outgrows the grid. Set dr_q_max >= %g, "
            "or omit it entirely to get that value by default.\n",
            name.c_str(),
            q_max_used,
            q_star_max,
            parent->GetQ().back(),
            parent->GetMass(),
            DrPsdSpecies::kQMaxMargin * q_star_max);
  }

  // DNCDMSpecies stores Gamma internally as Gamma/c in 1/Mpc (dncdm_species.cpp:121).
  // Everything that talks ABOUT it -- the input key, the campaign, the sizing rule in
  // kTableMaxDlnaDefault -- speaks km/s/Mpc, so convert back once, here.
  const double gamma_kms_mpc = parent->Gamma() * (_c_ / 1.e3);

  auto composite = std::make_unique<DNCDMInvSpecies>(std::move(parent),
                                                     std::move(fermion),
                                                     std::move(boson),
                                                     cfg,
                                                     ctx.pba,
                                                     ctx.bgm);
  // Reduced-moment daughters (0 = off, the shipped default). See
  // DrPsdSpecies::SetReducedBasis for what the representation is and why its basis is
  // frozen; 6 is the measured working point (closure p90 <= 0.29 at Gamma = 1e7 and 1e8,
  // and frozen-at-6 is as accurate as live-at-4).
  const int red = in.get_or("dr_reduced_moments", 0);
  class_test_severe(red < 0 || red > DrPsdSpecies::kMaxReducedMoments,
                    "species '%s': dr_reduced_moments (=%d) must be 0 (off) .. %d",
                    name.c_str(),
                    red,
                    DrPsdSpecies::kMaxReducedMoments);
  // SEVERE, like the range check above it: dr_reduced_moments names a discretisation, not
  // a cosmological parameter, so a sampler can never vary its way out of this. class_test
  // here would raise CosmoComputationError and merely REJECT THE POINT, which in a
  // shooting or MCMC run means the message is swallowed and the fit quietly loses samples.
  class_test_severe(red == 1,
                    "species '%s': dr_reduced_moments = 1 cannot represent an equilibrium "
                    "perturbation (that space is two-dimensional: one temperature and one "
                    "chemical potential), so it breaks detailed balance by construction. Use "
                    ">= 2, and see the design note for why 6 is the measured working point",
                    name.c_str());
  composite->set_reduced_moments(red);
  // Reduced collision table ln-a spacing. See kTableMaxDlnaDefault for the measured
  // sizing rule and for why the Gamma dependence is the caller's to apply.
  const auto dlna_set = in.get<double>("dr_table_max_dlna");
  const double dlna   = dlna_set.value_or(DNCDMInvSpecies::kTableMaxDlnaDefault);
  class_test_severe(!(dlna > 0.),
                    "species '%s': dr_table_max_dlna (=%g) must be positive -- it is a ln-a "
                    "spacing, and the table sub-divides each background row until its own "
                    "spacing meets it",
                    name.c_str(),
                    dlna);
  composite->set_table_max_dlna(dlna);
  if (red > 0 && !dlna_set)
    WarnTableSpacingUnset(name, gamma_kms_mpc);
  composite->background_verbose_ = ctx.pfc->get_or("background_verbose", 0);
  // Resolved here and not in ProcessBackgroundTable, which runs before any perturbation
  // layout exists.  Mirrors what AddCouplingDerivs computes per call:
  // min(l_max_dncdm_col, the hierarchies' own l_max).
  if (ctx.ppr != nullptr)
    composite->set_reduced_table_l_max(std::min(ctx.ppr->l_max_dncdm_col, ctx.ppr->l_max_ncdm));
  return {name, std::move(composite)};
}

// ─────────────────────────────────────────────────────────────────────────────
// Constructor
// ─────────────────────────────────────────────────────────────────────────────

DNCDMInvSpecies::DNCDMInvSpecies(std::unique_ptr<DNCDMSpecies> parent,
                                 std::unique_ptr<DrPsdSpecies> fermion,
                                 std::unique_ptr<DrPsdSpecies> boson,
                                 DecayTransitionKernel::Config cfg,
                                 const background* pba,
                                 const BackgroundModule* bgm)
    : CompositeSpecies(parent->name(), BaseSpecies::EnergyType::Other), pba_(pba), bgm_(bgm) {
  parent_  = parent.get();
  fermion_ = fermion.get();
  boson_   = boson.get();
  cfg_     = cfg;
  parent_->SetCollisionOwned(true);  // before any index registration (design §3 table)

  // Pin both daughters' T to the parent's (#385 kernel-unit consistency): the
  // kappa boundary conversion (KappaStoredToBare) and the kernel's shared-q
  // convention only hold if every species reports occupation at the SAME
  // dimensionless q = p/T. Direct-built test composites (BuildCompositeDirect)
  // must get this too, so it lives here rather than in Create's guards.
  fermion_->PinTemperature(parent_->GetT());
  boson_->PinTemperature(parent_->GetT());

  // Grid views (single-grid parent: q_/dq_; daughters own their log-trapezoid
  // background grids). The pointed-to buffers outlive the kernel (owned by the
  // children, never resized after construction).
  const GridView parent_view{parent_->GetQ().data(), parent_->dq().data(), parent_->q_size()};
  const GridView fermion_view{fermion_->q_bg().data(),
                              fermion_->dq_bg().data(),
                              static_cast<int>(fermion_->q_bg().size())};
  const GridView boson_view{boson_->q_bg().data(),
                            boson_->dq_bg().data(),
                            static_cast<int>(boson_->q_bg().size())};
  // cfg_, not the cfg parameter: the two are equal except that cfg_ carries the
  // diagnostic overrides read above, and CreatePerturbScratch already builds the per-k
  // perturbation kernels from cfg_. Both legs must see ONE Config, or f-bar and delta-f
  // evolve under different effective rates -- i.e. two different physical models.
  kernel_ = std::make_unique<DecayTransitionKernel>(parent_view,
                                                    fermion_view,
                                                    boson_view,
                                                    fermion_->statistics(),
                                                    boson_->statistics(),
                                                    cfg_);

  // The PERTURBATION kernels are NOT built here. They are per-k, built by
  // CreatePerturbScratch on these same grids and Config, because the kernel caches
  // one (a, f-bar) transition geometry and a rolling Legendre recurrence -- state
  // that belongs to a single k-mode integration. Keeping a shared one here is what
  // made threads>1 abort with the operator's strict-ascending-l contract.
  // Child order is the child_layouts contract: kParent, kFermion, kBoson (#358).
  children_.push_back(std::move(parent));
  children_.push_back(std::move(fermion));
  children_.push_back(std::move(boson));

  df_H_.assign(parent_->q_size(), 0.);
  df_l_.assign(fermion_->q_bg().size(), 0.);
  df_phi_.assign(boson_->q_bg().size(), 0.);
  fH_bare_.assign(parent_->q_size(), 0.);
  diag_H_.assign(parent_->q_size(), 0.);
  diag_l_.assign(fermion_->q_bg().size(), 0.);
  diag_phi_.assign(boson_->q_bg().size(), 0.);

  // The perturbation gather/scatter buffers are NOT members any more: they are
  // per-k, allocated by CreatePerturbScratch and owned by the workspace.
}

void DNCDMInvSpecies::SetBackgroundModule(const BackgroundModule* bgm) {
  bgm_ = bgm;
  CompositeSpecies::SetBackgroundModule(bgm);
}

// ─────────────────────────────────────────────────────────────────────────────
// Background RHS: kernel source into the parent+daughter f-slots
// ─────────────────────────────────────────────────────────────────────────────

void DNCDMInvSpecies::ApplyKernelBackgroundDerivs(double a,
                                                  const double* y,
                                                  double* dy,
                                                  double a_prime_over_a) {
  const double m     = parent_->GetMass();
  const double Gamma = parent_->Gamma();
  const int Np       = parent_->q_size();
  // Kernel boundary (#385): the kernel's (1±f) coefficients and cubic-cancelling
  // band factor physically require BARE per-dof occupation; the parent's f-slot
  // is STORED occupation (deg_H-suppressed). Convert in (kappa*), then the
  // write-back converts the kernel's bare rate back to stored units (/kappa) —
  // the daughters already store bare occupation, so their legs are untouched.
  //
  // The parent's INTEGRATION slot additionally carries DNCDMSpecies::kFScale (the
  // published ln f column does not), so this leg — and only this leg — folds that
  // in. KappaStoredToBare() itself stays column-semantics, because that is what the
  // perturbation gather and the tests read.
  const double kappa = KappaStoredToBare() / DNCDMSpecies::kFScale;
  for (int i = 0; i < Np; ++i)
    fH_bare_[i] = kappa * y[parent_->bi_f_parent_index() + i];
  kernel_->ComputeBackgroundDerivs(a,
                                   m,
                                   Gamma,
                                   fH_bare_.data(),
                                   &y[fermion_->bi_f_index()],
                                   &y[boson_->bi_f_index()],
                                   df_H_.data(),
                                   df_l_.data(),
                                   df_phi_.data(),
                                   a_prime_over_a);
  // Assign (=) the kernel source: the parent's decay-only diagonal term is not
  // registered in collision-owned mode, and massless daughters have no comoving
  // dilution term, so the kernel source IS the full f-RHS for all three.
  for (int i = 0; i < Np; ++i)
    dy[parent_->bi_f_parent_index() + i] = df_H_[i] / kappa;
  const int Nf = static_cast<int>(fermion_->q_bg().size());
  for (int i = 0; i < Nf; ++i)
    dy[fermion_->bi_f_index() + i] = df_l_[i];
  const int Nb = static_cast<int>(boson_->q_bg().size());
  for (int i = 0; i < Nb; ++i)
    dy[boson_->bi_f_index() + i] = df_phi_[i];
}

void DNCDMInvSpecies::BackgroundDerivs(double tau,
                                       const double* y,
                                       double* dy,
                                       const double* pvecback) {
  CompositeSpecies::BackgroundDerivs(tau, y, dy, pvecback);  // parent/daughters: no-op sources
  // a'/a drives the kernel's stiffness cap; the background and perturbation legs
  // MUST pass the same quantity or f-bar and delta-f would evolve under different
  // effective rates, i.e. two different physical models.
  const double a_bg = pvecback[bgm_->index_bg_a_];
  ApplyKernelBackgroundDerivs(a_bg, y, dy, a_bg * pvecback[bgm_->index_bg_H_]);
}

void DNCDMInvSpecies::BackgroundDerivsDiagonal(double tau,
                                               const double* y,
                                               double* diag,
                                               const double* pvecback) {
  CompositeSpecies::BackgroundDerivsDiagonal(tau, y, diag, pvecback);
  const double a_bg = pvecback[bgm_->index_bg_a_];

  // The kernel caches its transitions from the last PrepareTransitions, and the
  // evolver always calls the derivs before the diagonal at the same state, so this
  // could read the cache. It re-prepares anyway: a diagonal silently taken from a
  // DIFFERENT state than the RHS it is paired with is exactly the failure mode that
  // stays stable while converging to the wrong attractor, and the guarantee is worth
  // more than one kernel pass per step (the evolver takes several RHS evaluations
  // per step regardless).
  const double m     = parent_->GetMass();
  const double Gamma = parent_->Gamma();
  const int Np       = parent_->q_size();
  const double kappa = KappaStoredToBare() / DNCDMSpecies::kFScale;
  for (int i = 0; i < Np; ++i)
    fH_bare_[i] = kappa * y[parent_->bi_f_parent_index() + i];
  kernel_->ComputeBackgroundDerivs(a_bg,
                                   m,
                                   Gamma,
                                   fH_bare_.data(),
                                   &y[fermion_->bi_f_index()],
                                   &y[boson_->bi_f_index()],
                                   df_H_.data(),
                                   df_l_.data(),
                                   df_phi_.data(),
                                   a_bg * pvecback[bgm_->index_bg_H_]);
  kernel_->CollisionDiagonal(diag_H_.data(), diag_l_.data(), diag_phi_.data());

  // d(dy_i)/d(y_i) is INVARIANT under the parent's stored/bare rescaling: the write
  // back divides the rate by the same kappa the read multiplies the occupation by,
  // so the two cancel. The daughters store bare occupation and need no conversion.
  for (int i = 0; i < Np; ++i)
    diag[parent_->bi_f_parent_index() + i] += diag_H_[i];
  const int Nf = static_cast<int>(fermion_->q_bg().size());
  for (int i = 0; i < Nf; ++i)
    diag[fermion_->bi_f_index() + i] += diag_l_[i];
  const int Nb = static_cast<int>(boson_->q_bg().size());
  for (int i = 0; i < Nb; ++i)
    diag[boson_->bi_f_index() + i] += diag_phi_[i];
}

DecayTransitionKernel::Moments DNCDMInvSpecies::ConservationMoments(double a,
                                                                    const double* dy) const {
  // dy's parent leg is in the integration slot's units (ApplyKernelBackgroundDerivs
  // divides the kernel's bare rate by kappa AND by kFScale on write-back);
  // ComputeMoments needs the kernel's own BARE units, so re-scale a scratch copy by
  // the same factor here (#385 kappa, plus the kFScale unit change). Test-only (see
  // the header doc), so a local allocation is fine.
  const double kappa = KappaStoredToBare() / DNCDMSpecies::kFScale;
  const int Np       = parent_->q_size();
  std::vector<double> dyH_bare(Np);
  for (int i = 0; i < Np; ++i)
    dyH_bare[i] = kappa * dy[parent_->bi_f_parent_index() + i];
  return kernel_->ComputeMoments(a,
                                 parent_->GetMass(),
                                 dyH_bare.data(),
                                 &dy[fermion_->bi_f_index()],
                                 &dy[boson_->bi_f_index()]);
}

void DNCDMInvSpecies::BuildReducedBasis(const double* background_table,
                                        int n_rows,
                                        int row_stride) {
  const int Np = parent_->q_size(), Nf = fermion_->q_size(), Nb = boson_->q_size();
  const GridView parent_view{parent_->GetQ().data(), parent_->dq().data(), Np};
  const GridView fermion_view{fermion_->q().data(), fermion_->dq().data(), Nf};
  const GridView boson_view{boson_->q().data(), boson_->dq().data(), Nb};
  const GridView fermion_bg_view{fermion_->q_bg().data(),
                                 fermion_->dq_bg().data(),
                                 static_cast<int>(fermion_->q_bg().size())};
  const GridView boson_bg_view{boson_->q_bg().data(),
                               boson_->dq_bg().data(),
                               static_cast<int>(boson_->q_bg().size())};
  const double m     = parent_->GetMass();
  const double Gamma = parent_->Gamma();
  const double kappa = KappaStoredToBare();
  const int pf       = parent_->bg_lnf_index();
  const int fbase    = fermion_->bg_f_index();
  const int bbase    = boson_->bg_f_index();

  reduced_basis_kernel_        = std::make_unique<DecayTransitionKernel>(parent_view,
                                                                         fermion_view,
                                                                         boson_view,
                                                                         fermion_->statistics(),
                                                                         boson_->statistics(),
                                                                         cfg_,
                                                                         fermion_bg_view,
                                                                         boson_bg_view);
  DecayTransitionKernel& probe = *reduced_basis_kernel_;
  std::vector<double> fH(Np), diagH(Np), diagl(Nf), diagp(Nb);

  // Find the DECAY EPOCH: the row where the parent's comoving number falls fastest.
  //
  // NOT the row where the collision rate is largest against aH: the parent's diagonal
  // tends to a*Gamma while aH falls, so rate/aH grows monotonically and peaks at a = 1
  // -- long after the parent is extinct (rho_H/rho_tot ~ 1e-60) and the collision has
  // nothing left to act on. That criterion freezes the basis at a = 1 on an empty
  // daughter spectrum. n_H falls through 1/e where the decay actually happens, which
  // is the window the reduction has to be accurate in.
  auto comoving_n = [&](const double* bg) {
    double n = 0.;
    for (int i = 0; i < Np; ++i) {
      fH[i]  = kappa * std::exp(bg[pf + i]);
      n     += parent_->dq()[i] * parent_->GetQ()[i] * parent_->GetQ()[i] * fH[i];
    }
    return n;
  };
  const double n0 = comoving_n(background_table);
  class_test_severe(!(n0 > 0.),
                    "species '%s': dr_reduced_moments is set but the parent's comoving number "
                    "is zero at the first background row",
                    name().c_str());

  // The first row where the parent has decayed to 1/e of its initial comoving number.
  // Deliberately a THRESHOLD CROSSING, not an extremum of any rate: both rate-based
  // criteria tried here picked late times, where the parent's occupation sits on
  // kFParentFloor and every log-derivative of it is floor noise (the run warns about
  // exactly that). n_H falls monotonically through 1/e once and only once.
  int best_row = -1;
  for (int row = 0; row < n_rows; ++row) {
    if (comoving_n(background_table + (size_t) row * row_stride) < n0 / _E_) {
      best_row = row;
      break;
    }
  }
  class_test_severe(best_row < 0,
                    "species '%s': dr_reduced_moments is set but the parent never decays to 1/e "
                    "of its initial comoving number, so there is no decay epoch to freeze a "
                    "basis at. Is Gamma too small for this run's a-range?",
                    name().c_str());
  const double best_rate = n0;

  const double* bg = background_table + (size_t) best_row * row_stride;
  const double a   = bg[bgm_->index_bg_a_];
  for (int i = 0; i < Np; ++i)
    fH[i] = kappa * std::exp(bg[pf + i]);
  // The daughters' PSDs on their STATE grid (pvecback carries them on the finer
  // background grid; bg_index is the exact subsample).
  std::vector<double> fl(Nf), fphi(Nb);
  for (int j = 0; j < Nf; ++j)
    fl[j] = bg[fbase + fermion_->bg_index(j)];
  for (int j = 0; j < Nb; ++j)
    fphi[j] = bg[bbase + boson_->bg_index(j)];

  reduced_op_ = std::make_unique<ReducedCollisionOperator>(probe,
                                                           parent_view,
                                                           fermion_view,
                                                           boson_view,
                                                           fermion_->statistics(),
                                                           boson_->statistics(),
                                                           ReducedCollisionOperator::Config{
                                                               reduced_moments_});
  reduced_op_->SetBackground(a, m, fH.data(), fl.data(), fphi.data());

  fermion_->SetReducedBasis(reduced_moments_,
                            reduced_op_->psi(0),
                            reduced_op_->gnorm(0),
                            reduced_op_->q_ref(0),
                            reduced_op_->alpha0(0));
  boson_->SetReducedBasis(reduced_moments_,
                          reduced_op_->psi(1),
                          reduced_op_->gnorm(1),
                          reduced_op_->q_ref(1),
                          reduced_op_->alpha0(1));

  if (background_verbose_ > 1)
    fprintf(stderr,
            "[dncdm-reduced] %s: %d moments per daughter, basis frozen at a = %.4e "
            "(decay epoch: n_H = n_ini/e, n_ini = %.4g); daughter perturbation variables "
            "%d -> %d per multipole\n",
            name().c_str(),
            reduced_moments_,
            a,
            best_rate,
            Nf + Nb,
            2 * reduced_moments_);

  // The basis weight is f-bar(1 -/+ f-bar), which is only a weight while the occupation
  // stays in range. The gather does not enforce Pauli blocking and DrPsdSpecies floors f
  // but never caps it at 1, so a fermion bin CAN come back above saturation -- and a
  // negative weight is silent: it flips the sign of the reconstruction across the
  // saturated core without failing anything. ReducedCollisionOperator clamps it (see
  // EntropyWeight); this is the line that says whether the clamp ever had to.
  const char* dname[2] = {"nu_l", "phi"};
  for (int d = 0; d < 2; ++d) {
    const auto c = reduced_op_->weight_clamp(d);
    if (c.n_clamped > 0)
      fprintf(stderr,
              "[dncdm-reduced] %s: WARNING -- %s occupation outside [0, 1] in %d of %d bins "
              "at the freeze row (peak f = %.6f). The entropy weight is clamped to "
              "saturation there, so those bins carry no basis support. This is a "
              "background-solve defect, not a reduction one.\n",
              name().c_str(),
              dname[d],
              c.n_clamped,
              (d == 0) ? Nf : Nb,
              c.f_max);
    else if (background_verbose_ > 1)
      fprintf(stderr,
              "[dncdm-reduced] %s: %s peak occupation at the freeze row = %.6f, no clamping\n",
              name().c_str(),
              dname[d],
              c.f_max);
  }
}

void DNCDMInvSpecies::BuildReducedTable(const double* background_table,
                                        int n_rows,
                                        int row_stride) {
  const int Np = parent_->q_size(), Nf = fermion_->q_size(), Nb = boson_->q_size();
  const GridView parent_view{parent_->GetQ().data(), parent_->dq().data(), Np};
  const GridView fermion_view{fermion_->q().data(), fermion_->dq().data(), Nf};
  const GridView boson_view{boson_->q().data(), boson_->dq().data(), Nb};
  const GridView fermion_bg_view{fermion_->q_bg().data(),
                                 fermion_->dq_bg().data(),
                                 static_cast<int>(fermion_->q_bg().size())};
  const GridView boson_bg_view{boson_->q_bg().data(),
                               boson_->dq_bg().data(),
                               static_cast<int>(boson_->q_bg().size())};
  const double m     = parent_->GetMass();
  const double Gamma = parent_->Gamma();
  const double kappa = KappaStoredToBare();
  const int pf       = parent_->bg_lnf_index();
  const int fbase    = fermion_->bg_f_index();
  const int bbase    = boson_->bg_f_index();
  const int l_max    = reduced_table_l_max_;
  const int S        = Np + 2 * reduced_moments_;

  // Pass 1: bound the window, by the kernel's own diagonal against the expansion rate.
  // Only the FIRST and LAST row are needed, because the table has to be contiguous in a.
  DecayTransitionKernel probe(parent_view,
                              fermion_view,
                              boson_view,
                              fermion_->statistics(),
                              boson_->statistics(),
                              cfg_,
                              fermion_bg_view,
                              boson_bg_view);
  std::vector<double> fH(Np), diagH(Np), diagl(Nf), diagp(Nb);
  int row_lo = -1;
  for (int row = 0; row < n_rows; ++row) {
    const double* bg = background_table + (size_t) row * row_stride;
    const double a   = bg[bgm_->index_bg_a_];
    const double aH  = a * bg[bgm_->index_bg_H_];
    for (int i = 0; i < Np; ++i)
      fH[i] = kappa * std::exp(bg[pf + i]);
    probe.PrepareTransitions(a, m, Gamma, fH.data(), &bg[fbase], &bg[bbase], aH);
    probe.CollisionDiagonal(diagH.data(), diagl.data(), diagp.data());
    double rate = 0.;
    for (double d : diagH)
      rate = std::fmax(rate, std::fabs(d));
    for (double d : diagl)
      rate = std::fmax(rate, std::fabs(d));
    for (double d : diagp)
      rate = std::fmax(rate, std::fabs(d));
    // ONE condition, on the EARLY side only. The window runs from the first active row
    // to the end of the table.
    //
    // There WAS a late cut, on the parent's remaining comoving number (n_H > 1e-6 n_ini),
    // on the reasoning that "the parent has to exist for its collision to move anything".
    // That reasoning is wrong, and it is the bug this comment replaces. The operator does
    // not only move the parent: with inverse decays it couples the two daughter sectors to
    // each other, and that coupling is what keeps the reduced moment system dissipative.
    // Truncate it -- Apply() contributes exactly ZERO outside the tabulated range -- and
    // the moments are left with an undamped growing mode that runs to z = 0.
    //
    // Measured, Gamma = 1e6, m = 0.3, six moments, max|phi| (0.4734 is correct): the
    // sector tracks the matrix-free path to 1e-4 until a = 2.6e-2, then diverges, reaching
    // 2.1e27 by z = 0. The divergence begins ~3x PAST the old window's late edge, in the
    // region the table had truncated away. Removing the late cut alone restores
    // max|phi| = 0.4734 and C_2^TT = 1.4777e-10 at every background stepsize tested.
    //
    // Why it stayed hidden: where the last qualifying row lands depends on the background
    // table's row spacing, so the damage is a function of `back_integration_stepsize` --
    // a global precision knob this species does not own. At the 7e-3 default this code was
    // developed against, the cut fell late enough to be harmless; when #397 moved that
    // default to 0.07 every reduced run diverged, while the exact scheme (no table, no
    // window) and the matrix-free reduced path (no window) stayed correct.
    //
    // The cost of not truncating is bounded and small: the tail from the decay epoch to
    // a = 1 is ~2 decades, which at the default spacing is 65 further rows -- 21.6 MB and
    // 0.14 s against 14.6 MB and 0.10 s truncated (Gamma = 1e6, 28x28 blocks). The blocks
    // out there are near zero and cost only storage, which is the same trade the rows
    // inside the window already make.
    if (rate > kTableActiveRate * aH) {
      row_lo = row;
      break;
    }
  }
  // The late edge is the end of the table, unconditionally -- not "the last row that also
  // passed the test above". Deriving it from the rate would silently reinstate a late cut
  // the moment the ratio dipped back under the threshold anywhere, which is the exact
  // failure this function was rewritten to remove; the argument that it never dips (the
  // parent's diagonal tends to a*Gamma while aH falls, so rate/aH rises monotonically and
  // is still rising at a = 1) is a property of the physics, not something the code should
  // be relying on to stay dissipative.
  const int row_hi = n_rows - 1;
  if (row_lo < 0) {
    // Not gated: falling back to the matrix-free path is a silent 10x, and a user who
    // asked for the table needs to know they did not get it.
    fprintf(stderr,
            "[dncdm-reduced] %s: the collision is never active above %.0e*aH on any "
            "background row, so there is no window to tabulate; falling back to the "
            "matrix-free path\n",
            name().c_str(),
            kTableActiveRate);
    return;
  }

  // Pass 2: assemble every row in [row_lo, row_hi], one chunk per task. The kernel is
  // stateful (a cached transition geometry and a rolling multipole recurrence), so each
  // task builds its own; GridView is non-owning, so that duplicates only the per-node
  // scratch. Chunks are contiguous, so each task writes a slice the merge below just
  // concatenates; every row costs the same S assemblies, so there is no imbalance to
  // spread by interleaving.
  // The tabulation is interpolated in ln a, so its accuracy is set by the SPACING of the
  // points it holds -- and those were the background integrator's own rows, i.e.
  // `back_integration_stepsize`, a knob this table does not own and cannot see. At the
  // 7e-3 this code was developed against that spacing was fine; at the 0.07 of #397 it is
  // not, and the residue shows up where the operator varies fastest (Gamma = 1e9,
  // m = 0.06: max|phi| 101 against 0.4734 correct, the last cell still wrong once the
  // truncation above is fixed). So the table sub-divides each background interval until
  // its own spacing meets `dr_table_max_dlna`, interpolating the kernel's INPUTS -- the parent's
  // ln f, the daughters' f, a and H -- rather than the assembled operator. n_sub is global
  // and taken from the widest interval, so a background table uniform in ln a stays uniform
  // here and ReducedOperatorTable::LowerRow keeps its O(1) lookup.
  const int n_win = row_hi - row_lo + 1;
  double dlna_row = 0.;
  for (int row = row_lo; row < row_hi; ++row) {
    const double a0 = background_table[(size_t) row * row_stride + bgm_->index_bg_a_];
    const double a1 = background_table[(size_t) (row + 1) * row_stride + bgm_->index_bg_a_];
    dlna_row        = std::fmax(dlna_row, std::log(a1 / a0));
  }
  const int n_sub =
      std::max(1, static_cast<int>(std::ceil(dlna_row / table_max_dlna_ - kTableDlnaSlack)));
  const long long n_pts = (long long) (n_win - 1) * n_sub + 1;
  const int n_threads   = std::max(1,
                                   pba_ != nullptr ? static_cast<int>(pba_->number_of_threads) : 1);
  const int n_chunks    = static_cast<int>(std::min<long long>(n_pts, 4LL * n_threads));
  const size_t block    = (size_t) (l_max + 1) * S * S;
  std::vector<std::vector<double>> chunk_M(n_chunks);
  std::vector<std::vector<double>> chunk_a(n_chunks);

  const auto t0 = std::chrono::steady_clock::now();
  auto do_chunk = [&](int c) {
    const long long lo = (long long) c * n_pts / n_chunks;
    const long long hi = (long long) (c + 1) * n_pts / n_chunks;
    // Both final sizes are known exactly here, and chunk_M is the big one: at the heaviest
    // cell the table is 277 MB, so letting it grow by reallocation would recopy the
    // assembled blocks log2(n) times and transiently hold ~1.5x a chunk on top of what the
    // finished chunks already occupy -- while every chunk is live, since they are merged
    // only after the last one returns.
    chunk_a[c].reserve((size_t) (hi - lo));
    chunk_M[c].reserve((size_t) (hi - lo) * block);
    DecayTransitionKernel k2(parent_view,
                             fermion_view,
                             boson_view,
                             fermion_->statistics(),
                             boson_->statistics(),
                             cfg_,
                             fermion_bg_view,
                             boson_bg_view);
    ReducedCollisionOperator op(k2,
                                parent_view,
                                fermion_view,
                                boson_view,
                                fermion_->statistics(),
                                boson_->statistics(),
                                ReducedCollisionOperator::Config{reduced_moments_});
    std::vector<double> f(Np), f_l(Nf), f_phi(Nb), M;
    for (long long p = lo; p < hi; ++p) {
      const int row    = row_lo + static_cast<int>(p / n_sub);
      const int sub    = static_cast<int>(p % n_sub);
      const double t   = static_cast<double>(sub) / n_sub;
      const double* bg = background_table + (size_t) row * row_stride;
      // sub == 0 reads the row itself, so the last row -- which has no successor -- is only
      // ever reached there and bg1 is never formed past the end of the table.
      const double* bg1 = (sub == 0) ? bg : background_table + (size_t) (row + 1) * row_stride;
      double a, aH;
      const double* f_l_src;
      const double* f_phi_src;
      if (sub == 0) {
        a  = bg[bgm_->index_bg_a_];
        aH = a * bg[bgm_->index_bg_H_];
        for (int i = 0; i < Np; ++i)
          f[i] = kappa * std::exp(bg[pf + i]);
        f_l_src   = &bg[fbase];
        f_phi_src = &bg[bbase];
      }
      else {
        const double lna0 = std::log(bg[bgm_->index_bg_a_]);
        const double lnH0 = std::log(bg[bgm_->index_bg_H_]);
        a                 = std::exp(lna0 + t * (std::log(bg1[bgm_->index_bg_a_]) - lna0));
        aH                = a * std::exp(lnH0 + t * (std::log(bg1[bgm_->index_bg_H_]) - lnH0));
        // The parent is stored as ln f, so linear here is GEOMETRIC in f, which is what an
        // occupation that falls many decades across the window needs. The daughters are
        // stored as f and are built up from zero, where geometric would be wrong.
        for (int i = 0; i < Np; ++i)
          f[i] = kappa * std::exp(bg[pf + i] + t * (bg1[pf + i] - bg[pf + i]));
        for (int j = 0; j < Nf; ++j)
          f_l[j] = bg[fbase + j] + t * (bg1[fbase + j] - bg[fbase + j]);
        for (int j = 0; j < Nb; ++j)
          f_phi[j] = bg[bbase + j] + t * (bg1[bbase + j] - bg[bbase + j]);
        f_l_src   = f_l.data();
        f_phi_src = f_phi.data();
      }
      k2.PrepareTransitions(a, m, Gamma, f.data(), f_l_src, f_phi_src, aH);

      // The BASIS is the frozen one -- SetBackground here would rebuild it per row and
      // silently make every row's moments mean something different. Only the operator
      // varies with a; the coordinates must not.
      op.AdoptFrozenBasis(a, m, f.data(), *reduced_op_);
      op.Assemble(l_max, M);

      chunk_a[c].push_back(a);
      chunk_M[c].insert(chunk_M[c].end(), M.begin(), M.end());
    }
  };

  if (n_chunks > 1 && n_threads > 1) {
    Tools::TaskSystem ts(n_threads);
    std::vector<std::future<void>> fut;
    fut.reserve(n_chunks);
    for (int c = 0; c < n_chunks; ++c)
      fut.push_back(ts.AsyncTask([&, c]() { do_chunk(c); }));
    for (auto& f : fut)
      f.get();
  }
  else {
    for (int c = 0; c < n_chunks; ++c)
      do_chunk(c);
  }

  reduced_table_.Reset(S, l_max);
  for (int c = 0; c < n_chunks; ++c)
    for (size_t r = 0; r < chunk_a[c].size(); ++r)
      reduced_table_.AddRow(chunk_a[c][r],
                            std::vector<double>(chunk_M[c].begin() + r * block,
                                                chunk_M[c].begin() + (r + 1) * block));

  // The table is built and installed above; all that is left is to report it.
  if (background_verbose_ <= 1)
    return;
  const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  fprintf(stderr,
          "[dncdm-reduced] %s: tabulated M~_l at %d points across %d of %d background rows "
          "(dln a %.2e, %dx sub-divided; a = %.4e .. %.4e), %d x %d per multipole, %.1f MB, "
          "precompute %.2f s on %d threads. The k-loop no longer sweeps the %d-point "
          "daughter grid.\n",
          name().c_str(),
          reduced_table_.n_rows(),
          n_win,
          n_rows,
          dlna_row / n_sub,
          n_sub,
          background_table[(size_t) row_lo * row_stride + bgm_->index_bg_a_],
          background_table[(size_t) row_hi * row_stride + bgm_->index_bg_a_],
          S,
          S,
          reduced_table_.bytes() / 1048576.,
          secs,
          n_threads,
          Nf);
}

void DNCDMInvSpecies::ProcessBackgroundTable(const double* background_table,
                                             int n_rows,
                                             int row_stride,
                                             const double* /*z_table*/) {
  if (reduced_moments_ == 0 || bgm_ == nullptr)
    return;
  // ORDER MATTERS: the table is the reduction of the operator ONTO the frozen basis,
  // so the basis has to exist first.
  BuildReducedBasis(background_table, n_rows, row_stride);
  // reduced_table_l_max_ comes from the precision structure at Create time, so it is
  // still -1 for a composite built directly by a test (no ppr). Tabulating then would
  // walk the whole two-pass window scan to assemble zero-sized blocks; the RHS falls
  // through to the matrix-free path on an empty table anyway.
  if (reduced_table_l_max_ >= 0)
    BuildReducedTable(background_table, n_rows, row_stride);
}

// ─────────────────────────────────────────────────────────────────────────────
// Perturbations: per-ℓ collision coupling (design §4.3)
// ─────────────────────────────────────────────────────────────────────────────

std::unique_ptr<BaseSpecies::PerturbScratch> DNCDMInvSpecies::CreatePerturbScratch() const {
  auto s = std::make_unique<Scratch>();
  // A private kernel on the same grids and Config as the background one. GridView is
  // non-owning, so this duplicates the per-node scratch only, never the grids.
  const GridView parent_view{parent_->GetQ().data(), parent_->dq().data(), parent_->q_size()};
  const GridView fermion_pt_view{fermion_->q().data(), fermion_->dq().data(), fermion_->q_size()};
  const GridView boson_pt_view{boson_->q().data(), boson_->dq().data(), boson_->q_size()};
  const GridView fermion_bg_view{fermion_->q_bg().data(),
                                 fermion_->dq_bg().data(),
                                 static_cast<int>(fermion_->q_bg().size())};
  const GridView boson_bg_view{boson_->q_bg().data(),
                               boson_->dq_bg().data(),
                               static_cast<int>(boson_->q_bg().size())};
  s->kernel_pt = std::make_unique<DecayTransitionKernel>(parent_view,
                                                         fermion_pt_view,
                                                         boson_pt_view,
                                                         fermion_->statistics(),
                                                         boson_->statistics(),
                                                         cfg_,
                                                         fermion_bg_view,
                                                         boson_bg_view);
  s->fH_gather.assign(parent_->q_size(), 0.);
  s->F_H.assign(parent_->q_size(), 0.);
  s->dF_H.assign(parent_->q_size(), 0.);
  s->F_l.assign(fermion_->q_size(), 0.);
  s->dF_l.assign(fermion_->q_size(), 0.);
  s->F_phi.assign(boson_->q_size(), 0.);
  s->dF_phi.assign(boson_->q_size(), 0.);
  // Sized on the PERTURBATION grids, matching what kernel_pt writes: CollisionDiagonal
  // loops parent_.n / fermion_.n / boson_.n, which for this kernel are the pt views.
  s->diag_H_pt.assign(parent_->q_size(), 0.);
  s->diag_l_pt.assign(fermion_->q_size(), 0.);
  s->diag_phi_pt.assign(boson_->q_size(), 0.);
  return s;
}

void DNCDMInvSpecies::ApplyKernelPerturbDerivs(double a,
                                               int l_max_col,
                                               const double* pvecback,
                                               const NCDMBaseSpecies::PerturbLayout& p_lay,
                                               const NCDMBaseSpecies::PerturbLayout& f_lay,
                                               const NCDMBaseSpecies::PerturbLayout& b_lay,
                                               const double* y,
                                               double* dy,
                                               Scratch& scratch,
                                               double k) const {
  const double m              = parent_->GetMass();
  const double Gamma          = parent_->Gamma();
  const int Np                = parent_->q_size();
  DecayTransitionKernel& kern = *scratch.kernel_pt;
  // PERTURBATION-grid sizes for the daughters (== the layouts' q_size); the parent's
  // grid is untouched by dr_bg_refine.
  const int Nf = fermion_->q_size();
  const int Nb = boson_->q_size();

  // THE TABULATED PATH. M~_l is k-independent, so it was assembled once per background
  // row and is replayed here as a dense mat-vec. This skips PrepareTransitions, the
  // reconstruction onto the quadrature grid, the kernel sweep and the projection back --
  // the table already IS P M R.
  //
  // The guard is on the table's l_max as well as its existence: it is resolved from ppr
  // at Create time, before any layout exists, so if a layout ever asks for more
  // multipoles than were tabulated this falls through to the matrix-free path below
  // rather than reading past the block. That costs speed, never correctness.
  if (tabulated()) {
    const int L = std::min({p_lay.l_max, f_lay.l_max, b_lay.l_max, l_max_col});
    if (L <= reduced_table_.l_max()) {
      const int S        = Np + 2 * reduced_moments_;
      const double kappa = KappaStoredToBare();
      // Sized lazily, like the other Scratch buffers: S depends on the layout, which
      // does not exist at CreatePerturbScratch time.
      if (static_cast<int>(scratch.red_z.size()) < S) {
        scratch.red_z.resize(S);
        scratch.red_dz.resize(S);
      }
      for (int l = 0; l <= L; ++l) {
        // The parent leg crosses the kappa unit boundary (#385); the daughters' moments
        // are bare on both sides. kappa cancels on the parent's own diagonal, but not on
        // the off-diagonal blocks, so it has to be applied here and undone below.
        for (int i = 0; i < Np; ++i)
          scratch.red_z[i] = kappa * y[p_lay.index_per_q[i] + l];
        for (int j = 0; j < reduced_moments_; ++j) {
          scratch.red_z[Np + j]                    = y[f_lay.index_per_q[j] + l];
          scratch.red_z[Np + reduced_moments_ + j] = y[b_lay.index_per_q[j] + l];
        }
        std::fill(scratch.red_dz.begin(), scratch.red_dz.begin() + S, 0.);
        reduced_table_.Apply(a, l, scratch.red_z.data(), scratch.red_dz.data());
        for (int i = 0; i < Np; ++i)
          dy[p_lay.index_per_q[i] + l] += scratch.red_dz[i] / kappa;
        for (int j = 0; j < reduced_moments_; ++j) {
          dy[f_lay.index_per_q[j] + l] += scratch.red_dz[Np + j];
          dy[b_lay.index_per_q[j] + l] += scratch.red_dz[Np + reduced_moments_ + j];
        }
      }
      return;
    }
  }

  // Gather the background PSDs from the live background columns: parent publishes
  // ln f̄_H (exp back, then kappa-converted to BARE occupation, #385), daughters
  // publish bare f̄ directly. These are the SAME columns the children's
  // PerturbDerivs read, so the collision coefficients and the background the
  // hierarchies stream on are the same numbers.
  const int pf       = parent_->bg_lnf_index();
  const double kappa = KappaStoredToBare();
  for (int i = 0; i < Np; ++i)
    scratch.fH_gather[i] = kappa * std::exp(pvecback[pf + i]);
  // Daughters: hand the kernel the FULL background columns. It was built
  // with the fine grids as its background-gather views, so it samples f_l, f_phi at
  // the node momenta on the grid the background was integrated on, while still
  // depositing the collision onto the coarse PT grid the hierarchy lives on. No
  // subsampling copy is needed (or wanted: subsampling is what threw away the
  // injection front — dr_psd_species.h explains why splining that cliff is not an
  // option either, and this reads it at full resolution instead).
  const int fbase = fermion_->bg_f_index();
  const int bbase = boson_->bg_f_index();

  // One transition build per RHS; the operator reuses the node scratch for every ℓ.
  // Same cap as the background leg (see BackgroundDerivs). bgm_ is null when the
  // operator is driven directly by dncdm_inv_test with a stub pvecback; 0 then
  // disables the cap, which is what those exact-identity assertions require.
  const double a_prime_over_a = (bgm_ != nullptr) ? a * pvecback[bgm_->index_bg_H_] : 0.;
  kern.PrepareTransitions(a,
                          m,
                          Gamma,
                          scratch.fH_gather.data(),
                          &pvecback[fbase],
                          &pvecback[bbase],
                          a_prime_over_a);

  // ONE representation for all three species (#386): every leg carries the
  // UNNORMALIZED F = f̄Ψ = δf, which is what the kernel operator speaks natively.
  //
  // The parent used to be the exception — it evolved the normalized Ψ_H, so this
  // loop converted in (F_H = f̄_H Ψ_H) and out (/f̄_H) and carried a
  // −(f̄̇_H/f̄_H)Ψ_H dilution term to compensate. That is only well-conditioned
  // while f̄_H stays finite, and it does not: the parent's high-q bins decay past
  // 1e-37, at which point f̄̇_H/f̄_H — which with inverse decays on is
  // (K/ε)(f̄_l f̄_φ/f̄_H − 1), not a decay rate — reached 8.8e+5 / Mpc, and the two
  // terms of the write-back agreed to six digits so their difference was pure
  // roundoff. The solver was integrating noise: h collapsed to 8e-8 Mpc at
  // Γ=1e5 while Ψ_H ran off to 3e+5. In F-space none of that exists. The damping
  // coefficient is the physical −K/ε ≤ aΓ, the inverse-decay source
  // (K/ε)(f̄_φ F_l + f̄_l F_φ) is bounded, the dilution term cancels structurally
  // (the daughters' derivation in dr_psd_species.cpp, applied to the parent), and
  // F = 0 says exactly what it should: no perturbation to a distribution that is
  // no longer there.
  //
  // Beyond l_max_dncdm_col there is now nothing left to add — the dilution term
  // was the only all-ℓ contribution — so the loop simply stops there.
  const int L  = std::min({p_lay.l_max, f_lay.l_max, b_lay.l_max, l_max_col});
  const int L1 = L + 1;

  /* ONE kernel sweep over all multipoles instead of L+1 of them. The operator is 80% of
     total runtime and was re-streaming its whole node structure per multipole; inverting
     to node-outer/l-inner makes each node's data resident across the hierarchy. The
     gather/scatter below is the price: (Np+Nf+Nb)*L1 copies per RHS, negligible against
     the node work it saves. Results are unchanged -- the accumulation order per
     (bin, l) is still ascending i then n (see ApplyPerturbationOperatorAllL). */
  auto fit = [](std::vector<double>& v, size_t want) {
    if (v.size() < want)
      v.assign(want, 0.);
  };
  fit(scratch.FH_all, (size_t) Np * L1);
  fit(scratch.Fl_all, (size_t) Nf * L1);
  fit(scratch.Fphi_all, (size_t) Nb * L1);
  fit(scratch.dFH_all, (size_t) Np * L1);
  fit(scratch.dFl_all, (size_t) Nf * L1);
  fit(scratch.dFphi_all, (size_t) Nb * L1);

  // Same kappa boundary as before (#385): the parent's slot is STORED occupation, the
  // kernel speaks BARE. Daughters need no conversion in either direction.
  for (int i = 0; i < Np; ++i)
    for (int l = 0; l <= L; ++l)
      scratch.FH_all[(size_t) i * L1 + l] = kappa * y[p_lay.index_per_q[i] + l];
  // The daughters. In the REDUCED representation their y-slots hold moments, not grid
  // occupations, so they are reconstructed onto the quadrature grid here (R), swept by
  // the ordinary kernel below, and projected back (P) after it. That is the whole of
  // dF = P M R z, applied MATRIX-FREE: one kernel sweep, exactly as the exact scheme
  // pays. Assembling M~ per RHS would instead cost size() sweeps -- ~28x more than the
  // operator it replaces -- and precomputing it over the background table needs an
  // interpolation table in ln a. The saving here is therefore NOT in the RHS; it is that
  // the ODE carries 2*n_moments daughter variables per multipole instead of Nf + Nb.
  if (reduced()) {
    fit(scratch.red_m, (size_t) reduced_moments_);
    fit(scratch.red_grid, (size_t) std::max(Nf, Nb));
    for (int d = 0; d < 2; ++d) {
      const auto& lay  = (d == 0) ? f_lay : b_lay;
      const int n_grid = (d == 0) ? Nf : Nb;
      double* F        = (d == 0) ? scratch.Fl_all.data() : scratch.Fphi_all.data();
      for (int l = 0; l <= L; ++l) {
        for (int j = 0; j < reduced_moments_; ++j)
          scratch.red_m[j] = y[lay.index_per_q[j] + l];
        reduced_op_->Reconstruct(d, scratch.red_m.data(), scratch.red_grid.data());
        for (int i = 0; i < n_grid; ++i)
          F[(size_t) i * L1 + l] = scratch.red_grid[i];
      }
    }
  }
  else {
    for (int j = 0; j < Nf; ++j)
      for (int l = 0; l <= L; ++l)
        scratch.Fl_all[(size_t) j * L1 + l] = y[f_lay.index_per_q[j] + l];
    for (int k = 0; k < Nb; ++k)
      for (int l = 0; l <= L; ++l)
        scratch.Fphi_all[(size_t) k * L1 + l] = y[b_lay.index_per_q[k] + l];
  }

  kern.ApplyPerturbationOperatorAllL(L,
                                     L1,
                                     scratch.FH_all.data(),
                                     scratch.Fl_all.data(),
                                     scratch.Fphi_all.data(),
                                     scratch.dFH_all.data(),
                                     scratch.dFl_all.data(),
                                     scratch.dFphi_all.data());

  for (int l = 0; l <= L; ++l) {
    // Parent: the operator output verbatim, back to STORED units. No 1/f̄_H, no
    // floor guard, no dilution term (see the block comment above).
    for (int i = 0; i < Np; ++i)
      dy[p_lay.index_per_q[i] + l] += scratch.dFH_all[(size_t) i * L1 + l] / kappa;
    // Daughters: same F-space collision, already in their own units. An empty bin
    // MUST receive it — that is how it fills.
    if (reduced()) {
      // ...and the P of P M R: project the grid-resolved output back onto the moments
      // the state actually carries. Because P retains the number and energy functionals
      // exactly, whatever conservation the kernel's output has survives this step.
      for (int d = 0; d < 2; ++d) {
        const auto& lay  = (d == 0) ? f_lay : b_lay;
        const int n_grid = (d == 0) ? Nf : Nb;
        const double* dF = (d == 0) ? scratch.dFl_all.data() : scratch.dFphi_all.data();
        for (int i = 0; i < n_grid; ++i)
          scratch.red_grid[i] = dF[(size_t) i * L1 + l];
        reduced_op_->Project(d, scratch.red_grid.data(), scratch.red_m.data());
        for (int j = 0; j < reduced_moments_; ++j)
          dy[lay.index_per_q[j] + l] += scratch.red_m[j];
      }
    }
    else {
      for (int j = 0; j < Nf; ++j)
        dy[f_lay.index_per_q[j] + l] += scratch.dFl_all[(size_t) j * L1 + l];
      for (int k = 0; k < Nb; ++k)
        dy[b_lay.index_per_q[k] + l] += scratch.dFphi_all[(size_t) k * L1 + l];
    }
  }
}

void DNCDMInvSpecies::ApplyKernelPerturbDiagonal(double a,
                                                 int l_max_col,
                                                 const double* pvecback,
                                                 const NCDMBaseSpecies::PerturbLayout& p_lay,
                                                 const NCDMBaseSpecies::PerturbLayout& f_lay,
                                                 const NCDMBaseSpecies::PerturbLayout& b_lay,
                                                 double* diag,
                                                 Scratch& scratch) const {
  DecayTransitionKernel& kern = *scratch.kernel_pt;
  const double m              = parent_->GetMass();
  const double Gamma          = parent_->Gamma();
  const int Np = parent_->q_size(), Nf = fermion_->q_size(), Nb = boson_->q_size();

  // Identical gather to ApplyKernelPerturbDerivs: same columns, same kappa, so the
  // diagonal linearises the same discrete network the RHS evaluates.
  const int pf       = parent_->bg_lnf_index();
  const double kappa = KappaStoredToBare();
  for (int i = 0; i < Np; ++i)
    scratch.fH_gather[i] = kappa * std::exp(pvecback[pf + i]);
  const int fbase             = fermion_->bg_f_index();
  const int bbase             = boson_->bg_f_index();
  const double a_prime_over_a = (bgm_ != nullptr) ? a * pvecback[bgm_->index_bg_H_] : 0.;

  const int L = std::min({p_lay.l_max, f_lay.l_max, b_lay.l_max, l_max_col});

  // TABULATED: read the diagonal straight off the interpolated matrix. Doing this from
  // the kernel instead would force a PrepareTransitions per step, which is most of the
  // work the tabulation exists to remove -- and this diagonal is the better one anyway:
  // it is the actual diagonal of the operator the ODE sees, where the matrix-free path
  // can only offer each moment a Rayleigh quotient of the per-bin diagonal.
  //
  // The kappa unit boundary cancels on a DIAGONAL entry (dy = dF/kappa and F = kappa*y),
  // so the parent's entries need no conversion here, unlike in the mat-vec above.
  if (tabulated() && L <= reduced_table_.l_max()) {
    const int S = Np + 2 * reduced_moments_;
    if (static_cast<int>(scratch.red_dz.size()) < S)
      scratch.red_dz.resize(S);
    for (int l = 0; l <= L; ++l) {
      reduced_table_.Diagonal(a, l, scratch.red_dz.data());
      for (int i = 0; i < Np; ++i)
        diag[p_lay.index_per_q[i] + l] += scratch.red_dz[i];
      for (int j = 0; j < reduced_moments_; ++j) {
        diag[f_lay.index_per_q[j] + l] += scratch.red_dz[Np + j];
        diag[b_lay.index_per_q[j] + l] += scratch.red_dz[Np + reduced_moments_ + j];
      }
    }
    return;
  }

  kern.PrepareTransitions(a,
                          m,
                          Gamma,
                          scratch.fH_gather.data(),
                          &pvecback[fbase],
                          &pvecback[bbase],
                          a_prime_over_a);
  kern.CollisionDiagonal(scratch.diag_H_pt.data(),
                         scratch.diag_l_pt.data(),
                         scratch.diag_phi_pt.data());

  // d(dy_i)/d(y_i) is INVARIANT under the parent's stored/bare rescaling -- the
  // write-back divides the rate by the same kappa the read multiplies the
  // occupation by, so the two cancel. Same argument as BackgroundDerivsDiagonal.
  // The daughters store bare occupation and need no conversion either way.

  // In the reduced representation the daughters have no per-bin diagonal to hand out --
  // their state is moments, and the collision does not act diagonally on them. What the
  // exponential evolver wants is a representative rate per variable, so each moment gets
  // the RAYLEIGH QUOTIENT of the per-bin diagonal in that moment's own measure:
  //     diag_j = Sum dmu psi_j^2 diag(q) / Sum dmu psi_j^2,   dmu = dq q^2 w.
  // This is an integrating-factor choice, not a physical claim: ETD stays convergent for
  // any diagonal, and a wrong one costs step size rather than correctness. Using the
  // per-bin values would be simply ill-typed, and using zero would hand the daughters'
  // loss rate -- which is the stiff part after the parent goes extinct -- to the explicit
  // part of the splitting.
  double red_diag[2][DrPsdSpecies::kMaxReducedMoments] = {{0.}, {0.}};
  if (reduced()) {
    for (int d = 0; d < 2; ++d) {
      const DrPsdSpecies* sp = (d == 0) ? fermion_ : boson_;
      const double* dg       = (d == 0) ? scratch.diag_l_pt.data() : scratch.diag_phi_pt.data();
      const int n_grid       = sp->q_size();
      const std::vector<double>& psi = reduced_op_->psi(d);
      const std::vector<double>& w   = reduced_op_->weight(d);
      for (int j = 0; j < reduced_moments_; ++j) {
        double num = 0., den = 0.;
        for (int i = 0; i < n_grid; ++i) {
          const double p2   = psi[(size_t) j * n_grid + i] * psi[(size_t) j * n_grid + i];
          const double dmu  = sp->dq()[i] * sp->q()[i] * sp->q()[i] * w[i] * p2;
          num              += dmu * dg[i];
          den              += dmu;
        }
        red_diag[d][j] = (den > 0.) ? num / den : 0.;
      }
    }
  }

  for (int l = 0; l <= L; ++l) {
    for (int i = 0; i < Np; ++i)
      diag[p_lay.index_per_q[i] + l] += scratch.diag_H_pt[i];
    if (reduced()) {
      for (int j = 0; j < reduced_moments_; ++j) {
        diag[f_lay.index_per_q[j] + l] += red_diag[0][j];
        diag[b_lay.index_per_q[j] + l] += red_diag[1][j];
      }
    }
    else {
      for (int j = 0; j < Nf; ++j)
        diag[f_lay.index_per_q[j] + l] += scratch.diag_l_pt[j];
      for (int k = 0; k < Nb; ++k)
        diag[b_lay.index_per_q[k] + l] += scratch.diag_phi_pt[k];
    }
  }
}

void DNCDMInvSpecies::AddCouplingDerivs(double /*tau*/,
                                        const double* y,
                                        double* dy,
                                        const perturb_parameters_and_workspace& ppaw) const {
  const perturb_workspace* ppw    = ppaw.ppw;
  const precision* ppr            = ppaw.perturbations_module->GetPrecision();
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  const BaseSpecies::PerturbLayout& my        = *ppw->pv->species_layouts[collection_index_];
  const NCDMBaseSpecies::PerturbLayout& p_lay = parent_layout(my);
  const NCDMBaseSpecies::PerturbLayout& f_lay = fermion_layout(my);
  const NCDMBaseSpecies::PerturbLayout& b_lay = boson_layout(my);

  // l_max_dncdm_col truncates the collision. Clamped to the layout's own l_max here
  // rather than validated at input: the hierarchy length is what the indexing can
  // actually reach, and a collision multipole above it has nothing to write to.
  int l_max_col = ppr->l_max_dncdm_col;
  if (l_max_col > p_lay.l_max)
    l_max_col = p_lay.l_max;

  // Per-k working memory (see Scratch): one private kernel + gather buffers per
  // k-mode, so concurrent modes cannot corrupt each other's transition geometry
  // or Legendre recurrence.
  auto* scratch = static_cast<Scratch*>(ppw->species_scratch[collection_index_].get());
  ApplyKernelPerturbDerivs(ctx.a,
                           l_max_col,
                           ppw->pvecback.data(),
                           p_lay,
                           f_lay,
                           b_lay,
                           y,
                           dy,
                           *scratch,
                           ctx.k);
}

void DNCDMInvSpecies::PerturbDerivsDiagonal(const BaseSpecies::PerturbLayout& /*layout*/,
                                            double /*tau*/,
                                            const double* /*y*/,
                                            double* diag,
                                            const perturb_parameters_and_workspace& ppaw) const {
  // Children report no diagonal of their own (free-streaming is pure l-coupling,
  // which has none), so unlike the background leg there is no
  // CompositeSpecies::PerturbDerivsDiagonal call to chain here -- the coupling IS
  // the whole diagonal. Adding the base call anyway would be a no-op, but it would
  // also imply the children have something to say, and they do not.
  const perturb_workspace* ppw = ppaw.ppw;
  const precision* ppr         = ppaw.perturbations_module->GetPrecision();

  const BaseSpecies::PerturbLayout& my        = *ppw->pv->species_layouts[collection_index_];
  const NCDMBaseSpecies::PerturbLayout& p_lay = parent_layout(my);
  const NCDMBaseSpecies::PerturbLayout& f_lay = fermion_layout(my);
  const NCDMBaseSpecies::PerturbLayout& b_lay = boson_layout(my);

  int l_max_col = ppr->l_max_dncdm_col;
  if (l_max_col > p_lay.l_max)
    l_max_col = p_lay.l_max;

  // `a` comes from the pvecback the DIAGONAL callback just interpolated, NOT from
  // ppw->scalar_ctx: scalar_ctx.a is written by perturb_derivs_member, and on a
  // REJECTED step the last RHS ran at an ETDRK4 stage point ahead of t, so ctx.a
  // would be a(t+h) paired with pvecback at t. That mismatch is the failure mode the
  // module callback's comment warns about -- it stays stable and quietly degrades the
  // step. Reading both from the same vector makes it impossible by construction.
  const double* pvecback = ppw->pvecback.data();
  const double a         = pvecback[bgm_->index_bg_a_];

  auto* scratch = static_cast<Scratch*>(ppw->species_scratch[collection_index_].get());
  ApplyKernelPerturbDiagonal(a, l_max_col, pvecback, p_lay, f_lay, b_lay, diag, *scratch);
}

// ─────────────────────────────────────────────────────────────────────────────
// Background column output (parent + both daughters)
// ─────────────────────────────────────────────────────────────────────────────

void DNCDMInvSpecies::WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const {
  parent_->WriteBackgroundColumnTitles(w);
  fermion_->WriteBackgroundColumnTitles(w);
  boson_->WriteBackgroundColumnTitles(w);
  w.Add("bose_occupancy_" + name(), 0.);
}

void DNCDMInvSpecies::WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const {
  parent_->WriteBackgroundData(pvecback, w);
  fermion_->WriteBackgroundData(pvecback, w);
  boson_->WriteBackgroundData(pvecback, w);
  // max_q f̄_φ, the neutrino-lasing order parameter (thesis §6.4): the Bose runaway
  // is driven by 1+f̄_φ ≈ f̄_φ, so f̄_φ ≳ 1 is the condition for the (1±f) terms to
  // dominate. The per-bin f_dr_*_phi[i] columns already carry the full PSD; this is
  // the single scalar that makes the onset visible without post-processing 100
  // columns, and its dependence on dr_q_min is the whole q_min-convergence question.
  const double* f_phi = pvecback + boson_->bg_f_index();
  const int N         = static_cast<int>(boson_->q_bg().size());
  double f_max        = f_phi[0];
  for (int i = 1; i < N; ++i)
    f_max = std::max(f_max, f_phi[i]);
  w.Add("bose_occupancy_" + name(), f_max);
}

// ─────────────────────────────────────────────────────────────────────────────
// Closure + shooter hooks (combined parent+fermion+boson reserve, mirrors DNCDM_DR)
// ─────────────────────────────────────────────────────────────────────────────

double DNCDMInvSpecies::GetOmega0() const {
  if (parent_->Omega_dncdmdr_pending().has_value())
    return *parent_->Omega_dncdmdr_pending();
  return CompositeSpecies::GetOmega0();  // pre-shooting discovery fallback
}

std::vector<ShootingTarget> DNCDMInvSpecies::GetShootingTargets() const {
  if (parent_->InitialAbundanceMode()) {
    const double seed = parent_->Omega_dncdmdr_pending().value_or(
        parent_->Omega_ini_pending().value_or(0.1));
    return {{name() + ".Omega_dncdmdr_fixedpoint", name() + ".Omega_dncdmdr", seed}};
  }
  if (parent_->Omega_dncdmdr_pending().has_value()) {
    return {{name() + ".Omega_dncdmdr", name() + ".deg", *parent_->Omega_dncdmdr_pending()}};
  }
  return {};
}

void DNCDMInvSpecies::ComputeShootingGuess(const SpeciesBuildContext& ctx,
                                           std::vector<double>& guess,
                                           std::vector<double>& dxdy) const {
  if (parent_->InitialAbundanceMode()) {
    double seed = 0.1;
    if (parent_->Omega_ini_pending().has_value())
      seed = *parent_->Omega_ini_pending();
    else if (parent_->Neff_ini_pending().has_value() && ctx.pba)
      seed = *parent_->Neff_ini_pending() * 7. / 8. * std::pow(4. / 11., 4. / 3.) *
             ctx.pba->Omega0_g;
    guess.push_back(seed);
    dxdy.push_back(1.0);
    return;
  }
  if (!parent_->Omega_dncdmdr_pending().has_value())
    return;
  auto [g, d] = parent_->DegGuessFromOmegaToday(ctx, *parent_->Omega_dncdmdr_pending());
  guess.push_back(g);
  dxdy.push_back(d);
}

double DNCDMInvSpecies::ComputeShootingResidual(const ShootingResidualContext& ctx,
                                                const ShootingTarget& target) const {
  const double* bg      = ctx.bg_today;
  const double H0       = ctx.pba->H0;
  const double combined = (parent_->Rho(bg) + fermion_->Rho(bg) + boson_->Rho(bg)) / (H0 * H0);
  if (target.target_name == name() + ".Omega_dncdmdr_fixedpoint")
    return -combined + GetOmega0();
  return combined - target.target_value;
}

void DNCDMInvSpecies::PrintVariables(PerturbColumnWriter& w,
                                     const BaseSpecies::PerturbLayout* base,
                                     double tau,
                                     const double* y,
                                     const PerturbationsModule& mod,
                                     const perturb_workspace* ppw) const {
  // Chain first (children currently publish nothing) so a child that later gains
  // PrintVariables is picked up without editing this; then the sector aggregates
  // below, which no child can form on its own.
  CompositeSpecies::PrintVariables(w, base, tau, y, mod, ppw);
  // a⁴Π_νφ ≡ a⁴·Σ_{i∈{H,l,φ}} (ρ̄_i+p̄_i)σ_i — the total anisotropic stress of the
  // decaying sector, eq. (6.19) of the thesis / the Ma&Bertschinger source of
  // k²(φ−ψ) = 12πGa²Π. Each child already forms (ρ+p)σ in its own StressEnergy
  // (the parent from its normalized Ψ moments, the daughters from their
  // unnormalized F moments), so this only sums the three and multiplies by a⁴.
  // Densities are in CLASS's 8πG/3 units and the thesis absorbs the overall
  // normalization into its constant C, so only the TIME DEPENDENCE is physical:
  // a⁴ removes the free radiation redshift, leaving a constant-amplitude
  // free-streaming oscillation against which interaction-induced isotropisation
  // shows up as a decaying envelope.
  //
  // Composite-owned (mirrors Type3Species::PrintVariables): the children are not
  // separately registered in the collection, so their layouts are reached through
  // this composite's nested layout, never through their own collection_index_.
  double pi_H = 0., pi_l = 0., pi_phi = 0.;
  // The sector's OTHER two moments, on the same footing as a4Pi: δ is ρ-weighted,
  // θ and σ are (ρ+p)-weighted, all three summed over {H, l, φ}. Assembled here
  // rather than downstream from three published ratios because the sum has to be
  // taken on the UNNORMALIZED moments: once the parent is extinct δ_H = δρ_H/ρ_H is
  // a ratio of two vanishing numbers (measured: ρ_H/ρ_sector = 1e-23 at a = 0.1 for
  // Γ=1e5) and any weighting applied to the published ratio inherits that noise,
  // while Σδρ / Σρ never forms it. The per-component columns are published too, so
  // the decomposition stays checkable — they are the ones to distrust late.
  double d_sec = 0., t_sec = 0., s_sec = 0.;
  double d_H = 0., d_l = 0., d_phi = 0., t_H = 0., t_l = 0., t_phi = 0.;
  if (!w.IsTitleMode()) {
    const perturb_vector* pv = ppw->pv.get();
    const double* pvecback   = ppw->pvecback.data();
    const double a           = pvecback[mod.GetBackgroundModule()->index_bg_a_];
    const double a4          = a * a * a * a;
    const auto& my           = static_cast<const CompositeSpecies::PerturbLayout&>(*base);
    const auto se_H   = parent_->StressEnergy(*my.child_layouts[kParent], pv, y, pvecback, ppw);
    const auto se_l   = fermion_->StressEnergy(*my.child_layouts[kFermion], pv, y, pvecback, ppw);
    const auto se_phi = boson_->StressEnergy(*my.child_layouts[kBoson], pv, y, pvecback, ppw);
    pi_H              = a4 * se_H.rho_plus_p_shear;
    pi_l              = a4 * se_l.rho_plus_p_shear;
    pi_phi            = a4 * se_phi.rho_plus_p_shear;

    StressEnergyContribution tot  = se_H;
    tot                          += se_l;
    tot                          += se_phi;
    auto ratio       = [](double num, double den) { return (den > 0.) ? num / den : 0.; };
    const double rpp = tot.rho + tot.p;
    d_sec            = ratio(tot.delta_rho, tot.rho);
    t_sec            = ratio(tot.rho_plus_p_theta, rpp);
    s_sec            = ratio(tot.rho_plus_p_shear, rpp);
    d_H              = ratio(se_H.delta_rho, se_H.rho);
    d_l              = ratio(se_l.delta_rho, se_l.rho);
    d_phi            = ratio(se_phi.delta_rho, se_phi.rho);
    t_H              = ratio(se_H.rho_plus_p_theta, se_H.rho + se_H.p);
    t_l              = ratio(se_l.rho_plus_p_theta, se_l.rho + se_l.p);
    t_phi            = ratio(se_phi.rho_plus_p_theta, se_phi.rho + se_phi.p);

    // Synchronous → Newtonian, exactly as NCDMSpecies::PrintVariables does it for a
    // stand-alone ncdm. Without this the decaying sector's δ and θ would not be
    // comparable with the `delta_<name>` / `theta_<name>` columns of the
    // self-interacting comparison runs, which are gauge-corrected. σ is
    // gauge-invariant and is left alone. The correction is linear in each species'
    // (1+w_i), so applying it to the ρ-weighted sum with the sector's own w is the
    // same as correcting each component first and then summing.
    if (mod.GetPerturbs()->gauge == possible_gauges::synchronous) {
      const double H     = pvecback[mod.GetBackgroundModule()->index_bg_H_];
      const double k     = ppw->scalar_ctx.k;
      const double alpha = ppw->pvecmetric[ppw->index_mt_alpha];
      auto to_newt       = [&](double& d, double& t, double rho, double p) {
        d -= 3. * a * H * (1. + ((rho > 0.) ? p / rho : 0.)) * alpha;
        t += k * k * alpha;
      };
      to_newt(d_sec, t_sec, tot.rho, tot.p);
      to_newt(d_H, t_H, se_H.rho, se_H.p);
      to_newt(d_l, t_l, se_l.rho, se_l.p);
      to_newt(d_phi, t_phi, se_phi.rho, se_phi.p);
    }
  }
  // Total first (the measured quantity), then the decomposition: the thesis' own
  // footnote flags that the a⁴Π_νφ notation hides a parent contribution which is
  // "often negligible" — with the three columns side by side that is checkable
  // rather than assumed.
  w.Add("a4Pi_" + name(), pi_H + pi_l + pi_phi, true);
  w.Add("a4Pi_" + name() + "_H", pi_H, true);
  w.Add("a4Pi_" + name() + "_l", pi_l, true);
  w.Add("a4Pi_" + name() + "_phi", pi_phi, true);

  w.Add("delta_" + name(), d_sec, true);
  w.Add("theta_" + name(), t_sec, true);
  w.Add("shear_" + name(), s_sec, true);
  w.Add("delta_" + name() + "_H", d_H, true);
  w.Add("delta_" + name() + "_l", d_l, true);
  w.Add("delta_" + name() + "_phi", d_phi, true);
  w.Add("theta_" + name() + "_H", t_H, true);
  w.Add("theta_" + name() + "_l", t_l, true);
  w.Add("theta_" + name() + "_phi", t_phi, true);

  // ── Transport rate: the COLLISIONAL part of d(a⁴Π)/dτ, and the parent's γ ─────
  //
  // a⁴Π is simultaneously BUILT UP by free streaming and BROKEN DOWN by collisions,
  // so its total logarithmic derivative is not the transport rate -- in a no-decay
  // reference run that expression is nonzero and entirely Γ-independent. The
  // collision term is emitted directly instead: (ρ+p)σ is LINEAR AND HOMOGENEOUS in
  // the perturbation moments (each child's StressEnergy contracts y against fixed
  // background weights), so evaluating the SAME contraction on a vector holding only
  // the kernel's deposit gives d(a⁴Π)/dτ|_coll exactly, with no differencing and no
  // second copy of the contraction to drift out of sync with the a4Pi columns above.
  //
  // The isotropisation rate then follows from the emitted columns as
  //     Γ_T = -a4PiDot_coll / (a · a4Pi),
  // the 1/a converting the conformal-time rate to the convention in which
  // d(a⁴Π)/dτ = -aΓ_T·a⁴Π. ⚠ Γ_T is a small residual of two much larger terms, so
  // read it as an emitted quantity to be checked, not as a converged observable.
  // Sampling-point cost only; never touched on the RHS.
  double pidot_H = 0., pidot_l = 0., pidot_phi = 0., gamma_H = 0.;
  if (!w.IsTitleMode()) {
    const perturb_vector* pv = ppw->pv.get();
    const double* pvecback   = ppw->pvecback.data();
    const double a           = pvecback[mod.GetBackgroundModule()->index_bg_a_];
    const double a4          = a * a * a * a;
    const auto& my           = static_cast<const CompositeSpecies::PerturbLayout&>(*base);
    auto* scratch            = static_cast<Scratch*>(ppw->species_scratch[collection_index_].get());
    if (scratch != nullptr && scratch->kernel_pt != nullptr) {
      const auto& p_lay = parent_layout(my);
      int l_max_col     = mod.GetPrecision()->l_max_dncdm_col;
      if (l_max_col > p_lay.l_max)
        l_max_col = p_lay.l_max;
      // Zeroed every call: ApplyKernelPerturbDerivs ACCUMULATES into dy (it is
      // written to be one contributor among several on the real RHS), so a stale
      // buffer would silently sum sampling points together.
      scratch->dy_coll.assign(pv->pt_size, 0.);
      ApplyKernelPerturbDerivs(a,
                               l_max_col,
                               pvecback,
                               p_lay,
                               fermion_layout(my),
                               boson_layout(my),
                               y,
                               scratch->dy_coll.data(),
                               *scratch);
      const double* dyc = scratch->dy_coll.data();
      // The DAUGHTERS' StressEnergy contracts the y it is HANDED, so evaluating it on
      // the collision vector gives their Π-rate directly.
      pidot_l = a4 * fermion_->StressEnergy(*my.child_layouts[kFermion], pv, dyc, pvecback, ppw)
                         .rho_plus_p_shear;
      pidot_phi =
          a4 *
          boson_->StressEnergy(*my.child_layouts[kBoson], pv, dyc, pvecback, ppw).rho_plus_p_shear;
      // The PARENT'S DOES NOT: DNCDMSpecies::RescaledPerturbations reads
      // ppw->pv->y[index_pt] directly (dncdm_species.cpp), ignoring the y argument.
      // That is harmless on the RHS, where the two are the same vector, but it means
      // parent_->StressEnergy(..., dy_coll, ...) would silently return the ordinary
      // Π_H rather than its rate.
      //
      // Rescale instead. Π_H = norm · Σ_q dq q⁴/ε · y₂ with a y-INDEPENDENT norm (the
      // 2/3, (ρ+p), the exp(lnN) rescaling and the ρ+p normalisation all factor out),
      // so the same norm carries the collision vector and the ratio of the two raw
      // ℓ=2 contractions converts one into the other. The rescaling factor cancels
      // between numerator and denominator, so it need not be recomputed here, and
      // there is no second copy of the normalisation to drift out of sync.
      // Weights are the collision-owned ones; the composite always sets that mode on
      // its parent (SetCollisionOwned(true) in this constructor).
      double num = 0., den = 0.;
      for (int i = 0; i < parent_->q_size(); ++i) {
        const double q    = parent_->GetQ()[i];
        const double eps  = std::sqrt(q * q + a * a * parent_->GetMass() * parent_->GetMass());
        const double wq   = parent_->dq()[i] * q * q * q * q / eps;
        const int i2      = p_lay.index_per_q[i] + 2;
        num              += wq * dyc[i2];
        den              += wq * pv->y[i2];
      }
      // pi_H/den is the y-independent norm and stays well conditioned even at a zero
      // crossing of den, because pi_H carries the same vanishing factor.
      pidot_H = (den != 0.) ? (pi_H / den) * num : 0.;
    }
    // Parent's number-weighted mean Lorentz factor <eps>/(a·M). Emitted rather than
    // reconstructed downstream, which would mean re-deriving the parent's q grid and
    // its stored-vs-bare f convention outside the code that owns them.
    const double M = parent_->GetMass();
    if (M > 0. && a > 0.) {
      const auto& q  = parent_->GetQ();
      const auto& dq = parent_->dq();
      double num = 0., den = 0.;
      for (int i = 0; i < parent_->q_size(); ++i) {
        const double f   = std::exp(pvecback[parent_->bg_lnf_index() + i]);
        const double w0  = dq[i] * q[i] * q[i] * f;
        num             += w0 * std::sqrt(q[i] * q[i] + a * a * M * M);
        den             += w0;
      }
      gamma_H = (den > 0.) ? num / (den * a * M) : 0.;
    }
  }
  w.Add("a4PiDot_coll_" + name(), pidot_H + pidot_l + pidot_phi, true);
  w.Add("a4PiDot_coll_" + name() + "_H", pidot_H, true);
  w.Add("a4PiDot_coll_" + name() + "_l", pidot_l, true);
  w.Add("a4PiDot_coll_" + name() + "_phi", pidot_phi, true);
  w.Add("gamma_" + name(), gamma_H, true);

  // The parent's own background density and the UNNORMALIZED moments of its
  // perturbation. rho_H is what makes the sector's contribution to the metric
  // physically negligible once the parent has decayed; drho_H is the moment that
  // actually enters that metric source. Publishing them side by side is the only
  // way to see whether a large delta_H = drho_H/rho_H is a harmless ratio of two
  // vanishing numbers or is genuinely driving the metric — the two look identical
  // in any normalized output.
  double rho_H = 0., drho_H = 0., rpp_shear_H = 0.;
  if (!w.IsTitleMode()) {
    const perturb_vector* pv = ppw->pv.get();
    const double* pvecback   = ppw->pvecback.data();
    const auto& my           = static_cast<const CompositeSpecies::PerturbLayout&>(*base);
    const auto se = parent_->StressEnergy(*my.child_layouts[kParent], pv, y, pvecback, ppw);
    rho_H         = se.rho;
    drho_H        = se.delta_rho;
    rpp_shear_H   = se.rho_plus_p_shear;
  }
  w.Add("rho_" + name() + "_H", rho_H, true);
  w.Add("drho_" + name() + "_H", drho_H, true);
  w.Add("rppshear_" + name() + "_H", rpp_shear_H, true);
}
