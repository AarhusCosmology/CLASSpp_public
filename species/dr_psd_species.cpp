#include "dr_psd_species.h"

#include <algorithm>
#include <cmath>

#include "arrays.h"
#include "background_module.h"
#include "errors.h"
#include "perturbations_module.h"
#include "precision.h"
#include "species/species_input.h"

// ─────────────────────────────────────────────────────────────────────────────
// Constructor
// ─────────────────────────────────────────────────────────────────────────────

DrPsdSpecies::DrPsdSpecies(FileContent* pfc,
                           const std::string& instance_name,
                           const NcdmSettings& settings,
                           const background* pba,
                           const BackgroundModule* bgm,
                           Statistics stat,
                           double kinematic_q_max)
    : NCDMBaseSpecies(instance_name, EnergyType::Other, pfc, instance_name, settings, DeferInit{}),
      pba_(pba), stat_(stat) {
  bgm_ = bgm;
  SpeciesInput input(pfc, instance_name);

  // Massless daughter: mass zero, T_dr = T_cmb convention (dimensionless q = p/T).
  m_in_eV_ = 0.;
  M_       = 0.;
  T_       = input.get_or("dr_T", 1.0);   // T_dr / T_cmb
  ksi_     = input.get_or("dr_ksi", 0.);  // chemical potential mu/T

  // Spin dof: fermion g = 2, boson g = 1. The g_H·g_l/g_φ = 2 boson dof ratio is
  // folded into the kernel deposit (design §4.1), not here.
  // Bare-occupation storage (kernel boundary, #385) => factor carries g/(2*pi)^3,
  // so deg means physical dof (g=2 fermion = particle+antiparticle, matching a
  // deg=1 standard NCDM neutrino; g=1 boson).
  SetDegAndFactor((stat_ == Statistics::Fermion ? 2.0 : 1.0) / pow(2 * _PI_, 3));

  // Grid controls. dr_q_min/dr_q_max/dr_N_q define the BACKGROUND grid; the
  // perturbation grid is its exact dr_bg_refine-subsample (no independent endpoints —
  // sharing them is what keeps both grids covering the same decay band).
  q_min_bg_ = input.get_or("dr_q_min", q_min_bg_);
  N_q_bg_   = input.get_or("dr_N_q", N_q_bg_);
  refine_   = input.get_or("dr_bg_refine", refine_);

  // dr_q_max: an explicit value ALWAYS wins (the composite warns once, and only
  // warns, if it sits below the emission band -- an undersized grid is sometimes a
  // deliberate cost choice and the run still has to be reproducible). Absent the key,
  // a composite-built daughter is sized to its own kinematics rather than to the
  // standalone default: the band top is a*m/2 for a non-relativistic parent, which
  // for a 1 eV parent at T = 0.71611 is q ~ 2980 at a = 1, ten times the fixed
  // default. Under-sizing does not degrade gracefully -- TwoBinSplit clamps the
  // deposit into the top bin, keeping number and losing energy in proportion to
  // q_max/q* (measured: the daughters received 40% of the parent's energy loss at
  // z=3, 20% at z=1, and the sector's comoving energy came out 39% low at Gamma=1e3).
  //
  // Note the cost asymmetry, which is why this default is safe for log sampling and
  // blunt for linear: on a LOG grid a ten-fold q_max is one more decade of bins, so
  // widening is nearly free (measured flat, 124 s -> 121 s over q_max 300 -> 1e4). On
  // a LINEAR grid the same widening coarsens every cell by the same factor, since
  // dr_N_q is deliberately NOT adjusted here -- the bin count is the user's to choose
  // and silently inflating it would change the ODE dimension behind their back.
  if (auto explicit_q_max = input.get<double>("dr_q_max")) {
    q_max_bg_ = *explicit_q_max;
  }
  else if (kinematic_q_max > 0.) {
    q_max_bg_ = kQMaxMargin * kinematic_q_max;
  }

  class_test(!(q_min_bg_ > 0.) || !(q_max_bg_ > q_min_bg_) || N_q_bg_ < 2,
             "species '%s': require 0 < dr_q_min < dr_q_max and dr_N_q >= 2",
             instance_name.c_str());
  class_test(refine_ < 1,
             "species '%s': dr_bg_refine (=%d) must be >= 1 (1 = perturbation grid identical "
             "to the background grid)",
             instance_name.c_str(),
             refine_);

  // Momentum-grid family. "log" (default) is the repo's log-trapezoid; "linear" is
  // quadrature.c's qm_trapz, present so this species can be run on the SAME grid as
  // the reference implementations (arXiv:2011.01502's background figures use a
  // uniform 240-bin grid to q_max=22) and so the historical uniform-bin neutrino-
  // lasing calculation is reachable at all. Numeric, not structural: an unknown
  // value rejects the point rather than aborting the run.
  {
    const std::string sampling = input.get_or("dr_q_sampling", std::string("log"));
    class_test(sampling != "log" && sampling != "linear",
               "species '%s': dr_q_sampling (='%s') must be 'log' or 'linear'",
               instance_name.c_str(),
               sampling.c_str());
    q_sampling_ = (sampling == "linear") ? QSampling::Linear : QSampling::Log;
  }
  // Divisibility of the exact subsample (see bg_index): the log grid counts
  // INTERVALS from its endpoint at index 0, the linear grid counts CELLS from its
  // implicit zero, so the two differ by one. Numeric (not structural): it couples
  // two sampler-varyable integers, so a bad combination rejects the point rather
  // than aborting the run.
  if (q_sampling_ == QSampling::Linear) {
    class_test(N_q_bg_ % refine_ != 0,
               "species '%s': with dr_q_sampling=linear, dr_bg_refine (=%d) must divide "
               "dr_N_q (=%d) exactly. The perturbation grid is the subsample "
               "q_pt[j] = q_bg[(j+1)*dr_bg_refine - 1] = (j+1)*dr_bg_refine*h, i.e. N_bg "
               "cells of width h regrouped into N_bg/dr_bg_refine cells of width "
               "dr_bg_refine*h. Choose dr_N_q = dr_bg_refine*n",
               instance_name.c_str(),
               refine_,
               N_q_bg_);
  }
  else {
    class_test((N_q_bg_ - 1) % refine_ != 0,
               "species '%s': dr_bg_refine (=%d) must divide dr_N_q-1 (=%d) exactly. The "
               "perturbation grid is the exact subsample q_pt[j] = q_bg[j*dr_bg_refine], and N "
               "points span N-1 intervals, so the requirement is N_bg-1 = dr_bg_refine*(N_pt-1) "
               "(NOT N_bg = dr_bg_refine*N_pt). Choose dr_N_q = dr_bg_refine*n + 1",
               instance_name.c_str(),
               refine_,
               N_q_bg_ - 1);
  }

  // Initial occupation amplitude, PER DAUGHTER. This is the paper's A_n
  // (arXiv:2011.01502; `ncdm_psd_parameters` in the reference implementation):
  // f(q, a_ini) = A · f_eq(q), with f_eq the equilibrium PSD for this statistics.
  //
  // The daughters are constructed by the composite under ONE instance name, so
  // every other dr_* key is necessarily shared between them — which is correct:
  // the kernel requires parent and both daughters to report occupation on the SAME
  // dimensionless q axis at the same T (see the dr_T pinning guard in
  // DNCDMInvSpecies::Create), so per-daughter grids or temperatures would only let
  // a user build a physically invalid sector. The initial abundance is the one
  // quantity that legitimately differs per daughter, and it is exactly the freedom
  // the published background figures use: A = [0, 1, 1] means phi starts EMPTY
  // while nu_l starts at a full Fermi-Dirac. That configuration is unreachable if
  // both daughters share a single seeding flag, so it gets statistics-qualified
  // keys — resolved off this instance's own statistics, no composite plumbing.
  {
    const char* key = (stat_ == Statistics::Fermion) ? "dr_f_ini_l" : "dr_f_ini_phi";
    if (auto explicit_A = input.get<double>(key)) {
      f_ini_ = *explicit_A;
      class_test(!(f_ini_ >= 0.),
                 "species '%s': %s (=%g) must be >= 0 (it multiplies the equilibrium PSD)",
                 instance_name.c_str(),
                 key,
                 f_ini_);
    }
    else {
      // Back-compatible default: any thermal key or Omega_ini seeds BOTH daughters
      // at full equilibrium, which is what this species did before the per-daughter
      // keys existed. Absent all of them, daughters start empty.
      f_ini_ = (input.get<double>("Omega_ini").has_value() ||
                input.get<double>("dr_T").has_value() || input.get<double>("dr_ksi").has_value())
                   ? 1.
                   : 0.;
    }
  }

  if (q_sampling_ == QSampling::Linear)
    BuildLinearTrapezoidGrids();
  else
    BuildLogTrapezoidGrids();
}

void DrPsdSpecies::PinTemperature(double T) {
  T_ = T;
  SetDegAndFactor(GetDeg());
}

Statistics DrPsdSpecies::ParseStatistics(FileContent* pfc, const std::string& instance_name) {
  SpeciesInput input(pfc, instance_name);
  const std::string s = input.get_or<std::string>("statistics", "fermion");
  if (s == "fermion")
    return Statistics::Fermion;
  if (s == "boson")
    return Statistics::Boson;
  class_stop_severe("species '%s': statistics must be 'fermion' or 'boson' (got '%s')",
                    instance_name.c_str(),
                    s.c_str());
}

std::vector<Named> DrPsdSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  const auto instances = ctx.pfc->instances_with("type", std::string(kTypeName));
  std::vector<Named> result;
  result.reserve(instances.size());
  for (const auto& name : instances) {
    (void) ctx.pfc->get<std::string>(name + ".type");  // mark consumed
    const Statistics stat = ParseStatistics(ctx.pfc, name);
    result.push_back({name,
                      std::make_unique<DrPsdSpecies>(ctx.pfc,
                                                     name,
                                                     *ctx.ncdm_settings,
                                                     ctx.pba,
                                                     ctx.bgm,
                                                     stat)});
  }
  return result;
}

// ─────────────────────────────────────────────────────────────────────────────
// Grid + thermal seed
// ─────────────────────────────────────────────────────────────────────────────

void DrPsdSpecies::BuildLogTrapezoidGrids() {
  const int N    = N_q_bg_;
  const double h = std::log(q_max_bg_ / q_min_bg_) / (N - 1);
  q_bg_.resize(N);
  dq_bg_.resize(N);
  for (int i = 0; i < N; ++i)
    q_bg_[i] = q_min_bg_ * std::exp(i * h);
  // Log-trapezoid cell widths (dq = q du), half-weight endpoints (thesis §3).
  for (int i = 0; i < N; ++i) {
    const double wq = (i == 0 || i == N - 1) ? 0.5 : 1.0;
    dq_bg_[i]       = wq * h * q_bg_[i];
  }
  // Perturbation grid: EXACT integer subsample, q_[j] read out of q_bg_[j·refine_].
  // Not a spline (background quantities ring across the daughter's ~87-decade
  // injection cliff — see ComputeBackground) and not an independently log-spaced
  // grid (two exp() formulas agree only to rounding; direct indexing makes a PT
  // point literally a BG array entry, which is what lets every background column be
  // read at PT index j as BG column j·refine_).
  // N points span N−1 intervals ⇒ N_bg−1 = refine_·(N_pt−1), enforced in the ctor.
  const int N_pt    = (N - 1) / refine_ + 1;
  const double h_pt = refine_ * h;  // same log-trapezoid convention, coarser step
  N_q_              = N_pt;
  q_.resize(N_pt);
  dq_.resize(N_pt);
  for (int j = 0; j < N_pt; ++j) {
    q_[j]           = q_bg_[j * refine_];
    const double wq = (j == 0 || j == N_pt - 1) ? 0.5 : 1.0;
    dq_[j]          = wq * h_pt * q_[j];
  }
  w_.assign(N_pt, 0.);
  dlnf0_dlnq_.assign(N_pt, 0.);
  w_bg_.assign(N, 0.);
}

void DrPsdSpecies::BuildLinearTrapezoidGrids() {
  // quadrature.c `qm_trapz`, qmin==0 branch: trapezoid over [0, q_max] with q=0
  // carried as an implicit zero-weight point (q²f vanishes there), so the FIRST
  // interior point takes full weight h and only the last is halved. Byte-identical
  // to both reference implementations and to np.linspace(0, q_max, N+1)[1:].
  // dr_q_min is deliberately unused: the grid is anchored at zero and its first
  // point is q_max/N by construction.
  const int N    = N_q_bg_;
  const double h = q_max_bg_ / N;
  q_bg_.resize(N);
  dq_bg_.resize(N);
  for (int i = 0; i < N; ++i) {
    q_bg_[i]  = h + i * h;
    dq_bg_[i] = (i == N - 1) ? 0.5 * h : h;
  }
  q_bg_[N - 1] = q_max_bg_;  // make the endpoint exact, as the reference codes do

  // Perturbation grid: EXACT integer subsample taken from the TOP,
  // q_[j] = q_bg_[(j+1)*refine_ - 1] = (j+1)*(refine_*h). Anchoring on the implicit
  // zero (q_bg_[j*refine_], as the log grid does) would leave a first cell of width
  // h among cells of width refine_*h — not a trapezoid rule; this offset makes the
  // coarse grid a canonical qm_trapz grid in its own right, with step refine_*h and
  // the same q_max. N divisible by refine_ is enforced in the ctor.
  const int N_pt    = N / refine_;
  const double h_pt = refine_ * h;
  N_q_              = N_pt;
  q_.resize(N_pt);
  dq_.resize(N_pt);
  for (int j = 0; j < N_pt; ++j) {
    q_[j]  = q_bg_[(j + 1) * refine_ - 1];
    dq_[j] = (j == N_pt - 1) ? 0.5 * h_pt : h_pt;
  }

  w_.assign(N_q_, 0.);
  dlnf0_dlnq_.assign(N_q_, 0.);
  w_bg_.assign(N, 0.);
}

double DrPsdSpecies::ThermalF0(double q) const {
  const double x = q - ksi_;
  if (stat_ == Statistics::Fermion)
    return 1.0 / (std::exp(x) + 1.0);
  // Bose: regular as long as x > 0 (ksi < q_min); clamp defensively.
  const double denom = std::expm1(x);
  return (denom > kFFloor) ? 1.0 / denom : 1.0 / kFFloor;
}

// ─────────────────────────────────────────────────────────────────────────────
// Background
// ─────────────────────────────────────────────────────────────────────────────

void DrPsdSpecies::RegisterBackgroundIndices(int& index_bg) {
  index_bg_number_    = index_bg++;
  index_bg_rho_       = index_bg++;
  index_bg_p_         = index_bg++;
  index_bg_pseudo_p_  = index_bg++;
  index_bg_f_         = index_bg;
  index_bg           += N_q_bg_;
  index_bg_dfdlnq_    = index_bg;
  index_bg           += N_q_bg_;
}

void DrPsdSpecies::RegisterIntegrationIndices(int& index_bi) {
  index_bi_f_  = index_bi;
  index_bi    += N_q_bg_;
}

void DrPsdSpecies::SetBackgroundInitialConditions(const BackgroundICContext& ctx) {
  for (int i = 0; i < N_q_bg_; ++i) {
    // Seed at the positivity floor (not a large subtracted offset): the per-bin f̄
    // feeds ∂f̄/∂lnq, the kernel's band factor / (1±f) coefficients and ρ,p, so f̄
    // must be the TRUE occupation, not `bi_f − kFSeed`. A large seed (kFSeed=1e-10)
    // recovered by subtraction loses all precision once the true occupation falls
    // below ulp(seed) ≈2e-26 (high-q daughter bins), leaving f̄ as rounding noise.
    // Floored directly, no subtraction.
    // f_ini_ is the paper's per-daughter amplitude A (0 = starts empty).
    const double f0                           = (f_ini_ > 0.) ? f_ini_ * ThermalF0(q_bg_[i]) : 0.;
    ctx.pvecback_integration[index_bi_f_ + i] = f0 + kFFloor;
  }
}

void DrPsdSpecies::ComputeBackground(double a, const double* pvecback_B, double* pvecback) {
  const int N = N_q_bg_;
  std::vector<double> f_df(2 * N), ddf(N), lnq(N);
  for (int i = 0; i < N; ++i) {
    // Floor at kFFloor (gains dominate; f stays strictly positive). No seed
    // subtraction: f̄ is the true occupation (see the IC seed).
    const double f_raw = pvecback_B[index_bi_f_ + i];
    // The clamp below is what makes a negative occupation INVISIBLE downstream, so
    // report it here or it cannot be observed at all: every published column, and
    // every diagnostic built on one, sees kFFloor and not the negative value. This
    // is the detector for the non-Metzler loss (decay_transition_kernel.h,
    // DecayTransitionKernel::LumpTheta) -- without it, "no negatives" only means
    // "the clamp ran".
    // Gated on StoringBackgroundTable so the evolver's trial states, which ARE
    // allowed to be unphysical, do not trip it.
    if (f_raw < 0. && bgm_ != nullptr && bgm_->StoringBackgroundTable()) {
      // Record the WORST excursion relative to the PSD's own scale on that row, not
      // merely the sign. A bin sitting at -1e-120 while the peak is 1e-1 is the
      // integrator's round-off against a floored zero and is physically irrelevant;
      // a bin at -1e-3 of the peak is the non-Metzler loss actually running the
      // state negative. Reporting only "a negative occurred" cannot tell those apart,
      // and the difference is the whole positivity question (see LumpTheta).
      double peak = 0.;
      for (int j = 0; j < N; ++j)
        peak = std::max(peak, pvecback_B[index_bi_f_ + j]);
      const double rel = (peak > 0.) ? -f_raw / peak : 0.;
      if (rel > worst_negative_rel_) {
        worst_negative_rel_ = rel;
        // Announce only once per DECADE of severity, so a run that merely grazes
        // round-off stays quiet while one that is genuinely running the state
        // negative reports each time it gets worse.
        if (rel > 1e-8 && rel > 10. * announced_negative_rel_) {
          announced_negative_rel_ = rel;
          fprintf(stderr,
                  "WARNING: species '%s': daughter occupation NEGATIVE at %.2e of the "
                  "PSD peak (bin %d, q=%g, a=%g, f=%g). Clamped to %g downstream, so "
                  "nothing else will show it.\n",
                  name().c_str(),
                  rel,
                  i,
                  q_bg_[i],
                  a,
                  f_raw,
                  kFFloor);
        }
      }
    }
    const double f            = std::max(f_raw, kFFloor);
    w_bg_[i]                  = f * dq_bg_[i];
    pvecback[index_bg_f_ + i] = f;
    f_df[i]                   = f;
    lnq[i]                    = std::log(q_bg_[i]);
  }
  // Live ∂f̄/∂lnq via an f–ln q spline. Splining f̄ ITSELF (not ln f̄, the DNCDM/NCDM
  // seam) is load-bearing for a species sourced from zero: a decay-filled daughter
  // PSD is a narrow bump sitting on the kFFloor plateau, i.e. an ~87-DECADE cliff in
  // one bin. In ln f̄ that cliff makes ∂lnf̄/∂lnq ~ 1e2-1e3 and rings into the
  // physically-occupied bins, so f̄·∂lnf̄/∂lnq overstates ∂f̄/∂lnq by ~2 orders right
  // where the metric drives the F-hierarchy — and the bump MOVES as the decay band
  // shifts with a, so the ringing pattern jumps from step to step and the RHS is
  // effectively discontinuous in time (ndf15 step-size collapse, rkdp45 ringing).
  // On f̄ itself the floor is simply ≈0: the edge is an ordinary "falls to zero" and
  // ∂f̄/∂lnq stays bounded by the bump scale. This is the same normalized-quantity
  // ill-conditioning that motivates evolving F instead of Ψ (see the perturbation
  // section below) — ∂lnf̄/∂lnq is exactly such a normalized quantity.
  array_spline_table_lines(lnq.data(), N, f_df.data(), 1, ddf.data(), _SPLINE_EST_DERIV_);
  array_derive_spline(lnq.data(), N, f_df.data(), ddf.data(), 1, 0, N);
  for (int i = 0; i < N; ++i)
    pvecback[index_bg_dfdlnq_ + i] = f_df[N + i];

  double number, rho, p, pseudo_p;
  ComputeMomenta(1. / a - 1., &number, &rho, &p, nullptr, &pseudo_p);  // M_ = 0 ⇒ ε = q
  pvecback[index_bg_number_]   = number;
  pvecback[index_bg_rho_]      = rho;
  pvecback[index_bg_p_]        = p;
  pvecback[index_bg_pseudo_p_] = pseudo_p;
}

void DrPsdSpecies::WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const {
  const std::string nm = ColumnLabel();
  w.Add("(.)number_dr_" + nm, 0.);
  w.Add("(.)rho_dr_" + nm, 0.);
  w.Add("(.)p_dr_" + nm, 0.);
  for (int i = 0; i < N_q_bg_; ++i)
    w.Add("f_dr_" + nm + "[" + std::to_string(i) + "]", 0.);
}

void DrPsdSpecies::WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const {
  const std::string nm = ColumnLabel();
  w.Add("(.)number_dr_" + nm, Number(pvecback));
  w.Add("(.)rho_dr_" + nm, Rho(pvecback));
  w.Add("(.)p_dr_" + nm, P(pvecback));
  for (int i = 0; i < N_q_bg_; ++i)
    w.Add("f_dr_" + nm + "[" + std::to_string(i) + "]", pvecback[index_bg_f_ + i]);
}

// ─────────────────────────────────────────────────────────────────────────────
// Perturbations — UNNORMALIZED F ≡ f̄·Ψ = δf hierarchy (NOT the standard CLASS
// normalized Ψ = δf/f̄ that DNCDMSpecies uses).
//
// WHY (load-bearing, not a style choice): a daughter is sourced from EXACTLY zero
// abundance, so its f̄̇/f̄ is smooth but legitimately enormous at early times
// (measured 2.1e9 at a=3.3e-14, decaying to ~10 by a=6e-6). Both the −(f̄̇/f̄)Ψ
// dilution term and the collision operator's 1/f̄ normalization divide by that, so
// any residual Ψ blows the hierarchy up long before the ratio settles (the
// inverse-decay NaN). Substituting F = f̄Ψ into the Ψ-hierarchy and using
// f̄·(∂lnf̄/∂lnq) = ∂f̄/∂lnq exactly, the f̄̇Ψ terms cancel IDENTICALLY:
//   Ḟ + i(qk/ε)μF + (∂f̄/∂lnq)[η̇ − μ²(ḣ+6η̇)/2] = (dF/dτ)_C
// — no f̄̇/f̄ term exists to divide by, and F=0 is the exact statement "no
// perturbation to an empty distribution". The free-streaming ℓ,ℓ±1 propagation
// coefficients are linear and representation-independent, so they are UNCHANGED
// from the Ψ form (identical in FORM to DNCDMSpecies with M_ = 0, ε = q ⇒ qk/ε = k);
// only the metric-source coefficient changes (dlnf0/dlnq → f̄·dlnf0/dlnq) and the
// moment weights drop their f̄ (already baked into the state variable).
// ─────────────────────────────────────────────────────────────────────────────

void DrPsdSpecies::SetReducedBasis(int n_moments,
                                   const std::vector<double>& psi,
                                   const std::vector<double>& gnorm,
                                   double q_ref,
                                   double alpha0) {
  reduced_moments_ = n_moments;
  red_psi_         = psi;
  red_gnorm_       = gnorm;
  red_q_ref_       = q_ref;
  red_alpha0_      = alpha0;
}

void DrPsdSpecies::ReducedDriver(const double* pvecback, double* D) const {
  // The metric driver, projected: D_j = Sum dq q^2 psi_j (d f-bar / d ln q). Read on the
  // BACKGROUND grid through bg_index, exactly as the q-resolved path does -- the exact
  // integer subsample guarantees the mapped column is the value AT q_[iq].
  const int nq = q_size();
  for (int j = 0; j < reduced_moments_; ++j)
    D[j] = 0.;
  for (int iq = 0; iq < nq; ++iq) {
    const double base = dq_[iq] * q_[iq] * q_[iq] * pvecback[index_bg_dfdlnq_ + bg_index(iq)];
    for (int j = 0; j < reduced_moments_; ++j)
      D[j] += base * red_psi_[(size_t) j * nq + iq];
  }
}

void DrPsdSpecies::ReducedProject(const double* fq, double* m) const {
  const int nq = q_size();
  for (int j = 0; j < reduced_moments_; ++j)
    m[j] = 0.;
  for (int iq = 0; iq < nq; ++iq) {
    const double base = dq_[iq] * q_[iq] * q_[iq] * fq[iq];
    for (int j = 0; j < reduced_moments_; ++j)
      m[j] += base * red_psi_[(size_t) j * nq + iq];
  }
}

void DrPsdSpecies::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                               perturb_vector* /*pv*/,
                                               const precision* ppr,
                                               int& index_pt,
                                               const perturb_workspace* /*ppw*/,
                                               int /*gauge*/) {
  auto& layout = static_cast<NCDMBaseSpecies::PerturbLayout&>(base);

  // No fluid approximation (design §5 excludes dr_psd): always full hierarchy.
  layout.l_max = ppr->l_max_ncdm;
  // In the reduced representation the "q" index counts MOMENTS, not grid points. That
  // is the whole cost saving, and it is safe because eps = q makes the free-streaming
  // coefficient q-independent, so every moment obeys the same hierarchy exactly.
  layout.q_size = reduced() ? reduced_moments_ : q_size();

  layout.index_per_q.clear();
  layout.index_per_q.reserve(layout.q_size);
  for (int iq = 0; iq < layout.q_size; ++iq)
    layout.index_per_q.push_back(index_pt + iq * (layout.l_max + 1));

  index_pt += layout.total_size();
}

void DrPsdSpecies::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                                 double /*tau*/,
                                 const double* y,
                                 double* dy,
                                 const perturb_parameters_and_workspace& ppaw) const {
  const auto& layout = static_cast<const NCDMBaseSpecies::PerturbLayout&>(base);
  if (layout.q_size <= 0 || layout.index_per_q.empty())
    return;

  const perturb_workspace* ppw    = ppaw.ppw;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;
  const double* s_l               = ppw->s_l.data();
  const double k                  = ctx.k;
  const double a2                 = ctx.a2;
  const double metric_continuity  = ctx.metric_continuity;
  const double metric_euler       = ctx.metric_euler;
  const double metric_shear       = ctx.metric_shear;
  const double cotKgen            = ctx.cotKgen;

  const double* pvecback = ppw->pvecback.data();
  const int lmax         = layout.l_max;

  if (reduced()) {
    // The reduced hierarchy is the SAME hierarchy. eps = q makes qk/eps exactly k, so
    // the transport carries no q and every moment free-streams identically -- there is
    // no closure error here at all, only in the collision. What changes is the metric
    // driver: the projection of d f-bar / d ln q onto the frozen basis, in place of its
    // pointwise value at q_[iq].
    double D[kMaxReducedMoments] = {0.};
    ReducedDriver(pvecback, D);
    for (int j = 0; j < layout.q_size; ++j) {
      const int idx = layout.index_per_q[j];
      dy[idx]       = -k * y[idx + 1] + metric_continuity * D[j] / 3.;
      dy[idx + 1]   = k / 3. * (y[idx] - 2. * s_l[2] * y[idx + 2]) - metric_euler / (3. * k) * D[j];
      dy[idx + 2]   = k / 5. * (2. * s_l[2] * y[idx + 1] - 3. * s_l[3] * y[idx + 3]) -
                      s_l[2] * metric_shear * 2. / 15. * D[j];
      for (int l = 3; l < lmax; ++l)
        dy[idx + l] = k / (2. * l + 1.) *
                      (l * s_l[l] * y[idx + l - 1] - (l + 1.) * s_l[l + 1] * y[idx + l + 1]);
      dy[idx + lmax] = k * y[idx + lmax - 1] - (1. + lmax) * k * cotKgen * y[idx + lmax];
    }
    return;
  }

  for (int iq = 0; iq < layout.q_size; ++iq) {
    const double q = q_[iq];
    // F-space metric driver: ∂f̄/∂lnq, NOT the dimensionless ∂lnf̄/∂lnq the
    // Ψ-hierarchy uses. Splined from f̄ directly (see ComputeBackground) — never
    // reconstructed as f̄·∂lnf̄/∂lnq, which the floor plateau corrupts.
    // iq is a PT-grid index; the column is BG-grid, so map it (the exact subsample
    // guarantees the mapped column is the value AT q_[iq], not near it).
    const double df0_dlnq = pvecback[index_bg_dfdlnq_ + bg_index(iq)];

    const double epsilon        = std::sqrt(q * q + a2 * M_ * M_);  // = q
    const double qk_div_epsilon = k * q / epsilon;
    const int idx               = layout.index_per_q[iq];

    dy[idx]     = -qk_div_epsilon * y[idx + 1] + metric_continuity * df0_dlnq / 3.;
    dy[idx + 1] = qk_div_epsilon / 3. * (y[idx] - 2. * s_l[2] * y[idx + 2]) -
                  epsilon * metric_euler / (3. * q * k) * df0_dlnq;
    dy[idx + 2] = qk_div_epsilon / 5. * (2. * s_l[2] * y[idx + 1] - 3. * s_l[3] * y[idx + 3]) -
                  s_l[2] * metric_shear * 2. / 15. * df0_dlnq;
    for (int l = 3; l < lmax; ++l)
      dy[idx + l] = qk_div_epsilon / (2. * l + 1.) *
                    (l * s_l[l] * y[idx + l - 1] - (l + 1.) * s_l[l + 1] * y[idx + l + 1]);
    dy[idx + lmax] = qk_div_epsilon * y[idx + lmax - 1] - (1. + lmax) * k * cotKgen * y[idx + lmax];
  }
}

void DrPsdSpecies::ApplyInitialConditions(const BaseSpecies::PerturbLayout& base,
                                          double* y,
                                          const PerturbIcContext& ctx) {
  const auto& layout = static_cast<const NCDMBaseSpecies::PerturbLayout&>(base);
  if (layout.q_size <= 0 || layout.index_per_q.empty())
    return;

  // F = f̄·Ψ_adiabatic, which the ∂f̄/∂lnq driver already carries — the standard
  // adiabatic seed with ∂lnf̄/∂lnq swapped for ∂f̄/∂lnq. No empty-bin special case is
  // needed (and none is wanted): on an empty bin ∂f̄/∂lnq ≈ 0, so F = 0 falls out —
  // "no perturbation to a distribution that is empty anyway", which is exactly what
  // the old Ψ = 0 branch faked.
  const double* pvecback = ctx.ppw->pvecback.data();

  if (reduced()) {
    // Same adiabatic seed, projected. eps/q = 1 for a massless daughter, so the l = 1
    // coefficient loses its only q-dependence and the seed is a pure multiple of D_j.
    double D[kMaxReducedMoments] = {0.};
    ReducedDriver(pvecback, D);
    for (int j = 0; j < layout.q_size; ++j) {
      const int idx  = layout.index_per_q[j];
      const int lmax = layout.l_max;
      y[idx + 0]     = -0.25 * ctx.delta_ur * D[j];
      if (lmax >= 1)
        y[idx + 1] = -1. / 3. / ctx.k * ctx.theta_ur * D[j];
      if (lmax >= 2)
        y[idx + 2] = -0.5 * ctx.shear_ur * D[j];
      if (lmax >= 3)
        y[idx + 3] = -0.25 * ctx.l3_ur * D[j];
    }
    return;
  }

  for (int index_q = 0; index_q < layout.q_size; ++index_q) {
    const int idx  = layout.index_per_q[index_q];
    const int lmax = layout.l_max;

    const double q       = q_[index_q];
    const double epsilon = std::sqrt(q * q + ctx.a * ctx.a * M_ * M_);  // = q
    // PT index -> BG column (see PerturbDerivs).
    const double df0_dlnq = pvecback[index_bg_dfdlnq_ + bg_index(index_q)];  // ∂f̄/∂lnq

    y[idx + 0] = -0.25 * ctx.delta_ur * df0_dlnq;
    if (lmax >= 1)
      y[idx + 1] = -epsilon / 3. / q / ctx.k * ctx.theta_ur * df0_dlnq;
    if (lmax >= 2)
      y[idx + 2] = -0.5 * ctx.shear_ur * df0_dlnq;
    if (lmax >= 3)
      y[idx + 3] = -0.25 * ctx.l3_ur * df0_dlnq;
  }
}

BaseSpecies::StressEnergyContribution DrPsdSpecies::StressEnergy(
    const BaseSpecies::PerturbLayout& base,
    const perturb_vector* /*pv*/,
    const double* y,
    const double* pvecback,
    const perturb_workspace* ppw) const {
  const auto& layout = static_cast<const NCDMBaseSpecies::PerturbLayout&>(base);
  StressEnergyContribution se;
  se.rho = pvecback[index_bg_rho_];
  se.p   = pvecback[index_bg_p_];
  if (layout.q_size <= 0 || layout.index_per_q.empty())
    return se;

  const double k  = ppw->scalar_ctx.k;
  const double a2 = ppw->scalar_ctx.a2;

  // UNNORMALIZED weights: w0 = dq only. The state IS δf, so the f̄ that the
  // normalized-Ψ convention carries in the weight is already baked into y.
  double s_drho = 0., s_theta = 0., s_dp = 0., s_shear = 0.;

  if (reduced()) {
    // For a massless daughter eps = q, so ALL FOUR sums carry the same weight q^3 --
    // they are the ENERGY moment at l = 0, 1, 2 and nothing else. That is why the metric
    // never needs the daughter's q-resolution, which is the design's whole premise
    // (section 2.2), and it is exact rather than an approximation.
    s_drho              = ReducedEnergyAtL(layout, y, 0);
    s_dp                = s_drho;
    s_theta             = ReducedEnergyAtL(layout, y, 1);
    s_shear             = ReducedEnergyAtL(layout, y, 2);
    const double fac    = factor_ / (a2 * a2);
    se.delta_rho        = fac * s_drho;
    se.rho_plus_p_theta = fac * k * s_theta;
    se.delta_p          = fac * s_dp / 3.;
    se.rho_plus_p_shear = fac * (2. / 3.) * s_shear;
    return se;
  }

  for (int iq = 0; iq < layout.q_size; ++iq) {
    const double q        = q_[iq];
    const double q2       = q * q;
    const double epsilon  = std::sqrt(q2 + a2 * M_ * M_);  // = q
    const int idx         = layout.index_per_q[iq];
    const double w0       = dq_[iq];
    s_drho               += w0 * q2 * epsilon * y[idx];
    s_theta              += w0 * q2 * q * y[idx + 1];
    s_dp                 += w0 * q2 * q2 / epsilon * y[idx];
    s_shear              += w0 * q2 * q2 / epsilon * y[idx + 2];
  }
  const double fac    = factor_ / (a2 * a2);
  se.delta_rho        = fac * s_drho;
  se.rho_plus_p_theta = fac * k * s_theta;
  se.delta_p          = fac * s_dp / 3.;
  se.rho_plus_p_shear = fac * (2. / 3.) * s_shear;
  return se;
}

void DrPsdSpecies::FillSources(const BaseSpecies::PerturbLayout& base,
                               const double* y,
                               const double* /*dy*/,
                               PerturbSourceContext& ctx) const {
  if (ctx.index_md != ctx.p_mod->index_md_scalars_)
    return;
  if (index_tp_delta_ < 0 && index_tp_theta_ < 0)
    return;

  const auto& layout     = static_cast<const NCDMBaseSpecies::PerturbLayout&>(base);
  const double* pvecback = ctx.ppw->pvecback.data();
  const double rho       = Rho(pvecback);
  const double p         = P(pvecback);
  const double a2        = ctx.a2;
  const double k         = ctx.k;

  double s_drho = 0., s_theta = 0.;
  if (reduced() && layout.q_size > 0 && !layout.index_per_q.empty()) {
    s_drho  = ReducedEnergyAtL(layout, y, 0);
    s_theta = ReducedEnergyAtL(layout, y, 1);
  }
  else if (layout.q_size > 0 && !layout.index_per_q.empty()) {
    for (int iq = 0; iq < layout.q_size; ++iq) {
      const double q        = q_[iq];
      const double q2       = q * q;
      const double epsilon  = std::sqrt(q2 + a2 * M_ * M_);  // = q
      const int idx         = layout.index_per_q[iq];
      const double w0       = dq_[iq];  // unnormalized state: y IS δf (see StressEnergy)
      s_drho               += w0 * q2 * epsilon * y[idx];
      s_theta              += w0 * q2 * q * y[idx + 1];
    }
  }
  const double fac = factor_ / (a2 * a2);

  if (index_tp_delta_ >= 0) {
    const double delta = (rho > 0.) ? fac * s_drho / rho : 0.;
    const double src =
        (rho > 0.) ? delta - RhoPrimeOverRho(pvecback, ctx.a_prime_over_a) * ctx.theta_over_k2 : 0.;
    ctx.p_mod->SetSourceValue(ctx.index_md,
                              ctx.index_ic,
                              index_tp_delta_,
                              ctx.index_tau,
                              ctx.index_k,
                              src);
  }
  if (index_tp_theta_ >= 0) {
    const double rho_plus_p = rho + p;
    const double theta      = (rho_plus_p > 0.) ? fac * k * s_theta / rho_plus_p : 0.;
    const double src        = (rho_plus_p > 0.) ? theta + ctx.theta_shift : 0.;
    ctx.p_mod->SetSourceValue(ctx.index_md,
                              ctx.index_ic,
                              index_tp_theta_,
                              ctx.index_tau,
                              ctx.index_k,
                              src);
  }
}

void DrPsdSpecies::WriteOutputColumns(PerturbColumnWriter& w,
                                      const PerturbationsModule& mod,
                                      file_format fmt,
                                      BaseSpecies::TransferColumnSection section) const {
  if (fmt != file_format::class_format)
    return;
  const perturbs* ppt = mod.GetPerturbs();
  if (section != TransferColumnSection::velocity && ppt->has_density_transfers)
    w.Add("d_dr_" + ColumnLabel(), index_tp_delta_, index_tp_delta_ >= 0);
  if (section != TransferColumnSection::density && ppt->has_velocity_transfers)
    w.Add("t_dr_" + ColumnLabel(), index_tp_theta_, index_tp_theta_ >= 0);
}
