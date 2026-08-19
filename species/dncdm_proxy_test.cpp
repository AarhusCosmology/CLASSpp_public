// Unit test for DNCDMProxySpecies (dr_representation = proxy): the two closed-form
// pieces of the relaxation-time closure, checked against arXiv:2203.09075 rather
// than against this code's own output, plus the dr_rta_form selector.
//
// Chen, Oldengott, Pierobon & Wong give 𝓕 twice -- once exactly (eq. 14, in terms
// of the incomplete gamma function) and once as the two-branch approximation they
// actually integrate, with stated accuracies. TransportRate ships the second. So
// the exact form is reimplemented here, independently, and the shipped branches are
// held to the paper's own error bars. That is the check a transcription bug in a
// seven-term series fails and a self-comparison would not.
#include <cassert>
#include <cmath>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "background.h"
#include "parser.h"
#include "precision.h"
#include "species/dncdm_proxy_species.h"
#include "species/ncdm_base_species.h"
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

static background MakeBackground() {
  background pba{};
  pba.H0       = 2.2e-4;
  pba.Omega0_g = 5.4e-5;  // finite radiation so the ctor's k_rad guess is regular
  return pba;
}

// A proxy instance at the shipped defaults. Deliberately small grids: nothing here
// integrates anything, the tests only read back parsed configuration.
static void SetProxyBase(FileContent& fc) {
  fc.set("dncdm1.type", "ncdm_decay_dr");
  fc.set("dncdm1.dr_representation", "proxy");
  fc.set("dncdm1.m", "1.0");
  fc.set("dncdm1.Gamma", "10");
  fc.set("dncdm1.Omega_dncdmdr", "0.001");
  fc.set("dncdm1.quadrature_strategy", "3");
  fc.set("dncdm1.momenta_bins", "16");
  fc.set("dncdm1.dr_N_q", "16");
}

static std::unique_ptr<DNCDMProxySpecies> BuildProxy(FileContent& fc,
                                                     background& pba,
                                                     NcdmSettings& settings) {
  auto parent = std::make_unique<DNCDMSpecies>(&fc, "dncdm1", settings, &pba, nullptr);
  // Defaults, but non-null: Create reads ppr->l_max_ncdm / l_max_dr to check them
  // against the RTA closure's fixed multipole scratch, and the production path always
  // supplies a precision. Leaving it null here made the fixture unlike production.
  static precision ppr_defaults{};
  SpeciesBuildContext ctx{};
  ctx.pfc           = &fc;
  ctx.pba           = &pba;
  ctx.ncdm_settings = &settings;
  ctx.ppr           = &ppr_defaults;
  ctx.bgm           = nullptr;
  Named n           = DNCDMProxySpecies::Create(std::move(parent), ctx);
  auto* raw         = dynamic_cast<DNCDMProxySpecies*>(n.species.get());
  assert(raw != nullptr);
  n.species.release();
  return std::unique_ptr<DNCDMProxySpecies>(raw);
}

// ─────────────────────────────────────────────────────────────────────────────
// Reference implementations, written from the paper and used by nothing else.
// ─────────────────────────────────────────────────────────────────────────────

/** E_1(x) = Γ(0,x), the exponential integral. Abramowitz & Stegun 5.1.11 (series)
 *  below x = 1 and 5.1.22 (Lentz continued fraction) above it -- the standard split,
 *  because the alternating series loses all its digits to cancellation once x > 1.
 *  Independent of anything in the shipped translation unit. */
static double E1(double x) {
  assert(x > 0.);
  constexpr double kEulerGamma = 0.5772156649015329;
  if (x <= 1.) {
    // E_1(x) = -gamma - ln x + sum_{n>=1} (-1)^(n+1) x^n / (n n!)
    double sum = 0., term = 1.;
    for (int n = 1; n <= 60; ++n) {
      term *= -x / n;
      sum  -= term / n;
    }
    return -kEulerGamma - std::log(x) + sum;
  }
  // E_1(x) = e^-x / (x + 1 - 1/(x + 3 - 4/(x + 5 - ...))), by modified Lentz.
  constexpr double kTiny = 1e-300;
  double b = x + 1., c = 1. / kTiny, d = 1. / b, h = d;
  for (int i = 1; i <= 200; ++i) {
    const double an   = -1. * i * i;
    b                += 2.;
    d                 = 1. / (an * d + b);
    c                 = b + an / c;
    const double del  = c * d;
    h                *= del;
    if (std::fabs(del - 1.) < 1e-15)
      break;
  }
  return h * std::exp(-x);
}

/** arXiv:2203.09075 eq. (14) verbatim: 𝓕(x) = ½ e^-x [ -1 + x - e^x (x²-2) Γ(0,x) ]. */
static double CurlyFExact(double x) {
  return 0.5 * std::exp(-x) * (-1. + x - std::exp(x) * (x * x - 2.) * E1(x));
}

static int failures = 0;

static void Check(bool ok, const char* what) {
  if (!ok) {
    std::fprintf(stderr, "FAIL: %s\n", what);
    ++failures;
  }
}

static void CheckClose(double got, double want, double rtol, const char* what) {
  const double denom = std::max(std::fabs(want), 1e-300);
  const double rel   = std::fabs(got - want) / denom;
  if (!(rel <= rtol)) {
    std::fprintf(stderr,
                 "FAIL: %s: got %.10g, want %.10g, rel %.3g > %.3g\n",
                 what,
                 got,
                 want,
                 rel,
                 rtol);
    ++failures;
  }
}

int main() {
  // ── 𝓕: the shipped two-branch approximation against eq. (14) ───────────────
  //
  // The paper's own error bars, and they are the assertion: "<7% for x <= 2" on the
  // small-x series, "<10% at 2 < x < 17" on the tail. A mistyped coefficient in the
  // series breaks the first at small x; a wrong tail coefficient breaks the second.
  {
    const double small[] = {1e-8, 1e-6, 1e-4, 1e-3, 0.01, 0.1, 0.3, 0.5, 1.0, 1.5, 1.9, 2.0};
    for (double x : small) {
      char msg[128];
      std::snprintf(msg, sizeof(msg), "CurlyF small-x branch at x = %g", x);
      CheckClose(DNCDMProxySpecies::CurlyF(x), CurlyFExact(x), 0.07, msg);
    }
    const double large[] = {2.0001, 2.5, 3., 5., 8., 12., 16.9};
    for (double x : large) {
      char msg[128];
      std::snprintf(msg, sizeof(msg), "CurlyF large-x branch at x = %g", x);
      CheckClose(DNCDMProxySpecies::CurlyF(x), CurlyFExact(x), 0.10, msg);
    }
  }

  // 𝓕 is positive and decreasing over the range the rate is ever evaluated on, and
  // 𝓕 ~ -ln x as x -> 0 (the paper: "for x ~ 1e-10 -> 0.1, 𝓕 evaluates to ~20 -> 1").
  // The rate carries x^5 𝓕(x), so the log divergence is integrable and must NOT be
  // clamped away -- a floor here would silently flatten the small-x rate.
  {
    // The paper's two stated magnitudes are round numbers, so they are checked as
    // such: eq. (14) gives 21.96 and 1.407 at these points (the series agrees), and
    // the paper rounds them to "~20" and "~1".
    Check(DNCDMProxySpecies::CurlyF(1e-10) > 15. && DNCDMProxySpecies::CurlyF(1e-10) < 25.,
          "CurlyF(1e-10) is the paper's ~20");
    Check(DNCDMProxySpecies::CurlyF(0.1) > 1. && DNCDMProxySpecies::CurlyF(0.1) < 2.,
          "CurlyF(0.1) is the paper's ~1");
    double prev = DNCDMProxySpecies::CurlyF(1e-6);
    for (double x = 1e-5; x < 20.; x *= 1.3) {
      const double cur = DNCDMProxySpecies::CurlyF(x);
      Check(cur > 0. && cur < prev, "CurlyF is positive and strictly decreasing");
      prev = cur;
    }
    // Guarded, not evaluated, outside the domain: x <= 0 has no 𝓕, and past x = 40
    // the parent is so non-relativistic that inverse decay is shut.
    Check(DNCDMProxySpecies::CurlyF(0.) == 0., "CurlyF(0) is guarded to zero");
    Check(DNCDMProxySpecies::CurlyF(-1.) == 0., "CurlyF(negative) is guarded to zero");
    Check(DNCDMProxySpecies::CurlyF(41.) == 0., "CurlyF past the cutoff is zero");
  }

  // ── α_ℓ: eq. (16) ─────────────────────────────────────────────────────────
  //
  // α_0 = α_1 = 0 is the load-bearing part -- the ℓ <= 1 collision integrals vanish
  // identically by energy/momentum conservation, so a nonzero α there would source
  // the Einstein equations with a term that must not exist.
  {
    CheckClose(DNCDMProxySpecies::AlphaL(0), 0., 0., "AlphaL(0) = 0");
    CheckClose(DNCDMProxySpecies::AlphaL(1), 0., 0., "AlphaL(1) = 0");
    CheckClose(DNCDMProxySpecies::AlphaL(2), 1., 1e-14, "AlphaL(2) = 1");
    for (int l = 0; l <= 20; ++l) {
      const double x    = static_cast<double>(l);
      const double want = (3. * x * x * x * x + 2. * x * x * x - 11. * x * x + 6. * x) / 32.;
      char msg[64];
      std::snprintf(msg, sizeof(msg), "AlphaL(%d) matches eq. (16)", l);
      CheckClose(DNCDMProxySpecies::AlphaL(l), want, 1e-14, msg);
    }
    for (int l = 2; l < 20; ++l)
      Check(DNCDMProxySpecies::AlphaL(l + 1) > DNCDMProxySpecies::AlphaL(l),
            "AlphaL increases with l above the conservation floor");
  }

  // ── dr_rta_form ───────────────────────────────────────────────────────────
  {
    NcdmSettings settings = TestSettings();
    background pba        = MakeBackground();

    FileContent fc;
    SetProxyBase(fc);
    Check(BuildProxy(fc, pba, settings)->rta_form() == DNCDMProxySpecies::RtaForm::kPowers,
          "dr_rta_form defaults to the fitted two-power form");

    FileContent fc_powers;
    SetProxyBase(fc_powers);
    fc_powers.set("dncdm1.dr_rta_form", "powers");
    Check(BuildProxy(fc_powers, pba, settings)->rta_form() == DNCDMProxySpecies::RtaForm::kPowers,
          "dr_rta_form = powers selects the fitted form");

    FileContent fc_copw;
    SetProxyBase(fc_copw);
    fc_copw.set("dncdm1.dr_rta_form", "copw");
    Check(BuildProxy(fc_copw, pba, settings)->rta_form() == DNCDMProxySpecies::RtaForm::kCOPW,
          "dr_rta_form = copw selects the analytic form");

    // A misspelt form must not silently fall back to the default: the two forms
    // differ by more than the effect being measured, so a typo that ran anyway
    // would be indistinguishable from a physics result.
    FileContent fc_bad;
    SetProxyBase(fc_bad);
    fc_bad.set("dncdm1.dr_rta_form", "paper");
    bool threw = false;
    try {
      BuildProxy(fc_bad, pba, settings);
    }
    catch (const std::invalid_argument&) {
      threw = true;
    }
    Check(threw, "an unknown dr_rta_form is rejected rather than defaulted");
  }

  // ── TransportRate under copw IS eq. (13) ──────────────────────────────────
  //
  // The formula is short enough that the only way to get it wrong is a stray factor,
  // which is exactly what a stray factor looks like: nothing crashes and the rate is
  // merely wrong. So it is checked against the closed form on a hand-filled
  // background row, with the scale factor passed in so no BackgroundModule is needed.
  //
  //   Gamma_T,l = alpha_l * a * Gamma^0 * (rho_H/rho_sec) * (1/12) * X^5 * F(X)
  //
  // The defaults are load-bearing here: eq. (13) has no free normalisation, so this
  // only holds if dr_rta_C5 defaults to 1 under copw rather than to the fitted
  // 2.5228, and if no gamma^-3 term is added on top.
  {
    NcdmSettings settings = TestSettings();
    background pba        = MakeBackground();
    FileContent fc;
    SetProxyBase(fc);
    fc.set("dncdm1.dr_rta_form", "copw");
    auto proxy = BuildProxy(fc, pba, settings);

    int index_bg = 0;
    proxy->RegisterBackgroundIndices(index_bg);
    std::vector<double> pvecback(index_bg + 8, 0.);

    // A row on which every guard in TransportRate passes: a live parent holding a
    // decent share of the sector, and a mildly out-of-balance background.
    //
    // GetMass() is m_nuH/T_0 (~6e3 here), not the eV mass, so the row is built from
    // it rather than from round numbers: a is chosen to put X = a*m at 1, where 𝓕 is
    // O(1) and the seam is nearby, and the comoving number is chosen to give a
    // definite mean Lorentz factor, since gamma = rho_H/(number*m) identically.
    const double m       = proxy->parent().GetMass();
    const double G0      = proxy->parent().Gamma();
    const double a       = 1.0 / m;  // => X = 1
    const double X       = a * m;
    const double rho_H   = 3.0;
    const double rho_sec = 10.0;
    const double gamma   = 3.0;
    const double number  = rho_H / (gamma * m);

    pvecback[proxy->parent().bg_rho_index()]    = rho_H;
    pvecback[proxy->parent().bg_number_index()] = number;
    pvecback[proxy->bg_rho_sec_index()]         = rho_sec;
    pvecback[proxy->bg_eps_ne_index()]          = 0.5;

    // If this trips, the row above is wrong, not the code: gamma < 1 is unphysical
    // and TransportRate clamps it, which would make the comparison meaningless.
    CheckClose((a * rho_H / number) / X,
               gamma,
               1e-12,
               "the hand-filled row has the intended gamma");

    for (int l : {2, 3, 5, 9}) {
      const double want = DNCDMProxySpecies::AlphaL(l) * a * G0 * (rho_H / rho_sec) * (1. / 12.) *
                          X * X * X * X * X * DNCDMProxySpecies::CurlyF(X);
      char msg[96];
      std::snprintf(msg, sizeof(msg), "copw TransportRate(l=%d) is eq. (13)", l);
      CheckClose(proxy->TransportRate(l, a, pvecback.data()), want, 1e-12, msg);
    }
    // l <= 1 vanishes identically, and alpha_l already encodes that.
    for (int l : {0, 1})
      CheckClose(proxy->TransportRate(l, a, pvecback.data()), 0., 0., "copw rate is 0 for l <= 1");

    // The fitted form must NOT agree: if these coincide, the selector is not wired
    // through and every comparison made with it would be vacuous.
    FileContent fcp;
    SetProxyBase(fcp);
    auto powers = BuildProxy(fcp, pba, settings);
    int idx_p   = 0;
    powers->RegisterBackgroundIndices(idx_p);
    std::vector<double> pv_p(idx_p + 8, 0.);
    pv_p[powers->parent().bg_rho_index()]    = rho_H;
    pv_p[powers->parent().bg_number_index()] = number;
    pv_p[powers->bg_rho_sec_index()]         = rho_sec;
    pv_p[powers->bg_eps_ne_index()]          = 0.5;
    Check(std::fabs(powers->TransportRate(4, a, pv_p.data()) -
                    proxy->TransportRate(4, a, pvecback.data())) > 1e-30,
          "the two forms give different rates on the same row");
  }

  if (failures > 0) {
    std::fprintf(stderr, "dncdm_proxy_test: %d failure(s)\n", failures);
    return 1;
  }
  std::printf("dncdm_proxy_test: all checks passed\n");
  return 0;
}
