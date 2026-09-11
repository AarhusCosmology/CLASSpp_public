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

// As SetProxyBase, but the parent's normalisation comes from `deg` rather than
// Omega_dncdmdr: channel multiplicity on the parent side is only meaningful when the
// parent's weight counts SPECIES, which is the neutrino case the campaign runs.
static void SetProxyBaseDeg(FileContent& fc, const char* deg) {
  fc.set("dncdm1.type", "ncdm_decay_dr");
  fc.set("dncdm1.dr_representation", "proxy");
  fc.set("dncdm1.m", "1.0");
  fc.set("dncdm1.Gamma", "10");
  fc.set("dncdm1.T", "0.71611");
  fc.set("dncdm1.deg", deg);
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

    // dr_rta_alpha applies to BOTH forms -- AlphaFor multiplies the gamma^-5 term on
    // the common return path, whichever branch built it. That is deliberate (see
    // AlphaForm), so it is asserted rather than left to be rediscovered: copw out of
    // the box is the paper because the DEFAULT is the paper's quartic, and copw with
    // `integrated` is the same eq. (13) carrying the corrected l-dependence.
    {
      FileContent fc_ai;
      SetProxyBase(fc_ai);
      fc_ai.set("dncdm1.dr_rta_form", "copw");
      fc_ai.set("dncdm1.dr_rta_alpha", "integrated");
      auto copw_int = BuildProxy(fc_ai, pba, settings);
      // A SECOND instance registers its OWN background indices, so it needs its own
      // row: handing it the vector filled for `proxy` reads slots that mean something
      // else and every guard in TransportRate returns 0 -- silently, and the test
      // would then be asserting against a zero it had produced itself.
      int idx_ai = 0;
      copw_int->RegisterBackgroundIndices(idx_ai);
      std::vector<double> pv_ai(idx_ai + 8, 0.);
      pv_ai[copw_int->parent().bg_rho_index()]    = rho_H;
      pv_ai[copw_int->parent().bg_number_index()] = number;
      pv_ai[copw_int->bg_rho_sec_index()]         = rho_sec;
      pv_ai[copw_int->bg_eps_ne_index()]          = 0.5;
      for (int l : {2, 5, 9}) {
        const double base = a * G0 * (rho_H / rho_sec) * (1. / 12.) * X * X * X * X * X *
                            DNCDMProxySpecies::CurlyF(X);
        char msg[112];
        std::snprintf(msg, sizeof(msg), "copw + integrated is eq. (13) x alpha_eff (l=%d)", l);
        CheckClose(copw_int->TransportRate(l, a, pvecback.data()),
                   DNCDMProxySpecies::AlphaLEff(l, X) * base,
                   1e-12,
                   msg);
      }
      // ... and at l = 2 the two forms coincide, because alpha_2 = 1 either way.
      CheckClose(copw_int->TransportRate(2, a, pvecback.data()),
                 proxy->TransportRate(2, a, pvecback.data()),
                 1e-12,
                 "copw l=2 is unaffected by dr_rta_alpha");
    }

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

  // ── alpha_l^eff: the tabulated, momentum-integrated l-dependence ───────────
  //
  // The table is generated, so what is checked here are the properties that make it
  // usable at all -- not its values, which have no independent oracle.
  {
    // Normalised to l = 2 at every X, by construction. If this drifts, C5 has been
    // silently rescaled and every calibrated amplitude is wrong.
    for (double X : {1e-6, 1e-3, 0.01, 0.3, 1.0, 10.}) {
      char msg[128];
      std::snprintf(msg, sizeof(msg), "AlphaLEff(2, %g) == 1", X);
      CheckClose(DNCDMProxySpecies::AlphaLEff(2, X), 1.0, 1e-12, msg);
    }
    // Monotonically increasing in l wherever the term matters (X <~ 1), and strictly
    // BELOW the published quartic everywhere -- that is the whole content of the form.
    for (double X : {1e-6, 1e-3, 0.01, 0.1}) {
      double prev = 1.0;
      for (int l = 3; l <= 17; ++l) {
        const double v = DNCDMProxySpecies::AlphaLEff(l, X);
        char msg[128];
        std::snprintf(msg, sizeof(msg), "AlphaLEff increasing in l at X = %g, l = %d", X, l);
        Check(v > prev, msg);
        std::snprintf(msg, sizeof(msg), "AlphaLEff < AlphaL at X = %g, l = %d", X, l);
        Check(v < DNCDMProxySpecies::AlphaL(l), msg);
        prev = v;
      }
    }
    // Falling in X at fixed l: less boost, less able to wipe fine angular structure.
    Check(DNCDMProxySpecies::AlphaLEff(17, 1e-3) > DNCDMProxySpecies::AlphaLEff(17, 0.1),
          "AlphaLEff(17, .) decreasing in X");
    // Clamped, not extrapolated, outside the table.
    CheckClose(DNCDMProxySpecies::AlphaLEff(17, 1e-12),
               DNCDMProxySpecies::AlphaLEff(17, 1e-6),
               1e-12,
               "AlphaLEff clamps below");
    CheckClose(DNCDMProxySpecies::AlphaLEff(17, 1e6),
               DNCDMProxySpecies::AlphaLEff(17, 31.6),
               1e-2,
               "AlphaLEff clamps above");
    // Outside the tabulated l range it must fall back to the quartic, not to garbage.
    CheckClose(DNCDMProxySpecies::AlphaLEff(25, 0.01),
               DNCDMProxySpecies::AlphaL(25),
               1e-12,
               "AlphaLEff falls back to AlphaL above l_max");
    // The published values it is meant to replace, at the campaign's own X.
    CheckClose(DNCDMProxySpecies::AlphaL(17), 8041.0, 1e-9, "AlphaL(17) = 8041");
    Check(DNCDMProxySpecies::AlphaLEff(17, 0.01) < 0.5 * DNCDMProxySpecies::AlphaL(17),
          "AlphaLEff(17) is less than half the quartic at X = 0.01");
  }

  // ── the structured form's two shape functions ─────────────────────────────
  {
    // Phi_4 IS curly-F, evaluated exactly instead of through the paper's two-branch
    // approximation. CurlyFExact above is an independent implementation (it is what
    // the CurlyF tolerance checks are made against), so this pins the closed form
    // against something not derived from it.
    for (double X : {1e-4, 1e-2, 0.1, 0.3, 1.0, 2.0, 5.0, 12.0}) {
      char msg[128];
      std::snprintf(msg, sizeof(msg), "PhiNLO(%g) is curly-F exactly", X);
      CheckClose(DNCDMProxySpecies::PhiNLO(X), CurlyFExact(X), 1e-6, msg);
    }
    // Phi_2(0) = 1 is what makes X^3 Phi_2 -> X^3 propto gamma^-3: the
    // Hannestad-Raffelt asymptote is RECOVERED, not assumed. If this drifts the LO
    // term no longer has the limit it is named for.
    CheckClose(DNCDMProxySpecies::PhiLO(0.), 1.0, 1e-12, "PhiLO(0) = 1");
    CheckClose(DNCDMProxySpecies::PhiLO(1e-8), 1.0, 1e-7, "PhiLO -> 1 as X -> 0");
    // Closed form against values from an independent quadrature of the defining
    // integral (scipy, 8 digits): (1+X)e^-X - X^2 Gamma(0,X).
    const double xs[] = {0.1, 0.3, 1.0, 2.0, 5.0, 10.0};
    // mpmath, 30 digits, of the DEFINING integral -- not of the closed form, so this
    // is an independent check rather than a restatement. Quoted to 10 digits because
    // at 4 the reference itself was what failed the tolerance.
    const double want[] = {9.7709192026e-01,
                           8.8155278824e-01,
                           5.1637494795e-01,
                           2.1040380688e-01,
                           1.1720292213e-02,
                           8.3702334419e-05};
    for (int i = 0; i < 6; ++i) {
      char msg[128];
      std::snprintf(msg, sizeof(msg), "PhiLO(%g) against quadrature", xs[i]);
      CheckClose(DNCDMProxySpecies::PhiLO(xs[i]), want[i], 1e-9, msg);
    }
    // BOTH die once the parent is non-relativistic -- the whole point of the form.
    // A gamma-power model saturates here instead and isotropises forever.
    double pl = DNCDMProxySpecies::PhiLO(1.), pn = DNCDMProxySpecies::PhiNLO(1.);
    for (double X : {2., 4., 8., 16.}) {
      const double l2 = DNCDMProxySpecies::PhiLO(X), n2 = DNCDMProxySpecies::PhiNLO(X);
      char msg[128];
      std::snprintf(msg, sizeof(msg), "PhiLO decreasing at X = %g", X);
      Check(l2 < pl && l2 > 0., msg);
      std::snprintf(msg, sizeof(msg), "PhiNLO decreasing at X = %g", X);
      Check(n2 < pn && n2 > 0., msg);
      pl = l2;
      pn = n2;
    }
    Check(DNCDMProxySpecies::PhiLO(30.) < 1e-12 * DNCDMProxySpecies::PhiLO(1.),
          "PhiLO has collapsed by X = 30");
    // And the form is reachable.
    NcdmSettings settings = TestSettings();
    background pba        = MakeBackground();
    FileContent fcs;
    SetProxyBase(fcs);
    fcs.set("dncdm1.dr_rta_form", "structured");
    auto st = BuildProxy(fcs, pba, settings);
    Check(st->rta_form() == DNCDMProxySpecies::RtaForm::kStructured,
          "dr_rta_form = structured selects the structured form");

    // The ASSEMBLY, on a hand-filled row -- the Phi's being right does not make
    // (1/12) X^3 Phi_2 and (1/12) X^5 Phi_4 right, and this is the piece a refit
    // would silently invalidate.
    int idx_s = 0;
    st->RegisterBackgroundIndices(idx_s);
    std::vector<double> pv_s(idx_s + 8, 0.);
    const double aS = 1e-3, rhoS = 3.0, secS = 10.0, numS = 2.0, epsS = 0.25;
    pv_s[st->parent().bg_rho_index()]    = rhoS;
    pv_s[st->parent().bg_number_index()] = numS;
    pv_s[st->bg_rho_sec_index()]         = secS;
    pv_s[st->bg_eps_ne_index()]          = epsS;
    const double XS                      = aS * st->parent().GetMass();
    const double G0S                     = st->parent().Gamma();
    for (int l : {2, 4, 7}) {
      const double b3   = l * (l + 1.) / 6.;
      const double r3   = 0.0848 * std::pow(epsS, 0.5) * (1. / 12.) * XS * XS * XS *
                          DNCDMProxySpecies::PhiLO(XS);
      const double r5   = 0.5283 * (1. / 12.) * XS * XS * XS * XS * XS *
                          DNCDMProxySpecies::PhiNLO(XS);
      const double want = aS * G0S * (rhoS / secS) * (b3 * r3 + DNCDMProxySpecies::AlphaL(l) * r5);
      char msg[112];
      std::snprintf(msg, sizeof(msg), "structured TransportRate(l=%d) assembly", l);
      CheckClose(st->TransportRate(l, aS, pv_s.data()), want, 1e-12, msg);
    }
  }

  // ── beta_l: both forms normalised at l = 2, and distinct above ────────────
  {
    NcdmSettings settings = TestSettings();
    background pba        = MakeBackground();
    FileContent fcb;
    SetProxyBase(fcb);
    fcb.set("dncdm1.dr_rta_beta", "legs");
    auto legs = BuildProxy(fcb, pba, settings);
    FileContent fcl;
    SetProxyBase(fcl);
    auto leg = BuildProxy(fcl, pba, settings);

    int i1 = 0, i2 = 0;
    legs->RegisterBackgroundIndices(i1);
    leg->RegisterBackgroundIndices(i2);
    std::vector<double> p1(i1 + 8, 0.), p2(i2 + 8, 0.);
    const double aB = 2e-3, rB = 3.0, sB = 10.0, nB = 2.0;
    for (auto* pr : {&p1, &p2}) {
      auto& sp                             = (pr == &p1) ? *legs : *leg;
      (*pr)[sp.parent().bg_rho_index()]    = rB;
      (*pr)[sp.parent().bg_number_index()] = nB;
      (*pr)[sp.bg_rho_sec_index()]         = sB;
      (*pr)[sp.bg_eps_ne_index()]          = 0.9;  // far from balance: LO term lives
    }
    // Both are 1 at l = 2 by construction, so C3 does not move with the choice --
    // if that drifts, every calibrated amplitude silently changes meaning.
    CheckClose(legs->TransportRate(2, aB, p1.data()),
               leg->TransportRate(2, aB, p2.data()),
               1e-12,
               "beta forms agree at l = 2");
    // ... and differ above it, in the direction the O(mu) coefficient says:
    // l(l+1)/2 - 2 = 8 at l = 4 against l(l+1)/6 = 10/3.
    Check(legs->TransportRate(4, aB, p1.data()) > leg->TransportRate(4, aB, p2.data()),
          "dr_rta_beta = legs is steeper than legendre at l = 4");
  }

  // ── channel multiplicity: the weights, and the guard that keeps them honest ──
  //
  // Two separate things carry a multiplicity, and the whole design rests on not
  // confusing them (design 2026-09-11-dncdm-scenario-b-multiplicity):
  //   * the kernel's channel COUNT, which scales rates;
  //   * the species' momentum-integration WEIGHT, which scales rho/n/Pi.
  // The daughter has no `deg` of its own, so `dr_n_daughter` has to supply its weight;
  // the parent's weight is `deg`, which already works (measured: deg = 2 doubles
  // rho_dncdm1 exactly). What must not be spellable is the two disagreeing -- a parent
  // weighted 2 whose daughter is fed at the 1-parent rate runs, errors nothing, and
  // silently violates sector energy conservation.
  {
    background pba  = MakeBackground();
    NcdmSettings st = TestSettings();
    const int Nd    = 16;
    std::vector<double> f(Nd, 0.3);
    const double aD = 1e-3;

    FileContent f1{};
    SetProxyBaseDeg(f1, "1");
    auto one = BuildProxy(f1, pba, st);

    FileContent f2{};
    SetProxyBaseDeg(f2, "1");
    f2.set("dncdm1.dr_n_daughter", "2");
    auto two = BuildProxy(f2, pba, st);

    // The fermion daughter is n_daughter species; the boson is always one.
    CheckClose(two->DaughterRho(f.data(), aD, false),
               2. * one->DaughterRho(f.data(), aD, false),
               1e-14,
               "dr_n_daughter = 2 doubles the fermion daughter's rho");
    CheckClose(two->DaughterRho(f.data(), aD, true),
               one->DaughterRho(f.data(), aD, true),
               1e-14,
               "dr_n_daughter does not change the boson's rho (phi is one species)");

    // PARENT SIDE: rejected outright, at any deg. `deg` is overloaded and cannot carry
    // a species count.
    //
    // In the Omega_dncdmdr shooting path deg is a PSD AMPLITUDE -- a genuinely diluted
    // population, where f_bare = deg * f0 is the physical occupation and the kernel
    // boundary (KappaStoredToBare, which carries GetDeg()) is right to scale it. A
    // species COUNT is the opposite: n identical species share ONE per-dof occupation,
    // and scaling it corrupts the (1 -+ f) blocking terms, which are nonlinear in f.
    //
    // Measured, proxy background, m = 0.06, Gamma = 1e7, sector a^4 rho drift across the
    // decay epoch: 0.0043% at (deg 1, n_p 1), 0.0207% at (deg 2, n_p 1), and 3.64% at
    // (deg 2, n_p 2) -- an 850x leak on exactly the combination a naive guard would
    // bless. rho_phi, which starts empty and is therefore pure decay product, comes out
    // 7.4x at deg = 2 rather than the 2x a pure weight would give: proof the occupation
    // itself scaled.
    //
    // The daughter side has no such problem and is NOT rejected: n_daughter is applied
    // to rho/n outside the kernel and never reaches the boundary (measured drift 0.003%,
    // i.e. the baseline).
    for (const char* deg : {"1", "2"}) {
      FileContent fx{};
      SetProxyBaseDeg(fx, deg);
      fx.set("dncdm1.dr_n_parent", "2");
      bool rejected = false;
      try {
        auto bad = BuildProxy(fx, pba, st);
      }
      catch (const std::exception&) {
        rejected = true;
      }
      Check(rejected, "dr_n_parent > 1 is rejected (deg cannot carry a species count)");
    }
    // n_parent = 1 stays spellable, so the key is not simply dead.
    FileContent f5{};
    SetProxyBaseDeg(f5, "1");
    f5.set("dncdm1.dr_n_parent", "1");
    bool ok1 = true;
    try {
      auto fine = BuildProxy(f5, pba, st);
    }
    catch (const std::exception&) {
      ok1 = false;
    }
    Check(ok1, "dr_n_parent = 1 is accepted");
  }

  if (failures > 0) {
    std::fprintf(stderr, "dncdm_proxy_test: %d failure(s)\n", failures);
    return 1;
  }
  std::printf("dncdm_proxy_test: all checks passed\n");
  return 0;
}
