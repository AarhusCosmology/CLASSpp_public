#include "greybody_ncdm_species.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "background_module.h"
#include "bisection.h"
#include "errors.h"
#include "species/ncdm_family.h"
#include "species/species_input.h"

namespace {
constexpr double kGreyBodyLeadingFactor = 2.0;  // the leading 2 in f0

// Riemann zeta via the Borwein eta (alternating-series) algorithm. Vendored
// because libc++ lacks std::riemann_zeta; used uniformly on all compilers so
// results are toolchain-independent. Valid for s != 1.
double riemann_zeta(double s) {
  // Target near-machine-epsilon accuracy: the Borwein error scales like
  // (3 + sqrt(8))^(-n), so n = 21 - floor(s) suffices for s <= 20. Clamp to 6
  // (ample for large s, where zeta(s) -> 1 exponentially fast).
  int n = (s > 0.) ? static_cast<int>(21 - s) : 21;
  if (n < 6)
    n = 6;

  // Borwein d_k coefficients via the e_j auxiliary array.
  std::vector<double> e(n + 1);
  double bnj = 1.0;
  e[n]       = bnj;
  for (int j = n - 1; j >= 0; --j) {
    bnj  *= (j + 1.) / (n - j);
    e[j]  = e[j + 1] + bnj;
  }

  double S1 = 0., S2 = 0.;
  for (int k = 1; k <= n; ++k) {
    if ((k - 1) % 2 == 0)
      S1 += 1. / std::pow(k, s);
    else
      S1 -= 1. / std::pow(k, s);
  }
  for (int k = n + 1; k <= 2 * n; ++k) {
    if ((k - 1) % 2 == 0)
      S2 += e[k - n] / std::pow(k, s);
    else
      S2 -= e[k - n] / std::pow(k, s);
  }
  return (S1 + S2 / std::pow(2.0, n)) / (1. - std::pow(2.0, 1. - s));
}

// J_alpha = (1 - 2^-(alpha+1)) zeta(alpha+2), with Taylor fallback near alpha=-1.
double J_alpha(double alpha) {
  if (std::fabs(alpha + 1.) < 1e-3) {
    double j = alpha + 1., j2 = j * j, j3 = j2 * j, j4 = j3 * j, j5 = j4 * j, j6 = j5 * j;
    return 0.693147180559945309 + 0.159868903742430972 * j - 0.0326862962794492996 * j2 +
           0.00156899170541551496 * j3 + 0.000749872421120475325 * j4 -
           0.000204290897136748191 * j5 + 0.0000231747166672581420 * j6;
  }
  double A = std::pow(2., alpha + 1.);
  return (1. - 1. / A) * riemann_zeta(alpha + 2.);
}

// r_M(alpha) = M2*M4/M3^2 as a function of alpha only.
double r_M(double alpha) {
  // At alpha = -1 the prefactor 4A^2-5A+1 vanishes (A=1) while zeta(alpha+2)
  // diverges (zeta(1)); the 0*inf limit is finite. Reachable as a bisect
  // midpoint when lo+hi rounds to -2 exactly.
  if (alpha == -1.) {
    return 216. * std::log(2.) * riemann_zeta(3.) / std::pow(_PI_, 4);
  }
  double A = std::pow(2., alpha + 1.);
  return ((4. * A * A - 5. * A + 1.) * (alpha + 3.) * riemann_zeta(alpha + 2.) *
          riemann_zeta(alpha + 4.)) /
         (std::pow(2. * A - 1., 2) * (alpha + 2.) * std::pow(riemann_zeta(alpha + 3.), 2));
}

// r_M(0) = 21 zeta(2) zeta(4) / (18 zeta(3)^2): the moment ratio at alpha = 0,
// the anchor of the rational initial guess in solve_alpha.
constexpr double kRMatZero = 1.43748132827993497652;

// Bracket the root of r_M(alpha) = r, then bisect with tools/bisection.h.
// PRECONDITION: r > 1 (Cauchy-Schwarz: M2*M4 >= M3^2). r_M is monotonic,
// spanning (1, inf) as alpha runs from +inf down to -2, so every r > 1 has a
// unique solution; r <= 1 is unphysical (and would make the initial guess and
// bracket search ill-defined), so we reject it up front.
double solve_alpha(double r) {
  class_test(r <= 1., "grey-body: moment ratio r = %g must exceed 1 (requires M2*M4 >= M3^2)", r);
  auto f = [&](double a) { return r_M(a) - r; };
  // Initial guess from inverting the rational approximation
  // r_M(alpha) ~ (alpha + 2*kRMatZero) / (alpha + 2) (exact at alpha=0, -> 1 as
  // alpha -> inf): alpha ~ 2*(kRMatZero - r)/(r - 1).
  double a    = 2. * (kRMatZero - r) / (r - 1.);
  double left = a, right = a, lim = 1.;
  while (f(left) < 0.) {  // walk left toward alpha = -2
    lim   *= 0.1;
    right  = left;
    left   = -2. + lim;
  }
  while (f(right) > 0.) {  // walk right
    left  = right;
    right = std::fabs(right) * 10.;
  }
  // f(left) >= 0, f(right) <= 0; bisect_value treats pred(mid)==true as the hi side.
  class_test(f(left) * f(right) > 0., "grey-body: failed to bracket alpha for r=%g", r);
  return bisect_value(left, right, 1e-13, [&](double mid) { return f(mid) < 0.; });
}
}  // namespace

namespace greybody {

GreyBodyParams GreyBodyParams::FromDirect(double alpha, double x, double q0) {
  // The PSD q^(alpha-1)/(exp(alpha*x*q)+1) and its GB-Laguerre quadrature
  // (which divides by alpha*x) require strictly positive parameters; reject
  // non-physical direct input rather than silently producing inf/nan.
  class_test(alpha <= 0. || x <= 0. || q0 <= 0.,
             "grey-body: direct parameters alpha, x, q0 must all be > 0");
  GreyBodyParams p;
  p.alpha_         = alpha;
  p.x_             = x;
  p.q0_            = q0;
  p.x_times_alpha_ = x * alpha;
  p.alpham1_logq0_ = (alpha - 1.) * std::log(q0);
  return p;
}

GreyBodyParams GreyBodyParams::FromMoments(double r, double M2, double M3) {
  // x divides by M3 and the normalization takes log(M2); both moments of a
  // positive PSD are strictly positive, so reject non-positive input up front.
  class_test(M2 <= 0. || M3 <= 0., "grey-body: moments M2 and M3 must both be > 0");
  double alpha = solve_alpha(r);
  double A     = std::pow(2., alpha + 1.);
  double x     = ((2. * A - 1.) * M2 * riemann_zeta(alpha + 3.) * std::tgamma(alpha + 3.)) /
                 (2. * (A - 1.) * alpha * M3 * riemann_zeta(alpha + 2.) * std::tgamma(alpha + 2.));
  double x_times_alpha = x * alpha;
  double alpham1_logq0 = std::log(J_alpha(alpha)) + std::lgamma(alpha + 2.) -
                         (alpha + 2.) * std::log(x_times_alpha) - std::log(M2);
  double q0            = std::min(1e100, std::exp(alpham1_logq0 / (alpha - 1.)));

  GreyBodyParams p;
  p.alpha_         = alpha;
  p.x_             = x;
  p.q0_            = q0;
  p.x_times_alpha_ = x_times_alpha;
  p.alpham1_logq0_ = alpham1_logq0;
  return p;
}

void GreyBodyParams::moments(double& M2, double& M3, double& M4) const {
  // q0^(1-alpha) == exp(-alpham1_logq0_) by definition of alpham1_logq0_ =
  // (alpha-1)*log(q0). We use the latter (the same quantity the PSD evaluates)
  // so the moments stay consistent with the PSD even when q0_ itself is the
  // clamped/degenerate value reported near the alpha -> 1 (Fermi-Dirac) limit,
  // where q0 is not identifiable.
  auto Mn = [&](int n) {
    double s = n + alpha_;
    return std::exp(-alpham1_logq0_) * (1. - std::pow(2., 1. - s)) * std::tgamma(s) *
           riemann_zeta(s) / std::pow(x_times_alpha_, s);
  };
  M2 = Mn(2);
  M3 = Mn(3);
  M4 = Mn(4);
}

}  // namespace greybody

GreyBodyNCDMSpecies::GreyBodyNCDMSpecies(FileContent* pfc,
                                         const std::string& instance_name,
                                         const NcdmSettings& settings,
                                         const background* pba,
                                         const BackgroundModule* bgm)
    : NCDMSpecies(pfc, instance_name, settings, pba, bgm, NCDMBaseSpecies::DeferInit{}) {
  SpeciesInput input(pfc, instance_name);

  auto mode_opt = input.get<std::string>("parameterization");
  class_test_severe(!mode_opt,
                    "grey-body species '%s': required field 'parameterization' is missing (set "
                    "it to 'direct' or 'moments')",
                    instance_name.c_str());
  const std::string& mode = *mode_opt;
  class_test_severe(mode != "direct" && mode != "moments",
                    "grey-body species '%s': parameterization must be 'direct' or 'moments'",
                    instance_name.c_str());

  if (mode == "direct") {
    auto alpha = input.get<double>("alpha");
    auto q0    = input.get<double>("q0");
    auto x     = input.get<double>("x");
    class_test_severe(!(alpha && q0 && x),
                      "grey-body species '%s': direct mode requires alpha, q0, x",
                      instance_name.c_str());
    gb_ = greybody::GreyBodyParams::FromDirect(*alpha, *x, *q0);
  }
  else {
    auto r  = input.get<double>("r");
    auto M2 = input.get<double>("M2");
    auto M3 = input.get<double>("M3");
    class_test_severe(!(r && M2 && M3),
                      "grey-body species '%s': moments mode requires r, M2, M3",
                      instance_name.c_str());
    gb_ = greybody::GreyBodyParams::FromMoments(*r, *M2, *M3);

    // Optional mass-via-M2 convenience: m = m_M2 / M2.
    if (auto m_M2 = input.get<double>("m_M2")) {
      m_in_eV_ = *m_M2 / *M2;
    }
    class_test(m_in_eV_ != 0.0 && m_in_eV_ < 0.01,
               "grey-body species '%s': inferred mass < 0.01 eV; moment inversion is unreliable "
               "in the low-mass limit",
               instance_name.c_str());
  }
  gb_ready_ = true;

  // Object is now fully configured: build quadrature with the grey-body PSD
  // (dispatched via the virtual EvaluatePsdAnalytic) and resolve mass/Omega.
  BuildQuadratureAndMass(settings);
  ResolveMassOmegaClosure(settings);
}

double GreyBodyNCDMSpecies::EvaluatePsdAnalytic(double q) const {
  if (!gb_ready_) {
    // Invariant: BuildQuadratureAndMass (the only caller path into here) is
    // deferred until gb_ is configured, so this is unreachable. Fail loudly
    // rather than silently returning the wrong (Fermi-Dirac) PSD.
    throw std::logic_error("GreyBodyNCDMSpecies::EvaluatePsdAnalytic called before gb_ready_");
  }
  // Log-space evaluation for numerical stability at large parameters. The
  // denominator log(exp(t)+1) is computed as the stable softplus t + log1p(e^-t)
  // (t = alpha*x*q >= 0 here, since all quadrature nodes q are positive) so the
  // exponential never overflows for large t.
  double t    = gb_.x_times_alpha() * q;
  double logf = (gb_.alpha() - 1.) * std::log(q) - gb_.alpham1_logq0() -
                (t + std::log1p(std::exp(-t)));
  return kGreyBodyLeadingFactor / std::pow(2. * _PI_, 3) * std::exp(logf);
}

quadrature_method GreyBodyNCDMSpecies::DefaultQuadratureStrategy() const {
  return qm_GB_Laguerre;
}

void GreyBodyNCDMSpecies::FillQuadratureParams(GBQuadParams& p) const {
  p.active        = true;
  p.alpha         = gb_.alpha();
  p.x_times_alpha = gb_.x_times_alpha();
}

void GreyBodyNCDMSpecies::GreyBodyMoments(double& M2, double& M3, double& M4) const {
  // Report the moments of the *actual* sampled distribution (q_bg_/w_bg_, whose
  // weights already fold in f0), not the analytic gb_.moments(): this reflects
  // what the background integration used, including quadrature resolution.
  const std::vector<double>& qv = q_bg_;
  const std::vector<double>& wv = w_bg_;
  M2 = M3 = M4 = 0.;
  for (size_t i = 0; i < qv.size(); ++i) {
    double q2  = qv[i] * qv[i];
    M2        += q2 * wv[i];
    M3        += q2 * qv[i] * wv[i];
    M4        += q2 * q2 * wv[i];
  }
  double pref  = 0.5 * std::pow(2. * _PI_, 3);
  M2          *= pref;
  M3          *= pref;
  M4          *= pref;
}

std::optional<double> GreyBodyNCDMSpecies::GetParam(const std::string& name) const {
  if (name == "alpha")
    return gb_.alpha();
  if (name == "q0")
    return gb_.q0();
  if (name == "x")
    return gb_.x();
  if (name == "M_2" || name == "M_3" || name == "M_4") {
    double M2, M3, M4;
    GreyBodyMoments(M2, M3, M4);
    if (name == "M_2")
      return M2;
    if (name == "M_3")
      return M3;
    return M4;
  }
  return BaseSpecies::GetParam(name);
}

std::vector<Named> GreyBodyNCDMSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  return CreateAllNcdmInstances<GreyBodyNCDMSpecies>(ctx);
}
