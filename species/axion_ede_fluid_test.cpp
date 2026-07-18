#include "axion_ede_fluid.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

#include "background.h"
#include "parser.h"
#include "species_build_context.h"

// Log-spaced trapezoid quadrature of 3(1+w)/a on [a, 1] for the reference integral.
static double QuadIntegral(const AxionEDEFluid& fld, double a) {
  const int N     = 200000;  // log-spaced, plenty for 1e-14 -> 1
  const double la = std::log(a);
  const double h  = -la / N;
  auto f          = [&](double lna) {
    double w, dw, integ;
    fld.ComputeWFld(std::exp(lna), &w, &dw, &integ);
    return 3. * (1. + w);
  };
  double s = 0.5 * (f(la) + f(0.));
  for (int i = 1; i < N; ++i)
    s += f(la + i * h);
  return s * h;
}

int main() {
  background pba{};
  pba.H0 = 2.2e-4;

  const double a_c = std::pow(10., -3.5), nu = 1., w_i = -1., n = 3.;
  const double w_f = AxionEDEFluid::WFinal(n);
  assert(std::fabs(w_f - 0.5) < 1e-15);

  AxionEDEFluid fld(pba, /*omega0=*/1e-6, a_c, n, nu, w_i, w_f, /*theta_i=*/2.8);

  // ── w(a): sigmoid limits and midpoint.
  double w, dw, integ;
  fld.ComputeWFld(1e-8, &w, &dw, &integ);
  assert(std::fabs(w - w_i) < 1e-10);
  fld.ComputeWFld(1.0, &w, &dw, &integ);
  assert(std::fabs(w - w_f) < 1e-3);  // a >> a_c
  fld.ComputeWFld(a_c, &w, &dw, &integ);
  assert(std::fabs(w - 0.5 * (w_i + w_f)) < 1e-12);  // x = 1 at a = a_c

  // ── a = 0 exactly must not NaN (BackgroundModule::background_checks calls it).
  fld.ComputeWFld(0., &w, &dw, &integ);
  assert(w == w_i);
  assert(dw == 0.);
  assert(std::isfinite(integ) || integ == 0.);

  // ── dw/da vs central finite difference over the transition.
  for (double a : {a_c / 10., a_c / 2., a_c, 2. * a_c, 10. * a_c}) {
    const double h = a * 1e-6;
    double wp, wm, d_;
    fld.ComputeWFld(a + h, &wp, &d_, &integ);
    fld.ComputeWFld(a - h, &wm, &d_, &integ);
    const double fd = (wp - wm) / (2. * h);
    fld.ComputeWFld(a, &w, &dw, &integ);
    assert(std::fabs(dw - fd) < 1e-4 * std::fabs(fd) + 1e-12);
  }

  // ── Closed-form integral vs quadrature at several epochs.
  for (double a : {1e-10, 1e-5, a_c, 1e-2, 0.5}) {
    fld.ComputeWFld(a, &w, &dw, &integ);
    const double q = QuadIntegral(fld, a);
    assert(std::fabs(integ - q) < 1e-4 * std::fabs(q) + 1e-8);
  }

  // ── Omega0 <-> Omega_ac round trip.
  const double om_ac = 0.05;
  const double om0   = AxionEDEFluid::OmegaZeroFromOmegaAc(om_ac, a_c, nu, w_i, w_f);
  assert(std::fabs(AxionEDEFluid::OmegaAcFromOmegaZero(om0, a_c, nu, w_i, w_f) - om_ac) <
         1e-12 * om_ac);
  // And it matches direct density evolution: rho(a_c)/rho(1) = exp(integral(a_c)).
  fld.ComputeWFld(a_c, &w, &dw, &integ);
  assert(std::fabs(om0 * std::exp(integ) - om_ac) < 1e-10 * om_ac);

  // ── cs2 limits: k -> infinity gives 1; k -> 0 gives w_f = (n-1)/(n+1);
  //    monotonic in k^2. (omega_axion_ must be set: use the test hook.)
  fld.SetOmegaAxionForTest(1e2 * pba.H0);
  const double a_test = 1e-2;
  assert(std::fabs(fld.Cs2(1e30, a_test) - 1.) < 1e-6);
  assert(std::fabs(fld.Cs2(1e-30, a_test) - (n - 1.) / (n + 1.)) < 1e-6);
  double prev = 0.;
  for (double k2 : {1e-8, 1e-4, 1., 1e4, 1e8}) {
    const double c = fld.Cs2(k2, a_test);
    assert(c >= prev - 1e-15 && c <= 1. + 1e-15);
    prev = c;
  }

  // ════ CreateAll wiring ════
  auto make_base_fc = [] {
    FileContent fc;
    fc.set("fluid_equation_of_state", "pheno_axion");
    fc.set("n_pheno_axion", "3");
    fc.set("log10_axion_ac", "-3.5");
    fc.set("Theta_initial_fld", "2.8");
    return fc;
  };
  auto run_create = [&pba](FileContent& fc) {
    SpeciesBuildContext ctx{};
    ctx.pfc = &fc;
    ctx.pba = &pba;
    return FluidSpecies::CreateAll(ctx);
  };

  // Valid: Omega_fld_ac path.
  {
    FileContent fc = make_base_fc();
    fc.set("Omega_fld_ac", "0.05");
    auto result = run_create(fc);
    assert(result.size() == 1);
    auto* ax = dynamic_cast<AxionEDEFluid*>(result[0].species.get());
    assert(ax != nullptr);
    const double ac  = std::pow(10., -3.5);
    const double om0 = AxionEDEFluid::OmegaZeroFromOmegaAc(0.05, ac, 1., -1., 0.5);
    assert(std::fabs(ax->GetOmega0() - om0) < 1e-12 * om0);
    assert(std::fabs(ax->a_c() - ac) < 1e-15);
    assert(std::fabs(ax->n_axion() - 3.) < 1e-15);
  }

  // Valid: direct Omega_fld path.
  {
    FileContent fc = make_base_fc();
    fc.set("Omega_fld", "1e-6");
    auto result = run_create(fc);
    assert(result.size() == 1);
    assert(std::fabs(result[0].species->GetOmega0() - 1e-6) < 1e-18);
  }

  // Valid: w_fld_f instead of n (n derived as (1+w_f)/(1-w_f)).
  {
    FileContent fc;
    fc.set("fluid_equation_of_state", "pheno_axion");
    fc.set("w_fld_f", "0.5");
    fc.set("log10_axion_ac", "-3.5");
    fc.set("Theta_initial_fld", "2.8");
    fc.set("Omega_fld_ac", "0.05");
    auto result = run_create(fc);
    assert(result.size() == 1);
    auto* ax = dynamic_cast<AxionEDEFluid*>(result[0].species.get());
    assert(std::fabs(ax->n_axion() - 3.) < 1e-12);
  }

  // Rejections.
  auto expect_create_throw = [&](FileContent fc) {
    bool threw = false;
    try {
      run_create(fc);
    }
    catch (const std::invalid_argument&) {
      threw = true;
    }
    assert(threw);
  };
  // Structural checks abort (invalid_argument); checks on a parsed numeric value's
  // range are computation rejections (runtime_error) so a sampler can reject the
  // point and keep the chain alive.
  auto expect_create_throw_computation = [&](FileContent fc) {
    bool threw = false;
    try {
      run_create(fc);
    }
    catch (const std::runtime_error&) {
      threw = true;
    }
    assert(threw);
  };
  {  // both n and w_fld_f
    FileContent fc = make_base_fc();
    fc.set("Omega_fld_ac", "0.05");
    fc.set("w_fld_f", "0.5");
    expect_create_throw(std::move(fc));
  }
  {  // missing Theta_initial_fld
    FileContent fc;
    fc.set("fluid_equation_of_state", "pheno_axion");
    fc.set("n_pheno_axion", "3");
    fc.set("log10_axion_ac", "-3.5");
    fc.set("Omega_fld_ac", "0.05");
    expect_create_throw(std::move(fc));
  }
  {  // both a_c and log10_axion_ac
    FileContent fc = make_base_fc();
    fc.set("Omega_fld_ac", "0.05");
    fc.set("a_c", "3e-4");
    expect_create_throw(std::move(fc));
  }
  {  // two density keys
    FileContent fc = make_base_fc();
    fc.set("Omega_fld_ac", "0.05");
    fc.set("Omega_fld", "1e-6");
    expect_create_throw(std::move(fc));
  }
  {  // no density key
    FileContent fc = make_base_fc();
    expect_create_throw(std::move(fc));
  }
  {  // explicit use_ppf = yes
    FileContent fc = make_base_fc();
    fc.set("Omega_fld_ac", "0.05");
    fc.set("use_ppf", "yes");
    expect_create_throw(std::move(fc));
  }
  {  // cs2_fld makes no sense here
    FileContent fc = make_base_fc();
    fc.set("Omega_fld_ac", "0.05");
    fc.set("cs2_fld", "1.0");
    expect_create_throw(std::move(fc));
  }
  {  // Theta out of range
    FileContent fc = make_base_fc();
    fc.set("Omega_fld_ac", "0.05");
    fc.set("Theta_initial_fld", "3.15");
    expect_create_throw_computation(std::move(fc));
  }
  {  // closure-budget override is rejected for pheno_axion
    FileContent fc = make_base_fc();
    fc.set("Omega_fld_ac", "0.05");
    SpeciesBuildContext ctx{};
    ctx.pfc                     = &fc;
    ctx.pba                     = &pba;
    ctx.omega0_closure_override = 0.7;
    bool threw                  = false;
    try {
      FluidSpecies::CreateAll(ctx);
    }
    catch (const std::invalid_argument&) {
      threw = true;
    }
    assert(threw);
  }
  {  // nu_fld must be > 0
    FileContent fc = make_base_fc();
    fc.set("Omega_fld_ac", "0.05");
    fc.set("nu_fld", "0");
    expect_create_throw_computation(std::move(fc));
  }
  {  // n_pheno_axion must be >= 1
    FileContent fc = make_base_fc();
    fc.set("Omega_fld_ac", "0.05");
    fc.set("n_pheno_axion", "0.5");
    expect_create_throw_computation(std::move(fc));
  }
  {  // w_fld_f must lie in [0, 1)
    FileContent fc;
    fc.set("fluid_equation_of_state", "pheno_axion");
    fc.set("w_fld_f", "1.0");
    fc.set("log10_axion_ac", "-3.5");
    fc.set("Theta_initial_fld", "2.8");
    fc.set("Omega_fld_ac", "0.05");
    expect_create_throw_computation(std::move(fc));
  }
  {  // negative w_fld_f would give derived n = (1+w_f)/(1-w_f) < 1
    FileContent fc;
    fc.set("fluid_equation_of_state", "pheno_axion");
    fc.set("w_fld_f", "-0.5");
    fc.set("log10_axion_ac", "-3.5");
    fc.set("Theta_initial_fld", "2.8");
    fc.set("Omega_fld_ac", "0.05");
    expect_create_throw_computation(std::move(fc));
  }
  {  // w_fld_f (= 0.5 from n = 3) must exceed w_fld_i
    FileContent fc = make_base_fc();
    fc.set("Omega_fld_ac", "0.05");
    fc.set("w_fld_i", "0.9");
    expect_create_throw_computation(std::move(fc));
  }
  {  // w_fld_i < -1 would cross the phantom divide, contradicting
     // AxionEDEFluid::ReachesPhantomDivide() == false
    FileContent fc = make_base_fc();
    fc.set("Omega_fld_ac", "0.05");
    fc.set("w_fld_i", "-1.5");
    expect_create_throw_computation(std::move(fc));
  }
  {  // a_c must lie in (0, 1)
    FileContent fc;
    fc.set("fluid_equation_of_state", "pheno_axion");
    fc.set("n_pheno_axion", "3");
    fc.set("a_c", "1.5");
    fc.set("Theta_initial_fld", "2.8");
    fc.set("Omega_fld_ac", "0.05");
    expect_create_throw_computation(std::move(fc));
  }
  {  // Omega_fld_ac must be positive
    FileContent fc = make_base_fc();
    fc.set("Omega_fld_ac", "-0.05");
    expect_create_throw_computation(std::move(fc));
  }
  {  // Omega_fld must be positive
    FileContent fc = make_base_fc();
    fc.set("Omega_fld", "-1e-6");
    expect_create_throw_computation(std::move(fc));
  }

  // ════ Phantom-divide + HyRec-CPL species hooks ════
  {
    // AxionEDEFluid: never reaches the divide; refuses a CPL pair for HyRec.
    assert(!fld.ReachesPhantomDivide());
    double w0h = -2., wah = -2.;
    assert(!fld.HyrecCplApproximation(&w0h, &wah));

    // Base fluid, CLP w0 = -1, wa = 0 (exact divide contact): reaches it.
    FluidSpecies lambda_like(pba, 0.7, CLP, -1., 0., 1., 0.);
    assert(lambda_like.ReachesPhantomDivide());

    // Base fluid, CLP w0 = -0.9, wa = 0: never reaches -1.
    FluidSpecies quint(pba, 0.7, CLP, -0.9, 0., 1., 0.);
    assert(!quint.ReachesPhantomDivide());
    // Its HyRec CPL tangent reproduces (w0, wa) exactly.
    assert(quint.HyrecCplApproximation(&w0h, &wah));
    assert(std::fabs(w0h + 0.9) < 1e-15);
    assert(std::fabs(wah - 0.) < 1e-15);

    // Base fluid, CLP crossing (w(0) = w0+wa = -1.4, w(1) = -0.9): reaches it.
    FluidSpecies crosser(pba, 0.7, CLP, -0.9, -0.5, 1., 0.);
    assert(crosser.ReachesPhantomDivide());
  }

  std::printf("axion EDE fluid tests passed\n");
  return 0;
}
