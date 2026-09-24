// The Jordan-frame scalar-tensor species: its factory's guards and shooting target,
// and the background closure it hands the BackgroundModule, checked against the
// field equations it solves.
// Design: docs/superpowers/specs/2026-09-24-scalar-tensor-gravity-design.md
#include <cassert>
#include <cmath>
#include <exception>
#include <functional>
#include <string>
#include <vector>

#include "background.h"
#include "parser.h"
#include "scalar_tensor_species.h"
#include "species_build_context.h"

namespace {

SpeciesBuildContext MakeCtx(FileContent& fc, const background& pba) {
  SpeciesBuildContext ctx{};
  ctx.pfc = &fc;
  ctx.pba = &pba;
  return ctx;
}

bool Throws(const std::function<void()>& f) {
  try {
    f();
  }
  catch (const std::exception&) {
    return true;
  }
  return false;
}

// A valid EMG-like input, to which a test adds or overrides keys.
FileContent Emg() {
  FileContent fc;
  fc.set("xi_st", "0.3");
  fc.set("phi_ini_st", "0.5");
  fc.set("lambda_st", "3e-105");
  return fc;
}

// The built species with its background slots and ODE slots registered from 0.
struct Built {
  std::unique_ptr<BaseSpecies> owner;
  ScalarTensorSpecies* st = nullptr;
  int bg_size             = 0;
  int bi_size             = 0;
};

Built Build(FileContent fc, const background& pba) {
  auto ctx    = MakeCtx(fc, pba);
  auto result = ScalarTensorSpecies::CreateAll(ctx);
  assert(result.size() == 1);
  Built b;
  b.owner = std::move(result[0].species);
  b.st    = static_cast<ScalarTensorSpecies*>(b.owner.get());
  b.st->RegisterBackgroundIndices(b.bg_size);
  b.st->RegisterIntegrationIndices(b.bi_size);
  return b;
}

}  // namespace

int main() {
  background pba{};
  pba.H0 = 2.2e-4;

  // ── Absent xi_st: no species.
  {
    FileContent fc;
    auto ctx = MakeCtx(fc, pba);
    assert(ScalarTensorSpecies::CreateAll(ctx).empty());
  }

  // ── phi_ini_st is required.
  {
    FileContent fc;
    fc.set("xi_st", "0.3");
    auto ctx = MakeCtx(fc, pba);
    assert(Throws([&] { ScalarTensorSpecies::CreateAll(ctx); }));
  }

  // ── The restrictions fail loudly rather than run the wrong equations.
  {
    const std::vector<std::pair<std::string, std::string>> refused = {
        {"gauge", "newtonian"},
        {"modes", "s,t"},
        {"modes", "s, v"},
        {"ic", "ad,cdi"},
    };
    for (const auto& [key, value] : refused) {
      FileContent fc = Emg();
      fc.set(key, value);
      auto ctx = MakeCtx(fc, pba);
      assert(Throws([&] { ScalarTensorSpecies::CreateAll(ctx); }));
    }
    background curved = pba;
    curved.Omega0_k   = 0.01;
    FileContent fc    = Emg();
    auto ctx          = MakeCtx(fc, curved);
    assert(Throws([&] { ScalarTensorSpecies::CreateAll(ctx); }));
  }

  // ── The allowed settings pass.
  {
    FileContent fc = Emg();
    fc.set("gauge", "synchronous");
    fc.set("modes", "s");
    fc.set("ic", "ad");
    auto ctx = MakeCtx(fc, pba);
    assert(ScalarTensorSpecies::CreateAll(ctx).size() == 1);
  }

  // ── Closure: shoot the potential's constant until the effective density today
  //    vanishes; a given Omega_Lambda_st (the resolved value) switches shooting off.
  {
    Built b = Build(Emg(), pba);
    assert(b.st->GetOmega0() == 0.);
    const auto targets = b.st->GetShootingTargets();
    assert(targets.size() == 1);
    assert(targets[0].unknown_param == "Omega_Lambda_st");
    assert(targets[0].target_value == 0.);

    FileContent fc = Emg();
    fc.set("Omega_Lambda_st", "0.01");
    Built given = Build(fc, pba);
    assert(given.st->GetShootingTargets().empty());
    assert(std::fabs(given.st->V(0.) / (3. * pba.H0 * pba.H0 * 0.01) - 1.) < 1e-15);
  }

  // ── The potential: V = 3 H0^2 Omega_Lambda + lambda M_pl^2 phi^4 / 4, M_pl in 1/Mpc.
  {
    Built b              = Build(Emg(), pba);
    const double M_pl    = 3.8093e56;  // reduced Planck mass in 1/Mpc
    const double V_expct = 3e-105 * M_pl * M_pl * std::pow(0.5, 4) / 4.;
    assert(std::fabs(b.st->V(0.5) / V_expct - 1.) < 1e-3);
    assert(std::fabs(b.st->F(0.5) - (1. + 0.3 * 0.25)) < 1e-15);
  }

  // ── GR limit: xi = 0, no potential, a frozen field. The closure adds nothing,
  //    and GR's H = sqrt(rho_tot) is the H the species solved for.
  {
    FileContent fc;
    fc.set("xi_st", "0");
    fc.set("phi_ini_st", "0.4");
    fc.set("Omega_Lambda_st", "0");
    Built b = Build(fc, pba);
    std::vector<double> y(b.bi_size, 0.), back(b.bg_size, 0.);
    const double a          = 1e-3;
    y[b.st->bi_phi_index()] = 0.4;
    b.st->ComputeBackground(a, y.data(), back.data());
    double rho = 2.5, p = 0.8, rho_r = 2.4, rho_m = 0.1;
    b.st->CloseFriedmann(a, 0., y.data(), back.data(), rho, p, rho_r, rho_m);
    assert(std::fabs(rho - 2.5) < 1e-15 && std::fabs(p - 0.8) < 1e-15);
    assert(rho_r == 2.4 && rho_m == 0.1);
    assert(std::fabs(b.st->Rho(back.data())) < 1e-15);
  }

  // ── A generic state: the H, H-dot and phi-double-dot the closure implies solve
  //    the Jordan-frame equations (M_pl = 1, rho in CLASS units = rho_phys / 3):
  //      (00)  f H^2 + f_dot H = rho_bare
  //      (ij)  -f (2 H_dot + 3 H^2) = 3 p_bare + f_ddot + 2 H f_dot
  //      (KG)  box phi - V_phi + f_phi R / 2 = 0,  R = 6 (H_dot + 2 H^2)
  {
    FileContent fc = Emg();
    fc.set("Omega_Lambda_st", "1.4e6");  // V_Lambda ~ 0.2 in the units of this state
    Built b = Build(fc, pba);
    std::vector<double> y(b.bi_size, 0.), back(b.bg_size, 0.);
    const double a = 0.5, phi = 0.5, phi_prime = 0.3;
    y[b.st->bi_phi_index()]       = phi;
    y[b.st->bi_phi_prime_index()] = phi_prime;
    b.st->ComputeBackground(a, y.data(), back.data());
    const double rho_rest = 1.0, p_rest = 0.2;
    double rho = rho_rest, p = p_rest, rho_r = 0.6, rho_m = 0.4;
    b.st->CloseFriedmann(a, 0., y.data(), back.data(), rho, p, rho_r, rho_m);

    const double H     = std::sqrt(rho);  // GR's line, now exact
    const double H_dot = -(3. * p + 3. * H * H) / 2.;
    const double f = b.st->F(phi), f_phi = 0.6 * phi, f_phiphi = 0.6;
    const double h = 1e-6;
    const double V = b.st->V(phi), V_phi = (b.st->V(phi + h) - b.st->V(phi - h)) / (2. * h);
    const double phi_dot  = phi_prime / a;
    const double phi_ddot = (b.st->PhiPrimePrime(back.data()) - a * H * phi_prime) / (a * a);
    const double f_dot    = f_phi * phi_dot;
    const double f_ddot   = f_phi * phi_ddot + f_phiphi * phi_dot * phi_dot;
    const double rho_bare = rho_rest + (phi_dot * phi_dot / 2. + V) / 3.;
    const double p_bare   = p_rest + (phi_dot * phi_dot / 2. - V) / 3.;

    const double e00 = f * H * H + f_dot * H - rho_bare;
    const double eij = -f * (2. * H_dot + 3. * H * H) - (3. * p_bare + f_ddot + 2. * H * f_dot);
    const double R   = 6. * (H_dot + 2. * H * H);
    const double eKG = -(phi_ddot + 3. * H * phi_dot) - V_phi + f_phi * R / 2.;
    assert(std::fabs(e00) < 1e-12 * rho_bare);
    assert(std::fabs(eij) < 1e-12 * rho_bare);
    assert(std::fabs(eKG) < 1e-6 * std::fabs(V_phi));  // V_phi by finite difference
    // Radiation and matter enter the expansion with G/f: their buckets scale by 1/f.
    assert(std::fabs(rho_r - 0.6 / f) < 1e-15 && std::fabs(rho_m - 0.4 / f) < 1e-15);
  }

  // ── The effective Planck mass must stay positive.
  {
    FileContent fc;
    fc.set("xi_st", "-1");
    fc.set("phi_ini_st", "1.5");
    fc.set("Omega_Lambda_st", "0");
    Built b = Build(fc, pba);
    std::vector<double> y(b.bi_size, 0.), back(b.bg_size, 0.);
    y[b.st->bi_phi_index()] = 1.5;  // f = 1 - 2.25 < 0
    b.st->ComputeBackground(1e-3, y.data(), back.data());
    double rho = 1., p = 0.3, rho_r = 0.9, rho_m = 0.1;
    assert(Throws(
        [&] { b.st->CloseFriedmann(1e-3, 0., y.data(), back.data(), rho, p, rho_r, rho_m); }));
  }

  // ── No expanding solution of (00): a bare density that is not positive must not
  //    come back as H or -H through GR's sqrt(rho_tot). Both roots complex (phi' = 0),
  //    and both roots negative (f_dot > 0).
  {
    for (const double phi_prime : {0., 0.5}) {
      FileContent fc;
      fc.set("xi_st", "0.3");
      fc.set("phi_ini_st", "0.5");
      fc.set("Omega_Lambda_st", "0");
      Built b = Build(fc, pba);
      std::vector<double> y(b.bi_size, 0.), back(b.bg_size, 0.);
      y[b.st->bi_phi_index()]       = 0.5;
      y[b.st->bi_phi_prime_index()] = phi_prime;
      b.st->ComputeBackground(0.5, y.data(), back.data());
      double rho = -1., p = 0., rho_r = 0., rho_m = 0.;
      assert(Throws(
          [&] { b.st->CloseFriedmann(0.5, 0., y.data(), back.data(), rho, p, rho_r, rho_m); }));
    }
  }

  // ── The pre-solve density query assumes GR; this species cannot answer it.
  {
    Built b = Build(Emg(), pba);
    assert(Throws([&] { b.st->BackgroundDensityOverH0Sq(1e-4, pba.H0); }));
  }

  return 0;
}
