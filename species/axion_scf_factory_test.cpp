#include <cassert>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

#include "background.h"
#include "parser.h"
#include "scalar_field.h"
#include "species_build_context.h"

static SpeciesBuildContext MakeCtx(FileContent& fc, const background& pba) {
  SpeciesBuildContext ctx{};
  ctx.pfc = &fc;
  ctx.pba = &pba;
  return ctx;
}

int main() {
  background pba{};
  pba.H0 = 2.2e-4;

  // ── Valid axion config: shooting on m, frozen ICs, params = [m, f, n, Theta].
  {
    FileContent fc;
    fc.set("Omega_scf", "0.05");
    fc.set("scf_potential", "axion");
    fc.set("f_axion", "0.5");
    fc.set("n_axion", "2");
    fc.set("Theta_initial_scf", "2.0");
    auto ctx    = MakeCtx(fc, pba);
    auto result = ScalarFieldSpecies::CreateAll(ctx);
    assert(result.size() == 1);
    auto& scf = static_cast<ScalarFieldSpecies&>(*result[0].species);
    assert(!scf.attractor_ic_scf());
    assert(std::fabs(scf.phi_ini_scf() - 2.0 * 0.5) < 1e-15);
    assert(scf.phi_prime_ini_scf() == 0.);
    assert(scf.scf_tuning_index() == 0);
    assert(scf.scf_parameters().size() == 4);
    assert(scf.scf_parameters()[1] == 0.5);
    assert(scf.scf_parameters()[2] == 2.0);
    assert(scf.scf_parameters()[3] == 2.0);
    assert(scf.scf_parameters()[0] > 0.);  // frozen-field guess seeded
    assert(scf.GetShootingTargets().size() == 1);
    assert(scf.GetShootingTargets()[0].target_name == "Omega_scf");
    assert(scf.GetShootingTargets()[0].unknown_param == "scf_shooting_parameter");
  }

  // ── m_axion, when given, is used as the Newton seed.
  {
    FileContent fc;
    fc.set("Omega_scf", "0.05");
    fc.set("scf_potential", "axion");
    fc.set("m_axion", "3.3e-3");
    fc.set("f_axion", "0.5");
    fc.set("Theta_initial_scf", "2.0");
    auto ctx    = MakeCtx(fc, pba);
    auto result = ScalarFieldSpecies::CreateAll(ctx);
    assert(result.size() == 1);
    auto& scf = static_cast<ScalarFieldSpecies&>(*result[0].species);
    assert(std::fabs(scf.scf_parameters()[0] - 3.3e-3) < 1e-18);
    assert(scf.scf_parameters()[2] == 1.0);  // n_axion defaults to 1
    assert(scf.GetShootingTargets().size() == 1);
  }

  // ── scf_shooting_parameter (resolved value from DoShooting) disables shooting.
  {
    FileContent fc;
    fc.set("Omega_scf", "0.05");
    fc.set("scf_potential", "axion");
    fc.set("f_axion", "0.5");
    fc.set("Theta_initial_scf", "2.0");
    fc.set("scf_shooting_parameter", "4.4e-3");
    auto ctx    = MakeCtx(fc, pba);
    auto result = ScalarFieldSpecies::CreateAll(ctx);
    assert(result.size() == 1);
    auto& scf = static_cast<ScalarFieldSpecies&>(*result[0].species);
    assert(std::fabs(scf.scf_parameters()[0] - 4.4e-3) < 1e-18);
    assert(scf.GetShootingTargets().empty());
  }

  // ── Closure mode: override supplies Omega, shooting still armed.
  {
    FileContent fc;
    fc.set("scf_potential", "axion");
    fc.set("f_axion", "0.5");
    fc.set("Theta_initial_scf", "2.0");
    auto ctx                    = MakeCtx(fc, pba);
    ctx.omega0_closure_override = 0.68;
    auto result                 = ScalarFieldSpecies::CreateAll(ctx);
    assert(result.size() == 1);
    auto& scf = static_cast<ScalarFieldSpecies&>(*result[0].species);
    assert(scf.GetShootingTargets().size() == 1);
    assert(std::fabs(scf.GetShootingTargets()[0].target_value - 0.68) < 1e-15);
  }

  // ── Rejections: attractor requested; scf_parameters given; missing keys; bad ranges.
  auto expect_throw = [&pba](void (*mut)(FileContent&)) {
    FileContent fc;
    fc.set("Omega_scf", "0.05");
    fc.set("scf_potential", "axion");
    fc.set("f_axion", "0.5");
    fc.set("Theta_initial_scf", "2.0");
    mut(fc);
    auto ctx   = MakeCtx(fc, pba);
    bool threw = false;
    try {
      ScalarFieldSpecies::CreateAll(ctx);
    }
    catch (const std::invalid_argument&) {
      threw = true;
    }
    assert(threw);
  };
  expect_throw([](FileContent& fc) { fc.set("attractor_ic_scf", "yes"); });
  expect_throw([](FileContent& fc) { fc.set("scf_parameters", "1.0, 2.0"); });

  // Checks on a parsed numeric value's range are computation rejections
  // (runtime_error) so a sampler can reject the point and keep the chain alive.
  auto expect_throw_computation = [&pba](void (*mut)(FileContent&)) {
    FileContent fc;
    fc.set("Omega_scf", "0.05");
    fc.set("scf_potential", "axion");
    fc.set("f_axion", "0.5");
    fc.set("Theta_initial_scf", "2.0");
    mut(fc);
    auto ctx   = MakeCtx(fc, pba);
    bool threw = false;
    try {
      ScalarFieldSpecies::CreateAll(ctx);
    }
    catch (const std::runtime_error&) {
      threw = true;
    }
    assert(threw);
  };
  expect_throw_computation([](FileContent& fc) { fc.set("f_axion", "-0.5"); });
  expect_throw_computation([](FileContent& fc) { fc.set("n_axion", "0.5"); });
  expect_throw_computation([](FileContent& fc) { fc.set("Theta_initial_scf", "3.5"); });
  expect_throw_computation([](FileContent& fc) { fc.set("Theta_initial_scf", "0"); });

  // Missing f_axion / Theta_initial_scf.
  for (const char* missing : {"f_axion", "Theta_initial_scf"}) {
    FileContent fc;
    fc.set("Omega_scf", "0.05");
    fc.set("scf_potential", "axion");
    if (std::string(missing) != "f_axion")
      fc.set("f_axion", "0.5");
    if (std::string(missing) != "Theta_initial_scf")
      fc.set("Theta_initial_scf", "2.0");
    auto ctx   = MakeCtx(fc, pba);
    bool threw = false;
    try {
      ScalarFieldSpecies::CreateAll(ctx);
    }
    catch (const std::invalid_argument&) {
      threw = true;
    }
    assert(threw);
  }

  // ── Unknown scf_potential value still errors.
  {
    FileContent fc;
    fc.set("Omega_scf", "0.05");
    fc.set("scf_potential", "banana");
    auto ctx   = MakeCtx(fc, pba);
    bool threw = false;
    try {
      ScalarFieldSpecies::CreateAll(ctx);
    }
    catch (const std::invalid_argument&) {
      threw = true;
    }
    assert(threw);
  }

  std::printf("axion scf factory tests passed\n");
  return 0;
}
