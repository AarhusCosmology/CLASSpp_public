// Multi-instance scalar fields (#346): every "<N>.type = scalar_field" builds one
// independent field under the key N, beside the legacy *_scf field ("ScalarField"),
// which alone may close the budget. A field's output columns carry N (phi_N, d_N), so
// a name that repeats another column's title (b: (.)rho_b) is rejected where the
// titles are made.
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "background.h"
#include "background_module.h"
#include "common.h"
#include "cosmology.h"
#include "input_module.h"
#include "parser.h"
#include "perturbations_module.h"
#include "scalar_field.h"
#include "species_build_context.h"

namespace {

SpeciesBuildContext MakeCtx(FileContent& fc, const background& pba) {
  SpeciesBuildContext ctx{};
  ctx.pfc = &fc;
  ctx.pba = &pba;
  return ctx;
}

const ScalarFieldSpecies& Field(const std::vector<Named>& built, const std::string& key) {
  for (const auto& e : built) {
    if (e.key == key) {
      return static_cast<const ScalarFieldSpecies&>(*e.species);
    }
  }
  throw std::logic_error("no field '" + key + "'");
}

FileContent Quiet() {
  FileContent fc;
  fc.set("output", "");
  fc.set("input_verbose", "0");
  fc.set("background_verbose", "0");
  fc.set("write warnings", "no");
  fc.set("write parameters", "no");
  return fc;
}

template <class E, class F>
bool Throws(const F& f) {
  try {
    f();
  }
  catch (const E&) {
    return true;
  }
  return false;
}

void test_factory(const background& pba) {
  // ── Two dot instances beside the legacy field: three fields, three keys ──
  {
    FileContent fc;
    fc.set("Omega_scf", "0.01");
    fc.set("scf_parameters", "10, 0, 0, 0");
    fc.set("q.type", "scalar_field");
    fc.set("q.Omega", "0.02");
    fc.set("q.parameters", "12, 0, 0, 0");
    fc.set("ax.type", "scalar_field");
    fc.set("ax.Omega", "0.03");
    fc.set("ax.potential", "axion");
    fc.set("ax.f", "0.5");
    fc.set("ax.Theta_initial", "2.0");
    auto ctx         = MakeCtx(fc, pba);
    const auto built = ScalarFieldSpecies::CreateAll(ctx);
    assert(built.size() == 3);

    const auto& legacy = Field(built, "ScalarField");
    assert(legacy.name() == "ScalarField");
    assert(legacy.column_tag() == "scf");
    assert(legacy.GetOmega0() == 0.01);
    assert(legacy.GetShootingTargets().size() == 1);
    assert(legacy.GetShootingTargets()[0].unknown_param == "scf_shooting_parameter");

    const auto& q = Field(built, "q");
    assert(q.name() == "q");
    assert(q.column_tag() == "q");
    assert(q.GetOmega0() == 0.02);
    assert(q.attractor_ic_scf());
    assert(q.GetShootingTargets().size() == 1);
    assert(q.GetShootingTargets()[0].target_name == "q.Omega");
    assert(q.GetShootingTargets()[0].unknown_param == "q.shooting_parameter");
    assert(q.GetShootingTargets()[0].target_value == 0.02);

    const auto& ax = Field(built, "ax");
    assert(!ax.attractor_ic_scf());
    assert(std::fabs(ax.phi_ini_scf() - 2.0 * 0.5) < 1e-15);
    assert(ax.scf_parameters().size() == 4);
    assert(ax.scf_parameters()[1] == 0.5);
    assert(ax.scf_parameters()[2] == 1.0);  // n defaults to 1
    assert(ax.GetShootingTargets()[0].unknown_param == "ax.shooting_parameter");

    // Every dot key the fields use was read.
    for (const auto& key : fc.unread_parameters()) {
      assert(key.rfind("q.", 0) != 0 && key.rfind("ax.", 0) != 0);
    }
  }

  // ── A resolved shooting_parameter fills the tuning slot and disarms shooting ──
  {
    FileContent fc;
    fc.set("q.type", "scalar_field");
    fc.set("q.Omega", "0.02");
    fc.set("q.parameters", "12, 0, 0, 0");
    fc.set("q.tuning_index", "0");
    fc.set("q.shooting_parameter", "7.5");
    auto ctx         = MakeCtx(fc, pba);
    const auto built = ScalarFieldSpecies::CreateAll(ctx);
    assert(built.size() == 1);
    const auto& q = Field(built, "q");
    assert(q.scf_parameters()[0] == 7.5);
    assert(q.GetShootingTargets().empty());
  }

  // ── Non-attractor ICs come from the last two parameters, per instance ──
  {
    FileContent fc;
    fc.set("q.type", "scalar_field");
    fc.set("q.Omega", "0.02");
    fc.set("q.parameters", "12, 0, 0, 0, 0.3, 0.004");
    fc.set("q.attractor_ic", "no");
    auto ctx         = MakeCtx(fc, pba);
    const auto built = ScalarFieldSpecies::CreateAll(ctx);
    const auto& q    = Field(built, "q");
    assert(!q.attractor_ic_scf());
    assert(q.phi_ini_scf() == 0.3);
    assert(q.phi_prime_ini_scf() == 0.004);
  }

  // ── Pass 2 (closure override) builds the legacy field alone ──
  {
    FileContent fc;
    fc.set("Omega_scf", "-1");
    fc.set("scf_parameters", "10, 0, 0, 0");
    fc.set("q.type", "scalar_field");
    fc.set("q.Omega", "0.02");
    fc.set("q.parameters", "12, 0, 0, 0");
    auto ctx                    = MakeCtx(fc, pba);
    ctx.omega0_closure_override = 0.68;
    const auto built            = ScalarFieldSpecies::CreateAll(ctx);
    assert(built.size() == 1);
    assert(built[0].key == "ScalarField");
    assert(built[0].species->GetOmega0() == 0.68);
  }

  // ── Pass 1 with a legacy closure field builds the dot instances only ──
  {
    FileContent fc;
    fc.set("Omega_scf", "-1");
    fc.set("q.type", "scalar_field");
    fc.set("q.Omega", "0.02");
    fc.set("q.parameters", "12, 0, 0, 0");
    auto ctx         = MakeCtx(fc, pba);
    const auto built = ScalarFieldSpecies::CreateAll(ctx);
    assert(built.size() == 1);
    assert(built[0].key == "q");
  }

  // ── Omega = 0 builds no field, as for the legacy key ──
  {
    FileContent fc;
    fc.set("q.type", "scalar_field");
    fc.set("q.Omega", "0");
    auto ctx = MakeCtx(fc, pba);
    assert(ScalarFieldSpecies::CreateAll(ctx).empty());
  }

  // ── Rejections ──
  auto build = [&pba](void (*mut)(FileContent&)) {
    return [&pba, mut]() {
      FileContent fc;
      fc.set("q.type", "scalar_field");
      fc.set("q.Omega", "0.02");
      fc.set("q.parameters", "12, 0, 0, 0");
      mut(fc);
      auto ctx = MakeCtx(fc, pba);
      ScalarFieldSpecies::CreateAll(ctx);
    };
  };
  // Structural: missing Omega, reserved name, unknown potential, bad tuning index.
  assert(Throws<std::invalid_argument>([&pba]() {
    FileContent fc;
    fc.set("q.type", "scalar_field");
    fc.set("q.parameters", "12, 0, 0, 0");
    auto ctx = MakeCtx(fc, pba);
    ScalarFieldSpecies::CreateAll(ctx);
  }));
  assert(Throws<std::invalid_argument>(build([](FileContent& fc) {
    fc.set("ScalarField.type", "scalar_field");
    fc.set("ScalarField.Omega", "0.01");
  })));
  assert(Throws<std::invalid_argument>(build([](FileContent& fc) {
    fc.set("scf.type", "scalar_field");
    fc.set("scf.Omega", "0.01");
  })));
  assert(Throws<std::invalid_argument>(
      build([](FileContent& fc) { fc.set("q.potential", "cosine"); })));
  assert(
      Throws<std::invalid_argument>(build([](FileContent& fc) { fc.set("q.tuning_index", "4"); })));
  // A value a sampler may vary: a computation rejection, not an abort.
  assert(Throws<std::runtime_error>(build([](FileContent& fc) { fc.set("q.Omega", "-0.1"); })));
}

void test_input_module() {
  // ── Two independent fields beside a Lambda that closes the budget ──
  {
    FileContent fc = Quiet();
    fc.set("q1.type", "scalar_field");
    fc.set("q1.Omega", "0.01");
    fc.set("q1.parameters", "10, 0, 0, 0");
    fc.set("q2.type", "scalar_field");
    fc.set("q2.Omega", "0.02");
    fc.set("q2.parameters", "12, 0, 0, 0");
    InputModule input(fc);
    assert(input.all_species_.contains("q1"));
    assert(input.all_species_.contains("q2"));
    assert(!input.all_species_.contains("ScalarField"));
  }

  // ── The legacy field closes the budget; a dot field sits beside it and is
  //    counted in the budget the closure fills ──
  {
    FileContent fc = Quiet();
    fc.set("Omega_Lambda", "0");
    fc.set("Omega_fld", "0");
    fc.set("Omega_scf", "-1");
    fc.set("scf_parameters", "10, 0, 0, 0");
    fc.set("q.type", "scalar_field");
    fc.set("q.Omega", "0.1");
    fc.set("q.parameters", "12, 0, 0, 0");
    InputModule input(fc);
    assert(input.all_species_.contains("q"));
    assert(input.all_species_.contains("ScalarField"));
    double omega_sum = 0.;
    for (const auto& [key, sp] : input.all_species_) {
      omega_sum += sp->GetOmega0();
    }
    assert(std::fabs(omega_sum + input.background_.Omega0_k - 1.) < 1e-12);
  }
}

bool Has(const std::string& haystack, const std::string& needle) {
  return haystack.find(needle) != std::string::npos;
}

FileContent WithField(FileContent fc, const std::string& name) {
  fc.set(name + ".type", "scalar_field");
  fc.set(name + ".Omega", "0.02");
  fc.set(name + ".parameters", "10, 0, 0, 0");
  return fc;
}

std::string BackgroundTitles(FileContent fc, std::vector<std::string>* notes) {
  Cosmology cosmo(fc);
  std::string titles;
  cosmo.GetBackgroundModule()->background_output_titles(titles);
  if (notes != nullptr) {
    *notes = cosmo.GetBackgroundModule()->background_output_notes();
  }
  return titles;
}

void test_column_titles() {
  assert(first_duplicate_title("").empty());
  assert(first_duplicate_title("z\ta\tb\t").empty());
  assert(first_duplicate_title("z\ta\tb\ta\t") == "a");

  // ── Two fields write their own columns, and the units note once ──
  {
    std::vector<std::string> notes;
    const std::string titles = BackgroundTitles(WithField(WithField(Quiet(), "q1"), "q2"), &notes);
    assert(Has(titles, "\tphi_q1\t") && Has(titles, "\tphi_q2\t"));
    assert(notes.size() == 1);
    assert(Has(notes[0], "reduced Planck mass"));
  }
  // ── No scalar field, no note ──
  {
    std::vector<std::string> notes;
    BackgroundTitles(Quiet(), &notes);
    assert(notes.empty());
  }
  // ── b repeats the baryons' (.)rho_b, tot the module's (.)rho_tot ──
  for (const char* name : {"b", "tot"}) {
    const FileContent fc = WithField(Quiet(), name);
    assert(Throws<std::invalid_argument>([&fc]() { BackgroundTitles(fc, nullptr); }));
  }
  // ── The k_output_values columns: cdm repeats CDM's delta_cdm ──
  FileContent k_out = Quiet();
  k_out.set("output", "mTk");
  k_out.set("P_k_max_1/Mpc", "1e-5");  // one k mode beside the k_output one
  k_out.set("k_output_values", "0.01");
  k_out.set("perturbations_verbose", "0");
  {
    FileContent fc = WithField(k_out, "cdm");
    assert(Throws<std::invalid_argument>([&fc]() {
      Cosmology cosmo(fc);
      cosmo.GetPerturbationsModule();
    }));
  }
  // ── The transfer-function columns: tot repeats d_tot (its delta_tot is no clash) ──
  {
    FileContent fc = WithField(k_out, "tot");
    Cosmology cosmo(fc);
    const auto& pt = cosmo.GetPerturbationsModule();
    std::string titles;
    assert(Throws<std::invalid_argument>(
        [&pt, &titles]() { pt->perturb_output_titles(file_format::class_format, titles); }));
  }
}

}  // namespace

int main() {
  background pba{};
  pba.H0 = 2.2e-4;
  test_factory(pba);
  test_input_module();
  test_column_titles();
  return 0;
}
