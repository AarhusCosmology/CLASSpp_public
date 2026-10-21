// What a species cannot evolve, it says: SupportsNewtonianGauge() and
// SupportsIsocurvatureInitialConditions(). The perturbation module asks every species
// once, before it integrates anything, and refuses the run with the species' name.
//
// Before, the module named the composites it knew to be adiabatic-only or synchronous-only,
// and Type-3 and NEDE parsed the `gauge` string again at construction, each with its own
// spelling rules. The parser's ppt->gauge is now the only reading of that key.
#include <cassert>
#include <cstdio>
#include <stdexcept>
#include <string>

#include "cosmology.h"

namespace {

FileContent Quiet() {
  FileContent fc;
  fc.set("output", "mPk");
  fc.set("P_k_max_1/Mpc", "0.1");
  fc.set("input_verbose", "0");
  fc.set("background_verbose", "0");
  fc.set("thermodynamics_verbose", "0");
  fc.set("perturbations_verbose", "0");
  fc.set("write warnings", "no");
  fc.set("write parameters", "no");
  return fc;
}

/* The message of the severe error (std::invalid_argument) the perturbation module throws,
   or "" if it was built. */
std::string PerturbationsMessage(FileContent fc) {
  try {
    Cosmology cosmology{fc};
    cosmology.GetPerturbationsModule();
  }
  catch (const std::invalid_argument& e) {
    return e.what();
  }
  return "";
}

bool Has(const std::string& haystack, const std::string& needle) {
  return haystack.find(needle) != std::string::npos;
}

FileContent IdmDr() {
  FileContent fc = Quiet();
  fc.set("Omega_idm_dr", "0.01");
  fc.set("xi_idr", "0.5");
  fc.set("a_idm_dr", "1e3");
  return fc;
}

FileContent Drmd() {
  FileContent fc = Quiet();
  fc.set("f_idm_drmd", "0.1");
  fc.set("delta_Neff_drmd", "0.5");
  fc.set("z_stop", "1e4");
  fc.set("G_over_aH_drmd_ini", "1.0");
  return fc;
}

FileContent Dcdm() {
  FileContent fc = Quiet();
  fc.set("Omega_dcdmdr", "0.12");
  fc.set("Gamma_dcdm", "10.0");
  return fc;
}

FileContent Type3() {
  FileContent fc = Quiet();
  fc.set("Omega_fld", "0");
  fc.set("Omega_scf", "0.1");
  fc.set("attractor_ic_scf", "yes");
  fc.set("scf_parameters", "10.0, 0.0, 0.0, 0.0, 100.0, 0.0");
  fc.set("scf_veta", "0.2");
  return fc;
}

FileContent Nede() {
  FileContent fc = Quiet();
  fc.set("f_NEDE", "0.1");
  fc.set("m_NEDE", "316.2278");
  return fc;
}

void ExpectIsocurvatureRefused(FileContent fc, const std::string& species) {
  fc.set("ic", "ad,cdi");
  const std::string message = PerturbationsMessage(fc);
  if (!(Has(message, species) && Has(message, "adiabatic"))) {
    std::printf("%s with cdi: got \"%s\"\n", species.c_str(), message.c_str());
    assert(false);
  }
  std::printf("%s refuses isocurvature initial conditions: ok\n", species.c_str());
}

void ExpectNewtonianRefused(FileContent fc, const std::string& species, const char* gauge) {
  fc.set("gauge", gauge);
  const std::string message = PerturbationsMessage(fc);
  if (!(Has(message, species) && Has(message, "synchronous gauge only"))) {
    std::printf("%s with gauge = %s: got \"%s\"\n", species.c_str(), gauge, message.c_str());
    assert(false);
  }
  std::printf("%s refuses gauge = %s: ok\n", species.c_str(), gauge);
}

/* A switch in tensor modes reaches every species' CopyPerturbationsAcrossSwitch, but only
   scalar workspaces hold species approximations: DRMD must not ask for its regime there. */
void TestDrmdRunsWithTensors() {
  FileContent fc = Drmd();
  fc.set("output", "tCl,mPk");
  fc.set("modes", "s,t");
  fc.set("r", "0.1");
  fc.set("l_max_scalars", "200");
  fc.set("l_max_tensors", "200");
  const std::string message = PerturbationsMessage(fc);
  if (!message.empty()) {
    std::printf("IDM_DRMD_IDR_DRMD with tensor modes: got \"%s\"\n", message.c_str());
    assert(false);
  }
  std::printf("IDM_DRMD_IDR_DRMD with tensor modes: ok\n");
}

}  // namespace

int main() {
  ExpectIsocurvatureRefused(Dcdm(), "DCDM_DR");
  ExpectIsocurvatureRefused(IdmDr(), "IDM_DR_IDR");
  ExpectIsocurvatureRefused(Drmd(), "IDM_DRMD_IDR_DRMD");

  ExpectNewtonianRefused(Drmd(), "IDM_DRMD_IDR_DRMD", "newtonian");
  ExpectNewtonianRefused(Nede(), "NEDE", "newtonian");
  // Every spelling the parser reads as Newtonian; Type-3's own matcher once let "new" through.
  for (const char* gauge : {"newtonian", "Newtonian", "new", "newt"}) {
    ExpectNewtonianRefused(Type3(), "CDM_SCF_Momentum", gauge);
  }

  TestDrmdRunsWithTensors();

  // The synchronous-gauge, adiabatic run of a refusing species is not refused.
  assert(PerturbationsMessage(Type3()).empty());
  std::printf("CDM_SCF_Momentum in synchronous gauge: ok\n");
  return 0;
}
