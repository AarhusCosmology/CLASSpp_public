/** @file input_module.h Documented includes for input module */

#ifndef __INPUT__
#define __INPUT__

#include <memory>
#include <string>
#include <vector>

#include "../species/base_species.h"
#include "../species/species_build_context.h"
#include "../species/species_collection.h"
#include "background.h"
#include "common.h"
#include "lensing.h"
#include "nonlinear.h"
#include "output.h"
#include "parser.h"
#include "perturbations.h"
#include "primordial.h"
#include "quadrature.h"
#include "spectra.h"
#include "thermodynamics.h"
#include "transfer.h"

class InputModule {
 public:
  InputModule(FileContent& fc);
  static int file_content_from_arguments(int argc, char** argv, FileContent& fc);

  /** Resolve any per-species / theta_s shooting targets by root-finding, returning a fully
   *  resolved module. No targets (or already inside a shooting context) → returns the input
   *  unchanged. Called lazily by Cosmology::GetInputModule. */
  static InputModulePtr DoShooting(InputModulePtr input_module);

  FileContent& file_content_;
  precision precision_;
  background background_;
  thermo thermodynamics_;
  perturbs perturbations_;
  transfers transfers_;
  primordial primordial_;
  spectra spectra_;
  nonlinear nonlinear_;
  lensing lensing_;
  output output_;
  /** All cosmological species, constructed at end of InputModule ctor. */
  SpeciesCollection all_species_;
  /** Pre-resolved Omega budget for coupled species (CDM, IDM_DR/IDR,
   *  IDM_DRMD/IDR_DRMD, DCDM_DR). Populated by ReadCoupledCluster during
   *  ReadContext and passed to ConstructSpecies via SpeciesBuildContext. */
  SpeciesOmegaBudget omega_budget_;
  /** Raw intermediates from phase-i parsing that coupled factories need for
   *  physics construction. Populated by ReadCoupledCluster alongside
   *  omega_budget_ and threaded through SpeciesBuildContext. */
  CoupledClusterInputs coupled_inputs_;

 private:
  // Hook-based shooting (used by DoShooting). theta_s is identified by
  // target_name == "100*theta_s" / unknown_param == "h"; all targets are scalar.
  struct ShootingWorkspace {
    explicit ShootingWorkspace(FileContent& fc_ref) : fc(fc_ref) {}
    FileContent& fc;
    std::vector<ShootingTarget> targets;  // theta_s first (if present), then species in lex order
    // Parallel to targets: the all_species_ key owning each target ("" for the module-level
    // theta_s). Lets ShootingResidual route each slot to its species' ComputeShootingResidual
    // with the authoritative target (so the species never re-derives it from the file content).
    std::vector<std::string> target_species_keys;
  };
  static int ShootingResidual(double* x, int x_size, void* pworkspace, double* output);

  void ConstructSpecies();

  /** Parse the coupled-species cluster (CDM, IDM_DR/IDR, IDM_DRMD/IDR_DRMD,
   *  DCDM_DR) into omega_budget_ and coupled_inputs_. Resolves f_idm_dr/f_idm_drmd
   *  CDM subtractions, IDR derivation from T_idr*Omega0_g, and the synchronous-gauge
   *  CDM minimum. Species-specific physics params (T_idr, l_max_idr, f_idm_drmd,
   *  delta_Neff_drmd, z_stop, G_over_aH_drmd, Gamma_dcdm, Omega_ini_dcdm) are
   *  owned by their species and parsed in each species' CreateAll; this
   *  function also stores them in coupled_inputs_ so factories need not re-parse.
   *  Called from ReadContext. */
  int ReadCoupledCluster();

  void ReadContext();          // phase i: inputs needed to build species
  void ReadDerived();          // phase iii: everything else + species-dependent reads
  void WriteParameterFiles();  // read/unread parameter dump (runs after ReadDerived)
  int input_read_precisions();
};

/* macro for reading parameter values with routines from the parser */
#define class_read_double(name, destination)        \
  do {                                              \
    parser_read_double(pfc, name, &param1, &flag1); \
    if (flag1 == _TRUE_)                            \
      destination = param1;                         \
  } while (0);

#define class_read_int(name, destination)       \
  do {                                          \
    parser_read_int(pfc, name, &int1, &flag1);  \
    if (flag1 == _TRUE_)                        \
      destination = (typeof(destination)) int1; \
  } while (0);

/* `destination` must be a std::string (the assignment below relies on it); legacy
 * char[] fields are not supported by this macro. `string1` is a std::string local. */
#define class_read_string(name, destination)        \
  do {                                              \
    parser_read_string(pfc, name, string1, &flag1); \
    if (flag1 == _TRUE_)                            \
      destination = string1;                        \
  } while (0);

#define class_read_double_one_of_two(name1, name2, destination)               \
  do {                                                                        \
    parser_read_double(pfc, name1, &param1, &flag1);                          \
    parser_read_double(pfc, name2, &param2, &flag2);                          \
    class_test((flag1 == _TRUE_) && (flag2 == _TRUE_),                        \
               "In input file, you can only enter one of %s, %s, choose one", \
               name1,                                                         \
               name2);                                                        \
    if (flag1 == _TRUE_)                                                      \
      destination = param1;                                                   \
    if (flag2 == _TRUE_)                                                      \
      destination = param2;                                                   \
  } while (0);

#define class_at_least_two_of_three(a, b, c)                              \
  ((a == _TRUE_) && (b == _TRUE_)) || ((a == _TRUE_) && (c == _TRUE_)) || \
      ((b == _TRUE_) && (c == _TRUE_))

#define class_at_least_two_of_four(a, b, c, d)                                \
  ((a == _TRUE_) && (b == _TRUE_)) || ((a == _TRUE_) && (c == _TRUE_)) ||     \
      ((a == _TRUE_) && (d == _TRUE_)) || ((b == _TRUE_) && (c == _TRUE_)) || \
      ((b == _TRUE_) && (d == _TRUE_)) || ((c == _TRUE_) && (d == _TRUE_))

#define class_any_nonzero_four(a, b, c, d) ((a) != 0. || (b) != 0. || (c) != 0. || (d) != 0.)

#define class_all_nonzero_four(a, b, c, d) ((a) != 0. && (b) != 0. && (c) != 0. && (d) != 0.)

#define class_none_of_three(a, b, c) (a == _FALSE_) && (b == _FALSE_) && (c == _FALSE_)

#endif
