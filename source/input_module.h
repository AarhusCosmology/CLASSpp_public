/** @file input_module.h Input parsing and cosmological-species construction. */

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
  static void file_content_from_arguments(int argc, char** argv, FileContent& fc);

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
  static void ShootingResidual(double* x, int x_size, void* pworkspace, double* output);

  void ConstructSpecies();

  /** Parse the coupled-species cluster (CDM, IDM_DR/IDR, IDM_DRMD/IDR_DRMD,
   *  DCDM_DR) into omega_budget_ and coupled_inputs_. Resolves f_idm_dr/f_idm_drmd
   *  CDM subtractions, IDR derivation from T_idr*Omega0_g, and the synchronous-gauge
   *  CDM minimum. Species-specific physics params (T_idr, l_max_idr, f_idm_drmd,
   *  delta_Neff_drmd, z_stop, G_over_aH_drmd, Gamma_dcdm, Omega_ini_dcdm) are
   *  owned by their species and parsed in each species' CreateAll; this
   *  function also stores them in coupled_inputs_ so factories need not re-parse.
   *  Called from ReadContext. */
  void ReadCoupledCluster();

  void ReadContext();                // phase i: inputs needed to build species
  void ReadDerived();                // phase iii: everything else + species-dependent reads
  void SelectPerturbationEvolver();  // phase iv: opt into the explicit evolver
  void WriteParameterFiles();        // read/unread parameter dump (runs after ReadDerived)
  void input_read_precisions();
};

#endif
