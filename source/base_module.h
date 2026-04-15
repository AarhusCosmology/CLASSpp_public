#ifndef BASE_MODULE_H
#define BASE_MODULE_H

/* class modules */
#include <map>
#include <memory>
#include <string>

#include "../species/base_species.h"
#include "background.h"
#include "common.h"
#include "input_module.h"
#include "lensing.h"
#include "nonlinear.h"
#include "output.h"
#include "perturbations.h"
#include "primordial.h"
#include "spectra.h"
#include "thermodynamics.h"
#include "transfer.h"

class BaseSpecies;
class NCDMSpecies;

class BaseModule {
 public:
  BaseModule(InputModulePtr input_module)
      : ncdm_(input_module->ncdm_), dr_(input_module->dr_),
        all_species_(input_module->all_species_), ppr(&input_module->precision_),
        pba(&input_module->background_), pth(&input_module->thermodynamics_),
        ppt(&input_module->perturbations_), ppm(&input_module->primordial_),
        pnl(&input_module->nonlinear_), ptr(&input_module->transfers_),
        psp(&input_module->spectra_), ple(&input_module->lensing_), pop(&input_module->output_) {
    input_module_     = std::move(input_module);
    error_message_[0] = '\n';
  }
  BaseModule(const BaseModule&) = delete;

  mutable ErrorMsg error_message_;

 public:
  const std::shared_ptr<NonColdDarkMatter> ncdm_;
  const std::shared_ptr<DarkRadiation> dr_;

  /** Const map of all cosmological species, keyed by name.
   *  Use .at("CDM") -- never operator[] -- to preserve const safety. */
  const std::map<std::string, std::unique_ptr<BaseSpecies>>& all_species_;

 protected:
  InputModulePtr input_module_;

  const precision* const ppr;
  const background* const pba;
  const thermo* const pth;
  const perturbs* const ppt;
  const primordial* const ppm;
  const nonlinear* const pnl;
  const transfers* const ptr;
  const spectra* const psp;
  const lensing* const ple;
  const output* const pop;
};

#endif  //BASE_MODULE_H
