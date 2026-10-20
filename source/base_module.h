#ifndef BASE_MODULE_H
#define BASE_MODULE_H

/* class modules */
#include <memory>
#include <string>

#include "../species/base_species.h"
#include "../species/species_collection.h"
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
      // Declaration order: all_species_ reads the argument before input_module_ takes it;
      // everything after reads input_module_.
      : all_species_(input_module->all_species_), input_module_(std::move(input_module)),
        ppr(&input_module_->precision_), pba(&input_module_->background_),
        pth(&input_module_->thermodynamics_), ppt(&input_module_->perturbations_),
        ppm(&input_module_->primordial_), pnl(&input_module_->nonlinear_),
        ptr(&input_module_->transfers_), psp(&input_module_->spectra_),
        ple(&input_module_->lensing_), pop(&input_module_->output_) {}
  BaseModule(const BaseModule&) = delete;

  /** Const reference to the species owned by InputModule. Use .at("CDM")
   *  (throws if absent) or .find("CDM") (returns pointer-or-nullptr). */
  const SpeciesCollection& all_species_;

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
