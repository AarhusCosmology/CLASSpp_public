#ifndef COSMOLOGY_H
#define COSMOLOGY_H

#include <memory>

#include "bessel_table_cache.h"
#include "input_module.h"

class Cosmology {
 public:
  /** bessel_cache: the transfer module's Bessel table, kept across cosmologies by a caller that
      computes many of them (classy). Null gives this cosmology a cache of its own. */
  Cosmology(FileContent& fc, std::shared_ptr<BesselTableCache> bessel_cache = nullptr)
      : input_module_ptr_(InputModulePtr(new InputModule(fc))),
        bessel_cache_(bessel_cache ? std::move(bessel_cache)
                                   : std::make_shared<BesselTableCache>()) {
    InputModule::RejectUnreadParameters(fc);
  }
  Cosmology(std::unique_ptr<InputModule> input_module)
      : input_module_ptr_(std::move(input_module)),
        bessel_cache_(std::make_shared<BesselTableCache>()) {}

  InputModulePtr& GetInputModule();
  BackgroundModulePtr& GetBackgroundModule();
  ThermodynamicsModulePtr& GetThermodynamicsModule();
  PerturbationsModulePtr& GetPerturbationsModule();
  PrimordialModulePtr& GetPrimordialModule();
  NonlinearModulePtr& GetNonlinearModule();
  TransferModulePtr& GetTransferModule();
  SpectraModulePtr& GetSpectraModule();
  LensingModulePtr& GetLensingModule();

 private:
  bool shot_ = false;  // GetInputModule runs DoShooting once, on first use
  InputModulePtr input_module_ptr_;
  BackgroundModulePtr background_module_ptr_;
  ThermodynamicsModulePtr thermodynamics_module_ptr_;
  PerturbationsModulePtr perturbations_module_ptr_;
  PrimordialModulePtr primordial_module_ptr_;
  NonlinearModulePtr nonlinear_module_ptr_;
  TransferModulePtr transfer_module_ptr_;
  SpectraModulePtr spectra_module_ptr_;
  LensingModulePtr lensing_module_ptr_;
  std::shared_ptr<BesselTableCache> bessel_cache_;
};
#endif  //COSMOLOGY_H
