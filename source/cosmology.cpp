#include "cosmology.h"
#include "background_module.h"
#include "thermodynamics_module.h"
#include "perturbations_module.h"
#include "primordial_module.h"
#include "nonlinear_module.h"
#include "transfer_module.h"
#include "spectra_module.h"
#include "lensing_module.h"
#include "output_module.h"

InputModulePtr& Cosmology::GetInputModule() {
  return input_module_ptr_;
}

BackgroundModulePtr& Cosmology::GetBackgroundModule() {
  if (!background_module_ptr_) {
    background_module_ptr_ = BackgroundModulePtr(new BackgroundModule(input_module_ptr_));
  }
  return background_module_ptr_;
}

ThermodynamicsModulePtr& Cosmology::GetThermodynamicsModule() {
  if (!thermodynamics_module_ptr_) {
    thermodynamics_module_ptr_ = ThermodynamicsModulePtr(new ThermodynamicsModule(input_module_ptr_, GetBackgroundModule()));
  }
  return thermodynamics_module_ptr_;
}

PerturbationsModulePtr& Cosmology::GetPerturbationsModule() {
  if (!perturbations_module_ptr_) {
    perturbations_module_ptr_ = PerturbationsModulePtr(new PerturbationsModule(input_module_ptr_, GetBackgroundModule(), GetThermodynamicsModule()));
  }
  return perturbations_module_ptr_;
}

PrimordialModulePtr& Cosmology::GetPrimordialModule() {
  if (!primordial_module_ptr_) {
    /** If sigma8 was input, compute local pm and nl module here, compute sigma8, update As and continue*/
    if (input_module_ptr_->primordial_.sigma8 > 0){
      auto pm = PrimordialModulePtr(new PrimordialModule(input_module_ptr_, GetPerturbationsModule()));
      auto nl = NonlinearModule(input_module_ptr_, GetBackgroundModule(), GetPerturbationsModule(), pm);
      double sigma8 = 0;
      if (nl.has_pk_m_ == _TRUE_) {
        sigma8 = nl.sigma8_[nl.index_pk_m_]; 
      }
      else if (nl.has_pk_cb_ == _TRUE_) {
        sigma8 = nl.sigma8_[nl.index_pk_cb_];
      }
      else {
        throw std::invalid_argument("No valid power spectrum found in nonlinear module for calculating sigma8.");
      }
      const_cast<primordial*>(&input_module_ptr_->primordial_)->A_s *= pow(input_module_ptr_->primordial_.sigma8 / sigma8, 2);
    }
    primordial_module_ptr_ = PrimordialModulePtr(new PrimordialModule(input_module_ptr_, GetPerturbationsModule()));
  }
  return primordial_module_ptr_;
}

NonlinearModulePtr& Cosmology::GetNonlinearModule() {
  if (!nonlinear_module_ptr_) {
    nonlinear_module_ptr_ = NonlinearModulePtr(new NonlinearModule(input_module_ptr_, GetBackgroundModule(), GetPerturbationsModule(), GetPrimordialModule()));
  }
  return nonlinear_module_ptr_;
}

TransferModulePtr& Cosmology::GetTransferModule() {
  if (!transfer_module_ptr_) {
    transfer_module_ptr_ = TransferModulePtr(new TransferModule(input_module_ptr_, GetBackgroundModule(), GetThermodynamicsModule(), GetPerturbationsModule(), GetNonlinearModule()));
  }
  return transfer_module_ptr_;
}

SpectraModulePtr& Cosmology::GetSpectraModule() {
  if (!spectra_module_ptr_) {
    spectra_module_ptr_ = SpectraModulePtr(new SpectraModule(input_module_ptr_, GetPerturbationsModule(), GetPrimordialModule(), GetNonlinearModule(), GetTransferModule()));
  }
  return spectra_module_ptr_;
}

LensingModulePtr& Cosmology::GetLensingModule() {
  if (!lensing_module_ptr_) {
    lensing_module_ptr_ = LensingModulePtr(new LensingModule(input_module_ptr_, GetSpectraModule()));
  }
  return lensing_module_ptr_;
}
