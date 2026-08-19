#pragma once

#include "fluid.h"

/**
 * PPF (Parametrised Post-Friedmann) dark-energy fluid. Same background as
 * FluidSpecies; the perturbations use the Gamma closure (Hu 0808.3125) instead
 * of delta/theta. At most one PPF fluid may exist (PPF is defined relative to
 * the whole universe).
 */
class PpfFluid : public FluidSpecies {
 public:
  struct PerturbLayout : FluidSpecies::PerturbLayout {
    int idx_Gamma = -1;
  };
  std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
    return std::make_unique<PerturbLayout>();
  }

  PpfFluid(const background& pba,
           double omega0_fld,
           equation_of_state fluid_eos,
           double w0_fld,
           double wa_fld,
           double cs2_fld,
           double Omega_EDE,
           double c_gamma_over_c_fld);

  void RegisterPerturbationIndices(BaseSpecies::PerturbLayout& layout,
                                   perturb_vector* pv,
                                   const precision* ppr,
                                   int& index_pt,
                                   const perturb_workspace* ppw,
                                   int gauge) override;
  void PerturbDerivs(const BaseSpecies::PerturbLayout& layout,
                     double tau,
                     const double* y,
                     double* dy,
                     const perturb_parameters_and_workspace& ppaw) const override;
  void FillSources(const BaseSpecies::PerturbLayout& layout,
                   const double* y,
                   const double* dy,
                   PerturbSourceContext& ctx) const override;
  void PrintVariables(PerturbColumnWriter& writer,
                      const BaseSpecies::PerturbLayout* base,
                      double tau,
                      const double* y,
                      const PerturbationsModule& mod,
                      const perturb_workspace* ppw) const override;
  void CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_layout,
                                     const BaseSpecies::PerturbLayout& new_layout,
                                     const double* old_y,
                                     double* new_y,
                                     const PerturbSwitchContext& ctx) const override;

  void ComputePpf(double k,
                  double a,
                  double a_prime_over_a,
                  const precision* ppr,
                  const double* y,
                  perturb_workspace* ppw) const;

 private:
  double c_gamma_over_c_fld_ = 0.4;
};
