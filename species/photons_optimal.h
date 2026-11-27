#pragma once
#include <algorithm>

#include "photons.h"

/** Photons in the optimal basis of Tram & Lesgourgues (1305.3261): temperature multipoles
 *  F_l and polarisation multipoles G_l, from l = 0, in every mode. Exact in flat space, and for
 *  temperature also in curved space; for polarisation in curved space an approximation (see
 *  PhotonsSpecies). The default (`hierarchy = optimal`). */
class PhotonsOptimalSpecies final : public PhotonsSpecies {
 public:
  using PhotonsSpecies::PhotonsSpecies;

  int MaxMultipole(const precision* ppr, bool tensors) const override {
    return tensors ? std::max(ppr->l_max_g_ten, ppr->l_max_pol_g_ten)
                   : std::max(ppr->l_max_g, ppr->l_max_pol_g);
  }

  /** Scalar, vector and tensor modes alike: idx_delta, idx_theta and idx_shear hold F_0,
   *  F_1 (as theta) and F_2 (as shear); idx_l3 is F_3, followed by F_4 .. F_l_max; idx_pol0
   *  .. idx_pol3 are G_0 .. G_3, followed by G_4 .. G_l_max_pol. */
  struct PerturbLayout : PhotonsSpecies::PerturbLayout {
    int idx_l3   = -1;
    int idx_pol0 = -1, idx_pol1 = -1, idx_pol2 = -1, idx_pol3 = -1;
    int l_max     = -1;
    int l_max_pol = -1;
  };

  std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
    return std::make_unique<PerturbLayout>();
  }

  /** The polarisation combination \f$ P^{(2)} \f$ that couples the temperature
   *  and polarisation tensor hierarchies. Shared so that PerturbTensorDerivs,
   *  the source, and the tight-coupling closure cannot drift apart. */
  static double TensorP2(double F0, double shear, double F4, double G0, double G2, double G4) {
    return -1. / _SQRT6_ *
           (1. / 10. * F0 + 2. / 7. * shear + 3. / 70. * F4 - 3. / 5. * G0 + 6. / 7. * G2 -
            3. / 70. * G4);
  }

  void RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                   perturb_vector* pv,
                                   const precision* ppr,
                                   int& index_pt,
                                   const perturb_workspace* ppw,
                                   int gauge) override;
  void RegisterVectorPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                         perturb_vector* pv,
                                         const precision* ppr,
                                         int& index_pt,
                                         const perturb_workspace* ppw,
                                         int gauge) override;
  void RegisterTensorPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                         perturb_vector* pv,
                                         const precision* ppr,
                                         int& index_pt,
                                         const perturb_workspace* ppw,
                                         int gauge) override;

  void PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                     double tau,
                     const double* y,
                     double* dy,
                     const perturb_parameters_and_workspace& ppaw) const override;
  void PerturbVectorDerivs(const BaseSpecies::PerturbLayout& base,
                           double tau,
                           const double* y,
                           double* dy,
                           const perturb_parameters_and_workspace& ppaw) const override;
  void PerturbTensorDerivs(const BaseSpecies::PerturbLayout& base,
                           double tau,
                           const double* y,
                           double* dy,
                           const perturb_parameters_and_workspace& ppaw) const override;

  /** Photon temperature l>=3 and polarization l=1,3+ are not needed in source
   *  evaluation when both rsa and tca are off. Called only in scalar mode
   *  (the dispatch loop in the module is guarded by if (_scalars_)). */
  void MarkUsedInSources(const BaseSpecies::PerturbLayout& base,
                         const perturb_workspace* ppw,
                         int* used_in_sources) const override;

  /** In tensor mode, photon temperature l=0,2,4 and pol l=0,2,4 are needed;
   *  mark higher multipoles as unused when both rsa and tca are off. */
  void MarkTensorUsedInSources(const BaseSpecies::PerturbLayout& base,
                               const perturb_workspace* ppw,
                               int* used_in_sources) const override;

  void CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_base,
                                     const BaseSpecies::PerturbLayout& new_base,
                                     const double* old_y,
                                     double* new_y,
                                     const PerturbSwitchContext& ctx) const override;

  void StartScalarHierarchy(const PhotonsSpecies::PerturbLayout& base,
                            double* y,
                            double k,
                            double dkappa,
                            const perturb_workspace* ppw) const override;
  double ScalarPolarisationSource(const PhotonsSpecies::PerturbLayout& base,
                                  const double* y,
                                  const perturb_workspace* ppw) const override;
  void StartTensorHierarchy(const PhotonsSpecies::PerturbLayout& base,
                            double* y,
                            const TensorTcaClosure& tca) const override;
  double TensorPolarisationSource(const PhotonsSpecies::PerturbLayout& base,
                                  const double* y) const override;
  double TensorGwSource(const PhotonsSpecies::PerturbLayout& base,
                        const double* y,
                        double a2,
                        double rho_g) const override;
  void AddTensorPerturbationOutput(OutputTable& row,
                                   const PhotonsSpecies::PerturbLayout& base,
                                   const double* y,
                                   const TensorTcaClosure* tca) const override;

 protected:
  void AddPolarisationOutput(OutputTable& row,
                             const PhotonsSpecies::PerturbLayout& base,
                             const double* y,
                             double k,
                             double dkappa,
                             const perturb_workspace* ppw) const override;
};
