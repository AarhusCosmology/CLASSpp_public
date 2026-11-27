#pragma once
#include <algorithm>
#include <vector>

#include "photons.h"

/** Photons in the total angular momentum (TAM) basis: temperature multipoles Theta_l^(m), and
 *  E_l^(m), B_l^(m) polarisation multipoles from l = 2 (Hu & White astro-ph/9702170; Hu, Seljak,
 *  White & Zaldarriaga astro-ph/9709066), truncated with eqs. (9) and (10) of Pitrou, Pereira &
 *  Lesgourgues 2005.12119. Exact for polarisation in curved space, where the optimal basis is an
 *  approximation: errors of order |Omega_k| per cent in scalar EE and 50 |Omega_k| per cent in
 *  tensor EE and BB at low l (2005.12119). Chosen with `hierarchy = tam`.
 *
 *  Cyril Pitrou's work: the equations, coupling coefficients, closure and tight-coupling
 *  relations here are his implementation in class_public v3.4.0, ported to this species.
 *
 *  delta_g, theta_g and shear_g keep the optimal normalisation (PhotonsSpecies); above them
 *  Theta_l = (2l+1)/4 F_l. Tensor modes run two multipoles beyond l_max_g_ten and
 *  l_max_pol_g_ten, as class_public does: Theta_2 .. Theta_{L+2} are as many multipoles as
 *  F_0 .. F_L, so each hierarchy is as long as the optimal basis's, and B is a third one.
 *  Truncation alone, TAM at L+2 then lies between the optimal basis at L+1 and L+2 (tensor
 *  EE at l <= 30, L = 5: 0.6-0.8% against 2.3-2.9%); TAM at L+1 matches it at L. No vector
 *  modes yet.
 *
 *  Design: docs/superpowers/specs/2026-11-27-tam-photon-hierarchy-design.md */
class PhotonsTamSpecies final : public PhotonsSpecies {
 public:
  using PhotonsSpecies::PhotonsSpecies;

  /** Tensor hierarchies start at l = 2: two more multipoles make each as long as the optimal
   *  basis's, which starts at l = 0. */
  static constexpr int kTensorExtraMultipoles = 2;

  int MaxMultipole(const precision* ppr, bool tensors) const override {
    return tensors ? std::max(ppr->l_max_g_ten, ppr->l_max_pol_g_ten) + kTensorExtraMultipoles
                   : std::max(ppr->l_max_g, ppr->l_max_pol_g);
  }

  /** Scalar modes: idx_delta, idx_theta, idx_shear (PhotonsSpecies); idx_l3 is Theta_3, followed
   *  by Theta_4 .. Theta_l_max; idx_e2 is E_2, followed by E_3 .. E_l_max_pol. No B for m = 0.
   *  Tensor modes: idx_t2 is Theta_2^(2), followed by Theta_3 .. Theta_l_max; idx_e2 and idx_b2
   *  are E_2^(2) and B_2^(2), each followed by l = 3 .. l_max_pol. */
  struct PerturbLayout : PhotonsSpecies::PerturbLayout {
    int idx_l3    = -1;
    int idx_t2    = -1;
    int idx_e2    = -1;
    int idx_b2    = -1;
    int l_max     = -1;
    int l_max_pol = -1;
  };

  std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
    return std::make_unique<PerturbLayout>();
  }

  /** The free-streaming coefficients of one k-mode and one m, eq. (5) of 2005.12119:
   *  {}_s kappa_l^m = sqrt((l^2 - m^2)(l^2 - s^2)/l^2) sqrt(q^2 - K l^2), q^2 = k^2 + (1+m) K,
   *  zero where l <= m or l <= s. For m = 0, {}_0 kappa_l^0 = k l s_l. */
  struct Coefficients {
    double q = 0.;
    std::vector<double> kappa0; /**< {}_0 kappa_l^m, l = 0 .. l_max */
    std::vector<double> kappa2; /**< {}_2 kappa_l^m, l = 0 .. l_max */
  };
  static double Kappa(int s, int m, int l, double q2, double K);
  static void FillCoefficients(Coefficients& c, double k, double K, int m, int l_max);

  /** Per-k scratch: the coefficients of the k-mode being integrated, filled on first use. */
  struct Scratch : BaseSpecies::PerturbScratch {
    Coefficients c;
    double k = -1.;
    int m    = -1;
  };
  std::unique_ptr<BaseSpecies::PerturbScratch> CreatePerturbScratch() const override {
    return std::make_unique<Scratch>();
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

  /** The sources need delta_g, theta_g, shear_g and E_2 (scalars), Theta_2 and E_2 (tensors). */
  void MarkUsedInSources(const BaseSpecies::PerturbLayout& base,
                         const perturb_workspace* ppw,
                         int* used_in_sources) const override;
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

 private:
  /** This k-mode's coefficients for m, from the workspace's scratch (filled once per k). */
  const Coefficients& CoefficientsFor(const perturb_workspace* ppw,
                                      double k,
                                      int m,
                                      int l_max) const;
};
