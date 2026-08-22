#pragma once
#include <string_view>

#include "../species/base_species.h"
#include "background.h"
#include "constants.h"
#include "perturbations.h"
#include "species_build_context.h"

/** Photons: rho ~ a^{-4}. Boltzmann hierarchy with l_max_g multipoles. */
class PhotonsSpecies : public BaseSpecies {
 public:
  /** Benchmarked against ndf15 across the standard parameter range; see
   *  BaseSpecies::SupportsExplicitPerturbationEvolver(). */
  bool SupportsExplicitPerturbationEvolver() const override {
    return true;
  }
  static constexpr std::string_view kTypeName = "photons";

  explicit PhotonsSpecies(const background& pba)
      : BaseSpecies("Photons", EnergyType::Radiation), pba_(pba) {}

  // ── PerturbLayout ──────────────────────────────────────────────────────────

  struct PerturbLayout : BaseSpecies::PerturbLayout {
    // Scalar / vector / tensor — same struct holds whichever mode this pv is for.
    int idx_delta = -1, idx_theta = -1, idx_shear = -1, idx_l3 = -1;
    int idx_pol0 = -1, idx_pol1 = -1, idx_pol2 = -1, idx_pol3 = -1;
    int l_max     = -1;
    int l_max_pol = -1;
  };

  std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
    return std::make_unique<PerturbLayout>();
  }

  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);

  double GetOmega0() const override {
    return pba_.Omega0_g;
  }

  void RegisterBackgroundIndices(int& index_bg) override {
    index_bg_rho_ = index_bg++;
  }

  void ComputeBackground(double a, const double* /*pvecback_B*/, double* pvecback) override {
    pvecback[index_bg_rho_] = pba_.Omega0_g * pba_.H0 * pba_.H0 / (a * a * a * a);
  }

  double Rho(const double* pvecback) const override {
    return pvecback[index_bg_rho_];
  }
  double P(const double* pvecback) const override {
    return pvecback[index_bg_rho_] / 3.;
  }
  double PPrime(double a,
                double H,
                const double* /*pvecback_B*/,
                const double* pvecback) const override {
    return a * H * (-4. / 3. * pvecback[index_bg_rho_]);
  }

  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override {
    w.Add("(.)rho_g", 0.);
  }
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override {
    w.Add("(.)rho_g", pvecback[index_bg_rho_]);
  }

  // ── Perturbation index registration ────────────────────────────────────────

  void RegisterTransferSourceIndices(int& index_tp, const SourceRequestContext& ctx) override;

  void RegisterPerturbationIndices(BaseSpecies::PerturbLayout& layout,
                                   perturb_vector* pv,
                                   const precision* ppr,
                                   int& index_pt,
                                   const perturb_workspace* ppw,
                                   int gauge) override;

  void RegisterVectorPerturbationIndices(BaseSpecies::PerturbLayout& layout,
                                         perturb_vector* pv,
                                         const precision* ppr,
                                         int& index_pt,
                                         const perturb_workspace* ppw,
                                         int gauge) override;

  void RegisterTensorPerturbationIndices(BaseSpecies::PerturbLayout& layout,
                                         perturb_vector* pv,
                                         const precision* ppr,
                                         int& index_pt,
                                         const perturb_workspace* ppw,
                                         int gauge) override;

  // ── PerturbDerivs ──────────────────────────────────────────────────────────

  void PerturbDerivs(const BaseSpecies::PerturbLayout& layout,
                     double tau,
                     const double* y,
                     double* dy,
                     const perturb_parameters_and_workspace& ppaw) const override;

  void PerturbVectorDerivs(const BaseSpecies::PerturbLayout& layout,
                           double tau,
                           const double* y,
                           double* dy,
                           const perturb_parameters_and_workspace& ppaw) const override;

  void PerturbTensorDerivs(const BaseSpecies::PerturbLayout& layout,
                           double tau,
                           const double* y,
                           double* dy,
                           const perturb_parameters_and_workspace& ppaw) const override;

  // ── Tensor tight coupling ──────────────────────────────────────────────────

  /** The polarisation combination \f$ P^{(2)} \f$ that couples the temperature
   *  and polarisation tensor hierarchies. Shared so that PerturbTensorDerivs,
   *  the source, and the tight-coupling closure cannot drift apart. */
  static double TensorP2(double F0, double shear, double F4, double G0, double G2, double G4) {
    return -1. / _SQRT6_ *
           (1. / 10. * F0 + 2. / 7. * shear + 3. / 70. * F4 - 3. / 5. * G0 + 6. / 7. * G2 -
            3. / 70. * G4);
  }

  /** Leading-order tensor tight-coupling state, in the multipoles the hierarchy
   *  evolves. Everything above l = 0 is higher order in 1/kappa' and stays zero. */
  struct TensorTcaClosure {
    double F0; /**< photon temperature monopole (delta_g) */
    double G0; /**< photon polarisation monopole (pol0_g) */
    double P2; /**< the combination above, evaluated on this state */
  };

  /** \f$ D = d/d\tau (\dot h / \kappa') \f$, the drift the first-order closure
   *  needs, from the quantities a call site already has. */
  static double TensorTcaDrift(double gwdot, double dkappa, double gw_second, double ddkappa) {
    return gw_second / dkappa - gwdot * ddkappa / (dkappa * dkappa);
  }

  /** Solve the tensor hierarchy's stiff balance for (F_0, G_0, P^{(2)}).
   *
   *  \param gwdot   \f$ \dot h \f$, the gravitational-wave velocity.
   *  \param dkappa  \f$ \kappa' \f$, the Thomson scattering rate.
   *  \param drift   \f$ D \f$ from TensorTcaDrift(). Zero -- the default -- gives
   *                 the leading-order closure, for call sites with no time
   *                 derivative available. */
  static TensorTcaClosure TensorTightCoupling(double gwdot, double dkappa, double drift = 0.);

  // ── Source filling and initial conditions ──────────────────────────────────

  void FillSources(const BaseSpecies::PerturbLayout& layout,
                   const double* y,
                   const double* dy,
                   PerturbSourceContext& ctx) const override;

  void ApplyInitialConditions(const BaseSpecies::PerturbLayout& layout,
                              double* y,
                              const PerturbIcContext& ctx) override;

  void PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& layout,
                                     double* y,
                                     const PerturbIcContext& ctx) override;

  // ── MarkUsedInSources ──────────────────────────────────────────────────────

  /** Photon temperature l>=3 and polarization l=1,3+ are not needed in source
   *  evaluation when both rsa and tca are off. Called only in scalar mode
   *  (the dispatch loop in the module is guarded by if (_scalars_)). */
  void MarkUsedInSources(const BaseSpecies::PerturbLayout& layout,
                         const perturb_workspace* ppw,
                         int* used_in_sources) const override;

  /** In tensor mode, photon temperature l=0,2,4 and pol l=0,2,4 are needed;
   *  mark higher multipoles as unused when both rsa and tca are off. */
  void MarkTensorUsedInSources(const BaseSpecies::PerturbLayout& layout,
                               const perturb_workspace* ppw,
                               int* used_in_sources) const override;

  // ── Switch-copy hook ────────────────────────────────────────────────────────

  void CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_layout,
                                     const BaseSpecies::PerturbLayout& new_layout,
                                     const double* old_y,
                                     double* new_y,
                                     const PerturbSwitchContext& ctx) const override;

  void WriteOutputColumns(
      PerturbColumnWriter& writer,
      const PerturbationsModule& mod,
      file_format fmt,
      TransferColumnSection section = TransferColumnSection::all) const override;
  void PrintVariables(PerturbColumnWriter& writer,
                      const BaseSpecies::PerturbLayout* base,
                      double tau,
                      const double* y,
                      const PerturbationsModule& mod,
                      const perturb_workspace* ppw) const override;

  // ── Stress-energy observables ───────────────────────────────────────────────

  StressEnergyContribution StressEnergy(const BaseSpecies::PerturbLayout& layout,
                                        const perturb_vector* pv,
                                        const double* y,
                                        const double* pvecback,
                                        const perturb_workspace* ppw) const override;

  void SetThermodynamicsModule(const ThermodynamicsModule* thm) override {
    thm_ = thm;
  }
  void SetPerturbs(const perturbs* ppt) override {
    ppt_ = ppt;
  }

 private:
  const background& pba_;
  const ThermodynamicsModule* thm_ = nullptr;
  const perturbs* ppt_             = nullptr;

  int index_tp_delta_ = -1;  // #309 transfer-source slot
  int index_tp_theta_ = -1;
};
