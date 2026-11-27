#pragma once
#include <algorithm>
#include <string_view>

#include "../species/base_species.h"
#include "background.h"
#include "constants.h"
#include "perturbations.h"
#include "species_build_context.h"

/** Photons: rho ~ a^{-4}, and a Boltzmann hierarchy for temperature and polarisation.
 *
 *  The hierarchy can be written in two angular bases, each a leaf of this class:
 *  PhotonsOptimalSpecies (Tram & Lesgourgues 1305.3261, the default) and the total angular
 *  momentum basis (Hu, Seljak, White & Zaldarriaga astro-ph/9709066, exact for polarisation
 *  in curved space; Cyril Pitrou's implementation in class_public v3.4.0). The input key
 *  `hierarchy` picks one, and CreateAll builds exactly that one.
 *
 *  Both bases keep delta_g, theta_g and shear_g in the same normalisation. Everything that
 *  reads only those, and the background, lives here. Every hook whose answer depends on the
 *  basis is pure here (`override = 0`), so that a leaf which forgets one does not compile,
 *  instead of running the other basis's code on its own layout.
 *
 *  Design: docs/superpowers/specs/2026-11-27-tam-photon-hierarchy-design.md */
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

  /** What both bases share. In scalar modes these are delta_g, theta_g and shear_g, at -1 when
   *  not evolved (radiation streaming; shear_g also under tight coupling). Each leaf derives
   *  its own layout from this one. */
  struct PerturbLayout : BaseSpecies::PerturbLayout {
    int idx_delta = -1, idx_theta = -1, idx_shear = -1;
  };

  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);

  // ── Background ─────────────────────────────────────────────────────────────

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

  void AddBackgroundOutput(OutputTable& row, const double* pvecback) const override {
    row.Set("(.)rho_g", pvecback[index_bg_rho_]);
  }

  // ── Tensor tight coupling (the physical state, the same in both bases) ─────

  /** Leading-order tensor tight-coupling state, written as the optimal basis's monopoles: the
   *  other multipoles of that basis are higher order in 1/kappa' and vanish. A leaf in another
   *  basis maps this state onto its own variables (StartTensorHierarchy). */
  struct TensorTcaClosure {
    double F0; /**< photon temperature monopole F_0^(2) */
    double G0; /**< photon polarisation monopole G_0^(2) */
    double P2; /**< the polarisation source P^(2) on this state */
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

  // ── Basis-dependent hooks of BaseSpecies: each leaf provides them ──────────

  int MaxMultipole(const precision* ppr, bool tensors) const override              = 0;
  std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override = 0;

  void RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                   perturb_vector* pv,
                                   const precision* ppr,
                                   int& index_pt,
                                   const perturb_workspace* ppw,
                                   int gauge) override       = 0;
  void RegisterVectorPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                         perturb_vector* pv,
                                         const precision* ppr,
                                         int& index_pt,
                                         const perturb_workspace* ppw,
                                         int gauge) override = 0;
  void RegisterTensorPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                         perturb_vector* pv,
                                         const precision* ppr,
                                         int& index_pt,
                                         const perturb_workspace* ppw,
                                         int gauge) override = 0;

  void PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                     double tau,
                     const double* y,
                     double* dy,
                     const perturb_parameters_and_workspace& ppaw) const override       = 0;
  void PerturbVectorDerivs(const BaseSpecies::PerturbLayout& base,
                           double tau,
                           const double* y,
                           double* dy,
                           const perturb_parameters_and_workspace& ppaw) const override = 0;
  void PerturbTensorDerivs(const BaseSpecies::PerturbLayout& base,
                           double tau,
                           const double* y,
                           double* dy,
                           const perturb_parameters_and_workspace& ppaw) const override = 0;

  void MarkUsedInSources(const BaseSpecies::PerturbLayout& base,
                         const perturb_workspace* ppw,
                         int* used_in_sources) const override       = 0;
  void MarkTensorUsedInSources(const BaseSpecies::PerturbLayout& base,
                               const perturb_workspace* ppw,
                               int* used_in_sources) const override = 0;

  void CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_base,
                                     const BaseSpecies::PerturbLayout& new_base,
                                     const double* old_y,
                                     double* new_y,
                                     const PerturbSwitchContext& ctx) const override = 0;

  // ── Photon-only hooks the perturbations module calls ───────────────────────

  /** Scalar modes, tight coupling ending: set every evolved multipole beyond theta_g from the
   *  tight-coupling solution (shear_g = ppw->tca_shear_g). */
  virtual void StartScalarHierarchy(const PerturbLayout& layout,
                                    double* y,
                                    double k,
                                    double dkappa,
                                    const perturb_workspace* ppw) const = 0;

  /** Scalar modes: the polarisation source \f$ P^{(0)} \f$ from the evolved hierarchy (tight
   *  coupling and radiation streaming off). */
  virtual double ScalarPolarisationSource(const PerturbLayout& layout,
                                          const double* y,
                                          const perturb_workspace* ppw) const = 0;

  /** Tensor modes, tight coupling ending: start the hierarchy on the closure's state. */
  virtual void StartTensorHierarchy(const PerturbLayout& layout,
                                    double* y,
                                    const TensorTcaClosure& tca) const = 0;

  /** Tensor modes: \f$ P^{(2)} \f$ from the evolved hierarchy. */
  virtual double TensorPolarisationSource(const PerturbLayout& layout, const double* y) const = 0;

  /** Tensor modes: the photons' contribution to the gravitational-wave source (their
   *  anisotropic stress, times a^2) from the evolved hierarchy. */
  virtual double TensorGwSource(const PerturbLayout& layout,
                                const double* y,
                                double a2,
                                double rho_g) const = 0;

  /** Tensor modes: the photon columns of the perturbation output. From the closure's state when
   *  `tca` is given (tight coupling on), from the hierarchy when it is evolved, zero otherwise
   *  (radiation streaming). */
  virtual void AddTensorPerturbationOutput(OutputTable& row,
                                           const PerturbLayout& layout,
                                           const double* y,
                                           const TensorTcaClosure* tca) const = 0;

  // ── Shared: read only delta_g, theta_g, shear_g ────────────────────────────

  void RegisterTransferSourceIndices(int& index_tp, const SourceRequestContext& ctx) override;

  void FillSources(const BaseSpecies::PerturbLayout& base,
                   const double* y,
                   const double* dy,
                   PerturbSourceContext& ctx) const override;

  void ApplyInitialConditions(const BaseSpecies::PerturbLayout& base,
                              double* y,
                              const PerturbIcContext& ctx) override;

  void PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                     double* y,
                                     const PerturbIcContext& ctx) override;

  StressEnergyContribution StressEnergy(const BaseSpecies::PerturbLayout& base,
                                        const perturb_vector* pv,
                                        const double* y,
                                        const double* pvecback,
                                        const perturb_workspace* ppw) const override;

  void AddTransferOutput(OutputTable& row,
                         const double* tk,
                         const PerturbationsModule& mod,
                         file_format fmt,
                         TransferColumnSection section = TransferColumnSection::all) const override;

  /** delta_g, theta_g and shear_g (in the Newtonian gauge), then the leaf's polarisation
   *  columns. */
  void AddPerturbationOutput(OutputTable& row,
                             const BaseSpecies::PerturbLayout& base,
                             double tau,
                             const double* y,
                             const PerturbationsModule& mod,
                             const perturb_workspace* ppw) const override;

  void SetThermodynamicsModule(const ThermodynamicsModule* thm) override {
    thm_ = thm;
  }
  void SetPerturbs(const perturbs* ppt) override {
    ppt_ = ppt;
  }

 protected:
  /** Scalar modes, the equations both bases share: delta_g' always, and theta_g' (which
   *  reads shear_g, or the tight-coupling shear). Returns true under tight coupling, when
   *  nothing else is evolved. The layout must have delta_g (radiation streaming off). */
  bool DensityAndVelocityDerivs(const PerturbLayout& layout,
                                const double* y,
                                double* dy,
                                const perturb_parameters_and_workspace& ppaw) const;

  /** The leaf's scalar polarisation columns of the perturbation output, after delta_g, theta_g
   *  and shear_g. Under radiation streaming they are zero; under tight coupling they come from
   *  ppw->tca_shear_g. */
  virtual void AddPolarisationOutput(OutputTable& row,
                                     const PerturbLayout& layout,
                                     const double* y,
                                     double k,
                                     double dkappa,
                                     const perturb_workspace* ppw) const = 0;

  const background& pba_;
  const ThermodynamicsModule* thm_ = nullptr;
  const perturbs* ppt_             = nullptr;

 private:
  int index_tp_delta_ = -1;  // #309 transfer-source slot
  int index_tp_theta_ = -1;
};
