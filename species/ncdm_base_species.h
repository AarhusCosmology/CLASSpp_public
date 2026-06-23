#pragma once
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "../species/base_species.h"
#include "background.h"
#include "parser.h"
#include "quadrature.h"

class BackgroundModule;

#ifndef _zeta3_
#define _zeta3_ \
  1.2020569031595942853997381615114499907649862923404988817922 /**< for quadrature test function */
#endif
#ifndef _zeta5_
#define _zeta5_ \
  1.0369277551433699263313654864570341680570809195019128119741 /**< for quadrature test function */
#endif
#ifndef _PSD_DERIVATIVE_EXP_MIN_
#define _PSD_DERIVATIVE_EXP_MIN_ -30 /**< for ncdm, for accurate computation of dlnf0/dlnq */
#endif
#ifndef _PSD_DERIVATIVE_EXP_MAX_
#define _PSD_DERIVATIVE_EXP_MAX_ 2 /**< for ncdm, for accurate computation of dlnf0/dlnq */
#endif

struct NcdmSettings {
  double h;
  double T_cmb;
  double tol_ncdm;
  double tol_ncdm_bg;
  double tol_M_ncdm;
};

/**
 * Abstract base for all NCDM flavors. Owns per-species quadrature,
 * distribution function, and thermodynamic parameters. Subclasses
 * implement the perturbation interface.
 */
class NCDMBaseSpecies : public BaseSpecies {
 public:
  // ── PerturbLayout ─────────────────────────────────────────────────────────
  /**
   * Per-pv layout for all NCDM-family species.
   * Holds l_max, q_size, and per-q index offsets into pv->y.
   * Owned by perturb_vector::species_layouts (one per pv per species).
   */
  struct PerturbLayout : BaseSpecies::PerturbLayout {
    int l_max  = -1;
    int q_size = -1;
    std::vector<int> index_per_q;  // absolute offsets into pv->y, one per momentum bin
    int total_size() const {
      return q_size * (l_max + 1);
    }
  };

  std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
    return std::make_unique<PerturbLayout>();
  }

  // ── Public accessors ──────────────────────────────────────────────────────
  double GetOmega0() const override;
  double GetNeff(double z) const;
  double GetMassInElectronvolt() const {
    return m_in_eV_;
  }
  double GetDeg() const {
    return deg_;
  }
  /**
   * NCDM (and DNCDM) are "warm" matter — they participate in delta_m but
   * NOT in the delta_cb cold-matter tally. EnergyType::Other → both base
   * defaults would return false, so both explicit overrides are needed.
   */
  bool ClustersAsMatter() const override {
    return true;
  }
  bool IsColdMatterSpecies() const override {
    return false;
  }

  double GetIni(double a, double a_today, double tol_ncdm_initial_w) const;

  /** NCDM may need to start integration earlier so it is relativistic at a_ini. */
  double BackgroundAIni(double a_proposed, double a_today, double tol) const override {
    return GetIni(a_proposed, a_today, tol);
  }

  double GetRescalingFactor(const double* lnf_array) const;
  virtual std::tuple<double, double> GetRescaledParameters(double a, const double* lnf_array) const;

  // Called from NCDMSpecies/DNCDMSpecies background methods
  int ComputeMomenta(
      double z, double* n, double* rho, double* p, double* drho_dM, double* pseudo_p) const;

  int q_size() const {
    return static_cast<int>(q_.size());
  }
  int q_size_bg() const {
    return static_cast<int>(q_bg_.size());
  }

  // Public accessors for quadrature/distribution data (used by module code)
  const std::vector<double>& q() const {
    return q_;
  }
  const std::vector<double>& w() const {
    return w_;
  }
  const std::vector<double>& dlnf0_dlnq() const {
    return dlnf0_dlnq_;
  }
  double M() const {
    return M_;
  }
  double factor() const {
    return factor_;
  }

  void PrintNeffInfo() const override;
  void PrintMassInfo() const override;

  double NeutrinoOmega0() const override {
    return GetOmega0();
  }
  double NeffContribution(double z) const override {
    return GetNeff(z);
  }
  double TensorMasslessRelativisticRho(const double* pvecback) const override {
    return 3. * P(pvecback);
  }
  void CheckUltraRelativisticAtIc(const double* pvecback, double tol) const override;
  bool IsUltraRelativisticAtIc(const double* pvecback, double tol) const override;
  void WarnIfTooHeavyForHalofit(double m_ev_threshold) const override;
  void PrintOmegaInfo() const;

  void SetDeg_from_Omega_ini(double z_ini, double H0, double Omega_ini);

  // Background overrides (shared by NCDMSpecies and DNCDMSpecies)
  void SetBackgroundModule(const BackgroundModule* bgm) override {
    bgm_ = bgm;
  }
  void SetPerturbs(const perturbs* ppt) override {
    ppt_ = ppt;
  }

  void RegisterTransferSourceIndices(int& index_tp, const SourceRequestContext& ctx) override;

  // ── Layout-based tensor register (shared by NCDM and DNCDM) ──────────────
  void RegisterTensorPerturbationIndices(BaseSpecies::PerturbLayout& layout,
                                         perturb_vector* pv,
                                         const precision* ppr,
                                         int& index_pt,
                                         const perturb_workspace* ppw,
                                         int gauge) override;

  // ── Layout-based tensor Boltzmann hierarchy ───────────────────────────────
  /** Tensor-mode Boltzmann hierarchy for NCDM species.
   *  Subclasses must provide GetDlnf0Dlnq to specialise the driving
   *  term (static table for NCDMSpecies; pvecback array for DNCDMSpecies). */
  void PerturbTensorDerivs(const BaseSpecies::PerturbLayout& layout,
                           double tau,
                           const double* y,
                           double* dy,
                           const perturb_parameters_and_workspace& ppaw) const override;

  // ── Tensor GW source contribution ────────────────────────────────────────
  /** Contributes this species' anisotropic-stress integral to ppw->gw_source.
   *  Uses the species's per-pv NCDMBaseSpecies::PerturbLayout slot. */
  void ContributeTensorGwSource(const BaseSpecies::PerturbLayout& layout,
                                double a,
                                double a_today,
                                const double* y,
                                perturb_workspace* ppw) const override;

  // ── MarkUsedInSources ────────────────────────────────────────────────────
  /** NCDM multipoles l > 2 do not enter source functions. */
  void MarkUsedInSources(const BaseSpecies::PerturbLayout& layout,
                         const perturb_workspace* ppw,
                         int* used_in_sources) const override;

  void PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& layout,
                                     double* y,
                                     const PerturbIcContext& ctx) override;

 protected:
  int index_tp_delta_ = -1;  // #309 transfer-source slot (this ncdm flavor)
  int index_tp_theta_ = -1;

  /** Return d ln f_0 / d ln q for momentum bin iq.
   *  NCDMSpecies reads a static table; DNCDMSpecies reads pvecback.
   *  Used by PerturbTensorDerivs. */
  virtual double GetDlnf0Dlnq(int iq, const double* pvecback) const = 0;

  /** Return the quadrature weight w_0[iq] for the GW source integral.
   *  NCDMBaseSpecies returns the static w_[iq]; DNCDMSpecies overrides. */
  virtual double GetW0ForGwSource(int iq, const double* pvecback) const;

  /** Analytic PSD f0(q) for the non-file case. Base: Fermi-Dirac with chemical
   *  potential ksi_. GreyBodyNCDMSpecies overrides this. */
  virtual double EvaluatePsdAnalytic(double q) const;

  /** Default quadrature strategy when the user leaves quadrature_strategy_ at
   *  qm_auto. Base: qm_auto (unchanged). GreyBodyNCDMSpecies returns a
   *  GB-specific method. Wired into InitQuadrature in Task 5.1; not consulted
   *  yet. */
  virtual quadrature_method DefaultQuadratureStrategy() const {
    return qm_auto;
  }

  /** Populate grey-body quadrature parameters. Base: leaves p.active == false. */
  virtual void FillQuadratureParams(GBQuadParams& /*p*/) const {}

  // Constructor: reads parameters per-instance via SpeciesInput (dot-syntax).
  NCDMBaseSpecies(std::string name,
                  EnergyType energy_type,
                  FileContent* pfc,
                  const std::string& instance_name,
                  const NcdmSettings& settings);

  /** Tag type selecting the deferred-init constructor: reads parameters but
   *  does NOT build quadrature or mass, so a subclass can set up an overridden
   *  PSD first and then call BuildQuadratureAndMass itself. */
  struct DeferInit {};

  NCDMBaseSpecies(std::string name,
                  EnergyType energy_type,
                  FileContent* pfc,
                  const std::string& instance_name,
                  const NcdmSettings& settings,
                  DeferInit);

  /** Build perturbation/background quadrature and the dimensionless mass M_.
   *  Calls InitQuadrature (which dispatches the virtual PSD) then computes M_
   *  from m_in_eV_. Call once, after the object is fully constructed. */
  void BuildQuadratureAndMass(const NcdmSettings& settings);

  void InitQuadrature(const NcdmSettings& settings);

  void SetOmega0(double Omega0, double h);
  void SetDegAndFactor(double deg);

  // Momenta variant with variable degeneracy (used by SetDeg_from_Omega_ini):
  int ComputeMomentaDeg(double deg,
                        double z,
                        double* n,
                        double* rho,
                        double* p,
                        double* drho_ddeg,
                        double* pseudo_p) const;

  // Infer M from Omega:
  double MFromOmega(double H0, double Omega0, double tol_M_ncdm) const;

  // Compute dq[i] = w_bg_[i] / f0(q_bg_[i]) for decay_dr species.
  // Must be called after InitQuadrature (i.e., at end of subclass constructor).
  std::vector<double> ComputeDq() const;

  // ── Quadrature — perturbation sampling ───────────────────────────────────
  std::vector<double> q_;
  std::vector<double> w_;
  std::vector<double> dlnf0_dlnq_;

  // ── Quadrature — background sampling ─────────────────────────────────────
  std::vector<double> q_bg_;
  std::vector<double> w_bg_;

  double factor_ = 0.;

  // ── Species thermodynamic parameters ─────────────────────────────────────
  double m_in_eV_ = 0.;
  double M_       = 0.;  // dimensionless mass: m / (k_B * T_ncdm)
  double deg_     = 1.;
  double T_       = 0.71611;  // T_ncdm / T_cmb
  double ksi_     = 0.;       // mu / T_ncdm

  const BackgroundModule* bgm_ = nullptr;
  double T_cmb_                = 0.;
  double h_                    = 0.;

 private:
  const perturbs* ppt_ = nullptr;

  struct DistributionParams {
    const NCDMBaseSpecies* sp = nullptr;
    // For file-based PSD interpolation:
    int tablesize = 0;
    std::vector<double> q, f0, d2f0;
    int last_index = 0;
  };

  void ReadParametersByInstance(FileContent* pfc,
                                const std::string& instance_name,
                                const NcdmSettings& settings);
  void InitDistribution(FileContent* pfc, int species_index);

  int ComputeMomentaMass(double M,
                         double z,
                         double* n,
                         double* rho,
                         double* p,
                         double* drho_dM,
                         double* pseudo_p) const;

  static int DistributionFunction(void* params, double q, double* f0);
  static int TestFunction(void* params, double q, double* test);

  int quadrature_strategy_ = 0;
  int input_q_size_        = 5;
  double qmax_             = 15.;
  std::vector<double> psd_parameters_;
  bool got_file_ = false;
  std::string psd_file_;

  double rho_nu_rel_ = 0.;
  double Omega0_     = 0.;
  double omega0_     = 0.;
};
