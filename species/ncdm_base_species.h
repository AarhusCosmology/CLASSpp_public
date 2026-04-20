#pragma once
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "../species/base_species.h"
#include "background.h"
#include "parser.h"

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
  // ── Public accessors ──────────────────────────────────────────────────────
  double GetOmega0() const;
  double GetNeff(double z) const;
  double GetMassInElectronvolt() const {
    return m_in_eV_;
  }
  double GetDeg() const {
    return deg_;
  }
  double GetIni(double a, double a_today, double tol_ncdm_initial_w) const;
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

  void PrintNeffInfo() const;
  void PrintMassInfo() const;
  void PrintOmegaInfo() const;

  // Background overrides (shared by NCDMSpecies and DNCDMSpecies)
  void SetBackgroundModule(const BackgroundModule* bgm) override {
    bgm_ = bgm;
  }

 protected:
  // Constructor for standard NCDM species (reads from *_ncdm_standard parameter lists)
  NCDMBaseSpecies(std::string name,
                  EnergyType energy_type,
                  FileContent* pfc,
                  int species_index,
                  const NcdmSettings& settings);

  // Constructor for decay_dr NCDM species (reads from *_ncdm_decay_dr parameter lists)
  // dncdm_index: 0-based index among decay_dr species
  // Pass is_decay_dr=true to select the decay_dr parameter reading path
  NCDMBaseSpecies(std::string name,
                  EnergyType energy_type,
                  FileContent* pfc,
                  int dncdm_index,
                  const NcdmSettings& settings,
                  bool is_decay_dr);

  void SetOmega0(double Omega0, double h);
  void SetDegAndFactor(double deg);
  void SetDeg_from_Omega_ini(double z_ini, double H0, double Omega_ini);

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
  struct DistributionParams {
    const NCDMBaseSpecies* sp = nullptr;
    // For file-based PSD interpolation:
    int tablesize = 0;
    std::vector<double> q, f0, d2f0;
    int last_index = 0;
  };

  void ReadParameters(FileContent* pfc, int species_index, const NcdmSettings& settings);
  void ReadDecayDrParameters(FileContent* pfc, int dncdm_index, const NcdmSettings& settings);
  void InitQuadrature(const NcdmSettings& settings);
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

  int quadrature_strategy_ = 3;
  int input_q_size_        = -1;
  double qmax_             = 15.;
  std::vector<double> psd_parameters_;
  bool got_file_ = false;
  std::string psd_file_;

  double rho_nu_rel_ = 0.;
  double Omega0_     = 0.;
  double omega0_     = 0.;

  mutable ErrorMsg error_message_;
};
