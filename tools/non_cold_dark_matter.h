#ifndef NON_COLD_DARK_MATTER_H
#define NON_COLD_DARK_MATTER_H

#include "common.h"
#include "parser.h"
#include "input_module.h"
#include "arrays.h"
#include "quadrature.h"
#include <exception>

#include <memory>

#define _zeta3_ 1.2020569031595942853997381615114499907649862923404988817922 /**< for quandrature test function */
#define _zeta5_ 1.0369277551433699263313654864570341680570809195019128119741 /**< for quandrature test function */
#define _PSD_DERIVATIVE_EXP_MIN_ -30 /**< for ncdm, for accurate computation of dlnf0/dlnq, q step is varied in range specified by these parameters */
#define _PSD_DERIVATIVE_EXP_MAX_ 2  /**< for ncdm, for accurate computation of dlnf0/dlnq, q step is varied in range specified by these parameters */

struct NcdmSettings {
  double h;
  double T_cmb;
  double tol_ncdm;
  double tol_ncdm_bg;
  double tol_M_ncdm;
};

/**
 * Extra properties needed for decaying species
 * Currently only works with NCDM->DR; a different type of struct should be made for the future NCDM->NCDM (or they can be derived from a common struct)
 */

struct DecayDRProperties {
  std::vector<double> dq;
  double Gamma;
  
  int dr_id;     // Index of DR species that this DNCDM species is sourcing
  int dncdm_id;  // Index of ncdm_decay_dr species that this instance describes
  int q_offset;

  int quadrature_strategy;
};

enum class NCDMType { standard, decay_dr }; /** Contains implemented ncdm types */

class NonColdDarkMatter {
public:
  static std::shared_ptr<NonColdDarkMatter> Create(FileContent* pfc, const NcdmSettings&);
  ~NonColdDarkMatter();
  int background_ncdm_momenta(int n_ncdm, double z, double* n, double* rho, double* p, double* drho_dM, double* pseudo_p) const;
  int background_ncdm_momenta_deg(int n_ncdm, double deg, double z, double T_cmb, double* n, double* rho, double* p, double* drho_ddeg, double* pseudo_p) const;
  double GetOmega0() const;
  double GetNeff(double z) const;
  double GetMassInElectronvolt(int n_ncdm) const;
  void PrintNeffInfo() const;
  void PrintMassInfo() const;
  void PrintOmegaInfo() const;
  void SetBackgroundWeight(int n_ncdm, int q_index, double weight);
  void SetOmega0(int n_ncdm, double Omega0, double h);
  void SetDegAndFactor(int n_ncdm, double deg, double T_cmb);
  void SetDeg_from_Omega_ini(int n_ncdm, double z_ini, double H0, double Omega_ini, double T_cmb);
  double GetIni(double a, double a_today, double tol_ncdm_initial_w) const;
  double GetDeg(int n_ncdm) const;
  double GetRescalingFactor(int n_ncdm, double* pvecback_begin) const;
  std::tuple<double, double> GetRescaledParameters(int n_ncdm, double a, double* pvecback_begin);

  std::vector<NCDMType> ncdm_types_; /** Contains information about the types of each ncdm species */
  
  int N_ncdm_ = 0;
  int N_ncdm_standard_ = 0;
  int N_ncdm_decay_dr_ = 0;
  double* m_ncdm_in_eV_ = nullptr;
  double* M_ncdm_ = nullptr;
  int* q_size_ncdm_ = nullptr;    /**< Size of the q_ncdm arrays */
  double* factor_ncdm_ = nullptr; /**< List of normalization factors for calculating energy density etc.*/
  double** q_ncdm_ = nullptr;     /**< Pointers to vectors of perturbation sampling in q */
  double** w_ncdm_ = nullptr;     /**< Pointers to vectors of corresponding quadrature weights w */
  double** dlnf0_dlnq_ncdm_ = nullptr; /**< Pointers to vectors of logarithmic derivatives of p-s-d */
  double* Omega_dncdmdr_ = nullptr;
  
  std::map<int, DecayDRProperties> decay_dr_map_; /**< Map holding decay-wdm-to-dr-specific properties */
  int q_total_size_dncdm_; /**< Total amount of q-bins in DNCDM species */
  
  mutable ErrorMsg error_message_;

private:

  NonColdDarkMatter(FileContent* pfc, const NcdmSettings&);
  int background_ncdm_init(FileContent* pfc, const NcdmSettings&);
  int background_ncdm_momenta_mass(int n_ncdm, double M, double z, double* n, double* rho, double* p, double* drho_dM, double* pseudo_p) const;
  double background_ncdm_M_from_Omega(int n_ncdm, double H0, double Omega0, double tol_M_ncdm);
  static int background_ncdm_distribution(void* pba, double q, double* f0);
  static int background_ncdm_test_function(void* pba, double q, double* test);

  void SafeFree(double* pointer);
  void SafeFree(int* pointer);
  void SafeFree(char* pointer);
  void SafeFree(double** double_pointer);

  const double deg_ncdm_default_ = 1.0;
  const double T_ncdm_default_ = 0.71611; /* this value gives m/omega = 93.14 eV b*/
  const double T_dncdm_default_ = T_ncdm_default_;
  const double ksi_ncdm_default = 0.;

  double rho_nu_rel_;

  double* deg_ncdm_ = nullptr;
  double* T_ncdm_ = nullptr;
  double* ksi_ncdm_ = nullptr;
  double* ncdm_psd_parameters_ = nullptr;

  int* got_files_ = nullptr;
  char* ncdm_psd_files_ = nullptr;
  int* ncdm_quadrature_strategy_ = nullptr; /**< Vector of integers according to quadrature strategy. */
  int* ncdm_input_q_size_ = nullptr; /**< Vector of numbers of q bins */
  double* ncdm_qmax_ = nullptr;   /**< Vector of maximum value of q */
  double* Omega0_ncdm_ = nullptr;
  double* omega0_ncdm_ = nullptr;

  int* q_size_ncdm_bg_ = nullptr; /**< Size of the q_ncdm_bg arrays */
  double** q_ncdm_bg_ = nullptr;  /**< Pointers to vectors of background sampling in q */
  double** w_ncdm_bg_ = nullptr;  /**< Pointers to vectors of corresponding quadrature weights w */
};

/**
 * temporary parameters and workspace passed to phase space distribution function
 */

struct background_parameters_for_distributions {
  background_parameters_for_distributions(const NonColdDarkMatter* ncdm_input, int n_ncdm)
  : ncdm(ncdm_input)
  , n_ncdm(n_ncdm) {}
  ~background_parameters_for_distributions() {
    if (q) free(q);
    if (f0) free(f0);
    if (d2f0) free(d2f0);
  }

  /* structures containing fixed input parameters (indices, ...) */
  const NonColdDarkMatter* const ncdm;
/* Index of current distribution function */
  int n_ncdm;

  /* Used for interpolating in file of tabulated p-s-d: */
  int tablesize = 0;
  double* q = nullptr;
  double* f0 = nullptr;
  double* d2f0 = nullptr;
  int last_index = 0;

};

class NoNcdmRequested : public std::exception {
  std::string message;
public:
  NoNcdmRequested(const char* error) : message(error) {}
  virtual const char* what() const throw () {
         return message.c_str();
  }
};


#endif //NON_COLD_DARK_MATTER_H
