#include "ncdm_base_species.h"

#include <cfloat>  // DBL_EPSILON, DBL_MAX
#include <cmath>
#include <cstring>

#include "background_module.h"
#include "perturbations.h"
#include "perturbations_module.h"
#include "species_input.h"

// ─────────────────────────────────────────────────────────────────────────────
// Constructors + BuildQuadratureAndMass
// ─────────────────────────────────────────────────────────────────────────────

NCDMBaseSpecies::NCDMBaseSpecies(std::string name,
                                 EnergyType energy_type,
                                 FileContent* pfc,
                                 const std::string& instance_name,
                                 const NcdmSettings& settings)
    : BaseSpecies(std::move(name), energy_type), T_cmb_(settings.T_cmb), h_(settings.h) {
  ReadParametersByInstance(pfc, instance_name, settings);
  BuildQuadratureAndMass(settings);
}

NCDMBaseSpecies::NCDMBaseSpecies(std::string name,
                                 EnergyType energy_type,
                                 FileContent* pfc,
                                 const std::string& instance_name,
                                 const NcdmSettings& settings,
                                 DeferInit)
    : BaseSpecies(std::move(name), energy_type), T_cmb_(settings.T_cmb), h_(settings.h) {
  ReadParametersByInstance(pfc, instance_name, settings);
  // Caller (subclass) is responsible for BuildQuadratureAndMass after setting up its PSD.
}

void NCDMBaseSpecies::BuildQuadratureAndMass(const NcdmSettings& settings) {
  InitQuadrature(settings);

  if (m_in_eV_ != 0.0) {
    M_ = m_in_eV_ / _k_B_ * _eV_ / T_ / T_cmb_;
    double rho_ncdm;
    ComputeMomenta(0., nullptr, &rho_ncdm, nullptr, nullptr, nullptr);
  }
  else {
    M_ = 0.;
  }
}

// ─────────────────────────────────────────────────────────────────────────────
// ReadParametersByInstance — read parameters for one species via SpeciesInput
// ─────────────────────────────────────────────────────────────────────────────

void NCDMBaseSpecies::ReadParametersByInstance(FileContent* pfc,
                                               const std::string& instance_name,
                                               const NcdmSettings& settings) {
  SpeciesInput input(pfc, instance_name);

  // Common scalar fields. read_double / read_int leave the destination unchanged
  // on miss, so each member keeps its in-class default. (Defaults that were
  // previously fixed-up here, like T_ncdm_default = 0.71611, must already be set
  // as in-class member initialisers — see Task 4.)
  input.read_double("m", m_in_eV_);
  input.read_double("T", T_);
  input.read_double("deg", deg_);
  input.read_double("Omega", Omega0_);
  input.read_double("omega", omega0_);
  input.read_double("ksi", ksi_);
  input.read_int("quadrature_strategy", quadrature_strategy_);
  input.read_int("momenta_bins", input_q_size_);
  input.read_double("max_q", qmax_);

  // psd_parameters: variable-length list (single instance, comma-separated values).
  std::vector<double> psd_params;
  if (input.read_list_of_doubles("psd_parameters", psd_params)) {
    psd_parameters_ = std::move(psd_params);
  }

  // PSD file: a single-instance flag and an optional filename.
  int use_psd_file = 0;
  input.read_int("use_psd_file", use_psd_file);
  got_file_ = (use_psd_file != 0);
  if (got_file_) {
    std::string filename;
    if (!input.read_string("psd_filename", filename)) {
      throw std::invalid_argument("species '" + instance_name +
                                  "': use_psd_file=1 but psd_filename is missing");
    }
    psd_file_ = std::move(filename);
  }

  // Resolve omega/Omega conflict (matches existing semantics).
  if (omega0_ != 0.0) {
    if (Omega0_ != 0.0) {
      throw std::invalid_argument("species '" + instance_name +
                                  "': both Omega and omega specified — choose one");
    }
    Omega0_ = omega0_ / settings.h / settings.h;
  }
  else {
    omega0_ = Omega0_ * settings.h * settings.h;
  }

  // Ultra-relativistic default: if both Omega and m are zero, give a tiny mass.
  if ((Omega0_ == 0.0) && (m_in_eV_ == 0.0)) {
    m_in_eV_ = 1.e-5;
  }
}

// ─────────────────────────────────────────────────────────────────────────────
// InitQuadrature
// ─────────────────────────────────────────────────────────────────────────────

void NCDMBaseSpecies::InitQuadrature(const NcdmSettings& settings) {
  DistributionParams pbadist;
  pbadist.sp        = this;
  pbadist.tablesize = 0;

  // If file-based PSD, read file and build spline table
  if (got_file_) {
    FILE* psdfile = fopen(psd_file_.c_str(), "r");
    class_test(psdfile == nullptr, error_message_, "Could not open file %s!", psd_file_.c_str());
    double tmp1, tmp2;
    int row = 0;
    while (fscanf(psdfile, "%lf %lf", &tmp1, &tmp2) == 2)
      ++row;
    rewind(psdfile);
    pbadist.tablesize = row;
    class_test(pbadist.tablesize < 2,
               error_message_,
               "PSD file '%s' has %d data row(s); need at least 2 for spline interpolation",
               psd_file_.c_str(),
               row);

    pbadist.q.resize(pbadist.tablesize);
    pbadist.f0.resize(pbadist.tablesize);
    pbadist.d2f0.resize(pbadist.tablesize);
    for (row = 0; row < pbadist.tablesize; row++) {
      fscanf(psdfile, "%lf %lf", &pbadist.q[row], &pbadist.f0[row]);
    }
    fclose(psdfile);
    class_call(array_spline_table_lines(pbadist.q.data(),
                                        pbadist.tablesize,
                                        pbadist.f0.data(),
                                        1,
                                        pbadist.d2f0.data(),
                                        _SPLINE_EST_DERIV_,
                                        error_message_),
               error_message_,
               error_message_);
  }

  // Grey-body subclasses fill these and select a GB-specific default strategy;
  // for all other species gb.active stays false and strategy == quadrature_strategy_.
  GBQuadParams gb;
  FillQuadratureParams(gb);
  int strategy = quadrature_strategy_;
  if (strategy == qm_auto && gb.active) {
    strategy = DefaultQuadratureStrategy();
  }

  // Handle perturbation quadrature
  if (strategy == qm_auto) {
    q_.resize(_QUADRATURE_MAX_);
    w_.resize(_QUADRATURE_MAX_);
    int q_size;
    class_call(get_qsampling(q_.data(),
                             w_.data(),
                             &q_size,
                             _QUADRATURE_MAX_,
                             settings.tol_ncdm,
                             pbadist.q.data(),
                             pbadist.tablesize,
                             TestFunction,
                             DistributionFunction,
                             &pbadist,
                             error_message_),
               error_message_,
               error_message_);
    q_.resize(q_size);
    w_.resize(q_size);

    q_bg_.resize(_QUADRATURE_MAX_BG_);
    w_bg_.resize(_QUADRATURE_MAX_BG_);
    int q_size_bg;
    class_call(get_qsampling(q_bg_.data(),
                             w_bg_.data(),
                             &q_size_bg,
                             _QUADRATURE_MAX_BG_,
                             settings.tol_ncdm_bg,
                             pbadist.q.data(),
                             pbadist.tablesize,
                             TestFunction,
                             DistributionFunction,
                             &pbadist,
                             error_message_),
               error_message_,
               error_message_);
    q_bg_.resize(q_size_bg);
    w_bg_.resize(q_size_bg);
  }
  else {
    // Manual sampling
    q_.resize(input_q_size_);
    w_.resize(input_q_size_);
    std::vector<double> dq_dummy(input_q_size_);
    class_call(get_qsampling_manual(q_.data(),
                                    w_.data(),
                                    dq_dummy.data(),
                                    input_q_size_,
                                    qmax_,
                                    (enum quadrature_method) strategy,
                                    pbadist.q.data(),
                                    pbadist.tablesize,
                                    DistributionFunction,
                                    &pbadist,
                                    gb,
                                    error_message_),
               error_message_,
               error_message_);
    q_bg_ = q_;
    w_bg_ = w_;
  }

  // Compute dlnf0_dlnq_
  dlnf0_dlnq_.resize(q_.size());
  for (int index_q = 0; index_q < static_cast<int>(q_.size()); index_q++) {
    double q = q_[index_q];
    double f0;
    class_call(DistributionFunction(&pbadist, q, &f0), error_message_, error_message_);

    double dq = 1., f0m2 = 0., f0p2 = 0.;
    for (int tolexp = _PSD_DERIVATIVE_EXP_MIN_; tolexp < _PSD_DERIVATIVE_EXP_MAX_; tolexp++) {
      if (index_q == 0)
        dq = MIN((0.5 - DBL_EPSILON) * q, 2 * exp(tolexp) * (q_[index_q + 1] - q));
      else if (index_q == static_cast<int>(q_.size()) - 1)
        dq = exp(tolexp) * 2.0 * (q_[index_q] - q_[index_q - 1]);
      else
        dq = exp(tolexp) * (q_[index_q + 1] - q_[index_q - 1]);

      class_call(DistributionFunction(&pbadist, q - 2 * dq, &f0m2), error_message_, error_message_);
      class_call(DistributionFunction(&pbadist, q + 2 * dq, &f0p2), error_message_, error_message_);
      if (fabs((f0p2 - f0m2) / f0) > sqrt(DBL_EPSILON))
        break;
    }

    double f0m1, f0p1;
    class_call(DistributionFunction(&pbadist, q - dq, &f0m1), error_message_, error_message_);
    class_call(DistributionFunction(&pbadist, q + dq, &f0p1), error_message_, error_message_);
    double df0dq = (+f0m2 - 8 * f0m1 + 8 * f0p1 - f0p2) / 12.0 / dq;
    if (fabs(f0) == 0.)
      dlnf0_dlnq_[index_q] = -q;
    else
      dlnf0_dlnq_[index_q] = q / f0 * df0dq;
  }

  // Compute factor_
  factor_ = deg_ * 4 * _PI_ * pow(T_cmb_ * T_ * _k_B_, 4) * 8 * _PI_ * _G_ / 3. /
            pow(_h_P_ / 2. / _PI_, 3) / pow(_c_, 7) * _Mpc_over_m_ * _Mpc_over_m_;

  // Compute rho_nu_rel_ (relativistic reference density)
  rho_nu_rel_ = 56.0 / 45.0 * pow(_PI_, 6) * pow(4.0 / 11.0, 4.0 / 3.0) * _G_ / pow(_h_P_, 3) /
                pow(_c_, 7) * pow(_Mpc_over_m_, 2) * pow(T_cmb_ * _k_B_, 4);
}

// ─────────────────────────────────────────────────────────────────────────────
// InitDistribution — file reading done inline in InitQuadrature
// ─────────────────────────────────────────────────────────────────────────────

void NCDMBaseSpecies::InitDistribution(FileContent* /*pfc*/, int /*species_index*/) {}

// ─────────────────────────────────────────────────────────────────────────────
// EvaluatePsdAnalytic — virtual, overridable Fermi-Dirac PSD
// ─────────────────────────────────────────────────────────────────────────────

double NCDMBaseSpecies::EvaluatePsdAnalytic(double q) const {
  // Fermi-Dirac with chemical potential.
  return 1.0 / pow(2 * _PI_, 3) * (1. / (exp(q - ksi_) + 1.) + 1. / (exp(q + ksi_) + 1.));
}

// ─────────────────────────────────────────────────────────────────────────────
// DistributionFunction — static callback
// ─────────────────────────────────────────────────────────────────────────────

int NCDMBaseSpecies::DistributionFunction(void* params, double q, double* f0) {
  auto* p                   = static_cast<DistributionParams*>(params);
  const NCDMBaseSpecies* sp = p->sp;

  if (p->tablesize > 0) {
    int lastidx = p->tablesize - 1;
    if (q < p->q[0]) {
      *f0 = p->f0[0];
    }
    else if (q > p->q[lastidx]) {
      double qlast   = p->q[lastidx];
      double f0last  = p->f0[lastidx];
      double dqlast  = qlast - p->q[lastidx - 1];
      double df0last = f0last - p->f0[lastidx - 1];
      *f0            = f0last * exp(-(qlast - q) * df0last / f0last / dqlast);
    }
    else {
      class_call(array_interpolate_spline(p->q.data(),
                                          p->tablesize,
                                          p->f0.data(),
                                          p->d2f0.data(),
                                          1,
                                          q,
                                          &p->last_index,
                                          f0,
                                          1,
                                          const_cast<char*>(sp->error_message_)),
                 const_cast<char*>(sp->error_message_),
                 const_cast<char*>(sp->error_message_));
    }
  }
  else {
    *f0 = sp->EvaluatePsdAnalytic(q);
  }
  return _SUCCESS_;
}

// ─────────────────────────────────────────────────────────────────────────────
// TestFunction — static callback
// ─────────────────────────────────────────────────────────────────────────────

int NCDMBaseSpecies::TestFunction(void* /*params*/, double q, double* test) {
  double c = 2.0 / (3.0 * _zeta3_);
  double d = 120.0 / (7.0 * pow(_PI_, 4));
  double e = 2.0 / (45.0 * _zeta5_);
  *test    = pow(2.0 * _PI_, 3) / 6.0 * (c * q * q - d * q * q * q - e * q * q * q * q);
  return _SUCCESS_;
}

// ─────────────────────────────────────────────────────────────────────────────
// ComputeMomentaMass — background momenta integration at given M
// ─────────────────────────────────────────────────────────────────────────────

int NCDMBaseSpecies::ComputeMomentaMass(double M,
                                        double z,
                                        double* n,
                                        double* rho,
                                        double* p,
                                        double* drho_dM,
                                        double* pseudo_p) const {
  double factor2     = factor_ * pow(1 + z, 4);
  const double* qvec = q_bg_.data();
  const double* wvec = w_bg_.data();
  int qsize          = static_cast<int>(q_bg_.size());

  if (n != nullptr)
    *n = 0.;
  if (rho != nullptr)
    *rho = 0.;
  if (p != nullptr)
    *p = 0.;
  if (drho_dM != nullptr)
    *drho_dM = 0.;
  if (pseudo_p != nullptr)
    *pseudo_p = 0.;

  for (int index_q = 0; index_q < qsize; index_q++) {
    double q2      = qvec[index_q] * qvec[index_q];
    double epsilon = sqrt(q2 + M * M / (1. + z) / (1. + z));
    if (n != nullptr)
      *n += q2 * wvec[index_q];
    if (rho != nullptr)
      *rho += q2 * epsilon * wvec[index_q];
    if (p != nullptr)
      *p += q2 * q2 / 3. / epsilon * wvec[index_q];
    if (drho_dM != nullptr)
      *drho_dM += q2 * M / (1. + z) / (1. + z) / epsilon * wvec[index_q];
    if (pseudo_p != nullptr)
      *pseudo_p += pow(q2 / epsilon, 3) / 3. * wvec[index_q];
  }

  if (n != nullptr)
    *n *= factor2 / (1. + z);
  if (rho != nullptr)
    *rho *= factor2;
  if (p != nullptr)
    *p *= factor2;
  if (drho_dM != nullptr)
    *drho_dM *= factor2;
  if (pseudo_p != nullptr)
    *pseudo_p *= factor2;

  return _SUCCESS_;
}

// ─────────────────────────────────────────────────────────────────────────────
// ComputeMomenta — uses stored M_
// ─────────────────────────────────────────────────────────────────────────────

int NCDMBaseSpecies::ComputeMomenta(
    double z, double* n, double* rho, double* p, double* drho_dM, double* pseudo_p) const {
  return ComputeMomentaMass(M_, z, n, rho, p, drho_dM, pseudo_p);
}

// ─────────────────────────────────────────────────────────────────────────────
// ComputeMomentaDeg — same as ComputeMomentaMass but with variable degeneracy
// ─────────────────────────────────────────────────────────────────────────────

int NCDMBaseSpecies::ComputeMomentaDeg(double deg,
                                       double z,
                                       double* n,
                                       double* rho,
                                       double* p,
                                       double* drho_ddeg,
                                       double* pseudo_p) const {
  double factor2     = deg * 4 * _PI_ * pow(T_cmb_ * T_ * _k_B_, 4) * 8 * _PI_ * _G_ / 3. /
                       pow(_h_P_ / 2. / _PI_, 3) / pow(_c_, 7) * _Mpc_over_m_ * _Mpc_over_m_ *
                       pow(1 + z, 4);
  const double* qvec = q_bg_.data();
  const double* wvec = w_bg_.data();
  int qsize          = static_cast<int>(q_bg_.size());

  if (n != nullptr)
    *n = 0.;
  if (rho != nullptr)
    *rho = 0.;
  if (p != nullptr)
    *p = 0.;
  if (pseudo_p != nullptr)
    *pseudo_p = 0.;

  for (int index_q = 0; index_q < qsize; index_q++) {
    double q2      = qvec[index_q] * qvec[index_q];
    double epsilon = sqrt(q2 + M_ * M_ / (1. + z) / (1. + z));
    if (n != nullptr)
      *n += q2 * wvec[index_q];
    if (rho != nullptr)
      *rho += q2 * epsilon * wvec[index_q];
    if (p != nullptr)
      *p += q2 * q2 / 3. / epsilon * wvec[index_q];
    if (pseudo_p != nullptr)
      *pseudo_p += pow(q2 / epsilon, 3) / 3.0 * wvec[index_q];
  }

  if (n != nullptr)
    *n *= factor2 / (1. + z);
  if (rho != nullptr)
    *rho *= factor2;
  if (p != nullptr)
    *p *= factor2;
  if (pseudo_p != nullptr)
    *pseudo_p *= factor2;

  // rho is linear in deg
  if (drho_ddeg != nullptr) {
    *drho_ddeg = (rho != nullptr) ? (*rho / deg) : 0.;
    if (rho == nullptr) {
      // Compute rho separately if needed
      double rho_tmp = 0.;
      for (int index_q = 0; index_q < qsize; index_q++) {
        double q2       = qvec[index_q] * qvec[index_q];
        double epsilon  = sqrt(q2 + M_ * M_ / (1. + z) / (1. + z));
        rho_tmp        += q2 * epsilon * wvec[index_q];
      }
      *drho_ddeg = rho_tmp * factor2 / deg;
    }
  }

  return _SUCCESS_;
}

// ─────────────────────────────────────────────────────────────────────────────
// MFromOmega — Newton iteration to find dimensionless mass from Omega
// ─────────────────────────────────────────────────────────────────────────────

double NCDMBaseSpecies::MFromOmega(double H0, double Omega0, double tol_M_ncdm) const {
  int maxiter = 50;
  double rho0 = H0 * H0 * Omega0;
  double M    = 0;
  double rho;
  double n;
  ComputeMomentaMass(M, 0., &n, &rho, nullptr, nullptr, nullptr);

  class_test(rho0 < rho,
             const_cast<char*>(error_message_),
             "The value of Omega for this species, %g, is less than for a massless species! "
             "It should be at least %g. Check your input.",
             Omega0,
             Omega0 * rho / rho0);

  M = rho0 / n;  // zeroth order guess (strict NR limit)
  for (int iter = 1; iter <= maxiter; iter++) {
    double drhodM;
    ComputeMomentaMass(M, 0., nullptr, &rho, nullptr, &drhodM, nullptr);
    double deltaM = (rho0 - rho) / drhodM;
    if ((M + deltaM) < 0.0) {
      deltaM = -M / 2.0;
    }
    M += deltaM;
    if (fabs(deltaM / M) < tol_M_ncdm) {
      return M;
    }
  }
  ThrowRuntimeError("Newton iteration could not converge on a mass for some reason.");
  return 0.;
}

// ─────────────────────────────────────────────────────────────────────────────
// GetOmega0, GetNeff, GetIni, GetRescalingFactor, GetRescaledParameters
// ─────────────────────────────────────────────────────────────────────────────

double NCDMBaseSpecies::GetOmega0() const {
  return Omega0_;
}

double NCDMBaseSpecies::GetNeff(double z) const {
  double rho_ncdm_rel;
  ComputeMomentaMass(0., z, nullptr, &rho_ncdm_rel, nullptr, nullptr, nullptr);
  return rho_ncdm_rel / rho_nu_rel_;
}

double NCDMBaseSpecies::GetIni(double a, double a_today, double tol_ncdm_initial_w) const {
  double rho_ncdm, p_ncdm;
  bool converged = false;
  for (int counter = 0; counter < _MAX_IT_; counter++) {
    ComputeMomenta(a_today / a - 1.0, nullptr, &rho_ncdm, &p_ncdm, nullptr, nullptr);
    if (fabs(p_ncdm / rho_ncdm - 1. / 3.) <= tol_ncdm_initial_w) {
      converged = true;
      break;
    }
    a *= _SCALE_BACK_;
  }
  class_test(converged == false,
             error_message_,
             "Search for initial scale factor a such that ncdm species is relativistic failed.");
  return a;
}

double NCDMBaseSpecies::GetRescalingFactor(const double* lnf_array) const {
  double lnN = DBL_MAX;
  for (int index_q = 0; index_q < static_cast<int>(q_.size()); index_q++) {
    lnN = std::min(lnN, -lnf_array[index_q]);
  }
  return lnN + 400.;
}

std::tuple<double, double> NCDMBaseSpecies::GetRescaledParameters(double a,
                                                                  const double* lnf_array) const {
  // dq_ (decay_dr quadrature step) is only in DNCDMSpecies; base is a no-op.
  (void) a;
  (void) lnf_array;
  return {0., 0.};
}

// ─────────────────────────────────────────────────────────────────────────────
// ComputeDq — dq_[i] = w_bg_[i] / f0(q_bg_[i]) for decay_dr species
// ─────────────────────────────────────────────────────────────────────────────

std::vector<double> NCDMBaseSpecies::ComputeDq() const {
  DistributionParams pbadist;
  pbadist.sp        = this;
  pbadist.tablesize = 0;
  // (file-based PSD would need re-reading; DNCDM always uses Fermi-Dirac so tablesize==0)

  const int qsize = static_cast<int>(q_bg_.size());
  std::vector<double> dq(qsize);
  for (int i = 0; i < qsize; i++) {
    double f0 = 0.;
    DistributionFunction(&pbadist, q_bg_[i], &f0);
    dq[i] = (f0 != 0.) ? w_bg_[i] / f0 : 0.;
  }
  return dq;
}

// ─────────────────────────────────────────────────────────────────────────────
// SetOmega0, SetDegAndFactor, SetDeg_from_Omega_ini
// ─────────────────────────────────────────────────────────────────────────────

void NCDMBaseSpecies::SetOmega0(double Omega0, double h) {
  Omega0_ = Omega0;
  omega0_ = Omega0 * h * h;
}

void NCDMBaseSpecies::SetDegAndFactor(double deg) {
  deg_    = deg;
  factor_ = deg_ * 4 * _PI_ * pow(T_cmb_ * T_ * _k_B_, 4) * 8 * _PI_ * _G_ / 3. /
            pow(_h_P_ / 2. / _PI_, 3) / pow(_c_, 7) * _Mpc_over_m_ * _Mpc_over_m_;
}

void NCDMBaseSpecies::SetDeg_from_Omega_ini(double z_ini, double H0, double Omega_ini) {
  double rho_deg1;
  ComputeMomentaDeg(1.0, z_ini, nullptr, &rho_deg1, nullptr, nullptr, nullptr);
  double Omega_deg1 = rho_deg1 * pow(1 + z_ini, -4.0) / H0 / H0;
  deg_              = Omega_ini / Omega_deg1;
  factor_           = deg_ * 4 * _PI_ * pow(T_cmb_ * T_ * _k_B_, 4) * 8 * _PI_ * _G_ / 3. /
                      pow(_h_P_ / 2. / _PI_, 3) / pow(_c_, 7) * _Mpc_over_m_ * _Mpc_over_m_;
}

// ─────────────────────────────────────────────────────────────────────────────
// Print methods
// ─────────────────────────────────────────────────────────────────────────────

void NCDMBaseSpecies::PrintNeffInfo() const {
  if (got_file_) {
    printf(" -> ncdm species read from file %s\n", psd_file_.c_str());
  }
  double rho_ncdm_rel;
  ComputeMomentaMass(0., 0., nullptr, &rho_ncdm_rel, nullptr, nullptr, nullptr);
  printf(
      " -> ncdm species sampled with %d (resp. %d) points for purpose of background "
      "(resp. perturbation) integration. In the relativistic limit it gives Delta N_eff = %g\n",
      static_cast<int>(q_bg_.size()),
      static_cast<int>(q_.size()),
      rho_ncdm_rel / rho_nu_rel_);
}

void NCDMBaseSpecies::PrintMassInfo() const {
  printf(" -> non-cold dark matter species has m_i = %e eV (so m_i / omega_i = %e eV)\n",
         m_in_eV_,
         m_in_eV_ * deg_ / omega0_);
}

void NCDMBaseSpecies::PrintOmegaInfo() const {
  printf("-> %-26s Omega = %-15g , omega = %-15g\n", "Neutrino Species", Omega0_, omega0_);
}

// ─────────────────────────────────────────────────────────────────────────────
// RegisterTensorPerturbationIndices (layout-based)
// Reads ppr->l_max_ncdm and q_size() directly; no caller pre-setup required.
// Thread-safe: layout/pv are per-thread.
// ─────────────────────────────────────────────────────────────────────────────

void NCDMBaseSpecies::RegisterTensorPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                        perturb_vector* /*pv*/,
                                                        const precision* ppr,
                                                        int& index_pt,
                                                        const perturb_workspace* /*ppw*/,
                                                        int /*gauge*/) {
  auto& layout = static_cast<PerturbLayout&>(base);

  /* NCDM tensor slots are only reserved when tensor_method == tm_exact.
     Other tensor methods (photons-only, massless approximation) do not
     need NCDM tensor multipoles. */
  if (ppt_ == nullptr || ppt_->tensor_method != tm_exact)
    return;

  layout.l_max  = ppr->l_max_ncdm;
  layout.q_size = q_size();

  layout.index_per_q.clear();
  layout.index_per_q.reserve(layout.q_size);
  for (int iq = 0; iq < layout.q_size; ++iq)
    layout.index_per_q.push_back(index_pt + iq * (layout.l_max + 1));

  index_pt += layout.total_size();
}

// ─────────────────────────────────────────────────────────────────────────────
// MarkUsedInSources
// NCDM multipoles l > 2 do not enter source functions.
// ─────────────────────────────────────────────────────────────────────────────

void NCDMBaseSpecies::MarkUsedInSources(const BaseSpecies::PerturbLayout& base,
                                        const perturb_workspace* /*ppw*/,
                                        int* used_in_sources) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  /* NCDMBaseSpecies scalar mask: multipoles l > 2 do not enter source functions. */
  if (layout.q_size < 0)
    return;
  for (int iq = 0; iq < layout.q_size; ++iq) {
    const int b = layout.index_per_q[iq];
    for (int l = 0; l <= layout.l_max; ++l) {
      if (l > 2)
        used_in_sources[b + l] = _FALSE_;
    }
  }
}

// ─────────────────────────────────────────────────────────────────────────────
// PerturbSynchronousToNewtonian
// Per-q synchronous→Newtonian gauge shift for all NCDM-family species.
// delta_ur/theta_ur are the shared free-streaming radiation IC, gauge-shifted
// as radiation (delta: -4ℋα; theta: +k²α). shear_ur/l3_ur are gauge-invariant.
// ─────────────────────────────────────────────────────────────────────────────

void NCDMBaseSpecies::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                                    double* y,
                                                    const PerturbIcContext& ctx) {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  if (layout.q_size <= 0 || layout.index_per_q.empty())
    return;
  const double* pvecback = ctx.ppw->pvecback;
  const double delta_ur  = ctx.delta_ur - 4. * ctx.a_prime_over_a * ctx.alpha;
  const double theta_ur  = ctx.theta_ur + ctx.k * ctx.k * ctx.alpha;
  const int lmax         = layout.l_max;
  for (int index_q = 0; index_q < layout.q_size; ++index_q) {
    const int idx           = layout.index_per_q[index_q];
    const double q          = q_[index_q];
    const double epsilon    = std::sqrt(q * q + ctx.a * ctx.a * M_ * M_);
    const double dlnf0_dlnq = GetDlnf0Dlnq(index_q, pvecback);
    y[idx + 0]              = -0.25 * delta_ur * dlnf0_dlnq;
    if (lmax >= 1)
      y[idx + 1] = -epsilon / 3. / q / ctx.k * theta_ur * dlnf0_dlnq;
    if (lmax >= 2)
      y[idx + 2] = -0.5 * ctx.shear_ur * dlnf0_dlnq;
    if (lmax >= 3)
      y[idx + 3] = -0.25 * ctx.l3_ur * dlnf0_dlnq;
  }
}

// ─────────────────────────────────────────────────────────────────────────────
// PerturbTensorDerivs (layout-based)
// Ports the per-species tensor Boltzmann hierarchy from perturb_derivs_member.
// GetDlnf0Dlnq is supplied by the concrete subclass:
//   NCDMSpecies → static dlnf0_dlnq_[iq]
//   DNCDMSpecies → pvecback[bg_dlnfdlnq_index + iq]
// ─────────────────────────────────────────────────────────────────────────────

void NCDMBaseSpecies::PerturbTensorDerivs(const BaseSpecies::PerturbLayout& base,
                                          double /*tau*/,
                                          const double* y,
                                          double* dy,
                                          const perturb_parameters_and_workspace& ppaw) {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  if (layout.q_size <= 0 || layout.index_per_q.empty())
    return;

  const perturb_workspace* ppw = ppaw.ppw;
  const perturb_vector* pv     = ppw->pv;
  const double* pvecback       = ppw->pvecback;
  const double* s_l            = ppw->s_l;
  const double k               = ppaw.k;
  const double a               = pvecback[bgm_->index_bg_a_];
  const double a2              = a * a;
  const double cotKgen         = ppw->cotKgen;
  const int lmax               = layout.l_max;

  for (int iq = 0; iq < layout.q_size; ++iq) {
    const double q              = q_[iq];
    const double dlnf0_dlnq     = GetDlnf0Dlnq(iq, pvecback);
    const double epsilon        = std::sqrt(q * q + a2 * M_ * M_);
    const double qk_div_epsilon = k * q / epsilon;
    const int ncdm_idx          = layout.index_per_q[iq];

    dy[ncdm_idx] = -qk_div_epsilon * y[ncdm_idx + 1] -
                   0.25 * _SQRT6_ * y[pv->index_pt_gwdot] * dlnf0_dlnq;

    for (int l = 1; l < lmax; ++l) {
      dy[ncdm_idx + l] = qk_div_epsilon / (2. * l + 1.0) *
                         (l * s_l[l] * y[ncdm_idx + (l - 1)] -
                          (l + 1.) * s_l[l + 1] * y[ncdm_idx + (l + 1)]);
    }

    dy[ncdm_idx + lmax] = qk_div_epsilon * y[ncdm_idx + lmax - 1] -
                          (1. + lmax) * k * cotKgen * y[ncdm_idx + lmax];
  }
}

// ─────────────────────────────────────────────────────────────────────────────
// GetW0ForGwSource  (base-class version: static quadrature weight)
// DNCDMSpecies overrides to use pvecback[bg_lnf_index + iq].
// ─────────────────────────────────────────────────────────────────────────────

double NCDMBaseSpecies::GetW0ForGwSource(int iq, const double* /*pvecback*/) const {
  return w_[iq];
}

// ─────────────────────────────────────────────────────────────────────────────
// ContributeTensorGwSource
// Ports the per-species GW source integral from perturb_total_stress_energy.
// Uses the species's NCDMBaseSpecies::PerturbLayout slot from ppw->pv.
// ─────────────────────────────────────────────────────────────────────────────

void NCDMBaseSpecies::ContributeTensorGwSource(const BaseSpecies::PerturbLayout& base,
                                               double a,
                                               double a_today,
                                               const double* y,
                                               perturb_workspace* ppw) const {
  const auto& layout     = static_cast<const PerturbLayout&>(base);
  const double* pvecback = ppw->pvecback;
  const double a2        = a * a;
  const double factor    = factor_ * std::pow(a_today / a, 4);

  double gwncdm = 0.;
  for (int iq = 0; iq < layout.q_size; ++iq) {
    const int idx        = layout.index_per_q[iq];
    const double w0      = GetW0ForGwSource(iq, pvecback);
    const double q       = q_[iq];
    const double q2      = q * q;
    const double epsilon = std::sqrt(q2 + M_ * M_ * a2);

    gwncdm += q2 * q2 / epsilon * w0 *
              (1. / 15. * y[idx] + 2. / 21. * y[idx + 2] + 1. / 35. * y[idx + 4]);
  }

  gwncdm         *= -_SQRT6_ * 4 * a2 * factor;
  ppw->gw_source += gwncdm;
}
