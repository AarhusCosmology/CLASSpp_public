#include "ncdm_base_species.h"

#include <cfloat>  // DBL_EPSILON, DBL_MAX
#include <cstring>

namespace {

int readDoubleList(
    FileContent* pfc, const char* name, std::vector<double>& values, int* found, ErrorMsg errmsg) {
  try {
    *found = pfc->read_list_of_doubles(name, values) ? _TRUE_ : _FALSE_;
  }
  catch (const std::exception& e) {
    class_stop(errmsg, "%s", e.what());
  }
  return _SUCCESS_;
}

int readIntegerList(
    FileContent* pfc, const char* name, std::vector<int>& values, int* found, ErrorMsg errmsg) {
  try {
    *found = pfc->read_list_of_integers(name, values) ? _TRUE_ : _FALSE_;
  }
  catch (const std::exception& e) {
    class_stop(errmsg, "%s", e.what());
  }
  return _SUCCESS_;
}

int readStringList(FileContent* pfc,
                   const char* name,
                   std::vector<std::string>& values,
                   int* found,
                   ErrorMsg errmsg) {
  try {
    *found = pfc->read_list_of_strings(name, values) ? _TRUE_ : _FALSE_;
  }
  catch (const std::exception& e) {
    class_stop(errmsg, "%s", e.what());
  }
  return _SUCCESS_;
}

}  // namespace

// ─────────────────────────────────────────────────────────────────────────────
// Constructor
// ─────────────────────────────────────────────────────────────────────────────

NCDMBaseSpecies::NCDMBaseSpecies(std::string name,
                                 EnergyType energy_type,
                                 FileContent* pfc,
                                 int species_index,
                                 const NcdmSettings& settings,
                                 const std::string& suffix)
    : BaseSpecies(std::move(name), energy_type), T_cmb_(settings.T_cmb), h_(settings.h) {
  ReadParameters(pfc, species_index, suffix, settings);
}

NCDMBaseSpecies::NCDMBaseSpecies(std::string name,
                                 EnergyType energy_type,
                                 FileContent* pfc,
                                 int dncdm_index,
                                 const NcdmSettings& settings,
                                 bool /*is_decay_dr*/)
    : BaseSpecies(std::move(name), energy_type), T_cmb_(settings.T_cmb), h_(settings.h) {
  ReadDecayDrParameters(pfc, dncdm_index, settings);
}

// ─────────────────────────────────────────────────────────────────────────────
// ReadDecayDrParameters — read parameters for one decay_dr species
// ─────────────────────────────────────────────────────────────────────────────

void NCDMBaseSpecies::ReadDecayDrParameters(FileContent* pfc,
                                            int dncdm_index,
                                            const NcdmSettings& settings) {
  char* errmsg = error_message_;

  // Read N_ncdm_decay_dr so we know how many entries to expect
  int N_ncdm_decay_dr = 0;
  int flag;
  class_call(parser_read_int(pfc, "N_ncdm_decay_dr", &N_ncdm_decay_dr, &flag, errmsg),
             errmsg,
             errmsg);

  // Helper: read a list of ints, fill with default if not found, return element [dncdm_index]
  auto read_int_elem = [&](const char* key, int default_value) -> int {
    int found;
    std::vector<int> vals;
    readIntegerList(pfc, key, vals, &found, errmsg);
    if (found == _FALSE_) {
      return default_value;
    }
    if (dncdm_index < static_cast<int>(vals.size())) {
      return vals[dncdm_index];
    }
    return default_value;
  };

  auto read_double_elem = [&](const char* key, double default_value) -> double {
    int found;
    std::vector<double> vals;
    readDoubleList(pfc, key, vals, &found, errmsg);
    if (found == _FALSE_) {
      return default_value;
    }
    if (dncdm_index < static_cast<int>(vals.size())) {
      return vals[dncdm_index];
    }
    return default_value;
  };

  const double T_dncdm_default = 0.71611;

  quadrature_strategy_ = read_int_elem("quadrature_strategy_ncdm_decay_dr", 0);
  input_q_size_        = read_int_elem("N_momentum_bins_ncdm_decay_dr", 5);
  qmax_                = read_double_elem("maximum_q_ncdm_decay_dr", 15.0);
  T_                   = read_double_elem("T_ncdm_decay_dr", T_dncdm_default);
  ksi_                 = read_double_elem("ksi_ncdm_decay_dr", 0.0);
  m_in_eV_             = read_double_elem("m_ncdm_decay_dr", 1.0);
  deg_                 = read_double_elem("deg_ncdm_decay_dr", 1.0);

  // No file-based PSD for decay_dr (got_file_ stays false, psd_parameters_ stays empty)
  got_file_ = false;

  // Initialize quadrature (uses quadrature_strategy_, input_q_size_, qmax_, ksi_)
  InitQuadrature(settings);

  // Compute M_ from m_in_eV_: M = m / (k_B * T_ncdm * T_cmb)
  double H0 = settings.h * 1.e5 / _c_;
  if (m_in_eV_ != 0.0) {
    M_ = m_in_eV_ / _k_B_ * _eV_ / T_ / T_cmb_;
    double rho_ncdm;
    ComputeMomenta(0., nullptr, &rho_ncdm, nullptr, nullptr, nullptr);
    // For DNCDM, Omega will be set later (via SetDegAndFactor / SetDeg_from_Omega_ini)
    // so we don't compute Omega0_ here from rho.
    // We DO keep M_ from m_in_eV_.
  }
  else {
    M_ = 0.;
  }
}

// ─────────────────────────────────────────────────────────────────────────────
// ReadParameters — extract single-species data from the parameter file
// ─────────────────────────────────────────────────────────────────────────────

void NCDMBaseSpecies::ReadParameters(FileContent* pfc,
                                     int species_index,
                                     const std::string& suffix,
                                     const NcdmSettings& settings) {
  char* errmsg = error_message_;

  // 1. Determine expected number of species based on suffix
  int N_expected    = 0;
  std::string var_N = "N_ncdm" + suffix;
  int flag1, flag2 = _FALSE_;

  class_call(parser_read_int(pfc, var_N.c_str(), &N_expected, &flag1, errmsg), errmsg, errmsg);

  // Maintain backward compatibility for standard NCDM missing the suffix
  if (suffix == "_standard") {
    int int2;
    class_call(parser_read_int(pfc, "N_ncdm", &int2, &flag2, errmsg), errmsg, errmsg);
    if ((flag1 == _TRUE_) && (flag2 == _TRUE_)) {
      throw std::invalid_argument(
          "In input file, you can only enter one of N_ncdm_standard and N_ncdm, choose one");
    }
    if (flag2 == _TRUE_)
      N_expected = int2;
  }

  const double deg_ncdm_default = 1.0;
  const double T_ncdm_default   = 0.71611;
  const double ksi_ncdm_default = 0.;

  // 2. Dynamic Lambda helpers
  auto read_list_of_ints = [&](const std::string& base_name,
                               const std::string& legacy_name,
                               std::vector<int>& output,
                               int default_value) {
    int flg1, flg2 = _FALSE_, entries_read;
    std::vector<int> raw, raw2;
    std::string varname = base_name + suffix;

    class_call(readIntegerList(pfc, varname.c_str(), raw, &flg1, errmsg), errmsg, errmsg);

    // Only look for deprecated names if we are parsing standard NCDM
    if (suffix == "_standard" && !legacy_name.empty()) {
      class_call(readIntegerList(pfc, legacy_name.c_str(), raw2, &flg2, errmsg), errmsg, errmsg);
    }

    if ((flg1 == _TRUE_) && (flg2 == _TRUE_)) {
      throw std::invalid_argument("In input file, you can only enter one of " + varname + " and " +
                                  legacy_name + ", choose one");
    }
    else if ((flg1 == _TRUE_) || (flg2 == _TRUE_)) {
      entries_read = static_cast<int>((flg1 == _TRUE_) ? raw.size() : raw2.size());
      if (entries_read != N_expected) {
        throw std::invalid_argument(
            "Number of entries in " + varname +
            " does not match the expected number: " + std::to_string(N_expected));
      }
      output = (flg1 == _TRUE_) ? raw : raw2;
    }
    else {
      output.assign(N_expected, default_value);
    }
    return 0;
  };

  auto read_list_of_doubles = [&](const std::string& base_name,
                                  const std::string& legacy_name,
                                  std::vector<double>& output,
                                  double default_value) {
    int flg1, flg2 = _FALSE_, entries_read;
    std::vector<double> raw, raw2;
    std::string varname = base_name + suffix;

    class_call(readDoubleList(pfc, varname.c_str(), raw, &flg1, errmsg), errmsg, errmsg);

    if (suffix == "_standard" && !legacy_name.empty()) {
      class_call(readDoubleList(pfc, legacy_name.c_str(), raw2, &flg2, errmsg), errmsg, errmsg);
    }

    if ((flg1 == _TRUE_) && (flg2 == _TRUE_)) {
      throw std::invalid_argument("In input file, you can only enter one of " + varname + " and " +
                                  legacy_name + ", choose one");
    }
    else if ((flg1 == _TRUE_) || (flg2 == _TRUE_)) {
      entries_read = static_cast<int>((flg1 == _TRUE_) ? raw.size() : raw2.size());
      if (entries_read != N_expected) {
        throw std::invalid_argument(
            "Number of entries in " + varname +
            " does not match the expected number: " + std::to_string(N_expected));
      }
      output = (flg1 == _TRUE_) ? raw : raw2;
    }
    else {
      output.assign(N_expected, default_value);
    }
    return 0;
  };

  // 3. Extract Element [species_index] using clean, suffix-agnostic base strings
  std::vector<int> quadrature_strategies;
  read_list_of_ints("quadrature_strategy_ncdm", "Quadrature strategy", quadrature_strategies, 0);
  quadrature_strategy_ = quadrature_strategies[species_index];

  std::vector<int> input_q_sizes;
  read_list_of_ints("N_momentum_bins_ncdm", "Number of momentum bins", input_q_sizes, 5);
  input_q_size_ = input_q_sizes[species_index];

  std::vector<double> qmaxes;
  read_list_of_doubles("maximum_q_ncdm", "Maximum q", qmaxes, 15.0);
  qmax_ = qmaxes[species_index];

  std::vector<double> T_ncdm_list;
  read_list_of_doubles("T_ncdm", "T_ncdm", T_ncdm_list, T_ncdm_default);
  T_ = T_ncdm_list[species_index];

  std::vector<double> ksi_ncdm_list;
  read_list_of_doubles("ksi_ncdm", "ksi_ncdm", ksi_ncdm_list, ksi_ncdm_default);
  ksi_ = ksi_ncdm_list[species_index];

  std::vector<double> deg_ncdm_list;
  read_list_of_doubles("deg_ncdm", "deg_ncdm", deg_ncdm_list, deg_ncdm_default);
  deg_ = deg_ncdm_list[species_index];

  std::vector<double> m_ncdm_list;
  read_list_of_doubles("m_ncdm", "m_ncdm", m_ncdm_list, 0.0);
  m_in_eV_ = m_ncdm_list[species_index];

  std::vector<double> Omega0_ncdm_list;
  read_list_of_doubles("Omega_ncdm", "Omega_ncdm", Omega0_ncdm_list, 0.0);
  Omega0_ = Omega0_ncdm_list[species_index];

  std::vector<double> omega0_ncdm_list;
  read_list_of_doubles("omega_ncdm", "omega_ncdm", omega0_ncdm_list, 0.0);
  omega0_ = omega0_ncdm_list[species_index];

  // Resolve omega/Omega conflict and set defaults
  if (omega0_ != 0.0) {
    class_test(Omega0_ != 0.0,
               errmsg,
               "Nonzero values for both Omega and omega for ncdm species %d are specified!",
               species_index);
    Omega0_ = omega0_ / settings.h / settings.h;
  }
  else {
    omega0_ = Omega0_ * settings.h * settings.h;
  }

  if ((Omega0_ == 0.0) && (m_in_eV_ == 0.0)) {
    m_in_eV_ = 1.e-5;  // ultra-relativistic default
  }

  // Read PSD file flags
  std::vector<int> got_files;
  {
    int found = _FALSE_;
    // We check for "use_ncdm_psd_files" + suffix. Fall back to global if standard.
    std::string psd_flag_name = "use_ncdm_psd_files" + (suffix == "_standard" ? "" : suffix);
    class_call(readIntegerList(pfc, psd_flag_name.c_str(), got_files, &found, errmsg),
               errmsg,
               errmsg);

    if (found == _FALSE_) {
      got_files.assign(N_expected, _FALSE_);
    }
    class_test(static_cast<int>(got_files.size()) != N_expected,
               errmsg,
               "Number of entries in %s does not match expected number, %d.",
               psd_flag_name.c_str(),
               N_expected);
  }

  got_file_ = (got_files[species_index] == _TRUE_);
  if (got_file_) {
    int filenum = 0;
    for (int n = 0; n < species_index; n++)
      if (got_files[n] == _TRUE_)
        filenum++;

    int fileentries = 0;
    for (int n = 0; n < N_expected; n++)
      if (got_files[n] == _TRUE_)
        fileentries++;

    std::vector<std::string> raw_psd_files;
    int flag_files;
    std::string psd_file_name = "ncdm_psd_filenames" + (suffix == "_standard" ? "" : suffix);
    class_call(readStringList(pfc, psd_file_name.c_str(), raw_psd_files, &flag_files, errmsg),
               errmsg,
               errmsg);

    class_test(flag_files == _FALSE_,
               errmsg,
               "Input use_ncdm_psd_files is found, but no filenames found!");
    class_test(static_cast<int>(raw_psd_files.size()) != fileentries,
               errmsg,
               "Number of filenames found, %d, does not match number of _TRUE_ values in "
               "use_ncdm_psd_files, %d",
               static_cast<int>(raw_psd_files.size()),
               fileentries);
    psd_file_ = raw_psd_files[filenum];
  }

  // Read optional PSD parameters
  {
    int found              = _FALSE_;
    std::string psd_params = "ncdm_psd_parameters" + (suffix == "_standard" ? "" : suffix);
    class_call(readDoubleList(pfc, psd_params.c_str(), psd_parameters_, &found, errmsg),
               errmsg,
               errmsg);
    if (found == _FALSE_)
      psd_parameters_.clear();
  }

  // The rest remains unchanged: InitQuadrature(settings) and mass evaluations
  InitQuadrature(settings);

  double H0 = settings.h * 1.e5 / _c_;
  if (m_in_eV_ != 0.0) {
    M_ = m_in_eV_ / _k_B_ * _eV_ / T_ / T_cmb_;
    double rho_ncdm;
    ComputeMomenta(0., nullptr, &rho_ncdm, nullptr, nullptr, nullptr);
    // Overwrite behavior for pure decay products can be managed by the child class if needed,
    // but standard initialization flows normally here.
    if (Omega0_ == 0.0) {
      Omega0_ = rho_ncdm / H0 / H0;
      omega0_ = Omega0_ * settings.h * settings.h;
    }
    else {
      double fnu_factor  = H0 * H0 * Omega0_ / rho_ncdm;
      factor_           *= fnu_factor;
      deg_              *= fnu_factor;
    }
  }
  else {
    M_       = MFromOmega(H0, Omega0_, settings.tol_M_ncdm);
    m_in_eV_ = _k_B_ / _eV_ * T_ * M_ * T_cmb_;
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

  // Handle perturbation quadrature
  if (quadrature_strategy_ == qm_auto) {
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
                                    (enum quadrature_method) quadrature_strategy_,
                                    pbadist.q.data(),
                                    pbadist.tablesize,
                                    DistributionFunction,
                                    &pbadist,
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
// DistributionFunction — static callback
// ─────────────────────────────────────────────────────────────────────────────

int NCDMBaseSpecies::DistributionFunction(void* params, double q, double* f0) {
  auto* p                   = static_cast<DistributionParams*>(params);
  const NCDMBaseSpecies* sp = p->sp;
  double ksi                = sp->ksi_;

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
    // Fermi-Dirac with chemical potential
    *f0 = 1.0 / pow(2 * _PI_, 3) * (1. / (exp(q - ksi) + 1.) + 1. / (exp(q + ksi) + 1.));
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
