/** @file lensing.c Documented lensing module
 *
 * Simon Prunet and Julien Lesgourgues, 6.12.2010
 *
 * This module computes the lensed temperature and polarization
 * anisotropy power spectra \f$ C_l^{X}, P(k), ... \f$'s given the
 * unlensed temperature, polarization and lensing potential spectra.
 *
 * Follows Challinor and Lewis full-sky method, astro-ph/0502425
 *
 * The following functions can be called from other modules:
 *
 * -# lensing_init() at the beginning (but after spectra_init())
 * -# lensing_cl_at_l() at any time for computing Cl_lensed at any l
 * -# lensing_free() at the end
 */

#include "lensing_module.h"

#include <numeric>

#include "spectra_module.h"
#include "thread_pool.h"

/**
 * Anisotropy power spectra \f$ C_l\f$'s for all types, modes and initial conditions.
 * SO FAR: ONLY SCALAR
 *
 * This routine evaluates all the lensed \f$ C_l\f$'s at a given value of l by
 * picking it in the pre-computed table. When relevant, it also
 * sums over all initial conditions for each mode, and over all modes.
 *
 * This function can be called from whatever module at whatever time,
 * provided that lensing_init() has been called before, and
 * lensing_free() has not been called yet.
 *
 */

LensingModule::LensingModule(InputModulePtr input_module, SpectraModulePtr spectra_module)
    : BaseModule(std::move(input_module)), spectra_module_(std::move(spectra_module)) {
  lensing_init();
}

LensingModule::~LensingModule() {
  lensing_free();
}

std::map<std::string, std::vector<double>> LensingModule::cl_output(int lmax) const {
  ThrowRuntimeErrorIf((lmax > l_lensed_max_) || (lmax < 0),
                      "Error: lmax = %d is outside the allowed range [0, %d]\n",
                      lmax,
                      l_lensed_max_);
  std::vector<int> l_values(lmax + 1 - 2);
  std::iota(l_values.begin(), l_values.end(), 2);
  return cl_output_at_l_values(l_values);
}

std::map<std::string, std::vector<double>> LensingModule::cl_output_computed() const {
  std::vector<int> l_values;
  for (int index_l = 0; index_l < l_size_; ++index_l) {
    double ell = l_[index_l];
    if (ell > l_lensed_max_) {
      break;
    }
    l_values.push_back(ell);
  }
  return cl_output_at_l_values(l_values);
}

std::map<std::string, std::vector<double>> LensingModule::cl_output_at_l_values(
    const std::vector<int>& l_values) const {
  ThrowRuntimeErrorIf(ple->has_lensed_cls == _FALSE_,
                      "No lensed Cls was computed, adjust your inputs.\n");

  std::map<std::string, int> index_map;

  if (has_tt_)
    index_map["tt"] = index_lt_tt_;
  if (has_ee_)
    index_map["ee"] = index_lt_ee_;
  if (has_te_)
    index_map["te"] = index_lt_te_;
  if (has_bb_)
    index_map["bb"] = index_lt_bb_;
  if (has_pp_)
    index_map["pp"] = index_lt_pp_;
  if (has_tp_)
    index_map["tp"] = index_lt_tp_;

  // Create vectors for fast iteration in nested loop below.
  std::vector<std::vector<double>> data_vectors;
  std::vector<int> indices;
  std::vector<std::string> keys;

  for (const auto& element : index_map) {
    data_vectors.push_back(std::vector<double>(l_values.size() + 2, 0.0));
    indices.push_back(element.second);
    keys.push_back(element.first);
  }

  std::vector<double> cl_lensed(lt_size_);
  for (int index_l = 0; index_l < l_values.size(); ++index_l) {
    int status = lensing_cl_at_l(l_values[index_l], cl_lensed.data());
    ThrowRuntimeErrorIf(status != _SUCCESS_,
                        "Error in LensingModule::cl_output: %s",
                        error_message_);
    for (int i = 0; i < data_vectors.size(); ++i) {
      data_vectors[i][index_l + 2] = cl_lensed[indices[i]];
    }
  }
  // Now move vectors into map. We could have created the vectors inside the map directly, but that would
  // lead to many unnecessary map-lookups in the l-loop above.
  std::map<std::string, std::vector<double>> output;
  for (int i = 0; i < data_vectors.size(); ++i) {
    output[keys[i]] = std::move(data_vectors[i]);
  }
  std::vector<double> ell = {0.0, 1.0};
  ell.insert(ell.end(), l_values.begin(), l_values.end());
  output["ell"] = std::move(ell);
  return output;
}

int LensingModule::lensing_cl_at_l(int l, double* cl_lensed) const {
  class_test(l > l_lensed_max_,
             error_message_,
             "you asked for lensed Cls at l=%d, they were computed only up to l=%d, you should "
             "increase l_max_scalars or decrease the precision parameter delta_l_max",
             l,
             l_lensed_max_);

  int last_index;
  class_call(array_interpolate_spline(const_cast<double*>(l_.data()),
                                      l_size_,
                                      const_cast<double*>(cl_lens_.data()),
                                      const_cast<double*>(ddcl_lens_.data()),
                                      lt_size_,
                                      l,
                                      &last_index,
                                      cl_lensed,
                                      lt_size_,
                                      error_message_),
             error_message_,
             error_message_);

  /* set to zero for the types such that l<l_max */
  for (int index_lt = 0; index_lt < lt_size_; index_lt++)
    if ((int) l > l_max_lt_[index_lt])
      cl_lensed[index_lt] = 0.;

  return _SUCCESS_;
}

/**
 * This routine initializes the lensing structure (in particular,
 * computes table of lensed anisotropy spectra \f$ C_l^{X} \f$)
 *
 * @return the error status
 */

int LensingModule::lensing_init() {
  /** Summary: */
  /** - Define local variables */

  double* sqrt1;
  double* sqrt2;
  double* sqrt3;
  double* sqrt4;
  double* sqrt5;

  /** - check that we really want to compute at least one spectrum */

  if (ple->has_lensed_cls == _FALSE_) {
    if (ple->lensing_verbose > 0)
      printf("No lensing requested. Lensing module skipped.\n");
    return _SUCCESS_;
  }
  else {
    if (ple->lensing_verbose > 0) {
      printf("Computing lensed spectra ");
      if (ppr->accurate_lensing == _TRUE_)
        printf("(accurate mode)\n");
      else
        printf("(fast mode)\n");
    }
  }

  /** - initialize indices and allocate some of the arrays in the
      lensing structure */

  class_call(lensing_indices(), error_message_, error_message_);

  /** - put all precision variables hare; will be stored later in precision structure */
  /** - Last element in \f$ \mu \f$ will be for \f$ \mu=1 \f$, needed for sigma2.
      The rest will be chosen as roots of a Gauss-Legendre quadrature **/

  int num_mu;
  if (ppr->accurate_lensing == _TRUE_) {
    num_mu  = (l_unlensed_max_ + ppr->num_mu_minus_lmax); /* Must be even ?? CHECK */
    num_mu += num_mu % 2;                                 /* Force it to be even */
  }
  else {
    /* Integrate correlation function difference on [0,pi/16] */
    num_mu = (l_unlensed_max_ * 2) / 16;
  }
  /** - allocate array of \f$ \mu \f$ values, as well as quadrature weights */

  std::vector<double> mu(num_mu);
  /* Reserve last element of mu for mu=1, needed for sigma2 */
  mu[num_mu - 1] = 1.0;

  std::vector<double> w8(num_mu - 1);

  if (ppr->accurate_lensing == _TRUE_) {
    class_call(quadrature_gauss_legendre(mu.data(),
                                         w8.data(),
                                         num_mu - 1,
                                         ppr->tol_gauss_legendre,
                                         error_message_),
               error_message_,
               error_message_);
  }
  else { /* Crude integration on [0,pi/16]: Riemann sum on theta */

    double delta_theta = _PI_ / 16. / (double) (num_mu - 1);
    for (int index_mu = 0; index_mu < num_mu - 1; index_mu++) {
      double theta = (index_mu + 1) * delta_theta;
      mu[index_mu] = cos(theta);
      w8[index_mu] = sin(theta) * delta_theta; /* We integrate on mu */
    }
  }

  /** - Compute \f$ d^l_{mm'} (\mu) \f$*/

  int icount = 0;
  std::vector<double*> d00(num_mu);
  std::vector<double*> d11(num_mu);
  std::vector<double*> d1m1(num_mu);
  std::vector<double*> d2m2(num_mu);
  icount += 4 * num_mu * (l_unlensed_max_ + 1);

  std::vector<double*> d20, d3m1, d4m2;
  if (has_te_ == _TRUE_) {
    d20.resize(num_mu);
    d3m1.resize(num_mu);
    d4m2.resize(num_mu);
    icount += 3 * num_mu * (l_unlensed_max_ + 1);
  }

  std::vector<double*> d22, d31, d3m3, d40, d4m4;
  if (has_ee_ == _TRUE_ || has_bb_ == _TRUE_) {
    d22.resize(num_mu);
    d31.resize(num_mu);
    d3m3.resize(num_mu);
    d40.resize(num_mu);
    d4m4.resize(num_mu);
    icount += 5 * num_mu * (l_unlensed_max_ + 1);
  }

  icount += 5 * (l_unlensed_max_ + 1); /* for arrays sqrt1[l] to sqrt5[l] */

  /** - Allocate main contiguous buffer **/
  std::vector<double> buf_dxx(icount);

  icount = 0;
  for (int index_mu = 0; index_mu < num_mu; index_mu++) {
    d00[index_mu]  = &(buf_dxx[icount + (index_mu + 0 * num_mu) * (l_unlensed_max_ + 1)]);
    d11[index_mu]  = &(buf_dxx[icount + (index_mu + 1 * num_mu) * (l_unlensed_max_ + 1)]);
    d1m1[index_mu] = &(buf_dxx[icount + (index_mu + 2 * num_mu) * (l_unlensed_max_ + 1)]);
    d2m2[index_mu] = &(buf_dxx[icount + (index_mu + 3 * num_mu) * (l_unlensed_max_ + 1)]);
  }
  icount += 4 * num_mu * (l_unlensed_max_ + 1);

  if (has_te_ == _TRUE_) {
    for (int index_mu = 0; index_mu < num_mu; index_mu++) {
      d20[index_mu]  = &(buf_dxx[icount + (index_mu + 0 * num_mu) * (l_unlensed_max_ + 1)]);
      d3m1[index_mu] = &(buf_dxx[icount + (index_mu + 1 * num_mu) * (l_unlensed_max_ + 1)]);
      d4m2[index_mu] = &(buf_dxx[icount + (index_mu + 2 * num_mu) * (l_unlensed_max_ + 1)]);
    }
    icount += 3 * num_mu * (l_unlensed_max_ + 1);
  }

  if (has_ee_ == _TRUE_ || has_bb_ == _TRUE_) {
    for (int index_mu = 0; index_mu < num_mu; index_mu++) {
      d22[index_mu]  = &(buf_dxx[icount + (index_mu + 0 * num_mu) * (l_unlensed_max_ + 1)]);
      d31[index_mu]  = &(buf_dxx[icount + (index_mu + 1 * num_mu) * (l_unlensed_max_ + 1)]);
      d3m3[index_mu] = &(buf_dxx[icount + (index_mu + 2 * num_mu) * (l_unlensed_max_ + 1)]);
      d40[index_mu]  = &(buf_dxx[icount + (index_mu + 3 * num_mu) * (l_unlensed_max_ + 1)]);
      d4m4[index_mu] = &(buf_dxx[icount + (index_mu + 4 * num_mu) * (l_unlensed_max_ + 1)]);
    }
    icount += 5 * num_mu * (l_unlensed_max_ + 1);
  }

  sqrt1   = &(buf_dxx[icount]);
  icount += l_unlensed_max_ + 1;
  sqrt2   = &(buf_dxx[icount]);
  icount += l_unlensed_max_ + 1;
  sqrt3   = &(buf_dxx[icount]);
  icount += l_unlensed_max_ + 1;
  sqrt4   = &(buf_dxx[icount]);
  icount += l_unlensed_max_ + 1;
  sqrt5   = &(buf_dxx[icount]);
  Tools::TaskSystem task_system(pba->number_of_threads);
  std::vector<std::future<int>> dXX_tasks;

  dXX_tasks.push_back(task_system.AsyncTask(
      [&]() { return lensing_d00(mu.data(), num_mu, l_unlensed_max_, d00.data()); }));

  dXX_tasks.push_back(task_system.AsyncTask(
      [&]() { return lensing_d11(mu.data(), num_mu, l_unlensed_max_, d11.data()); }));

  dXX_tasks.push_back(task_system.AsyncTask(
      [&]() { return lensing_d1m1(mu.data(), num_mu, l_unlensed_max_, d1m1.data()); }));

  dXX_tasks.push_back(task_system.AsyncTask(
      [&]() { return lensing_d2m2(mu.data(), num_mu, l_unlensed_max_, d2m2.data()); }));

  if (has_te_ == _TRUE_) {
    dXX_tasks.push_back(task_system.AsyncTask(
        [&]() { return lensing_d20(mu.data(), num_mu, l_unlensed_max_, d20.data()); }));
    dXX_tasks.push_back(task_system.AsyncTask(
        [&]() { return lensing_d3m1(mu.data(), num_mu, l_unlensed_max_, d3m1.data()); }));
    dXX_tasks.push_back(task_system.AsyncTask(
        [&]() { return lensing_d4m2(mu.data(), num_mu, l_unlensed_max_, d4m2.data()); }));
  }

  if (has_ee_ == _TRUE_ || has_bb_ == _TRUE_) {
    dXX_tasks.push_back(task_system.AsyncTask(
        [&]() { return lensing_d22(mu.data(), num_mu, l_unlensed_max_, d22.data()); }));
    dXX_tasks.push_back(task_system.AsyncTask(
        [&]() { return lensing_d31(mu.data(), num_mu, l_unlensed_max_, d31.data()); }));
    dXX_tasks.push_back(task_system.AsyncTask(
        [&]() { return lensing_d3m3(mu.data(), num_mu, l_unlensed_max_, d3m3.data()); }));
    dXX_tasks.push_back(task_system.AsyncTask(
        [&]() { return lensing_d40(mu.data(), num_mu, l_unlensed_max_, d40.data()); }));
    dXX_tasks.push_back(task_system.AsyncTask(
        [&]() { return lensing_d4m4(mu.data(), num_mu, l_unlensed_max_, d4m4.data()); }));
  }

  for (auto& task : dXX_tasks) {
    class_call(task.get(), error_message_, error_message_);
  }

  /** - compute \f$ Cgl(\mu)\f$, \f$ Cgl2(\mu) \f$ and sigma2(\f$\mu\f$) */

  std::vector<double> Cgl(num_mu);

  std::vector<double> Cgl2(num_mu);

  std::vector<double> sigma2(num_mu - 1);

  std::vector<double> cl_unlensed(spectra_module_->ct_size_);

  /** - Locally store unlensed temperature \f$ cl_{tt}\f$ and potential \f$ cl_{pp}\f$ spectra **/
  std::vector<double> cl_tt(l_unlensed_max_ + 1);
  std::vector<double> cl_te;
  if (has_te_ == _TRUE_) {
    cl_te.resize(l_unlensed_max_ + 1);
  }
  std::vector<double> cl_ee, cl_bb;
  if (has_ee_ == _TRUE_ || has_bb_ == _TRUE_) {
    cl_ee.resize(l_unlensed_max_ + 1);
    cl_bb.resize(l_unlensed_max_ + 1);
  }
  std::vector<double> cl_pp(l_unlensed_max_ + 1);

  std::vector<std::vector<double>> cl_md_ic_storage(spectra_module_->md_size_);
  std::vector<std::vector<double>> cl_md_storage(spectra_module_->md_size_);
  std::vector<double*> cl_md(spectra_module_->md_size_, nullptr);
  std::vector<double*> cl_md_ic(spectra_module_->md_size_, nullptr);

  for (int index_md = 0; index_md < spectra_module_->md_size_; index_md++) {
    if (spectra_module_->md_size_ > 1) {
      cl_md_storage[index_md].resize(spectra_module_->ct_size_);
      cl_md[index_md] = cl_md_storage[index_md].data();
    }

    if (spectra_module_->ic_size_[index_md] > 1) {
      cl_md_ic_storage[index_md].resize(spectra_module_->ic_ic_size_[index_md] *
                                        spectra_module_->ct_size_);
      cl_md_ic[index_md] = cl_md_ic_storage[index_md].data();
    }
  }

  for (int l = 2; l <= l_unlensed_max_; l++) {
    class_call(spectra_module_->spectra_cl_at_l(l,
                                                cl_unlensed.data(),
                                                cl_md.data(),
                                                cl_md_ic.data()),
               psp->error_message,
               error_message_);
    cl_tt[l] = cl_unlensed[index_lt_tt_];
    cl_pp[l] = cl_unlensed[index_lt_pp_];
    if (has_te_ == _TRUE_) {
      cl_te[l] = cl_unlensed[index_lt_te_];
    }
    if (has_ee_ == _TRUE_ || has_bb_ == _TRUE_) {
      cl_ee[l] = cl_unlensed[index_lt_ee_];
      cl_bb[l] = cl_unlensed[index_lt_bb_];
    }
  }

  /** - Compute sigma2\f$(\mu)\f$ and Cgl2(\f$\mu\f$) **/
  for (int index_mu = 0; index_mu < num_mu; index_mu++) {
    Cgl[index_mu]  = 0;
    Cgl2[index_mu] = 0;

    for (int l = 2; l <= l_unlensed_max_; l++) {
      Cgl[index_mu] += (2. * l + 1.) * l * (l + 1.) * cl_pp[l] * d11[index_mu][l];

      Cgl2[index_mu] += (2. * l + 1.) * l * (l + 1.) * cl_pp[l] * d1m1[index_mu][l];
    }

    Cgl[index_mu]  /= 4. * _PI_;
    Cgl2[index_mu] /= 4. * _PI_;
  }

  for (int index_mu = 0; index_mu < num_mu - 1; index_mu++) {
    /* Cgl(1.0) - Cgl(mu) */
    sigma2[index_mu] = Cgl[num_mu - 1] - Cgl[index_mu];
  }

  /** - compute ksi, ksi+, ksi-, ksiX */

  /** - --> ksi is for TT **/
  std::vector<double> ksi;
  if (has_tt_ == _TRUE_) {
    ksi.assign(num_mu - 1, 0.0);
  }

  /** - --> ksiX is for TE **/
  std::vector<double> ksiX;
  if (has_te_ == _TRUE_) {
    ksiX.assign(num_mu - 1, 0.0);
  }

  /** - --> ksip, ksim for EE, BB **/
  std::vector<double> ksip, ksim;
  if (has_ee_ == _TRUE_ || has_bb_ == _TRUE_) {
    ksip.assign(num_mu - 1, 0.0);
    ksim.assign(num_mu - 1, 0.0);
  }

  for (int l = 2; l <= l_unlensed_max_; l++) {
    double ll = (double) l;
    sqrt1[l]  = sqrt((ll + 2) * (ll + 1) * ll * (ll - 1));
    sqrt2[l]  = sqrt((ll + 2) * (ll - 1));
    sqrt3[l]  = sqrt((ll + 3) * (ll - 2));
    sqrt4[l]  = sqrt((ll + 4) * (ll + 3) * (ll - 2.) * (ll - 3));
    sqrt5[l]  = sqrt(ll * (ll + 1));
  }

  for (int index_mu = 0; index_mu < num_mu - 1; index_mu++) {
    for (int l = 2; l <= l_unlensed_max_; l++) {
      double ll = (double) l;

      double fac  = ll * (ll + 1) / 4.;
      double fac1 = (2 * ll + 1) / (4. * _PI_);

      /* In the following we will keep terms of the form (sigma2)^k*(Cgl2)^m
         with k+m <= 2 */

      double X_000  = exp(-fac * sigma2[index_mu]);
      double X_p000 = -fac * X_000;
      /* X_220 = 0.25*sqrt1[l] * exp(-(fac-0.5)*sigma2[index_mu]); */
      double X_220 = 0.25 * sqrt1[l] * X_000; /* Order 0 */
      /* next 5 lines useless, but avoid compiler warning 'may be used uninitialized' */
      double X_242  = 0.;
      double X_132  = 0.;
      double X_121  = 0.;
      double X_p022 = 0.;
      double X_022  = 0.;

      if (has_te_ == _TRUE_ || has_ee_ == _TRUE_ || has_bb_ == _TRUE_) {
        /* X_022 = exp(-(fac-1.)*sigma2[index_mu]); */
        X_022  = X_000 * (1 + sigma2[index_mu] * (1 + 0.5 * sigma2[index_mu])); /* Order 2 */
        X_p022 = -(fac - 1.) * X_022; /* Old versions were missing the
        minus sign in this line, which introduced a very small error
        on the high-l C_l^TE lensed spectrum [credits for bug fix:
        Selim Hotinli] */

        /* X_242 = 0.25*sqrt4[l] * exp(-(fac-5./2.)*sigma2[index_mu]); */
        X_242 = 0.25 * sqrt4[l] * X_000; /* Order 0 */
        if (has_ee_ == _TRUE_ || has_bb_ == _TRUE_) {
          /* X_121 = - 0.5*sqrt2[l] * exp(-(fac-2./3.)*sigma2[index_mu]);
             X_132 = - 0.5*sqrt3[l] * exp(-(fac-5./3.)*sigma2[index_mu]); */
          X_121 = -0.5 * sqrt2[l] * X_000 * (1 + 2. / 3. * sigma2[index_mu]); /* Order 1 */
          X_132 = -0.5 * sqrt3[l] * X_000 * (1 + 5. / 3. * sigma2[index_mu]); /* Order 1 */
        }
      }

      if (has_tt_ == _TRUE_) {
        double res = fac1 * cl_tt[l];

        double lens = (X_000 * X_000 * d00[index_mu][l] +
                       X_p000 * X_p000 * d1m1[index_mu][l] * Cgl2[index_mu] * 8. / (ll * (ll + 1)) +
                       (X_p000 * X_p000 * d00[index_mu][l] + X_220 * X_220 * d2m2[index_mu][l]) *
                           Cgl2[index_mu] * Cgl2[index_mu]);
        if (ppr->accurate_lensing == _FALSE_) {
          /* Remove unlensed correlation function */
          lens -= d00[index_mu][l];
        }
        res           *= lens;
        ksi[index_mu] += res;
      }

      if (has_te_ == _TRUE_) {
        double resX = fac1 * cl_te[l];

        double lens = (X_022 * X_000 * d20[index_mu][l] +
                       Cgl2[index_mu] * 2. * X_p000 / sqrt5[l] *
                           (X_121 * d11[index_mu][l] + X_132 * d3m1[index_mu][l]) +
                       0.5 * Cgl2[index_mu] * Cgl2[index_mu] *
                           ((2. * X_p022 * X_p000 + X_220 * X_220) * d20[index_mu][l] +
                            X_220 * X_242 * d4m2[index_mu][l]));
        if (ppr->accurate_lensing == _FALSE_) {
          lens -= d20[index_mu][l];
        }
        resX           *= lens;
        ksiX[index_mu] += resX;
      }

      if (has_ee_ == _TRUE_ || has_bb_ == _TRUE_) {
        double resp = fac1 * (cl_ee[l] + cl_bb[l]);
        double resm = fac1 * (cl_ee[l] - cl_bb[l]);

        double lensp =
            (X_022 * X_022 * d22[index_mu][l] +
             2. * Cgl2[index_mu] * X_132 * X_121 * d31[index_mu][l] +
             Cgl2[index_mu] * Cgl2[index_mu] *
                 (X_p022 * X_p022 * d22[index_mu][l] + X_242 * X_220 * d40[index_mu][l]));

        double lensm = (X_022 * X_022 * d2m2[index_mu][l] +
                        Cgl2[index_mu] * (X_121 * X_121 * d1m1[index_mu][l] +
                                          X_132 * X_132 * d3m3[index_mu][l]) +
                        0.5 * Cgl2[index_mu] * Cgl2[index_mu] *
                            (2. * X_p022 * X_p022 * d2m2[index_mu][l] +
                             X_220 * X_220 * d00[index_mu][l] + X_242 * X_242 * d4m4[index_mu][l]));
        if (ppr->accurate_lensing == _FALSE_) {
          lensp -= d22[index_mu][l];
          lensm -= d2m2[index_mu][l];
        }
        resp           *= lensp;
        resm           *= lensm;
        ksip[index_mu] += resp;
        ksim[index_mu] += resm;
      }
    }
  }

  /** - compute lensed \f$ C_l\f$'s by integration */
  if (has_tt_ == _TRUE_) {
    class_call(lensing_lensed_cl_tt(ksi.data(), d00.data(), w8.data(), num_mu - 1),
               error_message_,
               error_message_);
    if (ppr->accurate_lensing == _FALSE_) {
      class_call(lensing_addback_cl_tt(cl_tt.data()), error_message_, error_message_);
    }
  }

  if (has_te_ == _TRUE_) {
    class_call(lensing_lensed_cl_te(ksiX.data(), d20.data(), w8.data(), num_mu - 1),
               error_message_,
               error_message_);
    if (ppr->accurate_lensing == _FALSE_) {
      class_call(lensing_addback_cl_te(cl_te.data()), error_message_, error_message_);
    }
  }

  if (has_ee_ == _TRUE_ || has_bb_ == _TRUE_) {
    class_call(lensing_lensed_cl_ee_bb(ksip.data(),
                                       ksim.data(),
                                       d22.data(),
                                       d2m2.data(),
                                       w8.data(),
                                       num_mu - 1),
               error_message_,
               error_message_);
    if (ppr->accurate_lensing == _FALSE_) {
      class_call(lensing_addback_cl_ee_bb(cl_ee.data(), cl_bb.data()),
                 error_message_,
                 error_message_);
    }
  }

  /** - spline computed \f$ C_l\f$'s in view of interpolation */

  class_call(array_spline_table_lines(l_.data(),
                                      l_size_,
                                      cl_lens_.data(),
                                      lt_size_,
                                      ddcl_lens_.data(),
                                      _SPLINE_EST_DERIV_,
                                      error_message_),
             error_message_,
             error_message_);

  /** - Exit **/

  return _SUCCESS_;
}

/**
 * This routine frees all the memory space allocated by lensing_init().
 *
 * To be called at the end of each run, only when no further calls to
 * lensing_cl_at_l() are needed.
 *
 * @return the error status
 */

int LensingModule::lensing_free() {
  return _SUCCESS_;
}

/**
 * This routine defines indices and allocates tables in the lensing structure
 *
 * @return the error status
 */

int LensingModule::lensing_indices() {
  /* indices of all Cl types (lensed and unlensed) */

  if (spectra_module_->has_tt_ == _TRUE_) {
    has_tt_      = _TRUE_;
    index_lt_tt_ = spectra_module_->index_ct_tt_;
  }
  else {
    has_tt_ = _FALSE_;
  }

  if (spectra_module_->has_ee_ == _TRUE_) {
    has_ee_      = _TRUE_;
    index_lt_ee_ = spectra_module_->index_ct_ee_;
  }
  else {
    has_ee_ = _FALSE_;
  }

  if (spectra_module_->has_te_ == _TRUE_) {
    has_te_      = _TRUE_;
    index_lt_te_ = spectra_module_->index_ct_te_;
  }
  else {
    has_te_ = _FALSE_;
  }

  if (spectra_module_->has_bb_ == _TRUE_) {
    has_bb_      = _TRUE_;
    index_lt_bb_ = spectra_module_->index_ct_bb_;
  }
  else {
    has_bb_ = _FALSE_;
  }

  if (spectra_module_->has_pp_ == _TRUE_) {
    has_pp_      = _TRUE_;
    index_lt_pp_ = spectra_module_->index_ct_pp_;
  }
  else {
    has_pp_ = _FALSE_;
  }

  if (spectra_module_->has_tp_ == _TRUE_) {
    has_tp_      = _TRUE_;
    index_lt_tp_ = spectra_module_->index_ct_tp_;
  }
  else {
    has_tp_ = _FALSE_;
  }

  if (spectra_module_->has_dd_ == _TRUE_) {
    has_dd_      = _TRUE_;
    index_lt_dd_ = spectra_module_->index_ct_dd_;
  }
  else {
    has_dd_ = _FALSE_;
  }

  if (spectra_module_->has_td_ == _TRUE_) {
    has_td_      = _TRUE_;
    index_lt_td_ = spectra_module_->index_ct_td_;
  }
  else {
    has_td_ = _FALSE_;
  }

  if (spectra_module_->has_ll_ == _TRUE_) {
    has_ll_      = _TRUE_;
    index_lt_ll_ = spectra_module_->index_ct_ll_;
  }
  else {
    has_ll_ = _FALSE_;
  }

  if (spectra_module_->has_tl_ == _TRUE_) {
    has_tl_      = _TRUE_;
    index_lt_tl_ = spectra_module_->index_ct_tl_;
  }
  else {
    has_tl_ = _FALSE_;
  }

  lt_size_ = spectra_module_->ct_size_;

  /* number of multipoles */

  l_unlensed_max_ = spectra_module_->l_max_tot_;

  l_lensed_max_ = l_unlensed_max_ - ppr->delta_l_max;

  int index_l;
  for (index_l = 0;
       (index_l < spectra_module_->l_size_max_) && (spectra_module_->l_[index_l] <= l_lensed_max_);
       index_l++)
    ;

  if (index_l < spectra_module_->l_size_max_)
    index_l++; /* one more point in order to be able to interpolate till l_lensed_max_ */

  l_size_ = index_l + 1;

  l_.resize(l_size_);

  for (index_l = 0; index_l < l_size_; index_l++) {
    l_[index_l] = spectra_module_->l_[index_l];
  }

  /* allocate table where results will be stored */

  cl_lens_.resize(l_size_ * lt_size_);

  ddcl_lens_.resize(l_size_ * lt_size_);

  /* fill with unlensed cls */

  std::vector<std::vector<double>> cl_md_ic_storage(spectra_module_->md_size_);
  std::vector<std::vector<double>> cl_md_storage(spectra_module_->md_size_);
  std::vector<double*> cl_md_ptrs(spectra_module_->md_size_, nullptr);
  std::vector<double*> cl_md_ic_ptrs(spectra_module_->md_size_, nullptr);

  for (int index_md = 0; index_md < spectra_module_->md_size_; index_md++) {
    if (spectra_module_->md_size_ > 1) {
      cl_md_storage[index_md].resize(spectra_module_->ct_size_);
      cl_md_ptrs[index_md] = cl_md_storage[index_md].data();
    }

    if (spectra_module_->ic_size_[index_md] > 1) {
      cl_md_ic_storage[index_md].resize(spectra_module_->ic_ic_size_[index_md] *
                                        spectra_module_->ct_size_);
      cl_md_ic_ptrs[index_md] = cl_md_ic_storage[index_md].data();
    }
  }

  for (index_l = 0; index_l < l_size_; index_l++) {
    class_call(spectra_module_->spectra_cl_at_l(l_[index_l],
                                                &(cl_lens_[index_l * lt_size_]),
                                                cl_md_ptrs.data(),
                                                cl_md_ic_ptrs.data()),
               psp->error_message,
               error_message_);
  }

  /* we want to output Cl_lensed up to the same l_max as Cl_unlensed
     (even if a number delta_l_max of extra values of l have been used
     internally for more accurate results). Notable exception to the
     above rule: ClBB_lensed(scalars) must be outputed at least up to the same l_max as
     ClEE_unlensed(scalars) (since ClBB_unlensed is null for scalars)
  */

  l_max_lt_.resize(lt_size_);
  for (int index_lt = 0; index_lt < lt_size_; index_lt++) {
    l_max_lt_[index_lt] = 0;
    for (int index_md = 0; index_md < spectra_module_->md_size_; index_md++) {
      l_max_lt_[index_lt] = MAX(l_max_lt_[index_lt],
                                spectra_module_->l_max_ct_[index_md][index_lt]);

      if ((has_bb_ == _TRUE_) && (has_ee_ == _TRUE_) && (index_lt == index_lt_bb_)) {
        l_max_lt_[index_lt] = MAX(l_max_lt_[index_lt],
                                  spectra_module_->l_max_ct_[index_md][index_lt_ee_]);
      }
    }
  }

  return _SUCCESS_;
}

/**
 * This routine computes the lensed power spectra by Gaussian quadrature
 *
 * @param ksi  Input: Lensed correlation function (ksi[index_mu])
 * @param d00  Input: Legendre polynomials (\f$ d^l_{00}\f$[l][index_mu])
 * @param w8   Input: Legendre quadrature weights (w8[index_mu])
 * @param nmu  Input: Number of quadrature points (0<=index_mu<=nmu)
 * @return the error status
 */

int LensingModule::lensing_lensed_cl_tt(double* ksi, double** d00, double* w8, int nmu) {
  /** Integration by Gauss-Legendre quadrature. **/

  for (int index_l = 0; index_l < l_size_; index_l++) {
    double cle = 0;
    for (int imu = 0; imu < nmu; imu++) {
      cle += ksi[imu] * d00[imu][(int) l_[index_l]] * w8[imu]; /* loop could be optimized */
    }
    cl_lens_[index_l * lt_size_ + index_lt_tt_] = cle * 2.0 * _PI_;
  }

  return _SUCCESS_;
}

/**
 * This routine adds back the unlensed \f$ cl_{tt}\f$ power spectrum
 * Used in case of fast (and BB inaccurate) integration of
 * correlation functions.
 *
 * @return the error status
 */

int LensingModule::lensing_addback_cl_tt(double* cl_tt) {
  for (int index_l = 0; index_l < l_size_; index_l++) {
    int l                                        = (int) l_[index_l];
    cl_lens_[index_l * lt_size_ + index_lt_tt_] += cl_tt[l];
  }
  return _SUCCESS_;
}

/**
 * This routine computes the lensed power spectra by Gaussian quadrature
 *
 * @param ksiX Input: Lensed correlation function (ksiX[index_mu])
 * @param d20  Input: Wigner d-function (\f$ d^l_{20}\f$[l][index_mu])
 * @param w8   Input: Legendre quadrature weights (w8[index_mu])
 * @param nmu  Input: Number of quadrature points (0<=index_mu<=nmu)
 * @return the error status
 */

int LensingModule::lensing_lensed_cl_te(double* ksiX, double** d20, double* w8, int nmu) {
  /** Integration by Gauss-Legendre quadrature. **/

  for (int index_l = 0; index_l < l_size_; index_l++) {
    double clte = 0;
    for (int imu = 0; imu < nmu; imu++) {
      clte += ksiX[imu] * d20[imu][(int) l_[index_l]] * w8[imu]; /* loop could be optimized */
    }
    cl_lens_[index_l * lt_size_ + index_lt_te_] = clte * 2.0 * _PI_;
  }

  return _SUCCESS_;
}

/**
 * This routine adds back the unlensed \f$ cl_{te}\f$ power spectrum
 * Used in case of fast (and BB inaccurate) integration of
 * correlation functions.
 *
 * @param cl_te Input: Array of unlensed power spectrum
 * @return the error status
 */

int LensingModule::lensing_addback_cl_te(double* cl_te) {
  for (int index_l = 0; index_l < l_size_; index_l++) {
    int l                                        = (int) l_[index_l];
    cl_lens_[index_l * lt_size_ + index_lt_te_] += cl_te[l];
  }
  return _SUCCESS_;
}

/**
 * This routine computes the lensed power spectra by Gaussian quadrature
 *
 * @param ksip Input: Lensed correlation function (ksi+[index_mu])
 * @param ksim Input: Lensed correlation function (ksi-[index_mu])
 * @param d22  Input: Wigner d-function (\f$ d^l_{22}\f$[l][index_mu])
 * @param d2m2 Input: Wigner d-function (\f$ d^l_{2-2}\f$[l][index_mu])
 * @param w8   Input: Legendre quadrature weights (w8[index_mu])
 * @param nmu  Input: Number of quadrature points (0<=index_mu<=nmu)
 * @return the error status
 */

int LensingModule::lensing_lensed_cl_ee_bb(
    double* ksip, double* ksim, double** d22, double** d2m2, double* w8, int nmu) {
  /** Integration by Gauss-Legendre quadrature. **/

  for (int index_l = 0; index_l < l_size_; index_l++) {
    double clp = 0, clm = 0;
    for (int imu = 0; imu < nmu; imu++) {
      clp += ksip[imu] * d22[imu][(int) l_[index_l]] * w8[imu];  /* loop could be optimized */
      clm += ksim[imu] * d2m2[imu][(int) l_[index_l]] * w8[imu]; /* loop could be optimized */
    }
    cl_lens_[index_l * lt_size_ + index_lt_ee_] = (clp + clm) * _PI_;
    cl_lens_[index_l * lt_size_ + index_lt_bb_] = (clp - clm) * _PI_;
  }

  return _SUCCESS_;
}

/**
 * This routine adds back the unlensed \f$ cl_{ee}\f$, \f$ cl_{bb}\f$ power spectra
 * Used in case of fast (and BB inaccurate) integration of
 * correlation functions.
 *
 * @param cl_ee Input: Array of unlensed power spectrum
 * @param cl_bb Input: Array of unlensed power spectrum
 * @return the error status
 */

int LensingModule::lensing_addback_cl_ee_bb(double* cl_ee, double* cl_bb) {
  for (int index_l = 0; index_l < l_size_; index_l++) {
    int l                                        = (int) l_[index_l];
    cl_lens_[index_l * lt_size_ + index_lt_ee_] += cl_ee[l];
    cl_lens_[index_l * lt_size_ + index_lt_bb_] += cl_bb[l];
  }
  return _SUCCESS_;
}

/**
 * This routine computes the d00 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d00    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int LensingModule::lensing_d00(double* mu, int num_mu, int lmax, double** d00) {
  std::vector<double> fac1(lmax), fac2(lmax), fac3(lmax);

  for (int l = 1; l < lmax; l++) {
    double ll = (double) l;
    fac1[l]   = sqrt((2 * ll + 3) / (2 * ll + 1)) * (2 * ll + 1) / (ll + 1);
    fac2[l]   = sqrt((2 * ll + 3) / (2 * ll - 1)) * ll / (ll + 1);
    fac3[l]   = sqrt(2. / (2 * ll + 3));
  }

  for (int index_mu = 0; index_mu < num_mu; index_mu++) {
    double dlm1      = 1.0 / sqrt(2.); /* l=0 */
    d00[index_mu][0] = dlm1 * sqrt(2.);
    double dl        = mu[index_mu] * sqrt(3. / 2.); /*l=1*/
    d00[index_mu][1] = dl * sqrt(2. / 3.);
    for (int l = 1; l < lmax; l++) {
      /* sqrt((2l+1)/2)*d00 recurrence, supposed to be more stable */
      double dlp1          = fac1[l] * mu[index_mu] * dl - fac2[l] * dlm1;
      d00[index_mu][l + 1] = dlp1 * fac3[l];
      dlm1                 = dl;
      dl                   = dlp1;
    }
  }
  return _SUCCESS_;
}

/**
 * This routine computes the d11 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d11    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int LensingModule::lensing_d11(double* mu, int num_mu, int lmax, double** d11) {
  std::vector<double> fac1(lmax), fac2(lmax), fac3(lmax), fac4(lmax);
  for (int l = 2; l < lmax; l++) {
    double ll = (double) l;
    fac1[l]   = sqrt((2 * ll + 3) / (2 * ll + 1)) * (ll + 1) * (2 * ll + 1) / (ll * (ll + 2));
    fac2[l]   = 1.0 / (ll * (ll + 1.));
    fac3[l] = sqrt((2 * ll + 3) / (2 * ll - 1)) * (ll - 1) * (ll + 1) / (ll * (ll + 2)) * (ll + 1) /
              ll;
    fac4[l] = sqrt(2. / (2 * ll + 3));
  }

  for (int index_mu = 0; index_mu < num_mu; index_mu++) {
    d11[index_mu][0] = 0;
    double dlm1      = (1.0 + mu[index_mu]) / 2. * sqrt(3. / 2.); /*l=1*/
    d11[index_mu][1] = dlm1 * sqrt(2. / 3.);
    double dl = (1.0 + mu[index_mu]) / 2. * (2.0 * mu[index_mu] - 1.0) * sqrt(5. / 2.); /*l=2*/
    d11[index_mu][2] = dl * sqrt(2. / 5.);
    for (int l = 2; l < lmax; l++) {
      /* sqrt((2l+1)/2)*d11 recurrence, supposed to be more stable */
      double dlp1          = fac1[l] * (mu[index_mu] - fac2[l]) * dl - fac3[l] * dlm1;
      d11[index_mu][l + 1] = dlp1 * fac4[l];
      dlm1                 = dl;
      dl                   = dlp1;
    }
  }
  return _SUCCESS_;
}

/**
 * This routine computes the d1m1 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d1m1    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int LensingModule::lensing_d1m1(double* mu, int num_mu, int lmax, double** d1m1) {
  std::vector<double> fac1(lmax), fac2(lmax), fac3(lmax), fac4(lmax);
  for (int l = 2; l < lmax; l++) {
    double ll = (double) l;
    fac1[l]   = sqrt((2 * ll + 3) / (2 * ll + 1)) * (ll + 1) * (2 * ll + 1) / (ll * (ll + 2));
    fac2[l]   = 1.0 / (ll * (ll + 1.));
    fac3[l] = sqrt((2 * ll + 3) / (2 * ll - 1)) * (ll - 1) * (ll + 1) / (ll * (ll + 2)) * (ll + 1) /
              ll;
    fac4[l] = sqrt(2. / (2 * ll + 3));
  }

  for (int index_mu = 0; index_mu < num_mu; index_mu++) {
    d1m1[index_mu][0] = 0;
    double dlm1       = (1.0 - mu[index_mu]) / 2. * sqrt(3. / 2.); /*l=1*/
    d1m1[index_mu][1] = dlm1 * sqrt(2. / 3.);
    double dl = (1.0 - mu[index_mu]) / 2. * (2.0 * mu[index_mu] + 1.0) * sqrt(5. / 2.); /*l=2*/
    d1m1[index_mu][2] = dl * sqrt(2. / 5.);
    for (int l = 2; l < lmax; l++) {
      /* sqrt((2l+1)/2)*d1m1 recurrence, supposed to be more stable */
      double dlp1           = fac1[l] * (mu[index_mu] + fac2[l]) * dl - fac3[l] * dlm1;
      d1m1[index_mu][l + 1] = dlp1 * fac4[l];
      dlm1                  = dl;
      dl                    = dlp1;
    }
  }
  return _SUCCESS_;
}

/**
 * This routine computes the d2m2 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d2m2   Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int LensingModule::lensing_d2m2(double* mu, int num_mu, int lmax, double** d2m2) {
  std::vector<double> fac1(lmax), fac2(lmax), fac3(lmax), fac4(lmax);
  for (int l = 2; l < lmax; l++) {
    double ll = (double) l;
    fac1[l]   = sqrt((2 * ll + 3) / (2 * ll + 1)) * (ll + 1) * (2 * ll + 1) / ((ll - 1) * (ll + 3));
    fac2[l]   = 4.0 / (ll * (ll + 1));
    fac3[l]   = sqrt((2 * ll + 3) / (2 * ll - 1)) * (ll - 2) * (ll + 2) / ((ll - 1) * (ll + 3)) *
                (ll + 1) / ll;
    fac4[l]   = sqrt(2. / (2 * ll + 3));
  }

  for (int index_mu = 0; index_mu < num_mu; index_mu++) {
    d2m2[index_mu][0] = 0;
    double dlm1       = 0.; /*l=1*/
    d2m2[index_mu][1] = 0;
    double dl         = (1.0 - mu[index_mu]) * (1.0 - mu[index_mu]) / 4. * sqrt(5. / 2.); /*l=2*/
    d2m2[index_mu][2] = dl * sqrt(2. / 5.);
    for (int l = 2; l < lmax; l++) {
      /* sqrt((2l+1)/2)*d2m2 recurrence, supposed to be more stable */
      double dlp1           = fac1[l] * (mu[index_mu] + fac2[l]) * dl - fac3[l] * dlm1;
      d2m2[index_mu][l + 1] = dlp1 * fac4[l];
      dlm1                  = dl;
      dl                    = dlp1;
    }
  }
  return _SUCCESS_;
}

/**
 * This routine computes the d22 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d22    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int LensingModule::lensing_d22(double* mu, int num_mu, int lmax, double** d22) {
  std::vector<double> fac1(lmax), fac2(lmax), fac3(lmax), fac4(lmax);
  for (int l = 2; l < lmax; l++) {
    double ll = (double) l;
    fac1[l]   = sqrt((2 * ll + 3) / (2 * ll + 1)) * (ll + 1) * (2 * ll + 1) / ((ll - 1) * (ll + 3));
    fac2[l]   = 4.0 / (ll * (ll + 1));
    fac3[l]   = sqrt((2 * ll + 3) / (2 * ll - 1)) * (ll - 2) * (ll + 2) / ((ll - 1) * (ll + 3)) *
                (ll + 1) / ll;
    fac4[l]   = sqrt(2. / (2 * ll + 3));
  }

  for (int index_mu = 0; index_mu < num_mu; index_mu++) {
    d22[index_mu][0] = 0;
    double dlm1      = 0.; /*l=1*/
    d22[index_mu][1] = 0;
    double dl        = (1.0 + mu[index_mu]) * (1.0 + mu[index_mu]) / 4. * sqrt(5. / 2.); /*l=2*/
    d22[index_mu][2] = dl * sqrt(2. / 5.);
    for (int l = 2; l < lmax; l++) {
      /* sqrt((2l+1)/2)*d22 recurrence, supposed to be more stable */
      double dlp1          = fac1[l] * (mu[index_mu] - fac2[l]) * dl - fac3[l] * dlm1;
      d22[index_mu][l + 1] = dlp1 * fac4[l];
      dlm1                 = dl;
      dl                   = dlp1;
    }
  }
  return _SUCCESS_;
}

/**
 * This routine computes the d20 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d20    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int LensingModule::lensing_d20(double* mu, int num_mu, int lmax, double** d20) {
  std::vector<double> fac1(lmax), fac3(lmax), fac4(lmax);
  for (int l = 2; l < lmax; l++) {
    double ll = (double) l;
    fac1[l]   = sqrt((2 * ll + 3) * (2 * ll + 1) / ((ll - 1) * (ll + 3)));
    fac3[l]   = sqrt((2 * ll + 3) * (ll - 2) * (ll + 2) / ((2 * ll - 1) * (ll - 1) * (ll + 3)));
    fac4[l]   = sqrt(2. / (2 * ll + 3));
  }

  for (int index_mu = 0; index_mu < num_mu; index_mu++) {
    d20[index_mu][0] = 0;
    double dlm1      = 0.; /*l=1*/
    d20[index_mu][1] = 0;
    double dl        = sqrt(15.) / 4. * (1 - mu[index_mu] * mu[index_mu]); /*l=2*/
    d20[index_mu][2] = dl * sqrt(2. / 5.);
    for (int l = 2; l < lmax; l++) {
      /* sqrt((2l+1)/2)*d22 recurrence, supposed to be more stable */
      double dlp1          = fac1[l] * mu[index_mu] * dl - fac3[l] * dlm1;
      d20[index_mu][l + 1] = dlp1 * fac4[l];
      dlm1                 = dl;
      dl                   = dlp1;
    }
  }
  return _SUCCESS_;
}

/**
 * This routine computes the d31 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d31    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int LensingModule::lensing_d31(double* mu, int num_mu, int lmax, double** d31) {
  std::vector<double> fac1(lmax), fac2(lmax), fac3(lmax), fac4(lmax);
  for (int l = 3; l < lmax; l++) {
    double ll = (double) l;
    fac1[l] = sqrt((2 * ll + 3) * (2 * ll + 1) / ((ll - 2) * (ll + 4) * ll * (ll + 2))) * (ll + 1);
    fac2[l] = 3.0 / (ll * (ll + 1));
    fac3[l] = sqrt((2 * ll + 3) / (2 * ll - 1) * (ll - 3) * (ll + 3) * (ll - 1) * (ll + 1) /
                   ((ll - 2) * (ll + 4) * ll * (ll + 2))) *
              (ll + 1) / ll;
    fac4[l] = sqrt(2. / (2 * ll + 3));
  }

  for (int index_mu = 0; index_mu < num_mu; index_mu++) {
    d31[index_mu][0] = 0;
    d31[index_mu][1] = 0;
    double dlm1      = 0.; /*l=2*/
    d31[index_mu][2] = 0;
    double dl = sqrt(105. / 2.) * (1 + mu[index_mu]) * (1 + mu[index_mu]) * (1 - mu[index_mu]) /
                8.; /*l=3*/
    d31[index_mu][3] = dl * sqrt(2. / 7.);
    for (int l = 3; l < lmax; l++) {
      /* sqrt((2l+1)/2)*d22 recurrence, supposed to be more stable */
      double dlp1          = fac1[l] * (mu[index_mu] - fac2[l]) * dl - fac3[l] * dlm1;
      d31[index_mu][l + 1] = dlp1 * fac4[l];
      dlm1                 = dl;
      dl                   = dlp1;
    }
  }
  return _SUCCESS_;
}

/**
 * This routine computes the d3m1 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d3m1   Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int LensingModule::lensing_d3m1(double* mu, int num_mu, int lmax, double** d3m1) {
  std::vector<double> fac1(lmax), fac2(lmax), fac3(lmax), fac4(lmax);
  for (int l = 3; l < lmax; l++) {
    double ll = (double) l;
    fac1[l] = sqrt((2 * ll + 3) * (2 * ll + 1) / ((ll - 2) * (ll + 4) * ll * (ll + 2))) * (ll + 1);
    fac2[l] = 3.0 / (ll * (ll + 1));
    fac3[l] = sqrt((2 * ll + 3) / (2 * ll - 1) * (ll - 3) * (ll + 3) * (ll - 1) * (ll + 1) /
                   ((ll - 2) * (ll + 4) * ll * (ll + 2))) *
              (ll + 1) / ll;
    fac4[l] = sqrt(2. / (2 * ll + 3));
  }

  for (int index_mu = 0; index_mu < num_mu; index_mu++) {
    d3m1[index_mu][0] = 0;
    d3m1[index_mu][1] = 0;
    double dlm1       = 0.; /*l=2*/
    d3m1[index_mu][2] = 0;
    double dl = sqrt(105. / 2.) * (1 + mu[index_mu]) * (1 - mu[index_mu]) * (1 - mu[index_mu]) /
                8.; /*l=3*/
    d3m1[index_mu][3] = dl * sqrt(2. / 7.);
    for (int l = 3; l < lmax; l++) {
      /* sqrt((2l+1)/2)*d22 recurrence, supposed to be more stable */
      double dlp1           = fac1[l] * (mu[index_mu] + fac2[l]) * dl - fac3[l] * dlm1;
      d3m1[index_mu][l + 1] = dlp1 * fac4[l];
      dlm1                  = dl;
      dl                    = dlp1;
    }
  }
  return _SUCCESS_;
}

/**
 * This routine computes the d3m3 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d3m3   Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int LensingModule::lensing_d3m3(double* mu, int num_mu, int lmax, double** d3m3) {
  std::vector<double> fac1(lmax), fac2(lmax), fac3(lmax), fac4(lmax);
  for (int l = 3; l < lmax; l++) {
    double ll = (double) l;
    fac1[l]   = sqrt((2 * ll + 3) * (2 * ll + 1)) * (ll + 1) / ((ll - 2) * (ll + 4));
    fac2[l]   = 9.0 / (ll * (ll + 1));
    fac3[l]   = sqrt((2 * ll + 3) / (2 * ll - 1)) * (ll - 3) * (ll + 3) * (l + 1) /
                ((ll - 2) * (ll + 4) * ll);
    fac4[l]   = sqrt(2. / (2 * ll + 3));
  }

  for (int index_mu = 0; index_mu < num_mu; index_mu++) {
    d3m3[index_mu][0] = 0;
    d3m3[index_mu][1] = 0;
    double dlm1       = 0.; /*l=2*/
    d3m3[index_mu][2] = 0;
    double dl = sqrt(7. / 2.) * (1 - mu[index_mu]) * (1 - mu[index_mu]) * (1 - mu[index_mu]) /
                8.; /*l=3*/
    d3m3[index_mu][3] = dl * sqrt(2. / 7.);
    for (int l = 3; l < lmax; l++) {
      /* sqrt((2l+1)/2)*d22 recurrence, supposed to be more stable */
      double dlp1           = fac1[l] * (mu[index_mu] + fac2[l]) * dl - fac3[l] * dlm1;
      d3m3[index_mu][l + 1] = dlp1 * fac4[l];
      dlm1                  = dl;
      dl                    = dlp1;
    }
  }
  return _SUCCESS_;
}

/**
 * This routine computes the d40 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d40    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int LensingModule::lensing_d40(double* mu, int num_mu, int lmax, double** d40) {
  std::vector<double> fac1(lmax), fac3(lmax), fac4(lmax);
  for (int l = 4; l < lmax; l++) {
    double ll = (double) l;
    fac1[l]   = sqrt((2 * ll + 3) * (2 * ll + 1) / ((ll - 3) * (ll + 5)));
    fac3[l]   = sqrt((2 * ll + 3) * (ll - 4) * (ll + 4) / ((2 * ll - 1) * (ll - 3) * (ll + 5)));
    fac4[l]   = sqrt(2. / (2 * ll + 3));
  }

  for (int index_mu = 0; index_mu < num_mu; index_mu++) {
    d40[index_mu][0] = 0;
    d40[index_mu][1] = 0;
    d40[index_mu][2] = 0;
    double dlm1      = 0.; /*l=3*/
    d40[index_mu][3] = 0;
    double dl        = sqrt(315.) * (1 + mu[index_mu]) * (1 + mu[index_mu]) * (1 - mu[index_mu]) *
                       (1 - mu[index_mu]) / 16.; /*l=4*/
    d40[index_mu][4] = dl * sqrt(2. / 9.);
    for (int l = 4; l < lmax; l++) {
      /* sqrt((2l+1)/2)*d22 recurrence, supposed to be more stable */
      double dlp1          = fac1[l] * mu[index_mu] * dl - fac3[l] * dlm1;
      d40[index_mu][l + 1] = dlp1 * fac4[l];
      dlm1                 = dl;
      dl                   = dlp1;
    }
  }
  return _SUCCESS_;
}

/**
 * This routine computes the d4m2 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d4m2   Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int LensingModule::lensing_d4m2(double* mu, int num_mu, int lmax, double** d4m2) {
  std::vector<double> fac1(lmax), fac2(lmax), fac3(lmax), fac4(lmax);
  for (int l = 4; l < lmax; l++) {
    double ll = (double) l;
    fac1[l]   = sqrt((2 * ll + 3) * (2 * ll + 1) / ((ll - 3) * (ll + 5) * (ll - 1) * (ll + 3))) *
                (ll + 1.);
    fac2[l]   = 8. / (ll * (ll + 1));
    fac3[l]   = sqrt((2 * ll + 3) * (ll - 4) * (ll + 4) * (ll - 2) * (ll + 2) /
                     ((2 * ll - 1) * (ll - 3) * (ll + 5) * (ll - 1) * (ll + 3))) *
                (ll + 1) / ll;
    fac4[l]   = sqrt(2. / (2 * ll + 3));
  }

  for (int index_mu = 0; index_mu < num_mu; index_mu++) {
    d4m2[index_mu][0] = 0;
    d4m2[index_mu][1] = 0;
    d4m2[index_mu][2] = 0;
    double dlm1       = 0.; /*l=3*/
    d4m2[index_mu][3] = 0;
    double dl         = sqrt(126.) * (1 + mu[index_mu]) * (1 - mu[index_mu]) * (1 - mu[index_mu]) *
                        (1 - mu[index_mu]) / 16.; /*l=4*/
    d4m2[index_mu][4] = dl * sqrt(2. / 9.);
    for (int l = 4; l < lmax; l++) {
      /* sqrt((2l+1)/2)*d22 recurrence, supposed to be more stable */
      double dlp1           = fac1[l] * (mu[index_mu] + fac2[l]) * dl - fac3[l] * dlm1;
      d4m2[index_mu][l + 1] = dlp1 * fac4[l];
      dlm1                  = dl;
      dl                    = dlp1;
    }
  }
  return _SUCCESS_;
}

/**
 * This routine computes the d4m4 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d4m4   Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

int LensingModule::lensing_d4m4(double* mu, int num_mu, int lmax, double** d4m4) {
  std::vector<double> fac1(lmax), fac2(lmax), fac3(lmax), fac4(lmax);
  for (int l = 4; l < lmax; l++) {
    double ll = (double) l;
    fac1[l]   = sqrt((2 * ll + 3) * (2 * ll + 1)) * (ll + 1) / ((ll - 3) * (ll + 5));
    fac2[l]   = 16. / (ll * (ll + 1));
    fac3[l]   = sqrt((2 * ll + 3) / (2 * ll - 1)) * (ll - 4) * (ll + 4) * (ll + 1) /
                ((ll - 3) * (ll + 5) * ll);
    fac4[l]   = sqrt(2. / (2 * ll + 3));
  }

  for (int index_mu = 0; index_mu < num_mu; index_mu++) {
    d4m4[index_mu][0] = 0;
    d4m4[index_mu][1] = 0;
    d4m4[index_mu][2] = 0;
    double dlm1       = 0.; /*l=3*/
    d4m4[index_mu][3] = 0;
    double dl = sqrt(9. / 2.) * (1 - mu[index_mu]) * (1 - mu[index_mu]) * (1 - mu[index_mu]) *
                (1 - mu[index_mu]) / 16.; /*l=4*/
    d4m4[index_mu][4] = dl * sqrt(2. / 9.);
    for (int l = 4; l < lmax; l++) {
      /* sqrt((2l+1)/2)*d22 recurrence, supposed to be more stable */
      double dlp1           = fac1[l] * (mu[index_mu] + fac2[l]) * dl - fac3[l] * dlm1;
      d4m4[index_mu][l + 1] = dlp1 * fac4[l];
      dlm1                  = dl;
      dl                    = dlp1;
    }
  }
  return _SUCCESS_;
}
