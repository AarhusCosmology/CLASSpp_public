/** @file lensing_module.cpp
 * Computes lensed CMB power spectra from unlensed spectra using the full-sky
 * method of Challinor and Lewis for the lifetime of a LensingModule instance.
 */

#include "lensing_module.h"

#include <algorithm>
#include <array>
#include <numeric>

#include "errors.h"
#include "spectra_module.h"
#include "thread_pool.h"
#include "wigner_d.h"

/**
 * Anisotropy power spectra \f$ C_l\f$'s for all types, modes and initial conditions.
 * SO FAR: ONLY SCALAR
 *
 * This routine evaluates all the lensed \f$ C_l\f$'s at a given value of l by
 * picking it in the pre-computed table. When relevant, it also
 * sums over all initial conditions for each mode, and over all modes.
 *
 * This function can be called from whatever module at whatever time,
 * valid for the lifetime of the module instance.
 *
 */

LensingModule::LensingModule(InputModulePtr input_module, SpectraModulePtr spectra_module)
    : BaseModule(std::move(input_module)), spectra_module_(std::move(spectra_module)) {
  lensing_init();
}

LensingModule::~LensingModule() = default;

std::map<std::string, std::vector<double>> LensingModule::cl_output(int lmax) const {
  class_test_severe(!ple->has_lensed_cls,
                    "No lensed Cls were computed. Add 'lCl' to the list of outputs and set "
                    "'lensing' to 'yes'.");
  class_test_severe((lmax > l_lensed_max_) || (lmax < 0),
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
  class_test_severe(!ple->has_lensed_cls, "No lensed Cls was computed, adjust your inputs.\n");

  std::map<std::string, int> index_map;

  if (has_tt_) {
    index_map["tt"] = index_lt_tt_;
  }
  if (has_ee_) {
    index_map["ee"] = index_lt_ee_;
  }
  if (has_te_) {
    index_map["te"] = index_lt_te_;
  }
  if (has_bb_) {
    index_map["bb"] = index_lt_bb_;
  }
  if (has_pp_) {
    index_map["pp"] = index_lt_pp_;
  }
  if (has_tp_) {
    index_map["tp"] = index_lt_tp_;
  }

  // Create vectors for fast iteration in nested loop below.
  std::vector<std::vector<double>> data_vectors;
  std::vector<int> indices;
  std::vector<std::string> keys;

  for (const auto& element : index_map) {
    data_vectors.emplace_back(l_values.size() + 2, 0.0);
    indices.push_back(element.second);
    keys.push_back(element.first);
  }

  std::vector<double> cl_lensed(lt_size_);
  for (int index_l = 0; index_l < l_values.size(); ++index_l) {
    lensing_cl_at_l(l_values[index_l], cl_lensed.data());
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

void LensingModule::lensing_cl_at_l(int l, double* cl_lensed) const {
  class_test(l > l_lensed_max_,
             "you asked for lensed Cls at l=%d, they were computed only up to l=%d, you should "
             "increase l_max_scalars or decrease the precision parameter delta_l_max",
             l,
             l_lensed_max_);

  int last_index;
  array_interpolate_spline(l_.data(),
                           l_size_,
                           cl_lens_.data(),
                           ddcl_lens_.data(),
                           lt_size_,
                           l,
                           &last_index,
                           cl_lensed,
                           lt_size_);
  for (const auto& [index_lt, spline] : acoustic_) {
    cl_lensed[index_lt] = spline(l);
  }

  /* set to zero for the types such that l<l_max */
  for (int index_lt = 0; index_lt < lt_size_; index_lt++) {
    if (l > l_max_lt_[index_lt]) {
      cl_lensed[index_lt] = 0.;
    }
  }
}

/**
 * This routine initializes the lensing structure (in particular,
 * computes table of lensed anisotropy spectra \f$ C_l^{X} \f$)
 *
 * @return the error status
 */

void LensingModule::lensing_init() {
  /** Summary: */

  /** - check that we really want to compute at least one spectrum */

  if (!ple->has_lensed_cls) {
    if (ple->lensing_verbose > 0) {
      printf("No lensing requested. Lensing module skipped.\n");
    }
    return;
  }
  if (ple->lensing_verbose > 0) {
    printf("Computing lensed spectra ");
    if (ppr->accurate_lensing) {
      printf("(accurate mode)\n");
    }
    else {
      printf("(fast mode)\n");
    }
  }

  /** - initialize indices and allocate some of the arrays in the
      lensing structure */

  lensing_indices();

  /** - put all precision variables hare; will be stored later in precision structure */
  /** - Last element in \f$ \mu \f$ will be for \f$ \mu=1 \f$, needed for sigma2.
      The rest will be chosen as roots of a Gauss-Legendre quadrature **/

  int num_mu;
  if (ppr->accurate_lensing) {
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

  if (ppr->accurate_lensing) {
    quadrature_gauss_legendre(mu.data(), w8.data(), num_mu - 1, ppr->tol_gauss_legendre);
  }
  else { /* Crude integration on [0,pi/16]: Riemann sum on theta */

    double delta_theta = _PI_ / 16. / static_cast<double>(num_mu - 1);
    for (int index_mu = 0; index_mu < num_mu - 1; index_mu++) {
      double theta = (index_mu + 1) * delta_theta;
      mu[index_mu] = cos(theta);
      w8[index_mu] = sin(theta) * delta_theta; /* We integrate on mu */
    }
  }

  /** - Locally store unlensed temperature \f$ cl_{tt}\f$ and potential \f$ cl_{pp}\f$ spectra **/
  std::vector<double> cl_unlensed(spectra_module_->ct_size_);
  std::vector<double> cl_tt(l_unlensed_max_ + 1);
  std::vector<double> cl_te;
  if (has_te_) {
    cl_te.resize(l_unlensed_max_ + 1);
  }
  std::vector<double> cl_ee;
  std::vector<double> cl_bb;
  if (has_ee_ || has_bb_) {
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
    spectra_module_->spectra_cl_at_l(l, cl_unlensed.data(), cl_md.data(), cl_md_ic.data());
    /* lensed EE and BB can be requested without TT (output = pCl,lCl): index_lt_tt_ is then unset */
    if (has_tt_) {
      cl_tt[l] = cl_unlensed[index_lt_tt_];
    }
    cl_pp[l] = cl_unlensed[index_lt_pp_];
    if (has_te_) {
      cl_te[l] = cl_unlensed[index_lt_te_];
    }
    if (has_ee_ || has_bb_) {
      cl_ee[l] = cl_unlensed[index_lt_ee_];
      cl_bb[l] = cl_unlensed[index_lt_bb_];
    }
  }

  std::vector<double> sqrt1(l_unlensed_max_ + 1);
  std::vector<double> sqrt2(l_unlensed_max_ + 1);
  std::vector<double> sqrt3(l_unlensed_max_ + 1);
  std::vector<double> sqrt4(l_unlensed_max_ + 1);
  std::vector<double> sqrt5(l_unlensed_max_ + 1);
  for (int l = 2; l <= l_unlensed_max_; l++) {
    double ll = static_cast<double>(l);
    sqrt1[l]  = sqrt((ll + 2) * (ll + 1) * ll * (ll - 1));
    sqrt2[l]  = sqrt((ll + 2) * (ll - 1));
    sqrt3[l]  = sqrt((ll + 3) * (ll - 2));
    sqrt4[l]  = sqrt((ll + 4) * (ll + 3) * (ll - 2.) * (ll - 3));
    sqrt5[l]  = sqrt(ll * (ll + 1));
  }

  /** - \f$ Cgl(\mu=1) \f$, the point sigma2\f$(\mu)\f$ = Cgl(1) - Cgl(\f$\mu\f$) is measured from */
  double Cgl_at_1 = 0.;
  {
    wigner::Nodes mu_at_1{};
    mu_at_1.fill(mu[num_mu - 1]);
    std::vector<double> d11_at_1(static_cast<size_t>(l_unlensed_max_ + 1) * wigner::kChunk);
    wigner::D11(mu_at_1, l_unlensed_max_, d11_at_1.data());
    for (int l = 2; l <= l_unlensed_max_; l++) {
      Cgl_at_1 += (2. * l + 1.) * l * (l + 1.) * cl_pp[l] *
                  d11_at_1[static_cast<size_t>(l) * wigner::kChunk];
    }
    Cgl_at_1 /= 4. * _PI_;
  }

  /** - Lensed correlation functions and their transform back to multipoles, streamed over \f$ \mu \f$.

      Both steps need every Wigner \f$ d^l_{mm'}(\mu) \f$ at every node and multipole, for up to
      twelve (m,m') pairs. Tabulated all at once that is 12 x num_mu x l_max doubles -- 12 GB at
      the l_max ~ 11000 a high-l CMB likelihood needs, with accurate_lensing's num_mu ~ l_max.
      Neither step needs more than a few nodes at a time, so the nodes are cut into chunks of
      wigner::kChunk: a task tabulates the d-functions of its chunk only, as d[l][node], forms the
      correlation functions there, and adds the chunk's share of the Gauss-Legendre (or Riemann)
      sum into that chunk's own row of partial sums. The partial sums are then added in chunk
      order. Chunks, not threads, are the unit of summation, so the result does not depend on the
      number of threads. Every loop runs over the nodes of the chunk innermost, so that the
      recurrences and the sums over l are vectorised across nodes; each node's sums still add
      their terms in increasing l. */

  const int num_nodes  = num_mu - 1; /* the last mu (=1) only serves Cgl_at_1 */
  const int num_chunks = (num_nodes + wigner::kChunk - 1) / wigner::kChunk;
  const int row_size   = l_unlensed_max_ + 1;

  std::vector<double> partial_tt;
  std::vector<double> partial_te;
  std::vector<double> partial_p;
  std::vector<double> partial_m;
  if (has_tt_) {
    partial_tt.assign(static_cast<size_t>(num_chunks) * l_size_, 0.0);
  }
  if (has_te_) {
    partial_te.assign(static_cast<size_t>(num_chunks) * l_size_, 0.0);
  }
  if (has_ee_ || has_bb_) {
    partial_p.assign(static_cast<size_t>(num_chunks) * l_size_, 0.0);
    partial_m.assign(static_cast<size_t>(num_chunks) * l_size_, 0.0);
  }

  const int num_d = 4 + (has_te_ ? 3 : 0) + ((has_ee_ || has_bb_) ? 5 : 0);

  /* One task per thread, each walking its own stride of chunks with one reused buffer. */
  auto process_chunks = [&](int first_chunk, int stride) {
    std::vector<double> buf(static_cast<size_t>(num_d) * row_size * wigner::kChunk);
    int next_d  = 0;
    auto table  = [&]() { return &buf[static_cast<size_t>(next_d++) * row_size * wigner::kChunk]; };
    double* d00 = table();
    double* d11 = table();
    double* d1m1 = table();
    double* d2m2 = table();
    double* d20  = nullptr;
    double* d3m1 = nullptr;
    double* d4m2 = nullptr;
    double* d22  = nullptr;
    double* d31  = nullptr;
    double* d3m3 = nullptr;
    double* d40  = nullptr;
    double* d4m4 = nullptr;
    if (has_te_) {
      d20  = table();
      d3m1 = table();
      d4m2 = table();
    }
    if (has_ee_ || has_bb_) {
      d22  = table();
      d31  = table();
      d3m3 = table();
      d40  = table();
      d4m4 = table();
    }
    /* element (l, node) of a table */
    auto at = [](const double* d, int l, int j) {
      return d[static_cast<size_t>(l) * wigner::kChunk + j];
    };

    for (int index_chunk = first_chunk; index_chunk < num_chunks; index_chunk += stride) {
      const int mu_first = index_chunk * wigner::kChunk;
      const int n        = std::min(wigner::kChunk, num_nodes - mu_first);

      /* the last chunk is padded with copies of its last node, whose results are not used */
      wigner::Nodes mu_chunk{};
      for (int j = 0; j < wigner::kChunk; j++) {
        mu_chunk[j] = mu[mu_first + std::min(j, n - 1)];
      }

      /** - --> Wigner d-functions of this chunk, computed by recurrence in l */
      wigner::D00(mu_chunk, l_unlensed_max_, d00);
      wigner::D11(mu_chunk, l_unlensed_max_, d11);
      wigner::D1m1(mu_chunk, l_unlensed_max_, d1m1);
      wigner::D2m2(mu_chunk, l_unlensed_max_, d2m2);
      if (has_te_) {
        wigner::D20(mu_chunk, l_unlensed_max_, d20);
        wigner::D3m1(mu_chunk, l_unlensed_max_, d3m1);
        wigner::D4m2(mu_chunk, l_unlensed_max_, d4m2);
      }
      if (has_ee_ || has_bb_) {
        wigner::D22(mu_chunk, l_unlensed_max_, d22);
        wigner::D31(mu_chunk, l_unlensed_max_, d31);
        wigner::D3m3(mu_chunk, l_unlensed_max_, d3m3);
        wigner::D40(mu_chunk, l_unlensed_max_, d40);
        wigner::D4m4(mu_chunk, l_unlensed_max_, d4m4);
      }

      /** - --> Cgl(\f$\mu\f$), Cgl2(\f$\mu\f$) and sigma2(\f$\mu\f$) at the nodes of the chunk */
      wigner::Nodes Cgl{};
      wigner::Nodes Cgl2{};
      for (int l = 2; l <= l_unlensed_max_; l++) {
        for (int j = 0; j < n; j++) {
          Cgl[j]  += (2. * l + 1.) * l * (l + 1.) * cl_pp[l] * at(d11, l, j);
          Cgl2[j] += (2. * l + 1.) * l * (l + 1.) * cl_pp[l] * at(d1m1, l, j);
        }
      }
      wigner::Nodes sigma2{};
      for (int j = 0; j < n; j++) {
        Cgl[j]    /= 4. * _PI_;
        Cgl2[j]   /= 4. * _PI_;
        sigma2[j]  = Cgl_at_1 - Cgl[j];
      }

      /** - --> ksi (TT), ksiX (TE), ksip and ksim (EE, BB) at the nodes of the chunk */
      wigner::Nodes ksi{};
      wigner::Nodes ksiX{};
      wigner::Nodes ksip{};
      wigner::Nodes ksim{};

      for (int l = 2; l <= l_unlensed_max_; l++) {
        double ll = static_cast<double>(l);

        double fac  = ll * (ll + 1) / 4.;
        double fac1 = (2 * ll + 1) / (4. * _PI_);

        for (int j = 0; j < n; j++) {
          /* In the following we will keep terms of the form (sigma2)^k*(Cgl2)^m
             with k+m <= 2 */

          double X_000  = exp(-fac * sigma2[j]);
          double X_p000 = -fac * X_000;
          /* X_220 = 0.25*sqrt1[l] * exp(-(fac-0.5)*sigma2); */
          double X_220 = 0.25 * sqrt1[l] * X_000; /* Order 0 */
          /* next 5 lines useless, but avoid compiler warning 'may be used uninitialized' */
          double X_242  = 0.;
          double X_132  = 0.;
          double X_121  = 0.;
          double X_p022 = 0.;
          double X_022  = 0.;

          if (has_te_ || has_ee_ || has_bb_) {
            /* X_022 = exp(-(fac-1.)*sigma2); */
            X_022  = X_000 * (1 + sigma2[j] * (1 + 0.5 * sigma2[j])); /* Order 2 */
            X_p022 = -(fac - 1.) * X_022; /* Old versions were missing the
            minus sign in this line, which introduced a very small error
            on the high-l C_l^TE lensed spectrum [credits for bug fix:
            Selim Hotinli] */

            /* X_242 = 0.25*sqrt4[l] * exp(-(fac-5./2.)*sigma2); */
            X_242 = 0.25 * sqrt4[l] * X_000; /* Order 0 */
            if (has_ee_ || has_bb_) {
              /* X_121 = - 0.5*sqrt2[l] * exp(-(fac-2./3.)*sigma2);
                 X_132 = - 0.5*sqrt3[l] * exp(-(fac-5./3.)*sigma2); */
              X_121 = -0.5 * sqrt2[l] * X_000 * (1 + 2. / 3. * sigma2[j]); /* Order 1 */
              X_132 = -0.5 * sqrt3[l] * X_000 * (1 + 5. / 3. * sigma2[j]); /* Order 1 */
            }
          }

          if (has_tt_) {
            double res = fac1 * cl_tt[l];

            double lens = (X_000 * X_000 * at(d00, l, j) +
                           X_p000 * X_p000 * at(d1m1, l, j) * Cgl2[j] * 8. / (ll * (ll + 1)) +
                           (X_p000 * X_p000 * at(d00, l, j) + X_220 * X_220 * at(d2m2, l, j)) *
                               Cgl2[j] * Cgl2[j]);
            if (!ppr->accurate_lensing) {
              /* Remove unlensed correlation function */
              lens -= at(d00, l, j);
            }
            res    *= lens;
            ksi[j] += res;
          }

          if (has_te_) {
            double resX = fac1 * cl_te[l];

            double lens = (X_022 * X_000 * at(d20, l, j) +
                           Cgl2[j] * 2. * X_p000 / sqrt5[l] *
                               (X_121 * at(d11, l, j) + X_132 * at(d3m1, l, j)) +
                           0.5 * Cgl2[j] * Cgl2[j] *
                               ((2. * X_p022 * X_p000 + X_220 * X_220) * at(d20, l, j) +
                                X_220 * X_242 * at(d4m2, l, j)));
            if (!ppr->accurate_lensing) {
              lens -= at(d20, l, j);
            }
            resX    *= lens;
            ksiX[j] += resX;
          }

          if (has_ee_ || has_bb_) {
            double resp = fac1 * (cl_ee[l] + cl_bb[l]);
            double resm = fac1 * (cl_ee[l] - cl_bb[l]);

            double lensp = (X_022 * X_022 * at(d22, l, j) +
                            2. * Cgl2[j] * X_132 * X_121 * at(d31, l, j) +
                            Cgl2[j] * Cgl2[j] *
                                (X_p022 * X_p022 * at(d22, l, j) + X_242 * X_220 * at(d40, l, j)));

            double lensm = (X_022 * X_022 * at(d2m2, l, j) +
                            Cgl2[j] *
                                (X_121 * X_121 * at(d1m1, l, j) + X_132 * X_132 * at(d3m3, l, j)) +
                            0.5 * Cgl2[j] * Cgl2[j] *
                                (2. * X_p022 * X_p022 * at(d2m2, l, j) +
                                 X_220 * X_220 * at(d00, l, j) + X_242 * X_242 * at(d4m4, l, j)));
            if (!ppr->accurate_lensing) {
              lensp -= at(d22, l, j);
              lensm -= at(d2m2, l, j);
            }
            resp    *= lensp;
            resm    *= lensm;
            ksip[j] += resp;
            ksim[j] += resm;
          }
        }
      }

      /** - --> this chunk's share of the transform back to the lensed \f$ C_l\f$'s */
      const size_t offset = static_cast<size_t>(index_chunk) * l_size_;
      for (int index_l = 0; index_l < l_size_; index_l++) {
        const int l = static_cast<int>(l_[index_l]);
        for (int j = 0; j < n; j++) {
          const int index_mu = mu_first + j;
          if (has_tt_) {
            partial_tt[offset + index_l] += ksi[j] * at(d00, l, j) * w8[index_mu];
          }
          if (has_te_) {
            partial_te[offset + index_l] += ksiX[j] * at(d20, l, j) * w8[index_mu];
          }
          if (has_ee_ || has_bb_) {
            partial_p[offset + index_l] += ksip[j] * at(d22, l, j) * w8[index_mu];
            partial_m[offset + index_l] += ksim[j] * at(d2m2, l, j) * w8[index_mu];
          }
        }
      }
    }
  };

  const int num_tasks = std::max(1, std::min(static_cast<int>(pba->number_of_threads), num_chunks));
  Tools::TaskSystem task_system(num_tasks);
  std::vector<std::future<void>> tasks;
  for (int index_task = 0; index_task < num_tasks; index_task++) {
    tasks.push_back(
        task_system.AsyncTask([&, index_task]() { process_chunks(index_task, num_tasks); }));
  }
  for (auto& task : tasks) {
    task.get();
  }

  /** - compute lensed \f$ C_l\f$'s by adding up the chunks' partial sums, in chunk order */
  for (int index_l = 0; index_l < l_size_; index_l++) {
    double cltt = 0.;
    double clte = 0.;
    double clp  = 0.;
    double clm  = 0.;
    for (int index_chunk = 0; index_chunk < num_chunks; index_chunk++) {
      const size_t index = static_cast<size_t>(index_chunk) * l_size_ + index_l;
      if (has_tt_) {
        cltt += partial_tt[index];
      }
      if (has_te_) {
        clte += partial_te[index];
      }
      if (has_ee_ || has_bb_) {
        clp += partial_p[index];
        clm += partial_m[index];
      }
    }
    if (has_tt_) {
      cl_lens_[index_l * lt_size_ + index_lt_tt_] = cltt * 2.0 * _PI_;
    }
    if (has_te_) {
      cl_lens_[index_l * lt_size_ + index_lt_te_] = clte * 2.0 * _PI_;
    }
    if (has_ee_ || has_bb_) {
      cl_lens_[index_l * lt_size_ + index_lt_ee_] = (clp + clm) * _PI_;
      cl_lens_[index_l * lt_size_ + index_lt_bb_] = (clp - clm) * _PI_;
    }
  }

  if (!ppr->accurate_lensing) {
    if (has_tt_) {
      lensing_addback_cl_tt(cl_tt.data());
    }
    if (has_te_) {
      lensing_addback_cl_te(cl_te.data());
    }
    if (has_ee_ || has_bb_) {
      lensing_addback_cl_ee_bb(cl_ee.data(), cl_bb.data());
    }
  }

  /** - spline computed \f$ C_l\f$'s in view of interpolation */

  array_spline_table_lines(l_.data(),
                           l_size_,
                           cl_lens_.data(),
                           lt_size_,
                           ddcl_lens_.data(),
                           _SPLINE_EST_DERIV_);

  /** - the lensed TT, EE and TE still oscillate with the acoustic peaks: interpolated like the
      unlensed ones (SpectraModule::spectra_cls) */
  for (const int index_lt : {index_lt_tt_, index_lt_ee_, index_lt_te_}) {
    if ((index_lt >= 0) && AcousticSpline::Applies(l_.data(),
                                                   l_size_,
                                                   cl_lens_.data() + index_lt,
                                                   lt_size_,
                                                   spectra_module_->acoustic_omega_)) {
      acoustic_.emplace_back(index_lt,
                             AcousticSpline(l_.data(),
                                            l_size_,
                                            cl_lens_.data() + index_lt,
                                            lt_size_,
                                            spectra_module_->acoustic_omega_));
    }
  }

  /** - Exit **/
}

/**
 * This routine defines indices and allocates tables in the lensing structure
 *
 * @return the error status
 */

void LensingModule::lensing_indices() {
  /* indices of all Cl types (lensed and unlensed) */

  if (spectra_module_->has_tt_) {
    has_tt_      = true;
    index_lt_tt_ = spectra_module_->index_ct_tt_;
  }
  else {
    has_tt_ = false;
  }

  if (spectra_module_->has_ee_) {
    has_ee_      = true;
    index_lt_ee_ = spectra_module_->index_ct_ee_;
  }
  else {
    has_ee_ = false;
  }

  if (spectra_module_->has_te_) {
    has_te_      = true;
    index_lt_te_ = spectra_module_->index_ct_te_;
  }
  else {
    has_te_ = false;
  }

  if (spectra_module_->has_bb_) {
    has_bb_      = true;
    index_lt_bb_ = spectra_module_->index_ct_bb_;
  }
  else {
    has_bb_ = false;
  }

  if (spectra_module_->has_pp_) {
    has_pp_      = true;
    index_lt_pp_ = spectra_module_->index_ct_pp_;
  }
  else {
    has_pp_ = false;
  }

  if (spectra_module_->has_tp_) {
    has_tp_      = true;
    index_lt_tp_ = spectra_module_->index_ct_tp_;
  }
  else {
    has_tp_ = false;
  }

  if (spectra_module_->has_dd_) {
    has_dd_      = true;
    index_lt_dd_ = spectra_module_->index_ct_dd_;
  }
  else {
    has_dd_ = false;
  }

  if (spectra_module_->has_td_) {
    has_td_      = true;
    index_lt_td_ = spectra_module_->index_ct_td_;
  }
  else {
    has_td_ = false;
  }

  if (spectra_module_->has_ll_) {
    has_ll_      = true;
    index_lt_ll_ = spectra_module_->index_ct_ll_;
  }
  else {
    has_ll_ = false;
  }

  if (spectra_module_->has_tl_) {
    has_tl_      = true;
    index_lt_tl_ = spectra_module_->index_ct_tl_;
  }
  else {
    has_tl_ = false;
  }

  lt_size_ = spectra_module_->ct_size_;

  /* number of multipoles */

  l_unlensed_max_ = spectra_module_->l_max_tot_;

  l_lensed_max_ = l_unlensed_max_ - ppr->delta_l_max;

  int index_l;
  for (index_l = 0;
       (index_l < spectra_module_->l_size_max_) && (spectra_module_->l_[index_l] <= l_lensed_max_);
       index_l++) {
    ;
  }

  if (index_l < spectra_module_->l_size_max_) {
    index_l++; /* one more point in order to be able to interpolate till l_lensed_max_ */
  }

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
    spectra_module_->spectra_cl_at_l(l_[index_l],
                                     &(cl_lens_[index_l * lt_size_]),
                                     cl_md_ptrs.data(),
                                     cl_md_ic_ptrs.data());
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
      l_max_lt_[index_lt] = std::max(l_max_lt_[index_lt],
                                     spectra_module_->l_max_ct_[index_md][index_lt]);

      if ((has_bb_) && (has_ee_) && (index_lt == index_lt_bb_)) {
        l_max_lt_[index_lt] = std::max(l_max_lt_[index_lt],
                                       spectra_module_->l_max_ct_[index_md][index_lt_ee_]);
      }
    }
  }
}

/**
 * This routine adds back the unlensed \f$ cl_{tt}\f$ power spectrum
 * Used in case of fast (and BB inaccurate) integration of
 * correlation functions.
 *
 * @return the error status
 */

void LensingModule::lensing_addback_cl_tt(const double* cl_tt) {
  for (int index_l = 0; index_l < l_size_; index_l++) {
    int l                                        = static_cast<int>(l_[index_l]);
    cl_lens_[index_l * lt_size_ + index_lt_tt_] += cl_tt[l];
  }
}

/**
 * This routine adds back the unlensed \f$ cl_{te}\f$ power spectrum
 * Used in case of fast (and BB inaccurate) integration of
 * correlation functions.
 *
 * @param cl_te Input: Array of unlensed power spectrum
 * @return the error status
 */

void LensingModule::lensing_addback_cl_te(const double* cl_te) {
  for (int index_l = 0; index_l < l_size_; index_l++) {
    int l                                        = static_cast<int>(l_[index_l]);
    cl_lens_[index_l * lt_size_ + index_lt_te_] += cl_te[l];
  }
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

void LensingModule::lensing_addback_cl_ee_bb(const double* cl_ee, const double* cl_bb) {
  for (int index_l = 0; index_l < l_size_; index_l++) {
    int l                                        = static_cast<int>(l_[index_l]);
    cl_lens_[index_l * lt_size_ + index_lt_ee_] += cl_ee[l];
    cl_lens_[index_l * lt_size_ + index_lt_bb_] += cl_bb[l];
  }
}
