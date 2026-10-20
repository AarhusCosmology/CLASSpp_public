/** @file transfer_module.cpp
 * Computes harmonic-space transfer functions for the lifetime of a
 * TransferModule instance.
 */

#include "transfer_module.h"

#include <algorithm>

#include "background_module.h"
#include "nonlinear_module.h"
#include "perturbations_module.h"
#include "thermodynamics_module.h"
#include "thread_pool.h"

TransferModule::TransferModule(InputModulePtr input_module,
                               BackgroundModulePtr background_module,
                               ThermodynamicsModulePtr thermodynamics_module,
                               PerturbationsModulePtr perturbations_module,
                               NonlinearModulePtr nonlinear_module)
    : BaseModule(std::move(input_module)), background_module_(std::move(background_module)),
      thermodynamics_module_(std::move(thermodynamics_module)),
      perturbations_module_(std::move(perturbations_module)),
      nonlinear_module_(std::move(nonlinear_module)) {
  transfer_init();
}

TransferModule::~TransferModule() {}

/**
 * Transfer function \f$ \Delta_l^{X} (q) \f$ at a given wavenumber q.
 *
 * For a given mode (scalar, vector, tensor), initial condition, type
 * (temperature, polarization, lensing, etc) and multipole, computes
 * the transfer function for an arbitrary value of q by interpolating
 * between pre-computed values of q. This
 * function can be called from whatever module at whatever time,
 * valid for the lifetime of the module instance.
 *
 * Wavenumbers are called q in this module and k in the perturbation
 * module. In flat universes k=q. In non-flat universes q and k differ
 * through \f$ q2 = k2 + K(1+m)\f$, where m=0,1,2 for scalar, vector,
 * tensor. q should be used throughout the transfer module, excepted
 * when interpolating or manipulating the source functions S(k,tau)
 * calculated in the perturbation module: for a given value of q, this
 * should be done at the corresponding k(q).
 *
 * @param index_md   Input: index of requested mode
 * @param index_ic   Input: index of requested initial condition
 * @param index_tt   Input: index of requested type
 * @param index_l    Input: index of requested multipole
 * @param q          Input: any wavenumber
 * @param transfer_function Output: transfer function
 * @return the error status
 */

void TransferModule::transfer_functions_at_q(
    int index_md, int index_ic, int index_tt, int index_l, double q, double* transfer_function) {
  /** Summary: */

  /** - interpolate in pre-computed table using array_interpolate_two() */
  array_interpolate_two(q_.data(),
                        1,
                        0,
                        transfer_[index_md].data() +
                            ((index_ic * tt_size_[index_md] + index_tt) * l_size_[index_md] +
                             index_l) *
                                q_size_,
                        q_size_,
                        q,
                        transfer_function,
                        1);
}

/**
 * This routine initializes the transfers structure, (in particular,
 * computes table of transfer functions \f$ \Delta_l^{X} (q) \f$)
 *
 * Main steps:
 *
 * - initialize all indices in the transfers structure
 *   and allocate all its arrays using transfer_indices_of_transfers().
 *
 * - for each thread (in case of parallel run), initialize the fields of a memory zone called the transfer_workspace with transfer_workspace_init()
 *
 * - loop over q values. For each q, compute the Bessel functions if needed with transfer_update_HIS(), and defer the calculation of all transfer functions to transfer_compute_for_each_q()
 * - for each thread, the workspace is destroyed automatically at scope exit
 *
 * @return the error status
 */

void TransferModule::transfer_init() {
  /** Summary: */

  /** - define local variables */

  /* conformal time today */
  double tau0;
  /* conformal time at recombination */
  double tau_rec;
  /* order of magnitude of the oscillation period of transfer functions */
  double q_period;

  /** - check whether any spectrum in harmonic space (i.e., any \f$C_l\f$'s) is actually requested */

  if (!ppt->has_cls) {
    has_cls_             = false;
    do_lcmb_full_limber_ = false;
    if (ptr->transfer_verbose > 0) {
      printf("No harmonic space transfer functions to compute. Transfer module skipped.\n");
    }
    return;
  }
  else {
    has_cls_ = true;
  }

  if (ptr->transfer_verbose > 0) {
    fprintf(stdout, "Computing transfers\n");
  }

  do_lcmb_full_limber_ = ppt->has_cl_cmb_lensing_potential && ppt->want_lcmb_full_limber;

  /** - get number of modes (scalars, tensors...) */

  md_size_ = perturbations_module_->md_size_;

  /** - get conformal age / recombination time
      from background / thermodynamics structures
      (only place where these structures are used in this module) */

  tau0    = background_module_->conformal_age_;
  tau_rec = thermodynamics_module_->tau_rec_;

  /** - order of magnitude of the oscillation period of transfer functions */

  q_period = 2. * _PI_ / (tau0 - tau_rec) * thermodynamics_module_->angular_rescaling_;

  /** - initialize all indices in the transfers structure and
      allocate all its arrays using transfer_indices_of_transfers() */

  transfer_indices_of_transfers(q_period, pba->K, pba->sgnK);

  /** - warn if the closed-space Limber delta is incomplete: past the equator
      the far turning point at pi R - chi_L falls inside the sources and is not
      accounted for (issue #423). */

  /* The CMB lensing source is tabulated from recombination to today, so the
     largest distance it reaches is tau0 - tau_rec, not tau0. */
  if (ppt->has_cl_cmb_lensing_potential &&
      LimberSourceCrossesEquator(tau0 - tau_rec, pba->sgnK, pba->K) && ptr->transfer_verbose > 0) {
    fprintf(stdout,
            " -> [WARNING:] the CMB lensing sources reach past the equator of this closed "
            "universe (tau0-tau_rec=%.4g exceeds (pi/2)/sqrt(K)=%.4g).\n    The Limber "
            "approximation used for the CMB lensing potential above l_switch_limber=%g places "
            "its delta at the near turning point only,\n    and the second turning point at "
            "pi/sqrt(K) - chi now falls inside the sources, so C_l^phiphi is systematically "
            "incomplete there.\n    This needs |Omega_k| greater than about 0.25; reduce the "
            "curvature or raise l_switch_limber to integrate exactly instead.\n",
            tau0 - tau_rec,
            0.5 * _PI_ / sqrt(pba->K),
            ppr->l_switch_limber);
  }

  /** - copy sources to a local array sources (in fact, only the pointers are copied, not the data), and eventually apply non-linear corrections to the sources */

  std::vector<std::vector<double*>> sources_storage(md_size_);
  std::vector<std::vector<std::vector<double>>> owned_sources(md_size_);
  transfer_perturbation_copy_sources_and_nl_corrections(sources_storage, owned_sources);

  /** - spline all the sources passed by the perturbation module with respect to k (in order to interpolate later at a given value of k) */

  std::vector<std::vector<double*>> sources_spline_storage(md_size_);
  std::vector<std::vector<std::vector<double>>> spline_owned_storage(md_size_);
  transfer_perturbation_source_spline(sources_storage,
                                      sources_spline_storage,
                                      spline_owned_storage);

  /** - eventually read the selection and evolution functions */

  transfer_global_selection_read();

  /** - define once for all the time sampling of each type, and its
      window function W(tau) on this sampling */

  transfer_selection_times();
  transfer_time_samplings(tau_rec);
  transfer_precompute_windows(tau_rec);

  /** - compute flat spherical bessel functions */

  double xmax = q_[q_size_ - 1] * tau0;
  if (pba->sgnK == -1 && index_q_flat_approximation_ < q_size_) {
    xmax *= sqrt(pba->sgnK * pba->K) / q_[q_size_ - 1] * (l_[l_size_max_ - 1] + 1) /
            asinh((l_[l_size_max_ - 1] + 1) / q_[q_size_ - 1] * sqrt(pba->sgnK * pba->K)) * 1.01;
  }

  HyperInterpStruct BIS(0,
                        1.,
                        l_size_max_,
                        l_.data(),
                        ppr->hyper_x_min,
                        xmax,
                        ppr->hyper_sampling_flat,
                        l_[l_size_max_ - 1] + 1,
                        ppr->hyper_phi_min_abs);

  std::vector<double**> sources(md_size_);
  std::vector<double**> sources_spline(md_size_);
  for (int index_md = 0; index_md < md_size_; ++index_md) {
    sources[index_md]        = sources_storage[index_md].data();
    sources_spline[index_md] = sources_spline_storage[index_md].data();
  }

  /** - loop over all wavenumbers (parallelized), in the normal scheme and,
      if needed, in the full Limber scheme */

  Tools::TaskSystem task_system(pba->number_of_threads);
  std::vector<std::future<void>> future_output;
  for (int pass = 0; pass < (do_lcmb_full_limber_ ? 2 : 1); pass++) {
    const bool use_full_limber = (pass == 1);
    const int q_size           = use_full_limber ? q_size_limber_ : q_size_;
    for (int index_q = 0; index_q < q_size; index_q++) {
      future_output.push_back(task_system.AsyncTask([this,
                                                     use_full_limber,
                                                     q_size,
                                                     tau_rec,
                                                     tau0,
                                                     index_q,
                                                     &sources,
                                                     &sources_spline,
                                                     &BIS]() {
        struct transfer_workspace ptw = {};
        transfer_workspace_init(&ptw,
                                pba->K,
                                pba->sgnK,
                                tau0 - thermodynamics_module_->tau_cut_,
                                &BIS);

        if (ptr->transfer_verbose > 2) {
          printf("Compute %stransfer for wavenumber [%d/%d]\n",
                 use_full_limber ? "full-Limber " : "",
                 index_q,
                 q_size - 1);
        }

        /* Update interpolation structure (not needed in the full Limber scheme): */
        if (!use_full_limber) {
          transfer_update_HIS(&ptw, index_q, tau0);
        }

        transfer_compute_for_each_q(index_q,
                                    tau_rec,
                                    sources.data(),
                                    sources_spline.data(),
                                    &ptw,
                                    use_full_limber);
      }));
    }
  }

  for (std::future<void>& future : future_output) {
    future.get();
  }
  future_output.clear();
  /** - finally, temporary source tables are cleaned up by RAII */
}

/**
 * This routine defines all indices and allocates all tables
 * in the transfers structure
 *
 * Compute list of (k, l) values, allocate and fill corresponding
 * arrays in the transfers structure. Allocate the array of transfer
 * function tables.
 *
 * @param q_period Input: order of magnitude of the oscillation period of transfer functions
 * @param K        Input: spatial curvature (in absolute value)
 * @param sgnK     Input: spatial curvature sign (open/closed/flat)
 * @return the error status
 */

void TransferModule::transfer_indices_of_transfers(double q_period, double K, int sgnK) {
  /** Summary: */

  /** - define local variables */

  int index_md;
  int index_tt;
  int index_tt_common;

  /** - define indices for transfer types */

  tt_size_.resize(md_size_);

  /** - type indices common to scalars and tensors */

  index_tt = 0;

  class_define_index(index_tt_t2_, ppt->has_cl_cmb_temperature, index_tt, 1);
  class_define_index(index_tt_e_, ppt->has_cl_cmb_polarization, index_tt, 1);

  index_tt_common = index_tt;

  /** - type indices for scalars */

  if (ppt->has_scalars) {
    index_tt = index_tt_common;

    class_define_index(index_tt_t0_, ppt->has_cl_cmb_temperature, index_tt, 1);
    class_define_index(index_tt_t1_, ppt->has_cl_cmb_temperature, index_tt, 1);
    class_define_index(index_tt_lcmb_, ppt->has_cl_cmb_lensing_potential, index_tt, 1);
    class_define_index(index_tt_density_, ppt->has_nc_density, index_tt, ppt->selection_num);
    class_define_index(index_tt_rsd_, ppt->has_nc_rsd, index_tt, ppt->selection_num);
    class_define_index(index_tt_d0_, ppt->has_nc_rsd, index_tt, ppt->selection_num);
    class_define_index(index_tt_d1_, ppt->has_nc_rsd, index_tt, ppt->selection_num);
    class_define_index(index_tt_nc_lens_, ppt->has_nc_lens, index_tt, ppt->selection_num);
    class_define_index(index_tt_nc_g1_, ppt->has_nc_gr, index_tt, ppt->selection_num);
    class_define_index(index_tt_nc_g2_, ppt->has_nc_gr, index_tt, ppt->selection_num);
    class_define_index(index_tt_nc_g3_, ppt->has_nc_gr, index_tt, ppt->selection_num);
    class_define_index(index_tt_nc_g4_, ppt->has_nc_gr, index_tt, ppt->selection_num);
    class_define_index(index_tt_nc_g5_, ppt->has_nc_gr, index_tt, ppt->selection_num);
    class_define_index(index_tt_lensing_,
                       ppt->has_cl_lensing_potential,
                       index_tt,
                       ppt->selection_num);

    tt_size_[perturbations_module_->index_md_scalars_] = index_tt;
  }

  /** - type indices for vectors */

  if (ppt->has_vectors) {
    index_tt = index_tt_common;

    class_define_index(index_tt_t1_, ppt->has_cl_cmb_temperature, index_tt, 1);
    class_define_index(index_tt_b_, ppt->has_cl_cmb_polarization, index_tt, 1);

    tt_size_[perturbations_module_->index_md_vectors_] = index_tt;
  }

  /** - type indices for tensors */

  if (ppt->has_tensors) {
    index_tt = index_tt_common;

    class_define_index(index_tt_b_, ppt->has_cl_cmb_polarization, index_tt, 1);

    tt_size_[perturbations_module_->index_md_tensors_] = index_tt;
  }

  /** - allocate arrays of (k, l) values and transfer functions */

  /* number of l values for each mode and type,
     l_size_tt[index_md][index_tt], and maximized for each mode,
     l_size[index_md] */

  l_size_.resize(md_size_);

  l_size_tt_.resize(md_size_);

  for (index_md = 0; index_md < md_size_; index_md++) {
    l_size_tt_[index_md].resize(tt_size_[index_md]);
  }

  /* array (of array) of transfer functions for each mode, transfer[index_md] */

  transfer_.resize(md_size_);

  /** - get q values using transfer_get_q_list() */

  transfer_get_q_list(q_period, K, sgnK);

  if (do_lcmb_full_limber_) {
    transfer_get_q_limber_list(K, sgnK);
  }

  /** - get k values using transfer_get_k_list() */
  transfer_get_k_list(K);

  /** - define once for all the properties of each transfer type */
  transfer_type_properties();

  /** - get l values using transfer_get_l_list() */
  transfer_get_l_list();

  /** - loop over modes (scalar, etc). For each mode: */

  for (index_md = 0; index_md < md_size_; index_md++) {
    /** - allocate arrays of transfer functions, (transfer_[index_md])[index_ic][index_tt][index_l][index_k] */
    transfer_[index_md].resize(perturbations_module_->ic_size_[index_md] * tt_size_[index_md] *
                               l_size_[index_md] * q_size_);
  }

  if (do_lcmb_full_limber_) {
    transfer_limber_.resize(md_size_);
    const int index_md = perturbations_module_->index_md_scalars_;
    transfer_limber_[index_md].resize(perturbations_module_->ic_size_[index_md] *
                                      tt_size_[index_md] * l_size_[index_md] * q_size_limber_);
  }
}

void TransferModule::transfer_perturbation_copy_sources_and_nl_corrections(
    std::vector<std::vector<double*>>& sources,
    std::vector<std::vector<std::vector<double>>>& owned_sources) {
  for (int index_md = 0; index_md < md_size_; index_md++) {
    sources[index_md].resize(perturbations_module_->ic_size_[index_md] *
                             perturbations_module_->tp_size_[index_md]);
    owned_sources[index_md].resize(perturbations_module_->ic_size_[index_md] *
                                   perturbations_module_->tp_size_[index_md]);

    for (int index_ic = 0; index_ic < perturbations_module_->ic_size_[index_md]; index_ic++) {
      for (int index_tp = 0; index_tp < perturbations_module_->tp_size_[index_md]; index_tp++) {
        if ((pnl->method != nl_none) && (_scalarsEXT_) &&
            (((perturbations_module_->has_source_delta_m_) &&
              (index_tp == perturbations_module_->index_tp_delta_m_)) ||
             ((perturbations_module_->has_source_delta_cb_) &&
              (index_tp == perturbations_module_->index_tp_delta_cb_)) ||
             ((perturbations_module_->has_source_theta_m_) &&
              (index_tp == perturbations_module_->index_tp_theta_m_)) ||
             ((perturbations_module_->has_source_theta_cb_) &&
              (index_tp == perturbations_module_->index_tp_theta_cb_)) ||
             ((perturbations_module_->has_source_phi_) &&
              (index_tp == perturbations_module_->index_tp_phi_)) ||
             ((perturbations_module_->has_source_phi_prime_) &&
              (index_tp == perturbations_module_->index_tp_phi_prime_)) ||
             ((perturbations_module_->has_source_phi_plus_psi_) &&
              (index_tp == perturbations_module_->index_tp_phi_plus_psi_)) ||
             ((perturbations_module_->has_source_psi_) &&
              (index_tp == perturbations_module_->index_tp_psi_)))) {
          owned_sources[index_md][index_ic * perturbations_module_->tp_size_[index_md] + index_tp]
              .resize(perturbations_module_->k_size_[index_md] * perturbations_module_->tau_size_);
          sources[index_md][index_ic * perturbations_module_->tp_size_[index_md] + index_tp] =
              owned_sources[index_md]
                           [index_ic * perturbations_module_->tp_size_[index_md] + index_tp]
                               .data();

          for (int index_tau = 0; index_tau < perturbations_module_->tau_size_; index_tau++) {
            for (int index_k = 0; index_k < perturbations_module_->k_size_[index_md]; index_k++) {
              if (((perturbations_module_->has_source_delta_cb_) &&
                   (index_tp == perturbations_module_->index_tp_delta_cb_)) ||
                  ((perturbations_module_->has_source_theta_cb_) &&
                   (index_tp == perturbations_module_->index_tp_theta_cb_))) {
                sources[index_md][index_ic * perturbations_module_->tp_size_[index_md] +
                                  index_tp][index_tau * perturbations_module_->k_size_[index_md] +
                                            index_k] =
                    perturbations_module_
                        ->sources_[index_md]
                                  [index_ic * perturbations_module_->tp_size_[index_md] + index_tp]
                                  [index_tau * perturbations_module_->k_size_[index_md] + index_k] *
                    nonlinear_module_
                        ->nl_corr_density_[nonlinear_module_->index_pk_cb_]
                                          [index_tau * perturbations_module_->k_size_[index_md] +
                                           index_k];
              }
              else {
                sources[index_md][index_ic * perturbations_module_->tp_size_[index_md] +
                                  index_tp][index_tau * perturbations_module_->k_size_[index_md] +
                                            index_k] =
                    perturbations_module_
                        ->sources_[index_md]
                                  [index_ic * perturbations_module_->tp_size_[index_md] + index_tp]
                                  [index_tau * perturbations_module_->k_size_[index_md] + index_k] *
                    nonlinear_module_
                        ->nl_corr_density_[nonlinear_module_->index_pk_m_]
                                          [index_tau * perturbations_module_->k_size_[index_md] +
                                           index_k];
              }
            }
          }
        }
        else {
          // The table mixes writable local copies (nonlinear case above) with
          // module-owned tables used read-only; the cast is confined here.
          sources[index_md][index_ic * perturbations_module_->tp_size_[index_md] + index_tp] =
              const_cast<double*>(
                  perturbations_module_
                      ->sources_[index_md]
                                [index_ic * perturbations_module_->tp_size_[index_md] + index_tp]
                      .data());
        }
      }
    }
  }
}

void TransferModule::transfer_perturbation_source_spline(
    const std::vector<std::vector<double*>>& sources,
    std::vector<std::vector<double*>>& sources_spline_ptrs,
    std::vector<std::vector<std::vector<double>>>& sources_spline_storage) {
  for (int index_md = 0; index_md < md_size_; index_md++) {
    sources_spline_ptrs[index_md].resize(perturbations_module_->ic_size_[index_md] *
                                         perturbations_module_->tp_size_[index_md]);
    sources_spline_storage[index_md].resize(perturbations_module_->ic_size_[index_md] *
                                            perturbations_module_->tp_size_[index_md]);

    for (int index_ic = 0; index_ic < perturbations_module_->ic_size_[index_md]; index_ic++) {
      for (int index_tp = 0; index_tp < perturbations_module_->tp_size_[index_md]; index_tp++) {
        sources_spline_storage[index_md][index_ic * perturbations_module_->tp_size_[index_md] +
                                         index_tp]
            .resize(perturbations_module_->k_size_[index_md] * perturbations_module_->tau_size_);
        sources_spline_ptrs[index_md][index_ic * perturbations_module_->tp_size_[index_md] +
                                      index_tp] =
            sources_spline_storage[index_md]
                                  [index_ic * perturbations_module_->tp_size_[index_md] + index_tp]
                                      .data();

        array_spline_table_columns2(
            perturbations_module_->k_[index_md].data(),
            perturbations_module_->k_size_[index_md],
            sources[index_md][index_ic * perturbations_module_->tp_size_[index_md] + index_tp],
            perturbations_module_->tau_size_,
            sources_spline_ptrs[index_md]
                               [index_ic * perturbations_module_->tp_size_[index_md] + index_tp],
            _SPLINE_EST_DERIV_);
      }
    }
  }
}

/**
 * This routine defines the number and values of multipoles l for all modes.
 *
 * @return the error status
 */

void TransferModule::transfer_get_l_list() {
  int index_l;
  int l_max = 0;
  int increment;
  /** Summary: */

  /* check that largests need value of l_max */

  if (ppt->has_cls) {
    if (ppt->has_scalars) {
      if ((ppt->has_cl_cmb_temperature) || (ppt->has_cl_cmb_polarization) ||
          (ppt->has_cl_cmb_lensing_potential)) {
        l_max = std::max(ppt->l_scalar_max, l_max);
      }

      if ((ppt->has_cl_lensing_potential) || (ppt->has_cl_number_count)) {
        l_max = std::max(ppt->l_lss_max, l_max);
      }
    }

    if (ppt->has_vectors) {
      l_max = std::max(ppt->l_vector_max, l_max);
    }

    if (ppt->has_tensors) {
      l_max = std::max(ppt->l_tensor_max, l_max);
    }
  }

  /** - allocate and fill l array */

  /** - start from l = 2 and increase with logarithmic step */

  l_.clear();
  l_.push_back(2);
  increment = std::max(static_cast<int>(
                           l_.back() *
                           (pow(ppr->l_logstep, thermodynamics_module_->angular_rescaling_) - 1.)),
                       1);

  while (((l_.back() + increment) < l_max) &&
         (increment < ppr->l_linstep * thermodynamics_module_->angular_rescaling_)) {
    l_.push_back(l_.back() + increment);
    increment =
        std::max(static_cast<int>(
                     l_.back() *
                     (pow(ppr->l_logstep, thermodynamics_module_->angular_rescaling_) - 1.)),
                 1);
  }

  /** - when the logarithmic step becomes larger than some linear step,
      stick to this linear step till l_max */

  increment = ppr->l_linstep * thermodynamics_module_->angular_rescaling_;

  while ((l_.back() + increment) <= l_max) {
    l_.push_back(l_.back() + increment);
  }

  /** - last value set to exactly l_max */

  if (l_.back() != l_max) {
    l_.push_back(l_max);
  }

  l_size_max_ = static_cast<int>(l_.size());

  /* for each mode and type, find relevant size of l array,
     l_size_tt[index_md][index_tt] (since for some modes and types
     l_max can be smaller). Also, maximize this size for each mode to
     find l_size[index_md]. */

  for (int index_md = 0; index_md < perturbations_module_->md_size_; index_md++) {
    l_size_[index_md] = 0;

    for (int index_tt = 0; index_tt < tt_size_[index_md]; index_tt++) {
      l_max = types_[index_md][index_tt].l_max;

      class_test(l_max > l_[l_size_max_ - 1],
                 "For mode %d, type %d, asked for l_max=%d greater than in Bessel table where "
                 "l_max=%d",
                 index_md,
                 index_tt,
                 l_max,
                 l_[l_size_max_ - 1]);

      index_l = 0;
      while (l_[index_l] < l_max) {
        index_l++;
      }
      l_size_tt_[index_md][index_tt] = index_l + 1;

      if (l_size_tt_[index_md][index_tt] < l_size_max_) {
        l_size_tt_[index_md][index_tt]++;
      }
      if (l_size_tt_[index_md][index_tt] < l_size_max_) {
        l_size_tt_[index_md][index_tt]++;
      }

      l_size_[index_md] = std::max(l_size_[index_md], l_size_tt_[index_md][index_tt]);
    }
  }
}

/**
 * This routine defines the number and values of wavenumbers q for
 * each mode (goes smoothly from logarithmic step for small q's to
 * linear step for large q's).
 *
 * @param q_period Input: order of magnitude of the oscillation period of transfer functions
 * @param K        Input: spatial curvature (in absolute value)
 * @param sgnK     Input: spatial curvature sign (open/closed/flat)
 * @return the error status
 */

void TransferModule::transfer_get_q_list(double q_period, double K, int sgnK) {
  int index_q;
  double q;
  double q_min = 0.;
  double q_max = 0.;
  double q_step;
  int nu;
  int nu_min;
  int nu_proposed;
  double q_approximation;
  double last_step = 0.;
  int last_index   = 0;
  double q_logstep_spline;
  double q_logstep_trapzd;

  /* first and last value in flat case*/

  if (sgnK == 0) {
    q_min = perturbations_module_->k_min_;

    q_max = 0.;
    for (int index_md = 0; index_md < perturbations_module_->md_size_; index_md++) {
      q_max = std::max(q_max,
                       perturbations_module_
                           ->k_[index_md][perturbations_module_->k_size_cl_[index_md] - 1]);
    }

    K = 0;
  }

  /* first and last value in open case*/

  else if (sgnK == -1) {
    q_min = sqrt(perturbations_module_->k_min_ * perturbations_module_->k_min_ + K);

    double k_max = 0.;
    for (int index_md = 0; index_md < perturbations_module_->md_size_; index_md++) {
      k_max = std::max(k_max,
                       perturbations_module_
                           ->k_[index_md][perturbations_module_->k_size_cl_[index_md] - 1]);
    }

    q_max = sqrt(k_max * k_max + K);
    if (ppt->has_vectors) {
      q_max = std::min(q_max, sqrt(k_max * k_max + 2. * K));
    }
    if (ppt->has_tensors) {
      q_max = std::min(q_max, sqrt(k_max * k_max + 3. * K));
    }
  }

  /* first and last value in closed case*/

  else if (sgnK == 1) {
    nu_min = 3;
    q_min  = nu_min * sqrt(K);

    q_max = 0.;
    for (int index_md = 0; index_md < perturbations_module_->md_size_; index_md++) {
      q_max = std::max(q_max,
                       perturbations_module_
                           ->k_[index_md][perturbations_module_->k_size_cl_[index_md] - 1]);
    }
  }

  /* adjust the parameter governing the log step size to curvature */

  q_logstep_spline = ppr->q_logstep_spline /
                     pow(thermodynamics_module_->angular_rescaling_, ppr->q_logstep_open);
  q_logstep_trapzd = ppr->q_logstep_trapzd;

  /* assign the first value before starting the loop */

  q_.clear();
  q_.push_back(q_min);
  nu = 3;

  /* loop over the values */

  while (q_.back() < q_max) {
    index_q = static_cast<int>(q_.size());

    /* cap at 2^25 (~3*10^7) values so that pathological precision
       parameters fail cleanly instead of sampling forever */
    class_test(index_q >= (1 << 25), "buggy q-list definition (q_size=%d)", index_q);

    /* step size formula in flat/open case. Step goes gradually from
       logarithmic to linear:

       - in the small q limit, it is logarithmic with: (delta q / q) =
       q_period * q_logstep_spline

       - in the large q limit, it is linear with: (delta q) = q_period
       * ppr->q_linstep
       */

    if (sgnK <= 0) {
      q = q_[index_q - 1] + q_period * ppr->q_linstep * q_[index_q - 1] /
                                (q_[index_q - 1] + ppr->q_linstep / q_logstep_spline);
    }

    /* step size formula in closed case. Same thing excepted that:

       - in the small q limit, the logarithmic step is reduced, being
       given by q_logstep_trapzd, and values are rounded to integer
       values of nu=q/sqrt(K). This happens as long as
       nu<nu_flat_approximation

       - for nu>nu_flat_approximation, the step gradually catches up
       the same expression as in the flat/open case, and there is no
       need to round up to integer nu's.
    */

    else {
      if (nu < static_cast<int>(ppr->hyper_flat_approximation_nu)) {
        q = q_[index_q - 1] + q_period * ppr->q_linstep * q_[index_q - 1] /
                                  (q_[index_q - 1] + ppr->q_linstep / q_logstep_trapzd);

        nu_proposed = static_cast<int>(q / sqrt(K));
        if (nu_proposed <= nu + 1) {
          nu = nu + 1;
        }
        else {
          nu = nu_proposed;
        }

        q          = nu * sqrt(K);
        last_step  = q - q_[index_q - 1];
        last_index = index_q + 1;
      }
      else {
        q_step = q_period * ppr->q_linstep * q_[index_q - 1] /
                 (q_[index_q - 1] + ppr->q_linstep / q_logstep_spline);

        if (index_q - last_index < static_cast<int>(ppr->q_numstep_transition)) {
          q = q_[index_q - 1] +
              (1 - static_cast<double>(index_q - last_index) / ppr->q_numstep_transition) *
                  last_step +
              static_cast<double>(index_q - last_index) / ppr->q_numstep_transition * q_step;
        }
        else {
          q = q_[index_q - 1] + q_step;
        }
      }
    }

    q_.push_back(q);
  }

  /* drop the last point if we overshot q_max, then infer total number of values */

  if (q_.back() > q_max) {
    q_.pop_back();
  }

  q_size_ = static_cast<int>(q_.size());

  class_test(q_size_ < 2, "buggy q-list definition");

  /* in curved universe, check at which index the flat rescaling
     approximation will start being used */

  if (sgnK != 0) {
    q_approximation = ppr->hyper_flat_approximation_nu * sqrt(sgnK * K);
    if (q_approximation < q_[q_size_ - 1]) {
      for (index_q_flat_approximation_ = 0; index_q_flat_approximation_ < q_size_ - 1;
           index_q_flat_approximation_++) {
        if (q_[index_q_flat_approximation_] > q_approximation) {
          break;
        }
      }
      if (ptr->transfer_verbose > 1) {
        printf(
            "Flat bessel approximation spares hyperspherical bessel computations for %d "
            "wavenumbers over a total of %d\n",
            q_size_ - index_q_flat_approximation_,
            q_size_);
      }
    }
    else {
      /* All q values are below the flat approximation threshold, so
         no flat Bessel approximation is used for any q point. Set
         index past the end of the q grid so no valid q index enters
         the flat-approximation branch. Note: consumers that use this
         as an array index (e.g. spectra_compute_cl) must clamp to
         < q_size_. */
      index_q_flat_approximation_ = q_size_;
      if (ptr->transfer_verbose > 1) {
        printf(
            "Cannot use flat bessel approximation, since q_max = %.5e < q_start_approximation = "
            "%.5e\n",
            q_[q_size_ - 1],
            q_approximation);
      }
    }
  }
}

/**
 * Define the coarse logarithmic q-grid used only by the full-Limber CMB
 * lensing calculation. Its upper limit follows the complete perturbation
 * source grid, not k_size_cl.
 */
/**
 * Lower anchor of the full-Limber CMB lensing grid.
 *
 * @param k_min      Input: smallest wavenumber solved by the perturbation module
 * @param k_min_flat Input: the flat truncation scale, k_min_tau0/tau0
 * @param sgnK       Input: sign of the spatial curvature
 * @param K          Input: spatial curvature
 * @return the q value at which the geometric grid starts
 */
double TransferModule::LimberGridQMin(double k_min, double k_min_flat, int sgnK, double K) {
  if (sgnK == 0) {
    return k_min;
  }

  /* Closed: the lowest scalar mode is nu = 3. Never anchor below the scale at
     which the flat run truncates, or the grid slides with K instead of
     converging to the flat one. Written as a pure clamp of the existing
     3 sqrt(K), so that every curvature above the crossover is untouched. */
  if (sgnK == 1) {
    return std::max(3. * sqrt(K), k_min_flat);
  }

  /* Open: sqrt(k_min^2 + K) with k_min = sqrt(-K + k_min_flat^2) already
     collapses to the flat anchor, so this converges on its own. */
  return sqrt(k_min * k_min + K);
}

void TransferModule::transfer_get_q_limber_list(double K, int sgnK) {
  double q_min = 0.;
  double q_max = 0.;

  class_test(ppr->q_logstep_limber <= 1.,
             "q_logstep_limber=%e should be strictly larger than one",
             ppr->q_logstep_limber);

  q_min = LimberGridQMin(perturbations_module_->k_min_,
                         ppr->k_min_tau0 / background_module_->conformal_age_,
                         sgnK,
                         K);

  if (sgnK == 0) {
    for (int index_md = 0; index_md < perturbations_module_->md_size_; index_md++) {
      q_max = std::max(q_max,
                       perturbations_module_
                           ->k_[index_md][perturbations_module_->k_size_[index_md] - 1]);
    }
  }
  else if (sgnK == -1) {
    double k_max = 0.;
    for (int index_md = 0; index_md < perturbations_module_->md_size_; index_md++) {
      k_max = std::max(k_max,
                       perturbations_module_
                           ->k_[index_md][perturbations_module_->k_size_[index_md] - 1]);
    }
    q_max = sqrt(k_max * k_max + K);
    if (ppt->has_vectors) {
      q_max = std::min(q_max, sqrt(k_max * k_max + 2. * K));
    }
    if (ppt->has_tensors) {
      q_max = std::min(q_max, sqrt(k_max * k_max + 3. * K));
    }
  }
  else {
    for (int index_md = 0; index_md < perturbations_module_->md_size_; index_md++) {
      q_max = std::max(q_max,
                       perturbations_module_
                           ->k_[index_md][perturbations_module_->k_size_[index_md] - 1]);
    }
  }

  class_test(q_max <= q_min, "full-Limber q range is invalid: q_min=%e, q_max=%e", q_min, q_max);

  q_size_limber_ = static_cast<int>(log(q_max / q_min) / log(ppr->q_logstep_limber)) + 1;
  class_test(q_size_limber_ < 2,
             "full-Limber q grid contains fewer than two points: q_min=%e, q_max=%e",
             q_min,
             q_max);

  q_limber_.resize(q_size_limber_);
  q_limber_[0] = q_min;
  for (int index_q = 1; index_q < q_size_limber_; index_q++) {
    q_limber_[index_q] = ppr->q_logstep_limber * q_limber_[index_q - 1];
  }
}

/**
 * This routine infers from the q values a list of corresponding k
 * values for each mode.
 *
 * @param K       Input: spatial curvature
 * @return the error status
 */

void TransferModule::transfer_get_k_list(double K) {
  double m = 0.;

  k_.resize(md_size_);
  if (do_lcmb_full_limber_) {
    k_limber_.resize(md_size_);
  }

  for (int index_md = 0; index_md < md_size_; index_md++) {
    k_[index_md].resize(q_size_);
    if (do_lcmb_full_limber_) {
      k_limber_[index_md].resize(q_size_limber_);
    }

    if (_scalarsEXT_) {
      m = 0.;
    }
    if (_vectorsEXT_) {
      m = 1.;
    }
    if (_tensorsEXT_) {
      m = 2.;
    }

    for (int index_q = 0; index_q < q_size_; index_q++) {
      k_[index_md][index_q] = sqrt(q_[index_q] * q_[index_q] - K * (m + 1.));
    }
    if (do_lcmb_full_limber_) {
      for (int index_q = 0; index_q < q_size_limber_; index_q++) {
        k_limber_[index_md][index_q] = sqrt(q_limber_[index_q] * q_limber_[index_q] - K * (m + 1.));
      }
    }

    if (k_[index_md][0] < perturbations_module_->k_[index_md][0]) {
      /* If k_[index_md][0] < perturbations_module_->k_[index_md][0] at the level of rounding,
         adjust first value of k_list to avoid interpolation errors: */
      if ((perturbations_module_->k_[index_md][0] - k_[index_md][0]) < 10. * DBL_EPSILON) {
        k_[index_md][0] = perturbations_module_->k_[index_md][0];
      }
      else {
        class_stop(
            "bug in k_list calculation: in perturbation module k_min=%e, in transfer module "
            "k_min[mode=%d]=%e, interpolation impossible",
            perturbations_module_->k_[0][0],
            index_md,
            k_[index_md][0]);
      }
    }

    if (do_lcmb_full_limber_ && k_limber_[index_md][0] < perturbations_module_->k_[index_md][0]) {
      if ((perturbations_module_->k_[index_md][0] - k_limber_[index_md][0]) < 10. * DBL_EPSILON) {
        k_limber_[index_md][0] = perturbations_module_->k_[index_md][0];
      }
      else {
        class_stop(
            "bug in full-Limber k-list calculation: in perturbation module k_min=%e, in "
            "transfer module k_min[mode=%d]=%e, interpolation impossible",
            perturbations_module_->k_[0][0],
            index_md,
            k_limber_[index_md][0]);
      }
    }

    class_test(k_[index_md][q_size_ - 1] >
                   perturbations_module_->k_[0][perturbations_module_->k_size_cl_[0] - 1],
               "bug in k_list calculation: in perturbation module k_max=%e, in transfer module "
               "k_max[mode=%d]=%e, interpolation impossible",
               perturbations_module_->k_[0][perturbations_module_->k_size_cl_[0]],
               index_md,
               k_[index_md][q_size_ - 1]);

    if (do_lcmb_full_limber_) {
      class_test(k_limber_[index_md][q_size_limber_ - 1] >
                     perturbations_module_->k_[0][perturbations_module_->k_size_[0] - 1],
                 "bug in full-Limber k-list calculation: in perturbation module k_max=%e, in "
                 "transfer "
                 "module k_max[mode=%d]=%e, interpolation impossible",
                 perturbations_module_->k_[0][perturbations_module_->k_size_[0] - 1],
                 index_md,
                 k_limber_[index_md][q_size_limber_ - 1]);
    }
  }
}

/**
 * Define once for all the properties of each transfer type: perturbation
 * source, redshift bin, largest multipole, radial function, window, power
 * of k, Limber threshold and neglect rules. Every per-type decision in
 * this module reads these properties.
 */

void TransferModule::transfer_type_properties() {
  const int nbin = ppt->selection_num;

  types_.assign(md_size_, {});

  for (int index_md = 0; index_md < md_size_; index_md++) {
    types_[index_md].resize(tt_size_[index_md]);

    for (int index_tt = 0; index_tt < tt_size_[index_md]; index_tt++) {
      TransferType& type = types_[index_md][index_tt];

      /* is index_tt one of the nbin types starting at index_first? Then set its term and bin. */
      auto in_bins = [&](int index_first, bool flag, TransferTerm term) {
        if (flag && (index_tt >= index_first) && (index_tt < index_first + nbin)) {
          type.term = term;
          type.bin  = index_tt - index_first;
        }
      };

      if (_scalarsEXT_) {
        if (ppt->has_cl_cmb_temperature) {
          if (index_tt == index_tt_t0_) {
            type.term = TransferTerm::cmb_t0;
          }
          if (index_tt == index_tt_t1_) {
            type.term = TransferTerm::cmb_t1;
          }
          if (index_tt == index_tt_t2_) {
            type.term = TransferTerm::cmb_t2;
          }
        }
        if (ppt->has_cl_cmb_polarization && (index_tt == index_tt_e_)) {
          type.term = TransferTerm::cmb_e;
        }
        if (ppt->has_cl_cmb_lensing_potential && (index_tt == index_tt_lcmb_)) {
          type.term = TransferTerm::cmb_lensing;
        }
        in_bins(index_tt_density_, ppt->has_nc_density, TransferTerm::nc_density);
        in_bins(index_tt_rsd_, ppt->has_nc_rsd, TransferTerm::nc_rsd);
        in_bins(index_tt_d0_, ppt->has_nc_rsd, TransferTerm::nc_d0);
        in_bins(index_tt_d1_, ppt->has_nc_rsd, TransferTerm::nc_d1);
        in_bins(index_tt_nc_g1_, ppt->has_nc_gr, TransferTerm::nc_g1);
        in_bins(index_tt_nc_g2_, ppt->has_nc_gr, TransferTerm::nc_g2);
        in_bins(index_tt_nc_g3_, ppt->has_nc_gr, TransferTerm::nc_g3);
        in_bins(index_tt_nc_lens_, ppt->has_nc_lens, TransferTerm::nc_lens);
        in_bins(index_tt_nc_g4_, ppt->has_nc_gr, TransferTerm::nc_g4);
        in_bins(index_tt_nc_g5_, ppt->has_nc_gr, TransferTerm::nc_g5);
        in_bins(index_tt_lensing_, ppt->has_cl_lensing_potential, TransferTerm::galaxy_lensing);
      }
      else {
        if (ppt->has_cl_cmb_temperature) {
          if (_vectorsEXT_ && (index_tt == index_tt_t1_)) {
            type.term = TransferTerm::cmb_t1;
          }
          if (index_tt == index_tt_t2_) {
            type.term = TransferTerm::cmb_t2;
          }
        }
        if (ppt->has_cl_cmb_polarization) {
          if (index_tt == index_tt_e_) {
            type.term = TransferTerm::cmb_e;
          }
          if (index_tt == index_tt_b_) {
            type.term = TransferTerm::cmb_b;
          }
        }
      }

      /* the transfer function can be neglected for l < (q-delta_k_neglect) r_a(tau_rec) */
      auto cmb = [&](int index_tp, radial_function_type radial, double delta_k, bool late_cut) {
        type.index_tp        = index_tp;
        type.radial_type     = radial;
        type.delta_k_neglect = delta_k;
        type.late_source_cut = late_cut;
      };
      /* number counts and galaxy lensing: the Limber threshold scales with the mean redshift of the bin */
      auto lss = [&](int index_tp,
                     radial_function_type radial,
                     TransferWindow window,
                     double k_power,
                     double l_limber_over_z) {
        type.index_tp    = index_tp;
        type.radial_type = radial;
        type.window_type = window;
        type.k_power     = k_power;
        type.l_limber    = l_limber_over_z * ppt->selection_mean[type.bin];
        type.l_max       = ppt->l_lss_max;
      };
      const auto* pt        = perturbations_module_.get();
      const double nc_local = ppr->l_switch_limber_for_nc_local_over_z;
      const double nc_los   = ppr->l_switch_limber_for_nc_los_over_z;
      const auto selection  = TransferWindow::selection;
      const auto integrated = TransferWindow::integrated;

      if (_scalarsEXT_) {
        type.l_max = ppt->l_scalar_max;
        switch (type.term) {
          /* sources at late times can be neglected at high l for CMB
             types, excepted T0 which contains the late ISW effect */
          case TransferTerm::cmb_t0:
            cmb(pt->index_tp_t0_, SCALAR_TEMPERATURE_0, ppr->transfer_neglect_delta_k_S_t0, false);
            break;
          case TransferTerm::cmb_t1:
            cmb(pt->index_tp_t1_, SCALAR_TEMPERATURE_1, ppr->transfer_neglect_delta_k_S_t1, true);
            break;
          case TransferTerm::cmb_t2:
            cmb(pt->index_tp_t2_, SCALAR_TEMPERATURE_2, ppr->transfer_neglect_delta_k_S_t2, true);
            break;
          case TransferTerm::cmb_e:
            cmb(pt->index_tp_p_, SCALAR_POLARISATION_E, ppr->transfer_neglect_delta_k_S_e, true);
            break;
          case TransferTerm::cmb_b:
            break;
          /* CMB lensing potential: an integrated type with support [tau_rec, tau0] */
          case TransferTerm::cmb_lensing:
            type.index_tp    = pt->index_tp_phi_plus_psi_;
            type.window_type = TransferWindow::cmb_lensing;
            type.k_power     = ptr->lcmb_tilt;
            type.k_pivot     = ptr->lcmb_pivot;
            type.l_limber    = floor(ppr->l_switch_limber) + 1.;
            type.full_limber = do_lcmb_full_limber_;
            break;
          /* use here delta_cb rather than delta_m if density number counts calculated only for cold dark matter + baryon */
          /* (this important comment is referenced in a WARNING message in perturbations.c) */
          case TransferTerm::nc_density:
            lss(pt->index_tp_delta_m_, SCALAR_TEMPERATURE_0, selection, 0., nc_local);
            break;
          case TransferTerm::nc_rsd:
            lss(pt->index_tp_theta_m_, NC_RSD, selection, 0., nc_local);
            break;
          /* factors 1/k^2 and 1/k from the original CLASSgal paper (arXiv:1307.1459) */
          case TransferTerm::nc_d0:
            lss(pt->index_tp_theta_m_, SCALAR_TEMPERATURE_0, selection, -2., nc_local);
            break;
          case TransferTerm::nc_d1:
            lss(pt->index_tp_theta_m_, SCALAR_TEMPERATURE_1, selection, -1., nc_local);
            break;
          case TransferTerm::nc_g1:
            lss(pt->index_tp_psi_, SCALAR_TEMPERATURE_0, selection, 0., nc_local);
            break;
          case TransferTerm::nc_g2:
            lss(pt->index_tp_phi_, SCALAR_TEMPERATURE_0, selection, 0., nc_local);
            break;
          case TransferTerm::nc_g3:
            lss(pt->index_tp_phi_prime_, SCALAR_TEMPERATURE_0, selection, 0., nc_local);
            break;
          case TransferTerm::nc_lens:
            lss(pt->index_tp_phi_plus_psi_, SCALAR_TEMPERATURE_0, integrated, 0., nc_los);
            break;
          case TransferTerm::nc_g4:
            lss(pt->index_tp_phi_plus_psi_, SCALAR_TEMPERATURE_0, integrated, 0., nc_los);
            break;
          /* factor k from the chi derivative d/dchi j_l(k chi) = k j_l'(k chi). Even if G5 is
             integrated along the line-of-sight, it has the Limber threshold of local terms,
             because it involves the derivative of the Bessel function */
          case TransferTerm::nc_g5:
            lss(pt->index_tp_phi_plus_psi_, SCALAR_TEMPERATURE_1, integrated, 1., nc_local);
            break;
          case TransferTerm::galaxy_lensing:
            lss(pt->index_tp_phi_plus_psi_, SCALAR_TEMPERATURE_0, integrated, 0., nc_los);
            break;
        }

        /* number counts with Dirac selection functions never use the Limber approximation */
        if ((type.bin >= 0) && (ppt->selection == dirac) &&
            (type.term != TransferTerm::galaxy_lensing)) {
          type.l_limber = _HUGE_;
        }
      }

      if (_vectorsEXT_) {
        type.l_max = ppt->l_vector_max;
        switch (type.term) {
          case TransferTerm::cmb_t1:
            cmb(pt->index_tp_t1_, VECTOR_TEMPERATURE_1, ppr->transfer_neglect_delta_k_V_t1, true);
            break;
          case TransferTerm::cmb_t2:
            cmb(pt->index_tp_t2_, VECTOR_TEMPERATURE_2, ppr->transfer_neglect_delta_k_V_t2, true);
            break;
          case TransferTerm::cmb_e:
            cmb(pt->index_tp_p_, VECTOR_POLARISATION_E, ppr->transfer_neglect_delta_k_V_e, true);
            break;
          case TransferTerm::cmb_b:
            cmb(pt->index_tp_p_, VECTOR_POLARISATION_B, ppr->transfer_neglect_delta_k_V_b, true);
            break;
          default:
            break;
        }
      }

      /* the late time sources of tensor temperature can never be neglected */
      if (_tensorsEXT_) {
        type.l_max = ppt->l_tensor_max;
        switch (type.term) {
          case TransferTerm::cmb_t2:
            cmb(pt->index_tp_t2_, TENSOR_TEMPERATURE_2, ppr->transfer_neglect_delta_k_T_t2, false);
            break;
          case TransferTerm::cmb_e:
            cmb(pt->index_tp_p_, TENSOR_POLARISATION_E, ppr->transfer_neglect_delta_k_T_e, true);
            break;
          case TransferTerm::cmb_b:
            cmb(pt->index_tp_p_, TENSOR_POLARISATION_B, ppr->transfer_neglect_delta_k_T_b, true);
            break;
          default:
            break;
        }
      }

      class_test(type.index_tp < 0,
                 "no perturbation source associated to transfer type %d of mode %d",
                 index_tt,
                 index_md);
    }
  }
}

/**
 * For each selection function, compute once for all the min, mean and
 * max values of conformal time (associated to the min, mean and max
 * values of redshift specified by the user)
 */

void TransferModule::transfer_selection_times() {
  const int nbin = ppt->selection_num;
  selection_tau_min_.assign(nbin, 0.);
  selection_tau_mean_.assign(nbin, 0.);
  selection_tau_max_.assign(nbin, 0.);

  for (int bin = 0; bin < nbin; bin++) {
    /* redshift interval of the bin; the few lines below should be consistent with their
       counterpart in the input module */
    double half_width = 0.;
    if (ppt->selection == gaussian) {
      half_width = ppt->selection_width[bin] * ppr->selection_cut_at_sigma;
    }
    if (ppt->selection == tophat) {
      half_width = (1. + ppr->selection_cut_at_sigma * ppr->selection_tophat_edge) *
                   ppt->selection_width[bin];
    }
    double z_high = ppt->selection_mean[bin] + half_width;
    double z_low  = std::max(ppt->selection_mean[bin] - half_width, 0.);
    if (ppt->selection == tabulated) {
      z_high = ptr->selection_z_max[bin];
      z_low  = ptr->selection_z_min[bin];
    }

    /* lower edge of time interval, higher edge, and central value */
    background_module_->background_tau_of_z(z_high, &selection_tau_min_[bin]);
    background_module_->background_tau_of_z(z_low, &selection_tau_max_[bin]);
    background_module_->background_tau_of_z(std::max(ppt->selection_mean[bin], 0.),
                                            &selection_tau_mean_[bin]);
  }
}

/**
 * Define once for all the time sampling of the line-of-sight integral of
 * each transfer type. Types with an identical sampling share it. For each
 * sampling, store the values of (tau0-tau), the weights of the trapezoidal
 * rule, and the indices and weights for resampling the perturbation
 * sources by linear interpolation.
 *
 * The sampling must catch both the time dependence of the transfer source
 * W(tau) S(k, tau) and the oscillations of the Bessel functions. CMB
 * types use the times of the perturbation module. Number count and galaxy
 * lensing types use evenly spaced times: at least ppr->selection_sampling,
 * and enough to sample the oscillations of \f$ j_l(k (\tau_0-\tau)) \f$
 * (period \f$ \Delta \tau = 2\pi/k \f$) for the largest relevant wavenumber,
 * k_max = l_limber/(tau0-tau_mean) for local types and
 * l_limber/((tau0-tau_mean)/2) for integrated ones.
 *
 * @param tau_rec Input: recombination time
 */

void TransferModule::transfer_time_samplings(double tau_rec) {
  const double tau0                 = background_module_->conformal_age_;
  const std::vector<double>& tau_pt = perturbations_module_->tau_sampling_;
  const int tau_size_pt             = perturbations_module_->tau_size_;

  /* description of a time sampling: two types share a sampling if these are equal */
  struct Spec {
    TransferSampling kind;
    double tau_min;
    double tau_max;
    int tau_size;
    double log_step;
    bool operator==(const Spec& o) const {
      return kind == o.kind && tau_min == o.tau_min && tau_max == o.tau_max &&
             tau_size == o.tau_size && log_step == o.log_step;
    }
  };

  class_test(ppr->selection_sampling_bessel_min <= 0.,
             "selection_sampling_bessel_min=%e should be strictly positive (it is a minimum number "
             "of "
             "sampled points per radian of the Bessel functions)",
             ppr->selection_sampling_bessel_min);

  samplings_.assign(md_size_, {});
  tau_size_max_  = 0;
  tt_per_ts_max_ = 0;

  for (int index_md = 0; index_md < md_size_; index_md++) {
    std::vector<Spec> specs;

    for (int index_tt = 0; index_tt < tt_size_[index_md]; index_tt++) {
      TransferType& type = types_[index_md][index_tt];
      const int bin      = type.bin;
      Spec spec{TransferSampling::from_perturbations, 0., tau0, 0, 0.};

      switch (type.window_type) {
        /* CMB temperature and polarisation: all times of the perturbation module */
        case TransferWindow::none:
          break;

        /* CMB lensing potential: logarithmic late-time sampling after recombination (class_public #267) */
        case TransferWindow::cmb_lensing:
          spec = {TransferSampling::logarithmic,
                  tau_rec,
                  (1. - ppr->cmb_lensing_sampling_fraction) * tau0,
                  0,
                  ppr->cmb_lensing_sampling_factor /
                      std::min(ppr->l_switch_limber, static_cast<double>(ppt->l_scalar_max))};
          class_test(ppr->cmb_lensing_sampling_factor <= 0.,
                     "cmb_lensing_sampling_factor=%e should be strictly positive",
                     ppr->cmb_lensing_sampling_factor);
          class_test((ppr->cmb_lensing_sampling_fraction <= 0.) ||
                         (ppr->cmb_lensing_sampling_fraction >= 1.),
                     "cmb_lensing_sampling_fraction=%e should be strictly between 0 and 1",
                     ppr->cmb_lensing_sampling_fraction);
          class_test(
              ppr->cmb_lensing_sampling_fraction * spec.log_step <
                  10. * ppr->smallest_allowed_variation,
              "cmb_lensing_sampling_fraction=%e and/or cmb_lensing_sampling_factor=%e are too "
              "small, and/or l_switch_limber=%g and/or l_scalar_max=%d too big: the logarithmic "
              "step sampling the CMB lensing source would reach machine precision",
              ppr->cmb_lensing_sampling_fraction,
              ppr->cmb_lensing_sampling_factor,
              ppr->l_switch_limber,
              ppt->l_scalar_max);
          break;

        /* local number count contributions: the support of the selection function of the bin,
           with one point for a Dirac selection function */
        case TransferWindow::selection: {
          spec = {TransferSampling::uniform,
                  selection_tau_min_[bin],
                  selection_tau_max_[bin],
                  1,
                  0.};
          if (spec.tau_min != spec.tau_max) {
            const int l_limber = ppr->l_switch_limber_for_nc_local_over_z *
                                 ppt->selection_mean[bin];
            spec.tau_size =
                std::max(static_cast<int>(ppr->selection_sampling),
                         static_cast<int>(static_cast<int>((spec.tau_max - spec.tau_min) /
                                                           ((tau0 - selection_tau_mean_[bin]) /
                                                            std::min(l_limber, ppt->l_lss_max))) *
                                          ppr->selection_sampling_bessel));
          }
          break;
        }

        /* number count contributions integrated along the line of sight, and galaxy
           lensing: from the beginning of the selection function of the bin till today */
        case TransferWindow::integrated: {
          spec               = {TransferSampling::uniform, selection_tau_min_[bin], tau0, 0, 0.};
          const bool g5      = (type.term == TransferTerm::nc_g5);
          const int l_limber = (g5 ? ppr->l_switch_limber_for_nc_local_over_z
                                   : ppr->l_switch_limber_for_nc_los_over_z) *
                               ppt->selection_mean[bin];
          spec.tau_size =
              std::max(static_cast<int>(ppr->selection_sampling),
                       static_cast<int>(static_cast<int>((tau0 - spec.tau_min) /
                                                         ((tau0 - selection_tau_mean_[bin]) / 2. /
                                                          std::min(l_limber, ppt->l_lss_max))) *
                                        (g5 ? ppr->selection_sampling_bessel
                                            : ppr->selection_sampling_bessel_los)));
          break;
        }
      }

      if (spec.kind == TransferSampling::from_perturbations) {
        spec.tau_size = static_cast<int>(
            tau_pt.end() - std::upper_bound(tau_pt.begin(), tau_pt.end(), spec.tau_min));
      }

      /* reuse the sampling of a previous type if identical */
      auto found    = std::find(specs.begin(), specs.end(), spec);
      type.index_ts = static_cast<int>(found - specs.begin());
      if (found == specs.end()) {
        specs.push_back(spec);
      }
    }

    /* fill each sampling */
    samplings_[index_md].resize(specs.size());

    for (size_t index_ts = 0; index_ts < specs.size(); index_ts++) {
      const Spec& spec               = specs[index_ts];
      TransferTimeSampling& sampling = samplings_[index_md][index_ts];
      std::vector<double>& t0mt      = sampling.tau0_minus_tau;

      switch (spec.kind) {
        /* a subset of the perturbation times: no interpolation needed */
        case TransferSampling::from_perturbations:
          for (int index_tau_pt = tau_size_pt - spec.tau_size; index_tau_pt < tau_size_pt;
               index_tau_pt++) {
            t0mt.push_back(tau0 - tau_pt[index_tau_pt]);
            sampling.index_tau_pt.push_back(std::min(index_tau_pt, tau_size_pt - 2));
            sampling.weight_tau_pt.push_back((index_tau_pt == tau_size_pt - 1) ? 1. : 0.);
          }
          break;

        /* evenly spaced times */
        case TransferSampling::uniform:
          if (spec.tau_size == 1) {
            t0mt.push_back(tau0 - spec.tau_min);
          }
          else {
            for (int index_tau = 0; index_tau < spec.tau_size - 1; index_tau++) {
              t0mt.push_back(tau0 - spec.tau_min -
                             (static_cast<double>(index_tau)) /
                                 (static_cast<double>(spec.tau_size) - 1.) *
                                 (spec.tau_max - spec.tau_min));
            }
            t0mt.push_back(tau0 - spec.tau_max);
          }
          break;

        case TransferSampling::logarithmic:
          t0mt = transfer_logarithmic_sampling(spec.tau_min, spec.tau_max, spec.log_step);
          break;
      }

      /* bracketing perturbation times by bisection */
      const int tau_size = static_cast<int>(t0mt.size());
      if (spec.kind != TransferSampling::from_perturbations) {
        for (int index_tau = 0; index_tau < tau_size; index_tau++) {
          const double tau = tau0 - t0mt[index_tau];
          int inf;
          array_search_bisect(tau_size_pt, tau_pt.data(), tau, &inf);
          sampling.index_tau_pt.push_back(inf);
          sampling.weight_tau_pt.push_back((tau - tau_pt[inf]) / (tau_pt[inf + 1] - tau_pt[inf]));
        }
      }

      /* for evenly spaced times: largest q at which the Bessel functions j_l(q(tau0-tau))
         are sampled with at least selection_sampling_bessel_min points per radian. (The
         logarithmic sampling is built to resolve them down to tau0.) */
      if ((spec.kind == TransferSampling::uniform) && (tau_size > 1)) {
        sampling.q_max_sampled = (tau_size - 1.) / (spec.tau_max - spec.tau_min) /
                                 ppr->selection_sampling_bessel_min;
      }

      sampling.w_trapz.resize(tau_size);
      array_trapezoidal_mweights(t0mt.data(), tau_size, sampling.w_trapz.data());
      tau_size_max_ = std::max(tau_size_max_, tau_size);
    }

    for (int index_tt = 0; index_tt < tt_size_[index_md]; index_tt++) {
      std::vector<int>& types = samplings_[index_md][types_[index_md][index_tt].index_ts].types;
      types.push_back(index_tt);
      tt_per_ts_max_ = std::max(tt_per_ts_max_, static_cast<int>(types.size()));
    }
  }
}

/**
 * Time sampling of the CMB lensing source (class_public #267): keep the
 * perturbation times after tau_min as long as their step is smaller than
 * (tau0-tau)*log_step; then continue with the logarithmic step
 * Delta tau = (tau0-tau)*log_step until tau_max; then end at tau0. The
 * logarithmic steps resolve the 1/chi lensing kernel and the Bessel
 * functions down to small distances, which the perturbation sampling,
 * designed for the time dependence of the sources, does not.
 *
 * @param tau_min  Input: sampled times are larger than tau_min
 * @param tau_max  Input: largest time sampled with logarithmic steps
 * @param log_step Input: logarithmic step Delta tau/(tau0-tau)
 * @return values of (tau0-tau) in decreasing order
 */

std::vector<double> TransferModule::transfer_logarithmic_sampling(double tau_min,
                                                                  double tau_max,
                                                                  double log_step) {
  const double tau0                 = background_module_->conformal_age_;
  const std::vector<double>& tau_pt = perturbations_module_->tau_sampling_;
  const int tau_size_pt             = perturbations_module_->tau_size_;
  std::vector<double> tau0_minus_tau;

  /* throw away times before tau_min */
  int index_tau_pt = static_cast<int>(std::upper_bound(tau_pt.begin(), tau_pt.end(), tau_min) -
                                      tau_pt.begin());
  class_test((index_tau_pt < 1) || (index_tau_pt >= tau_size_pt),
             "the perturbation sources should be sampled on both sides of tau=%e (index %d of %d)",
             tau_min,
             index_tau_pt,
             tau_size_pt);

  /* use the time steps of the perturbation module as long as they are smaller than (tau0-tau)*log_step */
  while ((index_tau_pt < tau_size_pt) && (tau_pt[index_tau_pt] - tau_pt[index_tau_pt - 1] <
                                          (tau0 - tau_pt[index_tau_pt - 1]) * log_step)) {
    tau0_minus_tau.push_back(tau0 - tau_pt[index_tau_pt]);
    index_tau_pt++;
  }

  /* if we have not reached tau0 already, use the logarithmic step until tau = tau_max, then tau0
     (replacing the last value if it overshoots tau0) */
  if (index_tau_pt < tau_size_pt) {
    double tau = tau_pt[index_tau_pt - 1];
    while (tau < tau_max) {
      tau += (tau0 - tau) * log_step;
      tau0_minus_tau.push_back(tau0 - tau);
    }
    if (tau0_minus_tau.back() > 0.) {
      tau0_minus_tau.push_back(0.);
    }
    else {
      tau0_minus_tau.back() = 0.;
    }
  }

  return tau0_minus_tau;
}

/**
 * Compute the transfer functions \f$ \Delta_l^{X} (q) \f$ for a given
 * wavenumber q, for all modes, initial conditions, types and multipoles.
 *
 * For each mode, interpolate all the needed perturbation sources at k(q).
 * Then, for each time sampling, build the transfer sources of all the types
 * sharing it with transfer_sources(). For each l, choose the method for each
 * of these types and initial conditions with transfer_select_method(),
 * interpolate the needed Bessel functions once for all of them with
 * transfer_bessel_interpolate(), and compute the transfer functions with
 * transfer_limber() or transfer_integrate().
 *
 * @param index_q             Input: index of wavenumber in q_ or q_limber_
 * @param tau_rec             Input: recombination time
 * @param pert_sources        Input: perturbation sources
 * @param pert_sources_spline Input: second derivatives of the perturbation sources with respect to k
 * @param ptw                 Input: workspace
 * @param use_full_limber     Input: whether index_q refers to the full-Limber grid q_limber_
 */

void TransferModule::transfer_compute_for_each_q(int index_q,
                                                 double tau_rec,
                                                 double** const* pert_sources,
                                                 double** const* pert_sources_spline,
                                                 struct transfer_workspace* ptw,
                                                 bool use_full_limber) {
  const double ra_rec = (background_module_->conformal_age_ - tau_rec) *
                        thermodynamics_module_->angular_rescaling_;

  for (int index_md = 0; index_md < md_size_; index_md++) {
    const std::vector<TransferType>& types = types_[index_md];
    const int ic_size                      = perturbations_module_->ic_size_[index_md];
    const int tp_size                      = perturbations_module_->tp_size_[index_md];

    /* in the full Limber scheme, only the types with the full_limber flag are computed */
    auto computed = [&](int index_tt) { return !use_full_limber || types[index_tt].full_limber; };

    const double q = use_full_limber ? q_limber_[index_q] : q_[index_q];
    const double k = use_full_limber ? k_limber_[index_md][index_q] : k_[index_md][index_q];
    const double k_max =
        use_full_limber
            ? perturbations_module_->k_[index_md][perturbations_module_->k_size_[index_md] - 1]
            : perturbations_module_->k_[index_md][perturbations_module_->k_size_cl_[index_md] - 1];
    std::vector<double>& transfer = use_full_limber ? transfer_limber_[index_md]
                                                    : transfer_[index_md];
    const int q_size              = use_full_limber ? q_size_limber_ : q_size_;

    /* if we reached q_max for this mode, there is nothing to be done */
    if (k > k_max) {
      continue;
    }

    /** - interpolate at k all the perturbation sources needed by this mode */

    std::fill(ptw->is_interpolated.begin(), ptw->is_interpolated.end(), 0);
    for (int index_ic = 0; index_ic < ic_size; index_ic++) {
      for (int index_tt = 0; index_tt < tt_size_[index_md]; index_tt++) {
        const int index_source = index_ic * tp_size + types[index_tt].index_tp;
        if (!computed(index_tt) || ptw->is_interpolated[index_source]) {
          continue;
        }
        ptw->interpolated_sources[index_source].resize(perturbations_module_->tau_size_);
        transfer_interpolate_sources(k,
                                     index_md,
                                     pert_sources[index_md][index_source],
                                     pert_sources_spline[index_md][index_source],
                                     ptw->interpolated_sources[index_source].data());
        ptw->is_interpolated[index_source] = 1;
      }
    }

    /** - loop over time samplings */

    for (const TransferTimeSampling& sampling : samplings_[index_md]) {
      const int tt_size_ts = static_cast<int>(sampling.types.size());
      if (std::none_of(sampling.types.begin(), sampling.types.end(), computed)) {
        continue;
      }

      const std::vector<double>& tau0_minus_tau = sampling.tau0_minus_tau;
      const int tau_size                        = static_cast<int>(tau0_minus_tau.size());

      /* radial coordinates r(tau0-tau) and a few related trigonometric functions,
         and the largest q at which the Bessel functions are tabulated */
      double q_max_bessel = 0.;
      if (!use_full_limber) {
        transfer_radial_coordinates(ptw, tau0_minus_tau.data(), tau_size, k);
        if (ptw->sgnK == 0) {
          q_max_bessel = ptw->pBIS->x[ptw->pBIS->x_size - 1] / tau0_minus_tau[0];
        }
        else {
          q_max_bessel = q_[q_size_ - 1];
        }
      }

      /* transfer sources W(tau) S(k,tau) (k/k_pivot)^k_power of the types using this sampling */
      for (int index_ic = 0; index_ic < ic_size; index_ic++) {
        for (int n = 0; n < tt_size_ts; n++) {
          const int index_tt = sampling.types[n];
          if (computed(index_tt)) {
            transfer_sources(ptw->interpolated_sources[index_ic * tp_size +
                                                       types[index_tt].index_tp]
                                 .data(),
                             k,
                             index_md,
                             index_tt,
                             ptw->sources.data() + (index_ic * tt_per_ts_max_ + n) * tau_size_max_);
          }
        }
      }

      for (int index_l = 0; index_l < l_size_[index_md]; index_l++) {
        const double l = static_cast<double>(l_[index_l]);

        /* choose the method for each (initial condition, type) using this sampling, and find
           which Bessel functions are needed, and at how many times */
        const double tau0_minus_tau_min_bessel =
            use_full_limber ? 0. : transfer_bessel_cut(ptw, index_q, index_l, k);
        int bessel_size = 0;
        BesselNeeds needs{false, false, false};

        for (int index_ic = 0; index_ic < ic_size; index_ic++) {
          for (int n = 0; n < tt_size_ts; n++) {
            const int index_tt   = sampling.types[n];
            const int index_n    = index_ic * tt_per_ts_max_ + n;
            ptw->method[index_n] = transfer_select_method(ptw,
                                                          index_q,
                                                          index_md,
                                                          index_tt,
                                                          index_l,
                                                          q,
                                                          ra_rec,
                                                          q_max_bessel,
                                                          sampling,
                                                          ptw->sources.data() +
                                                              index_n * tau_size_max_,
                                                          use_full_limber,
                                                          tau0_minus_tau_min_bessel,
                                                          &ptw->index_tau_max[index_n],
                                                          &ptw->bessel_truncation[index_n]);
            if (ptw->method[index_n] == TransferMethod::line_of_sight) {
              bessel_size                  = std::max(bessel_size, ptw->index_tau_max[index_n] + 1);
              const BesselNeeds type_needs = transfer_radial_function_needs(
                  types[index_tt].radial_type);
              needs.Phi   = needs.Phi || type_needs.Phi;
              needs.dPhi  = needs.dPhi || type_needs.dPhi;
              needs.d2Phi = needs.d2Phi || type_needs.d2Phi;
            }
          }
        }

        /* interpolate these Bessel functions once for all types and initial conditions */
        if (bessel_size > 0) {
          transfer_bessel_interpolate(ptw, index_q, index_l, tau_size, bessel_size, needs);
        }

        /* compute and store the transfer functions */
        for (int index_ic = 0; index_ic < ic_size; index_ic++) {
          for (int n = 0; n < tt_size_ts; n++) {
            const int index_tt    = sampling.types[n];
            const int index_n     = index_ic * tt_per_ts_max_ + n;
            const double* sources = ptw->sources.data() + index_n * tau_size_max_;

            double transfer_function = 0.;
            if (ptw->method[index_n] == TransferMethod::limber) {
              transfer_limber(ptw,
                              l,
                              q,
                              types[index_tt].radial_type,
                              tau0_minus_tau,
                              sources,
                              &transfer_function);
            }
            else if (ptw->method[index_n] == TransferMethod::line_of_sight) {
              transfer_integrate(ptw,
                                 index_l,
                                 k,
                                 types[index_tt].radial_type,
                                 sampling,
                                 sources,
                                 ptw->index_tau_max[index_n],
                                 ptw->bessel_truncation[index_n],
                                 tau0_minus_tau_min_bessel,
                                 &transfer_function);
            }

            if (computed(index_tt)) {
              transfer[((index_ic * tt_size_[index_md] + index_tt) * l_size_[index_md] + index_l) *
                           q_size +
                       index_q] = transfer_function;
            }
          }
        }
      } /* end of loop over l */
    } /* end of loop over time samplings */
  } /* end of loop over modes */
}

/**
 * Compute the radial coordinates (argument of the Bessel functions) of a
 * time sampling, and some trigonometric functions of them.
 *
 * @param ptw            Input/Output: workspace (results stored there)
 * @param tau0_minus_tau Input: sampled values of (tau0-tau)
 * @param tau_size       Input: number of sampled values
 * @param k              Input: wavenumber
 */

void TransferModule::transfer_radial_coordinates(struct transfer_workspace* ptw,
                                                 const double* tau0_minus_tau,
                                                 int tau_size,
                                                 double k) {
  double sqrt_absK = 0.;

  switch (ptw->sgnK) {
    case 1:
      sqrt_absK = sqrt(ptw->K);
      for (int index_tau = 0; index_tau < tau_size; index_tau++) {
        ptw->chi[index_tau]     = sqrt_absK * tau0_minus_tau[index_tau];
        ptw->cscKgen[index_tau] = sqrt_absK / k / sin(ptw->chi[index_tau]);
        ptw->cotKgen[index_tau] = ptw->cscKgen[index_tau] * cos(ptw->chi[index_tau]);
      }
      break;
    case 0:
      for (int index_tau = 0; index_tau < tau_size; index_tau++) {
        ptw->chi[index_tau]     = k * tau0_minus_tau[index_tau];
        ptw->cscKgen[index_tau] = 1.0 / ptw->chi[index_tau];
        ptw->cotKgen[index_tau] = 1.0 / ptw->chi[index_tau];
      }
      /* fully-reversed chi array, used by the flat-space Bessel interpolation */
      for (int j = 0; j < tau_size; j++) {
        ptw->chi_full_reverse[j] = ptw->chi[tau_size - 1 - j];
      }
      break;
    case -1:
      sqrt_absK = sqrt(-ptw->K);
      for (int index_tau = 0; index_tau < tau_size; index_tau++) {
        ptw->chi[index_tau]     = sqrt_absK * tau0_minus_tau[index_tau];
        ptw->cscKgen[index_tau] = sqrt_absK / k / sinh(ptw->chi[index_tau]);
        ptw->cotKgen[index_tau] = ptw->cscKgen[index_tau] * cosh(ptw->chi[index_tau]);
      }
      break;
  }
}

/**
 * This routine interpolates sources \f$ S(k, \tau) \f$ for each mode,
 * initial condition and type (of perturbation module), to get them at
 * the right values of k, using the spline interpolation method.
 *
 * @param k                     Input: wavenumber
 * @param index_md              Input: index of mode
 * @param pert_source           Input: array of sources
 * @param pert_source_spline    Input: array of second derivative of sources
 * @param interpolated_sources  Output: array of interpolated sources (filled here but allocated in transfer_init() to avoid numerous reallocation)
 * @return the error status
 */

void TransferModule::transfer_interpolate_sources(
    double k,
    int index_md,
    const double*
        pert_source, /* array with argument pert_source[index_tau*perturbations_module_->k_size_[index_md]+index_k] (must be allocated) */
    const double*
        pert_source_spline, /* array with argument pert_source_spline[index_tau*perturbations_module_->k_size_[index_md]+index_k] (must be allocated) */
    double*
        interpolated_sources /* array with argument interpolated_sources[index_tau] (must be allocated) */
) {
  /** Summary: */

  /** - define local variables */

  /* index running on k values in the original source array */
  int index_k;

  /* index running on time */
  int index_tau;

  /* variables used for spline interpolation algorithm */
  double h;
  double a;
  double b;

  /** - interpolate at each k value using the usual
      spline interpolation algorithm. */

  index_k = 0;
  h       = perturbations_module_->k_[index_md][index_k + 1] -
            perturbations_module_->k_[index_md][index_k];

  while (((index_k + 1) < perturbations_module_->k_size_[index_md]) &&
         (perturbations_module_->k_[index_md][index_k + 1] < k)) {
    index_k++;
    h = perturbations_module_->k_[index_md][index_k + 1] -
        perturbations_module_->k_[index_md][index_k];
  }

  class_test(h == 0., "stop to avoid division by zero");

  b = (k - perturbations_module_->k_[index_md][index_k]) / h;
  a = 1. - b;

  for (index_tau = 0; index_tau < perturbations_module_->tau_size_; index_tau++) {
    interpolated_sources[index_tau] =
        a * pert_source[index_tau * perturbations_module_->k_size_[index_md] + index_k] +
        b * pert_source[index_tau * perturbations_module_->k_size_[index_md] + index_k + 1] +
        ((a * a * a - a) *
             pert_source_spline[index_tau * perturbations_module_->k_size_[index_md] + index_k] +
         (b * b * b - b) * pert_source_spline[index_tau * perturbations_module_->k_size_[index_md] +
                                              index_k + 1]) *
            h * h / 6.0;
  }
}

/**
 * The code makes a distinction between "perturbation sources"
 * (e.g. gravitational potential) and "transfer sources" (e.g. total
 * density fluctuations, obtained through the Poisson equation, and
 * observed with a given selection function).
 *
 * This routine computes the transfer source of a type on its time
 * sampling, \f$ S(k,\tau) W(\tau) (k/k_{pivot})^{k_{power}} \f$, by
 * resampling the interpolated perturbation source with the precomputed
 * linear-interpolation weights of the sampling.
 *
 * @param interpolated_sources  Input: perturbation source at k, on the times of the perturbation module
 * @param k                     Input: wavenumber
 * @param index_md              Input: index of mode
 * @param index_tt              Input: index of type
 * @param sources               Output: transfer source
 */

void TransferModule::transfer_sources(
    const double* interpolated_sources, double k, int index_md, int index_tt, double* sources) {
  const TransferType& type             = types_[index_md][index_tt];
  const TransferTimeSampling& sampling = samplings_[index_md][type.index_ts];
  const int tau_size                   = static_cast<int>(sampling.tau0_minus_tau.size());
  const int* index_tau_pt              = sampling.index_tau_pt.data();
  const double* weight_tau_pt          = sampling.weight_tau_pt.data();
  const double* window                 = type.window.data();

  const double k_factor = (type.k_power != 0.) ? pow(k / type.k_pivot, type.k_power) : 1.;

  for (int index_tau = 0; index_tau < tau_size; index_tau++) {
    const int i        = index_tau_pt[index_tau];
    const double w     = weight_tau_pt[index_tau];
    sources[index_tau] = (interpolated_sources[i] * (1. - w) + w * interpolated_sources[i + 1]) *
                         window[index_tau] * k_factor;
  }
}

/**
 * Arbitrarily normalized selection function dN/dz(z,bin)
 *
 * @param bin                   Input: redshift bin number
 * @param z                     Input: one value of redshift
 * @param selection             Output: pointer to selection function
 * @return the error status
 */

void TransferModule::transfer_selection_function(int bin, double z, double* selection) {
  double x;
  double dNdz;
  double dln_dNdz_dz;
  int last_index;

  /* trivial dirac case */
  if (ppt->selection == dirac) {
    *selection = 1.;

    return;
  }

  /* tabulated case: spline interpolation of the selection function of the bin, set to zero
     outside of the tabulated range (and wherever the spline would become negative) */
  if (ppt->selection == tabulated) {
    const std::vector<double>& z_file = ptr->selection_file_z;
    if ((z < z_file.front()) || (z > z_file.back())) {
      *selection = 0.;
    }
    else {
      array_interpolate_spline(z_file.data(),
                               static_cast<int>(z_file.size()),
                               ptr->selection_file_dNdz[bin].data(),
                               selection_file_ddNdz_[bin].data(),
                               1,
                               z,
                               &last_index,
                               selection,
                               1);
      *selection = std::max(*selection, 0.);
    }

    return;
  }

  /* difference between z and the bin center (we can take the absolute
     value as long as all selection functions are symmetric around
     x=0) */
  x = fabs(z - ppt->selection_mean[bin]);

  /* gaussian case (the function is anyway normalized later
     automatically, but could not resist to normalize it already
     here) */
  if (ppt->selection == gaussian) {
    *selection = exp(-0.5 * pow(x / ppt->selection_width[bin], 2)) / ppt->selection_width[bin] /
                 sqrt(2. * _PI_);

    if ((ptr->has_nz_file) || (ptr->has_nz_analytic)) {
      if (ptr->has_nz_file) {
        class_test((z < nz_z_[0]) || (z > nz_z_[nz_size_ - 1]),
                   "Your input file for the selection function only covers the redshift range [%f "
                   ": %f]. However, your input for the selection function requires z=%f",
                   nz_z_[0],
                   nz_z_[nz_size_ - 1],
                   z);

        array_interpolate_spline(nz_z_.data(),
                                 nz_size_,
                                 nz_nz_.data(),
                                 nz_ddnz_.data(),
                                 1,
                                 z,
                                 &last_index,
                                 &dNdz,
                                 1);
      }
      else {
        transfer_dNdz_analytic(z, &dNdz, &dln_dNdz_dz);
      }

      *selection *= dNdz;
    }

    return;
  }

  /* top-hat case, with smoothed edges. The problem with sharp edges
     is that the final result will be affected by random
     noise. Indeed, the values of k at which the transfer functions
     Delta_l(k) are sampled will never coincide with the actual edges
     of the true transfer function (computed with or even without the
     Limber approximation). Hence the integral Cl=\int dk
     Delta_l(k)**2 (...) will be imprecise and will fluctuate randomly
     with the resolution along k. With smooth edges, the problem is
     solved, and the final Cls become mildly dependent on the
     resolution along k. */

  if (ppt->selection == tophat) {
    /* selection function, centered on z=mean (i.e. on x=0), equal to
       one around x=0, with tanh step centered on x=width, of width
       delta x = 0.1*width
    */
    *selection = (1. - tanh((x - ppt->selection_width[bin]) /
                            (ppr->selection_tophat_edge * ppt->selection_width[bin]))) /
                 2.;

    if ((ptr->has_nz_file) || (ptr->has_nz_analytic)) {
      if (ptr->has_nz_file) {
        array_interpolate_spline(nz_z_.data(),
                                 nz_size_,
                                 nz_nz_.data(),
                                 nz_ddnz_.data(),
                                 1,
                                 z,
                                 &last_index,
                                 &dNdz,
                                 1);
      }
      else {
        transfer_dNdz_analytic(z, &dNdz, &dln_dNdz_dz);
      }

      *selection *= dNdz;
    }

    return;
  }

  /* get here only if selection type was recognized */
  class_stop("invalid choice of selection function");
}

/**
 * Analytic form for dNdz distribution, from arXiv:1004.4640
 *
 * @param z            Input: redshift
 * @param dNdz         Output: density per redshift, dN/dZ
 * @param dln_dNdz_dz  Output: dln(dN/dz)/dz, used optionally for the source evolution
 * @return the error status
 */

void TransferModule::transfer_dNdz_analytic(double z, double* dNdz, double* dln_dNdz_dz) {
  /* Implement here your favorite analytic ansatz for the selection
     function. Typical function for photometric sample: dN/dz =
     (z/z0)^alpha exp[-(z/z0)^beta]. Then: dln(dN/dz)/dz = (alpha -
     beta*(z/z0)^beta)/z. In principle, one is free to use different
     ansatz for the selection function and the evolution
     function. Since the selection function uses only dN/dz, while the
     evolution uses only dln(dN/dz)/dz, it is possible to use
     different functions for dN/dz and dln(dN/dz)/dz */

  double z0;
  double alpha;
  double beta;

  //Euclid IST dNdz, do not change this!
  z0    = 0.9 / pow(2., 1. / 2.);
  alpha = 2.0;
  beta  = 1.5;

  *dNdz = pow(z / z0, alpha) * exp(-pow(z / z0, beta));

  *dln_dNdz_dz = (alpha - pow(z / z0, beta) * beta) / z;
}

/**
 * Compute and normalize selection function for a set of time values
 *
 * @param selection             Output: normalized selection function
 * @param tau0_minus_tau        Input: values of (tau0-tau) at which source are sample
 * @param w_trapz               Input: trapezoidal weights for integration over tau
 * @param tau_size              Input: size of previous two arrays
 * @param pvecback              Input: allocated array of background values
 * @param tau0                  Input: time today
 * @param bin                   Input: redshift bin number
 * @return the error status
 */

void TransferModule::transfer_selection_compute(double* selection,
                                                const double* tau0_minus_tau,
                                                const double* w_trapz,
                                                int tau_size,
                                                double* pvecback,
                                                double tau0,
                                                int bin) {
  /* running value of time */
  double tau;

  /* used for normalizing the selection to one */
  double norm;

  /* used for calling background_at_tau() */
  int last_index;

  /* running value of redshift */
  double z;

  if (tau_size > 1) {
    /* loop over time */
    for (int index_tau = 0; index_tau < tau_size; index_tau++) {
      /* running value of time */
      tau = tau0 - tau0_minus_tau[index_tau];

      /* get background quantities at this time */
      background_module_->background_at_tau(tau,
                                            pba->long_info,
                                            pba->inter_normal,
                                            &last_index,
                                            pvecback);

      /* infer redshift */
      z = 1. / pvecback[background_module_->index_bg_a_] - 1.;

      /* get corresponding dN/dz(z,bin) */
      transfer_selection_function(bin, z, &(selection[index_tau]));

      /* get corresponding dN/dtau = dN/dz * dz/dtau = dN/dz * H */
      selection[index_tau] *= pvecback[background_module_->index_bg_H_];
    }

    /* compute norm = \int W(tau) dtau */
    array_trapezoidal_integral(selection, tau_size, w_trapz, &norm);

    /* divide W by norm so that \int W(tau) dtau = 1 */
    for (int index_tau = 0; index_tau < tau_size; index_tau++) {
      selection[index_tau] /= norm;
    }
  }

  /* trivial case: dirac distribution */
  else {
    selection[0] = 1.;
  }
}

/**
 * Choose how to compute the transfer function of a type at a given (q, l):
 * zero (above l_max, neglected, or no overlap between the source and the
 * Bessel function), Limber approximation, or line-of-sight integral. For
 * the latter, also find the last sampled time contributing to the integral.
 *
 * @param ptw                       Input: workspace
 * @param index_q                   Input: index of wavenumber
 * @param index_md                  Input: index of mode
 * @param index_tt                  Input: index of type
 * @param index_l                   Input: index of multipole
 * @param q                         Input: wavenumber
 * @param ra_rec                    Input: angular diameter distance to recombination
 * @param q_max_bessel              Input: largest q at which the Bessel functions are tabulated
 * @param sampling                  Input: time sampling of the type
 * @param sources                   Input: transfer source of the type
 * @param use_full_limber           Input: whether we are on the full-Limber grid
 * @param tau0_minus_tau_min_bessel Input: (tau0-tau) below which the Bessel function is negligible (see transfer_bessel_cut())
 * @param index_tau_max             Output: for line_of_sight, last contributing sampled time
 * @param bessel_truncation         Output: for line_of_sight, is the integral truncated by the Bessel cut?
 * @return the method
 */

TransferMethod TransferModule::transfer_select_method(struct transfer_workspace* ptw,
                                                      int index_q,
                                                      int index_md,
                                                      int index_tt,
                                                      int index_l,
                                                      double q,
                                                      double ra_rec,
                                                      double q_max_bessel,
                                                      const TransferTimeSampling& sampling,
                                                      const double* sources,
                                                      bool use_full_limber,
                                                      double tau0_minus_tau_min_bessel,
                                                      int* index_tau_max,
                                                      char* bessel_truncation) {
  const TransferType& type                  = types_[index_md][index_tt];
  const std::vector<double>& tau0_minus_tau = sampling.tau0_minus_tau;
  const int tau_size                        = static_cast<int>(tau0_minus_tau.size());
  const double l                            = static_cast<double>(l_[index_l]);

  *index_tau_max     = -1;
  *bessel_truncation = 0;

  /** - zero above the l_max of the type, and in the full Limber scheme for types not computed there */
  if ((index_l >= l_size_tt_[index_md][index_tt]) || (use_full_limber && !type.full_limber)) {
    return TransferMethod::zero;
  }

  /** - neglect the transfer function when l is much smaller than k*tau0 */
  if (l < (q - type.delta_k_neglect) * ra_rec) {
    return TransferMethod::zero;
  }

  /** - for K>0 (closed), transfer functions only defined for l<nu */
  if ((ptw->sgnK == 1) && (l_[index_l] >= static_cast<int>(q / sqrt(ptw->K) + 0.2))) {
    return TransferMethod::zero;
  }

  /** - in the full Limber scheme, always use the Limber approximation */
  if (use_full_limber) {
    return TransferMethod::limber;
  }

  /** - in curved space, hyperspherical Bessel functions are not computed above some l */
  if ((ptw->sgnK != 0) && (index_q < index_q_flat_approximation_) && (index_l >= ptw->HIS.l_size)) {
    return TransferMethod::zero;
  }

  if (ptr->transfer_verbose > 3) {
    printf("Compute transfer for l=%d type=%d\n", static_cast<int>(l), index_tt);
  }

  /** - Limber approximation above the largest q at which the Bessel functions are tabulated,
      above the largest q at which the time sampling resolves their oscillations (otherwise
      the line-of-sight integral aliases), and above the l threshold of the type */
  //TBC: in principle the Limber condition should be adapted to account for curvature effects
  if ((q > q_max_bessel) || (q > sampling.q_max_sampled) || (l >= type.l_limber)) {
    return TransferMethod::limber;
  }

  /** - otherwise, line-of-sight integral */
  const double cut = tau0_minus_tau_min_bessel;

  /** - if there is no overlap between the region in which bessels and sources are non-zero, return zero */
  if (cut >= tau0_minus_tau[0]) {
    return TransferMethod::zero;
  }

  /** - index of the last point in the overlapping region. After this step, it can be as small
      as zero, but not negative. Remember it, to know if the truncation of the integral is due to
      the Bessel function and not to the source. */
  int index = tau_size - 1;
  while (tau0_minus_tau[index] < cut) {
    index--;
  }
  const int index_tau_max_Bessel = index;

  /** - the source function can vanish at large \f$ \tau \f$. Check if further points can be eliminated. */
  while (sources[index] == 0.) {
    index--;
    if (index < 0) {
      return TransferMethod::zero;
    }
  }

  /** - neglect late time CMB sources when l is above threshold */
  if (type.late_source_cut &&
      (l > ppr->transfer_neglect_late_source * thermodynamics_module_->angular_rescaling_)) {
    while (tau0_minus_tau[index] < ptw->tau0_minus_tau_cut) {
      index--;
      if (index < 0) {
        return TransferMethod::zero;
      }
    }
  }

  *index_tau_max     = index;
  *bessel_truncation = (index != tau_size - 1) && (index == index_tau_max_Bessel);
  return TransferMethod::line_of_sight;
}

/**
 * This routine computes the transfer functions \f$ \Delta_l^{X} (k) \f$
 * for a given type, multipole l and wavenumber k, by convolving the
 * transfer source with the radial function (Bessel function or a
 * combination of its derivatives) of the type.
 *
 * @param ptw                       Input: workspace, with the Bessel functions interpolated by transfer_bessel_interpolate()
 * @param index_l                   Input: index of multipole
 * @param k                         Input: wavenumber
 * @param radial_type               Input: type of radial function
 * @param sampling                  Input: time sampling of the transfer source
 * @param sources                   Input: transfer source
 * @param index_tau_max             Input: last sampled time contributing to the integral
 * @param bessel_truncation         Input: whether the integral is truncated by the Bessel cut
 * @param tau0_minus_tau_min_bessel Input: (tau0-tau) below which the Bessel function is negligible
 * @param trsf                      Output: transfer function \f$ \Delta_l(k) \f$
 */

void TransferModule::transfer_integrate(struct transfer_workspace* ptw,
                                        int index_l,
                                        double k,
                                        radial_function_type radial_type,
                                        const TransferTimeSampling& sampling,
                                        const double* sources,
                                        int index_tau_max,
                                        bool bessel_truncation,
                                        double tau0_minus_tau_min_bessel,
                                        double* trsf) {
  double* radial_function = ptw->radial_function.data();

  /** - Compute the radial function from the Bessel functions interpolated before: */
  transfer_radial_function(ptw, k, index_l, index_tau_max + 1, radial_function, radial_type);

  /** - Now we do most of the convolution integral: */
  array_trapezoidal_convolution(sources,
                                radial_function,
                                index_tau_max + 1,
                                sampling.w_trapz.data(),
                                trsf);

  /** - This integral is correct for the case where no truncation has
      occurred. If it has been truncated at some index_tau_max because
      f[index_tau_max+1]==0, it is still correct. The 'mistake' in using
      the wrong weight w_trapz[index_tau_max] is exactly compensated by the
      triangle we miss. However, for the Bessel cut off, we must subtract the
      wrong triangle and add the correct triangle. */
  if (bessel_truncation) {
    *trsf -= 0.5 * (sampling.tau0_minus_tau[index_tau_max + 1] - tau0_minus_tau_min_bessel) *
             radial_function[index_tau_max] * sources[index_tau_max];
  }
}

/**
 * Locate the Limber delta and the curvature factor that belongs to it.
 *
 * @param l     Input: multipole
 * @param q     Input: wavenumber
 * @param sgnK  Input: sign of the spatial curvature
 * @param K     Input: spatial curvature
 * @return the evaluation distance tau0-tau and the amplitude at that distance
 */
TransferModule::LimberPoint TransferModule::LimberDeltaPoint(double l,
                                                             double q,
                                                             int sgnK,
                                                             double K) {
  /* The Langer-corrected multipole, in every geometry: the delta of the Limber
     prescription sits at sin_K(chi_L) = L/q, and its weight carries
     cos_K(chi_L)^(-1/2) evaluated at that same L. Using the classical turning
     point sqrt(l(l+1)) in one geometry and l+1/2 in another made the lensing
     output discontinuous across Omega_k = 0 (issue #423). */
  const double L = l + 0.5;

  LimberPoint point{0., 1.};

  if (sgnK == 0) {
    point.tau0_minus_tau = L / q;
  }
  else if (sgnK == 1) {
    /* nu = q/sqrt(K) is an integer >= l+1 in a closed universe, so L/nu < 1
       and the near turning point always exists. */
    point.tau0_minus_tau = asin(L / q * sqrt(K)) / sqrt(K);
    point.amplitude      = pow(1. - K * L * L / q / q, -1. / 4.);
  }
  else if (sgnK == -1) {
    point.tau0_minus_tau = asinh(L / q * sqrt(-K)) / sqrt(-K);
    point.amplitude      = pow(1. - K * L * L / q / q, -1. / 4.);
  }

  return point;
}

/**
 * Does the source support reach past the equator of a closed universe?
 *
 * @param tau0_minus_tau_max Input: largest distance at which sources are held
 * @param sgnK               Input: sign of the spatial curvature
 * @param K                  Input: spatial curvature
 * @return whether the far turning point can fall inside the sources
 */
bool TransferModule::LimberSourceCrossesEquator(double tau0_minus_tau_max, int sgnK, double K) {
  if (sgnK != 1) {
    return false;
  }
  return tau0_minus_tau_max > 0.5 * _PI_ / sqrt(K);
}

/**
 * This routine computes the transfer functions \f$ \Delta_l^{X} (k) \f$)
 * for a given type, multipole l and wavenumber k, by using the Limber
 * approximation, i.e by evaluating the transfer source at a single value
 * of tau (the Bessel function being approximated as a Dirac distribution).
 *
 * @param ptw            Input: pointer to transfer workspace structure
 * @param l              Input: multipole
 * @param q              Input: wavenumber
 * @param radial_type    Input: type of radial (Bessel) functions to convolve with
 * @param tau0_minus_tau Input: time sampling of the transfer source
 * @param sources        Input: transfer source
 * @param trsf           Output: transfer function \f$ \Delta_l(k) \f$
 */

void TransferModule::transfer_limber(struct transfer_workspace* ptw,
                                     double l,
                                     double q,
                                     radial_function_type radial_type,
                                     const std::vector<double>& tau0_minus_tau,
                                     const double* sources,
                                     double* trsf) {
  /** Summary: */

  /** - define local variables */

  /* interpolated source and its derivatives at this value */
  double S;
  double Sp;
  double Sm;

  const int tau_size = static_cast<int>(tau0_minus_tau.size());

  double tau0_minus_tau_limber = 0.;
  double IPhiFlat              = 0.;

  if (radial_type == SCALAR_TEMPERATURE_0) {
    /** - get k, l and infer tau such that k(tau0-tau)=l+1/2;
        check that tau is in appropriate range */

    const LimberPoint limber = LimberDeltaPoint(l, q, ptw->sgnK, ptw->K);
    tau0_minus_tau_limber    = limber.tau0_minus_tau;

    if ((tau0_minus_tau_limber > tau0_minus_tau[0]) ||
        (tau0_minus_tau_limber < tau0_minus_tau[tau_size - 1])) {
      *trsf = 0.;
      return;
    }

    transfer_limber_interpolate(tau0_minus_tau.data(),
                                sources,
                                tau_size,
                                tau0_minus_tau_limber,
                                &S);

    /** - get transfer = source * \f$ \sqrt{\pi/(2l+1)}/q \f$
        = source*[tau0-tau] * \f$ \sqrt{\pi/(2l+1)}/(l+1/2)\f$
    */

    IPhiFlat = sqrt(_PI_ / (2. * l)) * (1. - 0.25 / l + 1. / 32. / (l * l));

    *trsf = IPhiFlat * S;

    if (ptw->sgnK == 0) {
      /* Identical to the curved expression below (amplitude is 1 and
         tau0_minus_tau_limber * q is l + 1/2), but spelled without the
         round trip so flat space stays bit-identical. */
      *trsf /= (l + 0.5);
    }
    else {
      *trsf *= limber.amplitude / (tau0_minus_tau_limber * q);
    }
  }

  else if (radial_type == SCALAR_TEMPERATURE_1) {
    if (((l + 1.5) / q > tau0_minus_tau[0]) || ((l - 0.5) / q < tau0_minus_tau[tau_size - 1])) {
      *trsf = 0.;
      return;
    }

    transfer_limber_interpolate(tau0_minus_tau.data(), sources, tau_size, (l + 1.5) / q, &Sp);

    transfer_limber_interpolate(tau0_minus_tau.data(), sources, tau_size, (l - 0.5) / q, &Sm);

    *trsf = -sqrt(_PI_ / (2. * l + 3.)) * Sp / (l + 1.5) * (l + 1.) / (2. * l + 1) +
            sqrt(_PI_ / (2. * l - 1.)) * Sm / (l - 0.5) * l / (2. * l + 1.);
  }

  else if (radial_type == NC_RSD) {
    if (((l + 2.5) / q > tau0_minus_tau[0]) || ((l - 1.5) / q < tau0_minus_tau[tau_size - 1])) {
      *trsf = 0.;
      return;
    }

    transfer_limber_interpolate(tau0_minus_tau.data(), sources, tau_size, (l + 2.5) / q, &Sp);

    transfer_limber_interpolate(tau0_minus_tau.data(), sources, tau_size, (l - 1.5) / q, &Sm);

    transfer_limber_interpolate(tau0_minus_tau.data(), sources, tau_size, (l + 0.5) / q, &S);

    *trsf =
        sqrt(_PI_ / (2. * l + 5.)) * Sp / (l + 2.5) * l * (l + 2.) / (2. * l + 1.) / (2. * l + 3.) -
        sqrt(_PI_ / (2. * l + 1.)) * S / (l + 0.5) * l / (2. * l + 1.) *
            (l / (2. * l - 1.) + (l + 1.) / (2. * l + 3.)) +
        sqrt(_PI_ / (2. * l - 3.)) * Sm / (l - 1.5) * l * (l - 1.) / (2. * l + 1.) / (2. * l - 1.);
  }

  else {
    class_stop("Limber approximation has not been coded for the radial_type of index %d\n",
               radial_type);
  }
}

void TransferModule::transfer_limber_interpolate(const double* tau0_minus_tau,
                                                 const double* sources,
                                                 int tau_size,
                                                 double tau0_minus_tau_limber,
                                                 double* S) {
  int index_tau;
  double dS;
  double ddS;

  /** - find  bracketing indices by bisection: the first index_tau with
      tau0_minus_tau[index_tau] <= tau0_minus_tau_limber (the times are in decreasing
      order). index_tau must be at least 1 (so that index_tau-1 is at least 0)
      and at most tau_size-2 (so that index_tau+1 is at most tau_size-1).
  */
  index_tau = static_cast<int>(std::partition_point(tau0_minus_tau + 1,
                                                    tau0_minus_tau + std::max(tau_size - 2, 1),
                                                    [tau0_minus_tau_limber](double t) {
                                                      return t > tau0_minus_tau_limber;
                                                    }) -
                               tau0_minus_tau);

  /** - interpolate by fitting a polynomial of order two; get source
      and its first two derivatives. Note that we are not
      interpolating S, but the product S*(tau0-tau). Indeed this
      product is regular in tau=tau0, while S alone diverges for
      lensing. */

  /* the case where the last of the three point is the edge (tau0=tau) must be treated separately, see below */
  if (index_tau < tau_size - 2) {
    array_interpolate_parabola(tau0_minus_tau[index_tau - 1],
                               tau0_minus_tau[index_tau],
                               tau0_minus_tau[index_tau + 1],
                               tau0_minus_tau_limber,
                               sources[index_tau - 1] * tau0_minus_tau[index_tau - 1],
                               sources[index_tau] * tau0_minus_tau[index_tau],
                               sources[index_tau + 1] * tau0_minus_tau[index_tau + 1],
                               S,
                               &dS,
                               &ddS);
  }

  /* in this case, we have stored a zero for sources[index_k*tau_size+index_tau+1]. But we can use in very good approximation the fact that S*(tau0-tau) is constant near tau=tau0 and replace sources[index_k*tau_size+index_tau+1]*tau0_minus_tau[index_tau+1] by sources[index_k*tau_size+index_tau]*tau0_minus_tau[index_tau] */
  else {
    array_interpolate_parabola(tau0_minus_tau[index_tau - 1],
                               tau0_minus_tau[index_tau],
                               tau0_minus_tau[index_tau + 1],
                               tau0_minus_tau_limber,
                               sources[index_tau - 1] * tau0_minus_tau[index_tau - 1],
                               sources[index_tau] * tau0_minus_tau[index_tau],
                               sources[index_tau] * tau0_minus_tau[index_tau],
                               S,
                               &dS,
                               &ddS);
  }
}

/**
 * Which of the Bessel function \f$ \Phi_l \f$ and its first two derivatives
 * enter the radial function of a given type.
 *
 * @param radial_type Input: type of radial function
 * @return the needed functions
 */

TransferModule::BesselNeeds TransferModule::transfer_radial_function_needs(
    radial_function_type radial_type) {
  switch (radial_type) {
    case SCALAR_TEMPERATURE_1:
      return {false, true, false};
    case SCALAR_TEMPERATURE_2:
      return {true, false, true};
    case VECTOR_TEMPERATURE_2:
    case VECTOR_POLARISATION_E:
    case TENSOR_POLARISATION_B:
      return {true, true, false};
    case TENSOR_POLARISATION_E:
      return {true, true, true};
    case NC_RSD:
      return {false, false, true};
    default:
      return {true, false, false};
  }
}

/**
 * Minimum value of \f$ (\tau_0-\tau) \f$ at which \f$ j_l(k[\tau_0-\tau]) \f$ is
 * known, given that \f$ j_l(x) \f$ is sampled above some finite value
 * \f$ x_{\min} \f$ (below which it can be approximated by zero).
 *
 * @param ptw      Input: workspace
 * @param index_q  Input: index of wavenumber
 * @param index_l  Input: index of multipole
 * @param k        Input: wavenumber
 * @return the minimum value of (tau0-tau); _HUGE_ where the hyperspherical Bessel functions are not computed
 */

double TransferModule::transfer_bessel_cut(struct transfer_workspace* ptw,
                                           int index_q,
                                           int index_l,
                                           double k) {
  const double l = static_cast<double>(l_[index_l]);

  if (ptw->sgnK == 0) {
    return ptw->pBIS->chi_at_phimin[index_l] /
           k; /* segmentation fault impossible, checked before that k != 0 */
  }

  if (index_q < index_q_flat_approximation_) {
    /* above HIS.l_size, transfer_select_method() returns zero anyway */
    if (index_l >= ptw->HIS.l_size) {
      return _HUGE_;
    }
    return ptw->HIS.chi_at_phimin[index_l] / sqrt(ptw->sgnK * ptw->K);
  }

  const double x_turning_point =
      (ptw->sgnK == 1) ? asin(sqrt(l * (l + 1.)) / q_[index_q] * sqrt(ptw->sgnK * ptw->K))
                       : asinh(sqrt(l * (l + 1.)) / q_[index_q] * sqrt(ptw->sgnK * ptw->K));
  return ptw->pBIS->chi_at_phimin[index_l] / sqrt(ptw->sgnK * ptw->K) * x_turning_point /
         sqrt(l * (l + 1.));
}

/**
 * Interpolate the Bessel functions (hyperspherical Bessel functions in
 * non-flat space) and the requested derivatives at the first x_size radial
 * coordinates of the current time sampling (stored in the workspace by
 * transfer_radial_coordinates()), once for all the types and initial
 * conditions using this sampling. They are stored in the workspace in
 * reverse order (order of increasing argument), and combined into the
 * radial function of each type by transfer_radial_function().
 *
 * @param ptw         Input/Output: workspace
 * @param index_q     Input: index of wavenumber
 * @param index_l     Input: index of multipole
 * @param tau_size    Input: number of times of the current time sampling
 * @param x_size      Input: number of times at which the functions are needed
 * @param needs       Input: which functions are needed
 */

void TransferModule::transfer_bessel_interpolate(struct transfer_workspace* ptw,
                                                 int index_q,
                                                 int index_l,
                                                 int tau_size,
                                                 int x_size,
                                                 BesselNeeds needs) {
  HyperInterpStruct* pHIS;
  const double* chi        = ptw->chi.data();
  double* chireverse       = ptw->chireverse.data();
  double* rescale_function = ptw->rescale_function.data();
  const double K           = ptw->K;
  const double l           = static_cast<double>(l_[index_l]);
  double nu                = 0.;
  double chi_tp            = 0.;
  double rescale_argument;
  double rescale_amplitude;
  using HermiteInterpolation =
      void (*)(const HyperInterpStruct*, int, int, const double*, double*, double*, double*);
  HermiteInterpolation interpolate_Phi;
  HermiteInterpolation interpolate_dPhi;
  HermiteInterpolation interpolate_PhidPhi;
  HermiteInterpolation interpolate_Phid2Phi;
  HermiteInterpolation interpolate_PhidPhid2Phi;
  enum Hermite_Interpolation_Order HIorder;

  if (ptw->sgnK == 0) {
    pHIS              = ptw->pBIS;
    rescale_argument  = 1.;
    rescale_amplitude = 1.;
    HIorder           = HERMITE4;
  }
  else if (index_q < index_q_flat_approximation_) {
    pHIS              = &(ptw->HIS);
    rescale_argument  = 1.;
    rescale_amplitude = 1.;
    HIorder           = HERMITE6;
  }
  else {
    pHIS = ptw->pBIS;
    if (ptw->sgnK == 1) {
      nu     = q_[index_q] / sqrt(K);
      chi_tp = asin(sqrt(l * (l + 1.)) / nu);
    }
    else {
      nu     = q_[index_q] / sqrt(-K);
      chi_tp = asinh(sqrt(l * (l + 1.)) / nu);
    }
    rescale_argument  = sqrt(l * (l + 1.)) / chi_tp;
    rescale_amplitude = pow(1. - K * l * (l + 1.) / q_[index_q] / q_[index_q], -1. / 12.);
    HIorder           = HERMITE4;
  }

  ptw->bessel_size      = x_size;
  ptw->rescale_argument = rescale_argument;

  if (HIorder == HERMITE4) {
    interpolate_Phi          = hyperspherical_Hermite_interpolation<4, true, false, false>;
    interpolate_dPhi         = hyperspherical_Hermite_interpolation<4, false, true, false>;
    interpolate_PhidPhi      = hyperspherical_Hermite_interpolation<4, true, true, false>;
    interpolate_Phid2Phi     = hyperspherical_Hermite_interpolation<4, true, false, true>;
    interpolate_PhidPhid2Phi = hyperspherical_Hermite_interpolation<4, true, true, true>;
  }
  else {
    interpolate_Phi          = hyperspherical_Hermite_interpolation<6, true, false, false>;
    interpolate_dPhi         = hyperspherical_Hermite_interpolation<6, false, true, false>;
    interpolate_PhidPhi      = hyperspherical_Hermite_interpolation<6, true, true, false>;
    interpolate_Phid2Phi     = hyperspherical_Hermite_interpolation<6, true, false, true>;
    interpolate_PhidPhid2Phi = hyperspherical_Hermite_interpolation<6, true, true, true>;
  }

  //Reverse chi
  if (ptw->sgnK == 0) {
    chireverse = ptw->chi_full_reverse.data() + (tau_size - x_size);
    for (int j = 0; j < x_size; j++) {
      rescale_function[j] = 1.;
    }
  }
  else {
    for (int j = 0; j < x_size; j++) {
      chireverse[j] = chi[x_size - 1 - j] * rescale_argument;
      if (rescale_amplitude == 1.) {
        rescale_function[j] = 1.;
      }
      else {
        const double x = chireverse[j] / rescale_argument;
        if (ptw->sgnK == 1) {
          rescale_function[j] = std::min(rescale_amplitude *
                                             (1 + 0.34 * atan(l / nu) * (x - chi_tp) +
                                              2.00 * pow(atan(l / nu) * (x - chi_tp), 2)),
                                         x / sin(x));
        }
        else {
          rescale_function[j] = std::max(rescale_amplitude *
                                             (1 - 0.38 * atan(l / nu) * (x - chi_tp) +
                                              0.40 * pow(atan(l / nu) * (x - chi_tp), 2)),
                                         x / sinh(x));
        }
      }
    }
  }

  class_test((pHIS->x[pHIS->x_size - 1] < chireverse[x_size - 1]) && (ptw->sgnK != 1),
             "Bessels need to be interpolated at %e, outside the range in which they have been "
             "computed (<%e). Increase their x_max.",
             chireverse[x_size - 1],
             pHIS->x[pHIS->x_size - 1]);

  /* interpolate only the needed functions (Phi comes for free with d2Phi) */
  double* Phi   = ptw->Phi.data();
  double* dPhi  = ptw->dPhi.data();
  double* d2Phi = ptw->d2Phi.data();
  if (needs.d2Phi) {
    if (needs.dPhi) {
      interpolate_PhidPhid2Phi(pHIS, x_size, index_l, chireverse, Phi, dPhi, d2Phi);
    }
    else {
      interpolate_Phid2Phi(pHIS, x_size, index_l, chireverse, Phi, nullptr, d2Phi);
    }
  }
  else if (needs.Phi) {
    if (needs.dPhi) {
      interpolate_PhidPhi(pHIS, x_size, index_l, chireverse, Phi, dPhi, nullptr);
    }
    else {
      interpolate_Phi(pHIS, x_size, index_l, chireverse, Phi, nullptr, nullptr);
    }
  }
  else {
    interpolate_dPhi(pHIS, x_size, index_l, chireverse, nullptr, dPhi, nullptr);
  }
}

/**
 * Compute the radial function (Bessel function or combination of its
 * derivatives) of a given type at the first x_size times of the current
 * time sampling, from the Bessel functions interpolated by
 * transfer_bessel_interpolate() for at least x_size times.
 *
 * @param ptw             Input: workspace (radial coordinates and Bessel functions)
 * @param k               Input: wavenumber
 * @param index_l         Input: index of multipole
 * @param x_size          Input: number of values
 * @param radial_function Output: radial function, in the order of the time sampling
 * @param radial_type     Input: type of radial function
 */

void TransferModule::transfer_radial_function(struct transfer_workspace* ptw,
                                              double k,
                                              int index_l,
                                              int x_size,
                                              double* radial_function,
                                              radial_function_type radial_type) {
  const double* cscKgen = ptw->cscKgen.data();
  const double* cotKgen = ptw->cotKgen.data();
  /* the Bessel functions are stored in reverse order: the value at index i of the time sampling is at n-i */
  const double* Phi              = ptw->Phi.data();
  const double* dPhi             = ptw->dPhi.data();
  const double* d2Phi            = ptw->d2Phi.data();
  const double* rescale_function = ptw->rescale_function.data();
  const int n                    = ptw->bessel_size - 1;
  const double rescale_argument  = ptw->rescale_argument;
  const double K                 = ptw->K;
  const double k2                = k * k;
  const double l                 = static_cast<double>(l_[index_l]);
  double factor;
  double s0;
  double s2;
  double ssqrt3;
  double si;
  double ssqrt2;
  double ssqrt2i;

  class_test(x_size > ptw->bessel_size,
             "radial function needed at %d times, but Bessel functions interpolated at only %d "
             "times",
             x_size,
             ptw->bessel_size);

  /* sqrt(|K|)/k; in flat space, the choice consistent with chi=k*(tau0-tau) and nu=1 */
  const double sqrt_absK_over_k = (ptw->sgnK == 0) ? 1.0 : sqrt(ptw->sgnK * K) / k;
  const double absK_over_k2     = sqrt_absK_over_k * sqrt_absK_over_k;

  switch (radial_type) {
    case SCALAR_TEMPERATURE_0:
      for (int i = 0; i < x_size; i++) {
        radial_function[i] = Phi[n - i] * rescale_function[n - i];
      }
      break;
    case SCALAR_TEMPERATURE_1:
      for (int i = 0; i < x_size; i++) {
        radial_function[i] = sqrt_absK_over_k * dPhi[n - i] * rescale_argument *
                             rescale_function[n - i];
      }
      break;
    case SCALAR_TEMPERATURE_2:
      s2     = sqrt(1.0 - 3.0 * K / k2);
      factor = 1.0 / (2.0 * s2);
      for (int i = 0; i < x_size; i++) {
        radial_function[i] =
            factor *
            (3 * absK_over_k2 * d2Phi[n - i] * rescale_argument * rescale_argument + Phi[n - i]) *
            rescale_function[n - i];
      }
      break;
    case SCALAR_POLARISATION_E:
      s2     = sqrt(1.0 - 3.0 * K / k2);
      factor = sqrt(3.0 / 8.0 * (l + 2.0) * (l + 1.0) * l * (l - 1.0)) / s2;
      for (int i = 0; i < x_size; i++) {
        radial_function[i] = factor * cscKgen[i] * cscKgen[i] * Phi[n - i] *
                             rescale_function[n - i];
      }
      break;
    case VECTOR_TEMPERATURE_1:
      s0     = sqrt(1.0 + K / k2);
      factor = sqrt(0.5 * l * (l + 1)) / s0;
      for (int i = 0; i < x_size; i++) {
        radial_function[i] = factor * cscKgen[i] * Phi[n - i] * rescale_function[n - i];
      }
      break;
    case VECTOR_TEMPERATURE_2:
      s0     = sqrt(1.0 + K / k2);
      ssqrt3 = sqrt(1.0 - 2.0 * K / k2);
      factor = sqrt(1.5 * l * (l + 1)) / s0 / ssqrt3;
      for (int i = 0; i < x_size; i++) {
        radial_function[i] = factor * cscKgen[i] *
                             (sqrt_absK_over_k * dPhi[n - i] * rescale_argument -
                              cotKgen[i] * Phi[n - i]) *
                             rescale_function[n - i];
      }
      break;
    case VECTOR_POLARISATION_E:
      s0     = sqrt(1.0 + K / k2);
      ssqrt3 = sqrt(1.0 - 2.0 * K / k2);
      factor = 0.5 * sqrt((l - 1.0) * (l + 2.0)) / s0 / ssqrt3;
      for (int i = 0; i < x_size; i++) {
        radial_function[i] = factor * cscKgen[i] *
                             (cotKgen[i] * Phi[n - i] +
                              sqrt_absK_over_k * dPhi[n - i] * rescale_argument) *
                             rescale_function[n - i];
      }
      break;
    case VECTOR_POLARISATION_B:
      s0     = sqrt(1.0 + K / k2);
      ssqrt3 = sqrt(1.0 - 2.0 * K / k2);
      si     = sqrt(1.0 + 2.0 * K / k2);
      factor = 0.5 * sqrt((l - 1.0) * (l + 2.0)) * si / s0 / ssqrt3;
      for (int i = 0; i < x_size; i++) {
        radial_function[i] = factor * cscKgen[i] * Phi[n - i] * rescale_function[n - i];
      }
      break;
    case TENSOR_TEMPERATURE_2:
      ssqrt2 = sqrt(1.0 - 1.0 * K / k2);
      si     = sqrt(1.0 + 2.0 * K / k2);
      factor = sqrt(3.0 / 8.0 * (l + 2.0) * (l + 1.0) * l * (l - 1.0)) / si / ssqrt2;
      for (int i = 0; i < x_size; i++) {
        radial_function[i] = factor * cscKgen[i] * cscKgen[i] * Phi[n - i] *
                             rescale_function[n - i];
      }
      break;
    case TENSOR_POLARISATION_E:
      ssqrt2 = sqrt(1.0 - 1.0 * K / k2);
      si     = sqrt(1.0 + 2.0 * K / k2);
      factor = 0.25 / si / ssqrt2;
      for (int i = 0; i < x_size; i++) {
        radial_function[i] = factor *
                             (absK_over_k2 * d2Phi[n - i] * rescale_argument * rescale_argument +
                              4.0 * cotKgen[i] * sqrt_absK_over_k * dPhi[n - i] * rescale_argument -
                              (1.0 + 4 * K / k2 - 2.0 * cotKgen[i] * cotKgen[i]) * Phi[n - i]) *
                             rescale_function[n - i];
      }
      break;
    case TENSOR_POLARISATION_B:
      ssqrt2i = sqrt(1.0 + 3.0 * K / k2);
      ssqrt2  = sqrt(1.0 - 1.0 * K / k2);
      si      = sqrt(1.0 + 2.0 * K / k2);
      factor  = 0.5 * ssqrt2i / ssqrt2 / si;
      for (int i = 0; i < x_size; i++) {
        radial_function[i] = factor *
                             (sqrt_absK_over_k * dPhi[n - i] * rescale_argument +
                              2.0 * cotKgen[i] * Phi[n - i]) *
                             rescale_function[n - i];
      }
      break;
    case NC_RSD:
      for (int i = 0; i < x_size; i++) {
        radial_function[i] = absK_over_k2 * d2Phi[n - i] * rescale_argument * rescale_argument *
                             rescale_function[n - i];
      }
      // Note: in previous line there was a missing factor absK_over_k2 until version 2.4.3. Credits Francesco Montanari.
      break;
  }
}

/* for reading the global selection function (ie the one multiplying the selection function of each
   bin), the evolution function, and splining the selection functions of selection=tabulated */

void TransferModule::transfer_global_selection_read() {
  /* for selection=tabulated, spline the selection function of each bin (read in the input module) */
  if (ppt->selection == tabulated) {
    selection_file_ddNdz_.assign(ptr->selection_file_dNdz.size(), {});
    for (size_t bin = 0; bin < ptr->selection_file_dNdz.size(); bin++) {
      selection_file_ddNdz_[bin].resize(ptr->selection_file_z.size());
      array_spline_table_lines(ptr->selection_file_z.data(),
                               static_cast<int>(ptr->selection_file_z.size()),
                               ptr->selection_file_dNdz[bin].data(),
                               1,
                               selection_file_ddNdz_[bin].data(),
                               _SPLINE_EST_DERIV_);
    }
  }

  /* for reading selection function */
  FILE* input_file;
  int row;
  int status;
  double tmp1;
  double tmp2;

  nz_size_ = 0;

  if (ptr->has_nz_file) {
    input_file = fopen(ptr->nz_file_name.c_str(), "r");
    class_test(input_file == nullptr, "Could not open file %s!", ptr->nz_file_name.c_str());

    /* Find size of table */
    for (row = 0, status = 2; status == 2; row++) {
      status = fscanf(input_file, "%lf %lf", &tmp1, &tmp2);
    }
    rewind(input_file);
    nz_size_ = row - 1;

    /* Allocate room for interpolation table */
    nz_z_.resize(nz_size_);
    nz_nz_.resize(nz_size_);
    nz_ddnz_.resize(nz_size_);

    for (row = 0; row < nz_size_; row++) {
      status = fscanf(input_file, "%lf %lf", &nz_z_[row], &nz_nz_[row]);
    }
    fclose(input_file);

    /* Call spline interpolation: */
    array_spline_table_lines(nz_z_.data(),
                             nz_size_,
                             nz_nz_.data(),
                             1,
                             nz_ddnz_.data(),
                             _SPLINE_EST_DERIV_);
  }

  nz_evo_size_ = 0;

  if (ptr->has_nz_evo_file) {
    input_file = fopen(ptr->nz_evo_file_name.c_str(), "r");
    class_test(input_file == nullptr, "Could not open file %s!", ptr->nz_evo_file_name.c_str());

    /* Find size of table */
    for (row = 0, status = 2; status == 2; row++) {
      status = fscanf(input_file, "%lf %lf", &tmp1, &tmp2);
    }
    rewind(input_file);
    nz_evo_size_ = row - 1;

    /* Allocate room for interpolation table */
    nz_evo_z_.resize(nz_evo_size_);
    nz_evo_nz_.resize(nz_evo_size_);
    nz_evo_dlog_nz_.resize(nz_evo_size_);
    nz_evo_dd_dlog_nz_.resize(nz_evo_size_);

    for (row = 0; row < nz_evo_size_; row++) {
      status = fscanf(input_file, "%lf %lf", &nz_evo_z_[row], &nz_evo_nz_[row]);
    }
    fclose(input_file);

    /* infer dlog(dN/dz)/dz from dN/dz */
    nz_evo_dlog_nz_[0] = (log(nz_evo_nz_[1]) - log(nz_evo_nz_[0])) / (nz_evo_z_[1] - nz_evo_z_[0]);
    for (row = 1; row < nz_evo_size_ - 1; row++) {
      nz_evo_dlog_nz_[row] = (log(nz_evo_nz_[row + 1]) - log(nz_evo_nz_[row - 1])) /
                             (nz_evo_z_[row + 1] - nz_evo_z_[row - 1]);
    }
    nz_evo_dlog_nz_[nz_evo_size_ - 1] = (log(nz_evo_nz_[nz_evo_size_ - 1]) -
                                         log(nz_evo_nz_[nz_evo_size_ - 2])) /
                                        (nz_evo_z_[nz_evo_size_ - 1] - nz_evo_z_[nz_evo_size_ - 2]);

    /* Call spline interpolation: */
    array_spline_table_lines(nz_evo_z_.data(),
                             nz_evo_size_,
                             nz_evo_dlog_nz_.data(),
                             1,
                             nz_evo_dd_dlog_nz_.data(),
                             _SPLINE_EST_DERIV_);
  }
};

void TransferModule::transfer_workspace_init(struct transfer_workspace* ptw,
                                             double K,
                                             int sgnK,
                                             double tau0_minus_tau_cut,
                                             HyperInterpStruct* pBIS) {
  ptw->pBIS               = pBIS;
  ptw->K                  = K;
  ptw->sgnK               = sgnK;
  ptw->tau0_minus_tau_cut = tau0_minus_tau_cut;

  int ic_tp_max = 0;
  int ic_max    = 0;
  for (int index_md = 0; index_md < md_size_; index_md++) {
    ic_tp_max = std::max(ic_tp_max,
                         perturbations_module_->ic_size_[index_md] *
                             perturbations_module_->tp_size_[index_md]);
    ic_max    = std::max(ic_max, perturbations_module_->ic_size_[index_md]);
  }
  const int n_max = ic_max * tt_per_ts_max_;

  ptw->interpolated_sources.resize(ic_tp_max);
  ptw->is_interpolated.resize(ic_tp_max);
  ptw->sources.resize(n_max * tau_size_max_);
  ptw->method.resize(n_max);
  ptw->index_tau_max.resize(n_max);
  ptw->bessel_truncation.resize(n_max);
  ptw->chi.resize(tau_size_max_);
  ptw->cscKgen.resize(tau_size_max_);
  ptw->cotKgen.resize(tau_size_max_);
  ptw->Phi.resize(tau_size_max_);
  ptw->dPhi.resize(tau_size_max_);
  ptw->d2Phi.resize(tau_size_max_);
  ptw->chireverse.resize(tau_size_max_);
  ptw->rescale_function.resize(tau_size_max_);
  ptw->radial_function.resize(tau_size_max_);
  ptw->chi_full_reverse.resize(tau_size_max_);
}

void TransferModule::transfer_update_HIS(struct transfer_workspace* ptw, int index_q, double tau0) {
  double nu;
  double new_nu;
  int int_nu;
  double xmin;
  double xmax;
  double sampling;
  double phiminabs;
  double xtol;
  double sqrt_absK;
  int l_size_max;
  int index_l_left;
  int index_l_right;

  if ((ptw->sgnK != 0) && (index_q < index_q_flat_approximation_)) {
    xmin = ppr->hyper_x_min;

    sqrt_absK = sqrt(ptw->sgnK * ptw->K);

    xmax = sqrt_absK * tau0;
    nu   = q_[index_q] / sqrt_absK;

    if (ptw->sgnK == 1) {
      xmax = std::min(xmax, _PI_ / 2.0 - ppr->hyper_x_min);  //We only need solution on [0;pi/2]

      int_nu = static_cast<int>(nu + 0.2);
      new_nu = static_cast<double>(int_nu);
      class_test(nu - new_nu > 1.e-6,
                 "problem in q list definition in closed case for index_q=%d, nu=%e, nu-int(nu)=%e",
                 index_q,
                 nu,
                 nu - new_nu);
      nu = new_nu;
    }

    if (nu > ppr->hyper_nu_sampling_step) {
      sampling = ppr->hyper_sampling_curved_high_nu;
    }
    else {
      sampling = ppr->hyper_sampling_curved_low_nu;
    }

    /* find the highest value of l such that x_nonzero < xmax = sqrt(|K|) tau0. That will be l_max. */
    l_size_max = l_size_max_;
    if (ptw->sgnK == 1) {
      while (static_cast<double>(l_[l_size_max - 1]) >= nu) {
        l_size_max--;
      }
    }

    if (ptw->sgnK == -1) {
      xtol      = ppr->hyper_x_tol;
      phiminabs = ppr->hyper_phi_min_abs;

      /* First try to find lmax using fast approximation: */
      index_l_left  = 0;
      index_l_right = l_size_max - 1;
      transfer_get_lmax(hyperspherical_get_xmin_from_approx,
                        ptw->sgnK,
                        nu,
                        l_.data(),
                        l_size_max,
                        phiminabs,
                        xmax,
                        xtol,
                        &index_l_left,
                        &index_l_right);

      /* Now use WKB approximation to eventually modify borders: */
      transfer_get_lmax(hyperspherical_get_xmin_from_Airy,
                        ptw->sgnK,
                        nu,
                        l_.data(),
                        l_size_max,
                        phiminabs,
                        xmax,
                        xtol,
                        &index_l_left,
                        &index_l_right);
      l_size_max = index_l_right + 1;
    }

    class_test(nu <= 0.,
               "nu=%e when index_q=%d, q=%e, K=%e, sqrt(|K|)=%e; instead nu should always be "
               "strictly positive",
               nu,
               index_q,
               q_[index_q],
               ptw->K,
               sqrt_absK);

    ptw->HIS = HyperInterpStruct(ptw->sgnK,
                                 nu,
                                 l_size_max,
                                 l_.data(),
                                 xmin,
                                 xmax,
                                 sampling,
                                 l_[l_size_max - 1] + 1,
                                 ppr->hyper_phi_min_abs);
  }
}

void TransferModule::transfer_get_lmax(
    void (*get_xmin_generic)(
        int sgnK, int l, double nu, double xtol, double phiminabs, double* x_nonzero, int* fevals),
    int sgnK,
    double nu,
    int* lvec,
    int lsize,
    double phiminabs,
    double xmax,
    double xtol,
    int* index_l_left,
    int* index_l_right) {
  double x_nonzero;
  int fevals = 0;
  int index_l_mid;
  int multiplier;
  int right_boundary_checked = false;
  get_xmin_generic(sgnK, lvec[0], nu, xtol, phiminabs, &x_nonzero, &fevals);
  if (x_nonzero >= xmax) {
    //x at left boundary is already larger than xmax.
    *index_l_right = std::max(lsize - 1, 1);
    return;
  }
  get_xmin_generic(sgnK, lvec[lsize - 1], nu, xtol, phiminabs, &x_nonzero, &fevals);

  if (x_nonzero < xmax) {
    //All Bessels are relevant
    *index_l_left = std::max(0, (lsize - 2));
    return;
  }
  /* Hunt for left boundary: */
  for (multiplier = 1;; multiplier *= 5) {
    get_xmin_generic(sgnK, lvec[*index_l_left], nu, xtol, phiminabs, &x_nonzero, &fevals);
    if (x_nonzero <= xmax) {
      //Boundary found
      break;
    }
    else {
      //We can use current index_l_left as index_l_right:
      *index_l_right         = *index_l_left;
      right_boundary_checked = true;
    }
    //Update index_l_left:
    *index_l_left = (*index_l_left) - multiplier;
    if (*index_l_left <= 0) {
      *index_l_left = 0;
      break;
    }
  }
  /* If not found, hunt for right boundary: */
  if (!right_boundary_checked) {
    for (multiplier = 1;; multiplier *= 5) {
      get_xmin_generic(sgnK, lvec[*index_l_right], nu, xtol, phiminabs, &x_nonzero, &fevals);
      if (x_nonzero >= xmax) {
        //Boundary found
        break;
      }
      else {
        //We can use current index_l_right as index_l_left:
        *index_l_left = *index_l_right;
      }
      //Update index_l_right:
      *index_l_right = (*index_l_right) + multiplier;
      if (*index_l_right >= (lsize - 1)) {
        *index_l_right = lsize - 1;
        break;
      }
    }
  }
  /* binary search */
  while (((*index_l_right) - (*index_l_left)) > 1) {
    index_l_mid = static_cast<int>(0.5 * ((*index_l_right) + (*index_l_left)));
    get_xmin_generic(sgnK, lvec[index_l_mid], nu, xtol, phiminabs, &x_nonzero, &fevals);
    if (x_nonzero < xmax) {
      *index_l_left = index_l_mid;
    }
    else {
      *index_l_right = index_l_mid;
    }
  }
}

/**
 * Compute once for all the window function W(tau) of each type on its
 * time sampling:
 *
 * - 1 for CMB temperature and polarisation;
 * - the lensing kernel for the CMB lensing potential;
 * - for number count contributions local in a bin, the selection
 *   function combined with the corresponding prefactor (e.g. 1/aH);
 * - for galaxy lensing and number count contributions integrated along
 *   the line of sight (labeled lens, g4, g5 in the original CLASSgal
 *   paper), the selection function convolved with the corresponding
 *   kernel.
 *
 * Factors of k are added later, through (k/k_pivot)^k_power.
 *
 * @param tau_rec               Input: recombination time
 */

void TransferModule::transfer_precompute_windows(double tau_rec) {
  /* conformal time today */
  const double tau0 = background_module_->conformal_age_;

  /* for calling background_at_tau */
  int last_index;
  std::vector<double> pvecback(background_module_->bg_size_);

  /* geometrical quantities */
  double sinKgen_source         = 0.;
  double sinKgen_source_to_lens = 0.;
  double cotKgen_source         = 0.;
  double cscKgen_lens           = 0.;

  /* source evolution factor */
  double f_evo = 0.;

  for (int index_md = 0; index_md < md_size_; index_md++) {
    for (int index_tt = 0; index_tt < tt_size_[index_md]; index_tt++) {
      TransferType& type                        = types_[index_md][index_tt];
      const TransferTimeSampling& sampling      = samplings_[index_md][type.index_ts];
      const std::vector<double>& tau0_minus_tau = sampling.tau0_minus_tau;
      const int tau_size                        = static_cast<int>(tau0_minus_tau.size());
      const int bin                             = type.bin;
      const double mag_bias       = (bin >= 0) ? ptr->selection_magnification_bias[bin] : 0.;
      std::vector<double>& window = type.window;

      window.assign(tau_size, 1.);

      switch (type.window_type) {
        /* CMB temperature and polarisation: no window */
        case TransferWindow::none:
          break;

        /* lensing source =  - W(tau) (phi(k,tau) + psi(k,tau)) Heaviside(tau-tau_rec)
           with
           psi,phi = metric perturbation in newtonian gauge (phi+psi = Phi_A-Phi_H of Bardeen)
           W = (tau-tau_rec)/(tau_0-tau)/(tau_0-tau_rec)
           H(x) = Heaviside
           (in tau = tau_0, set source = 0 to avoid division by zero;
           regulated anyway by Bessel). Possibly rescaled by hand by lcmb_rescale. */
        case TransferWindow::cmb_lensing:
          for (int index_tau = 0; index_tau < tau_size; index_tau++) {
            const double tau = tau0 - tau0_minus_tau[index_tau];
            double rescaling = 0.;
            if (index_tau < tau_size - 1) {
              switch (pba->sgnK) {
                case 1:
                  rescaling = sqrt(pba->K) * sin((tau_rec - tau) * sqrt(pba->K)) /
                              sin((tau0 - tau) * sqrt(pba->K)) /
                              sin((tau0 - tau_rec) * sqrt(pba->K));
                  break;
                case 0:
                  rescaling = (tau_rec - tau) / (tau0 - tau) / (tau0 - tau_rec);
                  break;
                case -1:
                  rescaling = sqrt(-pba->K) * sinh((tau_rec - tau) * sqrt(-pba->K)) /
                              sinh((tau0 - tau) * sqrt(-pba->K)) /
                              sinh((tau0 - tau_rec) * sqrt(-pba->K));
                  break;
              }
              // Note: until 2.4.3 there was a bug here: the curvature effects had been omitted.
            }
            window[index_tau] = rescaling * ptr->lcmb_rescale;
          }
          break;

        /* number count contributions local in the bin */
        case TransferWindow::selection: {
          /* values of selection function at sampled values of tau */
          std::vector<double> selection(tau_size);
          transfer_selection_compute(selection.data(),
                                     tau0_minus_tau.data(),
                                     sampling.w_trapz.data(),
                                     tau_size,
                                     pvecback.data(),
                                     tau0,
                                     bin);

          for (int index_tau = 0; index_tau < tau_size; index_tau++) {
            /* conformal time */
            const double tau = tau0 - tau0_minus_tau[index_tau];

            /* geometrical quantity */
            switch (pba->sgnK) {
              case 1:
                cotKgen_source = sqrt(pba->K) * cos(tau0_minus_tau[index_tau] * sqrt(pba->K)) /
                                 sin(tau0_minus_tau[index_tau] * sqrt(pba->K));
                break;
              case 0:
                cotKgen_source = 1. / (tau0_minus_tau[index_tau]);
                break;
              case -1:
                cotKgen_source = sqrt(-pba->K) * cosh(tau0_minus_tau[index_tau] * sqrt(-pba->K)) /
                                 sinh(tau0_minus_tau[index_tau] * sqrt(-pba->K));
                break;
            }

            /* corresponding background quantities */
            background_module_->background_at_tau(tau,
                                                  pba->long_info,
                                                  pba->inter_normal,
                                                  &last_index,
                                                  pvecback.data());
            const double a       = pvecback[background_module_->index_bg_a_];
            const double H       = pvecback[background_module_->index_bg_H_];
            const double H_prime = pvecback[background_module_->index_bg_H_prime_];

            /* Source evolution, used by nCl doppler and nCl gravity terms */
            if ((type.term == TransferTerm::nc_d0) || (type.term == TransferTerm::nc_d1) ||
                (type.term == TransferTerm::nc_g2)) {
              transfer_f_evo(pvecback.data(), last_index, cotKgen_source, &f_evo);
              /* Error in old CLASS 2.6.3 : Number count evolution did not respect curvature */
            }

            double rescaling = 0.;
            switch (type.term) {
              /* matter density source =  [- (dz/dtau) W(z)] * delta_m(k,tau)
                 = W(tau) delta_m(k,tau)
                 with
                 delta_m = total matter perturbation (defined in gauge-independent way, see arXiv 1307.1459)
                 W(z) = redshift space selection function = dN/dz
                 W(tau) = same wrt conformal time = dN/dtau
                 (in tau = tau_0, set source = 0 to avoid division by zero;
                 regulated anyway by Bessel). */
              case TransferTerm::nc_density:
                rescaling = ptr->selection_bias[bin] * selection[index_tau];
                break;
              /* redshift space distortion source = - [- (dz/dtau) W(z)] * (k/H) * theta(k,tau) */
              case TransferTerm::nc_rsd:
                rescaling = selection[index_tau] / H / a;
                break;
              case TransferTerm::nc_d0:
                rescaling = (f_evo - 3.) * selection[index_tau] * H * a;
                break;
              case TransferTerm::nc_d1:
                rescaling = selection[index_tau] *
                            (1. + H_prime / a / H / H +
                             (2. - 5. * mag_bias)
                                 // /tau0_minus_tau[index_tau] // in flat space
                                 * cotKgen_source  // in general case
                                 / a / H +
                             5. * mag_bias - f_evo);
                break;
              case TransferTerm::nc_g1:
                rescaling = selection[index_tau];
                break;
              case TransferTerm::nc_g2:
                rescaling = -selection[index_tau] *
                            (3. + H_prime / a / H / H +
                             (2. - 5. * mag_bias)
                                 // /tau0_minus_tau[index_tau]  // in flat space
                                 * cotKgen_source  // in general case
                                 / a / H -
                             f_evo);
                break;
              case TransferTerm::nc_g3:
                rescaling = selection[index_tau] / a / H;
                break;
              default:
                break;
            }

            window[index_tau] = rescaling;
          }
          break;
        }

        /* galaxy lensing and number count contributions integrated along the line of sight */
        case TransferWindow::integrated: {
          /* time sampling of the source selection function: evenly spaced over the bin
             (just its centre for a Dirac selection function) */
          const int tau_sources_size = (ppt->selection == dirac)
                                           ? 1
                                           : static_cast<int>(ppr->selection_sampling);
          std::vector<double> tau0_minus_tau_lensing_sources(tau_sources_size);
          std::vector<double> w_trapz_lensing_sources(tau_sources_size);
          std::vector<double> selection(tau_sources_size);

          if (tau_sources_size == 1) {
            tau0_minus_tau_lensing_sources[0] = tau0 - selection_tau_mean_[bin];
          }
          else {
            for (int index_tau = 0; index_tau < tau_sources_size - 1; index_tau++) {
              tau0_minus_tau_lensing_sources[index_tau] =
                  tau0 - selection_tau_min_[bin] -
                  (static_cast<double>(index_tau)) / (static_cast<double>(tau_sources_size) - 1.) *
                      (selection_tau_max_[bin] - selection_tau_min_[bin]);
            }
            tau0_minus_tau_lensing_sources[tau_sources_size - 1] = tau0 - selection_tau_max_[bin];
          }

          /* Compute trapezoidal weights for integration over tau */
          array_trapezoidal_mweights(tau0_minus_tau_lensing_sources.data(),
                                     tau_sources_size,
                                     w_trapz_lensing_sources.data());

          /* compute values of selection function at sampled values of tau */
          transfer_selection_compute(selection.data(),
                                     tau0_minus_tau_lensing_sources.data(),
                                     w_trapz_lensing_sources.data(),
                                     tau_sources_size,
                                     pvecback.data(),
                                     tau0,
                                     bin);

          /* loop over time and rescale */
          for (int index_tau = 0; index_tau < tau_size; index_tau++) {
            /* lensing source =  - W(tau) (phi(k,tau) + psi(k,tau)) Heaviside(tau-tau_rec)
               with
               psi,phi = metric perturbation in newtonian gauge (phi+psi = Phi_A-Phi_H of Bardeen)
               W = (tau-tau_rec)/(tau_0-tau)/(tau_0-tau_rec)
               H(x) = Heaviside
               (in tau = tau_0, set source = 0 to avoid division by zero;
               regulated anyway by Bessel).
            */
            double rescaling = 0.;

            if (index_tau < tau_size - 1) {
              for (int index_tau_sources = 0; index_tau_sources < tau_sources_size;
                   index_tau_sources++) {
                const double t0mt_source = tau0_minus_tau_lensing_sources[index_tau_sources];
                switch (pba->sgnK) {
                  case 1:
                    sinKgen_source         = sin(t0mt_source * sqrt(pba->K)) / sqrt(pba->K);
                    sinKgen_source_to_lens = sin((tau0_minus_tau[index_tau] - t0mt_source) *
                                                 sqrt(pba->K)) /
                                             sqrt(pba->K);
                    cotKgen_source         = cos(t0mt_source * sqrt(pba->K)) / sinKgen_source;
                    cscKgen_lens = sqrt(pba->K) / sin(sqrt(pba->K) * tau0_minus_tau[index_tau]);
                    break;
                  case 0:
                    sinKgen_source         = t0mt_source;
                    sinKgen_source_to_lens = (tau0_minus_tau[index_tau] - t0mt_source);
                    cotKgen_source         = 1. / (t0mt_source);
                    cscKgen_lens           = 1. / (tau0_minus_tau[index_tau]);
                    break;
                  case -1:
                    sinKgen_source         = sinh(t0mt_source * sqrt(-pba->K)) / sqrt(-pba->K);
                    sinKgen_source_to_lens = sinh((tau0_minus_tau[index_tau] - t0mt_source) *
                                                  sqrt(-pba->K)) /
                                             sqrt(-pba->K);
                    cotKgen_source         = cosh(t0mt_source * sqrt(-pba->K)) / sinKgen_source;
                    cscKgen_lens = sqrt(-pba->K) / sinh(sqrt(-pba->K) * tau0_minus_tau[index_tau]);
                    break;
                }

                /* condition for excluding from the sum the sources located in z=zero */
                if ((t0mt_source > 0.) && (t0mt_source - tau0_minus_tau[index_tau] > 0.)) {
                  const double sel = selection[index_tau_sources];
                  const double w   = w_trapz_lensing_sources[index_tau_sources];

                  switch (type.term) {
                    case TransferTerm::galaxy_lensing:
                      rescaling += sinKgen_source_to_lens * cscKgen_lens / sinKgen_source * sel * w;
                      break;
                    case TransferTerm::nc_lens:
                      rescaling -= (2. - 5. * mag_bias) / 2. * sinKgen_source_to_lens *
                                   cscKgen_lens / sinKgen_source * sel * w;
                      break;
                    case TransferTerm::nc_g4:
                      rescaling += (2. - 5. * mag_bias) * cotKgen_source * sel * w;
                      break;
                    case TransferTerm::nc_g5: {
                      /* background quantities and source evolution at time tau_lensing_source */
                      background_module_->background_at_tau(tau0 - t0mt_source,
                                                            pba->long_info,
                                                            pba->inter_normal,
                                                            &last_index,
                                                            pvecback.data());
                      transfer_f_evo(pvecback.data(), last_index, cotKgen_source, &f_evo);
                      const double a       = pvecback[background_module_->index_bg_a_];
                      const double H       = pvecback[background_module_->index_bg_H_];
                      const double H_prime = pvecback[background_module_->index_bg_H_prime_];
                      rescaling += (1. + H_prime / a / H / H +
                                    (2. - 5. * mag_bias) * cotKgen_source / a / H + 5. * mag_bias -
                                    f_evo) *
                                   sel * w;
                      break;
                    }
                    default:
                      break;
                  }
                }
              }
            }

            /* Finally store integrated result for later use */
            window[index_tau] = rescaling;
          }
          break;
        }
      }
    }
  }
}

void TransferModule::transfer_f_evo(
    double* pvecback,
    int last_index,
    double cotKgen, /* Should be FILLED with values of corresponding time */
    double* f_evo) {
  /* Allocate temporary variables for calculation of f_evo */
  double z;
  double dNdz;
  double dln_dNdz_dz;
  double temp_f_evo;

  if ((ptr->has_nz_evo_file) || (ptr->has_nz_evo_analytic)) {
    temp_f_evo = 2. / pvecback[background_module_->index_bg_H_] /
                     pvecback[background_module_->index_bg_a_] * cotKgen +
                 pvecback[background_module_->index_bg_H_prime_] /
                     pvecback[background_module_->index_bg_H_] /
                     pvecback[background_module_->index_bg_H_] /
                     pvecback[background_module_->index_bg_a_];

    z = 1. / pvecback[background_module_->index_bg_a_] - 1.;

    if (ptr->has_nz_evo_file) {
      class_test((z < nz_evo_z_[0]) || (z > nz_evo_z_[nz_evo_size_ - 1]),
                 "Your input file for the selection function only covers the redshift range [%f : "
                 "%f]. However, your input for the selection function requires z=%f",
                 nz_evo_z_[0],
                 nz_evo_z_[nz_evo_size_ - 1],
                 z);

      array_interpolate_spline(nz_evo_z_.data(),
                               nz_evo_size_,
                               nz_evo_dlog_nz_.data(),
                               nz_evo_dd_dlog_nz_.data(),
                               1,
                               z,
                               &last_index,
                               &dln_dNdz_dz,
                               1);
    }
    else {
      transfer_dNdz_analytic(z, &dNdz, &dln_dNdz_dz);
    }

    temp_f_evo -= dln_dNdz_dz / pvecback[background_module_->index_bg_a_];
  }
  else {
    temp_f_evo = 0.;
  }

  /* after obtaining f_evo, store it in output */
  *f_evo = temp_f_evo;
}
