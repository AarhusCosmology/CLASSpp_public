//
//  dark_radiation.cpp
//  ppCLASS
//
//  Created by Emil Brinch Holm on 29/12/2020.
//  Copyright © 2020 Aarhus University. All rights reserved.
//

#include "dark_radiation.h"

std::shared_ptr<DarkRadiation> DarkRadiation::Create(FileContent* pfc, std::vector<SourceType> source_types, std::vector<DRType> dr_types, std::vector<double> deg, double T_cmb) {
  try {
    return std::shared_ptr<DarkRadiation>(new DarkRadiation(pfc, source_types, dr_types, deg, T_cmb));
  }
  catch (NoDrRequested& error) {
    return nullptr;
  }
}

/**
 *  The constructor reads from input files and creates quadrature abscissas.
 *
 *
 *
 *
*/

DarkRadiation::DarkRadiation(FileContent* pfc, std::vector<SourceType> source_types, std::vector<DRType> dr_types, std::vector<double> deg, double T_cmb)
: N_species_(static_cast<int>(source_types.size())), source_types_(source_types), dr_types_(dr_types), deg_(deg) {
  if (this->Init(pfc, T_cmb) == _FAILURE_) {
    throw std::runtime_error(error_message_);
  }
}

/**
 *
 *
 *
 *
 *
*/

int DarkRadiation::Init(FileContent* pfc, double T_cmb) {
  if (source_types_.size() == 0) {
    // Throw error to let DarkRadiation::Create know to continue with no DR module
    throw NoDrRequested("No DR species requested; continuing without a DR module.");
  }
  
  /* Get number of source types */
  N_dcdm_ = 0;
  N_dncdm_ = 0;
  for (SourceType source : source_types_) {
    if (source == SourceType::dcdm) {
      ++N_dcdm_;
    } else if (source == SourceType::dncdm) {
      ++N_dncdm_;
    }
  }
  
  int int1, flag1, n, entries_read;
  int* temp_list;
  char* errmsg = error_message_;
  /* Quadrature strategy */
  class_read_list_of_integers_or_default("Quadrature strategy DR", temp_list, quadrature_strategy_default_, 1);
  quadrature_strategy_ = temp_list[0];
  
  /* Maximum q */
  class_read_list_of_integers_or_default("Maximum q DR", temp_list, qmax_default_, 1);
  qmax_ = temp_list[0];
  
  /* Number of momentum bins, same for all species */
  class_read_list_of_integers_or_default("N_q_dr", temp_list, N_q_default_, 1);
  N_q_ = temp_list[0];
  
  q_.resize(N_q_);
  dq_.resize(N_q_);
  w_.resize(N_q_);
  
  /* Handle q-sampling
     Currently the same for background and perturbations
  */
  if (quadrature_strategy_ == qm_auto) {
    throw std::runtime_error("Automatic DR quadrature not yet implemented. Assign the variable Quadrature strategy DR.");
  }
  else {
    /* Manual q-sampling */
    class_call(get_qsampling_manual(q_.data(),
                                    w_.data(),
                                    dq_.data(),
                                    N_q_,
                                    qmax_,
                                    (enum quadrature_method) quadrature_strategy_,
                                    NULL,
                                    0,
                                    InitialDistribution,
                                    NULL,
                                    error_message_),
               error_message_, error_message_);
    /*
    printf("Created the following quadrature scheme for dark radiation:\n");
    for (int n = 0; n < N_q_; n++) {
      printf("q = %g, dq = %g, w = %g \n", q_[n], dq_[n], w_[n]);
    }
    */
  }
  
  // Prepare member vectors
  cumulative_q_index_.push_back(0);
  for (int n_dr = 0; n_dr < N_species_; n_dr++) {
    if (source_types_[n_dr] == SourceType::dncdm) {
      // Assume all dncdm sourced DRs have inverse decays if has_inv is true
      // Here, we make custom initial conditions that allow the fermionic daughter to have a Fermi-Dirac initial distribution
      //   with all other initial distribution being zero. Thus, we currently don't use InitialDistribution for anything.
      int has_inv = 0;
      class_read_int("Inverse decay", has_inv);
      int initial_fermionic_pop = 0;
      class_read_int("Initial daughter population", initial_fermionic_pop);

        // If no inv, we have one species
        std::vector<double> w_temp(N_q_);
          for (int index_q = 0; index_q < N_q_; index_q++) {
            if (initial_fermionic_pop == _FALSE_) {
              w_temp[index_q] = 0.;
            }
            else if (initial_fermionic_pop == _TRUE_) {
              double q = q_[index_q];
              double f0;
              InitialDistribution(NULL, q, &f0);
              w_temp[index_q] = f0*dq_[index_q];
            }
          }
        w_species_.push_back(w_temp);

      /*
      printf("Made the following initial population:\n");
      if (dr_types_[n_dr] == DRType::fermion) printf("Fermionic daughter:\n");
      if (dr_types_[n_dr] == DRType::boson) printf("Bosonic daughter:\n");
      for (int index_q = 0; index_q < N_q_; index_q++) {
        printf("w_(q=%g) = %g \n", q_[index_q], w_species_.back()[index_q]);
      }
      */

      cumulative_q_index_.push_back(cumulative_q_index_.back() + N_q_); /* All species have same q_size */

      // Only the dncdm-sourced DR carry a degeneracy factor
      double T_dr_over_T_cmb = 0.71611; // Same as T_ncdm_default
      factor_.push_back(deg_[n_dr]*4*_PI_*pow(T_dr_over_T_cmb*T_cmb*_k_B_, 4)*8*_PI_*_G_/3./pow(_h_P_/2./_PI_, 3)/pow(_c_, 7)*_Mpc_over_m_*_Mpc_over_m_);

      /*std::vector<double> w_temp(N_q_);
      for (int index_q = 0; index_q < N_q_; index_q++) {
        // Initialize all single-source weights identically
        w_temp[index_q] = w_[index_q]/(1.0*N_dncdm_);
      }
      w_species_.push_back(w_temp);*/
    }


    // Maybe needs fixing with initial populations; double check this!
    rho_species_.push_back(0.);
  }
  return _SUCCESS_;
}

/**
 *
 *
 *
 *
 *
*/

int DarkRadiation::InitialDistribution(void * params, double q, double* f0) {
  // Here, we should put a Taylor expansion around small a of the DR EoM to include the small amount of decay that has already occurred at a_ini, otherwise all weights are identically zero and quadrature does not make sense!
  double f_FD = 1.0/pow(2*_PI_, 3)*2./(exp(q) + 1.);

  // *f0 = 0.0*q;
  *f0 = f_FD;
  return _SUCCESS_;
}

/**
 *
 *
 *
 *
 *
*/

void DarkRadiation::IntegrateDistribution(double z, double* number, double* rho, double* p, int index_dr) {
  std::vector<double> w_vec(N_q_);
  
  // Which species to sum over?
  if (index_dr == 42) { // Sum over total distribution function, this is default
    for (int index_q = 0; index_q < N_q_; index_q++) {
      w_vec[index_q] = w_[index_q];
    }
  } else {
    for (int index_q = 0; index_q < N_q_; index_q++) {
      w_vec[index_q] = w_species_[index_dr][index_q];
    }
  }
  
  if (number != NULL) *number = 0.;
  if ((rho != NULL) || (p != NULL)) *rho = 0.;
  if (p != NULL) *p = 0.;
  
  /* Sum over quadrature abscissas
     the q^2 is from d^3 q volume element
  */
  for (int index_q = 0; index_q < N_q_; index_q++) {
    if (number != NULL) *number += pow(q_[index_q],2)*w_vec[index_q];
    if ((rho != NULL) || (p != NULL)) *rho += pow(q_[index_q],3)*w_vec[index_q]; // DR is massless, so epsilon = q
  }
  if (p != NULL) *p = *rho/3.; // Equation of state for a massless particle is p = rho/3
  
  /* Adjust normalization */
  // The correct factor needs double checking...
  double integral_factor = factor_[index_dr]*pow(1. + z, 4);
  if (number != NULL) *number *= integral_factor/(1. + z);
  if (rho != NULL) *rho *= integral_factor;
  if (p != NULL) *p *= integral_factor;
}

/**
 *
 *
 *
 *
 *
*/

double DarkRadiation::ComputeMeanMomentum(int index_dr) {
  std::vector<double> w_vec(N_q_);
  double zero_test = 0.;
  // Which species to sum over?
  if (index_dr == 42) { // Sum over total distribution function, this is default
    for (int index_q = 0; index_q < N_q_; index_q++) {
      w_vec[index_q] = w_[index_q];
      zero_test += w_[index_q];
    }
  } else {
    for (int index_q = 0; index_q < N_q_; index_q++) {
      w_vec[index_q] = w_species_[index_dr][index_q];
      zero_test += w_species_[index_dr][index_q];
    }
  }
  if (zero_test == 0.) {
    // No DR particles, possibly since no decay has occurred yet
    return 0.;
  }
  double num = 0.;
  double denom = 0.;
  for (int index_q = 0; index_q < N_q_; index_q++) {
    double q = q_[index_q];
    num += w_vec[index_q]*pow(q, 3);
    denom += w_vec[index_q]*pow(q, 2);
  }
  return num/denom;
}
