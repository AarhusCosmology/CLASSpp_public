//
//  dark_radiation.cpp
//  ppCLASS
//
//  Created by Emil Brinch Holm on 29/12/2020.
//  Copyright © 2020 Aarhus University. All rights reserved.
//

#include "dark_radiation.h"

std::shared_ptr<DarkRadiation> DarkRadiation::Create(FileContent* pfc, std::vector<SourceType> source_types, double T_cmb) {
  if (source_types.size() > 0) {
    try {
        return std::shared_ptr<DarkRadiation>(new DarkRadiation(pfc, source_types, T_cmb));
    }
    catch (std::exception& error) {
      printf("Could not create DarkRadiation class:\n %s", error.what());
      return nullptr;
    }
  }
  else {
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

DarkRadiation::DarkRadiation(FileContent* pfc, std::vector<SourceType> source_types, double T_cmb)
: N_species_(source_types.size()), source_types_(source_types) {
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
  // Do we want separate quadrature for each species? Or the same quadrature as the ncdm species?
  
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
    rho_species_.push_back(0.);
    
    if (source_types_[n_dr] == SourceType::dncdm) {
      cumulative_q_index_.push_back(cumulative_q_index_.back() + N_q_); /* All species have same q_size */
      std::vector<double> w_temp(N_q_);
      for (int index_q = 0; index_q < N_q_; index_q++) {
        // Initialize all single-source weights identically
        w_temp[index_q] = w_[index_q]/(1.0*N_dncdm_);
      }
      w_species_.push_back(w_temp);
    }
  }
  double T_dr_over_T_cmb = 0.71611; // Same as T_ncdm_default
  
  // Conventions copied from the NonColdDarkMatter factor
  // TO-DO: Double check!
  factor_ = 4*_PI_*pow(T_dr_over_T_cmb*T_cmb*_k_B_, 4)*8*_PI_*_G_/3./pow(_h_P_/2./_PI_, 3)/pow(_c_, 7)*_Mpc_over_m_*_Mpc_over_m_;
  
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

  *f0 = 0.0*q;
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
  double integral_factor = factor_*pow(1. + z, 4);
  if (number != NULL) *number *= integral_factor/(1. + z);
  if (rho != NULL) *rho *= integral_factor;
  if (p != NULL) *p *= integral_factor;
}
