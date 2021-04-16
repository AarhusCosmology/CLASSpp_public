//
//  dark_radiation.h
//  ppCLASS
//
//  Created by Emil Brinch Holm on 29/12/2020.
//  Copyright © 2020 Aarhus University. All rights reserved.
//

#ifndef dark_radiation_h
#define dark_radiation_h

#include <stdio.h>
#include "common.h"
#include "parser.h"
#include "quadrature.h"
#include "input_module.h"



class DarkRadiation {
public:
  enum class SourceType { dcdm, dncdm };

  static std::shared_ptr<DarkRadiation> Create(FileContent* pfc, std::vector<SourceType> source_types, double T_cmb);
  void IntegrateDistribution(double z, double* number, double* rho, double* p, int index_dr = 42);
  double ComputeMeanMomentum(double z, int index_dr = 42);

  std::vector<double> q_;
  std::vector<double> dq_;
  std::vector<double> w_; /* Total weights */
  std::vector<std::vector<double> > w_species_; /* Individual weights from separate sources */
  
  std::vector<SourceType> source_types_;
  
  std::vector<int> cumulative_q_index_;
  std::vector<double> rho_species_; /* Energy from each source channel */
  
  double factor_; /* Factor multiplying integrals over distribution functions */
  int N_species_; /* Amount of channels that source DR */
  int N_dcdm_;
  int N_dncdm_; /* Amount of dncdm sources */
  int N_q_; /* Number of momentum bins */

  mutable ErrorMsg error_message_;

private:
  DarkRadiation(FileContent* pfc, std::vector<SourceType> source_types, double T_cmb);
  int Init(FileContent* pfc, double T_cmb);
  static int InitialDistribution(void * params, double q, double* f0);
  
  const int N_q_default_ = 5;
  int quadrature_strategy_;
  const int quadrature_strategy_default_ = 3;
  double qmax_;
  const double qmax_default_ = 15.;

};



#endif /* dark_radiation_h */
