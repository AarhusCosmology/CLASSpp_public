/** @file recfast_model.h RECFAST's effective three-level atom. */

#pragma once

#include "precision.h"
#include "recombination_model.h"
#include "thermodynamics.h"

/**
 * RECFAST 1.5's hydrogen and helium ionization derivatives: Peebles' effective
 * three-level atom with the fudge factor and the Gaussian corrections of
 * version 1.5, plus the version 1.4 helium triplet correction.
 *
 * The code is the one the thermodynamics module carried inline before HYREC-2
 * arrived, moved here unchanged so that the two recombination models sit behind
 * one interface.
 */

class RecfastModel : public RecombinationModel {
 public:
  RecfastModel(const precision* ppr, const recombination* preco) : ppr_(ppr), preco_(preco) {}

  IonisationDerivatives Derivatives(const RecombinationState& state,
                                    const EnergyDeposition& dep,
                                    double energy_rate) const override;

  const char* Name() const override {
    return "RECFAST";
  }

 private:
  const precision* const ppr_;
  const recombination* const preco_;
};
