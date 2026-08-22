/** @file hyrec_model.h HYREC-2 as a supplier of ionization derivatives. */

#pragma once

#include <string>

#include "recombination_model.h"

/* Opaque here: HYREC-2's headers are C and are included only in hyrec_model.cpp,
   which owns the extern "C" seam. Keeping them out of this header matters -- they
   define macros such as RECFAST and FULL that would leak into every translation
   unit that touches a recombination model. */
struct HyrecData;

/**
 * HYREC-2 (Ali-Haimoud, Hirata & Lee; arXiv:2007.14114), driven as a derivative
 * supplier rather than as a history builder.
 *
 * This is class_public v3's arrangement: CLASS integrates x_H, x_He and T_mat
 * with its own evolver and its own background, and calls HYREC-2 only for the
 * two ionization derivatives. Two things follow.
 *
 * First, HYREC-2 never computes an expansion rate -- `rec_HubbleRate` is not
 * reached -- so the CPL dark-energy reconstruction it would otherwise use, and
 * the massive neutrino that reconstruction double-counted, are simply out of the
 * picture (issues #396 and #369). Whatever CLASS's background says, including
 * scalar-field dark energy, is what recombination sees.
 *
 * Second, the model is SWIFT, HYREC-2's own default and a fit to its full
 * radiative-transfer calculation. FULL is unreachable this way: it needs the
 * radiation-field history indexed on HYREC-2's uniform lna grid, which an
 * adaptive evolver taking rejected steps cannot supply. HYREC-2's own guards
 * still fall back to PEEBLES where the effective-rate tables run out.
 */

class HyrecModel : public RecombinationModel {
 public:
  /**
   * @param path      directory holding HYREC-2's data tables
   * @param T_cmb     CMB temperature today in K
   * @param obh2      Omega_b h^2
   * @param ocbh2     Omega h^2 of NON-free-streaming matter -- baryons and cold
   *                  dark matter, excluding massive neutrinos. It enters SWIFT's
   *                  fiducial-difference correction, and including neutrinos here
   *                  would reintroduce #396's double-count in a different slot.
   * @param okh2      Omega_k h^2
   * @param Neff      effective number of neutrino species
   * @param YHe       primordial helium mass fraction
   * @param fHe       helium abundance by number, relative to hydrogen
   * @param nH0       hydrogen number density today in m^-3
   * @param z_initial redshift at which CLASS starts the recombination integration
   */
  HyrecModel(const std::string& path,
             double T_cmb,
             double obh2,
             double ocbh2,
             double okh2,
             double Neff,
             double YHe,
             double fHe,
             double nH0,
             double z_initial);
  ~HyrecModel() override;

  HyrecModel(const HyrecModel&)            = delete;
  HyrecModel& operator=(const HyrecModel&) = delete;

  IonisationDerivatives Derivatives(const RecombinationState& state,
                                    const EnergyDeposition& dep,
                                    double energy_rate) const override;

  const char* Name() const override;

 private:
  HyrecData* data_;
  std::string path_; /**< kept alive: HYREC-2 stores the pointer, not a copy */
  double fHe_;
};
