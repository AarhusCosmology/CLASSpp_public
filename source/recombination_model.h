/** @file recombination_model.h The atomic physics behind dx_H/dz and dx_He/dz. */

#pragma once

#include "energy_deposition.h"

/**
 * The fine-structure constant and the electron mass at one redshift, relative to
 * today, and the power laws through which they enter the thermal history (Hart &
 * Chluba, arXiv:1705.03925). {1, 1} is standard physics.
 */
struct FundamentalConstants {
  double alpha = 1.; /**< alpha(z)/alpha_0 */
  double me    = 1.; /**< m_e(z)/m_e,0 */

  /** Atomic energy levels, so every ionization temperature: alpha^2 m_e. The atomic
   *  physics applies it by dividing the temperature instead. */
  double TemperatureRescale() const {
    return alpha * alpha * me;
  }
  /** Thomson cross section: alpha^2/m_e^2. */
  double ThomsonRescale() const {
    return alpha * alpha / (me * me);
  }
  /** Compton heating of the baryons, sigma_T/m_e: alpha^2/m_e^3. */
  double ComptonRescale() const {
    return alpha * alpha / (me * me * me);
  }
  /** Saha equation, left over once T has been divided by TemperatureRescale():
   *  alpha^3 m_e^3. */
  double SahaRescale() const {
    return alpha * alpha * alpha * me * me * me;
  }
};

/**
 * State handed to a recombination model at one point of the integration.
 *
 * Everything a model needs and nothing it should be choosing for itself: the
 * expansion rate comes from CLASS's background, not from any internal
 * reconstruction, and the phase flags are owned by the driver rather than
 * re-derived from the state (mirroring the perturbations module).
 */

struct RecombinationState {
  double z;     /**< redshift */
  double x_H;   /**< HII fraction, relative to total hydrogen */
  double x_He;  /**< HeII fraction, relative to total helium */
  double x;     /**< x_H + f_He x_He, the free-electron fraction per hydrogen nucleus */
  double n_H;   /**< hydrogen number density in m^-3 */
  double H;     /**< Hubble rate in s^-1 */
  double T_mat; /**< baryon temperature in K */
  double T_rad; /**< photon temperature in K */

  bool hydrogen_frozen; /**< hydrogen follows its Saha branch; report dx_H/dz = 0 */
  bool helium_ode;      /**< helium is being evolved rather than held at Saha equilibrium */

  FundamentalConstants constants; /**< alpha and m_e at z */
};

/** Derivatives with respect to redshift of the two ionization fractions. */

struct IonisationDerivatives {
  double dx_H_dz;
  double dx_He_dz;
};

/**
 * A recombination model supplies the hydrogen and helium ionization derivatives;
 * CLASS integrates them, evolves the baryon temperature, and owns the phase
 * switching. Both derivatives come from one call because every model shares
 * expensive intermediates between them.
 */

class RecombinationModel {
 public:
  virtual ~RecombinationModel() = default;

  virtual IonisationDerivatives Derivatives(const RecombinationState& state,
                                            const EnergyDeposition& dep,
                                            double energy_rate) const = 0;

  /** Name for verbose output. */
  virtual const char* Name() const = 0;
};
