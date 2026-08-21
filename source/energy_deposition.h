/** @file energy_deposition.h How injected energy is split between deposition channels. */

#ifndef __ENERGY_DEPOSITION__
#define __ENERGY_DEPOSITION__

/**
 * List of possible energy-deposition functions chi_c(x).
 *
 * Energy injected by decaying or annihilating dark matter does not all end up in
 * the same place: some ionizes hydrogen or helium, some excites the Lyman-alpha
 * transition, some heats the gas, and some is left below the Lyman-alpha
 * threshold and does nothing. These are the fractions chi_c(x) of the *deposited*
 * energy that reach each channel, as a function of the ionization fraction.
 */

enum energy_deposition_function {
  deposition_Galli_2013, /**< tabulated fractions computed by Galli et al. 2013 (arXiv:1306.0563) */
  deposition_legacy      /**< the two analytic fits CLASS carried inline before the table was
                              available here; provided so pre-existing results stay reproducible */
};

/**
 * Fractions of the deposited energy reaching each channel, each in [0,1].
 *
 * They are NOT a partition of unity. Galli et al. 2013 computes every channel
 * separately, and the published rows sum to between 0.993 and 1.019 across the
 * table. So do not read one channel as "whatever the others did not take", and
 * do not rely on the set being exactly energy-conserving.
 */

struct EnergyDeposition {
  double heat;   /**< heating of the baryon gas */
  double ion_H;  /**< ionization of hydrogen */
  double ion_He; /**< ionization of helium */
  double lya;    /**< Lyman-alpha excitation of hydrogen */
  double lowE;   /**< photons below 10.2 eV, which deposit nothing */
};

/**
 * Evaluate the deposition fractions at ionization fraction x.
 *
 * x is CLASS's x = x_H + f_He x_He, so it exceeds 1 while helium is still
 * ionized; the fractions saturate to pure heating there. Values are clamped to
 * [0,1] so that spline ringing across the shoulder near x = 1 can never produce a
 * negative fraction.
 */

EnergyDeposition energy_deposition_fractions(energy_deposition_function which, double x);

#endif
