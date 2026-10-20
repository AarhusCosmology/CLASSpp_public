#include "hyrec_model.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <utility>

#include "constants.h"
#include "errors.h"

/* HYREC-2 is C, and this is the only place in CLASS that includes its headers.
   history.h pulls in hyrectools.h, hydrogen.h and energy_injection.h. */
extern "C" {
#include "helium.h"
#include "history.h"
}

/* HYREC-2's HYREC_DATA is an anonymous-struct typedef, so it cannot be forward
   declared. Wrap it. */
struct HyrecData {
  HYREC_DATA d;
};

namespace {

/* Ionization and Lyman-alpha energies of hydrogen, in eV. The Lyman-alpha
   transition is close to 3/4 of the ionization energy (1s -> 2p, E ~ 1/n^2). */
constexpr double kE_H_ion = 13.5984336478;
constexpr double kE_H_lya = 10.1988356821;

/* HYREC-2's helium rate calibrated to CosmoRec, in the form CAMB 2 uses for RECFAST (Lewis,
   arXiv:2607.14854 eq. 23) but refitted to HYREC-2: dx_HeII/dlna is multiplied by
   F = 1 + (a0 + a1 u + a2 u^2) / (1 + d2 u^2 + d4 u^4), u = (z_scale - z0) / w, where z_scale is
   the redshift at which a 2.7255 K sky has today's photon temperature (in atomic units, so it
   follows He I recombination under varying alpha and m_e). F peaks at 1.65 near z = 1780, where
   CosmoRec's 584 A photons, diffusing into the H I wings, recombine helium fastest relative to
   HYREC-2's escape probability. No window: F -> 1 away from the peak (1.005 at z = 3000), so the
   rate stays continuous. Fitted at one cosmology to direct CosmoRec (Lewis' reference settings);
   docs/superpowers/specs/2026-10-08-hyrec-helium-cosmorec-design.md has the fit and its transfer. */
double HeliumRateCorrection(double T_rad_atomic_eV) {
  constexpr double a0  = 0.635977;
  constexpr double a1  = 0.298024;
  constexpr double a2  = 0.239789;
  constexpr double d2  = 2.91828;
  constexpr double d4  = 3.31521;
  constexpr double z0  = 1751.94;
  constexpr double w   = 286.606;
  const double z_scale = T_rad_atomic_eV / (kBoltz * 2.7255) - 1.;
  const double u       = (z_scale - z0) / w;
  const double u2      = u * u;
  return 1. + (a0 + a1 * u + a2 * u2) / (1. + d2 * u2 + d4 * u2 * u2);
}

}  // namespace

HyrecModel::HyrecModel(std::string path,
                       double T_cmb,
                       double obh2,
                       double ocbh2,
                       double okh2,
                       double Neff,
                       double YHe,
                       double fHe,
                       double nH0,
                       double z_initial,
                       bool He_rate_correction)
    : data_(new HyrecData), path_(std::move(path)), fHe_(fHe),
      He_rate_correction_(He_rate_correction) {
  /* hyrec_allocate reads the effective-rate and SWIFT-fit tables from this
     directory, so it has to be set first. HYREC-2 keeps the pointer, hence the
     copy that outlives the call. */
  data_->d.path_to_hyrec = path_.data();

  /* Integrating down to z = 0: CLASS asks for derivatives across its whole
     recombination range. Only the output-array sizes depend on these bounds --
     we never ask HYREC-2 to build or interpolate a history. */
  hyrec_allocate(&data_->d, z_initial, 0.);
  class_test(data_->d.error != 0, "HYREC-2 failed to initialise: %s", data_->d.error_message);

  REC_COSMOPARAMS* cosmo = data_->d.cosmo;
  cosmo->T0              = T_cmb;
  cosmo->obh2            = obh2;
  cosmo->ocbh2           = ocbh2;
  cosmo->okh2            = okh2;
  cosmo->Neff            = Neff;
  cosmo->YHe             = YHe;
  cosmo->fHe             = fHe;
  cosmo->nH0             = nH0 * 1e-6; /* HYREC-2 wants cm^-3 */

  /* Varying fundamental constants: set per call in Derivatives(). */
  cosmo->fsR = 1.;
  cosmo->meR = 1.;

  /* Set for completeness. Nothing on the derivative path reads these: the
     expansion rate comes from CLASS's background, and the dark-energy and
     neutrino parameters exist only to feed rec_HubbleRate, which is never
     called. */
  cosmo->h      = 0.;
  cosmo->odeh2  = 0.;
  cosmo->orh2   = 0.;
  cosmo->onuh2  = 0.;
  cosmo->w0     = -1.;
  cosmo->wa     = 0.;
  cosmo->Nur    = Neff;
  cosmo->Nmnu   = 0.;
  cosmo->mnu[0] = cosmo->mnu[1] = cosmo->mnu[2] = 0.;
  cosmo->dlna                                   = (MODEL == SWIFT ? DLNA_SWIFT : DLNA_HYREC);
  cosmo->nz                                     = data_->d.Nz;

  /* Energy injection reaches HYREC-2 as already-deposited rates, filled per call
     in Derivatives(). CLASS owns the injection recipe and the channel split, so
     HYREC-2's own energy_injection.c never runs. */
  INJ_PARAMS* inj = cosmo->inj_params;
  std::memset(inj, 0, sizeof(INJ_PARAMS));
  inj->on_the_spot = 1;
}

HyrecModel::~HyrecModel() {
  hyrec_free(&data_->d);
  delete data_;
}

const char* HyrecModel::Name() const {
  return "HYREC-2 (" HYREC_VERSION ", SWIFT)";
}

IonisationDerivatives HyrecModel::Derivatives(const RecombinationState& state,
                                              const EnergyDeposition& dep,
                                              double energy_rate) const {
  IonisationDerivatives result = {0., 0.};

  const double nH_cgs   = state.n_H * 1e-6;
  const double T_mat_eV = state.T_mat * kBoltz;
  const double T_rad_eV = state.T_rad * kBoltz;
  INJ_PARAMS* inj       = data_->d.cosmo->inj_params;

  /* HYREC-2 applies varying alpha and m_e itself, throughout hydrogen.c and helium.c;
     they are z-dependent, so they are published on every call. */
  data_->d.cosmo->fsR = state.constants.alpha;
  data_->d.cosmo->meR = state.constants.me;

  /* HYREC-2 wants deposition rates per hydrogen atom per second, and applies the
     branching between ionization and Lyman-alpha escape itself -- unlike RECFAST,
     which needs the Peebles factor C to do the same job. */
  const double lya_rate = (energy_rate > 0.) ? energy_rate * dep.lya / (state.n_H * kE_H_lya * _eV_)
                                             : 0.;

  /************/
  /* hydrogen */
  /************/

  if (!state.hydrogen_frozen) {
    inj->ion    = (energy_rate > 0.) ? energy_rate * dep.ion_H / (state.n_H * kE_H_ion * _eV_) : 0.;
    inj->exclya = lya_rate;

    /* Outside the effective-rate tables HYREC-2 has no data for its multilevel
       atom, and asking anyway is a hard error rather than a graceful fallback.
       class_public gets away without the upper guard because its phase structure
       starts the hydrogen ODE below TR_MAX; ours starts wherever
       recfast_x_H0_trigger puts it, which is a precision parameter. Peebles is
       the right answer at both ends: above TR_MAX hydrogen is fully ionized, and
       below TR_MIN recombination has long finished. */
    /* The tables are in atomic units of temperature, which move with alpha^2 m_e. */
    const double T_rad_atomic = T_rad_eV / state.constants.TemperatureRescale();
    const bool outside_tables = (T_rad_atomic <= TR_MIN || T_rad_atomic >= TR_MAX ||
                                 state.T_mat / state.T_rad <= T_RATIO_MIN);
    const int model           = outside_tables ? PEEBLES : MODEL;

    /* iz indexes the radiation-field history, which only FULL mode keeps. */
    const double dxHIIdlna = rec_dxHIIdlna(&data_->d,
                                           model,
                                           state.x,
                                           state.x_H,
                                           nH_cgs,
                                           state.H,
                                           T_mat_eV,
                                           T_rad_eV,
                                           0,
                                           state.z);
    class_test(data_->d.error != 0, "HYREC-2 failed in rec_dxHIIdlna: %s", data_->d.error_message);

    result.dx_H_dz = -dxHIIdlna / (1. + state.z);
  }

  /************/
  /* helium   */
  /************/

  /* HYREC-2 counts HeII relative to hydrogen where CLASS counts it relative to
     total helium, hence the f_He on the way in and out. */
  const double xHeII = state.x_He * fHe_;

  /* HYREC-2's helium routine covers HeII -> HeI only, and needs neutral helium to
     exist: its escape probability goes as 1/xHeI. CLASS's analytic phase runs
     from HeIII (x_He = 2) down through the second Saha transition, which is
     outside that range, and holds helium at Saha rather than evolving it -- so
     the derivative is unused there anyway, entering only a negligible correction
     term in dT_mat/dz.

     Below XHEII_MIN HYREC-2 stops tracking helium and takes it to be entirely
     neutral (history.c). Here it relaxes to zero at kRelax e-folds per e-fold of
     the scale factor. Holding it at XHEII_MIN instead, as this wrapper used to,
     left 1e-6 free electrons per hydrogen in the late tail, and made XHEII_MIN a
     trap: with a zero derivative below it, a stiff step that landed there -- even
     below zero -- was self-consistent and stuck. The two rates are blended over
     [XHEII_MIN, 2 XHEII_MIN] so the derivative stays continuous, and HYREC-2 is
     never asked about helium at the low temperatures where its Boltzmann factors
     overflow. */
  if (state.helium_ode && xHeII < fHe_) {
    constexpr double kRelax = 100.;
    const double relax      = -kRelax * xHeII;
    double dxHeIIdlna       = relax;

    if (xHeII >= XHEII_MIN) {
      inj->ion = (energy_rate > 0.) ? energy_rate * dep.ion_He / (state.n_H * kE_H_ion * _eV_) : 0.;
      inj->exclya = lya_rate;

      double hyrec = rec_helium_dxHeIIdlna(&data_->d, state.z, 1. - state.x_H, xHeII, state.H);
      class_test(data_->d.error != 0,
                 "HYREC-2 failed in rec_helium_dxHeIIdlna: %s",
                 data_->d.error_message);

      if (He_rate_correction_) {
        hyrec *= HeliumRateCorrection(T_rad_eV / state.constants.TemperatureRescale());
      }

      const double weight = std::min(1., (xHeII - XHEII_MIN) / XHEII_MIN);
      dxHeIIdlna          = weight * hyrec + (1. - weight) * relax;
    }

    result.dx_He_dz = -dxHeIIdlna / (1. + state.z) / fHe_;
  }

  return result;
}
