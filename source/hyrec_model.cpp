#include "hyrec_model.h"

#include <cstdlib>
#include <cstring>

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

}  // namespace

HyrecModel::HyrecModel(const std::string& path,
                       double T_cmb,
                       double obh2,
                       double ocbh2,
                       double okh2,
                       double Neff,
                       double YHe,
                       double fHe,
                       double nH0,
                       double z_initial)
    : data_(new HyrecData), path_(path), fHe_(fHe) {
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

  /* No varying fundamental constants in CLASS++. */
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
    const bool outside_tables = (T_rad_eV <= TR_MIN || T_rad_eV >= TR_MAX ||
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
     term in dT_mat/dz. Below XHEII_MIN helium has finished recombining and
     HYREC-2 stops tracking it. */
  if (state.helium_ode && xHeII >= XHEII_MIN && xHeII < fHe_) {
    inj->ion = (energy_rate > 0.) ? energy_rate * dep.ion_He / (state.n_H * kE_H_ion * _eV_) : 0.;
    inj->exclya = lya_rate;

    const double dxHeIIdlna =
        rec_helium_dxHeIIdlna(&data_->d, state.z, 1. - state.x_H, xHeII, state.H);
    class_test(data_->d.error != 0,
               "HYREC-2 failed in rec_helium_dxHeIIdlna: %s",
               data_->d.error_message);

    result.dx_He_dz = -dxHeIIdlna / (1. + state.z) / fHe_;
  }

  return result;
}
