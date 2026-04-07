#include "idm_dr_idr_species.h"
#include "background_module.h"
#include "perturbations_module.h"
#include "thermodynamics_module.h"
#include "thermodynamics.h"

IDM_DR_IDR_Species::IDM_DR_IDR_Species(const background& pba)
  : CompositeSpecies("IDM_DR_IDR", BaseSpecies::EnergyType::Other)
  , pba_(pba)
{
  auto idm = std::make_unique<IDM_DRSpecies>(pba);
  auto idr = std::make_unique<IDRSpecies>(pba);
  idm_dr_ = idm.get();
  idr_    = idr.get();
  children_.push_back(std::move(idm));
  children_.push_back(std::move(idr));
}

void IDM_DR_IDR_Species::AddCouplingDerivs(double /*tau*/, const double* y, double* dy,
                                            const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace*    ppw = ppaw.ppw;
  const perturb_vector*       pv  = ppw->pv;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  auto* pth_mod  = ppaw.perturbations_module->GetThermodynamicsModule().get();
  const double* pvecback   = ppw->pvecback;
  const double* pvecthermo = ppw->pvecthermo;
  auto* ppt = ppaw.perturbations_module->GetPerturbs();

  const double dmu_idm_dr = pvecthermo[pth_mod->index_th_dmu_idm_dr_];
  const double rho_idm_dr = idm_dr_->Rho(pvecback);
  const double rho_idr    = idr_->Rho(pvecback);

  // No coupling terms if either sector is absent or coupling rate vanishes.
  // This also guards the Sinv = rho_idr/rho_idm_dr and 1/dmu_idm_dr divisions.
  if (rho_idm_dr <= 0. || rho_idr <= 0. || dmu_idm_dr <= 0.) return;

  const double Sinv         = 4./3. * rho_idr / rho_idm_dr;
  const double theta_idm_dr = y[pv->index_pt_theta_idm_dr];

  // Under RSA for IDR: IDR Boltzmann hierarchy is not evolved; use RSA-approximated
  // theta_idr to keep the IDM_DR drag correct. Do not write to IDR equations.
  if (ppw->approx[ppw->index_ap_rsa_idr] == (int)rsa_idr_on) {
    const double theta_idr = ppw->rsa_theta_idr;
    dy[pv->index_pt_theta_idm_dr] -= Sinv * dmu_idm_dr * (theta_idm_dr - theta_idr)
                                    - ctx.k2 * pvecthermo[pth_mod->index_th_cidm_dr2_]
                                      * y[pv->index_pt_delta_idm_dr];
    return;
  }

  const double theta_idr = (pv->index_pt_theta_idr >= 0)
                           ? y[pv->index_pt_theta_idr] : 0.;

  if (ppw->approx[ppw->index_ap_tca_idm_dr] == (int)tca_idm_dr_off) {
    const thermo* pth = pth_mod->GetThermodynamics();
    const double dmu_idr = pth->b_idr / pth->a_idm_dr
                         * pba_.Omega0_idr / pba_.Omega0_idm_dr * dmu_idm_dr;

    // IDM_DR velocity coupling
    dy[pv->index_pt_theta_idm_dr] -= Sinv * dmu_idm_dr * (theta_idm_dr - theta_idr)
                                    - ctx.k2 * pvecthermo[pth_mod->index_th_cidm_dr2_]
                                      * y[pv->index_pt_delta_idm_dr];

    // IDR velocity coupling
    if (ctx.idr_nature == idr_free_streaming) {
      dy[pv->index_pt_theta_idr] += dmu_idm_dr * (theta_idm_dr - theta_idr);

      // IDR Compton collision terms in hierarchy l>=2
      const int l_max = pv->l_max_idr;
      for (int l = 2; l <= l_max; l++) {
        dy[pv->index_pt_delta_idr + l] -=
            (ppt->alpha_idm_dr[l-2] * dmu_idm_dr + ppt->beta_idr[l-2] * dmu_idr)
            * y[pv->index_pt_delta_idr + l];
      }
    }
  } else {
    // TCA on: compute tca_shear and tca_slip locally
    const double delta_idr = (pv->index_pt_delta_idr >= 0) ? y[pv->index_pt_delta_idr] : 0.;

    const double tca_shear_idm_dr =
        0.5 * 8./15. / dmu_idm_dr / ppt->alpha_idm_dr[0]
        * (theta_idm_dr + ctx.metric_shear);

    const double tca_slip_idm_dr =
        (pth_mod->GetThermodynamics()->nindex_idm_dr - 2./(1.+Sinv))
        * ctx.a_prime_over_a * (theta_idm_dr - theta_idr)
        + 1./(1.+Sinv) / dmu_idm_dr
        * (-ctx.a_prime_over_a * theta_idm_dr
           + ctx.k2 * pvecthermo[pth_mod->index_th_cidm_dr2_]
             * y[pv->index_pt_delta_idm_dr]
           + ctx.k2 * Sinv * (delta_idr/4. - tca_shear_idm_dr));

    // ASSIGN (=): TCA replaces the free-streaming velocity written by the children
    dy[pv->index_pt_theta_idm_dr] =
        1./(1.+Sinv) * (-ctx.a_prime_over_a*theta_idm_dr
                        + ctx.k2*pvecthermo[pth_mod->index_th_cidm_dr2_]
                          *y[pv->index_pt_delta_idm_dr]
                        + ctx.k2*Sinv*(delta_idr/4. - tca_shear_idm_dr))
        + ctx.metric_euler + Sinv/(1.+Sinv)*tca_slip_idm_dr;

    dy[pv->index_pt_theta_idr] =
        1./(1.+Sinv) * (-ctx.a_prime_over_a*theta_idm_dr
                        + ctx.k2*pvecthermo[pth_mod->index_th_cidm_dr2_]
                          *y[pv->index_pt_delta_idm_dr]
                        + ctx.k2*Sinv*(delta_idr/4. - tca_shear_idm_dr))
        + ctx.metric_euler - 1./(1.+Sinv)*tca_slip_idm_dr;
  }
}
