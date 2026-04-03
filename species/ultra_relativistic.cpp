#include "ultra_relativistic.h"
#include "perturbations.h"
#include "perturbations_module.h"
#include "background.h"

UltraRelativisticSpecies::UltraRelativisticSpecies(const background& pba)
    : BaseSpecies("UR", EnergyType::Radiation),
      Omega0_ur_(pba.Omega0_ur), H0_(pba.H0) {}

// ── Background ────────────────────────────────────────────────────────────────

void UltraRelativisticSpecies::RegisterBackgroundIndices(int& index_bg) {
  class_define_index(index_bg_rho_, _TRUE_, index_bg, 1);
}

void UltraRelativisticSpecies::ComputeBackground(double a_rel,
                                                   const double* /*pvecback_B*/,
                                                   double* pvecback) {
  const double a4 = a_rel * a_rel * a_rel * a_rel;
  pvecback[index_bg_rho_] = Omega0_ur_ * H0_ * H0_ / a4;
}

double UltraRelativisticSpecies::Rho(const double* pvecback) const {
  return pvecback[index_bg_rho_];
}

double UltraRelativisticSpecies::P(const double* pvecback) const {
  return pvecback[index_bg_rho_] / 3.;
}

double UltraRelativisticSpecies::DpDloga(const double* pvecback) const {
  return -4./3. * pvecback[index_bg_rho_];
}

// ── Perturbations ─────────────────────────────────────────────────────────────

void UltraRelativisticSpecies::RegisterPerturbationIndices(
    perturb_vector* pv, int& index_pt,
    const perturb_workspace* ppw, int /*gauge*/) {

  /* If radiation streaming approximation is on, no UR variables are integrated.
     Sentinel -1 signals RSA to Delta/Theta/DeltaP/RhoPlusPShear. */
  if (ppw->approx[ppw->index_ap_rsa] == (int)rsa_on) {
    pv->index_pt_delta_ur = -1;
    pv->index_pt_theta_ur = -1;
    pv->index_pt_shear_ur = -1;
    pv->index_pt_l3_ur    = -1;
    return;
  }

  class_define_index(pv->index_pt_delta_ur, _TRUE_, index_pt, 1);
  class_define_index(pv->index_pt_theta_ur, _TRUE_, index_pt, 1);
  class_define_index(pv->index_pt_shear_ur, _TRUE_, index_pt, 1);

  if (ppw->approx[ppw->index_ap_ufa] == (int)ufa_off) {
    if (pv->l_max_ur >= 3) {
      class_define_index(pv->index_pt_l3_ur, _TRUE_, index_pt, pv->l_max_ur - 2);
    }
  }
}

void UltraRelativisticSpecies::PerturbDerivs(
    double tau, const double* y, double* dy,
    const perturb_parameters_and_workspace& ppaw) {

  const perturb_workspace* ppw    = ppaw.ppw;
  const perturb_vector* pv        = ppw->pv;
  const PerturbScalarContext& ctx  = ppw->scalar_ctx;
  const double* s_l               = ppw->s_l;
  const double  k                  = ppaw.k;
  const double  k2                 = k * k;
  const double  cotKgen            = ctx.cotKgen;
  const double  s2_squared         = ctx.s2_squared;
  const double  a_prime_over_a     = ctx.a_prime_over_a;
  const double  metric_continuity  = ctx.metric_continuity;
  const double  metric_euler       = ctx.metric_euler;
  const double  metric_shear       = ctx.metric_shear;
  const double  metric_ufa_class   = ctx.metric_ufa_class;

  const perturbs*   ppt = ppaw.perturbations_module->GetPerturbs();
  const background* pba = ppaw.perturbations_module->GetBackground();
  const precision*  ppr = ppaw.perturbations_module->GetPrecision();

  /* RSA active: UR quantities handled analytically, nothing to integrate */
  if (ppw->approx[ppw->index_ap_rsa] == (int)rsa_on) return;

  int l;
  const double three_ceff2 = ppt->three_ceff2_ur;
  const double three_cvis2 = ppt->three_cvis2_ur;

  /* ── density ─────────────────────────────────────────────────────────── */
  dy[pv->index_pt_delta_ur] =
      -4./3. * (y[pv->index_pt_theta_ur] + metric_continuity)
      + (1. - three_ceff2) * a_prime_over_a *
        (y[pv->index_pt_delta_ur] + 4. * a_prime_over_a * y[pv->index_pt_theta_ur] / k2);

  /* ── velocity ────────────────────────────────────────────────────────── */
  dy[pv->index_pt_theta_ur] =
      k2 * (three_ceff2 * y[pv->index_pt_delta_ur] / 4.
            - s2_squared * y[pv->index_pt_shear_ur])
      + metric_euler
      - (1. - three_ceff2) * a_prime_over_a * y[pv->index_pt_theta_ur];

  if (ppw->approx[ppw->index_ap_ufa] == (int)ufa_off) {
    /* ── exact shear ──────────────────────────────────────────────────── */
    dy[pv->index_pt_shear_ur] =
        0.5 * (8./15. * (y[pv->index_pt_theta_ur] + metric_shear)
               - 3./5. * k * s_l[3] / s_l[2] * y[pv->index_pt_shear_ur + 1]
               - (1. - three_cvis2) * 8./15.
                 * (y[pv->index_pt_theta_ur] + metric_shear));

    /* ── l = 3 ────────────────────────────────────────────────────────── */
    l = 3;
    dy[pv->index_pt_l3_ur] =
        k / (2.*l + 1.) *
        (l * 2. * s_l[l] * s_l[2] * y[pv->index_pt_shear_ur]
         - (l + 1.) * s_l[l + 1] * y[pv->index_pt_l3_ur + 1]);

    /* ── l = 4 .. l_max_ur - 1 ────────────────────────────────────────── */
    for (l = 4; l < pv->l_max_ur; l++) {
      dy[pv->index_pt_delta_ur + l] =
          k / (2.*l + 1) *
          (l * s_l[l] * y[pv->index_pt_delta_ur + l - 1]
           - (l + 1.) * s_l[l + 1] * y[pv->index_pt_delta_ur + l + 1]);
    }

    /* ── l = l_max_ur  (truncation) ───────────────────────────────────── */
    l = pv->l_max_ur;
    dy[pv->index_pt_delta_ur + l] =
        k * (s_l[l] * y[pv->index_pt_delta_ur + l - 1]
             - (1. + l) * cotKgen * y[pv->index_pt_delta_ur + l]);

    /* ── optional self-interaction (G_eff_ur) damping ─────────────────── */
    if (ppt->G_eff_ur > 0.) {
      const double a = ctx.a;
      double taudot = pow(a, -4)
                    * pow(pow(4./11., 1./3.) * pba->T_cmb * _k_B_, 5)
                    * pow(ppt->G_eff_ur / (1e12 * _eV_ * _eV_), 2)
                    * (2. * _PI_ / _h_P_) / _c_ * _Mpc_over_m_;
      taudot = std::min(taudot, a_prime_over_a * 1e9);
      if (taudot > 0.) {
        const double alpha_RTA[5] = {0.40, 0.43, 0.46, 0.47, 0.48};
        for (l = 2; l <= pv->l_max_ur; l++) {
          int alpha_index = std::min(4, l - 2);
          dy[pv->index_pt_delta_ur + l] -=
              alpha_RTA[alpha_index] * taudot * y[pv->index_pt_delta_ur + l];
        }
      }
    }

  } else {
    /* ── UFA: fluid approximation for shear only ──────────────────────── */
    const int method = ppr->ur_fluid_approximation;

    if (method == (int)ufa_mb) {
      dy[pv->index_pt_shear_ur] =
          -3. / tau * y[pv->index_pt_shear_ur]
          + 2./3. * (y[pv->index_pt_theta_ur] + metric_shear);
    } else if (method == (int)ufa_hu) {
      dy[pv->index_pt_shear_ur] =
          -3. * a_prime_over_a * y[pv->index_pt_shear_ur]
          + 2./3. * (y[pv->index_pt_theta_ur] + metric_shear);
    } else { /* ufa_CLASS (default) */
      dy[pv->index_pt_shear_ur] =
          -3. / tau * y[pv->index_pt_shear_ur]
          + 2./3. * (y[pv->index_pt_theta_ur] + metric_ufa_class);
    }

    /* optional G_eff_ur damping on shear */
    if (ppt->G_eff_ur > 0.) {
      const double a = ctx.a;
      double taudot = pow(a, -4)
                    * pow(pow(4./11., 1./3.) * pba->T_cmb * _k_B_, 5)
                    * pow(ppt->G_eff_ur / (1e12 * _eV_ * _eV_), 2)
                    * (2. * _PI_ / _h_P_) / _c_ * _Mpc_over_m_;
      taudot = std::min(taudot, a_prime_over_a * 1e9);
      dy[pv->index_pt_shear_ur] -= 0.40 * taudot * y[pv->index_pt_shear_ur];
    }
  }
}

double UltraRelativisticSpecies::Delta(const perturb_vector* pv, const double* y,
                                        const double* /*pvecback*/) const {
  return (pv->index_pt_delta_ur >= 0) ? y[pv->index_pt_delta_ur] : 0.;
}

double UltraRelativisticSpecies::Theta(const perturb_vector* pv, const double* y,
                                        const double* /*pvecback*/) const {
  return (pv->index_pt_theta_ur >= 0) ? y[pv->index_pt_theta_ur] : 0.;
}

double UltraRelativisticSpecies::DeltaP(const perturb_vector* pv, const double* y,
                                         const double* pvecback) const {
  /* delta_p = c_s^2 * delta_rho = (1/3) * rho * delta */
  return (pv->index_pt_delta_ur >= 0) ? Rho(pvecback) * y[pv->index_pt_delta_ur] / 3. : 0.;
}

double UltraRelativisticSpecies::RhoPlusPShear(const perturb_vector* pv, const double* y,
                                                const double* pvecback) const {
  /* (rho + p) * sigma = 4/3 * rho_ur * shear_ur */
  return (pv->index_pt_shear_ur >= 0) ? 4./3. * Rho(pvecback) * y[pv->index_pt_shear_ur] : 0.;
}
