#include "recfast_model.h"

#include <algorithm>
#include <cmath>

IonisationDerivatives RecfastModel::Derivatives(const RecombinationState& state,
                                                const EnergyDeposition& dep,
                                                double energy_rate) const {
  IonisationDerivatives result = {0., 0.};

  /* Names the moved RECFAST code uses. */
  const double z                = state.z;
  const double x_H              = state.x_H;
  const double x_He             = state.x_He;
  const double x                = state.x;
  const double n                = state.n_H;
  const double n_He             = preco_->fHe * n;
  const double Tmat             = state.T_mat;
  const double Hz               = state.H;
  const bool hydrogen_frozen    = state.hydrogen_frozen;
  const bool helium_corrections = state.helium_ode;

  double Rdown = 1.e-19 * _a_PPB_ * pow((Tmat / 1.e4), _b_PPB_) /
                 (1. + _c_PPB_ * pow((Tmat / 1.e4), _d_PPB_));
  double Rup   = Rdown * pow((preco_->CR * Tmat), 1.5) * exp(-preco_->CDB / Tmat);

  double sq_0     = sqrt(Tmat / _T_0_);
  double sq_1     = sqrt(Tmat / _T_1_);
  double Rdown_He = _a_VF_ /
                    (sq_0 * pow((1. + sq_0), (1. - _b_VF_)) * pow((1. + sq_1), (1. + _b_VF_)));
  double Rup_He   = 4. * Rdown_He * pow((preco_->CR * Tmat), 1.5) * exp(-preco_->CDB_He / Tmat);
  double K        = preco_->CK / Hz;

  /* following is from recfast 1.5 */

  if (ppr_->recfast_Hswitch)
    K *= 1. +
         ppr_->recfast_AGauss1 *
             exp(-pow((log(1. + z) - ppr_->recfast_zGauss1) / ppr_->recfast_wGauss1, 2)) +
         ppr_->recfast_AGauss2 *
             exp(-pow((log(1. + z) - ppr_->recfast_zGauss2) / ppr_->recfast_wGauss2, 2));

  /* end of new recfast 1.5 piece */

  /* following is from recfast 1.4 */

  double Rdown_trip = _a_trip_ / (sq_0 * pow((1. + sq_0), (1. - _b_trip_)) *
                                  pow((1. + sq_1), (1. + _b_trip_)));
  double Rup_trip   = Rdown_trip * exp(-_h_P_ * _c_ * _L_He2St_ion_ / (_k_B_ * Tmat)) *
                      pow(preco_->CR * Tmat, 1.5) * 4. / 3.;

  int Heflag = 0;
  if (helium_corrections && (x_He >= 5.e-9))
    Heflag = ppr_->recfast_Heswitch;

  double CfHe_t = 0.;

  double K_He = 0.;
  if (Heflag == 0)
    K_He = preco_->CK_He / Hz;
  else {
    double tauHe_s = _A2P_s_ * preco_->CK_He * 3. * n_He * (1. - x_He) / Hz;
    double pHe_s   = (1. - exp(-tauHe_s)) / tauHe_s;
    K_He           = 1. / (_A2P_s_ * pHe_s * 3. * n_He * (1. - x_He));

    double Doppler = 0.;
    double pb      = 0.;
    double qb      = 0.;
    double AHcon   = 0.;

    /*    if (((Heflag == 2) || (Heflag >= 5)) && (x_H < 0.99999)) { */
    if (((Heflag == 2) || (Heflag >= 5)) &&
        (x_H < 0.9999999)) { /* threshold changed by Antony Lewis in 2008 to get smoother Helium */

      Doppler          = 2. * _k_B_ * Tmat / (_m_H_ * _not4_ * _c_ * _c_);
      Doppler          = _c_ * _L_He_2p_ * sqrt(Doppler);
      double gamma_2Ps = 3. * _A2P_s_ * preco_->fHe * (1. - x_He) * _c_ * _c_ /
                         (sqrt(_PI_) * _sigma_He_2Ps_ * 8. * _PI_ * Doppler * (1. - x_H)) /
                         pow(_c_ * _L_He_2p_, 2);
      pb               = 0.36;
      qb               = ppr_->recfast_fudge_He;
      AHcon            = _A2P_s_ / (1. + pb * pow(gamma_2Ps, qb));
      K_He             = 1. / ((_A2P_s_ * pHe_s + AHcon) * 3. * n_He * (1. - x_He));
    }

    if (Heflag >= 3) {
      double tauHe_t = _A2P_t_ * n_He * (1. - x_He) * 3. / (8. * _PI_ * Hz * pow(_L_He_2Pt_, 3));
      double pHe_t   = (1. - exp(-tauHe_t)) / tauHe_t;
      double CL_PSt  = _h_P_ * _c_ * (_L_He_2Pt_ - _L_He_2St_) / _k_B_;
      if ((Heflag == 3) || (Heflag == 5) || (x_H >= 0.99999)) {
        CfHe_t = _A2P_t_ * pHe_t * exp(-CL_PSt / Tmat);
        CfHe_t = CfHe_t / (Rup_trip + CfHe_t);
      }
      else {
        Doppler          = 2. * _k_B_ * Tmat / (_m_H_ * _not4_ * _c_ * _c_);
        Doppler          = _c_ * _L_He_2Pt_ * sqrt(Doppler);
        double gamma_2Pt = 3. * _A2P_t_ * preco_->fHe * (1. - x_He) * _c_ * _c_ /
                           (sqrt(_PI_) * _sigma_He_2Pt_ * 8. * _PI_ * Doppler * (1. - x_H)) /
                           pow(_c_ * _L_He_2Pt_, 2);
        pb               = 0.66;
        qb               = 0.9;
        AHcon            = _A2P_t_ / (1. + pb * pow(gamma_2Pt, qb)) / 3.;
        CfHe_t           = (_A2P_t_ * pHe_t + AHcon) * exp(-CL_PSt / Tmat);
        CfHe_t           = CfHe_t / (Rup_trip + CfHe_t);
      }
    }
  }

  /* end of new recfast 1.4 piece */

  /************/
  /* hydrogen */
  /************/

  if (hydrogen_frozen)
    result.dx_H_dz = 0.;
  else {
    /* Peebles' coefficient (approximated as one when the Hydrogen
       ionization fraction is very close to one) */
    double C = 0.;
    if (x_H < ppr_->recfast_x_H0_trigger2) {
      C = (1. + K * _Lambda_ * n * (1. - x_H)) /
          (1. / preco_->fu + K * _Lambda_ * n * (1. - x_H) / preco_->fu + K * Rup * n * (1. - x_H));
    }
    else {
      C = 1.;
    }

    /* For DM annihilation: fractions of the injected energy going into hydrogen
       ionization and into Lyman-alpha excitation. The escape probability C splits
       the excitation channel the way RECFAST's three-level atom sees it. */

    /* evolution of hydrogen ionisation fraction: */

    // JL: test for debugginf reio_inter
    //fprintf(stdout,"%e  %e  %e  %e\n",z,Tmat,K*_Lambda_*n,K*Rup*n);

    result.dx_H_dz = (x * x_H * n * Rdown - Rup * (1. - x_H) * exp(-preco_->CL / Tmat)) * C /
                         (Hz * (1. + z)) /* Peeble's equation with fudged factors */
                     -
                     energy_rate / n * (dep.ion_H / _L_H_ion_ + dep.lya * (1. - C) / _L_H_alpha_) /
                         (_h_P_ * _c_ * Hz *
                          (1. + z)); /* energy injection (neglect fraction going to helium) */
  }

  /************/
  /* helium   */
  /************/

  if (x_He < 1.e-15)
    result.dx_He_dz = 0.;
  else {
    double He_Boltz = exp(std::min(680., preco_->Bfact / Tmat));

    /* equations modified to take into account energy injection from dark matter */
    //C_He=(1. + K_He*_Lambda_He_*n_He*(1.-x_He)*He_Boltz)/(1. + K_He*(_Lambda_He_+Rup_He)*n_He*(1.-x_He)*He_Boltz);

    result.dx_He_dz =
        ((x * x_He * n * Rdown_He - Rup_He * (1. - x_He) * exp(-preco_->CL_He / Tmat)) *
         (1. + K_He * _Lambda_He_ * n_He * (1. - x_He) * He_Boltz)) /
        (Hz * (1 + z) *
         (1. +
          K_He * (_Lambda_He_ + Rup_He) * n_He * (1. - x_He) *
              He_Boltz)); /* in case of energy injection due to DM, we neglect the contribution to helium ionization */

    /* following is from recfast 1.4 */
    /* this correction is not self-consistent when there is energy injection  from dark matter, and leads to nan's  at small redshift (unimportant when reionization takes over before that redshift) */

    if (Heflag >= 3)
      result.dx_He_dz = result.dx_He_dz + (x * x_He * n * Rdown_trip -
                                           (1. - x_He) * 3. * Rup_trip *
                                               exp(-_h_P_ * _c_ * _L_He_2St_ / (_k_B_ * Tmat))) *
                                              CfHe_t / (Hz * (1. + z));

    /* end of new recfast 1.4 piece */
  }

  return result;
}
