// Unit test for AxionNCDMSpecies (ncdm_axion): construction guards, the
// gstar_dec -> T mapping, Bose-Einstein PSD wiring, analytic moments, and an
// A/B check against ncdm_standard with the same PSD tabulated to file.
#include <cassert>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

#include "background.h"
#include "background_module.h"
#include "cosmology.h"
#include "parser.h"
#include "species/axion_ncdm_species.h"
#include "species/ncdm_species.h"
#include "species/species_build_context.h"

static NcdmSettings TestSettings() {
  NcdmSettings s{};
  s.h           = 0.67;
  s.T_cmb       = 2.7255;
  s.tol_ncdm    = 1e-3;
  s.tol_ncdm_bg = 1e-5;
  s.tol_M_ncdm  = 1e-8;
  return s;
}

static FileContent BaseFc() {
  FileContent fc;
  fc.set("ax.type", "ncdm_axion");
  fc.set("ax.m", "1.0");
  fc.set("ax.T", "0.39");
  return fc;
}

static bool Throws(FileContent& fc) {
  background pba{};
  pba.H0 = 2.2e-4;
  try {
    AxionNCDMSpecies sp(&fc, "ax", TestSettings(), &pba, nullptr);
    return false;
  }
  catch (const std::invalid_argument&) {
    return true;
  }
}

// Value checks (T <= 0, gstar_dec below today's g*S, ksi != 0): these reject a
// numeric parameter, so the species code raises them via the plain (runtime_error)
// severity macros — the sampler should reject the point, not abort the chain.
static bool ThrowsComputation(FileContent& fc) {
  background pba{};
  pba.H0 = 2.2e-4;
  try {
    AxionNCDMSpecies sp(&fc, "ax", TestSettings(), &pba, nullptr);
    return false;
  }
  catch (const std::runtime_error&) {
    return true;
  }
}

int main() {
  background pba{};
  pba.H0 = 2.2e-4;

  // ── Construction with explicit T; quadrature built; Omega resolved from m ──
  {
    FileContent fc = BaseFc();
    AxionNCDMSpecies sp(&fc, "ax", TestSettings(), &pba, nullptr);
    assert(sp.GetMassInElectronvolt() == 1.0);
    assert(sp.GetDeg() == 1.0);
    assert(sp.q_size() > 0 && sp.q_size_bg() > 0);
    assert(sp.GetOmega0() > 0.);
  }

  // ── gstar_dec -> T mapping: g*S = 10.75 must reproduce T = (43/11/10.75)^(1/3).
  //    T_ is protected, so probe it through the (T^4-scaling) massless Neff. ──
  {
    FileContent fg;
    fg.set("ax.type", "ncdm_axion");
    fg.set("ax.gstar_dec", "10.75");
    AxionNCDMSpecies sp_g(&fg, "ax", TestSettings(), &pba, nullptr);

    char t_buf[32];
    std::snprintf(t_buf, sizeof(t_buf), "%.16e", std::cbrt((43. / 11.) / 10.75));
    FileContent ft;
    ft.set("ax.type", "ncdm_axion");
    ft.set("ax.T", t_buf);
    AxionNCDMSpecies sp_t(&ft, "ax", TestSettings(), &pba, nullptr);

    assert(std::fabs(sp_g.GetNeff(0.) / sp_t.GetNeff(0.) - 1.) < 1e-12);
  }

  // ── Constructor guards ──
  {
    FileContent fc = BaseFc();
    fc.set("ax.gstar_dec", "30");  // both T and gstar_dec
    assert(Throws(fc));
  }
  {
    FileContent fc;
    fc.set("ax.type", "ncdm_axion");
    fc.set("ax.m", "1.0");  // neither T nor gstar_dec
    assert(Throws(fc));
  }
  {
    FileContent fc;
    fc.set("ax.type", "ncdm_axion");
    fc.set("ax.gstar_dec", "3.0");  // below today's g*S = 43/11
    assert(ThrowsComputation(fc));
  }
  {
    FileContent fc = BaseFc();
    fc.set("ax.ksi", "0.1");  // chemical potential unsupported
    assert(ThrowsComputation(fc));
  }
  {
    FileContent fc = BaseFc();
    fc.set("ax.use_psd_file", "1");  // file PSD contradicts the built-in BE PSD
    fc.set("ax.psd_filename", "psd_FD_single.dat");
    assert(Throws(fc));
  }
  {
    FileContent fc;
    fc.set("ax.type", "ncdm_axion");
    fc.set("ax.m", "1.0");
    fc.set("ax.T", "0.0");  // non-positive temperature
    assert(ThrowsComputation(fc));
  }

  // ── CreateAll factory: one instance per dot-syntax block of this type ──
  {
    FileContent fc        = BaseFc();
    NcdmSettings settings = TestSettings();
    SpeciesBuildContext ctx{};
    ctx.pfc           = &fc;
    ctx.pba           = &pba;
    ctx.ncdm_settings = &settings;
    auto result       = AxionNCDMSpecies::CreateAll(ctx);
    assert(result.size() == 1);
    assert(result[0].key == "ax");
    assert(dynamic_cast<AxionNCDMSpecies*>(result[0].species.get()) != nullptr);
  }

  // ── Omega-driven closure: with Omega given (no m), the base derives M via
  //    MFromOmega and back-computes m. Re-deriving Omega from the resulting
  //    species must reproduce the input. ──
  {
    FileContent fc;
    fc.set("ax.type", "ncdm_axion");
    fc.set("ax.Omega", "0.003");
    fc.set("ax.T", "0.39");
    AxionNCDMSpecies sp(&fc, "ax", TestSettings(), &pba, nullptr);
    assert(std::fabs(sp.GetOmega0() / 0.003 - 1.) < 1e-6);
    assert(sp.GetMassInElectronvolt() > 0.);
    double rho = 0.;
    sp.ComputeMomenta(0., nullptr, &rho, nullptr, nullptr, nullptr);
    const double H0 = 0.67 * 1.e5 / _c_;
    assert(std::fabs(rho / (H0 * H0) / 0.003 - 1.) < 1e-6);
  }

  // ── Massless limit: Neff = (4/7) (11/4)^(4/3) T^4 (one bosonic dof).
  //    GetNeff integrates at M = 0 exactly, so this pins the BE energy integral
  //    pi^4/15 against the code's rho_nu_rel_ normalization. ──
  {
    FileContent fc;
    fc.set("ax.type", "ncdm_axion");
    fc.set("ax.T", "0.39");
    AxionNCDMSpecies sp(&fc, "ax", TestSettings(), &pba, nullptr);
    const double expected = 4. / 7. * std::pow(11. / 4., 4. / 3.) * std::pow(0.39, 4);
    assert(std::fabs(sp.GetNeff(0.) / expected - 1.) < 1e-4);
  }

  // ── Statistics + dof normalization in one number:
  //    rho_BE(1 dof) / rho_FD(2 dof, ncdm_standard deg=1) = (pi^4/15)/(2 * 7pi^4/120) = 4/7 ──
  {
    FileContent fa;
    fa.set("ax.type", "ncdm_axion");
    fa.set("ax.T", "0.5");
    AxionNCDMSpecies ax(&fa, "ax", TestSettings(), &pba, nullptr);
    FileContent fn;
    fn.set("nu.type", "ncdm_standard");
    fn.set("nu.T", "0.5");
    NCDMSpecies nu(&fn, "nu", TestSettings(), &pba, nullptr);
    assert(std::fabs(ax.GetNeff(0.) / nu.GetNeff(0.) - 4. / 7.) < 1e-4);
  }

  // ── Omega h^2 against the analytic non-relativistic hot-relic value:
  //    m = 1 eV, T = 0.39 -> M ~ 1.1e4, so rho0 = n0 * m with
  //    n-integral 2*zeta3 (Bose-Einstein) to O((T/m)^2) ~ 1e-8.
  //    Omega h^2 = rho_class * (c/1e5)^2 since H0_class = h*1e5/c [1/Mpc]. ──
  {
    FileContent fc = BaseFc();
    AxionNCDMSpecies sp(&fc, "ax", TestSettings(), &pba, nullptr);
    const double T = 0.39, T_cmb = 2.7255, m_eV = 1.0;
    const double factor   = 4. * _PI_ * std::pow(T_cmb * T * _k_B_, 4) * 8. * _PI_ * _G_ / 3. /
                            std::pow(_h_P_ / 2. / _PI_, 3) / std::pow(_c_, 7) * _Mpc_over_m_ *
                            _Mpc_over_m_;
    const double M        = m_eV * _eV_ / (_k_B_ * T * T_cmb);
    const double rho0     = factor * (2. * _zeta3_ / std::pow(2. * _PI_, 3)) * M;
    const double expected = rho0 * (_c_ / 1.e5) * (_c_ / 1.e5);  // Omega h^2
    const double got      = sp.GetOmega0() * 0.67 * 0.67;
    assert(std::fabs(got / expected - 1.) < 1e-3);
  }

  // ── Internal consistency of the momenta integrals: rho/n -> M at z=0 for M >> 1 ──
  {
    FileContent fc = BaseFc();
    AxionNCDMSpecies sp(&fc, "ax", TestSettings(), &pba, nullptr);
    double n = 0., rho = 0.;
    sp.ComputeMomenta(0., &n, &rho, nullptr, nullptr, nullptr);
    assert(std::fabs(rho / n / sp.M() - 1.) < 1e-6);
  }

  // ── PSD wiring: dlnf0/dlnq = -q/(1 - e^-q) (Bose-Einstein) at every
  //    perturbation-grid node; -> -1 as q -> 0, -> -q at large q. ──
  {
    FileContent fc = BaseFc();
    AxionNCDMSpecies sp(&fc, "ax", TestSettings(), &pba, nullptr);
    const auto& q = sp.q();
    const auto& d = sp.dlnf0_dlnq();
    assert(q.size() == d.size());
    for (size_t i = 0; i < q.size(); ++i) {
      const double expected = -q[i] / (1. - std::exp(-q[i]));
      assert(std::fabs(d[i] - expected) < 1e-4 * std::fabs(expected));
    }
  }

  // ── Perturbation-grid quadrature: sum_i w_i q_i^3 = Int q^3 f0 dq = (pi^4/15)/(2pi)^3.
  //    (The w_ weights fold in f0; tol_ncdm = 1e-3 sets the grid accuracy.) ──
  {
    FileContent fc = BaseFc();
    AxionNCDMSpecies sp(&fc, "ax", TestSettings(), &pba, nullptr);
    const auto& q = sp.q();
    const auto& w = sp.w();
    double I3     = 0.;
    for (size_t i = 0; i < q.size(); ++i) {
      I3 += q[i] * q[i] * q[i] * w[i];
    }
    const double expected = std::pow(_PI_, 4) / 15. / std::pow(2. * _PI_, 3);
    assert(std::fabs(I3 / expected - 1.) < 1e-3);
  }

  // ── A/B: ncdm_standard + the same BE PSD tabulated to file must reproduce
  //    ncdm_axion. Log-spaced grid so the spline resolves the 1/q rise; the
  //    constant extrapolation below q_min contributes O(q_min^3) ~ 1e-12. ──
  {
    const std::string psd_path = "output/axion_be_psd_test.dat";
    {
      FILE* f = std::fopen(psd_path.c_str(), "w");
      assert(f != nullptr);
      const double kPref = 1.0 / std::pow(2. * _PI_, 3);
      const int n_rows   = 5000;
      const double q_lo = 1e-4, q_hi = 30.;
      for (int i = 0; i < n_rows; ++i) {
        const double q = q_lo * std::pow(q_hi / q_lo, i / (double) (n_rows - 1));
        std::fprintf(f, "%.16e %.16e\n", q, kPref / std::expm1(q));
      }
      std::fclose(f);
    }
    FileContent fa = BaseFc();
    AxionNCDMSpecies ax(&fa, "ax", TestSettings(), &pba, nullptr);
    FileContent fn;
    fn.set("nu.type", "ncdm_standard");
    fn.set("nu.m", "1.0");
    fn.set("nu.T", "0.39");
    fn.set("nu.use_psd_file", "1");
    fn.set("nu.psd_filename", psd_path);
    NCDMSpecies nu(&fn, "nu", TestSettings(), &pba, nullptr);
    for (double z : {0., 100., 3000.}) {
      double rho_ax = 0., rho_nu = 0.;
      ax.ComputeMomenta(z, nullptr, &rho_ax, nullptr, nullptr, nullptr);
      nu.ComputeMomenta(z, nullptr, &rho_nu, nullptr, nullptr, nullptr);
      assert(std::fabs(rho_ax / rho_nu - 1.) < 1e-4);
    }
    assert(std::fabs(ax.GetOmega0() / nu.GetOmega0() - 1.) < 1e-4);
    std::remove(psd_path.c_str());
  }

  // ── Full pipeline: ini-driven construction, background + spectra end-to-end.
  //    Omega_ax must match the same NR analytic value used above (m=1eV, T=0.39). ──
  {
    FileContent fc;
    fc.set("h", "0.67");
    fc.set("omega_b", "0.022");
    fc.set("omega_cdm", "0.12");
    fc.set("ax.type", "ncdm_axion");
    fc.set("ax.m", "1.0");
    fc.set("ax.T", "0.39");
    fc.set("output", "tCl,mPk");
    fc.set("l_max_scalars", "600");  // keep the smoke test fast
    fc.set("P_k_max_h/Mpc", "1.");
    Cosmology cosmo(fc);
    auto& bgm = cosmo.GetBackgroundModule();

    const double T = 0.39, T_cmb_v = 2.7255, m_eV = 1.0, h = 0.67;
    const double factor      = 4. * _PI_ * std::pow(T_cmb_v * T * _k_B_, 4) * 8. * _PI_ * _G_ / 3. /
                               std::pow(_h_P_ / 2. / _PI_, 3) / std::pow(_c_, 7) * _Mpc_over_m_ *
                               _Mpc_over_m_;
    const double M           = m_eV * _eV_ / (_k_B_ * T * T_cmb_v);
    const double expected_om = factor * (2. * _zeta3_ / std::pow(2. * _PI_, 3)) * M * (_c_ / 1.e5) *
                               (_c_ / 1.e5) / (h * h);
    const double om_ax       = bgm->GetOmega0Species("ax");
    assert(std::fabs(om_ax / expected_om - 1.) < 2e-3);

    auto& spm = cosmo.GetSpectraModule();
    (void) spm;  // reaching here without a throw = the full pipeline ran
    std::printf("full ncdm_axion pipeline ran, Omega_ax = %.4e\n", om_ax);
  }

  std::printf("axion ncdm test passed\n");
  return 0;
}
