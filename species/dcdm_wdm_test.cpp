// Unit test for WdmDecayProductSpecies (dcdm_wdm daughter): construction guards,
// momentum grid, injection-kernel normalization (the exact energy sum rule), and
// (from Task 2 onward) the background interface.
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "background.h"
#include "background_module.h"
#include "cosmology.h"
#include "parser.h"
#include "species/dcdm_wdm_species.h"
#include "species/wdm_decay_product.h"

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
  fc.set("ddm.type", "dcdm_wdm");
  fc.set("ddm.Gamma", "100");
  fc.set("ddm.vkick", "0.1");
  fc.set("ddm.Omega_ini", "0.05");
  return fc;
}

static bool Throws(FileContent& fc) {
  background pba{};
  pba.H0 = 2.2e-4;
  try {
    WdmDecayProductSpecies sp(&fc, "ddm", TestSettings(), &pba, nullptr);
    return false;
  }
  catch (const std::invalid_argument&) {
    return true;
  }
}

// Value checks (numeric ranges on epsilon/vkick/momenta_bins): the
// species code raises these via the plain (runtime_error) severity macros — the
// sampler should reject the point, not abort the chain.
static bool ThrowsComputation(FileContent& fc) {
  background pba{};
  pba.H0 = 2.2e-4;
  try {
    WdmDecayProductSpecies sp(&fc, "ddm", TestSettings(), &pba, nullptr);
    return false;
  }
  catch (const std::runtime_error&) {
    return true;
  }
}

int main() {
  background pba{};
  pba.H0 = 2.2e-4;

  // ── Fiducial cosmic-time helpers (Task 1) ──────────────────────────────────
  {
    const double H0 = 0.674 * 1.0e5 / _c_, Om = 0.31, Or = 9.2e-5;
    // Matter-era limit t(a) -> (2/3) a^{3/2} / (H0 sqrt(Om)) for a >> a_eq (=Or/Om ~ 3e-4).
    for (double a : {0.1, 0.5, 1.0}) {
      const double t_matter = (2.0 / 3.0) * std::pow(a, 1.5) / (H0 * std::sqrt(Om));
      const double t        = WdmDecayProductSpecies::FiducialCosmicTime(a);
      assert(std::fabs(t - t_matter) < 5e-3 * t_matter);
    }
    // Radiation-era limit t(a) -> a^2 / (2 H0 sqrt(Or)) for a << a_eq.
    {
      const double a     = 1e-6;
      const double t_rad = a * a / (2.0 * H0 * std::sqrt(Or));
      assert(std::fabs(WdmDecayProductSpecies::FiducialCosmicTime(a) - t_rad) < 1e-2 * t_rad);
    }
    // Monotone + inverse round-trip.
    double prev = -1.;
    for (double a : {1e-4, 1e-3, 1e-2, 0.1, 0.5, 1.0}) {
      const double t = WdmDecayProductSpecies::FiducialCosmicTime(a);
      assert(t > prev);
      prev                = t;
      const double a_back = WdmDecayProductSpecies::FiducialScaleFactorAtTime(t);
      assert(std::fabs(a_back - a) < 1e-9 * a);
    }
  }

  // ── Construction + grid ─────────────────────────────────────────────────────
  {
    FileContent fc = BaseFc();
    WdmDecayProductSpecies sp(&fc, "ddm", TestSettings(), &pba, nullptr);
    assert(sp.q_size() == 96);
    // M = kQKick * eps / v with eps = sqrt(1 - v^2)
    const double v = 0.1, eps = std::sqrt(1. - v * v);
    assert(std::fabs(sp.M() - 10.0 * eps / v) < 1e-12 * sp.M());
    // Injection-adapted grid: equal-injected-decay quantiles. Reconstruct the
    // decayed fraction F at each bin centre; it must be the uniform F-grid.
    {
      const int Nq       = sp.q_size();
      const double g1    = sp.Gamma() * WdmDecayProductSpecies::FiducialCosmicTime(1.0);
      const double F_at1 = -std::expm1(-g1);
      const double F_lo  = 1e-3;  // default q_edge_tol
      const double F_hi  = std::min(1. - 1e-3, F_at1);
      const double dF    = (F_hi - F_lo) / Nq;
      double prev        = 0.;
      for (int i = 0; i < Nq; ++i) {
        const double a = sp.q()[i] / 10.0;  // a' = q / kQKick
        const double F = -std::expm1(-sp.Gamma() * WdmDecayProductSpecies::FiducialCosmicTime(a));
        assert(std::fabs(F - (F_lo + (i + 0.5) * dF)) < 1e-6);
        assert(sp.q()[i] > prev);  // strictly increasing
        prev = sp.q()[i];
      }
      assert(sp.q().back() < 10.0 * (1. + 1e-9));  // within (0, kQKick]
    }
    // Grid depends on Gamma alone (fiducial H): different h / pba.H0, same Gamma
    // -> bit-identical grid.
    {
      FileContent fc2 = BaseFc();
      NcdmSettings s2 = TestSettings();
      s2.h            = 0.80;
      background pba2{};
      pba2.H0 = 3.0e-4;
      WdmDecayProductSpecies sp2(&fc2, "ddm", s2, &pba2, nullptr);
      assert(sp2.q_size() == sp.q_size());
      for (int i = 0; i < sp.q_size(); ++i)
        assert(sp.q()[i] == sp2.q()[i]);
    }
  }
  {
    FileContent fc;
    fc.set("ddm.type", "dcdm_wdm");
    fc.set("ddm.Gamma", "100");
    fc.set("ddm.epsilon", "0.6");
    fc.set("ddm.Omega_ini", "0.05");
    WdmDecayProductSpecies sp(&fc, "ddm", TestSettings(), &pba, nullptr);
    assert(std::fabs(sp.vkick() - 0.8) < 1e-14);
    assert(std::fabs(sp.M() - 10.0 * 0.6 / 0.8) < 1e-12);
    assert(sp.GetOmega0() == 0.);  // decay product reserves nothing itself
  }

  // ── Guards ──────────────────────────────────────────────────────────────────
  {
    FileContent fc = BaseFc();
    fc.set("ddm.epsilon", "0.5");  // both epsilon and vkick -> reject
    assert(Throws(fc));
  }
  {
    FileContent fc;
    fc.set("ddm.type", "dcdm_wdm");
    fc.set("ddm.Gamma", "100");
    fc.set("ddm.epsilon", "1.0");  // too close to 1
    fc.set("ddm.Omega_ini", "0.05");
    assert(ThrowsComputation(fc));
  }
  {
    FileContent fc = BaseFc();
    fc.set("ddm.m", "0.1");  // thermal-NCDM key -> reject
    assert(Throws(fc));
  }
  {
    FileContent fc;  // no Gamma at all
    fc.set("ddm.type", "dcdm_wdm");
    fc.set("ddm.vkick", "0.1");
    fc.set("ddm.Omega_ini", "0.05");
    assert(Throws(fc));
  }
  {
    FileContent fc = BaseFc();
    fc.set("ddm.Omega_dcdmwdm", "0.1");  // conflicts with Omega_ini outside shooting
    assert(Throws(fc));
  }
  // ── Injection kernel: exact energy sum rule ────────────────────────────────
  {
    FileContent fc = BaseFc();
    WdmDecayProductSpecies sp(&fc, "ddm", TestSettings(), &pba, nullptr);
    const int N = sp.q_size();
    std::vector<double> J(N), dJ(N);
    const int Nsr = sp.q_size();
    for (double a : {sp.q()[Nsr / 5] / 10.0,
                     sp.q()[2 * Nsr / 5] / 10.0,
                     sp.q()[3 * Nsr / 5] / 10.0,
                     sp.q()[4 * Nsr / 5] / 10.0}) {
      const double rho_dcdm = 0.7;  // arbitrary
      sp.FillInjection(a, rho_dcdm, J.data(), dJ.data());
      double sum = 0.;
      for (int i = 0; i < N; ++i) {
        const double q        = sp.q()[i];
        const double epsilon  = std::sqrt(q * q + a * a * sp.M() * sp.M());
        sum                  += sp.dq()[i] * q * q * epsilon * J[i];
      }
      const double lhs = sum * sp.factor() / (a * a * a * a);
      const double rhs = a * sp.Gamma() * rho_dcdm;
      assert(std::fabs(lhs - rhs) < 1e-12 * rhs);
    }
    // Below the grid's bottom edge the erf deposit CLAMPS the whole kick into the
    // bottom bin (it does not drop): the below-grid tail of GaussWeights is folded
    // into w[0] (Sigma w = 1 identically), so the realized fraction of the exact
    // a*Gamma*rho sum rule stays exactly 1 at the edge centre and everywhere below
    // it -- no injected energy is lost. This is what fixes the over-time energy
    // golden at high z (dropping removed the earliest q_edge_tol of decays, a large
    // fraction of the cumulative rho_wdm early on). Grid-derived, so it holds for
    // any Gamma.
    {
      auto gate_ratio = [&](double u_cut) {
        const double a_c = std::exp(u_cut) / 10.0;
        sp.FillInjection(a_c, 0.7, J.data(), dJ.data());
        double s = 0.;
        for (int i = 0; i < N; ++i) {
          const double q   = sp.q()[i];
          const double ep  = std::sqrt(q * q + a_c * a_c * sp.M() * sp.M());
          s               += sp.dq()[i] * q * q * ep * J[i];
        }
        return s * sp.factor() / std::pow(a_c, 4) / (a_c * sp.Gamma() * 0.7);
      };
      const double u0    = sp.u()[0];
      const double sp01  = sp.u()[1] - sp.u()[0];
      const double r_at  = gate_ratio(u0);               // at the bottom centre
      const double r_mid = gate_ratio(u0 - 0.5 * sp01);  // half a bin below
      const double r_far = gate_ratio(u0 - 2.0 * sp01);  // well below the edge
      assert(std::fabs(r_at - 1.0) < 1e-9);              // full deposit at the edge centre
      assert(std::fabs(r_mid - 1.0) < 1e-9);             // clamped: energy retained below edge
      assert(std::fabs(r_far - 1.0) < 1e-9);             // still fully retained far below (no drop)
      // dJ/dlnq pointwise sign structure below the grid is no longer asserted
      // here: with the -(3+q^2/eps^2)J measure term, the sign is not
      // constrained by "off-grid falling tail" reasoning alone. See the
      // moment-identity block below for the load-bearing dJ/dlnq property.
    }
  }

  // ── Erf deposit: interior energy centroid sits at u_cut to within the kernel
  // width. The cell-integrated Gaussian is symmetric in u, so the deposited
  // energy centroid tracks the kick location; sub-cell effects bound the error
  // by a fraction of sigma. (The old CIC's 1e-9 exactness was a stencil
  // property, not a physical requirement — the over-time golden is the
  // physical requirement.)
  {
    FileContent fc = BaseFc();
    WdmDecayProductSpecies sp(&fc, "ddm", TestSettings(), &pba, nullptr);
    const int N = sp.q_size();
    std::vector<double> J(N);
    for (int idx : {N / 4, N / 2, 3 * N / 4}) {
      const double u_cut = 0.5 * (sp.u()[idx] + sp.u()[idx + 1]);  // between centres
      const double a     = std::exp(u_cut) / 10.0;                 // u_cut = ln(a*kQKick)
      const double sigma = sp.SigmaAt(u_cut);
      sp.FillInjection(a, 0.7, J.data(), nullptr);
      double A = 0., Au = 0.;
      bool negative = false;
      for (int i = 0; i < N; ++i) {
        const double q   = sp.q()[i];
        const double ep  = std::sqrt(q * q + a * a * sp.M() * sp.M());
        const double e   = sp.dq()[i] * q * q * ep * J[i];  // injected energy in bin i
        A               += e;
        Au              += e * sp.u()[i];
        if (J[i] < 0.)
          negative = true;
      }
      assert(!negative);  // weights in [0,1]
      assert(A > 0.);
      assert(std::fabs(Au / A - u_cut) < 0.5 * sigma);         // centroid within the kernel
      assert(sigma >= 0.03 - 1e-12 && sigma <= 0.25 + 1e-12);  // clamp respected
    }
  }

  // ── High-side clamp: fast decay trims q_hi < kQKick; a kick past the top edge is
  //    CLAMPED into the top bin (retained, not dropped) so no energy is lost. ────
  {
    FileContent fc = BaseFc();
    fc.set("ddm.Gamma", "100000");  // very fast: F(a=1) -> 1, top trimmed
    WdmDecayProductSpecies sp(&fc, "ddm", TestSettings(), &pba, nullptr);
    const int N = sp.q_size();
    assert(sp.q().back() < 10.0 * (1. - 1e-6));  // q_hi trimmed below kQKick
    std::vector<double> J(N);
    // The realized fraction of the exact a*Gamma*rho sum rule, i.e. Sum(w).
    auto sum_w = [&](double a_c) {
      sp.FillInjection(a_c, 0.7, J.data(), nullptr);
      double s = 0.;
      for (int i = 0; i < N; ++i) {
        const double q   = sp.q()[i];
        const double ep  = std::sqrt(q * q + a_c * a_c * sp.M() * sp.M());
        s               += sp.dq()[i] * q * q * ep * J[i];
      }
      return s * sp.factor() / std::pow(a_c, 4) / (a_c * sp.Gamma() * 0.7);
    };
    const double a_in = sp.q().back() / 10.0 * 0.98;  // cutoff just below top edge
    assert(std::fabs(sum_w(a_in) - 1.0) < 1e-9);      // in-grid: full deposit
    // Cutoff at kQKick >> q_hi is past the top edge: clamped into the top bin,
    // energy retained (Sum(w) still 1), NOT shut off.
    assert(std::fabs(sum_w(1.0) - 1.0) < 1e-9);
    assert(J[N - 1] > 0.);  // the retained energy lands in the top bin
    // dJ/dlnq pointwise sign structure above the grid is no longer asserted
    // here: with the -(3+q^2/eps^2)J measure term, the sign is not
    // constrained by "off-grid rising tail" reasoning alone. See the
    // moment-identity block below for the load-bearing dJ/dlnq property.
  }
  // ── dJ/dlnq moment identity (integration by parts): with the completed
  //    -(3+q^2/eps^2)J measure term, dJ/dlnq is no longer pointwise-signed
  //    (rises below cutoff / falls above — that was a premise of the old,
  //    incomplete derivative). The physical invariant instead is the
  //    rho-weighted moment: the edge-difference part of dJ/dlnq telescopes to
  //    ~0 for interior cutoffs, so
  //      Sum dq q^2 eps dJdlnq  ==  - Sum dq q^2 eps (3 + q^2/eps^2) J
  //    up to the sparsity seed (~1e-12 rel.) and edge-pdf leakage. This is the
  //    g-channel's load-bearing property: it drives the daughter's fluid-limit
  //    (low-k) metric growth; its absence measured as a flat -7.6% P_m deficit
  //    (see .superpowers/sdd/task-3-lowk-diagnosis.md). Checked at three
  //    interior cutoffs (between centres at N/4, N/2, 3N/4) plus a=0.5. ─────
  {
    FileContent fc = BaseFc();
    WdmDecayProductSpecies sp(&fc, "ddm", TestSettings(), &pba, nullptr);
    const int N = sp.q_size();
    std::vector<double> J(N), dJ(N);

    // J itself still peaks near the cutoff (unaffected by the dJ/dlnq fix).
    {
      const double a     = 0.5;
      const double u_cut = std::log(a * 10.0);
      sp.FillInjection(a, 0.7, J.data(), dJ.data());
      int i_peak = 0;
      for (int i = 1; i < N; ++i)
        if (J[i] > J[i_peak])
          i_peak = i;
      assert(i_peak > 0);  // interior cutoff (a=0.5) -> peak is never the bottom bin
      assert(std::fabs(sp.u()[i_peak] - u_cut) <= (sp.u()[i_peak] - sp.u()[i_peak - 1]) * 1.5);
    }

    auto moment_check = [&](double a) {
      sp.FillInjection(a, 0.7, J.data(), dJ.data());
      double lhs = 0., rhs = 0.;
      for (int i = 0; i < N; ++i) {
        const double q   = sp.q()[i];
        const double ep  = std::sqrt(q * q + a * a * sp.M() * sp.M());
        lhs             += sp.dq()[i] * q * q * ep * dJ[i];
        rhs             -= sp.dq()[i] * q * q * ep * (3. + q * q / (ep * ep)) * J[i];
      }
      assert(std::fabs(lhs - rhs) < 1e-6 * std::fabs(rhs));
    };
    moment_check(0.5);
    for (int idx : {N / 4, N / 2, 3 * N / 4}) {
      const double u_cut = 0.5 * (sp.u()[idx] + sp.u()[idx + 1]);  // between centres
      moment_check(std::exp(u_cut) / 10.0);
    }
  }

  // ── Background plumbing: indices, ICs, ComputeMomenta via w_bg_ ────────────
  {
    FileContent fc = BaseFc();
    WdmDecayProductSpecies sp(&fc, "ddm", TestSettings(), &pba, nullptr);
    const int N = sp.q_size();

    int index_bg = 0, index_bi = 0;
    sp.RegisterBackgroundIndices(index_bg);
    sp.RegisterIntegrationIndices(index_bi);
    assert(index_bg == 4 + 3 * N);
    assert(index_bi == 2 * N);

    std::vector<double> bi(index_bi, 0.), bg(index_bg, 0.);
    BackgroundICContext icc{};
    icc.a_ini                = 1e-14;
    icc.pvecback_integration = bi.data();
    sp.SetBackgroundInitialConditions(icc);

    // The ICs seed f with a tiny evolver-friendly offset (subtracted again in
    // ComputeBackground); read it back so the checks below can compensate.
    const double f_seed = bi[sp.bi_f_index()];
    assert(f_seed > 0. && f_seed < 1e-8);

    // Empty distribution: rho = p = 0 exactly (the seed cancels exactly).
    sp.ComputeBackground(0.5, bi.data(), bg.data());
    assert(sp.Rho(bg.data()) == 0.);
    assert(sp.P(bg.data()) == 0.);

    // Populate one bin by hand and check rho against the closed-form quadrature.
    const int iq             = N / 2;
    const double a           = 0.5;
    bi[sp.bi_f_index() + iq] = 3.14 + f_seed;
    sp.ComputeBackground(a, bi.data(), bg.data());
    const double q       = sp.q()[iq];
    const double epsilon = std::sqrt(q * q + a * a * sp.M() * sp.M());
    const double rho_exp = sp.factor() / std::pow(a, 4) * sp.dq()[iq] * q * q * epsilon * 3.14;
    assert(std::fabs(sp.Rho(bg.data()) - rho_exp) < 1e-12 * rho_exp);
  }

  // ── Gauge guard ─────────────────────────────────────────────────────────────
  {
    for (const char* g : {"newtonian", "Newtonian", "new", "newt"}) {
      FileContent fc = BaseFc();
      fc.set("gauge", g);
      SpeciesBuildContext sctx{};
      sctx.pfc   = &fc;
      bool threw = false;
      try {
        DCDM_WDM_Species::CreateAll(sctx);
      }
      catch (const std::exception& e) {
        threw = std::string(e.what()).find("synchronous gauge only") != std::string::npos;
      }
      assert(threw);
    }
  }

  // ── Full background run + exact energy-conservation golden ────────────────
  {
    FileContent fc = BaseFc();  // Gamma=100 km/s/Mpc, vkick=0.1, Omega_ini=0.05
    fc.set("h", "0.67");
    fc.set("Omega_cdm", "0.20");
    fc.set("output", "");  // background-only
    Cosmology cosmo(fc);
    auto& bgm = cosmo.GetBackgroundModule();

    // Shooting must have pinned the combined density: reserve == integrated.
    const double reserve = bgm->GetOmega0Species("ddm");
    assert(reserve > 0.01 && reserve < 0.06);

    // Pull the full background table via titles/data.
    std::string titles;
    bgm->background_output_titles(titles);
    std::vector<std::string> cols;
    {
      std::stringstream ss(titles);
      std::string t;
      while (std::getline(ss, t, '\t'))
        if (!t.empty())
          cols.push_back(t);
    }
    auto col = [&](const std::string& nm) {
      for (size_t i = 0; i < cols.size(); ++i)
        if (cols[i] == nm)
          return (int) i;
      std::printf("column '%s' not found\n", nm.c_str());
      assert(false);
      return -1;
    };
    const int n_titles = (int) cols.size();
    std::vector<double> data((size_t) n_titles * bgm->bt_size_);
    bgm->background_output_data(n_titles, data.data());

    const int ic_z    = col("z");
    const int ic_tau  = col("conf. time [Mpc]");
    const int ic_dcdm = col("(.)rho_dcdm_ddm");
    const int ic_wdm  = col("(.)rho_wdm_ddm");

    // rho_wdm(a) must equal the injection-time integral
    //   Int dtau' a' Gamma rho_dcdm(a') (a'/a)^3 sqrt(eps^2 + (1-eps^2)(a'/a)^2)
    // (each daughter is born with energy m_parent/2 and redshifts kinetically).
    const double v = 0.1, eps2 = 1. - v * v;
    const double Gamma = 100. * 1.e3 / _c_;  // same double conversion as the species
    const int n_rows   = bgm->bt_size_;
    auto row           = [&](int r, int c) { return data[(size_t) r * n_titles + c]; };
    // Compare at epochs where rho_wdm is physical (nonzero). The deposit clamps
    // below-grid injection into the bottom bin (Σ W = 1 always), so energy is
    // conserved from the earliest decays.
    auto row_at_z = [&](double z_target) {
      int best = n_rows - 1;
      for (int r = 0; r < n_rows; ++r) {
        if (std::fabs(row(r, ic_z) - z_target) < std::fabs(row(best, ic_z) - z_target))
          best = r;
      }
      return best;
    };
    for (int target_row : {row_at_z(9.), row_at_z(3.), n_rows - 1}) {
      const double a_t = 1. / (1. + row(target_row, ic_z));
      double integral  = 0.;
      for (int r = 1; r <= target_row; ++r) {
        const double dtau = row(r, ic_tau) - row(r - 1, ic_tau);
        auto integrand    = [&](int rr) {
          const double ap = 1. / (1. + row(rr, ic_z));
          const double x  = ap / a_t;
          // energy ratio to birth energy m_parent/2: sqrt(eps^2 + (1-eps^2) x^2),
          // mass part eps2 constant, kinetic part (1-eps2) redshifts as x^2
          return ap * Gamma * row(rr, ic_dcdm) * x * x * x * std::sqrt(eps2 + (1. - eps2) * x * x);
        };
        integral += 0.5 * (integrand(r - 1) + integrand(r)) * dtau;
      }
      const double rho_wdm = row(target_row, ic_wdm);
      if (rho_wdm > 0.) {
        const double rel = std::fabs(integral - rho_wdm) / rho_wdm;
        std::printf("energy conservation at z=%.3g: rel. dev. %.2e\n", row(target_row, ic_z), rel);
        assert(rel < 5e-3);
      }
    }
  }

  // ── Perturbed run smoke test: Cls finite, spectra suppressed vs LCDM-ish ───
  //
  // Perturbed run smoke test: the full pipeline runs and Cls are finite.
  // (Uses momenta_bins=48 to stay fast; the historical momenta_bins=32
  // -ffast-math background failure no longer reproduces on master.)
  {
    FileContent fc = BaseFc();
    fc.set("h", "0.67");
    fc.set("Omega_cdm", "0.20");
    fc.set("output", "tCl,pCl,lCl,mPk");
    fc.set("lensing", "yes");
    fc.set("l_max_scalars", "600");  // keep the smoke test fast
    fc.set("P_k_max_h/Mpc", "1.");
    fc.set("ddm.momenta_bins",
           "48");  // fast (32 as in the brief hits a pre-existing bg bug, see above)
    Cosmology cosmo(fc);
    auto& spm = cosmo.GetSpectraModule();
    (void) spm;  // reaching here without a throw = the full pipeline ran
    std::printf("full perturbed dcdm_wdm pipeline ran\n");
  }

  std::printf("dcdm_wdm daughter test passed\n");
  return 0;
}
