// Unit test for DrPsdSpecies (dr_psd massless daughter): construction + grid,
// background plumbing (index counts, closed-form Rho, live dlnf/dlnq spline),
// free-streaming radiation redshift, the statistics guard, the #372
// CopyPerturbationsAcrossSwitch regression, and background-column-name
// disambiguation between the fermion/boson daughters sharing an instance name.
#include <cassert>
#include <cmath>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "background.h"
#include "decay_transition_kernel.h"  // Statistics
#include "parser.h"
#include "species/background_ic_context.h"
#include "species/dr_psd_species.h"
#include "species/ncdm_base_species.h"
#include "species/perturb_source_context.h"
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

// Construct a standalone dr_psd directly (bypasses the factory) for plumbing
// checks; `extra`-keyed FileContent lets a caller add e.g. an abundance key.
static std::unique_ptr<DrPsdSpecies> MakeDrPsd(FileContent& fc, Statistics stat) {
  background pba{};
  pba.H0 = 2.2e-4;
  return std::make_unique<DrPsdSpecies>(&fc, "drx", TestSettings(), &pba, nullptr, stat);
}

static FileContent BaseFc() {
  FileContent fc;
  fc.set("drx.type", "dr_psd");
  return fc;
}

// Distinct nonzero marker per (iq, l) slot (mirrors dncdm_switch_copy_test.cpp).
static double Marker(int iq, int l) {
  return 1. + 100. * iq + l;
}

static void FillNcdmLayout(NCDMBaseSpecies::PerturbLayout& l, int q_size, int l_max, int idx0) {
  l.q_size = q_size;
  l.l_max  = l_max;
  l.index_per_q.clear();
  for (int iq = 0; iq < q_size; ++iq)
    l.index_per_q.push_back(idx0 + iq * (l_max + 1));
}

int main() {
  // ── Construction + grid ────────────────────────────────────────────────────
  {
    FileContent fc = BaseFc();
    auto sp        = MakeDrPsd(fc, Statistics::Boson);
    assert(sp->statistics() == Statistics::Boson);
    assert(sp->M() == 0.);  // massless
    const int N = sp->q_size();
    assert(N >= 2);
    // Default dr_bg_refine == 1: the perturbation grid IS the background grid.
    assert(sp->refine() == 1);
    assert((int) sp->q_bg().size() == N);
    assert((int) sp->q().size() == N);
    assert((int) sp->dq_bg().size() == N);
    for (int i = 0; i < N; ++i) {
      assert(sp->q()[i] == sp->q_bg()[i]);    // bit-for-bit, not "close"
      assert(sp->dq()[i] == sp->dq_bg()[i]);  // refine=1 must be a true no-op
    }
    // log spacing: strictly increasing q, constant ratio q_{i+1}/q_i.
    const double r0 = sp->q_bg()[1] / sp->q_bg()[0];
    for (int i = 0; i + 1 < N; ++i) {
      assert(sp->q_bg()[i + 1] > sp->q_bg()[i]);
      assert(std::fabs(sp->q_bg()[i + 1] / sp->q_bg()[i] - r0) < 1e-10 * r0);
    }
  }

  // ── dr_bg_refine: PT grid is an EXACT integer subsample of the BG grid ──────
  // The perturbation hierarchy costs q_size x (l_max+1) variables per k-mode, so it
  // must be able to run coarser than the background grid the injection cliff needs.
  // The coarse grid is built by DIRECT INDEXING q_bg_[j*m] -- never a fresh
  // log-spaced formula, never a spline of background columns -- so a PT-grid point
  // is literally the same array entry (and the same background column) as a BG-grid
  // point. Points span N-1 intervals, hence N_bg-1 = m*(N_pt-1).
  {
    FileContent fc = BaseFc();
    fc.set("drx.dr_N_q", "100");  // 99 intervals = 9 * 11
    fc.set("drx.dr_bg_refine", "9");
    auto sp = MakeDrPsd(fc, Statistics::Fermion);
    assert(sp->refine() == 9);
    assert((int) sp->q_bg().size() == 100);
    assert(sp->q_size() == 12);  // (100-1)/9 + 1, NOT 100/9
    assert((int) sp->dq().size() == 12);

    for (int j = 0; j < sp->q_size(); ++j)
      assert(sp->q()[j] == sp->q_bg()[j * 9]);  // bit-for-bit: same array entry
    assert(sp->q().front() == sp->q_bg().front());
    assert(sp->q().back() == sp->q_bg().back());  // endpoints coincide -> same band cover

    // Same half-weight-endpoint log-trapezoid convention as dq_bg_, at h_pt = 9*h_bg.
    const double h_bg = std::log(sp->q_bg().back() / sp->q_bg().front()) / 99.;
    const double h_pt = 9. * h_bg;
    for (int j = 0; j < sp->q_size(); ++j) {
      const double wq     = (j == 0 || j == sp->q_size() - 1) ? 0.5 : 1.0;
      const double expect = wq * h_pt * sp->q()[j];
      assert(std::fabs(sp->dq()[j] - expect) <= 1e-13 * expect);
    }
    // Both grids integrate d(ln q) over the SAME range: sum dq/q = ln(qmax/qmin).
    double s_pt = 0., s_bg = 0.;
    for (int j = 0; j < sp->q_size(); ++j)
      s_pt += sp->dq()[j] / sp->q()[j];
    for (int i = 0; i < (int) sp->q_bg().size(); ++i)
      s_bg += sp->dq_bg()[i] / sp->q_bg()[i];
    assert(std::fabs(s_pt - s_bg) < 1e-12 * s_bg);
  }

  // ── dr_q_sampling = linear: quadrature.c qm_trapz, byte-compatible with the
  // reference implementations ────────────────────────────────────────────────
  // The published arXiv:2011.01502 background figures were produced on
  // `np.linspace(0, q_max, N+1)[1:]` with q_max=22, N=240 -- i.e. quadrature.c's
  // qm_trapz qmin==0 branch. Reproducing that grid EXACTLY is what makes a
  // like-for-like cross-code comparison possible, so pin it here rather than
  // merely asserting "uniform".
  {
    FileContent fc = BaseFc();
    fc.set("drx.dr_q_sampling", "linear");
    fc.set("drx.dr_q_max", "22");
    fc.set("drx.dr_N_q", "240");
    auto sp = MakeDrPsd(fc, Statistics::Fermion);

    const int N    = 240;
    const double h = 22. / N;
    assert(sp->q_size() == N);
    assert((int) sp->q_bg().size() == N);

    // q_i = (i+1)*h, with the endpoint made exact.
    for (int i = 0; i < N; ++i) {
      const double expect = (i == N - 1) ? 22. : h + i * h;
      assert(std::fabs(sp->q_bg()[i] - expect) <= 1e-12 * expect);
    }
    assert(sp->q_bg().back() == 22.);  // exact, as the reference codes force

    // dq_i = h, LAST point half-weighted only: q=0 is carried as an implicit
    // zero-weight point, so the first interior point keeps full weight. (This is
    // the asymmetry that distinguishes qm_trapz from a plain closed trapezoid --
    // getting it wrong shifts every moment by O(h).)
    for (int i = 0; i < N; ++i) {
      const double expect = (i == N - 1) ? 0.5 * h : h;
      assert(std::fabs(sp->dq_bg()[i] - expect) <= 1e-13 * h);
    }

    // refine defaults to 1, so the PT grid IS the BG grid.
    assert(sp->refine() == 1);
    for (int i = 0; i < N; ++i) {
      assert(sp->q()[i] == sp->q_bg()[i]);  // bit-for-bit
      assert(sp->dq()[i] == sp->dq_bg()[i]);
      assert(sp->bg_index(i) == i);
    }
    // Weights integrate dq over [0, q_max]: sum dq = q_max - h/2.
    double s = 0.;
    for (int i = 0; i < N; ++i)
      s += sp->dq_bg()[i];
    assert(std::fabs(s - (22. - 0.5 * h)) < 1e-12 * 22.);
  }

  // ── linear grid with dr_bg_refine > 1 ──────────────────────────────────────
  // The coarse grid must be a canonical qm_trapz grid IN ITS OWN RIGHT, of step
  // refine*h. That is what forces the subsample to be taken from the top,
  // q_pt[j] = q_bg[(j+1)*refine - 1]: anchoring on the implicit zero the way the
  // log grid does (q_bg[j*refine]) would leave a first cell of width h among cells
  // of width refine*h, which is not a trapezoid rule and shifts every moment.
  {
    FileContent fc = BaseFc();
    fc.set("drx.dr_q_sampling", "linear");
    fc.set("drx.dr_q_max", "22");
    fc.set("drx.dr_N_q", "240");
    fc.set("drx.dr_bg_refine", "6");
    auto sp = MakeDrPsd(fc, Statistics::Fermion);

    const int N_pt    = 240 / 6;
    const double h_pt = 22. / N_pt;
    assert(sp->q_size() == N_pt);
    for (int j = 0; j < N_pt; ++j) {
      // Canonical qm_trapz at the COARSE step: q_j = (j+1)*h_pt, only the last
      // point half-weighted.
      const double expect_q = (j == N_pt - 1) ? 22. : h_pt + j * h_pt;
      assert(std::fabs(sp->q()[j] - expect_q) <= 1e-12 * expect_q);
      const double expect_dq = (j == N_pt - 1) ? 0.5 * h_pt : h_pt;
      assert(std::fabs(sp->dq()[j] - expect_dq) <= 1e-13 * h_pt);
      // and still an EXACT subsample, so background columns stay directly indexable
      // (never interpolated across the injection cliff).
      assert(sp->bg_index(j) == (j + 1) * 6 - 1);
      assert(sp->q()[j] == sp->q_bg()[sp->bg_index(j)]);  // bit-for-bit
    }
    assert(sp->q().back() == 22.);              // endpoint inherited exactly
    assert(sp->bg_index(N_pt - 1) == 240 - 1);  // last PT point IS last BG point
    double s = 0.;
    for (int j = 0; j < N_pt; ++j)
      s += sp->dq()[j];
    assert(std::fabs(s - (22. - 0.5 * h_pt)) < 1e-12 * 22.);
  }

  // Divisibility is sampling-dependent, because the log grid counts INTERVALS from
  // its endpoint while the linear grid counts CELLS from its implicit zero: log
  // needs refine | N-1, linear needs refine | N. An unknown sampling name is also
  // rejected. All numeric (sampler-varyable) -> runtime_error, not invalid_argument.
  {
    auto ctor_throws_runtime = [](const char* sampling, const char* refine) {
      FileContent fc = BaseFc();
      fc.set("drx.dr_q_sampling", sampling);
      fc.set("drx.dr_N_q", "100");
      fc.set("drx.dr_bg_refine", refine);
      try {
        MakeDrPsd(fc, Statistics::Fermion);
      }
      catch (const std::runtime_error&) {
        return true;
      }
      catch (...) {
        return false;
      }
      return false;
    };
    // N = 100: log wants refine | 99, linear wants refine | 100.
    assert(ctor_throws_runtime("linear", "9"));     // 9 does not divide 100
    assert(!ctor_throws_runtime("linear", "10"));   // 10 does
    assert(ctor_throws_runtime("log", "10"));       // 10 does not divide 99
    assert(!ctor_throws_runtime("log", "9"));       // 9 does
    assert(ctor_throws_runtime("quadratic", "1"));  // unknown sampling name
    assert(!ctor_throws_runtime("linear", "1"));
  }

  // ── per-daughter initial occupation: the published A = [0, 1, 1] ───────────
  // arXiv:2011.01502's background figures start phi EMPTY and nu_l at a full
  // Fermi-Dirac. Both daughters are built by the composite under one instance
  // name, so before dr_f_ini_l / dr_f_ini_phi existed this was unreachable: the
  // seeding flag was shared and could only be all-empty or all-thermal.
  {
    auto seed_max = [](FileContent& fc, Statistics stat) {
      auto sp      = MakeDrPsd(fc, stat);
      int index_bg = 0, index_bi = 0;
      sp->RegisterBackgroundIndices(index_bg);
      sp->RegisterIntegrationIndices(index_bi);
      std::vector<double> bi(index_bi, 0.);
      BackgroundICContext icc{};
      icc.pvecback_integration = bi.data();
      sp->SetBackgroundInitialConditions(icc);
      double m = 0.;
      for (int i = 0; i < sp->q_size(); ++i)
        m = std::max(m, bi[sp->bi_f_index() + i]);
      return m;
    };

    // A = [phi:0, nu_l:1] -- asymmetric, the configuration that matters.
    {
      FileContent fc = BaseFc();
      fc.set("drx.dr_f_ini_l", "1");
      fc.set("drx.dr_f_ini_phi", "0");
      const double f_l   = seed_max(fc, Statistics::Fermion);
      const double f_phi = seed_max(fc, Statistics::Boson);
      assert(f_l > 1e-3);           // nu_l thermally populated
      assert(f_phi <= 2. * 1e-30);  // phi empty (only the positivity floor)
      assert(f_l > 1e20 * f_phi);   // and unambiguously so
    }
    // The amplitude is linear: A=0.5 halves the seed at every q.
    {
      FileContent fc1 = BaseFc();
      fc1.set("drx.dr_f_ini_l", "1");
      FileContent fc2 = BaseFc();
      fc2.set("drx.dr_f_ini_l", "0.5");
      const double a1 = seed_max(fc1, Statistics::Fermion);
      const double a2 = seed_max(fc2, Statistics::Fermion);
      assert(std::fabs(a2 / a1 - 0.5) < 1e-9);
    }
    // Back-compat: a legacy thermal key still seeds BOTH daughters at full
    // equilibrium, and no key at all still leaves both empty.
    {
      FileContent fc = BaseFc();
      fc.set("drx.dr_ksi", "0");
      assert(seed_max(fc, Statistics::Fermion) > 1e-3);
      assert(seed_max(fc, Statistics::Boson) > 1e-3);
      FileContent bare = BaseFc();
      assert(seed_max(bare, Statistics::Fermion) <= 2. * 1e-30);
      assert(seed_max(bare, Statistics::Boson) <= 2. * 1e-30);
    }
    // Negative amplitude is a rejectable numeric error, not a structural one.
    {
      FileContent fc = BaseFc();
      fc.set("drx.dr_f_ini_l", "-1");
      bool threw = false;
      try {
        MakeDrPsd(fc, Statistics::Fermion);
      }
      catch (const std::runtime_error&) {
        threw = true;
      }
      assert(threw);
    }
  }

  // ── dr_bg_refine guards (class_test -> runtime_error, NOT invalid_argument) ──
  // The divisibility requirement couples two numeric inputs a sampler can vary, so
  // it is a rejectable computation error, not a structural (severe) one.
  {
    auto ctor_throws_runtime = [](const char* n_q, const char* refine) {
      FileContent fc = BaseFc();
      fc.set("drx.dr_N_q", n_q);
      fc.set("drx.dr_bg_refine", refine);
      try {
        MakeDrPsd(fc, Statistics::Fermion);
      }
      catch (const std::runtime_error&) {
        return true;
      }
      catch (...) {
        return false;  // invalid_argument (severe) would be the wrong channel
      }
      return false;
    };
    assert(ctor_throws_runtime("100", "7"));  // 99 % 7 != 0
    assert(ctor_throws_runtime("100", "0"));  // refine must be >= 1
    assert(ctor_throws_runtime("100", "-3"));
    assert(!ctor_throws_runtime("100", "11"));  // 99 = 11*9 -> N_pt = 10
    assert(!ctor_throws_runtime("97", "8"));    // 96 = 8*12 -> N_pt = 13
  }

  // ── Background plumbing: index counts, IC seed, closed-form Rho ─────────────
  {
    FileContent fc = BaseFc();
    auto sp        = MakeDrPsd(fc, Statistics::Fermion);
    const int N    = sp->q_size();

    int index_bg = 0, index_bi = 0;
    sp->RegisterBackgroundIndices(index_bg);
    sp->RegisterIntegrationIndices(index_bi);
    assert(index_bg == 4 + 2 * N);  // number,rho,p,pseudo_p + f(N) + dlnfdlnq(N)
    assert(index_bi == N);          // f(N) only

    std::vector<double> bi(index_bi, 0.), bg(index_bg, 0.);
    BackgroundICContext icc{};
    icc.a_ini                = 1e-14;
    icc.pvecback_integration = bi.data();
    sp->SetBackgroundInitialConditions(icc);
    const double f_seed = bi[sp->bi_f_index()];
    assert(f_seed > 0. && f_seed < 1e-8);  // empty daughter starts at the evolver seed

    // Populate one bg bin by hand and match rho to the closed-form quadrature
    // rho = factor/a^4 · dq q² ε f with ε = q (massless).
    const int iq              = N / 2;
    const double a            = 0.5;
    bi[sp->bi_f_index() + iq] = 3.14 + f_seed;
    sp->ComputeBackground(a, bi.data(), bg.data());
    const double q       = sp->q_bg()[iq];
    const double rho_exp = sp->factor() / std::pow(a, 4) * sp->dq_bg()[iq] * q * q * q * 3.14;
    assert(std::fabs(sp->Rho(bg.data()) - rho_exp) < 1e-10 * rho_exp);
    // Massless equation of state: p == rho/3.
    assert(std::fabs(sp->P(bg.data()) - sp->Rho(bg.data()) / 3.) < 1e-10 * sp->Rho(bg.data()));
  }

  // ── ∂f̄/∂lnq: pure power-law f ∝ q^p splines to p·q^p at interior bins ───────
  // The published column is ∂f̄/∂lnq (splined from f̄ itself, NOT f̄·∂lnf̄/∂lnq — the
  // ln form rings across the kFFloor plateau of a decay-filled PSD), so the
  // power-law check is dimensional: ∂(q^p)/∂lnq = p·q^p. Tolerance is RELATIVE
  // because q^p spans decades across the grid. Only the INTERIOR is asserted: the
  // ln f̄ form was exact for a power law (ln f̄ is linear in ln q), the f̄ form carries
  // the usual _SPLINE_EST_DERIV_ end-condition bias, which decays geometrically off
  // each edge (1.5e-2, 5.1e-3, 1.8e-3, 5.7e-4, 2.1e-4, then ≤5e-5). Harmless where
  // it lands: the daughter grid edges are the empty tails of a decay-filled PSD.
  {
    FileContent fc = BaseFc();
    auto sp        = MakeDrPsd(fc, Statistics::Fermion);
    const int N    = sp->q_size();

    int index_bg = 0, index_bi = 0;
    sp->RegisterBackgroundIndices(index_bg);
    sp->RegisterIntegrationIndices(index_bi);
    std::vector<double> bi(index_bi, 0.), bg(index_bg, 0.);
    BackgroundICContext icc{};
    icc.pvecback_integration = bi.data();
    sp->SetBackgroundInitialConditions(icc);
    const double f_seed = bi[sp->bi_f_index()];

    const double p = -2.0;
    for (int i = 0; i < N; ++i)
      bi[sp->bi_f_index() + i] = std::pow(sp->q_bg()[i], p) + f_seed;
    sp->ComputeBackground(0.5, bi.data(), bg.data());
    double worst = 0.;
    for (int i = 5; i < N - 5; ++i) {
      const double expect = p * std::pow(sp->q_bg()[i], p);
      worst = std::max(worst,
                       std::fabs(bg[sp->bg_dfdlnq_index() + i] - expect) / std::fabs(expect));
    }
    std::printf("  dfdlnq q^%g: worst interior rel err = %.2e (need < 1e-4)\n", p, worst);
    assert(worst < 1e-4);
  }

  // ── Free-streaming: thermal seed, massless ρ ∝ a⁻⁴ radiation redshift ───────
  // v1 unit-level check (see DEVIATION note in the file header comment): a
  // standalone dr_psd reserves no sector (GetOmega0()==0), so a fully-closed
  // standalone cosmology is ill-posed by design (the composite reserves it, §5).
  // The physical invariant — massless free-streaming redshifts as radiation — is
  // asserted directly: a static comoving f gives ρ·a⁴ = const and p == ρ/3.
  {
    FileContent fc = BaseFc();
    fc.set("drx.dr_T", "1.0");  // request a thermal initial abundance (BE)
    auto sp     = MakeDrPsd(fc, Statistics::Boson);
    const int N = sp->q_size();

    int index_bg = 0, index_bi = 0;
    sp->RegisterBackgroundIndices(index_bg);
    sp->RegisterIntegrationIndices(index_bi);
    std::vector<double> bi(index_bi, 0.), bg(index_bg, 0.);
    BackgroundICContext icc{};
    icc.pvecback_integration = bi.data();
    sp->SetBackgroundInitialConditions(icc);
    // Thermal seed populated the bins well above the bare evolver seed.
    double max_f = 0.;
    for (int i = 0; i < N; ++i)
      max_f = std::max(max_f, bi[sp->bi_f_index() + i]);
    assert(max_f > 1e-3);

    double rho_a4_ref = -1.;
    for (double a : {0.25, 0.5, 1.0}) {
      sp->ComputeBackground(a, bi.data(), bg.data());
      const double rho_a4 = sp->Rho(bg.data()) * std::pow(a, 4);
      if (rho_a4_ref < 0.)
        rho_a4_ref = rho_a4;
      assert(std::fabs(rho_a4 - rho_a4_ref) < 1e-10 * rho_a4_ref);  // ρ ∝ a⁻⁴
      assert(std::fabs(sp->P(bg.data()) - sp->Rho(bg.data()) / 3.) <
             1e-10 * sp->Rho(bg.data()));  // massless EoS
    }
  }

  // ── Statistics guard: statistics key must be fermion|boson ─────────────────
  {
    background pba{};
    pba.H0            = 2.2e-4;
    NcdmSettings sett = TestSettings();

    auto make_ctx = [&](FileContent& fc) {
      SpeciesBuildContext ctx{};
      ctx.pfc           = &fc;
      ctx.pba           = &pba;
      ctx.ncdm_settings = &sett;
      ctx.bgm           = nullptr;
      return ctx;
    };

    {  // valid fermion
      FileContent fc = BaseFc();
      fc.set("drx.statistics", "fermion");
      auto ctx = make_ctx(fc);
      auto vec = DrPsdSpecies::CreateAll(ctx);
      assert(vec.size() == 1);
    }
    {  // valid boson
      FileContent fc = BaseFc();
      fc.set("drx.statistics", "boson");
      auto ctx = make_ctx(fc);
      auto vec = DrPsdSpecies::CreateAll(ctx);
      assert(vec.size() == 1);
    }
    {  // garbage statistics -> class_test_severe (invalid_argument)
      FileContent fc = BaseFc();
      fc.set("drx.statistics", "tachyon");
      auto ctx   = make_ctx(fc);
      bool threw = false;
      try {
        DrPsdSpecies::CreateAll(ctx);
      }
      catch (const std::invalid_argument&) {
        threw = true;
      }
      assert(threw);
    }
  }

  // ── CopyPerturbationsAcrossSwitch regression (#372 class, Task 2.2) ─────────
  // The per-bin layout must migrate slot-by-slot when the pv is reallocated at
  // another species' approximation switch (inherits the NCDM-family base copy).
  {
    FileContent fc  = BaseFc();
    auto sp         = MakeDrPsd(fc, Statistics::Fermion);
    const int qn    = sp->q_size();
    const int l_max = 5;
    assert(qn > 0);

    NCDMBaseSpecies::PerturbLayout old_l, new_l;
    FillNcdmLayout(old_l, qn, l_max, /*idx0=*/7);
    FillNcdmLayout(new_l, qn, l_max, /*idx0=*/3);

    std::vector<double> old_y(7 + qn * (l_max + 1), 0.);
    std::vector<double> new_y(7 + qn * (l_max + 1), 0.);
    for (int iq = 0; iq < qn; ++iq)
      for (int l = 0; l <= l_max; ++l)
        old_y[old_l.index_per_q[iq] + l] = Marker(iq, l);

    const PerturbSwitchContext ctx{};
    sp->CopyPerturbationsAcrossSwitch(old_l, new_l, old_y.data(), new_y.data(), ctx);

    for (int iq = 0; iq < qn; ++iq)
      for (int l = 0; l <= l_max; ++l)
        assert(new_y[new_l.index_per_q[iq] + l] == Marker(iq, l));
    // Slots outside the species' layout stay untouched (still zero).
    for (int i = 0; i < 3; ++i)
      assert(new_y[i] == 0.);
    // Uninitialized layouts (shape -1) must be a safe no-op.
    NCDMBaseSpecies::PerturbLayout empty_old, empty_new;
    sp->CopyPerturbationsAcrossSwitch(empty_old, empty_new, old_y.data(), new_y.data(), ctx);
  }

  // ── Column-name disambiguation ──────────────────────────────────────────────
  // The composite constructs fermion and boson daughters under the SAME
  // instance name (they are not separately dot-addressable) -- their
  // background/output column names must not collide.
  {
    FileContent fc = BaseFc();
    auto sp_l      = MakeDrPsd(fc, Statistics::Fermion);
    auto sp_phi    = MakeDrPsd(fc, Statistics::Boson);

    std::string titles;
    BackgroundColumnWriter w(titles);
    sp_l->WriteBackgroundColumnTitles(w);
    sp_phi->WriteBackgroundColumnTitles(w);

    // No duplicate column name in the combined title string.
    std::vector<std::string> cols;
    size_t start = 0;
    for (size_t i = 0; i < titles.size(); ++i) {
      if (titles[i] == _DELIMITER_[0]) {
        cols.push_back(titles.substr(start, i - start));
        start = i + 1;
      }
    }
    for (size_t i = 0; i < cols.size(); ++i)
      for (size_t j = i + 1; j < cols.size(); ++j)
        assert(cols[i] != cols[j]);
    // Both statistics-specific density columns are present.
    assert(titles.find("rho_dr_drx_l") != std::string::npos);
    assert(titles.find("rho_dr_drx_phi") != std::string::npos);
  }

  // ── Absolute normalization (#385) ────────────────────────────────────────────
  // A standalone fermion daughter thermal at the standard neutrino temperature
  // T_nu/T_cmb = (4/11)^(1/3), deg left at its class default (physical dof=2,
  // "particle+antiparticle" -- see dr_psd_species.cpp:37), must reproduce
  // rho_nu_rel() -- CLASS's own closed-form reference density for ONE standard
  // massless neutrino -- to quadrature tolerance. Pre-fix (deg=2 with no
  // 1/(2*pi)^3) this is off by a large, constant, state-independent factor.
  {
    FileContent fc = BaseFc();
    fc.set("drx.dr_T", std::to_string(std::pow(4. / 11., 1. / 3.)));
    fc.set("drx.dr_q_min", "1e-3");
    fc.set("drx.dr_q_max", "50");
    fc.set("drx.dr_N_q", "400");
    auto sp     = MakeDrPsd(fc, Statistics::Fermion);
    const int N = sp->q_size();

    int index_bg = 0, index_bi = 0;
    sp->RegisterBackgroundIndices(index_bg);
    sp->RegisterIntegrationIndices(index_bi);
    std::vector<double> bi(index_bi, 0.), bg(index_bg, 0.);
    BackgroundICContext icc{};
    icc.pvecback_integration = bi.data();
    sp->SetBackgroundInitialConditions(icc);
    double max_f = 0.;
    for (int i = 0; i < N; ++i)
      max_f = std::max(max_f, bi[sp->bi_f_index() + i]);
    assert(max_f > 1e-3);  // thermal seed actually triggered (f_ini_ fallback == 1)

    sp->ComputeBackground(1.0, bi.data(), bg.data());
    const double rho = sp->Rho(bg.data());
    const double rel = std::fabs(rho - sp->rho_nu_rel()) / sp->rho_nu_rel();
    std::printf(
        "  absolute normalization: rho=%.6e  rho_nu_rel=%.6e  rel=%.3e  ratio=%.3f  (need "
        "rel < 1e-3)\n",
        rho,
        sp->rho_nu_rel(),
        rel,
        rho / sp->rho_nu_rel());
    assert(rel < 1e-3);
  }

  std::printf("dr_psd test passed\n");
  return 0;
}
