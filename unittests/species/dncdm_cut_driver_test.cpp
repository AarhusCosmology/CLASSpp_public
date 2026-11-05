// #496: the cut of a trapezoid grid (quadrature_strategy = 3) on the metric drivers of a
// decaying ncdm in collision-owned (F-space) mode. The grid integrates the parent's PSD
// truncated at q_max, whose d ln f/d ln q carries -f(q_max) delta(q - q_max): the last bin
// gets an extra driver -q_max/dq_last = -2N times its occupation f-bar, every other bin none.
// With y = 0 the derivatives are the drivers alone, so each bin's l = 0, 1, 2 can be checked
// against the formula, with ncdm_driver_normalisation off. With it on, the drivers' q-sum
// must be -3 (rho+p) on the grid (#497), which a spline of ln f on 20 bins misses. (Decay-only mode is covered end to end by the gauge test in
// python/test_class.py, since its Newtonian block needs a perturbations module.)

#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

#include "background.h"
#include "dncdm_species.h"
#include "ncdm_base_species.h"
#include "parser.h"
#include "perturbations.h"
#include "perturbations_module.h"

namespace {

std::unique_ptr<DNCDMSpecies> MakeTrapezoidDncdm(const background* pba, int bins, bool normalise) {
  FileContent fc;
  fc.set("dncdm1.type", "ncdm_decay_dr");
  fc.set("dncdm1.m", "1.0");
  fc.set("dncdm1.Gamma", "1e3");
  fc.set("dncdm1.Omega_dncdmdr", "0.001");
  fc.set("dncdm1.quadrature_strategy", "3");  // qm_trapz: stops at max_q = 15
  fc.set("dncdm1.momenta_bins", std::to_string(bins));
  NcdmSettings settings{};
  settings.h                 = 0.67556;
  settings.T_cmb             = 2.7255;
  settings.tol_ncdm          = 1.e-3;
  settings.tol_ncdm_bg       = 1.e-5;
  settings.tol_M_ncdm        = 1.e-5;
  settings.normalise_drivers = normalise;
  return std::make_unique<DNCDMSpecies>(&fc, "dncdm1", settings, pba, /*bgm=*/nullptr);
}

bool Close(double got, double expect) {
  return std::fabs(got - expect) <= 1e-12 * std::fabs(expect) + 1e-300;
}

}  // namespace

int Run(bool normalise) {
  background pba{};
  pba.H0 = 2.2e-4;

  const int bins = 20;
  auto sp        = MakeTrapezoidDncdm(&pba, bins, normalise);
  sp->SetCollisionOwned(true);
  const int N                   = sp->q_size();
  const std::vector<double>& q  = sp->GetQ();
  const std::vector<double>& dq = sp->dq();
  int rc                        = 0;

  // The cut of a trapezoid grid: -q_max/dq_last, with dq_last = h/2 and q_max = N h.
  const double cut = -q.back() / dq.back();
  if (N != bins || !Close(cut, -2. * N)) {
    std::fprintf(stderr, "FAIL: N = %d, cut = %g, expected -2N = %d\n", N, cut, -2 * N);
    rc = 1;
  }

  int index_bg = 0;
  sp->RegisterBackgroundIndices(index_bg);
  std::vector<double> pvecback(index_bg, 0.);
  for (int iq = 0; iq < N; ++iq) {
    pvecback[sp->bg_lnf_index() + iq]      = -q[iq] - 0.3;  // a decayed occupation
    pvecback[sp->bg_dlnfdlnq_index() + iq] = -q[iq];        // its smooth d ln f/d ln q
  }

  const int l_max = 6;
  NCDMBaseSpecies::PerturbLayout layout;
  layout.q_size = N;
  layout.l_max  = l_max;
  for (int iq = 0; iq < N; ++iq) {
    layout.index_per_q.push_back(iq * (l_max + 1));
  }
  std::vector<double> y(N * (l_max + 1), 0.);
  std::vector<double> dy(N * (l_max + 1), 0.);

  perturb_workspace ppw;
  ppw.pvecback                     = pvecback;
  ppw.s_l                          = std::vector<double>(l_max + 2, 1.);
  ppw.scalar_ctx.k                 = 0.01;
  ppw.scalar_ctx.a                 = 1e-2;
  ppw.scalar_ctx.a2                = 1e-4;
  ppw.scalar_ctx.metric_continuity = 0.3;  // all three nonzero: each driver is checked
  ppw.scalar_ctx.metric_euler      = 2e-5;
  ppw.scalar_ctx.metric_shear      = 0.2;
  perturb_parameters_and_workspace ppaw(nullptr);  // collision-owned mode needs no module
  ppaw.ppw = &ppw;
  ppaw.k   = ppw.scalar_ctx.k;

  sp->PerturbDerivs(layout, 1., y.data(), dy.data(), ppaw);

  const double k = ppw.scalar_ctx.k;
  const double M = sp->GetMass();
  if (normalise) {
    // sum_q dq q^2 eps (driver of bin q) = -3 sum_q dq f q^2 (eps + q^2 / 3 eps), the l = 0
    // driver being metric_continuity * driver / 3.
    double driver_sum = 0.;
    double bg_sum     = 0.;
    for (int iq = 0; iq < N; ++iq) {
      const double f    = std::exp(pvecback[sp->bg_lnf_index() + iq]);
      const double eps  = std::sqrt(q[iq] * q[iq] + ppw.scalar_ctx.a2 * M * M);
      driver_sum       += dq[iq] * q[iq] * q[iq] * eps * 3. * dy[layout.index_per_q[iq]] /
                          ppw.scalar_ctx.metric_continuity;
      bg_sum           += dq[iq] * f * q[iq] * q[iq] * (eps + q[iq] * q[iq] / (3. * eps));
    }
    if (!Close(driver_sum, -3. * bg_sum)) {
      std::fprintf(stderr,
                   "FAIL: normalised driver sum %.17g, expected -3(rho+p) = %.17g\n",
                   driver_sum,
                   -3. * bg_sum);
      rc = 1;
    }
    else {
      std::printf("dncdm drivers normalised: q-sum = -3(rho+p) on the grid: OK\n");
    }
    return rc;
  }
  for (int iq = 0; iq < N; ++iq) {
    const double f         = std::exp(pvecback[sp->bg_lnf_index() + iq]);
    const double dlnf      = pvecback[sp->bg_dlnfdlnq_index() + iq] + (iq == N - 1 ? cut : 0.);
    const double driver    = dlnf * f;
    const double epsilon   = std::sqrt(q[iq] * q[iq] + ppw.scalar_ctx.a2 * M * M);
    const int idx          = layout.index_per_q[iq];
    const double expect[3] = {ppw.scalar_ctx.metric_continuity * driver / 3.,
                              -epsilon * ppw.scalar_ctx.metric_euler / (3. * q[iq] * k) * driver,
                              -ppw.s_l[2] * ppw.scalar_ctx.metric_shear * 2. / 15. * driver};
    for (int l = 0; l < 3; ++l) {
      if (!Close(dy[idx + l], expect[l])) {
        std::fprintf(stderr,
                     "FAIL: bin %d, l = %d: dy = %.17g, expected %.17g\n",
                     iq,
                     l,
                     dy[idx + l],
                     expect[l]);
        rc = 1;
      }
    }
  }
  if (rc == 0) {
    std::printf("dncdm cut driver (collision-owned): OK, cut = %g on bin %d of %d\n",
                cut,
                N - 1,
                N);
  }
  return rc;
}

int main() {
  return Run(false) | Run(true);
}
