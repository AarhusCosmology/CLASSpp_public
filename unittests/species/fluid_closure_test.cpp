// A CPL fluid must close the background at a = 1: H(z = 0) is the H0 that was asked for, and
// rho_fld(a) is rho_fld,0 exp(integral of 3(1+w) dln a), whatever the background tolerance.
//
// rho_fld used to be integrated as an ODE from a_ini. For w0 + wa << -1 the density climbs by
// ~100 decades before today and the integrator's relative error lands in H(z = 0): at the
// default tol_background_integration = 1e-7, |H(0)/H0 - 1| was 2.6e-4 at (w0, wa) =
// (-0.58, -1.30), near DESI DR2's CMB+DESI best fit, and 6.2e-2 at (0.16, -3.64).
#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

#include "background.h"
#include "background_module.h"
#include "constants.h"
#include "cosmology.h"
#include "parser.h"

namespace {

struct Closure {
  double h_ratio_minus_one;  // H(z = 0) / H0 - 1
  double rho_rel_err_max;    // max over the table of |rho_fld / closed form - 1|
};

Closure Run(double w0, double wa, const char* use_ppf) {
  const double h = 0.6487;
  FileContent fc;
  fc.set("h", std::to_string(h));
  fc.set("omega_b", "0.02223");
  fc.set("omega_cdm", "0.1191");
  fc.set("Omega_Lambda", "0");
  fc.set("w0_fld", std::to_string(w0));
  fc.set("wa_fld", std::to_string(wa));
  fc.set("cs2_fld", "1");
  fc.set("use_ppf", use_ppf);
  fc.set("output", "");  // background only
  Cosmology cosmo(fc);
  auto& bgm = cosmo.GetBackgroundModule();

  std::string titles;
  bgm->background_output_titles(titles);
  std::vector<std::string> cols;
  {
    std::stringstream ss(titles);
    std::string t;
    while (std::getline(ss, t, '\t')) {
      if (!t.empty()) {
        cols.push_back(t);
      }
    }
  }
  auto col = [&](const std::string& nm) {
    for (size_t i = 0; i < cols.size(); ++i) {
      if (cols[i] == nm) {
        return static_cast<int>(i);
      }
    }
    std::printf("column '%s' not found\n", nm.c_str());
    assert(false);
    return -1;
  };
  const int n_titles = static_cast<int>(cols.size());
  const int n_rows   = bgm->bt_size_;
  std::vector<double> data(static_cast<size_t>(n_titles) * n_rows);
  bgm->background_output_data(n_titles, data.data());
  auto row = [&](int r, int c) { return data[static_cast<size_t>(r) * n_titles + c]; };

  const int ic_z    = col("z");
  const int ic_H    = col("H [1/Mpc]");
  const int ic_rho  = col("(.)rho_fld");
  const int ic_crit = col("(.)rho_crit");

  const double H0 = h * 1.e5 / _c_;
  Closure out{};
  for (int r = 0; r < n_rows; ++r) {
    if (row(r, ic_z) == 0.) {
      out.h_ratio_minus_one = row(r, ic_H) / H0 - 1.;
    }
  }
  // Closed form for CPL: rho(a)/rho0 = a^{-3(1+w0+wa)} exp(3 wa (a - 1)). Omega_fld,0 is
  // read off the table's own last row so the check is about the a-dependence alone.
  const double rho0 = row(n_rows - 1, ic_rho);
  assert(std::fabs(row(n_rows - 1, ic_z)) < 1e-12);
  for (int r = 0; r < n_rows; ++r) {
    const double a      = 1. / (1. + row(r, ic_z));
    const double closed = rho0 * std::pow(a, -3. * (1. + w0 + wa)) * std::exp(3. * wa * (a - 1.));
    if (closed > 1e-30 * row(r, ic_crit)) {  // skip the epochs where the fluid is negligible
      out.rho_rel_err_max = std::fmax(out.rho_rel_err_max, std::fabs(row(r, ic_rho) / closed - 1.));
    }
  }
  return out;
}

}  // namespace

int main() {
  struct Case {
    double w0, wa;
  };
  const std::array<Case, 3> cases = {{{-0.578, -1.298}, {0.159, -3.640}, {-0.9, 0.3}}};
  for (const char* ppf : {"yes", "no"}) {
    for (const auto& c : cases) {
      const Closure r = Run(c.w0, c.wa, ppf);
      std::printf(
          "w0 = %+.3f wa = %+.3f ppf = %-3s  H(0)/H0 - 1 = %+.2e  max |rho/closed - 1| = %.2e\n",
          c.w0,
          c.wa,
          ppf,
          r.h_ratio_minus_one,
          r.rho_rel_err_max);
      assert(std::fabs(r.h_ratio_minus_one) < 1e-10);
      assert(r.rho_rel_err_max < 1e-10);
    }
  }
  std::printf("fluid closure: OK\n");
  return 0;
}
