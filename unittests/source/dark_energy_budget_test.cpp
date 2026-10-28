// The background's dark energy (Omega_de, w_de), which halofit reads: the species of type
// DarkEnergy, and any other species while its pressure is negative. A species of type Other is
// split into radiation (3p) and matter (rho - 3p) only while its pressure is positive.
//
// Checked against the per-species columns of the background table: Lambda alone (w = -1), a
// CPL fluid (w = w_fld(a), not w_fld today), Lambda given with a fluid that closes the budget
// (both present, the density-weighted mean), and a quintessence field with Lambda, which
// counts as dark energy once its pressure turns negative and stays out of Omega_m then. Until 2026-10-27 Omega_m(z = 0) held rho - 3p of that field:
// 0.486 instead of 0.312 for Omega_scf = 0.1 at w = -0.25.
#include <cassert>
#include <cmath>
#include <cstdio>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "background.h"
#include "background_module.h"
#include "cosmology.h"
#include "parser.h"

namespace {

struct Table {
  std::map<std::string, int> col;
  std::vector<double> data;
  int n_titles = 0;
  int n_rows   = 0;
  double at(int r, const std::string& name) const {
    auto it = col.find(name);
    if (it == col.end()) {
      std::printf("column '%s' not found\n", name.c_str());
      assert(false);
    }
    return data[static_cast<size_t>(r) * n_titles + it->second];
  }
  double get(int r, const std::string& name) const {  // 0 for a species that is absent
    return col.count(name) ? at(r, name) : 0.;
  }
};

struct Check {
  double w_err_max      = 0.;  // max |w_de - sum p / sum rho| over the table
  double omega_err_max  = 0.;  // max |Omega_de - sum rho / rho_crit|
  double closure_maxdev = 0.;  // max |Omega_m + Omega_r + Omega_de - 1| (flat models)
  double w_today        = 0.;
  double omega_m_today  = 0.;
  double matter_today   = 0.;  // (rho_cdm + rho_b + non-relativistic ncdm) / rho_crit
  double lambda_today   = 0.;  // rho_lambda / rho_crit
  double fld_today      = 0.;  // rho_fld / rho_crit
  int rows_with_scf_de  = 0;   // rows where the field's pressure is negative
  int rows_with_scf_pos = 0;   // and positive
};

Check Run(const std::map<std::string, std::string>& extra) {
  FileContent fc;
  fc.set("h", "0.6781");
  fc.set("omega_b", "0.02238");
  fc.set("omega_cdm", "0.1201");
  fc.set("output", "");  // background only
  for (const auto& [k, v] : extra) {
    fc.set(k, v);
  }
  Cosmology cosmo(fc);
  auto& bgm = cosmo.GetBackgroundModule();

  Table t;
  std::string titles;
  bgm->background_output_titles(titles);
  {
    std::stringstream ss(titles);
    std::string s;
    while (std::getline(ss, s, '\t')) {
      if (!s.empty()) {
        t.col[s] = t.n_titles++;
      }
    }
  }
  t.n_rows = bgm->bt_size_;
  t.data.resize(static_cast<size_t>(t.n_titles) * t.n_rows);
  bgm->background_output_data(t.n_titles, t.data.data());

  Check out;
  for (int r = 0; r < t.n_rows; ++r) {
    /* the stored row behind this line of the output table */
    const double* pvecback = bgm->background_table_.data() + static_cast<size_t>(r) * bgm->bg_size_;
    const double rho_l     = t.get(r, "(.)rho_lambda");
    const double rho_f     = t.get(r, "(.)rho_fld");
    const double w_f       = t.col.count("(.)w_fld") ? t.at(r, "(.)w_fld") : 0.;
    const double rho_s     = t.get(r, "(.)rho_scf");
    const double p_s       = t.get(r, "(.)p_scf");
    double rho_de          = rho_l + rho_f;
    double p_de            = -rho_l + w_f * rho_f;
    if (p_s < 0.) {
      rho_de += rho_s;
      p_de   += p_s;
      ++out.rows_with_scf_de;
    }
    else if (rho_s > 0.) {
      ++out.rows_with_scf_pos;
    }
    const double rho_crit = t.at(r, "(.)rho_crit");
    const double Omega_de = pvecback[bgm->index_bg_Omega_de_];
    const double w_de     = pvecback[bgm->index_bg_w_de_];
    if (rho_de > 1e-6 * rho_crit) {  // where the dark energy is not negligible
      out.w_err_max = std::max(out.w_err_max, std::fabs(w_de - p_de / rho_de));
    }
    out.omega_err_max  = std::max(out.omega_err_max, std::fabs(Omega_de - rho_de / rho_crit));
    out.closure_maxdev = std::max(out.closure_maxdev,
                                  std::fabs(pvecback[bgm->index_bg_Omega_m_] +
                                            pvecback[bgm->index_bg_Omega_r_] + Omega_de - 1.));
    if (t.at(r, "z") == 0.) {
      out.w_today       = w_de;
      out.omega_m_today = pvecback[bgm->index_bg_Omega_m_];
      double rho_nr     = t.at(r, "(.)rho_cdm") + t.at(r, "(.)rho_b");
      out.matter_today  = rho_nr / rho_crit;
      out.lambda_today  = rho_l / rho_crit;
      out.fld_today     = rho_f / rho_crit;
    }
  }
  return out;
}

void Report(const char* label, const Check& c) {
  std::printf(
      "%-34s w_de err %.1e | Omega_de err %.1e | Omega_m+r+de-1 %.1e | w_de(0) %+.4f | "
      "Omega_m(0) %.4f (cdm+b %.4f) | field rows p<0 %d, p>0 %d\n",
      label,
      c.w_err_max,
      c.omega_err_max,
      c.closure_maxdev,
      c.w_today,
      c.omega_m_today,
      c.matter_today,
      c.rows_with_scf_de,
      c.rows_with_scf_pos);
}

}  // namespace

int main() {
  /* Lambda alone: w = -1 exactly, and Omega_m + Omega_r + Omega_de = 1 (flat, no ncdm). */
  {
    const Check c = Run({});
    Report("LCDM", c);
    assert(c.w_today == -1.);
    assert(c.w_err_max == 0.);
    assert(c.omega_err_max < 1e-12);
    assert(c.closure_maxdev < 1e-10);
  }

  /* A CPL fluid: w_de follows w_fld(a), not its value today (-0.45). */
  {
    const Check c = Run({{"Omega_Lambda", "0"},
                         {"w0_fld", "-0.45"},
                         {"wa_fld", "-1.6"},
                         {"cs2_fld", "1"},
                         {"use_ppf", "yes"}});
    Report("CPL fluid (-0.45, -1.6)", c);
    assert(std::fabs(c.w_today + 0.45) < 1e-10);
    assert(c.w_err_max < 1e-12);
    assert(c.omega_err_max < 1e-12);
    assert(c.closure_maxdev < 1e-10);
  }

  /* Lambda given with a fluid that closes the budget: both are there (until 2026-10-27 Lambda
     was built only as the closure species, and this gave a fluid of 0.69 and no Lambda), and
     w_de is the density-weighted mean of -1 and w_fld. */
  {
    const Check c = Run({{"Omega_Lambda", "0.3"},
                         {"w0_fld", "-0.8"},
                         {"wa_fld", "0.2"},
                         {"cs2_fld", "1"},
                         {"use_ppf", "no"}});
    Report("Lambda 0.3 + fluid (-0.8, 0.2)", c);
    assert(std::fabs(c.lambda_today - 0.3) < 1e-10);
    assert(c.fld_today > 0.38 && c.fld_today < 0.40);  // 1 - 0.3 - Omega_m - Omega_r
    assert(std::fabs(c.w_today - (-0.3 - 0.8 * c.fld_today) / (0.3 + c.fld_today)) < 1e-10);
    assert(c.w_err_max < 1e-12);
    assert(c.omega_err_max < 1e-12);
    assert(c.closure_maxdev < 1e-10);
  }

  /* A quintessence field with Lambda (the lvl2 scalar-field scenario): dark energy once its
     pressure is negative, matter and radiation before, and out of Omega_m today; w_de is the
     density-weighted mean of the two. */
  {
    const Check c = Run(
        {{"Omega_scf", "0.1"}, {"attractor_ic_scf", "yes"}, {"scf_parameters", "10, 0, 0, 0"}});
    Report("quintessence (Omega_scf 0.1) + Lambda", c);
    assert(c.rows_with_scf_de > 0 && c.rows_with_scf_pos > 0);
    assert(c.w_today > -0.95 && c.w_today < -0.85);  // (-0.59 + 0.1 x -0.245) / 0.69
    assert(c.w_err_max < 1e-10);
    assert(c.omega_err_max < 1e-10);
    assert(c.closure_maxdev < 1e-8);
    assert(std::fabs(c.omega_m_today - c.matter_today) < 1e-10);
  }

  std::printf("dark energy budget: all checks passed\n");
  return 0;
}
