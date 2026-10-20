/**
 * ppw->s_l must be sized from the species that actually own hierarchies, not
 * from a hard-coded list of species keys.
 *
 * Regression test for #421: a DNCDM composite owns a DarkRadiationSpecies
 * daughter running to l_max_dr, but is keyed by its instance name rather than
 * "DCDM_DR", so the key test that used to size s_l missed it and every
 * free-streaming loop read past the end of the array.
 */

#include <cstdio>
#include <memory>

#include "background.h"
#include "composite_species.h"
#include "dark_radiation_species.h"
#include "dcdm_dr_species.h"
#include "dncdm_dr_species.h"
#include "dncdm_species.h"
#include "parser.h"
#include "photons.h"
#include "precision.h"
#include "ultra_relativistic.h"

namespace {

int failures = 0;

void expect(bool ok, const char* what) {
  if (!ok) {
    std::fprintf(stderr, "FAIL: %s\n", what);
    ++failures;
  }
}

/* The invariant that matters: whatever a species' hierarchy will index, s_l must
   be at least that long. The free-streaming loops read s_l[l+1] for l < l_max, so
   the highest index touched is s_l[l_max] and s_l needs l_max + 1 entries. */
void test_dark_radiation_reports_l_max_dr() {
  precision ppr;
  ppr.l_max_dr = 80;
  background pba{};
  DarkRadiationSpecies dr("DR", &pba, /*bgm=*/nullptr);
  expect(dr.MaxMultipole(&ppr, /*tensors=*/false) == 80,
         "DarkRadiationSpecies must report l_max_dr for scalars");
  expect(dr.MaxMultipole(&ppr, /*tensors=*/true) == 0,
         "DarkRadiationSpecies carries no tensor hierarchy");
}

/* #421 itself: the composite must forward its children's requirement. Before the
   fix this species contributed nothing, s_l was sized 18, and the daughter
   indexed s_l[81]. */
void test_dncdm_composite_forwards_its_dr_daughter() {
  precision ppr;
  ppr.l_max_dr   = 80;
  ppr.l_max_ncdm = 10;
  background pba{};

  FileContent fc;
  fc.set("dncdm1.type", "ncdm_decay_dr");
  fc.set("dncdm1.m", "1.0");
  fc.set("dncdm1.Gamma", "1e3");
  fc.set("dncdm1.Omega_dncdmdr", "0.001");
  NcdmSettings settings{};
  settings.h           = 0.67556;
  settings.T_cmb       = 2.7255;
  settings.tol_ncdm    = 1.e-3;
  settings.tol_ncdm_bg = 1.e-5;
  settings.tol_M_ncdm  = 1.e-5;

  DNCDM_DR_Species
      dncdm_dr(std::make_unique<DNCDMSpecies>(&fc, "dncdm1", settings, &pba, /*bgm=*/nullptr),
               &pba,
               /*bgm=*/nullptr);

  expect(dncdm_dr.MaxMultipole(&ppr, /*tensors=*/false) >= ppr.l_max_dr,
         "DNCDM_DR_Species must forward its DarkRadiationSpecies daughter's l_max_dr");
}

void test_dcdm_dr_composite_forwards_too() {
  precision ppr;
  ppr.l_max_dr = 80;
  background pba{};
  DCDM_DR_Species dcdm_dr(&pba,
                          /*bgm=*/nullptr,
                          /*omega0_dcdmdr=*/0.05,
                          /*Gamma_dcdm=*/100.,
                          /*Omega_ini_dcdm=*/0.);
  expect(dcdm_dr.MaxMultipole(&ppr, /*tensors=*/false) >= ppr.l_max_dr,
         "DCDM_DR_Species must forward its DarkRadiationSpecies daughter's l_max_dr");
}

/* Photons are the only species whose requirement differs between modes, and the
   tensor hierarchy is the shorter one -- so a mode-blind implementation would
   over-allocate silently rather than fail. Pin both. */
void test_photons_are_mode_dependent() {
  precision ppr;
  ppr.l_max_g         = 12;
  ppr.l_max_pol_g     = 10;
  ppr.l_max_g_ten     = 5;
  ppr.l_max_pol_g_ten = 7;
  background pba{};
  PhotonsSpecies g(pba);
  expect(g.MaxMultipole(&ppr, /*tensors=*/false) == 12, "photons: scalars use l_max_g");
  expect(g.MaxMultipole(&ppr, /*tensors=*/true) == 7, "photons: tensors use l_max_pol_g_ten");
}

void test_ultra_relativistic_reports_l_max_ur() {
  precision ppr;
  ppr.l_max_ur = 33;
  background pba{};
  UltraRelativisticSpecies ur(pba, /*omega0_ur=*/1.7e-5);
  expect(ur.MaxMultipole(&ppr, /*tensors=*/false) == 33, "ur: scalars use l_max_ur");
  expect(ur.MaxMultipole(&ppr, /*tensors=*/true) == 33, "ur: tensors use l_max_ur too");
}

}  // namespace

int main() {
  test_dark_radiation_reports_l_max_dr();
  test_dncdm_composite_forwards_its_dr_daughter();
  test_dcdm_dr_composite_forwards_too();
  test_photons_are_mode_dependent();
  test_ultra_relativistic_reports_l_max_ur();

  if (failures != 0) {
    std::fprintf(stderr, "max_multipole_test: %d failure(s)\n", failures);
    return 1;
  }
  std::printf("max_multipole_test: all tests passed\n");
  return 0;
}
