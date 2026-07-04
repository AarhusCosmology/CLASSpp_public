#include <cstdio>
#include <memory>
#include <vector>

#include "background.h"
#include "cdm.h"
#include "composite_species.h"
#include "dark_radiation_species.h"
#include "dcdm.h"
#include "dcdm_dr_species.h"
#include "dncdm_dr_species.h"
#include "dncdm_species.h"
#include "idm_dr.h"
#include "idm_dr_idr_species.h"
#include "idm_drmd.h"
#include "idm_drmd_idr_drmd_species.h"
#include "idr.h"
#include "idr_drmd.h"
#include "ncdm_base_species.h"
#include "parser.h"
#include "scalar_field.h"
#include "type3_species.h"

namespace {

// Every composite's CreatePerturbLayout must return a layout deriving from
// CompositeSpecies::PerturbLayout with one correctly-typed sub-layout per child,
// in children_ order. The generic child-loops (RegisterPerturbationIndices,
// StressEnergy, DelegateTally, ...) static_cast to that type and index
// child_layouts, so a composite whose layout lacks child_layouts is UB (#358).
template <typename Child0Layout, typename Child1Layout>
int CheckCompositeLayout(const CompositeSpecies& sp, const char* what) {
  auto layout = sp.CreatePerturbLayout();
  auto* comp  = dynamic_cast<CompositeSpecies::PerturbLayout*>(layout.get());
  if (comp == nullptr) {
    std::fprintf(stderr,
                 "FAIL: %s layout does not derive from CompositeSpecies::PerturbLayout\n",
                 what);
    return 1;
  }
  if (comp->child_layouts.size() != 2) {
    std::fprintf(stderr,
                 "FAIL: %s expected 2 child_layouts, got %zu\n",
                 what,
                 comp->child_layouts.size());
    return 1;
  }
  if (dynamic_cast<Child0Layout*>(comp->child_layouts[0].get()) == nullptr) {
    std::fprintf(stderr, "FAIL: %s child_layouts[0] has the wrong dynamic type\n", what);
    return 1;
  }
  if (dynamic_cast<Child1Layout*>(comp->child_layouts[1].get()) == nullptr) {
    std::fprintf(stderr, "FAIL: %s child_layouts[1] has the wrong dynamic type\n", what);
    return 1;
  }
  return 0;
}

}  // namespace

int main() {
  background pba{};
  pba.H0                           = 1e-4;
  const std::vector<double> params = {1.22, 0.0, 0.0, 0.0};
  auto scf = std::make_unique<ScalarFieldSpecies>(pba,
                                                  /*omega0=*/0.7,
                                                  params,
                                                  /*tuning=*/0,
                                                  /*attractor=*/true,
                                                  /*phi_ini=*/1.,
                                                  /*phi_prime_ini=*/1.,
                                                  DefaultScalarFieldPotential(),
                                                  /*beta=*/-2.0);
  Type3Species t3(pba, /*omega0_cdm=*/0.12, std::move(scf));

  int rc  = 0;
  rc     |= CheckCompositeLayout<CDMSpecies::PerturbLayout,
                                 ScalarFieldSpecies::PerturbLayout>(t3, "Type3Species");

  IDM_DR_IDR_Species idm_dr_idr(pba,
                                /*omega0_idm_dr=*/0.05,
                                /*omega0_idr=*/0.02,
                                /*T_idr=*/1.,
                                /*l_max_idr=*/17,
                                /*a_idm_dr=*/1.e3,
                                /*nindex_idm_dr=*/4.,
                                /*m_idm=*/1.e11,
                                /*b_idr=*/0.,
                                /*idr_nature=*/0,
                                /*alpha_idm_dr=*/{1.5},
                                /*beta_idr=*/{1.5});
  rc |= CheckCompositeLayout<IDM_DRSpecies::PerturbLayout,
                             IDRSpecies::PerturbLayout>(idm_dr_idr, "IDM_DR_IDR_Species");

  IDM_DRMD_IDR_DRMD_Species idm_drmd(pba,
                                     /*omega0_idm_drmd=*/0.05,
                                     /*omega0_idr_drmd=*/0.02,
                                     /*f_idm_drmd=*/0.1,
                                     /*G_over_aH_drmd=*/1.,
                                     /*delta_Neff_drmd=*/0.,
                                     /*z_stop=*/1.e4);
  rc |= CheckCompositeLayout<IDM_DRMDSpecies::PerturbLayout,
                             IDR_DRMDSpecies::PerturbLayout>(idm_drmd, "IDM_DRMD_IDR_DRMD_Species");

  DCDM_DR_Species dcdm_dr(&pba,
                          /*bgm=*/nullptr,
                          /*omega0_dcdmdr=*/0.05,
                          /*Gamma_dcdm=*/100.,
                          /*Omega_ini_dcdm=*/0.);
  rc |= CheckCompositeLayout<DCDMSpecies::PerturbLayout,
                             DarkRadiationSpecies::PerturbLayout>(dcdm_dr, "DCDM_DR_Species");

  {
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
    rc |= CheckCompositeLayout<NCDMBaseSpecies::PerturbLayout,
                               DarkRadiationSpecies::PerturbLayout>(dncdm_dr, "DNCDM_DR_Species");
  }

  if (rc != 0)
    return 1;

  std::printf("composite layout test passed\n");
  return 0;
}
