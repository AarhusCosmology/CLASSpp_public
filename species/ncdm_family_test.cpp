// Unit test for the NCDM-family helpers (issues #376 and #378):
// SynthesiseNcdmFluidApproximation must gather <name>.fluid_approximation from
// every family type that consults index_ap_ncdmfa (not just ncdm_standard),
// reject the key on species that ignore it, and CreateAllNcdmInstances<T> must
// reproduce the shared per-instance factory loop.
#include "species/ncdm_family.h"

#include <cassert>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

#include "background.h"
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

static bool SynthesisThrows(FileContent& fc) {
  try {
    SynthesiseNcdmFluidApproximation(&fc);
    return false;
  }
  catch (const std::invalid_argument&) {
    return true;
  }
}

static bool Contains(const std::vector<std::string>& v, const std::string& s) {
  for (const auto& e : v) {
    if (e == s)
      return true;
  }
  return false;
}

int main() {
  // ── #376 repro: an axion-only instance must feed the family-wide key ──
  {
    FileContent fc;
    fc.set("ax.type", "ncdm_axion");
    fc.set("ax.fluid_approximation", "3");
    SynthesiseNcdmFluidApproximation(&fc);
    auto v = fc.get<std::string>("ncdm_fluid_approximation");
    assert(v && *v == "3");
  }

  // ── Every consumer type contributes (greybody, self-interacting, decay_dr) ──
  {
    FileContent fc;
    fc.set("gb.type", "ncdm_greybody");
    fc.set("gb.fluid_approximation", "2");
    fc.set("si.type", "ncdm_self_interacting");
    fc.set("si.fluid_approximation", "2");
    fc.set("dn.type", "ncdm_decay_dr");
    fc.set("dn.fluid_approximation", "2");
    SynthesiseNcdmFluidApproximation(&fc);
    auto v = fc.get<std::string>("ncdm_fluid_approximation");
    assert(v && *v == "2");
  }

  // ── Identity guard now spans the family: standard vs axion disagreeing ──
  {
    FileContent fc;
    fc.set("nu.type", "ncdm_standard");
    fc.set("nu.fluid_approximation", "2");
    fc.set("ax.type", "ncdm_axion");
    fc.set("ax.fluid_approximation", "3");
    assert(SynthesisThrows(fc));
  }

  // ── If any family instance sets the key, all of them must ──
  {
    FileContent fc;
    fc.set("nu.type", "ncdm_standard");
    fc.set("nu.fluid_approximation", "2");
    fc.set("ax.type", "ncdm_axion");
    assert(SynthesisThrows(fc));
  }

  // ── Agreement across types synthesizes once, no throw ──
  {
    FileContent fc;
    fc.set("nu.type", "ncdm_standard");
    fc.set("nu.fluid_approximation", "1");
    fc.set("ax.type", "ncdm_axion");
    fc.set("ax.fluid_approximation", "1");
    SynthesiseNcdmFluidApproximation(&fc);
    auto v = fc.get<std::string>("ncdm_fluid_approximation");
    assert(v && *v == "1");
  }

  // ── dcdm_wdm ignores the fluid approximation: dot key must be rejected ──
  {
    FileContent fc;
    fc.set("ddm.type", "dcdm_wdm");
    fc.set("ddm.fluid_approximation", "2");
    assert(SynthesisThrows(fc));
  }

  // ── Same rejection for any non-consumer species type ──
  {
    FileContent fc;
    fc.set("x.type", "fluid");
    fc.set("x.fluid_approximation", "2");
    assert(SynthesisThrows(fc));
  }

  // ── Conflicting explicit global + dot key still throws (existing rule) ──
  {
    FileContent fc;
    fc.set("ncdm_fluid_approximation", "2");
    fc.set("ax.type", "ncdm_axion");
    fc.set("ax.fluid_approximation", "3");
    assert(SynthesisThrows(fc));
  }

  // ── Matching explicit global + dot key is fine ──
  {
    FileContent fc;
    fc.set("ncdm_fluid_approximation", "3");
    fc.set("ax.type", "ncdm_axion");
    fc.set("ax.fluid_approximation", "3");
    SynthesiseNcdmFluidApproximation(&fc);
    auto v = fc.get<std::string>("ncdm_fluid_approximation");
    assert(v && *v == "3");
  }

  // ── No dot keys anywhere: no-op, nothing synthesized ──
  {
    FileContent fc;
    fc.set("nu.type", "ncdm_standard");
    SynthesiseNcdmFluidApproximation(&fc);
    assert(!fc.get<std::string>("ncdm_fluid_approximation"));
  }

  // ── A stray key with no matching <name>.type is not ours to judge:
  //    left unread for the generic unused-parameter warning ──
  {
    FileContent fc;
    fc.set("ghost.fluid_approximation", "1");
    SynthesiseNcdmFluidApproximation(&fc);
    assert(!fc.get<std::string>("ncdm_fluid_approximation"));
    assert(Contains(fc.unread_parameters(), "ghost.fluid_approximation"));
  }

  // ── CreateAllNcdmInstances<T>: the shared factory loop (#378) ──
  background pba{};
  pba.H0 = 2.2e-4;

  {
    FileContent fc;
    fc.set("nu1.type", "ncdm_standard");
    fc.set("nu1.m", "0.06");
    fc.set("nu1.T", "0.71611");
    fc.set("nu2.type", "ncdm_standard");
    fc.set("nu2.m", "0.1");
    fc.set("nu2.T", "0.71611");
    fc.set("other.type", "cdm");  // must not be picked up

    NcdmSettings settings = TestSettings();
    SpeciesBuildContext ctx{};
    ctx.pfc           = &fc;
    ctx.pba           = &pba;
    ctx.ncdm_settings = &settings;

    auto result = CreateAllNcdmInstances<NCDMSpecies>(ctx);
    assert(result.size() == 2);
    assert(result[0].key == "nu1");
    assert(result[1].key == "nu2");
    assert(result[0].species && result[1].species);
    assert(result[0].species->name() == "nu1");

    // The loop must mark each <instance>.type consumed (instances_with does not).
    auto unread = fc.unread_parameters();
    assert(!Contains(unread, "nu1.type"));
    assert(!Contains(unread, "nu2.type"));
    assert(Contains(unread, "other.type"));
  }

  {
    FileContent fc;
    fc.set("ax.type", "ncdm_axion");
    fc.set("ax.m", "1.0");
    fc.set("ax.T", "0.39");

    NcdmSettings settings = TestSettings();
    SpeciesBuildContext ctx{};
    ctx.pfc           = &fc;
    ctx.pba           = &pba;
    ctx.ncdm_settings = &settings;

    auto result = CreateAllNcdmInstances<AxionNCDMSpecies>(ctx);
    assert(result.size() == 1);
    assert(result[0].key == "ax");
    // Axion normalization: deg = 1 is one bosonic dof, so N_eff differs from
    // a fermionic ncdm_standard at the same T; just probe it exists and is > 0.
    assert(result[0].species->GetOmega0() > 0.);
  }

  printf("ncdm_family_test: all assertions passed\n");
  return 0;
}
