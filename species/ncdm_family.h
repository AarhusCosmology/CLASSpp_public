#pragma once

#include <memory>
#include <string>
#include <vector>

#include "parser.h"
#include "species/species_build_context.h"

/**
 * Shared per-instance CreateAll loop for the leaf NCDM-family factories
 * (ncdm_standard, ncdm_greybody, ncdm_axion, ncdm_self_interacting).
 * T must expose a static kTypeName and a
 * T(FileContent*, instance_name, NcdmSettings, background*, BackgroundModule*)
 * constructor. Factories with extra pre-steps (legacy transmute, legacy-key
 * rejection) run those first and then delegate here.
 */
template <class T>
std::vector<Named> CreateAllNcdmInstances(const SpeciesBuildContext& ctx) {
  const auto instances = ctx.pfc->instances_with("type", std::string(T::kTypeName));

  std::vector<Named> result;
  result.reserve(instances.size());
  for (const auto& name : instances) {
    // Mark <instance>.type as consumed (instances_with does not mark anything as read)
    (void) ctx.pfc->get<std::string>(name + ".type");

    auto sp = std::make_unique<T>(ctx.pfc, name, *ctx.ncdm_settings, ctx.pba, ctx.bgm);
    result.push_back({name, std::move(sp)});
  }
  return result;
}

/**
 * Synthesise the family-wide ncdm_fluid_approximation precision parameter from
 * per-instance <name>.fluid_approximation dot keys. Gathers every NCDM-family
 * type whose perturbation code consults index_ap_ncdmfa (not just
 * ncdm_standard), enforces the "identical for all instances" rule across the
 * whole family, and throws std::invalid_argument when the key appears on a
 * species type that ignores the fluid approximation (e.g. dcdm_wdm).
 */
void SynthesiseNcdmFluidApproximation(FileContent* pfc);
