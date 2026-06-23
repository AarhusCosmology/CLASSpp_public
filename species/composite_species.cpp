#include "composite_species.h"

void CompositeSpecies::RegisterBackgroundIndices(int& index_bg) {
  for (auto& child : children_)
    child->RegisterBackgroundIndices(index_bg);
}

void CompositeSpecies::RegisterIntegrationIndices(int& index_bi) {
  for (auto& child : children_)
    child->RegisterIntegrationIndices(index_bi);
}

void CompositeSpecies::RegisterTransferSourceIndices(int& index_tp,
                                                     const SourceRequestContext& ctx) {
  for (auto& child : children_)
    child->RegisterTransferSourceIndices(index_tp, ctx);
}

void CompositeSpecies::SetBackgroundModule(const BackgroundModule* bgm) {
  for (auto& child : children_)
    child->SetBackgroundModule(bgm);
}

void CompositeSpecies::SetThermodynamicsModule(const ThermodynamicsModule* thm) {
  for (auto& child : children_)
    child->SetThermodynamicsModule(thm);
}

void CompositeSpecies::SetPerturbs(const perturbs* ppt) {
  for (auto& child : children_)
    child->SetPerturbs(ppt);
}

void CompositeSpecies::SetBackgroundInitialConditions(const BackgroundICContext& ctx) {
  for (auto& child : children_)
    child->SetBackgroundInitialConditions(ctx);
}

void CompositeSpecies::ComputeBackground(double a_rel, const double* pvecback_B, double* pvecback) {
  for (auto& child : children_)
    child->ComputeBackground(a_rel, pvecback_B, pvecback);
}

void CompositeSpecies::BackgroundDerivs(double tau,
                                        const double* y,
                                        double* dy,
                                        const double* pvecback) {
  for (auto& child : children_)
    child->BackgroundDerivs(tau, y, dy, pvecback);
}

double CompositeSpecies::Rho(const double* pvecback) const {
  double rho = 0.;
  for (const auto& child : children_)
    rho += child->Rho(pvecback);
  return rho;
}

double CompositeSpecies::P(const double* pvecback) const {
  double p = 0.;
  for (const auto& child : children_)
    p += child->P(pvecback);
  return p;
}

double CompositeSpecies::PPrime(double a,
                                double H,
                                const double* pvecback_B,
                                const double* pvecback) const {
  double pp = 0.;
  for (const auto& child : children_)
    pp += child->PPrime(a, H, pvecback_B, pvecback);
  return pp;
}

double CompositeSpecies::FreestreamingRho(const double* pvecback) const {
  double rho = 0.;
  for (const auto& child : children_)
    rho += child->FreestreamingRho(pvecback);
  return rho;
}

void CompositeSpecies::AddCouplingDerivs(double /*tau*/,
                                         const double* /*y*/,
                                         double* /*dy*/,
                                         const perturb_parameters_and_workspace& /*ppaw*/) const {}

bool CompositeSpecies::ClustersAsMatter() const {
  for (const auto& child : children_)
    if (child->ClustersAsMatter())
      return true;
  return false;
}

bool CompositeSpecies::IsColdMatterSpecies() const {
  for (const auto& child : children_)
    if (child->IsColdMatterSpecies())
      return true;
  return false;
}
