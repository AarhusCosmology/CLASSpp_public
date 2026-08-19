#include "composite_species.h"

#include "perturbations.h"  // perturb_vector, perturb_workspace, contexts, precision

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

void CompositeSpecies::ComputeBackground(double a, const double* pvecback_B, double* pvecback) {
  for (auto& child : children_)
    child->ComputeBackground(a, pvecback_B, pvecback);
}

void CompositeSpecies::BackgroundDerivs(double tau,
                                        const double* y,
                                        double* dy,
                                        const double* pvecback) {
  for (auto& child : children_)
    child->BackgroundDerivs(tau, y, dy, pvecback);
}

void CompositeSpecies::BackgroundDerivsDiagonal(double tau,
                                                const double* y,
                                                double* diag,
                                                const double* pvecback) {
  for (auto& child : children_)
    child->BackgroundDerivsDiagonal(tau, y, diag, pvecback);
}

void CompositeSpecies::FinalizeBackground(double a,
                                          double H,
                                          const double* pvecback_B,
                                          double* pvecback) {
  for (auto& child : children_)
    child->FinalizeBackground(a, H, pvecback_B, pvecback);
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

std::unique_ptr<BaseSpecies::PerturbLayout> CompositeSpecies::CreatePerturbLayout() const {
  auto l = std::make_unique<PerturbLayout>();
  l->child_layouts.reserve(children_.size());
  for (const auto& c : children_)
    l->child_layouts.push_back(c->CreatePerturbLayout());
  return l;
}

void CompositeSpecies::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                   perturb_vector* pv,
                                                   const precision* ppr,
                                                   int& index_pt,
                                                   const perturb_workspace* ppw,
                                                   int gauge) {
  auto& my = static_cast<PerturbLayout&>(base);
  for (size_t i = 0; i < children_.size(); ++i)
    children_[i]->RegisterPerturbationIndices(*my.child_layouts[i], pv, ppr, index_pt, ppw, gauge);
}

void CompositeSpecies::ApplyInitialConditions(const BaseSpecies::PerturbLayout& base,
                                              double* y,
                                              const PerturbIcContext& ctx) {
  const auto& my = static_cast<const PerturbLayout&>(base);
  for (size_t i = 0; i < children_.size(); ++i)
    children_[i]->ApplyInitialConditions(*my.child_layouts[i], y, ctx);
}

BaseSpecies::StressEnergyContribution CompositeSpecies::StressEnergy(
    const BaseSpecies::PerturbLayout& base,
    const perturb_vector* pv,
    const double* y,
    const double* pvecback,
    const perturb_workspace* ppw) const {
  const auto& my = static_cast<const PerturbLayout&>(base);
  StressEnergyContribution se;
  for (size_t i = 0; i < children_.size(); ++i)
    se += children_[i]->StressEnergy(*my.child_layouts[i], pv, y, pvecback, ppw);
  return se;
}

void CompositeSpecies::DelegateTally(const BaseSpecies::PerturbLayout& base,
                                     const perturb_vector* pv,
                                     const double* y,
                                     const double* pvecback,
                                     const perturb_workspace* ppw,
                                     StressEnergyContribution& total,
                                     StressEnergyContribution& total_cold,
                                     StressEnergyContribution& total_warm) const {
  const auto& my = static_cast<const PerturbLayout&>(base);
  for (size_t i = 0; i < children_.size(); ++i)
    children_[i]->TallyStressEnergy(*my.child_layouts[i],
                                    pv,
                                    y,
                                    pvecback,
                                    ppw,
                                    total,
                                    total_cold,
                                    total_warm);
}

void CompositeSpecies::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                                     double tau,
                                     const double* y,
                                     double* dy,
                                     const perturb_parameters_and_workspace& ppaw) const {
  const auto& my = static_cast<const PerturbLayout&>(base);
  for (size_t i = 0; i < children_.size(); ++i)
    children_[i]->PerturbDerivs(*my.child_layouts[i], tau, y, dy, ppaw);
  AddCouplingDerivs(tau, y, dy, ppaw);  // two-phase contract: children first, then coupling
}

void CompositeSpecies::PerturbDerivsDiagonal(const BaseSpecies::PerturbLayout& base,
                                             double tau,
                                             const double* y,
                                             double* diag,
                                             const perturb_parameters_and_workspace& ppaw) const {
  const auto& my = static_cast<const PerturbLayout&>(base);
  for (size_t i = 0; i < children_.size(); ++i)
    children_[i]->PerturbDerivsDiagonal(*my.child_layouts[i], tau, y, diag, ppaw);
  /* No coupling phase here, unlike PerturbDerivs: a composite whose coupling has a
     diagonal overrides this outright (dncdm_inv does), and the base class has no
     AddCouplingDerivs counterpart to call. */
}

void CompositeSpecies::FillSources(const BaseSpecies::PerturbLayout& base,
                                   const double* y,
                                   const double* dy,
                                   PerturbSourceContext& ctx) const {
  const auto& my = static_cast<const PerturbLayout&>(base);
  for (size_t i = 0; i < children_.size(); ++i)
    children_[i]->FillSources(*my.child_layouts[i], y, dy, ctx);
}

void CompositeSpecies::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                                     double* y,
                                                     const PerturbIcContext& ctx) {
  const auto& my = static_cast<const PerturbLayout&>(base);
  for (size_t i = 0; i < children_.size(); ++i)
    children_[i]->PerturbSynchronousToNewtonian(*my.child_layouts[i], y, ctx);
}

void CompositeSpecies::CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_base,
                                                     const BaseSpecies::PerturbLayout& new_base,
                                                     const double* old_y,
                                                     double* new_y,
                                                     const PerturbSwitchContext& ctx) const {
  const auto& old_l = static_cast<const PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const PerturbLayout&>(new_base);
  for (size_t i = 0; i < children_.size(); ++i)
    children_[i]->CopyPerturbationsAcrossSwitch(*old_l.child_layouts[i],
                                                *new_l.child_layouts[i],
                                                old_y,
                                                new_y,
                                                ctx);
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

bool CompositeSpecies::HasWarmMatter() const {
  // Per-child scan: the composite's own OR-ed predicates can't answer this
  // (a cold parent + warm daughter would read ClustersAsMatter && IsCold).
  for (const auto& child : children_)
    if (child->HasWarmMatter())
      return true;
  return false;
}

void CompositeSpecies::FinalizeMatterClassification() {
  BaseSpecies::
      FinalizeMatterClassification();  // stamp this composite's own (unused but valid) bits
  for (auto& child : children_)
    child->FinalizeMatterClassification();  // children's bits drive DelegateTally (Task 3)
}
