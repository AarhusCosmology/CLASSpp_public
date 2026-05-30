#include "composite_species.h"

void CompositeSpecies::RegisterBackgroundIndices(int& index_bg) {
  for (auto& child : children_)
    child->RegisterBackgroundIndices(index_bg);
}

void CompositeSpecies::RegisterIntegrationIndices(int& index_bi) {
  for (auto& child : children_)
    child->RegisterIntegrationIndices(index_bi);
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

void CompositeSpecies::SetBackgroundInitialConditions(double a_rel, double* pvecback_integration) {
  for (auto& child : children_)
    child->SetBackgroundInitialConditions(a_rel, pvecback_integration);
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

double CompositeSpecies::DpDloga(const double* pvecback) const {
  double dp = 0.;
  for (const auto& child : children_)
    dp += child->DpDloga(pvecback);
  return dp;
}

double CompositeSpecies::FreestreamingRho(const double* pvecback) const {
  double rho = 0.;
  for (const auto& child : children_)
    rho += child->FreestreamingRho(pvecback);
  return rho;
}

void CompositeSpecies::PerturbDerivs(double tau,
                                     const double* y,
                                     double* dy,
                                     const perturb_parameters_and_workspace& ppaw) {
  for (auto& child : children_)
    child->PerturbDerivs(tau, y, dy, ppaw);
  AddCouplingDerivs(tau, y, dy, ppaw);
}

void CompositeSpecies::PerturbVectorDerivs(double tau,
                                           const double* y,
                                           double* dy,
                                           const perturb_parameters_and_workspace& ppaw) {
  for (auto& child : children_)
    child->PerturbVectorDerivs(tau, y, dy, ppaw);
}

void CompositeSpecies::PerturbTensorDerivs(double tau,
                                           const double* y,
                                           double* dy,
                                           const perturb_parameters_and_workspace& ppaw) {
  for (auto& child : children_)
    child->PerturbTensorDerivs(tau, y, dy, ppaw);
}

double CompositeSpecies::Delta(const perturb_vector* pv,
                               const double* y,
                               const double* pvecback,
                               const perturb_workspace* ppw) const {
  double rho_delta = 0., rho_total = 0.;
  for (const auto& child : children_) {
    const double rho  = child->Rho(pvecback);
    rho_delta        += rho * child->Delta(pv, y, pvecback, ppw);
    rho_total        += rho;
  }
  return (rho_total > 0.) ? rho_delta / rho_total : 0.;
}

double CompositeSpecies::Theta(const perturb_vector* pv,
                               const double* y,
                               const double* pvecback,
                               const perturb_workspace* ppw) const {
  double rho_plus_p_theta = 0., rho_plus_p_total = 0.;
  for (const auto& child : children_) {
    const double rho_plus_p  = child->Rho(pvecback) + child->P(pvecback);
    rho_plus_p_theta        += rho_plus_p * child->Theta(pv, y, pvecback, ppw);
    rho_plus_p_total        += rho_plus_p;
  }
  return (rho_plus_p_total > 0.) ? rho_plus_p_theta / rho_plus_p_total : 0.;
}

double CompositeSpecies::DeltaP(const perturb_vector* pv,
                                const double* y,
                                const double* pvecback,
                                const perturb_workspace* ppw) const {
  double dp = 0.;
  for (const auto& child : children_)
    dp += child->DeltaP(pv, y, pvecback, ppw);
  return dp;
}

void CompositeSpecies::AddCouplingDerivs(double /*tau*/,
                                         const double* /*y*/,
                                         double* /*dy*/,
                                         const perturb_parameters_and_workspace& /*ppaw*/) {}

double CompositeSpecies::RhoPlusPShear(const perturb_vector* pv,
                                       const double* y,
                                       const double* pvecback,
                                       const perturb_workspace* ppw) const {
  double s = 0.;
  for (const auto& child : children_)
    s += child->RhoPlusPShear(pv, y, pvecback, ppw);
  return s;
}

bool CompositeSpecies::IsMatterSpecies() const {
  for (const auto& child : children_)
    if (child->IsMatterSpecies())
      return true;
  return false;
}

bool CompositeSpecies::IsColdMatterSpecies() const {
  for (const auto& child : children_)
    if (child->IsColdMatterSpecies())
      return true;
  return false;
}

double CompositeSpecies::MatterRho(const double* pvecback) const {
  double r = 0.;
  for (const auto& child : children_)
    r += child->MatterRho(pvecback);
  return r;
}

double CompositeSpecies::MatterRhoDelta(const perturb_vector* pv,
                                        const double* y,
                                        const double* pvecback,
                                        const perturb_workspace* ppw) const {
  double rd = 0.;
  for (const auto& child : children_)
    rd += child->MatterRhoDelta(pv, y, pvecback, ppw);
  return rd;
}

double CompositeSpecies::MatterRhoPlusPTheta(const perturb_vector* pv,
                                             const double* y,
                                             const double* pvecback,
                                             const perturb_workspace* ppw) const {
  double rpt = 0.;
  for (const auto& child : children_)
    rpt += child->MatterRhoPlusPTheta(pv, y, pvecback, ppw);
  return rpt;
}

double CompositeSpecies::MatterRhoPlusP(const double* pvecback) const {
  double rp = 0.;
  for (const auto& child : children_)
    rp += child->MatterRhoPlusP(pvecback);
  return rp;
}
