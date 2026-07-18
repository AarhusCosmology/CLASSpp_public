#include <cassert>
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include "base_species.h"
#include "composite_species.h"

// Minimal concrete species: configurable name + energy_type, preset Rho.
// Implements only the 7 pure virtuals of BaseSpecies; Rho ignores pvecback.
class BudgetFake : public BaseSpecies {
 public:
  BudgetFake(std::string name, EnergyType et, double rho)
      : BaseSpecies(std::move(name), et), rho_(rho) {}
  void RegisterBackgroundIndices(int&) override {}
  void ComputeBackground(double, const double*, double*) override {}
  double Rho(const double*) const override {
    return rho_;
  }
  double P(const double*) const override {
    return 0.;
  }
  void PerturbDerivs(const PerturbLayout&,
                     double,
                     const double*,
                     double*,
                     const perturb_parameters_and_workspace&) const override {}
  StressEnergyContribution StressEnergy(const PerturbLayout&,
                                        const perturb_vector*,
                                        const double*,
                                        const double*,
                                        const perturb_workspace*) const override {
    return {};
  }
  double GetOmega0() const override {
    return 0.;
  }

 private:
  double rho_;
};

// CompositeSpecies::children_ is protected; this exposes AddChild for the test.
class TestComposite : public CompositeSpecies {
 public:
  TestComposite() : CompositeSpecies("test_composite", EnergyType::Other) {}
  void AddChild(std::unique_ptr<BaseSpecies> c) {
    children_.push_back(std::move(c));
  }
};

int main() {
  // 1) BudgetBucketOf maps energy_type for non-NCDM species.
  BudgetFake rad("g", BaseSpecies::EnergyType::Radiation, 1.0);
  BudgetFake mat("cdm", BaseSpecies::EnergyType::Matter, 1.0);
  BudgetFake de("lambda", BaseSpecies::EnergyType::DarkEnergy, 1.0);
  BudgetFake oth("scf", BaseSpecies::EnergyType::Other, 1.0);
  assert(BudgetBucketOf(rad) == BudgetBucket::Radiation);
  assert(BudgetBucketOf(mat) == BudgetBucket::NonRelativistic);
  assert(BudgetBucketOf(de) == BudgetBucket::Other);
  assert(BudgetBucketOf(oth) == BudgetBucket::Other);

  // 2) AppendBudgetLines: one line per child, omega = Rho/rho_crit,
  //    label = child name(), bucket = BudgetBucketOf(child). (Rho ignores the
  //    null pvecback here.)
  TestComposite comp;
  comp.AddChild(std::make_unique<BudgetFake>("DCDM", BaseSpecies::EnergyType::Matter, 2.0));
  comp.AddChild(std::make_unique<BudgetFake>("DR", BaseSpecies::EnergyType::Radiation, 3.0));
  std::vector<BudgetLine> lines;
  comp.AppendBudgetLines(/*pvecback_today=*/nullptr, /*rho_crit=*/10.0, lines);
  assert(lines.size() == 2);
  assert(lines[0].label == "DCDM");
  assert(lines[0].bucket == BudgetBucket::NonRelativistic);
  assert(std::fabs(lines[0].omega - 0.2) < 1e-12);
  assert(lines[1].label == "DR");
  assert(lines[1].bucket == BudgetBucket::Radiation);
  assert(std::fabs(lines[1].omega - 0.3) < 1e-12);

  printf("species_budget_test: all assertions passed\n");
  return 0;
}
