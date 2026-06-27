#include <cassert>
#include <cstdio>
#include <memory>

#include "background.h"
#include "cdm.h"
#include "composite_species.h"

// CompositeSpecies provides generic PerturbDerivs/StressEnergy and children_ is protected;
// this minimal concrete composite lets the unit test build one and add children.
class TestComposite : public CompositeSpecies {
 public:
  TestComposite() : CompositeSpecies("test_composite", EnergyType::Matter) {}
  void AddChild(std::unique_ptr<BaseSpecies> c) {
    children_.push_back(std::move(c));
  }
};

// Minimal plain species whose StressEnergy returns a preset contribution and whose
// matter classification is set explicitly. TallyStressEnergy passes pv/y/pvecback/ppw
// straight through to StressEnergy (never dereferencing them), so the nullptrs the
// tally test hands in are fine. Used to verify per-child bucketing/delegation.
class FakeSpecies : public BaseSpecies {
 public:
  FakeSpecies(BaseSpecies::StressEnergyContribution se, bool clusters, bool cold)
      : BaseSpecies("fake", EnergyType::Other), se_(se), clusters_(clusters), cold_(cold) {}
  void RegisterBackgroundIndices(int&) override {}
  void ComputeBackground(double, const double*, double*) override {}
  double Rho(const double*) const override {
    return se_.rho;
  }
  double P(const double*) const override {
    return se_.p;
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
    return se_;
  }
  double GetOmega0() const override {
    return 0.;
  }
  bool ClustersAsMatter() const override {
    return clusters_;
  }
  bool IsColdMatterSpecies() const override {
    return cold_;
  }

 private:
  BaseSpecies::StressEnergyContribution se_;
  bool clusters_;
  bool cold_;
};

int main() {
  background pba{};
  pba.H0 = 1e-4;

  // Cold dark matter: clusters as matter AND is cold. Cached values must match
  // the virtual predicates after FinalizeMatterClassification().
  CDMSpecies cdm(pba, 0.12, /*coupled=*/false);
  cdm.FinalizeMatterClassification();
  assert(cdm.ClustersAsMatterCached() == cdm.ClustersAsMatter());
  assert(cdm.IsColdCached() == cdm.IsColdMatterSpecies());
  assert(cdm.ClustersAsMatterCached() == true);
  assert(cdm.IsColdCached() == true);

  // ── Composite recursion: the reason FinalizeMatterClassification() is virtual ──
  TestComposite comp;
  auto cdm_child    = std::make_unique<CDMSpecies>(pba, 0.12, /*coupled=*/false);
  CDMSpecies* child = cdm_child.get();
  comp.AddChild(std::move(cdm_child));

  // Discriminator: the child cache starts at its false default; ONLY the composite's
  // recursive override flips it. (Delete the child loop in
  // CompositeSpecies::FinalizeMatterClassification and the post-finalize asserts fail.)
  assert(child->ClustersAsMatterCached() == false);
  assert(child->IsColdCached() == false);

  comp.FinalizeMatterClassification();

  assert(child->ClustersAsMatterCached() == child->ClustersAsMatter());  // recursion ran
  assert(child->IsColdCached() == child->IsColdMatterSpecies());
  assert(child->ClustersAsMatterCached() == true);
  assert(child->IsColdCached() == true);

  // The composite's own aggregate bits (the BaseSpecies:: line of the override)
  assert(comp.ClustersAsMatterCached() == comp.ClustersAsMatter());
  assert(comp.IsColdCached() == comp.IsColdMatterSpecies());
  assert(comp.ClustersAsMatterCached() == true);
  assert(comp.IsColdCached() == true);

  // ── Per-child matter tally (I-2): composite delegates the tally per child, and
  // each child folds into a bucket by its OWN cached classification, not the
  // composite's aggregate. The dark-energy (scf-like, w=-1) child clusters as
  // neither cold nor warm matter, so it must drop out of both matter buckets —
  // exactly the Type3 bug fix — while still contributing to the `total` (Einstein)
  // accumulator. ────────────────────────────────────────────────────────────────
  using SE = BaseSpecies::StressEnergyContribution;
  SE se_cold;
  se_cold.rho              = 2.0;
  se_cold.delta_rho        = 0.20;
  se_cold.rho_plus_p_theta = 0.020;
  SE se_warm;
  se_warm.rho              = 1.0;
  se_warm.p                = 0.10;
  se_warm.delta_rho        = 0.10;
  se_warm.rho_plus_p_theta = 0.010;
  SE se_de;
  se_de.rho       = 0.7;
  se_de.p         = -0.7;
  se_de.delta_rho = 0.05;  // w = -1, clusters as neither

  TestComposite tally_comp;
  tally_comp.AddChild(std::make_unique<FakeSpecies>(se_cold, /*clusters=*/true, /*cold=*/true));
  tally_comp.AddChild(std::make_unique<FakeSpecies>(se_warm, /*clusters=*/true, /*cold=*/false));
  tally_comp.AddChild(std::make_unique<FakeSpecies>(se_de, /*clusters=*/false, /*cold=*/false));
  tally_comp.FinalizeMatterClassification();  // stamp the child caches the tally reads

  auto layout = tally_comp.CreatePerturbLayout();
  SE total, total_cold, total_warm;
  tally_comp.TallyStressEnergy(*layout,
                               nullptr,
                               nullptr,
                               nullptr,
                               nullptr,
                               total,
                               total_cold,
                               total_warm);

  // `total` sums every child regardless of matter classification (it feeds Einstein).
  assert(total.rho == se_cold.rho + se_warm.rho + se_de.rho);
  assert(total.p == se_cold.p + se_warm.p + se_de.p);
  assert(total.delta_rho == se_cold.delta_rho + se_warm.delta_rho + se_de.delta_rho);
  // Cold bucket = the cold-clustering child ONLY (no rho-3P proxy; w=-1 child absent).
  assert(total_cold.rho == se_cold.rho);
  assert(total_cold.delta_rho == se_cold.delta_rho);
  assert(total_cold.p == se_cold.p);
  // Warm bucket = the warm-clustering child ONLY.
  assert(total_warm.rho == se_warm.rho);
  assert(total_warm.delta_rho == se_warm.delta_rho);
  assert(total_warm.p == se_warm.p);
  // Discriminator: the w=-1 child (rho=0.7) reached neither matter bucket. A wrong
  // composite-level bucketing would have folded all three into total_cold (3.7).
  assert(total_cold.rho + total_warm.rho == se_cold.rho + se_warm.rho);

  std::printf("composite classification test passed\n");
  return 0;
}
