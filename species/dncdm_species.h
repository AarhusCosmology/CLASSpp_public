#pragma once
#include <cmath>
#include <memory>
#include <optional>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include "../species/ncdm_base_species.h"
#include "../species/species_build_context.h"
#include "background.h"
#include "perturbations.h"

class BackgroundModule;

/**
 * Decaying Non-Cold Dark Matter (DNCDM).
 * Inherits NCDMBaseSpecies; owns per-species quadrature, distribution function,
 * Gamma decay rate, and dq volume elements absorbed from DecayDRProperties.
 */
class DNCDMSpecies : public NCDMBaseSpecies {
 public:
  static constexpr std::string_view kTypeName = "ncdm_decay_dr";

  // Reads all DNCDM-specific parameters from the dot-syntax instance
  // identified by instance_name (e.g. "dncdm1").
  DNCDMSpecies(FileContent* pfc,
               const std::string& instance_name,
               const NcdmSettings& settings,
               const background* pba,
               const BackgroundModule* bgm);

  // Accessors for deferred closure (used by CreateAll after construction)
  const std::optional<double>& Omega_ini_pending() const {
    return Omega_ini_pending_;
  }
  const std::optional<double>& Neff_ini_pending() const {
    return Neff_ini_pending_;
  }
  const std::optional<double>& Omega_dncdmdr_pending() const {
    return Omega_dncdmdr_pending_;
  }

  /** True iff this flavor is normalized by initial abundance (Omega_ini/omega_ini/Neff_ini) —
   *  the mode that needs the Omega_dncdmdr fixed-point shoot for closure (vs combined mode,
   *  which shoots deg). */
  bool InitialAbundanceMode() const {
    return Omega_ini_pending_.has_value() || Neff_ini_pending_.has_value();
  }

  /** Backfill this species' today density fraction Omega0_ from its integrated density
   *  at a=1. In combined/initial modes the matter child is normalized via `deg`, never
   *  SetOmega0, so GetOmega0() stays 0 — which drops it from fnu (GetOmega0NcdmTot) and
   *  the budget-print neutrino line. BackgroundModule calls this post-integration.
   *  (SetOmega0 is protected; this public wrapper lets the module trigger the backfill.) */
  void BackfillOmega0FromToday(const double* pvecback_today, double H0, double h) {
    SetOmega0(Rho(pvecback_today) / (H0 * H0), h);
  }

  // Compute a (deg_guess, dxdy) pair for Newton shooting that varies deg to hit
  // a target today-density Omega_target = (rho_dncdm + rho_dr) / H0^2 at z=0.
  // Ported from input_module.cpp:3759-3800 (single-flavor, no loop).
  std::pair<double, double> DegGuessFromOmegaToday(const SpeciesBuildContext& ctx,
                                                   double Omega_target) const;

  struct Named {
    std::string key;
    std::unique_ptr<DNCDMSpecies> species;
  };

  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);

  // ── Collision-owned mode (DNCDMInvSpecies composite, design §3 table) ─────
  // When true, the parent drops its ln f + separate-dlnfdlnq integration
  // variables for a single f-per-bin variable and its background RHS becomes a
  // no-op; the composite writes the kernel-supplied RHS into the f-slots. Every
  // background column and the entire perturbation surface are unchanged (both
  // modes publish index_bg_lnf_decay_dr1_ / index_bg_dlnfdlnq_decay_). Decay-only
  // runs never set this, so their code path is byte-for-byte unchanged.
  void SetCollisionOwned(bool v) {
    collision_owned_ = v;
  }
  int bi_f_parent_index() const {
    return index_bi_f_parent_;
  }

  // ── Background ──────────────────────────────────────────────────────────
  void RegisterBackgroundIndices(int& index_bg) override;
  void RegisterIntegrationIndices(int& index_bi) override;
  void SetBackgroundInitialConditions(const BackgroundICContext& ctx) override;
  void ComputeBackground(double a, const double* pvecback_B, double* pvecback) override;
  void BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) override;

  double Rho(const double* pvecback) const override {
    return pvecback[index_bg_rho_];
  }
  double P(const double* pvecback) const override {
    return pvecback[index_bg_p_];
  }
  double PPrime(double a,
                double H,
                const double* /*pvecback_B*/,
                const double* pvecback) const override {
    return a * H * (pvecback[index_bg_pseudo_p_] - 5. * pvecback[index_bg_p_]);
  }

  // ── Perturbations ────────────────────────────────────────────────────────

  // Layout-based scalar register (writes both layout and legacy pv arrays).
  void RegisterPerturbationIndices(BaseSpecies::PerturbLayout& layout,
                                   perturb_vector* pv,
                                   const precision* ppr,
                                   int& index_pt,
                                   const perturb_workspace* ppw,
                                   int gauge) override;

  // Scalar PerturbDerivs (called by DNCDM_DR_Species composite with my.dncdm).
  void PerturbDerivs(const BaseSpecies::PerturbLayout& layout,
                     double tau,
                     const double* y,
                     double* dy,
                     const perturb_parameters_and_workspace& ppaw) const override;

  void ApplyInitialConditions(const BaseSpecies::PerturbLayout& layout,
                              double* y,
                              const PerturbIcContext& ctx) override;

  /** Fused override: one RescaledPerturbations call for delta/theta/shear;
   *  DeltaP uses its own independent loop (different weights), matching the
   *  individual methods' per-term expressions and operand order exactly. */
  StressEnergyContribution StressEnergy(const BaseSpecies::PerturbLayout& layout,
                                        const perturb_vector* pv,
                                        const double* y,
                                        const double* pvecback,
                                        const perturb_workspace* ppw) const override;

  /** Transfer sources for the parent's own delta/theta. The slots have always been
   *  allocated (NCDMBaseSpecies::RegisterTransferSourceIndices) but were never filled
   *  and never named, so the decaying parent was the one member of its own sector
   *  absent from mTk/vTk output while both daughters were present — which made the
   *  sector's total delta unassemblable from a single run. Mirrors DrPsdSpecies. */
  void FillSources(const BaseSpecies::PerturbLayout& layout,
                   const double* y,
                   const double* dy,
                   PerturbSourceContext& ctx) const override;
  void WriteOutputColumns(
      PerturbColumnWriter& writer,
      const PerturbationsModule& mod,
      file_format fmt,
      TransferColumnSection section = TransferColumnSection::all) const override;

  bool IsFreestreaming() const override {
    return true;
  }
  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override;

  // ── Accessors for DNCDM_DR_Species coupling ───────────────────────────────
  int bg_number_index() const {
    return index_bg_number_;
  }
  int bg_pseudo_p_index() const {
    return index_bg_pseudo_p_;
  }
  int bg_lnf_index() const {
    return index_bg_lnf_decay_dr1_;
  }
  int bg_dlnfdlnq_index() const {
    return index_bg_dlnfdlnq_decay_;
  }
  int bg_dlnfdlnq_sep_index() const {
    return index_bg_dlnfdlnq_sep_;
  }
  int bi_lnf_index() const {
    return index_bi_lnf_decay_dr1_;
  }
  int bi_dlnfdlnq_sep_index() const {
    return index_bi_dlnfdlnq_separate_decay_;
  }

  double Gamma() const {
    return Gamma_;
  }
  const std::vector<double>& dq() const {
    return dq_;
  }
  double GetMass() const {
    return M_;
  }
  const std::vector<double>& GetQ() const {
    return q_;
  }

  // Override GetRescaledParameters to use dq_ (not standard w_ weights)
  std::tuple<double, double> GetRescaledParameters(double a,
                                                   const double* lnf_array) const override;

  /**
   * Returns rescaled (delta, theta, shear) for this decaying NCDM flavor.
   * Rescaling subtracts a common lnN from every lnf to prevent exp(lnf)
   * underflow near the precision floor; lnN cancels in the delta/theta/shear
   * ratios, so this is mathematically equivalent to the unrescaled form but
   * numerically stable.
   *
   * @param layout this species' per-pv NCDM layout (provides index_per_q).
   * @param a Scale factor at the current integration time.
   * @param k Fourier wavenumber of the perturbation mode.
   * @param ppw Perturbation workspace holding the current phase-space state.
   */
  std::tuple<double, double, double> RescaledPerturbations(
      const NCDMBaseSpecies::PerturbLayout& layout,
      double a,
      double k,
      const perturb_workspace* ppw) const;

 protected:
  double GetDlnf0Dlnq(int iq, const double* pvecback) const override {
    return pvecback[index_bg_dlnfdlnq_decay_ + iq];
  }
  double GetW0ForGwSource(int iq, const double* pvecback) const override {
    return dq_[iq] * std::exp(pvecback[index_bg_lnf_decay_dr1_ + iq]);
  }

 private:
  const background* pba_;

  // Deferred closure stash (instance-name constructor only).
  // Set when the user specifies Omega_ini/omega_ini or Neff_ini; CreateAll
  // applies SetDeg_from_Omega_ini once a_ini is available.
  std::optional<double> Omega_ini_pending_;
  std::optional<double> Neff_ini_pending_;
  // Today-density target: shoot deg so that (rho_dncdm+rho_dr)/H0^2 == this value at z=0.
  std::optional<double> Omega_dncdmdr_pending_;

  // Absorbed from DecayDRProperties
  double Gamma_ = 0.;
  std::vector<double> dq_;

  // Collision-owned mode (set by the DNCDMInvSpecies factory before any index
  // registration). kFParentFloor keeps ln f defined as gains repopulate f from 0.
  // In this mode the perturbation state variable is F = δf, not the normalized Ψ
  // (#386) — see PerturbDerivs / RescaledPerturbations.
  bool collision_owned_ = false;
  /** Soft floor added to the collision-owned parent's occupation before its ln f
   *  column is published (see ComputeBackground). It must sit comfortably ABOVE the
   *  integrator's own noise level, which after the kFScale rescale is
   *  `threshold/kFScale = (abstol/rtol)/kFScale` ~ 1e-159 at the default
   *  tol_background_integration -- i.e. the column should plateau while the state is
   *  still well resolved, never track it down into the noise. 1e-100 keeps ~59
   *  decades of margin, and a species at 1e-100 of its initial occupation is
   *  irrelevant to every observable, so there is nothing to gain by going lower and
   *  a real failure mode to be bought by it. */
  static constexpr double kFParentFloor = 1e-100;
  // Set from STORED table rows only (BackgroundModule::StoringBackgroundTable) --
  // never from the evolver's trial states, which are allowed to be unphysical.
  // See FloorReport.
  mutable bool floor_touched_   = false;
  mutable bool floor_left_      = false;
  mutable bool negative_f_rows_ = false;

 public:
  /** Units of the collision-owned parent's background f-slot: the integrator carries
   *  kFScale*f, not f.
   *
   *  Both solvers control error with max(threshold, |y|) against a HARD-CODED
   *  abstol = 1e-15 that tol_background_integration cannot move (rkdp45:
   *  threshold = abstol/rtol, so the effective absolute bound is abstol again;
   *  ndf15: threshold = abstol outright). This species' f starts near 0.5 and decays
   *  hundreds of decades below that, so unrescaled it is integrated as noise the
   *  moment it drops under ~1e-9 -- f chatters against zero, the published ln f
   *  column acquires square wells, and the perturbation module's cubic spline of that
   *  column turns P(k) into NaN.
   *
   *  Rescaling is EXACT, a change of units rather than an approximation, and it moves
   *  the noise floor below anything this species can reach while carrying energy. It
   *  works here because the decay is MONOTONE: ln f falls smoothly, reaches
   *  kFParentFloor and stays there, and a flat plateau is something a spline can
   *  represent. A configuration in which inverse decays genuinely REPOPULATE the
   *  parent from ~0 would come back off that plateau and needs the FLOOR addressed
   *  too, not just the scale.
   *
   *  Public because the composite converts at the kernel boundary; see
   *  DNCDMInvSpecies::ApplyKernelBackgroundDerivs. */
  static constexpr double kFScale = 1e150;

  /** Whether the collision-owned parent's occupation went under kFParentFloor during
   *  the background solve, and whether it later came back out.
   *
   *  `touched` alone is benign: the species has decayed into irrelevance and its
   *  published ln f column simply plateaus, which splines fine.
   *
   *  `left` is the one that matters -- inverse decays repopulated the parent from
   *  below the floor, so the plateau was hiding real physics AND the column now has a
   *  rise whose height is set by the floor rather than by the solution. The soft floor
   *  keeps such a recovery smooth, but it still says kFParentFloor is too high for the
   *  configuration being run.
   *
   *  `negative` is stronger than either: the integrator produced a NEGATIVE occupation
   *  on a stored row. Unphysical whatever the floor is, and clamped to zero only so
   *  the logarithm stays total.
   *
   *  ALL THREE are evaluated on STORED rows only. An earlier version tested every
   *  ComputeBackground call, including the trial states the evolver proposes and
   *  rejects, and its false positives condemned configurations whose accepted solution
   *  stayed dozens of decades clear of the floor. */
  struct FloorReport {
    bool touched  = false;
    bool left     = false;
    bool negative = false;
  };
  FloorReport BackgroundFloorReport() const {
    return {floor_touched_, floor_left_, negative_f_rows_};
  }

 private:
  // Background indices
  int index_bg_number_   = -1;
  int index_bg_pseudo_p_ = -1;

  int index_bi_lnf_decay_dr1_           = -1;
  int index_bi_dlnfdlnq_separate_decay_ = -1;
  int index_bi_f_parent_                = -1;  // collision-owned: single f-per-bin variable

  int index_bg_lnf_decay_dr1_  = -1;
  int index_bg_dlnfdlnq_decay_ = -1;
  int index_bg_dlnfdlnq_sep_   = -1;

  // Perturbation indices
  int index_pt_psi0_ = -1;
};
