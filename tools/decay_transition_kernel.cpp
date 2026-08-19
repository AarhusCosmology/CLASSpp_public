#include "decay_transition_kernel.h"

#include <cmath>

#include "errors.h"
#include "quadrature.h"

DecayTransitionKernel::DecayTransitionKernel(GridView parent,
                                             GridView fermion,
                                             GridView boson,
                                             Statistics fermion_stat,
                                             Statistics boson_stat,
                                             Config cfg,
                                             GridView fermion_bg,
                                             GridView boson_bg)
    : parent_(parent), fermion_(fermion), boson_(boson),
      fermion_bg_(fermion_bg.n > 0 ? fermion_bg : fermion),
      boson_bg_(boson_bg.n > 0 ? boson_bg : boson), fermion_stat_(fermion_stat),
      boson_stat_(boson_stat), cfg_(cfg) {
  separate_bg_grids_ = (fermion_bg_.q != fermion_.q) || (boson_bg_.q != boson_.q);

  // Worst-case nodes per parent: the union partition cuts the band at every fermion
  // cell edge AND every image of a boson cell edge inside it, so at most
  // (n_f + n_b + 1) sub-intervals, each carrying one midpoint node.
  node_stride_ = fermion_.n + boson_.n + 1;
  breakpoints_.assign(fermion_.n + boson_.n + 2, 0.);

  // All hot-path scratch sized once here so PrepareTransitions never allocates.
  const int total = parent_.n * node_stride_;
  node_off_.assign(parent_.n, 0);
  node_cnt_.assign(parent_.n, 0);
  for (int i = 0; i < parent_.n; ++i)
    node_off_[i] = i * node_stride_;  // node_cnt_ is set per call
  q2_star_.assign(total, 0.);
  q3_star_.assign(total, 0.);
  wn_.assign(total, 0.);
  lambda_s_.assign(total, 0.);
  fl_gather_.assign(total, 0.);
  fphi_gather_.assign(total, 0.);
  cos_a_.assign(total, 0.);
  cos_b_.assign(total, 0.);
  cos_g_.assign(total, 0.);
  Pa_prev_.assign(total, 0.);
  Pa_curr_.assign(total, 0.);
  Pb_prev_.assign(total, 0.);
  Pb_curr_.assign(total, 0.);
  Pg_prev_.assign(total, 0.);
  Pg_curr_.assign(total, 0.);
  j_fermion_.assign(total, 0);
  k_boson_.assign(total, 0);
  w_fermion_.assign(static_cast<size_t>(total) * DepositStencil::kWidth, 0.);
  w_boson_.assign(static_cast<size_t>(total) * DepositStencil::kWidth, 0.);
  eps1_.assign(parent_.n, 0.);
  df_bg_H_.assign(parent_.n, 0.);
  df_bg_l_.assign(fermion_.n, 0.);
  df_bg_phi_.assign(boson_.n, 0.);
  fl_state_.assign(fermion_.n, 0.);
  fphi_state_.assign(boson_.n, 0.);
  inv_meas_f_.assign(fermion_.n, 0.);
  for (int k = 0; k < fermion_.n; ++k)
    inv_meas_f_[k] = 1.0 / (fermion_.q[k] * fermion_.q[k] * fermion_.dq[k]);
  inv_meas_b_.assign(boson_.n, 0.);
  for (int k = 0; k < boson_.n; ++k)
    inv_meas_b_[k] = 1.0 / (boson_.q[k] * boson_.q[k] * boson_.dq[k]);
  theta_l_.assign(total, 1.);
  theta_phi_.assign(total, 1.);
  d_l_edge_.assign(static_cast<size_t>(total) * DepositStencil::kWidth, 0.);
  d_phi_edge_.assign(static_cast<size_t>(total) * DepositStencil::kWidth, 0.);
  dev_l_edge_.assign(static_cast<size_t>(total) * DepositStencil::kWidth, 0.);
  dev_phi_edge_.assign(static_cast<size_t>(total) * DepositStencil::kWidth, 0.);
  dth_l_edge_.assign(static_cast<size_t>(total) * DepositStencil::kWidth, 0.);
  dth_phi_edge_.assign(static_cast<size_t>(total) * DepositStencil::kWidth, 0.);
}

void DecayTransitionKernel::TwoBinSplit(
    const GridView& g, double q_star, int& j, double& lambda, bool& clamped, int hint) {
  clamped = false;
  if (q_star <= g.q[0]) {
    j       = 0;
    lambda  = 1.;  // full deposit into the bottom bin
    clamped = (q_star < g.q[0]);
    return;
  }
  if (q_star >= g.q[g.n - 1]) {
    j       = g.n - 2;
    lambda  = 0.;  // full deposit into the top bin (j+1)
    clamped = (q_star > g.q[g.n - 1]);
    return;
  }
  // Bracket g.q[lo] <= q_star < g.q[hi=lo+1]. Both endpoints are already excluded
  // above, so the invariant holds for the full range and for any sub-range the hunt
  // narrows it to; bisection then finds the SAME lo whatever the starting bracket.
  int lo = 0;
  int hi = g.n - 1;
  if (hint > 0 && hint < g.n - 1) {
    int inc = 1;
    if (q_star >= g.q[hint]) {  // hunt upward, doubling the stride
      lo = hint;
      hi = lo + 1;
      while (hi < g.n - 1 && q_star >= g.q[hi]) {
        lo    = hi;
        inc <<= 1;
        hi    = lo + inc;
      }
      if (hi > g.n - 1)
        hi = g.n - 1;
    }
    else {  // hunt downward
      hi = hint;
      lo = hi - 1;
      while (lo > 0 && q_star < g.q[lo]) {
        hi    = lo;
        inc <<= 1;
        lo    = hi - inc;
      }
      if (lo < 0)
        lo = 0;
    }
  }
  while (hi - lo > 1) {
    const int mid = (lo + hi) / 2;
    if (q_star < g.q[mid]) {
      hi = mid;
    }
    else {
      lo = mid;
    }
  }
  j      = lo;
  lambda = (g.q[j + 1] - q_star) / (g.q[j + 1] - g.q[j]);
}

void DecayTransitionKernel::BuildDeposit(const GridView& g,
                                         double q_star,
                                         DepositStencil& out,
                                         int* hint) {
  int j        = 0;
  double lam   = 0.;
  bool clamped = false;
  TwoBinSplit(g, q_star, j, lam, clamped, hint ? *hint : -1);
  if (hint)
    *hint = j;
  out.j0      = j;
  out.w[0]    = lam;
  out.w[1]    = 1.0 - lam;
  out.clamped = clamped;
}

void DecayTransitionKernel::PrepareTransitions(double a,
                                               double m,
                                               double Gamma,
                                               const double* f_H,
                                               const double* f_l,
                                               const double* f_phi,
                                               double a_prime_over_a) const {
  am2_     = a * a * m * m;
  k_coeff_ = a * a * m * Gamma;  // K = a^2 m Gamma (design 4.1)

  // Stiffness cap (Config::max_rate). The per-bin collision rate is K/eps1_i (see
  // the decay-only identity dF_H = -(K/eps1) F_H), in the same conformal 1/Mpc
  // units as a'/a, so the two are directly comparable. eps1 is smallest at the
  // smallest q, so the FASTEST bin sets the constraint:
  //
  //     K_eff = min(K,  max_rate * (a'/a) * eps_min)
  //
  // ONE global factor on K, so all three species' deposits scale together and the
  // exact conservation identities survive (every deposit is linear in K). Capping
  // per bin would rescale the legs of a single transition differently and destroy
  // them. Cheap: one extra pass over the parent grid, only when the cap is armed.
  if (cfg_.max_rate > 0. && a_prime_over_a > 0. && k_coeff_ > 0.) {
    double eps_min = std::sqrt(parent_.q[0] * parent_.q[0] + am2_);
    for (int i = 1; i < parent_.n; ++i)
      eps_min = std::min(eps_min, std::sqrt(parent_.q[i] * parent_.q[i] + am2_));
    const double k_max = cfg_.max_rate * a_prime_over_a * eps_min;
    k_coeff_           = std::min(k_coeff_, k_max);
  }
  clamped_energy_residual_      = 0.;
  split_energy_residual_        = 0.;
  split_energy_residual_pert_   = 0.;
  split_momentum_residual_pert_ = 0.;
  // The rolling Legendre recurrence restarts with every transition build: the next
  // ApplyPerturbationOperator call must be l = 0 (it seeds P_0 for every node).
  next_l_ = 0;
  // Cache the parent f for the operator's quantum-statistics coefficients (valid
  // until the next PrepareTransitions; the caller owns the buffer).
  fH_bg_ = f_H;

  // K = 0: SKIP THE NODE LAYOUT ENTIRELY.
  //
  // Every deposit, every operator coefficient and every diagonal entry is linear in
  // K, so at K = 0 they are all exactly zero and the several-hundred-node layout,
  // gather and deposit build below produce nothing but zeros. Reached for Gamma = 0.
  //
  // Skipping it rather than multiplying by zero is what makes this a COST saving and
  // not merely an arithmetic one: a profile puts ~97% of perturbation runtime in that
  // node walk, so zeroing the coefficients while still walking it saves ~4%.
  //
  // Emptying node_cnt_ is what makes the skip safe for the OTHER entry points without
  // adding an early return to each: ApplyPerturbationOperator(AllL),
  // CollisionDiagonal and PositivityDecomposition all reach node data only through
  // `for (n = node_off_[i]; n < node_off_[i] + node_cnt_[i]; ++n)`, and all three zero
  // their outputs before that loop. eps1_ is still filled because those three read it
  // OUTSIDE the node loop (in the discarded `k_coeff_ * dq1 * q1/eps1_[i]` prefactor),
  // and a stale or zero eps1_ would divide by zero there.
  if (k_coeff_ == 0.) {
    for (int i = 0; i < parent_.n; ++i) {
      eps1_[i]     = std::sqrt(parent_.q[i] * parent_.q[i] + am2_);
      node_cnt_[i] = 0;
      df_bg_H_[i]  = 0.;
    }
    for (int j = 0; j < fermion_.n; ++j)
      df_bg_l_[j] = 0.;
    for (int k = 0; k < boson_.n; ++k)
      df_bg_phi_[k] = 0.;
    return;
  }

  // Rolling deposit brackets. Within a parent bin the node momenta are monotone --
  // q2* ascends along the stratified partition, q3* = eps1 - q2* descends -- and
  // across bins the band edges move smoothly, so carrying these over the whole sweep
  // makes every BuildDeposit an O(1) hunt from the previous answer. Purely a search
  // accelerator; the stencils are unchanged.
  int hint_f = -1, hint_b = -1, hint_fbg = -1, hint_bbg = -1;
  for (int i = 0; i < parent_.n; ++i) {
    const double q1   = parent_.q[i];
    const double dq1  = parent_.dq[i];
    const double eps1 = std::sqrt(q1 * q1 + am2_);
    eps1_[i]          = eps1;
    // band(q1) = [(eps1-q1)/2, (eps1+q1)/2]; width = q1 exactly (qhi = qlo + q1),
    // computed this way so Sum wn_s = q1 with no catastrophic cancellation.
    const double qlo        = 0.5 * (eps1 - q1);
    const double half_width = 0.5 * q1;
    const double fH_i       = f_H[i];

    // Lay out this parent's nodes on the union partition, integrating the emission
    // band piece by piece (see the class doc's stratified-quadrature note).
    const int off = node_off_[i];
    int cnt       = 0;
    // Breakpoints strictly inside the band, in ascending order, produced by MERGING
    // two already-sorted streams rather than sorting -- the fermion cell edges q_k,
    // and the boson's images eps1 - q_k which descend in k. Note the boson edges land
    // inside the band for exactly the same k as the fermion's, since eps1 - qhi = qlo
    // and eps1 - qlo = qhi: the band is symmetric under q2 <-> q3.
    // The band is NOT clipped to the grids: nodes outside stay outside and the
    // deposit clamps them, so clamped_energy_residual and the total band weight
    // (= q1) both stay exact.
    const double qhi    = qlo + q1;
    int nbp             = 0;
    breakpoints_[nbp++] = qlo;
    int af              = 0;
    while (af < fermion_.n && fermion_.q[af] <= qlo)
      ++af;
    int ab = boson_.n - 1;
    while (ab >= 0 && eps1 - boson_.q[ab] <= qlo)
      --ab;
    while (true) {
      const bool have_f = (af < fermion_.n) && (fermion_.q[af] < qhi);
      const bool have_b = (ab >= 0) && (eps1 - boson_.q[ab] < qhi);
      if (!have_f && !have_b)
        break;
      double x;
      if (have_f && (!have_b || fermion_.q[af] <= eps1 - boson_.q[ab]))
        x = fermion_.q[af++];
      else
        x = eps1 - boson_.q[ab--];
      if (x <= breakpoints_[nbp - 1])
        continue;
      // Checked BEFORE the write, not after: a capacity guard that fires once the
      // overflow has already happened is not a guard. Capacity is n_f + n_b + 2 and
      // at most that many can be pushed (qlo, every edge of either grid inside the
      // band, qhi), so this can only fire on a logic error -- and must abort, because
      // dropping a breakpoint silently under-counts the decay rate while leaving
      // every per-transition identity satisfied, invisible to the conservation tests.
      class_test_severe(nbp >= static_cast<int>(breakpoints_.size()),
                        "stratified quadrature: parent bin %d overran the %d-breakpoint "
                        "scratch",
                        i,
                        static_cast<int>(breakpoints_.size()));
      breakpoints_[nbp++] = x;
    }
    if (qhi > breakpoints_[nbp - 1])
      breakpoints_[nbp++] = qhi;

    class_test_severe(nbp - 1 > node_stride_,
                      "stratified quadrature: parent bin %d needs %d nodes but the "
                      "stride is %d",
                      i,
                      nbp - 1,
                      node_stride_);
    // One midpoint node per sub-interval. The integrand is smooth across a cell, and
    // the midpoint rule there already matches an 8-node-per-cell reference to <5e-5.
    for (int c = 0; c + 1 < nbp; ++c) {
      const double ca = breakpoints_[c];
      const double hw = 0.5 * (breakpoints_[c + 1] - ca);
      if (!(hw > 0.))
        continue;
      q2_star_[off + cnt] = ca + hw;
      wn_[off + cnt]      = 2.0 * hw;
      ++cnt;
    }

    node_cnt_[i] = cnt;

    for (int n = off; n < off + cnt; ++n) {
      const double q2 = q2_star_[n];
      const double q3 = eps1 - q2;  // q2* + q3* = eps1 exactly (energy conservation per transition)
      q3_star_[n]     = q3;

      // Gather f with the two-bin split weights; the SAME weights scatter into
      // df_* (transpose consistency -> exact discrete conservation).
      DepositStencil sf;
      DepositStencil sb;
      BuildDeposit(fermion_, q2, sf, &hint_f);
      BuildDeposit(boson_, q3, sb, &hint_b);
      const int jf   = sf.j0;
      const int kb   = sb.j0;
      const bool clf = sf.clamped;
      const bool clb = sb.clamped;
      j_fermion_[n]  = jf;
      k_boson_[n]    = kb;
      for (int e = 0; e < DepositStencil::kWidth; ++e) {
        w_fermion_[DepositStencil::kWidth * n + e] = sf.w[e];
        w_boson_[DepositStencil::kWidth * n + e]   = sb.w[e];
      }
      // f is sampled on the BACKGROUND grid (== the state grid unless the caller
      // supplied separate views, in which case f_l/f_phi are sized to those). The
      // deposit indices above stay on the state grid: gather-of-F and scatter-of-dF
      // must remain transposes for the conservation identities, but WHERE the
      // background occupation is read is a free choice, and reading it coarsely
      // fakes the daughters' injection front (see the ctor doc).
      const auto gather = [](const DepositStencil& s, const double* f) {
        double acc = 0.;
        for (int e = 0; e < DepositStencil::kWidth; ++e)
          acc += s.w[e] * f[s.j0 + e];
        return acc;
      };
      double fl_g   = 0.;
      double fphi_g = 0.;
      if (separate_bg_grids_) {
        DepositStencil sfb;
        DepositStencil sbb;
        BuildDeposit(fermion_bg_, q2, sfb, &hint_fbg);
        BuildDeposit(boson_bg_, q3, sbb, &hint_bbg);
        fl_g   = gather(sfb, f_l);
        fphi_g = gather(sbb, f_phi);
      }
      else {
        fl_g   = gather(sf, f_l);
        fphi_g = gather(sb, f_phi);
      }
      fl_gather_[n]   = fl_g;
      fphi_gather_[n] = fphi_g;

      // Lambda = f_l f_phi (inv) - f_H (dec) + f_H (f_l - f_phi) (qs); drop the
      // inv term when inverse decays are off, the qs term when qs is off.
      double lam = -fH_i;
      if (cfg_.inverse_decays)
        lam += fl_g * fphi_g;
      if (cfg_.quantum_statistics)
        lam += fH_i * (fl_g - fphi_g);
      lambda_s_[n] = lam;

      // Emission-angle cosines (formalism-notes eq. 6.9, massless daughters).
      cos_a_[n] = (2.0 * eps1 * q2 - am2_) / (2.0 * q1 * q2);
      cos_b_[n] = (2.0 * eps1 * q3 - am2_) / (2.0 * q1 * q3);
      cos_g_[n] = 1.0 - am2_ / (2.0 * q2 * q3);

      // Off-grid clamp diagnostic: a clamped deposit preserves number but places
      // energy at the edge bin, not at q*. Book that mismatch (0 with covering
      // grids). Daughter numbers per transition: fermion -dN, boson -2 dN.
      if (clf || clb) {
        const double dN = k_coeff_ * dq1 * (q1 / eps1) * wn_[n] * lam;
        if (clf) {
          double e_dep = 0.;
          for (int e = 0; e < DepositStencil::kWidth; ++e)
            e_dep += sf.w[e] * fermion_.q[jf + e];
          clamped_energy_residual_ += (e_dep - q2) * (-dN);
        }
        if (clb) {
          double e_dep = 0.;
          for (int e = 0; e < DepositStencil::kWidth; ++e)
            e_dep += sb.w[e] * boson_.q[kb + e];
          clamped_energy_residual_ += (e_dep - q3) * (-2.0 * dN);
        }
      }
    }
  }

  // Background scatter into df_bg_* (cached here so ComputeBackgroundDerivs copies
  // it out AND the perturbation operator reuses the SAME rates for -f-dot/f).
  for (int i = 0; i < parent_.n; ++i)
    df_bg_H_[i] = 0.;
  for (int j = 0; j < fermion_.n; ++j)
    df_bg_l_[j] = 0.;
  for (int k = 0; k < boson_.n; ++k)
    df_bg_phi_[k] = 0.;

  // LUMPED DAUGHTER LOSS (positivity). Lambda is linear in each daughter's own
  // occupation, so split it per leg as
  //
  //     Lambda = Lambda0 + g * f_gathered,      g = dLambda/df_daughter,
  //
  // and evaluate the f-linear part at the DEPOSIT BIN's own occupation instead of
  // the gathered one. The gathered value IS sum_e w_e f_{k+e}, so summing the two
  // deposits reproduces dN exactly: the total number leaving each daughter, and
  // therefore the whole parent leg and every number identity, is UNCHANGED. What
  // changes is only where in q the energy leaves, by O(dq) per node, booked in
  // split_energy_residual_.
  //
  // Why it has to be this way: with the exact (lambda, 1-lambda) weights the loss
  // charged to bin k is proportional to its NEIGHBOUR's occupation as well as its
  // own, i.e. the operator has negative off-diagonals. At an injection front, where
  // f_{k+1}/f_k is enormous, an EMPTY bin is then billed for particles it does not
  // hold and the right-hand side drives it negative -- measured at Gamma=1e8, boson
  // bin 6: gain +1.52 against off-diagonal -1.59, state running to f = -1.13. That
  // is a property of the RHS, so no integrator can repair it; strict positivity
  // rejection, clamping f at 0, and evolving sqrt(f) were each tried and each
  // stalled the solver instead, all against this same wall.
  //
  // The lumped loss reads the occupation of the bin it deposits INTO, so when the
  // gather grid is finer than the deposit grid the caller's f arrays cannot supply
  // it. Resample once here rather than skipping the fix on that path: the state grid
  // is an exact integer subsample of the bg grid, so this is a lookup.
  const double* f_l_st   = f_l;
  const double* f_phi_st = f_phi;
  if (separate_bg_grids_) {
    for (int k = 0; k < fermion_.n; ++k) {
      int j      = 0;
      double lam = 0.;
      bool dummy = false;
      TwoBinSplit(fermion_bg_, fermion_.q[k], j, lam, dummy);
      fl_state_[k] = lam * f_l[j] + (1.0 - lam) * f_l[j + 1];
    }
    for (int k = 0; k < boson_.n; ++k) {
      int j      = 0;
      double lam = 0.;
      bool dummy = false;
      TwoBinSplit(boson_bg_, boson_.q[k], j, lam, dummy);
      fphi_state_[k] = lam * f_phi[j] + (1.0 - lam) * f_phi[j + 1];
    }
    f_l_st   = fl_state_.data();
    f_phi_st = fphi_state_.data();
  }
  // ── deposit sweep ────────────────────────────────────────────────────────
  // ~95% of PrepareTransitions.
  auto sweep = [&](int i, double* acc_l, double* acc_phi, double& split_acc, double& clamp_acc) {
    const double q1           = parent_.q[i];
    const double dq1          = parent_.dq[i];
    const double q1_over_eps1 = q1 / eps1_[i];
    const double fH_i         = fH_bg_[i];
    for (int n = node_off_[i]; n < node_off_[i] + node_cnt_[i]; ++n) {
      // Signed number rate for this transition (dN < 0 for net decay).
      // NOTE on the decay-only rung (inverse_decays = no, quantum_statistics = no):
      // there g_l = g_phi = 0, so lumping is algebraically a NO-OP and the rung's
      // physics is untouched -- but it is not bit-identical, because restructuring
      // the deposit into a two-edge loop lets the compiler vectorize and reassociate
      // it differently. Measured on that rung: max |dP/P| = 1.0e-5, median 3.9e-6,
      // i.e. pure last-bit drift amplified by the ODE, ~100x inside this branch's
      // 0.1% verification tolerance. Do not chase it back to MD5 equality; spelling
      // dN as a single expression again does not recover it.
      const double C  = k_coeff_ * dq1 * q1_over_eps1 * wn_[n];
      const double dN = C * lambda_s_[n];

      // Parent: dN/(q1^2 dq1) = (K/(eps1 q1)) wn Lambda (direct continuum form).
      df_bg_H_[i] += dN / (q1 * q1 * dq1);

      // dLambda/df_l and dLambda/df_phi at this node. g_l >= 0 always (absorption
      // plus Pauli blocking); g_phi = f_l - f_H changes sign, negative meaning Bose
      // enhancement outweighs absorption, in which case the term is a gain and is
      // lumped for the same reason -- one uniform O(dq) treatment, no branch, so
      // the RHS stays smooth in the state (a kink there is what stalled the clamp).
      const double g_l   = (cfg_.inverse_decays ? fphi_gather_[n] : 0.) +
                           (cfg_.quantum_statistics ? fH_i : 0.);
      const double g_phi = (cfg_.inverse_decays ? fl_gather_[n] : 0.) -
                           (cfg_.quantum_statistics ? fH_i : 0.);

      // Fermion daughter: split -dN into (jf, jf+1); 1/(q^2 dq) converts number
      // back to a PSD derivative at each deposit bin.
      //
      // Every lumped quantity is written as (gathered value) + (per-bin DEVIATION),
      // and the deviation is taken against the STATE-grid gather, never the
      // background-grid one. That is what makes sum_e w_e (deviation) vanish
      // identically and so leaves dN, the parent's leg and all number identities
      // untouched. Getting this wrong is not subtle in its effect but is easy to
      // write: referencing the finer bg gather here broke numH+numl to 1e-2 at
      // dr_bg_refine = 9 while looking perfect at refine = 1, where the two gathers
      // coincide. The ACCURACY of the total still comes from the bg gather, which is
      // what Lambda and g_l are built from; only the SPLIT between two bins uses
      // state-grid information, which is all the state grid can represent anyway.
      constexpr int W  = DepositStencil::kWidth;
      const int jf     = j_fermion_[n];
      constexpr int nf = DepositStencil::kWidth;
      const double* wf = &w_fermion_[W * n];
      double flg_st    = 0.;
      for (int e = 0; e < nf; ++e)
        flg_st += wf[e] * f_l_st[jf + e];
      // theta is the lumping FRACTION for this node (see LumpTheta): 1 lumps fully,
      // 0 uses the exact split. Cached because the perturbation operator, the
      // diagonal and the positivity decomposition all have to use the SAME blend.
      const double th_l = LumpTheta(g_l, wf, nf, jf, f_l_st, flg_st);
      theta_l_[n]       = th_l;
      for (int e = 0; e < nf; ++e) {
        const int k           = jf + e;
        const double d_l      = th_l * (f_l_st[k] - flg_st);
        const double dNe      = dN + C * g_l * d_l;
        acc_l[k]             += -wf[e] * dNe / (fermion_.q[k] * fermion_.q[k] * fermion_.dq[k]);
        split_acc            += -2.0 * wf[e] * fermion_.q[k] * (dNe - dN);
        d_l_edge_[W * n + e]  = d_l;  // ApplyPerturbationOperator rebuilds its
                                      // coefficients from this (Lambda_k = Lambda + g_l d_l)
        dev_l_edge_[W * n + e] = f_l_st[k] - flg_st;
        dth_l_edge_[W * n + e] = LumpDTheta(g_l, wf, nf, jf, f_l_st, flg_st, e);
      }

      // Boson daughter: x2 = g_H g_l / g_phi dof ratio, folded into df_phi.
      const int kb     = k_boson_[n];
      constexpr int nb = DepositStencil::kWidth;
      const double* wb = &w_boson_[W * n];
      double fphig_st  = 0.;
      for (int e = 0; e < nb; ++e)
        fphig_st += wb[e] * f_phi_st[kb + e];
      const double th_p = LumpTheta(g_phi, wb, nb, kb, f_phi_st, fphig_st);
      theta_phi_[n]     = th_p;
      for (int e = 0; e < nb; ++e) {
        const int k               = kb + e;
        const double d_phi        = th_p * (f_phi_st[k] - fphig_st);
        const double dNe          = dN + C * g_phi * d_phi;
        acc_phi[k]               += -2.0 * wb[e] * dNe / (boson_.q[k] * boson_.q[k] * boson_.dq[k]);
        split_energy_residual_   += -2.0 * wb[e] * boson_.q[k] * (dNe - dN);
        d_phi_edge_[W * n + e]    = d_phi;
        dev_phi_edge_[W * n + e]  = f_phi_st[k] - fphig_st;
        dth_phi_edge_[W * n + e]  = LumpDTheta(g_phi, wb, nb, kb, f_phi_st, fphig_st, e);
      }
    }
  };

  for (int i = 0; i < parent_.n; ++i)
    sweep(i, df_bg_l_.data(), df_bg_phi_.data(), split_energy_residual_, clamped_energy_residual_);
}

void DecayTransitionKernel::PositivityDecomposition(PositivityTerms& fermion_out,
                                                    PositivityTerms& boson_out) const {
  fermion_out.Resize(fermion_.n);
  boson_out.Resize(boson_.n);
  for (int i = 0; i < parent_.n; ++i) {
    const double q1   = parent_.q[i];
    const double dq1  = parent_.dq[i];
    const double fH_i = fH_bg_[i];
    for (int n = node_off_[i]; n < node_off_[i] + node_cnt_[i]; ++n) {
      const double C      = k_coeff_ * dq1 * (q1 / eps1_[i]) * wn_[n];
      const double fl_g   = fl_gather_[n];
      const double fphi_g = fphi_gather_[n];
      // dLambda/df_l and dLambda/df_phi at this node (Lambda is linear in each).
      const double g_l   = (cfg_.inverse_decays ? fphi_g : 0.) +
                           (cfg_.quantum_statistics ? fH_i : 0.);
      const double g_phi = (cfg_.inverse_decays ? fl_g : 0.) -
                           (cfg_.quantum_statistics ? fH_i : 0.);
      // The rest of Lambda, i.e. what survives when that daughter's f is set to zero.
      // Under the limiter only the theta fraction of the f-linear part is diagonal;
      // the (1-theta) fraction stays on the gathered value and is therefore part of
      // the gain as far as THIS bin is concerned. theta = 1 recovers the fully
      // lumped decomposition exactly.
      const double th_l     = theta_l_[n];
      const double th_p     = theta_phi_[n];
      const double lam0_l   = lambda_s_[n] - th_l * g_l * fl_g;
      const double lam0_phi = lambda_s_[n] - th_p * g_phi * fphi_g;

      // Mirrors the deposit in PrepareTransitions exactly, INCLUDING whether the
      // f-linear part is lumped: the coupling to the neighbour is gone (offdiag == 0
      // identically) and gain_k reduces to Lambda0, which is -f_H (1 + f_phi) for the
      // fermion and -f_H (1 - f_l) for the boson. So a negative RHS-at-zero can only
      // mean f_H < 0 or f_l > 1 -- an upstream unphysical state, not a discretisation
      // defect. That is the point: after the fix this diagnostic stops reporting its
      // own scheme and starts reporting the inputs. offdiag is kept in the struct
      // because it is what the DEFECT looked like, and a reader comparing against the
      // commit that introduced this needs the slot to exist.
      constexpr int W        = DepositStencil::kWidth;
      const int jf           = j_fermion_[n];
      const double* const wf = &w_fermion_[W * n];
      for (int e = 0; e < DepositStencil::kWidth; ++e) {
        const int k          = jf + e;
        const double meas    = fermion_.q[k] * fermion_.q[k] * fermion_.dq[k];
        const double pref    = -1.0 * C * wf[e] / meas;
        fermion_out.gain[k] += pref * lam0_l;
        fermion_out.diag[k] += pref * th_l * g_l;
      }
      const int kb           = k_boson_[n];
      const double* const wb = &w_boson_[W * n];
      for (int e = 0; e < DepositStencil::kWidth; ++e) {
        const int k        = kb + e;
        const double meas  = boson_.q[k] * boson_.q[k] * boson_.dq[k];
        const double pref  = -2.0 * C * wb[e] / meas;
        boson_out.gain[k] += pref * lam0_phi;
        boson_out.diag[k] += pref * th_p * g_phi;
      }
    }
  }
}

void DecayTransitionKernel::ApplyPerturbationOperatorAllL(int l_max,
                                                          int stride,
                                                          const double* __restrict F_H,
                                                          const double* __restrict F_l,
                                                          const double* __restrict F_phi,
                                                          double* __restrict dF_H,
                                                          double* __restrict dF_l,
                                                          double* __restrict dF_phi) const {
  // Same contract as the per-l form, stated once for the whole sweep: this must follow a
  // PrepareTransitions, and it consumes every multipole, so nothing may have advanced the
  // recurrence first. API misuse, not a numeric condition, hence severe.
  class_test_severe(next_l_ != 0,
                    "ApplyPerturbationOperatorAllL must start from l=0 (next_l_=%d): it "
                    "consumes l=0..l_max in one pass and cannot resume mid-hierarchy",
                    next_l_);
  next_l_ = l_max + 1;

  for (int i = 0; i < parent_.n; ++i)
    for (int l = 0; l <= l_max; ++l)
      dF_H[(size_t) i * stride + l] = 0.;
  for (int j = 0; j < fermion_.n; ++j)
    for (int l = 0; l <= l_max; ++l)
      dF_l[(size_t) j * stride + l] = 0.;
  for (int k = 0; k < boson_.n; ++k)
    for (int l = 0; l <= l_max; ++l)
      dF_phi[(size_t) k * stride + l] = 0.;

  // Booked into locals and added once at the end, so the accumulation order within each
  // residual matches the per-l form's (ascending i, then n) even though l now interleaves.
  double resid_e = 0., resid_m = 0.;

  // Per-node, per-multipole scratch. On the STACK deliberately: it is written and read
  // within one node's block, so heap members would add nothing but an aliasing question
  // and a thread-safety one (this object is per-workspace, but stack arrays are per-call
  // whatever the ownership turns out to be). kMaxSweepL is a generous bound on
  // l_max_dncdm_col, which defaults to 17 and is capped by l_max_ncdm; exceeding it is a
  // configuration error, not a numeric condition, hence severe.
  constexpr int kMaxSweepL = 128;
  class_test_severe(l_max > kMaxSweepL,
                    "l_max = %d exceeds the sweep's scratch bound %d: raise kMaxSweepL in "
                    "ApplyPerturbationOperatorAllL",
                    l_max,
                    kMaxSweepL);
  double Pa_l[kMaxSweepL + 1], Pb_l[kMaxSweepL + 1], Pg_l[kMaxSweepL + 1];
  double Fl_node_l[kMaxSweepL + 1], Fphi_node_l[kMaxSweepL + 1];
  double dthg_l_l[kMaxSweepL + 1], dthg_p_l[kMaxSweepL + 1];

  for (int i = 0; i < parent_.n; ++i) {
    const double q1            = parent_.q[i];
    const double dq1           = parent_.dq[i];
    const double fH_i          = fH_bg_[i];
    const double base          = k_coeff_ * dq1 * (q1 / eps1_[i]);
    const double inv_meas_p    = 1.0 / (q1 * q1 * dq1);
    const double* const FH_bin = F_H + (size_t) i * stride;
    double* const dFH_bin      = dF_H + (size_t) i * stride;
    for (int n = node_off_[i]; n < node_off_[i] + node_cnt_[i]; ++n) {
      // Everything in this block is l-INDEPENDENT and was previously recomputed once per
      // multipole. It is the whole point of the inversion.
      const double B           = base * wn_[n];
      const double fl_g        = fl_gather_[n];
      const double fphi_g      = fphi_gather_[n];
      constexpr int W          = DepositStencil::kWidth;
      const int jf             = j_fermion_[n];
      constexpr int nf         = DepositStencil::kWidth;
      const double* const wf   = &w_fermion_[W * n];
      const int kb             = k_boson_[n];
      constexpr int nb         = DepositStencil::kWidth;
      const double* const wb   = &w_boson_[W * n];
      const double* const dle  = &d_l_edge_[W * n];
      const double* const dpe  = &d_phi_edge_[W * n];
      const double th_l        = theta_l_[n];
      const double th_p        = theta_phi_[n];
      const double one_m_thl   = 1. - th_l;
      const double one_m_thp   = 1. - th_p;
      const double* const devl = &dev_l_edge_[W * n];
      const double* const devp = &dev_phi_edge_[W * n];
      const double* const dthl = &dth_l_edge_[W * n];
      const double* const dthp = &dth_phi_edge_[W * n];
      const double ca = cos_a_[n], cb = cos_b_[n], cg = cos_g_[n];

      double c_H   = -1.0;
      double c_l   = cfg_.inverse_decays ? fphi_g : 0.0;
      double c_phi = cfg_.inverse_decays ? fl_g : 0.0;
      if (cfg_.quantum_statistics) {
        c_H   += fl_g - fphi_g;
        c_l   += fH_i;
        c_phi += -fH_i;
      }

      // Legendre factors for the whole hierarchy, up front. The recurrence is
      // sequential in l and so is a vectorisation barrier wherever it sits; evaluating
      // it ONCE into three small arrays (18 doubles each, L1-resident) buys every loop
      // below a stride-1 body the compiler can put in vector registers. Same recurrence,
      // same order, same values as the rolling pair it replaces.
      Pa_l[0] = 1.0;
      Pb_l[0] = 1.0;
      Pg_l[0] = 1.0;
      if (l_max >= 1) {
        Pa_l[1] = ca;
        Pb_l[1] = cb;
        Pg_l[1] = cg;
      }
      for (int l = 2; l <= l_max; ++l) {
        const double c_cur = 2.0 * l - 1.0;
        const double c_prv = l - 1;
        const double inv_l = 1.0 / l;
        Pa_l[l]            = (c_cur * ca * Pa_l[l - 1] - c_prv * Pa_l[l - 2]) * inv_l;
        Pb_l[l]            = (c_cur * cb * Pb_l[l - 1] - c_prv * Pb_l[l - 2]) * inv_l;
        Pg_l[l]            = (c_cur * cg * Pg_l[l - 1] - c_prv * Pg_l[l - 2]) * inv_l;
      }

      // Node gathers, as one AXPY per deposit edge over contiguous l. Summing over e
      // in the same ascending order as before keeps these bit-identical.
      for (int l = 0; l <= l_max; ++l) {
        Fl_node_l[l]   = 0.;
        Fphi_node_l[l] = 0.;
        dthg_l_l[l]    = 0.;
        dthg_p_l[l]    = 0.;
      }
      for (int e = 0; e < nf; ++e) {
        const double c                     = wf[e];
        const double* const __restrict src = F_l + (size_t) (jf + e) * stride;
        for (int l = 0; l <= l_max; ++l)
          Fl_node_l[l] += c * src[l];
      }
      for (int e = 0; e < nb; ++e) {
        const double c                     = wb[e];
        const double* const __restrict src = F_phi + (size_t) (kb + e) * stride;
        for (int l = 0; l <= l_max; ++l)
          Fphi_node_l[l] += c * src[l];
      }
      // delta-theta gathers: one per (node, l), NOT one per deposit edge -- theta is a
      // property of the whole stencil, so dthg carries no e index.
      for (int t = 0; t < nf; ++t) {
        const double c                     = dthl[t];
        const double* const __restrict src = F_l + (size_t) (jf + t) * stride;
        for (int l = 0; l <= l_max; ++l)
          dthg_l_l[l] += c * src[l];
      }
      for (int t = 0; t < nb; ++t) {
        const double c                     = dthp[t];
        const double* const __restrict src = F_phi + (size_t) (kb + t) * stride;
        for (int l = 0; l <= l_max; ++l)
          dthg_p_l[l] += c * src[l];
      }

      // Parent leg.
      for (int l = 0; l <= l_max; ++l) {
        const double sH    = c_H * FH_bin[l];
        const double sl    = c_l * Fl_node_l[l];
        const double sphi  = c_phi * Fphi_node_l[l];
        dFH_bin[l]        += B * (sH + sl * Pa_l[l] + sphi * Pb_l[l]) * inv_meas_p;
      }

      // Daughter scatters. Everything that depends only on the edge is hoisted, so the
      // l body is pure stride-1 arithmetic on doubles.
      for (int e = 0; e < nf; ++e) {
        const int k                        = jf + e;
        const double dl                    = dle[e];
        const double cH_e                  = c_H + (cfg_.quantum_statistics ? dl : 0.);
        const double cphi_e                = c_phi + (cfg_.inverse_decays ? dl : 0.);
        const double dv                    = devl[e];
        const double pref                  = -wf[e] * inv_meas_f_[k];
        const double* const __restrict src = F_l + (size_t) k * stride;
        double* const __restrict dst       = dF_l + (size_t) k * stride;
        for (int l = 0; l <= l_max; ++l) {
          // theta blend on the daughter's own leg, matching the background's lumped
          // loss (see LumpTheta); th_l = 1 gives the fully lumped expression. theta
          // is a function of the state, so perturbing it perturbs Lambda_k too, which
          // is the dthg term.
          const double Fl_e    = th_l * src[l] + one_m_thl * Fl_node_l[l] + dv * dthg_l_l[l];
          const double dN_l_k  = B * (cH_e * FH_bin[l] * Pa_l[l] + c_l * Fl_e +
                                      cphi_e * Fphi_node_l[l] * Pg_l[l]);
          dst[l]              += pref * dN_l_k;
        }
      }
      for (int e = 0; e < nb; ++e) {
        const int k                        = kb + e;
        const double dphi                  = dpe[e];
        const double cH_e                  = c_H - (cfg_.quantum_statistics ? dphi : 0.);
        const double cl_e                  = c_l + (cfg_.inverse_decays ? dphi : 0.);
        const double dv                    = devp[e];
        const double pref                  = -2.0 * wb[e] * inv_meas_b_[k];
        const double* const __restrict src = F_phi + (size_t) k * stride;
        double* const __restrict dst       = dF_phi + (size_t) k * stride;
        for (int l = 0; l <= l_max; ++l) {
          const double Fphi_e    = th_p * src[l] + one_m_thp * Fphi_node_l[l] + dv * dthg_p_l[l];
          const double dN_phi_k  = B * (cH_e * FH_bin[l] * Pb_l[l] + cl_e * Fl_node_l[l] * Pg_l[l] +
                                        c_phi * Fphi_e);
          dst[l]                += pref * dN_phi_k;
        }
      }

      // Booked split residuals, which exist only at l = 0 (energy) and l = 1 (momentum).
      // Recomputing dN_*_k for those two multipoles costs 2 of 18 and keeps the scatter
      // loops above branch-free; accumulating straight into resid_e/resid_m preserves
      // the summation order (fermion edges, then boson edges, ascending i then n).
      for (int l = 0; l <= 1 && l <= l_max; ++l) {
        const double sH     = c_H * FH_bin[l];
        const double sl     = c_l * Fl_node_l[l];
        const double sphi   = c_phi * Fphi_node_l[l];
        const double dN_l   = B * (sH * Pa_l[l] + sl + sphi * Pg_l[l]);
        const double dN_phi = B * (sH * Pb_l[l] + sl * Pg_l[l] + sphi);
        double& R           = (l == 0) ? resid_e : resid_m;
        for (int e = 0; e < nf; ++e) {
          const int k          = jf + e;
          const double dl      = dle[e];
          const double cH_e    = c_H + (cfg_.quantum_statistics ? dl : 0.);
          const double cphi_e  = c_phi + (cfg_.inverse_decays ? dl : 0.);
          const double Fl_e    = th_l * F_l[(size_t) k * stride + l] + one_m_thl * Fl_node_l[l] +
                                 devl[e] * dthg_l_l[l];
          const double dN_l_k  = B * (cH_e * FH_bin[l] * Pa_l[l] + c_l * Fl_e +
                                      cphi_e * Fphi_node_l[l] * Pg_l[l]);
          R                   += -2.0 * wf[e] * fermion_.q[k] * (dN_l_k - dN_l);
        }
        for (int e = 0; e < nb; ++e) {
          const int k         = kb + e;
          const double dphi   = dpe[e];
          const double cH_e   = c_H - (cfg_.quantum_statistics ? dphi : 0.);
          const double cl_e   = c_l + (cfg_.inverse_decays ? dphi : 0.);
          const double Fphi_e = th_p * F_phi[(size_t) k * stride + l] + one_m_thp * Fphi_node_l[l] +
                                devp[e] * dthg_p_l[l];
          const double dN_phi_k  = B * (cH_e * FH_bin[l] * Pb_l[l] + cl_e * Fl_node_l[l] * Pg_l[l] +
                                        c_phi * Fphi_e);
          R                     += -2.0 * wb[e] * boson_.q[k] * (dN_phi_k - dN_phi);
        }
      }
    }
  }
  split_energy_residual_pert_   += resid_e;
  split_momentum_residual_pert_ += resid_m;
}

void DecayTransitionKernel::CollisionDiagonal(double* diag_H,
                                              double* diag_l,
                                              double* diag_phi) const {
  for (int i = 0; i < parent_.n; ++i)
    diag_H[i] = 0.;
  for (int j = 0; j < fermion_.n; ++j)
    diag_l[j] = 0.;
  for (int k = 0; k < boson_.n; ++k)
    diag_phi[k] = 0.;

  for (int i = 0; i < parent_.n; ++i) {
    const double q1    = parent_.q[i];
    const double dq1   = parent_.dq[i];
    const double fH_i  = fH_bg_[i];
    const double measH = q1 * q1 * dq1;
    for (int n = node_off_[i]; n < node_off_[i] + node_cnt_[i]; ++n) {
      const double C      = k_coeff_ * dq1 * (q1 / eps1_[i]) * wn_[n];
      const double fl_g   = fl_gather_[n];
      const double fphi_g = fphi_gather_[n];

      // d(Lambda)/d(f_H). The decay term contributes -1 always; quantum statistics
      // adds (f_l - f_phi), i.e. Pauli blocking minus Bose enhancement. Physical
      // states keep this <= 0 (f_l <= 1, f_phi >= 0), which is what lets the
      // exponential step treat it as a relaxation rate rather than a growth rate.
      const double dLam_dfH  = -1.0 + (cfg_.quantum_statistics ? (fl_g - fphi_g) : 0.);
      diag_H[i]             += C * dLam_dfH / measH;

      // Daughters: same coefficients the lumped deposit uses in PrepareTransitions,
      // so g_l / g_phi here are d(Lambda)/d(f_l) and d(Lambda)/d(f_phi).
      const double g_l   = (cfg_.inverse_decays ? fphi_g : 0.) +
                           (cfg_.quantum_statistics ? fH_i : 0.);
      const double g_phi = (cfg_.inverse_decays ? fl_g : 0.) -
                           (cfg_.quantum_statistics ? fH_i : 0.);

      // d(Lambda_k)/d(f_k) = g (theta + (1-theta) dG/df_k), where G is the gathered
      // occupation. Fully lumped (theta = 1) this is just g -- which is why the
      // shipped diagonal is a bare g and is exact. Under the limiter the (1-theta)
      // fraction still rides on the gather, so its derivative has to be carried:
      // dG/df_k is the deposit weight for the linear gather and the log-odds chain
      // factor for the balanced one. Getting this term wrong is not cosmetic --
      // this diagonal is what the ETD evolver integrates analytically.
      // ...and theta is itself a function of the state under the limiter, so the
      // (f_k - G) d(theta)/df_k term rides along too. Dropping it is a 4.1e-2 error
      // in the diagonal (measured), which is not acceptable in the operator ETD
      // integrates analytically.
      const double th_l      = theta_l_[n];
      const double th_p      = theta_phi_[n];
      constexpr int W        = DepositStencil::kWidth;
      const int jf           = j_fermion_[n];
      const double* const wf = &w_fermion_[W * n];
      constexpr int nf       = DepositStencil::kWidth;
      for (int e = 0; e < nf; ++e) {
        const int k       = jf + e;
        const double meas = fermion_.q[k] * fermion_.q[k] * fermion_.dq[k];
        const double dGdf = wf[e];
        // Read the deviation and d(theta)/df that PrepareTransitions ALREADY stored.
        // Recomputing them here would have to reproduce its STATE-grid gather, and
        // the tempting `fl_gather_[n]` is the BACKGROUND-grid one -- a distinction
        // that vanishes at dr_bg_refine = 1 and silently corrupts everything above it
        // (the same trap the lumped-loss comment above records).
        const double dth   = dth_l_edge_[W * n + e];
        const double dLam  = th_l + (1. - th_l) * dGdf + dev_l_edge_[W * n + e] * dth;
        diag_l[k]         += -1.0 * C * wf[e] * g_l * dLam / meas;
      }
      const int kb           = k_boson_[n];
      const double* const wb = &w_boson_[W * n];
      constexpr int nb       = DepositStencil::kWidth;
      for (int e = 0; e < nb; ++e) {
        const int k        = kb + e;
        const double meas  = boson_.q[k] * boson_.q[k] * boson_.dq[k];
        const double dGdf  = wb[e];
        const double dth   = dth_phi_edge_[W * n + e];
        const double dLam  = th_p + (1. - th_p) * dGdf + dev_phi_edge_[W * n + e] * dth;
        diag_phi[k]       += -2.0 * C * wb[e] * g_phi * dLam / meas;
      }
    }
  }
}

void DecayTransitionKernel::ComputeBackgroundDerivs(double a,
                                                    double m,
                                                    double Gamma,
                                                    const double* f_H,
                                                    const double* f_l,
                                                    const double* f_phi,
                                                    double* df_H,
                                                    double* df_l,
                                                    double* df_phi,
                                                    double a_prime_over_a) const {
  PrepareTransitions(a, m, Gamma, f_H, f_l, f_phi, a_prime_over_a);
  for (int i = 0; i < parent_.n; ++i)
    df_H[i] = df_bg_H_[i];
  for (int j = 0; j < fermion_.n; ++j)
    df_l[j] = df_bg_l_[j];
  for (int k = 0; k < boson_.n; ++k)
    df_phi[k] = df_bg_phi_[k];
}

DecayTransitionKernel::Moments DecayTransitionKernel::ComputeMoments(
    double a, double m, const double* df_H, const double* df_l, const double* df_phi) const {
  const double A = a * a * m * m;
  Moments M{0., 0., 0., 0., 0., 0.};
  for (int i = 0; i < parent_.n; ++i) {
    const double q    = parent_.q[i];
    const double w    = parent_.dq[i] * q * q;
    const double eps  = std::sqrt(q * q + A);
    M.N_H            += w * df_H[i];
    M.E_H            += w * eps * df_H[i];
  }
  for (int j = 0; j < fermion_.n; ++j) {
    const double q  = fermion_.q[j];
    const double w  = fermion_.dq[j] * q * q;
    M.N_l          += w * df_l[j];
    M.E_l          += w * q * df_l[j];  // massless daughter: eps = q
  }
  for (int k = 0; k < boson_.n; ++k) {
    const double q  = boson_.q[k];
    const double w  = boson_.dq[k] * q * q;
    M.N_phi        += w * df_phi[k];
    M.E_phi        += w * q * df_phi[k];
  }
  return M;
}

double DecayTransitionKernel::LegendreP(int l, double x) {
  // Upward recurrence (2l+1) x P_l = (l+1) P_{l+1} + l P_{l-1}; never a table.
  if (l <= 0)
    return 1.0;
  double pm1 = 1.0;  // P_0
  double p   = x;    // P_1
  for (int n = 1; n < l; ++n) {
    const double pnext = ((2.0 * n + 1.0) * x * p - n * pm1) / (n + 1.0);
    pm1                = p;
    p                  = pnext;
  }
  return p;
}

void DecayTransitionKernel::ApplyPerturbationOperator(int l,
                                                      const double* F_H,
                                                      const double* F_l,
                                                      const double* F_phi,
                                                      double* dF_H,
                                                      double* dF_l,
                                                      double* dF_phi) const {
  // dF_* accumulate the UNNORMALIZED F-space collision rate (dF_i/dtau)^(1)_{C,l}.
  // No 1/f_i anywhere: F_i = delta-f_i is the natural variable for a species that can
  // fill from zero abundance (see the header). The band factor is linearized as
  //   dLambda = c_H F_H + c_l F_l + c_phi F_phi,
  // c_i = dLambda/df_i evaluated at the node (c_H=-1+qs(f_l-f_phi), c_l=f_phi+qs f_H,
  // c_phi=f_l-qs f_H; the inverse f_l f_phi term is dropped when inverse off). When
  // projecting onto species i's multipole l, species j enters with P_l(cos angle_ij)
  // via the Legendre addition theorem (angle_iH<->cos_a for l, cos_b for phi;
  // angle_lphi<->cos_g). The same transpose two-bin split as the background scatter
  // makes l=0 (N/E) and l=1 (momentum: q1=q2* cos_a+q3* cos_b) conserve per node.
  //
  // CALL-ORDER CONTRACT: l must advance by exactly 1 per call, starting from 0 after
  // each PrepareTransitions. The angular factors are carried in the per-node rolling
  // (P_{l-2}, P_{l-1}) state (Pa_/Pb_/Pg_prev_,curr_) and advanced ONE recurrence step
  // per call, so a skipped or repeated l would silently evaluate the wrong multipole.
  // This is an API-misuse check, not a numeric one (no cosmological parameter can
  // reach it), hence severe.
  class_test_severe(l != next_l_,
                    "ApplyPerturbationOperator called out of order: expected l=%d, got l=%d "
                    "(must be called for l=0,1,2,... in strict sequence after each "
                    "PrepareTransitions)",
                    next_l_,
                    l);
  next_l_ = l + 1;
  // Only l = 0 (energy) and l = 1 (momentum) carry a conservation identity, so only
  // those two book a residual. Resolved once here, not per node.
  double* const resid = (l == 0)   ? &split_energy_residual_pert_
                        : (l == 1) ? &split_momentum_residual_pert_
                                   : nullptr;

  for (int i = 0; i < parent_.n; ++i)
    dF_H[i] = 0.;
  for (int j = 0; j < fermion_.n; ++j)
    dF_l[j] = 0.;
  for (int k = 0; k < boson_.n; ++k)
    dF_phi[k] = 0.;

  for (int i = 0; i < parent_.n; ++i) {
    const double q1   = parent_.q[i];
    const double dq1  = parent_.dq[i];
    const double fH_i = fH_bg_[i];
    const double base = k_coeff_ * dq1 * (q1 / eps1_[i]);  // B_s prefactor / wn
    const double FH_i = F_H[i];
    for (int n = node_off_[i]; n < node_off_[i] + node_cnt_[i]; ++n) {
      const double B      = base * wn_[n];  // number-rate prefactor for this node
      const double fl_g   = fl_gather_[n];
      const double fphi_g = fphi_gather_[n];

      // Gather the daughter delta-f at the nodes with the SAME two-bin split as f
      // (linear interpolation of F itself, the transpose of the deposit below).
      constexpr int W        = DepositStencil::kWidth;
      const int jf           = j_fermion_[n];
      constexpr int nf       = DepositStencil::kWidth;
      const double* const wf = &w_fermion_[W * n];
      double Fl_node         = 0.;
      for (int e = 0; e < nf; ++e)
        Fl_node += wf[e] * F_l[jf + e];
      const int kb           = k_boson_[n];
      constexpr int nb       = DepositStencil::kWidth;
      const double* const wb = &w_boson_[W * n];
      double Fphi_node       = 0.;
      for (int e = 0; e < nb; ++e)
        Fphi_node += wb[e] * F_phi[kb + e];

      // Linearized band-factor coefficients c_i = dLambda/df_i at this node.
      double c_H   = -1.0;
      double c_l   = cfg_.inverse_decays ? fphi_g : 0.0;
      double c_phi = cfg_.inverse_decays ? fl_g : 0.0;
      if (cfg_.quantum_statistics) {
        c_H   += fl_g - fphi_g;
        c_l   += fH_i;
        c_phi += -fH_i;
      }
      // Per-species source amplitudes s_i = c_i F_i (F-space perturbation of Lambda).
      const double sH   = c_H * FH_i;
      const double sl   = c_l * Fl_node;
      const double sphi = c_phi * Fphi_node;

      // Angular factors advanced by ONE recurrence step from the cached pair instead of
      // rebuilt from P_0: same relation, same operand order as LegendreP, so the value is
      // LegendreP(l, x) to within the last bit. (Not bit-for-bit: freed of LegendreP's
      // serial inner loop this loop vectorizes, and -ffast-math then contracts/reassociates
      // it differently -- a <=1 ulp difference in P_l. The l=0/l=1 conservation identities
      // and the high-l P_l check in decay_kernel_test.cpp bound the consequence.)
      double Pa;  // angle H<->l
      double Pb;  // angle H<->phi
      double Pg;  // angle l<->phi
      if (l == 0) {
        Pa = 1.0;  // P_0, LegendreP's l<=0 branch
        Pb = 1.0;
        Pg = 1.0;
      }
      else if (l == 1) {
        Pa = cos_a_[n];  // P_1 = x, LegendreP's loop seed
        Pb = cos_b_[n];
        Pg = cos_g_[n];
      }
      else {
        // LegendreP's step P_{n+1} = ((2n+1) x P_n - n P_{n-1})/(n+1) at n = l-1.
        const double c_cur = 2.0 * l - 1.0;
        const double c_prv = l - 1;
        Pa                 = (c_cur * cos_a_[n] * Pa_curr_[n] - c_prv * Pa_prev_[n]) / l;
        Pb                 = (c_cur * cos_b_[n] * Pb_curr_[n] - c_prv * Pb_prev_[n]) / l;
        Pg                 = (c_cur * cos_g_[n] * Pg_curr_[n] - c_prv * Pg_prev_[n]) / l;
      }
      Pa_prev_[n] = Pa_curr_[n];
      Pa_curr_[n] = Pa;
      Pb_prev_[n] = Pb_curr_[n];
      Pb_curr_[n] = Pb;
      Pg_prev_[n] = Pg_curr_[n];
      Pg_curr_[n] = Pg;

      // Parent deposit (project onto H's multipole): {H:1, l:Pa, phi:Pb}.
      const double dLambda_H  = sH + sl * Pa + sphi * Pb;
      dF_H[i]                += B * dLambda_H / (q1 * q1 * dq1);

      // Fermion deposit (project onto l's multipole): {H:Pa, l:1, phi:Pg}, LUMPED --
      // the exact Jacobian of the lumped background deposit, so f-bar and delta-f
      // evolve under the same discrete model. For the deposit into bin k the band
      // factor is that bin's own Lambda_k = Lambda0_l + g_l f_l[k], so F_l is read at
      // k rather than gathered, and the coefficients multiplying F_H and F_phi carry
      // f_l[k] (precomputed in PrepareTransitions: they are background-only).
      //
      // The l=0 NUMBER identities survive exactly, by the same telescoping that keeps
      // the background's: sum_e w_e of each per-edge coefficient returns the gathered
      // one, so sum_e w_e dN_l_k == dN_l. The q-WEIGHTED identities (l=0 energy,
      // l=1 momentum) pick up the same O(dq) misplacement as the background's energy,
      // booked in split_energy_residual_pert_ / split_momentum_residual_pert_.
      const double dN_l   = B * (sH * Pa + sl + sphi * Pg);
      const double dN_phi = B * (sH * Pb + sl * Pg + sphi);
      const double th_l   = theta_l_[n];
      const double th_p   = theta_phi_[n];
      for (int e = 0; e < nf; ++e) {
        const int k         = jf + e;
        const double dl     = d_l_edge_[W * n + e];
        const double cH_e   = c_H + (cfg_.quantum_statistics ? dl : 0.);
        const double cphi_e = c_phi + (cfg_.inverse_decays ? dl : 0.);
        // The daughter's OWN leg carries the same theta blend as the background's
        // lumped loss (see LumpTheta): theta F_l[k] + (1-theta) Fl_node. At
        // theta = 1 this is the shipped c_l F_l[k] exactly. Blending here rather
        // than only in d_l_edge_ is what keeps this the true Jacobian.
        double dthg_l = 0.;
        for (int t = 0; t < nf; ++t)
          dthg_l += dth_l_edge_[W * n + t] * F_l[jf + t];
        const double Fl_e = th_l * F_l[k] + (1. - th_l) * Fl_node + dev_l_edge_[W * n + e] * dthg_l;
        const double dN_l_k  = B * (cH_e * FH_i * Pa + c_l * Fl_e + cphi_e * Fphi_node * Pg);
        dF_l[k]             += -wf[e] * dN_l_k / (fermion_.q[k] * fermion_.q[k] * fermion_.dq[k]);
        if (resid)
          *resid += -2.0 * wf[e] * fermion_.q[k] * (dN_l_k - dN_l);
      }

      // Boson deposit (project onto phi's multipole): {H:Pb, l:Pg, phi:1}; x2 dof.
      for (int e = 0; e < nb; ++e) {
        const int k       = kb + e;
        const double dphi = d_phi_edge_[W * n + e];
        const double cH_e = c_H - (cfg_.quantum_statistics ? dphi : 0.);
        const double cl_e = c_l + (cfg_.inverse_decays ? dphi : 0.);
        double dthg_p     = 0.;
        for (int t = 0; t < nb; ++t)
          dthg_p += dth_phi_edge_[W * n + t] * F_phi[kb + t];
        const double Fphi_e   = th_p * F_phi[k] + (1. - th_p) * Fphi_node +
                                dev_phi_edge_[W * n + e] * dthg_p;
        const double dN_phi_k = B * (cH_e * FH_i * Pb + cl_e * Fl_node * Pg + c_phi * Fphi_e);
        dF_phi[k] += -2.0 * wb[e] * dN_phi_k / (boson_.q[k] * boson_.q[k] * boson_.dq[k]);
        if (resid)
          *resid += -2.0 * wb[e] * boson_.q[k] * (dN_phi_k - dN_phi);
      }
    }
  }
}

void DecayTransitionKernel::NodeKinematics(
    int i, int s, double& q2, double& q3, double& cos_a, double& cos_b, double& cos_g) const {
  const int n = node_off_[i] + s;
  q2          = q2_star_[n];
  q3          = q3_star_[n];
  cos_a       = cos_a_[n];
  cos_b       = cos_b_[n];
  cos_g       = cos_g_[n];
}
