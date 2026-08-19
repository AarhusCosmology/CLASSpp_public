#include "reduced_collision_operator.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>

namespace {

// ---------------------------------------------------------------- small LA --
// Everything here is O(n^3) on matrices of at most a few hundred rows, run once per
// background state, never on a hot path. CLASSpp has no BLAS/LAPACK, so these are
// written out; they are deliberately the simplest correct versions.

/** Cholesky S = L L^T, lower triangular, in place into `L` (row-major n x n).
 *  Throws if S is not positive definite -- which for the metric below means a
 *  daughter basis whose Gram matrix is singular, i.e. a weight with no support. */
void Cholesky(const std::vector<double>& S, int n, std::vector<double>& L) {
  L.assign((size_t) n * n, 0.);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j <= i; ++j) {
      double sum = S[(size_t) i * n + j];
      for (int k = 0; k < j; ++k)
        sum -= L[(size_t) i * n + k] * L[(size_t) j * n + k];
      if (i == j) {
        if (!(sum > 0.)) {
          char msg[192];
          std::snprintf(msg,
                        sizeof msg,
                        "ReducedCollisionOperator: matrix not positive definite at pivot %d/%d "
                        "(pivot = %.6e, diagonal = %.6e)",
                        i,
                        n,
                        sum,
                        S[(size_t) i * n + i]);
          throw std::runtime_error(msg);
        }
        L[(size_t) i * n + i] = std::sqrt(sum);
      }
      else {
        L[(size_t) i * n + j] = sum / L[(size_t) j * n + j];
      }
    }
  }
}

/** Solve L z = b for lower-triangular L (row-major n x n), in place on `b`. */
void ForwardSolve(const std::vector<double>& L, int n, double* b) {
  for (int i = 0; i < n; ++i) {
    double sum = b[i];
    for (int k = 0; k < i; ++k)
      sum -= L[(size_t) i * n + k] * b[k];
    b[i] = sum / L[(size_t) i * n + i];
  }
}

/** Largest eigenvalue of a real symmetric matrix by cyclic Jacobi rotations.
 *  `A` is row-major n x n and is destroyed. */
double MaxEigenvalueSymmetric(std::vector<double>& A, int n) {
  const int kMaxSweeps = 60;
  for (int sweep = 0; sweep < kMaxSweeps; ++sweep) {
    double off = 0.;
    for (int p = 0; p < n; ++p)
      for (int q = p + 1; q < n; ++q)
        off += A[(size_t) p * n + q] * A[(size_t) p * n + q];
    if (off <= 1e-30)
      break;
    for (int p = 0; p < n; ++p) {
      for (int q = p + 1; q < n; ++q) {
        const double apq = A[(size_t) p * n + q];
        if (std::fabs(apq) < 1e-300)
          continue;
        const double app   = A[(size_t) p * n + p];
        const double aqq   = A[(size_t) q * n + q];
        const double theta = 0.5 * (aqq - app) / apq;
        const double t     = (theta >= 0. ? 1. : -1.) /
                             (std::fabs(theta) + std::sqrt(theta * theta + 1.));
        const double c     = 1. / std::sqrt(t * t + 1.);
        const double s     = t * c;
        for (int k = 0; k < n; ++k) {
          const double akp      = A[(size_t) k * n + p];
          const double akq      = A[(size_t) k * n + q];
          A[(size_t) k * n + p] = c * akp - s * akq;
          A[(size_t) k * n + q] = s * akp + c * akq;
        }
        for (int k = 0; k < n; ++k) {
          const double apk      = A[(size_t) p * n + k];
          const double aqk      = A[(size_t) q * n + k];
          A[(size_t) p * n + k] = c * apk - s * aqk;
          A[(size_t) q * n + k] = s * apk + c * aqk;
        }
      }
    }
  }
  double lam = -1e300;
  for (int i = 0; i < n; ++i)
    lam = std::fmax(lam, A[(size_t) i * n + i]);
  return lam;
}

/** max Re lambda(M) <= largest eigenvalue of the S-symmetric part of M.
 *
 *  With S = L L^T the S-inner-product numerical range of M equals the ordinary
 *  numerical range of Mhat = L^T M L^-T, because for y = L^T x
 *      x^T S M x / x^T S x  =  y^T (L^T M L^-T) y / y^T y.
 *  So form Mhat, symmetrise, and take its largest eigenvalue. */
std::vector<double> MetricTransform(const double* M, const std::vector<double>& S, int n) {
  std::vector<double> L;
  Cholesky(S, n, L);

  // A = L^T M   (A_ij = sum_k L_ki M_kj, and L is lower so k >= i)
  std::vector<double> A((size_t) n * n, 0.);
  for (int i = 0; i < n; ++i)
    for (int k = i; k < n; ++k) {
      const double lki = L[(size_t) k * n + i];
      if (lki == 0.)
        continue;
      for (int j = 0; j < n; ++j)
        A[(size_t) i * n + j] += lki * M[(size_t) k * n + j];
    }

  // Mhat = A L^-T, i.e. Mhat L^T = A, i.e. L Mhat^T = A^T: forward-solve each row.
  std::vector<double> Mhat((size_t) n * n, 0.);
  std::vector<double> row(n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j)
      row[j] = A[(size_t) i * n + j];
    ForwardSolve(L, n, row.data());
    for (int j = 0; j < n; ++j)
      Mhat[(size_t) i * n + j] = row[j];
  }
  return Mhat;
}

double NumericalRangeBound(const double* M, const std::vector<double>& S, int n) {
  std::vector<double> Mhat = MetricTransform(M, S, n);
  std::vector<double> Sym((size_t) n * n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      Sym[(size_t) i * n + j] = 0.5 * (Mhat[(size_t) i * n + j] + Mhat[(size_t) j * n + i]);
  return MaxEigenvalueSymmetric(Sym, n);
}

/** f(1-f) for a fermion, f(1+f) for a boson -- and it MUST come out non-negative.
 *
 *  The discrete gather does not enforce Pauli blocking, so a fermion bin can overshoot
 *  saturation (the design note measures a peak f_l ~ 1.01, and DrPsdSpecies floors f at
 *  kFFloor but never caps it at 1). f(1-f) is then negative, and nothing downstream
 *  notices: dmu = dq q^2 w becomes a SIGNED measure, the Stieltjes recurrence
 *  orthogonalises in an indefinite form, and Reconstruct returns a df whose sign is
 *  flipped across the saturated core. BuildMetric's Cholesky would catch the indefinite
 *  metric, but the shipped species path only calls Reconstruct/Project and never builds
 *  one.
 *
 *  Saturation is the physical limit -- a bin at f = 1 has no capacity left -- so the
 *  overshoot is clamped to it, at the kernel's own kFMax so this operator's view of the
 *  background is the same one the kernel already uses (EtaOf/ClampF).
 *
 *  The lower end is deliberately NOT floored to kFMin: w = 0 on an unpopulated bin is
 *  the correct weight, and it is what makes the "basis weight has no support" guard in
 *  SetBackground fire on an empty daughter instead of silently building a basis out of
 *  1e-100. */
double EntropyWeight(double f, Statistics s) {
  const double x = (s == Statistics::Fermion) ? std::fmin(f, kFMax) : f;
  if (!(x > 0.))
    return 0.;
  return (s == Statistics::Fermion) ? x * (1. - x) : x * (1. + x);
}

/** Relative floor on the entropy weight when it appears in a DENOMINATOR.
 *
 *  A bin with w = 0 holds no thermodynamic state, so any perturbation placed there has
 *  infinite entropy norm. The exact operator's state vector does carry such bins (a
 *  daughter above the injection front); the reduced one never can. The floor is what
 *  makes the two numerical ranges comparable at all, so it MUST be the same constant on
 *  both sides -- a tighter floor on one side alone would silently break the containment
 *  the whole design rests on. */
constexpr double kWeightFloor = 1e-8;

double FlooredWeight(const std::vector<double>& w, size_t k) {
  double wmax = 0.;
  for (double x : w)
    wmax = std::fmax(wmax, x);
  const double fl = (wmax > 0.) ? kWeightFloor * wmax : 1.;
  return std::fmax(w[k], fl);
}

}  // namespace

// ------------------------------------------------------------- constructor --

ReducedCollisionOperator::ReducedCollisionOperator(const DecayTransitionKernel& kernel,
                                                   GridView parent,
                                                   GridView fermion,
                                                   GridView boson,
                                                   Statistics fermion_stat,
                                                   Statistics boson_stat,
                                                   Config cfg)
    : kernel_(&kernel), parent_(parent), fermion_(fermion), boson_(boson),
      fermion_stat_(fermion_stat), boson_stat_(boson_stat), cfg_(cfg) {
  if (cfg_.n_moments < 1)
    throw std::runtime_error("ReducedCollisionOperator: n_moments must be >= 1");
}

// -------------------------------------------------------------- the basis --

void ReducedCollisionOperator::SetBackground(
    double a, double m, const double* f_H, const double* f_l, const double* f_phi) {
  a_               = a;
  m_               = m;
  const double am2 = a * a * m * m;

  eps_.resize(parent_.n);
  w_H_.resize(parent_.n);
  for (int b = 0; b < parent_.n; ++b) {
    eps_[b] = std::sqrt(parent_.q[b] * parent_.q[b] + am2);
    w_H_[b] = EntropyWeight(f_H[b], Statistics::Fermion);
  }

  const GridView* g[2]     = {&fermion_, &boson_};
  const double* f[2]       = {f_l, f_phi};
  const Statistics stat[2] = {fermion_stat_, boson_stat_};
  const int nm             = cfg_.n_moments;

  for (int d = 0; d < 2; ++d) {
    Basis& B = basis_[d];
    B.w.resize(g[d]->n);
    B.clamp = WeightClamp{};
    for (int k = 0; k < g[d]->n; ++k) {
      // Report before weighting: what matters is that the BACKGROUND went outside the
      // range an occupation can occupy, not what this class then did about it.
      B.clamp.f_max = std::fmax(B.clamp.f_max, f[d][k]);
      if (f[d][k] < 0. || (stat[d] == Statistics::Fermion && f[d][k] > kFMax))
        ++B.clamp.n_clamped;
      B.w[k] = EntropyWeight(f[d][k], stat[d]);
    }

    // q_ref is the weight's own mean momentum, so the monomials (q/q_ref)^j stay
    // O(1) across the basis. The Gram matrix is a Hilbert-like moment matrix and its
    // conditioning is the only thing limiting how many moments can be carried; with
    // this scaling n_moments = 4 is comfortable and 6 is where it starts to bite.
    double s2 = 0., s3 = 0.;
    for (int k = 0; k < g[d]->n; ++k) {
      const double base  = g[d]->dq[k] * g[d]->q[k] * g[d]->q[k] * B.w[k];
      s2                += base;
      s3                += base * g[d]->q[k];
    }
    B.q_ref = (s2 > 0. && s3 > 0.) ? s3 / s2 : 1.;

    // Stieltjes: orthogonal polynomials in the discrete measure dmu_k = dq q^2 w.
    //     psi_0     = 1
    //     psi_1     = (t - alpha_0) psi_0
    //     psi_{j+1} = (t - alpha_j) psi_j - beta_j psi_{j-1}
    // with alpha_j = <t psi_j, psi_j>/<psi_j, psi_j> and beta_j = n_j / n_{j-1}.
    // The recurrence never forms t^j, which is exactly why it survives where the
    // monomial Gram does not.
    const int ng = g[d]->n;
    std::vector<double> dmu(ng), t(ng);
    for (int k = 0; k < ng; ++k) {
      dmu[k] = g[d]->dq[k] * g[d]->q[k] * g[d]->q[k] * B.w[k];
      t[k]   = g[d]->q[k] / B.q_ref;
    }
    B.psi.assign((size_t) nm * ng, 0.);
    B.gnorm.assign(nm, 0.);
    for (int k = 0; k < ng; ++k)
      B.psi[k] = 1.;
    for (int k = 0; k < ng; ++k)
      B.gnorm[0] += dmu[k];
    if (!(B.gnorm[0] > 0.))
      throw std::runtime_error("ReducedCollisionOperator: basis weight has no support");

    for (int j = 0; j + 1 < nm; ++j) {
      double num = 0.;
      for (int k = 0; k < ng; ++k)
        num += dmu[k] * t[k] * B.psi[(size_t) j * ng + k] * B.psi[(size_t) j * ng + k];
      const double alpha = num / B.gnorm[j];
      const double beta  = (j == 0) ? 0. : B.gnorm[j] / B.gnorm[j - 1];
      if (j == 0)
        B.alpha0 = alpha;
      for (int k = 0; k < ng; ++k) {
        const double prev                = (j == 0) ? 0. : B.psi[(size_t) (j - 1) * ng + k];
        B.psi[(size_t) (j + 1) * ng + k] = (t[k] - alpha) * B.psi[(size_t) j * ng + k] -
                                           beta * prev;
      }
      double n2 = 0.;
      for (int k = 0; k < ng; ++k)
        n2 += dmu[k] * B.psi[(size_t) (j + 1) * ng + k] * B.psi[(size_t) (j + 1) * ng + k];
      // The measure lives on `ng` points, so it supports at most `ng` orthogonal
      // polynomials -- and in practice far fewer, because the weight is concentrated.
      // A collapsing norm is the basis genuinely running out, not a numerical wobble,
      // and it must be reported rather than divided by.
      if (!(n2 > 1e-24 * B.gnorm[0])) {
        char msg[160];
        std::snprintf(msg,
                      sizeof msg,
                      "ReducedCollisionOperator: basis exhausted at n_moments = %d "
                      "(||psi_%d||^2 / ||psi_0||^2 = %.3e)",
                      j + 2,
                      j + 1,
                      n2 / B.gnorm[0]);
        throw std::runtime_error(msg);
      }
      B.gnorm[j + 1] = n2;
    }

    // The entropy norm induced on the moment coordinates is G^-1 H G^-1 with
    // H_jj' = Sum dq q^2 (w^2/w) t^(j+j'). The basis weight IS the entropy weight, so
    // H == G, and with G diagonal the metric collapses to diag(1/gnorm).
    B.metric.assign((size_t) nm * nm, 0.);
    for (int j = 0; j < nm; ++j)
      B.metric[(size_t) j * nm + j] = 1. / B.gnorm[j];
  }
}

void ReducedCollisionOperator::AdoptFrozenBasis(double a,
                                                double m,
                                                const double* f_H,
                                                const ReducedCollisionOperator& src) {
  if (src.cfg_.n_moments != cfg_.n_moments)
    throw std::runtime_error("ReducedCollisionOperator: AdoptFrozenBasis moment-count mismatch");
  a_ = a;
  m_ = m;
  // The parent quantities DO belong to this row -- eps enters the conservation weights
  // and w_H the metric, and both are functions of a. Only the daughter BASIS is frozen.
  const double am2 = a * a * m * m;
  eps_.resize(parent_.n);
  w_H_.resize(parent_.n);
  for (int b = 0; b < parent_.n; ++b) {
    eps_[b] = std::sqrt(parent_.q[b] * parent_.q[b] + am2);
    w_H_[b] = EntropyWeight(f_H[b], Statistics::Fermion);
  }
  basis_[0] = src.basis_[0];
  basis_[1] = src.basis_[1];
}

void ReducedCollisionOperator::Reconstruct(int daughter, const double* coeff, double* df) const {
  const Basis& B    = basis_[daughter];
  const GridView& g = (daughter == 0) ? fermion_ : boson_;
  const int nm      = cfg_.n_moments;

  // c = Ginv * coeff, and G is diagonal, so Project(Reconstruct(coeff)) == coeff by a
  // scaling rather than a linear solve.
  for (int k = 0; k < g.n; ++k) {
    double acc = 0.;
    for (int j = 0; j < nm; ++j)
      acc += coeff[j] / B.gnorm[j] * B.psi[(size_t) j * g.n + k];
    df[k] = B.w[k] * acc;
  }
}

void ReducedCollisionOperator::Project(int daughter, const double* df, double* coeff) const {
  const Basis& B    = basis_[daughter];
  const GridView& g = (daughter == 0) ? fermion_ : boson_;
  const int nm      = cfg_.n_moments;

  for (int j = 0; j < nm; ++j)
    coeff[j] = 0.;
  for (int k = 0; k < g.n; ++k) {
    const double base = g.dq[k] * g.q[k] * g.q[k] * df[k];
    for (int j = 0; j < nm; ++j)
      coeff[j] += base * B.psi[(size_t) j * g.n + k];
  }
}

// ----------------------------------------------------------- the assembly --

void ReducedCollisionOperator::Assemble(int l_max, std::vector<double>& out) const {
  const int S      = size();
  const int nm     = cfg_.n_moments;
  const int Np     = parent_.n;
  const int stride = l_max + 1;
  out.assign((size_t) (l_max + 1) * S * S, 0.);

  std::vector<double> FH((size_t) Np * stride), Fl((size_t) fermion_.n * stride),
      Fp((size_t) boson_.n * stride);
  std::vector<double> dH((size_t) Np * stride), dl((size_t) fermion_.n * stride),
      dp((size_t) boson_.n * stride);
  std::vector<double> recon(std::max(fermion_.n, boson_.n));
  std::vector<double> slice(std::max(fermion_.n, boson_.n)), mom(nm), coeff(nm);

  for (int col = 0; col < S; ++col) {
    std::fill(FH.begin(), FH.end(), 0.);
    std::fill(Fl.begin(), Fl.end(), 0.);
    std::fill(Fp.begin(), Fp.end(), 0.);

    // The operator is block diagonal in l, so placing the SAME basis vector in every
    // multipole slot returns column `col` of M_l for every l in a single sweep --
    // which is what makes the whole assembly cost ~S operator applications rather
    // than S*(l_max+1).
    if (col < Np) {
      for (int l = 0; l <= l_max; ++l)
        FH[(size_t) col * stride + l] = 1.;
    }
    else {
      const int d = (col - Np) < nm ? 0 : 1;
      const int j = (col - Np) - d * nm;
      std::fill(coeff.begin(), coeff.end(), 0.);
      coeff[j] = 1.;
      Reconstruct(d, coeff.data(), recon.data());
      double* F        = (d == 0) ? Fl.data() : Fp.data();
      const int n_grid = (d == 0) ? fermion_.n : boson_.n;
      for (int k = 0; k < n_grid; ++k)
        for (int l = 0; l <= l_max; ++l)
          F[(size_t) k * stride + l] = recon[k];
    }

    kernel_->ResetMultipoleRecurrence();
    kernel_->ApplyPerturbationOperatorAllL(l_max,
                                           stride,
                                           FH.data(),
                                           Fl.data(),
                                           Fp.data(),
                                           dH.data(),
                                           dl.data(),
                                           dp.data());

    for (int l = 0; l <= l_max; ++l) {
      double* Ml = out.data() + (size_t) l * S * S;
      for (int i = 0; i < Np; ++i)
        Ml[(size_t) i * S + col] = dH[(size_t) i * stride + l];
      for (int d = 0; d < 2; ++d) {
        const double* src = (d == 0) ? dl.data() : dp.data();
        const int n_grid  = (d == 0) ? fermion_.n : boson_.n;
        for (int k = 0; k < n_grid; ++k)
          slice[k] = src[(size_t) k * stride + l];
        Project(d, slice.data(), mom.data());
        for (int j = 0; j < nm; ++j)
          Ml[(size_t) (Np + d * nm + j) * S + col] = mom[j];
      }
    }
  }
}

void ReducedCollisionOperator::AssembleExact(int l_max, std::vector<double>& out) const {
  const int N      = full_size();
  const int Np     = parent_.n;
  const int Nf     = fermion_.n;
  const int Nb     = boson_.n;
  const int stride = l_max + 1;
  out.assign((size_t) (l_max + 1) * N * N, 0.);

  std::vector<double> FH((size_t) Np * stride), Fl((size_t) Nf * stride), Fp((size_t) Nb * stride);
  std::vector<double> dH((size_t) Np * stride), dl((size_t) Nf * stride), dp((size_t) Nb * stride);

  for (int col = 0; col < N; ++col) {
    std::fill(FH.begin(), FH.end(), 0.);
    std::fill(Fl.begin(), Fl.end(), 0.);
    std::fill(Fp.begin(), Fp.end(), 0.);
    double* F     = (col < Np) ? FH.data() : (col < Np + Nf ? Fl.data() : Fp.data());
    const int idx = (col < Np) ? col : (col < Np + Nf ? col - Np : col - Np - Nf);
    for (int l = 0; l <= l_max; ++l)
      F[(size_t) idx * stride + l] = 1.;

    kernel_->ResetMultipoleRecurrence();
    kernel_->ApplyPerturbationOperatorAllL(l_max,
                                           stride,
                                           FH.data(),
                                           Fl.data(),
                                           Fp.data(),
                                           dH.data(),
                                           dl.data(),
                                           dp.data());

    for (int l = 0; l <= l_max; ++l) {
      double* Ml = out.data() + (size_t) l * N * N;
      for (int i = 0; i < Np; ++i)
        Ml[(size_t) i * N + col] = dH[(size_t) i * stride + l];
      for (int i = 0; i < Nf; ++i)
        Ml[(size_t) (Np + i) * N + col] = dl[(size_t) i * stride + l];
      for (int i = 0; i < Nb; ++i)
        Ml[(size_t) (Np + Nf + i) * N + col] = dp[(size_t) i * stride + l];
    }
  }
}

// ---------------------------------------------------------- conservation --

namespace {

/** ||w^T M|| / (||w|| * ||M||_max-col), the scale-free size of a broken identity. */
double LeftNullResidual(const std::vector<double>& w, const double* M, int n) {
  double wnorm = 0.;
  for (double x : w)
    wnorm += x * x;
  wnorm = std::sqrt(wnorm);
  if (!(wnorm > 0.))
    return 0.;

  double res = 0., scale = 0.;
  for (int j = 0; j < n; ++j) {
    double s = 0., col = 0.;
    for (int i = 0; i < n; ++i) {
      s   += w[i] * M[(size_t) i * n + j];
      col += w[i] * w[i] * M[(size_t) i * n + j] * M[(size_t) i * n + j];
    }
    res   = std::fmax(res, std::fabs(s));
    scale = std::fmax(scale, std::sqrt(col));
  }
  return (scale > 0.) ? res / scale : 0.;
}

}  // namespace

ReducedCollisionOperator::ConservationResidual ReducedCollisionOperator::Conservation(
    const double* M_l) const {
  const int S  = size();
  const int Np = parent_.n;
  const int nm = cfg_.n_moments;
  ConservationResidual out;

  auto parent_weight = [&](std::vector<double>& w, double g, int power) {
    for (int b = 0; b < Np; ++b) {
      const double base = g * parent_.dq[b] * parent_.q[b] * parent_.q[b];
      w[b]              = (power == 0) ? base : (power == 1 ? base * eps_[b] : base * parent_.q[b]);
    }
  };

  // N_H + N_l = 0 : the fermion-leg number identity. m_0 IS the number moment.
  std::vector<double> w(S, 0.);
  parent_weight(w, 1., 0);
  w[Np]              = 1.;
  out.number_fermion = LeftNullResidual(w, M_l, S);

  // 2 N_H + N_phi = 0 : the boson leg carries the folded 2 dof.
  std::fill(w.begin(), w.end(), 0.);
  parent_weight(w, 2., 0);
  w[Np + nm]       = 1.;
  out.number_boson = LeftNullResidual(w, M_l, S);

  if (nm >= 2) {
    // 2 E_H + 2 E_l + E_phi = 0, with E_i = Sum dq q^2 eps_i df_i and eps = q for the
    // massless daughters. The test functions are orthogonal polynomials, not monomials,
    // so t = psi_1 + alpha0 psi_0 and the daughter energy is q_ref (m_1 + alpha0 m_0)
    // rather than q_ref m_1. Getting this wrong would report a broken identity on an
    // operator that conserves perfectly well.
    auto daughter_energy = [&](std::vector<double>& v, int d, double deg) {
      const int base  = Np + d * nm;
      v[base]        += deg * basis_[d].q_ref * basis_[d].alpha0;
      v[base + 1]    += deg * basis_[d].q_ref;
    };
    std::fill(w.begin(), w.end(), 0.);
    parent_weight(w, 2., 1);
    daughter_energy(w, 0, 2.);
    daughter_energy(w, 1, 1.);
    out.energy = LeftNullResidual(w, M_l, S);

    // Momentum, the l = 1 channel: the parent's weight is q rather than eps, the
    // daughters' is unchanged. This is the identity whose discrete violation owns the
    // spurious positive eigenvalue of the shipped operator.
    std::fill(w.begin(), w.end(), 0.);
    parent_weight(w, 2., 2);
    daughter_energy(w, 0, 2.);
    daughter_energy(w, 1, 1.);
    out.momentum = LeftNullResidual(w, M_l, S);
  }
  return out;
}

ReducedCollisionOperator::ConservationResidual ReducedCollisionOperator::ConservationExact(
    const double* M_l) const {
  const int N  = full_size();
  const int Np = parent_.n;
  const int Nf = fermion_.n;
  ConservationResidual out;

  auto fill = [&](std::vector<double>& w, double gH, int power, double gf, double gb) {
    std::fill(w.begin(), w.end(), 0.);
    for (int b = 0; b < Np; ++b) {
      const double base = gH * parent_.dq[b] * parent_.q[b] * parent_.q[b];
      w[b]              = (power == 0) ? base : (power == 1 ? base * eps_[b] : base * parent_.q[b]);
    }
    for (int k = 0; k < Nf; ++k) {
      const double base = gf * fermion_.dq[k] * fermion_.q[k] * fermion_.q[k];
      w[Np + k]         = (power == 0) ? base : base * fermion_.q[k];
    }
    for (int k = 0; k < boson_.n; ++k) {
      const double base = gb * boson_.dq[k] * boson_.q[k] * boson_.q[k];
      w[Np + Nf + k]    = (power == 0) ? base : base * boson_.q[k];
    }
  };

  std::vector<double> w(N, 0.);
  fill(w, 1., 0, 1., 0.);
  out.number_fermion = LeftNullResidual(w, M_l, N);
  fill(w, 2., 0, 0., 1.);
  out.number_boson = LeftNullResidual(w, M_l, N);
  fill(w, 2., 1, 2., 1.);
  out.energy = LeftNullResidual(w, M_l, N);
  fill(w, 2., 2, 2., 1.);
  out.momentum = LeftNullResidual(w, M_l, N);
  return out;
}

// ---------------------------------------------------------- dissipativity --

void ReducedCollisionOperator::BuildMetric(std::vector<double>& S) const {
  const int n  = size();
  const int Np = parent_.n;
  const int nm = cfg_.n_moments;
  S.assign((size_t) n * n, 0.);

  // The parent is q-resolved, so its metric is the plain entropy norm. Bins where the
  // parent has decayed away carry w_H -> 0 and would give an unbounded metric; they
  // are floored, which is a diagnostic compromise and not a proof (documented on
  // DissipativityBound). The reduced daughter blocks need no such floor: the
  // reconstruction only ever produces df proportional to w, so the Gram matrix is
  // exactly the induced metric and is finite by construction.
  for (int b = 0; b < Np; ++b)
    S[(size_t) b * n + b] = 2. * parent_.dq[b] * parent_.q[b] * parent_.q[b] /
                            FlooredWeight(w_H_, (size_t) b);

  const double g[2] = {2., 1.};
  for (int d = 0; d < 2; ++d)
    for (int i = 0; i < nm; ++i)
      for (int j = 0; j < nm; ++j)
        S[(size_t) (Np + d * nm + i) * n + (Np + d * nm + j)] =
            g[d] * basis_[d].metric[(size_t) i * nm + j];
}

void ReducedCollisionOperator::BuildMetricExact(std::vector<double>& S) const {
  const int n  = full_size();
  const int Np = parent_.n;
  const int Nf = fermion_.n;
  S.assign((size_t) n * n, 0.);

  auto block = [&](int base, const GridView& g, const std::vector<double>& w, double deg) {
    for (int k = 0; k < g.n; ++k)
      S[(size_t) (base + k) * n + (base + k)] = deg * g.dq[k] * g.q[k] * g.q[k] /
                                                FlooredWeight(w, (size_t) k);
  };
  block(0, parent_, w_H_, 2.);
  block(Np, fermion_, basis_[0].w, 2.);
  block(Np + Nf, boson_, basis_[1].w, 1.);
}

double ReducedCollisionOperator::DissipativityBound(const double* M_l) const {
  std::vector<double> S;
  BuildMetric(S);
  return NumericalRangeBound(M_l, S, size());
}

double ReducedCollisionOperator::DissipativityBoundExact(const double* M_l) const {
  std::vector<double> S;
  BuildMetricExact(S);
  return NumericalRangeBound(M_l, S, full_size());
}

// ------------------------------------------------------------ equilibrium --

void ReducedCollisionOperator::EquilibriumNullVectors(std::vector<double>& out) const {
  const int S  = size();
  const int Np = parent_.n;
  const int nm = cfg_.n_moments;
  out.assign((size_t) 3 * S, 0.);

  std::vector<double> df(std::max(fermion_.n, boson_.n)), mom(nm);
  auto put_daughter = [&](int vec, int d, bool linear) {
    const GridView& g = (d == 0) ? fermion_ : boson_;
    for (int k = 0; k < g.n; ++k)
      df[k] = basis_[d].w[k] * (linear ? g.q[k] : 1.);
    Project(d, df.data(), mom.data());
    for (int j = 0; j < nm; ++j)
      out[(size_t) vec * S + Np + d * nm + j] = mom[j];
  };
  auto put_parent = [&](int vec, bool linear) {
    for (int b = 0; b < Np; ++b)
      out[(size_t) vec * S + b] = w_H_[b] * (linear ? parent_.q[b] : 1.);
  };

  // mu_H = mu_l + mu_phi leaves two independent chemical-potential directions.
  put_parent(0, false);
  put_daughter(0, 0, false);
  put_parent(1, false);
  put_daughter(1, 1, false);
  // A common temperature shift: df_i = w_i * q for every species.
  put_parent(2, true);
  put_daughter(2, 0, true);
  put_daughter(2, 1, true);
}

void ReducedCollisionOperator::EquilibriumNullVectorsExact(std::vector<double>& out) const {
  const int N  = full_size();
  const int Np = parent_.n;
  const int Nf = fermion_.n;
  out.assign((size_t) 3 * N, 0.);

  auto put = [&](int vec, bool linear, bool with_f, bool with_b) {
    for (int b = 0; b < Np; ++b)
      out[(size_t) vec * N + b] = w_H_[b] * (linear ? parent_.q[b] : 1.);
    if (with_f)
      for (int k = 0; k < Nf; ++k)
        out[(size_t) vec * N + Np + k] = basis_[0].w[k] * (linear ? fermion_.q[k] : 1.);
    if (with_b)
      for (int k = 0; k < boson_.n; ++k)
        out[(size_t) vec * N + Np + Nf + k] = basis_[1].w[k] * (linear ? boson_.q[k] : 1.);
  };
  put(0, false, true, false);
  put(1, false, false, true);
  put(2, true, true, true);
}

// --------------------------------------------------------- the tabulation --

void ReducedOperatorTable::Reset(int size, int l_max) {
  size_  = size;
  l_max_ = l_max;
  a_.clear();
  ln_a_.clear();
  M_.clear();
}

void ReducedOperatorTable::AddRow(double a, const std::vector<double>& M) {
  const size_t block = (size_t) (l_max_ + 1) * size_ * size_;
  if (size_ <= 0 || l_max_ < 0)
    throw std::runtime_error("ReducedOperatorTable: AddRow before Reset");
  if (M.size() != block)
    throw std::runtime_error("ReducedOperatorTable: row has the wrong block size");
  // Ascending is not a convenience, it is what makes the binary search in Apply valid;
  // an out-of-order row would silently return the wrong operator at every a above it.
  if (!a_.empty() && !(a > a_.back()))
    throw std::runtime_error("ReducedOperatorTable: rows must be added in ascending a");
  if (!(a > 0.))
    throw std::runtime_error("ReducedOperatorTable: row at non-positive a");
  a_.push_back(a);
  ln_a_.push_back(std::log(a));
  // Uniformity in ln a, tracked as rows arrive so the lookup can skip the search.
  if (a_.size() == 2) {
    const double d = ln_a_[1] - ln_a_[0];
    inv_dlna_      = (d != 0.) ? 1. / d : 0.;
    uniform_ln_a_  = (d > 0.);
  }
  else if (a_.size() > 2 && uniform_ln_a_) {
    const size_t n = ln_a_.size();
    const double d = ln_a_[n - 1] - ln_a_[n - 2];
    if (std::fabs(d * inv_dlna_ - 1.) > kUniformTol)
      uniform_ln_a_ = false;
  }
  M_.insert(M_.end(), M.begin(), M.end());
}

int ReducedOperatorTable::LowerRow(double a) const {
  int lo = 0, hi = static_cast<int>(a_.size()) - 1;
  if (hi < 0 || a < a_[0])
    return -1;
  if (uniform_ln_a_) {
    // One division. Clamped rather than trusted: log() and the stored ln_a_ need not
    // round identically at a row boundary, and an index off by one past the end would
    // read the next row's block.
    int i = static_cast<int>((std::log(a) - ln_a_[0]) * inv_dlna_);
    if (i < 0)
      i = 0;
    if (i > hi)
      i = hi;
    while (i < hi && a_[i + 1] <= a)
      ++i;
    while (i > 0 && a_[i] > a)
      --i;
    return i;
  }
  while (lo < hi) {
    const int mid = (lo + hi + 1) / 2;
    if (a_[mid] <= a)
      lo = mid;
    else
      hi = mid - 1;
  }
  return lo;
}

void ReducedOperatorTable::Apply(double a, int l, const double* z, double* dz) const {
  if (a_.empty() || l < 0 || l > l_max_ || a < a_.front() || a > a_.back())
    return;  // outside the window the collision is off; add nothing

  const int i        = LowerRow(a);
  const size_t block = (size_t) (l_max_ + 1) * size_ * size_;
  const size_t off_l = (size_t) l * size_ * size_;
  const double* A    = M_.data() + (size_t) i * block + off_l;

  // At or beyond the last row there is nothing to blend with; the exact-row case is
  // also taken here, which is what makes a tabulated a reproduce its row bit for bit
  // rather than through a weight that happens to be 1.
  if (i + 1 >= static_cast<int>(a_.size()) || a == a_[i]) {
    for (int r = 0; r < size_; ++r) {
      double s = 0.;
      for (int c = 0; c < size_; ++c)
        s += A[(size_t) r * size_ + c] * z[c];
      dz[r] += s;
    }
    return;
  }

  const double* B = M_.data() + (size_t) (i + 1) * block + off_l;
  const double t  = (std::log(a) - ln_a_[i]) / (ln_a_[i + 1] - ln_a_[i]);
  const double u  = 1. - t;
  for (int r = 0; r < size_; ++r) {
    double s = 0.;
    for (int c = 0; c < size_; ++c)
      s += (u * A[(size_t) r * size_ + c] + t * B[(size_t) r * size_ + c]) * z[c];
    dz[r] += s;
  }
}

void ReducedOperatorTable::Diagonal(double a, int l, double* out) const {
  for (int r = 0; r < size_; ++r)
    out[r] = 0.;
  if (a_.empty() || l < 0 || l > l_max_ || a < a_.front() || a > a_.back())
    return;

  const int i        = LowerRow(a);
  const size_t block = (size_t) (l_max_ + 1) * size_ * size_;
  const size_t off_l = (size_t) l * size_ * size_;
  const double* A    = M_.data() + (size_t) i * block + off_l;
  if (i + 1 >= static_cast<int>(a_.size()) || a == a_[i]) {
    for (int r = 0; r < size_; ++r)
      out[r] = A[(size_t) r * size_ + r];
    return;
  }
  const double* B = M_.data() + (size_t) (i + 1) * block + off_l;
  const double t  = (std::log(a) - ln_a_[i]) / (ln_a_[i + 1] - ln_a_[i]);
  for (int r = 0; r < size_; ++r)
    out[r] = (1. - t) * A[(size_t) r * size_ + r] + t * B[(size_t) r * size_ + r];
}
