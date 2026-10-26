#ifndef BESSEL_TABLE_CACHE_H
#define BESSEL_TABLE_CACHE_H

#include <memory>
#include <mutex>
#include <vector>

#include "hyperspherical.h"

/**
 * Keeps the transfer module's table of flat-space spherical Bessel functions from one
 * cosmology to the next.
 *
 * The table depends on the multipoles, on the x range and on three precision parameters, but
 * on the cosmology only through how far in x it must reach. Built afresh for every cosmology,
 * it cost 24 ms of 0.85 s on one thread for LCDM to l = 3000, in a serial step ahead of the
 * threaded transfer loop. A sampler calls compute() thousands of times on one classy object,
 * which owns one cache (CAMB keeps its Bessel tables for the same reason).
 *
 * The cached table is reused only when it was built for the same multipoles and parameters
 * and reaches at least as far in x. Tables are built on a fixed spacing
 * (HyperInterpStruct::OnFixedSpacing), so a reused table holds exactly the nodes and values a
 * fresh one would. A computation therefore gives the same result whatever is cached.
 * docs/superpowers/specs/2026-10-26-bessel-table-cache-design.md
 */
class BesselTableCache {
 public:
  /** The table for these multipoles on [xmin, xmax]: the cached one if it serves, else a new
      one, which then replaces it. */
  std::shared_ptr<const HyperInterpStruct> Get(const std::vector<int>& l,
                                               double xmin,
                                               double xmax,
                                               double sampling,
                                               int l_WKB,
                                               double phiminabs) {
    std::scoped_lock guard(mutex_);
    const bool reusable = table_ && l == l_ && xmin == xmin_ && sampling == sampling_ &&
                          l_WKB == l_WKB_ && phiminabs == phiminabs_ && table_->x.back() >= xmax;
    if (!reusable) {
      /* With some headroom, so that the next cosmology, which may need a little more x, can
         reuse it too. The headroom adds nodes, it does not move them. */
      table_ = std::make_shared<const HyperInterpStruct>(
          HyperInterpStruct::OnFixedSpacing(0,
                                            1.,
                                            static_cast<int>(l.size()),
                                            l.data(),
                                            xmin,
                                            kHeadroom * xmax,
                                            sampling,
                                            l_WKB,
                                            phiminabs));
      l_         = l;
      xmin_      = xmin;
      sampling_  = sampling;
      l_WKB_     = l_WKB;
      phiminabs_ = phiminabs;
    }
    return table_;
  }

 private:
  static constexpr double kHeadroom = 1.05;

  std::mutex mutex_;
  std::shared_ptr<const HyperInterpStruct> table_;
  std::vector<int> l_;
  double xmin_      = 0.;
  double sampling_  = 0.;
  int l_WKB_        = 0;
  double phiminabs_ = 0.;
};

#endif  // BESSEL_TABLE_CACHE_H
