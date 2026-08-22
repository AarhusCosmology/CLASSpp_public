#include "evolver_erk.h"

#include <atomic>
#include <cmath>

/**
 * Out-of-line state for the templated explicit-RK driver in
 * evolver_erk_impl.h: the process-wide controller configuration and the
 * cross-thread step counters. Everything numerical lives in the header so the
 * tableaux stay constexpr at the point of use.
 */

namespace {
ErkControllerConfig g_config;

std::atomic<bool> g_stats_on{false};
bool g_hist_on = false;
std::atomic<long long> g_accepted{0};
std::atomic<long long> g_rejected{0};
std::atomic<long long> g_derivs{0};
std::atomic<long long> g_dense{0};
std::atomic<long long> g_exact{0};

std::atomic<long long> g_err_acc[ErkHistograms::kErrBins];
std::atomic<long long> g_err_rej[ErkHistograms::kErrBins];
std::atomic<long long> g_x_acc[ErkHistograms::kXBins];
std::atomic<long long> g_x_rej[ErkHistograms::kXBins];

int bin_of(double value, double lo, double step, int nbins) {
  if (!(value > 0.0))
    return 0;
  const int b = static_cast<int>((log10(value) - lo) / step);
  return b < 0 ? 0 : (b >= nbins ? nbins - 1 : b);
}
}  // namespace

void evolver_erk_configure(const ErkControllerConfig& config) {
  g_config = config;
}
const ErkControllerConfig& evolver_erk_config() {
  return g_config;
}

void evolver_erk_stats_enable(bool on) {
  g_stats_on.store(on, std::memory_order_relaxed);
}
bool evolver_erk_stats_enabled() {
  return g_stats_on.load(std::memory_order_relaxed);
}
void evolver_erk_stats_reset() {
  g_accepted.store(0, std::memory_order_relaxed);
  g_rejected.store(0, std::memory_order_relaxed);
  g_derivs.store(0, std::memory_order_relaxed);
  g_dense.store(0, std::memory_order_relaxed);
  g_exact.store(0, std::memory_order_relaxed);
  for (int i = 0; i < ErkHistograms::kErrBins; i++) {
    g_err_acc[i].store(0, std::memory_order_relaxed);
    g_err_rej[i].store(0, std::memory_order_relaxed);
  }
  for (int i = 0; i < ErkHistograms::kXBins; i++) {
    g_x_acc[i].store(0, std::memory_order_relaxed);
    g_x_rej[i].store(0, std::memory_order_relaxed);
  }
}
ErkStats evolver_erk_stats_get() {
  ErkStats s;
  s.steps_accepted = g_accepted.load(std::memory_order_relaxed);
  s.steps_rejected = g_rejected.load(std::memory_order_relaxed);
  s.derivs_calls   = g_derivs.load(std::memory_order_relaxed);
  s.dense_points   = g_dense.load(std::memory_order_relaxed);
  s.exact_points   = g_exact.load(std::memory_order_relaxed);
  return s;
}

void evolver_erk_histograms_enable(bool on) {
  g_hist_on = on;
}

ErkHistograms evolver_erk_histograms_get() {
  ErkHistograms h;
  for (int i = 0; i < ErkHistograms::kErrBins; i++) {
    h.err_accepted[i] = g_err_acc[i].load(std::memory_order_relaxed);
    h.err_rejected[i] = g_err_rej[i].load(std::memory_order_relaxed);
  }
  for (int i = 0; i < ErkHistograms::kXBins; i++) {
    h.x_accepted[i] = g_x_acc[i].load(std::memory_order_relaxed);
    h.x_rejected[i] = g_x_rej[i].load(std::memory_order_relaxed);
  }
  return h;
}

namespace erk_detail {
void CountDerivs() {
  g_derivs.fetch_add(1, std::memory_order_relaxed);
}
void CountAccepted(double x, double err_ratio) {
  g_accepted.fetch_add(1, std::memory_order_relaxed);
  if (!g_hist_on)
    return;
  g_err_acc
      [bin_of(err_ratio, ErkHistograms::kErrLo, ErkHistograms::kErrStep, ErkHistograms::kErrBins)]
          .fetch_add(1, std::memory_order_relaxed);
  g_x_acc[bin_of(fabs(x), ErkHistograms::kXLo, ErkHistograms::kXStep, ErkHistograms::kXBins)]
      .fetch_add(1, std::memory_order_relaxed);
}
void CountRejected(double x, double err_ratio) {
  g_rejected.fetch_add(1, std::memory_order_relaxed);
  if (!g_hist_on)
    return;
  g_err_rej
      [bin_of(err_ratio, ErkHistograms::kErrLo, ErkHistograms::kErrStep, ErkHistograms::kErrBins)]
          .fetch_add(1, std::memory_order_relaxed);
  g_x_rej[bin_of(fabs(x), ErkHistograms::kXLo, ErkHistograms::kXStep, ErkHistograms::kXBins)]
      .fetch_add(1, std::memory_order_relaxed);
}
void CountDense() {
  g_dense.fetch_add(1, std::memory_order_relaxed);
}
void CountExact() {
  g_exact.fetch_add(1, std::memory_order_relaxed);
}
}  // namespace erk_detail
