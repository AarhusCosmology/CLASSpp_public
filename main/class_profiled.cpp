/** @file class_profiled.cpp
 *
 * Per-module timing harness for CLASS++.
 *
 * Drives the lazy Cosmology pipeline one module at a time, timing each stage
 * with a steady clock. Repeats the whole run N times (fresh Cosmology each
 * loop) so we can report min / median / mean / stddev per module and judge
 * which modules actually dominate the wall time -- and whether a proposed
 * optimization has any chance of moving the needle.
 *
 * Usage:  ./class_profiled <file.ini> [num_loops]   (default num_loops = 5)
 *
 * Timings are pure wall time of each GetXxxModule() call. Because the modules
 * are memoized inside Cosmology and built in dependency order here, each call
 * does only its own work.
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "cosmology.h"
#include "input_module.h"

namespace {

using Clock = std::chrono::steady_clock;

double ms_since(const Clock::time_point& t0) {
  return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
}

struct Stats {
  double min, median, mean, stddev;
};

Stats summarize(std::vector<double> v) {
  Stats s{0, 0, 0, 0};
  if (v.empty())
    return s;
  std::sort(v.begin(), v.end());
  s.min          = v.front();
  const size_t n = v.size();
  s.median       = (n % 2 == 1) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
  double sum     = 0.0;
  for (double x : v)
    sum += x;
  s.mean     = sum / n;
  double acc = 0.0;
  for (double x : v)
    acc += (x - s.mean) * (x - s.mean);
  s.stddev = (n > 1) ? std::sqrt(acc / (n - 1)) : 0.0;
  return s;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    printf("Usage: %s <file.ini> [num_loops]\n", argv[0]);
    return _FAILURE_;
  }
  int num_loops = (argc >= 3) ? std::atoi(argv[2]) : 5;
  if (num_loops < 1)
    num_loops = 1;

  // Parse arguments once into a master FileContent; copy it per loop so each
  // run starts from an identical, all-unread parameter set.
  FileContent fc_master;
  {
    char* prof_argv[2] = {argv[0], argv[1]};
    InputModule::file_content_from_arguments(2, prof_argv, fc_master);
  }

  // Module stage labels, in dependency order.
  const std::vector<std::string> labels = {"Input(+shoot)",
                                           "Background",
                                           "Thermodynamics",
                                           "Perturbations",
                                           "Primordial",
                                           "Nonlinear",
                                           "Transfer",
                                           "Spectra",
                                           "Lensing"};
  const size_t n_stages                 = labels.size();
  std::vector<std::vector<double>> samples(n_stages);
  std::vector<double> totals;

  printf("Profiling '%s' over %d loop(s)...\n", argv[1], num_loops);

  for (int loop = 0; loop < num_loops; ++loop) {
    try {
      FileContent fc = fc_master;  // fresh, all-unread copy
      auto t_total   = Clock::now();
      std::vector<double> stage_ms(n_stages, 0.0);

      auto t = Clock::now();
      Cosmology cosmology{fc};
      cosmology.GetInputModule();  // triggers shooting
      stage_ms[0] = ms_since(t);

      t = Clock::now();
      cosmology.GetBackgroundModule();
      stage_ms[1] = ms_since(t);
      t           = Clock::now();
      cosmology.GetThermodynamicsModule();
      stage_ms[2] = ms_since(t);
      t           = Clock::now();
      cosmology.GetPerturbationsModule();
      stage_ms[3] = ms_since(t);
      t           = Clock::now();
      cosmology.GetPrimordialModule();
      stage_ms[4] = ms_since(t);
      t           = Clock::now();
      cosmology.GetNonlinearModule();
      stage_ms[5] = ms_since(t);
      t           = Clock::now();
      cosmology.GetTransferModule();
      stage_ms[6] = ms_since(t);
      t           = Clock::now();
      cosmology.GetSpectraModule();
      stage_ms[7] = ms_since(t);
      t           = Clock::now();
      cosmology.GetLensingModule();
      stage_ms[8] = ms_since(t);

      double total = ms_since(t_total);
      for (size_t i = 0; i < n_stages; ++i)
        samples[i].push_back(stage_ms[i]);
      totals.push_back(total);

      printf("  loop %2d/%d: total %8.2f ms\n", loop + 1, num_loops, total);
    }
    catch (const std::exception& e) {
      printf("\n\nError in loop %d:\n%s\n", loop + 1, e.what());
      return _FAILURE_;
    }
  }

  // Report.
  printf("\n");
  printf("%-16s %9s %9s %9s %9s %7s\n", "Module", "min", "median", "mean", "stddev", "%total");
  printf("%-16s %9s %9s %9s %9s %7s\n",
         "----------------",
         "--------",
         "--------",
         "--------",
         "--------",
         "------");
  Stats tot = summarize(totals);
  for (size_t i = 0; i < n_stages; ++i) {
    Stats s    = summarize(samples[i]);
    double pct = (tot.median > 0) ? 100.0 * s.median / tot.median : 0.0;
    printf("%-16s %9.2f %9.2f %9.2f %9.2f %6.1f%%\n",
           labels[i].c_str(),
           s.min,
           s.median,
           s.mean,
           s.stddev,
           pct);
  }
  printf("%-16s %9s %9s %9s %9s %7s\n",
         "----------------",
         "--------",
         "--------",
         "--------",
         "--------",
         "------");
  printf("%-16s %9.2f %9.2f %9.2f %9.2f %6.1f%%\n",
         "TOTAL",
         tot.min,
         tot.median,
         tot.mean,
         tot.stddev,
         100.0);

  return _SUCCESS_;
}
