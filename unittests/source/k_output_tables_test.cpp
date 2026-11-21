// The perturbations at k_output_values are recorded while the perturbations are integrated.
// perturb_init queues every initial condition and wavenumber at once, so the ad and cdi
// integrations of one k_output_values mode run on two threads at the same time. Until the
// tables were kept per initial condition, both appended to one buffer: with ic = ad,cdi and
// eight threads, ./class died with SIGTRAP in eight runs out of eight (#504, 2026-11-21).
//
// The table of a mode is the rows of each initial condition in turn, as class_public writes
// them (its loop over initial conditions is outside the parallel region). So the table must
// not depend on the number of threads, bit for bit, and its first rows are the ad-only table.
#include <cassert>
#include <cstddef>
#include <string>
#include <vector>

#include "cosmology.h"
#include "perturbations_module.h"

namespace {

std::vector<OutputTable> Run(const std::string& ic, int threads) {
  FileContent fc;
  fc.set("output", "mPk");
  fc.set("ic", ic);
  if (ic != "ad") {
    fc.set("f_cdi", "0.5");
  }
  fc.set("k_output_values", "0.01");
  fc.set("P_k_max_1/Mpc", "1e-5");  // a few modes beside the k_output_values one
  fc.set("threads", std::to_string(threads));
  fc.set("input_verbose", "0");
  fc.set("perturbations_verbose", "0");
  fc.set("write warnings", "no");
  fc.set("write parameters", "no");
  Cosmology cosmology{fc};
  const auto& pt = cosmology.GetPerturbationsModule();
  return pt->GetPerturbationTables(pt->index_md_scalars_);
}

bool SameRows(const OutputTable& a, const OutputTable& b, std::size_t rows) {
  if (a.names() != b.names()) {
    return false;
  }
  for (std::size_t c = 0; c < a.names().size(); ++c) {
    for (std::size_t r = 0; r < rows; ++r) {
      if (a.column(c)[r] != b.column(c)[r]) {
        return false;
      }
    }
  }
  return true;
}

}  // namespace

int main() {
  const std::vector<OutputTable> serial = Run("ad,cdi", 1);
  assert(serial.size() == 1);
  const OutputTable& table = serial[0];

  // ad, then cdi: the time goes back once, where the cdi rows start.
  const OutputColumn tau = table.column("tau [Mpc]");
  std::size_t restarts   = 0;
  std::size_t ad_rows    = 0;
  for (std::size_t r = 1; r < tau.size(); ++r) {
    if (tau[r] < tau[r - 1]) {
      ++restarts;
      ad_rows = r;
    }
  }
  assert(restarts == 1);

  const std::vector<OutputTable> ad = Run("ad", 1);
  assert(ad[0].rows() == ad_rows);
  assert(SameRows(ad[0], table, ad_rows));

  for (int attempt = 0; attempt < 4; ++attempt) {
    const std::vector<OutputTable> threaded = Run("ad,cdi", 8);
    assert(threaded[0].rows() == table.rows());
    assert(SameRows(threaded[0], table, table.rows()));
  }
  return 0;
}
