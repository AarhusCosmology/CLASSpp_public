#include "energy_deposition.h"

#include <algorithm>
#include <array>
#include <cmath>

#include "arrays.h"

namespace {

constexpr int kLines    = 18;
constexpr int kChannels = 5;

/* Energy fractions from DM annihilation or decay, computed by Galli et al. 2013
   (arXiv:1306.0563) and distributed with class_public as
   external/heating/Galli_et_al_2013.dat. The table is small enough to live here
   rather than as a data file needing an install rule and a path parameter.

   Columns, in file order: heat, lya, ion_H, ion_He, lowE (E < 10.2 eV).

   The last two rows carry the table past x = 1, which CLASS's
   x = x_H + f_He x_He reaches while helium is still ionized. */

// clang-format off
constexpr double kX[kLines] = {
    0.000000, 0.000100, 0.000300, 0.000500, 0.001000, 0.003000,
    0.005000, 0.010000, 0.030000, 0.050000, 0.100000, 0.300000,
    0.500000, 0.800000, 0.900000, 0.990000, 1.100000, 1.200000,
};

constexpr double kChi[kLines*kChannels] = {
    0.151800, 0.323526, 0.350798, 0.024367, 0.142202,
    0.151880, 0.323526, 0.350798, 0.024367, 0.142202,
    0.174825, 0.308840, 0.349058, 0.023397, 0.136233,
    0.188520, 0.302591, 0.345508, 0.023737, 0.133317,
    0.210027, 0.291280, 0.341822, 0.023134, 0.128324,
    0.258912, 0.269481, 0.327298, 0.023415, 0.118130,
    0.289871, 0.256105, 0.316798, 0.023029, 0.111856,
    0.338316, 0.238304, 0.301893, 0.021302, 0.103617,
    0.458621, 0.192119, 0.255925, 0.018550, 0.083760,
    0.531628, 0.165741, 0.228453, 0.016601, 0.072038,
    0.654816, 0.121705, 0.175739, 0.012852, 0.053388,
    0.849031, 0.052273, 0.083885, 0.007317, 0.023016,
    0.923644, 0.026168, 0.043901, 0.004914, 0.011310,
    0.975679, 0.007178, 0.013518, 0.003640, 0.003117,
    0.987026, 0.003234, 0.006406, 0.003223, 0.001417,
    0.995299, 0.000290, 0.001700, 0.002924, 0.000118,
    1.000000, 0.000000, 0.000000, 0.000000, 0.000000,
    1.000000, 0.000000, 0.000000, 0.000000, 0.000000,
};
// clang-format on

/* Second derivatives for the cubic spline through kChi. Built once on first use;
   the table is compile-time constant, so this needs no per-run state. */
struct GalliSpline {
  std::array<double, kLines * kChannels> ddchi;

  GalliSpline() {
    array_spline_table_lines(kX, kLines, kChi, kChannels, ddchi.data(), _SPLINE_NATURAL_);
  }
};

double clamp01(double v) {
  return std::min(std::max(v, 0.), 1.);
}

}  // namespace

EnergyDeposition energy_deposition_fractions(energy_deposition_function which, double x) {
  EnergyDeposition dep = {0., 0., 0., 0., 0.};

  if (which == deposition_legacy) {
    /* Two analytic fits by Vivian Poulin of columns 1 and 2 in Table V of
       Slatyer et al. 2013, as carried inline in the thermodynamics module before
       the Galli table was available here. They cover only heating and hydrogen
       ionization; RECFAST's energy-injection term reused the ionization
       coefficient for the Lyman-alpha channel, so that is what `lya` reports and
       old results reproduce exactly. */
    if (x < 1.) {
      dep.heat  = std::min(0.996857 * (1. - std::pow(1. - std::pow(x, 0.300134), 1.51035)), 1.0);
      dep.ion_H = 0.369202 * std::pow(1. - std::pow(x, 0.463929), 1.70237);
    }
    else {
      dep.heat  = 1.;
      dep.ion_H = 0.;
    }
    dep.lya = dep.ion_H;
    return dep;
  }

  static const GalliSpline spline;

  double chi[kChannels];
  int last_index;
  array_interpolate_spline(kX,
                           kLines,
                           kChi,
                           spline.ddchi.data(),
                           kChannels,
                           std::min(std::max(x, kX[0]), kX[kLines - 1]),
                           &last_index,
                           chi,
                           kChannels);

  /* A cubic spline through the shoulder near x = 1 rings slightly; clamping keeps
     every channel a physical fraction. */
  dep.heat   = clamp01(chi[0]);
  dep.lya    = clamp01(chi[1]);
  dep.ion_H  = clamp01(chi[2]);
  dep.ion_He = clamp01(chi[3]);
  dep.lowE   = clamp01(chi[4]);

  return dep;
}
