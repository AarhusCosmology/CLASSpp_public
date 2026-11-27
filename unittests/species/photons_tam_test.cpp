/* The total angular momentum (TAM) photon hierarchy (PhotonsTamSpecies): its coupling
   coefficients against their closed forms, its tensor tight-coupling state against the
   optimal basis's, and the input key that chooses it.

   - {}_s kappa_l^m = sqrt((l^2 - m^2)(l^2 - s^2)/l^2) sqrt(q^2 - K l^2), q^2 = k^2 + (1+m) K
     (2005.12119, eq. 5). For m = 0 the temperature coefficient is k l s_l, which is how the
     scalar temperature hierarchy uses it.
   - The tight-coupling state of the tensor hierarchy is a physical state: written in either
     basis it must give the same polarisation source P^(2) and the same anisotropic stress in
     the gravitational-wave source. At leading order it is Theta_2 = -4/3 h'/kappa',
     E_2 = sqrt(2/3) h'/kappa' (C. Pitrou, class_public v3.4.0).
   - `hierarchy = optimal | tam` builds one photon species of that basis; with
     `<name>.type = photons`, `<name>.hierarchy` is the same key; anything else is a severe
     input error. */

#include "photons_tam.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

#include "background.h"
#include "constants.h"
#include "cosmology.h"
#include "input_module.h"
#include "photons_optimal.h"

namespace {

bool Close(double a, double b, double rtol = 1e-12) {
  return std::fabs(a - b) <= rtol * std::max({std::fabs(a), std::fabs(b), 1e-300});
}

constexpr double kK = 0.05; /* Mpc^-1 */

// ── Coefficients ──────────────────────────────────────────────────────────────

void test_m0_temperature_coefficient_is_k_l_s_l() {
  for (const double K : {0., 2.5e-4, -2.5e-4}) {
    PhotonsTamSpecies::Coefficients c;
    PhotonsTamSpecies::FillCoefficients(c, kK, K, 0, 20);
    for (int l = 1; l <= 20; ++l) {
      const double s_l = std::sqrt(std::max(1. - K * (l * l - 1.) / (kK * kK), 0.));
      assert(Close(c.kappa0[l], kK * l * s_l));
    }
  }
}

void test_flat_coefficients() {
  for (const int m : {0, 2}) {
    PhotonsTamSpecies::Coefficients c;
    PhotonsTamSpecies::FillCoefficients(c, kK, 0., m, 20);
    assert(Close(c.q, kK));
    for (int l = 3; l <= 20; ++l) {
      assert(Close(c.kappa0[l], kK * std::sqrt(l * l - m * m)));
      assert(Close(c.kappa2[l], kK * std::sqrt((l * l - m * m) * (l * l - 4.)) / l));
    }
  }
}

void test_curved_coefficients() {
  const double K = -3e-4;
  const int m    = 2;
  PhotonsTamSpecies::Coefficients c;
  PhotonsTamSpecies::FillCoefficients(c, kK, K, m, 12);
  const double q2 = kK * kK + 3. * K;
  assert(Close(c.q, std::sqrt(q2)));
  const int l = 7;
  assert(Close(c.kappa2[l], std::sqrt((49. - 4.) * (49. - 4.) / 49. * (q2 - 49. * K))));
  assert(Close(c.kappa0[l], std::sqrt((49. - 4.) * (q2 - 49. * K))));
}

void test_coefficients_vanish_at_and_below_m_and_s() {
  PhotonsTamSpecies::Coefficients c;
  PhotonsTamSpecies::FillCoefficients(c, kK, 1e-4, 2, 8);
  for (int l = 0; l <= 2; ++l) {
    assert(c.kappa0[l] == 0. && c.kappa2[l] == 0.);
  }
  assert(PhotonsTamSpecies::Kappa(2, 0, 2, kK * kK, 0.) == 0.);
  assert(PhotonsTamSpecies::Kappa(0, 0, 0, kK * kK, 0.) == 0.);
}

// ── Tensor tight coupling ─────────────────────────────────────────────────────

constexpr double kGwdot  = 3.7e-5;
constexpr double kDkappa = 2.4e3;
constexpr double kDrift  = 1.3e-6;

PhotonsTamSpecies::PerturbLayout TamTensorLayout(int l_max, int l_max_pol) {
  PhotonsTamSpecies::PerturbLayout lay;
  lay.l_max     = l_max;
  lay.l_max_pol = l_max_pol;
  lay.idx_t2    = 0;
  lay.idx_e2    = lay.idx_t2 + l_max - 1;
  lay.idx_b2    = lay.idx_e2 + l_max_pol - 1;
  return lay;
}

PhotonsOptimalSpecies::PerturbLayout OptimalTensorLayout(int l_max, int l_max_pol) {
  PhotonsOptimalSpecies::PerturbLayout lay;
  lay.l_max     = l_max;
  lay.l_max_pol = l_max_pol;
  lay.idx_delta = 0;
  lay.idx_theta = 1;
  lay.idx_shear = 2;
  lay.idx_l3    = 3;
  lay.idx_pol0  = l_max + 1;
  lay.idx_pol1  = lay.idx_pol0 + 1;
  lay.idx_pol2  = lay.idx_pol0 + 2;
  lay.idx_pol3  = lay.idx_pol0 + 3;
  return lay;
}

void test_tensor_state_is_the_same_in_both_bases() {
  background pba{};
  const PhotonsTamSpecies tam(pba);
  const PhotonsOptimalSpecies opt(pba);
  for (const double drift : {0., kDrift}) {
    const PhotonsSpecies::TensorTcaClosure tca =
        PhotonsSpecies::TensorTightCoupling(kGwdot, kDkappa, drift);

    const auto tlay = TamTensorLayout(7, 7);
    std::vector<double> yt(tlay.idx_b2 + tlay.l_max_pol - 1, 0.);
    tam.StartTensorHierarchy(tlay, yt.data(), tca);

    const auto olay = OptimalTensorLayout(5, 5);
    std::vector<double> yo(olay.idx_pol0 + olay.l_max_pol + 1, 0.);
    opt.StartTensorHierarchy(olay, yo.data(), tca);

    /* the polarisation source and the stress: physical, so equal in both bases */
    assert(Close(tam.TensorPolarisationSource(tlay, yt.data()), tca.P2));
    assert(Close(opt.TensorPolarisationSource(olay, yo.data()), tca.P2));
    const double a2    = 3.2e-7;
    const double rho_g = 4.1e3;
    assert(Close(tam.TensorGwSource(tlay, yt.data(), a2, rho_g),
                 opt.TensorGwSource(olay, yo.data(), a2, rho_g)));

    /* B_2 and every higher multipole stay at zero */
    for (int i = 0; i < static_cast<int>(yt.size()); ++i) {
      if (i != tlay.idx_t2 && i != tlay.idx_e2) {
        assert(yt[i] == 0.);
      }
    }
  }
}

void test_tensor_state_at_leading_order() {
  background pba{};
  const PhotonsTamSpecies tam(pba);
  const auto lay = TamTensorLayout(7, 7);
  std::vector<double> y(lay.idx_b2 + lay.l_max_pol - 1, 0.);
  tam.StartTensorHierarchy(lay, y.data(), PhotonsSpecies::TensorTightCoupling(kGwdot, kDkappa));
  const double u = kGwdot / kDkappa;
  assert(Close(y[lay.idx_t2], -4. / 3. * u));
  assert(Close(y[lay.idx_e2], std::sqrt(2. / 3.) * u));
}

// ── The input key ─────────────────────────────────────────────────────────────

FileContent Quiet() {
  FileContent fc;
  fc.set("output", "");
  fc.set("input_verbose", "0");
  fc.set("background_verbose", "0");
  fc.set("thermodynamics_verbose", "0");
  fc.set("write warnings", "no");
  fc.set("write parameters", "no");
  return fc;
}

/* "optimal", "tam", or the message of the std::invalid_argument the build threw */
std::string BuiltBasis(FileContent fc) {
  try {
    Cosmology cosmology{fc};
    const BaseSpecies& g = cosmology.GetInputModule()->all_species_.photons();
    if (dynamic_cast<const PhotonsTamSpecies*>(&g) != nullptr) {
      return "tam";
    }
    if (dynamic_cast<const PhotonsOptimalSpecies*>(&g) != nullptr) {
      return "optimal";
    }
    return "neither";
  }
  catch (const std::invalid_argument& e) {
    return e.what();
  }
}

void test_hierarchy_key() {
  assert(BuiltBasis(Quiet()) == "optimal");

  FileContent opt = Quiet();
  opt.set("hierarchy", "optimal");
  assert(BuiltBasis(opt) == "optimal");

  for (const char* spelling : {"tam", "TAM"}) {
    FileContent fc = Quiet();
    fc.set("hierarchy", spelling);
    assert(BuiltBasis(fc) == "tam");
  }

  FileContent dot = Quiet();
  dot.set("g.type", "photons");
  dot.set("g.hierarchy", "tam");
  assert(BuiltBasis(dot) == "tam");

  FileContent bad = Quiet();
  bad.set("hierarchy", "tma");
  const std::string message = BuiltBasis(bad);
  std::printf("hierarchy = tma: %s\n", message.c_str());
  assert(message.find("hierarchy") != std::string::npos &&
         message.find("tma") != std::string::npos);
}

}  // namespace

int main() {
  test_m0_temperature_coefficient_is_k_l_s_l();
  test_flat_coefficients();
  test_curved_coefficients();
  test_coefficients_vanish_at_and_below_m_and_s();
  test_tensor_state_is_the_same_in_both_bases();
  test_tensor_state_at_leading_order();
  test_hierarchy_key();
  std::printf("photons TAM tests passed\n");
  return 0;
}
