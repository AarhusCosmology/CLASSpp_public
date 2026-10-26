/** @file wigner_d.cpp Wigner d-functions d^l_{mm'}(mu) for the lensing module, by recurrence in l.
 *
 * docs/superpowers/specs/2026-10-26-lensing-wigner-recursions-design.md
 */

#include "wigner_d.h"

#include <cmath>
#include <vector>

using wigner::kChunk;
using wigner::Nodes;

namespace {

/** How mu enters the recurrence of a d^l_{mm'}: alone, or shifted by +b_l or -b_l. Three
    spellings of the expression rather than a sign factor, so that each recurrence does
    exactly the arithmetic it did when it had a loop of its own. */
enum class WignerShift { kNone, kPlus, kMinus };

/** The recurrence of every d^l_{mm'} below, for kChunk nodes at once:
      D_{l+1} = a_l (mu [+ or - b_l]) D_l - c_l D_{l-1},   d^{l+1}_{mm'} = e_l D_{l+1},
    with D = sqrt((2l+1)/2) d, the quantity Kostelec & Rockmore (2003) recur on for stability.
    Each node's chain is serial in l and the nodes are independent, so l runs outside and the
    nodes inside: the nodes fill the vector lanes and the pipeline, and the table
    d[l * kChunk + node] is written contiguously. Run one node at a time with l inside,
    the chains were latency-bound at ~2 ns a value, two thirds of the lensing module.
    On entry dl and dlm1 hold D_l and D_{l-1} of each node at l = l_first. */
template <WignerShift Shift>
void WignerRecurrence(const std::array<double, kChunk>& mu,
                      int l_first,
                      int lmax,
                      const std::vector<double>& a,
                      const std::vector<double>& b,
                      const std::vector<double>& c,
                      const std::vector<double>& e,
                      std::array<double, kChunk>& dl,
                      std::array<double, kChunk>& dlm1,
                      double* d) {
  for (int l = l_first; l < lmax; l++) {
    const double al = a[l];
    const double bl = b[l];
    const double cl = c[l];
    const double el = e[l];
    double* row     = d + static_cast<size_t>(l + 1) * kChunk;
    for (int j = 0; j < kChunk; j++) {
      double dlp1 = 0.;
      if constexpr (Shift == WignerShift::kNone) {
        dlp1 = al * mu[j] * dl[j] - cl * dlm1[j];
      }
      else if constexpr (Shift == WignerShift::kPlus) {
        dlp1 = al * (mu[j] + bl) * dl[j] - cl * dlm1[j];
      }
      else {
        dlp1 = al * (mu[j] - bl) * dl[j] - cl * dlm1[j];
      }
      row[j]  = dlp1 * el;
      dlm1[j] = dl[j];
      dl[j]   = dlp1;
    }
  }
}

}  // namespace

namespace wigner {

/**
 * This routine computes the d00 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d00    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

void D00(const Nodes& mu, int lmax, double* d00) {
  std::vector<double> fac1(lmax);
  std::vector<double> fac2(lmax);
  std::vector<double> fac3(lmax);

  for (int l = 1; l < lmax; l++) {
    double ll = static_cast<double>(l);
    fac1[l]   = sqrt((2 * ll + 3) / (2 * ll + 1)) * (2 * ll + 1) / (ll + 1);
    fac2[l]   = sqrt((2 * ll + 3) / (2 * ll - 1)) * ll / (ll + 1);
    fac3[l]   = sqrt(2. / (2 * ll + 3));
  }

  Nodes dl{};
  Nodes dlm1{};
  for (int index_mu = 0; index_mu < kChunk; index_mu++) {
    dlm1[index_mu]             = 1.0 / sqrt(2.); /* l=0 */
    d00[0 * kChunk + index_mu] = dlm1[index_mu] * sqrt(2.);
    dl[index_mu]               = mu[index_mu] * sqrt(3. / 2.); /*l=1*/
    d00[1 * kChunk + index_mu] = dl[index_mu] * sqrt(2. / 3.);
  }
  /* sqrt((2l+1)/2)*d00 recurrence, supposed to be more stable */
  WignerRecurrence<WignerShift::kNone>(mu, 1, lmax, fac1, fac1, fac2, fac3, dl, dlm1, d00);
}

/**
 * This routine computes the d11 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d11    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

void D11(const Nodes& mu, int lmax, double* d11) {
  std::vector<double> fac1(lmax);
  std::vector<double> fac2(lmax);
  std::vector<double> fac3(lmax);
  std::vector<double> fac4(lmax);
  for (int l = 2; l < lmax; l++) {
    double ll = static_cast<double>(l);
    fac1[l]   = sqrt((2 * ll + 3) / (2 * ll + 1)) * (ll + 1) * (2 * ll + 1) / (ll * (ll + 2));
    fac2[l]   = 1.0 / (ll * (ll + 1.));
    fac3[l] = sqrt((2 * ll + 3) / (2 * ll - 1)) * (ll - 1) * (ll + 1) / (ll * (ll + 2)) * (ll + 1) /
              ll;
    fac4[l] = sqrt(2. / (2 * ll + 3));
  }

  Nodes dl{};
  Nodes dlm1{};
  for (int index_mu = 0; index_mu < kChunk; index_mu++) {
    d11[0 * kChunk + index_mu] = 0;
    dlm1[index_mu]             = (1.0 + mu[index_mu]) / 2. * sqrt(3. / 2.); /*l=1*/
    d11[1 * kChunk + index_mu] = dlm1[index_mu] * sqrt(2. / 3.);
    dl[index_mu] = (1.0 + mu[index_mu]) / 2. * (2.0 * mu[index_mu] - 1.0) * sqrt(5. / 2.); /*l=2*/
    d11[2 * kChunk + index_mu] = dl[index_mu] * sqrt(2. / 5.);
  }
  /* sqrt((2l+1)/2)*d11 recurrence, supposed to be more stable */
  WignerRecurrence<WignerShift::kMinus>(mu, 2, lmax, fac1, fac2, fac3, fac4, dl, dlm1, d11);
}

/**
 * This routine computes the d1m1 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d1m1    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

void D1m1(const Nodes& mu, int lmax, double* d1m1) {
  std::vector<double> fac1(lmax);
  std::vector<double> fac2(lmax);
  std::vector<double> fac3(lmax);
  std::vector<double> fac4(lmax);
  for (int l = 2; l < lmax; l++) {
    double ll = static_cast<double>(l);
    fac1[l]   = sqrt((2 * ll + 3) / (2 * ll + 1)) * (ll + 1) * (2 * ll + 1) / (ll * (ll + 2));
    fac2[l]   = 1.0 / (ll * (ll + 1.));
    fac3[l] = sqrt((2 * ll + 3) / (2 * ll - 1)) * (ll - 1) * (ll + 1) / (ll * (ll + 2)) * (ll + 1) /
              ll;
    fac4[l] = sqrt(2. / (2 * ll + 3));
  }

  Nodes dl{};
  Nodes dlm1{};
  for (int index_mu = 0; index_mu < kChunk; index_mu++) {
    d1m1[0 * kChunk + index_mu] = 0;
    dlm1[index_mu]              = (1.0 - mu[index_mu]) / 2. * sqrt(3. / 2.); /*l=1*/
    d1m1[1 * kChunk + index_mu] = dlm1[index_mu] * sqrt(2. / 3.);
    dl[index_mu] = (1.0 - mu[index_mu]) / 2. * (2.0 * mu[index_mu] + 1.0) * sqrt(5. / 2.); /*l=2*/
    d1m1[2 * kChunk + index_mu] = dl[index_mu] * sqrt(2. / 5.);
  }
  /* sqrt((2l+1)/2)*d1m1 recurrence, supposed to be more stable */
  WignerRecurrence<WignerShift::kPlus>(mu, 2, lmax, fac1, fac2, fac3, fac4, dl, dlm1, d1m1);
}

/**
 * This routine computes the d2m2 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d2m2   Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

void D2m2(const Nodes& mu, int lmax, double* d2m2) {
  std::vector<double> fac1(lmax);
  std::vector<double> fac2(lmax);
  std::vector<double> fac3(lmax);
  std::vector<double> fac4(lmax);
  for (int l = 2; l < lmax; l++) {
    double ll = static_cast<double>(l);
    fac1[l]   = sqrt((2 * ll + 3) / (2 * ll + 1)) * (ll + 1) * (2 * ll + 1) / ((ll - 1) * (ll + 3));
    fac2[l]   = 4.0 / (ll * (ll + 1));
    fac3[l]   = sqrt((2 * ll + 3) / (2 * ll - 1)) * (ll - 2) * (ll + 2) / ((ll - 1) * (ll + 3)) *
                (ll + 1) / ll;
    fac4[l]   = sqrt(2. / (2 * ll + 3));
  }

  Nodes dl{};
  Nodes dlm1{};
  for (int index_mu = 0; index_mu < kChunk; index_mu++) {
    d2m2[0 * kChunk + index_mu] = 0;
    dlm1[index_mu]              = 0.; /*l=1*/
    d2m2[1 * kChunk + index_mu] = 0;
    dl[index_mu] = (1.0 - mu[index_mu]) * (1.0 - mu[index_mu]) / 4. * sqrt(5. / 2.); /*l=2*/
    d2m2[2 * kChunk + index_mu] = dl[index_mu] * sqrt(2. / 5.);
  }
  /* sqrt((2l+1)/2)*d2m2 recurrence, supposed to be more stable */
  WignerRecurrence<WignerShift::kPlus>(mu, 2, lmax, fac1, fac2, fac3, fac4, dl, dlm1, d2m2);
}

/**
 * This routine computes the d22 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d22    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

void D22(const Nodes& mu, int lmax, double* d22) {
  std::vector<double> fac1(lmax);
  std::vector<double> fac2(lmax);
  std::vector<double> fac3(lmax);
  std::vector<double> fac4(lmax);
  for (int l = 2; l < lmax; l++) {
    double ll = static_cast<double>(l);
    fac1[l]   = sqrt((2 * ll + 3) / (2 * ll + 1)) * (ll + 1) * (2 * ll + 1) / ((ll - 1) * (ll + 3));
    fac2[l]   = 4.0 / (ll * (ll + 1));
    fac3[l]   = sqrt((2 * ll + 3) / (2 * ll - 1)) * (ll - 2) * (ll + 2) / ((ll - 1) * (ll + 3)) *
                (ll + 1) / ll;
    fac4[l]   = sqrt(2. / (2 * ll + 3));
  }

  Nodes dl{};
  Nodes dlm1{};
  for (int index_mu = 0; index_mu < kChunk; index_mu++) {
    d22[0 * kChunk + index_mu] = 0;
    dlm1[index_mu]             = 0.; /*l=1*/
    d22[1 * kChunk + index_mu] = 0;
    dl[index_mu] = (1.0 + mu[index_mu]) * (1.0 + mu[index_mu]) / 4. * sqrt(5. / 2.); /*l=2*/
    d22[2 * kChunk + index_mu] = dl[index_mu] * sqrt(2. / 5.);
  }
  /* sqrt((2l+1)/2)*d22 recurrence, supposed to be more stable */
  WignerRecurrence<WignerShift::kMinus>(mu, 2, lmax, fac1, fac2, fac3, fac4, dl, dlm1, d22);
}

/**
 * This routine computes the d20 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d20    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

void D20(const Nodes& mu, int lmax, double* d20) {
  std::vector<double> fac1(lmax);
  std::vector<double> fac3(lmax);
  std::vector<double> fac4(lmax);
  for (int l = 2; l < lmax; l++) {
    double ll = static_cast<double>(l);
    fac1[l]   = sqrt((2 * ll + 3) * (2 * ll + 1) / ((ll - 1) * (ll + 3)));
    fac3[l]   = sqrt((2 * ll + 3) * (ll - 2) * (ll + 2) / ((2 * ll - 1) * (ll - 1) * (ll + 3)));
    fac4[l]   = sqrt(2. / (2 * ll + 3));
  }

  Nodes dl{};
  Nodes dlm1{};
  for (int index_mu = 0; index_mu < kChunk; index_mu++) {
    d20[0 * kChunk + index_mu] = 0;
    dlm1[index_mu]             = 0.; /*l=1*/
    d20[1 * kChunk + index_mu] = 0;
    dl[index_mu]               = sqrt(15.) / 4. * (1 - mu[index_mu] * mu[index_mu]); /*l=2*/
    d20[2 * kChunk + index_mu] = dl[index_mu] * sqrt(2. / 5.);
  }
  /* sqrt((2l+1)/2)*d22 recurrence, supposed to be more stable */
  WignerRecurrence<WignerShift::kNone>(mu, 2, lmax, fac1, fac1, fac3, fac4, dl, dlm1, d20);
}

/**
 * This routine computes the d31 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d31    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

void D31(const Nodes& mu, int lmax, double* d31) {
  std::vector<double> fac1(lmax);
  std::vector<double> fac2(lmax);
  std::vector<double> fac3(lmax);
  std::vector<double> fac4(lmax);
  for (int l = 3; l < lmax; l++) {
    double ll = static_cast<double>(l);
    fac1[l] = sqrt((2 * ll + 3) * (2 * ll + 1) / ((ll - 2) * (ll + 4) * ll * (ll + 2))) * (ll + 1);
    fac2[l] = 3.0 / (ll * (ll + 1));
    fac3[l] = sqrt((2 * ll + 3) / (2 * ll - 1) * (ll - 3) * (ll + 3) * (ll - 1) * (ll + 1) /
                   ((ll - 2) * (ll + 4) * ll * (ll + 2))) *
              (ll + 1) / ll;
    fac4[l] = sqrt(2. / (2 * ll + 3));
  }

  Nodes dl{};
  Nodes dlm1{};
  for (int index_mu = 0; index_mu < kChunk; index_mu++) {
    d31[0 * kChunk + index_mu] = 0;
    d31[1 * kChunk + index_mu] = 0;
    dlm1[index_mu]             = 0.; /*l=2*/
    d31[2 * kChunk + index_mu] = 0;
    dl[index_mu] = sqrt(105. / 2.) * (1 + mu[index_mu]) * (1 + mu[index_mu]) * (1 - mu[index_mu]) /
                   8.; /*l=3*/
    d31[3 * kChunk + index_mu] = dl[index_mu] * sqrt(2. / 7.);
  }
  /* sqrt((2l+1)/2)*d22 recurrence, supposed to be more stable */
  WignerRecurrence<WignerShift::kMinus>(mu, 3, lmax, fac1, fac2, fac3, fac4, dl, dlm1, d31);
}

/**
 * This routine computes the d3m1 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d3m1   Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

void D3m1(const Nodes& mu, int lmax, double* d3m1) {
  std::vector<double> fac1(lmax);
  std::vector<double> fac2(lmax);
  std::vector<double> fac3(lmax);
  std::vector<double> fac4(lmax);
  for (int l = 3; l < lmax; l++) {
    double ll = static_cast<double>(l);
    fac1[l] = sqrt((2 * ll + 3) * (2 * ll + 1) / ((ll - 2) * (ll + 4) * ll * (ll + 2))) * (ll + 1);
    fac2[l] = 3.0 / (ll * (ll + 1));
    fac3[l] = sqrt((2 * ll + 3) / (2 * ll - 1) * (ll - 3) * (ll + 3) * (ll - 1) * (ll + 1) /
                   ((ll - 2) * (ll + 4) * ll * (ll + 2))) *
              (ll + 1) / ll;
    fac4[l] = sqrt(2. / (2 * ll + 3));
  }

  Nodes dl{};
  Nodes dlm1{};
  for (int index_mu = 0; index_mu < kChunk; index_mu++) {
    d3m1[0 * kChunk + index_mu] = 0;
    d3m1[1 * kChunk + index_mu] = 0;
    dlm1[index_mu]              = 0.; /*l=2*/
    d3m1[2 * kChunk + index_mu] = 0;
    dl[index_mu] = sqrt(105. / 2.) * (1 + mu[index_mu]) * (1 - mu[index_mu]) * (1 - mu[index_mu]) /
                   8.; /*l=3*/
    d3m1[3 * kChunk + index_mu] = dl[index_mu] * sqrt(2. / 7.);
  }
  /* sqrt((2l+1)/2)*d22 recurrence, supposed to be more stable */
  WignerRecurrence<WignerShift::kPlus>(mu, 3, lmax, fac1, fac2, fac3, fac4, dl, dlm1, d3m1);
}

/**
 * This routine computes the d3m3 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d3m3   Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

void D3m3(const Nodes& mu, int lmax, double* d3m3) {
  std::vector<double> fac1(lmax);
  std::vector<double> fac2(lmax);
  std::vector<double> fac3(lmax);
  std::vector<double> fac4(lmax);
  for (int l = 3; l < lmax; l++) {
    double ll = static_cast<double>(l);
    fac1[l]   = sqrt((2 * ll + 3) * (2 * ll + 1)) * (ll + 1) / ((ll - 2) * (ll + 4));
    fac2[l]   = 9.0 / (ll * (ll + 1));
    fac3[l]   = sqrt((2 * ll + 3) / (2 * ll - 1)) * (ll - 3) * (ll + 3) * (l + 1) /
                ((ll - 2) * (ll + 4) * ll);
    fac4[l]   = sqrt(2. / (2 * ll + 3));
  }

  Nodes dl{};
  Nodes dlm1{};
  for (int index_mu = 0; index_mu < kChunk; index_mu++) {
    d3m3[0 * kChunk + index_mu] = 0;
    d3m3[1 * kChunk + index_mu] = 0;
    dlm1[index_mu]              = 0.; /*l=2*/
    d3m3[2 * kChunk + index_mu] = 0;
    dl[index_mu] = sqrt(7. / 2.) * (1 - mu[index_mu]) * (1 - mu[index_mu]) * (1 - mu[index_mu]) /
                   8.; /*l=3*/
    d3m3[3 * kChunk + index_mu] = dl[index_mu] * sqrt(2. / 7.);
  }
  /* sqrt((2l+1)/2)*d22 recurrence, supposed to be more stable */
  WignerRecurrence<WignerShift::kPlus>(mu, 3, lmax, fac1, fac2, fac3, fac4, dl, dlm1, d3m3);
}

/**
 * This routine computes the d40 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d40    Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

void D40(const Nodes& mu, int lmax, double* d40) {
  std::vector<double> fac1(lmax);
  std::vector<double> fac3(lmax);
  std::vector<double> fac4(lmax);
  for (int l = 4; l < lmax; l++) {
    double ll = static_cast<double>(l);
    fac1[l]   = sqrt((2 * ll + 3) * (2 * ll + 1) / ((ll - 3) * (ll + 5)));
    fac3[l]   = sqrt((2 * ll + 3) * (ll - 4) * (ll + 4) / ((2 * ll - 1) * (ll - 3) * (ll + 5)));
    fac4[l]   = sqrt(2. / (2 * ll + 3));
  }

  Nodes dl{};
  Nodes dlm1{};
  for (int index_mu = 0; index_mu < kChunk; index_mu++) {
    d40[0 * kChunk + index_mu] = 0;
    d40[1 * kChunk + index_mu] = 0;
    d40[2 * kChunk + index_mu] = 0;
    dlm1[index_mu]             = 0.; /*l=3*/
    d40[3 * kChunk + index_mu] = 0;
    dl[index_mu] = sqrt(315.) * (1 + mu[index_mu]) * (1 + mu[index_mu]) * (1 - mu[index_mu]) *
                   (1 - mu[index_mu]) / 16.; /*l=4*/
    d40[4 * kChunk + index_mu] = dl[index_mu] * sqrt(2. / 9.);
  }
  /* sqrt((2l+1)/2)*d22 recurrence, supposed to be more stable */
  WignerRecurrence<WignerShift::kNone>(mu, 4, lmax, fac1, fac1, fac3, fac4, dl, dlm1, d40);
}

/**
 * This routine computes the d4m2 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d4m2   Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

void D4m2(const Nodes& mu, int lmax, double* d4m2) {
  std::vector<double> fac1(lmax);
  std::vector<double> fac2(lmax);
  std::vector<double> fac3(lmax);
  std::vector<double> fac4(lmax);
  for (int l = 4; l < lmax; l++) {
    double ll = static_cast<double>(l);
    fac1[l]   = sqrt((2 * ll + 3) * (2 * ll + 1) / ((ll - 3) * (ll + 5) * (ll - 1) * (ll + 3))) *
                (ll + 1.);
    fac2[l]   = 8. / (ll * (ll + 1));
    fac3[l]   = sqrt((2 * ll + 3) * (ll - 4) * (ll + 4) * (ll - 2) * (ll + 2) /
                     ((2 * ll - 1) * (ll - 3) * (ll + 5) * (ll - 1) * (ll + 3))) *
                (ll + 1) / ll;
    fac4[l]   = sqrt(2. / (2 * ll + 3));
  }

  Nodes dl{};
  Nodes dlm1{};
  for (int index_mu = 0; index_mu < kChunk; index_mu++) {
    d4m2[0 * kChunk + index_mu] = 0;
    d4m2[1 * kChunk + index_mu] = 0;
    d4m2[2 * kChunk + index_mu] = 0;
    dlm1[index_mu]              = 0.; /*l=3*/
    d4m2[3 * kChunk + index_mu] = 0;
    dl[index_mu] = sqrt(126.) * (1 + mu[index_mu]) * (1 - mu[index_mu]) * (1 - mu[index_mu]) *
                   (1 - mu[index_mu]) / 16.; /*l=4*/
    d4m2[4 * kChunk + index_mu] = dl[index_mu] * sqrt(2. / 9.);
  }
  /* sqrt((2l+1)/2)*d22 recurrence, supposed to be more stable */
  WignerRecurrence<WignerShift::kPlus>(mu, 4, lmax, fac1, fac2, fac3, fac4, dl, dlm1, d4m2);
}

/**
 * This routine computes the d4m4 term
 *
 * @param mu     Input: Vector of cos(beta) values
 * @param num_mu Input: Number of cos(beta) values
 * @param lmax   Input: maximum multipole
 * @param d4m4   Input/output: Result is stored here
 *
 * Wigner d-functions, computed by recurrence
 * actual recurrence on \f$ \sqrt{(2l+1)/2} d^l_{mm'} \f$ for stability
 * Formulae from Kostelec & Rockmore 2003
 **/

void D4m4(const Nodes& mu, int lmax, double* d4m4) {
  std::vector<double> fac1(lmax);
  std::vector<double> fac2(lmax);
  std::vector<double> fac3(lmax);
  std::vector<double> fac4(lmax);
  for (int l = 4; l < lmax; l++) {
    double ll = static_cast<double>(l);
    fac1[l]   = sqrt((2 * ll + 3) * (2 * ll + 1)) * (ll + 1) / ((ll - 3) * (ll + 5));
    fac2[l]   = 16. / (ll * (ll + 1));
    fac3[l]   = sqrt((2 * ll + 3) / (2 * ll - 1)) * (ll - 4) * (ll + 4) * (ll + 1) /
                ((ll - 3) * (ll + 5) * ll);
    fac4[l]   = sqrt(2. / (2 * ll + 3));
  }

  Nodes dl{};
  Nodes dlm1{};
  for (int index_mu = 0; index_mu < kChunk; index_mu++) {
    d4m4[0 * kChunk + index_mu] = 0;
    d4m4[1 * kChunk + index_mu] = 0;
    d4m4[2 * kChunk + index_mu] = 0;
    dlm1[index_mu]              = 0.; /*l=3*/
    d4m4[3 * kChunk + index_mu] = 0;
    dl[index_mu] = sqrt(9. / 2.) * (1 - mu[index_mu]) * (1 - mu[index_mu]) * (1 - mu[index_mu]) *
                   (1 - mu[index_mu]) / 16.; /*l=4*/
    d4m4[4 * kChunk + index_mu] = dl[index_mu] * sqrt(2. / 9.);
  }
  /* sqrt((2l+1)/2)*d22 recurrence, supposed to be more stable */
  WignerRecurrence<WignerShift::kPlus>(mu, 4, lmax, fac1, fac2, fac3, fac4, dl, dlm1, d4m4);
}

}  // namespace wigner
