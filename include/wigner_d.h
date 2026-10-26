#ifndef WIGNER_D_H
#define WIGNER_D_H

#include <array>

/**
 * Wigner d-functions d^l_{mm'}(mu) for l = 0..lmax, as the lensing module needs them for its
 * correlation functions, at kChunk nodes mu at once. The result is stored as
 * d[l * kChunk + node]: the recurrences in l run across the nodes, which keeps the pipeline and
 * the vector lanes busy, and the sums over l that use them run across the nodes as well.
 * Recurrences on sqrt((2l+1)/2) d^l_{mm'} after Kostelec & Rockmore (2003), for stability.
 * docs/superpowers/specs/2026-10-26-lensing-wigner-recursions-design.md
 */
namespace wigner {

/** nodes handled together */
constexpr int kChunk = 32;
using Nodes          = std::array<double, kChunk>;

/** d must hold (lmax + 1) * kChunk values */
void D00(const Nodes& mu, int lmax, double* d00);
void D11(const Nodes& mu, int lmax, double* d11);
void D1m1(const Nodes& mu, int lmax, double* d1m1);
void D2m2(const Nodes& mu, int lmax, double* d2m2);
void D22(const Nodes& mu, int lmax, double* d22);
void D20(const Nodes& mu, int lmax, double* d20);
void D31(const Nodes& mu, int lmax, double* d31);
void D3m1(const Nodes& mu, int lmax, double* d3m1);
void D3m3(const Nodes& mu, int lmax, double* d3m3);
void D40(const Nodes& mu, int lmax, double* d40);
void D4m2(const Nodes& mu, int lmax, double* d4m2);
void D4m4(const Nodes& mu, int lmax, double* d4m4);

}  // namespace wigner

#endif  // WIGNER_D_H
