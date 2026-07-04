#ifndef __QSS__
#define __QSS__

#define _MIN_NUMBER_OF_LAGUERRE_POINTS_ 5

/******************************************/
/* Quadrature Sampling Strategy for CLASS */
/* 10/12 2010                             */
/* Thomas Tram                            */
/******************************************/
#include <memory>
#include <vector>

#include "common.h"

enum quadrature_method {
  qm_auto,
  qm_Laguerre,
  qm_trapz_indefinite,
  qm_trapz,
  qm_GB_Laguerre,
  qm_FermiDirac,
};

/** Optional grey-body parameters threaded into manual quadrature. When
 *  `active` is false the GB-specific methods are not used. Defined here (not in
 *  a species header) so quadrature.cpp stays species-agnostic. Distinct from
 *  greybody::GreyBodyParams (the moment-solver class) — this is just the few
 *  coefficients the quadrature engine needs. */
struct GBQuadParams {
  bool active          = false;
  double alpha         = 0.;
  double x_times_alpha = 0.;
};

/* Structures for QSS */

typedef struct adaptive_integration_tree_node {
  /* binary tree node: */
  double I;              /* Estimate of integral */
  double err;            /* Estimated error */
  std::vector<double> x; /* Pointer to the abscissas of node */
  std::vector<double> w; /* Pointer to the corresponding weights */
  int leaf_childs;       /* Number of leafs under current node. 1 means that the node is a leaf. */
  /* Pointer to children: */
  std::unique_ptr<struct adaptive_integration_tree_node> left, right; /* Pointer to left child. */
} qss_node;

void get_qsampling(double* x,
                   double* w,
                   int* N,
                   int N_max,
                   double rtol,
                   double* qvec,
                   int qsiz,
                   void (*test)(void* params_for_function, double q, double* psi),
                   void (*function)(void* params_for_function, double q, double* f0),
                   void* params_for_function);
bool get_qsampling_manual(double* x,
                          double* w,
                          double* dq,
                          int N,
                          double qmax,
                          enum quadrature_method method,
                          double* qvec,
                          int qsiz,
                          void (*function)(void* params_for_function, double q, double* f0),
                          void* params_for_function,
                          const GBQuadParams& gb);
void sort_x_and_w(double* x, double* w, double* workx, double* workw, int startidx, int endidx);
void get_leaf_x_and_w(qss_node* node, int* ind, double* x, double* w, int isindefinite);
void reduce_tree(qss_node* node, int level);
void leaf_count(qss_node* node);
double get_integral(qss_node* node, int level);
void gk_adapt(std::unique_ptr<qss_node>& node,
              void (*test)(void* params_for_function, double q, double* psi),
              void (*function)(void* params_for_function, double q, double* f0),
              void* params_for_function,
              double tol,
              int treemode,
              double a,
              double b,
              int isindefinite);
void compute_Laguerre(
    double* x, double* w, int N, double alpha, double* b, double* c, int totalweight);

/** Compute the N-point Gaussian rule for the weight 1 / (exp(x) + 1) on
 *  [0, infinity). The rule is constructed from Laguerre modified moments
 *  using the modified-Chebyshev algorithm and Golub-Welsch diagonalization.
 *  It is reliable through order 17; higher orders are intentionally rejected
 *  rather than silently returning an ill-conditioned rule. */
bool compute_FermiDirac(double* x, double* w, int N);
void gk_quad(void (*test)(void* params_for_function, double q, double* psi),
             void (*function)(void* params_for_function, double q, double* f0),
             void* params_for_function,
             qss_node* node,
             double a,
             double b,
             int isindefinite);

void quadrature_gauss_legendre(double* mu, double* w8, int n, double tol);

#endif
