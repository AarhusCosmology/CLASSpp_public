#ifndef __EVO__
#define __EVO__
#include "common.h"
// #include "perturbations.h"
#include <memory>
#include <vector>

#include "evolver_options.h"
#include "sparse.h"
#define TINY 1e-50
/**************************************************************/

struct jacobian {
  /*Stuff for normal method: */
  std::vector<double*> dfdy;  /* row pointers into dfdy_data_vec */
  std::vector<double> jacvec; /*Stores experience gained from subsequent calls */
  std::vector<double*> LU;    /* row pointers into LU_data_vec */
  std::vector<double> LUw;
  std::vector<int> luidx;
  /*Sparse stuff:*/
  int use_sparse               = 0;
  int sparse_stuff_initialized = 0;
  int max_nonzero              = 0; /*Maximal number of non-zero entries to be considered sparse */
  int repeated_pattern         = 0;
  int trust_sparse =
      0; /* Number of times a pattern is repeated (actually included) before we trust it. */
  int has_grouping = 0;
  int has_pattern  = 0;
  int new_jacobian = 0; /* True if sp_ludcmp has not been run on the current jacobian. */
  int cnzmax       = 0;
  std::vector<int> col_group;  /* Column grouping. Groups go from 0 to max_group*/
  std::vector<int> col_wi;     /* Workarray for column grouping*/
  int max_group = 0;           /*Number of columngroups -1 */
  std::unique_ptr<sp_mat> spJ; /* Stores the matrix we want to decompose */
  std::vector<double> xjac;    /*Stores the values of the sparse jacobian. (Same pattern as spJ) */
  std::unique_ptr<sp_num> Numerical; /*Stores the LU decomposition.*/
  std::vector<int> Cp; /* Stores the column pointers of the spJ+spJ' sparsity pattern. */
  std::vector<int> Ci; /* Stores the row indices of the  spJ+spJ' sparsity pattern. */

  /* Contiguous flat backing for the 2D row-pointer arrays above */
  std::vector<double> dfdy_data_vec;
  std::vector<double> LU_data_vec;
};

struct numjac_workspace {
  std::vector<double> yscale;
  std::vector<double> del;
  std::vector<double> Difmax;
  std::vector<double> absFdelRm;
  std::vector<double> absFvalue;
  std::vector<double> absFvalueRm;
  std::vector<double> Fscale;
  std::vector<double> ffdel;
  std::vector<double> yydel;
  std::vector<double> tmp;

  std::vector<double*> ydel_Fdel; /* row pointers into ydel_Fdel_data_vec */

  std::vector<int> logj;
  std::vector<int> Rowmax;

  /* Contiguous flat backing for ydel_Fdel */
  std::vector<double> ydel_Fdel_data_vec;
};

void initialize_jacobian(struct jacobian* jac, int neq);
void uninitialize_jacobian(struct jacobian* jac);
void initialize_numjac_workspace(struct numjac_workspace* nj_ws, int neq);
void uninitialize_numjac_workspace(struct numjac_workspace* nj_ws);
void calc_C(struct jacobian* jac);
void interp_from_dif(double tinterp,
                     double tnew,
                     double* ynew,
                     double h,
                     double** dif,
                     int k,
                     double* yinterp,
                     double* ypinterp,
                     double* yppinterp,
                     int* mask,
                     int neq,
                     int output);
void new_linearisation(struct jacobian* jac, double hinvGak, int neq);
void adjust_stepsize(double** dif, double abshdivabshlast, int neq, int k);
void eqvec(double* datavec, double* emptyvec, int n);
void lubksb(double** a, int n, int* indx, double b[]);
bool ludcmp(double** a, int n, int* indx, double* d, double* vv);
void fzero_Newton(void (*func)(double* x, int x_size, void* param, double* F),
                  double* x_inout,
                  double* dxdF,
                  int x_size,
                  double tolx,
                  double tolF,
                  void* param,
                  int* fevals);

void numjac(void (*derivs)(double x, double* y, double* dy, void* parameters_and_workspace),
            double t,
            double* y,
            double* fval,
            struct jacobian* jac,
            struct numjac_workspace* nj_ws,
            double thresh,
            int neq,
            int* nfe,
            void* parameters_and_workspace_for_derivs);

/** Variable-order (1-5) BDF evolver with an adaptive step, for stiff systems.
 *
 *  Honours EvolverOptions::max_order and ::stats in addition to the fields every
 *  evolver honours. Setting anything else is an error rather than an ignore --
 *  see evolver_options.h. */
void evolver_ndf15(EvolverDerivs derivs,
                   double x_ini,
                   double x_final,
                   double* y_inout,
                   int neq,
                   void* parameters_and_workspace_for_derivs,
                   const EvolverOptions& options);

/**************************************************************/

#endif
