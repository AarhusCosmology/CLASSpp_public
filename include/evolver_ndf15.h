#ifndef __EVO__
#define __EVO__
#include "common.h"
// #include "perturbations.h"
#include <memory>
#include <vector>

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
  int use_sparse;
  int sparse_stuff_initialized;
  int max_nonzero; /*Maximal number of non-zero entries to be considered sparse */
  int repeated_pattern;
  int trust_sparse; /* Number of times a pattern is repeated (actually included) before we trust it. */
  int has_grouping;
  int has_pattern;
  int new_jacobian; /* True if sp_ludcmp has not been run on the current jacobian. */
  int cnzmax;
  std::vector<int> col_group;  /* Column grouping. Groups go from 0 to max_group*/
  std::vector<int> col_wi;     /* Workarray for column grouping*/
  int max_group;               /*Number of columngroups -1 */
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

/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

int initialize_jacobian(struct jacobian* jac, int neq);
int uninitialize_jacobian(struct jacobian* jac);
int initialize_numjac_workspace(struct numjac_workspace* nj_ws, int neq);
int uninitialize_numjac_workspace(struct numjac_workspace* nj_ws);
int calc_C(struct jacobian* jac);
int interp_from_dif(double tinterp,
                    double tnew,
                    double* ynew,
                    double h,
                    double** dif,
                    int k,
                    double* yinterp,
                    double* ypinterp,
                    double* yppinterp,
                    int* index,
                    int neq,
                    int output);
int new_linearisation(struct jacobian* jac, double hinvGak, int neq);
int adjust_stepsize(double** dif, double abshdivabshlast, int neq, int k);
void eqvec(double* datavec, double* emptyvec, int n);
int lubksb(double** a, int n, int* indx, double b[]);
int ludcmp(double** a, int n, int* indx, double* d, double* vv);
int fzero_Newton(int (*func)(double* x, int x_size, void* param, double* F),
                 double* x_inout,
                 double* dxdF,
                 int x_size,
                 double tolx,
                 double tolF,
                 void* param,
                 int* fevals);

int numjac(int (*derivs)(double x, double* y, double* dy, void* parameters_and_workspace),
           double t,
           double* y,
           double* fval,
           struct jacobian* jac,
           struct numjac_workspace* nj_ws,
           double thresh,
           int neq,
           int* nfe,
           void* parameters_and_workspace_for_derivs);

int evolver_ndf15(
    int (*derivs)(double x, double* y, double* dy, void* parameters_and_workspace),
    double x_ini,
    double x_final,
    double* y_inout,
    int* used_in_output,
    int neq,
    void* parameters_and_workspace_for_derivs,
    double rtol,
    double minimum_variation,
    int (*timescale_and_approximation)(double x,
                                       void* parameters_and_workspace,
                                       double* timescales),
    double timestep_over_timescale,
    double* t_vec,
    int t_res,
    int (*output)(double x, double y[], double dy[], int index_x, void* parameters_and_workspace),
    int (*print_variables)(double x, double y[], double dy[], void* parameters_and_workspace));

#ifdef __cplusplus
}
#endif

/**************************************************************/

#endif
