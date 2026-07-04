/******************************************/
/* Quadrature Sampling Strategy for CLASS */
/* 10/12 2010                             */
/* Thomas Tram                            */
/******************************************/
#include "quadrature.h"

#include <algorithm>
#include <limits>
#define _N_COMB_LAG_ 16

namespace {

constexpr int kFermiDiracMaxOrder   = 17;
constexpr int kFermiDiracCohenTerms = 40;

double fermi_dirac_weight(double x) {
  return 1. / (std::exp(x) + 1.);
}

/**
 * Compute \int_0^infinity L_n(x)/(exp(x)+1) dx, where L_n is the ordinary
 * Laguerre polynomial. Expanding the Fermi-Dirac weight gives the alternating
 * series sum_{k>=1} (-1)^(k-1) (1-1/k)^n/k. Cohen--Rodriguez Villegas--Zagier
 * acceleration makes the evaluation accurate enough for the subsequent
 * modified-Chebyshev recurrence.
 */
double fermi_dirac_modified_moment(int n) {
  const long double root = 3.L + std::sqrt(8.L);
  long double d          = std::pow(root, kFermiDiracCohenTerms);
  d                      = 0.5L * (d + 1.L / d);
  long double b          = -1.L;
  long double c          = -d;
  long double sum        = 0.L;

  for (int k = 0; k < kFermiDiracCohenTerms; ++k) {
    c                        = b - c;
    const long double index  = static_cast<long double>(k + 1);
    const long double term   = std::pow(1.L - 1.L / index, n) / index;
    sum                     += c * term;
    b *= static_cast<long double>((k + kFermiDiracCohenTerms) * (k - kFermiDiracCohenTerms)) /
         ((static_cast<long double>(k) + 0.5L) * (k + 1));
  }
  return static_cast<double>(sum / d);
}

/** Implicit-QL diagonalization of a real symmetric tridiagonal matrix. The
 *  eigenvectors are stored by columns in eigenvectors. */
bool diagonalize_symmetric_tridiagonal(std::vector<double>& diagonal,
                                       std::vector<double>& subdiagonal,
                                       std::vector<double>& eigenvectors) {
  const int n = static_cast<int>(diagonal.size());
  if (n == 0)
    return false;

  // The QL iteration expects the subdiagonal shifted one element to the left.
  for (int i = 1; i < n; ++i)
    subdiagonal[i - 1] = subdiagonal[i];
  subdiagonal[n - 1] = 0.;

  constexpr int kMaxIterations = 100;
  for (int l = 0; l < n; ++l) {
    int iterations = 0;
    while (true) {
      int m = l;
      for (; m < n - 1; ++m) {
        const double scale = std::abs(diagonal[m]) + std::abs(diagonal[m + 1]);
        if (std::abs(subdiagonal[m]) <= std::numeric_limits<double>::epsilon() * scale)
          break;
      }
      if (m == l)
        break;
      if (++iterations > kMaxIterations)
        return false;

      double g = (diagonal[l + 1] - diagonal[l]) / (2. * subdiagonal[l]);
      double r = std::hypot(g, 1.);
      g        = diagonal[m] - diagonal[l] + subdiagonal[l] / (g + std::copysign(r, g));
      double s = 1.;
      double c = 1.;
      double p = 0.;

      int i = m - 1;
      for (; i >= l; --i) {
        const double f     = s * subdiagonal[i];
        const double b     = c * subdiagonal[i];
        r                  = std::hypot(f, g);
        subdiagonal[i + 1] = r;
        if (r == 0.) {
          diagonal[i + 1] -= p;
          subdiagonal[m]   = 0.;
          break;
        }
        s               = f / r;
        c               = g / r;
        g               = diagonal[i + 1] - p;
        r               = (diagonal[i] - g) * s + 2. * c * b;
        p               = s * r;
        diagonal[i + 1] = g + p;
        g               = c * r - b;

        for (int row = 0; row < n; ++row) {
          const double z                = eigenvectors[row * n + i + 1];
          eigenvectors[row * n + i + 1] = s * eigenvectors[row * n + i] + c * z;
          eigenvectors[row * n + i]     = c * eigenvectors[row * n + i] - s * z;
        }
      }
      if (r == 0. && i >= l)
        continue;
      diagonal[l]    -= p;
      subdiagonal[l]  = g;
      subdiagonal[m]  = 0.;
    }
  }

  // Golub-Welsch rules conventionally expose nodes in increasing order.
  for (int i = 0; i < n - 1; ++i) {
    int min_index = i;
    for (int j = i + 1; j < n; ++j) {
      if (diagonal[j] < diagonal[min_index])
        min_index = j;
    }
    if (min_index == i)
      continue;
    std::swap(diagonal[i], diagonal[min_index]);
    for (int row = 0; row < n; ++row)
      std::swap(eigenvectors[row * n + i], eigenvectors[row * n + min_index]);
  }
  return true;
}

}  // namespace

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
                          const GBQuadParams& gb) {
  double y, h, t;
  switch (method) {
    case (qm_auto):
      return false;
    case (qm_Laguerre): {
      /* Allocate storage for Laguerre coefficients: */
      std::vector<double> b(N);
      std::vector<double> c(N);
      compute_Laguerre(x, dq, N, 0.0, b.data(), c.data(), true);
      for (int i = 0; i < N; i++) {
        (*function)(params_for_function, x[i], &y);
        w[i] = dq[i] * y;
      }
      return true;
    }
    case (qm_FermiDirac): {
      std::vector<double> fd_weight(N);
      class_test(!compute_FermiDirac(x, fd_weight.data(), N),
                 "qm_FermiDirac supports 1 through %d momentum bins; the modified-moment "
                 "recurrence is ill-conditioned above this range",
                 kFermiDiracMaxOrder);
      for (int i = 0; i < N; ++i) {
        const double reference_weight = fermi_dirac_weight(x[i]);
        class_test(reference_weight == 0.,
                   "qm_FermiDirac produced an abscissa outside the representable Fermi-Dirac "
                   "weight range");
        dq[i] = fd_weight[i] / reference_weight;
        (*function)(params_for_function, x[i], &y);
        w[i] = dq[i] * y;
      }
      return true;
    }
    case (qm_GB_Laguerre): {
      // Generalized Gauss-Laguerre on the GL weight q^(alpha+1) e^(-alpha x q),
      // i.e. the integrand q^2 * f0 with the grey-body f0 ~ q^(alpha-1) e^(-alpha x q).
      // The rule divides nodes by gb.x_times_alpha, so it is only meaningful for a
      // grey-body species. Fail loudly if it was selected without valid GB params
      // (e.g. quadrature_strategy=qm_GB_Laguerre on a non-grey-body species) rather
      // than dividing by zero into NaN abscissas.
      class_test(!gb.active || gb.x_times_alpha <= 0.,
                 "qm_GB_Laguerre requires grey-body parameters (active alpha*x > 0); it was "
                 "selected for a species that did not provide them");
      std::vector<double> b(N);
      std::vector<double> c(N);
      compute_Laguerre(x, dq, N, gb.alpha + 1., b.data(), c.data(), true);
      double scaling = 1. / gb.x_times_alpha;
      for (int i = 0; i < N; i++) {
        dq[i] *= pow(x[i], -(gb.alpha + 1.)) * scaling;  // bare weight (uses original node)
        x[i]  *= scaling;                                // rescale node onto the alpha*x grid
        (*function)(params_for_function, x[i], &y);      // y = f0(x[i])
        w[i] = dq[i] * y;
      }
      return true;
    }
    case (qm_trapz):
      for (int i = 0; i < N; ++i) {
        /** Note that we count q=0 as an extra point with weight 0 */
        h    = qmax / N;
        x[i] = h + i * h;
        (*function)(params_for_function, x[i], &y);
        dq[i] = h;
        w[i]  = y * h;
        if (i == N - 1) {
          dq[i] *= 0.5;
          w[i]  *= 0.5;
        }
      }
      return true;
    case (qm_trapz_indefinite):
      /** We do the variable transformation q = 1/t-1. The trapezoidal rule is closed, but since the distribution function
	goes to zero in both limits, we can use an effectively N+2 rule simply by not using the exterior points. */
      for (int i = 0; i < N; ++i) {
        h    = 1.0 / (N + 1.0);
        t    = h + i * h;
        x[i] = 1.0 / t - 1.0;
        (*function)(params_for_function, x[i], &y);
        dq[i] = h / t / t;
        w[i]  = dq[i] * y;
      }
      return true;
  }
  return true;
}
void get_qsampling(double* x,
                   double* w,
                   int* N,
                   int N_max,
                   double rtol,
                   double* qvec,
                   int qsiz,
                   void (*test)(void* params_for_function, double q, double* psi),
                   void (*function)(void* params_for_function, double q, double* f0),
                   void* params_for_function) {
  /* This routine returns the fewest possible number of abscissas and weights under
     the requirement that a test function folded with the neutrino distribution function
     can be integrated to an accuracy of rtol. If the distribution function is Fermi-Dirac
     or close, a Laguerre quadrature formula is often the best choice.

     This function combines two completely different strategies: Adaptive Gauss-Kronrod
     quadrature and Laguerres quadrature formula. */

  int i, NL = 2, NR, level, Nadapt = 0, NLag, NLag_max, Nold = NL;
  int NFD, NFD_max;
  int adapt_converging = false, Laguerre_converging = false, FermiDirac_converging = false,
      combined_converging = false;
  double y, y1, y2, I, Igk, err, ILag, IFD;
  std::unique_ptr<qss_node> root      = nullptr;
  std::unique_ptr<qss_node> root_comb = nullptr;
  double I_comb, I_atzero, I_atinf, I_comb2;
  int N_comb = 0, N_comb_leg = 4;
  double a_comb, b_comb;
  double q_leg[4], w_leg[4];
  double q_lag[_N_COMB_LAG_], w_lag[_N_COMB_LAG_];
  std::string method_chosen;
  double qmin = 0., qmax = 0., qmaxm1 = 0.;
  std::vector<double> wcomb2;
  double delq;
  double Itot  = 0.0;
  int zeroskip = 0;

  /* Set roots and weights for Gauss-Legendre 4 point rule: */
  q_leg[3] = sqrt((3.0 + 2.0 * sqrt(6.0 / 5.0)) / 7.0);
  q_leg[2] = sqrt((3.0 - 2.0 * sqrt(6.0 / 5.0)) / 7.0);
  q_leg[1] = -q_leg[2];
  q_leg[0] = -q_leg[3];
  w_leg[3] = (18.0 - sqrt(30.0)) / 36.0;
  w_leg[2] = (18.0 + sqrt(30.0)) / 36.0;
  w_leg[1] = w_leg[2];
  w_leg[0] = w_leg[3];

  /* Allocate storage for Laguerre coefficients: */
  std::vector<double> b(N_max);
  std::vector<double> c(N_max);
  /* If a vector of q values has been passed, use it: */
  if ((qvec != nullptr) && (qsiz > 1)) {
    qmin   = qvec[0];
    qmax   = qvec[qsiz - 1];
    qmaxm1 = qvec[qsiz - 2];
  }
  else {
    qvec = nullptr;
  }

  /* First do the adaptive quadrature - this will also give the value of the integral: */
  gk_adapt(root, (*test), (*function), params_for_function, rtol * 1e-4, 1, 0.0, 1.0, true);
  /* Do a leaf count: */
  leaf_count(root.get());
  /* I can get the integral now: */
  I = get_integral(root.get(), 1);
  //printf("I = %le |%le, used points: %d\n",I,I-1.0,15*root->leaf_childs);
  Itot += I;

  /* Starting from the top, move down in levels until tolerance is met: */
  for (level = root->leaf_childs; level >= 1; level--) {
    Igk = get_integral(root.get(), level);
    err = I - Igk;
    if (fabs(err / Itot) < rtol)
      break;
  }
  if (level > 0) {
    /* Reduce tree to the found level:*/
    reduce_tree(root.get(), level);
    /* Count the new leafs: */
    leaf_count(root.get());
    /* I know know how many function evaluations is
       required by the adaptively sampled grid:*/
    Nadapt = 15 * root->leaf_childs;
    /* The adaptive routine could not recieve required precision
       using less than the required maximal number of points.*/
    if (Nadapt <= N_max)
      adapt_converging = true;
  }

  /* Combined adaptive quadrature and Laguerre rescaled quadrature: */
  if (qvec != nullptr) {
    /* Evaluate [0;qmin] using 4 point Gauss-Legendre rule: */
    (*function)(params_for_function, qmin, &y2);
    for (i = 0, I_atzero = 0.0; i < N_comb_leg; i++) {
      q_leg[i] = 0.5 * qmin * (1.0 + q_leg[i]);
      w_leg[i] = w_leg[i] * 0.5 * qmin * y2;
      (*test)(params_for_function, q_leg[i], &y);
      I_atzero += w_leg[i] * y;
    }

    /* Find asymptotic extrapolation:*/
    (*function)(params_for_function, qmaxm1, &y1);
    (*function)(params_for_function, qmax, &y2);

    b_comb = (y1 / y2 - 1.0) / (qmax - qmaxm1);
    b_comb = std::max(b_comb, 1e-100);
    //c_comb = -b_comb*qmax;
    a_comb = y2 * exp(b_comb * qmax);
    // printf("f(q) = %g*exp(-%g*q) \n",a_comb,b_comb);
    //(*function)(params_for_function,100,&y2);
    //printf("f(100) = %e ?= %e\n",y2,a_comb*exp(-b_comb*100));

    /* Evaluate tail using 6 point Laguerre: */
    compute_Laguerre(q_lag, w_lag, _N_COMB_LAG_, 0.0, b.data(), c.data(), true);
    for (i = 0, I_atinf = 0.0; i < _N_COMB_LAG_; i++) {
      w_lag[i] *= exp(-q_lag[i]);
      q_lag[i]  = qmax + q_lag[i] / b_comb;
      w_lag[i]  = a_comb / b_comb * exp(-b_comb * qmax) * w_lag[i];
      (*test)(params_for_function, q_lag[i], &y);
      I_atinf += w_lag[i] * y;
    }

    /* Do the adaptive quadrature - this will also give the main part of the integral: */
    gk_adapt(root_comb,
             (*test),
             (*function),
             params_for_function,
             rtol * 1e-2,
             1,
             qmin,
             qmax,
             false);
    /* Do a leaf count: */
    leaf_count(root_comb.get());
    /* Starting from the top, move down in levels until tolerance is met: */
    for (level = root_comb->leaf_childs; level >= 1; level--) {
      I_comb = get_integral(root_comb.get(), level);
      //printf("%le + %le + %le = %le | %le\n",
      //     I_atzero,I_atinf,I_comb,I_comb+I_atinf+I_atzero,I_comb+I_atinf+I_atzero-1.0);
      I_comb += (I_atinf + I_atzero);
      err     = I - I_comb;
      if (fabs(err / Itot) < rtol)
        break;
    }
    /* Reduce tree to the found level:*/
    if (level > 0) {
      reduce_tree(root_comb.get(), level);
      /* Count the new leafs: */
      leaf_count(root_comb.get());
      N_comb = 15 * root_comb->leaf_childs + _N_COMB_LAG_ + N_comb_leg;
      if (N_comb <= N_max)
        combined_converging = true;
    }

    /* Do the second combined quadrature: Same as above, but with trapezoidal rule
       using the given q grid:  */
    wcomb2.resize(qsiz);
    I_comb2 = 0.0;
    for (i = 0; i < qsiz; i++) {
      if (i == 0) {
        delq = qvec[1] - qvec[0];
      }
      else if (i == qsiz - 1) {
        delq = qvec[qsiz - 1] - qvec[qsiz - 2];
      }
      else {
        delq = qvec[i + 1] - qvec[i - 1];
      }
      (*function)(params_for_function, qvec[i], &y2);
      wcomb2[i] = 0.5 * y2 * delq;
      (*test)(params_for_function, qvec[i], &y);
      I_comb2 += wcomb2[i] * y;
    }
    I_comb2 += (I_atzero + I_atinf);
    err      = I - I_comb2;
    //    if(fabs(err/Itot)<rtol) combined2_converging= true;
    //printf("I_comb2 = %e, rerr = %e\n",I_comb2,fabs(err/I));
  }

  /* Search for the minimal Laguerre quadrature rule: */
  NLag_max = std::min(N_max, 80);
  for (NLag = NL; NLag <= NLag_max; NLag = std::min(NLag_max, NLag + 10)) {
    /* Evaluate integral: */
    compute_Laguerre(x, w, NLag, 0.0, b.data(), c.data(), true);
    ILag = 0.0;
    for (i = 0; i < NLag; i++) {
      (*test)(params_for_function, x[i], &y);
      (*function)(params_for_function, x[i], &y2);
      w[i] *= y2;
      ILag += y * w[i];
    }
    err = I - ILag;
    //fprintf(stderr,"\n Computing Laguerre, N=%d, I=%g and err=%g.\n",NLag,ILag,err);
    if (fabs(err / I) < rtol) {
      Laguerre_converging = true;
      break;
    }
    Nold = NLag;
    if (Nold == NLag_max)
      break;
  }

  if (Laguerre_converging) {
    /* We must refine NLag: */
    NL = Nold;
    NR = NLag;
    while ((NR - NL) > 1) {
      NLag = (NL + NR) / 2;
      /* Evaluate integral: */
      compute_Laguerre(x, w, NLag, 0.0, b.data(), c.data(), true);
      ILag = 0.0;
      for (i = 0; i < NLag; i++) {
        (*test)(params_for_function, x[i], &y);
        (*function)(params_for_function, x[i], &y2);
        w[i] *= y2;
        ILag += y * w[i];
      }
      err = I - ILag;
      //fprintf(stderr,"\n NLag=%d, rerr=%g.\n",NLag,fabs(err/I));
      if (fabs(err / Itot) < rtol) {
        NR = NLag;
      }
      else {
        NL = NLag;
      }
    }
  }

  /* Search for the minimal Fermi-Dirac Gaussian quadrature rule. The rule is
     constructed once per trial order from modified Laguerre moments; it is
     exact for the FD-weighted test polynomial through degree 2*N-1. As for
     Gauss-Laguerre, fold the actual PSD into the stored weights so this remains
     a valid importance rule for arbitrary analytic or tabulated PSDs. */
  NFD_max = std::min(N_max, kFermiDiracMaxOrder);
  for (NFD = 2; NFD <= NFD_max; ++NFD) {
    std::vector<double> x_fd(NFD);
    std::vector<double> w_fd(NFD);
    if (!compute_FermiDirac(x_fd.data(), w_fd.data(), NFD))
      break;

    IFD = 0.;
    for (i = 0; i < NFD; ++i) {
      (*test)(params_for_function, x_fd[i], &y);
      (*function)(params_for_function, x_fd[i], &y2);
      w_fd[i] *= y2 / fermi_dirac_weight(x_fd[i]);
      IFD     += y * w_fd[i];
    }
    err = I - IFD;
    if (fabs(err / Itot) < rtol) {
      FermiDirac_converging = true;
      break;
    }
  }

  /* Choose best method if both works: */
  *N = N_max;
  //Laguerre_converging = false;
  if (adapt_converging) {
    *N = Nadapt;
  }
  if (combined_converging) {
    if (N_comb <= *N) {
      *N               = N_comb;
      adapt_converging = false;
    }
    else {
      combined_converging = false;
    }
  }
  if (Laguerre_converging) {
    if (NLag <= *N) {
      *N                  = NLag;
      combined_converging = false;
      adapt_converging    = false;
    }
    else {
      Laguerre_converging = false;
    }
  }
  if (FermiDirac_converging) {
    if (NFD <= *N) {
      *N                  = NFD;
      Laguerre_converging = false;
      combined_converging = false;
      adapt_converging    = false;
    }
    else {
      FermiDirac_converging = false;
    }
  }
  //printf("N_adapt=%d, N_combined=%d at level=%d, Nlag=%d\n",Nadapt,N_comb,level,NLag);
  if (adapt_converging) {
    method_chosen = "Adaptive Gauss-Kronrod Quadrature";
    /* Gather weights and xvalues from tree: */
    i = Nadapt - 1;
    get_leaf_x_and_w(root.get(), &i, x, w, true);
  }
  else if (Laguerre_converging) {
    method_chosen = "Gauss-Laguerre Quadrature";
    /* x and w is already populated in this case. */
  }
  else if (FermiDirac_converging) {
    method_chosen = "Fermi-Dirac Gaussian Quadrature";
    class_test(!compute_FermiDirac(x, w, *N),
               "Could not reconstruct the selected Fermi-Dirac quadrature rule");
    for (i = 0; i < *N; ++i) {
      (*function)(params_for_function, x[i], &y);
      w[i] *= y / fermi_dirac_weight(x[i]);
    }
  }
  else if (combined_converging) {
    method_chosen = "Combined Quadrature";
    for (i = 0; i < N_comb_leg; i++) {
      x[i] = q_leg[i];
      w[i] = w_leg[i];
    }
    i = N_comb_leg;
    get_leaf_x_and_w(root_comb.get(), &i, x, w, false);
    //printf("from %d to %d\n",N_comb_leg,i);

    for (i = 0; i < _N_COMB_LAG_; i++) {
      x[N_comb - _N_COMB_LAG_ + i] = q_lag[i];
      w[N_comb - _N_COMB_LAG_ + i] = w_lag[i];
    }
  }
  else {
    /* Failed to converge! */
    class_stop(
        "get_qsampling fails to obtain a relative tolerance of %g as required using at most "
        "%d points. If the PSD is interpolated from a file, try increasing the resolution "
        "and the q-interval (qmin;qmax) if possible, or decrease tol_ncdm/tol_ncdm_bg. As a "
        "last resort, increase _QUADRATURE_MAX_/_QUADRATURE_MAX_BG_.",
        rtol,
        N_max);
  }
  /* Trim weights to avoid zero weights: */
  for (i = 0, zeroskip = 0; i < *N; i++) {
    for (; (i < *N) && (w[i + zeroskip] == 0.0); zeroskip++, (*N)--)
      ;
    x[i] = x[i + zeroskip];
    w[i] = w[i + zeroskip];
  }

  //printf("Chosen sampling: %s, with %d points.\n", method_chosen.c_str(), *N);
  //for(i=0; i<*N; i++) printf("(q,w) = (%g,%g)\n",x[i],w[i]);
  /* Deallocate tree handled by unique_ptr automatically. */
}

void sort_x_and_w(double* x, double* w, double* workx, double* workw, int startidx, int endidx) {
  int i, top = endidx, bot = startidx;
  double pivot;
  /* End recursion if only one element left in array: */
  if ((endidx - startidx) < 1) {
    return;
  }
  else {
    /*Copy x and w to workarray: */
    for (i = startidx; i <= endidx; i++) {
      workx[i] = x[i];
      workw[i] = w[i];
    }
    pivot = x[endidx];
    //printf("pivot chosen: x[%d] = %g\n",endidx,pivot);
    for (i = startidx; i < endidx; i++) {
      if (workx[i] <= pivot) {
        //printf("<--%g  ",workx[i]);
        x[bot]   = workx[i];
        w[bot++] = workw[i];
      }
      else {
        //printf("  %g-->",workx[i]);
        x[top]   = workx[i];
        w[top--] = workw[i];
      }
    }
    //printf("\n top=%d, bot=%d, left=%d, right=%d\n",top,bot,startidx,endidx);
    x[top] = pivot;
    w[top] = workw[endidx];
    /* Recursive call: */
    sort_x_and_w(x, w, workx, workw, startidx, bot - 1);
    sort_x_and_w(x, w, workx, workw, top + 1, endidx);
  }
}

void get_leaf_x_and_w(qss_node* node, int* ind, double* x, double* w, int isindefinite) {
  /* x and w should be exactly 15*root_node->leafchilds, and a leaf count should have
     been performed. Or perhaps I just use the fact that a leaf won't have children.
     Nah, let me use the leaf-count then. */
  int k;
  if (node->leaf_childs == 1) {
    for (k = 0; k < 15; k++) {
      x[*ind] = node->x[k];
      w[*ind] = node->w[k];
      if (isindefinite) {
        (*ind)--;
      }
      else {
        (*ind)++;
      }
    }
  }
  else {
    /* Do recursive call: */
    get_leaf_x_and_w(node->left.get(), ind, x, w, isindefinite);
    get_leaf_x_and_w(node->right.get(), ind, x, w, isindefinite);
  }
}

void reduce_tree(qss_node* node, int level) {
  /* Reduce the tree to a given level. Make all nodes with
     node->leaf_childs==level into leafs.
     If we call reduce_tree(root,1), nothing happens.*/
  if (node->leaf_childs == level) {
    node->left.reset();
    node->right.reset();
  }
  else if (node->leaf_childs > level) {
    /* else try to see if children nodes can be simplified: */
    reduce_tree(node->left.get(), level);
    reduce_tree(node->right.get(), level);
  }
  /* If called on a node which has leaf_childs<level, it does nothing. */
}

void leaf_count(qss_node* node) {
  /* Count the amount of leafs under a given node and write the number in the node. */
  /* We call recursively, until a node is a leaf - then we add the numbers on our
     way back:*/
  if (node->left) {
    /* This is not a leaf, do recursive call: */
    leaf_count(node->left.get());
    leaf_count(node->right.get());
    node->leaf_childs = node->left->leaf_childs + node->right->leaf_childs;
  }
  else {
    /* This is a leaf, by definition leaf_childs = 1: */
    node->leaf_childs = 1;
  }
}

double get_integral(qss_node* node, int level) {
  /* Traverse the tree and return the estimate of the integral at a given level.
     level 1 is the best estimate. */
  double IL, IR;
  /* An updated leaf_count is assumed. */
  if (node->leaf_childs <= level) {
    return node->I;
  }
  else {
    IL = get_integral(node->left.get(), level);
    IR = get_integral(node->right.get(), level);
    /* Combine the integrals: */
    return (IL + IR);
  }
}

void gk_adapt(std::unique_ptr<qss_node>& node,
              void (*test)(void* params_for_function, double q, double* psi),
              void (*function)(void* params_for_function, double q, double* f0),
              void* params_for_function,
              double tol,
              int treemode,
              double a,
              double b,
              int isindefinite) {
  /* Do adaptive Gauss-Kronrod quadrature, while building the
     recurrence tree. If treemode!=0, store x-values and weights aswell.
     At first call, a and b should be 0 and 1 if isdefinite. */
  double mid;
  /* Allocate current node: */
  node = std::make_unique<qss_node>();
  if (treemode != 0) {
    node->x.resize(15);
    node->w.resize(15);
  }

  gk_quad((*test), (*function), params_for_function, node.get(), a, b, isindefinite);
  if ((fabs(node->err / node->I) < tol) || (tol >= 1.0)) {
    /* Stop recursion and return. tol>=1.0 in case of I=0 infinite recursion */
    return;
  }
  else {
    /* Call gk_adapt recursively on children:*/
    mid = 0.5 * (a + b);
    //printf("<-%g,%g,%g,%g",mid,tol,(*node)->err,(*node)->I);
    gk_adapt(node->left,
             (*test),
             (*function),
             params_for_function,
             1.5 * tol,
             treemode,
             a,
             mid,
             isindefinite);
    //printf("%g->",mid);
    gk_adapt(node->right,
             (*test),
             (*function),
             params_for_function,
             1.5 * tol,
             treemode,
             mid,
             b,
             isindefinite);
    /* Update integral and error in this node and return: */
    /* Actually, it is more convenient just to keep the nodes own estimate of the
       integral for our purposes.
       (*node)->I = (*node)->left->I + (*node)->right->I;
       (*node)->err = sqrt(pow(node->left->err,2)+pow(node->right->err,2));
    */
  }
}

void compute_Laguerre(
    double* x, double* w, int N, double alpha, double* b, double* c, int totalweight) {
  int i, j, iter, maxiter = 10;
  double x0 = 0., r1, r2, ratio, d, logprod, logcc;
  double p0, p1 = 0.0, p2, dp0, dp1, dp2 = 0.0;
  double eps = 1e-14;
  /* Initialise recursion coefficients: */
  for (i = 0; i < N; i++) {
    b[i] = alpha + 2.0 * i + 1.0;
    c[i] = i * (alpha + i);
  }
  logprod = 0.0;
  for (i = 1; i < N; i++)
    logprod += log(c[i]);
  logcc = lgamma(alpha + 1) + logprod;

  /* Loop over roots: */
  for (i = 0; i < N; i++) {
    /* Estimate root: */
    if (i == 0) {
      x0 = (1.0 + alpha) * (3.0 + 0.92 * alpha) / (1.0 + 2.4 * N + 1.8 * alpha);
    }
    else if (i == 1) {
      x0 += (15.0 + 6.25 * alpha) / (1.0 + 0.9 * alpha + 2.5 * N);
    }
    else {
      r1     = (1.0 + 2.55 * (i - 1)) / (1.9 * (i - 1));
      r2     = 1.26 * (i - 1) * alpha / (1.0 + 3.5 * (i - 1));
      ratio  = (r1 + r2) / (1.0 + 0.3 * alpha);
      x0    += ratio * (x0 - x[i - 2]);
    }
    /* Refine root using Newtons method: */
    for (iter = 1; iter <= maxiter; iter++) {
      /* We need to find p2=L_N(x0), dp2=L'_N(x0) and
	 p1 = L_(N-1)(x0): */
      p1  = 1.0;
      dp1 = 0.0;
      p2  = x0 - alpha - 1.0;
      dp2 = 1.0;
      for (j = 1; j < N; j++) {
        p0  = p1;
        dp0 = dp1;
        p1  = p2;
        dp1 = dp2;
        p2  = (x0 - b[j]) * p1 - c[j] * p0;
        dp2 = (x0 - b[j]) * dp1 + p1 - c[j] * dp0;
      }
      /* New guess at root: */
      d   = p2 / dp2;
      x0 -= d;
      if (fabs(d) <= eps * (fabs(x0) + 1.0))
        break;
    }
    /* Okay, write root and weight: */
    x[i] = x0;

    if (totalweight)
      w[i] = exp(x0 + logcc - log(dp2 * p1));
    else
      w[i] = exp(logcc - log(dp2 * p1));
  }
}

bool compute_FermiDirac(double* x, double* w, int N) {
  if (N < 1 || N > kFermiDiracMaxOrder)
    return false;

  /* The reference monic Laguerre recurrence is
       L_(n+1) = (x - (2n+1)) L_n - n^2 L_(n-1).
     Modified Chebyshev obtains the recurrence of the FD orthogonal
     polynomials from nu_n = integral L_n(x)/(exp(x)+1) dx without ever
     forming ill-conditioned power moments. */
  std::vector<double> moments(2 * N);
  double factorial = 1.;
  for (int n = 0; n < 2 * N; ++n) {
    if (n > 0)
      factorial *= n;
    // fermi_dirac_modified_moment() uses the conventional Laguerre L_n
    // (leading coefficient (-1)^n/n!), while the recurrence below requires
    // the monic P_n=(-1)^n n! L_n.
    moments[n] = (n % 2 == 0 ? factorial : -factorial) * fermi_dirac_modified_moment(n);
  }

  if (!(moments[0] > 0.) || !std::isfinite(moments[0]))
    return false;

  std::vector<std::vector<double>> sigma(N, std::vector<double>(2 * N, 0.));
  std::vector<double> diagonal(N);
  std::vector<double> off_diagonal_squared(N, 0.);
  sigma[0]    = moments;
  diagonal[0] = 1. + moments[1] / moments[0];

  for (int n = 1; n < N; ++n) {
    for (int l = n; l <= 2 * N - n - 1; ++l) {
      const double laguerre_diagonal        = 2. * l + 1.;
      const double laguerre_subdiag_squared = static_cast<double>(l) * l;
      sigma[n][l] = sigma[n - 1][l + 1] - (diagonal[n - 1] - laguerre_diagonal) * sigma[n - 1][l] -
                    off_diagonal_squared[n - 1] * (n > 1 ? sigma[n - 2][l] : 0.) +
                    laguerre_subdiag_squared * sigma[n - 1][l - 1];
    }

    const double previous = sigma[n - 1][n - 1];
    const double current  = sigma[n][n];
    if (previous == 0. || current == 0. || !std::isfinite(previous) || !std::isfinite(current)) {
      return false;
    }

    diagonal[n] = (2. * n + 1.) - sigma[n - 1][n] / previous + sigma[n][n + 1] / current;
    off_diagonal_squared[n] = current / previous;
    if (!(off_diagonal_squared[n] > 0.) || !std::isfinite(diagonal[n]) ||
        !std::isfinite(off_diagonal_squared[n])) {
      return false;
    }
  }

  std::vector<double> subdiagonal(N, 0.);
  for (int n = 1; n < N; ++n)
    subdiagonal[n] = std::sqrt(off_diagonal_squared[n]);

  std::vector<double> eigenvectors(N * N, 0.);
  for (int n = 0; n < N; ++n)
    eigenvectors[n * N + n] = 1.;
  if (!diagonalize_symmetric_tridiagonal(diagonal, subdiagonal, eigenvectors))
    return false;

  for (int n = 0; n < N; ++n) {
    x[n] = diagonal[n];
    w[n] = moments[0] * eigenvectors[n] * eigenvectors[n];
    if (!(x[n] > 0.) || !(w[n] > 0.) || !std::isfinite(x[n]) || !std::isfinite(w[n]))
      return false;
  }
  return true;
}

void gk_quad(void (*test)(void* params_for_function, double q, double* psi),
             void (*function)(void* params_for_function, double q, double* f0),
             void* params_for_function,
             qss_node* node,
             double a,
             double b,
             int isindefinite) {
  const double z_k[15] = {-0.991455371120813,
                          -0.949107912342759,
                          -0.864864423359769,
                          -0.741531185599394,
                          -0.586087235467691,
                          -0.405845151377397,
                          -0.207784955007898,
                          0.0,
                          0.207784955007898,
                          0.405845151377397,
                          0.586087235467691,
                          0.741531185599394,
                          0.864864423359769,
                          0.949107912342759,
                          0.991455371120813};
  const double w_k[15] = {0.022935322010529,
                          0.063092092629979,
                          0.104790010322250,
                          0.140653259715525,
                          0.169004726639267,
                          0.190350578064785,
                          0.204432940075298,
                          0.209482141084728,
                          0.204432940075298,
                          0.190350578064785,
                          0.169004726639267,
                          0.140653259715525,
                          0.104790010322250,
                          0.063092092629979,
                          0.022935322010529};
  const double w_g[7]  = {0.129484966168870,
                          0.279705391489277,
                          0.381830050505119,
                          0.417959183673469,
                          0.381830050505119,
                          0.279705391489277,
                          0.129484966168870};
  int i, j;
  double x, wg, wk, t, Ik, Ig, y, y2;

  /* 	Loop through abscissas, transform the interval and form the Kronrod
     15 point estimate of the integral.
     Every second time we update the Gauss 7 point quadrature estimate. */

  Ik = 0.0;
  Ig = 0.0;
  for (i = 0; i < 15; i++) {
    /* Transform z into t in interval between a and b: */
    t = 0.5 * (a * (1 - z_k[i]) + b * (1 + z_k[i]));
    /* Modify weight such that it reflects the linear transformation above: */
    wk = 0.5 * (b - a) * w_k[i];
    if (isindefinite) {
      /* Transform t into x in interval between 0 and inf: */
      x = 1.0 / t - 1.0;
      /* Modify weight accordingly: */
      wk = wk / (t * t);
    }
    else {
      x = t;
    }
    (*test)(params_for_function, x, &y);
    (*function)(params_for_function, x, &y2);
    wk *= y2;
    /* Update Kronrod integral: */
    Ik += wk * y;
    /* If node->x and node->w is allocated, store values: */
    if (!node->x.empty())
      node->x[i] = x;
    if (!node->w.empty())
      node->w[i] = wk;
    /* If i is uneven, update Gauss integral: */
    if ((i % 2) == 1) {
      j = (i - 1) / 2;
      /* Transform weight according to linear transformation: */
      wg = 0.5 * (b - a) * w_g[j];
      if (isindefinite) {
        /* Transform weight according to non-linear transformation x = 1/t -1: */
        wg = wg / (t * t);
      }
      /* Update integral: */
      Ig += wg * y * y2;
    }
  }
  node->err = pow(200 * fabs(Ik - Ig), 1.5);
  node->I   = Ik;
}

/**
 * This routine computes the weights and abscissas of a Gauss-Legendre quadrature between -1 and 1
 *
 * @param mu     Input/output: Vector of cos(beta) values
 * @param w8     Input/output: Vector of quadrature weights
 * @param n      Input       : Number of quadrature points
 * @param tol    Input       : tolerance on each mu
 *
 * From Numerical recipes
 **/

void quadrature_gauss_legendre(double* mu, double* w8, int n, double tol) {
  int m, j, i, counter;
  double z1, z, pp, p3, p2, p1;

  m = (n + 1) / 2;
  for (i = 1; i <= m; i++) {
    z       = cos(_PI_ * ((double) i - 0.25) / ((double) n + 0.5));
    counter = 0;
    do {
      p1 = 1.0;
      p2 = 0.0;
      for (j = 1; j <= n; j++) {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
      }
      pp = n * (z * p1 - p2) / (z * z - 1.0);
      z1 = z;
      z  = z1 - p1 / pp;
      counter++;
      class_test(counter == _MAX_IT_,
                 "maximum number of iteration reached: increase either _MAX_IT_ or tol\n");
    } while (fabs(z - z1) > tol);
    mu[i - 1] = -z;
    mu[n - i] = z;
    w8[i - 1] = 2.0 / ((1.0 - z * z) * pp * pp);
    w8[n - i] = w8[i - 1];
  }
}
