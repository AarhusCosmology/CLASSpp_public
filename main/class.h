#ifndef __CLASS__
#define __CLASS__

/* standard libraries */
#include "float.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

/* tools for class */
#include "arrays.h"
#include "dei_rkck.h"
#include "growTable.h"
#include "parser.h"
#include "quadrature.h"

/* class modules */
#include "background.h"
#include "common.h"
#include "input_module.h"
#include "lensing.h"
#include "nonlinear.h"
#include "output.h"
#include "perturbations.h"
#include "primordial.h"
#include "spectra.h"
#include "thermodynamics.h"
#include "transfer.h"

class ClassConstants {
 public:
  static constexpr int sMAXTITLESTRINGLENGTH = _MAXTITLESTRINGLENGTH_;
  static constexpr int sFALSE                = _FALSE_;
  static constexpr int sARGUMENT_LENGTH_MAX  = _ARGUMENT_LENGTH_MAX_;
  static constexpr int sFAILURE              = _FAILURE_;
};

#endif
