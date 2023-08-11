#ifndef __CLASS__
#define __CLASS__

/* standard libraries */
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#include "float.h"

/* tools for class */
#include "quadrature.h"
#include "growTable.h"
#include "arrays.h"
#include "dei_rkck.h"
#include "parser.h"

/* class modules */
#include "common.h"
#include "input_module.h"
#include "background.h"
#include "thermodynamics.h"
#include "perturbations.h"
#include "primordial.h"
#include "nonlinear.h"
#include "transfer.h"
#include "spectra.h"
#include "lensing.h"
#include "output.h"

class ClassConstants {
public:
  static constexpr int sMAXTITLESTRINGLENGTH = _MAXTITLESTRINGLENGTH_;
  static constexpr int sFALSE = _FALSE_;
  static constexpr int sARGUMENT_LENGTH_MAX = _ARGUMENT_LENGTH_MAX_;
  static constexpr int sFAILURE = _FAILURE_;
};

#endif
