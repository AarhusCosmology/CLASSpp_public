/** @file common.h Generic libraries, parameters and functions used in the whole code. */

#include "float.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "svnversion.h"

#ifdef __cplusplus
#include <algorithm>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "exceptions.h"

class InputModule;
class BackgroundModule;
class ThermodynamicsModule;
class PerturbationsModule;
class PrimordialModule;
class NonlinearModule;
class TransferModule;
class SpectraModule;
class LensingModule;
class FileContent;

typedef std::shared_ptr<const InputModule> InputModulePtr;
typedef std::shared_ptr<const BackgroundModule> BackgroundModulePtr;
typedef std::shared_ptr<const ThermodynamicsModule> ThermodynamicsModulePtr;
typedef std::shared_ptr<const PerturbationsModule> PerturbationsModulePtr;
typedef std::shared_ptr<const PrimordialModule> PrimordialModulePtr;
typedef std::shared_ptr<const NonlinearModule> NonlinearModulePtr;
typedef std::shared_ptr<const TransferModule> TransferModulePtr;
typedef std::shared_ptr<const SpectraModule> SpectraModulePtr;
typedef std::shared_ptr<const LensingModule> LensingModulePtr;
#endif

#ifndef __COMMON__
#define __COMMON__

#include "constants.h"
#include "errors.h"
#ifdef __cplusplus
#include "precision.h"
#endif

#define _VERSION_ "v2.9.0"
/* @cond INCLUDE_WITH_DOXYGEN */

#define _MAX_IT_ \
  10000 /**< default maximum number of iterations in conditional loops (to avoid infinite loops) */

#define _QUADRATURE_MAX_ \
  250 /**< maximum allowed number of abssices in quadrature integral estimation */

#define _QUADRATURE_MAX_BG_ \
  800 /**< maximum allowed number of abssices in quadrature integral estimation */

#define _TOLVAR_ \
  100. /**< The minimum allowed variation is the machine precision times this number */

#define _HUGE_ 1.e99

#define _EPSILON_ 1.e-10

#define _OUTPUTPRECISION_ 12 /**< Number of significant digits in some output files */

#define _COLUMNWIDTH_ \
  24 /**< Must be at least _OUTPUTPRECISION_+8 for guaranteed fixed width columns */

#define _DELIMITER_ "\t" /**< character used for delimiting titles in the title strings */

#ifndef __CLASSDIR__
#define __CLASSDIR__ \
  "." /**< The directory of CLASS. This is set to the absolute path to the CLASS directory so this is just a failsafe. */
#endif

#define index_symmetric_matrix(i1, i2, N)            \
  (((i1) <= (i2))                                    \
       ? ((i2) + N * (i1) - ((i1) * ((i1) + 1)) / 2) \
       : ((i1) + N * (i2) -                          \
          ((i2) * ((i2) + 1)) /                      \
              2)) /**< assigns an index from 0 to [N(N+1)/2-1] to the coefficients M_{i1,i2} of an N*N symmetric matrix; useful for converting a symmetric matrix to a vector, without losing or double-counting any information */
/* @endcond */
// needed because of weird openmp bug on macosx lion...

/* macro for defining indices (usually one, sometimes a block) */
#define class_define_index(index, condition, running_index, number_of_indices) \
  {                                                                            \
    if (condition) {                                                           \
      index          = running_index;                                          \
      running_index += number_of_indices;                                      \
    }                                                                          \
  }

/* macros for writing formatted output */
#define class_fprintf_double(file, output, condition)                    \
  {                                                                      \
    if (condition == _TRUE_)                                             \
      fprintf(file, "%*.*e ", _COLUMNWIDTH_, _OUTPUTPRECISION_, output); \
  }

#define class_fprintf_double_or_default(file, output, condition, defaultvalue) \
  {                                                                            \
    if (condition == _TRUE_)                                                   \
      fprintf(file, "%*.*e ", _COLUMNWIDTH_, _OUTPUTPRECISION_, output);       \
    else                                                                       \
      fprintf(file, "%*.*e ", _COLUMNWIDTH_, _OUTPUTPRECISION_, defaultvalue); \
  }

#define class_fprintf_int(file, output, condition)                \
  {                                                               \
    if (condition == _TRUE_)                                      \
      fprintf(file,                                               \
              "%*d%*s ",                                          \
              std::max(0, _COLUMNWIDTH_ - _OUTPUTPRECISION_ - 5), \
              output,                                             \
              _OUTPUTPRECISION_ + 5,                              \
              " ");                                               \
  }

#define class_fprintf_columntitle(file, title, condition, colnum)            \
  {                                                                          \
    if (condition == _TRUE_)                                                 \
      fprintf(file,                                                          \
              "%*s%2d:%-*s ",                                                \
              std::max(0,                                                    \
                       std::min(_COLUMNWIDTH_ - _OUTPUTPRECISION_ - 6 - 3,   \
                                _COLUMNWIDTH_ - ((int) strlen(title)) - 3)), \
              "",                                                            \
              colnum++,                                                      \
              _OUTPUTPRECISION_ + 6,                                         \
              title);                                                        \
  }

#define class_store_columntitle(titlestring, title, condition) \
  {                                                            \
    if (condition == _TRUE_) {                                 \
      titlestring += title;                                    \
      titlestring += _DELIMITER_;                              \
    }                                                          \
  }

#define class_store_double(storage, value, condition, dataindex) \
  {                                                              \
    if (condition == _TRUE_)                                     \
      storage[dataindex++] = value;                              \
  }

#define class_store_double_or_default(storage, value, condition, dataindex, defaultvalue) \
  {                                                                                       \
    if (condition == _TRUE_)                                                              \
      storage[dataindex++] = value;                                                       \
    else                                                                                  \
      storage[dataindex++] = defaultvalue;                                                \
  }

#ifdef __cplusplus
/* Not in the extern "C" block: takes a std::string& (a C++ type), so it cannot
 * have C language linkage; guarded by __cplusplus so C translation units never
 * see this declaration. */
int get_number_of_titles(const std::string& titlestring);
#endif

#endif
