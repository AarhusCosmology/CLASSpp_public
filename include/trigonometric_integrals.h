/**
 * definitions for module trigonometric_integrals.c
 */

#ifndef __TRIGONOMETRIC_INTEGRALS__
#define __TRIGONOMETRIC_INTEGRALS__

#include "common.h"

/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

int cosine_integral(double x, double* Ci);

int sine_integral(double x, double* Si);

#ifdef __cplusplus
}
#endif

#endif
