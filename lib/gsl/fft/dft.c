#include <config.h.in>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/err/gsl_errno.h>
#include <gsl/complex/gsl_complex.h>

#include <gsl/gsl_dft_complex.h>
#include <gsl/gsl_dft_complex_float.h>

#include "complex_internal.h"

#define BASE_DOUBLE
#include "templates_on.h"
#include "dft_source.c"
#include "templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "templates_on.h"
#include "dft_source.c"
#include "templates_off.h"
#undef  BASE_FLOAT

