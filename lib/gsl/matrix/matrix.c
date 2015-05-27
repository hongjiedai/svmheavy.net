#include "stdafx.h"
#include <config.h.in>
#include <gsl/err/gsl_errno.h>
#include <gsl/matrix/gsl_matrix.h>

#define BASE_GSL_COMPLEX_LONG
#include "templates_on.h"
#ifndef HIDE_INLINE_STATIC
BASE
FUNCTION(gsl_matrix, get) (const TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	BASE zero = ZERO;

	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VAL("first index out of range", GSL_EINVAL, zero);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VAL("second index out of range", GSL_EINVAL, zero);
		}
	}
	return *(BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}

void
FUNCTION(gsl_matrix, set) (TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const BASE x)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("second index out of range", GSL_EINVAL);
		}
	}
	*(BASE *)(m->data + MULTIPLICITY * (i * m->tda + j)) = x;
}


BASE *
FUNCTION(gsl_matrix, ptr) (TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("second index out of range", GSL_EINVAL);
		}
	}
	return (BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}


const BASE *
FUNCTION(gsl_matrix, const_ptr) (const TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("second index out of range", GSL_EINVAL);
		}
	}
	return (const BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}
#endif

#include "templates_off.h"
#undef  BASE_GSL_COMPLEX_LONG

#define BASE_GSL_COMPLEX
#include "templates_on.h"
#ifndef HIDE_INLINE_STATIC
BASE
FUNCTION(gsl_matrix, get) (const TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	BASE zero = ZERO;

	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VAL("first index out of range", GSL_EINVAL, zero);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VAL("second index out of range", GSL_EINVAL, zero);
		}
	}
	return *(BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}

void
FUNCTION(gsl_matrix, set) (TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const BASE x)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("second index out of range", GSL_EINVAL);
		}
	}
	*(BASE *)(m->data + MULTIPLICITY * (i * m->tda + j)) = x;
}


BASE *
FUNCTION(gsl_matrix, ptr) (TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("second index out of range", GSL_EINVAL);
		}
	}
	return (BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}


const BASE *
FUNCTION(gsl_matrix, const_ptr) (const TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("second index out of range", GSL_EINVAL);
		}
	}
	return (const BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}
#endif

#include "templates_off.h"
#undef  BASE_GSL_COMPLEX

#define BASE_GSL_COMPLEX_FLOAT
#include "templates_on.h"
#ifndef HIDE_INLINE_STATIC
BASE
FUNCTION(gsl_matrix, get) (const TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	BASE zero = ZERO;

	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VAL("first index out of range", GSL_EINVAL, zero);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VAL("second index out of range", GSL_EINVAL, zero);
		}
	}
	return *(BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}

void
FUNCTION(gsl_matrix, set) (TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const BASE x)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("second index out of range", GSL_EINVAL);
		}
	}
	*(BASE *)(m->data + MULTIPLICITY * (i * m->tda + j)) = x;
}


BASE *
FUNCTION(gsl_matrix, ptr) (TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("second index out of range", GSL_EINVAL);
		}
	}
	return (BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}


const BASE *
FUNCTION(gsl_matrix, const_ptr) (const TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("second index out of range", GSL_EINVAL);
		}
	}
	return (const BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}
#endif

#include "templates_off.h"
#undef  BASE_GSL_COMPLEX_FLOAT

#define BASE_LONG_DOUBLE
#include "templates_on.h"
#ifndef HIDE_INLINE_STATIC
BASE
FUNCTION(gsl_matrix, get) (const TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	BASE zero = ZERO;

	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VAL("first index out of range", GSL_EINVAL, zero);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VAL("second index out of range", GSL_EINVAL, zero);
		}
	}
	return *(BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}

void
FUNCTION(gsl_matrix, set) (TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const BASE x)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("second index out of range", GSL_EINVAL);
		}
	}
	*(BASE *)(m->data + MULTIPLICITY * (i * m->tda + j)) = x;
}


BASE *
FUNCTION(gsl_matrix, ptr) (TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("second index out of range", GSL_EINVAL);
		}
	}
	return (BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}


const BASE *
FUNCTION(gsl_matrix, const_ptr) (const TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("second index out of range", GSL_EINVAL);
		}
	}
	return (const BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}
#endif

#include "templates_off.h"
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
#include "templates_on.h"
#ifndef HIDE_INLINE_STATIC
BASE
FUNCTION(gsl_matrix, get) (const TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	BASE zero = ZERO;

	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VAL("first index out of range", GSL_EINVAL, zero);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VAL("second index out of range", GSL_EINVAL, zero);
		}
	}
	return *(BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}

void
FUNCTION(gsl_matrix, set) (TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const BASE x)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("second index out of range", GSL_EINVAL);
		}
	}
	*(BASE *)(m->data + MULTIPLICITY * (i * m->tda + j)) = x;
}


BASE *
FUNCTION(gsl_matrix, ptr) (TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("second index out of range", GSL_EINVAL);
		}
	}
	return (BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}


const BASE *
FUNCTION(gsl_matrix, const_ptr) (const TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("second index out of range", GSL_EINVAL);
		}
	}
	return (const BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}
#endif

#include "templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "templates_on.h"
#ifndef HIDE_INLINE_STATIC
BASE
FUNCTION(gsl_matrix, get) (const TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	BASE zero = ZERO;

	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VAL("first index out of range", GSL_EINVAL, zero);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VAL("second index out of range", GSL_EINVAL, zero);
		}
	}
	return *(BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}

void
FUNCTION(gsl_matrix, set) (TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const BASE x)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("second index out of range", GSL_EINVAL);
		}
	}
	*(BASE *)(m->data + MULTIPLICITY * (i * m->tda + j)) = x;
}


BASE *
FUNCTION(gsl_matrix, ptr) (TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("second index out of range", GSL_EINVAL);
		}
	}
	return (BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}


const BASE *
FUNCTION(gsl_matrix, const_ptr) (const TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("second index out of range", GSL_EINVAL);
		}
	}
	return (const BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}
#endif

#include "templates_off.h"
#undef  BASE_FLOAT

#define BASE_ULONG
#include "templates_on.h"
#ifndef HIDE_INLINE_STATIC
BASE
FUNCTION(gsl_matrix, get) (const TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	BASE zero = ZERO;

	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VAL("first index out of range", GSL_EINVAL, zero);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VAL("second index out of range", GSL_EINVAL, zero);
		}
	}
	return *(BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}

void
FUNCTION(gsl_matrix, set) (TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const BASE x)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("second index out of range", GSL_EINVAL);
		}
	}
	*(BASE *)(m->data + MULTIPLICITY * (i * m->tda + j)) = x;
}


BASE *
FUNCTION(gsl_matrix, ptr) (TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("second index out of range", GSL_EINVAL);
		}
	}
	return (BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}


const BASE *
FUNCTION(gsl_matrix, const_ptr) (const TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("second index out of range", GSL_EINVAL);
		}
	}
	return (const BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}
#endif

#include "templates_off.h"
#undef  BASE_ULONG

#define BASE_LONG
#include "templates_on.h"
#ifndef HIDE_INLINE_STATIC
BASE
FUNCTION(gsl_matrix, get) (const TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	BASE zero = ZERO;

	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VAL("first index out of range", GSL_EINVAL, zero);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VAL("second index out of range", GSL_EINVAL, zero);
		}
	}
	return *(BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}

void
FUNCTION(gsl_matrix, set) (TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const BASE x)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("second index out of range", GSL_EINVAL);
		}
	}
	*(BASE *)(m->data + MULTIPLICITY * (i * m->tda + j)) = x;
}


BASE *
FUNCTION(gsl_matrix, ptr) (TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("second index out of range", GSL_EINVAL);
		}
	}
	return (BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}


const BASE *
FUNCTION(gsl_matrix, const_ptr) (const TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("second index out of range", GSL_EINVAL);
		}
	}
	return (const BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}
#endif

#include "templates_off.h"
#undef  BASE_LONG

#define BASE_UINT
#include "templates_on.h"
#ifndef HIDE_INLINE_STATIC
BASE
FUNCTION(gsl_matrix, get) (const TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	BASE zero = ZERO;

	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VAL("first index out of range", GSL_EINVAL, zero);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VAL("second index out of range", GSL_EINVAL, zero);
		}
	}
	return *(BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}

void
FUNCTION(gsl_matrix, set) (TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const BASE x)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("second index out of range", GSL_EINVAL);
		}
	}
	*(BASE *)(m->data + MULTIPLICITY * (i * m->tda + j)) = x;
}


BASE *
FUNCTION(gsl_matrix, ptr) (TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("second index out of range", GSL_EINVAL);
		}
	}
	return (BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}


const BASE *
FUNCTION(gsl_matrix, const_ptr) (const TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("second index out of range", GSL_EINVAL);
		}
	}
	return (const BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}
#endif

#include "templates_off.h"
#undef  BASE_UINT

#define BASE_INT
#include "templates_on.h"
#ifndef HIDE_INLINE_STATIC
BASE
FUNCTION(gsl_matrix, get) (const TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	BASE zero = ZERO;

	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VAL("first index out of range", GSL_EINVAL, zero);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VAL("second index out of range", GSL_EINVAL, zero);
		}
	}
	return *(BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}

void
FUNCTION(gsl_matrix, set) (TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const BASE x)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("second index out of range", GSL_EINVAL);
		}
	}
	*(BASE *)(m->data + MULTIPLICITY * (i * m->tda + j)) = x;
}


BASE *
FUNCTION(gsl_matrix, ptr) (TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("second index out of range", GSL_EINVAL);
		}
	}
	return (BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}


const BASE *
FUNCTION(gsl_matrix, const_ptr) (const TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("second index out of range", GSL_EINVAL);
		}
	}
	return (const BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}
#endif

#include "templates_off.h"
#undef  BASE_INT

#define BASE_USHORT
#include "templates_on.h"
#ifndef HIDE_INLINE_STATIC
BASE
FUNCTION(gsl_matrix, get) (const TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	BASE zero = ZERO;

	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VAL("first index out of range", GSL_EINVAL, zero);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VAL("second index out of range", GSL_EINVAL, zero);
		}
	}
	return *(BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}

void
FUNCTION(gsl_matrix, set) (TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const BASE x)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("second index out of range", GSL_EINVAL);
		}
	}
	*(BASE *)(m->data + MULTIPLICITY * (i * m->tda + j)) = x;
}


BASE *
FUNCTION(gsl_matrix, ptr) (TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("second index out of range", GSL_EINVAL);
		}
	}
	return (BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}


const BASE *
FUNCTION(gsl_matrix, const_ptr) (const TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("second index out of range", GSL_EINVAL);
		}
	}
	return (const BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}
#endif

#include "templates_off.h"
#undef  BASE_USHORT

#define BASE_SHORT
#include "templates_on.h"
#ifndef HIDE_INLINE_STATIC
BASE
FUNCTION(gsl_matrix, get) (const TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	BASE zero = ZERO;

	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VAL("first index out of range", GSL_EINVAL, zero);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VAL("second index out of range", GSL_EINVAL, zero);
		}
	}
	return *(BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}

void
FUNCTION(gsl_matrix, set) (TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const BASE x)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("second index out of range", GSL_EINVAL);
		}
	}
	*(BASE *)(m->data + MULTIPLICITY * (i * m->tda + j)) = x;
}


BASE *
FUNCTION(gsl_matrix, ptr) (TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("second index out of range", GSL_EINVAL);
		}
	}
	return (BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}


const BASE *
FUNCTION(gsl_matrix, const_ptr) (const TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("second index out of range", GSL_EINVAL);
		}
	}
	return (const BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}
#endif

#include "templates_off.h"
#undef  BASE_SHORT

#define BASE_UCHAR
#include "templates_on.h"
#ifndef HIDE_INLINE_STATIC
BASE
FUNCTION(gsl_matrix, get) (const TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	BASE zero = ZERO;

	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VAL("first index out of range", GSL_EINVAL, zero);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VAL("second index out of range", GSL_EINVAL, zero);
		}
	}
	return *(BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}

void
FUNCTION(gsl_matrix, set) (TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const BASE x)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("second index out of range", GSL_EINVAL);
		}
	}
	*(BASE *)(m->data + MULTIPLICITY * (i * m->tda + j)) = x;
}


BASE *
FUNCTION(gsl_matrix, ptr) (TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("second index out of range", GSL_EINVAL);
		}
	}
	return (BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}


const BASE *
FUNCTION(gsl_matrix, const_ptr) (const TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("second index out of range", GSL_EINVAL);
		}
	}
	return (const BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}
#endif

#include "templates_off.h"
#undef  BASE_UCHAR

#define BASE_CHAR
#include "templates_on.h"
#ifndef HIDE_INLINE_STATIC
BASE
FUNCTION(gsl_matrix, get) (const TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	BASE zero = ZERO;

	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VAL("first index out of range", GSL_EINVAL, zero);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VAL("second index out of range", GSL_EINVAL, zero);
		}
	}
	return *(BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}

void
FUNCTION(gsl_matrix, set) (TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const BASE x)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("second index out of range", GSL_EINVAL);
		}
	}
	*(BASE *)(m->data + MULTIPLICITY * (i * m->tda + j)) = x;
}


BASE *
FUNCTION(gsl_matrix, ptr) (TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("second index out of range", GSL_EINVAL);
		}
	}
	return (BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}


const BASE *
FUNCTION(gsl_matrix, const_ptr) (const TYPE(gsl_matrix) * m,
const size_t i, const size_t j)
{
	if (gsl_check_range)
	{
		if (i >= m->size1)        /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("first index out of range", GSL_EINVAL);
		}
		else if (j >= m->size2)   /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("second index out of range", GSL_EINVAL);
		}
	}
	return (const BASE *)(m->data + MULTIPLICITY * (i * m->tda + j));
}
#endif

#include "templates_off.h"
#undef  BASE_CHAR
