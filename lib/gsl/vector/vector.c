/* vector/vector.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */
#include "stdafx.h"
#include "../config.h.in"
#include <gsl/err/gsl_errno.h>
#include <gsl/vector/gsl_vector.h>

/* turn off range checking at runtime if zero */
int gsl_check_range = 1;

#define BASE_GSL_COMPLEX_LONG
#include "../templates_on.h"
#ifndef HIDE_INLINE_STATIC
BASE
FUNCTION(gsl_vector, get) (const TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			BASE zero = ZERO;
			GSL_ERROR_VAL("index out of range", GSL_EINVAL, zero);
		}
	}

	/* The following line is a generalization of return v->data[i] */

	return *(BASE *)(v->data + MULTIPLICITY * i * v->stride);
}

void
FUNCTION(gsl_vector, set) (TYPE(gsl_vector) * v, const size_t i, BASE x)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("index out of range", GSL_EINVAL);
		}
	}

	/* The following line is a generalization of v->data[i] = x */

	*(BASE *)(v->data + MULTIPLICITY * i * v->stride) = x;
}

BASE *
FUNCTION(gsl_vector, ptr) (TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("index out of range", GSL_EINVAL);
		}
	}

	return (BASE *)(v->data + MULTIPLICITY * i * v->stride);
}

const BASE *
FUNCTION(gsl_vector, const_ptr) (const TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("index out of range", GSL_EINVAL);
		}
	}

	/* The following line is a generalization of return v->data[i] */

	return (const BASE *)(v->data + MULTIPLICITY * i * v->stride);
}
#endif

#include "../templates_off.h"
#undef  BASE_GSL_COMPLEX_LONG

#define BASE_GSL_COMPLEX
#include "../templates_on.h"
#ifndef HIDE_INLINE_STATIC
BASE
FUNCTION(gsl_vector, get) (const TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			BASE zero = ZERO;
			GSL_ERROR_VAL("index out of range", GSL_EINVAL, zero);
		}
	}

	/* The following line is a generalization of return v->data[i] */

	return *(BASE *)(v->data + MULTIPLICITY * i * v->stride);
}

void
FUNCTION(gsl_vector, set) (TYPE(gsl_vector) * v, const size_t i, BASE x)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("index out of range", GSL_EINVAL);
		}
	}

	/* The following line is a generalization of v->data[i] = x */

	*(BASE *)(v->data + MULTIPLICITY * i * v->stride) = x;
}

BASE *
FUNCTION(gsl_vector, ptr) (TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("index out of range", GSL_EINVAL);
		}
	}

	return (BASE *)(v->data + MULTIPLICITY * i * v->stride);
}

const BASE *
FUNCTION(gsl_vector, const_ptr) (const TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("index out of range", GSL_EINVAL);
		}
	}

	/* The following line is a generalization of return v->data[i] */

	return (const BASE *)(v->data + MULTIPLICITY * i * v->stride);
}
#endif

#include "../templates_off.h"
#undef  BASE_GSL_COMPLEX

#define BASE_GSL_COMPLEX_FLOAT
#include "../templates_on.h"
#ifndef HIDE_INLINE_STATIC
BASE
FUNCTION(gsl_vector, get) (const TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			BASE zero = ZERO;
			GSL_ERROR_VAL("index out of range", GSL_EINVAL, zero);
		}
	}

	/* The following line is a generalization of return v->data[i] */

	return *(BASE *)(v->data + MULTIPLICITY * i * v->stride);
}

void
FUNCTION(gsl_vector, set) (TYPE(gsl_vector) * v, const size_t i, BASE x)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("index out of range", GSL_EINVAL);
		}
	}

	/* The following line is a generalization of v->data[i] = x */

	*(BASE *)(v->data + MULTIPLICITY * i * v->stride) = x;
}

BASE *
FUNCTION(gsl_vector, ptr) (TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("index out of range", GSL_EINVAL);
		}
	}

	return (BASE *)(v->data + MULTIPLICITY * i * v->stride);
}

const BASE *
FUNCTION(gsl_vector, const_ptr) (const TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("index out of range", GSL_EINVAL);
		}
	}

	/* The following line is a generalization of return v->data[i] */

	return (const BASE *)(v->data + MULTIPLICITY * i * v->stride);
}
#endif

#include "../templates_off.h"
#undef  BASE_GSL_COMPLEX_FLOAT

#define BASE_LONG_DOUBLE
#include "../templates_on.h"
#ifndef HIDE_INLINE_STATIC
BASE
FUNCTION(gsl_vector, get) (const TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			BASE zero = ZERO;
			GSL_ERROR_VAL("index out of range", GSL_EINVAL, zero);
		}
	}

	/* The following line is a generalization of return v->data[i] */

	return *(BASE *)(v->data + MULTIPLICITY * i * v->stride);
}

void
FUNCTION(gsl_vector, set) (TYPE(gsl_vector) * v, const size_t i, BASE x)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("index out of range", GSL_EINVAL);
		}
	}

	/* The following line is a generalization of v->data[i] = x */

	*(BASE *)(v->data + MULTIPLICITY * i * v->stride) = x;
}

BASE *
FUNCTION(gsl_vector, ptr) (TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("index out of range", GSL_EINVAL);
		}
	}

	return (BASE *)(v->data + MULTIPLICITY * i * v->stride);
}

const BASE *
FUNCTION(gsl_vector, const_ptr) (const TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("index out of range", GSL_EINVAL);
		}
	}

	/* The following line is a generalization of return v->data[i] */

	return (const BASE *)(v->data + MULTIPLICITY * i * v->stride);
}
#endif

#include "../templates_off.h"
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
#include "../templates_on.h"
#ifndef HIDE_INLINE_STATIC
BASE
FUNCTION(gsl_vector, get) (const TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			BASE zero = ZERO;
			GSL_ERROR_VAL("index out of range", GSL_EINVAL, zero);
		}
	}

	/* The following line is a generalization of return v->data[i] */

	return *(BASE *)(v->data + MULTIPLICITY * i * v->stride);
}

void
FUNCTION(gsl_vector, set) (TYPE(gsl_vector) * v, const size_t i, BASE x)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("index out of range", GSL_EINVAL);
		}
	}

	/* The following line is a generalization of v->data[i] = x */

	*(BASE *)(v->data + MULTIPLICITY * i * v->stride) = x;
}

BASE *
FUNCTION(gsl_vector, ptr) (TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("index out of range", GSL_EINVAL);
		}
	}

	return (BASE *)(v->data + MULTIPLICITY * i * v->stride);
}

const BASE *
FUNCTION(gsl_vector, const_ptr) (const TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("index out of range", GSL_EINVAL);
		}
	}

	/* The following line is a generalization of return v->data[i] */

	return (const BASE *)(v->data + MULTIPLICITY * i * v->stride);
}
#endif

#include "../templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "../templates_on.h"
#ifndef HIDE_INLINE_STATIC
BASE
FUNCTION(gsl_vector, get) (const TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			BASE zero = ZERO;
			GSL_ERROR_VAL("index out of range", GSL_EINVAL, zero);
		}
	}

	/* The following line is a generalization of return v->data[i] */

	return *(BASE *)(v->data + MULTIPLICITY * i * v->stride);
}

void
FUNCTION(gsl_vector, set) (TYPE(gsl_vector) * v, const size_t i, BASE x)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("index out of range", GSL_EINVAL);
		}
	}

	/* The following line is a generalization of v->data[i] = x */

	*(BASE *)(v->data + MULTIPLICITY * i * v->stride) = x;
}

BASE *
FUNCTION(gsl_vector, ptr) (TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("index out of range", GSL_EINVAL);
		}
	}

	return (BASE *)(v->data + MULTIPLICITY * i * v->stride);
}

const BASE *
FUNCTION(gsl_vector, const_ptr) (const TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("index out of range", GSL_EINVAL);
		}
	}

	/* The following line is a generalization of return v->data[i] */

	return (const BASE *)(v->data + MULTIPLICITY * i * v->stride);
}
#endif

#include "../templates_off.h"
#undef  BASE_FLOAT

#define BASE_ULONG
#include "../templates_on.h"
#ifndef HIDE_INLINE_STATIC
BASE
FUNCTION(gsl_vector, get) (const TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			BASE zero = ZERO;
			GSL_ERROR_VAL("index out of range", GSL_EINVAL, zero);
		}
	}

	/* The following line is a generalization of return v->data[i] */

	return *(BASE *)(v->data + MULTIPLICITY * i * v->stride);
}

void
FUNCTION(gsl_vector, set) (TYPE(gsl_vector) * v, const size_t i, BASE x)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("index out of range", GSL_EINVAL);
		}
	}

	/* The following line is a generalization of v->data[i] = x */

	*(BASE *)(v->data + MULTIPLICITY * i * v->stride) = x;
}

BASE *
FUNCTION(gsl_vector, ptr) (TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("index out of range", GSL_EINVAL);
		}
	}

	return (BASE *)(v->data + MULTIPLICITY * i * v->stride);
}

const BASE *
FUNCTION(gsl_vector, const_ptr) (const TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("index out of range", GSL_EINVAL);
		}
	}

	/* The following line is a generalization of return v->data[i] */

	return (const BASE *)(v->data + MULTIPLICITY * i * v->stride);
}
#endif

#include "../templates_off.h"
#undef  BASE_ULONG

#define BASE_LONG
#include "../templates_on.h"
#ifndef HIDE_INLINE_STATIC
BASE
FUNCTION(gsl_vector, get) (const TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			BASE zero = ZERO;
			GSL_ERROR_VAL("index out of range", GSL_EINVAL, zero);
		}
	}

	/* The following line is a generalization of return v->data[i] */

	return *(BASE *)(v->data + MULTIPLICITY * i * v->stride);
}

void
FUNCTION(gsl_vector, set) (TYPE(gsl_vector) * v, const size_t i, BASE x)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("index out of range", GSL_EINVAL);
		}
	}

	/* The following line is a generalization of v->data[i] = x */

	*(BASE *)(v->data + MULTIPLICITY * i * v->stride) = x;
}

BASE *
FUNCTION(gsl_vector, ptr) (TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("index out of range", GSL_EINVAL);
		}
	}

	return (BASE *)(v->data + MULTIPLICITY * i * v->stride);
}

const BASE *
FUNCTION(gsl_vector, const_ptr) (const TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("index out of range", GSL_EINVAL);
		}
	}

	/* The following line is a generalization of return v->data[i] */

	return (const BASE *)(v->data + MULTIPLICITY * i * v->stride);
}
#endif

#include "../templates_off.h"
#undef  BASE_LONG

#define BASE_UINT
#include "../templates_on.h"
#ifndef HIDE_INLINE_STATIC
BASE
FUNCTION(gsl_vector, get) (const TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			BASE zero = ZERO;
			GSL_ERROR_VAL("index out of range", GSL_EINVAL, zero);
		}
	}

	/* The following line is a generalization of return v->data[i] */

	return *(BASE *)(v->data + MULTIPLICITY * i * v->stride);
}

void
FUNCTION(gsl_vector, set) (TYPE(gsl_vector) * v, const size_t i, BASE x)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("index out of range", GSL_EINVAL);
		}
	}

	/* The following line is a generalization of v->data[i] = x */

	*(BASE *)(v->data + MULTIPLICITY * i * v->stride) = x;
}

BASE *
FUNCTION(gsl_vector, ptr) (TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("index out of range", GSL_EINVAL);
		}
	}

	return (BASE *)(v->data + MULTIPLICITY * i * v->stride);
}

const BASE *
FUNCTION(gsl_vector, const_ptr) (const TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("index out of range", GSL_EINVAL);
		}
	}

	/* The following line is a generalization of return v->data[i] */

	return (const BASE *)(v->data + MULTIPLICITY * i * v->stride);
}
#endif

#include "../templates_off.h"
#undef  BASE_UINT

#define BASE_INT
#include "../templates_on.h"
#ifndef HIDE_INLINE_STATIC
BASE
FUNCTION(gsl_vector, get) (const TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			BASE zero = ZERO;
			GSL_ERROR_VAL("index out of range", GSL_EINVAL, zero);
		}
	}

	/* The following line is a generalization of return v->data[i] */

	return *(BASE *)(v->data + MULTIPLICITY * i * v->stride);
}

void
FUNCTION(gsl_vector, set) (TYPE(gsl_vector) * v, const size_t i, BASE x)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("index out of range", GSL_EINVAL);
		}
	}

	/* The following line is a generalization of v->data[i] = x */

	*(BASE *)(v->data + MULTIPLICITY * i * v->stride) = x;
}

BASE *
FUNCTION(gsl_vector, ptr) (TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("index out of range", GSL_EINVAL);
		}
	}

	return (BASE *)(v->data + MULTIPLICITY * i * v->stride);
}

const BASE *
FUNCTION(gsl_vector, const_ptr) (const TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("index out of range", GSL_EINVAL);
		}
	}

	/* The following line is a generalization of return v->data[i] */

	return (const BASE *)(v->data + MULTIPLICITY * i * v->stride);
}
#endif

#include "../templates_off.h"
#undef  BASE_INT

#define BASE_USHORT
#include "../templates_on.h"
#ifndef HIDE_INLINE_STATIC
BASE
FUNCTION(gsl_vector, get) (const TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			BASE zero = ZERO;
			GSL_ERROR_VAL("index out of range", GSL_EINVAL, zero);
		}
	}

	/* The following line is a generalization of return v->data[i] */

	return *(BASE *)(v->data + MULTIPLICITY * i * v->stride);
}

void
FUNCTION(gsl_vector, set) (TYPE(gsl_vector) * v, const size_t i, BASE x)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("index out of range", GSL_EINVAL);
		}
	}

	/* The following line is a generalization of v->data[i] = x */

	*(BASE *)(v->data + MULTIPLICITY * i * v->stride) = x;
}

BASE *
FUNCTION(gsl_vector, ptr) (TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("index out of range", GSL_EINVAL);
		}
	}

	return (BASE *)(v->data + MULTIPLICITY * i * v->stride);
}

const BASE *
FUNCTION(gsl_vector, const_ptr) (const TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("index out of range", GSL_EINVAL);
		}
	}

	/* The following line is a generalization of return v->data[i] */

	return (const BASE *)(v->data + MULTIPLICITY * i * v->stride);
}
#endif

#include "../templates_off.h"
#undef  BASE_USHORT

#define BASE_SHORT
#include "../templates_on.h"
#ifndef HIDE_INLINE_STATIC
BASE
FUNCTION(gsl_vector, get) (const TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			BASE zero = ZERO;
			GSL_ERROR_VAL("index out of range", GSL_EINVAL, zero);
		}
	}

	/* The following line is a generalization of return v->data[i] */

	return *(BASE *)(v->data + MULTIPLICITY * i * v->stride);
}

void
FUNCTION(gsl_vector, set) (TYPE(gsl_vector) * v, const size_t i, BASE x)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("index out of range", GSL_EINVAL);
		}
	}

	/* The following line is a generalization of v->data[i] = x */

	*(BASE *)(v->data + MULTIPLICITY * i * v->stride) = x;
}

BASE *
FUNCTION(gsl_vector, ptr) (TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("index out of range", GSL_EINVAL);
		}
	}

	return (BASE *)(v->data + MULTIPLICITY * i * v->stride);
}

const BASE *
FUNCTION(gsl_vector, const_ptr) (const TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("index out of range", GSL_EINVAL);
		}
	}

	/* The following line is a generalization of return v->data[i] */

	return (const BASE *)(v->data + MULTIPLICITY * i * v->stride);
}
#endif

#include "../templates_off.h"
#undef  BASE_SHORT

#define BASE_UCHAR
#include "../templates_on.h"
#ifndef HIDE_INLINE_STATIC
BASE
FUNCTION(gsl_vector, get) (const TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			BASE zero = ZERO;
			GSL_ERROR_VAL("index out of range", GSL_EINVAL, zero);
		}
	}

	/* The following line is a generalization of return v->data[i] */

	return *(BASE *)(v->data + MULTIPLICITY * i * v->stride);
}

void
FUNCTION(gsl_vector, set) (TYPE(gsl_vector) * v, const size_t i, BASE x)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("index out of range", GSL_EINVAL);
		}
	}

	/* The following line is a generalization of v->data[i] = x */

	*(BASE *)(v->data + MULTIPLICITY * i * v->stride) = x;
}

BASE *
FUNCTION(gsl_vector, ptr) (TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("index out of range", GSL_EINVAL);
		}
	}

	return (BASE *)(v->data + MULTIPLICITY * i * v->stride);
}

const BASE *
FUNCTION(gsl_vector, const_ptr) (const TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("index out of range", GSL_EINVAL);
		}
	}

	/* The following line is a generalization of return v->data[i] */

	return (const BASE *)(v->data + MULTIPLICITY * i * v->stride);
}
#endif

#include "../templates_off.h"
#undef  BASE_UCHAR

#define BASE_CHAR
#include "../templates_on.h"
#ifndef HIDE_INLINE_STATIC
BASE
FUNCTION(gsl_vector, get) (const TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			BASE zero = ZERO;
			GSL_ERROR_VAL("index out of range", GSL_EINVAL, zero);
		}
	}

	/* The following line is a generalization of return v->data[i] */

	return *(BASE *)(v->data + MULTIPLICITY * i * v->stride);
}

void
FUNCTION(gsl_vector, set) (TYPE(gsl_vector) * v, const size_t i, BASE x)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VOID("index out of range", GSL_EINVAL);
		}
	}

	/* The following line is a generalization of v->data[i] = x */

	*(BASE *)(v->data + MULTIPLICITY * i * v->stride) = x;
}

BASE *
FUNCTION(gsl_vector, ptr) (TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("index out of range", GSL_EINVAL);
		}
	}

	return (BASE *)(v->data + MULTIPLICITY * i * v->stride);
}

const BASE *
FUNCTION(gsl_vector, const_ptr) (const TYPE(gsl_vector) * v, const size_t i)
{
	if (gsl_check_range)
	{
		if (i >= v->size)         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_NULL("index out of range", GSL_EINVAL);
		}
	}

	/* The following line is a generalization of return v->data[i] */

	return (const BASE *)(v->data + MULTIPLICITY * i * v->stride);
}
#endif
#include "../templates_off.h"
#undef  BASE_CHAR
