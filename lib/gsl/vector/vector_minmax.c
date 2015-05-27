#include "stdafx.h"
//#include "config.h.in"
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/vector/gsl_vector.h>

#define BASE_LONG_DOUBLE
#include "../templates_on.h"
BASE
FUNCTION(gsl_vector, max) (const TYPE(gsl_vector) * v)
{
	/* finds the largest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE max = v->data[0 * stride];
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x > max)
			max = x;
#ifdef FP
		if (isnan(x))
			return x;
#endif
	}

	return max;
}

BASE
FUNCTION(gsl_vector, min) (const TYPE(gsl_vector) * v)
{
	/* finds the smallest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE min = v->data[0 * stride];
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
			min = x;
#ifdef FP
		if (isnan(x))
			return x;
#endif
	}

	return min;
}

void
FUNCTION(gsl_vector, minmax) (const TYPE(gsl_vector) * v,
BASE * min_out,
BASE * max_out)
{
	/* finds the smallest and largest elements of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE max = v->data[0 * stride];
	BASE min = v->data[0 * stride];

	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
		{
			min = x;
		}
		if (x > max)
		{
			max = x;
		}
#ifdef FP
		if (isnan(x))
		{
			min = x;
			max = x;
			break;
		}
#endif
	}

	*min_out = min;
	*max_out = max;
}


size_t
FUNCTION(gsl_vector, max_index) (const TYPE(gsl_vector) * v)
{
	/* finds the largest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE max = v->data[0 * stride];
	size_t imax = 0;
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x > max)
		{
			max = x;
			imax = i;
		}
#ifdef FP
		if (isnan(x))
		{
			return i;
		}
#endif
	}

	return imax;
}

size_t
FUNCTION(gsl_vector, min_index) (const TYPE(gsl_vector) * v)
{
	/* finds the smallest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE min = v->data[0 * stride];
	size_t imin = 0;
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
		{
			min = x;
			imin = i;
		}
#ifdef FP
		if (isnan(x))
		{
			return i;
		}
#endif
	}

	return imin;
}


void
FUNCTION(gsl_vector, minmax_index) (const TYPE(gsl_vector) * v,
size_t * imin_out,
size_t * imax_out)
{
	/* finds the smallest and largest elements of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	size_t imin = 0, imax = 0;
	BASE max = v->data[0 * stride];
	BASE min = v->data[0 * stride];

	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
		{
			min = x;
			imin = i;
		}
		if (x > max)
		{
			max = x;
			imax = i;
		}
#ifdef FP
		if (isnan(x))
		{
			imin = i;
			imax = i;
			break;
		}
#endif
	}

	*imin_out = imin;
	*imax_out = imax;
}
#include "../templates_off.h"
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
#include "../templates_on.h"
BASE
FUNCTION(gsl_vector, max) (const TYPE(gsl_vector) * v)
{
	/* finds the largest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE max = v->data[0 * stride];
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x > max)
			max = x;
#ifdef FP
		if (isnan(x))
			return x;
#endif
	}

	return max;
}

BASE
FUNCTION(gsl_vector, min) (const TYPE(gsl_vector) * v)
{
	/* finds the smallest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE min = v->data[0 * stride];
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
			min = x;
#ifdef FP
		if (isnan(x))
			return x;
#endif
	}

	return min;
}

void
FUNCTION(gsl_vector, minmax) (const TYPE(gsl_vector) * v,
BASE * min_out,
BASE * max_out)
{
	/* finds the smallest and largest elements of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE max = v->data[0 * stride];
	BASE min = v->data[0 * stride];

	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
		{
			min = x;
		}
		if (x > max)
		{
			max = x;
		}
#ifdef FP
		if (isnan(x))
		{
			min = x;
			max = x;
			break;
		}
#endif
	}

	*min_out = min;
	*max_out = max;
}


size_t
FUNCTION(gsl_vector, max_index) (const TYPE(gsl_vector) * v)
{
	/* finds the largest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE max = v->data[0 * stride];
	size_t imax = 0;
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x > max)
		{
			max = x;
			imax = i;
		}
#ifdef FP
		if (isnan(x))
		{
			return i;
		}
#endif
	}

	return imax;
}

size_t
FUNCTION(gsl_vector, min_index) (const TYPE(gsl_vector) * v)
{
	/* finds the smallest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE min = v->data[0 * stride];
	size_t imin = 0;
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
		{
			min = x;
			imin = i;
		}
#ifdef FP
		if (isnan(x))
		{
			return i;
		}
#endif
	}

	return imin;
}


void
FUNCTION(gsl_vector, minmax_index) (const TYPE(gsl_vector) * v,
size_t * imin_out,
size_t * imax_out)
{
	/* finds the smallest and largest elements of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	size_t imin = 0, imax = 0;
	BASE max = v->data[0 * stride];
	BASE min = v->data[0 * stride];

	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
		{
			min = x;
			imin = i;
		}
		if (x > max)
		{
			max = x;
			imax = i;
		}
#ifdef FP
		if (isnan(x))
		{
			imin = i;
			imax = i;
			break;
		}
#endif
	}

	*imin_out = imin;
	*imax_out = imax;
}
#include "../templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "../templates_on.h"
BASE
FUNCTION(gsl_vector, max) (const TYPE(gsl_vector) * v)
{
	/* finds the largest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE max = v->data[0 * stride];
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x > max)
			max = x;
#ifdef FP
		if (isnan(x))
			return x;
#endif
	}

	return max;
}

BASE
FUNCTION(gsl_vector, min) (const TYPE(gsl_vector) * v)
{
	/* finds the smallest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE min = v->data[0 * stride];
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
			min = x;
#ifdef FP
		if (isnan(x))
			return x;
#endif
	}

	return min;
}

void
FUNCTION(gsl_vector, minmax) (const TYPE(gsl_vector) * v,
BASE * min_out,
BASE * max_out)
{
	/* finds the smallest and largest elements of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE max = v->data[0 * stride];
	BASE min = v->data[0 * stride];

	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
		{
			min = x;
		}
		if (x > max)
		{
			max = x;
		}
#ifdef FP
		if (isnan(x))
		{
			min = x;
			max = x;
			break;
		}
#endif
	}

	*min_out = min;
	*max_out = max;
}


size_t
FUNCTION(gsl_vector, max_index) (const TYPE(gsl_vector) * v)
{
	/* finds the largest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE max = v->data[0 * stride];
	size_t imax = 0;
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x > max)
		{
			max = x;
			imax = i;
		}
#ifdef FP
		if (isnan(x))
		{
			return i;
		}
#endif
	}

	return imax;
}

size_t
FUNCTION(gsl_vector, min_index) (const TYPE(gsl_vector) * v)
{
	/* finds the smallest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE min = v->data[0 * stride];
	size_t imin = 0;
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
		{
			min = x;
			imin = i;
		}
#ifdef FP
		if (isnan(x))
		{
			return i;
		}
#endif
	}

	return imin;
}


void
FUNCTION(gsl_vector, minmax_index) (const TYPE(gsl_vector) * v,
size_t * imin_out,
size_t * imax_out)
{
	/* finds the smallest and largest elements of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	size_t imin = 0, imax = 0;
	BASE max = v->data[0 * stride];
	BASE min = v->data[0 * stride];

	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
		{
			min = x;
			imin = i;
		}
		if (x > max)
		{
			max = x;
			imax = i;
		}
#ifdef FP
		if (isnan(x))
		{
			imin = i;
			imax = i;
			break;
		}
#endif
	}

	*imin_out = imin;
	*imax_out = imax;
}
#include "../templates_off.h"
#undef  BASE_FLOAT

#define BASE_ULONG
#include "../templates_on.h"
BASE
FUNCTION(gsl_vector, max) (const TYPE(gsl_vector) * v)
{
	/* finds the largest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE max = v->data[0 * stride];
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x > max)
			max = x;
#ifdef FP
		if (isnan(x))
			return x;
#endif
	}

	return max;
}

BASE
FUNCTION(gsl_vector, min) (const TYPE(gsl_vector) * v)
{
	/* finds the smallest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE min = v->data[0 * stride];
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
			min = x;
#ifdef FP
		if (isnan(x))
			return x;
#endif
	}

	return min;
}

void
FUNCTION(gsl_vector, minmax) (const TYPE(gsl_vector) * v,
BASE * min_out,
BASE * max_out)
{
	/* finds the smallest and largest elements of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE max = v->data[0 * stride];
	BASE min = v->data[0 * stride];

	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
		{
			min = x;
		}
		if (x > max)
		{
			max = x;
		}
#ifdef FP
		if (isnan(x))
		{
			min = x;
			max = x;
			break;
		}
#endif
	}

	*min_out = min;
	*max_out = max;
}


size_t
FUNCTION(gsl_vector, max_index) (const TYPE(gsl_vector) * v)
{
	/* finds the largest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE max = v->data[0 * stride];
	size_t imax = 0;
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x > max)
		{
			max = x;
			imax = i;
		}
#ifdef FP
		if (isnan(x))
		{
			return i;
		}
#endif
	}

	return imax;
}

size_t
FUNCTION(gsl_vector, min_index) (const TYPE(gsl_vector) * v)
{
	/* finds the smallest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE min = v->data[0 * stride];
	size_t imin = 0;
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
		{
			min = x;
			imin = i;
		}
#ifdef FP
		if (isnan(x))
		{
			return i;
		}
#endif
	}

	return imin;
}


void
FUNCTION(gsl_vector, minmax_index) (const TYPE(gsl_vector) * v,
size_t * imin_out,
size_t * imax_out)
{
	/* finds the smallest and largest elements of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	size_t imin = 0, imax = 0;
	BASE max = v->data[0 * stride];
	BASE min = v->data[0 * stride];

	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
		{
			min = x;
			imin = i;
		}
		if (x > max)
		{
			max = x;
			imax = i;
		}
#ifdef FP
		if (isnan(x))
		{
			imin = i;
			imax = i;
			break;
		}
#endif
	}

	*imin_out = imin;
	*imax_out = imax;
}
#include "../templates_off.h"
#undef  BASE_ULONG

#define BASE_LONG
#include "../templates_on.h"
BASE
FUNCTION(gsl_vector, max) (const TYPE(gsl_vector) * v)
{
	/* finds the largest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE max = v->data[0 * stride];
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x > max)
			max = x;
#ifdef FP
		if (isnan(x))
			return x;
#endif
	}

	return max;
}

BASE
FUNCTION(gsl_vector, min) (const TYPE(gsl_vector) * v)
{
	/* finds the smallest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE min = v->data[0 * stride];
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
			min = x;
#ifdef FP
		if (isnan(x))
			return x;
#endif
	}

	return min;
}

void
FUNCTION(gsl_vector, minmax) (const TYPE(gsl_vector) * v,
BASE * min_out,
BASE * max_out)
{
	/* finds the smallest and largest elements of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE max = v->data[0 * stride];
	BASE min = v->data[0 * stride];

	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
		{
			min = x;
		}
		if (x > max)
		{
			max = x;
		}
#ifdef FP
		if (isnan(x))
		{
			min = x;
			max = x;
			break;
		}
#endif
	}

	*min_out = min;
	*max_out = max;
}


size_t
FUNCTION(gsl_vector, max_index) (const TYPE(gsl_vector) * v)
{
	/* finds the largest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE max = v->data[0 * stride];
	size_t imax = 0;
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x > max)
		{
			max = x;
			imax = i;
		}
#ifdef FP
		if (isnan(x))
		{
			return i;
		}
#endif
	}

	return imax;
}

size_t
FUNCTION(gsl_vector, min_index) (const TYPE(gsl_vector) * v)
{
	/* finds the smallest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE min = v->data[0 * stride];
	size_t imin = 0;
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
		{
			min = x;
			imin = i;
		}
#ifdef FP
		if (isnan(x))
		{
			return i;
		}
#endif
	}

	return imin;
}


void
FUNCTION(gsl_vector, minmax_index) (const TYPE(gsl_vector) * v,
size_t * imin_out,
size_t * imax_out)
{
	/* finds the smallest and largest elements of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	size_t imin = 0, imax = 0;
	BASE max = v->data[0 * stride];
	BASE min = v->data[0 * stride];

	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
		{
			min = x;
			imin = i;
		}
		if (x > max)
		{
			max = x;
			imax = i;
		}
#ifdef FP
		if (isnan(x))
		{
			imin = i;
			imax = i;
			break;
		}
#endif
	}

	*imin_out = imin;
	*imax_out = imax;
}
#include "../templates_off.h"
#undef  BASE_LONG

#define BASE_UINT
#include "../templates_on.h"
BASE
FUNCTION(gsl_vector, max) (const TYPE(gsl_vector) * v)
{
	/* finds the largest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE max = v->data[0 * stride];
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x > max)
			max = x;
#ifdef FP
		if (isnan(x))
			return x;
#endif
	}

	return max;
}

BASE
FUNCTION(gsl_vector, min) (const TYPE(gsl_vector) * v)
{
	/* finds the smallest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE min = v->data[0 * stride];
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
			min = x;
#ifdef FP
		if (isnan(x))
			return x;
#endif
	}

	return min;
}

void
FUNCTION(gsl_vector, minmax) (const TYPE(gsl_vector) * v,
BASE * min_out,
BASE * max_out)
{
	/* finds the smallest and largest elements of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE max = v->data[0 * stride];
	BASE min = v->data[0 * stride];

	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
		{
			min = x;
		}
		if (x > max)
		{
			max = x;
		}
#ifdef FP
		if (isnan(x))
		{
			min = x;
			max = x;
			break;
		}
#endif
	}

	*min_out = min;
	*max_out = max;
}


size_t
FUNCTION(gsl_vector, max_index) (const TYPE(gsl_vector) * v)
{
	/* finds the largest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE max = v->data[0 * stride];
	size_t imax = 0;
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x > max)
		{
			max = x;
			imax = i;
		}
#ifdef FP
		if (isnan(x))
		{
			return i;
		}
#endif
	}

	return imax;
}

size_t
FUNCTION(gsl_vector, min_index) (const TYPE(gsl_vector) * v)
{
	/* finds the smallest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE min = v->data[0 * stride];
	size_t imin = 0;
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
		{
			min = x;
			imin = i;
		}
#ifdef FP
		if (isnan(x))
		{
			return i;
		}
#endif
	}

	return imin;
}


void
FUNCTION(gsl_vector, minmax_index) (const TYPE(gsl_vector) * v,
size_t * imin_out,
size_t * imax_out)
{
	/* finds the smallest and largest elements of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	size_t imin = 0, imax = 0;
	BASE max = v->data[0 * stride];
	BASE min = v->data[0 * stride];

	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
		{
			min = x;
			imin = i;
		}
		if (x > max)
		{
			max = x;
			imax = i;
		}
#ifdef FP
		if (isnan(x))
		{
			imin = i;
			imax = i;
			break;
		}
#endif
	}

	*imin_out = imin;
	*imax_out = imax;
}
#include "../templates_off.h"
#undef  BASE_UINT

#define BASE_INT
#include "../templates_on.h"
BASE
FUNCTION(gsl_vector, max) (const TYPE(gsl_vector) * v)
{
	/* finds the largest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE max = v->data[0 * stride];
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x > max)
			max = x;
#ifdef FP
		if (isnan(x))
			return x;
#endif
	}

	return max;
}

BASE
FUNCTION(gsl_vector, min) (const TYPE(gsl_vector) * v)
{
	/* finds the smallest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE min = v->data[0 * stride];
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
			min = x;
#ifdef FP
		if (isnan(x))
			return x;
#endif
	}

	return min;
}

void
FUNCTION(gsl_vector, minmax) (const TYPE(gsl_vector) * v,
BASE * min_out,
BASE * max_out)
{
	/* finds the smallest and largest elements of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE max = v->data[0 * stride];
	BASE min = v->data[0 * stride];

	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
		{
			min = x;
		}
		if (x > max)
		{
			max = x;
		}
#ifdef FP
		if (isnan(x))
		{
			min = x;
			max = x;
			break;
		}
#endif
	}

	*min_out = min;
	*max_out = max;
}


size_t
FUNCTION(gsl_vector, max_index) (const TYPE(gsl_vector) * v)
{
	/* finds the largest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE max = v->data[0 * stride];
	size_t imax = 0;
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x > max)
		{
			max = x;
			imax = i;
		}
#ifdef FP
		if (isnan(x))
		{
			return i;
		}
#endif
	}

	return imax;
}

size_t
FUNCTION(gsl_vector, min_index) (const TYPE(gsl_vector) * v)
{
	/* finds the smallest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE min = v->data[0 * stride];
	size_t imin = 0;
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
		{
			min = x;
			imin = i;
		}
#ifdef FP
		if (isnan(x))
		{
			return i;
		}
#endif
	}

	return imin;
}


void
FUNCTION(gsl_vector, minmax_index) (const TYPE(gsl_vector) * v,
size_t * imin_out,
size_t * imax_out)
{
	/* finds the smallest and largest elements of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	size_t imin = 0, imax = 0;
	BASE max = v->data[0 * stride];
	BASE min = v->data[0 * stride];

	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
		{
			min = x;
			imin = i;
		}
		if (x > max)
		{
			max = x;
			imax = i;
		}
#ifdef FP
		if (isnan(x))
		{
			imin = i;
			imax = i;
			break;
		}
#endif
	}

	*imin_out = imin;
	*imax_out = imax;
}
#include "../templates_off.h"
#undef  BASE_INT

#define BASE_USHORT
#include "../templates_on.h"
BASE
FUNCTION(gsl_vector, max) (const TYPE(gsl_vector) * v)
{
	/* finds the largest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE max = v->data[0 * stride];
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x > max)
			max = x;
#ifdef FP
		if (isnan(x))
			return x;
#endif
	}

	return max;
}

BASE
FUNCTION(gsl_vector, min) (const TYPE(gsl_vector) * v)
{
	/* finds the smallest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE min = v->data[0 * stride];
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
			min = x;
#ifdef FP
		if (isnan(x))
			return x;
#endif
	}

	return min;
}

void
FUNCTION(gsl_vector, minmax) (const TYPE(gsl_vector) * v,
BASE * min_out,
BASE * max_out)
{
	/* finds the smallest and largest elements of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE max = v->data[0 * stride];
	BASE min = v->data[0 * stride];

	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
		{
			min = x;
		}
		if (x > max)
		{
			max = x;
		}
#ifdef FP
		if (isnan(x))
		{
			min = x;
			max = x;
			break;
		}
#endif
	}

	*min_out = min;
	*max_out = max;
}


size_t
FUNCTION(gsl_vector, max_index) (const TYPE(gsl_vector) * v)
{
	/* finds the largest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE max = v->data[0 * stride];
	size_t imax = 0;
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x > max)
		{
			max = x;
			imax = i;
		}
#ifdef FP
		if (isnan(x))
		{
			return i;
		}
#endif
	}

	return imax;
}

size_t
FUNCTION(gsl_vector, min_index) (const TYPE(gsl_vector) * v)
{
	/* finds the smallest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE min = v->data[0 * stride];
	size_t imin = 0;
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
		{
			min = x;
			imin = i;
		}
#ifdef FP
		if (isnan(x))
		{
			return i;
		}
#endif
	}

	return imin;
}


void
FUNCTION(gsl_vector, minmax_index) (const TYPE(gsl_vector) * v,
size_t * imin_out,
size_t * imax_out)
{
	/* finds the smallest and largest elements of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	size_t imin = 0, imax = 0;
	BASE max = v->data[0 * stride];
	BASE min = v->data[0 * stride];

	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
		{
			min = x;
			imin = i;
		}
		if (x > max)
		{
			max = x;
			imax = i;
		}
#ifdef FP
		if (isnan(x))
		{
			imin = i;
			imax = i;
			break;
		}
#endif
	}

	*imin_out = imin;
	*imax_out = imax;
}
#include "../templates_off.h"
#undef  BASE_USHORT

#define BASE_SHORT
#include "../templates_on.h"
BASE
FUNCTION(gsl_vector, max) (const TYPE(gsl_vector) * v)
{
	/* finds the largest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE max = v->data[0 * stride];
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x > max)
			max = x;
#ifdef FP
		if (isnan(x))
			return x;
#endif
	}

	return max;
}

BASE
FUNCTION(gsl_vector, min) (const TYPE(gsl_vector) * v)
{
	/* finds the smallest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE min = v->data[0 * stride];
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
			min = x;
#ifdef FP
		if (isnan(x))
			return x;
#endif
	}

	return min;
}

void
FUNCTION(gsl_vector, minmax) (const TYPE(gsl_vector) * v,
BASE * min_out,
BASE * max_out)
{
	/* finds the smallest and largest elements of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE max = v->data[0 * stride];
	BASE min = v->data[0 * stride];

	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
		{
			min = x;
		}
		if (x > max)
		{
			max = x;
		}
#ifdef FP
		if (isnan(x))
		{
			min = x;
			max = x;
			break;
		}
#endif
	}

	*min_out = min;
	*max_out = max;
}


size_t
FUNCTION(gsl_vector, max_index) (const TYPE(gsl_vector) * v)
{
	/* finds the largest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE max = v->data[0 * stride];
	size_t imax = 0;
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x > max)
		{
			max = x;
			imax = i;
		}
#ifdef FP
		if (isnan(x))
		{
			return i;
		}
#endif
	}

	return imax;
}

size_t
FUNCTION(gsl_vector, min_index) (const TYPE(gsl_vector) * v)
{
	/* finds the smallest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE min = v->data[0 * stride];
	size_t imin = 0;
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
		{
			min = x;
			imin = i;
		}
#ifdef FP
		if (isnan(x))
		{
			return i;
		}
#endif
	}

	return imin;
}


void
FUNCTION(gsl_vector, minmax_index) (const TYPE(gsl_vector) * v,
size_t * imin_out,
size_t * imax_out)
{
	/* finds the smallest and largest elements of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	size_t imin = 0, imax = 0;
	BASE max = v->data[0 * stride];
	BASE min = v->data[0 * stride];

	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
		{
			min = x;
			imin = i;
		}
		if (x > max)
		{
			max = x;
			imax = i;
		}
#ifdef FP
		if (isnan(x))
		{
			imin = i;
			imax = i;
			break;
		}
#endif
	}

	*imin_out = imin;
	*imax_out = imax;
}
#include "../templates_off.h"
#undef  BASE_SHORT

#define BASE_UCHAR
#include "../templates_on.h"
BASE
FUNCTION(gsl_vector, max) (const TYPE(gsl_vector) * v)
{
	/* finds the largest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE max = v->data[0 * stride];
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x > max)
			max = x;
#ifdef FP
		if (isnan(x))
			return x;
#endif
	}

	return max;
}

BASE
FUNCTION(gsl_vector, min) (const TYPE(gsl_vector) * v)
{
	/* finds the smallest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE min = v->data[0 * stride];
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
			min = x;
#ifdef FP
		if (isnan(x))
			return x;
#endif
	}

	return min;
}

void
FUNCTION(gsl_vector, minmax) (const TYPE(gsl_vector) * v,
BASE * min_out,
BASE * max_out)
{
	/* finds the smallest and largest elements of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE max = v->data[0 * stride];
	BASE min = v->data[0 * stride];

	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
		{
			min = x;
		}
		if (x > max)
		{
			max = x;
		}
#ifdef FP
		if (isnan(x))
		{
			min = x;
			max = x;
			break;
		}
#endif
	}

	*min_out = min;
	*max_out = max;
}


size_t
FUNCTION(gsl_vector, max_index) (const TYPE(gsl_vector) * v)
{
	/* finds the largest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE max = v->data[0 * stride];
	size_t imax = 0;
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x > max)
		{
			max = x;
			imax = i;
		}
#ifdef FP
		if (isnan(x))
		{
			return i;
		}
#endif
	}

	return imax;
}

size_t
FUNCTION(gsl_vector, min_index) (const TYPE(gsl_vector) * v)
{
	/* finds the smallest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE min = v->data[0 * stride];
	size_t imin = 0;
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
		{
			min = x;
			imin = i;
		}
#ifdef FP
		if (isnan(x))
		{
			return i;
		}
#endif
	}

	return imin;
}


void
FUNCTION(gsl_vector, minmax_index) (const TYPE(gsl_vector) * v,
size_t * imin_out,
size_t * imax_out)
{
	/* finds the smallest and largest elements of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	size_t imin = 0, imax = 0;
	BASE max = v->data[0 * stride];
	BASE min = v->data[0 * stride];

	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
		{
			min = x;
			imin = i;
		}
		if (x > max)
		{
			max = x;
			imax = i;
		}
#ifdef FP
		if (isnan(x))
		{
			imin = i;
			imax = i;
			break;
		}
#endif
	}

	*imin_out = imin;
	*imax_out = imax;
}
#include "../templates_off.h"
#undef  BASE_UCHAR

#define BASE_CHAR
#include "../templates_on.h"
BASE
FUNCTION(gsl_vector, max) (const TYPE(gsl_vector) * v)
{
	/* finds the largest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE max = v->data[0 * stride];
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x > max)
			max = x;
#ifdef FP
		if (isnan(x))
			return x;
#endif
	}

	return max;
}

BASE
FUNCTION(gsl_vector, min) (const TYPE(gsl_vector) * v)
{
	/* finds the smallest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE min = v->data[0 * stride];
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
			min = x;
#ifdef FP
		if (isnan(x))
			return x;
#endif
	}

	return min;
}

void
FUNCTION(gsl_vector, minmax) (const TYPE(gsl_vector) * v,
BASE * min_out,
BASE * max_out)
{
	/* finds the smallest and largest elements of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE max = v->data[0 * stride];
	BASE min = v->data[0 * stride];

	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
		{
			min = x;
		}
		if (x > max)
		{
			max = x;
		}
#ifdef FP
		if (isnan(x))
		{
			min = x;
			max = x;
			break;
		}
#endif
	}

	*min_out = min;
	*max_out = max;
}


size_t
FUNCTION(gsl_vector, max_index) (const TYPE(gsl_vector) * v)
{
	/* finds the largest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE max = v->data[0 * stride];
	size_t imax = 0;
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x > max)
		{
			max = x;
			imax = i;
		}
#ifdef FP
		if (isnan(x))
		{
			return i;
		}
#endif
	}

	return imax;
}

size_t
FUNCTION(gsl_vector, min_index) (const TYPE(gsl_vector) * v)
{
	/* finds the smallest element of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	BASE min = v->data[0 * stride];
	size_t imin = 0;
	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
		{
			min = x;
			imin = i;
		}
#ifdef FP
		if (isnan(x))
		{
			return i;
		}
#endif
	}

	return imin;
}


void
FUNCTION(gsl_vector, minmax_index) (const TYPE(gsl_vector) * v,
size_t * imin_out,
size_t * imax_out)
{
	/* finds the smallest and largest elements of a vector */

	const size_t N = v->size;
	const size_t stride = v->stride;

	size_t imin = 0, imax = 0;
	BASE max = v->data[0 * stride];
	BASE min = v->data[0 * stride];

	size_t i;

	for (i = 0; i < N; i++)
	{
		BASE x = v->data[i*stride];
		if (x < min)
		{
			min = x;
			imin = i;
		}
		if (x > max)
		{
			max = x;
			imax = i;
		}
#ifdef FP
		if (isnan(x))
		{
			imin = i;
			imax = i;
			break;
		}
#endif
	}

	*imin_out = imin;
	*imax_out = imax;
}
#include "../templates_off.h"
#undef  BASE_CHAR


