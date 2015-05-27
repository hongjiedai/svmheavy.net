#include "stdafx.h"
#include "config.h.in"
#include <stdlib.h>
#include <gsl/vector/gsl_vector.h>

#define BASE_LONG_DOUBLE
#include "../templates_on.h"
int
FUNCTION(gsl_vector, add) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] += b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, sub) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] -= b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, mul) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] *= b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, div) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] /= b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, scale) (TYPE(gsl_vector) * a, const double x)
{
	const size_t N = a->size;
	const size_t stride = a->stride;

	size_t i;

	for (i = 0; i < N; i++)
	{
		a->data[i * stride] *= x;
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_vector, add_constant) (TYPE(gsl_vector) * a, const double x)
{
	const size_t N = a->size;
	const size_t stride = a->stride;

	size_t i;

	for (i = 0; i < N; i++)
	{
		a->data[i * stride] += x;
	}

	return GSL_SUCCESS;
}

#include "../templates_off.h"
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
#include "../templates_on.h"
int
FUNCTION(gsl_vector, add) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] += b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, sub) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] -= b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, mul) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] *= b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, div) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] /= b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, scale) (TYPE(gsl_vector) * a, const double x)
{
	const size_t N = a->size;
	const size_t stride = a->stride;

	size_t i;

	for (i = 0; i < N; i++)
	{
		a->data[i * stride] *= x;
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_vector, add_constant) (TYPE(gsl_vector) * a, const double x)
{
	const size_t N = a->size;
	const size_t stride = a->stride;

	size_t i;

	for (i = 0; i < N; i++)
	{
		a->data[i * stride] += x;
	}

	return GSL_SUCCESS;
}

#include "../templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "../templates_on.h"
int
FUNCTION(gsl_vector, add) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] += b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, sub) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] -= b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, mul) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] *= b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, div) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] /= b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, scale) (TYPE(gsl_vector) * a, const double x)
{
	const size_t N = a->size;
	const size_t stride = a->stride;

	size_t i;

	for (i = 0; i < N; i++)
	{
		a->data[i * stride] *= x;
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_vector, add_constant) (TYPE(gsl_vector) * a, const double x)
{
	const size_t N = a->size;
	const size_t stride = a->stride;

	size_t i;

	for (i = 0; i < N; i++)
	{
		a->data[i * stride] += x;
	}

	return GSL_SUCCESS;
}

#include "../templates_off.h"
#undef  BASE_FLOAT

#define BASE_ULONG
#include "../templates_on.h"
int
FUNCTION(gsl_vector, add) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] += b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, sub) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] -= b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, mul) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] *= b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, div) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] /= b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, scale) (TYPE(gsl_vector) * a, const double x)
{
	const size_t N = a->size;
	const size_t stride = a->stride;

	size_t i;

	for (i = 0; i < N; i++)
	{
		a->data[i * stride] *= x;
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_vector, add_constant) (TYPE(gsl_vector) * a, const double x)
{
	const size_t N = a->size;
	const size_t stride = a->stride;

	size_t i;

	for (i = 0; i < N; i++)
	{
		a->data[i * stride] += x;
	}

	return GSL_SUCCESS;
}

#include "../templates_off.h"
#undef  BASE_ULONG

#define BASE_LONG
#include "../templates_on.h"
int
FUNCTION(gsl_vector, add) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] += b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, sub) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] -= b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, mul) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] *= b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, div) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] /= b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, scale) (TYPE(gsl_vector) * a, const double x)
{
	const size_t N = a->size;
	const size_t stride = a->stride;

	size_t i;

	for (i = 0; i < N; i++)
	{
		a->data[i * stride] *= x;
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_vector, add_constant) (TYPE(gsl_vector) * a, const double x)
{
	const size_t N = a->size;
	const size_t stride = a->stride;

	size_t i;

	for (i = 0; i < N; i++)
	{
		a->data[i * stride] += x;
	}

	return GSL_SUCCESS;
}

#include "../templates_off.h"
#undef  BASE_LONG

#define BASE_UINT
#include "../templates_on.h"
int
FUNCTION(gsl_vector, add) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] += b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, sub) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] -= b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, mul) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] *= b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, div) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] /= b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, scale) (TYPE(gsl_vector) * a, const double x)
{
	const size_t N = a->size;
	const size_t stride = a->stride;

	size_t i;

	for (i = 0; i < N; i++)
	{
		a->data[i * stride] *= x;
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_vector, add_constant) (TYPE(gsl_vector) * a, const double x)
{
	const size_t N = a->size;
	const size_t stride = a->stride;

	size_t i;

	for (i = 0; i < N; i++)
	{
		a->data[i * stride] += x;
	}

	return GSL_SUCCESS;
}

#include "../templates_off.h"
#undef  BASE_UINT

#define BASE_INT
#include "../templates_on.h"
int
FUNCTION(gsl_vector, add) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] += b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, sub) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] -= b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, mul) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] *= b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, div) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] /= b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, scale) (TYPE(gsl_vector) * a, const double x)
{
	const size_t N = a->size;
	const size_t stride = a->stride;

	size_t i;

	for (i = 0; i < N; i++)
	{
		a->data[i * stride] *= x;
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_vector, add_constant) (TYPE(gsl_vector) * a, const double x)
{
	const size_t N = a->size;
	const size_t stride = a->stride;

	size_t i;

	for (i = 0; i < N; i++)
	{
		a->data[i * stride] += x;
	}

	return GSL_SUCCESS;
}

#include "../templates_off.h"
#undef  BASE_INT

#define BASE_USHORT
#include "../templates_on.h"
int
FUNCTION(gsl_vector, add) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] += b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, sub) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] -= b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, mul) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] *= b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, div) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] /= b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, scale) (TYPE(gsl_vector) * a, const double x)
{
	const size_t N = a->size;
	const size_t stride = a->stride;

	size_t i;

	for (i = 0; i < N; i++)
	{
		a->data[i * stride] *= x;
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_vector, add_constant) (TYPE(gsl_vector) * a, const double x)
{
	const size_t N = a->size;
	const size_t stride = a->stride;

	size_t i;

	for (i = 0; i < N; i++)
	{
		a->data[i * stride] += x;
	}

	return GSL_SUCCESS;
}

#include "../templates_off.h"
#undef  BASE_USHORT

#define BASE_SHORT
#include "../templates_on.h"
int
FUNCTION(gsl_vector, add) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] += b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, sub) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] -= b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, mul) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] *= b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, div) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] /= b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, scale) (TYPE(gsl_vector) * a, const double x)
{
	const size_t N = a->size;
	const size_t stride = a->stride;

	size_t i;

	for (i = 0; i < N; i++)
	{
		a->data[i * stride] *= x;
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_vector, add_constant) (TYPE(gsl_vector) * a, const double x)
{
	const size_t N = a->size;
	const size_t stride = a->stride;

	size_t i;

	for (i = 0; i < N; i++)
	{
		a->data[i * stride] += x;
	}

	return GSL_SUCCESS;
}

#include "../templates_off.h"
#undef  BASE_SHORT

#define BASE_UCHAR
#include "../templates_on.h"
int
FUNCTION(gsl_vector, add) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] += b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, sub) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] -= b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, mul) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] *= b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, div) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] /= b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, scale) (TYPE(gsl_vector) * a, const double x)
{
	const size_t N = a->size;
	const size_t stride = a->stride;

	size_t i;

	for (i = 0; i < N; i++)
	{
		a->data[i * stride] *= x;
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_vector, add_constant) (TYPE(gsl_vector) * a, const double x)
{
	const size_t N = a->size;
	const size_t stride = a->stride;

	size_t i;

	for (i = 0; i < N; i++)
	{
		a->data[i * stride] += x;
	}

	return GSL_SUCCESS;
}

#include "../templates_off.h"
#undef  BASE_UCHAR

#define BASE_CHAR
#include "../templates_on.h"
int
FUNCTION(gsl_vector, add) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] += b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, sub) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] -= b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, mul) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] *= b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, div) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
	const size_t N = a->size;

	if (b->size != N)
	{
		GSL_ERROR("vectors must have same length", GSL_EBADLEN);
	}
	else
	{
		const size_t stride_a = a->stride;
		const size_t stride_b = b->stride;

		size_t i;

		for (i = 0; i < N; i++)
		{
			a->data[i * stride_a] /= b->data[i * stride_b];
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_vector, scale) (TYPE(gsl_vector) * a, const double x)
{
	const size_t N = a->size;
	const size_t stride = a->stride;

	size_t i;

	for (i = 0; i < N; i++)
	{
		a->data[i * stride] *= x;
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_vector, add_constant) (TYPE(gsl_vector) * a, const double x)
{
	const size_t N = a->size;
	const size_t stride = a->stride;

	size_t i;

	for (i = 0; i < N; i++)
	{
		a->data[i * stride] += x;
	}

	return GSL_SUCCESS;
}
#include "../templates_off.h"
#undef  BASE_CHAR


