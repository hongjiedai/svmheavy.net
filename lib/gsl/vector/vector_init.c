#include "stdafx.h"
#include "config.h.in"
#include <stdlib.h>
#include <gsl/vector/gsl_vector.h>

#define BASE_GSL_COMPLEX_LONG
#include "../templates_on.h"
TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc) (const size_t n)
{
	TYPE(gsl_block) * block;
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	block = FUNCTION(gsl_block, alloc) (n);

	if (block == 0)
	{
		free(v);

		GSL_ERROR_VAL("failed to allocate space for block",
			GSL_ENOMEM, 0);
	}

	v->data = block->data;
	v->size = n;
	v->stride = 1;
	v->block = block;
	v->owner = 1;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, calloc) (const size_t n)
{
	size_t i;

	TYPE(gsl_vector) * v = FUNCTION(gsl_vector, alloc) (n);

	if (v == 0)
		return 0;

	/* initialize vector to zero */

	for (i = 0; i < MULTIPLICITY * n; i++)
	{
		v->data[i] = 0;
	}

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_from_block) (TYPE(gsl_block) * block,
const size_t offset,
const size_t n,
const size_t stride)
{
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer", GSL_EINVAL, 0);
	}

	if (block->size <= offset + (n - 1) * stride)
	{
		GSL_ERROR_VAL("vector would extend past end of block", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = block->data + MULTIPLICITY * offset;
	v->size = n;
	v->stride = stride;
	v->block = block;
	v->owner = 0;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_from_vector) (TYPE(gsl_vector) * w,
const size_t offset,
const size_t n,
const size_t stride)
{
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer", GSL_EINVAL, 0);
	}

	if (offset + (n - 1) * stride >= w->size)
	{
		GSL_ERROR_VAL("vector would extend past end of block", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = w->data + MULTIPLICITY * w->stride * offset;
	v->size = n;
	v->stride = stride * w->stride;
	v->block = w->block;
	v->owner = 0;

	return v;
}


void
FUNCTION(gsl_vector, free) (TYPE(gsl_vector) * v)
{
	if (v->owner)
	{
		FUNCTION(gsl_block, free) (v->block);
	}
	free(v);
}


void
FUNCTION(gsl_vector, set_all) (TYPE(gsl_vector) * v, BASE x)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;

	size_t i;

	for (i = 0; i < n; i++)
	{
		*(BASE *)(data + MULTIPLICITY * i * stride) = x;
	}
}

void
FUNCTION(gsl_vector, set_zero) (TYPE(gsl_vector) * v)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;
	const BASE zero = ZERO;

	size_t i;

	for (i = 0; i < n; i++)
	{
		*(BASE *)(data + MULTIPLICITY * i * stride) = zero;
	}
}

int
FUNCTION(gsl_vector, set_basis) (TYPE(gsl_vector) * v, size_t i)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;
	const BASE zero = ZERO;
	const BASE one = ONE;

	size_t k;

	if (i >= n)
	{
		GSL_ERROR("index out of range", GSL_EINVAL);
	}

	for (k = 0; k < n; k++)
	{
		*(BASE *)(data + MULTIPLICITY * k * stride) = zero;
	}

	*(BASE *)(data + MULTIPLICITY * i * stride) = one;

	return GSL_SUCCESS;
}
#include "../templates_off.h"
#undef  BASE_GSL_COMPLEX_LONG

#define BASE_GSL_COMPLEX
#include "../templates_on.h"
TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc) (const size_t n)
{
	TYPE(gsl_block) * block;
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	block = FUNCTION(gsl_block, alloc) (n);

	if (block == 0)
	{
		free(v);

		GSL_ERROR_VAL("failed to allocate space for block",
			GSL_ENOMEM, 0);
	}

	v->data = block->data;
	v->size = n;
	v->stride = 1;
	v->block = block;
	v->owner = 1;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, calloc) (const size_t n)
{
	size_t i;

	TYPE(gsl_vector) * v = FUNCTION(gsl_vector, alloc) (n);

	if (v == 0)
		return 0;

	/* initialize vector to zero */

	for (i = 0; i < MULTIPLICITY * n; i++)
	{
		v->data[i] = 0;
	}

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_from_block) (TYPE(gsl_block) * block,
const size_t offset,
const size_t n,
const size_t stride)
{
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer", GSL_EINVAL, 0);
	}

	if (block->size <= offset + (n - 1) * stride)
	{
		GSL_ERROR_VAL("vector would extend past end of block", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = block->data + MULTIPLICITY * offset;
	v->size = n;
	v->stride = stride;
	v->block = block;
	v->owner = 0;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_from_vector) (TYPE(gsl_vector) * w,
const size_t offset,
const size_t n,
const size_t stride)
{
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer", GSL_EINVAL, 0);
	}

	if (offset + (n - 1) * stride >= w->size)
	{
		GSL_ERROR_VAL("vector would extend past end of block", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = w->data + MULTIPLICITY * w->stride * offset;
	v->size = n;
	v->stride = stride * w->stride;
	v->block = w->block;
	v->owner = 0;

	return v;
}


void
FUNCTION(gsl_vector, free) (TYPE(gsl_vector) * v)
{
	if (v->owner)
	{
		FUNCTION(gsl_block, free) (v->block);
	}
	free(v);
}


void
FUNCTION(gsl_vector, set_all) (TYPE(gsl_vector) * v, BASE x)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;

	size_t i;

	for (i = 0; i < n; i++)
	{
		*(BASE *)(data + MULTIPLICITY * i * stride) = x;
	}
}

void
FUNCTION(gsl_vector, set_zero) (TYPE(gsl_vector) * v)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;
	const BASE zero = ZERO;

	size_t i;

	for (i = 0; i < n; i++)
	{
		*(BASE *)(data + MULTIPLICITY * i * stride) = zero;
	}
}

int
FUNCTION(gsl_vector, set_basis) (TYPE(gsl_vector) * v, size_t i)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;
	const BASE zero = ZERO;
	const BASE one = ONE;

	size_t k;

	if (i >= n)
	{
		GSL_ERROR("index out of range", GSL_EINVAL);
	}

	for (k = 0; k < n; k++)
	{
		*(BASE *)(data + MULTIPLICITY * k * stride) = zero;
	}

	*(BASE *)(data + MULTIPLICITY * i * stride) = one;

	return GSL_SUCCESS;
}
#include "../templates_off.h"
#undef  BASE_GSL_COMPLEX

#define BASE_GSL_COMPLEX_FLOAT
#include "../templates_on.h"
TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc) (const size_t n)
{
	TYPE(gsl_block) * block;
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	block = FUNCTION(gsl_block, alloc) (n);

	if (block == 0)
	{
		free(v);

		GSL_ERROR_VAL("failed to allocate space for block",
			GSL_ENOMEM, 0);
	}

	v->data = block->data;
	v->size = n;
	v->stride = 1;
	v->block = block;
	v->owner = 1;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, calloc) (const size_t n)
{
	size_t i;

	TYPE(gsl_vector) * v = FUNCTION(gsl_vector, alloc) (n);

	if (v == 0)
		return 0;

	/* initialize vector to zero */

	for (i = 0; i < MULTIPLICITY * n; i++)
	{
		v->data[i] = 0;
	}

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_from_block) (TYPE(gsl_block) * block,
const size_t offset,
const size_t n,
const size_t stride)
{
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer", GSL_EINVAL, 0);
	}

	if (block->size <= offset + (n - 1) * stride)
	{
		GSL_ERROR_VAL("vector would extend past end of block", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = block->data + MULTIPLICITY * offset;
	v->size = n;
	v->stride = stride;
	v->block = block;
	v->owner = 0;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_from_vector) (TYPE(gsl_vector) * w,
const size_t offset,
const size_t n,
const size_t stride)
{
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer", GSL_EINVAL, 0);
	}

	if (offset + (n - 1) * stride >= w->size)
	{
		GSL_ERROR_VAL("vector would extend past end of block", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = w->data + MULTIPLICITY * w->stride * offset;
	v->size = n;
	v->stride = stride * w->stride;
	v->block = w->block;
	v->owner = 0;

	return v;
}


void
FUNCTION(gsl_vector, free) (TYPE(gsl_vector) * v)
{
	if (v->owner)
	{
		FUNCTION(gsl_block, free) (v->block);
	}
	free(v);
}


void
FUNCTION(gsl_vector, set_all) (TYPE(gsl_vector) * v, BASE x)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;

	size_t i;

	for (i = 0; i < n; i++)
	{
		*(BASE *)(data + MULTIPLICITY * i * stride) = x;
	}
}

void
FUNCTION(gsl_vector, set_zero) (TYPE(gsl_vector) * v)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;
	const BASE zero = ZERO;

	size_t i;

	for (i = 0; i < n; i++)
	{
		*(BASE *)(data + MULTIPLICITY * i * stride) = zero;
	}
}

int
FUNCTION(gsl_vector, set_basis) (TYPE(gsl_vector) * v, size_t i)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;
	const BASE zero = ZERO;
	const BASE one = ONE;

	size_t k;

	if (i >= n)
	{
		GSL_ERROR("index out of range", GSL_EINVAL);
	}

	for (k = 0; k < n; k++)
	{
		*(BASE *)(data + MULTIPLICITY * k * stride) = zero;
	}

	*(BASE *)(data + MULTIPLICITY * i * stride) = one;

	return GSL_SUCCESS;
}
#include "../templates_off.h"
#undef  BASE_GSL_COMPLEX_FLOAT

#define BASE_LONG_DOUBLE
#include "../templates_on.h"
TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc) (const size_t n)
{
	TYPE(gsl_block) * block;
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	block = FUNCTION(gsl_block, alloc) (n);

	if (block == 0)
	{
		free(v);

		GSL_ERROR_VAL("failed to allocate space for block",
			GSL_ENOMEM, 0);
	}

	v->data = block->data;
	v->size = n;
	v->stride = 1;
	v->block = block;
	v->owner = 1;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, calloc) (const size_t n)
{
	size_t i;

	TYPE(gsl_vector) * v = FUNCTION(gsl_vector, alloc) (n);

	if (v == 0)
		return 0;

	/* initialize vector to zero */

	for (i = 0; i < MULTIPLICITY * n; i++)
	{
		v->data[i] = 0;
	}

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_from_block) (TYPE(gsl_block) * block,
const size_t offset,
const size_t n,
const size_t stride)
{
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer", GSL_EINVAL, 0);
	}

	if (block->size <= offset + (n - 1) * stride)
	{
		GSL_ERROR_VAL("vector would extend past end of block", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = block->data + MULTIPLICITY * offset;
	v->size = n;
	v->stride = stride;
	v->block = block;
	v->owner = 0;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_from_vector) (TYPE(gsl_vector) * w,
const size_t offset,
const size_t n,
const size_t stride)
{
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer", GSL_EINVAL, 0);
	}

	if (offset + (n - 1) * stride >= w->size)
	{
		GSL_ERROR_VAL("vector would extend past end of block", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = w->data + MULTIPLICITY * w->stride * offset;
	v->size = n;
	v->stride = stride * w->stride;
	v->block = w->block;
	v->owner = 0;

	return v;
}


void
FUNCTION(gsl_vector, free) (TYPE(gsl_vector) * v)
{
	if (v->owner)
	{
		FUNCTION(gsl_block, free) (v->block);
	}
	free(v);
}

void
FUNCTION(gsl_vector, set_all) (TYPE(gsl_vector) * v, BASE x)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;

	size_t i;

	for (i = 0; i < n; i++)
	{
		*(BASE *)(data + MULTIPLICITY * i * stride) = x;
	}
}

void
FUNCTION(gsl_vector, set_zero) (TYPE(gsl_vector) * v)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;
	const BASE zero = ZERO;

	size_t i;

	for (i = 0; i < n; i++)
	{
		*(BASE *)(data + MULTIPLICITY * i * stride) = zero;
	}
}

int
FUNCTION(gsl_vector, set_basis) (TYPE(gsl_vector) * v, size_t i)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;
	const BASE zero = ZERO;
	const BASE one = ONE;

	size_t k;

	if (i >= n)
	{
		GSL_ERROR("index out of range", GSL_EINVAL);
	}

	for (k = 0; k < n; k++)
	{
		*(BASE *)(data + MULTIPLICITY * k * stride) = zero;
	}

	*(BASE *)(data + MULTIPLICITY * i * stride) = one;

	return GSL_SUCCESS;
}
#include "../templates_off.h"
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
#include "../templates_on.h"
TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc) (const size_t n)
{
	TYPE(gsl_block) * block;
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	block = FUNCTION(gsl_block, alloc) (n);

	if (block == 0)
	{
		free(v);

		GSL_ERROR_VAL("failed to allocate space for block",
			GSL_ENOMEM, 0);
	}

	v->data = block->data;
	v->size = n;
	v->stride = 1;
	v->block = block;
	v->owner = 1;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, calloc) (const size_t n)
{
	size_t i;

	TYPE(gsl_vector) * v = FUNCTION(gsl_vector, alloc) (n);

	if (v == 0)
		return 0;

	/* initialize vector to zero */

	for (i = 0; i < MULTIPLICITY * n; i++)
	{
		v->data[i] = 0;
	}

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_from_block) (TYPE(gsl_block) * block,
const size_t offset,
const size_t n,
const size_t stride)
{
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer", GSL_EINVAL, 0);
	}

	if (block->size <= offset + (n - 1) * stride)
	{
		GSL_ERROR_VAL("vector would extend past end of block", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = block->data + MULTIPLICITY * offset;
	v->size = n;
	v->stride = stride;
	v->block = block;
	v->owner = 0;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_from_vector) (TYPE(gsl_vector) * w,
const size_t offset,
const size_t n,
const size_t stride)
{
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer", GSL_EINVAL, 0);
	}

	if (offset + (n - 1) * stride >= w->size)
	{
		GSL_ERROR_VAL("vector would extend past end of block", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = w->data + MULTIPLICITY * w->stride * offset;
	v->size = n;
	v->stride = stride * w->stride;
	v->block = w->block;
	v->owner = 0;

	return v;
}


void
FUNCTION(gsl_vector, free) (TYPE(gsl_vector) * v)
{
	if (v->owner)
	{
		FUNCTION(gsl_block, free) (v->block);
	}
	free(v);
}


void
FUNCTION(gsl_vector, set_all) (TYPE(gsl_vector) * v, BASE x)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;

	size_t i;

	for (i = 0; i < n; i++)
	{
		*(BASE *)(data + MULTIPLICITY * i * stride) = x;
	}
}

void
FUNCTION(gsl_vector, set_zero) (TYPE(gsl_vector) * v)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;
	const BASE zero = ZERO;

	size_t i;

	for (i = 0; i < n; i++)
	{
		*(BASE *)(data + MULTIPLICITY * i * stride) = zero;
	}
}

int
FUNCTION(gsl_vector, set_basis) (TYPE(gsl_vector) * v, size_t i)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;
	const BASE zero = ZERO;
	const BASE one = ONE;

	size_t k;

	if (i >= n)
	{
		GSL_ERROR("index out of range", GSL_EINVAL);
	}

	for (k = 0; k < n; k++)
	{
		*(BASE *)(data + MULTIPLICITY * k * stride) = zero;
	}

	*(BASE *)(data + MULTIPLICITY * i * stride) = one;

	return GSL_SUCCESS;
}
#include "../templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "../templates_on.h"
TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc) (const size_t n)
{
	TYPE(gsl_block) * block;
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	block = FUNCTION(gsl_block, alloc) (n);

	if (block == 0)
	{
		free(v);

		GSL_ERROR_VAL("failed to allocate space for block",
			GSL_ENOMEM, 0);
	}

	v->data = block->data;
	v->size = n;
	v->stride = 1;
	v->block = block;
	v->owner = 1;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, calloc) (const size_t n)
{
	size_t i;

	TYPE(gsl_vector) * v = FUNCTION(gsl_vector, alloc) (n);

	if (v == 0)
		return 0;

	/* initialize vector to zero */

	for (i = 0; i < MULTIPLICITY * n; i++)
	{
		v->data[i] = 0;
	}

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_from_block) (TYPE(gsl_block) * block,
const size_t offset,
const size_t n,
const size_t stride)
{
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer", GSL_EINVAL, 0);
	}

	if (block->size <= offset + (n - 1) * stride)
	{
		GSL_ERROR_VAL("vector would extend past end of block", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = block->data + MULTIPLICITY * offset;
	v->size = n;
	v->stride = stride;
	v->block = block;
	v->owner = 0;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_from_vector) (TYPE(gsl_vector) * w,
const size_t offset,
const size_t n,
const size_t stride)
{
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer", GSL_EINVAL, 0);
	}

	if (offset + (n - 1) * stride >= w->size)
	{
		GSL_ERROR_VAL("vector would extend past end of block", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = w->data + MULTIPLICITY * w->stride * offset;
	v->size = n;
	v->stride = stride * w->stride;
	v->block = w->block;
	v->owner = 0;

	return v;
}


void
FUNCTION(gsl_vector, free) (TYPE(gsl_vector) * v)
{
	if (v->owner)
	{
		FUNCTION(gsl_block, free) (v->block);
	}
	free(v);
}


void
FUNCTION(gsl_vector, set_all) (TYPE(gsl_vector) * v, BASE x)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;

	size_t i;

	for (i = 0; i < n; i++)
	{
		*(BASE *)(data + MULTIPLICITY * i * stride) = x;
	}
}

void
FUNCTION(gsl_vector, set_zero) (TYPE(gsl_vector) * v)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;
	const BASE zero = ZERO;

	size_t i;

	for (i = 0; i < n; i++)
	{
		*(BASE *)(data + MULTIPLICITY * i * stride) = zero;
	}
}

int
FUNCTION(gsl_vector, set_basis) (TYPE(gsl_vector) * v, size_t i)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;
	const BASE zero = ZERO;
	const BASE one = ONE;

	size_t k;

	if (i >= n)
	{
		GSL_ERROR("index out of range", GSL_EINVAL);
	}

	for (k = 0; k < n; k++)
	{
		*(BASE *)(data + MULTIPLICITY * k * stride) = zero;
	}

	*(BASE *)(data + MULTIPLICITY * i * stride) = one;

	return GSL_SUCCESS;
}
#include "../templates_off.h"
#undef  BASE_FLOAT

#define BASE_ULONG
#include "../templates_on.h"
TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc) (const size_t n)
{
	TYPE(gsl_block) * block;
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	block = FUNCTION(gsl_block, alloc) (n);

	if (block == 0)
	{
		free(v);

		GSL_ERROR_VAL("failed to allocate space for block",
			GSL_ENOMEM, 0);
	}

	v->data = block->data;
	v->size = n;
	v->stride = 1;
	v->block = block;
	v->owner = 1;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, calloc) (const size_t n)
{
	size_t i;

	TYPE(gsl_vector) * v = FUNCTION(gsl_vector, alloc) (n);

	if (v == 0)
		return 0;

	/* initialize vector to zero */

	for (i = 0; i < MULTIPLICITY * n; i++)
	{
		v->data[i] = 0;
	}

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_from_block) (TYPE(gsl_block) * block,
const size_t offset,
const size_t n,
const size_t stride)
{
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer", GSL_EINVAL, 0);
	}

	if (block->size <= offset + (n - 1) * stride)
	{
		GSL_ERROR_VAL("vector would extend past end of block", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = block->data + MULTIPLICITY * offset;
	v->size = n;
	v->stride = stride;
	v->block = block;
	v->owner = 0;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_from_vector) (TYPE(gsl_vector) * w,
const size_t offset,
const size_t n,
const size_t stride)
{
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer", GSL_EINVAL, 0);
	}

	if (offset + (n - 1) * stride >= w->size)
	{
		GSL_ERROR_VAL("vector would extend past end of block", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = w->data + MULTIPLICITY * w->stride * offset;
	v->size = n;
	v->stride = stride * w->stride;
	v->block = w->block;
	v->owner = 0;

	return v;
}


void
FUNCTION(gsl_vector, free) (TYPE(gsl_vector) * v)
{
	if (v->owner)
	{
		FUNCTION(gsl_block, free) (v->block);
	}
	free(v);
}


void
FUNCTION(gsl_vector, set_all) (TYPE(gsl_vector) * v, BASE x)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;

	size_t i;

	for (i = 0; i < n; i++)
	{
		*(BASE *)(data + MULTIPLICITY * i * stride) = x;
	}
}

void
FUNCTION(gsl_vector, set_zero) (TYPE(gsl_vector) * v)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;
	const BASE zero = ZERO;

	size_t i;

	for (i = 0; i < n; i++)
	{
		*(BASE *)(data + MULTIPLICITY * i * stride) = zero;
	}
}

int
FUNCTION(gsl_vector, set_basis) (TYPE(gsl_vector) * v, size_t i)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;
	const BASE zero = ZERO;
	const BASE one = ONE;

	size_t k;

	if (i >= n)
	{
		GSL_ERROR("index out of range", GSL_EINVAL);
	}

	for (k = 0; k < n; k++)
	{
		*(BASE *)(data + MULTIPLICITY * k * stride) = zero;
	}

	*(BASE *)(data + MULTIPLICITY * i * stride) = one;

	return GSL_SUCCESS;
}
#include "../templates_off.h"
#undef  BASE_ULONG

#define BASE_LONG
#include "../templates_on.h"
TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc) (const size_t n)
{
	TYPE(gsl_block) * block;
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	block = FUNCTION(gsl_block, alloc) (n);

	if (block == 0)
	{
		free(v);

		GSL_ERROR_VAL("failed to allocate space for block",
			GSL_ENOMEM, 0);
	}

	v->data = block->data;
	v->size = n;
	v->stride = 1;
	v->block = block;
	v->owner = 1;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, calloc) (const size_t n)
{
	size_t i;

	TYPE(gsl_vector) * v = FUNCTION(gsl_vector, alloc) (n);

	if (v == 0)
		return 0;

	/* initialize vector to zero */

	for (i = 0; i < MULTIPLICITY * n; i++)
	{
		v->data[i] = 0;
	}

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_from_block) (TYPE(gsl_block) * block,
const size_t offset,
const size_t n,
const size_t stride)
{
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer", GSL_EINVAL, 0);
	}

	if (block->size <= offset + (n - 1) * stride)
	{
		GSL_ERROR_VAL("vector would extend past end of block", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = block->data + MULTIPLICITY * offset;
	v->size = n;
	v->stride = stride;
	v->block = block;
	v->owner = 0;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_from_vector) (TYPE(gsl_vector) * w,
const size_t offset,
const size_t n,
const size_t stride)
{
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer", GSL_EINVAL, 0);
	}

	if (offset + (n - 1) * stride >= w->size)
	{
		GSL_ERROR_VAL("vector would extend past end of block", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = w->data + MULTIPLICITY * w->stride * offset;
	v->size = n;
	v->stride = stride * w->stride;
	v->block = w->block;
	v->owner = 0;

	return v;
}


void
FUNCTION(gsl_vector, free) (TYPE(gsl_vector) * v)
{
	if (v->owner)
	{
		FUNCTION(gsl_block, free) (v->block);
	}
	free(v);
}


void
FUNCTION(gsl_vector, set_all) (TYPE(gsl_vector) * v, BASE x)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;

	size_t i;

	for (i = 0; i < n; i++)
	{
		*(BASE *)(data + MULTIPLICITY * i * stride) = x;
	}
}

void
FUNCTION(gsl_vector, set_zero) (TYPE(gsl_vector) * v)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;
	const BASE zero = ZERO;

	size_t i;

	for (i = 0; i < n; i++)
	{
		*(BASE *)(data + MULTIPLICITY * i * stride) = zero;
	}
}

int
FUNCTION(gsl_vector, set_basis) (TYPE(gsl_vector) * v, size_t i)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;
	const BASE zero = ZERO;
	const BASE one = ONE;

	size_t k;

	if (i >= n)
	{
		GSL_ERROR("index out of range", GSL_EINVAL);
	}

	for (k = 0; k < n; k++)
	{
		*(BASE *)(data + MULTIPLICITY * k * stride) = zero;
	}

	*(BASE *)(data + MULTIPLICITY * i * stride) = one;

	return GSL_SUCCESS;
}
#include "../templates_off.h"
#undef  BASE_LONG

#define BASE_UINT
#include "../templates_on.h"
TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc) (const size_t n)
{
	TYPE(gsl_block) * block;
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	block = FUNCTION(gsl_block, alloc) (n);

	if (block == 0)
	{
		free(v);

		GSL_ERROR_VAL("failed to allocate space for block",
			GSL_ENOMEM, 0);
	}

	v->data = block->data;
	v->size = n;
	v->stride = 1;
	v->block = block;
	v->owner = 1;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, calloc) (const size_t n)
{
	size_t i;

	TYPE(gsl_vector) * v = FUNCTION(gsl_vector, alloc) (n);

	if (v == 0)
		return 0;

	/* initialize vector to zero */

	for (i = 0; i < MULTIPLICITY * n; i++)
	{
		v->data[i] = 0;
	}

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_from_block) (TYPE(gsl_block) * block,
const size_t offset,
const size_t n,
const size_t stride)
{
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer", GSL_EINVAL, 0);
	}

	if (block->size <= offset + (n - 1) * stride)
	{
		GSL_ERROR_VAL("vector would extend past end of block", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = block->data + MULTIPLICITY * offset;
	v->size = n;
	v->stride = stride;
	v->block = block;
	v->owner = 0;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_from_vector) (TYPE(gsl_vector) * w,
const size_t offset,
const size_t n,
const size_t stride)
{
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer", GSL_EINVAL, 0);
	}

	if (offset + (n - 1) * stride >= w->size)
	{
		GSL_ERROR_VAL("vector would extend past end of block", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = w->data + MULTIPLICITY * w->stride * offset;
	v->size = n;
	v->stride = stride * w->stride;
	v->block = w->block;
	v->owner = 0;

	return v;
}


void
FUNCTION(gsl_vector, free) (TYPE(gsl_vector) * v)
{
	if (v->owner)
	{
		FUNCTION(gsl_block, free) (v->block);
	}
	free(v);
}


void
FUNCTION(gsl_vector, set_all) (TYPE(gsl_vector) * v, BASE x)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;

	size_t i;

	for (i = 0; i < n; i++)
	{
		*(BASE *)(data + MULTIPLICITY * i * stride) = x;
	}
}

void
FUNCTION(gsl_vector, set_zero) (TYPE(gsl_vector) * v)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;
	const BASE zero = ZERO;

	size_t i;

	for (i = 0; i < n; i++)
	{
		*(BASE *)(data + MULTIPLICITY * i * stride) = zero;
	}
}

int
FUNCTION(gsl_vector, set_basis) (TYPE(gsl_vector) * v, size_t i)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;
	const BASE zero = ZERO;
	const BASE one = ONE;

	size_t k;

	if (i >= n)
	{
		GSL_ERROR("index out of range", GSL_EINVAL);
	}

	for (k = 0; k < n; k++)
	{
		*(BASE *)(data + MULTIPLICITY * k * stride) = zero;
	}

	*(BASE *)(data + MULTIPLICITY * i * stride) = one;

	return GSL_SUCCESS;
}
#include "../templates_off.h"
#undef  BASE_UINT

#define BASE_INT
#include "../templates_on.h"
TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc) (const size_t n)
{
	TYPE(gsl_block) * block;
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	block = FUNCTION(gsl_block, alloc) (n);

	if (block == 0)
	{
		free(v);

		GSL_ERROR_VAL("failed to allocate space for block",
			GSL_ENOMEM, 0);
	}

	v->data = block->data;
	v->size = n;
	v->stride = 1;
	v->block = block;
	v->owner = 1;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, calloc) (const size_t n)
{
	size_t i;

	TYPE(gsl_vector) * v = FUNCTION(gsl_vector, alloc) (n);

	if (v == 0)
		return 0;

	/* initialize vector to zero */

	for (i = 0; i < MULTIPLICITY * n; i++)
	{
		v->data[i] = 0;
	}

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_from_block) (TYPE(gsl_block) * block,
const size_t offset,
const size_t n,
const size_t stride)
{
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer", GSL_EINVAL, 0);
	}

	if (block->size <= offset + (n - 1) * stride)
	{
		GSL_ERROR_VAL("vector would extend past end of block", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = block->data + MULTIPLICITY * offset;
	v->size = n;
	v->stride = stride;
	v->block = block;
	v->owner = 0;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_from_vector) (TYPE(gsl_vector) * w,
const size_t offset,
const size_t n,
const size_t stride)
{
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer", GSL_EINVAL, 0);
	}

	if (offset + (n - 1) * stride >= w->size)
	{
		GSL_ERROR_VAL("vector would extend past end of block", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = w->data + MULTIPLICITY * w->stride * offset;
	v->size = n;
	v->stride = stride * w->stride;
	v->block = w->block;
	v->owner = 0;

	return v;
}


void
FUNCTION(gsl_vector, free) (TYPE(gsl_vector) * v)
{
	if (v->owner)
	{
		FUNCTION(gsl_block, free) (v->block);
	}
	free(v);
}


void
FUNCTION(gsl_vector, set_all) (TYPE(gsl_vector) * v, BASE x)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;

	size_t i;

	for (i = 0; i < n; i++)
	{
		*(BASE *)(data + MULTIPLICITY * i * stride) = x;
	}
}

void
FUNCTION(gsl_vector, set_zero) (TYPE(gsl_vector) * v)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;
	const BASE zero = ZERO;

	size_t i;

	for (i = 0; i < n; i++)
	{
		*(BASE *)(data + MULTIPLICITY * i * stride) = zero;
	}
}

int
FUNCTION(gsl_vector, set_basis) (TYPE(gsl_vector) * v, size_t i)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;
	const BASE zero = ZERO;
	const BASE one = ONE;

	size_t k;

	if (i >= n)
	{
		GSL_ERROR("index out of range", GSL_EINVAL);
	}

	for (k = 0; k < n; k++)
	{
		*(BASE *)(data + MULTIPLICITY * k * stride) = zero;
	}

	*(BASE *)(data + MULTIPLICITY * i * stride) = one;

	return GSL_SUCCESS;
}
#include "../templates_off.h"
#undef  BASE_INT

#define BASE_USHORT
#include "../templates_on.h"
TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc) (const size_t n)
{
	TYPE(gsl_block) * block;
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	block = FUNCTION(gsl_block, alloc) (n);

	if (block == 0)
	{
		free(v);

		GSL_ERROR_VAL("failed to allocate space for block",
			GSL_ENOMEM, 0);
	}

	v->data = block->data;
	v->size = n;
	v->stride = 1;
	v->block = block;
	v->owner = 1;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, calloc) (const size_t n)
{
	size_t i;

	TYPE(gsl_vector) * v = FUNCTION(gsl_vector, alloc) (n);

	if (v == 0)
		return 0;

	/* initialize vector to zero */

	for (i = 0; i < MULTIPLICITY * n; i++)
	{
		v->data[i] = 0;
	}

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_from_block) (TYPE(gsl_block) * block,
const size_t offset,
const size_t n,
const size_t stride)
{
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer", GSL_EINVAL, 0);
	}

	if (block->size <= offset + (n - 1) * stride)
	{
		GSL_ERROR_VAL("vector would extend past end of block", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = block->data + MULTIPLICITY * offset;
	v->size = n;
	v->stride = stride;
	v->block = block;
	v->owner = 0;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_from_vector) (TYPE(gsl_vector) * w,
const size_t offset,
const size_t n,
const size_t stride)
{
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer", GSL_EINVAL, 0);
	}

	if (offset + (n - 1) * stride >= w->size)
	{
		GSL_ERROR_VAL("vector would extend past end of block", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = w->data + MULTIPLICITY * w->stride * offset;
	v->size = n;
	v->stride = stride * w->stride;
	v->block = w->block;
	v->owner = 0;

	return v;
}


void
FUNCTION(gsl_vector, free) (TYPE(gsl_vector) * v)
{
	if (v->owner)
	{
		FUNCTION(gsl_block, free) (v->block);
	}
	free(v);
}


void
FUNCTION(gsl_vector, set_all) (TYPE(gsl_vector) * v, BASE x)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;

	size_t i;

	for (i = 0; i < n; i++)
	{
		*(BASE *)(data + MULTIPLICITY * i * stride) = x;
	}
}

void
FUNCTION(gsl_vector, set_zero) (TYPE(gsl_vector) * v)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;
	const BASE zero = ZERO;

	size_t i;

	for (i = 0; i < n; i++)
	{
		*(BASE *)(data + MULTIPLICITY * i * stride) = zero;
	}
}

int
FUNCTION(gsl_vector, set_basis) (TYPE(gsl_vector) * v, size_t i)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;
	const BASE zero = ZERO;
	const BASE one = ONE;

	size_t k;

	if (i >= n)
	{
		GSL_ERROR("index out of range", GSL_EINVAL);
	}

	for (k = 0; k < n; k++)
	{
		*(BASE *)(data + MULTIPLICITY * k * stride) = zero;
	}

	*(BASE *)(data + MULTIPLICITY * i * stride) = one;

	return GSL_SUCCESS;
}
#include "../templates_off.h"
#undef  BASE_USHORT

#define BASE_SHORT
#include "../templates_on.h"
TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc) (const size_t n)
{
	TYPE(gsl_block) * block;
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	block = FUNCTION(gsl_block, alloc) (n);

	if (block == 0)
	{
		free(v);

		GSL_ERROR_VAL("failed to allocate space for block",
			GSL_ENOMEM, 0);
	}

	v->data = block->data;
	v->size = n;
	v->stride = 1;
	v->block = block;
	v->owner = 1;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, calloc) (const size_t n)
{
	size_t i;

	TYPE(gsl_vector) * v = FUNCTION(gsl_vector, alloc) (n);

	if (v == 0)
		return 0;

	/* initialize vector to zero */

	for (i = 0; i < MULTIPLICITY * n; i++)
	{
		v->data[i] = 0;
	}

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_from_block) (TYPE(gsl_block) * block,
const size_t offset,
const size_t n,
const size_t stride)
{
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer", GSL_EINVAL, 0);
	}

	if (block->size <= offset + (n - 1) * stride)
	{
		GSL_ERROR_VAL("vector would extend past end of block", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = block->data + MULTIPLICITY * offset;
	v->size = n;
	v->stride = stride;
	v->block = block;
	v->owner = 0;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_from_vector) (TYPE(gsl_vector) * w,
const size_t offset,
const size_t n,
const size_t stride)
{
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer", GSL_EINVAL, 0);
	}

	if (offset + (n - 1) * stride >= w->size)
	{
		GSL_ERROR_VAL("vector would extend past end of block", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = w->data + MULTIPLICITY * w->stride * offset;
	v->size = n;
	v->stride = stride * w->stride;
	v->block = w->block;
	v->owner = 0;

	return v;
}


void
FUNCTION(gsl_vector, free) (TYPE(gsl_vector) * v)
{
	if (v->owner)
	{
		FUNCTION(gsl_block, free) (v->block);
	}
	free(v);
}


void
FUNCTION(gsl_vector, set_all) (TYPE(gsl_vector) * v, BASE x)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;

	size_t i;

	for (i = 0; i < n; i++)
	{
		*(BASE *)(data + MULTIPLICITY * i * stride) = x;
	}
}

void
FUNCTION(gsl_vector, set_zero) (TYPE(gsl_vector) * v)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;
	const BASE zero = ZERO;

	size_t i;

	for (i = 0; i < n; i++)
	{
		*(BASE *)(data + MULTIPLICITY * i * stride) = zero;
	}
}

int
FUNCTION(gsl_vector, set_basis) (TYPE(gsl_vector) * v, size_t i)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;
	const BASE zero = ZERO;
	const BASE one = ONE;

	size_t k;

	if (i >= n)
	{
		GSL_ERROR("index out of range", GSL_EINVAL);
	}

	for (k = 0; k < n; k++)
	{
		*(BASE *)(data + MULTIPLICITY * k * stride) = zero;
	}

	*(BASE *)(data + MULTIPLICITY * i * stride) = one;

	return GSL_SUCCESS;
}
#include "../templates_off.h"
#undef  BASE_SHORT

#define BASE_UCHAR
#include "../templates_on.h"
TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc) (const size_t n)
{
	TYPE(gsl_block) * block;
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	block = FUNCTION(gsl_block, alloc) (n);

	if (block == 0)
	{
		free(v);

		GSL_ERROR_VAL("failed to allocate space for block",
			GSL_ENOMEM, 0);
	}

	v->data = block->data;
	v->size = n;
	v->stride = 1;
	v->block = block;
	v->owner = 1;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, calloc) (const size_t n)
{
	size_t i;

	TYPE(gsl_vector) * v = FUNCTION(gsl_vector, alloc) (n);

	if (v == 0)
		return 0;

	/* initialize vector to zero */

	for (i = 0; i < MULTIPLICITY * n; i++)
	{
		v->data[i] = 0;
	}

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_from_block) (TYPE(gsl_block) * block,
const size_t offset,
const size_t n,
const size_t stride)
{
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer", GSL_EINVAL, 0);
	}

	if (block->size <= offset + (n - 1) * stride)
	{
		GSL_ERROR_VAL("vector would extend past end of block", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = block->data + MULTIPLICITY * offset;
	v->size = n;
	v->stride = stride;
	v->block = block;
	v->owner = 0;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_from_vector) (TYPE(gsl_vector) * w,
const size_t offset,
const size_t n,
const size_t stride)
{
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer", GSL_EINVAL, 0);
	}

	if (offset + (n - 1) * stride >= w->size)
	{
		GSL_ERROR_VAL("vector would extend past end of block", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = w->data + MULTIPLICITY * w->stride * offset;
	v->size = n;
	v->stride = stride * w->stride;
	v->block = w->block;
	v->owner = 0;

	return v;
}


void
FUNCTION(gsl_vector, free) (TYPE(gsl_vector) * v)
{
	if (v->owner)
	{
		FUNCTION(gsl_block, free) (v->block);
	}
	free(v);
}


void
FUNCTION(gsl_vector, set_all) (TYPE(gsl_vector) * v, BASE x)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;

	size_t i;

	for (i = 0; i < n; i++)
	{
		*(BASE *)(data + MULTIPLICITY * i * stride) = x;
	}
}

void
FUNCTION(gsl_vector, set_zero) (TYPE(gsl_vector) * v)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;
	const BASE zero = ZERO;

	size_t i;

	for (i = 0; i < n; i++)
	{
		*(BASE *)(data + MULTIPLICITY * i * stride) = zero;
	}
}

int
FUNCTION(gsl_vector, set_basis) (TYPE(gsl_vector) * v, size_t i)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;
	const BASE zero = ZERO;
	const BASE one = ONE;

	size_t k;

	if (i >= n)
	{
		GSL_ERROR("index out of range", GSL_EINVAL);
	}

	for (k = 0; k < n; k++)
	{
		*(BASE *)(data + MULTIPLICITY * k * stride) = zero;
	}

	*(BASE *)(data + MULTIPLICITY * i * stride) = one;

	return GSL_SUCCESS;
}
#include "../templates_off.h"
#undef  BASE_UCHAR

#define BASE_CHAR
#include "../templates_on.h"
TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc) (const size_t n)
{
	TYPE(gsl_block) * block;
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	block = FUNCTION(gsl_block, alloc) (n);

	if (block == 0)
	{
		free(v);

		GSL_ERROR_VAL("failed to allocate space for block",
			GSL_ENOMEM, 0);
	}

	v->data = block->data;
	v->size = n;
	v->stride = 1;
	v->block = block;
	v->owner = 1;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, calloc) (const size_t n)
{
	size_t i;

	TYPE(gsl_vector) * v = FUNCTION(gsl_vector, alloc) (n);

	if (v == 0)
		return 0;

	/* initialize vector to zero */

	for (i = 0; i < MULTIPLICITY * n; i++)
	{
		v->data[i] = 0;
	}

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_from_block) (TYPE(gsl_block) * block,
const size_t offset,
const size_t n,
const size_t stride)
{
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer", GSL_EINVAL, 0);
	}

	if (block->size <= offset + (n - 1) * stride)
	{
		GSL_ERROR_VAL("vector would extend past end of block", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = block->data + MULTIPLICITY * offset;
	v->size = n;
	v->stride = stride;
	v->block = block;
	v->owner = 0;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_from_vector) (TYPE(gsl_vector) * w,
const size_t offset,
const size_t n,
const size_t stride)
{
	TYPE(gsl_vector) * v;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, 0);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer", GSL_EINVAL, 0);
	}

	if (offset + (n - 1) * stride >= w->size)
	{
		GSL_ERROR_VAL("vector would extend past end of block", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = w->data + MULTIPLICITY * w->stride * offset;
	v->size = n;
	v->stride = stride * w->stride;
	v->block = w->block;
	v->owner = 0;

	return v;
}


void
FUNCTION(gsl_vector, free) (TYPE(gsl_vector) * v)
{
	if (v->owner)
	{
		FUNCTION(gsl_block, free) (v->block);
	}
	free(v);
}


void
FUNCTION(gsl_vector, set_all) (TYPE(gsl_vector) * v, BASE x)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;

	size_t i;

	for (i = 0; i < n; i++)
	{
		*(BASE *)(data + MULTIPLICITY * i * stride) = x;
	}
}

void
FUNCTION(gsl_vector, set_zero) (TYPE(gsl_vector) * v)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;
	const BASE zero = ZERO;

	size_t i;

	for (i = 0; i < n; i++)
	{
		*(BASE *)(data + MULTIPLICITY * i * stride) = zero;
	}
}

int
FUNCTION(gsl_vector, set_basis) (TYPE(gsl_vector) * v, size_t i)
{
	ATOMIC * const data = v->data;
	const size_t n = v->size;
	const size_t stride = v->stride;
	const BASE zero = ZERO;
	const BASE one = ONE;

	size_t k;

	if (i >= n)
	{
		GSL_ERROR("index out of range", GSL_EINVAL);
	}

	for (k = 0; k < n; k++)
	{
		*(BASE *)(data + MULTIPLICITY * k * stride) = zero;
	}

	*(BASE *)(data + MULTIPLICITY * i * stride) = one;

	return GSL_SUCCESS;
}
#include "../templates_off.h"
#undef  BASE_CHAR
