#include "stdafx.h"
#include <config.h.in>
#include <stdlib.h>
#include <gsl/block/gsl_block.h>

#define BASE_GSL_COMPLEX_LONG
#include "templates_on.h"
TYPE(gsl_block) *
FUNCTION(gsl_block, alloc) (const size_t n)
{
	TYPE(gsl_block) * b;

	if (n == 0)
	{
		GSL_ERROR_VAL("block length n must be positive integer",
			GSL_EINVAL, 0);
	}
	
	b = (TYPE(gsl_block) *) malloc(sizeof(TYPE(gsl_block)));

	if (b == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for block struct",
			GSL_ENOMEM, 0);
	}

	b->data = (ATOMIC *)malloc(MULTIPLICITY * n * sizeof(ATOMIC));

	if (b->data == 0)
	{
		free(b);         /* exception in constructor, avoid memory leak */

		GSL_ERROR_VAL("failed to allocate space for block data",
			GSL_ENOMEM, 0);
	}

	b->size = n;

	return b;
}

TYPE(gsl_block) *
FUNCTION(gsl_block, calloc) (const size_t n)
{
	size_t i;

	TYPE(gsl_block) * b = FUNCTION(gsl_block, alloc) (n);

	if (b == 0)
		return 0;

	/* initialize block to zero */

	for (i = 0; i < MULTIPLICITY * n; i++)
	{
		b->data[i] = 0;
	}

	return b;
}

void
FUNCTION(gsl_block, free) (TYPE(gsl_block) * b)
{
	free(b->data);
	free(b);
}

#include "templates_off.h"
#undef  BASE_GSL_COMPLEX_LONG

#define BASE_GSL_COMPLEX
#include "templates_on.h"
TYPE(gsl_block) *
FUNCTION(gsl_block, alloc) (const size_t n)
{
	TYPE(gsl_block) * b;

	if (n == 0)
	{
		GSL_ERROR_VAL("block length n must be positive integer",
			GSL_EINVAL, 0);
	}

	b = (TYPE(gsl_block) *) malloc(sizeof(TYPE(gsl_block)));

	if (b == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for block struct",
			GSL_ENOMEM, 0);
	}

	b->data = (ATOMIC *)malloc(MULTIPLICITY * n * sizeof(ATOMIC));

	if (b->data == 0)
	{
		free(b);         /* exception in constructor, avoid memory leak */

		GSL_ERROR_VAL("failed to allocate space for block data",
			GSL_ENOMEM, 0);
	}

	b->size = n;

	return b;
}

TYPE(gsl_block) *
FUNCTION(gsl_block, calloc) (const size_t n)
{
	size_t i;

	TYPE(gsl_block) * b = FUNCTION(gsl_block, alloc) (n);

	if (b == 0)
		return 0;

	/* initialize block to zero */

	for (i = 0; i < MULTIPLICITY * n; i++)
	{
		b->data[i] = 0;
	}

	return b;
}

void
FUNCTION(gsl_block, free) (TYPE(gsl_block) * b)
{
	free(b->data);
	free(b);
}

#include "templates_off.h"
#undef  BASE_GSL_COMPLEX

#define BASE_GSL_COMPLEX_FLOAT
#include "templates_on.h"
TYPE(gsl_block) *
FUNCTION(gsl_block, alloc) (const size_t n)
{
	TYPE(gsl_block) * b;

	if (n == 0)
	{
		GSL_ERROR_VAL("block length n must be positive integer",
			GSL_EINVAL, 0);
	}

	b = (TYPE(gsl_block) *) malloc(sizeof(TYPE(gsl_block)));

	if (b == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for block struct",
			GSL_ENOMEM, 0);
	}

	b->data = (ATOMIC *)malloc(MULTIPLICITY * n * sizeof(ATOMIC));

	if (b->data == 0)
	{
		free(b);         /* exception in constructor, avoid memory leak */

		GSL_ERROR_VAL("failed to allocate space for block data",
			GSL_ENOMEM, 0);
	}

	b->size = n;

	return b;
}

TYPE(gsl_block) *
FUNCTION(gsl_block, calloc) (const size_t n)
{
	size_t i;

	TYPE(gsl_block) * b = FUNCTION(gsl_block, alloc) (n);

	if (b == 0)
		return 0;

	/* initialize block to zero */

	for (i = 0; i < MULTIPLICITY * n; i++)
	{
		b->data[i] = 0;
	}

	return b;
}

void
FUNCTION(gsl_block, free) (TYPE(gsl_block) * b)
{
	free(b->data);
	free(b);
}

#include "templates_off.h"
#undef  BASE_GSL_COMPLEX_FLOAT

#define BASE_LONG_DOUBLE
#include "templates_on.h"
TYPE(gsl_block) *
FUNCTION(gsl_block, alloc) (const size_t n)
{
	TYPE(gsl_block) * b;

	if (n == 0)
	{
		GSL_ERROR_VAL("block length n must be positive integer",
			GSL_EINVAL, 0);
	}

	b = (TYPE(gsl_block) *) malloc(sizeof(TYPE(gsl_block)));

	if (b == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for block struct",
			GSL_ENOMEM, 0);
	}

	b->data = (ATOMIC *)malloc(MULTIPLICITY * n * sizeof(ATOMIC));

	if (b->data == 0)
	{
		free(b);         /* exception in constructor, avoid memory leak */

		GSL_ERROR_VAL("failed to allocate space for block data",
			GSL_ENOMEM, 0);
	}

	b->size = n;

	return b;
}

TYPE(gsl_block) *
FUNCTION(gsl_block, calloc) (const size_t n)
{
	size_t i;

	TYPE(gsl_block) * b = FUNCTION(gsl_block, alloc) (n);

	if (b == 0)
		return 0;

	/* initialize block to zero */

	for (i = 0; i < MULTIPLICITY * n; i++)
	{
		b->data[i] = 0;
	}

	return b;
}

void
FUNCTION(gsl_block, free) (TYPE(gsl_block) * b)
{
	free(b->data);
	free(b);
}

#include "templates_off.h"
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
#include "templates_on.h"
TYPE(gsl_block) *
FUNCTION(gsl_block, alloc) (const size_t n)
{
	TYPE(gsl_block) * b;

	if (n == 0)
	{
		GSL_ERROR_VAL("block length n must be positive integer",
			GSL_EINVAL, 0);
	}

	b = (TYPE(gsl_block) *) malloc(sizeof(TYPE(gsl_block)));

	if (b == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for block struct",
			GSL_ENOMEM, 0);
	}

	b->data = (ATOMIC *)malloc(MULTIPLICITY * n * sizeof(ATOMIC));

	if (b->data == 0)
	{
		free(b);         /* exception in constructor, avoid memory leak */

		GSL_ERROR_VAL("failed to allocate space for block data",
			GSL_ENOMEM, 0);
	}

	b->size = n;

	return b;
}

TYPE(gsl_block) *
FUNCTION(gsl_block, calloc) (const size_t n)
{
	size_t i;

	TYPE(gsl_block) * b = FUNCTION(gsl_block, alloc) (n);

	if (b == 0)
		return 0;

	/* initialize block to zero */

	for (i = 0; i < MULTIPLICITY * n; i++)
	{
		b->data[i] = 0;
	}

	return b;
}

void
FUNCTION(gsl_block, free) (TYPE(gsl_block) * b)
{
	free(b->data);
	free(b);
}

#include "templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "templates_on.h"
TYPE(gsl_block) *
FUNCTION(gsl_block, alloc) (const size_t n)
{
	TYPE(gsl_block) * b;

	if (n == 0)
	{
		GSL_ERROR_VAL("block length n must be positive integer",
			GSL_EINVAL, 0);
	}

	b = (TYPE(gsl_block) *) malloc(sizeof(TYPE(gsl_block)));

	if (b == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for block struct",
			GSL_ENOMEM, 0);
	}

	b->data = (ATOMIC *)malloc(MULTIPLICITY * n * sizeof(ATOMIC));

	if (b->data == 0)
	{
		free(b);         /* exception in constructor, avoid memory leak */

		GSL_ERROR_VAL("failed to allocate space for block data",
			GSL_ENOMEM, 0);
	}

	b->size = n;

	return b;
}

TYPE(gsl_block) *
FUNCTION(gsl_block, calloc) (const size_t n)
{
	size_t i;

	TYPE(gsl_block) * b = FUNCTION(gsl_block, alloc) (n);

	if (b == 0)
		return 0;

	/* initialize block to zero */

	for (i = 0; i < MULTIPLICITY * n; i++)
	{
		b->data[i] = 0;
	}

	return b;
}

void
FUNCTION(gsl_block, free) (TYPE(gsl_block) * b)
{
	free(b->data);
	free(b);
}

#include "templates_off.h"
#undef  BASE_FLOAT

#define BASE_ULONG
#include "templates_on.h"
TYPE(gsl_block) *
FUNCTION(gsl_block, alloc) (const size_t n)
{
	TYPE(gsl_block) * b;

	if (n == 0)
	{
		GSL_ERROR_VAL("block length n must be positive integer",
			GSL_EINVAL, 0);
	}

	b = (TYPE(gsl_block) *) malloc(sizeof(TYPE(gsl_block)));

	if (b == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for block struct",
			GSL_ENOMEM, 0);
	}

	b->data = (ATOMIC *)malloc(MULTIPLICITY * n * sizeof(ATOMIC));

	if (b->data == 0)
	{
		free(b);         /* exception in constructor, avoid memory leak */

		GSL_ERROR_VAL("failed to allocate space for block data",
			GSL_ENOMEM, 0);
	}

	b->size = n;

	return b;
}

TYPE(gsl_block) *
FUNCTION(gsl_block, calloc) (const size_t n)
{
	size_t i;

	TYPE(gsl_block) * b = FUNCTION(gsl_block, alloc) (n);

	if (b == 0)
		return 0;

	/* initialize block to zero */

	for (i = 0; i < MULTIPLICITY * n; i++)
	{
		b->data[i] = 0;
	}

	return b;
}

void
FUNCTION(gsl_block, free) (TYPE(gsl_block) * b)
{
	free(b->data);
	free(b);
}

#include "templates_off.h"
#undef  BASE_ULONG

#define BASE_LONG
#include "templates_on.h"
TYPE(gsl_block) *
FUNCTION(gsl_block, alloc) (const size_t n)
{
	TYPE(gsl_block) * b;

	if (n == 0)
	{
		GSL_ERROR_VAL("block length n must be positive integer",
			GSL_EINVAL, 0);
	}

	b = (TYPE(gsl_block) *) malloc(sizeof(TYPE(gsl_block)));

	if (b == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for block struct",
			GSL_ENOMEM, 0);
	}

	b->data = (ATOMIC *)malloc(MULTIPLICITY * n * sizeof(ATOMIC));

	if (b->data == 0)
	{
		free(b);         /* exception in constructor, avoid memory leak */

		GSL_ERROR_VAL("failed to allocate space for block data",
			GSL_ENOMEM, 0);
	}

	b->size = n;

	return b;
}

TYPE(gsl_block) *
FUNCTION(gsl_block, calloc) (const size_t n)
{
	size_t i;

	TYPE(gsl_block) * b = FUNCTION(gsl_block, alloc) (n);

	if (b == 0)
		return 0;

	/* initialize block to zero */

	for (i = 0; i < MULTIPLICITY * n; i++)
	{
		b->data[i] = 0;
	}

	return b;
}

void
FUNCTION(gsl_block, free) (TYPE(gsl_block) * b)
{
	free(b->data);
	free(b);
}

#include "templates_off.h"
#undef  BASE_LONG

#define BASE_UINT
#include "templates_on.h"
TYPE(gsl_block) *
FUNCTION(gsl_block, alloc) (const size_t n)
{
	TYPE(gsl_block) * b;

	if (n == 0)
	{
		GSL_ERROR_VAL("block length n must be positive integer",
			GSL_EINVAL, 0);
	}

	b = (TYPE(gsl_block) *) malloc(sizeof(TYPE(gsl_block)));

	if (b == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for block struct",
			GSL_ENOMEM, 0);
	}

	b->data = (ATOMIC *)malloc(MULTIPLICITY * n * sizeof(ATOMIC));

	if (b->data == 0)
	{
		free(b);         /* exception in constructor, avoid memory leak */

		GSL_ERROR_VAL("failed to allocate space for block data",
			GSL_ENOMEM, 0);
	}

	b->size = n;

	return b;
}

TYPE(gsl_block) *
FUNCTION(gsl_block, calloc) (const size_t n)
{
	size_t i;

	TYPE(gsl_block) * b = FUNCTION(gsl_block, alloc) (n);

	if (b == 0)
		return 0;

	/* initialize block to zero */

	for (i = 0; i < MULTIPLICITY * n; i++)
	{
		b->data[i] = 0;
	}

	return b;
}

void
FUNCTION(gsl_block, free) (TYPE(gsl_block) * b)
{
	free(b->data);
	free(b);
}

#include "templates_off.h"
#undef  BASE_UINT

#define BASE_INT
#include "templates_on.h"
TYPE(gsl_block) *
FUNCTION(gsl_block, alloc) (const size_t n)
{
	TYPE(gsl_block) * b;

	if (n == 0)
	{
		GSL_ERROR_VAL("block length n must be positive integer",
			GSL_EINVAL, 0);
	}

	b = (TYPE(gsl_block) *) malloc(sizeof(TYPE(gsl_block)));

	if (b == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for block struct",
			GSL_ENOMEM, 0);
	}

	b->data = (ATOMIC *)malloc(MULTIPLICITY * n * sizeof(ATOMIC));

	if (b->data == 0)
	{
		free(b);         /* exception in constructor, avoid memory leak */

		GSL_ERROR_VAL("failed to allocate space for block data",
			GSL_ENOMEM, 0);
	}

	b->size = n;

	return b;
}

TYPE(gsl_block) *
FUNCTION(gsl_block, calloc) (const size_t n)
{
	size_t i;

	TYPE(gsl_block) * b = FUNCTION(gsl_block, alloc) (n);

	if (b == 0)
		return 0;

	/* initialize block to zero */

	for (i = 0; i < MULTIPLICITY * n; i++)
	{
		b->data[i] = 0;
	}

	return b;
}

void
FUNCTION(gsl_block, free) (TYPE(gsl_block) * b)
{
	free(b->data);
	free(b);
}

#include "templates_off.h"
#undef  BASE_INT

#define BASE_USHORT
#include "templates_on.h"
TYPE(gsl_block) *
FUNCTION(gsl_block, alloc) (const size_t n)
{
	TYPE(gsl_block) * b;

	if (n == 0)
	{
		GSL_ERROR_VAL("block length n must be positive integer",
			GSL_EINVAL, 0);
	}

	b = (TYPE(gsl_block) *) malloc(sizeof(TYPE(gsl_block)));

	if (b == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for block struct",
			GSL_ENOMEM, 0);
	}

	b->data = (ATOMIC *)malloc(MULTIPLICITY * n * sizeof(ATOMIC));

	if (b->data == 0)
	{
		free(b);         /* exception in constructor, avoid memory leak */

		GSL_ERROR_VAL("failed to allocate space for block data",
			GSL_ENOMEM, 0);
	}

	b->size = n;

	return b;
}

TYPE(gsl_block) *
FUNCTION(gsl_block, calloc) (const size_t n)
{
	size_t i;

	TYPE(gsl_block) * b = FUNCTION(gsl_block, alloc) (n);

	if (b == 0)
		return 0;

	/* initialize block to zero */

	for (i = 0; i < MULTIPLICITY * n; i++)
	{
		b->data[i] = 0;
	}

	return b;
}

void
FUNCTION(gsl_block, free) (TYPE(gsl_block) * b)
{
	free(b->data);
	free(b);
}

#include "templates_off.h"
#undef  BASE_USHORT

#define BASE_SHORT
#include "templates_on.h"
TYPE(gsl_block) *
FUNCTION(gsl_block, alloc) (const size_t n)
{
	TYPE(gsl_block) * b;

	if (n == 0)
	{
		GSL_ERROR_VAL("block length n must be positive integer",
			GSL_EINVAL, 0);
	}

	b = (TYPE(gsl_block) *) malloc(sizeof(TYPE(gsl_block)));

	if (b == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for block struct",
			GSL_ENOMEM, 0);
	}

	b->data = (ATOMIC *)malloc(MULTIPLICITY * n * sizeof(ATOMIC));

	if (b->data == 0)
	{
		free(b);         /* exception in constructor, avoid memory leak */

		GSL_ERROR_VAL("failed to allocate space for block data",
			GSL_ENOMEM, 0);
	}

	b->size = n;

	return b;
}

TYPE(gsl_block) *
FUNCTION(gsl_block, calloc) (const size_t n)
{
	size_t i;

	TYPE(gsl_block) * b = FUNCTION(gsl_block, alloc) (n);

	if (b == 0)
		return 0;

	/* initialize block to zero */

	for (i = 0; i < MULTIPLICITY * n; i++)
	{
		b->data[i] = 0;
	}

	return b;
}

void
FUNCTION(gsl_block, free) (TYPE(gsl_block) * b)
{
	free(b->data);
	free(b);
}

#include "templates_off.h"
#undef  BASE_SHORT

#define BASE_UCHAR
#include "templates_on.h"
TYPE(gsl_block) *
FUNCTION(gsl_block, alloc) (const size_t n)
{
	TYPE(gsl_block) * b;

	if (n == 0)
	{
		GSL_ERROR_VAL("block length n must be positive integer",
			GSL_EINVAL, 0);
	}

	b = (TYPE(gsl_block) *) malloc(sizeof(TYPE(gsl_block)));

	if (b == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for block struct",
			GSL_ENOMEM, 0);
	}

	b->data = (ATOMIC *)malloc(MULTIPLICITY * n * sizeof(ATOMIC));

	if (b->data == 0)
	{
		free(b);         /* exception in constructor, avoid memory leak */

		GSL_ERROR_VAL("failed to allocate space for block data",
			GSL_ENOMEM, 0);
	}

	b->size = n;

	return b;
}

TYPE(gsl_block) *
FUNCTION(gsl_block, calloc) (const size_t n)
{
	size_t i;

	TYPE(gsl_block) * b = FUNCTION(gsl_block, alloc) (n);

	if (b == 0)
		return 0;

	/* initialize block to zero */

	for (i = 0; i < MULTIPLICITY * n; i++)
	{
		b->data[i] = 0;
	}

	return b;
}

void
FUNCTION(gsl_block, free) (TYPE(gsl_block) * b)
{
	free(b->data);
	free(b);
}

#include "templates_off.h"
#undef  BASE_UCHAR

#define BASE_CHAR
#include "templates_on.h"
TYPE(gsl_block) *
FUNCTION(gsl_block, alloc) (const size_t n)
{
	TYPE(gsl_block) * b;

	if (n == 0)
	{
		GSL_ERROR_VAL("block length n must be positive integer",
			GSL_EINVAL, 0);
	}

	b = (TYPE(gsl_block) *) malloc(sizeof(TYPE(gsl_block)));

	if (b == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for block struct",
			GSL_ENOMEM, 0);
	}

	b->data = (ATOMIC *)malloc(MULTIPLICITY * n * sizeof(ATOMIC));

	if (b->data == 0)
	{
		free(b);         /* exception in constructor, avoid memory leak */

		GSL_ERROR_VAL("failed to allocate space for block data",
			GSL_ENOMEM, 0);
	}

	b->size = n;

	return b;
}

TYPE(gsl_block) *
FUNCTION(gsl_block, calloc) (const size_t n)
{
	size_t i;

	TYPE(gsl_block) * b = FUNCTION(gsl_block, alloc) (n);

	if (b == 0)
		return 0;

	/* initialize block to zero */

	for (i = 0; i < MULTIPLICITY * n; i++)
	{
		b->data[i] = 0;
	}

	return b;
}

void
FUNCTION(gsl_block, free) (TYPE(gsl_block) * b)
{
	free(b->data);
	free(b);
}

#include "templates_off.h"
#undef  BASE_CHAR
