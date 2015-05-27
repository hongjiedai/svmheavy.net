#include "stdafx.h"
#include "../config.h.in"
#include <gsl/vector/gsl_vector.h>
#include <gsl/err/gsl_errno.h>

#define BASE_GSL_COMPLEX_LONG
#include "../templates_on.h"
int
FUNCTION(gsl_vector, memcpy) (TYPE(gsl_vector) * dest,
const TYPE(gsl_vector) * src)
{
	const size_t src_size = src->size;
	const size_t dest_size = dest->size;

	if (src_size != dest_size)
	{
		GSL_ERROR("vector lengths are not equal", GSL_EBADLEN);
	}

  {
	  const size_t src_stride = src->stride;
	  const size_t dest_stride = dest->stride;
	  size_t j;

	  for (j = 0; j < src_size; j++)
	  {
		  size_t k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  dest->data[MULTIPLICITY * dest_stride * j + k]
				  = src->data[MULTIPLICITY * src_stride * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}
#include "../templates_off.h"
#undef  BASE_GSL_COMPLEX_LONG

#define BASE_GSL_COMPLEX
#include "../templates_on.h"
int
FUNCTION(gsl_vector, memcpy) (TYPE(gsl_vector) * dest,
const TYPE(gsl_vector) * src)
{
	const size_t src_size = src->size;
	const size_t dest_size = dest->size;

	if (src_size != dest_size)
	{
		GSL_ERROR("vector lengths are not equal", GSL_EBADLEN);
	}

  {
	  const size_t src_stride = src->stride;
	  const size_t dest_stride = dest->stride;
	  size_t j;

	  for (j = 0; j < src_size; j++)
	  {
		  size_t k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  dest->data[MULTIPLICITY * dest_stride * j + k]
				  = src->data[MULTIPLICITY * src_stride * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}
#include "../templates_off.h"
#undef  BASE_GSL_COMPLEX

#define BASE_GSL_COMPLEX_FLOAT
#include "../templates_on.h"
int
FUNCTION(gsl_vector, memcpy) (TYPE(gsl_vector) * dest,
const TYPE(gsl_vector) * src)
{
	const size_t src_size = src->size;
	const size_t dest_size = dest->size;

	if (src_size != dest_size)
	{
		GSL_ERROR("vector lengths are not equal", GSL_EBADLEN);
	}

  {
	  const size_t src_stride = src->stride;
	  const size_t dest_stride = dest->stride;
	  size_t j;

	  for (j = 0; j < src_size; j++)
	  {
		  size_t k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  dest->data[MULTIPLICITY * dest_stride * j + k]
				  = src->data[MULTIPLICITY * src_stride * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}
#include "../templates_off.h"
#undef  BASE_GSL_COMPLEX_FLOAT

#define BASE_LONG_DOUBLE
#include "../templates_on.h"
int
FUNCTION(gsl_vector, memcpy) (TYPE(gsl_vector) * dest,
const TYPE(gsl_vector) * src)
{
	const size_t src_size = src->size;
	const size_t dest_size = dest->size;

	if (src_size != dest_size)
	{
		GSL_ERROR("vector lengths are not equal", GSL_EBADLEN);
	}

  {
	  const size_t src_stride = src->stride;
	  const size_t dest_stride = dest->stride;
	  size_t j;

	  for (j = 0; j < src_size; j++)
	  {
		  size_t k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  dest->data[MULTIPLICITY * dest_stride * j + k]
				  = src->data[MULTIPLICITY * src_stride * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}
#include "../templates_off.h"
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
#include "../templates_on.h"
int
FUNCTION(gsl_vector, memcpy) (TYPE(gsl_vector) * dest,
const TYPE(gsl_vector) * src)
{
	const size_t src_size = src->size;
	const size_t dest_size = dest->size;

	if (src_size != dest_size)
	{
		GSL_ERROR("vector lengths are not equal", GSL_EBADLEN);
	}

  {
	  const size_t src_stride = src->stride;
	  const size_t dest_stride = dest->stride;
	  size_t j;

	  for (j = 0; j < src_size; j++)
	  {
		  size_t k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  dest->data[MULTIPLICITY * dest_stride * j + k]
				  = src->data[MULTIPLICITY * src_stride * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}
#include "../templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "../templates_on.h"
int
FUNCTION(gsl_vector, memcpy) (TYPE(gsl_vector) * dest,
const TYPE(gsl_vector) * src)
{
	const size_t src_size = src->size;
	const size_t dest_size = dest->size;

	if (src_size != dest_size)
	{
		GSL_ERROR("vector lengths are not equal", GSL_EBADLEN);
	}

  {
	  const size_t src_stride = src->stride;
	  const size_t dest_stride = dest->stride;
	  size_t j;

	  for (j = 0; j < src_size; j++)
	  {
		  size_t k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  dest->data[MULTIPLICITY * dest_stride * j + k]
				  = src->data[MULTIPLICITY * src_stride * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}
#include "../templates_off.h"
#undef  BASE_FLOAT

#define BASE_ULONG
#include "../templates_on.h"
int
FUNCTION(gsl_vector, memcpy) (TYPE(gsl_vector) * dest,
const TYPE(gsl_vector) * src)
{
	const size_t src_size = src->size;
	const size_t dest_size = dest->size;

	if (src_size != dest_size)
	{
		GSL_ERROR("vector lengths are not equal", GSL_EBADLEN);
	}

  {
	  const size_t src_stride = src->stride;
	  const size_t dest_stride = dest->stride;
	  size_t j;

	  for (j = 0; j < src_size; j++)
	  {
		  size_t k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  dest->data[MULTIPLICITY * dest_stride * j + k]
				  = src->data[MULTIPLICITY * src_stride * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}
#include "../templates_off.h"
#undef  BASE_ULONG

#define BASE_LONG
#include "../templates_on.h"
int
FUNCTION(gsl_vector, memcpy) (TYPE(gsl_vector) * dest,
const TYPE(gsl_vector) * src)
{
	const size_t src_size = src->size;
	const size_t dest_size = dest->size;

	if (src_size != dest_size)
	{
		GSL_ERROR("vector lengths are not equal", GSL_EBADLEN);
	}

  {
	  const size_t src_stride = src->stride;
	  const size_t dest_stride = dest->stride;
	  size_t j;

	  for (j = 0; j < src_size; j++)
	  {
		  size_t k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  dest->data[MULTIPLICITY * dest_stride * j + k]
				  = src->data[MULTIPLICITY * src_stride * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}
#include "../templates_off.h"
#undef  BASE_LONG

#define BASE_UINT
#include "../templates_on.h"
int
FUNCTION(gsl_vector, memcpy) (TYPE(gsl_vector) * dest,
const TYPE(gsl_vector) * src)
{
	const size_t src_size = src->size;
	const size_t dest_size = dest->size;

	if (src_size != dest_size)
	{
		GSL_ERROR("vector lengths are not equal", GSL_EBADLEN);
	}

  {
	  const size_t src_stride = src->stride;
	  const size_t dest_stride = dest->stride;
	  size_t j;

	  for (j = 0; j < src_size; j++)
	  {
		  size_t k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  dest->data[MULTIPLICITY * dest_stride * j + k]
				  = src->data[MULTIPLICITY * src_stride * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}
#include "../templates_off.h"
#undef  BASE_UINT

#define BASE_INT
#include "../templates_on.h"
int
FUNCTION(gsl_vector, memcpy) (TYPE(gsl_vector) * dest,
const TYPE(gsl_vector) * src)
{
	const size_t src_size = src->size;
	const size_t dest_size = dest->size;

	if (src_size != dest_size)
	{
		GSL_ERROR("vector lengths are not equal", GSL_EBADLEN);
	}

  {
	  const size_t src_stride = src->stride;
	  const size_t dest_stride = dest->stride;
	  size_t j;

	  for (j = 0; j < src_size; j++)
	  {
		  size_t k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  dest->data[MULTIPLICITY * dest_stride * j + k]
				  = src->data[MULTIPLICITY * src_stride * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}
#include "../templates_off.h"
#undef  BASE_INT

#define BASE_USHORT
#include "../templates_on.h"
int
FUNCTION(gsl_vector, memcpy) (TYPE(gsl_vector) * dest,
const TYPE(gsl_vector) * src)
{
	const size_t src_size = src->size;
	const size_t dest_size = dest->size;

	if (src_size != dest_size)
	{
		GSL_ERROR("vector lengths are not equal", GSL_EBADLEN);
	}

  {
	  const size_t src_stride = src->stride;
	  const size_t dest_stride = dest->stride;
	  size_t j;

	  for (j = 0; j < src_size; j++)
	  {
		  size_t k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  dest->data[MULTIPLICITY * dest_stride * j + k]
				  = src->data[MULTIPLICITY * src_stride * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}
#include "../templates_off.h"
#undef  BASE_USHORT

#define BASE_SHORT
#include "../templates_on.h"
int
FUNCTION(gsl_vector, memcpy) (TYPE(gsl_vector) * dest,
const TYPE(gsl_vector) * src)
{
	const size_t src_size = src->size;
	const size_t dest_size = dest->size;

	if (src_size != dest_size)
	{
		GSL_ERROR("vector lengths are not equal", GSL_EBADLEN);
	}

  {
	  const size_t src_stride = src->stride;
	  const size_t dest_stride = dest->stride;
	  size_t j;

	  for (j = 0; j < src_size; j++)
	  {
		  size_t k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  dest->data[MULTIPLICITY * dest_stride * j + k]
				  = src->data[MULTIPLICITY * src_stride * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}
#include "../templates_off.h"
#undef  BASE_SHORT

#define BASE_UCHAR
#include "../templates_on.h"
int
FUNCTION(gsl_vector, memcpy) (TYPE(gsl_vector) * dest,
const TYPE(gsl_vector) * src)
{
	const size_t src_size = src->size;
	const size_t dest_size = dest->size;

	if (src_size != dest_size)
	{
		GSL_ERROR("vector lengths are not equal", GSL_EBADLEN);
	}

  {
	  const size_t src_stride = src->stride;
	  const size_t dest_stride = dest->stride;
	  size_t j;

	  for (j = 0; j < src_size; j++)
	  {
		  size_t k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  dest->data[MULTIPLICITY * dest_stride * j + k]
				  = src->data[MULTIPLICITY * src_stride * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}
#include "../templates_off.h"
#undef  BASE_UCHAR

#define BASE_CHAR
#include "../templates_on.h"
int
FUNCTION(gsl_vector, memcpy) (TYPE(gsl_vector) * dest,
const TYPE(gsl_vector) * src)
{
	const size_t src_size = src->size;
	const size_t dest_size = dest->size;

	if (src_size != dest_size)
	{
		GSL_ERROR("vector lengths are not equal", GSL_EBADLEN);
	}

  {
	  const size_t src_stride = src->stride;
	  const size_t dest_stride = dest->stride;
	  size_t j;

	  for (j = 0; j < src_size; j++)
	  {
		  size_t k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  dest->data[MULTIPLICITY * dest_stride * j + k]
				  = src->data[MULTIPLICITY * src_stride * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}
#include "../templates_off.h"
#undef  BASE_CHAR
