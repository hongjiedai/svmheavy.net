#include "stdafx.h"
#include "config.h.in"
#include <stdlib.h>
#include <gsl/vector/gsl_vector.h>

#include "view.h"

#define BASE_GSL_COMPLEX_LONG
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector_with_stride) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t stride, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) * stride >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride * stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_GSL_COMPLEX_LONG

#define BASE_GSL_COMPLEX
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector_with_stride) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t stride, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) * stride >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride * stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_GSL_COMPLEX

#define BASE_GSL_COMPLEX_FLOAT
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector_with_stride) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t stride, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) * stride >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride * stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_GSL_COMPLEX_FLOAT

#define BASE_LONG_DOUBLE
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector_with_stride) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t stride, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) * stride >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride * stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector_with_stride) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t stride, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) * stride >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride * stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector_with_stride) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t stride, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) * stride >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride * stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_FLOAT

#define BASE_ULONG
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector_with_stride) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t stride, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) * stride >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride * stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_ULONG

#define BASE_LONG
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector_with_stride) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t stride, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) * stride >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride * stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_LONG

#define BASE_UINT
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector_with_stride) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t stride, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) * stride >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride * stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_UINT

#define BASE_INT
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector_with_stride) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t stride, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) * stride >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride * stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_INT

#define BASE_USHORT
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector_with_stride) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t stride, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) * stride >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride * stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_USHORT

#define BASE_SHORT
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector_with_stride) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t stride, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) * stride >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride * stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_SHORT

#define BASE_UCHAR
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector_with_stride) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t stride, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) * stride >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride * stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_UCHAR

#define BASE_CHAR
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector_with_stride) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t stride, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) * stride >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride * stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_CHAR

#define USE_QUALIFIER
#define QUALIFIER const

#define BASE_GSL_COMPLEX_LONG
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector_with_stride) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t stride, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) * stride >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride * stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_GSL_COMPLEX_LONG

#define BASE_GSL_COMPLEX
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector_with_stride) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t stride, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) * stride >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride * stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_GSL_COMPLEX

#define BASE_GSL_COMPLEX_FLOAT
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector_with_stride) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t stride, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) * stride >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride * stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_GSL_COMPLEX_FLOAT

#define BASE_LONG_DOUBLE
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector_with_stride) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t stride, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) * stride >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride * stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector_with_stride) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t stride, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) * stride >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride * stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector_with_stride) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t stride, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) * stride >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride * stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_FLOAT

#define BASE_ULONG
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector_with_stride) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t stride, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) * stride >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride * stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_ULONG

#define BASE_LONG
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector_with_stride) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t stride, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) * stride >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride * stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_LONG

#define BASE_UINT
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector_with_stride) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t stride, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) * stride >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride * stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_UINT

#define BASE_INT
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector_with_stride) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t stride, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) * stride >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride * stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_INT

#define BASE_USHORT
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector_with_stride) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t stride, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) * stride >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride * stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_USHORT

#define BASE_SHORT
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector_with_stride) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t stride, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) * stride >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride * stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_SHORT

#define BASE_UCHAR
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector_with_stride) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t stride, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) * stride >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride * stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_UCHAR

#define BASE_CHAR
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, subvector_with_stride) (QUALIFIED_TYPE(gsl_vector) * v, size_t offset, size_t stride, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

	if (stride == 0)
	{
		GSL_ERROR_VAL("stride must be positive integer",
			GSL_EINVAL, view);
	}

	if (offset + (n - 1) * stride >= v->size)
	{
		GSL_ERROR_VAL("view would extend past end of vector",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) s = NULL_VECTOR;

	  s.data = v->data + MULTIPLICITY * v->stride * offset;
	  s.size = n;
	  s.stride = v->stride * stride;
	  s.block = v->block;
	  s.owner = 0;

	  view.vector = s;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_CHAR
