#include "stdafx.h"
#include "config.h.in"
#include <stdlib.h>
#include <gsl/vector/gsl_vector.h>

#include "view.h"

#define BASE_GSL_COMPLEX_LONG
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array) (QUALIFIER ATOMIC * base, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = 1;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array_with_stride) (QUALIFIER ATOMIC * base,
size_t stride,
size_t n)
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

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = stride;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_GSL_COMPLEX_LONG

#define BASE_GSL_COMPLEX
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array) (QUALIFIER ATOMIC * base, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = 1;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array_with_stride) (QUALIFIER ATOMIC * base,
size_t stride,
size_t n)
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

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = stride;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_GSL_COMPLEX

#define BASE_GSL_COMPLEX_FLOAT
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array) (QUALIFIER ATOMIC * base, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = 1;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array_with_stride) (QUALIFIER ATOMIC * base,
size_t stride,
size_t n)
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

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = stride;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_GSL_COMPLEX_FLOAT

#define BASE_LONG_DOUBLE
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array) (QUALIFIER ATOMIC * base, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = 1;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array_with_stride) (QUALIFIER ATOMIC * base,
size_t stride,
size_t n)
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

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = stride;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array) (QUALIFIER ATOMIC * base, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = 1;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array_with_stride) (QUALIFIER ATOMIC * base,
size_t stride,
size_t n)
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

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = stride;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array) (QUALIFIER ATOMIC * base, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = 1;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array_with_stride) (QUALIFIER ATOMIC * base,
size_t stride,
size_t n)
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

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = stride;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_FLOAT

#define BASE_ULONG
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array) (QUALIFIER ATOMIC * base, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = 1;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array_with_stride) (QUALIFIER ATOMIC * base,
size_t stride,
size_t n)
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

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = stride;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_ULONG

#define BASE_LONG
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array) (QUALIFIER ATOMIC * base, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = 1;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array_with_stride) (QUALIFIER ATOMIC * base,
size_t stride,
size_t n)
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

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = stride;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_LONG

#define BASE_UINT
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array) (QUALIFIER ATOMIC * base, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = 1;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array_with_stride) (QUALIFIER ATOMIC * base,
size_t stride,
size_t n)
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

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = stride;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_UINT

#define BASE_INT
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array) (QUALIFIER ATOMIC * base, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = 1;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array_with_stride) (QUALIFIER ATOMIC * base,
size_t stride,
size_t n)
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

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = stride;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_INT

#define BASE_USHORT
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array) (QUALIFIER ATOMIC * base, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = 1;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array_with_stride) (QUALIFIER ATOMIC * base,
size_t stride,
size_t n)
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

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = stride;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_USHORT

#define BASE_SHORT
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array) (QUALIFIER ATOMIC * base, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = 1;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array_with_stride) (QUALIFIER ATOMIC * base,
size_t stride,
size_t n)
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

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = stride;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_SHORT

#define BASE_UCHAR
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array) (QUALIFIER ATOMIC * base, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = 1;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array_with_stride) (QUALIFIER ATOMIC * base,
size_t stride,
size_t n)
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

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = stride;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_UCHAR

#define BASE_CHAR
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array) (QUALIFIER ATOMIC * base, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = 1;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array_with_stride) (QUALIFIER ATOMIC * base,
size_t stride,
size_t n)
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

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = stride;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
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
FUNCTION(gsl_vector, view_array) (QUALIFIER ATOMIC * base, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = 1;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array_with_stride) (QUALIFIER ATOMIC * base,
size_t stride,
size_t n)
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

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = stride;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_GSL_COMPLEX_LONG

#define BASE_GSL_COMPLEX
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array) (QUALIFIER ATOMIC * base, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = 1;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array_with_stride) (QUALIFIER ATOMIC * base,
size_t stride,
size_t n)
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

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = stride;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_GSL_COMPLEX

#define BASE_GSL_COMPLEX_FLOAT
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array) (QUALIFIER ATOMIC * base, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = 1;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array_with_stride) (QUALIFIER ATOMIC * base,
size_t stride,
size_t n)
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

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = stride;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_GSL_COMPLEX_FLOAT

#define BASE_LONG_DOUBLE
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array) (QUALIFIER ATOMIC * base, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = 1;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array_with_stride) (QUALIFIER ATOMIC * base,
size_t stride,
size_t n)
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

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = stride;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array) (QUALIFIER ATOMIC * base, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = 1;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array_with_stride) (QUALIFIER ATOMIC * base,
size_t stride,
size_t n)
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

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = stride;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array) (QUALIFIER ATOMIC * base, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = 1;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array_with_stride) (QUALIFIER ATOMIC * base,
size_t stride,
size_t n)
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

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = stride;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_FLOAT

#define BASE_ULONG
#include "../templates_on.h"
#include "view_source.c"
#include "../templates_off.h"
#undef  BASE_ULONG

#define BASE_LONG
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array) (QUALIFIER ATOMIC * base, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = 1;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array_with_stride) (QUALIFIER ATOMIC * base,
size_t stride,
size_t n)
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

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = stride;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_LONG

#define BASE_UINT
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array) (QUALIFIER ATOMIC * base, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = 1;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array_with_stride) (QUALIFIER ATOMIC * base,
size_t stride,
size_t n)
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

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = stride;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_UINT

#define BASE_INT
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array) (QUALIFIER ATOMIC * base, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = 1;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array_with_stride) (QUALIFIER ATOMIC * base,
size_t stride,
size_t n)
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

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = stride;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_INT

#define BASE_USHORT
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array) (QUALIFIER ATOMIC * base, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = 1;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array_with_stride) (QUALIFIER ATOMIC * base,
size_t stride,
size_t n)
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

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = stride;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_USHORT

#define BASE_SHORT
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array) (QUALIFIER ATOMIC * base, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = 1;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array_with_stride) (QUALIFIER ATOMIC * base,
size_t stride,
size_t n)
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

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = stride;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_SHORT

#define BASE_UCHAR
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array) (QUALIFIER ATOMIC * base, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = 1;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array_with_stride) (QUALIFIER ATOMIC * base,
size_t stride,
size_t n)
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

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = stride;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_UCHAR

#define BASE_CHAR
#include "../templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array) (QUALIFIER ATOMIC * base, size_t n)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (n == 0)
	{
		GSL_ERROR_VAL("vector length n must be positive integer",
			GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = 1;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, view_array_with_stride) (QUALIFIER ATOMIC * base,
size_t stride,
size_t n)
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

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = (ATOMIC *)base;
	  v.size = n;
	  v.stride = stride;
	  v.block = 0;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

#include "../templates_off.h"
#undef  BASE_CHAR
