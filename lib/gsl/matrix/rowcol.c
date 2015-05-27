#include "stdafx.h"
#include <config.h.in>
#include <gsl/gsl_math.h>
#include <gsl/matrix/gsl_matrix.h>
#include <gsl/vector/gsl_vector.h>
#include <gsl/err/gsl_errno.h>

#include "view.h"

#define BASE_GSL_COMPLEX_LONG
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, row) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t i)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + i * MULTIPLICITY * m->tda;
	  v.size = m->size2;
	  v.stride = 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, column) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t j)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + j * MULTIPLICITY;
	  v.size = m->size1;
	  v.stride = m->tda;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, diagonal) (QUALIFIED_TYPE(gsl_matrix) * m)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	TYPE(gsl_vector) v = NULL_VECTOR;
	v.data = m->data;
	v.size = GSL_MIN(m->size1, m->size2);
	v.stride = m->tda + 1;
	v.block = m->block;
	v.owner = 0;

	view.vector = v;
	return view;
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, subdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (k >= m->size1)
	{
		GSL_ERROR_VAL("subdiagonal index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY * m->tda;
	  v.size = GSL_MIN(m->size1 - k, m->size2);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, superdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;


	if (k >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY;
	  v.size = GSL_MIN(m->size1, m->size2 - k);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_GSL_COMPLEX_LONG

#define BASE_GSL_COMPLEX
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, row) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t i)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + i * MULTIPLICITY * m->tda;
	  v.size = m->size2;
	  v.stride = 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, column) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t j)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + j * MULTIPLICITY;
	  v.size = m->size1;
	  v.stride = m->tda;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, diagonal) (QUALIFIED_TYPE(gsl_matrix) * m)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	TYPE(gsl_vector) v = NULL_VECTOR;
	v.data = m->data;
	v.size = GSL_MIN(m->size1, m->size2);
	v.stride = m->tda + 1;
	v.block = m->block;
	v.owner = 0;

	view.vector = v;
	return view;
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, subdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (k >= m->size1)
	{
		GSL_ERROR_VAL("subdiagonal index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY * m->tda;
	  v.size = GSL_MIN(m->size1 - k, m->size2);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, superdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;


	if (k >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY;
	  v.size = GSL_MIN(m->size1, m->size2 - k);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_GSL_COMPLEX

#define BASE_GSL_COMPLEX_FLOAT
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, row) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t i)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + i * MULTIPLICITY * m->tda;
	  v.size = m->size2;
	  v.stride = 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, column) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t j)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + j * MULTIPLICITY;
	  v.size = m->size1;
	  v.stride = m->tda;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, diagonal) (QUALIFIED_TYPE(gsl_matrix) * m)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	TYPE(gsl_vector) v = NULL_VECTOR;
	v.data = m->data;
	v.size = GSL_MIN(m->size1, m->size2);
	v.stride = m->tda + 1;
	v.block = m->block;
	v.owner = 0;

	view.vector = v;
	return view;
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, subdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (k >= m->size1)
	{
		GSL_ERROR_VAL("subdiagonal index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY * m->tda;
	  v.size = GSL_MIN(m->size1 - k, m->size2);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, superdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;


	if (k >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY;
	  v.size = GSL_MIN(m->size1, m->size2 - k);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_GSL_COMPLEX_FLOAT

#define BASE_LONG_DOUBLE
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, row) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t i)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + i * MULTIPLICITY * m->tda;
	  v.size = m->size2;
	  v.stride = 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, column) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t j)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + j * MULTIPLICITY;
	  v.size = m->size1;
	  v.stride = m->tda;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, diagonal) (QUALIFIED_TYPE(gsl_matrix) * m)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	TYPE(gsl_vector) v = NULL_VECTOR;
	v.data = m->data;
	v.size = GSL_MIN(m->size1, m->size2);
	v.stride = m->tda + 1;
	v.block = m->block;
	v.owner = 0;

	view.vector = v;
	return view;
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, subdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (k >= m->size1)
	{
		GSL_ERROR_VAL("subdiagonal index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY * m->tda;
	  v.size = GSL_MIN(m->size1 - k, m->size2);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, superdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;


	if (k >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY;
	  v.size = GSL_MIN(m->size1, m->size2 - k);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, row) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t i)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + i * MULTIPLICITY * m->tda;
	  v.size = m->size2;
	  v.stride = 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, column) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t j)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + j * MULTIPLICITY;
	  v.size = m->size1;
	  v.stride = m->tda;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, diagonal) (QUALIFIED_TYPE(gsl_matrix) * m)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	TYPE(gsl_vector) v = NULL_VECTOR;
	v.data = m->data;
	v.size = GSL_MIN(m->size1, m->size2);
	v.stride = m->tda + 1;
	v.block = m->block;
	v.owner = 0;

	view.vector = v;
	return view;
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, subdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (k >= m->size1)
	{
		GSL_ERROR_VAL("subdiagonal index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY * m->tda;
	  v.size = GSL_MIN(m->size1 - k, m->size2);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, superdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;


	if (k >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY;
	  v.size = GSL_MIN(m->size1, m->size2 - k);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, row) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t i)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + i * MULTIPLICITY * m->tda;
	  v.size = m->size2;
	  v.stride = 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, column) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t j)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + j * MULTIPLICITY;
	  v.size = m->size1;
	  v.stride = m->tda;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, diagonal) (QUALIFIED_TYPE(gsl_matrix) * m)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	TYPE(gsl_vector) v = NULL_VECTOR;
	v.data = m->data;
	v.size = GSL_MIN(m->size1, m->size2);
	v.stride = m->tda + 1;
	v.block = m->block;
	v.owner = 0;

	view.vector = v;
	return view;
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, subdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (k >= m->size1)
	{
		GSL_ERROR_VAL("subdiagonal index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY * m->tda;
	  v.size = GSL_MIN(m->size1 - k, m->size2);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, superdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;


	if (k >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY;
	  v.size = GSL_MIN(m->size1, m->size2 - k);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_FLOAT

#define BASE_ULONG
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, row) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t i)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + i * MULTIPLICITY * m->tda;
	  v.size = m->size2;
	  v.stride = 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, column) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t j)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + j * MULTIPLICITY;
	  v.size = m->size1;
	  v.stride = m->tda;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, diagonal) (QUALIFIED_TYPE(gsl_matrix) * m)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	TYPE(gsl_vector) v = NULL_VECTOR;
	v.data = m->data;
	v.size = GSL_MIN(m->size1, m->size2);
	v.stride = m->tda + 1;
	v.block = m->block;
	v.owner = 0;

	view.vector = v;
	return view;
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, subdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (k >= m->size1)
	{
		GSL_ERROR_VAL("subdiagonal index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY * m->tda;
	  v.size = GSL_MIN(m->size1 - k, m->size2);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, superdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;


	if (k >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY;
	  v.size = GSL_MIN(m->size1, m->size2 - k);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_ULONG

#define BASE_LONG
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, row) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t i)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + i * MULTIPLICITY * m->tda;
	  v.size = m->size2;
	  v.stride = 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, column) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t j)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + j * MULTIPLICITY;
	  v.size = m->size1;
	  v.stride = m->tda;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, diagonal) (QUALIFIED_TYPE(gsl_matrix) * m)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	TYPE(gsl_vector) v = NULL_VECTOR;
	v.data = m->data;
	v.size = GSL_MIN(m->size1, m->size2);
	v.stride = m->tda + 1;
	v.block = m->block;
	v.owner = 0;

	view.vector = v;
	return view;
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, subdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (k >= m->size1)
	{
		GSL_ERROR_VAL("subdiagonal index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY * m->tda;
	  v.size = GSL_MIN(m->size1 - k, m->size2);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, superdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;


	if (k >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY;
	  v.size = GSL_MIN(m->size1, m->size2 - k);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_LONG

#define BASE_UINT
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, row) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t i)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + i * MULTIPLICITY * m->tda;
	  v.size = m->size2;
	  v.stride = 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, column) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t j)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + j * MULTIPLICITY;
	  v.size = m->size1;
	  v.stride = m->tda;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, diagonal) (QUALIFIED_TYPE(gsl_matrix) * m)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	TYPE(gsl_vector) v = NULL_VECTOR;
	v.data = m->data;
	v.size = GSL_MIN(m->size1, m->size2);
	v.stride = m->tda + 1;
	v.block = m->block;
	v.owner = 0;

	view.vector = v;
	return view;
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, subdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (k >= m->size1)
	{
		GSL_ERROR_VAL("subdiagonal index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY * m->tda;
	  v.size = GSL_MIN(m->size1 - k, m->size2);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, superdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;


	if (k >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY;
	  v.size = GSL_MIN(m->size1, m->size2 - k);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_UINT

#define BASE_INT
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, row) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t i)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + i * MULTIPLICITY * m->tda;
	  v.size = m->size2;
	  v.stride = 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, column) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t j)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + j * MULTIPLICITY;
	  v.size = m->size1;
	  v.stride = m->tda;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, diagonal) (QUALIFIED_TYPE(gsl_matrix) * m)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	TYPE(gsl_vector) v = NULL_VECTOR;
	v.data = m->data;
	v.size = GSL_MIN(m->size1, m->size2);
	v.stride = m->tda + 1;
	v.block = m->block;
	v.owner = 0;

	view.vector = v;
	return view;
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, subdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (k >= m->size1)
	{
		GSL_ERROR_VAL("subdiagonal index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY * m->tda;
	  v.size = GSL_MIN(m->size1 - k, m->size2);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, superdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;


	if (k >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY;
	  v.size = GSL_MIN(m->size1, m->size2 - k);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_INT

#define BASE_USHORT
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, row) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t i)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + i * MULTIPLICITY * m->tda;
	  v.size = m->size2;
	  v.stride = 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, column) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t j)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + j * MULTIPLICITY;
	  v.size = m->size1;
	  v.stride = m->tda;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, diagonal) (QUALIFIED_TYPE(gsl_matrix) * m)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	TYPE(gsl_vector) v = NULL_VECTOR;
	v.data = m->data;
	v.size = GSL_MIN(m->size1, m->size2);
	v.stride = m->tda + 1;
	v.block = m->block;
	v.owner = 0;

	view.vector = v;
	return view;
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, subdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (k >= m->size1)
	{
		GSL_ERROR_VAL("subdiagonal index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY * m->tda;
	  v.size = GSL_MIN(m->size1 - k, m->size2);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, superdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;


	if (k >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY;
	  v.size = GSL_MIN(m->size1, m->size2 - k);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_USHORT

#define BASE_SHORT
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, row) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t i)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + i * MULTIPLICITY * m->tda;
	  v.size = m->size2;
	  v.stride = 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, column) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t j)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + j * MULTIPLICITY;
	  v.size = m->size1;
	  v.stride = m->tda;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, diagonal) (QUALIFIED_TYPE(gsl_matrix) * m)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	TYPE(gsl_vector) v = NULL_VECTOR;
	v.data = m->data;
	v.size = GSL_MIN(m->size1, m->size2);
	v.stride = m->tda + 1;
	v.block = m->block;
	v.owner = 0;

	view.vector = v;
	return view;
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, subdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (k >= m->size1)
	{
		GSL_ERROR_VAL("subdiagonal index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY * m->tda;
	  v.size = GSL_MIN(m->size1 - k, m->size2);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, superdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;


	if (k >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY;
	  v.size = GSL_MIN(m->size1, m->size2 - k);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_SHORT

#define BASE_UCHAR
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, row) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t i)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + i * MULTIPLICITY * m->tda;
	  v.size = m->size2;
	  v.stride = 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, column) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t j)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + j * MULTIPLICITY;
	  v.size = m->size1;
	  v.stride = m->tda;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, diagonal) (QUALIFIED_TYPE(gsl_matrix) * m)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	TYPE(gsl_vector) v = NULL_VECTOR;
	v.data = m->data;
	v.size = GSL_MIN(m->size1, m->size2);
	v.stride = m->tda + 1;
	v.block = m->block;
	v.owner = 0;

	view.vector = v;
	return view;
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, subdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (k >= m->size1)
	{
		GSL_ERROR_VAL("subdiagonal index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY * m->tda;
	  v.size = GSL_MIN(m->size1 - k, m->size2);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, superdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;


	if (k >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY;
	  v.size = GSL_MIN(m->size1, m->size2 - k);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_UCHAR

#define BASE_CHAR
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, row) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t i)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + i * MULTIPLICITY * m->tda;
	  v.size = m->size2;
	  v.stride = 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, column) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t j)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + j * MULTIPLICITY;
	  v.size = m->size1;
	  v.stride = m->tda;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, diagonal) (QUALIFIED_TYPE(gsl_matrix) * m)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	TYPE(gsl_vector) v = NULL_VECTOR;
	v.data = m->data;
	v.size = GSL_MIN(m->size1, m->size2);
	v.stride = m->tda + 1;
	v.block = m->block;
	v.owner = 0;

	view.vector = v;
	return view;
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, subdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (k >= m->size1)
	{
		GSL_ERROR_VAL("subdiagonal index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY * m->tda;
	  v.size = GSL_MIN(m->size1 - k, m->size2);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, superdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;


	if (k >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY;
	  v.size = GSL_MIN(m->size1, m->size2 - k);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_CHAR

#define USE_QUALIFIER
#define QUALIFIER const

#define BASE_GSL_COMPLEX_LONG
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, row) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t i)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + i * MULTIPLICITY * m->tda;
	  v.size = m->size2;
	  v.stride = 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, column) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t j)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + j * MULTIPLICITY;
	  v.size = m->size1;
	  v.stride = m->tda;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, diagonal) (QUALIFIED_TYPE(gsl_matrix) * m)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	TYPE(gsl_vector) v = NULL_VECTOR;
	v.data = m->data;
	v.size = GSL_MIN(m->size1, m->size2);
	v.stride = m->tda + 1;
	v.block = m->block;
	v.owner = 0;

	view.vector = v;
	return view;
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, subdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (k >= m->size1)
	{
		GSL_ERROR_VAL("subdiagonal index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY * m->tda;
	  v.size = GSL_MIN(m->size1 - k, m->size2);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, superdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;


	if (k >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY;
	  v.size = GSL_MIN(m->size1, m->size2 - k);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_GSL_COMPLEX_LONG

#define BASE_GSL_COMPLEX
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, row) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t i)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + i * MULTIPLICITY * m->tda;
	  v.size = m->size2;
	  v.stride = 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, column) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t j)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + j * MULTIPLICITY;
	  v.size = m->size1;
	  v.stride = m->tda;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, diagonal) (QUALIFIED_TYPE(gsl_matrix) * m)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	TYPE(gsl_vector) v = NULL_VECTOR;
	v.data = m->data;
	v.size = GSL_MIN(m->size1, m->size2);
	v.stride = m->tda + 1;
	v.block = m->block;
	v.owner = 0;

	view.vector = v;
	return view;
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, subdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (k >= m->size1)
	{
		GSL_ERROR_VAL("subdiagonal index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY * m->tda;
	  v.size = GSL_MIN(m->size1 - k, m->size2);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, superdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;


	if (k >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY;
	  v.size = GSL_MIN(m->size1, m->size2 - k);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_GSL_COMPLEX

#define BASE_GSL_COMPLEX_FLOAT
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, row) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t i)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + i * MULTIPLICITY * m->tda;
	  v.size = m->size2;
	  v.stride = 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, column) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t j)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + j * MULTIPLICITY;
	  v.size = m->size1;
	  v.stride = m->tda;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, diagonal) (QUALIFIED_TYPE(gsl_matrix) * m)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	TYPE(gsl_vector) v = NULL_VECTOR;
	v.data = m->data;
	v.size = GSL_MIN(m->size1, m->size2);
	v.stride = m->tda + 1;
	v.block = m->block;
	v.owner = 0;

	view.vector = v;
	return view;
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, subdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (k >= m->size1)
	{
		GSL_ERROR_VAL("subdiagonal index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY * m->tda;
	  v.size = GSL_MIN(m->size1 - k, m->size2);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, superdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;


	if (k >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY;
	  v.size = GSL_MIN(m->size1, m->size2 - k);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_GSL_COMPLEX_FLOAT

#define BASE_LONG_DOUBLE
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, row) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t i)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + i * MULTIPLICITY * m->tda;
	  v.size = m->size2;
	  v.stride = 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, column) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t j)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + j * MULTIPLICITY;
	  v.size = m->size1;
	  v.stride = m->tda;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, diagonal) (QUALIFIED_TYPE(gsl_matrix) * m)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	TYPE(gsl_vector) v = NULL_VECTOR;
	v.data = m->data;
	v.size = GSL_MIN(m->size1, m->size2);
	v.stride = m->tda + 1;
	v.block = m->block;
	v.owner = 0;

	view.vector = v;
	return view;
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, subdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (k >= m->size1)
	{
		GSL_ERROR_VAL("subdiagonal index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY * m->tda;
	  v.size = GSL_MIN(m->size1 - k, m->size2);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, superdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;


	if (k >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY;
	  v.size = GSL_MIN(m->size1, m->size2 - k);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, row) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t i)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + i * MULTIPLICITY * m->tda;
	  v.size = m->size2;
	  v.stride = 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, column) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t j)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + j * MULTIPLICITY;
	  v.size = m->size1;
	  v.stride = m->tda;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, diagonal) (QUALIFIED_TYPE(gsl_matrix) * m)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	TYPE(gsl_vector) v = NULL_VECTOR;
	v.data = m->data;
	v.size = GSL_MIN(m->size1, m->size2);
	v.stride = m->tda + 1;
	v.block = m->block;
	v.owner = 0;

	view.vector = v;
	return view;
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, subdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (k >= m->size1)
	{
		GSL_ERROR_VAL("subdiagonal index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY * m->tda;
	  v.size = GSL_MIN(m->size1 - k, m->size2);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, superdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;


	if (k >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY;
	  v.size = GSL_MIN(m->size1, m->size2 - k);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, row) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t i)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + i * MULTIPLICITY * m->tda;
	  v.size = m->size2;
	  v.stride = 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, column) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t j)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + j * MULTIPLICITY;
	  v.size = m->size1;
	  v.stride = m->tda;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, diagonal) (QUALIFIED_TYPE(gsl_matrix) * m)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	TYPE(gsl_vector) v = NULL_VECTOR;
	v.data = m->data;
	v.size = GSL_MIN(m->size1, m->size2);
	v.stride = m->tda + 1;
	v.block = m->block;
	v.owner = 0;

	view.vector = v;
	return view;
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, subdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (k >= m->size1)
	{
		GSL_ERROR_VAL("subdiagonal index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY * m->tda;
	  v.size = GSL_MIN(m->size1 - k, m->size2);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, superdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;


	if (k >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY;
	  v.size = GSL_MIN(m->size1, m->size2 - k);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_FLOAT

#define BASE_ULONG
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, row) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t i)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + i * MULTIPLICITY * m->tda;
	  v.size = m->size2;
	  v.stride = 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, column) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t j)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + j * MULTIPLICITY;
	  v.size = m->size1;
	  v.stride = m->tda;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, diagonal) (QUALIFIED_TYPE(gsl_matrix) * m)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	TYPE(gsl_vector) v = NULL_VECTOR;
	v.data = m->data;
	v.size = GSL_MIN(m->size1, m->size2);
	v.stride = m->tda + 1;
	v.block = m->block;
	v.owner = 0;

	view.vector = v;
	return view;
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, subdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (k >= m->size1)
	{
		GSL_ERROR_VAL("subdiagonal index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY * m->tda;
	  v.size = GSL_MIN(m->size1 - k, m->size2);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, superdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;


	if (k >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY;
	  v.size = GSL_MIN(m->size1, m->size2 - k);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_ULONG

#define BASE_LONG
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, row) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t i)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + i * MULTIPLICITY * m->tda;
	  v.size = m->size2;
	  v.stride = 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, column) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t j)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + j * MULTIPLICITY;
	  v.size = m->size1;
	  v.stride = m->tda;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, diagonal) (QUALIFIED_TYPE(gsl_matrix) * m)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	TYPE(gsl_vector) v = NULL_VECTOR;
	v.data = m->data;
	v.size = GSL_MIN(m->size1, m->size2);
	v.stride = m->tda + 1;
	v.block = m->block;
	v.owner = 0;

	view.vector = v;
	return view;
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, subdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (k >= m->size1)
	{
		GSL_ERROR_VAL("subdiagonal index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY * m->tda;
	  v.size = GSL_MIN(m->size1 - k, m->size2);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, superdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;


	if (k >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY;
	  v.size = GSL_MIN(m->size1, m->size2 - k);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_LONG

#define BASE_UINT
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, row) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t i)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + i * MULTIPLICITY * m->tda;
	  v.size = m->size2;
	  v.stride = 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, column) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t j)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + j * MULTIPLICITY;
	  v.size = m->size1;
	  v.stride = m->tda;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, diagonal) (QUALIFIED_TYPE(gsl_matrix) * m)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	TYPE(gsl_vector) v = NULL_VECTOR;
	v.data = m->data;
	v.size = GSL_MIN(m->size1, m->size2);
	v.stride = m->tda + 1;
	v.block = m->block;
	v.owner = 0;

	view.vector = v;
	return view;
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, subdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (k >= m->size1)
	{
		GSL_ERROR_VAL("subdiagonal index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY * m->tda;
	  v.size = GSL_MIN(m->size1 - k, m->size2);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, superdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;


	if (k >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY;
	  v.size = GSL_MIN(m->size1, m->size2 - k);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_UINT

#define BASE_INT
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, row) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t i)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + i * MULTIPLICITY * m->tda;
	  v.size = m->size2;
	  v.stride = 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, column) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t j)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + j * MULTIPLICITY;
	  v.size = m->size1;
	  v.stride = m->tda;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, diagonal) (QUALIFIED_TYPE(gsl_matrix) * m)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	TYPE(gsl_vector) v = NULL_VECTOR;
	v.data = m->data;
	v.size = GSL_MIN(m->size1, m->size2);
	v.stride = m->tda + 1;
	v.block = m->block;
	v.owner = 0;

	view.vector = v;
	return view;
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, subdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (k >= m->size1)
	{
		GSL_ERROR_VAL("subdiagonal index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY * m->tda;
	  v.size = GSL_MIN(m->size1 - k, m->size2);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, superdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;


	if (k >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY;
	  v.size = GSL_MIN(m->size1, m->size2 - k);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_INT

#define BASE_USHORT
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, row) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t i)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + i * MULTIPLICITY * m->tda;
	  v.size = m->size2;
	  v.stride = 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, column) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t j)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + j * MULTIPLICITY;
	  v.size = m->size1;
	  v.stride = m->tda;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, diagonal) (QUALIFIED_TYPE(gsl_matrix) * m)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	TYPE(gsl_vector) v = NULL_VECTOR;
	v.data = m->data;
	v.size = GSL_MIN(m->size1, m->size2);
	v.stride = m->tda + 1;
	v.block = m->block;
	v.owner = 0;

	view.vector = v;
	return view;
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, subdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (k >= m->size1)
	{
		GSL_ERROR_VAL("subdiagonal index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY * m->tda;
	  v.size = GSL_MIN(m->size1 - k, m->size2);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, superdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;


	if (k >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY;
	  v.size = GSL_MIN(m->size1, m->size2 - k);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_USHORT

#define BASE_SHORT
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, row) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t i)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + i * MULTIPLICITY * m->tda;
	  v.size = m->size2;
	  v.stride = 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, column) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t j)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + j * MULTIPLICITY;
	  v.size = m->size1;
	  v.stride = m->tda;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, diagonal) (QUALIFIED_TYPE(gsl_matrix) * m)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	TYPE(gsl_vector) v = NULL_VECTOR;
	v.data = m->data;
	v.size = GSL_MIN(m->size1, m->size2);
	v.stride = m->tda + 1;
	v.block = m->block;
	v.owner = 0;

	view.vector = v;
	return view;
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, subdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (k >= m->size1)
	{
		GSL_ERROR_VAL("subdiagonal index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY * m->tda;
	  v.size = GSL_MIN(m->size1 - k, m->size2);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, superdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;


	if (k >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY;
	  v.size = GSL_MIN(m->size1, m->size2 - k);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_SHORT

#define BASE_UCHAR
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, row) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t i)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + i * MULTIPLICITY * m->tda;
	  v.size = m->size2;
	  v.stride = 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, column) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t j)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + j * MULTIPLICITY;
	  v.size = m->size1;
	  v.stride = m->tda;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, diagonal) (QUALIFIED_TYPE(gsl_matrix) * m)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	TYPE(gsl_vector) v = NULL_VECTOR;
	v.data = m->data;
	v.size = GSL_MIN(m->size1, m->size2);
	v.stride = m->tda + 1;
	v.block = m->block;
	v.owner = 0;

	view.vector = v;
	return view;
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, subdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (k >= m->size1)
	{
		GSL_ERROR_VAL("subdiagonal index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY * m->tda;
	  v.size = GSL_MIN(m->size1 - k, m->size2);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, superdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;


	if (k >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY;
	  v.size = GSL_MIN(m->size1, m->size2 - k);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_UCHAR

#define BASE_CHAR
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, row) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t i)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + i * MULTIPLICITY * m->tda;
	  v.size = m->size2;
	  v.stride = 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, column) (QUALIFIED_TYPE(gsl_matrix) * m, const size_t j)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + j * MULTIPLICITY;
	  v.size = m->size1;
	  v.stride = m->tda;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, diagonal) (QUALIFIED_TYPE(gsl_matrix) * m)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	TYPE(gsl_vector) v = NULL_VECTOR;
	v.data = m->data;
	v.size = GSL_MIN(m->size1, m->size2);
	v.stride = m->tda + 1;
	v.block = m->block;
	v.owner = 0;

	view.vector = v;
	return view;
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, subdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;

	if (k >= m->size1)
	{
		GSL_ERROR_VAL("subdiagonal index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY * m->tda;
	  v.size = GSL_MIN(m->size1 - k, m->size2);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}

QUALIFIED_VIEW(_gsl_vector, view)
FUNCTION(gsl_matrix, superdiagonal) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t k)
{
	QUALIFIED_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;


	if (k >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_vector) v = NULL_VECTOR;

	  v.data = m->data + k * MULTIPLICITY;
	  v.size = GSL_MIN(m->size1, m->size2 - k);
	  v.stride = m->tda + 1;
	  v.block = m->block;
	  v.owner = 0;

	  view.vector = v;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_CHAR
