#include "stdafx.h"
#include <config.h.in>
#include <gsl/gsl_math.h>
#include <gsl/matrix/gsl_matrix.h>
#include <gsl/vector/gsl_vector.h>
#include <gsl/err/gsl_errno.h>

#include "view.h"

#define BASE_GSL_COMPLEX_LONG
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_matrix, view)
FUNCTION(gsl_matrix, submatrix) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const size_t n1, const size_t n2)
{
	QUALIFIED_VIEW(_gsl_matrix, view) view = NULL_MATRIX_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}
	else if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}
	else if (n1 == 0)
	{
		GSL_ERROR_VAL("first dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (n2 == 0)
	{
		GSL_ERROR_VAL("second dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (i + n1 > m->size1)
	{
		GSL_ERROR_VAL("first dimension overflows matrix", GSL_EINVAL, view);
	}
	else if (j + n2 > m->size2)
	{
		GSL_ERROR_VAL("second dimension overflows matrix", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_matrix) s = NULL_MATRIX;

	  s.data = m->data + MULTIPLICITY * (i * m->tda + j);
	  s.size1 = n1;
	  s.size2 = n2;
	  s.tda = m->tda;
	  s.block = m->block;
	  s.owner = 0;

	  view.matrix = s;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_GSL_COMPLEX_LONG

#define BASE_GSL_COMPLEX
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_matrix, view)
FUNCTION(gsl_matrix, submatrix) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const size_t n1, const size_t n2)
{
	QUALIFIED_VIEW(_gsl_matrix, view) view = NULL_MATRIX_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}
	else if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}
	else if (n1 == 0)
	{
		GSL_ERROR_VAL("first dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (n2 == 0)
	{
		GSL_ERROR_VAL("second dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (i + n1 > m->size1)
	{
		GSL_ERROR_VAL("first dimension overflows matrix", GSL_EINVAL, view);
	}
	else if (j + n2 > m->size2)
	{
		GSL_ERROR_VAL("second dimension overflows matrix", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_matrix) s = NULL_MATRIX;

	  s.data = m->data + MULTIPLICITY * (i * m->tda + j);
	  s.size1 = n1;
	  s.size2 = n2;
	  s.tda = m->tda;
	  s.block = m->block;
	  s.owner = 0;

	  view.matrix = s;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_GSL_COMPLEX

#define BASE_GSL_COMPLEX_FLOAT
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_matrix, view)
FUNCTION(gsl_matrix, submatrix) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const size_t n1, const size_t n2)
{
	QUALIFIED_VIEW(_gsl_matrix, view) view = NULL_MATRIX_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}
	else if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}
	else if (n1 == 0)
	{
		GSL_ERROR_VAL("first dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (n2 == 0)
	{
		GSL_ERROR_VAL("second dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (i + n1 > m->size1)
	{
		GSL_ERROR_VAL("first dimension overflows matrix", GSL_EINVAL, view);
	}
	else if (j + n2 > m->size2)
	{
		GSL_ERROR_VAL("second dimension overflows matrix", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_matrix) s = NULL_MATRIX;

	  s.data = m->data + MULTIPLICITY * (i * m->tda + j);
	  s.size1 = n1;
	  s.size2 = n2;
	  s.tda = m->tda;
	  s.block = m->block;
	  s.owner = 0;

	  view.matrix = s;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_GSL_COMPLEX_FLOAT

#define BASE_LONG_DOUBLE
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_matrix, view)
FUNCTION(gsl_matrix, submatrix) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const size_t n1, const size_t n2)
{
	QUALIFIED_VIEW(_gsl_matrix, view) view = NULL_MATRIX_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}
	else if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}
	else if (n1 == 0)
	{
		GSL_ERROR_VAL("first dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (n2 == 0)
	{
		GSL_ERROR_VAL("second dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (i + n1 > m->size1)
	{
		GSL_ERROR_VAL("first dimension overflows matrix", GSL_EINVAL, view);
	}
	else if (j + n2 > m->size2)
	{
		GSL_ERROR_VAL("second dimension overflows matrix", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_matrix) s = NULL_MATRIX;

	  s.data = m->data + MULTIPLICITY * (i * m->tda + j);
	  s.size1 = n1;
	  s.size2 = n2;
	  s.tda = m->tda;
	  s.block = m->block;
	  s.owner = 0;

	  view.matrix = s;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_matrix, view)
FUNCTION(gsl_matrix, submatrix) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const size_t n1, const size_t n2)
{
	QUALIFIED_VIEW(_gsl_matrix, view) view = NULL_MATRIX_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}
	else if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}
	else if (n1 == 0)
	{
		GSL_ERROR_VAL("first dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (n2 == 0)
	{
		GSL_ERROR_VAL("second dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (i + n1 > m->size1)
	{
		GSL_ERROR_VAL("first dimension overflows matrix", GSL_EINVAL, view);
	}
	else if (j + n2 > m->size2)
	{
		GSL_ERROR_VAL("second dimension overflows matrix", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_matrix) s = NULL_MATRIX;

	  s.data = m->data + MULTIPLICITY * (i * m->tda + j);
	  s.size1 = n1;
	  s.size2 = n2;
	  s.tda = m->tda;
	  s.block = m->block;
	  s.owner = 0;

	  view.matrix = s;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_matrix, view)
FUNCTION(gsl_matrix, submatrix) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const size_t n1, const size_t n2)
{
	QUALIFIED_VIEW(_gsl_matrix, view) view = NULL_MATRIX_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}
	else if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}
	else if (n1 == 0)
	{
		GSL_ERROR_VAL("first dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (n2 == 0)
	{
		GSL_ERROR_VAL("second dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (i + n1 > m->size1)
	{
		GSL_ERROR_VAL("first dimension overflows matrix", GSL_EINVAL, view);
	}
	else if (j + n2 > m->size2)
	{
		GSL_ERROR_VAL("second dimension overflows matrix", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_matrix) s = NULL_MATRIX;

	  s.data = m->data + MULTIPLICITY * (i * m->tda + j);
	  s.size1 = n1;
	  s.size2 = n2;
	  s.tda = m->tda;
	  s.block = m->block;
	  s.owner = 0;

	  view.matrix = s;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_FLOAT

#define BASE_ULONG
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_matrix, view)
FUNCTION(gsl_matrix, submatrix) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const size_t n1, const size_t n2)
{
	QUALIFIED_VIEW(_gsl_matrix, view) view = NULL_MATRIX_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}
	else if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}
	else if (n1 == 0)
	{
		GSL_ERROR_VAL("first dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (n2 == 0)
	{
		GSL_ERROR_VAL("second dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (i + n1 > m->size1)
	{
		GSL_ERROR_VAL("first dimension overflows matrix", GSL_EINVAL, view);
	}
	else if (j + n2 > m->size2)
	{
		GSL_ERROR_VAL("second dimension overflows matrix", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_matrix) s = NULL_MATRIX;

	  s.data = m->data + MULTIPLICITY * (i * m->tda + j);
	  s.size1 = n1;
	  s.size2 = n2;
	  s.tda = m->tda;
	  s.block = m->block;
	  s.owner = 0;

	  view.matrix = s;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_ULONG

#define BASE_LONG
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_matrix, view)
FUNCTION(gsl_matrix, submatrix) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const size_t n1, const size_t n2)
{
	QUALIFIED_VIEW(_gsl_matrix, view) view = NULL_MATRIX_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}
	else if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}
	else if (n1 == 0)
	{
		GSL_ERROR_VAL("first dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (n2 == 0)
	{
		GSL_ERROR_VAL("second dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (i + n1 > m->size1)
	{
		GSL_ERROR_VAL("first dimension overflows matrix", GSL_EINVAL, view);
	}
	else if (j + n2 > m->size2)
	{
		GSL_ERROR_VAL("second dimension overflows matrix", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_matrix) s = NULL_MATRIX;

	  s.data = m->data + MULTIPLICITY * (i * m->tda + j);
	  s.size1 = n1;
	  s.size2 = n2;
	  s.tda = m->tda;
	  s.block = m->block;
	  s.owner = 0;

	  view.matrix = s;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_LONG

#define BASE_UINT
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_matrix, view)
FUNCTION(gsl_matrix, submatrix) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const size_t n1, const size_t n2)
{
	QUALIFIED_VIEW(_gsl_matrix, view) view = NULL_MATRIX_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}
	else if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}
	else if (n1 == 0)
	{
		GSL_ERROR_VAL("first dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (n2 == 0)
	{
		GSL_ERROR_VAL("second dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (i + n1 > m->size1)
	{
		GSL_ERROR_VAL("first dimension overflows matrix", GSL_EINVAL, view);
	}
	else if (j + n2 > m->size2)
	{
		GSL_ERROR_VAL("second dimension overflows matrix", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_matrix) s = NULL_MATRIX;

	  s.data = m->data + MULTIPLICITY * (i * m->tda + j);
	  s.size1 = n1;
	  s.size2 = n2;
	  s.tda = m->tda;
	  s.block = m->block;
	  s.owner = 0;

	  view.matrix = s;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_UINT

#define BASE_INT
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_matrix, view)
FUNCTION(gsl_matrix, submatrix) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const size_t n1, const size_t n2)
{
	QUALIFIED_VIEW(_gsl_matrix, view) view = NULL_MATRIX_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}
	else if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}
	else if (n1 == 0)
	{
		GSL_ERROR_VAL("first dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (n2 == 0)
	{
		GSL_ERROR_VAL("second dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (i + n1 > m->size1)
	{
		GSL_ERROR_VAL("first dimension overflows matrix", GSL_EINVAL, view);
	}
	else if (j + n2 > m->size2)
	{
		GSL_ERROR_VAL("second dimension overflows matrix", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_matrix) s = NULL_MATRIX;

	  s.data = m->data + MULTIPLICITY * (i * m->tda + j);
	  s.size1 = n1;
	  s.size2 = n2;
	  s.tda = m->tda;
	  s.block = m->block;
	  s.owner = 0;

	  view.matrix = s;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_INT

#define BASE_USHORT
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_matrix, view)
FUNCTION(gsl_matrix, submatrix) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const size_t n1, const size_t n2)
{
	QUALIFIED_VIEW(_gsl_matrix, view) view = NULL_MATRIX_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}
	else if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}
	else if (n1 == 0)
	{
		GSL_ERROR_VAL("first dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (n2 == 0)
	{
		GSL_ERROR_VAL("second dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (i + n1 > m->size1)
	{
		GSL_ERROR_VAL("first dimension overflows matrix", GSL_EINVAL, view);
	}
	else if (j + n2 > m->size2)
	{
		GSL_ERROR_VAL("second dimension overflows matrix", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_matrix) s = NULL_MATRIX;

	  s.data = m->data + MULTIPLICITY * (i * m->tda + j);
	  s.size1 = n1;
	  s.size2 = n2;
	  s.tda = m->tda;
	  s.block = m->block;
	  s.owner = 0;

	  view.matrix = s;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_USHORT

#define BASE_SHORT
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_matrix, view)
FUNCTION(gsl_matrix, submatrix) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const size_t n1, const size_t n2)
{
	QUALIFIED_VIEW(_gsl_matrix, view) view = NULL_MATRIX_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}
	else if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}
	else if (n1 == 0)
	{
		GSL_ERROR_VAL("first dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (n2 == 0)
	{
		GSL_ERROR_VAL("second dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (i + n1 > m->size1)
	{
		GSL_ERROR_VAL("first dimension overflows matrix", GSL_EINVAL, view);
	}
	else if (j + n2 > m->size2)
	{
		GSL_ERROR_VAL("second dimension overflows matrix", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_matrix) s = NULL_MATRIX;

	  s.data = m->data + MULTIPLICITY * (i * m->tda + j);
	  s.size1 = n1;
	  s.size2 = n2;
	  s.tda = m->tda;
	  s.block = m->block;
	  s.owner = 0;

	  view.matrix = s;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_SHORT

#define BASE_UCHAR
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_matrix, view)
FUNCTION(gsl_matrix, submatrix) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const size_t n1, const size_t n2)
{
	QUALIFIED_VIEW(_gsl_matrix, view) view = NULL_MATRIX_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}
	else if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}
	else if (n1 == 0)
	{
		GSL_ERROR_VAL("first dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (n2 == 0)
	{
		GSL_ERROR_VAL("second dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (i + n1 > m->size1)
	{
		GSL_ERROR_VAL("first dimension overflows matrix", GSL_EINVAL, view);
	}
	else if (j + n2 > m->size2)
	{
		GSL_ERROR_VAL("second dimension overflows matrix", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_matrix) s = NULL_MATRIX;

	  s.data = m->data + MULTIPLICITY * (i * m->tda + j);
	  s.size1 = n1;
	  s.size2 = n2;
	  s.tda = m->tda;
	  s.block = m->block;
	  s.owner = 0;

	  view.matrix = s;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_UCHAR

#define BASE_CHAR
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_matrix, view)
FUNCTION(gsl_matrix, submatrix) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const size_t n1, const size_t n2)
{
	QUALIFIED_VIEW(_gsl_matrix, view) view = NULL_MATRIX_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}
	else if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}
	else if (n1 == 0)
	{
		GSL_ERROR_VAL("first dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (n2 == 0)
	{
		GSL_ERROR_VAL("second dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (i + n1 > m->size1)
	{
		GSL_ERROR_VAL("first dimension overflows matrix", GSL_EINVAL, view);
	}
	else if (j + n2 > m->size2)
	{
		GSL_ERROR_VAL("second dimension overflows matrix", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_matrix) s = NULL_MATRIX;

	  s.data = m->data + MULTIPLICITY * (i * m->tda + j);
	  s.size1 = n1;
	  s.size2 = n2;
	  s.tda = m->tda;
	  s.block = m->block;
	  s.owner = 0;

	  view.matrix = s;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_CHAR

#define USE_QUALIFIER
#define QUALIFIER const

#define BASE_GSL_COMPLEX_LONG
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_matrix, view)
FUNCTION(gsl_matrix, submatrix) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const size_t n1, const size_t n2)
{
	QUALIFIED_VIEW(_gsl_matrix, view) view = NULL_MATRIX_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}
	else if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}
	else if (n1 == 0)
	{
		GSL_ERROR_VAL("first dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (n2 == 0)
	{
		GSL_ERROR_VAL("second dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (i + n1 > m->size1)
	{
		GSL_ERROR_VAL("first dimension overflows matrix", GSL_EINVAL, view);
	}
	else if (j + n2 > m->size2)
	{
		GSL_ERROR_VAL("second dimension overflows matrix", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_matrix) s = NULL_MATRIX;

	  s.data = m->data + MULTIPLICITY * (i * m->tda + j);
	  s.size1 = n1;
	  s.size2 = n2;
	  s.tda = m->tda;
	  s.block = m->block;
	  s.owner = 0;

	  view.matrix = s;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_GSL_COMPLEX_LONG

#define BASE_GSL_COMPLEX
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_matrix, view)
FUNCTION(gsl_matrix, submatrix) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const size_t n1, const size_t n2)
{
	QUALIFIED_VIEW(_gsl_matrix, view) view = NULL_MATRIX_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}
	else if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}
	else if (n1 == 0)
	{
		GSL_ERROR_VAL("first dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (n2 == 0)
	{
		GSL_ERROR_VAL("second dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (i + n1 > m->size1)
	{
		GSL_ERROR_VAL("first dimension overflows matrix", GSL_EINVAL, view);
	}
	else if (j + n2 > m->size2)
	{
		GSL_ERROR_VAL("second dimension overflows matrix", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_matrix) s = NULL_MATRIX;

	  s.data = m->data + MULTIPLICITY * (i * m->tda + j);
	  s.size1 = n1;
	  s.size2 = n2;
	  s.tda = m->tda;
	  s.block = m->block;
	  s.owner = 0;

	  view.matrix = s;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_GSL_COMPLEX

#define BASE_GSL_COMPLEX_FLOAT
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_matrix, view)
FUNCTION(gsl_matrix, submatrix) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const size_t n1, const size_t n2)
{
	QUALIFIED_VIEW(_gsl_matrix, view) view = NULL_MATRIX_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}
	else if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}
	else if (n1 == 0)
	{
		GSL_ERROR_VAL("first dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (n2 == 0)
	{
		GSL_ERROR_VAL("second dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (i + n1 > m->size1)
	{
		GSL_ERROR_VAL("first dimension overflows matrix", GSL_EINVAL, view);
	}
	else if (j + n2 > m->size2)
	{
		GSL_ERROR_VAL("second dimension overflows matrix", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_matrix) s = NULL_MATRIX;

	  s.data = m->data + MULTIPLICITY * (i * m->tda + j);
	  s.size1 = n1;
	  s.size2 = n2;
	  s.tda = m->tda;
	  s.block = m->block;
	  s.owner = 0;

	  view.matrix = s;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_GSL_COMPLEX_FLOAT

#define BASE_LONG_DOUBLE
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_matrix, view)
FUNCTION(gsl_matrix, submatrix) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const size_t n1, const size_t n2)
{
	QUALIFIED_VIEW(_gsl_matrix, view) view = NULL_MATRIX_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}
	else if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}
	else if (n1 == 0)
	{
		GSL_ERROR_VAL("first dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (n2 == 0)
	{
		GSL_ERROR_VAL("second dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (i + n1 > m->size1)
	{
		GSL_ERROR_VAL("first dimension overflows matrix", GSL_EINVAL, view);
	}
	else if (j + n2 > m->size2)
	{
		GSL_ERROR_VAL("second dimension overflows matrix", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_matrix) s = NULL_MATRIX;

	  s.data = m->data + MULTIPLICITY * (i * m->tda + j);
	  s.size1 = n1;
	  s.size2 = n2;
	  s.tda = m->tda;
	  s.block = m->block;
	  s.owner = 0;

	  view.matrix = s;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_matrix, view)
FUNCTION(gsl_matrix, submatrix) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const size_t n1, const size_t n2)
{
	QUALIFIED_VIEW(_gsl_matrix, view) view = NULL_MATRIX_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}
	else if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}
	else if (n1 == 0)
	{
		GSL_ERROR_VAL("first dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (n2 == 0)
	{
		GSL_ERROR_VAL("second dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (i + n1 > m->size1)
	{
		GSL_ERROR_VAL("first dimension overflows matrix", GSL_EINVAL, view);
	}
	else if (j + n2 > m->size2)
	{
		GSL_ERROR_VAL("second dimension overflows matrix", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_matrix) s = NULL_MATRIX;

	  s.data = m->data + MULTIPLICITY * (i * m->tda + j);
	  s.size1 = n1;
	  s.size2 = n2;
	  s.tda = m->tda;
	  s.block = m->block;
	  s.owner = 0;

	  view.matrix = s;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_matrix, view)
FUNCTION(gsl_matrix, submatrix) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const size_t n1, const size_t n2)
{
	QUALIFIED_VIEW(_gsl_matrix, view) view = NULL_MATRIX_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}
	else if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}
	else if (n1 == 0)
	{
		GSL_ERROR_VAL("first dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (n2 == 0)
	{
		GSL_ERROR_VAL("second dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (i + n1 > m->size1)
	{
		GSL_ERROR_VAL("first dimension overflows matrix", GSL_EINVAL, view);
	}
	else if (j + n2 > m->size2)
	{
		GSL_ERROR_VAL("second dimension overflows matrix", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_matrix) s = NULL_MATRIX;

	  s.data = m->data + MULTIPLICITY * (i * m->tda + j);
	  s.size1 = n1;
	  s.size2 = n2;
	  s.tda = m->tda;
	  s.block = m->block;
	  s.owner = 0;

	  view.matrix = s;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_FLOAT

#define BASE_ULONG
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_matrix, view)
FUNCTION(gsl_matrix, submatrix) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const size_t n1, const size_t n2)
{
	QUALIFIED_VIEW(_gsl_matrix, view) view = NULL_MATRIX_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}
	else if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}
	else if (n1 == 0)
	{
		GSL_ERROR_VAL("first dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (n2 == 0)
	{
		GSL_ERROR_VAL("second dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (i + n1 > m->size1)
	{
		GSL_ERROR_VAL("first dimension overflows matrix", GSL_EINVAL, view);
	}
	else if (j + n2 > m->size2)
	{
		GSL_ERROR_VAL("second dimension overflows matrix", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_matrix) s = NULL_MATRIX;

	  s.data = m->data + MULTIPLICITY * (i * m->tda + j);
	  s.size1 = n1;
	  s.size2 = n2;
	  s.tda = m->tda;
	  s.block = m->block;
	  s.owner = 0;

	  view.matrix = s;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_ULONG

#define BASE_LONG
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_matrix, view)
FUNCTION(gsl_matrix, submatrix) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const size_t n1, const size_t n2)
{
	QUALIFIED_VIEW(_gsl_matrix, view) view = NULL_MATRIX_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}
	else if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}
	else if (n1 == 0)
	{
		GSL_ERROR_VAL("first dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (n2 == 0)
	{
		GSL_ERROR_VAL("second dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (i + n1 > m->size1)
	{
		GSL_ERROR_VAL("first dimension overflows matrix", GSL_EINVAL, view);
	}
	else if (j + n2 > m->size2)
	{
		GSL_ERROR_VAL("second dimension overflows matrix", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_matrix) s = NULL_MATRIX;

	  s.data = m->data + MULTIPLICITY * (i * m->tda + j);
	  s.size1 = n1;
	  s.size2 = n2;
	  s.tda = m->tda;
	  s.block = m->block;
	  s.owner = 0;

	  view.matrix = s;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_LONG

#define BASE_UINT
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_matrix, view)
FUNCTION(gsl_matrix, submatrix) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const size_t n1, const size_t n2)
{
	QUALIFIED_VIEW(_gsl_matrix, view) view = NULL_MATRIX_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}
	else if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}
	else if (n1 == 0)
	{
		GSL_ERROR_VAL("first dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (n2 == 0)
	{
		GSL_ERROR_VAL("second dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (i + n1 > m->size1)
	{
		GSL_ERROR_VAL("first dimension overflows matrix", GSL_EINVAL, view);
	}
	else if (j + n2 > m->size2)
	{
		GSL_ERROR_VAL("second dimension overflows matrix", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_matrix) s = NULL_MATRIX;

	  s.data = m->data + MULTIPLICITY * (i * m->tda + j);
	  s.size1 = n1;
	  s.size2 = n2;
	  s.tda = m->tda;
	  s.block = m->block;
	  s.owner = 0;

	  view.matrix = s;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_UINT

#define BASE_INT
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_matrix, view)
FUNCTION(gsl_matrix, submatrix) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const size_t n1, const size_t n2)
{
	QUALIFIED_VIEW(_gsl_matrix, view) view = NULL_MATRIX_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}
	else if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}
	else if (n1 == 0)
	{
		GSL_ERROR_VAL("first dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (n2 == 0)
	{
		GSL_ERROR_VAL("second dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (i + n1 > m->size1)
	{
		GSL_ERROR_VAL("first dimension overflows matrix", GSL_EINVAL, view);
	}
	else if (j + n2 > m->size2)
	{
		GSL_ERROR_VAL("second dimension overflows matrix", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_matrix) s = NULL_MATRIX;

	  s.data = m->data + MULTIPLICITY * (i * m->tda + j);
	  s.size1 = n1;
	  s.size2 = n2;
	  s.tda = m->tda;
	  s.block = m->block;
	  s.owner = 0;

	  view.matrix = s;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_INT

#define BASE_USHORT
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_matrix, view)
FUNCTION(gsl_matrix, submatrix) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const size_t n1, const size_t n2)
{
	QUALIFIED_VIEW(_gsl_matrix, view) view = NULL_MATRIX_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}
	else if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}
	else if (n1 == 0)
	{
		GSL_ERROR_VAL("first dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (n2 == 0)
	{
		GSL_ERROR_VAL("second dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (i + n1 > m->size1)
	{
		GSL_ERROR_VAL("first dimension overflows matrix", GSL_EINVAL, view);
	}
	else if (j + n2 > m->size2)
	{
		GSL_ERROR_VAL("second dimension overflows matrix", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_matrix) s = NULL_MATRIX;

	  s.data = m->data + MULTIPLICITY * (i * m->tda + j);
	  s.size1 = n1;
	  s.size2 = n2;
	  s.tda = m->tda;
	  s.block = m->block;
	  s.owner = 0;

	  view.matrix = s;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_USHORT

#define BASE_SHORT
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_matrix, view)
FUNCTION(gsl_matrix, submatrix) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const size_t n1, const size_t n2)
{
	QUALIFIED_VIEW(_gsl_matrix, view) view = NULL_MATRIX_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}
	else if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}
	else if (n1 == 0)
	{
		GSL_ERROR_VAL("first dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (n2 == 0)
	{
		GSL_ERROR_VAL("second dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (i + n1 > m->size1)
	{
		GSL_ERROR_VAL("first dimension overflows matrix", GSL_EINVAL, view);
	}
	else if (j + n2 > m->size2)
	{
		GSL_ERROR_VAL("second dimension overflows matrix", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_matrix) s = NULL_MATRIX;

	  s.data = m->data + MULTIPLICITY * (i * m->tda + j);
	  s.size1 = n1;
	  s.size2 = n2;
	  s.tda = m->tda;
	  s.block = m->block;
	  s.owner = 0;

	  view.matrix = s;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_SHORT

#define BASE_UCHAR
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_matrix, view)
FUNCTION(gsl_matrix, submatrix) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const size_t n1, const size_t n2)
{
	QUALIFIED_VIEW(_gsl_matrix, view) view = NULL_MATRIX_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}
	else if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}
	else if (n1 == 0)
	{
		GSL_ERROR_VAL("first dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (n2 == 0)
	{
		GSL_ERROR_VAL("second dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (i + n1 > m->size1)
	{
		GSL_ERROR_VAL("first dimension overflows matrix", GSL_EINVAL, view);
	}
	else if (j + n2 > m->size2)
	{
		GSL_ERROR_VAL("second dimension overflows matrix", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_matrix) s = NULL_MATRIX;

	  s.data = m->data + MULTIPLICITY * (i * m->tda + j);
	  s.size1 = n1;
	  s.size2 = n2;
	  s.tda = m->tda;
	  s.block = m->block;
	  s.owner = 0;

	  view.matrix = s;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_UCHAR

#define BASE_CHAR
#include "templates_on.h"
QUALIFIED_VIEW(_gsl_matrix, view)
FUNCTION(gsl_matrix, submatrix) (QUALIFIED_TYPE(gsl_matrix) * m,
const size_t i, const size_t j,
const size_t n1, const size_t n2)
{
	QUALIFIED_VIEW(_gsl_matrix, view) view = NULL_MATRIX_VIEW;

	if (i >= m->size1)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, view);
	}
	else if (j >= m->size2)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, view);
	}
	else if (n1 == 0)
	{
		GSL_ERROR_VAL("first dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (n2 == 0)
	{
		GSL_ERROR_VAL("second dimension must be non-zero", GSL_EINVAL, view);
	}
	else if (i + n1 > m->size1)
	{
		GSL_ERROR_VAL("first dimension overflows matrix", GSL_EINVAL, view);
	}
	else if (j + n2 > m->size2)
	{
		GSL_ERROR_VAL("second dimension overflows matrix", GSL_EINVAL, view);
	}

  {
	  TYPE(gsl_matrix) s = NULL_MATRIX;

	  s.data = m->data + MULTIPLICITY * (i * m->tda + j);
	  s.size1 = n1;
	  s.size2 = n2;
	  s.tda = m->tda;
	  s.block = m->block;
	  s.owner = 0;

	  view.matrix = s;
	  return view;
  }
}


#include "templates_off.h"
#undef  BASE_CHAR
