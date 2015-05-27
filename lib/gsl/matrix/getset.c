#include "stdafx.h"
#include <config.h.in>
#include <gsl/err/gsl_errno.h>
#include <gsl/matrix/gsl_matrix.h>
#include <gsl/vector/gsl_vector.h>

#define BASE_GSL_COMPLEX_LONG
#include "templates_on.h"
int
FUNCTION(gsl_matrix, get_row) (TYPE(gsl_vector) * v,
const TYPE(gsl_matrix) * m,
const size_t i)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (i >= M)
	{
		GSL_ERROR("row index is out of range", GSL_EINVAL);
	}

	if (v->size != N)
	{
		GSL_ERROR("matrix row size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  ATOMIC *v_data = v->data;
	  const ATOMIC *row_data = m->data + MULTIPLICITY * i * tda;
	  const size_t stride = v->stride;
	  size_t j;

	  for (j = 0; j < N; j++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  v_data[MULTIPLICITY * stride * j + k]
				  = row_data[MULTIPLICITY * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, get_col) (TYPE(gsl_vector) * v,
const TYPE(gsl_matrix) * m,
const size_t j)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (j >= N)
	{
		GSL_ERROR("column index is out of range", GSL_EINVAL);
	}

	if (v->size != M)
	{
		GSL_ERROR("matrix column size and vector length are not equal",
			GSL_EBADLEN);
	}


  {
	  ATOMIC *v_data = v->data;
	  const ATOMIC *column_data = m->data + MULTIPLICITY * j;
	  const size_t stride = v->stride;
	  size_t i;

	  for (i = 0; i < M; i++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  v_data[stride * MULTIPLICITY * i + k] =
				  column_data[MULTIPLICITY * i * tda + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, set_row) (TYPE(gsl_matrix) * m,
const size_t i,
const TYPE(gsl_vector) * v)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (i >= M)
	{
		GSL_ERROR("row index is out of range", GSL_EINVAL);
	}

	if (v->size != N)
	{
		GSL_ERROR("matrix row size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  const ATOMIC *v_data = v->data;
	  ATOMIC *row_data = m->data + MULTIPLICITY * i * tda;
	  const size_t stride = v->stride;
	  size_t j;

	  for (j = 0; j < N; j++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  row_data[MULTIPLICITY*j + k]
				  = v_data[MULTIPLICITY * stride * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, set_col) (TYPE(gsl_matrix) * m,
const size_t j,
const TYPE(gsl_vector) * v)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (j >= N)
	{
		GSL_ERROR("column index is out of range", GSL_EINVAL);
	}

	if (v->size != M)
	{
		GSL_ERROR("matrix column size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  const ATOMIC *v_data = v->data;
	  ATOMIC *column_data = m->data + MULTIPLICITY * j;
	  const size_t stride = v->stride;
	  size_t i;

	  for (i = 0; i < M; i++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  column_data[MULTIPLICITY * i * tda + k]
				  = v_data[MULTIPLICITY * stride * i + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}


TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_row_from_matrix) (TYPE(gsl_matrix) * m,
const size_t i)
{
	TYPE(gsl_vector) * v;

	const size_t M = m->size1;

	if (i >= M)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = m->data + MULTIPLICITY * i * m->tda;
	v->size = m->size2;
	v->stride = 1;
	v->block = 0;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_col_from_matrix) (TYPE(gsl_matrix) * m,
const size_t j)
{
	TYPE(gsl_vector) * v;

	const size_t N = m->size2;

	if (j >= N)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = m->data + MULTIPLICITY * j;
	v->size = m->size1;
	v->stride = m->tda;
	v->block = 0;

	return v;
}


#include "templates_off.h"
#undef  BASE_GSL_COMPLEX_LONG

#define BASE_GSL_COMPLEX
#include "templates_on.h"
int
FUNCTION(gsl_matrix, get_row) (TYPE(gsl_vector) * v,
const TYPE(gsl_matrix) * m,
const size_t i)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (i >= M)
	{
		GSL_ERROR("row index is out of range", GSL_EINVAL);
	}

	if (v->size != N)
	{
		GSL_ERROR("matrix row size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  ATOMIC *v_data = v->data;
	  const ATOMIC *row_data = m->data + MULTIPLICITY * i * tda;
	  const size_t stride = v->stride;
	  size_t j;

	  for (j = 0; j < N; j++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  v_data[MULTIPLICITY * stride * j + k]
				  = row_data[MULTIPLICITY * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, get_col) (TYPE(gsl_vector) * v,
const TYPE(gsl_matrix) * m,
const size_t j)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (j >= N)
	{
		GSL_ERROR("column index is out of range", GSL_EINVAL);
	}

	if (v->size != M)
	{
		GSL_ERROR("matrix column size and vector length are not equal",
			GSL_EBADLEN);
	}


  {
	  ATOMIC *v_data = v->data;
	  const ATOMIC *column_data = m->data + MULTIPLICITY * j;
	  const size_t stride = v->stride;
	  size_t i;

	  for (i = 0; i < M; i++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  v_data[stride * MULTIPLICITY * i + k] =
				  column_data[MULTIPLICITY * i * tda + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, set_row) (TYPE(gsl_matrix) * m,
const size_t i,
const TYPE(gsl_vector) * v)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (i >= M)
	{
		GSL_ERROR("row index is out of range", GSL_EINVAL);
	}

	if (v->size != N)
	{
		GSL_ERROR("matrix row size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  const ATOMIC *v_data = v->data;
	  ATOMIC *row_data = m->data + MULTIPLICITY * i * tda;
	  const size_t stride = v->stride;
	  size_t j;

	  for (j = 0; j < N; j++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  row_data[MULTIPLICITY*j + k]
				  = v_data[MULTIPLICITY * stride * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, set_col) (TYPE(gsl_matrix) * m,
const size_t j,
const TYPE(gsl_vector) * v)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (j >= N)
	{
		GSL_ERROR("column index is out of range", GSL_EINVAL);
	}

	if (v->size != M)
	{
		GSL_ERROR("matrix column size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  const ATOMIC *v_data = v->data;
	  ATOMIC *column_data = m->data + MULTIPLICITY * j;
	  const size_t stride = v->stride;
	  size_t i;

	  for (i = 0; i < M; i++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  column_data[MULTIPLICITY * i * tda + k]
				  = v_data[MULTIPLICITY * stride * i + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}


TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_row_from_matrix) (TYPE(gsl_matrix) * m,
const size_t i)
{
	TYPE(gsl_vector) * v;

	const size_t M = m->size1;

	if (i >= M)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = m->data + MULTIPLICITY * i * m->tda;
	v->size = m->size2;
	v->stride = 1;
	v->block = 0;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_col_from_matrix) (TYPE(gsl_matrix) * m,
const size_t j)
{
	TYPE(gsl_vector) * v;

	const size_t N = m->size2;

	if (j >= N)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = m->data + MULTIPLICITY * j;
	v->size = m->size1;
	v->stride = m->tda;
	v->block = 0;

	return v;
}


#include "templates_off.h"
#undef  BASE_GSL_COMPLEX

#define BASE_GSL_COMPLEX_FLOAT
#include "templates_on.h"
int
FUNCTION(gsl_matrix, get_row) (TYPE(gsl_vector) * v,
const TYPE(gsl_matrix) * m,
const size_t i)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (i >= M)
	{
		GSL_ERROR("row index is out of range", GSL_EINVAL);
	}

	if (v->size != N)
	{
		GSL_ERROR("matrix row size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  ATOMIC *v_data = v->data;
	  const ATOMIC *row_data = m->data + MULTIPLICITY * i * tda;
	  const size_t stride = v->stride;
	  size_t j;

	  for (j = 0; j < N; j++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  v_data[MULTIPLICITY * stride * j + k]
				  = row_data[MULTIPLICITY * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, get_col) (TYPE(gsl_vector) * v,
const TYPE(gsl_matrix) * m,
const size_t j)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (j >= N)
	{
		GSL_ERROR("column index is out of range", GSL_EINVAL);
	}

	if (v->size != M)
	{
		GSL_ERROR("matrix column size and vector length are not equal",
			GSL_EBADLEN);
	}


  {
	  ATOMIC *v_data = v->data;
	  const ATOMIC *column_data = m->data + MULTIPLICITY * j;
	  const size_t stride = v->stride;
	  size_t i;

	  for (i = 0; i < M; i++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  v_data[stride * MULTIPLICITY * i + k] =
				  column_data[MULTIPLICITY * i * tda + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, set_row) (TYPE(gsl_matrix) * m,
const size_t i,
const TYPE(gsl_vector) * v)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (i >= M)
	{
		GSL_ERROR("row index is out of range", GSL_EINVAL);
	}

	if (v->size != N)
	{
		GSL_ERROR("matrix row size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  const ATOMIC *v_data = v->data;
	  ATOMIC *row_data = m->data + MULTIPLICITY * i * tda;
	  const size_t stride = v->stride;
	  size_t j;

	  for (j = 0; j < N; j++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  row_data[MULTIPLICITY*j + k]
				  = v_data[MULTIPLICITY * stride * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, set_col) (TYPE(gsl_matrix) * m,
const size_t j,
const TYPE(gsl_vector) * v)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (j >= N)
	{
		GSL_ERROR("column index is out of range", GSL_EINVAL);
	}

	if (v->size != M)
	{
		GSL_ERROR("matrix column size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  const ATOMIC *v_data = v->data;
	  ATOMIC *column_data = m->data + MULTIPLICITY * j;
	  const size_t stride = v->stride;
	  size_t i;

	  for (i = 0; i < M; i++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  column_data[MULTIPLICITY * i * tda + k]
				  = v_data[MULTIPLICITY * stride * i + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}


TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_row_from_matrix) (TYPE(gsl_matrix) * m,
const size_t i)
{
	TYPE(gsl_vector) * v;

	const size_t M = m->size1;

	if (i >= M)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = m->data + MULTIPLICITY * i * m->tda;
	v->size = m->size2;
	v->stride = 1;
	v->block = 0;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_col_from_matrix) (TYPE(gsl_matrix) * m,
const size_t j)
{
	TYPE(gsl_vector) * v;

	const size_t N = m->size2;

	if (j >= N)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = m->data + MULTIPLICITY * j;
	v->size = m->size1;
	v->stride = m->tda;
	v->block = 0;

	return v;
}


#include "templates_off.h"
#undef  BASE_GSL_COMPLEX_FLOAT

#define BASE_LONG_DOUBLE
#include "templates_on.h"
int
FUNCTION(gsl_matrix, get_row) (TYPE(gsl_vector) * v,
const TYPE(gsl_matrix) * m,
const size_t i)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (i >= M)
	{
		GSL_ERROR("row index is out of range", GSL_EINVAL);
	}

	if (v->size != N)
	{
		GSL_ERROR("matrix row size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  ATOMIC *v_data = v->data;
	  const ATOMIC *row_data = m->data + MULTIPLICITY * i * tda;
	  const size_t stride = v->stride;
	  size_t j;

	  for (j = 0; j < N; j++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  v_data[MULTIPLICITY * stride * j + k]
				  = row_data[MULTIPLICITY * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, get_col) (TYPE(gsl_vector) * v,
const TYPE(gsl_matrix) * m,
const size_t j)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (j >= N)
	{
		GSL_ERROR("column index is out of range", GSL_EINVAL);
	}

	if (v->size != M)
	{
		GSL_ERROR("matrix column size and vector length are not equal",
			GSL_EBADLEN);
	}


  {
	  ATOMIC *v_data = v->data;
	  const ATOMIC *column_data = m->data + MULTIPLICITY * j;
	  const size_t stride = v->stride;
	  size_t i;

	  for (i = 0; i < M; i++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  v_data[stride * MULTIPLICITY * i + k] =
				  column_data[MULTIPLICITY * i * tda + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, set_row) (TYPE(gsl_matrix) * m,
const size_t i,
const TYPE(gsl_vector) * v)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (i >= M)
	{
		GSL_ERROR("row index is out of range", GSL_EINVAL);
	}

	if (v->size != N)
	{
		GSL_ERROR("matrix row size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  const ATOMIC *v_data = v->data;
	  ATOMIC *row_data = m->data + MULTIPLICITY * i * tda;
	  const size_t stride = v->stride;
	  size_t j;

	  for (j = 0; j < N; j++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  row_data[MULTIPLICITY*j + k]
				  = v_data[MULTIPLICITY * stride * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, set_col) (TYPE(gsl_matrix) * m,
const size_t j,
const TYPE(gsl_vector) * v)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (j >= N)
	{
		GSL_ERROR("column index is out of range", GSL_EINVAL);
	}

	if (v->size != M)
	{
		GSL_ERROR("matrix column size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  const ATOMIC *v_data = v->data;
	  ATOMIC *column_data = m->data + MULTIPLICITY * j;
	  const size_t stride = v->stride;
	  size_t i;

	  for (i = 0; i < M; i++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  column_data[MULTIPLICITY * i * tda + k]
				  = v_data[MULTIPLICITY * stride * i + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}


TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_row_from_matrix) (TYPE(gsl_matrix) * m,
const size_t i)
{
	TYPE(gsl_vector) * v;

	const size_t M = m->size1;

	if (i >= M)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = m->data + MULTIPLICITY * i * m->tda;
	v->size = m->size2;
	v->stride = 1;
	v->block = 0;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_col_from_matrix) (TYPE(gsl_matrix) * m,
const size_t j)
{
	TYPE(gsl_vector) * v;

	const size_t N = m->size2;

	if (j >= N)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = m->data + MULTIPLICITY * j;
	v->size = m->size1;
	v->stride = m->tda;
	v->block = 0;

	return v;
}


#include "templates_off.h"
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
#include "templates_on.h"
int
FUNCTION(gsl_matrix, get_row) (TYPE(gsl_vector) * v,
const TYPE(gsl_matrix) * m,
const size_t i)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (i >= M)
	{
		GSL_ERROR("row index is out of range", GSL_EINVAL);
	}

	if (v->size != N)
	{
		GSL_ERROR("matrix row size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  ATOMIC *v_data = v->data;
	  const ATOMIC *row_data = m->data + MULTIPLICITY * i * tda;
	  const size_t stride = v->stride;
	  size_t j;

	  for (j = 0; j < N; j++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  v_data[MULTIPLICITY * stride * j + k]
				  = row_data[MULTIPLICITY * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, get_col) (TYPE(gsl_vector) * v,
const TYPE(gsl_matrix) * m,
const size_t j)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (j >= N)
	{
		GSL_ERROR("column index is out of range", GSL_EINVAL);
	}

	if (v->size != M)
	{
		GSL_ERROR("matrix column size and vector length are not equal",
			GSL_EBADLEN);
	}


  {
	  ATOMIC *v_data = v->data;
	  const ATOMIC *column_data = m->data + MULTIPLICITY * j;
	  const size_t stride = v->stride;
	  size_t i;

	  for (i = 0; i < M; i++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  v_data[stride * MULTIPLICITY * i + k] =
				  column_data[MULTIPLICITY * i * tda + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, set_row) (TYPE(gsl_matrix) * m,
const size_t i,
const TYPE(gsl_vector) * v)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (i >= M)
	{
		GSL_ERROR("row index is out of range", GSL_EINVAL);
	}

	if (v->size != N)
	{
		GSL_ERROR("matrix row size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  const ATOMIC *v_data = v->data;
	  ATOMIC *row_data = m->data + MULTIPLICITY * i * tda;
	  const size_t stride = v->stride;
	  size_t j;

	  for (j = 0; j < N; j++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  row_data[MULTIPLICITY*j + k]
				  = v_data[MULTIPLICITY * stride * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, set_col) (TYPE(gsl_matrix) * m,
const size_t j,
const TYPE(gsl_vector) * v)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (j >= N)
	{
		GSL_ERROR("column index is out of range", GSL_EINVAL);
	}

	if (v->size != M)
	{
		GSL_ERROR("matrix column size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  const ATOMIC *v_data = v->data;
	  ATOMIC *column_data = m->data + MULTIPLICITY * j;
	  const size_t stride = v->stride;
	  size_t i;

	  for (i = 0; i < M; i++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  column_data[MULTIPLICITY * i * tda + k]
				  = v_data[MULTIPLICITY * stride * i + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}


TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_row_from_matrix) (TYPE(gsl_matrix) * m,
const size_t i)
{
	TYPE(gsl_vector) * v;

	const size_t M = m->size1;

	if (i >= M)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = m->data + MULTIPLICITY * i * m->tda;
	v->size = m->size2;
	v->stride = 1;
	v->block = 0;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_col_from_matrix) (TYPE(gsl_matrix) * m,
const size_t j)
{
	TYPE(gsl_vector) * v;

	const size_t N = m->size2;

	if (j >= N)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = m->data + MULTIPLICITY * j;
	v->size = m->size1;
	v->stride = m->tda;
	v->block = 0;

	return v;
}


#include "templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "templates_on.h"
int
FUNCTION(gsl_matrix, get_row) (TYPE(gsl_vector) * v,
const TYPE(gsl_matrix) * m,
const size_t i)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (i >= M)
	{
		GSL_ERROR("row index is out of range", GSL_EINVAL);
	}

	if (v->size != N)
	{
		GSL_ERROR("matrix row size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  ATOMIC *v_data = v->data;
	  const ATOMIC *row_data = m->data + MULTIPLICITY * i * tda;
	  const size_t stride = v->stride;
	  size_t j;

	  for (j = 0; j < N; j++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  v_data[MULTIPLICITY * stride * j + k]
				  = row_data[MULTIPLICITY * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, get_col) (TYPE(gsl_vector) * v,
const TYPE(gsl_matrix) * m,
const size_t j)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (j >= N)
	{
		GSL_ERROR("column index is out of range", GSL_EINVAL);
	}

	if (v->size != M)
	{
		GSL_ERROR("matrix column size and vector length are not equal",
			GSL_EBADLEN);
	}


  {
	  ATOMIC *v_data = v->data;
	  const ATOMIC *column_data = m->data + MULTIPLICITY * j;
	  const size_t stride = v->stride;
	  size_t i;

	  for (i = 0; i < M; i++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  v_data[stride * MULTIPLICITY * i + k] =
				  column_data[MULTIPLICITY * i * tda + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, set_row) (TYPE(gsl_matrix) * m,
const size_t i,
const TYPE(gsl_vector) * v)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (i >= M)
	{
		GSL_ERROR("row index is out of range", GSL_EINVAL);
	}

	if (v->size != N)
	{
		GSL_ERROR("matrix row size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  const ATOMIC *v_data = v->data;
	  ATOMIC *row_data = m->data + MULTIPLICITY * i * tda;
	  const size_t stride = v->stride;
	  size_t j;

	  for (j = 0; j < N; j++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  row_data[MULTIPLICITY*j + k]
				  = v_data[MULTIPLICITY * stride * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, set_col) (TYPE(gsl_matrix) * m,
const size_t j,
const TYPE(gsl_vector) * v)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (j >= N)
	{
		GSL_ERROR("column index is out of range", GSL_EINVAL);
	}

	if (v->size != M)
	{
		GSL_ERROR("matrix column size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  const ATOMIC *v_data = v->data;
	  ATOMIC *column_data = m->data + MULTIPLICITY * j;
	  const size_t stride = v->stride;
	  size_t i;

	  for (i = 0; i < M; i++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  column_data[MULTIPLICITY * i * tda + k]
				  = v_data[MULTIPLICITY * stride * i + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}


TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_row_from_matrix) (TYPE(gsl_matrix) * m,
const size_t i)
{
	TYPE(gsl_vector) * v;

	const size_t M = m->size1;

	if (i >= M)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = m->data + MULTIPLICITY * i * m->tda;
	v->size = m->size2;
	v->stride = 1;
	v->block = 0;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_col_from_matrix) (TYPE(gsl_matrix) * m,
const size_t j)
{
	TYPE(gsl_vector) * v;

	const size_t N = m->size2;

	if (j >= N)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = m->data + MULTIPLICITY * j;
	v->size = m->size1;
	v->stride = m->tda;
	v->block = 0;

	return v;
}


#include "templates_off.h"
#undef  BASE_FLOAT

#define BASE_ULONG
#include "templates_on.h"
int
FUNCTION(gsl_matrix, get_row) (TYPE(gsl_vector) * v,
const TYPE(gsl_matrix) * m,
const size_t i)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (i >= M)
	{
		GSL_ERROR("row index is out of range", GSL_EINVAL);
	}

	if (v->size != N)
	{
		GSL_ERROR("matrix row size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  ATOMIC *v_data = v->data;
	  const ATOMIC *row_data = m->data + MULTIPLICITY * i * tda;
	  const size_t stride = v->stride;
	  size_t j;

	  for (j = 0; j < N; j++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  v_data[MULTIPLICITY * stride * j + k]
				  = row_data[MULTIPLICITY * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, get_col) (TYPE(gsl_vector) * v,
const TYPE(gsl_matrix) * m,
const size_t j)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (j >= N)
	{
		GSL_ERROR("column index is out of range", GSL_EINVAL);
	}

	if (v->size != M)
	{
		GSL_ERROR("matrix column size and vector length are not equal",
			GSL_EBADLEN);
	}


  {
	  ATOMIC *v_data = v->data;
	  const ATOMIC *column_data = m->data + MULTIPLICITY * j;
	  const size_t stride = v->stride;
	  size_t i;

	  for (i = 0; i < M; i++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  v_data[stride * MULTIPLICITY * i + k] =
				  column_data[MULTIPLICITY * i * tda + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, set_row) (TYPE(gsl_matrix) * m,
const size_t i,
const TYPE(gsl_vector) * v)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (i >= M)
	{
		GSL_ERROR("row index is out of range", GSL_EINVAL);
	}

	if (v->size != N)
	{
		GSL_ERROR("matrix row size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  const ATOMIC *v_data = v->data;
	  ATOMIC *row_data = m->data + MULTIPLICITY * i * tda;
	  const size_t stride = v->stride;
	  size_t j;

	  for (j = 0; j < N; j++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  row_data[MULTIPLICITY*j + k]
				  = v_data[MULTIPLICITY * stride * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, set_col) (TYPE(gsl_matrix) * m,
const size_t j,
const TYPE(gsl_vector) * v)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (j >= N)
	{
		GSL_ERROR("column index is out of range", GSL_EINVAL);
	}

	if (v->size != M)
	{
		GSL_ERROR("matrix column size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  const ATOMIC *v_data = v->data;
	  ATOMIC *column_data = m->data + MULTIPLICITY * j;
	  const size_t stride = v->stride;
	  size_t i;

	  for (i = 0; i < M; i++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  column_data[MULTIPLICITY * i * tda + k]
				  = v_data[MULTIPLICITY * stride * i + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}


TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_row_from_matrix) (TYPE(gsl_matrix) * m,
const size_t i)
{
	TYPE(gsl_vector) * v;

	const size_t M = m->size1;

	if (i >= M)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = m->data + MULTIPLICITY * i * m->tda;
	v->size = m->size2;
	v->stride = 1;
	v->block = 0;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_col_from_matrix) (TYPE(gsl_matrix) * m,
const size_t j)
{
	TYPE(gsl_vector) * v;

	const size_t N = m->size2;

	if (j >= N)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = m->data + MULTIPLICITY * j;
	v->size = m->size1;
	v->stride = m->tda;
	v->block = 0;

	return v;
}


#include "templates_off.h"
#undef  BASE_ULONG

#define BASE_LONG
#include "templates_on.h"
int
FUNCTION(gsl_matrix, get_row) (TYPE(gsl_vector) * v,
const TYPE(gsl_matrix) * m,
const size_t i)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (i >= M)
	{
		GSL_ERROR("row index is out of range", GSL_EINVAL);
	}

	if (v->size != N)
	{
		GSL_ERROR("matrix row size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  ATOMIC *v_data = v->data;
	  const ATOMIC *row_data = m->data + MULTIPLICITY * i * tda;
	  const size_t stride = v->stride;
	  size_t j;

	  for (j = 0; j < N; j++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  v_data[MULTIPLICITY * stride * j + k]
				  = row_data[MULTIPLICITY * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, get_col) (TYPE(gsl_vector) * v,
const TYPE(gsl_matrix) * m,
const size_t j)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (j >= N)
	{
		GSL_ERROR("column index is out of range", GSL_EINVAL);
	}

	if (v->size != M)
	{
		GSL_ERROR("matrix column size and vector length are not equal",
			GSL_EBADLEN);
	}


  {
	  ATOMIC *v_data = v->data;
	  const ATOMIC *column_data = m->data + MULTIPLICITY * j;
	  const size_t stride = v->stride;
	  size_t i;

	  for (i = 0; i < M; i++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  v_data[stride * MULTIPLICITY * i + k] =
				  column_data[MULTIPLICITY * i * tda + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, set_row) (TYPE(gsl_matrix) * m,
const size_t i,
const TYPE(gsl_vector) * v)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (i >= M)
	{
		GSL_ERROR("row index is out of range", GSL_EINVAL);
	}

	if (v->size != N)
	{
		GSL_ERROR("matrix row size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  const ATOMIC *v_data = v->data;
	  ATOMIC *row_data = m->data + MULTIPLICITY * i * tda;
	  const size_t stride = v->stride;
	  size_t j;

	  for (j = 0; j < N; j++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  row_data[MULTIPLICITY*j + k]
				  = v_data[MULTIPLICITY * stride * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, set_col) (TYPE(gsl_matrix) * m,
const size_t j,
const TYPE(gsl_vector) * v)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (j >= N)
	{
		GSL_ERROR("column index is out of range", GSL_EINVAL);
	}

	if (v->size != M)
	{
		GSL_ERROR("matrix column size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  const ATOMIC *v_data = v->data;
	  ATOMIC *column_data = m->data + MULTIPLICITY * j;
	  const size_t stride = v->stride;
	  size_t i;

	  for (i = 0; i < M; i++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  column_data[MULTIPLICITY * i * tda + k]
				  = v_data[MULTIPLICITY * stride * i + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}


TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_row_from_matrix) (TYPE(gsl_matrix) * m,
const size_t i)
{
	TYPE(gsl_vector) * v;

	const size_t M = m->size1;

	if (i >= M)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = m->data + MULTIPLICITY * i * m->tda;
	v->size = m->size2;
	v->stride = 1;
	v->block = 0;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_col_from_matrix) (TYPE(gsl_matrix) * m,
const size_t j)
{
	TYPE(gsl_vector) * v;

	const size_t N = m->size2;

	if (j >= N)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = m->data + MULTIPLICITY * j;
	v->size = m->size1;
	v->stride = m->tda;
	v->block = 0;

	return v;
}


#include "templates_off.h"
#undef  BASE_LONG

#define BASE_UINT
#include "templates_on.h"
int
FUNCTION(gsl_matrix, get_row) (TYPE(gsl_vector) * v,
const TYPE(gsl_matrix) * m,
const size_t i)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (i >= M)
	{
		GSL_ERROR("row index is out of range", GSL_EINVAL);
	}

	if (v->size != N)
	{
		GSL_ERROR("matrix row size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  ATOMIC *v_data = v->data;
	  const ATOMIC *row_data = m->data + MULTIPLICITY * i * tda;
	  const size_t stride = v->stride;
	  size_t j;

	  for (j = 0; j < N; j++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  v_data[MULTIPLICITY * stride * j + k]
				  = row_data[MULTIPLICITY * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, get_col) (TYPE(gsl_vector) * v,
const TYPE(gsl_matrix) * m,
const size_t j)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (j >= N)
	{
		GSL_ERROR("column index is out of range", GSL_EINVAL);
	}

	if (v->size != M)
	{
		GSL_ERROR("matrix column size and vector length are not equal",
			GSL_EBADLEN);
	}


  {
	  ATOMIC *v_data = v->data;
	  const ATOMIC *column_data = m->data + MULTIPLICITY * j;
	  const size_t stride = v->stride;
	  size_t i;

	  for (i = 0; i < M; i++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  v_data[stride * MULTIPLICITY * i + k] =
				  column_data[MULTIPLICITY * i * tda + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, set_row) (TYPE(gsl_matrix) * m,
const size_t i,
const TYPE(gsl_vector) * v)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (i >= M)
	{
		GSL_ERROR("row index is out of range", GSL_EINVAL);
	}

	if (v->size != N)
	{
		GSL_ERROR("matrix row size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  const ATOMIC *v_data = v->data;
	  ATOMIC *row_data = m->data + MULTIPLICITY * i * tda;
	  const size_t stride = v->stride;
	  size_t j;

	  for (j = 0; j < N; j++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  row_data[MULTIPLICITY*j + k]
				  = v_data[MULTIPLICITY * stride * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, set_col) (TYPE(gsl_matrix) * m,
const size_t j,
const TYPE(gsl_vector) * v)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (j >= N)
	{
		GSL_ERROR("column index is out of range", GSL_EINVAL);
	}

	if (v->size != M)
	{
		GSL_ERROR("matrix column size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  const ATOMIC *v_data = v->data;
	  ATOMIC *column_data = m->data + MULTIPLICITY * j;
	  const size_t stride = v->stride;
	  size_t i;

	  for (i = 0; i < M; i++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  column_data[MULTIPLICITY * i * tda + k]
				  = v_data[MULTIPLICITY * stride * i + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}


TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_row_from_matrix) (TYPE(gsl_matrix) * m,
const size_t i)
{
	TYPE(gsl_vector) * v;

	const size_t M = m->size1;

	if (i >= M)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = m->data + MULTIPLICITY * i * m->tda;
	v->size = m->size2;
	v->stride = 1;
	v->block = 0;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_col_from_matrix) (TYPE(gsl_matrix) * m,
const size_t j)
{
	TYPE(gsl_vector) * v;

	const size_t N = m->size2;

	if (j >= N)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = m->data + MULTIPLICITY * j;
	v->size = m->size1;
	v->stride = m->tda;
	v->block = 0;

	return v;
}


#include "templates_off.h"
#undef  BASE_UINT

#define BASE_INT
#include "templates_on.h"
int
FUNCTION(gsl_matrix, get_row) (TYPE(gsl_vector) * v,
const TYPE(gsl_matrix) * m,
const size_t i)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (i >= M)
	{
		GSL_ERROR("row index is out of range", GSL_EINVAL);
	}

	if (v->size != N)
	{
		GSL_ERROR("matrix row size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  ATOMIC *v_data = v->data;
	  const ATOMIC *row_data = m->data + MULTIPLICITY * i * tda;
	  const size_t stride = v->stride;
	  size_t j;

	  for (j = 0; j < N; j++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  v_data[MULTIPLICITY * stride * j + k]
				  = row_data[MULTIPLICITY * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, get_col) (TYPE(gsl_vector) * v,
const TYPE(gsl_matrix) * m,
const size_t j)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (j >= N)
	{
		GSL_ERROR("column index is out of range", GSL_EINVAL);
	}

	if (v->size != M)
	{
		GSL_ERROR("matrix column size and vector length are not equal",
			GSL_EBADLEN);
	}


  {
	  ATOMIC *v_data = v->data;
	  const ATOMIC *column_data = m->data + MULTIPLICITY * j;
	  const size_t stride = v->stride;
	  size_t i;

	  for (i = 0; i < M; i++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  v_data[stride * MULTIPLICITY * i + k] =
				  column_data[MULTIPLICITY * i * tda + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, set_row) (TYPE(gsl_matrix) * m,
const size_t i,
const TYPE(gsl_vector) * v)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (i >= M)
	{
		GSL_ERROR("row index is out of range", GSL_EINVAL);
	}

	if (v->size != N)
	{
		GSL_ERROR("matrix row size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  const ATOMIC *v_data = v->data;
	  ATOMIC *row_data = m->data + MULTIPLICITY * i * tda;
	  const size_t stride = v->stride;
	  size_t j;

	  for (j = 0; j < N; j++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  row_data[MULTIPLICITY*j + k]
				  = v_data[MULTIPLICITY * stride * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, set_col) (TYPE(gsl_matrix) * m,
const size_t j,
const TYPE(gsl_vector) * v)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (j >= N)
	{
		GSL_ERROR("column index is out of range", GSL_EINVAL);
	}

	if (v->size != M)
	{
		GSL_ERROR("matrix column size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  const ATOMIC *v_data = v->data;
	  ATOMIC *column_data = m->data + MULTIPLICITY * j;
	  const size_t stride = v->stride;
	  size_t i;

	  for (i = 0; i < M; i++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  column_data[MULTIPLICITY * i * tda + k]
				  = v_data[MULTIPLICITY * stride * i + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}


TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_row_from_matrix) (TYPE(gsl_matrix) * m,
const size_t i)
{
	TYPE(gsl_vector) * v;

	const size_t M = m->size1;

	if (i >= M)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = m->data + MULTIPLICITY * i * m->tda;
	v->size = m->size2;
	v->stride = 1;
	v->block = 0;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_col_from_matrix) (TYPE(gsl_matrix) * m,
const size_t j)
{
	TYPE(gsl_vector) * v;

	const size_t N = m->size2;

	if (j >= N)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = m->data + MULTIPLICITY * j;
	v->size = m->size1;
	v->stride = m->tda;
	v->block = 0;

	return v;
}


#include "templates_off.h"
#undef  BASE_INT

#define BASE_USHORT
#include "templates_on.h"
int
FUNCTION(gsl_matrix, get_row) (TYPE(gsl_vector) * v,
const TYPE(gsl_matrix) * m,
const size_t i)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (i >= M)
	{
		GSL_ERROR("row index is out of range", GSL_EINVAL);
	}

	if (v->size != N)
	{
		GSL_ERROR("matrix row size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  ATOMIC *v_data = v->data;
	  const ATOMIC *row_data = m->data + MULTIPLICITY * i * tda;
	  const size_t stride = v->stride;
	  size_t j;

	  for (j = 0; j < N; j++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  v_data[MULTIPLICITY * stride * j + k]
				  = row_data[MULTIPLICITY * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, get_col) (TYPE(gsl_vector) * v,
const TYPE(gsl_matrix) * m,
const size_t j)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (j >= N)
	{
		GSL_ERROR("column index is out of range", GSL_EINVAL);
	}

	if (v->size != M)
	{
		GSL_ERROR("matrix column size and vector length are not equal",
			GSL_EBADLEN);
	}


  {
	  ATOMIC *v_data = v->data;
	  const ATOMIC *column_data = m->data + MULTIPLICITY * j;
	  const size_t stride = v->stride;
	  size_t i;

	  for (i = 0; i < M; i++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  v_data[stride * MULTIPLICITY * i + k] =
				  column_data[MULTIPLICITY * i * tda + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, set_row) (TYPE(gsl_matrix) * m,
const size_t i,
const TYPE(gsl_vector) * v)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (i >= M)
	{
		GSL_ERROR("row index is out of range", GSL_EINVAL);
	}

	if (v->size != N)
	{
		GSL_ERROR("matrix row size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  const ATOMIC *v_data = v->data;
	  ATOMIC *row_data = m->data + MULTIPLICITY * i * tda;
	  const size_t stride = v->stride;
	  size_t j;

	  for (j = 0; j < N; j++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  row_data[MULTIPLICITY*j + k]
				  = v_data[MULTIPLICITY * stride * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, set_col) (TYPE(gsl_matrix) * m,
const size_t j,
const TYPE(gsl_vector) * v)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (j >= N)
	{
		GSL_ERROR("column index is out of range", GSL_EINVAL);
	}

	if (v->size != M)
	{
		GSL_ERROR("matrix column size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  const ATOMIC *v_data = v->data;
	  ATOMIC *column_data = m->data + MULTIPLICITY * j;
	  const size_t stride = v->stride;
	  size_t i;

	  for (i = 0; i < M; i++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  column_data[MULTIPLICITY * i * tda + k]
				  = v_data[MULTIPLICITY * stride * i + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}


TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_row_from_matrix) (TYPE(gsl_matrix) * m,
const size_t i)
{
	TYPE(gsl_vector) * v;

	const size_t M = m->size1;

	if (i >= M)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = m->data + MULTIPLICITY * i * m->tda;
	v->size = m->size2;
	v->stride = 1;
	v->block = 0;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_col_from_matrix) (TYPE(gsl_matrix) * m,
const size_t j)
{
	TYPE(gsl_vector) * v;

	const size_t N = m->size2;

	if (j >= N)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = m->data + MULTIPLICITY * j;
	v->size = m->size1;
	v->stride = m->tda;
	v->block = 0;

	return v;
}


#include "templates_off.h"
#undef  BASE_USHORT

#define BASE_SHORT
#include "templates_on.h"
int
FUNCTION(gsl_matrix, get_row) (TYPE(gsl_vector) * v,
const TYPE(gsl_matrix) * m,
const size_t i)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (i >= M)
	{
		GSL_ERROR("row index is out of range", GSL_EINVAL);
	}

	if (v->size != N)
	{
		GSL_ERROR("matrix row size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  ATOMIC *v_data = v->data;
	  const ATOMIC *row_data = m->data + MULTIPLICITY * i * tda;
	  const size_t stride = v->stride;
	  size_t j;

	  for (j = 0; j < N; j++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  v_data[MULTIPLICITY * stride * j + k]
				  = row_data[MULTIPLICITY * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, get_col) (TYPE(gsl_vector) * v,
const TYPE(gsl_matrix) * m,
const size_t j)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (j >= N)
	{
		GSL_ERROR("column index is out of range", GSL_EINVAL);
	}

	if (v->size != M)
	{
		GSL_ERROR("matrix column size and vector length are not equal",
			GSL_EBADLEN);
	}


  {
	  ATOMIC *v_data = v->data;
	  const ATOMIC *column_data = m->data + MULTIPLICITY * j;
	  const size_t stride = v->stride;
	  size_t i;

	  for (i = 0; i < M; i++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  v_data[stride * MULTIPLICITY * i + k] =
				  column_data[MULTIPLICITY * i * tda + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, set_row) (TYPE(gsl_matrix) * m,
const size_t i,
const TYPE(gsl_vector) * v)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (i >= M)
	{
		GSL_ERROR("row index is out of range", GSL_EINVAL);
	}

	if (v->size != N)
	{
		GSL_ERROR("matrix row size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  const ATOMIC *v_data = v->data;
	  ATOMIC *row_data = m->data + MULTIPLICITY * i * tda;
	  const size_t stride = v->stride;
	  size_t j;

	  for (j = 0; j < N; j++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  row_data[MULTIPLICITY*j + k]
				  = v_data[MULTIPLICITY * stride * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, set_col) (TYPE(gsl_matrix) * m,
const size_t j,
const TYPE(gsl_vector) * v)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (j >= N)
	{
		GSL_ERROR("column index is out of range", GSL_EINVAL);
	}

	if (v->size != M)
	{
		GSL_ERROR("matrix column size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  const ATOMIC *v_data = v->data;
	  ATOMIC *column_data = m->data + MULTIPLICITY * j;
	  const size_t stride = v->stride;
	  size_t i;

	  for (i = 0; i < M; i++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  column_data[MULTIPLICITY * i * tda + k]
				  = v_data[MULTIPLICITY * stride * i + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}


TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_row_from_matrix) (TYPE(gsl_matrix) * m,
const size_t i)
{
	TYPE(gsl_vector) * v;

	const size_t M = m->size1;

	if (i >= M)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = m->data + MULTIPLICITY * i * m->tda;
	v->size = m->size2;
	v->stride = 1;
	v->block = 0;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_col_from_matrix) (TYPE(gsl_matrix) * m,
const size_t j)
{
	TYPE(gsl_vector) * v;

	const size_t N = m->size2;

	if (j >= N)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = m->data + MULTIPLICITY * j;
	v->size = m->size1;
	v->stride = m->tda;
	v->block = 0;

	return v;
}


#include "templates_off.h"
#undef  BASE_SHORT

#define BASE_UCHAR
#include "templates_on.h"
int
FUNCTION(gsl_matrix, get_row) (TYPE(gsl_vector) * v,
const TYPE(gsl_matrix) * m,
const size_t i)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (i >= M)
	{
		GSL_ERROR("row index is out of range", GSL_EINVAL);
	}

	if (v->size != N)
	{
		GSL_ERROR("matrix row size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  ATOMIC *v_data = v->data;
	  const ATOMIC *row_data = m->data + MULTIPLICITY * i * tda;
	  const size_t stride = v->stride;
	  size_t j;

	  for (j = 0; j < N; j++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  v_data[MULTIPLICITY * stride * j + k]
				  = row_data[MULTIPLICITY * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, get_col) (TYPE(gsl_vector) * v,
const TYPE(gsl_matrix) * m,
const size_t j)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (j >= N)
	{
		GSL_ERROR("column index is out of range", GSL_EINVAL);
	}

	if (v->size != M)
	{
		GSL_ERROR("matrix column size and vector length are not equal",
			GSL_EBADLEN);
	}


  {
	  ATOMIC *v_data = v->data;
	  const ATOMIC *column_data = m->data + MULTIPLICITY * j;
	  const size_t stride = v->stride;
	  size_t i;

	  for (i = 0; i < M; i++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  v_data[stride * MULTIPLICITY * i + k] =
				  column_data[MULTIPLICITY * i * tda + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, set_row) (TYPE(gsl_matrix) * m,
const size_t i,
const TYPE(gsl_vector) * v)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (i >= M)
	{
		GSL_ERROR("row index is out of range", GSL_EINVAL);
	}

	if (v->size != N)
	{
		GSL_ERROR("matrix row size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  const ATOMIC *v_data = v->data;
	  ATOMIC *row_data = m->data + MULTIPLICITY * i * tda;
	  const size_t stride = v->stride;
	  size_t j;

	  for (j = 0; j < N; j++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  row_data[MULTIPLICITY*j + k]
				  = v_data[MULTIPLICITY * stride * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, set_col) (TYPE(gsl_matrix) * m,
const size_t j,
const TYPE(gsl_vector) * v)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (j >= N)
	{
		GSL_ERROR("column index is out of range", GSL_EINVAL);
	}

	if (v->size != M)
	{
		GSL_ERROR("matrix column size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  const ATOMIC *v_data = v->data;
	  ATOMIC *column_data = m->data + MULTIPLICITY * j;
	  const size_t stride = v->stride;
	  size_t i;

	  for (i = 0; i < M; i++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  column_data[MULTIPLICITY * i * tda + k]
				  = v_data[MULTIPLICITY * stride * i + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}


TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_row_from_matrix) (TYPE(gsl_matrix) * m,
const size_t i)
{
	TYPE(gsl_vector) * v;

	const size_t M = m->size1;

	if (i >= M)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = m->data + MULTIPLICITY * i * m->tda;
	v->size = m->size2;
	v->stride = 1;
	v->block = 0;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_col_from_matrix) (TYPE(gsl_matrix) * m,
const size_t j)
{
	TYPE(gsl_vector) * v;

	const size_t N = m->size2;

	if (j >= N)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = m->data + MULTIPLICITY * j;
	v->size = m->size1;
	v->stride = m->tda;
	v->block = 0;

	return v;
}


#include "templates_off.h"
#undef  BASE_UCHAR

#define BASE_CHAR
#include "templates_on.h"
int
FUNCTION(gsl_matrix, get_row) (TYPE(gsl_vector) * v,
const TYPE(gsl_matrix) * m,
const size_t i)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (i >= M)
	{
		GSL_ERROR("row index is out of range", GSL_EINVAL);
	}

	if (v->size != N)
	{
		GSL_ERROR("matrix row size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  ATOMIC *v_data = v->data;
	  const ATOMIC *row_data = m->data + MULTIPLICITY * i * tda;
	  const size_t stride = v->stride;
	  size_t j;

	  for (j = 0; j < N; j++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  v_data[MULTIPLICITY * stride * j + k]
				  = row_data[MULTIPLICITY * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, get_col) (TYPE(gsl_vector) * v,
const TYPE(gsl_matrix) * m,
const size_t j)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (j >= N)
	{
		GSL_ERROR("column index is out of range", GSL_EINVAL);
	}

	if (v->size != M)
	{
		GSL_ERROR("matrix column size and vector length are not equal",
			GSL_EBADLEN);
	}


  {
	  ATOMIC *v_data = v->data;
	  const ATOMIC *column_data = m->data + MULTIPLICITY * j;
	  const size_t stride = v->stride;
	  size_t i;

	  for (i = 0; i < M; i++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  v_data[stride * MULTIPLICITY * i + k] =
				  column_data[MULTIPLICITY * i * tda + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, set_row) (TYPE(gsl_matrix) * m,
const size_t i,
const TYPE(gsl_vector) * v)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (i >= M)
	{
		GSL_ERROR("row index is out of range", GSL_EINVAL);
	}

	if (v->size != N)
	{
		GSL_ERROR("matrix row size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  const ATOMIC *v_data = v->data;
	  ATOMIC *row_data = m->data + MULTIPLICITY * i * tda;
	  const size_t stride = v->stride;
	  size_t j;

	  for (j = 0; j < N; j++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  row_data[MULTIPLICITY*j + k]
				  = v_data[MULTIPLICITY * stride * j + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, set_col) (TYPE(gsl_matrix) * m,
const size_t j,
const TYPE(gsl_vector) * v)
{
	const size_t M = m->size1;
	const size_t N = m->size2;
	const size_t tda = m->tda;

	if (j >= N)
	{
		GSL_ERROR("column index is out of range", GSL_EINVAL);
	}

	if (v->size != M)
	{
		GSL_ERROR("matrix column size and vector length are not equal",
			GSL_EBADLEN);
	}

  {
	  const ATOMIC *v_data = v->data;
	  ATOMIC *column_data = m->data + MULTIPLICITY * j;
	  const size_t stride = v->stride;
	  size_t i;

	  for (i = 0; i < M; i++)
	  {
		  unsigned int k;

		  for (k = 0; k < MULTIPLICITY; k++)
		  {
			  column_data[MULTIPLICITY * i * tda + k]
				  = v_data[MULTIPLICITY * stride * i + k];
		  }
	  }
  }

	return GSL_SUCCESS;
}


TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_row_from_matrix) (TYPE(gsl_matrix) * m,
const size_t i)
{
	TYPE(gsl_vector) * v;

	const size_t M = m->size1;

	if (i >= M)
	{
		GSL_ERROR_VAL("row index is out of range", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = m->data + MULTIPLICITY * i * m->tda;
	v->size = m->size2;
	v->stride = 1;
	v->block = 0;

	return v;
}

TYPE(gsl_vector) *
FUNCTION(gsl_vector, alloc_col_from_matrix) (TYPE(gsl_matrix) * m,
const size_t j)
{
	TYPE(gsl_vector) * v;

	const size_t N = m->size2;

	if (j >= N)
	{
		GSL_ERROR_VAL("column index is out of range", GSL_EINVAL, 0);
	}

	v = (TYPE(gsl_vector) *) malloc(sizeof(TYPE(gsl_vector)));

	if (v == 0)
	{
		GSL_ERROR_VAL("failed to allocate space for vector struct",
			GSL_ENOMEM, 0);
	}

	v->data = m->data + MULTIPLICITY * j;
	v->size = m->size1;
	v->stride = m->tda;
	v->block = 0;

	return v;
}


#include "templates_off.h"
#undef  BASE_CHAR
