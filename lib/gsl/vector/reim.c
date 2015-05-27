#include "stdafx.h"
#include "config.h.in"
#include <stdlib.h>
#include <gsl/vector/gsl_vector.h>

#include "view.h"

#define BASE_GSL_COMPLEX_LONG
#include "../templates_on.h"
QUALIFIED_REAL_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, real) (QUALIFIED_TYPE(gsl_vector) * v)
{
	REAL_TYPE(gsl_vector) s = NULL_VECTOR;

	s.data = v->data;
	s.size = v->size;
	s.stride = MULTIPLICITY * v->stride;
	s.block = 0;  /* FIXME: should be v->block, but cannot point to
				  block of different type */
	s.owner = 0;

	{
		QUALIFIED_REAL_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;
		view.vector = s;
		return view;
	}
}

QUALIFIED_REAL_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, imag) (QUALIFIED_TYPE(gsl_vector) * v)
{
	REAL_TYPE(gsl_vector) s = NULL_VECTOR;

	s.data = v->data + 1;
	s.size = v->size;
	s.stride = MULTIPLICITY * v->stride;
	s.block = 0;  /* FIXME: cannot point to block of different type */
	s.owner = 0;

	{
		QUALIFIED_REAL_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;
		view.vector = s;
		return view;
	}
}


#include "../templates_off.h"
#undef  BASE_GSL_COMPLEX_LONG

#define BASE_GSL_COMPLEX
#include "../templates_on.h"
QUALIFIED_REAL_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, real) (QUALIFIED_TYPE(gsl_vector) * v)
{
	REAL_TYPE(gsl_vector) s = NULL_VECTOR;

	s.data = v->data;
	s.size = v->size;
	s.stride = MULTIPLICITY * v->stride;
	s.block = 0;  /* FIXME: should be v->block, but cannot point to
				  block of different type */
	s.owner = 0;

	{
		QUALIFIED_REAL_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;
		view.vector = s;
		return view;
	}
}

QUALIFIED_REAL_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, imag) (QUALIFIED_TYPE(gsl_vector) * v)
{
	REAL_TYPE(gsl_vector) s = NULL_VECTOR;

	s.data = v->data + 1;
	s.size = v->size;
	s.stride = MULTIPLICITY * v->stride;
	s.block = 0;  /* FIXME: cannot point to block of different type */
	s.owner = 0;

	{
		QUALIFIED_REAL_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;
		view.vector = s;
		return view;
	}
}


#include "../templates_off.h"
#undef  BASE_GSL_COMPLEX

#define BASE_GSL_COMPLEX_FLOAT
#include "../templates_on.h"
QUALIFIED_REAL_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, real) (QUALIFIED_TYPE(gsl_vector) * v)
{
	REAL_TYPE(gsl_vector) s = NULL_VECTOR;

	s.data = v->data;
	s.size = v->size;
	s.stride = MULTIPLICITY * v->stride;
	s.block = 0;  /* FIXME: should be v->block, but cannot point to
				  block of different type */
	s.owner = 0;

	{
		QUALIFIED_REAL_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;
		view.vector = s;
		return view;
	}
}

QUALIFIED_REAL_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, imag) (QUALIFIED_TYPE(gsl_vector) * v)
{
	REAL_TYPE(gsl_vector) s = NULL_VECTOR;

	s.data = v->data + 1;
	s.size = v->size;
	s.stride = MULTIPLICITY * v->stride;
	s.block = 0;  /* FIXME: cannot point to block of different type */
	s.owner = 0;

	{
		QUALIFIED_REAL_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;
		view.vector = s;
		return view;
	}
}


#include "../templates_off.h"
#undef  BASE_GSL_COMPLEX_FLOAT

#define USE_QUALIFIER
#define QUALIFIER const

#define BASE_GSL_COMPLEX_LONG
#include "../templates_on.h"
QUALIFIED_REAL_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, real) (QUALIFIED_TYPE(gsl_vector) * v)
{
	REAL_TYPE(gsl_vector) s = NULL_VECTOR;

	s.data = v->data;
	s.size = v->size;
	s.stride = MULTIPLICITY * v->stride;
	s.block = 0;  /* FIXME: should be v->block, but cannot point to
				  block of different type */
	s.owner = 0;

	{
		QUALIFIED_REAL_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;
		view.vector = s;
		return view;
	}
}

QUALIFIED_REAL_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, imag) (QUALIFIED_TYPE(gsl_vector) * v)
{
	REAL_TYPE(gsl_vector) s = NULL_VECTOR;

	s.data = v->data + 1;
	s.size = v->size;
	s.stride = MULTIPLICITY * v->stride;
	s.block = 0;  /* FIXME: cannot point to block of different type */
	s.owner = 0;

	{
		QUALIFIED_REAL_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;
		view.vector = s;
		return view;
	}
}


#include "../templates_off.h"
#undef  BASE_GSL_COMPLEX_LONG

#define BASE_GSL_COMPLEX
#include "../templates_on.h"
QUALIFIED_REAL_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, real) (QUALIFIED_TYPE(gsl_vector) * v)
{
	REAL_TYPE(gsl_vector) s = NULL_VECTOR;

	s.data = v->data;
	s.size = v->size;
	s.stride = MULTIPLICITY * v->stride;
	s.block = 0;  /* FIXME: should be v->block, but cannot point to
				  block of different type */
	s.owner = 0;

	{
		QUALIFIED_REAL_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;
		view.vector = s;
		return view;
	}
}

QUALIFIED_REAL_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, imag) (QUALIFIED_TYPE(gsl_vector) * v)
{
	REAL_TYPE(gsl_vector) s = NULL_VECTOR;

	s.data = v->data + 1;
	s.size = v->size;
	s.stride = MULTIPLICITY * v->stride;
	s.block = 0;  /* FIXME: cannot point to block of different type */
	s.owner = 0;

	{
		QUALIFIED_REAL_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;
		view.vector = s;
		return view;
	}
}


#include "../templates_off.h"
#undef  BASE_GSL_COMPLEX

#define BASE_GSL_COMPLEX_FLOAT
#include "../templates_on.h"
QUALIFIED_REAL_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, real) (QUALIFIED_TYPE(gsl_vector) * v)
{
	REAL_TYPE(gsl_vector) s = NULL_VECTOR;

	s.data = v->data;
	s.size = v->size;
	s.stride = MULTIPLICITY * v->stride;
	s.block = 0;  /* FIXME: should be v->block, but cannot point to
				  block of different type */
	s.owner = 0;

	{
		QUALIFIED_REAL_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;
		view.vector = s;
		return view;
	}
}

QUALIFIED_REAL_VIEW(_gsl_vector, view)
FUNCTION(gsl_vector, imag) (QUALIFIED_TYPE(gsl_vector) * v)
{
	REAL_TYPE(gsl_vector) s = NULL_VECTOR;

	s.data = v->data + 1;
	s.size = v->size;
	s.stride = MULTIPLICITY * v->stride;
	s.block = 0;  /* FIXME: cannot point to block of different type */
	s.owner = 0;

	{
		QUALIFIED_REAL_VIEW(_gsl_vector, view) view = NULL_VECTOR_VIEW;
		view.vector = s;
		return view;
	}
}


#include "../templates_off.h"
#undef  BASE_GSL_COMPLEX_FLOAT
