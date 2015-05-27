#include "stdafx.h"
#include <config.h.in>
#include <gsl/matrix/gsl_matrix.h>
#include <gsl/err/gsl_errno.h>

#define BASE_GSL_COMPLEX_LONG
#include "templates_on.h"
int
FUNCTION(gsl_matrix, isnull) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] != 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


int
FUNCTION(gsl_matrix, ispos) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] <= 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


int
FUNCTION(gsl_matrix, isneg) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] >= 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


#include "templates_off.h"
#undef  BASE_GSL_COMPLEX_LONG

#define BASE_GSL_COMPLEX
#include "templates_on.h"
int
FUNCTION(gsl_matrix, isnull) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] != 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


int
FUNCTION(gsl_matrix, ispos) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] <= 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


int
FUNCTION(gsl_matrix, isneg) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] >= 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


#include "templates_off.h"
#undef  BASE_GSL_COMPLEX

#define BASE_GSL_COMPLEX_FLOAT
#include "templates_on.h"
int
FUNCTION(gsl_matrix, isnull) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] != 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


int
FUNCTION(gsl_matrix, ispos) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] <= 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


int
FUNCTION(gsl_matrix, isneg) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] >= 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


#include "templates_off.h"
#undef  BASE_GSL_COMPLEX_FLOAT

#define BASE_LONG_DOUBLE
#include "templates_on.h"
int
FUNCTION(gsl_matrix, isnull) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] != 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


int
FUNCTION(gsl_matrix, ispos) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] <= 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


int
FUNCTION(gsl_matrix, isneg) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] >= 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


#include "templates_off.h"
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
#include "templates_on.h"
int
FUNCTION(gsl_matrix, isnull) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] != 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


int
FUNCTION(gsl_matrix, ispos) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] <= 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


int
FUNCTION(gsl_matrix, isneg) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] >= 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


#include "templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "templates_on.h"
int
FUNCTION(gsl_matrix, isnull) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] != 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


int
FUNCTION(gsl_matrix, ispos) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] <= 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


int
FUNCTION(gsl_matrix, isneg) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] >= 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


#include "templates_off.h"
#undef  BASE_FLOAT

#define BASE_ULONG
#include "templates_on.h"
int
FUNCTION(gsl_matrix, isnull) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] != 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


int
FUNCTION(gsl_matrix, ispos) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] <= 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


int
FUNCTION(gsl_matrix, isneg) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] >= 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


#include "templates_off.h"
#undef  BASE_ULONG

#define BASE_LONG
#include "templates_on.h"
int
FUNCTION(gsl_matrix, isnull) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] != 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


int
FUNCTION(gsl_matrix, ispos) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] <= 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


int
FUNCTION(gsl_matrix, isneg) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] >= 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


#include "templates_off.h"
#undef  BASE_LONG

#define BASE_UINT
#include "templates_on.h"
int
FUNCTION(gsl_matrix, isnull) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] != 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


int
FUNCTION(gsl_matrix, ispos) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] <= 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


int
FUNCTION(gsl_matrix, isneg) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] >= 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


#include "templates_off.h"
#undef  BASE_UINT

#define BASE_INT
#include "templates_on.h"
int
FUNCTION(gsl_matrix, isnull) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] != 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


int
FUNCTION(gsl_matrix, ispos) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] <= 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


int
FUNCTION(gsl_matrix, isneg) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] >= 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


#include "templates_off.h"
#undef  BASE_INT

#define BASE_USHORT
#include "templates_on.h"
int
FUNCTION(gsl_matrix, isnull) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] != 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


int
FUNCTION(gsl_matrix, ispos) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] <= 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


int
FUNCTION(gsl_matrix, isneg) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] >= 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


#include "templates_off.h"
#undef  BASE_USHORT

#define BASE_SHORT
#include "templates_on.h"
int
FUNCTION(gsl_matrix, isnull) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] != 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


int
FUNCTION(gsl_matrix, ispos) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] <= 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


int
FUNCTION(gsl_matrix, isneg) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] >= 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


#include "templates_off.h"
#undef  BASE_SHORT

#define BASE_UCHAR
#include "templates_on.h"
int
FUNCTION(gsl_matrix, isnull) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] != 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


int
FUNCTION(gsl_matrix, ispos) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] <= 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


int
FUNCTION(gsl_matrix, isneg) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] >= 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


#include "templates_off.h"
#undef  BASE_UCHAR

#define BASE_CHAR
#include "templates_on.h"
int
FUNCTION(gsl_matrix, isnull) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] != 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


int
FUNCTION(gsl_matrix, ispos) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] <= 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


int
FUNCTION(gsl_matrix, isneg) (const TYPE(gsl_matrix) * m)
{
	const size_t size1 = m->size1;
	const size_t size2 = m->size2;
	const size_t tda = m->tda;

	size_t i, j, k;

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			for (k = 0; k < MULTIPLICITY; k++)
			{
				if (m->data[(i * tda + j) * MULTIPLICITY + k] >= 0.0)
				{
					return 0;
				}
			}
		}
	}

	return 1;
}


#include "templates_off.h"
#undef  BASE_CHAR
