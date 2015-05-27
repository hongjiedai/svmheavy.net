#include "stdafx.h"
#include <config.h.in>
#include <gsl/err/gsl_errno.h>
#include <gsl/vector/gsl_vector.h>
#include <gsl/permutation/gsl_permute.h>
#include <gsl/permutation/gsl_permute_vector.h>

#define BASE_GSL_COMPLEX_LONG
#include "templates_on.h"
int
TYPE(gsl_permute) (const size_t * p, ATOMIC * data, const size_t stride, const size_t n)
{
	size_t i, k, pk;

	for (i = 0; i < n; i++)
	{
		k = p[i];

		while (k > i)
			k = p[k];

		if (k < i)
			continue;

		/* Now have k == i, i.e the least in its cycle */

		pk = p[k];

		if (pk == i)
			continue;

		/* shuffle the elements of the cycle */

		{
			unsigned int a;

			ATOMIC t[MULTIPLICITY];

			for (a = 0; a < MULTIPLICITY; a++)
				t[a] = data[i*stride*MULTIPLICITY + a];

			while (pk != i)
			{
				for (a = 0; a < MULTIPLICITY; a++)
				{
					ATOMIC r1 = data[pk*stride*MULTIPLICITY + a];
					data[k*stride*MULTIPLICITY + a] = r1;
				}
				k = pk;
				pk = p[k];
			};

			for (a = 0; a < MULTIPLICITY; a++)
				data[k*stride*MULTIPLICITY + a] = t[a];
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_permute, inverse) (const size_t * p, ATOMIC * data, const size_t stride, const size_t n)
{
	size_t i, k, pk;

	for (i = 0; i < n; i++)
	{
		k = p[i];

		while (k > i)
			k = p[k];

		if (k < i)
			continue;

		/* Now have k == i, i.e the least in its cycle */

		pk = p[k];

		if (pk == i)
			continue;

		/* shuffle the elements of the cycle in the inverse direction */

		{
			unsigned int a;

			ATOMIC t[MULTIPLICITY];

			for (a = 0; a < MULTIPLICITY; a++)
				t[a] = data[k*stride*MULTIPLICITY + a];

			while (pk != i)
			{
				for (a = 0; a < MULTIPLICITY; a++)
				{
					ATOMIC r1 = data[pk*stride*MULTIPLICITY + a];
					data[pk*stride*MULTIPLICITY + a] = t[a];
					t[a] = r1;
				}

				k = pk;
				pk = p[k];
			};

			for (a = 0; a < MULTIPLICITY; a++)
				data[pk*stride*MULTIPLICITY + a] = t[a];
		}
	}

	return GSL_SUCCESS;
}


int
TYPE(gsl_permute_vector) (const gsl_permutation * p, TYPE(gsl_vector) * v)
{
	if (v->size != p->size)
	{
		GSL_ERROR("vector and permutation must be the same length", GSL_EBADLEN);
	}

	TYPE(gsl_permute) (p->data, v->data, v->stride, v->size);

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_permute_vector, inverse) (const gsl_permutation * p, TYPE(gsl_vector) * v)
{
	if (v->size != p->size)
	{
		GSL_ERROR("vector and permutation must be the same length", GSL_EBADLEN);
	}

	FUNCTION(gsl_permute, inverse) (p->data, v->data, v->stride, v->size);

	return GSL_SUCCESS;
}

#include "templates_off.h"
#undef  BASE_GSL_COMPLEX_LONG

#define BASE_GSL_COMPLEX
#include "templates_on.h"
int
TYPE(gsl_permute) (const size_t * p, ATOMIC * data, const size_t stride, const size_t n)
{
	size_t i, k, pk;

	for (i = 0; i < n; i++)
	{
		k = p[i];

		while (k > i)
			k = p[k];

		if (k < i)
			continue;

		/* Now have k == i, i.e the least in its cycle */

		pk = p[k];

		if (pk == i)
			continue;

		/* shuffle the elements of the cycle */

		{
			unsigned int a;

			ATOMIC t[MULTIPLICITY];

			for (a = 0; a < MULTIPLICITY; a++)
				t[a] = data[i*stride*MULTIPLICITY + a];

			while (pk != i)
			{
				for (a = 0; a < MULTIPLICITY; a++)
				{
					ATOMIC r1 = data[pk*stride*MULTIPLICITY + a];
					data[k*stride*MULTIPLICITY + a] = r1;
				}
				k = pk;
				pk = p[k];
			};

			for (a = 0; a < MULTIPLICITY; a++)
				data[k*stride*MULTIPLICITY + a] = t[a];
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_permute, inverse) (const size_t * p, ATOMIC * data, const size_t stride, const size_t n)
{
	size_t i, k, pk;

	for (i = 0; i < n; i++)
	{
		k = p[i];

		while (k > i)
			k = p[k];

		if (k < i)
			continue;

		/* Now have k == i, i.e the least in its cycle */

		pk = p[k];

		if (pk == i)
			continue;

		/* shuffle the elements of the cycle in the inverse direction */

		{
			unsigned int a;

			ATOMIC t[MULTIPLICITY];

			for (a = 0; a < MULTIPLICITY; a++)
				t[a] = data[k*stride*MULTIPLICITY + a];

			while (pk != i)
			{
				for (a = 0; a < MULTIPLICITY; a++)
				{
					ATOMIC r1 = data[pk*stride*MULTIPLICITY + a];
					data[pk*stride*MULTIPLICITY + a] = t[a];
					t[a] = r1;
				}

				k = pk;
				pk = p[k];
			};

			for (a = 0; a < MULTIPLICITY; a++)
				data[pk*stride*MULTIPLICITY + a] = t[a];
		}
	}

	return GSL_SUCCESS;
}


int
TYPE(gsl_permute_vector) (const gsl_permutation * p, TYPE(gsl_vector) * v)
{
	if (v->size != p->size)
	{
		GSL_ERROR("vector and permutation must be the same length", GSL_EBADLEN);
	}

	TYPE(gsl_permute) (p->data, v->data, v->stride, v->size);

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_permute_vector, inverse) (const gsl_permutation * p, TYPE(gsl_vector) * v)
{
	if (v->size != p->size)
	{
		GSL_ERROR("vector and permutation must be the same length", GSL_EBADLEN);
	}

	FUNCTION(gsl_permute, inverse) (p->data, v->data, v->stride, v->size);

	return GSL_SUCCESS;
}

#include "templates_off.h"
#undef  BASE_GSL_COMPLEX

#define BASE_GSL_COMPLEX_FLOAT
#include "templates_on.h"
int
TYPE(gsl_permute) (const size_t * p, ATOMIC * data, const size_t stride, const size_t n)
{
	size_t i, k, pk;

	for (i = 0; i < n; i++)
	{
		k = p[i];

		while (k > i)
			k = p[k];

		if (k < i)
			continue;

		/* Now have k == i, i.e the least in its cycle */

		pk = p[k];

		if (pk == i)
			continue;

		/* shuffle the elements of the cycle */

		{
			unsigned int a;

			ATOMIC t[MULTIPLICITY];

			for (a = 0; a < MULTIPLICITY; a++)
				t[a] = data[i*stride*MULTIPLICITY + a];

			while (pk != i)
			{
				for (a = 0; a < MULTIPLICITY; a++)
				{
					ATOMIC r1 = data[pk*stride*MULTIPLICITY + a];
					data[k*stride*MULTIPLICITY + a] = r1;
				}
				k = pk;
				pk = p[k];
			};

			for (a = 0; a < MULTIPLICITY; a++)
				data[k*stride*MULTIPLICITY + a] = t[a];
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_permute, inverse) (const size_t * p, ATOMIC * data, const size_t stride, const size_t n)
{
	size_t i, k, pk;

	for (i = 0; i < n; i++)
	{
		k = p[i];

		while (k > i)
			k = p[k];

		if (k < i)
			continue;

		/* Now have k == i, i.e the least in its cycle */

		pk = p[k];

		if (pk == i)
			continue;

		/* shuffle the elements of the cycle in the inverse direction */

		{
			unsigned int a;

			ATOMIC t[MULTIPLICITY];

			for (a = 0; a < MULTIPLICITY; a++)
				t[a] = data[k*stride*MULTIPLICITY + a];

			while (pk != i)
			{
				for (a = 0; a < MULTIPLICITY; a++)
				{
					ATOMIC r1 = data[pk*stride*MULTIPLICITY + a];
					data[pk*stride*MULTIPLICITY + a] = t[a];
					t[a] = r1;
				}

				k = pk;
				pk = p[k];
			};

			for (a = 0; a < MULTIPLICITY; a++)
				data[pk*stride*MULTIPLICITY + a] = t[a];
		}
	}

	return GSL_SUCCESS;
}


int
TYPE(gsl_permute_vector) (const gsl_permutation * p, TYPE(gsl_vector) * v)
{
	if (v->size != p->size)
	{
		GSL_ERROR("vector and permutation must be the same length", GSL_EBADLEN);
	}

	TYPE(gsl_permute) (p->data, v->data, v->stride, v->size);

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_permute_vector, inverse) (const gsl_permutation * p, TYPE(gsl_vector) * v)
{
	if (v->size != p->size)
	{
		GSL_ERROR("vector and permutation must be the same length", GSL_EBADLEN);
	}

	FUNCTION(gsl_permute, inverse) (p->data, v->data, v->stride, v->size);

	return GSL_SUCCESS;
}

#include "templates_off.h"
#undef  BASE_GSL_COMPLEX_FLOAT

#define BASE_LONG_DOUBLE
#include "templates_on.h"
int
TYPE(gsl_permute) (const size_t * p, ATOMIC * data, const size_t stride, const size_t n)
{
	size_t i, k, pk;

	for (i = 0; i < n; i++)
	{
		k = p[i];

		while (k > i)
			k = p[k];

		if (k < i)
			continue;

		/* Now have k == i, i.e the least in its cycle */

		pk = p[k];

		if (pk == i)
			continue;

		/* shuffle the elements of the cycle */

		{
			unsigned int a;

			ATOMIC t[MULTIPLICITY];

			for (a = 0; a < MULTIPLICITY; a++)
				t[a] = data[i*stride*MULTIPLICITY + a];

			while (pk != i)
			{
				for (a = 0; a < MULTIPLICITY; a++)
				{
					ATOMIC r1 = data[pk*stride*MULTIPLICITY + a];
					data[k*stride*MULTIPLICITY + a] = r1;
				}
				k = pk;
				pk = p[k];
			};

			for (a = 0; a < MULTIPLICITY; a++)
				data[k*stride*MULTIPLICITY + a] = t[a];
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_permute, inverse) (const size_t * p, ATOMIC * data, const size_t stride, const size_t n)
{
	size_t i, k, pk;

	for (i = 0; i < n; i++)
	{
		k = p[i];

		while (k > i)
			k = p[k];

		if (k < i)
			continue;

		/* Now have k == i, i.e the least in its cycle */

		pk = p[k];

		if (pk == i)
			continue;

		/* shuffle the elements of the cycle in the inverse direction */

		{
			unsigned int a;

			ATOMIC t[MULTIPLICITY];

			for (a = 0; a < MULTIPLICITY; a++)
				t[a] = data[k*stride*MULTIPLICITY + a];

			while (pk != i)
			{
				for (a = 0; a < MULTIPLICITY; a++)
				{
					ATOMIC r1 = data[pk*stride*MULTIPLICITY + a];
					data[pk*stride*MULTIPLICITY + a] = t[a];
					t[a] = r1;
				}

				k = pk;
				pk = p[k];
			};

			for (a = 0; a < MULTIPLICITY; a++)
				data[pk*stride*MULTIPLICITY + a] = t[a];
		}
	}

	return GSL_SUCCESS;
}


int
TYPE(gsl_permute_vector) (const gsl_permutation * p, TYPE(gsl_vector) * v)
{
	if (v->size != p->size)
	{
		GSL_ERROR("vector and permutation must be the same length", GSL_EBADLEN);
	}

	TYPE(gsl_permute) (p->data, v->data, v->stride, v->size);

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_permute_vector, inverse) (const gsl_permutation * p, TYPE(gsl_vector) * v)
{
	if (v->size != p->size)
	{
		GSL_ERROR("vector and permutation must be the same length", GSL_EBADLEN);
	}

	FUNCTION(gsl_permute, inverse) (p->data, v->data, v->stride, v->size);

	return GSL_SUCCESS;
}

#include "templates_off.h"
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
#include "templates_on.h"
int
TYPE(gsl_permute) (const size_t * p, ATOMIC * data, const size_t stride, const size_t n)
{
	size_t i, k, pk;

	for (i = 0; i < n; i++)
	{
		k = p[i];

		while (k > i)
			k = p[k];

		if (k < i)
			continue;

		/* Now have k == i, i.e the least in its cycle */

		pk = p[k];

		if (pk == i)
			continue;

		/* shuffle the elements of the cycle */

		{
			unsigned int a;

			ATOMIC t[MULTIPLICITY];

			for (a = 0; a < MULTIPLICITY; a++)
				t[a] = data[i*stride*MULTIPLICITY + a];

			while (pk != i)
			{
				for (a = 0; a < MULTIPLICITY; a++)
				{
					ATOMIC r1 = data[pk*stride*MULTIPLICITY + a];
					data[k*stride*MULTIPLICITY + a] = r1;
				}
				k = pk;
				pk = p[k];
			};

			for (a = 0; a < MULTIPLICITY; a++)
				data[k*stride*MULTIPLICITY + a] = t[a];
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_permute, inverse) (const size_t * p, ATOMIC * data, const size_t stride, const size_t n)
{
	size_t i, k, pk;

	for (i = 0; i < n; i++)
	{
		k = p[i];

		while (k > i)
			k = p[k];

		if (k < i)
			continue;

		/* Now have k == i, i.e the least in its cycle */

		pk = p[k];

		if (pk == i)
			continue;

		/* shuffle the elements of the cycle in the inverse direction */

		{
			unsigned int a;

			ATOMIC t[MULTIPLICITY];

			for (a = 0; a < MULTIPLICITY; a++)
				t[a] = data[k*stride*MULTIPLICITY + a];

			while (pk != i)
			{
				for (a = 0; a < MULTIPLICITY; a++)
				{
					ATOMIC r1 = data[pk*stride*MULTIPLICITY + a];
					data[pk*stride*MULTIPLICITY + a] = t[a];
					t[a] = r1;
				}

				k = pk;
				pk = p[k];
			};

			for (a = 0; a < MULTIPLICITY; a++)
				data[pk*stride*MULTIPLICITY + a] = t[a];
		}
	}

	return GSL_SUCCESS;
}


int
TYPE(gsl_permute_vector) (const gsl_permutation * p, TYPE(gsl_vector) * v)
{
	if (v->size != p->size)
	{
		GSL_ERROR("vector and permutation must be the same length", GSL_EBADLEN);
	}

	TYPE(gsl_permute) (p->data, v->data, v->stride, v->size);

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_permute_vector, inverse) (const gsl_permutation * p, TYPE(gsl_vector) * v)
{
	if (v->size != p->size)
	{
		GSL_ERROR("vector and permutation must be the same length", GSL_EBADLEN);
	}

	FUNCTION(gsl_permute, inverse) (p->data, v->data, v->stride, v->size);

	return GSL_SUCCESS;
}

#include "templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "templates_on.h"
int
TYPE(gsl_permute) (const size_t * p, ATOMIC * data, const size_t stride, const size_t n)
{
	size_t i, k, pk;

	for (i = 0; i < n; i++)
	{
		k = p[i];

		while (k > i)
			k = p[k];

		if (k < i)
			continue;

		/* Now have k == i, i.e the least in its cycle */

		pk = p[k];

		if (pk == i)
			continue;

		/* shuffle the elements of the cycle */

		{
			unsigned int a;

			ATOMIC t[MULTIPLICITY];

			for (a = 0; a < MULTIPLICITY; a++)
				t[a] = data[i*stride*MULTIPLICITY + a];

			while (pk != i)
			{
				for (a = 0; a < MULTIPLICITY; a++)
				{
					ATOMIC r1 = data[pk*stride*MULTIPLICITY + a];
					data[k*stride*MULTIPLICITY + a] = r1;
				}
				k = pk;
				pk = p[k];
			};

			for (a = 0; a < MULTIPLICITY; a++)
				data[k*stride*MULTIPLICITY + a] = t[a];
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_permute, inverse) (const size_t * p, ATOMIC * data, const size_t stride, const size_t n)
{
	size_t i, k, pk;

	for (i = 0; i < n; i++)
	{
		k = p[i];

		while (k > i)
			k = p[k];

		if (k < i)
			continue;

		/* Now have k == i, i.e the least in its cycle */

		pk = p[k];

		if (pk == i)
			continue;

		/* shuffle the elements of the cycle in the inverse direction */

		{
			unsigned int a;

			ATOMIC t[MULTIPLICITY];

			for (a = 0; a < MULTIPLICITY; a++)
				t[a] = data[k*stride*MULTIPLICITY + a];

			while (pk != i)
			{
				for (a = 0; a < MULTIPLICITY; a++)
				{
					ATOMIC r1 = data[pk*stride*MULTIPLICITY + a];
					data[pk*stride*MULTIPLICITY + a] = t[a];
					t[a] = r1;
				}

				k = pk;
				pk = p[k];
			};

			for (a = 0; a < MULTIPLICITY; a++)
				data[pk*stride*MULTIPLICITY + a] = t[a];
		}
	}

	return GSL_SUCCESS;
}


int
TYPE(gsl_permute_vector) (const gsl_permutation * p, TYPE(gsl_vector) * v)
{
	if (v->size != p->size)
	{
		GSL_ERROR("vector and permutation must be the same length", GSL_EBADLEN);
	}

	TYPE(gsl_permute) (p->data, v->data, v->stride, v->size);

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_permute_vector, inverse) (const gsl_permutation * p, TYPE(gsl_vector) * v)
{
	if (v->size != p->size)
	{
		GSL_ERROR("vector and permutation must be the same length", GSL_EBADLEN);
	}

	FUNCTION(gsl_permute, inverse) (p->data, v->data, v->stride, v->size);

	return GSL_SUCCESS;
}

#include "templates_off.h"
#undef  BASE_FLOAT

#define BASE_ULONG
#include "templates_on.h"
int
TYPE(gsl_permute) (const size_t * p, ATOMIC * data, const size_t stride, const size_t n)
{
	size_t i, k, pk;

	for (i = 0; i < n; i++)
	{
		k = p[i];

		while (k > i)
			k = p[k];

		if (k < i)
			continue;

		/* Now have k == i, i.e the least in its cycle */

		pk = p[k];

		if (pk == i)
			continue;

		/* shuffle the elements of the cycle */

		{
			unsigned int a;

			ATOMIC t[MULTIPLICITY];

			for (a = 0; a < MULTIPLICITY; a++)
				t[a] = data[i*stride*MULTIPLICITY + a];

			while (pk != i)
			{
				for (a = 0; a < MULTIPLICITY; a++)
				{
					ATOMIC r1 = data[pk*stride*MULTIPLICITY + a];
					data[k*stride*MULTIPLICITY + a] = r1;
				}
				k = pk;
				pk = p[k];
			};

			for (a = 0; a < MULTIPLICITY; a++)
				data[k*stride*MULTIPLICITY + a] = t[a];
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_permute, inverse) (const size_t * p, ATOMIC * data, const size_t stride, const size_t n)
{
	size_t i, k, pk;

	for (i = 0; i < n; i++)
	{
		k = p[i];

		while (k > i)
			k = p[k];

		if (k < i)
			continue;

		/* Now have k == i, i.e the least in its cycle */

		pk = p[k];

		if (pk == i)
			continue;

		/* shuffle the elements of the cycle in the inverse direction */

		{
			unsigned int a;

			ATOMIC t[MULTIPLICITY];

			for (a = 0; a < MULTIPLICITY; a++)
				t[a] = data[k*stride*MULTIPLICITY + a];

			while (pk != i)
			{
				for (a = 0; a < MULTIPLICITY; a++)
				{
					ATOMIC r1 = data[pk*stride*MULTIPLICITY + a];
					data[pk*stride*MULTIPLICITY + a] = t[a];
					t[a] = r1;
				}

				k = pk;
				pk = p[k];
			};

			for (a = 0; a < MULTIPLICITY; a++)
				data[pk*stride*MULTIPLICITY + a] = t[a];
		}
	}

	return GSL_SUCCESS;
}


int
TYPE(gsl_permute_vector) (const gsl_permutation * p, TYPE(gsl_vector) * v)
{
	if (v->size != p->size)
	{
		GSL_ERROR("vector and permutation must be the same length", GSL_EBADLEN);
	}

	TYPE(gsl_permute) (p->data, v->data, v->stride, v->size);

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_permute_vector, inverse) (const gsl_permutation * p, TYPE(gsl_vector) * v)
{
	if (v->size != p->size)
	{
		GSL_ERROR("vector and permutation must be the same length", GSL_EBADLEN);
	}

	FUNCTION(gsl_permute, inverse) (p->data, v->data, v->stride, v->size);

	return GSL_SUCCESS;
}

#include "templates_off.h"
#undef  BASE_ULONG

#define BASE_LONG
#include "templates_on.h"
int
TYPE(gsl_permute) (const size_t * p, ATOMIC * data, const size_t stride, const size_t n)
{
	size_t i, k, pk;

	for (i = 0; i < n; i++)
	{
		k = p[i];

		while (k > i)
			k = p[k];

		if (k < i)
			continue;

		/* Now have k == i, i.e the least in its cycle */

		pk = p[k];

		if (pk == i)
			continue;

		/* shuffle the elements of the cycle */

		{
			unsigned int a;

			ATOMIC t[MULTIPLICITY];

			for (a = 0; a < MULTIPLICITY; a++)
				t[a] = data[i*stride*MULTIPLICITY + a];

			while (pk != i)
			{
				for (a = 0; a < MULTIPLICITY; a++)
				{
					ATOMIC r1 = data[pk*stride*MULTIPLICITY + a];
					data[k*stride*MULTIPLICITY + a] = r1;
				}
				k = pk;
				pk = p[k];
			};

			for (a = 0; a < MULTIPLICITY; a++)
				data[k*stride*MULTIPLICITY + a] = t[a];
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_permute, inverse) (const size_t * p, ATOMIC * data, const size_t stride, const size_t n)
{
	size_t i, k, pk;

	for (i = 0; i < n; i++)
	{
		k = p[i];

		while (k > i)
			k = p[k];

		if (k < i)
			continue;

		/* Now have k == i, i.e the least in its cycle */

		pk = p[k];

		if (pk == i)
			continue;

		/* shuffle the elements of the cycle in the inverse direction */

		{
			unsigned int a;

			ATOMIC t[MULTIPLICITY];

			for (a = 0; a < MULTIPLICITY; a++)
				t[a] = data[k*stride*MULTIPLICITY + a];

			while (pk != i)
			{
				for (a = 0; a < MULTIPLICITY; a++)
				{
					ATOMIC r1 = data[pk*stride*MULTIPLICITY + a];
					data[pk*stride*MULTIPLICITY + a] = t[a];
					t[a] = r1;
				}

				k = pk;
				pk = p[k];
			};

			for (a = 0; a < MULTIPLICITY; a++)
				data[pk*stride*MULTIPLICITY + a] = t[a];
		}
	}

	return GSL_SUCCESS;
}


int
TYPE(gsl_permute_vector) (const gsl_permutation * p, TYPE(gsl_vector) * v)
{
	if (v->size != p->size)
	{
		GSL_ERROR("vector and permutation must be the same length", GSL_EBADLEN);
	}

	TYPE(gsl_permute) (p->data, v->data, v->stride, v->size);

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_permute_vector, inverse) (const gsl_permutation * p, TYPE(gsl_vector) * v)
{
	if (v->size != p->size)
	{
		GSL_ERROR("vector and permutation must be the same length", GSL_EBADLEN);
	}

	FUNCTION(gsl_permute, inverse) (p->data, v->data, v->stride, v->size);

	return GSL_SUCCESS;
}

#include "templates_off.h"
#undef  BASE_LONG

#define BASE_UINT
#include "templates_on.h"
int
TYPE(gsl_permute) (const size_t * p, ATOMIC * data, const size_t stride, const size_t n)
{
	size_t i, k, pk;

	for (i = 0; i < n; i++)
	{
		k = p[i];

		while (k > i)
			k = p[k];

		if (k < i)
			continue;

		/* Now have k == i, i.e the least in its cycle */

		pk = p[k];

		if (pk == i)
			continue;

		/* shuffle the elements of the cycle */

		{
			unsigned int a;

			ATOMIC t[MULTIPLICITY];

			for (a = 0; a < MULTIPLICITY; a++)
				t[a] = data[i*stride*MULTIPLICITY + a];

			while (pk != i)
			{
				for (a = 0; a < MULTIPLICITY; a++)
				{
					ATOMIC r1 = data[pk*stride*MULTIPLICITY + a];
					data[k*stride*MULTIPLICITY + a] = r1;
				}
				k = pk;
				pk = p[k];
			};

			for (a = 0; a < MULTIPLICITY; a++)
				data[k*stride*MULTIPLICITY + a] = t[a];
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_permute, inverse) (const size_t * p, ATOMIC * data, const size_t stride, const size_t n)
{
	size_t i, k, pk;

	for (i = 0; i < n; i++)
	{
		k = p[i];

		while (k > i)
			k = p[k];

		if (k < i)
			continue;

		/* Now have k == i, i.e the least in its cycle */

		pk = p[k];

		if (pk == i)
			continue;

		/* shuffle the elements of the cycle in the inverse direction */

		{
			unsigned int a;

			ATOMIC t[MULTIPLICITY];

			for (a = 0; a < MULTIPLICITY; a++)
				t[a] = data[k*stride*MULTIPLICITY + a];

			while (pk != i)
			{
				for (a = 0; a < MULTIPLICITY; a++)
				{
					ATOMIC r1 = data[pk*stride*MULTIPLICITY + a];
					data[pk*stride*MULTIPLICITY + a] = t[a];
					t[a] = r1;
				}

				k = pk;
				pk = p[k];
			};

			for (a = 0; a < MULTIPLICITY; a++)
				data[pk*stride*MULTIPLICITY + a] = t[a];
		}
	}

	return GSL_SUCCESS;
}


int
TYPE(gsl_permute_vector) (const gsl_permutation * p, TYPE(gsl_vector) * v)
{
	if (v->size != p->size)
	{
		GSL_ERROR("vector and permutation must be the same length", GSL_EBADLEN);
	}

	TYPE(gsl_permute) (p->data, v->data, v->stride, v->size);

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_permute_vector, inverse) (const gsl_permutation * p, TYPE(gsl_vector) * v)
{
	if (v->size != p->size)
	{
		GSL_ERROR("vector and permutation must be the same length", GSL_EBADLEN);
	}

	FUNCTION(gsl_permute, inverse) (p->data, v->data, v->stride, v->size);

	return GSL_SUCCESS;
}

#include "templates_off.h"
#undef  BASE_UINT

#define BASE_INT
#include "templates_on.h"
int
TYPE(gsl_permute) (const size_t * p, ATOMIC * data, const size_t stride, const size_t n)
{
	size_t i, k, pk;

	for (i = 0; i < n; i++)
	{
		k = p[i];

		while (k > i)
			k = p[k];

		if (k < i)
			continue;

		/* Now have k == i, i.e the least in its cycle */

		pk = p[k];

		if (pk == i)
			continue;

		/* shuffle the elements of the cycle */

		{
			unsigned int a;

			ATOMIC t[MULTIPLICITY];

			for (a = 0; a < MULTIPLICITY; a++)
				t[a] = data[i*stride*MULTIPLICITY + a];

			while (pk != i)
			{
				for (a = 0; a < MULTIPLICITY; a++)
				{
					ATOMIC r1 = data[pk*stride*MULTIPLICITY + a];
					data[k*stride*MULTIPLICITY + a] = r1;
				}
				k = pk;
				pk = p[k];
			};

			for (a = 0; a < MULTIPLICITY; a++)
				data[k*stride*MULTIPLICITY + a] = t[a];
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_permute, inverse) (const size_t * p, ATOMIC * data, const size_t stride, const size_t n)
{
	size_t i, k, pk;

	for (i = 0; i < n; i++)
	{
		k = p[i];

		while (k > i)
			k = p[k];

		if (k < i)
			continue;

		/* Now have k == i, i.e the least in its cycle */

		pk = p[k];

		if (pk == i)
			continue;

		/* shuffle the elements of the cycle in the inverse direction */

		{
			unsigned int a;

			ATOMIC t[MULTIPLICITY];

			for (a = 0; a < MULTIPLICITY; a++)
				t[a] = data[k*stride*MULTIPLICITY + a];

			while (pk != i)
			{
				for (a = 0; a < MULTIPLICITY; a++)
				{
					ATOMIC r1 = data[pk*stride*MULTIPLICITY + a];
					data[pk*stride*MULTIPLICITY + a] = t[a];
					t[a] = r1;
				}

				k = pk;
				pk = p[k];
			};

			for (a = 0; a < MULTIPLICITY; a++)
				data[pk*stride*MULTIPLICITY + a] = t[a];
		}
	}

	return GSL_SUCCESS;
}


int
TYPE(gsl_permute_vector) (const gsl_permutation * p, TYPE(gsl_vector) * v)
{
	if (v->size != p->size)
	{
		GSL_ERROR("vector and permutation must be the same length", GSL_EBADLEN);
	}

	TYPE(gsl_permute) (p->data, v->data, v->stride, v->size);

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_permute_vector, inverse) (const gsl_permutation * p, TYPE(gsl_vector) * v)
{
	if (v->size != p->size)
	{
		GSL_ERROR("vector and permutation must be the same length", GSL_EBADLEN);
	}

	FUNCTION(gsl_permute, inverse) (p->data, v->data, v->stride, v->size);

	return GSL_SUCCESS;
}

#include "templates_off.h"
#undef  BASE_INT

#define BASE_USHORT
#include "templates_on.h"
int
TYPE(gsl_permute) (const size_t * p, ATOMIC * data, const size_t stride, const size_t n)
{
	size_t i, k, pk;

	for (i = 0; i < n; i++)
	{
		k = p[i];

		while (k > i)
			k = p[k];

		if (k < i)
			continue;

		/* Now have k == i, i.e the least in its cycle */

		pk = p[k];

		if (pk == i)
			continue;

		/* shuffle the elements of the cycle */

		{
			unsigned int a;

			ATOMIC t[MULTIPLICITY];

			for (a = 0; a < MULTIPLICITY; a++)
				t[a] = data[i*stride*MULTIPLICITY + a];

			while (pk != i)
			{
				for (a = 0; a < MULTIPLICITY; a++)
				{
					ATOMIC r1 = data[pk*stride*MULTIPLICITY + a];
					data[k*stride*MULTIPLICITY + a] = r1;
				}
				k = pk;
				pk = p[k];
			};

			for (a = 0; a < MULTIPLICITY; a++)
				data[k*stride*MULTIPLICITY + a] = t[a];
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_permute, inverse) (const size_t * p, ATOMIC * data, const size_t stride, const size_t n)
{
	size_t i, k, pk;

	for (i = 0; i < n; i++)
	{
		k = p[i];

		while (k > i)
			k = p[k];

		if (k < i)
			continue;

		/* Now have k == i, i.e the least in its cycle */

		pk = p[k];

		if (pk == i)
			continue;

		/* shuffle the elements of the cycle in the inverse direction */

		{
			unsigned int a;

			ATOMIC t[MULTIPLICITY];

			for (a = 0; a < MULTIPLICITY; a++)
				t[a] = data[k*stride*MULTIPLICITY + a];

			while (pk != i)
			{
				for (a = 0; a < MULTIPLICITY; a++)
				{
					ATOMIC r1 = data[pk*stride*MULTIPLICITY + a];
					data[pk*stride*MULTIPLICITY + a] = t[a];
					t[a] = r1;
				}

				k = pk;
				pk = p[k];
			};

			for (a = 0; a < MULTIPLICITY; a++)
				data[pk*stride*MULTIPLICITY + a] = t[a];
		}
	}

	return GSL_SUCCESS;
}


int
TYPE(gsl_permute_vector) (const gsl_permutation * p, TYPE(gsl_vector) * v)
{
	if (v->size != p->size)
	{
		GSL_ERROR("vector and permutation must be the same length", GSL_EBADLEN);
	}

	TYPE(gsl_permute) (p->data, v->data, v->stride, v->size);

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_permute_vector, inverse) (const gsl_permutation * p, TYPE(gsl_vector) * v)
{
	if (v->size != p->size)
	{
		GSL_ERROR("vector and permutation must be the same length", GSL_EBADLEN);
	}

	FUNCTION(gsl_permute, inverse) (p->data, v->data, v->stride, v->size);

	return GSL_SUCCESS;
}

#include "templates_off.h"
#undef  BASE_USHORT

#define BASE_SHORT
#include "templates_on.h"
int
TYPE(gsl_permute) (const size_t * p, ATOMIC * data, const size_t stride, const size_t n)
{
	size_t i, k, pk;

	for (i = 0; i < n; i++)
	{
		k = p[i];

		while (k > i)
			k = p[k];

		if (k < i)
			continue;

		/* Now have k == i, i.e the least in its cycle */

		pk = p[k];

		if (pk == i)
			continue;

		/* shuffle the elements of the cycle */

		{
			unsigned int a;

			ATOMIC t[MULTIPLICITY];

			for (a = 0; a < MULTIPLICITY; a++)
				t[a] = data[i*stride*MULTIPLICITY + a];

			while (pk != i)
			{
				for (a = 0; a < MULTIPLICITY; a++)
				{
					ATOMIC r1 = data[pk*stride*MULTIPLICITY + a];
					data[k*stride*MULTIPLICITY + a] = r1;
				}
				k = pk;
				pk = p[k];
			};

			for (a = 0; a < MULTIPLICITY; a++)
				data[k*stride*MULTIPLICITY + a] = t[a];
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_permute, inverse) (const size_t * p, ATOMIC * data, const size_t stride, const size_t n)
{
	size_t i, k, pk;

	for (i = 0; i < n; i++)
	{
		k = p[i];

		while (k > i)
			k = p[k];

		if (k < i)
			continue;

		/* Now have k == i, i.e the least in its cycle */

		pk = p[k];

		if (pk == i)
			continue;

		/* shuffle the elements of the cycle in the inverse direction */

		{
			unsigned int a;

			ATOMIC t[MULTIPLICITY];

			for (a = 0; a < MULTIPLICITY; a++)
				t[a] = data[k*stride*MULTIPLICITY + a];

			while (pk != i)
			{
				for (a = 0; a < MULTIPLICITY; a++)
				{
					ATOMIC r1 = data[pk*stride*MULTIPLICITY + a];
					data[pk*stride*MULTIPLICITY + a] = t[a];
					t[a] = r1;
				}

				k = pk;
				pk = p[k];
			};

			for (a = 0; a < MULTIPLICITY; a++)
				data[pk*stride*MULTIPLICITY + a] = t[a];
		}
	}

	return GSL_SUCCESS;
}


int
TYPE(gsl_permute_vector) (const gsl_permutation * p, TYPE(gsl_vector) * v)
{
	if (v->size != p->size)
	{
		GSL_ERROR("vector and permutation must be the same length", GSL_EBADLEN);
	}

	TYPE(gsl_permute) (p->data, v->data, v->stride, v->size);

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_permute_vector, inverse) (const gsl_permutation * p, TYPE(gsl_vector) * v)
{
	if (v->size != p->size)
	{
		GSL_ERROR("vector and permutation must be the same length", GSL_EBADLEN);
	}

	FUNCTION(gsl_permute, inverse) (p->data, v->data, v->stride, v->size);

	return GSL_SUCCESS;
}

#include "templates_off.h"
#undef  BASE_SHORT

#define BASE_UCHAR
#include "templates_on.h"
int
TYPE(gsl_permute) (const size_t * p, ATOMIC * data, const size_t stride, const size_t n)
{
	size_t i, k, pk;

	for (i = 0; i < n; i++)
	{
		k = p[i];

		while (k > i)
			k = p[k];

		if (k < i)
			continue;

		/* Now have k == i, i.e the least in its cycle */

		pk = p[k];

		if (pk == i)
			continue;

		/* shuffle the elements of the cycle */

		{
			unsigned int a;

			ATOMIC t[MULTIPLICITY];

			for (a = 0; a < MULTIPLICITY; a++)
				t[a] = data[i*stride*MULTIPLICITY + a];

			while (pk != i)
			{
				for (a = 0; a < MULTIPLICITY; a++)
				{
					ATOMIC r1 = data[pk*stride*MULTIPLICITY + a];
					data[k*stride*MULTIPLICITY + a] = r1;
				}
				k = pk;
				pk = p[k];
			};

			for (a = 0; a < MULTIPLICITY; a++)
				data[k*stride*MULTIPLICITY + a] = t[a];
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_permute, inverse) (const size_t * p, ATOMIC * data, const size_t stride, const size_t n)
{
	size_t i, k, pk;

	for (i = 0; i < n; i++)
	{
		k = p[i];

		while (k > i)
			k = p[k];

		if (k < i)
			continue;

		/* Now have k == i, i.e the least in its cycle */

		pk = p[k];

		if (pk == i)
			continue;

		/* shuffle the elements of the cycle in the inverse direction */

		{
			unsigned int a;

			ATOMIC t[MULTIPLICITY];

			for (a = 0; a < MULTIPLICITY; a++)
				t[a] = data[k*stride*MULTIPLICITY + a];

			while (pk != i)
			{
				for (a = 0; a < MULTIPLICITY; a++)
				{
					ATOMIC r1 = data[pk*stride*MULTIPLICITY + a];
					data[pk*stride*MULTIPLICITY + a] = t[a];
					t[a] = r1;
				}

				k = pk;
				pk = p[k];
			};

			for (a = 0; a < MULTIPLICITY; a++)
				data[pk*stride*MULTIPLICITY + a] = t[a];
		}
	}

	return GSL_SUCCESS;
}


int
TYPE(gsl_permute_vector) (const gsl_permutation * p, TYPE(gsl_vector) * v)
{
	if (v->size != p->size)
	{
		GSL_ERROR("vector and permutation must be the same length", GSL_EBADLEN);
	}

	TYPE(gsl_permute) (p->data, v->data, v->stride, v->size);

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_permute_vector, inverse) (const gsl_permutation * p, TYPE(gsl_vector) * v)
{
	if (v->size != p->size)
	{
		GSL_ERROR("vector and permutation must be the same length", GSL_EBADLEN);
	}

	FUNCTION(gsl_permute, inverse) (p->data, v->data, v->stride, v->size);

	return GSL_SUCCESS;
}

#include "templates_off.h"
#undef  BASE_UCHAR

#define BASE_CHAR
#include "templates_on.h"
int
TYPE(gsl_permute) (const size_t * p, ATOMIC * data, const size_t stride, const size_t n)
{
	size_t i, k, pk;

	for (i = 0; i < n; i++)
	{
		k = p[i];

		while (k > i)
			k = p[k];

		if (k < i)
			continue;

		/* Now have k == i, i.e the least in its cycle */

		pk = p[k];

		if (pk == i)
			continue;

		/* shuffle the elements of the cycle */

		{
			unsigned int a;

			ATOMIC t[MULTIPLICITY];

			for (a = 0; a < MULTIPLICITY; a++)
				t[a] = data[i*stride*MULTIPLICITY + a];

			while (pk != i)
			{
				for (a = 0; a < MULTIPLICITY; a++)
				{
					ATOMIC r1 = data[pk*stride*MULTIPLICITY + a];
					data[k*stride*MULTIPLICITY + a] = r1;
				}
				k = pk;
				pk = p[k];
			};

			for (a = 0; a < MULTIPLICITY; a++)
				data[k*stride*MULTIPLICITY + a] = t[a];
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_permute, inverse) (const size_t * p, ATOMIC * data, const size_t stride, const size_t n)
{
	size_t i, k, pk;

	for (i = 0; i < n; i++)
	{
		k = p[i];

		while (k > i)
			k = p[k];

		if (k < i)
			continue;

		/* Now have k == i, i.e the least in its cycle */

		pk = p[k];

		if (pk == i)
			continue;

		/* shuffle the elements of the cycle in the inverse direction */

		{
			unsigned int a;

			ATOMIC t[MULTIPLICITY];

			for (a = 0; a < MULTIPLICITY; a++)
				t[a] = data[k*stride*MULTIPLICITY + a];

			while (pk != i)
			{
				for (a = 0; a < MULTIPLICITY; a++)
				{
					ATOMIC r1 = data[pk*stride*MULTIPLICITY + a];
					data[pk*stride*MULTIPLICITY + a] = t[a];
					t[a] = r1;
				}

				k = pk;
				pk = p[k];
			};

			for (a = 0; a < MULTIPLICITY; a++)
				data[pk*stride*MULTIPLICITY + a] = t[a];
		}
	}

	return GSL_SUCCESS;
}


int
TYPE(gsl_permute_vector) (const gsl_permutation * p, TYPE(gsl_vector) * v)
{
	if (v->size != p->size)
	{
		GSL_ERROR("vector and permutation must be the same length", GSL_EBADLEN);
	}

	TYPE(gsl_permute) (p->data, v->data, v->stride, v->size);

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_permute_vector, inverse) (const gsl_permutation * p, TYPE(gsl_vector) * v)
{
	if (v->size != p->size)
	{
		GSL_ERROR("vector and permutation must be the same length", GSL_EBADLEN);
	}

	FUNCTION(gsl_permute, inverse) (p->data, v->data, v->stride, v->size);

	return GSL_SUCCESS;
}

#include "templates_off.h"
#undef  BASE_CHAR
