#include "stdafx.h"
#include <config.h.in>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/matrix/gsl_matrix.h>

#define BASE_GSL_COMPLEX_LONG
#include "templates_on.h"
int
FUNCTION(gsl_matrix, add) (TYPE(gsl_matrix) * a,
const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				const size_t aij = 2 * (i * tda_a + j);
				const size_t bij = 2 * (i * tda_b + j);

				a->data[aij] += b->data[bij];
				a->data[aij + 1] += b->data[bij + 1];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, sub) (TYPE(gsl_matrix) * a,
const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				const size_t aij = 2 * (i * tda_a + j);
				const size_t bij = 2 * (i * tda_b + j);

				a->data[aij] -= b->data[bij];
				a->data[aij + 1] -= b->data[bij + 1];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, mul_elements) (TYPE(gsl_matrix) * a,
const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				const size_t aij = 2 * (i * tda_a + j);
				const size_t bij = 2 * (i * tda_b + j);

				ATOMIC ar = a->data[aij];
				ATOMIC ai = a->data[aij + 1];

				ATOMIC br = b->data[bij];
				ATOMIC bi = b->data[bij + 1];

				a->data[aij] = ar * br - ai * bi;
				a->data[aij + 1] = ar * bi + ai * br;
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, div_elements) (TYPE(gsl_matrix) * a,
const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				const size_t aij = 2 * (i * tda_a + j);
				const size_t bij = 2 * (i * tda_b + j);

				ATOMIC ar = a->data[aij];
				ATOMIC ai = a->data[aij + 1];

				ATOMIC br = b->data[bij];
				ATOMIC bi = b->data[bij + 1];

				ATOMIC s = 1.0 / hypot(br, bi);

				ATOMIC sbr = s * br;
				ATOMIC sbi = s * bi;

				a->data[aij] = (ar * sbr + ai * sbi) * s;
				a->data[aij + 1] = (ai * sbr - ar * sbi) * s;
			}
		}

		return GSL_SUCCESS;
	}
}

int FUNCTION(gsl_matrix, scale) (TYPE(gsl_matrix) * a, const BASE x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;

	size_t i, j;

	ATOMIC xr = GSL_REAL(x);
	ATOMIC xi = GSL_IMAG(x);

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			const size_t aij = 2 * (i * tda + j);

			ATOMIC ar = a->data[aij];
			ATOMIC ai = a->data[aij + 1];

			a->data[aij] = ar * xr - ai * xi;
			a->data[aij + 1] = ar * xi + ai * xr;
		}
	}

	return GSL_SUCCESS;
}

int FUNCTION(gsl_matrix, add_constant) (TYPE(gsl_matrix) * a, const BASE x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;

	size_t i, j;

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			a->data[2 * (i * tda + j)] += GSL_REAL(x);
			a->data[2 * (i * tda + j) + 1] += GSL_IMAG(x);
		}
	}

	return GSL_SUCCESS;
}


int FUNCTION(gsl_matrix, add_diagonal) (TYPE(gsl_matrix) * a, const BASE x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;
	const size_t loop_lim = (M < N ? M : N);
	size_t i;
	for (i = 0; i < loop_lim; i++)
	{
		a->data[2 * (i * tda + i)] += GSL_REAL(x);
		a->data[2 * (i * tda + i) + 1] += GSL_IMAG(x);
	}

	return GSL_SUCCESS;
}

#include "templates_off.h"
#undef  BASE_GSL_COMPLEX_LONG

#define BASE_GSL_COMPLEX
#include "templates_on.h"
int
FUNCTION(gsl_matrix, add) (TYPE(gsl_matrix) * a,
const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				const size_t aij = 2 * (i * tda_a + j);
				const size_t bij = 2 * (i * tda_b + j);

				a->data[aij] += b->data[bij];
				a->data[aij + 1] += b->data[bij + 1];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, sub) (TYPE(gsl_matrix) * a,
const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				const size_t aij = 2 * (i * tda_a + j);
				const size_t bij = 2 * (i * tda_b + j);

				a->data[aij] -= b->data[bij];
				a->data[aij + 1] -= b->data[bij + 1];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, mul_elements) (TYPE(gsl_matrix) * a,
const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				const size_t aij = 2 * (i * tda_a + j);
				const size_t bij = 2 * (i * tda_b + j);

				ATOMIC ar = a->data[aij];
				ATOMIC ai = a->data[aij + 1];

				ATOMIC br = b->data[bij];
				ATOMIC bi = b->data[bij + 1];

				a->data[aij] = ar * br - ai * bi;
				a->data[aij + 1] = ar * bi + ai * br;
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, div_elements) (TYPE(gsl_matrix) * a,
const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				const size_t aij = 2 * (i * tda_a + j);
				const size_t bij = 2 * (i * tda_b + j);

				ATOMIC ar = a->data[aij];
				ATOMIC ai = a->data[aij + 1];

				ATOMIC br = b->data[bij];
				ATOMIC bi = b->data[bij + 1];

				ATOMIC s = 1.0 / hypot(br, bi);

				ATOMIC sbr = s * br;
				ATOMIC sbi = s * bi;

				a->data[aij] = (ar * sbr + ai * sbi) * s;
				a->data[aij + 1] = (ai * sbr - ar * sbi) * s;
			}
		}

		return GSL_SUCCESS;
	}
}

int FUNCTION(gsl_matrix, scale) (TYPE(gsl_matrix) * a, const BASE x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;

	size_t i, j;

	ATOMIC xr = GSL_REAL(x);
	ATOMIC xi = GSL_IMAG(x);

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			const size_t aij = 2 * (i * tda + j);

			ATOMIC ar = a->data[aij];
			ATOMIC ai = a->data[aij + 1];

			a->data[aij] = ar * xr - ai * xi;
			a->data[aij + 1] = ar * xi + ai * xr;
		}
	}

	return GSL_SUCCESS;
}

int FUNCTION(gsl_matrix, add_constant) (TYPE(gsl_matrix) * a, const BASE x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;

	size_t i, j;

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			a->data[2 * (i * tda + j)] += GSL_REAL(x);
			a->data[2 * (i * tda + j) + 1] += GSL_IMAG(x);
		}
	}

	return GSL_SUCCESS;
}


int FUNCTION(gsl_matrix, add_diagonal) (TYPE(gsl_matrix) * a, const BASE x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;
	const size_t loop_lim = (M < N ? M : N);
	size_t i;
	for (i = 0; i < loop_lim; i++)
	{
		a->data[2 * (i * tda + i)] += GSL_REAL(x);
		a->data[2 * (i * tda + i) + 1] += GSL_IMAG(x);
	}

	return GSL_SUCCESS;
}

#include "templates_off.h"
#undef  BASE_GSL_COMPLEX

#define BASE_GSL_COMPLEX_FLOAT
#include "templates_on.h"
int
FUNCTION(gsl_matrix, add) (TYPE(gsl_matrix) * a,
const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				const size_t aij = 2 * (i * tda_a + j);
				const size_t bij = 2 * (i * tda_b + j);

				a->data[aij] += b->data[bij];
				a->data[aij + 1] += b->data[bij + 1];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, sub) (TYPE(gsl_matrix) * a,
const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				const size_t aij = 2 * (i * tda_a + j);
				const size_t bij = 2 * (i * tda_b + j);

				a->data[aij] -= b->data[bij];
				a->data[aij + 1] -= b->data[bij + 1];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, mul_elements) (TYPE(gsl_matrix) * a,
const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				const size_t aij = 2 * (i * tda_a + j);
				const size_t bij = 2 * (i * tda_b + j);

				ATOMIC ar = a->data[aij];
				ATOMIC ai = a->data[aij + 1];

				ATOMIC br = b->data[bij];
				ATOMIC bi = b->data[bij + 1];

				a->data[aij] = ar * br - ai * bi;
				a->data[aij + 1] = ar * bi + ai * br;
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, div_elements) (TYPE(gsl_matrix) * a,
const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				const size_t aij = 2 * (i * tda_a + j);
				const size_t bij = 2 * (i * tda_b + j);

				ATOMIC ar = a->data[aij];
				ATOMIC ai = a->data[aij + 1];

				ATOMIC br = b->data[bij];
				ATOMIC bi = b->data[bij + 1];

				ATOMIC s = 1.0 / hypot(br, bi);

				ATOMIC sbr = s * br;
				ATOMIC sbi = s * bi;

				a->data[aij] = (ar * sbr + ai * sbi) * s;
				a->data[aij + 1] = (ai * sbr - ar * sbi) * s;
			}
		}

		return GSL_SUCCESS;
	}
}

int FUNCTION(gsl_matrix, scale) (TYPE(gsl_matrix) * a, const BASE x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;

	size_t i, j;

	ATOMIC xr = GSL_REAL(x);
	ATOMIC xi = GSL_IMAG(x);

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			const size_t aij = 2 * (i * tda + j);

			ATOMIC ar = a->data[aij];
			ATOMIC ai = a->data[aij + 1];

			a->data[aij] = ar * xr - ai * xi;
			a->data[aij + 1] = ar * xi + ai * xr;
		}
	}

	return GSL_SUCCESS;
}

int FUNCTION(gsl_matrix, add_constant) (TYPE(gsl_matrix) * a, const BASE x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;

	size_t i, j;

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			a->data[2 * (i * tda + j)] += GSL_REAL(x);
			a->data[2 * (i * tda + j) + 1] += GSL_IMAG(x);
		}
	}

	return GSL_SUCCESS;
}


int FUNCTION(gsl_matrix, add_diagonal) (TYPE(gsl_matrix) * a, const BASE x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;
	const size_t loop_lim = (M < N ? M : N);
	size_t i;
	for (i = 0; i < loop_lim; i++)
	{
		a->data[2 * (i * tda + i)] += GSL_REAL(x);
		a->data[2 * (i * tda + i) + 1] += GSL_IMAG(x);
	}

	return GSL_SUCCESS;
}

#include "templates_off.h"
#undef  BASE_GSL_COMPLEX_FLOAT

#define BASE_LONG_DOUBLE
#include "templates_on.h"
int
FUNCTION(gsl_matrix, add) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] += b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, sub) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] -= b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, mul_elements) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] *= b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, div_elements) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] /= b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, scale) (TYPE(gsl_matrix) * a, const double x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;

	size_t i, j;

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			a->data[i * tda + j] *= x;
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, add_constant) (TYPE(gsl_matrix) * a, const double x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;

	size_t i, j;

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			a->data[i * tda + j] += x;
		}
	}

	return GSL_SUCCESS;
}


int
FUNCTION(gsl_matrix, add_diagonal) (TYPE(gsl_matrix) * a, const double x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;
	const size_t loop_lim = (M < N ? M : N);
	size_t i;
	for (i = 0; i < loop_lim; i++)
	{
		a->data[i * tda + i] += x;
	}

	return GSL_SUCCESS;
}

#include "templates_off.h"
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
#include "templates_on.h"
int
FUNCTION(gsl_matrix, add) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] += b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, sub) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] -= b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, mul_elements) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] *= b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, div_elements) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] /= b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, scale) (TYPE(gsl_matrix) * a, const double x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;

	size_t i, j;

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			a->data[i * tda + j] *= x;
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, add_constant) (TYPE(gsl_matrix) * a, const double x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;

	size_t i, j;

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			a->data[i * tda + j] += x;
		}
	}

	return GSL_SUCCESS;
}


int
FUNCTION(gsl_matrix, add_diagonal) (TYPE(gsl_matrix) * a, const double x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;
	const size_t loop_lim = (M < N ? M : N);
	size_t i;
	for (i = 0; i < loop_lim; i++)
	{
		a->data[i * tda + i] += x;
	}

	return GSL_SUCCESS;
}

#include "templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "templates_on.h"
int
FUNCTION(gsl_matrix, add) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] += b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, sub) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] -= b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, mul_elements) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] *= b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, div_elements) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] /= b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, scale) (TYPE(gsl_matrix) * a, const double x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;

	size_t i, j;

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			a->data[i * tda + j] *= x;
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, add_constant) (TYPE(gsl_matrix) * a, const double x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;

	size_t i, j;

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			a->data[i * tda + j] += x;
		}
	}

	return GSL_SUCCESS;
}


int
FUNCTION(gsl_matrix, add_diagonal) (TYPE(gsl_matrix) * a, const double x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;
	const size_t loop_lim = (M < N ? M : N);
	size_t i;
	for (i = 0; i < loop_lim; i++)
	{
		a->data[i * tda + i] += x;
	}

	return GSL_SUCCESS;
}

#include "templates_off.h"
#undef  BASE_FLOAT

#define BASE_ULONG
#include "templates_on.h"
int
FUNCTION(gsl_matrix, add) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] += b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, sub) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] -= b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, mul_elements) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] *= b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, div_elements) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] /= b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, scale) (TYPE(gsl_matrix) * a, const double x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;

	size_t i, j;

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			a->data[i * tda + j] *= x;
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, add_constant) (TYPE(gsl_matrix) * a, const double x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;

	size_t i, j;

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			a->data[i * tda + j] += x;
		}
	}

	return GSL_SUCCESS;
}


int
FUNCTION(gsl_matrix, add_diagonal) (TYPE(gsl_matrix) * a, const double x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;
	const size_t loop_lim = (M < N ? M : N);
	size_t i;
	for (i = 0; i < loop_lim; i++)
	{
		a->data[i * tda + i] += x;
	}

	return GSL_SUCCESS;
}

#include "templates_off.h"
#undef  BASE_ULONG

#define BASE_LONG
#include "templates_on.h"
int
FUNCTION(gsl_matrix, add) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] += b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, sub) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] -= b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, mul_elements) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] *= b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, div_elements) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] /= b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, scale) (TYPE(gsl_matrix) * a, const double x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;

	size_t i, j;

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			a->data[i * tda + j] *= x;
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, add_constant) (TYPE(gsl_matrix) * a, const double x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;

	size_t i, j;

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			a->data[i * tda + j] += x;
		}
	}

	return GSL_SUCCESS;
}


int
FUNCTION(gsl_matrix, add_diagonal) (TYPE(gsl_matrix) * a, const double x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;
	const size_t loop_lim = (M < N ? M : N);
	size_t i;
	for (i = 0; i < loop_lim; i++)
	{
		a->data[i * tda + i] += x;
	}

	return GSL_SUCCESS;
}

#include "templates_off.h"
#undef  BASE_LONG

#define BASE_UINT
#include "templates_on.h"
int
FUNCTION(gsl_matrix, add) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] += b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, sub) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] -= b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, mul_elements) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] *= b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, div_elements) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] /= b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, scale) (TYPE(gsl_matrix) * a, const double x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;

	size_t i, j;

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			a->data[i * tda + j] *= x;
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, add_constant) (TYPE(gsl_matrix) * a, const double x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;

	size_t i, j;

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			a->data[i * tda + j] += x;
		}
	}

	return GSL_SUCCESS;
}


int
FUNCTION(gsl_matrix, add_diagonal) (TYPE(gsl_matrix) * a, const double x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;
	const size_t loop_lim = (M < N ? M : N);
	size_t i;
	for (i = 0; i < loop_lim; i++)
	{
		a->data[i * tda + i] += x;
	}

	return GSL_SUCCESS;
}

#include "templates_off.h"
#undef  BASE_UINT

#define BASE_INT
#include "templates_on.h"
int
FUNCTION(gsl_matrix, add) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] += b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, sub) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] -= b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, mul_elements) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] *= b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, div_elements) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] /= b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, scale) (TYPE(gsl_matrix) * a, const double x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;

	size_t i, j;

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			a->data[i * tda + j] *= x;
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, add_constant) (TYPE(gsl_matrix) * a, const double x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;

	size_t i, j;

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			a->data[i * tda + j] += x;
		}
	}

	return GSL_SUCCESS;
}


int
FUNCTION(gsl_matrix, add_diagonal) (TYPE(gsl_matrix) * a, const double x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;
	const size_t loop_lim = (M < N ? M : N);
	size_t i;
	for (i = 0; i < loop_lim; i++)
	{
		a->data[i * tda + i] += x;
	}

	return GSL_SUCCESS;
}

#include "templates_off.h"
#undef  BASE_INT

#define BASE_USHORT
#include "templates_on.h"
int
FUNCTION(gsl_matrix, add) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] += b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, sub) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] -= b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, mul_elements) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] *= b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, div_elements) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] /= b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, scale) (TYPE(gsl_matrix) * a, const double x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;

	size_t i, j;

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			a->data[i * tda + j] *= x;
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, add_constant) (TYPE(gsl_matrix) * a, const double x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;

	size_t i, j;

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			a->data[i * tda + j] += x;
		}
	}

	return GSL_SUCCESS;
}


int
FUNCTION(gsl_matrix, add_diagonal) (TYPE(gsl_matrix) * a, const double x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;
	const size_t loop_lim = (M < N ? M : N);
	size_t i;
	for (i = 0; i < loop_lim; i++)
	{
		a->data[i * tda + i] += x;
	}

	return GSL_SUCCESS;
}

#include "templates_off.h"
#undef  BASE_USHORT

#define BASE_SHORT
#include "templates_on.h"
int
FUNCTION(gsl_matrix, add) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] += b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, sub) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] -= b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, mul_elements) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] *= b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, div_elements) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] /= b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, scale) (TYPE(gsl_matrix) * a, const double x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;

	size_t i, j;

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			a->data[i * tda + j] *= x;
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, add_constant) (TYPE(gsl_matrix) * a, const double x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;

	size_t i, j;

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			a->data[i * tda + j] += x;
		}
	}

	return GSL_SUCCESS;
}


int
FUNCTION(gsl_matrix, add_diagonal) (TYPE(gsl_matrix) * a, const double x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;
	const size_t loop_lim = (M < N ? M : N);
	size_t i;
	for (i = 0; i < loop_lim; i++)
	{
		a->data[i * tda + i] += x;
	}

	return GSL_SUCCESS;
}

#include "templates_off.h"
#undef  BASE_SHORT

#define BASE_UCHAR
#include "templates_on.h"
int
FUNCTION(gsl_matrix, add) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] += b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, sub) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] -= b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, mul_elements) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] *= b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, div_elements) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] /= b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, scale) (TYPE(gsl_matrix) * a, const double x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;

	size_t i, j;

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			a->data[i * tda + j] *= x;
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, add_constant) (TYPE(gsl_matrix) * a, const double x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;

	size_t i, j;

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			a->data[i * tda + j] += x;
		}
	}

	return GSL_SUCCESS;
}


int
FUNCTION(gsl_matrix, add_diagonal) (TYPE(gsl_matrix) * a, const double x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;
	const size_t loop_lim = (M < N ? M : N);
	size_t i;
	for (i = 0; i < loop_lim; i++)
	{
		a->data[i * tda + i] += x;
	}

	return GSL_SUCCESS;
}

#include "templates_off.h"
#undef  BASE_UCHAR

#define BASE_CHAR
#include "templates_on.h"
int
FUNCTION(gsl_matrix, add) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] += b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, sub) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] -= b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, mul_elements) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] *= b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, div_elements) (TYPE(gsl_matrix) * a, const TYPE(gsl_matrix) * b)
{
	const size_t M = a->size1;
	const size_t N = a->size2;

	if (b->size1 != M || b->size2 != N)
	{
		GSL_ERROR("matrices must have same dimensions", GSL_EBADLEN);
	}
	else
	{
		const size_t tda_a = a->tda;
		const size_t tda_b = b->tda;

		size_t i, j;

		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				a->data[i * tda_a + j] /= b->data[i * tda_b + j];
			}
		}

		return GSL_SUCCESS;
	}
}

int
FUNCTION(gsl_matrix, scale) (TYPE(gsl_matrix) * a, const double x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;

	size_t i, j;

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			a->data[i * tda + j] *= x;
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_matrix, add_constant) (TYPE(gsl_matrix) * a, const double x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;

	size_t i, j;

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			a->data[i * tda + j] += x;
		}
	}

	return GSL_SUCCESS;
}


int
FUNCTION(gsl_matrix, add_diagonal) (TYPE(gsl_matrix) * a, const double x)
{
	const size_t M = a->size1;
	const size_t N = a->size2;
	const size_t tda = a->tda;
	const size_t loop_lim = (M < N ? M : N);
	size_t i;
	for (i = 0; i < loop_lim; i++)
	{
		a->data[i * tda + i] += x;
	}

	return GSL_SUCCESS;
}

#include "templates_off.h"
#undef  BASE_CHAR


