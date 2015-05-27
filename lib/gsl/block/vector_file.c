#include "stdafx.h"
#include <config.h.in>
#include <stdio.h>
#include <gsl/err/gsl_errno.h>
#include <gsl/block/gsl_block.h>

#define BASE_GSL_COMPLEX_LONG
#include "templates_on.h"
int
FUNCTION (gsl_block, fread) (FILE * stream, TYPE(gsl_block) * b)
{
	size_t n = b->size ;

	ATOMIC * data = b->data ;
	
	size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

	if (items != n)
	{
		GSL_ERROR ("fread failed", GSL_EFAILED);
	}

	return 0;
}

int
FUNCTION (gsl_block, fwrite) (FILE * stream, const TYPE(gsl_block) * b)
{
	size_t n = b->size ;

	ATOMIC * data = b->data ;

	size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

	if (items != n)
	{
		GSL_ERROR ("fwrite failed", GSL_EFAILED);
	}

	return 0;
}

int
FUNCTION (gsl_block, raw_fread) (FILE * stream, ATOMIC * data, 
const size_t n, const size_t stride)
{
	if (stride == 1)
	{
		size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

		if (items != n)
		{
			GSL_ERROR ("fread failed", GSL_EFAILED);
		}
	}
	else
	{
		size_t i;

		for (i = 0; i < n; i++)
		{
			size_t item = fread (data + MULTIPLICITY * i * stride,
				MULTIPLICITY * sizeof (ATOMIC), 1, stream);
			if (item != 1)
			{
				GSL_ERROR ("fread failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION (gsl_block, raw_fwrite) (FILE * stream, const ATOMIC * data,
const size_t n, const size_t stride)
{

	if (stride == 1)
	{
		size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

		if (items != n)
		{
			GSL_ERROR ("fwrite failed", GSL_EFAILED);
		}
	}
	else
	{
		size_t i;

		for (i = 0; i < n; i++)
		{
			size_t item = fwrite (data + MULTIPLICITY * i * stride,
				MULTIPLICITY * sizeof (ATOMIC),
				1, stream);
			if (item != 1)
			{
				GSL_ERROR ("fwrite failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}
#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)

int
FUNCTION(gsl_block, fprintf) (FILE * stream, const TYPE(gsl_block) * b, const char *format)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		int status;

		for (k = 0; k < MULTIPLICITY; k++)
		{
			if (k > 0)
			{
				status = putc(' ', stream);

				if (status == EOF)
				{
					GSL_ERROR("putc failed", GSL_EFAILED);
				}
			}
			status = fprintf(stream,
				format,
				data[MULTIPLICITY * i + k]);
			if (status < 0)
			{
				GSL_ERROR("fprintf failed", GSL_EFAILED);
			}
		}

		status = putc('\n', stream);

		if (status == EOF)
		{
			GSL_ERROR("putc failed", GSL_EFAILED);
		}
	}

	return 0;
}

int
FUNCTION(gsl_block, fscanf) (FILE * stream, TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		for (k = 0; k < MULTIPLICITY; k++)
		{
			ATOMIC_IO tmp;

			int status = fscanf(stream, IN_FORMAT, &tmp);

			data[MULTIPLICITY * i + k] = tmp;


			if (status != 1)
			{
				GSL_ERROR("fscanf failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}


int
FUNCTION(gsl_block, raw_fprintf) (FILE * stream,
const ATOMIC * data,
const size_t n,
const size_t stride,
const char *format)
{
	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		int status;

		for (k = 0; k < MULTIPLICITY; k++)
		{
			if (k > 0)
			{
				status = putc(' ', stream);

				if (status == EOF)
				{
					GSL_ERROR("putc failed", GSL_EFAILED);
				}
			}
			status = fprintf(stream,
				format,
				data[MULTIPLICITY * i * stride + k]);
			if (status < 0)
			{
				GSL_ERROR("fprintf failed", GSL_EFAILED);
			}
		}

		status = putc('\n', stream);

		if (status == EOF)
		{
			GSL_ERROR("putc failed", GSL_EFAILED);
		}
	}

	return 0;
}

int
FUNCTION(gsl_block, raw_fscanf) (FILE * stream,
ATOMIC * data,
const size_t n,
const size_t stride)
{
	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		for (k = 0; k < MULTIPLICITY; k++)
		{
			ATOMIC_IO tmp;

			int status = fscanf(stream, IN_FORMAT, &tmp);

			data[MULTIPLICITY * i * stride + k] = tmp;

			if (status != 1)
				GSL_ERROR("fscanf failed", GSL_EFAILED);
		}
	}

	return GSL_SUCCESS;
}

#endif

#include "templates_off.h"
#undef  BASE_GSL_COMPLEX_LONG

#define BASE_GSL_COMPLEX
#include "templates_on.h"
int
FUNCTION(gsl_block, fread) (FILE * stream, TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t items = fread(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

	if (items != n)
	{
		GSL_ERROR("fread failed", GSL_EFAILED);
	}

	return 0;
}

int
FUNCTION(gsl_block, fwrite) (FILE * stream, const TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t items = fwrite(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

	if (items != n)
	{
		GSL_ERROR("fwrite failed", GSL_EFAILED);
	}

	return 0;
}

int
FUNCTION(gsl_block, raw_fread) (FILE * stream, ATOMIC * data,
const size_t n, const size_t stride)
{
	if (stride == 1)
	{
		size_t items = fread(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

		if (items != n)
		{
			GSL_ERROR("fread failed", GSL_EFAILED);
		}
	}
	else
	{
		size_t i;

		for (i = 0; i < n; i++)
		{
			size_t item = fread(data + MULTIPLICITY * i * stride,
				MULTIPLICITY * sizeof(ATOMIC), 1, stream);
			if (item != 1)
			{
				GSL_ERROR("fread failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_block, raw_fwrite) (FILE * stream, const ATOMIC * data,
const size_t n, const size_t stride)
{

	if (stride == 1)
	{
		size_t items = fwrite(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

		if (items != n)
		{
			GSL_ERROR("fwrite failed", GSL_EFAILED);
		}
	}
	else
	{
		size_t i;

		for (i = 0; i < n; i++)
		{
			size_t item = fwrite(data + MULTIPLICITY * i * stride,
				MULTIPLICITY * sizeof(ATOMIC),
				1, stream);
			if (item != 1)
			{
				GSL_ERROR("fwrite failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)

int
FUNCTION(gsl_block, fprintf) (FILE * stream, const TYPE(gsl_block) * b, const char *format)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		int status;

		for (k = 0; k < MULTIPLICITY; k++)
		{
			if (k > 0)
			{
				status = putc(' ', stream);

				if (status == EOF)
				{
					GSL_ERROR("putc failed", GSL_EFAILED);
				}
			}
			status = fprintf(stream,
				format,
				data[MULTIPLICITY * i + k]);
			if (status < 0)
			{
				GSL_ERROR("fprintf failed", GSL_EFAILED);
			}
		}

		status = putc('\n', stream);

		if (status == EOF)
		{
			GSL_ERROR("putc failed", GSL_EFAILED);
		}
	}

	return 0;
}

int
FUNCTION(gsl_block, fscanf) (FILE * stream, TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		for (k = 0; k < MULTIPLICITY; k++)
		{
			ATOMIC_IO tmp;

			int status = fscanf(stream, IN_FORMAT, &tmp);

			data[MULTIPLICITY * i + k] = tmp;


			if (status != 1)
			{
				GSL_ERROR("fscanf failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}


int
FUNCTION(gsl_block, raw_fprintf) (FILE * stream,
const ATOMIC * data,
const size_t n,
const size_t stride,
const char *format)
{
	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		int status;

		for (k = 0; k < MULTIPLICITY; k++)
		{
			if (k > 0)
			{
				status = putc(' ', stream);

				if (status == EOF)
				{
					GSL_ERROR("putc failed", GSL_EFAILED);
				}
			}
			status = fprintf(stream,
				format,
				data[MULTIPLICITY * i * stride + k]);
			if (status < 0)
			{
				GSL_ERROR("fprintf failed", GSL_EFAILED);
			}
		}

		status = putc('\n', stream);

		if (status == EOF)
		{
			GSL_ERROR("putc failed", GSL_EFAILED);
		}
	}

	return 0;
}

int
FUNCTION(gsl_block, raw_fscanf) (FILE * stream,
ATOMIC * data,
const size_t n,
const size_t stride)
{
	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		for (k = 0; k < MULTIPLICITY; k++)
		{
			ATOMIC_IO tmp;

			int status = fscanf(stream, IN_FORMAT, &tmp);

			data[MULTIPLICITY * i * stride + k] = tmp;

			if (status != 1)
				GSL_ERROR("fscanf failed", GSL_EFAILED);
		}
	}

	return GSL_SUCCESS;
}

#endif

#include "templates_off.h"
#undef  BASE_GSL_COMPLEX

#define BASE_GSL_COMPLEX_FLOAT
#include "templates_on.h"
int
FUNCTION(gsl_block, fread) (FILE * stream, TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t items = fread(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

	if (items != n)
	{
		GSL_ERROR("fread failed", GSL_EFAILED);
	}

	return 0;
}

int
FUNCTION(gsl_block, fwrite) (FILE * stream, const TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t items = fwrite(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

	if (items != n)
	{
		GSL_ERROR("fwrite failed", GSL_EFAILED);
	}

	return 0;
}

int
FUNCTION(gsl_block, raw_fread) (FILE * stream, ATOMIC * data,
const size_t n, const size_t stride)
{
	if (stride == 1)
	{
		size_t items = fread(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

		if (items != n)
		{
			GSL_ERROR("fread failed", GSL_EFAILED);
		}
	}
	else
	{
		size_t i;

		for (i = 0; i < n; i++)
		{
			size_t item = fread(data + MULTIPLICITY * i * stride,
				MULTIPLICITY * sizeof(ATOMIC), 1, stream);
			if (item != 1)
			{
				GSL_ERROR("fread failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_block, raw_fwrite) (FILE * stream, const ATOMIC * data,
const size_t n, const size_t stride)
{

	if (stride == 1)
	{
		size_t items = fwrite(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

		if (items != n)
		{
			GSL_ERROR("fwrite failed", GSL_EFAILED);
		}
	}
	else
	{
		size_t i;

		for (i = 0; i < n; i++)
		{
			size_t item = fwrite(data + MULTIPLICITY * i * stride,
				MULTIPLICITY * sizeof(ATOMIC),
				1, stream);
			if (item != 1)
			{
				GSL_ERROR("fwrite failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)

int
FUNCTION(gsl_block, fprintf) (FILE * stream, const TYPE(gsl_block) * b, const char *format)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		int status;

		for (k = 0; k < MULTIPLICITY; k++)
		{
			if (k > 0)
			{
				status = putc(' ', stream);

				if (status == EOF)
				{
					GSL_ERROR("putc failed", GSL_EFAILED);
				}
			}
			status = fprintf(stream,
				format,
				data[MULTIPLICITY * i + k]);
			if (status < 0)
			{
				GSL_ERROR("fprintf failed", GSL_EFAILED);
			}
		}

		status = putc('\n', stream);

		if (status == EOF)
		{
			GSL_ERROR("putc failed", GSL_EFAILED);
		}
	}

	return 0;
}

int
FUNCTION(gsl_block, fscanf) (FILE * stream, TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		for (k = 0; k < MULTIPLICITY; k++)
		{
			ATOMIC_IO tmp;

			int status = fscanf(stream, IN_FORMAT, &tmp);

			data[MULTIPLICITY * i + k] = tmp;


			if (status != 1)
			{
				GSL_ERROR("fscanf failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}


int
FUNCTION(gsl_block, raw_fprintf) (FILE * stream,
const ATOMIC * data,
const size_t n,
const size_t stride,
const char *format)
{
	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		int status;

		for (k = 0; k < MULTIPLICITY; k++)
		{
			if (k > 0)
			{
				status = putc(' ', stream);

				if (status == EOF)
				{
					GSL_ERROR("putc failed", GSL_EFAILED);
				}
			}
			status = fprintf(stream,
				format,
				data[MULTIPLICITY * i * stride + k]);
			if (status < 0)
			{
				GSL_ERROR("fprintf failed", GSL_EFAILED);
			}
		}

		status = putc('\n', stream);

		if (status == EOF)
		{
			GSL_ERROR("putc failed", GSL_EFAILED);
		}
	}

	return 0;
}

int
FUNCTION(gsl_block, raw_fscanf) (FILE * stream,
ATOMIC * data,
const size_t n,
const size_t stride)
{
	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		for (k = 0; k < MULTIPLICITY; k++)
		{
			ATOMIC_IO tmp;

			int status = fscanf(stream, IN_FORMAT, &tmp);

			data[MULTIPLICITY * i * stride + k] = tmp;

			if (status != 1)
				GSL_ERROR("fscanf failed", GSL_EFAILED);
		}
	}

	return GSL_SUCCESS;
}

#endif

#include "templates_off.h"
#undef  BASE_GSL_COMPLEX_FLOAT

#define BASE_LONG_DOUBLE
#include "templates_on.h"
int
FUNCTION (gsl_block, fread) (FILE * stream, TYPE(gsl_block) * b)
{
	size_t n = b->size ;

	ATOMIC * data = b->data ;

	size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

	if (items != n)
	{
		GSL_ERROR ("fread failed", GSL_EFAILED);
	}

	return 0;
}

int
FUNCTION (gsl_block, fwrite) (FILE * stream, const TYPE(gsl_block) * b)
{
	size_t n = b->size ;

	ATOMIC * data = b->data ;

	size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

	if (items != n)
	{
		GSL_ERROR ("fwrite failed", GSL_EFAILED);
	}

	return 0;
}

int
FUNCTION (gsl_block, raw_fread) (FILE * stream, ATOMIC * data, 
const size_t n, const size_t stride)
{
	if (stride == 1)
	{
		size_t items = fread (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

		if (items != n)
		{
			GSL_ERROR ("fread failed", GSL_EFAILED);
		}
	}
	else
	{
		size_t i;

		for (i = 0; i < n; i++)
		{
			size_t item = fread (data + MULTIPLICITY * i * stride,
				MULTIPLICITY * sizeof (ATOMIC), 1, stream);
			if (item != 1)
			{
				GSL_ERROR ("fread failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION (gsl_block, raw_fwrite) (FILE * stream, const ATOMIC * data,
const size_t n, const size_t stride)
{

	if (stride == 1)
	{
		size_t items = fwrite (data, MULTIPLICITY * sizeof (ATOMIC), n, stream);

		if (items != n)
		{
			GSL_ERROR ("fwrite failed", GSL_EFAILED);
		}
	}
	else
	{
		size_t i;

		for (i = 0; i < n; i++)
		{
			size_t item = fwrite (data + MULTIPLICITY * i * stride,
				MULTIPLICITY * sizeof (ATOMIC),
				1, stream);
			if (item != 1)
			{
				GSL_ERROR ("fwrite failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)

int
FUNCTION(gsl_block, fprintf) (FILE * stream, const TYPE(gsl_block) * b, const char *format)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		int status;

		for (k = 0; k < MULTIPLICITY; k++)
		{
			if (k > 0)
			{
				status = putc(' ', stream);

				if (status == EOF)
				{
					GSL_ERROR("putc failed", GSL_EFAILED);
				}
			}
			status = fprintf(stream,
				format,
				data[MULTIPLICITY * i + k]);
			if (status < 0)
			{
				GSL_ERROR("fprintf failed", GSL_EFAILED);
			}
		}

		status = putc('\n', stream);

		if (status == EOF)
		{
			GSL_ERROR("putc failed", GSL_EFAILED);
		}
	}

	return 0;
}

int
FUNCTION(gsl_block, fscanf) (FILE * stream, TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		for (k = 0; k < MULTIPLICITY; k++)
		{
			ATOMIC_IO tmp;

			int status = fscanf(stream, IN_FORMAT, &tmp);

			data[MULTIPLICITY * i + k] = tmp;


			if (status != 1)
			{
				GSL_ERROR("fscanf failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}


int
FUNCTION(gsl_block, raw_fprintf) (FILE * stream,
const ATOMIC * data,
const size_t n,
const size_t stride,
const char *format)
{
	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		int status;

		for (k = 0; k < MULTIPLICITY; k++)
		{
			if (k > 0)
			{
				status = putc(' ', stream);

				if (status == EOF)
				{
					GSL_ERROR("putc failed", GSL_EFAILED);
				}
			}
			status = fprintf(stream,
				format,
				data[MULTIPLICITY * i * stride + k]);
			if (status < 0)
			{
				GSL_ERROR("fprintf failed", GSL_EFAILED);
			}
		}

		status = putc('\n', stream);

		if (status == EOF)
		{
			GSL_ERROR("putc failed", GSL_EFAILED);
		}
	}

	return 0;
}

int
FUNCTION(gsl_block, raw_fscanf) (FILE * stream,
ATOMIC * data,
const size_t n,
const size_t stride)
{
	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		for (k = 0; k < MULTIPLICITY; k++)
		{
			ATOMIC_IO tmp;

			int status = fscanf(stream, IN_FORMAT, &tmp);

			data[MULTIPLICITY * i * stride + k] = tmp;

			if (status != 1)
				GSL_ERROR("fscanf failed", GSL_EFAILED);
		}
	}

	return GSL_SUCCESS;
}

#endif

#include "templates_off.h"
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
#include "templates_on.h"
int
FUNCTION(gsl_block, fread) (FILE * stream, TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t items = fread(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

	if (items != n)
	{
		GSL_ERROR("fread failed", GSL_EFAILED);
	}

	return 0;
}

int
FUNCTION(gsl_block, fwrite) (FILE * stream, const TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t items = fwrite(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

	if (items != n)
	{
		GSL_ERROR("fwrite failed", GSL_EFAILED);
	}

	return 0;
}

int
FUNCTION(gsl_block, raw_fread) (FILE * stream, ATOMIC * data,
const size_t n, const size_t stride)
{
	if (stride == 1)
	{
		size_t items = fread(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

		if (items != n)
		{
			GSL_ERROR("fread failed", GSL_EFAILED);
		}
	}
	else
	{
		size_t i;

		for (i = 0; i < n; i++)
		{
			size_t item = fread(data + MULTIPLICITY * i * stride,
				MULTIPLICITY * sizeof(ATOMIC), 1, stream);
			if (item != 1)
			{
				GSL_ERROR("fread failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_block, raw_fwrite) (FILE * stream, const ATOMIC * data,
const size_t n, const size_t stride)
{

	if (stride == 1)
	{
		size_t items = fwrite(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

		if (items != n)
		{
			GSL_ERROR("fwrite failed", GSL_EFAILED);
		}
	}
	else
	{
		size_t i;

		for (i = 0; i < n; i++)
		{
			size_t item = fwrite(data + MULTIPLICITY * i * stride,
				MULTIPLICITY * sizeof(ATOMIC),
				1, stream);
			if (item != 1)
			{
				GSL_ERROR("fwrite failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)

int
FUNCTION(gsl_block, fprintf) (FILE * stream, const TYPE(gsl_block) * b, const char *format)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		int status;

		for (k = 0; k < MULTIPLICITY; k++)
		{
			if (k > 0)
			{
				status = putc(' ', stream);

				if (status == EOF)
				{
					GSL_ERROR("putc failed", GSL_EFAILED);
				}
			}
			status = fprintf(stream,
				format,
				data[MULTIPLICITY * i + k]);
			if (status < 0)
			{
				GSL_ERROR("fprintf failed", GSL_EFAILED);
			}
		}

		status = putc('\n', stream);

		if (status == EOF)
		{
			GSL_ERROR("putc failed", GSL_EFAILED);
		}
	}

	return 0;
}

int
FUNCTION(gsl_block, fscanf) (FILE * stream, TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		for (k = 0; k < MULTIPLICITY; k++)
		{
			ATOMIC_IO tmp;

			int status = fscanf(stream, IN_FORMAT, &tmp);

			data[MULTIPLICITY * i + k] = tmp;


			if (status != 1)
			{
				GSL_ERROR("fscanf failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}


int
FUNCTION(gsl_block, raw_fprintf) (FILE * stream,
const ATOMIC * data,
const size_t n,
const size_t stride,
const char *format)
{
	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		int status;

		for (k = 0; k < MULTIPLICITY; k++)
		{
			if (k > 0)
			{
				status = putc(' ', stream);

				if (status == EOF)
				{
					GSL_ERROR("putc failed", GSL_EFAILED);
				}
			}
			status = fprintf(stream,
				format,
				data[MULTIPLICITY * i * stride + k]);
			if (status < 0)
			{
				GSL_ERROR("fprintf failed", GSL_EFAILED);
			}
		}

		status = putc('\n', stream);

		if (status == EOF)
		{
			GSL_ERROR("putc failed", GSL_EFAILED);
		}
	}

	return 0;
}

int
FUNCTION(gsl_block, raw_fscanf) (FILE * stream,
ATOMIC * data,
const size_t n,
const size_t stride)
{
	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		for (k = 0; k < MULTIPLICITY; k++)
		{
			ATOMIC_IO tmp;

			int status = fscanf(stream, IN_FORMAT, &tmp);

			data[MULTIPLICITY * i * stride + k] = tmp;

			if (status != 1)
				GSL_ERROR("fscanf failed", GSL_EFAILED);
		}
	}

	return GSL_SUCCESS;
}

#endif

#include "templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "templates_on.h"
int
FUNCTION(gsl_block, fread) (FILE * stream, TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t items = fread(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

	if (items != n)
	{
		GSL_ERROR("fread failed", GSL_EFAILED);
	}

	return 0;
}

int
FUNCTION(gsl_block, fwrite) (FILE * stream, const TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t items = fwrite(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

	if (items != n)
	{
		GSL_ERROR("fwrite failed", GSL_EFAILED);
	}

	return 0;
}

int
FUNCTION(gsl_block, raw_fread) (FILE * stream, ATOMIC * data,
const size_t n, const size_t stride)
{
	if (stride == 1)
	{
		size_t items = fread(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

		if (items != n)
		{
			GSL_ERROR("fread failed", GSL_EFAILED);
		}
	}
	else
	{
		size_t i;

		for (i = 0; i < n; i++)
		{
			size_t item = fread(data + MULTIPLICITY * i * stride,
				MULTIPLICITY * sizeof(ATOMIC), 1, stream);
			if (item != 1)
			{
				GSL_ERROR("fread failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_block, raw_fwrite) (FILE * stream, const ATOMIC * data,
const size_t n, const size_t stride)
{

	if (stride == 1)
	{
		size_t items = fwrite(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

		if (items != n)
		{
			GSL_ERROR("fwrite failed", GSL_EFAILED);
		}
	}
	else
	{
		size_t i;

		for (i = 0; i < n; i++)
		{
			size_t item = fwrite(data + MULTIPLICITY * i * stride,
				MULTIPLICITY * sizeof(ATOMIC),
				1, stream);
			if (item != 1)
			{
				GSL_ERROR("fwrite failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)

int
FUNCTION(gsl_block, fprintf) (FILE * stream, const TYPE(gsl_block) * b, const char *format)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		int status;

		for (k = 0; k < MULTIPLICITY; k++)
		{
			if (k > 0)
			{
				status = putc(' ', stream);

				if (status == EOF)
				{
					GSL_ERROR("putc failed", GSL_EFAILED);
				}
			}
			status = fprintf(stream,
				format,
				data[MULTIPLICITY * i + k]);
			if (status < 0)
			{
				GSL_ERROR("fprintf failed", GSL_EFAILED);
			}
		}

		status = putc('\n', stream);

		if (status == EOF)
		{
			GSL_ERROR("putc failed", GSL_EFAILED);
		}
	}

	return 0;
}

int
FUNCTION(gsl_block, fscanf) (FILE * stream, TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		for (k = 0; k < MULTIPLICITY; k++)
		{
			ATOMIC_IO tmp;

			int status = fscanf(stream, IN_FORMAT, &tmp);

			data[MULTIPLICITY * i + k] = tmp;


			if (status != 1)
			{
				GSL_ERROR("fscanf failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}


int
FUNCTION(gsl_block, raw_fprintf) (FILE * stream,
const ATOMIC * data,
const size_t n,
const size_t stride,
const char *format)
{
	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		int status;

		for (k = 0; k < MULTIPLICITY; k++)
		{
			if (k > 0)
			{
				status = putc(' ', stream);

				if (status == EOF)
				{
					GSL_ERROR("putc failed", GSL_EFAILED);
				}
			}
			status = fprintf(stream,
				format,
				data[MULTIPLICITY * i * stride + k]);
			if (status < 0)
			{
				GSL_ERROR("fprintf failed", GSL_EFAILED);
			}
		}

		status = putc('\n', stream);

		if (status == EOF)
		{
			GSL_ERROR("putc failed", GSL_EFAILED);
		}
	}

	return 0;
}

int
FUNCTION(gsl_block, raw_fscanf) (FILE * stream,
ATOMIC * data,
const size_t n,
const size_t stride)
{
	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		for (k = 0; k < MULTIPLICITY; k++)
		{
			ATOMIC_IO tmp;

			int status = fscanf(stream, IN_FORMAT, &tmp);

			data[MULTIPLICITY * i * stride + k] = tmp;

			if (status != 1)
				GSL_ERROR("fscanf failed", GSL_EFAILED);
		}
	}

	return GSL_SUCCESS;
}

#endif

#include "templates_off.h"
#undef  BASE_FLOAT

#define BASE_ULONG
#include "templates_on.h"
int
FUNCTION(gsl_block, fread) (FILE * stream, TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t items = fread(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

	if (items != n)
	{
		GSL_ERROR("fread failed", GSL_EFAILED);
	}

	return 0;
}

int
FUNCTION(gsl_block, fwrite) (FILE * stream, const TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t items = fwrite(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

	if (items != n)
	{
		GSL_ERROR("fwrite failed", GSL_EFAILED);
	}

	return 0;
}

int
FUNCTION(gsl_block, raw_fread) (FILE * stream, ATOMIC * data,
const size_t n, const size_t stride)
{
	if (stride == 1)
	{
		size_t items = fread(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

		if (items != n)
		{
			GSL_ERROR("fread failed", GSL_EFAILED);
		}
	}
	else
	{
		size_t i;

		for (i = 0; i < n; i++)
		{
			size_t item = fread(data + MULTIPLICITY * i * stride,
				MULTIPLICITY * sizeof(ATOMIC), 1, stream);
			if (item != 1)
			{
				GSL_ERROR("fread failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_block, raw_fwrite) (FILE * stream, const ATOMIC * data,
const size_t n, const size_t stride)
{

	if (stride == 1)
	{
		size_t items = fwrite(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

		if (items != n)
		{
			GSL_ERROR("fwrite failed", GSL_EFAILED);
		}
	}
	else
	{
		size_t i;

		for (i = 0; i < n; i++)
		{
			size_t item = fwrite(data + MULTIPLICITY * i * stride,
				MULTIPLICITY * sizeof(ATOMIC),
				1, stream);
			if (item != 1)
			{
				GSL_ERROR("fwrite failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)

int
FUNCTION(gsl_block, fprintf) (FILE * stream, const TYPE(gsl_block) * b, const char *format)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		int status;

		for (k = 0; k < MULTIPLICITY; k++)
		{
			if (k > 0)
			{
				status = putc(' ', stream);

				if (status == EOF)
				{
					GSL_ERROR("putc failed", GSL_EFAILED);
				}
			}
			status = fprintf(stream,
				format,
				data[MULTIPLICITY * i + k]);
			if (status < 0)
			{
				GSL_ERROR("fprintf failed", GSL_EFAILED);
			}
		}

		status = putc('\n', stream);

		if (status == EOF)
		{
			GSL_ERROR("putc failed", GSL_EFAILED);
		}
	}

	return 0;
}

int
FUNCTION(gsl_block, fscanf) (FILE * stream, TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		for (k = 0; k < MULTIPLICITY; k++)
		{
			ATOMIC_IO tmp;

			int status = fscanf(stream, IN_FORMAT, &tmp);

			data[MULTIPLICITY * i + k] = tmp;


			if (status != 1)
			{
				GSL_ERROR("fscanf failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}


int
FUNCTION(gsl_block, raw_fprintf) (FILE * stream,
const ATOMIC * data,
const size_t n,
const size_t stride,
const char *format)
{
	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		int status;

		for (k = 0; k < MULTIPLICITY; k++)
		{
			if (k > 0)
			{
				status = putc(' ', stream);

				if (status == EOF)
				{
					GSL_ERROR("putc failed", GSL_EFAILED);
				}
			}
			status = fprintf(stream,
				format,
				data[MULTIPLICITY * i * stride + k]);
			if (status < 0)
			{
				GSL_ERROR("fprintf failed", GSL_EFAILED);
			}
		}

		status = putc('\n', stream);

		if (status == EOF)
		{
			GSL_ERROR("putc failed", GSL_EFAILED);
		}
	}

	return 0;
}

int
FUNCTION(gsl_block, raw_fscanf) (FILE * stream,
ATOMIC * data,
const size_t n,
const size_t stride)
{
	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		for (k = 0; k < MULTIPLICITY; k++)
		{
			ATOMIC_IO tmp;

			int status = fscanf(stream, IN_FORMAT, &tmp);

			data[MULTIPLICITY * i * stride + k] = tmp;

			if (status != 1)
				GSL_ERROR("fscanf failed", GSL_EFAILED);
		}
	}

	return GSL_SUCCESS;
}

#endif

#include "templates_off.h"
#undef  BASE_ULONG

#define BASE_LONG
#include "templates_on.h"
int
FUNCTION(gsl_block, fread) (FILE * stream, TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t items = fread(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

	if (items != n)
	{
		GSL_ERROR("fread failed", GSL_EFAILED);
	}

	return 0;
}

int
FUNCTION(gsl_block, fwrite) (FILE * stream, const TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t items = fwrite(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

	if (items != n)
	{
		GSL_ERROR("fwrite failed", GSL_EFAILED);
	}

	return 0;
}

int
FUNCTION(gsl_block, raw_fread) (FILE * stream, ATOMIC * data,
const size_t n, const size_t stride)
{
	if (stride == 1)
	{
		size_t items = fread(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

		if (items != n)
		{
			GSL_ERROR("fread failed", GSL_EFAILED);
		}
	}
	else
	{
		size_t i;

		for (i = 0; i < n; i++)
		{
			size_t item = fread(data + MULTIPLICITY * i * stride,
				MULTIPLICITY * sizeof(ATOMIC), 1, stream);
			if (item != 1)
			{
				GSL_ERROR("fread failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_block, raw_fwrite) (FILE * stream, const ATOMIC * data,
const size_t n, const size_t stride)
{

	if (stride == 1)
	{
		size_t items = fwrite(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

		if (items != n)
		{
			GSL_ERROR("fwrite failed", GSL_EFAILED);
		}
	}
	else
	{
		size_t i;

		for (i = 0; i < n; i++)
		{
			size_t item = fwrite(data + MULTIPLICITY * i * stride,
				MULTIPLICITY * sizeof(ATOMIC),
				1, stream);
			if (item != 1)
			{
				GSL_ERROR("fwrite failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)

int
FUNCTION(gsl_block, fprintf) (FILE * stream, const TYPE(gsl_block) * b, const char *format)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		int status;

		for (k = 0; k < MULTIPLICITY; k++)
		{
			if (k > 0)
			{
				status = putc(' ', stream);

				if (status == EOF)
				{
					GSL_ERROR("putc failed", GSL_EFAILED);
				}
			}
			status = fprintf(stream,
				format,
				data[MULTIPLICITY * i + k]);
			if (status < 0)
			{
				GSL_ERROR("fprintf failed", GSL_EFAILED);
			}
		}

		status = putc('\n', stream);

		if (status == EOF)
		{
			GSL_ERROR("putc failed", GSL_EFAILED);
		}
	}

	return 0;
}

int
FUNCTION(gsl_block, fscanf) (FILE * stream, TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		for (k = 0; k < MULTIPLICITY; k++)
		{
			ATOMIC_IO tmp;

			int status = fscanf(stream, IN_FORMAT, &tmp);

			data[MULTIPLICITY * i + k] = tmp;


			if (status != 1)
			{
				GSL_ERROR("fscanf failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}


int
FUNCTION(gsl_block, raw_fprintf) (FILE * stream,
const ATOMIC * data,
const size_t n,
const size_t stride,
const char *format)
{
	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		int status;

		for (k = 0; k < MULTIPLICITY; k++)
		{
			if (k > 0)
			{
				status = putc(' ', stream);

				if (status == EOF)
				{
					GSL_ERROR("putc failed", GSL_EFAILED);
				}
			}
			status = fprintf(stream,
				format,
				data[MULTIPLICITY * i * stride + k]);
			if (status < 0)
			{
				GSL_ERROR("fprintf failed", GSL_EFAILED);
			}
		}

		status = putc('\n', stream);

		if (status == EOF)
		{
			GSL_ERROR("putc failed", GSL_EFAILED);
		}
	}

	return 0;
}

int
FUNCTION(gsl_block, raw_fscanf) (FILE * stream,
ATOMIC * data,
const size_t n,
const size_t stride)
{
	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		for (k = 0; k < MULTIPLICITY; k++)
		{
			ATOMIC_IO tmp;

			int status = fscanf(stream, IN_FORMAT, &tmp);

			data[MULTIPLICITY * i * stride + k] = tmp;

			if (status != 1)
				GSL_ERROR("fscanf failed", GSL_EFAILED);
		}
	}

	return GSL_SUCCESS;
}

#endif

#include "templates_off.h"
#undef  BASE_LONG

#define BASE_UINT
#include "templates_on.h"
int
FUNCTION(gsl_block, fread) (FILE * stream, TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t items = fread(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

	if (items != n)
	{
		GSL_ERROR("fread failed", GSL_EFAILED);
	}

	return 0;
}

int
FUNCTION(gsl_block, fwrite) (FILE * stream, const TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t items = fwrite(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

	if (items != n)
	{
		GSL_ERROR("fwrite failed", GSL_EFAILED);
	}

	return 0;
}

int
FUNCTION(gsl_block, raw_fread) (FILE * stream, ATOMIC * data,
const size_t n, const size_t stride)
{
	if (stride == 1)
	{
		size_t items = fread(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

		if (items != n)
		{
			GSL_ERROR("fread failed", GSL_EFAILED);
		}
	}
	else
	{
		size_t i;

		for (i = 0; i < n; i++)
		{
			size_t item = fread(data + MULTIPLICITY * i * stride,
				MULTIPLICITY * sizeof(ATOMIC), 1, stream);
			if (item != 1)
			{
				GSL_ERROR("fread failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_block, raw_fwrite) (FILE * stream, const ATOMIC * data,
const size_t n, const size_t stride)
{

	if (stride == 1)
	{
		size_t items = fwrite(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

		if (items != n)
		{
			GSL_ERROR("fwrite failed", GSL_EFAILED);
		}
	}
	else
	{
		size_t i;

		for (i = 0; i < n; i++)
		{
			size_t item = fwrite(data + MULTIPLICITY * i * stride,
				MULTIPLICITY * sizeof(ATOMIC),
				1, stream);
			if (item != 1)
			{
				GSL_ERROR("fwrite failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)

int
FUNCTION(gsl_block, fprintf) (FILE * stream, const TYPE(gsl_block) * b, const char *format)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		int status;

		for (k = 0; k < MULTIPLICITY; k++)
		{
			if (k > 0)
			{
				status = putc(' ', stream);

				if (status == EOF)
				{
					GSL_ERROR("putc failed", GSL_EFAILED);
				}
			}
			status = fprintf(stream,
				format,
				data[MULTIPLICITY * i + k]);
			if (status < 0)
			{
				GSL_ERROR("fprintf failed", GSL_EFAILED);
			}
		}

		status = putc('\n', stream);

		if (status == EOF)
		{
			GSL_ERROR("putc failed", GSL_EFAILED);
		}
	}

	return 0;
}

int
FUNCTION(gsl_block, fscanf) (FILE * stream, TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		for (k = 0; k < MULTIPLICITY; k++)
		{
			ATOMIC_IO tmp;

			int status = fscanf(stream, IN_FORMAT, &tmp);

			data[MULTIPLICITY * i + k] = tmp;


			if (status != 1)
			{
				GSL_ERROR("fscanf failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}


int
FUNCTION(gsl_block, raw_fprintf) (FILE * stream,
const ATOMIC * data,
const size_t n,
const size_t stride,
const char *format)
{
	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		int status;

		for (k = 0; k < MULTIPLICITY; k++)
		{
			if (k > 0)
			{
				status = putc(' ', stream);

				if (status == EOF)
				{
					GSL_ERROR("putc failed", GSL_EFAILED);
				}
			}
			status = fprintf(stream,
				format,
				data[MULTIPLICITY * i * stride + k]);
			if (status < 0)
			{
				GSL_ERROR("fprintf failed", GSL_EFAILED);
			}
		}

		status = putc('\n', stream);

		if (status == EOF)
		{
			GSL_ERROR("putc failed", GSL_EFAILED);
		}
	}

	return 0;
}

int
FUNCTION(gsl_block, raw_fscanf) (FILE * stream,
ATOMIC * data,
const size_t n,
const size_t stride)
{
	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		for (k = 0; k < MULTIPLICITY; k++)
		{
			ATOMIC_IO tmp;

			int status = fscanf(stream, IN_FORMAT, &tmp);

			data[MULTIPLICITY * i * stride + k] = tmp;

			if (status != 1)
				GSL_ERROR("fscanf failed", GSL_EFAILED);
		}
	}

	return GSL_SUCCESS;
}

#endif

#include "templates_off.h"
#undef  BASE_UINT

#define BASE_INT
#include "templates_on.h"
int
FUNCTION(gsl_block, fread) (FILE * stream, TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t items = fread(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

	if (items != n)
	{
		GSL_ERROR("fread failed", GSL_EFAILED);
	}

	return 0;
}

int
FUNCTION(gsl_block, fwrite) (FILE * stream, const TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t items = fwrite(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

	if (items != n)
	{
		GSL_ERROR("fwrite failed", GSL_EFAILED);
	}

	return 0;
}

int
FUNCTION(gsl_block, raw_fread) (FILE * stream, ATOMIC * data,
const size_t n, const size_t stride)
{
	if (stride == 1)
	{
		size_t items = fread(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

		if (items != n)
		{
			GSL_ERROR("fread failed", GSL_EFAILED);
		}
	}
	else
	{
		size_t i;

		for (i = 0; i < n; i++)
		{
			size_t item = fread(data + MULTIPLICITY * i * stride,
				MULTIPLICITY * sizeof(ATOMIC), 1, stream);
			if (item != 1)
			{
				GSL_ERROR("fread failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_block, raw_fwrite) (FILE * stream, const ATOMIC * data,
const size_t n, const size_t stride)
{

	if (stride == 1)
	{
		size_t items = fwrite(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

		if (items != n)
		{
			GSL_ERROR("fwrite failed", GSL_EFAILED);
		}
	}
	else
	{
		size_t i;

		for (i = 0; i < n; i++)
		{
			size_t item = fwrite(data + MULTIPLICITY * i * stride,
				MULTIPLICITY * sizeof(ATOMIC),
				1, stream);
			if (item != 1)
			{
				GSL_ERROR("fwrite failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)

int
FUNCTION(gsl_block, fprintf) (FILE * stream, const TYPE(gsl_block) * b, const char *format)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		int status;

		for (k = 0; k < MULTIPLICITY; k++)
		{
			if (k > 0)
			{
				status = putc(' ', stream);

				if (status == EOF)
				{
					GSL_ERROR("putc failed", GSL_EFAILED);
				}
			}
			status = fprintf(stream,
				format,
				data[MULTIPLICITY * i + k]);
			if (status < 0)
			{
				GSL_ERROR("fprintf failed", GSL_EFAILED);
			}
		}

		status = putc('\n', stream);

		if (status == EOF)
		{
			GSL_ERROR("putc failed", GSL_EFAILED);
		}
	}

	return 0;
}

int
FUNCTION(gsl_block, fscanf) (FILE * stream, TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		for (k = 0; k < MULTIPLICITY; k++)
		{
			ATOMIC_IO tmp;

			int status = fscanf(stream, IN_FORMAT, &tmp);

			data[MULTIPLICITY * i + k] = tmp;


			if (status != 1)
			{
				GSL_ERROR("fscanf failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}


int
FUNCTION(gsl_block, raw_fprintf) (FILE * stream,
const ATOMIC * data,
const size_t n,
const size_t stride,
const char *format)
{
	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		int status;

		for (k = 0; k < MULTIPLICITY; k++)
		{
			if (k > 0)
			{
				status = putc(' ', stream);

				if (status == EOF)
				{
					GSL_ERROR("putc failed", GSL_EFAILED);
				}
			}
			status = fprintf(stream,
				format,
				data[MULTIPLICITY * i * stride + k]);
			if (status < 0)
			{
				GSL_ERROR("fprintf failed", GSL_EFAILED);
			}
		}

		status = putc('\n', stream);

		if (status == EOF)
		{
			GSL_ERROR("putc failed", GSL_EFAILED);
		}
	}

	return 0;
}

int
FUNCTION(gsl_block, raw_fscanf) (FILE * stream,
ATOMIC * data,
const size_t n,
const size_t stride)
{
	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		for (k = 0; k < MULTIPLICITY; k++)
		{
			ATOMIC_IO tmp;

			int status = fscanf(stream, IN_FORMAT, &tmp);

			data[MULTIPLICITY * i * stride + k] = tmp;

			if (status != 1)
				GSL_ERROR("fscanf failed", GSL_EFAILED);
		}
	}

	return GSL_SUCCESS;
}

#endif

#include "templates_off.h"
#undef  BASE_INT

#define BASE_USHORT
#include "templates_on.h"
int
FUNCTION(gsl_block, fread) (FILE * stream, TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t items = fread(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

	if (items != n)
	{
		GSL_ERROR("fread failed", GSL_EFAILED);
	}

	return 0;
}

int
FUNCTION(gsl_block, fwrite) (FILE * stream, const TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t items = fwrite(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

	if (items != n)
	{
		GSL_ERROR("fwrite failed", GSL_EFAILED);
	}

	return 0;
}

int
FUNCTION(gsl_block, raw_fread) (FILE * stream, ATOMIC * data,
const size_t n, const size_t stride)
{
	if (stride == 1)
	{
		size_t items = fread(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

		if (items != n)
		{
			GSL_ERROR("fread failed", GSL_EFAILED);
		}
	}
	else
	{
		size_t i;

		for (i = 0; i < n; i++)
		{
			size_t item = fread(data + MULTIPLICITY * i * stride,
				MULTIPLICITY * sizeof(ATOMIC), 1, stream);
			if (item != 1)
			{
				GSL_ERROR("fread failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_block, raw_fwrite) (FILE * stream, const ATOMIC * data,
const size_t n, const size_t stride)
{

	if (stride == 1)
	{
		size_t items = fwrite(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

		if (items != n)
		{
			GSL_ERROR("fwrite failed", GSL_EFAILED);
		}
	}
	else
	{
		size_t i;

		for (i = 0; i < n; i++)
		{
			size_t item = fwrite(data + MULTIPLICITY * i * stride,
				MULTIPLICITY * sizeof(ATOMIC),
				1, stream);
			if (item != 1)
			{
				GSL_ERROR("fwrite failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)

int
FUNCTION(gsl_block, fprintf) (FILE * stream, const TYPE(gsl_block) * b, const char *format)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		int status;

		for (k = 0; k < MULTIPLICITY; k++)
		{
			if (k > 0)
			{
				status = putc(' ', stream);

				if (status == EOF)
				{
					GSL_ERROR("putc failed", GSL_EFAILED);
				}
			}
			status = fprintf(stream,
				format,
				data[MULTIPLICITY * i + k]);
			if (status < 0)
			{
				GSL_ERROR("fprintf failed", GSL_EFAILED);
			}
		}

		status = putc('\n', stream);

		if (status == EOF)
		{
			GSL_ERROR("putc failed", GSL_EFAILED);
		}
	}

	return 0;
}

int
FUNCTION(gsl_block, fscanf) (FILE * stream, TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		for (k = 0; k < MULTIPLICITY; k++)
		{
			ATOMIC_IO tmp;

			int status = fscanf(stream, IN_FORMAT, &tmp);

			data[MULTIPLICITY * i + k] = tmp;


			if (status != 1)
			{
				GSL_ERROR("fscanf failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}


int
FUNCTION(gsl_block, raw_fprintf) (FILE * stream,
const ATOMIC * data,
const size_t n,
const size_t stride,
const char *format)
{
	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		int status;

		for (k = 0; k < MULTIPLICITY; k++)
		{
			if (k > 0)
			{
				status = putc(' ', stream);

				if (status == EOF)
				{
					GSL_ERROR("putc failed", GSL_EFAILED);
				}
			}
			status = fprintf(stream,
				format,
				data[MULTIPLICITY * i * stride + k]);
			if (status < 0)
			{
				GSL_ERROR("fprintf failed", GSL_EFAILED);
			}
		}

		status = putc('\n', stream);

		if (status == EOF)
		{
			GSL_ERROR("putc failed", GSL_EFAILED);
		}
	}

	return 0;
}

int
FUNCTION(gsl_block, raw_fscanf) (FILE * stream,
ATOMIC * data,
const size_t n,
const size_t stride)
{
	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		for (k = 0; k < MULTIPLICITY; k++)
		{
			ATOMIC_IO tmp;

			int status = fscanf(stream, IN_FORMAT, &tmp);

			data[MULTIPLICITY * i * stride + k] = tmp;

			if (status != 1)
				GSL_ERROR("fscanf failed", GSL_EFAILED);
		}
	}

	return GSL_SUCCESS;
}

#endif

#include "templates_off.h"
#undef  BASE_USHORT

#define BASE_SHORT
#include "templates_on.h"
int
FUNCTION(gsl_block, fread) (FILE * stream, TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t items = fread(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

	if (items != n)
	{
		GSL_ERROR("fread failed", GSL_EFAILED);
	}

	return 0;
}

int
FUNCTION(gsl_block, fwrite) (FILE * stream, const TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t items = fwrite(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

	if (items != n)
	{
		GSL_ERROR("fwrite failed", GSL_EFAILED);
	}

	return 0;
}

int
FUNCTION(gsl_block, raw_fread) (FILE * stream, ATOMIC * data,
const size_t n, const size_t stride)
{
	if (stride == 1)
	{
		size_t items = fread(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

		if (items != n)
		{
			GSL_ERROR("fread failed", GSL_EFAILED);
		}
	}
	else
	{
		size_t i;

		for (i = 0; i < n; i++)
		{
			size_t item = fread(data + MULTIPLICITY * i * stride,
				MULTIPLICITY * sizeof(ATOMIC), 1, stream);
			if (item != 1)
			{
				GSL_ERROR("fread failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_block, raw_fwrite) (FILE * stream, const ATOMIC * data,
const size_t n, const size_t stride)
{

	if (stride == 1)
	{
		size_t items = fwrite(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

		if (items != n)
		{
			GSL_ERROR("fwrite failed", GSL_EFAILED);
		}
	}
	else
	{
		size_t i;

		for (i = 0; i < n; i++)
		{
			size_t item = fwrite(data + MULTIPLICITY * i * stride,
				MULTIPLICITY * sizeof(ATOMIC),
				1, stream);
			if (item != 1)
			{
				GSL_ERROR("fwrite failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)

int
FUNCTION(gsl_block, fprintf) (FILE * stream, const TYPE(gsl_block) * b, const char *format)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		int status;

		for (k = 0; k < MULTIPLICITY; k++)
		{
			if (k > 0)
			{
				status = putc(' ', stream);

				if (status == EOF)
				{
					GSL_ERROR("putc failed", GSL_EFAILED);
				}
			}
			status = fprintf(stream,
				format,
				data[MULTIPLICITY * i + k]);
			if (status < 0)
			{
				GSL_ERROR("fprintf failed", GSL_EFAILED);
			}
		}

		status = putc('\n', stream);

		if (status == EOF)
		{
			GSL_ERROR("putc failed", GSL_EFAILED);
		}
	}

	return 0;
}

int
FUNCTION(gsl_block, fscanf) (FILE * stream, TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		for (k = 0; k < MULTIPLICITY; k++)
		{
			ATOMIC_IO tmp;

			int status = fscanf(stream, IN_FORMAT, &tmp);

			data[MULTIPLICITY * i + k] = tmp;


			if (status != 1)
			{
				GSL_ERROR("fscanf failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}


int
FUNCTION(gsl_block, raw_fprintf) (FILE * stream,
const ATOMIC * data,
const size_t n,
const size_t stride,
const char *format)
{
	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		int status;

		for (k = 0; k < MULTIPLICITY; k++)
		{
			if (k > 0)
			{
				status = putc(' ', stream);

				if (status == EOF)
				{
					GSL_ERROR("putc failed", GSL_EFAILED);
				}
			}
			status = fprintf(stream,
				format,
				data[MULTIPLICITY * i * stride + k]);
			if (status < 0)
			{
				GSL_ERROR("fprintf failed", GSL_EFAILED);
			}
		}

		status = putc('\n', stream);

		if (status == EOF)
		{
			GSL_ERROR("putc failed", GSL_EFAILED);
		}
	}

	return 0;
}

int
FUNCTION(gsl_block, raw_fscanf) (FILE * stream,
ATOMIC * data,
const size_t n,
const size_t stride)
{
	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		for (k = 0; k < MULTIPLICITY; k++)
		{
			ATOMIC_IO tmp;

			int status = fscanf(stream, IN_FORMAT, &tmp);

			data[MULTIPLICITY * i * stride + k] = tmp;

			if (status != 1)
				GSL_ERROR("fscanf failed", GSL_EFAILED);
		}
	}

	return GSL_SUCCESS;
}

#endif

#include "templates_off.h"
#undef  BASE_SHORT

#define BASE_UCHAR
#include "templates_on.h"
int
FUNCTION(gsl_block, fread) (FILE * stream, TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t items = fread(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

	if (items != n)
	{
		GSL_ERROR("fread failed", GSL_EFAILED);
	}

	return 0;
}

int
FUNCTION(gsl_block, fwrite) (FILE * stream, const TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t items = fwrite(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

	if (items != n)
	{
		GSL_ERROR("fwrite failed", GSL_EFAILED);
	}

	return 0;
}

int
FUNCTION(gsl_block, raw_fread) (FILE * stream, ATOMIC * data,
const size_t n, const size_t stride)
{
	if (stride == 1)
	{
		size_t items = fread(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

		if (items != n)
		{
			GSL_ERROR("fread failed", GSL_EFAILED);
		}
	}
	else
	{
		size_t i;

		for (i = 0; i < n; i++)
		{
			size_t item = fread(data + MULTIPLICITY * i * stride,
				MULTIPLICITY * sizeof(ATOMIC), 1, stream);
			if (item != 1)
			{
				GSL_ERROR("fread failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_block, raw_fwrite) (FILE * stream, const ATOMIC * data,
const size_t n, const size_t stride)
{

	if (stride == 1)
	{
		size_t items = fwrite(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

		if (items != n)
		{
			GSL_ERROR("fwrite failed", GSL_EFAILED);
		}
	}
	else
	{
		size_t i;

		for (i = 0; i < n; i++)
		{
			size_t item = fwrite(data + MULTIPLICITY * i * stride,
				MULTIPLICITY * sizeof(ATOMIC),
				1, stream);
			if (item != 1)
			{
				GSL_ERROR("fwrite failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)

int
FUNCTION(gsl_block, fprintf) (FILE * stream, const TYPE(gsl_block) * b, const char *format)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		int status;

		for (k = 0; k < MULTIPLICITY; k++)
		{
			if (k > 0)
			{
				status = putc(' ', stream);

				if (status == EOF)
				{
					GSL_ERROR("putc failed", GSL_EFAILED);
				}
			}
			status = fprintf(stream,
				format,
				data[MULTIPLICITY * i + k]);
			if (status < 0)
			{
				GSL_ERROR("fprintf failed", GSL_EFAILED);
			}
		}

		status = putc('\n', stream);

		if (status == EOF)
		{
			GSL_ERROR("putc failed", GSL_EFAILED);
		}
	}

	return 0;
}

int
FUNCTION(gsl_block, fscanf) (FILE * stream, TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		for (k = 0; k < MULTIPLICITY; k++)
		{
			ATOMIC_IO tmp;

			int status = fscanf(stream, IN_FORMAT, &tmp);

			data[MULTIPLICITY * i + k] = tmp;


			if (status != 1)
			{
				GSL_ERROR("fscanf failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}


int
FUNCTION(gsl_block, raw_fprintf) (FILE * stream,
const ATOMIC * data,
const size_t n,
const size_t stride,
const char *format)
{
	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		int status;

		for (k = 0; k < MULTIPLICITY; k++)
		{
			if (k > 0)
			{
				status = putc(' ', stream);

				if (status == EOF)
				{
					GSL_ERROR("putc failed", GSL_EFAILED);
				}
			}
			status = fprintf(stream,
				format,
				data[MULTIPLICITY * i * stride + k]);
			if (status < 0)
			{
				GSL_ERROR("fprintf failed", GSL_EFAILED);
			}
		}

		status = putc('\n', stream);

		if (status == EOF)
		{
			GSL_ERROR("putc failed", GSL_EFAILED);
		}
	}

	return 0;
}

int
FUNCTION(gsl_block, raw_fscanf) (FILE * stream,
ATOMIC * data,
const size_t n,
const size_t stride)
{
	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		for (k = 0; k < MULTIPLICITY; k++)
		{
			ATOMIC_IO tmp;

			int status = fscanf(stream, IN_FORMAT, &tmp);

			data[MULTIPLICITY * i * stride + k] = tmp;

			if (status != 1)
				GSL_ERROR("fscanf failed", GSL_EFAILED);
		}
	}

	return GSL_SUCCESS;
}

#endif

#include "templates_off.h"
#undef  BASE_UCHAR

#define BASE_CHAR
#include "templates_on.h"
int
FUNCTION(gsl_block, fread) (FILE * stream, TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t items = fread(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

	if (items != n)
	{
		GSL_ERROR("fread failed", GSL_EFAILED);
	}

	return 0;
}

int
FUNCTION(gsl_block, fwrite) (FILE * stream, const TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t items = fwrite(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

	if (items != n)
	{
		GSL_ERROR("fwrite failed", GSL_EFAILED);
	}

	return 0;
}

int
FUNCTION(gsl_block, raw_fread) (FILE * stream, ATOMIC * data,
const size_t n, const size_t stride)
{
	if (stride == 1)
	{
		size_t items = fread(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

		if (items != n)
		{
			GSL_ERROR("fread failed", GSL_EFAILED);
		}
	}
	else
	{
		size_t i;

		for (i = 0; i < n; i++)
		{
			size_t item = fread(data + MULTIPLICITY * i * stride,
				MULTIPLICITY * sizeof(ATOMIC), 1, stream);
			if (item != 1)
			{
				GSL_ERROR("fread failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}

int
FUNCTION(gsl_block, raw_fwrite) (FILE * stream, const ATOMIC * data,
const size_t n, const size_t stride)
{

	if (stride == 1)
	{
		size_t items = fwrite(data, MULTIPLICITY * sizeof(ATOMIC), n, stream);

		if (items != n)
		{
			GSL_ERROR("fwrite failed", GSL_EFAILED);
		}
	}
	else
	{
		size_t i;

		for (i = 0; i < n; i++)
		{
			size_t item = fwrite(data + MULTIPLICITY * i * stride,
				MULTIPLICITY * sizeof(ATOMIC),
				1, stream);
			if (item != 1)
			{
				GSL_ERROR("fwrite failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)

int
FUNCTION(gsl_block, fprintf) (FILE * stream, const TYPE(gsl_block) * b, const char *format)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		int status;

		for (k = 0; k < MULTIPLICITY; k++)
		{
			if (k > 0)
			{
				status = putc(' ', stream);

				if (status == EOF)
				{
					GSL_ERROR("putc failed", GSL_EFAILED);
				}
			}
			status = fprintf(stream,
				format,
				data[MULTIPLICITY * i + k]);
			if (status < 0)
			{
				GSL_ERROR("fprintf failed", GSL_EFAILED);
			}
		}

		status = putc('\n', stream);

		if (status == EOF)
		{
			GSL_ERROR("putc failed", GSL_EFAILED);
		}
	}

	return 0;
}

int
FUNCTION(gsl_block, fscanf) (FILE * stream, TYPE(gsl_block) * b)
{
	size_t n = b->size;

	ATOMIC * data = b->data;

	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		for (k = 0; k < MULTIPLICITY; k++)
		{
			ATOMIC_IO tmp;

			int status = fscanf(stream, IN_FORMAT, &tmp);

			data[MULTIPLICITY * i + k] = tmp;


			if (status != 1)
			{
				GSL_ERROR("fscanf failed", GSL_EFAILED);
			}
		}
	}

	return GSL_SUCCESS;
}


int
FUNCTION(gsl_block, raw_fprintf) (FILE * stream,
const ATOMIC * data,
const size_t n,
const size_t stride,
const char *format)
{
	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		int status;

		for (k = 0; k < MULTIPLICITY; k++)
		{
			if (k > 0)
			{
				status = putc(' ', stream);

				if (status == EOF)
				{
					GSL_ERROR("putc failed", GSL_EFAILED);
				}
			}
			status = fprintf(stream,
				format,
				data[MULTIPLICITY * i * stride + k]);
			if (status < 0)
			{
				GSL_ERROR("fprintf failed", GSL_EFAILED);
			}
		}

		status = putc('\n', stream);

		if (status == EOF)
		{
			GSL_ERROR("putc failed", GSL_EFAILED);
		}
	}

	return 0;
}

int
FUNCTION(gsl_block, raw_fscanf) (FILE * stream,
ATOMIC * data,
const size_t n,
const size_t stride)
{
	size_t i;

	for (i = 0; i < n; i++)
	{
		int k;
		for (k = 0; k < MULTIPLICITY; k++)
		{
			ATOMIC_IO tmp;

			int status = fscanf(stream, IN_FORMAT, &tmp);

			data[MULTIPLICITY * i * stride + k] = tmp;

			if (status != 1)
				GSL_ERROR("fscanf failed", GSL_EFAILED);
		}
	}

	return GSL_SUCCESS;
}

#endif

#include "templates_off.h"
#undef  BASE_CHAR
