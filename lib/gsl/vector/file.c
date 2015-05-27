#include "stdafx.h"
#include "config.h.in"
#include <stdio.h>
#include <gsl/err/gsl_errno.h>
#include <gsl/block/gsl_block.h>
#include <gsl/vector/gsl_vector.h>

#define BASE_GSL_COMPLEX_LONG
#include "../templates_on.h"
int
FUNCTION(gsl_vector, fread) (FILE * stream, TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fread) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}

int
FUNCTION(gsl_vector, fwrite) (FILE * stream, const TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fwrite) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)
int
FUNCTION(gsl_vector, fprintf) (FILE * stream, const TYPE(gsl_vector) * v,
const char *format)
{
	int status = FUNCTION(gsl_block, raw_fprintf) (stream,
		v->data,
		v->size,
		v->stride,
		format);
	return status;
}

int
FUNCTION(gsl_vector, fscanf) (FILE * stream, TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fscanf) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}
#endif
#include "../templates_off.h"
#undef  BASE_GSL_COMPLEX_LONG

#define BASE_GSL_COMPLEX
#include "../templates_on.h"
int
FUNCTION(gsl_vector, fread) (FILE * stream, TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fread) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}

int
FUNCTION(gsl_vector, fwrite) (FILE * stream, const TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fwrite) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)
int
FUNCTION(gsl_vector, fprintf) (FILE * stream, const TYPE(gsl_vector) * v,
const char *format)
{
	int status = FUNCTION(gsl_block, raw_fprintf) (stream,
		v->data,
		v->size,
		v->stride,
		format);
	return status;
}

int
FUNCTION(gsl_vector, fscanf) (FILE * stream, TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fscanf) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}
#endif
#include "../templates_off.h"
#undef  BASE_GSL_COMPLEX

#define BASE_GSL_COMPLEX_FLOAT
#include "../templates_on.h"
int
FUNCTION(gsl_vector, fread) (FILE * stream, TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fread) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}

int
FUNCTION(gsl_vector, fwrite) (FILE * stream, const TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fwrite) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)
int
FUNCTION(gsl_vector, fprintf) (FILE * stream, const TYPE(gsl_vector) * v,
const char *format)
{
	int status = FUNCTION(gsl_block, raw_fprintf) (stream,
		v->data,
		v->size,
		v->stride,
		format);
	return status;
}

int
FUNCTION(gsl_vector, fscanf) (FILE * stream, TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fscanf) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}
#endif
#include "../templates_off.h"
#undef  BASE_GSL_COMPLEX_FLOAT

#define BASE_LONG_DOUBLE
#include "../templates_on.h"
int
FUNCTION(gsl_vector, fread) (FILE * stream, TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fread) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}

int
FUNCTION(gsl_vector, fwrite) (FILE * stream, const TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fwrite) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)
int
FUNCTION(gsl_vector, fprintf) (FILE * stream, const TYPE(gsl_vector) * v,
const char *format)
{
	int status = FUNCTION(gsl_block, raw_fprintf) (stream,
		v->data,
		v->size,
		v->stride,
		format);
	return status;
}

int
FUNCTION(gsl_vector, fscanf) (FILE * stream, TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fscanf) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}
#endif
#include "../templates_off.h"
#undef  BASE_LONG_DOUBLE

#define BASE_DOUBLE
#include "../templates_on.h"
int
FUNCTION(gsl_vector, fread) (FILE * stream, TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fread) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}

int
FUNCTION(gsl_vector, fwrite) (FILE * stream, const TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fwrite) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)
int
FUNCTION(gsl_vector, fprintf) (FILE * stream, const TYPE(gsl_vector) * v,
const char *format)
{
	int status = FUNCTION(gsl_block, raw_fprintf) (stream,
		v->data,
		v->size,
		v->stride,
		format);
	return status;
}

int
FUNCTION(gsl_vector, fscanf) (FILE * stream, TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fscanf) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}
#endif
#include "../templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "../templates_on.h"
int
FUNCTION(gsl_vector, fread) (FILE * stream, TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fread) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}

int
FUNCTION(gsl_vector, fwrite) (FILE * stream, const TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fwrite) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)
int
FUNCTION(gsl_vector, fprintf) (FILE * stream, const TYPE(gsl_vector) * v,
const char *format)
{
	int status = FUNCTION(gsl_block, raw_fprintf) (stream,
		v->data,
		v->size,
		v->stride,
		format);
	return status;
}

int
FUNCTION(gsl_vector, fscanf) (FILE * stream, TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fscanf) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}
#endif
#include "../templates_off.h"
#undef  BASE_FLOAT

#define BASE_ULONG
#include "../templates_on.h"
int
FUNCTION(gsl_vector, fread) (FILE * stream, TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fread) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}

int
FUNCTION(gsl_vector, fwrite) (FILE * stream, const TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fwrite) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)
int
FUNCTION(gsl_vector, fprintf) (FILE * stream, const TYPE(gsl_vector) * v,
const char *format)
{
	int status = FUNCTION(gsl_block, raw_fprintf) (stream,
		v->data,
		v->size,
		v->stride,
		format);
	return status;
}

int
FUNCTION(gsl_vector, fscanf) (FILE * stream, TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fscanf) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}
#endif
#include "../templates_off.h"
#undef  BASE_ULONG

#define BASE_LONG
#include "../templates_on.h"
int
FUNCTION(gsl_vector, fread) (FILE * stream, TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fread) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}

int
FUNCTION(gsl_vector, fwrite) (FILE * stream, const TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fwrite) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)
int
FUNCTION(gsl_vector, fprintf) (FILE * stream, const TYPE(gsl_vector) * v,
const char *format)
{
	int status = FUNCTION(gsl_block, raw_fprintf) (stream,
		v->data,
		v->size,
		v->stride,
		format);
	return status;
}

int
FUNCTION(gsl_vector, fscanf) (FILE * stream, TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fscanf) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}
#endif
#include "../templates_off.h"
#undef  BASE_LONG

#define BASE_UINT
#include "../templates_on.h"
int
FUNCTION(gsl_vector, fread) (FILE * stream, TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fread) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}

int
FUNCTION(gsl_vector, fwrite) (FILE * stream, const TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fwrite) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)
int
FUNCTION(gsl_vector, fprintf) (FILE * stream, const TYPE(gsl_vector) * v,
const char *format)
{
	int status = FUNCTION(gsl_block, raw_fprintf) (stream,
		v->data,
		v->size,
		v->stride,
		format);
	return status;
}

int
FUNCTION(gsl_vector, fscanf) (FILE * stream, TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fscanf) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}
#endif
#include "../templates_off.h"
#undef  BASE_UINT

#define BASE_INT
#include "../templates_on.h"
int
FUNCTION(gsl_vector, fread) (FILE * stream, TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fread) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}

int
FUNCTION(gsl_vector, fwrite) (FILE * stream, const TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fwrite) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)
int
FUNCTION(gsl_vector, fprintf) (FILE * stream, const TYPE(gsl_vector) * v,
const char *format)
{
	int status = FUNCTION(gsl_block, raw_fprintf) (stream,
		v->data,
		v->size,
		v->stride,
		format);
	return status;
}

int
FUNCTION(gsl_vector, fscanf) (FILE * stream, TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fscanf) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}
#endif
#include "../templates_off.h"
#undef  BASE_INT

#define BASE_USHORT
#include "../templates_on.h"
int
FUNCTION(gsl_vector, fread) (FILE * stream, TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fread) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}

int
FUNCTION(gsl_vector, fwrite) (FILE * stream, const TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fwrite) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)
int
FUNCTION(gsl_vector, fprintf) (FILE * stream, const TYPE(gsl_vector) * v,
const char *format)
{
	int status = FUNCTION(gsl_block, raw_fprintf) (stream,
		v->data,
		v->size,
		v->stride,
		format);
	return status;
}

int
FUNCTION(gsl_vector, fscanf) (FILE * stream, TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fscanf) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}
#endif
#include "../templates_off.h"
#undef  BASE_USHORT

#define BASE_SHORT
#include "../templates_on.h"
int
FUNCTION(gsl_vector, fread) (FILE * stream, TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fread) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}

int
FUNCTION(gsl_vector, fwrite) (FILE * stream, const TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fwrite) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)
int
FUNCTION(gsl_vector, fprintf) (FILE * stream, const TYPE(gsl_vector) * v,
const char *format)
{
	int status = FUNCTION(gsl_block, raw_fprintf) (stream,
		v->data,
		v->size,
		v->stride,
		format);
	return status;
}

int
FUNCTION(gsl_vector, fscanf) (FILE * stream, TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fscanf) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}
#endif
#include "../templates_off.h"
#undef  BASE_SHORT

#define BASE_UCHAR
#include "../templates_on.h"
int
FUNCTION(gsl_vector, fread) (FILE * stream, TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fread) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}

int
FUNCTION(gsl_vector, fwrite) (FILE * stream, const TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fwrite) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)
int
FUNCTION(gsl_vector, fprintf) (FILE * stream, const TYPE(gsl_vector) * v,
const char *format)
{
	int status = FUNCTION(gsl_block, raw_fprintf) (stream,
		v->data,
		v->size,
		v->stride,
		format);
	return status;
}

int
FUNCTION(gsl_vector, fscanf) (FILE * stream, TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fscanf) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}
#endif
#include "../templates_off.h"
#undef  BASE_UCHAR

#define BASE_CHAR
#include "../templates_on.h"
int
FUNCTION(gsl_vector, fread) (FILE * stream, TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fread) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}

int
FUNCTION(gsl_vector, fwrite) (FILE * stream, const TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fwrite) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}

#if !(USES_LONGDOUBLE && !HAVE_PRINTF_LONGDOUBLE)
int
FUNCTION(gsl_vector, fprintf) (FILE * stream, const TYPE(gsl_vector) * v,
const char *format)
{
	int status = FUNCTION(gsl_block, raw_fprintf) (stream,
		v->data,
		v->size,
		v->stride,
		format);
	return status;
}

int
FUNCTION(gsl_vector, fscanf) (FILE * stream, TYPE(gsl_vector) * v)
{
	int status = FUNCTION(gsl_block, raw_fscanf) (stream,
		v->data,
		v->size,
		v->stride);
	return status;
}
#endif
#include "../templates_off.h"
#undef  BASE_CHAR
