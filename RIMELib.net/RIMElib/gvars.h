
/*
 *  RIMElib: RuntIme Mathematical Equation Library (data types)
 *  Copyright (C) 2004  Alistair Shilton
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */


#ifndef _gvars_h
#define _gvars_h

/*********************************************************************

Generic number structure
========================

The aim of this module is to provide a dynamically self-typing data
type for mathematics.  Thus a variable of type GenVar can be 0,1, an
integer, a real number, +/- infinity, finite indeterminant, infinite
indeterminant, an error code or even "unknown" (ie. the answer depends
on information we don't have), depending on the situation.  This type
will change dynamically as operations are carried out on the variable,
always taking the simplest possible appropriate form.

In detail, possible forms the data may take are:

DATA_IS_ZERO:      value is zero.
DATA_IS_ONE:       value is one.
DATA_IS_INTEGER:   value is an integer (type GVAR_INT), IntegerValue.
DATA_IS_REAL:      value is a real number (type double), RealValue.
DATA_IS_POS_INFTY: value is positive infinity.
DATA_IS_NEG_INFTY: value is negative infinity.
DATA_IS_FIN_INDET: value is indeterminate, but known to be finite.
DATA_IS_INF_INDET: value is indeterminate and may be infinite.
DATA_IS_ERROR:     error has occured, type ErrorCode.
DATA_IS_UNCALCED:  value cannot be determined, information lacking.



Error handling
==============

If an error occurs, the datatype will be set to DATA_IS_ERROR.  The
error code is ErrorCode - each bit of this is a specific error, so
multiple error codes can be defined by ORing different errors.  Errors
are defined below, and have the form GVAR_.



Constructors
============

The following constructors are to be used for the various types:

GenVar ResultZero(void)
GenVar ResultOne(void)
GenVar ResultInteger(GVAR_INT what)
GenVar ResultReal(double what)
GenVar ResultPosInfnty(void)
GenVar ResultNegInfnty(void)
GenVar ResultFinIndet(void)
GenVar ResultInfIndet(void)
GenVar ResultError(int errorCode)
GenVar ResultUncalced(void)

NB: integer and real type constructors do bounds checking on the
    number, so there is no need to worry about NaNs and infs that
    might be given to these constructors.


Data conversion functions
=========================

The datatype can be cast to double and GVAR_INT types using the
following functions:

double getReal(GenVar a)
GVAR_INT getInteger(GenVar a)

Notes: - getReal will return NaN if conversion is not possible.
       - getInteger will return GVAR_INT_MAX if conver not possible.
       - getInteger will return GVAR_INT_MAX-1 if result is inf or too
         large for type GVAR_INT, or greater than GVAR_INT_MAX-1.
       - getInteger will return GVAR_INT_MIN-1 if result is -inf or
         too small for type GVAR_INT, or less than GVAR_INT_MIN+1.


Type testing functions
======================

The following functions allow testing of data type.  They return 0
(false) unless given conditions are met:

int isFiniteReal(a):    returns 1 if number is finite real.
int isFiniteInteger(a): returns 1 if number is finite integer.

int isReal(a):    returns 1 if number is real, including infinite.
int isInteger(a): returns 1 if number is integer, including infinite.

int isError(a):    returns 1 if number indicates error.
int isUncalced(a): returns 1 if number is unknown.

NB: Indeterminant numbers are classed as real (which they are - they
    just aren't known).  They are not integers.


Data comparison functions
=========================

The following relational operators are defined.  They return 0 (false)
unless given conditions are met:

int isaequalb(a,b): returns 1 if a == b, a and b real.

int isalessb(a,b):         returns 1 if a <  b, a and b real.
int isagreaterb(a,b):      returns 1 if a >  b, a and b real.
int isalessequalb(a,b):    returns 1 if a <= b, a and b real.
int isagreaterequalb(a,b): returns 1 if a >= b, a and b real.

Notes: - if the comparison is meaningless (eg. a is an error, or
         uncalced) then 0 (false) will be returned.
       - we assume: inf == inf, -inf == -inf, inf != -inf, inf !< inf,
         -inf !< -inf, -inf < inf, inf !< -inf.
       - remember that data is simplified.  So, if a is of type real
         then we can assume it isn't an integer, and likewise if a is
         an integer we can assume that it isn't 0 or 1.


Basic arithmetic operations
===========================

The unary additive and multiplicative inversion operators are defined
as follows:

GenVar negGenVar(GenVar a): additive inversion (negation)
GenVar invGenVar(GenVar a): multiplicative inversion

Also, we have the usual binary operations:

GenVar addGenVar(GenVar a, GenVar b)  - addition a+b
GenVar subGenVar(GenVar a, GenVar b)  - subtraction a-b
GenVar mulGenVar(GenVar a, GenVar b)  - multiplication a*b
GenVar divGenVar(GenVar a, GenVar b)  - division a/b
GenVar idivGenVar(GenVar a, GenVar b) - integer only division a/b
GenVar modGenVar(GenVar a, GenVar b)  - mod a%b

Notes: - If one or both of a, b is an error then the result will copy
         this error (if both are errors, it will combine both errors).
       - If one is an error and the other uncalced, the result will be
         an error.
       - inf/0 = inf
       - -inf/0 = -inf
       - the inverse of 0 is given as inf (we assume inv(a) = 1/a).
         One side-effect of this is that -inf, inverted twice, will
         give inf, not -inf as may be expected (ie. 0+ and 0- not
         supported).
       - i%0 = 0 by definition here
       - inf%j assumes that inf is integer inf.  Hence inf%0 = 0,
         -inf%0 = 0, inf%i = i, -inf%i = -i.
       - i%inf = i%-inf = i by definition here
       - divGenVar acting on two integers will, in general, give a
         real result.  use idivGenVar to indicate strict integer
         division.


Number parsing and deparsing
============================

Number format for parsing/deparsing is thus:

- A GVAR_INT conforming to the usual c printf standards "%L".
- A double conforming to the usual c printf standards "%1.15le".
- "inf" to indicate infinity.
- "-inf" to indicate negative infinity.
- "FiniteIndet" to indicate a finite indefinite number.
- "InfintIndet" to indicate an infinite indefinite number.
- "t" followed by a number to indicate an error.
- "uc" to indicated uncalculated variable.

The following functions are provided for parsing:

GenVar parseNumber(char *express, long len)
char *deparseNumber(GenVar a)
char *deparseNumber_cstruct(GenVar a)

Upon failure, deparse will return NULL.


Helper Functons
===============

GenVar combineErrors(GenVar a, GenVar b): assuming both a and b are
     errors, this will return an error with the error types
     appropriately combined (by ORing).

int combineDataDescr(GenVar a, GenVar b): Returns the
     combined types of a and b (in order a_b).  This is useful for
     functions of two variables.


Overflow safe binary operations
===============================

To prevent overflow problems, integer and real addition and
multiplication should be done by passing the relevant ints or doubles
to the following functions.  These do type-safe add or mult and
returns the result in GenVar form.  An integer overflow will return a
real result (no overflow, some rounding possible) - no overflow will
return integer type.  A real overflow will give positive or negative
infinity, which may cause problems.

GenVar safeIntAdd(GVAR_INT a, GVAR_INT b) - overflow safe int add.
GenVar safeIntMul(GVAR_INT a, GVAR_INT b) - overflow safe int mul.

GenVar safeRealAdd(double a, double b) - overflow safe real addition
GenVar safeRealMul(double a, double b) - overflow safe real multiplication

*********************************************************************/

#include <limits.h>

#ifndef LONG_LONG_MAX
#define LONG_LONG_MAX LLONG_MAX
#endif
#ifndef LONG_LONG_MIN
#define LONG_LONG_MIN LLONG_MIN
#endif

#define GVAR_INT        long long
#define GVAR_INT_MAX    LONG_LONG_MAX
#define GVAR_INT_MIN    LONG_LONG_MIN

/*
   I was using enum here originally, but C++ refused to grep ints to enums,
   and string conversion was too painful, so I've reverted to good old
   macros.
*/

#define DATA_IS_ZERO            0
#define DATA_IS_ONE             1
#define DATA_IS_INTEGER         2
#define DATA_IS_REAL            3
#define DATA_IS_POS_INFTY       4
#define DATA_IS_NEG_INFTY       5
#define DATA_IS_FIN_INDET       6
#define DATA_IS_INF_INDET       7
#define DATA_IS_ERROR           8
#define DATA_IS_UNCALCED        9

#define DATA_ZERO_ZERO                  0
#define DATA_ONE_ZERO                   1
#define DATA_INTEGER_ZERO               2
#define DATA_REAL_ZERO                  3
#define DATA_POS_INFTY_ZERO             4
#define DATA_NEG_INFTY_ZERO             5
#define DATA_FIN_INDET_ZERO             6
#define DATA_INF_INDET_ZERO             7
#define DATA_ERROR_ZERO                 8
#define DATA_UNCALCED_ZERO              9
#define DATA_ZERO_ONE                   10
#define DATA_ONE_ONE                    11
#define DATA_INTEGER_ONE                12
#define DATA_REAL_ONE                   13
#define DATA_POS_INFTY_ONE              14
#define DATA_NEG_INFTY_ONE              15
#define DATA_FIN_INDET_ONE              16
#define DATA_INF_INDET_ONE              17
#define DATA_ERROR_ONE                  18
#define DATA_UNCALCED_ONE               19
#define DATA_ZERO_INTEGER               20
#define DATA_ONE_INTEGER                21
#define DATA_INTEGER_INTEGER            22
#define DATA_REAL_INTEGER               23
#define DATA_POS_INFTY_INTEGER          24
#define DATA_NEG_INFTY_INTEGER          25
#define DATA_FIN_INDET_INTEGER          26
#define DATA_INF_INDET_INTEGER          27
#define DATA_ERROR_INTEGER              28
#define DATA_UNCALCED_INTEGER           29
#define DATA_ZERO_REAL                  30
#define DATA_ONE_REAL                   31
#define DATA_INTEGER_REAL               32
#define DATA_REAL_REAL                  33
#define DATA_POS_INFTY_REAL             34
#define DATA_NEG_INFTY_REAL             35
#define DATA_FIN_INDET_REAL             36
#define DATA_INF_INDET_REAL             37
#define DATA_ERROR_REAL                 38
#define DATA_UNCALCED_REAL              39
#define DATA_ZERO_POS_INFTY             40
#define DATA_ONE_POS_INFTY              41
#define DATA_INTEGER_POS_INFTY          42
#define DATA_REAL_POS_INFTY             43
#define DATA_POS_INFTY_POS_INFTY        44
#define DATA_NEG_INFTY_POS_INFTY        45
#define DATA_FIN_INDET_POS_INFTY        46
#define DATA_INF_INDET_POS_INFTY        47
#define DATA_ERROR_POS_INFTY            48
#define DATA_UNCALCED_POS_INFTY         49
#define DATA_ZERO_NEG_INFTY             50
#define DATA_ONE_NEG_INFTY              51
#define DATA_INTEGER_NEG_INFTY          52
#define DATA_REAL_NEG_INFTY             53
#define DATA_POS_INFTY_NEG_INFTY        54
#define DATA_NEG_INFTY_NEG_INFTY        55
#define DATA_FIN_INDET_NEG_INFTY        56
#define DATA_INF_INDET_NEG_INFTY        57
#define DATA_ERROR_NEG_INFTY            58
#define DATA_UNCALCED_NEG_INFTY         59
#define DATA_ZERO_FIN_INDET             60
#define DATA_ONE_FIN_INDET              61
#define DATA_INTEGER_FIN_INDET          62
#define DATA_REAL_FIN_INDET             63
#define DATA_POS_INFTY_FIN_INDET        64
#define DATA_NEG_INFTY_FIN_INDET        65
#define DATA_FIN_INDET_FIN_INDET        66
#define DATA_INF_INDET_FIN_INDET        67
#define DATA_ERROR_FIN_INDET            68
#define DATA_UNCALCED_FIN_INDET         69
#define DATA_ZERO_INF_INDET             70
#define DATA_ONE_INF_INDET              71
#define DATA_INTEGER_INF_INDET          72
#define DATA_REAL_INF_INDET             73
#define DATA_POS_INFTY_INF_INDET        74
#define DATA_NEG_INFTY_INF_INDET        75
#define DATA_FIN_INDET_INF_INDET        76
#define DATA_INF_INDET_INF_INDET        77
#define DATA_ERROR_INF_INDET            78
#define DATA_UNCALCED_INF_INDET         79
#define DATA_ZERO_ERROR                 80
#define DATA_ONE_ERROR                  81
#define DATA_INTEGER_ERROR              82
#define DATA_REAL_ERROR                 83
#define DATA_POS_INFTY_ERROR            84
#define DATA_NEG_INFTY_ERROR            85
#define DATA_FIN_INDET_ERROR            86
#define DATA_INF_INDET_ERROR            87
#define DATA_ERROR_ERROR                88
#define DATA_UNCALCED_ERROR             89
#define DATA_ZERO_UNCALCED              90
#define DATA_ONE_UNCALCED               91
#define DATA_INTEGER_UNCALCED           92
#define DATA_REAL_UNCALCED              93
#define DATA_POS_INFTY_UNCALCED         94
#define DATA_NEG_INFTY_UNCALCED         95
#define DATA_FIN_INDET_UNCALCED         96
#define DATA_INF_INDET_UNCALCED         97
#define DATA_ERROR_UNCALCED             98
#define DATA_UNCALCED_UNCALCED          99



#define GVAR_EUNKNOWN   0x00000000000000001 /* Unknown error */
#define GVAR_ENANDET    0x00000000000000002 /* NaN detected during simplification of real */
#define GVAR_ESYNTAX    0x00000000000000004 /* Syntax error */
#define GVAR_ETYPEERR   0x00000000000000008 /* Datatype unknown error */
#define GVAR_EBADLOCAL  0x00000000000000010 /* Bad local variable */
#define GVAR_EILLINTEG  0x00000000000000020 /* Ill-defined (backward) integral step size */
#define GVAR_c          0x00000000000000040
#define GVAR_d          0x00000000000000080
#define GVAR_e          0x00000000000000100
#define GVAR_f          0x00000000000000200
#define GVAR_g          0x00000000000000400
#define GVAR_h          0x00000000000000800
#define GVAR_i          0x00000000000001000
#define GVAR_j          0x00000000000002000
#define GVAR_k          0x00000000000004000
#define GVAR_l          0x00000000000008000
#define GVAR_m          0x00000000000010000
#define GVAR_n          0x00000000000020000
#define GVAR_o          0x00000000000040000
#define GVAR_p          0x00000000000080000
#define GVAR_q          0x00000000000100000
#define GVAR_r          0x00000000000200000
#define GVAR_s          0x00000000000400000
#define GVAR_t          0x00000000000800000
#define GVAR_GSL_FAILUR 0x00000000001000000 /* gsl: generic failure */
#define GVAR_GSL_CONTIN 0x00000000002000000 /* gsl: iteration has not converged */
#define GVAR_GSL_EDOM   0x00000000004000000 /* gsl: input domain error  e.g sqrt(-1) */
#define GVAR_GSL_ERANGE 0x00000000008000000 /* gsl: output range error, e.g. exp(1e100) */
#define GVAR_GSL_EFAULT 0x00000000010000000 /* gsl: invalid pointer */
#define GVAR_GSL_EINVAL 0x00000000020000000 /* gsl: invalid argument supplied by user */
#define GVAR_GSL_EFAILE 0x00000000040000000 /* gsl: generic failure */
#define GVAR_GSL_EFACTO 0x00000000080000000 /* gsl: factorization failed */
#define GVAR_GSL_ESANIT 0x00000000100000000 /* gsl: sanity check failed - shouldn't happen */
#define GVAR_GSL_ENOMEM 0x00000000200000000 /* gsl: malloc failed */
#define GVAR_GSL_EBADFN 0x00000000400000000 /* gsl: problem with user-supplied function */
#define GVAR_GSL_ERNWAY 0x00000000800000000 /* gsl: iterative process is out of control */
#define GVAR_GSL_EMAXIT 0x00000001000000000 /* gsl: exceeded max number of iterations */
#define GVAR_GSL_EZRDIV 0x00000002000000000 /* gsl: tried to divide by zero */
#define GVAR_GSL_EBADTL 0x00000004000000000 /* gsl: user specified an invalid tolerance */
#define GVAR_GSL_ETOL   0x00000008000000000 /* gsl: failed to reach the specified tolerance */
#define GVAR_GSL_EUNDRF 0x00000010000000000 /* gsl: underflow */
#define GVAR_GSL_EOVRFL 0x00000020000000000 /* gsl: overflow  */
#define GVAR_GSL_ELOSS  0x00000040000000000 /* gsl: loss of accuracy */
#define GVAR_GSL_EROUND 0x00000080000000000 /* gsl: failed because of roundoff error */
#define GVAR_GSL_EBADLE 0x00000100000000000 /* gsl: matrix, vector lengths are not conformant */
#define GVAR_GSL_ENOTSQ 0x00000200000000000 /* gsl: matrix not square */
#define GVAR_GSL_ESING  0x00000400000000000 /* gsl: apparent singularity detected */
#define GVAR_GSL_EDIVER 0x00000800000000000 /* gsl: integral or series is divergent */
#define GVAR_GSL_EUNSUP 0x00001000000000000 /* gsl: requested feature is not supported by the hardware */
#define GVAR_GSL_EUNIMP 0x00002000000000000 /* gsl: requested feature not (yet) implemented */
#define GVAR_GSL_ECACHE 0x00004000000000000 /* gsl: cache limit exceeded */
#define GVAR_GSL_ETABLE 0x00008000000000000 /* gsl: table limit exceeded */
#define GVAR_GSL_ENOPRO 0x00010000000000000 /* gsl: iteration is not making progress towards solution */
#define GVAR_GSL_ENOPRX 0x00020000000000000 /* gsl: jacobian evaluations are not improving the solution */
#define GVAR_GSL_ETOLF  0x00040000000000000 /* gsl: cannot reach the specified tolerance in F */
#define GVAR_GSL_ETOLX  0x00080000000000000 /* gsl: cannot reach the specified tolerance in X */
#define GVAR_GSL_ETOLG  0x00100000000000000 /* gsl: cannot reach the specified tolerance in gradient */
#define GVAR_GSL_EOF    0x00200000000000000 /* gsl: end of file */
#define GVAR_u          0x00400000000000000
#define GVAR_v          0x00800000000000000
#define GVAR_w          0x01000000000000000
#define GVAR_x          0x02000000000000000
#define GVAR_y          0x04000000000000000


typedef struct
{
    int DataType;

    GVAR_INT           IntegerValue;
    double             RealValue;
    unsigned long long ErrorCode;
}
GenVar;


GenVar ResultZero(void);
GenVar ResultOne(void);
GenVar ResultInteger(GVAR_INT what);
GenVar ResultReal(double what);
GenVar ResultPosInfnty(void);
GenVar ResultNegInfnty(void);
GenVar ResultFinIndet(void);
GenVar ResultInfIndet(void);
GenVar ResultError(unsigned long long errorCode);
GenVar ResultUncalced(void);


double getReal(GenVar a);
GVAR_INT getInteger(GenVar a);


int isFiniteReal(GenVar a);
int isFiniteInteger(GenVar a);

int isReal(GenVar a);
int isInteger(GenVar a);

int isError(GenVar a);
int isUncalced(GenVar a);


int isaequalb(GenVar a, GenVar b);

int isalessb(GenVar a, GenVar b);
int isagreaterb(GenVar a, GenVar b);
int isalessequalb(GenVar a, GenVar b);
int isagreaterequalb(GenVar a, GenVar b);


GenVar addGenVar(GenVar a, GenVar b);
GenVar subGenVar(GenVar a, GenVar b);
GenVar mulGenVar(GenVar a, GenVar b);
GenVar divGenVar(GenVar a, GenVar b);
GenVar idivGenVar(GenVar a, GenVar b);
GenVar modGenVar(GenVar a, GenVar b);
GenVar negGenVar(GenVar a);
GenVar invGenVar(GenVar a);


char *deparseNumber(GenVar a);
GenVar parseNumber(char *express, long len);

char *deparseNumber_cstruct(GenVar a);


GenVar combineErrors(GenVar a, GenVar b);
int combineDataDescr(GenVar a, GenVar b);


unsigned long long convGSLErr(int gsl_errno);


GenVar safeIntAdd(GVAR_INT a, GVAR_INT b);
GenVar safeIntMul(GVAR_INT a, GVAR_INT b);

GenVar safeRealAdd(double a, double b);
GenVar safeRealMul(double a, double b);


#endif
