#include "stdafx.h"
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


#include <math.h>
#include <malloc.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/err/gsl_errno.h>

#include "gvars.h"


#define MAX_NUM_LENGTH  512 /* numbers must print < this length in chars */


/*

Data simplification function.  Will return argument converted to its
"simplest" form, so:

- if DataType = DATA_IS_REAL and RealValue is 0, convert to DATA_IS_ZERO.
- if DataType = DATA_IS_REAL and RealValue is 1, convert to DATA_IS_ONE.
- if DataType = DATA_IS_REAL and RealValue is an integer, convert to data
  of type DATA_IS_INTEGER.
- if DataType = DATA_IS_INTEGER and IntegerValue is 0, conv to DATA_IS_ZERO.
- if DataType = DATA_IS_INTEGER and IntegerValue is 1, conv to DATA_IS_ONE.

Finally, integers outside the range [GVAR_INT_MIN+1,GVAR_INT_MAX-1] are
converted to reals, and reals outside the range (-GSL_DBL_MAX,GSL_DBL_MAX)
will be converted to "infinity" (positive or negative) as is appropriate.

*/

GenVar simplifyGenVar(GenVar a);



/*
Constructors for various data types.
*/

GenVar ResultZero(void)
{
    GenVar result;

    result.DataType = DATA_IS_ZERO;

    result.IntegerValue = 0;
    result.RealValue    = 0;
    result.ErrorCode    = 0;

    return result;
}

GenVar ResultOne(void)
{
    GenVar result;

    result.DataType = DATA_IS_ONE;

    result.IntegerValue = 0;
    result.RealValue    = 0;
    result.ErrorCode    = 0;

    return result;
}

GenVar ResultInteger(GVAR_INT what)
{
    GenVar result;

    result.DataType = DATA_IS_INTEGER;

    result.IntegerValue = what;
    result.RealValue    = 0;
    result.ErrorCode    = 0;

    return simplifyGenVar(result);
}

GenVar ResultReal(double what)
{
    GenVar result;

    result.DataType = DATA_IS_REAL;

    result.IntegerValue = 0;
    result.RealValue    = what;
    result.ErrorCode    = 0;

    return simplifyGenVar(result);
}

GenVar ResultPosInfnty(void)
{
    GenVar result;

    result.DataType = DATA_IS_POS_INFTY;

    result.IntegerValue = 0;
    result.RealValue    = 0;
    result.ErrorCode    = 0;

    return result;
}

GenVar ResultNegInfnty(void)
{
    GenVar result;

    result.DataType = DATA_IS_NEG_INFTY;

    result.IntegerValue = 0;
    result.RealValue    = 0;
    result.ErrorCode    = 0;

    return result;
}

GenVar ResultFinIndet(void)
{
    GenVar result;

    result.DataType = DATA_IS_FIN_INDET;

    result.IntegerValue = 0;
    result.RealValue    = 0;
    result.ErrorCode    = 0;

    return result;
}

GenVar ResultInfIndet(void)
{
    GenVar result;

    result.DataType = DATA_IS_INF_INDET;

    result.IntegerValue = 0;
    result.RealValue    = 0;
    result.ErrorCode    = 0;

    return result;
}

GenVar ResultError(unsigned long long errorCode)
{
    GenVar result;

    result.DataType = DATA_IS_ERROR;

    result.IntegerValue = 0;
    result.RealValue    = 0;
    result.ErrorCode    = errorCode;

    return result;
}

GenVar ResultUncalced(void)
{
    GenVar result;

    result.DataType = DATA_IS_UNCALCED;

    result.IntegerValue = 0;
    result.RealValue    = 0;
    result.ErrorCode    = 0;

    return result;
}




/*
Data conversion functions.
*/

double getReal(GenVar a)
{
    double result;

    switch ( a.DataType )
    {
        case DATA_IS_ZERO:
        {
            result = 0;

            break;
        }

        case DATA_IS_ONE:
        {
            result = 1;

            break;
        }

        case DATA_IS_INTEGER:
        {
            result = (double) (a.IntegerValue);

            break;
        }

        case DATA_IS_REAL:
        {
            result = a.RealValue;

            break;
        }

        case DATA_IS_POS_INFTY:
        {
            result = GSL_POSINF;

            break;
        }

        case DATA_IS_NEG_INFTY:
        {
            result = GSL_NEGINF;

            break;
        }

        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        case DATA_IS_UNCALCED:
        case DATA_IS_ERROR:
        default:
        {
            result = GSL_NAN;

            break;
        }
    }

    return result;
}

GVAR_INT getInteger(GenVar a)
{
    GVAR_INT result;

    switch ( a.DataType )
    {
        case DATA_IS_ZERO:
        {
            result = 0;

            break;
        }

        case DATA_IS_ONE:
        {
            result = 1;

            break;
        }

        case DATA_IS_INTEGER:
        {
            result = (a.IntegerValue);

            break;
        }

        case DATA_IS_REAL:
        {
            if ( (a.RealValue) > GVAR_INT_MAX-1 )
            {
                result = GVAR_INT_MAX-1;
            }

            else
            {
                if ( (a.RealValue) < GVAR_INT_MIN+1 )
                {
                    result = GVAR_INT_MIN+1;
                }

                else
                {
                    result = (GVAR_INT) (a.RealValue);
                }
            }

            break;
        }

        case DATA_IS_POS_INFTY:
        {
            result = GVAR_INT_MAX-1;

            break;
        }

        case DATA_IS_NEG_INFTY:
        {
            result = GVAR_INT_MIN+1;

            break;
        }

        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        case DATA_IS_UNCALCED:
        case DATA_IS_ERROR:
        default:
        {
            result = GVAR_INT_MAX;

            /*
               FIXME: this is a hack.  It is needed because there is no
                      NaN or inf type for GVAR_INTs in c.
            */


            break;
        }
    }

    return result;
}





/*
Data testing functions (all return 0 by default):

Programming notes: remember that data is simplified.  So, if a is of type
real then we can assume it isn't an integer, and likewise if a is an integer
we can assume that it isn't 0 or 1.
*/

int isReal(GenVar a)
{
    int result;

    switch ( a.DataType )
    {
        case DATA_IS_ZERO:
        case DATA_IS_ONE:
        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        case DATA_IS_POS_INFTY:
        case DATA_IS_NEG_INFTY:
        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        {
            result = 1;

            break;
        }

        case DATA_IS_UNCALCED:
        case DATA_IS_ERROR:
        default:
        {
            result = 0;

            break;
        }
    }

    return result;
}

int isInteger(GenVar a)
{
    int result;

    switch ( a.DataType )
    {
        case DATA_IS_ZERO:
        case DATA_IS_ONE:
        case DATA_IS_INTEGER:
        {
            result = 1;

            break;
        }

        case DATA_IS_REAL:
        {
            result = 0;

            break;
        }

        case DATA_IS_POS_INFTY:
        case DATA_IS_NEG_INFTY:
        {
            result = 1;

            break;
        }

        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = 0;

            break;
        }
    }

    return result;
}

int isFiniteReal(GenVar a)
{
    int result;

    switch ( a.DataType )
    {
        case DATA_IS_ZERO:
        case DATA_IS_ONE:
        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        {
            result = 1;

            break;
        }

        case DATA_IS_POS_INFTY:
        case DATA_IS_NEG_INFTY:
        {
            result = 0;

            break;
        }

        case DATA_IS_FIN_INDET:
        {
            result = 1;

            break;
        }

        case DATA_IS_INF_INDET:
        case DATA_IS_UNCALCED:
        case DATA_IS_ERROR:
        default:
        {
            result = 0;

            break;
        }
    }

    return result;
}

int isFiniteInteger(GenVar a)
{
    int result;

    switch ( a.DataType )
    {
        case DATA_IS_ZERO:
        case DATA_IS_ONE:
        case DATA_IS_INTEGER:
        {
            result = 1;

            break;
        }

        case DATA_IS_REAL:
        case DATA_IS_POS_INFTY:
        case DATA_IS_NEG_INFTY:
        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        case DATA_IS_UNCALCED:
        case DATA_IS_ERROR:
        default:
        {
            result = 0;

            break;
        }
    }

    return result;
}

int isError(GenVar a)
{
    int result;

    switch ( a.DataType )
    {
        case DATA_IS_ZERO:
        case DATA_IS_ONE:
        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        case DATA_IS_POS_INFTY:
        case DATA_IS_NEG_INFTY:
        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        {
            result = 0;

            break;
        }

        case DATA_IS_ERROR:
        {
            result = 1;

            break;
        }

        case DATA_IS_UNCALCED:
        default:
        {
            result = 0;

            break;
        }
    }

    return result;
}

int isUncalced(GenVar a)
{
    int result;

    switch ( a.DataType )
    {
        case DATA_IS_UNCALCED:
        {
            result = 1;

            break;
        }

        case DATA_IS_INTEGER:
        case DATA_IS_ZERO:
        case DATA_IS_ONE:
        case DATA_IS_POS_INFTY:
        case DATA_IS_NEG_INFTY:
        case DATA_IS_REAL:
        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        case DATA_IS_ERROR:
        default:
        {
            result = 0;

            break;
        }
    }

    return result;
}

int isaequalb(GenVar a, GenVar b)
{
    int result;

    switch ( combineDataDescr(a,b) )
    {
        case DATA_ZERO_ZERO:
        case DATA_ONE_ONE:
        case DATA_POS_INFTY_POS_INFTY:
        case DATA_NEG_INFTY_NEG_INFTY:
        {
            result = 1;

            break;
        }

        case DATA_INTEGER_INTEGER:
        {
            result = 0;

            if ( a.IntegerValue == b.IntegerValue )
            {
                result = 1;
            }

            break;
        }

        case DATA_REAL_REAL:
        {
            result = 0;

            if ( a.RealValue == b.RealValue )
            {
                result = 1;
            }

            break;
        }

        default:
        {
            result = 0;

            break;
        }
    }

    return result;
}

int isalessb(GenVar a, GenVar b)
{
    int result;

    switch ( combineDataDescr(a,b) )
    {
        case DATA_ZERO_ONE:
        case DATA_ZERO_POS_INFTY:
        case DATA_NEG_INFTY_ZERO:
        case DATA_ONE_POS_INFTY:
        case DATA_NEG_INFTY_ONE:
        case DATA_INTEGER_POS_INFTY:
        case DATA_NEG_INFTY_INTEGER:
        case DATA_REAL_POS_INFTY:
        case DATA_NEG_INFTY_REAL:
        case DATA_NEG_INFTY_POS_INFTY:
        {
            result = 1;

            break;
        }

        case DATA_INTEGER_INTEGER:
        case DATA_INTEGER_ZERO:
        case DATA_INTEGER_ONE:
        case DATA_ONE_INTEGER:
        case DATA_ZERO_INTEGER:
        {
            result = 0;

            if ( getInteger(a) < getInteger(b) )
            {
                result = 1;
            }

            break;
        }

        case DATA_REAL_REAL:
        case DATA_REAL_ZERO:
        case DATA_ZERO_REAL:
        case DATA_REAL_ONE:
        case DATA_ONE_REAL:
        case DATA_REAL_INTEGER:
        case DATA_INTEGER_REAL:
        {
            result = 0;

            if ( getReal(a) < getReal(b) )
            {
                result = 1;
            }

            break;
        }

        default:
        {
            result = 0;

            break;
        }
    }

    return result;
}

int isagreaterb(GenVar a, GenVar b)
{
    return isalessb(b,a);
}

int isalessequalb(GenVar a, GenVar b)
{
    return isalessb(a,b) | isaequalb(a,b);
}

int isagreaterequalb(GenVar a, GenVar b)
{
    return isagreaterequalb(b,a);
}






/*
Basic arithmetic operations.
*/

GenVar addGenVar(GenVar a, GenVar b)
{
    GenVar result;

    switch ( combineDataDescr(a,b) )
    {
        case DATA_ONE_ZERO:
        case DATA_INTEGER_ZERO:
        case DATA_REAL_ZERO:
        case DATA_POS_INFTY_ZERO:
        case DATA_NEG_INFTY_ZERO:
        case DATA_FIN_INDET_ZERO:
        case DATA_INF_INDET_ZERO:
        case DATA_ERROR_ZERO:
        case DATA_ERROR_ONE:
        case DATA_ERROR_INTEGER:
        case DATA_ERROR_REAL:
        case DATA_ERROR_POS_INFTY:
        case DATA_ERROR_NEG_INFTY:
        case DATA_ERROR_FIN_INDET:
        case DATA_ERROR_INF_INDET:
        case DATA_ERROR_UNCALCED:
        case DATA_UNCALCED_ZERO:
        case DATA_UNCALCED_ONE:
        case DATA_UNCALCED_INTEGER:
        case DATA_UNCALCED_REAL:
        case DATA_UNCALCED_POS_INFTY:
        case DATA_UNCALCED_NEG_INFTY:
        case DATA_UNCALCED_FIN_INDET:
        case DATA_UNCALCED_INF_INDET:
        case DATA_UNCALCED_UNCALCED:
        case DATA_POS_INFTY_POS_INFTY:
        case DATA_POS_INFTY_ONE:
        case DATA_POS_INFTY_INTEGER:
        case DATA_POS_INFTY_REAL:
        case DATA_POS_INFTY_FIN_INDET:
        case DATA_NEG_INFTY_NEG_INFTY:
        case DATA_NEG_INFTY_ONE:
        case DATA_NEG_INFTY_INTEGER:
        case DATA_NEG_INFTY_REAL:
        case DATA_NEG_INFTY_FIN_INDET:
        case DATA_FIN_INDET_FIN_INDET:
        case DATA_FIN_INDET_ONE:
        case DATA_FIN_INDET_INTEGER:
        case DATA_FIN_INDET_REAL:
        case DATA_INF_INDET_INF_INDET:
        case DATA_INF_INDET_ONE:
        case DATA_INF_INDET_INTEGER:
        case DATA_INF_INDET_REAL:
        case DATA_INF_INDET_POS_INFTY:
        case DATA_INF_INDET_NEG_INFTY:
        case DATA_INF_INDET_FIN_INDET:
        {
            result = a;

            break;
        }

        case DATA_ZERO_ZERO:
        case DATA_ZERO_ONE:
        case DATA_ZERO_INTEGER:
        case DATA_ZERO_REAL:
        case DATA_ZERO_POS_INFTY:
        case DATA_ZERO_NEG_INFTY:
        case DATA_ZERO_FIN_INDET:
        case DATA_ZERO_INF_INDET:
        case DATA_ZERO_ERROR:
        case DATA_ONE_ERROR:
        case DATA_INTEGER_ERROR:
        case DATA_REAL_ERROR:
        case DATA_POS_INFTY_ERROR:
        case DATA_NEG_INFTY_ERROR:
        case DATA_FIN_INDET_ERROR:
        case DATA_INF_INDET_ERROR:
        case DATA_UNCALCED_ERROR:
        case DATA_ZERO_UNCALCED:
        case DATA_ONE_UNCALCED:
        case DATA_INTEGER_UNCALCED:
        case DATA_REAL_UNCALCED:
        case DATA_POS_INFTY_UNCALCED:
        case DATA_NEG_INFTY_UNCALCED:
        case DATA_FIN_INDET_UNCALCED:
        case DATA_INF_INDET_UNCALCED:
        case DATA_ONE_POS_INFTY:
        case DATA_INTEGER_POS_INFTY:
        case DATA_REAL_POS_INFTY:
        case DATA_FIN_INDET_POS_INFTY:
        case DATA_ONE_NEG_INFTY:
        case DATA_INTEGER_NEG_INFTY:
        case DATA_REAL_NEG_INFTY:
        case DATA_FIN_INDET_NEG_INFTY:
        case DATA_ONE_FIN_INDET:
        case DATA_INTEGER_FIN_INDET:
        case DATA_REAL_FIN_INDET:
        case DATA_ONE_INF_INDET:
        case DATA_INTEGER_INF_INDET:
        case DATA_REAL_INF_INDET:
        case DATA_POS_INFTY_INF_INDET:
        case DATA_NEG_INFTY_INF_INDET:
        case DATA_FIN_INDET_INF_INDET:
        {
            result = b;

            break;
        }

        case DATA_ONE_ONE:
        case DATA_INTEGER_INTEGER:
        case DATA_INTEGER_ONE:
        case DATA_ONE_INTEGER:
        {
            result = safeIntAdd(getInteger(a),getInteger(b));

            break;
        }

        case DATA_REAL_REAL:
        case DATA_REAL_ONE:
        case DATA_ONE_REAL:
        case DATA_REAL_INTEGER:
        case DATA_INTEGER_REAL:
        {
            result = safeRealAdd(getReal(a),getReal(b));

            break;
        }

        case DATA_NEG_INFTY_POS_INFTY:
        case DATA_POS_INFTY_NEG_INFTY:
        {
            result = ResultInfIndet();

            break;
        }

        case DATA_ERROR_ERROR:
        {
            result = combineErrors(a,b);

            break;
        }

        default:
        {
            result = ResultError(GVAR_ETYPEERR);

            break;
        }
    }

    return simplifyGenVar(result);
}

GenVar subGenVar(GenVar a, GenVar b)
{
    return addGenVar(a,negGenVar(b));
}


GenVar mulGenVar(GenVar a, GenVar b)
{
    GenVar result;

    switch ( combineDataDescr(a,b) )
    {
        case DATA_ZERO_ZERO:
        case DATA_ZERO_ONE:
        case DATA_ZERO_INTEGER:
        case DATA_ZERO_REAL:
        case DATA_ZERO_UNCALCED:
        case DATA_ZERO_FIN_INDET:
        case DATA_ONE_ONE:
        case DATA_INTEGER_ONE:
        case DATA_REAL_ONE:
        case DATA_POS_INFTY_ONE:
        case DATA_NEG_INFTY_ONE:
        case DATA_FIN_INDET_ONE:
        case DATA_INF_INDET_ONE:
        case DATA_FIN_INDET_FIN_INDET:
        case DATA_INF_INDET_INF_INDET:
        case DATA_INF_INDET_FIN_INDET:
        case DATA_FIN_INDET_INTEGER:
        case DATA_INF_INDET_INTEGER:
        case DATA_FIN_INDET_REAL:
        case DATA_INF_INDET_REAL:
        case DATA_ERROR_ZERO:
        case DATA_ERROR_ONE:
        case DATA_ERROR_INTEGER:
        case DATA_ERROR_REAL:
        case DATA_ERROR_POS_INFTY:
        case DATA_ERROR_NEG_INFTY:
        case DATA_ERROR_FIN_INDET:
        case DATA_ERROR_INF_INDET:
        case DATA_ERROR_UNCALCED:
        case DATA_UNCALCED_ONE:
        case DATA_UNCALCED_INTEGER:
        case DATA_UNCALCED_REAL:
        case DATA_UNCALCED_POS_INFTY:
        case DATA_UNCALCED_NEG_INFTY:
        case DATA_UNCALCED_FIN_INDET:
        case DATA_UNCALCED_INF_INDET:
        case DATA_UNCALCED_UNCALCED:
        {
            result = a;

            break;
        }

        case DATA_ONE_ZERO:
        case DATA_INTEGER_ZERO:
        case DATA_REAL_ZERO:
        case DATA_FIN_INDET_ZERO:
        case DATA_ONE_INTEGER:
        case DATA_ONE_REAL:
        case DATA_ONE_POS_INFTY:
        case DATA_ONE_NEG_INFTY:
        case DATA_ONE_FIN_INDET:
        case DATA_ONE_INF_INDET:
        case DATA_INTEGER_FIN_INDET:
        case DATA_INTEGER_INF_INDET:
        case DATA_REAL_FIN_INDET:
        case DATA_REAL_INF_INDET:
        case DATA_ZERO_ERROR:
        case DATA_ONE_ERROR:
        case DATA_INTEGER_ERROR:
        case DATA_REAL_ERROR:
        case DATA_FIN_INDET_INF_INDET:
        case DATA_POS_INFTY_ERROR:
        case DATA_NEG_INFTY_ERROR:
        case DATA_FIN_INDET_ERROR:
        case DATA_INF_INDET_ERROR:
        case DATA_ONE_UNCALCED:
        case DATA_INTEGER_UNCALCED:
        case DATA_REAL_UNCALCED:
        case DATA_POS_INFTY_UNCALCED:
        case DATA_NEG_INFTY_UNCALCED:
        case DATA_FIN_INDET_UNCALCED:
        case DATA_INF_INDET_UNCALCED:
        case DATA_UNCALCED_ERROR:
        case DATA_UNCALCED_ZERO:
        {
            result = b;

            break;
        }

        case DATA_INTEGER_INTEGER:
        {
            result = safeIntMul(a.IntegerValue,b.IntegerValue);

            break;
        }

        case DATA_REAL_REAL:
        case DATA_REAL_INTEGER:
        case DATA_INTEGER_REAL:
        {
            result = safeRealMul(getReal(a),getReal(b));

            break;
        }

        case DATA_POS_INFTY_POS_INFTY:
        case DATA_NEG_INFTY_NEG_INFTY:
        {
            result = ResultPosInfnty();

            break;
        }

        case DATA_NEG_INFTY_POS_INFTY:
        case DATA_POS_INFTY_NEG_INFTY:
        {
            result = ResultNegInfnty();

            break;
        }

        case DATA_POS_INFTY_ZERO:
        case DATA_ZERO_POS_INFTY:
        case DATA_NEG_INFTY_ZERO:
        case DATA_ZERO_NEG_INFTY:
        case DATA_INF_INDET_ZERO:
        case DATA_ZERO_INF_INDET:
        case DATA_FIN_INDET_POS_INFTY:
        case DATA_POS_INFTY_FIN_INDET:
        case DATA_INF_INDET_POS_INFTY:
        case DATA_POS_INFTY_INF_INDET:
        case DATA_FIN_INDET_NEG_INFTY:
        case DATA_NEG_INFTY_FIN_INDET:
        case DATA_INF_INDET_NEG_INFTY:
        case DATA_NEG_INFTY_INF_INDET:
        {
            result = ResultInfIndet();

            break;
        }

        case DATA_POS_INFTY_INTEGER:
        {
            if ( getInteger(b) > 0 )
            {
                result = ResultPosInfnty();
            }

            else
            {
                result = ResultNegInfnty();
            }

            break;
        }

        case DATA_POS_INFTY_REAL:
        {
            if ( getReal(b) > 0 )
            {
                result = ResultPosInfnty();
            }

            else
            {
                result = ResultNegInfnty();
            }

            break;
        }

        case DATA_INTEGER_POS_INFTY:
        {
            if ( getInteger(a) > 0 )
            {
                result = ResultPosInfnty();
            }

            else
            {
                result = ResultNegInfnty();
            }

            break;
        }

        case DATA_REAL_POS_INFTY:
        {
            if ( getReal(a) > 0 )
            {
                result = ResultPosInfnty();
            }

            else
            {
                result = ResultNegInfnty();
            }

            break;
        }

        case DATA_NEG_INFTY_INTEGER:
        {
            if ( getInteger(b) < 0 )
            {
                result = ResultPosInfnty();
            }

            else
            {
                result = ResultNegInfnty();
            }

            break;
        }

        case DATA_NEG_INFTY_REAL:
        {
            if ( getReal(b) < 0 )
            {
                result = ResultPosInfnty();
            }

            else
            {
                result = ResultNegInfnty();
            }

            break;
        }

        case DATA_INTEGER_NEG_INFTY:
        {
            if ( getInteger(a) < 0 )
            {
                result = ResultPosInfnty();
            }

            else
            {
                result = ResultNegInfnty();
            }

            break;
        }

        case DATA_REAL_NEG_INFTY:
        {
            if ( getReal(a) < 0 )
            {
                result = ResultPosInfnty();
            }

            else
            {
                result = ResultNegInfnty();
            }

            break;
        }

        case DATA_ERROR_ERROR:
        {
            result = combineErrors(a,b);

            break;
        }

        default:
        {
            result = ResultError(GVAR_ETYPEERR);

            break;
        }
    }

    return simplifyGenVar(result);
}

GenVar divGenVar(GenVar a, GenVar b)
{
    GenVar result;

    switch ( combineDataDescr(a,b) )
    {
        case DATA_ONE_ONE:
        case DATA_ZERO_ONE:
        case DATA_ZERO_INTEGER:
        case DATA_ZERO_REAL:
        case DATA_POS_INFTY_ZERO:
        case DATA_ZERO_POS_INFTY:
        case DATA_NEG_INFTY_ZERO:
        case DATA_ZERO_NEG_INFTY:
        case DATA_INTEGER_ONE:
        case DATA_REAL_ONE:
        case DATA_POS_INFTY_ONE:
        case DATA_NEG_INFTY_ONE:
        case DATA_FIN_INDET_ONE:
        case DATA_FIN_INDET_INTEGER:
        case DATA_FIN_INDET_REAL:
        case DATA_INF_INDET_ONE:
        case DATA_INF_INDET_INTEGER:
        case DATA_INF_INDET_REAL:
        case DATA_ERROR_ZERO:
        case DATA_ERROR_ONE:
        case DATA_ERROR_INTEGER:
        case DATA_ERROR_REAL:
        case DATA_ERROR_POS_INFTY:
        case DATA_ERROR_NEG_INFTY:
        case DATA_ERROR_FIN_INDET:
        case DATA_ERROR_INF_INDET:
        case DATA_UNCALCED_ZERO:
        case DATA_UNCALCED_ONE:
        case DATA_UNCALCED_INTEGER:
        case DATA_UNCALCED_REAL:
        case DATA_UNCALCED_POS_INFTY:
        case DATA_UNCALCED_NEG_INFTY:
        case DATA_UNCALCED_FIN_INDET:
        case DATA_UNCALCED_INF_INDET:
        case DATA_ERROR_UNCALCED:
        case DATA_UNCALCED_UNCALCED:
        case DATA_ZERO_UNCALCED:
        case DATA_INF_INDET_UNCALCED:
        case DATA_INF_INDET_INF_INDET:
        case DATA_INF_INDET_ZERO:
        case DATA_INF_INDET_FIN_INDET:
        case DATA_INF_INDET_POS_INFTY:
        case DATA_INF_INDET_NEG_INFTY:
        {
            result = a;

            break;
        }

        case DATA_ZERO_ERROR:
        case DATA_ONE_ERROR:
        case DATA_INTEGER_ERROR:
        case DATA_REAL_ERROR:
        case DATA_POS_INFTY_ERROR:
        case DATA_NEG_INFTY_ERROR:
        case DATA_FIN_INDET_ERROR:
        case DATA_INF_INDET_ERROR:
        case DATA_ONE_UNCALCED:
        case DATA_INTEGER_UNCALCED:
        case DATA_REAL_UNCALCED:
        case DATA_POS_INFTY_UNCALCED:
        case DATA_NEG_INFTY_UNCALCED:
        case DATA_FIN_INDET_UNCALCED:
        case DATA_UNCALCED_ERROR:
        case DATA_ZERO_INF_INDET:
        case DATA_ONE_INF_INDET:
        case DATA_INTEGER_INF_INDET:
        case DATA_REAL_INF_INDET:
        case DATA_POS_INFTY_INF_INDET:
        case DATA_NEG_INFTY_INF_INDET:
        case DATA_FIN_INDET_INF_INDET:
        {
            result = b;

            break;
        }

        case DATA_INTEGER_INTEGER:
        case DATA_ONE_INTEGER:
        {
            if ( (getInteger(a)%getInteger(b)) == 0 )
            {
                result = ResultInteger(getInteger(a)/getInteger(b));
            }

            else
            {
                result = ResultReal(getReal(a)/getReal(b));
            }

            break;
        }

        case DATA_REAL_REAL:
        case DATA_REAL_INTEGER:
        case DATA_INTEGER_REAL:
        case DATA_ONE_REAL:
        {
            result = ResultReal(getReal(a)/getReal(b));

            break;
        }

        case DATA_POS_INFTY_POS_INFTY:
        case DATA_NEG_INFTY_NEG_INFTY:
        case DATA_NEG_INFTY_POS_INFTY:
        case DATA_POS_INFTY_NEG_INFTY:
        case DATA_FIN_INDET_FIN_INDET:
        case DATA_POS_INFTY_FIN_INDET:
        case DATA_NEG_INFTY_FIN_INDET:
        case DATA_FIN_INDET_ZERO:
        case DATA_ZERO_FIN_INDET:
        case DATA_ONE_FIN_INDET:
        case DATA_INTEGER_FIN_INDET:
        case DATA_REAL_FIN_INDET:
        case DATA_ZERO_ZERO:
        {
            result = ResultInfIndet();

            break;
        }

        case DATA_ONE_ZERO:
        {
            result = ResultPosInfnty();

            break;
        }

        case DATA_INTEGER_ZERO:
        {
            if ( getInteger(a) > 0 )
            {
                result = ResultPosInfnty();
            }

            else
            {
                result = ResultNegInfnty();
            }

            break;
        }

        case DATA_REAL_ZERO:
        {
            if ( getReal(a) > 0 )
            {
                result = ResultPosInfnty();
            }

            else
            {
                result = ResultNegInfnty();
            }

            break;
        }

        case DATA_ONE_POS_INFTY:
        case DATA_ONE_NEG_INFTY:
        case DATA_INTEGER_POS_INFTY:
        case DATA_INTEGER_NEG_INFTY:
        case DATA_REAL_POS_INFTY:
        case DATA_REAL_NEG_INFTY:
        case DATA_FIN_INDET_POS_INFTY:
        case DATA_FIN_INDET_NEG_INFTY:
        {
            result = ResultZero();

            break;
        }

        case DATA_POS_INFTY_INTEGER:
        {
            if ( getInteger(b) > 0 )
            {
                result = ResultPosInfnty();
            }

            else
            {
                result = ResultNegInfnty();
            }

            break;
        }

        case DATA_POS_INFTY_REAL:
        {
            if ( getReal(b) > 0 )
            {
                result = ResultPosInfnty();
            }

            else
            {
                result = ResultNegInfnty();
            }

            break;
        }

        case DATA_NEG_INFTY_INTEGER:
        {
            if ( getInteger(b) < 0 )
            {
                result = ResultPosInfnty();
            }

            else
            {
                result = ResultNegInfnty();
            }

            break;
        }

        case DATA_NEG_INFTY_REAL:
        {
            if ( getReal(b) < 0 )
            {
                result = ResultPosInfnty();
            }

            else
            {
                result = ResultNegInfnty();
            }

            break;
        }

        case DATA_ERROR_ERROR:
        {
            result = combineErrors(a,b);

            break;
        }

        default:
        {
            result = ResultError(GVAR_ETYPEERR);

            break;
        }
    }

    return simplifyGenVar(result);
}

GenVar idivGenVar(GenVar a, GenVar b)
{
    GenVar result;

    switch ( combineDataDescr(a,b) )
    {
        case DATA_ONE_ONE:
        case DATA_POS_INFTY_ZERO:
        case DATA_NEG_INFTY_ZERO:
        case DATA_ZERO_ONE:
        case DATA_ZERO_INTEGER:
        case DATA_ZERO_POS_INFTY:
        case DATA_ZERO_NEG_INFTY:
        case DATA_INTEGER_ONE:
        case DATA_POS_INFTY_ONE:
        case DATA_NEG_INFTY_ONE:
        case DATA_ERROR_ZERO:
        case DATA_ERROR_ONE:
        case DATA_ERROR_INTEGER:
        case DATA_ERROR_POS_INFTY:
        case DATA_ERROR_NEG_INFTY:
        case DATA_UNCALCED_ZERO:
        case DATA_UNCALCED_ONE:
        case DATA_UNCALCED_INTEGER:
        case DATA_UNCALCED_POS_INFTY:
        case DATA_UNCALCED_NEG_INFTY:
        case DATA_ERROR_UNCALCED:
        case DATA_UNCALCED_UNCALCED:
        {
            result = a;

            break;
        }

        case DATA_ZERO_ERROR:
        case DATA_ONE_ERROR:
        case DATA_INTEGER_ERROR:
        case DATA_POS_INFTY_ERROR:
        case DATA_NEG_INFTY_ERROR:
        case DATA_ONE_UNCALCED:
        case DATA_INTEGER_UNCALCED:
        case DATA_POS_INFTY_UNCALCED:
        case DATA_NEG_INFTY_UNCALCED:
        case DATA_UNCALCED_ERROR:
        {
            result = b;

            break;
        }

        case DATA_INTEGER_INTEGER:
        case DATA_ONE_INTEGER:
        {
            result = ResultInteger(getInteger(a)/getInteger(b));

            break;
        }

        case DATA_ONE_ZERO:
        {
            result = ResultPosInfnty();

            break;
        }

        case DATA_INTEGER_ZERO:
        {
            if ( getInteger(a) > 0 )
            {
                result = ResultPosInfnty();
            }

            else
            {
                result = ResultNegInfnty();
            }

            break;
        }

        case DATA_ONE_POS_INFTY:
        case DATA_ONE_NEG_INFTY:
        case DATA_INTEGER_POS_INFTY:
        case DATA_INTEGER_NEG_INFTY:
        {
            result = ResultZero();

            break;
        }

        case DATA_POS_INFTY_INTEGER:
        {
            if ( getInteger(b) > 0 )
            {
                result = ResultPosInfnty();
            }

            else
            {
                result = ResultNegInfnty();
            }

            break;
        }

        case DATA_NEG_INFTY_INTEGER:
        {
            if ( getInteger(b) < 0 )
            {
                result = ResultPosInfnty();
            }

            else
            {
                result = ResultNegInfnty();
            }

            break;
        }

        case DATA_ERROR_ERROR:
        {
            result = combineErrors(a,b);

            break;
        }

        case DATA_ZERO_REAL:
        case DATA_REAL_ONE:
        case DATA_FIN_INDET_ONE:
        case DATA_INF_INDET_ONE:
        case DATA_FIN_INDET_INTEGER:
        case DATA_INF_INDET_INTEGER:
        case DATA_FIN_INDET_REAL:
        case DATA_INF_INDET_REAL:
        case DATA_REAL_REAL:
        case DATA_REAL_INTEGER:
        case DATA_INTEGER_REAL:
        case DATA_ONE_REAL:
        case DATA_FIN_INDET_FIN_INDET:
        case DATA_INF_INDET_INF_INDET:
        case DATA_FIN_INDET_ZERO:
        case DATA_ZERO_FIN_INDET:
        case DATA_INF_INDET_ZERO:
        case DATA_ZERO_INF_INDET:
        case DATA_ONE_FIN_INDET:
        case DATA_ONE_INF_INDET:
        case DATA_INTEGER_FIN_INDET:
        case DATA_INTEGER_INF_INDET:
        case DATA_REAL_FIN_INDET:
        case DATA_REAL_INF_INDET:
        case DATA_POS_INFTY_FIN_INDET:
        case DATA_INF_INDET_POS_INFTY:
        case DATA_POS_INFTY_INF_INDET:
        case DATA_NEG_INFTY_FIN_INDET:
        case DATA_INF_INDET_NEG_INFTY:
        case DATA_NEG_INFTY_INF_INDET:
        case DATA_INF_INDET_FIN_INDET:
        case DATA_FIN_INDET_INF_INDET:
        case DATA_REAL_ZERO:
        case DATA_REAL_POS_INFTY:
        case DATA_REAL_NEG_INFTY:
        case DATA_FIN_INDET_POS_INFTY:
        case DATA_FIN_INDET_NEG_INFTY:
        case DATA_POS_INFTY_REAL:
        case DATA_NEG_INFTY_REAL:
        case DATA_UNCALCED_REAL:
        case DATA_UNCALCED_FIN_INDET:
        case DATA_UNCALCED_INF_INDET:
        case DATA_REAL_UNCALCED:
        case DATA_FIN_INDET_UNCALCED:
        case DATA_INF_INDET_UNCALCED:
        case DATA_ZERO_ZERO:
        case DATA_ZERO_UNCALCED:
        case DATA_POS_INFTY_POS_INFTY:
        case DATA_NEG_INFTY_NEG_INFTY:
        case DATA_NEG_INFTY_POS_INFTY:
        case DATA_POS_INFTY_NEG_INFTY:
        {
            result = ResultError(GVAR_GSL_EDOM);

            break;
        }

        case DATA_ERROR_REAL:
        case DATA_ERROR_FIN_INDET:
        case DATA_ERROR_INF_INDET:
        {
            result = combineErrors(a,ResultError(GVAR_GSL_EDOM));

            break;
        }

        case DATA_REAL_ERROR:
        case DATA_FIN_INDET_ERROR:
        case DATA_INF_INDET_ERROR:
        {
            result = combineErrors(ResultError(GVAR_GSL_EDOM),b);

            break;
        }

        default:
        {
            result = ResultError(GVAR_ETYPEERR);

            break;
        }
    }

    return simplifyGenVar(result);
}

GenVar modGenVar(GenVar a, GenVar b)
{
    GenVar result;

    switch ( combineDataDescr(a,b) )
    {
        case DATA_ZERO_POS_INFTY:
        case DATA_ONE_POS_INFTY:
        case DATA_INTEGER_POS_INFTY:
        case DATA_ZERO_NEG_INFTY:
        case DATA_ONE_NEG_INFTY:
        case DATA_INTEGER_NEG_INFTY:
        case DATA_ERROR_ZERO:
        case DATA_ERROR_ONE:
        case DATA_ERROR_INTEGER:
        case DATA_ERROR_POS_INFTY:
        case DATA_ERROR_NEG_INFTY:
        case DATA_ERROR_UNCALCED:
        case DATA_UNCALCED_ZERO:
        case DATA_UNCALCED_ONE:
        case DATA_UNCALCED_INTEGER:
        case DATA_UNCALCED_POS_INFTY:
        case DATA_UNCALCED_NEG_INFTY:
        case DATA_UNCALCED_UNCALCED:
        {
            result = a;

            break;
        }

        case DATA_ZERO_ZERO:
        case DATA_ONE_ZERO:
        case DATA_INTEGER_ZERO:
        case DATA_POS_INFTY_ZERO:
        case DATA_POS_INFTY_ONE:
        case DATA_POS_INFTY_INTEGER:
        case DATA_ZERO_ERROR:
        case DATA_ONE_ERROR:
        case DATA_INTEGER_ERROR:
        case DATA_POS_INFTY_ERROR:
        case DATA_NEG_INFTY_ERROR:
        case DATA_UNCALCED_ERROR:
        case DATA_ZERO_UNCALCED:
        case DATA_ONE_UNCALCED:
        case DATA_INTEGER_UNCALCED:
        case DATA_POS_INFTY_UNCALCED:
        case DATA_NEG_INFTY_UNCALCED:
        {
            result = b;

            break;
        }

        case DATA_NEG_INFTY_ZERO:
        case DATA_NEG_INFTY_ONE:
        case DATA_NEG_INFTY_INTEGER:
        {
            result = negGenVar(b);

            break;
        }

        case DATA_ZERO_ONE:
        case DATA_ZERO_INTEGER:
        case DATA_ONE_ONE:
        case DATA_ONE_INTEGER:
        case DATA_INTEGER_ONE:
        case DATA_INTEGER_INTEGER:
        {
            result = ResultInteger(getInteger(a)%getInteger(b));

            break;
        }

        case DATA_POS_INFTY_POS_INFTY:
        case DATA_NEG_INFTY_NEG_INFTY:
        case DATA_NEG_INFTY_POS_INFTY:
        case DATA_POS_INFTY_NEG_INFTY:
        {
            result = ResultInfIndet();

            break;
        }

        case DATA_ZERO_REAL:
        case DATA_REAL_ONE:
        case DATA_FIN_INDET_ONE:
        case DATA_INF_INDET_ONE:
        case DATA_FIN_INDET_INTEGER:
        case DATA_INF_INDET_INTEGER:
        case DATA_FIN_INDET_REAL:
        case DATA_INF_INDET_REAL:
        case DATA_UNCALCED_REAL:
        case DATA_UNCALCED_FIN_INDET:
        case DATA_UNCALCED_INF_INDET:
        case DATA_REAL_UNCALCED:
        case DATA_FIN_INDET_UNCALCED:
        case DATA_INF_INDET_UNCALCED:
        case DATA_REAL_REAL:
        case DATA_REAL_INTEGER:
        case DATA_INTEGER_REAL:
        case DATA_ONE_REAL:
        case DATA_FIN_INDET_FIN_INDET:
        case DATA_INF_INDET_INF_INDET:
        case DATA_FIN_INDET_ZERO:
        case DATA_ZERO_FIN_INDET:
        case DATA_INF_INDET_ZERO:
        case DATA_ZERO_INF_INDET:
        case DATA_ONE_FIN_INDET:
        case DATA_ONE_INF_INDET:
        case DATA_INTEGER_FIN_INDET:
        case DATA_INTEGER_INF_INDET:
        case DATA_REAL_FIN_INDET:
        case DATA_REAL_INF_INDET:
        case DATA_POS_INFTY_FIN_INDET:
        case DATA_INF_INDET_POS_INFTY:
        case DATA_POS_INFTY_INF_INDET:
        case DATA_NEG_INFTY_FIN_INDET:
        case DATA_INF_INDET_NEG_INFTY:
        case DATA_NEG_INFTY_INF_INDET:
        case DATA_INF_INDET_FIN_INDET:
        case DATA_FIN_INDET_INF_INDET:
        case DATA_REAL_ZERO:
        case DATA_REAL_POS_INFTY:
        case DATA_REAL_NEG_INFTY:
        case DATA_FIN_INDET_POS_INFTY:
        case DATA_FIN_INDET_NEG_INFTY:
        case DATA_POS_INFTY_REAL:
        case DATA_NEG_INFTY_REAL:
        {
            result = ResultError(GVAR_GSL_EDOM);

            break;
        }

        case DATA_ERROR_REAL:
        case DATA_ERROR_FIN_INDET:
        case DATA_ERROR_INF_INDET:
        {
            result = combineErrors(a,ResultError(GVAR_GSL_EDOM));

            break;
        }

        case DATA_REAL_ERROR:
        case DATA_FIN_INDET_ERROR:
        case DATA_INF_INDET_ERROR:
        {
            result = combineErrors(ResultError(GVAR_GSL_EDOM),b);

            break;
        }

        case DATA_ERROR_ERROR:
        {
            result = combineErrors(a,b);

            break;
        }

        default:
        {
            result = ResultError(GVAR_ETYPEERR);

            break;
        }
    }

    return simplifyGenVar(result);
}

GenVar negGenVar(GenVar a)
{
    GenVar result;

    switch ( a.DataType )
    {
        case DATA_IS_ZERO:
        {
            result = a;

            break;
        }

        case DATA_IS_ONE:
        {
            result = ResultInteger(-1);

            break;
        }

        case DATA_IS_INTEGER:
        {
            result = ResultInteger(-(a.IntegerValue));

            break;
        }

        case DATA_IS_REAL:
        {
            result = ResultReal(-(a.RealValue));

            break;
        }

        case DATA_IS_POS_INFTY:
        {
            result = ResultNegInfnty();

            break;
        }

        case DATA_IS_NEG_INFTY:
        {
            result = ResultPosInfnty();

            break;
        }

        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        case DATA_IS_UNCALCED:
        case DATA_IS_ERROR:
        {
            result = a;

            break;
        }

        default:
        {
            result = ResultError(GVAR_ETYPEERR);

            break;
        }
    }

    return simplifyGenVar(result);
}

GenVar invGenVar(GenVar a)
{
    GenVar result;

    switch ( a.DataType )
    {
        case DATA_IS_ZERO:
        {
            result = ResultPosInfnty();

            break;
        }

        case DATA_IS_ONE:
        {
            result = a;

            break;
        }

        case DATA_IS_INTEGER:
        {
            result = ResultReal(1/getReal(a));

            break;
        }

        case DATA_IS_REAL:
        {
            result = ResultReal(1/(a.RealValue));

            break;
        }

        case DATA_IS_POS_INFTY:
        case DATA_IS_NEG_INFTY:
        {
            result = ResultZero();

            break;
        }

        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        {
            result = ResultInfIndet();

            break;
        }

        case DATA_IS_UNCALCED:
        case DATA_IS_ERROR:
        {
            result = a;

            break;
        }

        default:
        {
            result = ResultError(GVAR_ETYPEERR);

            break;
        }
    }

    return simplifyGenVar(result);
}





/*
Number parsing and deparsing.
*/

char *deparseNumber(GenVar a)
{
    char *result;

    if ( ( result = (char *) malloc((MAX_NUM_LENGTH+1)*sizeof(char)) ) == NULL )
    {
        return NULL;
    }

    switch ( a.DataType )
    {
        case DATA_IS_ZERO:
        {
            sprintf(result,"0");

            break;
        }

        case DATA_IS_ONE:
        {
            sprintf(result,"1");

            break;
        }

        case DATA_IS_INTEGER:
        {
            if ( getInteger(a) > 0 )
            {
                sprintf(result,"%Lu",((unsigned GVAR_INT) getInteger(a)));
            }

            else
            {
                sprintf(result,"-%Lu",((unsigned GVAR_INT) -getInteger(a)));
            }

            break;
        }

        case DATA_IS_REAL:
        {
            sprintf(result,"%1.15le",getReal(a));

            break;
        }

        case DATA_IS_POS_INFTY:
        {
            sprintf(result,"inf");

            break;
        }

        case DATA_IS_NEG_INFTY:
        {
            sprintf(result,"-inf");

            break;
        }

        case DATA_IS_FIN_INDET:
        {
            sprintf(result,"FiniteIndet");

            break;
        }

        case DATA_IS_INF_INDET:
        {
            sprintf(result,"InfintIndet");

            break;
        }

        case DATA_IS_UNCALCED:
        {
            sprintf(result,"uc");

            break;
        }

        case DATA_IS_ERROR:
        {
            sprintf(result,"t%Lu",a.ErrorCode);

            break;
        }

        default:
        {
            sprintf(result,"t%Lu",((unsigned GVAR_INT) GVAR_ESYNTAX));

            break;
        }
    }

    return result;
}

char *deparseNumber_cstruct(GenVar a)
{
    char *result;

    if ( ( result = (char *) malloc((MAX_NUM_LENGTH+1)*sizeof(char)) ) == NULL )
    {
        return NULL;
    }

    if ( getInteger(a) >= 0 )
    {
        sprintf(result,"{ %d , %Lu , %1.15le , %Lu }",a.DataType,((unsigned GVAR_INT) getInteger(a)),a.RealValue,a.ErrorCode);
    }

    else
    {
        sprintf(result,"{ %d , -%Lu , %1.15le , %Lu }",a.DataType,((unsigned GVAR_INT) -getInteger(a)),a.RealValue,a.ErrorCode);
    }


    return result;
}

GenVar parseNumber(char *express, long len)
{
    GenVar result;
    long i;
    double divisor;
    GenVar exponent;
    unsigned long long whaterr;

    if ( len > 0 )
    {
        switch ( express[0] )
        {
            case '-':
            {
                result = negGenVar(parseNumber(express+1,len-1));

                break;
            }

            case 'i':
            {
                if ( ( len == strlen("inf") ) && ( strcmp(express,"inf") == 0 ) )
                {
                    result = ResultPosInfnty();
                }

                else
                {
                    result = ResultError(GVAR_ESYNTAX);
                }

                break;
            }

            case 'F':
            {
                if ( ( len == strlen("FiniteIndet") ) && ( strcmp(express,"FiniteIndet") == 0 ) )
                {
                    result = ResultFinIndet();
                }

                else
                {
                    result = ResultError(GVAR_ESYNTAX);
                }

                break;
            }

            case 'I':
            {
                if ( ( len == strlen("InfintIndet") ) && ( strcmp(express,"InfintIndet") == 0 ) )
                {
                    result = ResultInfIndet();
                }

                else
                {
                    result = ResultError(GVAR_ESYNTAX);
                }

                break;
            }

            case '0': case '1': case '2': case '3': case '4':
            case '5': case '6': case '7': case '8': case '9':
            case '.':
            {
                result = ResultZero();

                i = 0;

                while ( i < len )
                {
                    switch ( express[i] )
                    {
                        case '0': { result = addGenVar(mulGenVar(result,ResultInteger(10)),ResultInteger(0)); break; }
                        case '1': { result = addGenVar(mulGenVar(result,ResultInteger(10)),ResultInteger(1)); break; }
                        case '2': { result = addGenVar(mulGenVar(result,ResultInteger(10)),ResultInteger(2)); break; }
                        case '3': { result = addGenVar(mulGenVar(result,ResultInteger(10)),ResultInteger(3)); break; }
                        case '4': { result = addGenVar(mulGenVar(result,ResultInteger(10)),ResultInteger(4)); break; }
                        case '5': { result = addGenVar(mulGenVar(result,ResultInteger(10)),ResultInteger(5)); break; }
                        case '6': { result = addGenVar(mulGenVar(result,ResultInteger(10)),ResultInteger(6)); break; }
                        case '7': { result = addGenVar(mulGenVar(result,ResultInteger(10)),ResultInteger(7)); break; }
                        case '8': { result = addGenVar(mulGenVar(result,ResultInteger(10)),ResultInteger(8)); break; }
                        case '9': { result = addGenVar(mulGenVar(result,ResultInteger(10)),ResultInteger(9)); break; }

                        case '.':
                        case 'e':
                        case 'E':
                        {
                            goto quick_entry_real;

                            break;
                        }

                        default:
                        {
                            result = ResultError(GVAR_ESYNTAX);

                            goto quick_exit_point;

                            break;
                        }
                    }

                    i++;
                }

                goto quick_exit_point;

                quick_entry_real:

                if ( express[i] == '.' )
                {
                    i++;

                    divisor = 10;

                    while ( i < len )
                    {
                        switch ( express[i] )
                        {
                            case '0': { result = addGenVar(result,ResultReal(0/divisor)); break; }
                            case '1': { result = addGenVar(result,ResultReal(1/divisor)); break; }
                            case '2': { result = addGenVar(result,ResultReal(2/divisor)); break; }
                            case '3': { result = addGenVar(result,ResultReal(3/divisor)); break; }
                            case '4': { result = addGenVar(result,ResultReal(4/divisor)); break; }
                            case '5': { result = addGenVar(result,ResultReal(5/divisor)); break; }
                            case '6': { result = addGenVar(result,ResultReal(6/divisor)); break; }
                            case '7': { result = addGenVar(result,ResultReal(7/divisor)); break; }
                            case '8': { result = addGenVar(result,ResultReal(8/divisor)); break; }
                            case '9': { result = addGenVar(result,ResultReal(9/divisor)); break; }

                            case 'e':
                            case 'E':
                            {
                                goto quick_entry_exponent;

                                break;
                            }

                            default:
                            {
                                result = ResultError(GVAR_ESYNTAX);

                                goto quick_exit_point;

                                break;
                            }
                        }

                        i++;
                        divisor *= 10;
                    }
                }

                else
                {
                    quick_entry_exponent:

                    i++;

                    exponent = parseNumber(express+i,len-i);

                    switch ( exponent.DataType )
                    {
                        case DATA_IS_ZERO:
                        case DATA_IS_ONE:
                        case DATA_IS_INTEGER:
                        case DATA_IS_REAL:
                        {
                            result = mulGenVar(result,ResultReal(pow(10,getReal(exponent))));

                            break;
                        }

                        case DATA_IS_POS_INFTY:
                        case DATA_IS_NEG_INFTY:
                        case DATA_IS_FIN_INDET:
                        case DATA_IS_INF_INDET:
                        default:
                        {
                            result = ResultError(GVAR_ESYNTAX);

                            break;
                        }
                    }
                }

                quick_exit_point:

                break;
            }

            case 'u':
            {
                if ( ( len == strlen("uc") ) && ( strcmp(express,"uc") == 0 ) )
                {
                    result = ResultUncalced();
                }

                else
                {
                    result = ResultError(GVAR_ESYNTAX);
                }

                break;
            }

            case 't':
            {
                sscanf(express+1,"%Lu",&whaterr);

                result = parseNumber(express+1,len-1);

                if ( isInteger(result) )
                {
                    result = ResultError(getInteger(result));
                }

                else
                {
                    result = ResultError(GVAR_ESYNTAX);
                }

                break;
            }

            default:
            {
                result = ResultError(GVAR_ESYNTAX);

                break;
            }
        }
    }

    else
    {
        result = ResultError(GVAR_ESYNTAX);
    }

    return simplifyGenVar(result);
}






/*
Data simplification function.
*/

GenVar simplifyGenVar(GenVar a)
{
    GenVar result;

    switch ( a.DataType )
    {
        case DATA_IS_ZERO:
        case DATA_IS_ONE:
        {
            result = a;

            break;
        }

        case DATA_IS_INTEGER:
        {
            if ( a.IntegerValue == 0 )
            {
                result = ResultZero();
            }

            else
            {
                if ( a.IntegerValue == 1 )
                {
                    result = ResultOne();
                }

                else
                {
                    if ( ( a.IntegerValue >= GVAR_INT_MAX ) || 
                         ( a.IntegerValue <= GVAR_INT_MIN )     )
                    {
                        result = ResultReal((double) (a.IntegerValue));
                    }

                    else
                    {
                        result = a;
                    }
                }
            }

            break;
        }

        case DATA_IS_REAL:
        {
            if ( gsl_isnan(a.RealValue) )
            {
                result = ResultError(GVAR_ENANDET);
            }

            else
            {
                if ( a.RealValue == 0 )
                {
                    result = ResultZero();
                }

                else
                {
                    if ( a.RealValue == 1 )
                    {
                        result = ResultOne();
                    }

                    else
                    {
                        if ( ( ((double) ((GVAR_INT) (a.RealValue))) == a.RealValue ) &&
                             ( ((GVAR_INT) (a.RealValue)) < GVAR_INT_MAX ) &&
                             ( ((GVAR_INT) (a.RealValue)) > GVAR_INT_MIN )    )
                        {
                            result = simplifyGenVar(ResultInteger((GVAR_INT) (a.RealValue)));
                        }

                        else
                        {
                            switch ( gsl_isinf(a.RealValue) )
                            {
                                case +1:
                                {
                                    result = ResultPosInfnty();

                                    break;
                                }

                                case -1:
                                {
                                    result = ResultNegInfnty();

                                    break;
                                }

                                default:
                                {
                                    if ( ( a.RealValue >= GSL_DBL_MAX  ) ||
                                         ( a.RealValue <= -GSL_DBL_MAX )    )
                                    {
                                        if ( a.RealValue > 1 )
                                        {
                                            result = ResultPosInfnty();
                                        }

                                        else
                                        {
                                            result = ResultNegInfnty();
                                        }
                                    }

                                    else
                                    {
                                        result = a;
                                    }

                                    break;
                                }
                            }
                        }
                    }
                }
            }

            break;
        }

        case DATA_IS_POS_INFTY:
        case DATA_IS_NEG_INFTY:
        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        case DATA_IS_UNCALCED:
        case DATA_IS_ERROR:
        {
            result = a;

            break;
        }

        default:
        {
            result = ResultError(GVAR_ETYPEERR);

            break;
        }
    }

    return result;
}



int combineDataDescr(GenVar a, GenVar b)
{
    int result;

    switch ( a.DataType )
    {
        case DATA_IS_ZERO:
        {
            switch ( b.DataType )
            {
                case DATA_IS_ZERO:      { result = DATA_ZERO_ZERO;      break; }
                case DATA_IS_ONE:       { result = DATA_ZERO_ONE;       break; }
                case DATA_IS_INTEGER:   { result = DATA_ZERO_INTEGER;   break; }
                case DATA_IS_REAL:      { result = DATA_ZERO_REAL;      break; }
                case DATA_IS_POS_INFTY: { result = DATA_ZERO_POS_INFTY; break; }
                case DATA_IS_NEG_INFTY: { result = DATA_ZERO_NEG_INFTY; break; }
                case DATA_IS_FIN_INDET: { result = DATA_ZERO_FIN_INDET; break; }
                case DATA_IS_INF_INDET: { result = DATA_ZERO_INF_INDET; break; }
                case DATA_IS_UNCALCED:  { result = DATA_ZERO_UNCALCED;  break; }
                case DATA_IS_ERROR:     { result = DATA_ZERO_ERROR;     break; }
                default:                { result = DATA_ZERO_ERROR;     break; }
            }

            break;
        }

        case DATA_IS_ONE:
        {
            switch ( b.DataType )
            {
                case DATA_IS_ZERO:      { result = DATA_ONE_ZERO;      break; }
                case DATA_IS_ONE:       { result = DATA_ONE_ONE;       break; }
                case DATA_IS_INTEGER:   { result = DATA_ONE_INTEGER;   break; }
                case DATA_IS_REAL:      { result = DATA_ONE_REAL;      break; }
                case DATA_IS_POS_INFTY: { result = DATA_ONE_POS_INFTY; break; }
                case DATA_IS_NEG_INFTY: { result = DATA_ONE_NEG_INFTY; break; }
                case DATA_IS_FIN_INDET: { result = DATA_ONE_FIN_INDET; break; }
                case DATA_IS_INF_INDET: { result = DATA_ONE_INF_INDET; break; }
                case DATA_IS_UNCALCED:  { result = DATA_ONE_UNCALCED;  break; }
                case DATA_IS_ERROR:     { result = DATA_ONE_ERROR;     break; }
                default:                { result = DATA_ONE_ERROR;     break; }
            }

            break;
        }

        case DATA_IS_INTEGER:
        {
            switch ( b.DataType )
            {
                case DATA_IS_ZERO:      { result = DATA_INTEGER_ZERO;      break; }
                case DATA_IS_ONE:       { result = DATA_INTEGER_ONE;       break; }
                case DATA_IS_INTEGER:   { result = DATA_INTEGER_INTEGER;   break; }
                case DATA_IS_REAL:      { result = DATA_INTEGER_REAL;      break; }
                case DATA_IS_POS_INFTY: { result = DATA_INTEGER_POS_INFTY; break; }
                case DATA_IS_NEG_INFTY: { result = DATA_INTEGER_NEG_INFTY; break; }
                case DATA_IS_FIN_INDET: { result = DATA_INTEGER_FIN_INDET; break; }
                case DATA_IS_INF_INDET: { result = DATA_INTEGER_INF_INDET; break; }
                case DATA_IS_UNCALCED:  { result = DATA_INTEGER_UNCALCED;  break; }
                case DATA_IS_ERROR:     { result = DATA_INTEGER_ERROR;     break; }
                default:                { result = DATA_INTEGER_ERROR;     break; }
            }

            break;
        }

        case DATA_IS_REAL:
        {
            switch ( b.DataType )
            {
                case DATA_IS_ZERO:      { result = DATA_REAL_ZERO;      break; }
                case DATA_IS_ONE:       { result = DATA_REAL_ONE;       break; }
                case DATA_IS_INTEGER:   { result = DATA_REAL_INTEGER;   break; }
                case DATA_IS_REAL:      { result = DATA_REAL_REAL;      break; }
                case DATA_IS_POS_INFTY: { result = DATA_REAL_POS_INFTY; break; }
                case DATA_IS_NEG_INFTY: { result = DATA_REAL_NEG_INFTY; break; }
                case DATA_IS_FIN_INDET: { result = DATA_REAL_FIN_INDET; break; }
                case DATA_IS_INF_INDET: { result = DATA_REAL_INF_INDET; break; }
                case DATA_IS_UNCALCED:  { result = DATA_REAL_UNCALCED;  break; }
                case DATA_IS_ERROR:     { result = DATA_REAL_ERROR;     break; }
                default:                { result = DATA_REAL_ERROR;     break; }
            }

            break;
        }

        case DATA_IS_POS_INFTY:
        {
            switch ( b.DataType )
            {
                case DATA_IS_ZERO:      { result = DATA_POS_INFTY_ZERO;      break; }
                case DATA_IS_ONE:       { result = DATA_POS_INFTY_ONE;       break; }
                case DATA_IS_INTEGER:   { result = DATA_POS_INFTY_INTEGER;   break; }
                case DATA_IS_REAL:      { result = DATA_POS_INFTY_REAL;      break; }
                case DATA_IS_POS_INFTY: { result = DATA_POS_INFTY_POS_INFTY; break; }
                case DATA_IS_NEG_INFTY: { result = DATA_POS_INFTY_NEG_INFTY; break; }
                case DATA_IS_FIN_INDET: { result = DATA_POS_INFTY_FIN_INDET; break; }
                case DATA_IS_INF_INDET: { result = DATA_POS_INFTY_INF_INDET; break; }
                case DATA_IS_UNCALCED:  { result = DATA_POS_INFTY_UNCALCED;  break; }
                case DATA_IS_ERROR:     { result = DATA_POS_INFTY_ERROR;     break; }
                default:                { result = DATA_POS_INFTY_ERROR;     break; }
            }

            break;
        }

        case DATA_IS_NEG_INFTY:
        {
            switch ( b.DataType )
            {
                case DATA_IS_ZERO:      { result = DATA_NEG_INFTY_ZERO;      break; }
                case DATA_IS_ONE:       { result = DATA_NEG_INFTY_ONE;       break; }
                case DATA_IS_INTEGER:   { result = DATA_NEG_INFTY_INTEGER;   break; }
                case DATA_IS_REAL:      { result = DATA_NEG_INFTY_REAL;      break; }
                case DATA_IS_POS_INFTY: { result = DATA_NEG_INFTY_POS_INFTY; break; }
                case DATA_IS_NEG_INFTY: { result = DATA_NEG_INFTY_NEG_INFTY; break; }
                case DATA_IS_FIN_INDET: { result = DATA_NEG_INFTY_FIN_INDET; break; }
                case DATA_IS_INF_INDET: { result = DATA_NEG_INFTY_INF_INDET; break; }
                case DATA_IS_UNCALCED:  { result = DATA_NEG_INFTY_UNCALCED;  break; }
                case DATA_IS_ERROR:     { result = DATA_NEG_INFTY_ERROR;     break; }
                default:                { result = DATA_NEG_INFTY_ERROR;     break; }
            }

            break;
        }

        case DATA_IS_FIN_INDET:
        {
            switch ( b.DataType )
            {
                case DATA_IS_ZERO:      { result = DATA_FIN_INDET_ZERO;      break; }
                case DATA_IS_ONE:       { result = DATA_FIN_INDET_ONE;       break; }
                case DATA_IS_INTEGER:   { result = DATA_FIN_INDET_INTEGER;   break; }
                case DATA_IS_REAL:      { result = DATA_FIN_INDET_REAL;      break; }
                case DATA_IS_POS_INFTY: { result = DATA_FIN_INDET_POS_INFTY; break; }
                case DATA_IS_NEG_INFTY: { result = DATA_FIN_INDET_NEG_INFTY; break; }
                case DATA_IS_FIN_INDET: { result = DATA_FIN_INDET_FIN_INDET; break; }
                case DATA_IS_INF_INDET: { result = DATA_FIN_INDET_INF_INDET; break; }
                case DATA_IS_UNCALCED:  { result = DATA_FIN_INDET_UNCALCED;  break; }
                case DATA_IS_ERROR:     { result = DATA_FIN_INDET_ERROR;     break; }
                default:                { result = DATA_FIN_INDET_ERROR;     break; }
            }

            break;
        }

        case DATA_IS_INF_INDET:
        {
            switch ( b.DataType )
            {
                case DATA_IS_ZERO:      { result = DATA_INF_INDET_ZERO;      break; }
                case DATA_IS_ONE:       { result = DATA_INF_INDET_ONE;       break; }
                case DATA_IS_INTEGER:   { result = DATA_INF_INDET_INTEGER;   break; }
                case DATA_IS_REAL:      { result = DATA_INF_INDET_REAL;      break; }
                case DATA_IS_POS_INFTY: { result = DATA_INF_INDET_POS_INFTY; break; }
                case DATA_IS_NEG_INFTY: { result = DATA_INF_INDET_NEG_INFTY; break; }
                case DATA_IS_FIN_INDET: { result = DATA_INF_INDET_FIN_INDET; break; }
                case DATA_IS_INF_INDET: { result = DATA_INF_INDET_INF_INDET; break; }
                case DATA_IS_UNCALCED:  { result = DATA_INF_INDET_UNCALCED;  break; }
                case DATA_IS_ERROR:     { result = DATA_INF_INDET_ERROR;     break; }
                default:                { result = DATA_INF_INDET_ERROR;     break; }
            }

            break;
        }

        case DATA_IS_UNCALCED:
        {
            switch ( b.DataType )
            {
                case DATA_IS_ZERO:      { result = DATA_UNCALCED_ZERO;      break; }
                case DATA_IS_ONE:       { result = DATA_UNCALCED_ONE;       break; }
                case DATA_IS_INTEGER:   { result = DATA_UNCALCED_INTEGER;   break; }
                case DATA_IS_REAL:      { result = DATA_UNCALCED_REAL;      break; }
                case DATA_IS_POS_INFTY: { result = DATA_UNCALCED_POS_INFTY; break; }
                case DATA_IS_NEG_INFTY: { result = DATA_UNCALCED_NEG_INFTY; break; }
                case DATA_IS_FIN_INDET: { result = DATA_UNCALCED_FIN_INDET; break; }
                case DATA_IS_INF_INDET: { result = DATA_UNCALCED_INF_INDET; break; }
                case DATA_IS_UNCALCED:  { result = DATA_UNCALCED_UNCALCED;  break; }
                case DATA_IS_ERROR:     { result = DATA_UNCALCED_ERROR;     break; }
                default:                { result = DATA_UNCALCED_ERROR;     break; }
            }

            break;
        }

        case DATA_IS_ERROR:
        default:
        {
            switch ( b.DataType )
            {
                case DATA_IS_ZERO:      { result = DATA_ERROR_ZERO;      break; }
                case DATA_IS_ONE:       { result = DATA_ERROR_ONE;       break; }
                case DATA_IS_INTEGER:   { result = DATA_ERROR_INTEGER;   break; }
                case DATA_IS_REAL:      { result = DATA_ERROR_REAL;      break; }
                case DATA_IS_POS_INFTY: { result = DATA_ERROR_POS_INFTY; break; }
                case DATA_IS_NEG_INFTY: { result = DATA_ERROR_NEG_INFTY; break; }
                case DATA_IS_FIN_INDET: { result = DATA_ERROR_FIN_INDET; break; }
                case DATA_IS_INF_INDET: { result = DATA_ERROR_INF_INDET; break; }
                case DATA_IS_UNCALCED:  { result = DATA_ERROR_UNCALCED;  break; }
                case DATA_IS_ERROR:     { result = DATA_ERROR_ERROR;     break; }
                default:                { result = DATA_ERROR_ERROR;     break; }
            }

            break;
        }
    }

    return result;
}


GenVar safeIntAdd(GVAR_INT a, GVAR_INT b)
{
    GenVar result;

    if ( ( a > 0 ) && ( b > 0 ) )
    {
        if ( a <= GVAR_INT_MAX-b )
        {
            result = ResultInteger(a+b);
        }

        else
        {
            result = safeRealAdd((double) a,(double) b);
        }
    }

    else
    {
        if ( ( a < 0 ) && ( b < 0 ) )
        {
            if ( a >= GVAR_INT_MIN-b )
            {
                result = ResultInteger(a+b);
            }

            else
            {
                result = safeRealAdd((double) a,(double) b);
            }
        }

        else
        {
            result = ResultInteger(a+b);
        }
    }

    return simplifyGenVar(result);
}

GenVar safeIntMul(GVAR_INT a, GVAR_INT b)
{
    GVAR_INT sign_result;
    unsigned GVAR_INT mag_a;
    unsigned GVAR_INT mag_b;
    unsigned GVAR_INT mag_res;
    unsigned GVAR_INT a_lsb;
    unsigned GVAR_INT a_msb;
    unsigned GVAR_INT b_lsb;
    unsigned GVAR_INT b_msb;
    unsigned GVAR_INT ares_lsb;
    unsigned GVAR_INT ares_msb;
    unsigned GVAR_INT ares_lsb_over;
    unsigned GVAR_INT ares_msb_over;
    unsigned GVAR_INT bres_lsb;
    unsigned GVAR_INT bres_msb;
    unsigned GVAR_INT bres_lsb_over;
    unsigned GVAR_INT bres_msb_over;
    unsigned GVAR_INT final_over;
    unsigned GVAR_INT lsb_mask;

    /*
       Do special cases first.
    */

    if ( ( a == 0 ) || ( b == 0 ) )
    {
        return ResultZero();
    }

    if ( a == 1 )
    {
        return ResultInteger(b);
    }

    if ( b == 1 )
    {
        return ResultInteger(a);
    }

    /*
       Don't bother if it's too close to call.
    */

    if ( ( a >= GVAR_INT_MAX ) || ( a <= GVAR_INT_MIN ) || 
         ( b >= GVAR_INT_MAX ) || ( b <= GVAR_INT_MIN )    )
    {
        return safeRealMul((double) a,(double) b);
    }

    /*
       Assumptions: 1 byte = 8 bits
                    2 <= |a| <= min(GVAR_INT_MAX,-GVAR_INT_MAX)
                    2 <= |b| <= min(GVAR_INT_MAX,-GVAR_INT_MAX)
                    |GVAR_INT_MAX+GVAR_INT_MIN| <= 1
                    sizeof returns # bytes
    */

    lsb_mask = (((unsigned GVAR_INT) 0x01)>>(4*sizeof(unsigned GVAR_INT)))-1;

    sign_result = 1;

    if ( a < 0 ) { sign_result *= -1; mag_a = (unsigned GVAR_INT) -a; } else { mag_a = (unsigned GVAR_INT) a; }
    if ( b < 0 ) { sign_result *= -1; mag_b = (unsigned GVAR_INT) -b; } else { mag_b = (unsigned GVAR_INT) b; }

    a_lsb = mag_a & lsb_mask;
    b_lsb = mag_b & lsb_mask;

    a_msb = ( mag_a >> (4*sizeof(unsigned GVAR_INT)) ) & lsb_mask;
    b_msb = ( mag_b >> (4*sizeof(unsigned GVAR_INT)) ) & lsb_mask;

    /*
                      b_msb    b_lsb
                        *      a_lsb
                     -----------------
       ares_msb_over ares_msb ares_lsb

       overflow if ares_msb_over != 0
    */

    ares_lsb = a_lsb * b_lsb;
    ares_msb = a_lsb * b_msb;

    ares_lsb_over = ( ares_lsb >> (4*sizeof(unsigned GVAR_INT)) ) & lsb_mask;
    ares_lsb      &= lsb_mask;

    ares_msb += ares_lsb_over;

    ares_msb_over = ( ares_msb >> (4*sizeof(unsigned GVAR_INT)) ) & lsb_mask;
    ares_msb      &= lsb_mask;

    if ( ares_msb_over != 0 )
    {
        return safeRealMul((double) a,(double) b);
    }

    /*
                               b_msb    b_lsb
                        *      a_msb
                     -------------------------
       bres_msb_over bres_msb bres_lsb   0

       overflow if bres_msb_over != 0
       overflow if bres_msb      != 0
    */

    bres_lsb = a_msb * b_lsb;
    bres_msb = a_msb * b_msb;

    bres_lsb_over = ( bres_lsb >> (4*sizeof(unsigned GVAR_INT)) ) & lsb_mask;
    bres_lsb      &= lsb_mask;

    bres_msb += bres_lsb_over;

    bres_msb_over = ( bres_msb >> (4*sizeof(unsigned GVAR_INT)) ) & lsb_mask;
    bres_msb      &= lsb_mask;

    if ( ( bres_msb_over != 0 ) || ( bres_msb != 0 ) )
    {
        return safeRealMul((double) a,(double) b);
    }

    /*
                   ares_msb ares_lsb
                 + bres_lsb   0
                   -----------------
       final_over    ...    ares_lsb

       overflow if final_over > 0
    */

    bres_lsb += ares_msb;

    final_over = ( bres_lsb >> (4*sizeof(unsigned GVAR_INT)) ) & lsb_mask;

    if ( final_over != 0 )
    {
        return safeRealMul((double) a,(double) b);
    }

    /*
       OK, mag_res is the result - does it overflow GVAR_INT?
    */

    mag_res = mag_a * mag_b;

    if ( ( mag_res >= GVAR_INT_MAX ) && ( sign_result == +1 ) )
    {
        return ResultReal((double) mag_res);
    }

    if ( ( mag_res >= GVAR_INT_MAX ) && ( sign_result == -1 ) )
    {
        return ResultReal(-((double) mag_res));
    }

    return ResultInteger(sign_result*((GVAR_INT) mag_res));
}


GenVar safeRealAdd(double a, double b)
{
    GenVar result;

    if ( ( a > 0 ) && ( b > 0 ) )
    {
        if ( a <= GSL_DBL_MAX-b )
        {
            result = ResultReal(a+b);
        }

        else
        {
            result = ResultPosInfnty();
        }
    }

    else
    {
        if ( ( a < 0 ) && ( b < 0 ) )
        {
            if ( a >= -GSL_DBL_MAX-b )
            {
                result = ResultReal(a+b);
            }

            else
            {
                result = ResultNegInfnty();
            }
        }

        else
        {
            result = ResultReal(a+b);
        }
    }

    return simplifyGenVar(result);
}

GenVar safeRealMul(double a, double b)
{
    GenVar result;
    double bas_a;
    double bas_b;
    double bas_res;
    double bas_max;
    int exp_a;
    int exp_b;
    int exp_res;
    int exp_max;

    if ( ( a == 0 ) || ( b == 0 ) )
    {
        result = ResultZero();
    }

    else
    {
        bas_max = gsl_frexp(GSL_DBL_MAX,&exp_max);

        bas_a = gsl_frexp(a,&exp_a);
        bas_b = gsl_frexp(b,&exp_b);

        bas_res = bas_a * bas_b;
        exp_res = exp_a + exp_b;

        bas_a = gsl_frexp(bas_res,&exp_a);

        bas_res =  bas_a;
        exp_res += exp_a;

        if ( ( exp_res > exp_max ) || ( ( exp_res == exp_max ) && ( GSL_MAX(bas_res,-bas_res) >= bas_max ) ) )
        {
            result = ResultPosInfnty();
        }

        else
        {
            result = ResultReal(gsl_ldexp(bas_res,exp_res));
        }
    }

    return result;
}



unsigned long long convGSLErr(int gsl_errno)
{
    unsigned long long result = 0;

    switch ( gsl_errno )
    {
        case GSL_FAILURE:  { result = (unsigned long long) GVAR_GSL_FAILUR; break; }
        case GSL_CONTINUE: { result = (unsigned long long) GVAR_GSL_CONTIN; break; }
        case GSL_EDOM:     { result = (unsigned long long) GVAR_GSL_EDOM;   break; }
        case GSL_ERANGE:   { result = (unsigned long long) GVAR_GSL_ERANGE; break; }
        case GSL_EFAULT:   { result = (unsigned long long) GVAR_GSL_EFAULT; break; }
        case GSL_EINVAL:   { result = (unsigned long long) GVAR_GSL_EINVAL; break; }
        case GSL_EFAILED:  { result = (unsigned long long) GVAR_GSL_EFAILE; break; }
        case GSL_EFACTOR:  { result = (unsigned long long) GVAR_GSL_EFACTO; break; }
        case GSL_ESANITY:  { result = (unsigned long long) GVAR_GSL_ESANIT; break; }
        case GSL_ENOMEM:   { result = (unsigned long long) GVAR_GSL_ENOMEM; break; }
        case GSL_EBADFUNC: { result = (unsigned long long) GVAR_GSL_EBADFN; break; }
        case GSL_ERUNAWAY: { result = (unsigned long long) GVAR_GSL_ERNWAY; break; }
        case GSL_EMAXITER: { result = (unsigned long long) GVAR_GSL_EMAXIT; break; }
        case GSL_EZERODIV: { result = (unsigned long long) GVAR_GSL_EZRDIV; break; }
        case GSL_EBADTOL:  { result = (unsigned long long) GVAR_GSL_EBADTL; break; }
        case GSL_ETOL:     { result = (unsigned long long) GVAR_GSL_ETOL;   break; }
        case GSL_EUNDRFLW: { result = (unsigned long long) GVAR_GSL_EUNDRF; break; }
        case GSL_EOVRFLW:  { result = (unsigned long long) GVAR_GSL_EOVRFL; break; }
        case GSL_ELOSS:    { result = (unsigned long long) GVAR_GSL_ELOSS;  break; }
        case GSL_EROUND:   { result = (unsigned long long) GVAR_GSL_EROUND; break; }
        case GSL_EBADLEN:  { result = (unsigned long long) GVAR_GSL_EBADLE; break; }
        case GSL_ENOTSQR:  { result = (unsigned long long) GVAR_GSL_ENOTSQ; break; }
        case GSL_ESING:    { result = (unsigned long long) GVAR_GSL_ESING;  break; }
        case GSL_EDIVERGE: { result = (unsigned long long) GVAR_GSL_EDIVER; break; }
        case GSL_EUNSUP:   { result = (unsigned long long) GVAR_GSL_EUNSUP; break; }
        case GSL_EUNIMPL:  { result = (unsigned long long) GVAR_GSL_EUNIMP; break; }
        case GSL_ECACHE:   { result = (unsigned long long) GVAR_GSL_ECACHE; break; }
        case GSL_ETABLE:   { result = (unsigned long long) GVAR_GSL_ETABLE; break; }
        case GSL_ENOPROG:  { result = (unsigned long long) GVAR_GSL_ENOPRO; break; }
        case GSL_ENOPROGJ: { result = (unsigned long long) GVAR_GSL_ENOPRX; break; }
        case GSL_ETOLF:    { result = (unsigned long long) GVAR_GSL_ETOLF;  break; }
        case GSL_ETOLX:    { result = (unsigned long long) GVAR_GSL_ETOLX;  break; }
        case GSL_ETOLG:    { result = (unsigned long long) GVAR_GSL_ETOLG;  break; }
        case GSL_EOF:      { result = (unsigned long long) GVAR_GSL_EOF;    break; }

        default:
        {
            break;
        }
    }

    return result;
}

GenVar combineErrors(GenVar a, GenVar b)
{
    a.ErrorCode |= b.ErrorCode;

    return a;
}


