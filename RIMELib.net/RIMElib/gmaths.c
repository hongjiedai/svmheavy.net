#include "stdafx.h"
/*
 *  RIMElib: RuntIme Mathematical Equation Library (underlying functions)
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

#include <gsl/specfunc/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <gsl/err/gsl_errno.h>

#include "gmaths.h"
#include "gvars.h"

GenVar Gen_kronDelta(GenVar ii, GenVar jj)
{
    GenVar result;

    switch ( combineDataDescr(ii,jj) )
    {
        case DATA_ERROR_ZERO:
        case DATA_ERROR_ONE:
        case DATA_ERROR_INTEGER:
        case DATA_ERROR_UNCALCED:
        case DATA_UNCALCED_ZERO:
        case DATA_UNCALCED_ONE:
        case DATA_UNCALCED_INTEGER:
        case DATA_UNCALCED_UNCALCED:
        {
            result = ii;

            break;
        }

        case DATA_ZERO_ERROR:
        case DATA_ONE_ERROR:
        case DATA_INTEGER_ERROR:
        case DATA_UNCALCED_ERROR:
        case DATA_ZERO_UNCALCED:
        case DATA_ONE_UNCALCED:
        case DATA_INTEGER_UNCALCED:
        {
            result = jj;

            break;
        }

        case DATA_ZERO_ZERO:
        case DATA_ONE_ZERO:
        case DATA_INTEGER_ZERO:
        case DATA_ZERO_ONE:
        case DATA_ONE_ONE:
        case DATA_INTEGER_ONE:
        case DATA_ZERO_INTEGER:
        case DATA_ONE_INTEGER:
        case DATA_INTEGER_INTEGER:
        {
            if ( isaequalb(ii,jj) )
            {
                result = ResultOne();
            }

            else
            {
                result = ResultZero();
            }

            break;
        }

        case DATA_ERROR_REAL:
        case DATA_ERROR_POS_INFTY:
        case DATA_ERROR_NEG_INFTY:
        case DATA_ERROR_FIN_INDET:
        case DATA_ERROR_INF_INDET:
        {
            result = combineErrors(ii,ResultError(GVAR_GSL_EDOM));

            break;
        }

        case DATA_REAL_ERROR:
        case DATA_POS_INFTY_ERROR:
        case DATA_NEG_INFTY_ERROR:
        case DATA_FIN_INDET_ERROR:
        case DATA_INF_INDET_ERROR:
        {
            result = combineErrors(ResultError(GVAR_GSL_EDOM),jj);

            break;
        }

        case DATA_UNCALCED_REAL:
        case DATA_UNCALCED_POS_INFTY:
        case DATA_UNCALCED_NEG_INFTY:
        case DATA_UNCALCED_FIN_INDET:
        case DATA_UNCALCED_INF_INDET:
        case DATA_REAL_UNCALCED:
        case DATA_POS_INFTY_UNCALCED:
        case DATA_NEG_INFTY_UNCALCED:
        case DATA_FIN_INDET_UNCALCED:
        case DATA_INF_INDET_UNCALCED:
        case DATA_REAL_ZERO:
        case DATA_POS_INFTY_ZERO:
        case DATA_NEG_INFTY_ZERO:
        case DATA_FIN_INDET_ZERO:
        case DATA_INF_INDET_ZERO:
        case DATA_REAL_ONE:
        case DATA_POS_INFTY_ONE:
        case DATA_NEG_INFTY_ONE:
        case DATA_FIN_INDET_ONE:
        case DATA_INF_INDET_ONE:
        case DATA_REAL_INTEGER:
        case DATA_POS_INFTY_INTEGER:
        case DATA_NEG_INFTY_INTEGER:
        case DATA_FIN_INDET_INTEGER:
        case DATA_INF_INDET_INTEGER:
        case DATA_ZERO_REAL:
        case DATA_ONE_REAL:
        case DATA_INTEGER_REAL:
        case DATA_REAL_REAL:
        case DATA_POS_INFTY_REAL:
        case DATA_NEG_INFTY_REAL:
        case DATA_FIN_INDET_REAL:
        case DATA_INF_INDET_REAL:
        case DATA_ZERO_POS_INFTY:
        case DATA_ONE_POS_INFTY:
        case DATA_INTEGER_POS_INFTY:
        case DATA_REAL_POS_INFTY:
        case DATA_POS_INFTY_POS_INFTY:
        case DATA_NEG_INFTY_POS_INFTY:
        case DATA_FIN_INDET_POS_INFTY:
        case DATA_INF_INDET_POS_INFTY:
        case DATA_ZERO_NEG_INFTY:
        case DATA_ONE_NEG_INFTY:
        case DATA_INTEGER_NEG_INFTY:
        case DATA_REAL_NEG_INFTY:
        case DATA_POS_INFTY_NEG_INFTY:
        case DATA_NEG_INFTY_NEG_INFTY:
        case DATA_FIN_INDET_NEG_INFTY:
        case DATA_INF_INDET_NEG_INFTY:
        case DATA_ZERO_FIN_INDET:
        case DATA_ONE_FIN_INDET:
        case DATA_INTEGER_FIN_INDET:
        case DATA_REAL_FIN_INDET:
        case DATA_POS_INFTY_FIN_INDET:
        case DATA_NEG_INFTY_FIN_INDET:
        case DATA_FIN_INDET_FIN_INDET:
        case DATA_INF_INDET_FIN_INDET:
        case DATA_ZERO_INF_INDET:
        case DATA_ONE_INF_INDET:
        case DATA_INTEGER_INF_INDET:
        case DATA_REAL_INF_INDET:
        case DATA_POS_INFTY_INF_INDET:
        case DATA_NEG_INFTY_INF_INDET:
        case DATA_FIN_INDET_INF_INDET:
        case DATA_INF_INDET_INF_INDET:
        {
            result = ResultError(GVAR_GSL_EDOM);

            break;
        }

        case DATA_ERROR_ERROR:
        default:
        {
            result = combineErrors(ii,jj);

            break;
        }
    }

    return result;
}

GenVar Gen_diracDelta(GenVar xx, GenVar yy)
{
    GenVar result;

    switch ( combineDataDescr(xx,yy) )
    {
        case DATA_ERROR_ZERO:
        case DATA_ERROR_ONE:
        case DATA_ERROR_INTEGER:
        case DATA_ERROR_REAL:
        case DATA_ERROR_UNCALCED:
        case DATA_UNCALCED_ZERO:
        case DATA_UNCALCED_ONE:
        case DATA_UNCALCED_INTEGER:
        case DATA_UNCALCED_REAL:
        case DATA_UNCALCED_UNCALCED:
        {
            result = xx;

            break;
        }

        case DATA_ZERO_ERROR:
        case DATA_ONE_ERROR:
        case DATA_INTEGER_ERROR:
        case DATA_REAL_ERROR:
        case DATA_UNCALCED_ERROR:
        case DATA_ZERO_UNCALCED:
        case DATA_ONE_UNCALCED:
        case DATA_INTEGER_UNCALCED:
        case DATA_REAL_UNCALCED:
        {
            result = yy;

            break;
        }

        case DATA_ZERO_ZERO:
        case DATA_ONE_ZERO:
        case DATA_INTEGER_ZERO:
        case DATA_REAL_ZERO:
        case DATA_ZERO_ONE:
        case DATA_ONE_ONE:
        case DATA_INTEGER_ONE:
        case DATA_REAL_ONE:
        case DATA_ZERO_INTEGER:
        case DATA_ONE_INTEGER:
        case DATA_INTEGER_INTEGER:
        case DATA_REAL_INTEGER:
        case DATA_ZERO_REAL:
        case DATA_ONE_REAL:
        case DATA_INTEGER_REAL:
        case DATA_REAL_REAL:
        {
            if ( isaequalb(xx,yy) )
            {
                result = ResultPosInfnty();
            }

            else
            {
                result = ResultZero();
            }

            break;
        }

        case DATA_ERROR_POS_INFTY:
        case DATA_ERROR_NEG_INFTY:
        case DATA_ERROR_FIN_INDET:
        case DATA_ERROR_INF_INDET:
        {
            result = combineErrors(xx,ResultError(GVAR_GSL_EDOM));

            break;
        }

        case DATA_POS_INFTY_ERROR:
        case DATA_NEG_INFTY_ERROR:
        case DATA_FIN_INDET_ERROR:
        case DATA_INF_INDET_ERROR:
        {
            result = combineErrors(ResultError(GVAR_GSL_EDOM),yy);

            break;
        }

        case DATA_UNCALCED_POS_INFTY:
        case DATA_UNCALCED_NEG_INFTY:
        case DATA_UNCALCED_FIN_INDET:
        case DATA_UNCALCED_INF_INDET:
        case DATA_POS_INFTY_UNCALCED:
        case DATA_NEG_INFTY_UNCALCED:
        case DATA_FIN_INDET_UNCALCED:
        case DATA_INF_INDET_UNCALCED:
        case DATA_POS_INFTY_ZERO:
        case DATA_NEG_INFTY_ZERO:
        case DATA_FIN_INDET_ZERO:
        case DATA_INF_INDET_ZERO:
        case DATA_POS_INFTY_ONE:
        case DATA_NEG_INFTY_ONE:
        case DATA_FIN_INDET_ONE:
        case DATA_INF_INDET_ONE:
        case DATA_POS_INFTY_INTEGER:
        case DATA_NEG_INFTY_INTEGER:
        case DATA_FIN_INDET_INTEGER:
        case DATA_INF_INDET_INTEGER:
        case DATA_POS_INFTY_REAL:
        case DATA_NEG_INFTY_REAL:
        case DATA_FIN_INDET_REAL:
        case DATA_INF_INDET_REAL:
        case DATA_ZERO_POS_INFTY:
        case DATA_ONE_POS_INFTY:
        case DATA_INTEGER_POS_INFTY:
        case DATA_REAL_POS_INFTY:
        case DATA_POS_INFTY_POS_INFTY:
        case DATA_NEG_INFTY_POS_INFTY:
        case DATA_FIN_INDET_POS_INFTY:
        case DATA_INF_INDET_POS_INFTY:
        case DATA_ZERO_NEG_INFTY:
        case DATA_ONE_NEG_INFTY:
        case DATA_INTEGER_NEG_INFTY:
        case DATA_REAL_NEG_INFTY:
        case DATA_POS_INFTY_NEG_INFTY:
        case DATA_NEG_INFTY_NEG_INFTY:
        case DATA_FIN_INDET_NEG_INFTY:
        case DATA_INF_INDET_NEG_INFTY:
        case DATA_ZERO_FIN_INDET:
        case DATA_ONE_FIN_INDET:
        case DATA_INTEGER_FIN_INDET:
        case DATA_REAL_FIN_INDET:
        case DATA_POS_INFTY_FIN_INDET:
        case DATA_NEG_INFTY_FIN_INDET:
        case DATA_FIN_INDET_FIN_INDET:
        case DATA_INF_INDET_FIN_INDET:
        case DATA_ZERO_INF_INDET:
        case DATA_ONE_INF_INDET:
        case DATA_INTEGER_INF_INDET:
        case DATA_REAL_INF_INDET:
        case DATA_POS_INFTY_INF_INDET:
        case DATA_NEG_INFTY_INF_INDET:
        case DATA_FIN_INDET_INF_INDET:
        case DATA_INF_INDET_INF_INDET:
        {
            result = ResultError(GVAR_GSL_EDOM);

            break;
        }

        case DATA_ERROR_ERROR:
        default:
        {
            result = combineErrors(xx,yy);

            break;
        }
    }

    return result;
}

GenVar Gen_finDiracDelta(GenVar xx, GenVar yy)
{
    GenVar result;

    switch ( combineDataDescr(xx,yy) )
    {
        case DATA_ERROR_ZERO:
        case DATA_ERROR_ONE:
        case DATA_ERROR_INTEGER:
        case DATA_ERROR_REAL:
        case DATA_ERROR_UNCALCED:
        case DATA_UNCALCED_ZERO:
        case DATA_UNCALCED_ONE:
        case DATA_UNCALCED_INTEGER:
        case DATA_UNCALCED_REAL:
        case DATA_UNCALCED_UNCALCED:
        {
            result = xx;

            break;
        }

        case DATA_ZERO_ERROR:
        case DATA_ONE_ERROR:
        case DATA_INTEGER_ERROR:
        case DATA_REAL_ERROR:
        case DATA_UNCALCED_ERROR:
        case DATA_ZERO_UNCALCED:
        case DATA_ONE_UNCALCED:
        case DATA_INTEGER_UNCALCED:
        case DATA_REAL_UNCALCED:
        {
            result = yy;

            break;
        }

        case DATA_ZERO_ZERO:
        case DATA_ONE_ZERO:
        case DATA_INTEGER_ZERO:
        case DATA_REAL_ZERO:
        case DATA_ZERO_ONE:
        case DATA_ONE_ONE:
        case DATA_INTEGER_ONE:
        case DATA_REAL_ONE:
        case DATA_ZERO_INTEGER:
        case DATA_ONE_INTEGER:
        case DATA_INTEGER_INTEGER:
        case DATA_REAL_INTEGER:
        case DATA_ZERO_REAL:
        case DATA_ONE_REAL:
        case DATA_INTEGER_REAL:
        case DATA_REAL_REAL:
        {
            if ( isaequalb(xx,yy) )
            {
                result = ResultOne();
            }

            else
            {
                result = ResultZero();
            }

            break;
        }

        case DATA_ERROR_POS_INFTY:
        case DATA_ERROR_NEG_INFTY:
        case DATA_ERROR_FIN_INDET:
        case DATA_ERROR_INF_INDET:
        {
            result = combineErrors(xx,ResultError(GVAR_GSL_EDOM));

            break;
        }

        case DATA_POS_INFTY_ERROR:
        case DATA_NEG_INFTY_ERROR:
        case DATA_FIN_INDET_ERROR:
        case DATA_INF_INDET_ERROR:
        {
            result = combineErrors(ResultError(GVAR_GSL_EDOM),yy);

            break;
        }

        case DATA_UNCALCED_POS_INFTY:
        case DATA_UNCALCED_NEG_INFTY:
        case DATA_UNCALCED_FIN_INDET:
        case DATA_UNCALCED_INF_INDET:
        case DATA_POS_INFTY_UNCALCED:
        case DATA_NEG_INFTY_UNCALCED:
        case DATA_FIN_INDET_UNCALCED:
        case DATA_INF_INDET_UNCALCED:
        case DATA_POS_INFTY_ZERO:
        case DATA_NEG_INFTY_ZERO:
        case DATA_FIN_INDET_ZERO:
        case DATA_INF_INDET_ZERO:
        case DATA_POS_INFTY_ONE:
        case DATA_NEG_INFTY_ONE:
        case DATA_FIN_INDET_ONE:
        case DATA_INF_INDET_ONE:
        case DATA_POS_INFTY_INTEGER:
        case DATA_NEG_INFTY_INTEGER:
        case DATA_FIN_INDET_INTEGER:
        case DATA_INF_INDET_INTEGER:
        case DATA_POS_INFTY_REAL:
        case DATA_NEG_INFTY_REAL:
        case DATA_FIN_INDET_REAL:
        case DATA_INF_INDET_REAL:
        case DATA_ZERO_POS_INFTY:
        case DATA_ONE_POS_INFTY:
        case DATA_INTEGER_POS_INFTY:
        case DATA_REAL_POS_INFTY:
        case DATA_POS_INFTY_POS_INFTY:
        case DATA_NEG_INFTY_POS_INFTY:
        case DATA_FIN_INDET_POS_INFTY:
        case DATA_INF_INDET_POS_INFTY:
        case DATA_ZERO_NEG_INFTY:
        case DATA_ONE_NEG_INFTY:
        case DATA_INTEGER_NEG_INFTY:
        case DATA_REAL_NEG_INFTY:
        case DATA_POS_INFTY_NEG_INFTY:
        case DATA_NEG_INFTY_NEG_INFTY:
        case DATA_FIN_INDET_NEG_INFTY:
        case DATA_INF_INDET_NEG_INFTY:
        case DATA_ZERO_FIN_INDET:
        case DATA_ONE_FIN_INDET:
        case DATA_INTEGER_FIN_INDET:
        case DATA_REAL_FIN_INDET:
        case DATA_POS_INFTY_FIN_INDET:
        case DATA_NEG_INFTY_FIN_INDET:
        case DATA_FIN_INDET_FIN_INDET:
        case DATA_INF_INDET_FIN_INDET:
        case DATA_ZERO_INF_INDET:
        case DATA_ONE_INF_INDET:
        case DATA_INTEGER_INF_INDET:
        case DATA_REAL_INF_INDET:
        case DATA_POS_INFTY_INF_INDET:
        case DATA_NEG_INFTY_INF_INDET:
        case DATA_FIN_INDET_INF_INDET:
        case DATA_INF_INDET_INF_INDET:
        {
            result = ResultError(GVAR_GSL_EDOM);

            break;
        }

        case DATA_ERROR_ERROR:
        default:
        {
            result = combineErrors(xx,yy);

            break;
        }
    }

    return result;
}

GenVar Gen_max(GenVar xx, GenVar yy)
{
    GenVar result;

    switch ( combineDataDescr(xx,yy) )
    {
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
        case DATA_UNCALCED_NEG_INFTY:
        case DATA_UNCALCED_FIN_INDET:
        case DATA_UNCALCED_INF_INDET:
        case DATA_UNCALCED_UNCALCED:
        case DATA_FIN_INDET_ZERO:
        case DATA_FIN_INDET_ONE:
        case DATA_FIN_INDET_INTEGER:
        case DATA_FIN_INDET_REAL:
        case DATA_FIN_INDET_NEG_INFTY:
        case DATA_FIN_INDET_FIN_INDET:
        case DATA_POS_INFTY_FIN_INDET:
        case DATA_INF_INDET_FIN_INDET:
        case DATA_INF_INDET_ZERO:
        case DATA_INF_INDET_ONE:
        case DATA_INF_INDET_INTEGER:
        case DATA_INF_INDET_REAL:
        case DATA_INF_INDET_NEG_INFTY:
        case DATA_POS_INFTY_INF_INDET:
        case DATA_POS_INFTY_UNCALCED:
        {
            result = xx;

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
        case DATA_UNCALCED_ERROR:
        case DATA_ZERO_UNCALCED:
        case DATA_ONE_UNCALCED:
        case DATA_INTEGER_UNCALCED:
        case DATA_REAL_UNCALCED:
        case DATA_NEG_INFTY_UNCALCED:
        case DATA_FIN_INDET_UNCALCED:
        case DATA_INF_INDET_UNCALCED:
        case DATA_FIN_INDET_POS_INFTY:
        case DATA_ZERO_FIN_INDET:
        case DATA_ONE_FIN_INDET:
        case DATA_INTEGER_FIN_INDET:
        case DATA_REAL_FIN_INDET:
        case DATA_NEG_INFTY_FIN_INDET:
        case DATA_INF_INDET_POS_INFTY:
        case DATA_ZERO_INF_INDET:
        case DATA_ONE_INF_INDET:
        case DATA_INTEGER_INF_INDET:
        case DATA_REAL_INF_INDET:
        case DATA_NEG_INFTY_INF_INDET:
        case DATA_FIN_INDET_INF_INDET:
        case DATA_INF_INDET_INF_INDET:
        case DATA_UNCALCED_POS_INFTY:
        {
            result = yy;

            break;
        }

        case DATA_ZERO_ZERO:
        case DATA_ONE_ZERO:
        case DATA_INTEGER_ZERO:
        case DATA_REAL_ZERO:
        case DATA_POS_INFTY_ZERO:
        case DATA_NEG_INFTY_ZERO:
        case DATA_ZERO_ONE:
        case DATA_ONE_ONE:
        case DATA_INTEGER_ONE:
        case DATA_REAL_ONE:
        case DATA_POS_INFTY_ONE:
        case DATA_NEG_INFTY_ONE:
        case DATA_ZERO_INTEGER:
        case DATA_ONE_INTEGER:
        case DATA_INTEGER_INTEGER:
        case DATA_REAL_INTEGER:
        case DATA_POS_INFTY_INTEGER:
        case DATA_NEG_INFTY_INTEGER:
        case DATA_ZERO_REAL:
        case DATA_ONE_REAL:
        case DATA_INTEGER_REAL:
        case DATA_REAL_REAL:
        case DATA_POS_INFTY_REAL:
        case DATA_NEG_INFTY_REAL:
        case DATA_ZERO_POS_INFTY:
        case DATA_ONE_POS_INFTY:
        case DATA_INTEGER_POS_INFTY:
        case DATA_REAL_POS_INFTY:
        case DATA_POS_INFTY_POS_INFTY:
        case DATA_NEG_INFTY_POS_INFTY:
        case DATA_ZERO_NEG_INFTY:
        case DATA_ONE_NEG_INFTY:
        case DATA_INTEGER_NEG_INFTY:
        case DATA_REAL_NEG_INFTY:
        case DATA_POS_INFTY_NEG_INFTY:
        case DATA_NEG_INFTY_NEG_INFTY:
        {
            if ( isalessb(xx,yy) )
            {
                result = yy;
            }

            else
            {
                result = xx;
            }

            break;
        }

        case DATA_ERROR_ERROR:
        default:
        {
            result = combineErrors(xx,yy);

            break;
        }
    }

    return result;
}

GenVar Gen_min(GenVar xx, GenVar yy)
{
    return negGenVar(Gen_max(negGenVar(xx),negGenVar(yy)));
}

GenVar Gen_abs(GenVar xx)
{
    GenVar result;

    switch ( xx.DataType )
    {
        case DATA_IS_ZERO:
        case DATA_IS_ONE:
        {
            result = xx;

            break;
        }

        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        {
            /* this is quite deliberate - remember overflows */

            if ( getReal(xx) < 0 )
            {
                result = negGenVar(xx);
            }

            else
            {
                result = xx;
            }

            break;
        }

        case DATA_IS_POS_INFTY:
        case DATA_IS_NEG_INFTY:
        {
            result = ResultPosInfnty();

            break;
        }

        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        case DATA_IS_UNCALCED:
        case DATA_IS_ERROR:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_sgn(GenVar xx)
{
    GenVar result;

    switch ( xx.DataType )
    {
        case DATA_IS_ZERO:
        {
            result = ResultZero();

            break;
        }

        case DATA_IS_ONE:
        {
            result = ResultOne();

            break;
        }

        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        {
            if ( getReal(xx) > 0 )
            {
                result = ResultOne();
            }

            else
            {
                result = ResultInteger(-1);
            }

            break;
        }

        case DATA_IS_POS_INFTY:
        {
            result = ResultOne();

            break;
        }

        case DATA_IS_NEG_INFTY:
        {
            result = ResultInteger(-1);

            break;
        }

        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        {
            result = ResultFinIndet();

            break;
        }

        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_rint(GenVar xx)
{
    GenVar result;
    double tempa;
    double tempb;
    double tempc;

    switch ( xx.DataType )
    {
        case DATA_IS_ZERO:
        case DATA_IS_ONE:
        case DATA_IS_INTEGER:
        {
            result = xx;

            break;
        }

        case DATA_IS_REAL:
        {
            tempa = getReal(xx);

            tempb = ceil(tempa);
            tempc = floor(tempa);

            if ( (tempa-tempc) <= (tempb-tempa) )
            {
                if ( ( tempc < GVAR_INT_MAX ) && ( tempc > GVAR_INT_MIN ) )
                {
                    result = ResultInteger(((GVAR_INT) tempc));
                }

                else
                {
                    result = ResultReal(tempc);
                }
            }

            else
            {
                if ( ( tempb < GVAR_INT_MAX ) && ( tempb > GVAR_INT_MIN ) )
                {
                    result = ResultInteger(((GVAR_INT) tempb));
                }

                else
                {
                    result = ResultReal(tempb);
                }
            }

            break;
        }

        case DATA_IS_POS_INFTY:
        case DATA_IS_NEG_INFTY:
        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_ceil(GenVar xx)
{
    GenVar result;
    double tempa;

    switch ( xx.DataType )
    {
        case DATA_IS_ZERO:
        case DATA_IS_ONE:
        case DATA_IS_INTEGER:
        {
            result = xx;

            break;
        }

        case DATA_IS_REAL:
        {
            tempa = ceil(getReal(xx));

            if ( ( tempa < GVAR_INT_MAX ) && ( tempa > GVAR_INT_MIN ) )
            {
                result = ResultInteger((GVAR_INT) tempa);
            }

            else
            {
                result = ResultReal(tempa);
            }

            break;
        }

        case DATA_IS_POS_INFTY:
        case DATA_IS_NEG_INFTY:
        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_floor(GenVar xx)
{
    GenVar result;
    double tempa;

    switch ( xx.DataType )
    {
        case DATA_IS_ZERO:
        case DATA_IS_ONE:
        case DATA_IS_INTEGER:
        {
            result = xx;

            break;
        }

        case DATA_IS_REAL:
        {
            tempa = floor(getReal(xx));

            if ( ( tempa < GVAR_INT_MAX ) && ( tempa > GVAR_INT_MIN ) )
            {
                result = ResultInteger((GVAR_INT) tempa);
            }

            else
            {
                result = ResultReal(tempa);
            }

            break;
        }

        case DATA_IS_POS_INFTY:
        case DATA_IS_NEG_INFTY:
        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_perm(GenVar ii, GenVar jj)
{
    GenVar result;
    GVAR_INT i;
    GVAR_INT j;
    GVAR_INT k;

    switch ( combineDataDescr(ii,jj) )
    {
        case DATA_ERROR_ZERO:
        case DATA_ERROR_ONE:
        case DATA_ERROR_INTEGER:
        case DATA_ERROR_UNCALCED:
        case DATA_UNCALCED_ZERO:
        case DATA_UNCALCED_ONE:
        case DATA_UNCALCED_INTEGER:
        case DATA_UNCALCED_UNCALCED:
        {
            result = ii;

            break;
        }

        case DATA_ZERO_ERROR:
        case DATA_ONE_ERROR:
        case DATA_INTEGER_ERROR:
        case DATA_UNCALCED_ERROR:
        case DATA_ZERO_UNCALCED:
        case DATA_ONE_UNCALCED:
        case DATA_INTEGER_UNCALCED:
        {
            result = jj;

            break;
        }

        case DATA_ZERO_ZERO:
        case DATA_ONE_ZERO:
        case DATA_INTEGER_ZERO:
        case DATA_ZERO_ONE:
        case DATA_ONE_ONE:
        case DATA_INTEGER_ONE:
        case DATA_ZERO_INTEGER:
        case DATA_ONE_INTEGER:
        case DATA_INTEGER_INTEGER:
        {
            i = getInteger(ii);
            j = getInteger(jj);

            if ( ( i >= 0 ) && ( j >= 0 ) && ( i >= j ) )
            {
                result = ResultOne();

                /*
                   Recall: perm(0,0) = 1
                           perm(i,0) = i!
                           perm(i,j) = i!/j! (assuming j > i)

                   BUT: overflows could happen here easily, so
                        need to use multiplication in GenVars to
                        ensure overflows are handled correctly
                        (ie. by conversion to double).
                */

                if ( i > j )
                {
                    for ( k = j+1 ; k <= i ; k++ )
                    {
                        result = mulGenVar(result,ResultInteger(k));
                    }
                }
            }

            else
            {
                result = ResultError(GVAR_GSL_EDOM);
            }

            break;
        }

        case DATA_ERROR_REAL:
        case DATA_ERROR_POS_INFTY:
        case DATA_ERROR_NEG_INFTY:
        case DATA_ERROR_FIN_INDET:
        case DATA_ERROR_INF_INDET:
        {
            result = combineErrors(ii,ResultError(GVAR_GSL_EDOM));

            break;
        }

        case DATA_REAL_ERROR:
        case DATA_POS_INFTY_ERROR:
        case DATA_NEG_INFTY_ERROR:
        case DATA_FIN_INDET_ERROR:
        case DATA_INF_INDET_ERROR:
        {
            result = combineErrors(ResultError(GVAR_GSL_EDOM),jj);

            break;
        }

        case DATA_UNCALCED_REAL:
        case DATA_UNCALCED_POS_INFTY:
        case DATA_UNCALCED_NEG_INFTY:
        case DATA_UNCALCED_FIN_INDET:
        case DATA_UNCALCED_INF_INDET:
        case DATA_REAL_UNCALCED:
        case DATA_POS_INFTY_UNCALCED:
        case DATA_NEG_INFTY_UNCALCED:
        case DATA_FIN_INDET_UNCALCED:
        case DATA_INF_INDET_UNCALCED:
        case DATA_REAL_ZERO:
        case DATA_POS_INFTY_ZERO:
        case DATA_NEG_INFTY_ZERO:
        case DATA_FIN_INDET_ZERO:
        case DATA_INF_INDET_ZERO:
        case DATA_REAL_ONE:
        case DATA_POS_INFTY_ONE:
        case DATA_NEG_INFTY_ONE:
        case DATA_FIN_INDET_ONE:
        case DATA_INF_INDET_ONE:
        case DATA_REAL_INTEGER:
        case DATA_POS_INFTY_INTEGER:
        case DATA_NEG_INFTY_INTEGER:
        case DATA_FIN_INDET_INTEGER:
        case DATA_INF_INDET_INTEGER:
        case DATA_ZERO_REAL:
        case DATA_ONE_REAL:
        case DATA_INTEGER_REAL:
        case DATA_REAL_REAL:
        case DATA_POS_INFTY_REAL:
        case DATA_NEG_INFTY_REAL:
        case DATA_FIN_INDET_REAL:
        case DATA_INF_INDET_REAL:
        case DATA_ZERO_POS_INFTY:
        case DATA_ONE_POS_INFTY:
        case DATA_INTEGER_POS_INFTY:
        case DATA_REAL_POS_INFTY:
        case DATA_POS_INFTY_POS_INFTY:
        case DATA_NEG_INFTY_POS_INFTY:
        case DATA_FIN_INDET_POS_INFTY:
        case DATA_INF_INDET_POS_INFTY:
        case DATA_ZERO_NEG_INFTY:
        case DATA_ONE_NEG_INFTY:
        case DATA_INTEGER_NEG_INFTY:
        case DATA_REAL_NEG_INFTY:
        case DATA_POS_INFTY_NEG_INFTY:
        case DATA_NEG_INFTY_NEG_INFTY:
        case DATA_FIN_INDET_NEG_INFTY:
        case DATA_INF_INDET_NEG_INFTY:
        case DATA_ZERO_FIN_INDET:
        case DATA_ONE_FIN_INDET:
        case DATA_INTEGER_FIN_INDET:
        case DATA_REAL_FIN_INDET:
        case DATA_POS_INFTY_FIN_INDET:
        case DATA_NEG_INFTY_FIN_INDET:
        case DATA_FIN_INDET_FIN_INDET:
        case DATA_INF_INDET_FIN_INDET:
        case DATA_ZERO_INF_INDET:
        case DATA_ONE_INF_INDET:
        case DATA_INTEGER_INF_INDET:
        case DATA_REAL_INF_INDET:
        case DATA_POS_INFTY_INF_INDET:
        case DATA_NEG_INFTY_INF_INDET:
        case DATA_FIN_INDET_INF_INDET:
        case DATA_INF_INDET_INF_INDET:
        {
            result = ResultError(GVAR_GSL_EDOM);

            break;
        }

        case DATA_ERROR_ERROR:
        default:
        {
            result = combineErrors(ii,jj);

            break;
        }
    }

    return result;
}

GenVar Gen_comb(GenVar ii, GenVar jj)
{
    GenVar result;
    GVAR_INT i;
    GVAR_INT n;
    GVAR_INT r;

    switch ( combineDataDescr(ii,jj) )
    {
        case DATA_ERROR_ZERO:
        case DATA_ERROR_ONE:
        case DATA_ERROR_INTEGER:
        case DATA_ERROR_UNCALCED:
        case DATA_UNCALCED_ZERO:
        case DATA_UNCALCED_ONE:
        case DATA_UNCALCED_INTEGER:
        case DATA_UNCALCED_UNCALCED:
        {
            result = ii;

            break;
        }

        case DATA_ZERO_ERROR:
        case DATA_ONE_ERROR:
        case DATA_INTEGER_ERROR:
        case DATA_UNCALCED_ERROR:
        case DATA_ZERO_UNCALCED:
        case DATA_ONE_UNCALCED:
        case DATA_INTEGER_UNCALCED:
        {
            result = jj;

            break;
        }

        case DATA_ZERO_ZERO:
        case DATA_ONE_ZERO:
        case DATA_INTEGER_ZERO:
        case DATA_ZERO_ONE:
        case DATA_ONE_ONE:
        case DATA_INTEGER_ONE:
        case DATA_ZERO_INTEGER:
        case DATA_ONE_INTEGER:
        case DATA_INTEGER_INTEGER:
        {
            n = getInteger(ii);
            r = getInteger(jj);

            if ( ( n >= 0 ) && ( r >= 0 ) && ( n >= r ) )
            {
                /*
                   Recall: comb(0,0) = 1
                           comb(n,0) = 1
                           comb(n,r) = n!/(r!(n-r)!) (assuming r > n)

                   BUT: overflows could happen here easily, so
                        need to use multiplication in GenVars to
                        ensure overflows are handled correctly
                        (ie. by conversion to double).
                */

                result = ResultOne();

                if ( ( r > 0 ) && ( n != r ) )
                {
                    /*
                       Want r as large as possible
                    */

                    if ( (n-r) > r )
                    {
                        r = n-r;
                    }

                    for ( i = r+1 ; i <= n ; i++ )
                    {
                        result = mulGenVar(result,ResultInteger(i));
                        result = divGenVar(result,ResultInteger(i-r));
                    }
                }
            }

            else
            {
                result = ResultError(GVAR_GSL_EDOM);
            }

            break;
        }

        case DATA_ERROR_REAL:
        case DATA_ERROR_POS_INFTY:
        case DATA_ERROR_NEG_INFTY:
        case DATA_ERROR_FIN_INDET:
        case DATA_ERROR_INF_INDET:
        {
            result = combineErrors(ii,ResultError(GVAR_GSL_EDOM));

            break;
        }

        case DATA_REAL_ERROR:
        case DATA_POS_INFTY_ERROR:
        case DATA_NEG_INFTY_ERROR:
        case DATA_FIN_INDET_ERROR:
        case DATA_INF_INDET_ERROR:
        {
            result = combineErrors(ResultError(GVAR_GSL_EDOM),jj);

            break;
        }

        case DATA_UNCALCED_REAL:
        case DATA_UNCALCED_POS_INFTY:
        case DATA_UNCALCED_NEG_INFTY:
        case DATA_UNCALCED_FIN_INDET:
        case DATA_UNCALCED_INF_INDET:
        case DATA_REAL_UNCALCED:
        case DATA_POS_INFTY_UNCALCED:
        case DATA_NEG_INFTY_UNCALCED:
        case DATA_FIN_INDET_UNCALCED:
        case DATA_INF_INDET_UNCALCED:
        case DATA_REAL_ZERO:
        case DATA_POS_INFTY_ZERO:
        case DATA_NEG_INFTY_ZERO:
        case DATA_FIN_INDET_ZERO:
        case DATA_INF_INDET_ZERO:
        case DATA_REAL_ONE:
        case DATA_POS_INFTY_ONE:
        case DATA_NEG_INFTY_ONE:
        case DATA_FIN_INDET_ONE:
        case DATA_INF_INDET_ONE:
        case DATA_REAL_INTEGER:
        case DATA_POS_INFTY_INTEGER:
        case DATA_NEG_INFTY_INTEGER:
        case DATA_FIN_INDET_INTEGER:
        case DATA_INF_INDET_INTEGER:
        case DATA_ZERO_REAL:
        case DATA_ONE_REAL:
        case DATA_INTEGER_REAL:
        case DATA_REAL_REAL:
        case DATA_POS_INFTY_REAL:
        case DATA_NEG_INFTY_REAL:
        case DATA_FIN_INDET_REAL:
        case DATA_INF_INDET_REAL:
        case DATA_ZERO_POS_INFTY:
        case DATA_ONE_POS_INFTY:
        case DATA_INTEGER_POS_INFTY:
        case DATA_REAL_POS_INFTY:
        case DATA_POS_INFTY_POS_INFTY:
        case DATA_NEG_INFTY_POS_INFTY:
        case DATA_FIN_INDET_POS_INFTY:
        case DATA_INF_INDET_POS_INFTY:
        case DATA_ZERO_NEG_INFTY:
        case DATA_ONE_NEG_INFTY:
        case DATA_INTEGER_NEG_INFTY:
        case DATA_REAL_NEG_INFTY:
        case DATA_POS_INFTY_NEG_INFTY:
        case DATA_NEG_INFTY_NEG_INFTY:
        case DATA_FIN_INDET_NEG_INFTY:
        case DATA_INF_INDET_NEG_INFTY:
        case DATA_ZERO_FIN_INDET:
        case DATA_ONE_FIN_INDET:
        case DATA_INTEGER_FIN_INDET:
        case DATA_REAL_FIN_INDET:
        case DATA_POS_INFTY_FIN_INDET:
        case DATA_NEG_INFTY_FIN_INDET:
        case DATA_FIN_INDET_FIN_INDET:
        case DATA_INF_INDET_FIN_INDET:
        case DATA_ZERO_INF_INDET:
        case DATA_ONE_INF_INDET:
        case DATA_INTEGER_INF_INDET:
        case DATA_REAL_INF_INDET:
        case DATA_POS_INFTY_INF_INDET:
        case DATA_NEG_INFTY_INF_INDET:
        case DATA_FIN_INDET_INF_INDET:
        case DATA_INF_INDET_INF_INDET:
        {
            result = ResultError(GVAR_GSL_EDOM);

            break;
        }

        case DATA_ERROR_ERROR:
        default:
        {
            result = combineErrors(ii,jj);

            break;
        }
    }

    return result;
}

GenVar Gen_fact(GenVar ii)
{
    GenVar result;
    GVAR_INT i;
    GVAR_INT j;

    switch ( ii.DataType )
    {
        case DATA_IS_ONE:
        case DATA_IS_ZERO:
        case DATA_IS_INTEGER:
        {
            i = getInteger(ii);

            if ( i >= 0 )
            {
                /*
                   Beware of overflow.
                */

                result = ResultOne();

                if ( i > 1 )
                {
                    for ( j = 2 ; j <= i ; j++ )
                    {
                         result = mulGenVar(result,ResultInteger(j));
                    }
                }
            }

            else
            {
                result = ResultError(GVAR_GSL_EDOM);
            }

            break;
        }

        case DATA_IS_REAL:
        case DATA_IS_POS_INFTY:
        case DATA_IS_NEG_INFTY:
        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        {
            result = ResultError(GVAR_GSL_EDOM);

            break;
        }

        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = ii;

            break;
        }
    }

    return result;
}

GenVar Gen_pow(GenVar xx, GenVar yy)
{
    GenVar result;
    unsigned GVAR_INT i,j;
    GVAR_INT k;

    /*
       Notes: - pow(0,0)        = 1
              - pow(+-inf,0)    = 1
              - pow(0,+-inf)    = 0
              - pow(x,-inf)     = 0
              - pow(+-inf,-inf) = 0
    */

    switch ( combineDataDescr(xx,yy) )
    {
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
        case DATA_ONE_ONE:
        case DATA_INTEGER_ONE:
        case DATA_REAL_ONE:
        case DATA_POS_INFTY_ONE:
        case DATA_NEG_INFTY_ONE:
        case DATA_FIN_INDET_ONE:
        case DATA_INF_INDET_ONE:
        case DATA_INF_INDET_INTEGER:
        case DATA_POS_INFTY_POS_INFTY:
        {
            result = xx;

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
        case DATA_UNCALCED_ERROR:
        case DATA_ZERO_UNCALCED:
        case DATA_ONE_UNCALCED:
        case DATA_INTEGER_UNCALCED:
        case DATA_REAL_UNCALCED:
        case DATA_POS_INFTY_UNCALCED:
        case DATA_NEG_INFTY_UNCALCED:
        case DATA_FIN_INDET_UNCALCED:
        case DATA_INF_INDET_UNCALCED:
        {
            result = yy;

            break;
        }

        case DATA_POS_INFTY_INTEGER:
        case DATA_POS_INFTY_REAL:
        {
            if ( getReal(yy) > 0 )
            {
                result = ResultPosInfnty();
            }

            else
            {
                result = ResultZero();
            }

            break;
        }

        case DATA_ZERO_ZERO:
        case DATA_ONE_ZERO:
        case DATA_INTEGER_ZERO:
        case DATA_REAL_ZERO:
        case DATA_POS_INFTY_ZERO:
        case DATA_NEG_INFTY_ZERO:
        case DATA_FIN_INDET_ZERO:
        case DATA_INF_INDET_ZERO:
        case DATA_UNCALCED_ZERO:
        case DATA_ONE_INTEGER:
        case DATA_ONE_REAL:
        case DATA_ONE_POS_INFTY:
        case DATA_ONE_NEG_INFTY:
        case DATA_ONE_FIN_INDET:
        case DATA_ONE_INF_INDET:
        {
            result = ResultOne();

            break;
        }

        case DATA_ZERO_ONE:
        case DATA_ZERO_POS_INFTY:
        case DATA_POS_INFTY_NEG_INFTY:
        case DATA_NEG_INFTY_NEG_INFTY:
        {
            result = ResultZero();

            break;
        }

        case DATA_ZERO_INTEGER:
        case DATA_ZERO_REAL:
        {
            if ( getReal(yy) > 0 )
            {
                result = ResultZero();
            }

            else
            {
                result = ResultPosInfnty();
            }

            break;
        }

        case DATA_ZERO_NEG_INFTY:
        {
            result = ResultPosInfnty();

            break;
        }

        case DATA_INTEGER_INTEGER:
        case DATA_REAL_INTEGER:
        {
            /*
               For finite integers, leave the work to mulGenVar - this
               will keep type where possible, track overflows and
               generally do the hard, annoying work.
            */

            result = ResultOne();

            k = getInteger(yy);

            if ( k > 0 )
            {
                j = (unsigned GVAR_INT) k;

                for ( i = 1 ; i <= j ; i++ )
                {
                     result = mulGenVar(result,xx);
                }
            }

            else
            {
                j = (unsigned GVAR_INT) -k;

                for ( i = 1 ; i <= j ; i++ )
                {
                    result = mulGenVar(result,xx);
                }

                result = invGenVar(result);
            }

            break;
        }

        case DATA_INTEGER_REAL:
        case DATA_REAL_REAL:
        {
            if ( getReal(xx) > 0 )
            {
                result = ResultReal(pow(getReal(xx),getReal(yy)));
            }

            else
            {
                result = ResultError(GVAR_GSL_EDOM);
            }

            break;
        }

        case DATA_NEG_INFTY_INTEGER:
        {
            if ( getInteger(yy)%2 )
            {
                result = ResultNegInfnty();
            }

            else
            {
                result = ResultPosInfnty();
            }

            break;
        }

        case DATA_NEG_INFTY_REAL:
        case DATA_NEG_INFTY_FIN_INDET:
        case DATA_NEG_INFTY_INF_INDET:
        {
            result = ResultError(GVAR_GSL_EDOM);

            break;
        }

        case DATA_FIN_INDET_INTEGER:
        {
            if ( getInteger(yy) > 0 )
            {
                result = ResultFinIndet();
            }

            else
            {
                result = ResultInfIndet();
            }

            break;
        }

        case DATA_INF_INDET_REAL:
        case DATA_FIN_INDET_REAL:
        case DATA_FIN_INDET_FIN_INDET:
        case DATA_INF_INDET_FIN_INDET:
        case DATA_FIN_INDET_INF_INDET:
        case DATA_INF_INDET_INF_INDET:
        {
            result = ResultError(GVAR_GSL_EDOM);

            break;
        }

        case DATA_NEG_INFTY_POS_INFTY:
        case DATA_FIN_INDET_POS_INFTY:
        case DATA_INF_INDET_POS_INFTY:
        case DATA_INTEGER_INF_INDET:
        case DATA_REAL_INF_INDET:
        case DATA_POS_INFTY_FIN_INDET:
        case DATA_POS_INFTY_INF_INDET:
        case DATA_ZERO_FIN_INDET:
        case DATA_ZERO_INF_INDET:
        case DATA_FIN_INDET_NEG_INFTY:
        case DATA_INF_INDET_NEG_INFTY:
        {
            result = ResultInfIndet();

            break;
        }

        case DATA_INTEGER_POS_INFTY:
        {
            if ( getInteger(xx) > 0 )
            {
                /* we know that xx != 1 */

                result = ResultPosInfnty();

                break;
            }

            else
            {
                if ( getInteger(xx) == -1 )
                {
                    result = ResultFinIndet();
                }

                else
                {
                    result = ResultInfIndet();
                }
            }

            break;
        }

        case DATA_REAL_POS_INFTY:
        {
            if ( getReal(xx) > 0 )
            {
                if ( getReal(xx) < 1 )
                {
                    result = ResultZero();
                }

                else
                {
                    result = ResultPosInfnty();
                }
            }

            else
            {
                if ( getReal(xx) > -1 )
                {
                    result = ResultZero();
                }

                else
                {
                    result = ResultInfIndet();
                }
            }

            break;
        }

        case DATA_INTEGER_NEG_INFTY:
        {
            if ( getInteger(xx) == -1 )
            {
                result = ResultFinIndet();
            }

            else
            {
                result = ResultZero();
            }

            break;
        }

        case DATA_REAL_NEG_INFTY:
        {
            if ( ( getReal(xx) > -1 ) && ( getReal(xx) < 1 ) )
            {
                if ( getReal(xx) < 0 )
                {
                    result = ResultInfIndet();
                }

                else
                {
                    result = ResultPosInfnty();
                }
            }

            else
            {
                result = ResultZero();
            }

            break;
        }

        case DATA_INTEGER_FIN_INDET:
        case DATA_REAL_FIN_INDET:
        {
            result = ResultFinIndet();

            break;
        }

        case DATA_ERROR_ERROR:
        default:
        {
            result = combineErrors(xx,yy);

            break;
        }
    }

    return result;
}

GenVar Gen_sqrt(GenVar xx)
{
    GenVar result;

    switch ( xx.DataType )
    {
        case DATA_IS_ZERO:
        case DATA_IS_ONE:
        {
            result = xx;

            break;
        }

        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        {
            if ( getReal(xx) > 0 )
            {
                result = ResultReal(sqrt(getReal(xx)));
            }

            else
            {
                result = ResultError(GVAR_GSL_EDOM);
            }

            break;
        }

        case DATA_IS_POS_INFTY:
        {
            result = xx;

            break;
        }

        case DATA_IS_NEG_INFTY:
        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        {
            result = ResultError(GVAR_GSL_EDOM);

            break;
        }

        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_sin(GenVar xx)
{
    GenVar result;

    switch ( xx.DataType )
    {
        case DATA_IS_ZERO:
        {
            result = ResultZero();

            break;
        }

        case DATA_IS_ONE:
        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        {
            result = ResultReal(sin(getReal(xx)));

            break;
        }

        case DATA_IS_POS_INFTY:
        case DATA_IS_NEG_INFTY:
        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        {
            result = ResultFinIndet();

            break;
        }

        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_cos(GenVar xx)
{
    GenVar result;

    switch ( xx.DataType )
    {
        case DATA_IS_ZERO:
        {
            result = ResultOne();

            break;
        }

        case DATA_IS_ONE:
        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        {
            result = ResultReal(cos(getReal(xx)));

            break;
        }

        case DATA_IS_POS_INFTY:
        case DATA_IS_NEG_INFTY:
        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        {
            result = ResultFinIndet();

            break;
        }

        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_tan(GenVar xx)
{
    GenVar result;

    switch ( xx.DataType )
    {
        case DATA_IS_ZERO:
        {
            result = ResultZero();

            break;
        }

        case DATA_IS_ONE:
        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        {
            result = ResultReal(tan(getReal(xx)));

            if ( !isFiniteReal(result) )
            {
                result = ResultInfIndet();
            }

            break;
        }

        case DATA_IS_POS_INFTY:
        case DATA_IS_NEG_INFTY:
        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        {
            result = ResultInfIndet();

            break;
        }

        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_cosec(GenVar xx)
{
    GenVar result;

    /*
       Can't use shorthand, as it stuffs up the infinite indeterminants.
    */

    switch ( xx.DataType )
    {
        case DATA_IS_ZERO:
        {
            result = ResultInfIndet();

            break;
        }

        case DATA_IS_ONE:
        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        {
            result = ResultReal(sin(getReal(xx)));

            if ( result.DataType == DATA_IS_ZERO )
            {
                result = ResultInfIndet();
            }

            else
            {
                result = invGenVar(result);
            }

            break;
        }

        case DATA_IS_POS_INFTY:
        case DATA_IS_NEG_INFTY:
        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        {
            result = ResultInfIndet();

            break;
        }

        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_sec(GenVar xx)
{
    GenVar result;

    /*
       Can't use shorthand, as it stuffs up the infinite indeterminants.
    */

    switch ( xx.DataType )
    {
        case DATA_IS_ZERO:
        {
            result = ResultOne();

            break;
        }

        case DATA_IS_ONE:
        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        {
            result = ResultReal(cos(getReal(xx)));

            if ( result.DataType == DATA_IS_ZERO )
            {
                result = ResultInfIndet();
            }

            else
            {
                result = invGenVar(result);
            }

            break;
        }

        case DATA_IS_POS_INFTY:
        case DATA_IS_NEG_INFTY:
        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        {
            result = ResultInfIndet();

            break;
        }

        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_cot(GenVar xx)
{
    GenVar result;

    /*
       Can't use shorthand, as it stuffs up the infinite indeterminants.
       (and zeroes in this case).
    */

    switch ( xx.DataType )
    {
        case DATA_IS_ZERO:
        {
            result = ResultInfIndet();

            break;
        }

        case DATA_IS_ONE:
        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        {
            result = ResultReal(tan(getReal(xx)));

            if ( !isFiniteReal(result) )
            {
                result = ResultZero();
            }

            else
            {
                if ( result.DataType == DATA_IS_ZERO )
                {
                    result = ResultInfIndet();
                }

                else
                {
                    result = invGenVar(result);
                }
            }

            break;
        }

        case DATA_IS_POS_INFTY:
        case DATA_IS_NEG_INFTY:
        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        {
            result = ResultInfIndet();

            break;
        }

        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_asin(GenVar xx)
{
    GenVar result;

    switch ( xx.DataType )
    {
        case DATA_IS_ZERO:
        {
            result = ResultZero();

            break;
        }

        case DATA_IS_ONE:
        {
            result = ResultReal(M_PI/2);

            break;
        }

        case DATA_IS_INTEGER:
        {
            if ( getInteger(xx) == -1 )
            {
                result = ResultReal(-M_PI/2);
            }

            else
            {
                result = ResultError(GVAR_GSL_EDOM);
            }

            break;
        }

        case DATA_IS_REAL:
        {
            if ( ( getReal(xx) < 1 ) && ( getReal(xx) > -1 ) )
            {
                result = ResultReal(asin(getReal(xx)));
            }

            else
            {
                result = ResultError(GVAR_GSL_EDOM);
            }

            break;
        }

        case DATA_IS_POS_INFTY:
        case DATA_IS_NEG_INFTY:
        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        {
            result = ResultError(GVAR_GSL_EDOM);

            break;
        }

        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_acos(GenVar xx)
{
    GenVar result;

    switch ( xx.DataType )
    {
        case DATA_IS_ZERO:
        {
            result = ResultReal(M_PI/2);

            break;
        }

        case DATA_IS_ONE:
        {
            result = ResultZero();

            break;
        }

        case DATA_IS_INTEGER:
        {
            if ( getInteger(xx) == -1 )
            {
                result = ResultReal(M_PI);
            }

            else
            {
                result = ResultError(GVAR_GSL_EDOM);
            }

            break;
        }

        case DATA_IS_REAL:
        {
            if ( ( getReal(xx) < 1 ) && ( getReal(xx) > -1 ) )
            {
                result = ResultReal(acos(getReal(xx)));
            }

            else
            {
                result = ResultError(GVAR_GSL_EDOM);
            }

            break;
        }

        case DATA_IS_POS_INFTY:
        case DATA_IS_NEG_INFTY:
        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        {
            result = ResultError(GVAR_GSL_EDOM);

            break;
        }

        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_atan(GenVar xx)
{
    GenVar result;

    switch ( xx.DataType )
    {
        case DATA_IS_ZERO:
        {
            result = ResultZero();

            break;
        }

        case DATA_IS_ONE:
        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        {
            result = ResultReal(atan(getReal(xx)));

            break;
        }

        case DATA_IS_POS_INFTY:
        {
            result = ResultReal(M_PI/2);

            break;
        }

        case DATA_IS_NEG_INFTY:
        {
            result = ResultReal(-M_PI/2);

            break;
        }

        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        {
            result = ResultFinIndet();

            break;
        }

        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_acosec(GenVar xx)
{
    return Gen_asin(invGenVar(xx));
}

GenVar Gen_asec(GenVar xx)
{
    return Gen_acos(invGenVar(xx));
}

GenVar Gen_acot(GenVar xx)
{
    GenVar result;

    switch ( xx.DataType )
    {
        case DATA_IS_ZERO:
        {
            result = ResultFinIndet();

            break;
        }

        case DATA_IS_ONE:
        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        {
            result = ResultReal(atan(1/getReal(xx)));

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
            result = ResultFinIndet();

            break;
        }

        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_sinh(GenVar xx)
{
    GenVar result;

    switch ( xx.DataType )
    {
        case DATA_IS_ZERO:
        {
            result = ResultZero();

            break;
        }

        case DATA_IS_ONE:
        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        {
            result = ResultReal(sinh(getReal(xx)));

            break;
        }

        case DATA_IS_POS_INFTY:
        case DATA_IS_NEG_INFTY:
        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_cosh(GenVar xx)
{
    GenVar result;

    switch ( xx.DataType )
    {
        case DATA_IS_ZERO:
        {
            result = ResultOne();

            break;
        }

        case DATA_IS_ONE:
        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        {
            result = ResultReal(cosh(getReal(xx)));

            break;
        }

        case DATA_IS_POS_INFTY:
        case DATA_IS_NEG_INFTY:
        {
            result = ResultPosInfnty();

            break;
        }

        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_tanh(GenVar xx)
{
    GenVar result;

    switch ( xx.DataType )
    {
        case DATA_IS_ZERO:
        {
            result = ResultZero();

            break;
        }

        case DATA_IS_ONE:
        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        {
            result = ResultReal(tanh(getReal(xx)));

            break;
        }

        case DATA_IS_POS_INFTY:
        {
            result = ResultOne();

            break;
        }

        case DATA_IS_NEG_INFTY:
        {
            result = ResultInteger(-1);

            break;
        }

        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        {
            result = ResultFinIndet();

            break;
        }

        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_cosech(GenVar xx)
{
    GenVar result;

    switch ( xx.DataType )
    {
        case DATA_IS_ZERO:
        {
            result = ResultInfIndet();

            break;
        }

        case DATA_IS_ONE:
        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        {
            result = invGenVar(ResultReal(sinh(getReal(xx))));

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

        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}


GenVar Gen_sech(GenVar xx)
{
    return invGenVar(Gen_cosh(xx));
}

GenVar Gen_coth(GenVar xx)
{
    GenVar result;

    switch ( xx.DataType )
    {
        case DATA_IS_ZERO:
        {
            result = ResultInfIndet();

            break;
        }

        case DATA_IS_ONE:
        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        {
            result = invGenVar(ResultReal(tanh(getReal(xx))));

            break;
        }

        case DATA_IS_POS_INFTY:
        {
            result = ResultOne();

            break;
        }

        case DATA_IS_NEG_INFTY:
        {
            result = ResultInteger(-1);

            break;
        }

        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        {
            result = ResultInfIndet();

            break;
        }

        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_asinh(GenVar xx)
{
    GenVar result;

    switch ( xx.DataType )
    {
        case DATA_IS_ZERO:
        {
            result = ResultZero();

            break;
        }

        case DATA_IS_ONE:
        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        {
            result = ResultReal(gsl_asinh(getReal(xx)));

            break;
        }

        case DATA_IS_POS_INFTY:
        case DATA_IS_NEG_INFTY:
        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_acosh(GenVar xx)
{
    GenVar result;

    switch ( xx.DataType )
    {
        case DATA_IS_ZERO:
        {
            result = ResultError(GVAR_GSL_EDOM);

            break;
        }

        case DATA_IS_ONE:
        {
            result = ResultZero();

            break;
        }

        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        {
            if ( getReal(xx) > 1 )
            {
                result = ResultReal(gsl_acosh(getReal(xx)));
            }

            else
            {
                result = ResultError(GVAR_GSL_EDOM);
            }

            break;
        }

        case DATA_IS_POS_INFTY:
        {
            result = ResultPosInfnty();

            break;
        }

        case DATA_IS_NEG_INFTY:
        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        {
            result = ResultError(GVAR_GSL_EDOM);

            break;
        }

        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_atanh(GenVar xx)
{
    GenVar result;
    double temp;

    switch ( xx.DataType )
    {
        case DATA_IS_ZERO:
        {
            result = ResultZero();

            break;
        }

        case DATA_IS_ONE:
        {
            result = ResultPosInfnty();

            break;
        }

        case DATA_IS_INTEGER:
        {
            if ( getInteger(xx) == -1 )
            {
                result = ResultNegInfnty();
            }

            else
            {
                result = ResultError(GVAR_GSL_EDOM);
            }

            break;
        }

        case DATA_IS_REAL:
        {
            temp = getReal(xx);

            if ( ( temp < 1 ) && ( temp > -1 ) )
            {
                result = ResultReal(gsl_atanh(temp));
            }

            else
            {
                result = ResultError(GVAR_GSL_EDOM);
            }

            break;
        }

        case DATA_IS_POS_INFTY:
        case DATA_IS_NEG_INFTY:
        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        {
            result = ResultError(GVAR_GSL_EDOM);

            break;
        }

        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_acosech(GenVar xx)
{
    GenVar result;
    double temp;

    switch ( xx.DataType )
    {
        case DATA_IS_ZERO:
        {
            result = ResultInfIndet();

            break;
        }

        case DATA_IS_ONE:
        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        {
            temp = getReal(xx);

            result = invGenVar(ResultReal(gsl_asinh(temp)));

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

        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_asech(GenVar xx)
{
    GenVar result;
    double temp;

    switch ( xx.DataType )
    {
        case DATA_IS_ZERO:
        {
            result = ResultError(GVAR_GSL_EDOM);

            break;
        }

        case DATA_IS_ONE:
        {
            result = ResultPosInfnty();

            break;
        }

        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        {
            temp = getReal(xx);

            if ( temp > 1 )
            {
                result = invGenVar(ResultReal(gsl_acosh(temp)));
            }

            else
            {
                result = ResultError(GVAR_GSL_EDOM);
            }

            break;
        }

        case DATA_IS_POS_INFTY:
        case DATA_IS_NEG_INFTY:
        {
            result = ResultError(GVAR_GSL_EDOM);

            break;
        }

        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        {
            result = ResultError(GVAR_GSL_EDOM);

            break;
        }

        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_acoth(GenVar xx)
{
    return Gen_atanh(invGenVar(xx));
}

GenVar Gen_sinc(GenVar xx)
{
    GenVar result;

    switch ( xx.DataType )
    {
        case DATA_IS_ZERO:
        {
            result = ResultOne();

            break;
        }

        case DATA_IS_ONE:
        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        {
            result = ResultReal(gsl_sf_sinc(getReal(xx)));

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
            result = ResultFinIndet();

            break;
        }

        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_gd(GenVar xx)
{
    GenVar result;

    switch ( xx.DataType )
    {
        case DATA_IS_ZERO:
        {
            result = ResultZero();

            break;
        }

        case DATA_IS_ONE:
        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        {
            result = ResultReal(((2*atan(exp(getReal(xx))))-(M_PI/2)));

            break;
        }

        case DATA_IS_POS_INFTY:
        {
            result = ResultReal(M_PI/2);

            break;
        }

        case DATA_IS_NEG_INFTY:
        {
            result = ResultReal(-M_PI/2);

            break;
        }

        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        {
            result = ResultFinIndet();

            break;
        }

        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_agd(GenVar xx)
{
    GenVar result;

    switch ( xx.DataType )
    {
        case DATA_IS_ZERO:
        case DATA_IS_ONE:
        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        {
            if ( ( getReal(xx) <= M_PI/2 ) && ( getReal(xx) >= -M_PI/2 ) )
            {
                if ( getReal(xx) == -M_PI/2 )
                {
                    result = ResultPosInfnty();
                }

                else
                {
                    if ( getReal(xx) == M_PI/2 )
                    {
                        result = ResultNegInfnty();
                    }

                    else
                    {
                        result = ResultReal(log(tan((getReal(xx)/2)+(M_PI/4))));
                    }
                }
            }

            else
            {
                result = ResultError(GVAR_GSL_EDOM);
            }

            break;
        }

        case DATA_IS_POS_INFTY:
        case DATA_IS_NEG_INFTY:
        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        {
            result = ResultError(GVAR_GSL_EDOM);

            break;
        }

        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_exp(GenVar xx)
{
    GenVar result;
    gsl_sf_result gres;
    int ires;

    switch ( xx.DataType )
    {
        case DATA_IS_ZERO:
        {
            result = ResultOne();

            break;
        }

        case DATA_IS_ONE:
        {
            result = ResultReal(M_E);

            break;
        }

        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        {
            switch ( ( ires = gsl_sf_exp_e(getReal(xx),&gres) ) )
            {
                case GSL_SUCCESS:
                {
                    result = ResultReal(gres.val);

                    break;
                }

                case GSL_EUNDRFLW:
                case GSL_EOVRFLW:
                case GSL_ERANGE:
                {
                    if ( getReal(xx) > 0 )
                    {
                        result = ResultPosInfnty();
                    }

                    else
                    {
                        result = ResultZero();
                    }

                    break;
                }

                default:
                {
                    result = ResultError(convGSLErr(ires));

                    break;
                }
            }

            break;
        }

        case DATA_IS_POS_INFTY:
        {
            result = ResultPosInfnty();

            break;
        }

        case DATA_IS_NEG_INFTY:
        {
            result = ResultZero();

            break;
        }

        case DATA_IS_FIN_INDET:
        {
            result = ResultFinIndet();

            break;
        }

        case DATA_IS_INF_INDET:
        {
            result = ResultInfIndet();

            break;
        }

        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_tenup(GenVar xx)
{
    GenVar result;

    switch ( xx.DataType )
    {
        case DATA_IS_ZERO:
        {
            result = ResultOne();

            break;
        }

        case DATA_IS_ONE:
        {
            result = ResultInteger(10);

            break;
        }

        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        {
            result = Gen_pow(ResultInteger(10),xx);

            break;
        }

        case DATA_IS_POS_INFTY:
        {
            result = ResultPosInfnty();

            break;
        }

        case DATA_IS_NEG_INFTY:
        {
            result = ResultZero();

            break;
        }

        case DATA_IS_FIN_INDET:
        {
            result = ResultFinIndet();

            break;
        }

        case DATA_IS_INF_INDET:
        {
            result = ResultInfIndet();

            break;
        }

        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_log(GenVar xx)
{
    GenVar result;
    gsl_sf_result gres;
    int ires;

    switch ( xx.DataType )
    {
        case DATA_IS_ZERO:
        {
            result = ResultNegInfnty();

            break;
        }

        case DATA_IS_ONE:
        {
            result = ResultZero();

            break;
        }

        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        {
            switch ( ( ires = gsl_sf_log_e(getReal(xx),&gres) ) )
            {
                case GSL_SUCCESS:
                {
                    result = ResultReal(gres.val);

                    break;
                }

                case GSL_EUNDRFLW:
                case GSL_EOVRFLW:
                case GSL_ERANGE:
                {
                    if ( getReal(xx) > 1 )
                    {
                        result = ResultPosInfnty();
                    }

                    else
                    {
                        result = ResultZero();
                    }

                    break;
                }

                default:
                {
                    result = ResultError(convGSLErr(ires));

                    break;
                }
            }

            break;
        }

        case DATA_IS_POS_INFTY:
        {
            result = ResultPosInfnty();

            break;
        }

        case DATA_IS_NEG_INFTY:
        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        {
            result = ResultError(GVAR_GSL_EDOM);

            break;
        }

        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_logX(GenVar xx)
{
    return Gen_log(Gen_pow(xx,ResultReal(0.4342944819032518)));
}

GenVar Gen_gamma(GenVar aa)
{
    GenVar result;
    gsl_sf_result gres;
    int ires;

    switch ( aa.DataType )
    {
        case DATA_IS_ZERO:
        {
            result = ResultInfIndet();

            break;
        }

        case DATA_IS_ONE:
        {
            result = ResultOne();

            break;
        }

        case DATA_IS_INTEGER:
        {
            if ( getInteger(aa) > 0 )
            {
                result = Gen_fact(ResultInteger(getInteger(aa)-1));
            }

            else
            {
                result = ResultInfIndet();
            }

            break;
        }

        case DATA_IS_REAL:
        {
            switch ( ( ires = gsl_sf_gamma_e(getReal(aa),&gres) ) )
            {
                case GSL_SUCCESS:
                {
                    result = ResultReal(gres.val);

                    break;
                }

                case GSL_EUNDRFLW:
                {
                    result = ResultZero();

                    break;
                }

                case GSL_EOVRFLW:
                case GSL_ERANGE:
                {
                    if ( getReal(aa) > 0 )
                    {
                        result = ResultPosInfnty();
                    }

                    else
                    {
                        result = ResultInfIndet();
                    }

                    break;
                }

                default:
                {
                    result = ResultError(convGSLErr(ires));

                    break;
                }
            }

            break;
        }

        case DATA_IS_POS_INFTY:
        {
            result = ResultPosInfnty();

            break;
        }

        case DATA_IS_NEG_INFTY:
        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        {
            result = ResultInfIndet();

            break;
        }

        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = aa;

            break;
        }
    }

    return result;
}

GenVar Gen_lngamma(GenVar aa)
{
    GenVar result;
    gsl_sf_result gres;
    int ires;

    switch ( aa.DataType )
    {
        case DATA_IS_ZERO:
        {
            result = ResultPosInfnty();

            break;
        }

        case DATA_IS_ONE:
        {
            result = ResultZero();

            break;
        }

        case DATA_IS_INTEGER:
        {
            if ( getInteger(aa) < 0 )
            {
                result = ResultPosInfnty();
            }

            else
            {
                switch ( ( ires = gsl_sf_lngamma_e(getReal(aa),&gres) ) )
                {
                    case GSL_SUCCESS:
                    {
                        result = ResultReal(gres.val);

                        break;
                    }

                    case GSL_EUNDRFLW:
                    {
                        result = ResultZero();

                        break;
                    }

                    case GSL_EOVRFLW:
                    case GSL_ERANGE:
                    {
                        result = ResultPosInfnty();

                        break;
                    }

                    default:
                    {
                        result = ResultError(convGSLErr(ires));

                        break;
                    }
                }
            }

            break;
        }

        case DATA_IS_REAL:
        {
            switch ( ( ires = gsl_sf_lngamma_e(getReal(aa),&gres) ) )
            {
                case GSL_SUCCESS:
                {
                    result = ResultReal(gres.val);

                    break;
                }

                case GSL_EUNDRFLW:
                {
                    result = ResultZero();

                    break;
                }

                case GSL_EOVRFLW:
                case GSL_ERANGE:
                {
                    /*
                       FIXME: there are two possible sources for overflow
                              here.  First, if |gamma(aa)| is sufficiently
                              close to zero then -inf will result.  The case
                              done here, though, assumes |gamma(aa)| is too
                              large, resulting in an overflow to +inf.
                              Need to deal with both possibilities.
                    */

                    result = ResultPosInfnty();

                    break;
                }

                default:
                {
                    result = ResultError(convGSLErr(ires));

                    break;
                }
            }

            break;
        }

        case DATA_IS_POS_INFTY:
        {
            result = ResultPosInfnty();

            break;
        }

        case DATA_IS_NEG_INFTY:
        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        {
            result = ResultInfIndet();

            break;
        }

        case DATA_IS_ERROR:
        case DATA_IS_UNCALCED:
        default:
        {
            result = aa;

            break;
        }
    }

    return result;
}

GenVar Gen_psi(GenVar xx)
{
    GenVar result;
    gsl_sf_result gres;
    int ires;

    switch ( xx.DataType )
    {
        case DATA_IS_ZERO:
        case DATA_IS_ONE:
        case DATA_IS_INTEGER:
        {
            if ( getInteger(xx) <= 0 )
            {
                /*
                   the log gamma has points of inflexion at negative
                   integers and 0, where it starts approaching +inf
                   and then continues down from +inf.  Overflow will
                   only happen near hear.
                */

                result = ResultNegInfnty();
            }

            else
            {
                switch ( ( ires = gsl_sf_psi_int_e(getInteger(xx),&gres) ) )
                {
                    case GSL_SUCCESS:
                    {
                        result = ResultReal(gres.val);

                        break;
                    }

                    case GSL_EUNDRFLW:
                    {
                        result = ResultZero();

                        break;
                    }

                    case GSL_EOVRFLW:
                    case GSL_ERANGE:
                    {
                        result = ResultPosInfnty();

                        break;
                    }

                    default:
                    {
                        result = ResultError(convGSLErr(ires));

                        break;
                    }
                }
            }

            break;
        }

        case DATA_IS_REAL:
        {
            switch ( ( ires = gsl_sf_psi_e(getReal(xx),&gres) ) )
            {
                case GSL_SUCCESS:
                {
                    result = ResultReal(gres.val);

                    break;
                }

                case GSL_EUNDRFLW:
                {
                    result = ResultZero();

                    break;
                }

                case GSL_EOVRFLW:
                case GSL_ERANGE:
                {
                    if ( getReal(xx) > 0 )
                    {
                        result = ResultPosInfnty();
                    }

                    else
                    {
                        /*
                           the log gamma has points of inflexion at negative
                           integers and 0, where it starts approaching +inf
                           and then continues down from +inf.  Overflow will
                           only happen near hear.
                        */

                        /*
                           FIXME: there are two possible sources for overflow
                               lngamma.  First, if |gamma(aa)| is sufficiently
                               close to zero then -inf will result.  The case
                               done here, though, assumes |gamma(aa)| is too
                               large, resulting in an overflow to +inf.
                               Need to deal with both possibilities.
                        */

                        result = ResultNegInfnty();
                    }

                    break;
                }

                default:
                {
                    result = ResultError(convGSLErr(ires));

                    break;
                }
            }

            break;
        }

        case DATA_IS_POS_INFTY:
        {
            result = ResultPosInfnty();

            break;
        }

        case DATA_IS_NEG_INFTY:
        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        {
            result = ResultInfIndet();

            break;
        }

        case DATA_IS_UNCALCED:
        case DATA_IS_ERROR:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_psi_n(GenVar ii, GenVar xx)
{
    GenVar result;
    gsl_sf_result gres;
    int ires;

    if ( isFiniteInteger(ii) && ( getInteger(ii) >= 0 ) )
    {
        switch ( xx.DataType )
        {
            case DATA_IS_ZERO:
            case DATA_IS_ONE:
            case DATA_IS_INTEGER:
            case DATA_IS_REAL:
            {
                if ( getReal(xx) <= 0 )
                {
                    result = ResultError(GVAR_GSL_EDOM);
                }

                else
                {
                    switch ( ( ires = gsl_sf_psi_n_e(getInteger(ii),getReal(xx),&gres) ) )
                    {
                        case GSL_SUCCESS:
                        {
                            result = ResultReal(gres.val);

                            break;
                        }

                        default:
                        {
                            result = ResultError(convGSLErr(ires));

                            break;
                        }
                    }
                }

                break;
            }

            case DATA_IS_POS_INFTY:
            {
                if ( getInteger(ii) == 0 )
                {
                    result = ResultPosInfnty();
                }

                else
                {
                    result = ResultError(GVAR_GSL_EDOM);
                }

                break;
            }

            case DATA_IS_NEG_INFTY:
            case DATA_IS_FIN_INDET:
            case DATA_IS_INF_INDET:
            {
                result = ResultError(GVAR_GSL_EDOM);

                break;
            }

            case DATA_IS_UNCALCED:
            case DATA_IS_ERROR:
            default:
            {
                result = xx;

                break;
            }
        }
    }

    else
    {
        result = ResultError(GVAR_GSL_EDOM);
    }

    return result;
}

GenVar Gen_gami(GenVar a, GenVar x)
{
    return subGenVar(Gen_gamma(a),Gen_gamic(a,x));
}

GenVar Gen_gamic(GenVar a, GenVar x)
{
    GenVar result;
    gsl_sf_result gres;
    int ires;

    /*
       FIXME: this does not account for all limits and stuff.  Nor does
              it correctly allow overflow to goto infinity.
    */

    switch ( combineDataDescr(a,x) )
    {
        case DATA_ERROR_ZERO:
        case DATA_ERROR_ONE:
        case DATA_ERROR_POS_INFTY:
        case DATA_ERROR_UNCALCED:
        case DATA_UNCALCED_ZERO:
        case DATA_UNCALCED_ONE:
        case DATA_UNCALCED_POS_INFTY:
        case DATA_UNCALCED_UNCALCED:
        {
            result = a;

            break;
        }

        case DATA_ONE_UNCALCED:
        case DATA_ONE_ERROR:
        case DATA_ZERO_UNCALCED:
        case DATA_ZERO_ERROR:
        case DATA_POS_INFTY_ERROR:
        case DATA_NEG_INFTY_ERROR:
        case DATA_FIN_INDET_ERROR:
        case DATA_INF_INDET_ERROR:
        case DATA_POS_INFTY_UNCALCED:
        case DATA_NEG_INFTY_UNCALCED:
        case DATA_FIN_INDET_UNCALCED:
        case DATA_INF_INDET_UNCALCED:
        case DATA_UNCALCED_ERROR:
        case DATA_INTEGER_UNCALCED:
        case DATA_INTEGER_ERROR:
        case DATA_REAL_UNCALCED:
        case DATA_REAL_ERROR:
        {
            result = x;

            break;
        }

        case DATA_ZERO_ZERO:
        case DATA_ONE_ZERO:
        case DATA_INTEGER_ZERO:
        case DATA_REAL_ZERO:
        case DATA_POS_INFTY_ZERO:
        case DATA_NEG_INFTY_ZERO:
        case DATA_FIN_INDET_ZERO:
        case DATA_INF_INDET_ZERO:
        {
            result = Gen_gamma(a);

            break;
        }

        case DATA_ZERO_POS_INFTY:
        case DATA_ONE_POS_INFTY:
        case DATA_INTEGER_POS_INFTY:
        case DATA_REAL_POS_INFTY:
        case DATA_NEG_INFTY_POS_INFTY:
        case DATA_FIN_INDET_POS_INFTY:
        case DATA_NEG_INFTY_ONE:
        {
            result = ResultZero();

            break;
        }

        case DATA_POS_INFTY_POS_INFTY:
        case DATA_INF_INDET_POS_INFTY:
        {
            result = ResultInfIndet();

            break;
        }

        case DATA_ZERO_ONE:
        case DATA_ZERO_INTEGER:
        case DATA_ZERO_REAL:
        case DATA_ONE_ONE:
        case DATA_ONE_INTEGER:
        case DATA_ONE_REAL:
        case DATA_INTEGER_ONE:
        case DATA_INTEGER_INTEGER:
        case DATA_INTEGER_REAL:
        case DATA_REAL_ONE:
        case DATA_REAL_INTEGER:
        case DATA_REAL_REAL:
        {
            if ( getReal(x) < 0 )
            {
                result = ResultError(GVAR_GSL_EDOM);
            }

            else
            {
                switch ( ( ires = gsl_sf_gamma_inc_e(getReal(a),getReal(x),&gres) ) )
                {
                    case GSL_SUCCESS:
                    {
                        result = ResultReal(gres.val);

                        break;
                    }

                    default:
                    {
                        result = ResultError(convGSLErr(ires));

                        break;
                    }
                }
            }

            break;
        }

        case DATA_ZERO_NEG_INFTY:
        case DATA_ZERO_FIN_INDET:
        case DATA_ZERO_INF_INDET:
        case DATA_ONE_NEG_INFTY:
        case DATA_ONE_FIN_INDET:
        case DATA_ONE_INF_INDET:
        case DATA_INTEGER_NEG_INFTY:
        case DATA_INTEGER_FIN_INDET:
        case DATA_INTEGER_INF_INDET:
        case DATA_REAL_NEG_INFTY:
        case DATA_REAL_FIN_INDET:
        case DATA_REAL_INF_INDET:
        case DATA_POS_INFTY_NEG_INFTY:
        case DATA_POS_INFTY_FIN_INDET:
        case DATA_POS_INFTY_INF_INDET:
        case DATA_NEG_INFTY_NEG_INFTY:
        case DATA_NEG_INFTY_FIN_INDET:
        case DATA_NEG_INFTY_INF_INDET:
        case DATA_FIN_INDET_NEG_INFTY:
        case DATA_FIN_INDET_FIN_INDET:
        case DATA_FIN_INDET_INF_INDET:
        case DATA_INF_INDET_NEG_INFTY:
        case DATA_INF_INDET_FIN_INDET:
        case DATA_INF_INDET_INF_INDET:
        case DATA_UNCALCED_NEG_INFTY:
        case DATA_UNCALCED_FIN_INDET:
        case DATA_UNCALCED_INF_INDET:
        {
            result = ResultError(GVAR_GSL_EDOM);

            break;
        }

        case DATA_POS_INFTY_ONE:
        case DATA_POS_INFTY_INTEGER:
        case DATA_POS_INFTY_REAL:
        {
            if ( getReal(x) < 0 )
            {
                result = ResultError(GVAR_GSL_EDOM);
            }

            else
            {
                result = ResultPosInfnty();
            }

            break;
        }

        case DATA_NEG_INFTY_INTEGER:
        case DATA_NEG_INFTY_REAL:
        {
            if ( getReal(x) < 0 )
            {
                result = ResultError(GVAR_GSL_EDOM);
            }

            else
            {
                result = ResultZero();
            }

            break;
        }

        case DATA_FIN_INDET_ONE:
        case DATA_FIN_INDET_INTEGER:
        case DATA_FIN_INDET_REAL:
        case DATA_INF_INDET_ONE:
        case DATA_INF_INDET_INTEGER:
        case DATA_INF_INDET_REAL:
        {
            if ( getReal(x) < 0 )
            {
                result = ResultError(GVAR_GSL_EDOM);
            }

            else
            {
                result = ResultInfIndet();
            }

            break;
        }

        case DATA_UNCALCED_INTEGER:
        case DATA_UNCALCED_REAL:
        {
            if ( getReal(x) < 0 )
            {
                result = ResultError(GVAR_GSL_EDOM);
            }

            else
            {
                result = a;
            }

            break;
        }

        case DATA_ERROR_INTEGER:
        case DATA_ERROR_REAL:
        {
            if ( getReal(x) < 0 )
            {
                result = combineErrors(a,ResultError(GVAR_GSL_EDOM));
            }

            else
            {
                result = a;
            }

            break;
        }

        case DATA_ERROR_NEG_INFTY:
        case DATA_ERROR_FIN_INDET:
        case DATA_ERROR_INF_INDET:
        {
            result = combineErrors(a,ResultError(GVAR_GSL_EDOM));

            break;
        }

        case DATA_ERROR_ERROR:
        default:
        {
            result = combineErrors(a,x);

            break;
        }
    }

    return result;
}

GenVar Gen_erf(GenVar xx)
{
    GenVar result;
    gsl_sf_result gres;
    int ires;

    switch ( xx.DataType )
    {
        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        case DATA_IS_ONE:
        case DATA_IS_ZERO:
        {
            switch ( ( ires = gsl_sf_erf_e(getReal(xx),&gres) ) )
            {
                case GSL_SUCCESS:
                {
                    result = ResultReal(gres.val);

                    break;
                }

                case GSL_EUNDRFLW:
                case GSL_EOVRFLW:
                case GSL_ERANGE:
                {
                    if ( getReal(xx) > 0 )
                    {
                        result = ResultOne();
                    }

                    else
                    {
                        result = ResultInteger(-1);
                    }

                    break;
                }

                default:
                {
                    result = ResultError(convGSLErr(ires));

                    break;
                }
            }

            break;
        }

        case DATA_IS_POS_INFTY:
        {
            result = ResultOne();

            break;
        }

        case DATA_IS_NEG_INFTY:
        {
            result = ResultInteger(-1);

            break;
        }

        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        {
            result = ResultFinIndet();

            break;
        }

        case DATA_IS_UNCALCED:
        case DATA_IS_ERROR:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_erfc(GenVar xx)
{
    GenVar result;
    gsl_sf_result gres;
    int ires;

    switch ( xx.DataType )
    {
        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        case DATA_IS_ONE:
        case DATA_IS_ZERO:
        {
            switch ( ( ires = gsl_sf_erfc_e(getReal(xx),&gres) ) )
            {
                case GSL_SUCCESS:
                {
                    result = ResultReal(gres.val);

                    break;
                }

                case GSL_EUNDRFLW:
                case GSL_EOVRFLW:
                case GSL_ERANGE:
                {
                    if ( getReal(xx) > 0 )
                    {
                        result = ResultZero();
                    }

                    else
                    {
                        result = ResultInteger(2);
                    }

                    break;
                }

                default:
                {
                    result = ResultError(convGSLErr(ires));

                    break;
                }
            }

            break;
        }

        case DATA_IS_POS_INFTY:
        {
            result = ResultZero();

            break;
        }

        case DATA_IS_NEG_INFTY:
        {
            result = ResultInteger(2);

            break;
        }

        case DATA_IS_FIN_INDET:
        case DATA_IS_INF_INDET:
        {
            result = ResultFinIndet();

            break;
        }

        case DATA_IS_UNCALCED:
        case DATA_IS_ERROR:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_normDistr(GenVar xx)
{
    GenVar result;

    switch ( xx.DataType )
    {
        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        case DATA_IS_ONE:
        case DATA_IS_ZERO:
        {
            result = mulGenVar(ResultReal(0.398942280401),Gen_exp(negGenVar(mulGenVar(divGenVar(xx,ResultInteger(2)),xx))));

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
            result = ResultFinIndet();

            break;
        }

        case DATA_IS_UNCALCED:
        case DATA_IS_ERROR:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

GenVar Gen_polyDistr(GenVar nn, GenVar xx)
{
    GenVar result;
    GenVar tempa;
    GenVar tempb;
    GenVar tempc;
    GenVar cp;
    GenVar cp_prime;

    switch ( combineDataDescr(nn,xx) )
    {
        case DATA_ZERO_ZERO:
        case DATA_ZERO_ONE:
        case DATA_ZERO_INTEGER:
        case DATA_ZERO_REAL:
        case DATA_ZERO_POS_INFTY:
        case DATA_ZERO_NEG_INFTY:
        case DATA_ZERO_FIN_INDET:
        case DATA_ZERO_INF_INDET:
        case DATA_ZERO_UNCALCED:
        {
            result = ResultError(GVAR_GSL_EDOM);

            break;
        }

        case DATA_ZERO_ERROR:
        {
            result = combineErrors(ResultError(GVAR_GSL_EDOM),xx);

            break;
        }

        case DATA_ONE_ZERO:
        case DATA_ONE_ONE:
        case DATA_ONE_INTEGER:
        case DATA_ONE_REAL:
        {
            solve_polydist:

            tempa = Gen_gamma(divGenVar(ResultOne(),nn));
            tempb = Gen_gamma(divGenVar(ResultInteger(3),nn));

            tempc = Gen_sqrt(divGenVar(tempb,tempa));

            cp       = divGenVar(mulGenVar(nn,tempc),mulGenVar(ResultInteger(2),tempa));
            cp_prime = Gen_pow(tempc,nn);

            xx = Gen_abs(xx);

            result = mulGenVar(cp,Gen_exp(negGenVar(mulGenVar(cp_prime,Gen_pow(xx,nn)))));

            break;
        }

        case DATA_ONE_POS_INFTY:
        case DATA_ONE_NEG_INFTY:
        {
            result = ResultZero();

            break;
        }

        case DATA_ONE_FIN_INDET:
        case DATA_ONE_INF_INDET:
        {
            result = ResultFinIndet();

            break;
        }

        case DATA_ONE_UNCALCED:
        case DATA_ONE_ERROR:
        {
            result = xx;

            break;
        }

        case DATA_INTEGER_ZERO:
        case DATA_INTEGER_ONE:
        case DATA_INTEGER_INTEGER:
        case DATA_INTEGER_REAL:
        case DATA_REAL_ZERO:
        case DATA_REAL_ONE:
        case DATA_REAL_INTEGER:
        case DATA_REAL_REAL:
        {
            if ( getReal(nn) < 0 )
            {
                result = ResultError(GVAR_GSL_EDOM);
            }

            else
            {
                goto solve_polydist;
            }

            break;
        }

        case DATA_INTEGER_POS_INFTY:
        case DATA_INTEGER_NEG_INFTY:
        case DATA_REAL_POS_INFTY:
        case DATA_REAL_NEG_INFTY:
        {
            if ( getReal(nn) < 0 )
            {
                result = ResultError(GVAR_GSL_EDOM);
            }

            else
            {
                result = ResultZero();
            }

            break;
        }

        case DATA_INTEGER_FIN_INDET:
        case DATA_INTEGER_INF_INDET:
        case DATA_REAL_FIN_INDET:
        case DATA_REAL_INF_INDET:
        {
            if ( getReal(nn) < 0 )
            {
                result = ResultError(GVAR_GSL_EDOM);
            }

            else
            {
                result = ResultFinIndet();
            }

            break;
        }

        case DATA_INTEGER_UNCALCED:
        case DATA_REAL_UNCALCED:
        {
            if ( getReal(nn) < 0 )
            {
                result = ResultError(GVAR_GSL_EDOM);
            }

            else
            {
                result = xx;
            }

            break;
        }

        case DATA_INTEGER_ERROR:
        case DATA_REAL_ERROR:
        {
            if ( getReal(nn) < 0 )
            {
                result = combineErrors(ResultError(GVAR_GSL_EDOM),xx);
            }

            else
            {
                result = xx;
            }

            break;
        }

        case DATA_POS_INFTY_ZERO:
        case DATA_POS_INFTY_ONE:
        case DATA_POS_INFTY_INTEGER:
        case DATA_POS_INFTY_REAL:
        case DATA_POS_INFTY_POS_INFTY:
        case DATA_POS_INFTY_NEG_INFTY:
        {
            result = Gen_diracDelta(ResultZero(),xx);

            break;
        }

        case DATA_POS_INFTY_FIN_INDET:
        case DATA_POS_INFTY_INF_INDET:
        {
            result = ResultInfIndet();

            break;
        }

        case DATA_POS_INFTY_UNCALCED:
        case DATA_POS_INFTY_ERROR:
        {
            result = xx;

            break;
        }

        case DATA_NEG_INFTY_ZERO:
        case DATA_NEG_INFTY_ONE:
        case DATA_NEG_INFTY_INTEGER:
        case DATA_NEG_INFTY_REAL:
        case DATA_NEG_INFTY_POS_INFTY:
        case DATA_NEG_INFTY_NEG_INFTY:
        case DATA_NEG_INFTY_FIN_INDET:
        case DATA_NEG_INFTY_INF_INDET:
        case DATA_NEG_INFTY_UNCALCED:
        case DATA_FIN_INDET_ZERO:
        case DATA_FIN_INDET_ONE:
        case DATA_FIN_INDET_INTEGER:
        case DATA_FIN_INDET_REAL:
        case DATA_FIN_INDET_POS_INFTY:
        case DATA_FIN_INDET_NEG_INFTY:
        case DATA_FIN_INDET_FIN_INDET:
        case DATA_FIN_INDET_INF_INDET:
        case DATA_FIN_INDET_UNCALCED:
        case DATA_INF_INDET_ZERO:
        case DATA_INF_INDET_ONE:
        case DATA_INF_INDET_INTEGER:
        case DATA_INF_INDET_REAL:
        case DATA_INF_INDET_POS_INFTY:
        case DATA_INF_INDET_NEG_INFTY:
        case DATA_INF_INDET_FIN_INDET:
        case DATA_INF_INDET_INF_INDET:
        case DATA_INF_INDET_UNCALCED:
        {
            result = ResultError(GVAR_GSL_EDOM);

            break;
        }

        case DATA_NEG_INFTY_ERROR:
        case DATA_FIN_INDET_ERROR:
        case DATA_INF_INDET_ERROR:
        {
            result = combineErrors(ResultError(GVAR_GSL_EDOM),xx);

            break;
        }

        case DATA_UNCALCED_ZERO:
        case DATA_UNCALCED_ONE:
        case DATA_UNCALCED_INTEGER:
        case DATA_UNCALCED_REAL:
        case DATA_UNCALCED_POS_INFTY:
        case DATA_UNCALCED_NEG_INFTY:
        case DATA_UNCALCED_FIN_INDET:
        case DATA_UNCALCED_INF_INDET:
        case DATA_UNCALCED_UNCALCED:
        {
            result = nn;

            break;
        }

        case DATA_UNCALCED_ERROR:
        {
            result = xx;

            break;
        }

        case DATA_ERROR_ZERO:
        case DATA_ERROR_ONE:
        case DATA_ERROR_INTEGER:
        case DATA_ERROR_REAL:
        case DATA_ERROR_POS_INFTY:
        case DATA_ERROR_NEG_INFTY:
        case DATA_ERROR_FIN_INDET:
        case DATA_ERROR_INF_INDET:
        case DATA_ERROR_UNCALCED:
        {
            result = nn;

            break;
        }

        case DATA_ERROR_ERROR:
        default:
        {
            result = combineErrors(nn,xx);

            break;
        }
    }

    return result;
}

GenVar Gen_dawson(GenVar xx)
{
    GenVar result;
    gsl_sf_result gres;
    int ires;

    switch ( xx.DataType )
    {
        case DATA_IS_INTEGER:
        case DATA_IS_REAL:
        case DATA_IS_ONE:
        case DATA_IS_ZERO:
        {
            switch ( ( ires = gsl_sf_dawson_e(getReal(xx),&gres) ) )
            {
                case GSL_SUCCESS:
                {
                    result = ResultReal(gres.val);

                    break;
                }

                case GSL_EUNDRFLW:
                {
                    result = ResultZero();

                    break;
                }

                default:
                {
                    result = ResultError(convGSLErr(ires));

                    break;
                }
            }

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
            result = ResultFinIndet();

            break;
        }

        case DATA_IS_UNCALCED:
        case DATA_IS_ERROR:
        default:
        {
            result = xx;

            break;
        }
    }

    return result;
}

