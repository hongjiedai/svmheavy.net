#include "stdafx.h"
/*
 *  SVMheavy - Another SVM Library
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


//
// Sparse vector functions.
//
// Written by: Alistair Shilton
//             Melbourne University
//




#include <iostream>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "sparsevector.h"
#include "svdefs.h"
#include "search.h"
#include "c_double.h"
#include "vector.h"


double dot_prod_sparse(const fVECTOR &x, const fVECTOR &y)
{
    long i,j,k,l,m,n;
    double result;

    result = 0;

    k = (long) x.get_offset_element(0);
    l = (long) y.get_offset_element(0);

    if ( ( k >= 1 ) && ( l >= 1 ) )
    {
        i = 1;
        j = 1;

        while ( ( i <= k ) && ( j <= l ) )
        {
            m = (long) x.get_offset_element((2*i)-1);
            n = (long) y.get_offset_element((2*j)-1);

            if ( n == m )
            {
                result += ( ((double) (x.get_offset_element(2*i))) * ((double) (y.get_offset_element(2*j))) );

                i++;
                j++;
            }

            else
            {
                if ( n > m )
                {
                    i++;
                }

                else
                {
                    j++;
                }
            }
        }
    }

    return result;
}

double dif_meas_sparse(const fVECTOR &x, const fVECTOR &y)
{
    // Note that (x-y)'(x-y) = x'x + y'y - 2x'y

    long i,j,k,l,m,n;
    int i_done;
    int j_done;
    double result;

    result = 0;

    k = (long) x.get_offset_element(0);
    l = (long) y.get_offset_element(0);

    if ( ( k >= 1 ) && ( l >= 1 ) )
    {
        i = 1;
        j = 1;

        i_done = 0;
        j_done = 0;

        // This loop is basically the dot product, but the x'y terms are now
        // scaled by -2, and x'x and y'y terms have been added.

        while ( ( i <= k ) && ( j <= l ) )
        {
            m = (long) x.get_offset_element((2*i)-1);
            n = (long) y.get_offset_element((2*j)-1);

            if ( !i_done )
            {
                result += ( ((double) (x.get_offset_element(2*i))) * ((double) (x.get_offset_element(2*i))) );

                i_done = 1;
            }

            if ( !j_done )
            {
                result += ( ((double) (y.get_offset_element(2*j))) * ((double) (y.get_offset_element(2*j))) );

                j_done = 1;
            }

            if ( n == m )
            {
                result -= ( 2 * ((double) (x.get_offset_element(2*i))) * ((double) (y.get_offset_element(2*j))) );

                i++;
                j++;

                i_done = 0;
                j_done = 0;
            }

            else
            {
                if ( n > m )
                {
                    i++;

                    i_done = 0;
                }

                else
                {
                    j++;

                    j_done = 0;
                }
            }
        }

        // Now do the remainder of the x'x and y'y terms

        if ( i_done )
        {
            i++;
        }

        if ( j_done )
        {
            j++;
        }

        while ( i <= k )
        {
            result += ( ((double) (x.get_offset_element(2*i))) * ((double) (x.get_offset_element(2*i))) );

            i++;
        }

        while ( j <= l )
        {
            result += ( ((double) (y.get_offset_element(2*j))) * ((double) (y.get_offset_element(2*j))) );

            j++;
        }
    }

    return result;
}

double eucldist_sparse(const fVECTOR &x, const fVECTOR &y)
{
    return sqrt(dif_meas_sparse(x,y));
}


fVECTOR &vector_scale_sparse(fVECTOR &x, double a)
{
    long dim;
    long i;

    dim = (long) x[0];

    if ( dim > 0 )
    {
        for ( i = 1 ; i <= dim ; i++ )
        {
            x[2*i] *= a;
        }
    }

    return x;
}

fVECTOR &vector_add_sparse(fVECTOR &x, const fVECTOR &y)
{
    long dimi,dimj;
    long i,j,n,m;

    dimi = (long) x.get_offset_element(0);
    dimj = (long) y.get_offset_element(0);

    if ( dimj > 0 )
    {
        i = 1;

        for ( j = 1 ; j <= dimj ; j++ )
        {
            if ( i <= dimi )
            {
                m = (long) x.get_offset_element((2*i)-1);
                n = (long) y.get_offset_element((2*j)-1);

                while ( n > m )
                {
                    i++;

                    if ( i > dimi )
                    {
                        goto addmethod;
                    }

                    m = (long) x.get_offset_element((2*i)-1);
                }

                if ( n == m )
                {
                    x[2*i] += y.get_offset_element(2*j);

                    i++;
                }

                else
                {
                    x.addend();
                    x.addend();

                    x.bswap(i+1,(2*dimi)+3);
                    x.bswap(i+2,(2*dimi)+2);

                    x[(2*i)+1] = (double) n;
                    x[(2*i)+2] = y.get_offset_element(2*j);

                    i++;
                    dimi++;
                }
            }

            else
            {
                addmethod:

                x.addend();
                x.addend();

                x[(2*i)+1] = y.get_offset_element((2*j)-1);
                x[(2*i)+2] = y.get_offset_element(2*j);

                i++;
                dimi++;
            }
        }
    }

    return x;
}

fVECTOR &vector_sub_sparse(fVECTOR &x, const fVECTOR &y)
{
    long dimi,dimj;
    long i,j,n,m;

    dimi = (long) x.get_offset_element(0);
    dimj = (long) y.get_offset_element(0);

    if ( dimj > 0 )
    {
        i = 1;

        for ( j = 1 ; j <= dimj ; j++ )
        {
            if ( i <= dimi )
            {
                m = (long) x.get_offset_element((2*i)-1);
                n = (long) y.get_offset_element((2*j)-1);

                while ( n > m )
                {
                    i++;

                    if ( i > dimi )
                    {
                        goto addmethod;
                    }

                    m = (long) x.get_offset_element((2*i)-1);
                }

                if ( n == m )
                {
                    x[2*i] -= y.get_offset_element(2*j);

                    i++;
                }

                else
                {
                    x.addend();
                    x.addend();

                    x.bswap(i+1,(2*dimi)+3);
                    x.bswap(i+2,(2*dimi)+2);

                    x[(2*i)+1] = (double) n;
                    x[(2*i)+2] = -(y.get_offset_element(2*j));

                    i++;
                    dimi++;
                }
            }

            else
            {
                addmethod:

                x.addend();
                x.addend();

                x[(2*i)+1] = y.get_offset_element((2*j)-1);
                x[(2*i)+2] = -(y.get_offset_element(2*j));

                i++;
                dimi++;
            }
        }
    }

    return x;
}

fVECTOR &convert_standard_format(char *input, fVECTOR &result)
{
    long d = 0;
    long len = 0;
    char *inrun;
    char *xpoint;
    long i,j;

    while ( input[len] != '\0' )
    {
        if ( input[len] == ':' )
        {
            d++;
        }

        len++;
    }

    result.trim_to_size((2*d)+1);

    result[0] = (L_DOUBLE) d;

    inrun = input;

    if ( d > 0 )
    {
		// TODO delete marker and value?
		// long marker[d];
		// L_DOUBLE values[d];
        long *marker = new long[d];
        L_DOUBLE *values = new L_DOUBLE[d];
        long mtemp;
        L_DOUBLE vtemp;

        /*
           Get the data
        */

        for ( i = 1 ; i <= d ; i++ )
        {
            /*
               First off, seek the :
            */

            xpoint = inrun;

            while ( *inrun != ':' )
            {
                inrun++;
            }

            /*
               Get the element number
            */

            *inrun = '\0';
            sscanf(xpoint,"%ld",&(marker[i-1]));
            *inrun = ':';

            inrun++;

            /*
               Get the value
            */

            double temp;

            sscanf(inrun,"%lf",&temp);

            values[i-1] = temp;

            /*
               Find the whitespace
            */

            if ( i < d )
            {
                while ( !isspace(*inrun) )
                {
                    inrun++;
                }
            }

            /*
               And cross it
            */

            while ( isspace(*inrun) )
            {
                inrun++;
            }
        }

        /*
           Sort the data
        */

        if ( d > 1 )
        {
            for ( i = 1 ; i < d ; i++ )
            {
                for ( j = i+1 ; j <= d ; j++ )
                {
                    if ( marker[j-1] < marker[i-1] )
                    {
                        mtemp = marker[i-1];
                        vtemp = values[i-1];

                        marker[i-1] = marker[j-1];
                        values[i-1] = values[j-1];

                        marker[j-1] = mtemp;
                        values[j-1] = vtemp;
                    }
                }
            }
        }

        /*
           Store the data
        */

        for ( i = 1 ; i <= d ; i++ )
        {
            result[(2*i)-1] = (L_DOUBLE) marker[i-1];
            result[(2*i)]   = values[i-1];
        }
    }

    return result;
}

fVECTOR &convert_full_format(char *buffer, fVECTOR &result)
{
    double temp;
    long dim;

    while ( isspace(*buffer) )
    {
        buffer++;
    }

    char *tbuff;

    tbuff = buffer;
    dim = 0;

    while ( *tbuff != '\0' )
    {
        while ( !isspace(*tbuff) )
        {
            tbuff++;
        }

        while ( isspace(*tbuff) )
        {
            tbuff++;
        }

        dim++;
    }

    result.trim_to_size(dim);

    tbuff = buffer;
    dim = 0;

    while ( *tbuff != '\0' )
    {
        sscanf(tbuff,"%lf",&temp);
        result[dim] = temp;

        while ( !isspace(*tbuff) )
        {
            tbuff++;
        }

        while ( isspace(*tbuff) )
        {
            tbuff++;
        }

        dim++;
    }

    return result;
}

f_fVECTOR &desparsifier(f_fVECTOR &input, f_fVECTOR &result)
{
    long dim;
    long maxdim;
    long numvects;
    long i,j;

    numvects = input.get_effective_size();

    /*
       Calculate the max dimension of data
    */

    maxdim = 0;

    if ( numvects > 0 )
    {
        for ( i = 1 ; i <= numvects ; i++ )
        {
            dim = (long) (input[i-1])[0];

            if ( dim > maxdim )
            {
                maxdim = dim;
            }
        }
    }

    /*
       Construct the desparsed version
    */

    result.trim_to_size(numvects);

    if ( numvects > 0 )
    {
        for ( i = 1 ; i <= numvects ; i++ )
        {
            (result[i-1]).trim_to_size(maxdim);
            (result[i-1]) = 0.0;

            dim = (long) (input[i-1])[0];

            if ( dim > 0 )
            {
                for ( j = 1 ; j <= dim ; j++ )
                {
                    (result[i-1])[((long) (input[i-1])[(2*j)-1])-1] = (input[i-1])[2*j];
                }
            }
        }
    }

    return result;
}



double dot_prod_sparse_opt(L_DOUBLE *x, L_DOUBLE *y)
{
    long i,j,k,l,m,n;
    double result;

    result = 0;

    k = (long) x[0];
    l = (long) y[0];

    if ( ( k >= 1 ) && ( l >= 1 ) )
    {
        i = 1;
        j = 1;

        while ( ( i <= k ) && ( j <= l ) )
        {
            m = (long) x[(2*i)-1];
            n = (long) y[(2*j)-1];

            if ( n == m )
            {
                result += ( ((double) (x[2*i])) * ((double) (y[2*j])) );

                i++;
                j++;
            }

            else
            {
                if ( n > m )
                {
                    i++;
                }

                else
                {
                    j++;
                }
            }
        }
    }

    return result;
}

double dif_meas_sparse_opt(L_DOUBLE *x, L_DOUBLE *y)
{
    // Note that (x-y)'(x-y) = x'x + y'y - 2x'y

    long i,j,k,l,m,n;
    int i_done;
    int j_done;
    double result;

    result = 0;

    k = (long) x[0];
    l = (long) y[0];

    if ( ( k >= 1 ) && ( l >= 1 ) )
    {
        i = 1;
        j = 1;

        i_done = 0;
        j_done = 0;

        // This loop is basically the dot product, but the x'y terms are now
        // scaled by -2, and x'x and y'y terms have been added.

        while ( ( i <= k ) && ( j <= l ) )
        {
            m = (long) x[(2*i)-1];
            n = (long) y[(2*j)-1];

            if ( !i_done )
            {
                result += ( ((double) x[2*i]) * ((double) x[2*i]) );

                i_done = 1;
            }

            if ( !j_done )
            {
                result += ( ((double) y[2*j]) * ((double) y[2*j]) );

                j_done = 1;
            }

            if ( n == m )
            {
                result -= ( 2 * ((double) x[2*i]) * ((double) y[2*j]) );

                i++;
                j++;

                i_done = 0;
                j_done = 0;
            }

            else
            {
                if ( n > m )
                {
                    i++;

                    i_done = 0;
                }

                else
                {
                    j++;

                    j_done = 0;
                }
            }
        }

        // Now do the remainder of the x'x and y'y terms

        if ( i_done )
        {
            i++;
        }

        if ( j_done )
        {
            j++;
        }

        while ( i <= k )
        {
            result += ( ((double) x[2*i]) * ((double) x[2*i]) );

            i++;
        }

        while ( j <= l )
        {
            result += ( ((double) y[2*j]) * ((double) y[2*j]) );

            j++;
        }
    }

    return result;
}

double eucldist_sparse_opt(L_DOUBLE *x, L_DOUBLE *y)
{
    return sqrt(dif_meas_sparse_opt(x,y));
}



