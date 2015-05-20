
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
// Levi civita function.
//
// Written by: Alistair Shilton
//             Melbourne University
//

#include "levicivita.h"


//
// levi_civita: calculate levi-civita terms.
//              dim = # terms in expansion (>= 1)
//              terms = some permutation of {1,2,...,dim}
//                      (or a subset repeated to make dim terms)
//

int levi_civita(long dim, long *terms)
{
    long i;
    long term_one_marker = 0;
    int result = +1;

    if ( dim == 0 )
    {
        return 1;
    }

    // are we done?

    for ( i = 1 ; i <= dim ; i++ )
    {
        if ( terms[i-1] != i )
        {
            goto not_ordered;
        }
    }

    return 1;

    not_ordered:

    if ( dim == 1 )
    {
        return 0;
    }

    // find the first term
    // if we can't find it, return zero

    for ( i = 1 ; i <= dim ; i++ )
    {
        // have we found a repeated term?
        // if so, return 0 - no need to proceed further

        if ( terms[i-1] < 1 )
        {
            return 0;
        }

        // get out if we have found term 1

        if ( terms[i-1] == 1 )
        {
            term_one_marker = i;

            goto found_one;
        }
    }

    // term 1 does not exist, so result is 0

    return 0;

    found_one:

    // if not the first term, negate result.

    if ( term_one_marker != 1 )
    {
        result = -1;
    }

    // swap to first term

    i                        = terms[0];
    terms[0]                 = terms[term_one_marker-1];
    terms[term_one_marker-1] = i;

    // decrement the remainder of the sequence

    for ( i = 2 ; i <= dim ; i++ )
    {
        terms[i-1]--;
    }

    // recurse

    result *= levi_civita(dim-1,terms+1);

    // get rid of increment

    for ( i = 2 ; i <= dim ; i++ )
    {
        terms[i-1]++;
    }

    // swap to first term back

    i                        = terms[0];
    terms[0]                 = terms[term_one_marker-1];
    terms[term_one_marker-1] = i;

    return result;
}

