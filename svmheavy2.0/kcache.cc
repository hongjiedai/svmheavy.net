#include "stdafx.h"

/*
 *  SVMheavy - Another SVM Library
 *  Copyright (C) 2005  Alistair Shilton
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
// Kernel cache class
//
// Written by: Alistair Shilton
//             Melbourne University
//

#include "kcache.h"
#include "svdata.h"
#include "svflags.h"
#include "svdata.h"
#include "kernel.h"
#include "common/vector.h"
#include "common/matrix.h"
#include "common/factor.h"
#include <limits.h>

#define PRINT_CACHE_STATE                                               \
{                                                                       \
    long __i__;                                                         \
    Klink *__x__;                                                       \
                                                                        \
    __x__ = first_element;                                              \
                                                                        \
    for ( __i__ = 1 ; __i__ <= maxrows ; __i__++ )                      \
    {                                                                   \
        if ( __x__ == last_element )                                    \
        {                                                               \
            std::cerr << "+";                                           \
        }                                                               \
                                                                        \
        std::cerr << __x__->row_ident << ", ";                          \
                                                                        \
        __x__ = __x__->next;                                            \
    }                                                                   \
                                                                        \
    std::cerr << "\n";                                                  \
                                                                        \
    if ( __x__ != NULL )                                                \
    {                                                                   \
        std::cerr << "NULL ending missing\n";                           \
                                                                        \
        exit(1);                                                        \
    }                                                                   \
}

#define MIN_MAXROWS 5
#define MAXROWS_STRESS 5

Klink::Klink()
{
    kernel_row      = NULL;
    fast_kernel_row = NULL;
    row_ident       = 0;

    next = NULL;
    prev = NULL;

    return;
}

Klink::Klink(const fVECTOR &_kernel_row, long _row_ident)
{
    kernel_row      = _kernel_row;
    fast_kernel_row = NULL;
    row_ident       = _row_ident;

    next = NULL;
    prev = NULL;

    return;
}

Kcache::Kcache()
{
    #ifdef DEBCACHE
    std::cerr << "Call to constructor\n";
    #endif

    trainsize  = 0;
    maxrows    = 0;
    rowdim     = 0;
    memsize    = 0;
    min_rowdim = 0;

    first_element = NULL;
    last_element  = NULL;

    return;
}

void Kcache::reset(void)
{
    #ifdef DEBCACHE
    std::cerr << "Call to reset\n";
    #endif

    Klink *now;
    Klink *then;

    now = first_element;

    while ( now != NULL )
    {
        then = now->next;
        REPDEL(now);
        now = then;
    }

    trainsize  = 0;
    maxrows    = 0;
    rowdim     = 0;
    memsize    = 0;
    min_rowdim = 0;

    first_element = NULL;
    last_element  = NULL;

    return;
}

void Kcache::setmemsize(SVdata *caller, long _memsize, long _min_rowdim)
{
    #ifdef DEBCACHE
    std::cerr << "Call to setmemsize(" << _memsize << "," << _min_rowdim << ")\n";
    #endif

    THROW_ASSERT( _memsize > 0 );
    THROW_ASSERT( _min_rowdim > 0 );

    long i;
    Klink *pos_ptr;

    if ( memsize == 0 )
    {
        memsize    = _memsize;
        min_rowdim = _min_rowdim;

        /*
           Set all the relevant values.
        */

        trainsize = caller->N;

        rowdim = ( trainsize > min_rowdim ) ? trainsize : min_rowdim;

        if ( ( maxrows = (memsize*1024*1024)/(rowdim*sizeof(L_DOUBLE)) ) < MIN_MAXROWS )
        {
            maxrows = MIN_MAXROWS;
        }

        #ifdef DEBSTRESS
        maxrows = MAXROWS_STRESS;
        #endif

        /*
           Allocate the vectors.
        */

        fVECTOR temp('c',rowdim);

        REPNEW(first_element,Klink(temp,0));

        pos_ptr = first_element;

        if ( maxrows >= 2 )
        {
            for ( i = 2 ; i <= maxrows ; i++ )
            {
                REPNEW(pos_ptr->next,Klink(temp,0));

                (pos_ptr->next)->prev = pos_ptr;

                pos_ptr = pos_ptr->next;
            }
        }

        last_element = pos_ptr;

        #ifdef DEBCACHE
        std::cerr << "end of first go setmemsize:\n";
        PRINT_CACHE_STATE;
        std::cerr << "Finished\n\n";
        #endif
    }

    else
    {
        THROW_ASSERT( caller->N >= trainsize-1 );
        THROW_ASSERT( caller->N <= trainsize+1 );

        long new_trainsize;
        long new_rowdim;
        long new_maxrows;

        /*
           Assumption - (implicit) call afer remove operation.  The
           data has all been changed correctly.  ==OR== only the
           memsize stuff has changed.  In either case, easy.

           Failure of assumption if a single point added to end dealt with
           shortly.
        */

        memsize    = _memsize;
        min_rowdim = _min_rowdim;

        new_trainsize = caller->N;

        new_rowdim = ( new_trainsize > min_rowdim ) ? new_trainsize : min_rowdim;

        if ( ( new_maxrows = (memsize*1024*1024)/(new_rowdim*sizeof(L_DOUBLE)) ) < MIN_MAXROWS )
        {
            new_maxrows = MIN_MAXROWS;
        }

        #ifdef DEBSTRESS
        new_maxrows = MAXROWS_STRESS;
        #endif

        if ( new_rowdim != rowdim )
        {
            pos_ptr = first_element;

            while ( pos_ptr != NULL )
            {
                (pos_ptr->kernel_row).trim_to_size(new_rowdim);

                pos_ptr = pos_ptr->next;
            }
        }

        rowdim = new_rowdim;

        if ( new_maxrows > maxrows )
        {
            fVECTOR temp('c',rowdim);

            pos_ptr = last_element;

            for ( i = maxrows+1 ; i <= new_maxrows ; i++ )
            {
                REPNEW(pos_ptr->next,Klink(temp,0));

                (pos_ptr->next)->prev = pos_ptr;

                pos_ptr = pos_ptr->next;
            }

            last_element = pos_ptr;
        }

        else if ( new_maxrows < maxrows )
        {
            Klink *temp;

            pos_ptr = first_element;
            temp    = first_element;

            for ( i = 1 ; i <= new_maxrows ; i++ )
            {
                temp = pos_ptr;
                pos_ptr = pos_ptr->next;
            }

            last_element = temp;
            last_element->next = NULL;

            for ( i = new_maxrows+1 ; i <= maxrows ; i++ )
            {
                temp = pos_ptr;
                pos_ptr = pos_ptr->next;

                REPDEL(temp);
            }
        }

        maxrows = new_maxrows;

        if ( new_trainsize > trainsize )
        {
            /*
               Assumption - call after a single point was added to the
               training set.  Need to update the relevant element of all
               used rows
            */

            if ( first_element->row_ident != 0 )
            {
                pos_ptr = first_element;

                while ( pos_ptr != NULL )
                {
                    (pos_ptr->kernel_row)[new_trainsize-1] = (caller->K).kernel((caller->x)[(pos_ptr->row_ident)-1],(caller->x)[new_trainsize-1],NULL,NULL,0.0,0.0,(pos_ptr->row_ident),new_trainsize,(caller->z)[(pos_ptr->row_ident)-1],(caller->z)[new_trainsize-1]);

                    pos_ptr = pos_ptr->next;

                    if ( pos_ptr->row_ident == 0 )
                    {
                        break;
                    }
                }
            }
        }

        trainsize = new_trainsize;

        #ifdef DEBCACHE
        std::cerr << "finished later call to setmemsize:\n";
        PRINT_CACHE_STATE;
        std::cerr << "Finished\n\n";
        #endif
    }

    return;
}

L_DOUBLE Kcache::getval(SVdata *caller, long numi, long numj)
{
    #ifdef DEBCACHE
    std::cerr << "Call to getval(" << numi << "," << numj << ")\n";
    PRINT_CACHE_STATE;
    #endif

    THROW_ASSERT( memsize > 0 );
    THROW_ASSERT( numi >= 1 );
    THROW_ASSERT( numi <= trainsize );
    THROW_ASSERT( numj >= 1 );
    THROW_ASSERT( numj <= trainsize );

    Klink *pos_ptr;
    long j;

    if ( first_element->row_ident != 0 )
    {
        pos_ptr = first_element;

        while ( pos_ptr != NULL )
        {
            if ( pos_ptr->row_ident == numi )
            {
                break;
            }

            else if ( pos_ptr->row_ident == numj )
            {
                j    = numi;
                numi = numj;
                numj = j;

                break;
            }

            else if ( pos_ptr->row_ident == 0 )
            {
                break;
            }

            pos_ptr = pos_ptr->next;
        }
    }

    #ifdef DEBCACHE
    L_DOUBLE result;
    result = (*getrow(caller,numi))[numj-1];
    std::cerr << "result = " << result << "\n";
    return result;
    #endif

    return (*getrow(caller,numi))[numj-1];
}

L_DOUBLE Kcache::getval_statcache(SVdata *caller, long numi, long numj)
{
    #ifdef DEBCACHE
    std::cerr << "Call to getval_statcache(" << numi << "," << numj << ")\n";
    PRINT_CACHE_STATE;
    #endif

    THROW_ASSERT( memsize > 0 );
    THROW_ASSERT( numi >= 1 );
    THROW_ASSERT( numi <= trainsize );
    THROW_ASSERT( numj >= 1 );
    THROW_ASSERT( numj <= trainsize );

    Klink *pos_ptr;
    int isincache = 0;
    L_DOUBLE result = 0.0;

    if ( first_element->row_ident != 0 )
    {
        pos_ptr = first_element;

        while ( pos_ptr != NULL )
        {
            if ( pos_ptr->row_ident == numi )
            {
                isincache = 1;
                result    = (pos_ptr->kernel_row)[numj-1];

                break;
            }

            else if ( pos_ptr->row_ident == numj )
            {
                isincache = 1;
                result    = (pos_ptr->kernel_row)[numi-1];

                break;
            }

            else if ( pos_ptr->row_ident == 0 )
            {
                break;
            }

            pos_ptr = pos_ptr->next;
        }
    }

    if ( !isincache )
    {
        result = (caller->K).kernel((caller->x)[numi-1],(caller->x)[numj-1],NULL,NULL,0.0,0.0,numi,numj,(caller->z)[numi-1],(caller->z)[numj-1]);
    }

    #ifdef DEBCACHE
    std::cerr << "result = " << result << "\n";
    #endif

    return result;
}

fVECTOR *Kcache::getrow(SVdata *caller, long numi)
{
    #ifdef DEBCACHE
    std::cerr << "Call to getrow(" << numi << ")\n";
    PRINT_CACHE_STATE;
    #endif

    THROW_ASSERT( memsize > 0 );
    THROW_ASSERT( numi >= 1 );
    THROW_ASSERT( numi <= trainsize );

    long j;
    Klink *pos_ptr = NULL;

    if ( first_element->row_ident != 0 )
    {
        pos_ptr = first_element;

        while ( pos_ptr != NULL )
        {
            if ( pos_ptr->row_ident == numi )
            {
                /*
                   Put the row at the top of the list
                */

                if ( pos_ptr == last_element )
                {
                    #ifdef DEBCACHE
                    std::cerr << "found last one\n";
                    #endif

                    last_element = last_element->prev;
                    last_element->next = NULL;

                    pos_ptr->next = first_element;
                    pos_ptr->prev = NULL;

                    first_element->prev = pos_ptr;
                    first_element = pos_ptr;
                }

                else if ( pos_ptr != first_element )
                {
                    #ifdef DEBCACHE
                    std::cerr << "found middle one\n";
                    #endif

                    (pos_ptr->prev)->next = pos_ptr->next;
                    (pos_ptr->next)->prev = pos_ptr->prev;

                    pos_ptr->next = first_element;
                    pos_ptr->prev = NULL;

                    first_element->prev = pos_ptr;
                    first_element = pos_ptr;
                }

                break;
            }

            else if ( pos_ptr->row_ident == 0 )
            {
                pos_ptr = NULL;

                break;
            }

            pos_ptr = pos_ptr->next;
        }
    }

    /*
       If not found, construct the relevant row.
    */

    if ( pos_ptr == NULL )
    {
        #ifdef DEBCACHE
        std::cerr << "construction happens\n";
        #endif

        pos_ptr = last_element;

        last_element = last_element->prev;
        last_element->next = NULL;

        pos_ptr->next = first_element;
        pos_ptr->prev = NULL;

        first_element->prev = pos_ptr;
        first_element = pos_ptr;

        first_element->row_ident = numi;

        for ( j = 1 ; j <= trainsize ; j++ )
        {
            (first_element->kernel_row)[j-1] = (caller->K).kernel((caller->x)[numi-1],(caller->x)[j-1],NULL,NULL,0.0,0.0,numi,j,(caller->z)[numi-1],(caller->z)[j-1]);
        }
    }

    return &(first_element->kernel_row);
}

void Kcache::remove(SVdata *caller, long num)
{
    #ifdef DEBCACHE
    std::cerr << "Call to remove(" << num << ")\n";
    #endif

    THROW_ASSERT( memsize > 0 );
    THROW_ASSERT( num >= 1 );
    THROW_ASSERT( num <= trainsize );

    Klink *pos_ptr = NULL;
    Klink *pos_hold = NULL;

    /*
       See if it is in the cache, and if found remove it.  Also swap removed
       bit to end otherwise.
    */

    if ( first_element->row_ident != 0 )
    {
        pos_ptr = first_element;

        while ( pos_ptr != NULL )
        {
            if ( pos_ptr->row_ident == num )
            {
                pos_ptr->row_ident = 0;

                if ( pos_ptr == first_element )
                {
                    first_element = first_element->next;
                    first_element->prev = NULL;

                    pos_ptr->next = NULL;
                    pos_ptr->prev = last_element;

                    last_element->next = pos_ptr;
                    last_element = pos_ptr;

                    pos_ptr = first_element;
                }

                else if ( pos_ptr != last_element )
                {
                    pos_hold = pos_ptr->next;

                    (pos_ptr->prev)->next = pos_ptr->next;
                    (pos_ptr->next)->prev = pos_ptr->prev;

                    pos_ptr->next = NULL;
                    pos_ptr->prev = last_element;

                    last_element->next = pos_ptr;
                    last_element = pos_ptr;

                    pos_ptr = pos_hold;
                }
            }

            /* deliberately not an ELSE here */

            if ( pos_ptr->row_ident == 0 )
            {
                break;
            }

            (pos_ptr->kernel_row).fswap(num,trainsize);

            pos_ptr = pos_ptr->next;
        }
    }

    /*
       Finally, call setmemsize
    */

    setmemsize(caller,memsize,min_rowdim);

    return;
}

void Kcache::flush(void)
{
    #ifdef DEBCACHE
    std::cerr << "Call to flush()\n";
    #endif

    THROW_ASSERT( memsize > 0 );

    Klink *pos_ptr;

    if ( first_element->row_ident != 0 )
    {
        pos_ptr = first_element;

        while ( pos_ptr != NULL )
        {
            if ( pos_ptr->row_ident == 0 )
            {
                break;
            }

            pos_ptr->row_ident = 0;

            pos_ptr = pos_ptr->next;
        }
    }

    return;
}






























void Kcache::enter_opt(SVdata *caller)
{
    Klink *pos_ptr;
    long i;

    REPNEWB(fast_lookup,Klink *,trainsize);

    for ( i = 0 ; i < trainsize ; i++ )
    {
        fast_lookup[i] = last_element;
    }

    pos_ptr = first_element;

    while ( pos_ptr != NULL )
    {
        pos_ptr->fast_kernel_row = &((pos_ptr->kernel_row)[0]);

        if ( pos_ptr->row_ident != 0 )
        {
            fast_lookup[(pos_ptr->row_ident)-1] = pos_ptr;
        }

        pos_ptr = pos_ptr->next;
    }

    fast_kernel_fn = (caller->K).get_fast_kern(&fast_uc,&fast_uv,(double **) (&fast_covw));

    return;
}

void Kcache::exit_opt(void)
{
    REPDELB(fast_lookup);

    return;
}

L_DOUBLE Kcache::opt_getval_statcache(SVdata *caller, long numi, long numj)
{
    #ifdef DEBCACHE
    std::cerr << "Call to getval_statcache(" << numi << "," << numj << ")\n";
    PRINT_CACHE_STATE;
    #endif

    THROW_ASSERT( memsize > 0 );
    THROW_ASSERT( numi >= 1 );
    THROW_ASSERT( numi <= trainsize );
    THROW_ASSERT( numj >= 1 );
    THROW_ASSERT( numj <= trainsize );

    L_DOUBLE result = 0.0;

    #ifdef USEMINORCACHE
    if ( (fast_lookup[numi-1])->row_ident == numi )
    {
        result = ((fast_lookup[numi-1])->fast_kernel_row)[numj-1];
    }

    else if ( (fast_lookup[numj-1])->row_ident == numj )
    {
        result = ((fast_lookup[numi-1])->fast_kernel_row)[numi-1];
    }

    else
    {
        if ( (caller->fast_x) == NULL )
        {
            result = (caller->K).kernel((caller->x)[numi-1],(caller->x)[numj-1],NULL,NULL,0.0,0.0,numi,numj,(caller->z)[numi-1],(caller->z)[numj-1]);
        }

        else
        {
            if ( fast_kernel_fn == NULL )
            {
                result = (caller->K).kernel(*((caller->fast_x)[numi-1]),*((caller->fast_x)[numj-1]),(caller->fast_opt_x)[numi-1],(caller->fast_opt_x)[numj-1],0.0,0.0,numi,numj,(caller->fast_z)[numi-1],(caller->fast_z)[numj-1]);
            }

            else
            {
                result = fast_kernel_fn((caller->fast_opt_x)[numi-1],(caller->fast_opt_x)[numj-1],fast_uc,fast_uv,fast_covw);
            }
        }
    }
    #endif

    #ifndef USEMINORCACHE
    if ( (caller->fast_x) == NULL )
    {
        result = (caller->K).kernel((caller->x)[numi-1],(caller->x)[numj-1],NULL,NULL,0.0,0.0,numi,numj,(caller->z)[numi-1],(caller->z)[numj-1]);
    }

    else
    {
        if ( fast_kernel_fn == NULL )
        {
            result = (caller->K).kernel(*((caller->fast_x)[numi-1]),*((caller->fast_x)[numj-1]),(caller->fast_opt_x)[numi-1],(caller->fast_opt_x)[numj-1],0.0,0.0,numi,numj,(caller->fast_z)[numi-1],(caller->fast_z)[numj-1]);
        }

        else
        {
            result = fast_kernel_fn((caller->fast_opt_x)[numi-1],
            (caller->fast_opt_x)[numj-1],fast_uc,fast_uv,(double *) fast_covw);
        }
    }
    #endif

    #ifdef DEBCACHE
    std::cerr << "result = " << result << "\n";
    #endif

    return result;
}

L_DOUBLE *Kcache::opt_getrow(SVdata *caller, long numi)
{
    #ifdef DEBCACHE
    std::cerr << "Call to getrow(" << numi << ")\n";
    PRINT_CACHE_STATE;
    #endif

    THROW_ASSERT( memsize > 0 );
    THROW_ASSERT( numi >= 1 );
    THROW_ASSERT( numi <= trainsize );

    long j;

    if ( (fast_lookup[numi-1])->row_ident == numi )
    {
        /*
           Put the row at the top of the list
        */

        if ( (fast_lookup[numi-1]) == last_element )
        {
            #ifdef DEBCACHE
            std::cerr << "found last one\n";
            #endif

            last_element = last_element->prev;
            last_element->next = NULL;

            (fast_lookup[numi-1])->next = first_element;
            (fast_lookup[numi-1])->prev = NULL;

            first_element->prev = (fast_lookup[numi-1]);
            first_element = (fast_lookup[numi-1]);
        }

        else if ( (fast_lookup[numi-1]) != first_element )
        {
            #ifdef DEBCACHE
            std::cerr << "found middle one\n";
            #endif

            ((fast_lookup[numi-1])->prev)->next = (fast_lookup[numi-1])->next;
            ((fast_lookup[numi-1])->next)->prev = (fast_lookup[numi-1])->prev;

            (fast_lookup[numi-1])->next = first_element;
            (fast_lookup[numi-1])->prev = NULL;

            first_element->prev = (fast_lookup[numi-1]);
            first_element = (fast_lookup[numi-1]);
        }
    }

    else
    {
        #ifdef DEBCACHE
        std::cerr << "construction happens\n";
        #endif

        (fast_lookup[numi-1]) = last_element;

        last_element = last_element->prev;
        last_element->next = NULL;

        (fast_lookup[numi-1])->next = first_element;
        (fast_lookup[numi-1])->prev = NULL;

        first_element->prev = (fast_lookup[numi-1]);
        first_element= (fast_lookup[numi-1]);

        first_element->row_ident = numi;

        fVECTOR **fast_x;
        L_DOUBLE **fast_opt_x;
        L_DOUBLE *fast_z;

        fast_x     = caller->fast_x;
        fast_opt_x = caller->fast_opt_x;
        fast_z     = caller->fast_z;

        if ( fast_x == NULL )
        {
            for ( j = 1 ; j <= trainsize ; j++ )
            {
                (first_element->fast_kernel_row)[j-1] = (caller->K).kernel((caller->x)[numi-1],(caller->x)[j-1],NULL,NULL,0.0,0.0,numi,j,fast_z[numi-1],fast_z[j-1]);
            }
        }

        else
        {
            if ( fast_kernel_fn == NULL )
            {
                for ( j = 1 ; j <= trainsize ; j++ )
                {
                    (first_element->fast_kernel_row)[j-1] = (caller->K).kernel(*(fast_x[numi-1]),*(fast_x[j-1]),fast_opt_x[numi-1],fast_opt_x[j-1],0.0,0.0,numi,j,fast_z[numi-1],fast_z[j-1]);
                }
            }

            else
            {
                for ( j = 1 ; j <= trainsize ; j++ )
                {
                    (first_element->fast_kernel_row)[j-1] = fast_kernel_fn(fast_opt_x[numi-1],fast_opt_x[j-1],fast_uc,fast_uv,(double *) fast_covw);
                }
            }
        }
    }

    return first_element->fast_kernel_row;
}

L_DOUBLE Kcache::opt_getval_statcache_fast(SVdata *caller, long numi, long numj)
{
    #ifdef DEBCACHE
    std::cerr << "Call to getval_statcache(" << numi << "," << numj << ")\n";
    PRINT_CACHE_STATE;
    #endif

    THROW_ASSERT( memsize > 0 );
    THROW_ASSERT( numi >= 1 );
    THROW_ASSERT( numi <= trainsize );
    THROW_ASSERT( numj >= 1 );
    THROW_ASSERT( numj <= trainsize );

    L_DOUBLE result = 0.0;

    #ifdef USEMINORCACHE
    if ( (fast_lookup[numi-1])->row_ident == numi )
    {
        result = ((fast_lookup[numi-1])->fast_kernel_row)[numj-1];
    }

    else if ( (fast_lookup[numj-1])->row_ident == numj )
    {
        result = ((fast_lookup[numi-1])->fast_kernel_row)[numi-1];
    }

    else
    {
        result = fast_kernel_fn((caller->fast_opt_x)[numi-1],(caller->fast_opt_x)[numj-1],fast_uc,fast_uv,fast_covw);
    }
    #endif

    #ifndef USEMINORCACHE
    result = fast_kernel_fn((caller->fast_opt_x)[numi-1],(caller->fast_opt_x)[numj-1],fast_uc,fast_uv,(double *) fast_covw);
    #endif

    #ifdef DEBCACHE
    std::cerr << "result = " << result << "\n";
    #endif

    return result;
}

L_DOUBLE *Kcache::opt_getrow_fast(SVdata *caller, long numi)
{
    #ifdef DEBCACHE
    std::cerr << "Call to getrow(" << numi << ")\n";
    PRINT_CACHE_STATE;
    #endif

    THROW_ASSERT( memsize > 0 );
    THROW_ASSERT( numi >= 1 );
    THROW_ASSERT( numi <= trainsize );

    long j;

    if ( (fast_lookup[numi-1])->row_ident == numi )
    {
        /*
           Put the row at the top of the list
        */

        if ( (fast_lookup[numi-1]) == last_element )
        {
            #ifdef DEBCACHE
            std::cerr << "found last one\n";
            #endif

            last_element = last_element->prev;
            last_element->next = NULL;

            (fast_lookup[numi-1])->next = first_element;
            (fast_lookup[numi-1])->prev = NULL;

            first_element->prev = (fast_lookup[numi-1]);
            first_element = (fast_lookup[numi-1]);
        }

        else if ( (fast_lookup[numi-1]) != first_element )
        {
            #ifdef DEBCACHE
            std::cerr << "found middle one\n";
            #endif

            ((fast_lookup[numi-1])->prev)->next = (fast_lookup[numi-1])->next;
            ((fast_lookup[numi-1])->next)->prev = (fast_lookup[numi-1])->prev;

            (fast_lookup[numi-1])->next = first_element;
            (fast_lookup[numi-1])->prev = NULL;

            first_element->prev = (fast_lookup[numi-1]);
            first_element = (fast_lookup[numi-1]);
        }
    }

    else
    {
        #ifdef DEBCACHE
        std::cerr << "construction happens\n";
        #endif

        (fast_lookup[numi-1]) = last_element;

        last_element = last_element->prev;
        last_element->next = NULL;

        (fast_lookup[numi-1])->next = first_element;
        (fast_lookup[numi-1])->prev = NULL;

        first_element->prev = (fast_lookup[numi-1]);
        first_element= (fast_lookup[numi-1]);

        first_element->row_ident = numi;

        fVECTOR **fast_x;
        L_DOUBLE **fast_opt_x;
        L_DOUBLE *fast_z;

        fast_x     = caller->fast_x;
        fast_opt_x = caller->fast_opt_x;
        fast_z     = caller->fast_z;

        for ( j = 1 ; j <= trainsize ; j++ )
        {
            (first_element->fast_kernel_row)[j-1] = fast_kernel_fn(fast_opt_x[numi-1],fast_opt_x[j-1],fast_uc,fast_uv,(double *) fast_covw);
        }
    }

    return first_element->fast_kernel_row;
}

