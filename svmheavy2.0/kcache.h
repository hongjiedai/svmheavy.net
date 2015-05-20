
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

#ifndef _kcache_h
#define _kcache_h

#include "svflags.h"
#include "kernel.h"
#include "common/vector.h"
#include "common/matrix.h"
#include "common/factor.h"

class Kcache;
class Klink;
class SVdata;

class Klink
{
    /*
       kernel_row = the kernel row.
       row_ident  = row number (or zero if not used).
    */

    public:

    Klink();
    Klink(const fVECTOR &_kernel_row, long _row_ident);

    fVECTOR kernel_row;
    L_DOUBLE *fast_kernel_row;
    long row_ident;

    Klink *next;
    Klink *prev;
};

class Kcache
{
    public:

    /*
       Constructor:

       - allocate all vectors etc.  Must be followed by a called to
         setmemsize before use (by default, the constructor sets
         everything to 0, and vectors empty).
    */

    Kcache();

    /*
       getval(caller,numi,numj):

       - returns the value for elements numi,numj.  If the element is stored
         (which can be told easily by looking for numi (or numj) in the
         row_idents vector) then will simply return the value.

         If the value is not stored, the relevant row will be constructed
         and filled out, and the relevant element returned.  If memory is
         already full before this then the oldest row (last in the linked
         list) will be overwritten.

         When constructing a new row, if possible kernel values will be
         taken from existing rows instead of re-calculating them.

         numi is used for getrow calls if possible.

       getval_statcache(caller,numi,numj):

       - return the value for elements numi,numj.  If the element is stored
         then return this.  Otherwise, calculate and return.  The cache
         itself will not be affected by this operation (unlike getval).

       getrow(caller,numi):

       - like getval, but returns (vector) row numi.

       remove(caller,num)

       - if row num is currently stored, then remove it (actually a null
         operation on kernel_rows, but row_idents will be written with a
         zero). This assumes that remove has been called because a row/
         column has been removed from the training set.  Hence not only will
         the row be removed, but the elements of the other kernel_rows
         corresponding this this will be removed (and those elements lying
         after this moved back one).

         It is assumed that remove will be called after the relevant data
         has been removed from caller.  An implicit call will be made to
         setmemsize after this function is complete.

       flush()

       - remove all kernel_rows currently in the cache by overwriting
         row_idents with zero.

       setmemsize(caller,memsize,min_rowdim)

       - this function does two things.  Firstly, rowdim is recalculated,
         and is the larger of min_rowdim and the training set size.  If this
         value has changed then the dimension of all elements of kernel_rows
         will be increased to suit (junk left as padding).  If the training
         set size has increased then (it is assumed that the new training
         data will be at the end of the training set in caller, and already
         added) the relevant kernel values will be calculated.  If the
         training set size has decreased, no action is taken (it is assumed
         in this case that this function has been called by the remove
         function, as described above).

         Secondly, numrows is calculated.  This is the number of rows of
         size rowdim which can be fitted into memsize MB of memory
         (roughly, anyhow - there will be some bookkeeping overhead which
         is ignored).  If numrows is less than the size of kernel_rows then
         the oldest rows will be removed so that the values coincide.
         Similarly, if numrows exceeds the size of kernel_rows then blank
         (filler) rows with row_idents set to zero will be added.

         It is of course assumed here that the size of the training set will
         only change in increments of 1, and that memsize and min_rowdim will
         only change by small amounts (if at all).  If this is not the case
         then things may go badly, as memory size may temporarily increase
         substantially (the exception is if this is the first call).

       reset()

       - puts cache back into state pre doing anything.
    */

    L_DOUBLE getval(SVdata *caller, long numi, long numj);
    L_DOUBLE getval_statcache(SVdata *caller, long numi, long numj);
    fVECTOR *getrow(SVdata *caller, long numi);
    void remove(SVdata *caller, long num);
    void flush(void);
    void setmemsize(SVdata *caller, long _memsize, long _min_rowdim);
    void reset(void);

    /*
       Optimisations
    */

    void enter_opt(SVdata *caller);
    void exit_opt(void);
    L_DOUBLE opt_getval_statcache(SVdata *caller, long numi, long numj);
    L_DOUBLE *opt_getrow(SVdata *caller, long numi);
    L_DOUBLE opt_getval_statcache_fast(SVdata *caller, long numi, long numj);
    L_DOUBLE *opt_getrow_fast(SVdata *caller, long numi);

    private:

    /*
       The kernel data is as follows:

       first_element = first (most recently accessed) kernel cache row.
       last_element  = last kernel cache row (may be unused).

       trainsize  = the size of the current training set.  Kept b/c it is
                    necessary to know when it has changed in the setmemsize
                    function.
       numrows    = the number of rows in the cache (including empty rows).
       rowdim     = size of each row vector in the cache.  If rowdim exceeds
                    trainsize then random padding will be present.
       memsize    = size of memory (in MB) which the cache must fit into.
       min_rowdim = minimum allowable dimension of rows in the cache.
    */

    Klink *first_element;
    Klink *last_element;

    long trainsize;
    long maxrows;
    long rowdim;
    long memsize;
    long min_rowdim;

    /*
       Optimisations
    */

    opt_kern_ptr fast_kernel_fn;
    L_DOUBLE *fast_uv;
    L_DOUBLE *fast_uc;
    L_DOUBLE *fast_covw;
    Klink **fast_lookup;
};

#endif
