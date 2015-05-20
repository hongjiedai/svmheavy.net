
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
// Support vector machine - regression class
//
// Written by: Alistair Shilton
//             Melbourne University
//



#include <iostream>
#include <math.h>
#include "svm_regress.h"
#include "svdata.h"
#include "svoptim.h"
#include "kernel.h"
#include "svflags.h"
#include "common/vector.h"
#include "common/sparsevector.h"
#include "common/matrix.h"
#include "common/factor.h"
#include "common/svdefs.h"
#include "common/outfilt.h"
#include "svm_thread.h"

#define SET_RAW_SVDATA(xopttype,xraw,xcachesize,xNexpect,xd2cbufflen)   \
{                                                                       \
    int svflags;                                                        \
                                                                        \
    CLEAR_FLAGS(svflags);                                               \
                                                                        \
    switch ( (xopttype) )                                               \
    {                                                                   \
        case SVM_OPT_ACTIVEST_FULLECACHE_FULLGCACHE_INVFACT:            \
        {                                                               \
            MAKE_OPT_ACTIVE_SET(svflags);                               \
            MAKE_E_CACHE_FULL(svflags);                                 \
            MAKE_KERN_CACHE_FULL(svflags);                              \
            MAKE_H_FACT_INVE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_ACTIVEST_SUPPECACHE_FULLGCACHE_INVFACT:            \
        {                                                               \
            MAKE_OPT_ACTIVE_SET(svflags);                               \
            MAKE_E_CACHE_SUPP(svflags);                                 \
            MAKE_KERN_CACHE_FULL(svflags);                              \
            MAKE_H_FACT_INVE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_ACTIVEST_FREEECACHE_FULLGCACHE_INVFACT:            \
        {                                                               \
            MAKE_OPT_ACTIVE_SET(svflags);                               \
            MAKE_E_CACHE_FREE(svflags);                                 \
            MAKE_KERN_CACHE_FULL(svflags);                              \
            MAKE_H_FACT_INVE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_ACTIVEST_EMPTECACHE_FULLGCACHE_INVFACT:            \
        {                                                               \
            MAKE_OPT_ACTIVE_SET(svflags);                               \
            MAKE_E_CACHE_NONE(svflags);                                 \
            MAKE_KERN_CACHE_FULL(svflags);                              \
            MAKE_H_FACT_INVE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_ACTIVEST_FULLECACHE_FREEGCACHE_INVFACT:            \
        {                                                               \
            MAKE_OPT_ACTIVE_SET(svflags);                               \
            MAKE_E_CACHE_FULL(svflags);                                 \
            MAKE_KERN_CACHE_FREE(svflags);                              \
            MAKE_H_FACT_INVE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_ACTIVEST_SUPPECACHE_FREEGCACHE_INVFACT:            \
        {                                                               \
            MAKE_OPT_ACTIVE_SET(svflags);                               \
            MAKE_E_CACHE_SUPP(svflags);                                 \
            MAKE_KERN_CACHE_FREE(svflags);                              \
            MAKE_H_FACT_INVE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_ACTIVEST_FREEECACHE_FREEGCACHE_INVFACT:            \
        {                                                               \
            MAKE_OPT_ACTIVE_SET(svflags);                               \
            MAKE_E_CACHE_FREE(svflags);                                 \
            MAKE_KERN_CACHE_FREE(svflags);                              \
            MAKE_H_FACT_INVE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_ACTIVEST_EMPTECACHE_FREEGCACHE_INVFACT:            \
        {                                                               \
            MAKE_OPT_ACTIVE_SET(svflags);                               \
            MAKE_E_CACHE_NONE(svflags);                                 \
            MAKE_KERN_CACHE_FREE(svflags);                              \
            MAKE_H_FACT_INVE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_ACTIVEST_FULLECACHE_FULLGCACHE_CHOLFACT:           \
        {                                                               \
            MAKE_OPT_ACTIVE_SET(svflags);                               \
            MAKE_E_CACHE_FULL(svflags);                                 \
            MAKE_KERN_CACHE_FULL(svflags);                              \
            MAKE_H_FACT_CHOL(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_ACTIVEST_SUPPECACHE_FULLGCACHE_CHOLFACT:           \
        {                                                               \
            MAKE_OPT_ACTIVE_SET(svflags);                               \
            MAKE_E_CACHE_SUPP(svflags);                                 \
            MAKE_KERN_CACHE_FULL(svflags);                              \
            MAKE_H_FACT_CHOL(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_ACTIVEST_FREEECACHE_FULLGCACHE_CHOLFACT:           \
        {                                                               \
            MAKE_OPT_ACTIVE_SET(svflags);                               \
            MAKE_E_CACHE_FREE(svflags);                                 \
            MAKE_KERN_CACHE_FULL(svflags);                              \
            MAKE_H_FACT_CHOL(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_ACTIVEST_EMPTECACHE_FULLGCACHE_CHOLFACT:           \
        {                                                               \
            MAKE_OPT_ACTIVE_SET(svflags);                               \
            MAKE_E_CACHE_NONE(svflags);                                 \
            MAKE_KERN_CACHE_FULL(svflags);                              \
            MAKE_H_FACT_CHOL(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_ACTIVEST_FULLECACHE_FREEGCACHE_CHOLFACT:           \
        {                                                               \
            MAKE_OPT_ACTIVE_SET(svflags);                               \
            MAKE_E_CACHE_FULL(svflags);                                 \
            MAKE_KERN_CACHE_FREE(svflags);                              \
            MAKE_H_FACT_CHOL(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_ACTIVEST_SUPPECACHE_FREEGCACHE_CHOLFACT:           \
        {                                                               \
            MAKE_OPT_ACTIVE_SET(svflags);                               \
            MAKE_E_CACHE_SUPP(svflags);                                 \
            MAKE_KERN_CACHE_FREE(svflags);                              \
            MAKE_H_FACT_CHOL(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_ACTIVEST_FREEECACHE_FREEGCACHE_CHOLFACT:           \
        {                                                               \
            MAKE_OPT_ACTIVE_SET(svflags);                               \
            MAKE_E_CACHE_FREE(svflags);                                 \
            MAKE_KERN_CACHE_FREE(svflags);                              \
            MAKE_H_FACT_CHOL(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_ACTIVEST_EMPTECACHE_FREEGCACHE_CHOLFACT:           \
        {                                                               \
            MAKE_OPT_ACTIVE_SET(svflags);                               \
            MAKE_E_CACHE_NONE(svflags);                                 \
            MAKE_KERN_CACHE_FREE(svflags);                              \
            MAKE_H_FACT_CHOL(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_PLATTSMO_FULLECACHE_FULLGCACHE:                    \
        {                                                               \
            MAKE_OPT_PLATT_SMO(svflags);                                \
            MAKE_E_CACHE_FULL(svflags);                                 \
            MAKE_KERN_CACHE_FULL(svflags);                              \
            MAKE_H_FACT_NONE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_PLATTSMO_SUPPECACHE_FULLGCACHE:                    \
        {                                                               \
            MAKE_OPT_PLATT_SMO(svflags);                                \
            MAKE_E_CACHE_SUPP(svflags);                                 \
            MAKE_KERN_CACHE_FULL(svflags);                              \
            MAKE_H_FACT_NONE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_PLATTSMO_FREEECACHE_FULLGCACHE:                    \
        {                                                               \
            MAKE_OPT_PLATT_SMO(svflags);                                \
            MAKE_E_CACHE_FREE(svflags);                                 \
            MAKE_KERN_CACHE_FULL(svflags);                              \
            MAKE_H_FACT_NONE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_PLATTSMO_EMPTECACHE_FULLGCACHE:                    \
        {                                                               \
            MAKE_OPT_PLATT_SMO(svflags);                                \
            MAKE_E_CACHE_NONE(svflags);                                 \
            MAKE_KERN_CACHE_FULL(svflags);                              \
            MAKE_H_FACT_NONE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_PLATTSMO_FULLECACHE_FREEGCACHE:                    \
        {                                                               \
            MAKE_OPT_PLATT_SMO(svflags);                                \
            MAKE_E_CACHE_FULL(svflags);                                 \
            MAKE_KERN_CACHE_FREE(svflags);                              \
            MAKE_H_FACT_NONE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_PLATTSMO_SUPPECACHE_FREEGCACHE:                    \
        {                                                               \
            MAKE_OPT_PLATT_SMO(svflags);                                \
            MAKE_E_CACHE_SUPP(svflags);                                 \
            MAKE_KERN_CACHE_FREE(svflags);                              \
            MAKE_H_FACT_NONE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_PLATTSMO_FREEECACHE_FREEGCACHE:                    \
        {                                                               \
            MAKE_OPT_PLATT_SMO(svflags);                                \
            MAKE_E_CACHE_FREE(svflags);                                 \
            MAKE_KERN_CACHE_FREE(svflags);                              \
            MAKE_H_FACT_NONE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_PLATTSMO_EMPTECACHE_FREEGCACHE:                    \
        {                                                               \
            MAKE_OPT_PLATT_SMO(svflags);                                \
            MAKE_E_CACHE_NONE(svflags);                                 \
            MAKE_KERN_CACHE_FREE(svflags);                              \
            MAKE_H_FACT_NONE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_PLATTSMO_FULLECACHE_DYNAGCACHE:                    \
        {                                                               \
            MAKE_OPT_PLATT_SMO(svflags);                                \
            MAKE_E_CACHE_FULL(svflags);                                 \
            MAKE_KERN_CACHE_DYNA(svflags);                              \
            MAKE_H_FACT_NONE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_PLATTSMO_SUPPECACHE_DYNAGCACHE:                    \
        {                                                               \
            MAKE_OPT_PLATT_SMO(svflags);                                \
            MAKE_E_CACHE_SUPP(svflags);                                 \
            MAKE_KERN_CACHE_DYNA(svflags);                              \
            MAKE_H_FACT_NONE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_PLATTSMO_FREEECACHE_DYNAGCACHE:                    \
        {                                                               \
            MAKE_OPT_PLATT_SMO(svflags);                                \
            MAKE_E_CACHE_FREE(svflags);                                 \
            MAKE_KERN_CACHE_DYNA(svflags);                              \
            MAKE_H_FACT_NONE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_PLATTSMO_EMPTECACHE_DYNAGCACHE:                    \
        {                                                               \
            MAKE_OPT_PLATT_SMO(svflags);                                \
            MAKE_E_CACHE_NONE(svflags);                                 \
            MAKE_KERN_CACHE_DYNA(svflags);                              \
            MAKE_H_FACT_NONE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_PLATTSMO_FULLECACHE_EMPTGCACHE:                    \
        {                                                               \
            MAKE_OPT_PLATT_SMO(svflags);                                \
            MAKE_E_CACHE_FULL(svflags);                                 \
            MAKE_KERN_CACHE_NONE(svflags);                              \
            MAKE_H_FACT_NONE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_PLATTSMO_SUPPECACHE_EMPTGCACHE:                    \
        {                                                               \
            MAKE_OPT_PLATT_SMO(svflags);                                \
            MAKE_E_CACHE_SUPP(svflags);                                 \
            MAKE_KERN_CACHE_NONE(svflags);                              \
            MAKE_H_FACT_NONE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_PLATTSMO_FREEECACHE_EMPTGCACHE:                    \
        {                                                               \
            MAKE_OPT_PLATT_SMO(svflags);                                \
            MAKE_E_CACHE_FREE(svflags);                                 \
            MAKE_KERN_CACHE_NONE(svflags);                              \
            MAKE_H_FACT_NONE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_PLATTSMO_EMPTECACHE_EMPTGCACHE:                    \
        {                                                               \
            MAKE_OPT_PLATT_SMO(svflags);                                \
            MAKE_E_CACHE_NONE(svflags);                                 \
            MAKE_KERN_CACHE_NONE(svflags);                              \
            MAKE_H_FACT_NONE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_D2CALGOR_FULLECACHE_FULLGCACHE:                    \
        {                                                               \
            MAKE_OPT_DANIEL_D2C(svflags);                               \
            MAKE_E_CACHE_FULL(svflags);                                 \
            MAKE_KERN_CACHE_FULL(svflags);                              \
            MAKE_H_FACT_NONE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_D2CALGOR_SUPPECACHE_FULLGCACHE:                    \
        {                                                               \
            MAKE_OPT_DANIEL_D2C(svflags);                               \
            MAKE_E_CACHE_SUPP(svflags);                                 \
            MAKE_KERN_CACHE_FULL(svflags);                              \
            MAKE_H_FACT_NONE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_D2CALGOR_FREEECACHE_FULLGCACHE:                    \
        {                                                               \
            MAKE_OPT_DANIEL_D2C(svflags);                               \
            MAKE_E_CACHE_FREE(svflags);                                 \
            MAKE_KERN_CACHE_FULL(svflags);                              \
            MAKE_H_FACT_NONE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_D2CALGOR_EMPTECACHE_FULLGCACHE:                    \
        {                                                               \
            MAKE_OPT_DANIEL_D2C(svflags);                               \
            MAKE_E_CACHE_NONE(svflags);                                 \
            MAKE_KERN_CACHE_FULL(svflags);                              \
            MAKE_H_FACT_NONE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_D2CALGOR_FULLECACHE_FREEGCACHE:                    \
        {                                                               \
            MAKE_OPT_DANIEL_D2C(svflags);                               \
            MAKE_E_CACHE_FULL(svflags);                                 \
            MAKE_KERN_CACHE_FREE(svflags);                              \
            MAKE_H_FACT_NONE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_D2CALGOR_SUPPECACHE_FREEGCACHE:                    \
        {                                                               \
            MAKE_OPT_DANIEL_D2C(svflags);                               \
            MAKE_E_CACHE_SUPP(svflags);                                 \
            MAKE_KERN_CACHE_FREE(svflags);                              \
            MAKE_H_FACT_NONE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_D2CALGOR_FREEECACHE_FREEGCACHE:                    \
        {                                                               \
            MAKE_OPT_DANIEL_D2C(svflags);                               \
            MAKE_E_CACHE_FREE(svflags);                                 \
            MAKE_KERN_CACHE_FREE(svflags);                              \
            MAKE_H_FACT_NONE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_D2CALGOR_EMPTECACHE_FREEGCACHE:                    \
        {                                                               \
            MAKE_OPT_DANIEL_D2C(svflags);                               \
            MAKE_E_CACHE_NONE(svflags);                                 \
            MAKE_KERN_CACHE_FREE(svflags);                              \
            MAKE_H_FACT_NONE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_D2CALGOR_FULLECACHE_DYNAGCACHE:                    \
        {                                                               \
            MAKE_OPT_DANIEL_D2C(svflags);                               \
            MAKE_E_CACHE_FULL(svflags);                                 \
            MAKE_KERN_CACHE_DYNA(svflags);                              \
            MAKE_H_FACT_NONE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_D2CALGOR_SUPPECACHE_DYNAGCACHE:                    \
        {                                                               \
            MAKE_OPT_DANIEL_D2C(svflags);                               \
            MAKE_E_CACHE_SUPP(svflags);                                 \
            MAKE_KERN_CACHE_DYNA(svflags);                              \
            MAKE_H_FACT_NONE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_D2CALGOR_FREEECACHE_DYNAGCACHE:                    \
        {                                                               \
            MAKE_OPT_DANIEL_D2C(svflags);                               \
            MAKE_E_CACHE_FREE(svflags);                                 \
            MAKE_KERN_CACHE_DYNA(svflags);                              \
            MAKE_H_FACT_NONE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_D2CALGOR_EMPTECACHE_DYNAGCACHE:                    \
        {                                                               \
            MAKE_OPT_DANIEL_D2C(svflags);                               \
            MAKE_E_CACHE_NONE(svflags);                                 \
            MAKE_KERN_CACHE_DYNA(svflags);                              \
            MAKE_H_FACT_NONE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_D2CALGOR_FULLECACHE_EMPTGCACHE:                    \
        {                                                               \
            MAKE_OPT_DANIEL_D2C(svflags);                               \
            MAKE_E_CACHE_FULL(svflags);                                 \
            MAKE_KERN_CACHE_NONE(svflags);                              \
            MAKE_H_FACT_NONE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_D2CALGOR_SUPPECACHE_EMPTGCACHE:                    \
        {                                                               \
            MAKE_OPT_DANIEL_D2C(svflags);                               \
            MAKE_E_CACHE_SUPP(svflags);                                 \
            MAKE_KERN_CACHE_NONE(svflags);                              \
            MAKE_H_FACT_NONE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_D2CALGOR_FREEECACHE_EMPTGCACHE:                    \
        {                                                               \
            MAKE_OPT_DANIEL_D2C(svflags);                               \
            MAKE_E_CACHE_FREE(svflags);                                 \
            MAKE_KERN_CACHE_NONE(svflags);                              \
            MAKE_H_FACT_NONE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SVM_OPT_D2CALGOR_EMPTECACHE_EMPTGCACHE:                    \
        {                                                               \
            MAKE_OPT_DANIEL_D2C(svflags);                               \
            MAKE_E_CACHE_NONE(svflags);                                 \
            MAKE_KERN_CACHE_NONE(svflags);                              \
            MAKE_H_FACT_NONE(svflags);                                  \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        default:                                                        \
        {                                                               \
            L_THROW(0);                                                 \
                                                                        \
            break;                                                      \
        }                                                               \
    }                                                                   \
                                                                        \
                                                                        \
                                                                        \
    if ( fix_bias )                                                     \
    {                                                                   \
        REPNEW((xraw),SVdata(_kern,default_bias,svflags,xcachesize,xNexpect,xd2cbufflen)); \
    }                                                                   \
                                                                        \
    else                                                                \
    {                                                                   \
        REPNEW((xraw),SVdata(_kern,svflags,xcachesize,xNexpect,xd2cbufflen)); \
    }                                                                   \
}

#define INIT_ASYNCH                             \
{                                               \
    async_exit_flag  = 0;                       \
    back_opt_flag    = 0;                       \
    loop_count       = 0;                       \
    loop_count_accum = 0;                       \
                                                \
    J = 0.0;                                    \
}

#define LOCK_FN_START           \
int back_opt_flag_temp = 0;     \
{                               \
    if ( is_back_opt() )        \
    {                           \
        back_opt_off();         \
        back_opt_flag_temp = 1; \
    }                           \
}
#define LOCK_FN_END             \
{                               \
    if ( back_opt_flag_temp )   \
    {                           \
        back_opt_on();          \
    }                           \
}

#define EXLOCK_FN_START(what)   \
int back_opt_flag_temp = 0;     \
{                               \
    if ( what.is_back_opt() )   \
    {                           \
        what.back_opt_off();    \
        back_opt_flag_temp = 1; \
    }                           \
}
#define EXLOCK_FN_END(what)     \
{                               \
    if ( back_opt_flag_temp )   \
    {                           \
        what.back_opt_on();     \
    }                           \
}



//
// Stream IO overloading
//

std::ostream &operator<<(std::ostream &output, SVM_regress &dumpee)
{
    EXLOCK_FN_START(dumpee);

    THROW_ASSERT(dumpee.raw != NULL);

    output << "r\n";

    output << "===========SVM regression data===========\n\n";

    output << "Optim type:  " << dumpee.opt_type    << "\n";
    output << "Risk type:   " << dumpee.risk_type   << "\n";
    output << "Tube shrink: " << dumpee.tube_shrink << "\n\n";

    output << "C/N: " << fixnum(dumpee.C)             << "\n";
    output << "C+:  " << fixnum(dumpee.C_plus_scale)  << "\n";
    output << "C-:  " << fixnum(dumpee.C_minus_scale) << "\n\n";

    output << "E:  " << fixnum(dumpee.E)             << "\n";
    output << "E+: " << fixnum(dumpee.E_plus_scale)  << "\n";
    output << "E-: " << fixnum(dumpee.E_minus_scale) << "\n\n";

    output << "nu:  " << fixnum(dumpee.nu) << "\n\n";

    output << "sol_tol: " << dumpee.sol_tol << "\n";
    output << "epochs:  " << dumpee.epochs  << "\n\n";

    SVdata_summary temp(*(dumpee.raw));

    output << "SVM data:\n\n" << temp << "\n\n";

    EXLOCK_FN_END(dumpee);

    return output;
}

std::istream &operator>>(std::istream &input, SVM_regress &source)
{
    EXLOCK_FN_START(source);

    wait_dummy howzat;
    L_DOUBLE default_bias;

    std::ostream *where_to;
    int echo_level;                                                       
    int fix_bias;



    where_to   = source.where_to;                                         
    echo_level = source.echo_level;                                       
                                                                       
    source.prime_io(NULL,0);                                                
                                                                       
    if ( where_to == NULL )                                               
    {                                                                     
        input >> howzat; input >> source.opt_type;
        input >> howzat; input >> source.risk_type;
        input >> howzat; input >> source.tube_shrink;

        input >> howzat; input >> source.C;
        input >> howzat; input >> source.C_plus_scale;
        input >> howzat; input >> source.C_minus_scale;

        input >> howzat; input >> source.E;
        input >> howzat; input >> source.E_plus_scale;
        input >> howzat; input >> source.E_minus_scale;

        input >> howzat; input >> source.nu;

        input >> howzat; input >> source.sol_tol;
        input >> howzat; input >> source.epochs;

        if ( source.raw != NULL )
        {
            REPDEL((source.raw));
        }

        REPNEW(source.raw,SVdata(DEFAULT_SVFLAGS,1,DEFAULT_TRAININGSIZE,DEFAULT_D2C_BUFF_LEN));

        input >> howzat; input >> *(source.raw);
    }
                                                                       
    else                                                                  
    {
        kernel_dicto _kern;

        *where_to << "Optimisation type (244 default): "; input >> source.opt_type; if ( echo_level ) { *where_to << source.opt_type << "\n"; }
        *where_to << "Fixed bias (0 = no, 1 = yes): ";    input >> fix_bias;        if ( echo_level ) { *where_to << fix_bias        << "\n"; }

        if ( fix_bias )
        {
            *where_to << "Bias: "; input >> default_bias; if ( echo_level ) { *where_to << default_bias << "\n"; }
        }

        else
        {
            default_bias = 0.0;
        }

        *where_to << "Empirical risk type (0 = linear, 1 = quadratic): "; input >> source.risk_type; if ( echo_level ) { *where_to << source.risk_type << "\n"; }
        *where_to << "Tube shrinking (0 = no, 1 = quadratic): "; input >> source.tube_shrink; if ( echo_level ) { *where_to << source.tube_shrink << "\n"; }

        if ( source.tube_shrink )
        {
            *where_to << "nu: "; input >> source.nu; if ( echo_level ) { *where_to << source.nu << "\n"; }
        }

        else
        {
            source.nu = 1.0;
        }

        *where_to << "\n";
        *where_to << "C/N: "; input >> source.C;             if ( echo_level ) { *where_to << source.C             << "\n"; }
        *where_to << "C+:  "; input >> source.C_plus_scale;  if ( echo_level ) { *where_to << source.C_plus_scale  << "\n"; }
        *where_to << "C-:  "; input >> source.C_minus_scale; if ( echo_level ) { *where_to << source.C_minus_scale << "\n"; }

        *where_to << "\n";
        *where_to << "E:  "; input >> source.E;             if ( echo_level ) { *where_to << source.E             << "\n"; }
        *where_to << "E+: "; input >> source.E_plus_scale;  if ( echo_level ) { *where_to << source.E_plus_scale  << "\n"; }
        *where_to << "E-: "; input >> source.E_minus_scale; if ( echo_level ) { *where_to << source.E_minus_scale << "\n"; }

        *where_to << "\n";
        *where_to << "solution tolerance: "; input >> source.sol_tol; if ( echo_level ) { *where_to << source.sol_tol << "\n"; }
        *where_to << "max iterations: "; input >> source.epochs;      if ( echo_level ) { *where_to << source.epochs  << "\n"; }

        *where_to << "\n";
        *where_to << "Kernel:\n\n";

        _kern.prime_io(where_to,echo_level);

        input >> _kern;

        if ( source.raw != NULL )
        {
            REPDEL((source.raw));
        }

        SET_RAW_SVDATA(source.opt_type,source.raw,DEFAULT_CACHE_SIZE,DEFAULT_TRAININGSIZE,DEFAULT_D2C_BUFF_LEN);
    }                                                                     

    EXLOCK_FN_END(source);

    return input;
}




//
// Constructor:
//

SVM_regress::SVM_regress()
{
    INIT_ASYNCH;

    raw = NULL;

    return;
}

SVM_regress::SVM_regress(const SVM_regress &source)
{
    INIT_ASYNCH;

    opt_type    = source.opt_type;
    risk_type   = source.risk_type;
    tube_shrink = source.tube_shrink;

    C             = source.C;
    C_plus_scale  = source.C_plus_scale;
    C_minus_scale = source.C_minus_scale;

    E             = source.E;
    E_plus_scale  = source.E_plus_scale;
    E_minus_scale = source.E_minus_scale;

    nu   = source.nu;

    sol_tol = source.sol_tol;
    epochs  = source.epochs;

    REPNEW(raw,SVdata(*(source.raw)));

    return;
}

SVM_regress::SVM_regress(const kernel_dicto &_kern,
                         int _risk_type,
                         int _tube_shrink,
                         L_DOUBLE _cee,
                         L_DOUBLE _C_plus_scale,
                         L_DOUBLE _C_minus_scale,
                         L_DOUBLE _E,
                         L_DOUBLE _E_plus_scale,
                         L_DOUBLE _E_minus_scale,
                         L_DOUBLE _nu,
                         int fix_bias,
                         L_DOUBLE default_bias,
                         int _opt_type,
                         long _epochs,
                         double _sol_tol,
                         long cachesize,
                         long Nexpect,
                         long d2cbufflen)
{
    INIT_ASYNCH;

    opt_type    = _opt_type;
    risk_type   = _risk_type;
    tube_shrink = _tube_shrink;

    C             = _cee;
    C_plus_scale  = _C_plus_scale;
    C_minus_scale = _C_minus_scale;

    E             = _E;
    E_plus_scale  = _E_plus_scale;
    E_minus_scale = _E_minus_scale;

    nu   = _nu;

    sol_tol = _sol_tol;
    epochs  = _epochs;

    SET_RAW_SVDATA(opt_type,raw,cachesize,Nexpect,d2cbufflen);

    return;
}

SVM_regress &SVM_regress::operator=(const SVM_regress &source)
{
    INIT_ASYNCH;

    if ( raw != NULL )
    {
        REPDEL(raw);
    }

    opt_type    = source.opt_type;
    risk_type   = source.risk_type;
    tube_shrink = source.tube_shrink;

    C             = source.C;
    C_plus_scale  = source.C_plus_scale;
    C_minus_scale = source.C_minus_scale;

    E             = source.E;
    E_plus_scale  = source.E_plus_scale;
    E_minus_scale = source.E_minus_scale;

    nu   = source.nu;

    sol_tol = source.sol_tol;
    epochs  = source.epochs;

    THROW_ASSERT( (source.raw) != NULL );

    REPNEW(raw,SVdata(*(source.raw)));

    return *this;
}




//
// Destructor:
//

SVM_regress::~SVM_regress()
{
    back_opt_off();

    if ( raw != NULL )
    {
        REPDEL(raw);
    }

    return;
}




//
// Asynchronous operation functions.
//

int SVM_regress::back_opt_on(void)
{
    int result = 0;

    if ( !back_opt_flag )
    {
        if ( !( result = svmthread_create(background_opt_callback_regress,(void *) this) ) )
        {
            loop_count    = 1;
            back_opt_flag = 1;
        }
    }

    return result;
}

void SVM_regress::back_opt_off(void)
{
    if ( back_opt_flag )
    {
        while ( loop_count )
        {
            async_exit_flag = 1;
        }

        async_exit_flag = 0;
        back_opt_flag   = 0;
    }

    return;
}

int SVM_regress::is_back_opt(void) const
{
    return back_opt_flag;
}

unsigned long SVM_regress::get_iter(void) const
{
    return loop_count;
}

unsigned long SVM_regress::get_iter_accum(void) const
{
    return loop_count_accum;
}

void SVM_regress::set_iter_accum(unsigned long itval)
{
    loop_count_accum = itval;

    return;
}

L_DOUBLE SVM_regress::get_J(void) const
{
    return J;
}

void SVM_regress::set_J(L_DOUBLE Jval)
{
    J = Jval;

    return;
}




//
// Training set modification functions.
//

void SVM_regress::add_point(L_DOUBLE z, fVECTOR &x, L_DOUBLE t, L_DOUBLE t_star, L_DOUBLE eps, L_DOUBLE eps_star, long id)
{
    LOCK_FN_START;

    THROW_ASSERT(raw != NULL);

    if ( tube_shrink )
    {
        if ( risk_type )
        {
            raw->addpoint(3,id,x,z,0.0,0.0,MAX_UPPER_BOUND,-MAX_UPPER_BOUND,1.0/(C*C_plus_scale*t),1.0/(C*C_minus_scale*t_star),eps/(C*nu),eps_star/(C*nu));
        }
                                      
        else                          
        {
            raw->addpoint(3,id,x,z,0.0,0.0,C*C_plus_scale*t,-C*C_minus_scale*t_star,0.0,0.0,eps/(C*nu),eps_star/(C*nu));
        }
    }

    else
    {
        if ( risk_type )
        {
            raw->addpoint(3,id,x,z,E*E_plus_scale*eps,E*E_minus_scale*eps_star,MAX_UPPER_BOUND,-MAX_UPPER_BOUND,1.0/(C*C_plus_scale*t),1.0/(C*C_minus_scale*t_star),0.0,0.0);
        }

        else
        {
            raw->addpoint(3,id,x,z,E*E_plus_scale*eps,E*E_minus_scale*eps_star,C*C_plus_scale*t,-C*C_minus_scale*t_star,0.0,0.0,0.0,0.0);
        }
    }

    LOCK_FN_END;

    return;
}

void SVM_regress::add_points(fVECTOR &z, f_fVECTOR &x, fVECTOR &t, fVECTOR &t_star, fVECTOR &eps, fVECTOR &eps_star, iVECTOR &id)
{
    long i;
    long numpoints;

    numpoints = z.get_effective_size();

    THROW_ASSERT(numpoints >= 0);
    THROW_ASSERT(numpoints == x.get_real_size());
    THROW_ASSERT(numpoints == t.get_effective_size());
    THROW_ASSERT(numpoints == t_star.get_effective_size());
    THROW_ASSERT(numpoints == eps.get_effective_size());
    THROW_ASSERT(numpoints == eps_star.get_effective_size());
    THROW_ASSERT(numpoints == id.get_effective_size());

    if ( numpoints >= 1 )
    {
        for ( i = 1 ; i <= numpoints ; i++ )
        {
            add_point(z[i-1],x[i-1],t[i-1],t_star[i-1],eps[i-1],eps_star[i-1],id[i-1]);
        }
    }

    return;
}

void SVM_regress::add_point_lower_bound(L_DOUBLE z, fVECTOR &x, L_DOUBLE t, L_DOUBLE eps, long id)
{
    LOCK_FN_START;

    THROW_ASSERT(raw != NULL);

    if ( tube_shrink )
    {
        if ( risk_type )
        {
            raw->addpoint(2,id,x,z,0.0,0.0,MAX_UPPER_BOUND,0.0,1.0/(C*C_plus_scale*t),0.0,eps/(C*nu),0.0);
        }
                                      
        else                          
        {
            raw->addpoint(2,id,x,z,0.0,0.0,C*C_plus_scale*t,0.0,0.0,0.0,eps/(C*nu),0.0);
        }
    }

    else
    {
        if ( risk_type )
        {
            raw->addpoint(2,id,x,z,E*E_plus_scale*eps,0.0,MAX_UPPER_BOUND,0.0,1.0/(C*C_plus_scale*t),0.0,0.0,0.0);
        }

        else
        {
            raw->addpoint(2,id,x,z,E*E_plus_scale*eps,0.0,C*C_plus_scale*t,0.0,0.0,0.0,0.0,0.0);
        }
    }

    LOCK_FN_END;

    return;
}

void SVM_regress::add_points_lower_bound(fVECTOR &z, f_fVECTOR &x, fVECTOR &t, fVECTOR &eps, iVECTOR &id)
{
    long i;
    long numpoints;

    numpoints = z.get_effective_size();

    THROW_ASSERT(numpoints >= 0);
    THROW_ASSERT(numpoints == x.get_real_size());
    THROW_ASSERT(numpoints == t.get_effective_size());
    THROW_ASSERT(numpoints == eps.get_effective_size());
    THROW_ASSERT(numpoints == id.get_effective_size());

    if ( numpoints >= 1 )
    {
        for ( i = 1 ; i <= numpoints ; i++ )
        {
            add_point_lower_bound(z[i-1],x[i-1],t[i-1],eps[i-1],id[i-1]);
        }
    }

    return;
}

void SVM_regress::add_point_upper_bound(L_DOUBLE z, fVECTOR &x, L_DOUBLE t_star, L_DOUBLE eps_star, long id)
{
    LOCK_FN_START;

    THROW_ASSERT(raw != NULL);

    if ( tube_shrink )
    {
        if ( risk_type )
        {
            raw->addpoint(1,id,x,z,0.0,0.0,0.0,-MAX_UPPER_BOUND,0.0,1.0/(C*C_minus_scale*t_star),0.0,eps_star/(C*nu));
        }
                                      
        else                          
        {
            raw->addpoint(1,id,x,z,0.0,0.0,0.0,C*C_minus_scale*t_star,0.0,0.0,0.0,eps_star/(C*nu));
        }
    }

    else
    {
        if ( risk_type )
        {
            raw->addpoint(1,id,x,z,0.0,E*E_minus_scale*eps_star,0.0,-MAX_UPPER_BOUND,0.0,1.0/(C*C_minus_scale*t_star),0.0,0.0);
        }

        else
        {
            raw->addpoint(1,id,x,z,0.0,E*E_minus_scale*eps_star,0.0,-C*C_minus_scale*t_star,0.0,0.0,0.0,0.0);
        }
    }

    LOCK_FN_END;

    return;
}

void SVM_regress::add_points_upper_bound(fVECTOR &z, f_fVECTOR &x, fVECTOR &t_star, fVECTOR &eps_star, iVECTOR &id)
{
    long i;
    long numpoints;

    numpoints = z.get_effective_size();

    THROW_ASSERT(numpoints >= 0);
    THROW_ASSERT(numpoints == x.get_real_size());
    THROW_ASSERT(numpoints == t_star.get_effective_size());
    THROW_ASSERT(numpoints == eps_star.get_effective_size());
    THROW_ASSERT(numpoints == id.get_effective_size());

    if ( numpoints >= 1 )
    {
        for ( i = 1 ; i <= numpoints ; i++ )
        {
            add_point_upper_bound(z[i-1],x[i-1],t_star[i-1],eps_star[i-1],id[i-1]);
        }
    }

    return;
}

void SVM_regress::del_point(long id)
{
    LOCK_FN_START;

    THROW_ASSERT(raw != NULL);

    long which;

    which = raw->find_marker(id);

    if ( which )
    {
        raw->delpoint(which);
    }

    LOCK_FN_END;

    return;
}

void SVM_regress::del_points(iVECTOR id)
{
    long i;
    long numpoints;

    numpoints = id.get_effective_size();

    THROW_ASSERT(numpoints >= 0);

    if ( numpoints >= 1 )
    {
        for ( i = 1 ; i <= numpoints ; i++ )
        {
            del_point(id[i-1]);
        }
    }

    return;
}




//
// Grid helper functions.
//

void SVM_regress::set_z(long i, L_DOUBLE z)
{
    THROW_ASSERT( i >= 0 );
    THROW_ASSERT( i <= get_N() );

    raw->set_z(i,z);

    return;
}




//
// Training parameter modification functions:
//

void SVM_regress::set_C(L_DOUBLE C_new)
{
    scale_C(C_new/C);

    return;
}

void SVM_regress::set_E(L_DOUBLE E_new)
{
    scale_E(E_new/E);

    return;
}

void SVM_regress::set_nu(L_DOUBLE nu_new)
{
    scale_nu(nu_new/nu);

    return;
}

void SVM_regress::set_C_plus_scale(L_DOUBLE C_new)
{
    scale_C_plus_scale(C_new/C_plus_scale);

    return;
}

void SVM_regress::set_C_minus_scale(L_DOUBLE C_new)
{
    scale_C_minus_scale(C_new/C_minus_scale);

    return;
}

void SVM_regress::set_E_plus_scale(L_DOUBLE E_new)
{
    scale_E_plus_scale(E_new/E_plus_scale);

    return;
}

void SVM_regress::set_E_minus_scale(L_DOUBLE E_new)
{
    scale_E_minus_scale(E_new/E_minus_scale);

    return;
}

void SVM_regress::scale_C(L_DOUBLE C_scale)
{
    LOCK_FN_START;

    THROW_ASSERT(raw != NULL);

    if ( risk_type )
    {
        raw->scale_gamma_both(1.0/C_scale);
    }

    else
    {
        raw->scale_hv(C_scale);
    }

    if ( tube_shrink )
    {
        raw->scale_mu_both(1.0/C_scale);
    }

    C *= C_scale;

    LOCK_FN_END;

    return;
}

void SVM_regress::scale_E(L_DOUBLE E_scale)
{
    LOCK_FN_START;

    THROW_ASSERT(raw != NULL);

    raw->scale_rho_both(E_scale);

    E *= E_scale;

    LOCK_FN_END;

    return;
}

void SVM_regress::scale_nu(L_DOUBLE nu_scale)
{
    LOCK_FN_START;

    THROW_ASSERT(raw != NULL);

    raw->scale_mu_both(1.0/nu_scale);

    nu *= nu_scale;

    LOCK_FN_END;

    return;
}

void SVM_regress::scale_C_plus_scale(L_DOUBLE C_scale)
{
    LOCK_FN_START;

    THROW_ASSERT(raw != NULL);

    if ( risk_type )
    {
        raw->scale_gamma(1.0/C_scale);
    }

    else
    {
        raw->scale_h(C_scale);
    }

    if ( tube_shrink )
    {
        raw->scale_mu(1.0/C_scale);
    }

    C_plus_scale *= C_scale;

    LOCK_FN_END;

    return;
}

void SVM_regress::scale_C_minus_scale(L_DOUBLE C_scale)
{
    LOCK_FN_START;

    THROW_ASSERT(raw != NULL);

    if ( risk_type )
    {
        raw->scale_gamma_star(1.0/C_scale);
    }

    else
    {
        raw->scale_v(C_scale);
    }

    if ( tube_shrink )
    {
        raw->scale_mu_star(1.0/C_scale);
    }

    C_minus_scale *= C_scale;

    LOCK_FN_END;

    return;
}

void SVM_regress::scale_E_plus_scale(L_DOUBLE E_scale)
{
    LOCK_FN_START;

    THROW_ASSERT(raw != NULL);

    raw->scale_rho(E_scale);

    E_plus_scale *= E_scale;

    LOCK_FN_END;

    return;
}

void SVM_regress::scale_E_minus_scale(L_DOUBLE E_scale)
{
    LOCK_FN_START;

    THROW_ASSERT(raw != NULL);

    raw->scale_rho_star(E_scale);

    E_minus_scale *= E_scale;

    LOCK_FN_END;

    return;
}

void SVM_regress::set_kernel(kernel_dicto &kerndict)
{
    LOCK_FN_START;

    THROW_ASSERT(raw != NULL);

    raw->set_kernel(kerndict);

    LOCK_FN_END;

    return;
}

void SVM_regress::set_epochs(long _epochs)
{
    epochs = _epochs;

    return;
}




//
// Information functions
//

L_DOUBLE SVM_regress::get_C(void)
{
    THROW_ASSERT(raw != NULL);

    return C;
}

L_DOUBLE SVM_regress::get_C_plus_scale(void)
{
    THROW_ASSERT(raw != NULL);

    return C_plus_scale;
}

L_DOUBLE SVM_regress::get_C_minus_scale(void)
{
    THROW_ASSERT(raw != NULL);

    return C_minus_scale;
}

L_DOUBLE SVM_regress::get_E(void)
{
    THROW_ASSERT(raw != NULL);

    return E;
}

L_DOUBLE SVM_regress::get_E_plus_scale(void)
{
    THROW_ASSERT(raw != NULL);

    return E_minus_scale;
}

L_DOUBLE SVM_regress::get_E_minus_scale(void)
{
    THROW_ASSERT(raw != NULL);

    return E_plus_scale;
}

L_DOUBLE SVM_regress::get_nu(void)
{
    THROW_ASSERT(raw != NULL);

    return nu;
}

kernel_dicto &SVM_regress::get_kernel(void)
{
    THROW_ASSERT(raw != NULL);

    return raw->K;
}

fVECTOR &SVM_regress::get_x(long i)
{
    THROW_ASSERT(raw != NULL);
    THROW_ASSERT( i >= 0 );
    THROW_ASSERT( i <= get_N() );

    return (raw->x)[i-1];
}




//
// Optimise: Optimise the SVM.
//

unsigned long SVM_regress::optimise(void)
{
    LOCK_FN_START;

    unsigned long result;

    result = optimise_unsafe();

    LOCK_FN_END;

    return result;
}

unsigned long SVM_regress::optimise_unsafe(void)
{
    THROW_ASSERT(raw != NULL);

    async_exit_flag = 0;

    return solve(*raw,sol_tol,epochs,&async_exit_flag,&loop_count,&loop_count_accum,&J);
}

unsigned long SVM_regress::optimise(L_DOUBLE temp_sol_tol)
{
    LOCK_FN_START;

    THROW_ASSERT(raw != NULL);

    async_exit_flag = 0;

    unsigned long result;
    result = solve(*raw,temp_sol_tol,epochs,&async_exit_flag,&loop_count,&loop_count_accum,&J);

    LOCK_FN_END;

    return result;
}

unsigned long SVM_regress::optimise_from_zero(L_DOUBLE temp_sol_tol, long MaxItersPerGeneration)
{
    LOCK_FN_START;

    long i;
    long j;
    long k;

    THROW_ASSERT(raw != NULL);

    i = get_N_Z()+get_N_L()+get_N_U()+1;
    j = get_N();

    if ( i <= j )
    {
        //
        // Observation: The active set found by an SVM will tend to be
        //              much the same, regardless of the kernel used
        //              (within reason).  However, the actual alpha
        //              corresponding to this is not constant, so it
        //              seems reasonable to retain the support vectors
        //              and throw away alpha.
        //

        fVECTOR alpha_pivot('n',j-i+1);

        for ( k = i ; k <= j ; k++ )
        {
            alpha_pivot[k-i] = -(raw->get_alpha_pivot(k));
        }

        if ( raw->fix_bias )
        {
            raw->step_general_pivotx(i,j,alpha_pivot,0.0);
        }

        else
        {
            raw->step_general_pivotx(i,j,alpha_pivot,-(raw->b));
        }
    }

    async_exit_flag = 0;

    unsigned long result;

    result = solve(*raw,temp_sol_tol,MaxItersPerGeneration,&async_exit_flag,&loop_count,&loop_count_accum,&J);

    LOCK_FN_END;

    return result;
}

void SVM_regress::removeNonSupports(void)
{
    LOCK_FN_START;

    long i,j;

    if ( ( j = get_N() ) > 0 )
    {
        for ( i = j ; i >= 1 ; i-- )
        {
            if ( (raw->tau)[i-1] == 0 )
            {
                raw->delpoint(i);
            }
        }
    }

    LOCK_FN_END;

    return;
}




//
// test_point returns g(x)
// cross_test_point returns dg(x)/dx
//

L_DOUBLE SVM_regress::test_point(const fVECTOR &x)
{
    LOCK_FN_START;

    THROW_ASSERT(raw != NULL);

    L_DOUBLE result;
    fVECTOR loc_x(x);
    result = raw->calc_g(loc_x);

    LOCK_FN_END;

    return result;
}

L_DOUBLE SVM_regress::cross_test_point(const fVECTOR &x)
{
    LOCK_FN_START;

    THROW_ASSERT(raw != NULL);

    L_DOUBLE result;
    result = raw->calc_g_grad(x);

    LOCK_FN_END;

    return result;
}
















//
// Get information.
//

long SVM_regress::get_N_S(void)
{
    THROW_ASSERT(raw != NULL);

    return (raw->N_L)+(raw->N_U)+(raw->N_F);
}

long SVM_regress::get_N(void)
{
    THROW_ASSERT(raw != NULL);

    return raw->N;
}

long SVM_regress::get_N_Z(void)
{
    THROW_ASSERT(raw != NULL);

    return raw->N_Z;
}

long SVM_regress::get_N_L(void)
{
    THROW_ASSERT(raw != NULL);

    return raw->N_L;
}

long SVM_regress::get_N_U(void)
{
    THROW_ASSERT(raw != NULL);

    return raw->N_U;
}

long SVM_regress::get_N_F(void)
{
    THROW_ASSERT(raw != NULL);

    return raw->N_F;
}

long SVM_regress::get_N_FP(void)
{
    THROW_ASSERT(raw != NULL);

    return raw->N_FP;
}

long SVM_regress::get_N_FN(void)
{
    THROW_ASSERT(raw != NULL);

    return raw->N_FN;
}

long SVM_regress::get_max_id(void)
{
    THROW_ASSERT(raw != NULL);

    long i,j,k;

    j = (raw->ident).get_effective_size();
    k = 0;

    if ( j > 0 )
    {
        k = (raw->ident)[0];

        for ( i = 1 ; i <= j ; i++ )
        {
            if ( (raw->ident)[i-1] >= k )
            {
                k = (raw->ident)[i-1];
            }
        }
    }

    return k;
}

fVECTOR SVM_regress::get_alpha(void)
{
    LOCK_FN_START;

    fVECTOR result('x',get_N());
    long i,N;

    N = get_N();

    if ( N > 0 )
    {
        for ( i = 0 ; i < N ; i++ )
        {
            result[i] = (raw->alpha)[i];
        }
    }

    LOCK_FN_END;

    return result;
}

iVECTOR SVM_regress::get_tau(void)
{
    LOCK_FN_START;

    iVECTOR result('x',get_N());
    long i,N;

    N = get_N();

    if ( N > 0 )
    {
        for ( i = 0 ; i < N ; i++ )
        {
            result[i] = (raw->tau)[i];
        }
    }

    LOCK_FN_END;

    return result;
}




//
// IO stuff
//

void SVM_regress::prime_io(std::ostream *_where_to, int _echo_level)
{
    where_to   = _where_to;
    echo_level = _echo_level;

    return;
}

void *background_opt_callback_regress(void *caller)
{
    ((SVM_regress *) caller)->optimise_unsafe();

    ((SVM_regress *) caller)->loop_count = 0;

    return caller;
}
