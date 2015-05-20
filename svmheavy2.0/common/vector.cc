
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
// Vector class.
//
// Written by: Alistair Shilton
//             Melbourne University
//




#include <iostream>
#include <string.h>

#include "vector.h"
#include "svdefs.h"
#include "search.h"
#include "c_double.h"
#include "friends.h"
#include "outfilt.h"


const fVECTOR &fixnum(const fVECTOR &whatever);
const fVECTOR &fixnum(const fVECTOR &whatever)
{
    return whatever;
}

//
// Macros
//

//
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=                                                                +=+=
// +=+=                      OPERATOR OVERLOADS                        +=+=
// +=+=                                                                +=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
//

#define UNARY_VECTOR_OPERATION_CONST_RES(what_type,loc_zero,what_left,what_right) \
{                                                                       \
    long i;                                                             \
    long size;                                                          \
    long _lstart,_lend;                                                 \
                                                                        \
    switch ( left_op.v_type )                                           \
    {                                                                   \
        case NORMAL_VECTOR:                                             \
        {                                                               \
            result.make_normal(DO_FORCE);                               \
                                                                        \
            size = left_op.get_effective_size();                        \
                                                                        \
            if ( size > 0 )                                             \
            {                                                           \
                result.pad_vector(size);                                \
                                                                        \
                _lstart = result.offset_start_effective;                \
                _lend   = result.offset_end_effective;                  \
                                                                        \
                for ( i = _lstart ; i <= _lend ; i++ )                  \
                {                                                       \
                    (result.vect)[i-1] = what_left (left_op.get_offset_element(i-_lstart)) what_right; \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case ZERO_VECTOR:                                               \
        {                                                               \
            result.make_scalar(DO_FORCE);                               \
                                                                        \
            result[-1] = what_left (left_op.get_offset_element(0)) what_right; \
                                                                        \
            if ( result[-1] == loc_zero )                               \
            {                                                           \
                result.make_zero(DO_FORCE);                             \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SCALAR_VECTOR:                                             \
        {                                                               \
            result.make_scalar(DO_FORCE);                               \
                                                                        \
            result[-1] = what_left (left_op.get_offset_element(-1)) what_right; \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SELECT_VECTOR:                                             \
        {                                                               \
            result.make_select(DO_FORCE);                               \
            result.set_sel_elm(left_op.get_sel_elm());                  \
                                                                        \
            what_type tempa;                                            \
            what_type tempb;                                            \
                                                                        \
            tempa = what_left (left_op.get_offset_element(-1)) what_right; \
            tempb = what_left (left_op.get_offset_element(42+(left_op.get_sel_elm()))) what_right; \
                                                                        \
            if ( tempb == loc_zero )                                    \
            {                                                           \
                result[-1] = tempa;                                     \
            }                                                           \
                                                                        \
            else                                                        \
            {                                                           \
                THROW_ASSERTB(1,tempa == tempb);                        \
                                                                        \
                result.make_scalar(DO_FORCE);                           \
                                                                        \
                result[-1] = tempa;                                     \
                                                                        \
                if ( result[-1] == loc_zero )                           \
                {                                                       \
                    result.make_zero(DO_FORCE);                         \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        default:                                                        \
        {                                                               \
            L_THROW(2);                                                 \
                                                                        \
            break;                                                      \
        }                                                               \
    }                                                                   \
}

#define UNARY_VECTOR_OPERATION_VARIA(what_type,loc_zero,what_left,what_right) \
{                                                                       \
    long i;                                                             \
    long size;                                                          \
    long _lstart,_lend;                                                 \
                                                                        \
    switch ( left_op.v_type )                                           \
    {                                                                   \
        case NORMAL_VECTOR:                                             \
        {                                                               \
            size = left_op.get_effective_size();                        \
                                                                        \
            if ( size > 0 )                                             \
            {                                                           \
                _lstart = left_op.offset_start_effective;               \
                _lend   = left_op.offset_end_effective;                 \
                                                                        \
                for ( i = _lstart ; i <= _lend ; i++ )                  \
                {                                                       \
                    what_left (left_op.vect)[i-1] what_right;           \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case ZERO_VECTOR:                                               \
        {                                                               \
            left_op.make_scalar(DO_FORCE);                              \
                                                                        \
            left_op[-1] = loc_zero;                                     \
                                                                        \
            what_left left_op[-1] what_right;                           \
                                                                        \
            if ( left_op[-1] == loc_zero )                              \
            {                                                           \
                left_op.make_zero(DO_FORCE);                            \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SCALAR_VECTOR:                                             \
        {                                                               \
            what_left left_op[-1] what_right;                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SELECT_VECTOR:                                             \
        {                                                               \
            what_type tempa;                                            \
            what_type tempb;                                            \
                                                                        \
            tempa = left_op[-1];                                        \
            tempb = loc_zero;                                           \
                                                                        \
            what_left tempa what_right;                                 \
            what_left tempb what_right;                                 \
                                                                        \
            if ( tempb == loc_zero )                                    \
            {                                                           \
                left_op[-1] = tempa;                                    \
            }                                                           \
                                                                        \
            else                                                        \
            {                                                           \
                THROW_ASSERTB(2,tempa == tempb);                        \
                                                                        \
                left_op.make_scalar(DO_FORCE);                          \
                                                                        \
                left_op[-1] = tempa;                                    \
                                                                        \
                if ( left_op[-1] == loc_zero )                          \
                {                                                       \
                    left_op.make_zero(DO_FORCE);                        \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        default:                                                        \
        {                                                               \
            L_THROW(2);                                                 \
                                                                        \
            break;                                                      \
        }                                                               \
    }                                                                   \
}

#define UNARY_VECTOR_OPERATION_VARIA_RES(what_type,loc_zero,what_left,what_right) \
{                                                                       \
    long i;                                                             \
    long size;                                                          \
    long _lstart,_lend;                                                 \
                                                                        \
    result = left_op;                                                   \
                                                                        \
    switch ( left_op.v_type )                                           \
    {                                                                   \
        case NORMAL_VECTOR:                                             \
        {                                                               \
            size = left_op.get_effective_size();                        \
                                                                        \
            if ( size > 0 )                                             \
            {                                                           \
                _lstart = left_op.offset_start_effective;               \
                _lend   = left_op.offset_end_effective;                 \
                                                                        \
                for ( i = _lstart ; i <= _lend ; i++ )                  \
                {                                                       \
                    (result.vect)[i-_lstart] = what_left (left_op.vect)[i-1] what_right; \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case ZERO_VECTOR:                                               \
        {                                                               \
            left_op.make_scalar(DO_FORCE);                              \
            result.make_scalar(DO_FORCE);                               \
                                                                        \
            left_op[-1] = loc_zero;                                     \
            result[-1] = loc_zero;                                      \
                                                                        \
            result[-1] = what_left left_op[-1] what_right;              \
                                                                        \
            if ( left_op[-1] == loc_zero )                              \
            {                                                           \
                left_op.make_zero(DO_FORCE);                            \
            }                                                           \
                                                                        \
            if ( result[-1] == loc_zero )                               \
            {                                                           \
                result.make_zero(DO_FORCE);                             \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SCALAR_VECTOR:                                             \
        {                                                               \
            result[-1] = what_left left_op[-1] what_right;              \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SELECT_VECTOR:                                             \
        {                                                               \
            what_type tempa;                                            \
            what_type tempb;                                            \
            what_type tempc;                                            \
            what_type tempd;                                            \
                                                                        \
            tempa = left_op[-1];                                        \
            tempb = loc_zero;                                           \
                                                                        \
            tempc = what_left tempa what_right;                         \
            tempd = what_left tempb what_right;                         \
                                                                        \
            if ( tempb == loc_zero )                                    \
            {                                                           \
                left_op[-1] = tempa;                                    \
            }                                                           \
                                                                        \
            else                                                        \
            {                                                           \
                THROW_ASSERTB(3,tempa == tempb);                        \
                                                                        \
                left_op.make_scalar(DO_FORCE);                          \
                                                                        \
                left_op[-1] = tempa;                                    \
                                                                        \
                if ( left_op[-1] == loc_zero )                          \
                {                                                       \
                    left_op.make_zero(DO_FORCE);                        \
                }                                                       \
            }                                                           \
                                                                        \
            if ( tempd == loc_zero )                                    \
            {                                                           \
                result[-1] = tempc;                                     \
            }                                                           \
                                                                        \
            else                                                        \
            {                                                           \
                THROW_ASSERTB(80,tempc == tempd);                       \
                                                                        \
                result.make_scalar(DO_FORCE);                           \
                                                                        \
                result[-1] = tempc;                                     \
                                                                        \
                if ( result[-1] == loc_zero )                           \
                {                                                       \
                    result.make_zero(DO_FORCE);                         \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        default:                                                        \
        {                                                               \
            L_THROW(3);                                                 \
                                                                        \
            break;                                                      \
        }                                                               \
    }                                                                   \
}

#define BINARY_VECTOR_OPERATION_CONST_CONST_RES(what_type,loc_zero,what_mid)      \
{                                                                       \
    long i;                                                             \
    long le_size;                                                       \
    long re_size;                                                       \
                                                                        \
    switch ( left_op.v_type )                                           \
    {                                                                   \
        case NORMAL_VECTOR:                                             \
        {                                                               \
            result = left_op;                                           \
                                                                        \
            switch ( right_op.v_type )                                  \
            {                                                           \
                case NORMAL_VECTOR:                                     \
                case SCALAR_VECTOR:                                     \
                {                                                       \
                    le_size = left_op.get_effective_size();             \
                                                                        \
                    if ( le_size > 0 )                                  \
                    {                                                   \
                        for ( i = 1 ; i <= le_size ; i++ )              \
                        {                                               \
                            (result.vect)[i-1] = (left_op.get_offset_element(i-1)) what_mid (right_op.get_offset_element(i-1)); \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case ZERO_VECTOR:                                       \
                {                                                       \
                    break;                                              \
                }                                                       \
                                                                        \
                case SELECT_VECTOR:                                     \
                {                                                       \
                    result[(right_op.get_sel_elm())-1] = (left_op.get_offset_element((right_op.get_sel_elm())-1)) what_mid (right_op.get_offset_element((right_op.get_sel_elm())-1)); \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(1);                                         \
                                                                        \
                    break;                                              \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case ZERO_VECTOR:                                               \
        {                                                               \
            result = right_op;                                          \
                                                                        \
            switch ( right_op.v_type )                                  \
            {                                                           \
                case NORMAL_VECTOR:                                     \
                {                                                       \
                    re_size = right_op.get_effective_size();            \
                                                                        \
                    if ( re_size > 0 )                                  \
                    {                                                   \
                        for ( i = 1 ; i <= re_size ; i++ )              \
                        {                                               \
                            (result.vect)[i-1] = loc_zero what_mid (right_op.get_offset_element(i-1)); \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case ZERO_VECTOR:                                       \
                {                                                       \
                    break;                                              \
                }                                                       \
                                                                        \
                case SCALAR_VECTOR:                                     \
                case SELECT_VECTOR:                                     \
                {                                                       \
                    result[-1] = loc_zero what_mid (right_op.get_offset_element(-1)); \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(2);                                         \
                                                                        \
                    break;                                              \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SCALAR_VECTOR:                                             \
        {                                                               \
            switch ( right_op.v_type )                                  \
            {                                                           \
                case NORMAL_VECTOR:                                     \
                {                                                       \
                    result = right_op;                                  \
                                                                        \
                    re_size = right_op.get_effective_size();            \
                                                                        \
                    if ( re_size > 0 )                                  \
                    {                                                   \
                        for ( i = 1 ; i <= re_size ; i++ )              \
                        {                                               \
                            (result.vect)[i-1] = (left_op.get_offset_element(i-1)) what_mid (right_op.get_offset_element(i-1)); \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case ZERO_VECTOR:                                       \
                {                                                       \
                    result = left_op;                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case SCALAR_VECTOR:                                     \
                {                                                       \
                    result = left_op;                                   \
                                                                        \
                    result[-1] = (left_op.get_offset_element(-1)) what_mid (right_op.get_offset_element(-1)); \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case SELECT_VECTOR:                                     \
                {                                                       \
                    what_type tempa;                                    \
                    what_type tempb;                                    \
                                                                        \
                    tempa = (left_op.get_offset_element(-1)) what_mid (right_op.get_offset_element((right_op.get_sel_elm())-1)); \
                    tempb = (left_op.get_offset_element(-1)) what_mid loc_zero; \
                                                                        \
                    if ( tempb == loc_zero )                            \
                    {                                                   \
                        result = right_op;                              \
                                                                        \
                        result[-1] = tempa;                             \
                    }                                                   \
                                                                        \
                    else                                                \
                    {                                                   \
                        result = left_op;                               \
                                                                        \
                        THROW_ASSERTB(4,tempa == tempb);                \
                                                                        \
                        result[-1] = tempa;                             \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    result = right_op;                                  \
                                                                        \
                    L_THROW(4);                                         \
                                                                        \
                    break;                                              \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SELECT_VECTOR:                                             \
        {                                                               \
            switch ( right_op.v_type )                                  \
            {                                                           \
                case NORMAL_VECTOR:                                     \
                {                                                       \
                    result = right_op;                                  \
                                                                        \
                    re_size = right_op.get_effective_size();            \
                                                                        \
                    if ( re_size > 0 )                                  \
                    {                                                   \
                        for ( i = 1 ; i <= re_size ; i++ )              \
                        {                                               \
                            (result.vect)[i-1] = (left_op.get_offset_element(i-1)) what_mid (right_op.get_offset_element(i-1)); \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case ZERO_VECTOR:                                       \
                {                                                       \
                    result = left_op;                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case SCALAR_VECTOR:                                     \
                {                                                       \
                    what_type tempa;                                    \
                    what_type tempb;                                    \
                                                                        \
                    tempa = (left_op.get_offset_element((left_op.get_sel_elm())-1)) what_mid (right_op.get_offset_element(-1)); \
                    tempb = loc_zero what_mid (right_op.get_offset_element(-1)); \
                                                                        \
                    if ( tempb == loc_zero )                            \
                    {                                                   \
                        result = left_op;                               \
                                                                        \
                        result[-1] = tempa;                             \
                    }                                                   \
                                                                        \
                    else                                                \
                    {                                                   \
                        result = right_op;                              \
                                                                        \
                        THROW_ASSERTB(5,tempa == tempb);                \
                                                                        \
                        result[-1] = tempa;                             \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case SELECT_VECTOR:                                     \
                {                                                       \
                    if ( left_op.get_sel_elm() !=                       \
                         right_op.get_sel_elm()    )                    \
                    {                                                   \
                        what_type tempb;                                \
                        what_type tempc;                                \
                        what_type tempd;                                \
                                                                        \
                        tempb = loc_zero what_mid (right_op.get_offset_element(right_op.get_sel_elm()-1)); \
                        tempc = (left_op.get_offset_element((left_op.get_sel_elm())-1)) what_mid loc_zero; \
                        tempd = loc_zero what_mid loc_zero;             \
                                                                        \
                        if ( tempb == loc_zero )                        \
                        {                                               \
                            result = left_op;                           \
                                                                        \
                            if ( tempc == loc_zero )                    \
                            {                                           \
                                THROW_ASSERTB(6,tempd == loc_zero);     \
                                                                        \
                                result.make_zero(DO_FORCE);             \
                            }                                           \
                                                                        \
                            else                                        \
                            {                                           \
                                THROW_ASSERTB(7,tempd == loc_zero);     \
                                                                        \
                                result[-1] = tempc;                     \
                            }                                           \
                        }                                               \
                                                                        \
                        else                                            \
                        {                                               \
                            result = right_op;                          \
                                                                        \
                            if ( tempc == loc_zero )                    \
                            {                                           \
                                THROW_ASSERTB(8,tempd == loc_zero);     \
                                                                        \
                                result[-1] = tempb;                     \
                            }                                           \
                                                                        \
                            else                                        \
                            {                                           \
                                THROW_ASSERTB(9,(tempd != loc_zero) && (tempb == tempc) && (tempb == tempd)); \
                                                                        \
                                result.make_scalar(DO_FORCE);           \
                                                                        \
                                result[-1] = tempb;                     \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    else                                                \
                    {                                                   \
                        result = left_op;                               \
                                                                        \
                        what_type tempa;                                \
                        what_type tempb;                                \
                                                                        \
                        tempa = (left_op.get_offset_element((left_op.get_sel_elm())-1)) what_mid (right_op.get_offset_element((right_op.get_sel_elm())-1)); \
                        tempb = loc_zero what_mid loc_zero;             \
                                                                        \
                        if ( tempb == loc_zero )                        \
                        {                                               \
                            result[-1] = tempa;                         \
                        }                                               \
                                                                        \
                        else                                            \
                        {                                               \
                            THROW_ASSERTB(10,tempa == tempb);           \
                                                                        \
                            result.make_scalar(DO_FORCE);               \
                                                                        \
                            result[-1] = tempa;                         \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    result = left_op;                                   \
                                                                        \
                    L_THROW(11);                                        \
                                                                        \
                    break;                                              \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        default:                                                        \
        {                                                               \
            result = left_op;                                           \
                                                                        \
            L_THROW(12);                                                \
                                                                        \
            break;                                                      \
        }                                                               \
    }                                                                   \
}

#define BINARY_VECTOR_OPERATION_VARIA_CONST(what_type,loc_zero,what_mid) \
{                                                                       \
    long i;                                                             \
    long le_size;                                                       \
    long re_size;                                                       \
    long _lstart,_lend;                                                 \
                                                                        \
    switch ( left_op.v_type )                                           \
    {                                                                   \
        case NORMAL_VECTOR:                                             \
        {                                                               \
            switch ( right_op.v_type )                                  \
            {                                                           \
                case NORMAL_VECTOR:                                     \
                case SCALAR_VECTOR:                                     \
                {                                                       \
                    le_size = left_op.get_effective_size();             \
                                                                        \
                    if ( le_size > 0 )                                  \
                    {                                                   \
                        _lstart = left_op.offset_start_effective;       \
                        _lend   = left_op.offset_end_effective;         \
                                                                        \
                        for ( i = _lstart ; i <= _lend ; i++ )          \
                        {                                               \
                            (left_op.vect)[i-1] what_mid (right_op.get_offset_element(i-_lstart)); \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case SELECT_VECTOR:                                     \
                {                                                       \
                    left_op[(right_op.get_sel_elm())-1] what_mid (right_op.get_offset_element((right_op.get_sel_elm())-1)); \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case ZERO_VECTOR:                                       \
                {                                                       \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(1);                                         \
                                                                        \
                    break;                                              \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case ZERO_VECTOR:                                               \
        {                                                               \
            switch ( right_op.v_type )                                  \
            {                                                           \
                case NORMAL_VECTOR:                                     \
                {                                                       \
                    re_size = right_op.get_effective_size();            \
                                                                        \
                    left_op.make_normal(DO_FORCE);                      \
                    left_op.pad_vector(re_size);                        \
                                                                        \
                    if ( re_size > 0 )                                  \
                    {                                                   \
                        for ( i = 1 ; i <= re_size ; i++ )              \
                        {                                               \
                            left_op[i-1] what_mid (right_op.get_offset_element(i-1)); \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case ZERO_VECTOR:                                       \
                {                                                       \
                    break;                                              \
                }                                                       \
                                                                        \
                case SCALAR_VECTOR:                                     \
                {                                                       \
                    left_op.make_scalar(DO_FORCE);                      \
                                                                        \
                    left_op[-1] what_mid (right_op.get_offset_element(-1)); \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case SELECT_VECTOR:                                     \
                {                                                       \
                    left_op.make_select(DO_FORCE);                      \
                    left_op.set_sel_elm(right_op.get_sel_elm());        \
                                                                        \
                    left_op.vector_m_diag = loc_zero;                   \
                                                                        \
                    left_op[-1] what_mid (right_op.get_offset_element(-1)); \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(3);                                         \
                                                                        \
                    break;                                              \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SCALAR_VECTOR:                                             \
        {                                                               \
            switch ( right_op.v_type )                                  \
            {                                                           \
                case NORMAL_VECTOR:                                     \
                {                                                       \
                    what_type temp;                                     \
                                                                        \
                    temp = left_op[-1];                                 \
                                                                        \
                    re_size = right_op.get_effective_size();            \
                                                                        \
                    left_op.make_normal(DO_FORCE);                      \
                    left_op.pad_vector(re_size);                        \
                                                                        \
                    if ( re_size > 0 )                                  \
                    {                                                   \
                        for ( i = 1 ; i <= re_size ; i++ )              \
                        {                                               \
                            left_op[i-1] = temp;                        \
                                                                        \
                            left_op[i-1] what_mid (right_op.get_offset_element(i-1)); \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case ZERO_VECTOR:                                       \
                {                                                       \
                    break;                                              \
                }                                                       \
                                                                        \
                case SCALAR_VECTOR:                                     \
                {                                                       \
                    left_op[-1] what_mid (right_op.get_offset_element(-1)); \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case SELECT_VECTOR:                                     \
                {                                                       \
                    what_type tempa;                                    \
                    what_type tempb;                                    \
                                                                        \
                    tempa = left_op[-1];                                \
                    tempb = left_op[-1];                                \
                                                                        \
                    tempa what_mid (right_op.get_offset_element((left_op.get_sel_elm())-1)); \
                    tempb what_mid (right_op.get_offset_element((left_op.get_sel_elm())  )); \
                                                                        \
                    if ( tempb == loc_zero )                            \
                    {                                                   \
                        left_op.make_select(DO_FORCE);                  \
                        left_op.set_sel_elm(right_op.get_sel_elm());    \
                                                                        \
                        left_op[-1] = tempa;                            \
                    }                                                   \
                                                                        \
                    else                                                \
                    {                                                   \
                        THROW_ASSERTB(11,tempa == tempb);               \
                                                                        \
                        left_op.make_scalar(DO_FORCE);                  \
                                                                        \
                        left_op[-1] = tempa;                            \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(5);                                         \
                                                                        \
                    break;                                              \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SELECT_VECTOR:                                             \
        {                                                               \
            switch ( right_op.v_type )                                  \
            {                                                           \
                case NORMAL_VECTOR:                                     \
                {                                                       \
                    long where;                                         \
                    what_type temp;                                     \
                                                                        \
                    temp  = left_op[-1];                                \
                    where = left_op.get_sel_elm();                      \
                                                                        \
                    re_size = right_op.get_effective_size();            \
                                                                        \
                    left_op.make_normal(DO_FORCE);                      \
                    left_op.pad_vector(re_size);                        \
                                                                        \
                    if ( re_size > 0 )                                  \
                    {                                                   \
                        for ( i = 1 ; i <= re_size ; i++ )              \
                        {                                               \
                            if ( i == where )                           \
                            {                                           \
                                left_op[i-1] = temp;                    \
                            }                                           \
                                                                        \
                            left_op[i-1] what_mid (right_op.get_offset_element(i-1)); \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case ZERO_VECTOR:                                       \
                {                                                       \
                    break;                                              \
                }                                                       \
                                                                        \
                case SCALAR_VECTOR:                                     \
                {                                                       \
                    what_type tempa;                                    \
                    what_type tempb;                                    \
                                                                        \
                    tempa = left_op[-1];                                \
                    tempb = loc_zero;                                   \
                                                                        \
                    tempa what_mid (right_op.get_offset_element(-1));   \
                    tempb what_mid (right_op.get_offset_element(-1));   \
                                                                        \
                    if ( tempb == loc_zero )                            \
                    {                                                   \
                        left_op[-1] = tempa;                            \
                    }                                                   \
                                                                        \
                    else                                                \
                    {                                                   \
                        THROW_ASSERTB(12,tempa == tempb);               \
                                                                        \
                        left_op.make_scalar(DO_FORCE);                  \
                                                                        \
                        left_op[-1] = tempa;                            \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case SELECT_VECTOR:                                     \
                {                                                       \
                    if ( left_op.get_sel_elm() !=                       \
                         right_op.get_sel_elm()    )                    \
                    {                                                   \
                        what_type tempb;                                \
                        what_type tempc;                                \
                        what_type tempd;                                \
                                                                        \
                        tempb = left_op[left_op.get_sel_elm()];         \
                        tempc = left_op[left_op.get_sel_elm()-1];       \
                        tempd = left_op[left_op.get_sel_elm()];         \
                                                                        \
                        tempb what_mid (right_op.get_offset_element((right_op.get_sel_elm())-1)); \
                        tempc what_mid (right_op.get_offset_element((right_op.get_sel_elm())  )); \
                        tempd what_mid (right_op.get_offset_element((right_op.get_sel_elm())  )); \
                                                                        \
                        if ( tempb == loc_zero )                        \
                        {                                               \
                            if ( tempc == loc_zero )                    \
                            {                                           \
                                THROW_ASSERTB(13,tempd == loc_zero);    \
                                                                        \
                                left_op.make_zero(DO_FORCE);            \
                            }                                           \
                                                                        \
                            else                                        \
                            {                                           \
                                THROW_ASSERTB(14,tempd == loc_zero);    \
                                                                        \
                                left_op[-1] = tempc;                    \
                            }                                           \
                        }                                               \
                                                                        \
                        else                                            \
                        {                                               \
                            if ( tempc == loc_zero )                    \
                            {                                           \
                                THROW_ASSERTB(15,tempd == loc_zero);    \
                                                                        \
                                left_op.set_sel_elm(right_op.get_sel_elm()); \
                                                                        \
                                left_op[-1] = tempb;                    \
                            }                                           \
                                                                        \
                            else                                        \
                            {                                           \
                                THROW_ASSERTB(16,tempd != loc_zero);    \
                                                                        \
                                if ( ( tempb == tempc ) &&              \
                                     ( tempb == tempd )    )            \
                                {                                       \
                                    left_op.make_scalar(DO_FORCE);      \
                                                                        \
                                    left_op[-1] = tempb;                \
                                }                                       \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    else                                                \
                    {                                                   \
                        what_type tempa;                                \
                        what_type tempb;                                \
                                                                        \
                        tempa = left_op[-1];                            \
                        tempb = loc_zero;                               \
                                                                        \
                        tempa what_mid (right_op.get_offset_element((right_op.get_sel_elm())-1)); \
                        tempb what_mid (right_op.get_offset_element((right_op.get_sel_elm())  )); \
                                                                        \
                        if ( tempb == loc_zero )                        \
                        {                                               \
                            left_op[-1] = tempa;                        \
                        }                                               \
                                                                        \
                        else                                            \
                        {                                               \
                            THROW_ASSERTB(17,tempa == tempb);           \
                                                                        \
                            left_op.make_scalar(DO_FORCE);              \
                                                                        \
                            left_op[-1] = tempa;                        \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(13);                                        \
                                                                        \
                    break;                                              \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        default:                                                        \
        {                                                               \
            L_THROW(14);                                                \
                                                                        \
            break;                                                      \
        }                                                               \
    }                                                                   \
}





//
// Generic vector and matrix relational operators
//
// what_a: either || or && - this defines how individual relationships
//         are concatenated to form complete result
// what_b: actual relation between individual elements
// what_c: what result starts as - generally, 0 if what_a is || or 1 if
//         what_a is &&
//

#define VECTOR_VECTOR_RELATION(loc_zero,what_a,what_b,what_c)           \
{                                                                       \
    long i;                                                             \
    long lsize;                                                         \
    long rsize;                                                         \
                                                                        \
    result = what_c;                                                    \
                                                                        \
    switch ( left_op.v_type )                                           \
    {                                                                   \
        case NORMAL_VECTOR:                                             \
        {                                                               \
            lsize = left_op.get_effective_size();                       \
                                                                        \
            if ( ( right_op.get_effective_size() >  0     ) &&          \
                 ( right_op.get_effective_size() != lsize )    )        \
            {                                                           \
                result = 0;                                             \
            }                                                           \
                                                                        \
            else                                                        \
            {                                                           \
                if ( lsize > 0 )                                        \
                {                                                       \
                    for ( i = 1 ; i <=  lsize ; i++ )                   \
                    {                                                   \
                        result =                                        \
                        ( result what_a                                 \
                        ( left_op.get_offset_element(i-1) what_b        \
                          right_op.get_offset_element(i-1) ) );         \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case ZERO_VECTOR:                                               \
        {                                                               \
            switch ( right_op.v_type )                                  \
            {                                                           \
                case NORMAL_VECTOR:                                     \
                {                                                       \
                    rsize = right_op.get_effective_size();              \
                                                                        \
                    if ( rsize > 0 )                                    \
                    {                                                   \
                        for ( i = 1 ; i <=  rsize ; i++ )               \
                        {                                               \
                            result =                                    \
                            ( result what_a                             \
                            ( left_op.get_offset_element(i-1) what_b    \
                              right_op.get_offset_element(i-1) ) );     \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case ZERO_VECTOR:                                       \
                {                                                       \
                    i = 1;                                              \
                                                                        \
                    result =                                            \
                    ( result what_a                                     \
                    ( left_op.get_offset_element(i-1) what_b            \
                      right_op.get_offset_element(i-1) ) );             \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case SCALAR_VECTOR:                                     \
                {                                                       \
                    i = 1;                                              \
                                                                        \
                    result =                                            \
                    ( result what_a                                     \
                    ( left_op.get_offset_element(i-1) what_b            \
                      right_op.get_offset_element(i-1) ) );             \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case SELECT_VECTOR:                                     \
                {                                                       \
                    i = right_op.get_sel_elm();                         \
                                                                        \
                    result =                                            \
                    ( result what_a                                     \
                    ( left_op.get_offset_element(i-1) what_b            \
                      right_op.get_offset_element(i-1) ) );             \
                                                                        \
                    i = right_op.get_sel_elm()+1;                       \
                                                                        \
                    result =                                            \
                    ( result what_a                                     \
                    ( left_op.get_offset_element(i-1) what_b            \
                      right_op.get_offset_element(i-1) ) );             \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(1);                                         \
                                                                        \
                    break;                                              \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SCALAR_VECTOR:                                             \
        {                                                               \
            switch ( right_op.v_type )                                  \
            {                                                           \
                case NORMAL_VECTOR:                                     \
                {                                                       \
                    rsize = right_op.get_effective_size();              \
                                                                        \
                    if ( rsize > 0 )                                    \
                    {                                                   \
                        for ( i = 1 ; i <=  rsize ; i++ )               \
                        {                                               \
                            result =                                    \
                            ( result what_a                             \
                            ( left_op.get_offset_element(i-1) what_b    \
                              right_op.get_offset_element(i-1) ) );     \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case ZERO_VECTOR:                                       \
                {                                                       \
                    i = 1;                                              \
                                                                        \
                    result =                                            \
                    ( result what_a                                     \
                    ( left_op.get_offset_element(i-1) what_b            \
                      right_op.get_offset_element(i-1) ) );             \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case SCALAR_VECTOR:                                     \
                {                                                       \
                    i = 1;                                              \
                                                                        \
                    result =                                            \
                    ( result what_a                                     \
                    ( left_op.get_offset_element(i-1) what_b            \
                      right_op.get_offset_element(i-1) ) );             \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case SELECT_VECTOR:                                     \
                {                                                       \
                    i = right_op.get_sel_elm();                         \
                                                                        \
                    result =                                            \
                    ( result what_a                                     \
                    ( left_op.get_offset_element(i-1) what_b            \
                      right_op.get_offset_element(i-1) ) );             \
                                                                        \
                    i = right_op.get_sel_elm()+1;                       \
                                                                        \
                    result =                                            \
                    ( result what_a                                     \
                    ( left_op.get_offset_element(i-1) what_b            \
                      right_op.get_offset_element(i-1) ) );             \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(2);                                         \
                                                                        \
                    break;                                              \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SELECT_VECTOR:                                             \
        {                                                               \
            switch ( right_op.v_type )                                  \
            {                                                           \
                case NORMAL_VECTOR:                                     \
                {                                                       \
                    rsize = right_op.get_effective_size();              \
                                                                        \
                    if ( rsize > 0 )                                    \
                    {                                                   \
                        for ( i = 1 ; i <=  rsize ; i++ )               \
                        {                                               \
                            result =                                    \
                            ( result what_a                             \
                            ( left_op.get_offset_element(i-1) what_b    \
                              right_op.get_offset_element(i-1) ) );     \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case ZERO_VECTOR:                                       \
                {                                                       \
                    i = left_op.get_sel_elm();                          \
                                                                        \
                    result =                                            \
                    ( result what_a                                     \
                    ( left_op.get_offset_element(i-1) what_b            \
                      right_op.get_offset_element(i-1) ) );             \
                                                                        \
                    i = left_op.get_sel_elm() + 1;                      \
                                                                        \
                    result =                                            \
                    ( result what_a                                     \
                    ( left_op.get_offset_element(i-1) what_b            \
                      right_op.get_offset_element(i-1) ) );             \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case SCALAR_VECTOR:                                     \
                {                                                       \
                    i = left_op.get_sel_elm();                          \
                                                                        \
                    result =                                            \
                    ( result what_a                                     \
                    ( left_op.get_offset_element(i-1) what_b            \
                      right_op.get_offset_element(i-1) ) );             \
                                                                        \
                    i = left_op.get_sel_elm() + 1;                      \
                                                                        \
                    result =                                            \
                    ( result what_a                                     \
                    ( left_op.get_offset_element(i-1) what_b            \
                      right_op.get_offset_element(i-1) ) );             \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case SELECT_VECTOR:                                     \
                {                                                       \
                    i = left_op.get_sel_elm();                          \
                                                                        \
                    result =                                            \
                    ( result what_a                                     \
                    ( left_op.get_offset_element(i-1) what_b            \
                      right_op.get_offset_element(i-1) ) );             \
                                                                        \
                    i = left_op.get_sel_elm() + 1;                      \
                                                                        \
                    result =                                            \
                    ( result what_a                                     \
                    ( left_op.get_offset_element(i-1) what_b            \
                      right_op.get_offset_element(i-1) ) );             \
                                                                        \
                    i = right_op.get_sel_elm();                         \
                                                                        \
                    result =                                            \
                    ( result what_a                                     \
                    ( left_op.get_offset_element(i-1) what_b            \
                      right_op.get_offset_element(i-1) ) );             \
                                                                        \
                    i = right_op.get_sel_elm() + 1;                     \
                                                                        \
                    result =                                            \
                    ( result what_a                                     \
                    ( left_op.get_offset_element(i-1) what_b            \
                      right_op.get_offset_element(i-1) ) );             \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(3);                                         \
                                                                        \
                    break;                                              \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        default:                                                        \
        {                                                               \
            L_THROW(4);                                                 \
                                                                        \
            break;                                                      \
        }                                                               \
    }                                                                   \
}


#define VECTOR_SCALAR_RELATION(what_a,what_b,what_c,loc_zero)           \
{                                                                       \
    long i;                                                             \
    long lsize;                                                         \
                                                                        \
    result = what_c;                                                    \
                                                                        \
    switch ( left_op.v_type )                                           \
    {                                                                   \
        case NORMAL_VECTOR:                                             \
        {                                                               \
            lsize = left_op.get_effective_size();                       \
                                                                        \
            if ( lsize > 0 )                                            \
            {                                                           \
                for ( i = 1 ; i <=  lsize ; i++ )                       \
                {                                                       \
                    result =                                            \
                    ( result what_a                                     \
                    ( left_op.get_offset_element(i-1) what_b right_op ) ); \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case ZERO_VECTOR:                                               \
        {                                                               \
            result = ( result what_a ( loc_zero what_b right_op ) );    \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SCALAR_VECTOR:                                             \
        {                                                               \
            result = ( result what_a ( left_op.get_offset_element(-1) what_b right_op ) ); \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SELECT_VECTOR:                                             \
        {                                                               \
            result = ( result what_a ( loc_zero what_b right_op ) );    \
            result = ( result what_a ( left_op.get_offset_element(-1) what_b right_op ) ); \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        default:                                                        \
        {                                                               \
            L_THROW(4);                                                 \
                                                                        \
            break;                                                      \
        }                                                               \
    }                                                                   \
}



//
// Stream IO macros
//

#define VECTOR_OUT                                                      \
{                                                                       \
    long i;                                                             \
                                                                        \
    if ( source.size == source.get_effective_size() )                   \
    {                                                                   \
        output << "size: " << source.size << "\n";                      \
    }                                                                   \
                                                                        \
    else                                                                \
    {                                                                   \
        output << "modified size (-10*size)-10: ";                      \
        output << -(10+(10*(source.size))) << "\n";                     \
    }                                                                   \
                                                                        \
    switch ( source.v_type )                                            \
    {                                                                   \
        case NORMAL_VECTOR:                                             \
        {                                                               \
            if ( source.size > 0 )                                      \
            {                                                           \
                if ( source.size != source.get_effective_size() )       \
                {                                                       \
                    output << "offset start: " << source.offset_start_effective << "\n";   \
                    output << "offset end:   " << source.offset_end_effective   << "\n\n"; \
                }                                                       \
                                                                        \
                output << "vector \n\n";                                \
                                                                        \
                for ( i = 1 ; i <= source.size ; i++ )                  \
                {                                                       \
                    output << i << ": " << fixnum((source.vect)[i-1]) << "\n";  \
                }                                                       \
            }                                                           \
                                                                        \
            output << "\n";                                             \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case ZERO_VECTOR:                                               \
        {                                                               \
            break;                                                      \
        }                                                               \
                                                                        \
        case SCALAR_VECTOR:                                             \
        {                                                               \
            output << "scalar: " << fixnum(source.get_offset_element(-1)) << "\n\n"; \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SELECT_VECTOR:                                             \
        {                                                               \
            output << "selector: " << source.get_sel_elm()                  << "\n";   \
            output << "scalar:   " << fixnum(source.get_offset_element(-1)) << "\n\n"; \
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
}

#define VECTOR_OUT_INTLONG                                              \
{                                                                       \
    long i;                                                             \
                                                                        \
    if ( source.size == source.get_effective_size() )                   \
    {                                                                   \
        output << "size: " << source.size << "\n";                      \
    }                                                                   \
                                                                        \
    else                                                                \
    {                                                                   \
        output << "modified size (-10*size)-10: ";                      \
        output << -(10+(10*(source.size))) << "\n";                     \
    }                                                                   \
                                                                        \
    switch ( source.v_type )                                            \
    {                                                                   \
        case NORMAL_VECTOR:                                             \
        {                                                               \
            if ( source.size > 0 )                                      \
            {                                                           \
                if ( source.size != source.get_effective_size() )       \
                {                                                       \
                    output << "offset start: " << source.offset_start_effective << "\n";   \
                    output << "offset end:   " << source.offset_end_effective   << "\n\n"; \
                }                                                       \
                                                                        \
                output << "vector \n\n";                                \
                                                                        \
                for ( i = 1 ; i <= source.size ; i++ )                  \
                {                                                       \
                    output << i << ": " << (source.vect)[i-1] << "\n";  \
                }                                                       \
            }                                                           \
                                                                        \
            output << "\n";                                             \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case ZERO_VECTOR:                                               \
        {                                                               \
            break;                                                      \
        }                                                               \
                                                                        \
        case SCALAR_VECTOR:                                             \
        {                                                               \
            output << "scalar: " << source.get_offset_element(-1) << "\n\n"; \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SELECT_VECTOR:                                             \
        {                                                               \
            output << "selector: " << source.get_sel_elm()          << "\n";   \
            output << "scalar:   " << source.get_offset_element(-1) << "\n\n"; \
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
}

#define VECTOR_IN                                                       \
{                                                                       \
    long i;                                                             \
    wait_dummy howzat;                                                  \
    int modified_size = 0;                                              \
                                                                        \
    int size;                                                           \
    int v_type;                                                         \
                                                                        \
    long offset_start_effective;                                        \
    long offset_end_effective;                                          \
                                                                        \
    std::ostream *where_to;                                             \
    int echo_level;                                                     \
                                                                        \
    where_to   = dest.where_to;                                         \
    echo_level = dest.echo_level;                                       \
                                                                        \
    dest.prime_io(NULL,0);                                              \
                                                                        \
    if ( where_to == NULL )                                             \
    {                                                                   \
        dest.make_zero(DO_FORCE);                                       \
                                                                        \
        input >> howzat; input >> size;                                 \
                                                                        \
        switch ( size )                                                 \
        {                                                               \
            case ZERO_VECTOR_SIZE:                                      \
            {                                                           \
                v_type = ZERO_VECTOR;                                   \
                                                                        \
                break;                                                  \
            }                                                           \
                                                                        \
            case SCALAR_VECTOR_SIZE:                                    \
            {                                                           \
                v_type = SCALAR_VECTOR;                                 \
                                                                        \
                break;                                                  \
            }                                                           \
                                                                        \
            case SELECT_VECTOR_SIZE:                                    \
            {                                                           \
                v_type = SELECT_VECTOR;                                 \
                                                                        \
                break;                                                  \
            }                                                           \
                                                                        \
            default:                                                    \
            {                                                           \
                if ( size <= -10 )                                      \
                {                                                       \
                    modified_size = 1;                                  \
                                                                        \
                    size += 10;                                         \
                    size /= -10;                                        \
                }                                                       \
                                                                        \
                THROW_ASSERTB(18,size >= 0);                            \
                                                                        \
                v_type = NORMAL_VECTOR;                                 \
                                                                        \
                break;                                                  \
            }                                                           \
        }                                                               \
                                                                        \
        switch ( v_type )                                               \
        {                                                               \
            case NORMAL_VECTOR:                                         \
            {                                                           \
                dest.make_normal(DO_FORCE);                             \
                                                                        \
                if ( size > 0 )                                         \
                {                                                       \
                    dest.pad_vector(size);                              \
                                                                        \
                    if ( modified_size )                                \
                    {                                                   \
                        input >> howzat; input >> offset_start_effective; \
                        input >> howzat; input >> offset_end_effective; \
                    }                                                   \
                                                                        \
                    else                                                \
                    {                                                   \
                        offset_start_effective = 1;                     \
                        offset_end_effective   = size;                  \
                    }                                                   \
                                                                        \
                    for ( i = 1 ; i <= dest.size ; i++ )                \
                    {                                                   \
                        input >> howzat; input >> dest[i-1];            \
                    }                                                   \
                                                                        \
                    dest.set_offset_start(offset_start_effective);      \
                    dest.set_offset_end(offset_end_effective);          \
                }                                                       \
                                                                        \
                break;                                                  \
            }                                                           \
                                                                        \
            case ZERO_VECTOR:                                           \
            {                                                           \
                dest.make_zero(DO_FORCE);                               \
                                                                        \
                break;                                                  \
            }                                                           \
                                                                        \
            case SCALAR_VECTOR:                                         \
            {                                                           \
                dest.make_scalar(DO_FORCE);                             \
                                                                        \
                input >> howzat; input >> dest.vector_m_diag;           \
                                                                        \
                break;                                                  \
            }                                                           \
                                                                        \
            case SELECT_VECTOR:                                         \
            {                                                           \
                dest.make_select(DO_FORCE);                             \
                                                                        \
                input >> howzat; input >> dest.sel_elm;                 \
                input >> howzat; input >> dest.vector_m_diag;           \
                                                                        \
                break;                                                  \
            }                                                           \
                                                                        \
            default:                                                    \
            {                                                           \
                L_THROW(0);                                             \
                                                                        \
                break;                                                  \
            }                                                           \
        }                                                               \
    }                                                                   \
                                                                        \
    else                                                                \
    {                                                                   \
        dest.make_zero(DO_FORCE);                                       \
                                                                        \
        *where_to << "Vector sizes: >= 0 gives normal vector.\n";       \
        *where_to << "              -1 gives zero vector.\n";           \
        *where_to << "              -2 gives scalar vector.\n";         \
        *where_to << "              -3 gives select vector.\n\n";       \
                                                                        \
        *where_to << "Vector size: ";                                   \
        input >> size;                                                  \
        if ( echo_level ) { *where_to << size << "\n"; }                \
                                                                        \
        switch ( size )                                                 \
        {                                                               \
            case ZERO_VECTOR_SIZE:                                      \
            {                                                           \
                v_type = ZERO_VECTOR;                                   \
                                                                        \
                break;                                                  \
            }                                                           \
                                                                        \
            case SCALAR_VECTOR_SIZE:                                    \
            {                                                           \
                v_type = SCALAR_VECTOR;                                 \
                                                                        \
                break;                                                  \
            }                                                           \
                                                                        \
            case SELECT_VECTOR_SIZE:                                    \
            {                                                           \
                v_type = SELECT_VECTOR;                                 \
                                                                        \
                break;                                                  \
            }                                                           \
                                                                        \
            default:                                                    \
            {                                                           \
                THROW_ASSERTB(19,size >= 0);                            \
                                                                        \
                v_type = NORMAL_VECTOR;                                 \
                                                                        \
                break;                                                  \
            }                                                           \
        }                                                               \
                                                                        \
        switch ( v_type )                                               \
        {                                                               \
            case NORMAL_VECTOR:                                         \
            {                                                           \
                dest.make_normal(DO_FORCE);                             \
                                                                        \
                if ( size > 0 )                                         \
                {                                                       \
                    dest.pad_vector(size);                              \
                                                                        \
                    offset_start_effective = 1;                         \
                    offset_end_effective   = size;                      \
                                                                        \
                    *where_to << "Enter vector:\n";                     \
                                                                        \
                    for ( i = 1 ; i <= dest.size ; i++ )                \
                    {                                                   \
                        input >> dest[i-1];                             \
                        if ( echo_level ) { *where_to << dest[i-1] << "\n"; } \
                    }                                                   \
                                                                        \
                    dest.set_offset_start(offset_start_effective);      \
                    dest.set_offset_end(offset_end_effective);          \
                }                                                       \
                                                                        \
                break;                                                  \
            }                                                           \
                                                                        \
            case ZERO_VECTOR:                                           \
            {                                                           \
                dest.make_zero(DO_FORCE);                               \
                                                                        \
                break;                                                  \
            }                                                           \
                                                                        \
            case SCALAR_VECTOR:                                         \
            {                                                           \
                dest.make_scalar(DO_FORCE);                             \
                                                                        \
                *where_to << "Scalar: ";                                \
                input >> dest.vector_m_diag;                            \
                if ( echo_level ) { *where_to << dest.vector_m_diag << "\n"; } \
                                                                        \
                break;                                                  \
            }                                                           \
                                                                        \
            case SELECT_VECTOR:                                         \
            {                                                           \
                dest.make_select(DO_FORCE);                             \
                                                                        \
                *where_to << "Which element is non-zero: ";             \
                input >> dest.sel_elm;                                  \
                if ( echo_level ) { *where_to << dest.sel_elm << "\n"; }\
                *where_to << "Scalar: ";                                \
                input >> dest.vector_m_diag;                            \
                if ( echo_level ) { *where_to << dest.vector_m_diag << "\n"; } \
                                                                        \
                break;                                                  \
            }                                                           \
                                                                        \
            default:                                                    \
            {                                                           \
                L_THROW(0);                                             \
                                                                        \
                break;                                                  \
            }                                                           \
        }                                                               \
    }                                                                   \
}

#define VECTOR_IN_PROP_PRIME                                            \
{                                                                       \
    long i;                                                             \
    wait_dummy howzat;                                                  \
    int modified_size = 0;                                              \
                                                                        \
    int size;                                                           \
    int v_type;                                                         \
                                                                        \
    long offset_start_effective;                                        \
    long offset_end_effective;                                          \
                                                                        \
    std::ostream *where_to;                                             \
    int echo_level;                                                     \
                                                                        \
    where_to   = dest.where_to;                                         \
    echo_level = dest.echo_level;                                       \
                                                                        \
    dest.prime_io(NULL,0);                                              \
                                                                        \
    if ( where_to == NULL )                                             \
    {                                                                   \
        dest.make_zero(DO_FORCE);                                       \
                                                                        \
        input >> howzat; input >> size;                                 \
                                                                        \
        switch ( size )                                                 \
        {                                                               \
            case ZERO_VECTOR_SIZE:                                      \
            {                                                           \
                v_type = ZERO_VECTOR;                                   \
                                                                        \
                break;                                                  \
            }                                                           \
                                                                        \
            case SCALAR_VECTOR_SIZE:                                    \
            {                                                           \
                v_type = SCALAR_VECTOR;                                 \
                                                                        \
                break;                                                  \
            }                                                           \
                                                                        \
            case SELECT_VECTOR_SIZE:                                    \
            {                                                           \
                v_type = SELECT_VECTOR;                                 \
                                                                        \
                break;                                                  \
            }                                                           \
                                                                        \
            default:                                                    \
            {                                                           \
                if ( size <= -10 )                                      \
                {                                                       \
                    modified_size = 1;                                  \
                                                                        \
                    size += 10;                                         \
                    size /= -10;                                        \
                }                                                       \
                                                                        \
                THROW_ASSERTB(20,size >= 0);                            \
                                                                        \
                v_type = NORMAL_VECTOR;                                 \
                                                                        \
                break;                                                  \
            }                                                           \
        }                                                               \
                                                                        \
        switch ( v_type )                                               \
        {                                                               \
            case NORMAL_VECTOR:                                         \
            {                                                           \
                dest.make_normal(DO_FORCE);                             \
                                                                        \
                if ( size > 0 )                                         \
                {                                                       \
                    dest.pad_vector(size);                              \
                                                                        \
                    if ( modified_size )                                \
                    {                                                   \
                        input >> howzat; input >> offset_start_effective; \
                        input >> howzat; input >> offset_end_effective; \
                    }                                                   \
                                                                        \
                    else                                                \
                    {                                                   \
                        offset_start_effective = 1;                     \
                        offset_end_effective   = size;                  \
                    }                                                   \
                                                                        \
                    for ( i = 1 ; i <= dest.size ; i++ )                \
                    {                                                   \
                        input >> howzat; input >> dest[i-1];            \
                    }                                                   \
                                                                        \
                    dest.set_offset_start(offset_start_effective);      \
                    dest.set_offset_end(offset_end_effective);          \
                }                                                       \
                                                                        \
                break;                                                  \
            }                                                           \
                                                                        \
            case ZERO_VECTOR:                                           \
            {                                                           \
                dest.make_zero(DO_FORCE);                               \
                                                                        \
                break;                                                  \
            }                                                           \
                                                                        \
            case SCALAR_VECTOR:                                         \
            {                                                           \
                dest.make_scalar(DO_FORCE);                             \
                                                                        \
                input >> howzat; input >> dest.vector_m_diag;           \
                                                                        \
                break;                                                  \
            }                                                           \
                                                                        \
            case SELECT_VECTOR:                                         \
            {                                                           \
                dest.make_select(DO_FORCE);                             \
                                                                        \
                input >> howzat; input >> dest.sel_elm;                 \
                input >> howzat; input >> dest.vector_m_diag;           \
                                                                        \
                break;                                                  \
            }                                                           \
                                                                        \
            default:                                                    \
            {                                                           \
                L_THROW(0);                                             \
                                                                        \
                break;                                                  \
            }                                                           \
        }                                                               \
    }                                                                   \
                                                                        \
    else                                                                \
    {                                                                   \
        dest.make_zero(DO_FORCE);                                       \
                                                                        \
        *where_to << "Vector sizes: >= 0 gives normal vector.\n";       \
        *where_to << "              -1 gives zero vector.\n";           \
        *where_to << "              -2 gives scalar vector.\n";         \
        *where_to << "              -3 gives select vector.\n\n";       \
                                                                        \
        *where_to << "Vector size: ";                                   \
        input >> size;                                                  \
        if ( echo_level ) { *where_to << size << "\n"; }                \
                                                                        \
        switch ( size )                                                 \
        {                                                               \
            case ZERO_VECTOR_SIZE:                                      \
            {                                                           \
                v_type = ZERO_VECTOR;                                   \
                                                                        \
                break;                                                  \
            }                                                           \
                                                                        \
            case SCALAR_VECTOR_SIZE:                                    \
            {                                                           \
                v_type = SCALAR_VECTOR;                                 \
                                                                        \
                break;                                                  \
            }                                                           \
                                                                        \
            case SELECT_VECTOR_SIZE:                                    \
            {                                                           \
                v_type = SELECT_VECTOR;                                 \
                                                                        \
                break;                                                  \
            }                                                           \
                                                                        \
            default:                                                    \
            {                                                           \
                THROW_ASSERTB(21,size >= 0);                            \
                                                                        \
                v_type = NORMAL_VECTOR;                                 \
                                                                        \
                break;                                                  \
            }                                                           \
        }                                                               \
                                                                        \
        switch ( v_type )                                               \
        {                                                               \
            case NORMAL_VECTOR:                                         \
            {                                                           \
                dest.make_normal(DO_FORCE);                             \
                                                                        \
                if ( size > 0 )                                         \
                {                                                       \
                    dest.pad_vector(size);                              \
                                                                        \
                    offset_start_effective = 1;                         \
                    offset_end_effective   = size;                      \
                                                                        \
                    *where_to << "Enter vector:\n";                     \
                                                                        \
                    for ( i = 1 ; i <= dest.size ; i++ )                \
                    {                                                   \
                        (dest[i-1]).prime_io(where_to,echo_level);      \
                        input >> dest[i-1];                             \
                        if ( echo_level ) { *where_to << dest[i-1] << "\n"; } \
                    }                                                   \
                                                                        \
                    dest.set_offset_start(offset_start_effective);      \
                    dest.set_offset_end(offset_end_effective);          \
                }                                                       \
                                                                        \
                break;                                                  \
            }                                                           \
                                                                        \
            case ZERO_VECTOR:                                           \
            {                                                           \
                dest.make_zero(DO_FORCE);                               \
                                                                        \
                break;                                                  \
            }                                                           \
                                                                        \
            case SCALAR_VECTOR:                                         \
            {                                                           \
                dest.make_scalar(DO_FORCE);                             \
                                                                        \
                *where_to << "Scalar: ";                                \
                (dest.vector_m_diag).prime_io(where_to,echo_level);     \
                input >> dest.vector_m_diag;                            \
                if ( echo_level ) { *where_to << dest.vector_m_diag << "\n"; } \
                                                                        \
                break;                                                  \
            }                                                           \
                                                                        \
            case SELECT_VECTOR:                                         \
            {                                                           \
                dest.make_select(DO_FORCE);                             \
                                                                        \
                *where_to << "Which element is non-zero: ";             \
                input >> dest.sel_elm;                                  \
                if ( echo_level ) { *where_to << dest.sel_elm << "\n"; } \
                *where_to << "Scalar: ";                                \
                (dest.vector_m_diag).prime_io(where_to,echo_level);     \
                input >> dest.vector_m_diag;                            \
                if ( echo_level ) { *where_to << dest.vector_m_diag << "\n"; } \
                                                                        \
                break;                                                  \
            }                                                           \
                                                                        \
            default:                                                    \
            {                                                           \
                L_THROW(0);                                             \
                                                                        \
                break;                                                  \
            }                                                           \
        }                                                               \
    }                                                                   \
}







#define VECTOR_MULT_VECTOR                                              \
{                                                                       \
    long i;                                                             \
    long le_size;                                                       \
    long re_size;                                                       \
                                                                        \
    switch ( left_op.v_type )                                           \
    {                                                                   \
        case NORMAL_VECTOR:                                             \
        {                                                               \
            switch ( right_op.v_type )                                  \
            {                                                           \
                case NORMAL_VECTOR:                                     \
                {                                                       \
                    le_size = left_op.get_effective_size();             \
                    re_size = right_op.get_effective_size();            \
                                                                        \
                    THROW_ASSERTB(24,le_size == re_size);               \
                                                                        \
                    if ( le_size > 0 )                                  \
                    {                                                   \
                        for ( i = 1 ; i <= le_size ; i++ )              \
                        {                                               \
                            result += ((left_op.get_offset_element(i-1))*(right_op.get_offset_element(i-1))); \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case ZERO_VECTOR:                                       \
                {                                                       \
                    break;                                              \
                }                                                       \
                                                                        \
                case SCALAR_VECTOR:                                     \
                {                                                       \
                    le_size = left_op.get_effective_size();             \
                                                                        \
                    if ( le_size > 0 )                                  \
                    {                                                   \
                        for ( i = 1 ; i <= le_size ; i++ )              \
                        {                                               \
                            result += ((left_op.get_offset_element(i-1))*(right_op.get_offset_element(i-1))); \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case SELECT_VECTOR:                                     \
                {                                                       \
                    i = right_op.get_sel_elm();                         \
                                                                        \
                    result += ((left_op.get_offset_element(i-1))*(right_op.get_offset_element(i-1))); \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(5);                                         \
                                                                        \
                    break;                                              \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case ZERO_VECTOR:                                               \
        {                                                               \
            break;                                                      \
        }                                                               \
                                                                        \
        case SCALAR_VECTOR:                                             \
        {                                                               \
            switch ( right_op.v_type )                                  \
            {                                                           \
                case NORMAL_VECTOR:                                     \
                {                                                       \
                    re_size = right_op.get_effective_size();            \
                                                                        \
                    if ( re_size > 0 )                                  \
                    {                                                   \
                        for ( i = 1 ; i <= re_size ; i++ )              \
                        {                                               \
                            result += ((left_op.get_offset_element(i-1))*(right_op.get_offset_element(i-1))); \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case ZERO_VECTOR:                                       \
                {                                                       \
                    break;                                              \
                }                                                       \
                                                                        \
                case SCALAR_VECTOR:                                     \
                {                                                       \
                    L_THROW(55);                                        \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case SELECT_VECTOR:                                     \
                {                                                       \
                    i = right_op.get_sel_elm();                         \
                                                                        \
                    result += ((left_op.get_offset_element(i-1))*(right_op.get_offset_element(i-1))); \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(7);                                         \
                                                                        \
                    break;                                              \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SELECT_VECTOR:                                             \
        {                                                               \
            switch ( right_op.v_type )                                  \
            {                                                           \
                case NORMAL_VECTOR:                                     \
                case SCALAR_VECTOR:                                     \
                {                                                       \
                    i = left_op.get_sel_elm();                          \
                                                                        \
                    result += ((left_op.get_offset_element(i-1))*(right_op.get_offset_element(i-1))); \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case ZERO_VECTOR:                                       \
                {                                                       \
                    break;                                              \
                }                                                       \
                                                                        \
                case SELECT_VECTOR:                                     \
                {                                                       \
                    if ( left_op.get_sel_elm() == right_op.get_sel_elm() ) \
                    {                                                   \
                        i = left_op.get_sel_elm();                      \
                                                                        \
                        result += ((left_op.get_offset_element(i-1))*(right_op.get_offset_element(i-1))); \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(8);                                         \
                                                                        \
                    break;                                              \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        default:                                                        \
        {                                                               \
            L_THROW(9);                                                 \
                                                                        \
            break;                                                      \
        }                                                               \
    }                                                                   \
}



#define VECTOR_VOID_CONSTRUCTOR(loc_zero)                               \
{                                                                       \
    size = ZERO_VECTOR_SIZE;                                            \
    vector_m_diag = loc_zero;                                           \
    vect = NULL;                                                        \
    v_type = ZERO_VECTOR;                                               \
    where_to = NULL;                                                    \
    echo_level = 0;                                                     \
}

#define GENERIC_VECTOR_CONSTRUCTOR(what_type,loc_zero,loc_one)          \
{                                                                       \
    VECTOR_VOID_CONSTRUCTOR(loc_zero);                                  \
                                                                        \
    long i;                                                             \
                                                                        \
    where_to = NULL;                                                    \
    echo_level = 0;                                                     \
                                                                        \
    size = _size;                                                       \
                                                                        \
    switch ( size )                                                     \
    {                                                                   \
        case ZERO_VECTOR_SIZE:   v_type = ZERO_VECTOR;   break;         \
        case SCALAR_VECTOR_SIZE: v_type = SCALAR_VECTOR; break;         \
        case SELECT_VECTOR_SIZE: v_type = SELECT_VECTOR; break;         \
                                                                        \
        default:                                                        \
        {                                                               \
            THROW_ASSERTB(25,size >= 0);                                \
                                                                        \
            v_type = NORMAL_VECTOR;                                     \
                                                                        \
            break;                                                      \
        }                                                               \
    }                                                                   \
                                                                        \
    vector_m_diag = loc_zero;                                           \
                                                                        \
    switch ( v_type )                                                   \
    {                                                                   \
        case ZERO_VECTOR:                                               \
        {                                                               \
            break;                                                      \
        }                                                               \
                                                                        \
        case SCALAR_VECTOR:                                             \
        {                                                               \
            vector_m_diag = loc_one;                                    \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SELECT_VECTOR:                                             \
        {                                                               \
            vector_m_diag = loc_one;                                    \
            sel_elm = 1;                                                \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        default:                                                        \
        {                                                               \
            alloc_lookahead = size+10;                                  \
                                                                        \
            REPNEWB(vect,what_type,(size+10));                          \
                                                                        \
            reset_offsets();                                            \
                                                                        \
            if ( size > 0 )                                             \
            {                                                           \
                _size = size;                                           \
                size  = 0;                                              \
                                                                        \
                for ( i = 1 ; i <= _size ; i++ )                        \
                {                                                       \
                    addend();                                           \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
    }                                                                   \
}

#define GENERIC_VECTOR_COPY_CONSTRUCTOR(what_type,loc_zero)             \
{                                                                       \
    VECTOR_VOID_CONSTRUCTOR(loc_zero);                                  \
                                                                        \
    long jjj;                                                           \
                                                                        \
    where_to = NULL;                                                    \
    echo_level = 0;                                                     \
                                                                        \
    vector_m_diag = source.vector_m_diag;                               \
    sel_elm = source.sel_elm;                                           \
                                                                        \
    v_type = source.v_type;                                             \
                                                                        \
    size = source.size;                                                 \
                                                                        \
    switch ( v_type )                                                   \
    {                                                                   \
        case NORMAL_VECTOR:                                             \
        {                                                               \
            alloc_lookahead = (source.size)+10;                         \
                                                                        \
            REPNEWB(vect,what_type,((source.size)+10));                 \
                                                                        \
            size = 0;                                                   \
                                                                        \
            reset_offsets();                                            \
                                                                        \
            if ( source.size > 0 )                                      \
            {                                                           \
                for ( jjj = 1 ; jjj <= source.size ; jjj++ )            \
                {                                                       \
                    addend(((source.vect)[jjj-1]));                     \
                }                                                       \
            }                                                           \
                                                                        \
            offset_start_effective = source.offset_start_effective;     \
            offset_end_effective   = source.offset_end_effective;       \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case ZERO_VECTOR:                                               \
        case SCALAR_VECTOR:                                             \
        case SELECT_VECTOR:                                             \
        {                                                               \
            break;                                                      \
        }                                                               \
                                                                        \
        default:                                                        \
        {                                                               \
            L_THROW(1);                                                 \
                                                                        \
            break;                                                      \
        }                                                               \
    }                                                                   \
}

/*FIXME: the vect != NULL test is a hack - it should not be needed*/

#define VECTOR_DESTRUCT                                                 \
{                                                                       \
    if ( ( v_type == NORMAL_VECTOR ) && size+alloc_lookahead > 0 )      \
    {                                                                   \
        if ( vect != NULL )                                             \
        {                                                               \
            REPDELB(vect);                                              \
        }                                                               \
                                                                        \
        vect = NULL;                                                    \
    }                                                                   \
}

#define VECTOR_MAKE_NORMAL(what_type,loc_look_size)                     \
{                                                                       \
    if ( v_type != NORMAL_VECTOR )                                      \
    {                                                                   \
        THROW_ASSERTB(26,force);                                        \
                                                                        \
        alloc_lookahead = loc_look_size;                                \
        REPNEWB(vect,what_type,loc_look_size);                          \
                                                                        \
        size = 0;                                                       \
                                                                        \
        v_type = NORMAL_VECTOR;                                         \
                                                                        \
        reset_offsets();                                                \
    }                                                                   \
}

#define VECTOR_MAKE_ZERO                                                \
{                                                                       \
    if ( v_type != ZERO_VECTOR )                                        \
    {                                                                   \
        THROW_ASSERTB(27,force);                                        \
                                                                        \
        if ( ( v_type == NORMAL_VECTOR ) && size+alloc_lookahead > 0 )  \
        {                                                               \
            REPDELB(vect);                                              \
                                                                        \
            vect = NULL;                                                \
        }                                                               \
                                                                        \
        size = ZERO_VECTOR_SIZE;                                        \
                                                                        \
        v_type = ZERO_VECTOR;                                           \
    }                                                                   \
}

#define VECTOR_MAKE_SCALAR(loc_one)                                     \
{                                                                       \
    if ( v_type != SCALAR_VECTOR )                                      \
    {                                                                   \
        THROW_ASSERTB(28,force);                                        \
                                                                        \
        if ( ( v_type == NORMAL_VECTOR ) && size+alloc_lookahead > 0 )  \
        {                                                               \
            REPDELB(vect);                                              \
                                                                        \
            vect = NULL;                                                \
        }                                                               \
                                                                        \
        size = SCALAR_VECTOR_SIZE;                                      \
                                                                        \
        v_type = SCALAR_VECTOR;                                         \
                                                                        \
        vector_m_diag = loc_one;                                        \
    }                                                                   \
}

#define VECTOR_MAKE_SELECT(loc_one)                                     \
{                                                                       \
    if ( v_type != SELECT_VECTOR )                                      \
    {                                                                   \
        THROW_ASSERTB(29,force);                                        \
                                                                        \
        if ( ( v_type == NORMAL_VECTOR ) && size+alloc_lookahead > 0 )  \
        {                                                               \
            REPDELB(vect);                                              \
                                                                        \
            vect = NULL;                                                \
        }                                                               \
                                                                        \
        size = SELECT_VECTOR_SIZE;                                      \
                                                                        \
        v_type = SELECT_VECTOR;                                         \
                                                                        \
        vector_m_diag = loc_one;                                        \
                                                                        \
        sel_elm = 1;                                                    \
    }                                                                   \
}

#define VECTOR_ASSIGNMENT(what_type,loc_zero,loc_zero_right,loc_look_size) \
{                                                                       \
    long jsq;                                                           \
    long jjj;                                                           \
    long res_size_left;                                                 \
    long res_size_right;                                                \
                                                                        \
    if ( ((void *) this) != ((void *) (&right_op)) )                    \
    {                                                                   \
        {                                                               \
            {                                                           \
                switch ( v_type )                                       \
                {                                                       \
                    case NORMAL_VECTOR:                                 \
                    {                                                   \
                        normal_here:                                    \
                                                                        \
                        jsq = right_op.get_effective_size();            \
                                                                        \
                        if ( size == 0 )                                \
                        {                                               \
                            if ( jsq < 0 )                              \
                            {                                           \
                                switch ( right_op.v_type )              \
                                {                                       \
                                    case NORMAL_VECTOR:                 \
                                    {                                   \
                                        L_THROW(1);                     \
                                                                        \
                                        break;                          \
                                    }                                   \
                                                                        \
                                    case ZERO_VECTOR:                   \
                                    {                                   \
                                        make_zero(DO_FORCE);            \
                                                                        \
                                        break;                          \
                                    }                                   \
                                                                        \
                                    case SCALAR_VECTOR:                 \
                                    {                                   \
                                        make_scalar(DO_FORCE);          \
                                                                        \
                                        break;                          \
                                    }                                   \
                                                                        \
                                    case SELECT_VECTOR:                 \
                                    {                                   \
                                        make_select(DO_FORCE);          \
                                                                        \
                                        sel_elm = right_op.get_sel_elm(); \
                                                                        \
                                        break;                          \
                                    }                                   \
                                                                        \
                                    default:                            \
                                    {                                   \
                                        L_THROW(2);                     \
                                                                        \
                                        break;                          \
                                    }                                   \
                                }                                       \
                            }                                           \
                                                                        \
                            else                                        \
                            {                                           \
                                if ( jsq > 0 )                          \
                                {                                       \
                                    if ( alloc_lookahead < jsq )        \
                                    {                                   \
                                        if ( alloc_lookahead > 0 )      \
                                        {                               \
                                            REPDELB(vect);              \
                                                                        \
                                            vect = NULL;                \
                                        }                               \
                                                                        \
                                        size = jsq;                     \
                                                                        \
                                        alloc_lookahead = loc_look_size; \
                                                                        \
                                        REPNEWB(vect,what_type,(size+loc_look_size)); \
                                    }                                   \
                                                                        \
                                    else                                \
                                    {                                   \
                                        alloc_lookahead -= jsq;         \
                                                                        \
                                        size = jsq;                     \
                                    }                                   \
                                                                        \
                                    reset_offsets();                    \
                                }                                       \
                                                                        \
                                else                                    \
                                {                                       \
                                    if ( alloc_lookahead > 0 )          \
                                    {                                   \
                                        REPDELB(vect);                  \
                                                                        \
                                        vect = NULL;                    \
                                    }                                   \
                                }                                       \
                            }                                           \
                        }                                               \
                                                                        \
                        break;                                          \
                    }                                                   \
                                                                        \
                    case ZERO_VECTOR:                                   \
                    case SCALAR_VECTOR:                                 \
                    case SELECT_VECTOR:                                 \
                    {                                                   \
                        switch ( right_op.v_type )                      \
                        {                                               \
                            case NORMAL_VECTOR:                         \
                            {                                           \
                                make_normal(DO_FORCE);                  \
                                                                        \
                                goto normal_here;                       \
                                                                        \
                                break;                                  \
                            }                                           \
                                                                        \
                            case ZERO_VECTOR:                           \
                            {                                           \
                                make_zero(DO_FORCE);                    \
                                                                        \
                                break;                                  \
                            }                                           \
                                                                        \
                            case SCALAR_VECTOR:                         \
                            {                                           \
                                make_scalar(DO_FORCE);                  \
                                                                        \
                                break;                                  \
                            }                                           \
                                                                        \
                            case SELECT_VECTOR:                         \
                            {                                           \
                                make_select(DO_FORCE);                  \
                                                                        \
                                sel_elm = right_op.get_sel_elm();       \
                                                                        \
                                break;                                  \
                            }                                           \
                                                                        \
                            default:                                    \
                            {                                           \
                                L_THROW(4);                             \
                                                                        \
                                break;                                  \
                            }                                           \
                        }                                               \
                                                                        \
                        break;                                          \
                    }                                                   \
                                                                        \
                    default:                                            \
                    {                                                   \
                        L_THROW(7);                                     \
                                                                        \
                        break;                                          \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            switch ( v_type )                                           \
            {                                                           \
                case NORMAL_VECTOR:                                     \
                {                                                       \
                    res_size_left  = get_effective_size();              \
                    res_size_right = right_op.get_effective_size();     \
                                                                        \
                    if ( res_size_right < 0 )                           \
                    {                                                   \
                        res_size_right = res_size_left;                 \
                    }                                                   \
                                                                        \
                    THROW_ASSERTB(30,res_size_left == res_size_right);  \
                                                                        \
                    if ( res_size_left > 0 )                            \
                    {                                                   \
                        for ( jjj = 1 ; jjj <= res_size_left ; jjj++ )  \
                        {                                               \
                            (*this)[jjj-1] = right_op.get_offset_element(jjj-1); \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case ZERO_VECTOR:                                       \
                {                                                       \
                    switch ( right_op.v_type )                          \
                    {                                                   \
                        case NORMAL_VECTOR:                             \
                        {                                               \
                            L_THROW(3);                                 \
                                                                        \
                            break;                                      \
                        }                                               \
                                                                        \
                        case ZERO_VECTOR:                               \
                        {                                               \
                            break;                                      \
                        }                                               \
                                                                        \
                        case SCALAR_VECTOR:                             \
                        case SELECT_VECTOR:                             \
                        {                                               \
                            THROW_ASSERTB(31,right_op.vector_m_diag == loc_zero_right); \
                                                                        \
                            break;                                      \
                        }                                               \
                                                                        \
                        default:                                        \
                        {                                               \
                            L_THROW(5);                                 \
                                                                        \
                            break;                                      \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case SCALAR_VECTOR:                                     \
                {                                                       \
                    switch ( right_op.v_type )                          \
                    {                                                   \
                        case NORMAL_VECTOR:                             \
                        {                                               \
                            L_THROW(6);                                 \
                                                                        \
                            break;                                      \
                        }                                               \
                                                                        \
                        case ZERO_VECTOR:                               \
                        {                                               \
                            vector_m_diag = loc_zero;                   \
                                                                        \
                            break;                                      \
                        }                                               \
                                                                        \
                        case SCALAR_VECTOR:                             \
                        {                                               \
                            (*this)[-1] = right_op.get_offset_element(-1); \
                                                                        \
                            break;                                      \
                        }                                               \
                                                                        \
                        case SELECT_VECTOR:                             \
                        {                                               \
                            THROW_ASSERTB(32,right_op.vector_m_diag == loc_zero_right); \
                                                                        \
                            vector_m_diag = loc_zero;                   \
                                                                        \
                            break;                                      \
                        }                                               \
                                                                        \
                        default:                                        \
                        {                                               \
                            L_THROW(8);                                 \
                                                                        \
                            break;                                      \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case SELECT_VECTOR:                                     \
                {                                                       \
                    switch ( right_op.v_type )                          \
                    {                                                   \
                        case NORMAL_VECTOR:                             \
                        {                                               \
                            L_THROW(9);                                 \
                                                                        \
                            break;                                      \
                        }                                               \
                                                                        \
                        case ZERO_VECTOR:                               \
                        {                                               \
                            vector_m_diag = loc_zero;                   \
                                                                        \
                            break;                                      \
                        }                                               \
                                                                        \
                        case SCALAR_VECTOR:                             \
                        {                                               \
                            THROW_ASSERTB(33,right_op.vector_m_diag == loc_zero_right); \
                                                                        \
                            vector_m_diag = loc_zero;                   \
                                                                        \
                            break;                                      \
                        }                                               \
                                                                        \
                        case SELECT_VECTOR:                             \
                        {                                               \
                            if ( get_sel_elm() == right_op.get_sel_elm() ) \
                            {                                           \
                                (*this)[-1] = right_op.get_offset_element(-1); \
                            }                                           \
                                                                        \
                            else                                        \
                            {                                           \
                                THROW_ASSERTB(34,right_op.vector_m_diag == loc_zero_right); \
                                                                        \
                                vector_m_diag = loc_zero;               \
                            }                                           \
                                                                        \
                            break;                                      \
                        }                                               \
                                                                        \
                        default:                                        \
                        {                                               \
                            L_THROW(12);                                \
                                                                        \
                            break;                                      \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(13);                                        \
                                                                        \
                    break;                                              \
                }                                                       \
            }                                                           \
        }                                                               \
    }                                                                   \
}

#define VECTOR_SCALAR_ASSIGNMENT                                        \
{                                                                       \
    long jjj;                                                           \
    long res_size_left;                                                 \
                                                                        \
    switch ( v_type )                                                   \
    {                                                                   \
        case NORMAL_VECTOR:                                             \
        {                                                               \
            if ( size == 0 )                                            \
            {                                                           \
                goto scalar_case;                                       \
            }                                                           \
                                                                        \
            res_size_left = get_effective_size();                       \
                                                                        \
            if ( res_size_left > 0 )                                    \
            {                                                           \
                for ( jjj = 1 ; jjj <= res_size_left ; jjj++ )          \
                {                                                       \
                    (*this)[jjj-1] = right_op;                          \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case ZERO_VECTOR:                                               \
        case SCALAR_VECTOR:                                             \
        case SELECT_VECTOR:                                             \
        {                                                               \
            scalar_case:                                                \
                                                                        \
            make_scalar(DO_FORCE);                                      \
                                                                        \
            (*this)[-1] = right_op;                                     \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        default:                                                        \
        {                                                               \
            L_THROW(7);                                                 \
                                                                        \
            break;                                                      \
        }                                                               \
    }                                                                   \
}

#define VECTOR_GET_OFFSET_START                                         \
{                                                                       \
    switch ( v_type )                                                   \
    {                                                                   \
        case NORMAL_VECTOR:                                             \
        {                                                               \
            FN_EXIT_POINT offset_start_effective;                       \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case ZERO_VECTOR:                                               \
        case SCALAR_VECTOR:                                             \
        case SELECT_VECTOR:                                             \
        {                                                               \
            FN_EXIT_POINT 1;                                            \
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
}

#define VECTOR_GET_EFFECTIVE_SIZE                                       \
{                                                                       \
    switch ( v_type )                                                   \
    {                                                                   \
        case NORMAL_VECTOR:                                             \
        {                                                               \
            FN_EXIT_POINT offset_end_effective-offset_start_effective+1; \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case ZERO_VECTOR:                                               \
        case SCALAR_VECTOR:                                             \
        case SELECT_VECTOR:                                             \
        {                                                               \
            FN_EXIT_POINT size;                                         \
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
}

#define VECTOR_REF_DEREF(loc_static_zero_ref,loc_static_m_diag_ref)     \
{                                                                       \
    if ( j_elm == -1 )                                                  \
    {                                                                   \
        FN_EXIT_POINT vector_m_diag;                                    \
    }                                                                   \
                                                                        \
    switch ( v_type )                                                   \
    {                                                                   \
        case NORMAL_VECTOR:                                             \
        {                                                               \
            THROW_ASSERTB(35,j_elm >= 0);                               \
            THROW_ASSERTB(36,j_elm <= (get_effective_size()-1));        \
                                                                        \
            FN_EXIT_POINT vect[j_elm+offset_start_effective-1];         \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SCALAR_VECTOR:                                             \
        {                                                               \
            THROW_ASSERTB(37,j_elm >= 0);                               \
                                                                        \
            FN_EXIT_POINT loc_static_m_diag_ref;                        \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case ZERO_VECTOR:                                               \
        {                                                               \
            THROW_ASSERTB(38,j_elm >= 0);                               \
                                                                        \
            FN_EXIT_POINT loc_static_zero_ref;                          \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SELECT_VECTOR:                                             \
        {                                                               \
            THROW_ASSERTB(39,j_elm >= 0);                               \
                                                                        \
            if ( j_elm == sel_elm-1 )                                   \
            {                                                           \
                FN_EXIT_POINT vector_m_diag;                            \
            }                                                           \
                                                                        \
            else                                                        \
            {                                                           \
                FN_EXIT_POINT loc_static_zero_ref;                      \
            }                                                           \
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
}

#define VECTOR_GET_OFFSET_ELEMENT(loc_zero)                             \
{                                                                       \
    if ( j_elm == -1 )                                                  \
    {                                                                   \
        FN_EXIT_POINT vector_m_diag;                                    \
    }                                                                   \
                                                                        \
    switch ( v_type )                                                   \
    {                                                                   \
        case NORMAL_VECTOR:                                             \
        {                                                               \
            THROW_ASSERTB(40,j_elm >= 0);                               \
            THROW_ASSERTB(41,j_elm <= (get_effective_size()-1));        \
                                                                        \
            FN_EXIT_POINT vect[j_elm+offset_start_effective-1];         \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case ZERO_VECTOR:                                               \
        {                                                               \
            THROW_ASSERTB(42,j_elm >= 0);                               \
                                                                        \
            FN_EXIT_POINT loc_zero;                                     \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SCALAR_VECTOR:                                             \
        {                                                               \
            THROW_ASSERTB(43,j_elm >= 0);                               \
                                                                        \
            FN_EXIT_POINT vector_m_diag;                                \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SELECT_VECTOR:                                             \
        {                                                               \
            THROW_ASSERTB(44,j_elm >= 0);                               \
                                                                        \
            if ( j_elm == sel_elm-1 )                                   \
            {                                                           \
                FN_EXIT_POINT vector_m_diag;                            \
            }                                                           \
                                                                        \
            else                                                        \
            {                                                           \
                FN_EXIT_POINT loc_zero;                                 \
            }                                                           \
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
}

#define VECTOR_BSWAP(what_type)                                         \
{                                                                       \
    long k;                                                             \
                                                                        \
    THROW_ASSERTB(45,i >= 1);                                           \
    THROW_ASSERTB(46,j >= 1);                                           \
                                                                        \
    if ( j < i )                                                        \
    {                                                                   \
        k = j;                                                          \
        j = i;                                                          \
        i = k;                                                          \
    }                                                                   \
                                                                        \
    switch ( v_type )                                                   \
    {                                                                   \
        case NORMAL_VECTOR:                                             \
        {                                                               \
            THROW_ASSERTB(47,i <= size);                                \
            THROW_ASSERTB(48,j <= size);                                \
                                                                        \
            if ( i != j )                                               \
            {                                                           \
                what_type temp;                                         \
                                                                        \
                temp = (*this)[j-1];                                    \
                                                                        \
                for ( k = j ; k > i ; k-- )                             \
                {                                                       \
                    (*this)[k-1] = (*this)[k-2];                        \
                }                                                       \
                                                                        \
                (*this)[i-1] = temp;                                    \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case ZERO_VECTOR:                                               \
        case SCALAR_VECTOR:                                             \
        {                                                               \
            break;                                                      \
        }                                                               \
                                                                        \
        case SELECT_VECTOR:                                             \
        {                                                               \
            if ( i != j )                                               \
            {                                                           \
                if ( j == get_sel_elm() )                               \
                {                                                       \
                    set_sel_elm(i);                                     \
                }                                                       \
                                                                        \
                else                                                    \
                {                                                       \
                    if ( ( get_sel_elm() < j ) && ( get_sel_elm() >= i ) ) \
                    {                                                   \
                        set_sel_elm(get_sel_elm()+1);                   \
                    }                                                   \
                }                                                       \
            }                                                           \
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
}

#define VECTOR_FSWAP(what_type)                                         \
{                                                                       \
    long k;                                                             \
                                                                        \
    THROW_ASSERTB(49,i >= 1);                                           \
    THROW_ASSERTB(50,j >= 1);                                           \
                                                                        \
    if ( j < i )                                                        \
    {                                                                   \
        k = j;                                                          \
        j = i;                                                          \
        i = k;                                                          \
    }                                                                   \
                                                                        \
    switch ( v_type )                                                   \
    {                                                                   \
        case NORMAL_VECTOR:                                             \
        {                                                               \
            THROW_ASSERTB(51,i <= size);                                \
            THROW_ASSERTB(52,j <= size);                                \
                                                                        \
            if ( i != j )                                               \
            {                                                           \
                what_type temp;                                         \
                                                                        \
                temp = (*this)[i-1];                                    \
                                                                        \
                for ( k = i ; k < j ; k++ )                             \
                {                                                       \
                    (*this)[k-1] = (*this)[k];                          \
                }                                                       \
                                                                        \
                (*this)[j-1] = temp;                                    \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case ZERO_VECTOR:                                               \
        case SCALAR_VECTOR:                                             \
        {                                                               \
            break;                                                      \
        }                                                               \
                                                                        \
        case SELECT_VECTOR:                                             \
        {                                                               \
            if ( i != j )                                               \
            {                                                           \
                if ( i == get_sel_elm() )                               \
                {                                                       \
                    set_sel_elm(j);                                     \
                }                                                       \
                                                                        \
                else                                                    \
                {                                                       \
                    if ( ( get_sel_elm() <= j ) && ( get_sel_elm() > i ) ) \
                    {                                                   \
                        set_sel_elm(get_sel_elm()-1);                   \
                    }                                                   \
                }                                                       \
            }                                                           \
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
}

#define VECTOR_SQUARESWAP(what_type)                                    \
{                                                                       \
    long k;                                                             \
                                                                        \
    THROW_ASSERTB(53,i >= 1);                                           \
    THROW_ASSERTB(54,j >= 1);                                           \
                                                                        \
    if ( j < i )                                                        \
    {                                                                   \
        k = i;                                                          \
        i = j;                                                          \
        j = k;                                                          \
    }                                                                   \
                                                                        \
    switch ( v_type )                                                   \
    {                                                                   \
        case NORMAL_VECTOR:                                             \
        {                                                               \
            THROW_ASSERTB(55,i <= size);                                \
            THROW_ASSERTB(56,j <= size);                                \
                                                                        \
            if ( i != j )                                               \
            {                                                           \
                what_type temp;                                         \
                                                                        \
                temp         = (*this)[i-1];                            \
                (*this)[i-1] = (*this)[j-1];                            \
                (*this)[j-1] = temp;                                    \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case ZERO_VECTOR:                                               \
        case SCALAR_VECTOR:                                             \
        {                                                               \
            break;                                                      \
        }                                                               \
                                                                        \
        case SELECT_VECTOR:                                             \
        {                                                               \
            if ( i != j )                                               \
            {                                                           \
                if ( i == get_sel_elm() )                               \
                {                                                       \
                    set_sel_elm(j);                                     \
                }                                                       \
                                                                        \
                else                                                    \
                {                                                       \
                    if ( j == get_sel_elm() )                           \
                    {                                                   \
                        set_sel_elm(i);                                 \
                    }                                                   \
                }                                                       \
            }                                                           \
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
}

#define VECTOR_ADDSTART_VOID(loc_zero)                                  \
{                                                                       \
    switch ( v_type )                                                   \
    {                                                                   \
        case NORMAL_VECTOR:                                             \
        {                                                               \
            addstart(loc_zero);                                         \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case ZERO_VECTOR:                                               \
        {                                                               \
            break;                                                      \
        }                                                               \
                                                                        \
        case SCALAR_VECTOR:                                             \
        {                                                               \
            THROW_ASSERTB(57,vector_m_diag == loc_zero);                \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SELECT_VECTOR:                                             \
        {                                                               \
            set_sel_elm(get_sel_elm()+1);                               \
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
}

#define VECTOR_ADDSTART_WHATTYPE(what_type,loc_zero,loc_look_size)      \
{                                                                       \
    what_type *_vect;                                                   \
    long i;                                                             \
                                                                        \
    switch ( v_type )                                                   \
    {                                                                   \
        case NORMAL_VECTOR:                                             \
        {                                                               \
            if ( alloc_lookahead == 0 )                                 \
            {                                                           \
                REPNEWB(_vect,what_type,(size+1+loc_look_size));        \
                                                                        \
                alloc_lookahead = 1 + loc_look_size;                    \
                                                                        \
                if ( size > 0 )                                         \
                {                                                       \
                    for ( i = 1 ; i <= size ; i++ )                     \
                    {                                                   \
                        _vect[i-1] = vect[i-1];                         \
                    }                                                   \
                                                                        \
                    REPDELB(vect);                                      \
                                                                        \
                    vect = NULL;                                        \
                }                                                       \
                                                                        \
                vect = _vect;                                           \
            }                                                           \
                                                                        \
            offset_end_effective++;                                     \
                                                                        \
            if ( offset_start_effective > 1 )                           \
            {                                                           \
                offset_start_effective++;                               \
            }                                                           \
                                                                        \
            if ( size > 0 )                                             \
            {                                                           \
                for ( i = size ; i >= 1 ; i-- )                         \
                {                                                       \
                    vect[i] = vect[i-1];                                \
                }                                                       \
            }                                                           \
                                                                        \
            size++;                                                     \
                                                                        \
            alloc_lookahead--;                                          \
                                                                        \
            vect[0] = what;                                             \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case ZERO_VECTOR:                                               \
        {                                                               \
            THROW_ASSERTB(1,what == loc_zero);                          \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SCALAR_VECTOR:                                             \
        {                                                               \
            THROW_ASSERTB(1,what == vector_m_diag);                     \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SELECT_VECTOR:                                             \
        {                                                               \
            THROW_ASSERTB(59,what == loc_zero);                         \
                                                                        \
            set_sel_elm(get_sel_elm()+1);                               \
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
}

#define VECTOR_ADDSTART_VECTOR(loc_zero)                                \
{                                                                       \
    long i;                                                             \
                                                                        \
    switch ( v_type )                                                   \
    {                                                                   \
        case NORMAL_VECTOR:                                             \
        {                                                               \
            if ( size > 0 )                                             \
            {                                                           \
                switch ( what.v_type )                                  \
                {                                                       \
                    case NORMAL_VECTOR:                                 \
                    {                                                   \
                        if ( what.get_effective_size() > 0 )            \
                        {                                               \
                            for ( i = what.get_effective_size() ; i >= 1 ; i-- ) \
                            {                                           \
                                addstart(what.get_offset_element(i-1)); \
                            }                                           \
                        }                                               \
                                                                        \
                        break;                                          \
                    }                                                   \
                                                                        \
                    case ZERO_VECTOR:                                   \
                    case SCALAR_VECTOR:                                 \
                    case SELECT_VECTOR:                                 \
                    {                                                   \
                        L_THROW(0);                                     \
                                                                        \
                        break;                                          \
                    }                                                   \
                                                                        \
                    default:                                            \
                    {                                                   \
                        L_THROW(0);                                     \
                                                                        \
                        break;                                          \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            else                                                        \
            {                                                           \
                (*this) = what;                                         \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case ZERO_VECTOR:                                               \
        {                                                               \
            switch ( what.v_type )                                      \
            {                                                           \
                case NORMAL_VECTOR:                                     \
                {                                                       \
                    if ( what.get_effective_size() > 0 )                \
                    {                                                   \
                        for ( i = what.get_effective_size() ; i >= 1 ; i-- ) \
                        {                                               \
                            addstart(what.get_offset_element(i-1));     \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case ZERO_VECTOR:                                       \
                {                                                       \
                    break;                                              \
                }                                                       \
                                                                        \
                case SCALAR_VECTOR:                                     \
                {                                                       \
                    THROW_ASSERTB(60,what.vector_m_diag == loc_zero);   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case SELECT_VECTOR:                                     \
                {                                                       \
                    if ( what.vector_m_diag != loc_zero )               \
                    {                                                   \
                        make_select(DO_FORCE);                          \
                                                                        \
                        set_sel_elm(what.get_sel_elm());                \
                                                                        \
                        vector_m_diag = what.vector_m_diag;             \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(0);                                         \
                                                                        \
                    break;                                              \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SCALAR_VECTOR:                                             \
        {                                                               \
            switch ( what.v_type )                                      \
            {                                                           \
                case NORMAL_VECTOR:                                     \
                {                                                       \
                    if ( what.get_effective_size() > 0 )                \
                    {                                                   \
                        for ( i = what.get_effective_size() ; i >= 1 ; i-- ) \
                        {                                               \
                            addstart(what.get_offset_element(i-1));     \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case ZERO_VECTOR:                                       \
                {                                                       \
                    THROW_ASSERTB(61,vector_m_diag == loc_zero);        \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case SCALAR_VECTOR:                                     \
                {                                                       \
                    THROW_ASSERTB(62,what.vector_m_diag == vector_m_diag); \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case SELECT_VECTOR:                                     \
                {                                                       \
                    if ( ( what.vector_m_diag != loc_zero ) || ( vector_m_diag != loc_zero ) ) \
                    {                                                   \
                        THROW_ASSERTB(63,vector_m_diag == loc_zero);    \
                                                                        \
                        make_select(DO_FORCE);                          \
                                                                        \
                        set_sel_elm(what.get_sel_elm());                \
                                                                        \
                        vector_m_diag = what.vector_m_diag;             \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(0);                                         \
                                                                        \
                    break;                                              \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SELECT_VECTOR:                                             \
        {                                                               \
            switch ( what.v_type )                                      \
            {                                                           \
                case NORMAL_VECTOR:                                     \
                {                                                       \
                    if ( what.get_effective_size() > 0 )                \
                    {                                                   \
                        for ( i = what.get_effective_size() ; i >= 1 ; i-- ) \
                        {                                               \
                            addstart(what.get_offset_element(i-1));     \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case ZERO_VECTOR:                                       \
                {                                                       \
                    THROW_ASSERTB(64,vector_m_diag == loc_zero);        \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case SCALAR_VECTOR:                                     \
                {                                                       \
                    THROW_ASSERTB(65,( vector_m_diag == loc_zero ) && ( what.vector_m_diag == loc_zero )); \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case SELECT_VECTOR:                                     \
                {                                                       \
                    THROW_ASSERTB(81,vector_m_diag == loc_zero);        \
                                                                        \
                    set_sel_elm(what.get_sel_elm());                    \
                                                                        \
                    vector_m_diag = what.vector_m_diag;                 \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(0);                                         \
                                                                        \
                    break;                                              \
                }                                                       \
            }                                                           \
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
}

#define VECTOR_ADDEND_VOID(loc_zero)                                    \
{                                                                       \
    switch ( v_type )                                                   \
    {                                                                   \
        case NORMAL_VECTOR:                                             \
        {                                                               \
            addend(loc_zero);                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case ZERO_VECTOR:                                               \
        {                                                               \
            break;                                                      \
        }                                                               \
                                                                        \
        case SCALAR_VECTOR:                                             \
        {                                                               \
            THROW_ASSERTB(82,vector_m_diag == loc_zero);                \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SELECT_VECTOR:                                             \
        {                                                               \
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
}

#define VECTOR_ADDEND_WHATTYPE(what_type,loc_zero,loc_look_size)        \
{                                                                       \
    what_type *_vect;                                                   \
    long i;                                                             \
                                                                        \
    switch ( v_type )                                                   \
    {                                                                   \
        case NORMAL_VECTOR:                                             \
        {                                                               \
            if ( alloc_lookahead == 0 )                                 \
            {                                                           \
                REPNEWB(_vect,what_type,(size+1+loc_look_size));        \
                                                                        \
                alloc_lookahead = 1 + loc_look_size;                    \
                                                                        \
                if ( size > 0 )                                         \
                {                                                       \
                    for ( i = 1 ; i <= size ; i++ )                     \
                    {                                                   \
                        _vect[i-1] = vect[i-1];                         \
                    }                                                   \
                                                                        \
                    REPDELB(vect);                                      \
                                                                        \
                    vect = NULL;                                        \
                }                                                       \
                                                                        \
                vect = _vect;                                           \
            }                                                           \
                                                                        \
            if ( offset_end_effective == size )                         \
            {                                                           \
                offset_end_effective++;                                 \
            }                                                           \
                                                                        \
            size++;                                                     \
                                                                        \
            alloc_lookahead--;                                          \
                                                                        \
            (vect[size-1]) = (what_type) what;                          \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case ZERO_VECTOR:                                               \
        {                                                               \
            THROW_ASSERTB(84,what == loc_zero);                         \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SCALAR_VECTOR:                                             \
        {                                                               \
            THROW_ASSERTB(85,what == vector_m_diag);                    \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SELECT_VECTOR:                                             \
        {                                                               \
            THROW_ASSERTB(85,what == loc_zero);                         \
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
}

#define VECTOR_ADDEND_VECTOR(loc_zero)                                  \
{                                                                       \
    long i;                                                             \
                                                                        \
    switch ( v_type )                                                   \
    {                                                                   \
        case NORMAL_VECTOR:                                             \
        {                                                               \
            if ( size > 0 )                                             \
            {                                                           \
                switch ( what.v_type )                                  \
                {                                                       \
                    case NORMAL_VECTOR:                                 \
                    {                                                   \
                        if ( what.get_effective_size() > 0 )            \
                        {                                               \
                            for ( i = 1 ; i <= what.get_effective_size() ; i++ ) \
                            {                                           \
                                addend(what.get_offset_element(i-1));   \
                            }                                           \
                        }                                               \
                                                                        \
                        break;                                          \
                    }                                                   \
                                                                        \
                    case ZERO_VECTOR:                                   \
                    {                                                   \
                        L_THROW(0);                                     \
                                                                        \
                        break;                                          \
                    }                                                   \
                                                                        \
                    case SCALAR_VECTOR:                                 \
                    {                                                   \
                        L_THROW(0);                                     \
                                                                        \
                        break;                                          \
                    }                                                   \
                                                                        \
                    case SELECT_VECTOR:                                 \
                    {                                                   \
                        L_THROW(0);                                     \
                                                                        \
                        break;                                          \
                    }                                                   \
                                                                        \
                    default:                                            \
                    {                                                   \
                        L_THROW(0);                                     \
                                                                        \
                        break;                                          \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            else                                                        \
            {                                                           \
                (*this) = what;                                         \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case ZERO_VECTOR:                                               \
        {                                                               \
            switch ( what.v_type )                                      \
            {                                                           \
                case NORMAL_VECTOR:                                     \
                {                                                       \
                    if ( what.get_effective_size() > 0 )                \
                    {                                                   \
                        for ( i = 1 ; i <= what.get_effective_size() ; i++ ) \
                        {                                               \
                            addend(what.get_offset_element(i-1));       \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case ZERO_VECTOR:                                       \
                {                                                       \
                    break;                                              \
                }                                                       \
                                                                        \
                case SCALAR_VECTOR:                                     \
                {                                                       \
                    THROW_ASSERTB(86,what.vector_m_diag == loc_zero);   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case SELECT_VECTOR:                                     \
                {                                                       \
                    THROW_ASSERTB(87,what.vector_m_diag == loc_zero);   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(0);                                         \
                                                                        \
                    break;                                              \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SCALAR_VECTOR:                                             \
        {                                                               \
            switch ( what.v_type )                                      \
            {                                                           \
                case NORMAL_VECTOR:                                     \
                {                                                       \
                    if ( what.get_effective_size() > 0 )                \
                    {                                                   \
                        for ( i = 1 ; i <= what.get_effective_size() ; i++ ) \
                        {                                               \
                            addend(what.get_offset_element(i-1));       \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case ZERO_VECTOR:                                       \
                {                                                       \
                    THROW_ASSERTB(88,vector_m_diag == loc_zero);        \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case SCALAR_VECTOR:                                     \
                {                                                       \
                    THROW_ASSERTB(89,what.vector_m_diag == vector_m_diag); \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case SELECT_VECTOR:                                     \
                {                                                       \
                    THROW_ASSERTB(90,( what.vector_m_diag == loc_zero ) && ( vector_m_diag == loc_zero )); \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(0);                                         \
                                                                        \
                    break;                                              \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case SELECT_VECTOR:                                             \
        {                                                               \
            switch ( what.v_type )                                      \
            {                                                           \
                case NORMAL_VECTOR:                                     \
                {                                                       \
                    if ( what.get_effective_size() > 0 )                \
                    {                                                   \
                        for ( i = 1 ; i <= what.get_effective_size() ; i++ ) \
                        {                                               \
                            addend(what.get_offset_element(i-1));       \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case ZERO_VECTOR:                                       \
                {                                                       \
                    break;                                              \
                }                                                       \
                                                                        \
                case SCALAR_VECTOR:                                     \
                {                                                       \
                    THROW_ASSERTB(91,what.vector_m_diag == loc_zero);   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case SELECT_VECTOR:                                     \
                {                                                       \
                    THROW_ASSERTB(92,what.vector_m_diag == loc_zero);   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(0);                                         \
                                                                        \
                    break;                                              \
                }                                                       \
            }                                                           \
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
}

#define VECTOR_REMOVE                                                   \
{                                                                       \
    long i;                                                             \
                                                                        \
    THROW_ASSERTB(66,which >= 1);                                       \
                                                                        \
    switch ( v_type )                                                   \
    {                                                                   \
        case NORMAL_VECTOR:                                             \
        {                                                               \
            THROW_ASSERTB(67,which <= size);                            \
                                                                        \
            if ( which < size )                                         \
            {                                                           \
                for ( i = which+1 ; i <= size ; i++ )                   \
                {                                                       \
                    vect[i-2] = vect[i-1];                              \
                }                                                       \
            }                                                           \
                                                                        \
            if ( offset_end_effective >= which )                        \
            {                                                           \
                offset_end_effective--;                                 \
            }                                                           \
                                                                        \
            if ( offset_start_effective > which )                       \
            {                                                           \
                offset_start_effective--;                               \
            }                                                           \
                                                                        \
            size--;                                                     \
                                                                        \
            alloc_lookahead++;                                          \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case ZERO_VECTOR:                                               \
        {                                                               \
            break;                                                      \
        }                                                               \
                                                                        \
        case SCALAR_VECTOR:                                             \
        {                                                               \
            break;                                                      \
        }                                                               \
                                                                        \
        case SELECT_VECTOR:                                             \
        {                                                               \
            if ( which < get_sel_elm() )                                \
            {                                                           \
                set_sel_elm(get_sel_elm()-1);                           \
            }                                                           \
                                                                        \
            else                                                        \
            {                                                           \
                if ( which == get_sel_elm() )                           \
                {                                                       \
                    make_zero(DO_FORCE);                                \
                }                                                       \
            }                                                           \
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
}

#define VECTOR_TRIM                                                     \
{                                                                       \
    long i;                                                             \
                                                                        \
    reset_offsets();                                                    \
                                                                        \
    THROW_ASSERTB(200,which >= 0);                                      \
                                                                        \
    switch ( v_type )                                                   \
    {                                                                   \
        case NORMAL_VECTOR:                                             \
        {                                                               \
            if ( size > which )                                         \
            {                                                           \
                for ( i = size ; i > which ; i-- )                      \
                {                                                       \
                    remove(i);                                          \
                }                                                       \
            }                                                           \
                                                                        \
            else                                                        \
            {                                                           \
                if ( size < which )                                     \
                {                                                       \
                    for ( i = size+1 ; i <= which ; i++ )               \
                    {                                                   \
                        addend();                                       \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case ZERO_VECTOR:                                               \
        case SCALAR_VECTOR:                                             \
        case SELECT_VECTOR:                                             \
        {                                                               \
            make_normal(DO_FORCE);                                      \
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
}

#define PAD_VECTOR                                                      \
{                                                                       \
    long i;                                                             \
                                                                        \
    if ( _size == 0 )                                                   \
    {                                                                   \
        FN_EXIT_POINT;                                                  \
    }                                                                   \
                                                                        \
    THROW_ASSERTB(68,_size >= 0);                                       \
                                                                        \
    switch ( v_type )                                                   \
    {                                                                   \
        case ZERO_VECTOR:                                               \
        case SCALAR_VECTOR:                                             \
        case SELECT_VECTOR:                                             \
        {                                                               \
            make_normal(DO_FORCE);                                      \
                                                                        \
            MACRO_COMMENT("there is intentionally not break statement"); \
            MACRO_COMMENT("here, and ordering is quite deliberately the"); \
            MACRO_COMMENT("reverse of what is usually used.");          \
        }                                                               \
                                                                        \
        case NORMAL_VECTOR:                                             \
        {                                                               \
            MACRO_COMMENT("Already know _size > 0");                    \
                                                                        \
            for ( i = 1 ; i <= _size ; i++ )                            \
            {                                                           \
                addend();                                               \
            }                                                           \
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
}






























//=====================================================================
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//=====================================================================
//=====================================================================
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//==*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*==
//===*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*===
//=====================================================================


//
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=                                                                +=+=
// +=+=                      Actual code starts here                   +=+=
// +=+=                                                                +=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
//




//
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=                                                                +=+=
// +=+=                      OPERATOR OVERLOADS                        +=+=
// +=+=                                                                +=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
//

//
// Stream operator overloading
//

std::ostream &operator<<(std::ostream &output, const fVECTOR   &source)
{
    FN_ENTRY_POINT

    VECTOR_OUT;

    FN_EXIT_POINT output;
}

std::istream &operator>>(std::istream &input, fVECTOR   &dest)
{
    FN_ENTRY_POINT

    VECTOR_IN;

    FN_EXIT_POINT input;
}




//
// Mathematical operator overloading
//

// + posation - unary, return rvalue
// - negation - unary, return rvalue

fVECTOR   operator+ (const fVECTOR  &left_op)
{
    FN_ENTRY_POINT

    fVECTOR result;

    UNARY_VECTOR_OPERATION_CONST_RES(L_DOUBLE,0.0,+,);

    FN_EXIT_POINT result;
}


fVECTOR   operator- (const fVECTOR  &left_op)
{
    FN_ENTRY_POINT

    fVECTOR result;

    UNARY_VECTOR_OPERATION_CONST_RES(L_DOUBLE,0.0,-,);

    FN_EXIT_POINT result;
}




// + addition       - binary, return rvalue
// - subtraction    - binary, return rvalue

fVECTOR   operator+ (const fVECTOR  &left_op, const fVECTOR  &right_op)
{
    FN_ENTRY_POINT

    fVECTOR result;

    BINARY_VECTOR_OPERATION_CONST_CONST_RES(L_DOUBLE,0.0,+);

    FN_EXIT_POINT result;
}

fVECTOR   operator+ (const fVECTOR  &left_op, const   double &right_op)
{
    FN_ENTRY_POINT

    fVECTOR result;

    UNARY_VECTOR_OPERATION_CONST_RES(L_DOUBLE,0.0,,+right_op);

    FN_EXIT_POINT result;
}

fVECTOR   operator+ (const fVECTOR  &left_op, const c_double &right_op)
{
    FN_ENTRY_POINT

    fVECTOR result;

    UNARY_VECTOR_OPERATION_CONST_RES(L_DOUBLE,0.0,,+right_op);

    FN_EXIT_POINT result;
}

fVECTOR   operator+ (const   double &right_op, const fVECTOR  &left_op)
{
    FN_ENTRY_POINT

    fVECTOR result;

    UNARY_VECTOR_OPERATION_CONST_RES(L_DOUBLE,0.0,right_op+,);

    FN_EXIT_POINT result;
}

fVECTOR   operator+ (const c_double &right_op, const fVECTOR  &left_op)
{
    FN_ENTRY_POINT

    fVECTOR result;

    UNARY_VECTOR_OPERATION_CONST_RES(L_DOUBLE,0.0,right_op+,);

    FN_EXIT_POINT result;
}


fVECTOR   operator- (const fVECTOR  &left_op, const fVECTOR  &right_op)
{
    FN_ENTRY_POINT

    fVECTOR result;

    BINARY_VECTOR_OPERATION_CONST_CONST_RES(L_DOUBLE,0.0,-);

    FN_EXIT_POINT result;
}

fVECTOR   operator- (const fVECTOR  &left_op, const   double &right_op)
{
    FN_ENTRY_POINT

    fVECTOR result;

    UNARY_VECTOR_OPERATION_CONST_RES(L_DOUBLE,0.0,,-right_op);

    FN_EXIT_POINT result;
}

fVECTOR   operator- (const fVECTOR  &left_op, const c_double &right_op)
{
    FN_ENTRY_POINT

    fVECTOR result;

    UNARY_VECTOR_OPERATION_CONST_RES(L_DOUBLE,0.0,,-right_op);

    FN_EXIT_POINT result;
}

fVECTOR   operator- (const   double &right_op, const fVECTOR  &left_op)
{
    FN_ENTRY_POINT

    fVECTOR result;

    UNARY_VECTOR_OPERATION_CONST_RES(L_DOUBLE,0.0,right_op-,);

    FN_EXIT_POINT result;
}

fVECTOR   operator- (const c_double &right_op, const fVECTOR  &left_op)
{
    FN_ENTRY_POINT

    fVECTOR result;

    UNARY_VECTOR_OPERATION_CONST_RES(L_DOUBLE,0.0,right_op-,);

    FN_EXIT_POINT result;
}




// += additive       assignment - binary, return lvalue
// -= subtractive    assignment - binary, return lvalue

fVECTOR  &operator+=(      fVECTOR  &left_op, const fVECTOR  &right_op)
{
    FN_ENTRY_POINT

    BINARY_VECTOR_OPERATION_VARIA_CONST(L_DOUBLE,0.0,+=);

    FN_EXIT_POINT left_op;
}

fVECTOR  &operator+=(      fVECTOR  &left_op, const   double &right_op)
{
    FN_ENTRY_POINT

    UNARY_VECTOR_OPERATION_VARIA(L_DOUBLE,0.0,,+=right_op);

    FN_EXIT_POINT left_op;
}

fVECTOR  &operator+=(      fVECTOR  &left_op, const c_double &right_op)
{
    FN_ENTRY_POINT

    UNARY_VECTOR_OPERATION_VARIA(L_DOUBLE,0.0,,+=right_op);

    FN_EXIT_POINT left_op;
}


fVECTOR  &operator-=(      fVECTOR  &left_op, const fVECTOR  &right_op)
{
    FN_ENTRY_POINT

    BINARY_VECTOR_OPERATION_VARIA_CONST(L_DOUBLE,0.0,-=);

    FN_EXIT_POINT left_op;
}

fVECTOR  &operator-=(      fVECTOR  &left_op, const   double &right_op)
{
    FN_ENTRY_POINT

    UNARY_VECTOR_OPERATION_VARIA(L_DOUBLE,0.0,,-=right_op);

    FN_EXIT_POINT left_op;
}

fVECTOR  &operator-=(      fVECTOR  &left_op, const c_double &right_op)
{
    FN_ENTRY_POINT

    UNARY_VECTOR_OPERATION_VARIA(L_DOUBLE,0.0,,-=right_op);

    FN_EXIT_POINT left_op;
}




// * multiplication - binary, return rvalue
// / division       - binary, return rvalue

fVECTOR   operator* (const fVECTOR  &left_op, const   double &right_op)
{
    FN_ENTRY_POINT

    fVECTOR result;

    UNARY_VECTOR_OPERATION_CONST_RES(L_DOUBLE,0.0,,*right_op);

    FN_EXIT_POINT result;
}

fVECTOR   operator* (const fVECTOR  &left_op, const c_double &right_op)
{
    FN_ENTRY_POINT

    fVECTOR result;

    UNARY_VECTOR_OPERATION_CONST_RES(L_DOUBLE,0.0,,*right_op);

    FN_EXIT_POINT result;
}

fVECTOR   operator* (const   double &right_op, const fVECTOR  &left_op)
{
    FN_ENTRY_POINT

    fVECTOR result;

    UNARY_VECTOR_OPERATION_CONST_RES(L_DOUBLE,0.0,right_op*,);

    FN_EXIT_POINT result;
}

fVECTOR   operator* (const c_double &right_op, const fVECTOR  &left_op)
{
    FN_ENTRY_POINT

    fVECTOR result;

    UNARY_VECTOR_OPERATION_CONST_RES(L_DOUBLE,0.0,right_op*,);

    FN_EXIT_POINT result;
}



fVECTOR   operator/ (const fVECTOR  &left_op, const   double &right_op)
{
    FN_ENTRY_POINT

    fVECTOR result;

    UNARY_VECTOR_OPERATION_CONST_RES(L_DOUBLE,0.0,,/right_op);

    FN_EXIT_POINT result;
}

fVECTOR   operator/ (const fVECTOR  &left_op, const c_double &right_op)
{
    FN_ENTRY_POINT

    fVECTOR result;

    UNARY_VECTOR_OPERATION_CONST_RES(L_DOUBLE,0.0,,/right_op);

    FN_EXIT_POINT result;
}






// * multiplication - binary, return rvalue

L_DOUBLE  operator* (const fVECTOR  &left_op, const fVECTOR  &right_op)
{
    FN_ENTRY_POINT

    L_DOUBLE result;

    result = 0.0;

    VECTOR_MULT_VECTOR;

    FN_EXIT_POINT result;
}




// *= multiplicative assignment - binary, return lvalue
// /= divisive       assignment - binary, return lvalue

fVECTOR  &operator*=(      fVECTOR  &left_op, const   double &right_op)
{
    FN_ENTRY_POINT

    UNARY_VECTOR_OPERATION_VARIA(L_DOUBLE,0.0,,*=right_op);

    FN_EXIT_POINT left_op;
}

fVECTOR  &operator*=(      fVECTOR  &left_op, const c_double &right_op)
{
    FN_ENTRY_POINT

    UNARY_VECTOR_OPERATION_VARIA(L_DOUBLE,0.0,,*=right_op);

    FN_EXIT_POINT left_op;
}


fVECTOR  &operator/=(      fVECTOR  &left_op, const   double &right_op)
{
    FN_ENTRY_POINT

    UNARY_VECTOR_OPERATION_VARIA(L_DOUBLE,0.0,,/=right_op);

    FN_EXIT_POINT left_op;
}

fVECTOR  &operator/=(      fVECTOR  &left_op, const c_double &right_op)
{
    FN_ENTRY_POINT

    UNARY_VECTOR_OPERATION_VARIA(L_DOUBLE,0.0,,/=right_op);

    FN_EXIT_POINT left_op;
}







//
// Relational operator overloading
//

// == equivalence

int operator==(const fVECTOR  &left_op, const fVECTOR  &right_op)
{
    FN_ENTRY_POINT

    int result;

    VECTOR_VECTOR_RELATION(0.0,&&,==,1);

    FN_EXIT_POINT result;
}

int operator==(const fVECTOR  &left_op, const   double &right_op)
{
    int result;

    VECTOR_SCALAR_RELATION(&&,==,1,0.0);

    FN_EXIT_POINT result;
}

int operator==(const fVECTOR  &left_op, const c_double &right_op)
{
    FN_ENTRY_POINT

    int result;

    VECTOR_SCALAR_RELATION(&&,==,1,0.0);

    FN_EXIT_POINT result;
}

int operator==(const   double &left_op, const fVECTOR  &right_op)
{
    FN_ENTRY_POINT

    int result;

    result = ( right_op == left_op );

    FN_EXIT_POINT result;
}

int operator==(const c_double &left_op, const fVECTOR  &right_op)
{
    FN_ENTRY_POINT

    int result;

    result = ( right_op == left_op );

    FN_EXIT_POINT result;
}





// != inequivalence

int operator!=(const fVECTOR  &left_op, const fVECTOR  &right_op)
{
    FN_ENTRY_POINT

    int result;

    VECTOR_VECTOR_RELATION(0.0,||,!=,0);

    FN_EXIT_POINT result;
}

int operator!=(const fVECTOR  &left_op, const   double &right_op)
{
    FN_ENTRY_POINT

    int result;

    VECTOR_SCALAR_RELATION(||,!=,0,0.0);

    FN_EXIT_POINT result;
}

int operator!=(const fVECTOR  &left_op, const c_double &right_op)
{
    FN_ENTRY_POINT

    int result;

    VECTOR_SCALAR_RELATION(||,!=,0,0.0);

    FN_EXIT_POINT result;
}

int operator!=(const   double &left_op, const fVECTOR  &right_op)
{
    FN_ENTRY_POINT

    int result;

    result = ( right_op != left_op );

    FN_EXIT_POINT result;
}

int operator!=(const c_double &left_op, const fVECTOR  &right_op)
{
    FN_ENTRY_POINT

    int result;

    result = ( right_op != left_op );

    FN_EXIT_POINT result;
}




// < left than

int operator< (const fVECTOR  &left_op, const fVECTOR  &right_op)
{
    FN_ENTRY_POINT

    int result;

    VECTOR_VECTOR_RELATION(0.0,&&,< ,1);

    FN_EXIT_POINT result;
}

int operator< (const fVECTOR  &left_op, const   double &right_op)
{
    FN_ENTRY_POINT

    int result;

    VECTOR_SCALAR_RELATION(&&,< ,1,0.0);

    FN_EXIT_POINT result;
}

int operator< (const fVECTOR  &left_op, const c_double &right_op)
{
    FN_ENTRY_POINT

    int result;

    VECTOR_SCALAR_RELATION(&&,< ,1,0.0);

    FN_EXIT_POINT result;
}

int operator< (const   double &left_op, const fVECTOR  &right_op)
{
    FN_ENTRY_POINT

    int result;

    result = ( right_op > left_op );

    FN_EXIT_POINT result;
}

int operator< (const c_double &left_op, const fVECTOR  &right_op)
{
    FN_ENTRY_POINT

    int result;

    result = ( right_op > left_op );

    FN_EXIT_POINT result;
}




// <= less than or equal to

int operator<=(const fVECTOR  &left_op, const fVECTOR  &right_op)
{
    FN_ENTRY_POINT

    int result;

    VECTOR_VECTOR_RELATION(0.0,&&,<=,1);

    FN_EXIT_POINT result;
}

int operator<=(const fVECTOR  &left_op, const   double &right_op)
{
    FN_ENTRY_POINT

    int result;

    VECTOR_SCALAR_RELATION(&&,<=,1,0.0);

    FN_EXIT_POINT result;
}

int operator<=(const fVECTOR  &left_op, const c_double &right_op)
{
    FN_ENTRY_POINT

    int result;

    VECTOR_SCALAR_RELATION(&&,<=,1,0.0);

    FN_EXIT_POINT result;
}

int operator<=(const   double &left_op, const fVECTOR  &right_op)
{
    FN_ENTRY_POINT

    int result;

    result = ( right_op >= left_op );

    FN_EXIT_POINT result;
}

int operator<=(const c_double &left_op, const fVECTOR  &right_op)
{
    FN_ENTRY_POINT

    int result;

    result = ( right_op >= left_op );

    FN_EXIT_POINT result;
}




// > greater than

int operator> (const fVECTOR  &left_op, const fVECTOR  &right_op)
{
    FN_ENTRY_POINT

    int result;

    VECTOR_VECTOR_RELATION(0.0,&&,> ,1);

    FN_EXIT_POINT result;
}

int operator> (const fVECTOR  &left_op, const   double &right_op)
{
    FN_ENTRY_POINT

    int result;

    VECTOR_SCALAR_RELATION(&&,> ,1,0.0);

    FN_EXIT_POINT result;
}

int operator> (const fVECTOR  &left_op, const c_double &right_op)
{
    FN_ENTRY_POINT

    int result;

    VECTOR_SCALAR_RELATION(&&,> ,1,0.0);

    FN_EXIT_POINT result;
}

int operator> (const   double &left_op, const fVECTOR  &right_op)
{
    FN_ENTRY_POINT

    int result;

    result = ( right_op < left_op );

    FN_EXIT_POINT result;
}

int operator> (const c_double &left_op, const fVECTOR  &right_op)
{
    FN_ENTRY_POINT

    int result;

    result = ( right_op < left_op );

    FN_EXIT_POINT result;
}




// >= greater than or equal to

int operator>= (const fVECTOR  &left_op, const fVECTOR  &right_op)
{
    FN_ENTRY_POINT

    int result;

    VECTOR_VECTOR_RELATION(0.0,&&,>= ,1);

    FN_EXIT_POINT result;
}

int operator>= (const fVECTOR  &left_op, const   double &right_op)
{
    FN_ENTRY_POINT

    int result;

    VECTOR_SCALAR_RELATION(&&,>= ,1,0.0);

    FN_EXIT_POINT result;
}

int operator>= (const fVECTOR  &left_op, const c_double &right_op)
{
    FN_ENTRY_POINT

    int result;

    VECTOR_SCALAR_RELATION(&&,>= ,1,0.0);

    FN_EXIT_POINT result;
}

int operator>= (const   double &left_op, const fVECTOR  &right_op)
{
    FN_ENTRY_POINT

    int result;

    result = ( right_op <= left_op );

    FN_EXIT_POINT result;
}

int operator>= (const c_double &left_op, const fVECTOR  &right_op)
{
    FN_ENTRY_POINT

    int result;

    result = ( right_op <= left_op );

    FN_EXIT_POINT result;
}
















std::ostream &operator<<(std::ostream &output, const fRef_pair &source)
{
    FN_ENTRY_POINT

    output << "vector(" << source.element_num_row << "," << source.element_num << ") = " << fixnum(source.element_val) << "\n";

    FN_EXIT_POINT output;
}
















//
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=                                                                +=+=
// +=+=                       VECTOR CLASS                             +=+=
// +=+=                                                                +=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
//

//
// Constructor and destructor:
//

fVECTOR::fVECTOR()
{
    FN_ENTRY_POINT

    VECTOR_VOID_CONSTRUCTOR(0.0);

    FN_EXIT_POINT;
}

fVECTOR::fVECTOR(char, long _size)
{
    FN_ENTRY_POINT

    GENERIC_VECTOR_CONSTRUCTOR(L_DOUBLE,0.0,1.0);

    FN_EXIT_POINT;
}

fVECTOR::fVECTOR(const fVECTOR &source)
{
    FN_ENTRY_POINT

    GENERIC_VECTOR_COPY_CONSTRUCTOR(L_DOUBLE,0.0);

    FN_EXIT_POINT;
}

fVECTOR::~fVECTOR()
{
    FN_ENTRY_POINT

    VECTOR_DESTRUCT;

    FN_EXIT_POINT;
}




//
// Negation:
//

void fVECTOR::negate_it(void)
{
    FN_ENTRY_POINT

    (*this) *= -1.0;

    FN_EXIT_POINT;
}





//
// Test and modify vector type.
//

void fVECTOR::make_normal(int force)
{
    FN_ENTRY_POINT

    VECTOR_MAKE_NORMAL(L_DOUBLE,LOOK_SIZE);

    FN_EXIT_POINT;

    force = 0;
}

void fVECTOR::make_zero(int force)
{
    FN_ENTRY_POINT

    VECTOR_MAKE_ZERO;

    FN_EXIT_POINT;

    force = 0;
}

void fVECTOR::make_scalar(int force)
{
    FN_ENTRY_POINT

    VECTOR_MAKE_SCALAR(1.0);

    FN_EXIT_POINT;

    force = 0;
}

void fVECTOR::make_select(int force)
{
    FN_ENTRY_POINT

    VECTOR_MAKE_SELECT(1.0);

    FN_EXIT_POINT;

    force = 0;
}

int fVECTOR::is_normal(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT ( v_type == NORMAL_VECTOR );
}

int fVECTOR::is_zero(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT ( v_type == ZERO_VECTOR );
}

int fVECTOR::is_scalar(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT ( v_type == SCALAR_VECTOR );
}

int fVECTOR::is_select(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT ( v_type == SELECT_VECTOR );
}

int fVECTOR::get_v_type(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT v_type;
}





//
// Overloaded assignment operators.
//

//
// Method: - coerce the v_types to compatible types.  A normal vector of
//           size 0 will be coerced to the same type as the source,
//           otherwise it will remain normal with size unchanged.
//         - do the assignment.
//

fVECTOR &fVECTOR::operator=(const fVECTOR &right_op)
{
    FN_ENTRY_POINT

    VECTOR_ASSIGNMENT(L_DOUBLE,0.0,0.0,LOOK_SIZE);

    FN_EXIT_POINT (*this);
}

fVECTOR &fVECTOR::operator=(const   double &right_op)
{
    FN_ENTRY_POINT

    VECTOR_SCALAR_ASSIGNMENT;

    FN_EXIT_POINT (*this);
}

fVECTOR &fVECTOR::operator=(const c_double &right_op)
{
    FN_ENTRY_POINT

    VECTOR_SCALAR_ASSIGNMENT;

    FN_EXIT_POINT (*this);
}




//
// Set and get offsets:
//

void fVECTOR::reset_offsets(void)
{
    FN_ENTRY_POINT

    offset_start_effective = 1;
    offset_end_effective   = size;

    FN_EXIT_POINT;
}

void fVECTOR::set_offset_start(long start)
{
    FN_ENTRY_POINT

    set_offsets(start,offset_end_effective);

    FN_EXIT_POINT;
}

void fVECTOR::set_offset_end(long end)
{
    FN_ENTRY_POINT

    set_offsets(get_offset_start(),end);

    FN_EXIT_POINT;
}

void fVECTOR::set_offsets(long start, long end)
{
    FN_ENTRY_POINT

    THROW_ASSERT(start >= 1);
    THROW_ASSERT(start <= (end+1));
    THROW_ASSERT(end >= (start-1));
    THROW_ASSERT(end <= size);

    offset_start_effective = start;
    offset_end_effective   = end;

    FN_EXIT_POINT;
}

long fVECTOR::get_offset_start(void) const
{
    FN_ENTRY_POINT

    VECTOR_GET_OFFSET_START;

    FN_EXIT_POINT 1;
}

long fVECTOR::get_offset_end(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT offset_end_effective;
}

long fVECTOR::get_effective_size(void) const
{
    FN_ENTRY_POINT

    VECTOR_GET_EFFECTIVE_SIZE;

    // unreachable code

    FN_EXIT_POINT 0;
}

long fVECTOR::get_real_size(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT size;
}




//
// Get standard (T *) version of vector, assuming that this is
// possible in finite dimensions.
//

L_DOUBLE *fVECTOR::operator()(void) const
{
    FN_ENTRY_POINT

    L_DOUBLE *result = NULL;
    long i;

    switch ( v_type )
    {
        case NORMAL_VECTOR:
        {
            REPNEWB(result,L_DOUBLE,(get_effective_size()+1));

            if ( get_effective_size() > 0 )
            {
                for ( i = 1 ; i <= get_effective_size() ; i++ )
                {
                    result[i-1] = get_offset_element(i-1);
                }
            }

            break;
        }

        case ZERO_VECTOR:
        {
            REPNEWB(result,L_DOUBLE,1);

            result[0] = get_offset_element(0);

            break;
        }

        case SCALAR_VECTOR:
        case SELECT_VECTOR:
        {
            REPNEWB(result,L_DOUBLE,1);

            result[0] = get_offset_element(-1);

            break;
        }

        default:
        {
            L_THROW(0);

            break;
        }
    }

    FN_EXIT_POINT result;
}





//
// Array element overloading
//

L_DOUBLE &fVECTOR::operator[](long j_elm)
{
    FN_ENTRY_POINT

    VECTOR_REF_DEREF(l_double_get_static_zero_ref(),l_double_get_static_ref(vector_m_diag));

    // unreachable code

    FN_EXIT_POINT l_double_get_static_zero_ref();
}

L_DOUBLE &fVECTOR::operator()(long j_elm)
{
    FN_ENTRY_POINT

    VECTOR_REF_DEREF(l_double_get_static_zero_ref(),l_double_get_static_ref(vector_m_diag));

    // unreachable code

    FN_EXIT_POINT l_double_get_static_zero_ref();
}




//
// Offset element functions:
//

L_DOUBLE fVECTOR::get_offset_element(long j_elm) const
{
    FN_ENTRY_POINT

    VECTOR_GET_OFFSET_ELEMENT(0.0);

    // unreachable code

    FN_EXIT_POINT 0.0;
}





//
// Misc operator function
//

fRef_pair fVECTOR::operator()(char, Vect_operation what) const
{
    FN_ENTRY_POINT

    THROW_ASSERT((what == FIND_MAX) || (what == FIND_MIN));

    fRef_pair result;
    long jjj;
    long e_size;
    L_DOUBLE temp;

    result.element_num_row = size;
    result.element_val     = 0.0;

    switch ( v_type )
    {
        case NORMAL_VECTOR:
        {
            e_size = get_effective_size();

            result.element_val = 0.0;
            result.element_num = -1;

            if ( what == FIND_MAX )
            {
                if ( e_size > 0 )
                {
                    result.element_num = 1;
                    result.element_val = get_offset_element(1-1);

                    if ( e_size > 1 )
                    {
                        for ( jjj = 2 ; jjj <= e_size ; jjj++ )
                        {
                            if ((result.element_val)< (temp=get_offset_element(jjj-1)))
                            {
                                result.element_num = jjj;
                                result.element_val = temp;
                            }
                        }
                    }
                }
            }

            else
            {
                if ( e_size > 0 )
                {
                    result.element_num = 1;
                    result.element_val = get_offset_element(1-1);

                    if ( e_size > 1 )
                    {
                        for ( jjj = 2 ; jjj <= e_size ; jjj++ )
                        {
                            if ((result.element_val)>(temp=get_offset_element(jjj-1)))
                            {
                                result.element_num = jjj;
                                result.element_val = temp;
                            }
                        }
                    }
                }
            }

            break;
        }

        case ZERO_VECTOR:
        {
            result.element_num = 1;
            result.element_val = 0.0;

            break;
        }

        case SCALAR_VECTOR:
        {
            result.element_num = 1;
            result.element_val = get_offset_element(-1);

            break;
        }

        case SELECT_VECTOR:
        {
            if ( what == FIND_MAX )
            {
                if ( get_offset_element(-1) > 0.0 )
                {
                    result.element_num = get_sel_elm();
                    result.element_val = get_offset_element(-1);
                }

                else
                {
                    result.element_num = get_sel_elm()+1;
                    result.element_val = 0.0;
                }
            }

            else
            {
                if ( get_offset_element(-1) < 0.0 )
                {
                    result.element_num = get_sel_elm();
                    result.element_val = get_offset_element(-1);
                }

                else
                {
                    result.element_num = get_sel_elm()+1;
                    result.element_val = 0.0;
                }
            }

            break;
        }

        default:
        {
            L_THROW(5);

            break;
        }
    }

    FN_EXIT_POINT result;
}




//
// Various swap functions
//

void fVECTOR::bswap(long i, long j)
{
    FN_ENTRY_POINT

    VECTOR_BSWAP(L_DOUBLE);

    FN_EXIT_POINT;
}

void fVECTOR::fswap(long i, long j)
{
    FN_ENTRY_POINT

    VECTOR_FSWAP(L_DOUBLE);

    FN_EXIT_POINT;
}

void fVECTOR::squareswap(long i, long j)
{
    FN_ENTRY_POINT

    VECTOR_SQUARESWAP(L_DOUBLE);

    FN_EXIT_POINT;
}




//
// Add and remove elements
//
// addstart: add element to start of vector
// addend:   add element to end of vector
// remove:   remove element which from vector
//

void fVECTOR::addstart(void)
{
    FN_ENTRY_POINT

    VECTOR_ADDSTART_VOID(0.0);

    FN_EXIT_POINT;
}

void fVECTOR::addstart(const   double &what)
{
    FN_ENTRY_POINT

    VECTOR_ADDSTART_WHATTYPE(L_DOUBLE,0.0,LOOK_SIZE);

    FN_EXIT_POINT;
}

void fVECTOR::addstart(const c_double &what)
{
    FN_ENTRY_POINT

    VECTOR_ADDSTART_WHATTYPE(L_DOUBLE,0.0,LOOK_SIZE);

    FN_EXIT_POINT;
}

void fVECTOR::addstart(const fVECTOR &what)
{
    FN_ENTRY_POINT

    VECTOR_ADDSTART_VECTOR(0.0);

    FN_EXIT_POINT;
}

void fVECTOR::addend(void)
{
    FN_ENTRY_POINT

    VECTOR_ADDEND_VOID(0.0);

    FN_EXIT_POINT;
}

void fVECTOR::addend(const   double &what)
{
    FN_ENTRY_POINT

    VECTOR_ADDEND_WHATTYPE(L_DOUBLE,0.0,LOOK_SIZE);

    FN_EXIT_POINT;
}

void fVECTOR::addend(const c_double &what)
{
    FN_ENTRY_POINT

    VECTOR_ADDEND_WHATTYPE(L_DOUBLE,0.0,LOOK_SIZE);

    FN_EXIT_POINT;
}

void fVECTOR::addend(const fVECTOR &what)
{
    FN_ENTRY_POINT

    VECTOR_ADDEND_VECTOR(0.0);

    FN_EXIT_POINT;
}

void fVECTOR::remove(long which)
{
    FN_ENTRY_POINT

    VECTOR_REMOVE;

    FN_EXIT_POINT;
}

void fVECTOR::trim_to_size(long which)
{
    FN_ENTRY_POINT

    VECTOR_TRIM;

    FN_EXIT_POINT;
}




//
// Miscellaneous:
//

void fVECTOR::pad_vector(long _size)
{
    FN_ENTRY_POINT

    PAD_VECTOR;

    FN_EXIT_POINT;
}




//
// Get and set selector element.
//

void fVECTOR::set_sel_elm(long what)
{
    FN_ENTRY_POINT

    if ( what != get_sel_elm() )
    {
        THROW_ASSERT(what > 0);

        sel_elm = what;
    }

    FN_EXIT_POINT;
}

long fVECTOR::get_sel_elm(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT sel_elm;
}




//
// return the vect pointer.
//

L_DOUBLE *fVECTOR::get_direct_ref(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT vect;
}



//
// IO stuff
//

void fVECTOR::prime_io(std::ostream *_where_to, int _echo_level)
{
    FN_ENTRY_POINT

    where_to   = _where_to;
    echo_level = _echo_level;

    FN_EXIT_POINT;
}










//
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=                                                                +=+=
// +=+=                 VECTOR OF VECTORS CLASS (slow version)         +=+=
// +=+=                                                                +=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
//

//
// Stream operator overloading
//

std::ostream &operator<<(std::ostream &output, const slow_f_fVECTOR   &source)
{
    FN_ENTRY_POINT

    VECTOR_OUT;

    FN_EXIT_POINT output;
}

std::istream &operator>>(std::istream &input, slow_f_fVECTOR   &dest)
{
    FN_ENTRY_POINT

    VECTOR_IN_PROP_PRIME;

    FN_EXIT_POINT input;
}

//
// Constructor and destructor:
//

slow_f_fVECTOR::slow_f_fVECTOR()
{
    FN_ENTRY_POINT

    fVECTOR loc_zero;

    VECTOR_VOID_CONSTRUCTOR(loc_zero);

    FN_EXIT_POINT;
}

slow_f_fVECTOR::slow_f_fVECTOR(char, long _size)
{
    FN_ENTRY_POINT

    fVECTOR loc_zero;
    fVECTOR loc_one;

    GENERIC_VECTOR_CONSTRUCTOR(fVECTOR,loc_zero,loc_one);

    FN_EXIT_POINT;
}

slow_f_fVECTOR::slow_f_fVECTOR(const slow_f_fVECTOR &source)
{
    FN_ENTRY_POINT

    fVECTOR loc_zero;

    GENERIC_VECTOR_COPY_CONSTRUCTOR(fVECTOR,loc_zero);

    FN_EXIT_POINT;
}

slow_f_fVECTOR::~slow_f_fVECTOR()
{
    FN_ENTRY_POINT

    VECTOR_DESTRUCT;

    FN_EXIT_POINT;
}





//
// Test and modify vector type.
//

void slow_f_fVECTOR::make_normal(int force)
{
    FN_ENTRY_POINT

    VECTOR_MAKE_NORMAL(fVECTOR,LOOK_SIZE);

    FN_EXIT_POINT;

    force = 0;
}

void slow_f_fVECTOR::make_zero(int force)
{
    FN_ENTRY_POINT

    VECTOR_MAKE_ZERO;

    FN_EXIT_POINT;

    force = 0;
}

void slow_f_fVECTOR::make_scalar(int force)
{
    FN_ENTRY_POINT

    fVECTOR loc_one;

    VECTOR_MAKE_SCALAR(loc_one);

    FN_EXIT_POINT;

    force = 0;
}

void slow_f_fVECTOR::make_select(int force)
{
    FN_ENTRY_POINT

    fVECTOR loc_one;

    VECTOR_MAKE_SELECT(loc_one);

    FN_EXIT_POINT;

    force = 0;
}

int slow_f_fVECTOR::is_normal(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT ( v_type == NORMAL_VECTOR );
}

int slow_f_fVECTOR::is_zero(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT ( v_type == ZERO_VECTOR );
}

int slow_f_fVECTOR::is_scalar(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT ( v_type == SCALAR_VECTOR );
}

int slow_f_fVECTOR::is_select(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT ( v_type == SELECT_VECTOR );
}

int slow_f_fVECTOR::get_v_type(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT v_type;
}





//
// Overloaded assignment operators.
//

slow_f_fVECTOR &slow_f_fVECTOR::operator=(const slow_f_fVECTOR &right_op)
{
    FN_ENTRY_POINT

    fVECTOR loc_zero;

    VECTOR_ASSIGNMENT(fVECTOR,loc_zero,loc_zero,LOOK_SIZE);

    FN_EXIT_POINT (*this);
}




//
// Set and get offsets:
//

void slow_f_fVECTOR::reset_offsets(void)
{
    FN_ENTRY_POINT

    offset_start_effective = 1;
    offset_end_effective   = size;

    FN_EXIT_POINT;
}

void slow_f_fVECTOR::set_offset_start(long start)
{
    FN_ENTRY_POINT

    set_offsets(start,offset_end_effective);

    FN_EXIT_POINT;
}

void slow_f_fVECTOR::set_offset_end(long end)
{
    FN_ENTRY_POINT

    set_offsets(get_offset_start(),end);

    FN_EXIT_POINT;
}

void slow_f_fVECTOR::set_offsets(long start, long end)
{
    FN_ENTRY_POINT

    THROW_ASSERT(start >= 1);
    THROW_ASSERT(start <= (end+1));
    THROW_ASSERT(end >= (start-1));
    THROW_ASSERT(end <= size);

    offset_start_effective = start;
    offset_end_effective   = end;

    FN_EXIT_POINT;
}

long slow_f_fVECTOR::get_offset_start(void) const
{
    FN_ENTRY_POINT

    VECTOR_GET_OFFSET_START;

    FN_EXIT_POINT 1;
}

long slow_f_fVECTOR::get_offset_end(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT offset_end_effective;
}

long slow_f_fVECTOR::get_effective_size(void) const
{
    FN_ENTRY_POINT

    VECTOR_GET_EFFECTIVE_SIZE;

    // unreachable code

    FN_EXIT_POINT 0;
}

long slow_f_fVECTOR::get_real_size(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT size;
}




//
// Array element overloading
//

fVECTOR &slow_f_fVECTOR::operator[](long j_elm)
{
    FN_ENTRY_POINT

    VECTOR_REF_DEREF(fvector_get_static_zero_ref(),fvector_get_static_ref(vector_m_diag));

    // unreachable code

    FN_EXIT_POINT fvector_get_static_zero_ref();
}

fVECTOR &slow_f_fVECTOR::operator()(long j_elm)
{
    FN_ENTRY_POINT

    VECTOR_REF_DEREF(fvector_get_static_zero_ref(),fvector_get_static_ref(vector_m_diag));

    // unreachable code

    FN_EXIT_POINT fvector_get_static_zero_ref();
}




//
// Offset element functions:
//

fVECTOR slow_f_fVECTOR::get_offset_element(long j_elm) const
{
    FN_ENTRY_POINT

    fVECTOR loc_zero;

    VECTOR_GET_OFFSET_ELEMENT(loc_zero);

    // unreachable code

    FN_EXIT_POINT loc_zero;
}




//
// Various swap functions
//

void slow_f_fVECTOR::bswap(long i, long j)
{
    FN_ENTRY_POINT

    VECTOR_BSWAP(fVECTOR);

    FN_EXIT_POINT;
}

void slow_f_fVECTOR::fswap(long i, long j)
{
    FN_ENTRY_POINT

    VECTOR_FSWAP(fVECTOR);

    FN_EXIT_POINT;
}

void slow_f_fVECTOR::squareswap(long i, long j)
{
    FN_ENTRY_POINT

    VECTOR_SQUARESWAP(fVECTOR);

    FN_EXIT_POINT;
}




//
// Add and remove elements
//
// addstart: add element to start of vector
// addend:   add element to end of vector
// remove:   remove element which from vector
//

void slow_f_fVECTOR::addstart(void)
{
    FN_ENTRY_POINT

    fVECTOR loc_zero;

    VECTOR_ADDSTART_VOID(loc_zero);

    FN_EXIT_POINT;
}

void slow_f_fVECTOR::addstart(const fVECTOR &what)
{
    FN_ENTRY_POINT

    fVECTOR loc_zero;

    VECTOR_ADDSTART_WHATTYPE(fVECTOR,loc_zero,LOOK_SIZE);

    FN_EXIT_POINT;
}

void slow_f_fVECTOR::addstart(const slow_f_fVECTOR &what)
{
    FN_ENTRY_POINT

    fVECTOR loc_zero;

    VECTOR_ADDSTART_VECTOR(loc_zero);

    FN_EXIT_POINT;
}

void slow_f_fVECTOR::addend(void)
{
    FN_ENTRY_POINT

    fVECTOR loc_zero;

    VECTOR_ADDEND_VOID(loc_zero);

    FN_EXIT_POINT;
}

void slow_f_fVECTOR::addend(const fVECTOR &what)
{
    FN_ENTRY_POINT

    fVECTOR loc_zero;

    VECTOR_ADDEND_WHATTYPE(fVECTOR,loc_zero,LOOK_SIZE);

    FN_EXIT_POINT;
}

void slow_f_fVECTOR::addend(const slow_f_fVECTOR &what)
{
    FN_ENTRY_POINT

    fVECTOR loc_zero;

    VECTOR_ADDEND_VECTOR(loc_zero);

    FN_EXIT_POINT;
}

void slow_f_fVECTOR::remove(long which)
{
    FN_ENTRY_POINT

    VECTOR_REMOVE;

    FN_EXIT_POINT;
}

void slow_f_fVECTOR::trim_to_size(long which)
{
    FN_ENTRY_POINT

    VECTOR_TRIM;

    FN_EXIT_POINT;
}




//
// Miscellaneous:
//

void slow_f_fVECTOR::pad_vector(long _size)
{
    FN_ENTRY_POINT

    PAD_VECTOR;

    FN_EXIT_POINT;
}




//
// Get and set selector element.
//

void slow_f_fVECTOR::set_sel_elm(long what)
{
    FN_ENTRY_POINT

    if ( what != get_sel_elm() )
    {
        THROW_ASSERT(what > 0);

        sel_elm = what;
    }

    FN_EXIT_POINT;
}

long slow_f_fVECTOR::get_sel_elm(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT sel_elm;
}




//
// return the vect pointer.
//

fVECTOR *slow_f_fVECTOR::get_direct_ref(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT vect;
}



//
// IO stuff
//

void slow_f_fVECTOR::prime_io(std::ostream *_where_to, int _echo_level)
{
    FN_ENTRY_POINT

    where_to   = _where_to;
    echo_level = _echo_level;

    FN_EXIT_POINT;
}









//
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=                                                                +=+=
// +=+=                 VECTOR OF int     CLASS                        +=+=
// +=+=                                                                +=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
//

//
// Stream operator overloading
//

std::ostream &operator<<(std::ostream &output, const iVECTOR   &source)
{
    FN_ENTRY_POINT

    VECTOR_OUT_INTLONG;

    FN_EXIT_POINT output;
}

std::istream &operator>>(std::istream &input, iVECTOR   &dest)
{
    FN_ENTRY_POINT

    VECTOR_IN;

    FN_EXIT_POINT input;
}

//
// Constructor and destructor:
//

iVECTOR::iVECTOR()
{
    FN_ENTRY_POINT

    VECTOR_VOID_CONSTRUCTOR(0);

    FN_EXIT_POINT;
}

iVECTOR::iVECTOR(char, long _size)
{
    FN_ENTRY_POINT

    GENERIC_VECTOR_CONSTRUCTOR(int,0,1);

    FN_EXIT_POINT;
}

iVECTOR::iVECTOR(const iVECTOR &source)
{
    FN_ENTRY_POINT

    GENERIC_VECTOR_COPY_CONSTRUCTOR(int,0);

    FN_EXIT_POINT;
}

iVECTOR::~iVECTOR()
{
    FN_ENTRY_POINT

    VECTOR_DESTRUCT;

    FN_EXIT_POINT;
}





//
// Test and modify vector type.
//

void iVECTOR::make_normal(int force)
{
    FN_ENTRY_POINT

    VECTOR_MAKE_NORMAL(int,LOOK_SIZE);

    FN_EXIT_POINT;

    force = 0;
}

void iVECTOR::make_zero(int force)
{
    FN_ENTRY_POINT

    VECTOR_MAKE_ZERO;

    FN_EXIT_POINT;

    force = 0;
}

void iVECTOR::make_scalar(int force)
{
    FN_ENTRY_POINT

    VECTOR_MAKE_SCALAR(1);

    FN_EXIT_POINT;

    force = 0;
}

void iVECTOR::make_select(int force)
{
    FN_ENTRY_POINT

    VECTOR_MAKE_SELECT(1);

    FN_EXIT_POINT;

    force = 0;
}

int iVECTOR::is_normal(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT ( v_type == NORMAL_VECTOR );
}

int iVECTOR::is_zero(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT ( v_type == ZERO_VECTOR );
}

int iVECTOR::is_scalar(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT ( v_type == SCALAR_VECTOR );
}

int iVECTOR::is_select(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT ( v_type == SELECT_VECTOR );
}

int iVECTOR::get_v_type(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT v_type;
}





//
// Overloaded assignment operators.
//

iVECTOR &iVECTOR::operator=(const iVECTOR &right_op)
{
    FN_ENTRY_POINT

    VECTOR_ASSIGNMENT(int,0,0,LOOK_SIZE);

    FN_EXIT_POINT (*this);
}




//
// Set and get offsets:
//

void iVECTOR::reset_offsets(void)
{
    FN_ENTRY_POINT

    offset_start_effective = 1;
    offset_end_effective   = size;

    FN_EXIT_POINT;
}

void iVECTOR::set_offset_start(long start)
{
    FN_ENTRY_POINT

    set_offsets(start,offset_end_effective);

    FN_EXIT_POINT;
}

void iVECTOR::set_offset_end(long end)
{
    FN_ENTRY_POINT

    set_offsets(get_offset_start(),end);

    FN_EXIT_POINT;
}

void iVECTOR::set_offsets(long start, long end)
{
    FN_ENTRY_POINT

    THROW_ASSERT(start >= 1);
    THROW_ASSERT(start <= (end+1));
    THROW_ASSERT(end >= (start-1));
    THROW_ASSERT(end <= size);

    offset_start_effective = start;
    offset_end_effective   = end;

    FN_EXIT_POINT;
}

long iVECTOR::get_offset_start(void) const
{
    FN_ENTRY_POINT

    VECTOR_GET_OFFSET_START;

    FN_EXIT_POINT 1;
}

long iVECTOR::get_offset_end(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT offset_end_effective;
}

long iVECTOR::get_effective_size(void) const
{
    FN_ENTRY_POINT

    VECTOR_GET_EFFECTIVE_SIZE;

    // unreachable code

    FN_EXIT_POINT 0;
}

long iVECTOR::get_real_size(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT size;
}



//
// Get standard (T *) version of vector, assuming that this is
// possible in finite dimensions.
//

int *iVECTOR::operator()(void) const
{
    FN_ENTRY_POINT

    int *result = NULL;
    long i;

    switch ( v_type )
    {
        case NORMAL_VECTOR:
        {
            REPNEWB(result,int,(get_effective_size()+1));

            if ( get_effective_size() > 0 )
            {
                for ( i = 1 ; i <= get_effective_size() ; i++ )
                {
                    result[i-1] = get_offset_element(i-1);
                }
            }

            break;
        }

        case ZERO_VECTOR:
        {
            REPNEWB(result,int,1);

            result[0] = get_offset_element(0);

            break;
        }

        case SCALAR_VECTOR:
        case SELECT_VECTOR:
        {
            REPNEWB(result,int,1);

            result[0] = get_offset_element(-1);

            break;
        }

        default:
        {
            L_THROW(0);

            break;
        }
    }

    FN_EXIT_POINT result;
}




//
// Array element overloading
//

int &iVECTOR::operator[](long j_elm)
{
    FN_ENTRY_POINT

    VECTOR_REF_DEREF(int_get_static_zero_ref(),int_get_static_ref(vector_m_diag));

    // unreachable code

    FN_EXIT_POINT int_get_static_zero_ref();
}

int &iVECTOR::operator()(long j_elm)
{
    FN_ENTRY_POINT

    VECTOR_REF_DEREF(int_get_static_zero_ref(),int_get_static_ref(vector_m_diag));

    // unreachable code

    FN_EXIT_POINT int_get_static_zero_ref();
}




//
// Offset element functions:
//

int iVECTOR::get_offset_element(long j_elm) const
{
    FN_ENTRY_POINT

    VECTOR_GET_OFFSET_ELEMENT(0);

    // unreachable code

    FN_EXIT_POINT 0;
}




//
// Various swap functions
//

void iVECTOR::bswap(long i, long j)
{
    FN_ENTRY_POINT

    VECTOR_BSWAP(int);

    FN_EXIT_POINT;
}

void iVECTOR::fswap(long i, long j)
{
    FN_ENTRY_POINT

    VECTOR_FSWAP(int);

    FN_EXIT_POINT;
}

void iVECTOR::squareswap(long i, long j)
{
    FN_ENTRY_POINT

    VECTOR_SQUARESWAP(int);

    FN_EXIT_POINT;
}




//
// Add and remove elements
//
// addstart: add element to start of vector
// addend:   add element to end of vector
// remove:   remove element which from vector
//

void iVECTOR::addstart(void)
{
    FN_ENTRY_POINT

    VECTOR_ADDSTART_VOID(0);

    FN_EXIT_POINT;
}

void iVECTOR::addstart(const int &what)
{
    FN_ENTRY_POINT

    VECTOR_ADDSTART_WHATTYPE(int,0,LOOK_SIZE);

    FN_EXIT_POINT;
}

void iVECTOR::addstart(const iVECTOR &what)
{
    FN_ENTRY_POINT

    VECTOR_ADDSTART_VECTOR(0);

    FN_EXIT_POINT;
}

void iVECTOR::addend(void)
{
    FN_ENTRY_POINT

    VECTOR_ADDEND_VOID(0);

    FN_EXIT_POINT;
}

void iVECTOR::addend(const int &what)
{
    FN_ENTRY_POINT

    VECTOR_ADDEND_WHATTYPE(int,0,LOOK_SIZE);

    FN_EXIT_POINT;
}

void iVECTOR::addend(const iVECTOR &what)
{
    FN_ENTRY_POINT

    VECTOR_ADDEND_VECTOR(0);

    FN_EXIT_POINT;
}

void iVECTOR::remove(long which)
{
    FN_ENTRY_POINT

    VECTOR_REMOVE;

    FN_EXIT_POINT;
}

void iVECTOR::trim_to_size(long which)
{
    FN_ENTRY_POINT

    VECTOR_TRIM;

    FN_EXIT_POINT;
}




//
// Miscellaneous:
//

void iVECTOR::pad_vector(long _size)
{
    FN_ENTRY_POINT

    PAD_VECTOR;

    FN_EXIT_POINT;
}




//
// Get and set selector element.
//

void iVECTOR::set_sel_elm(long what)
{
    FN_ENTRY_POINT

    if ( what != get_sel_elm() )
    {
        THROW_ASSERT(what > 0);

        sel_elm = what;
    }

    FN_EXIT_POINT;
}

long iVECTOR::get_sel_elm(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT sel_elm;
}





//
// IO stuff
//

void iVECTOR::prime_io(std::ostream *_where_to, int _echo_level)
{
    FN_ENTRY_POINT

    where_to   = _where_to;
    echo_level = _echo_level;

    FN_EXIT_POINT;
}






//
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=                                                                +=+=
// +=+=                 VECTOR OF longs   CLASS                        +=+=
// +=+=                                                                +=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
//

//
// Stream operator overloading
//

std::ostream &operator<<(std::ostream &output, const lVECTOR   &source)
{
    FN_ENTRY_POINT

    VECTOR_OUT_INTLONG;

    FN_EXIT_POINT output;
}

std::istream &operator>>(std::istream &input, lVECTOR   &dest)
{
    FN_ENTRY_POINT

    VECTOR_IN;

    FN_EXIT_POINT input;
}


//
// Constructor and destructor:
//

lVECTOR::lVECTOR()
{
    FN_ENTRY_POINT

    VECTOR_VOID_CONSTRUCTOR(0);

    FN_EXIT_POINT;
}

lVECTOR::lVECTOR(char, long _size)
{
    FN_ENTRY_POINT

    GENERIC_VECTOR_CONSTRUCTOR(long,0,1);

    FN_EXIT_POINT;
}

lVECTOR::lVECTOR(const lVECTOR &source)
{
    FN_ENTRY_POINT

    GENERIC_VECTOR_COPY_CONSTRUCTOR(long,0);

    FN_EXIT_POINT;
}

lVECTOR::~lVECTOR()
{
    FN_ENTRY_POINT

    VECTOR_DESTRUCT;

    FN_EXIT_POINT;
}





//
// Test and modify vector type.
//

void lVECTOR::make_normal(int force)
{
    FN_ENTRY_POINT

    VECTOR_MAKE_NORMAL(long,LOOK_SIZE);

    FN_EXIT_POINT;

    force = 0;
}

void lVECTOR::make_zero(int force)
{
    FN_ENTRY_POINT

    VECTOR_MAKE_ZERO;

    FN_EXIT_POINT;

    force = 0;
}

void lVECTOR::make_scalar(int force)
{
    FN_ENTRY_POINT

    VECTOR_MAKE_SCALAR(1);

    FN_EXIT_POINT;

    force = 0;
}

void lVECTOR::make_select(int force)
{
    FN_ENTRY_POINT

    VECTOR_MAKE_SELECT(1);

    FN_EXIT_POINT;

    force = 0;
}

int lVECTOR::is_normal(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT ( v_type == NORMAL_VECTOR );
}

int lVECTOR::is_zero(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT ( v_type == ZERO_VECTOR );
}

int lVECTOR::is_scalar(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT ( v_type == SCALAR_VECTOR );
}

int lVECTOR::is_select(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT ( v_type == SELECT_VECTOR );
}

int lVECTOR::get_v_type(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT v_type;
}





//
// Overloaded assignment operators.
//

lVECTOR &lVECTOR::operator=(const lVECTOR &right_op)
{
    FN_ENTRY_POINT

    VECTOR_ASSIGNMENT(long,0,0,LOOK_SIZE);

    FN_EXIT_POINT (*this);
}




//
// Set and get offsets:
//

void lVECTOR::reset_offsets(void)
{
    FN_ENTRY_POINT

    offset_start_effective = 1;
    offset_end_effective   = size;

    FN_EXIT_POINT;
}

void lVECTOR::set_offset_start(long start)
{
    FN_ENTRY_POINT

    set_offsets(start,offset_end_effective);

    FN_EXIT_POINT;
}

void lVECTOR::set_offset_end(long end)
{
    FN_ENTRY_POINT

    set_offsets(get_offset_start(),end);

    FN_EXIT_POINT;
}

void lVECTOR::set_offsets(long start, long end)
{
    FN_ENTRY_POINT

    THROW_ASSERT(start >= 1);
    THROW_ASSERT(start <= (end+1));
    THROW_ASSERT(end >= (start-1));
    THROW_ASSERT(end <= size);

    offset_start_effective = start;
    offset_end_effective   = end;

    FN_EXIT_POINT;
}

long lVECTOR::get_offset_start(void) const
{
    FN_ENTRY_POINT

    VECTOR_GET_OFFSET_START;

    FN_EXIT_POINT 1;
}

long lVECTOR::get_offset_end(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT offset_end_effective;
}

long lVECTOR::get_effective_size(void) const
{
    FN_ENTRY_POINT

    VECTOR_GET_EFFECTIVE_SIZE;

    // unreachable code

    FN_EXIT_POINT 0;
}

long lVECTOR::get_real_size(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT size;
}



//
// Get standard (T *) version of vector, assuming that this is
// possible in finite dimensions.
//

long *lVECTOR::operator()(void) const
{
    FN_ENTRY_POINT

    long *result = NULL;
    long i;

    switch ( v_type )
    {
        case NORMAL_VECTOR:
        {
            REPNEWB(result,long,(get_effective_size()+1));

            if ( get_effective_size() > 0 )
            {
                for ( i = 1 ; i <= get_effective_size() ; i++ )
                {
                    result[i-1] = get_offset_element(i-1);
                }
            }

            break;
        }

        case ZERO_VECTOR:
        {
            REPNEWB(result,long,1);

            result[0] = get_offset_element(0);

            break;
        }

        case SCALAR_VECTOR:
        case SELECT_VECTOR:
        {
            REPNEWB(result,long,1);

            result[0] = get_offset_element(-1);

            break;
        }

        default:
        {
            L_THROW(0);

            break;
        }
    }

    FN_EXIT_POINT result;
}




//
// Array element overloading
//

long &lVECTOR::operator[](long j_elm)
{
    FN_ENTRY_POINT

    VECTOR_REF_DEREF(long_get_static_zero_ref(),long_get_static_ref(vector_m_diag));

    // unreachable code

    FN_EXIT_POINT long_get_static_zero_ref();
}

long &lVECTOR::operator()(long j_elm)
{
    FN_ENTRY_POINT

    VECTOR_REF_DEREF(long_get_static_zero_ref(),long_get_static_ref(vector_m_diag));

    // unreachable code

    FN_EXIT_POINT long_get_static_zero_ref();
}




//
// Offset element functions:
//

long lVECTOR::get_offset_element(long j_elm) const
{
    FN_ENTRY_POINT

    VECTOR_GET_OFFSET_ELEMENT(0);

    // unreachable code

    FN_EXIT_POINT 0;
}




//
// Various swap functions
//

void lVECTOR::bswap(long i, long j)
{
    FN_ENTRY_POINT

    VECTOR_BSWAP(long);

    FN_EXIT_POINT;
}

void lVECTOR::fswap(long i, long j)
{
    FN_ENTRY_POINT

    VECTOR_FSWAP(long);

    FN_EXIT_POINT;
}

void lVECTOR::squareswap(long i, long j)
{
    FN_ENTRY_POINT

    VECTOR_SQUARESWAP(long);

    FN_EXIT_POINT;
}




//
// Add and remove elements
//
// addstart: add element to start of vector
// addend:   add element to end of vector
// remove:   remove element which from vector
//

void lVECTOR::addstart(void)
{
    FN_ENTRY_POINT

    VECTOR_ADDSTART_VOID(0);

    FN_EXIT_POINT;
}

void lVECTOR::addstart(const long &what)
{
    FN_ENTRY_POINT

    VECTOR_ADDSTART_WHATTYPE(long,0,LOOK_SIZE);

    FN_EXIT_POINT;
}

void lVECTOR::addstart(const lVECTOR &what)
{
    FN_ENTRY_POINT

    VECTOR_ADDSTART_VECTOR(0);

    FN_EXIT_POINT;
}

void lVECTOR::addend(void)
{
    FN_ENTRY_POINT

    VECTOR_ADDEND_VOID(0);

    FN_EXIT_POINT;
}

void lVECTOR::addend(const long &what)
{
    FN_ENTRY_POINT

    VECTOR_ADDEND_WHATTYPE(long,0,LOOK_SIZE);

    FN_EXIT_POINT;
}

void lVECTOR::addend(const lVECTOR &what)
{
    FN_ENTRY_POINT

    VECTOR_ADDEND_VECTOR(0);

    FN_EXIT_POINT;
}

void lVECTOR::remove(long which)
{
    FN_ENTRY_POINT

    VECTOR_REMOVE;

    FN_EXIT_POINT;
}

void lVECTOR::trim_to_size(long which)
{
    FN_ENTRY_POINT

    VECTOR_TRIM;

    FN_EXIT_POINT;
}




//
// Miscellaneous:
//

void lVECTOR::pad_vector(long _size)
{
    FN_ENTRY_POINT

    PAD_VECTOR;

    FN_EXIT_POINT;
}




//
// Get and set selector element.
//

void lVECTOR::set_sel_elm(long what)
{
    FN_ENTRY_POINT

    if ( what != get_sel_elm() )
    {
        THROW_ASSERT(what > 0);

        sel_elm = what;
    }

    FN_EXIT_POINT;
}

long lVECTOR::get_sel_elm(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT sel_elm;
}




//
// IO stuff
//

void lVECTOR::prime_io(std::ostream *_where_to, int _echo_level)
{
    FN_ENTRY_POINT

    where_to   = _where_to;
    echo_level = _echo_level;

    FN_EXIT_POINT;
}









//
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=                                                                +=+=
// +=+=                 VECTOR OF void *'s CLASS                       +=+=
// +=+=                                                                +=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
//



//
// Stream operator overloading
//

std::ostream &operator<<(std::ostream &output, const vpVECTOR   &source)
{
    FN_ENTRY_POINT

    long i;                                                             
                                                                        
    if ( source.size == source.get_effective_size() )        
    {                                                                   
        output << "size: " << source.size << "\n";           
    }                                                                   
                                                                        
    else                                                                
    {                                                                   
        output << "modified size (-10*size)-10: ";                      
        output << -(10+(10*(source.size))) << "\n";          
    }                                                                   
                                                                        
    switch ( source.v_type )                                      
    {                                                                   
        case NORMAL_VECTOR:                                             
        {                                                               
            if ( source.size > 0 )                           
            {                                                           
                if ( source.size != source.get_effective_size() ) 
                {                                                       
                    output << "offset start: " << source.offset_start_effective << "\n";   
                    output << "offset end:   " << source.offset_end_effective   << "\n\n"; 
                }                                                       
                                                                        
                output << "vector \n\n";                                
                                                                        
                THROW_ASSERTB(61,(source.void_point_printer) != NULL);  
                                                                        
                for ( i = 1 ; i <= source.size ; i++ )       
                {                                                       
                    output << i << ": ";                                
                    (*(source.void_point_printer))(output,(source.vect)[i-1]); 
                    output << "\n";                                     
                }                                                       
            }                                                           
                                                                        
            output << "\n";                                             
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case ZERO_VECTOR:                                               
        {                                                               
            break;                                                      
        }                                                               
                                                                        
        case SCALAR_VECTOR:                                             
        {
            L_THROW(0);
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case SELECT_VECTOR:                                             
        {                                                               
            L_THROW(0);
                                                                        
            break;                                                      
        }                                                               
                                                                        
        default:                                                        
        {                                                               
            L_THROW(0);                                                 
                                                                        
            break;                                                      
        }                                                               
    }                                                                   

    FN_EXIT_POINT output;
}

std::istream &operator>>(std::istream &input, vpVECTOR   &dest)
{
    FN_ENTRY_POINT

    long i;                                                             
    wait_dummy howzat;                                                  
    int modified_size = 0;                                              
                                                                        
    int size;                                                           
    int v_type;                                                         
                                                                        
    long offset_start_effective;                                        
    long offset_end_effective;                                          
                                                                        
    std::ostream *where_to;                                             
    int echo_level;                                                     
                                                                        
    where_to   = dest.where_to;                                         
    echo_level = dest.echo_level;                                       
                                                                        
    dest.prime_io(NULL,0);                                              
                                                                        
    if ( where_to == NULL )                                             
    {                                                                   
        dest.make_zero(DO_FORCE);                                       
                                                                        
        input >> howzat; input >> size;                                 
                                                                        
        switch ( size )                                                 
        {                                                               
            case ZERO_VECTOR_SIZE:                                      
            {                                                           
                v_type = ZERO_VECTOR;
                                                                        
                break;                                                  
            }                                                           
                                                                        
            case SCALAR_VECTOR_SIZE:                                    
            {                                                           
                L_THROW(0);
                                                                        
                break;                                                  
            }                                                           
                                                                        
            case SELECT_VECTOR_SIZE:                                    
            {                                                           
                L_THROW(0);
                                                                        
                break;                                                  
            }                                                           
                                                                        
            default:                                                    
            {                                                           
                if ( size <= -10 )                                      
                {                                                       
                    modified_size = 1;                                  
                                                                        
                    size += 10;                                         
                                                                        
                    THROW_ASSERTB(123,(size%10) == 0);                  
                                                                        
                    size /= -10;                                        
                }                                                       
                                                                        
                THROW_ASSERTB(22,size >= 0);                            
                                                                        
                v_type = NORMAL_VECTOR;                                 
                                                                        
                break;                                                  
            }                                                           
        }                                                               
                                                                        
        switch ( v_type )                                               
        {                                                               
            case NORMAL_VECTOR:                                         
            {                                                           
                dest.make_normal(DO_FORCE);                             
                                                                        
                if ( size > 0 )                                         
                {                                                       
                    dest.pad_vector(size);                              
                                                                        
                    if ( modified_size )                                
                    {                                                   
                        input >> howzat; input >> offset_start_effective; 
                        input >> howzat; input >> offset_end_effective; 
                    }                                                   
                                                                        
                    else                                                
                    {                                                   
                        offset_start_effective = 1;                     
                        offset_end_effective   = size;                  
                    }                                                   
                                                                        
                    THROW_ASSERTB(53,(dest.void_point_scanner) != NULL); 
                                                                        
                    for ( i = 1 ; i <= dest.size ; i++ )                
                    {                                                   
                        input >> howzat;                                
                        (*(dest.void_point_scanner))(input,dest[i-1],where_to,echo_level); 
                    }                                                   
                                                                        
                    dest.set_offset_start(offset_start_effective);      
                    dest.set_offset_end(offset_end_effective);          
                }                                                       
                                                                        
                break;                                                  
            }                                                           
                                                                        
            case ZERO_VECTOR:                                           
            {                                                           
                dest.make_zero(DO_FORCE);                               
                                                                        
                break;                                                  
            }                                                           
                                                                        
            case SCALAR_VECTOR:                                         
            {                                                           
                L_THROW(0);
                                                                        
                break;                                                  
            }                                                           
                                                                        
            case SELECT_VECTOR:                                         
            {                                                           
                L_THROW(0);
                                                                        
                break;                                                  
            }                                                           
                                                                        
            default:                                                    
            {                                                           
                L_THROW(0);                                             
                                                                        
                break;                                                  
            }                                                           
        }                                                               
    }                                                                   
                                                                        
    else                                                                
    {                                                                   
        dest.make_zero(DO_FORCE);                                       
                                                                        
        *where_to << "Vector sizes: >= 0 gives normal vector.\n";       
        *where_to << "              -1 gives zero vector.\n\n";           
                                                                        
        *where_to << "Vector size: ";                                   
        input >> size;                                                  
        if ( echo_level ) { *where_to << size << "\n"; }                
                                                                        
        switch ( size )                                                 
        {                                                               
            case ZERO_VECTOR_SIZE:                                      
            {                                                           
                v_type = ZERO_VECTOR;                                   
                                                                        
                break;                                                  
            }                                                           
                                                                        
            case SCALAR_VECTOR_SIZE:                                    
            {                                                           
                L_THROW(0);
                                                                        
                break;                                                  
            }                                                           
                                                                        
            case SELECT_VECTOR_SIZE:                                    
            {                                                           
                L_THROW(0);
                                                                        
                break;                                                  
            }                                                           
                                                                        
            default:                                                    
            {                                                           
                THROW_ASSERTB(23,size >= 0);                            
                                                                        
                v_type = NORMAL_VECTOR;                                 
                                                                        
                break;                                                  
            }                                                           
        }                                                               
                                                                        
        switch ( v_type )                                               
        {                                                               
            case NORMAL_VECTOR:                                         
            {                                                           
                dest.make_normal(DO_FORCE);                             
                                                                        
                if ( size > 0 )                                         
                {                                                       
                    dest.pad_vector(size);                              
                                                                        
                    offset_start_effective = 1;                         
                    offset_end_effective   = size;                      
                                                                        
                    *where_to << "Enter vector:\n";                     
                                                                        
                    THROW_ASSERTB(83,(dest.void_point_scanner) != NULL); 
                                                                        
                    for ( i = 1 ; i <= dest.size ; i++ )                
                    {                                                   
                        (*(dest.void_point_scanner))(input,dest[i-1],where_to,echo_level); 
                    }                                                   
                                                                        
                    dest.set_offset_start(offset_start_effective);      
                    dest.set_offset_end(offset_end_effective);          
                }                                                       
                                                                        
                break;                                                  
            }                                                           
                                                                        
            case ZERO_VECTOR:                                           
            {                                                           
                dest.make_zero(DO_FORCE);                               
                                                                        
                break;                                                  
            }                                                           
                                                                        
            case SCALAR_VECTOR:                                         
            {                                                           
                L_THROW(0);
                                                                        
                break;                                                  
            }                                                           
                                                                        
            case SELECT_VECTOR:                                         
            {                                                           
                L_THROW(0);
                                                                        
                break;                                                  
            }                                                           
                                                                        
            default:                                                    
            {                                                           
                L_THROW(0);                                             
                                                                        
                break;                                                  
            }                                                           
        }                                                               
    }                                                                   

    FN_EXIT_POINT input;
}




void_point_ref::void_point_ref()
{
    FN_ENTRY_POINT

    void_point = NULL;

    FN_EXIT_POINT;
}

void_point_ref::void_point_ref(const void_point_ref &source)
{
    FN_ENTRY_POINT

    void_point = source.void_point;

    FN_EXIT_POINT;
}

void_point_ref::~void_point_ref()
{
    FN_ENTRY_POINT

    FN_EXIT_POINT;
}

void_point_ref &void_point_ref::operator=(const void_point_ref &source)
{
    FN_ENTRY_POINT

    void_point = source.void_point;

    FN_EXIT_POINT *this;
}

int operator==(const void_point_ref &left_op, const void_point_ref &right_op)
{
    FN_ENTRY_POINT

    FN_EXIT_POINT ( left_op.void_point == right_op.void_point );
}

int operator!=(const void_point_ref &left_op, const void_point_ref &right_op)
{
    FN_ENTRY_POINT

    FN_EXIT_POINT ( left_op.void_point != right_op.void_point );
}

void fvector_void_point_printer(std::ostream &output, const void_point_ref &what)
{
    FN_ENTRY_POINT

    THROW_ASSERT(what.void_point != NULL);

    output << *((fVECTOR *) (what.void_point));

    FN_EXIT_POINT;
}

void fvector_void_point_scanner(std::istream &input, void_point_ref &what, std::ostream *where_to, int echo_level)
{
    FN_ENTRY_POINT

    fVECTOR *temp = NULL;

    if ( what.void_point == NULL )
    {
        REPNEW(temp,fVECTOR);

        (what.void_point) = (void *) temp;
    }

    ((fVECTOR *) (what.void_point))->prime_io(where_to,echo_level);

    input >> *((fVECTOR *) (what.void_point));

    FN_EXIT_POINT;
}

void throw_void_point_printer(std::ostream &output, const void_point_ref &what)
{
    FN_ENTRY_POINT

    L_THROW(0);

    // OK - remove this warning ;)

    output << *((fVECTOR *) (what.void_point));

    FN_EXIT_POINT;
}

void throw_void_point_scanner(std::istream &input, void_point_ref &what, std::ostream *where_to, int echo_level)
{
    FN_ENTRY_POINT

    L_THROW(0);

    // Baaad waaarning...

    input >> *((fVECTOR *) (what.void_point));

    where_to   = NULL;
    echo_level = 0;

    FN_EXIT_POINT;
}




//
// Constructor and destructor:
//

vpVECTOR::vpVECTOR()
{
    FN_ENTRY_POINT

    void_point_ref loc_zero;

    VECTOR_VOID_CONSTRUCTOR(loc_zero);

    void_point_printer = throw_void_point_printer;
    void_point_scanner = throw_void_point_scanner;

    FN_EXIT_POINT;
}

vpVECTOR::vpVECTOR(char, long _size)
{
    FN_ENTRY_POINT

    THROW_ASSERT(_size != SCALAR_VECTOR_SIZE);
    THROW_ASSERT(_size != SELECT_VECTOR_SIZE);

    void_point_ref loc_zero;
    void_point_ref loc_one;

    GENERIC_VECTOR_CONSTRUCTOR(void_point_ref,loc_zero,loc_one);

    void_point_printer = throw_void_point_printer;
    void_point_scanner = throw_void_point_scanner;

    FN_EXIT_POINT;
}

vpVECTOR::~vpVECTOR()
{
    FN_ENTRY_POINT

    VECTOR_DESTRUCT;

    FN_EXIT_POINT;
}





//
// Test and modify vector type.
//

void vpVECTOR::make_normal(int force)
{
    FN_ENTRY_POINT

    VECTOR_MAKE_NORMAL(void_point_ref,LOOK_SIZE);

    FN_EXIT_POINT;

    force = 0;
}

void vpVECTOR::make_zero(int force)
{
    FN_ENTRY_POINT

    VECTOR_MAKE_ZERO;

    FN_EXIT_POINT;

    force = 0;
}

int vpVECTOR::is_normal(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT ( v_type == NORMAL_VECTOR );
}

int vpVECTOR::is_zero(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT ( v_type == ZERO_VECTOR );
}

int vpVECTOR::get_v_type(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT v_type;
}





//
// Set and get offsets:
//

void vpVECTOR::reset_offsets(void)
{
    FN_ENTRY_POINT

    offset_start_effective = 1;
    offset_end_effective   = size;

    FN_EXIT_POINT;
}

void vpVECTOR::set_offset_start(long start)
{
    FN_ENTRY_POINT

    set_offsets(start,offset_end_effective);

    FN_EXIT_POINT;
}

void vpVECTOR::set_offset_end(long end)
{
    FN_ENTRY_POINT

    set_offsets(get_offset_start(),end);

    FN_EXIT_POINT;
}

void vpVECTOR::set_offsets(long start, long end)
{
    FN_ENTRY_POINT

    THROW_ASSERT(start >= 1);
    THROW_ASSERT(start <= (end+1));
    THROW_ASSERT(end >= (start-1));
    THROW_ASSERT(end <= size);

    offset_start_effective = start;
    offset_end_effective   = end;

    FN_EXIT_POINT;
}

long vpVECTOR::get_offset_start(void) const
{
    FN_ENTRY_POINT

    VECTOR_GET_OFFSET_START;

    FN_EXIT_POINT 1;
}

long vpVECTOR::get_offset_end(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT offset_end_effective;
}

long vpVECTOR::get_effective_size(void) const
{
    FN_ENTRY_POINT

    VECTOR_GET_EFFECTIVE_SIZE;

    // unreachable code

    FN_EXIT_POINT 0;
}

long vpVECTOR::get_real_size(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT size;
}




//
// Array element overloading
//

void_point_ref &vpVECTOR::operator[](long j_elm)
{
    FN_ENTRY_POINT

    VECTOR_REF_DEREF(void_point_ref_get_static_zero_ref(),void_point_ref_get_static_ref(vector_m_diag));

    // unreachable code

    FN_EXIT_POINT void_point_ref_get_static_zero_ref();
}

void_point_ref &vpVECTOR::operator()(long j_elm)
{
    FN_ENTRY_POINT

    VECTOR_REF_DEREF(void_point_ref_get_static_zero_ref(),void_point_ref_get_static_ref(vector_m_diag));

    // unreachable code

    FN_EXIT_POINT void_point_ref_get_static_zero_ref();
}




//
// Offset element functions:
//

void_point_ref vpVECTOR::get_offset_element(long j_elm) const
{
    FN_ENTRY_POINT

    void_point_ref loc_zero;

    VECTOR_GET_OFFSET_ELEMENT(loc_zero);

    // unreachable code

    FN_EXIT_POINT loc_zero;
}




//
// Various swap functions
//

void vpVECTOR::bswap(long i, long j)
{
    FN_ENTRY_POINT

    VECTOR_BSWAP(void_point_ref);

    FN_EXIT_POINT;
}

void vpVECTOR::fswap(long i, long j)
{
    FN_ENTRY_POINT

    VECTOR_FSWAP(void_point_ref);

    FN_EXIT_POINT;
}

void vpVECTOR::squareswap(long i, long j)
{
    FN_ENTRY_POINT

    VECTOR_SQUARESWAP(void_point_ref);

    FN_EXIT_POINT;
}




//
// Add and remove elements
//
// addstart: add element to start of vector
// addend:   add element to end of vector
// remove:   remove element which from vector
//

void vpVECTOR::addstart(void)
{
    FN_ENTRY_POINT

    void_point_ref loc_zero;

    VECTOR_ADDSTART_VOID(loc_zero);

    FN_EXIT_POINT;
}

void vpVECTOR::addstart(const void_point_ref &what)
{
    FN_ENTRY_POINT

    void_point_ref loc_zero;

    VECTOR_ADDSTART_WHATTYPE(void_point_ref,loc_zero,LOOK_SIZE);

    FN_EXIT_POINT;
}

void vpVECTOR::addend(void)
{
    FN_ENTRY_POINT

    void_point_ref loc_zero;

    VECTOR_ADDEND_VOID(loc_zero);

    FN_EXIT_POINT;
}

void vpVECTOR::addend(const void_point_ref &what)
{
    FN_ENTRY_POINT

    void_point_ref loc_zero;

    VECTOR_ADDEND_WHATTYPE(void_point_ref,loc_zero,LOOK_SIZE);

    FN_EXIT_POINT;
}

void vpVECTOR::remove(long which)
{
    FN_ENTRY_POINT

    VECTOR_REMOVE;

    FN_EXIT_POINT;
}

void vpVECTOR::trim_to_size(long which)
{
    FN_ENTRY_POINT

    VECTOR_TRIM;

    FN_EXIT_POINT;
}




//
// Miscellaneous:
//

void vpVECTOR::pad_vector(long _size)
{
    FN_ENTRY_POINT

    PAD_VECTOR;

    FN_EXIT_POINT;
}




//
// Get and set selector element.
//

void vpVECTOR::set_sel_elm(long what)
{
    FN_ENTRY_POINT

    if ( what != get_sel_elm() )
    {
        THROW_ASSERT(what > 0);

        sel_elm = what;
    }

    FN_EXIT_POINT;
}

long vpVECTOR::get_sel_elm(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT sel_elm;
}




//
// return the vect pointer.
//

void_point_ref *vpVECTOR::get_direct_ref(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT vect;
}


//
// IO stuff
//

void vpVECTOR::prime_io(std::ostream *_where_to, int _echo_level)
{
    FN_ENTRY_POINT

    where_to   = _where_to;
    echo_level = _echo_level;

    FN_EXIT_POINT;
}

void vpVECTOR::set_io_funs(void (*_void_point_printer)(std::ostream &output,const void_point_ref &what), void (*_void_point_scanner)(std::istream &input, void_point_ref &what, std::ostream *where_to, int echo_level))
{
    FN_ENTRY_POINT

    void_point_printer = _void_point_printer;
    void_point_scanner = _void_point_scanner;

    FN_EXIT_POINT;
}





//
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=                                                                +=+=
// +=+=                 VECTOR OF VECTORS CLASS (fast version)         +=+=
// +=+=                                                                +=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
//

//
// Stream operator overloading
//

std::ostream &operator<<(std::ostream &output, const f_fVECTOR &source)
{
    FN_ENTRY_POINT

    output << source.void_base;

    FN_EXIT_POINT output;
}

std::istream &operator>>(std::istream &input, f_fVECTOR &dest)
{
    FN_ENTRY_POINT

    long i;
    long size;

    size = dest.get_real_size();

    if ( size > 0 )
    {
        for ( i = size ; i >= 1 ; i-- )
        {
            dest.remove(i);
        }
    }

    (dest.void_base).prime_io(dest.where_to,dest.echo_level);

    dest.where_to   = NULL;
    dest.echo_level = 0;

    input >> dest.void_base;

    FN_EXIT_POINT input;
}

//
// Constructor and destructor:
//

f_fVECTOR::f_fVECTOR()
{
    FN_ENTRY_POINT

    where_to   = NULL;
    echo_level = 0;

    void_base.make_normal(DO_FORCE);
    void_base.set_io_funs(fvector_void_point_printer,fvector_void_point_scanner);

    FN_EXIT_POINT;
}

f_fVECTOR::f_fVECTOR(char dummy, long _size)
{
    FN_ENTRY_POINT

    long i;
    fVECTOR loc_zero;

    THROW_ASSERT(_size >= 0);

    where_to   = NULL;
    echo_level = 0;

    void_base.make_normal(DO_FORCE);
    void_base.set_io_funs(fvector_void_point_printer,fvector_void_point_scanner);

    if ( _size > 0 )
    {
        for ( i = 1 ; i <= _size ; i++ )
        {
            addend(loc_zero);
        }
    }

    FN_EXIT_POINT;

    dummy = 'x';
}

f_fVECTOR::f_fVECTOR(const f_fVECTOR &source)
{
    FN_ENTRY_POINT

    long i;

    where_to   = NULL;
    echo_level = 0;

    void_base.make_normal(DO_FORCE);
    void_base.set_io_funs(fvector_void_point_printer,fvector_void_point_scanner);

    if ( source.get_real_size() > 0 )
    {
        for ( i = 1 ; i <= source.get_real_size() ; i++ )
        {
            addend(source.get_offset_element(i-1));
        }
    }

    FN_EXIT_POINT;
}

f_fVECTOR::~f_fVECTOR()
{
    FN_ENTRY_POINT

    long i;
    long size;

    size = get_real_size();

    if ( size > 0 )
    {
        for ( i = size ; i >= 1 ; i-- )
        {
            remove(i);
        }
    }

    FN_EXIT_POINT;
}




//
// Overloaded assignment operators.
//

f_fVECTOR &f_fVECTOR::operator=(const f_fVECTOR &right_op)
{
    FN_ENTRY_POINT

    long i;
    long lsize = get_real_size();
    long rsize = right_op.get_real_size();
    fVECTOR loc_zero;

    if ( ((void *) this) != ((void *) (&right_op)) )
    {
        THROW_ASSERT(lsize >= 0);
        THROW_ASSERT(rsize >= 0);

        if ( lsize == 0 )
        {
            if ( rsize > 0 )
            {
                for ( i = 1 ; i <= rsize ; i++ )
                {
                    addend(loc_zero);
                }
            }
        }

        THROW_ASSERT(lsize == rsize);

        if ( rsize > 0 )
        {
            for ( i = 1 ; i <= rsize ; i++ )
            {
                (*this)[i-1] = right_op.get_offset_element(i-1);
            }
        }
    }

    FN_EXIT_POINT *this;
}




//
// Information function
//

long f_fVECTOR::get_real_size(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT void_base.get_real_size();
}

long f_fVECTOR::get_effective_size(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT void_base.get_real_size();
}




//
// Array element overloading
//

fVECTOR &f_fVECTOR::operator[](long j_elm)
{
    FN_ENTRY_POINT

    void *temp;
    fVECTOR *tempb;

    temp = (void_base[j_elm]).void_point;

    THROW_ASSERT(temp != NULL);

    tempb = (fVECTOR *) temp;

    FN_EXIT_POINT *tempb;
}

fVECTOR &f_fVECTOR::operator()(long j_elm)
{
    FN_ENTRY_POINT

    void *temp;
    fVECTOR *tempb;

    temp = (void_base[j_elm]).void_point;

    THROW_ASSERT(temp != NULL);

    tempb = (fVECTOR *) temp;

    FN_EXIT_POINT *tempb;
}




//
// Offset element functions:
//

fVECTOR f_fVECTOR::get_offset_element(long j_elm) const
{
    FN_ENTRY_POINT

    void *temp;
    fVECTOR *tempb;

    temp = (void_base.get_offset_element(j_elm)).void_point;

    THROW_ASSERT(temp != NULL);

    tempb = (fVECTOR *) temp;

    FN_EXIT_POINT *tempb;
}




//
// Various swap functions
//

void f_fVECTOR::bswap(long i, long j)
{
    FN_ENTRY_POINT

    void_base.bswap(i,j);

    FN_EXIT_POINT;
}

void f_fVECTOR::fswap(long i, long j)
{
    FN_ENTRY_POINT

    void_base.fswap(i,j);

    FN_EXIT_POINT;
}

void f_fVECTOR::squareswap(long i, long j)
{
    FN_ENTRY_POINT

    void_base.squareswap(i,j);

    FN_EXIT_POINT;
}




//
// Add and remove elements
//
// addstart: add element to start of vector
// addend:   add element to end of vector
// remove:   remove element which from vector
//

void f_fVECTOR::addstart(const fVECTOR &what)
{
    FN_ENTRY_POINT

    void_point_ref temp;

    REPNEW((temp.void_point),fVECTOR(what));

    void_base.addstart(temp);

    FN_EXIT_POINT;
}

void f_fVECTOR::addend(const fVECTOR &what)
{
    FN_ENTRY_POINT

    void_point_ref temp;

    REPNEW((temp.void_point),fVECTOR(what));

    void_base.addend(temp);

    FN_EXIT_POINT;
}

void f_fVECTOR::remove(long which)
{
    FN_ENTRY_POINT

    REPDEL(((fVECTOR *) (void_base[which-1]).void_point));

    (void_base[which-1]).void_point = NULL;

    void_base.remove(which);

    FN_EXIT_POINT;
}

void f_fVECTOR::trim_to_size(long which)
{
    FN_ENTRY_POINT

    fVECTOR loc_zero;                                                   
    long i;                                                             
    long size = void_base.get_real_size();                              
                                                                        
    THROW_ASSERTB(200,which >= 0);                                      
                                                                        
    if ( size > which )                                                 
    {                                                                   
        for ( i = size ; i > which ; i-- )                              
        {                                                               
            remove(i);                                                  
        }                                                               
    }                                                                   
                                                                        
    else                                                                
    {                                                                   
        if ( size < which )                                             
        {                                                               
            for ( i = size+1 ; i <= which ; i++ )                       
            {                                                           
                addend(loc_zero);                                       
            }                                                           
        }                                                               
    }                                                                   

    FN_EXIT_POINT;
}



//
// IO stuff
//

void f_fVECTOR::prime_io(std::ostream *_where_to, int _echo_level)
{
    FN_ENTRY_POINT

    where_to   = _where_to;
    echo_level = _echo_level;

    FN_EXIT_POINT;
}































//
// Hackery:
//
// Problem: want to return a reference to an element that is constant
// (cannot be changed) in a situation where it is assumed that change will
// not be attempted but, stylistically, a reference is still nicer.  For
// example, if a vector x is a select vector with element i nonzero, and
// x[j-1] (j != i) is asked for, this is a constant, but for simplicity
// it is nice to still be able to use this notation.
//
// Solution: these functions contain static variables.  References may be
// returned to these statics instead of the real thing.  Certain assumptions
// must be made - if too many statics are hanging around in an expression,
// the circular buffer containing all the different values will overflow and
// the solution will be unpredictable.  But then, I tend to avoid expression
// that are so absurdly long (see vector.h for definition of absurdity ;).
// In essence, as long as the expression contains < 100 references to doubles
// no trouble will occur.
//
// Error checking (defining NDEBUG or NO_REF_BUF_CHECK overrides this):
// Two copies of all static data is kept in two separate circular buffers.
// Only one of these is used to provide references to data.  When this
// function is called, each element of this buffer is compared with the
// corresponding element of the second buffer, and if differences are found
// an exception will occur.  Differences will occur if the reference is used
// in a manner which is stylistically non-constant (ie. it is assigned to),
// and so indicates at runtime that the constantness of an element has not
// been respected.
//

L_DOUBLE &l_double_get_static_zero_ref(void)
{
    FN_ENTRY_POINT

    static L_DOUBLE temp_zero = 0.0;

    THROW_ASSERT(temp_zero == 0.0);

    #ifdef DEBUG_STATIC_ZERO
    std::cerr << "static zero reference here\n";
    #endif

    FN_EXIT_POINT temp_zero;
}

L_DOUBLE &l_double_get_static_ref(L_DOUBLE what)
{
    FN_ENTRY_POINT

    static L_DOUBLE *temp_buf = NULL;
    static L_DOUBLE *temp_buf_cmp = NULL;
    static long buf_point = 1;
    long where;
    long i,j,k;

    if ( temp_buf == NULL )
    {
        // Allocate memory

        ZREPNEWB(temp_buf,L_DOUBLE,TEMP_BUF_SIZE_L_DOUBLE);
        ZREPNEWB(temp_buf_cmp,L_DOUBLE,TEMP_BUF_SIZE_L_DOUBLE);

        // Fill look-back buffers

        for ( i = TEMP_BUF_SIZE_L_DOUBLE-CHECK_BUF_DEPTH_L_DOUBLE+1 ; i <= TEMP_BUF_SIZE_L_DOUBLE ; i++ )
        {
            temp_buf[i-1] = 0.0;
            temp_buf_cmp[i-1] = 0.0;
        }

        buf_point = 1;
    }

    // Compare lookback buffers

    j = buf_point-CHECK_BUF_DEPTH_L_DOUBLE;
    k = buf_point-1;

    if ( j <= 0 )
    {
        j += TEMP_BUF_SIZE_L_DOUBLE;
    }

    if ( k <= 0 )
    {
        k += TEMP_BUF_SIZE_L_DOUBLE;
    }

    #ifndef NDEBUG
    #ifndef NO_REF_BUF_CHECK
    {
        if ( j > k )
        {
            for ( i = j ; i <= TEMP_BUF_SIZE_L_DOUBLE ; i++ )
            {
                if ( temp_buf[i-1] != temp_buf_cmp[i-1] )
                {
                    L_THROW(0);
                }
            }

            j = 1;
        }

        for ( i = j ; i <= k ; i++ )
        {
            if ( temp_buf[i-1] != temp_buf_cmp[i-1] )
            {
                L_THROW(0);
            }
        }
    }
    #endif
    #endif

    // Save position

    where = buf_point;

    // Write new data

    temp_buf[where-1] = what;
    temp_buf_cmp[where-1] = what;

    // Update position

    buf_point++;

    if ( buf_point > TEMP_BUF_SIZE_L_DOUBLE )
    {
        buf_point = 1;
    }

    // return reference

    FN_EXIT_POINT temp_buf[where-1];
}

int &int_get_static_zero_ref(void)
{
    FN_ENTRY_POINT

    static int temp_zero = 0;

    THROW_ASSERT(temp_zero == 0);

    FN_EXIT_POINT temp_zero;
}

int &int_get_static_ref(int what)
{
    FN_ENTRY_POINT

    static int *temp_buf = NULL;
    static int *temp_buf_cmp = NULL;
    static long buf_point = 1;
    long where;
    long i,j,k;

    if ( temp_buf == NULL )
    {
        // Allocate memory

        ZREPNEWB(temp_buf,int,TEMP_BUF_SIZE_INT);
        ZREPNEWB(temp_buf_cmp,int,TEMP_BUF_SIZE_INT);

        // Fill look-back buffers

        for ( i = TEMP_BUF_SIZE_INT-CHECK_BUF_DEPTH_INT+1 ; i <= TEMP_BUF_SIZE_INT ; i++ )
        {
            temp_buf[i-1] = 0;
            temp_buf_cmp[i-1] = 0;
        }

        buf_point = 1;
    }

    // Compare lookback buffers

    j = buf_point-CHECK_BUF_DEPTH_INT;
    k = buf_point-1;

    if ( j <= 0 )
    {
        j += TEMP_BUF_SIZE_INT;
    }

    if ( k <= 0 )
    {
        k += TEMP_BUF_SIZE_INT;
    }

    #ifndef NDEBUG
    #ifndef NO_REF_BUF_CHECK
    {
        if ( j > k )
        {
            for ( i = j ; i <= TEMP_BUF_SIZE_INT ; i++ )
            {
                if ( temp_buf[i-1] != temp_buf_cmp[i-1] )
                {
                    L_THROW(0);
                }
            }

            j = 1;
        }

        for ( i = j ; i <= k ; i++ )
        {
            if ( temp_buf[i-1] != temp_buf_cmp[i-1] )
            {
                L_THROW(0);
            }
        }
    }
    #endif
    #endif

    // Save position

    where = buf_point;

    // Write new data

    temp_buf[where-1] = what;
    temp_buf_cmp[where-1] = what;

    // Update position

    buf_point++;

    if ( buf_point > TEMP_BUF_SIZE_INT )
    {
        buf_point = 1;
    }

    // return reference

    FN_EXIT_POINT temp_buf[where-1];
}






fVECTOR &fvector_get_static_zero_ref(void)
{
    FN_ENTRY_POINT

    static fVECTOR temp_zero;

    THROW_ASSERT(temp_zero.is_zero());

    FN_EXIT_POINT temp_zero;
}

fVECTOR &fvector_get_static_ref(fVECTOR what)
{
    FN_ENTRY_POINT

    static fVECTOR *temp_buf = NULL;
    static fVECTOR *temp_buf_cmp = NULL;
    static long buf_point = 1;
    long where;
    long i,j,k;

    if ( temp_buf == NULL )
    {
        // Allocate memory

        ZREPNEWB(temp_buf,fVECTOR,TEMP_BUF_SIZE_FVECTOR);
        ZREPNEWB(temp_buf_cmp,fVECTOR,TEMP_BUF_SIZE_FVECTOR);

        // Fill look-back buffers

        for ( i = TEMP_BUF_SIZE_FVECTOR-CHECK_BUF_DEPTH_FVECTOR+1 ; i <= TEMP_BUF_SIZE_FVECTOR ; i++ )
        {
            temp_buf[i-1] = 0;
            temp_buf_cmp[i-1] = 0;
        }

        buf_point = 1;
    }

    // Compare lookback buffers

    j = buf_point-CHECK_BUF_DEPTH_FVECTOR;
    k = buf_point-1;

    if ( j <= 0 )
    {
        j += TEMP_BUF_SIZE_FVECTOR;
    }

    if ( k <= 0 )
    {
        k += TEMP_BUF_SIZE_FVECTOR;
    }

    #ifndef NDEBUG
    #ifndef NO_REF_BUF_CHECK
    {
        if ( j > k )
        {
            for ( i = j ; i <= TEMP_BUF_SIZE_FVECTOR ; i++ )
            {
                if ( temp_buf[i-1] != temp_buf_cmp[i-1] )
                {
                    L_THROW(0);
                }
            }

            j = 1;
        }

        for ( i = j ; i <= k ; i++ )
        {
            if ( temp_buf[i-1] != temp_buf_cmp[i-1] )
            {
                L_THROW(0);
            }
        }
    }
    #endif
    #endif

    // Save position

    where = buf_point;

    // Write new data

    temp_buf[where-1] = what;
    temp_buf_cmp[where-1] = what;

    // Update position

    buf_point++;

    if ( buf_point > TEMP_BUF_SIZE_FVECTOR )
    {
        buf_point = 1;
    }

    // return reference

    FN_EXIT_POINT temp_buf[where-1];
}




long &long_get_static_zero_ref(void)
{
    FN_ENTRY_POINT

    static long temp_zero = 0;

    THROW_ASSERT(temp_zero == 0);

    FN_EXIT_POINT temp_zero;
}

long &long_get_static_ref(long what)
{
    FN_ENTRY_POINT

    static long *temp_buf = NULL;
    static long *temp_buf_cmp = NULL;
    static long buf_point = 1;
    long where;
    long i,j,k;

    if ( temp_buf == NULL )
    {
        // Allocate memory

        ZREPNEWB(temp_buf,long,TEMP_BUF_SIZE_LONG);
        ZREPNEWB(temp_buf_cmp,long,TEMP_BUF_SIZE_LONG);

        // Fill look-back buffers

        for ( i = TEMP_BUF_SIZE_LONG-CHECK_BUF_DEPTH_LONG+1 ; i <= TEMP_BUF_SIZE_LONG ; i++ )
        {
            temp_buf[i-1] = 0;
            temp_buf_cmp[i-1] = 0;
        }

        buf_point = 1;
    }

    // Compare lookback buffers

    j = buf_point-CHECK_BUF_DEPTH_LONG;
    k = buf_point-1;

    if ( j <= 0 )
    {
        j += TEMP_BUF_SIZE_LONG;
    }

    if ( k <= 0 )
    {
        k += TEMP_BUF_SIZE_LONG;
    }

    #ifndef NDEBUG
    #ifndef NO_REF_BUF_CHECK
    {
        if ( j > k )
        {
            for ( i = j ; i <= TEMP_BUF_SIZE_LONG ; i++ )
            {
                if ( temp_buf[i-1] != temp_buf_cmp[i-1] )
                {
                    L_THROW(0);
                }
            }

            j = 1;
        }

        for ( i = j ; i <= k ; i++ )
        {
            if ( temp_buf[i-1] != temp_buf_cmp[i-1] )
            {
                L_THROW(0);
            }
        }
    }
    #endif
    #endif

    // Save position

    where = buf_point;

    // Write new data

    temp_buf[where-1] = what;
    temp_buf_cmp[where-1] = what;

    // Update position

    buf_point++;

    if ( buf_point > TEMP_BUF_SIZE_LONG )
    {
        buf_point = 1;
    }

    // return reference

    FN_EXIT_POINT temp_buf[where-1];
}






void_point_ref &void_point_ref_get_static_zero_ref(void)
{
    FN_ENTRY_POINT

    static void_point_ref temp_zero;

    THROW_ASSERT(temp_zero.void_point == NULL);

    FN_EXIT_POINT temp_zero;
}

void_point_ref &void_point_ref_get_static_ref(void_point_ref what)
{
    FN_ENTRY_POINT

    L_THROW(0);

    static void_point_ref *temp_buf = NULL;
    static void_point_ref *temp_buf_cmp = NULL;
    static long buf_point = 1;
    long where;
    long i,j,k;

    if ( temp_buf == NULL )
    {
        // Allocate memory

        ZREPNEWB(temp_buf,void_point_ref,TEMP_BUF_SIZE_VPR);
        ZREPNEWB(temp_buf_cmp,void_point_ref,TEMP_BUF_SIZE_VPR);

        // Fill look-back buffers

        for ( i = TEMP_BUF_SIZE_VPR-CHECK_BUF_DEPTH_VPR+1 ; i <= TEMP_BUF_SIZE_VPR ; i++ )
        {
            (temp_buf[i-1]).void_point     = NULL;
            (temp_buf_cmp[i-1]).void_point = NULL;
        }

        buf_point = 1;
    }

    // Compare lookback buffers

    j = buf_point-CHECK_BUF_DEPTH_VPR;
    k = buf_point-1;

    if ( j <= 0 )
    {
        j += TEMP_BUF_SIZE_VPR;
    }

    if ( k <= 0 )
    {
        k += TEMP_BUF_SIZE_VPR;
    }

    #ifndef NDEBUG
    #ifndef NO_REF_BUF_CHECK
    {
        if ( j > k )
        {
            for ( i = j ; i <= TEMP_BUF_SIZE_VPR ; i++ )
            {
                if ( (temp_buf[i-1]).void_point != (temp_buf_cmp[i-1]).void_point )
                {
                    L_THROW(0);
                }
            }

            j = 1;
        }

        for ( i = j ; i <= k ; i++ )
        {
            if ( (temp_buf[i-1]).void_point != (temp_buf_cmp[i-1]).void_point )
            {
                L_THROW(0);
            }
        }
    }
    #endif
    #endif

    // Save position

    where = buf_point;

    // Write new data

    (temp_buf[where-1]).void_point     = what.void_point;
    (temp_buf_cmp[where-1]).void_point = what.void_point;

    // Update position

    buf_point++;

    if ( buf_point > TEMP_BUF_SIZE_VPR )
    {
        buf_point = 1;
    }

    // return reference

    FN_EXIT_POINT temp_buf[where-1];
}






