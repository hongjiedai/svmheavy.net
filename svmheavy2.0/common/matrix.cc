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
// Matrix class.
//
// Written by: Alistair Shilton
//             Melbourne University
//




#include <iostream>
#include <string.h>

#include "matrix.h"
#include "svdefs.h"
#include "search.h"
#include "levicivita.h"
#include "c_double.h"
#include "vector.h"
#include "friends.h"
#include "outfilt.h"


//
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=                                                                +=+=
// +=+=                      OPERATOR OVERLOAD MACROS                  +=+=
// +=+=                                                                +=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
//

#define UNARY_MATRIX_OPERATION_CONST_RES(what_left,what_right)          \
{                                                                       \
    long i,j;                                                           \
    long size_row;                                                      \
    long size_col;                                                      \
    long offset;                                                        \
                                                                        \
    result = left_op;                                                   \
                                                                        \
    MACRO_COMMENT("Make result as symmetric as possible.");             \
                                                                        \
    if ( left_op.get_e_type() == E_SYMM )                               \
    {                                                                   \
        result.make_symmetric(DO_FORCE);                                \
    }                                                                   \
                                                                        \
    if ( ( left_op.get_e_type() == E_EMPTY          ) &&                \
         ( left_op.get_m_type() == SYMMETRIC_MATRIX )    )              \
    {                                                                   \
        result.make_symmetric(DO_FORCE);                                \
    }                                                                   \
                                                                        \
    offset = left_op.get_e_offs();                                      \
                                                                        \
    size_row = left_op.get_effective_height();                          \
    size_col = left_op.get_effective_width();                           \
                                                                        \
    switch ( left_op.get_e_type() )                                     \
    {                                                                   \
        case E_EMPTY:                                                   \
        {                                                               \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_ZERO:                                                    \
        {                                                               \
            L_DOUBLE tempa;                                             \
                                                                        \
            tempa = what_left 0.0 what_right;                           \
                                                                        \
            if ( tempa != 0.0 )                                         \
            {                                                           \
                for ( i = 1 ; i <= size_row ; i++ )                     \
                {                                                       \
                    for ( j = 1 ; j <= size_col ; j++ )                 \
                    {                                                   \
                        result(i-1,j-1) = tempa;                        \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_DIAG:                                                    \
        {                                                               \
            L_DOUBLE tempa;                                             \
                                                                        \
            tempa = what_left 0.0 what_right;                           \
                                                                        \
            if ( tempa != 0.0 )                                         \
            {                                                           \
                for ( i = 1 ; i <= size_row ; i++ )                     \
                {                                                       \
                    for ( j = 1 ; j <= size_col ; j++ )                 \
                    {                                                   \
                        if ( i-offset == j )                            \
                        {                                               \
                            result(i-1,j-1) = what_left (left_op.get_offset_element(i-1,j-1)) what_right; \
                        }                                               \
                                                                        \
                        else                                            \
                        {                                               \
                            result(i-1,j-1) = tempa;                    \
                        }                                               \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            else                                                        \
            {                                                           \
                for ( i = 1 ; i <= size_row ; i++ )                     \
                {                                                       \
                    for ( j = 1 ; j <= size_col ; j++ )                 \
                    {                                                   \
                        if ( i-offset == j )                            \
                        {                                               \
                            result(i-1,j-1) = what_left (left_op.get_offset_element(i-1,j-1)) what_right; \
                        }                                               \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_LOWER:                                                   \
        {                                                               \
            L_DOUBLE tempa;                                             \
                                                                        \
            tempa = what_left 0.0 what_right;                           \
                                                                        \
            if ( tempa != 0.0 )                                         \
            {                                                           \
                for ( i = 1 ; i <= size_row ; i++ )                     \
                {                                                       \
                    for ( j = 1 ; j <= size_col ; j++ )                 \
                    {                                                   \
                        if ( i-offset >= j )                            \
                        {                                               \
                            result(i-1,j-1) = what_left (left_op.get_offset_element(i-1,j-1)) what_right; \
                        }                                               \
                                                                        \
                        else                                            \
                        {                                               \
                            result(i-1,j-1) = tempa;                    \
                        }                                               \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            else                                                        \
            {                                                           \
                for ( i = 1 ; i <= size_row ; i++ )                     \
                {                                                       \
                    for ( j = 1 ; j <= size_col ; j++ )                 \
                    {                                                   \
                        if ( i-offset >= j )                            \
                        {                                               \
                            result(i-1,j-1) = what_left (left_op.get_offset_element(i-1,j-1)) what_right; \
                        }                                               \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_UPPER:                                                   \
        {                                                               \
            L_DOUBLE tempa;                                             \
                                                                        \
            tempa = what_left 0.0 what_right;                           \
                                                                        \
            if ( tempa != 0.0 )                                         \
            {                                                           \
                for ( i = 1 ; i <= size_row ; i++ )                     \
                {                                                       \
                    for ( j = 1 ; j <= size_col ; j++ )                 \
                    {                                                   \
                        if ( i-offset <= j )                            \
                        {                                               \
                            result(i-1,j-1) = what_left (left_op.get_offset_element(i-1,j-1)) what_right; \
                        }                                               \
                                                                        \
                        else                                            \
                        {                                               \
                            result(i-1,j-1) = tempa;                    \
                        }                                               \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            else                                                        \
            {                                                           \
                for ( i = 1 ; i <= size_row ; i++ )                     \
                {                                                       \
                    for ( j = 1 ; j <= size_col ; j++ )                 \
                    {                                                   \
                        if ( i-offset <= j )                            \
                        {                                               \
                            result(i-1,j-1) = what_left (left_op.get_offset_element(i-1,j-1)) what_right; \
                        }                                               \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_SYMM:                                                    \
        {                                                               \
            for ( i = 1 ; i <= size_row ; i++ )                         \
            {                                                           \
                for ( j = 1 ; j <= i ; j++ )                            \
                {                                                       \
                    result(i-1,j-1) = what_left (left_op.get_offset_element(i-1,j-1)) what_right; \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_ASYMM:                                                   \
        {                                                               \
            for ( i = 1 ; i <= size_row ; i++ )                         \
            {                                                           \
                for ( j = 1 ; j <= size_col ; j++ )                     \
                {                                                       \
                    result(i-1,j-1) = what_left (left_op.get_offset_element(i-1,j-1)) what_right; \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_VZ:                                                      \
        {                                                               \
            L_DOUBLE tempa;                                             \
                                                                        \
            tempa = what_left 0.0 what_right;                           \
                                                                        \
            THROW_ASSERTB(10,tempa == 0.0);                             \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_VS:                                                      \
        {                                                               \
            L_DOUBLE tempa;                                             \
                                                                        \
            tempa = what_left 0.0 what_right;                           \
                                                                        \
            THROW_ASSERTB(11,tempa == 0.0);                             \
                                                                        \
            result(-1,-1) = what_left left_op.get_offset_element(-1,-1) what_right; \
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

#define UNARY_MATRIX_OPERATION_VARIA(what_left,what_right)              \
{                                                                       \
    long i,j;                                                           \
    long size_row;                                                      \
    long size_col;                                                      \
    long offset;                                                        \
                                                                        \
    left_op.fix_symm();                                                 \
                                                                        \
    offset = left_op.get_e_offs();                                      \
                                                                        \
    size_row = left_op.get_effective_height();                          \
    size_col = left_op.get_effective_width();                           \
                                                                        \
    switch ( left_op.get_e_type() )                                     \
    {                                                                   \
        case E_EMPTY:                                                   \
        {                                                               \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_ZERO:                                                    \
        {                                                               \
            L_DOUBLE tempb = 0.0;                                       \
                                                                        \
            what_left tempb what_right;                                 \
                                                                        \
            if ( tempb != 0.0 )                                         \
            {                                                           \
                for ( i = 1 ; i <= size_row ; i++ )                     \
                {                                                       \
                    for ( j = 1 ; j <= size_col ; j++ )                 \
                    {                                                   \
                        left_op(i-1,j-1) = tempb;                       \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_DIAG:                                                    \
        {                                                               \
            L_DOUBLE tempb = 0.0;                                       \
                                                                        \
            what_left tempb what_right;                                 \
                                                                        \
            if ( tempb != 0.0 )                                         \
            {                                                           \
                for ( i = 1 ; i <= size_row ; i++ )                     \
                {                                                       \
                    for ( j = 1 ; j <= size_col ; j++ )                 \
                    {                                                   \
                        if ( i-offset == j )                            \
                        {                                               \
                            what_left left_op(i-1,j-1) what_right;      \
                        }                                               \
                                                                        \
                        else                                            \
                        {                                               \
                            left_op(i-1,j-1) = tempb;                   \
                        }                                               \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            else                                                        \
            {                                                           \
                for ( i = 1 ; i <= size_row ; i++ )                     \
                {                                                       \
                    for ( j = 1 ; j <= size_col ; j++ )                 \
                    {                                                   \
                        if ( i-offset == j )                            \
                        {                                               \
                            what_left left_op(i-1,j-1) what_right;      \
                        }                                               \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_LOWER:                                                   \
        {                                                               \
            L_DOUBLE tempb = 0.0;                                       \
                                                                        \
            what_left tempb what_right;                                 \
                                                                        \
            if ( tempb != 0.0 )                                         \
            {                                                           \
                for ( i = 1 ; i <= size_row ; i++ )                     \
                {                                                       \
                    for ( j = 1 ; j <= size_col ; j++ )                 \
                    {                                                   \
                        if ( i-offset >= j )                            \
                        {                                               \
                            what_left left_op(i-1,j-1) what_right;      \
                        }                                               \
                                                                        \
                        else                                            \
                        {                                               \
                            left_op(i-1,j-1) = tempb;                   \
                        }                                               \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            else                                                        \
            {                                                           \
                for ( i = 1 ; i <= size_row ; i++ )                     \
                {                                                       \
                    for ( j = 1 ; j <= size_col ; j++ )                 \
                    {                                                   \
                        if ( i-offset >= j )                            \
                        {                                               \
                            what_left left_op(i-1,j-1) what_right;      \
                        }                                               \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_UPPER:                                                   \
        {                                                               \
            L_DOUBLE tempb = 0.0;                                       \
                                                                        \
            what_left tempb what_right;                                 \
                                                                        \
            if ( tempb != 0.0 )                                         \
            {                                                           \
                for ( i = 1 ; i <= size_row ; i++ )                     \
                {                                                       \
                    for ( j = 1 ; j <= size_col ; j++ )                 \
                    {                                                   \
                        if ( i-offset <= j )                            \
                        {                                               \
                            what_left left_op(i-1,j-1) what_right;      \
                        }                                               \
                                                                        \
                        else                                            \
                        {                                               \
                            left_op(i-1,j-1) = tempb;                   \
                        }                                               \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            else                                                        \
            {                                                           \
                for ( i = 1 ; i <= size_row ; i++ )                     \
                {                                                       \
                    for ( j = 1 ; j <= size_col ; j++ )                 \
                    {                                                   \
                        if ( i-offset <= j )                            \
                        {                                               \
                            what_left left_op(i-1,j-1) what_right;      \
                        }                                               \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_SYMM:                                                    \
        {                                                               \
            for ( i = 1 ; i <= size_row ; i++ )                         \
            {                                                           \
                for ( j = 1 ; j <= i ; j++ )                            \
                {                                                       \
                    what_left left_op(i-1,j-1) what_right;              \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_ASYMM:                                                   \
        {                                                               \
            for ( i = 1 ; i <= size_row ; i++ )                         \
            {                                                           \
                for ( j = 1 ; j <= size_col ; j++ )                     \
                {                                                       \
                    what_left left_op(i-1,j-1) what_right;              \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_VZ:                                                      \
        {                                                               \
            L_DOUBLE tempb = 0.0;                                       \
                                                                        \
            what_left tempb what_right;                                 \
                                                                        \
            THROW_ASSERTB(10,tempb == 0.0);                             \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_VS:                                                      \
        {                                                               \
            L_DOUBLE tempb = 0.0;                                       \
                                                                        \
            what_left tempb what_right;                                 \
                                                                        \
            THROW_ASSERTB(11,tempb == 0.0);                             \
                                                                        \
            what_left left_op(-1,-1) what_right;                        \
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

#define UNARY_MATRIX_OPERATION_VARIA_RES(what_left,what_right)          \
{                                                                       \
    long i,j;                                                           \
    long size_row;                                                      \
    long size_col;                                                      \
    long offset;                                                        \
                                                                        \
    left_op.fix_symm();                                                 \
                                                                        \
    result = left_op;                                                   \
                                                                        \
    if ( left_op.get_e_type() == E_SYMM )                               \
    {                                                                   \
        result.make_symmetric(DO_FORCE);                                \
    }                                                                   \
                                                                        \
    if ( ( left_op.get_e_type() == E_EMPTY          ) &&                \
         ( left_op.get_m_type() == SYMMETRIC_MATRIX )    )              \
    {                                                                   \
        result.make_symmetric(DO_FORCE);                                \
    }                                                                   \
                                                                        \
    offset = left_op.get_e_offs();                                      \
                                                                        \
    size_row = left_op.get_effective_height();                          \
    size_col = left_op.get_effective_width();                           \
                                                                        \
    switch ( left_op.get_e_type() )                                     \
    {                                                                   \
        case E_EMPTY:                                                   \
        {                                                               \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_ZERO:                                                    \
        {                                                               \
            L_DOUBLE tempa;                                             \
            L_DOUBLE tempb = 0.0;                                       \
                                                                        \
            tempa = what_left tempb what_right;                         \
                                                                        \
            if ( ( tempa != 0.0 ) && ( tempb != 0.0 ) )                 \
            {                                                           \
                for ( i = 1 ; i <= size_row ; i++ )                     \
                {                                                       \
                    for ( j = 1 ; j <= size_col ; j++ )                 \
                    {                                                   \
                        result(i-1,j-1)  = tempa;                       \
                        left_op(i-1,j-1) = tempb;                       \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            else                                                        \
            {                                                           \
                if ( tempa != 0.0 )                                     \
                {                                                       \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            result(i-1,j-1) = tempa;                    \
                        }                                               \
                    }                                                   \
                }                                                       \
                                                                        \
                else                                                    \
                {                                                       \
                    if ( tempb != 0.0 )                                 \
                    {                                                   \
                        for ( i = 1 ; i <= size_row ; i++ )             \
                        {                                               \
                            for ( j = 1 ; j <= size_col ; j++ )         \
                            {                                           \
                                left_op(i-1,j-1) = tempb;               \
                            }                                           \
                        }                                               \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_DIAG:                                                    \
        {                                                               \
            L_DOUBLE tempa;                                             \
            L_DOUBLE tempb = 0.0;                                       \
                                                                        \
            tempa = what_left tempb what_right;                         \
                                                                        \
            if ( ( tempa != 0.0 ) && ( tempb != 0.0 ) )                 \
            {                                                           \
                for ( i = 1 ; i <= size_row ; i++ )                     \
                {                                                       \
                    for ( j = 1 ; j <= size_col ; j++ )                 \
                    {                                                   \
                        if ( i-offset == j )                            \
                        {                                               \
                            result(i-1,j-1) = what_left left_op(i-1,j-1) what_right; \
                        }                                               \
                                                                        \
                        else                                            \
                        {                                               \
                            result(i-1,j-1)  = tempa;                   \
                            left_op(i-1,j-1) = tempb;                   \
                        }                                               \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            else                                                        \
            {                                                           \
                if ( tempa != 0.0 )                                     \
                {                                                       \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-offset == j )                        \
                            {                                           \
                                result(i-1,j-1) = what_left left_op(i-1,j-1) what_right; \
                            }                                           \
                                                                        \
                            else                                        \
                            {                                           \
                                result(i-1,j-1) = tempa;                \
                            }                                           \
                        }                                               \
                    }                                                   \
                }                                                       \
                                                                        \
                else                                                    \
                {                                                       \
                    if ( tempb != 0.0 )                                 \
                    {                                                   \
                        for ( i = 1 ; i <= size_row ; i++ )             \
                        {                                               \
                            for ( j = 1 ; j <= size_col ; j++ )         \
                            {                                           \
                                if ( i-offset == j )                    \
                                {                                       \
                                    result(i-1,j-1) = what_left left_op(i-1,j-1) what_right; \
                                }                                       \
                                                                        \
                                else                                    \
                                {                                       \
                                    left_op(i-1,j-1) = tempb;           \
                                }                                       \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    else                                                \
                    {                                                   \
                        for ( i = 1 ; i <= size_row ; i++ )             \
                        {                                               \
                            for ( j = 1 ; j <= size_col ; j++ )         \
                            {                                           \
                                if ( i-offset == j )                    \
                                {                                       \
                                    result(i-1,j-1) = what_left left_op(i-1,j-1) what_right; \
                                }                                       \
                            }                                           \
                        }                                               \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_LOWER:                                                   \
        {                                                               \
            L_DOUBLE tempa;                                             \
            L_DOUBLE tempb = 0.0;                                       \
                                                                        \
            tempa = what_left tempb what_right;                         \
                                                                        \
            if ( ( tempa != 0.0 ) && ( tempb != 0.0 ) )                 \
            {                                                           \
                for ( i = 1 ; i <= size_row ; i++ )                     \
                {                                                       \
                    for ( j = 1 ; j <= size_col ; j++ )                 \
                    {                                                   \
                        if ( i-offset >= j )                            \
                        {                                               \
                            result(i-1,j-1) = what_left left_op(i-1,j-1) what_right; \
                        }                                               \
                                                                        \
                        else                                            \
                        {                                               \
                            result(i-1,j-1)  = tempa;                   \
                            left_op(i-1,j-1) = tempb;                   \
                        }                                               \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            else                                                        \
            {                                                           \
                if ( tempa != 0.0 )                                     \
                {                                                       \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-offset >= j )                        \
                            {                                           \
                                result(i-1,j-1) = what_left left_op(i-1,j-1) what_right; \
                            }                                           \
                                                                        \
                            else                                        \
                            {                                           \
                                result(i-1,j-1) = tempa;                \
                            }                                           \
                        }                                               \
                    }                                                   \
                }                                                       \
                                                                        \
                else                                                    \
                {                                                       \
                    if ( tempb != 0.0 )                                 \
                    {                                                   \
                        for ( i = 1 ; i <= size_row ; i++ )             \
                        {                                               \
                            for ( j = 1 ; j <= size_col ; j++ )         \
                            {                                           \
                                if ( i-offset >= j )                    \
                                {                                       \
                                    result(i-1,j-1) = what_left left_op(i-1,j-1) what_right; \
                                }                                       \
                                                                        \
                                else                                    \
                                {                                       \
                                    left_op(i-1,j-1) = tempb;           \
                                }                                       \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    else                                                \
                    {                                                   \
                        for ( i = 1 ; i <= size_row ; i++ )             \
                        {                                               \
                            for ( j = 1 ; j <= size_col ; j++ )         \
                            {                                           \
                                if ( i-offset >= j )                    \
                                {                                       \
                                    result(i-1,j-1) = what_left left_op(i-1,j-1) what_right; \
                                }                                       \
                            }                                           \
                        }                                               \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_UPPER:                                                   \
        {                                                               \
            L_DOUBLE tempa;                                             \
            L_DOUBLE tempb = 0.0;                                       \
                                                                        \
            tempa = what_left tempb what_right;                         \
                                                                        \
            if ( ( tempa != 0.0 ) && ( tempb != 0.0 ) )                 \
            {                                                           \
                for ( i = 1 ; i <= size_row ; i++ )                     \
                {                                                       \
                    for ( j = 1 ; j <= size_col ; j++ )                 \
                    {                                                   \
                        if ( i-offset <= j )                            \
                        {                                               \
                            result(i-1,j-1) = what_left left_op(i-1,j-1) what_right; \
                        }                                               \
                                                                        \
                        else                                            \
                        {                                               \
                            result(i-1,j-1)  = tempa;                   \
                            left_op(i-1,j-1) = tempb;                   \
                        }                                               \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            else                                                        \
            {                                                           \
                if ( tempa != 0.0 )                                     \
                {                                                       \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-offset <= j )                        \
                            {                                           \
                                result(i-1,j-1) = what_left left_op(i-1,j-1) what_right; \
                            }                                           \
                                                                        \
                            else                                        \
                            {                                           \
                                result(i-1,j-1) = tempa;                \
                            }                                           \
                        }                                               \
                    }                                                   \
                }                                                       \
                                                                        \
                else                                                    \
                {                                                       \
                    if ( tempb != 0.0 )                                 \
                    {                                                   \
                        for ( i = 1 ; i <= size_row ; i++ )             \
                        {                                               \
                            for ( j = 1 ; j <= size_col ; j++ )         \
                            {                                           \
                                if ( i-offset <= j )                    \
                                {                                       \
                                    result(i-1,j-1) = what_left left_op(i-1,j-1) what_right; \
                                }                                       \
                                                                        \
                                else                                    \
                                {                                       \
                                    left_op(i-1,j-1) = tempb;           \
                                }                                       \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    else                                                \
                    {                                                   \
                        for ( i = 1 ; i <= size_row ; i++ )             \
                        {                                               \
                            for ( j = 1 ; j <= size_col ; j++ )         \
                            {                                           \
                                if ( i-offset <= j )                    \
                                {                                       \
                                    result(i-1,j-1) = what_left left_op(i-1,j-1) what_right; \
                                }                                       \
                            }                                           \
                        }                                               \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_SYMM:                                                    \
        {                                                               \
            for ( i = 1 ; i <= size_row ; i++ )                         \
            {                                                           \
                for ( j = 1 ; j <= i ; j++ )                            \
                {                                                       \
                   result(i-1,j-1) = what_left left_op(i-1,j-1) what_right; \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_ASYMM:                                                   \
        {                                                               \
            for ( i = 1 ; i <= size_row ; i++ )                         \
            {                                                           \
                for ( j = 1 ; j <= size_col ; j++ )                     \
                {                                                       \
                   result(i-1,j-1) = what_left left_op(i-1,j-1) what_right; \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_VZ:                                                      \
        {                                                               \
            L_DOUBLE tempa;                                             \
            L_DOUBLE tempb = 0.0;                                       \
                                                                        \
            tempa = what_left tempb what_right;                         \
                                                                        \
            THROW_ASSERTB(10,( tempa == 0.0 ) || ( tempb == 0.0 ));     \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_VS:                                                      \
        {                                                               \
            L_DOUBLE tempa;                                             \
            L_DOUBLE tempb = 0.0;                                       \
                                                                        \
            tempa = what_left tempb what_right;                         \
                                                                        \
            THROW_ASSERTB(11,( tempa == 0.0 ) || ( tempb == 0.0 ));     \
                                                                        \
            result(-1,-1) = what_left left_op(-1,-1) what_right;        \
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

#define BINARY_MATRIX_OPERATION_CONST_CONST_RES(what_mid)               \
{                                                                       \
    long i,j;                                                           \
    long size_row;                                                      \
    long size_col;                                                      \
    long roffset;                                                       \
                                                                        \
    result = left_op;                                                   \
                                                                        \
    if ( left_op.get_e_type() == E_SYMM )                               \
    {                                                                   \
        result.make_symmetric(DO_FORCE);                                \
    }                                                                   \
                                                                        \
    if ( ( left_op.get_e_type() == E_EMPTY          ) &&                \
         ( left_op.get_m_type() == SYMMETRIC_MATRIX )    )              \
    {                                                                   \
        result.make_symmetric(DO_FORCE);                                \
    }                                                                   \
                                                                        \
    switch ( left_op.get_e_type() )                                     \
    {                                                                   \
        case E_EMPTY:                                                   \
        {                                                               \
            switch ( right_op.get_e_type() )                            \
            {                                                           \
                case E_EMPTY:                                           \
                case E_VZ:                                              \
                case E_VS:                                              \
                {                                                       \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_ZERO:                                            \
                case E_DIAG:                                            \
                case E_LOWER:                                           \
                case E_UPPER:                                           \
                case E_SYMM:                                            \
                case E_ASYMM:                                           \
                {                                                       \
                    L_THROW(1);                                         \
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
        case E_ZERO:                                                    \
        case E_DIAG:                                                    \
        case E_LOWER:                                                   \
        case E_UPPER:                                                   \
        case E_ASYMM:                                                   \
        {                                                               \
            switch ( right_op.get_e_type() )                            \
            {                                                           \
                case E_EMPTY:                                           \
                {                                                       \
                    L_THROW(5);                                         \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_ZERO:                                            \
                case E_VZ:                                              \
                {                                                       \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_DIAG:                                            \
                {                                                       \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-roffset == j )                       \
                            {                                           \
                                result(i-1,j-1) = (left_op.get_offset_element(i-1,j-1)) what_mid (right_op.get_offset_element(i-1,j-1)); \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_LOWER:                                           \
                {                                                       \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-roffset >= j )                       \
                            {                                           \
                                result(i-1,j-1) = (left_op.get_offset_element(i-1,j-1)) what_mid (right_op.get_offset_element(i-1,j-1)); \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_UPPER:                                           \
                {                                                       \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-roffset <= j )                       \
                            {                                           \
                                result(i-1,j-1) = (left_op.get_offset_element(i-1,j-1)) what_mid (right_op.get_offset_element(i-1,j-1)); \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_SYMM:                                            \
                case E_ASYMM:                                           \
                {                                                       \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            result(i-1,j-1) = (left_op.get_offset_element(i-1,j-1)) what_mid (right_op.get_offset_element(i-1,j-1)); \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_VS:                                              \
                {                                                       \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i == j )                               \
                            {                                           \
                                result(i-1,j-1) = (left_op.get_offset_element(i-1,j-1)) what_mid (right_op.get_offset_element(i-1,j-1)); \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(6);                                         \
                                                                        \
                    break;                                              \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_SYMM:                                                    \
        {                                                               \
            switch ( right_op.get_e_type() )                            \
            {                                                           \
                case E_EMPTY:                                           \
                {                                                       \
                    L_THROW(11);                                        \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_ZERO:                                            \
                case E_VZ:                                              \
                {                                                       \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_SYMM:                                            \
                {                                                       \
                    size_row = left_op.get_effective_height();          \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= i ; j++ )                    \
                        {                                               \
                            result(i-1,j-1) = (left_op.get_offset_element(i-1,j-1)) what_mid (right_op.get_offset_element(i-1,j-1)); \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_VS:                                              \
                {                                                       \
                    size_row = left_op.get_effective_height();          \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= i ; j++ )                    \
                        {                                               \
                            if ( i == j )                               \
                            {                                           \
                                result(i-1,j-1) = (left_op.get_offset_element(i-1,j-1)) what_mid (right_op.get_offset_element(i-1,j-1)); \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_DIAG:                                            \
                {                                                       \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    if ( roffset != 0 )                                 \
                    {                                                   \
                        result.make_asymmetric(DO_FORCE);               \
                    }                                                   \
                                                                        \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-roffset == j )                       \
                            {                                           \
                                result(i-1,j-1) = (left_op.get_offset_element(i-1,j-1)) what_mid (right_op.get_offset_element(i-1,j-1)); \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_LOWER:                                           \
                {                                                       \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    result.make_asymmetric(DO_FORCE);                   \
                                                                        \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-roffset >= j )                       \
                            {                                           \
                                result(i-1,j-1) = (left_op.get_offset_element(i-1,j-1)) what_mid (right_op.get_offset_element(i-1,j-1)); \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_UPPER:                                           \
                {                                                       \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    result.make_asymmetric(DO_FORCE);                   \
                                                                        \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-roffset <= j )                       \
                            {                                           \
                                result(i-1,j-1) = (left_op.get_offset_element(i-1,j-1)) what_mid (right_op.get_offset_element(i-1,j-1)); \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_ASYMM:                                           \
                {                                                       \
                    result.make_asymmetric(DO_FORCE);                   \
                                                                        \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            result(i-1,j-1) = (left_op.get_offset_element(i-1,j-1)) what_mid (right_op.get_offset_element(i-1,j-1)); \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(12);                                        \
                                                                        \
                    break;                                              \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_VZ:                                                      \
        {                                                               \
            switch ( right_op.get_e_type() )                            \
            {                                                           \
                case E_EMPTY:                                           \
                {                                                       \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_ZERO:                                            \
                {                                                       \
                    result.make_diagonal(DO_FORCE);                     \
                                                                        \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    if ( size_row > size_col )                          \
                    {                                                   \
                        result.pad_matrix(size_row);                    \
                    }                                                   \
                                                                        \
                    else                                                \
                    {                                                   \
                        result.pad_matrix(size_col);                    \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_DIAG:                                            \
                {                                                       \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    result.make_diagonal(DO_FORCE);                     \
                                                                        \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    if ( size_row > size_col )                          \
                    {                                                   \
                        result.pad_matrix(size_row);                    \
                    }                                                   \
                                                                        \
                    else                                                \
                    {                                                   \
                        result.pad_matrix(size_col);                    \
                    }                                                   \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-roffset == j )                       \
                            {                                           \
                                result(i-1,j-1) = (left_op.get_offset_element(i-1,j-1)) what_mid (right_op.get_offset_element(i-1,j-1)); \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_LOWER:                                           \
                {                                                       \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    result.make_diagonal(DO_FORCE);                     \
                                                                        \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    if ( size_row > size_col )                          \
                    {                                                   \
                        result.pad_matrix(size_row);                    \
                    }                                                   \
                                                                        \
                    else                                                \
                    {                                                   \
                        result.pad_matrix(size_col);                    \
                    }                                                   \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-roffset >= j )                       \
                            {                                           \
                                result(i-1,j-1) = (left_op.get_offset_element(i-1,j-1)) what_mid (right_op.get_offset_element(i-1,j-1)); \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_UPPER:                                           \
                {                                                       \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    result.make_diagonal(DO_FORCE);                     \
                                                                        \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    if ( size_row > size_col )                          \
                    {                                                   \
                        result.pad_matrix(size_row);                    \
                    }                                                   \
                                                                        \
                    else                                                \
                    {                                                   \
                        result.pad_matrix(size_col);                    \
                    }                                                   \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-roffset <= j )                       \
                            {                                           \
                                result(i-1,j-1) = (left_op.get_offset_element(i-1,j-1)) what_mid (right_op.get_offset_element(i-1,j-1)); \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_SYMM:                                            \
                {                                                       \
                    result.make_symmetric(DO_FORCE);                    \
                                                                        \
                    size_row = right_op.get_effective_height();         \
                                                                        \
                    result.pad_matrix(size_row);                        \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= i ; j++ )                    \
                        {                                               \
                            result(i-1,j-1) = (left_op.get_offset_element(i-1,j-1)) what_mid (right_op.get_offset_element(i-1,j-1)); \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_ASYMM:                                           \
                {                                                       \
                    result.make_diagonal(DO_FORCE);                     \
                                                                        \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    if ( size_row > size_col )                          \
                    {                                                   \
                        result.pad_matrix(size_row);                    \
                    }                                                   \
                                                                        \
                    else                                                \
                    {                                                   \
                        result.pad_matrix(size_col);                    \
                    }                                                   \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            result(i-1,j-1) = (left_op.get_offset_element(i-1,j-1)) what_mid (right_op.get_offset_element(i-1,j-1)); \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_VZ:                                              \
                {                                                       \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_VS:                                              \
                {                                                       \
                    result.make_scalar(DO_FORCE);                       \
                                                                        \
                    result(-1,-1) = (left_op.get_offset_element(0,0)) what_mid (right_op.get_offset_element(-1,-1)); \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(17);                                        \
                                                                        \
                    break;                                              \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_VS:                                                      \
        {                                                               \
            switch ( right_op.get_e_type() )                            \
            {                                                           \
                case E_EMPTY:                                           \
                {                                                       \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_ZERO:                                            \
                {                                                       \
                    L_DOUBLE tempa;                                     \
                                                                        \
                    tempa = left_op.get_offset_element(-1,-1);          \
                                                                        \
                    result.make_diagonal(DO_FORCE);                     \
                                                                        \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    if ( size_row > size_col )                          \
                    {                                                   \
                        result.pad_matrix(size_row);                    \
                                                                        \
                        for ( i = 1 ; i <= size_row ; i++ )             \
                        {                                               \
                            result(i-1,i-1) = tempa;                    \
                        }                                               \
                    }                                                   \
                                                                        \
                    else                                                \
                    {                                                   \
                        result.pad_matrix(size_col);                    \
                                                                        \
                        for ( i = 1 ; i <= size_col ; i++ )             \
                        {                                               \
                            result(i-1,i-1) = tempa;                    \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_DIAG:                                            \
                {                                                       \
                    L_DOUBLE tempa;                                     \
                                                                        \
                    tempa = left_op.get_offset_element(-1,-1);          \
                                                                        \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    result.make_diagonal(DO_FORCE);                     \
                                                                        \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    if ( size_row > size_col )                          \
                    {                                                   \
                        result.pad_matrix(size_row);                    \
                                                                        \
                        for ( i = 1 ; i <= size_row ; i++ )             \
                        {                                               \
                            result(i-1,i-1) = tempa;                    \
                        }                                               \
                    }                                                   \
                                                                        \
                    else                                                \
                    {                                                   \
                        result.pad_matrix(size_col);                    \
                                                                        \
                        for ( i = 1 ; i <= size_col ; i++ )             \
                        {                                               \
                            result(i-1,i-1) = tempa;                    \
                        }                                               \
                    }                                                   \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-roffset == j )                       \
                            {                                           \
                                result(i-1,j-1) = (left_op.get_offset_element(i-1,j-1)) what_mid (right_op.get_offset_element(i-1,j-1)); \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_LOWER:                                           \
                {                                                       \
                    L_DOUBLE tempa;                                     \
                                                                        \
                    tempa = left_op.get_offset_element(-1,-1);          \
                                                                        \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    result.make_diagonal(DO_FORCE);                     \
                                                                        \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    if ( size_row > size_col )                          \
                    {                                                   \
                        result.pad_matrix(size_row);                    \
                                                                        \
                        for ( i = 1 ; i <= size_row ; i++ )             \
                        {                                               \
                            result(i-1,i-1) = tempa;                    \
                        }                                               \
                    }                                                   \
                                                                        \
                    else                                                \
                    {                                                   \
                        result.pad_matrix(size_col);                    \
                                                                        \
                        for ( i = 1 ; i <= size_col ; i++ )             \
                        {                                               \
                            result(i-1,i-1) = tempa;                    \
                        }                                               \
                    }                                                   \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-roffset >= j )                       \
                            {                                           \
                                result(i-1,j-1) = (left_op.get_offset_element(i-1,j-1)) what_mid (right_op.get_offset_element(i-1,j-1)); \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_UPPER:                                           \
                {                                                       \
                    L_DOUBLE tempa;                                     \
                                                                        \
                    tempa = left_op.get_offset_element(-1,-1);          \
                                                                        \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    result.make_diagonal(DO_FORCE);                     \
                                                                        \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    if ( size_row > size_col )                          \
                    {                                                   \
                        result.pad_matrix(size_row);                    \
                                                                        \
                        for ( i = 1 ; i <= size_row ; i++ )             \
                        {                                               \
                            result(i-1,i-1) = tempa;                    \
                        }                                               \
                    }                                                   \
                                                                        \
                    else                                                \
                    {                                                   \
                        result.pad_matrix(size_col);                    \
                                                                        \
                        for ( i = 1 ; i <= size_col ; i++ )             \
                        {                                               \
                            result(i-1,i-1) = tempa;                    \
                        }                                               \
                    }                                                   \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-roffset <= j )                       \
                            {                                           \
                                result(i-1,j-1) = (left_op.get_offset_element(i-1,j-1)) what_mid (right_op.get_offset_element(i-1,j-1)); \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_SYMM:                                            \
                {                                                       \
                    L_DOUBLE tempa;                                     \
                                                                        \
                    tempa = left_op.get_offset_element(-1,-1);          \
                                                                        \
                    result.make_symmetric(DO_FORCE);                    \
                                                                        \
                    size_row = right_op.get_effective_height();         \
                                                                        \
                    result.pad_matrix(size_row);                        \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        result(i-1,i-1) = tempa;                        \
                    }                                                   \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= i ; j++ )                    \
                        {                                               \
                            result(i-1,j-1) = (left_op.get_offset_element(i-1,j-1)) what_mid (right_op.get_offset_element(i-1,j-1)); \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_ASYMM:                                           \
                {                                                       \
                    L_DOUBLE tempa;                                     \
                                                                        \
                    tempa = left_op.get_offset_element(-1,-1);          \
                                                                        \
                    result.make_diagonal(DO_FORCE);                     \
                                                                        \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    if ( size_row > size_col )                          \
                    {                                                   \
                        result.pad_matrix(size_row);                    \
                                                                        \
                        for ( i = 1 ; i <= size_row ; i++ )             \
                        {                                               \
                            result(i-1,i-1) = tempa;                    \
                        }                                               \
                    }                                                   \
                                                                        \
                    else                                                \
                    {                                                   \
                        result.pad_matrix(size_col);                    \
                                                                        \
                        for ( i = 1 ; i <= size_col ; i++ )             \
                        {                                               \
                            result(i-1,i-1) = tempa;                    \
                        }                                               \
                    }                                                   \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            result(i-1,j-1) = (left_op.get_offset_element(i-1,j-1)) what_mid (right_op.get_offset_element(i-1,j-1)); \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_VZ:                                              \
                {                                                       \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_VS:                                              \
                {                                                       \
                    result(-1,-1) = (left_op.get_offset_element(-1,-1)) what_mid (right_op.get_offset_element(-1,-1)); \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(17);                                        \
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
            L_THROW(21);                                                \
                                                                        \
            break;                                                      \
        }                                                               \
    }                                                                   \
}

#define BINARY_MATRIX_OPERATION_VARIA_CONST(what_mid)                   \
{                                                                       \
    long i,j;                                                           \
    long size_row;                                                      \
    long size_col;                                                      \
    long roffset;                                                       \
                                                                        \
    switch ( left_op.get_e_type() )                                     \
    {                                                                   \
        case E_EMPTY:                                                   \
        {                                                               \
            switch ( right_op.get_e_type() )                            \
            {                                                           \
                case E_EMPTY:                                           \
                case E_VZ:                                              \
                case E_VS:                                              \
                {                                                       \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_ZERO:                                            \
                case E_DIAG:                                            \
                case E_LOWER:                                           \
                case E_UPPER:                                           \
                case E_SYMM:                                            \
                case E_ASYMM:                                           \
                {                                                       \
                    L_THROW(1);                                         \
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
        case E_ZERO:                                                    \
        case E_DIAG:                                                    \
        case E_LOWER:                                                   \
        case E_UPPER:                                                   \
        case E_ASYMM:                                                   \
        {                                                               \
            switch ( right_op.get_e_type() )                            \
            {                                                           \
                case E_EMPTY:                                           \
                {                                                       \
                    L_THROW(5);                                         \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_ZERO:                                            \
                case E_VZ:                                              \
                {                                                       \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_DIAG:                                            \
                {                                                       \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-roffset == j )                       \
                            {                                           \
                                left_op(i-1,j-1) what_mid (right_op.get_offset_element(i-1,j-1)); \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_LOWER:                                           \
                {                                                       \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-roffset >= j )                       \
                            {                                           \
                                left_op(i-1,j-1) what_mid (right_op.get_offset_element(i-1,j-1)); \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_UPPER:                                           \
                {                                                       \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-roffset <= j )                       \
                            {                                           \
                                left_op(i-1,j-1) what_mid (right_op.get_offset_element(i-1,j-1)); \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_SYMM:                                            \
                case E_ASYMM:                                           \
                {                                                       \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            left_op(i-1,j-1) what_mid (right_op.get_offset_element(i-1,j-1)); \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_VS:                                              \
                {                                                       \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i == j )                               \
                            {                                           \
                                left_op(i-1,j-1) what_mid (right_op.get_offset_element(i-1,j-1)); \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(6);                                         \
                                                                        \
                    break;                                              \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_SYMM:                                                    \
        {                                                               \
            switch ( right_op.get_e_type() )                            \
            {                                                           \
                case E_EMPTY:                                           \
                {                                                       \
                    L_THROW(11);                                        \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_ZERO:                                            \
                case E_VZ:                                              \
                {                                                       \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_SYMM:                                            \
                {                                                       \
                    size_row = left_op.get_effective_height();          \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= i ; j++ )                    \
                        {                                               \
                            left_op(i-1,j-1) what_mid (right_op.get_offset_element(i-1,j-1)); \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_VS:                                              \
                {                                                       \
                    size_row = left_op.get_effective_height();          \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= i ; j++ )                    \
                        {                                               \
                            if ( i == j )                               \
                            {                                           \
                                left_op(i-1,j-1) what_mid (right_op.get_offset_element(i-1,j-1)); \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_DIAG:                                            \
                {                                                       \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    if ( roffset != 0 )                                 \
                    {                                                   \
                        left_op.make_asymmetric(DO_FORCE);              \
                    }                                                   \
                                                                        \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-roffset == j )                       \
                            {                                           \
                                left_op(i-1,j-1) what_mid (right_op.get_offset_element(i-1,j-1)); \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_LOWER:                                           \
                {                                                       \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    left_op.make_asymmetric(DO_FORCE);                  \
                                                                        \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-roffset >= j )                       \
                            {                                           \
                                left_op(i-1,j-1) what_mid (right_op.get_offset_element(i-1,j-1)); \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_UPPER:                                           \
                {                                                       \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    left_op.make_asymmetric(DO_FORCE);                  \
                                                                        \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-roffset <= j )                       \
                            {                                           \
                                left_op(i-1,j-1) what_mid (right_op.get_offset_element(i-1,j-1)); \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_ASYMM:                                           \
                {                                                       \
                    left_op.make_asymmetric(DO_FORCE);                  \
                                                                        \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            left_op(i-1,j-1) what_mid (right_op.get_offset_element(i-1,j-1)); \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(12);                                        \
                                                                        \
                    break;                                              \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_VZ:                                                      \
        {                                                               \
            switch ( right_op.get_e_type() )                            \
            {                                                           \
                case E_EMPTY:                                           \
                {                                                       \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_ZERO:                                            \
                {                                                       \
                    left_op.make_diagonal(DO_FORCE);                    \
                                                                        \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    if ( size_row > size_col )                          \
                    {                                                   \
                        left_op.pad_matrix(size_row);                   \
                    }                                                   \
                                                                        \
                    else                                                \
                    {                                                   \
                        left_op.pad_matrix(size_col);                   \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_DIAG:                                            \
                {                                                       \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    left_op.make_diagonal(DO_FORCE);                    \
                                                                        \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    if ( size_row > size_col )                          \
                    {                                                   \
                        left_op.pad_matrix(size_row);                   \
                    }                                                   \
                                                                        \
                    else                                                \
                    {                                                   \
                        left_op.pad_matrix(size_col);                   \
                    }                                                   \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-roffset == j )                       \
                            {                                           \
                                left_op(i-1,j-1) what_mid (right_op.get_offset_element(i-1,j-1)); \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_LOWER:                                           \
                {                                                       \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    left_op.make_diagonal(DO_FORCE);                    \
                                                                        \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    if ( size_row > size_col )                          \
                    {                                                   \
                        left_op.pad_matrix(size_row);                   \
                    }                                                   \
                                                                        \
                    else                                                \
                    {                                                   \
                        left_op.pad_matrix(size_col);                   \
                    }                                                   \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-roffset >= j )                       \
                            {                                           \
                                left_op(i-1,j-1) what_mid (right_op.get_offset_element(i-1,j-1)); \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_UPPER:                                           \
                {                                                       \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    left_op.make_diagonal(DO_FORCE);                    \
                                                                        \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    if ( size_row > size_col )                          \
                    {                                                   \
                        left_op.pad_matrix(size_row);                   \
                    }                                                   \
                                                                        \
                    else                                                \
                    {                                                   \
                        left_op.pad_matrix(size_col);                   \
                    }                                                   \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-roffset <= j )                       \
                            {                                           \
                                left_op(i-1,j-1) what_mid (right_op.get_offset_element(i-1,j-1)); \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_SYMM:                                            \
                {                                                       \
                    left_op.make_symmetric(DO_FORCE);                   \
                                                                        \
                    size_row = right_op.get_effective_height();         \
                                                                        \
                    left_op.pad_matrix(size_row);                       \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= i ; j++ )                    \
                        {                                               \
                            left_op(i-1,j-1) what_mid (right_op.get_offset_element(i-1,j-1)); \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_ASYMM:                                           \
                {                                                       \
                    left_op.make_diagonal(DO_FORCE);                    \
                                                                        \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    if ( size_row > size_col )                          \
                    {                                                   \
                        left_op.pad_matrix(size_row);                   \
                    }                                                   \
                                                                        \
                    else                                                \
                    {                                                   \
                        left_op.pad_matrix(size_col);                   \
                    }                                                   \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            left_op(i-1,j-1) what_mid (right_op.get_offset_element(i-1,j-1)); \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_VZ:                                              \
                {                                                       \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_VS:                                              \
                {                                                       \
                    left_op.make_scalar(DO_FORCE);                      \
                                                                        \
                    left_op(-1,-1) what_mid (right_op.get_offset_element(-1,-1)); \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(17);                                        \
                                                                        \
                    break;                                              \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_VS:                                                      \
        {                                                               \
            switch ( right_op.get_e_type() )                            \
            {                                                           \
                case E_EMPTY:                                           \
                {                                                       \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_ZERO:                                            \
                {                                                       \
                    L_DOUBLE tempa;                                     \
                                                                        \
                    tempa = left_op(-1,-1);                             \
                                                                        \
                    left_op.make_diagonal(DO_FORCE);                    \
                                                                        \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    if ( size_row > size_col )                          \
                    {                                                   \
                        left_op.pad_matrix(size_row);                   \
                                                                        \
                        for ( i = 1 ; i <= size_row ; i++ )             \
                        {                                               \
                            left_op(i-1,i-1) = tempa;                   \
                        }                                               \
                    }                                                   \
                                                                        \
                    else                                                \
                    {                                                   \
                        left_op.pad_matrix(size_col);                   \
                                                                        \
                        for ( i = 1 ; i <= size_col ; i++ )             \
                        {                                               \
                            left_op(i-1,i-1) = tempa;                   \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_DIAG:                                            \
                {                                                       \
                    L_DOUBLE tempa;                                     \
                                                                        \
                    tempa = left_op(-1,-1);                             \
                                                                        \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    left_op.make_diagonal(DO_FORCE);                    \
                                                                        \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    if ( size_row > size_col )                          \
                    {                                                   \
                        left_op.pad_matrix(size_row);                   \
                                                                        \
                        for ( i = 1 ; i <= size_row ; i++ )             \
                        {                                               \
                            left_op(i-1,i-1) = tempa;                   \
                        }                                               \
                    }                                                   \
                                                                        \
                    else                                                \
                    {                                                   \
                        left_op.pad_matrix(size_col);                   \
                                                                        \
                        for ( i = 1 ; i <= size_col ; i++ )             \
                        {                                               \
                            left_op(i-1,i-1) = tempa;                   \
                        }                                               \
                    }                                                   \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-roffset == j )                       \
                            {                                           \
                                left_op(i-1,j-1) what_mid (right_op.get_offset_element(i-1,j-1)); \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_LOWER:                                           \
                {                                                       \
                    L_DOUBLE tempa;                                     \
                                                                        \
                    tempa = left_op(-1,-1);                             \
                                                                        \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    left_op.make_diagonal(DO_FORCE);                    \
                                                                        \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    if ( size_row > size_col )                          \
                    {                                                   \
                        left_op.pad_matrix(size_row);                   \
                                                                        \
                        for ( i = 1 ; i <= size_row ; i++ )             \
                        {                                               \
                            left_op(i-1,i-1) = tempa;                   \
                        }                                               \
                    }                                                   \
                                                                        \
                    else                                                \
                    {                                                   \
                        left_op.pad_matrix(size_col);                   \
                                                                        \
                        for ( i = 1 ; i <= size_col ; i++ )             \
                        {                                               \
                            left_op(i-1,i-1) = tempa;                   \
                        }                                               \
                    }                                                   \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-roffset >= j )                       \
                            {                                           \
                                left_op(i-1,j-1) what_mid (right_op.get_offset_element(i-1,j-1)); \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_UPPER:                                           \
                {                                                       \
                    L_DOUBLE tempa;                                     \
                                                                        \
                    tempa = left_op(-1,-1);                             \
                                                                        \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    left_op.make_diagonal(DO_FORCE);                    \
                                                                        \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    if ( size_row > size_col )                          \
                    {                                                   \
                        left_op.pad_matrix(size_row);                   \
                                                                        \
                        for ( i = 1 ; i <= size_row ; i++ )             \
                        {                                               \
                            left_op(i-1,i-1) = tempa;                   \
                        }                                               \
                    }                                                   \
                                                                        \
                    else                                                \
                    {                                                   \
                        left_op.pad_matrix(size_col);                   \
                                                                        \
                        for ( i = 1 ; i <= size_col ; i++ )             \
                        {                                               \
                            left_op(i-1,i-1) = tempa;                   \
                        }                                               \
                    }                                                   \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-roffset <= j )                       \
                            {                                           \
                                left_op(i-1,j-1) what_mid (right_op.get_offset_element(i-1,j-1)); \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_SYMM:                                            \
                {                                                       \
                    L_DOUBLE tempa;                                     \
                                                                        \
                    tempa = left_op(-1,-1);                             \
                                                                        \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    left_op.make_diagonal(DO_FORCE);                    \
                                                                        \
                    size_row = right_op.get_effective_height();         \
                                                                        \
                    left_op.pad_matrix(size_row);                       \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        left_op(i-1,i-1) = tempa;                       \
                    }                                                   \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_row ; j++ )             \
                        {                                               \
                            left_op(i-1,j-1) what_mid (right_op.get_offset_element(i-1,j-1)); \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_ASYMM:                                           \
                {                                                       \
                    L_DOUBLE tempa;                                     \
                                                                        \
                    tempa = left_op(-1,-1);                             \
                                                                        \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    left_op.make_diagonal(DO_FORCE);                    \
                                                                        \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    if ( size_row > size_col )                          \
                    {                                                   \
                        left_op.pad_matrix(size_row);                   \
                                                                        \
                        for ( i = 1 ; i <= size_row ; i++ )             \
                        {                                               \
                            left_op(i-1,i-1) = tempa;                   \
                        }                                               \
                    }                                                   \
                                                                        \
                    else                                                \
                    {                                                   \
                        left_op.pad_matrix(size_col);                   \
                                                                        \
                        for ( i = 1 ; i <= size_col ; i++ )             \
                        {                                               \
                            left_op(i-1,i-1) = tempa;                   \
                        }                                               \
                    }                                                   \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            left_op(i-1,j-1) what_mid (right_op.get_offset_element(i-1,j-1)); \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_VZ:                                              \
                {                                                       \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_VS:                                              \
                {                                                       \
                    left_op(-1,-1) what_mid (right_op.get_offset_element(-1,-1)); \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(17);                                        \
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
            L_THROW(21);                                                \
                                                                        \
            break;                                                      \
        }                                                               \
    }                                                                   \
}




//
// Generic scalar, vector and matrix relational operators
//
// what_a: either || or && - this defines how individual relationships
//         are concatenated to form complete result
// what_b: actual relation between individual elements
// what_c: what result starts as - generally, 0 if what_a is || or 1 if
//         what_a is &&
//

#define MATRIX_MATRIX_RELATION(what_a,what_b,what_c)                    \
{                                                                       \
    long i,j;                                                           \
    long size_row;                                                      \
    long size_col;                                                      \
    long loffset;                                                       \
    long roffset;                                                       \
                                                                        \
    result = what_c;                                                    \
                                                                        \
    switch ( left_op.get_e_type() )                                     \
    {                                                                   \
        case E_EMPTY:                                                   \
        {                                                               \
            switch ( right_op.get_e_type() )                            \
            {                                                           \
                case E_EMPTY:                                           \
                case E_VZ:                                              \
                case E_VS:                                              \
                {                                                       \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_ZERO:                                            \
                case E_DIAG:                                            \
                case E_LOWER:                                           \
                case E_UPPER:                                           \
                case E_SYMM:                                            \
                case E_ASYMM:                                           \
                {                                                       \
                    L_THROW(1);                                         \
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
        case E_ZERO:                                                    \
        {                                                               \
            switch ( right_op.get_e_type() )                            \
            {                                                           \
                case E_EMPTY:                                           \
                {                                                       \
                    L_THROW(3);                                         \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_ZERO:                                            \
                {                                                       \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_DIAG:                                            \
                {                                                       \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-roffset == j )                       \
                            {                                           \
                                result =                                \
                                ( result what_a                         \
                                ( left_op.get_offset_element(i-1,j-1) what_b \
                                  right_op.get_offset_element(i-1,j-1) ) );  \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_LOWER:                                           \
                {                                                       \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-roffset >= j )                       \
                            {                                           \
                                result =                                \
                                ( result what_a                         \
                                ( left_op.get_offset_element(i-1,j-1) what_b \
                                  right_op.get_offset_element(i-1,j-1) ) );  \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_UPPER:                                           \
                {                                                       \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-roffset <= j )                       \
                            {                                           \
                                result =                                \
                                ( result what_a                         \
                                ( left_op.get_offset_element(i-1,j-1) what_b \
                                  right_op.get_offset_element(i-1,j-1) ) );  \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_SYMM:                                            \
                {                                                       \
                    size_row = right_op.get_effective_height();         \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= i ; j++ )                    \
                        {                                               \
                            result =                                    \
                            ( result what_a                             \
                            ( left_op.get_offset_element(i-1,j-1) what_b \
                              right_op.get_offset_element(i-1,j-1) ) );  \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_ASYMM:                                           \
                {                                                       \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            result =                                    \
                            ( result what_a                             \
                            ( left_op.get_offset_element(i-1,j-1) what_b \
                              right_op.get_offset_element(i-1,j-1) ) );  \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_VZ:                                              \
                {                                                       \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_VS:                                              \
                {                                                       \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                    result = ( result what_a ( 0.0 what_b right_op.get_offset_element(-1,-1) ) ); \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(4);                                         \
                                                                        \
                    break;                                              \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_DIAG:                                                    \
        {                                                               \
            switch ( right_op.get_e_type() )                            \
            {                                                           \
                case E_EMPTY:                                           \
                {                                                       \
                    L_THROW(5);                                         \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_ZERO:                                            \
                {                                                       \
                    loffset = left_op.get_e_offs();                     \
                                                                        \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-loffset == j )                       \
                            {                                           \
                                result =                                \
                                ( result what_a                         \
                                ( left_op.get_offset_element(i-1,j-1) what_b \
                                  right_op.get_offset_element(i-1,j-1) ) );  \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_DIAG:                                            \
                {                                                       \
                    loffset = left_op.get_e_offs();                     \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( ( i-loffset == j ) ||                  \
                                 ( i-roffset == j )    )                \
                            {                                           \
                                result =                                \
                                ( result what_a                         \
                                ( left_op.get_offset_element(i-1,j-1) what_b \
                                  right_op.get_offset_element(i-1,j-1) ) );  \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_LOWER:                                           \
                {                                                       \
                    loffset = left_op.get_e_offs();                     \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( ( i-loffset == j ) ||                  \
                                 ( i-roffset >= j )    )                \
                            {                                           \
                                result =                                \
                                ( result what_a                         \
                                ( left_op.get_offset_element(i-1,j-1) what_b \
                                  right_op.get_offset_element(i-1,j-1) ) );  \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_UPPER:                                           \
                {                                                       \
                    loffset = left_op.get_e_offs();                     \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( ( i-loffset == j ) ||                  \
                                 ( i-roffset <= j )    )                \
                            {                                           \
                                result =                                \
                                ( result what_a                         \
                                ( left_op.get_offset_element(i-1,j-1) what_b \
                                  right_op.get_offset_element(i-1,j-1) ) );  \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_SYMM:                                            \
                case E_ASYMM:                                           \
                {                                                       \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            result =                                    \
                            ( result what_a                             \
                            ( left_op.get_offset_element(i-1,j-1) what_b \
                              right_op.get_offset_element(i-1,j-1) ) );  \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_VZ:                                              \
                {                                                       \
                    loffset = left_op.get_e_offs();                     \
                                                                        \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-loffset == j )                       \
                            {                                           \
                                result =                                \
                                ( result what_a                         \
                                ( left_op.get_offset_element(i-1,j-1) what_b \
                                  right_op.get_offset_element(i-1,j-1) ) );  \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_VS:                                              \
                {                                                       \
                    loffset = left_op.get_e_offs();                     \
                                                                        \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( ( i-loffset == j ) ||                  \
                                 ( i         == j )    )                \
                            {                                           \
                                result =                                \
                                ( result what_a                         \
                                ( left_op.get_offset_element(i-1,j-1) what_b \
                                  right_op.get_offset_element(i-1,j-1) ) );  \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(6);                                         \
                                                                        \
                    break;                                              \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_LOWER:                                                   \
        {                                                               \
            switch ( right_op.get_e_type() )                            \
            {                                                           \
                case E_EMPTY:                                           \
                {                                                       \
                    L_THROW(7);                                         \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_ZERO:                                            \
                {                                                       \
                    loffset = left_op.get_e_offs();                     \
                                                                        \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-loffset >= j )                       \
                            {                                           \
                                result =                                \
                                ( result what_a                         \
                                ( left_op.get_offset_element(i-1,j-1) what_b \
                                  right_op.get_offset_element(i-1,j-1) ) );  \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_DIAG:                                            \
                {                                                       \
                    loffset = left_op.get_e_offs();                     \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( ( i-loffset >= j ) ||                  \
                                 ( i-roffset == j )    )                \
                            {                                           \
                                result =                                \
                                ( result what_a                         \
                                ( left_op.get_offset_element(i-1,j-1) what_b \
                                  right_op.get_offset_element(i-1,j-1) ) );  \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_LOWER:                                           \
                {                                                       \
                    loffset = left_op.get_e_offs();                     \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( ( i-loffset >= j ) ||                  \
                                 ( i-roffset >= j )    )                \
                            {                                           \
                                result =                                \
                                ( result what_a                         \
                                ( left_op.get_offset_element(i-1,j-1) what_b \
                                  right_op.get_offset_element(i-1,j-1) ) );  \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_UPPER:                                           \
                {                                                       \
                    loffset = left_op.get_e_offs();                     \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( ( i-loffset >= j ) ||                  \
                                 ( i-roffset <= j )    )                \
                            {                                           \
                                result =                                \
                                ( result what_a                         \
                                ( left_op.get_offset_element(i-1,j-1) what_b \
                                  right_op.get_offset_element(i-1,j-1) ) );  \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_SYMM:                                            \
                case E_ASYMM:                                           \
                {                                                       \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            result =                                    \
                            ( result what_a                             \
                            ( left_op.get_offset_element(i-1,j-1) what_b \
                              right_op.get_offset_element(i-1,j-1) ) );  \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_VZ:                                              \
                {                                                       \
                    loffset = left_op.get_e_offs();                     \
                                                                        \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-loffset >= j )                       \
                            {                                           \
                                result =                                \
                                ( result what_a                         \
                                ( left_op.get_offset_element(i-1,j-1) what_b \
                                  right_op.get_offset_element(i-1,j-1) ) );  \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_VS:                                              \
                {                                                       \
                    loffset = left_op.get_e_offs();                     \
                                                                        \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( ( i-loffset >= j ) ||                  \
                                 ( i         == j )    )                \
                            {                                           \
                                result =                                \
                                ( result what_a                         \
                                ( left_op.get_offset_element(i-1,j-1) what_b \
                                  right_op.get_offset_element(i-1,j-1) ) );  \
                            }                                           \
                        }                                               \
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
        case E_UPPER:                                                   \
        {                                                               \
            switch ( right_op.get_e_type() )                            \
            {                                                           \
                case E_EMPTY:                                           \
                {                                                       \
                    L_THROW(9);                                         \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_ZERO:                                            \
                {                                                       \
                    loffset = left_op.get_e_offs();                     \
                                                                        \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-loffset <= j )                       \
                            {                                           \
                                result =                                \
                                ( result what_a                         \
                                ( left_op.get_offset_element(i-1,j-1) what_b \
                                  right_op.get_offset_element(i-1,j-1) ) );  \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_DIAG:                                            \
                {                                                       \
                    loffset = left_op.get_e_offs();                     \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( ( i-loffset <= j ) ||                  \
                                 ( i-roffset == j )    )                \
                            {                                           \
                                result =                                \
                                ( result what_a                         \
                                ( left_op.get_offset_element(i-1,j-1) what_b \
                                  right_op.get_offset_element(i-1,j-1) ) );  \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_LOWER:                                           \
                {                                                       \
                    loffset = left_op.get_e_offs();                     \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( ( i-loffset <= j ) ||                  \
                                 ( i-roffset >= j )    )                \
                            {                                           \
                                result =                                \
                                ( result what_a                         \
                                ( left_op.get_offset_element(i-1,j-1) what_b \
                                  right_op.get_offset_element(i-1,j-1) ) );  \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_UPPER:                                           \
                {                                                       \
                    loffset = left_op.get_e_offs();                     \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( ( i-loffset <= j ) ||                  \
                                 ( i-roffset <= j )    )                \
                            {                                           \
                                result =                                \
                                ( result what_a                         \
                                ( left_op.get_offset_element(i-1,j-1) what_b \
                                  right_op.get_offset_element(i-1,j-1) ) );  \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_SYMM:                                            \
                case E_ASYMM:                                           \
                {                                                       \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            result =                                    \
                            ( result what_a                             \
                            ( left_op.get_offset_element(i-1,j-1) what_b \
                              right_op.get_offset_element(i-1,j-1) ) );  \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_VZ:                                              \
                {                                                       \
                    loffset = left_op.get_e_offs();                     \
                                                                        \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-loffset <= j )                       \
                            {                                           \
                                result =                                \
                                ( result what_a                         \
                                ( left_op.get_offset_element(i-1,j-1) what_b \
                                  right_op.get_offset_element(i-1,j-1) ) );  \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_VS:                                              \
                {                                                       \
                    loffset = left_op.get_e_offs();                     \
                                                                        \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( ( i-loffset <= j ) ||                  \
                                 ( i         == j )    )                \
                            {                                           \
                                result =                                \
                                ( result what_a                         \
                                ( left_op.get_offset_element(i-1,j-1) what_b \
                                  right_op.get_offset_element(i-1,j-1) ) );  \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(10);                                        \
                                                                        \
                    break;                                              \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_SYMM:                                                    \
        {                                                               \
            switch ( right_op.get_e_type() )                            \
            {                                                           \
                case E_EMPTY:                                           \
                {                                                       \
                    L_THROW(11);                                        \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_ZERO:                                            \
                case E_SYMM:                                            \
                case E_VZ:                                              \
                case E_VS:                                              \
                {                                                       \
                    size_row = left_op.get_effective_height();          \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= i ; j++ )                    \
                        {                                               \
                            result =                                    \
                            ( result what_a                             \
                            ( left_op.get_offset_element(i-1,j-1) what_b \
                              right_op.get_offset_element(i-1,j-1) ) );  \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_DIAG:                                            \
                {                                                       \
                    if ( right_op.get_e_offs() == 0 )                   \
                    {                                                   \
                        size_row = left_op.get_effective_height();      \
                                                                        \
                        for ( i = 1 ; i <= size_row ; i++ )             \
                        {                                               \
                            for ( j = 1 ; j <= i ; j++ )                \
                            {                                           \
                                result =                                \
                                ( result what_a                         \
                                ( left_op.get_offset_element(i-1,j-1) what_b \
                                  right_op.get_offset_element(i-1,j-1) ) );  \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    else                                                \
                    {                                                   \
                        size_row = left_op.get_effective_height();      \
                        size_col = left_op.get_effective_width();       \
                                                                        \
                        for ( i = 1 ; i <= size_row ; i++ )             \
                        {                                               \
                            for ( j = 1 ; j <= size_col ; j++ )         \
                            {                                           \
                                result =                                \
                                ( result what_a                         \
                                ( left_op.get_offset_element(i-1,j-1) what_b \
                                  right_op.get_offset_element(i-1,j-1) ) );  \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_LOWER:                                           \
                case E_UPPER:                                           \
                case E_ASYMM:                                           \
                {                                                       \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            result =                                    \
                            ( result what_a                             \
                            ( left_op.get_offset_element(i-1,j-1) what_b \
                              right_op.get_offset_element(i-1,j-1) ) );  \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(12);                                        \
                                                                        \
                    break;                                              \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_ASYMM:                                                   \
        {                                                               \
            switch ( right_op.get_e_type() )                            \
            {                                                           \
                case E_EMPTY:                                           \
                {                                                       \
                    L_THROW(13);                                        \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_ZERO:                                            \
                case E_SYMM:                                            \
                case E_DIAG:                                            \
                case E_LOWER:                                           \
                case E_UPPER:                                           \
                case E_ASYMM:                                           \
                case E_VZ:                                              \
                case E_VS:                                              \
                {                                                       \
                    size_row = left_op.get_effective_height();          \
                    size_col = left_op.get_effective_width();           \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            result =                                    \
                            ( result what_a                             \
                            ( left_op.get_offset_element(i-1,j-1) what_b \
                              right_op.get_offset_element(i-1,j-1) ) );  \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(14);                                        \
                                                                        \
                    break;                                              \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_VZ:                                                      \
        {                                                               \
            switch ( right_op.get_e_type() )                            \
            {                                                           \
                case E_EMPTY:                                           \
                {                                                       \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_ZERO:                                            \
                {                                                       \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_DIAG:                                            \
                {                                                       \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-roffset == j )                       \
                            {                                           \
                                result =                                \
                                ( result what_a                         \
                                ( left_op.get_offset_element(i-1,j-1) what_b \
                                  right_op.get_offset_element(i-1,j-1) ) );  \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_LOWER:                                           \
                {                                                       \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-roffset >= j )                       \
                            {                                           \
                                result =                                \
                                ( result what_a                         \
                                ( left_op.get_offset_element(i-1,j-1) what_b \
                                  right_op.get_offset_element(i-1,j-1) ) );  \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_UPPER:                                           \
                {                                                       \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i-roffset <= j )                       \
                            {                                           \
                                result =                                \
                                ( result what_a                         \
                                ( left_op.get_offset_element(i-1,j-1) what_b \
                                  right_op.get_offset_element(i-1,j-1) ) );  \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_SYMM:                                            \
                {                                                       \
                    size_row = right_op.get_effective_height();         \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= i ; j++ )                    \
                        {                                               \
                            result =                                    \
                            ( result what_a                             \
                            ( left_op.get_offset_element(i-1,j-1) what_b \
                              right_op.get_offset_element(i-1,j-1) ) );  \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_ASYMM:                                           \
                {                                                       \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            result =                                    \
                            ( result what_a                             \
                            ( left_op.get_offset_element(i-1,j-1) what_b \
                              right_op.get_offset_element(i-1,j-1) ) );  \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_VZ:                                              \
                {                                                       \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_VS:                                              \
                {                                                       \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    result =                                            \
                    ( result what_a                                     \
                    ( left_op.get_offset_element(0,0) what_b            \
                      right_op.get_offset_element(-1,-1) ) );           \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(15);                                        \
                                                                        \
                    break;                                              \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_VS:                                                      \
        {                                                               \
            switch ( right_op.get_e_type() )                            \
            {                                                           \
                case E_EMPTY:                                           \
                {                                                       \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_ZERO:                                            \
                {                                                       \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( i == j )                               \
                            {                                           \
                                result =                                \
                                ( result what_a                         \
                                ( left_op.get_offset_element(i-1,j-1) what_b \
                                  right_op.get_offset_element(i-1,j-1) ) );  \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_DIAG:                                            \
                {                                                       \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( ( i         == j ) ||                  \
                                 ( i-roffset == j )    )                \
                            {                                           \
                                result =                                \
                                ( result what_a                         \
                                ( left_op.get_offset_element(i-1,j-1) what_b \
                                  right_op.get_offset_element(i-1,j-1) ) );  \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_LOWER:                                           \
                {                                                       \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( ( i         == j ) ||                  \
                                 ( i-roffset >= j )    )                \
                            {                                           \
                                result =                                \
                                ( result what_a                         \
                                ( left_op.get_offset_element(i-1,j-1) what_b \
                                  right_op.get_offset_element(i-1,j-1) ) );  \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_UPPER:                                           \
                {                                                       \
                    roffset = right_op.get_e_offs();                    \
                                                                        \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            if ( ( i         == j ) ||                  \
                                 ( i-roffset <= j )    )                \
                            {                                           \
                                result =                                \
                                ( result what_a                         \
                                ( left_op.get_offset_element(i-1,j-1) what_b \
                                  right_op.get_offset_element(i-1,j-1) ) );  \
                            }                                           \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_SYMM:                                            \
                {                                                       \
                    size_row = right_op.get_effective_height();         \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= i ; j++ )                    \
                        {                                               \
                            result =                                    \
                            ( result what_a                             \
                            ( left_op.get_offset_element(i-1,j-1) what_b \
                              right_op.get_offset_element(i-1,j-1) ) );  \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_ASYMM:                                           \
                {                                                       \
                    size_row = right_op.get_effective_height();         \
                    size_col = right_op.get_effective_width();          \
                                                                        \
                    for ( i = 1 ; i <= size_row ; i++ )                 \
                    {                                                   \
                        for ( j = 1 ; j <= size_col ; j++ )             \
                        {                                               \
                            result =                                    \
                            ( result what_a                             \
                            ( left_op.get_offset_element(i-1,j-1) what_b \
                              right_op.get_offset_element(i-1,j-1) ) );  \
                        }                                               \
                    }                                                   \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_VZ:                                              \
                {                                                       \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    result =                                            \
                    ( result what_a                                     \
                    ( left_op.get_offset_element(-1,-1) what_b          \
                      right_op.get_offset_element(0,0) ) );             \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                case E_VS:                                              \
                {                                                       \
                    result = ( result what_a ( 0.0 what_b 0.0 ) );      \
                                                                        \
                    result =                                            \
                    ( result what_a                                     \
                    ( left_op.get_offset_element(-1,-1) what_b          \
                      right_op.get_offset_element(-1,-1) ) );           \
                                                                        \
                                                                        \
                    break;                                              \
                }                                                       \
                                                                        \
                default:                                                \
                {                                                       \
                    L_THROW(18);                                        \
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
            L_THROW(19);                                                \
                                                                        \
            break;                                                      \
        }                                                               \
    }                                                                   \
}

#define MATRIX_SCALAR_RELATION(what_a,what_b,what_c)                    \
{                                                                       \
    long i,j;                                                           \
    long size_row;                                                      \
    long size_col;                                                      \
    long loffset;                                                       \
    long roffset;                                                       \
                                                                        \
    result = what_c;                                                    \
                                                                        \
    switch ( left_op.get_e_type() )                                     \
    {                                                                   \
        case E_EMPTY:                                                   \
        {                                                               \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_ZERO:                                                    \
        {                                                               \
            result = ( result what_a ( 0.0 what_b right_op ) );         \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_DIAG:                                                    \
        {                                                               \
            loffset = left_op.get_e_offs();                             \
                                                                        \
            size_row = left_op.get_effective_height();                  \
            size_col = left_op.get_effective_width();                   \
                                                                        \
            result = ( result what_a ( 0.0 what_b right_op ) );         \
                                                                        \
            for ( i = 1 ; i <= size_row ; i++ )                         \
            {                                                           \
                for ( j = 1 ; j <= size_col ; j++ )                     \
                {                                                       \
                    if ( i-loffset == j )                               \
                    {                                                   \
                        result =                                        \
                        ( result what_a                                 \
                        ( left_op.get_offset_element(i-1,j-1) what_b right_op ) ); \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_LOWER:                                                   \
        {                                                               \
            loffset = left_op.get_e_offs();                             \
                                                                        \
            size_row = left_op.get_effective_height();                  \
            size_col = left_op.get_effective_width();                   \
                                                                        \
            result = ( result what_a ( 0.0 what_b right_op ) );         \
                                                                        \
            for ( i = 1 ; i <= size_row ; i++ )                         \
            {                                                           \
                for ( j = 1 ; j <= size_col ; j++ )                     \
                {                                                       \
                    if ( i-loffset >= j )                               \
                    {                                                   \
                        result =                                        \
                        ( result what_a                                 \
                        ( left_op.get_offset_element(i-1,j-1) what_b right_op ) ); \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_UPPER:                                                   \
        {                                                               \
            loffset = left_op.get_e_offs();                             \
                                                                        \
            size_row = left_op.get_effective_height();                  \
            size_col = left_op.get_effective_width();                   \
                                                                        \
            result = ( result what_a ( 0.0 what_b right_op ) );         \
                                                                        \
            for ( i = 1 ; i <= size_row ; i++ )                         \
            {                                                           \
                for ( j = 1 ; j <= size_col ; j++ )                     \
                {                                                       \
                    if ( i-loffset <= j )                               \
                    {                                                   \
                        result =                                        \
                        ( result what_a                                 \
                        ( left_op.get_offset_element(i-1,j-1) what_b right_op ) ); \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_SYMM:                                                    \
        {                                                               \
            size_row = left_op.get_effective_height();                  \
                                                                        \
            for ( i = 1 ; i <= size_row ; i++ )                         \
            {                                                           \
                for ( j = 1 ; j <= i ; j++ )                            \
                {                                                       \
                    result =                                            \
                    ( result what_a                                     \
                    ( left_op.get_offset_element(i-1,j-1) what_b right_op ) ); \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_ASYMM:                                                   \
        {                                                               \
            size_row = left_op.get_effective_height();                  \
            size_col = left_op.get_effective_width();                   \
                                                                        \
            for ( i = 1 ; i <= size_row ; i++ )                         \
            {                                                           \
                for ( j = 1 ; j <= size_col ; j++ )                     \
                {                                                       \
                    result =                                            \
                    ( result what_a                                     \
                    ( left_op.get_offset_element(i-1,j-1) what_b right_op ) ); \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_VZ:                                                      \
        {                                                               \
            result = ( result what_a ( 0.0 what_b right_op ) );         \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_VS:                                                      \
        {                                                               \
            result = ( result what_a ( 0.0 what_b right_op ) );         \
            result = ( result what_a ( left_op.get_offset_element(-1,-1) what_b right_op ) ); \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        default:                                                        \
        {                                                               \
            L_THROW(19);                                                \
                                                                        \
            break;                                                      \
        }                                                               \
    }                                                                   \
}







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

std::ostream &operator<<(std::ostream &output, const fMATRIX   &source)
{
    FN_ENTRY_POINT

    long i,j;

    output << "size:   " << source.size   << "\n";
    output << "m_type: " << source.m_type << "\n\n";

    if ( source.size > 0 )
    {
        output << "start_row: " << source.offset_start_effective_row << "\n";
        output << "start_col: " << source.offset_start_effective_col << "\n\n";

        output << "end_row: " << source.offset_end_effective_row << "\n";
        output << "end_col: " << source.offset_end_effective_col << "\n\n";

        output << "Matrix:\n";

        for ( i = 1 ; i <= source.size ; i++ )
        {
            for ( j = 1 ; j <= source.size ; j++ )
            {
                output << fixnum(source.getelm(i,j)) << " ";
            }

            output << "\n";
        }
    }

    if ( source.size < 0 )
    {
        output << "Matrix:\n";

        output << fixnum(source.m_diag) << "\n";
    }

    FN_EXIT_POINT output;
}

std::istream &operator>>(std::istream &input, fMATRIX   &dest)
{
    FN_ENTRY_POINT

    long i,j;
    L_DOUBLE dummy;
    wait_dummy howzat;
    long size;
    long m_type;

    std::ostream *where_to;
    int echo_level;

    where_to   = dest.where_to;
    echo_level = dest.echo_level;

    dest.prime_io(NULL,0);

    if ( where_to == NULL )
    {
        dest.make_zero(DO_FORCE);

        input >> howzat; input >> size;
        input >> howzat; input >> m_type;

        switch ( m_type )
        {
            case ASYMMETRIC_MATRIX:
            {
                dest.make_asymmetric(DO_FORCE);

                break;
            }

            case SYMMETRIC_MATRIX:
            {
                dest.make_symmetric(DO_FORCE);

                break;
            }

            case LOWER_TRIANGULAR_MATRIX:
            {
                dest.make_lower_triangular(DO_FORCE);

                break;
            }

            case UPPER_TRIANGULAR_MATRIX:
            {
                dest.make_upper_triangular(DO_FORCE);

                break;
            }

            case DIAGONAL_MATRIX:
            {
                dest.make_diagonal(DO_FORCE);

                break;
            }

            case ZERO_MATRIX:
            {
                dest.make_zero(DO_FORCE);

                break;
            }

            case SCALAR_MATRIX:
            {
                dest.make_scalar(DO_FORCE);

                break;
            }

            default:
            {
                L_THROW(0);

                break;
            }
        }

        if ( size > 0 )
        {
            dest.pad_matrix(size);

            input >> howzat; input >> dest.offset_start_effective_row;
            input >> howzat; input >> dest.offset_start_effective_col;

            input >> howzat; input >> dest.offset_end_effective_row;
            input >> howzat; input >> dest.offset_end_effective_col;

            input >> howzat;

            for ( i = 1 ; i <= size ; i++ )
            {
                for ( j = 1 ; j <= size ; j++ )
                {
                    if ( i == j )
                    {
                        input >> dest.getref(i,j);
                    }

                    else
                    {
                        if ( i > j )
                        {
                            if ( ( m_type != UPPER_TRIANGULAR_MATRIX ) &&
                                 ( m_type != DIAGONAL_MATRIX         )    )
                            {
                                input >> dest.getref(i,j);
                            }

                            else
                            {
                                input >> dummy;
                            }
                        }

                        else
                        {
                            if ( ( m_type != LOWER_TRIANGULAR_MATRIX ) &&
                                 ( m_type != DIAGONAL_MATRIX         )    )
                            {
                                input >> dest.getref(i,j);
                            }

                            else
                            {
                                input >> dummy;
                            }
                        }
                    }
                }
            }
        }

        if ( size < 0 )
        {
            input >> howzat; input >> dest.m_diag;
        }
    }

    else
    {
        dest.make_zero(DO_FORCE);

        *where_to << "Matrix types: 0 = Asymmetric matrix.\n";
        *where_to << "              1 = Symmetric matrix.\n";
        *where_to << "              2 = Lower triangular matrix.\n";
        *where_to << "              3 = Upper triangular matrix.\n";
        *where_to << "              4 = Diagonal matrix.\n";
        *where_to << "              5 = Scalar matrix.\n";
        *where_to << "              6 = Zero matrix.\n\n";

        *where_to << "Matrix type: ";
        input >> m_type;
        if ( echo_level ) { *where_to << m_type << "\n"; }

        switch ( m_type )
        {
            case ASYMMETRIC_MATRIX:
            {
                *where_to << "Size: ";
                input >> size;
                if ( echo_level ) { *where_to << size << "\n"; }

                dest.make_asymmetric(DO_FORCE);

                break;
            }

            case SYMMETRIC_MATRIX:
            {
                *where_to << "Size: ";
                input >> size;
                if ( echo_level ) { *where_to << size << "\n"; }

                dest.make_symmetric(DO_FORCE);

                break;
            }

            case LOWER_TRIANGULAR_MATRIX:
            {
                *where_to << "Size: ";
                input >> size;
                if ( echo_level ) { *where_to << size << "\n"; }

                dest.make_lower_triangular(DO_FORCE);

                break;
            }

            case UPPER_TRIANGULAR_MATRIX:
            {
                *where_to << "Size: ";
                input >> size;
                if ( echo_level ) { *where_to << size << "\n"; }

                dest.make_upper_triangular(DO_FORCE);

                break;
            }

            case DIAGONAL_MATRIX:
            {
                *where_to << "Size: ";
                input >> size;
                if ( echo_level ) { *where_to << size << "\n"; }

                dest.make_diagonal(DO_FORCE);

                break;
            }

            case SCALAR_MATRIX:
            {
                size = SCALAR_MATRIX_SIZE;

                dest.make_scalar(DO_FORCE);

                break;
            }

            case ZERO_MATRIX:
            {
                size = ZERO_MATRIX_SIZE;

                dest.make_zero(DO_FORCE);

                break;
            }

            default:
            {
                L_THROW(0);

                break;
            }
        }

        if ( size > 0 )
        {
            dest.pad_matrix(size);

            dest.offset_start_effective_row = 1;
            dest.offset_start_effective_col = 1;

            dest.offset_end_effective_row = size;
            dest.offset_end_effective_col = size;

            switch ( m_type )
            {
                case ASYMMETRIC_MATRIX:
                {
                    for ( i = 1 ; i <= size ; i++ )
                    {
                        *where_to << "Row " << i << " (1 - " << size << "): ";

                        for ( j = 1 ; j <= size ; j++ )
                        {
                            input >> dest.getref(i,j);
                            if ( echo_level ) { *where_to << dest.getref(i,j) << " "; }
                        }

                        if ( echo_level ) { *where_to << "\n"; }
                    }

                    break;
                }

                case SYMMETRIC_MATRIX:
                {
                    for ( i = 1 ; i <= size ; i++ )
                    {
                        *where_to << "Row " << i << " (1 - " << i << "): ";

                        for ( j = 1 ; j <= i ; j++ )
                        {
                            input >> dest.getref(i,j);
                            if ( echo_level ) { *where_to << dest.getref(i,j) << " "; }
                        }

                        if ( echo_level ) { *where_to << "\n"; }
                    }

                    break;
                }

                case LOWER_TRIANGULAR_MATRIX:
                {
                    for ( i = 1 ; i <= size ; i++ )
                    {
                        *where_to << "Row " << i << " (1 - " << i << "): ";

                        for ( j = 1 ; j <= i ; j++ )
                        {
                            input >> dest.getref(i,j);
                            if ( echo_level ) { *where_to << dest.getref(i,j) << " "; }
                        }

                        if ( echo_level ) { *where_to << "\n"; }
                    }

                    break;
                }

                case UPPER_TRIANGULAR_MATRIX:
                {
                    for ( i = 1 ; i <= size ; i++ )
                    {
                        *where_to << "Row " << i << " (" << i << " - " << size << "): ";

                        for ( j = i ; j <= size ; j++ )
                        {
                            input >> dest.getref(i,j);
                            if ( echo_level ) { *where_to << dest.getref(i,j) << " "; }
                        }

                        if ( echo_level ) { *where_to << "\n"; }
                    }

                    break;
                }

                case DIAGONAL_MATRIX:
                {
                    for ( i = 1 ; i <= size ; i++ )
                    {
                        *where_to << "Element (" << i << "," << i << "): ";

                        input >> dest.getref(i,i);
                        if ( echo_level ) { *where_to << dest.getref(i,i) << "\n"; }
                    }

                    break;
                }

                case SCALAR_MATRIX:
                case ZERO_MATRIX:
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

        if ( size == SCALAR_MATRIX_SIZE )
        {
            *where_to << "Scalar: ";
            input >> dest.m_diag;
            if ( echo_level ) { *where_to << dest.m_diag << "\n"; }
        }
    }

    FN_EXIT_POINT input;
}




//
// Mathematical operator overloading
//

// + posation - unary, return rvalue
// - negation - unary, return rvalue

fMATRIX   operator+ (const fMATRIX  &left_op)
{
    FN_ENTRY_POINT

    fMATRIX result;

    UNARY_MATRIX_OPERATION_CONST_RES(+,);

    FN_EXIT_POINT result;
}


fMATRIX   operator- (const fMATRIX  &left_op)
{
    FN_ENTRY_POINT

    fMATRIX result;

    UNARY_MATRIX_OPERATION_CONST_RES(-,);

    FN_EXIT_POINT result;
}




// + addition       - binary, return rvalue
// - subtraction    - binary, return rvalue

fMATRIX   operator+ (const fMATRIX  &left_op, const fMATRIX  &right_op)
{
    FN_ENTRY_POINT

    fMATRIX result;

    BINARY_MATRIX_OPERATION_CONST_CONST_RES(+);

    FN_EXIT_POINT result;
}


fMATRIX   operator- (const fMATRIX  &left_op, const fMATRIX  &right_op)
{
    FN_ENTRY_POINT

    fMATRIX result;

    BINARY_MATRIX_OPERATION_CONST_CONST_RES(-);

    FN_EXIT_POINT result;
}




// += additive       assignment - binary, return lvalue
// -= subtractive    assignment - binary, return lvalue

fMATRIX  &operator+=(      fMATRIX  &left_op, const fMATRIX  &right_op)
{
    FN_ENTRY_POINT

    BINARY_MATRIX_OPERATION_VARIA_CONST(+=);

    FN_EXIT_POINT left_op;
}


fMATRIX  &operator-=(      fMATRIX  &left_op, const fMATRIX  &right_op)
{
    FN_ENTRY_POINT

    BINARY_MATRIX_OPERATION_VARIA_CONST(-=);

    FN_EXIT_POINT left_op;
}




// * multiplication - binary, return rvalue
// / division       - binary, return rvalue

fMATRIX   operator* (const fMATRIX  &left_op, const   double &right_op)
{
    FN_ENTRY_POINT

    fMATRIX result;

    UNARY_MATRIX_OPERATION_CONST_RES(,*right_op);

    FN_EXIT_POINT result;
}

fMATRIX   operator* (const fMATRIX  &left_op, const c_double &right_op)
{
    FN_ENTRY_POINT

    fMATRIX result;

    UNARY_MATRIX_OPERATION_CONST_RES(,*right_op);

    FN_EXIT_POINT result;
}


fMATRIX   operator* (const   double &right_op, const fMATRIX  &left_op)
{
    FN_ENTRY_POINT

    fMATRIX result;

    UNARY_MATRIX_OPERATION_CONST_RES(right_op*,);

    FN_EXIT_POINT result;
}

fMATRIX   operator* (const c_double &right_op, const fMATRIX  &left_op)
{
    FN_ENTRY_POINT

    fMATRIX result;

    UNARY_MATRIX_OPERATION_CONST_RES(right_op*,);

    FN_EXIT_POINT result;
}



fMATRIX   operator/ (const fMATRIX  &left_op, const   double &right_op)
{
    FN_ENTRY_POINT

    fMATRIX result;

    UNARY_MATRIX_OPERATION_CONST_RES(,/right_op);

    FN_EXIT_POINT result;
}

fMATRIX   operator/ (const fMATRIX  &left_op, const c_double &right_op)
{
    FN_ENTRY_POINT

    fMATRIX result;

    UNARY_MATRIX_OPERATION_CONST_RES(,/right_op);

    FN_EXIT_POINT result;
}







// * multiplication - binary, FN_EXIT_POINT rvalue

fMATRIX   operator* (const fMATRIX  &left_op, const fMATRIX  &right_op)
{
    FN_ENTRY_POINT

    fMATRIX result;

    long i,j,k;                                                         
    long height;                                                        
    long width;                                                         
    long length;                                                        
    long size;                                                          
                                                                        
    switch ( left_op.get_e_type() )                                     
    {                                                                   
        case E_EMPTY:                                                   
        {                                                               
            switch ( right_op.get_e_type() )                            
            {                                                           
                case E_EMPTY:                                           
                case E_VZ:                                              
                case E_VS:                                              
                {                                                       
                    result.make_diagonal(DO_FORCE);                     
                                                                        
                    height = 0;                                         
                    width  = 0;                                         
                                                                        
                    size = 0;                                           
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case E_ZERO:                                            
                case E_DIAG:                                            
                case E_LOWER:                                           
                case E_UPPER:                                           
                case E_SYMM:                                            
                case E_ASYMM:                                           
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
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case E_ZERO:                                                    
        {                                                               
            switch ( right_op.get_e_type() )                            
            {                                                           
                case E_EMPTY:                                           
                {                                                       
                    L_THROW(0);                                         
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case E_ZERO:                                            
                case E_DIAG:                                            
                case E_LOWER:                                           
                case E_UPPER:                                           
                case E_SYMM:                                            
                case E_ASYMM:                                           
                {                                                       
                    THROW_ASSERT((left_op.get_effective_width()) == (right_op.get_effective_height())); 
                                                                        
                    height = left_op.get_effective_height();            
                    width  = right_op.get_effective_width();            
                                                                        
                    size = height;                                      
                                                                        
                    if ( width > height )                               
                    {                                                   
                        size = width;                                   
                    }                                                   
                                                                        
                    result.make_diagonal(DO_FORCE);                     
                    result.pad_matrix(size);                            
                                                                        
                    result.set_offset_end_at_row(height);               
                    result.set_offset_end_at_col(width);                
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case E_VZ:                                              
                case E_VS:                                              
                {                                                       
                    height = left_op.get_effective_height();            
                    width  = left_op.get_effective_width();             
                                                                        
                    size = height;                                      
                                                                        
                    result.make_diagonal(DO_FORCE);                     
                    result.pad_matrix(size);                            
                                                                        
                    result.set_offset_end_at_row(height);               
                    result.set_offset_end_at_col(width);                
                                                                        
                    break;                                              
                }                                                       
                                                                        
                default:                                                
                {                                                       
                    L_THROW(0);                                         
                                                                        
                    break;                                              
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case E_DIAG:                                                    
        case E_LOWER:                                                   
        case E_UPPER:                                                   
        case E_ASYMM:                                                   
        {                                                               
            switch ( right_op.get_e_type() )                            
            {                                                           
                case E_EMPTY:                                           
                {                                                       
                    L_THROW(0);                                         
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case E_ZERO:                                            
                {                                                       
                    height = left_op.get_effective_height();            
                    width  = right_op.get_effective_width();            
                    length = left_op.get_effective_width();             
                                                                        
                    THROW_ASSERT(right_op.get_effective_height() == length);  
                                                                        
                    size = height;                                      
                                                                        
                    if ( width > height )                               
                    {                                                   
                        size = width;                                   
                    }                                                   
                                                                        
                    result.make_diagonal(DO_FORCE);                     
                    result.pad_matrix(size);                            
                                                                        
                    result.set_offset_end_at_row(height);               
                    result.set_offset_end_at_col(width);                
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case E_DIAG:                                            
                case E_LOWER:                                           
                case E_UPPER:                                           
                case E_SYMM:                                            
                case E_ASYMM:                                           
                {                                                       
                    height = left_op.get_effective_height();            
                    width  = right_op.get_effective_width();            
                    length = left_op.get_effective_width();             
                                                                        
                    THROW_ASSERT(right_op.get_effective_height() == length);  
                                                                        
                    size = height;                                      
                                                                        
                    if ( width > height )                               
                    {                                                   
                        size = width;                                   
                    }                                                   
                                                                        
                    result.make_diagonal(DO_FORCE);                     
                    result.pad_matrix(size);                            
                                                                        
                    result.set_offset_end_at_row(height);               
                    result.set_offset_end_at_col(width);                
                                                                        
                    for ( i = 1 ; i <= height ; i++ )                   
                    {                                                   
                        for ( j = 1 ; j <= width ; j++ )                
                        {                                               
                            for ( k = 1 ; k <= length ; k++ )           
                            {                                           
                                if ( k == 1 )                           
                                {                                       
                                    result(i-1,j-1) = (left_op.get_offset_element(i-1,k-1) * right_op.get_offset_element(k-1,j-1)); 
                                }                                       
                                                                        
                                else                                    
                                {                                       
                                    result(i-1,j-1) += (left_op.get_offset_element(i-1,k-1) * right_op.get_offset_element(k-1,j-1)); 
                                }                                       
                            }                                           
                        }                                               
                    }                                                   
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case E_VZ:                                              
                {                                                       
                    height = left_op.get_effective_height();            
                    width  = left_op.get_effective_width();             
                                                                        
                    size = height;                                      
                                                                        
                    if ( width > height )                               
                    {                                                   
                        size = width;                                   
                    }                                                   
                                                                        
                    result.make_diagonal(DO_FORCE);                     
                    result.pad_matrix(size);                            
                                                                        
                    result.set_offset_end_at_row(height);               
                    result.set_offset_end_at_col(width);                
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case E_VS:                                              
                {                                                       
                    height = left_op.get_effective_height();            
                    width  = left_op.get_effective_width();             
                                                                        
                    size = height;                                      
                                                                        
                    if ( width > height )                               
                    {                                                   
                        size = width;                                   
                    }                                                   
                                                                        
                    result = (left_op * right_op.get_offset_element(-1,-1)); 
                                                                        
                    result.set_offset_end_at_row(height);               
                    result.set_offset_end_at_col(width);                
                                                                        
                    break;                                              
                }                                                       
                                                                        
                default:                                                
                {                                                       
                    L_THROW(0);                                         
                                                                        
                    break;                                              
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case E_SYMM:                                                    
        {                                                               
            switch ( right_op.get_e_type() )                            
            {                                                           
                case E_EMPTY:                                           
                {                                                       
                    L_THROW(0);                                         
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case E_ZERO:                                            
                {                                                       
                    THROW_ASSERT((left_op.get_effective_width()) == (right_op.get_effective_height())); 
                                                                        
                    height = left_op.get_effective_height();            
                    width  = right_op.get_effective_width();            
                                                                        
                    size = height;                                      
                                                                        
                    if ( width > height )                               
                    {                                                   
                        size = width;                                   
                    }                                                   
                                                                        
                    result.make_diagonal(DO_FORCE);                     
                    result.pad_matrix(size);                            
                                                                        
                    result.set_offset_end_at_row(height);               
                    result.set_offset_end_at_col(width);                
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case E_DIAG:                                            
                case E_LOWER:                                           
                case E_UPPER:                                           
                case E_ASYMM:                                           
                {                                                       
                    height = left_op.get_effective_height();            
                    width  = right_op.get_effective_width();            
                    length = left_op.get_effective_width();             
                                                                        
                    THROW_ASSERT(right_op.get_effective_height() == length);  
                                                                        
                    size = height;                                      
                                                                        
                    if ( width > height )                               
                    {                                                   
                        size = width;                                   
                    }                                                   
                                                                        
                    result.make_diagonal(DO_FORCE);                     
                    result.pad_matrix(size);                            
                                                                        
                    result.set_offset_end_at_row(height);               
                    result.set_offset_end_at_col(width);                
                                                                        
                    for ( i = 1 ; i <= height ; i++ )                   
                    {                                                   
                        for ( j = 1 ; j <= width ; j++ )                
                        {                                               
                            for ( k = 1 ; k <= length ; k++ )           
                            {                                           
                                if ( k == 1 )                           
                                {                                       
                                    result(i-1,j-1) = left_op.get_offset_element(i-1,k-1) * right_op.get_offset_element(k-1,j-1); 
                                }                                       
                                                                        
                                else                                    
                                {                                       
                                    result(i-1,j-1) += left_op.get_offset_element(i-1,k-1) * right_op.get_offset_element(k-1,j-1); 
                                }                                       
                            }                                           
                        }                                               
                    }                                                   
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case E_SYMM:                                            
                {                                                       
                    height = left_op.get_effective_height();            
                    width  = right_op.get_effective_width();            
                    length = left_op.get_effective_width();             
                                                                        
                    THROW_ASSERT(right_op.get_effective_height() == length);  
                                                                        
                    size = height;                                      
                                                                        
                    if ( width > height )                               
                    {                                                   
                        size = width;                                   
                    }                                                   
                                                                        
                    result.make_symmetric(DO_FORCE);                    
                    result.pad_matrix(size);                            
                                                                        
                    result.set_offset_end_at_row(height);               
                    result.set_offset_end_at_col(width);                
                                                                        
                    for ( i = 1 ; i <= height ; i++ )                   
                    {                                                   
                        for ( j = 1 ; j <= i ; j++ )                    
                        {                                               
                            for ( k = 1 ; k <= length ; k++ )           
                            {                                           
                                if ( k == 1 )                           
                                {                                       
                                    result(i-1,j-1) = left_op.get_offset_element(i-1,k-1) * right_op.get_offset_element(k-1,j-1); 
                                }                                       
                                                                        
                                else                                    
                                {                                       
                                    result(i-1,j-1) += left_op.get_offset_element(i-1,k-1) * right_op.get_offset_element(k-1,j-1); 
                                }                                       
                            }                                           
                        }                                               
                    }                                                   
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case E_VZ:                                              
                {                                                       
                    height = left_op.get_effective_height();            
                    width  = left_op.get_effective_width();             
                                                                        
                    size = height;                                      
                                                                        
                    if ( width > height )                               
                    {                                                   
                        size = width;                                   
                    }                                                   
                                                                        
                    result.make_symmetric(DO_FORCE);                    
                    result.pad_matrix(size);                            
                                                                        
                    result.set_offset_end_at_row(height);               
                    result.set_offset_end_at_col(width);                
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case E_VS:                                              
                {                                                       
                    height = left_op.get_effective_height();            
                    width  = left_op.get_effective_width();             
                                                                        
                    size = height;                                      
                                                                        
                    if ( width > height )                               
                    {                                                   
                        size = width;                                   
                    }                                                   
                                                                        
                    result = ( left_op * right_op.get_offset_element(-1,-1) ); 
                                                                        
                    result.set_offset_end_at_row(height);               
                    result.set_offset_end_at_col(width);                
                                                                        
                    break;                                              
                }                                                       
                                                                        
                default:                                                
                {                                                       
                    L_THROW(0);                                         
                                                                        
                    break;                                              
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case E_VZ:                                                      
        {                                                               
            switch ( right_op.get_e_type() )                            
            {                                                           
                case E_EMPTY:                                           
                {                                                       
                    result.make_diagonal(DO_FORCE);                     
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case E_ZERO:                                            
                case E_DIAG:                                            
                case E_LOWER:                                           
                case E_UPPER:                                           
                case E_ASYMM:                                           
                {                                                       
                    height = right_op.get_effective_height();           
                    width  = right_op.get_effective_width();            
                                                                        
                    size = height;                                      
                                                                        
                    if ( width > height )                               
                    {                                                   
                        size = width;                                   
                    }                                                   
                                                                        
                    result.make_diagonal(DO_FORCE);                     
                    result.pad_matrix(size);                            
                                                                        
                    result.set_offset_end_at_row(height);               
                    result.set_offset_end_at_col(width);                
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case E_SYMM:                                            
                {                                                       
                    height = right_op.get_effective_height();           
                    width  = right_op.get_effective_width();            
                                                                        
                    size = height;                                      
                                                                        
                    if ( width > height )                               
                    {                                                   
                        size = width;                                   
                    }                                                   
                                                                        
                    result.make_symmetric(DO_FORCE);                    
                    result.pad_matrix(size);                            
                                                                        
                    result.set_offset_end_at_row(height);               
                    result.set_offset_end_at_col(width);                
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case E_VZ:                                              
                case E_VS:                                              
                {                                                       
                    result.make_zero(DO_FORCE);                         
                                                                        
                    break;                                              
                }                                                       
                                                                        
                default:                                                
                {                                                       
                    L_THROW(0);                                         
                                                                        
                    break;                                              
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case E_VS:                                                      
        {                                                               
            switch ( right_op.get_e_type() )                            
            {                                                           
                case E_EMPTY:                                           
                {                                                       
                    result.make_diagonal(DO_FORCE);                     
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case E_ZERO:                                            
                case E_DIAG:                                            
                case E_LOWER:                                           
                case E_UPPER:                                           
                case E_SYMM:                                            
                case E_ASYMM:                                           
                {                                                       
                    height = right_op.get_effective_height();           
                    width  = right_op.get_effective_width();            
                                                                        
                    size = height;                                      
                                                                        
                    if ( width > height )                               
                    {                                                   
                        size = width;                                   
                    }                                                   
                                                                        
                    result = (right_op * (left_op.get_offset_element(-1,-1)) ); 
                                                                        
                    result.set_offset_end_at_row(height);               
                    result.set_offset_end_at_col(width);                
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case E_VZ:                                              
                {                                                       
                    result.make_zero(DO_FORCE);                         
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case E_VS:                                              
                {                                                       
                    result = ( left_op * right_op.get_offset_element(-1,-1) ); 
                                                                        
                    break;                                              
                }                                                       
                                                                        
                default:                                                
                {                                                       
                    L_THROW(0);                                         
                                                                        
                    break;                                              
                }                                                       
            }                                                           
                                                                        
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

fVECTOR   operator* (const fMATRIX  &left_op, const fVECTOR  &right_op)
{
    FN_ENTRY_POINT

    fVECTOR result;

    long i,j;                                                           
    long size;                                                          
    long length;                                                        
                                                                        
    switch ( left_op.get_e_type() )                                     
    {                                                                   
        case E_EMPTY:                                                   
        {                                                               
            switch ( right_op.get_v_type() )                            
            {                                                           
                case NORMAL_VECTOR:                                     
                {                                                       
                    THROW_ASSERT(right_op.get_effective_size() == 0);         
                                                                        
                    result.make_normal(DO_FORCE);                       
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case ZERO_VECTOR:                                       
                case SCALAR_VECTOR:                                     
                case SELECT_VECTOR:                                     
                {                                                       
                    result.make_normal(DO_FORCE);                       
                                                                        
                    break;                                              
                }                                                       
                                                                        
                default:                                                
                {                                                       
                    L_THROW(0);                                         
                                                                        
                    break;                                              
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case E_ZERO:                                                    
        {                                                               
            switch ( right_op.get_v_type() )                            
            {                                                           
                case NORMAL_VECTOR:                                     
                {                                                       
                    THROW_ASSERT((left_op.get_effective_width()) == (right_op.get_effective_size())); 
                                                                        
                    result.make_normal(DO_FORCE);                       
                                                                        
                    result.pad_vector(left_op.get_effective_height());  
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case ZERO_VECTOR:                                       
                case SCALAR_VECTOR:                                     
                case SELECT_VECTOR:                                     
                {                                                       
                    result.make_normal(DO_FORCE);                       
                                                                        
                    result.pad_vector(left_op.get_effective_height());  
                                                                        
                    break;                                              
                }                                                       
                                                                        
                default:                                                
                {                                                       
                    L_THROW(0);                                         
                                                                        
                    break;                                              
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case E_DIAG:                                                    
        case E_LOWER:                                                   
        case E_UPPER:                                                   
        case E_SYMM:                                                    
        case E_ASYMM:                                                   
        {                                                               
            switch ( right_op.get_v_type() )                            
            {                                                           
                case NORMAL_VECTOR:                                     
                {                                                       
                    THROW_ASSERT((left_op.get_effective_width()) == (right_op.get_effective_size())); 
                                                                        
                    result.make_normal(DO_FORCE);                       
                                                                        
                    result.pad_vector(left_op.get_effective_height());  
                                                                        
                    size = left_op.get_effective_height();              
                    length = left_op.get_effective_width();             
                                                                        
                    for ( i = 1 ; i <= size ; i++ )                     
                    {                                                   
                        for ( j = 1 ; j <= length ; j++ )               
                        {                                               
                            if ( j == 1 )                               
                            {                                           
                                result[i-1] = left_op.get_offset_element(i-1,j-1) * right_op.get_offset_element(j-1); 
                            }                                           
                                                                        
                            else                                        
                            {                                           
                                result[i-1] += left_op.get_offset_element(i-1,j-1) * right_op.get_offset_element(j-1); 
                            }                                           
                        }                                               
                    }                                                   
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case SCALAR_VECTOR:                                     
                {                                                       
                    result.make_normal(DO_FORCE);                       
                                                                        
                    result.pad_vector(left_op.get_effective_height());  
                                                                        
                    size = left_op.get_effective_height();              
                    length = left_op.get_effective_width();             
                                                                        
                    for ( i = 1 ; i <= size ; i++ )                     
                    {                                                   
                        for ( j = 1 ; j <= length ; j++ )               
                        {                                               
                            if ( j == 1 )                               
                            {                                           
                                result[i-1] = left_op.get_offset_element(i-1,j-1) * right_op.get_offset_element(j-1); 
                            }                                           
                                                                        
                            else                                        
                            {                                           
                                result[i-1] += left_op.get_offset_element(i-1,j-1) * right_op.get_offset_element(j-1); 
                            }                                           
                        }                                               
                    }                                                   
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case SELECT_VECTOR:                                     
                {                                                       
                    result =  left_op((right_op.get_sel_elm())-1);      
                    result *= right_op.get_offset_element((right_op.get_sel_elm())-1); 
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case ZERO_VECTOR:                                       
                {                                                       
                    result.make_normal(DO_FORCE);                       
                                                                        
                    result.pad_vector(left_op.get_effective_height());  
                                                                        
                    break;                                              
                }                                                       
                                                                        
                default:                                                
                {                                                       
                    L_THROW(0);                                         
                                                                        
                    break;                                              
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case E_VZ:                                                      
        {                                                               
            switch ( right_op.get_v_type() )                            
            {                                                           
                case NORMAL_VECTOR:                                     
                case ZERO_VECTOR:                                       
                case SCALAR_VECTOR:                                     
                case SELECT_VECTOR:                                     
                {                                                       
                    result = left_op.get_offset_element(0,0) * right_op; 
                                                                        
                    break;                                              
                }                                                       
                                                                        
                default:                                                
                {                                                       
                    L_THROW(0);                                         
                                                                        
                    break;                                              
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case E_VS:                                                      
        {                                                               
            switch ( right_op.get_v_type() )                            
            {                                                           
                case NORMAL_VECTOR:                                     
                case ZERO_VECTOR:                                       
                case SCALAR_VECTOR:                                     
                case SELECT_VECTOR:                                     
                {                                                       
                    result = left_op.get_offset_element(0,0) * right_op; 
                                                                        
                    break;                                              
                }                                                       
                                                                        
                default:                                                
                {                                                       
                    L_THROW(0);                                         
                                                                        
                    break;                                              
                }                                                       
            }                                                           
                                                                        
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

fVECTOR   operator* (const fVECTOR  &left_op, const fMATRIX  &right_op)
{
    FN_ENTRY_POINT

    fVECTOR result;

    long i,j;                                                           
    long size;                                                          
    long length;                                                        
                                                                        
    switch ( right_op.get_e_type() )                                    
    {                                                                   
        case E_EMPTY:                                                   
        {                                                               
            switch ( left_op.get_v_type() )                             
            {                                                           
                case NORMAL_VECTOR:                                     
                {                                                       
                    THROW_ASSERT(left_op.get_effective_size() == 0);          
                                                                        
                    result.make_normal(DO_FORCE);                       
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case ZERO_VECTOR:                                       
                case SCALAR_VECTOR:                                     
                case SELECT_VECTOR:                                     
                {                                                       
                    result.make_normal(DO_FORCE);                       
                                                                        
                    break;                                              
                }                                                       
                                                                        
                default:                                                
                {                                                       
                    L_THROW(0);                                         
                                                                        
                    break;                                              
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case E_ZERO:                                                    
        {                                                               
            switch ( left_op.get_v_type() )                             
            {                                                           
                case NORMAL_VECTOR:                                     
                {                                                       
                    THROW_ASSERT((right_op.get_effective_height()) == (left_op.get_effective_size())); 
                                                                        
                    result.make_normal(DO_FORCE);                       
                                                                        
                    result.pad_vector(right_op.get_effective_width());  
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case ZERO_VECTOR:                                       
                case SCALAR_VECTOR:                                     
                case SELECT_VECTOR:                                     
                {                                                       
                    result.make_normal(DO_FORCE);                       
                                                                        
                    result.pad_vector(right_op.get_effective_width());  
                                                                        
                    break;                                              
                }                                                       
                                                                        
                default:                                                
                {                                                       
                    L_THROW(0);                                         
                                                                        
                    break;                                              
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case E_DIAG:                                                    
        case E_LOWER:                                                   
        case E_UPPER:                                                   
        case E_SYMM:                                                    
        case E_ASYMM:                                                   
        {                                                               
            size = right_op.get_effective_width();                      
            length = right_op.get_effective_height();                   
                                                                        
            switch ( left_op.get_v_type() )                             
            {                                                           
                case NORMAL_VECTOR:                                     
                {                                                       
                    THROW_ASSERT((right_op.get_effective_height()) == (left_op.get_effective_size())); 
                                                                        
                    result.make_normal(DO_FORCE);                       
                                                                        
                    result.pad_vector(right_op.get_effective_width());  
                                                                        
                    for ( i = 1 ; i <= size ; i++ )                     
                    {                                                   
                        for ( j = 1 ; j <= length ; j++ )               
                        {                                               
                            if ( j == 1 )                               
                            {                                           
                                result[i-1] = left_op.get_offset_element(j-1) * right_op.get_offset_element(j-1,i-1); 
                            }                                           
                                                                        
                            else                                        
                            {                                           
                                result[i-1] += left_op.get_offset_element(j-1) * right_op.get_offset_element(j-1,i-1); 
                            }                                           
                        }                                               
                    }                                                   
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case SCALAR_VECTOR:                                     
                {                                                       
                    result.make_normal(DO_FORCE);                       
                                                                        
                    result.pad_vector(right_op.get_effective_width());  
                                                                        
                    for ( i = 1 ; i <= size ; i++ )                     
                    {                                                   
                        for ( j = 1 ; j <= length ; j++ )               
                        {                                               
                            if ( j == 1 )                               
                            {                                           
                                result[i-1] = left_op.get_offset_element(j-1) * right_op.get_offset_element(j-1,i-1); 
                            }                                           
                                                                        
                            else                                        
                            {                                           
                                result[i-1] += left_op.get_offset_element(j-1) * right_op.get_offset_element(j-1,i-1); 
                            }                                           
                        }                                               
                    }                                                   
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case SELECT_VECTOR:                                     
                {                                                       
                    result = left_op.get_offset_element((left_op.get_sel_elm())-1) * right_op[(left_op.get_sel_elm())-1]; 
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case ZERO_VECTOR:                                       
                {                                                       
                    result.make_normal(DO_FORCE);                       
                                                                        
                    result.pad_vector(right_op.get_effective_width());  
                                                                        
                    break;                                              
                }                                                       
                                                                        
                default:                                                
                {                                                       
                    L_THROW(0);                                         
                                                                        
                    break;                                              
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case E_VZ:                                                      
        {                                                               
            switch ( left_op.get_v_type() )                             
            {                                                           
                case NORMAL_VECTOR:                                     
                case ZERO_VECTOR:                                       
                case SCALAR_VECTOR:                                     
                case SELECT_VECTOR:                                     
                {                                                       
                    result = left_op * right_op.get_offset_element(0,0); 
                                                                        
                    break;                                              
                }                                                       
                                                                        
                default:                                                
                {                                                       
                    L_THROW(0);                                         
                                                                        
                    break;                                              
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case E_VS:                                                      
        {                                                               
            switch ( left_op.get_v_type() )                             
            {                                                           
                case NORMAL_VECTOR:                                     
                case ZERO_VECTOR:                                       
                case SCALAR_VECTOR:                                     
                case SELECT_VECTOR:                                     
                {                                                       
                    result = left_op * right_op.get_offset_element(0,0); 
                                                                        
                    break;                                              
                }                                                       
                                                                        
                default:                                                
                {                                                       
                    L_THROW(0);                                         
                                                                        
                    break;                                              
                }                                                       
            }                                                           
                                                                        
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





// *= multiplicative assignment - binary, return lvalue
// /= divisive       assignment - binary, return lvalue

fMATRIX  &operator*=(      fMATRIX  &left_op, const   double &right_op)
{
    FN_ENTRY_POINT

    UNARY_MATRIX_OPERATION_VARIA(,*=right_op);

    FN_EXIT_POINT left_op;
}

fMATRIX  &operator*=(      fMATRIX  &left_op, const c_double &right_op)
{
    FN_ENTRY_POINT

    UNARY_MATRIX_OPERATION_VARIA(,*=right_op);

    FN_EXIT_POINT left_op;
}

fMATRIX  &assign_all(      fMATRIX  &left_op, const   double &right_op)
{
    FN_ENTRY_POINT

    UNARY_MATRIX_OPERATION_VARIA(,=right_op);

    FN_EXIT_POINT left_op;
}

fMATRIX  &assign_all(      fMATRIX  &left_op, const c_double &right_op)
{
    FN_ENTRY_POINT

    UNARY_MATRIX_OPERATION_VARIA(,=right_op);

    FN_EXIT_POINT left_op;
}


fMATRIX  &operator/=(      fMATRIX  &left_op, const   double &right_op)
{
    FN_ENTRY_POINT

    UNARY_MATRIX_OPERATION_VARIA(,/=right_op);

    FN_EXIT_POINT left_op;
}

fMATRIX  &operator/=(      fMATRIX  &left_op, const c_double &right_op)
{
    FN_ENTRY_POINT

    UNARY_MATRIX_OPERATION_VARIA(,/=right_op);

    FN_EXIT_POINT left_op;
}




// *= multiplicative assignment - binary, return lvalue

fMATRIX  &operator*=(      fMATRIX  &left_op, const fMATRIX  &right_op)
{
    FN_ENTRY_POINT

    fMATRIX temp;

    temp = left_op * right_op;

    left_op = temp;

    FN_EXIT_POINT left_op;
}

fVECTOR  &operator*=(      fVECTOR  &left_op, const fMATRIX  &right_op)
{
    FN_ENTRY_POINT

    fVECTOR temp;

    temp = left_op * right_op;

    left_op = temp;

    FN_EXIT_POINT left_op;
}




//
// Relational operator overloading
//

// == equivalence

int operator==(const fMATRIX  &left_op, const fMATRIX  &right_op)
{
    FN_ENTRY_POINT

    int result;

    MATRIX_MATRIX_RELATION(&&,==,1);

    FN_EXIT_POINT result;
}






// != inequivalence

int operator!=(const fMATRIX  &left_op, const fMATRIX  &right_op)
{
    FN_ENTRY_POINT

    int result;

    MATRIX_MATRIX_RELATION(||,!=,0);

    FN_EXIT_POINT result;
}




// < left than

int operator< (const fMATRIX  &left_op, const fMATRIX  &right_op)
{
    FN_ENTRY_POINT

    int result;

    MATRIX_MATRIX_RELATION(&&,< ,1);

    FN_EXIT_POINT result;
}




// <= less than or equal to

int operator<=(const fMATRIX  &left_op, const fMATRIX  &right_op)
{
    FN_ENTRY_POINT

    int result;

    MATRIX_MATRIX_RELATION(&&,<=,1);

    FN_EXIT_POINT result;
}




// > greater than

int operator> (const fMATRIX  &left_op, const fMATRIX  &right_op)
{
    FN_ENTRY_POINT

    int result;

    MATRIX_MATRIX_RELATION(&&,> ,1);

    FN_EXIT_POINT result;
}




// >= greater than or equal to

int operator>= (const fMATRIX  &left_op, const fMATRIX  &right_op)
{
    FN_ENTRY_POINT

    int result;

    MATRIX_MATRIX_RELATION(&&,>= ,1);

    FN_EXIT_POINT result;
}













//
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=                                                                +=+=
// +=+=                      MATRIX CLASSES                            +=+=
// +=+=                                                                +=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
//


//
// Constructor and destructor:
//

fMATRIX::fMATRIX()
{
    FN_ENTRY_POINT

    where_to   = NULL;
    echo_level = 0;

    m_diag = 0.0;
    m_type = ZERO_MATRIX;
    size   = ZERO_MATRIX_SIZE;

    reset_offsets();

    all_it = NULL;
    matrix_struct_lookup.make_normal(DO_FORCE);

    FN_EXIT_POINT;
}

fMATRIX::fMATRIX(char, long _size, int _m_type)
{
    FN_ENTRY_POINT

    long i;

    where_to   = NULL;
    echo_level = 0;

    m_diag = 0.0;
    m_type = _m_type;
    size   = 0;

    reset_offsets();

    all_it = NULL;
    matrix_struct_lookup.make_normal(DO_FORCE);

    switch ( m_type )
    {
        case ASYMMETRIC_MATRIX:
        case SYMMETRIC_MATRIX:
        case LOWER_TRIANGULAR_MATRIX:
        case UPPER_TRIANGULAR_MATRIX:
        case DIAGONAL_MATRIX:
        {
            THROW_ASSERT(_size >= 0);

            reset_offsets();

            if ( _size > 0 )
            {
                for ( i = 1 ; i <= _size ; i++ )
                {
                    addend();
                }
            }

            break;
        }

        case SCALAR_MATRIX:
        {
            size   = SCALAR_MATRIX_SIZE;
            m_diag = 1.0;

            break;
        }

        case ZERO_MATRIX:
        {
            size = ZERO_MATRIX_SIZE;

            break;
        }

        default:
        {
            L_THROW(0);

            break;
        }
    }

    FN_EXIT_POINT;
}

fMATRIX::fMATRIX(const fMATRIX &source)
{
    FN_ENTRY_POINT

    l_xMatrix *naut;

    long i,j;
    long start_row;                                                     
    long start_col;                                                     
    long end_row;                                                       
    long end_col;                                                       
    long rsize;                                                         
                                                                        
    where_to   = NULL;                                                    
    echo_level = 0;                                                     
                                                                        
    m_diag = 0.0;                                                  
    size   = source.size;                                               
    m_type = source.m_type;                                             
                                                                        
    start_row = ( offset_start_effective_row = source.offset_start_effective_row ); 
    start_col = ( offset_start_effective_col = source.offset_start_effective_col ); 
    end_row   = ( offset_end_effective_row   = source.offset_end_effective_row   ); 
    end_col   = ( offset_end_effective_col   = source.offset_end_effective_col   ); 
                                                                        
    m_diag = source.m_diag;                                             
                                                                        
    all_it = NULL;                                                      
    matrix_struct_lookup.make_normal(DO_FORCE);
                                                                        
    if ( size > 0 )                                                     
    {                                                                   
        rsize = size;                                                   
        size  = 0;                                                      
                                                                        
        reset_offsets();                                                

        switch ( m_type )                                               
        {                                                               
            case ASYMMETRIC_MATRIX:                                     
            {                                                           
                for ( i = 1 ; i <= rsize ; i++ )                        
                {                                                       
                    addend();                                           
                                                                        
                    naut = GET_XMATRIX_N(i);                            
                                                                        
                    (naut->diag) = source.getelm(i,i);                  
                                                                        
                    ((naut->row)[i-1]) = (naut->diag);                  
                    ((naut->col)[i-1]) = (naut->diag);                  
                                                                        
                    if ( i > 1 )                                        
                    {                                                   
                        for ( j = 1 ; j < i ; j++ )                     
                        {                                               
                            ((naut->row)[j-1]) = source.getelm(i,j);    
                            ((naut->col)[j-1]) = source.getelm(j,i);    
                        }                                               
                    }                                                   
                }                                                       
                                                                        
                break;                                                  
            }                                                           
                                                                        
            case SYMMETRIC_MATRIX:                                      
            case LOWER_TRIANGULAR_MATRIX:                               
            {                                                           
                for ( i = 1 ; i <= rsize ; i++ )                        
                {                                                       
                    addend();                                           
                                                                        
                    naut = GET_XMATRIX_N(i);                            
                                                                        
                    (naut->diag) = source.getelm(i,i);                  
                                                                        
                    ((naut->row)[i-1]) = (naut->diag);                  
                                                                        
                    if ( i > 1 )                                        
                    {                                                   
                        for ( j = 1 ; j < i ; j++ )                     
                        {                                               
                            ((naut->row)[j-1]) = source.getelm(i,j);    
                        }                                               
                    }                                                   
                }                                                       
                                                                        
                break;                                                  
            }                                                           
                                                                        
            case UPPER_TRIANGULAR_MATRIX:                               
            {                                                           
                for ( i = 1 ; i <= rsize ; i++ )                        
                {                                                       
                    addend();                                           
                                                                        
                    naut = GET_XMATRIX_N(i);                            
                                                                        
                    (naut->diag) = source.getelm(i,i);                  
                                                                        
                    ((naut->col)[i-1]) = (naut->diag);                  
                                                                        
                    if ( i > 1 )                                        
                    {                                                   
                        for ( j = 1 ; j < i ; j++ )                     
                        {                                               
                            ((naut->col)[j-1]) = source.getelm(j,i);    
                        }                                               
                    }                                                   
                }                                                       
                                                                        
                break;                                                  
            }                                                           
                                                                        
            case DIAGONAL_MATRIX:                                       
            {                                                           
                for ( i = 1 ; i <= rsize ; i++ )                        
                {                                                       
                    addend();                                           
                                                                        
                    naut = GET_XMATRIX_N(i);                            
                                                                        
                    (naut->diag) = source.getelm(i,i);                  
                }                                                       
                                                                        
                break;                                                  
            }                                                           
                                                                        
            case SCALAR_MATRIX:                                         
            case ZERO_MATRIX:                                           
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
                                                                        
        offset_start_effective_row = start_row;                         
        offset_start_effective_col = start_col;                         
        offset_end_effective_row   = end_row;                           
        offset_end_effective_col   = end_col;
    }                                                                   
                                                                        
    FN_EXIT_POINT;
}

fMATRIX::fMATRIX(const fVECTOR &source)
{
    FN_ENTRY_POINT

    vector_constructor_subfunction(source);

    FN_EXIT_POINT;
}

fMATRIX::fMATRIX(const   double &x)
{
    FN_ENTRY_POINT

    fVECTOR source('x',SCALAR_VECTOR_SIZE);

    source[-1] = x;

    vector_constructor_subfunction(source);

    FN_EXIT_POINT;
}

fMATRIX::fMATRIX(const c_double &x)
{
    FN_ENTRY_POINT

    fVECTOR source('x',SCALAR_VECTOR_SIZE);

    source[-1] = x;

    vector_constructor_subfunction(source);

    FN_EXIT_POINT;
}


void fMATRIX::vector_constructor_subfunction(const fVECTOR &source)
{
    FN_ENTRY_POINT

    long i;                                                             
    long rsize;                                                         
                                                                        
    where_to   = NULL;                                                    
    echo_level = 0;                                                     
                                                                        
    m_diag = 0.0;

    all_it = NULL;
    matrix_struct_lookup.make_normal(DO_FORCE);

    switch ( source.get_v_type() )                                      
    {                                                                   
        case NORMAL_VECTOR:                                             
        {                                                               
            size   = source.get_effective_size();                       
            m_type = DIAGONAL_MATRIX;                                   
                                                                        
            rsize = size;                                               
            size  = 0;                                                  
                                                                        
            reset_offsets();                                            
                                                                        
            if ( rsize > 0 )                                            
            {                                                           
                for ( i = 1 ; i <= rsize ; i++ )                        
                {                                                       
                    addend(source.get_offset_element(i-1));             
                }                                                       
            }                                                           
                                                                        
            reset_offsets();                                            
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case ZERO_VECTOR:                                               
        {                                                               
            size   = ZERO_MATRIX_SIZE;                                  
            m_type = ZERO_MATRIX;                                       
                                                                        
            reset_offsets();                                            
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case SCALAR_VECTOR:                                             
        {                                                               
            size   = SCALAR_MATRIX_SIZE;                                
            m_type = SCALAR_MATRIX;                                     
                                                                        
            m_diag = source.vector_m_diag;                              
                                                                        
            reset_offsets();                                            
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case SELECT_VECTOR:                                             
        {                                                               
            L_THROW(1);                                                 
                                                                        
            break;                                                      
        }                                                               
                                                                        
        default:                                                        
        {                                                               
            L_THROW(2);                                                 
                                                                        
            break;                                                      
        }                                                               
    }                                                                   

    FN_EXIT_POINT;
}

fMATRIX::~fMATRIX()
{
    FN_ENTRY_POINT

    l_xMatrix *Xnow;
    l_xMatrix *Xthen;

    long i;                                                             
                                                                        
    reset_offsets();                                                    
                                                                        
    if ( size <= 0 ) { FN_EXIT_POINT; }                                        
                                                                        
    Xnow = all_it;                                                      
                                                                        
    for ( i = 1 ; i <= size ; i++ )                                     
    {                                                                   
        Xthen = Xnow->next;                                             
                                                                        
        switch ( m_type )                                               
        {                                                               
            case ASYMMETRIC_MATRIX:                                     
            {                                                           
                REPDELB(Xnow->row);                                     
                REPDELB(Xnow->col);

                Xnow->row = NULL;
                Xnow->col = NULL;
                                                                        
                break;                                                  
            }                                                           
                                                                        
            case SYMMETRIC_MATRIX:                                      
            case LOWER_TRIANGULAR_MATRIX:                               
            {                                                           
                REPDELB(Xnow->row);

                Xnow->row = NULL;
                                                                        
                break;                                                  
            }                                                           
                                                                        
            case UPPER_TRIANGULAR_MATRIX:                               
            {                                                           
                REPDELB(Xnow->col);

                Xnow->col = NULL;
                                                                        
                break;                                                  
            }                                                           
                                                                        
            case DIAGONAL_MATRIX:                                       
            {                                                           
                break;                                                  
            }                                                           
                                                                        
            case SCALAR_MATRIX:                                         
            case ZERO_MATRIX:                                           
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
                                                                        
        REPDEL(Xnow);

        Xnow = Xthen;                                                   
    }                                                                   

    FN_EXIT_POINT;
}




//
// Overloaded assignment for matrices
//

fMATRIX &fMATRIX::operator=(const fMATRIX &left_op)
{
    FN_ENTRY_POINT

    long left_size_row;
    long left_size_col;
    long size_row;
    long size_col;
    long offset;
    long left_offset;
    long left_size_max;

    if ( ((void *) this) != ((void *) (&left_op)) )
    {
        {
            fix_symm();

            left_size_row = left_op.get_effective_height();
            left_size_col = left_op.get_effective_width();

            // If size is zero, make it consistent

            if ( ( left_size_row == 0 ) || ( left_size_col == 0 ) )
            {
                left_size_row = 0;
                left_size_col = 0;
            }

            left_offset = left_op.get_e_offs();

            // If we have an empty, zero or scalar matrix, we need to
            // re-size it to the source size (assuming it to be square).

            switch ( size )
            {
                case 0:
                case SCALAR_MATRIX_SIZE:
                case ZERO_MATRIX_SIZE:
                {
                    // We need to make a SQUARE matrix big enought to hold
                    // the source matrix (don't worry about the possibility
                    // that the source is zero or scalar - this will be
                    // taken care of in any case).

                    if ( left_size_row > left_size_col )
                    {
                        left_size_max = left_size_row;
                    }

                    else
                    {
                        left_size_max = left_size_col;
                    }

                    // Fix the symmetry of the matrix based on the source

                    switch ( left_op.get_e_type() )
                    {
                        case E_EMPTY:
                        {
                            make_diagonal(DO_FORCE);

                            reset_offsets();

                            break;
                        }

                        case E_ZERO:
                        {
                            make_diagonal(DO_FORCE);

                            pad_matrix(left_size_max);

                            reset_offsets();

                            // allow for non-squareness

                            set_offset_start(1,1);
                            set_offset_end(left_size_row,left_size_col);

                            break;
                        }

                        case E_DIAG:
                        {
                            if ( left_offset == 0 )
                            {
                                make_diagonal(DO_FORCE);
                            }

                            else
                            {
                                if ( left_offset > 0 )
                                {
                                    make_lower_triangular(DO_FORCE);
                                }

                                else
                                {
                                    make_upper_triangular(DO_FORCE);
                                }
                            }

                            pad_matrix(left_size_max);

                            reset_offsets();

                            set_offset_start(1,1);
                            set_offset_end(left_size_row,left_size_col);

                            break;
                        }

                        case E_LOWER:                                   
                        {                                               
                            if ( left_offset >= 0 )                     
                            {                                           
                                make_lower_triangular(DO_FORCE);        
                            }                                           
                                                                        
                            else                                        
                            {                                           
                                make_asymmetric(DO_FORCE);              
                            }                                           
                                                                        
                            pad_matrix(left_size_max);                  
                                                                        
                            reset_offsets();                            
                                                                        
                            set_offset_start(1,1);                      
                            set_offset_end(left_size_row,left_size_col); 
                                                                        
                            break;                                      
                        }                                               
                                                                        
                        case E_UPPER:                                   
                        {                                               
                            if ( left_offset <= 0 )                     
                            {                                           
                                make_upper_triangular(DO_FORCE);        
                            }                                           
                                                                        
                            else                                        
                            {                                           
                                make_asymmetric(DO_FORCE);              
                            }                                           
                                                                        
                            pad_matrix(left_size_max);                  
                                                                        
                            reset_offsets();                            
                                                                        
                            set_offset_start(1,1);                      
                            set_offset_end(left_size_row,left_size_col); 
                                                                        
                            break;                                      
                        }                                               
                                                                        
                        case E_SYMM:                                    
                        {                                               
                            make_symmetric(DO_FORCE);                   
                                                                        
                            pad_matrix(left_size_max);                  
                                                                        
                            reset_offsets();                            
                                                                        
                            set_offset_start(1,1);                      
                            set_offset_end(left_size_row,left_size_col); 
                                                                        
                            break;                                      
                        }                                               
                                                                        
                        case E_ASYMM:                                   
                        {                                               
                            if ( left_size_row > 1 )                    
                            {                                           
                                make_asymmetric(DO_FORCE);              
                            }                                           
                                                                        
                            pad_matrix(left_size_max);                  
                                                                        
                            reset_offsets();                            
                                                                        
                            set_offset_start(1,1);                      
                            set_offset_end(left_size_row,left_size_col); 
                                                                        
                            break;                                      
                        }                                               
                                                                        
                        case E_VZ:                                      
                        {                                               
                            make_zero(DO_FORCE);                        
                                                                        
                            m_diag = 0.0;                          
                                                                        
                            break;                                      
                        }                                               
                                                                        
                        case E_VS:                                      
                        {                                               
                            make_scalar(DO_FORCE);                      
                                                                        
                            m_diag = left_op.m_diag;                    
                                                                        
                            break;                                      
                        }                                               
                                                                        
                        default:                                        
                        {                                               
                            L_THROW(0);                                 
                                                                        
                            break;                                      
                        }                                               
                    }                                                   
                                                                        
                    break;                                              
                }                                                       
                                                                        
                default:                                                
                {
                    THROW_ASSERT(size >= 0);
                                                                        
                    break;                                              
                }                                                       
            }                                                           
                                                                        
            size_row = get_effective_height();                          
            size_col = get_effective_width();                           
                                                                        
            if ( ( size_row == 0 ) || ( size_col == 0 ) )               
            {                                                           
                size_row = 0;                                           
                size_col = 0;                                           
            }                                                           
                                                                        
            offset = get_e_offs();

            // Fix the symmetry if required, and also ensure that
            // things match up correctly.

            // There is no need to be TOO pedantic here - the overwrites
            // through the scalar references will deal with all cases
            // except for the conversion of symmetric matrices to
            // assymetric ones.

            #ifndef NDEBUG
            {
                if ( ( left_op.get_e_type() != E_VZ ) &&
                     ( left_op.get_e_type() != E_VS )    )
                {                                                           
                    THROW_ASSERT(( size_row == left_size_row ) && ( size_col == left_size_col ));
                }
            }
            #endif
                                                                        
            switch ( get_e_type() )                                     
            {                                                           
                case E_EMPTY:                                           
                {                                                       
                    switch ( left_op.get_e_type() )                     
                    {                                                   
                        case E_EMPTY:                                   
                        {                                               
                            break;                                      
                        }                                               
                                                                        
                        case E_ZERO:                                    
                        case E_DIAG:                                    
                        case E_LOWER:                                   
                        case E_UPPER:                                   
                        case E_SYMM:                                    
                        case E_ASYMM:                                   
                        {                                               
                            L_THROW(0);                                 
                                                                        
                            break;                                      
                        }                                               
                                                                        
                        case E_VZ:                                      
                        case E_VS:                                      
                        {                                               
                            break;                                      
                        }                                               
                                                                        
                        default:                                        
                        {                                               
                            L_THROW(0);                                 
                                                                        
                            break;                                      
                        }                                               
                    }                                                   
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case E_ZERO:                                            
                case E_DIAG:                                            
                case E_LOWER:                                           
                case E_UPPER:                                           
                case E_ASYMM:                                           
                {                                                       
                    switch ( left_op.get_e_type() )                     
                    {                                                   
                        case E_EMPTY:                                   
                        {                                               
                            L_THROW(0);                                 
                                                                        
                            break;                                      
                        }                                               
                                                                        
                        case E_ZERO:                                    
                        case E_DIAG:                                    
                        case E_LOWER:                                   
                        case E_UPPER:                                   
                        case E_SYMM:                                    
                        case E_ASYMM:                                   
                        case E_VZ:                                      
                        case E_VS:                                      
                        {                                               
                            break;                                      
                        }                                               
                                                                        
                        default:                                        
                        {                                               
                            L_THROW(0);                                 
                                                                        
                            break;                                      
                        }                                               
                    }                                                   
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case E_SYMM:                                            
                {                                                       
                    switch ( left_op.get_e_type() )                     
                    {                                                   
                        case E_EMPTY:                                   
                        {                                               
                            L_THROW(0);                                 
                                                                        
                            break;                                      
                        }                                               
                                                                        
                        case E_ZERO:                                    
                        case E_DIAG:                                    
                        case E_LOWER:                                   
                        case E_UPPER:                                   
                        case E_ASYMM:                                   
                        case E_VZ:                                      
                        case E_VS:                                      
                        {                                               
                            make_asymmetric();                          
                                                                        
                            break;                                      
                        }                                               
                                                                        
                        case E_SYMM:                                    
                        {                                               
                            break;                                      
                        }                                               
                                                                        
                        default:                                        
                        {                                               
                            L_THROW(0);                                 
                                                                        
                            break;                                      
                        }                                               
                    }                                                   
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case E_VZ:                                              
                {                                                       
                    switch ( left_op.get_e_type() )                     
                    {                                                   
                        case E_EMPTY:                                   
                        case E_ZERO:                                    
                        case E_DIAG:                                    
                        case E_LOWER:                                   
                        case E_UPPER:                                   
                        case E_SYMM:                                    
                        case E_ASYMM:                                   
                        {                                               
                            L_THROW(0);                                 
                                                                        
                            break;                                      
                        }                                               
                                                                        
                        case E_VZ:                                      
                        {                                               
                            break;                                      
                        }                                               
                                                                        
                        case E_VS:                                      
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
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case E_VS:                                              
                {                                                       
                    switch ( left_op.get_e_type() )                     
                    {                                                   
                        case E_EMPTY:                                   
                        case E_ZERO:                                    
                        case E_DIAG:                                    
                        case E_LOWER:                                   
                        case E_UPPER:                                   
                        case E_SYMM:                                    
                        case E_ASYMM:                                   
                        {                                               
                            L_THROW(0);                                 
                                                                        
                            break;                                      
                        }                                               
                                                                        
                        case E_VZ:                                      
                        {                                               
                            L_THROW(0);                                 
                                                                        
                            break;                                      
                        }                                               
                                                                        
                        case E_VS:                                      
                        {                                               
                            break;                                      
                        }                                               
                                                                        
                        default:                                        
                        {                                               
                            L_THROW(0);                                 
                                                                        
                            break;                                      
                        }                                               
                    }                                                   
                                                                        
                    break;                                              
                }                                                       
                                                                        
                default:                                                
                {                                                       
                    L_THROW(0);                                         
                                                                        
                    break;                                              
                }                                                       
            }

            // Finally, actually do the copy.

            assign_all(*this,0.0);

            (*this) += left_op;
        }
    }

    FN_EXIT_POINT (*this);
}

fMATRIX &fMATRIX::operator=(const fVECTOR &left_op)
{
    FN_ENTRY_POINT

    fMATRIX temp(left_op);

    (*this) = temp;

    FN_EXIT_POINT (*this);
}

fMATRIX &fMATRIX::operator=(const   double &left_op)
{
    FN_ENTRY_POINT

    fMATRIX temp(SCALAR_MATRIX_SIZE,SCALAR_MATRIX);

    temp.m_diag = left_op;

    (*this) = temp;

    FN_EXIT_POINT (*this);
}

fMATRIX &fMATRIX::operator=(const c_double &left_op)
{
    FN_ENTRY_POINT

    fMATRIX temp(SCALAR_MATRIX_SIZE,SCALAR_MATRIX);

    temp.m_diag = left_op;

    (*this) = temp;

    FN_EXIT_POINT (*this);
}




//
// Assignment from indexed matrix.
//

fMATRIX &fMATRIX::assign_from_index(const fMATRIX &left_op, const lVECTOR &b)
{
    FN_ENTRY_POINT

    long left_size;
    long left_offset;
    long result_size;
    long i,j;
    L_DOUBLE temp;

    if ( ((void *) this) != ((void *) (&left_op)) )
    {
        {
            THROW_ASSERT(!is_not_eff_square());

            fix_symm();

            left_size   = b.get_effective_size();
            left_offset = left_op.get_e_offs();

            // If we have an empty, zero or scalar matrix, we need to
            // re-size it to the source size (assuming it to be square).

            switch ( size )
            {
                case 0:
                case SCALAR_MATRIX_SIZE:
                case ZERO_MATRIX_SIZE:
                {
                    // Fix the symmetry of the matrix based on the source

                    switch ( left_op.get_e_type() )
                    {
                        case E_EMPTY:
                        {
                            make_diagonal(DO_FORCE);

                            reset_offsets();

                            break;
                        }

                        case E_ZERO:
                        {
                            make_diagonal(DO_FORCE);

                            pad_matrix(left_size);

                            reset_offsets();

                            break;
                        }

                        case E_DIAG:
                        {
                            if ( left_offset == 0 )
                            {
                                make_diagonal(DO_FORCE);
                            }

                            else
                            {
                                if ( left_offset > 0 )
                                {
                                    make_lower_triangular(DO_FORCE);
                                }

                                else
                                {
                                    make_upper_triangular(DO_FORCE);
                                }
                            }

                            pad_matrix(left_size);

                            reset_offsets();

                            break;
                        }

                        case E_LOWER:                                   
                        {                                               
                            if ( left_offset >= 0 )                     
                            {                                           
                                make_lower_triangular(DO_FORCE);        
                            }                                           
                                                                        
                            else                                        
                            {                                           
                                make_asymmetric(DO_FORCE);              
                            }                                           
                                                                        
                            pad_matrix(left_size);
                                                                        
                            reset_offsets();                            
                                                                        
                            break;                                      
                        }                                               
                                                                        
                        case E_UPPER:                                   
                        {                                               
                            if ( left_offset <= 0 )                     
                            {                                           
                                make_upper_triangular(DO_FORCE);        
                            }                                           
                                                                        
                            else                                        
                            {                                           
                                make_asymmetric(DO_FORCE);              
                            }                                           
                                                                        
                            pad_matrix(left_size);                  
                                                                        
                            reset_offsets();                            
                                                                        
                            break;                                      
                        }                                               
                                                                        
                        case E_SYMM:                                    
                        {                                               
                            make_symmetric(DO_FORCE);                   
                                                                        
                            pad_matrix(left_size);                  
                                                                        
                            reset_offsets();                            
                                                                        
                            break;                                      
                        }                                               
                                                                        
                        case E_ASYMM:                                   
                        {                                               
                            if ( left_size > 1 )                    
                            {                                           
                                make_asymmetric(DO_FORCE);              
                            }                                           
                                                                        
                            pad_matrix(left_size);                  
                                                                        
                            reset_offsets();                            
                                                                        
                            break;                                      
                        }                                               
                                                                        
                        case E_VZ:                                      
                        {                                               
                            make_zero(DO_FORCE);                        
                                                                        
                            m_diag = 0.0;                          
                                                                        
                            break;                                      
                        }                                               
                                                                        
                        case E_VS:                                      
                        {                                               
                            make_scalar(DO_FORCE);                      
                                                                        
                            m_diag = left_op.m_diag;                    
                                                                        
                            break;                                      
                        }                                               
                                                                        
                        default:                                        
                        {                                               
                            L_THROW(0);                                 
                                                                        
                            break;                                      
                        }                                               
                    }                                                   
                                                                        
                    break;                                              
                }                                                       
                                                                        
                default:                                                
                {
                    THROW_ASSERT(size > 0);
                                                                        
                    break;                                              
                }                                                       
            }                                                           
                                                                        
            result_size = get_effective_height();                          

            // Fix the symmetry if required, and also ensure that
            // things match up correctly.

            THROW_ASSERT( ( left_op.get_e_type() == E_VZ ) || ( left_op.get_e_type() == E_VS ) || ( result_size == left_size ) );
                                                                        
            switch ( get_e_type() )                                     
            {                                                           
                case E_EMPTY:                                           
                {                                                       
                    switch ( left_op.get_e_type() )                     
                    {                                                   
                        case E_EMPTY:                                   
                        {                                               
                            break;                                      
                        }                                               
                                                                        
                        case E_ZERO:
                        case E_DIAG:                                    
                        case E_LOWER:                                   
                        case E_UPPER:                                   
                        case E_SYMM:                                    
                        case E_ASYMM:                                   
                        {                                               
                            L_THROW(0);                                 
                                                                        
                            break;                                      
                        }

                        case E_VZ:                                      
                        case E_VS:
                        {
                            break;
                        }
                                                                        
                        default:                                        
                        {                                               
                            L_THROW(0);                                 
                                                                        
                            break;                                      
                        }                                               
                    }                                                   
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case E_ZERO:                                            
                case E_DIAG:                                            
                case E_LOWER:                                           
                case E_UPPER:                                           
                case E_ASYMM:                                           
                {
                    switch ( left_op.get_e_type() )                     
                    {                                                   
                        case E_EMPTY:                                   
                        {                                               
                            L_THROW(0);                                 
                                                                        
                            break;                                      
                        }                                               
                                                                        
                        case E_ZERO:
                        case E_DIAG:
                        case E_LOWER:                                   
                        case E_UPPER:                                   
                        case E_SYMM:                                    
                        case E_ASYMM:                                   
                        {
                            THROW_ASSERT(left_size == result_size);

                            for ( i = 1 ; i <= result_size ; i++ )
                            {
                                for ( j = 1 ; j <= result_size ; j++ )
                                {
                                    (*this)(i-1,j-1) = left_op.get_offset_element((b.get_offset_element(i-1))-1,(b.get_offset_element(j-1))-1);
                                }
                            }

                            break;
                        }

                        case E_VZ:
                        {
                            for ( i = 1 ; i <= result_size ; i++ )
                            {
                                for ( j = 1 ; j <= result_size ; j++ )
                                {
                                    (*this)(i-1,j-1) = 0.0;
                                }
                            }

                            break;
                        }

                        case E_VS:                                      
                        {
                            temp = left_op.get_offset_element(-1,-1);

                            for ( i = 1 ; i <= result_size ; i++ )
                            {
                                for ( j = 1 ; j <= result_size ; j++ )
                                {
                                    (*this)(i-1,j-1) = 0.0;
                                }

                                (*this)(i-1,i-1) = temp;
                            }

                            break;
                        }

                        default:                                        
                        {                                               
                            L_THROW(0);                                 
                                                                        
                            break;                                      
                        }                                               
                    }                                                   
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case E_SYMM:                                            
                {
                    switch ( left_op.get_e_type() )                     
                    {                                                   
                        case E_EMPTY:                                   
                        {                                               
                            L_THROW(0);                                 
                                                                        
                            break;                                      
                        }                                               
                                                                        
                        case E_ZERO:
                        case E_DIAG:
                        case E_LOWER:                                   
                        case E_UPPER:                                   
                        case E_ASYMM:                                   
                        {
                            THROW_ASSERT(left_size == result_size);

                            make_asymmetric();

                            for ( i = 1 ; i <= result_size ; i++ )
                            {
                                for ( j = 1 ; j <= result_size ; j++ )
                                {
                                    (*this)(i-1,j-1) = left_op.get_offset_element((b.get_offset_element(i-1))-1,(b.get_offset_element(j-1))-1);
                                }
                            }

                            break;
                        }

                        case E_SYMM:                                    
                        {
                            THROW_ASSERT(left_size == result_size);

                            for ( i = 1 ; i <= result_size ; i++ )
                            {
                                for ( j = 1 ; j <= i ; j++ )
                                {
                                    (*this)(i-1,j-1) = left_op.get_offset_element((b.get_offset_element(i-1))-1,(b.get_offset_element(j-1))-1);
                                }
                            }

                            break;
                        }

                        case E_VZ:
                        {
                            for ( i = 1 ; i <= result_size ; i++ )
                            {
                                for ( j = 1 ; j <= result_size ; j++ )
                                {
                                    (*this)(i-1,j-1) = 0.0;
                                }
                            }

                            break;
                        }

                        case E_VS:                                      
                        {
                            temp = left_op.get_offset_element(-1,-1);

                            for ( i = 1 ; i <= result_size ; i++ )
                            {
                                for ( j = 1 ; j <= result_size ; j++ )
                                {
                                    (*this)(i-1,j-1) = 0.0;
                                }

                                (*this)(i-1,i-1) = temp;
                            }

                            break;
                        }

                        default:                                        
                        {                                               
                            L_THROW(0);                                 
                                                                        
                            break;                                      
                        }                                               
                    }                                                   
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case E_VZ:                                              
                {                                                       
                    switch ( left_op.get_e_type() )                     
                    {                                                   
                        case E_EMPTY:                                   
                        case E_ZERO:                                    
                        case E_DIAG:                                    
                        case E_LOWER:                                   
                        case E_UPPER:                                   
                        case E_SYMM:                                    
                        case E_ASYMM:                                   
                        {                                               
                            L_THROW(0);                                 
                                                                        
                            break;                                      
                        }                                               
                                                                        
                        case E_VZ:                                      
                        {                                               
                            break;                                      
                        }                                               
                                                                        
                        case E_VS:                                      
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
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case E_VS:                                              
                {                                                       
                    switch ( left_op.get_e_type() )                     
                    {                                                   
                        case E_EMPTY:                                   
                        case E_ZERO:                                    
                        case E_DIAG:                                    
                        case E_LOWER:                                   
                        case E_UPPER:                                   
                        case E_SYMM:                                    
                        case E_ASYMM:                                   
                        {                                               
                            L_THROW(0);                                 
                                                                        
                            break;                                      
                        }                                               
                                                                        
                        case E_VZ:                                      
                        {                                               
                            L_THROW(0);                                 
                                                                        
                            break;                                      
                        }                                               
                                                                        
                        case E_VS:                                      
                        {
                            (*this)(-1,-1) = left_op.get_offset_element(-1,-1);

                            break;                                      
                        }                                               
                                                                        
                        default:                                        
                        {                                               
                            L_THROW(0);                                 
                                                                        
                            break;                                      
                        }                                               
                    }                                                   
                                                                        
                    break;                                              
                }                                                       
                                                                        
                default:                                                
                {                                                       
                    L_THROW(0);                                         
                                                                        
                    break;                                              
                }                                                       
            }
        }
    }

    FN_EXIT_POINT (*this);
}




//
// Set and get offsets - or reset to defaults (whole matrix, but does
// not effect attributes).
//

void fMATRIX::reset_offsets(void)
{
    FN_ENTRY_POINT

    offset_start_effective_row = 1;
    offset_start_effective_col = 1;

    offset_end_effective_row = size;
    offset_end_effective_col = size;

    FN_EXIT_POINT;
}

void fMATRIX::set_offset_start_at_row(long start_row)
{
    FN_ENTRY_POINT

    set_offsets(start_row,get_offset_start_at_col(),get_offset_end_at_row(),get_offset_end_at_col());

    FN_EXIT_POINT;
}

void fMATRIX::set_offset_start_at_col(long start_col)
{
    FN_ENTRY_POINT

    set_offsets(get_offset_start_at_row(),start_col,get_offset_end_at_row(),get_offset_end_at_col());

    FN_EXIT_POINT;
}

void fMATRIX::set_offset_end_at_row(long end_row)
{
    FN_ENTRY_POINT

    set_offsets(get_offset_start_at_row(),get_offset_start_at_col(),end_row,get_offset_end_at_col());

    FN_EXIT_POINT;
}

void fMATRIX::set_offset_end_at_col(long end_col)
{
    FN_ENTRY_POINT

    set_offsets(get_offset_start_at_row(),get_offset_start_at_col(),get_offset_end_at_row(),end_col);

    FN_EXIT_POINT;
}

void fMATRIX::set_offset_start(long start_row, long start_col)
{
    FN_ENTRY_POINT

    set_offsets(start_row,start_col,get_offset_end_at_row(),get_offset_end_at_col());

    FN_EXIT_POINT;
}

void fMATRIX::set_offset_end(long end_row, long end_col)
{
    FN_ENTRY_POINT

    set_offsets(get_offset_start_at_row(),get_offset_start_at_col(),end_row,end_col);

    FN_EXIT_POINT;
}

void fMATRIX::set_offsets(long start_row, long start_col, long end_row, long end_col)
{
    FN_ENTRY_POINT

    switch ( m_type )                                                   
    {                                                                   
        case ASYMMETRIC_MATRIX:                                         
        case SYMMETRIC_MATRIX:                                          
        case LOWER_TRIANGULAR_MATRIX:                                   
        case UPPER_TRIANGULAR_MATRIX:                                   
        case DIAGONAL_MATRIX:                                           
        {                                                               
            THROW_ASSERT(start_row >= 1);                                     
            THROW_ASSERT(start_row <= (end_row+1));                           
            THROW_ASSERT(end_row >= (start_row-1));                           
            THROW_ASSERT(end_row <= size);                                    
            THROW_ASSERT(start_col >= 1);                                     
            THROW_ASSERT(start_col <= (end_col+1));                           
            THROW_ASSERT(end_col >= (start_col-1));                           
            THROW_ASSERT(end_col <= size);                                    
                                                                        
            offset_start_effective_row = start_row;                     
            offset_start_effective_col = start_col;                     
                                                                        
            offset_end_effective_row = end_row;                         
            offset_end_effective_col = end_col;                         
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case SCALAR_MATRIX:                                             
        case ZERO_MATRIX:                                               
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

    FN_EXIT_POINT;
}

long fMATRIX::get_offset_start_at_row(void) const
{
    FN_ENTRY_POINT

    switch ( m_type )                                                   
    {                                                                   
        case ASYMMETRIC_MATRIX:                                         
        case SYMMETRIC_MATRIX:                                          
        case LOWER_TRIANGULAR_MATRIX:                                   
        case UPPER_TRIANGULAR_MATRIX:                                   
        case DIAGONAL_MATRIX:                                           
        {                                                               
            FN_EXIT_POINT offset_start_effective_row;                              
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case SCALAR_MATRIX:                                             
        case ZERO_MATRIX:                                               
        {                                                               
            FN_EXIT_POINT 1;                                                   
                                                                        
            break;                                                      
        }                                                               
                                                                        
        default:                                                        
        {                                                               
            L_THROW(0);                                                 
                                                                        
            break;                                                      
        }                                                               
    }                                                                   

    // unreachable code

    FN_EXIT_POINT 1;
}

long fMATRIX::get_offset_start_at_col(void) const
{
    FN_ENTRY_POINT

    switch ( m_type )                                                   
    {                                                                   
        case ASYMMETRIC_MATRIX:                                         
        case SYMMETRIC_MATRIX:                                          
        case LOWER_TRIANGULAR_MATRIX:                                   
        case UPPER_TRIANGULAR_MATRIX:                                   
        case DIAGONAL_MATRIX:                                           
        {                                                               
            FN_EXIT_POINT offset_start_effective_col;
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case SCALAR_MATRIX:                                             
        case ZERO_MATRIX:                                               
        {                                                               
            FN_EXIT_POINT 1;                                                   
                                                                        
            break;                                                      
        }                                                               
                                                                        
        default:                                                        
        {                                                               
            L_THROW(0);                                                 
                                                                        
            break;                                                      
        }                                                               
    }                                                                   

    // unreachable code

    FN_EXIT_POINT 1;
}

long fMATRIX::get_offset_end_at_row(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT offset_end_effective_row;
}

long fMATRIX::get_offset_end_at_col(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT offset_end_effective_col;
}

long fMATRIX::get_effective_height(void) const
{
    FN_ENTRY_POINT

    switch ( size )
    {
        case SCALAR_MATRIX_SIZE:
        case ZERO_MATRIX_SIZE:
        {
            FN_EXIT_POINT size;

            break;
        }

        default:
        {
            break;
        }
    }

    FN_EXIT_POINT offset_end_effective_row-offset_start_effective_row+1;
}

long fMATRIX::get_effective_width(void) const
{
    FN_ENTRY_POINT

    switch ( size )
    {
        case SCALAR_MATRIX_SIZE:
        case ZERO_MATRIX_SIZE:
        {
            FN_EXIT_POINT size;

            break;
        }

        default:
        {
            break;
        }
    }

    FN_EXIT_POINT offset_end_effective_col-offset_start_effective_col+1;
}

long fMATRIX::get_real_size(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT size;
}




//
// Modify matrix symmetry attribuates - not that not all types can be
// converted to just any other type.
//
// force: 0 - normal.  Throw an error is change will cause non-zero
//            elements to go to zero.
//        1 - force elements off symmetry to zero if necessary
//

void fMATRIX::make_asymmetric(int force)
{
    FN_ENTRY_POINT

    l_xMatrix *naut;

    long i,j;

    switch ( m_type )                                                   
    {                                                                   
        case ASYMMETRIC_MATRIX:                                         
        {                                                               
            break;                                                      
        }                                                               
                                                                        
        case SYMMETRIC_MATRIX:                                          
        {                                                               
            if ( size > 0 )                                             
            {                                                           
                naut = all_it;                                          
                                                                        
                for ( i = 1 ; i <= size ; i++ )                         
                {                                                       
                    REPNEWB((naut->col),L_DOUBLE,i);                   
                                                                        
                    if ( i > 1 )                                        
                    {                                                   
                        for ( j = 1 ; j < i ; j++ )                     
                        {                                               
                            (naut->col)[j-1] = (naut->row)[j-1];        
                        }                                               
                    }                                                   
                                                                        
                    (naut->col)[i-1] = naut->diag;                      
                                                                        
                    naut = naut->next;                                  
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case LOWER_TRIANGULAR_MATRIX:                                   
        {                                                               
            if ( size > 0 )                                             
            {                                                           
                naut = all_it;                                          
                                                                        
                for ( i = 1 ; i <= size ; i++ )                         
                {                                                       
                    REPNEWB((naut->col),L_DOUBLE,i);                   
                                                                        
                    if ( i > 1 )                                        
                    {                                                   
                        for ( j = 1 ; j < i ; j++ )                     
                        {                                               
                            (naut->col)[j-1] = 0.0;                
                        }                                               
                    }                                                   
                                                                        
                    (naut->col)[i-1] = naut->diag;                      
                                                                        
                    naut = naut->next;                                  
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case UPPER_TRIANGULAR_MATRIX:                                   
        {                                                               
            if ( size > 0 )                                             
            {                                                           
                naut = all_it;                                          
                                                                        
                for ( i = 1 ; i <= size ; i++ )                         
                {                                                       
                    REPNEWB((naut->row),L_DOUBLE,i);                   
                                                                        
                    if ( i > 1 )                                        
                    {                                                   
                        for ( j = 1 ; j < i ; j++ )                     
                        {                                               
                            (naut->row)[j-1] = 0.0;                
                        }                                               
                    }                                                   
                                                                        
                    (naut->row)[i-1] = naut->diag;                      
                                                                        
                    naut = naut->next;                                  
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case DIAGONAL_MATRIX:                                           
        {                                                               
            if ( size > 0 )                                             
            {                                                           
                naut = all_it;                                          
                                                                        
                for ( i = 1 ; i <= size ; i++ )                         
                {                                                       
                    REPNEWB((naut->row),L_DOUBLE,i);                   
                    REPNEWB((naut->col),L_DOUBLE,i);                   
                                                                        
                    if ( i > 1 )                                        
                    {                                                   
                        for ( j = 1 ; j < i ; j++ )                     
                        {                                               
                            (naut->row)[j-1] = 0.0;                
                            (naut->col)[j-1] = 0.0;                
                        }                                               
                    }                                                   
                                                                        
                    (naut->row)[i-1] = naut->diag;                      
                    (naut->col)[i-1] = naut->diag;                      
                                                                        
                    naut = naut->next;                                  
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case SCALAR_MATRIX:                                             
        case ZERO_MATRIX:                                               
        {
            THROW_ASSERT(force);

            size   = 0;
            m_type = ASYMMETRIC_MATRIX;                             
                                                                        
            all_it = NULL;                                          
                                                                        
            reset_offsets();                                        
                                                                        
            break;                                                      
        }                                                               
                                                                        
        default:                                                        
        {                                                               
            L_THROW(0);                                                 
                                                                        
            break;                                                      
        }                                                               
    }                                                                   
                                                                        
    m_type = ASYMMETRIC_MATRIX;                                         

    FN_EXIT_POINT;

    force = 0;
}

void fMATRIX::make_symmetric(int force)
{
    FN_ENTRY_POINT

    l_xMatrix *naut;

    long i,j;                                                           
    int throw_flag = 0;                                                 
                                                                        
    switch ( m_type )                                                   
    {                                                                   
        case ASYMMETRIC_MATRIX:                                         
        {                                                               
            if ( size > 0 )                                             
            {                                                           
                naut = all_it;                                          
                                                                        
                for ( i = 1 ; i <= size ; i++ )                         
                {                                                       
                    if ( ( i > 1 ) && ( force == 0 ) )                  
                    {                                                   
                        for ( j = 1 ; j < i ; j++ )                     
                        {                                               
                            if ( (naut->row)[j-1] != (naut->col)[j-1] ) 
                            {                                           
                                throw_flag = 1;                         
                            }                                           
                        }                                               
                    }                                                   
                                                                        
                    REPDELB(naut->col);

                    naut->col = naut->row;                              
                                                                        
                    naut = naut->next;                                  
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case SYMMETRIC_MATRIX:                                          
        {                                                               
            break;                                                      
        }                                                               
                                                                        
        case LOWER_TRIANGULAR_MATRIX:                                   
        {                                                               
            if ( size > 0 )                                             
            {                                                           
                naut = all_it;                                          
                                                                        
                for ( i = 1 ; i <= size ; i++ )                         
                {                                                       
                    if ( ( i > 1 ) && ( force == 0 ) )                  
                    {                                                   
                        for ( j = 1 ; j < i ; j++ )                     
                        {                                               
                            if ( (naut->row)[j-1] != 0.0 )         
                            {                                           
                                throw_flag = 1;                         
                            }                                           
                        }                                               
                    }                                                   
                                                                        
                    naut->col = naut->row;                              
                                                                        
                    naut = naut->next;                                  
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case UPPER_TRIANGULAR_MATRIX:                                   
        {                                                               
            if ( size > 0 )                                             
            {                                                           
                naut = all_it;                                          
                                                                        
                for ( i = 1 ; i <= size ; i++ )                         
                {                                                       
                    if ( ( i > 1 ) && ( force == 0 ) )                  
                    {                                                   
                        for ( j = 1 ; j < i ; j++ )                     
                        {                                               
                            if ( (naut->col)[j-1] != 0.0 )         
                            {                                           
                                throw_flag = 1;                         
                            }                                           
                        }                                               
                    }                                                   
                                                                        
                    naut->row = naut->col;                              
                                                                        
                    naut = naut->next;                                  
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case DIAGONAL_MATRIX:                                           
        {                                                               
            if ( size > 0 )                                             
            {                                                           
                naut = all_it;                                          
                                                                        
                for ( i = 1 ; i <= size ; i++ )                         
                {                                                       
                    REPNEWB((naut->row),L_DOUBLE,i);                   
                                                                        
                    naut->col = naut->row;                              
                                                                        
                    if ( i > 1 )                                        
                    {                                                   
                        for ( j = 1 ; j < i ; j++ )                     
                        {                                               
                            (naut->row)[j-1] = 0.0;                
                        }                                               
                    }                                                   
                                                                        
                    (naut->row)[i-1] = naut->diag;                      
                                                                        
                    naut = naut->next;                                  
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case SCALAR_MATRIX:                                             
        case ZERO_MATRIX:                                               
        {
            THROW_ASSERT(force);

            size   = 0;
            m_type = SYMMETRIC_MATRIX;                              

            all_it = NULL;
                                                                        
            reset_offsets();                                        
                                                                        
            break;                                                      
        }                                                               
                                                                        
        default:                                                        
        {                                                               
            L_THROW(0);                                                 
                                                                        
            break;                                                      
        }                                                               
    }                                                                   
                                                                        
    m_type = SYMMETRIC_MATRIX;

    THROW_ASSERT(!throw_flag);

    FN_EXIT_POINT;
}

void fMATRIX::make_lower_triangular(int force)
{
    FN_ENTRY_POINT

    l_xMatrix *naut;

    long i,j;                                                           
    int throw_flag = 0;                                                 
                                                                        
    switch ( m_type )                                                   
    {                                                                   
        case ASYMMETRIC_MATRIX:                                         
        {                                                               
            if ( size > 0 )                                             
            {                                                           
                naut = all_it;                                          
                                                                        
                for ( i = 1 ; i <= size ; i++ )                         
                {                                                       
                    if ( ( i > 1 ) && ( force == 0 ) )                  
                    {                                                   
                        for ( j = 1 ; j < i ; j++ )                     
                        {                                               
                            if ( (naut->col)[j-1] != 0.0 )         
                            {                                           
                                throw_flag = 1;                         
                            }                                           
                        }                                               
                    }                                                   
                                                                        
                    REPDELB(naut->col);                                 
                                                                        
                    naut->col = NULL;                                   
                                                                        
                    naut = naut->next;                                  
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case SYMMETRIC_MATRIX:                                          
        {                                                               
            if ( size > 0 )                                             
            {                                                           
                naut = all_it;                                          
                                                                        
                for ( i = 1 ; i <= size ; i++ )                         
                {                                                       
                    if ( ( i > 1 ) && ( force == 0 ) )                  
                    {                                                   
                        for ( j = 1 ; j < i ; j++ )                     
                        {                                               
                            if ( (naut->col)[j-1] != 0.0 )         
                            {                                           
                                throw_flag = 1;                         
                            }                                           
                        }                                               
                    }                                                   
                                                                        
                    naut->col = NULL;                                   
                                                                        
                    naut = naut->next;                                  
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case LOWER_TRIANGULAR_MATRIX:                                   
        {                                                               
            break;                                                      
        }                                                               
                                                                        
        case UPPER_TRIANGULAR_MATRIX:                                   
        {                                                               
            if ( size > 0 )                                             
            {                                                           
                naut = all_it;                                          
                                                                        
                for ( i = 1 ; i <= size ; i++ )                         
                {                                                       
                    if ( ( i > 1 ) && ( force == 0 ) )                  
                    {                                                   
                        for ( j = 1 ; j < i ; j++ )                     
                        {                                               
                            if ( (naut->col)[j-1] != 0.0 )         
                            {                                           
                                throw_flag = 1;                         
                            }                                           
                        }                                               
                    }                                                   
                                                                        
                    naut->row = naut->col;                              
                    naut->col = NULL;                                   
                                                                        
                    naut = naut->next;                                  
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case DIAGONAL_MATRIX:                                           
        {                                                               
            if ( size > 0 )                                             
            {                                                           
                naut = all_it;                                          
                                                                        
                for ( i = 1 ; i <= size ; i++ )                         
                {                                                       
                    REPNEWB((naut->row),L_DOUBLE,i);                   
                                                                        
                    if ( i > 1 )                                        
                    {                                                   
                        for ( j = 1 ; j < i ; j++ )                     
                        {                                               
                            (naut->row)[j-1] = 0.0;                
                        }                                               
                    }                                                   
                                                                        
                    (naut->row)[i-1] = naut->diag;                      
                                                                        
                    naut = naut->next;                                  
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case SCALAR_MATRIX:                                             
        case ZERO_MATRIX:                                               
        {
            THROW_ASSERT(force);

            size   = 0;
            m_type = LOWER_TRIANGULAR_MATRIX;                       
                                                                        
            all_it = NULL;                                          
                                                                        
            reset_offsets();                                        
                                                                        
            break;                                                      
        }                                                               
                                                                        
        default:                                                        
        {                                                               
            L_THROW(0);                                                 
                                                                        
            break;                                                      
        }                                                               
    }                                                                   
                                                                        
    m_type = LOWER_TRIANGULAR_MATRIX;                                   
                                                                        
    THROW_ASSERT(!throw_flag);

    FN_EXIT_POINT;
}

void fMATRIX::make_upper_triangular(int force)
{
    FN_ENTRY_POINT

    l_xMatrix *naut;

    long i,j;                                                           
    int throw_flag = 0;                                                 
                                                                        
    switch ( m_type )                                                   
    {                                                                   
        case ASYMMETRIC_MATRIX:                                         
        {                                                               
            if ( size > 0 )                                             
            {                                                           
                naut = all_it;                                          
                                                                        
                for ( i = 1 ; i <= size ; i++ )                         
                {                                                       
                    if ( ( i > 1 ) && ( force == 0 ) )                  
                    {                                                   
                        for ( j = 1 ; j < i ; j++ )                     
                        {                                               
                            if ( (naut->row)[j-1] != 0.0 )         
                            {                                           
                                throw_flag = 1;                         
                            }                                           
                        }                                               
                    }                                                   
                                                                        
                    REPDELB(naut->row);                                 
                                                                        
                    naut->row = NULL;                                   
                                                                        
                    naut = naut->next;                                  
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case SYMMETRIC_MATRIX:                                          
        {                                                               
            if ( size > 0 )                                             
            {                                                           
                naut = all_it;                                          
                                                                        
                for ( i = 1 ; i <= size ; i++ )                         
                {                                                       
                    if ( ( i > 1 ) && ( force == 0 ) )                  
                    {                                                   
                        for ( j = 1 ; j < i ; j++ )                     
                        {                                               
                            if ( (naut->row)[j-1] != 0.0 )         
                            {                                           
                                throw_flag = 1;                         
                            }                                           
                        }                                               
                    }                                                   
                                                                        
                    naut->row = NULL;                                   
                                                                        
                    naut = naut->next;                                  
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case LOWER_TRIANGULAR_MATRIX:                                   
        {                                                               
            if ( size > 0 )                                             
            {                                                           
                naut = all_it;                                          
                                                                        
                for ( i = 1 ; i <= size ; i++ )                         
                {                                                       
                    if ( ( i > 1 ) && ( force == 0 ) )                  
                    {                                                   
                        for ( j = 1 ; j < i ; j++ )                     
                        {                                               
                            if ( (naut->row)[j-1] != 0.0 )         
                            {                                           
                                throw_flag = 1;                         
                            }                                           
                        }                                               
                    }                                                   
                                                                        
                    naut->col = naut->row;                              
                    naut->row = NULL;                                   
                                                                        
                    naut = naut->next;                                  
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case UPPER_TRIANGULAR_MATRIX:                                   
        {                                                               
            break;                                                      
        }                                                               
                                                                        
        case DIAGONAL_MATRIX:                                           
        {                                                               
            if ( size > 0 )                                             
            {                                                           
                naut = all_it;                                          
                                                                        
                for ( i = 1 ; i <= size ; i++ )                         
                {                                                       
                    REPNEWB((naut->col),L_DOUBLE,i);                   
                                                                        
                    if ( i > 1 )                                        
                    {                                                   
                        for ( j = 1 ; j < i ; j++ )                     
                        {                                               
                            (naut->col)[j-1] = 0.0;                
                        }                                               
                    }                                                   
                                                                        
                    (naut->col)[i-1] = naut->diag;                      
                                                                        
                    naut = naut->next;                                  
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case SCALAR_MATRIX:                                             
        case ZERO_MATRIX:                                               
        {                                                               
            THROW_ASSERT(force);

            size   = 0;                                             
            m_type = UPPER_TRIANGULAR_MATRIX;                       
                                                                        
            all_it = NULL;                                          
                                                                        
            reset_offsets();                                        
                                                                        
            break;                                                      
        }                                                               
                                                                        
        default:                                                        
        {                                                               
            L_THROW(0);                                                 
                                                                        
            break;                                                      
        }                                                               
    }                                                                   
                                                                        
    m_type = UPPER_TRIANGULAR_MATRIX;                                   
                                                                        
    THROW_ASSERT(!throw_flag);

    FN_EXIT_POINT;
}

void fMATRIX::make_diagonal(int force)
{
    FN_ENTRY_POINT

    l_xMatrix *naut;

    long i,j;                                                           
    int throw_flag = 0;                                                 
                                                                        
    switch ( m_type )                                                   
    {                                                                   
        case ASYMMETRIC_MATRIX:                                         
        {                                                               
            if ( size > 0 )                                             
            {                                                           
                naut = all_it;                                          
                                                                        
                for ( i = 1 ; i <= size ; i++ )                         
                {                                                       
                    if ( ( i > 1 ) && ( force == 0 ) )                  
                    {                                                   
                        for ( j = 1 ; j < i ; j++ )                     
                        {                                               
                            if ( ( (naut->row)[j-1] != 0.0 ) ||    
                                 ( (naut->col)[j-1] != 0.0 )    )  
                            {                                           
                                throw_flag = 1;                         
                            }                                           
                        }                                               
                    }                                                   
                                                                        
                    REPDELB(naut->row);                                 
                    REPDELB(naut->col);                                 
                                                                        
                    naut->row = NULL;                                   
                    naut->col = NULL;                                   
                                                                        
                    naut = naut->next;                                  
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case SYMMETRIC_MATRIX:                                          
        {                                                               
            if ( size > 0 )                                             
            {                                                           
                naut = all_it;                                          
                                                                        
                for ( i = 1 ; i <= size ; i++ )                         
                {                                                       
                    if ( ( i > 1 ) && ( force == 0 ) )                  
                    {                                                   
                        for ( j = 1 ; j < i ; j++ )                     
                        {                                               
                            if ( (naut->row)[j-1] != 0.0 )         
                            {                                           
                                throw_flag = 1;                         
                            }                                           
                        }                                               
                    }                                                   
                                                                        
                    REPDELB(naut->row);                                 
                                                                        
                    naut->row = NULL;                                   
                    naut->col = NULL;                                   
                                                                        
                    naut = naut->next;                                  
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case LOWER_TRIANGULAR_MATRIX:                                   
        {                                                               
            if ( size > 0 )                                             
            {                                                           
                naut = all_it;                                          
                                                                        
                for ( i = 1 ; i <= size ; i++ )                         
                {                                                       
                    if ( ( i > 1 ) && ( force == 0 ) )                  
                    {                                                   
                        for ( j = 1 ; j < i ; j++ )                     
                        {                                               
                            if ( (naut->row)[j-1] != 0.0 )         
                            {                                           
                                throw_flag = 1;                         
                            }                                           
                        }                                               
                    }                                                   
                                                                        
                    REPDELB(naut->row);                                 
                                                                        
                    naut->row = NULL;                                   
                                                                        
                    naut = naut->next;                                  
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case UPPER_TRIANGULAR_MATRIX:                                   
        {                                                               
            if ( size > 0 )                                             
            {                                                           
                naut = all_it;                                          
                                                                        
                for ( i = 1 ; i <= size ; i++ )                         
                {                                                       
                    if ( ( i > 1 ) && ( force == 0 ) )                  
                    {                                                   
                        for ( j = 1 ; j < i ; j++ )                     
                        {                                               
                            if ( (naut->col)[j-1] != 0.0 )         
                            {                                           
                                throw_flag = 1;                         
                            }                                           
                        }                                               
                    }                                                   
                                                                        
                    REPDELB(naut->col);                                 
                                                                        
                    naut->col = NULL;                                   
                                                                        
                    naut = naut->next;                                  
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case DIAGONAL_MATRIX:                                           
        {                                                               
            break;                                                      
        }                                                               
                                                                        
        case SCALAR_MATRIX:                                             
        case ZERO_MATRIX:                                               
        {                                                               
            THROW_ASSERT(force);

            size   = 0;
            m_type = DIAGONAL_MATRIX;                               
                                                                        
            all_it = NULL;                                          
                                                                        
            reset_offsets();                                        
                                                                        
            break;                                                      
        }                                                               
                                                                        
        default:                                                        
        {                                                               
            L_THROW(0);                                                 
                                                                        
            break;                                                      
        }                                                               
    }                                                                   
                                                                        
    m_type = DIAGONAL_MATRIX;                                           
                                                                        
    THROW_ASSERT(!throw_flag);

    FN_EXIT_POINT;
}

void fMATRIX::make_scalar(int force)
{
    FN_ENTRY_POINT

    long i;                                                             
                                                                        
    switch ( m_type )                                                   
    {                                                                   
        case ASYMMETRIC_MATRIX:                                         
        case SYMMETRIC_MATRIX:                                          
        case LOWER_TRIANGULAR_MATRIX:                                   
        case UPPER_TRIANGULAR_MATRIX:                                   
        case DIAGONAL_MATRIX:                                           
        {                                                               
            THROW_ASSERT(force);

            if ( size > 0 )
            {                                                       
                for ( i = size ; i >= 1 ; i-- )                     
                {                                                   
                    remove(i);                                      
                }                                                   
            }                                                       
                                                                        
            size   = SCALAR_MATRIX_SIZE;                            
            m_type = SCALAR_MATRIX;                                 
                                                                        
            all_it = NULL;                                          
                                                                        
            reset_offsets();                                        
        }                                                               
                                                                        
        case ZERO_MATRIX:                                               
        {                                                               
            THROW_ASSERT(force);

            size   = SCALAR_MATRIX_SIZE;
            m_type = SCALAR_MATRIX;                                 
                                                                        
            all_it = NULL;                                          
                                                                        
            reset_offsets();                                        
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case SCALAR_MATRIX:                                             
        {                                                               
            break;                                                      
        }                                                               
                                                                        
        default:                                                        
        {                                                               
            L_THROW(0);                                                 
                                                                        
            break;                                                      
        }                                                               
    }                                                                   
                                                                        
    m_diag = 1.0;                                                   
                                                                        
    m_type = SCALAR_MATRIX;                                             

    FN_EXIT_POINT;

    force = 0;
}

void fMATRIX::make_zero(int force)
{
    FN_ENTRY_POINT

    long i;                                                             
                                                                        
    switch ( m_type )                                                   
    {                                                                   
        case ASYMMETRIC_MATRIX:                                         
        case SYMMETRIC_MATRIX:                                          
        case LOWER_TRIANGULAR_MATRIX:                                   
        case UPPER_TRIANGULAR_MATRIX:                                   
        case DIAGONAL_MATRIX:                                           
        {                                                               
            THROW_ASSERT(force);

            if ( size > 0 )                                         
            {                                                       
                for ( i = size ; i >= 1 ; i-- )                     
                {                                                   
                    remove(i);                                      
                }                                                   
            }                                                       
                                                                        
            size   = ZERO_MATRIX_SIZE;                              
            m_type = ZERO_MATRIX;                                   
                                                                        
            all_it = NULL;                                          
                                                                        
            reset_offsets();                                        
        }                                                               
                                                                        
        case SCALAR_MATRIX:                                             
        {                                                               
            THROW_ASSERT(force);

            size   = ZERO_MATRIX_SIZE;                              
            m_type = ZERO_MATRIX;                                   
                                                                        
            all_it = NULL;                                          
                                                                        
            reset_offsets();                                        
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case ZERO_MATRIX:                                               
        {                                                               
            break;                                                      
        }                                                               
                                                                        
        default:                                                        
        {                                                               
            L_THROW(0);                                                 
                                                                        
            break;                                                      
        }                                                               
    }                                                                   
                                                                        
    m_diag = 0.0;                                                  
                                                                        
    m_type = ZERO_MATRIX;                                               

    FN_EXIT_POINT;

    force = 0;
}

int fMATRIX::is_asymmetric(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT ( m_type == ASYMMETRIC_MATRIX );
}

int fMATRIX::is_symmetric(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT ( m_type == SYMMETRIC_MATRIX );
}

int fMATRIX::is_lower_triangular(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT ( m_type == LOWER_TRIANGULAR_MATRIX );
}

int fMATRIX::is_upper_triangular(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT ( m_type == UPPER_TRIANGULAR_MATRIX );
}

int fMATRIX::is_diagonal(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT ( m_type == DIAGONAL_MATRIX );
}

int fMATRIX::is_scalar(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT ( m_type == SCALAR_MATRIX );
}

int fMATRIX::is_zero(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT ( m_type == ZERO_MATRIX );
}

int fMATRIX::get_m_type(void) const
{
    FN_ENTRY_POINT

    FN_EXIT_POINT m_type;
}





//
// Transposition:
//

void fMATRIX::transpose(void)
{
    FN_ENTRY_POINT

    l_xMatrix *naut;

    long i;                                                             
    L_DOUBLE *temp;                                                    
                                                                        
    switch ( m_type )                                                   
    {                                                                   
        case SYMMETRIC_MATRIX:
        {
            break;
        }

        case ASYMMETRIC_MATRIX:                                         
        case LOWER_TRIANGULAR_MATRIX:                                   
        case UPPER_TRIANGULAR_MATRIX:                                   
        case DIAGONAL_MATRIX:                                           
        {                                                               
            if ( size > 0 )                                             
            {                                                           
                naut = all_it;                                          
                                                                        
                for ( i = 1 ; i <= size ; i++ )                         
                {                                                       
                    if ( m_type != DIAGONAL_MATRIX )                    
                    {                                                   
                        temp      = naut->row;                          
                        naut->row = naut->col;                          
                        naut->col = temp;                               
                    }                                                   
                                                                        
                    naut = naut->next;                                  
                }                                                       
            }                                                           
                                                                        
            if ( m_type == LOWER_TRIANGULAR_MATRIX )                    
            {                                                           
                m_type = UPPER_TRIANGULAR_MATRIX;                       
            }                                                           
                                                                        
            else                                                        
            {                                                           
                if ( m_type == UPPER_TRIANGULAR_MATRIX )                
                {                                                       
                    m_type = LOWER_TRIANGULAR_MATRIX;                   
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case SCALAR_MATRIX:                                             
        case ZERO_MATRIX:                                               
        {                                                               
            break;                                                      
        }                                                               
                                                                        
        default:                                                        
        {                                                               
            L_THROW(0);                                                 
                                                                        
            break;                                                      
        }                                                               
    }                                                                   

    FN_EXIT_POINT;
}




//
// Negation:
//

void fMATRIX::negate_it(void)
{
    FN_ENTRY_POINT

    (*this) *= -1.0;

    FN_EXIT_POINT;
}




//
// Misc info functions:
//
// is_not_eff_square(): returns 1 is effective matrix is not square,
//                      0 otherwise.
//

int fMATRIX::is_not_eff_square(void) const
{
    FN_ENTRY_POINT

    if ( get_effective_width() != get_effective_height() )
    {
        switch ( m_type )
        {
            case ASYMMETRIC_MATRIX:
            case SYMMETRIC_MATRIX:
            case LOWER_TRIANGULAR_MATRIX:
            case UPPER_TRIANGULAR_MATRIX:
            case DIAGONAL_MATRIX:
            {
                FN_EXIT_POINT 1;

                break;
            }

            case SCALAR_MATRIX:
            case ZERO_MATRIX:
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

    FN_EXIT_POINT 0;
}




//
// Get standard (T **) version of matrix, assuming that this is
// possible in finite dimensions.
//

L_DOUBLE **fMATRIX::operator()(void) const
{
    FN_ENTRY_POINT

    L_DOUBLE **result = NULL;
    long i,j;

    switch ( get_m_type() )
    {
        case ASYMMETRIC_MATRIX:
        case SYMMETRIC_MATRIX:
        case LOWER_TRIANGULAR_MATRIX:
        case UPPER_TRIANGULAR_MATRIX:
        case DIAGONAL_MATRIX:
        {
            REPNEWB(result,L_DOUBLE *,get_effective_height());

            if ( get_effective_height() > 0 )
            {
                for ( i = 1 ; i <= get_effective_height() ; i++ )
                {
                    REPNEWB(result[i-1],L_DOUBLE,get_effective_width());

                    for ( j = 1 ; j <= get_effective_width() ; j++ )
                    {
                        result[i-1][j-1] = get_offset_element(i-1,j-1);
                    }
                }
            }

            break;
        }

        case SCALAR_MATRIX:
        case ZERO_MATRIX:
        {
            REPNEWB(result,L_DOUBLE *,1);
            REPNEWB(result[0],L_DOUBLE,1);

            result[0][0] = get_offset_element(0,0);

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
// Array element overloading:
//

fVECTOR fMATRIX::operator[](long i) const
{
    FN_ENTRY_POINT

    fVECTOR result;

    long j,k;                                                           
                                                                        
    switch ( m_type )                                                   
    {                                                                   
        case ASYMMETRIC_MATRIX:                                         
        case SYMMETRIC_MATRIX:                                          
        case LOWER_TRIANGULAR_MATRIX:                                   
        case UPPER_TRIANGULAR_MATRIX:                                   
        case DIAGONAL_MATRIX:                                           
        {                                                               
            result.make_normal(DO_FORCE);                               
                                                                        
            if ( i == -1 )                                              
            {                                                           
                k = get_effective_height();                             
                                                                        
                if ( get_effective_width() < k )                        
                {                                                       
                    k = get_effective_width();                          
                }                                                       
                                                                        
                if ( k > 0 )                                            
                {                                                       
                    for ( j = 1 ; j <= k ; j++ )                        
                    {                                                   
                        result.addend(get_offset_element(j-1,j-1));     
                    }                                                   
                }                                                       
            }                                                           
                                                                        
            else                                                        
            {                                                           
                THROW_ASSERT(i >= 0);                                         
                THROW_ASSERT(i <= (get_effective_height()-1));                
                                                                        
                k = get_effective_width();                              
                                                                        
                if ( k > 0 )                                            
                {                                                       
                    for ( j = 1 ; j <= k ; j++ )                        
                    {                                                   
                        result.addend(get_offset_element(i,j-1));       
                    }                                                   
                }                                                       
            }                                                           
                                                                        
            result.vector_m_diag = m_diag;                              
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case SCALAR_MATRIX:                                             
        {                                                               
            if ( i == -1 )                                              
            {                                                           
                result.make_scalar(DO_FORCE);                           
            }                                                           
                                                                        
            else                                                        
            {                                                           
                THROW_ASSERT(i >= 0);                                         
                                                                        
                result.make_select(DO_FORCE);                           
                                                                        
                result.set_sel_elm(i+1);                                
            }                                                           
                                                                        
            result.vector_m_diag = m_diag;                              
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case ZERO_MATRIX:                                               
        {                                                               
            result.make_zero(DO_FORCE);                                 
                                                                        
            result.vector_m_diag = m_diag;                              
                                                                        
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

fVECTOR fMATRIX::operator()(long i) const
{
    FN_ENTRY_POINT

    fVECTOR result;

    long j,k;                                                           
                                                                        
    switch ( m_type )                                                   
    {                                                                   
        case ASYMMETRIC_MATRIX:                                         
        case SYMMETRIC_MATRIX:                                          
        case LOWER_TRIANGULAR_MATRIX:                                   
        case UPPER_TRIANGULAR_MATRIX:                                   
        case DIAGONAL_MATRIX:                                           
        {                                                               
            result.make_normal(DO_FORCE);                               
                                                                        
            if ( i == -1 )                                              
            {                                                           
                k = get_effective_height();                             
                                                                        
                if ( get_effective_width() < k )                        
                {                                                       
                    k = get_effective_width();                          
                }                                                       
                                                                        
                if ( k > 0 )                                            
                {                                                       
                    for ( j = 1 ; j <= k ; j++ )                        
                    {                                                   
                        result.addend(get_offset_element(j-1,j-1));     
                    }                                                   
                }                                                       
            }                                                           
                                                                        
            else                                                        
            {                                                           
                THROW_ASSERT(i >= 0);                                         
                THROW_ASSERT(i <= (get_effective_width()-1));                 
                                                                        
                k = get_effective_height();                             
                                                                        
                if ( k > 0 )                                            
                {                                                       
                    for ( j = 1 ; j <= k ; j++ )                        
                    {                                                   
                        result.addend(get_offset_element(j-1,i));       
                    }                                                   
                }                                                       
            }                                                           
                                                                        
            result.vector_m_diag = m_diag;                              
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case SCALAR_MATRIX:                                             
        {                                                               
            if ( i == -1 )                                              
            {                                                           
                result.make_scalar(DO_FORCE);                           
            }                                                           
                                                                        
            else                                                        
            {                                                           
                THROW_ASSERT(i >= 0);                                         
                                                                        
                result.make_select(DO_FORCE);                           
                                                                        
                result.set_sel_elm(i+1);                                
            }                                                           
                                                                        
            result.vector_m_diag = m_diag;                              
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case ZERO_MATRIX:                                               
        {                                                               
            result.make_zero(DO_FORCE);                                 
                                                                        
            result.vector_m_diag = m_diag;                              
                                                                        
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

L_DOUBLE &fMATRIX::operator()(long i, long j)
{
    FN_ENTRY_POINT

    long ii;
    long jj;

    ii = i;
    jj = j;

    if ( ii >= 0 ) { ii += get_offset_start_at_row(); }
    if ( jj >= 0 ) { jj += get_offset_start_at_col(); }

    FN_EXIT_POINT getref(ii,jj,1);
}




//
// get_offset_elm(i-1,j-1): acts like this[i-1][j-1] (with all the
//                          offsetting and crap).
//

L_DOUBLE fMATRIX::get_offset_element(long i, long j) const
{
    FN_ENTRY_POINT

    long ii;
    long jj;

    ii = i;
    jj = j;

    if ( ii >= 0 ) { ii += get_offset_start_at_row(); }
    if ( jj >= 0 ) { jj += get_offset_start_at_col(); }

    FN_EXIT_POINT getelm(ii,jj);
}




//
// Misc properties function:
//
// ARG_DIAG_PROD - product of diagonals of effective matrix
// ARG_TRACE     - sum of diagonals of effective matrix
// ARG_DET       - determinant of effective matrix
//
// Note: these will throw an exception if the matrix is not square.
//

L_DOUBLE fMATRIX::operator()(Matr_operation what) const
{
    FN_ENTRY_POINT

    L_DOUBLE result;
    long eff_size;
    long i;
    long *j;
    L_DOUBLE temp;
    int levi_result;
    int repeat_flag = 1;
    int zero_flag = 1;

    result = 0.0;

    THROW_ASSERT( ( what == ARG_DIAG_PROD ) || ( what == ARG_TRACE ) || ( what == ARG_DET ) );

    switch ( m_type )
    {
        case ASYMMETRIC_MATRIX:
        case SYMMETRIC_MATRIX:
        case LOWER_TRIANGULAR_MATRIX:
        case UPPER_TRIANGULAR_MATRIX:
        case DIAGONAL_MATRIX:
        {
            eff_size = get_effective_height();

            if ( ( eff_size == 0 ) || ( get_effective_width() == 0 ) )
            {
                if ( what == ARG_DIAG_PROD )
                {
                    result = 1.0;

                    FN_EXIT_POINT result;
                }

                if ( what == ARG_TRACE )
                {
                    result = 0.0;

                    FN_EXIT_POINT result;
                }

                if ( what == ARG_DET )
                {
                    result = 1.0;

                    FN_EXIT_POINT result;
                }
            }

            else
            {
                THROW_ASSERT(!is_not_eff_square());

                if ( what == ARG_DIAG_PROD )
                {
                    result = get_offset_element(1-1,1-1);

                    if ( eff_size > 1 )
                    {
                        for ( i = 2 ; i <= eff_size ; i++ )
                        {
                            result *= get_offset_element(i-1,i-1);
                        }
                    }

                    FN_EXIT_POINT result;
                }

                if ( what == ARG_TRACE )
                {
                    result = get_offset_element(1-1,1-1);

                    if ( eff_size > 1 )
                    {
                        for ( i = 2 ; i <= eff_size ; i++ )
                        {
                            result += get_offset_element(i-1,i-1);
                        }
                    }

                    FN_EXIT_POINT result;
                }

                if ( what == ARG_DET )
                {
                    REPNEWB(j,long,eff_size);

                    result = 0.0;

                    for ( i = 1 ; i <= eff_size ; i++ )
                    {
                        j[i-1] = i;
                    }

                    while ( repeat_flag )
                    {
                        if ( ( levi_result = levi_civita(eff_size,j) ) != 0 )
                        {
                            temp = get_offset_element(1-1,j[1-1]-1);

                            if ( eff_size > 1 )
                            {
                                for ( i = 2 ; i <= eff_size ; i++ )
                                {
                                    temp *= get_offset_element(i-1,j[i-1]-1);
                                }
                            }

                            if ( levi_result == +1 )
                            {
                                if ( zero_flag )
                                {
                                    result  = temp;
                                }

                                else
                                {
                                    result += temp;
                                }
                            }

                            else
                            {
                                if ( zero_flag )
                                {
                                    result  = -temp;
                                }

                                else
                                {
                                    result -= temp;
                                }
                            }

                            zero_flag = 0;
                        }

                        repeat_flag = 0;

                        for ( i = eff_size ; i >= 1 ; i-- )
                        {
                            if ( repeat_flag == 0 )
                            {
                                j[i-1]++;

                                if ( j[i-1] <= eff_size )
                                {
                                    repeat_flag = 1;
                                }

                                else
                                {
                                    j[i-1] = 1;
                                }
                            }
                        }
                    }

                    REPDELB(j);

                    j = NULL;
                }
            }

            break;
        }

        case SCALAR_MATRIX:
        {
            if ( what == ARG_TRACE )
            {
                result = 0.0;
            }

            else
            {
                result = m_diag;
            }

            break;
        }

        case ZERO_MATRIX:
        {
            result = 0.0;

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
// Rank update operators:
//

void fMATRIX::operator()(Matr_operation what, const fVECTOR &x, const   double &c)
{
    FN_ENTRY_POINT

    c_double temp;

    temp = c;

    (*this)(what,x,temp);

    FN_EXIT_POINT;
}

void fMATRIX::operator()(Matr_operation what, const fVECTOR &x, const c_double &c)
{
    FN_ENTRY_POINT

    long i,j;
    long eff_size;
    L_DOUBLE temp;
    L_DOUBLE tempb;

    THROW_ASSERT(what == ARG_RANK_ONE);

    if ( ( x != 0.0 ) && !(x.is_select()) && ( c != 0.0 ) )
    {
        switch ( get_e_type() )
        {
            case E_EMPTY:
            {
                THROW_ASSERT(!is_not_eff_square());
                THROW_ASSERT(x.get_effective_size() <= 0);

                break;
            }

            case E_ZERO:
            case E_DIAG:
            case E_UPPER:
            case E_LOWER:
            case E_ASYMM:
            {
                asym_rank_one_up:

                fix_symm();

                eff_size = get_effective_height();

                THROW_ASSERT(!is_not_eff_square());
                THROW_ASSERT((x.get_effective_size() < 0) || (x.get_effective_size() == eff_size));

                for ( i = 1 ; i <= eff_size ; i++ )
                {
                    for ( j = 1 ; j <= eff_size ; j++ )
                    {
                        temp  = x.get_offset_element(i-1);
                        temp *= c;
                        temp *= x.get_offset_element(j-1);

                        (*this)(i-1,j-1) += temp;
                    }
                }

                break;
            }

            case E_SYMM:
            {
                fix_symm();

                if ( get_e_type() != E_SYMM )
                {
                    goto asym_rank_one_up;
                }

                eff_size = get_effective_height();

                THROW_ASSERT(!is_not_eff_square());
                THROW_ASSERT((x.get_effective_size() < 0) || (x.get_effective_size() == eff_size));

                for ( i = 1 ; i <= eff_size ; i++ )
                {
                    for ( j = 1 ; j <= i ; j++ )
                    {
                        temp  = x.get_offset_element(i-1);
                        temp *= c;
                        temp *= x.get_offset_element(j-1);

                        (*this)(i-1,j-1) += temp;
                    }
                }

                break;
            }

            case E_VS:
            case E_VZ:
            {
                L_THROW(13);

                break;
            }

            default:
            {
                L_THROW(14);

                break;
            }
        }
    }

    else
    {
        if ( ( x != 0.0 ) && ( c != 0.0 ) )
        {
            // the only possibility here is that x.is_select() == 1

            fix_symm();

            tempb = x.get_offset_element((x.get_sel_elm())-1);

            temp  = tempb;
            temp *= c;
            temp *= tempb;

            (*this)((x.get_sel_elm())-1,(x.get_sel_elm())-1) += temp;
        }
    }

    FN_EXIT_POINT;

    what = ARG_DIAG_PROD;
}




//
// Misc operator function
//

fRef_pair fMATRIX::operator()(Matr_operation what, char, char) const
{
    FN_ENTRY_POINT

    fRef_pair result;
    long i,j;
    long e_size_row;
    long e_size_col;
    L_DOUBLE temp;

    THROW_ASSERT( ( what == ARG__MIN ) || ( what == ARG__MAX ) );

    temp = 0.0;
    result.element_val = 0.0;

    switch ( m_type )
    {
        case ASYMMETRIC_MATRIX:
        case SYMMETRIC_MATRIX:
        case LOWER_TRIANGULAR_MATRIX:
        case UPPER_TRIANGULAR_MATRIX:
        case DIAGONAL_MATRIX:
        {
            e_size_row = get_effective_height();
            e_size_col = get_effective_width();

            result.element_val = 0.0;

            result.element_num     = -1;
            result.element_num_row = -1;

            if ( what == ARG__MAX )
            {
                if ( ( e_size_row > 0 ) && ( e_size_col > 0 ) )
                {
                    for ( i = 1 ; i <= e_size_row ; i++ )
                    {
                        for ( j = 1 ; j <= e_size_col ; j++ )
                        {
                            if ( ( i == 1 ) && ( j == 1 ) )
                            {
                                result.element_num     = j;
                                result.element_num_row = i;

                                result.element_val = get_offset_element(i-1,j-1);
                            }

                            else
                            {
                                if ( (result.element_val) < (temp=get_offset_element(i-1,j-1)) )
                                {
                                    result.element_num     = j;
                                    result.element_num_row = i;

                                    result.element_val = temp;
                                }
                            }
                        }
                    }
                }
            }

            else
            {
                if ( ( e_size_row > 0 ) && ( e_size_col > 0 ) )
                {
                    for ( i = 1 ; i <= e_size_row ; i++ )
                    {
                        for ( j = 1 ; j <= e_size_col ; j++ )
                        {
                            if ( ( i == 1 ) && ( j == 1 ) )
                            {
                                result.element_num     = j;
                                result.element_num_row = i;

                                result.element_val = get_offset_element(i-1,j-1);
                            }

                            else
                            {
                                if ( (result.element_val) > (temp=get_offset_element(i-1,j-1)) )
                                {
                                    result.element_num     = j;
                                    result.element_num_row = i;

                                    result.element_val = temp;
                                }
                            }
                        }
                    }
                }
            }

            break;
        }

        case ZERO_MATRIX:
        {
            result.element_num     = 1;
            result.element_num_row = 1;

            result.element_val = 0.0;

            break;
        }

        case SCALAR_MATRIX:
        {
            if ( what == ARG__MAX )
            {
                if ( m_diag > (result.element_val) )
                {
                    result.element_num     = 1;
                    result.element_num_row = 1;

                    result.element_val = m_diag;
                }

                else
                {
                    result.element_num     = 1;
                    result.element_num_row = 2;

                    result.element_val = 0.0;
                }
            }

            else
            {
                if ( m_diag < (result.element_val) )
                {
                    result.element_num     = 1;
                    result.element_num_row = 1;

                    result.element_val = m_diag;
                }

                else
                {
                    result.element_num     = 1;
                    result.element_num_row = 2;

                    result.element_val = 0.0;
                }
            }

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
// Inversion operators:
//

//
// Generic inversion macro
//
// Given bxx = an array copy of b, overwrite bxx with y such that Gy = b
// (y'G = b' if ARG_INVERT_RIGHT used), where G is the present matrix and y
// and b have definite size e_size >= 0.
//
// Requires: e_size >= 0
//           bxx = array of size e_size containing b
//           declared longs i,j,k
// 

//
// Notes: - diagonal matrices are only invertable if the diagonal values
//          actually lie on the diagonal
//        - y'G = b' is the same as G'y = b, to do this case we just deal
//          with transposed values.
//        - the triangular cases are dealt with using standard foward
//          elimination / backward substitution.
//        - in the general (square) case, we do the following:
//          - make a row-wise copy of G in an array, so that it can be
//            pivoted easily
//          - moving down the matrix, for each row
//            - switch around do the diagonal on this row is non-zero
//            - divide through so this diagonal is one (divide G_row and bxx)
//              (optimisation note - actually, only need divide to the right
//              of the diagonal)
//            - For all rows below the present one, subtract the relevant
//              multiple of the present row to ensure that all elements
//              below the present diagonal are zero
//              (optimisation note - once again, only actually do maths to
//              the right of the present diagonal)
//          - this leaves us with an upper triangular matrix with 1's on
//            the diagonal.
//            (optimisation note - actually, there will be rubbish both on
//            and below the diagonal, but this is to be ignored as it is not
//            needed).
//          - From here on, we just deal with it like a triangular matrix,
//            with the added knowledge that all diagonals are 1.
//          - When G is scalar or zero, the method is trivial.
//

#define INV_GENERIC                                                     \
{                                                                       \
    L_DOUBLE whew;                                                      \
                                                                        \
    switch ( get_e_type() )                                             \
    {                                                                   \
        case E_EMPTY:                                                   \
        {                                                               \
            L_THROW(0);                                                 \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_ZERO:                                                    \
        {                                                               \
            throw INV_ERR;                                              \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_DIAG:                                                    \
        {                                                               \
            if ( get_e_offs() != 0 )                                    \
            {                                                           \
                throw INV_ERR;                                          \
            }                                                           \
                                                                        \
            for ( i = 1 ; i <= e_size ; i++ )                           \
            {                                                           \
                if ( ( i > z_start ) && ( i <= e_size-z_end ) )         \
                {                                                       \
                    whew = get_offset_element(i-1,i-1);                 \
                                                                        \
                    NZ_ASSERT(whew,zt);                                 \
                                                                        \
                    bxx[i-1] /= whew;                                   \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_LOWER:                                                   \
        {                                                               \
            if ( get_e_offs() > 0 )                                     \
            {                                                           \
                throw INV_ERR;                                          \
            }                                                           \
                                                                        \
            if ( get_e_offs() < 0 )                                     \
            {                                                           \
                goto generic_inversion;                                 \
            }                                                           \
                                                                        \
            if ( inv_type == INV_IS_LEFT )                              \
            {                                                           \
                for ( i = 1 ; i <= e_size ; i++ )                       \
                {                                                       \
                    if ( i > z_start )                                  \
                    {                                                   \
                        if ( i > 1+z_start )                            \
                        {                                               \
                            for ( j = 1+z_start ; j <= i-1 ; j++ )      \
                            {                                           \
                                bxx[i-1] -= ( (bxx[j-1]) * get_offset_element(i-1,j-1) ); \
                            }                                           \
                        }                                               \
                                                                        \
                        whew = get_offset_element(i-1,i-1);             \
                                                                        \
                        NZ_ASSERT(whew,zt);                             \
                                                                        \
                        bxx[i-1] /= whew;                               \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            else                                                        \
            {                                                           \
                for ( i = e_size ; i >= 1 ; i-- )                       \
                {                                                       \
                    if ( i <= e_size-z_end )                            \
                    {                                                   \
                        if ( i < e_size-z_end )                         \
                        {                                               \
                            for ( j=i+1 ; j <= e_size-z_end ; j++ )     \
                            {                                           \
                                bxx[i-1] -= ( (bxx[j-1]) * get_offset_element(j-1,i-1) ); \
                            }                                           \
                        }                                               \
                                                                        \
                        whew = get_offset_element(i-1,i-1);             \
                                                                        \
                        NZ_ASSERT(whew,zt);                             \
                                                                        \
                        bxx[i-1] /= whew;                               \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_UPPER:                                                   \
        {                                                               \
            if ( get_e_offs() < 0 )                                     \
            {                                                           \
                throw INV_ERR;                                          \
            }                                                           \
                                                                        \
            if ( get_e_offs() > 0 )                                     \
            {                                                           \
                goto generic_inversion;                                 \
            }                                                           \
                                                                        \
            if ( inv_type == INV_IS_LEFT )                              \
            {                                                           \
                for ( i = e_size ; i >= 1 ; i-- )                       \
                {                                                       \
                    if ( i <= e_size-z_end )                            \
                    {                                                   \
                        if ( i < e_size-z_end )                         \
                        {                                               \
                            for ( j=i+1 ; j <= e_size-z_end ; j++ )     \
                            {                                           \
                                bxx[i-1] -= ( (bxx[j-1]) * get_offset_element(i-1,j-1) ); \
                            }                                           \
                        }                                               \
                                                                        \
                        whew = get_offset_element(i-1,i-1);             \
                                                                        \
                        NZ_ASSERT(whew,zt);                             \
                                                                        \
                        bxx[i-1] /= whew;                               \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            else                                                        \
            {                                                           \
                for ( i = 1 ; i <= e_size ; i++ )                       \
                {                                                       \
                    if ( i > z_start )                                  \
                    {                                                   \
                        if ( i > 1+z_start )                            \
                        {                                               \
                            for ( j = 1+z_start ; j <= i-1 ; j++ )      \
                            {                                           \
                                bxx[i-1] -= ( (bxx[j-1]) * get_offset_element(j-1,i-1) ); \
                            }                                           \
                        }                                               \
                                                                        \
                        whew = get_offset_element(i-1,i-1);             \
                                                                        \
                        NZ_ASSERT(whew,zt);                             \
                                                                        \
                        bxx[i-1] /= whew;                               \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_SYMM:                                                    \
        case E_ASYMM:                                                   \
        {                                                               \
            generic_inversion:                                          \
                                                                        \
            L_DOUBLE **axx;                                             \
                                                                        \
            REPNEWB(axx,L_DOUBLE *,e_size);                             \
                                                                        \
            for ( i = 1 ; i <= e_size ; i++ )                           \
            {                                                           \
                REPNEWB(axx[i-1],L_DOUBLE,e_size);                      \
            }                                                           \
                                                                        \
            if ( inv_type == INV_IS_LEFT )                              \
            {                                                           \
                for ( i = 1 ; i <= e_size ; i++ )                       \
                {                                                       \
                    for ( j = 1 ; j <= e_size ; j++ )                   \
                    {                                                   \
                        axx[i-1][j-1] = get_offset_element(i-1,j-1);    \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            else                                                        \
            {                                                           \
                for ( i = 1 ; i <= e_size ; i++ )                       \
                {                                                       \
                    for ( j = 1 ; j <= e_size ; j++ )                   \
                    {                                                   \
                        axx[i-1][j-1] = get_offset_element(j-1,i-1);    \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            for ( i = 1 ; i <= e_size ; i++ )                           \
            {                                                           \
                j = i;                                                  \
                                                                        \
                rep_point_a_left:                                       \
                                                                        \
                try                                                     \
                {                                                       \
                    NZ_ASSERT(axx[j-1][i-1],zt);                        \
                }                                                       \
                                                                        \
                catch ( int err_num )                                   \
                {                                                       \
                    switch ( err_num )                                  \
                    {                                                   \
                        case INV_ERR:                                   \
                        {                                               \
                            if ( j < e_size )                           \
                            {                                           \
                                j++;                                    \
                                                                        \
                                goto rep_point_a_left;                  \
                            }                                           \
                        }                                               \
                                                                        \
                        default:                                        \
                        {                                               \
                            throw err_num;                              \
                                                                        \
                            break;                                      \
                        }                                               \
                    }                                                   \
                }                                                       \
                                                                        \
                if ( i != j )                                           \
                {                                                       \
                    tempa    = axx[i-1];                                \
                    axx[i-1] = axx[j-1];                                \
                    axx[j-1] = tempa;                                   \
                                                                        \
                    tempb    = bxx[i-1];                                \
                    bxx[i-1] = bxx[j-1];                                \
                    bxx[j-1] = tempb;                                   \
                }                                                       \
                                                                        \
                bxx[i-1] /= axx[i-1][i-1];                              \
                                                                        \
                if ( i < e_size )                                       \
                {                                                       \
                    for ( k = i+1 ; k <= e_size ; k++ )                 \
                    {                                                   \
                        axx[i-1][k-1] /= axx[i-1][i-1];                 \
                    }                                                   \
                                                                        \
                    for ( k = i+1 ; k <= e_size ; k++ )                 \
                    {                                                   \
                        bxx[k-1] -= ( (bxx[i-1]) * axx[k-1][i-1] );     \
                                                                        \
                        for ( l = i+1 ; l <= e_size ; l++ )             \
                        {                                               \
                            axx[k-1][l-1] -= ( axx[i-1][l-1] * axx[k-1][i-1] ); \
                        }                                               \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            if ( e_size > 1 )                                           \
            {                                                           \
                for ( i = e_size-1 ; i >= 1 ; i-- )                     \
                {                                                       \
                    for ( j = i+1 ; j <= e_size ; j++ )                 \
                    {                                                   \
                        bxx[i-1] -= ( bxx[j-1] * axx[i-1][j-1] );       \
                    }                                                   \
                }                                                       \
            }                                                           \
                                                                        \
            for ( i = 1 ; i <= e_size ; i++ )                           \
            {                                                           \
                REPDELB(axx[i-1]);                                      \
                                                                        \
                axx[i-1] = NULL;                                        \
            }                                                           \
                                                                        \
            REPDELB(axx);                                               \
                                                                        \
            axx = NULL;                                                 \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_VS:                                                      \
        {                                                               \
            NZ_ASSERT(m_diag,zt);                                       \
                                                                        \
            for ( i = 1 ; i <= e_size ; i++ )                           \
            {                                                           \
                bxx[i-1] /= m_diag;                                     \
            }                                                           \
                                                                        \
            break;                                                      \
        }                                                               \
                                                                        \
        case E_VZ:                                                      \
        {                                                               \
            throw INV_ERR;                                              \
                                                                        \
            break;                                                      \
        }                                                               \
    }                                                                   \
}


void fMATRIX::operator()(Matr_operation what, fVECTOR &y, const fVECTOR &b, long z_start, long z_end, double zt) const
{
    FN_ENTRY_POINT

    long i,j,k,l;
    long e_size;
    L_DOUBLE *bxx;
    int inv_type;
    long b_size;
    long y_size;

    if ( what == ARG_INVERT_LEFT )
    {
        inv_type = INV_IS_LEFT;
    }

    else
    {
        THROW_ASSERT(what == ARG_INVERT_RIGHT);

        inv_type = INV_IS_RIGHT;
    }

    b_size = b.get_effective_size();
    y_size = y.get_effective_size();
    e_size = get_effective_height();

    // check validity of arguments

    #ifndef NDEBUG
    switch ( get_e_type() )
    {
        case E_EMPTY:
        {
            THROW_ASSERT(b_size <= 0);
            THROW_ASSERT(y_size == 0);

            break;
        }

        case E_ZERO:
        case E_DIAG:
        case E_LOWER:
        case E_UPPER:
        case E_SYMM:
        case E_ASYMM:
        {
            THROW_ASSERT(!is_not_eff_square());
            THROW_ASSERT((b_size < 0) || (e_size == b_size));
            THROW_ASSERT(e_size == y_size);

            break;
        }

        case E_VS:
        case E_VZ:
        {
            break;
        }

        default:
        {
            L_THROW(0);

            break;
        }
    }
    #endif

    // Fix and check z_start and z_end

    switch ( b.get_v_type() )
    {
        case ZERO_VECTOR:
        case NORMAL_VECTOR:
        {
            break;
        }

        case SELECT_VECTOR:
        {
            THROW_ASSERT(b.get_sel_elm() >= (1+z_start));
            THROW_ASSERT(b.get_sel_elm() <= (e_size-z_end));

            z_start = b.get_sel_elm() - 1;
            z_end   = e_size - b.get_sel_elm();

            break;
        }

        case SCALAR_VECTOR:
        {
            THROW_ASSERT(z_start == 0);
            THROW_ASSERT(z_end == 0);

            break;
        }

        default:
        {
            L_THROW(0);

            break;
        }
    }

    // Do inversion, method dependant on form of y

    switch ( y.get_v_type() )
    {
        case ZERO_VECTOR:
        case SCALAR_VECTOR:
        case SELECT_VECTOR:
        {
            // In this case, we may only validly do the inversion if the
            // form of G is abstract

            switch ( get_e_type() )
            {
                case E_EMPTY:
                case E_ZERO:
                case E_DIAG:
                case E_LOWER:
                case E_UPPER:
                case E_SYMM:
                case E_ASYMM:
                {
                    L_THROW(0);

                    break;
                }

                case E_VS:
                {
                    y = b;

                    if ( b.get_v_type() != ZERO_VECTOR )
                    {
                        NZ_ASSERT(m_diag,zt);

                        y /= m_diag;
                    }

                    break;
                }

                case E_VZ:
                {
                    y = b;

                    if ( b.get_v_type() != ZERO_VECTOR )
                    {
                        throw INV_ERR;
                    }

                    break;
                }
            }

            break;
        }

        case NORMAL_VECTOR:
        {
            // Either the answer is zero, or we must do the full inversion

            if ( b.get_v_type() == ZERO_VECTOR )
            {
                y = b;
            }

            else
            {
                if ( e_size != 0 )
                {
                    // Allocate memory and setup problem variables

                    REPNEWB(bxx,L_DOUBLE,e_size);

                    for ( i = 1 ; i <= e_size ; i++ )
                    {
                        bxx[i-1] = b.get_offset_element(i-1);
                    }

                    // Do the inversion

                    L_DOUBLE *tempa;
                    L_DOUBLE tempb;

                    INV_GENERIC;

                    // Copy result and free memory.

                    for ( i = 1 ; i <= e_size ; i++ )
                    {
                        y[i-1] = bxx[i-1];
                    }

                    REPDELB(bxx);

                    bxx = NULL;
                }

                else
                {
                    y = b;
                }
            }

            break;
        }

        default:
        {
            L_THROW(0);

            break;
        }
    }

    FN_EXIT_POINT;
}

void fMATRIX::operator()(Matr_operation what, fMATRIX &y, const fMATRIX &b, double zt)
{
    FN_ENTRY_POINT

    // Undocumented actuality: this function also does rank-n updates

    if ( what != ARG_RANK_N )
    {
        long z_start = 0;
        long z_end = 0;
        long i,j,k,l;
        long e_size;
        fVECTOR *bxx;
        int inv_type;
        long b_size;
        long y_size;
        long y_type;

        THROW_ASSERT(!is_not_eff_square());
        THROW_ASSERT(b.is_not_eff_square() == y.is_not_eff_square());

        e_size = get_effective_height();

        if ( what == ARG_INVERT_LEFT )
        {
            inv_type = INV_IS_LEFT;

            b_size = b.get_effective_width();
            y_size = y.get_effective_width();

            j = b.get_effective_height();

            if ( b_size < 0 )
            {
                b_size = y_size;
                j = e_size;
            }

            THROW_ASSERT(b_size == y_size);

            if ( e_size >= 0 )
            {
                THROW_ASSERT(e_size == j);
                THROW_ASSERT(e_size == y.get_effective_height());
            }
        }

        else
        {
            THROW_ASSERT(what == ARG_INVERT_RIGHT);

            inv_type = INV_IS_RIGHT;

            b_size = b.get_effective_height();
            y_size = y.get_effective_height();

            j = b.get_effective_width();

            if ( b_size < 0 )
            {
                b_size = y_size;
                j = e_size;
            }

            THROW_ASSERT(b_size == y_size);

            if ( e_size >= 0 )
            {
                THROW_ASSERT(e_size == j);
                THROW_ASSERT(e_size == y.get_effective_width());
            }
        }

        // Do inversion, method dependant on form of y

        y_type = y.get_e_type();

        switch ( y_type )
        {
            case E_EMPTY:
            case E_ZERO:
            case E_DIAG:
            case E_LOWER:
            case E_UPPER:
            case E_SYMM:
            case E_ASYMM:
            {
                if ( e_size > 0 )
                {
                    REPNEWB(bxx,fVECTOR,e_size);

                    if ( inv_type == INV_IS_LEFT )
                    {
                        for ( i = 1 ; i <= e_size ; i++ )
                        {
                            (bxx[i-1]).pad_vector(b_size);

                            bxx[i-1] = b[i-1];
                        }
                    }

                    else
                    {
                        for ( i = 1 ; i <= e_size ; i++ )
                        {
                            (bxx[i-1]).pad_vector(b_size);

                            bxx[i-1] = b(i-1);
                        }
                    }

                    L_DOUBLE *tempa;
                    fVECTOR tempb('x',e_size);

                    INV_GENERIC;

                    if ( inv_type == INV_IS_LEFT )
                    {
                        for ( i = 1 ; i <= e_size ; i++ )
                        {
                            for ( j = 1 ; j <= y_size ; j++ )
                            {
                                y(i-1,j-1) = (bxx[i-1])[j-1];
                            }
                        }
                    }

                    else
                    {
                        for ( i = 1 ; i <= y_size ; i++ )
                        {
                            for ( j = 1 ; j <= e_size ; j++ )
                            {
                                y(i-1,j-1) = (bxx[j-1])[i-1];
                            }
                        }
                    }

                    REPDELB(bxx);

                    bxx = NULL;
                }

                else
                {
                    if ( e_size == 0 )
                    {
                        y = b;
                    }

                    else
                    {
                        if ( is_scalar() )
                        {
                            NZ_ASSERT(m_diag,zt);

                            y = b;
                            y /= m_diag;
                        }

                        if ( is_zero() )
                        {
                            throw INV_ERR;
                        }
                    }
                }

                break;
            }

            case E_VS:
            case E_VZ:
            {
                // In this case, we may only validly do the inversion if the
                // form of G is abstract

                switch ( get_e_type() )
                {
                    case E_EMPTY:
                    case E_ZERO:
                    case E_DIAG:
                    case E_LOWER:
                    case E_UPPER:
                    case E_SYMM:
                    case E_ASYMM:
                    {
                        L_THROW(0);

                        break;
                    }

                    case E_VS:
                    {
                        y = b;

                        if ( b.get_e_type() != E_VZ )
                        {
                            NZ_ASSERT(m_diag,zt);

                            y /= m_diag;
                        }

                        break;
                    }

                    case E_VZ:
                    {
                        y = b;

                        if ( b.get_e_type() != E_VZ )
                        {
                            throw INV_ERR;
                        }

                        break;
                    }
                }

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
        // and here is the rank-n bit

        // G := G + y.b.y'

        long i,j,k,n;
        long eff_size;
        L_DOUBLE ax;
        fVECTOR bx;

        // because this is an unofficial, make some assumptions, viz:
        //
        // sizes line up ok
        // no matrices are zero
        // both matrices are standard
        // Y is a symmetric matrix

        eff_size = get_effective_height();
        n = b.get_effective_height();

        if ( ( eff_size > 0 ) && ( n > 0 ) )
        {
            fix_symm();

            for ( i = 1 ; i <= n ; i++ )
            {
                bx = y(i-1);
                ax = b.get_offset_element(i-1,i-1);

                (*this)(ARG_RANK_ONE,bx,ax);

                if ( n > 1 )
                {
                    for ( j = 1 ; j <= n-1 ; j++ )
                    {
                        for ( k = j+1 ; k <= n ; k++ )
                        {
                            ax = 0.5;

                            bx  = y(j-1);
                            bx *= b.get_offset_element(j-1,k-1);
                            bx += y(k-1);

                            (*this)(ARG_RANK_ONE,bx,ax);
                        }
                    }

                    for ( j = 1 ; j <= n-1 ; j++ )
                    {
                        for ( k = j+1 ; k <= n ; k++ )
                        {
                            ax = -0.5;

                            bx  = y(j-1);
                            bx *= b.get_offset_element(j-1,k-1);
                            bx -= y(k-1);

                            (*this)(ARG_RANK_ONE,bx,ax);
                        }
                    }
                }
            }
        }
    }

    FN_EXIT_POINT;
}




//
// Matrix pivoting/resizing:
//

void fMATRIX::squareswap(long i, long j)
{
    FN_ENTRY_POINT

    l_xMatrix *Vi;
    l_xMatrix *Vj;
    l_xMatrix *Vm;

    long m;                                                             
    L_DOUBLE x;                                                        
                                                                        
    if ( i == j )                                                       
    {                                                                   
        FN_EXIT_POINT;                                                         
    }                                                                   
                                                                        
    THROW_ASSERT(i >= 1);                                                     
    THROW_ASSERT(i <= size);                                                  
    THROW_ASSERT(j >= 1);                                                     
    THROW_ASSERT(j <= size);                                                  
                                                                        
    if ( i > j )                                                        
    {                                                                   
        m = j;                                                          
        j = i;                                                          
        i = m;                                                          
    }                                                                   
                                                                        
    switch ( m_type )                                                   
    {                                                                   
        case ASYMMETRIC_MATRIX:                                         
        {                                                               
            // i-1   [ A  o  Q  t  R  ]    [ A  t  Q  o  R ] i-1
            // 1     [ b' f  p' s  u' ]    [ d' m  k' h  r'] 1
            // j-i-1 [ C  g  J  q  S  ] -> [ C  q  J  g  S ] j-i-1
            // 1     [ d' h  k' m  r' ]    [ b' s  p' f  u'] 1
            // N-j   [ E  i  L  n  P  ]    [ E  n  L  i  P ] N-j
                                                                        
            Vi = GET_XMATRIX_N(i);                                      
            Vj = GET_XMATRIX_N(j);                                      
                                                                        
            x        = Vi->diag;                                        
            Vi->diag = Vj->diag;                                        
            Vj->diag = x;                                               
                                                                        
            if ( j-i > 1 )                                              
            {                                                           
                Vm = Vi->next;                                          
                                                                        
                for ( m = i+1 ; m < j ; m++ )                           
                {                                                       
                    x              = (Vm->row)[i-1];                    
                    (Vm->row)[i-1] = (Vj->col)[m-1];                    
                    (Vj->col)[m-1] = x;                                 
                                                                        
                    x              = (Vm->col)[i-1];                    
                    (Vm->col)[i-1] = (Vj->row)[m-1];                    
                    (Vj->row)[m-1] = x;                                 
                                                                        
                    Vm = Vm->next;                                      
                }                                                       
            }                                                           
                                                                        
            if ( i > 1 )                                                
            {                                                           
                for ( m = 1 ; m < i ; m++ )                             
                {                                                       
                    x              = (Vi->row)[m-1];                    
                    (Vi->row)[m-1] = (Vj->row)[m-1];                    
                    (Vj->row)[m-1] = x;                                 
                                                                        
                    x              = (Vi->col)[m-1];                    
                    (Vi->col)[m-1] = (Vj->col)[m-1];                    
                    (Vj->col)[m-1] = x;                                 
                }                                                       
            }                                                           
                                                                        
            x              = (Vi->row)[i-1];                            
            (Vi->row)[i-1] = (Vj->row)[j-1];                            
            (Vj->row)[j-1] = x;                                         
                                                                        
            x              = (Vi->col)[i-1];                            
            (Vi->col)[i-1] = (Vj->col)[j-1];                            
            (Vj->col)[j-1] = x;                                         
                                                                        
            x              = (Vj->row)[i-1];                            
            (Vj->row)[i-1] = (Vj->col)[i-1];                            
            (Vj->col)[i-1] = x;                                         
                                                                        
            if ( j < size )                                             
            {                                                           
                Vm = Vj->next;                                          
                                                                        
                for ( m = (j+1) ; m <= size ; m++ )                     
                {                                                       
                    x              = (Vm->row)[i-1];                    
                    (Vm->row)[i-1] = (Vm->row)[j-1];                    
                    (Vm->row)[j-1] = x;                                 
                                                                        
                    x              = (Vm->col)[i-1];                    
                    (Vm->col)[i-1] = (Vm->col)[j-1];                    
                    (Vm->col)[j-1] = x;                                 
                                                                        
                    Vm = Vm->next;                                      
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case SYMMETRIC_MATRIX:                                          
        case LOWER_TRIANGULAR_MATRIX:                                   
        {                                                               
            // i-1   [ A  '* '* '* '* ]    [ A  '* '*  '* '* ] i-1
            // 1     [ b' f  '* '* '* ]    [ d' m  '*  '* '* ] 1
            // j-i-1 [ C  g  J  '* '* ] -> [ C  k* J   '* '* ] j-i-1
            // 1     [ d' h  k' m  '* ]    [ b' h* g'* f  '* ] 1
            // N-j   [ E  i  L  n  P  ]    [ E  n  L   i  P  ] N-j
                                                                        
            Vi = GET_XMATRIX_N(i);                                      
            Vj = GET_XMATRIX_N(j);                                      
                                                                        
            x        = Vi->diag;                                        
            Vi->diag = Vj->diag;                                        
            Vj->diag = x;                                               
                                                                        
            if ( i > 1 )                                                
            {                                                           
                for ( m = 1 ; m < i ; m++ )                             
                {                                                       
                    x              = (Vi->row)[m-1];                    
                    (Vi->row)[m-1] = (Vj->row)[m-1];                    
                    (Vj->row)[m-1] = x;                                 
                }                                                       
            }                                                           
                                                                        
            x              = (Vi->row)[i-1];                            
            (Vi->row)[i-1] = (Vj->row)[j-1];                            
            (Vj->row)[j-1] = x;                                         
                                                                        
            if ( j > (i+1) )                                            
            {                                                           
                Vm = Vi->next;                                          
                                                                        
                for ( m = (i+1) ; m < j ; m++ )                         
                {                                                       
                    x              = (Vj->row)[m-1];                    
                    (Vj->row)[m-1] = (Vm->row)[i-1];                    
                    (Vm->row)[i-1] = x;                                 
                                                                        
                    Vm = Vm->next;                                      
                }                                                       
            }                                                           
                                                                        
            if ( j < size )                                             
            {                                                           
                Vm = Vj->next;                                          
                                                                        
                for ( m = (j+1) ; m <= size ; m++ )                     
                {                                                       
                    x              = (Vm->row)[i-1];                    
                    (Vm->row)[i-1] = (Vm->row)[j-1];                    
                    (Vm->row)[j-1] = x;                                 
                                                                        
                    Vm = Vm->next;                                      
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case UPPER_TRIANGULAR_MATRIX:                                   
        {                                                               
            Vi = GET_XMATRIX_N(i);                                      
            Vj = GET_XMATRIX_N(j);                                      
                                                                        
            x        = Vi->diag;                                        
            Vi->diag = Vj->diag;                                        
            Vj->diag = x;                                               
                                                                        
            if ( i > 1 )                                                
            {                                                           
                for ( m = 1 ; m < i ; m++ )                             
                {                                                       
                    x              = (Vi->col)[m-1];                    
                    (Vi->col)[m-1] = (Vj->col)[m-1];                    
                    (Vj->col)[m-1] = x;                                 
                }                                                       
            }                                                           
                                                                        
            x              = (Vi->col)[i-1];                            
            (Vi->col)[i-1] = (Vj->col)[j-1];                            
            (Vj->col)[j-1] = x;                                         
                                                                        
            if ( j > (i+1) )                                            
            {                                                           
                Vm = Vi->next;                                          
                                                                        
                for ( m = (i+1) ; m < j ; m++ )                         
                {                                                       
                    x              = (Vj->col)[m-1];                    
                    (Vj->col)[m-1] = (Vm->col)[i-1];                    
                    (Vm->col)[i-1] = x;                                 
                                                                        
                    Vm = Vm->next;                                      
                }                                                       
            }                                                           
                                                                        
            if ( j < size )                                             
            {                                                           
                Vm = Vj->next;                                          
                                                                        
                for ( m = (j+1) ; m <= size ; m++ )                     
                {                                                       
                    x              = (Vm->col)[i-1];                    
                    (Vm->col)[i-1] = (Vm->col)[j-1];                    
                    (Vm->col)[j-1] = x;                                 
                                                                        
                    Vm = Vm->next;                                      
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case DIAGONAL_MATRIX:                                           
        {                                                               
            // i-1   [ A             ]    [ A             ] i-1
            // 1     [    f          ]    [    m          ] 1
            // j-i-1 [       J       ] -> [       J       ] j-i-1
            // 1     [          m    ]    [          f    ] 1
            // N-j   [             P ]    [             P ] N-j
                                                                        
            Vi = GET_XMATRIX_N(i);                                      
            Vj = GET_XMATRIX_N(j);                                      
                                                                        
            x        = Vi->diag;                                        
            Vi->diag = Vj->diag;                                        
            Vj->diag = x;                                               
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case SCALAR_MATRIX:                                             
        case ZERO_MATRIX:                                               
        {                                                               
            break;                                                      
        }                                                               
                                                                        
        default:                                                        
        {                                                               
            L_THROW(0);                                                 
                                                                        
            break;                                                      
        }                                                               
    }                                                                   

    FN_EXIT_POINT;
}

void fMATRIX::bswap(long i, long j)
{
    FN_ENTRY_POINT

    l_xMatrix *Vi;
    l_xMatrix *Vj;
    l_xMatrix *Vm;
    l_xMatrix *Vn;

    long m,n;                                                           
    L_DOUBLE x;                                                        
    L_DOUBLE y;                                                        
                                                                        
    if ( i == j )                                                       
    {                                                                   
        FN_EXIT_POINT;                                                         
    }                                                                   
                                                                        
    THROW_ASSERT(i >= 1);                                                     
    THROW_ASSERT(i <= size);                                                  
    THROW_ASSERT(j >= 1);                                                     
    THROW_ASSERT(j <= size);                                                  
                                                                        
    if ( j > i )                                                        
    {                                                                   
        m = j;                                                          
        j = i;                                                          
        i = m;                                                          
    }                                                                   
                                                                        
    switch ( m_type )                                                   
    {                                                                   
        case ASYMMETRIC_MATRIX:                                         
        {                                                               
            // j-1 [ A  K  l  M ]    [ A  l  K  M ] j-1
            // i-j [ B  E  n  O ] -> [ c' h  f' p'] 1
            // 1   [ c' f' h  p']    [ B  n  E  O ] i-j
            // N-i [ D  G  i  J ]    [ D  i  G  J ] N-i
                                                                        
            Vj = GET_XMATRIX_N(j);                                      
            Vi = GET_XMATRIX_N(i);                                      
                                                                        
            Vn = Vi;                                                    
                                                                        
            x = Vn->diag;                                               
                                                                        
            for ( m = i ; m > j ; m-- )                                 
            {                                                           
                Vn->diag = (Vn->prev)->diag;                            
                                                                        
                Vn = Vn->prev;                                          
            }                                                           
                                                                        
            Vn->diag = x;                                               
                                                                        
            for ( m = i-1 ; m >= j ; m-- )                              
            {                                                           
                x = (Vi->col)[m-1];                                     
                                                                        
                Vn = Vi;                                                
                                                                        
                if ( m > j )                                            
                {                                                       
                    for ( n = m ; n > j ; n-- )                         
                    {                                                   
                        (Vn->col)[n-1] = ((Vn->prev)->col)[n-2];        
                                                                        
                        Vn = Vn->prev;                                  
                    }                                                   
                }                                                       
                                                                        
                (Vn->col)[j-1] = (Vi->row)[j+((i-1)-m)-1];              
                                                                        
                Vn = Vi;                                                
                                                                        
                if ( m < i-1 )                                          
                {                                                       
                    for ( n = j+((i-1)-m) ; n > j ; n-- )               
                    {                                                   
                        (Vn->row)[n-1] = ((Vn->prev)->row)[n-2];        
                                                                        
                        Vn = Vn->prev;                                  
                    }                                                   
                }                                                       
                                                                        
                (Vn->row)[j-1] = x;                                     
            }                                                           
                                                                        
            if ( j > 1 )                                                
            {                                                           
                for ( m = 1 ; m < j ; m++ )                             
                {                                                       
                    Vn = Vi;                                            
                                                                        
                    x = (Vn->row)[m-1];                                 
                    y = (Vn->col)[m-1];                                 
                                                                        
                    for ( n = i ; n > j ; n-- )                         
                    {                                                   
                        (Vn->row)[m-1] = ((Vn->prev)->row)[m-1];        
                        (Vn->col)[m-1] = ((Vn->prev)->col)[m-1];        
                                                                        
                        Vn = Vn->prev;                                  
                    }                                                   
                                                                        
                    (Vn->row)[m-1] = x;                                 
                    (Vn->col)[m-1] = y;                                 
                }                                                       
            }                                                           
                                                                        
            if ( i < size )                                             
            {                                                           
                Vm = Vi->next;                                          
                                                                        
                for ( m = (i+1) ; m <= size ; m++ )                     
                {                                                       
                    x = (Vm->row)[i-1];                                 
                    y = (Vm->col)[i-1];                                 
                                                                        
                    for ( n = i ; n > j ; n-- )                         
                    {                                                   
                        (Vm->row)[n-1] = (Vm->row)[n-2];                
                        (Vm->col)[n-1] = (Vm->col)[n-2];                
                    }                                                   
                                                                        
                    (Vm->row)[j-1] = x;                                 
                    (Vm->col)[j-1] = y;                                 
                                                                        
                    Vm = Vm->next;                                      
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case SYMMETRIC_MATRIX:                                          
        case LOWER_TRIANGULAR_MATRIX:                                   
        {                                                               
            // j-1 [ A  '* '* '* ]    [ A  '* '* '* ] j-1
            // i-j [ B  E  '* '* ] -> [ c' h  '* '* ] 1
            // 1   [ c' f' h  '* ]    [ B  f* E  '* ] i-j
            // N-i [ D  G  i  J  ]    [ D  i  G  J  ] N-i
                                                                        
            Vj = GET_XMATRIX_N(j);                                      
            Vi = GET_XMATRIX_N(i);                                      
                                                                        
            Vn = Vi;                                                    
                                                                        
            x = Vn->diag;                                               
                                                                        
            for ( m = i ; m > j ; m-- )                                 
            {                                                           
                Vn->diag = (Vn->prev)->diag;                            
                                                                        
                Vn = Vn->prev;                                          
            }                                                           
                                                                        
            Vn->diag = x;                                               
                                                                        
            if ( j > 1 )                                                
            {                                                           
                for ( m = 1 ; m < j ; m++ )                             
                {                                                       
                    Vn = Vi;                                            
                                                                        
                    x = (Vn->row)[m-1];                                 
                                                                        
                    for ( n = i ; n > j ; n-- )                         
                    {                                                   
                        (Vn->row)[m-1] = ((Vn->prev)->row)[m-1];        
                                                                        
                        Vn = Vn->prev;                                  
                    }                                                   
                                                                        
                    (Vn->row)[m-1] = x;                                 
                }                                                       
            }                                                           
                                                                        
            for ( m = i ; m >= j ; m-- )                                
            {                                                           
                Vn = Vi;                                                
                                                                        
                if ( m > j )                                            
                {                                                       
                    x = (Vn->row)[m-1];                                 
                                                                        
                    for ( n = m ; n > j ; n-- )                         
                    {                                                   
                        (Vn->row)[n-1] = ((Vn->prev)->row)[n-2];        
                                                                        
                        Vn = Vn->prev;                                  
                    }                                                   
                                                                        
                    (Vn->row)[j-1] = x;                                 
                }                                                       
            }                                                           
                                                                        
            Vm = Vj->next;  m = j+1;                                    
            Vn = Vi;        n = i;                                      
                                                                        
            if ( n > m )                                                
            {                                                           
                repeat_a:                                               
                                                                        
                x              = (Vm->row)[j-1];                        
                (Vm->row)[j-1] = (Vn->row)[j-1];                        
                (Vn->row)[j-1] = x;                                     
                                                                        
                Vm = Vm->next;   m++;                                   
                Vn = Vn->prev;   n--;                                   
                                                                        
                if ( n > m ) { goto repeat_a; }                         
            }                                                           
                                                                        
            if ( i < size )                                             
            {                                                           
                Vm = Vi->next;                                          
                                                                        
                for ( m = (i+1) ; m <= size ; m++ )                     
                {                                                       
                    x = (Vm->row)[i-1];                                 
                                                                        
                    for ( n = i ; n > j ; n-- )                         
                    {                                                   
                        (Vm->row)[n-1] = (Vm->row)[n-2];                
                    }                                                   
                                                                        
                    (Vm->row)[j-1] = x;                                 
                                                                        
                    Vm = Vm->next;                                      
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case UPPER_TRIANGULAR_MATRIX:                                   
        {                                                               
            Vj = GET_XMATRIX_N(j);                                      
            Vi = GET_XMATRIX_N(i);                                      
                                                                        
            Vn = Vi;                                                    
                                                                        
            x = Vn->diag;                                               
                                                                        
            for ( m = i ; m > j ; m-- )                                 
            {                                                           
                Vn->diag = (Vn->prev)->diag;                            
                                                                        
                Vn = Vn->prev;                                          
            }                                                           
                                                                        
            Vn->diag = x;                                               
                                                                        
            if ( j > 1 )                                                
            {                                                           
                for ( m = 1 ; m < j ; m++ )                             
                {                                                       
                    Vn = Vi;                                            
                                                                        
                    x = (Vn->col)[m-1];                                 
                                                                        
                    for ( n = i ; n > j ; n-- )                         
                    {                                                   
                        (Vn->col)[m-1] = ((Vn->prev)->col)[m-1];        
                                                                        
                        Vn = Vn->prev;                                  
                    }                                                   
                                                                        
                    (Vn->col)[m-1] = x;                                 
                }                                                       
            }                                                           
                                                                        
            for ( m = i ; m >= j ; m-- )                                
            {                                                           
                Vn = Vi;                                                
                                                                        
                if ( m > j )                                            
                {                                                       
                    x = (Vn->col)[m-1];                                 
                                                                        
                    for ( n = m ; n > j ; n-- )                         
                    {                                                   
                        (Vn->col)[n-1] = ((Vn->prev)->col)[n-2];        
                                                                        
                        Vn = Vn->prev;                                  
                    }                                                   
                                                                        
                    (Vn->col)[j-1] = x;                                 
                }                                                       
            }                                                           
                                                                        
            Vm = Vj->next;  m = j+1;                                    
            Vn = Vi;        n = i;                                      
                                                                        
            if ( n > m )                                                
            {                                                           
                repeat_b:                                               
                                                                        
                x              = (Vm->col)[j-1];                        
                (Vm->col)[j-1] = (Vn->col)[j-1];                        
                (Vn->col)[j-1] = x;                                     
                                                                        
                Vm = Vm->next;   m++;                                   
                Vn = Vn->prev;   n--;                                   
                                                                        
                if ( n > m ) { goto repeat_b; }                         
            }                                                           
                                                                        
            if ( i < size )                                             
            {                                                           
                Vm = Vi->next;                                          
                                                                        
                for ( m = (i+1) ; m <= size ; m++ )                     
                {                                                       
                    x = (Vm->col)[i-1];                                 
                                                                        
                    for ( n = i ; n > j ; n-- )                         
                    {                                                   
                        (Vm->col)[n-1] = (Vm->col)[n-2];                
                    }                                                   
                                                                        
                    (Vm->col)[j-1] = x;                                 
                                                                        
                    Vm = Vm->next;                                      
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case DIAGONAL_MATRIX:                                           
        {                                                               
            // j-1 [ A          ]    [ A          ] j-1
            // i-j [    E       ] -> [    h       ] 1  
            // 1   [       h    ]    [       E    ] i-j
            // N-i [          J ]    [          J ] N-i
                                                                        
            Vj = GET_XMATRIX_N(j);                                      
            Vi = GET_XMATRIX_N(i);                                      
                                                                        
            Vn = Vi;                                                    
                                                                        
            x = Vn->diag;                                               
                                                                        
            for ( m = i ; m > j ; m-- )                                 
            {                                                           
                Vn->diag = (Vn->prev)->diag;                            
                                                                        
                Vn = Vn->prev;                                          
            }                                                           
                                                                        
            Vn->diag = x;                                               
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case SCALAR_MATRIX:                                             
        case ZERO_MATRIX:                                               
        {                                                               
            break;                                                      
        }                                                               
                                                                        
        default:                                                        
        {                                                               
            L_THROW(0);                                                 
                                                                        
            break;                                                      
        }                                                               
    }                                                                   

    FN_EXIT_POINT;
}

void fMATRIX::fswap(long i, long j)
{
    FN_ENTRY_POINT

    l_xMatrix *Vi;
    l_xMatrix *Vj;
    l_xMatrix *Vm;
    l_xMatrix *Vn;
    l_xMatrix *Vk;

    long m,n;                                                           
    L_DOUBLE x;                                                        
    L_DOUBLE y;                                                        
                                                                        
    if ( i == j )                                                       
    {                                                                   
        FN_EXIT_POINT;                                                         
    }                                                                   
                                                                        
    THROW_ASSERT(i >= 1);                                                     
    THROW_ASSERT(i <= size);                                                  
    THROW_ASSERT(j >= 1);                                                     
    THROW_ASSERT(j <= size);                                                  
                                                                        
    if ( i > j )                                                        
    {                                                                   
        m = j;                                                          
        j = i;                                                          
        i = m;                                                          
    }                                                                   
                                                                        
    switch ( m_type )                                                   
    {                                                                   
        case ASYMMETRIC_MATRIX:                                         
        {                                                               
            // i-1 [ A  l  M  N ]    [ A  M  l  N ] i-1
            // 1   [ b' e  o' p'] -> [ C  H  f  Q ] j-i
            // j-i [ C  f  H  Q ]    [ b' o' e  p'] 1  
            // N-j [ D  g  J  K ]    [ D  J  g  K ] N-j
                                                                        
            Vi = GET_XMATRIX_N(i);                                      
            Vj = GET_XMATRIX_N(j);                                      
                                                                        
            Vn = Vi;                                                    
                                                                        
            x = Vn->diag;                                               
                                                                        
            for ( m = i ; m < j ; m++ )                                 
            {                                                           
                Vn->diag = (Vn->next)->diag;                            
                                                                        
                Vn = Vn->next;                                          
            }                                                           
                                                                        
            Vn->diag = x;                                               
                                                                        
            Vm = Vi->next;                                              
            Vk = Vj;                                                    
                                                                        
            for ( m = i+1 ; m <= j ; m++ )                              
            {                                                           
                x = (Vm->row)[i-1];                                     
                                                                        
                Vn = Vm;                                                
                                                                        
                if ( m < j )                                            
                {                                                       
                    for ( n = m ; n < j ; n++ )                         
                    {                                                   
                        (Vn->row)[(n-m+i)-1] = ((Vn->next)->row)[(n-m+i)]; 
                                                                        
                        Vn = Vn->next;                                  
                    }                                                   
                }                                                       
                                                                        
                (Vn->row)[(j-m+i)-1] = (Vk->col)[i-1];                  
                                                                        
                Vn = Vk;                                                
                                                                        
                if ( m > i+1 )                                          
                {                                                       
                    for ( n = i ; n < m-1 ; n++ )                       
                    {                                                   
                        (Vn->col)[n-1] = ((Vn->next)->col)[n];          
                                                                        
                        Vn = Vn->next;                                  
                    }                                                   
                }                                                       
                                                                        
                (Vn->col)[(m-1)-1] = x;                                 
                                                                        
                Vm = Vm->next;                                          
                Vk = Vk->prev;                                          
            }                                                           
                                                                        
            if ( i > 1 )                                                
            {                                                           
                for ( m = 1 ; m < i ; m++ )                             
                {                                                       
                    Vn = Vi;                                            
                                                                        
                    x = (Vn->row)[m-1];                                 
                    y = (Vn->col)[m-1];                                 
                                                                        
                    for ( n = i ; n < j ; n++ )                         
                    {                                                   
                        (Vn->row)[m-1] = ((Vn->next)->row)[m-1];        
                        (Vn->col)[m-1] = ((Vn->next)->col)[m-1];        
                                                                        
                        Vn = Vn->next;                                  
                    }                                                   
                                                                        
                    (Vn->row)[m-1] = x;                                 
                    (Vn->col)[m-1] = y;                                 
                }                                                       
            }                                                           
                                                                        
            if ( j < size )                                             
            {                                                           
                Vm = Vj->next;                                          
                                                                        
                for ( m = (j+1) ; m <= size ; m++ )                     
                {                                                       
                    x = (Vm->row)[i-1];                                 
                    y = (Vm->col)[i-1];                                 
                                                                        
                    for ( n = i ; n < j ; n++ )                         
                    {                                                   
                        (Vm->row)[n-1] = (Vm->row)[n];                  
                        (Vm->col)[n-1] = (Vm->col)[n];                  
                    }                                                   
                                                                        
                    (Vm->row)[j-1] = x;                                 
                    (Vm->col)[j-1] = y;                                 
                                                                        
                    Vm = Vm->next;                                      
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case SYMMETRIC_MATRIX:                                          
        case LOWER_TRIANGULAR_MATRIX:                                   
        {                                                               
            // i-1 [ A  '* '* '* ]    [ A  '*  '* '* ] i-1
            // 1   [ b' e  '* '* ] -> [ C  H   '* '* ] j-i
            // j-i [ C  f  H  '* ]    [ b' f'* e  '* ] 1  
            // N-j [ D  g  J  K  ]    [ D  J   g  K  ] N-j
                                                                        
            Vi = GET_XMATRIX_N(i);                                      
            Vj = GET_XMATRIX_N(j);                                      
                                                                        
            Vn = Vi;                                                    
                                                                        
            x = Vn->diag;                                               
                                                                        
            for ( m = i ; m < j ; m++ )                                 
            {                                                           
                Vn->diag = (Vn->next)->diag;                            
                                                                        
                Vn = Vn->next;                                          
            }                                                           
                                                                        
            Vn->diag = x;                                               
                                                                        
            if ( i > 1 )                                                
            {                                                           
                for ( m = 1 ; m < i ; m++ )                             
                {                                                       
                    Vn = Vi;                                            
                                                                        
                    x = (Vn->row)[m-1];                                 
                                                                        
                    for ( n = i ; n < j ; n++ )                         
                    {                                                   
                        (Vn->row)[m-1] = ((Vn->next)->row)[m-1];        
                                                                        
                        Vn = Vn->next;                                  
                    }                                                   
                                                                        
                    (Vn->row)[m-1] = x;                                 
                }                                                       
            }                                                           
                                                                        
            Vm = Vi;                                                    
                                                                        
            for ( m = i ; m <= j ; m++ )                                
            {                                                           
                Vn = Vm;                                                
                                                                        
                if ( m < j )                                            
                {                                                       
                    x = (Vn->row)[i-1];                                 
                                                                        
                    for ( n = m ; n < j ; n++ )                         
                    {                                                   
                        (Vn->row)[n-(m-i)-1] = ((Vn->next)->row)[n-(m-i)]; 
                                                                        
                        Vn = Vn->next;                                  
                    }                                                   
                                                                        
                    (Vn->row)[j-(m-i)-1] = x;                           
                }                                                       
                                                                        
                if ( m < j )                                            
                {                                                       
                    Vm = Vm->next;                                      
                }                                                       
            }                                                           
                                                                        
            m = i;                                                      
            n = j-1;                                                    
                                                                        
            if ( n > m )                                                
            {                                                           
                repeat_a:                                               
                                                                        
                x              = (Vj->row)[n-1];                        
                (Vj->row)[n-1] = (Vj->row)[m-1];                        
                (Vj->row)[m-1] = x;                                     
                                                                        
                m++;                                                    
                n--;                                                    
                                                                        
                if ( n > m ) { goto repeat_a; }                         
            }                                                           
                                                                        
            if ( j < size )                                             
            {                                                           
                Vm = Vj->next;                                          
                                                                        
                for ( m = (j+1) ; m <= size ; m++ )                     
                {                                                       
                    x = (Vm->row)[i-1];                                 
                                                                        
                    for ( n = i ; n < j ; n++ )                         
                    {                                                   
                        (Vm->row)[n-1] = (Vm->row)[n];                  
                    }                                                   
                                                                        
                    (Vm->row)[j-1] = x;                                 
                                                                        
                    Vm = Vm->next;                                      
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case UPPER_TRIANGULAR_MATRIX:                                   
        {                                                               
            Vi = GET_XMATRIX_N(i);                                      
            Vj = GET_XMATRIX_N(j);                                      
                                                                        
            Vn = Vi;                                                    
                                                                        
            x = Vn->diag;                                               
                                                                        
            for ( m = i ; m < j ; m++ )                                 
            {                                                           
                Vn->diag = (Vn->next)->diag;                            
                                                                        
                Vn = Vn->next;                                          
            }                                                           
                                                                        
            Vn->diag = x;                                               
                                                                        
            if ( i > 1 )                                                
            {                                                           
                for ( m = 1 ; m < i ; m++ )                             
                {                                                       
                    Vn = Vi;                                            
                                                                        
                    x = (Vn->col)[m-1];                                 
                                                                        
                    for ( n = i ; n < j ; n++ )                         
                    {                                                   
                        (Vn->col)[m-1] = ((Vn->next)->col)[m-1];        
                                                                        
                        Vn = Vn->next;                                  
                    }                                                   
                                                                        
                    (Vn->col)[m-1] = x;                                 
                }                                                       
            }                                                           
                                                                        
            Vm = Vi;                                                    
                                                                        
            for ( m = i ; m <= j ; m++ )                                
            {                                                           
                Vn = Vm;                                                
                                                                        
                if ( m < j )                                            
                {                                                       
                    x = (Vn->col)[i-1];                                 
                                                                        
                    for ( n = m ; n < j ; n++ )                         
                    {                                                   
                        (Vn->col)[n-(m-i)-1] = ((Vn->next)->col)[n-(m-i)]; 
                                                                        
                        Vn = Vn->next;                                  
                    }                                                   
                                                                        
                    (Vn->col)[j-(m-i)-1] = x;                           
                }                                                       
                                                                        
                if ( m < j )                                            
                {                                                       
                    Vm = Vm->next;                                      
                }                                                       
            }                                                           
                                                                        
            m = i;                                                      
            n = j-1;                                                    
                                                                        
            if ( n > m )                                                
            {                                                           
                repeat_b:                                               
                                                                        
                x              = (Vj->col)[n-1];                        
                (Vj->col)[n-1] = (Vj->col)[m-1];                        
                (Vj->col)[m-1] = x;                                     
                                                                        
                m++;                                                    
                n--;                                                    
                                                                        
                if ( n > m ) { goto repeat_b; }                         
            }                                                           
                                                                        
            if ( j < size )                                             
            {                                                           
                Vm = Vj->next;                                          
                                                                        
                for ( m = (j+1) ; m <= size ; m++ )                     
                {                                                       
                    x = (Vm->col)[i-1];                                 
                                                                        
                    for ( n = i ; n < j ; n++ )                         
                    {                                                   
                        (Vm->col)[n-1] = (Vm->col)[n];                  
                    }                                                   
                                                                        
                    (Vm->col)[j-1] = x;                                 
                                                                        
                    Vm = Vm->next;                                      
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case DIAGONAL_MATRIX:                                           
        {                                                               
            // i-1 [ A          ]    [ A          ] i-1
            // 1   [    e       ] -> [    H       ] j-i
            // j-i [       H    ]    [       e    ] 1  
            // N-j [          K ]    [          K ] N-j
                                                                        
            Vi = GET_XMATRIX_N(i);                                      
            Vj = GET_XMATRIX_N(j);                                      
                                                                        
            Vn = Vi;                                                    
                                                                        
            x = Vn->diag;                                               
                                                                        
            for ( m = i ; m < j ; m++ )                                 
            {                                                           
                Vn->diag = (Vn->next)->diag;                            
                                                                        
                Vn = Vn->next;                                          
            }                                                           
                                                                        
            Vn->diag = x;                                               
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case SCALAR_MATRIX:                                             
        case ZERO_MATRIX:                                               
        {                                                               
            break;                                                      
        }                                                               
                                                                        
        default:                                                        
        {                                                               
            L_THROW(0);                                                 
                                                                        
            break;                                                      
        }                                                               
    }                                                                   

    FN_EXIT_POINT;
}

void fMATRIX::addend(void)
{
    FN_ENTRY_POINT

    switch ( m_type )                                                   
    {                                                                   
        case ASYMMETRIC_MATRIX:                                         
        {                                                               
            fVECTOR zero_vec_ere;                                   
                                                                        
            addend(0.0,zero_vec_ere,zero_vec_ere);                 
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case SYMMETRIC_MATRIX:                                          
        case LOWER_TRIANGULAR_MATRIX:                                   
        case UPPER_TRIANGULAR_MATRIX:                                   
        {                                                               
            fVECTOR zero_vec_ere;                                   
                                                                        
            addend(0.0,zero_vec_ere);                              
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case DIAGONAL_MATRIX:                                           
        {                                                               
            addend(0.0);                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case SCALAR_MATRIX:                                             
        {                                                               
            THROW_ASSERT(m_diag == 0.0);
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case ZERO_MATRIX:                                               
        {                                                               
            break;                                                      
        }                                                               
                                                                        
        default:                                                        
        {                                                               
            L_THROW(0);                                                 
                                                                        
            break;                                                      
        }                                                               
    }                                                                   


    FN_EXIT_POINT;
}

void fMATRIX::addend(const   double &b)
{
    FN_ENTRY_POINT

    l_xMatrix *here;

    switch ( m_type )                                                   
    {                                                                   
        case ASYMMETRIC_MATRIX:                                         
        {                                                               
            fVECTOR zero_vec_ere;                                   
                                                                        
            addend(b,zero_vec_ere,zero_vec_ere);                        
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case SYMMETRIC_MATRIX:                                          
        case LOWER_TRIANGULAR_MATRIX:                                   
        case UPPER_TRIANGULAR_MATRIX:                                   
        {                                                               
            fVECTOR zero_vec_ere;                                   
                                                                        
            addend(b,zero_vec_ere);                                     
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case DIAGONAL_MATRIX:                                           
        {                                                               
            size++;                                                     
                                                                        
            if ( size == 1 )                                            
            {                                                           
                REPNEW(all_it,l_xMatrix);                                  
                                                                        
                all_it->prev = NULL;                                    
                all_it->next = NULL;                                    
                                                                        
                reset_offsets();                                        
                                                                        
                here = all_it;                                          
            }                                                           
                                                                        
            else                                                        
            {                                                           
                here = GET_XMATRIX_N((size-1));                           
                                                                        
                REPNEW(here->next,l_xMatrix);                              
                                                                        
                (here->next)->prev = here;                              
                (here->next)->next = NULL;                              
                                                                        
                here = here->next;                                      
                                                                        
                if ( offset_end_effective_row == size-1 ) { offset_end_effective_row++; } 
                if ( offset_end_effective_col == size-1 ) { offset_end_effective_col++; } 
            }

            matrix_struct_lookup.addend();
			GET_XMATRIX_NS(size) = (void *)here;
                                                                        
            here->diag = b;                                             
                                                                        
            here->row = NULL;                                           
            here->col = NULL;                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case SCALAR_MATRIX:                                             
        {                                                               
            THROW_ASSERT(b == m_diag);
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case ZERO_MATRIX:                                               
        {                                                               
            THROW_ASSERT(b == 0.0);
                                                                        
            break;                                                      
        }                                                               
                                                                        
        default:                                                        
        {                                                               
            L_THROW(0);                                                 
                                                                        
            break;                                                      
        }                                                               
    }                                                                   

    FN_EXIT_POINT;
}

void fMATRIX::addend(const   double &b, const fVECTOR &x)
{
    FN_ENTRY_POINT

    l_xMatrix *here;

    long i;                                                             
                                                                        
    switch ( m_type )                                                   
    {                                                                   
        case ASYMMETRIC_MATRIX:                                         
        {                                                               
            fVECTOR y(x);                                           
                                                                        
            addend(b,x,y);                                              
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case DIAGONAL_MATRIX:                                           
        {                                                               
            make_asymmetric();                                          
                                                                        
            addend(b,x,x);                                              
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case SYMMETRIC_MATRIX:                                          
        {
            THROW_ASSERT((x.get_effective_size() < 0) || (x.get_effective_size() == size+1));
                                                                        
            size++;                                                     
                                                                        
            if ( size == 1 )                                            
            {                                                           
                REPNEW(all_it,l_xMatrix);                                  
                                                                        
                all_it->prev = NULL;                                    
                all_it->next = NULL;                                    
                                                                        
                reset_offsets();                                        
                                                                        
                here = all_it;                                          
            }                                                           
                                                                        
            else                                                        
            {                                                           
                here = GET_XMATRIX_N((size-1));                           
                                                                        
                REPNEW(here->next,l_xMatrix);                              
                                                                        
                (here->next)->prev = here;                              
                (here->next)->next = NULL;                              
                                                                        
                here = here->next;                                      
                                                                        
                if ( offset_end_effective_row == size-1 )               
                {                                                       
                    offset_end_effective_row++;                         
                }                                                       
                                                                        
                if ( offset_end_effective_col == size-1 )               
                {                                                       
                    offset_end_effective_col++;                         
                }                                                       
            }                                                           
                                                                        
            matrix_struct_lookup.addend();
			GET_XMATRIX_NS(size) = (void *)here;
                                                                        
            here->diag = b;                                             
                                                                        
            REPNEWB(here->row,L_DOUBLE,size);                          
                                                                        
            for ( i = 1 ; i <= size ; i++ )                             
            {                                                           
                (here->row)[i-1] = x.get_offset_element(i-1);           
            }                                                           
                                                                        
            here->col = here->row;                                      
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case LOWER_TRIANGULAR_MATRIX:                                   
        {                                                               
            THROW_ASSERT((x.get_effective_size() < 0) || (x.get_effective_size() == size+1));
                                                                        
            size++;                                                     
                                                                        
            if ( size == 1 )                                            
            {                                                           
                REPNEW(all_it,l_xMatrix);                                  
                                                                        
                all_it->prev = NULL;                                    
                all_it->next = NULL;                                    
                                                                        
                reset_offsets();                                        
                                                                        
                here = all_it;                                          
            }                                                           
                                                                        
            else                                                        
            {                                                           
                here = GET_XMATRIX_N((size-1));                           
                                                                        
                REPNEW(here->next,l_xMatrix);                              
                                                                        
                (here->next)->prev = here;                              
                (here->next)->next = NULL;                              
                                                                        
                here = here->next;                                      
                                                                        
                if ( offset_end_effective_row == size-1 ) { offset_end_effective_row++; } 
                if ( offset_end_effective_col == size-1 ) { offset_end_effective_col++; } 
            }                                                           
                                                                        
            matrix_struct_lookup.addend();
			GET_XMATRIX_NS(size) = (void *)here;
                                                                        
            here->diag = b;                                             
                                                                        
            REPNEWB(here->row,L_DOUBLE,size);                          
                                                                        
            for ( i = 1 ; i <= size ; i++ )                             
            {                                                           
                (here->row)[i-1] = x.get_offset_element(i-1);           
            }                                                           
                                                                        
            here->col = NULL;                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case UPPER_TRIANGULAR_MATRIX:                                   
        {                                                               
            THROW_ASSERT((x.get_effective_size() < 0) || (x.get_effective_size() == size+1));
                                                                        
            size++;                                                     
                                                                        
            if ( size == 1 )                                            
            {                                                           
                REPNEW(all_it,l_xMatrix);                                  
                                                                        
                all_it->prev = NULL;                                    
                all_it->next = NULL;                                    
                                                                        
                reset_offsets();                                        
                                                                        
                here = all_it;                                          
            }                                                           
                                                                        
            else                                                        
            {                                                           
                here = GET_XMATRIX_N((size-1));                           
                                                                        
                REPNEW(here->next,l_xMatrix);                              
                                                                        
                (here->next)->prev = here;                              
                (here->next)->next = NULL;                              
                                                                        
                here = here->next;                                      
                                                                        
                if ( offset_end_effective_row == size-1 ) { offset_end_effective_row++; } 
                if ( offset_end_effective_col == size-1 ) { offset_end_effective_col++; } 
            }                                                           
                                                                        
            matrix_struct_lookup.addend();
			GET_XMATRIX_NS(size) = (void *)here;
                                                                        
            here->diag = b;                                             
                                                                        
            REPNEWB(here->col,L_DOUBLE,size);                          
                                                                        
            for ( i = 1 ; i <= size ; i++ )                             
            {                                                           
                (here->col)[i-1] = x.get_offset_element(i-1);           
            }                                                           
                                                                        
            here->row = NULL;                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case SCALAR_MATRIX:                                             
        case ZERO_MATRIX:                                               
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

    FN_EXIT_POINT;
}

void fMATRIX::addend(const   double &b, const fVECTOR &x, const fVECTOR &y)
{
    FN_ENTRY_POINT

    l_xMatrix *here;

    long i;                                                             
                                                                        
    switch ( m_type )                                                   
    {                                                                   
        case SYMMETRIC_MATRIX:                                          
        case LOWER_TRIANGULAR_MATRIX:                                   
        case UPPER_TRIANGULAR_MATRIX:                                   
        case DIAGONAL_MATRIX:                                           
        {                                                               
            make_asymmetric();                                          
        }                                                               
                                                                        
        case ASYMMETRIC_MATRIX:                                         
        {                                                               
            THROW_ASSERT((x.get_effective_size() < 0) || (x.get_effective_size() == size+1));
            THROW_ASSERT((y.get_effective_size() < 0) || (y.get_effective_size() == size+1));
                                                                        
            size++;                                                     
                                                                        
            if ( size == 1 )                                            
            {                                                           
                REPNEW(all_it,l_xMatrix);                                  
                                                                        
                all_it->prev = NULL;                                    
                all_it->next = NULL;                                    
                                                                        
                reset_offsets();                                        
                                                                        
                here = all_it;                                          
            }                                                           
                                                                        
            else                                                        
            {                                                           
                here = GET_XMATRIX_N((size-1));                           
                                                                        
                REPNEW(here->next,l_xMatrix);                              
                                                                        
                (here->next)->prev = here;                              
                (here->next)->next = NULL;                              
                                                                        
                if ( offset_end_effective_row == size-1 )               
                {                                                       
                    offset_end_effective_row++;                         
                }                                                       
                                                                        
                if ( offset_end_effective_col == size-1 )               
                {                                                       
                    offset_end_effective_col++;                         
                }                                                       
                                                                        
                here = here->next;                                      
            }                                                           
                                                                        
            matrix_struct_lookup.addend();
			GET_XMATRIX_NS(size) = (void *)here;
                                                                        
            here->diag = b;                                             
                                                                        
            REPNEWB(here->row,L_DOUBLE,size);                          
            REPNEWB(here->col,L_DOUBLE,size);
                                                                        
            for ( i = 1 ; i <= size ; i++ )                             
            {                                                           
                (here->row)[i-1] = x.get_offset_element(i-1);           
                (here->col)[i-1] = y.get_offset_element(i-1);           
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case SCALAR_MATRIX:                                             
        case ZERO_MATRIX:                                               
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

    FN_EXIT_POINT;
}

void fMATRIX::addend(const c_double &b)
{
    FN_ENTRY_POINT

    double temp;

    temp = b;

    addend(temp);

    FN_EXIT_POINT;
}

void fMATRIX::addend(const c_double &b, const fVECTOR &x)
{
    FN_ENTRY_POINT

    double temp;

    temp = b;

    addend(temp,x);

    FN_EXIT_POINT;
}

void fMATRIX::addend(const c_double &b, const fVECTOR &x, const fVECTOR &y)
{
    FN_ENTRY_POINT

    double temp;

    temp = b;

    addend(temp,x,y);

    FN_EXIT_POINT;
}

void fMATRIX::addstart(void)
{
    FN_ENTRY_POINT

    addend();
    bswap(size,1);

    FN_EXIT_POINT;
}

void fMATRIX::addstart(const   double &b)
{
    FN_ENTRY_POINT

    addend(b);
    bswap(size,1);

    FN_EXIT_POINT;
}

void fMATRIX::addstart(const   double &b, const fVECTOR &x)
{
    FN_ENTRY_POINT

    fVECTOR tempx;

    tempx = x;
    tempx.fswap(1,size+1);
    addend(b,tempx);
    bswap(size,1);

    FN_EXIT_POINT;
}

void fMATRIX::addstart(const   double &b, const fVECTOR &x, const fVECTOR &y)
{
    FN_ENTRY_POINT

    fVECTOR tempx;
    fVECTOR tempy;

    tempx = x;
    tempy = y;
    tempx.fswap(1,size+1);
    tempy.fswap(1,size+1);
    addend(b,tempx,tempy);
    bswap(size,1);

    FN_EXIT_POINT;
}

void fMATRIX::addstart(const c_double &b)
{
    FN_ENTRY_POINT

    addend(b);
    bswap(size,1);

    FN_EXIT_POINT;
}

void fMATRIX::addstart(const c_double &b, const fVECTOR &x)
{
    FN_ENTRY_POINT

    fVECTOR tempx;

    tempx = x;
    tempx.fswap(1,size+1);
    addend(b,tempx);
    bswap(size,1);

    FN_EXIT_POINT;
}

void fMATRIX::addstart(const c_double &b, const fVECTOR &x, const fVECTOR &y)
{
    FN_ENTRY_POINT

    fVECTOR tempx;
    fVECTOR tempy;

    tempx = x;
    tempy = y;
    tempx.fswap(1,size+1);
    tempy.fswap(1,size+1);
    addend(b,tempx,tempy);
    bswap(size,1);

    FN_EXIT_POINT;
}

void fMATRIX::remove(long i)
{
    FN_ENTRY_POINT

    l_xMatrix *temp;

    switch ( m_type )                                                   
    {                                                                   
        case ASYMMETRIC_MATRIX:                                         
        case SYMMETRIC_MATRIX:                                          
        case LOWER_TRIANGULAR_MATRIX:                                   
        case UPPER_TRIANGULAR_MATRIX:                                   
        case DIAGONAL_MATRIX:                                           
        {
            THROW_ASSERT(i >= 1);
            THROW_ASSERT(i <= size);
                                                                        
            fswap(i,size);                                              
                                                                        
            temp = GET_XMATRIX_N(size);                                 
                                                                        
            if ( i < offset_start_effective_col )                       
            {                                                           
                offset_start_effective_col--;                           
                offset_end_effective_col--;                             
            }                                                           
                                                                        
            else                                                        
            {                                                           
                if ( i <= offset_end_effective_col )                    
                {                                                       
                    offset_end_effective_col--;                         
                }                                                       
            }                                                           
                                                                        
            if ( i < offset_start_effective_row )                       
            {                                                           
                offset_start_effective_row--;                           
                offset_end_effective_row--;                             
            }                                                           
                                                                        
            else                                                        
            {                                                           
                if ( i <= offset_end_effective_row )                    
                {                                                       
                    offset_end_effective_row--;                         
                }                                                       
            }

            matrix_struct_lookup.remove(size);
                                                                        
            size--;                                                     
                                                                        
            if ( size == 0 )                                            
            {                                                           
                all_it = NULL;                                          
                                                                        
                reset_offsets();                                        
            }                                                           
                                                                        
            else                                                        
            {                                                           
                (temp->prev)->next = NULL;                              
            }                                                           
                                                                        
            switch ( m_type )                                           
            {                                                           
                case ASYMMETRIC_MATRIX:                                 
                {                                                       
                    REPDELB((temp->row));                               
                    REPDELB((temp->col));

                    temp->row = NULL;
                    temp->col = NULL;
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case SYMMETRIC_MATRIX:                                  
                case LOWER_TRIANGULAR_MATRIX:                           
                {                                                       
                    REPDELB((temp->row));

                    temp->row = NULL;
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case UPPER_TRIANGULAR_MATRIX:                           
                {                                                       
                    REPDELB((temp->col));

                    temp->col = NULL;
                                                                        
                    break;                                              
                }                                                       
                                                                        
                case DIAGONAL_MATRIX:                                   
                {                                                       
                    break;                                              
                }                                                       
                                                                        
                case SCALAR_MATRIX:                                     
                case ZERO_MATRIX:                                       
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
                                                                        
            REPDEL(temp);

            temp = NULL;
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case SCALAR_MATRIX:                                             
        case ZERO_MATRIX:                                               
        {                                                               
            break;                                                      
        }                                                               
                                                                        
        default:                                                        
        {                                                               
            L_THROW(0);                                                 
                                                                        
            break;                                                      
        }                                                               
    }                                                                   

    FN_EXIT_POINT;
}




//
// Misc:
//

void fMATRIX::pad_matrix(long _size)
{
    FN_ENTRY_POINT

    long i;                                                             
                                                                        
    if ( _size == 0 )                                                   
    {                                                                   
        FN_EXIT_POINT;                                                         
    }

    THROW_ASSERT(_size >= 0);
                                                                        
    switch ( get_m_type() )                                             
    {                                                                   
        case ZERO_MATRIX:                                               
        case SCALAR_MATRIX:                                             
        {                                                               
            make_diagonal(DO_FORCE);                                    
        }                                                               
                                                                        
        case ASYMMETRIC_MATRIX:                                         
        case SYMMETRIC_MATRIX:                                          
        case LOWER_TRIANGULAR_MATRIX:                                   
        case UPPER_TRIANGULAR_MATRIX:                                   
        case DIAGONAL_MATRIX:                                           
        {                                                               
            if ( _size > 0 )                                            
            {                                                           
                for ( i = 1 ; i <= _size ; i++ )                        
                {                                                       
                    addend();                                           
                }                                                       
            }                                                           
                                                                        
            break;                                                      
        }                                                               
                                                                        
        default:                                                        
        {                                                               
            L_THROW(0);                                                 
                                                                        
            break;                                                      
        }                                                               
    }                                                                   

    FN_EXIT_POINT;
}

void fMATRIX::invert_matrix(double zt)
{
    FN_ENTRY_POINT

    fMATRIX tempa;
    fMATRIX tempb('r',SCALAR_MATRIX_SIZE,SCALAR_MATRIX);

    long temp_off_start_row;
    long temp_off_start_col;
    long temp_off_end_row;
    long temp_off_end_col;

    long i;

    if ( size > 0 )
    {
        switch ( m_type )
        {
            case SCALAR_MATRIX:
            {
                break;
            }

            case ZERO_MATRIX:
            {
                throw INV_ERR;

                break;
            }

            case ASYMMETRIC_MATRIX:
            {
                THROW_ASSERT(!is_not_eff_square());

                temp_off_start_row = get_offset_start_at_row();
                temp_off_start_col = get_offset_start_at_col();

                temp_off_end_row = get_offset_end_at_row();
                temp_off_end_col = get_offset_end_at_col();

                reset_offsets();

                tempa = (*this);
                tempb(-1,-1) = 1.0;

                (*this)(ARG_INVERT_LEFT,tempa,tempb,zt);

                (*this) = tempa;

                set_offsets(temp_off_start_row,temp_off_start_col,temp_off_end_row,temp_off_end_col);

                break;
            }

            case SYMMETRIC_MATRIX:
            {
                THROW_ASSERT(!is_not_eff_square());

                temp_off_start_row = get_offset_start_at_row();
                temp_off_start_col = get_offset_start_at_col();

                temp_off_end_row = get_offset_end_at_row();
                temp_off_end_col = get_offset_end_at_col();

                reset_offsets();

                tempa = (*this);
                tempb(-1,-1) = 1.0;

                (*this)(ARG_INVERT_LEFT,tempa,tempb,zt);

                (*this) = tempa;

                make_symmetric(DO_FORCE);

                set_offsets(temp_off_start_row,temp_off_start_col,temp_off_end_row,temp_off_end_col);

                break;
            }

            case LOWER_TRIANGULAR_MATRIX:
            {
                THROW_ASSERT(!is_not_eff_square());

                temp_off_start_row = get_offset_start_at_row();
                temp_off_start_col = get_offset_start_at_col();

                temp_off_end_row = get_offset_end_at_row();
                temp_off_end_col = get_offset_end_at_col();

                reset_offsets();

                tempa = (*this);
                tempb(-1,-1) = 1.0;

                (*this)(ARG_INVERT_LEFT,tempa,tempb,zt);

                (*this) = tempa;

                make_lower_triangular(DO_FORCE);

                set_offsets(temp_off_start_row,temp_off_start_col,temp_off_end_row,temp_off_end_col);

                break;
            }

            case UPPER_TRIANGULAR_MATRIX:
            {
                THROW_ASSERT(!is_not_eff_square());

                temp_off_start_row = get_offset_start_at_row();
                temp_off_start_col = get_offset_start_at_col();

                temp_off_end_row = get_offset_end_at_row();
                temp_off_end_col = get_offset_end_at_col();

                reset_offsets();

                tempa = (*this);
                tempb(-1,-1) = 1.0;

                (*this)(ARG_INVERT_LEFT,tempa,tempb,zt);

                (*this) = tempa;

                make_upper_triangular(DO_FORCE);

                set_offsets(temp_off_start_row,temp_off_start_col,temp_off_end_row,temp_off_end_col);

                break;
            }

            case DIAGONAL_MATRIX:
            {
                THROW_ASSERT(!is_not_eff_square());

                temp_off_start_row = get_offset_start_at_row();
                temp_off_start_col = get_offset_start_at_col();

                temp_off_end_row = get_offset_end_at_row();
                temp_off_end_col = get_offset_end_at_col();

                reset_offsets();

                for ( i = 1 ; i <= size ; i++ )
                {
                    NZ_ASSERT((*this)(i-1,i-1),zt);

                    (*this)(i-1,i-1) = ( 1.0 / (*this)(i-1,i-1) );
                }

                set_offsets(temp_off_start_row,temp_off_start_col,temp_off_end_row,temp_off_end_col);

                break;
            }

            default:
            {
                L_THROW(0);

                break;
            }
        }
    }

    FN_EXIT_POINT;
}




//
// Information about effective matrix:
//

int fMATRIX::get_e_type(void) const
{
    FN_ENTRY_POINT

    int result = 0;

    if ( ( get_effective_height() > 0 ) &&                              
         ( get_effective_width() > 0 )    )                             
    {                                                                   
        switch ( m_type )                                               
        {                                                               
            case ASYMMETRIC_MATRIX:                                     
            {                                                           
                result = E_ASYMM;                                       
                                                                        
                break;                                                  
            }                                                           
                                                                        
            case SYMMETRIC_MATRIX:                                      
            {                                                           
                if ( is_not_eff_square() )                              
                {                                                       
                    result = E_ASYMM;                                   
                }                                                       
                                                                        
                else                                                    
                {                                                       
                    if ( get_offset_start_at_row() != get_offset_start_at_col() ) 
                    {                                                   
                        result = E_ASYMM;                               
                    }                                                   
                                                                        
                    else                                                
                    {                                                   
                        result = E_SYMM;                                
                    }                                                   
                }                                                       
                                                                        
                break;                                                  
            }                                                           
                                                                        
            case LOWER_TRIANGULAR_MATRIX:                               
            {                                                           
                if ( get_offset_end_at_row() < get_offset_start_at_col() )    
                {                                                       
                    result = E_ZERO;                                    
                }                                                       
                                                                        
                else                                                    
                {                                                       
                    if ( get_offset_start_at_row() < get_offset_end_at_col() ) 
                    {                                                   
                        result = E_LOWER;                               
                    }                                                   
                                                                        
                    else                                                
                    {                                                   
                        result = E_ASYMM;                               
                    }                                                   
                }                                                       
                                                                        
                break;                                                  
            }                                                           
                                                                        
            case UPPER_TRIANGULAR_MATRIX:                               
            {                                                           
                if ( get_offset_start_at_row() > get_offset_end_at_col() )    
                {                                                       
                    result = E_ZERO;                                    
                }                                                       
                                                                        
                else                                                    
                {                                                       
                    if ( get_offset_end_at_row() > get_offset_start_at_col() ) 
                    {                                                   
                        result = E_UPPER;                               
                    }                                                   
                                                                        
                    else                                                
                    {                                                   
                        result = E_ASYMM;                               
                    }                                                   
                }                                                       
                                                                        
                break;                                                  
            }                                                           
                                                                        
            case DIAGONAL_MATRIX:                                       
            {                                                           
                if ( ( get_offset_end_at_row() < get_offset_start_at_col() ) ||   
                     ( get_offset_start_at_row() > get_offset_end_at_col() )    ) 
                {                                                       
                    result = E_ZERO;                                    
                }                                                       
                                                                        
                else                                                    
                {                                                       
                    if ( ( get_effective_height() > 1 ) ||              
                         ( get_effective_width() > 1 )    )             
                    {                                                   
                        result = E_DIAG;                                
                    }                                                   
                                                                        
                    else                                                
                    {                                                   
                        result = E_ASYMM;                               
                    }                                                   
                }                                                       
                                                                        
                break;                                                  
            }                                                           
                                                                        
            case ZERO_MATRIX:                                           
            {                                                           
                result = E_VZ;                                          
                                                                        
                break;                                                  
            }                                                           
                                                                        
            case SCALAR_MATRIX:                                         
            {                                                           
                result = E_VS;                                          
                                                                        
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
        switch ( m_type )                                               
        {                                                               
            case ASYMMETRIC_MATRIX:                                     
            case SYMMETRIC_MATRIX:                                      
            case LOWER_TRIANGULAR_MATRIX:                               
            case UPPER_TRIANGULAR_MATRIX:                               
            case DIAGONAL_MATRIX:                                       
            {                                                           
                result = E_EMPTY;                                       
                                                                        
                break;                                                  
            }                                                           
                                                                        
            case ZERO_MATRIX:                                           
            {                                                           
                result = E_VZ;                                          
                                                                        
                break;                                                  
            }                                                           
                                                                        
            case SCALAR_MATRIX:                                         
            {                                                           
                result = E_VS;                                          
                                                                        
                break;                                                  
            }                                                           
                                                                        
            default:                                                    
            {                                                           
                L_THROW(0);                                             
                                                                        
                break;                                                  
            }                                                           
        }                                                               
    }                                                                   

    FN_EXIT_POINT result;
}

long fMATRIX::get_e_offs(void) const
{
    FN_ENTRY_POINT

    long result = 0;

    switch ( get_e_type() )                                             
    {                                                                   
        case E_EMPTY:                                                   
        case E_VZ:                                                      
        case E_VS:                                                      
        {                                                               
            result = 0;                                                 
                                                                        
            break;                                                      
        }                                                               
                                                                        
        case E_ZERO:                                                    
        case E_SYMM:                                                    
        case E_ASYMM:                                                   
        case E_DIAG:                                                    
        case E_LOWER:                                                   
        case E_UPPER:                                                   
        {                                                               
            result = get_offset_start_at_col() - get_offset_start_at_row();   
                                                                        
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

void fMATRIX::fix_symm(void)
{
    FN_ENTRY_POINT

    switch ( get_e_type() )
    {
        case E_ASYMM:
        {
            if ( is_symmetric() )
            {
                make_asymmetric();
            }

            break;
        }

        case E_EMPTY:
        case E_ZERO:
        case E_SYMM:
        case E_DIAG:
        case E_LOWER:
        case E_UPPER:
        case E_VZ:
        case E_VS:
        {
            break;
        }

        default:
        {
            L_THROW(0);

            break;
        }
    }

    FN_EXIT_POINT;
}




//
// IO stuff
//

void fMATRIX::prime_io(std::ostream *_where_to, int _echo_level)
{
    FN_ENTRY_POINT

    where_to   = _where_to;
    echo_level = _echo_level;

    FN_EXIT_POINT;
}




//
// Raw data operations:
//

L_DOUBLE fMATRIX::getelm(long i, long j) const
{
    FN_ENTRY_POINT

    if ( ( i >= 1 ) && ( j >= 1 ) )                                     
    {                                                                   
        if ( ( m_type == SCALAR_MATRIX ) || ( m_type == ZERO_MATRIX ) ) 
        {                                                               
            if ( ( i == j ) && ( m_type == SCALAR_MATRIX ) )            
            {                                                           
                FN_EXIT_POINT m_diag;                                          
            }                                                           
                                                                        
            else                                                        
            {                                                           
                FN_EXIT_POINT 0.0;                                        
            }                                                           
        }                                                               
                                                                        
        else                                                            
        {                                                               
            THROW_ASSERT(i >= 1);                                             
            THROW_ASSERT(i <= size);                                          
            THROW_ASSERT(j >= 1);                                             
            THROW_ASSERT(j <= size);                                          
                                                                        
            if ( i == j )                                               
            {                                                           
                FN_EXIT_POINT GET_MATRIX_N(i)->diag;
            }                                                           
                                                                        
            else                                                        
            {                                                           
                if ( i > j )                                            
                {                                                       
                    switch ( m_type )                                   
                    {                                                   
                        case ASYMMETRIC_MATRIX:                         
                        case SYMMETRIC_MATRIX:                          
                        case LOWER_TRIANGULAR_MATRIX:                   
                        {                                               
                            FN_EXIT_POINT (GET_MATRIX_N(i)->row)[j-1];        
                                                                        
                            break;                                      
                        }                                               
                                                                        
                        case UPPER_TRIANGULAR_MATRIX:                   
                        case DIAGONAL_MATRIX:                           
                        {                                               
                            FN_EXIT_POINT 0.0;                            
                                                                        
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
                    switch ( m_type )                                   
                    {                                                   
                        case ASYMMETRIC_MATRIX:                         
                        case UPPER_TRIANGULAR_MATRIX:                   
                        case SYMMETRIC_MATRIX:                          
                        {                                               
                            FN_EXIT_POINT (GET_MATRIX_N(j)->col)[i-1];                    
                                                                        
                            break;                                      
                        }                                               
                                                                        
                        case LOWER_TRIANGULAR_MATRIX:                   
                        case DIAGONAL_MATRIX:                           
                        {                                               
                            FN_EXIT_POINT 0.0;                            
                                                                        
                            break;                                      
                        }                                               
                                                                        
                        default:                                        
                        {                                               
                            L_THROW(0);                                 
                                                                        
                            break;                                      
                        }                                               
                    }                                                   
                }                                                       
            }                                                           
        }                                                               
    }

    THROW_ASSERT((i == -1) && (j == -1));

    FN_EXIT_POINT m_diag;
}

L_DOUBLE &fMATRIX::getref(long i, long j, int slack)
{
    FN_ENTRY_POINT

    l_xMatrix *temp;

    if ( ( i >= 1 ) && ( j >= 1 ) )                                     
    {                                                                   
        if ( ( m_type == SCALAR_MATRIX ) || ( m_type == ZERO_MATRIX ) ) 
        {                                                               
            THROW_ASSERT(slack);

            if ( ( i == j ) && ( m_type == SCALAR_MATRIX ) )        
            {                                                       
                FN_EXIT_POINT l_double_get_static_ref(m_diag);                               
            }                                                       
                                                                        
            else                                                    
            {                                                       
                FN_EXIT_POINT l_double_get_static_zero_ref();                                 
            }                                                       
        }                                                               
                                                                        
        else                                                            
        {                                                               
            THROW_ASSERT(i >= 1);                                             
            THROW_ASSERT(i <= size);                                          
            THROW_ASSERT(j >= 1);                                             
            THROW_ASSERT(j <= size);                                          
                                                                        
            if ( i == j )                                               
            {                                                           
                FN_EXIT_POINT GET_XMATRIX_N(i)->diag;                          
            }                                                           
                                                                        
            else                                                        
            {                                                           
                if ( i > j )                                            
                {                                                       
                    switch ( m_type )                                   
                    {                                                   
                        case ASYMMETRIC_MATRIX:                         
                        case SYMMETRIC_MATRIX:                          
                        case LOWER_TRIANGULAR_MATRIX:                   
                        {                                               
                            FN_EXIT_POINT (GET_XMATRIX_N(i)->row)[j-1];        
                                                                        
                            break;                                      
                        }                                               
                                                                        
                        case UPPER_TRIANGULAR_MATRIX:                   
                        {                                               
                            make_asymmetric();                          
                                                                        
                            FN_EXIT_POINT getref(i,j);                         
                                                                        
                            break;                                      
                        }                                               
                                                                        
                        case DIAGONAL_MATRIX:                           
                        {                                               
                            make_lower_triangular();                    
                                                                        
                            FN_EXIT_POINT getref(i,j);                         
                                                                        
                            break;                                      
                        }                                               
                                                                        
                                                                        
                        default:                                        
                        {                                               
                            L_THROW(1);                                 
                                                                        
                            break;                                      
                        }                                               
                    }                                                   
                }                                                       
                                                                        
                else                                                    
                {                                                       
                    switch ( m_type )                                   
                    {                                                   
                        case ASYMMETRIC_MATRIX:                         
                        case UPPER_TRIANGULAR_MATRIX:                   
                        {                                               
                            FN_EXIT_POINT (GET_XMATRIX_N(j)->col)[i-1];        
                                                                        
                            break;                                      
                        }                                               
                                                                        
                        case SYMMETRIC_MATRIX:                          
                        {                                               
                            temp = GET_XMATRIX_N(j);                    
                                                                        
                            FN_EXIT_POINT (temp->col)[i-1];                    
                                                                        
                            break;                                      
                        }                                               
                                                                        
                        case LOWER_TRIANGULAR_MATRIX:                   
                        {                                               
                            make_asymmetric();                          
                                                                        
                            FN_EXIT_POINT getref(i,j);                         
                                                                        
                            break;                                      
                        }                                               
                                                                        
                        case DIAGONAL_MATRIX:                           
                        {                                               
                            make_upper_triangular();                    
                                                                        
                            FN_EXIT_POINT getref(i,j);                         
                                                                        
                            break;                                      
                        }                                               
                                                                        
                        default:                                        
                        {                                               
                            L_THROW(2);                                 
                                                                        
                            break;                                      
                        }                                               
                    }                                                   
                }                                                       
            }                                                           
        }                                                               
    }

    THROW_ASSERT((i == -1) && (j == -1));

    FN_EXIT_POINT m_diag;

    slack = 0;
}

void fMATRIX::setelm(const   double &x, long i, long j)
{
    FN_ENTRY_POINT

    getref(i,j) = x;

    FN_EXIT_POINT;
}

void fMATRIX::setelm(const c_double &x, long i, long j)
{
    FN_ENTRY_POINT

    getref(i,j) = x;

    FN_EXIT_POINT;
}

