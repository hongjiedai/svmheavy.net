
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
// Indexed matrix factorisation class.
//
// Written by: Alistair Shilton
//             Melbourne University
//




#include <iostream>
#include <string.h>

#include "factor.h"
#include "search.h"
#include "svdefs.h"
#include "c_double.h"
#include "vector.h"
#include "matrix.h"
#include "friends.h"
#include "outfilt.h"


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

std::ostream &operator<<(std::ostream &output, const fMATRFact &source)
{
    FN_ENTRY_POINT

    output << "nbad: " << source.nbad << "\n\n";

    output << "Primary swap flag: " << source.primary_swap_flag << "\n";
    output << "Case four flag:    " << source.case_four_flag    << "\n\n";

    output << "s_valid: " << source.s_valid << "\n";
    output << "s:       " << source.s       << "\n\n";

    output << "L: " << source.L << "\n";

    output << "Zero tolerance: " << source.zt << "\n\n";

    output << "Indexing: " << source.use_index << "\n";

    output << "f_form: " << source.f_form << "\n";
    output << "f_type: " << source.f_type << "\n\n";

    FN_EXIT_POINT output;
}


std::istream &operator>>(std::istream &input, fMATRFact &dest)
{
    FN_ENTRY_POINT

    wait_dummy howzat;

    input >> howzat; input >> dest.nbad;

    input >> howzat; input >> dest.primary_swap_flag;
    input >> howzat; input >> dest.case_four_flag;

    input >> howzat; input >> dest.s_valid;
    input >> howzat; input >> dest.s;

    input >> howzat; input >> dest.L;

    input >> howzat; input >> dest.zt;

    input >> howzat; input >> dest.use_index;

    input >> howzat; input >> dest.f_form;
    input >> howzat; input >> dest.f_type;

    FN_EXIT_POINT input;
}





//
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=                                                                +=+=
// +=+=              MATRIX FACTORISATION CLASSES                      +=+=
// +=+=                                                                +=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
//

//
// Constuctors:
//

fMATRFact::fMATRFact()
{
    FN_ENTRY_POINT

    L.make_zero(DO_FORCE);
    L.make_symmetric(DO_FORCE);

    f_form = FORM_ONE;
    f_type = DEFAULT_FACT;

    case_four_flag    = 0;
    primary_swap_flag = 0;

    zt = 0;

    fact_none_size = 0;

    if ( f_type == FACT_NONE )
    {
        nbad = -1;
    }

    else
    {
        nbad = 0;

        fix_symmetry();
    }

    FN_EXIT_POINT;
}

fMATRFact::fMATRFact(int _use_index, const lVECTOR &Gindex, const fMATRIX &G_real, double _zt, int _f_type)
{
    FN_ENTRY_POINT

    use_index = _use_index;

    THROW_ASSERT(!(G_real.is_not_eff_square()));

    f_form = FORM_ONE;
    f_type = _f_type;

    case_four_flag    = 0;
    primary_swap_flag = 0;

    zt = _zt;

    fact_none_size = 0;

    if ( f_type == FACT_NONE )
    {
        if ( use_index )
        {
            fact_none_size = Gindex.get_effective_size();
        }

        else
        {
            fact_none_size = G_real.get_effective_height();
        }

        nbad = -1;
    }

    else
    {
        nbad = 0;

        fix_symmetry();

        L.make_zero(DO_FORCE);
        L.make_symmetric(DO_FORCE);
        L.pad_matrix(Gindex.get_effective_size());
        overwrite_L(Gindex,G_real);

        refresh_fact();
        xfact(Gindex,G_real);
        fix_symmetry();
    }

    FN_EXIT_POINT;
}

fMATRFact::fMATRFact(int _use_index, const lVECTOR &Gindex, const fMATRIX &G_real, const fVECTOR &d_real, double _zt, int _f_type)
{
    FN_ENTRY_POINT

    use_index = _use_index;

    THROW_ASSERT(!(G_real.is_not_eff_square()));

    f_form = FORM_TWO;
    f_type = _f_type;

    case_four_flag    = 0;
    primary_swap_flag = 0;

    zt = _zt;

    fact_none_size = 0;

    if ( f_type == FACT_NONE )
    {
        if ( use_index )
        {
            fact_none_size = Gindex.get_effective_size();
        }

        else
        {
            fact_none_size = G_real.get_effective_height();
        }

        nbad = -1;
    }

    else
    {
        nbad = 0;

        L.make_zero(DO_FORCE);
        L.make_symmetric(DO_FORCE);
        L.pad_matrix(Gindex.get_effective_size());
        overwrite_L(Gindex,G_real);

        refresh_fact(d_real);
        xfact(Gindex,G_real,d_real);
        fix_symmetry();
    }

    FN_EXIT_POINT;
}

fMATRFact::fMATRFact(const fMATRFact &source)
{
    FN_ENTRY_POINT

    use_index = source.use_index;

    L.make_zero(DO_FORCE);
    s.make_zero(DO_FORCE);

    case_four_flag    = source.case_four_flag;
    primary_swap_flag = source.primary_swap_flag;

    s_valid = source.s_valid;

    L = source.L;
    s = source.s;

    nbad = source.nbad;

    zt = source.zt;

    f_form = source.f_form;
    f_type = source.f_type;

    fact_none_size = source.fact_none_size;

    if ( f_type != FACT_NONE )
    {
        fix_symmetry();
    }

    FN_EXIT_POINT;
}




//
// Overwrite assignment operator
//

fMATRFact &fMATRFact::operator=(const fMATRFact &source)
{
    FN_ENTRY_POINT

    use_index = source.use_index;

    s.make_zero(DO_FORCE);
    L.make_zero(DO_FORCE);

    nbad = source.nbad;

    case_four_flag    = source.case_four_flag;
    primary_swap_flag = source.primary_swap_flag;

    s_valid = source.s_valid;
    s       = source.s;

    L = source.L;

    zt = source.zt;

    f_form = source.f_form;
    f_type = source.f_type;

    fact_none_size = source.fact_none_size;

    if ( f_type != FACT_NONE )
    {
        fix_symmetry();
    }

    FN_EXIT_POINT (*this);
}




//
// Rank-1 updates:
//

long fMATRFact::rankone(const lVECTOR &Gindex, const fVECTOR &ax, const L_DOUBLE &b, const fMATRIX &G_real, long z_start, long z_end)
{
    FN_ENTRY_POINT

    long result;

    THROW_ASSERT((z_start+z_end) <= get_size());
    THROW_ASSERT((ax.get_effective_size() < 0) || (ax.get_effective_size() == get_size()));

    if ( f_type == FACT_NONE )
    {
        result = -1;
    }

    else
    {
        if ( !( ( ax.is_zero()                ) ||
                ( b == 0.0                    ) ||
                ( z_start+z_end == get_size() )    ) )
        {
            switch ( f_form )
            {
                case FORM_ONE:
                {
                    fVECTOR a('x',get_size());

                    a = ax;

                    enter_fixvect_normal_fn(a);
                    xrankone(Gindex,a,b,G_real,0,z_start,z_end);

                    break;
                }

                case FORM_TWO:
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

        fix_symmetry();

        result = xfact(Gindex,G_real);
    }

    FN_EXIT_POINT result;
}

long fMATRFact::rankone(const lVECTOR &Gindex, const fVECTOR &ax, const L_DOUBLE &b, const fMATRIX &G_real, const fVECTOR &d_real, long z_start, long z_end)
{
    FN_ENTRY_POINT

    long result;

    if ( f_type == FACT_NONE )
    {
        result = -1;
    }

    else
    {
        switch ( f_form )
        {
            case FORM_ONE:
            {
                rankone(Gindex,ax,b,G_real,z_start,z_end);

                break;
            }

            case FORM_TWO:
            {
                THROW_ASSERT((z_start+z_end) <= get_size());
                THROW_ASSERT((ax.get_effective_size() < 0) || (ax.get_effective_size() == get_size()));

                if ( !( ( ax.is_zero()                ) ||
                        ( b == 0.0                    ) ||
                        ( z_start+z_end == get_size() )    ) )
                {
                    fVECTOR a('x',get_size());

                    a = ax;

                    a.addstart();

                    if ( primary_swap_flag )
                    {
                        switch ( z_start )
                        {
                            case 0:
                            {
                                if ( z_end == get_size()-1 )
                                {
                                    z_end--;
                                    z_start++;
                                }

                                break;
                            }

                            case 1:
                            {
                                if ( z_end == get_size()-2 )
                                {
                                    z_end++;
                                    z_start--;
                                }

                                else
                                {
                                    z_start--;
                                }

                                break;
                            }

                            default:
                            {
                                break;
                            }
                        }
                    }

                    enter_fixvect_normal_fn(a);
                    xrankone(Gindex,a,b,G_real,d_real,0,z_start,z_end);
                }

                break;
            }

            default:
            {
                L_THROW(4);

                break;
            }
        }

        fix_symmetry();

        result = xfact(Gindex,G_real,d_real);
    }

    FN_EXIT_POINT result;
}

long fMATRFact::pert_one(const lVECTOR &Gindex, long i, const L_DOUBLE &b, const fMATRIX &G_real)
{
    FN_ENTRY_POINT

    long result;

    if ( f_type == FACT_NONE )
    {
        result = -1;
    }

    else
    {
        switch ( f_form )
        {
            case FORM_ONE:
            {
                fVECTOR a('x',SELECT_VECTOR_SIZE);

                a.set_sel_elm(i);

                rankone(Gindex,a,b,G_real,i-1,get_size()-i);

                break;
            }

            case FORM_TWO:
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

        fix_symmetry();

        result = nbad;
    }

    FN_EXIT_POINT result;
}

long fMATRFact::pert_one(const lVECTOR &Gindex, long i, const L_DOUBLE &b, const fMATRIX &G_real, const fVECTOR &d_real)
{
    long result;

    if ( f_type == FACT_NONE )
    {
        result = -1;
    }

    else
    {
        switch ( f_form )
        {
            case FORM_ONE:
            {
                pert_one(Gindex,i,b,G_real);

                break;
            }

            case FORM_TWO:
            {
                fVECTOR a('x',SELECT_VECTOR_SIZE);

                a.set_sel_elm(i);

                rankone(Gindex,a,b,G_real,d_real,i-1,get_size()-i);

                break;
            }

            default:
            {
                L_THROW(0);

                break;
            }
        }

        fix_symmetry();

        result = nbad;
    }

    FN_EXIT_POINT result;
}




//
// Inversion:
//

long fMATRFact::minverse(const lVECTOR &Gindex, fVECTOR &ax, const fVECTOR &bx, const fMATRIX &G_real, long z_start, long z_end)
{
    FN_ENTRY_POINT

    THROW_ASSERT((z_start+z_end) <= (get_size()-nbad));
    THROW_ASSERT(ax.get_effective_size() == (get_size()-nbad));
    THROW_ASSERT((bx.get_effective_size() < 0) || (bx.get_effective_size() == (get_size()-nbad)));

    long result;

    if ( f_type == FACT_NONE )
    {
        fMATRFact temp(use_index,Gindex,G_real,zt,DEFAULT_FACT);

        temp.minverse(Gindex,ax,bx,G_real,z_start,z_end);

        result = -1;
    }

    else
    {
        switch ( f_form )
        {
            case FORM_ONE:
            {
                fVECTOR a('x',get_size()-nbad);
                fVECTOR b('x',get_size()-nbad);

                b = bx;

                enter_fixvect_normal_fn(b);

                if ( get_size()-nbad > 0 )
                {
                    if ( z_start+z_end < get_size()-nbad )
                    {
                        switch ( f_type )
                        {
                            case FACT_INVERSE:
                            {
                                b.set_offset_start(1+z_start);
                                b.set_offset_end(get_size()-nbad-z_end);

                                L.set_offset_start_at_row(1);
                                L.set_offset_start_at_col(1+z_start);
                                L.set_offset_end_at_row(get_size()-nbad);
                                L.set_offset_end_at_col(get_size()-nbad-z_end);

                                a = ( L * b );

                                L.reset_offsets();
                                b.reset_offsets();

                                break;
                            }

                            case FACT_CHOL:
                            {
                                fVECTOR cx('x',get_size()-nbad);

                                {
                                    if ( z_start >= 1 )
                                    {
                                        cx.set_offset_start(1);
                                        cx.set_offset_end(z_start);

                                        cx = 0.0;

                                        cx.reset_offsets();
                                    }

                                    b.set_offset_start(1+z_start);
                                    b.set_offset_end(get_size()-nbad-z_end);

                                    cx.set_offset_start(1+z_start);
                                    cx.set_offset_end(get_size()-nbad-z_end);

                                    L.set_offset_start_at_row(1+z_start);
                                    L.set_offset_start_at_col(1+z_start);
                                    L.set_offset_end_at_row(get_size()-nbad-z_end);
                                    L.set_offset_end_at_col(get_size()-nbad-z_end);

                                    L(ARG_INVERT_LEFT,cx,b);

                                    L.reset_offsets();
                                    cx.reset_offsets();
                                    b.reset_offsets();

                                    if ( z_end >= 1 )
                                    {
                                        fVECTOR dx('x',z_end);

                                        b.set_offset_start(1+z_start);
                                        b.set_offset_end(get_size()-nbad-z_end);

                                        L.set_offset_start_at_row(get_size()-nbad-z_end+1);
                                        L.set_offset_start_at_col(1+z_start);
                                        L.set_offset_end_at_row(get_size()-nbad);
                                        L.set_offset_end_at_col(get_size()-nbad-z_end);

                                        dx = -( L * b );

                                        L.reset_offsets();
                                        b.reset_offsets();

                                        cx.set_offset_start(get_size()-nbad-z_end+1);
                                        cx.set_offset_end(get_size()-nbad);

                                        L.set_offset_start_at_row(get_size()-nbad-z_end+1);
                                        L.set_offset_start_at_col(get_size()-nbad-z_end+1);
                                        L.set_offset_end_at_row(get_size()-nbad);
                                        L.set_offset_end_at_col(get_size()-nbad);

                                        L(ARG_INVERT_LEFT,cx,dx);

                                        L.reset_offsets();
                                        cx.reset_offsets();
                                    }
                                }

                                {
                                    cx.set_offset_start(1+z_start);
                                    cx.set_offset_end(get_size()-nbad);

                                    a.set_offset_start(1+z_start);
                                    a.set_offset_end(get_size()-nbad);

                                    L.set_offset_start_at_row(1+z_start);
                                    L.set_offset_start_at_col(1+z_start);
                                    L.set_offset_end_at_row(get_size()-nbad);
                                    L.set_offset_end_at_col(get_size()-nbad);

                                    L(ARG_INVERT_RIGHT,a,cx);

                                    L.reset_offsets();
                                    a.reset_offsets();
                                    cx.reset_offsets();

                                    if ( z_start >= 1 )
                                    {
                                        fVECTOR dx('x',z_start);

                                        a.set_offset_start(1+z_start);
                                        a.set_offset_end(get_size()-nbad);

                                        L.set_offset_start_at_row(1+z_start);
                                        L.set_offset_start_at_col(1);
                                        L.set_offset_end_at_row(get_size()-nbad);
                                        L.set_offset_end_at_col(z_start);

                                        dx = -( a * L );

                                        L.reset_offsets();
                                        a.reset_offsets();

                                        a.set_offset_start(1);
                                        a.set_offset_end(z_start);

                                        L.set_offset_start_at_row(1);
                                        L.set_offset_start_at_col(1);
                                        L.set_offset_end_at_row(z_start);
                                        L.set_offset_end_at_col(z_start);

                                        L(ARG_INVERT_RIGHT,a,dx);

                                        L.reset_offsets();
                                        a.reset_offsets();
                                    }
                                }

                                break;
                            }

                            case FACT_NONE:
                            default:
                            {
                                L_THROW(4);

                                break;
                            }
                        }
                    }

                    else
                    {
                        a = 0.0;
                    }
                }

                exit_fixvect_normal_fn(a);

                ax = a;

                break;
            }

            case FORM_TWO:
            {
                L_THROW(0);

                break;
            }

            default:
            {
                L_THROW(6);

                break;
            }
        }

        result = nbad;
    }

    FN_EXIT_POINT result;
}

long fMATRFact::minverse(const lVECTOR &Gindex, fVECTOR &ax, L_DOUBLE &ay, const fVECTOR &bx, const L_DOUBLE &by, const fMATRIX &G_real, const fVECTOR &d_real, long z_start, long z_end, int elm_one_zero)
{
    FN_ENTRY_POINT

    long result;

    if ( f_type == FACT_NONE )
    {
        fMATRFact temp(use_index,Gindex,G_real,d_real,zt,DEFAULT_FACT);

        temp.minverse(Gindex,ax,ay,bx,by,G_real,d_real,z_start,z_end,elm_one_zero);

        result = -1;
    }

    else
    {
        switch ( f_form )
        {
            case FORM_ONE:
            {
                ay = 0.0;

                minverse(Gindex,ax,bx,G_real);

                break;
            }

            case FORM_TWO:
            {
                THROW_ASSERT((z_start+z_end) <= (get_size()-nbad));
                THROW_ASSERT(ax.get_effective_size() == (get_size()-nbad));
                THROW_ASSERT((bx.get_effective_size() < 0) || (bx.get_effective_size() == (get_size()-nbad)));

                fVECTOR a('x',get_size()-nbad+1);
                fVECTOR b('x',get_size()-nbad);

                b = bx;

                b.addstart(by);

                enter_fixvect_normal_fn(b);

                if ( primary_swap_flag )
                {
                    switch ( z_start )
                    {
                        case 0:
                        {
                            if ( z_end == get_size()-1 )
                            {
                                z_end--;
                                z_start++;
                            }

                            break;
                        }

                        case 1:
                        {
                            if ( z_end == get_size()-2 )
                            {
                                z_end++;
                                z_start--;
                            }

                            else
                            {
                                z_start--;
                            }

                            break;
                        }

                        default:
                        {
                            break;
                        }
                    }
                }

                if ( get_size()-nbad+1 > 0 )
                {
                    if ( ( z_start+z_end < get_size()-nbad ) ||
                         ( !(elm_one_zero)                 )    )
                    {
                        switch ( f_type )
                        {
                            case FACT_INVERSE:
                            {
                                if ( !elm_one_zero )
                                {
                                    b.set_offset_start(1);
                                    b.set_offset_end(1);

                                    L.set_offset_start_at_row(1);
                                    L.set_offset_start_at_col(1);
                                    L.set_offset_end_at_row(get_size()-nbad+1);
                                    L.set_offset_end_at_col(1);

                                    a = ( L * b );

                                    L.reset_offsets();
                                    b.reset_offsets();
                                }

                                else
                                {
                                    a = 0.0;
                                }

                                if ( z_start+z_end < get_size()-nbad )
                                {
                                    b.set_offset_start(1+z_start+1);
                                    b.set_offset_end(get_size()-nbad+1-z_end);

                                    L.set_offset_start_at_row(1);
                                    L.set_offset_start_at_col(1+z_start+1);
                                    L.set_offset_end_at_row(get_size()-nbad+1);
                                    L.set_offset_end_at_col(get_size()-nbad+1-z_end);

                                    a += ( L * b );

                                    L.reset_offsets();
                                    b.reset_offsets();
                                }

                                break;
                            }

                            case FACT_CHOL:
                            {
                                if ( case_four_flag )
                                {
                                    L_DOUBLE d_b_conj;

                                    d_b_conj = L(1,0);

                                    if ( elm_one_zero )
                                    {
                                        if ( z_start == 0 )
                                        {
                                            a[0] = 0.0;
                                            a[1] = ( d_b_conj * b[0] );
                                        }

                                        else
                                        {
                                            // We don't need this case, so the
                                            // code will never be reached - but
                                            // it is kept for consistency.

                                            a[0] = 0.0;
                                            a[1] = 0.0;
                                        }
                                    }

                                    else
                                    {
                                        if ( z_start == 0 )
                                        {
                                            a[0] = ( L(1,0)   * b[1] );
                                            a[1] = ( d_b_conj * b[0] );
                                        }

                                        else
                                        {
                                            a[0] = ( L(1,0) * b[1] );
                                            a[1] = 0.0;
                                        }
                                    }
                                }

                                else
                                {
                                    if ( z_start+z_end == get_size()-nbad )
                                    {
                                        // We may assume that elm_one_zero == 0

                                        fVECTOR cx('x',get_size()-nbad+1);

                                        {
                                            cx[0] = 0.0;
                                            cx[1] = ( b[0] / L(1,1) );

                                            if ( get_size()-nbad >= 2 )
                                            {
                                                fVECTOR dx('x',get_size()-nbad-1);

                                                cx.set_offset_start(2);
                                                cx.set_offset_end(2);

                                                L.set_offset_start_at_row(2+1);
                                                L.set_offset_start_at_col(2);
                                                L.set_offset_end_at_row(get_size()-nbad+1);
                                                L.set_offset_end_at_col(2);

                                                dx = ( L * cx );

                                                L.reset_offsets();
                                                cx.reset_offsets();

                                                cx.set_offset_start(2+1);
                                                cx.set_offset_end(get_size()-nbad+1);

                                                L.set_offset_start_at_row(2+1);
                                                L.set_offset_start_at_col(2+1);
                                                L.set_offset_end_at_row(get_size()-nbad+1);
                                                L.set_offset_end_at_col(get_size()-nbad+1);

                                                L(ARG_INVERT_LEFT,cx,dx);

                                                L.reset_offsets();
                                                cx.reset_offsets();
                                            }
                                        }

                                        cx[1] = -cx[1];

                                        {
                                            fVECTOR dx('x',1);

                                            cx.set_offset_start(2);
                                            cx.set_offset_end(get_size()-nbad+1);

                                            a.set_offset_start(2);
                                            a.set_offset_end(get_size()-nbad+1);

                                            L.set_offset_start_at_row(2);
                                            L.set_offset_start_at_col(2);
                                            L.set_offset_end_at_row(get_size()-nbad+1);
                                            L.set_offset_end_at_col(get_size()-nbad+1);

                                            L(ARG_INVERT_RIGHT,a,cx);

                                            L.reset_offsets();
                                            a.reset_offsets();
                                            cx.reset_offsets();

                                            a.set_offset_start(2);
                                            a.set_offset_end(get_size()-nbad+1);

                                            L.set_offset_start_at_row(2);
                                            L.set_offset_start_at_col(1);
                                            L.set_offset_end_at_row(get_size()-nbad+1);
                                            L.set_offset_end_at_col(1);

                                            dx = -( a * L );

                                            L.reset_offsets();
                                            a.reset_offsets();

                                            a.set_offset_start(1);
                                            a.set_offset_end(1);

                                            L.set_offset_start_at_row(1);
                                            L.set_offset_start_at_col(1);
                                            L.set_offset_end_at_row(1);
                                            L.set_offset_end_at_col(1);

                                            L(ARG_INVERT_RIGHT,a,dx);

                                            L.reset_offsets();
                                            a.reset_offsets();
                                        }
                                    }

                                    else
                                    {
                                        if ( z_start == 0 )
                                        {
                                            if ( !(elm_one_zero) )
                                            {
                                                fVECTOR cx('x',get_size()-nbad+1);

                                                {
                                                    b.set_offset_start(1);
                                                    b.set_offset_end(get_size()-nbad+1-z_end);

                                                    cx.set_offset_start(1);
                                                    cx.set_offset_end(get_size()-nbad+1-z_end);

                                                    L.set_offset_start_at_row(1);
                                                    L.set_offset_start_at_col(1);
                                                    L.set_offset_end_at_row(get_size()-nbad+1-z_end);
                                                    L.set_offset_end_at_col(get_size()-nbad+1-z_end);

                                                    L(ARG_INVERT_LEFT,cx,b);

                                                    L.reset_offsets();
                                                    cx.reset_offsets();
                                                    b.reset_offsets();

                                                    if ( z_end > 0 )
                                                    {
                                                        fVECTOR dx('x',z_end);

                                                        cx.set_offset_start(1);
                                                        cx.set_offset_end(get_size()-nbad+1-z_end);

                                                        L.set_offset_start_at_row(get_size()-nbad+1-z_end+1);
                                                        L.set_offset_start_at_col(1);
                                                        L.set_offset_end_at_row(get_size()-nbad+1);
                                                        L.set_offset_end_at_col(get_size()-nbad+1-z_end);

                                                        dx = -( L * cx );

                                                        L.reset_offsets();
                                                        cx.reset_offsets();

                                                        cx.set_offset_start(get_size()-nbad+1-z_end+1);
                                                        cx.set_offset_end(get_size()-nbad+1);

                                                        L.set_offset_start_at_row(get_size()-nbad+1-z_end+1);
                                                        L.set_offset_start_at_col(get_size()-nbad+1-z_end+1);
                                                        L.set_offset_end_at_row(get_size()-nbad+1);
                                                        L.set_offset_end_at_col(get_size()-nbad+1);

                                                        L(ARG_INVERT_LEFT,cx,dx);

                                                        L.reset_offsets();
                                                        cx.reset_offsets();
                                                    }
                                                }

                                                cx[1] = -cx[1];

                                                {
                                                    cx.set_offset_start(1);
                                                    cx.set_offset_end(get_size()-nbad+1);

                                                    a.set_offset_start(1);
                                                    a.set_offset_end(get_size()-nbad+1);

                                                    L.set_offset_start_at_row(1);
                                                    L.set_offset_start_at_col(1);
                                                    L.set_offset_end_at_row(get_size()-nbad+1);
                                                    L.set_offset_end_at_col(get_size()-nbad+1);

                                                    L(ARG_INVERT_RIGHT,a,cx);

                                                    L.reset_offsets();
                                                    a.reset_offsets();
                                                    cx.reset_offsets();
                                                }
                                            }

                                            else
                                            {
                                                fVECTOR cx('x',(get_size()-nbad+1));

                                                {
                                                    cx[0] = (              b[0]   / L(0,0) );
                                                    cx[1] = ( -( L(1,0) * cx[0] ) / L(1,1) );

                                                    if ( z_end < get_size()-nbad-1 )
                                                    {
                                                        fVECTOR dx('x',(get_size()-nbad-1-z_end));

                                                        b.set_offset_start(3);
                                                        b.set_offset_end(get_size()-nbad+1-z_end);

                                                        dx = b;

                                                        b.reset_offsets();

                                                        cx.set_offset_start(1);
                                                        cx.set_offset_end(2);

                                                        L.set_offset_start_at_row(3);
                                                        L.set_offset_start_at_col(1);
                                                        L.set_offset_end_at_row(get_size()-nbad+1-z_end);
                                                        L.set_offset_end_at_col(2);

                                                        dx -= ( L * cx );

                                                        L.reset_offsets();
                                                        cx.reset_offsets();

                                                        cx.set_offset_start(3);
                                                        cx.set_offset_end(get_size()-nbad+1-z_end);

                                                        L.set_offset_start_at_row(3);
                                                        L.set_offset_start_at_col(3);
                                                        L.set_offset_end_at_row(get_size()-nbad+1-z_end);
                                                        L.set_offset_end_at_col(get_size()-nbad+1-z_end);

                                                        L(ARG_INVERT_LEFT,cx,dx);

                                                        L.reset_offsets();
                                                        cx.reset_offsets();
                                                    }

                                                    if ( z_end > 0 )
                                                    {
                                                        fVECTOR dx('x',z_end);

                                                        cx.set_offset_start(1);
                                                        cx.set_offset_end(get_size()-nbad+1-z_end);

                                                        L.set_offset_start_at_row(get_size()-nbad+1-z_end+1);
                                                        L.set_offset_start_at_col(1);
                                                        L.set_offset_end_at_row(get_size()-nbad+1);
                                                        L.set_offset_end_at_col(get_size()-nbad+1-z_end);

                                                        dx = -( L * cx );

                                                        L.reset_offsets();
                                                        cx.reset_offsets();

                                                        cx.set_offset_start(get_size()-nbad+1-z_end+1);
                                                        cx.set_offset_end(get_size()-nbad+1);

                                                        L.set_offset_start_at_row(get_size()-nbad+1-z_end+1);
                                                        L.set_offset_start_at_col(get_size()-nbad+1-z_end+1);
                                                        L.set_offset_end_at_row(get_size()-nbad+1);
                                                        L.set_offset_end_at_col(get_size()-nbad+1);

                                                        L(ARG_INVERT_LEFT,cx,dx);

                                                        L.reset_offsets();
                                                        cx.reset_offsets();
                                                    }
                                                }

                                                cx[1] = -cx[1];

                                                {
                                                    cx.set_offset_start(1);
                                                    cx.set_offset_end(get_size()-nbad+1);

                                                    a.set_offset_start(1);
                                                    a.set_offset_end(get_size()-nbad+1);

                                                    L.set_offset_start_at_row(1);
                                                    L.set_offset_start_at_col(1);
                                                    L.set_offset_end_at_row(get_size()-nbad+1);
                                                    L.set_offset_end_at_col(get_size()-nbad+1);

                                                    L(ARG_INVERT_RIGHT,a,cx);

                                                    L.reset_offsets();
                                                    a.reset_offsets();
                                                    cx.reset_offsets();
                                                }
                                            }
                                        }

                                        else
                                        {
                                            if ( !(elm_one_zero) )
                                            {
                                                fVECTOR cx('x',get_size()-nbad+1);

                                                {
                                                    cx[0] = 0.0;
                                                    cx[1] = b[1] / L(1,1);

                                                    if ( z_start > 1 )
                                                    {
                                                        fVECTOR dx('x',z_start-1);

                                                        cx.set_offset_start(2);
                                                        cx.set_offset_end(2);

                                                        L.set_offset_start_at_row(3);
                                                        L.set_offset_start_at_col(2);
                                                        L.set_offset_end_at_row(z_start+1);
                                                        L.set_offset_end_at_col(2);

                                                        dx = -( L * cx );

                                                        L.reset_offsets();
                                                        cx.reset_offsets();

                                                        cx.set_offset_start(3);
                                                        cx.set_offset_end(z_start+1);

                                                        L.set_offset_start_at_row(3);
                                                        L.set_offset_start_at_col(3);
                                                        L.set_offset_end_at_row(z_start+1);
                                                        L.set_offset_end_at_col(z_start+1);

                                                        L(ARG_INVERT_LEFT,cx,dx);

                                                        L.reset_offsets();
                                                        cx.reset_offsets();
                                                    }

                                                    {
                                                        fVECTOR dx('x',get_size()-nbad-z_start-z_end);

                                                        b.set_offset_start(1+z_start+1);
                                                        b.set_offset_end(get_size()-nbad+1-z_end);

                                                        dx = b;

                                                        b.reset_offsets();

                                                        cx.set_offset_start(1+1);
                                                        cx.set_offset_end(1+z_start);

                                                        L.set_offset_start_at_row(1+z_start+1);
                                                        L.set_offset_start_at_col(1+1);
                                                        L.set_offset_end_at_row(get_size()-nbad+1-z_end);
                                                        L.set_offset_end_at_col(1+z_start);

                                                        dx -= (L * cx);

                                                        L.reset_offsets();
                                                        cx.reset_offsets();

                                                        cx.set_offset_start(1+z_start+1);
                                                        cx.set_offset_end(get_size()-nbad+1-z_end);

                                                        L.set_offset_start_at_row(1+z_start+1);
                                                        L.set_offset_start_at_col(1+z_start+1);
                                                        L.set_offset_end_at_row(get_size()-nbad+1-z_end);
                                                        L.set_offset_end_at_col(get_size()-nbad+1-z_end);

                                                        L(ARG_INVERT_LEFT,cx,dx);

                                                        L.reset_offsets();
                                                        cx.reset_offsets();
                                                    }

                                                    if ( z_end > 0 )
                                                    {
                                                        fVECTOR dx('x',z_end);

                                                        cx.set_offset_start(1+1);
                                                        cx.set_offset_end(get_size()-nbad+1-z_end);

                                                        L.set_offset_start_at_row(get_size()-nbad+1-z_end+1);
                                                        L.set_offset_start_at_col(1+1);
                                                        L.set_offset_end_at_row(get_size()-nbad+1);
                                                        L.set_offset_end_at_col(get_size()-nbad+1-z_end);

                                                        dx = -( L * cx );

                                                        L.reset_offsets();
                                                        cx.reset_offsets();

                                                        cx.set_offset_start(get_size()-nbad+1-z_end+1);
                                                        cx.set_offset_end(get_size()-nbad+1);

                                                        L.set_offset_start_at_row(get_size()-nbad+1-z_end+1);
                                                        L.set_offset_start_at_col(get_size()-nbad+1-z_end+1);
                                                        L.set_offset_end_at_row(get_size()-nbad+1);
                                                        L.set_offset_end_at_col(get_size()-nbad+1);

                                                        L(ARG_INVERT_LEFT,cx,dx);

                                                        L.reset_offsets();
                                                        cx.reset_offsets();
                                                    }
                                                }

                                                cx[1] = -cx[1];

                                                {
                                                    fVECTOR dx('x',1);

                                                    cx.set_offset_start(1+1);
                                                    cx.set_offset_end(get_size()-nbad+1);

                                                    a.set_offset_start(1+1);
                                                    a.set_offset_end(get_size()-nbad+1);

                                                    L.set_offset_start_at_row(1+1);
                                                    L.set_offset_start_at_col(1+1);
                                                    L.set_offset_end_at_row(get_size()-nbad+1);
                                                    L.set_offset_end_at_col(get_size()-nbad+1);

                                                    L(ARG_INVERT_RIGHT,a,cx);

                                                    L.reset_offsets();
                                                    a.reset_offsets();
                                                    cx.reset_offsets();

                                                    a.set_offset_start(2);
                                                    a.set_offset_end(get_size()-nbad+1);

                                                    L.set_offset_start_at_row(2);
                                                    L.set_offset_start_at_col(1);
                                                    L.set_offset_end_at_row(get_size()-nbad+1);
                                                    L.set_offset_end_at_col(1);

                                                    dx = -( a * L );

                                                    L.reset_offsets();
                                                    a.reset_offsets();

                                                    a.set_offset_start(1);
                                                    a.set_offset_end(1);

                                                    L.set_offset_start_at_row(1);
                                                    L.set_offset_start_at_col(1);
                                                    L.set_offset_end_at_row(1);
                                                    L.set_offset_end_at_col(1);

                                                    L(ARG_INVERT_RIGHT,a,dx);

                                                    L.reset_offsets();
                                                    a.reset_offsets();
                                                }
                                            }

                                            else
                                            {
                                                fVECTOR cx('x',get_size()-nbad+1);

                                                {
                                                    {
                                                        cx.set_offset_start(1);
                                                        cx.set_offset_end(z_start+1);

                                                        cx = 0.0;

                                                        cx.reset_offsets();
                                                    }

                                                    {
                                                        b.set_offset_start(z_start+2);
                                                        b.set_offset_end(get_size()-nbad+1-z_end);

                                                        cx.set_offset_start(z_start+2);
                                                        cx.set_offset_end(get_size()-nbad+1-z_end);

                                                        L.set_offset_start_at_row(z_start+2);
                                                        L.set_offset_start_at_col(z_start+2);
                                                        L.set_offset_end_at_row(get_size()-nbad+1-z_end);
                                                        L.set_offset_end_at_col(get_size()-nbad+1-z_end);

                                                        L(ARG_INVERT_LEFT,cx,b);

                                                        L.reset_offsets();
                                                        cx.reset_offsets();
                                                        b.reset_offsets();
                                                    }

                                                    if ( z_end > 0 )
                                                    {
                                                        fVECTOR dx('x',z_end);

                                                        cx.set_offset_start(z_start+2);
                                                        cx.set_offset_end(get_size()-nbad+1-z_end);

                                                        L.set_offset_start_at_row(get_size()-nbad+1-z_end+1);
                                                        L.set_offset_start_at_col(z_start+2);
                                                        L.set_offset_end_at_row(get_size()-nbad+1);
                                                        L.set_offset_end_at_col(get_size()-nbad+1-z_end);

                                                        dx = -( L * cx );

                                                        L.reset_offsets();
                                                        cx.reset_offsets();

                                                        cx.set_offset_start(get_size()-nbad+1-z_end+1);
                                                        cx.set_offset_end(get_size()-nbad+1);

                                                        L.set_offset_start_at_row(get_size()-nbad+1-z_end+1);
                                                        L.set_offset_start_at_col(get_size()-nbad+1-z_end+1);
                                                        L.set_offset_end_at_row(get_size()-nbad+1);
                                                        L.set_offset_end_at_col(get_size()-nbad+1);

                                                        L(ARG_INVERT_LEFT,cx,dx);

                                                        L.reset_offsets();
                                                        cx.reset_offsets();
                                                    }
                                                }

                                                {
                                                    fVECTOR dx('x',z_start+1);

                                                    cx.set_offset_start(z_start+2);
                                                    cx.set_offset_end(get_size()-nbad+1);

                                                    a.set_offset_start(z_start+2);
                                                    a.set_offset_end(get_size()-nbad+1);

                                                    L.set_offset_start_at_row(z_start+2);
                                                    L.set_offset_start_at_col(z_start+2);
                                                    L.set_offset_end_at_row(get_size()-nbad+1);
                                                    L.set_offset_end_at_col(get_size()-nbad+1);

                                                    L(ARG_INVERT_RIGHT,a,cx);

                                                    L.reset_offsets();
                                                    a.reset_offsets();
                                                    cx.reset_offsets();

                                                    a.set_offset_start(z_start+2);
                                                    a.set_offset_end(get_size()-nbad+1);

                                                    L.set_offset_start_at_row(z_start+2);
                                                    L.set_offset_start_at_col(1);
                                                    L.set_offset_end_at_row(get_size()-nbad+1);
                                                    L.set_offset_end_at_col(z_start+1);

                                                    dx = -( a * L );

                                                    L.reset_offsets();
                                                    a.reset_offsets();

                                                    a.set_offset_start(1);
                                                    a.set_offset_end(z_start+1);

                                                    L.set_offset_start_at_row(1);
                                                    L.set_offset_start_at_col(1);
                                                    L.set_offset_end_at_row(z_start+1);
                                                    L.set_offset_end_at_col(z_start+1);

                                                    L(ARG_INVERT_RIGHT,a,dx);

                                                    L.reset_offsets();
                                                    a.reset_offsets();
                                                }
                                            }
                                        }
                                    }
                                }

                                break;
                            }

                            case FACT_NONE:
                            default:
                            {
                                L_THROW(5);

                                break;
                            }
                        }
                    }

                    else
                    {
                        a = 0.0;
                    }
                }

                exit_fixvect_normal_fn(a);

                ay = a[0];

                a.remove(1);

                ax = a;

                break;
            }

            default:
            {
                L_THROW(6);

                break;
            }
        }

        result = nbad;
    }

    FN_EXIT_POINT result;
}

long fMATRFact::near_invert(const lVECTOR &Gindex, fVECTOR &a, const fMATRIX &G_real)
{
    FN_ENTRY_POINT

    THROW_ASSERT(a.get_effective_size() == (get_size()-nbad));
    THROW_ASSERT(nbad > 0);

    long result;

    if ( f_type == FACT_NONE )
    {
        fMATRFact temp(use_index,Gindex,G_real,zt,DEFAULT_FACT);

        temp.near_invert(Gindex,a,G_real);

        result = temp.get_nbad();
    }

    else
    {
        if ( get_size()-nbad > 0 )
        {
            switch ( f_form )
            {
                case FORM_ONE:
                {
                    if ( !s_valid )
                    {
                        switch ( f_type )
                        {
                            case FACT_INVERSE:
                            {
                                fVECTOR b;

                                L.set_offset_start_at_row(1);
                                L.set_offset_start_at_col(get_size()-nbad+1);
                                L.set_offset_end_at_row(get_size()-nbad);
                                L.set_offset_end_at_col(get_size()-nbad+1);

                                b = L(0);

                                L.reset_offsets();

                                L.set_offset_start_at_row(1);
                                L.set_offset_start_at_col(1);
                                L.set_offset_end_at_row(get_size()-nbad);
                                L.set_offset_end_at_col(get_size()-nbad);

                                s = ( L * b );

                                L.reset_offsets();

                                s_valid = 1;

                                break;
                            }

                            case FACT_CHOL:
                            {
                                fVECTOR b;

                                L.set_offset_start_at_row(get_size()-nbad+1);
                                L.set_offset_start_at_col(1);
                                L.set_offset_end_at_row(get_size()-nbad+1);
                                L.set_offset_end_at_col(get_size()-nbad);

                                b = L[0];

                                L.reset_offsets();

                                L.set_offset_start_at_row(1);
                                L.set_offset_start_at_col(1);
                                L.set_offset_end_at_row(get_size()-nbad);
                                L.set_offset_end_at_col(get_size()-nbad);

                                s.pad_vector(b.get_effective_size());

                                L(ARG_INVERT_RIGHT,s,b);

                                L.reset_offsets();

                                s_valid = 1;

                                break;
                            }

                            case FACT_NONE:
                            default:
                            {
                                L_THROW(0);

                                break;
                            }
                        }
                    }

                    a = s;

                    exit_fixvect_normal_fn(a);

                    break;
                }

                case FORM_TWO:
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

        result = nbad;
    }

    FN_EXIT_POINT result;
}

long fMATRFact::near_invert(const lVECTOR &Gindex, fVECTOR &ax, L_DOUBLE &ay, const fMATRIX &G_real, const fVECTOR &d_real)
{
    FN_ENTRY_POINT

    long result;

    if ( f_type == FACT_NONE )
    {
        fMATRFact temp(use_index,Gindex,G_real,d_real,zt,DEFAULT_FACT);

        temp.near_invert(Gindex,ax,ay,G_real,d_real);

        result = temp.get_nbad();
    }

    else
    {
        switch ( f_form )
        {
            case FORM_ONE:
            {
                ay = 0.0;

                near_invert(Gindex,ax,G_real);

                break;
            }

            case FORM_TWO:
            {
                THROW_ASSERT(ax.get_effective_size() == (get_size()-nbad));
                THROW_ASSERT(nbad > 0);

                ay = 0.0;

                if ( get_size()-nbad > 0 )
                {
                    fVECTOR a;

                    a = ax;

                    a.addstart(ay);

                    if ( !s_valid )
                    {
                        switch ( f_type )
                        {
                            case FACT_INVERSE:
                            {
                                fVECTOR b;

                                L.set_offset_start_at_row(1);
                                L.set_offset_start_at_col(get_size()-nbad+2);
                                L.set_offset_end_at_row(get_size()-nbad+1);
                                L.set_offset_end_at_col(get_size()-nbad+2);

                                b = L(0);

                                L.reset_offsets();

                                L.set_offset_start_at_row(1);
                                L.set_offset_start_at_col(1);
                                L.set_offset_end_at_row(get_size()-nbad+1);
                                L.set_offset_end_at_col(get_size()-nbad+1);

                                s = ( L * b );

                                L.reset_offsets();

                                s_valid = 1;

                                a = s;

                                break;
                            }

                            case FACT_CHOL:
                            {
                                if ( case_four_flag )
                                {
                                    fVECTOR b;
                                    L_DOUBLE d_b_conj;

                                    L.set_offset_start_at_row(3);
                                    L.set_offset_start_at_col(1);
                                    L.set_offset_end_at_row(3);
                                    L.set_offset_end_at_col(2);

                                    b = L[0];

                                    L.reset_offsets();

                                    d_b_conj = L(1,0);

                                    s[0] = ( L(1,0)   * b[1] );
                                    s[1] = ( d_b_conj * b[0] );
                                }

                                else
                                {
                                    fVECTOR b;

                                    L.set_offset_start_at_row(get_size()-nbad+2);
                                    L.set_offset_start_at_col(1);
                                    L.set_offset_end_at_row(get_size()-nbad+2);
                                    L.set_offset_end_at_col(get_size()-nbad+1);

                                    b = L[0];

                                    L.reset_offsets();

                                    L.set_offset_start_at_row(1);
                                    L.set_offset_start_at_col(1);
                                    L.set_offset_end_at_row(get_size()-nbad+1);
                                    L.set_offset_end_at_col(get_size()-nbad+1);

                                    s.pad_vector(b.get_effective_size());

                                    L(ARG_INVERT_RIGHT,s,b);

                                    L.reset_offsets();
                                }

                                s_valid = 1;

                                a = s;

                                break;
                            }

                            case FACT_NONE:
                            default:
                            {
                                L_THROW(0);

                                break;
                            }
                        }
                    }

                    exit_fixvect_normal_fn(a);

                    ay = a[0];

                    a.remove(1);

                    ax = a;
                }

                break;
            }

            default:
            {
                L_THROW(3);

                break;
            }
        }

        result = nbad;
    }

    FN_EXIT_POINT result;
}




//
// Matrix manipulations:
//

long fMATRFact::fact_addend(const lVECTOR &Gindex, const fVECTOR &ax, const fMATRIX &G_real)
{
    FN_ENTRY_POINT

    long result;

    if ( f_type == FACT_NONE )
    {
        result = -1;

        fact_none_size++;
    }

    else
    {
        switch ( f_form )
        {
            case FORM_ONE:
            {
                fVECTOR a('x',get_size()+1);

                THROW_ASSERT((ax.get_effective_size() < 0) || (ax.get_effective_size() == (get_size()+1)));

                a = ax;

                // We need to put L to the correct size temporarily so that
                // get_size() will FN_EXIT_POINT the correct value when called.

                L.addend();
                enter_fixvect_normal_fn(a);
                L.remove(get_size());

                switch ( f_type )
                {
                    case FACT_INVERSE:
                    {
                        L.addend(a[get_size()],a);

                        nbad++;

                        break;
                    }

                    case FACT_CHOL:
                    {
                        if ( nbad < get_size() )
                        {
                            fVECTOR cx('x',get_size()-nbad);

                            a.set_offset_start(1);
                            a.set_offset_end(get_size()-nbad);

                            L.set_offset_start_at_row(1);
                            L.set_offset_start_at_col(1);
                            L.set_offset_end_at_row(get_size()-nbad);
                            L.set_offset_end_at_col(get_size()-nbad);

                            L(ARG_INVERT_LEFT,cx,a);

                            a = cx;

                            L.reset_offsets();
                            a.reset_offsets();
                        }

                        L.addend(a[get_size()],a);

                        nbad++;

                        break;
                    }

                    case FACT_NONE:
                    default:
                    {
                        L_THROW(2);

                        break;
                    }
                }

                break;
            }

            case FORM_TWO:
            {
                L_THROW(3);

                break;
            }

            default:
            {
                L_THROW(4);

                break;
            }
        }

        fix_symmetry();

        result = xfact(Gindex,G_real);
    }

    FN_EXIT_POINT result;
}

long fMATRFact::fact_addend(const lVECTOR &Gindex, const fVECTOR &ax, const L_DOUBLE &d, const fMATRIX &G_real, const fVECTOR &d_real)
{
    FN_ENTRY_POINT

    long result;

    if ( f_type == FACT_NONE )
    {
        result = -1;

        fact_none_size++;
    }

    else
    {
        switch ( f_form )
        {
            case FORM_ONE:
            {
                fact_addend(Gindex,ax,G_real);

                break;
            }

            case FORM_TWO:
            {
                THROW_ASSERT((ax.get_effective_size() < 0) || (ax.get_effective_size() == (get_size()+1)));

                fVECTOR a('x',get_size()+1);

                a = ax;

                a.addstart(d);

                // We need to put L to the correct size temporarily so that
                // get_size() will return the correct value when called.

                L.addend();
                enter_fixvect_normal_fn(a);
                L.remove(get_size()+1);

                switch ( f_type )
                {
                    case FACT_INVERSE:
                    {
                        switch ( get_size() )
                        {
                            case 0:
                            {
                                // Both G and L have size get_size()+1 (L had d
                                // at start, which will be re-added shortly.  G
                                // has an extra row not in factorisation).

                                overwrite_L(Gindex,G_real); // modded

                                refresh_fact(d_real);

                                break;
                            }

                            default:
                            {
                                L.addend(a[get_size()+1],a);

                                nbad++;

                                break;
                            }
                        }

                        break;
                    }

                    case FACT_CHOL:
                    {
                        switch ( get_size() )
                        {
                            case 0:
                            case 1:
                            {
                                // Both G and L have size get_size()+1 (L had d
                                // at start, which will be re-added shortly.  G
                                // has an extra row not in factorisation).

                                overwrite_L(Gindex,G_real); // modded

                                refresh_fact(d_real);

                                break;
                            }

                            default:
                            {
                                if ( !case_four_flag )
                                {
                                    if ( nbad < get_size() )
                                    {
                                        fVECTOR cx('x',get_size()-nbad+1);

                                        a.set_offset_start(1);
                                        a.set_offset_end(get_size()-nbad+1);

                                        L.set_offset_start_at_row(1);
                                        L.set_offset_start_at_col(1);
                                        L.set_offset_end_at_row(get_size()-nbad+1);
                                        L.set_offset_end_at_col(get_size()-nbad+1);

                                        L(ARG_INVERT_LEFT,cx,a);

                                        cx[1] = -cx[1];

                                        a = cx;

                                        L.reset_offsets();
                                        a.reset_offsets();
                                    }
                                }

                                L.addend(a[get_size()+1],a);

                                nbad++;

                                break;
                            }
                        }

                        break;
                    }

                    case FACT_NONE:
                    default:
                    {
                        L_THROW(4);

                        break;
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

        fix_symmetry();

        result = xfact(Gindex,G_real,d_real);
    }

    FN_EXIT_POINT result;
}

long fMATRFact::fact_shrink(const lVECTOR &Gindex, long ix, const fMATRIX &G_real)
{
    FN_ENTRY_POINT

    long result;

    if ( f_type == FACT_NONE )
    {
        result = -1;

        fact_none_size--;
    }

    else
    {
        long i,j,k,l;
        long refresh_size = 0;
        int we_failed = 0;

        switch ( f_form )
        {
            case FORM_ONE:
            {
                i = ix;

                THROW_ASSERT(i >= 1);
                THROW_ASSERT(i <= get_size());

                switch ( f_type )
                {
                    case FACT_INVERSE:
                    {
                        // Let's fix up s first off

                        if ( i <= get_size()-nbad+1 )
                        {
                            s.make_zero(DO_FORCE);
                            s_valid = 0;
                        }

                        fVECTOR bx;
                        L_DOUBLE ax;

                        switch ( get_size() )
                        {
                            case 1:
                            {
                                // H no longer exists

                                L.remove(1);

                                nbad = 0;

                                break;
                            }

                            default:
                            {
                                // We have two possibilities.  Either the row/
                                // column is in the inverted part of the matrix
                                // or it is not.

                                if ( i <= get_size()-nbad )
                                {
                                    // case 1. The row/column is in the inverted
                                    // part of the matrix.

                                    // a. if there is only the inverted part is
                                    // 1*1 then things are really simple.

                                    if ( get_size()-nbad == 1 )
                                    {
                                        L.remove(1);
                                    }

                                    else
                                    {
                                        // b. there are more elements in the
                                        // inverted part that need to be updated.

                                        // Now, assuming the diagonal being
                                        // removed is non-zero, the update will
                                        // be simple.  Attempt this - if an
                                        // inversion fails then we must try
                                        // something a little more complex.

                                        // It is convenient to position it
                                        // at the end of the inverted part
                                        // of the matrix n (temporarily).

                                        L.fswap(i,get_size()-nbad);

                                        try
                                        {
                                            L.set_offset_start_at_row(1);
                                            L.set_offset_start_at_col(1);
                                            L.set_offset_end_at_row(get_size()-nbad);
                                            L.set_offset_end_at_col(get_size()-nbad);

                                            bx = L((get_size()-nbad)-1);

                                            L.reset_offsets();

                                            NZ_ASSERT(bx[(get_size()-nbad)-1],zt);
                                            ax = -1.0/bx[(get_size()-nbad)-1];

                                            // OK - if the operation was going to
                                            // fail, it would have failed on the
                                            // previous line, so we can safely
                                            // finish the operation now.

                                            // NOTE - we need not worry about
                                            // what will happen to the row of
                                            // L that we are removing.

                                            L.set_offset_start_at_row(1);
                                            L.set_offset_start_at_col(1);
                                            L.set_offset_end_at_row(get_size()-nbad);
                                            L.set_offset_end_at_col(get_size()-nbad);

                                            L(ARG_RANK_ONE,bx,ax);

                                            L.reset_offsets();

                                            refresh_size = 0;
                                        }

                                        catch ( int err_num )
                                        {
                                            switch ( err_num )
                                            {
                                                case INV_ERR:
                                                {
                                                    // If we have reached this
                                                    // point then we know that
                                                    // bx cannot be inverted, so
                                                    // we cannot proceed in the
                                                    // usual fashion.

                                                    // The following flag is
                                                    // used to indicate if the
                                                    // operation has completely
                                                    // failed and it is
                                                    // necessary to re-start from
                                                    // scratch.

                                                    we_failed = 1;

                                                    if ( (get_size()-nbad)/2 >= 2 )
                                                    {
                                                        // If we have 4 or more
                                                        // elements in the
                                                        // inverted portion of
                                                        // L, then we will
                                                        // attempt to find a
                                                        // slightly larger non-
                                                        // singular block that
                                                        // may be removed.  The
                                                        // following code tries
                                                        // progressively larger
                                                        // square blocks (size
                                                        // i*i) on the bottom
                                                        // right of the inverted
                                                        // upper left of H

                                                        for ( j = 2 ; j <= (get_size()-nbad)/2 ; j++ )
                                                        {
                                                            try
                                                            {
                                                                fMATRIX tempb;

                                                                L.set_offset_start_at_row(get_size()-nbad-j+1);
                                                                L.set_offset_start_at_col(get_size()-nbad-j+1);
                                                                L.set_offset_end_at_row(get_size()-nbad);
                                                                L.set_offset_end_at_col(get_size()-nbad);

                                                                tempb = -L;

                                                                tempb.invert_matrix(zt);

                                                                L.reset_offsets();

                                                                // OK - to have
                                                                // got to this
                                                                // point, we must
                                                                // have found a
                                                                // non-singular
                                                                // square. We
                                                                // continue with
                                                                // a rank-n
                                                                // analogy of
                                                                // the method we
                                                                // failed to use
                                                                // previously.

                                                                fMATRIX tempc;

                                                                L.set_offset_start_at_col(get_size()-nbad-j+1);
                                                                L.set_offset_start_at_row(1);
                                                                L.set_offset_end_at_col(get_size()-nbad);
                                                                L.set_offset_end_at_row(get_size()-nbad-j);

                                                                tempc = L;




                                                                L.reset_offsets();

                                                                L.set_offset_start_at_row(1);
                                                                L.set_offset_start_at_col(1);
                                                                L.set_offset_end_at_row(get_size()-nbad-j);
                                                                L.set_offset_end_at_col(get_size()-nbad-j);

                                                                L(ARG_RANK_N,tempc,tempb);

                                                                L.reset_offsets();

                                                                // Register our
                                                                // success.

                                                                we_failed = 0;

                                                                refresh_size = j-1;

                                                                goto get_out_a;

                                                                break;
                                                            }

                                                            catch ( int x_num )
                                                            {
                                                                switch ( x_num )
                                                                {
                                                                    case INV_ERR:
                                                                    {
                                                                        // OK -
                                                                        // that
                                                                        // one
                                                                        // failed
                                                                        // - on
                                                                        // to the
                                                                        // next
                                                                        // one.

                                                                        L.reset_offsets();

                                                                        break;
                                                                    }

                                                                    default:
                                                                    {
                                                                        throw x_num;

                                                                        break;
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }

                                                    get_out_a:

                                                    if ( we_failed )
                                                    {
                                                        refresh_size = get_size()-nbad-1;
                                                    }

                                                    break;
                                                }

                                                default:
                                                {
                                                    throw err_num;

                                                    break;
                                                }
                                            }
                                        }

                                        // Having done what we came to do, we
                                        // can re-order things

                                        L.bswap(i,get_size()-nbad);

                                        // refresh_size is the size of the block
                                        // that has been removed from the
                                        // inverted part of the factorisation,
                                        // not including (at present) the row/
                                        // column being removed from H
                                        // altogether.  We now add this if
                                        // needed.

                                        l = get_size()-nbad;

                                        if ( i >= l+1-refresh_size )
                                        {
                                            refresh_size++;
                                        }

                                        // Overwrite the relevant part and update
                                        // nbad as we go.

                                        if ( refresh_size >= 1 )
                                        {
                                            for ( j = l+1-refresh_size ; j <= l ; j++ )
                                            {
                                                long jj;
                                                long kk;

                                                if ( j != i )
                                                {
                                                    for ( k = 1 ; k <= j ; k++ )
                                                    {
                                                        if ( k != i )
                                                        {
                                                            jj = convert_to_virtual_pos(j);
                                                            kk = convert_to_virtual_pos(k);

                                                            if ( jj > ix ) { jj--; }
                                                            if ( kk > ix ) { kk--; }

                                                            L(j-1,k-1) = get_deindexed_value(Gindex,G_real,jj,kk);
                                                        }
                                                    }
                                                }

                                                nbad++;
                                            }
                                        }

                                        // Lastly, we remove the relevant rows
                                        // and columns from the matrix and its
                                        // diagonal mirror.

                                        L.remove(i);

                                        // Do or do not, there is no try.
                                        // - yoda
                                    }
                                }

                                else
                                {
                                    // case 2. The row/column is not in the
                                    // inverted part of the matrix.

                                    L.remove(i);

                                    nbad--;
                                }

                                break;
                            }
                        }

                        break;
                    }

                    case FACT_CHOL:
                    {
                        if ( i <= get_size()-nbad+1 )
                        {
                            s.make_zero(DO_FORCE);
                            s_valid = 0;
                        }

                        switch ( get_size() )
                        {
                            case 0:
                            {
                                L_THROW(0);

                                break;
                            }

                            case 1:
                            {
                                L.remove(1);

                                nbad = 0;

                                break;
                            }

                            default:
                            {
                                // We have two possibilities.  Either the row/
                                // column is in the inverted part of the matrix
                                // or it is not.

                                if ( i <= get_size()-nbad )
                                {
                                    // case 1. The row/column is in the inverted
                                    // part of the matrix.

                                    // a. if there is only the inverted part is
                                    // 1*1 then things are really simple.

                                    if ( get_size()-nbad == 1 )
                                    {
                                        L.remove(1);
                                    }

                                    // b. there are more elements in the inverted
                                    // part that need to be updated.

                                    else
                                    {
                                        if ( i == get_size()-nbad )
                                        {
                                            L.remove(i);
                                        }

                                        else
                                        {
                                            fVECTOR qx('x',get_size());

                                            L.set_offset_start_at_row(1);
                                            L.set_offset_start_at_col(i);
                                            L.set_offset_end_at_row(get_size());
                                            L.set_offset_end_at_col(i);

                                            qx = L(0);

                                            L.reset_offsets();

                                            L.remove(i);
                                            qx.remove(i);

                                            xrankone(Gindex,qx,1.0,G_real,1,i-1,0);
                                        }
                                    }
                                }

                                else
                                {
                                    // case 2. The row/column is not in the
                                    // inverted part of the matrix.

                                    L.remove(i);

                                    nbad--;
                                }

                                break;
                            }
                        }

                        break;
                    }

                    case FACT_NONE:
                    default:
                    {
                        L_THROW(0);

                        break;
                    }
                }

                break;
            }

            case FORM_TWO:
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

        fix_symmetry();

        result = xfact(Gindex,G_real);
    }

    FN_EXIT_POINT result;
}

long fMATRFact::fact_shrink(const lVECTOR &Gindex, long ix, const fMATRIX &G_real, const fVECTOR &d_real)
{
    FN_ENTRY_POINT

    long result;

    if ( f_type == FACT_NONE )
    {
        result = -1;

        fact_none_size--;
    }

    else
    {
        long i,j,k,l;
        long refresh_size = 0;
        int we_failed = 0;
        fVECTOR d_temp;

        switch ( f_form )
        {
            case FORM_ONE:
            {
                fact_shrink(Gindex,ix,G_real);

                break;
            }

            case FORM_TWO:
            {
                i = convert_to_actual_pos(ix);

                THROW_ASSERT(i >= 1);
                THROW_ASSERT(i <= get_size());

                switch ( f_type )
                {
                    case FACT_INVERSE:
                    {
                        // Let's fix up s first off

                        if ( i <= get_size()-nbad+1 )
                        {
                            s.make_zero(DO_FORCE);
                            s_valid = 0;
                        }

                        fVECTOR bx;
                        L_DOUBLE ax;

                        switch ( get_size() )
                        {
                            case 1:
                            {
                                // H no longer exists

                                L.remove(2);

                                L(0,0) = 0.0;

                                nbad = 0;

                                break;
                            }

                            default:
                            {
                                // We have two possibilities.  Either the row/
                                // column is in the inverted part of the matrix
                                // or it is not.

                                if ( i <= get_size()-nbad )
                                {
                                    // case 1. The row/column is in the inverted
                                    // part of the matrix.

                                    // a. if there is only the inverted part is
                                    // 1*1 then things are really simple.

                                    if ( get_size()-nbad == 1 )
                                    {
                                        L.remove(2);

                                        L(0,0) = 0.0;
                                    }

                                    else
                                    {
                                        // b. there are more elements in the
                                        // inverted part that need to be updated.

                                        // Now, assuming the diagonal being
                                        // removed is non-zero, the update will
                                        // be simple.  Attempt this - if an
                                        // inversion fails then we must try
                                        // something a little more complex.

                                        // It is convenient to position it
                                        // at the end of the inverted part
                                        // of the matrix n (temporarily).

                                        L.fswap(i+1,get_size()-nbad+1);

                                        try
                                        {
                                            L.set_offset_start_at_row(1);
                                            L.set_offset_start_at_col(1);
                                            L.set_offset_end_at_row(get_size()-nbad+1);
                                            L.set_offset_end_at_col(get_size()-nbad+1);

                                            bx = L(get_size()-nbad);

                                            L.reset_offsets();

                                            NZ_ASSERT(bx[get_size()-nbad],zt);
                                            ax = -1.0/bx[get_size()-nbad];

                                            // OK - if the operation was going to
                                            // fail, it would have failed on the
                                            // previous line, so we can safely
                                            // finish the operation now.

                                            // NOTE - we need not worry about
                                            // what will happen to the row of
                                            // L that we are removing.

                                            L.set_offset_start_at_row(1);
                                            L.set_offset_start_at_col(1);
                                            L.set_offset_end_at_row(get_size()-nbad+1);
                                            L.set_offset_end_at_col(get_size()-nbad+1);

                                            L(ARG_RANK_ONE,bx,ax);

                                            L.reset_offsets();

                                            refresh_size = 0;
                                        }

                                        catch ( int err_num )
                                        {
                                            switch ( err_num )
                                            {
                                                case INV_ERR:
                                                {
                                                    // If we have reached this
                                                    // point then we know that
                                                    // bx cannot be inverted, so
                                                    // we cannot proceed in the
                                                    // usual fashion.

                                                    // The following flag is
                                                    // used to indicate if the
                                                    // operation has completely
                                                    // failed and it is
                                                    // necessary to re-start from
                                                    // scratch.

                                                    we_failed = 1;

                                                    if ( (get_size()-nbad)/2 >= 2 )
                                                    {
                                                        // If we have 3 or more
                                                        // elements in the
                                                        // inverted portion of
                                                        // L, then we will
                                                        // attempt to find a
                                                        // slightly larger non-
                                                        // singular block that
                                                        // may be removed.  The
                                                        // following code tries
                                                        // progressively larger
                                                        // square blocks (size
                                                        // i*i) on the bottom
                                                        // right of the inverted
                                                        // upper left of H
                                                        // (except for a
                                                        // get_size()-nbad
                                                        // square, in which case
                                                        // we may as well just
                                                        // start from scratch).

                                                        for ( j = 2 ; j <= (get_size()-nbad)/2 ; j++ )
                                                        {
                                                            try
                                                            {
                                                                fMATRIX tempb;

                                                                L.set_offset_start_at_row(get_size()-nbad-j+2);
                                                                L.set_offset_start_at_col(get_size()-nbad-j+2);
                                                                L.set_offset_end_at_row(get_size()-nbad+1);
                                                                L.set_offset_end_at_col(get_size()-nbad+1);

                                                                tempb = -L;

                                                                tempb.invert_matrix(zt);

                                                                L.reset_offsets();

                                                                // OK - to have
                                                                // got to this
                                                                // point, we must
                                                                // have found a
                                                                // non-singular
                                                                // square. We
                                                                // continue with
                                                                // a rank-n
                                                                // analogy of
                                                                // the method we
                                                                // failed to use
                                                                // previously.

                                                                fMATRIX tempc;

                                                                L.set_offset_start_at_col(get_size()-nbad-j+2);
                                                                L.set_offset_start_at_row(1);
                                                                L.set_offset_end_at_col(get_size()-nbad+1);
                                                                L.set_offset_end_at_row(get_size()-nbad-j+1);

                                                                tempc = L;




                                                                L.reset_offsets();

                                                                L.set_offset_start_at_row(1);
                                                                L.set_offset_start_at_col(1);
                                                                L.set_offset_end_at_row(get_size()-nbad-j+1);
                                                                L.set_offset_end_at_col(get_size()-nbad-j+1);

                                                                L(ARG_RANK_N,tempc,tempb);

                                                                L.reset_offsets();

                                                                // Register our
                                                                // success.

                                                                we_failed = 0;

                                                                refresh_size = j-1;

                                                                goto get_out_a;

                                                                break;
                                                            }

                                                            catch ( int x_num )
                                                            {
                                                                switch ( x_num )
                                                                {
                                                                    case INV_ERR:
                                                                    {
                                                                        // OK -
                                                                        // that
                                                                        // one
                                                                        // failed
                                                                        // - on
                                                                        // to the
                                                                        // next
                                                                        // one.

                                                                        L.reset_offsets();

                                                                        break;
                                                                    }

                                                                    default:
                                                                    {
                                                                        throw x_num;

                                                                        break;
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }

                                                    get_out_a:

                                                    if ( we_failed )
                                                    {
                                                        refresh_size = get_size()-nbad-1;
                                                    }

                                                    break;
                                                }

                                                default:
                                                {
                                                    throw err_num;

                                                    break;
                                                }
                                            }
                                        }

                                        // Having done what we came to do, we
                                        // can re-order things

                                        L.bswap(i+1,get_size()-nbad+1);

                                        // refresh_size is the size of the block
                                        // that has been removed from the
                                        // inverted part of the factorisation,
                                        // not including (at present) the row/
                                        // column being removed from H
                                        // altogether.  We now add this if
                                        // needed.

                                        l = get_size()-nbad;

                                        if ( i >= l+1-refresh_size )
                                        {
                                            refresh_size++;
                                        }

                                        // Overwrite the relevant part and update
                                        // nbad as we go.

                                        if ( refresh_size >= 1 )
                                        {
                                            long jj;
                                            long kk;

                                            if ( we_failed )
                                            {
                                                L(0,0) = 0.0;
                                            }

                                            for ( j = l+1-refresh_size+1 ; j <= l+1 ; j++ )
                                            {
                                                if ( j != i )
                                                {
                                                    for ( k = 1+1 ; k <= j ; k++ )
                                                    {
                                                        if ( k != i )
                                                        {
                                                            jj = convert_to_virtual_pos(j-1);
                                                            kk = convert_to_virtual_pos(k-1);

                                                            if ( jj > ix ) { jj--; }
                                                            if ( kk > ix ) { kk--; }

                                                            L(j-1,k-1) = get_deindexed_value(Gindex,G_real,jj,kk);

                                                            // The bizarre dec
                                                            // on j and k in
                                                            // the call here
                                                            // is intentional as
                                                            // j and k refer to
                                                            // a H position here,
                                                            // which is not what
                                                            // we need when
                                                            // reffing G_real.
                                                        }
                                                    }
                                                }

                                                nbad++;
                                            }
                                        }

                                        // Lastly, we remove the relevant rows
                                        // and columns from the matrix and its
                                        // diagonal mirror.

                                        L.remove(i+1);

                                        // --- Insert wise words here. ---
                                    }
                                }

                                else
                                {
                                    // case 2. The row/column is not in the
                                    // inverted part of the matrix.

                                    L.remove(i+1);

                                    nbad--;
                                }

                                break;
                            }
                        }

                        break;
                    }

                    case FACT_CHOL:
                    {
                        if ( i <= get_size()-nbad+1 )
                        {
                            s_valid = 0;

                            s.make_zero(DO_FORCE);
                        }

                        switch ( get_size() )
                        {
                            case 1:
                            {
                                L.remove(2);

                                L(0,0) = 0.0;

                                nbad = 0;

                                break;
                            }

                            case 2:
                            {
                                long d_size = get_size()-1;

                                // L is being shrunk - currently it has size
                                // dsize+2 (with d at start), after op will have
                                // size dsize+1 (with d at end).  We need to
                                // overwrite L with current content (G of size
                                // dsize+1, d at start) before shrinking by 1.

                                L.remove(d_size+2);
                                overwrite_L(Gindex,G_real); // modded
 
                                fVECTOR d_loc('x',d_size+1);

                                d_loc = d_real;

                                L.remove(ix);
                                d_loc.remove(ix);

                                refresh_fact(d_loc);

                                break;
                            }

                            default:
                            {
                                // We have two possibilities.  Either the row/
                                // column is in the inverted part of the matrix
                                // or it is not.

                                if ( i <= get_size()-nbad )
                                {
                                    // case 1. The row/column is in the inverted
                                    // part of the matrix.

                                    switch ( i )
                                    {
                                        case 1:
                                        {
                                            // NB - there is probably a more
                                            // optimal and neat way of doing
                                            // this bit, but this is close
                                            // enough for now.

                                            long d_size = get_size()-1;


                                            L.remove(d_size+2);
                                            overwrite_L(Gindex,G_real); // modded

                                            fVECTOR d_loc('x',d_size+1);

                                            d_loc = d_real;

                                            L.remove(ix);
                                            d_loc.remove(ix);

                                            refresh_fact(d_loc);

                                            break;
                                        }

                                        default:
                                        {
                                            if ( i == get_size()-nbad )
                                            {
                                                L.remove(i+1);
                                            }

                                            else
                                            {
                                                fVECTOR qx('x',get_size()+1);

                                                L.set_offset_start_at_row(1);
                                                L.set_offset_start_at_col(i+1);
                                                L.set_offset_end_at_row(get_size()+1);
                                                L.set_offset_end_at_col(i+1);

                                                qx = L(0);

                                                L.reset_offsets();

                                                L.remove(i+1);
                                                qx.remove(i);

                                                xrankone(Gindex,qx,1.0,G_real,d_real,1,i-1,0);
                                            }

                                            break;
                                        }
                                    }
                                }

                                else
                                {
                                    // case 2. The row/column is not in the
                                    // inverted part of the matrix.

                                    // Note that we CANNOT have i == 1

                                    L.remove(i+1);

                                    nbad--;
                                }

                                break;
                            }
                        }

                        break;
                    }

                    case FACT_NONE:
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

        fix_symmetry();

        result = xfact(Gindex,G_real,d_real);
    }

    FN_EXIT_POINT result;
}

//
// Information functions
//
// get_size() - get size - size of G
// get_nbad() - get nbad - 0 if non-singular, >= 1 if singular
//
// get_case_four_flag()    - get case_four_flag
// get_primary_swap_flag() - get primary_swap_flag
//
// get_f_form() - get f_sign
// get_f_type() - get f_type
//

long fMATRFact::get_size(void)
{
    FN_ENTRY_POINT

    long result;

    if ( f_type == FACT_NONE )
    {
        result = fact_none_size;
    }

    else
    {
        switch ( f_form )
        {
            case FORM_ONE:
            {
                result = L.get_real_size();

                break;
            }

            case FORM_TWO:
            {
                result = L.get_real_size() - 1;

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

long fMATRFact::get_nbad(void)
{
    FN_ENTRY_POINT

    if ( f_type == FACT_NONE )
    {
        return -1;
    }

    FN_EXIT_POINT nbad;
}

int fMATRFact::get_case_four_flag(void)
{
    FN_ENTRY_POINT

    int result;

    switch ( f_form )
    {
        case FORM_ONE:
        {
            result = 0;

            break;
        }

        case FORM_TWO:
        {
            result = case_four_flag;

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

int fMATRFact::get_primary_swap_flag(void)
{
    FN_ENTRY_POINT

    int result;

    switch ( f_form )
    {
        case FORM_ONE:
        {
            result = 0;

            break;
        }

        case FORM_TWO:
        {
            result = primary_swap_flag;

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

int fMATRFact::get_f_type(void)
{
    FN_ENTRY_POINT

    FN_EXIT_POINT f_type;
}

int fMATRFact::get_f_form(void)
{
    FN_ENTRY_POINT

    FN_EXIT_POINT f_form;
}

fMATRIX fMATRFact::test_fact(void)
{
    FN_ENTRY_POINT

    if ( f_type == FACT_NONE )
    {
        L_THROW(0);
    }

    fMATRIX result;
    long i;

    switch ( f_form )
    {
        case FORM_ONE:
        {
            result.make_symmetric(DO_FORCE);
            result.pad_matrix(get_size());

            if ( get_size() > 0 )
            {
                switch ( f_type )
                {
                    case FACT_INVERSE:
                    {
                        result = L;

                        if ( get_size()-nbad > 0 )
                        {
                            result.set_offset_start_at_row(1);
                            result.set_offset_start_at_col(1);
                            result.set_offset_end_at_row(get_size()-nbad);
                            result.set_offset_end_at_col(get_size()-nbad);

                            L.set_offset_start_at_row(1);
                            L.set_offset_start_at_col(1);
                            L.set_offset_end_at_row(get_size()-nbad);
                            L.set_offset_end_at_col(get_size()-nbad);

                            result.invert_matrix();

                            L.reset_offsets();

                            result.reset_offsets();
                        }

                        break;
                    }

                    case FACT_CHOL:
                    {
                        result = L;
                        result.make_symmetric(DO_FORCE);

                        if ( get_size()-nbad > 0 )
                        {
                            fMATRIX M;

                            L.set_offset_start_at_row(1);
                            L.set_offset_start_at_col(1);
                            L.set_offset_end_at_row(get_size()-nbad);
                            L.set_offset_end_at_col(get_size()-nbad);

                            M = L;

                            M.transpose();

                            L.reset_offsets();

                            result.set_offset_start_at_row(1);
                            result.set_offset_start_at_col(1);
                            result.set_offset_end_at_row(get_size()-nbad);
                            result.set_offset_end_at_col(get_size()-nbad);

                            L.set_offset_start_at_row(1);
                            L.set_offset_start_at_col(1);
                            L.set_offset_end_at_row(get_size()-nbad);
                            L.set_offset_end_at_col(get_size()-nbad);

                            result =  L;
                            result *= M;

                            result.make_symmetric(DO_FORCE);

                            L.reset_offsets();
                            result.reset_offsets();

                            if ( nbad > 0 )
                            {
                                result.set_offset_start_at_row(get_size()-nbad+1);
                                result.set_offset_start_at_col(1);
                                result.set_offset_end_at_row(get_size());
                                result.set_offset_end_at_col(get_size()-nbad);

                                L.set_offset_start_at_row(get_size()-nbad+1);
                                L.set_offset_start_at_col(1);
                                L.set_offset_end_at_row(get_size());
                                L.set_offset_end_at_col(get_size()-nbad);

                                result =  L;
                                result *= M;

                                result.make_symmetric(DO_FORCE);

                                L.reset_offsets();
                                result.reset_offsets();
                            }
                        }

                        break;
                    }

                    case FACT_NONE:
                    default:
                    {
                        L_THROW(0);

                        break;
                    }
                }
            }

            break;
        }

        case FORM_TWO:
        {
            result.make_symmetric(DO_FORCE);
            result.pad_matrix(get_size()+1);

            if ( get_size() > 0 )
            {
                switch ( f_type )
                {
                    case FACT_INVERSE:
                    {
                        result = L;

                        if ( get_size()-nbad > 0 )
                        {
                            result.set_offset_start_at_row(1);
                            result.set_offset_start_at_col(1);
                            result.set_offset_end_at_row(get_size()-nbad+1);
                            result.set_offset_end_at_col(get_size()-nbad+1);

                            L.set_offset_start_at_row(1);
                            L.set_offset_start_at_col(1);
                            L.set_offset_end_at_row(get_size()-nbad+1);
                            L.set_offset_end_at_col(get_size()-nbad+1);

                            result.invert_matrix();

                            L.reset_offsets();

                            result.reset_offsets();
                        }

                        break;
                    }

                    case FACT_CHOL:
                    {
                        result = L;

                        if ( get_size()-nbad > 0 )
                        {
                            if ( !case_four_flag )
                            {
                                fMATRIX M;

                                L.set_offset_start_at_row(1);
                                L.set_offset_start_at_col(1);
                                L.set_offset_end_at_row(get_size()-nbad+1);
                                L.set_offset_end_at_col(get_size()-nbad+1);

                                M = L;

                                M.transpose();

                                for ( i = 1 ; i <= get_size()-nbad+1 ; i++ )
                                {
                                    M(1,i-1) *= -1.0;
                                }

                                L.reset_offsets();

                                result.set_offset_start_at_row(1);
                                result.set_offset_start_at_col(1);
                                result.set_offset_end_at_row(get_size()-nbad+1);
                                result.set_offset_end_at_col(get_size()-nbad+1);

                                L.set_offset_start_at_row(1);
                                L.set_offset_start_at_col(1);
                                L.set_offset_end_at_row(get_size()-nbad+1);
                                L.set_offset_end_at_col(get_size()-nbad+1);

                                result =  L;
                                result *= M;

                                result.make_symmetric(DO_FORCE);

                                L.reset_offsets();
                                result.reset_offsets();

                                if ( nbad > 0 )
                                {
                                    result.set_offset_start_at_row(get_size()-nbad+2);
                                    result.set_offset_start_at_col(1);
                                    result.set_offset_end_at_row(get_size()+1);
                                    result.set_offset_end_at_col(get_size()-nbad+1);

                                    L.set_offset_start_at_row(get_size()-nbad+2);
                                    L.set_offset_start_at_col(1);
                                    L.set_offset_end_at_row(get_size()+1);
                                    L.set_offset_end_at_col(get_size()-nbad+1);

                                    result =  L;
                                    result *= M;

                                    result.make_symmetric(DO_FORCE);

                                    L.reset_offsets();
                                    result.reset_offsets();
                                }
                            }
                        }

                        result.make_symmetric(DO_FORCE);

                        result.squareswap(1,2);

                        if ( primary_swap_flag )
                        {
                            result.squareswap(2,3);
                        }

                        break;
                    }

                    case FACT_NONE:
                    default:
                    {
                        L_THROW(0);

                        break;
                    }
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
// Internal functions:
//

long fMATRFact::xfact(const lVECTOR &Gindex, const fMATRIX &G_real)
{
    FN_ENTRY_POINT

    if ( f_type == FACT_NONE )
    {
        L_THROW(0);
    }

    long i,j;

    switch ( f_form )
    {
        case FORM_ONE:
        {
            repetition_point:

            if ( nbad > 0 )
            {
                switch ( f_type )
                {
                    case FACT_INVERSE:
                    {
                        if ( get_size() > 0 )
                        {
                            if ( get_size() == nbad )
                            {
                                s.make_zero(DO_FORCE);
                                s_valid = 0;

                                L_DOUBLE f;

                                try
                                {
                                    // Now comes the crunch

                                    f = L(0,0);
                                    NZ_ASSERT(f,zt);
                                    f = 1.0/f;

                                    // If we got past that, we can continue

                                    L(0,0) = f;

                                    nbad--;
                                }

                                catch ( int err_num )
                                {
                                    switch ( err_num )
                                    {
                                        case INV_ERR:
                                        {
                                            goto exit_point;

                                            break;
                                        }

                                        default:
                                        {
                                            throw err_num;

                                            break;
                                        }
                                    }
                                }
                            }

                            else
                            {
                                L_DOUBLE f;
                                fVECTOR q;
                                fVECTOR b;

                                L.set_offset_start_at_row(get_size()-nbad+1);
                                L.set_offset_start_at_col(1);
                                L.set_offset_end_at_row(get_size()-nbad+1);
                                L.set_offset_end_at_col(get_size()-nbad);

                                b = L[0];

                                L.reset_offsets();

                                if ( s_valid )
                                {
                                    q = s;
                                }

                                else
                                {
                                    L.set_offset_start_at_row(1);
                                    L.set_offset_start_at_col(1);
                                    L.set_offset_end_at_row(get_size()-nbad);
                                    L.set_offset_end_at_col(get_size()-nbad);

                                    q = ( b * L );

                                    L.reset_offsets();
                                }

                                f = ( L(get_size()-nbad,get_size()-nbad) - ( q * b ) );

                                try
                                {
                                    // Now comes the crunch

                                    NZ_ASSERT(f,zt);
                                    f = 1.0/f;

                                    // If we got past that, we can continue

                                    L(get_size()-nbad,get_size()-nbad) = f;

                                    f *= -1.0;

                                    L.set_offset_start_at_row(get_size()-nbad+1);
                                    L.set_offset_start_at_col(1);
                                    L.set_offset_end_at_row(get_size()-nbad+1);
                                    L.set_offset_end_at_col(get_size()-nbad);

                                    for ( j = 1 ; j <= L.get_effective_width() ; j++ )
                                    {
                                        L(0,j-1) = ( f * q.get_offset_element(j-1) );
                                    }

                                    L.reset_offsets();

                                    L.set_offset_start_at_row(1);
                                    L.set_offset_start_at_col(1);
                                    L.set_offset_end_at_row(get_size()-nbad);
                                    L.set_offset_end_at_col(get_size()-nbad);

                                    L(ARG_RANK_ONE,q,-f);

                                    L.reset_offsets();

                                    nbad--;

                                    s.make_zero(DO_FORCE);
                                    s_valid = 0;
                                }

                                catch ( int err_num )
                                {
                                    switch ( err_num )
                                    {
                                        case INV_ERR:
                                        {
                                            s.make_zero(DO_FORCE);
                                            s_valid = 0;

                                            s = q;
                                            s_valid = 1;

                                            goto exit_point;

                                            break;
                                        }

                                        default:
                                        {
                                            throw err_num;

                                            break;
                                        }
                                    }
                                }
                            }
                        }

                        else
                        {
                            s.make_zero(DO_FORCE);
                            s_valid = 0;
                        }

                        break;
                    }

                    case FACT_CHOL:
                    {
                        if ( get_size() > 0 )
                        {
                            if ( get_size() == nbad )
                            {
                                s.make_zero(DO_FORCE);
                                s_valid = 0;

                                L_DOUBLE f;

                                try
                                {
                                    f = L(0,0);
                                    POS_ASSERT(L(0,0),zt);
                                    f = sqrt(f);

                                    // success - we have reached here

                                    L(0,0) = f;

                                    // and now the flow down.

                                    f = 1.0/f;

                                    if ( nbad > 1 )
                                    {
                                        for ( i = 2 ; i <= get_size() ; i++ )
                                        {
                                            L(i-1,get_size()-nbad) *= f;
                                        }
                                    }

                                    nbad--;
                                }

                                catch ( int err_num )
                                {
                                    switch ( err_num )
                                    {
                                        case SQRT_ERR:
                                        {
                                            goto exit_point;

                                            break;
                                        }

                                        default:
                                        {
                                            throw err_num;

                                            break;
                                        }
                                    }
                                }
                            }

                            else
                            {
                                fVECTOR b;
                                L_DOUBLE f;
                                L_DOUBLE g;

                                L.set_offset_start_at_row(get_size()-nbad+1);
                                L.set_offset_start_at_col(1);
                                L.set_offset_end_at_row(get_size()-nbad+1);
                                L.set_offset_end_at_col(get_size()-nbad);

                                b = L[0];

                                f = -( L[0] * b );

                                L.reset_offsets();

                                f += L(get_size()-nbad,get_size()-nbad);

                                try
                                {
                                    POS_ASSERT(f,zt);
                                    f = sqrt(f);

                                    // success - we have reached here

                                    L(get_size()-nbad,get_size()-nbad) = f;

                                    if ( nbad > 1 )
                                    {
                                        // and now the flow down.

                                        f = 1.0/f;

                                        for ( i = get_size()-nbad+2 ; i <= get_size() ; i++ )
                                        {
                                            g = L(i-1,get_size()-nbad);

                                            L.set_offset_start_at_row(i);
                                            L.set_offset_start_at_col(1);
                                            L.set_offset_end_at_row(i);
                                            L.set_offset_end_at_col(get_size()-nbad);

                                            g -= ( L[0] * b );

                                            L.reset_offsets();

                                            g *= f;

                                            L(i-1,get_size()-nbad) = g;
                                        }
                                    }

                                    nbad--;

                                    s.make_zero(DO_FORCE);
                                    s_valid = 0;
                                }

                                catch ( int err_num )
                                {
                                    switch ( err_num )
                                    {
                                        case SQRT_ERR:
                                        {
                                            goto exit_point;

                                            break;
                                        }

                                        default:
                                        {
                                            throw err_num;

                                            break;
                                        }
                                    }
                                }
                            }
                        }

                        else
                        {
                            s.make_zero(DO_FORCE);
                            s_valid = 0;
                        }

                        break;
                    }

                    case FACT_NONE:
                    default:
                    {
                        L_THROW(0);

                        break;
                    }
                }

                if ( nbad != 0 )
                {
                    goto repetition_point;
                }
            }

            exit_point:

            break;
        }

        case FORM_TWO:
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

    fix_symmetry();

    FN_EXIT_POINT nbad;

    L_DOUBLE tempa;
    long tempb;

    tempa = G_real.get_offset_element(-1,-1);
    tempb = Gindex.get_offset_element(-1);
}

long fMATRFact::xfact(const lVECTOR &Gindex, const fMATRIX &G_real, const fVECTOR &d_real)
{
    FN_ENTRY_POINT

    if ( f_type == FACT_NONE )
    {
        L_THROW(0);
    }

    long i,j;

    switch ( f_form )
    {
        case FORM_ONE:
        {
            L_THROW(0);

            break;
        }

        case FORM_TWO:
        {
            repetition_point:

            if ( nbad > 0 )
            {
                switch ( f_type )
                {
                    case FACT_INVERSE:
                    {
                        if ( get_size() > 0 )
                        {
                            if ( get_size() == nbad )
                            {
                                s.make_zero(DO_FORCE);
                                s_valid = 0;


                                L.remove(get_size()+1);
                                overwrite_L(Gindex,G_real); // modded

                                refresh_fact(d_real);
                            }

                            else
                            {
                                L_DOUBLE f;
                                fVECTOR q;
                                fVECTOR b;

                                L.set_offset_start_at_row(get_size()-nbad+2);
                                L.set_offset_start_at_col(1);
                                L.set_offset_end_at_row(get_size()-nbad+2);
                                L.set_offset_end_at_col(get_size()-nbad+1);

                                b = L[0];

                                L.reset_offsets();

                                if ( s_valid )
                                {
                                    q = s;
                                }

                                else
                                {
                                    L.set_offset_start_at_row(1);
                                    L.set_offset_start_at_col(1);
                                    L.set_offset_end_at_row(get_size()-nbad+1);
                                    L.set_offset_end_at_col(get_size()-nbad+1);

                                    q = ( b * L );

                                    L.reset_offsets();
                                }

                                f = ( L(get_size()-nbad+1,get_size()-nbad+1) - ( q * b ) );

                                try
                                {
                                    // Now comes the crunch

                                    NZ_ASSERT(f,zt);
                                    f = 1.0/f;

                                    // If we got past that, we can continue

                                    L(get_size()-nbad+1,get_size()-nbad+1) = f;

                                    f *= -1.0;

                                    L.set_offset_start_at_row(get_size()-nbad+2);
                                    L.set_offset_start_at_col(1);
                                    L.set_offset_end_at_row(get_size()-nbad+2);
                                    L.set_offset_end_at_col(get_size()-nbad+1);

                                    for ( j = 1 ; j <= L.get_effective_width() ; j++ )
                                    {
                                        L(0,j-1) = ( f * q.get_offset_element(j-1) );
                                    }

                                    L.reset_offsets();

                                    L.set_offset_start_at_row(1);
                                    L.set_offset_start_at_col(1);
                                    L.set_offset_end_at_row(get_size()-nbad+1);
                                    L.set_offset_end_at_col(get_size()-nbad+1);

                                    L(ARG_RANK_ONE,q,-f);

                                    L.reset_offsets();

                                    nbad--;

                                    s.make_zero(DO_FORCE);
                                    s_valid = 0;
                                }

                                catch ( int err_num )
                                {
                                    switch ( err_num )
                                    {
                                        case INV_ERR:
                                        {
                                            s.make_zero(DO_FORCE);
                                            s_valid = 0;

                                            s = q;
                                            s_valid = 1;

                                            goto exit_point;

                                            break;
                                        }

                                        default:
                                        {
                                            throw err_num;

                                            break;
                                        }
                                    }
                                }
                            }
                        }

                        else
                        {
                            s.make_zero(DO_FORCE);
                            s_valid = 0;
                        }

                        break;
                    }

                    case FACT_CHOL:
                    {
                        if ( case_four_flag )
                        {
                            goto exit_point;
                        }

                        if ( get_size() > 0 )
                        {
                            if ( get_size() == nbad )
                            {
                                s.make_zero(DO_FORCE);
                                s_valid = 0;


                                L.remove(get_size()+1);
                                overwrite_L(Gindex,G_real); // modded

                                refresh_fact(d_real);
                            }

                            else
                            {
                                fVECTOR dhm;
                                L_DOUBLE q;
                                L_DOUBLE r;

                                L.set_offset_start_at_row(get_size()-nbad+2);
                                L.set_offset_start_at_col(1);
                                L.set_offset_end_at_row(get_size()-nbad+2);
                                L.set_offset_end_at_col(get_size()-nbad+1);

                                dhm = L[0];

                                dhm[1] *= -1.0;

                                q = -( L[0] * dhm );

                                L.reset_offsets();

                                q += L(get_size()-nbad+1,get_size()-nbad+1);

                                try
                                {
                                    POS_ASSERT(q,zt);
                                    q = sqrt(q);

                                    // success - we have reached here

                                    L(get_size()-nbad+1,get_size()-nbad+1) = q;

                                    // and now the flow down.

                                    q = 1.0/q;

                                    if ( nbad > 1 )
                                    {
                                        for ( i = get_size()-nbad+3 ; i <= get_size()+1 ; i++ )
                                        {
                                            r = L(i-1,get_size()-nbad+1);

                                            L.set_offset_start_at_row(i);
                                            L.set_offset_start_at_col(1);
                                            L.set_offset_end_at_row(i);
                                            L.set_offset_end_at_col(get_size()-nbad+1);

                                            r -= ( L[0] * dhm );

                                            L.reset_offsets();

                                            r *= q;

                                            L(i-1,get_size()-nbad+1) = r;
                                        }
                                    }

                                    nbad--;

                                    s.make_zero(DO_FORCE);
                                    s_valid = 0;
                                }

                                catch ( int err_num )
                                {
                                    switch ( err_num )
                                    {
                                        case SQRT_ERR:
                                        {
                                            goto exit_point;

                                            break;
                                        }

                                        default:
                                        {
                                            throw err_num;

                                            break;
                                        }
                                    }
                                }
                            }
                        }

                        else
                        {
                            s.make_zero(DO_FORCE);
                            s_valid = 0;
                        }

                        break;
                    }

                    case FACT_NONE:
                    default:
                    {
                        L_THROW(0);

                        break;
                    }
                }

                if ( nbad != 0 )
                {
                    goto repetition_point;
                }
            }

            exit_point:

            break;
        }

        default:
        {
            L_THROW(0);

            break;
        }
    }

    fix_symmetry();

    FN_EXIT_POINT nbad;
}

long fMATRFact::xrankone(const lVECTOR &Gindex, const fVECTOR &ax, const L_DOUBLE &bx, const fMATRIX &G_real, int hold_off, long z_start, long z_end)
{
    FN_ENTRY_POINT

    if ( f_type == FACT_NONE )
    {
        L_THROW(0);
    }

    long i,j,k;
    long z_lap;
    long size;
    long refresh_size;
    int fresh_start;

    // Matrix inversion lemma:
    //
    // inv( A - B.D.C ) = inv(A) + inv(A).B.inv( inv(D) - C.inv(A).B ).C.inv(A)
    // inv( A + B.D.C ) = inv(A) - inv(A).B.inv( inv(D) + C.inv(A).B ).C.inv(A)

    // Assume a is a normal vector

    from_the_top:

    size = get_size();

    THROW_ASSERT(size == ax.get_effective_size());

    if ( ( bx != 0.0 ) && ( size-z_start-z_end > 0 ) )
    {
        fVECTOR a('x',size);
        L_DOUBLE b;

        a = ax;
        b = bx;

        switch ( f_form )
        {
            case FORM_ONE:
            {
                switch ( f_type )
                {
                    case FACT_INVERSE:
                    {
                        s.make_zero(DO_FORCE);
                        s_valid = 0;

                        // Step 1: update the factorised corner

                        if ( get_size()-nbad > z_start )
                        {
                            z_lap = nbad;

                            if ( z_end > nbad )
                            {
                                z_lap = z_end;
                            }

                            fVECTOR q;
                            L_DOUBLE f;

                            L.set_offset_start_at_row(1);
                            L.set_offset_start_at_col(1+z_start);
                            L.set_offset_end_at_row(get_size()-nbad);
                            L.set_offset_end_at_col(get_size()-z_lap);

                            a.set_offset_start(1+z_start);
                            a.set_offset_end(get_size()-z_lap);

                            q = ( L * a );

                            a.reset_offsets();
                            L.reset_offsets();

                            a.set_offset_start(1);
                            a.set_offset_end(get_size()-nbad);

                            f = ( (1.0/b) + (q*a) );

                            a.reset_offsets();

                            try
                            {
                                NZ_ASSERT(f,zt);
                                f = -1.0/f;

                                L.set_offset_start_at_row(1);
                                L.set_offset_start_at_col(1);
                                L.set_offset_end_at_row(get_size()-nbad);
                                L.set_offset_end_at_col(get_size()-nbad);

                                L(ARG_RANK_ONE,q,f);

                                L.reset_offsets();
                            }

                            catch ( int err_num )
                            {
                                switch ( err_num )
                                {
                                    case INV_ERR:
                                    {
                                        refresh_size = get_size()-nbad;
                                        fresh_start  = 1;

                                        if ( (get_size()-nbad)/2 >= 2 )
                                        {
                                            for ( j = 1 ; j <= (get_size()-nbad)/2 ; j++ )
                                            {
                                                try
                                                {
                                                    fMATRIX tempb;

                                                    L.set_offset_start_at_row(get_size()-nbad-j+1);
                                                    L.set_offset_start_at_col(get_size()-nbad-j+1);
                                                    L.set_offset_end_at_row(get_size()-nbad);
                                                    L.set_offset_end_at_col(get_size()-nbad);

                                                    tempb = -L;
                                                    tempb.invert_matrix(zt);

                                                    L.reset_offsets();

                                                    fMATRIX tempc;

                                                    L.set_offset_start_at_col(get_size()-nbad-j+1);
                                                    L.set_offset_start_at_row(1);
                                                    L.set_offset_end_at_col(get_size()-nbad);
                                                    L.set_offset_end_at_row(get_size()-nbad-j);

                                                    tempc = L;




                                                    L.reset_offsets();

                                                    L.set_offset_start_at_row(1);
                                                    L.set_offset_start_at_col(1);
                                                    L.set_offset_end_at_row(get_size()-nbad-j);
                                                    L.set_offset_end_at_col(get_size()-nbad-j);

                                                    L(ARG_RANK_N,tempc,tempb);

                                                    L.reset_offsets();

                                                    refresh_size = j;
                                                    fresh_start  = 0;

                                                    goto get_out_a;

                                                    break;
                                                }

                                                catch ( int err_num )
                                                {
                                                    switch ( err_num )
                                                    {
                                                        case INV_ERR:
                                                        {
                                                            L.reset_offsets();

                                                            break;
                                                        }

                                                        default:
                                                        {
                                                            throw err_num;

                                                            break;
                                                        }
                                                    }
                                                }
                                            }
                                        }

                                        get_out_a:

                                        if ( refresh_size >= 1 )
                                        {
                                            for ( j = get_size()-nbad+1-refresh_size ; j <= get_size()-nbad ; j++ )
                                            {
                                                for ( k = 1 ; k <= j ; k++ )
                                                {
                                                    L(j-1,k-1) = get_reordered_value(Gindex,G_real,j,k);
                                                }

                                                nbad++;
                                            }
                                        }

                                        if ( fresh_start )
                                        {
                                            // In this case size(L) = size(G)

                                            overwrite_L(Gindex,G_real); // modded

                                            refresh_fact();
                                        }

                                        goto from_the_top;

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

                        // Step 2. update unfactorised corner

                        if ( ( nbad > z_end ) && ( hold_off == 0 ) )
                        {
                            L_DOUBLE temp;
                            L_DOUBLE tempb;

                            for ( i = get_size()-nbad+1 ; i <= get_size()-z_end ; i++ )
                            {
                                for ( j = 1 ; j <= i ; j++ )
                                {
                                    tempb = a[j-1];

                                    temp  = a[i-1];
                                    temp *= b;
                                    temp *= tempb;

                                    L(i-1,j-1) += temp;
                                }
                            }
                        }

                        break;
                    }

                    case FACT_CHOL:
                    {
                        s.make_zero(DO_FORCE);
                        s_valid = 0;

                        L_DOUBLE x;
                        L_DOUBLE y;
                        L_DOUBLE alpha;
                        L_DOUBLE beta;
                        L_DOUBLE gamma;
                        L_DOUBLE epsilon;

                        // Step 1: update the factorised corner

                        if ( b >= 0.0 )
                        {
                            a *= sqrt(b);

                            b /= b;
                        }

                        else
                        {
                            a *= sqrt(-b);

                            b /= -b;
                        }

                        if ( get_size()-nbad > z_start )
                        {
                            fVECTOR aaz;

                            aaz = a;

                            for ( i = 1+z_start ; i <= get_size()-nbad ; i++ )
                            {
                                x =  L(i-1,i-1);
                                x *= L(i-1,i-1);

                                y =  aaz[i-1];
                                y *= aaz[i-1];

                                if ( b > 0.0 )
                                {
                                    x += y;
                                }

                                else
                                {
                                    x -= y;
                                }

                                try
                                {
                                    POS_ASSERT(x,zt);
                                    x = sqrt(x);

                                    alpha =  aaz[i-1];
                                    alpha /= x;

                                    beta =  L(i-1,i-1);
                                    beta /= x;

                                    if ( i < get_size() )
                                    {
                                        for ( j = i+1 ; j <= get_size() ; j++ )
                                        {
                                            gamma   = aaz[j-1];
                                            epsilon = L(j-1,i-1);

                                            aaz[j-1] =  ( alpha * epsilon );
                                            aaz[j-1] -= ( beta  * gamma   );

                                            if ( b >= 0.0 )
                                            {
                                                L(j-1,i-1) =  ( beta  * epsilon );
                                                L(j-1,i-1) += ( alpha * gamma   );
                                            }

                                            else
                                            {
                                                L(j-1,i-1) =  ( beta  * epsilon );
                                                L(j-1,i-1) -= ( alpha * gamma   );
                                            }
                                        }
                                    }

                                    L(i-1,i-1) = x;
                                }

                                catch ( int err_num )
                                {
                                    switch ( err_num )
                                    {
                                        case SQRT_ERR:
                                        {
                                            refresh_size = get_size()-nbad-i+1;

                                            for ( j = get_size()-nbad+1-refresh_size ; j <= get_size()-nbad ; j++ )
                                            {
                                                for ( k = get_size()-nbad+1-refresh_size ; k <= j ; k++ )
                                                {
                                                    L(j-1,k-1) = get_reordered_value(Gindex,G_real,j,k);
                                                }

                                                nbad++;
                                            }

                                            goto keep_going;

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

                        // Step 2. update unfactorised corner

                        keep_going:

                        if ( ( nbad > z_end ) && ( hold_off == 0 ) )
                        {
                            L_DOUBLE temp;
                            L_DOUBLE tempb;

                            for ( i = get_size()-nbad+1 ; i <= get_size()-z_end ; i++ )
                            {
                                for ( j = get_size()-nbad+1 ; j <= i ; j++ )
                                {
                                    tempb = a[j-1];

                                    temp  = a[i-1];
                                    temp *= b;
                                    temp *= tempb;

                                    L(i-1,j-1) += temp;
                                }
                            }
                        }

                        break;
                    }

                    case FACT_NONE:
                    default:
                    {
                        L_THROW(0);

                        break;
                    }
                }

                break;
            }

            case FORM_TWO:
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

    fix_symmetry();

    FN_EXIT_POINT nbad;
}

long fMATRFact::xrankone(const lVECTOR &Gindex, const fVECTOR &ax, const L_DOUBLE &bx, const fMATRIX &G_real, const fVECTOR &d_real, int hold_off, long z_start, long z_end)
{
    FN_ENTRY_POINT

    if ( f_type == FACT_NONE )
    {
        L_THROW(0);
    }

    long i,j,k;
    long size;
    long refresh_size;
    int fresh_start;
    long z_lap;

    // Matrix inversion lemma:
    //
    // inv( A - B.D.C ) = inv(A) + inv(A).B.inv( inv(D) - C.inv(A).B ).C.inv(A)
    // inv( A + B.D.C ) = inv(A) - inv(A).B.inv( inv(D) + C.inv(A).B ).C.inv(A)

    // Assume a is a normal vector

    from_the_top:

    size  = get_size();

    THROW_ASSERT((size+1) == ax.get_effective_size());

    if ( ( bx != 0.0 ) && ( size-z_start-z_end > 0 ) )
    {
        fVECTOR a('x',size+1);
        L_DOUBLE b;

        a = ax;
        b = bx;

        switch ( f_form )
        {
            case FORM_ONE:
            {
                L_THROW(0);

                break;
            }

            case FORM_TWO:
            {
                switch ( f_type )
                {
                    case FACT_INVERSE:
                    {
                        s.make_zero(DO_FORCE);
                        s_valid = 0;

                        // Step 1: update the factorised corner

                        if ( get_size()-nbad > z_start )
                        {
                            z_lap = nbad;

                            if ( z_end > nbad )
                            {
                                z_lap = z_end;
                            }

                            fVECTOR q;
                            L_DOUBLE f;

                            L.set_offset_start_at_row(1);
                            L.set_offset_start_at_col(1+z_start);
                            L.set_offset_end_at_row(get_size()-nbad+1);
                            L.set_offset_end_at_col(get_size()-z_lap+1);

                            a.set_offset_start(1+z_start);
                            a.set_offset_end(get_size()-z_lap+1);

                            q = ( L * a );

                            a.reset_offsets();
                            L.reset_offsets();

                            a.set_offset_start(1);
                            a.set_offset_end(get_size()-nbad+1);

                            f = ( (1.0/b) + (q*a) );

                            a.reset_offsets();

                            try
                            {
                                NZ_ASSERT(f,zt);
                                f = -1.0/f;

                                L.set_offset_start_at_row(1);
                                L.set_offset_start_at_col(1);
                                L.set_offset_end_at_row(get_size()-nbad+1);
                                L.set_offset_end_at_col(get_size()-nbad+1);

                                L(ARG_RANK_ONE,q,f);

                                L.reset_offsets();
                            }

                            catch ( int err_num )
                            {
                                switch ( err_num )
                                {
                                    case INV_ERR:
                                    {
                                        refresh_size = get_size()-nbad+1;
                                        fresh_start  = 1;

                                        if ( (get_size()-nbad)/2 >= 2 )
                                        {
                                            for ( j = 1 ; j <= (get_size()-nbad)/2 ; j++ )
                                            {
                                                try
                                                {
                                                    fMATRIX tempb;

                                                    L.set_offset_start_at_row(get_size()-nbad-j+2);
                                                    L.set_offset_start_at_col(get_size()-nbad-j+2);
                                                    L.set_offset_end_at_row(get_size()-nbad+1);
                                                    L.set_offset_end_at_col(get_size()-nbad+1);

                                                    tempb = -L;
                                                    tempb.invert_matrix(zt);

                                                    L.reset_offsets();

                                                    fMATRIX tempc;

                                                    L.set_offset_start_at_col(get_size()-nbad-j+2);
                                                    L.set_offset_start_at_row(1);
                                                    L.set_offset_end_at_col(get_size()-nbad+1);
                                                    L.set_offset_end_at_row(get_size()-nbad-j+1);

                                                    tempc = L;




                                                    L.reset_offsets();

                                                    L.set_offset_start_at_row(1);
                                                    L.set_offset_start_at_col(1);
                                                    L.set_offset_end_at_row(get_size()-nbad-j+1);
                                                    L.set_offset_end_at_col(get_size()-nbad-j+1);

                                                    L(ARG_RANK_N,tempc,tempb);

                                                    L.reset_offsets();

                                                    refresh_size = j;
                                                    fresh_start  = 0;

                                                    goto get_out_a;

                                                    break;
                                                }

                                                catch ( int err_num )
                                                {
                                                    switch ( err_num )
                                                    {
                                                        case INV_ERR:
                                                        {
                                                            L.reset_offsets();

                                                            break;
                                                        }

                                                        default:
                                                        {
                                                            throw err_num;

                                                            break;
                                                        }
                                                    }
                                                }
                                            }
                                        }

                                        get_out_a:

                                        if ( refresh_size >= 1 )
                                        {
                                            long mbad;

                                            mbad = nbad;

                                            if ( fresh_start )
                                            {
                                                L(0,0) = 0.0;
                                            }

                                            for ( j = (get_size()-mbad)+1-refresh_size+1 ; j <= (get_size()-mbad)+1 ; j++ )
                                            {
                                                for ( k = 1+1 ; k <= j ; k++ )
                                                {
                                                    L(j-1,k-1) = get_reordered_value(Gindex,G_real,j-1,k-1);
                                                }

                                                nbad++;
                                            }
                                        }

                                        if ( fresh_start )
                                        {
                                            long d_size = get_size();


                                            L.remove(d_size+1);
                                            overwrite_L(Gindex,G_real); // modded

                                            exit_fixvect_normal_fn(a);

                                            a.remove(1);

                                            L(ARG_RANK_ONE,a,b);

                                            refresh_fact(d_real);

                                            goto quick_exit;
                                        }

                                        goto from_the_top;

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

                        // Step 2. update unfactorised corner

                        if ( ( nbad > z_end ) && ( hold_off == 0 ) )
                        {
                            L_DOUBLE temp;
                            L_DOUBLE tempb;

                            for ( i = get_size()-nbad+2 ; i <= get_size()-z_end+1 ; i++ )
                            {
                                for ( j = 1 ; j <= i ; j++ )
                                {
                                    tempb = a[j-1];

                                    temp  = a[i-1];
                                    temp *= b;
                                    temp *= tempb;

                                    L(i-1,j-1) += temp;
                                }
                            }
                        }

                        break;
                    }

                    case FACT_CHOL:
                    {
                        if ( case_four_flag )
                        {
                            long d_size = get_size();


                            L.remove(d_size+1);
                            overwrite_L(Gindex,G_real); // modded

                            exit_fixvect_normal_fn(a);

                            a.remove(1);

                            L(ARG_RANK_ONE,a,b);

                            refresh_fact(d_real);

                            goto quick_exit;
                        }

                        else
                        {
                            int b_sign;

                            s.make_zero(DO_FORCE);
                            s_valid = 0;

                            L_DOUBLE x;
                            L_DOUBLE y;
                            L_DOUBLE alpha;
                            L_DOUBLE beta;
                            L_DOUBLE gamma;
                            L_DOUBLE epsilon;

                            long s_point;

                            s_point = 1;

                            if ( z_start > 0 )
                            {
                                s_point = 2+z_start;
                            }

                            // Step 1: update the factorised corner

                            if ( b > 0.0 )
                            {
                                a *= sqrt(b);
                                b /= b;

                                b_sign = +1;
                            }

                            else
                            {
                                a *= sqrt(-b);
                                b /= -b;

                                b_sign = -1;
                            }

                            if ( s_point <= get_size()-nbad+1 )
                            {
                                fVECTOR aaz;

                                aaz = a;

                                for ( i = s_point ; i <= get_size()-nbad+1 ; i++ )
                                {
                                    x =  L(i-1,i-1);
                                    x *= L(i-1,i-1);

                                    y =  aaz[i-1];
                                    y *= aaz[i-1];

                                    if ( b_sign == +1 )
                                    {
                                        if ( i != 2 )
                                        {
                                            x += y;
                                        }

                                        else
                                        {
                                            x -= y;
                                        }
                                    }

                                    else
                                    {
                                        if ( i != 2 )
                                        {
                                            x -= y;
                                        }

                                        else
                                        {
                                            x += y;
                                        }
                                    }

                                    try
                                    {
                                        POS_ASSERT(x,zt);
                                        x = sqrt(x);

                                        alpha =  aaz[i-1];
                                        alpha /= x;

                                        beta =  L(i-1,i-1);
                                        beta /= x;

                                        if ( i < get_size()+1 )
                                        {
                                            for ( j = i+1 ; j <= get_size()+1 ; j++ )
                                            {
                                                gamma   = aaz[j-1];
                                                epsilon = L(j-1,i-1);

                                                aaz[j-1] =  ( alpha * epsilon );
                                                aaz[j-1] -= ( beta  * gamma   );

                                                if ( b_sign == +1 )
                                                {
                                                    if ( i != 2 )
                                                    {
                                                        L(j-1,i-1) =  ( beta  * epsilon );
                                                        L(j-1,i-1) += ( alpha * gamma   );
                                                    }

                                                    else
                                                    {
                                                        L(j-1,i-1) =  ( beta  * epsilon );
                                                        L(j-1,i-1) -= ( alpha * gamma   );
                                                    }
                                                }

                                                else
                                                {
                                                    if ( i != 2 )
                                                    {
                                                        L(j-1,i-1) =  ( beta  * epsilon );
                                                        L(j-1,i-1) -= ( alpha * gamma   );
                                                    }

                                                    else
                                                    {
                                                        L(j-1,i-1) =  ( beta  * epsilon );
                                                        L(j-1,i-1) += ( alpha * gamma   );
                                                    }
                                                }
                                            }
                                        }

                                        L(i-1,i-1) = x;
                                    }

                                    catch ( int err_num )
                                    {
                                        switch ( err_num )
                                        {
                                            case SQRT_ERR:
                                            {
                                                if ( i <= 2 )
                                                {
                                                    long d_size = get_size();


                                                    L.remove(d_size+1);
                                                    overwrite_L(Gindex,G_real); // modded

                                                    exit_fixvect_normal_fn(a);

                                                    a.remove(1);

                                                    L(ARG_RANK_ONE,a,b);

                                                    refresh_fact(d_real);

                                                    goto quick_exit;
                                                }

                                                else
                                                {
                                                    nbad = 0;

                                                    for ( j = i ; j <= get_size()+1 ; j++ )
                                                    {
                                                        for ( k = i ; k <= j ; k++ )
                                                        {
                                                            L(j-1,k-1) = get_reordered_value(Gindex,G_real,j-1,k-1);
                                                        }

                                                        nbad++;
                                                    }
                                                }

                                                goto keep_going;

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

                            // Step 2. update unfactorised corner

                            keep_going:

                            if ( ( nbad > z_end ) && ( hold_off == 0 ) )
                            {
                                L_DOUBLE temp;
                                L_DOUBLE tempb;

                                for ( i = get_size()-nbad+2 ; i <= get_size()-z_end+1 ; i++ )
                                {
                                    for ( j = get_size()-nbad+2 ; j <= i ; j++ )
                                    {
                                        tempb = a[j-1];

                                        temp  = a[i-1];
                                        temp *= b;
                                        temp *= tempb;

                                        L(i-1,j-1) += temp;
                                    }
                                }
                            }
                        }

                        break;
                    }

                    case FACT_NONE:
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

    quick_exit:

    fix_symmetry();

    FN_EXIT_POINT nbad;
}

long fMATRFact::refresh_fact(void)
{
    FN_ENTRY_POINT

    if ( f_type == FACT_NONE )
    {
        L_THROW(0);
    }

    switch ( f_form )
    {
        case FORM_ONE:
        {
            s.make_zero(DO_FORCE);
            s_valid = 0;

            switch ( f_type )
            {
                case FACT_INVERSE:
                {
                    L.make_symmetric(DO_FORCE);

                    nbad = get_size();

                    break;
                }

                case FACT_CHOL:
                {
                    L.make_symmetric(DO_FORCE);
                    L.make_lower_triangular(DO_FORCE);

                    nbad = get_size();

                    break;
                }

                case FACT_NONE:
                default:
                {
                    L_THROW(0);

                    break;
                }
            }

            break;
        }

        case FORM_TWO:
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

    FN_EXIT_POINT nbad;
}

long fMATRFact::refresh_fact(const fVECTOR &d_real)
{
    FN_ENTRY_POINT

    if ( f_type == FACT_NONE )
    {
        L_THROW(0);
    }

    long i;
    int is_sqrt_done;

    switch ( f_form )
    {
        case FORM_ONE:
        {
            L_THROW(0);

            break;
        }

        case FORM_TWO:
        {
            s.make_zero(DO_FORCE);
            s_valid = 0;

            switch ( f_type )
            {
                case FACT_INVERSE:
                {
                    L.make_symmetric(DO_FORCE);

                    nbad = get_size();

                    break;
                }

                case FACT_CHOL:
                {
                    L.make_symmetric(DO_FORCE);
                    L.make_lower_triangular(DO_FORCE);

                    nbad = get_size();

                    break;
                }

                case FACT_NONE:
                default:
                {
                    L_THROW(0);

                    break;
                }
            }

            case_four_flag    = 0;
            primary_swap_flag = 0;

            if ( L.get_effective_height() > 0 )
            {
                switch ( f_type )
                {
                    case FACT_INVERSE:
                    {
                        fVECTOR loc_d('x',L.get_effective_height());

                        loc_d = d_real;
                        loc_d.addstart();

                        L.addstart(loc_d[0],loc_d);

                        loc_d.make_zero(DO_FORCE);
                        loc_d.make_normal(DO_FORCE);

                        nbad = get_size()-1;

                        L_DOUBLE tempa;
                        L_DOUBLE tempb;
                        L_DOUBLE tempc;

                        tempa =  L(1,0);
                        tempb =  1.0/tempa;
                        tempc =  -tempb;
                        tempc *= L(1,1);
                        tempc *= tempb;

                        L(0,0) = tempc;
                        L(1,0) = tempb;
                        L(1,1) = 0.0;

                        break;
                    }

                    case FACT_CHOL:
                    {
                        fVECTOR loc_d('x',L.get_effective_height());

                        loc_d = d_real;
                        loc_d.addstart();

                        L.addstart(loc_d[0],loc_d);

                        loc_d.make_zero(DO_FORCE);
                        loc_d.make_normal(DO_FORCE);

                        nbad = get_size()-1;

                        L_DOUBLE temp;
                        L_DOUBLE tempb;

                        temp = 0.0;

                        L.make_lower_triangular(DO_FORCE);

                        // We need to check for case_four_flag and
                        // primary_swap_flag start conditions.

                        L.squareswap(1,2);

                        is_sqrt_done = 0;

                        try_again:

                        try
                        {
                            if ( !is_sqrt_done )
                            {
                                POS_ASSERT(L(0,0),zt);
                                temp = sqrt(L(0,0));
                            }

                            is_sqrt_done = 0;

                            // If we have got this far, things are OK.

                            tempb = 1.0/temp;

                            L(0,0) = temp;
                            L(1,0) = L(1,0) * tempb;
                            L(1,1) = sqrt(L(1,0)*L(1,0));

                            nbad = get_size()-1;

                            // Flow down

                            L_DOUBLE tempe;
                            L_DOUBLE tempf;
                            L_DOUBLE tempg;
                            L_DOUBLE temph;

                            tempf =  1.0/L(1,1);
                            tempe =  -tempf;
                            tempf *= L(1,0);
                            tempf *= tempb;

                            for ( i = 2 ; i <= get_size() ; i++ )
                            {
                                tempg = L(i,0);
                                temph = L(i,1);

                                L(i,0) = ( tempg * tempb );
                                L(i,1) = ( tempg * tempf )
                                       + ( temph * tempe );
                            }
                        }

                        catch ( int err_num )
                        {
                            switch ( err_num )
                            {
                                case SQRT_ERR:
                                {
                                    if ( get_size() == 1 )
                                    {
                                        // case 4

                                        nbad = 0;

                                        case_four_flag = 1;

                                        break;
                                    }

                                    else
                                    {
                                        // Try switching (primary_swap_flag)

                                        try
                                        {
                                            POS_ASSERT(L(2,2),zt);
                                            temp = sqrt(L(2,2));

                                            // If we have got this far,
                                            // things are OK.

                                            primary_swap_flag = 1;

                                            L.squareswap(1,3);

                                            is_sqrt_done = 1;

                                            goto try_again;
                                        }

                                        catch ( int x_err )
                                        {
                                            switch ( x_err )
                                            {
                                                case SQRT_ERR:
                                                {
                                                    // case 4

                                                    nbad = get_size()-1;

                                                    case_four_flag = 1;

                                                    break;
                                                }

                                                default:
                                                {
                                                    throw x_err;

                                                    break;
                                                }
                                            }
                                        }

                                        break;
                                    }

                                    default:
                                    {
                                        throw err_num;

                                        break;
                                    }
                                }
                            }
                        }

                        break;
                    }

                    case FACT_NONE:
                    default:
                    {
                        L_THROW(0);

                        break;
                    }
                }
            }

            else
            {
                switch ( f_type )
                {
                    case FACT_INVERSE:
                    case FACT_CHOL:
                    {
                        L.addstart();

                        nbad = 0;

                        break;
                    }

                    case FACT_NONE:
                    default:
                    {
                        L_THROW(0);

                        break;
                    }
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

    fix_symmetry();

    FN_EXIT_POINT nbad;
}




//
// Reality <-> illusion conversion functions:
//

long fMATRFact::enter_fixvect_normal_fn(fVECTOR &a)
{
    FN_ENTRY_POINT

    switch ( f_form )
    {
        case FORM_ONE:
        {
            break;
        }

        case FORM_TWO:
        {
            switch ( f_type )
            {
                case FACT_NONE:
                case FACT_INVERSE:
                {
                    break;
                }

                case FACT_CHOL:
                {
                    if ( primary_swap_flag )
                    {
                        a.squareswap(2,3);
                    }

                    a.squareswap(1,2);

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

    FN_EXIT_POINT nbad;
}

long fMATRFact::exit_fixvect_normal_fn(fVECTOR &a)
{
    FN_ENTRY_POINT

    switch ( f_form )
    {
        case FORM_ONE:
        {
            break;
        }

        case FORM_TWO:
        {
            switch ( f_type )
            {
                case FACT_NONE:
                case FACT_INVERSE:
                {
                    break;
                }

                case FACT_CHOL:
                {
                    a.squareswap(1,2);

                    if ( primary_swap_flag )
                    {
                        a.squareswap(2,3);
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

        default:
        {
            L_THROW(0);

            break;
        }
    }

    FN_EXIT_POINT nbad;
}

long fMATRFact::convert_to_actual_pos(long virtual_pos)
{
    FN_ENTRY_POINT

    long result;

    switch ( f_form )
    {
        case FORM_ONE:
        {
            result = virtual_pos;

            break;
        }

        case FORM_TWO:
        {
            switch ( f_type )
            {
                case FACT_NONE:
                case FACT_INVERSE:
                {
                    result = virtual_pos;

                    break;
                }

                case FACT_CHOL:
                {
                    if ( primary_swap_flag )
                    {
                        switch ( virtual_pos )
                        {
                            case 1:
                            {
                                result = 2;

                                break;
                            }

                            case 2:
                            {
                                result = 1;

                                break;
                            }

                            default:
                            {
                                result = virtual_pos;

                                break;
                            }
                        }
                    }

                    else
                    {
                        result = virtual_pos;
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

        default:
        {
            L_THROW(0);

            break;
        }
    }

    FN_EXIT_POINT result;
}

long fMATRFact::convert_to_virtual_pos(long actual_pos)
{
    FN_ENTRY_POINT

    long result;

    switch ( f_form )
    {
        case FORM_ONE:
        {
            result = actual_pos;

            break;
        }

        case FORM_TWO:
        {
            switch ( f_type )
            {
                case FACT_NONE:
                case FACT_INVERSE:
                {
                    result = actual_pos;

                    break;
                }

                case FACT_CHOL:
                {
                    if ( primary_swap_flag )
                    {
                        switch ( actual_pos )
                        {
                            case 1:
                            {
                                result = 2;

                                break;
                            }

                            case 2:
                            {
                                result = 1;

                                break;
                            }

                            default:
                            {
                                result = actual_pos;

                                break;
                            }
                        }
                    }

                    else
                    {
                        result = actual_pos;
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

        default:
        {
            L_THROW(0);

            break;
        }
    }

    FN_EXIT_POINT result;
}

L_DOUBLE fMATRFact::get_reordered_value(const fVECTOR &what, long i)
{
    FN_ENTRY_POINT

    FN_EXIT_POINT what.get_offset_element(convert_to_virtual_pos(i)-1);
}

L_DOUBLE fMATRFact::get_reordered_value(const lVECTOR &Gindex, const fMATRIX &G_real, long i, long j)
{
    FN_ENTRY_POINT

    FN_EXIT_POINT get_deindexed_value(Gindex,G_real,convert_to_virtual_pos(i),convert_to_virtual_pos(j));
}

L_DOUBLE fMATRFact::get_deindexed_value(const lVECTOR &Gindex, const fMATRIX &G_real, long i, long j)
{
    FN_ENTRY_POINT

    L_DOUBLE result;

    if ( use_index )
    {
        result = G_real.get_offset_element((Gindex.get_offset_element(i-1))-1,(Gindex.get_offset_element(j-1))-1);
    }

    else
    {
        result = G_real.get_offset_element(i-1,j-1);
    }

    FN_EXIT_POINT result;
}

void fMATRFact::overwrite_L(const lVECTOR &Gindex, const fMATRIX &G_real)
{
    FN_ENTRY_POINT

    if ( f_type == FACT_NONE )
    {
        L_THROW(0);
    }

    long i,j,ssiizzee;

    ssiizzee = Gindex.get_effective_size();

    if ( ssiizzee > 0 )
    {
        for ( i = 1 ; i <= ssiizzee ; i++ )
        {
            for ( j = 1 ; j <= i ; j++ )
            {
                L(i-1,j-1) = get_deindexed_value(Gindex,G_real,i,j);
            }
        }
    }

    fix_symmetry();

    FN_EXIT_POINT;
}



//
// Symmetry quick fix
//

void fMATRFact::fix_symmetry(void)
{
    FN_ENTRY_POINT

    switch ( f_type )
    {
        case FACT_INVERSE:
        {
            L.make_symmetric(DO_FORCE);

            break;
        }

        case FACT_CHOL:
        {
            L.make_lower_triangular(DO_FORCE);

            break;
        }

        case FACT_NONE:
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






