
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

#ifndef _factor_h
#define _factor_h

#include <iostream>
#include <string.h>


#include "search.h"
#include "svdefs.h"
#include "c_double.h"
#include "vector.h"
#include "matrix.h"
#include "friends.h"

// ========================================================================

class fMATRFact;

#define DO_FORCE        1
#define DESTROY_VECT    1

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

std::ostream &operator<<(std::ostream &output, const fMATRFact &source);
std::istream &operator>>(std::istream &input, fMATRFact &dest);


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
// These deal with a specific problem.  Suppose we have a square positive
// definite matrix G (symmetric under transposition), and we define either:
//
// form 1. H = [ G ]
//
// form 2. H = [ 0 d' ]
//             [ d G  ]
//
// ( size = size(G) )
//
// we want to be able to efficiently do the operation:
//
// a = inv(R).b
//
// where R = H if H is non-singular, or:
//
//     a = inv(R).b
// and a = inv(R).x (optimised operation)
//
// if H is singular, where:
//
// form 1: H  = [ Gu  Go' ] size(Gu)
//              [ Go  Gl  ] size(Gl)
//         G  = [ Gu  Go' ] size(Gu)
//              [ Go  Gl  ] size(Gl)
//         R  = [ Gu  ] size(Gu)
//         Go = [ x'  ] 1
//              [ Gol ] size(Gl)-1
//
// form 2: H  = [ 0   du' dl' ] 1
//              [ du  Gu  Go' ] size(Gu)
//              [ dl  Go  Gl  ] size(Gl)
//         G  = [ Gu  Go' ] size(Gu)
//              [ Go  Gl  ] size(Gl)
//         R  = [ 0   du' ] size(Gu)
//              [ du  Gu  ] 1
//         Go = [ xl' ] 1
//              [ Gol ] size(Gl)-1
//         dl = [ dm  ] 1
//              [ dll ] size(Gl)-1
//         x  = [ dm  ] 1
//              [ xl  ] size(Gu)
//
// and Gu is as large as possible.
//


//
// =====================
// Form of factorisation
// =====================
//
// exact:
//
// FACT_INVERSE: store as an inverse factorisation
// FACT_CHOL:    store as a cholesky factorisation
// FACT_NONE:    no factorisation kept, all operations done at run-time
//
// Storage: either L = some diagonal of appropriate size (form 1) or
//          a diagonal with the diag(L) = d (form 2).
//

#define FACT_INVERSE    1
#define FACT_CHOL       2
#define FACT_NONE       3

#define FORM_ONE        1
#define FORM_TWO        2

#define FACT_START_INDIC        -10

#define DEFAULT_FACT    FACT_CHOL


//
// =============
// Indexing of G
// =============
//
// If desired, the matrix G may be indexed.  That is, when accessing element
// G(i-1,j-1) of matrix G, instead element G(Gindex[i-1]-1,Gindex[j-1]-1) is
// accessed.  This allows for virtual pivoting - that is, G is a large
// matrix that is never shuffled, but when accessed in this way it looks
// like it has been shuffled.
//
// This is controlled by use_index.  If this is 0, no indexing is used.
// Otherwise, indexing will be used.  The only functions so affected are the
// function get_indexed_value, directly or indirectly (most commonly through
// overwrite_L).
//


class fMATRFact
{
    ALL_FRIENDS

    friend class fVECTOR;
    friend class fMATRIX;

    public:

    //
    // Constuctors:
    //
    // G      - matrix to be factorised.
    // zt     - zero tolerance.
    // d      - d vector, if used.
    // f_type - factorisation type, either FACT_INVERSE, FACT_CHOL or
    //          FACT_NONE.
    //
    // f_form = FORM_ONE: does not include d argument.
    // f_form = FORM_TWO: does include d argument.
    //

    fMATRFact();

    fMATRFact(int use_index, const lVECTOR &Gindex, const fMATRIX &G,                   double _zt, int _f_type);
    fMATRFact(int use_index, const lVECTOR &Gindex, const fMATRIX &G, const fVECTOR &d, double _zt, int _f_type);
    fMATRFact(const fMATRFact &source);

    //
    // Overwrite assignment operator
    //

    fMATRFact &operator=(const fMATRFact &source);

    //
    // Rank-1 updates:
    //
    // Perform a rank-one update on G and hence fix the factorisation.
    //
    // X := X + aba^
    //
    //     [ 0  ] z_start
    // a = [ ax ] get_size()-z_start-z_end
    //     [ 0  ] z_end
    //
    // ^ has been used to indicate transposition.
    //
    //     [ G_11 G_21^ G_31^ ] z_start
    // G = [ G_21 X     G_32^ ] size(G)-z_end-z_start
    //     [ G_31 G_32  G_33  ] z_end
    //
    // NOTE: 1. for b<0, nbad may increase, decrease or stay constant
    //       2. for b>0, nbad may only decrease or stay constant
    //
    // pert_one(i,s) is for the special case:
    //
    // b = s
    // a = [ 0 ..... 0 1 0 ..... 0 ]'
    //                 ^
    //            ith element
    //
    // equivalently, z_start = i-1, z_end = get_size()-i, a = 1, b = s
    //
    // Return value: all return nbad, or -1 if nbad not known (FACT_NONE)
    //
    // NB. - G_real and d_real must be presented in the form they take
    //       BEFORE any adjustments are made.
    //

    long rankone(const lVECTOR &Gindex, const fVECTOR &a, const L_DOUBLE &b, const fMATRIX &G_real,                        long z_start = 0, long z_end = 0);
    long rankone(const lVECTOR &Gindex, const fVECTOR &a, const L_DOUBLE &b, const fMATRIX &G_real, const fVECTOR &d_real, long z_start = 0, long z_end = 0);

    long pert_one(const lVECTOR &Gindex, long i, const L_DOUBLE &b, const fMATRIX &G_real                       );
    long pert_one(const lVECTOR &Gindex, long i, const L_DOUBLE &b, const fMATRIX &G_real, const fVECTOR &d_real);

    //
    // Inversion:
    //
    //     Compute a where R.a = b, where a and b have the form:
    //
    //         [ {by} ] {1 - present if f_form = FORM_TWO}
    //     b = [ 0    ] z_start
    //         [ bz   ] size(b)-z_start-z_end-{1} = size(bx)-z_start-z_end
    //         [ 0    ] z_end
    //
    //          [ 0  ] z_start
    //     bx = [ bz ] size(bx)-z_start-z_end - is what is given
    //          [ 0  ] z_end
    //
    //     a = [ {ay} ] {1 - present if f_form = FORM_TWO}
    //         [ ax   ] size(b)-{1} = size(bx)
    //
    //     ax = [ ax ] size(bx) - is what is written to
    //
    //     thus enabling certain optimisations to occur in the calculation.
    //     Both z_start and z_end arguments are optional.  The elm_one_zero
    //     argument is used with FORM_TWO to indicate that by is zero if
    //     elm_one_zero == 1, or non-zero if elm_one_zero == 0.
    //
    // Optimal inversion:
    //
    //     This calculates a, where R.a = x   (ie. it returns s)
    //
    // Return value: all return nbad, or -1 if nbad not known (FACT_NONE)
    //

    long minverse(const lVECTOR &Gindex, fVECTOR &ax,               const fVECTOR &bx,                     const fMATRIX &G_real,                        long z_start = 0, long z_end = 0                      );
    long minverse(const lVECTOR &Gindex, fVECTOR &ax, L_DOUBLE &ay, const fVECTOR &bx, const L_DOUBLE &by, const fMATRIX &G_real, const fVECTOR &d_real, long z_start = 0, long z_end = 0, int elm_one_zero = 0);

    long near_invert(const lVECTOR &Gindex, fVECTOR &ax,               const fMATRIX &G_real                       );
    long near_invert(const lVECTOR &Gindex, fVECTOR &ax, L_DOUBLE &ay, const fMATRIX &G_real, const fVECTOR &d_real);

    //
    // Matrix manipulations:
    //
    // fact_addend: add row/column to end of G.  G_real is assumed to be
    //              the non-factorised G with the relevant row/col added
    //              prior to calling this function.  Likewise for d_real.
    //
    // fact_shrink: remove row/column i from G (and d).  G_real is assumed
    //              to be the non-factorised G without the relevant row/col
    //              removed as yet.  Likewise for d_real.
    //
    // Return value: all return nbad, or -1 if nbad not known (FACT_NONE)
    //
    // NB: - all must make appropriate adjustments to stored offsets before
    //       calling end_normal_fn().
    //

    long fact_addend(const lVECTOR &Gindex, const fVECTOR &ax,                     const fMATRIX &G_real                       );
    long fact_addend(const lVECTOR &Gindex, const fVECTOR &ax, const L_DOUBLE &dx, const fMATRIX &G_real, const fVECTOR &d_real);

    long fact_shrink(const lVECTOR &Gindex, long i, const fMATRIX &G_real                       );
    long fact_shrink(const lVECTOR &Gindex, long i, const fMATRIX &G_real, const fVECTOR &d_real);

    //
    // Information functions
    //
    // get_size() - get size - size of G
    // get_nbad() - get nbad - 0 if non-singular, >= 1 if singular
    //
    // get_case_four_flag()    - get case_four_flag
    // get_primary_swap_flag() - get primary_swap_flag
    //
    // get_f_type() - get f_type (factorisation type)
    // get_f_form() - get f_form (is d included)
    //
    // test_fact() - calculates G_real based on what is present in the
    //               current factorisation.
    //

    long get_size(void);
    long get_nbad(void);

    int get_case_four_flag(void);
    int get_primary_swap_flag(void);

    int get_f_type(void);
    int get_f_form(void);

    fMATRIX test_fact(void);  // Note - ordering will be indexed

    private:

    //
    // nbad - size(Gl)
    //
    // case_four_flag    - described below.  Normally 0, but 1 indicates a
    //                     special case for f_type = FACT_CHOL
    // primary_swap_flag - described below.  Only used if f_type = FACT_CHOL
    //
    // s       - inverse factorisation only. s = inv(R).x - this is
    //           calculated during the inverse update part of the algorithm,
    //           but may or may not exist at any time (non-existance
    //           indicated by s_valid flag).  It's ordering is internal
    //           (L order).  Only relevant for inverse factorisation.
    // s_valid - 0 is normal, 1 indicates that s is up-to-date (exists).
    //
    // L - both factorisations.
    //
    // zt - zero tolerance.  Positive is taken to mean >= zt, negative
    //      means <= -zt and zero means < zt and > -zt.
    //
    // f_form - form of H, discussed previously.
    // f_type - type of factorisation, discussed previously.
    //

    long nbad;

    int case_four_flag;
    int primary_swap_flag;

    int s_valid;
    fVECTOR s;

    fMATRIX L;

    double zt;

    int f_form;
    int f_type;

    int use_index;

    //
    // If FACT_NONE then the following is used to store the size of the
    // (apparent) factorisation.
    //

    long fact_none_size;

    //
    // Internal functions:
    //
    // xfact:    attempt to extend the factorisation as far as possible.
    //
    // xrankone: internal version of the rank-1 update function.  This
    //           works much the same as the previous version, except that
    //           it updates H (not G.  H ordering as given below) and has an
    //           additional argument, namely:
    //
    //           - hold_off: 0 - normal operation
    //                       1 - only update the factorised portion of H
    //
    //           form 1: H = H
    //
    //                   size(a) == size(H)
    //
    //                       [ 0  ] z_start
    //                   a = [ ax ] size(G)-z_start-z_end
    //                       [ 0  ] z_end
    //
    //           form 2: H = H
    //
    //                   size(a) == size(H)
    //
    //                       [ 0  ] 1
    //                       [ 0  ] 1
    //                   a = [ 0  ] z_start-1              z_start >= 1
    //                       [ ax ] size(G)-z_start-z_end
    //                       [ 0  ] z_end
    //
    //                       [ ax1 ] 1
    //                       [ 0   ] 1
    //                   a = [ ax2 ] size(G)-z_end         z_start = 0
    //                       [ 0   ] z_end
    //
    // refresh_fact: Assuming that L has just been overwritten with G
    //               (not H) and likewise d_store with d, without reference
    //               to primary_swap_flag or case_four_flag, reset all
    //               relevant stuff to start factorisation from scratch
    //               (DOES NOT call xfact).
    //

    long xfact(const lVECTOR &Gindex, const fMATRIX &G_real                       );
    long xfact(const lVECTOR &Gindex, const fMATRIX &G_real, const fVECTOR &d_real);
    long xrankone(const lVECTOR &Gindex, const fVECTOR &a, const L_DOUBLE &b, const fMATRIX &G_real,                        int hold_off = 0, long z_start = 0, long z_end = 0);
    long xrankone(const lVECTOR &Gindex, const fVECTOR &a, const L_DOUBLE &b, const fMATRIX &G_real, const fVECTOR &d_real, int hold_off = 0, long z_start = 0, long z_end = 0);
    long refresh_fact(void);
    long refresh_fact(const fVECTOR &d_real);

    //
    // Reality <-> illusion conversion functions:
    //
    // enter_fixvect_normal_fn: if f_form == FORM_TWO then do the following,
    //                          in the order given:
    //                          a - if primary_swap_flag == 1 then swap
    //                              elements 2 and 3 in the vector.
    //                          b - if f_type = FACT_CHOL then swap elements
    //                              1 and 2 in the vector (even if get_size()
    //                              returns 0 - this behaviour is important
    //                              for a special case of addend).
    // exit_fixvect_normal_fn:  if f_form == FORM_TWO then do the following,
    //                          in the order given:
    //                          a - if f_type = FACT_CHOL then swap elements
    //                              2 and 1 in the vector.
    //                          b - if primary_swap_flag == 1 then swap
    //                              elements 3 and 2 in the vector.
    //                          exception: if get_size() == 0 then DO NOTHING
    //
    // convert_to_actual_pos:  convert position in visible G to actual
    //                         (stored as a part of H) G (but not actually
    //                         wrt H).
    // convert_to_virtual_pos: the above in reverse.
    //
    // get_reordered_value: when a variable is constant, it can't be
    //                      re-ordered, so to get an element (indexed from
    //                      1 {,1} to n {,n}, use this function.
    //                      get_reordered_value(v,i)   gives v[imod-1]
    //                      get_reordered_value(G,i,j) gives G[imod-1][jmod-1]
    //
    // get_deindexed_value: get element G_real(i-1,j-1), allowing for the
    //                      (potential) use of the index Gindex.  That is,
    //                      if indexing is off, just return G_real(i-1,j-1),
    //                      otherwise G_real(Gindex[i-1]-1,Gindex[j-1]-1).
    // overwrite_L:         overwrites the L matrix with G_real.

    long enter_fixvect_normal_fn(fVECTOR &a);
    long exit_fixvect_normal_fn(fVECTOR &a);

    long convert_to_actual_pos(long virtual_pos);
    long convert_to_virtual_pos(long actual_pos);

    L_DOUBLE get_reordered_value(const fVECTOR &d_real, long i);
    L_DOUBLE get_reordered_value(const lVECTOR &Gindex, const fMATRIX &G_real, long i, long j);

    L_DOUBLE get_deindexed_value(const lVECTOR &Gindex, const fMATRIX &G_real, long i, long j);
    void overwrite_L(const lVECTOR &Gindex, const fMATRIX &G_real);

    //
    // Symmetry quick fix
    //

    void fix_symmetry(void);
};




//
//                            =======================
//                            Factorisation specifics
//                            =======================
//
// =======================
// Cholesky factorisations
// =======================
//
// Note: Q'.q = s
//
// =============================
// Form 1 Cholesky factorisation
// =============================
//
// case 1
// ======
//
// If size == 0 then no H will exist.
//
// Stored as:
//
// L = empty matrix
//
// Relevant part
//
// Q is undefined.
// R is undefined.
// q is undefined.
// r is undefined.
//
// nbad == 0
//
// case 2
// ======
//
// Any positive definite symmetric matrix G may be written as follows:
//
// H = [ G_a ]
//   = [ L_a ].[ L_a ]^
//
// where: upper_left(G_a) >= zt*zt
//        diag(L_a) >= zt
//        L_a is lower triangular.
//
// Stored as:
//
// L = [ L_a ] size
//
// Relevant part:
//
// Q = [ L_a ] size
// R = [ G_a ] size
// q is undefined.
// r is undefined.
//
// nbad = 0
//
// case 3
// ======
//
// If G is positive semidefinite symmetric, then we may always write it as
// follows:
//
//     [ [ G_a ]              [ g_d'^ G_e^ ]     ] size-nbad
// H = [                                         ]
//     [ [ g_d' ]             [ g_m   g_n^ ]     ] 1
//     [ [ G_e  ]             [ g_n   G_o  ]     ] nbad-1
//
//     [ [ L_a  ].[ L_a ]^    [ L_a ].[ l_d' ]^ ]
//     [                              [ L_e  ]  ]
//   = [                                        ]
//     [ [ l_d' ].[ L_a ]^    [ g_m   g_n^ ]    ]
//     [ [ L_e  ]             [ g_n   G_o  ]    ]
//
// where: upper_left(G_a) >= zt*zt  (unless nbad == size, ie G_a is 0 size)
//        diag(L_a) >= zt           (unless nbad == size, ie L_a is 0 size)
//        L_a is lower triangular   (unless nbad == size, ie L_a is 0 size)
//        (g_m - l_d'.l_d'^) < zt*zt
//
// Stored as:
//
//     [ [ L_a  ]                       ] size-nbad
// L = [                                ]
//     [ [ l_d' ]    [ g_m            ] ] 1
//     [ [ L_e  ]    [ g_n lower(G_o) ] ] nbad-1
//
// Relevant part:
//
// Q = [ L_a ] size-nbad
// R = [ G_a ] size-nbad
// q = [ l_d ] size-nbad
// r = [ g_d ] size-nbad
//
// 1 <= nbad <= size
//
// case 4
// ======
//
// nbad = FACT_START_INDIC
//
// This indicates that factorisation has not been started yet.  It is a
// temporary state that only occurs between the constructor start and
// when the constructor calls xfact.
//
// =============================
// Form 2 Cholesky factorisation
// =============================
//
// If G is positive semidefinite symmetric and H is defined by
//
// H = [ 0  d^ ] 1
//     [ d  G  ] size
//
// where d has no zero elements.
//
// case 1
// ======
//
// If size = 0 then H will be as follows:
//
// H = [ 0 ] 1
//
// Stored as:
//
// L = [ 0 ] 1
//
// Relevant part:
//
// Q is undefined.
// R is undefined.
// q is undefined.
// r is undefined.
//
// nbad == 0
//
// case 2
// ======
//
// If G is positive definite symmetric then we may write H as follows:
//
//     [ l_a  0    0    ].[ 1 0  0  ].[ l_a  0    0    ]^ 1
// H = [ l_b  l_f  0'   ] [ 0 -1 0' ] [ l_b  l_f  0'   ]  1
//     [ l_c  l_g  L_j  ] [ 0 0  I  ] [ l_c  l_g  L_j  ]  size-1
//
//     [ g_a  d_b* g_c^ ] 1
// H = [ d_b  0    d_g^ ] 1
//     [ g_c  d_g  G_j  ] size-1
//
// where: g_a >= zt*zt
//        l_a = sqrt(g_a) >= zt
//        l_b = d_b*calc_inv(l_a)
//        l_f = l_b*l_b
//        diag(L_j) >= zt         (unless size == 1, ie. L_j is size 0)
//        L_j is lower triangular (unless size == 1, ie. L_j is size 0)
//
// special case: if l_f < zt then we must ignore this fact and continue
//               regardless.  With luck, this will not cause divide by
//               zero errors, but it is a good reason for making sure that
//               the matrix G is appropriately scaled (to approx 1 on all
//               diagonal elements).
//
// Stored as:
//
//     [ l_a            ] 1
// L = [ l_b  l_f       ] 1
//     [ l_c  l_g  L_j  ] size-1
//
// Relevant part:
//
//     [ l_a            ] 1
// Q = [ l_b  l_f       ] 1
//     [ l_c  l_g  L_j  ] size-1
//
//     [ g_a  d_b* g_c^ ] 1
// R = [ d_b  0    d_g^ ] 1
//     [ g_c  d_g  G_j  ] size-1
//
// q is undefined.
// r is undefined.
//
// nbad == 0
//
// case 3
// ======
//
// If G is positive semidefinite symmetric and sqrt(g_a) >= zt then we may
// write H as follows:
//
//     [ [ l_a  0    0    ].[ 1 0  0  ].[ l_a  0    0    ]^                                      ]
//     [ [ l_b  l_f  0'   ] [ 0 -1 0' ] [ l_b  l_f  0'   ]      "lower left transpose"           ]
//     [ [ l_c  l_g  L_j  ] [ 0 0  I  ] [ l_c  l_g  L_j  ]                                       ]
// H = [                                                                                         ]
//     [ [ l_d  l_h  l_k' ].[ 1 0  0  ].[ l_a  0    0    ]^     [ g_m  g_n^ ]                    ]
//     [ [ l_e  l_i  L_l  ] [ 0 -1 0' ] [ l_b  l_f  0'   ]      [ g_n  G_o  ]                    ]
//     [                    [ 0 0  I  ] [ l_c  l_g  L_j  ]                                       ]
//
//     [ [ g_a  d_b* g_c^ ]      [ g_d*   g_e^ ] ] 1
//     [ [ d_b  0    d_g^ ]      [ d_h*   d_i^ ] ] 1
// H = [ [ g_c  d_g  G_j  ]      [ g_k'^  G_l^ ] ] size-1-nbad
//     [                                       ]
//     [ [ g_d  d_h  g_k' ]      [ g_m    g_n^ ] ] 1
//     [ [ g_e  d_i  G_l  ]      [ g_n    G_o  ] ] nbad-1
//
// where: g_a >= zt*zt
//        l_a = sqrt(g_a) >= zt
//        l_b = d_b*calc_inv(l_a)
//        l_f = l_b*l_b
//        diag(L_j) >= zt         (unless nbad == size-1, ie. L_j is size 0)
//        L_j is lower triangular (unless nbad == size-1, ie. L_j is size 0)
//        (g_m - l_d*l_d + l_h*l_h - l_k'.l_k'^) < zt*zt
//
// special case: if l_f < zt then we must ignore this fact and continue
//               regardless.  With luck, this will not cause divide by
//               zero errors, but it is a good reason for making sure that
//               the matrix G is appropriately scaled (to approx 1 on all
//               diagonal elements).
//
// Stored as:
//
//     [ [ l_a            ]                       ] 1
//     [ [ l_b  l_f       ]                       ] 1
// L = [ [ l_c  l_g  L_j  ]                       ] size-1-nbad
//     [                                          ]
//     [ [ l_d  l_h  l_k' ]      [ g_m          ] ] 1
//     [ [ l_e  l_i  L_l  ]      [ g_n  lt(G_o) ] ] nbad-1
//
// Relevant part:
//
//     [ [ l_a            ] ] 1
// Q = [ [ l_b  l_f       ] ] 1
//     [ [ l_c  l_g  L_j  ] ] size-1-nbad
//
//     [ [ g_a  d_b* g_c^ ] ] 1
// R = [ [ d_b  0    d_g^ ] ] 1
//     [ [ g_c  d_g  G_j  ] ] size-1-nbad
//
//     [ l_d ] 1
// q = [ l_h ] 1
//     [ l_k ] size-1-nbad
//
//     [ g_d ] 1
// r = [ d_h ] 1
//     [ g_k ] size-1-nbad
//
// 1 <= nbad <= size-1
//
// case 4 (case_four_flag)
// =======================
//
// G is positive semidefinite symmetric and has the form:
//
//     [ g_a       g_d*  g_e^ ] 1
// G = [                      ]
//     [ g_d       g_m   g_n^ ] 1
//     [ g_e       g_n   G_o  ] size-2 (nbad == size-1)
//
//     [ [ g_a  d_b* ]       [ g_d*  g_e^ ] ] 1
//     [ [ d_b  0    ]       [ d_h*  d_i^ ] ] 1
// H = [                                    ]
//     [ [ g_d  d_h  ]       [ g_m   g_n^ ] ] 1
//     [ [ g_e  d_i  ]       [ g_n   G_o  ] ] size-2 (nbad == size-1)
//
// where: g_a < zt*zt
//        g_m < zt*zt ( or size == 1 )
//
// consequently: g_d <= g_a < zt*zt
//               g_e <= g_a < zt*zt
//               g_d <= g_m < zt*zt
//               g_n <= g_m < zt*zt
//
// We cannot write this as a cholesky factorisation.  However, we can note
// that relevant part is:
//
// R = [ [ g_a  d_b* ] ] 1
//     [ [ d_b  0    ] ] 1
//
// Furthermore, as d_b = +-1, we can readily see that:
//
// inv(R) = [ 0   d_b* ] 1
//          [ d_b -g_a ] 1
//
// r = [ g_d ] 1
//     [ d_h ] 1
//
// This case will be indicated by case_four_flag being set to 1.
//
// nbad = size-1
//
// case 5
// ======
//
// nbad = FACT_START_INDIC
//
// This indicates that factorisation has not been started yet.  It is a
// temporary state that only occurs between the constructor start and
// when the constructor calls xfact.
//
// primary_swap_flag
// =================
//
// During operation, we may occasionally come accross the case:
//
//     [ g_a       g_d*  g_e^ ] 1
// G = [                      ]
//     [ g_d       g_m   g_n^ ] 1
//     [ g_e       g_n   G_o  ] size-2
//
// where: g_a <  zt*zt
//        g_m >= zt*zt
//
// In this case, we may have a significant sized non-singular block, but
// because g_a is zero we cannot factorise it as such.  The workaround is
// to use the primary_swap_flag, which is defined as follows:
//
// primary_swap_flag = 0 - normal
//                   = 1 - row/column 1 and row/column 2 of all stored
//                         G ordered (1 and 3 of non G-ordered) variables
//                         have been switched.  Hence the fixvect
//                         functions must correct for this.
//
// So, if we arrive at the above situation at some point, we call the
// pri_swi() function which will invert the state of the primary_swap_flag
// and hence the ordering of all internal data.  The effect on G will be to
// convert it as follows:
//
//     [ g_a       g_d*  g_e^ ]           [ g_m       g_d*  g_n^ ]
// G = [                      ]  :->  G = [                      ]
//     [ g_d       g_m   g_n^ ]           [ g_d       g_a   g_e^ ]
//     [ g_e       g_n   G_o  ]           [ g_n       g_e   G_o  ]
//
// which circumvents the above difficulty.
//
//
//
//
//
//
//
//
//
//
//
//
// ======================
// Inverse Factorisations
// ======================
//
// Note: s = Q.r
//
// ============================
// Form 1 Inverse Factorisation
// ============================
//
// case 1
// ======
//
// If size == 0 then no G will exist.
//
// Stored as:
//
// L = empty matrix
//
// Relevant part
//
// R is undefined.
// r is undefined.
//
// nbad == 0
//
// case 2
// ======
//
// Any positive definite symmetric matrix G may be written as follows:
//
// H = [ G_a ]  size
//   = inv(L_a) size
//
// where: L_a is symmetric.
//
// Stored as:
//
// L = [ L_a ] size
//
// Relevant part:
//
// Q = [ L_a ] size
// R = [ G_a ] size
// r is undefined.
//
// nbad == 0
//
// case 3
// ======
//
// If G is positive semidefinite symmetric, then we may always write it as
// follows:
//
//     [ [ G_a  ]    [ g_d*  G_e^ ] ] size-nbad
// H = [                            ]
//     [ [ g_d' ]    [ g_m   g_n^ ] ] 1
//     [ [ G_e  ]    [ g_n   G_o  ] ] nbad-1
//
//     [ inv(L_a)    [ g_d*  G_e^ ] ] size-nbad
//   = [                            ]
//     [ [ g_d' ]    [ g_m   g_n^ ] ] 1
//     [ [ G_e  ]    [ g_n   G_o  ] ] nbad-1
//
// where: L_a is symmetric     (unless size == nbad, ie. L_a is size 0)
//        h_d'.s'^ - zt < h_m < h_d'.s'^ + zt.
//
// Stored as:
//
//     [ [ L_a  ]    [ g_d*  G_e^ ] ] size-nbad
// L = [                            ]
//     [ [ g_d' ]    [ g_m   g_n^ ] ] 1
//     [ [ G_e  ]    [ g_n   G_o  ] ] nbad-1
//
// Relevant part:
//
// Q = [ L_a ] size-nbad
// R = [ G_a ] size-nbad
// r = [ g_d ] size-nbad
//
// 1 <= nbad <= size
//
// case 4
// ======
//
// nbad = FACT_START_INDIC
//
// This indicates that factorisation has not been started yet.  It is a
// temporary state that only occurs between the constructor start and
// when the constructor calls xfact.
//
// ============================
// Form 2 Inverse Factorisation
// ============================
//
// If G is positive semidefinite symmetric and H is defined by
//
// H = [ 0  d^ ] 1
//     [ d  G  ] size
//
// where d has no zero elements.
//
// case 1
// ======
//
// If size = 0 then H will be as follows:
//
// H = [ 0 ] 1
//
// Stored as:
//
// L = [ 0 ] 1
//
// Relevant part:
//
// Q is undefined.
// R is undefined.
// r is undefined.
//
// nbad == 0
//
// case 2
// ======
//
// If G is positive definite symmetric then we may write H as follows:
//
// H = inv([ l_f  l_b^ ]) 1
//         [ l_b  L_a  ]  size
//
//   = [ 0    d_b^ ] 1
//     [ d_b  G_a  ] size
//
// where: L_a is symmetric
//
// Stored as:
//
// L = [ l_f  l_b^ ] 1
//     [ l_b  L_a  ] size
//
// Relevant part:
//
// Q = [ l_f  l_b' ] 1
//     [ l_b  L_a  ] size
//
// R = [ 0    d_b' ] 1
//     [ d_b  G_a  ] size
//
// r is undefined.
//
// nbad == 0
//
// case 3
// ======
//
// If G is positive semidefinite symmetric then we may always write H as
// follows:
//
//     [ inv([ l_f  l_b^ ])    [ d_h*  d_i^ ] ] 1
//     [     [ l_b  L_a  ]     [ g_d'^ G_e^ ] ] size-nbad
// H = [                                      ]
//     [ [ d_h  g_d' ]         [ g_m   g_n^ ] ] 1
//     [ [ d_i  G_e  ]         [ g_n   G_o  ] ] nbad-1
//
//     [ [ 0    d_b' ]         [ d_h*   d_i^ ] ] 1
//     [ [ d_b  G_a  ]         [ g_d'^  G_e^ ] ] size-nbad
//   = [                                       ]
//     [ [ d_h  g_d' ]         [ g_m    g_n^ ] ] 1
//     [ [ d_i  G_e  ]         [ g_n    G_o  ] ] nbad-1
//
// where: L_a is symmetric
//        [ d_h g_d' ].s'^ - zt < g_m < [ d_h g_d' ].s'^ + zt.
//
// Stored as:
//
//     [ [ l_f  l_b' ]    [ d_h*  d_i^ ] ] 1
//     [ [ l_b  L_a  ]    [ g_d'^ G_e^ ] ] size-nbad
// L = [                                ]
//     [ [ d_h  g_d' ]    [ g_m   g_n^ ] ] 1
//     [ [ d_i  G_e  ]    [ g_n   G_o  ] ] nbad-1
//
// Relevant part:
//
// Q = [ l_f  l_b^ ]
//     [ l_b  L_a  ]
//
// R = [ 0    d_b^ ]
//     [ d_b  G_a  ]
//
// r = [ d_h ]
//     [ g_d ]
//
// 1 <= nbad <= size-1
//
// case 4
// ======
//
// nbad = FACT_START_INDIC
//
// This indicates that factorisation has not been started yet.  It is a
// temporary state that only occurs between the constructor start and
// when the constructor calls xfact.
//



#endif
