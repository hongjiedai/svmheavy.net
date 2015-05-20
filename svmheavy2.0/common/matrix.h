
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

#ifndef _matrix_h
#define _matrix_h

#include <iostream>
#include <string.h>


#include "svdefs.h"
#include "search.h"
#include "levicivita.h"
#include "c_double.h"
#include "vector.h"
#include "friends.h"

// ========================================================================

class fRef_pair;
class fVECTOR;
class fMATRIX;
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

std::ostream &operator<<(std::ostream &output, const fMATRIX   &source);
std::istream &operator>>(std::istream &input, fMATRIX   &dest);

//
// Mathematical operator overloading
//

// + posation - unary, return rvalue
// - negation - unary, return rvalue

fMATRIX   operator+ (const fMATRIX  &left_op);
fMATRIX   operator- (const fMATRIX  &left_op);

// + addition       - binary, return rvalue
// - subtraction    - binary, return rvalue

fMATRIX   operator+ (const fMATRIX  &left_op, const fMATRIX  &right_op);
fMATRIX   operator- (const fMATRIX  &left_op, const fMATRIX  &right_op);

// += additive       assignment - binary, return lvalue
// -= subtractive    assignment - binary, return lvalue

fMATRIX  &operator+=(      fMATRIX  &left_op, const fMATRIX  &right_op);
fMATRIX  &operator-=(      fMATRIX  &left_op, const fMATRIX  &right_op);

// * multiplication - binary, return rvalue
// / division       - binary, return rvalue

fMATRIX   operator* (const fMATRIX  &left_op, const   double &right_op);
fMATRIX   operator* (const fMATRIX  &left_op, const c_double &right_op);

fMATRIX   operator* (const   double &left_op, const fMATRIX  &right_op);
fMATRIX   operator* (const c_double &left_op, const fMATRIX  &right_op);

fMATRIX   operator/ (const fMATRIX  &left_op, const   double &right_op);
fMATRIX   operator/ (const fMATRIX  &left_op, const c_double &right_op);

// * multiplication - binary, return rvalue

fMATRIX   operator* (const fMATRIX  &left_op, const fMATRIX  &right_op);
fVECTOR   operator* (const fMATRIX  &left_op, const fVECTOR  &right_op);
fVECTOR   operator* (const fVECTOR  &left_op, const fMATRIX  &right_op);

// *= multiplicative assignment - binary, return lvalue
// /= divisive       assignment - binary, return lvalue

fMATRIX  &operator*=(      fMATRIX  &left_op, const   double &right_op);
fMATRIX  &operator*=(      fMATRIX  &left_op, const c_double &right_op);

fMATRIX  &assign_all(      fMATRIX  &left_op, const   double &right_op);
fMATRIX  &assign_all(      fMATRIX  &left_op, const c_double &right_op);

fMATRIX  &operator/=(      fMATRIX  &left_op, const   double &right_op);
fMATRIX  &operator/=(      fMATRIX  &left_op, const c_double &right_op);

// *= multiplicative assignment - binary, return lvalue

fMATRIX  &operator*=(      fMATRIX  &left_op, const fMATRIX  &right_op);
fVECTOR  &operator*=(      fVECTOR  &left_op, const fMATRIX  &right_op);

//
// Relational operator overloading
//

int operator==(const fMATRIX  &left_op, const fMATRIX  &right_op);
int operator!=(const fMATRIX  &left_op, const fMATRIX  &right_op);
int operator< (const fMATRIX  &left_op, const fMATRIX  &right_op);
int operator<=(const fMATRIX  &left_op, const fMATRIX  &right_op);
int operator> (const fMATRIX  &left_op, const fMATRIX  &right_op);
int operator>=(const fMATRIX  &left_op, const fMATRIX  &right_op);




//
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=                                                                +=+=
// +=+=                       MATRIX CLASS                             +=+=
// +=+=                                                                +=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
//


//
// A matrix consists of a 2-d array of elements of some type, as well as
// a size (for bound checking) and also offset start and offset end markers
// for both rows and columns.
//
// For all overloaded operators (but not normal functions, like _swap), the
// matrix being refered to is assumed to start at row
// offset_start_effective_row and column offset_start_effective_col and
// end at row offset_end_effective_row and column offset_end_effective_col
// - the exception is stream io, which prints everything.
//
// Otherwise, references are to the actual matrices without offsets.
//
// The effect of this is that mathematical operations will be effected, but
// control operations will not be.
//
// If offset_start == offset_end+1 then a matrix of zero size is assumed for
// overloaded operations.
//

#define ASYMMETRIC_MATRIX            0 // generic asymmetric matrix
#define SYMMETRIC_MATRIX             1 // matrix is symmetric about diagonal
#define LOWER_TRIANGULAR_MATRIX      2 // matrix is lower triangular
#define UPPER_TRIANGULAR_MATRIX      3 // matrix is upper triangular
#define DIAGONAL_MATRIX              4 // matrix is diagonal
#define SCALAR_MATRIX                5 // matrix is the identity matrix * scalar
#define ZERO_MATRIX                  6 // matrix is zero


#define ASYMETRIC_MATRIX_SIZE        0
#define SYMETRIC_MATRIX_SIZE         0
#define LOWER_TRIANGULAR_MATRIX_SIZE 0
#define UPPER_TRIANGULAR_MATRIX_SIZE 0
#define DIAGONAL_MATRIX_SIZE         0
#define SCALAR_MATRIX_SIZE           -1
#define ZERO_MATRIX_SIZE             -2


#define E_EMPTY         1
#define E_ZERO          2
#define E_DIAG          3
#define E_LOWER         4
#define E_UPPER         5
#define E_SYMM          6
#define E_ASYMM         7
#define E_VZ            8
#define E_VS            9


#define INV_IS_LEFT     0
#define INV_IS_RIGHT    1


enum Matr_operation {ARG_DIAG_PROD,ARG_TRACE,ARG_DET,ARG_RANK_ONE,ARG_RANK_N,ARG_INVERT_LEFT,ARG_INVERT_RIGHT,ARG__MIN,ARG__MAX};


//
// Symmetric matrix storage struct - see below for details.
//

class l_xMatrix
{
    public:

    L_DOUBLE *row; // pointer to this row
    L_DOUBLE *col; // pointer to this column
    L_DOUBLE diag; // diagonal element

    l_xMatrix *next;  // pointer to next row matrix
    l_xMatrix *prev;  // pointer to previous row matrix
};








//
// Common zero test and throw if zero within tolerance funtion
//

#define NZ_ASSERT(a,zt)                                                 \
{                                                                       \
    if ( ( (a) <= (zt) ) && ( (a) >= -(zt) ) )                          \
    {                                                                   \
        throw INV_ERR;                                                  \
    }                                                                   \
}

#define POS_ASSERT(a,zt)                                                \
{                                                                       \
    if ( (a) <= (zt) )                                                  \
    {                                                                   \
        throw SQRT_ERR;                                                 \
    }                                                                   \
}



//
// Non-symmetric matrix type
// =========================
//
// m_type = ASYMMETRIC_MATRIX
//
// Each row/column is specified by an element in the xMatrix doubly linked
// list.  This contains a pointer to the leftmost edge of the row containing
// i elements (where i is the number of the row/column), a pointer to the
// topmost edge of the column containing i elements and a diagonal.  ie.
//
//             column pointers
//               \/ \/    \/
//          -> [ rcd c   .. c   ]
//  row     -> [ r   rcd .. c   ]
// pointers :  [ :   :      :   ]
//          -> [ r   r   .. rcd ]
//
// Note that the diagonal is repeated three times - once in the row, once
// in the column, and once as a diagonal element.  The diagonal element is
// the important one - this is used.  The other two may become different
// because of lazy coding.
//
// Symmetric Matrix type
// =====================
//
// m_type = SYMMETRIC_MATRIX
//
// These are the same as non-symmetric, but row and column pointers point
// to the same memory.
//
//          -> [ r"d "   .. "   ]
//  row     -> [ r"  r"d .. "   ]
// pointers :  [ :   :      :   ]
//          -> [ r   r   .. r"d ]
//
// move->row = move->col
//
// Lower Triangular Matrix type
// ============================
//
// m_type = LOWER_TRIANGULAR_MATRIX
//
// Essentially the same as symmetric, except that column pointers are NULL
// and columns are assumed to be zero above the diagonal.
//
// Upper Triangular Matrix type
// ============================
//
// m_type = UPPER_TRIANGULAR_MATRIX
//
// Essentially the same as symmetric, except that row pointers are NULL
// and rows are assumed to be zero above the diagonal.
//
// Diagonal Matrix type
// ====================
//
// m_type = DIAGONAL_MATRIX
//
// Both row and column pointers are NULL, so there are only diagonal
// elements.
//
// Scalar Matrix type
// ==================
//
// m_type = SCALAR_MATRIX
//
// Matrix does not actually store anything in this case.  It's reports its
// size as SCALAR_MATRIX_SIZE - references to diagonal elements will return m_diag
// while references to offdiagonals will report zero.  Default value for
// m_diag is one.
//
// Zero Matrix type
// ================
//
// m_type = ZERO_MATRIX
//
// Matrix does not actually store anything in this case.  It's reports its
// size as ZERO_MATRIX_SIZE - all references to the matrix return zero.
//
//
//
// Conversion between types
// ========================
//
// Various operations may force the matrix to cease to be of a particular
// symmetry and become another automatically.  Forced conversion may be
// achieved only between compatible types.  Symmetric matrices have a certain
// "persistence", namely that changes to the lower triangle are reflected by
// the upper diagonal in many cases, rather than being converted to
// asymmetric type.
//
// NB: scalar and zero matrices are of fixed type, and no conversion may be
//     applied to these (unless forced using DO_FORCE).
//

class kern_MATRIX;

class fMATRIX
{
    ALL_FRIENDS

    friend class fVECTOR;
    friend class fMATRFact;

    friend class kern_MATRIX;

    public:

    //
    // Constructor and destructor:
    //
    // fMATRIX({_size{,_type}}) - Construct a matrix of ZERO type if no
    //                            args used.  Otherwise, create a 0 matrix
    //                            of the given size and type.
    //
    // fMATRIX(source) - construct a matrix that is a the same as source
    //                   (the complete matrix, not just the effective
    //                   matrix).  If source is a vector then the matrix
    //                   so created will be diagonal and have the same size
    //                   as the EFFECTIVE vector.  Construction with a scalar
    //                   will create a SCALAR_MATRIX matrix with m_diag
    //                   source.
    //
    // ~fMATRIX() - destructor.
    //

    fMATRIX();
    fMATRIX(char dummy, long _size = ZERO_MATRIX_SIZE, int _m_type = ZERO_MATRIX);

    fMATRIX(const fMATRIX  &source);
    fMATRIX(const fVECTOR  &source);
    fMATRIX(const   double &source);
    fMATRIX(const c_double &source);

    ~fMATRIX();

    //
    // Overloaded assignment for matrices
    //

    fMATRIX &operator=(const fMATRIX  &right_op);
    fMATRIX &operator=(const fVECTOR  &right_op);
    fMATRIX &operator=(const   double &right_op);
    fMATRIX &operator=(const c_double &right_op);

    //
    // Assignment from indexed matrix.
    //
    // A.assign_from_index(B,b) means A(i-1,j-1) = B(b[i-1]-1,c[j-1]-1)
    // for all i,j in the (square) relevant submatrix A.
    //
    // Noting that the elements of B are selected in a symmetric manner from
    // B, it follows that symmetry re-assignment will be essentially the
    // same as for overloaded assignment above.
    //

    fMATRIX &assign_from_index(const fMATRIX &B, const lVECTOR &b);

    //
    // Set and get offsets - or reset to defaults (whole matrix, but does
    // not effect attributes).
    //

    void reset_offsets(void);

    void set_offset_start_at_row(long start_row);
    void set_offset_start_at_col(long start_col);

    void set_offset_end_at_row(long end_row);
    void set_offset_end_at_col(long end_col);

    void set_offset_start(long start_row, long start_col);
    void set_offset_end(long end_row, long end_col);

    void set_offsets(long start_row, long start_col, long end_row, long end_col);

    long get_offset_start_at_row(void) const;
    long get_offset_start_at_col(void) const;
    long get_offset_end_at_row(void)   const;
    long get_offset_end_at_col(void)   const;

    long get_effective_height(void) const;
    long get_effective_width(void) const;

    long get_real_size(void) const;

    //
    // Modify matrix symmetry attribuates - note that not all types can be
    // converted to just any other type.  Like the previous set of functions,
    // ignores any transposion of the matrix during operation.
    //
    // force: 0 - normal.  Throw an error is change will cause non-zero
    //            elements to go to zero.
    //        1 - force elements off symmetry to zero if necessary
    //

    void make_asymmetric(int force = 0);
    void make_symmetric(int force = 0);
    void make_lower_triangular(int force = 0);
    void make_upper_triangular(int force = 0);
    void make_diagonal(int force = 0);
    void make_scalar(int force = 0);
    void make_zero(int force = 0);

    int is_asymmetric(void)       const;
    int is_symmetric(void)        const;
    int is_lower_triangular(void) const;
    int is_upper_triangular(void) const;
    int is_diagonal(void)         const;
    int is_scalar(void)           const;
    int is_zero(void)             const;

    int get_m_type(void) const;

    //
    // Tranposition (whole matrix, not just effective matrix):
    //

    void transpose(void);

    //
    // Negation:
    //

    void negate_it(void);

    //
    // Misc info functions:
    //
    // is_not_eff_square(): returns 1 is effective matrix is not square,
    //                      0 otherwise.
    //

    int is_not_eff_square(void) const;

    //
    // Get standard (T **) version of matrix, assuming that this is
    // possible in finite dimensions.
    //

    L_DOUBLE **operator()(void) const;

    //
    // Array element overloading.  It works as follows (x is a matrix, y
    // is a vector, z is of class T):
    //
    // x[i-1] is a copy of row i of the effective matrix
    // x(i-1) is a copy of col i of the effective matrix
    //
    // x[-1] is a copy of the diagonal of the effective matrix
    // x(-1) is a copy of the diagonal of the effective matrix
    //
    // x(i-1,j-1) is a reference to element i,j of the effective matrix
    // x( -1,j-1) is a reference to element j,j of the effective matrix
    // x(i-1, -1) is a reference to m_diag
    // x( -1, -1) is a reference to m_diag
    //
    // Symmetry re-assignment: If an element of a diagonal or triangular
    // matrix which would otherwise be zero is overwritten with a non-zero
    // value then the symmetry will be re-assigned to either triangular
    // (for a diagonal matrix) or asymmetric.
    // If an element of a symmetric matrix is overwritten then symmetry
    // will not be re-assigned.  Instead, a mirror assignment will occur
    // to ensure that the symmetry is kept.
    // For zero and scalar ops, any write to something other than m_diag
    // may (or may not) result in a throw, and no change will be made to
    // the element in question.
    //
    // Note: x[i-1][j-1] = 10 WILL NOT WORK.
    //

    fVECTOR operator[](long i) const;
    fVECTOR operator()(long i) const;

    L_DOUBLE &operator()(long i, long j);

    //
    // get_offset_elm(i-1,j-1): acts like this(i-1,j-1) (with all the
    //                          offsetting and crap).
    //

    L_DOUBLE get_offset_element(long i, long j) const;

    //
    // Misc properties operator:
    //
    // ARG_DIAG_PROD - product of diagonals of effective matrix
    // ARG_TRACE     - sum of diagonals of effective matrix
    // ARG_DET       - determinant of effective matrix
    //
    // Note: these will throw an exception if the matrix is not square.
    //

    L_DOUBLE operator()(Matr_operation what) const;

    //
    // Rank update operators:
    //
    // operation = ARG_RANK_ONE - do a rank one update the effective matrix
    //
    // rank one: M(ARG_RANK_ONE, x, c) - M = M + xcx'
    //
    // notes: x = m dimensional vector
    //        c = scalar
    //

    void operator()(Matr_operation what, const fVECTOR &x, const   double &c);
    void operator()(Matr_operation what, const fVECTOR &x, const c_double &c);

    //
    // Misc operator function
    //
    // operation = ARG_MAX - get largest element
    //             ARG_MIN - get smallest element
    //

    fRef_pair operator()(Matr_operation what, char dummya, char dummyb) const;

    //
    // Inversion operators:
    //
    // Suppose we are given b (vector or matrix), and using the present
    // effective matrix (G), these functions calculate y (vector or matrix)
    // when:
    //
    // ARG_INVERT_LEFT:  G.y  = b
    // ARG_INVERT_RIGHT: y'.G = b'
    //
    // The arguments z_start and z_end define the form of b, such that:
    //
    //     [ 0  ] z_start
    // b = [ bx ] size(b)-z_start-z_end
    //     [ 0  ] z_end
    //
    // zt defines the zero tolerance of the operation.
    //

    void operator()(Matr_operation what, fVECTOR &y, const fVECTOR &b, long z_start = 0, long z_end = 0, double zt = 0) const;
    void operator()(Matr_operation what, fMATRIX &y, const fMATRIX &b, double zt = 0);

    //
    // Matrix pivoting/resizing:
    //
    // Notes: all of these operations ignore lower / upper triangular
    //        states and treat the matrix as either purely {a-}symmetric
    //        and square.
    //
    // squareswap(i,j) - swap row/column i/j
    //
    //                   Asymmetric:
    //
    //                   i-1   [ A  o  Q  t  R  ]    [ A  t  Q  o  R ] i-1
    //                   1     [ b' f  p' s  u' ]    [ d' m  k' h  r'] 1
    //                   j-i-1 [ C  g  J  q  S  ] -> [ C  q  J  g  S ] j-i-1
    //                   1     [ d' h  k' m  r' ]    [ b' s  p' f  u'] 1
    //                   N-j   [ E  i  L  n  P  ]    [ E  n  L  i  P ] N-j
    //
    //                   Symmetric:
    //
    //                   i-1   [ A  "* "* "* "* ]    [ A  "* "*  "* "* ] i-1
    //                   1     [ b' f  "* "* "* ]    [ d' m  "*  "* "* ] 1
    //                   j-i-1 [ C  g  J  "* "* ] -> [ C  k* J   "* "* ] j-i-1
    //                   1     [ d' h  k' m  "* ]    [ b' h* g'* f  "* ] 1
    //                   N-j   [ E  i  L  n  P  ]    [ E  n  L   i  P  ] N-j
    //
    //                   Triangular: Like symmetric - the upper (lower)
    //                               mirror is assumed to exist for the
    //                               duration of the operation.
    //
    //                   Diagonal:
    //
    //                   i-1   [ A             ]    [ A             ] i-1
    //                   1     [    f          ]    [    m          ] 1
    //                   j-i-1 [       J       ] -> [       J       ] j-i-1
    //                   1     [          m    ]    [          f    ] 1
    //                   N-j   [             P ]    [             P ] N-j
    //
    // bswap(i,j) - move row/column i back to position j
    //
    //              Asymmetric: j-1 [ A  K  l  M ]    [ A  l  K  M ] j-1
    //                          i-j [ B  E  n  O ] -> [ c' h  f' p'] 1
    //                          1   [ c' f' h  p']    [ B  n  E  O ] i-j
    //                          N-i [ D  G  i  J ]    [ D  i  G  J ] N-i
    //
    //              Symmetric: j-1 [ A  "* "* "* ]    [ A  "* "* "* ] j-1
    //                         i-j [ B  E  "* "* ] -> [ c' h  "* "* ] 1
    //                         1   [ c' f' h  "* ]    [ B  f* E  "* ] i-j
    //                         N-i [ D  G  i  J  ]    [ D  i  G  J  ] N-i
    //
    //              Triangular: Like symmetric - the upper (lower)
    //                          mirror is assumed to exist for the
    //                          duration of the operation.
    //
    //              Diagonal: j-1 [ A          ]    [ A          ] j-1
    //                        i-j [    E       ] -> [    h       ] 1
    //                        1   [       h    ]    [       E    ] i-j
    //                        N-i [          J ]    [          J ] N-i
    //
    // fswap(i,j) - move row/column i forward to position j
    //
    //              Asymmetric: i-1 [ A  l  M  N ]    [ A  M  l  N ] i-1
    //                          1   [ b' e  o' p'] -> [ C  H  f  Q ] j-i
    //                          j-i [ C  f  H  Q ]    [ b' o' e  p'] 1
    //                          N-j [ D  g  J  K ]    [ D  J  g  K ] N-j
    //
    //              Symmetric: i-1 [ A  "* "* "* ]    [ A  "*  "* "* ] i-1
    //                         1   [ b' e  "* "* ] -> [ C  H   "* "* ] j-i
    //                         j-i [ C  f  H  "* ]    [ b' f'* e  "* ] 1
    //                         N-j [ D  g  J  K  ]    [ D  J   g  K  ] N-j
    //
    //              Triangular: Like symmetric - the upper (lower)
    //                          mirror is assumed to exist for the
    //                          duration of the operation.
    //
    //              Diagonal: i-1 [ A          ]    [ A          ] i-1
    //                        1   [    e       ] -> [    H       ] j-i
    //                        j-i [       H    ]    [       e    ] 1
    //                        N-j [          K ]    [          K ] N-j
    //
    //
    // addend(c{,x{,y}}) - add new row/dolumn to end of matrix
    //
    //                     Asymmetric: M = [ Mold' d ]  x    = [ b ]
    //                                     [ b'    c ]         [ c ]
    //                                                  y    = [ d ]
    //                                                         [ c ]
    //                                                  diag =   c
    //
    //                     Symmetric: M = [ Mold b* ]  x    = [ b ]
    //                                    [  b'  c  ]         [ c ]
    //                                                diag =   c
    //
    //                     Lower triangular: M = [ Mold   ]  x    = [ b ]
    //                                           [  b'  c ]         [ c ]
    //                                                       diag =   c
    //
    //                     Upper triangular: M = [ Mold b ]  x    = [ b ]
    //                                           [      c ]         [ c ]
    //                                                       diag =   c
    //
    //                     Diagonal: M = [ Mold   ]  diag = c
    //                                   [      c ]
    //
    // addstart(b{,x{,y}}) - add new row/column to start of matrix
    //                       (non-optimal, uses callback to addend)
    //
    //                       Asymmetric: M = [ c  d'  ]  x    = [ c ]
    //                                       [ b Mold ]         [ b ]
    //                                                   y    = [ c ]
    //                                                          [ d ]
    //                                                   diag =   c
    //
    //                       Symmetric: M = [ c b'*  ]  x    = [ c ]
    //                                      [ b Mold ]         [ b ]
    //                                                  diag =   c
    //
    //                       Lower triangular: M = [ b      ]  x    = [ c ]
    //                                             [ c Mold ]         [ b ]
    //                                                         diag =   c
    //
    //                       Upper triangular: M = [ c  b'  ]  x    = [ c ]
    //                                             [   Mold ]         [ b ]
    //                                                         diag =   c
    //
    //                       Diagonal: M = [ c      ]  diag = c
    //                                     [   Mold ]
    //
    //
    // remove(rc)  - remove row/column rc from matrix (and submatrices)
    //
    // NOTE: - addstart and addend do not use pointer assignment.
    //       - same only add offset part of vector.
    //       - same check that b is same in x, y and b.
    //

    void squareswap(long i, long j);
    void bswap(long, long);
    void fswap(long, long);

    void addend(void);
    void addend(const   double &c);
    void addend(const   double &c, const fVECTOR &x);
    void addend(const   double &c, const fVECTOR &x, const fVECTOR &y);
    void addend(const c_double &c);
    void addend(const c_double &c, const fVECTOR &x);
    void addend(const c_double &c, const fVECTOR &x, const fVECTOR &y);

    void addstart(void);
    void addstart(const   double &c);
    void addstart(const   double &c, const fVECTOR &x);
    void addstart(const   double &c, const fVECTOR &x, const fVECTOR &y);
    void addstart(const c_double &c);
    void addstart(const c_double &c, const fVECTOR &x);
    void addstart(const c_double &c, const fVECTOR &x, const fVECTOR &y);

    void remove(long i);

    //
    // Misc:
    //
    // pad_matrix:    add some zero row/columns onto the end.
    // invert_matrix: invert the whole (not effective) matrix.
    //

    void pad_matrix(long _size);
    void invert_matrix(double zt = 0);

    //
    // Information about effective matrix.
    //
    // The form of the actual is not necessarily the same and the form of
    // the effective matrix.  The following functions provide functionality
    // similar to get_v_type et al, but for the effective matrix.
    //
    // get_e_type: get the type of the effective matrix (see below).
    // get_e_offs: get the offset of the effective matrix (see below).
    //
    // fix_symm: Fix symmetry for modifying ONLY the non-zero effective part
    //           of the matrix in a symmetric way.  Note that if v_type ==
    //           SYMMETRIC_MATRIX but the effective matrix type is E_ASYMM then
    //           this function will convert to v_type == ASYMMETRIC_MATRIX.
    //
    // In the following docs, all matrix references are effective unless
    // otherwise stated.  Zero refers to elements that are zero because of
    // the v_type of the matrix.
    //
    // Types:  E_EMPTY: matrix is empty.  v_type could be anything but
    //                  ZERO_MATRIX or SCALAR_MATRIX.  offset = 0.
    //         E_ZERO:  all elements of the matrix are zero, but the v_type
    //                  is not ZERO_MATRIX.  v_type may be DIAGONAL_MATRIX,
    //                  LOWER_TRIANGULAR_MATRIX or UPPER_TRIANGULAR_MATRIX.
    //                  offset as usual.
    //         E_DIAG:  matrix has a diagonal somewhere, with zero elements
    //                  on the upper right or lower left.  v_type ==
    //                  DIAGONAL_MATRIX.  offset is dependant on where the
    //                  diagonal hits the upper left corner.  offset = 0 if
    //                  is hits the corner, positive if it hits the left
    //                  side (offset = offset_start_col-offset_start_row)
    //                  and negative if it hits the top.
    //         E_LOWER: matrix is lower triangular such that there are
    //                  at least some zeroes in the upper right.  v_type ==
    //                  LOWER_TRIANGULAR_MATRIX and offset set by uppermost
    //                  diagonal as per E_DIAG.
    //         E_UPPER: matrix is upper triangular such that there are
    //                  at least some zeroes in the lower left.   v_type ==
    //                  UPPER_TRIANGULAR_MATRIX and offset set by lowermost
    //                  diagonal as per E_DIAG.
    //         E_SYMM:  matrix is symmetric and square (in this case we
    //                  must have v_type == SYMMETRIC_MATRIX).  offset = 0.
    //         E_ASYMM: matrix is asymmetric.  v_type could be anything but
    //                  ZERO_MATRIX or SCALAR_MATRIX.  offset as usual.
    //         E_VZ:    v_type == ZERO_MATRIX.  offset = 0.
    //         E_VS:    v_type == SCALAR_MATRIX.  offset = 0.
    //
    // NB - E_DIAG, E_LOWER and E_UPPER all include at least 1 zero element.
    //      Otherwise, they default to E_ASYMM.
    //

    int get_e_type(void) const;
    long get_e_offs(void) const;
    void fix_symm(void);

    //
    // IO stuff - stream IO can be either standard (non-informative) or
    // interactive (question/answer).  By calling this function, the output
    // stream for streamed input can be set, putting streamed input into
    // interactive mode for the next stream input only (it then reverts
    // to standard).
    //
    // where_to is set to NULL for standard.
    // echo_level can be set non-zero to make stream repeat everything.
    //

    void prime_io(std::ostream *_where_to,int _echo_level);
    private:
    std::ostream *where_to;
    int echo_level;

    private:

    //
    // Basic information:
    //
    // size:   size of matrix (zero if matrix empty)
    // m_type: type of matrix (see above)
    // m_diag: scalar element if matrix has type SCALAR_MATRIX
    //         default value is one.
    //

    long size;
    int m_type;

    L_DOUBLE m_diag;

    //
    // Offsetting information:
    //
    // offset_start_effective: row or column where the matrix starts
    // offset_end_effective:   number of rows/columns left off end
    //

    long offset_start_effective_row;
    long offset_start_effective_col;

    long offset_end_effective_row;
    long offset_end_effective_col;

    //
    // Matrix contents:
    //
    // all_it: matrix struct
    //

    l_xMatrix *all_it;
    vpVECTOR matrix_struct_lookup;

    //
    // Raw data operations:
    //
    // getelm(i,j): get copy of element i,j of actual matrix.  If symmetry
    //              means that this element does not exist then return zero.
    //              Conjugation compensation is automatic if upper
    //              triangular symmetric.
    //
    // getref(i,j): like getelm, but gets a reference to a point.  To do
    //              this, if the element does not exist it breaks symmetry
    //              (if possible - throw if not) and then returns the
    //              relevant reference to the point in question.
    //
    // setelm(x,i,j): set element i,j of actual matrix to x.  If symmetry
    //                means that this element does not exist, then symmetry
    //                will be adjusted.  EXCEPTION: if the matrix is
    //                symmetric, it will be assumed to remain so, so any
    //                assignment will be equivalent to two (symmetric)
    //                assignments.  This is transparent, as in the symmetric
    //                case both row and col point to the same memory.
    //                Conjugation compensation is automatic if upper
    //                triangular in this case.
    //
    // slack: 0 = normal (throw if a ref to an unchangeable is requested)
    //        1 = slack (don't throw - just return a reference which may
    //            or may not be checked for changing later.  Note that a
    //            throw may still occur if a change is found to have occured
    //            later, but this is not guaranteed, and if no change is made
    //            then no throws will occur).
    //

    L_DOUBLE getelm(long i, long j) const;
    L_DOUBLE &getref(long i, long j, int slack = 0);
    void setelm(const   double &x, long i, long j);
    void setelm(const c_double &x, long i, long j);

    //
    // Constructor subfunction.
    //

    void vector_constructor_subfunction(const fVECTOR &source);
};


//
// ACCELERATION MACRO: The following macro gets the pointer to element
//                     n in the all_it linked list.
//
// GET_XMATRIX_N - marginally faster, not const
// GET_MATRIX_N  - slower, const preserving
//

#define GET_XMATRIX_N(___qwerty_n) ((l_xMatrix *) ((matrix_struct_lookup[(___qwerty_n)-1]).void_point))
#define GET_MATRIX_N(___qwerty_n)  ((l_xMatrix *) ((matrix_struct_lookup.get_offset_element((___qwerty_n)-1)).void_point))



#endif
