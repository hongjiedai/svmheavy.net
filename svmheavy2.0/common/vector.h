
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

#ifndef _vector_h
#define _vector_h

#include <iostream>
#include <string.h>


#include "svdefs.h"
#include "search.h"
#include "c_double.h"
#include "friends.h"


// ========================================================================

class fRef_pair;
class fVECTOR;
class slow_f_fVECTOR;
class vpVECTOR;
class f_fVECTOR;
class iVECTOR;
class lVECTOR;
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

std::ostream &operator<<(std::ostream &output, const fVECTOR &source);
std::istream &operator>>(std::istream &input,        fVECTOR &dest  );

//
// Mathematical operator overloading
//

// + posation - unary, return rvalue
// - negation - unary, return rvalue

fVECTOR   operator+ (const fVECTOR  &left_op);
fVECTOR   operator- (const fVECTOR  &left_op);

// + addition       - binary, return rvalue
// - subtraction    - binary, return rvalue

fVECTOR   operator+ (const fVECTOR  &left_op, const fVECTOR  &right_op);
fVECTOR   operator+ (const fVECTOR  &left_op, const   double &right_op);
fVECTOR   operator+ (const fVECTOR  &left_op, const c_double &right_op);
fVECTOR   operator+ (const   double &left_op, const fVECTOR  &right_op);
fVECTOR   operator+ (const c_double &left_op, const fVECTOR  &right_op);

fVECTOR   operator- (const fVECTOR  &left_op, const fVECTOR  &right_op);
fVECTOR   operator- (const fVECTOR  &left_op, const   double &right_op);
fVECTOR   operator- (const fVECTOR  &left_op, const c_double &right_op);
fVECTOR   operator- (const   double &left_op, const fVECTOR  &right_op);
fVECTOR   operator- (const c_double &left_op, const fVECTOR  &right_op);

// += additive       assignment - binary, return lvalue
// -= subtractive    assignment - binary, return lvalue

fVECTOR  &operator+=(      fVECTOR  &left_op, const fVECTOR  &right_op);
fVECTOR  &operator+=(      fVECTOR  &left_op, const   double &right_op);
fVECTOR  &operator+=(      fVECTOR  &left_op, const c_double &right_op);

fVECTOR  &operator-=(      fVECTOR  &left_op, const fVECTOR  &right_op);
fVECTOR  &operator-=(      fVECTOR  &left_op, const   double &right_op);
fVECTOR  &operator-=(      fVECTOR  &left_op, const c_double &right_op);

// * multiplication - binary, return rvalue
// / division       - binary, return rvalue

fVECTOR   operator* (const fVECTOR  &left_op, const   double &right_op);
fVECTOR   operator* (const fVECTOR  &left_op, const c_double &right_op);

fVECTOR   operator* (const   double &left_op, const fVECTOR  &right_op);
fVECTOR   operator* (const c_double &left_op, const fVECTOR  &right_op);

fVECTOR   operator/ (const fVECTOR  &left_op, const   double &right_op);
fVECTOR   operator/ (const fVECTOR  &left_op, const c_double &right_op);

// * multiplication - binary, return rvalue

L_DOUBLE  operator* (const fVECTOR  &left_op, const fVECTOR  &right_op);

// *= multiplicative assignment - binary, return lvalue
// /= divisive       assignment - binary, return lvalue

fVECTOR  &operator*=(      fVECTOR  &left_op, const   double &right_op);
fVECTOR  &operator*=(      fVECTOR  &left_op, const c_double &right_op);

fVECTOR  &operator/=(      fVECTOR  &left_op, const   double &right_op);
fVECTOR  &operator/=(      fVECTOR  &left_op, const c_double &right_op);

//
// Relational operator overloading
//

// == equivalence

int operator==(const fVECTOR  &left_op, const fVECTOR  &right_op);

int operator==(const fVECTOR  &left_op, const   double &right_op);
int operator==(const fVECTOR  &left_op, const c_double &right_op);
int operator==(const   double &left_op, const fVECTOR  &right_op);
int operator==(const c_double &left_op, const fVECTOR  &right_op);

// != inequivalence

int operator!=(const fVECTOR  &left_op, const fVECTOR  &right_op);

int operator!=(const fVECTOR  &left_op, const   double &right_op);
int operator!=(const fVECTOR  &left_op, const c_double &right_op);
int operator!=(const   double &left_op, const fVECTOR  &right_op);
int operator!=(const c_double &left_op, const fVECTOR  &right_op);

// < less than

int operator< (const fVECTOR  &left_op, const fVECTOR  &right_op);

int operator< (const fVECTOR  &left_op, const   double &right_op);
int operator< (const fVECTOR  &left_op, const c_double &right_op);
int operator< (const   double &left_op, const fVECTOR  &right_op);
int operator< (const c_double &left_op, const fVECTOR  &right_op);

// <= less than or equal to

int operator<=(const fVECTOR  &left_op, const fVECTOR  &right_op);

int operator<=(const fVECTOR  &left_op, const   double &right_op);
int operator<=(const fVECTOR  &left_op, const c_double &right_op);
int operator<=(const   double &left_op, const fVECTOR  &right_op);
int operator<=(const c_double &left_op, const fVECTOR  &right_op);

// > greater than

int operator> (const fVECTOR  &left_op, const fVECTOR  &right_op);

int operator> (const fVECTOR  &left_op, const   double &right_op);
int operator> (const fVECTOR  &left_op, const c_double &right_op);
int operator> (const   double &left_op, const fVECTOR  &right_op);
int operator> (const c_double &left_op, const fVECTOR  &right_op);

// >= greater than or equal to

int operator>=(const fVECTOR  &left_op, const fVECTOR  &right_op);

int operator>=(const fVECTOR  &left_op, const   double &right_op);
int operator>=(const fVECTOR  &left_op, const c_double &right_op);
int operator>=(const   double &left_op, const fVECTOR  &right_op);
int operator>=(const c_double &left_op, const fVECTOR  &right_op);





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
// A vector consists of a 1-d array of elements of some type, as well as
// a size (for bound checking) and offset start and offset end markers.
// Hence, using overloaded operators, you can do things like safe matrix
// multiplication and stuff.
//
// For all overloaded operators (but not normal functions, like _swap), the
// vector being refered to is assumed to start at offset_start_effective and
// end at offset_end_effective - the exception is stream io, which prints
// everything.
//
// Otherwise, references are to the actual vectors without offsets.
//
// The effect of this is that mathematical operations will be effected, but
// control operations will not be.
//
// If offset_start == offset_end+1 then a vector of zero size is assumed for
// overloaded operations.
//
// Special types:
// ==============
//
// NORMAL_VECTOR: standard vector type.
// ZERO_VECTOR:   zero vector of indeterminate length.
// SCALAR_VECTOR: scalar vector of indeterminate length.
// SELECT_VECTOR: only a given element will be non-zero, indeterminate len.
//

#define NORMAL_VECTOR   0
#define ZERO_VECTOR     1
#define SCALAR_VECTOR   2
#define SELECT_VECTOR   3

//
// Related values for size.
//

#define NORMAL_VECTOR_SIZE      0
#define ZERO_VECTOR_SIZE        -1
#define SCALAR_VECTOR_SIZE      -2
#define SELECT_VECTOR_SIZE      -3


//
// Ref_pair
// ========
//
// Return structure which gives a value and it's position in the vector.
// Note that the position is relative to the offset start - ie. if the
// value is the first element in the offset vector, it will be marked as 1,
// not as offset start.
//
// This is used by the overloaded operations x(ARG_MAX) and x(ARG_MIN),
// which, respectively, return the largest and smallest elements in a vector.
//
// If the there IS no vector, then the return value will be element_num = -1
// and element_val = 0 (or possibly m_diag).
//

std::ostream &operator<<(std::ostream &output, const fRef_pair &source);

class fRef_pair
{
    friend std::ostream &operator<<(std::ostream &output, const fRef_pair &source);

    public:

    long element_num;
    L_DOUBLE element_val;

    long element_num_row;

    // for matrices, element_num gives the column number and element_num_row
    // the row number
};

enum Vect_operation {FIND_MAX,FIND_MIN};

//
// If elements are added to the vector that make it exceed the allocated
// size, then memory will be allocated to current size + LOOK_SIZE
//

#ifndef LOOK_SIZE
#define LOOK_SIZE       220
#endif




//
// Vector class
//

class fVECTOR
{
    ALL_FRIENDS

    friend class fMATRIX;
    friend class fMATRFact;

    public:

    //
    // Constructor and destructor:
    //
    // fVECTOR(_size ) - create a vector of the give type / size.
    //                   default is a ZERO_VECTOR.
    // fVECTOR(source) - copy constructor
    //
    // ~fVECTOR() - destructor
    //

    fVECTOR();
    fVECTOR(char dummy, long _size);
    fVECTOR(const fVECTOR  &source);

    ~fVECTOR();

    //
    // Negation:
    //

    void negate_it(void);

    //
    // Test and modify vector type.
    //

    void make_normal(int force = 0);
    void make_zero(int force = 0);
    void make_scalar(int force = 0);
    void make_select(int force = 0);

    int is_normal(void) const;
    int is_zero(void)   const;
    int is_scalar(void) const;
    int is_select(void) const;

    int get_v_type(void) const;

    //
    // Overloaded assignment operators.
    //
    // Notes: How this assignment works depends on what type this vector
    //        is. Hence we have four cases, differentiated by v_type:
    //
    //        NORMAL_VECTOR: in this case, the vector type will not change.
    //                       The elements themselves will be set according
    //                       to right_op, and effective sizes must align if
    //                       right_op is also of type NORMAL_VECTOR.  The
    //                       exception is when size == 0, in which case sizes
    //                       need not align, and if right_op is a
    //                       NORMAL_VECTOR with non-zero size the size of
    //                       this will be changed to the size of right_op
    //                       before assignment.  If size is zero and right_op
    //                       is !NORMAL_VECTOR then v_type will be changed.
    //        ZERO_VECTOR:    / type will be changed to whatever right_op is
    //        SCALAR_VECTOR: -| and size changed as may be necessary prior
    //        SELECT_VECTOR:  \ to assignment.
    //
    // The non-vectorial forms treat the scalar argument like a vector of
    // type SCALAR_VECTOR and assign to this accordingly.
    //

    fVECTOR &operator=(const fVECTOR  &right_op);
    fVECTOR &operator=(const   double &right_op);
    fVECTOR &operator=(const c_double &right_op);

    //
    // Set and get offsets - or reset to defaults (whole vector).  Defaults
    // are start at 1 and finish at size.  If end < start then the matrix
    // will be treated as having size 0.
    //

    void reset_offsets(void);

    void set_offset_start(long start);
    void set_offset_end(long end);

    void set_offsets(long start, long end);

    long get_offset_start(void) const;
    long get_offset_end(void)   const;

    long get_effective_size(void) const;
    long get_real_size(void)      const;

    //
    // Get standard (L_DOUBLE *) version of vector, assuming that this is
    // possible in finite dimensions.  Otherwise, throw an exception.
    //

    L_DOUBLE *operator()(void) const;

    //
    // Array element overloading
    //
    // -1: gives m_diag
    //

    L_DOUBLE &operator[](long j_elm);
    L_DOUBLE &operator()(long j_elm);

    //
    // Offset element functions:
    //
    // get_offset_element: this returns the value of the element referred to
    //                     by i+offset_start-1 (1 to size).  It returns a
    //                     *value*, NOT a reference.  -1 means m_diag.
    //                     Note: j_elm ranges from 0 to size-1, or -1.
    //

    L_DOUBLE get_offset_element(long j_elm) const;

    //
    // Misc operator function
    //

    fRef_pair operator()(char dummy, Vect_operation what) const;

    //
    // Various swap functions
    //
    // bswap: [ c ] (i-1)        [ c ] (i-1)
    //        [ d ] (j-i)  ->    [ e ] (1)
    //        [ e ] (1)          [ d ] (j-i)
    //        [ f ] (...)        [ f ] (...)
    //
    // fswap: [ c ] (i-1)        [ c ] (i-1)
    //        [ e ] (1)    ->    [ d ] (j-i)
    //        [ d ] (j-i)        [ e ] (1)
    //        [ f ] (...)        [ f ] (...)
    //
    // squareswap: [ c ] (i-1)        [ c ] (i-1)
    //             [ d ] (1)          [ f ] (1)
    //             [ e ] (j-i)  ->    [ e ] (j-i)
    //             [ f ] (1)          [ d ] (1)
    //             [ g ] (...)        [ g ] (...)
    //

    void bswap(long i, long j);
    void fswap(long i, long j);
    void squareswap(long i, long j);

    //
    // Add and remove elements
    //
    // addstart: add element to start of vector
    // addend:   add element to end of vector
    // remove:   remove element which from vector
    //
    // if no argument is given, zero is assumed.
    //

    void addstart(void);
    void addstart(const   double &what);
    void addstart(const c_double &what);
    void addstart(const fVECTOR  &what);

    void addend(void);
    void addend(const   double &what);
    void addend(const c_double &what);
    void addend(const fVECTOR  &what);

    void remove(long which);

    //
    // trim_to_size: make vector normal, with size what >= 0
    //

    void trim_to_size(long what);

    //
    // Miscellaneous:
    //
    // pad_vector: if size == 0 or type is not NORMAL_VECTOR, this will
    //             make type NORMAL_VECTOR and size = _size with all
    //             elements zero.  Otherwise, _size zeros are added to
    //             the end of the vector.
    //

    void pad_vector(long _size);

    //
    // Get and set selector element.
    //

    void set_sel_elm(long);
    long get_sel_elm(void) const;

    //
    // Return the vect pointer.
    //
    // Warning: no type checking is done, memory cannot be deleted, does
    //          not allow for offsets, should not be used.
    //

    L_DOUBLE *get_direct_ref(void) const;

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
    // Vector data:
    //
    // v_type: vector type
    //
    // size: size of vector.
    //
    // offset_start_effective: effective start.
    // offset_end_effective:   effective end.
    //
    // vect: contains actual data.
    //
    // m_diag: either scalar element or select element.
    //
    // sel_elm: selector element.
    //
    // alloc_lookahead: If addstart and addend are used, then vect may
    //                  actually have more than size T *'s allocated to it.
    //                  The actual number of T *'s allocated is size
    //                  size + alloc_lookahead.
    //

    long size;
    long offset_start_effective;
    long offset_end_effective;
    L_DOUBLE *vect;
    L_DOUBLE vector_m_diag;
    long sel_elm;
    long alloc_lookahead;
    int v_type;
};






//
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=                                                                +=+=
// +=+=                      VECTOR OF VECTORS (slow version)          +=+=
// +=+=                                                                +=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
//

//
// Stream operator overloading
//

std::ostream &operator<<(std::ostream &output, const slow_f_fVECTOR &source);
std::istream &operator>>(std::istream &input,        slow_f_fVECTOR &dest  );

//
// Vector class
//

class slow_f_fVECTOR
{
    ALL_FRIENDS

    public:

    //
    // Constructor and destructor:
    //
    // slow_f_fVECTOR(_size ) - create a vector of the give type / size.
    //                          default is a ZERO_VECTOR.
    // slow_f_fVECTOR(source) - copy constructor
    //
    // ~slow_f_fVECTOR() - destructor
    //

    slow_f_fVECTOR();
    slow_f_fVECTOR(char dummy, long _size);
    slow_f_fVECTOR(const slow_f_fVECTOR  &source);

    ~slow_f_fVECTOR();

    //
    // Test and modify vector type.
    //

    void make_normal(int force = 0);
    void make_zero(int force = 0);
    void make_scalar(int force = 0);
    void make_select(int force = 0);

    int is_normal(void) const;
    int is_zero(void)   const;
    int is_scalar(void) const;
    int is_select(void) const;

    int get_v_type(void) const;

    //
    // Overloaded assignment operators.
    //

    slow_f_fVECTOR &operator=(const slow_f_fVECTOR  &right_op);

    //
    // Set and get offsets - or reset to defaults (whole vector).  Defaults
    // are start at 1 and finish at size.  If end < start then the matrix
    // will be treated as having size 0.
    //

    void reset_offsets(void);

    void set_offset_start(long start);
    void set_offset_end(long end);

    void set_offsets(long start, long end);

    long get_offset_start(void) const;
    long get_offset_end(void)   const;

    long get_effective_size(void) const;
    long get_real_size(void)      const;

    //
    // Array element overloading
    //
    // -1: gives m_diag
    //

    fVECTOR &operator[](long j_elm);
    fVECTOR &operator()(long j_elm);

    //
    // Offset element functions:
    //

    fVECTOR get_offset_element(long j_elm) const;

    //
    // Various swap functions
    //

    void bswap(long i, long j);
    void fswap(long i, long j);
    void squareswap(long i, long j);

    //
    // Add and remove elements
    //

    void addstart(void);
    void addstart(const fVECTOR &what);
    void addstart(const slow_f_fVECTOR  &what);

    void addend(void);
    void addend(const fVECTOR &what);
    void addend(const slow_f_fVECTOR  &what);

    void remove(long which);

    //
    // trim_to_size: make vector normal, with size what >= 0
    //

    void trim_to_size(long what);

    //
    // Miscellaneous:
    //

    void pad_vector(long _size);

    //
    // Get and set selector element.
    //

    void set_sel_elm(long);
    long get_sel_elm(void) const;

    //
    // Return the vect pointer.
    //
    // Warning: no type checking is done, memory cannot be deleted, does
    //          not allow for offsets, should not be used.
    //

    fVECTOR *get_direct_ref(void) const;

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
    // Vector data:
    //
    // v_type: vector type
    //
    // size: size of vector.
    //
    // offset_start_effective: effective start.
    // offset_end_effective:   effective end.
    //
    // vect: contains actual data.
    //
    // m_diag: either scalar element or select element.
    //
    // sel_elm: selector element.
    //
    // alloc_lookahead: If addstart and addend are used, then vect may
    //                  actually have more than size T *'s allocated to it.
    //                  The actual number of T *'s allocated is size
    //                  size + alloc_lookahead.
    //

    long size;
    long offset_start_effective;
    long offset_end_effective;
    fVECTOR *vect;
    fVECTOR vector_m_diag;
    long sel_elm;
    long alloc_lookahead;
    int v_type;
};





//
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=                                                                +=+=
// +=+=                      VECTOR OF void *'s                        +=+=
// +=+=                                                                +=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
//

//
// Reason: a vector of vectors, implemented as above (slow_f_fVECTOR) is
//         (1) *extremely* slow and (2) requires that all vector elements
//         be the same size and type (for example, during fswap and bswap).
//         These problems may be overcome using a vector of pointers to
//         fVECTOR.  For simplicity and generality, I have chosen to do this
//         in two stages - firstly, implement a vector of pointers (this
//         class) that can hold anything, and then implement f_fVECTOR as
//         a wrapper to a vector of pointers.  The wrapper takes care of
//         details like memory management (making sure pointers all point
//         to something appropriate, freeing memory etc.) while the vector
//         of pointers does all the "vectorial" stuff.
//

//
// There seems to be no void * reference type in C++ - but we need one here,
// so this is it:
//

class void_point_ref
{
    public:

    void_point_ref();
    void_point_ref(const void_point_ref &source);

    ~void_point_ref();

    void_point_ref &operator=(const void_point_ref &source);

    void *void_point;
};

int operator==(const void_point_ref &left_op, const void_point_ref &right_op);
int operator!=(const void_point_ref &left_op, const void_point_ref &right_op);

void fvector_void_point_printer(std::ostream &output, const void_point_ref &what);
void fvector_void_point_scanner(std::istream &input ,       void_point_ref &what, std::ostream *where_to, int echo_level);

void throw_void_point_printer(std::ostream &output, const void_point_ref &what);
void throw_void_point_scanner(std::istream &input ,       void_point_ref &what, std::ostream *where_to, int echo_level);

//
// Stream operator overloading
//

std::ostream &operator<<(std::ostream &output, const vpVECTOR &source);
std::istream &operator>>(std::istream &input,        vpVECTOR &dest  );

//
// Vector class
//

class vpVECTOR
{
    ALL_FRIENDS

    public:

    //
    // Constructor and destructor:
    //
    // vpVECTOR(_size)  - create a vector of the give type / size.
    //                    default is a ZERO_VECTOR.
    // vpVECTOR(source) - copy contructor.
    //
    // ~vpVECTOR() - destructor
    //

    vpVECTOR();
    vpVECTOR(char dummy, long _size);

    ~vpVECTOR();

    //
    // Test and modify vector type.
    //

    void make_normal(int force = 0);
    void make_zero(int force = 0);

    int is_normal(void) const;
    int is_zero(void)   const;

    int get_v_type(void) const;

    //
    // Set and get offsets - or reset to defaults (whole vector).  Defaults
    // are start at 1 and finish at size.  If end < start then the matrix
    // will be treated as having size 0.
    //

    void reset_offsets(void);

    void set_offset_start(long start);
    void set_offset_end(long end);

    void set_offsets(long start, long end);

    long get_offset_start(void) const;
    long get_offset_end(void)   const;

    long get_effective_size(void) const;
    long get_real_size(void)      const;

    //
    // Array element overloading.
    //
    // -1: gives m_diag
    //

    void_point_ref &operator[](long j_elm);
    void_point_ref &operator()(long j_elm);

    //
    // Offset element functions
    //

    void_point_ref get_offset_element(long j_elm) const;

    //
    // Various swap functions
    //

    void bswap(long i, long j);
    void fswap(long i, long j);
    void squareswap(long i, long j);

    //
    // Add and remove elements - note that actually freeing any memory
    // associated with the pointer must be done externally.
    //

    void addstart(void);
    void addstart(const void_point_ref &what);

    void addend(void);
    void addend(const void_point_ref &what);

    void remove(long which);

    //
    // trim_to_size: make vector normal, with size what >= 0
    //

    void trim_to_size(long what);

    //
    // Miscellaneous:
    //

    void pad_vector(long _size);

    //
    // Get and set selector element.
    //

    void set_sel_elm(long);
    long get_sel_elm(void) const;

    //
    // Return the vect pointer.
    //
    // Warning: no type checking is done, memory cannot be deleted, does
    //          not allow for offsets, should not be used.
    //

    void_point_ref *get_direct_ref(void) const;

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
    void set_io_funs(void (*void_point_printer)(std::ostream &output,const void_point_ref &what), void (*void_point_scanner)(std::istream &input, void_point_ref &what, std::ostream *where_to, int echo_level));
    private:
    std::ostream *where_to;
    int echo_level;

    private:

    //
    // Vector data:
    //
    // v_type: vector type
    //
    // size: size of vector.
    //
    // offset_start_effective: effective start.
    // offset_end_effective:   effective end.
    //
    // vect: contains actual data.
    //
    // m_diag: either scalar element or select element.
    //
    // sel_elm: selector element.
    //
    // alloc_lookahead: If addstart and addend are used, then vect may
    //                  actually have more than size T *'s allocated to it.
    //                  The actual number of T *'s allocated is size
    //                  size + alloc_lookahead.
    //

    long size;
    long offset_start_effective;
    long offset_end_effective;
    void_point_ref *vect;
    void_point_ref vector_m_diag;
    long sel_elm;
    long alloc_lookahead;
    int v_type;

    //
    // Stream io helpers:
    //
    // Stream io doesn't know how to handle void (obviously).  These
    // functions are called (if set - calling a function at NULL is a
    // bad idea (tm)) to handle this stuff.
    //
    // NOTES: - void_point_printer must throw an exception for NULL pointer.
    //        - void_point_scanner must allocate for NULL pointer.
    //

    void (*void_point_printer)(std::ostream &output, const void_point_ref &what);
    void (*void_point_scanner)(std::istream &input ,       void_point_ref &what, std::ostream *where_to, int echo_level);
};





//
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=                                                                +=+=
// +=+=                      VECTOR OF VECTORS (fast but crude version)+=+=
// +=+=                                                                +=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
//

//
// Nearly functionally equivalent to the slow version, but is quicker
// (obviously), uses pointer hackery and can handle vector elements of
// different dimensions and types with significantly more applomb.
//
// NB: this class deviates from the standard macro approach.
//

//
// Stream operator overloading
//

std::ostream &operator<<(std::ostream &output, const f_fVECTOR &source);
std::istream &operator>>(std::istream &input,        f_fVECTOR &dest  );

//
// Vector class
//

class f_fVECTOR
{
    ALL_FRIENDS

    public:

    //
    // Constructor and destructor:
    //
    // f_fVECTOR(_size)  - create a vector of a given size (normal only).
    // f_fVECTOR(source) - copy constructor
    //
    // ~f_fVECTOR() - destructor
    //

    f_fVECTOR();
    f_fVECTOR(char dummy, long _size);
    f_fVECTOR(const f_fVECTOR &source);

    ~f_fVECTOR();

    //
    // Overloaded assignment operators.
    //

    f_fVECTOR &operator=(const f_fVECTOR &right_op);

    //
    // Information function
    //

    long get_effective_size(void) const;
    long get_real_size(void) const;

    //
    // Array element overloading
    //

    fVECTOR &operator[](long j_elm);
    fVECTOR &operator()(long j_elm);
    fVECTOR get_offset_element(long j_elm) const;

    //
    // Various swap functions
    //

    void bswap(long i, long j);
    void fswap(long i, long j);
    void squareswap(long i, long j);

    //
    // Add and remove elements
    //

    void addstart(const fVECTOR &what);
    void addend(const fVECTOR &what);
    void remove(long which);

    //
    // trim_to_size: make vector normal, with size what >= 0
    //

    void trim_to_size(long what);

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
    // A void pointer base takes care of things.
    //

    vpVECTOR void_base;
};





//
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=                                                                +=+=
// +=+=                      VECTOR OF ints                            +=+=
// +=+=                                                                +=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
//

//
// Stream operator overloading
//

std::ostream &operator<<(std::ostream &output, const iVECTOR &source);
std::istream &operator>>(std::istream &input,        iVECTOR &dest  );

//
// Vector class
//

class iVECTOR
{
    ALL_FRIENDS

    friend class fMATRIX;
    friend class fMATRFact;

    public:

    //
    // Constructor and destructor:
    //
    // iVECTOR(_size ) - create a vector of the give type / size.
    //                   default is a ZERO_VECTOR.
    // iVECTOR(source) - copy constructor
    //
    // ~iVECTOR() - destructor
    //

    iVECTOR();
    iVECTOR(char dummy, long _size);
    iVECTOR(const iVECTOR  &source);

    ~iVECTOR();

    //
    // Test and modify vector type.
    //

    void make_normal(int force = 0);
    void make_zero(int force = 0);
    void make_scalar(int force = 0);
    void make_select(int force = 0);

    int is_normal(void) const;
    int is_zero(void)   const;
    int is_scalar(void) const;
    int is_select(void) const;

    int get_v_type(void) const;

    //
    // Overloaded assignment operators.
    //

    iVECTOR &operator=(const iVECTOR  &right_op);

    //
    // Set and get offsets - or reset to defaults (whole vector).  Defaults
    // are start at 1 and finish at size.  If end < start then the matrix
    // will be treated as having size 0.
    //

    void reset_offsets(void);

    void set_offset_start(long start);
    void set_offset_end(long end);

    void set_offsets(long start, long end);

    long get_offset_start(void) const;
    long get_offset_end(void)   const;

    long get_effective_size(void) const;
    long get_real_size(void)      const;

    //
    // Get standard (long *) version of vector, assuming that this is
    // possible in finite dimensions.  Otherwise, throw an exception.
    //

    int *operator()(void) const;

    //
    // Array element overloading
    //
    // -1: gives m_diag
    //

    int &operator[](long j_elm);
    int &operator()(long j_elm);

    //
    // Offset element functions:
    //

    int get_offset_element(long j_elm) const;

    //
    // Various swap functions
    //

    void bswap(long i, long j);
    void fswap(long i, long j);
    void squareswap(long i, long j);

    //
    // Add and remove elements
    //

    void addstart(void);
    void addstart(const int &what);
    void addstart(const iVECTOR  &what);

    void addend(void);
    void addend(const int &what);
    void addend(const iVECTOR  &what);

    void remove(long which);

    //
    // trim_to_size: make vector normal, with size what >= 0
    //

    void trim_to_size(long what);

    //
    // Miscellaneous:
    //

    void pad_vector(long _size);

    //
    // Get and set selector element.
    //

    void set_sel_elm(long);
    long get_sel_elm(void) const;

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
    // Vector data:
    //
    // v_type: vector type
    //
    // size: size of vector.
    //
    // offset_start_effective: effective start.
    // offset_end_effective:   effective end.
    //
    // vect: contains actual data.
    //
    // m_diag: either scalar element or select element.
    //
    // sel_elm: selector element.
    //
    // alloc_lookahead: If addstart and addend are used, then vect may
    //                  actually have more than size T *'s allocated to it.
    //                  The actual number of T *'s allocated is size
    //                  size + alloc_lookahead.
    //

    long size;
    long offset_start_effective;
    long offset_end_effective;
    int *vect;
    int vector_m_diag;
    long sel_elm;
    long alloc_lookahead;
    int v_type;
};





//
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=                                                                +=+=
// +=+=                      VECTOR OF longs                           +=+=
// +=+=                                                                +=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
//

//
// Stream operator overloading
//

std::ostream &operator<<(std::ostream &output, const lVECTOR &source);
std::istream &operator>>(std::istream &input,        lVECTOR &dest  );

//
// Vector class
//

class lVECTOR
{
    ALL_FRIENDS

    friend class fMATRIX;
    friend class fMATRFact;

    public:

    //
    // Constructor and destructor:
    //
    // lVECTOR(_size ) - create a vector of the give type / size.
    //                   default is a ZERO_VECTOR.
    // lVECTOR(source) - copy constructor
    //
    // ~lVECTOR() - destructor
    //

    lVECTOR();
    lVECTOR(char dummy, long _size);
    lVECTOR(const lVECTOR  &source);

    ~lVECTOR();

    //
    // Test and modify vector type.
    //

    void make_normal(int force = 0);
    void make_zero(int force = 0);
    void make_scalar(int force = 0);
    void make_select(int force = 0);

    int is_normal(void) const;
    int is_zero(void)   const;
    int is_scalar(void) const;
    int is_select(void) const;

    int get_v_type(void) const;

    //
    // Overloaded assignment operators.
    //

    lVECTOR &operator=(const lVECTOR  &right_op);

    //
    // Set and get offsets - or reset to defaults (whole vector).  Defaults
    // are start at 1 and finish at size.  If end < start then the matrix
    // will be treated as having size 0.
    //

    void reset_offsets(void);

    void set_offset_start(long start);
    void set_offset_end(long end);

    void set_offsets(long start, long end);

    long get_offset_start(void) const;
    long get_offset_end(void)   const;

    long get_effective_size(void) const;
    long get_real_size(void)      const;

    //
    // Get standard (long *) version of vector, assuming that this is
    // possible in finite dimensions.  Otherwise, throw an exception.
    //

    long *operator()(void) const;

    //
    // Array element overloading
    //
    // -1: gives m_diag
    //

    long &operator[](long j_elm);
    long &operator()(long j_elm);

    //
    // Offset element functions:
    //

    long get_offset_element(long j_elm) const;

    //
    // Various swap functions
    //

    void bswap(long i, long j);
    void fswap(long i, long j);
    void squareswap(long i, long j);

    //
    // Add and remove elements
    //

    void addstart(void);
    void addstart(const long &what);
    void addstart(const lVECTOR  &what);

    void addend(void);
    void addend(const long &what);
    void addend(const lVECTOR  &what);

    void remove(long which);

    //
    // trim_to_size: make vector normal, with size what >= 0
    //

    void trim_to_size(long what);

    //
    // Miscellaneous:
    //

    void pad_vector(long _size);

    //
    // Get and set selector element.
    //

    void set_sel_elm(long);
    long get_sel_elm(void) const;

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
    // Vector data:
    //
    // v_type: vector type
    //
    // size: size of vector.
    //
    // offset_start_effective: effective start.
    // offset_end_effective:   effective end.
    //
    // vect: contains actual data.
    //
    // m_diag: either scalar element or select element.
    //
    // sel_elm: selector element.
    //
    // alloc_lookahead: If addstart and addend are used, then vect may
    //                  actually have more than size T *'s allocated to it.
    //                  The actual number of T *'s allocated is size
    //                  size + alloc_lookahead.
    //

    long size;
    long offset_start_effective;
    long offset_end_effective;
    long *vect;
    long vector_m_diag;
    long sel_elm;
    long alloc_lookahead;
    int v_type;
};







//
// Hackery:
//
// On occassion, we need to return references to elements of matrices or
// vectors that don't actually exist in memory and cannot be changed.  For
// example, the statement a = G(1,1), where G is a scalar type matrix, is
// a perfectly legal statement, whereas G(1,1) = a for the same matrix is
// not.  In this case, we use these functions to return references to static
// variables.  Upon the next call to these functions, the variable may be
// checked to see if it has changed, and a throw may occur if it has.
//

#define TEMP_BUF_SIZE_L_DOUBLE          10000
#define CHECK_BUF_DEPTH_L_DOUBLE        2
#define TEMP_BUF_SIZE_FVECTOR           100
#define CHECK_BUF_DEPTH_FVECTOR         2
#define TEMP_BUF_SIZE_INT               10000
#define CHECK_BUF_DEPTH_INT             2
#define TEMP_BUF_SIZE_LONG              10000
#define CHECK_BUF_DEPTH_LONG            2
#define TEMP_BUF_SIZE_VPR               10000
#define CHECK_BUF_DEPTH_VPR             2
L_DOUBLE &l_double_get_static_zero_ref(void);
L_DOUBLE &l_double_get_static_ref(L_DOUBLE what);
fVECTOR &fvector_get_static_zero_ref(void);
fVECTOR &fvector_get_static_ref(fVECTOR what);
int &int_get_static_zero_ref(void);
int &int_get_static_ref(int what);
long &long_get_static_zero_ref(void);
long &long_get_static_ref(long what);
void_point_ref &void_point_ref_get_static_zero_ref(void);
void_point_ref &void_point_ref_get_static_ref(void_point_ref what);


#endif
