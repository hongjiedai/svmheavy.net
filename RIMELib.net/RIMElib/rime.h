
/*
 *  RIMElib: RuntIme Mathematical Equation Library (C++ wrapper)
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



/*********************************************************************

See erime.h for details.


RuntIme Mathematical Equation Library (C++ wrapper)
===================================================

Written by: Alistair Shilton
            Melbourne University
            apsh@ecr.mu.oz.au

Date of most recent changes: 13/3/2008

The following is intended for parsing, evaluating, differentiating and
integrating arbitrary (text based) equations.  See erime.h for details.


Equation string format
======================

See erime.h


Local Variables, Placeholders and substitution markers.
=======================================================

See erime.h


Functions: parsing, deparsing, copying and deleting
===================================================

By default, MA_node objects are constructed to "0".  Various other
constructions are possible - see class definition below.  In these
constructors the integer insist_finite controls the operation of the
object.  If this is 0 then any result from an evaluation will be
returned as is, even if it is NaN or +/- inf.  Setting insist_finite
causes an exception to be thrown if an answer is not a finite real.

The assignment operator is overloaded in the obvious fashion.



Functions: manipulation
=======================

The usual manipulations described in erime.h are available, namely:

void simplify(void);

void replace(..., long n, long m);
void replaceRow(..., long n);
void replaceCol(..., long m);

where details on ... is given below (it is an equation in some format,
eg a string).  fast_, naive_ and fast_naive versions are also
available, and act as would be expected from erime.h



Functions: differentiation and integration
==========================================

The following functions are available:

MA_node diff(long n, long m);
MA_node diffRow(long n, long k, long l);
MA_node diffCol(long i, long j, long m);
MA_node diffMatrix(long i, long j, long k, long l);

MA_node *diffp(long n, long m);
MA_node *diffpRow(long n, long k, long l);
MA_node *diffpCol(long i, long j, long m);
MA_node *diffpMatrix(long i, long j, long k, long l);

as well as the usual fast_ prefixed forms.  The act as would be
expected.  The difference between diff and diffp forms is that the
former will pass the result on the stack ("good" C++ style, fine for
small equations, very slow for large equations) whereas the latter
will allocate once and return a pointer, thus avoiding the repetitive
implicit calling of constructors and such whenever the result gets
passed around.

  

Functions: evaluation
=====================

Evaluation is done via overloading of the () operator (so, if your
MA_node is called fancy, you would evaluate using fancy(x)).  The
following basically follow the standards of erime.h:

double operator()(double **variables) const;
double operator()(void) const;
double operator()(double x) const;
double operator()(double x, double y) const;
double operator()(double x, double y, double z) const;
double operator()(double (*DoubleArgEvaluationFunction)(void *argContents, long i, long j), void *arguments) const;


*********************************************************************/


#ifndef rime_h
#define rime_h

#include <iostream>

//
// To avoid certain repetitive warnings (which are important in context,
// but not to all the other headers that just happen to include this one),
// unless actually compiling this library, just tell the compiler that
// MathsNode is "some struct".  See fixme.txt
//

#ifndef _erime_h
struct mathsNode;
typedef mathsNode MathsNode;
#include "gsl/gsl_math.h"
#include "gsl/err/gsl_errno.h"
#define RIME_NAN GSL_NAN
#endif

int rime_finite(double num);

#define MATHS_ZEROPARSEERR      EERR_MAX_ERROR-1
#define MATHS_DUPLICATEERR      EERR_MAX_ERROR-2
#define MATHS_PARSEERR          EERR_MAX_ERROR-3
#define MATHS_DIFFERR           EERR_MAX_ERROR-4
#define MATHS_NOCONTERR         EERR_MAX_ERROR-5
#define MATHS_CALCERR           EERR_MAX_ERROR-6
#define MATHS_DEPARSEERR        EERR_MAX_ERROR-7
#define MATHS_DECONSTRUCTERR    EERR_MAX_ERROR-8
#define MATHS_DECONSTRUCTERRB   EERR_MAX_ERROR-8


//
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=                                                                +=+=
// +=+=                      MATHS FUNCTION CLASS                      +=+=
// +=+=                                                                +=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
// +=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=
//

class MA_node
{
    friend std::ostream &operator<<(std::ostream &output, const MA_node &dumpee);
    friend std::istream &operator>>(std::istream &input, MA_node &source);

    public:

    //
    // Constructors and destructors.
    //

    MA_node();
    MA_node(int insist_finite);
    MA_node(const double &what, int insist_finite = 0);
    MA_node(const MA_node &source);
    MA_node(const char *express, int insist_finite = 0);

    ~MA_node();

    //
    // Assignment function:
    //

    MA_node &operator=(const double &what);
    MA_node &operator=(const MA_node &source);
    MA_node &operator=(char *express);

    //
    // Equation manipulation functions:
    //

    void simplify(void);

    void replace(const double &with, long n, long m);
    void replace(const MA_node &with, long n, long m);
    void replace(char *express, long n, long m);

    void replaceRow(const double &with, long n);
    void replaceRow(const MA_node &with, long n);
    void replaceRow(char *express, long n);

    void replaceCol(const double &with, long m);
    void replaceCol(const MA_node &with, long m);
    void replaceCol(char *express, long m);

    void fast_replace(const double &with, long n, long m);
    void fast_replace(const MA_node &with, long n, long m);
    void fast_replace(char *express, long n, long m);

    void fast_replaceRow(const double &with, long n);
    void fast_replaceRow(const MA_node &with, long n);
    void fast_replaceRow(char *express, long n);

    void fast_replaceCol(const double &with, long m);
    void fast_replaceCol(const MA_node &with, long m);
    void fast_replaceCol(char *express, long m);

    void naive_replace(const double &with, long n, long m);
    void naive_replace(const MA_node &with, long n, long m);
    void naive_replace(char *express, long n, long m);

    void naive_replaceRow(const double &with, long n);
    void naive_replaceRow(const MA_node &with, long n);
    void naive_replaceRow(char *express, long n);

    void naive_replaceCol(const double &with, long m);
    void naive_replaceCol(const MA_node &with, long m);
    void naive_replaceCol(char *express, long m);

    void fast_naive_replace(const double &with, long n, long m);
    void fast_naive_replace(const MA_node &with, long n, long m);
    void fast_naive_replace(char *express, long n, long m);

    void fast_naive_replaceRow(const double &with, long n);
    void fast_naive_replaceRow(const MA_node &with, long n);
    void fast_naive_replaceRow(char *express, long n);

    void fast_naive_replaceCol(const double &with, long m);
    void fast_naive_replaceCol(const MA_node &with, long m);
    void fast_naive_replaceCol(char *express, long m);

    //
    // Partial differentiation function:
    //

    MA_node diff(long n, long m);
    MA_node diffRow(long n, long k, long l);
    MA_node diffCol(long i, long j, long m);
    MA_node diffMatrix(long i, long j, long k, long l);

    MA_node *diffp(long n, long m);
    MA_node *diffpRow(long n, long k, long l);
    MA_node *diffpCol(long i, long j, long m);
    MA_node *diffpMatrix(long i, long j, long k, long l);

    MA_node fast_diff(long n, long m);
    MA_node fast_diffRow(long n, long k, long l);
    MA_node fast_diffCol(long i, long j, long m);
    MA_node fast_diffMatrix(long i, long j, long k, long l);

    MA_node *fast_diffp(long n, long m);
    MA_node *fast_diffpRow(long n, long k, long l);
    MA_node *fast_diffpCol(long i, long j, long m);
    MA_node *fast_diffpMatrix(long i, long j, long k, long l);

    //
    // Evaluation function:
    //

    double operator()(double **variables) const;
    double operator()(void) const;
    double operator()(double x) const;
    double operator()(double x, double y) const;
    double operator()(double x, double y, double z) const;
    double operator()(double (*DoubleArgEvaluationFunction)(void *argContents, long i, long j), void *arguments) const;

    //
    // Mode control.
    //

    void set_insist_finite();
    void clear_insist_finite();

    private:

    MathsNode *content;

    //
    // If the following is set, results that are NaN or inf will cause
    // an exception.
    //

    int insist_finite;
};








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

std::ostream &operator<<(std::ostream &output, const MA_node    &source);
std::istream &operator>>(std::istream &input, MA_node    &dest);

//
// Mathematical operator overloading
//

// + posation - unary, return rvalue
// - negation - unary, return rvalue

MA_node    operator+ (const MA_node   &left_op);
MA_node    operator- (const MA_node   &left_op);

// + addition       - binary, return rvalue
// - subtraction    - binary, return rvalue

MA_node    operator+ (const MA_node   &left_op, const MA_node   &right_op);
MA_node    operator- (const MA_node   &left_op, const MA_node   &right_op);

// += additive       assignment - binary, return lvalue
// -= subtractive    assignment - binary, return lvalue

MA_node   &operator+=(      MA_node   &left_op, const MA_node   &right_op);
MA_node   &operator-=(      MA_node   &left_op, const MA_node   &right_op);

// * multiplication - binary, return rvalue
// / division       - binary, return rvalue

MA_node    operator* (const MA_node   &left_op, const MA_node   &right_op);
MA_node    operator* (const MA_node   &left_op, const   double  &right_op);
MA_node    operator* (const   double  &left_op, const MA_node   &right_op);

MA_node    operator/ (const MA_node   &left_op, const MA_node   &right_op);
MA_node    operator/ (const MA_node   &left_op, const   double  &right_op);
MA_node    operator/ (const   double  &left_op, const MA_node   &right_op);

// *= multiplicative assignment - binary, return lvalue
// /= divisive       assignment - binary, return lvalue

MA_node   &operator*=(      MA_node   &left_op, const MA_node   &right_op);
MA_node   &operator*=(      MA_node   &left_op, const   double  &right_op);

MA_node   &operator/=(      MA_node   &left_op, const MA_node   &right_op);
MA_node   &operator/=(      MA_node   &left_op, const   double  &right_op);




#endif

