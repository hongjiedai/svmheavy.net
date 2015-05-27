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
// double replacement for counting flops
//
// Written by: Alistair Shilton
//             Melbourne University
//

#include <iostream>
#include <math.h>

#include "c_double.h"
#include "svdefs.h"

c_double::c_double()
{
    value = 0;

    return;
}

c_double::c_double(const int &source)
{
    

    value = source;

    return;
}

c_double::c_double(const long &source)
{
    

    value = source;

    return;
}

c_double::c_double(const double &source)
{
    

    value = source;

    return;
}

c_double::c_double(const c_double &source)
{
    

    value = source.value;

    return;
}

c_double &c_double::operator=(const int &source)
{
    

    value = source;

    return *this;
}

c_double &c_double::operator=(const long &source)
{
    

    value = source;

    return *this;
}

c_double &c_double::operator=(const double &source)
{
    

    value = source;

    return *this;
}

c_double &c_double::operator=(const c_double &source)
{
    

    value = source.value;

    return *this;
}

c_double::operator double() const
{
    

    return value;
}

std::ostream &operator<<(std::ostream &output, const c_double &source)
{
    

    output << source.value;

    return output;
}

std::istream &operator>>(std::istream &input ,       c_double &destin)
{
    

    input  >> destin.value;

    return input;
}

// + posation - unary, return rvalue
// - negation - unary, return rvalue

c_double  operator+ (const c_double &left_op)
{
    

    c_double result;

    result.value = +(left_op.value);

    return result;
}

c_double  operator- (const c_double &left_op)
{
    

    c_double result;

    result.value = -(left_op.value);

    return result;
}

// + addition       - binary, return rvalue
// - subtraction    - binary, return rvalue

c_double  operator+ (const c_double &left_op, const c_double &right_op)
{
    

    #ifdef DO__FLOPS
    REG_f_plus;
    #endif

    c_double result;

    result.value = ( left_op.value + right_op.value );

    return result;
}

c_double  operator+ (const c_double &left_op, const   double &right_op)
{
    

    #ifdef DO__FLOPS
    REG_f_plus;
    #endif

    c_double result;

    result.value = ( left_op.value + right_op );

    return result;
}

c_double  operator+ (const   double &left_op, const c_double &right_op)
{
    

    #ifdef DO__FLOPS
    REG_f_plus;
    #endif

    c_double result;

    result.value = ( left_op + right_op.value );

    return result;
}

c_double  operator- (const c_double &left_op, const c_double &right_op)
{
    

    #ifdef DO__FLOPS
    REG_f_minus;
    #endif

    c_double result;

    result.value = ( left_op.value - right_op.value );

    return result;
}

c_double  operator- (const c_double &left_op, const   double &right_op)
{
    

    #ifdef DO__FLOPS
    REG_f_minus;
    #endif

    c_double result;

    result.value = ( left_op.value - right_op );

    return result;
}

c_double  operator- (const   double &left_op, const c_double &right_op)
{
    

    #ifdef DO__FLOPS
    REG_f_minus;
    #endif

    c_double result;

    result.value = ( left_op - right_op.value );

    return result;
}

// += additive       assignment - binary, return lvalue
// -= subtractive    assignment - binary, return lvalue

c_double &operator+=(      c_double &left_op, const c_double &right_op)
{
    

    #ifdef DO__FLOPS
    REG_f_plus;
    #endif

    left_op.value += right_op.value;

    return left_op;
}

c_double &operator+=(      c_double &left_op, const   double &right_op)
{
    

    #ifdef DO__FLOPS
    REG_f_plus;
    #endif

    left_op.value += right_op;

    return left_op;
}

  double &operator+=(        double &left_op, const c_double &right_op)
{
    

    #ifdef DO__FLOPS
    REG_f_plus;
    #endif

    left_op += right_op.value;

    return left_op;
}

c_double &operator-=(      c_double &left_op, const c_double &right_op)
{
    

    #ifdef DO__FLOPS
    REG_f_minus;
    #endif

    left_op.value -= right_op.value;

    return left_op;
}

c_double &operator-=(      c_double &left_op, const   double &right_op)
{
    

    #ifdef DO__FLOPS
    REG_f_minus;
    #endif

    left_op.value -= right_op;

    return left_op;
}

  double &operator-=(        double &left_op, const c_double &right_op)
{
    

    #ifdef DO__FLOPS
    REG_f_minus;
    #endif

    left_op -= right_op.value;

    return left_op;
}

// * multiplication - binary, return rvalue
// / division       - binary, return rvalue

c_double  operator* (const c_double &left_op, const c_double &right_op)
{
    

    #ifdef DO__FLOPS
    REG_f_mult;
    #endif

    c_double result;

    result.value = ( left_op.value * right_op.value );

    return result;
}

c_double  operator* (const c_double &left_op, const   double &right_op)
{
    

    #ifdef DO__FLOPS
    REG_f_mult;
    #endif

    c_double result;

    result.value = ( left_op.value * right_op );

    return result;
}

c_double  operator* (const   double &left_op, const c_double &right_op)
{
    

    #ifdef DO__FLOPS
    REG_f_mult;
    #endif

    c_double result;

    result.value = ( left_op * right_op.value );

    return result;
}

c_double  operator/ (const c_double &left_op, const c_double &right_op)
{
    

    #ifdef DO__FLOPS
    REG_f_div;
    #endif

    c_double result;

    result.value = ( left_op.value / right_op.value );

    return result;
}

c_double  operator/ (const c_double &left_op, const   double &right_op)
{
    

    #ifdef DO__FLOPS
    REG_f_div;
    #endif

    c_double result;

    result.value = ( left_op.value / right_op );

    return result;
}

c_double  operator/ (const   double &left_op, const c_double &right_op)
{
    

    #ifdef DO__FLOPS
    REG_f_div;
    #endif

    c_double result;

    result.value = ( left_op / right_op.value );

    return result;
}

// *= multiplicative assignment - binary, return lvalue
// /= divisive       assignment - binary, return lvalue

c_double &operator*=(      c_double &left_op, const c_double &right_op)
{
    

    #ifdef DO__FLOPS
    REG_f_mult;
    #endif

    left_op.value *= right_op.value;

    return left_op;
}

c_double &operator*=(      c_double &left_op, const   double &right_op)
{
    

    #ifdef DO__FLOPS
    REG_f_mult;
    #endif

    left_op.value *= right_op;

    return left_op;
}

  double &operator*=(        double &left_op, const c_double &right_op)
{
    

    #ifdef DO__FLOPS
    REG_f_mult;
    #endif

    left_op *= right_op.value;

    return left_op;
}

c_double &operator/=(      c_double &left_op, const c_double &right_op)
{
    

    #ifdef DO__FLOPS
    REG_f_div;
    #endif

    left_op.value /= right_op.value;

    return left_op;
}

c_double &operator/=(      c_double &left_op, const   double &right_op)
{
    

    #ifdef DO__FLOPS
    REG_f_div;
    #endif

    left_op.value /= right_op;

    return left_op;
}

  double &operator/=(        double &left_op, const c_double &right_op)
{
    

    #ifdef DO__FLOPS
    REG_f_div;
    #endif

    left_op /= right_op.value;

    return left_op;
}

//
// Relational operator overloading
//

// == equivalence

int operator==(const c_double &left_op, const c_double &right_op)
{
    

    return ( (left_op.value) == (right_op.value) );
}

int operator==(const c_double &left_op, const   double &right_op)
{
    

    return ( (left_op.value) == (right_op) );
}

int operator==(const   double &left_op, const c_double &right_op)
{
    

    return ( (left_op) == (right_op.value) );
}

// != inequivalence

int operator!=(const c_double &left_op, const c_double &right_op)
{
    

    return ( (left_op.value) != (right_op.value) );
}

int operator!=(const c_double &left_op, const   double &right_op)
{
    

    return ( (left_op.value) != (right_op) );
}

int operator!=(const   double &left_op, const c_double &right_op)
{
    

    return ( (left_op) != (right_op.value) );
}

// < less than

int operator< (const c_double &left_op, const c_double &right_op)
{
    

    return ( (left_op.value) < (right_op.value) );
}

int operator< (const c_double &left_op, const   double &right_op)
{
    

    return ( (left_op.value) < (right_op) );
}

int operator< (const   double &left_op, const c_double &right_op)
{
    

    return ( (left_op) < (right_op.value) );
}

// <= less than or equal to

int operator<=(const c_double &left_op, const c_double &right_op)
{
    

    return ( (left_op.value) <= (right_op.value) );
}

int operator<=(const c_double &left_op, const   double &right_op)
{
    

    return ( (left_op.value) <= (right_op) );
}

int operator<=(const   double &left_op, const c_double &right_op)
{
    

    return ( (left_op) <= (right_op.value) );
}

// > greater than

int operator> (const c_double &left_op, const c_double &right_op)
{
    

    return ( (left_op.value) > (right_op.value) );
}

int operator> (const c_double &left_op, const   double &right_op)
{
    

    return ( (left_op.value) > (right_op) );
}

int operator> (const   double &left_op, const c_double &right_op)
{
    

    return ( (left_op) > (right_op.value) );
}

// >= greater than or equal to

int operator>=(const c_double &left_op, const c_double &right_op)
{
    

    return ( (left_op.value) >= (right_op.value) );
}

int operator>=(const c_double &left_op, const   double &right_op)
{
    

    return ( (left_op.value) >= (right_op) );
}

int operator>=(const   double &left_op, const c_double &right_op)
{
    

    return ( (left_op) >= (right_op.value) );
}


// fpu operations

c_double pow(c_double x,c_double y)
{
    

    #ifdef DO__FLOPS
    REG_f_fpu;
    #endif

    c_double result;

    result = pow((double) x,(double) y);

    return result;
}

c_double pow(c_double x,double y)
{
    

    #ifdef DO__FLOPS
    REG_f_fpu;
    #endif

    c_double result;

    result = pow((double) x,y);

    return result;
}

c_double pow(double x,c_double y)
{
    

    #ifdef DO__FLOPS
    REG_f_fpu;
    #endif

    c_double result;

    result = pow(x,(double) y);

    return result;
}

c_double sin(c_double x)
{
    

    #ifdef DO__FLOPS
    REG_f_fpu;
    #endif

    c_double result;

    result = sin((double) x);

    return result;
    
}

c_double cos(c_double x)
{
    

    #ifdef DO__FLOPS
    REG_f_fpu;
    #endif

    c_double result;

    result = cos((double) x);

    return result;
}

c_double tan(c_double x)
{
    

    #ifdef DO__FLOPS
    REG_f_fpu;
    #endif

    c_double result;

    result = tan((double) x);

    return result;
}

c_double asin(c_double x)
{
    

    #ifdef DO__FLOPS
    REG_f_fpu;
    #endif

    c_double result;

    result = asin((double) x);

    return result;
}

c_double acos(c_double x)
{
    

    #ifdef DO__FLOPS
    REG_f_fpu;
    #endif

    c_double result;

    result = acos((double) x);

    return result;
}

c_double atan(c_double x)
{
    

    #ifdef DO__FLOPS
    REG_f_fpu;
    #endif

    c_double result;

    result = atan((double) x);

    return result;
}

c_double sinh(c_double x)
{
    

    #ifdef DO__FLOPS
    REG_f_fpu;
    #endif

    c_double result;

    result = sinh((double) x);

    return result;
}

c_double cosh(c_double x)
{
    

    #ifdef DO__FLOPS
    REG_f_fpu;
    #endif

    c_double result;

    result = cosh((double) x);

    return result;
}

c_double tanh(c_double x)
{
    

    #ifdef DO__FLOPS
    REG_f_fpu;
    #endif

    c_double result;

    result = tanh((double) x);

    return result;
}

c_double asinh(c_double x)
{
    

    #ifdef DO__FLOPS
    REG_f_fpu;
    #endif

    c_double result;

    result = asinh((double) x);

    return result;
}

c_double acosh(c_double x)
{
    

    #ifdef DO__FLOPS
    REG_f_fpu;
    #endif

    c_double result;

    result = acosh((double) x);

    return result;
}

c_double atanh(c_double x)
{
    

    #ifdef DO__FLOPS
    REG_f_fpu;
    #endif

    c_double result;

    result = atanh((double) x);

    return result;
}

c_double erf(c_double x)
{
    

    #ifdef DO__FLOPS
    REG_f_fpu;
    #endif

    c_double result;

    result = erf((double) x);

    return result;
}

c_double exp(c_double x)
{
    

    #ifdef DO__FLOPS
    REG_f_fpu;
    #endif

    c_double result;

    result = exp((double) x);

    return result;
}

c_double expm1(c_double x)
{
    

    #ifdef DO__FLOPS
    REG_f_fpu;
    #endif

    c_double result;

    result = expm1((double) x);

    return result;
}

c_double log(c_double x)
{
    

    #ifdef DO__FLOPS
    REG_f_fpu;
    #endif

    c_double result;

    result = log((double) x);

    return result;
}

c_double log1p(c_double x)
{
    

    #ifdef DO__FLOPS
    REG_f_fpu;
    #endif

    c_double result;

    result = log1p((double) x);

    return result;
}

c_double log10(c_double x)
{
    

    #ifdef DO__FLOPS
    REG_f_fpu;
    #endif

    c_double result;

    result = log10((double) x);

    return result;
}

c_double sqrt(c_double x)
{
    

    #ifdef DO__FLOPS
    REG_f_sqrt;
    #endif

    c_double result;

    result = sqrt((double) x);

    return result;
}

c_double gamma(c_double x)
{
    

    #ifdef DO__FLOPS
    REG_f_fpu;
    #endif

    c_double result;

    result = gamma((double) x);

    return result;
}

c_double jN(int x, c_double y)
{
    

    #ifdef DO__FLOPS
    REG_f_fpu;
    #endif

    c_double result;

    result = jN(x,(double) y);

    return result;
}

c_double yN(int x, c_double y)
{
    

    #ifdef DO__FLOPS
    REG_f_fpu;
    #endif

    c_double result;

    result = yN(x,(double) y);

    return result;
}




