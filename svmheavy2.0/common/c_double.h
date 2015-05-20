
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

#ifndef _c_double_h
#define _c_double_h

#include <iostream>
#include <math.h>

#include "svdefs.h"

//
// Count double type - used to count flops transparently
//

class c_double
{
    friend std::ostream &operator<<(std::ostream &output, const c_double &source);
    friend std::istream &operator>>(std::istream &input,  c_double       &destin);

    friend c_double  operator+ (const c_double &left_op);
    friend c_double  operator- (const c_double &left_op);
    friend c_double  operator+ (const c_double &left_op, const c_double &right_op);
    friend c_double  operator+ (const c_double &left_op, const   double &right_op);
    friend c_double  operator+ (const   double &left_op, const c_double &right_op);
    friend c_double  operator- (const c_double &left_op, const c_double &right_op);
    friend c_double  operator- (const c_double &left_op, const   double &right_op);
    friend c_double  operator- (const   double &left_op, const c_double &right_op);
    friend c_double &operator+=(      c_double &left_op, const c_double &right_op);
    friend c_double &operator+=(      c_double &left_op, const   double &right_op);
    friend   double &operator+=(        double &left_op, const c_double &right_op);
    friend c_double &operator-=(      c_double &left_op, const c_double &right_op);
    friend c_double &operator-=(      c_double &left_op, const   double &right_op);
    friend   double &operator-=(        double &left_op, const c_double &right_op);
    friend c_double  operator* (const c_double &left_op, const c_double &right_op);
    friend c_double  operator* (const c_double &left_op, const   double &right_op);
    friend c_double  operator* (const   double &left_op, const c_double &right_op);
    friend c_double  operator/ (const c_double &left_op, const c_double &right_op);
    friend c_double  operator/ (const c_double &left_op, const   double &right_op);
    friend c_double  operator/ (const   double &left_op, const c_double &right_op);
    friend c_double &operator*=(      c_double &left_op, const c_double &right_op);
    friend c_double &operator*=(      c_double &left_op, const   double &right_op);
    friend   double &operator*=(        double &left_op, const c_double &right_op);
    friend c_double &operator/=(      c_double &left_op, const c_double &right_op);
    friend c_double &operator/=(      c_double &left_op, const   double &right_op);
    friend   double &operator/=(        double &left_op, const c_double &right_op);

    friend int operator==(const c_double &left_op, const c_double &right_op);
    friend int operator==(const c_double &left_op, const   double &right_op);
    friend int operator==(const   double &left_op, const c_double &right_op);
    friend int operator!=(const c_double &left_op, const c_double &right_op);
    friend int operator!=(const c_double &left_op, const   double &right_op);
    friend int operator!=(const   double &left_op, const c_double &right_op);
    friend int operator< (const c_double &left_op, const c_double &right_op);
    friend int operator< (const c_double &left_op, const   double &right_op);
    friend int operator< (const   double &left_op, const c_double &right_op);
    friend int operator<=(const c_double &left_op, const c_double &right_op);
    friend int operator<=(const c_double &left_op, const   double &right_op);
    friend int operator<=(const   double &left_op, const c_double &right_op);
    friend int operator> (const c_double &left_op, const c_double &right_op);
    friend int operator> (const c_double &left_op, const   double &right_op);
    friend int operator> (const   double &left_op, const c_double &right_op);
    friend int operator>=(const c_double &left_op, const c_double &right_op);
    friend int operator>=(const c_double &left_op, const   double &right_op);
    friend int operator>=(const   double &left_op, const c_double &right_op);

    public:

    c_double();
    c_double(const int &);
    c_double(const long &);
    c_double(const double &);
    c_double(const c_double &);

    c_double &operator=(const int &);
    c_double &operator=(const long &);
    c_double &operator=(const double &);
    c_double &operator=(const c_double &);

    operator double() const;

    private:

    double value;
};

std::ostream &operator<<(std::ostream &output, const c_double &source);
std::istream &operator>>(std::istream &input ,       c_double &destin);

// + posation - unary, return rvalue
// - negation - unary, return rvalue

c_double  operator+ (const c_double &left_op); // ign
c_double  operator- (const c_double &left_op); // ign

// + addition       - binary, return rvalue
// - subtraction    - binary, return rvalue

c_double  operator+ (const c_double &left_op, const c_double &right_op);
c_double  operator+ (const c_double &left_op, const   double &right_op);
c_double  operator+ (const   double &left_op, const c_double &right_op);

c_double  operator- (const c_double &left_op, const c_double &right_op);
c_double  operator- (const c_double &left_op, const   double &right_op);
c_double  operator- (const   double &left_op, const c_double &right_op);

// += additive       assignment - binary, return lvalue
// -= subtractive    assignment - binary, return lvalue

c_double &operator+=(      c_double &left_op, const c_double &right_op);
c_double &operator+=(      c_double &left_op, const   double &right_op);
  double &operator+=(        double &left_op, const c_double &right_op);

c_double &operator-=(      c_double &left_op, const c_double &right_op);
c_double &operator-=(      c_double &left_op, const   double &right_op);
  double &operator-=(        double &left_op, const c_double &right_op);

// * multiplication - binary, return rvalue
// / division       - binary, return rvalue

c_double  operator* (const c_double &left_op, const c_double &right_op);
c_double  operator* (const c_double &left_op, const   double &right_op);
c_double  operator* (const   double &left_op, const c_double &right_op);

c_double  operator/ (const c_double &left_op, const c_double &right_op);
c_double  operator/ (const c_double &left_op, const   double &right_op);
c_double  operator/ (const   double &left_op, const c_double &right_op);

// *= multiplicative assignment - binary, return lvalue
// /= divisive       assignment - binary, return lvalue

c_double &operator*=(      c_double &left_op, const c_double &right_op);
c_double &operator*=(      c_double &left_op, const   double &right_op);
  double &operator*=(        double &left_op, const c_double &right_op);

c_double &operator/=(      c_double &left_op, const c_double &right_op);
c_double &operator/=(      c_double &left_op, const   double &right_op);
  double &operator/=(        double &left_op, const c_double &right_op);

//
// Relational operator overloading
//

// == equivalence

int operator==(const c_double &left_op, const c_double &right_op);
int operator==(const c_double &left_op, const   double &right_op);
int operator==(const   double &left_op, const c_double &right_op);

// != inequivalence

int operator!=(const c_double &left_op, const c_double &right_op);
int operator!=(const c_double &left_op, const   double &right_op);
int operator!=(const   double &left_op, const c_double &right_op);

// < less than

int operator< (const c_double &left_op, const c_double &right_op);
int operator< (const c_double &left_op, const   double &right_op);
int operator< (const   double &left_op, const c_double &right_op);

// <= less than or equal to

int operator<=(const c_double &left_op, const c_double &right_op);
int operator<=(const c_double &left_op, const   double &right_op);
int operator<=(const   double &left_op, const c_double &right_op);

// > greater than

int operator> (const c_double &left_op, const c_double &right_op);
int operator> (const c_double &left_op, const   double &right_op);
int operator> (const   double &left_op, const c_double &right_op);

// >= greater than or equal to

int operator>=(const c_double &left_op, const c_double &right_op);
int operator>=(const c_double &left_op, const   double &right_op);
int operator>=(const   double &left_op, const c_double &right_op);


// fpu operations

c_double pow(c_double x,c_double y);
c_double pow(c_double x,double y);
c_double pow(double x,c_double y);
c_double sin(c_double x);
c_double cos(c_double x);
c_double tan(c_double x);
c_double asin(c_double x);
c_double acos(c_double x);
c_double atan(c_double x);
c_double sinh(c_double x);
c_double cosh(c_double x);
c_double tanh(c_double x);
c_double asinh(c_double x);
c_double acosh(c_double x);
c_double atanh(c_double x);
c_double erf(c_double x);
c_double exp(c_double x);
c_double expm1(c_double x);
c_double log(c_double x);
c_double log1p(c_double x);
c_double log10(c_double x);
c_double sqrt(c_double x);
c_double gamma(c_double x);
c_double jN(int x, c_double y);
c_double yN(int x, c_double y);



#endif
