#include "stdafx.h"
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

//
// C++ wrapper for Runtime Mathematical Equation Library.  See erime.h
// for details on operation of this library.
//
// Written by: Alistair Shilton
//             Melbourne University
//             apsh@ecr.mu.oz.au
//
// Date of completion: 13/3/2008
//


#include <iostream>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <malloc.h>
#include "gsl/gsl_math.h"
#include "gsl/err/gsl_errno.h"

extern "C"
{
#include "erime.h"
}

#include "rime.h"

#define DEFAULT_STDIN_LEN       256000

//
// Constructors and destructors.
//

MA_node::MA_node()
{
    insist_finite = 0;

    if ( ( content = getZero() ) == NULL )
    {
        throw MATHS_ZEROPARSEERR;
    }

    return;
}

MA_node::MA_node(int _insist_finite)
{
    insist_finite = _insist_finite;

    if ( ( content = getZero() ) == NULL )
    {
        throw MATHS_ZEROPARSEERR;
    }

    return;
}

MA_node::MA_node(const MA_node &source)
{
    insist_finite = source.insist_finite;

    if ( ( content = copyEquation(source.content) ) == NULL )
    {
        throw MATHS_DUPLICATEERR;
    }

    return;
}

MA_node::MA_node(const char *express, int _insist_finite)
{
    char *copexpress = new char[strlen(express)+1];;

    strcpy(copexpress,express);

    insist_finite = _insist_finite;

    if ( ( content = parseMaths(copexpress) ) == NULL )
    {
        throw MATHS_PARSEERR;
    }

    delete copexpress;

    return;
}

MA_node::MA_node(const double &what, int _insist_finite)
{
    insist_finite = _insist_finite;

    if ( ( content = getZero() ) == NULL )
    {
        throw MATHS_ZEROPARSEERR;
    }

    *this = what;

    return;
}

MA_node::~MA_node()
{
    if ( content != NULL )
    {
        delEqn(content);
    }

    return;
}




//
// Assignment functions:
//

MA_node &MA_node::operator=(const double &what)
{
    char dummy[50];

    if ( content != NULL )
    {
        delEqn(content);
    }

    sprintf(dummy,"%f",what);

    if ( ( content = parseMaths(dummy) ) == NULL )
    {
        throw MATHS_PARSEERR;
    }

    return (*this);
}

MA_node &MA_node::operator=(const MA_node &source)
{
    if ( this != &source )
    {
        if ( content != NULL )
        {
            delEqn(content);
        }

        if ( ( content = copyEquation(source.content) ) == NULL )
        {
            throw MATHS_DUPLICATEERR;
        }

        insist_finite = source.insist_finite;
    }

    return (*this);
}

MA_node &MA_node::operator=(char *express)
{
    if ( content != NULL )
    {
        delEqn(content);
    }

    if ( ( content = parseMaths(express) ) == NULL )
    {
        throw MATHS_PARSEERR;
    }

    return (*this);
}




//
// Simplification function:
//

void MA_node::simplify(void)
{
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = simplifyEqn(content) ) )
    {
        throw tamp;
    }

    return;
}




//
// Substitution function:
//

void MA_node::replace(const double &with, long n, long m)
{
    int tamp;
    MA_node temp(with);

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = subVar(content,temp.content,n,m,0) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::replace(const MA_node &with, long n, long m)
{
    // Allow that we may be replacing recursively.

    MA_node temp(with);
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = subVar(content,temp.content,n,m,0) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::replace(char *express, long n, long m)
{
    MA_node temp(express);
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = subVar(content,temp.content,n,m,0) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::replaceRow(const double &with, long n)
{
    MA_node temp(with);
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = subVarRow(content,temp.content,n,0) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::replaceRow(const MA_node &with, long n)
{
    // Allow that we may be replacing with itself (recursive)

    MA_node temp(with);
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = subVarRow(content,temp.content,n,0) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::replaceRow(char *express, long n)
{
    MA_node temp(express);
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = subVarRow(content,temp.content,n,0) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::replaceCol(const double &with, long m)
{
    MA_node temp(with);
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = subVarCol(content,temp.content,m,0) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::replaceCol(const MA_node &with, long m)
{
    // Allow that we may be replacing with itself (recursive)

    MA_node temp(with);
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = subVarCol(content,temp.content,m,0) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::replaceCol(char *express, long m)
{
    MA_node temp(express);
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = subVarCol(content,temp.content,m,0) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::fast_replace(const double &with, long n, long m)
{
    int tamp;
    MA_node temp(with);

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = fast_subVar(content,temp.content,n,m,0) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::fast_replace(const MA_node &with, long n, long m)
{
    // Allow that we may be replacing with itself (recursive)

    MA_node temp(with);
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = fast_subVar(content,temp.content,n,m,0) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::fast_replace(char *express, long n, long m)
{
    MA_node temp(express);
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = fast_subVar(content,temp.content,n,m,0) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::fast_replaceRow(const double &with, long n)
{
    MA_node temp(with);
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = fast_subVarRow(content,temp.content,n,0) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::fast_replaceRow(const MA_node &with, long n)
{
    // Allow that we may be replacing with itself (recursive)

    MA_node temp(with);
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = fast_subVarRow(content,temp.content,n,0) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::fast_replaceRow(char *express, long n)
{
    MA_node temp(express);
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = fast_subVarRow(content,temp.content,n,0) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::fast_replaceCol(const double &with, long m)
{
    MA_node temp(with);
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = fast_subVarCol(content,temp.content,m,0) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::fast_replaceCol(const MA_node &with, long m)
{
    // Allow that we may be replacing with itself (recursive)

    MA_node temp(with);
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = fast_subVarCol(content,temp.content,m,0) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::fast_replaceCol(char *express, long m)
{
    MA_node temp(express);
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = fast_subVarCol(content,temp.content,m,0) ) )
    {
        throw tamp;
    }

    return;
}




//
// Naive substitution function:
//

void MA_node::naive_replace(const double &with, long n, long m)
{
    int tamp;
    MA_node temp(with);

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = subVar(content,temp.content,n,m,1) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::naive_replace(const MA_node &with, long n, long m)
{
    // Allow that we may be replacing with itself (recursive)

    MA_node temp(with);
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = subVar(content,temp.content,n,m,1) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::naive_replace(char *express, long n, long m)
{
    MA_node temp(express);
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = subVar(content,temp.content,n,m,1) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::naive_replaceRow(const double &with, long n)
{
    MA_node temp(with);
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = subVarRow(content,temp.content,n,1) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::naive_replaceRow(const MA_node &with, long n)
{
    // Allow that we may be replacing with itself (recursive)

    MA_node temp(with);
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = subVarRow(content,temp.content,n,1) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::naive_replaceRow(char *express, long n)
{
    MA_node temp(express);
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = subVarRow(content,temp.content,n,1) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::naive_replaceCol(const double &with, long m)
{
    MA_node temp(with);
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = subVarCol(content,temp.content,m,1) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::naive_replaceCol(const MA_node &with, long m)
{
    // Allow that we may be replacing with itself (recursive)

    MA_node temp(with);
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = subVarCol(content,temp.content,m,1) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::naive_replaceCol(char *express, long m)
{
    MA_node temp(express);
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = subVarCol(content,temp.content,m,1) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::fast_naive_replace(const double &with, long n, long m)
{
    int tamp;
    MA_node temp(with);

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = fast_subVar(content,temp.content,n,m,1) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::fast_naive_replace(const MA_node &with, long n, long m)
{
    // Allow that we may be replacing with itself (recursive)

    MA_node temp(with);
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = fast_subVar(content,temp.content,n,m,1) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::fast_naive_replace(char *express, long n, long m)
{
    MA_node temp(express);
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = fast_subVar(content,temp.content,n,m,1) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::fast_naive_replaceRow(const double &with, long n)
{
    MA_node temp(with);
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = fast_subVarRow(content,temp.content,n,1) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::fast_naive_replaceRow(const MA_node &with, long n)
{
    // Allow that we may be replacing with itself (recursive)

    MA_node temp(with);
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = fast_subVarRow(content,temp.content,n,1) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::fast_naive_replaceRow(char *express, long n)
{
    MA_node temp(express);
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = fast_subVarRow(content,temp.content,n,1) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::fast_naive_replaceCol(const double &with, long m)
{
    MA_node temp(with);
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = fast_subVarCol(content,temp.content,m,1) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::fast_naive_replaceCol(const MA_node &with, long m)
{
    // Allow that we may be replacing with itself (recursive)

    MA_node temp(with);
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = fast_subVarCol(content,temp.content,m,1) ) )
    {
        throw tamp;
    }

    return;
}

void MA_node::fast_naive_replaceCol(char *express, long m)
{
    MA_node temp(express);
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = fast_subVarCol(content,temp.content,m,1) ) )
    {
        throw tamp;
    }

    return;
}





//
// Partial differentiation function:
//

MA_node MA_node::diff(long n, long m)
{
    MA_node result;

    if ( result.content != NULL )
    {
        delEqn(result.content);
    }

    if ( ( result.content = makeDeriv(content,n,m) ) == NULL )
    {
        throw MATHS_DIFFERR;
    }

    return result;
}

MA_node MA_node::diffRow(long n, long k, long l)
{
    MA_node result;

    if ( result.content != NULL )
    {
        delEqn(result.content);
    }

    if ( ( result.content = RowMakeDeriv(content,n,k,l) ) == NULL )
    {
        throw MATHS_DIFFERR;
    }

    return result;
}

MA_node MA_node::diffCol(long i, long j, long m)
{
    MA_node result;

    if ( result.content != NULL )
    {
        delEqn(result.content);
    }

    if ( ( result.content = ColMakeDeriv(content,i,j,m) ) == NULL )
    {
        throw MATHS_DIFFERR;
    }

    return result;
}

MA_node MA_node::diffMatrix(long i, long j, long k, long l)
{
    MA_node result;

    if ( result.content != NULL )
    {
        delEqn(result.content);
    }

    if ( ( result.content = MatrixMakeDeriv(content,i,j,k,l) ) == NULL )
    {
        throw MATHS_DIFFERR;
    }

    return result;
}

MA_node *MA_node::diffp(long n, long m)
{
    MA_node *result;

    result = new MA_node;

    if ( result->content != NULL )
    {
        delEqn(result->content);
    }

    if ( ( result->content = makeDeriv(content,n,m) ) == NULL )
    {
        throw MATHS_DIFFERR;
    }

    return result;
}

MA_node *MA_node::diffpRow(long n, long k, long l)
{
    MA_node *result;

    result = new MA_node;

    if ( result->content != NULL )
    {
        delEqn(result->content);
    }

    if ( ( result->content = RowMakeDeriv(content,n,k,l) ) == NULL )
    {
        throw MATHS_DIFFERR;
    }

    return result;
}

MA_node *MA_node::diffpCol(long i, long j, long m)
{
    MA_node *result;

    result = new MA_node;

    if ( result->content != NULL )
    {
        delEqn(result->content);
    }

    if ( ( result->content = ColMakeDeriv(content,i,j,m) ) == NULL )
    {
        throw MATHS_DIFFERR;
    }

    return result;
}

MA_node *MA_node::diffpMatrix(long i, long j, long k, long l)
{
    MA_node *result;

    result = new MA_node;

    if ( result->content != NULL )
    {
        delEqn(result->content);
    }

    if ( ( result->content = MatrixMakeDeriv(content,i,j,k,l) ) == NULL )
    {
        throw MATHS_DIFFERR;
    }

    return result;
}

MA_node MA_node::fast_diff(long n, long m)
{
    MA_node result;

    if ( result.content != NULL )
    {
        delEqn(result.content);
    }

    if ( ( result.content = fast_makeDeriv(content,n,m) ) == NULL )
    {
        throw MATHS_DIFFERR;
    }

    return result;
}

MA_node MA_node::fast_diffRow(long n, long k, long l)
{
    MA_node result;

    if ( result.content != NULL )
    {
        delEqn(result.content);
    }

    if ( ( result.content = fast_RowMakeDeriv(content,n,k,l) ) == NULL )
    {
        throw MATHS_DIFFERR;
    }

    return result;
}

MA_node MA_node::fast_diffCol(long i, long j, long m)
{
    MA_node result;

    if ( result.content != NULL )
    {
        delEqn(result.content);
    }

    if ( ( result.content = fast_ColMakeDeriv(content,i,j,m) ) == NULL )
    {
        throw MATHS_DIFFERR;
    }

    return result;
}

MA_node MA_node::fast_diffMatrix(long i, long j, long k, long l)
{
    MA_node result;

    if ( result.content != NULL )
    {
        delEqn(result.content);
    }

    if ( ( result.content = fast_MatrixMakeDeriv(content,i,j,k,l) ) == NULL )
    {
        throw MATHS_DIFFERR;
    }

    return result;
}

MA_node *MA_node::fast_diffp(long n, long m)
{
    MA_node *result;

    result = new MA_node;

    if ( result->content != NULL )
    {
        delEqn(result->content);
    }

    if ( ( result->content = fast_makeDeriv(content,n,m) ) == NULL )
    {
        throw MATHS_DIFFERR;
    }

    return result;
}

MA_node *MA_node::fast_diffpRow(long n, long k, long l)
{
    MA_node *result;

    result = new MA_node;

    if ( result->content != NULL )
    {
        delEqn(result->content);
    }

    if ( ( result->content = fast_RowMakeDeriv(content,n,k,l) ) == NULL )
    {
        throw MATHS_DIFFERR;
    }

    return result;
}

MA_node *MA_node::fast_diffpCol(long i, long j, long m)
{
    MA_node *result;

    result = new MA_node;

    if ( result->content != NULL )
    {
        delEqn(result->content);
    }

    if ( ( result->content = fast_ColMakeDeriv(content,i,j,m) ) == NULL )
    {
        throw MATHS_DIFFERR;
    }

    return result;
}

MA_node *MA_node::fast_diffpMatrix(long i, long j, long k, long l)
{
    MA_node *result;

    result = new MA_node;

    if ( result->content != NULL )
    {
        delEqn(result->content);
    }

    if ( ( result->content = fast_MatrixMakeDeriv(content,i,j,k,l) ) == NULL )
    {
        throw MATHS_DIFFERR;
    }

    return result;
}




//
// Evaluate function:
//

double MA_node::operator()(double **variables) const
{
    double result;
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = evaluateEqnMatrix_e(&result,content,variables) ) )
    {
        throw tamp;
    }

    if ( !gsl_finite(result) && insist_finite )
    {
        throw MATHS_CALCERR;
    }

    return result;
}

double MA_node::operator()(double x) const
{
    double result;
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = evaluateEqnList_e(&result,content,1,x) ) )
    {
        throw tamp;
    }

    if ( !gsl_finite(result) && insist_finite )
    {
        throw MATHS_CALCERR;
    }

    return result;
}

double MA_node::operator()(void) const
{
    double result;
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = evaluateEqnNull_e(&result,content) ) )
    {
        throw tamp;
    }

    if ( !gsl_finite(result) && insist_finite )
    {
        throw MATHS_CALCERR;
    }

    return result;
}

double MA_node::operator()(double x, double y) const
{
    double result;
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = evaluateEqnList_e(&result,content,2,x,y) ) )
    {
        throw tamp;
    }

    if ( !gsl_finite(result) && insist_finite )
    {
        throw MATHS_CALCERR;
    }

    return result;
}

double MA_node::operator()(double x, double y, double z) const
{
    double result;
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = evaluateEqnList_e(&result,content,3,x,y,z) ) )
    {
        throw tamp;
    }

    if ( !gsl_finite(result) && insist_finite )
    {
        throw MATHS_CALCERR;
    }

    return result;
}

double MA_node::operator()(double (*DoubleArgEvaluationFunction)(void *argContents, long i, long j), void *arguments) const
{
    double result;
    int tamp;

    if ( content == NULL )
    {
        throw MATHS_NOCONTERR;
    }

    if ( ( tamp = evaluateEqnGeneral_e(&result,content,DoubleArgEvaluationFunction,arguments) ) )
    {
        throw tamp;
    }

    if ( !gsl_finite(result) && insist_finite )
    {
        throw MATHS_CALCERR;
    }

    return result;
}




void MA_node::set_insist_finite()
{
    insist_finite = 1;

    return;
}

void MA_node::clear_insist_finite()
{
    insist_finite = 0;

    return;
}





// Stream io

std::ostream &operator<<(std::ostream &output, const MA_node &dumpee)
{
    char *temp;

    if ( dumpee.content != NULL )
    {
        temp = deparseMaths(dumpee.content);

        if ( temp != NULL )
        {
            output << "#" << strlen(temp) << ": " << temp << "\n";

            // use free, as this string would have been allocated with new

            free(temp);
        }

        else
        {
            throw MATHS_DEPARSEERR;
        }
    }

    return output;
}

std::istream &operator>>(std::istream &input, MA_node &source)
{
    char temp[DEFAULT_STDIN_LEN+1];
    char *tempb;
    long len;

    // We allow two possible formats here.  The simplest is just a string
    // that is the equation (without spaces of course), which we read and
    // then parse.  This is only intended for small equations (specifically,
    // equations of length <= DEFAULT_STDIN_LEN.
    //
    // The second format (which is output by default) can deal with equations
    // of arbitrary length.  This format has #n:, followed by a string of
    // length n that is the equation, without spaces (the #n: string must be
    // of length no more than DEFAULT_STDIN_LEN, but this shouldn't be too
    // restrictive).

    if ( source.content != NULL )
    {
        delEqn(source.content);
    }

    input >> temp;

    if ( temp[0] != '#' )
    {
        // Case 1: string is just an equation.

        if ( ( source.content = parseMaths(temp) ) == NULL )
        {
            throw MATHS_PARSEERR;
        }
    }

    else
    {
        // Case 2: string has length specifier.  Prune #n: to #n and then
        //         read n using sscanf on the incremented string pointer.

        temp[strlen(temp)-1] = '\0';

        sscanf(temp+1,"%ld",&len);

        tempb = new char[len+1];

        input >> tempb;

        if ( ( source.content = parseMaths(tempb) ) == NULL )
        {
            throw MATHS_PARSEERR;
        }

        delete[] tempb;
    }

    return input;
}




// + posation - unary, return rvalue
// - negation - unary, return rvalue

MA_node  operator+ (const MA_node &left_op)
{
    return ( left_op * +1.0 );
}

MA_node  operator- (const MA_node &left_op)
{
    return ( left_op * -1.0 );
}




// + addition    - binary, return rvalue
// - subtraction - binary, return rvalue

MA_node  operator+ (const MA_node          &left_op, const MA_node          &right_op)
{
    MA_node temp("var(10000,10000)+var(20000,20000)");

    temp.replace(left_op,10000,10000);
    temp.replace(right_op,20000,20000);

    return temp;
}

MA_node  operator- (const MA_node          &left_op, const MA_node          &right_op)
{
    MA_node temp("var(100000,100000)-var(200000,200000)");

    temp.replace(left_op,100000,100000);
    temp.replace(right_op,200000,200000);

    return temp;
}




// += additive       assignment - binary, return lvalue
// -= subtractive    assignment - binary, return lvalue

MA_node &operator+=(      MA_node &left_op, const MA_node          &right_op)
{
    left_op = ( left_op + right_op );

    return left_op;
}

MA_node &operator-=(      MA_node &left_op, const MA_node          &right_op)
{
    left_op = ( left_op - right_op );

    return left_op;
}




// * multiplication - binary, return rvalue
// / division       - binary, return rvalue

MA_node  operator* (const MA_node          &left_op, const MA_node          &right_op)
{
    MA_node temp("var(100000,100000)*var(200000,200000)");

    temp.replace(left_op,100000,100000);
    temp.replace(right_op,200000,200000);

    return temp;
}

MA_node  operator* (const MA_node          &left_op, const   double         &right_op)
{
    MA_node temp(right_op);

    return ( left_op * temp );
}

MA_node  operator* (const   double         &left_op, const MA_node          &right_op)
{
    MA_node temp(left_op);

    return ( temp * right_op );
}

MA_node  operator/ (const MA_node          &left_op, const MA_node          &right_op)
{
    MA_node temp("var(100000,100000)/var(200000,200000)");

    temp.replace(left_op,100000,100000);
    temp.replace(right_op,200000,200000);

    return temp;
}

MA_node  operator/ (const MA_node          &left_op, const   double         &right_op)
{
    MA_node temp(right_op);

    return ( left_op / temp );
}

MA_node  operator/ (const   double         &left_op, const MA_node          &right_op)
{
    MA_node temp(left_op);

    return ( temp / right_op );
}




// *= multiplicative assignment - binary, return lvalue
// /= divisive       assignment - binary, return lvalue

MA_node &operator*=(      MA_node &left_op, const MA_node          &right_op)
{
    left_op = ( left_op * right_op );

    return left_op;
}

MA_node &operator*=(      MA_node &left_op, const   double         &right_op)
{
    left_op = ( left_op * right_op );

    return left_op;
}

MA_node &operator/=(      MA_node &left_op, const MA_node          &right_op)
{
    left_op = ( left_op / right_op );

    return left_op;
}

MA_node &operator/=(      MA_node &left_op, const   double         &right_op)
{
    left_op = ( left_op / right_op );

    return left_op;
}



int rime_finite(double num)
{
    return gsl_finite(num);
}
