
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
// Sparse vector functions.
//
// Written by: Alistair Shilton
//             Melbourne University
//

#ifndef _sparsevector_h
#define _sparsevector_h

#include <iostream>
#include <string.h>
#include <math.h>



#include "svdefs.h"
#include "search.h"
#include "c_double.h"
#include "vector.h"



/*

Sparse vectors are supported via normal vectors having a specific format.
A sparse vector x with d non-zero elements must have format:

x = [ d  ]
    [ i1 ]
    [ v1 ]
    [ i2 ]
    [ v2 ]
    [ :  ]
    [ :  ]
    [ id ]
    [ vd ]

where ij is the number of the jth non-zero element, and vj is the value.
Hence the dimension of x must be 2d+2 (an even number).  To optimise speed
of operations, elements must be ordered, ie: i1 < i2 < ... < id



Functions:
==========

dot_prod_sparse(x,y) = x'y
dif_meas_sparse(x,y) = ||x-y||^2
eucldist_sparse(x,y) = ||x-y||

vector_scale(x,a): x := ax

vector_add(x,y): x := x+y
vector_sub(x,y): x := x-y


Constructions
=============

convert_standard_format: converts a string in "standard" format, ie:
   i1:v1 i2:v2 ... in:vn
   to a sparse vector satisfying the above format
desparsifier: convert a set of sparse vectors to nonsparse vectors
   with consistent dimension (which is the minimum dimension that
   can store all the sparse vectors).

*/


double dot_prod_sparse(const fVECTOR &x, const fVECTOR &y);
double dif_meas_sparse(const fVECTOR &x, const fVECTOR &y);
double eucldist_sparse(const fVECTOR &x, const fVECTOR &y);

double dot_prod_sparse_opt(L_DOUBLE *x, L_DOUBLE *y);
double dif_meas_sparse_opt(L_DOUBLE *x, L_DOUBLE *y);
double eucldist_sparse_opt(L_DOUBLE *x, L_DOUBLE *y);

fVECTOR &vector_scale_sparse(fVECTOR &x, double a);

fVECTOR &vector_add_sparse(fVECTOR &x, const fVECTOR &y);
fVECTOR &vector_sub_sparse(fVECTOR &x, const fVECTOR &y);


fVECTOR &convert_standard_format(char *input, fVECTOR &result);
fVECTOR &convert_full_format(char *input, fVECTOR &result);
f_fVECTOR &desparsifier(f_fVECTOR &input, f_fVECTOR &result);


#endif
