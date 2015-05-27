
/*
 *  SVMheavy - Another SVM Library
 *  Copyright (C) 2005  Alistair Shilton
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
// Mercer Kernel functions
//
// Written by: Alistair Shilton
//             Melbourne University
//
// Date: 13-12-00 condensed from all over the place
//       6-12-01  fixed a few minor blunders
//       25-12-01 rewrote using maths.h parser
//       21-06-04 rewrote
//       31-10-04 optimisations done
//       **-09-05 more optimisations done
//

#ifndef _kernel_h
#define _kernel_h

#include <iostream>

#include "common/svdefs.h"
#include "common/vector.h"
#include "common/matrix.h"
#include "common/friends.h"

#define LINEAR_KERNEL           "var(5,12)"


/******************************************************************/
/*                                                                */
/* Simple kernel class                                            */
/*                                                                */
/******************************************************************/

class basic_kernel
{
    friend std::ostream &operator<<(std::ostream &output, const basic_kernel &dumpee);              \
    friend std::istream &operator>>(std::istream &input, basic_kernel &source);                     \

    //
    // Sparse Vector Format
    // ====================
    //
    // Sparse vectors with d non-zero elements must have the following
    // format:
    //
    // x = [ d  ]
    //     [ i1 ]
    //     [ v1 ]
    //     [ i2 ]
    //     [ v2 ]
    //     [ :  ]
    //     [ :  ]
    //     [ id ]
    //     [ vd ]
    //
    // where ij is the number of the jth non-zero element, and vj is the
    // value.  Hence the dimension of x must be 2d+2 (an even number).
    // Must have i1 < i2 < ... < id
    //
    //
    // Variable assignment
    // ===================
    //
    // The following variables may be used in the original expression:
    //
    // var(1,n) = xm vector (simplified out - not usable for sparse vectors)
    // var(2,n) = ym vector (simplified out - not usable for sparse vectors)
    //
    // var(3,n) = user variables (variable)
    // var(4,n) = user constants (fixed)
    //
    // var(5,1)  = xd
    // var(5,2)  = yd
    // var(5,3)  = xi
    // var(5,4)  = yi
    // var(5,5)  = xz
    // var(5,6)  = yz
    // var(5,7)  = covw (fixed if fix_covw == 1).
    // var(5,8)  = min(dim x,dim y)
    // var(5,9)  = max(dim x,dim y)
    // var(5,10) = reserved for differentiation
    // var(5,11) = reserved for differentiation
    // var(5,12) = x'y          (simplified out)
    // var(5,13) = (x-y)'(x-y)  (simplified out)
    // var(5,14) = ||x-y||      (simplified out)
    // var(5,15) = dim x
    // var(5,16) = dim y
    //
    // After initial processing, these will be expanded out (internally) to
    // the following:
    //
    // var(6,n) = x vector (not usable for sparse vectors)
    // var(7,n) = y vector (not usable for sparse vectors)
    //
    // var(8,n) = xc (mean correction) vector (fixed if fix_mean == 1).
    // var(9,n) = not used.
    //
    // var(9+m,n) = Covariance matrix (lower tri) (fixed if fix_covar == 1).
    //
    // The expansion itself works as follows:
    //
    // xm = covw.Ex.(x-xc)    (= covw.x for sparse vectors)
    // ym = covw.Ex.(y-xc)    (= covw.y for sparse vectors)
    //
    // where: xc = a vector
    //        Ex = a lower triangular matrix.
    //
    // where .* is the elementwise multiplication operation.
    //
    // NB: local variables var(0,n) where 1 <= n <= 12 are reserved.
    //
    // Sparse vectors note: var(6,n) and var(7,n) give access to the
    // relevant x or y structure - decoding this is your job.
    //
    // Other factors
    // =============
    //
    // scale_prod: if this is set then Ex/dim is used rather than Ex
    //
    // Derivatives
    // ===========
    //
    // The general derivative is w.r.t. var(var(5,10),var(5,11))
    //

    public:

    //
    // Constructors and Destructors
    // ============================
    //
    // default: uncorrected fixed linear kernel, covw = 1, scale_prod = 0
    // copy: copy the source
    // complete: all variables must be provided.
    //

    basic_kernel();
    basic_kernel(const basic_kernel &what);
    basic_kernel(char *kern_def, const fVECTOR &_uc, const fVECTOR &_uv, double _covw, const fVECTOR &_mean, const fMATRIX &_covar, int _fix_covw, int _fix_mean, int _fix_covar, int _scale_prod, int _sparse_kernel);
    basic_kernel(int kern_number, const fVECTOR &_uc, const fVECTOR &_uv, double _covw, const fVECTOR &_mean, const fMATRIX &_covar, int _fix_covw, int _fix_mean, int _fix_covar, int _scale_prod, int _sparse_kernel);

    basic_kernel &operator=(const basic_kernel &source);

    //
    // Kernel evaluation functions
    // ===========================
    //
    // Calculate the relevant result (after calculating the relevant
    // derivatives if this is the first call).
    //

    double kernel(const fVECTOR &x, const fVECTOR &y, L_DOUBLE *x_fast = NULL, L_DOUBLE *y_fast = NULL, double xd = 1, double yd = 1, long xi = -1, long yi = -1, double xz = 0, double yz = 0);
    double cross_kernel(const fVECTOR &x, const fVECTOR &y, double xd = 1, double yd = 1, long xi = -1, long yi = -1, double xz = 0, double yz = 0);
    basic_kernel kernel_deriv(const fVECTOR &x, const fVECTOR &y, double xd = 1, double yd = 1, long xi = -1, long yi = -1, double xz = 0, double yz = 0);

    //
    // Kernel manipulation functions
    // =============================
    //
    // The following functions operate on uv, covw, mean and covar, unless
    // the relevant fix_ flag is set.
    //
    // scale_kernel: scales variables.
    // add_kernel: assuming that what_scale is essentially the same kernel
    //             as this (ie. kernel_def is the same, any variables with
    //             the fix_ flag set are the same) then add a scaled version
    //             of the what variables to the local variables.
    //

    void scale_kernel(double scale);
    void add_kernel(const basic_kernel &what, double scale);

    //
    // Information
    // ===========
    //

    int is_sparse(void) const;

    //
    // IO stuff
    // ========
    //

    void prime_io(std::ostream *_where_to,int _echo_level);
    private:
    std::ostream *where_to;
    int echo_level;

    //
    // Variables
    // =========
    //

    MA_node kernel_def;
    MA_node cross_kernel_def;
    MA_node kernel_deriv_def;

    int is_cross_kernel_valid;
    int is_kernel_deriv_valid;

    int scale_prod;
    int sparse_kernel;

    public:
    fVECTOR uc;
    fVECTOR uv;

    double  covw;
    private:
    fVECTOR mean;
    fMATRIX covar;

    int fix_covw;
    int fix_mean;
    int fix_covar;

    //
    // Fast Kernel Functionality
    // =========================
    //
    // is_opt: indicates whether an optimalised kernel function is available
    //         If no, is_opt = 0.  Otherwise, is_opt indicates the built in
    //         kernel function.
    //

    public:
    int is_opt;
    double (*optim_kernel)(const fVECTOR &x, const fVECTOR &y, double xd, double yd, long xi, long yi, double xz, double yz, const fVECTOR &uc, const fVECTOR &uv, double covw, const fVECTOR &mean, const fMATRIX &covar);
    double (*optim_kernel_fast)(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw);
    private:

    //
    // Kernel control functions
    // ========================
    //

    void simplify_kernel(void);

    void calc_cross_kernel(long dim);
    void calc_kernel_deriv(void);
};




//
// Foldback function
// =================
//
// When MA_node needs to find a variable value, it will call back to this
// function to get the result.  argContents must actually point to an array
// of cast void pointers, which are, in order:
//
// argContents[1-1]  = NULL
// argContents[2-1]  = NULL
// argContents[3-1]  = uv
// argContents[4-1]  = uc
// argContents[5-1]  = an array of doubles that contain values of var(5,i)
// argContents[6-1]  = x
// argContents[7-1]  = y
// argContents[8-1]  = mean
// argContents[9-1]  = NULL
// argContents[10-1] = covar
// argContents[11-1] = scale_prod
//

double DoubleArgEvaluationFunctionKernel(void *argContents, long i, long j);



//
// Stream io
// =========
//

std::ostream &operator<<(std::ostream &output, const basic_kernel &dumpee);
std::istream &operator>>(std::istream &input,        basic_kernel &source);


//
// Operator overloading
// ====================
//
// These are just like the kernel manipulation functions described
// previously - the same provisos apply.
//
// +  posation                  - unary,  return rvalue
// -  negation                  - unary,  return rvalue
// +  addition                  - binary, return rvalue
// -  subtraction               - binary, return rvalue
// *  multiplication            - binary, return rvalue
// /  division                  - binary, return rvalue
// += additive       assignment - binary, return lvalue
// -= subtractive    assignment - binary, return lvalue
// *= multiplicative assignment - binary, return lvalue
// /= divisive       assignment - binary, return lvalue
//

basic_kernel  operator+ (const basic_kernel &left_op);
basic_kernel  operator- (const basic_kernel &left_op);
basic_kernel  operator+ (const basic_kernel &left_op, const basic_kernel &right_op);
basic_kernel  operator- (const basic_kernel &left_op, const basic_kernel &right_op);
basic_kernel  operator* (const basic_kernel &left_op, const double       &right_op);
basic_kernel  operator* (const double       &left_op, const basic_kernel &right_op);
basic_kernel  operator/ (const basic_kernel &left_op, const double       &right_op);
basic_kernel &operator+=(      basic_kernel &left_op, const basic_kernel &right_op);
basic_kernel &operator-=(      basic_kernel &left_op, const basic_kernel &right_op);
basic_kernel &operator*=(      basic_kernel &left_op, const double       &right_op);
basic_kernel &operator/=(      basic_kernel &left_op, const double       &right_op);





typedef double(*opt_kern_ptr)(L_DOUBLE *x, L_DOUBLE *y, L_DOUBLE *uc, L_DOUBLE *uv, double *covw);

/******************************************************************/
/*                                                                */
/* Kernel dictionary class.                                       */
/*                                                                */
/******************************************************************/

class kernel_dicto
{
    friend std::ostream &operator<<(std::ostream &output, const kernel_dicto &dumpee);
    friend std::istream &operator>>(std::istream &input, kernel_dicto &source);

    public:

    //
    // Constructors and Destructors
    // ============================
    //
    // default: uncorrected fixed linear kernel, covw = 1, scale_prod = 0
    // copy: copy the source
    // complete: all variables must be provided.
    //

    kernel_dicto();
    kernel_dicto(const kernel_dicto &what);
    kernel_dicto(long _dsize, int _fix_weights, double *weights, char **kern_def, fVECTOR **_uc, fVECTOR **_uv, double *_covw, fVECTOR **_mean, fMATRIX **_covar, int *_fix_covw, int *_fix_mean, int *_fix_covar, int *_scale_prod, int _sparse_kernel);
    kernel_dicto(long _dsize, int _fix_weights, double *weights, int *kern_number, fVECTOR **_uc, fVECTOR **_uv, double *_covw, fVECTOR **_mean, fMATRIX **_covar, int *_fix_covw, int *_fix_mean, int *_fix_covar, int *_scale_prod, int _sparse_kernel);
    kernel_dicto(double weight, int kern_number, const fVECTOR &_uc, const fVECTOR &_uv, double _covw, const fVECTOR &_mean, const fMATRIX &_covar, int _fix_covw = 1, int _fix_mean = 1, int _fix_covar = 1, int _scale_prod = 0, int _sparse_kernel = 0);

    kernel_dicto &operator=(const kernel_dicto &source);

    ~kernel_dicto();

    //
    // Normalisation
    // =============
    //
    // make kweights'.kweights = normal
    //

    void normalise(double normal);

    //
    // Kernel evaluation functions
    // ===========================
    //
    // Calculate the relevant result.
    //

    double kernel(const fVECTOR &x, const fVECTOR &y, L_DOUBLE *x_fast = NULL, L_DOUBLE *y_fast = NULL, double xd = 1.0, double yd = 1.0, long xi = -1, long yi = -1, double xz = 0, double yz = 0);
    double cross_kernel(const fVECTOR &x, const fVECTOR &y, double xd = 1.0, double yd = 1.0, long xi = -1, long yi = -1, double xz = 0, double yz = 0);
    kernel_dicto kernel_deriv(const fVECTOR &x, const fVECTOR &y, double xd = 1.0, double yd = 1.0, long xi = -1, long yi = -1, double xz = 0, double yz = 0);

    //
    // Optimisations
    // =============
    //
    // The fancy stuff here is really only necessary for very obscure cases
    // wherein kernel dictionaries or obscure kernel functions are used.  In
    // general, all that is needed is a very simple optimisation function.
    // The following provides that.
    //
    // To use, call get_fast_kern(uc,uv,covw).  If NULL is returned, then no
    // optimised version exists, so revert to using the usual kernel
    // function.  Otherwise, a pointer to the optimised function will be
    // returned.  uc, uv and covw must point to pointers which will be used
    // to store the relevant pointers for used when calling the optimised
    // kernel function.
    //

    opt_kern_ptr get_fast_kern(L_DOUBLE **uc, L_DOUBLE **uv, double **covw);

    //
    // Elementwise access functions
    // ============================
    //
    // While we can treat the kernel dictionary as one "element", we may
    // also wish to treat it as a set of separate kernels.  The following
    // functions allow this by providing access to individual elements of
    // the dictionary.
    //

    basic_kernel &kernel_element(long i);
    L_DOUBLE &weight_element(long j);
    long get_dsize(void);
    fVECTOR &getweights(void);

    //
    // Kernel manipulation functions
    // =============================
    //

    void scale_kernel(double scale);
    void add_kernel(const kernel_dicto &what, double scale);

    //
    // Information
    // ===========
    //

    int is_sparse(void) const;

    //
    // IO stuff
    // ========
    //

    void prime_io(std::ostream *_where_to,int _echo_level);
    private:
    std::ostream *where_to;
    int echo_level;

    public:

    //
    // data
    // ====
    //
    // dsize: number of kernels in dictionary.
    //
    // kern: kernels
    // kweights: corresponding weights
    //

    long dsize;

    basic_kernel **kern;
    fVECTOR kweights;

    int fix_weights;
};



//
// Stream io
// =========
//

std::ostream &operator<<(std::ostream &output, const kernel_dicto &dumpee);
std::istream &operator>>(std::istream &input,        kernel_dicto &source);


//
// Operator overloading
// ====================
//
// These are just like the kernel manipulation functions described
// previously - the same provisos apply.
//
// +  posation                  - unary,  return rvalue
// -  negation                  - unary,  return rvalue
// +  addition                  - binary, return rvalue
// -  subtraction               - binary, return rvalue
// *  multiplication            - binary, return rvalue
// /  division                  - binary, return rvalue
// += additive       assignment - binary, return lvalue
// -= subtractive    assignment - binary, return lvalue
// *= multiplicative assignment - binary, return lvalue
// /= divisive       assignment - binary, return lvalue
//

kernel_dicto  operator+ (const kernel_dicto &left_op);
kernel_dicto  operator- (const kernel_dicto &left_op);
kernel_dicto  operator+ (const kernel_dicto &left_op, const kernel_dicto &right_op);
kernel_dicto  operator- (const kernel_dicto &left_op, const kernel_dicto &right_op);
kernel_dicto  operator* (const kernel_dicto &left_op, const double       &right_op);
kernel_dicto  operator* (const double       &left_op, const kernel_dicto &right_op);
kernel_dicto  operator/ (const kernel_dicto &left_op, const double       &right_op);
kernel_dicto &operator+=(      kernel_dicto &left_op, const kernel_dicto &right_op);
kernel_dicto &operator-=(      kernel_dicto &left_op, const kernel_dicto &right_op);
kernel_dicto &operator*=(      kernel_dicto &left_op, const double       &right_op);
kernel_dicto &operator/=(      kernel_dicto &left_op, const double       &right_op);



#endif
