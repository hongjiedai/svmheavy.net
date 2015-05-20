
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
// Support vector machine - regression class
//
// Written by: Alistair Shilton
//             Melbourne University
//


#ifndef _svm_regress_h
#define _svm_regress_h

#include <iostream>
#include <limits.h>
#include "svoptim.h"
#include "svflags.h"
#include "kernel.h"
#include "common/svdefs.h"
#include "common/vector.h"

class SVM_regress;
class SVdata;

#define DEFAULT_SOL_TOL         1e-3
#define MAX_UPPER_BOUND         1e200
#define DEFAULT_CACHE_SIZE      64
#define DEFAULT_TRAININGSIZE    500
#define DEFAULT_D2C_BUFF_LEN    10

void *background_opt_callback_regress(void *caller);

//
// Regression class wrapper.
//

std::ostream &operator<<(std::ostream &output, SVM_regress &dumpee);
std::istream &operator>>(std::istream &input,  SVM_regress &source);


class SVM_regress
{
    friend std::ostream &operator<<(std::ostream &, SVM_regress &);
    friend std::istream &operator>>(std::istream &, SVM_regress &);

    friend void *background_opt_callback_regress(void *caller);

    public:

    //
    // Constructor.
    //
    // kern: kernel dictionary.
    //
    // risk_type: 0 = linear empirical misclassification error
    //            1 = quadratic empirical misclassification error
    //
    // tube_shrink: 0 = no tube shrinking
    //              1 = quadratic tube shrinking.
    //
    // C:             Misclassification cost.
    // C_plus_scale:  scale for C for d=+1 vectors.
    // C_minus_scale: scale for C for d=-1 vectors.
    //
    // E:             epsilon value.
    // E_plus_scale:  epsilon scale for overshoot.
    // E_minus_scale: epsilon scale for undershoot.
    //
    // nu: quadratic tube shrinking term.
    //
    // fix_bias: 0 = variable bias (standard)
    //           1 = fixed bias
    //
    // default_bias: b value used if fix_bias == 1
    //
    // opt_type: see above SVM_OPT_* macros for details.
    //
    // epochs:     maximum number of epochs for any given algorithm.
    // sol_tol:    solution tolerance requirement for gradient descent/SMO
    // cachesize:  if dynamic caching used, this is its size (in MB)
    // Nexpect:    expected training set size (used by dynamic caching)
    // d2cbufflen: length of error buffer used by d2c optimiser
    //

    SVM_regress();
    SVM_regress(const SVM_regress &source);
    SVM_regress(const kernel_dicto &_kern,
                int risk_type = 0,
                int tube_shrink = 0,
                L_DOUBLE C = 1.0,
                L_DOUBLE C_plus_scale = 1.0,
                L_DOUBLE C_minus_scale = 1.0,
                L_DOUBLE E = 1.0,
                L_DOUBLE E_plus_scale = 1.0,
                L_DOUBLE E_minus_scale = 1.0,
                L_DOUBLE nu = 0.1,
                int fix_bias = 0,
                L_DOUBLE default_bias = 0.0,
                int opt_type = 1,
                long epochs = LONG_MAX,
                double sol_tol = DEFAULT_SOL_TOL,
                long cachesize = DEFAULT_CACHE_SIZE,
                long Nexpect = DEFAULT_TRAININGSIZE,
                long d2cbufflen = DEFAULT_D2C_BUFF_LEN);

    SVM_regress &operator=(const SVM_regress &what);

    //
    // Destructor.
    //

    ~SVM_regress();
   
    //
    // Background optimisation options.
    //
    // On some multi-tasking architectures, it is possible for SVMheavy to
    // optimise the SVM as a background process, so the program using the
    // SVM_regress can add data, change parameters, find g(x) simultaneously
    // with the optimisation process (although naturally, g(x) will not
    // become optimal until optimisation is complete).  This ability is
    // controlled by the following functions.
    //
    // back_opt_on    - turn on background optimisation.  Will return 0 if
    //                  successful, or an error code otherwise.
    // back_opt_off   - turn off background optimisation.
    // is_back_opt    - returns 1 if background optimisation is turned on,
    //                  or 0 otherwise.
    // get_iter       - returns the contents of the iteration counter.  This
    //                  will be 0 if optimisation is not in process, nonzero
    //                  otherwise.
    // get_iter_accum - returns the number of iterations (total) taken by the
    //                  all background optimisation processes since the 
    //                  object was constructed (or set_iter_accum was last
    //                  called).
    // set_iter_accum - set the cumulative iteration counter.
    // get_J          - gets the value of the objective function.  This is
    //                  not necessarily accurate - it is only ever reset on
    //                  construction, will not be changed by addition or
    //                  removal of data or modification of training
    //                  parameters, and is only kept by the variable bias
    //                  SMO and variable bias D2C optimisers.
    // set_J          - set the "value" (reported) of the objective function
    //                  reported by get_J.
    //

    int back_opt_on(void);
    void back_opt_off(void);
    int is_back_opt(void) const;
    unsigned long get_iter(void) const;
    unsigned long get_iter_accum(void) const;
    void set_iter_accum(unsigned long itval);
    L_DOUBLE get_J(void) const;
    void set_J(L_DOUBLE Jval);

    //
    // Training set modification functions.
    //
    // Add/remove/modify individual training points:
    //
    // add_point:  Add point to the training set.
    // add_points: Add multiple points to the training set.
    // del_point:  Delete point with marker id from training set.
    // del_points: Delete multiple points.
    //

    void add_point(L_DOUBLE z, fVECTOR &x, L_DOUBLE t, L_DOUBLE t_star, L_DOUBLE eps, L_DOUBLE eps_star, long id);
    void add_points(fVECTOR &z, f_fVECTOR &x, fVECTOR &t, fVECTOR &t_star, fVECTOR &eps, fVECTOR &eps_star, iVECTOR &id);

    void add_point_lower_bound(L_DOUBLE z, fVECTOR &x, L_DOUBLE t, L_DOUBLE eps, long id);
    void add_points_lower_bound(fVECTOR &z, f_fVECTOR &x, fVECTOR &t, fVECTOR &eps, iVECTOR &id);

    void add_point_upper_bound(L_DOUBLE z, fVECTOR &x, L_DOUBLE t_star, L_DOUBLE eps_star, long id);
    void add_points_upper_bound(fVECTOR &z, f_fVECTOR &x, fVECTOR &t_star, fVECTOR &eps_star, iVECTOR &id);

    void del_point(long id);
    void del_points(iVECTOR id);

    //
    // Grid helper functions.
    //
    // set_z(z) - sets "z".  For this case, z is simply the gradient offset
    //            term, negated
    //

    void set_z(long which, L_DOUBLE z);

    //
    // Training parameter modification functions:
    //

    void set_C(L_DOUBLE C_new);
    void set_E(L_DOUBLE E_new);
    void set_nu(L_DOUBLE nu_new);
    void set_C_plus_scale(L_DOUBLE C_new);
    void set_C_minus_scale(L_DOUBLE C_new);
    void set_E_plus_scale(L_DOUBLE E_new);
    void set_E_minus_scale(L_DOUBLE E_new);

    void scale_C(L_DOUBLE C_scale);
    void scale_E(L_DOUBLE E_scale);
    void scale_nu(L_DOUBLE nu_scale);
    void scale_C_plus_scale(L_DOUBLE C_scale);
    void scale_C_minus_scale(L_DOUBLE C_scale);
    void scale_E_plus_scale(L_DOUBLE E_scale);
    void scale_E_minus_scale(L_DOUBLE E_scale);

    void set_kernel(kernel_dicto &kerndict);
    void set_epochs(long _epochs);

    //
    // Information functions
    //

    L_DOUBLE get_C(void);
    L_DOUBLE get_C_plus_scale(void);
    L_DOUBLE get_C_minus_scale(void);
    L_DOUBLE get_E(void);
    L_DOUBLE get_E_plus_scale(void);
    L_DOUBLE get_E_minus_scale(void);
    L_DOUBLE get_nu(void);
    kernel_dicto &get_kernel(void);
    fVECTOR &get_x(long i);

    //
    // Optimise: Optimise the SVM, returning the iteration count.
    //

    unsigned long optimise(void);
    unsigned long optimise(L_DOUBLE temp_sol_tol);
    unsigned long optimise_from_zero(L_DOUBLE temp_sol_tol, long MaxItersPerGeneration);

    void removeNonSupports(void);

    //
    // Use SVM: test_point returns g(x)
    //          cross_test_poing returns dg(x)/dx
    //

    L_DOUBLE test_point(const fVECTOR &x);
    L_DOUBLE cross_test_point(const fVECTOR &x);

    //
    // Get information.
    //

    long get_N(void);
    long get_N_Z(void);
    long get_N_L(void);
    long get_N_U(void);
    long get_N_F(void);
    long get_N_FP(void);
    long get_N_FN(void);
    long get_N_S(void);
    long get_max_id(void);
    fVECTOR get_alpha(void);
    iVECTOR get_tau(void);







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

    SVdata *raw;

    int risk_type;
    int tube_shrink;
    int opt_type;

    L_DOUBLE C;
    L_DOUBLE C_plus_scale;
    L_DOUBLE C_minus_scale;

    L_DOUBLE E;
    L_DOUBLE E_plus_scale;
    L_DOUBLE E_minus_scale;

    L_DOUBLE nu;

    L_DOUBLE sol_tol;
    unsigned long epochs;

    unsigned long optimise_unsafe(void);

    volatile int async_exit_flag;
    volatile int back_opt_flag;
    volatile unsigned long loop_count;
    volatile unsigned long loop_count_accum;
    JTYPE J;
};

#endif
